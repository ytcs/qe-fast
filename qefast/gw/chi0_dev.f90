! Copyright (C) 2026 The qe-fast authors. Portions derive from Quantum ESPRESSO; see NOTICE.
! SPDX-License-Identifier: GPL-2.0-or-later
!
! chi0_dev -- device path of chi0.x and gamma_ii.x (CUDA builds only). The stages of the host path, on plain arrays:
!   pairs   rho_t(G) for all conduction bands of one (k, v) by one batched cuFFT, kept in a device buffer
!   bin     buffered transitions sorted by energy bin on the host, one cuBLAS ZGEMM per occupied bin into A_b
!   kk      chi0 on the energy grid: zero-padded convolution and correlation along the bins, strided 1D cuFFT
!   inv     eps = 1 - v^1/2 chi0 v^1/2 and eps^-1, batched LU and inverse (cuBLAS), returned in chunks
! gamma_ii.x keeps the plane-wave coefficients of every k on the device (read once) and builds the periodic parts
! there; per q it turns eps^-1 into Im W in place (imw) and contracts it with the pair densities of the hot
! carriers (contract): one cuBLAS ZGEMM per energy of the grid over the pairs that interpolate at it.
! The host paths in chi0.f90 and gamma_ii.f90 are the reference for every stage.
!
MODULE chi0_dev
  USE cudafor
  USE cufft
  USE cublas_v2
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 200)
  INTEGER, PARAMETER :: TPB = 256
  PUBLIC :: dev_init, dev_add_k, dev_flush, dev_kk, dev_inv_chunk, dev_free, dev_cap
  PUBLIC :: dev_free_gb, dev_wfc_init, dev_wfc_put, dev_wfc_free, dev_add_kq, dev_imw, dev_contract
  !
  INTEGER, SAVE :: nr(3), nrx(3), nnr, npol, nc, ng, nb, cap, nt, ntot
  INTEGER, SAVE :: plan3 = 0, nplan3 = 0
  TYPE(cublasHandle), SAVE :: hb
  INTEGER, DEVICE, ALLOCATABLE, SAVE :: d_idx(:)
  COMPLEX(dp), DEVICE, ALLOCATABLE, SAVE :: d_uk(:,:,:), d_uq(:,:,:), d_f(:,:), d_r(:,:), d_a(:,:), d_chi(:,:)
  REAL(dp), ALLOCATABLE, SAVE :: h_x(:)
  ! gamma_ii.x: coefficients of every k (npwx npol, nbnd, nks), their Miller indices, cuFFT plans by batch size
  INTEGER, SAVE :: npwc = 0, nbndc = 0, nksc = 0
  INTEGER, DEVICE, ALLOCATABLE, SAVE :: d_mill(:,:,:), d_ngk(:)
  COMPLEX(dp), DEVICE, ALLOCATABLE, SAVE :: d_coef(:,:,:), d_bf(:,:,:), d_bi(:,:,:), d_ff(:,:), d_phi(:,:)
  INTEGER, SAVE :: pl_batch(8) = 0, pl_plan(8) = 0
  !
CONTAINS
  !
  SUBROUTINE chk(istat, what)
    INTEGER, INTENT(IN) :: istat
    CHARACTER(*), INTENT(IN) :: what
    IF (istat /= 0) THEN
      WRITE(*, '(a, a, a, i0)') 'chi0_dev: ', what, ' failed, status ', istat
      STOP 1
    ENDIF
  END SUBROUTINE chk
  !
  INTEGER FUNCTION dev_cap(nnr_, npol_, nv_, nc_, ng_, nb_, nw_, reserve_gb)
    !! transitions that the device buffer can hold after the fixed arrays, or 0 if even those do not fit
    INTEGER, INTENT(IN) :: nnr_, npol_, nv_, nc_, ng_, nb_, nw_
    REAL(dp), INTENT(IN) :: reserve_gb
    INTEGER(8) :: freeb, totb, fixed, per
    INTEGER :: istat
    istat = cudaMemGetInfo(freeb, totb)
    fixed = 16_8 * nnr_ * (npol_ * (nv_ + nc_) + nc_) + 16_8 * ng_ * ng_ * (nb_ + nw_)
    per = 16_8 * ng_ * 5             ! rho, and at flush two weighted copies with two entries per transition
    dev_cap = INT(MAX(0_8, (freeb - fixed - INT(reserve_gb * 1.0d9, 8)) / per))
    dev_cap = MIN(dev_cap, 4 * 1024 * 1024)
  END FUNCTION dev_cap
  !
  SUBROUTINE dev_init(nr_, nrx_, nnr_, npol_, nv_, nc_, ng_, idx_, nb_, nw_, cap_)
    !! idx_(ig): position of +G on the FFT box (1-based); rho is taken from the unnormalised inverse transform / N
    INTEGER, INTENT(IN) :: nr_(3), nrx_(3), nnr_, npol_, nv_, nc_, ng_, idx_(ng_), nb_, nw_, cap_
    INTEGER :: istat, n(3), emb(3)
    nr = nr_; nrx = nrx_; nnr = nnr_; npol = npol_; nc = nc_; ng = ng_; nb = nb_; cap = cap_; nt = 0; ntot = 0
    istat = cublasCreate(hb); CALL chk(istat, 'cublasCreate')
    ALLOCATE(d_idx(ng), d_uk(nnr, npol, nv_), d_uq(nnr, npol, nc), d_f(nnr, nc), d_r(ng, cap), h_x(cap))
    ALLOCATE(d_a(ng * ng, nb), d_chi(ng * ng, nw_))
    d_idx = idx_
    d_a = (0.0_dp, 0.0_dp)
    n = (/ nr(3), nr(2), nr(1) /)
    emb = (/ nrx(3), nrx(2), nrx(1) /)
    istat = cufftPlanMany(plan3, 3, n, emb, 1, nnr, emb, 1, nnr, CUFFT_Z2Z, nc); CALL chk(istat, 'cufftPlanMany 3d')
  END SUBROUTINE dev_init
  !
  ATTRIBUTES(GLOBAL) SUBROUTINE k_prod(f, uq, uk, nnr, npol, nc, nvk, iv)
    ! f(r, c) = sum_s conj(uq(r, s, c)) uk(r, s, iv)
    INTEGER, VALUE :: nnr, npol, nc, nvk, iv
    COMPLEX(dp) :: f(nnr, nc), uq(nnr, npol, nc), uk(nnr, npol, nvk)
    INTEGER :: r, c, s
    COMPLEX(dp) :: t
    r = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    c = blockIdx%y
    IF (r > nnr) RETURN
    t = (0.0_dp, 0.0_dp)
    DO s = 1, npol
      t = t + CONJG(uq(r, s, c)) * uk(r, s, iv)
    ENDDO
    f(r, c) = t
  END SUBROUTINE k_prod
  !
  ATTRIBUTES(GLOBAL) SUBROUTINE k_gather(rr, f, idx, ng, nnr, nc, off, scal)
    ! rr(ig, off + c) = f(idx(ig), c) * scal
    INTEGER, VALUE :: ng, nnr, nc, off
    REAL(dp), VALUE :: scal
    COMPLEX(dp) :: rr(ng, *), f(nnr, nc)
    INTEGER :: idx(ng)
    INTEGER :: ig, c
    ig = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    c = blockIdx%y
    IF (ig > ng) RETURN
    rr(ig, off + c) = f(idx(ig), c) * scal
  END SUBROUTINE k_gather
  !
  SUBROUTINE dev_add_k(uk, uq, nv, x)
    !! uk(:,:,v) (k) and uq(:,:,c) (k+q) on the FFT box; x(c, v) = (e_c - e_v) / de, the bin coordinate
    INTEGER, INTENT(IN) :: nv
    COMPLEX(dp), INTENT(IN) :: uk(nnr, npol, nv), uq(nnr, npol, nc)
    REAL(dp), INTENT(IN) :: x(nc, nv)
    d_uk = uk
    d_uq = uq
    CALL add_pairs(nv, x)
  END SUBROUTINE dev_add_k
  !
  SUBROUTINE add_pairs(nv, x)
    !! the transitions of d_uk(:,:,1:nv) -> d_uq into the buffer (flushed into A_b when full)
    INTEGER, INTENT(IN) :: nv
    REAL(dp), INTENT(IN) :: x(nc, nv)
    INTEGER :: iv, istat
    TYPE(dim3) :: g
    DO iv = 1, nv
      IF (nt + nc > cap) CALL dev_flush()
      g = dim3((nnr + TPB - 1) / TPB, nc, 1)
      CALL k_prod<<<g, TPB>>>(d_f, d_uq, d_uk, nnr, npol, nc, nv, iv)
      istat = cufftExecZ2Z(plan3, d_f, d_f, CUFFT_INVERSE); CALL chk(istat, 'cufftExecZ2Z pairs')
      g = dim3((ng + TPB - 1) / TPB, nc, 1)
      CALL k_gather<<<g, TPB>>>(d_r, d_f, d_idx, ng, nnr, nc, nt, 1.0_dp / (nr(1) * nr(2) * nr(3)))
      h_x(nt + 1:nt + nc) = x(:, iv)
      nt = nt + nc
    ENDDO
  END SUBROUTINE add_pairs
  !
  ATTRIBUTES(GLOBAL) SUBROUTINE k_weighted(xm, ym, rr, ent, wt, ng, ne)
    ! xm(:, e) = conj(rr(:, ent(e))) wt(e),  ym(:, e) = rr(:, ent(e))
    INTEGER, VALUE :: ng, ne
    COMPLEX(dp) :: xm(ng, ne), ym(ng, ne), rr(ng, *)
    INTEGER :: ent(ne)
    REAL(dp) :: wt(ne)
    INTEGER :: ig, e
    INTEGER(8) :: t
    t = INT(blockIdx%x - 1, 8) * blockDim%x + threadIdx%x            ! 1D grid: ne can exceed the y limit
    IF (t > INT(ng, 8) * ne) RETURN
    e = INT((t - 1) / ng) + 1
    ig = INT(t - INT(e - 1, 8) * ng)
    xm(ig, e) = CONJG(rr(ig, ent(e))) * wt(e)
    ym(ig, e) = rr(ig, ent(e))
  END SUBROUTINE k_weighted
  !
  SUBROUTINE dev_flush()
    !! bin the buffered transitions: each goes to bins floor(x)+1 and floor(x)+2 with weights 1-f and f (as the host)
    INTEGER, ALLOCATABLE :: cnt(:), off(:), ent(:), pos(:)
    REAL(dp), ALLOCATABLE :: wt(:)
    INTEGER, DEVICE, ALLOCATABLE :: d_ent(:)
    REAL(dp), DEVICE, ALLOCATABLE :: d_wt(:)
    COMPLEX(dp), DEVICE, ALLOCATABLE :: d_xm(:,:), d_ym(:,:)
    INTEGER :: t, j0, ne, b, istat, e
    REAL(dp) :: f
    TYPE(dim3) :: g
    IF (nt == 0) RETURN
    ALLOCATE(cnt(nb + 1), off(nb + 1), pos(nb + 1))
    cnt = 0
    DO t = 1, nt
      j0 = FLOOR(h_x(t)) + 1
      IF (j0 + 1 <= nb) THEN
        cnt(j0) = cnt(j0) + 1; cnt(j0 + 1) = cnt(j0 + 1) + 1
      ENDIF
    ENDDO
    off(1) = 0
    DO b = 2, nb + 1
      off(b) = off(b - 1) + cnt(b - 1)
    ENDDO
    ne = off(nb + 1)
    ALLOCATE(ent(MAX(ne, 1)), wt(MAX(ne, 1)))
    pos = off
    DO t = 1, nt
      j0 = FLOOR(h_x(t)) + 1
      f = h_x(t) - FLOOR(h_x(t))
      IF (j0 + 1 <= nb) THEN
        pos(j0) = pos(j0) + 1; ent(pos(j0)) = t; wt(pos(j0)) = 1.0_dp - f
        pos(j0 + 1) = pos(j0 + 1) + 1; ent(pos(j0 + 1)) = t; wt(pos(j0 + 1)) = f
      ENDIF
    ENDDO
    IF (ne > 0) THEN
      ALLOCATE(d_ent(ne), d_wt(ne), d_xm(ng, ne), d_ym(ng, ne))
      d_ent = ent(1:ne); d_wt = wt(1:ne)
      g = dim3(INT((INT(ng, 8) * ne + TPB - 1) / TPB), 1, 1)
      CALL k_weighted<<<g, TPB>>>(d_xm, d_ym, d_r, d_ent, d_wt, ng, ne)
      DO b = 1, nb
        IF (cnt(b) == 0) CYCLE
        e = off(b) + 1
        istat = cublasZgemm(hb, CUBLAS_OP_N, CUBLAS_OP_T, ng, ng, cnt(b), (1.0_dp, 0.0_dp), d_xm(1, e), ng, &
                            d_ym(1, e), ng, (1.0_dp, 0.0_dp), d_a(1, b), ng)
        CALL chk(istat, 'cublasZgemm bin')
      ENDDO
      DEALLOCATE(d_ent, d_wt, d_xm, d_ym)
    ENDIF
    ntot = ntot + nt
    nt = 0
    DEALLOCATE(cnt, off, pos, ent, wt)
  END SUBROUTINE dev_flush
  !
  ATTRIBUTES(GLOBAL) SUBROUTINE k_load(p, a, ncol, ldp, l, nb, c0, m, rev)
    ! p(c, j) = a(c0 + c - 1, j) for j <= nb (reversed in j if rev), 0 up to l; p has leading dimension ldp
    INTEGER, VALUE :: ncol, ldp, l, nb, c0, m, rev
    COMPLEX(dp) :: p(ldp, l), a(ncol, nb)
    INTEGER :: c, j
    c = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    j = blockIdx%y
    IF (c > m) RETURN
    IF (j > nb) THEN
      p(c, j) = (0.0_dp, 0.0_dp)
    ELSE IF (rev == 1) THEN
      p(c, j) = a(c0 + c - 1, nb + 1 - j)
    ELSE
      p(c, j) = a(c0 + c - 1, j)
    ENDIF
  END SUBROUTINE k_load
  !
  ATTRIBUTES(GLOBAL) SUBROUTINE k_mulrow(p, kf, ldp, l, m)
    INTEGER, VALUE :: ldp, l, m
    COMPLEX(dp) :: p(ldp, l), kf(l)
    INTEGER :: c, j
    c = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    j = blockIdx%y
    IF (c > m) RETURN
    p(c, j) = p(c, j) * kf(j)
  END SUBROUTINE k_mulrow
  !
  ATTRIBUTES(GLOBAL) SUBROUTINE k_store(chi, p, ncol, ldp, l, nb, nw, c0, m, sgn, scal)
    ! chi(c0 + c - 1, i) (+)= sgn * scal * p(c, nb + i - 1)
    INTEGER, VALUE :: ncol, ldp, l, nb, nw, c0, m
    REAL(dp), VALUE :: sgn, scal
    COMPLEX(dp) :: chi(ncol, nw), p(ldp, l)
    INTEGER :: c, i
    c = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    i = blockIdx%y
    IF (c > m) RETURN
    IF (sgn > 0.0_dp) THEN
      chi(c0 + c - 1, i) = scal * p(c, nb + i - 1)
    ELSE
      chi(c0 + c - 1, i) = chi(c0 + c - 1, i) - scal * p(c, nb + i - 1)
    ENDIF
  END SUBROUTINE k_store
  !
  SUBROUTINE dev_kk(pref, de_ha, eta_ha, nw)
    !! d_chi(:, i) = pref * sum_b d_a(:, b) [K(i - b) - K(i + b)],  K(m) = 1/(m de + i eta), by strided 1D cuFFTs
    REAL(dp), INTENT(IN) :: pref, de_ha, eta_ha
    INTEGER, INTENT(IN) :: nw
    COMPLEX(dp), ALLOCATABLE :: k1(:), k2(:)
    COMPLEX(dp), DEVICE, ALLOCATABLE :: d_k1(:), d_k2(:), d_p(:,:)
    INTEGER :: l, jj, ncol, mb, c0, m, istat, p1, pb
    INTEGER(8) :: freeb, totb
    TYPE(dim3) :: g
    CALL dev_flush()
    DEALLOCATE(d_uk, d_uq, d_f, d_r)                  ! pair stage done: room for the convolution blocks
    l = 1
    DO WHILE (l < 2 * nb + nw)
      l = 2 * l
    ENDDO
    ALLOCATE(k1(l), k2(l), d_k1(l), d_k2(l))
    k1 = (0.0_dp, 0.0_dp); k2 = (0.0_dp, 0.0_dp)
    DO jj = 1, nw + nb - 1
      k1(jj) = 1.0_dp / CMPLX(((jj - 1) - (nb - 1)) * de_ha, eta_ha, KIND=dp)
      k2(jj) = 1.0_dp / CMPLX((jj - 1) * de_ha, eta_ha, KIND=dp)
    ENDDO
    d_k1 = k1; d_k2 = k2
    istat = cufftPlan1d(p1, l, CUFFT_Z2Z, 1); CALL chk(istat, 'cufftPlan1d')
    istat = cufftExecZ2Z(p1, d_k1, d_k1, CUFFT_FORWARD); CALL chk(istat, 'k1 fft')
    istat = cufftExecZ2Z(p1, d_k2, d_k2, CUFFT_FORWARD); CALL chk(istat, 'k2 fft')
    istat = cufftDestroy(p1)
    ncol = ng * ng
    istat = cudaMemGetInfo(freeb, totb)
    mb = INT(MIN(INT(ncol, 8), MAX(1_8, (freeb / 3) / (16_8 * l))))      ! columns per block
    ALLOCATE(d_p(mb, l))
    DO c0 = 1, ncol, mb
      m = MIN(mb, ncol - c0 + 1)
      ! batched 1D transforms along j (stride mb, batch distance 1)
      istat = cufftPlanMany(pb, 1, (/ l /), (/ l /), mb, 1, (/ l /), mb, 1, CUFFT_Z2Z, m)
      CALL chk(istat, 'cufftPlanMany 1d')
      ! convolution term
      g = dim3((m + TPB - 1) / TPB, l, 1)
      CALL k_load<<<g, TPB>>>(d_p, d_a, ncol, mb, l, nb, c0, m, 0)
      istat = cufftExecZ2Z(pb, d_p, d_p, CUFFT_FORWARD); CALL chk(istat, 'kk fwd')
      CALL k_mulrow<<<g, TPB>>>(d_p, d_k1, mb, l, m)
      istat = cufftExecZ2Z(pb, d_p, d_p, CUFFT_INVERSE); CALL chk(istat, 'kk inv')
      g = dim3((m + TPB - 1) / TPB, nw, 1)
      CALL k_store<<<g, TPB>>>(d_chi, d_p, ncol, mb, l, nb, nw, c0, m, 1.0_dp, pref / l)
      ! correlation term
      g = dim3((m + TPB - 1) / TPB, l, 1)
      CALL k_load<<<g, TPB>>>(d_p, d_a, ncol, mb, l, nb, c0, m, 1)
      istat = cufftExecZ2Z(pb, d_p, d_p, CUFFT_FORWARD); CALL chk(istat, 'kk fwd 2')
      CALL k_mulrow<<<g, TPB>>>(d_p, d_k2, mb, l, m)
      istat = cufftExecZ2Z(pb, d_p, d_p, CUFFT_INVERSE); CALL chk(istat, 'kk inv 2')
      g = dim3((m + TPB - 1) / TPB, nw, 1)
      CALL k_store<<<g, TPB>>>(d_chi, d_p, ncol, mb, l, nb, nw, c0, m, -1.0_dp, pref / l)
      istat = cufftDestroy(pb)
    ENDDO
    DEALLOCATE(d_p, d_k1, d_k2, k1, k2)
    DEALLOCATE(d_a)
  END SUBROUTINE dev_kk
  !
  ATTRIBUTES(GLOBAL) SUBROUTINE k_eps(e, chi, sv, ng, i0, n, ncol)
    ! e(:, :, w) = delta - sv chi(:, i0 + w - 1) sv,  and an identity right-hand side is built by the caller
    INTEGER, VALUE :: ng, i0, n, ncol
    COMPLEX(dp) :: e(ng, ng, n), chi(ncol, *)
    REAL(dp) :: sv(ng)
    INTEGER :: i, j, w
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    j = blockIdx%y
    w = blockIdx%z
    IF (i > ng) RETURN
    e(i, j, w) = - sv(i) * chi(i + (j - 1) * ng, i0 + w - 1) * sv(j)
    IF (i == j) e(i, j, w) = e(i, j, w) + 1.0_dp
  END SUBROUTINE k_eps
  !
  SUBROUTINE dev_inv_chunk(sv, i0, n, eps_h, einv_h)
    !! eps and eps^-1 for energies i0 .. i0+n-1, returned to the host; batched LU and inverse (cuBLAS)
    REAL(dp), INTENT(IN) :: sv(ng)
    INTEGER, INTENT(IN) :: i0, n
    COMPLEX(dp), INTENT(OUT) :: eps_h(ng, ng, n), einv_h(ng, ng, n)
    REAL(dp), DEVICE, ALLOCATABLE :: d_sv(:)
    COMPLEX(dp), DEVICE, ALLOCATABLE, TARGET :: d_e(:,:,:), d_c(:,:,:)
    INTEGER, DEVICE, ALLOCATABLE :: d_piv(:), d_info(:)
    TYPE(c_devptr), ALLOCATABLE :: pa(:), pc(:)
    TYPE(c_devptr), DEVICE, ALLOCATABLE :: d_pa(:), d_pc(:)
    INTEGER, ALLOCATABLE :: info(:)
    INTEGER :: w, istat
    TYPE(dim3) :: g
    ALLOCATE(d_sv(ng), d_e(ng, ng, n), d_c(ng, ng, n), d_piv(ng * n), d_info(n), pa(n), pc(n), d_pa(n), d_pc(n), info(n))
    d_sv = sv
    g = dim3((ng + TPB - 1) / TPB, ng, n)
    CALL k_eps<<<g, TPB>>>(d_e, d_chi, d_sv, ng, i0, n, ng * ng)
    eps_h = d_e
    DO w = 1, n
      pa(w) = c_devloc(d_e(1, 1, w)); pc(w) = c_devloc(d_c(1, 1, w))
    ENDDO
    d_pa = pa; d_pc = pc
    istat = cublasZgetrfBatched(hb, ng, d_pa, ng, d_piv, d_info, n); CALL chk(istat, 'getrfBatched')
    info = d_info
    IF (ANY(info /= 0)) CALL chk(MAXVAL(ABS(info)), 'singular eps (getrfBatched info)')
    istat = cublasZgetriBatched(hb, ng, d_pa, ng, d_piv, d_pc, ng, d_info, n); CALL chk(istat, 'getriBatched')
    einv_h = d_c
    DEALLOCATE(d_sv, d_e, d_c, d_piv, d_info, pa, pc, d_pa, d_pc, info)
  END SUBROUTINE dev_inv_chunk
  !
  SUBROUTINE dev_free()
    INTEGER :: istat
    IF (ALLOCATED(d_chi)) DEALLOCATE(d_chi)
    IF (ALLOCATED(d_a)) DEALLOCATE(d_a)
    IF (ALLOCATED(d_uk)) DEALLOCATE(d_uk, d_uq, d_f, d_r)
    IF (ALLOCATED(d_idx)) DEALLOCATE(d_idx, h_x)
    IF (ALLOCATED(d_bf)) DEALLOCATE(d_bf, d_bi, d_ff, d_phi)
    istat = cufftDestroy(plan3)
    istat = cublasDestroy(hb)
  END SUBROUTINE dev_free
  !
  ! ---------------------------------------------------------------- gamma_ii.x
  !
  REAL(dp) FUNCTION dev_free_gb()
    INTEGER(8) :: freeb, totb
    INTEGER :: istat
    istat = cudaMemGetInfo(freeb, totb)
    dev_free_gb = freeb / 1.0d9
  END FUNCTION dev_free_gb
  !
  SUBROUTINE dev_wfc_init(npw_, nbnd_, nks_)
    !! room for the coefficients of every k: npw_ = npwx * npol rows
    INTEGER, INTENT(IN) :: npw_, nbnd_, nks_
    npwc = npw_; nbndc = nbnd_; nksc = nks_
    ALLOCATE(d_coef(npwc, nbndc, nksc), d_mill(3, npwc, nksc), d_ngk(nksc))
  END SUBROUTINE dev_wfc_init
  !
  SUBROUTINE dev_wfc_put(ik, coef, mill, ngk)
    !! coefficients of k-point ik (evc layout: spin s of G ig at ig + (s-1) npwx) and the Miller index of each G
    INTEGER, INTENT(IN) :: ik, ngk, mill(3, npwc)
    COMPLEX(dp), INTENT(IN) :: coef(npwc, nbndc)
    d_coef(:, :, ik) = coef
    d_mill(:, :, ik) = mill
    d_ngk(ik) = ngk
  END SUBROUTINE dev_wfc_put
  !
  SUBROUTINE dev_wfc_free()
    INTEGER :: i, istat
    IF (ALLOCATED(d_coef)) DEALLOCATE(d_coef, d_mill, d_ngk)
    DO i = 1, SIZE(pl_batch)
      IF (pl_batch(i) > 0) istat = cufftDestroy(pl_plan(i))
    ENDDO
    pl_batch = 0
  END SUBROUTINE dev_wfc_free
  !
  INTEGER FUNCTION plan_for(nbatch)
    !! batched 3D inverse-transform plan over nbatch boxes (cached)
    INTEGER, INTENT(IN) :: nbatch
    INTEGER :: i, istat, n(3), emb(3)
    DO i = 1, SIZE(pl_batch)
      IF (pl_batch(i) == nbatch) THEN
        plan_for = pl_plan(i); RETURN
      ENDIF
    ENDDO
    i = MINLOC(pl_batch, 1)
    IF (pl_batch(i) > 0) istat = cufftDestroy(pl_plan(i))
    n = (/ nr(3), nr(2), nr(1) /)
    emb = (/ nrx(3), nrx(2), nrx(1) /)
    istat = cufftPlanMany(pl_plan(i), 3, n, emb, 1, nnr, emb, 1, nnr, CUFFT_Z2Z, nbatch)
    CALL chk(istat, 'cufftPlanMany boxes')
    pl_batch(i) = nbatch
    plan_for = pl_plan(i)
  END FUNCTION plan_for
  !
  ATTRIBUTES(GLOBAL) SUBROUTINE k_zero(u, n)
    INTEGER(8), VALUE :: n
    COMPLEX(dp) :: u(n)
    INTEGER(8) :: t
    t = INT(blockIdx%x - 1, 8) * blockDim%x + threadIdx%x
    IF (t <= n) u(t) = (0.0_dp, 0.0_dp)
  END SUBROUTINE k_zero
  !
  ATTRIBUTES(GLOBAL) SUBROUTINE k_scatter(u, coef, mill, ngk, ik, n1, nbd, s1, s2, s3, nr1, nr2, nr3, nrx1, nrx2, &
                                          nnr, npol, npw, nbnd, nks)
    ! u(box(G - shift), s, b) = coef(G + (s-1) npw/npol, n1 + b - 1, ik)
    INTEGER, VALUE :: ik, n1, nbd, s1, s2, s3, nr1, nr2, nr3, nrx1, nrx2, nnr, npol, npw, nbnd, nks
    COMPLEX(dp) :: u(nnr, npol, nbd), coef(npw, nbnd, nks)
    INTEGER :: mill(3, npw, nks), ngk(nks)
    INTEGER :: ig, b, s, pos
    ig = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    b = blockIdx%y
    IF (ig > ngk(ik)) RETURN
    pos = 1 + MODULO(mill(1, ig, ik) - s1, nr1) + nrx1 * (MODULO(mill(2, ig, ik) - s2, nr2) + &
          nrx2 * MODULO(mill(3, ig, ik) - s3, nr3))
    DO s = 1, npol
      u(pos, s, b) = coef(ig + (s - 1) * (npw / npol), n1 + b - 1, ik)
    ENDDO
  END SUBROUTINE k_scatter
  !
  SUBROUTINE dev_box(ik, n1, nbd, shift, u)
    !! u(:, :, b) = periodic part of band n1 + b - 1 at k-point ik on the box, coefficients moved by -shift
    INTEGER, INTENT(IN) :: ik, n1, nbd, shift(3)
    COMPLEX(dp), DEVICE :: u(nnr, npol, nbd)
    INTEGER :: istat
    INTEGER(8) :: ntot8
    TYPE(dim3) :: g
    ntot8 = INT(nnr, 8) * npol * nbd
    CALL k_zero<<<INT((ntot8 + TPB - 1) / TPB), TPB>>>(u, ntot8)
    g = dim3((npwc / npol + TPB - 1) / TPB, nbd, 1)
    CALL k_scatter<<<g, TPB>>>(u, d_coef, d_mill, d_ngk, ik, n1, nbd, shift(1), shift(2), shift(3), nr(1), nr(2), &
                               nr(3), nrx(1), nrx(2), nnr, npol, npwc, nbndc, nksc)
    istat = cufftExecZ2Z(plan_for(npol * nbd), u, u, CUFFT_INVERSE); CALL chk(istat, 'cufftExecZ2Z boxes')
  END SUBROUTINE dev_box
  !
  SUBROUTINE dev_add_kq(ik, jk, g0, nv, x)
    !! the transitions v (k = ik) -> c (k + q = jk + g0) from the coefficients on the device; x(c, v) as dev_add_k
    INTEGER, INTENT(IN) :: ik, jk, g0(3), nv
    REAL(dp), INTENT(IN) :: x(nc, nv)
    CALL dev_box(ik, 1, nv, (/0, 0, 0/), d_uk)
    CALL dev_box(jk, nv + 1, nc, g0, d_uq)
    CALL add_pairs(nv, x)
  END SUBROUTINE dev_add_kq
  !
  ATTRIBUTES(GLOBAL) SUBROUTINE k_imw(chi, c, sv, ng, i0, n, ncol)
    ! chi(:, i0 + w - 1) = v^1/2 [(c - c^H) / 2i] v^1/2,  c = eps^-1 at that energy
    INTEGER, VALUE :: ng, i0, n, ncol
    COMPLEX(dp) :: chi(ncol, *), c(ng, ng, n)
    REAL(dp) :: sv(ng)
    INTEGER :: i, j, w
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    j = blockIdx%y
    w = blockIdx%z
    IF (i > ng) RETURN
    chi(i + (j - 1) * ng, i0 + w - 1) = sv(i) * sv(j) * (c(i, j, w) - CONJG(c(j, i, w))) * (0.0_dp, -0.5_dp)
  END SUBROUTINE k_imw
  !
  SUBROUTINE dev_imw(sv, nw)
    !! after dev_kk: d_chi(:, w) <- Im W(w) for every energy (batched LU and inverse in chunks)
    REAL(dp), INTENT(IN) :: sv(ng)
    INTEGER, INTENT(IN) :: nw
    REAL(dp), DEVICE, ALLOCATABLE :: d_sv(:)
    COMPLEX(dp), DEVICE, ALLOCATABLE, TARGET :: d_e(:,:,:), d_c(:,:,:)
    INTEGER, DEVICE, ALLOCATABLE :: d_piv(:), d_info(:)
    TYPE(c_devptr), ALLOCATABLE :: pa(:), pc(:)
    TYPE(c_devptr), DEVICE, ALLOCATABLE :: d_pa(:), d_pc(:)
    INTEGER, ALLOCATABLE :: info(:)
    INTEGER :: w, istat, nch, i0, n
    INTEGER(8) :: freeb, totb
    TYPE(dim3) :: g
    istat = cudaMemGetInfo(freeb, totb)
    nch = INT(MAX(1_8, MIN(INT(nw, 8), (freeb / 2) / (32_8 * ng * ng + 16_8 * ng + 64))))
    ALLOCATE(d_sv(ng), d_e(ng, ng, nch), d_c(ng, ng, nch), d_piv(ng * nch), d_info(nch), pa(nch), pc(nch), &
             d_pa(nch), d_pc(nch), info(nch))
    d_sv = sv
    DO w = 1, nch
      pa(w) = c_devloc(d_e(1, 1, w)); pc(w) = c_devloc(d_c(1, 1, w))
    ENDDO
    d_pa = pa; d_pc = pc
    DO i0 = 1, nw, nch
      n = MIN(nch, nw - i0 + 1)
      g = dim3((ng + TPB - 1) / TPB, ng, n)
      CALL k_eps<<<g, TPB>>>(d_e, d_chi, d_sv, ng, i0, n, ng * ng)
      istat = cublasZgetrfBatched(hb, ng, d_pa, ng, d_piv, d_info, n); CALL chk(istat, 'getrfBatched')
      info(1:n) = d_info(1:n)
      IF (ANY(info(1:n) /= 0)) CALL chk(MAXVAL(ABS(info(1:n))), 'singular eps (getrfBatched info)')
      istat = cublasZgetriBatched(hb, ng, d_pa, ng, d_piv, d_pc, ng, d_info, n); CALL chk(istat, 'getriBatched')
      CALL k_imw<<<g, TPB>>>(d_chi, d_c, d_sv, ng, i0, n, ng * ng)
    ENDDO
    DEALLOCATE(d_sv, d_e, d_c, d_piv, d_info, pa, pc, d_pa, d_pc, info)
  END SUBROUTINE dev_imw
  !
  ATTRIBUTES(GLOBAL) SUBROUTINE k_cols(x, phi, col, ng, ne)
    ! x(:, e) = phi(:, col(e))
    INTEGER, VALUE :: ng, ne
    COMPLEX(dp) :: x(ng, ne), phi(ng, *)
    INTEGER :: col(ne)
    INTEGER :: ig, e
    INTEGER(8) :: t
    t = INT(blockIdx%x - 1, 8) * blockDim%x + threadIdx%x
    IF (t > INT(ng, 8) * ne) RETURN
    e = INT((t - 1) / ng) + 1
    ig = INT(t - INT(e - 1, 8) * ng)
    x(ig, e) = phi(ig, col(e))
  END SUBROUTINE k_cols
  !
  ATTRIBUTES(GLOBAL) SUBROUTINE k_dot(out, x, y, wt, ng, ne)
    ! out(e) = wt(e) Re sum_G conj(x(G, e)) y(G, e); one block of 128 threads per entry
    INTEGER, VALUE :: ng, ne
    REAL(dp) :: out(ne), wt(ne)
    COMPLEX(dp) :: x(ng, ne), y(ng, ne)
    REAL(dp), SHARED :: s(128)
    INTEGER :: e, i, t, h
    e = blockIdx%x
    t = threadIdx%x
    s(t) = 0.0_dp
    DO i = t, ng, 128
      s(t) = s(t) + REAL(CONJG(x(i, e)) * y(i, e), dp)
    ENDDO
    CALL syncthreads()
    h = 64
    DO WHILE (h >= 1)
      IF (t <= h) s(t) = s(t) + s(t + h)
      CALL syncthreads()
      h = h / 2
    ENDDO
    IF (t == 1) out(e) = wt(e) * s(1)
  END SUBROUTINE k_dot
  !
  SUBROUTINE dev_contract(ik, jk, g0, nv, ni, nw, ne, col, wt, cnt, off, out)
    !! after dev_imw: phi(G, b, n') = <initial b (jk + g0)| e^{i(q+G)r} |final n' (ik)> for b = 1..ni (bands nv + b)
    !! and n' = 1..nc (bands nv + n'); entries e = off(w) + 1 .. off(w) + cnt(w) interpolate at grid energy w:
    !! out(e) = wt(e) phi_col(e)^H Im W(w) phi_col(e), col(e) = (n' - 1) ni + b
    INTEGER, INTENT(IN) :: ik, jk, g0(3), nv, ni, nw, ne, col(ne), cnt(nw), off(nw)
    REAL(dp), INTENT(IN) :: wt(ne)
    REAL(dp), INTENT(OUT) :: out(ne)
    INTEGER, DEVICE, ALLOCATABLE :: d_col(:)
    REAL(dp), DEVICE, ALLOCATABLE :: d_wt(:), d_out(:)
    COMPLEX(dp), DEVICE, ALLOCATABLE :: d_x(:,:), d_y(:,:)
    INTEGER :: jn, istat, w, pl
    TYPE(dim3) :: g
    IF (.NOT. ALLOCATED(d_bf)) ALLOCATE(d_bf(nnr, npol, nc), d_bi(nnr, npol, ni), d_ff(nnr, ni), d_phi(ng, ni * nc))
    CALL dev_box(ik, nv + 1, nc, (/0, 0, 0/), d_bf)
    CALL dev_box(jk, nv + 1, ni, g0, d_bi)
    pl = plan_for(ni)
    DO jn = 1, nc
      g = dim3((nnr + TPB - 1) / TPB, ni, 1)
      CALL k_prod<<<g, TPB>>>(d_ff, d_bi, d_bf, nnr, npol, ni, nc, jn)
      istat = cufftExecZ2Z(pl, d_ff, d_ff, CUFFT_INVERSE); CALL chk(istat, 'cufftExecZ2Z contraction pairs')
      g = dim3((ng + TPB - 1) / TPB, ni, 1)
      CALL k_gather<<<g, TPB>>>(d_phi, d_ff, d_idx, ng, nnr, ni, (jn - 1) * ni, 1.0_dp / (nr(1) * nr(2) * nr(3)))
    ENDDO
    IF (ne == 0) RETURN
    ALLOCATE(d_col(ne), d_wt(ne), d_out(ne), d_x(ng, ne), d_y(ng, ne))
    d_col = col; d_wt = wt
    CALL k_cols<<<INT((INT(ng, 8) * ne + TPB - 1) / TPB), TPB>>>(d_x, d_phi, d_col, ng, ne)
    DO w = 1, nw
      IF (cnt(w) == 0) CYCLE
      istat = cublasZgemm(hb, CUBLAS_OP_N, CUBLAS_OP_N, ng, cnt(w), ng, (1.0_dp, 0.0_dp), d_chi(1, w), ng, &
                          d_x(1, off(w) + 1), ng, (0.0_dp, 0.0_dp), d_y(1, off(w) + 1), ng)
      CALL chk(istat, 'cublasZgemm contraction')
    ENDDO
    CALL k_dot<<<ne, 128>>>(d_out, d_x, d_y, d_wt, ng, ne)
    out = d_out
    DEALLOCATE(d_col, d_wt, d_out, d_x, d_y)
  END SUBROUTINE dev_contract
  !
END MODULE chi0_dev
