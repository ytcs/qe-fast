! Copyright (C) 2026 The qe-fast authors. Portions derive from Quantum ESPRESSO; see NOTICE.
! SPDX-License-Identifier: GPL-2.0-or-later
!
! chi0_dev -- device path of chi0.x (CUDA builds only). The stages of the host path, on plain arrays:
!   pairs   rho_t(G) for all conduction bands of one (k, v) by one batched cuFFT, kept in a device buffer
!   bin     buffered transitions sorted by energy bin on the host, one cuBLAS ZGEMM per occupied bin into A_b
!   kk      chi0 on the energy grid: zero-padded convolution and correlation along the bins, strided 1D cuFFT
!   inv     eps = 1 - v^1/2 chi0 v^1/2 and eps^-1 per energy, cuSOLVER LU, returned to the host in chunks
! The host path in chi0.f90 is the reference for every stage.
!
MODULE chi0_dev
  USE cudafor
  USE cufft
  USE cublas_v2
  USE cusolverDn
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 200)
  INTEGER, PARAMETER :: TPB = 256
  PUBLIC :: dev_init, dev_add_k, dev_flush, dev_kk, dev_inv_chunk, dev_free, dev_cap
  !
  INTEGER, SAVE :: nr(3), nrx(3), nnr, npol, nc, ng, nb, cap, nt, ntot
  INTEGER, SAVE :: plan3 = 0, nplan3 = 0
  TYPE(cublasHandle), SAVE :: hb
  TYPE(cusolverDnHandle), SAVE :: hs
  INTEGER, DEVICE, ALLOCATABLE, SAVE :: d_idx(:)
  COMPLEX(dp), DEVICE, ALLOCATABLE, SAVE :: d_uk(:,:,:), d_uq(:,:,:), d_f(:,:), d_r(:,:), d_a(:,:), d_chi(:,:)
  REAL(dp), ALLOCATABLE, SAVE :: h_x(:)
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
    per = 16_8 * ng_ * 3                                    ! rho, and the two weighted copies at flush
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
    istat = cusolverDnCreate(hs); CALL chk(istat, 'cusolverDnCreate')
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
    INTEGER :: iv, istat
    TYPE(dim3) :: g
    d_uk = uk
    d_uq = uq
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
  END SUBROUTINE dev_add_k
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
    !! eps and eps^-1 for energies i0 .. i0+n-1, returned to the host (eps_h only its full matrices, for the head)
    REAL(dp), INTENT(IN) :: sv(ng)
    INTEGER, INTENT(IN) :: i0, n
    COMPLEX(dp), INTENT(OUT) :: eps_h(ng, ng, n), einv_h(ng, ng, n)
    REAL(dp), DEVICE, ALLOCATABLE :: d_sv(:)
    COMPLEX(dp), DEVICE, ALLOCATABLE :: d_e(:,:,:), d_b(:,:), d_w(:)
    INTEGER, DEVICE, ALLOCATABLE :: d_piv(:), d_info(:)
    COMPLEX(dp), ALLOCATABLE :: ident(:,:)
    INTEGER :: w, lw, istat, i
    TYPE(dim3) :: g
    ALLOCATE(d_sv(ng), d_e(ng, ng, n), d_b(ng, ng), d_piv(ng), d_info(1), ident(ng, ng))
    d_sv = sv
    g = dim3((ng + TPB - 1) / TPB, ng, n)
    CALL k_eps<<<g, TPB>>>(d_e, d_chi, d_sv, ng, i0, n, ng * ng)
    eps_h = d_e
    istat = cusolverDnZgetrf_bufferSize(hs, ng, ng, d_e(:, :, 1), ng, lw); CALL chk(istat, 'getrf_bufferSize')
    ALLOCATE(d_w(MAX(lw, 1)))
    ident = (0.0_dp, 0.0_dp)
    DO i = 1, ng
      ident(i, i) = (1.0_dp, 0.0_dp)
    ENDDO
    DO w = 1, n
      istat = cusolverDnZgetrf(hs, ng, ng, d_e(:, :, w), ng, d_w, d_piv, d_info); CALL chk(istat, 'getrf')
      d_b = ident
      istat = cusolverDnZgetrs(hs, CUBLAS_OP_N, ng, ng, d_e(:, :, w), ng, d_piv, d_b, ng, d_info)
      CALL chk(istat, 'getrs')
      einv_h(:, :, w) = d_b
    ENDDO
    DEALLOCATE(d_sv, d_e, d_b, d_piv, d_info, d_w, ident)
  END SUBROUTINE dev_inv_chunk
  !
  SUBROUTINE dev_free()
    INTEGER :: istat
    IF (ALLOCATED(d_chi)) DEALLOCATE(d_chi)
    IF (ALLOCATED(d_a)) DEALLOCATE(d_a)
    IF (ALLOCATED(d_uk)) DEALLOCATE(d_uk, d_uq, d_f, d_r, d_idx, h_x)
    istat = cufftDestroy(plan3)
    istat = cublasDestroy(hb)
    istat = cusolverDnDestroy(hs)
  END SUBROUTINE dev_free
  !
END MODULE chi0_dev
