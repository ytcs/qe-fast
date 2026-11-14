! Copyright (C) 2026 The qe-fast authors. Portions derive from Quantum ESPRESSO; see NOTICE.
! SPDX-License-Identifier: GPL-2.0-or-later
!
! gw_common -- the host stages shared by chi0.x and gamma_ii.x (issue #1): k + q on the grid, the local-field
! G basis of one q, periodic parts on the FFT box, pair densities, binning in energy, Kramers-Kronig by FFT,
! eps^-1 per energy. Serial FFT box, one MPI process.
!
MODULE gw_common
  USE kinds,     ONLY : DP
  USE constants, ONLY : rytoev
  IMPLICIT NONE
  PRIVATE
  REAL(DP), PARAMETER, PUBLIC :: ha = 2.0_DP * rytoev            ! eV per Hartree
  PUBLIC :: kq_grid, gbasis_q, box_index, wfc_r, pair_row, bin_add, kk_fft, eps_invert, mem_avail_gb
  !
CONTAINS
  !
  SUBROUTINE kq_grid(kcr, qcr, kq, g0)
    !! k + q = k_kq + g0 (crystal) for every k of the grid kcr(3, nk)
    REAL(DP), INTENT(IN) :: kcr(:,:), qcr(3)
    INTEGER, INTENT(OUT) :: kq(:), g0(:,:)
    INTEGER :: ik, jk
    DO ik = 1, SIZE(kcr, 2)
      kq(ik) = 0
      DO jk = 1, SIZE(kcr, 2)
        IF (MAXVAL(ABS(kcr(:, ik) + qcr - kcr(:, jk) - NINT(kcr(:, ik) + qcr - kcr(:, jk)))) < 1.0d-6) THEN
          kq(ik) = jk; g0(:, ik) = NINT(kcr(:, ik) + qcr - kcr(:, jk)); EXIT
        ENDIF
      ENDDO
      IF (kq(ik) == 0) CALL errore('gw_common', 'k + q is not on the k grid; choose q commensurate with the grid', ik)
    ENDDO
  END SUBROUTINE kq_grid
  !
  SUBROUTINE gbasis_q(qcr, must, ecut, ng, gm, q2, idx_mg, ihead)
    !! local-field basis of q (crystal): G with |q+G|^2 <= ecut (Ry = bohr^-2), plus the G 'must'; sorted by |q+G|^2
    !! (stable, rounded to 1e-10 bohr^-2: ties in enumeration order, as chi0_ref.py); idx_mg = where fwfft puts -G;
    !! ihead = position of 'must'
    USE cell_base, ONLY : bg, tpiba2
    USE gvect,     ONLY : mill, ngm
    USE fft_base,  ONLY : dffts
    REAL(DP), INTENT(IN) :: qcr(3), ecut
    INTEGER, INTENT(IN) :: must(3)
    INTEGER, INTENT(OUT) :: ng, ihead
    INTEGER, ALLOCATABLE, INTENT(OUT) :: gm(:,:), idx_mg(:)
    REAL(DP), ALLOCATABLE, INTENT(OUT) :: q2(:)
    INTEGER :: n, i1, i2, i3, cnt, k, l, ti(3)
    REAL(DP) :: v(3), t
    INTEGER, ALLOCATABLE :: tmp(:,:)
    REAL(DP), ALLOCATABLE :: tq(:)
    n = CEILING(SQRT(ecut / tpiba2) / MINVAL(SQRT(SUM(bg**2, 1)))) + 2
    ALLOCATE(tmp(3, (2 * n + 1)**3), tq((2 * n + 1)**3))
    cnt = 0
    DO i1 = -n, n
      DO i2 = -n, n
        DO i3 = -n, n
          v = (qcr(1) + i1) * bg(:, 1) + (qcr(2) + i2) * bg(:, 2) + (qcr(3) + i3) * bg(:, 3)
          t = SUM(v**2) * tpiba2
          IF (t <= ecut .OR. ALL((/i1, i2, i3/) == must)) THEN
            cnt = cnt + 1; tmp(:, cnt) = (/i1, i2, i3/); tq(cnt) = t
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    DO k = 2, cnt
      t = tq(k); ti = tmp(:, k); l = k - 1
      DO WHILE (l >= 1)
        IF (ANINT(tq(l) * 1.0d10) <= ANINT(t * 1.0d10)) EXIT
        tq(l + 1) = tq(l); tmp(:, l + 1) = tmp(:, l); l = l - 1
      ENDDO
      tq(l + 1) = t; tmp(:, l + 1) = ti
    ENDDO
    ng = cnt
    ALLOCATE(gm(3, cnt), q2(cnt), idx_mg(cnt))
    gm = tmp(:, 1:cnt); q2 = tq(1:cnt)
    ihead = 0
    DO k = 1, cnt
      ! where fwfft puts -G: through QE's own map (dffts%nl of the G-vector index), not a box formula, so that it
      ! holds for every FFT layout QE builds
      idx_mg(k) = 0
      DO l = 1, ngm
        IF (ALL(mill(:, l) == -gm(:, k))) THEN
          idx_mg(k) = dffts%nl(l); EXIT
        ENDIF
      ENDDO
      IF (idx_mg(k) == 0) CALL errore('gw_common', 'local-field G outside the density G sphere', k)
      IF (ALL(gm(:, k) == must)) ihead = k
    ENDDO
    DEALLOCATE(tmp, tq)
  END SUBROUTINE gbasis_q
  !
  INTEGER FUNCTION box_index(m)
    !! position of the Miller index m on the (serial) smooth FFT box, the layout of dffts%nl
    USE fft_base, ONLY : dffts
    INTEGER, INTENT(IN) :: m(3)
    box_index = 1 + MODULO(m(1), dffts%nr1) + dffts%nr1x * (MODULO(m(2), dffts%nr2) + &
                dffts%nr2x * MODULO(m(3), dffts%nr3))
  END FUNCTION box_index
  !
  SUBROUTINE wfc_r(ikk, n1, n2, shift, u)
    !! u(:, :, n) = periodic part of bands n1..n2 at k-point ikk on the FFT box, coefficients moved by -shift
    USE io_files,       ONLY : restart_dir
    USE pw_restart_new, ONLY : read_collected_wfc
    USE wavefunctions,  ONLY : evc
    USE wvfct,          ONLY : npwx
    USE klist,          ONLY : ngk, igk_k
    USE gvect,          ONLY : mill
    USE noncollin_module, ONLY : npol
    USE fft_base,       ONLY : dffts
    USE fft_interfaces, ONLY : invfft
    INTEGER, INTENT(IN) :: ikk, n1, n2, shift(3)
    COMPLEX(DP), INTENT(OUT) :: u(:,:,:)
    INTEGER :: n, igg, s, pos
    CALL read_collected_wfc(restart_dir(), ikk, evc)
    DO n = n1, n2
      DO s = 1, npol
        u(:, s, n - n1 + 1) = (0.0_DP, 0.0_DP)
        DO igg = 1, ngk(ikk)
          pos = box_index(mill(:, igk_k(igg, ikk)) - shift)
          u(pos, s, n - n1 + 1) = evc(igg + (s - 1) * npwx, n)
        ENDDO
        ! full-box transform: 'Wave' transforms only the columns of the wavefunction sphere, and the coefficients
        ! moved by -shift (k + q outside the first zone) can fall outside them
        CALL invfft('Rho', u(:, s, n - n1 + 1), dffts)
      ENDDO
    ENDDO
  END SUBROUTINE wfc_r
  !
  SUBROUTINE pair_row(uket, ubra, nbra, idx_mg, rho, f)
    !! rho(G, b) = <bra_b| e^{i(q+G)r} |ket> for one ket and bras 1..nbra on the box (bra at k + q, moved into the
    !! zone of k); f is a work array of the box size
    USE noncollin_module, ONLY : npol
    USE fft_base,       ONLY : dffts
    USE fft_interfaces, ONLY : fwfft
    INTEGER, INTENT(IN) :: nbra, idx_mg(:)
    COMPLEX(DP), INTENT(IN) :: uket(:,:), ubra(:,:,:)
    COMPLEX(DP), INTENT(OUT) :: rho(:,:), f(:)
    INTEGER :: ib, ig
    DO ib = 1, nbra
      f(:) = CONJG(ubra(:, 1, ib)) * uket(:, 1)
      IF (npol == 2) f(:) = f(:) + CONJG(ubra(:, 2, ib)) * uket(:, 2)
      CALL fwfft('Rho', f, dffts)                     ! (1/N) sum_r f e^{-iGr}: rho(G) sits at -G
      DO ig = 1, SIZE(idx_mg)
        rho(ig, ib) = f(idx_mg(ig))
      ENDDO
    ENDDO
  END SUBROUTINE pair_row
  !
  SUBROUTINE bin_add(a, rho, n, x)
    !! a(:, :, b) += sum_t w_tb rho_t^* rho_t^T over the transitions t = 1..n of energy x_t (in units of the bin
    !! width), linear binning into the bins floor(x) + 1 and floor(x) + 2; transitions past the last bin are dropped
    COMPLEX(DP), INTENT(INOUT) :: a(:,:,:)
    COMPLEX(DP), INTENT(IN) :: rho(:,:)
    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN) :: x(:)
    INTEGER :: ng, nb, b, i, j, t, np
    INTEGER, ALLOCATABLE :: tb(:), list(:)
    REAL(DP), ALLOCATABLE :: tw(:)
    COMPLEX(DP), ALLOCATABLE :: xm(:,:), ym(:,:)
    ng = SIZE(a, 1); nb = SIZE(a, 3)
    ALLOCATE(tb(n), list(n), tw(n), xm(ng, n), ym(ng, n))
    DO b = 1, 2
      DO t = 1, n
        tb(t) = FLOOR(x(t)) + 1 + (b - 1)
        IF (b == 1) THEN
          tw(t) = 1.0_DP - (x(t) - FLOOR(x(t)))
        ELSE
          tw(t) = x(t) - FLOOR(x(t))
        ENDIF
      ENDDO
      DO i = 1, nb
        np = 0
        DO t = 1, n
          IF (tb(t) == i .AND. FLOOR(x(t)) + 2 <= nb) THEN
            np = np + 1; list(np) = t
          ENDIF
        ENDDO
        IF (np == 0) CYCLE
        DO j = 1, np
          xm(:, j) = CONJG(rho(:, list(j))) * tw(list(j))
          ym(:, j) = rho(:, list(j))
        ENDDO
        CALL ZGEMM('N', 'T', ng, ng, np, (1.0_DP, 0.0_DP), xm, ng, ym, ng, (1.0_DP, 0.0_DP), a(1, 1, i), ng)
      ENDDO
    ENDDO
    DEALLOCATE(tb, list, tw, xm, ym)
  END SUBROUTINE bin_add
  !
  SUBROUTINE kk_fft(af, cf, ncol, nb, nw, de, eta, pref)
    !! cf(:, i) = pref sum_b af(:, b) [K(i - b) - K(i + b)],  K(m) = 1/(m de + i eta) (de, eta in eV; K in 1/Hartree),
    !! on the uniform grid w_i = (i-1) de: a zero-padded convolution and a correlation along b, by QE's 1D FFTs
    USE fft_scalar, ONLY : cft_1z
    INTEGER, INTENT(IN) :: ncol, nb, nw
    REAL(DP), INTENT(IN) :: de, eta, pref
    COMPLEX(DP), INTENT(IN) :: af(ncol, nb)
    COMPLEX(DP), INTENT(OUT) :: cf(ncol, nw)
    INTEGER :: l, m, c0, nbl, c, jj, mlast
    COMPLEX(DP), ALLOCATABLE :: k1(:), k2(:), k1f(:), k2f(:), buf(:), bf(:), d(:), df(:)
    REAL(DP) :: sf, sb, sf1, scal
    l = 1
    DO WHILE (l < 2 * nb + nw)
      l = 2 * l
    ENDDO
    ALLOCATE(k1(l), k2(l), k1f(l), k2f(l), d(l), df(l))
    k1 = (0.0_DP, 0.0_DP); k2 = (0.0_DP, 0.0_DP)
    DO jj = 1, nw + nb - 1
      m = (jj - 1) - (nb - 1)                         ! w_i - D_b = (i - b) de
      k1(jj) = 1.0_DP / CMPLX(m * de / ha, eta / ha, KIND=DP)
      m = jj - 1                                      ! w_i + D_b = (i + b) de
      k2(jj) = 1.0_DP / CMPLX(m * de / ha, eta / ha, KIND=DP)
    ENDDO
    CALL cft_1z(k1, 1, l, l, -1, k1f)
    CALL cft_1z(k2, 1, l, l, -1, k2f)
    d = (0.0_DP, 0.0_DP); d(1) = (1.0_DP, 0.0_DP)
    CALL cft_1z(d, 1, l, l, -1, df); sf1 = REAL(df(1))      ! forward scaling of the single (kernel) transforms
    nbl = MAX(1, MIN(ncol, INT(2.0d6 / l)))          ! two block buffers of <= ~32 MB each (column c at (c-1) l)
    ALLOCATE(buf(l * nbl), bf(l * nbl))
    mlast = -1
    DO c0 = 1, ncol, nbl
      m = MIN(nbl, ncol - c0 + 1)
      IF (m /= mlast) THEN
        ! the scalings of this FFT backend, measured with the batch size of the block: QE's internal backend (GPU
        ! builds) scales a batched call differently from a single transform, and applies the forward scaling to
        ! its input array in place (reset before each call)
        buf = (0.0_DP, 0.0_DP); buf(1) = (1.0_DP, 0.0_DP)
        CALL cft_1z(buf, m, l, l, -1, bf); sf = REAL(bf(1))
        buf = (0.0_DP, 0.0_DP); buf(1) = (1.0_DP, 0.0_DP)
        CALL cft_1z(buf, m, l, l, +1, bf); sb = REAL(bf(1))
        scal = pref / (sb * sf * sf1 * l)
        mlast = m
      ENDIF
      ! convolution term
      buf = (0.0_DP, 0.0_DP)
      DO c = 1, m
        buf((c - 1) * l + 1:(c - 1) * l + nb) = af(c0 + c - 1, 1:nb)
      ENDDO
      CALL cft_1z(buf, m, l, l, -1, bf)
      DO c = 1, m
        bf((c - 1) * l + 1:c * l) = bf((c - 1) * l + 1:c * l) * k1f(:)
      ENDDO
      CALL cft_1z(bf, m, l, l, +1, buf)
      DO c = 1, m
        cf(c0 + c - 1, 1:nw) = buf((c - 1) * l + nb:(c - 1) * l + nb + nw - 1)
      ENDDO
      ! correlation term (A reversed in b)
      buf = (0.0_DP, 0.0_DP)
      DO c = 1, m
        buf((c - 1) * l + 1:(c - 1) * l + nb) = af(c0 + c - 1, nb:1:-1)
      ENDDO
      CALL cft_1z(buf, m, l, l, -1, bf)
      DO c = 1, m
        bf((c - 1) * l + 1:c * l) = bf((c - 1) * l + 1:c * l) * k2f(:)
      ENDDO
      CALL cft_1z(bf, m, l, l, +1, buf)
      DO c = 1, m
        cf(c0 + c - 1, 1:nw) = (cf(c0 + c - 1, 1:nw) - buf((c - 1) * l + nb:(c - 1) * l + nb + nw - 1)) * scal
      ENDDO
    ENDDO
    DEALLOCATE(k1, k2, k1f, k2f, d, df, buf, bf)
  END SUBROUTINE kk_fft
  !
  SUBROUTINE eps_invert(sv, chiw, epsm, ehead, ihead)
    !! epsm = eps^-1 at one energy, eps = 1 - v^1/2 chi0 v^1/2; ehead = eps(ihead, ihead) (0 if ihead = 0)
    REAL(DP), INTENT(IN) :: sv(:)
    COMPLEX(DP), INTENT(IN) :: chiw(:,:)
    COMPLEX(DP), INTENT(OUT) :: epsm(:,:), ehead
    INTEGER, INTENT(IN) :: ihead
    INTEGER :: ng, i, j, info
    INTEGER, ALLOCATABLE :: ipiv(:)
    COMPLEX(DP), ALLOCATABLE :: work(:)
    ng = SIZE(sv)
    ALLOCATE(ipiv(ng), work(64 * ng))
    DO j = 1, ng
      DO i = 1, ng
        epsm(i, j) = - sv(i) * chiw(i, j) * sv(j)
      ENDDO
      epsm(j, j) = epsm(j, j) + 1.0_DP
    ENDDO
    ehead = (0.0_DP, 0.0_DP)
    IF (ihead > 0) ehead = epsm(ihead, ihead)
    CALL ZGETRF(ng, ng, epsm, ng, ipiv, info)
    IF (info /= 0) CALL errore('gw_common', 'singular eps', info)
    CALL ZGETRI(ng, epsm, ng, ipiv, work, SIZE(work), info)
    DEALLOCATE(ipiv, work)
  END SUBROUTINE eps_invert
  !
  REAL(DP) FUNCTION mem_avail_gb()
    !! MemAvailable from /proc/meminfo (GB); 2 GB if it cannot be read
    CHARACTER(LEN=128) :: line
    INTEGER :: u, ios2
    REAL(DP) :: kb
    mem_avail_gb = 2.0_DP
    OPEN(NEWUNIT=u, FILE='/proc/meminfo', STATUS='old', ACTION='read', IOSTAT=ios2)
    IF (ios2 /= 0) RETURN
    DO
      READ(u, '(a)', IOSTAT=ios2) line
      IF (ios2 /= 0) EXIT
      IF (line(1:13) == 'MemAvailable:') THEN
        READ(line(14:), *, IOSTAT=ios2) kb
        IF (ios2 == 0) mem_avail_gb = kb / 1.0d6
        EXIT
      ENDIF
    ENDDO
    CLOSE(u)
  END FUNCTION mem_avail_gb
  !
END MODULE gw_common
