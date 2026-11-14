! Copyright (C) 2026 The qe-fast authors. Portions derive from Quantum ESPRESSO; see NOTICE.
! SPDX-License-Identifier: GPL-2.0-or-later
!
! chi0.x -- inverse dielectric matrix eps^-1_GG'(q, w) by sum over bands (issue #1). Host path; the stages are
! those of qefast/gw/chi0_ref.py, which is the reference this program is checked against:
!
!   rho_t(G)    = <c, k+q| e^{i(q+G)r} |v, k>                   pair densities (QE FFTs)
!   A_b(G,G')   = sum_t w_tb rho_t(G)^* rho_t(G')                 transitions binned in energy (one ZGEMM per bin)
!   chi0(w)     = spin/(Nk Omega) sum_b A_b [1/(w - D_b + i eta) - 1/(w + D_b + i eta)]   (time reversal),
!                 on the uniform grid a convolution and a correlation in b, done by FFT
!   eps         = 1 - v^1/2 chi0 v^1/2,  v = 4 pi / |q+G|^2;  eps^-1 by LU per energy
!
! Supported: norm-conserving, fixed occupations, no spin polarisation or magnetisation, q /= 0, a pw.x nscf on a
! full k grid without symmetry (k + q on the grid), one MPI process. Hartree units inside, eV in the input/output.
!
! Input namelist &inputchi0: prefix, outdir, q(3) (2 pi / alat, Cartesian), ecut_eps (Ry), emax, de, eta (eV),
! nbnd (0: all), fileout, filfull ('' : no full matrix), max_mem (GB).
!
PROGRAM chi0
  USE kinds,            ONLY : DP
  USE constants,        ONLY : rytoev, pi, tpi
  USE io_files,         ONLY : prefix, tmp_dir, restart_dir
  USE io_global,        ONLY : ionode, ionode_id, stdout
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : mp_startup
  USE mp_images,        ONLY : intra_image_comm
  USE mp_pools,         ONLY : npool, nproc_pool
  USE environment,      ONLY : environment_start, environment_end
  USE wvfct,            ONLY : nbnd_tot => nbnd, npwx, et
  USE klist,            ONLY : xk, nks, ngk, igk_k, nelec, lgauss, ltetra
  USE wavefunctions,    ONLY : evc
  USE gvect,            ONLY : mill, ngm
  USE cell_base,        ONLY : at, bg, tpiba2, omega
  USE fft_base,         ONLY : dffts
  USE fft_interfaces,   ONLY : invfft, fwfft
  USE noncollin_module, ONLY : npol, noncolin, domag
  USE lsda_mod,         ONLY : lsda
  USE uspp,             ONLY : okvan
  USE control_flags,    ONLY : gamma_only
  USE pw_restart_new,   ONLY : read_collected_wfc
  USE qef_config,       ONLY : qef_note
#if defined(__CUDA)
  USE chi0_dev,         ONLY : dev_cap
#endif
  IMPLICIT NONE
  REAL(DP), PARAMETER :: ha = 2.0_DP * rytoev            ! eV per Hartree
  CHARACTER(LEN=256) :: outdir, fileout, filfull
  REAL(DP) :: q(3), ecut_eps, emax, de, eta, max_mem
  INTEGER :: nbnd, ios
  LOGICAL :: device
  CHARACTER(LEN=8) :: dbg
  NAMELIST / inputchi0 / prefix, outdir, q, ecut_eps, emax, de, eta, nbnd, fileout, filfull, max_mem, device
  !
  INTEGER :: nv, nc, ik, jk, iv, ic, ig, i, j, nb, nw, wc, i0, iw, nwc, ng_max, nt, b, iq_head, np, cap
  INTEGER, ALLOCATABLE :: kq(:), g0(:,:), gm(:,:), idx_mg(:), tb(:), list(:)
  REAL(DP) :: qcr(3), qc(3), x, spin, need, ahmin, ahmax
  REAL(DP), ALLOCATABLE :: kcr(:,:), q2(:), sv(:), wgrid(:), tw(:)
  COMPLEX(DP), ALLOCATABLE :: a(:,:,:), uk(:,:,:), uq(:,:,:), f(:), rho(:,:), xm(:,:), ym(:,:), &
                              chi(:,:,:), epsm(:,:), work(:), eqq(:), einvqq(:), ah(:,:)
  INTEGER, ALLOCATABLE :: ipiv(:)
  REAL(DP), ALLOCATABLE :: rwork(:), evals(:)
  LOGICAL :: needwf = .TRUE.
  !
  CALL mp_startup()
  CALL environment_start('CHI0')
  prefix = 'pwscf'; outdir = './'; q = 0.0_DP; ecut_eps = 4.0_DP; emax = 60.0_DP; de = 0.02_DP; eta = 0.1_DP
  nbnd = 0; fileout = 'chi0.dat'; filfull = ' '; max_mem = 2.0_DP; device = .TRUE.
  IF (ionode) THEN
    CALL input_from_file()
    READ(5, inputchi0, IOSTAT=ios)
  ENDIF
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  IF (ios /= 0) CALL errore('chi0', 'reading namelist inputchi0', ABS(ios))
  CALL mp_bcast(prefix, ionode_id, intra_image_comm); CALL mp_bcast(outdir, ionode_id, intra_image_comm)
  CALL mp_bcast(q, ionode_id, intra_image_comm); CALL mp_bcast(ecut_eps, ionode_id, intra_image_comm)
  CALL mp_bcast(emax, ionode_id, intra_image_comm); CALL mp_bcast(de, ionode_id, intra_image_comm)
  CALL mp_bcast(eta, ionode_id, intra_image_comm); CALL mp_bcast(nbnd, ionode_id, intra_image_comm)
  CALL mp_bcast(fileout, ionode_id, intra_image_comm); CALL mp_bcast(filfull, ionode_id, intra_image_comm)
  CALL mp_bcast(max_mem, ionode_id, intra_image_comm); CALL mp_bcast(device, ionode_id, intra_image_comm)
#if !defined(__CUDA)
  device = .FALSE.
#endif
  tmp_dir = TRIM(outdir) // '/'
  CALL read_file_new(needwf)
  !
  ! ---- the supported case
  IF (npool * nproc_pool > 1) CALL errore('chi0', 'host path: run with one MPI process for now', 1)
  IF (okvan) CALL errore('chi0', 'ultrasoft or PAW pseudopotentials are not supported (norm-conserving only)', 1)
  IF (lgauss .OR. ltetra) CALL errore('chi0', 'only fixed occupations (insulators) are supported', 1)
  IF (lsda) CALL errore('chi0', 'spin-polarised (lsda) runs are not supported', 1)
  IF (noncolin .AND. domag) CALL errore('chi0', 'magnetic runs are not supported (time reversal is assumed)', 1)
  IF (gamma_only) CALL errore('chi0', 'gamma_only runs are not supported', 1)
  IF (dffts%nnr /= dffts%nr1x * dffts%nr2x * dffts%nr3x) CALL errore('chi0', 'distributed FFT', 1)
  IF (nbnd <= 0 .OR. nbnd > nbnd_tot) nbnd = nbnd_tot
  IF (noncolin) THEN
    nv = NINT(nelec); spin = 1.0_DP
  ELSE
    nv = NINT(nelec / 2.0_DP); spin = 2.0_DP
  ENDIF
  nc = nbnd - nv
  IF (nc < 1) CALL errore('chi0', 'no empty bands', 1)
  !
  ! ---- q, k + q on the grid
  ALLOCATE(kcr(3, nks), kq(nks), g0(3, nks))
  DO ik = 1, nks
    kcr(:, ik) = MATMUL(xk(:, ik), at)                 ! crystal coordinates (at: columns a_i, alat units)
  ENDDO
  qcr = MATMUL(q, at)
  qc = CEILING(qcr - 0.5_DP - 1.0d-9)                 ! G_Q, crystal
  qcr = qcr - qc
  IF (MAXVAL(ABS(qcr)) < 1.0d-8) CALL errore('chi0', 'q = 0 needs the head and wings treatment; not supported', 1)
  DO ik = 1, nks
    kq(ik) = 0
    DO jk = 1, nks
      IF (MAXVAL(ABS(kcr(:, ik) + qcr - kcr(:, jk) - NINT(kcr(:, ik) + qcr - kcr(:, jk)))) < 1.0d-6) THEN
        kq(ik) = jk; g0(:, ik) = NINT(kcr(:, ik) + qcr - kcr(:, jk)); EXIT
      ENDIF
    ENDDO
    IF (kq(ik) == 0) CALL errore('chi0', 'k + q is not on the k grid; choose q commensurate with the grid', ik)
  ENDDO
  !
  ! ---- local-field basis |q+G|^2 <= ecut_eps (Ry = bohr^-2), sorted by |q+G| (stable), G_Q always in
  CALL gbasis()
  ALLOCATE(sv(ng_max))
  sv = SQRT(4.0_DP * pi / q2)
  WRITE(stdout, '(/5x,"chi0: ",I0," k, ",I0," bands (",I0," occupied), FFT ",3I4,", ",I0," G (|q+G|^2 <= ", &
       & F6.2," Ry), q crystal ",3F9.4," + G ",3I3)') nks, nbnd, nv, dffts%nr1, dffts%nr2, dffts%nr3, ng_max, &
       ecut_eps, qcr, NINT(qc)
  !
  ! ---- memory: binned matrices + one w chunk (chi, eps, work); shrink the chunk, stop only if A alone is too big
  nb = CEILING(emax / de) + 2
  nw = 0
  DO WHILE (nw * de < emax - 1.0d-12)
    nw = nw + 1
  ENDDO
  wc = nw
  need = (16.0_DP * (nb + nw) * ng_max * ng_max + 6.4d7) / 1.0d9   ! A, chi on the grid, FFT block buffers
  WRITE(stdout, '(5x,"chi0: memory estimate ",F6.2," GB (",I0," bins and ",I0," energies x ",I0,"^2)")') &
       need, nb, nw, ng_max
  IF (need > max_mem) CALL errore('chi0', 'memory estimate exceeds max_mem; raise de, lower emax or ecut_eps', 1)
  !
  ALLOCATE(wgrid(nw), eqq(nw), einvqq(nw), ah(ng_max, ng_max), rwork(3 * ng_max), evals(ng_max), &
           work(64 * ng_max), epsm(ng_max, ng_max), ipiv(ng_max))
  DO iw = 1, nw
    wgrid(iw) = (iw - 1) * de
  ENDDO
#if defined(__CUDA)
  IF (device) THEN
    cap = dev_cap(dffts%nnr, npol, nv, nc, ng_max, nb, nw, 2.0_DP)
    IF (cap < nc) THEN
      CALL qef_note('chi0 on the GPU', .FALSE., 'device memory too small for one k-point; host path')
      device = .FALSE.
    ENDIF
  ENDIF
  IF (device) THEN
    CALL qef_note('chi0 on the GPU', .TRUE.)
    CALL device_path()
  ENDIF
#endif
  IF (.NOT. device) THEN
  ! ---- A_b: pair densities and binning
  ALLOCATE(a(ng_max, ng_max, nb), uk(dffts%nnr, npol, nv), uq(dffts%nnr, npol, nc), f(dffts%nnr), &
           rho(ng_max, nc), xm(ng_max, nc), ym(ng_max, nc), tb(nc), tw(nc), list(nc))
  a = (0.0_DP, 0.0_DP)
  DO ik = 1, nks
    jk = kq(ik)
    CALL start_clock('chi0_wfc')
    CALL wfc_r(ik, 1, nv, (/0, 0, 0/), uk)
    CALL wfc_r(jk, nv + 1, nbnd, g0(:, ik), uq)
    CALL stop_clock('chi0_wfc')
    DO iv = 1, nv
      CALL start_clock('chi0_pairs')
      DO ic = 1, nc
        f(:) = CONJG(uq(:, 1, ic)) * uk(:, 1, iv)
        IF (npol == 2) f(:) = f(:) + CONJG(uq(:, 2, ic)) * uk(:, 2, iv)
        CALL fwfft('Rho', f, dffts)                     ! (1/N) sum_r f e^{-iGr}: rho(G) sits at -G
        DO ig = 1, ng_max
          rho(ig, ic) = f(idx_mg(ig))
        ENDDO
      ENDDO
      CALL stop_clock('chi0_pairs')
      CALL start_clock('chi0_bin')
      ! linear binning of D = e_c(k+q) - e_v(k) into bins b and b+1
      DO b = 1, 2
        DO ic = 1, nc
          x = (et(nv + ic, jk) - et(iv, ik)) * 0.5_DP * ha / de
          tb(ic) = FLOOR(x) + 1 + (b - 1)                   ! 1-based bin
          IF (b == 1) THEN
            tw(ic) = 1.0_DP - (x - FLOOR(x))
          ELSE
            tw(ic) = x - FLOOR(x)
          ENDIF
        ENDDO
        DO i = 1, nb
          np = 0
          DO ic = 1, nc
            IF (tb(ic) == i .AND. FLOOR((et(nv + ic, jk) - et(iv, ik)) * 0.5_DP * ha / de) + 2 <= nb) THEN
              np = np + 1; list(np) = ic
            ENDIF
          ENDDO
          IF (np == 0) CYCLE
          DO j = 1, np
            xm(:, j) = CONJG(rho(:, list(j))) * tw(list(j))
            ym(:, j) = rho(:, list(j))
          ENDDO
          CALL ZGEMM('N', 'T', ng_max, ng_max, np, (1.0_DP, 0.0_DP), xm, ng_max, ym, ng_max, (1.0_DP, 0.0_DP), &
                     a(1, 1, i), ng_max)
        ENDDO
      ENDDO
      CALL stop_clock('chi0_bin')
    ENDDO
  ENDDO
  DEALLOCATE(uk, uq, f, rho, xm, ym, tb, tw, list)
  !
  ! ---- chi0 on the energy grid by FFT (Kramers-Kronig as a convolution and a correlation), then eps^-1 per energy
  ALLOCATE(chi(ng_max, ng_max, nw))
  CALL get_environment_variable('CHI0_DEBUG', dbg)
  IF (LEN_TRIM(dbg) > 0) WRITE(stdout, '(5x,"chi0 debug host: sum|A| ",ES14.6,", A(1,1,b) sum ",2ES14.6)') &
       SUM(ABS(a)), SUM(a(1, 1, :))
  CALL start_clock('chi0_kk')
  CALL kk_fft(a, chi, ng_max * ng_max)
  CALL stop_clock('chi0_kk')
  IF (LEN_TRIM(dbg) > 0) WRITE(stdout, '(5x,"chi0 debug host: sum|chi| ",ES14.6,", chi(1,1,1) ",2ES14.6)') &
       SUM(ABS(chi)), chi(1, 1, 1)
  DEALLOCATE(a)
  CALL open_full()
  CALL start_clock('chi0_inv')
  DO iw = 1, nw
    DO j = 1, ng_max
      DO i = 1, ng_max
        epsm(i, j) = - sv(i) * chi(i, j, iw) * sv(j)
      ENDDO
      epsm(j, j) = epsm(j, j) + 1.0_DP
    ENDDO
    eqq(iw) = epsm(iq_head, iq_head)
    CALL ZGETRF(ng_max, ng_max, epsm, ng_max, ipiv, ios)
    IF (ios /= 0) CALL errore('chi0', 'singular eps', ios)
    CALL ZGETRI(ng_max, epsm, ng_max, ipiv, work, SIZE(work), ios)
    CALL finish_w(iw, epsm)
  ENDDO
  CALL stop_clock('chi0_inv')
  ENDIF   ! host path
  IF (LEN_TRIM(filfull) > 0 .AND. ionode) CLOSE(91)
  !
  ! ---- output (columns as chi0_ref.py)
  IF (ionode) THEN
    OPEN(90, FILE=TRIM(fileout), STATUS='replace')
    WRITE(90, '(a)') '# w(eV)  Im eps_NLF(=turbo IPA)  L_NLF  L_RPA  Re eps_M(RPA)  Im eps_M(RPA)'
    DO iw = 1, nw
      WRITE(90, '(6ES20.10)') wgrid(iw), AIMAG(eqq(iw)), -AIMAG(1.0_DP / eqq(iw)), -AIMAG(einvqq(iw)), &
           REAL(1.0_DP / einvqq(iw)), AIMAG(1.0_DP / einvqq(iw))
    ENDDO
    CLOSE(90)
  ENDIF
  WRITE(stdout, '(5x,"chi0: anti-Hermitian part of -eps^-1 over w > 0: eigenvalues in [",ES10.2,",",ES10.2,"]")') &
       ahmin, ahmax
  WRITE(stdout, '(5x,"chi0: RPA loss peak ",F8.2," eV, Re eps_M(w=0) RPA ",F8.3,", no local fields ",F8.3)') &
       wgrid(MAXLOC(-AIMAG(einvqq), 1)), REAL(1.0_DP / einvqq(1)), REAL(eqq(1))
  CALL print_clock('chi0_wfc'); CALL print_clock('chi0_pairs'); CALL print_clock('chi0_bin')
  CALL print_clock('chi0_kk'); CALL print_clock('chi0_inv')
  CALL environment_end('CHI0')
  CALL stop_pp()
  !
CONTAINS
  !
  SUBROUTINE open_full()
    ahmin = HUGE(1.0_DP); ahmax = 0.0_DP
    IF (LEN_TRIM(filfull) > 0 .AND. ionode) THEN
      OPEN(91, FILE=TRIM(filfull), FORM='unformatted', ACCESS='stream', STATUS='replace')
      OPEN(92, FILE=TRIM(filfull) // '.G.txt', STATUS='replace')
      WRITE(92, '(a)') "# eps^-1(i, j, w) = eps^-1_{G_i G_j}(q, w), column-major, complex(8) per element, w outer;"
      WRITE(92, '(a)') "# chi_GG' = sum rho*(G) rho(G'), rho(G) = <c,k+q|e^{i(q+G)r}|v,k>"
      WRITE(92, '("# q crystal ",3F12.8,"  w = 0 .. ",F10.4," eV step ",F8.4,", ",I0," energies, ",I0," G")') &
           qcr, emax, de, nw, ng_max
      DO ig = 1, ng_max
        WRITE(92, '(I6,3I5,F16.8)') ig - 1, gm(:, ig), q2(ig)
      ENDDO
      CLOSE(92)
    ENDIF
  END SUBROUTINE open_full
  !
  SUBROUTINE finish_w(iw_, einvw)
    !! per energy: the head of eps^-1, the PSD check of the anti-Hermitian part of -eps^-1, the full matrix
    INTEGER, INTENT(IN) :: iw_
    COMPLEX(DP), INTENT(IN) :: einvw(ng_max, ng_max)
    INTEGER :: info
    einvqq(iw_) = einvw(iq_head, iq_head)
    IF (wgrid(iw_) > 0.0_DP) THEN
      ah = -(einvw - CONJG(TRANSPOSE(einvw))) / (0.0_DP, 2.0_DP)
      CALL ZHEEV('N', 'U', ng_max, ah, ng_max, evals, work, SIZE(work), rwork, info)
      ahmin = MIN(ahmin, evals(1)); ahmax = MAX(ahmax, evals(ng_max))
    ENDIF
    IF (LEN_TRIM(filfull) > 0 .AND. ionode) WRITE(91) einvw
  END SUBROUTINE finish_w
  !
#if defined(__CUDA)
  SUBROUTINE device_path()
    !! the same stages on the GPU (chi0_dev); per-energy results return to the host in chunks
    USE chi0_dev, ONLY : dev_init, dev_add_k, dev_kk, dev_inv_chunk, dev_free
    INTEGER :: iw0, nch, ikk, jkk, ivv, icc, ig2
    INTEGER, ALLOCATABLE :: idx_pg(:)
    REAL(DP), ALLOCATABLE :: xb(:,:)
    COMPLEX(DP), ALLOCATABLE :: ukk(:,:,:), uqq(:,:,:), epsc(:,:,:), einvc(:,:,:)
    ALLOCATE(idx_pg(ng_max), xb(nc, nv), ukk(dffts%nnr, npol, nv), uqq(dffts%nnr, npol, nc))
    DO ig2 = 1, ng_max
      idx_pg(ig2) = box_index(gm(:, ig2))           ! unnormalised inverse transform / N puts rho(G) at +G
    ENDDO
    CALL dev_init((/dffts%nr1, dffts%nr2, dffts%nr3/), (/dffts%nr1x, dffts%nr2x, dffts%nr3x/), dffts%nnr, npol, &
                  nv, nc, ng_max, idx_pg, nb, nw, cap)
    DO ikk = 1, nks
      jkk = kq(ikk)
      CALL start_clock('chi0_wfc')
      CALL wfc_r(ikk, 1, nv, (/0, 0, 0/), ukk)
      CALL wfc_r(jkk, nv + 1, nbnd, g0(:, ikk), uqq)
      CALL stop_clock('chi0_wfc')
      DO ivv = 1, nv
        DO icc = 1, nc
          xb(icc, ivv) = (et(nv + icc, jkk) - et(ivv, ikk)) * 0.5_DP * ha / de
        ENDDO
      ENDDO
      CALL start_clock('chi0_pairs')
      CALL dev_add_k(ukk, uqq, nv, xb)
      CALL stop_clock('chi0_pairs')
    ENDDO
    DEALLOCATE(ukk, uqq, xb, idx_pg)
    CALL start_clock('chi0_kk')
    CALL dev_kk(spin / (nks * omega), de / ha, eta / ha, nw)
    CALL stop_clock('chi0_kk')
    CALL open_full()
    nch = MAX(1, MIN(nw, INT(2.5d8 / (32.0_DP * ng_max * ng_max))))
    ALLOCATE(epsc(ng_max, ng_max, nch), einvc(ng_max, ng_max, nch))
    CALL start_clock('chi0_inv')
    DO iw0 = 1, nw, nch
      CALL dev_inv_chunk(sv, iw0, MIN(nch, nw - iw0 + 1), epsc, einvc)
      DO ikk = 1, MIN(nch, nw - iw0 + 1)
        eqq(iw0 + ikk - 1) = epsc(iq_head, iq_head, ikk)
        CALL finish_w(iw0 + ikk - 1, einvc(:, :, ikk))
      ENDDO
    ENDDO
    CALL stop_clock('chi0_inv')
    DEALLOCATE(epsc, einvc)
    CALL dev_free()
  END SUBROUTINE device_path
#endif
  !
  SUBROUTINE kk_fft(af, cf, ncol)
    !! cf(:, i) = spin/(Nk Omega) sum_b af(:, b) [K(i - b) - K(i + b)],  K(m) = 1/(m de + i eta) (Hartree), on the
    !! uniform grid w_i = (i-1) de: a zero-padded convolution and a correlation along b, by QE's 1D FFTs (cft_1z)
    USE fft_scalar, ONLY : cft_1z
    INTEGER, INTENT(IN) :: ncol
    COMPLEX(DP), INTENT(IN) :: af(ncol, nb)
    COMPLEX(DP), INTENT(OUT) :: cf(ncol, nw)
    INTEGER :: l, m, c0, nbl, c, jj
    COMPLEX(DP), ALLOCATABLE :: k1(:), k2(:), k1f(:), k2f(:), buf(:), bf(:), d(:), df(:)
    REAL(DP) :: sf, sb, scal
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
    ! the scalings of this FFT backend: forward and backward of a unit impulse
    d = (0.0_DP, 0.0_DP); d(1) = (1.0_DP, 0.0_DP)
    CALL cft_1z(d, 1, l, l, -1, df); sf = REAL(df(1))
    CALL cft_1z(d, 1, l, l, +1, df); sb = REAL(df(1))
    scal = spin / (nks * omega) / (sb * sf * sf * l)
    nbl = MAX(1, MIN(ncol, INT(2.0d6 / l)))          ! two block buffers of <= ~32 MB each (column c at (c-1) l)
    ALLOCATE(buf(l * nbl), bf(l * nbl))
    DO c0 = 1, ncol, nbl
      m = MIN(nbl, ncol - c0 + 1)
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
  INTEGER FUNCTION box_index(m)
    !! position of the Miller index m on the (serial) smooth FFT box, the layout of dffts%nl
    INTEGER, INTENT(IN) :: m(3)
    box_index = 1 + MODULO(m(1), dffts%nr1) + dffts%nr1x * (MODULO(m(2), dffts%nr2) + &
                dffts%nr2x * MODULO(m(3), dffts%nr3))
  END FUNCTION box_index
  !
  SUBROUTINE gbasis()
    INTEGER :: n, i1, i2, i3, cnt, k, l, ti(3)
    REAL(DP) :: v(3), t
    INTEGER, ALLOCATABLE :: tmp(:,:)
    REAL(DP), ALLOCATABLE :: tq(:)
    n = CEILING(SQRT(ecut_eps / tpiba2) / MINVAL(SQRT(SUM(bg**2, 1)))) + 2
    ALLOCATE(tmp(3, (2 * n + 1)**3), tq((2 * n + 1)**3))
    cnt = 0
    DO i1 = -n, n
      DO i2 = -n, n
        DO i3 = -n, n
          v = (qcr(1) + i1) * bg(:, 1) + (qcr(2) + i2) * bg(:, 2) + (qcr(3) + i3) * bg(:, 3)
          t = SUM(v**2) * tpiba2
          IF (t <= ecut_eps .OR. ALL((/i1, i2, i3/) == NINT(qc))) THEN
            cnt = cnt + 1; tmp(:, cnt) = (/i1, i2, i3/); tq(cnt) = t
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    ! stable insertion sort by |q+G|^2 rounded to 1e-10 bohr^-2, so that ties fall in enumeration order exactly as
    ! in chi0_ref.py (numpy's stable argsort of the same rounded values)
    DO k = 2, cnt
      t = tq(k); ti = tmp(:, k); l = k - 1
      DO WHILE (l >= 1)
        IF (ANINT(tq(l) * 1.0d10) <= ANINT(t * 1.0d10)) EXIT
        tq(l + 1) = tq(l); tmp(:, l + 1) = tmp(:, l); l = l - 1
      ENDDO
      tq(l + 1) = t; tmp(:, l + 1) = ti
    ENDDO
    ng_max = cnt
    ALLOCATE(gm(3, cnt), q2(cnt), idx_mg(cnt))
    gm = tmp(:, 1:cnt); q2 = tq(1:cnt)
    iq_head = 0
    DO k = 1, cnt
      ! where fwfft puts -G: through QE's own map (dffts%nl of the G-vector index), not a box formula, so that it
      ! holds for every FFT layout QE builds
      idx_mg(k) = 0
      DO l = 1, ngm
        IF (ALL(mill(:, l) == -gm(:, k))) THEN
          idx_mg(k) = dffts%nl(l); EXIT
        ENDIF
      ENDDO
      IF (idx_mg(k) == 0) CALL errore('chi0', 'local-field G outside the density G sphere', k)
      IF (ALL(gm(:, k) == NINT(qc))) iq_head = k
    ENDDO
    DEALLOCATE(tmp, tq)
  END SUBROUTINE gbasis
  !
  SUBROUTINE wfc_r(ikk, n1, n2, shift, u)
    !! u(:, :, n) = periodic part of bands n1..n2 at k-point ikk on the FFT box, coefficients moved by -shift
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
END PROGRAM chi0
