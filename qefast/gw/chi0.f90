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
! nbnd (0: all), fileout, filfull ('' : no full matrix), max_mem (GB, host path; 0: half the available), device (GPU path in CUDA builds), check_psd
! (number of energies, evenly spread, at which the anti-Hermitian part of -eps^-1 is checked; >= nw: all), dev_reserve
! (GB of device memory left free), dev_buffer (most transitions buffered on the device before a flush; 0: from memory),
! scissor (eV, shift of the empty bands).
!
PROGRAM chi0
  USE kinds,            ONLY : DP
  USE constants,        ONLY : rytoev, pi
  USE io_files,         ONLY : prefix, tmp_dir
  USE io_global,        ONLY : ionode, ionode_id, stdout
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : mp_startup
  USE mp_images,        ONLY : intra_image_comm
  USE mp_pools,         ONLY : npool, nproc_pool
  USE environment,      ONLY : environment_start, environment_end
  USE wvfct,            ONLY : nbnd_tot => nbnd, et
  USE klist,            ONLY : xk, nks, nelec, lgauss, ltetra
  USE cell_base,        ONLY : at, omega
  USE fft_base,         ONLY : dffts
  USE noncollin_module, ONLY : npol, noncolin, domag
  USE lsda_mod,         ONLY : lsda
  USE uspp,             ONLY : okvan
  USE control_flags,    ONLY : gamma_only
  USE qef_config,       ONLY : qef_note
  USE gw_common,        ONLY : ha, kq_grid, gbasis_q, box_index, wfc_r, pair_row, bin_add, kk_fft, eps_invert, &
                               mem_avail_gb
#if defined(__CUDA)
  USE chi0_dev,         ONLY : dev_cap
#endif
  IMPLICIT NONE
  CHARACTER(LEN=256) :: outdir, fileout, filfull
  REAL(DP) :: q(3), ecut_eps, emax, de, eta, max_mem
  INTEGER :: nbnd, ios
  LOGICAL :: device
  INTEGER :: check_psd
  REAL(DP) :: dev_reserve
  INTEGER :: dev_buffer
  REAL(DP) :: scissor
  NAMELIST / inputchi0 / prefix, outdir, q, ecut_eps, emax, de, eta, nbnd, fileout, filfull, max_mem, device, &
                         check_psd, dev_reserve, dev_buffer, scissor
  !
  INTEGER :: nv, nc, ik, jk, iv, ic, ig, nb, nw, iw, ng_max, iq_head, cap
  INTEGER, ALLOCATABLE :: kq(:), g0(:,:), gm(:,:), idx_mg(:)
  REAL(DP) :: qcr(3), qc(3), spin, need, ahmin, ahmax
  REAL(DP), ALLOCATABLE :: kcr(:,:), q2(:), sv(:), wgrid(:), tx(:)
  COMPLEX(DP), ALLOCATABLE :: a(:,:,:), uk(:,:,:), uq(:,:,:), f(:), rho(:,:), &
                              chi(:,:,:), epsm(:,:), work(:), eqq(:), einvqq(:), ah(:,:)
  REAL(DP), ALLOCATABLE :: rwork(:), evals(:)
  LOGICAL :: needwf = .TRUE.
  !
  CALL mp_startup()
  CALL environment_start('CHI0')
  prefix = 'pwscf'; outdir = './'; q = 0.0_DP; ecut_eps = 4.0_DP; emax = 60.0_DP; de = 0.02_DP; eta = 0.1_DP
  nbnd = 0; fileout = 'chi0.dat'; filfull = ' '; max_mem = 0.0_DP; device = .TRUE.; check_psd = 100
  dev_reserve = 2.0_DP; dev_buffer = 0; scissor = 0.0_DP
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
  CALL mp_bcast(check_psd, ionode_id, intra_image_comm); CALL mp_bcast(dev_reserve, ionode_id, intra_image_comm)
  CALL mp_bcast(dev_buffer, ionode_id, intra_image_comm); CALL mp_bcast(scissor, ionode_id, intra_image_comm)
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
  ! scissor: every empty band moved up by scissor (eV); wavefunctions and pair densities are unchanged
  IF (scissor /= 0.0_DP) et(nv + 1:, :) = et(nv + 1:, :) + scissor / rytoev
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
  CALL kq_grid(kcr, qcr, kq, g0)
  !
  ! ---- local-field basis |q+G|^2 <= ecut_eps (Ry = bohr^-2), sorted by |q+G| (stable), G_Q always in
  CALL gbasis_q(qcr, NINT(qc), ecut_eps, ng_max, gm, q2, idx_mg, iq_head)
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
  need = (16.0_DP * (nb + nw) * ng_max * ng_max + 6.4d7) / 1.0d9   ! A, chi on the grid, FFT block buffers
  !
  ALLOCATE(wgrid(nw), eqq(nw), einvqq(nw), ah(ng_max, ng_max), rwork(3 * ng_max), evals(ng_max), &
           work(64 * ng_max), epsm(ng_max, ng_max))
  DO iw = 1, nw
    wgrid(iw) = (iw - 1) * de
  ENDDO
#if defined(__CUDA)
  IF (device) THEN
    cap = dev_cap(dffts%nnr, npol, nv, nc, ng_max, nb, nw, dev_reserve)
    IF (dev_buffer > 0) cap = MIN(cap, MAX(dev_buffer, nc))
    WRITE(stdout, '(5x,"chi0: device transition buffer ",I0," (",F6.2," GB left free)")') cap, dev_reserve
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
  ! host memory: binned matrices and chi on the grid; max_mem defaults to half of the memory available now
  IF (max_mem <= 0.0_DP) max_mem = 0.5_DP * mem_avail_gb()
  WRITE(stdout, '(5x,"chi0: host memory estimate ",F6.2," GB (",I0," bins and ",I0," energies x ",I0,"^2), limit ", &
       & F6.2," GB")') need, nb, nw, ng_max, max_mem
  IF (need > max_mem) CALL errore('chi0', 'memory estimate exceeds max_mem; raise de, lower emax or ecut_eps', 1)
  ! ---- A_b: pair densities and binning
  ALLOCATE(a(ng_max, ng_max, nb), uk(dffts%nnr, npol, nv), uq(dffts%nnr, npol, nc), f(dffts%nnr), &
           rho(ng_max, nc), tx(nc))
  a = (0.0_DP, 0.0_DP)
  DO ik = 1, nks
    jk = kq(ik)
    CALL start_clock('chi0_wfc')
    CALL wfc_r(ik, 1, nv, (/0, 0, 0/), uk)
    CALL wfc_r(jk, nv + 1, nbnd, g0(:, ik), uq)
    CALL stop_clock('chi0_wfc')
    DO iv = 1, nv
      CALL start_clock('chi0_pairs')
      CALL pair_row(uk(:, :, iv), uq, nc, idx_mg, rho, f)
      CALL stop_clock('chi0_pairs')
      CALL start_clock('chi0_bin')
      ! linear binning of D = e_c(k+q) - e_v(k) into bins b and b+1
      DO ic = 1, nc
        tx(ic) = (et(nv + ic, jk) - et(iv, ik)) * 0.5_DP * ha / de
      ENDDO
      CALL bin_add(a, rho, nc, tx)
      CALL stop_clock('chi0_bin')
    ENDDO
  ENDDO
  DEALLOCATE(uk, uq, f, rho, tx)
  !
  ! ---- chi0 on the energy grid by FFT (Kramers-Kronig as a convolution and a correlation), then eps^-1 per energy
  ALLOCATE(chi(ng_max, ng_max, nw))
  CALL start_clock('chi0_kk')
  CALL kk_fft(a, chi, ng_max * ng_max, nb, nw, de, eta, spin / (nks * omega))
  CALL stop_clock('chi0_kk')
  DEALLOCATE(a)
  CALL open_full()
  CALL start_clock('chi0_inv')
  DO iw = 1, nw
    CALL eps_invert(sv, chi(:, :, iw), epsm, eqq(iw), iq_head)
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
    IF (wgrid(iw_) > 0.0_DP .AND. MOD(iw_ - 1, MAX(1, nw / MAX(1, check_psd))) == 0) THEN
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
END PROGRAM chi0
