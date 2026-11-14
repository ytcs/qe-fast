! Copyright (C) 2026 The qe-fast authors. Portions derive from Quantum ESPRESSO; see NOTICE.
! SPDX-License-Identifier: GPL-2.0-or-later
!
! gamma_ii.x -- impact-ionisation rate of a hot electron from the RPA screened interaction (book E-13, issue #1).
! Host path; the stages are those of qefast/gw/gamma_ref.py, the reference this program is checked against:
!
!   Im Sigma_nk = 1/(Nk Omega) sum_{q /= 0} sum_{n' empty} phi^H Im W(q, w) phi,   w = e_nk - e_n'k-q > 0
!   phi(G)      = <n, k| e^{i(q+G)r} |n', k-q>
!   Im W        = v^1/2 [(eps^-1 - eps^-1^H) / 2i] v^1/2,  eps^-1 as chi0.x (RPA, local fields) on a uniform grid
!   Gamma_ii    = 2 |Im Sigma| / hbar   (z = 1)
!
! Initial states: every conduction state up to emax_carrier above the conduction minimum; final states: the empty
! states between the conduction minimum and the initial energy. Im W is interpolated linearly between grid
! energies. The q = 0 term (head and wings of W) is left out; its weight is 1/Nk.
!
! Supported: as chi0.x (norm-conserving, fixed occupations, no spin polarisation or magnetisation, a pw.x nscf on
! a full k grid without symmetry, one MPI process).
!
! Input namelist &inputgii: prefix, outdir, emax_carrier (eV), ecut_eps (Ry), de, eta (eV), scissor (eV),
! fileout (n, k, E - E_cbm, Im Sigma, Gamma per initial state), filpart ('': none; the rate into every final state
! (n', k'), whose sum over final states is Gamma: the cascade draws the scattered carrier from it), max_mem (GB;
! 0: half the available).
!
PROGRAM gamma_ii
  USE kinds,            ONLY : DP
  USE constants,        ONLY : rytoev, pi
  USE io_files,         ONLY : prefix, tmp_dir
  USE io_global,        ONLY : ionode, ionode_id, stdout
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : mp_startup
  USE mp_images,        ONLY : intra_image_comm
  USE mp_pools,         ONLY : npool, nproc_pool
  USE environment,      ONLY : environment_start, environment_end
  USE wvfct,            ONLY : nbnd, et
  USE klist,            ONLY : xk, nks, nelec, lgauss, ltetra
  USE cell_base,        ONLY : at, omega
  USE fft_base,         ONLY : dffts
  USE noncollin_module, ONLY : npol, noncolin, domag
  USE lsda_mod,         ONLY : lsda
  USE uspp,             ONLY : okvan
  USE control_flags,    ONLY : gamma_only
  USE gw_common,        ONLY : ha, kq_grid, gbasis_q, wfc_r, pair_row, bin_add, kk_fft, eps_invert, mem_avail_gb
  IMPLICIT NONE
  REAL(DP), PARAMETER :: hbar_ha = 2.4188843265857d-17   ! s, hbar / E_h
  CHARACTER(LEN=256) :: outdir, fileout, filpart
  REAL(DP) :: emax_carrier, ecut_eps, de, eta, scissor, max_mem
  INTEGER :: ios
  NAMELIST / inputgii / prefix, outdir, emax_carrier, ecut_eps, de, eta, scissor, fileout, filpart, max_mem
  !
  INTEGER :: nv, nc, ni, nmax, ik, jk, iq, nq, iv, ic, n, jn, iw, nb, nw, ng, ihead, np, p, j0, cnt, nrow
  INTEGER, ALLOCATABLE :: kq(:), g0(:,:), gm(:,:), idx_mg(:), pn(:), pj(:), pg(:), sel(:)
  REAL(DP) :: qcr(3), spin, wmax, ecbm, evbm, need, w, x, fr, lo, smax
  REAL(DP), ALLOCATABLE :: kcr(:,:), q2(:), sv(:), e(:,:), tx(:), imsig(:,:), pf(:), pc(:), sw(:)
  COMPLEX(DP), ALLOCATABLE :: a(:,:,:), chi(:,:,:), epsm(:,:), uk(:,:,:), uq(:,:,:), f(:), rho(:,:), &
                              phi(:,:,:), xm(:,:), ym(:,:)
  COMPLEX(DP) :: eh
  LOGICAL :: needwf = .TRUE.
  !
  CALL mp_startup()
  CALL environment_start('GAMMA_II')
  prefix = 'pwscf'; outdir = './'; emax_carrier = 6.0_DP; ecut_eps = 4.0_DP; de = 0.05_DP; eta = 0.05_DP
  scissor = 0.0_DP; fileout = 'gamma_ii.dat'; filpart = ' '; max_mem = 0.0_DP
  IF (ionode) THEN
    CALL input_from_file()
    READ(5, inputgii, IOSTAT=ios)
  ENDIF
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  IF (ios /= 0) CALL errore('gamma_ii', 'reading namelist inputgii', ABS(ios))
  CALL mp_bcast(prefix, ionode_id, intra_image_comm); CALL mp_bcast(outdir, ionode_id, intra_image_comm)
  CALL mp_bcast(emax_carrier, ionode_id, intra_image_comm); CALL mp_bcast(ecut_eps, ionode_id, intra_image_comm)
  CALL mp_bcast(de, ionode_id, intra_image_comm); CALL mp_bcast(eta, ionode_id, intra_image_comm)
  CALL mp_bcast(scissor, ionode_id, intra_image_comm); CALL mp_bcast(fileout, ionode_id, intra_image_comm)
  CALL mp_bcast(filpart, ionode_id, intra_image_comm); CALL mp_bcast(max_mem, ionode_id, intra_image_comm)
  tmp_dir = TRIM(outdir) // '/'
  CALL read_file_new(needwf)
  !
  IF (npool * nproc_pool > 1) CALL errore('gamma_ii', 'run with one MPI process for now', 1)
  IF (okvan) CALL errore('gamma_ii', 'ultrasoft or PAW pseudopotentials are not supported (norm-conserving only)', 1)
  IF (lgauss .OR. ltetra) CALL errore('gamma_ii', 'only fixed occupations (insulators) are supported', 1)
  IF (lsda) CALL errore('gamma_ii', 'spin-polarised (lsda) runs are not supported', 1)
  IF (noncolin .AND. domag) CALL errore('gamma_ii', 'magnetic runs are not supported (time reversal is assumed)', 1)
  IF (gamma_only) CALL errore('gamma_ii', 'gamma_only runs are not supported', 1)
  IF (dffts%nnr /= dffts%nr1x * dffts%nr2x * dffts%nr3x) CALL errore('gamma_ii', 'distributed FFT', 1)
  IF (noncolin) THEN
    nv = NINT(nelec); spin = 1.0_DP
  ELSE
    nv = NINT(nelec / 2.0_DP); spin = 2.0_DP
  ENDIF
  nc = nbnd - nv
  IF (nc < 1) CALL errore('gamma_ii', 'no empty bands', 1)
  IF (scissor /= 0.0_DP) et(nv + 1:, :) = et(nv + 1:, :) + scissor / rytoev
  ALLOCATE(e(nbnd, nks))
  e = et(1:nbnd, 1:nks) * rytoev
  ecbm = MINVAL(e(nv + 1, :)); evbm = MAXVAL(e(nv, :))
  ! initial states: bands nv+1..nmax; the highest band must lie above the carrier cutoff everywhere
  nmax = nv
  DO WHILE (nmax < nbnd)
    IF (MINVAL(e(nmax + 1, :)) > ecbm + emax_carrier) EXIT
    nmax = nmax + 1
  ENDDO
  IF (nmax >= nbnd - 1 .OR. MINVAL(e(nbnd, :)) < ecbm + emax_carrier) &
    CALL errore('gamma_ii', 'not enough bands above the carrier cutoff; raise nbnd in the nscf', 1)
  ni = nmax - nv
  wmax = emax_carrier + 0.5_DP
  nb = CEILING(wmax / de) + 2
  nw = FLOOR(wmax / de + 1.0d-9) + 1                ! grid 0, de, .., wmax
  ALLOCATE(kcr(3, nks), kq(nks), g0(3, nks), imsig(nbnd, nks))
  DO ik = 1, nks
    kcr(:, ik) = MATMUL(xk(:, ik), at)
  ENDDO
  nq = 0
  DO iq = 1, nks
    IF (MAXVAL(ABS(kcr(:, iq) - NINT(kcr(:, iq)))) > 1.0d-8) nq = nq + 1
  ENDDO
  WRITE(stdout, '(/5x,"gamma_ii: ",I0," k, bands ",I0,"-",I0," as initial states, E_cbm ",F8.3," eV, gap ",F7.3, &
       & " eV, W on ",I0," energies to ",F6.2," eV; ",I0," q (q = 0 left out, weight ",F7.4,")")') &
       nks, nv + 1, nmax, ecbm, ecbm - evbm, nw, wmax, nq, 1.0_DP / nks
  IF (max_mem <= 0.0_DP) max_mem = 0.5_DP * mem_avail_gb()
  IF (LEN_TRIM(filpart) > 0 .AND. ionode) THEN
    OPEN(93, FILE=TRIM(filpart), STATUS='replace')
    WRITE(93, '(a)') '# n k nf kf Gamma(1/s) E_init-E_cbm(eV) E_final-E_cbm(eV): rate of (n,k) -> (nf,kf), q = k - kf'
  ENDIF
  imsig = 0.0_DP
  ALLOCATE(uk(dffts%nnr, npol, nbnd), uq(dffts%nnr, npol, nc), f(dffts%nnr), tx(nc))
  !
  cnt = 0
  DO iq = 1, nks
    IF (MAXVAL(ABS(kcr(:, iq) - NINT(kcr(:, iq)))) <= 1.0d-8) CYCLE
    cnt = cnt + 1
    qcr = kcr(:, iq) - CEILING(kcr(:, iq) - 0.5_DP - 1.0d-9)
    CALL gbasis_q(qcr, (/0, 0, 0/), ecut_eps, ng, gm, q2, idx_mg, ihead)
    need = (16.0_DP * (nb + nw) * ng * ng + 16.0_DP * ng * nc * (ni + 1) + 6.4d7) / 1.0d9
    IF (need > max_mem) CALL errore('gamma_ii', 'memory estimate exceeds max_mem; raise de or lower ecut_eps', 1)
    ALLOCATE(sv(ng), a(ng, ng, nb), rho(ng, nc))
    sv = SQRT(4.0_DP * pi / q2)
    CALL kq_grid(kcr, qcr, kq, g0)
    ! ---- A_b (valence -> conduction transitions at transfer q), then chi0 by FFT and Im W per energy
    a = (0.0_DP, 0.0_DP)
    DO ik = 1, nks
      jk = kq(ik)
      CALL start_clock('gii_wfc')
      CALL wfc_r(ik, 1, nv, (/0, 0, 0/), uk)
      CALL wfc_r(jk, nv + 1, nbnd, g0(:, ik), uq)
      CALL stop_clock('gii_wfc')
      DO iv = 1, nv
        CALL start_clock('gii_pairs')
        CALL pair_row(uk(:, :, iv), uq, nc, idx_mg, rho, f)
        CALL stop_clock('gii_pairs')
        CALL start_clock('gii_bin')
        DO ic = 1, nc
          tx(ic) = (et(nv + ic, jk) - et(iv, ik)) * 0.5_DP * ha / de
        ENDDO
        CALL bin_add(a, rho, nc, tx)
        CALL stop_clock('gii_bin')
      ENDDO
    ENDDO
    DEALLOCATE(rho)
    ALLOCATE(chi(ng, ng, nw), epsm(ng, ng))
    CALL start_clock('gii_kk')
    CALL kk_fft(a, chi, ng * ng, nb, nw, de, eta, spin / (nks * omega))
    CALL stop_clock('gii_kk')
    DEALLOCATE(a)
    CALL start_clock('gii_w')
    DO iw = 1, nw                                    ! chi(:, :, iw) <- Im W(w_iw), Hartree bohr^3
      CALL eps_invert(sv, chi(:, :, iw), epsm, eh, 0)
      chi(:, :, iw) = (epsm - CONJG(TRANSPOSE(epsm))) / (0.0_DP, 2.0_DP)
      DO jn = 1, ng
        chi(:, jn, iw) = sv(:) * chi(:, jn, iw) * sv(jn)
      ENDDO
    ENDDO
    CALL stop_clock('gii_w')
    DEALLOCATE(epsm)
    ! ---- contraction: final state (n', k' = ik), initial state (n, k = k' + q)
    ALLOCATE(phi(ng, ni, nc), pn(ni * nc), pj(ni * nc), pg(ni * nc), pf(ni * nc), pc(ni * nc), sel(ni * nc), &
             sw(ni * nc), xm(ng, ni * nc), ym(ng, ni * nc))
    DO ik = 1, nks
      jk = kq(ik)
      CALL start_clock('gii_wfc')
      CALL wfc_r(ik, nv + 1, nbnd, (/0, 0, 0/), uk)
      CALL wfc_r(jk, nv + 1, nmax, g0(:, ik), uq)
      CALL stop_clock('gii_wfc')
      CALL start_clock('gii_pairs')
      DO jn = 1, nc
        CALL pair_row(uk(:, :, jn), uq, ni, idx_mg, phi(:, :, jn), f)
      ENDDO
      CALL stop_clock('gii_pairs')
      CALL start_clock('gii_contract')
      np = 0
      DO n = nv + 1, nmax
        IF (e(n, jk) > ecbm + emax_carrier) CYCLE
        DO jn = 1, nc
          w = e(n, jk) - e(nv + jn, ik)
          IF (w <= 0.0_DP) CYCLE
          x = w / de; j0 = INT(x); fr = x - j0
          IF (j0 + 2 > nw) CYCLE
          np = np + 1; pn(np) = n; pj(np) = jn; pg(np) = j0 + 1; pf(np) = fr
        ENDDO
      ENDDO
      pc(1:np) = 0.0_DP
      DO iw = 1, nw
        nrow = 0
        DO p = 1, np
          IF (pg(p) == iw) THEN
            nrow = nrow + 1; sel(nrow) = p; sw(nrow) = 1.0_DP - pf(p)
          ELSEIF (pg(p) + 1 == iw) THEN
            nrow = nrow + 1; sel(nrow) = p; sw(nrow) = pf(p)
          ENDIF
        ENDDO
        IF (nrow == 0) CYCLE
        DO p = 1, nrow
          xm(:, p) = phi(:, pn(sel(p)) - nv, pj(sel(p)))
        ENDDO
        CALL ZGEMM('N', 'N', ng, nrow, ng, (1.0_DP, 0.0_DP), chi(1, 1, iw), ng, xm, ng, (0.0_DP, 0.0_DP), ym, ng)
        DO p = 1, nrow
          pc(sel(p)) = pc(sel(p)) + sw(p) * REAL(DOT_PRODUCT(xm(:, p), ym(:, p)), DP)
        ENDDO
      ENDDO
      DO p = 1, np
        imsig(pn(p), jk) = imsig(pn(p), jk) + pc(p)
        IF (LEN_TRIM(filpart) > 0 .AND. ionode) WRITE(93, '(4I6,ES16.8,2F10.4)') pn(p), jk, nv + pj(p), ik, &
             2.0_DP * ABS(pc(p) / (nks * omega)) / hbar_ha, e(pn(p), jk) - ecbm, e(nv + pj(p), ik) - ecbm
      ENDDO
      CALL stop_clock('gii_contract')
    ENDDO
    DEALLOCATE(phi, pn, pj, pg, pf, pc, sel, sw, xm, ym, chi, sv, gm, q2, idx_mg)
    WRITE(stdout, '(5x,"  q ",I0,"/",I0,"  NG ",I0)') cnt, nq, ng
    FLUSH(stdout)
  ENDDO
  IF (LEN_TRIM(filpart) > 0 .AND. ionode) CLOSE(93)
  imsig = imsig / (nks * omega)
  !
  ! ---- output (columns as gamma_ref.py)
  smax = -HUGE(1.0_DP)
  IF (ionode) THEN
    OPEN(90, FILE=TRIM(fileout), STATUS='replace')
    WRITE(90, '(a)') '# band k E-E_cbm(eV) ImSigma(eV) Gamma_ii(1/s), z = 1, q = 0 term left out'
    DO ik = 1, nks
      DO n = nv + 1, nmax
        IF (e(n, ik) > ecbm + emax_carrier) CYCLE
        WRITE(90, '(2I6,F10.4,2ES16.6)') n, ik, e(n, ik) - ecbm, imsig(n, ik) * ha, &
             2.0_DP * ABS(imsig(n, ik)) / hbar_ha
        smax = MAX(smax, imsig(n, ik) * ha)
      ENDDO
    ENDDO
    CLOSE(90)
  ENDIF
  WRITE(stdout, '(5x,"gamma_ii: largest Im Sigma ",ES10.2," eV (must be <= 0)")') smax
  lo = 0.0_DP
  DO WHILE (lo < emax_carrier - 1.0d-9)
    CALL bin_report(lo, lo + 0.5_DP)
    lo = lo + 0.5_DP
  ENDDO
  CALL print_clock('gii_wfc'); CALL print_clock('gii_pairs'); CALL print_clock('gii_bin')
  CALL print_clock('gii_kk'); CALL print_clock('gii_w'); CALL print_clock('gii_contract')
  CALL environment_end('GAMMA_II')
  CALL stop_pp()
  !
CONTAINS
  !
  SUBROUTINE bin_report(e0, e1)
    REAL(DP), INTENT(IN) :: e0, e1
    INTEGER :: m, kk, nn
    REAL(DP) :: s, g, gmx
    m = 0; s = 0.0_DP; gmx = 0.0_DP
    DO kk = 1, nks
      DO nn = nv + 1, nmax
        IF (e(nn, kk) > ecbm + emax_carrier) CYCLE
        IF (e(nn, kk) - ecbm < e0 .OR. e(nn, kk) - ecbm >= e1) CYCLE
        g = 2.0_DP * ABS(imsig(nn, kk)) / hbar_ha
        m = m + 1; s = s + g; gmx = MAX(gmx, g)
      ENDDO
    ENDDO
    IF (m > 0) WRITE(stdout, '(5x,"gamma_ii: E-E_cbm ",F4.1,"-",F4.1," eV: ",I4," states, Gamma mean ",ES9.2, &
         & " 1/s, max ",ES9.2)') e0, e1, m, s / m, gmx
  END SUBROUTINE bin_report
  !
END PROGRAM gamma_ii
