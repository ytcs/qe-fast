!
! Copyright (C) 2026 The qe-fast authors
! SPDX-License-Identifier: GPL-2.0-or-later
!
!----------------------------------------------------------------------------
MODULE ph_imom
  !----------------------------------------------------------------------------
  !! Inverse spectral moments of the phonon perturbation (input flag inverse_moments).
  !!
  !! For the perturbation dV_mu of pattern mu at the wavevector q, with |c> the empty and |v> the occupied
  !! Kohn-Sham states and Delta = e_c(k+q) - e_v(k),
  !!
  !!    I_n(mu,nu) = sum_k w_k sum_v sum_c <v|dV_mu^+|c><c|dV_nu|v> / Delta^n ,   n = 2, 3, 4,
  !!
  !! with the spin factor in w_k (sum_k w_k = 2). They come from two Sternheimer solves per occupied state,
  !!    psi1 = -(H - e_v)^-1 P_c dV psi_v ,   psi2 = (H - e_v)^-1 P_c psi1 ,
  !!    I_2 = <psi1|psi1>,  I_3 = <psi1|psi2>,  I_4 = <psi2|psi2>  (summed over k and v),
  !! so no empty states are needed. Two vertices are used: the self-consistent dV = dV_bare + dV_Hxc
  !! (the converged response of each irrep; at q = 0 without the macroscopic field) and the bare dV_bare
  !! (local and nonlocal pseudopotential). With the self-consistent vertex, 2 I_3 is the coefficient of omega^2 in
  !! the electronic part of the force constants on the imaginary frequency axis, for the adiabatic kernel of the
  !! run; with the bare vertex it is the plain Kohn-Sham pair sum.
  !!
  !! Each irrep stores only its converged dV_Hxc; when every irrep of q is done, the solves run per k for all
  !! modes at once. The matrices are symmetrized like the dynamical matrix and written for the star of q in the
  !! layout of the dynamical-matrix file (Cartesian components, Rydberg atomic units, Ry^(2-n) bohr^-2, numbers
  !! in exponent format) to
  !!    <fildyn>.i2s, .i3s, .i4s (self-consistent vertex) and <fildyn>.i2b, .i3b, .i4b (bare vertex).
  !! At q = 0 the sum over all atom pairs is the response to a rigid translation; there the self-consistent
  !! psi1 of the translation equals -P_c d/dr psi_v exactly (also with nonlocal pseudopotentials), which is
  !! checked and printed. Its 4 I_3 = 2 sum |<c|p|v>|^2 / Delta (Hartree) equals the number of valence electrons
  !! only for a local potential (f-sum rule; H2 with a local pseudopotential: 2.01 for 2). Nonlocal projectors
  !! make p differ from the velocity and break that identity (Si, SG15 ONCV, 16^3 k: 8.80 for 8).
  !!
  !! Scope: norm-conserving pseudopotentials, insulators, collinear unpolarized spin, no Hubbard U,
  !! ecutrho = 4 ecutwfc; every irrep of a q point must be computed in the same run.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  LOGICAL, PUBLIC :: inverse_moments = .FALSE.
  !
  COMPLEX(DP), ALLOCATABLE :: dvall(:,:)     ! dV_Hxc of every mode of the current q (dense grid)
  LOGICAL, ALLOCATABLE :: mode_done(:)
  !
  INTEGER, PARAMETER :: nkind = 4             ! psi1 scf, psi2 scf, psi1 bare, psi2 bare
  CHARACTER(LEN=4), PARAMETER :: ext(6) = (/ '.i2s', '.i3s', '.i4s', '.i2b', '.i3b', '.i4b' /)
  !
  PUBLIC :: imom_check, imom_irrep, imom_finish
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE imom_check()
    !------------------------------------------------------------------------
    !! Stops for the cases the moments do not handle.
    USE uspp,             ONLY : okvan
    USE noncollin_module, ONLY : noncolin
    USE lsda_mod,         ONLY : lsda
    USE klist,            ONLY : lgauss, ltetra
    USE ldaU,             ONLY : lda_plus_u
    USE control_lr,       ONLY : lgamma_gamma
    USE gvecs,            ONLY : doublegrid
    USE control_ph,       ONLY : lmultipole
    !
    IF (.NOT. inverse_moments) RETURN
    IF (okvan) CALL errore('imom_check', 'inverse_moments: norm-conserving pseudopotentials only', 1)
    IF (noncolin .OR. lsda) CALL errore('imom_check', 'inverse_moments: collinear, unpolarized only', 1)
    IF (lgauss .OR. ltetra) CALL errore('imom_check', 'inverse_moments: insulators only', 1)
    IF (lda_plus_u) CALL errore('imom_check', 'inverse_moments: no Hubbard U', 1)
    IF (lgamma_gamma) CALL errore('imom_check', 'inverse_moments: not with the Gamma-only trick', 1)
    IF (doublegrid) CALL errore('imom_check', 'inverse_moments: needs ecutrho = 4 ecutwfc', 1)
    IF (lmultipole) CALL errore('imom_check', 'inverse_moments: not with lmultipole', 1)
  END SUBROUTINE imom_check
  !
  !--------------------------------------------------------------------------
  SUBROUTINE imom_irrep(imode0, npe, dvscf)
    !------------------------------------------------------------------------
    !! After the irrep (modes imode0+1 .. imode0+npe) has converged: keep its dV_Hxc.
    USE ions_base,        ONLY : nat
    USE fft_base,         ONLY : dfftp
    USE noncollin_module, ONLY : nspin_mag
    !
    INTEGER, INTENT(IN) :: imode0, npe
    COMPLEX(DP), INTENT(IN) :: dvscf(dfftp%nnr, nspin_mag, npe)
    INTEGER :: ipert
    !
    IF (.NOT. inverse_moments) RETURN
    IF (.NOT. ALLOCATED(dvall)) THEN
       ALLOCATE (dvall(dfftp%nnr, 3 * nat), mode_done(3 * nat))
       dvall = (0.d0, 0.d0)
       mode_done = .FALSE.
    ENDIF
    DO ipert = 1, npe
       dvall(:, imode0 + ipert) = dvscf(:, 1, ipert)
       mode_done(imode0 + ipert) = .TRUE.
    ENDDO
  END SUBROUTINE imom_irrep
  !
  !--------------------------------------------------------------------------
  SUBROUTINE imom_finish()
    !------------------------------------------------------------------------
    !! Every mode of q done: psi1 and psi2 for both vertices, every k and every mode; the six matrices,
    !! symmetrized, written for the star of q.
    USE io_global,        ONLY : stdout, ionode
    USE ions_base,        ONLY : nat
    USE cell_base,        ONLY : at, bg, tpiba
    USE symm_base,        ONLY : s, invs, irt, nsym, t_rev
    USE klist,            ONLY : xk, wk, ngk, igk_k, nelec
    USE gvect,            ONLY : g
    USE wvfct,            ONLY : nbnd, npwx, et
    USE wavefunctions,    ONLY : evc
    USE noncollin_module, ONLY : npol
    USE lsda_mod,         ONLY : current_spin
    USE uspp,             ONLY : vkb
    USE uspp_init,        ONLY : init_us_2
    USE buffers,          ONLY : get_buffer
    USE mp,               ONLY : mp_sum
    USE mp_pools,         ONLY : inter_pool_comm, intra_pool_comm
    USE units_lr,         ONLY : iuwfc, lrwfc
    USE control_lr,       ONLY : lgamma, nbnd_occ, tr2_ph
    USE qpoint,           ONLY : xq, nksq, ikks, ikqs
    USE eqv,              ONLY : dpsi, dvpsi, evq
    USE modes,            ONLY : u
    USE phus,             ONLY : alphap
    USE lrus,             ONLY : becp1
    USE lr_symm_base,     ONLY : minus_q, irotmq, nsymq, rtau
    USE output,           ONLY : fildyn
    USE apply_dpot_mod,   ONLY : apply_dpot_bands
    !
    INTEGER :: ik, ikk, ikq, npw, npwq, mu, nu, nocc, i, iu, nq, isq(48), imq, nmodes, na, nb, ipol, &
               ivert, lter, nfail, ig, ib
    REAL(DP) :: sxq(3, 48), tr(6), thresh, anorm, tres(2, 3)
    LOGICAL :: conv_root
    REAL(DP), ALLOCATABLE :: h_diag(:,:)
    COMPLEX(DP), ALLOCATABLE :: psi(:,:,:,:), mat(:,:,:), dvb(:,:), aux(:,:), dvloc(:,:)
    COMPLEX(DP) :: cmu
    INTEGER, EXTERNAL :: find_free_unit
    EXTERNAL ch_psi_all, cg_psi
    !
    IF (.NOT. inverse_moments .OR. .NOT. ALLOCATED(dvall)) RETURN
    nmodes = 3 * nat
    IF (.NOT. ALL(mode_done)) THEN
       WRITE(stdout, '(/,5x,"inverse moments: not every irrep of this q was computed in this run; none written")')
       CALL cleanup()
       RETURN
    ENDIF
    CALL start_clock('imom')
    thresh = MIN(1.d-1 * SQRT(tr2_ph / npol), 1.d-2)
    ALLOCATE (psi(npwx*npol, nbnd, nmodes, nkind), mat(nmodes, nmodes, 6), h_diag(npwx*npol, nbnd), &
              dvb(npwx*npol, nbnd), aux(npwx*npol, nbnd), dvloc(SIZE(dvall, 1), 1))
    mat = (0.d0, 0.d0)
    tres = 0.d0
    nfail = 0
    current_spin = 1
    !
    DO ik = 1, nksq
       ikk = ikks(ik)
       ikq = ikqs(ik)
       npw = ngk(ikk)
       npwq = ngk(ikq)
       nocc = nbnd_occ(ikk)
       CALL get_buffer(evc, lrwfc, iuwfc, ikk)
       !$acc update device(evc)
       IF (.NOT. lgamma) THEN
          CALL get_buffer(evq, lrwfc, iuwfc, ikq)
          !$acc update device(evq)
       ENDIF
       CALL init_us_2(npwq, igk_k(1, ikq), xk(1, ikq), vkb, .true.)
       !$acc update host(vkb)
       CALL g2_kin(ikq)
       CALL h_prec(ik, evq, h_diag)
       !
       psi = (0.d0, 0.d0)
       DO mu = 1, nmodes
          CALL dvqpsi_us(ik, u(1, mu), .FALSE., becp1, alphap)     ! dV_bare psi_v -> dvpsi
          dvb = dvpsi
          DO ivert = 1, 2
             dvpsi = dvb
             IF (ivert == 1) THEN                                  ! + dV_Hxc psi_v
                aux = (0.d0, 0.d0)
                dvloc(:, 1) = dvall(:, mu)
                CALL apply_dpot_bands(ik, nocc, dvloc, evc, aux)
                dvpsi = dvpsi + aux
             ENDIF
             ! psi1 = -(H - e_v)^-1 P_c dV psi_v  (orthogonalize returns -P_c dvpsi)
             CALL solve(psi(:, :, mu, 2 * ivert - 1))
             ! psi2 = (H - e_v)^-1 P_c psi1: the solver returns -psi2 for the source -P_c psi1
             dvpsi = psi(:, :, mu, 2 * ivert - 1)
             CALL solve(psi(:, :, mu, 2 * ivert))
             psi(:, :, mu, 2 * ivert) = -psi(:, :, mu, 2 * ivert)
          ENDDO
       ENDDO
       !
       DO nu = 1, nmodes
          DO mu = 1, nmodes
             mat(mu, nu, 1) = mat(mu, nu, 1) + wk(ikk) * dot(psi(:, :, mu, 1), psi(:, :, nu, 1), nocc)
             mat(mu, nu, 2) = mat(mu, nu, 2) + wk(ikk) * dot(psi(:, :, mu, 1), psi(:, :, nu, 2), nocc)
             mat(mu, nu, 3) = mat(mu, nu, 3) + wk(ikk) * dot(psi(:, :, mu, 2), psi(:, :, nu, 2), nocc)
             mat(mu, nu, 4) = mat(mu, nu, 4) + wk(ikk) * dot(psi(:, :, mu, 3), psi(:, :, nu, 3), nocc)
             mat(mu, nu, 5) = mat(mu, nu, 5) + wk(ikk) * dot(psi(:, :, mu, 3), psi(:, :, nu, 4), nocc)
             mat(mu, nu, 6) = mat(mu, nu, 6) + wk(ikk) * dot(psi(:, :, mu, 4), psi(:, :, nu, 4), nocc)
          ENDDO
       ENDDO
       !
       IF (lgamma) THEN
          ! rigid translation along ipol: pattern coefficients c_mu = sum_i conj(u_i,mu) d_i with d_i = 1 on the
          ! ipol component of every atom; the self-consistent psi1 must equal -P_c d/dr_ipol psi_v
          DO ipol = 1, 3
             aux = (0.d0, 0.d0)
             DO mu = 1, nmodes
                cmu = (0.d0, 0.d0)
                DO na = 1, nat
                   cmu = cmu + CONJG(u(3 * (na - 1) + ipol, mu))
                ENDDO
                aux = aux + cmu * psi(:, :, mu, 1)
             ENDDO
             dvpsi = (0.d0, 0.d0)
             DO ib = 1, nocc
                DO ig = 1, npw
                   dvpsi(ig, ib) = CMPLX(0.d0, tpiba * (xk(ipol, ikk) + g(ipol, igk_k(ig, ikk))), KIND=DP) &
                                   * evc(ig, ib)
                ENDDO
             ENDDO
             CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq, .FALSE.)   ! -P_c d psi
             tres(1, ipol) = tres(1, ipol) + wk(ikk) * DBLE(dot(aux - dvpsi, aux - dvpsi, nocc))
             tres(2, ipol) = tres(2, ipol) + wk(ikk) * DBLE(dot(dvpsi, dvpsi, nocc))
          ENDDO
       ENDIF
    ENDDO
    CALL mp_sum(mat, intra_pool_comm)
    CALL mp_sum(mat, inter_pool_comm)
    CALL mp_sum(tres, intra_pool_comm)
    CALL mp_sum(tres, inter_pool_comm)
    CALL mp_sum(nfail, inter_pool_comm)
    DEALLOCATE (psi, h_diag, dvb, aux, dvloc)
    IF (nfail > 0) WRITE(stdout, '(5x,"inverse moments: ",i6," Sternheimer roots not converged to ",es9.2)') &
       nfail, thresh
    !
    ! Hermitian part (the k sum over the irreducible wedge is Hermitian only after symmetrization)
    DO i = 1, 6
       mat(:, :, i) = 0.5d0 * (mat(:, :, i) + CONJG(TRANSPOSE(mat(:, :, i))))
       CALL symdyn_munu_new(mat(:, :, i), u, xq, s, invs, rtau, irt, at, bg, nsymq, nat, irotmq, minus_q)
    ENDDO
    !
    ! summary: the trace over all Cartesian components and atoms, and at q = 0 the rigid translation
    DO i = 1, 6
       tr(i) = 0.d0
       DO mu = 1, nmodes
          tr(i) = tr(i) + DBLE(mat(mu, mu, i))
       ENDDO
    ENDDO
    WRITE(stdout, '(/,5x,"Inverse moments (Ry a.u.), trace over atoms and directions:")')
    WRITE(stdout, '(5x,"  scf vertex:  I2 ",es14.6,"  I3 ",es14.6,"  I4 ",es14.6)') tr(1:3)
    WRITE(stdout, '(5x,"  bare vertex: I2 ",es14.6,"  I3 ",es14.6,"  I4 ",es14.6)') tr(4:6)
    WRITE(stdout, '(5x,"  4*I3 per atom and direction (m_e): scf ",f12.5,"  bare ",f12.5)') &
       4.d0 * tr(2) / nmodes, 4.d0 * tr(5) / nmodes
    IF (lgamma) THEN
       DO ipol = 1, 3
          tr(1) = 0.d0
          tr(2) = 0.d0
          DO na = 1, nat
             DO nb = 1, nat
                tr(1) = tr(1) + DBLE(mat(3*(na-1)+ipol, 3*(nb-1)+ipol, 2))
                tr(2) = tr(2) + DBLE(mat(3*(na-1)+ipol, 3*(nb-1)+ipol, 5))
             ENDDO
          ENDDO
          WRITE(stdout, '(5x,"  rigid translation along ",i1,": 4*I3 = ",f12.5," (scf), ",f12.5, &
               &" (bare) m_e per cell (scf = ",f8.3," valence electrons for a local potential only)")') &
               ipol, 4.d0 * tr(1), 4.d0 * tr(2), nelec
          WRITE(stdout, '(5x,"  translation check |psi1 + P_c d psi| / |P_c d psi| = ",es10.3)') &
               SQRT(tres(1, ipol) / MAX(tres(2, ipol), 1.d-300))
       ENDDO
    ENDIF
    !
    CALL star_q1(xq, at, bg, nsym, s, invs, nq, sxq, isq, imq, .FALSE., t_rev)
    DO i = 1, 6
       iu = 6
       IF (ionode) THEN
          iu = find_free_unit()
          OPEN (UNIT=iu, FILE=TRIM(fildyn)//ext(i), STATUS='unknown', FORM='formatted')
          CALL write_old_dyn_mat_head(iu)
       ENDIF
       CALL write_star(mat(:, :, i), nq, sxq, isq, imq, iu)
       IF (ionode) CLOSE (UNIT=iu, STATUS='keep')
    ENDDO
    DEALLOCATE (mat)
    CALL cleanup()
    CALL stop_clock('imom')
    !
  CONTAINS
    !
    SUBROUTINE solve(x)
      !! x = (H - e_v + alpha P_v)^-1 (-P_c dvpsi) on the occupied bands of k
      COMPLEX(DP), INTENT(OUT) :: x(:,:)
      CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq, .FALSE.)
      dpsi = (0.d0, 0.d0)
      conv_root = .TRUE.
      CALL cgsolve_all(ch_psi_all, cg_psi, et(1, ikk), dvpsi, dpsi, h_diag, npwx, npwq, thresh, ik, &
                       lter, conv_root, anorm, nocc, npol)
      IF (.NOT. conv_root) nfail = nfail + 1
      x = dpsi
    END SUBROUTINE solve
    !
    COMPLEX(DP) FUNCTION dot(a, b, n)
      COMPLEX(DP), INTENT(IN) :: a(:,:), b(:,:)
      INTEGER, INTENT(IN) :: n
      COMPLEX(DP), EXTERNAL :: zdotc
      INTEGER :: jb
      dot = (0.d0, 0.d0)
      DO jb = 1, n
         dot = dot + zdotc(SIZE(a, 1), a(1, jb), 1, b(1, jb), 1)
      ENDDO
    END FUNCTION dot
    !
  END SUBROUTINE imom_finish
  !
  !--------------------------------------------------------------------------
  SUBROUTINE write_star(dyn, nq, sxq, isq, imq, iu)
    !------------------------------------------------------------------------
    !! The star of q of a matrix that transforms like the dynamical matrix (the rotation of q2qstar_ph),
    !! written in the layout of write_dyn_on_file with numbers in exponent format.
    USE ions_base,    ONLY : nat
    USE cell_base,    ONLY : at, bg
    USE symm_base,    ONLY : s, invs, irt, nsym, t_rev
    USE lr_symm_base, ONLY : rtau
    !
    COMPLEX(DP), INTENT(IN) :: dyn(3 * nat, 3 * nat)
    INTEGER, INTENT(IN) :: nq, isq(48), imq, iu
    REAL(DP), INTENT(IN) :: sxq(3, 48)
    COMPLEX(DP) :: phi(3, 3, nat, nat), phi2(3, 3, nat, nat)
    INTEGER :: na, nb, iq, nsq, isym
    !
    nsq = nsym / nq
    IF (nsq * nq /= nsym) CALL errore('write_star', 'wrong degeneracy', 1)
    CALL scompact_dyn(nat, dyn, phi)
    DO na = 1, nat
       DO nb = 1, nat
          CALL trntnsc(phi(1, 1, na, nb), at, bg, -1)
       ENDDO
    ENDDO
    IF (imq /= 0) THEN
       phi2 = (0.d0, 0.d0)
       isym = 1
       DO WHILE (isq(isym) /= imq)
          isym = isym + 1
       ENDDO
       CALL rotate_and_add_dyn_tr(phi, phi2, nat, isym, s, invs, t_rev, irt, rtau, sxq(1, imq))
       phi = 0.5d0 * (phi + CONJG(phi2))
    ENDIF
    DO iq = 1, nq
       phi2 = (0.d0, 0.d0)
       DO isym = 1, nsym
          IF (isq(isym) == iq) CALL rotate_and_add_dyn_tr(phi, phi2, nat, isym, s, invs, t_rev, irt, rtau, &
                                                          sxq(1, iq))
       ENDDO
       phi2 = phi2 / DBLE(nsq)
       DO na = 1, nat
          DO nb = 1, nat
             CALL trntnsc(phi2(1, 1, na, nb), at, bg, +1)
          ENDDO
       ENDDO
       CALL put(sxq(:, iq), phi2)
       IF (imq == 0) CALL put(-sxq(:, iq), CONJG(phi2))
    ENDDO
    !
  CONTAINS
    !
    SUBROUTINE put(xq, p)
      REAL(DP), INTENT(IN) :: xq(3)
      COMPLEX(DP), INTENT(IN) :: p(3, 3, nat, nat)
      INTEGER :: i, j, ka, kb
      WRITE(iu, '(/,5x,"Dynamical  Matrix in cartesian axes",//,5x,"q = ( ",3f14.9," ) ",/)') xq
      DO ka = 1, nat
         DO kb = 1, nat
            WRITE(iu, '(2i5)') ka, kb
            DO i = 1, 3
               WRITE(iu, '(3(2(es21.12e3,1x),2x))') (p(i, j, ka, kb), j = 1, 3)
            ENDDO
         ENDDO
      ENDDO
    END SUBROUTINE put
    !
  END SUBROUTINE write_star
  !
  !--------------------------------------------------------------------------
  SUBROUTINE cleanup()
    IF (ALLOCATED(dvall)) DEALLOCATE (dvall)
    IF (ALLOCATED(mode_done)) DEALLOCATE (mode_done)
  END SUBROUTINE cleanup
  !
END MODULE ph_imom
