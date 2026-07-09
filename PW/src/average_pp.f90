!
! Copyright (C) 2005-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE average_pp( ntyp )
  !----------------------------------------------------------------------------
  !! Spin-orbit pseudopotentials transformed into standard pseudopotentials.
  !
  USE kinds,            ONLY : DP
  USE atom,             ONLY : rgrid
  USE uspp_param,       ONLY : upf
  USE ldaU,             ONLY : lda_plus_u
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ntyp
  !! number of species
  !
  ! ... local variables
  !
  INTEGER :: nt, nb, nbe, ind, ind1, l
  REAL(DP) :: weight
  !
  !
  DO nt = 1, ntyp
     !
     IF ( upf(nt)%has_so ) THEN
        !
        IF ( upf(nt)%tvanp .or. lda_plus_u) &
                CALL errore( 'average_pp', 'Fully relativistic PPs, need spin-orbit calc. (lspinorb=.true.)', 1 )
        !
        DO nb = 1, upf(nt)%nbeta
           !
           l = upf(nt)%lll(nb)
           !
           IF ( l == 0 ) THEN
              !
              ! ... s channel: only j = 1/2 exists, weight = 1, nothing to do
              !
              weight = 1.0_DP
              !
           ELSE IF ( upf(nt)%jjj(nb) > REAL( l, DP ) ) THEN
              !
              ! ... j = l + 1/2  ->  (2j+1)/[2(2l+1)] = (l+1)/(2l+1)
              !
              weight = ( l + 1.0_DP ) / ( 2.0_DP * l + 1.0_DP )
              !
           ELSE
              !
              ! ... j = l - 1/2  ->  (2j+1)/[2(2l+1)] = l/(2l+1)
              !
              weight = l / ( 2.0_DP * l + 1.0_DP )
              !
           ENDIF
           !
           upf(nt)%dion(nb,nb) = weight * upf(nt)%dion(nb,nb)
           !
        ENDDO
        !
        nbe = 0
        !
        DO nb = 1, upf(nt)%nwfc
           !
           nbe = nbe + 1
           !
           IF ( upf(nt)%lchi(nb) /= 0 .AND. &
                ABS(upf(nt)%jchi(nb)-upf(nt)%lchi(nb)-0.5D0 ) < 1.D-7 ) &
              nbe = nbe - 1
           !
        ENDDO
        !
        upf(nt)%nwfc = nbe
        !
        nbe = 0
        !
        DO nb = 1, upf(nt)%nwfc
           !
           nbe = nbe + 1
           !
           l = upf(nt)%lchi(nbe)
           !
           IF ( l /= 0 ) THEN
              !
              IF (ABS(upf(nt)%jchi(nbe)-upf(nt)%lchi(nbe)+0.5d0) < 1.d-7) THEN
                 IF ( ABS(upf(nt)%jchi(nbe+1)-upf(nt)%lchi(nbe+1)-0.5d0) > &
                      1.d-7) CALL errore( 'average_pp', 'wrong chi functions', 3 )
                 ind = nbe+1
                 ind1 = nbe
              ELSE
                 IF ( ABS(upf(nt)%jchi(nbe+1)-upf(nt)%lchi(nbe+1)+0.5d0) > &
                      1.d-7) CALL errore( 'average_pp', 'wrong chi functions', 4 )
                 ind = nbe
                 ind1 = nbe+1
              ENDIF
              !
              upf(nt)%chi(1:rgrid(nt)%mesh,nb) = &
                 ((l+1.D0) * upf(nt)%chi(1:rgrid(nt)%mesh,ind)+ &
                   l * upf(nt)%chi(1:rgrid(nt)%mesh,ind1)) / ( 2.D0 * l + 1.D0 )
              !
              nbe = nbe + 1
              !
           ELSE
              !
              upf(nt)%chi(1:rgrid(nt)%mesh,nb) = upf(nt)%chi(1:rgrid(nt)%mesh,nbe)
              !
           ENDIF
           !
           upf(nt)%lchi(nb) = upf(nt)%lchi(nbe)
           !
        ENDDO
        !
     ENDIF
     !
     upf(nt)%has_so = .FALSE.
     !
  ENDDO
  !
END SUBROUTINE average_pp
