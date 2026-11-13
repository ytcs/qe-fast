! Copyright (C) 2026 The qe-fast authors.
! SPDX-License-Identifier: GPL-2.0-or-later
!
! qe-fast: external procedures that the d3q plugin calls and Quantum ESPRESSO 7.6 no longer provides
! as externals. Linked into qe_d3q by external/CMakeLists.txt; d3q itself is unmodified.
!
LOGICAL FUNCTION matches(string1, string2)
  ! QE 7.6 moved it into the upf_utils module.
  USE upf_utils, ONLY : upf_matches => matches
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: string1, string2
  matches = upf_matches(string1, string2)
END FUNCTION matches
!
LOGICAL FUNCTION imatches(string1, string2)
  ! Case-insensitive variant, also moved into upf_utils.
  USE upf_utils, ONLY : upf_imatches => imatches
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: string1, string2
  imatches = upf_imatches(string1, string2)
END FUNCTION imatches
!
SUBROUTINE sgam_ph(at, bg, nsym, s, irt, tau, rtau, nat, sym)
  ! QE 7.6 replaced sgam_ph by sgam_lr, which computes rtau for every one of the nsym operations.
  ! d3q passes sym = .true. for all of them, so the two agree; sym is accepted and not used.
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nsym, s(3,3,48), nat, irt(48,nat)
  REAL(DP), INTENT(IN) :: at(3,3), bg(3,3), tau(3,nat)
  LOGICAL, INTENT(IN) :: sym(nsym)
  REAL(DP), INTENT(OUT) :: rtau(3,48,nat)
  CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
END SUBROUTINE sgam_ph
