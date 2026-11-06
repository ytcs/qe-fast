!
! Copyright (C) 2026 The qe-fast authors
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE qef_config
  !--------------------------------------------------------------------------
  !! Central switch table of qe-fast.
  !!
  !! Every run-time switch is read through qef_getenv. The value of a switch is, in order:
  !!   1. the environment variable of that name, when it is set;
  !!   2. the entry of the active preset in the table below;
  !!   3. unset (status 1): the code's built-in default for that switch.
  !! The preset comes from QEF_PRESET (off, safe, fast, max) or, when that is unset or empty, from
  !! default_preset. Any other value of QEF_PRESET stops the run.
  !!
  !! qefast/SWITCHES.md gives, per switch: default, class, gate evidence and fallback. A switch that
  !! no preset sets and that is not a lever in its own right (test knobs, tuning values) has no row
  !! here; it is still read through qef_getenv and is listed in SWITCHES.md.
  !!
  !! No dependencies, so that UtilXlib and the Wannier90 copy can use it too.
  !
  IMPLICIT NONE
  PRIVATE
  !
  PUBLIC :: qef_getenv, qef_preset_name, qef_report
  !
  CHARACTER(LEN=*), PARAMETER :: default_preset = 'safe'
  !
  INTEGER, PARAMETER :: npreset = 4
  CHARACTER(LEN=4), PARAMETER :: preset_names(npreset) = (/ 'off ', 'safe', 'fast', 'max ' /)
  !
  TYPE :: qef_switch
    CHARACTER(LEN=28) :: name
    CHARACTER(LEN=4)  :: val_off, val_safe, val_fast, val_max   ! blank: the preset leaves the switch unset
  END TYPE qef_switch
  !
  INTEGER, PARAMETER :: nswitch = 38
  TYPE(qef_switch), PARAMETER :: table(nswitch) = (/ &
    ! name                                off   safe  fast  max
    ! engine: `off` is the engine with every lever off (QEF_STH=0 is the plain QE loop)
    qef_switch('QEF_STH',                 '1',  '1',  '1',  '1'), &
    qef_switch('QEF_STRICT',              '0',  '0',  '0',  '0'), &
    ! levers of ph.x
    qef_switch('QEF_STH_DVB',             '0',  '1',  '1',  '1'), &
    qef_switch('QEF_STH_FFTDX',           ' ',  '1',  '1',  '1'), &
    qef_switch('QEF_STH_FDXDP',           ' ',  '1',  '1',  '1'), &
    qef_switch('QEF_STH_OBATCH',          ' ',  '1',  '1',  '1'), &
    qef_switch('QEF_STH_QCACHE',          ' ',  '1',  '1',  '1'), &
    qef_switch('QEF_NSCF_SYMROT',         ' ',  '1',  '1',  '1'), &
    qef_switch('QEF_PHQ_FUSE',            ' ',  '1',  '1',  '1'), &
    qef_switch('QEF_PHQ_VKB',             ' ',  '1',  '1',  '1'), &
    qef_switch('QEF_DVXC_GPU',            ' ',  '1',  '1',  '1'), &
    qef_switch('QEF_STH_KDEV',            ' ',  '2',  '2',  '2'), &
    qef_switch('QEF_STH_PSIR_PART',       ' ',  '1',  '1',  '1'), &
    qef_switch('QEF_STH_DRHO2',           ' ',  '1',  '1',  '1'), &
    ! levers of wannier90 and EPW
    qef_switch('QEF_W90_NT',              '1',  '8',  '8',  '8'), &
    qef_switch('QEF_W90_AMN',             ' ',  '1',  '1',  '1'), &
    qef_switch('QEF_W90_NOHR',            ' ',  '1',  '1',  '1'), &
    qef_switch('QEF_EPW_FDX',             ' ',  '1',  '1',  '1'), &
    qef_switch('QEF_EPW_BEC',             ' ',  '1',  '1',  '1'), &
    qef_switch('QEF_EPW_STAR',            ' ',  '1',  '1',  '1'), &
    qef_switch('QEF_EPW_B2W',             ' ',  '1',  '1',  '1'), &
    qef_switch('QEF_EPW_IDLEWAIT',        ' ',  '1',  '1',  '1'), &
    qef_switch('QEF_EPW_NOVME',           ' ',  '1',  '1',  '1'), &
    ! not in safe: they move the Wannier gauge, and lambda with it. In fast and max as in their gated definitions
    qef_switch('QEF_PW_SYMROT',           ' ',  ' ',  '1',  '1'), &
    qef_switch('QEF_W90_WNT',             ' ',  ' ',  '8',  '8'), &
    qef_switch('QEF_W90_MMN',             ' ',  ' ',  '1',  '1'), &
    ! fast only: looser inner threshold with the deflation preconditioner
    qef_switch('QEF_THRESH_L',            ' ',  ' ',  '10', ' '), &
    qef_switch('QEF_STH_DEFL',            ' ',  ' ',  '1',  ' '), &
    ! max only: CBSR, which keeps its automatic memory and payoff check
    qef_switch('QEF_STH_CBSR',            ' ',  ' ',  ' ',  '1'), &
    qef_switch('QEF_STH_CBSR_AUTO',       ' ',  ' ',  ' ',  '1'), &
    qef_switch('QEF_STH_CBSR_XCH',        ' ',  ' ',  ' ',  '1'), &
    qef_switch('QEF_STH_CBSR_HOST',       ' ',  ' ',  ' ',  '1'), &
    ! opt-in in every preset
    qef_switch('QEF_EPW_NOWRITE',         ' ',  ' ',  ' ',  ' '), &
    qef_switch('QEF_STH_MP',              ' ',  ' ',  ' ',  ' '), &
    qef_switch('QEF_TR2_SCALE',           ' ',  ' ',  ' ',  ' '), &
    qef_switch('QEF_STH_THRESH_SCALE',    ' ',  ' ',  ' ',  ' '), &
    qef_switch('QEF_STH_MEMLOG',          ' ',  ' ',  ' ',  ' '), &
    qef_switch('QEF_STH_TRUERES',         ' ',  ' ',  ' ',  ' ') /)
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE qef_getenv(name, value, status)
    !------------------------------------------------------------------------
    !! GET_ENVIRONMENT_VARIABLE with the table as a second source: same arguments and same status
    !! (0 = has a value, 1 = unset).
    CHARACTER(LEN=*), INTENT(IN)  :: name
    CHARACTER(LEN=*), INTENT(OUT) :: value
    INTEGER, OPTIONAL, INTENT(OUT) :: status
    INTEGER :: ios, i
    CHARACTER(LEN=4) :: v
    !
    CALL GET_ENVIRONMENT_VARIABLE(name, value, STATUS=ios)
    IF (ios == 1) THEN
      value = ' '
      DO i = 1, nswitch
        IF (table(i)%name /= name) CYCLE
        SELECT CASE (preset_index())
        CASE (1); v = table(i)%val_off
        CASE (2); v = table(i)%val_safe
        CASE (3); v = table(i)%val_fast
        CASE DEFAULT; v = table(i)%val_max
        END SELECT
        IF (LEN_TRIM(v) > 0) THEN
          value = v
          ios = 0
        ENDIF
        EXIT
      ENDDO
    ENDIF
    IF (PRESENT(status)) status = ios
  END SUBROUTINE qef_getenv
  !
  !--------------------------------------------------------------------------
  CHARACTER(LEN=4) FUNCTION qef_preset_name()
    !------------------------------------------------------------------------
    qef_preset_name = preset_names(preset_index())
  END FUNCTION qef_preset_name
  !
  !--------------------------------------------------------------------------
  SUBROUTINE qef_report(unit)
    !------------------------------------------------------------------------
    !! The active preset and the switches of the table that the environment overrides, for the run's output.
    INTEGER, INTENT(IN) :: unit
    CHARACTER(LEN=32) :: s
    INTEGER :: i, ios
    !
    CALL GET_ENVIRONMENT_VARIABLE('QEF_PRESET', s, STATUS=ios)
    IF (ios == 0 .AND. LEN_TRIM(s) > 0) THEN
      WRITE(unit, '(/5X,"qe-fast preset: ",A," (QEF_PRESET)")') TRIM(qef_preset_name())
    ELSE
      WRITE(unit, '(/5X,"qe-fast preset: ",A," (built-in default)")') TRIM(qef_preset_name())
    ENDIF
    DO i = 1, nswitch
      CALL GET_ENVIRONMENT_VARIABLE(TRIM(table(i)%name), s, STATUS=ios)
      IF (ios == 0) WRITE(unit, '(5X,"qe-fast override: ",A,"=",A)') TRIM(table(i)%name), TRIM(s)
    ENDDO
  END SUBROUTINE qef_report
  !
  !--------------------------------------------------------------------------
  INTEGER FUNCTION preset_index()
    !------------------------------------------------------------------------
    CHARACTER(LEN=32) :: s
    INTEGER :: ios, i
    !
    CALL GET_ENVIRONMENT_VARIABLE('QEF_PRESET', s, STATUS=ios)
    IF (ios /= 0 .OR. LEN_TRIM(s) == 0) s = default_preset
    DO i = 1, npreset
      IF (TRIM(s) == TRIM(preset_names(i))) THEN
        preset_index = i
        RETURN
      ENDIF
    ENDDO
    ERROR STOP 'qe-fast: QEF_PRESET must be off, safe, fast or max'
  END FUNCTION preset_index
  !
END MODULE qef_config
