program test_diagh_2

    USE laxlib_parallel_include
    USE mp,            ONLY : mp_bcast
    USE mp_world,      ONLY : mp_world_start, mp_world_end, mpime, &
                              root, world_comm
    USE mp_bands_util, ONLY : me_bgrp, root_bgrp, intra_bgrp_comm
    USE tester
    USE test_helpers,  ONLY : hermitian, symmetric, &
                              verify_generalized_eigenpairs
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    !
    TYPE(tester_t) :: test
    INTEGER :: world_group = 0
    !
    CALL test%init()
    !
#if defined(__MPI)
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    !
    me_bgrp = mpime; root_bgrp=root; intra_bgrp_comm=world_comm
    !
    CALL complex_1(test)
    !
    CALL real_1(test)
    !
    CALL collect_results(test)
    !
    CALL mp_world_end()
    !
    IF (mpime .eq. 0) CALL test%print()
    !
  CONTAINS
  !
  SUBROUTINE complex_1(test)
    USE LAXlib
    implicit none
    !
    TYPE(tester_t) :: test
    !
    integer, parameter :: n = 256
    !! size of the matrix. Must stay above the LAPACK tridiagonalization block
    !! size, or the m < n path below never reaches the blocked
    !! back-transformation of the eigenvectors. A few hundred is enough and
    !! keeps the test a few seconds.
    integer, parameter :: m = n/2
    !! number of eigenvalues and eigenvectors to be computed, must be <= n
    complex(DP) :: h(n,n)
    complex(DP) :: h_save(n,n)
    real(DP)    :: e(n)
    complex(DP) :: v(n,n)
    real(DP)    :: e_save(n)
    complex(DP) :: s(n,n)
    real(DP)    :: max_residual, max_ortho
    integer :: j
    !
    CALL hermitian(n, h)
    !
    h_save = h
    !
    v = (0.d0, 0.d0)
    e = 0.d0
    !
    CALL diagh(  n, n, h, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    DO j = 1, n
       CALL test%assert_close( h(1:n, j), h_save(1:n, j))
    END DO
    !
    e_save = e
    !
    ! Test that calling again gives the same results
    v = (0.d0, 0.d0)
    e = 0.d0
    CALL diagh(  n, n, h, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    DO j = 1, n
       CALL test%assert_close( h(1:n, j), h_save(1:n, j))
    END DO
    !
    test%tolerance32=1.e-5
    test%tolerance64=1.d-14
    CALL test%assert_close( e, e_save)
    !
    ! Test subset of eigenvalues
    v = (0.d0, 0.d0)
    e = 0.d0
    CALL diagh(  n, m, h, e, v(:,1:m), me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    DO j = 1, n
       CALL test%assert_close( h(1:n, j), h_save(1:n, j))
    END DO
    !
    test%tolerance32=1.e-5
    ! Use a looser tolerance for subset eigenvalues since ZHEEVD and ZHEEVX
    ! use different algorithms and may produce slightly different results
    test%tolerance64=1.d-12
    CALL test%assert_close( e(1:m), e_save(1:m))
    !
    ! Check eigenvectors. s = identity reduces verify_generalized_eigenpairs
    ! to the standard problem H v = e v.
    !
    s = (0.d0, 0.d0)
    DO j = 1, n
       s(j,j) = (1.d0, 0.d0)
    END DO
    !
    test%tolerance64=1.d-10
    IF (me_bgrp == root_bgrp) THEN
       CALL verify_generalized_eigenpairs( n, m, h_save, s, n, e, v, &
                                           max_residual, max_ortho )
       CALL test%assert_close( max_residual, 0.d0 )
       CALL test%assert_close( max_ortho, 0.d0 )
    END IF
    !
  END SUBROUTINE complex_1
  !
  SUBROUTINE real_1(test)
    USE LAXlib
    implicit none
    !
    TYPE(tester_t) :: test
    !
    integer, parameter :: n = 256
    !! size of the matrix; see the note in complex_1
    integer, parameter :: m = n/2
    !! number of eigenvalues and eigenvectors to be computed, must be <= n
    real(DP) :: h(n,n)
    real(DP) :: h_save(n,n)
    real(DP) :: e(n)
    real(DP) :: v(n,n)
    real(DP) :: e_save(n)
    real(DP) :: s(n,n)
    real(DP) :: max_residual, max_ortho
    integer :: j
    !
    CALL symmetric(n, h)
    !
    h_save = h
    !
    v = 0.d0
    e = 0.d0
    !
    CALL diagh(  n, n, h, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    DO j = 1, n
       CALL test%assert_close( h(1:n, j), h_save(1:n, j))
    END DO
    !
    e_save = e
    !
    ! Test that calling again gives the same results
    v = 0.d0
    e = 0.d0
    CALL diagh(  n, n, h, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    DO j = 1, n
       CALL test%assert_close( h(1:n, j), h_save(1:n, j))
    END DO
    !
    test%tolerance32=1.e-5
    test%tolerance64=1.d-14
    CALL test%assert_close( e, e_save)
    !
    ! Test subset of eigenvalues
    v = 0.d0
    e = 0.d0
    CALL diagh(  n, m, h, e, v(:,1:m), me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    DO j = 1, n
       CALL test%assert_close( h(1:n, j), h_save(1:n, j))
    END DO
    !
    test%tolerance32=1.e-5
    ! Use a looser tolerance for subset eigenvalues since DSYEVD and DSYEVX
    ! use different algorithms and may produce slightly different results
    test%tolerance64=1.d-12
    CALL test%assert_close( e(1:m), e_save(1:m))
    !
    ! Check eigenvectors; see complex_1
    !
    s = 0.d0
    DO j = 1, n
       s(j,j) = 1.d0
    END DO
    !
    test%tolerance64=1.d-10
    IF (me_bgrp == root_bgrp) THEN
       CALL verify_generalized_eigenpairs( n, m, h_save, s, n, e, v, &
                                           max_residual, max_ortho )
       CALL test%assert_close( max_residual, 0.d0 )
       CALL test%assert_close( max_ortho, 0.d0 )
    END IF
    !
  END SUBROUTINE real_1
  !
end program test_diagh_2
