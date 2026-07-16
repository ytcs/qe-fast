!
! Copyright (C) 2009 D. Forrer and M. Pavone
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Z=55-86 contributed by Martin Andersson (2011)
!------------------------------------------------------------------------------
!
MODULE london_module
  !
  !! Module for Dispersion Correction.  
  !! [ V. Barone et al. J. Comp. Chem., 30, 934 (2009) ]  
  !! [ S. Grimme, J. Comp. Chem., 27, 1787 (2006) ].
  !
  USE kinds ,           ONLY : DP
  USE parameters ,      ONLY : nsx
  !
  IMPLICIT NONE
  !
  SAVE
  !
  !
  REAL(DP) :: C6_i(nsx)
  !! C6\_i(ntyp) : atomic C6 coefficient of each atom type
  REAL ( DP ) , ALLOCATABLE :: R_vdw( : )
  !! R\_vdw(ntyp) : Van der Waals Radii of each atom type
  REAL ( DP ) , ALLOCATABLE :: C6_ij( : , : )
  !! C6\_ij(ntyp,ntyp) : C6 coefficients of each atom type pair:
  !! \(\sqrt{C6i\cdot C6j}\)
  REAL ( DP ) , ALLOCATABLE :: R_sum( : , : )
  !! R\_sum(ntyp,ntyp) : sum of VdW radii
  REAL ( DP ) , ALLOCATABLE :: r( : , : )
  !! r(3,mxr) : ordered distance vectors
  REAL ( DP ) , ALLOCATABLE :: dist2 ( : )
  !! dist2(mxr) : ordered distances
  !
  !
  REAL ( DP ) , PUBLIC :: scal6=0._dp
  !! global scaling factor
  REAL ( DP ) , PUBLIC :: lon_rcut=0._dp
  !! public cut-off radius
  REAL ( DP ) , PUBLIC :: in_C6 ( nsx ) = -1.0_dp
  !! in\_C6(ntyp) : input (user) specified atomic C6 coefficients
  REAL ( DP ) , PUBLIC :: in_rvdw( nsx ) = -1.0_dp
  !! in\_rvdw(ntyp) : input (user) specified atomic vdw radii
  !
  INTEGER :: mxr
  !! max number of r ( see rgen )
  !
  REAL ( DP ) :: r_cut
  !! cut-off radius in alat units
  REAL ( DP ) :: beta = 20.0_DP
  !! damping function parameter 
  !
  CONTAINS
   !
   !---------------------------------------------------------------------------
   ! Initialize parameters
   !---------------------------------------------------------------------------
   !
   SUBROUTINE init_london ( )
      !
      !! Extract parameters from database and compute C6\_ij and R\_sum(i,j).
      !
      USE ions_base ,          ONLY : ntyp => nsp, &
                                      atom_label => atm
      !
      USE cell_base ,          ONLY : alat, omega
      !
      USE constants,           ONLY : eps16
      !
      USE io_global,           ONLY : ionode, ionode_id, stdout
      !
      USE mp,                  ONLY : mp_bcast
      USE mp_images,           ONLY : intra_image_comm
      !
      INTEGER, PARAMETER :: maxZ = 86
      REAL (DP) :: vdw_coeffs(2,maxZ)
      !
      ! vdw C6 and radii for the first 86 atoms for the DFTD2 method
      ! Data from the DFT-D2 section of the dftd3.f file found on S.Grimme's home page:
      ! http://www.thch.uni-bonn.de/tc/index.php?section=downloads&subsection=DFT-D3
      ! See also S. Grimme, J. Comp. Chem., 27, 1787 (2006)
      ! First  column: C6, converted to Ry*Bohr^6 units
      ! (in the paper: J*nm^6/mol, conversion factor: 1 J*nm^6/mol = 34.69 Ry*Bohr^6)
      ! Second column: radii, in Bohr (in the paper they are in A)
      ! 
      DATA vdw_coeffs / &
         4.857_DP,    1.892_DP,&
         2.775_DP,    1.912_DP,&
        55.853_DP,    1.559_DP,&
        55.853_DP,    2.661_DP,&
       108.584_DP,    2.806_DP,&
        60.710_DP,    2.744_DP,&
        42.670_DP,    2.640_DP,&
        24.284_DP,    2.536_DP,&
        26.018_DP,    2.432_DP,&
        21.855_DP,    2.349_DP,&
       198.087_DP,    2.162_DP,&
       198.087_DP,    2.578_DP,&
       374.319_DP,    3.097_DP,&
       320.200_DP,    3.243_DP,&
       271.980_DP,    3.222_DP,&
       193.230_DP,    3.180_DP,&
       175.885_DP,    3.097_DP,&
       159.927_DP,    3.014_DP,&
       374.666_DP,    2.806_DP,&
       374.666_DP,    2.785_DP,&
       374.666_DP,    2.952_DP,&
       374.666_DP,    2.952_DP,&
       374.666_DP,    2.952_DP,&
       374.666_DP,    2.952_DP,&
       374.666_DP,    2.952_DP,&
       374.666_DP,    2.952_DP,&
       374.666_DP,    2.952_DP,&
       374.666_DP,    2.952_DP,&
       374.666_DP,    2.952_DP,&
       374.666_DP,    2.952_DP,&
       589.405_DP,    3.118_DP,&
       593.221_DP,    3.264_DP,&
       567.896_DP,    3.326_DP,&
       438.498_DP,    3.347_DP,&
       432.600_DP,    3.305_DP,&
       416.642_DP,    3.264_DP,&
       855.833_DP,    3.076_DP,&
       855.833_DP,    3.035_DP,&
       855.833_DP,    3.097_DP,&
       855.833_DP,    3.097_DP,&
       855.833_DP,    3.097_DP,&
       855.833_DP,    3.097_DP,&
       855.833_DP,    3.097_DP,&
       855.833_DP,    3.097_DP,&
       855.833_DP,    3.097_DP,&
       855.833_DP,    3.097_DP,&
       855.833_DP,    3.097_DP,&
       855.833_DP,    3.097_DP,&
      1294.678_DP,    3.160_DP,&
      1342.899_DP,    3.409_DP,&
      1333.532_DP,    3.555_DP,&
      1101.101_DP,    3.575_DP,&
      1092.775_DP,    3.575_DP,&
      1040.391_DP,    3.555_DP,&
     10937.246_DP,    3.405_DP,&
      7874.678_DP,    3.330_DP,&
      6114.381_DP,    3.251_DP,&
      4880.348_DP,    3.313_DP,&
      4880.348_DP,    3.313_DP,&
      4880.348_DP,    3.313_DP,&
      4880.348_DP,    3.313_DP,&
      4880.348_DP,    3.313_DP,&
      4880.348_DP,    3.313_DP,&
      4880.348_DP,    3.313_DP,&
      4880.348_DP,    3.313_DP,&
      4880.348_DP,    3.313_DP,&
      4880.348_DP,    3.313_DP,&
      4880.348_DP,    3.313_DP,&
      4880.348_DP,    3.313_DP,&
      4880.348_DP,    3.313_DP,&
      4880.348_DP,    3.313_DP,&
      3646.454_DP,    3.378_DP,&
      2818.308_DP,    3.349_DP,&
      2818.308_DP,    3.349_DP,&
      2818.308_DP,    3.349_DP,&
      2818.308_DP,    3.349_DP,&
      2818.308_DP,    3.349_DP,&
      2818.308_DP,    3.349_DP,&
      2818.308_DP,    3.349_DP,&
      1990.022_DP,    3.322_DP,&
      1986.206_DP,    3.752_DP,&
      2191.161_DP,    3.673_DP,&
      2204.274_DP,    3.586_DP,&
      1917.830_DP,    3.789_DP,&
      1983.327_DP,    3.762_DP,&
      1964.906_DP,    3.636_DP/
      !
      INTEGER :: ilab , ata , atb , i
      ! local : counter of atom type
      !         ata , atb : counters of C6_ij matrix
      !         counter
      INTEGER, EXTERNAL :: atomic_number
      !!
      REAL ( DP ) :: R_0, C_0, e_cut , sls
      ! local : buffers
      !
      IF ( ALLOCATED ( C6_ij ) ) RETURN
      !
      ! here we allocate parameters
      !
      ALLOCATE ( C6_ij ( ntyp , ntyp ) , &
                 R_sum ( ntyp , ntyp )   )
      !
      IF ( ionode ) THEN
         !
         ! and some buffers on ionode
         !
         IF (.NOT. ALLOCATED( R_vdw) ) ALLOCATE ( R_vdw ( ntyp )   )
         !
         ! here we initialize parameters to unphysical values
         !
         C6_i  ( : )     = -1.d0
         R_vdw ( : )     = -1.d0
         C6_ij ( : , : ) = -1.d0
         R_sum ( : , : ) = -1.d0
         !
         DO ilab = 1 , ntyp
           !
           i = atomic_number ( atom_label ( ilab ) )
           IF ( i > 0 .AND. i < 87 ) THEN
              IF ( in_C6 (ilab) > -eps16 ) THEN
                 C6_i  ( ilab )  = in_C6 (ilab)
              ELSE
                 C6_i  ( ilab )  = vdw_coeffs(1,i)
              END IF
              IF ( in_rvdw (ilab) > 0.0_DP ) THEN
                R_vdw ( ilab )  = in_rvdw (ilab)
              ELSE
                R_vdw ( ilab )  = vdw_coeffs(2,i)
              END IF
           ELSE
             CALL errore ( ' init_london ' ,&
                           'atom ' // atom_label(ilab) //' not found ' , ilab )
           END IF
           !
         END DO
         !
         ! are there all the parameters we need?
         !
         DO ilab = 1 , ntyp
           !
           IF ( ( C6_i  ( ilab ) < -eps16 ) .or. &
                ( R_vdw ( ilab ) < 0.0_DP ) ) THEN
             !
             CALL errore ( ' init_london ' ,&
                           ' one or more parameters not found ' , 4 )
             !
           END IF
           !
         END DO
         !
         ! ...here we store C6_ij parameters of each pair of atom types
         !    into a square matrix C6_ij = sqrt ( C6_i * C6_j )
         !
         DO atb = 1 , ntyp
           !
           DO ata = 1 , ntyp
             !
             C6_ij ( ata , atb ) = sqrt ( C6_i ( ata ) * C6_i ( atb ) )
             !
             R_sum ( ata , atb ) = R_vdw ( ata ) + R_vdw ( atb )
             !
           END DO
           !
         END DO
         !
         ! ... cutoff radius in alat units
         !
         r_cut = lon_rcut / alat
         !
         ! ... define a gross maximum bound of the mxr size
         !
         mxr = 1 + INT ( ( 2 * ( lon_rcut + alat ) ) ** 3 / omega )
         !
      END IF
      !
      ! broadcast data to all processors
      !
      CALL mp_bcast ( C6_ij,  ionode_id, intra_image_comm )
      CALL mp_bcast ( R_sum,  ionode_id, intra_image_comm )
      CALL mp_bcast ( r_cut,  ionode_id, intra_image_comm )
      CALL mp_bcast ( mxr  ,  ionode_id, intra_image_comm )
      !
      ALLOCATE ( r ( 3 , mxr ) , dist2 ( mxr ) )
      !
      RETURN
      !
   END SUBROUTINE init_london
   !
   SUBROUTINE print_london
      !
      USE io_global,           ONLY : ionode, stdout
      USE ions_base ,          ONLY : ntyp => nsp, atom_label => atm
      INTEGER :: ata
      !
      IF ( ionode .AND. ALLOCATED(R_vdw) ) THEN
         ! do not print anything if DFT-D2 variables not set
         WRITE ( stdout ,'( /, 5X, "-------------------------------------------------" , &
                          & /, 5X, "Parameters for Dispersion (Grimme-D2) Correction:" , &
                          & /, 5X, "-------------------------------------------------" , &
                          & /, 5X, "  atom      VdW radius       C_6     " , / )' )
         DO ata = 1 , ntyp
            !
            WRITE (stdout , '( 8X, A3 , 6X , F7.3 , 6X , F9.3 )' ) &
                               atom_label ( ata ) , R_vdw ( ata ) , C6_i ( ata )
            !
         END DO
      END IF
      !
   END SUBROUTINE print_london
   !
   !---------------------------------------------------------------------------
   ! Compute dispersion energy
   !---------------------------------------------------------------------------
   !
   FUNCTION energy_london ( alat , nat , ityp , at , bg , tau )
    !
    !! Here we compute the dispersion contribution to the total energy:
    !! $$ E = -( C_6^{ij} / R_{ij}^6 ) f_\text{damp}( R_{ij} ) \text{scal6} $$
    !! where f_damp is the damping function:
    !! $$ f_\text{damp}(R_{ij}) = [1+\exp{-\text{beta}(R_{ij}/(R_i^0+R_j^0)-1 )}]^{-1} $$
    !! and \(\text{scal6}\) is a global scaling factor.
    !
    USE mp_images,    ONLY : me_image , nproc_image, intra_image_comm
    USE mp,           ONLY : mp_sum
    !
    INTEGER , INTENT ( IN ) :: nat
    !! number of atoms
    INTEGER , INTENT ( IN ) :: ityp ( nat ) 
    !! type of each atom
    REAL ( DP ) , INTENT ( IN ) :: alat
    !! the cell parameter
    REAL ( DP ) , INTENT ( IN ) :: tau (3, nat)
    !! atomic positions in alat units
    REAL ( DP ) , INTENT ( IN ) :: at ( 3 , 3 )
    !! direct lattice vectors
    REAL ( DP ) , INTENT ( IN ) :: bg ( 3 , 3 )
    !! reciprocal lattice vectors
    !
    INTEGER :: ata , atb , nrm , nr
    ! locals : 
    ! ata , atb : atom counters
    ! nrm :       actual number of vectors computed by rgen
    ! nr :        counter
    !
    INTEGER :: na_s, na_e, mykey
    ! locals :    parallelization stuff
    !
    REAL ( DP ) :: dist , f_damp , energy_london , dtau ( 3 ) , dist6
    ! locals:
    ! dist          : distance R_ij between the current pair of atoms
    ! f_damp        : damping function
    ! energy_london : the dispersion energy
    ! dtau          : output of rgen ( not used )
    ! dist6         : distance**6
    !
    !
    CALL start_clock('energy_london')
    energy_london = 0.d0
    !
    ! poor-man parallelization over atoms
    ! - if nproc_image=1   : na_s=1, na_e=nat, mykey=0
    ! - if nproc_image<=nat: each proc. calculates atoms na_s to na_e; mykey=0
    ! - if nproc_image>nat : each processor takes care of atom na_s=na_e;
    !   mykey labels how many times each atom appears (mykey=0 first time etc.)
    !
    CALL block_distribute( nat, me_image, nproc_image, na_s, na_e, mykey )
    !
    ! ... the dispersion energy
    !
    IF ( mykey == 0 ) THEN
       !
       DO ata = na_s, na_e
          !
          DO atb = 1 , nat
             !
             dtau ( : ) = tau ( : , ata ) - tau ( : , atb )
             !
             CALL rgen ( dtau, r_cut, mxr, at, bg, r, dist2, nrm )
             !
 !$omp parallel do private(nr,dist,dist6,f_damp) default(shared), reduction(+:energy_london)
             DO nr = 1 , nrm
                !
                dist  = alat * sqrt ( dist2 ( nr ) )
                dist6 = beta*( dist / ( R_sum( ityp(atb), ityp(ata) ) ) - 1.0_dp )
                !
                ! dist6 is used here as temporary variable to avoid computing
                ! e^-x for too large x (for x=40, e^-x=4*10^-18)
                !
                IF ( dist6 < 40.0_dp ) THEN
                   f_damp = 1.0_dp / ( 1.d0 + exp ( -dist6 ) )
                ELSE
                   f_damp = 1.0_dp
                END IF
                !
                dist6 = dist**6
                energy_london = energy_london - &
                     ( C6_ij ( ityp ( atb ) , ityp ( ata ) ) / dist6 ) * &
                     f_damp
                !
             END DO
!$omp end parallel do
             !
          END DO
          !
       END DO
       !
       energy_london = scal6 * 0.5d0 * energy_london
       !
    ENDIF
    !
    CALL mp_sum ( energy_london , intra_image_comm )
    !
    CALL stop_clock('energy_london')
    RETURN
    !
   END FUNCTION energy_london
   !
   !---------------------------------------------------------------------------
   !
   FUNCTION force_london ( alat , nat , ityp , at , bg , tau )
    !
    !! Compute dispersion forces acting on atoms.
    !
    USE mp_images,    ONLY : me_image , nproc_image , intra_image_comm
    USE mp,           ONLY : mp_sum
    !
    INTEGER , INTENT ( IN ) :: nat
    !! number of atoms
    INTEGER , INTENT ( IN ) :: ityp ( nat )
    !! type of each atom
    REAL ( DP ) , INTENT ( IN ) :: alat
    !! the cell parameter
    REAL ( DP ) , INTENT ( IN ) :: tau (3, nat)
    !! atomic positions in alat units
    REAL ( DP ) , INTENT ( IN ) :: at ( 3 , 3 )
    !! direct lattice vectors
    REAL ( DP ) , INTENT ( IN ) :: bg ( 3 , 3 )
    !! reciprocal lattice vectors
    !
    INTEGER :: ata , atb , nrm , nr , ipol
    ! locals :
    ! ata , atb : atom counters
    ! nrm       : actual number of vectors computed by rgen
    ! nr        : counter on neighbours shells
    ! ipol      : counter on coords
    !
    INTEGER :: na_s, na_e, mykey
    ! locals :    parallelization stuff
    !
    REAL ( DP ) :: dist , f_damp , dtau ( 3 ) , force_london ( 3 , nat ) , &
                   dist6 , dist7 , exparg , expval , par , fac , add, aux(3)
    ! locals :
    ! dist         : distance R_ij between the current pair of atoms
    ! f_damp       :  damping function
    ! dtau         :  \vec R_ij
    ! force_london : dispersion forces
    ! dist6        :  dist**6
    ! dist7        :  dist**7
    ! ...  and some buffers
    !
    CALL start_clock('force_london')
    !
    ! parallelization: divide atoms across processors of this image
    ! (different images have different atomic positions)
    ! see energy_london for explanations
    !
    CALL block_distribute( nat, me_image, nproc_image, na_s, na_e, mykey )
    !
    ! ... the dispersion forces
    !
    force_london ( : , : ) = 0.d0
    !
    IF ( mykey == 0 ) THEN
       !
       DO ata = na_s, na_e
          !
          DO atb = 1 , nat
             !
             IF ( ata /= atb ) THEN
                !
                dtau ( : ) = tau ( : , ata ) - tau ( : , atb )
                !
                ! generate neighbours shells
                !
                CALL rgen ( dtau, r_cut, mxr, at, bg, r, dist2, nrm )
                !
                ! compute forces
                !
                par = beta / ( R_sum ( ityp ( atb ) , ityp ( ata ) ) )
                !
                aux(:) = 0.d0
!$omp parallel do private(nr,dist,dist6,dist7,exparg,expval,fac,add,ipol) default(shared), reduction(+:aux)
                DO nr = 1 , nrm
                   !
                   dist  = alat * sqrt ( dist2 ( nr ) )
                   dist6 = dist ** 6
                   dist7 = dist6 * dist
                   !
                   exparg = - beta * ( dist / ( R_sum ( ityp(atb) , ityp(ata) ) ) - 1 )
                   expval = exp ( exparg )
                   !
                   fac = C6_ij ( ityp ( atb ) , ityp ( ata ) ) / dist6
                   add = 6.d0 / dist
                   !
                   DO ipol = 1 , 3
                      ! workaround for OpenMP on Intel compilers
                      aux(ipol) = aux(ipol) + &
                           ( scal6 / ( 1 + expval ) * fac * &
                           ( - par * expval / ( 1.d0 + expval ) + add ) * &
                           r ( ipol , nr ) * alat / dist )
                      !
                   END DO
                   !
                END DO
!$omp end parallel do 
                DO ipol = 1 , 3
                   force_london ( ipol , ata ) = force_london ( ipol , ata ) + aux(ipol)
                ENDDO
                !
             END IF
             !
          END DO
          !
       END DO
       !
    END IF
    !
    CALL mp_sum ( force_london , intra_image_comm )
    !
    CALL stop_clock('force_london')
    RETURN
    !
   END FUNCTION force_london
   !
   !
   !---------------------------------------------------------------------------
   !
   FUNCTION stres_london ( alat , nat , ityp , at , bg , tau , omega )
    !
    !! Compute dispersion contribution to the stress tensor.
    !
    USE mp_images,    ONLY : me_image , nproc_image , intra_image_comm
    USE mp,           ONLY : mp_sum
    !
    INTEGER , INTENT ( IN ) :: nat
    !! number of atoms
    INTEGER , INTENT ( IN ) :: ityp ( nat ) 
    !! type of each atom
    REAL ( DP ) , INTENT ( IN ) :: alat
    !! the cell parameter
    REAL ( DP ) , INTENT ( IN ) :: tau (3, nat)
    !! atomic positions in alat units
    REAL ( DP ) , INTENT ( IN ) :: omega
    !! unit cell volume
    REAL ( DP ) , INTENT ( IN ) :: at ( 3 , 3 )
    !! direct lattice vectors
    REAL ( DP ) , INTENT ( IN ) :: bg ( 3 , 3 )
    !! reciprocal lattice vectors
    !
    INTEGER :: ata , atb , nrm , nr , ipol , lpol , spol
    ! locals :
    ! ata , atb : atom counters
    ! nrm       : actual number of vectors computed by rgen
    ! nr        : counter on neighbours shells
    ! xpol      : coords counters ipol lpol spol
    !
    INTEGER :: na_s, na_e, mykey
    ! locals :    parallelization stuff
    !
    REAL ( DP ) :: dist , f_damp , dtau ( 3 ) , stres_london ( 3 , 3 ) , &
                   dist6 , dist7 , exparg , expval , par , fac , add
    ! locals:
    ! dist : distance R_ij of current pair of atoms
    ! f_damp       : damping function
    ! dtau         : \vec R_ij
    ! stres_london : dispersion contribution to stress tensor
    ! dist6        : dist**6
    ! dist7        : dist**7
    !       and some buffers
    !
    CALL start_clock('stres_london')
    !
    ! parallelization: divide atoms across processors of this image
    ! (different images have different atomic positions)
    ! see energy_london for explanations
    !
    CALL block_distribute( nat, me_image, nproc_image, na_s, na_e, mykey )
    !
    ! ... the dispersion stress tensor
    !
    stres_london ( : , : ) = 0.d0
    !
    IF ( mykey == 0 ) THEN
       !
       DO ata = na_s, na_e
          !
          DO atb = 1 , nat
             !
             dtau ( : ) = tau ( : , ata ) - tau ( : , atb )
             !
             ! generate neighbours shells
             !
             CALL rgen ( dtau, r_cut, mxr, at, bg, r, dist2, nrm )
             !
             ! compute stress
             !
             par = beta / ( R_sum ( ityp ( atb ) , ityp ( ata ) ) )
             !
             DO nr = 1 , nrm
                !
                dist  = alat * sqrt ( dist2 ( nr ) )
                dist6 = dist ** 6
                dist7 = dist6 * dist
                !
                exparg = - beta * ( dist / ( R_sum ( ityp ( atb ) , ityp ( ata ) ) ) - 1 )
                !
                expval = exp ( exparg )
                !
                fac = C6_ij ( ityp ( atb ) , ityp ( ata ) ) / dist6
                !
                add = 6.d0 / dist
                !
                DO ipol = 1 , 3
                   !
                   DO lpol = 1 , ipol
                      !
                      stres_london ( lpol , ipol ) = stres_london ( lpol , ipol ) + &
                           ( scal6 / ( 1 + expval ) * fac * &
                           ( - par * expval / ( 1.d0 + expval ) + add ) * &
                           r ( ipol , nr ) * alat / dist ) * &
                           r ( lpol , nr ) * alat
                      !
                   END DO
                   !
                END DO
                !
             END DO
             !
          END DO
          !
       END DO
       !
    END IF
    !
    DO ipol = 1 , 3
       !
       DO lpol = ipol + 1 , 3
          !
          stres_london ( lpol , ipol ) = stres_london ( ipol , lpol )
          !
       END DO
       !
    END DO
    !
    stres_london ( : , : ) = - stres_london ( : , : ) / ( 2.d0 * omega )
    !
    CALL mp_sum ( stres_london , intra_image_comm )
    !
    CALL stop_clock('stres_london')
    RETURN
    !
   END FUNCTION stres_london
   !
   !---------------------------------------------------------------------------
   !
   SUBROUTINE dealloca_london
   !
   !! Clean memory.
   !
   IF ( ALLOCATED ( R_vdw ) ) DEALLOCATE ( R_vdw )
   IF ( ALLOCATED ( C6_ij ) ) DEALLOCATE ( C6_ij )
   IF ( ALLOCATED ( R_sum ) ) DEALLOCATE ( R_sum )
   IF ( ALLOCATED ( r     ) ) DEALLOCATE ( r     )
   IF ( ALLOCATED ( dist2 ) ) DEALLOCATE ( dist2 )
   !
   RETURN
   !
   END SUBROUTINE dealloca_london
   !
END MODULE london_module
