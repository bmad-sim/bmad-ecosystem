!+
! Module: cesr_mod
!
! This module holds all the CESR specific type definitions.
!-

#include "CESR_platform.inc"

module cesr_mod

  use precision_def

!-------------------------------------------------------------------------
! CESR logical names

  integer, parameter :: q49aw$ = 101, q47aw$ = 102
  integer, parameter :: q47ae$ = 103, q49ae$ = 104
  integer, parameter :: q43aw$ = 105, q43ae$ = 106
  integer, parameter :: q08aw$ = 107

  integer, parameter :: h_sep_08w$ = 1, h_sep_45w$ = 2
  integer, parameter :: h_sep_45e$ = 3, h_sep_08e$ = 4
  integer, parameter :: v_sep_48w$ = 5, v_sep_48e$ = 6

  integer, parameter :: wig_w$ = 1, wig_e$ = 2

  integer, parameter :: rf_w1$ = 1, rf_w2$ = 2, rf_e1$ = 3, rf_e2$ = 4

  integer, parameter :: scir_tilt_w$ = 1, scir_tilt_sk_w$ = 2
  integer, parameter :: scir_tilt_e$ = 3, scir_tilt_sk_e$ = 4

!---------------------------------------------------------------
! the cesr_struct is used by bmad_to_cesr.
! the element ordering here is similar to what is used by cesrv

  type b_struct
    real(rp) beta_mid             ! beta at the midpoint
    real(rp) beta_ave             ! beta averaged over the element
    real(rp) pos_mid              ! position at midpoint
    real(rp) pos_inj              ! position in injection lattice
  end type

  type cesr_element_struct 
    character*16 name              ! bmad name
    type (b_struct)  x, y          ! beta's and positions
    character*12 db_node_name
    integer ix_db                  ! element index for data base node
    integer ix_ring                ! index to element in ring structure
  end type

  integer, parameter :: n_quad_maxx = 120
  integer, parameter :: n_oct_maxx = 4
  integer, parameter :: n_hbnd_maxx = 6
  integer, parameter :: n_rf_maxx = 4
  integer, parameter :: n_skew_sex_maxx = 20
  integer, parameter :: n_wig_maxx = 2
  integer, parameter :: n_sep_maxx = 6
  integer, parameter :: n_qadd_maxx = 7
  integer, parameter :: n_scir_cam_maxx = 10
  integer, parameter :: n_scir_quad_maxx = 4
  integer, parameter :: n_scir_tilt_maxx = 8

! %ix_cesr is a pointer to cesr_struct for given type

  type cesr_struct
    type (cesr_element_struct) quad_(0:120)
    type (cesr_element_struct) skew_quad_(0:120)
    type (cesr_element_struct) sex_(0:120)
    type (cesr_element_struct) det_(0:120)
    type (cesr_element_struct) skew_sex_(n_skew_sex_maxx)
    type (cesr_element_struct) oct_(n_oct_maxx)
    type (cesr_element_struct) rf_(n_rf_maxx)
    type (cesr_element_struct) wig_(n_wig_maxx)
    type (cesr_element_struct) sep_(n_sep_maxx)
    type (cesr_element_struct) h_steer_(0:120)
    type (cesr_element_struct) v_steer_(0:120)
    type (cesr_element_struct) solenoid       ! solenoid struct
    type (cesr_element_struct) scir_cam_rho_(n_scir_cam_maxx)
    type (cesr_element_struct) scir_tilt_(n_scir_tilt_maxx)
    integer ix_ip_l3                     ! pointer to IP_L3
    integer, pointer :: ix_cesr(:) => null() 
  end type

!---------------------------------------------------------------------------
! DB_STRUCT:                       
! This structrue holds info on the correspondence between CESR data base
! elements and a BMAD ring. Use BMAD_TO_DB to initialize this structure.
! For completeness there are, in addition, arrays for the detectors and
! wigglers, etc.

! for an individual element

  type db_element_struct
    character*16 bmad_name    ! bmad name of element
    real(rp) dvar_dcu             ! calibration factor
    real(rp) var_theory           ! theory var value
    real(rp) var_0                ! extrapolated var value at CU = 0
    integer ix_ring           ! index to element array in ring struct
    integer ix_attrib         ! index to element attribute
    integer ix_cesrv          ! index to cesr_struct arrays
    character*12 db_node_name ! node name ("CSR QUAD CUR")
    integer ix_db             ! element index for data base node (5 for Q05W)
    character*16 db_ele_name  ! element name
    integer cu_high_lim       ! high limit
    integer cu_low_lim        ! low limit
    integer cu_now            ! current CU
    integer ix                ! Integer for general use.
  end type              

! for pointing to a db_struct array

  type db_node_struct
    type (db_element_struct), pointer :: ptr(:) => null()
  end type
                                                       
! db_struct def

  integer, parameter :: n_csr_sqewsext_maxx = 4

  type db_struct
    type (db_element_struct) :: csr_quad_cur(98)
    type (db_element_struct) :: csr_qadd_cur(n_qadd_maxx)
    type (db_element_struct) :: csr_horz_cur(98)
    type (db_element_struct) :: csr_hbnd_cur(n_hbnd_maxx)
    type (db_element_struct) :: csr_vert_cur(98)
    type (db_element_struct) :: csr_hsp_volt(n_sep_maxx)
    type (db_element_struct) :: csr_vsp_volt(n_sep_maxx)
    type (db_element_struct) :: csr_sext_cur(98)
    type (db_element_struct) :: csr_octu_cur(n_oct_maxx)
    type (db_element_struct) :: csr_sqewquad(98)
    type (db_element_struct) :: csr_sqewsext(n_csr_sqewsext_maxx)
    type (db_element_struct) :: scir_quadcur(n_scir_quad_maxx)
    type (db_element_struct) :: scir_skqucur(n_scir_quad_maxx)
    type (db_element_struct) :: scir_vertcur(n_scir_quad_maxx)
    type (db_element_struct) :: scir_sksxcur(1)
    type (db_node_struct) :: node(15)  ! does not include stuff below
! in db but without corresponding BMAD element
    type (db_element_struct) :: scir_pos_stp(n_scir_cam_maxx)
    type (db_element_struct) :: scir_enc_cnt(n_scir_cam_maxx)
    type (db_element_struct) :: scir_pos_rd(3*n_scir_cam_maxx)
! non data base stuff
    type (db_element_struct) :: quad(0:120) ! combinded csr_quad_cur, etc.
    type (db_element_struct) :: detector(0:120)
    type (db_element_struct) :: wiggler(1:n_wig_maxx)
    type (db_element_struct) :: scir_cam_rho(n_scir_cam_maxx)
    type (db_element_struct) :: scir_tilt(n_scir_tilt_maxx)
  end type

! for the synchrotron

  type synch_db_struct
    type (db_element_struct) :: bend(1:100)
    type (db_element_struct) :: detector(1:200)
    type (db_element_struct) :: hkicker(1:200)
    type (db_element_struct) :: vkicker(1:200)
    type (db_element_struct) :: htable(1:200)
    type (db_element_struct) :: vtable(1:200)
  end type

!-----------------------------------------------------------------------------
! For butns.nnnnn files

  type detector_struct
    real(rp) x_orb, y_orb
    integer amp(4)
    integer type
    logical ok
  endtype

  type butns_struct
    character*40 lattice
    character*20 date
    character*72 comment(5)
    type (detector_struct) det(0:120)
    integer save_set
    integer file_num
    integer turn    ! turn number for injection data
  end type

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine IDENTIFY_DB_NODE (DB_NAME, DB, DP_PTR, OK, TYPE_ERR)
!
! Subroutine to find which array in DB is associated with DB_NAME.
!
! Modules Needed:
!   use bmad
!
! Input:
!   DB_NAME -- Character*12: Data base name (eg. "CSR HSP VOLT")
!   DB      -- Db_node_struct: Data base structure.
!  
! Output:
!   DB_PTR   -- Db_element_struct, pointer(:) : Pointer to, eg., DB%HSP_VOLT(:)
!   OK       -- Logical: Set True if everything ok
!   TYPE_ERR -- Logical: If True then error message will be typed if needed.
!-

subroutine identify_db_node (db_name, db, db_ptr, ok, type_err)

  implicit none

  type (db_struct), target :: db
  type (db_element_struct), pointer :: db_ptr(:)

  integer k

  character*(*) db_name

  logical ok, type_err

!

  do k = 1, size(db%node)
    if (db%node(k)%ptr(1)%db_node_name /= db_name) cycle
    db_ptr => db%node(k)%ptr
    ok = .true.
    return
  enddo

  if (type_err) print *, &
        'ERROR IN IDENTIFY_DB_NODE: CANNOT FIND DATABASE NODE: ', db_name
  ok = .false.
  return

end subroutine

end module
