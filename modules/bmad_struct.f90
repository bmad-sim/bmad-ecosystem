!+
! BMAD_STRUCT holds the structure definitions for BMAD routines.
!-

module bmad_struct

  use physical_constants

! The "regular" elements are in positions: 1 to RING.N_ELE_RING
! regular elements are:
!     1) Superimpose slaves:      SUPER_SLAVE$
!     2) Overlay slaves:          OVERLAY_SLAVE$
!     3) Container slaves:        CONTAINER_SLAVE$
!     5) Free:                    FREE$
!
! The "control" elements are in positions: RING.N_ELE_RING+1 to RING.N_ELE_MAX
!
! The control elements are:
!     1) Superimposed controllers:      SUPER_LORD$
!     2) Overlay controllers:           OVERLAY_LORD$
!     3) Component:                     COMPONENT_LORD$
!     5) Group controllers:             GROUP_LORD$
!-

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! IF YOU CHANGE THE RING STRUCTURE YOU MUST INCREASE THE VERSION NUMBER !
!
  integer, parameter :: bmad_inc_version$ = 44
!
! THIS IS USED BY BMAD_PARSER TO MAKE SURE DIGESTED FILES ARE OK.
!
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

! parameter def

  integer, parameter :: n_attrib_maxx = 44         ! \ max ix for attributes
  integer, parameter :: n_attrib_special_maxx = 50 ! / max ix including s$ etc.
  integer, parameter :: n_ele_maxx = 2200
  integer, parameter :: n_control_maxx = 5000
  integer, parameter :: n_comp_maxx = 50   ! Max number of components in a COIL

  integer, parameter :: n_key = 26, n_pole_maxx = 9

! Structure definitions

  type twiss_struct
    real beta, alpha, gamma, phi, eta, etap
    real mobius_beta, mobius_eta   ! Mobius effective beta and eta
    real sigma
  end type

  type ele_struct
    character*16 name              ! name of element
    integer key                    ! key value
    character*16 type              ! type name (see MAD documentation)
    character*16 alias             ! Another name
    real value(n_attrib_maxx)      ! attribute values
    real mat6(6,6)                 ! transpot matrix 
    type (twiss_struct)  x,y,z     ! Twiss parameters at end of element
    real c_mat(2,2)                ! 2x2 C coupling matrix
    real gamma_c                   ! gamma associated with C matrix
    real s                         ! longitudinal position at the end
    real x_position, y_position    ! Floor position of element
    real z_position                ! Elevation of element
    real theta_position            ! Floor orientation angle of element
    real phi_position              ! Angle of attack
    integer control_type           ! SUPER_SLAVE$, OVERLAY_LORD$, etc.
    character*16 attribute_name    ! Used by overlays
    integer ix_value               ! Pointer for attribute to control
    integer n_slave                ! Number of slaves
    integer ix1_slave              ! Start index for slave elements
    integer ix2_slave              ! Stop  index for slave elements
    integer n_lord                 ! Number of lords
    integer ic1_lord               ! Start index for lord elements
    integer ic2_lord               ! Stop  index for lord elements
    integer ix_pointer             ! Pointer for general use
    integer ixx                    ! Pointer for BMAD internal use
    integer iyy                    ! Pointer for BMAD internal use
    logical coupled                ! Is transport matrix coupled?
    logical mode_flip              ! Have the normal modes traded places?
    logical is_on                  ! For turning element on/off
    logical multipoles_on          ! For turning multipoles on/off
    logical nonzero_multipoles     ! 
  end type

  type control_struct
    integer ix_lord                ! index to lord element
    integer ix_slave               ! index to slave element
    integer ix_attrib              ! index of attribute controlled
    real coef                      ! control coefficient
  end type

! parameter and mode structures

  type param_struct
    integer particle              ! +1 = positrons, -1 = electrons
    real energy                   ! beam energy in GeV.
    real n_part                   ! Number of particles in a bunch
    integer symmetry              ! symmetry of the ring (e/w symm, etc.)
    logical z_decoupled           ! is z motion decoupled from x and y?
    real total_length             ! total_length of ring
    logical stable                ! is closed ring stable?
    real growth_rate              ! growth rate/turn if not stable
    logical aperture_limit_on     ! use apertures in tracking?
    logical lost                  ! for use in tracking
    integer ix_lost               ! If lost at what element?
    integer lattice_type          ! linac_lattice$, etc...
    integer ixx                   ! Integer for general use
  end type

  type mode_info_struct
    real tune      ! "fractional" tune in radians: 0 < tune < 2pi
    real emit      ! Emittance
    real chrom     ! Chromaticity
  end type

  type dummy_parameter_struct
    integer dummy(200)
  end type

! RING_STRUCT

  type ring_struct
    union
      map
        integer version               ! Version number
        character*16 name             ! Name of ring given by USE statement
        character*40 lattice          ! Lattice
        character*80 input_file_name  ! name of the lattice input file
        integer n_ele_ring            ! number of regular ring elements
        integer n_ele_symm            ! symmetry point for rings w/ symmetry
        integer n_ele_use             ! number of elements used
        integer n_ele_max             ! Index of last element used
        integer n_control_array       ! last index used in CONTROL_ array
        integer n_ic_array            ! last index used in IC_ array
        type (mode_info_struct)  x, y, z  ! tunes, etc.
        type (param_struct)      param ! parameters
      endmap
      map
        type (dummy_parameter_struct) parameters
      endmap
    endunion

    type (ele_struct)  ele_(0:n_ele_maxx)    ! Array of ring elements
    type (ele_struct)  ele_init              ! For use by any program
    type (control_struct)  control_(n_control_maxx)  ! control list
    integer ic_(n_control_maxx)                      ! index to %control_(:)
  end type

! structure for a particle

  type pos_vel_struct
    real pos, vel                            ! position and velocity
  end type

  type coord_struct
    union
      map
        type (pos_vel_struct)  x, y, z
      endmap
      map
        real vec(6)
      endmap
    endunion
  end type

! KEY value definitions
! Note: Null_element must be the last element in the list.
! Note: overlay$ == overlay_lord$ 

  integer, parameter :: drift$ = 1, sbend$ = 2, quadrupole$ = 3, group$ = 4
  integer, parameter :: sextupole$ = 5, overlay$ = 6, custom$ = 7
  integer, parameter :: solenoid$ = 8, rfcavity$ = 9
  integer, parameter :: elseparator$ = 10, beambeam$ = 11, wiggler$ = 12
  integer, parameter :: sol_quad$ = 13, marker$ = 14, kicker$ = 15
  integer, parameter :: hybrid$ = 16, octupole$ = 17, rbend$ = 18
  integer, parameter :: multipole$ = 19, coil$ = 20, loop$ = 21
  integer, parameter :: accel_sol$ = 22, define_energy$ = 23
  integer, parameter :: def_beam$ = 24, ab_multipole$ = 25
  integer, parameter :: null_ele$ = 27

  character*16 :: key_name(n_key+1) = (/ &
      'DRIFT        ', 'SBEND        ', 'QUADRUPOLE   ', 'GROUP        ', &
      'SEXTUPOLE    ', 'OVERLAY      ', 'CUSTOM       ', 'SOLENOID     ', &
      'RFCAVITY     ', 'ELSEPARATOR  ', 'BEAMBEAM     ', 'WIGGLER      ', &
      'SOL_QUAD     ', 'MARKER       ', 'KICKER       ', 'HYBRID       ', &
      'OCTUPOLE     ', 'RBEND        ', 'MULTIPOLE    ', 'COIL         ', &
      'LOOP         ', 'ACCEL_SOL    ', 'DEFINE_ENERGY', 'DEF_BEAM     ', &
      'AB_MULTIPOLE ', 'NULL_ELEMENT ', '             ' /)

! Attribute name logical definitions
! Note: The following attributes must have unique number assignments:
!     X_OFFSET, Y_OFFSET, X_PITCH, Y_PITCH

  integer, parameter :: x_beg_limit$=2, y_beg_limit$=3, b_x2$=4, &
          b_y2$=5, l_st2$=9, b_z$=10, l_st1$=11, s_st2$=12, s_st1$=13, &
          b_x1$=14, b_y1$=15    

  integer, parameter :: val1$=2, val2$=3, val3$=4, val4$=5, val5$=6, &
          val6$=7, val7$=8, val8$=9, val9$=10, val10$=11, val11$=12, &
          val12$=13, val13$=14, val14$=15, val15$=16, val16$=17

  integer, parameter :: l$=1
  integer, parameter :: tilt$=2, command$=2
  integer, parameter :: old_command$=3, angle$=3
  integer, parameter :: k1$=4, sig_x$=4, harmon$=4, h_displace$=4
  integer, parameter :: k2$=5, sig_y$=5, b_max$=5, v_displace$=5
  integer, parameter :: k3$=6, sig_z$=6, rf_wavelength$=6, rho$=6
  integer, parameter :: ks$=7, volt$=7, e1$=7, n_pole$=7, bbi_const$=7
  integer, parameter :: lag$=8, e2$=8, charge$=8, gap$=8, ap2$=8
  integer, parameter :: n_slice$=9, l_cord$=9, bp2$=9
  integer, parameter :: new_energy$=10, ap4$=10
  integer, parameter :: r2$=11, bp4$=11
  integer, parameter :: ri$=12, ap6$=12
  integer, parameter :: coef$=13, current$=13, bp6$=13
  integer, parameter :: roll$=14, r2i$=14, ap8$=14
  integer, parameter :: diameter$=15, bp8$=15
  integer, parameter :: x_pitch$=16
  integer, parameter :: y_pitch$=17
  integer, parameter :: hkick$=18
  integer, parameter :: vkick$=19
  integer, parameter :: x_offset$=20
  integer, parameter :: y_offset$=21
  integer, parameter :: x_limit$=22 
  integer, parameter :: y_limit$=23
  integer, parameter :: radius$=24

  integer, parameter :: a0$=25, b0$=26, a1$=27, b1$=28, a2$=29, b2$=30, &
             a3$=31, b3$=32, a4$=33, b4$=34, a5$=35, b5$=36, a6$=37, &
             b6$=38, a7$=39, b7$=40, a8$=41, b8$=42, a9$=43, b9$=44

  integer, parameter :: ix1_m$ = a0$, ix2_m$ = b9$

  integer, parameter :: k0l$=25, t0$=26, k1l$=27, t1$=28, k2l$=29, t2$=30, &
              k3l$=31, t3$=32, k4l$=33, t4$=34, k5l$=35, t5$=36, k6l$=37, &
              t6$=38, k7l$=39, t7$=40, k8l$=41, t8$=42, k9l$=43, t9$=44 

  integer, parameter :: type$ = 45   ! this is 1 greater than n_attrib_maxx
  integer, parameter :: alias$ = 46 
  integer, parameter :: start_edge$ = 47     ! special for groups
  integer, parameter :: end_edge$ = 48       ! special for groups
  integer, parameter :: s$ = 49              ! special for groups
  integer, parameter :: accordian_edge$ = 50 ! special for groups

  integer, parameter :: particle$ = 1
  integer, parameter :: energy$   = 2
  integer, parameter :: n_part$    = 3

  character*16 null_name / 'NULL' /

! electron/positron

  integer, parameter :: positron$ = +1
  integer, parameter :: electron$ = -1

  character*8 :: particle_name(-1:1) = (/ 'ELECTRON', '???', 'POSITRON' /)

! SYMMETRY etc., logical names

  integer, parameter :: no_symmetry$ = 0
  integer, parameter :: mobius_symmetry$ = 2
  integer, parameter :: ew_antisymmetry$ = -1  ! coupling elements are antisymmetric
  integer, parameter :: ew_symmetry$ = 3
  integer, parameter :: linac_lattice$ = 10
  integer, parameter :: linear_lattice$ = 11
  integer, parameter :: circular_lattice$ = 12

! logicals for MAKE_HYBIRD_RING

  logical, parameter :: z_decoupled$ = .true.,    no_z_decoupled$ = .false.
  logical, parameter :: remove_markers$ = .true., no_remove_markers$ = .false.

! control element logicals
! Note: overlay_lord$ == overlay$ !!

  integer, parameter :: free$ = 1, super_slave$ = 2, overlay_slave$ = 3
  integer, parameter :: group_lord$ = 4, super_lord$ = 5, overlay_lord$ = 6
  integer, parameter :: component_lord$ = 7, container_slave$ = 8

  character*16 :: control_name(9) = (/ &
            'FREE_ELEMENT   ', 'SUPER_SLAVE    ', 'OVERLAY_SLAVE  ', &
            'GROUP_LORD     ', 'SUPER_LORD     ', 'OVERLAY_LORD   ', &
            'COMPONENT_LORD ', 'CONTAINER_SLAVE', '               ' /)

! plane list

  integer, parameter :: x_plane$ = 1, y_plane$ = 2, n_plane$ = 3

  character*16 plane_name(4) / 'X', 'Y', 'N', ' ' /

! aperture structure

  integer, parameter :: ap_array_maxx = 100

  type aperture_struct
    real dE_E
    type (coord_struct)  closed_orbit
    real x_(ap_array_maxx), y_(0:ap_array_maxx)
    integer plane_(ap_array_maxx)
    integer ix_ring(ap_array_maxx)
    integer i_turn(ap_array_maxx)
  end type

  type track_input_struct
    integer n_turn
    real x_init, y_init
    integer n_xy_pts
    real e_max
    integer n_energy_pts
    real accuracy
  end type

! garbage$ is, for example, for subroutines that want to communicate to
! the calling program that a variable has not been set properly.

  integer, parameter :: garbage$ = -9876

!

  type amode_struct
    real emittance
    real synch_int(4:5)
    real j_damp               ! damping partition number
    real alpha_damp           ! damping per turn
    real chrom                ! Chromaticity
    real tune                 ! "Fractional" tune in radians
  end type

  type modes_struct
    real synch_int(3)
    real sig_e
    real sig_z
    real energy_loss
    type (amode_struct)  a, b, z
  end type

  integer, parameter :: bends$ = 201
  integer, parameter :: wigglers$ = 202
  integer, parameter :: all$ = 203

! other

  integer, parameter :: ok$              = 1
  integer, parameter :: in_stop_band$    = 2
  integer, parameter :: non_symplectic$  = 3
  integer, parameter :: unstable$        = 4

  character*16 :: status_name(5) = (/     'OK            ', &
      'IN_STOP_BAND  ', 'NON_SYMPLECTIC', 'UNSTABLE      ', '              ' /)

! common flags
! status structure

  type ring_master_struct
    character*8 name
    real*8 z, mag_len, phys_len
    character*16 comment
  end type

  type bmad_status_struct
    logical :: ok             = .true.
    logical :: type_out       = .true.
    logical :: sub_type_out   = .true.
    logical :: exit_on_error  = .true.
    integer :: status         = ok$
  end type

  type (bmad_status_struct) bmad_status

!-----------------------------------------------------------------------------
! For butns.nnnnn files

  type detector_struct
    real x_orb, y_orb
    integer amp(4)
  endtype

  type butns_struct
    character*40 lattice
    integer save_set
    character*20 date
    character*72 comment(5)
    integer file_num
    type (detector_struct) det(0:99)
  end type

! For 6x27 matrices

  integer, parameter :: x11$ = 7, x12$ = 8, x13$ = 9, x14$ = 10, x15$ = 11
  integer, parameter :: x16$ = 12, x22$ = 13, x23$ = 14, x24$ = 15
  integer, parameter :: x25$ = 16, x26$ = 17, x33$ = 18, x34$ = 19, x35$ = 20
  integer, parameter :: x36$ = 21, x44$ = 22, x45$ = 23, x46$ = 24
  integer, parameter :: x55$ = 25, x56$ = 26, x66$ = 27

  type mat627_struct
    real m(6,27)
  end type

!---------------------------------------------------------------
! the cesr_struct is used by bmad_to_cesr.
! the element ordering here is similar to what is used by cesrv

  type b_struct
    real beta_mid             ! beta at the midpoint
    real beta_ave             ! beta averaged over the element
    real pos_mid              ! position at midpoint
    real pos_inj              ! position in injection lattice
  end type

  type cesr_element_struct
    character*16 name              ! bmad name
    character*12 db_node_name
    integer ix_db                  ! element index for data base node
    integer ix_ring                ! index to element in ring structure
    type (b_struct)  x, y          ! beta's and positions
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
    integer ix_cesr(n_ele_maxx)  ! pointer to cesr struct for given type
  end type

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

!---------------------------------------------------------------------------
! DB_STRUCT:                       
! This structrue holds info on the correspondence between CESR data base
! elements and a BMAD ring. Use BMAD_TO_DB to initialize this structure.
! For completeness there are, in addition, arrays for the detectors and
! wigglers, etc.

! for an individual element

  type db_element_struct
    character*16 bmad_name    ! bmad name of element
    integer ix_ring           ! index to element array in ring struct
    integer ix_attrib         ! index to element attribute
    integer ix_cesrv          ! index to cesr_struct arrays
    character*12 db_node_name ! node name ("CSR QUAD CUR")
    integer ix_db             ! element index for data base node (5 for Q05W)
    character*16 db_ele_name  ! element name
    real dvar_dcu             ! calibration factor
    real var_theory           ! theory var value
    real var_0                ! extrapolated var value at CU = 0
    integer cu_high_lim       ! high limit
    integer cu_low_lim        ! low limit
    integer cu_now            ! current CU
    integer ix                ! Integer for general use.
  end type              

! for pointing to a db_struct array

  type db_node_struct
    type (db_element_struct), pointer :: ptr(:)
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
    type (db_element_struct) :: quad_cur(0:120) ! combinded csr_quad_cur, etc.
    type (db_element_struct) :: quad_tilt(0:99)
    type (db_element_struct) :: qadd_tilt(1:n_qadd_maxx)
    type (db_element_struct) :: detector(0:99)
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

!----------------------------------------------------------------------------
! Common stuff used by various subroutines

  type bmad_common_struct
    real factor                   ! used by track_runge_kutta
    character*16 func_type        ! used by track_runge_kutta
  end type

  type (bmad_common_struct) bmad_common

end module
