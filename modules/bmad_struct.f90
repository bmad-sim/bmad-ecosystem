!+
! BMAD_STRUCT holds the structure definitions for BMAD routines.
!-

#include "CESR_platform.inc"

module bmad_struct

  use twiss_mod
  use bmad_taylor_mod

  use tpsalie_analysis, only: genfield

! The "regular" elements are in positions: 1 to RING.N_ELE_RING
! regular elements are:
!     1) Superimpose slaves:      SUPER_SLAVE$
!     2) Overlay slaves:          OVERLAY_SLAVE$
!     3) Free:                    FREE$
!
! The "control" elements are in positions: RING.N_ELE_RING+1 to RING.N_ELE_MAX
!
! The control elements are:
!     1) Superimposed controllers:      SUPER_LORD$
!     2) Overlay controllers:           OVERLAY_LORD$
!     3) Group controllers:             GROUP_LORD$
!-

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! IF YOU CHANGE THE RING STRUCTURE YOU MUST INCREASE THE VERSION NUMBER !

  integer, parameter :: bmad_inc_version$ = 60

! THIS IS USED BY BMAD_PARSER TO MAKE SURE DIGESTED FILES ARE OK.
!
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

! parameter def

  integer, parameter :: n_ele_maxx = 15000
  integer, parameter :: n_control_maxx = 5000
  integer, parameter :: n_attrib_maxx = 34

! Structure for a particle's coordinates.
! Coordinates are with respect to the reference trajectory.
! See the MAD-8 Physical Methods Manual.
! Usually (but there are some exceptions):
!   vec(1) -- x position: horizontal .
!   vec(2) -- p_x normalized momentum: P_x/P_0
!   vec(3) -- y position: vertical.
!   vec(4) -- p_y normalized momentum: P_Y/P_0
!   vec(5) -- z position: = -c*delta_t
!   vec(6) -- p_z normalized energy deviation: Delta_E/E_0
! Here P_0 is the reference momentum and E_0 is the reference energy.
! P_X and P_Y are the momentum along the x and y-axes.

  type coord_struct
    real(rdef) vec(6)
  end type

! Wiggler structures

  integer,parameter :: hyper_y$ = 1, hyper_xy$ = 2, hyper_x$ = 3
  character*8, parameter :: wig_term_type_name(0:3) = (/ &
                  'Garbage ', 'Hyper_Y ', 'Hyper_XY', 'Hyper_X ' /)

  type wig_term_struct
    real(rdef) coef
    real(rdef) kx, ky, kz
    real(rdef) phi_z
    integer type      ! hyper_y$, hyper_xy$, or hyper_x$
  end type

! Ele_struct
! Remember: If this struct is changed you have to modify:
!     read_digested_bmad_file
!     write_digested_bmad_file

  type ele_struct
    character*16 name              ! name of element
    character*16 type              ! type name (see MAD documentation)
    character*16 alias             ! Another name
    character*16 attribute_name    ! Used by overlays
    type (twiss_struct)  x,y,z         ! Twiss parameters at end of element
    real(rdef) value(n_attrib_maxx)    ! attribute values
    real(rdef) gen0(6)                 ! constant part of the genfield map
    real(rdef) vec0(6)                 ! 0th order transport vector
    real(rdef) mat6(6,6)               ! 1st order transport matrix 
    real(rdef) c_mat(2,2)              ! 2x2 C coupling matrix
    real(rdef) gamma_c                 ! gamma associated with C matrix
    real(rdef) s                       ! longitudinal position at the end
    real(rdef) x_position, y_position  ! Floor position of element
    real(rdef) z_position              ! Elevation of element
    real(rdef) theta_position          ! Floor orientation angle of element
    real(rdef) phi_position            ! Angle of attack
    real(rdef), pointer             :: a(:), b(:)     ! multipoles
    real(rdef), pointer             :: const(:)       ! Working constants.
    character*200, pointer          :: descrip        ! For general use
    type (wig_term_struct), pointer :: wig_term(:)    ! Wiggler Coefs
    type (genfield), pointer        :: gen_field      ! For symp_map$
    type (taylor_struct)            :: taylor(6)      ! Taylor terms
    integer key                    ! key value
    integer sub_key                ! For wigglers: map_type$, periodic_type$
    integer pointer_init           ! Pointers initialized when 
                                   !   pointer_init = has_been_inited$
    integer control_type           ! SUPER_SLAVE$, OVERLAY_LORD$, etc.
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
    integer mat6_calc_method       ! bmad_standard$, taylor$, etc.
    integer tracking_method        ! bmad_standard$, taylor$, etc.
    integer num_steps              ! number of slices for DA_maps
    integer integration_order      ! For Etiennes' PTC: 2, 4, or 6.
    integer ptc_kind               ! For setting the ptc kind type.
    integer taylor_order           ! Order of the taylor series.
    logical symplectify            ! Symplectify mat6 matrices.
    logical mode_flip              ! Have the normal modes traded places?
    logical multipoles_on          ! For turning multipoles on/off
    logical exact_rad_int_calc     ! Exact radiation integral calculation?
    logical B_field_master         ! Calculate K from the B_field?
    logical is_on                  ! For turning element on/off.
  end type

! struct for element to element control

  type control_struct
    real(rdef) coef                ! control coefficient
    integer ix_lord                ! index to lord element
    integer ix_slave               ! index to slave element
    integer ix_attrib              ! index of attribute controlled
  end type

! parameter and mode structures

  type param_struct
    real(rdef) energy             ! USE BEAM_ENERGY INSTEAD ! energy in GeV.
    real(rdef) beam_energy        ! beam energy in eV
    real(rdef) n_part             ! Number of particles in a bunch
    real(rdef) total_length       ! total_length of ring
    real(rdef) growth_rate        ! growth rate/turn if not stable
    integer particle              ! +1 = positrons, -1 = electrons
    integer symmetry              ! symmetry of the ring (e/w symm, etc.)
    integer ix_lost               ! If lost at what element?
    integer lattice_type          ! linac_lattice$, etc...
    integer ixx                   ! Integer for general use
    logical stable                ! is closed ring stable?
    logical aperture_limit_on     ! use apertures in tracking?
    logical lost                  ! for use in tracking
    logical damping_on            ! Radiation damping calculation toggle.
    logical fluctuations_on       ! Radiation fluctuations calculation toggle.
  end type

  type mode_info_struct
    real(rdef) tune      ! "fractional" tune in radians: 0 < tune < 2pi
    real(rdef) emit      ! Emittance
    real(rdef) chrom     ! Chromaticity
  end type

  type dummy_parameter_struct
    integer dummy(250)
  end type

! RING_STRUCT
! Remember: If this struct is changed you have to modify:
!     read_digested_bmad_file
!     write_digested_bmad_file
!     transfer_ring_parameters

  type ring_struct
    character*16 name             ! Name of ring given by USE statement
    character*40 lattice          ! Lattice
    character*80 input_file_name  ! name of the lattice input file
    character*80 title            ! general title
    type (mode_info_struct)  x, y, z  ! tunes, etc.
    type (param_struct)      param ! parameters
    integer version               ! Version number
    integer n_ele_ring            ! number of regular ring elements
    integer n_ele_symm            ! symmetry point for rings w/ symmetry
    integer n_ele_use             ! number of elements used
    integer n_ele_max             ! Index of last element used
    integer n_control_array       ! last index used in CONTROL_ array
    integer n_ic_array            ! last index used in IC_ array
    integer input_taylor_order    ! As set in the input file
    type (ele_struct)  ele_(0:n_ele_maxx)    ! Array of ring elements
    type (ele_struct)  ele_init              ! For use by any program
    type (control_struct)  control_(n_control_maxx)  ! control list
    integer ic_(n_control_maxx)                      ! index to %control_(:)
  end type


!

  character*3, parameter :: coord_name(6) = &
                              (/ "X  ", "P_x", "Y  ", "P_y", "Z  ", "P_z" /)

! KEY value definitions
! Note: overlay$ == overlay_lord$ 

  integer, parameter :: drift$ = 1, sbend$ = 2, quadrupole$ = 3, group$ = 4
  integer, parameter :: sextupole$ = 5, overlay$ = 6, custom$ = 7, taylor$ = 8
  integer, parameter :: rfcavity$ = 9
  integer, parameter :: elseparator$ = 10, beambeam$ = 11, wiggler$ = 12
  integer, parameter :: sol_quad$ = 13, marker$ = 14, kicker$ = 15
  integer, parameter :: hybrid$ = 16, octupole$ = 17, rbend$ = 18
  integer, parameter :: multipole$ = 19, accel_sol$ = 20
  integer, parameter :: def_beam$ = 21, ab_multipole$ = 22, solenoid$ = 23
  integer, parameter :: patch$ = 24, lcavity$ = 25, def_parameter$ = 26
  integer, parameter :: null_ele$ = 27, init_ele$ = 28, hom$ = 29
  integer, parameter :: matrix$ = 30

  integer, parameter :: n_key = 30

  character*16 :: key_name(n_key+1) = (/ &
    'DRIFT        ', 'SBEND        ', 'QUADRUPOLE   ', 'GROUP        ', &
    'SEXTUPOLE    ', 'OVERLAY      ', 'CUSTOM       ', 'TAYLOR       ', &
    'RFCAVITY     ', 'ELSEPARATOR  ', 'BEAMBEAM     ', 'WIGGLER      ', &
    'SOL_QUAD     ', 'MARKER       ', 'KICKER       ', 'HYBRID       ', &
    'OCTUPOLE     ', 'RBEND        ', 'MULTIPOLE    ', 'ACCEL_SOL    ', &
    'DEF BEAM     ', 'AB_MULTIPOLE ', 'SOLENOID     ', 'PATCH        ', &
    'LCAVITY      ', 'DEF PARAMETER', 'NULL_ELEMENT ', 'INIT_ELEMENT ', &
    'HOM          ', 'MATRIX       ', '             ' /)

! Attribute name logical definitions
! Note: The following attributes must have unique number assignments:
!     L$, TILT$, X_PITCH$ and higher

  integer, parameter :: lattice_type$ = 1, symmetry$ = 2, taylor_order$ = 3
  integer, parameter :: beam_energy$ = 4 

  integer, parameter :: x_beg_limit$=2, y_beg_limit$=3, b_x2$=4, &
          b_y2$=5, l_st2$=9, b_z$=10, l_st1$=11, s_st2$=12, s_st1$=13, &
          b_x1$=14, b_y1$=15    

  integer, parameter :: val1$=3, val2$=4, val3$=5, val4$=6, val5$=7, &
          val6$=8, val7$=9, val8$=10, val9$=11, val10$=12, val11$=13, &
          val12$=14

  integer, parameter :: beta_x0$ = 2, alpha_x0$ = 3, beta_y0$ = 4, &
          alpha_y0$ = 5, beta_x1$ = 6, alpha_x1$ = 7, beta_y1$ = 8, &
          alpha_y1$ = 9, dphi_x$ = 10, dphi_y$ = 11

!  integer, parameter :: x_position$ = 2, y_position$ = 3, z_position$ = 4, &
!          theta_position$ = 5, phi_position$ = 6, psi_position$ = 7, &
!          beta_x$ = 8, beta_y$ = 9, alpha_x$ = 10, alpha_y$ = 11, &
!          eta_x$ = 12, eta_y$ = 13, etap_x$ = 14, etap_y$ = 15, &
!          phase_x$ = 16, phase_y$ = 17, &
!          c11$ = 18, c12$ = 19, c21$ = 20, c22$ = 21

  integer, parameter :: l$=1
  integer, parameter :: tilt$=2, command$=2
  integer, parameter :: old_command$=3, angle$=3
  integer, parameter :: k1$=4, sig_x$=4, harmon$=4, h_displace$=4
  integer, parameter :: k2$=5, sig_y$=5, b_max$=5, v_displace$=5, g$=5
  integer, parameter :: k3$=6, sig_z$=6, rf_wavelength$=6, delta_g$=6
  integer, parameter :: ks$=7, volt$=7, e1$=7, n_pole$=7, bbi_const$=7
  integer, parameter :: e2$=8, charge$=8, gap$=8
  integer, parameter :: n_slice$=9, l_chord$=9, l_pole$=9, rf_frequency$=9
  integer, parameter :: fint$=10, polarity$ = 10, gradient$=10
  integer, parameter :: fintx$=11, z_patch$ = 11, phi0$=11
  integer, parameter :: rho$ = 12
  integer, parameter :: hgap$=13, energy_start$=13
  integer, parameter :: coef$=14, current$=14, hgapx$=14
  integer, parameter :: roll$=15
  integer, parameter :: l_original$ = 16
  integer, parameter :: l_start$=17
  integer, parameter :: l_end$=18
  integer, parameter :: x_pitch$=19
  integer, parameter :: y_pitch$=20
  integer, parameter :: hkick$=21
  integer, parameter :: vkick$=22
  integer, parameter :: x_offset$=23
  integer, parameter :: y_offset$=24
  integer, parameter :: s_offset$=25, z_offset$ = 25
  integer, parameter :: dE_offset$ = 26
  integer, parameter :: x_limit$=27
  integer, parameter :: y_limit$=28
  integer, parameter :: aperture$=29
  integer, parameter :: radius$=30
  integer, parameter :: energy$=31  ! formally new_energy$
  integer, parameter :: rel_tol$ = 32
  integer, parameter :: abs_tol$ = 33
  integer, parameter :: B_field$ = 34, B_gradient$ = 34

  integer, parameter :: type$ = 35   ! this is 1 greater than n_attrib_maxx
  integer, parameter :: alias$ = 36 
  integer, parameter :: start_edge$ = 37     ! special for groups
  integer, parameter :: end_edge$ = 38       ! special for groups
  integer, parameter :: accordion_edge$ = 39 ! special for groups
  integer, parameter :: symmetric_edge$ = 40
  integer, parameter :: mat6_calc_method$ = 41
  integer, parameter :: tracking_method$  = 42
  integer, parameter :: num_steps$ = 43
  integer, parameter :: integration_order$ = 44
  integer, parameter :: term$ = 45
  integer, parameter :: ptc_kind$ = 46
  integer, parameter :: symplectify$ = 47
  integer, parameter :: descrip$ = 48
  integer, parameter :: is_on$ = 49

! Warning: No other attribute parameters can have indexes larger than A0$.
! That is: multipole arrays An, Bn, KnL, and Tn must have the largest indexes

  integer, parameter :: a0$  = 50, k0l$  = 50
  integer, parameter :: a1$  = 51, k1l$  = 51
  integer, parameter :: a2$  = 52, k2l$  = 52
  integer, parameter :: a3$  = 53, k3l$  = 53
  integer, parameter :: a4$  = 54, k4l$  = 54
  integer, parameter :: a5$  = 55, k5l$  = 55
  integer, parameter :: a6$  = 56, k6l$  = 56
  integer, parameter :: a7$  = 57, k7l$  = 57
  integer, parameter :: a8$  = 58, k8l$  = 58
  integer, parameter :: a9$  = 59, k9l$  = 59
  integer, parameter :: a10$ = 60, k10l$ = 60
  integer, parameter :: a11$ = 61, k11l$ = 61
  integer, parameter :: a12$ = 62, k12l$ = 62
  integer, parameter :: a13$ = 63, k13l$ = 63
  integer, parameter :: a14$ = 64, k14l$ = 64
  integer, parameter :: a15$ = 65, k15l$ = 65
  integer, parameter :: a16$ = 66, k16l$ = 66
  integer, parameter :: a17$ = 67, k17l$ = 67
  integer, parameter :: a18$ = 68, k18l$ = 68
  integer, parameter :: a19$ = 69, k19l$ = 69
  integer, parameter :: a20$ = 70, k20l$ = 70

  integer, parameter :: b0$  = 80, t0$  = 80
  integer, parameter :: b1$  = 81, t1$  = 81
  integer, parameter :: b2$  = 82, t2$  = 82
  integer, parameter :: b3$  = 83, t3$  = 83
  integer, parameter :: b4$  = 84, t4$  = 84
  integer, parameter :: b5$  = 85, t5$  = 85
  integer, parameter :: b6$  = 86, t6$  = 86
  integer, parameter :: b7$  = 87, t7$  = 87
  integer, parameter :: b8$  = 88, t8$  = 88
  integer, parameter :: b9$  = 89, t9$  = 89
  integer, parameter :: b10$ = 90, t10$ = 90
  integer, parameter :: b11$ = 91, t11$ = 91
  integer, parameter :: b12$ = 92, t12$ = 92
  integer, parameter :: b13$ = 93, t13$ = 93
  integer, parameter :: b14$ = 94, t14$ = 94
  integer, parameter :: b15$ = 95, t15$ = 95
  integer, parameter :: b16$ = 96, t16$ = 96
  integer, parameter :: b17$ = 97, t17$ = 97
  integer, parameter :: b18$ = 98, t18$ = 98
  integer, parameter :: b19$ = 99, t19$ = 99
  integer, parameter :: b20$ =100, t20$ =100 ! this is n_attrib_special_maxx 

  integer, parameter :: n_attrib_special_maxx = 100

  integer, parameter :: particle$ = 1
  integer, parameter :: n_part$    = 3

  character*16, parameter :: null_name = 'NULL' 
  character*16, parameter :: blank_name = ' '

! electron/positron

  integer, parameter :: proton$     = +2
  integer, parameter :: positron$   = +1
  integer, parameter :: electron$   = -1
  integer, parameter :: antiproton$ = -2

  character*16 :: particle_name(-2:2) = (/ 'ANTIPROTON', &
                  'ELECTRON  ', '???       ', 'POSITRON  ', 'PROTON    ' /)

  integer, parameter :: charge_of(-2:2) = (/ -1, -1, 0, 1, 1 /)
  real(rdef), parameter :: mass_of(-2:2) = (/ 0.938279e9, 0.511003e6, 0.0, &
                                              0.511003e9, 0.938279e9 /)

! SYMMETRY etc., logical names

  integer, parameter :: no_symmetry$ = 0
  integer, parameter :: mobius_symmetry$ = 2
  integer, parameter :: ew_antisymmetry$ = -1  ! coupling elements are antisymmetric
  integer, parameter :: ew_symmetry$ = 3
  integer, parameter :: linac_lattice$ = 10
  integer, parameter :: linear_lattice$ = 11
  integer, parameter :: circular_lattice$ = 12

! logicals for MAKE_HYBIRD_RING

  logical, parameter :: remove_markers$ = .true., no_remove_markers$ = .false.

! control element logicals
! Note: Must have overlay_lord$ == overlay$ !!

  integer, parameter :: free$ = 1, super_slave$ = 2, overlay_slave$ = 3
  integer, parameter :: group_lord$ = 4, super_lord$ = 5, overlay_lord$ = 6

  character*16 :: control_name(7) = (/ &
            'FREE_ELEMENT   ', 'SUPER_SLAVE    ', 'OVERLAY_SLAVE  ', &
            'GROUP_LORD     ', 'SUPER_LORD     ', 'OVERLAY_LORD   ', &
            '               ' /)

! plane list, etc

  integer, parameter :: x_plane$ = 1, y_plane$ = 2
  integer, parameter :: z_plane$ = 3, n_plane$ = 4

  character(16) :: plane_name(5) = (/ 'X', 'Y', 'Z', 'N', ' ' /)

  logical, parameter :: set$ = .true., unset$ = .false.
  integer, parameter :: has_been_inited$ = 123456  ! something random

! garbage$ is, for example, for subroutines that want to communicate to
! the calling program that a variable has not been set properly.

  integer, parameter :: garbage$ = -9876

! Note: custom$ = 7, and taylor$ = 8 are taken from the element key list.

  integer, parameter :: bmad_standard$ = 1, symp_lie_ptc$ = 2
  integer, parameter :: runge_kutta$ = 3 
  integer, parameter :: linear$ = 4, tracking$ = 5, symp_map$ = 6
  integer, parameter :: wiedemann$ = 9, symp_lie_bmad$ = 10, none$ = 11
  integer, parameter :: boris$ = 12, adaptive_boris$ = 13, order_2$ = 14

  character(16), parameter :: calc_method_name(0:14) = (/ &
      "GARBAGE!      ", "BMAD_Standard ", "Symp_Lie_PTC  ", "Runge_Kutta   ", &
      "Linear        ", "Tracking      ", "Symp_Map      ", "Custom        ", &
      "Taylor        ", "Wiedemann     ", "Symp_Lie_BMAD ", "None          ", &
      "Boris         ", "Adaptive_Boris", "Order_2       " /)

  integer, parameter :: map_type$ = 1, periodic_type$ = 2
  character(16), parameter :: sub_key_name(0:2) = (/ &
         "GARBAGE!     ", "Map          ", "Periodic     " /)

!

  type amode_struct
    real(rdef) emittance
    real(rdef) synch_int(4:5)
    real(rdef) j_damp               ! damping partition number
    real(rdef) alpha_damp           ! damping per turn
    real(rdef) chrom                ! Chromaticity
    real(rdef) tune                 ! "Fractional" tune in radians
  end type

  type modes_struct
    real(rdef) synch_int(3)
    real(rdef) sig_e
    real(rdef) sig_z
    real(rdef) e_loss
    type (amode_struct)  a, b, z
  end type

  integer, parameter :: bends$ = 201
  integer, parameter :: wigglers$ = 202
  integer, parameter :: all$ = 203

! common flags
! status structure

  type bmad_status_struct
    logical :: ok             = .true.
    logical :: type_out       = .true.
    logical :: sub_type_out   = .true.
    logical :: exit_on_error  = .true.
    integer :: status         = ok$
  end type

  type (bmad_status_struct), save :: bmad_status

!-----------------------------------------------------------------------------
! For butns.nnnnn files

  type detector_struct
    real(rdef) x_orb, y_orb
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
  end type

! For 6x27 matrices

  integer, parameter :: x11$ = 7, x12$ = 8, x13$ = 9, x14$ = 10, x15$ = 11
  integer, parameter :: x16$ = 12, x22$ = 13, x23$ = 14, x24$ = 15
  integer, parameter :: x25$ = 16, x26$ = 17, x33$ = 18, x34$ = 19, x35$ = 20
  integer, parameter :: x36$ = 21, x44$ = 22, x45$ = 23, x46$ = 24
  integer, parameter :: x55$ = 25, x56$ = 26, x66$ = 27

  type mat627_struct
    real(rdef) m(6,27)
  end type

!---------------------------------------------------------------
! the cesr_struct is used by bmad_to_cesr.
! the element ordering here is similar to what is used by cesrv

  type b_struct
    real(rdef) beta_mid             ! beta at the midpoint
    real(rdef) beta_ave             ! beta averaged over the element
    real(rdef) pos_mid              ! position at midpoint
    real(rdef) pos_inj              ! position in injection lattice
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
! Units

  integer, parameter :: radians$ = 1, degrees$ = 2, cycles$ = 3, kHz$ = 4
  character*8 ::frequency_units_name(4) = (/ &
            'Radians ', 'Degrees ', 'Cycles  ', 'kHz     ' /)

!---------------------------------------------------------------------------
! DB_STRUCT:                       
! This structrue holds info on the correspondence between CESR data base
! elements and a BMAD ring. Use BMAD_TO_DB to initialize this structure.
! For completeness there are, in addition, arrays for the detectors and
! wigglers, etc.

! for an individual element

  type db_element_struct
    character*16 bmad_name    ! bmad name of element
    real(rdef) dvar_dcu             ! calibration factor
    real(rdef) var_theory           ! theory var value
    real(rdef) var_0                ! extrapolated var value at CU = 0
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

! electric and magnetic fields

  type em_field_struct
    real(rdef) E(3)         ! electric field
    real(rdef) B(3)         ! magnetic field
  end type

!------------------------------------------------------------------------------
! common stuff

! %taylor_order_ptc is what ptc has been set to.
! %taylor_order is what the user wants.
! %max_aperture_limit is used when no limit is specified or when 
!   ring%param%aperture_limit_on = False.

  type bmad_com_struct
    type (coord_struct) :: d_orb  ! for the transfer_mat_from_tracking routine
    real(rdef) :: beam_energy = 0
    real(rdef) :: max_aperture_limit = 1e3    
    integer :: taylor_order = 3              ! 3rd order is default
    integer :: taylor_order_ptc = 0          ! 0 -> not yet set 
    logical :: taylor_order_set = .false.    ! Used by set_taylor_order
    integer :: real_8_map_init
    integer :: default_integ_order = 2
    integer :: default_num_steps = 1
    logical :: init_needed = .true.
  end type
  
  type (bmad_com_struct), save :: bmad_com

! multi_turn_func_common is for multi_turn_tracking_to_mat.

  type (coord_struct), pointer :: multi_turn_func_common(:)  

end module
