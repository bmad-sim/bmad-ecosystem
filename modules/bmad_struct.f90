!+
! Bmad_struct holds the structure definitions for Bmad routines.
!-

#include "CESR_platform.inc"

module bmad_struct

  use twiss_mod
  use bmad_taylor_mod

  use tpsalie_analysis, only: genfield

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! IF YOU CHANGE THE RING STRUCTURE YOU MUST INCREASE THE VERSION NUMBER !

  integer, parameter :: bmad_inc_version$ = 71

! THIS IS USED BY BMAD_PARSER TO MAKE SURE DIGESTED FILES ARE OK.
!
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

! parameter def

  integer, parameter :: n_attrib_maxx = 41

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
    real(rp) vec(6)
  end type

! Wiggler structures

  integer,parameter :: hyper_y$ = 1, hyper_xy$ = 2, hyper_x$ = 3
  character*8, parameter :: wig_term_type_name(0:3) = (/ &
                  'Garbage ', 'Hyper_Y ', 'Hyper_XY', 'Hyper_X ' /)

  type wig_term_struct
    real(rp) coef
    real(rp) kx, ky, kz
    real(rp) phi_z
    integer type      ! hyper_y$, hyper_xy$, or hyper_x$
  end type

! Wakefield structs

  type sr_wake_struct     ! Short-Range Wake struct
    real(rp) z            ! Distance behind the leading particle
    real(rp) long         ! Longitudinal wake in V/C/m
    real(rp) trans        ! Transverse wake in V/C/m^2
  end type

  type lr_wake_struct   ! Long-Range Wake struct (transverse only)
    real(rp) freq       ! frequency in Hz
    real(rp) kick       ! Strength in V/C/m^2
    integer  i_cell     ! Index of the cell the mode is in.
    real(rp) Q          ! Quality factor
  end type

  type wake_struct
    character(200), pointer :: sr_file => null()
    character(200), pointer :: lr_file => null()
    type (sr_wake_struct), pointer :: sr(:) => null()
    type (lr_wake_struct), pointer :: lr(:) => null()
  end type

! Local reference frame position with respect to the global (floor) coordinates

  type floor_position_struct
    real(rp) x, y, z            ! offset from origin
    real(rp) theta, phi, psi    ! angular orientation
  end type

! Ele_struct
! REMEMBER: If this struct is changed you have to modify:
!     read_digested_bmad_file
!     write_digested_bmad_file
!     init_ele (in bmad_utils_mod)

  type ele_struct
    character(16) name              ! name of element
    character(16) type              ! type name
    character(16) alias             ! Another name
    character(16) attribute_name    ! Used by overlays
    type (twiss_struct)  x,y,z       ! Twiss parameters at end of element
    type (floor_position_struct) position
    real(rp) value(n_attrib_maxx)    ! attribute values
    real(rp) gen0(6)                 ! constant part of the genfield map
    real(rp) vec0(6)                 ! 0th order transport vector
    real(rp) mat6(6,6)               ! 1st order transport matrix 
    real(rp) c_mat(2,2)              ! 2x2 C coupling matrix
    real(rp) gamma_c                 ! gamma associated with C matrix
    real(rp) s                       ! longitudinal position at the end
    real(rp), pointer :: r(:,:) => null()  ! For general use. Not used by Bmad.
    real(rp), pointer :: a(:) => null()              ! multipole
    real(rp), pointer :: b(:) => null()              ! multipoles
    real(rp), pointer :: const(:) => null()          ! Working constants.
    character(200), pointer :: descrip => null()     ! For general use
    type (genfield), pointer :: gen_field => null()  ! For symp_map$
    type (taylor_struct) :: taylor(6)                ! Taylor terms
    type (wig_term_struct), pointer :: wig_term(:) => null()   ! Wiggler Coefs
    type (wake_struct) wake        ! Wakefields
    integer key                    ! key value
    integer sub_key                ! For wigglers: map_type$, periodic_type$
    integer control_type           ! SUPER_SLAVE$, OVERLAY_LORD$, etc.
    integer ix_value               ! Pointer for attribute to control
    integer n_slave                ! Number of slaves
    integer ix1_slave              ! Start index for slave elements
    integer ix2_slave              ! Stop  index for slave elements
    integer n_lord                 ! Number of lords
    integer ic1_lord               ! Start index for lord elements
    integer ic2_lord               ! Stop  index for lord elements
    integer ix_pointer             ! For general use. Not used by Bmad.
    integer ixx                    ! Index for Bmad internal use
    integer iyy                    ! Index for Bmad internal use
    integer mat6_calc_method       ! bmad_standard$, taylor$, etc.
    integer tracking_method        ! bmad_standard$, taylor$, etc.
    integer field_calc             ! Used with Boris or Runge-Kutta integrators.
    integer num_steps              ! number of slices for DA_maps
    integer integration_ord        ! For Etiennes' PTC: 2, 4, or 6.
    integer ptc_kind               ! For setting the ptc kind type.
    integer taylor_order           ! Order of the taylor series.
    logical symplectify            ! Symplectify mat6 matrices.
    logical mode_flip              ! Have the normal modes traded places?
    logical multipoles_on          ! For turning multipoles on/off
    logical exact_rad_int_calc     ! Exact radiation integral calculation?
    logical field_master           ! Calculate strength from the field value?
    logical is_on                  ! For turning element on/off.
    logical internal_logic         ! For Bmad internal use only.
    logical logic                  ! For general use. Not used by Bmad.
  end type

! struct for element to element control

  type control_struct
    real(rp) coef                  ! control coefficient
    integer ix_lord                ! index to lord element
    integer ix_slave               ! index to slave element
    integer ix_attrib              ! index of attribute controlled
  end type

! parameter and mode structures

  type param_struct
    real(rp) garbage            ! Saved for future use.
    real(rp) n_part             ! Particles/bunch (for BeamBeam elements).
    real(rp) charge             ! Charge of a bunch (used by LCavities).
    real(rp) total_length       ! total_length of ring
    real(rp) growth_rate        ! growth rate/turn if not stable
    real(rp) t1_with_RF(6,6)    ! Full 1-turn matrix with RF on.
    real(rp) t1_no_RF(6,6)      ! Full 1-turn matrix with RF off.
    integer particle            ! +1 = positrons, -1 = electrons
    integer iy                  ! Not currently used.
    integer ix_lost             ! If lost at what element?
    integer lattice_type        ! linac_lattice$, etc...
    integer ixx                 ! Integer for general use
    logical stable              ! is closed ring stable?
    logical aperture_limit_on   ! use apertures in tracking?
    logical lost                ! for use in tracking
  end type

  type mode_info_struct
    real(rp) tune      ! "fractional" tune in radians: 0 < tune < 2pi
    real(rp) emit      ! Emittance
    real(rp) chrom     ! Chromaticity
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
    character(16) name               ! Name of ring given by USE statement
    character(40) lattice            ! Lattice
    character(200) input_file_name   ! Name of the lattice input file
    character(80) title              ! General title
    type (mode_info_struct) x, y, z  ! Tunes, etc.
    type (param_struct) param        ! Parameters
    integer version                  ! Version number
    integer n_ele_use                ! Number of regular ring elements
    integer n_ele_ring               ! OBSOLETE: Identical to n_ele_use.
    integer n_ele_max                ! Index of last element used
    integer n_ele_maxx               ! Index of last element allocated
    integer n_control_max            ! Last index used in CONTROL_ array
    integer n_ic_max                 ! Last index used in IC_ array
    integer input_taylor_order       ! As set in the input file
    type (ele_struct)  ele_init      ! For use by any program
    type (ele_struct), pointer ::  ele_(:) => null()        ! Array of elements
    type (control_struct), pointer :: control_(:) => null() ! control list
    integer, pointer :: ic_(:) => null()                ! index to %control_(:)
    real(rp), pointer :: beam_energy ! points to ring%ele_(0)%value(beam_energy$)
  end type

!

  character(3), parameter :: coord_name(6) = &
                              (/ "X  ", "P_x", "Y  ", "P_y", "Z  ", "P_z" /)

! KEY value definitions
! Note: sbend$ and rbend$ also used for sub_key

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
  integer, parameter :: matrix$ = 30, monitor$ = 31, instrument$ = 32
  integer, parameter :: hkicker$ = 33, vkicker$ = 34, rcollimator$ = 35
  integer, parameter :: ecollimator$ = 36, i_beam$ = 37

  integer, parameter :: n_key = 37

  character(16) :: key_name(n_key+1) = (/ &
    'DRIFT        ', 'SBEND        ', 'QUADRUPOLE   ', 'GROUP        ', &
    'SEXTUPOLE    ', 'OVERLAY      ', 'CUSTOM       ', 'TAYLOR       ', &
    'RFCAVITY     ', 'ELSEPARATOR  ', 'BEAMBEAM     ', 'WIGGLER      ', &
    'SOL_QUAD     ', 'MARKER       ', 'KICKER       ', 'HYBRID       ', &
    'OCTUPOLE     ', 'RBEND        ', 'MULTIPOLE    ', 'ACCEL_SOL    ', &
    'DEF BEAM     ', 'AB_MULTIPOLE ', 'SOLENOID     ', 'PATCH        ', &
    'LCAVITY      ', 'DEF PARAMETER', 'NULL_ELEMENT ', 'INIT_ELEMENT ', &
    'HOM          ', 'MATRIX       ', 'MONITOR      ', 'INSTRUMENT   ', &
    'HKICKER      ', 'VKICKER      ', 'RCOLLIMATOR  ', 'ECOLLIMATOR  ', &
    'I_BEAM       ', '             ' /)

! Attribute name logical definitions
! Note: The following attributes must have unique number assignments:
!     L$, TILT$, X_PITCH$ and higher

  integer, parameter :: lattice_type$ = 1, symmetry$ = 2, taylor_order$ = 3
  integer, parameter :: energy_gev$ = 4 

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
  integer, parameter :: old_command$=3, angle$=3, kick$=3
  integer, parameter :: k1$=4, sig_x$=4, harmon$=4, h_displace$=4, e_loss$=4
  integer, parameter :: k2$=5, sig_y$=5, b_max$=5, v_displace$=5, g$=5
  integer, parameter :: k3$=6, sig_z$=6, rf_wavelength$=6, delta_g$=6
  integer, parameter :: ks$=7, voltage$=7, e1$=7, n_pole$=7, bbi_const$=7
  integer, parameter :: e2$=8, charge$=8, gap$=8
  integer, parameter :: n_slice$=9, l_chord$=9, l_pole$=9, rf_frequency$=9
  integer, parameter :: fint$=10, polarity$=10, gradient$=10
  integer, parameter :: fintx$=11, z_patch$=11, phi0$=11
  integer, parameter :: rho$=12, s_center$=12
  integer, parameter :: hgap$=13, energy_start$=13
  integer, parameter :: coef$=14, current$=14, hgapx$=14, delta_e$=14
  integer, parameter :: roll$=15
  integer, parameter :: l_original$=16
  integer, parameter :: l_start$=17, h1$=17
  integer, parameter :: l_end$=18, h2$=18
  integer, parameter :: x_pitch$=19
  integer, parameter :: y_pitch$=20
  integer, parameter :: hkick$=21
  integer, parameter :: vkick$=22
  integer, parameter :: x_offset$=23
  integer, parameter :: y_offset$=24
  integer, parameter :: s_offset$=25, z_offset$=25
  integer, parameter :: dE_offset$=26, check_sum$=26
  integer, parameter :: x_limit$=27
  integer, parameter :: y_limit$=28
  integer, parameter :: aperture$=29
  integer, parameter :: radius$=30
  integer, parameter :: beam_energy$=31
  integer, parameter :: rel_tol$=32
  integer, parameter :: abs_tol$=33
  integer, parameter :: B_field$=34, B_gradient$=34, &
                             E_field$=34, E_gradient$=34
  integer, parameter :: tilt_tot$=35
  integer, parameter :: x_pitch_tot$=36
  integer, parameter :: y_pitch_tot$=37
  integer, parameter :: x_offset_tot$=38
  integer, parameter :: y_offset_tot$=39
  integer, parameter :: s_offset_tot$=40
  integer, parameter :: p0c$ = 41

  integer, parameter :: alias$ = 42  ! this is 1 greater than n_attrib_maxx
  integer, parameter :: start_edge$ = 43     
  integer, parameter :: end_edge$ = 44     
  integer, parameter :: accordion_edge$ = 45, sr_wake_file$ = 45 
  integer, parameter :: symmetric_edge$ = 46, lr_wake_file$ = 46
  integer, parameter :: mat6_calc_method$ = 47
  integer, parameter :: tracking_method$  = 48
  integer, parameter :: num_steps$ = 49
  integer, parameter :: integration_ord$ = 50
  integer, parameter :: term$ = 51
  integer, parameter :: ptc_kind$ = 52
  integer, parameter :: symplectify$ = 53
  integer, parameter :: descrip$ = 54
  integer, parameter :: is_on$ = 55
  integer, parameter :: field_calc$ = 56
  integer, parameter :: type$ = 57

! Warning: No other attribute parameters can have indexes larger than A0$.
! That is: multipole arrays An, Bn, KnL, and Tn must have the largest indexes

  integer, parameter :: a0$  =  60, k0l$  =  60
  integer, parameter :: a1$  =  61, k1l$  =  61
  integer, parameter :: a2$  =  62, k2l$  =  62
  integer, parameter :: a3$  =  63, k3l$  =  63
  integer, parameter :: a4$  =  64, k4l$  =  64
  integer, parameter :: a5$  =  65, k5l$  =  65
  integer, parameter :: a6$  =  66, k6l$  =  66
  integer, parameter :: a7$  =  67, k7l$  =  67
  integer, parameter :: a8$  =  68, k8l$  =  68
  integer, parameter :: a9$  =  69, k9l$  =  69
  integer, parameter :: a10$ =  70, k10l$ =  70
  integer, parameter :: a11$ =  71, k11l$ =  71
  integer, parameter :: a12$ =  72, k12l$ =  72
  integer, parameter :: a13$ =  73, k13l$ =  73
  integer, parameter :: a14$ =  74, k14l$ =  74
  integer, parameter :: a15$ =  75, k15l$ =  75
  integer, parameter :: a16$ =  76, k16l$ =  76
  integer, parameter :: a17$ =  77, k17l$ =  77
  integer, parameter :: a18$ =  78, k18l$ =  78
  integer, parameter :: a19$ =  79, k19l$ =  79
  integer, parameter :: a20$ =  80, k20l$ =  80

  integer, parameter :: b0$  =  90, t0$  =  90
  integer, parameter :: b1$  =  91, t1$  =  91
  integer, parameter :: b2$  =  92, t2$  =  92
  integer, parameter :: b3$  =  93, t3$  =  93
  integer, parameter :: b4$  =  94, t4$  =  94
  integer, parameter :: b5$  =  95, t5$  =  95
  integer, parameter :: b6$  =  96, t6$  =  96
  integer, parameter :: b7$  =  97, t7$  =  97
  integer, parameter :: b8$  =  98, t8$  =  98
  integer, parameter :: b9$  =  99, t9$  =  99
  integer, parameter :: b10$ = 100, t10$ = 100
  integer, parameter :: b11$ = 101, t11$ = 101
  integer, parameter :: b12$ = 102, t12$ = 102
  integer, parameter :: b13$ = 103, t13$ = 103
  integer, parameter :: b14$ = 104, t14$ = 104
  integer, parameter :: b15$ = 105, t15$ = 105
  integer, parameter :: b16$ = 106, t16$ = 106
  integer, parameter :: b17$ = 107, t17$ = 107
  integer, parameter :: b18$ = 108, t18$ = 108
  integer, parameter :: b19$ = 109, t19$ = 109
  integer, parameter :: b20$ = 110, t20$ = 110 ! this is n_attrib_special_maxx 

  integer, parameter :: n_attrib_special_maxx = 110

  integer, parameter :: particle$ = 1
  integer, parameter :: n_part$    = 3

  character(16), parameter :: null_name = 'NULL' 
  character(16), parameter :: blank_name = ' '

! electron/positron

  integer, parameter :: proton$     = +2
  integer, parameter :: positron$   = +1
  integer, parameter :: electron$   = -1
  integer, parameter :: antiproton$ = -2

  character(16) :: particle_name(-2:2) = (/ 'ANTIPROTON', &
                  'ELECTRON  ', '???       ', 'POSITRON  ', 'PROTON    ' /)

  integer, parameter :: charge_of(-2:2) = (/ -1, -1, 0, 1, 1 /)
  real(rp), parameter :: mass_of(-2:2) = (/ m_proton, m_electron, 0.0_rp, &
                                            m_electron, m_proton /)

! lattice logical names

  integer, parameter :: linac_lattice$ = 10
  integer, parameter :: linear_lattice$ = 11
  integer, parameter :: circular_lattice$ = 12

  character(16) :: lattice_type(10:12) = &
        (/ 'LINAC_LATTICE   ', 'LINEAR_LATTICE  ', 'CIRCULAR_LATTICE' /)

! logicals for MAKE_HYBIRD_RING

  logical, parameter :: remove_markers$ = .true., no_remove_markers$ = .false.

! control element logicals

  integer, parameter :: free$ = 1, super_slave$ = 2, overlay_slave$ = 3
  integer, parameter :: group_lord$ = 4, super_lord$ = 5, overlay_lord$ = 6
  integer, parameter :: i_beam_lord$ = 7

  character(16) :: control_name(8) = (/ &
            'FREE_ELEMENT   ', 'SUPER_SLAVE    ', 'OVERLAY_SLAVE  ', &
            'GROUP_LORD     ', 'SUPER_LORD     ', 'OVERLAY_LORD   ', &
            'I_BEAM_LORD    ', '               ' /)

! plane list, etc

  integer, parameter :: x_plane$ = 1, y_plane$ = 2
  integer, parameter :: z_plane$ = 3, n_plane$ = 4, s_plane$ = 5

  character(16) :: plane_name(6) = (/ 'X', 'Y', 'Z', 'N', 'S', ' ' /)

  logical, parameter :: set$ = .true., unset$ = .false.

! garbage$ is, for example, for subroutines that want to communicate to
! the calling program that a variable has not been set properly.

  integer, parameter :: garbage$ = -9876

! Note: custom$ = 7, and taylor$ = 8 are taken from the element key list.

  integer, parameter :: bmad_standard$ = 1, symp_lie_ptc$ = 2
  integer, parameter :: runge_kutta$ = 3 
  integer, parameter :: linear$ = 4, tracking$ = 5, symp_map$ = 6
  integer, parameter :: wiedemann$ = 9, symp_lie_bmad$ = 10, none$ = 11
  integer, parameter :: boris$ = 12, adaptive_boris$ = 13, mad$ = 14

  character(16), parameter :: calc_method_name(0:14) = (/ &
      "GARBAGE!      ", "Bmad_Standard ", "Symp_Lie_PTC  ", "Runge_Kutta   ", &
      "Linear        ", "Tracking      ", "Symp_Map      ", "Custom        ", &
      "Taylor        ", "Wiedemann     ", "Symp_Lie_Bmad ", "None          ", &
      "Boris         ", "Adaptive_Boris", "MAD           " /)

! sbend$ and rbend$ are from key definitions.

  integer, parameter :: map_type$ = 1, periodic_type$ = 3
  character(16), parameter :: sub_key_name(0:18) = (/ "GARBAGE!     ", &
         "Map          ", "SBend        ", "Periodic     ", "GARBAGE!     ", "GARBAGE!     ", &
         "GARBAGE!     ", "GARBAGE!     ", "GARBAGE!     ", "GARBAGE!     ", "GARBAGE!     ", &
         "GARBAGE!     ", "GARBAGE!     ", "GARBAGE!     ", "GARBAGE!     ", "GARBAGE!     ", &
         "GARBAGE!     ", "GARBAGE!     ", "RBend        " /)

! The linac_mode_struct is basically the synchrotron integrals with the
! energy factors thrown in. Useful for linacs.

  type amode_struct
    real(rp) emittance        ! Beam emittance
    real(rp) synch_int(4:5)   ! Synchrotron integrals
    real(rp) j_damp           ! damping partition number
    real(rp) alpha_damp       ! damping per turn
    real(rp) chrom            ! Chromaticity
    real(rp) tune             ! "Fractional" tune in radians
  end type

  type linac_mode_struct
    real(rp) i2_E4        ! Integral: g^2 * gamma^4
    real(rp) i3_E7        ! Integral: g^3 * gamma^7
    real(rp) i5a_E6       ! Integral: (g^3 * H_a) * gamma^6
    real(rp) i5b_E6       ! Integral: (g^3 * H_b) * gamma^6
    real(rp) sig_E1       ! Energy spread after 1 pass (eV)
    real(rp) emittance_a  ! a mode emittance at end of linac
    real(rp) emittance_b  ! b mode emittance at end of linac
  end type

  type modes_struct
    real(rp) synch_int(3)  ! Synchrotron integrals I1, I2, and I3
    real(rp) sigE_E        ! SigmaE/E
    real(rp) sig_z         ! Sigma_Z
    real(rp) e_loss        ! Energy loss / turn (eV)
    type (amode_struct)  a, b, z
    type (linac_mode_struct) lin
  end type

  integer, parameter :: bends$ = 201
  integer, parameter :: wigglers$ = 202
  integer, parameter :: all$ = 203

! common flags
! status structure

  type bmad_status_struct
    integer :: status         = ok$      ! Computation status 
    logical :: ok             = .true.   ! Error flag
    logical :: type_out       = .true.   ! Print error messages?
    logical :: sub_type_out   = .true.   ! 
    logical :: exit_on_error  = .true.   ! Exit program on error?
  end type

  type (bmad_status_struct), save :: bmad_status

!---------------------------------------------------------------------------
! Units

  integer, parameter :: radians$ = 1, degrees$ = 2, cycles$ = 3, kHz$ = 4
  character*8 ::frequency_units_name(4) = (/ &
            'Radians ', 'Degrees ', 'Cycles  ', 'kHz     ' /)

! Electric and magnetic fields.

  integer, parameter :: kick_field$ = 1, em_field$ = 2

  type em_field_struct
    real(rp) E(3)         ! electric field
    real(rp) B(3)         ! magnetic field
    real(rp) kick(3)      ! kick
    real(rp) dE(3,3)      ! electric field gradient
    real(rp) dB(3,3)      ! magnetic field gradient
    real(rp) dkick(3,3)   ! kick gradiant
    integer  type         ! kick_field$ or em_field$
  end type

!------------------------------------------------------------------------------
! common stuff

! %taylor_order_ptc is what ptc has been set to.
! %taylor_order is what the user wants.
! The reason why there are two taylor_orders is that the Taylor order of PTC
!   cannot be set until the energy is set so Bmad must sometimes cache the 
!   taylor order until the energy is known.
! %max_aperture_limit is used when no limit is specified or when 
!   ring%param%aperture_limit_on = False.

  type bmad_com_struct
    real(rp) :: d_orb(6) = 1e-5  ! for the make_mat6_tracking routine
    real(rp) :: max_aperture_limit = 1e3    
    real(rp) :: k_loss = 0                   ! Internal var for LCavities.
#if defined(CESR_F90_DOUBLE)
    real(rp) :: rel_tollerance = 1e-5
    real(rp) :: abs_tollerance = 1e-8
#else
    real(rp) :: rel_tollerance = 1e-3
    real(rp) :: abs_tollerance = 1e-6
#endif
    integer :: taylor_order = 3              ! 3rd order is default
    integer :: taylor_order_ptc = 0          ! 0 -> not yet set 
    logical :: taylor_order_set = .false.    ! Used by set_taylor_order
    integer :: real_8_map_init               ! See PTC doc.
    integer :: default_integ_order = 2       ! PTC integration order
    integer :: default_num_steps = 1         ! Number integration steps
    logical :: canonical_coords = .true.     ! Use (x, px) [not (x, x')]
    logical :: use_liar_lcavity = .false.    ! Liar like tracking?
    logical :: sr_wakes_on = .true.          ! Short range wakefields?
    logical :: lr_wakes_on = .true.          ! Long range wakefields
    logical :: mat6_track_symmetric = .true. ! symmetric offsets
  end type
  
  type (bmad_com_struct), save :: bmad_com

! multi_turn_func_common is for multi_turn_tracking_to_mat.

  type (coord_struct), pointer :: multi_turn_func_common(:) 

end module
