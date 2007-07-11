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

! IF YOU CHANGE THE lat_struct OR ANY ASSOCIATED STRUCTURES YOU MUST 
! INCREASE THE VERSION NUMBER !
! THIS IS USED BY BMAD_PARSER TO MAKE SURE DIGESTED FILES ARE OK.

integer, parameter :: bmad_inc_version$ = 83

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

! Size of ele%value(:) array

integer, parameter :: n_attrib_maxx = 55

! coordinate def

type coord_struct   ! coordinates at a single point
  real(rp) vec(6)   ! (x, p_x, y, p_y, z, p_z)
  complex(rp) spin(2) ! particle spin in spinor notation
end type

type orbit_struct                            ! an entire orbit.
  type (coord_struct), allocatable :: at(:)  ! coords "at" end of each element.
end type 

! Wiggler structures

integer, parameter :: hyper_y$ = 1, hyper_xy$ = 2, hyper_x$ = 3
character(8), parameter :: wig_term_type_name(0:3) = (/ &
                  'Garbage ', 'Hyper_Y ', 'Hyper_XY', 'Hyper_X ' /)

type wig_term_struct
  real(rp) coef
  real(rp) kx, ky, kz
  real(rp) phi_z
  integer type      ! hyper_y$, hyper_xy$, or hyper_x$
end type

! Wakefield structs...
! Each sr_wake_struct represents a point on the wake vs. z curve.

type sr_table_wake_struct    ! Tabular short-Range Wake struct
  real(rp) z            ! Distance behind the leading particle
  real(rp) long         ! Longitudinal wake in V/C/m
  real(rp) trans        ! Transverse wake in V/C/m^2
end type

type sr_mode_wake_struct  ! Psudo-mode short-Range Wake struct 
  real(rp) amp        ! Amplitude
  real(rp) damp       ! Dampling factor.
  real(rp) k          ! k factor
  real(rp) phi        ! Phase in radians
  real(rp) norm_sin   ! non-skew sin-like component of the wake
  real(rp) norm_cos   ! non-skew cos-like component of the wake
  real(rp) skew_sin   ! skew sin-like component of the wake
  real(rp) skew_cos   ! skew cos-like component of the wake
end type

! Each lr_wake_struct represents a different mode.
! A non-zero freq_spread attribute value will make freq different from freq_in.

type lr_wake_struct     ! Long-Range Wake struct.
  real(rp) freq         ! Actual Frequency in Hz.
  real(rp) freq_in      ! Input frequency in Hz.
  real(rp) R_over_Q     ! Strength in V/C/m^2.
  real(rp) Q            ! Quality factor.
  real(rp) angle        ! polarization angle (radians/2pi).
  real(rp) norm_sin     ! non-skew sin-like component of the wake.
  real(rp) norm_cos     ! non-skew cos-like component of the wake.
  real(rp) skew_sin     ! skew sin-like component of the wake.
  real(rp) skew_cos     ! skew cos-like component of the wake.
  integer m             ! Order (1 = dipole, 2 = quad, etc.)
  logical polarized     ! Polaraized mode?
end type

! Note: Bmad routines observe the following rule: 
!   All pointers within a wake_struct are assumed to be allocated.

type wake_struct
  character(200) :: sr_file = ' '
  character(200) :: lr_file = ' '
  type (sr_table_wake_struct), pointer :: sr_table(:) => null()
  type (sr_mode_wake_struct), pointer :: sr_mode_long(:) => null()
  type (sr_mode_wake_struct), pointer :: sr_mode_trans(:) => null()
  type (lr_wake_struct), pointer :: lr(:) => null()
  real(rp) :: z_sr_mode_max = 0   ! Max allowable z value sr_mode. 
end type

! Local reference frame position with respect to the global (floor) coordinates

type floor_position_struct
  real(rp) x, y, z            ! offset from origin
  real(rp) theta, phi, psi    ! angular orientation
end type

! Transverse space charge structure. This structure contains information
! about the beam as a whole.

type trans_space_charge_struct
  type (coord_struct) closed_orb   ! beam orbit
  real(rp) kick_const
  real(rp) sig_x
  real(rp) sig_y
  real(rp) phi      ! Rotation angle to go from lab frame to rotated frame.
  real(rp) sin_phi
  real(rp) cos_phi
  real(rp) sig_z
endtype    

type xy_disp_struct
  real(rp) eta, etap
end type

! Ele_struct:
! REMEMBER: If this struct is changed you have to:
!     Increase bmad_inc_version by 1.
!     Modify read_digested_bmad_file.
!     Modify write_digested_bmad_file.
!     Modify init_ele (in bmad_utils_mod).
!     Modify ele_equal_ele

type ele_struct
  character(40) name                     ! name of element.
  character(40) type                     ! type name.
  character(40) alias                    ! Another name.
  character(40) attribute_name           ! Used by overlays.
  type (twiss_struct)  a, b, z           ! Twiss parameters at end of element
  type (xy_disp_struct) x, y             ! Projected dispersions.
  type (floor_position_struct) floor     ! Global floor position at end of ele.
  real(rp) value(n_attrib_maxx)          ! attribute values.
  real(rp) closed_orb(6)                 ! For Bmad internal use only.
  real(rp) gen0(6)                       ! constant part of the genfield map.
  real(rp) vec0(6)                       ! 0th order transport vector.
  real(rp) mat6(6,6)                     ! 1st order transport matrix.
  real(rp) c_mat(2,2)                    ! 2x2 C coupling matrix
  real(rp) gamma_c                       ! gamma associated with C matrix
  real(rp) s                             ! longitudinal position at the end
  real(rp), pointer :: r(:,:) => null()  ! For general use. Not used by Bmad.
  real(rp), pointer :: a_pole(:) => null()         ! multipole
  real(rp), pointer :: b_pole(:) => null()         ! multipoles
  real(rp), pointer :: const(:) => null()          ! Working constants.
  character(200), pointer :: descrip => null()     ! For general use
  type (genfield), pointer :: gen_field => null()  ! For symp_map$
  type (taylor_struct) :: taylor(6)                ! Taylor terms
  type (wake_struct), pointer :: wake => null()    ! Wakefields
  type (wig_term_struct), pointer :: wig_term(:) => null()   ! Wiggler Coefs
  type (trans_space_charge_struct), pointer :: trans_sc => null()
  integer key                ! key value
  integer sub_key            ! For wigglers: map_type$, periodic_type$
  integer control_type       ! SUPER_SLAVE$, OVERLAY_LORD$, etc.
  integer ix_value           ! Pointer for attribute to control
  integer n_slave            ! Number of slaves
  integer ix1_slave          ! Start index for slave elements
  integer ix2_slave          ! Stop  index for slave elements
  integer n_lord             ! Number of lords
  integer ic1_lord           ! Start index for lord elements
  integer ic2_lord           ! Stop  index for lord elements
  integer ix_pointer         ! For general use. Not used by Bmad.
  integer ixx                ! Index for Bmad internal use
  integer ix_ele             ! Index in lat%ele(:) array
  integer mat6_calc_method   ! bmad_standard$, taylor$, etc.
  integer tracking_method    ! bmad_standard$, taylor$, etc.
  integer field_calc         ! Used with Boris, Runge-Kutta integrators.
  integer num_steps          ! number of slices for DA_maps
  integer integrator_order   ! For Etiennes' PTC: 2, 4, or 6.
  integer ptc_kind           ! For setting the ptc kind type.
  integer taylor_order       ! Order of the taylor series.
  integer aperture_at        ! Aperture location: exit_end$, ...
  integer coupler_at         ! Lcavity coupler location: exit_end$, ...
  logical symplectify        ! Symplectify mat6 matrices.
  logical mode_flip          ! Have the normal modes traded places?
  logical multipoles_on      ! For turning multipoles on/off
  logical map_with_offsets   ! Taylor map calculated with element offsets?
  logical field_master       ! Calculate strength from the field value?
  logical is_on              ! For turning element on/off.
  logical internal_logic     ! For Bmad internal use only.
  logical logic              ! For general use. Not used by Bmad.
  logical on_an_i_beam       ! Have an I_Beam overlay_lord?
  logical csr_calc_on        ! Coherent synchrotron radiation calculation
end type

! struct for element to element control

type control_struct
  real(rp) coef                  ! control coefficient
  integer ix_lord                ! index to lord element
  integer ix_slave               ! index to slave element
  integer ix_attrib              ! index of attribute controlled
end type

! parameter and mode structures

type lat_param_struct
  real(rp) n_part             ! Particles/bunch (for BeamBeam elements).
  real(rp) total_length       ! total_length of lat
  real(rp) growth_rate        ! growth rate/turn if not stable
  real(rp) t1_with_RF(6,6)    ! Full 1-turn matrix with RF on.
  real(rp) t1_no_RF(6,6)      ! Full 1-turn matrix with RF off.
  integer particle            ! positron$, electron$, etc.
  integer ix_lost             ! Index of element particle was lost at.
  integer end_lost_at         ! entrance_end$ or exit_end$
  integer lattice_type        ! linear_lattice$, etc...
  integer ixx                 ! Integer for general use
  logical stable              ! is closed lat stable?
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

! lat_struct
! Remember: If this struct is changed you have to modify:
!     read_digested_bmad_file
!     write_digested_bmad_file
!     transfer_lat_parameters
!     lat_equal_lat

type lat_struct
  character(40) name                  ! Name of lat given by USE statement
  character(40) lattice               ! Lattice
  character(200) input_file_name      ! Name of the lattice input file
  character(80) title                 ! General title
  type (mode_info_struct) a, b, z     ! Tunes, etc.
  type (lat_param_struct) param       ! Parameters
  type (ele_struct)  ele_init         ! For use by any program
  type (ele_struct), pointer ::  ele(:) => null()        ! Array of elements
  type (control_struct), pointer :: control(:) => null() ! control list
  type (coord_struct) beam_start      ! Starting coords
  integer version                     ! Version number
  integer n_ele_track                 ! Number of lat elements to track through.
  integer n_ele_max                   ! Index of last valid element in %ele(:) array
  integer n_control_max               ! Last index used in CONTROL_array
  integer n_ic_max                    ! Last index used in IC_array
  integer input_taylor_order          ! As set in the input file
  integer, pointer :: ic(:) => null() ! Index to %control(:)
  real(rp), pointer :: E_TOT          ! points to lat%ele(0)%value(E_TOT$)
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
integer, parameter :: match$ = 30, monitor$ = 31, instrument$ = 32
integer, parameter :: hkicker$ = 33, vkicker$ = 34, rcollimator$ = 35
integer, parameter :: ecollimator$ = 36, i_beam$ = 37, bend_sol_quad$ = 38
integer, parameter :: def_beam_start$ = 39

integer, parameter :: n_key = 39

character(16) :: key_name(n_key) = (/ &
    'DRIFT        ', 'SBEND        ', 'QUADRUPOLE   ', 'GROUP        ', &
    'SEXTUPOLE    ', 'OVERLAY      ', 'CUSTOM       ', 'TAYLOR       ', &
    'RFCAVITY     ', 'ELSEPARATOR  ', 'BEAMBEAM     ', 'WIGGLER      ', &
    'SOL_QUAD     ', 'MARKER       ', 'KICKER       ', 'HYBRID       ', &
    'OCTUPOLE     ', 'RBEND        ', 'MULTIPOLE    ', 'ACCEL_SOL    ', &
    'DEF_BEAM     ', 'AB_MULTIPOLE ', 'SOLENOID     ', 'PATCH        ', &
    'LCAVITY      ', 'DEF_PARAMETER', 'NULL_ELE     ', 'INIT_ELE     ', &
    'HOM          ', 'MATCH        ', 'MONITOR      ', 'INSTRUMENT   ', &
    'HKICKER      ', 'VKICKER      ', 'RCOLLIMATOR  ', 'ECOLLIMATOR  ', &
    'I_BEAM       ', 'BEND_SOL_QUAD', 'BEAM_START   ' /)

! Attribute name logical definitions
! Note: The following attributes must have unique number assignments:
!     L$, TILT$, X_PITCH$ and higher

integer, parameter :: particle$ = 1, n_part$   = 2, taylor_order$ = 3
integer, parameter :: energy_gev$ = 4, lattice_type$ = 5, symmetry$ = 6

integer, parameter :: x_beg_limit$=2, y_beg_limit$=3, b_x2$=4, &
          b_y2$=5, l_st2$=9, b_z$=10, l_st1$=11, s_st2$=12, s_st1$=13, &
          b_x1$=14, b_y1$=15    

integer, parameter :: val1$=3, val2$=4, val3$=5, val4$=6, val5$=7, &
          val6$=8, val7$=9, val8$=10, val9$=11, val10$=12, val11$=13, &
          val12$=15

integer, parameter :: beta_a0$ = 2, alpha_a0$ = 3, beta_b0$ = 4, &
          alpha_b0$ = 5, beta_a1$ = 6, alpha_a1$ = 7, beta_b1$ = 8, &
          alpha_b1$ = 9, dphi_a$ = 10, dphi_b$ = 11, &
          eta_a0$ = 12, etap_a0$ = 13, eta_b0$ = 14, etap_b0$ = 15, &
          eta_a1$ = 16, etap_a1$ = 17, eta_b1$ = 18, etap_b1$ = 19

integer, parameter :: x$ = 1, p_x$ = 2, y$ = 3, p_y$ = 4, z$ = 5, p_z$ = 6

!  integer, parameter :: x_position$ = 2, y_position$ = 3, z_position$ = 4, &
!          theta_position$ = 5, phi_position$ = 6, psi_position$ = 7, &
!          beta_a$ = 8, beta_b$ = 9, alpha_a$ = 10, alpha_b$ = 11, &
!          eta_a$ = 12, eta_b$ = 13, etap_a$ = 14, etap_b$ = 15, &
!          phase_x$ = 16, phase_y$ = 17, &
!          c11$ = 18, c12$ = 19, c21$ = 20, c22$ = 21

integer, parameter :: l$=1
integer, parameter :: tilt$=2, command$=2
integer, parameter :: old_command$=3, angle$=3, kick$=3, gradient_err$=3
integer, parameter :: k1$=4, sig_x$=4, harmon$=4, h_displace$=4, e_loss$=4
integer, parameter :: k2$=5, sig_y$=5, b_max$=5, v_displace$=5, phi0_err$=5
integer, parameter :: k3$=6, sig_z$=6, rf_wavelength$=6, g_err$=6
integer, parameter ::        dks_ds$=6, lrad$=6   ! lrad is for felv testing.
integer, parameter :: g$=7, ks$=7, voltage$=7, n_pole$=7, bbi_const$=7
integer, parameter :: e1$=8, charge$=8, gap$=8, dphi0$=8
integer, parameter :: n_slice$=9, e2$=9, rf_frequency$=9
integer, parameter :: fint$=10, polarity$=10, gradient$=10
integer, parameter :: fintx$=11, z_patch$=11, phi0$=11
integer, parameter :: rho$=12, s_center$=12, p0c_start$=12
integer, parameter :: hgap$=13, E_TOT_START$=13, x_patch$=13
integer, parameter :: coef$=14, current$=14, hgapx$=14, delta_e$=14, l_pole$=14
integer, parameter :: roll$=15, quad_tilt$=15, freq_spread$=15, x_ray_line_len$=15
integer, parameter :: l_original$=16, l_chord$=16, bend_tilt$=16
integer, parameter :: l_start$=17, h1$=17, x_quad$=17
integer, parameter :: l_end$=18, h2$=18, y_quad$=18
integer, parameter :: x_pitch$=19
integer, parameter :: y_pitch$=20
integer, parameter :: hkick$=21
integer, parameter :: vkick$=22
integer, parameter :: BL_hkick$=23  
integer, parameter :: BL_vkick$=24
integer, parameter :: x_offset$=25
integer, parameter :: y_offset$=26
integer, parameter :: s_offset$=27, z_offset$=27
integer, parameter :: dE_offset$=28, check_sum$=28, B_field_err$=28
integer, parameter :: x_limit$=29
integer, parameter :: y_limit$=30
integer, parameter :: aperture$=31
integer, parameter :: radius$=32
integer, parameter :: E_TOT$=33
integer, parameter :: rel_tol$=34
integer, parameter :: abs_tol$=35
integer, parameter :: B_field$=36, Bs_field$ = 36, E_field$=36
integer, parameter :: B_gradient$=37, B1_gradient$=37, B2_gradient$=37, B3_gradient$=37
integer, parameter :: E_gradient$=37, E1_gradient$=37, E2_gradient$=37, E3_gradient$=37
integer, parameter :: tilt_tot$=38
integer, parameter :: x_pitch_tot$=39
integer, parameter :: y_pitch_tot$=40
integer, parameter :: x_offset_tot$=41
integer, parameter :: y_offset_tot$=42
integer, parameter :: s_offset_tot$=43
integer, parameter :: p0c$ = 44
integer, parameter :: BL_kick$ = 45
integer, parameter :: coupler_strength$ = 46
integer, parameter :: coupler_phase$ = 47
integer, parameter :: coupler_angle$ = 48
integer, parameter :: kick_tilt$ = 49
integer, parameter :: ds_step$ = 50
integer, parameter :: ref_orb$ = 51 ! 51 through 54 are reserved for the reference orbit

integer, parameter :: symmetric_edge$ = 56 ! this is 1 + n_attrib_maxx
integer, parameter :: mat6_calc_method$ = 57
integer, parameter :: tracking_method$  = 58
integer, parameter :: num_steps$ = 59
integer, parameter :: integrator_order$ = 60
integer, parameter :: term$ = 61
integer, parameter :: ptc_kind$ = 62
integer, parameter :: symplectify$ = 63
integer, parameter :: descrip$ = 64
integer, parameter :: is_on$ = 65
integer, parameter :: field_calc$ = 66
integer, parameter :: type$ = 67
integer, parameter :: aperture_at$ = 68
integer, parameter :: ran_seed$ = 69
integer, parameter :: sr_wake_file$ = 70 
integer, parameter :: lr_wake_file$ = 71
integer, parameter :: alias$ =72
integer, parameter :: start_edge$ =73
integer, parameter :: end_edge$ =74
integer, parameter :: accordion_edge$ =75
integer, parameter :: lattice$ = 76
integer, parameter :: coupler_at$ = 77
integer, parameter :: map_with_offsets$ = 78
integer, parameter :: csr_calc_on$ = 79

integer, parameter :: a0$  =  80, k0l$  =  80
integer, parameter :: a20$ = 100, k20l$ = 100

integer, parameter :: b0$  = 110, t0$  = 110
integer, parameter :: b20$ = 130, t20$ = 130 

integer, parameter :: n_attrib_special_maxx = 130

character(40), parameter :: null_name = '!NULL' 
character(40), parameter :: blank_name = ' '

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

integer, parameter :: linear_lattice$ = 10
integer, parameter :: circular_lattice$ = 12

character(16) :: lattice_type(10:12) = &
        (/ 'LINEAR_LATTICE  ', 'GARBAGE!        ', 'CIRCULAR_LATTICE' /)

! logicals for MAKE_HYBIRD_lat

logical, parameter :: remove_markers$ = .true., no_remove_markers$ = .false.

! control element logicals

integer, parameter :: free$ = 1, super_slave$ = 2, overlay_slave$ = 3
integer, parameter :: group_lord$ = 4, super_lord$ = 5, overlay_lord$ = 6
integer, parameter :: i_beam_lord$ = 7, multipass_lord$ = 8, multipass_slave$ = 9

character(16) :: control_name(10) = (/ &
            'FREE_ELEMENT   ', 'SUPER_SLAVE    ', 'OVERLAY_SLAVE  ', &
            'GROUP_LORD     ', 'SUPER_LORD     ', 'OVERLAY_LORD   ', &
            'I_BEAM_LORD    ', 'MULTIPASS_LORD ', 'MULTIPASS_SLAVE', &
            '               ' /)

! plane list, etc

integer, parameter :: x_plane$ = 1, y_plane$ = 2
integer, parameter :: z_plane$ = 3, n_plane$ = 4, s_plane$ = 5

character(16) :: plane_name(6) = (/ 'X', 'Y', 'Z', 'N', 'S', ' ' /)

logical, parameter :: set$ = .true., unset$ = .false.

! garbage$ is, for example, for subroutines that want to communicate to
! the calling program that a variable has not been set properly.

integer, parameter :: int_garbage$ = -9876
real(rp), parameter :: real_garbage$ = -9876.5

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
character(16), parameter :: sub_key_name(0:18) = (/ "GARBAGE!  ", &
     "Map       ", "SBend     ", "Periodic  ", "GARBAGE!  ", "GARBAGE!  ", &
     "GARBAGE!  ", "GARBAGE!  ", "GARBAGE!  ", "GARBAGE!  ", "GARBAGE!  ", &
     "GARBAGE!  ", "GARBAGE!  ", "GARBAGE!  ", "GARBAGE!  ", "GARBAGE!  ", &
     "GARBAGE!  ", "GARBAGE!  ", "RBend     " /)

! ele%aperture_at logical definitions.

integer, parameter :: entrance_end$ = 1, exit_end$ = 2, both_ends$ = 3
integer, parameter :: no_end$ = 4
character(16), parameter :: element_end_name(0:4) = (/ "GARBAGE!    ", &
      "Entrance_End", "Exit_End    ", "Both_Ends   ", "No_End      " /)

! The linac_normal_mode_struct is basically the synchrotron integrals with the
! energy factors thrown in. Useful for linacs.

type anormal_mode_struct
  real(rp) emittance        ! Beam emittance
  real(rp) synch_int(4:5)   ! Synchrotron integrals
  real(rp) j_damp           ! damping partition number
  real(rp) alpha_damp       ! damping per turn
  real(rp) chrom            ! Chromaticity
  real(rp) tune             ! "Fractional" tune in radians
end type

type linac_normal_mode_struct
  real(rp) i2_E4            ! Integral: g^2 * gamma^4
  real(rp) i3_E7            ! Integral: g^3 * gamma^7
  real(rp) i5a_E6           ! Integral: (g^3 * H_a) * gamma^6
  real(rp) i5b_E6           ! Integral: (g^3 * H_b) * gamma^6
  real(rp) sig_E1           ! Energy spread after 1 pass (eV)
  real(rp) a_emittance_end  ! a mode emittance at end of linac
  real(rp) b_emittance_end  ! b mode emittance at end of linac
end type

type normal_modes_struct
  real(rp) synch_int(3)  ! Synchrotron integrals I1, I2, and I3
  real(rp) sigE_E        ! SigmaE/E
  real(rp) sig_z         ! Sigma_Z
  real(rp) e_loss        ! Energy loss / turn (eV)
  real(rp) pz_aperture   ! pz aperture limit
  type (anormal_mode_struct)  a, b, z
  type (linac_normal_mode_struct) lin
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
character(8) ::frequency_units_name(4) = (/ &
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

! Structures for saving the track through an element.
! track%pt(0:n)  goes from 0 to n = track%n_pt

type track_point_struct
  real(rp) s              ! Longitudinal distance from beginning of element.
  type (coord_struct) orb ! position of a point.
  real(rp) vec0(6)        ! 0th order part of xfer map from the beginning.
  real(rp) mat6(6,6)      ! 1st order part of xfer map (transfer matrix).
end type

type track_struct
  type (track_point_struct), pointer :: pt(:) => null() ! An array of track points. 
  real(rp) :: ds_save = 1e-3         ! min distance between points
  real(rp) :: step0 = 1e-3           ! Initial step size.
  real(rp) :: step_min = 1e-8        ! min step size to step below which
                                     !   track1_adaptive_boris gives up
  integer :: max_step = 10000        ! maximum number of steps allowed
  logical :: save_track = .false.    ! save orbit?
  integer :: n_pt                    ! upper bound of track%pt(0:n)
  integer :: n_bad
  integer :: n_ok
end type

!------------------------------------------------------------------------------
! misc

! A real_array_struct is just a pointer to a real number.
! This is used to construct arrays of reals.

type real_array_struct
  real(rp), pointer :: r
end type 

!------------------------------------------------------------------------------
! common stuff

! %taylor_order_ptc is what ptc has been set to.
! %taylor_order is what the user wants.
! The reason why there are two taylor_orders is that the Taylor order of PTC
!   cannot be set until the energy is set so Bmad must sometimes cache the 
!   taylor order until the energy is known.
! %max_aperture_limit is used when no limit is specified or when 
!   lat%param%aperture_limit_on = False.

type bmad_common_struct
  real(rp) :: max_aperture_limit = 1e3       ! Max Aperture.
  real(rp) :: d_orb(6)           = 1e-5      ! for the make_mat6_tracking routine.
  real(rp) :: grad_loss_sr_wake  = 0         ! Internal var for LCavities.
  real(rp) :: default_ds_step    = 0.2       ! Integration step size.  
#if defined(CESR_F90_DOUBLE)
  real(rp) :: rel_tolerance = 1e-5
  real(rp) :: abs_tolerance = 1e-8
#else
  real(rp) :: rel_tolerance = 1e-3
  real(rp) :: abs_tolerance = 1e-6
#endif
  integer :: taylor_order = 3                ! 3rd order is default
  integer :: default_integ_order = 2         ! PTC integration order
  logical :: canonical_coords = .true.       ! Use (x, px) [not (x, x')]
  logical :: use_liar_lcavity = .false.      ! Liar like tracking?
  logical :: sr_wakes_on = .true.            ! Short range wakefields?
  logical :: lr_wakes_on = .true.            ! Long range wakefields
  logical :: mat6_track_symmetric = .true.   ! symmetric offsets
  logical :: auto_bookkeeper = .true.        ! Automatic bookkeeping?
  logical :: trans_space_charge_on = .false. ! Space charge switch
  logical :: coherent_synch_rad_on = .false. ! csr 
  logical :: spin_tracking_on = .false.      ! spin tracking?
  logical :: radiation_damping_on = .false.       ! Damping toggle.
  logical :: radiation_fluctuations_on = .false.  ! Fluctuations toggle.
  logical :: compute_ref_energy = .true.          ! Enable recomputation?
end type
  
type (bmad_common_struct), save :: bmad_com

! multi_turn_func_common is for multi_turn_tracking_to_mat.

type (coord_struct), pointer :: multi_turn_func_common(:) 

end module
