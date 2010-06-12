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

! IF YOU CHANGE THE LAT_STRUCT OR ANY ASSOCIATED STRUCTURES YOU MUST 
! INCREASE THE VERSION NUMBER !!!
! THIS IS USED BY BMAD_PARSER TO MAKE SURE DIGESTED FILES ARE OK.

integer, parameter :: bmad_inc_version$ = 91

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

! Size of ele%value(:) array

integer, parameter :: n_attrib_maxx = 60

! coordinate def

type coord_struct            ! Particle coordinates at a single point
  real(rp) :: vec(6) = 0     ! (x, px, y, py, z, pz)
  complex(rp) :: spin(2) = 0 ! Spin in spinor notation
end type

type coord_array_struct
  type (coord_struct), allocatable :: orb(:)
end type

! Structure for defining beam pipe cross-sections at given longitudinal positions.
! A cross-section is defined by an array of beam_pipe_vertex_struct components. 
! Each beam_pipe_vertex_struct defines a point (vertex) on the pipe.
! Vertices are connected by straight lines or circular arcs.

type beam_pipe_vertex_struct
  real(rp) x, y          ! coordinates of the vertex.
  real(rp) :: radius = 0 ! Radius of arc. 0 => Straight line.
  real(rp) angle         ! angle of (x, y).
end type

! Coupling structure

type bpm_phase_coupling_struct
  real(rp) K_22a  ! In-phase y/x for a-mode oscillations.
  real(rp) K_12a  ! Out-of-phase y/x for a-mode oscillations.
  real(rp) K_11b  ! In-phase x/y for b-mode oscillations.
  real(rp) K_12b  ! Out-of-phase x/y for b-mode oscillations.
  real(rp) Cbar22_a ! Cbar22 as calculated from K_22a.
  real(rp) Cbar12_a ! Cbar12 as calculated from K_12a.
  real(rp) Cbar11_b ! Cbar11 as calculated from K_11b.
  real(rp) Cbar12_b ! Cbar12 as calculated from K_12b.
  real(rp) phi_a    ! a-mode betatron phase.
  real(rp) phi_b    ! b-mode betatron phase.
end type

! Wiggler structures

integer, parameter :: hyper_y$ = 1, hyper_xy$ = 2, hyper_x$ = 3
character(8), parameter :: wig_term_type_name(0:3) = [ &
                  'Garbage ', 'Hyper_Y ', 'Hyper_XY', 'Hyper_X ' ]

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
  real(rp) b_sin      ! non-skew sin-like component of the wake
  real(rp) b_cos      ! non-skew cos-like component of the wake
  real(rp) a_sin      ! skew sin-like component of the wake
  real(rp) a_cos      ! skew cos-like component of the wake
end type

! Each lr_wake_struct represents a different mode.
! A non-zero lr_freq_spread attribute value will make freq different from freq_in.

type lr_wake_struct    ! Long-Range Wake struct.
  real(rp) freq        ! Actual Frequency in Hz.
  real(rp) freq_in     ! Input frequency in Hz.
  real(rp) R_over_Q    ! Strength in V/C/m^2.
  real(rp) Q           ! Quality factor.
  real(rp) angle       ! polarization angle (radians/2pi).
  real(rp) b_sin       ! non-skew sin-like component of the wake.
  real(rp) b_cos       ! non-skew cos-like component of the wake.
  real(rp) a_sin       ! skew sin-like component of the wake.
  real(rp) a_cos       ! skew cos-like component of the wake.
  real(rp) t_ref       ! time reference value for computing the wake amplitude.
                       !  This is used to prevent value overflow with long trains.
  integer m            ! Order (1 = dipole, 2 = quad, etc.)
  logical polarized    ! Polaraized mode?
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

! Structure to be used for an array of pointers to elements.
! The id component is not set by any Bmad routines and can be used by 
! programs that handle multiple lattices to indicate which lattice 
! the element pointer is pointing to.

type ele_pointer_struct
  type (ele_struct), pointer :: ele => null()
  integer id          
end type

! Structure to hold the information of where an individual element is in the lattice.

type lat_ele_loc_struct
  integer ix_ele
  integer ix_branch
end type

! The mode3_struct is used for normal mode analysis of the full 6x6 transfer matrix.

type mode3_struct
  real(rp) v(6,6)
  type (twiss_struct) a, b, c
  type (twiss_struct) x, y
end type

! Ele_struct:
! Remember: If this struct is changed you have to:
!     Increase bmad_inc_version by 1.
!     Modify read_digested_bmad_file.
!     Modify write_digested_bmad_file.
!     Modify init_ele (in bmad_utils_mod).
!     Modify ele_equal_ele

type ele_struct
  character(40) name                   ! name of element.
  character(40) type                   ! type name.
  character(40) alias                  ! Another name.
  character(40) component_name         ! Used by overlays, multipass patch, etc.
  character(200), pointer :: descrip => null() ! Description string.
  type (twiss_struct)  a, b, z         ! Twiss parameters at end of element
  type (xy_disp_struct) x, y           ! Projected dispersions.
  type (floor_position_struct) floor   ! Global floor position at end of ele.
  type (mode3_struct), pointer :: mode3 => null()
  type (coord_struct) map_ref_orb_in   ! Ref orbit at entrance of element.
  type (coord_struct) map_ref_orb_out  ! Ref orbit at exit of element.
  type (genfield), pointer :: gen_field => null() ! For symp_map$
  type (taylor_struct) :: taylor(6)               ! Taylor terms
  type (wake_struct), pointer :: wake => null()   ! Wakefields
  type (wig_term_struct), pointer :: wig_term(:) => null()   ! Wiggler Coefs
  type (trans_space_charge_struct), pointer :: trans_sc => null()
  real(rp) value(n_attrib_maxx)      ! attribute values.
  real(rp) old_value(n_attrib_maxx)  ! Used to see if %value(:) array has changed.
  real(rp) gen0(6)                   ! constant part of the genfield map.
  real(rp) vec0(6)                   ! 0th order transport vector.
  real(rp) mat6(6,6)                 ! 1st order transport matrix.
  real(rp) c_mat(2,2)                ! 2x2 C coupling matrix
  real(rp) gamma_c                   ! gamma associated with C matrix
  real(rp) s                         ! longitudinal position at the exit end.
  real(rp) ref_time                  ! Time ref particle passes exit end.
  real(rp), pointer :: r(:,:) => null()           ! For general use. Not used by Bmad.
  real(rp), pointer :: a_pole(:) => null()        ! multipole
  real(rp), pointer :: b_pole(:) => null()        ! multipoles
  real(rp), pointer :: const(:) => null()         ! Working constants.
  integer key                ! key value
  integer sub_key            ! For wigglers: map_type$, periodic_type$
  integer ix_ele             ! Index in lat%branch(n)%ele(:) array [n = 0 <==> lat%ele(:)].
  integer ix_branch          ! Index in lat%branch(:) array [0 => In lat%ele(:)].
  integer ix_value           ! Overlays: Index of control attribute. 
  integer slave_status       ! super_slave$, etc.
  integer n_slave            ! Number of slaves
  integer ix1_slave          ! Start index for slave elements
  integer ix2_slave          ! Stop  index for slave elements
  integer lord_status        ! overlay_lord$, etc.
  integer n_lord             ! Number of lords
  integer ic1_lord           ! Start index for lord elements
  integer ic2_lord           ! Stop  index for lord elements
  integer ix_pointer         ! For general use. Not used by Bmad.
  integer ixx                ! Index for Bmad internal use
  integer mat6_calc_method   ! bmad_standard$, taylor$, etc.
  integer tracking_method    ! bmad_standard$, taylor$, etc.
  integer field_calc         ! Used with Boris, Runge-Kutta integrators.
  integer num_steps          ! number of slices for DA_maps
  integer integrator_order   ! For Etiennes' PTC: 2, 4, or 6.
  integer ref_orbit          ! Multipass ref orb: single_ref$, match_global_coords$, 
                             !    match_at_entrance$, match_at_exit$, patch_in$, patch_out$
  integer taylor_order       ! Order of the taylor series.
  integer aperture_at        ! Aperture location: exit_end$, ...
  integer aperture_type      ! rectangular$, elliptical$, or star_shape$
  logical symplectify        ! Symplectify mat6 matrices.
  logical mode_flip          ! Have the normal modes traded places?
  logical multipoles_on      ! For turning multipoles on/off
  logical map_with_offsets   ! Taylor map calculated with element offsets?
  logical field_master       ! Calculate strength from the field value?
  logical is_on              ! For turning element on/off.
  logical old_is_on          ! For saving the element on/off state.
  logical logic              ! For general use. Not used by Bmad.
  logical bmad_logic         ! For Bmad internal use only.
  logical on_a_girder        ! Have an Girder overlay_lord?
  logical csr_calc_on        ! Coherent synchrotron radiation calculation
  logical offset_moves_aperture  ! element offsets affects aperture?          
end type

! struct for element to element control

type control_struct
  real(rp) :: coef = 0           ! Control coefficient
  integer :: ix_lord = -1        ! Index to lord element
  integer :: ix_slave = -1       ! Index to slave element
  integer :: ix_branch = 0       ! Index to branch line of slave
  integer :: ix_attrib = 0       ! Index of attribute controlled
end type

! lattice parameter struct

type lat_param_struct
  real(rp) n_part             ! Particles/bunch (for BeamBeam elements).
  real(rp) total_length       ! total_length of lat
  real(rp) unstable_factor    ! growth rate/turn for circular lats. |orbit/limit| for linear lats.
  real(rp) t1_with_RF(6,6)    ! Full 1-turn matrix with RF on.
  real(rp) t1_no_RF(6,6)      ! Full 1-turn matrix with RF off.
  integer particle            ! positron$, electron$, etc.
  integer ix_lost             ! Index of element particle was lost at.
  integer end_lost_at         ! entrance_end$ or exit_end$
  integer plane_lost_at       ! x_plane$, y_plane$, z_plane$ (reversed direction).
  integer lattice_type        ! linear_lattice$, etc...
  integer ixx                 ! Integer for general use
  logical stable              ! is closed lat stable?
  logical aperture_limit_on   ! use apertures in tracking?
  logical lost                ! for use in tracking
end type

!

type mode_info_struct
  real(rp) tune      ! "fractional" tune in radians: 0 < tune < 2pi
  real(rp) emit      ! Emittance.
  real(rp) chrom     ! Chromaticity.
  real(rp) sigma     ! Beam size.
  real(rp) sigmap    ! Beam divergence.
end type

type branch_struct
  character(40) name
  integer key               ! photon_branch$, branch$, etc.
  integer ix_branch
  integer ix_from_branch    ! 0 => main lattice line
  integer ix_from_ele
  integer, pointer :: n_ele_track
  integer, pointer :: n_ele_max
  type (ele_struct), pointer :: ele(:) => null()
  type (lat_param_struct), pointer :: param => null()
end type

type dummy_parameter_struct
  integer dummy(250)
end type

! lat_struct
! Remember: If this struct is changed you have to modify:
!     Increase bmad_inc_version by 1.
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
  type (ele_struct), pointer ::  ele(:) => null()  ! Array of elements [=> branch(0)].
  type (branch_struct), allocatable :: branch(:)   ! Branch arrays
  type (control_struct), allocatable :: control(:) ! Control list
  type (coord_struct) beam_start      ! Starting coords
  integer version                     ! Version number
  integer n_ele_track                 ! Number of lat elements to track through.
  integer n_ele_max                   ! Index of last valid element in %ele(:) array
  integer n_control_max               ! Last index used in control_array
  integer n_ic_max                    ! Last index used in ic_array
  integer input_taylor_order          ! As set in the input file
  integer, allocatable :: ic(:)       ! Index to %control(:)
end type

!

character(2), parameter :: coord_name(6) = &
                              [ "X ", "Px", "Y ", "Py", "Z ", "Pz" ]

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
integer, parameter :: ecollimator$ = 36, girder$ = 37, bend_sol_quad$ = 38
integer, parameter :: def_beam_start$ = 39, photon_branch$ = 40
integer, parameter :: branch$ = 41, mirror$ = 42, crystal$ = 43
integer, parameter :: pipe$ = 44

integer, parameter :: n_key = 44

! "bend_sol_" is used to force the use of at least "bend_sol_q" in defining bend_sol_quad elements

character(16) :: key_name(n_key) = [ &
    'DRIFT        ', 'SBEND        ', 'QUADRUPOLE   ', 'GROUP        ', &
    'SEXTUPOLE    ', 'OVERLAY      ', 'CUSTOM       ', 'TAYLOR       ', &
    'RFCAVITY     ', 'ELSEPARATOR  ', 'BEAMBEAM     ', 'WIGGLER      ', &
    'SOL_QUAD     ', 'MARKER       ', 'KICKER       ', 'HYBRID       ', &
    'OCTUPOLE     ', 'RBEND        ', 'MULTIPOLE    ', 'BEND_SOL_    ', &
    'DEF_BEAM     ', 'AB_MULTIPOLE ', 'SOLENOID     ', 'PATCH        ', &
    'LCAVITY      ', 'DEF_PARAMETER', 'NULL_ELE     ', 'INIT_ELE     ', &
    'HOM          ', 'MATCH        ', 'MONITOR      ', 'INSTRUMENT   ', &
    'HKICKER      ', 'VKICKER      ', 'RCOLLIMATOR  ', 'ECOLLIMATOR  ', &
    'GIRDER       ', 'BEND_SOL_QUAD', 'BEAM_START   ', 'PHOTON_BRANCH', &
    'BRANCH       ', 'MIRROR       ', 'CRYSTAL      ', 'PIPE         ' ]

! These logical arrays get set in init_attribute_name_array and are used
! to sort elements that have kick or orientation attributes from elements that do not.
! The orientation attributes are: tilt, x/y/s_offset, x/y_pitch, and *_tot versions.
! Note: A solenoid does not formally have a tilt but has everything else.
! Rule: Any element that formally has some but not all orientation attributes is considered
!   internally to have all attributes so any routine can safely work with all the 
!   orientation attributes as a block.

logical has_hkick_attributes(n_key)
logical has_kick_attributes(n_key)
logical has_orientation_attributes(n_key)

! Element attribute name logical definitions

integer, parameter :: particle$ = 1, n_part$   = 2, taylor_order$ = 3
integer, parameter :: lattice_type$ = 5, symmetry$ = 6

integer, parameter :: val1$=3, val2$=4, val3$=5, val4$=6, val5$=7, &
          val6$=8, val7$=9, val8$=10, val9$=11, val10$=12, val11$=13, &
          val12$=15

integer, parameter :: beta_a0$ = 2, alpha_a0$ = 3, beta_b0$ = 4, &
          alpha_b0$ = 5, beta_a1$ = 6, alpha_a1$ = 7, beta_b1$ = 8, &
          alpha_b1$ = 9, dphi_a$ = 10, dphi_b$ = 11, &
          eta_x0$ = 12, etap_x0$ = 13, eta_y0$ = 14, etap_y0$ = 15, &
          eta_x1$ = 16, etap_x1$ = 17, eta_y1$ = 18, etap_y1$ = 19, &
          match_end$ = 20, &
          x0$ = 21, px0$ = 22, y0$ = 23, py0$ = 24, z0$ = 25, pz0$ = 26, &
          x1$ = 34, px1$ = 35, y1$ = 36, py1$ = 37, z1$ = 38, pz1$ = 39, &
          match_end_orbit$ = 40, c_11$ = 46, c_12$ = 47, c_21$ = 48, c_22$ = 49, gamma_c$ = 50 

integer, parameter :: x$ = 1, px$ = 2, y$ = 3, py$ = 4, z$ = 5, pz$ = 6

integer, parameter :: l$=1    ! Assumed unique. Do not overload.
integer, parameter :: tilt$=2, command$=2, ix_branch_to$=2
integer, parameter :: old_command$=3, angle$=3, kick$=3, gradient_err$=3, x_gain_err$=3
integer, parameter :: direction$=3, graze_angle$=3
integer, parameter :: k1$=4, sig_x$=4, harmon$=4, h_displace$=4, e_loss$=4, y_gain_err$=4
integer, parameter ::       graze_angle_err$ = 4
integer, parameter :: k2$=5, sig_y$=5, b_max$=5, v_displace$=5, phi0_err$=5, crunch$=5
integer, parameter ::       critical_angle$ = 5
integer, parameter :: k3$=6, sig_z$=6, rf_wavelength$=6, g_err$=6, noise$=6
integer, parameter ::       dks_ds$=6, lrad$=6   ! lrad -> felv testing.
integer, parameter :: g$=7, ks$=7, voltage$=7, n_pole$=7, bbi_const$=7, osc_amplitude$=7
integer, parameter ::       g_graze$ = 7
integer, parameter :: e1$=8, charge$=8, gap$=8, dphi0$=8, x_gain_calib$=8, g_trans$=8
integer, parameter :: n_slice$=9, e2$=9, rf_frequency$=9, y_gain_calib$=9
integer, parameter :: fint$=10, polarity$=10, gradient$=10, crunch_calib$=10
integer, parameter :: fintx$=11, z_patch$=11, phi0$=11, x_offset_calib$=11
integer, parameter :: rho$=12, s_center$=12, p0c_start$=12, y_offset_calib$=12
integer, parameter :: hgap$=13, e_tot_start$=13, tilt_calib$=13
integer, parameter :: coef$=14, current$=14, hgapx$=14, delta_e$=14, l_pole$=14
integer, parameter ::       de_eta_meas$=14
integer, parameter :: roll$=15, quad_tilt$=15, lr_freq_spread$=15, x_ray_line_len$=15
integer, parameter :: n_sample$=15, delta_ref_time$=15
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
integer, parameter :: s_offset$=27, z_offset$=27 ! Assumed unique. Do not overload further.
integer, parameter :: B_field_err$=28, BL_kick$ = 28
integer, parameter :: radius$=29
integer, parameter :: n_ref_pass$=30  ! Assumed unique. Do not overload.
integer, parameter :: tilt_err$=31    
integer, parameter :: p0c$=32         ! Assumed unique. Do not overload.
integer, parameter :: e_tot$=33       ! Assumed unique. Do not overload.
integer, parameter :: Bs_field$=34
integer, parameter :: B_field$=35, E_field$=35
integer, parameter :: B_gradient$=36, E_gradient$=36
integer, parameter :: B1_gradient$=37, E1_gradient$=37
integer, parameter :: B2_gradient$=38, E2_gradient$=38, patch_end$ = 38
integer, parameter :: B3_gradient$=39, E3_gradient$=39, translate_after$=39
integer, parameter :: tilt_tot$=40      
integer, parameter :: x_pitch_tot$=41   
integer, parameter :: y_pitch_tot$=42   
integer, parameter :: x_offset_tot$=43  
integer, parameter :: y_offset_tot$=44  
integer, parameter :: s_offset_tot$=45  
integer, parameter :: coupler_strength$ = 46, Pz_offset$ = 46
integer, parameter :: coupler_phase$ = 47
integer, parameter :: coupler_angle$ = 48
integer, parameter :: pole_radius$ = 49, coupler_at$ = 49
integer, parameter :: ds_step$ = 50
integer, parameter :: general1$ = 51   ! For general use
integer, parameter :: general2$ = 52   ! For general use
integer, parameter :: general3$ = 53   ! For general use
integer, parameter :: general4$ = 54   ! For general use
integer, parameter :: general5$ = 55   ! For general use
integer, parameter :: x1_limit$ = 56   ! Assumed unique. Do not overload.
integer, parameter :: x2_limit$ = 57   ! Assumed unique. Do not overload.
integer, parameter :: y1_limit$ = 58   ! Assumed unique. Do not overload.
integer, parameter :: y2_limit$ = 59   ! Assumed unique. Do not overload.
integer, parameter :: check_sum$ = 60  ! Assumed unique. Do not overload.

!! 61 = 1 + n_attrib_maxx

integer, parameter :: ref_orbit$ = 61, term$ = 61
integer, parameter :: x_position$ = 62
integer, parameter :: symplectify$ = 63, y_position$ = 63
integer, parameter :: descrip$ = 64, z_position$ = 64
integer, parameter :: is_on$ = 65, theta_position$ = 65
integer, parameter :: field_calc$ = 66, phi_position$ = 66
integer, parameter :: type$ = 67, psi_position$ = 67
integer, parameter :: aperture_at$ = 68, beta_a$ = 68
integer, parameter :: ran_seed$ = 69, beta_b$ = 69
integer, parameter :: sr_wake_file$ = 70, alpha_a$ = 70, ref_patch$ = 70
integer, parameter :: lr_wake_file$ = 71, alpha_b$ = 71
integer, parameter :: alias$ =72, eta_x$ = 72
integer, parameter :: start_edge$ =73, eta_y$ = 73
integer, parameter :: end_edge$ =74, etap_x$ = 74
integer, parameter :: accordion_edge$ =75, etap_y$ = 75
integer, parameter :: lattice$ = 76, phi_a$ = 76
integer, parameter :: aperture_type$ = 77, phi_b$ = 77
integer, parameter :: map_with_offsets$ = 78, cmat_11$ = 78
integer, parameter :: csr_calc_on$ = 79, cmat_12$ = 79
integer, parameter :: symmetric_edge$ = 80, cmat_21$ = 80
integer, parameter :: mat6_calc_method$ = 81, cmat_22$ = 81
integer, parameter :: tracking_method$  = 82, s_long$ = 82
integer, parameter :: num_steps$ = 83, ref_time$ = 83
integer, parameter :: integrator_order$ = 84
integer, parameter :: aperture$ = 85
integer, parameter :: x_limit$ = 86
integer, parameter :: y_limit$ = 87
integer, parameter :: offset_moves_aperture$ = 88
integer, parameter :: aperture_limit_on$ = 89

! superimpose$ through common_lord$ assumed unique (or need to modify bmad_parser_mod.f90).

integer, parameter :: superimpose$    = 90   
integer, parameter :: offset$         = 91
integer, parameter :: reference$      = 92
integer, parameter :: ele_beginning$  = 93
integer, parameter :: ele_center$     = 94
integer, parameter :: ele_end$        = 95
integer, parameter :: ref_beginning$  = 96
integer, parameter :: ref_center$     = 97
integer, parameter :: ref_end$        = 98
integer, parameter :: common_lord$    = 99

integer, parameter :: to$ = 100
integer, parameter :: field_master$ = 101
integer, parameter :: star_aperture$ = 102

integer, parameter :: a0$  = 110, k0l$  = 110
integer, parameter :: a20$ = 130, k20l$ = 130

integer, parameter :: b0$  = 140, t0$  = 140
integer, parameter :: b20$ = 160, t20$ = 160 

integer, parameter :: n_attrib_special_maxx = t20$

character(40), parameter :: null_name$ = '!NULL' 
character(40), parameter :: blank_name$ = ' '

! electron/positron

integer, parameter :: proton$     = +2
integer, parameter :: positron$   = +1
integer, parameter :: photon$     =  0
integer, parameter :: electron$   = -1
integer, parameter :: antiproton$ = -2

character(16) :: particle_name(-2:2) = [ 'ANTIPROTON', &
                'ELECTRON  ', 'PHOTON    ', 'POSITRON  ', 'PROTON    ' ]

real(rp), parameter :: charge_of(-2:2) = [-1, -1, 0, 1, 1] * e_charge
real(rp), parameter :: mass_of(-2:2) = [ m_proton, m_electron, 0.0_rp, &
                                            m_electron, m_proton ]

! lattice logical names

integer, parameter :: linear_lattice$ = 10
integer, parameter :: circular_lattice$ = 12

character(16) :: lattice_type(10:12) = &
        [ 'LINEAR_LATTICE  ', 'GARBAGE!        ', 'CIRCULAR_LATTICE' ]

! logicals for MAKE_HYBIRD_lat

logical, parameter :: remove_markers$ = .true., no_remove_markers$ = .false.

! control element logicals

integer, parameter :: free$ = 1, super_slave$ = 2, overlay_slave$ = 3
integer, parameter :: group_lord$ = 4, super_lord$ = 5, overlay_lord$ = 6
integer, parameter :: girder_lord$ = 7, multipass_lord$ = 8, multipass_slave$ = 9
integer, parameter :: not_a_lord$ = 10, group_slave$ = 11, patch_in_slave$ = 12

character(16) :: control_name(12) = [ &
            'FREE           ', 'SUPER_SLAVE    ', 'OVERLAY_SLAVE  ', &
            'GROUP_LORD     ', 'SUPER_LORD     ', 'OVERLAY_LORD   ', &
            'GIRDER_LORD    ', 'MULTIPASS_LORD ', 'MULTIPASS_SLAVE', &
            'NOT_A_LORD     ', 'GROUP_SLAVE    ', 'PATCH_IN_SLAVE ' ]

! plane list, etc

integer, parameter :: x_plane$ = 1, y_plane$ = 2
integer, parameter :: z_plane$ = 3, n_plane$ = 4, s_plane$ = 5

character(16) :: plane_name(6) = [ 'X', 'Y', 'Z', 'N', 'S', ' ' ]

logical, parameter :: set$ = .true., unset$ = .false.

! Note: custom$ = 7, and taylor$ = 8 are taken from the element key list.

integer, parameter :: bmad_standard$ = 1, symp_lie_ptc$ = 2
integer, parameter :: runge_kutta$ = 3 
integer, parameter :: linear$ = 4, tracking$ = 5, symp_map$ = 6
integer, parameter :: symp_lie_bmad$ = 10, no_method$ = 11
integer, parameter :: boris$ = 12, adaptive_boris$ = 13, mad$ = 14

character(16), parameter :: calc_method_name(0:14) = [ &
      "GARBAGE!      ", "Bmad_Standard ", "Symp_Lie_PTC  ", "Runge_Kutta   ", &
      "Linear        ", "Tracking      ", "Symp_Map      ", "Custom        ", &
      "Taylor        ", "GARBAGE!      ", "Symp_Lie_Bmad ", "No_Method     ", &
      "Boris         ", "Adaptive_Boris", "MAD           " ]

! sbend$ and rbend$ are from key definitions.

integer, parameter :: map_type$ = 1, periodic_type$ = 3
character(16), parameter :: sub_key_name(0:18) = [ "GARBAGE!  ", &
     "Map       ", "SBend     ", "Periodic  ", "GARBAGE!  ", "GARBAGE!  ", &
     "GARBAGE!  ", "GARBAGE!  ", "GARBAGE!  ", "GARBAGE!  ", "GARBAGE!  ", &
     "GARBAGE!  ", "GARBAGE!  ", "GARBAGE!  ", "GARBAGE!  ", "GARBAGE!  ", &
     "GARBAGE!  ", "GARBAGE!  ", "RBend     " ]

! ele%aperture_at logical definitions.

integer, parameter :: entrance_end$ = 1, exit_end$ = 2, both_ends$ = 3
integer, parameter :: no_end$ = 4
character(16), parameter :: element_end_name(0:4) = [ "GARBAGE!    ", &
      "Entrance_End", "Exit_End    ", "Both_Ends   ", "No_End      " ]

! ref_orbit values.

integer, parameter :: single_ref$ = 1, match_at_entrance$ = 2, match_at_exit$ = 3 
integer, parameter :: match_global_coords$ = 4, patch_in$ = 5, patch_out$ = 6
character(20), parameter :: ref_orbit_name(0:6) = [      "GARBAGE!           ", &
            "Single_Ref         ", "Match_at_Entrance  ", "Match_at_Exit      ", &
            "Match_Global_Coords", "Patch_In           ", "Patch_Out          " ]

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
  real(rp) synch_int(0:3) ! Synchrotron integrals I0, I1, I2, and I3
  real(rp) sigE_E         ! SigmaE/E
  real(rp) sig_z          ! Sigma_Z
  real(rp) e_loss         ! Energy loss / turn (eV)
  real(rp) pz_aperture    ! pz aperture limit
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
character(8) ::frequency_units_name(4) = [ &
            'Radians ', 'Degrees ', 'Cycles  ', 'kHz     ' ]

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

! This is for debugging radiation damping and fluctuations.

type synch_rad_common_struct
  type (ele_struct), pointer :: ele0    ! Previous element. For i5 calc.
  real(rp) :: scale = 1.0               ! used to scale the radiation
  real(rp) :: i2 = 0, i3 = 0            ! radiation integrals
  real(rp) :: i5a = 0, i5b = 0
  logical :: i_calc_on = .false.        ! For calculating i2 and i3    
end type

type (synch_rad_common_struct), save :: synch_rad_com

integer, parameter :: not_lost$ = -1

integer, parameter :: is_logical$ = 1, is_integer$ = 2, is_real$ = 3, is_name$ = 4

integer, parameter :: rectangular$ = 1, elliptical$ = 2, star_shape$ = 3
character(16) :: shape_name(0:3) = ['garbage!   ', 'Rectangular', 'Elliptical ', 'Star_Shape ']

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
  real(rp) :: default_ds_step    = 0.2_rp    ! Integration step size.  
  real(rp) :: significant_longitudinal_length = 1e-10 ! meter 
  real(rp) :: rel_tolerance = 1e-5
  real(rp) :: abs_tolerance = 1e-8
  real(rp) :: rel_tol_adaptive_tracking = 1e-6  ! Adaptive tracking relative tolerance.
  real(rp) :: abs_tol_adaptive_tracking = 1e-7  ! Adaptive tracking absolute tolerance.
  integer :: taylor_order = 3                ! 3rd order is default
  integer :: default_integ_order = 2         ! PTC integration order
  logical :: canonical_coords = .true.       ! NOT USED.
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
  logical :: conserve_taylor_maps = .true.        ! Enable bookkeeper to set
                                                  ! ele%map_with_offsets = F?
end type
  
type (bmad_common_struct), save :: bmad_com

! multi_turn_func_common is for multi_turn_tracking_to_mat.

type (coord_struct), pointer :: multi_turn_func_common(:) 

! This structure stores the radiation integrals for the individual elements except
! lin_norm_emittance_a and lin_norm_emittance_b are running sums.

type rad_int_common_struct
  real(rp), allocatable :: i0(:)          ! Noe: All arrays are indexed from 0
  real(rp), allocatable :: i1(:)          ! Noe: All arrays are indexed from 0
  real(rp), allocatable :: i2(:) 
  real(rp), allocatable :: i3(:) 
  real(rp), allocatable :: i4a(:)
  real(rp), allocatable :: i4b(:)
  real(rp), allocatable :: i5a(:) 
  real(rp), allocatable :: i5b(:) 
  real(rp), allocatable :: lin_i2_E4(:) 
  real(rp), allocatable :: lin_i3_E7(:) 
  real(rp), allocatable :: lin_i5a_E6(:) 
  real(rp), allocatable :: lin_i5b_E6(:) 
  real(rp), allocatable :: lin_norm_emit_a(:)  ! Running sum
  real(rp), allocatable :: lin_norm_emit_b(:)  ! Running sum
  real(rp), allocatable :: n_steps(:)      ! number of qromb steps needed
end type

end module
