!+
! Bmad_struct holds the structure definitions for Bmad routines.
!-

module bmad_struct

use bmad_taylor_mod
use random_mod
use twiss_mod
use basic_bmad_mod

use definition, only: genfield, fibre, layout

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! IF YOU CHANGE THE LAT_STRUCT OR ANY ASSOCIATED STRUCTURES YOU MUST 
! INCREASE THE VERSION NUMBER !!!
! THIS IS USED BY BMAD_PARSER TO MAKE SURE DIGESTED FILES ARE OK.

integer, parameter :: bmad_inc_version$ = 129

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

! num_ele_attrib$ is size of ele%value(:) array.

integer, parameter :: num_ele_attrib$ = 70

! Misc

integer, parameter :: off$ = 1, on$ = 2

! Diffraction

integer, parameter :: bragg_diffracted$ = 1, forward_diffracted$ = 2, undiffracted$ = 3
character(20), parameter :: ref_orbit_follows_name(0:3) = [character(20) :: 'GARBAGE!', &
                                             'Bragg_Diffracted', 'Forward_Diffracted', 'Undefracted']

! wall3d definitions.

integer, parameter :: anchor_beginning$ = 1, anchor_center$ = 2, anchor_end$ = 3
character(12), parameter :: anchor_pt_name(0:3) = ['GARBAGE! ', 'Beginning', 'Center   ', 'End      ']

! Note: upstream_end$ = entrance_end$ & downstream_end$ = exit_end$

! first_track_edge$ is the edge a particle enters the element at. 
! This edge will depend upon whether a particle is moving in +s or -s direction.
! Similarly, second_track_edge$ is the edge a particle leaves the element at.

integer, parameter :: entrance_end$ = 1, exit_end$ = 2, both_ends$ = 3, no_end$ = 4
integer, parameter :: continuous$ = 5, surface$ = 6
integer, parameter :: first_track_edge$ = 11, second_track_edge$ = 12

character(16), parameter :: aperture_at_name(0:6) = [ &
      'GARBAGE!     ', 'Entrance_End ', 'Exit_End     ', 'Both_Ends    ', &
      'No_End       ', 'Continuous   ', 'Surface      ']

character(16), parameter :: end_at_name(0:4) = [ &
      'GARBAGE!     ', 'Entrance_End ', 'Exit_End     ', 'Both_Ends    ', &
      'No_End       ']

integer, parameter :: upstream_end$ = 1, downstream_end$ = 2
integer, parameter :: inside$ = 3, center_pt$ = 3

character(16), parameter :: ref_pt_name(0:3) = [ &
      'GARBAGE!      ', 'Entrance_End  ', 'Exit_End      ', 'Center        ']

character(16), parameter :: location_name(0:3) = [ &
      'GARBAGE!      ', 'Upstream_End  ', 'Downstream_End', 'Inside        ']

! Structures for defining cross-sections of beam pipes and capillaries
! A cross-section is defined by an array v(:) of wall3d_section_vertex_structs.
! Each vertex v(i) defines a point on the pipe/capillary.
! Vertices are connected by straight lines, circular arcs, or ellipses.
! The radius and tilt values are for the arc from the preceding vertex to this one.
! For v(1), the radius and tilt values are for the arc between v(n) and v(1) where
!   n = upper bound of v(:) array.

integer, parameter :: normal$ = 1, clear$ = 2, opaque$ = 3, crotch$ = 4, crotch1$ = 5
integer, parameter :: crotch2$ = 6, leg1$ = 7, leg2$ = 8
character(16), parameter :: wall3d_section_type_name(8) = [ &
                                                   'Normal  ', 'Clear   ', 'Opaque  ', 'Crotch  ', &
                                                   'Crotch1 ', 'Crotch2 ', 'Leg1    ', 'Leg2    ']
type wall3d_vertex_struct
  real(rp) x, y             ! Coordinates of the vertex.
  real(rp) :: radius_x = 0  ! Radius of arc or ellipse x-axis half width. 0 => Straight line.
  real(rp) :: radius_y = 0  ! Ellipse y-axis half height. 
  real(rp) :: tilt = 0      ! Tilt of ellipse
  real(rp) angle            ! Angle of (x, y) point.
  real(rp) x0, y0           ! Center of ellipse
end type

! A beam pipe or capillary cross section is a collection of vertexes.
! Vertices are always ordered in increasing angle.

type wall3d_section_struct
  integer :: type = normal$                 ! normal$, clear$, opaque$, crotch$, crotch1$, leg2$, ...
  real(rp) :: s = 0                         ! Longitudinal position
  integer :: n_vertex_input = 0             ! Number of vertices specified by the user.
  integer :: ix_ele = 0                     ! index of lattice element containing section
  type (wall3d_vertex_struct), allocatable :: v(:)     ! Array of vertices
  real(rp) :: x0 = 0, y0 = 0                ! Center of section
  ! Section-to-section spline interpolation of the center of the section
  real(rp) :: dx0_ds = 0                    ! Center of wall derivative
  real(rp) :: dy0_ds = 0                    ! Center of wall derivative
  real(rp) :: x0_coef(0:3) = 0              ! Spline coefs for x-center
  real(rp) :: y0_coef(0:3) = 0              ! Spline coefs for y-center
  ! Section-to_section spline interpolation of the wall.
  real(rp) :: dr_ds = real_garbage$         ! derivative of wall radius 
  real(rp) :: p1_coef(3) = 0                ! Spline coefs for p0 function
  real(rp) :: p2_coef(3) = 0                ! Spline coefs for p1 function
end type

! If, say, %ele_anchor_pt = center$ then center of wall is at the center of the element.

type wall3d_struct
  integer :: n_link = 1                           ! For memory management of %section
  real(rp) :: thickness = 0                       ! For diffraction_plate elements
  character(20) :: clear_material = ''            !
  character(20) :: opaque_material = ''           !
  logical :: superimpose = .false.                ! Can overlap another wall
  integer :: ele_anchor_pt = anchor_beginning$    ! anchor_beginning$, anchor_center$, or anchor_end$
  type (wall3d_section_struct), allocatable :: section(:)
end type  

! plane list, etc

integer, parameter :: x_plane$ = 1, y_plane$ = 2
integer, parameter :: z_plane$ = 3, n_plane$ = 4, s_plane$ = 5

character(1), parameter :: plane_name(6) = ['X', 'Y', 'Z', 'N', 'S', ' ']

! coordinate def
! Use coord_state_name for getting the string representation of coord%state

integer, parameter :: moving_forward$ = -9

integer, parameter :: alive$ = 1, lost$ = 2
integer, parameter :: lost_neg_x_aperture$ = 3, lost_pos_x_aperture$ = 4 
integer, parameter :: lost_neg_y_aperture$ = 5, lost_pos_y_aperture$ = 6
integer, parameter :: lost_z_aperture$ = 7

type coord_struct                 ! Particle coordinates at a single point
  real(rp) :: vec(6) = 0          ! (x, px, y, py, z, pz)
  real(rp) :: s = 0               ! Longitudinal position 
  real(rp) :: t = 0               ! Absolute time (not relative to reference).
  complex(rp) :: spin(2) = 0      ! Spin in spinor notation
  real(rp) :: field(2) = 0        ! Photon E-field intensity (x,y).
  real(rp) :: phase(2) = 0        ! Photon E-field phase (x,y)
  real(rp) :: charge = 0          ! macro charge (Coul).
  real(rp) :: p0c = 0             ! For non-photons: Reference momentum.
                                  !     For photons: Photon momentum (not reference).
  real(rp) :: beta = -1           ! Velocity / c_light.
  integer :: ix_ele = -1          ! Index of element particle was tracked through.
                                  !   May be -1 if element is not associated with a lattice.
  integer :: state = not_set$     ! alive$, lost$, lost_neg_x_aperture$, etc.
  integer :: direction = 1        ! Longitudinal direction of motion
  integer :: species = not_set$   ! Species being tracked.
  integer :: location = upstream_end$  ! upstream_end$, inside$, or downstream_end$
end type

type coord_array_struct
  type (coord_struct), allocatable :: orb(:)
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

! Wiggler 

integer, parameter :: hyper_y$ = 1, hyper_xy$ = 2, hyper_x$ = 3
character(8), parameter :: wig_term_type_name(0:3) = ['Garbage ', 'Hyper_Y ', 'Hyper_XY', 'Hyper_X ']

! Single wiggler term

type wig_term_struct
  real(rp) coef
  real(rp) kx, ky, kz
  real(rp) phi_z
  integer type      ! hyper_y$, hyper_xy$, or hyper_x$
end type

! Wiggler field

type wig_struct
  integer :: n_link = 1                            ! For memory management of %term
  type (wig_term_struct), allocatable :: term(:)   ! Wiggler Coefs
end type

! Wakefield structs...
! Each sr_wake_struct represents a point on the wake vs. z curve.

type rf_wake_sr_table_struct    ! Tabular short-Range Wake struct
  real(rp) z            ! Distance behind the leading particle
  real(rp) long         ! Longitudinal wake in V/C/m
  real(rp) trans        ! Transverse wake in V/C/m^2
end type

type rf_wake_sr_mode_struct  ! Psudo-mode short-Range Wake struct 
  real(rp) amp        ! Amplitude
  real(rp) damp       ! Dampling factor.
  real(rp) k          ! k factor
  real(rp) phi        ! Phase in radians
  real(rp) b_sin      ! non-skew sin-like component of the wake
  real(rp) b_cos      ! non-skew cos-like component of the wake
  real(rp) a_sin      ! skew sin-like component of the wake
  real(rp) a_cos      ! skew cos-like component of the wake
end type

! Each rf_wake_lr_struct represents a different mode.
! A non-zero lr_freq_spread attribute value will make freq different from freq_in.

type rf_wake_lr_struct    ! Long-Range Wake struct.
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
!   All pointers within a rf_wake_struct are assumed to be allocated.

type rf_wake_struct
  character(200) :: sr_file = ' '
  character(200) :: lr_file = ' '
  type (rf_wake_sr_table_struct), pointer :: sr_table(:) => null()
  type (rf_wake_sr_mode_struct), pointer :: sr_mode_long(:) => null()
  type (rf_wake_sr_mode_struct), pointer :: sr_mode_trans(:) => null()
  type (rf_wake_lr_struct), pointer :: lr(:) => null()
  real(rp) :: z_sr_mode_max = 0   ! Max allowable z value sr_mode. 
end type

type em_field_map_term_struct
  complex(rp) :: e_coef = 0
  complex(rp) :: b_coef = 0
end type

type em_field_map_struct
  character(200) :: file = ''   ! Input file name. Used also as ID for instances. 
  integer :: n_link = 1         ! For memory management of this structure
  integer :: ele_anchor_pt = anchor_beginning$
                                ! anchor_beginning$, anchor_center$, or anchor_end$
  real(rp) :: dz = 0            ! Distance between sampled field points.
  type (em_field_map_term_struct), allocatable :: term(:)
end type

type em_field_grid_pt_struct
  complex(rp) :: E(3) = 0
  complex(rp) :: B(3) = 0
end type

type em_field_grid_struct
  character(200) :: file = ''   ! Input file name. Used also as ID for instances. 
  integer :: type = 0           ! Type of grid structure
  integer :: ele_anchor_pt = anchor_beginning$  
                                ! anchor_beginning$, anchor_center$, or anchor_end$
  integer :: n_link = 1         ! For memory management of this structure
  type (em_field_grid_pt_struct), allocatable :: pt(:,:,:)
  real(rp) :: dr(3) = 0   ! Grid spacing.
  real(rp) :: r0(3) = 0   ! Grid origin.
end type

! Electro-Magnetic mode structure
! Note: Unlike most everything else, to save on space, different ele%field%mode%grid
! and ele%field%mode%map pointers may point to the same memeory location. 
! This being the case, these components are never deallocated.
! Rule: If %map is associated then %map%term(:) will be allocated.
! Rule: If %grid is associated then %grid%pt(:,:,:) will be allocated.

type em_field_mode_struct
  integer m                     ! Mode varies as cos(m*phi - phi0_azimuth)
  integer :: harmonic = 0       ! Harmonic of fundamental
  real(rp) :: f_damp = 0        ! 1/Q damping factor 
  real(rp) :: dphi0_ref = 0     ! Mode oscillates as: twopi * (f * t + dphi0_ref)
  real(rp) :: stored_energy = 0 ! epsilon_0/2 * \int_vol |E|^2 [Joules]
  real(rp) :: phi0_azimuth = 0  ! Azimuthal orientation of mode.
  real(rp) :: field_scale = 1   ! Factor to scale the fields by
  integer :: master_scale = 0   ! Master scaling parameter in ele%value(:) array.
  type (em_field_map_struct), pointer :: map => null()
  type (em_field_grid_struct), pointer :: grid => null()
end type

! The RF field may be characterized by a collection of modes.

type em_fields_struct
  type (em_field_mode_struct), allocatable :: mode(:)
end type

! Local reference frame position with respect to the global (floor) coordinates

type floor_position_struct
  real(rp) :: r(3) = 0                      ! (x, y, z) offset from origin
  real(rp) :: theta = 0, phi = 0, psi = 0   ! angular orientation
end type

! Space charge structure. This structure contains information about the beam as a whole.

type space_charge_struct
  type (coord_struct) closed_orb   ! beam orbit
  real(rp) kick_const
  real(rp) sig_x
  real(rp) sig_y
  real(rp) phi      ! Rotation angle to go from lab frame to rotated frame.
  real(rp) sin_phi
  real(rp) cos_phi
  real(rp) sig_z
end type    

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

integer, parameter :: super_ok$ = 0, stale$ = 2
integer, parameter :: attribute_group$ = 1, control_group$ = 2, floor_position_group$ = 3
integer, parameter :: s_position_group$ = 4, ref_energy_group$ = 5, mat6_group$ = 6
integer, parameter :: rad_int_group$ = 7, all_groups$ = 8

! The bookkeeping_state_struct is used for keeping track of what bookkeeping has
! been done on an element. NOTE: The information in this structure is ignored if 
! bmad_com%auto_bookkeeper = True which is the default.
! See the Bmad manual for more details.

type bookkeeping_state_struct
  integer :: attributes = stale$      ! Element dependent attributes: super_ok$, ok$ or stale$
  integer :: control = stale$         ! Lord/slave bookkeeping status: super_ok$, ok$ or stale$ 
  integer :: floor_position = stale$  ! Global (floor) geometry: super_ok$, ok$ or stale$
  integer :: s_position = stale$      ! Longitudinal position & element length: super_ok$, ok$ or stale$
  integer :: ref_energy = stale$      ! Reference energy and ref time: super_ok$, ok$ or stale$
  integer :: mat6 = stale$            ! Linear transfer map status: super_ok$, ok$ or stale$
  integer :: rad_int = stale$         ! Radiation integrals cache status
  integer :: ptc = stale$             ! Associated PTC fibre (or layout) status.
end type

! radiation integral data cache

type rad_int_ele_cache_struct
  real(rp) orb0(6)        ! Reference orbit for the calculation
  real(rp) g2_0           ! g2 factor when orbit = %vec0
  real(rp) g3_0           ! g3 factor when orbit = %vec0
  real(rp) dg2_dorb(6)    ! variation of g2 with respect to orbit.
  real(rp) dg3_dorb(6)    ! Variation of g3 with respect to orbit.
  logical :: stale = .true.
end type 

! Structure for surfaces of mirrors, crystals, etc.
! Rule: This structure is always allocated in the ele_struct for elements that need it.

type surface_grid_pt_struct
  real(rp) x_pitch, y_pitch, x_pitch_rms, y_pitch_rms
end type

integer, parameter :: segmented$ = 2, h_misalign$ = 3
character(16), parameter :: surface_grid_type_name(0:3) = &
                          ['GARBAGE!  ', 'Off       ', 'Segmented ', 'H_Misalign']

type surface_grid_struct
  character(200) :: file = ''
  integer :: type = off$   ! or segmented$, or h_misalign$
  real(rp) :: dr(2) = 0, r0(2) = 0
  type (surface_grid_pt_struct), allocatable :: pt(:,:)
end type

! Scratch space for segmented surface calculations

type segmented_surface_struct
  integer :: ix = int_garbage$, iy = int_garbage$    ! Index of segment
  real(rp) :: x0 = 0, y0 = 0, z0 = 0         ! Center of segment
  real(rp) :: slope_x = 0, slope_y = 0  ! Slopes of segment
end type

! Main surface container structure

type photon_surface_struct
  type (surface_grid_struct) :: grid = surface_grid_struct('', off$, 0, 0, null())
  type (segmented_surface_struct) :: segment = segmented_surface_struct()
  real(rp) :: curvature_xy(0:6,0:6) = 0
  logical :: has_curvature = .false.
end type

! Ele_struct:
! Remember: If this struct is changed you have to:
!     Increase bmad_inc_version by 1.
!     Modify read_digested_bmad_file.
!     Modify write_digested_bmad_file.
!     Modify init_ele (in bmad_utils_mod).
!     Modify ele_equal_ele

type ele_struct
  character(40) name                           ! name of element.
  character(40) type                           ! type name.
  character(40) alias                          ! Another name.
  character(40) component_name                 ! Used by overlays, multipass patch, etc.
  character(200), pointer :: descrip => null() ! Description string.
  type (twiss_struct)  a, b, z                 ! Twiss parameters at end of element
  type (xy_disp_struct) x, y                   ! Projected dispersions.
  type (bookkeeping_state_struct) bookkeeping_state     ! Element attribute bookkeeping
  type (branch_struct), pointer :: branch => null()     ! Pointer to branch containing element.
  type (em_fields_struct), pointer :: em_field => null()! DC and AC E/M fields
  type (floor_position_struct) floor                    ! Reference position in global coords.
  type (ele_struct), pointer :: lord => null()          ! Pointer to a slice lord.
  type (mode3_struct), pointer :: mode3 => null()       ! 6D normal mode structure.
  type (fibre), pointer :: ptc_fibre => null()          ! PTC tracking.
  type (genfield), pointer :: ptc_genfield => null()    ! For symp_map$
  type (rad_int_ele_cache_struct), pointer :: rad_int_cache => null() 
                                                        ! Radiation integral calc cached values 
  type (rf_wake_struct), pointer :: rf_wake => null()   ! Wakes
  type (space_charge_struct), pointer :: space_charge => null()
  type (photon_surface_struct), pointer :: surface => null()
  type (taylor_struct) :: taylor(6)                     ! Taylor terms
  type (wall3d_struct), pointer :: wall3d => null()     ! Chamber or capillary wall
  type (wig_struct), pointer :: wig => null()    ! Wiggler field
  type (coord_struct) map_ref_orb_in     ! Transfer map ref orbit at entrance end of element.
  type (coord_struct) map_ref_orb_out    ! Transfer map ref orbit at exit end of element.
  type (coord_struct) time_ref_orb_in    ! Reference orbit at entrance end for ref_time calc.
  type (coord_struct) time_ref_orb_out   ! Reference orbit at exit end for ref_time calc.
  real(rp) value(num_ele_attrib$)                ! attribute values.
  real(rp) old_value(num_ele_attrib$)            ! Used to see if %value(:) array has changed.
  real(rp) gen0(6)                               ! constant part of the genfield map.
  real(rp) vec0(6)                               ! 0th order transport vector.
  real(rp) mat6(6,6)                             ! 1st order transport matrix.
  real(rp) c_mat(2,2)                            ! 2x2 C coupling matrix
  real(rp) gamma_c                               ! gamma associated with C matrix
  real(rp) s                                     ! longitudinal ref position at the exit end.
  real(rp) ref_time                              ! Time ref particle passes exit end.
  real(rp), pointer :: r(:,:,:) => null()        ! For general use. Not used by Bmad.
  real(rp), pointer :: a_pole(:) => null()       ! knl for multipole elements.
  real(rp), pointer :: b_pole(:) => null()       ! tilt for multipole elements.
  integer key                    ! key value 
  integer sub_key                ! For wigglers: map_type$, periodic_type$
  integer :: ix_ele = -1         ! Index in lat%branch(n)%ele(:) array [n = 0 <==> lat%ele(:)].
  integer ix_branch              ! Index in lat%branch(:) array [0 => In lat%ele(:)].
  integer ix_value               ! Overlays: Index of control attribute. 
  integer slave_status           ! super_slave$, etc.
  integer n_slave                ! Number of slaves
  integer ix1_slave              ! Start index for slave elements
  integer ix2_slave              ! Stop  index for slave elements
  integer lord_status            ! overlay_lord$, etc.
  integer n_lord                 ! Number of lords
  integer ic1_lord               ! Start index for lord elements
  integer ic2_lord               ! Stop  index for lord elements
  integer ix_pointer             ! For general use. Not used by Bmad.
  integer ixx, iyy               ! Index for Bmad internal use
  integer mat6_calc_method       ! bmad_standard$, taylor$, etc.
  integer tracking_method        ! bmad_standard$, taylor$, etc.
  integer spin_tracking_method   ! bmad_standard$, symp_lie_ptc$, etc.
  integer ptc_integration_type   ! drift_kick$, matrix_kick$, or ripken_kick$
  integer field_calc             ! bmad_standard$, grid$, refer_to_lords$, or custom$
  integer aperture_at            ! Aperture location: downstream_end$, ...
  integer aperture_type          ! rectangular$, elliptical$, wall3d$, or custom$
  integer orientation            ! -1 -> Element is longitudinally reversed. +1 -> Normal.
  logical symplectify            ! Symplectify mat6 matrices.
  logical mode_flip              ! Have the normal modes traded places?
  logical multipoles_on          ! For turning multipoles on/off
  logical scale_multipoles       ! Are ab_multipoles within other elements (EG: quads, etc.) 
                                 !   scaled by the strength of the element?
  logical map_with_offsets       ! Taylor map calculated with element offsets?
  logical field_master           ! Calculate strength from the field value?
  logical is_on                  ! For turning element on/off.
  logical old_is_on              ! For saving the element on/off state.
  logical logic                  ! For general use. Not used by Bmad.
  logical bmad_logic             ! For Bmad internal use only.
  logical csr_calc_on            ! Coherent synchrotron radiation calculation
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

! lat_param_struct should be called branch_param_struct [Present name is historical artifact.]

type lat_param_struct
  real(rp) :: n_part = 0                     ! Particles/bunch (for BeamBeam elements).
  real(rp) :: total_length = 0               ! total_length of branch
  real(rp) :: unstable_factor = 0            ! growth rate/turn for closed branches. 
                                             ! |orbit/limit| for open branches.
  real(rp) :: t1_with_RF(6,6) = 0            ! Full 1-turn matrix with RF on.
  real(rp) :: t1_no_RF(6,6) = 0              ! Full 1-turn matrix with RF off.
  real(rp) :: rel_tracking_charge = 1        ! Charge relative to referece charge
  integer :: particle = positron$            ! Reference particle: positron$, electron$, etc.
  integer :: geometry = 0                    ! open$ or closed$
  integer :: ixx = 0                         ! Integer for general use
  logical :: stable = .false.                ! is closed lat stable?
  logical :: aperture_limit_on = .true.      ! use apertures in tracking?
  logical :: reverse_time_tracking = .false. ! Internal variable. Do not set.  
  type (bookkeeping_state_struct) :: bookkeeping_state = bookkeeping_state_struct()
                                          ! Overall status for the branch.
end type

! Structure for linking a branch_struct with a collection of ptc layouts

type ptc_layout_pointer_struct
  type (layout), pointer :: ptr => null()
end type

type ptc_branch1_info_struct
  type (ptc_layout_pointer_struct), allocatable :: layout(:)
end type

!

type mode_info_struct
  real(rp) tune      ! "fractional" tune in radians: 0 < tune < 2pi
  real(rp) emit      ! Emittance.
  real(rp) chrom     ! Chromaticity.
  real(rp) sigma     ! Beam size.
  real(rp) sigmap    ! Beam divergence.
end type

type normal_form_struct
  type (taylor_struct) :: M(6)             ! One-turn taylor map: M = A o R o A_inv, R = exp(:h:)
  type (taylor_struct) :: A(6)             ! Map from Floquet -> Lab coordinates
  type (taylor_struct) :: A_inv(6)         ! Map from Lab -> Floquet coordinates
  type (taylor_struct) :: dhdj(6)          ! Nonlinear tune function operating on Floquet coordinates
  type (ele_struct), pointer :: ele_origin ! Element at which the on-turn map was created.
                                           ! See subroutine: normal_form_taylors
end type

!

type branch_struct
  character(40) :: name = ''
  integer :: ix_branch = -1        ! Index of this branch. 0 => Main branch
  integer :: ix_root_branch = -1   ! Root branch index for this machine.
  integer :: ix_from_branch = -1   ! -1 => Not connected/
  integer :: ix_from_ele = -1      ! Branch ele in from_branch index.
  integer, pointer :: n_ele_track => null()
  integer, pointer :: n_ele_max => null()
  type (lat_struct), pointer :: lat => null()
  type (mode_info_struct), pointer :: a => null(), b => null(), z => null()
  type (ele_struct), pointer :: ele(:) => null()
  type (lat_param_struct), pointer :: param => null()
  type (wall3d_struct), pointer :: wall3d => null()
  type (ptc_branch1_info_struct) ptc
  type (normal_form_struct) normal_form_with_rf, normal_form_no_rf
end type

integer, parameter :: opal$ = 1, impactt$ = 2
character(16) :: pre_tracker_name(0:2) = ['NONE   ', 'OPAL   ', 'IMPACTT']

type pre_tracker_struct
  integer :: who = 0   ! Can be opal$, or impactt$
  integer :: ix_ele_start = 0
  integer :: ix_ele_end = 0
  character(200) :: input_file = ''
end type

! lat_struct
! Remember: If this struct is changed you have to modify:
!     Increase bmad_inc_version by 1.
!     read_digested_bmad_file
!     write_digested_bmad_file
!     transfer_lat_parameters
!     lat_equal_lat

type lat_struct
  character(40) use_name                      ! Name of lat given by USE statement
  character(40) lattice                       ! Lattice
  character(200) input_file_name              ! Name of the lattice input file
  character(80) title                         ! General title
  character(60), allocatable :: attribute_alias(:)  ! Aliases for custom1$, etc.
  type (mode_info_struct) a, b, z             ! Tunes, etc.
  type (lat_param_struct) param               ! Parameters
  type (bookkeeping_state_struct) lord_state  ! lord bookkeeping status.
  type (ele_struct)  ele_init                 ! For use by any program
  type (ele_struct), pointer ::  ele(:) => null()  ! Array of elements [=> branch(0)].
  type (branch_struct), allocatable :: branch(:)   ! Branch arrays
  type (control_struct), allocatable :: control(:) ! Control list
  type (coord_struct) beam_start          ! Starting coords
  type (pre_tracker_struct) pre_tracker   ! For OPAL/IMPACT-T
  integer version                         ! Version number
  integer n_ele_track                     ! Number of lat elements to track through.
  integer n_ele_max                       ! Index of last valid element in %ele(:) array
  integer n_control_max                   ! Last index used in control_array
  integer n_ic_max                        ! Last index used in ic_array
  integer input_taylor_order              ! As set in the input file
  integer, allocatable :: ic(:)           ! Index to %control(:)
  logical absolute_time_tracking          ! Use absolute time in lcavity and rfcavity tracking?
  logical rf_auto_scale_phase             ! See rf_auto_scale_phase_and_amp routine.
  logical rf_auto_scale_amp               ! See rf_auto_scale_phase_and_amp routine.
  logical use_ptc_layout                  ! Use ptc layout for lattice
end type

character(2), parameter :: coord_name(6) = ['X ', 'Px', 'Y ', 'Py', 'Z ', 'Pz']

! KEY value definitions
! Note: sbend$ and rbend$ also used for sub_key

integer, parameter :: drift$ = 1, sbend$ = 2, quadrupole$ = 3, group$ = 4
integer, parameter :: sextupole$ = 5, overlay$ = 6, custom$ = 7, taylor$ = 8
integer, parameter :: rfcavity$ = 9
integer, parameter :: elseparator$ = 10, beambeam$ = 11, wiggler$ = 12
integer, parameter :: sol_quad$ = 13, marker$ = 14, kicker$ = 15
integer, parameter :: hybrid$ = 16, octupole$ = 17, rbend$ = 18
integer, parameter :: multipole$ = 19, key_dummy$ = 20
integer, parameter :: def_beam$ = 21, ab_multipole$ = 22, solenoid$ = 23
integer, parameter :: patch$ = 24, lcavity$ = 25, def_parameter$ = 26
integer, parameter :: null_ele$ = 27, init_ele$ = 28, hom$ = 29
integer, parameter :: match$ = 30, monitor$ = 31, instrument$ = 32
integer, parameter :: hkicker$ = 33, vkicker$ = 34, rcollimator$ = 35
integer, parameter :: ecollimator$ = 36, girder$ = 37, bend_sol_quad$ = 38
integer, parameter :: def_beam_start$ = 39, photon_branch$ = 40
integer, parameter :: branch$ = 41, mirror$ = 42, crystal$ = 43
integer, parameter :: pipe$ = 44, capillary$ = 45, multilayer_mirror$ = 46
integer, parameter :: e_gun$ = 47, em_field$ = 48, floor_shift$ = 49, fiducial$ = 50
integer, parameter :: undulator$ = 51, diffraction_plate$ = 52

integer, parameter :: n_key$ = 52

! "bend_sol_" is used to force the use of at least "bend_sol_q" in defining bend_sol_quad elements

character(40), parameter :: key_name(n_key$) = [ &
    'DRIFT            ', 'SBEND            ', 'QUADRUPOLE       ', 'GROUP            ', &
    'SEXTUPOLE        ', 'OVERLAY          ', 'CUSTOM           ', 'TAYLOR           ', &
    'RFCAVITY         ', 'ELSEPARATOR      ', 'BEAMBEAM         ', 'WIGGLER          ', &
    'SOL_QUAD         ', 'MARKER           ', 'KICKER           ', 'HYBRID           ', &
    'OCTUPOLE         ', 'RBEND            ', 'MULTIPOLE        ', 'BEND_SOL_        ', &
    'DEF_BEAM         ', 'AB_MULTIPOLE     ', 'SOLENOID         ', 'PATCH            ', &
    'LCAVITY          ', 'DEF_PARAMETER    ', 'NULL_ELE         ', 'INIT_ELE         ', &
    'HOM              ', 'MATCH            ', 'MONITOR          ', 'INSTRUMENT       ', &
    'HKICKER          ', 'VKICKER          ', 'RCOLLIMATOR      ', 'ECOLLIMATOR      ', &
    'GIRDER           ', 'BEND_SOL_QUAD    ', 'DEF_BEAM_START   ', 'PHOTON_BRANCH    ', &
    'BRANCH           ', 'MIRROR           ', 'CRYSTAL          ', 'PIPE             ', &
    'CAPILLARY        ', 'MULTILAYER_MIRROR', 'E_GUN            ', 'EM_FIELD         ', &
    'FLOOR_SHIFT      ', 'FIDUCIAL         ', 'UNDULATOR        ', 'DIFFRACTION_PLATE']

! These logical arrays get set in init_attribute_name_array and are used
! to sort elements that have kick or orientation attributes from elements that do not.
! The orientation attributes are: tilt, x/y/z_offset, x/y_pitch, and *_tot versions.
! Note: A solenoid does not formally have a tilt but has everything else.
! Rule: Any element that formally has some but not all orientation attributes is considered
!   internally to have all attributes so any routine can safely work with all the 
!   orientation attributes as a block.

logical has_hkick_attributes(n_key$)
logical has_kick_attributes(n_key$)

! Element attribute name logical definitions

integer, parameter :: n_part$ = 2, taylor_order$ = 3, particle$ = 14
integer, parameter :: geometry$ = 15, lattice_type$ = 16, symmetry$ = 6

integer, parameter :: val1$=3, val2$=4, val3$=5, val4$=6, val5$=7, &
          val6$=9, val7$=10, val8$=11, val9$=12, val10$=13, val11$=14, &
          val12$=15

integer, parameter :: beta_a0$ = 2, alpha_a0$ = 3, beta_b0$ = 4, &
          alpha_b0$ = 5, beta_a1$ = 6, alpha_a1$ = 7, beta_b1$ = 8, &
          alpha_b1$ = 9, dphi_a$ = 10, dphi_b$ = 11, &
          eta_x0$ = 12, etap_x0$ = 13, eta_y0$ = 14, etap_y0$ = 15, &
          eta_x1$ = 16, etap_x1$ = 17, eta_y1$ = 18, etap_y1$ = 19, &
          match_end$ = 20, &
          x0$ = 21, px0$ = 22, y0$ = 23, py0$ = 24, z0$ = 25, pz0$ = 26, &
          x1$ = 27, px1$ = 28, y1$ = 29, py1$ = 30, z1$ = 31, pz1$ = 32, &
          match_end_orbit$ = 33, c_11$ = 34, c_12$ = 35, c_21$ = 36, c_22$ = 37, gamma_c$ = 39 

integer, parameter :: x$ = 1, px$ = 2, y$ = 3, py$ = 4, z$ = 5, pz$ = 6
integer, parameter :: t$ = 8
integer, parameter :: field_x$ = 10, field_y$ = 11, phase_x$ = 12, phase_y$ = 13
integer, parameter :: e_photon$ = 14

integer, parameter :: x_beam_start$ = 1, px_beam_start$ = 2, y_beam_start$ = 3
integer, parameter :: py_beam_start$ = 4, z_beam_start$ = 5, pz_beam_start$ = 6
integer, parameter :: abs_time_start$ = 8

integer, parameter :: l$=1                          ! Assumed unique. Do not assign 1 to another attribute.
integer, parameter :: tilt$=2, command$=2, roll$=2  ! Important: tilt$ = roll$
integer, parameter :: ref_tilt$ = 3, rf_frequency$=3, direction$=3
integer, parameter :: old_command$=3, kick$=3, x_gain_err$=3
integer, parameter :: rf_frequency_err$=4, k1$=4, sig_x$=4, harmon$=4, h_displace$=4, y_gain_err$=4
integer, parameter :: critical_angle_factor$ = 4, tilt_corr$ = 4, ref_coordinates$ = 4
integer, parameter :: lr_freq_spread$=5, graze_angle$=5, k2$=5, sig_y$=5, b_max$=5, v_displace$=5
integer, parameter :: flexible$ = 5, crunch$=5, ref_orbit_follows$=5
integer, parameter :: gradient$=6, k3$=6, sig_z$=6, noise$=6, new_branch$ = 6
integer, parameter :: g$=6, bragg_angle_in$ = 6
integer, parameter :: g_err$=7, n_pole$=7, bbi_const$=7, osc_amplitude$=7
integer, parameter :: gradient_err$=7, critical_angle$ = 7
integer, parameter :: bragg_angle_out$ = 7, ix_to_branch$=7
integer, parameter :: rho$=8, voltage$=8, delta_e$ = 8
integer, parameter :: charge$=8, x_gain_calib$=8, ix_to_element$=8
integer, parameter :: d1_thickness$ = 9, voltage_err$=9, rel_tracking_charge$ = 9
integer, parameter :: l_chord$=9, ks$=9, n_slice$=9, y_gain_calib$=9, bragg_angle$=9
integer, parameter :: polarity$=10, crunch_calib$=10, alpha_angle$=10, d2_thickness$ = 10
integer, parameter :: e1$=10, e_loss$=10, dks_ds$=10, gap$=10
integer, parameter :: grad_loss_sr_wake$=11, ds_path_length$=11
integer, parameter :: e2$=11, x_offset_calib$=11, v1_unitcell$=11, psi_angle$=11
integer, parameter :: y_offset_calib$=12, v_unitcell$=12, v2_unitcell$=12
integer, parameter :: traveling_wave$ = 12
integer, parameter :: fint$=12, fintx$=13, hgap$=14, hgapx$=15, h1$=16, h2$=17
integer, parameter :: phi0$=13, tilt_calib$=13, f0_re$=13, f0_re1$=13
integer, parameter :: phi0_err$=14, coef$=14, current$=14, l_pole$=14
integer, parameter :: de_eta_meas$=14, f0_im$=14, f0_im1$ = 14
integer, parameter :: quad_tilt$=14, bend_tilt$=15, x_quad$=16, y_quad$=17
integer, parameter :: dphi0$=15, n_sample$=15, fh_re$=15, f0_re2$=15, origin_ele_ref_pt$=15
integer, parameter :: dphi0_ref$ = 16, fh_im$=16, f0_im2$=16, x_half_length$=16, dx_origin$= 16
integer, parameter :: dphi0_max$=17, ref_polarization$=17, y_half_length$=17, dy_origin$ = 17
integer, parameter :: dz_origin$ = 18
integer, parameter :: fringe_type$ = 18, floor_set$ = 18, ptc_dir$ = 18
integer, parameter :: kill_fringe$ = 19, dtheta_origin$ = 19
integer, parameter :: b_param$ = 19
integer, parameter :: l_hard_edge$ = 20, dphi_origin$ = 20, ref_cap_gamma$ = 20
integer, parameter :: field_scale$ = 21, dpsi_origin$ = 21, darwin_width_sigma$ = 21
integer, parameter :: angle$=22, n_cell$=22, x_ray_line_len$=22, darwin_width_pi$ = 22

integer, parameter :: x_pitch$ = 23
integer, parameter :: y_pitch$ = 24  
integer, parameter :: x_offset$ = 25
integer, parameter :: y_offset$ = 26 
integer, parameter :: z_offset$ = 27 ! Assumed unique. Do not overload further.
integer, parameter :: hkick$ = 28, d_spacing$ = 28, t_offset$ = 28
integer, parameter :: vkick$ = 29, l_x$ = 29
integer, parameter :: BL_hkick$ = 30, l_y$ = 30        ! l_y$ = l_x$ + 1
integer, parameter :: BL_vkick$ = 31, l_z$ = 31        ! l_z$ = l_x$ + 2
integer, parameter :: BL_kick$ = 32, coupler_at$ = 32
integer, parameter :: B_field$ = 33, E_field$ = 33, coupler_phase$ = 33
integer, parameter :: coupler_angle$ = 34, B_field_err$ = 34
integer, parameter :: coupler_strength$ = 35
integer, parameter :: B1_gradient$ = 35, E1_gradient$ = 35
integer, parameter :: B2_gradient$ = 36, E2_gradient$ = 36, h_x_norm$ = 36
integer, parameter :: B3_gradient$ = 37, E3_gradient$ = 37, h_y_norm$ = 37, ptc_field_geometry$ = 38
integer, parameter :: Bs_field$ = 38, e_tot_offset$ = 38, h_z_norm$ = 38
integer, parameter :: delta_ref_time$ = 39 ! Assumed unique Do not overload.
integer, parameter :: p0c_start$ = 40
integer, parameter :: e_tot_start$ = 41   
integer, parameter :: p0c$ = 42         ! Assumed unique. Do not overload.
integer, parameter :: e_tot$ = 43       ! Assumed unique. Do not overload.
integer, parameter :: x_pitch_tot$ = 44, no_end_marker$ = 44
integer, parameter :: y_pitch_tot$ = 45
integer, parameter :: x_offset_tot$ = 46
integer, parameter :: y_offset_tot$ = 47
integer, parameter :: z_offset_tot$ = 48
integer, parameter :: tilt_tot$ = 49, roll_tot$ = 49  ! Important: tilt_tot$ = roll_tot$
integer, parameter :: pole_radius$ = 50, ref_tilt_tot$ = 50
integer, parameter :: n_ref_pass$ = 51
integer, parameter :: radius$ = 52
integer, parameter :: ref_time_start$ = 53
integer, parameter :: thickness$ = 54, integrator_order$ = 54   ! For Etiennes' PTC: 2, 4, or 6.
integer, parameter :: num_steps$ = 55
integer, parameter :: ds_step$ = 56
integer, parameter :: lord_pad1$ = 57
integer, parameter :: lord_pad2$ = 58, ref_wavelength$ = 58
integer, parameter :: scratch$ = 59
integer, parameter :: custom_attribute1$ = 61   ! For general use
integer, parameter :: custom_attribute2$ = 62   ! For general use
integer, parameter :: custom_attribute3$ = 63   ! For general use
integer, parameter :: custom_attribute4$ = 64   ! For general use
integer, parameter :: custom_attribute5$ = 65, custom_attribute_max$ = 65   ! For general use
integer, parameter :: x1_limit$ = 66   ! Assumed unique. Do not overload.
integer, parameter :: x2_limit$ = 67   ! Assumed unique. Do not overload.
integer, parameter :: y1_limit$ = 68   ! Assumed unique. Do not overload.
integer, parameter :: y2_limit$ = 69   ! Assumed unique. Do not overload.
integer, parameter :: check_sum$ = 70  ! Assumed unique. Do not overload.

!! 71 = 1 + num_ele_attrib$

integer, parameter :: lr_wake_file$ = 71, alpha_b$ = 71, use_hard_edge_drifts$ = 71
integer, parameter :: alias$ =72, eta_x$ = 72, ptc_max_fringe_order$ = 72
integer, parameter :: start_edge$ =73, eta_y$ = 73
integer, parameter :: end_edge$ =74, etap_x$ = 74
integer, parameter :: accordion_edge$ =75, etap_y$ = 75
integer, parameter :: lattice$ = 76, phi_a$ = 76, diffraction_type$ = 76
integer, parameter :: aperture_type$ = 77, eta_z$ = 77
integer, parameter :: map_with_offsets$ = 78, cmat_11$ = 78, surface_attrib$ = 78
integer, parameter :: csr_calc_on$ = 79, cmat_12$ = 79
integer, parameter :: s_position$ = 80, cmat_21$ = 80
integer, parameter :: mat6_calc_method$ = 81, cmat_22$ = 81
integer, parameter :: tracking_method$  = 82, s_long$ = 82
integer, parameter :: ref_time$ = 83, ptc_integration_type$ = 83
integer, parameter :: spin_tracking_method$ = 84, eta_a$ = 84
integer, parameter :: aperture$ = 85, rf_auto_scale_amp$ = 85, etap_a$ = 85
integer, parameter :: x_limit$ = 86, absolute_time_tracking$ = 86, eta_b$ = 86
integer, parameter :: y_limit$ = 87, rf_auto_scale_phase$ = 87, etap_b$ = 87
integer, parameter :: offset_moves_aperture$ = 88
integer, parameter :: aperture_limit_on$ = 89

integer, parameter :: ptc_exact_misalign$ = 90
integer, parameter :: sr_wake_file$ = 90, alpha_a$ = 90
integer, parameter :: term$ = 91, use_ptc_layout$ = 91
integer, parameter :: x_position$ = 92, s_spline$ = 92, ptc_exact_model$ = 92
integer, parameter :: symplectify$ = 93, y_position$ = 93, n_slice_spline$ = 93
integer, parameter :: z_position$ = 94
integer, parameter :: is_on$ = 95, theta_position$ = 95
integer, parameter :: field_calc$ = 96, phi_position$ = 96
integer, parameter :: psi_position$ = 97
integer, parameter :: aperture_at$ = 98, beta_a$ = 98
integer, parameter :: ran_seed$ = 99, beta_b$ = 99, origin_ele$= 99

integer, parameter :: to_line$ = 100
integer, parameter :: field_master$ = 101, harmon_master$ = 101, to_element$ = 101
integer, parameter :: descrip$ = 102
integer, parameter :: scale_multipoles$ = 103
integer, parameter :: wall_attribute$ = 104  ! Do not confuse this with wall3d$
integer, parameter :: field$ = 105
integer, parameter :: phi_b$ = 106, crystal_type$ = 106
integer, parameter :: type$ = 107
integer, parameter :: ref_origin$ = 108
integer, parameter :: ele_origin$ = 109

! superimpose$ through create_jumbo_slave$ assumed unique (or need to modify bmad_parser_mod.f90).

integer, parameter :: superimpose$    = 110   
integer, parameter :: offset$         = 111
integer, parameter :: reference$      = 112
integer, parameter :: ele_beginning$  = 113
integer, parameter :: ele_center$     = 114
integer, parameter :: ele_end$        = 115
integer, parameter :: ref_beginning$  = 116
integer, parameter :: ref_center$     = 117
integer, parameter :: ref_end$        = 118
integer, parameter :: create_jumbo_slave$ = 119

integer, parameter :: a0$  = 120, k0l$  = 120
integer, parameter :: a20$ = 140, k20l$ = 140

integer, parameter :: b0$  = 150, t0$  = 150
integer, parameter :: b20$ = 170, t20$ = 170 

integer, parameter :: num_ele_attrib_extended$ = t20$

character(40), parameter :: null_name$ = '!NULL' 
character(40), parameter :: blank_name$ = ' '

! lattice logical names

integer, parameter :: open$ = 1, closed$ = 2

character(16), parameter :: lattice_type_name(0:2) = ['GARBAGE!        ', 'LINEAR_LATTICE  ', 'CIRCULAR_LATTICE']
character(16), parameter :: geometry_name(0:2) = ['GARBAGE!    ', 'OPEN        ', 'CLOSED      ']

! logicals for MAKE_HYBIRD_lat

logical, parameter :: remove_markers$ = .true., no_remove_markers$ = .false.

! control element logicals
! Idea: Combine girder_lord, overlay_lord and group_lord -> control_lord

integer, parameter :: free$ = 1, super_slave$ = 2, control_slave$ = 3
integer, parameter :: group_lord$ = 4, super_lord$ = 5, overlay_lord$ = 6
integer, parameter :: girder_lord$ = 7, multipass_lord$ = 8, multipass_slave$ = 9
integer, parameter :: not_a_lord$ = 10, slice_slave$ = 11, control_lord$ = 12

character(16), parameter :: control_name(12) = [ &
            'FREE           ', 'SUPER_SLAVE    ', 'CONTROL_SLAVE  ', 'GROUP_LORD     ', &
            'SUPER_LORD     ', 'OVERLAY_LORD   ', 'GIRDER_LORD    ', 'MULTIPASS_LORD ', &
            'MULTIPASS_SLAVE', 'NOT_A_LORD     ', 'SLICE_SLAVE    ', 'CONTROL_LORD   ']

logical, parameter :: set$ = .true., unset$ = .false.

! Note: custom$ = 7, and taylor$ = 8 are taken from the element key list.

integer, parameter :: bmad_standard$ = 1, symp_lie_ptc$ = 2
integer, parameter :: runge_kutta$ = 3 
integer, parameter :: linear$ = 4, tracking$ = 5, symp_map$ = 6
integer, parameter :: hard_edge_model$ = 9, symp_lie_bmad$ = 10, static$ = 11
integer, parameter :: boris$ = 12, mad$ = 14
integer, parameter :: time_runge_kutta$ = 15, custom2$ = 16
integer, parameter :: n_methods$ = 16

character(16), parameter :: tracking_method_name(0:n_methods$) = [ &
      'GARBAGE!        ', 'Bmad_Standard   ', 'Symp_Lie_PTC    ', 'Runge_Kutta     ', &
      'Linear          ', 'Garbage         ', 'Symp_Map        ', 'Custom          ', &
      'Taylor          ', 'Garbage         ', 'Symp_Lie_Bmad   ', 'Static          ', &
      'Boris           ', 'GARBAGE!        ', 'MAD             ', 'Time_Runge_Kutta', &
      'Custom2         ']

character(16), parameter :: spin_tracking_method_name(0:n_methods$) = [ &
      'GARBAGE!        ', 'Bmad_Standard   ', 'Symp_Lie_PTC    ', 'Garbage         ', &
      'Garbage         ', 'Tracking        ', 'Garbage         ', 'Custom          ', &
      'Garbage         ', 'Garbage         ', 'Symp_Lie_Bmad   ', 'Garbage         ', &
      'Garbage         ', 'GARBAGE!        ', 'Garbage         ', 'Garbage         ', &
      'Custom2         ']

character(16), parameter :: mat6_calc_method_name(0:n_methods$) = [ &
      'GARBAGE!        ', 'Bmad_Standard   ', 'Symp_Lie_PTC    ', 'Garbage         ', &
      'Linear          ', 'Tracking        ', 'Symp_Map        ', 'Custom          ', &
      'Taylor          ', 'Garbage         ', 'Symp_Lie_Bmad   ', 'Static          ', &
      'Garbage         ', 'GARBAGE!        ', 'MAD             ', 'Garbage        a', &
      'Custom2         ']

integer, parameter :: drift_kick$ = 1, matrix_kick$ = 2, ripken_kick$ = 3
character(16), parameter :: ptc_integration_type_name(0:3) = [&
         'GARBAGE!   ', 'Drift_Kick ', 'Matrix_Kick', 'Ripken_Kick']


! sbend$ and rbend$ are from key definitions.

integer, parameter :: map_type$ = 1, periodic_type$ = 3, const_ref_energy$ = 4, nonconst_ref_energy$ = 5
character(20), parameter :: sub_key_name(0:18) = ['GARBAGE!           ', 'Map                ', &
    'SBend              ', 'Periodic           ', 'Const_Ref_Energy   ', 'NonConst_Ref_Energy', &
    'GARBAGE!           ', 'GARBAGE!           ', 'GARBAGE!           ', 'GARBAGE!           ', &
    'GARBAGE!           ', 'GARBAGE!           ', 'GARBAGE!           ', 'GARBAGE!           ', &
    'GARBAGE!           ', 'GARBAGE!           ', 'GARBAGE!           ', 'GARBAGE!           ', &
    'RBend              ']

! field_calc names.
! Note: refer_to_lords is an "internal" value which is not valid for use in a lattice file.
!   The period in "Refer_to_Lords." is used to prevent sets in the lattice file.

integer, parameter :: grid$ = 2, map$ = 3, Refer_to_lords$ = 4
character(16), parameter :: field_calc_name(0:7) = &
    ['GARBAGE!       ', 'Bmad_Standard  ', 'Grid           ', 'Map            ', &
     'Refer_to_Lords.', 'GARBAGE!       ', 'GARBAGE!       ', 'Custom         ']

! Crystal sub_key values.

integer, parameter :: bragg$ = 1, laue$ = 2
character(8), parameter :: diffraction_type_name(0:2) = ['GARBAGE!', 'Bragg   ', 'Laue    ']

! The linac_normal_mode_struct is basically the synchrotron integrals with the
! energy factors thrown in. Useful for linacs.

type anormal_mode_struct
  real(rp) emittance        ! Beam emittance
  real(rp) synch_int(4:6)   ! Synchrotron integrals
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
  real(rp) rf_voltage     ! Total rfcavity voltage (eV)
  real(rp) pz_aperture    ! pz aperture limit
  type (anormal_mode_struct)  a, b, z
  type (linac_normal_mode_struct) lin
end type

integer, parameter :: bends$ = 201
integer, parameter :: wigglers$ = 202
integer, parameter :: all$ = 203

!---------------------------------------------------------------------------
! Units

integer, parameter :: radians$ = 1, degrees$ = 2, cycles$ = 3, kHz$ = 4
character(8), parameter :: frequency_units_name(4) = ['Radians ', 'Degrees ', 'Cycles  ', 'kHz     ']

! Electric and magnetic fields.

type em_field_struct
  real(rp) E(3)         ! electric field
  real(rp) B(3)         ! magnetic field
  real(rp) dE(3,3)      ! electric field gradient
  real(rp) dB(3,3)      ! magnetic field gradient
end type


! Grid of em_grid information
integer, parameter :: rotationally_symmetric_rz$ = 1
character(30), parameter :: em_grid_type_name(1) = [ 'rotationally_symmetric_rz     ' ]
integer, parameter :: em_grid_dimension(1) = [ 2 ] 


! Structures for saving the track through an element.
! track%pt(0:n)  goes from 0 to n = track%n_pt

type track_map_struct
  real(rp) vec0(6)        ! 0th order part of xfer map from the beginning.
  real(rp) mat6(6,6)      ! 1st order part of xfer map (transfer matrix).
end type

type track_struct
  type (coord_struct), allocatable :: orb(:)      ! An array of track points: %orb(0:) 
  type (em_field_struct), allocatable:: field(:)  ! An array of em fields: %field(0:) 
  type (track_map_struct), allocatable :: map(:)  ! An array of maps: %map(0:)
  real(rp) :: ds_save = 1e-3                      ! Min distance between points.
  integer :: n_pt                                 ! Actual track upper bound for %orb(0:) and  %map(0:)
  integer :: n_bad                                ! Number of bad steps when adaptive tracking is done.
  integer :: n_ok                                 ! Number of good steps when adaptive tracking is done.
end type

!------------------------------------------------------------------------------
! misc

! This is for debugging radiation damping and fluctuations.

type synch_rad_common_struct
  real(rp) :: scale = 1.0               ! used to scale the radiation
  real(rp) :: i2 = 0, i3 = 0            ! radiation integrals
  real(rp) :: i5a = 0, i5b = 0
  logical :: i_calc_on = .false.        ! For calculating i2 and i3    
end type

type (synch_rad_common_struct), save :: synch_rad_com

integer, parameter :: is_logical$ = 1, is_integer$ = 2, is_real$ = 3, is_switch$ = 4, is_string$ = 5

! Note: custom$ = 7 is taken from element key names.

integer, parameter :: rectangular$ = 1, elliptical$ = 2, wall3d$ = 3
character(16), parameter :: aperture_type_name(0:7) = &
                                    ['garbage!   ', 'Rectangular', 'Elliptical ', 'Wall3D     ', &
                                     'garbage!   ', 'garbage!   ', 'garbage!   ', 'Custom     ']

integer, parameter :: sigma_polarization$ = 1, pi_polarization$ = 2
character(20) :: polarization_name(0:2) = ['Garbage!          ', 'Sigma_polarization', 'Pi_polarization   ']

! fringe_type$ 

integer, parameter :: full_straight$ = 1, full_bend$ = 2, none$ = 3, basic_bend$ = 4
character(16), parameter :: fringe_type_name(0:4) = ['Garbage!      ', &
                'FULL_STRAIGHT ', 'FULL_BEND     ',  'NONE          ', 'BASIC_BEND    ']

! extra_parsing_info_struct is used by parsing routines.
! %deterministic:
!   0 = Not, 1 = Ran state on input deterministic, 2 = ran state deterministice
!   will be generated in exactly the same way every time?
! %ran_function_was_called:
!   Only set True when ran function is called with ran_determinisitc = 0.
! %determinisitc_ran_function_was_called:
!   Only set True when ran function is called with ran_determinisitc = 1.

type extra_parsing_info_struct
  type (random_state_struct) :: initial_state
  integer deterministic   
  logical ran_function_was_called 
  logical deterministic_ran_function_was_called 
  logical taylor_order_set, ptc_max_fringe_order_set, use_hard_edge_drifts_set
  integer taylor_order, ptc_max_fringe_order
  logical use_hard_edge_drifts
end type

! ptc_field_geometry for bends

integer, parameter :: sector$ = 1, straight$ = 2, true_rbend$ = 3
character(16), parameter :: ptc_field_geometry_name(0:3) = [ &
                              'Garbage!  ', 'Sector    ', 'straight  ', 'true_rbend']

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
  real(rp) :: d_orb(6)           = 1e-5      ! Orbit deltas for the mat6 via tracking calc.
  real(rp) :: default_ds_step    = 0.2_rp    ! Integration step size.  
  real(rp) :: significant_length = 1e-10     ! meter 
  real(rp) :: rel_tol_tracking = 1e-8
  real(rp) :: abs_tol_tracking = 1e-10
  real(rp) :: rel_tol_adaptive_tracking = 1e-8     ! Adaptive tracking relative tolerance.
  real(rp) :: abs_tol_adaptive_tracking = 1e-10    ! Adaptive tracking absolute tolerance.
  real(rp) :: init_ds_adaptive_tracking = 1e-3     ! Initial step size
  real(rp) :: min_ds_adaptive_tracking = 1e-8      ! Min step size.
  integer :: taylor_order = 3                      ! 3rd order is default
  integer :: default_integ_order = 2               ! PTC integration order. 
  integer :: ptc_max_fringe_order = 2              ! PTC max fringe order (2 => Quadrupole !). 
                                                   !   Must call set_ptc after changing.
  logical :: use_hard_edge_drifts = .true.         ! Insert drifts when tracking through cavity?
  logical :: sr_wakes_on = .true.                  ! Short range wakefields?
  logical :: lr_wakes_on = .true.                  ! Long range wakefields
  logical :: mat6_track_symmetric = .true.         ! symmetric offsets
  logical :: auto_bookkeeper = .true.              ! Automatic bookkeeping?
  logical :: space_charge_on = .false.             ! Space charge switch
  logical :: coherent_synch_rad_on = .false.       ! csr 
  logical :: spin_tracking_on = .false.            ! spin tracking?
  logical :: radiation_damping_on = .false.        ! Damping toggle.
  logical :: radiation_fluctuations_on = .false.   ! Fluctuations toggle.
  logical :: conserve_taylor_maps = .true.         ! Enable bookkeeper to set ele%map_with_offsets = F?
  logical :: photon_tracking_uses_field = .false.  ! False => Track intensity (photons are incoherent).
  logical :: absolute_time_tracking_default = .false.   ! Default for lat%absolute_time_tracking
  logical :: rf_auto_scale_phase_default = .true.       ! Default for lat%rf_auto_scale_phase
  logical :: rf_auto_scale_amp_default = .true.         ! Default for lat%rf_auto_scale_amp
  logical :: use_ptc_layout_default = .false.           ! Default for lat%use_ptc_layout
  logical :: debug = .false.                            ! Used for code debugging.
end type
  
type (bmad_common_struct), save :: bmad_com

! Some routines need to keep track of where elements are when elements are added or removed from
! the lattice. 

type ele_loc_ele_struct
  type (lat_ele_loc_struct), allocatable :: ele(:)
end type

type ele_loc_branch_struct
  type (ele_loc_ele_struct), allocatable :: branch(:)
end type

type (ele_loc_branch_struct), save, target :: ele_loc_com

! This structure stores the radiation integrals for an individual element except
! lin_norm_emit_a and lin_norm_emit_b are running sums from the beginning of the branch.
! See also rad_int_all_ele_struct.

type rad_int1_struct
  real(rp) :: i0   = 0
  real(rp) :: i1   = 0
  real(rp) :: i2   = 0
  real(rp) :: i3   = 0
  real(rp) :: i4a  = 0
  real(rp) :: i4b  = 0
  real(rp) :: i4z  = 0
  real(rp) :: i5a  = 0
  real(rp) :: i5b  = 0
  real(rp) :: i6b  = 0
  real(rp) :: lin_i2_E4   = 0 
  real(rp) :: lin_i3_E7   = 0
  real(rp) :: lin_i5a_E6  = 0
  real(rp) :: lin_i5b_E6  = 0
  real(rp) :: lin_norm_emit_a = 0 ! Running sum
  real(rp) :: lin_norm_emit_b = 0 ! Running sum
  real(rp) :: n_steps         = 0 ! number of qromb steps needed
end type

! Structure for an array of rad_int1_structs for a single branch

type rad_int_all_ele_struct
  type (rad_int1_struct), allocatable :: ele(:) ! Array is indexed from 0
end type

contains

!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
!+
! Function coord_state_name (coord_state) result (state_str)
!
! Routine to return the string representation of a coord%state state.
!
! Input:
!   coord_state -- integer: coord%state value
!
! Output:
!   state_str   -- character(20): String representation.
!-

function coord_state_name (coord_state) result (state_str)

implicit none

integer coord_state
character(20) state_str

!

select case (coord_state)
case (alive$);                 state_str = 'Alive'
case (lost$);                  state_str = 'Lost'
case (not_set$);               state_str = 'Not_Set'
Case (lost_neg_x_aperture$);   state_str = 'Lost_Neg_X_Aperture'
case (lost_pos_x_aperture$);   state_str = 'Lost_Pos_X_Aperture'
case (lost_neg_y_aperture$);   state_str = 'Lost_Neg_Y_Aperture'
case (lost_pos_y_aperture$);   state_str = 'Lost_Pos_Y_Aperture'
case (lost_z_aperture$);       state_str = 'Lost_Z_Aperture'
case default;                  state_str = 'UNKNOWN!'
end select

end function coord_state_name

end module
