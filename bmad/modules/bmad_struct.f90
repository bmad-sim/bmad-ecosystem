!+
! Bmad_struct holds the structure definitions for Bmad routines.
!-

module bmad_struct

use random_mod
use spline_mod
use sim_utils
use cubic_interpolation_mod

use ptc_spin, only: genfield, fibre, layout, c_damap, c_normal_form, c_taylor, probe_8, internal_state, c_quaternion

private next_in_branch

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! IF YOU CHANGE THE LAT_STRUCT OR ANY ASSOCIATED STRUCTURES YOU MUST INCREASE THE VERSION NUMBER !!!
! THIS IS USED BY BMAD_PARSER TO MAKE SURE DIGESTED FILES ARE OK.

integer, parameter :: bmad_inc_version$ = 334

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

integer, parameter :: none$ = 1

!   hard_ele     -- ele_struct, pointer: Points to element with the hard edge.
!                     Will be nullified if there is no hard edge.
!                     This will be track_ele (if it has a hard edge) unless track_ele is a super_slave.
!   s_edge_hard  -- Real(rp): S-position of next hard egde in hard_ele frame.
!   particle_at  -- Edge of hard_ele being tracked through. 

type fringe_field_info_struct
  type (ele_struct), pointer :: hard_ele => null()
  real(rp) :: s_edge_hard = 0
  real(rp) :: ds_edge = 0                     ! Distance from particle to edge in hard_ele frame.
  integer :: particle_at = -1                 ! first_track_edge$, second_track_edge$, or none$
  integer, pointer :: hard_location => null() ! Particle location wrt hard_ele. Points to element in location(:).
  integer, allocatable :: location(:)         ! Particle location in an element. entrance_end$, inside$, or exit_end$
                                              ! Elements in list are the tracking element or its lords.
  logical :: has_fringe = .false.             ! Has a fringe to worry about?
end type

!

integer, parameter :: n_pole_maxx = 21  ! maximum multipole order

integer, parameter :: old_control_var_offset$ = 1000  ! For indexing into ele%control%var(:) array
integer, parameter :: var_offset$ = 2000              ! Important: var_offset$ > old_control_var_offset$
integer, parameter :: n_var_max$ = 999                ! Maximum number of variables per controller.
integer, parameter :: taylor_offset$ = 1000000000     ! Taylor term index offset.

type expression_atom_struct
  character(40) :: name = ''
  integer :: type = 0   ! plus$, minum$, sin$, cos$, etc. To convert to string use: expression_op_name
  real(rp) :: value = 0
end type

!-------------------------------------------------------------------------
! Note: custom$ = 7, and taylor$ = 8 are taken from the element key list.

integer, parameter :: bmad_standard$ = 1, symp_lie_ptc$ = 2, runge_kutta$ = 3 
integer, parameter :: linear$ = 4, tracking$ = 5, time_runge_kutta$ = 6
integer, parameter :: fixed_step_runge_kutta$ = 9, symp_lie_bmad$ = 10, magnus$ = 11
integer, parameter :: Auto$ = 12, sprint$ = 12, fixed_step_time_runge_kutta$ = 13, mad$ = 14
integer, parameter :: transverse_kick$ = 3, spin_integration$ = 99

character(28), parameter :: tracking_method_name(0:14) = [character(28) :: &
      'GARBAGE!', 'Bmad_Standard',               'Symp_Lie_PTC',     'Runge_Kutta', &
      'Linear',   'GARBAGE!',                    'Time_Runge_Kutta', 'Custom', &
      'Taylor',   'Fixed_Step_Runge_Kutta',      'Symp_Lie_Bmad',    'GARBAGE!', &
      'GARBAGE!', 'Fixed_Step_Time_Runge_kutta', 'MAD']

character(16), parameter :: spin_tracking_method_name(0:12) = [ &
      'GARBAGE!        ', 'Off             ', 'Symp_Lie_PTC    ', 'Transverse_Kick ', &
      'GARBAGE!        ', 'Tracking        ', 'GARBAGE!        ', 'Custom          ', &
      'GARBAGE!        ', 'GARBAGE!        ', 'GARBAGE!        ', 'Magnus          ', &
      'Sprint          ']

character(24), parameter :: mat6_calc_method_name(0:14) = [character(24):: 'GARBAGE!', &
      'Bmad_Standard', 'Symp_Lie_PTC', 'GARBAGE!',  'Linear', 'Tracking', &
      'GARBAGE!', 'Custom', 'Taylor', 'GARBAGE!', 'Symp_Lie_Bmad', &
      'GARBAGE!', 'Auto', 'GARBAGE!', 'MAD']

integer, parameter :: drift_kick$ = 1, matrix_kick$ = 2, ripken_kick$ = 3
character(16), parameter :: ptc_integration_type_name(0:3) = [&
         'GARBAGE!   ', 'Drift_Kick ', 'Matrix_Kick', 'Ripken_Kick']

! sbend$ and rbend$ are from key definitions.

character(16), parameter :: sub_key_name(0:18) = ['GARBAGE!     ', 'GARBAGE!     ', &
    'SBend        ', 'GARBAGE!     ', 'GARBAGE!     ', 'GARBAGE!     ', &
    'GARBAGE!     ', 'GARBAGE!     ', 'GARBAGE!     ', 'GARBAGE!     ', &
    'GARBAGE!     ', 'GARBAGE!     ', 'GARBAGE!     ', 'GARBAGE!     ', &
    'GARBAGE!     ', 'GARBAGE!     ', 'GARBAGE!     ', 'GARBAGE!     ', &
    'RBend        ']

! ptc_field_geometry for bends

integer, parameter :: sector$ = 1, straight$ = 2
character(16), parameter :: ptc_field_geometry_name(0:2) = [ &
                              'Garbage!  ', 'Sector    ', 'Straight  ']

! field_calc names.
! Note: refer_to_lords is an "internal" value which is not valid for use in a lattice file.
!   The period in "Refer_to_Lords." is used to prevent sets in the lattice file.

integer, parameter :: fieldmap$ = 2, planar_model$ = 3, Refer_to_lords$ = 4, no_field$ = 5
integer, parameter :: helical_model$ = 6, soft_edge$ = 8
character(16), parameter :: field_calc_name(0:8) = &
    [character(16):: 'GARBAGE!', 'Bmad_Standard', 'FieldMap', 'Planar_Model', &
     'Refer_to_Lords.', 'No_Field', 'Helical_Model', 'Custom', 'Soft_edge']

! Distribution

integer, parameter :: uniform$ = 1, gaussian$ = 2, spherical$ = 3, curve$ = 4
character(12), parameter :: distribution_name(0:4) = [character(12):: 'GARBAGE!', &
                                                       'Uniform', 'Gaussian', 'Spherical', 'Curve']

! Control element logicals.
! Note: super_slave$ and multipass_slave$ are also used as possible settings of the 
! why_not_free argument in attribute_free(...).

integer, parameter :: ix_slice_slave$ = -2 ! Index to set slice_slave%ix_ele to.

integer, parameter :: minor_slave$ = 1, super_slave$ = 2, free$ = 3
integer, parameter :: group_lord$ = 4, super_lord$ = 5, overlay_lord$ = 6
integer, parameter :: girder_lord$ = 7, multipass_lord$ = 8, multipass_slave$ = 9
integer, parameter :: not_a_lord$ = 10, slice_slave$ = 11, control_lord$ = 12, ramper_lord$ = 13
integer, parameter :: governor$ = 14, field_lord$ = 15    ! governor$ = Union of overlay and group lords.

character(20), parameter :: control_name(13) = [character(20):: &
            'Minor_Slave', 'Super_Slave', 'Free', 'Group_Lord', &
            'Super_Lord', 'Overlay_Lord', 'Girder_Lord', 'Multipass_Lord ', &
            'Multipass_Slave', 'Not_a_Lord', 'Slice_Slave', 'Control_Lord', 'Ramper_Lord']

logical, parameter :: set$ = .true., unset$ = .false.

! Auto aperture is different from Rectangular aperture in that with Auto, aperture values are calculated by Bmad 
! when the lattice file is parsed. Auto is only available for Mask, Detector and Diffraction_plate elements.

integer, parameter :: auto_aperture$ = 1, rectangular$ = 2, elliptical$ = 3, wall3d$ = 5, custom_aperture$ = 7
integer, parameter :: lord_defined$ = 8

character(16), parameter :: aperture_type_name(0:8) = [character(16):: &
                               'garbage!   ', 'Auto       ', 'Rectangular', 'Elliptical ', &
                               'Surface    ', 'Wall3D     ', 'garbage!   ', 'Custom     ', 'Lord_Defined']

! fringe_type
! non-bend fringe type names are in the range fringe_type(1:n_non_bend_fringe_type$)

integer, parameter :: soft_edge_only$ = 2, hard_edge_only$ = 3, full$ = 4
integer, parameter :: sad_full$ = 5, linear_edge$ = 6, basic_bend$ = 7

character(16), parameter :: fringe_type_name(0:7) = [character(16):: 'Garbage!', &
                                   'None', 'Soft_Edge_Only', 'Hard_edge_only', 'Full', &
                                   'SAD_Full', 'Linear_Edge', 'Basic_Bend']

integer, parameter :: standing_wave$ = 1, traveling_wave$ = 2, ptc_standard$ = 3
character(16), parameter :: cavity_type_name(0:3) = ['Garbage!      ', 'Standing_Wave ', 'Traveling_Wave', 'PTC_Standard  ']


integer, parameter :: x_invariant$ = 1, multipole_symmetry$ = 2
character(16), parameter :: ptc_fringe_geometry_name(0:2) = ['Garbage!          ', 'x_invariant       ', 'multipole_symmetry']

integer, parameter :: control_var$ = 1, old_control_var$ = 2, all_control_var$ = 3, elec_multipole$ = 4

!

integer, parameter :: ok$ = 1, in_stop_band$ = 2, non_symplectic$ = 3, unstable$ = 4
integer, parameter :: unstable_a$ = 5, unstable_b$ = 6
integer, parameter :: xfer_mat_calc_failure$ = 7, twiss_propagate_failure$ = 8, no_closed_orbit$ = 9

character(24) :: matrix_status_name(9) = [character(24) :: 'OK', 'IN_STOP_BAND', 'NON_SYMPLECTIC', &
                       'UNSTABLE', 'UNSTABLE A-MODE', 'UNSTABLE B-MODE', 'XFER_MAT_CALC_FAILURE', &
                       'TWISS_PROPAGATE_FAILURE', 'NO_CLOSED_ORBIT']

type twiss_struct
  real(rp) :: beta = 0, alpha = 0, gamma = 0, phi = 0, eta = 0, etap = 0, deta_ds = 0
  real(rp) :: sigma = 0, sigma_p = 0, emit = 0, norm_emit = 0
  real(rp) :: dbeta_dpz = 0, dalpha_dpz = 0, deta_dpz = 0, detap_dpz = 0
end type

! Misc parameters

integer, parameter :: include_kicks$ = 1, short$ = 8

integer, parameter :: user_set$ = 0, first_pass$ = 1
character(12), parameter :: multipass_ref_energy_name(0:1) = [character(12):: 'User_Set', 'First_Pass']

integer, parameter :: highland$ = 2, lynch_dahl$ = 3
character(12), parameter :: scatter_method_name(3) = [character(12):: 'Off', 'Highland', 'Lynch_Dahl']

!-------------------------------------------------------------------------
! Structure for holding the photon reflection probability tables.
! Used for both smooth surface reflections and custom crystal reflections.

type interval1_coef_struct
  real(rp) c0, c1, n_exp
end type

type photon_reflect_table_struct
  real(rp), allocatable :: angle(:)               ! Vector of angle values for %p_reflect
  real(rp), allocatable :: energy(:)              ! Vector of energy values for %p_reflect
  type (interval1_coef_struct), allocatable :: int1(:)
  real(rp), allocatable :: p_reflect(:,:)         ! (angle, ev) probability. Log used for smooth surface reflection
  real(rp) :: max_energy = -1                     ! maximum energy for this table
  real(rp), allocatable :: p_reflect_scratch(:)   ! Scratch space
  real(rp), allocatable :: bragg_angle(:)         ! Bragg angle at energy values.
end type

! Each photon_reflect_reflect_table_array(:) represents a different surface type.
! photon_reflect_reflect_table_array(1) is initialized by the photon_reflection_init routine
! All others can be set by an outside programmer. 

type photon_reflect_surface_struct
  character(40) :: name = ''
  character(80) :: description = ''                    ! Descriptive name
  character(200) :: reflectivity_file = ''
  type (photon_reflect_table_struct), allocatable :: table(:)
  real(rp) :: surface_roughness_rms = 0       ! sigma in Dugan's notation
  real(rp) :: roughness_correlation_len = 0   ! T in Dugan's notation
  integer :: ix_surface = -1
end type

integer, parameter :: incoherent$ = 1, coherent$ = 2
character(16), parameter :: photon_type_name(1:2) = ['Incoherent', 'Coherent  ']

!-------------------------------------------------------------------------

integer, parameter :: ascii$ = 1, binary$ = 2, hdf5$ = 3, one_file$ = 4
integer, parameter :: old_ascii$ = 44    ! For testing purposes.

! num_ele_attrib$ is size of ele%value(:) array.

integer, parameter :: num_ele_attrib$ = 75
integer, parameter :: off$ = 1, on$ = 2

integer, parameter :: save_state$ = 3, restore_state$ = 4, off_and_save$ = 5

integer, parameter :: horizontally_pure$ = 2, vertically_pure$ = 3
character(20), parameter :: exact_multipoles_name(3) = [character(20):: 'Off', 'Horizontally_Pure', 'Vertically_Pure']

integer, parameter :: one_dim$ = 2, steady_state_3d$ = 3
character(20), parameter :: csr_method_name(3) = [character(20):: 'Off', '1_Dim', 'Steady_State_3D']

integer, parameter :: slice$ = 2, fft_3D$ = 3, cathode_fft_3d$ = 4
character(20), parameter :: space_charge_method_name(4) = [character(20):: 'Off', 'Slice', 'FFT_3D', 'Cathode_FFT_3D']

! Pauli matrices

complex(rp), parameter :: pauli_1(2,2) = reshape([(1,0), (0,0), (0,0), (1,0)], [2,2])
complex(rp), parameter :: pauli_x(2,2) = reshape([(0,0), (1,0), (1,0), (0,0)], [2,2])
complex(rp), parameter :: pauli_y(2,2) = reshape([(0,0), (0,1), (0,-1), (0,0)], [2,2])
complex(rp), parameter :: pauli_z(2,2) = reshape([(1,0), (0,0), (0,0), (-1,0)], [2,2])

type pauli_struct
  complex(rp) sigma(2,2)
end type

type (pauli_struct), parameter :: pauli(0:3) = [pauli_struct(pauli_1), pauli_struct(pauli_x), &
                                                pauli_struct(pauli_y), pauli_struct(pauli_z)]

! Structure for spin matching calculations.
! Naming follows Barber & Ripkin section 2.78 in the Handbook of Accelerator Physics and Engineering.

type spin_eigen_struct
  complex(rp) :: vec(8) = 0
  complex(rp) :: val = 0
end type

type spin_axis_struct
  real(rp) :: l(3) = 0         ! Transverse axis.
  real(rp) :: n0(3) = 0        ! Invariant spin axis on closed orbit.
  real(rp) :: m(3) = 0         ! Transverse axis.
end type

type spin_matching_struct
  type (spin_axis_struct) :: axis = spin_axis_struct()
  type (spin_eigen_struct) :: eigen(8) = spin_eigen_struct()
  real(rp) :: dn_dpz(3) = 0         ! Invariant spin derivative
  real(rp) :: alpha(6) = 0          ! Alpha vector
  real(rp) :: beta(6) = 0           ! Beta vector
  real(rp) :: orb0(6) = 0           ! Closed orbit 
  real(rp) :: M_1turn(8,8) = 0      ! 1-turn matrix
  real(rp) :: M_ele(8,8) = 0        ! Transfer matrix through element.
  real(rp) :: sq_ele(0:3) = 0, sq_1turn(0:3) = 0
  logical :: valid = .false.
end type

! Note: Polarization is not 1 when the spin_polar struct represents an ensamble of spins.
! Note: Bmad now uses a (Sx, Sy, Sz) spin representation with quaternions for calculations. 

type spin_polar_struct
  real(rp) :: polarization = 1
  real(rp) :: theta = 0  ! Spherical coords: Angle from z-axis.
  real(rp) :: phi   = 0  ! Spherical coords: Angle in (x,y) plane.
  real(rp) :: xi    = 0  ! Spinor phase angle (See Bmad manual). 
end type

! Structure holding a spin/orbital first order map

type spin_orbit_map1_struct
  real(rp) :: orb_mat(6,6) = 0      ! Orbital matrix
  real(rp) :: vec0(6) = 0           ! Orbital 0th order map: r_out = mat6 * r_in + vec0
  real(rp) :: spin_q(0:3,0:6) = 0   ! 0th and 1st order quaternion spin map
end type

! Structure for linear Invariant Spin Field map

type linear_isf1_struct
  real(rp) :: orb0(6) = 0          ! Closed orbit.
  real(rp) :: isf(0:3, 0:6) = 0    ! Linear ISF map at a given point.
  real(rp) :: s = 0                ! Offset from beginning of element.
  !!! real(rp) :: m_1turn(6,6) = 0   ! Orbital 1-turn matrix.
end type

! Holds ISF info for one Bmad element.

type linear_ele_isf_struct
  type (linear_isf1_struct), allocatable :: node(:)   ! Array per PTC integration node.
end type


!

real(rp), parameter :: x_unit_vec(3) = [1, 0, 0], y_unit_vec(3) = [0, 1, 0], z_unit_vec(3) = [0, 0, 1]

integer, parameter :: magnetic$ = 1, electric$ = 2, mixed$ = 3
character(8), parameter :: em_field_type_name(3) = ['Magnetic', 'Electric', 'Mixed   ']

! Diffraction

integer, parameter :: bragg_diffracted$ = 1, forward_diffracted$ = 2, undiffracted$ = 3
character(20), parameter :: ref_orbit_follows_name(0:3) = [character(20) :: 'GARBAGE!', &
                                             'Bragg_Diffracted', 'Forward_Diffracted', 'Undefracted']

integer, parameter :: reflection$ = 1, transmission$ = 2
character(16), parameter :: mode_name(0:2) = [character(16) :: 'GARBAGE!', 'Reflection', 'Transmission']

! wall3d definitions.

integer, parameter :: anchor_beginning$ = 1, anchor_center$ = 2, anchor_end$ = 3
character(12), parameter :: anchor_pt_name(0:3) = ['GARBAGE! ', 'Beginning', 'Center   ', 'End      ']

! Note: upstream_end$ = entrance_end$ & downstream_end$ = exit_end$ for ele with %orientation = 1.

! first_track_edge$ is the edge a particle enters the element at. 
! This edge will depend upon whether a particle is moving in +s or -s direction.
! Similarly, second_track_edge$ is the edge a particle leaves the element at.

integer, parameter :: none_pt$ = 4
integer, parameter :: entrance_end$ = 1, exit_end$ = 2, both_ends$ = 3, no_end$ = 4, no_aperture$ = 4, nowhere$ = 4
integer, parameter :: continuous$ = 5, surface$ = 6, wall_transition$ = 7

integer, parameter :: upstream_end$ = 1, downstream_end$ = 2
integer, parameter :: inside$ = 3, center_pt$ = 3, start_end$ = 99

integer, parameter :: first_track_edge$ = 11, second_track_edge$ = 12, in_between$ = 13 ! Must be different from upstream_end$, downstream_end$

character(16), parameter :: fiducial_pt_name(4) = [character(16):: &
      'Entrance_end', 'Exit_End', 'Center', 'None']

character(16), parameter :: aperture_at_name(0:8) = [character(16):: &
      'GARBAGE!       ', 'Entrance_End   ', 'Exit_End       ', 'Both_Ends      ', &
      'No_Aperture    ', 'Continuous     ', 'Surface        ', 'Wall_Transition', 'Lord_Defined']

character(16), parameter :: end_at_name(0:4) = [ &
      'GARBAGE!     ', 'Entrance_End ', 'Exit_End     ', 'Both_Ends    ', &
      'No_End       ']

character(16), parameter :: ref_coords_name(0:4) = [ &
      'GARBAGE!     ', 'Entrance_End ', 'Exit_End     ', 'GARBAGE!     ', &
      'No_End       ']

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

integer, parameter :: normal$ = 1, clear$ = 2, opaque$ = 3, wall_start$ = 9, wall_end$ = 10
character(16), parameter :: wall3d_section_type_name(10) = [ &
                     'Normal     ', 'Clear      ', 'Opaque     ', 'Garbage!   ', &
                     'Garbage!   ', 'Garbage!   ', 'Garbage!   ', 'Garbage!   ', &
                     'Wall_Start ', 'Wall_End   ']

! Note: Component order in wall3d_vertex_struct is important since sr3d_read_wall_file uses
! a namelist read to input vertex points

type wall3d_vertex_struct
  real(rp) :: x = 0, y = 0      ! Coordinates of the vertex.
  real(rp) :: radius_x = 0      ! Radius of arc or ellipse x-axis half width. 0 => Straight line.
  real(rp) :: radius_y = 0      ! Ellipse y-axis half height. 
  real(rp) :: tilt = 0          ! Tilt of ellipse
  real(rp) :: angle = 0         ! Angle of (x, y) point.
  real(rp) :: x0 = 0, y0 = 0    ! Center of ellipse
  integer :: type = normal$     ! No longer used.
end type

! A beam pipe or capillary cross-section is a collection of vertexes.
! Vertices are always ordered in increasing angle.
! Note: %surface is not saved in digested files.

integer, parameter :: absolute$ = 1, relative$ = 2, shifted_to_relative$ = 3

type wall3d_section_struct
  character(40) :: name = ''                ! Identifying name
  character(20) :: material = ''            ! Material.
  type (wall3d_vertex_struct), allocatable :: v(:)  ! Array of vertices. Always stored relative.
  type (photon_reflect_surface_struct), pointer :: surface => null()
                                            ! Surface reflectivity tables.
  integer :: type = normal$                 ! normal$, clear$, opaque$, wall_start$, wall_end$
  integer :: n_vertex_input = 0             ! Number of vertices specified by the user.
  integer :: ix_ele = 0                     ! index of lattice element containing section
  integer :: ix_branch = 0                  ! Index of branch lattice element is in.
  integer :: vertices_state = relative$     ! absolute$, or shifted_to_relative$. If set to absolute$ on input, 
                                            !   will be changed to shifted_to_relative$ by section initalizer.
  logical :: patch_in_region = .false.      ! Patch element exists between this section and previous one?
  real(rp) :: thickness = -1                ! Material thickness.
  real(rp) :: s = 0                         ! Longitudinal position
  real(rp) :: r0(2) = 0                     ! Center of section
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

integer, parameter :: chamber_wall$ = 1, mask_plate$ = 2
character(12), parameter :: wall3d_name(2) = [character(12) :: 'Chamber_Wall', 'Mask_Plate']

type wall3d_struct
  character(40) :: name = ''
  integer :: type = chamber_wall$                 ! or mask_plate$
  integer :: ix_wall3d = 0                        ! Index in branch%wall3d(:) array.
  integer :: n_link = 1                           ! For memory management of ele%wall3d
  real(rp) :: thickness = -1                      ! For diffraction_plate elements
  character(20) :: clear_material = ''            !
  character(20) :: opaque_material = ''           !
  logical :: superimpose = .false.                ! Can overlap another wall
  integer :: ele_anchor_pt = anchor_beginning$    ! anchor_beginning$, anchor_center$, or anchor_end$
  type (wall3d_section_struct), allocatable :: section(:) ! Indexed from 1.
end type

!

type taylor_term_struct
  real(rp) :: coef = 0
  integer :: expn(6) = 0
end type

type complex_taylor_term_struct
  complex(rp) :: coef
  integer :: expn(6)  
end type

! Note: the taylor_struct uses the Bmad standard (x, p_x, y, p_y, z, p_z) 
! the universal_taylor in Etienne's PTC uses (x, p_x, y, p_y, p_z, -c*t)
! %ref is the reference point about which the taylor expansion was made.

type taylor_struct
  real (rp) :: ref = 0
  type (taylor_term_struct), pointer :: term(:) => null()
end type

! Note: the complex_taylor_struct uses the Bmad standard (x, p_x, y, p_y, z, p_z) 
! the universal_complex_taylor in Etienne's PTC uses (x, p_x, y, p_y, p_z, -c*t)
! %ref is the reference point about which the complex_taylor expansion was made

type complex_taylor_struct
  complex (rp) :: ref = 0
  type (complex_taylor_term_struct), pointer :: term(:) => null()
end type

! plane list, etc

integer, parameter :: x_plane$ = 1, y_plane$ = 2
integer, parameter :: z_plane$ = 3, n_plane$ = 4, s_plane$ = 5

character(1), parameter :: plane_name(6) = ['X', 'Y', 'Z', 'N', 'S', ' ']
character(2), parameter :: field_plane_name(3) = ['Bx', 'By', 'Bz']

! coordinate def
! Use coord_state_name for getting the string representation of coord%state

integer, parameter :: moving_forward$ = -9
integer, parameter :: pre_born$ = 0    ! EG: before cathode emission. Conforms to OpenPMD standard.
integer, parameter :: alive$ = 1       ! Conforms to OpenPMD standard.
integer, parameter :: lost$ = 2
integer, parameter :: lost_neg_x$ = 3, lost_pos_x$ = 4  
integer, parameter :: lost_neg_y$ = 5, lost_pos_y$ = 6
integer, parameter :: lost_z$ = 7
integer, parameter :: lost_pz$ = 8  ! Particle "turned around" when not tracking with time_runge_kutta.

integer, parameter :: lost_neg_x_aperture$ = 3, lost_pos_x_aperture$ = 4   ! old names.
integer, parameter :: lost_neg_y_aperture$ = 5, lost_pos_y_aperture$ = 6
integer, parameter :: lost_z_aperture$ = 7
integer, parameter :: lost_pz_aperture$ = 8  ! Particle "turned around" when not tracking with time_runge_kutta.

! State_name is not the full list of coord%state possible settings! Missing is not_set$
character(12), parameter ::state_name(0:8) = [character(12):: 'Pre_Born', 'Alive', 'Lost', 'Lost_Neg_X', &
                                                 'Lost_Pos_X', 'Lost_Neg_Y', 'Lost_Pos_Y', 'Lost_Pz', 'Lost_Z']


real(rp), parameter :: vec0$(6) = 0

! The %location component gives the particle location with respect to the element being tracked through
! even if that element is a super_slave or slice_slave. For example, a particle at the beginning of a
! slice_element will have %location = upstream_end$ even though the slice is at the center of the
! encopassing lord element.

type coord_struct                 ! Particle coordinates at a single point
  real(rp) :: vec(6) = 0          ! (x, px, y, py, z, pz). Generally phase space for charged particles. See Bmad manual.
  real(rp) :: s = 0               ! Longitudinal position 
  real(qp) :: t = 0               ! Absolute time (not relative to reference). Note: Quad precision!
  real(rp) :: spin(3) = 0         ! Spin.
  real(rp) :: field(2) = 0        ! Photon E-field intensity (x,y).
  real(rp) :: phase(2) = 0        ! Photon E-field phase (x,y).
  real(rp) :: charge = 0          ! Macroparticle weight (which is different from particle species charge). 
                                  !   For some space charge calcs the weight is in Coulombs.
  real(rp) :: dt_ref = 0          ! Used in:
                                  !   * time tracking for computing z.
                                  !   * by coherent photons = path_length/c_light.
  real(rp) :: r = 0               ! For general use. Not used by Bmad. 
  real(rp) :: p0c = 0             ! For non-photons: Reference momentum.
                                  !     For photons: Photon momentum (not reference).
  real(rp) :: E_potential = 0     ! Potential energy.
  real(rp) :: beta = -1           ! Velocity / c_light.
  integer :: ix_ele = -1          ! Index of the lattice element the particle is in.
                                  !   May be -1 if element is not associated with a lattice.
  integer :: ix_branch = -1       ! Index of the lattice branch the particle is in.
  integer :: ix_turn = 0          ! Turn index for multiturn tracking.
  integer :: ix_user = -1         ! For general use, not used by Bmad.
  integer :: state = not_set$     ! alive$, lost$, lost_neg_x_aperture$, lost_pz$, etc.
  integer :: direction = 1        ! +1 or -1. Sign of longitudinal direction of motion (ds/dt).
                                  !  This is independent of the element orientation.
  integer :: time_dir = 1         ! +1 or -1. Time direction. -1 => Traveling backwards in time.
  integer :: species = not_set$   ! positron$, proton$, etc.  
  integer :: location = upstream_end$  ! upstream_end$, inside$, or downstream_end$
end type

type coord_array_struct
  type (coord_struct), allocatable :: orbit(:)
end type

real(rp), parameter :: no_misalignment$ = 1.0_rp

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

! Wakefield structs...

integer, parameter :: x_polarization$ = 2, y_polarization$ = 3, xy$ = 2
character(8), parameter :: sr_transverse_polarization_name(3) = ['None  ', 'X_Axis', 'Y_Axis']

integer, parameter :: leading$ = 2, trailing$ = 3
integer, parameter :: x_leading$ = 2, y_leading$ = 3, x_trailing$ = 4, y_trailing$ = 5
character(8), parameter :: sr_transverse_position_dep_name(3) = [character(8):: 'none', 'leading', 'trailing']
character(12), parameter :: sr_longitudinal_position_dep_name(5) = &
                [character(12):: 'none', 'x_leading', 'y_leading', 'x_trailing', 'y_trailing']
character(8), parameter :: sr_z_plane_name(5) = [character(8):: 'X', 'XY', 'Y', null_name$, 'Z']

type wake_sr_z_long_struct
  real(rp), allocatable :: w(:)                        ! Input single particle Wake. Indexed from 1.
  complex(rp), allocatable :: fw(:)                    ! Fourier transform of w.
  complex(rp), allocatable :: fbunch(:), w_out(:)      ! Scratch space.
  real(rp) :: dz = 0                                   ! Distance between points. If zero there is no wake.
  real(rp) :: z0 = 0                                   ! Wake extent is [-z0, z0].
  real(rp) :: smoothing_sigma = 0                      ! 0 => No smoothing.
  integer :: position_dependence = none$               ! Transverse: leading$, trailing$, none$
                                                       ! Longitudinal: x_leading$, ..., y_trailing$, none$
  logical :: time_based = .false.                      ! Was input time based?
end type

type wake_sr_mode_struct    ! Psudo-mode Short-range wake struct 
  real(rp) :: amp = 0       ! Amplitude
  real(rp) :: damp = 0      ! Dampling factor.
  real(rp) :: k = 0         ! k factor
  real(rp) :: phi = 0       ! Phase in radians/2pi
  real(rp) :: b_sin = 0     ! non-skew (x) sin-like component of the wake
  real(rp) :: b_cos = 0     ! non-skew (x) cos-like component of the wake
  real(rp) :: a_sin = 0     ! skew (y) sin-like component of the wake
  real(rp) :: a_cos = 0     ! skew (y) cos-like component of the wake
  integer :: polarization = none$            ! Transverse: none$, x_axis$, y_axis$. Not used for longitudinal.
  integer :: position_dependence = not_set$  ! Transverse: leading$, trailing$, none$
                                             ! Longitudinal: x_leading$, ..., y_trailing$, none$
end type

type wake_sr_struct  ! Psudo-mode short-Range Wake struct
  character(400) :: file = ''
  type (wake_sr_z_long_struct) :: z_long
  type (wake_sr_mode_struct), allocatable :: long(:)
  type (wake_sr_mode_struct), allocatable :: trans(:)
  real(rp) :: z_ref_long = 0      ! z reference value for computing the wake amplitude.
  real(rp) :: z_ref_trans = 0     !  This is used to prevent value overflow with long bunches.
  real(rp) :: z_max = 0           ! Max allowable z value. 0-> ignore
  real(rp) :: amp_scale = 1       ! Wake amplitude scale factor.
  real(rp) :: z_scale = 1         ! z-distance scale factor.
  logical :: scale_with_length = .true. ! Scale wake with element length?
end type

! Each wake_lr_struct represents a different mode.
! A non-zero lr_freq_spread attribute value will make freq different from freq_in.

type wake_lr_mode_struct    ! Long-Range Wake struct.
  real(rp) :: freq = 0            ! Actual Frequency in Hz.
  real(rp) :: freq_in = 0         ! Input frequency in Hz.
  real(rp) :: R_over_Q = 0        ! Strength in V/C/m^(2*m_mode).
  real(rp) :: Q = real_garbage$   ! Used for backwards compatability.
  real(rp) :: damp = 0            ! Damping factor = omega / 2 * Q = pi * freq / Q
  real(rp) :: phi = 0             ! Phase in radians/2pi.
  real(rp) :: angle = 0           ! polarization angle (radians/2pi).
  real(rp) :: b_sin = 0           ! non-skew sin-like component of the wake.
  real(rp) :: b_cos = 0           ! non-skew cos-like component of the wake.
  real(rp) :: a_sin = 0           ! skew sin-like component of the wake.
  real(rp) :: a_cos = 0           ! skew cos-like component of the wake.
  integer :: m = 0                ! Mode order (1 = dipole, 2 = quad, etc.)
  logical :: polarized = .false.  ! Polaraized mode?
end type

type wake_lr_struct
  character(400) :: file = ''
  type (wake_lr_mode_struct), allocatable :: mode(:)
  real(rp) :: t_ref = 0             ! time reference value for computing the wake amplitude.
                                    !  This is used to prevent value overflow with long trains.
  real(rp) :: freq_spread = 0       ! Random frequency spread of long range modes.
  real(rp) :: amp_scale = 1         ! Wake amplitude scale factor.
  real(rp) :: time_scale = 1        ! time scale factor.
  logical :: self_wake_on = .true.  ! Long range self-wake used in tracking?
end type

!

type wake_struct
  type (wake_sr_struct) :: sr = wake_sr_struct('', wake_sr_z_long_struct(), null(), null(), 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 1.0_rp, .true.) ! Short-range wake
  type (wake_lr_struct) :: lr = wake_lr_struct('', null(), 0.0_rp, 0.0_rp, 1.0_rp, 1.0_rp, .true.) ! Long-range wake
end type

! ac_kicker structure

type ac_kicker_time_struct
  real(rp) :: amp = 0
  real(rp) :: time = 0
  type (spline_struct) :: spline = spline_struct()
end type

type ac_kicker_freq_struct
  real(rp) :: f = 0
  real(rp) :: amp = 0
  real(rp) :: phi = 0
  integer :: rf_clock_harmonic = 0  ! When RF clock is used.
end type

type ac_kicker_struct
  type (ac_kicker_time_struct), allocatable :: amp_vs_time(:)
  type (ac_kicker_freq_struct), allocatable :: frequency(:)
end type

! Cartesian field decomposition.
! Note: Integer mapping of family and form parameters must be maintained for interface to PTC equivalent.

integer, parameter :: family_y$ = 1, family_x$ = 2, family_qu$ = 3, family_sq$ = 4
integer, parameter :: hyper_y$ = 1, hyper_xy$ = 2, hyper_x$ = 3

character(8), parameter :: cartesian_map_family_name(0:4) = [character(8):: 'GARBAGE!', 'Y', 'X', 'QU', 'SQ']
character(8), parameter :: cartesian_map_form_name(0:3) = ['GARBAGE!', 'Hyper_Y ', 'Hyper_XY', 'Hyper_X ']

type cartesian_map_term1_struct
  real(rp) :: coef = 0
  real(rp) :: kx = 0, ky = 0, kz = 0
  real(rp) :: x0 = 0, y0 = 0, phi_z = 0
  integer :: family = 0        ! family_x$, etc.
  integer :: form = 0          ! hyper_y$, etc.
end type

type cartesian_map_term_struct
  character(400) :: file = ''      ! Input file name. Used also as ID for instances. 
  integer :: n_link = 1            ! For memory management of %term
  type (cartesian_map_term1_struct), allocatable :: term(:)  
end type

type cartesian_map_struct
  real(rp) :: field_scale = 1        ! Factor to scale the fields by
  real(rp) :: r0(3) = 0              ! Field origin offset.
  integer :: master_parameter = 0    ! Master parameter in ele%value(:) array to use for scaling the field.
  integer :: ele_anchor_pt = anchor_beginning$  ! anchor_beginning$, anchor_center$, or anchor_end$
  integer :: field_type = magnetic$  ! or electric$
  type (cartesian_map_term_struct), pointer :: ptr => null()
end type

! Cylindrical map field 

type cylindrical_map_term1_struct
  complex(rp) :: e_coef = 0
  complex(rp) :: b_coef = 0
end type

type cylindrical_map_term_struct
  character(400) :: file = ''   ! Input file name. Used also as ID for instances. 
  integer :: n_link = 1         ! For memory management of this structure
  type (cylindrical_map_term1_struct), allocatable :: term(:)
end type

type cylindrical_map_struct
  integer :: m = 0                  ! Azimuthal Mode: varies as cos(m*phi - theta0_azimuth)
  integer :: harmonic = 0           ! Harmonic of fundamental
  real(rp) :: phi0_fieldmap = 0     ! Mode oscillates as: twopi * (f * t + phi0_fieldmap)
  real(rp) :: theta0_azimuth = 0    ! Azimuthal ((x, y) plane) orientation of mode.
  real(rp) :: field_scale = 1       ! Factor to scale the fields by
  integer :: master_parameter = 0   ! Master parameter in ele%value(:) array to use for scaling the field.
  integer :: ele_anchor_pt = anchor_beginning$  ! anchor_beginning$, anchor_center$, or anchor_end$
  real(rp) :: dz = 0                ! Distance between sampled field points.
  real(rp) :: r0(3) = 0             ! Field origin offset.
  type (cylindrical_map_term_struct), pointer :: ptr => null()
end type

! Generalized gradients

type gen_grad1_struct
  integer :: m = 0                      ! Azimuthal index
  integer :: sincos = 0                 ! sin$ or cos$
  integer :: n_deriv_max = -1           ! Max GG derivative
  ! The derivative matrix is extended to include the interpolating spline polynomial.
  real(rp), allocatable :: deriv(:,:)   ! Range: (iz0:iz1, 0:2*n_deriv_max+1)
end type  

type gen_grad_map_struct
  character(400) :: file = ''   ! Input file name. Used also as ID for instances. 
  type (gen_grad1_struct), allocatable :: gg(:)
  integer :: ele_anchor_pt = anchor_beginning$  ! anchor_beginning$, anchor_center$, or anchor_end$
  integer :: field_type = magnetic$  ! or electric$
  integer :: iz0 = int_garbage$      ! gg%deriv(iz0:iz1, :) lower bound.
  integer :: iz1 = int_garbage$      ! gg%deriv(iz0:iz1, :) upper bound.
  real(rp) :: dz = 0                 ! Point spacing.
  real(rp) :: r0(3) = 0              ! field origin relative to ele_anchor_pt.
  real(rp) :: field_scale = 1        ! Factor to scale the fields by
  integer :: master_parameter = 0    ! Master parameter in ele%value(:) array to use for scaling the field.
  logical :: curved_ref_frame = .false.
end type

! Grid field

type grid_field_pt1_struct
  complex(rp) :: E(3) = 0
  complex(rp) :: B(3) = 0
end type

type grid_field_pt_struct
  character(400) :: file = ''   ! Input file name. Used also as ID for instances. 
  integer :: n_link = 1         ! For memory management of this structure
  type (grid_field_pt1_struct), allocatable :: pt(:,:,:)
end type

type grid_field_struct
  integer :: geometry = 0             ! Type of grid: xyz$, or rotationally_symmetric_rz$
  integer :: harmonic = 0             ! Harmonic of fundamental for AC fields.
  real(rp) :: phi0_fieldmap = 0       ! Mode oscillates as: twopi * (f * t + phi0_fieldmap)
  real(rp) :: field_scale = 1         ! Factor to scale the fields by
  integer :: field_type = mixed$      ! or magnetic$ or electric$
  integer :: master_parameter = 0     ! Master parameter in ele%value(:) array to use for scaling the field.
  integer :: ele_anchor_pt = anchor_beginning$  ! anchor_beginning$, anchor_center$, or anchor_end$
  integer :: interpolation_order = 1  ! Possibilities are 1 or 3.
  real(rp) :: dr(3) = 0   ! Grid spacing.
  real(rp) :: r0(3) = 0   ! Field origin relative to ele_anchor_pt.
  logical :: curved_ref_frame = .false.
  type (grid_field_pt_struct), pointer :: ptr => null()
  type (bicubic_cmplx_coef_struct) bi_coef(4, 2, 3)    ! Save computed coefs for faster tracking
  type (tricubic_cmplx_coef_struct) tri_coef(4, 2, 3)  ! Save computed coefs for faster tracking
end type

! The Taylor field is a set of evenly spaced planes.
! For each plane a taylor series is used to calculate the field components.
! The 2-vector of %expn(2) is (x, y).

type em_taylor_term_struct
  real(rp) :: coef = 0
  integer :: expn(2) = 0
end type

type em_taylor_struct
  real (rp) :: ref = 0
  type (em_taylor_term_struct), allocatable :: term(:)
end type

! Gfortran bug: "field(3) = em_taylor_struct()" not accepted.

! Unit and S-Matrices

real(rp), parameter :: vec3_zero$(3) = 0, vec6_zero$(6) = 0
real(rp), parameter :: mat3_unit$(3,3) = reshape( [1, 0, 0, 0, 1, 0, 0, 0, 1], [3,3])
real(rp), parameter :: mat6_unit$(6,6) = reshape( [1,0,0,0,0,0, 0,1,0,0,0,0, &
                             0,0,1,0,0,0, 0,0,0,1,0,0, 0,0,0,0,1,0, 0,0,0,0,0,1], [6,6])
real(rp), parameter :: s6_unit$(6,6) = reshape( [0,-1,0,0,0,0, 1,0,0,0,0,0, &
                             0,0,0,-1,0,0, 0,0,1,0,0,0, 0,0,0,0,0,-1, 0,0,0,0,1,0], [6,6])

type floor_position_struct
  real(rp) :: r(3) = 0                        ! (x, y, z) offset from origin
  real(rp) :: w(3,3) =  mat3_unit$            ! W matrix. Columns are unit vectors of the frame axes.
  real(rp) :: theta = 0, phi = 0, psi = 0     ! angular orientation consistent with W matrix
end type

! Space charge structure. This structure contains information about the beam as a whole.

type high_energy_space_charge_struct
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
  real(rp) :: eta = 0, etap = 0, deta_ds = 0, sigma = 0, deta_dpz = 0, detap_dpz = 0
end type

! Structure to hold the information of where an individual element is in the lattice.
! Use pointer_to_ele(lat, ele_loc) to point to an element given a lat_ele_loc_struct instance.

type lat_ele_loc_struct
  integer :: ix_ele = -1
  integer :: ix_branch = 0
end type

! Used in lat_ele_order_branch_struct

type lat_ele_order1_struct
  integer :: ix_branch = -1     ! Branch index
  integer :: ix_order = -1      ! Order index. -1 -> Unique in lattice, 0 -> unique in branch.
end type

type lat_ele_order_array_struct
  type (lat_ele_order1_struct), allocatable :: ele(:)
end type

! Structure for holding the order index of elements in a lattice.
! For example, If an element named "QQ" is in lattice branch #2 and is the fifth element
! in the branch with name "QQ", the order index for this element is 5.
! That is, a unique name for this element is "2>>QQ##5". 
! Note: The super_lord elements, which live in branch #0, are associated.
! with the branch that its super_slaves live in.
! That is, if a marker is superimpsed upon "2>>QQ##5", the "QQ" lord element will be placed in the lord
! section of branch 0 but this lord's order index, associated branch, and unique name are not changed.
!
! Given:
!   type (lat_ele_order_struct) order
! Then for a given lattice element in branch ix_branch and with element index ix_ele:
!   order%branch(ix_branch)%ele(ix_ele)%ix_order   => Order index. -1 -> Unique in lattice, 0 -> unique in branch.
!   order%branch(ix_branch)%ele(ix_ele)%ix_branch  => Associated branch index

type lat_ele_order_struct
  type (lat_ele_order_array_struct), allocatable :: branch(:)
end type

! Structure to be used for an array of pointers to lattice elements.
! The id component is not set by any Bmad routines and can be used, for example, by programs that 
! handle multiple lattices to indicate which lattice the element pointer is pointing to.
! A pointer to an element in a lattice is not usable if the number of elements in the 
! lattice is modified (so that the lattice element array is reallocated). In this case,
! the %loc component is potentially useful (as long as the element pointed to does not move).

type ele_pointer_struct
  type (ele_struct), pointer :: ele => null()
  type (lat_ele_loc_struct) :: loc = lat_ele_loc_struct()
  integer :: id = -1                    ! For general use. Not used by Bmad.
end type

! Structure to be used for an array of pointers to branches.

type branch_pointer_struct
  type (branch_struct), pointer :: branch => null()
end type

! Structure to be used for an array of pointers to lattices.

type lat_pointer_struct
  type (lat_struct), pointer :: lat => null()
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
integer, parameter :: rad_int_group$ = 7, all_groups$ = 8, s_and_floor_position_group$ = 9

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
  logical :: has_misalign = .false.   ! Used to avoid unnecessary calls to offset_particle.
end type

! Cache multipole values in the element. 
! In one simulation 25% of the time was spent constructing multipole arrays.

type multipole_cache_struct
  ! a_kick_mag and b_kick_mag are for non-multipole parameters like k1, k2, hkick, etc.
  real(rp), allocatable :: a_pole_mag(:), b_pole_mag(:)
  real(rp), allocatable :: a_kick_mag(:), b_kick_mag(:)
  integer :: ix_pole_mag_max = -1, ix_kick_mag_max = -1
  logical :: mag_valid = .false.
  ! From elseparator hkick and vkick.
  real(rp), allocatable :: a_pole_elec(:), b_pole_elec(:)
  real(rp), allocatable :: a_kick_elec(:), b_kick_elec(:)
  integer :: ix_pole_elec_max = -1, ix_kick_elec_max = -1
  logical :: elec_valid = .false.
end type

! Radiation damping and stochastic maps

type rad_map_struct
  real(rp) :: ref_orb(6) = -1                 ! Reference point around which damp_mat is calculated.
  real(rp) :: damp_dmat(6,6) = 0              ! damp_correction = xfer_mat_with_damping - xfer_mat_without_damping.
  real(rp) :: xfer_damp_vec(6) = 0            ! Transfer map with damping 0th order vector.
  real(rp) :: xfer_damp_mat(6,6) = mat6_unit$ ! 1st order matrix: xfer_no_damp_mat + xfer_damp_correction.
  real(rp) :: stoc_mat(6,6) = 0               ! Stochastic variance or "kick" (Cholesky decomposed) matrix.
end type

type rad_map_ele_struct
  ! 6D emit. In this structure rm0%stoc_mat and rm1%stoc_mat are the kick matrices.
  type (rad_map_struct) rm0, rm1  ! Upstream half and downstream half matrices for an element.
  logical :: stale = .true.
end type

! Segmented grid

type surface_segmented_pt_struct
  real(rp) :: x0 = 0, y0 = 0, z0 = 0         ! Position at center
  real(rp) :: dz_dx = 0, dz_dy = 0           ! Slope at center
end type

type surface_segmented_struct
  logical :: active = .false.
  real(rp) :: dr(2) = 0, r0(2) = 0
  type (surface_segmented_pt_struct), allocatable :: pt(:,:) 
end type

! H_misalign grid

type surface_h_misalign_pt_struct
  real(rp) :: x0 = 0, y0 = 0                 ! Position at center
  real(rp) :: rot_y = 0, rot_t = 0, rot_y_rms = 0, rot_t_rms = 0  ! rot_t = x-rotation for Bragg and z-rotation for Laue.
end type

type surface_h_misalign_struct
  logical :: active = .false.
  real(rp) :: dr(2) = 0, r0(2) = 0
  type (surface_h_misalign_pt_struct), allocatable :: pt(:,:) 
end type

! displacement grid

type surface_displacement_pt_struct
  real(rp) :: x0 = 0, y0 = 0                 ! Position at center
  real(rp) :: z0 = 0, dz_dx = 0, dz_dy = 0, d2z_dxdy = 0
end type

type surface_displacement_struct
  logical :: active = .false.
  real(rp) :: dr(2) = 0, r0(2) = 0
  type (surface_displacement_pt_struct), allocatable :: pt(:,:) 
end type

! Photon statistics at a detector

type pixel_pt_struct
  integer(8) :: n_photon = 0
  complex(rp) :: E_x  = 0, E_y = 0
  real(rp) ::  intensity_x = 0, intensity_y = 0, intensity = 0
  real(rp) :: orbit(6) = 0            ! x, Vx/c, y, Vy/c, dummy, E - E_ref.
  real(rp) :: orbit_rms(6) = 0        ! RMS statistics.
  real(rp) :: init_orbit(6) = 0       ! Initial orbit at start of lattice statistics.
  real(rp) :: init_orbit_rms(6) = 0   ! Initial orbit at start of lattice RMS statistics.
end type

type pixel_detec_struct
  real(rp) :: dr(2) = 0, r0(2) = 0
  integer(8) :: n_track_tot = 0       ! How many photons were launched from source element.
  integer(8) :: n_hit_detec = 0         ! How many photons hit the detector.
  integer(8) :: n_hit_pixel = 0       ! How many photons hit the pixel grid of the detector.
  type (pixel_pt_struct), allocatable :: pt(:,:) ! Grid of pixels
end type

! Surface container structure

type surface_curvature_struct
  real(rp) :: xy(0:6,0:6) = 0
  real(rp) :: spherical = 0
  real(rp) :: elliptical(3) = 0             ! Total curvature = elliptical + spherical
  logical :: has_curvature = .false.        ! Dependent var. Will be set by Bmad
end type

! Target points are in element coordinates.

type target_point_struct
  real(rp) :: r(3) = 0   ! (x, y, z)
end type

type photon_target_struct
  integer :: type = off$                      ! or rectangular$
  integer :: n_corner = 0
  type (lat_ele_loc_struct) :: ele_loc = lat_ele_loc_struct()
  type (target_point_struct) :: corner(8) = target_point_struct()
  type (target_point_struct) :: center = target_point_struct()
end type

type photon_material_struct
  complex(rp) :: f0_m1 = 0, f0_m2 = 0 ! For multilayer_mirror only.
  complex(rp) :: f_0 = 0
  complex(rp) :: f_h = 0     ! Structure factor for H direction.
  complex(rp) :: f_hbar = 0  ! Structure factor for -H direction.
  complex(rp) :: f_hkl = 0   ! = sqrt(f_h * f_hbar)
  real(rp) :: h_norm(3) = 0  ! Normalized H vector for crystals.
  real(rp) :: l_ref(3) = 0   ! Crystal reference orbit displacement vector in element coords.
end type

integer, parameter :: polarized$ = 1, unpolarized$ = 2
character(8), parameter :: polarization_name(3) = [character(8):: 'SIGMA', 'PI', 'BOTH']

! photon_element_struct is an ele_struct component holding photon parameters

type photon_element_struct
  type (surface_curvature_struct) :: curvature = surface_curvature_struct()
  type (photon_target_struct) :: target = photon_target_struct()
  type (photon_material_struct) :: material = photon_material_struct()
  type (surface_segmented_struct) :: segmented = surface_segmented_struct(.false., 0, 0, null())
  type (surface_h_misalign_struct) :: h_misalign = surface_h_misalign_struct(.false., 0, 0, null())
  type (surface_displacement_struct) :: displacement = surface_displacement_struct(.false., 0, 0, null())
  type (pixel_detec_struct) :: pixel = pixel_detec_struct([0.0_rp, 0.0_rp], [0.0_rp, 0.0_rp], 0, 0, 0, null())
  integer :: reflectivity_table_type = not_set$
  type (photon_reflect_table_struct) reflectivity_table_sigma  ! If polarization is ignored use sigma table.
  type (photon_reflect_table_struct) reflectivity_table_pi
  type (spline_struct), allocatable :: init_energy_prob(:)  ! Initial energy probability density
  real(rp), allocatable :: integrated_init_energy_prob(:)
end type

!------------------------------------------------------------------------------
! Beam structures

type bunch_struct
  type (coord_struct), allocatable :: particle(:)
  integer, allocatable :: ix_z(:)  ! bunch%ix_z(1) is index of head particle, etc.
  real(rp) :: charge_tot = 0       ! Total charge in a bunch (Coul).
  real(rp) :: charge_live = 0      ! Charge of live particles (Coul).
  real(rp) :: z_center = 0         ! Longitudinal center of bunch at creation time. Note: Generally, z_center of 
                                   !   bunch #1 is 0 and z_center of the other bunches is negative.
  real(rp) :: t_center = 0         ! Center of bunch at creation time relative to head bunch.
  real(rp) :: t0 = real_garbage$   ! Used by track1_bunch_space_charge for tracking so particles have constant t.
  logical :: drift_between_t_and_s = .false.
                                   ! Drift (ignore any fields) instead of tracking to speed up the calculation?
                                   ! This can only be done under certain circumstances.
  integer :: ix_ele = 0            ! Nominal element bunch is at. But, EG, dead particles can be someplace else.
  integer :: ix_bunch = 0          ! Bunch index. Head bunch = 1, etc.
  integer :: ix_turn = 0           ! Turn index for long term tracking. ix_turn = 0 before end of first turn, etc.
  integer :: n_live = 0
  integer :: n_good = 0            ! Number of accepted steps when using adaptive step size control.
  integer :: n_bad = 0             ! Number of rejected steps when using adaptive step size control.
end type

type beam_struct
  type (bunch_struct), allocatable :: bunch(:)
end type

type ellipse_beam_init_struct
  integer :: part_per_ellipse = 0  ! number of particles per ellipse
  integer :: n_ellipse = 1         ! number of ellipses (>= 1)
  real(rp) :: sigma_cutoff = 0     ! sigma cutoff of the representation
end type

type kv_beam_init_struct
  integer :: part_per_phi(2) = 0    ! number of particles per angle variable.
  integer :: n_I2 = 0               ! number of I2
  real(rp) :: A = 0                 ! A = I1/e
end type

type grid_beam_init_struct
  integer :: n_x = 0         ! Number of columns.
  integer :: n_px = 0        ! Number of rows.
  real(rp) :: x_min = 0      ! Lower x limit.
  real(rp) :: x_max = 0      ! Upper x limit.
  real(rp) :: px_min = 0     ! Lower px limit.
  real(rp) :: px_max = 0     ! Upper px limit.
end type

type beam_init_struct
  character(400) :: position_file = ''                ! File with particle positions.
  character(16) :: distribution_type(3) = 'RAN_GAUSS' ! distribution type (in x-px, y-py, and z-pz planes)
                                             ! "ELLIPSE", "KV", "GRID", "FILE", "RAN_GAUSS" or "" = "RAN_GAUSS" 
  real(rp) :: spin(3) = 0                    ! Spin (x, y, z)
  type (ellipse_beam_init_struct) :: ellipse(3) = ellipse_beam_init_struct() ! Ellipse beam distribution
  type (kv_beam_init_struct) :: KV = kv_beam_init_struct()                   ! KV beam distribution
  type (grid_beam_init_struct) :: grid(3) = grid_beam_init_struct()          ! Grid beam distribution
  real(rp) :: center_jitter(6) = 0.0         ! Bunch center rms jitter
  real(rp) :: emit_jitter(2)   = 0.0         ! a and b bunch emittance rms jitter normalized to emittance
  real(rp) :: sig_z_jitter     = 0.0         ! bunch length RMS jitter
  real(rp) :: sig_pz_jitter    = 0.0         ! RMS pz spread jitter 
  integer :: n_particle = 0                  ! Number of particles per bunch.
  logical :: renorm_center = .true.          ! Renormalize centroid?
  logical :: renorm_sigma = .true.           ! Renormalize sigma?
  character(16) :: random_engine = 'pseudo'  ! Or 'quasi'. Random number engine to use. 
  character(16) :: random_gauss_converter = 'exact'  
                                             ! Or 'quick'. Uniform to gauss conversion method.
  real(rp) :: random_sigma_cutoff = -1       ! Cut-off in sigmas.
  real(rp) :: a_norm_emit = 0                ! a-mode normalized emittance (emit * beta * gamma)
  real(rp) :: b_norm_emit = 0                ! b-mode normalized emittance (emit * beta * gamma)
  real(rp) :: a_emit = 0                     ! a-mode emittance
  real(rp) :: b_emit = 0                     ! b-mode emittance
  real(rp) :: dPz_dz = 0                     ! Correlation of Pz with long position.
  real(rp) :: center(6) = 0                  ! Bench phase space center offset relative to reference.
  real(rp) :: t_offset = 0                   ! Time center offset
  real(rp) :: dt_bunch = 0                   ! Time between bunches.
  real(rp) :: sig_z = 0                      ! Z sigma in m.
  real(rp) :: sig_pz = 0                     ! pz sigma
  real(rp) :: bunch_charge = 0               ! charge (Coul) in a bunch.
  integer :: n_bunch = 0                     ! Number of bunches.
  integer :: ix_turn = 0                     ! Turn index used to adjust particles time if needed.
  character(16) :: species = ""              ! "positron", etc. "" => use referece particle.
  logical :: full_6D_coupling_calc = .false. ! Use V from 6x6 1-turn mat to match distribution?  
                                             !   Else use 4x4 1-turn mat used.
  logical :: use_particle_start = .false.    ! Use lat%particle_start instead of beam_init%center, %spin?
  logical :: use_t_coords = .false.          ! If true, the distributions will be taken as in t-coordinates  
  logical :: use_z_as_t   = .false.          ! Only used if  use_t_coords = .true.
                                             !   If true,  z describes the t distribution 
                                             !   If false, z describes the s distribution
  character(200) :: file_name = ''           ! OLD!! DO NOT USE!!
end type


character(8), parameter :: random_engine_name(2) = [character(8):: 'pseudo', 'quasi'] ! Case sensitive
character(8), parameter :: random_gauss_converter_name(2) = [character(8):: 'exact', 'quick'] ! Case sensitive
character(12), parameter :: beam_distribution_type_name(5) = [character(12):: &
                                    'Ellipse', 'KV', 'Grid', 'File', 'Ran_Gauss']

! The routines calc_bunch_params and calc_bunch_params_slice calculate bunch parameters.
! Note: If, for example, there is only one particle, %twiss_valid = False but centroid calc is OK.

type bunch_params_struct
  type (coord_struct) :: centroid = coord_struct()  ! Lab frame
  type (twiss_struct) :: x = twiss_struct(), y = twiss_struct(), z = twiss_struct() ! Projected Twiss parameters
  type (twiss_struct) :: a = twiss_struct(), b = twiss_struct(), c = twiss_struct() ! Normal mode twiss parameters
  real(rp) :: sigma(6,6) = 0             ! beam size matrix
  real(rp) :: rel_max(7) = 0             ! Max orbit relative to centroid. 7 -> time.
  real(rp) :: rel_min(7) = 0             ! Min orbit relative to_centroid. 7 -> time.
  real(rp) :: s = -1                     ! Longitudinal position.
  real(rp) :: t = -1                     ! Time.
  real(rp) :: sigma_t = 0                ! RMS of time spread.
  real(rp) :: charge_live = 0            ! Charge of all non-lost particle
  real(rp) :: charge_tot = 0             ! Charge of all particles.
  integer :: n_particle_tot = 0          ! Total number of particles
  integer :: n_particle_live = 0         ! Number of non-lost particles
  integer :: n_particle_lost_in_ele = 0  ! Number lost in element (not calculated by Bmad)
  integer :: n_good_steps = 0            ! Number of good steps (set when tracking with space charge)
  integer :: n_bad_steps = 0             ! Number of bad steps (set when tracking with space charge)
  integer :: ix_ele = -1                 ! Lattice element where params evaluated at.
  integer :: location = not_set$         ! Location in element: upstream_end$, inside$, or downstream_end$
  logical :: twiss_valid = .false.       ! Is the data here valid? Note: IF there is no energy
                                         !   variation (RF off) twiss_valid may be true but in
                                         !   this case the z-twiss will not be valid.
end type

! Bunch_track_struct is the bunch analogue of a particle track_struct

type bunch_track_struct
  type (bunch_params_struct), allocatable :: pt(:)     ! Array indexed from 0
  real(rp) :: ds_save = -1                             ! Min distance between points.
  integer :: n_pt = -1                                 ! Track upper bound
end type

!------------------------------------------------------------------------

! Converter structure for Probability(pc_out, r_out; pc_in, thickness)
! Note: Probability at r = 0 is zero since using polar coordinates.

type converter_prob_pc_r_struct
  real(rp), allocatable :: pc_out(:)        ! Grid pc_out values.
  real(rp), allocatable :: r(:)             ! Grid r_out values.
  real(rp), allocatable :: prob(:,:)        ! Probability grid.
  real(rp), allocatable :: spin_z(:,:)      ! Z polarization grid.
  ! Stuff below is calculated rather than read in from the lattice file.
  real(rp) pc_out_min, pc_out_max
  real(rp) :: integrated_prob = 0           ! Integrated probability over (pc_out, r) with restrictions factered in.
  real(rp), allocatable :: p_norm(:,:)      ! Normalized probability taking into account.
                                            !   angle_out_max, pc_out_min, and pc_out_max restrictions.
  real(rp), allocatable :: integ_pc_out(:)  ! Normalized probability integrated from min pc_out up.
  real(rp), allocatable :: integ_r(:,:)
  real(rp), allocatable :: integ_r_ave(:)
end type

! coef fits

type converter_dir_1D_struct
  real(rp) :: pc_out = 0          ! pc_out value at fit
  real(rp) :: poly(0:4) = 0       ! param(r) = Sum: poly(i) * r^i
end type

type converter_dir_2D_struct
  real(rp) :: k = 0
  real(rp) :: poly(0:3) = 0
end type

type converter_dir_coef_struct
  type (converter_dir_1D_struct), allocatable :: fit_1d_r(:)
  type (converter_dir_2D_struct) :: fit_2d_r
  type (converter_dir_2D_struct) :: fit_2d_pc
  real(rp) :: c0 = 0
end type

type converter_direction_out_struct
  type (converter_dir_coef_struct) :: beta
  type (converter_dir_coef_struct) :: alpha_x, alpha_y
  type (converter_dir_coef_struct) :: dxds_min, dxds_max, dyds_max
  type (converter_dir_coef_struct) :: c_x
end type

! Converter structure for a given incoming particle energy.

type converter_sub_distribution_struct
  real(rp) :: pc_in = -1
  real(rp) :: spin_in(3)
  type (converter_prob_pc_r_struct) :: prob_pc_r
  type (converter_direction_out_struct) :: dir_out
end type

! Foil structure

type material_struct
  integer :: species = not_set$
  integer :: number = int_garbage$                      ! Relative number
  real(rp) :: density = real_garbage$, density_used = real_garbage$
  real(rp) :: area_density = real_garbage$, area_density_used = real_garbage$
  real(rp) :: radiation_length = real_garbage$, radiation_length_used = real_garbage$
end type

type foil_struct
  type (material_struct), allocatable :: material(:)
end type

! Distribution of outgoing particles for a given thickness.

type converter_distribution_struct
  real(rp) :: thickness = -1
  type (converter_sub_distribution_struct), allocatable :: sub_dist(:) ! Distribution at various pc_in values.
end type

! Converter structure
! For things like a positron converter in a linac.

type converter_struct
  integer :: species_out = 0    ! Output species
  character(40) :: material_type = ''
  type (converter_distribution_struct), allocatable :: dist(:)  ! Distribution at various thicknesses 
end type

! Struct for element to element control.
! Note: Unlike the old days, not all attributes have an associated index. 
!   So %ix_attrib may be -1. Using pointer_to_attribute with %attribute will always work.

type control_struct
  real(rp) :: value = 0          ! Used by group, and overlay elements.
  real(rp), allocatable :: y_knot(:)
  type (expression_atom_struct), allocatable :: stack(:) ! Evaluation stack
  type (lat_ele_loc_struct) :: slave = lat_ele_loc_struct()
  type (lat_ele_loc_struct) :: lord = lat_ele_loc_struct()
  character(40) :: slave_name = ''   ! Name of slave.
  character(40) :: attribute = ''    ! Name of attribute controlled. Set to "FIELD_OVERLAPS" for field overlaps.
                                     !   Set to "INPUT" or "OUTPUT" for feedback slaves.
  integer :: ix_attrib = -1          ! Index of attribute controlled. See note above!
end type

type control_var1_struct
  character(40) :: name = ''
  real(rp) :: value = 0
  real(rp) :: old_value = 0
end type

type control_ramp1_struct
  real(rp), allocatable :: y_knot(:)
  type (expression_atom_struct), allocatable :: stack(:) ! Evaluation stack
  character(40) :: attribute = ''     ! Name of attribute controlled. Set to "FIELD_OVERLAPS" for field overlaps.
  character(40) :: slave_name = ''    ! Name of slave.
  logical :: is_controller = .false.  ! Is the slave a controller? If so bookkeeping is different.
end type

integer, parameter :: cubic$ = 3
character(8), parameter :: interpolation_name(4) = [character(8):: null_name$, 'null_name$', 'Cubic', 'Linear']

type ramper_lord_struct
  integer :: ix_ele = 0       ! Lord index
  integer :: ix_con = 0       ! Index in lord%control%ramp(:) array
  real(rp), pointer :: attrib_ptr => null()    ! Pointer to attribute in this element.
end type

type controller_struct
  type (control_var1_struct), allocatable :: var(:)
  type (control_ramp1_struct), allocatable :: ramp(:)            ! For ramper lord elements
  type (ramper_lord_struct), allocatable :: ramper_lord(:)       ! Ramper lord info for this slave
  real(rp), allocatable :: x_knot(:)
end type


!-------------------------------------------------------------------------
! Ele_struct:
! Remember: If this struct is changed you have to:
!     Increase bmad_inc_version by 1.
!     run scripts to regenerate cpp_bmad_interface library.
!     Modify (this is not a complete list):
!       read_digested_bmad_file
!       write_digested_bmad_file
!       parsing routines to read in modified/new parameters...
!       deallocate_ele_pointers
!       ele_equal_ele
!       type_ele
!       write_bmad_lattice_file
!       pointer_to_attribute
!       pointers_to_attribute
!       Bmad manual

type ele_struct
  character(40) :: name = '<Initialized>'                ! name of element.
  character(40) :: type = ''                             ! type name.
  character(40) :: alias = ''                            ! Another name.
  character(40) :: component_name = ''                   ! Used by overlays, multipass patch, etc.
  character(200), pointer :: descrip => null()           ! Description string.
  type (twiss_struct) :: a = twiss_struct()              ! Twiss parameters at end of element
  type (twiss_struct) :: b = twiss_struct()              ! Twiss parameters at end of element
  type (twiss_struct) :: z = twiss_struct()              ! Twiss parameters at end of element
  type (xy_disp_struct) :: x = xy_disp_struct()          ! Projected dispersions.
  type (xy_disp_struct) :: y = xy_disp_struct()          ! Projected dispersions.
  type (ac_kicker_struct), pointer :: ac_kick => null()  ! ac_kicker element parameters.
  type (bookkeeping_state_struct) :: bookkeeping_state = bookkeeping_state_struct() ! Attribute bookkeeping
  type (branch_struct), pointer :: branch => null()                      ! Pointer to branch containing element.
  type (controller_struct), pointer :: control => null()                 ! group & overlay variables.
  type (converter_struct), pointer :: converter => null()                ! EG: Positron converter in linac.
  type (foil_struct), pointer :: foil => null()
  type (ele_struct), pointer :: lord => null()                           ! Pointer to a slice lord.
  type (fibre), pointer :: ptc_fibre => null()                           ! PTC track corresponding to this ele.
  type (floor_position_struct) :: floor = floor_position_struct(vec3_zero$, mat3_unit$, 0.0_rp, 0.0_rp, 0.0_rp)
  type (high_energy_space_charge_struct), pointer :: high_energy_space_charge => null()
  type (mode3_struct), pointer :: mode3 => null()                        ! 6D normal mode structure.
  type (photon_element_struct), pointer :: photon => null()
  type (multipole_cache_struct), allocatable :: multipole_cache
  type (rad_map_ele_struct), pointer :: rad_map => null()                ! Radiation kick parameters
  ! Note: The reference orbits for spin and orbit Taylor maps are not necessarily the same.
  ! For example, Sprint spin Taylor maps can be with respect to the zero orbit independent of the orbital map.
  type (taylor_struct) :: taylor(6) = taylor_struct()                    ! Phase space Taylor map.
  real(rp) :: spin_taylor_ref_orb_in(6) = real_garbage$
  type (taylor_struct) :: spin_taylor(0:3) = taylor_struct()             ! Quaternion Spin Taylor map.
  type (wake_struct), pointer :: wake => null()                          ! Wakes
  type (wall3d_struct), pointer :: wall3d(:) => null()                   ! Chamber or capillary wall
  ! E/M field structs.
  type (cartesian_map_struct), pointer :: cartesian_map(:) => null()     ! Used to define E/M fields
  type (cylindrical_map_struct), pointer :: cylindrical_map(:) => null() ! Used to define E/M fields
  type (gen_grad_map_struct), pointer :: gen_grad_map(:) => null()       ! Used to define E/M fields.
  type (grid_field_struct), pointer :: grid_field(:) => null()           ! Used to define E/M fields.
  ! The difference between map_ref_orb and time_ref_orb is that map_ref_orb is the reference orbit for the
  ! 1st order spin/orbit map which, in general, is non-zero while time_ref_orb follows the reference particle which is
  ! generally the zero orbit (non-zero, for example, in the second slice of a sliced wiggler).
  type (coord_struct) :: map_ref_orb_in = coord_struct()       ! Entrance end transfer map ref orbit
  type (coord_struct) :: map_ref_orb_out = coord_struct()      ! Exit end transfer map ref orbit
  type (coord_struct) :: time_ref_orb_in = coord_struct()      ! Reference orbit at entrance end for ref_time calc.
  type (coord_struct) :: time_ref_orb_out = coord_struct()     ! Reference orbit at exit end for ref_time calc.
  real(rp) :: value(num_ele_attrib$) = 0                       ! attribute values.
  real(rp) :: old_value(num_ele_attrib$) = 0                   ! Used to see if %value(:) array has changed.
  ! Note: The reference orbit for spin/orbit matrices is %map_ref_orb_in/out
  real(rp) :: spin_q(0:3,0:6) = real_garbage$                  ! 0th and 1st order Spin transport quaternion.
  real(rp) :: vec0(6) = 0                                      ! 0th order transport vector.
  real(rp) :: mat6(6,6) = 0                                    ! 1st order transport matrix.
  real(rp) :: c_mat(2,2) = 0                                   ! 2x2 C coupling matrix
  real(rp) :: gamma_c = 1                                      ! gamma associated with C matrix
  real(rp) :: s_start = 0                                      ! longitudinal ref position at entrance_end
  real(rp) :: s = 0                                            ! longitudinal ref position at the exit end.
  real(rp) :: ref_time = 0                                     ! Time ref particle passes exit end.
  real(rp), pointer :: a_pole(:) => null()                     ! knl for multipole elements.
  real(rp), pointer :: b_pole(:) => null()                     ! tilt for multipole elements.
  real(rp), pointer :: a_pole_elec(:) => null()                ! Electrostatic multipoles. ksnl for multipole elements.
  real(rp), pointer :: b_pole_elec(:) => null()                ! Electrostatic multipoles.
  real(rp), pointer :: custom(:) => null()                     ! Custom attributes.
  real(rp), pointer :: r(:,:,:) => null()                      ! For general use. Not used by Bmad.
  integer :: key = 0                              ! Element class (quadrupole, etc.).
  integer :: sub_key = 0                          ! Records bend input type.
  integer :: ix_ele = -1                          ! Index in branch ele(0:) array. Set to ix_slice_slave$ = -2 for slice_slave$ elements.
  integer :: ix_branch = 0                        ! Index in lat%branch(:) array. Note: lat%ele => lat%branch(0).
  integer :: lord_status = not_a_lord$            ! Type of lord element this is. overlay_lord$, etc.
  integer :: n_slave = 0                          ! Number of slaves (except field overlap slaves) of this element.
  integer :: n_slave_field = 0                    ! Number of field slaves of this element.
  integer :: ix1_slave = 0                        ! Pointer index to this element's slaves.
  integer :: slave_status = free$                 ! Type of slave element this is. multipass_slave$, slice_slave$, etc.
  integer :: n_lord = 0                           ! Number of lords (except field overlap and ramper lords).
  integer :: n_lord_field = 0                     ! Number of field lords of this element.
  integer :: n_lord_ramper = 0                    ! Number of ramper lords.
  integer :: ic1_lord = 0                         ! Pointer index to this element's lords.
  integer :: ix_pointer = 0                       ! For general use. Not used by Bmad.
  integer :: ixx = 0, iyy = 0, izz = 0            ! Index for Bmad internal use.
  integer :: mat6_calc_method = bmad_standard$    ! taylor$, symp_lie_ptc$, etc.
  integer :: tracking_method = bmad_standard$     ! taylor$, linear$, etc.
  integer :: spin_tracking_method = tracking$     ! symp_lie_ptc$, etc.
  integer :: csr_method = off$                    ! or one_dim$ ("1_dim"), steady_state_3d$
  integer :: space_charge_method = off$           ! slice$, slice_longitudinal$, slice_transverse$, fft_3D$, cathode_fft_3d$
  integer :: ptc_integration_type = matrix_kick$  ! drift_kick$, matrix_kick$, or ripken_kick$
  integer :: field_calc = bmad_standard$          ! no_field$, fieldmap$, refer_to_lords$, or custom$
  integer :: aperture_at = exit_end$              ! Aperture location: entrance_end$, ...
  integer :: aperture_type = rectangular$         ! rectangular$, elliptical$, auto_aperture$, ...
  integer :: ref_species = not_set$               ! Reference species
  integer :: orientation = 1                 ! -1 -> Element is longitudinally reversed. +1 -> Normal.
  logical :: symplectify = .false.           ! Symplectify mat6 matrices.
  logical :: mode_flip = .false.             ! Have the normal modes traded places?
  logical :: multipoles_on = .true.          ! For turning multipoles on/off
  logical :: scale_multipoles = .true.       ! Are ab_multipoles within other elements (EG: quads, etc.) 
                                             !        scaled by the strength of the element?
  logical :: taylor_map_includes_offsets = .true. ! Taylor map calculated with element misalignments?
  logical :: field_master = .false.          ! Calculate strength from the field value?
  logical :: is_on = .true.                  ! For turning element on/off.
  logical :: logic = .false.                 ! For general use. Not used by Bmad (except during lattice parsing).
  logical :: bmad_logic = .false.            ! For Bmad internal use only.
  logical :: select = .false.                ! For Bmad internal use only.
  logical :: offset_moves_aperture = .false. ! element offsets affects aperture?
contains
  procedure next_in_branch
  !! final :: ele_finalizer
end type

! The lat_param_struct should be called the branch_param_struct [Present name is a historical artifact.]

type lat_param_struct
  real(rp) :: n_part = 0                       ! Particles/bunch (for BeamBeam elements).
  real(rp) :: total_length = 0                 ! total_length of branch. Warning: branch may not start at s = 0.
  real(rp) :: unstable_factor = 0              ! If positive: Growth rate/turn if unstable in closed branches or
                                               !   |orbit-aperture|/aperture if particle hits wall. Zero otherwise.
  real(rp) :: t1_with_RF(6,6) = 0              ! Full 1-turn matrix with RF on.
  real(rp) :: t1_no_RF(6,6) = 0                ! Full 1-turn matrix with RF off.
  real(rp) :: spin_tune = 0                    ! Closed orbit spin tune.
  integer :: particle = not_set$               ! Reference particle: positron$, electron$, etc.
                                               !     Call lattice_bookkeeper if this is changed.
  integer :: default_tracking_species = ref_particle$  ! Default particle type to use in tracking.
  integer :: geometry = 0                      ! open$ or closed$
  integer :: ixx = 0                           ! Integer for general use
  logical :: stable = .false.                  ! is closed lat stable?
  logical :: live_branch = .true.              ! Should tracking be done on the branch?
  real(rp) :: g1_integral = -1                 ! Approximate |g| (bending strength) integral of branch. 
  real(rp) :: g2_integral = -1                 ! Approximate g^2 integral of branch. 
  real(rp) :: g3_integral = -1                 ! Approximate g^2 integral of branch. 
  type (bookkeeping_state_struct) :: bookkeeping_state = bookkeeping_state_struct()
                                               ! Overall status for the branch.
  type (beam_init_struct) :: beam_init = beam_init_struct() ! For beam initialization.
end type

! Structure for linking a branch_struct with a collection of ptc layouts

type ptc_layout_pointer_struct
  type (layout), pointer :: ptr => null()
end type

type ptc_branch1_struct
  type (ptc_layout_pointer_struct), allocatable :: m_u_layout(:)
  type (layout), pointer :: m_t_layout => null()      ! Tracking layout.
end type

!

type mode_info_struct
  logical :: stable = .false.  ! Is the mode stable?
  real(rp) :: tune   = 0       ! "fractional" tune in radians
  real(rp) :: emit   = 0       ! Emittance (unnormalized).
  real(rp) :: chrom  = 0       ! Chromaticity.
  real(rp) :: sigma  = 0       ! Beam size.
  real(rp) :: sigmap = 0       ! Beam divergence.
end type

! Resonance Driving (RD) term

type resonance_h_struct
  character(6) :: id = ''   ! 6 digit ID. EG: '003100'
  complex(rp) :: c_val = 0  ! Resonance value
end type

! Normal form components using Bmad structures and PTC structures.
! M = A1 o c_inv o L exp(F.grad)Identity o c o A1_inv
! A1 and L are linear, and c maps to the phasor basis: h+ = x + i p, h- = x - i p
! See subroutines: normal_form_taylors and normal_form_complex_taylors

type bmad_normal_form_struct
  type (ele_struct), pointer :: ele_origin => null()  ! Element at which the on-turn map was created.
  type (taylor_struct) :: M(6) = taylor_struct()                  ! One-turn taylor map: M = A o N o A_inv, N = exp(:h:)
  type (taylor_struct) :: A(6) = taylor_struct()                  ! Map from Floquet -> Lab coordinates
  type (taylor_struct) :: A_inv(6) = taylor_struct()              ! Map from Lab -> Floquet coordinates
  type (taylor_struct) :: dhdj(6) = taylor_struct()               ! Nonlinear tune function operating on Floquet coordinates
  type (complex_taylor_struct) :: F(6) = complex_taylor_struct()  ! Vector field factorization in phasor basis:
  type (complex_taylor_struct) :: L(6) = complex_taylor_struct()  ! L component
  type (resonance_h_struct), allocatable :: h(:)
end type

type ptc_normal_form_struct
  type (ele_struct), pointer :: ele_origin => null()  ! Element at which the on-turn map was created.
  type (probe_8) one_turn_map                ! One turn map
  real(rp) orb0(6)                           ! Closed orbit at element.
  type (c_normal_form) normal_form           ! Complex normal form
  type (c_taylor) phase(3)                   ! Phase/chromaticity maps
  type (c_taylor) path_length                ! Path length map. Gives momentum compaction.
  type (c_taylor) spin_tune                  ! Amplitude dependent spin tune 
  type (c_quaternion) isf                    ! Invariant spin field in (x, px, ...) space.
  type (internal_state) state                ! PTC state
  logical :: valid_map = .false.
end type

!

type branch_struct
  character(40) :: name = ''       ! Name of line that defines the branch.
  integer :: ix_branch = -1        ! Index of this branch. 0 => Main branch
  integer :: ix_from_branch = -1   ! -1 => No creating fork element to this branch.
  integer :: ix_from_ele = -1      ! Index of creating fork element which forks to this branch.
  integer :: ix_to_ele = -1        ! Index of element in this branch that creating fork element forks to.
  integer :: n_ele_track
  integer :: n_ele_max
  type (lat_struct), pointer :: lat => null()
  type (mode_info_struct) :: a , b , z  ! Note: Tunes are the fractional part.
  type (ele_struct), pointer :: ele(:) => null()
  type (lat_param_struct) :: param 
  type (wall3d_struct), pointer :: wall3d(:) => null()
  type (ptc_branch1_struct) ptc              ! Pointer to layout. Note: ptc info not transferred with "branch1 = branch2" set.
end type

integer, parameter :: opal$ = 1, impactt$ = 2
character(16) :: pre_tracker_name(0:2) = ['NONE   ', 'OPAL   ', 'IMPACTT']

type pre_tracker_struct
  integer :: who = 0   ! Can be opal$, or impactt$
  integer :: ix_ele_start = 0
  integer :: ix_ele_end = 0
  character(400) :: input_file = ''
end type

! lat_struct
! Remember: If this struct is changed you have to modify (among other things):
!     Increase bmad_inc_version by 1.
!     read_digested_bmad_file
!     write_digested_bmad_file
!     transfer_lat_parameters
!     lat_equal_lat
! Rule: When lat2 = lat2, lat2%surface and lat1%surface will point to the same location.

type lat_struct
  character(200) :: use_name = ''                     ! Name of lat given by USE statement
  character(40) :: lattice = ''                       ! Lattice
  character(40) :: machine = ''                       ! Name of the machine the lattice is for ("LHC", etc).
  character(400) :: input_file_name = ''              ! Name of the lattice input file
  character(80) :: title = ''                         ! General title
  character(100), allocatable :: print_str(:)         ! Saved print statements.
  type (expression_atom_struct), allocatable :: constant(:)  ! Constants defined in the lattice
  type (mode_info_struct), pointer :: a => null(), b => null(), z => null() ! Tunes (fractional part), etc. 
  type (lat_param_struct), pointer :: param => null() ! Parameters
  type (bookkeeping_state_struct) lord_state          ! lord bookkeeping status.
  type (ele_struct) ele_init                          ! For use by any program
  type (ele_struct), pointer ::  ele(:) => null()     ! Array of elements [=> branch(0)].
  type (branch_struct), allocatable :: branch(:)      ! Branch(0:) array
  type (control_struct), allocatable :: control(:)    ! Control list
  type (coord_struct) particle_start                  ! Starting particle_coords.
  type (beam_init_struct) beam_init                   ! Beam initialization.
  type (pre_tracker_struct) pre_tracker               ! For OPAL/IMPACT-T
  type (nametable_struct) nametable                   ! For quick searching by element name.
  real(rp), allocatable :: custom(:)                  ! Custom attributes.
  integer :: version = -1                             ! Version number
  integer, pointer :: n_ele_track => null()           ! Number of lat elements to track through.
  integer, pointer :: n_ele_max => null()             ! Index of last valid element in %ele(:) array
  integer :: n_control_max = 0                        ! Last index used in control_array
  integer :: n_ic_max = 0                             ! Last index used in ic_array
  integer :: input_taylor_order = 0                   ! As set in the input file
  integer, allocatable :: ic(:)                       ! Index to %control(:) from slaves.
  integer :: photon_type = incoherent$                ! Or coherent$. For X-ray simulations.
  integer :: creation_hash = 0                        ! Set by bmad_parser. creation_hash will vary if 
                                                      !   any of the lattice files are modified.
  integer :: ramper_slave_bookkeeping = stale$
end type

character(2), parameter :: coord_name(6) = ['x ', 'px', 'y ', 'py', 'z ', 'pz']
character(2), parameter :: coord_name_cap(6) = ['X ', 'Px', 'Y ', 'Py', 'Z ', 'Pz']

! KEY value definitions
! Note: sbend$ and rbend$ also used for sub_key

integer, parameter :: drift$ = 1, sbend$ = 2, quadrupole$ = 3, group$ = 4, sextupole$ = 5
integer, parameter :: overlay$ = 6, custom$ = 7, taylor$ = 8, rfcavity$ = 9, elseparator$ = 10
integer, parameter :: beambeam$ = 11, wiggler$ = 12, sol_quad$ = 13, marker$ = 14, kicker$ = 15
integer, parameter :: hybrid$ = 16, octupole$ = 17, rbend$ = 18, multipole$ = 19, def_bmad_com$ = 20
integer, parameter :: def_mad_beam$ = 21, ab_multipole$ = 22, solenoid$ = 23, patch$ = 24, lcavity$ = 25
integer, parameter :: def_parameter$ = 26, null_ele$ = 27, beginning_ele$ = 28, def_line$ = 29
integer, parameter :: match$ = 30, monitor$ = 31, instrument$ = 32, hkicker$ = 33, vkicker$ = 34
integer, parameter :: rcollimator$ = 35, ecollimator$ = 36, girder$ = 37, converter$ = 38
integer, parameter :: def_particle_start$ = 39, photon_fork$ = 40, fork$ = 41, mirror$ = 42, crystal$ = 43
integer, parameter :: pipe$ = 44, capillary$ = 45, multilayer_mirror$ = 46, e_gun$ = 47, em_field$ = 48
integer, parameter :: floor_shift$ = 49, fiducial$ = 50, undulator$ = 51, diffraction_plate$ = 52
integer, parameter :: photon_init$ = 53, sample$ = 54, detector$ = 55, sad_mult$ = 56, mask$ = 57
integer, parameter :: ac_kicker$ = 58, lens$ = 59, def_space_charge_com$ = 60, crab_cavity$ = 61
integer, parameter :: ramper$ = 62, def_ptc_com$ = 63, rf_bend$ = 64, gkicker$ = 65, foil$ = 66
integer, parameter :: thick_multipole$ = 67, pickup$ = 68, feedback$ = 69, n_key$ = 69

! A "!" as the first character is to prevent name matching by the key_name_to_key_index routine.

character(20), parameter :: key_name(n_key$) = [ &
    'Drift             ', 'SBend             ', 'Quadrupole        ', 'Group             ', 'Sextupole         ', &
    'Overlay           ', 'Custom            ', 'Taylor            ', 'RFCavity          ', 'ELSeparator       ', &
    'BeamBeam          ', 'Wiggler           ', 'Sol_Quad          ', 'Marker            ', 'Kicker            ', &
    'Hybrid            ', 'Octupole          ', 'RBend             ', 'Multipole         ', '!Bmad_Com         ', &
    '!Mad_Beam         ', 'AB_multipole      ', 'Solenoid          ', 'Patch             ', 'Lcavity           ', &
    '!Parameter        ', 'Null_Ele          ', 'Beginning_Ele     ', '!Line             ', 'Match             ', &
    'Monitor           ', 'Instrument        ', 'HKicker           ', 'VKicker           ', 'RCollimator       ', &
    'ECollimator       ', 'Girder            ', 'Converter         ', '!Particle_Start   ', 'Photon_Fork       ', &
    'Fork              ', 'Mirror            ', 'Crystal           ', 'Pipe              ', 'Capillary         ', &
    'Multilayer_Mirror ', 'E_Gun             ', 'EM_Field          ', 'Floor_Shift       ', 'Fiducial          ', &
    'Undulator         ', 'Diffraction_Plate ', 'Photon_Init       ', 'Sample            ', 'Detector          ', &
    'Sad_Mult          ', 'Mask              ', 'AC_Kicker         ', 'Lens              ', '!Space_Charge_Com ', &
    'Crab_Cavity       ', 'Ramper            ', '!PTC_Com          ', 'RF_Bend           ', 'GKicker           ', &
    'Foil              ', 'Thick_Multipole   ', 'Pickup            ', 'Feedback          ']

! These logical arrays get set in init_attribute_name_array and are used
! to sort elements that have kick or orientation attributes from elements that do not.
! The orientation attributes are: tilt, x/y/z_offset, x/y_pitch, and *_tot versions.
! Note: A solenoid does not formally have a tilt but has everything else.
! Rule: Any element that formally has some but not all orientation attributes is considered
!   internally to have all attributes so any routine can safely work with all the 
!   orientation attributes as a block.

logical has_hkick_attributes(n_key$)
logical has_kick_attributes(n_key$)

!

integer, parameter :: standard$ = 1, match_twiss$ = 2, identity$ = 3, phase_trombone$ = 4
integer, parameter :: match_orbit$ = 2, zero$ = 3
character(16) :: matrix_name(4) = [character(16):: 'standard', 'match_twiss', 'identity', 'phase_trombone']
character(12) :: kick0_name(3) = [character(16):: 'standard', 'match_orbit', 'zero']

! Element attribute name logical definitions

integer, parameter :: val1$=19, val2$=20, val3$=21, val4$=22, val5$=23, &
          val6$=24, val7$=25, val8$=26, val9$=27, val10$=28, val11$=29, &
          val12$=30

! Note: The following must not be varied without varying appropriate code:
!   Range [beta_a0$, alpha_b0$] and [beta_a1$, alpha_b1$] hold all twiss.
!   Range [eta_x0$, etap_y0$] and [eta_x1$, etap_y1$] hold all dispersion.
!   Range [c11_mat0$, c22_mat0$] and [c11_mat1$, c22_mat1$] hold all C-matrix values
    
integer, parameter :: beta_a0$ = 2, alpha_a0$ = 3, beta_b0$ = 4, alpha_b0$ = 5
integer, parameter :: beta_a1$ = 6, alpha_a1$ = 7, beta_b1$ = 8, alpha_b1$ = 9
integer, parameter :: dphi_a$ = 10, dphi_b$ = 11
integer, parameter :: eta_x0$ = 12, etap_x0$ = 13, eta_y0$ = 14, etap_y0$ = 15
integer, parameter :: eta_x1$ = 16, etap_x1$ = 17, eta_y1$ = 18, etap_y1$ = 19
integer, parameter :: c11_mat0$ = 20, c12_mat0$ = 21, c21_mat0$ = 22, c22_mat0$ = 23, mode_flip0$ = 24
integer, parameter :: c11_mat1$ = 25, c12_mat1$ = 26, c21_mat1$ = 27, c22_mat1$ = 28, mode_flip1$ = 29

integer, parameter :: x0$ = 30, px0$ = 31, y0$ = 32, py0$ = 33, z0$ = 34, pz0$ = 35
integer, parameter :: x1$ = 36, px1$ = 37, y1$ = 38, py1$ = 39, z1$ = 40, pz1$ = 41
integer, parameter :: matrix$ = 42, kick0$ = 43, recalc$ = 44
integer, parameter :: delta_time$ = 48

integer, parameter :: x$ = 1, px$ = 2, y$ = 3, py$ = 4, z$ = 5, pz$ = 6
integer, parameter :: t$ = 8
integer, parameter :: field_x$ = 10, field_y$ = 11, phase_x$ = 12, phase_y$ = 13
integer, parameter :: e_photon$ = 9

integer, parameter :: e1$ = 19, e2$ = 20
integer, parameter :: fint$ = 21, fintx$ = 22, hgap$ = 23, hgapx$ = 24, h1$ = 25, h2$ = 26

integer, parameter :: radius$ = 3, focal_strength$ = 5

integer, parameter :: l$ = 1                          ! Assumed unique. Do not assign 1 to another attribute.
integer, parameter :: tilt$ = 2, roll$ = 2, n_part$ = 2, inherit_from_fork$ = 2 ! Important: tilt$ = roll$
integer, parameter :: ref_tilt$ = 3, direction$ = 3, repetition_frequency$ = 3, deta_ds_master$ = 3, &
                      kick$ = 3, x_gain_err$ = 3, taylor_order$ = 3, r_solenoid$ = 3, final_charge$ = 3
integer, parameter :: k1$ = 4, kx$ = 4, harmon$ = 4, h_displace$ = 4, y_gain_err$ = 4, s_twiss_ref$ = 4, &
                      critical_angle_factor$ = 4, tilt_corr$ = 4, ref_coords$ = 4, dt_max$ = 4
integer, parameter :: graze_angle$ = 5, k2$ = 5, b_max$ = 5, v_displace$ = 5, gradient_tot$ = 5, harmon_master$ = 5, &
                      ks$ = 5, flexible$ = 5, crunch$ = 5, ref_orbit_follows$ = 5, pc_out_min$ = 5
integer, parameter :: gradient$ = 6, k3$ = 6, noise$ = 6, new_branch$ = 6, ix_branch$ = 6, g_max$ = 6, &
                      g$ = 6, symmetry$ = 6, field_scale_factor$ = 6, pc_out_max$ = 6
integer, parameter :: dg$ = 7, bbi_const$ = 7, osc_amplitude$ = 7, ix_to_branch$ = 7, angle_out_max$ = 7, &
                      gradient_err$ = 7, critical_angle$ = 7, bragg_angle_in$ = 7, spin_dn_dpz_x$ = 7
integer, parameter :: delta_e_ref$ = 8, interpolation$ = 8, bragg_angle_out$ = 8, k1x$ = 8, spin_dn_dpz_y$ = 8, &
                      charge$ = 8, x_gain_calib$ = 8, ix_to_element$ = 8, voltage$ = 8, g_tot$ = 8
integer, parameter :: rho$ = 9, voltage_err$ = 9, bragg_angle$ = 9, k1y$ = 9, n_particle$ = 9, spin_dn_dpz_z$ = 9
integer, parameter :: fringe_type$ = 10, dbragg_angle_de$ = 10
integer, parameter :: fringe_at$ = 11, gang$ = 11, darwin_width_sigma$ = 11
integer, parameter :: darwin_width_pi$ = 12
integer, parameter :: spin_fringe_on$ = 13, pendellosung_period_sigma$ = 13
integer, parameter :: sig_x$ = 14, exact_multipoles$ = 14, pendellosung_period_pi$ = 14
integer, parameter :: sig_y$ = 15, graze_angle_in$ = 15, r0_elec$ = 15, rf_frequency$ = 15
integer, parameter :: sig_z$ = 16, graze_angle_out$ = 16, r0_mag$ = 16, rf_wavelength$ = 16
integer, parameter :: sig_vx$ = 17, static_linear_map$ = 17
! longitudinal_mode$ is near to rf_wavelength$ for type_ele to print rf_bucket_length near rf_wavelength$
integer, parameter :: sig_vy$ = 18, constant_ref_energy$ = 18, longitudinal_mode$ = 18
integer, parameter :: sig_e$ = 19, sig_pz$ = 19, autoscale_amplitude$ = 19
integer, parameter :: d1_thickness$ = 20, default_tracking_species$ = 20, autoscale_phase$ = 20, &
                      n_slice$ = 20, y_gain_calib$ = 20, sig_e2$ = 20
integer, parameter :: fb1$ = 21, polarity$ = 21, crunch_calib$ = 21, alpha_angle$ = 21, d2_thickness$ = 21, &
                      beta_a_strong$ = 21, beta_a_out$ = 21, e_loss$ = 21, gap$ = 21, spin_x$ = 21, &
                      E_center$ = 21, scatter_test$ = 21
integer, parameter :: fb2$ = 22, x_offset_calib$ = 22, v1_unitcell$ = 22, psi_angle$ = 22, cavity_type$ = 22, &
                      beta_b_strong$ = 22, beta_b_out$ = 22, spin_y$ = 22, E2_center$ = 22, n_period$ = 22, &
                      emit_fraction$ = 22, x1_edge$ = 22
integer, parameter :: y_offset_calib$ = 23, v_unitcell$ = 23, v2_unitcell$ = 23, spin_z$ = 23, l_period$ = 23, &
                      fq1$ = 23, alpha_a_strong$ = 23, alpha_a_out$ = 23, E2_probability$ = 23, phi0_max$ = 23, &
                      x2_edge$ = 23
integer, parameter :: fq2$ = 24, phi0$ = 24, tilt_calib$ = 24, E_center_relative_to_ref$ = 24, y1_edge$ = 24, &
                      alpha_b_strong$ = 24, alpha_b_out$ = 24, is_mosaic$ = 24, px_aperture_width2$ = 24
integer, parameter :: phi0_err$ = 25, current$ = 25, mosaic_thickness$ = 25, px_aperture_center$ = 25, &
                      eta_x_out$ = 25, quad_tilt$ = 25, de_eta_meas$ = 25, spatial_distribution$ = 25, &
                      y2_edge$ = 25, species_strong$ = 25
integer, parameter :: eta_y_out$ = 26, mode$ = 26, velocity_distribution$ = 26, py_aperture_width2$ = 26, &
                      phi0_multipass$ = 26, n_sample$ = 26, origin_ele_ref_pt$ = 26, mosaic_angle_rms_in_plane$ = 26, &
                      eps_step_scale$ = 26, E_tot_strong$ = 26, dthickness_dx$ = 26, bend_tilt$ = 26
integer, parameter :: etap_x_out$ = 27, phi0_autoscale$ = 27, dx_origin$ = 27, energy_distribution$ = 27, &
                      x_quad$ = 27, ds_photon_slice$ = 27, mosaic_angle_rms_out_plane$ = 27, &
                      py_aperture_center$ = 27, x_dispersion_err$ = 27, l_rectangle$ = 27, pc_strong$ = 27
integer, parameter :: etap_y_out$ = 28, dy_origin$ = 28, y_quad$ = 28, e_field_x$ = 28, &
                      y_dispersion_err$ = 28, z_aperture_width2$ = 28, user_sets_length$ = 28, &
                      rf_clock_harmonic$ = 28, b_field_tot$ = 28
integer, parameter :: upstream_coord_dir$ = 29, dz_origin$ = 29, mosaic_diffraction_num$ = 29, &
                      cmat_11$ = 29, field_autoscale$ = 29, l_sagitta$ = 29, e_field_y$ = 29, &
                      x_dispersion_calib$ = 29, z_aperture_center$ = 29, f_factor$ = 29
integer, parameter :: cmat_12$ = 30, dtheta_origin$ = 30, b_param$ = 30, l_chord$ = 30, &
                      downstream_coord_dir$ = 30, pz_aperture_width2$ = 30, y_dispersion_calib$ = 30, &
                      scale_field_to_one$ = 30, voltage_tot$ = 30, scatter_method$ = 30
integer, parameter :: cmat_21$ = 31, l_active$ = 31, dphi_origin$ = 31, split_id$ = 31, ref_cap_gamma$ = 31, &
                      l_soft_edge$ = 31, transverse_sigma_cut$ = 31, pz_aperture_center$ = 31, &
                      mean_excitation_energy$ = 31, fiducial_pt$ = 31
integer, parameter :: cmat_22$ = 32, dpsi_origin$ = 32, t_offset$ = 32, ds_slice$ = 32, use_reflectivity_table$ = 32, init_needed$ = 32
integer, parameter :: angle$ = 33, n_cell$ = 33, mode_flip$ = 33, crossing_time$ = 33, x_kick$ = 33
integer, parameter :: x_pitch$ = 34, px_kick$ = 34   ! Note: [x_kick$, px_kick$, ..., pz_kick$] must be in order.
integer, parameter :: y_pitch$ = 35, y_kick$ = 35
integer, parameter :: x_offset$ = 36, py_kick$ = 36
integer, parameter :: y_offset$ = 37, z_kick$ = 37
integer, parameter :: z_offset$ = 38, pz_kick$ = 38
integer, parameter :: hkick$ = 39, d_spacing$ = 39, x_offset_mult$ = 39, emittance_a$ = 39, crab_x1$ = 39
integer, parameter :: vkick$ = 40, y_offset_mult$ = 40, p0c_ref_init$ = 40, emittance_b$ = 40, crab_x2$ = 40
integer, parameter :: BL_hkick$ = 41, e_tot_ref_init$ = 41, emittance_z$ = 41, crab_x3$ = 41
integer, parameter :: BL_vkick$ = 42, crab_tilt$ = 42
integer, parameter :: BL_kick$ = 43, B_field$ = 43, E_field$ = 43, high_energy_space_charge_on$ = 43, crab_x4$=43
integer, parameter :: photon_type$ = 44, coupler_phase$ = 44, dB_field$ = 44, crab_x5$=44
integer, parameter :: lattice_type$ = 45, B1_gradient$ = 45, E1_gradient$ = 45, coupler_angle$ = 45
integer, parameter :: live_branch$ = 46, B2_gradient$ = 46, E2_gradient$ = 46, coupler_strength$ = 46
integer, parameter :: geometry$ = 47, coupler_at$ = 47, E_tot_offset$ = 47, ptc_canonical_coords$ = 47
integer, parameter :: B3_gradient$ = 48, E3_gradient$ = 48, ptc_fringe_geometry$ = 48, e_tot_set$ = 48
integer, parameter :: Bs_field$ = 49, p0c_set$ = 49, ptc_field_geometry$ = 49, delta_ref_time_user_set$ = 49
integer, parameter :: delta_ref_time$ = 50
integer, parameter :: p0c_start$ = 51
integer, parameter :: e_tot_start$ = 52
integer, parameter :: p0c$ = 53
integer, parameter :: e_tot$ = 54
integer, parameter :: x_pitch_tot$ = 55, no_end_marker$ = 55
integer, parameter :: y_pitch_tot$ = 56
integer, parameter :: x_offset_tot$ = 57
integer, parameter :: y_offset_tot$ = 58
integer, parameter :: z_offset_tot$ = 59
integer, parameter :: tilt_tot$ = 60, roll_tot$ = 60  ! Important: tilt_tot$ = roll_tot$
integer, parameter :: ref_tilt_tot$ = 61
integer, parameter :: multipass_ref_energy$ = 62
integer, parameter :: dispatch$ = 63
integer, parameter :: ref_time_start$ = 64
integer, parameter :: thickness$ = 65, integrator_order$ = 65   ! For Etiennes' PTC: 2, 4, 6, or 8.
integer, parameter :: num_steps$ = 66   ! Assumed unique by set_flags_for_changed_real_attribute
integer, parameter :: ds_step$ = 67     ! Assumed unique by set_flags_for_changed_real_attribute
integer, parameter :: csr_ds_step$ = 68
integer, parameter :: lord_pad1$ = 69
integer, parameter :: lord_pad2$ = 70, ref_wavelength$ = 70
integer, parameter :: x1_limit$ = 71
integer, parameter :: x2_limit$ = 72
integer, parameter :: y1_limit$ = 73
integer, parameter :: y2_limit$ = 74
integer, parameter :: check_sum$ = 75

!!    = 1 + num_ele_attrib$

integer, parameter :: distribution$ = 81
integer, parameter :: tt$ = 81, x_knot$ = 81
integer, parameter :: alias$  = 82, max_fringe_order$ = 82, eta_x$ = 82
integer, parameter :: electric_dipole_moment$ = 83, lr_self_wake_on$ = 83, x_ref$ = 83, species_out$ = 83
integer, parameter :: y_knot$ = 83, eta_y$ = 83, density$ = 83
integer, parameter :: lr_wake_file$ = 84, px_ref$ = 84, etap_x$ = 84, slave$ = 84, &
                      density_used$ = 84
integer, parameter :: lr_freq_spread$ = 85, y_ref$ = 85, etap_y$ = 85, &
                      area_density$ = 85, input_ele$ = 85
integer, parameter :: lattice$ = 86, phi_a$ = 86, multipoles_on$ = 86, py_ref$ = 86, &
                      area_density_used$ = 86, output_ele$ = 86
integer, parameter :: aperture_type$ = 87, eta_z$ = 87, machine$ = 87
integer, parameter :: taylor_map_includes_offsets$ = 88, pixel$ = 88, p88$ = 88, radiation_length$ = 88, deta_dpz_x$ = 88
integer, parameter :: csr_method$ = 89, var$ = 89, z_ref$ = 89, p89$ = 89, radiation_length_used$ = 89, deta_dpz_y$ = 89

integer, parameter :: pz_ref$ = 90, space_charge_method$ = 90, p90$ = 90, detap_dpz_x$ = 90
integer, parameter :: mat6_calc_method$ = 91, detap_dpz_y$ = 91
integer, parameter :: tracking_method$  = 92, s_long$ = 92
integer, parameter :: ref_time$ = 93, ptc_integration_type$ = 93
integer, parameter :: spin_tracking_method$ = 94, eta_a$ = 94
integer, parameter :: aperture$ = 95, etap_a$ = 95
integer, parameter :: x_limit$ = 96, absolute_time_tracking$ = 96, eta_b$ = 96
integer, parameter :: y_limit$ = 97, etap_b$ = 97
integer, parameter :: offset_moves_aperture$ = 98
integer, parameter :: aperture_limit_on$ = 99, alpha_a$ = 99, reflectivity_table$ = 99, energy_probability_curve$ = 99

integer, parameter :: exact_misalign$ = 100, physical_source$ = 100
integer, parameter :: sr_wake_file$ = 100, alpha_b$ = 100
integer, parameter :: term$ = 101, frequencies$ = 101, old_integrator$ = 101, curvature$ = 101
integer, parameter :: x_position$ = 102, exact_model$ = 102
integer, parameter :: symplectify$ = 103, y_position$ = 103, n_slice_spline$ = 103
integer, parameter :: z_position$ = 104, amp_vs_time$ = 104
integer, parameter :: is_on$ = 105, theta_position$ = 105, vertical_kick$ = 105
integer, parameter :: field_calc$ = 106, phi_position$ = 106
integer, parameter :: psi_position$ = 107, wall$ = 107
integer, parameter :: aperture_at$ = 108, beta_a$ = 108
integer, parameter :: ran_seed$ = 109, origin_ele$ = 109, beta_b$ = 109

! 

integer, parameter :: to_line$ = 110, field_overlaps$ = 110, dbeta_dpz_a$ = 110
integer, parameter :: field_master$ = 111, to_element$ = 111, dbeta_dpz_b$ = 111
integer, parameter :: descrip$ = 112
integer, parameter :: scale_multipoles$ = 113, dalpha_dpz_a$ = 113
integer, parameter :: sr_wake$ = 114, dalpha_dpz_b$ = 114
integer, parameter :: ref_orbit$ = 115, lr_wake$ = 115
integer, parameter :: phi_b$ = 116, crystal_type$ = 116, material_type$ = 116
integer, parameter :: type$ = 117
integer, parameter :: ref_origin$ = 118
integer, parameter :: ele_origin$ = 119

integer, parameter :: superimpose$     = 120   
integer, parameter :: super_offset$    = 121
integer, parameter :: reference$       = 122
integer, parameter :: cartesian_map$   = 123
integer, parameter :: cylindrical_map$ = 124
integer, parameter :: grid_field$      = 125
integer, parameter :: gen_grad_map$    = 126
integer, parameter :: create_jumbo_slave$ = 127

integer, parameter :: accordion_edge$  = 128
integer, parameter :: start_edge$  = 129
integer, parameter :: end_edge$  = 130
integer, parameter :: s_position$ = 131
integer, parameter :: ref_species$ = 132, particle$ = 132
integer, parameter :: wrap_superimpose$ = 133

integer, parameter :: a0$  = 140, a21$  = 161
integer, parameter :: b0$  = 162, b21$  = 183

integer, parameter :: k0l$ = 140, k21l$ = 161
integer, parameter :: t0$  = 162, t21$  = 183

integer, parameter :: k0sl$ = 190, k21sl$ = 211

integer, parameter :: a0_elec$ = 190, a21_elec$ = 211
integer, parameter :: b0_elec$ = 212, b21_elec$ = 233

integer, parameter :: custom_attribute0$ = b21_elec$
integer, parameter :: custom_attribute_num$ = 40
integer, parameter :: num_ele_attrib_extended$ = custom_attribute0$ + custom_attribute_num$

integer, parameter :: g_err$ = dg$   ! For backwards compatibility.
integer, parameter :: B_field_err$ = dB_field$  ! For backwards compatibility

character(40), parameter :: blank_name$ = ' '

! lattice logical names

integer, parameter :: open$ = 1, closed$ = 2

character(16), parameter :: lattice_type_name(0:2) = ['GARBAGE!        ', 'Linear_Lattice  ', 'Circular_Lattice']
character(16), parameter :: geometry_name(0:2) = ['GARBAGE!    ', 'Open        ', 'Closed      ']

! The linac_normal_mode_struct is basically the synchrotron integrals with the energy factors thrown in.
! Note: The b%emittance calc from radiation integrals will include the photon vertical opening angle in the calc.

type anormal_mode_struct
  real(rp) :: emittance = 0         ! Beam emittance (unnormalized). Includes vertical photon opening angle.
  real(rp) :: emittance_no_vert = 0 ! Unnormalized beam emittance without the vertical photon opening angle taken into account.
  real(rp) :: synch_int(4:6) = 0    ! Synchrotron integrals
  real(rp) :: j_damp = 0            ! damping partition number
  real(rp) :: alpha_damp = 0        ! damping per turn
  real(rp) :: chrom = 0             ! Chromaticity
  real(rp) :: tune = 0              ! "Fractional" tune in radians
end type

type linac_normal_mode_struct
  real(rp) :: i2_E4 = 0           ! Integral: g^2 * gamma^4
  real(rp) :: i3_E7 = 0           ! Integral: g^3 * gamma^7
  real(rp) :: i5a_E6 = 0          ! Integral: (g^3 * H_a) * gamma^6
  real(rp) :: i5b_E6 = 0          ! Integral: (g^3 * H_b) * gamma^6
  real(rp) :: sig_E1 = 0          ! Energy spread after 1 pass (eV)
  real(rp) :: a_emittance_end = 0 ! a mode emittance at end of linac
  real(rp) :: b_emittance_end = 0 ! b mode emittance at end of linac
end type

type normal_modes_struct
  real(rp) :: synch_int(0:3) = 0  ! Synchrotron integrals I0, I1, I2, and I3
  real(rp) :: sigE_E = 0          ! SigmaE/E
  real(rp) :: sig_z = 0           ! Sigma_Z
  real(rp) :: e_loss = 0          ! Energy loss / turn (eV)
  real(rp) :: rf_voltage = 0      ! Total rfcavity voltage (eV)
  real(rp) :: pz_aperture = 0     ! pz aperture limit. Used with Touschek calculations.
  real(rp) :: pz_average = 0      ! Average over branch due to damping.
  real(rp) :: momentum_compaction = 0
  real(rp) :: dpz_damp = 0        ! Change in pz without RF
  type (anormal_mode_struct) :: a = anormal_mode_struct()
  type (anormal_mode_struct) :: b = anormal_mode_struct()
  type (anormal_mode_struct) :: z = anormal_mode_struct()
  type (linac_normal_mode_struct) :: lin = linac_normal_mode_struct()
end type

integer, parameter :: bends$ = 201
integer, parameter :: wigglers$ = 202
integer, parameter :: all$ = 203
integer, parameter :: upstream$ = 1, downstream$ = 2

!---------------------------------------------------------------------------
! Units

integer, parameter :: radians$ = 1, degrees$ = 2, cycles$ = 3, radians_over_2pi$ = 3
character(8), parameter :: angle_units_name(4) = [character(8):: 'radians', 'degrees', 'cycles', 'radians_over_2pi']
character(8), parameter :: short_angle_units_name(4) = [character(8):: 'rad', 'deg', 'rad/2pi', 'rad/2pi']
real(rp) :: radians_to_angle_units(4) = [1.0_rp, 180/pi, 1/twopi, 1/twopi]

! Electric and magnetic fields.

type em_field_struct
  real(rp) :: E(3) = 0        ! electric field.
  real(rp) :: B(3) = 0        ! magnetic field.
  real(rp) :: dE(3,3) = 0     ! electric field gradient.
  real(rp) :: dB(3,3) = 0     ! magnetic field gradient.
  real(rp) :: phi = 0         ! Electric scalar potential.
  real(rp) :: phi_B = 0       ! Magnetic scalar potential.
  real(rp) :: A(3) = 0        ! Magnetic vector potential.
end type

! Grid of grid_field information
integer, parameter :: rotationally_symmetric_rz$ = 1, xyz$ = 2
character(30), parameter :: grid_field_geometry_name(0:2) = &
                          [character(30) :: 'GARBAGE!', 'rotationally_symmetric_rz', 'xyz']
integer, parameter :: grid_field_dimension(2) = [2, 3] 

! Structures for saving the track through an element.

type strong_beam_struct
  integer :: ix_slice = 0                   ! 0 -> at element center and not at slice.
  real(rp) :: x_center = 0, y_center = 0    ! Strong beam slice center.
  real(rp) :: x_sigma = 0, y_sigma = 0      ! Strong beam slice sigma.
  real(rp) :: dx = 0, dy = 0                ! Particle - beam slice distance.
end type

type track_point_struct
  real(rp) s_lab                                  ! Longitudinal lab coord with respect to the upstream end.
  real(rp) s_body                                 ! Longitudinal body coord with respect to the entrance end.
  type (coord_struct) orb                         ! An array of track points indexed from 0 (%orb(0:)).
  type (em_field_struct) field                    ! An array of em fields indexed from 0 (%field(0:)).
  type (strong_beam_struct) strong_beam           ! Strong beam info for beambeam element.
  real(rp) vec0(6)                                ! 0th order part of xfer map from the beginning.
  real(rp) mat6(6,6)                              ! 1st order part of xfer map (transfer matrix).
end type

! Valid track%orb(:) points in range 0:track%n_pt

type track_struct
  type (track_point_struct), allocatable :: pt(:) ! Array of track points indexed from 0.
  real(rp) :: ds_save = 1d-3  ! Min distance between points. Not positive => Save at all points.
  integer :: n_pt = -1        ! Track upper bound for %pt(0:) array.
  ! n_bad and n_ok are used by adaptive trackers to record the number of times the step length had to be shortened.
  integer :: n_bad = 0        ! Number of "bad" steps where the step length was shortened.
  integer :: n_ok = 0         ! Number of "good" steps where the step length was not shortened.
end type

!------------------------------------------------------------------------------
! Multipass structures
! Multipass_lord_info_struct gives complete information about a single
! multipass_lord and all its slaves ("from the top down").
! If the multipass_lord has super_lords as slaves, n_super will be the number of
! super_slaves per super_lord.
! slave(1:n_pass, 1:n_super_slave) is a matrix of slaves in the tracking lattice.
! If there are no super_lords then n_super_slave = 1 and
!      super_lord(1:n_pass) = slave (1:n_pass, 1)

type multipass_lord_info_struct
  type (ele_struct), pointer :: lord          ! Lord element
  integer n_pass           ! Number of passes (= number of slaves)
  integer n_super_slave    ! Number of super_slaves per super_lord.
  type (ele_pointer_struct), allocatable :: super_lord(:)  ! Super_lord list if they exist.
  type (ele_pointer_struct), allocatable :: slave(:,:)     ! Slaves list in tracking part.
end type

! Multipass_ele_info_struct gives information about a singe element in the lattice
! ("from the bottom up").

type multipass_ele_info_struct
  logical multipass     ! True if involved in multipass. False otherwise
  integer ix_pass       ! Pass number
  integer, allocatable :: ix_lord(:)   ! Pointers to lord(:) array
  integer, allocatable :: ix_super(:) ! Indexes to slave(ix_pass, super_slave%ix_ele) matrix
end type

type multipass_branch_info_struct
  type (multipass_ele_info_struct), allocatable :: ele(:)
end type

! %lord(i), i = 1, ..., n = number of multipass_lords in the lattice.
! %branch(i), i = 0, ..., n = ubound(lat%branch)

type multipass_all_info_struct
  type (multipass_lord_info_struct), allocatable :: lord(:)      ! Array of lords
  type (multipass_branch_info_struct), allocatable :: branch(:)
end type

! Dynamic aperture structures...

! Dynamic aperture data point structure

type aperture_point_struct
  real(rp) x, y     ! (x,y) aperture point with respect to the reference orbit.
  integer plane     ! plane determining loss
  integer ix_ele    ! ele index particle lost at
  integer i_turn    ! turn particle lost at
end type

! Input parameters for a dynamic aperture scan

type aperture_param_struct
  real(rp) :: min_angle = 0
  real(rp) :: max_angle = pi
  integer :: n_angle   = 9
  integer :: n_turn = 100             ! Number of turns a particle must survive.
  real(rp) :: x_init = 1e-3_rp        ! Initial x coordinate to start with for theta_xy = 0.
  real(rp) :: y_init = 1e-3_rp        ! Initial y coordinate to start with for theta_xy = pi/2.
  real(rp) :: rel_accuracy = 1e-2_rp  ! Relative resolution of bracketed aperture.
  real(rp) :: abs_accuracy = 1e-5_rp  ! Absolute resolution of bracketed aperture (meters).
  character(40) :: start_ele = ''     ! Element to start tracking at.
end type

! Structure for a single dynamic aperture scan over a set of angles.

type aperture_scan_struct
  type (aperture_point_struct), allocatable :: point(:) ! Set of aperture points at different angles.
  type (coord_struct) :: ref_orb             ! Ref orbit around which the scan is made.
  real(rp) :: pz_start = 0                   ! Starting pz.
end type

!-------------------------------------------------------------------------
! Space charge parameters

type space_charge_common_struct                   ! Common block for space charge
  real(rp) :: ds_track_step = 0                   ! CSR tracking step size
  real(rp) :: dt_track_step = 1d-12               ! Time Runge kutta initial step.
  real(rp) :: cathode_strength_cutoff = 0.01      ! Cutoff for the cathode field calc.
  real(rp) :: rel_tol_tracking = 1d-8             ! Relative tolerance for tracking.
  real(rp) :: abs_tol_tracking = 1d-10            ! Absolute tolerance for tracking.
  real(rp) :: beam_chamber_height = 0             ! Used in shielding calculation.
  real(rp) :: lsc_sigma_cutoff = 0.1              ! Cutoff for the 1-dim longitudinal SC calc.
                                                  !   If a bin sigma is < cutoff * sigma_ave then ignore.
  real(rp) :: particle_sigma_cutoff = -1          ! 3D SC calc cutoff for particles with (x,y,z) position far from the center.
                                                  !  Negative or zero means ignore.
  integer :: space_charge_mesh_size(3) = [32, 32, 64]  ! Gird size for fft_3d space charge calc.
  integer :: csr3d_mesh_size(3) = [32, 32, 64]         ! Gird size for CSR.
  integer :: n_bin = 0                            ! Number of bins used
  integer :: particle_bin_span = 2                ! Longitudinal particle length / dz_bin
  integer :: n_shield_images = 0                  ! Chamber wall shielding. 0 = no shielding.
  integer :: sc_min_in_bin = 10                   ! Minimum number of particles in a bin for sigmas to be valid.
  logical :: lsc_kick_transverse_dependence = .false.
  logical :: debug = .false.
  character(400) :: diagnostic_output_file = ''   ! If non-blank write a diagnostic (EG wake) file
end type

type (space_charge_common_struct), save, target :: space_charge_com

!------------------------------------------------------------------------------

type time_runge_kutta_common_struct
  integer :: num_steps_done = -1              ! Number of integration steps. Not used by Bmad. For external use.
  logical :: print_too_many_step_err = .true.
end type

type (time_runge_kutta_common_struct), save :: time_runge_kutta_com

!------------------------------------------------------------------------------

integer, parameter :: invalid_name$ = 0, is_logical$ = 1, is_integer$ = 2, is_real$ = 3, is_switch$ = 4, is_string$ = 5
integer, parameter :: is_struct$ = 6, unknown$ = 7

! For coords_floor_to_curvilinear status argument

integer, parameter :: patch_problem$ = 2, outside$ = 3, cannot_find$ = 4

! extra_parsing_info_struct is used by parsing routines.
! %undeterministic_ran_function_called: Only set True when a ran function is called with ran_seed = 0

type extra_parsing_info_struct
  type (random_state_struct) :: ran_state           = random_state_struct()
  integer :: ran_seed                               = 0   
  logical :: undeterministic_ran_function_called    = .false.
  ! Used with bmad_com
  logical :: d_orb_set                              = .false.
  logical :: max_aperture_limit_set                 = .false.
  logical :: default_ds_step_set                    = .false.
  logical :: significant_length_set                 = .false.
  logical :: rel_tol_tracking_set                   = .false.
  logical :: abs_tol_tracking_set                   = .false.
  logical :: rel_tol_adaptive_tracking_set          = .false.
  logical :: abs_tol_adaptive_tracking_set          = .false.
  logical :: init_ds_adaptive_tracking_set          = .false.
  logical :: min_ds_adaptive_tracking_set           = .false.
  logical :: fatal_ds_adaptive_tracking_set         = .false.
  logical :: synch_rad_scale_set                    = .false.
  logical :: autoscale_amp_abs_tol_set              = .false.
  logical :: autoscale_amp_rel_tol_set              = .false.
  logical :: autoscale_phase_tol_set                = .false.
  logical :: rf_phase_below_transition_ref_set      = .false.
  logical :: electric_dipole_moment_set             = .false.
  logical :: taylor_order_set                       = .false.
  logical :: runge_kutta_order_set                  = .false.
  logical :: default_integ_order_set                = .false.
  logical :: sr_wakes_on_set                        = .false.
  logical :: lr_wakes_on_set                        = .false.
  logical :: high_energy_space_charge_on_set        = .false.
  logical :: csr_and_space_charge_on_set            = .false.
  logical :: spin_tracking_on_set                   = .false.
  logical :: spin_sokolov_ternov_flipping_on_set    = .false.
  logical :: radiation_damping_on_set               = .false.
  logical :: radiation_zero_average_set             = .false.
  logical :: radiation_fluctuations_on_set          = .false.
  logical :: conserve_taylor_maps_set               = .false.
  logical :: absolute_time_tracking_set             = .false.
  logical :: absolute_time_ref_shift_set            = .false.
  logical :: convert_to_kinetic_momentum_set        = .false.
  logical :: aperture_limit_on_set                  = .false.
  logical :: sad_eps_scale_set                      = .false.
  logical :: sad_amp_max_set                        = .false.
  logical :: sad_n_div_max_set                      = .false.
  logical :: max_num_runge_kutta_step_set           = .false.
  logical :: debug_set                              = .false.
  ! Used with space_charge_com
  logical :: ds_track_step_set                      = .false.
  logical :: dt_track_step_set                      = .false.
  logical :: cathode_strength_cutoff_set            = .false.
  logical :: sc_rel_tol_tracking_set                = .false.  ! For: space_charge_com%rel_tol_tracking
  logical :: sc_abs_tol_tracking_set                = .false.  ! For: space_charge_com%abs_tol_tracking
  logical :: beam_chamber_height_set                = .false.
  logical :: lsc_sigma_cutoff_set                   = .false.
  logical :: particle_sigma_cutoff_set              = .false.
  logical :: space_charge_mesh_size_set             = .false.
  logical :: csr3d_mesh_size_set                    = .false.
  logical :: n_bin_set                              = .false.
  logical :: particle_bin_span_set                  = .false.
  logical :: n_shield_images_set                    = .false.
  logical :: sc_min_in_bin_set                      = .false.
  logical :: lsc_kick_transverse_dependence_set     = .false.
  logical :: sc_debug_set                           = .false.
  logical :: diagnostic_output_file_set             = .false.
  ! Used with ptc_com
  logical :: old_integrator_set                     = .false.
  logical :: use_orientation_patches_set            = .false.
  logical :: print_info_messages_set                = .false.
  logical :: max_fringe_order_set                   = .false.
  logical :: exact_model_set                        = .false.
  logical :: exact_misalign_set                     = .false.
  logical :: vertical_kick_set                      = .false.
  logical :: cut_factor_set                         = .false.
  logical :: translate_patch_drift_time_set         = .false.
end type

!------------------------------------------------------------------------------
! common stuff

! %max_aperture_limit is used when no limit is specified or when 
!   lat%param%aperture_limit_on = False.
! Remember: Change extra_parsing_info_struct if bmad_common_struct changed.
! abs_tol_tracking was changed from 1d-10 to 1d-11 since Tao needs the tighter tol for
! calculating derivatives. DCS 8/2020.

type bmad_common_struct
  real(rp) :: max_aperture_limit = 1d3                 ! Max Aperture.
  real(rp) :: d_orb(6)           = 1d-5                ! Orbit deltas for the mat6 via tracking calc.
  real(rp) :: default_ds_step    = 0.2_rp              ! Default integration step for eles without an explicit step calc.
  real(rp) :: significant_length = 1d-10               ! meter 
  real(rp) :: rel_tol_tracking = 1d-9                  ! Closed orbit relative tolerance.
  real(rp) :: abs_tol_tracking = 1d-12                 ! Closed orbit absolute tolerance.
  real(rp) :: rel_tol_adaptive_tracking = 1d-8         ! Runge-Kutta tracking relative tolerance.
  real(rp) :: abs_tol_adaptive_tracking = 1d-10        ! Runge-Kutta tracking absolute tolerance.
  real(rp) :: init_ds_adaptive_tracking = 1d-3         ! Initial step size
  real(rp) :: min_ds_adaptive_tracking = 0             ! Min step size to take.
  real(rp) :: fatal_ds_adaptive_tracking = 1d-8        ! If actual step size is below this particle is lost.
  real(rp) :: autoscale_amp_abs_tol = 0.1_rp           ! Autoscale absolute amplitude tolerance (eV).
  real(rp) :: autoscale_amp_rel_tol = 1d-6             ! Autoscale relative amplitude tolerance
  real(rp) :: autoscale_phase_tol = 1d-5               ! Autoscale phase tolerance.
  real(rp) :: electric_dipole_moment = 0               ! Particle's EDM. Call set_ptc to transfer value to PTC.
  real(rp) :: synch_rad_scale = 1.0_rp                 ! Synch radiation kick scale. 1 => normal, 0 => no kicks.
  real(rp) :: sad_eps_scale = 5.0d-3                   ! Used in sad_mult step length calc.
  real(rp) :: sad_amp_max = 5.0d-2                     ! Used in sad_mult step length calc.
  integer :: sad_n_div_max = 1000                      ! Used in sad_mult step length calc.
  integer :: taylor_order = 0                          ! Taylor order to use. 0 -> default = ptc_private%taylor_order_saved.
  integer :: runge_kutta_order = 4                     ! Runge Kutta order.
  integer :: default_integ_order = 2                   ! PTC integration order. 
  integer :: max_num_runge_kutta_step = 10000          ! Maximum number of RK steps before particle is considered lost.
  logical :: rf_phase_below_transition_ref = .false.   ! Autoscale uses below transition stable point for RFCavities?
  logical :: sr_wakes_on = .true.                      ! Short range wakefields?
  logical :: lr_wakes_on = .true.                      ! Long range wakefields
  logical :: auto_bookkeeper = .true.                  ! Automatic bookkeeping?
  logical :: high_energy_space_charge_on = .false.     ! High energy space charge effect switch.
  logical :: csr_and_space_charge_on = .false.         ! Space charge switch.
  logical :: spin_tracking_on = .false.                ! spin tracking?
  logical :: spin_sokolov_ternov_flipping_on = .false. ! Spin flipping during synchrotron radiation emission?
  logical :: radiation_damping_on = .false.            ! Radiation damping toggle.
  logical :: radiation_zero_average = .false.          ! Shift damping to be zero on the zero orbit to get rid of sawtooth?
  logical :: radiation_fluctuations_on = .false.       ! Radiation fluctuations toggle.
  logical :: conserve_taylor_maps = .true.             ! Enable bookkeeper to set ele%taylor_map_includes_offsets = F?
  logical :: absolute_time_tracking = .false.          ! Absolute or relative time tracking?
  logical :: absolute_time_ref_shift = .true.          ! Apply reference time shift when using absolute time tracking?
  logical :: convert_to_kinetic_momentum = .false.     ! Cancel kicks due to finite vector potential when doing symplectic tracking?
                                                       !   Set to True to test symp_lie_bmad against runge_kutta.
  logical :: aperture_limit_on = .true.                ! use apertures in tracking?
  logical :: debug = .false.                           ! Used for code debugging.
end type
  
type (bmad_common_struct), save, target :: bmad_com

! Bmad global private structure
! For communication between Bmad routines and Bmad based programs.

type bmad_private_struct
  logical :: random_on = .true.       ! Temporarily turned off, for example, with the closed orbit calc.
end type

type (bmad_private_struct), save, target :: bmad_private

! ptc_com common block.
! Setup in: set_ptc_com_pointers

type ptc_common_struct
  integer, pointer :: max_fringe_order  => null()  ! Points to PTC HIGHEST_FRINGE. 2 (default) => Quadrupole.
  integer, pointer :: old_integrator    => null()  ! Points to PTC OLD_INTEGRATOR. -1 = False, 1 = True.
  logical, pointer :: exact_model       => null()  ! Points to PTC EXACT_MODEL. Default True.
  logical, pointer :: exact_misalign    => null()  ! Points to PTC ALWAYS_EXACTMIS. Default True. Notice different names.
  real(rp), pointer :: vertical_kick    => null()  ! Points to PTC VERTICAL_KICK for 6D emittance calc. 0 => off, 1 => on (default).
  real(rp) :: cut_factor = 0.006                   ! Cut factor for PTC tracking
  logical :: print_step_warning = .false.          ! Print warning if element uses too many steps.
  ! Below is stuff that should not be set except by experts
  logical :: use_orientation_patches = .true.      ! offset, pitch, and tilt attributes are put in ptc patch?
  logical :: print_info_messages = .false.         ! Allow PTC to print informational messages (which can clutter the output)?
  logical :: translate_patch_drift_time = .true.   ! When a Bmad patch is translated to a PTC fibre, is the drift
                                                   !   time included in the translation?
end type

type (ptc_common_struct), save, target :: ptc_com, ptc_com_default

! PTC related parameters that should not be set by non-Bmad routines.
!
! %taylor_order_saved is what is used if bmad_com%taylor_order is not set ( =  0).
! %taylor_order_saved is initially 3.
! When parsing a lattice file, %taylor_order_saved will be set to the taylor order of the lattice.
!
! %taylor_order_ptc is used to keep track as to what has been actually set for PTC. 
! If different from bmad_com[taylor_order], and bmad_com[taylor_order] is nonzero, the Bmad bookkeeping will 
! set PTC's order to bmad_com[taylor_order] and then set %taylor_order_ptc to bmad_com[taylor_order].

type ptc_private_struct
  type (internal_state) :: base_state   ! Base PTC state. 
  real(rp) :: e_tot_set = 0
  integer :: taylor_order_ptc = 0       ! What has been set in PTC. 0 -> not yet set. See above.
  integer :: taylor_order_saved = 3     ! Default to use at startup.
  logical :: init_ptc_needed = .true.
  logical :: init_spin_needed = .true.
end type

type (ptc_private_struct), save, target :: ptc_private

real(rp), parameter :: small_rel_change$ = 1d-14

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
  real(rp) :: lin_sig_E       = 0 ! Running sum
  real(rp) :: n_steps         = 0 ! number of qromb steps needed
end type

! Structures for a lattice of rad_int1_structs

type rad_int_branch_struct
  type (rad_int1_struct), allocatable :: ele(:) ! Array is indexed from 0
end type

type rad_int_all_ele_struct
  type (rad_int_branch_struct), allocatable :: branch(:) ! Array is indexed from 0
end type

! openPMD Header information

type pmd_header_struct
  character(:), allocatable :: openPMD
  character(:), allocatable :: openPMDextension
  character(:), allocatable :: basePath
  character(:), allocatable :: particlesPath
  character(:), allocatable :: meshesPath
  character(:), allocatable :: author
  character(:), allocatable :: software
  character(:), allocatable :: softwareVersion
  character(:), allocatable :: date
  character(:), allocatable :: latticeFile
  character(:), allocatable :: latticeName
end type

!-------------------------------------------------------------------------------------
! Parameters for expression_mod

! The numeric$ category is for numeric constants [EG: "1.3d-5"].
! The constant$ category is for constants like pi.
! The variable$ category includes symbolic constants defined in a lattice file, lattice parameters, etc.
! The species$ category is for the species() function. 
! The species_const$ category is for particle species ('He3', etc).

integer, parameter :: end_stack$ = 0, plus$ = 1, minus$ = 2, times$ = 3, divide$ = 4
integer, parameter :: l_parens$ = 5, r_parens$ = 6, power$ = 7
integer, parameter :: unary_minus$ = 8, unary_plus$ = 9, no_delim$ = 10
integer, parameter :: sin$ = 11, cos$ = 12, tan$ = 13
integer, parameter :: asin$ = 14, acos$ = 15, atan$ = 16, abs$ = 17, sqrt$ = 18
integer, parameter :: log$ = 19, exp$ = 20, ran$ = 21, ran_gauss$ = 22, atan2$ = 23
integer, parameter :: factorial$ = 24, int$ = 25, nint$ = 26, floor$ = 27, ceiling$ = 28
integer, parameter :: numeric$ = 29, variable$ = 30
integer, parameter :: mass_of$ = 31, charge_of$ = 32, anomalous_moment_of$ = 33, species$ = 34, species_const$ = 35
integer, parameter :: sinc$ = 36, constant$ = 37, comma$ = 38, rms$ = 39, average$ = 40, sum$ = 41, l_func_parens$ = 42
integer, parameter :: arg_count$ = 43, antiparticle$ = 44, cot$ = 45, sec$ = 46, csc$ = 47, sign$ = 48
integer, parameter :: sinh$ = 49, cosh$ = 50, tanh$ = 51, coth$ = 52, asinh$ = 53, acosh$ = 54, atanh$ = 55, acoth$ = 56
integer, parameter :: min$ = 57, max$ = 58, modulo$ = 59

! Names beginning with "?!+" are place holders that will never match to anything in an expression string.
! Note: "min", "max", "rms" and "average" are not implemented in Bmad but is used by Tao.

character(20), parameter :: expression_op_name(59) = [character(20) :: '+', '-', '*', '/', &
                                    '(', ')', '^', '-', '+', '', 'sin', 'cos', 'tan', &
                                    'asin', 'acos', 'atan', 'abs', 'sqrt', 'log', 'exp', 'ran', &
                                    'ran_gauss', 'atan2', 'factorial', 'int', 'nint', 'floor', 'ceiling', &
                                    '?!+Numeric', '?!+Variable', 'mass_of', 'charge_of', 'anomalous_moment_of', &
                                    'species', '?!+Species', 'sinc', '?!+Constant', ',', 'rms', 'average', 'sum', &
                                    '(', '?!+Arg Count', 'antiparticle', 'cot', 'sec', 'csc', 'sign', &
                                    'sinh', 'cosh', 'tanh', 'coth', 'asinh', 'acosh', 'atanh', 'acoth', 'min', 'max', 'modulo']

integer, parameter :: expression_eval_level(59) = [1, 1, 2, 2, 0, 0, 4, 3, 3, -1, &
              9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, &
              9, 9, 9, 9, 0, 9, 9, 9, 0, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9]

contains

!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
!+
! Function next_in_branch (ix_branch) result (ele, next_ele)
!
! Routine to return a pointer to the next element in a given lattice branch.
! This routine is a type procedure of ele_struct. That is, call looks like:
!     type (ele_struct) ele
!     ...
!     ele%next_in_branch(0)
!
! This routine is meant as an iterator over all elements in a lattice branch including lord elements.
! See: pointer_to_next_ele for a routine that iterates over the tracking elements.
!
! Note: Elements with field overlap are not considered.
! 
! Note: The ix_branch argument is needed since some lord elements are shared by multiple branches.
!
! Input:
!   ix_branch     -- integer: Branch index.
!
! Output:
!   next_ele      -- ele_struct, pointer: Pointer to next element. Will be null after last element.
!-

function next_in_branch (ele, ix_branch) result (next_ele)

implicit none

class (ele_struct), target :: ele
type (ele_struct), pointer :: next_ele
type (branch_struct), pointer :: branch

integer ix_branch

! If the next element is in the tracking part of the lattice

branch => ele%branch
if (ele%ix_ele < branch%n_ele_track) then
  next_ele => branch%ele(ele%ix_ele+1)
  return
endif

! Point to next element which will be in the lord section

nullify(next_ele)

if (ele%ix_ele == branch%n_ele_track) then
  branch => branch%lat%branch(0)
  if (branch%n_ele_max == branch%n_ele_track) return
  next_ele => branch%ele(branch%n_ele_track+1)

elseif (ele%ix_ele == branch%n_ele_max) then
  return

else
  next_ele => branch%ele(ele%ix_ele+1)
endif

! Now iterate until we find an element associated with the branch

do
  if (check_this_lord_in_branch(next_ele)) return

  if (next_ele%ix_ele == branch%n_ele_max) then
    nullify(next_ele)
    return
  else
    next_ele => branch%ele(next_ele%ix_ele+1)
  endif
enddo

!---------------------------------------
contains

recursive function check_this_lord_in_branch(lord) result (is_in_branch)

type (ele_struct) lord
type (ele_struct), pointer :: slave
type (branch_struct), pointer :: branch

integer i
logical is_in_branch

!

is_in_branch = .true.

do i = 1, lord%n_slave
  slave => pointer_to_slave(lord, i)
  if (slave%ix_ele <= slave%branch%n_ele_track) then
    if (slave%ix_branch == ix_branch) return
  else
    if (check_this_lord_in_branch(slave)) return
  endif
enddo

is_in_branch = .false.

end function check_this_lord_in_branch

end function next_in_branch

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
!
! Output:
!   state_str   -- character(16): String representation.
!-

function coord_state_name (coord_state, one_word) result (state_str)

implicit none

integer coord_state
character(16) state_str
logical, optional :: one_word

!

select case (coord_state)
case (not_set$);     state_str = 'Not_Set'
case (pre_born$);    state_str = 'Pre_Born'
case (alive$);       state_str = 'Alive'
case (lost$);        state_str = 'Lost'
case (lost_neg_x$);  state_str = 'Lost_Neg_X'
case (lost_pos_x$);  state_str = 'Lost_Pos_X'
case (lost_neg_y$);  state_str = 'Lost_Neg_Y'
case (lost_pos_y$);  state_str = 'Lost_Pos_Y'
case (lost_pz$);     state_str = 'Lost_Pz'
case (lost_z$);      state_str = 'Lost_Z'
case default;        state_str = null_name$
end select

end function coord_state_name

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Function is_attribute (ix_attrib, which) result (is_attrib)
!
! Routine to determine if an attribute index corresponds to a control variable for overlys/groups.
!
! Input:
!   ix_attrib -- integer: Attribute index.
!   which     -- integer: control_var$, old_control_var$, all_control_var$, multipole$, elec_multipole$
!
! Output:
!   is_attrib -- logical: True if a control variable
!-

function is_attribute (ix_attrib, which) result (is_attrib)

integer ix_attrib, which
logical is_attrib

!

select case (which)
case (control_var$)
  is_attrib = (ix_attrib > var_offset$ .and. ix_attrib < var_offset$ + n_var_max$)
  
case (old_control_var$)
  is_attrib = (ix_attrib > old_control_var_offset$ .and. ix_attrib < old_control_var_offset$ + n_var_max$)

case (all_control_var$)
  is_attrib = ((ix_attrib > var_offset$ .and. ix_attrib < var_offset$ + n_var_max$) .or. &
               (ix_attrib > old_control_var_offset$ .and. ix_attrib < old_control_var_offset$ + n_var_max$))

case (multipole$)
  is_attrib = (ix_attrib >= k0l$ .and. ix_attrib <= t21$)

case (elec_multipole$)
  is_attrib = (ix_attrib >= a0_elec$ .and. ix_attrib <= b21_elec$)

end select

end function is_attribute

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Function pointer_to_slave (lord, ix_slave, control, lord_type, ix_lord_back, ix_control, ix_ic) result (slave_ptr)
!
! Function to point to a slave of a lord.
! Note: Ramper lords do not have any associated slaves (slaves are assigned dynamically at run time).
!
! If lord_type = all$ (the default) the range for ix_slave is:
!   1 to lord%n_slave                                 for "regular" slaves.
!   lord%n_slave+1 to lord%n_slave+lord%n_slave_field for field overlap slaves.
!
! If lord_type = field_lord$, only the field overlap slaves may be accessed and the range for ix_slave is:
!   1 to lord%n_slave_field  
!
! Also see:
!   pointer_to_lord
!   pointer_to_super_lord
!   pointer_to_ele
!   num_lords
!
! Input:
!   lord             -- ele_struct: Lord element
!   ix_slave         -- integer: Index of the slave in the list of slaves controled by the lord.. 
!   lord_type        -- integer, optional: See above.
!
! Output:
!   slave_ptr      -- ele_struct, pointer: Pointer to the slave.
!                       Nullified if there is an error.
!   control        -- control_struct, pointer, optional: Pointer to control info for this lord/slave relationship.
!                       Nullified if there is an error.
!   ix_lord_back   -- integer, optional: Index back to the lord. That is, pointer_to_lord(slave_ptr, ix_lord_back)
!                       will point back to the lord. Set to -1 if there is an error.
!   ix_control     -- integer, optional: Index in lat%control(:) array the control argument is at.
!   ix_ic          -- integer, optional: Index of the lat%ic(:) element associated with the control argument.
!-

function pointer_to_slave (lord, ix_slave, control, lord_type, ix_lord_back, ix_control, ix_ic) result (slave_ptr)

implicit none

type (ele_struct), target :: lord
type (control_struct), pointer, optional :: control
type (ele_struct), pointer :: slave_ptr
type (control_struct), pointer :: con
type (lat_struct), pointer :: lat

integer, optional :: ix_lord_back, lord_type, ix_control, ix_ic
integer i, ix, ix_slave, icon, ixs

!

ixs = ix_slave
if (present(ix_control)) ix_control = -1
if (present(ix_ic)) ix_ic = -1
if (present(ix_lord_back)) ix_lord_back = -1

if (integer_option(all$, lord_type) == field_lord$) ixs = ixs + lord%n_slave

if (ixs > lord%n_slave+lord%n_slave_field .or. ix_slave < 1) then
  nullify(slave_ptr)
  if (present(control)) nullify(control)
  return
endif

lat => lord%branch%lat
icon = lord%ix1_slave + ixs - 1
con => lat%control(icon)
slave_ptr => lat%branch(con%slave%ix_branch)%ele(con%slave%ix_ele)
if (present(control)) control => con
if (present(ix_control)) ix_control = icon

if (present(ix_ic) .or. present(ix_lord_back)) then
  do i = 1, slave_ptr%n_lord + slave_ptr%n_lord_field
    ix = slave_ptr%ic1_lord + i - 1
    if (lat%ic(ix) == icon) then
      if (present(ix_ic)) ix_ic = ix
      if (present(ix_lord_back)) ix_lord_back = i
      exit
    endif
  enddo
endif

end function pointer_to_slave

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine ele_finalizer(ele)
!
! Finalizer routine for ele_struct instances.
! NOTE: Not currently used.
!
! Input:
!   ele   -- ele_struct: Element to cleanup.
!
! Output:
!   ele   -- ele_struct: Element with pointers deallocated as needed.
!-

subroutine ele_finalizer(ele)

type (ele_struct) ele

!

end subroutine ele_finalizer

end module
