module synrad_struct

use bmad_interface

! The wall is specified by an array of points with straight lines (faces) between the points. 
! the face between point i-1 and the point i is associated with the structure for point i.
! If a face represents an opening in the wall (eg the openings for the chess beam 
! lines) then that face is call a phantom.

! The "triangle" point is constructed so that the (possibly curved) wall 
! associated with a wall point wall%pt(i) is containd within the triangle 
! wall%pt(i)%r_floor, wall%pt(i)%r_floor_tri and wall%pt(i-1)%r_floor.

type wall_pt_struct           ! struct for input points
  character(40) :: name = ''          ! name of element (sliding_joint, etc.)
  real(rp) :: x = 0, s = 0            ! position of wall point
  real(rp) :: r_floor(3) = 0          ! Floor position.
  real(rp) :: r_floor_tri(3) = 0      ! Triangle point.
  integer :: ix_pt = 0                ! Self index int wall%pt array
  integer :: n_seg = 0                ! how many segments it will be broken up into
  integer :: ix_seg = 0               ! index to seg(:) array. From ix_seg+1 to ix_seg+n_seg
  logical :: phantom = .false.        ! is the face an opening?
  logical :: next_to_alley = .false.  ! point is near an alley way? 
  logical :: linear_wall = .false.    ! Is associated wall section a linear line segment?
end type wall_pt_struct

! For the synch light power computation each face is broken up into segments.
! It is assumed that the power landing upon a segment is constant
! throughout the segment

type source_struct
  integer :: ix_ele = 0           ! element index at source
  real(rp) :: power_per_len = 0   ! Power from this source
  real(rp) :: s = 0               ! Longitudinal s position.
end type source_struct

! substruct for a segment 

type seg_power_struct           
  real(rp) :: power_tot = 0            ! total power on segment (Watts)
  real(rp) :: power_per_len = 0        ! power density (Watts / m)
  real(rp) :: power_per_area = 0       ! power density (Watts / m^2)
  real(rp) :: photons_per_sec = 0      ! flux hitting segment in photons per sec
  integer :: n_source = 0              ! number of source points
  type (source_struct) main_source  ! main source info for rays hitting this seg
end type seg_power_struct

type wall_seg_struct         ! segment struct
  integer ix_seg             ! Self index in wall%seg array.
  integer ix_pt              ! index to which point owns this segment
  real(rp) s, x              ! s, x position of the segment at the endpoint
  real(rp) r_floor(3)        ! Floor position at end.
  real(rp) r_floor_mid(3)    ! Floor position at midpoint.
  real(rp) theta             ! Angle of segment orientation in global coords
  real(rp) len               ! length of segment
  type (seg_power_struct) power
  integer ix_ele             ! lattice ele at s_mid
  type (twiss_struct) :: a,b ! twiss info at s_mid
end type wall_seg_struct

! A wall is just a collection of points and segments
! %seg(0) is a dummy segment meant to hold position information.

type wall_struct           ! either positive_x or negative_x side wall
  integer side             ! positive_x$ or negative_x$ 
  type (wall_pt_struct), allocatable :: pt(:)    ! Indexed from 0.
  type (wall_seg_struct), allocatable :: seg(:)  ! Indexed from 0.
  integer n_pt_max
  integer n_seg_max
  real(rp) seg_len_max
end type wall_struct

! Include both walls and ends in walls_struct
! (Radiation could escape out ends if lattice is not circular.)
type walls_struct
  type (wall_struct) :: positive_x_wall, negative_x_wall
  type (wall_seg_struct) :: start_end, exit_end
  integer lat_geometry     ! open$ or closed$
  real(rp) s_max
end type walls_struct

! The computation tracks a set of synchrotron light rays from their source to the wall.

type ray_struct       ! struct for a light ray
  integer direction       ! direction of travel, +1 = forward direction
  real(rp) track_len      ! length of the track from the start
  type (coord_struct) start, now  ! coords. %vec(6) willl be negative if direction = -1
  type (floor_position_struct) start_floor, now_floor
  integer ix_wall_pt      ! index of wall point where hit
  integer ix_seg_pt
  real(rp) p1_factor      ! factor for computing the power/length
  real(rp) p2_factor      ! factor for computing the power/area
  real(rp) r_seg          !
  real(rp) g_bend         ! g = |1/rho|  bending radius inverse at source point
  type (twiss_struct) x_twiss, y_twiss    ! twiss at source point
  integer wall_side       ! Which wall is hit: positive_x_side$, negative_x_side$, 
                          !   start_side$, or exit_side$
end type ray_struct

! misc stuff

type ele_power_struct       ! power from a lat element
  real(rp) at_wall              ! power hitting the wall
  real(rp) radiated             ! power radiated
end type ele_power_struct

type synrad_param_struct
  character(100) lat_file
  real(rp) epsilon_y     ! vertical emit
  real(rp) i_beam        ! beam current
  integer n_slice        ! # of slice per element or wiggler pole
  logical :: debug = .false.
  logical filter_phantom_photons
end type synrad_param_struct

type synrad_mode_struct
 type (normal_modes_struct) pos_mode
 type (normal_modes_struct) ele_mode
end type synrad_mode_struct

type outline_pt_struct
  character(16) name
  character(16) blueprint
  real(rp) s, x
  logical phantom
end type outline_pt_struct

type outline_struct
  character(16) name
  character(16) blueprint
  type (outline_pt_struct) in(200), out(200)
  integer n_out, n_in
  integer ix_out_slide, ix_in_slide
  logical has_alley
  logical zero_is_center
  real(rp) s_center
  logical overlay
end type outline_struct

type wall_list_struct
  character(16) name
  integer ix_outline
  real(rp) s, len
end type wall_list_struct

type concat_part_struct
  character(16) name
  integer direction
end type concat_part_struct

type concat_struct
  character(16) name
  type (concat_part_struct) part(50)
  logical overlay
end type concat_struct

character(20) :: wall_name(-2:2) = ['start_side     ', 'negative_x_side', '?????          ', &
                                      'positive_x_side', 'exit_side      ' ]
integer, parameter :: negative_x$ = -1, positive_x$ = 1, start_side$ = -2, exit_side$ = 2
integer, save :: ix_ray = 0

real(rp), parameter :: synrad_significant_length = 1d-8

end module
