module synrad3d_struct

use bmad_struct
use bmad_interface

! This structure defines a photon at a particular point.
! for vec(6): (x, y, s) is the local coordinate system with
!  s being the longitudinal position (s = 0 is the start of the lattice), 
!  and  x, y are the local transverse coords. See the Bmad manual for more details.
! Notice that vec(1)^2 + vec(3)^2 + vec(5)^2 = 1


type photon3d_coord_struct
  real(rp) vec(6)             ! Photon position: (x, vx/c, y, vy/c, s, vz/c)
  real(rp) energy             ! In eV
  real(rp) track_len          ! length of the track from the start
  integer ix_ele              ! index of lattice element we are in.
  integer ix_wall             ! Index to wall cross-section array
end type

! This structure defines the full track of the photon from start to finish
! %start -- Starting position.
! %old   -- Used by the tracking code. Not useful otherwise.
! %now   -- Present position. At the end of tracking, %now will be the final position.
! %reflect(:)      -- Records the positions at which the photon was reflected off the wall
!                       including the final position. The array bounds are:
!                       %reflect(0:%n_reflect+1)
! %intensity       -- Intensity of this macro-photon in Photons/(beam_particle*turn).
! %crossed_lat_end -- Did the photon cross from the end of the lattice to the beginning
!                       or vice versa?
! %hit_antechamber -- Did the photon hit the antechamber at the final position?
! %ix_photon       -- Photon index. The first photon generated has index 1, etc.
! %n_reflect       -- Number of reflections. %reflect

type photon3d_track_struct
  type (photon3d_coord_struct) start, old, now  ! coords:
  type (photon3d_coord_struct), allocatable :: reflect(:) ! Photon reflection points
  real(rp) intensity          ! Intensity of this macro-photon in Photons/(beam_particle*turn)
  logical :: crossed_lat_end = .false.     ! Photon crossed through the lattice beginning or end?
  logical :: hit_antechamber = .false.     
  integer ix_photon                        ! Photon index.
  integer n_reflect                        ! Number of reflections
end type

!--------------
! The wall is specified by an array of points at given s locations.
! The wall between point i-1 and i is associated with wall%pt(i) (see the photon3d_coord_struct).
! If there is an antechamber: width2_plus and width2_minus are the antechamber horizontal extent.
! With no antechamber: width2_plus and width2_minus specify beam stops.

type wall3d_pt_struct
  real(rp) s                      ! Longitudinal position.
  character(16) basic_shape       ! "elliptical" or "rectangular".
  real(rp) width2                 ! Half width ignoring antechamber.
  real(rp) height2                ! Half height ignoring antechamber.
  real(rp) ante_height2_plus      ! Antechamber half height on +x side of the wall
  real(rp) width2_plus            ! Distance from pipe center to +x side edge.
  real(rp) ante_height2_minus     ! Antechamber half height on -x side of the wall
  real(rp) width2_minus           ! Distance from pipe center -x side edge.
  real(rp) ante_x0_plus           ! Computed: x coord at +x antechamber opening.
  real(rp) ante_x0_minus          ! Computed: x coord at -x antechamber opening.
  real(rp) y0_plus                ! Computed: y coord at edge of +x beam stop.
  real(rp) y0_minus               ! Computed: y coord at edge of -x beam stop.
end type

! This is just an array of chamber cross-sections.

type wall3d_struct
  type (wall3d_pt_struct), allocatable :: pt(:)
  integer n_pt_max
end type

! Some parameters that can be set. 

type sr3d_params_struct
  real(rp) :: ds_track_step_max = 3     ! Maximum longitudinal distance in one photon "step".
  real(rp) :: dr_track_step_max = 0.1   ! Maximum tranverse distance in one photon "step".
  logical :: allow_reflections = .true. ! If False, terminate tracking when photon hits the wall.
  logical :: debug = .false.           
  logical :: stop_if_hit_antechamber = .false. 
end type

type (sr3d_params_struct), save :: sr3d_params

end module
