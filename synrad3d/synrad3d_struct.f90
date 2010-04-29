module synrad3d_struct

use bmad_struct
use bmad_interface
use twiss_and_track_mod

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

type photon3d_wall_hit_struct
  type (photon3d_coord_struct) before_reflect   ! Coords before reflection.
  type (photon3d_coord_struct) after_reflect    ! Coords after reflection.
  real(rp) dw_perp(3)                   ! Wall perpendicular vector
  real(rp) cos_perp                     ! Cosine of hit angle
  real(rp) reflectivity                 ! Reflectivity probability
end type

! This structure defines the full track of the photon from start to finish
! %start -- Starting position.
! %old   -- Used by the tracking code. Not useful otherwise.
! %now   -- Present position. At the end of tracking, %now will be the final position.
! %wall_hit(:)     -- Records the positions at which the photon hit the wall
!                       including the final position. The array bounds are: %hit(1:%n_hit) 
! %intensity       -- Intensity of this macro-photon in Photons/(beam_particle*turn).
! %crossed_lat_end -- Did the photon cross from the end of the lattice to the beginning
!                       or vice versa?
! %hit_antechamber -- Did the photon hit the antechamber at the final position?
! %ix_photon       -- Unfiltered photon index. 
! %ix_photon_generated -- The first photon generated has index 1, etc.
! %n_wall_hit      -- Number of wall hits.

type photon3d_track_struct
  type (photon3d_coord_struct) start, old, now  ! coords:
  real(rp) intensity          ! Intensity of this macro-photon in Photons/(beam_particle*turn)
  logical :: crossed_lat_end = .false.     ! Photon crossed through the lattice beginning or end?
  logical :: hit_antechamber = .false.     
  integer ix_photon                        ! Photon index.
  integer ix_photon_generated
  integer :: n_wall_hit = 0                ! Number of wall hits
end type

!--------------
! The wall is specified by an array of points at given s locations.
! The wall between point i-1 and i is associated with wall%pt(i) (see the photon3d_coord_struct).
! If there is an antechamber: width2_plus and width2_minus are the antechamber horizontal extent.
! With no antechamber: width2_plus and width2_minus specify beam stops.

type wall3d_pt_struct
  real(rp) s                      ! Longitudinal position.
  character(16) basic_shape       ! "elliptical", "rectangular", or "polygon"
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
  integer ix_polygon              ! Polygon index 
end type

! Structure used to hold polygon shapes

type wall3d_polygon_struct
  type (polygon_vertex_struct), allocatable :: v(:)
end type

! This is just an array of chamber cross-sections.

type wall3d_struct
  type (wall3d_pt_struct), allocatable :: pt(:)
  type (wall3d_polygon_struct), allocatable :: polygon(:)
  integer n_pt_max
end type

! Some parameters that can be set. 

type sr3d_params_struct
  real(rp) :: ds_track_step_max = 3     ! Maximum longitudinal distance in one photon "step".
  real(rp) :: dr_track_step_max = 0.1   ! Maximum tranverse distance in one photon "step".
  logical :: allow_reflections = .true. ! If False, terminate tracking when photon hits the wall.
  logical :: stop_if_hit_antechamber = .false. 
  logical :: allow_absorbtion = .true.  ! If False, do not allow photon to be adsorbed at wall. 
                                        !   Used for debugging.
  logical :: debug_on
end type

type (sr3d_params_struct), save :: sr3d_params

end module
