module synrad3d_struct

use bmad_struct
use bmad_interface

type photon3d_coord_struct
  real(rp) vec(6)             ! Position: (x, vx/c, y, vy/c, z, vz/c)
  real(rp) energy             ! In eV
  real(rp) track_len          ! length of the track from the start
  integer ix_ele              ! index of element we are in.
  integer ix_wall             ! Index to wall segment
end type

type photon3d_track_struct
  type (photon3d_coord_struct) start, old, now  ! coords:
  type (photon3d_coord_struct), allocatable :: reflect(:) ! Records coords at points where 
                                                          !   the photon reflects off of the wall.
  real(rp) intensity          ! Intensity of this macro-photon in Photons/(beam_particle*turn)
  logical :: crossed_lat_end = .false.     ! Photon crossed through the lattice beginning or end?
  logical :: hit_antechamber = .false.     
  integer ix_photon                        ! Photon index.
  integer n_reflect                        ! Number of reflections
end type

!--------------
! The wall is specified by an array of points at given s locations.
! The wall between point i-1 and i is associated with wall%pt(i).
! If there is an antechamber: width2_plus and width2_minus are the antechamber horizontal extent.
! With no antechamber: width2_plus and width2_minus specify beam stops.

type wall3d_pt_struct
  real(rp) s                      ! Longitudinal position.
  character(16) type              ! Elliptical or rectangular.
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
