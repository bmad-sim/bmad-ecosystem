module synrad3d_struct

use bmad_struct
use bmad_interface

type photon3d_coord_struct
  real(rp) vec(6)             ! Position: (x, vx/c, y, vy/c, z, vz/c)
  real(rp) energy             ! eV
  real(rp) track_len          ! length of the track from the start
  integer ix_ele              ! index of element we are in.
  integer ix_wall             ! Index to wall segment
end type

type photon3d_track_struct
  type (photon3d_coord_struct) start, old, now  ! coords
  type (photon3d_coord_struct), allocatable :: reflect(:)
  real(rp) intensity          ! Intensity of this macro-photon in Photons/(beam_particle*turn)
  logical :: crossed_lat_end = .false.     ! photon crossed through the lattice beginning or end?
  logical :: hit_antichamber = .false.
  integer ix_photon           ! Photon index.
  integer n_reflect
end type

!--------------
! The wall is specified by an array of points as given s locations.
! The wall between point i-1 and i is associated with 

type wall3d_pt_struct
  real(rp) s           ! Longitudinal position.
  character(16) type   ! elliptical or rectangular
  real(rp) width2      ! half width
  real(rp) height2     ! half height
  real(rp) antichamber_plus_x_height2
  real(rp) antichamber_minus_x_height2
end type

type wall3d_struct
  type (wall3d_pt_struct), allocatable :: pt(:)
  integer n_pt_max
end type

!

type synrad3d_params_struct
  real(rp) :: ds_track_step_max = 3
  real(rp) :: dr_track_step_max = 0.1
  logical :: allow_reflections = .true.
  logical :: debug = .false.
end type

type (synrad3d_params_struct), save :: synrad3d_params

end module
