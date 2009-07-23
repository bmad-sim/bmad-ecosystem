module synrad3d_struct

use bmad_struct
use bmad_interface

type photon3d_coord_struct
  real(rp) vec(6)             ! Position: (x, vx/c, y, vy/c, z, vz/c)
  real(rp) energy             ! eV
  real(rp) intensity          ! Watts
end type

type photon3d_track_struct
  integer ix_ele          ! index of element we are now tracking through
  integer ix_source       ! element index at source of ray
  real(rp) track_len      ! length of the track from the start
  type (synrad3d_coord_struct) start, old, now  ! coords
  logical crossed_end     ! photon crossed through the lattice end?
  integer n_reflect
end type

!--------------

type wall_3d_pt_struct
  real(rp) width2      ! half width
  real(rp) height2     ! half height
  real(rp) s           ! Longitudinal position.
end type

type wall_3d_struct
  type (wall_3d_pt_struct), allocatable :: pt(:)
  integer n_pt_max
end type

end module
