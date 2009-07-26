module synrad3d_struct

use bmad_struct
use bmad_interface


integer, parameter :: elliptical$ = 1, rectangular$ = 2

type photon3d_coord_struct
  real(rp) vec(6)             ! Position: (x, vx/c, y, vy/c, z, vz/c)
  real(rp) energy             ! eV
  real(rp) intensity          ! Watts
  real(rp) radius             ! Normalized transverse position of the photon 
                              ! relative to the wall. 1 => at wall.
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
  integer type         ! elliptical$ or rectangular$
  real(rp) width2      ! half width
  real(rp) height2     ! half height
  real(rp) s           ! Longitudinal position.
end type

type wall_3d_struct
  type (wall_3d_pt_struct), allocatable :: pt(:)
  integer n_pt_max
end type

!

type synrad3d_params_struct
  real(rp) :: ds_track_step_max = 10
  real(rp) :: dr_track_step_max = 0.1
end type

type (synrad3d_params_struct) synrad3d_params

end module
