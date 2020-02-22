!+
! subroutine init_ray (ray, branch, ix_ele, l_offset, orb, direction)
!
! subroutine to start a new synch radiation ray
!
! Modules needed:
!   use synrad_mod
!
! Input:
!   branch   -- branch_struct with twiss propagated and mat6s made
!   ix_ele -- integer: index of lattice element to start ray from
!   l_offset -- real(rp): offset along the length of the element 
!                         to use as a starting point for ray
!   orb(0:*) -- coord_struct: orbit of particles to use as 
!                             source of ray
!   direction -- integer: +1 for in direction of s
!                         -1 for against s
!
! Output:
!   ray    -- ray_struct: synch radiation ray with starting
!                         parameters set
!-

subroutine init_ray (ray, branch, ix_ele, l_offset, orb, direction)

use synrad_interface, except => init_ray
use em_field_mod

implicit none

type (branch_struct), target :: branch
type (coord_struct) :: orb(0:*)
type (coord_struct), save :: orb0
type (ele_struct), pointer :: ele0, ele
type (ele_struct), save :: runt_ele, d_ele
type (ray_struct) :: ray
type (em_field_struct) :: field

real(rp) l_offset, g(3)
real(rp), save :: l_start = 0
integer direction, ix_ele
integer, save :: ix_ele_old = -1
logical err_flag

! initialize the ray

ele0 => branch%ele(ix_ele-1)
ele => branch%ele(ix_ele)

! Get the ray's initial twiss values.
! Starting from the old saved position saves time when tracking though a map_type wiggler.

if (ix_ele_old /= ix_ele .or. l_offset < l_start) then
  call transfer_ele (ele0, runt_ele)
  l_start = 0
  orb0 = orb(ix_ele-1)
endif

call twiss_and_track_intra_ele (ele, branch%param, l_start, l_offset, &
                   .true., .true., orb0, orb0, runt_ele, runt_ele, compute_floor_coords = .true.)
ray%x_twiss = runt_ele%a
ray%y_twiss = runt_ele%b

l_start = l_offset
ix_ele_old = ix_ele

! set the ray's g_bend value (inverse bending radius at src pt) 

call g_bending_strength_from_em_field (ele, branch%param, l_offset, orb0, .false., g)
ray%g_bend = norm2(g)

ray%start%vec(1) = orb0%vec(1)
ray%start%vec(2) = direction * orb0%vec(2)
ray%start%vec(3) = orb0%vec(3)
ray%start%vec(4) = direction * orb0%vec(4)
ray%start%vec(5) = l_offset
ray%start%vec(6) = direction * sqrt (1 - ray%start%vec(2)**2 - ray%start%vec(4)**2)
ray%start%ix_ele = ix_ele
ray%start%s = ele0%s + l_offset

ray%now = ray%start

ray%direction = direction
ray%track_len = 0

ray%start_floor = coords_curvilinear_to_floor ([ray%start%vec(1), 0.0_rp, ray%start%s], ele%branch, err_flag)
ray%start_floor%theta = runt_ele%floor%theta + atan2(ray%start%vec(2), ray%start%vec(6))
ray%start_floor%theta = modulo2(ray%start_floor%theta, pi)

ray%now_floor = ray%start_floor

ix_ray = ix_ray + 1

end subroutine
