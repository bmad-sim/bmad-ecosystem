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

use synrad_struct
use synrad_interface, except => init_ray
use boris_mod

implicit none

type (branch_struct), target :: branch
type (coord_struct) :: orb(0:*)
type (coord_struct), save :: orb0, orb1, orb2
type (ele_struct), pointer :: ele0, ele
type (ele_struct), save :: runt_ele, d_ele
type (ray_struct) :: ray
type (em_field_struct) :: field

real(rp) l_offset, k_wig, g_max, l_small
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

orb1 = orb0

select case (ele%key)
case (sbend$)

  ! sbends are easy
  ray%g_bend = abs(ele%value(g$) + ele%value(g_err$))

! for quads or sol_quads, get the bending radius
! from the change in x' and y' over a small 
! distance in the element

case (quadrupole$, sol_quad$, sad_mult$, elseparator$)
  call transfer_ele (runt_ele, d_ele)
  l_small = 1e-2      ! something small
  d_ele%value(l$) = l_small
  d_ele%value(fringe_at$) = none$
  call make_mat6 (d_ele, branch%param, orb1, orb1)
  call track1 (orb1, d_ele, branch%param, orb2)
  orb2%vec = orb2%vec - orb1%vec
  ray%g_bend = sqrt(orb2%vec(2)**2 + orb2%vec(4)**2) / l_small

! wiggler, undulator

case (wiggler$, undulator$)

  if (ele%sub_key == periodic_type$) then

    ! for periodic wigglers, get the max g_bend from 
    ! the max B field of the wiggler, then scale it 
    ! by the cos of the position along the poles
    ! Note: assumes particles are relativistic!!
    k_wig = twopi * ele%value(n_pole$) / (2 * ele%value(l$))

    g_max = c_light * ele%value(b_max$) / (ele%value(p0c$))
    ray%g_bend = abs(g_max * cos (k_wig * l_offset))
    orb1%vec(2) = orb1%vec(2) + (g_max / k_wig) * sin (k_wig * l_offset)

  else  ! map type

    ! for mapped wigglers, find the B field at the source point
    ! and extract the g_bend
    ! Note: assumes particles are relativistic!!
    call em_field_calc (runt_ele, branch%param, l_offset, orb1, .false., field)

    ray%g_bend = sqrt(sum(field%b(1:2)**2)) * c_light / ele%value(p0c$)

  endif

case default

  print *, 'ERROR: UNKNOWN ELEMENT HERE ', ele%name

end select

ray%start%vec(1) = orb1%vec(1)
ray%start%vec(2) = direction * orb1%vec(2)
ray%start%vec(3) = orb1%vec(3)
ray%start%vec(4) = direction * orb1%vec(4)
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
