!+
! Module rad_int_common
!
! Module needed:
!   use rad_int_common
!-

#include "CESR_platform.inc"

module rad_int_common               

use ptc_interface_mod
use runge_kutta_mod

! The "cache" is for saving values for g, etc through an element to speed
! up the calculation.

type rad_int_track_point_struct
  real(rp) mat6(6,6)
  real(rp) vec0(6)
  type (coord_struct) ref_orb_in
  type (coord_struct) ref_orb_out
  real(rp) g_x0, g_y0     ! Additional g factors for wiggler/bends.
  real(rp) dgx_dx, dgx_dy   ! bending strength gradiant
  real(rp) dgy_dx, dgy_dy   ! bending strength gradiant
  real(rp) l_pole              
end type

type ele_cache_struct
  type (rad_int_track_point_struct), allocatable :: pt(:)
  real(rp) del_z
  integer ix_ele
end type

type rad_int_cache_struct
  type (ele_cache_struct), allocatable :: ele(:)
  integer, allocatable :: ix_ele(:)
  logical :: set = .false.   ! is being used?
end type

type rad_int_info_struct
  type (lat_struct), pointer :: lat
  type (coord_struct), pointer :: orbit(:)
  type (twiss_struct)  a, b
  type (ele_cache_struct), pointer :: cache_ele ! pointer to cache in use
  real(rp) eta_a(4), eta_b(4)
  real(rp) g, g2          ! bending strength (1/bending_radius)
  real(rp) g_x, g_y       ! components in x-y plane
  real(rp) dg2_x, dg2_y
  integer ix_ele
end type

! This structure stores the radiation integrals for the individual elements except
! lin_norm_emittance_a and lin_norm_emittance_b are running sums.

type rad_int_common_struct
  real(rp), allocatable :: i1(:)          ! Noe: All arrays are indexed from 0
  real(rp), allocatable :: i2(:) 
  real(rp), allocatable :: i3(:) 
  real(rp), allocatable :: i4a(:)
  real(rp), allocatable :: i4b(:)
  real(rp), allocatable :: i5a(:) 
  real(rp), allocatable :: i5b(:) 
  real(rp), allocatable :: lin_i2_E4(:) 
  real(rp), allocatable :: lin_i3_E7(:) 
  real(rp), allocatable :: lin_i5a_E6(:) 
  real(rp), allocatable :: lin_i5b_E6(:) 
  real(rp), allocatable :: lin_norm_emit_a(:)  ! Running sum
  real(rp), allocatable :: lin_norm_emit_b(:)  ! Running sum
  real(rp), allocatable :: n_steps(:)      ! number of qromb steps needed
end type

type (rad_int_common_struct), target, save :: ric
type (rad_int_cache_struct), target, save :: rad_int_cache_common(0:10)

!

type synch_rad_common_struct
  type (ele_struct), pointer :: ele0    ! Previous element. For i5 calc.
  real(rp) :: scale = 1.0               ! used to scale the radiation
  real(rp) :: i2 = 0, i3 = 0            ! radiation integrals
  real(rp) :: i5a = 0, i5b = 0
  logical :: i_calc_on = .false.        ! For calculating i2 and i3    
end type

type (synch_rad_common_struct), save :: synch_rad_com

real(rp), parameter :: rad_fluct_const = 1.3231 * r_e * h_bar_planck * &
                                                       c_light / (m_electron * e_charge)

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! Subroutine qromb_rad_int(do_int, pt, info, int_tot)
!
! Function to do integration using Romberg's method on the 7 radiation 
! integrals.
! This is a modified version of QROMB from Num. Rec.
! See the Num. Rec. book for further details.
!
! This routine is only meant to be called by radiation_integrals and
! is not meant for general use.
!
! There are up to 7 integrals that are calculated:
!          I1, I2, I3, I4a, I4b, I5a, I5b
! If do_int(1:7) is False for an integral that means that the integral has
! been calculated by the calling routine using a formula and therefore does
! not have to be done by this routine.
!-

subroutine qromb_rad_int (do_int, pt, info, int_tot)

use precision_def
use nrtype
use nr, only: polint

implicit none

type (ele_struct), pointer :: ele
type (coord_struct) start, end
type (rad_int_track_point_struct) pt
type (rad_int_info_struct) info

integer, parameter :: jmax = 14
integer j, j0, n, n_pts, ix_ele

real(rp) :: int_tot(7)
real(rp) :: eps_int, eps_sum
real(rp) :: ll, del_z, l_ref, z_pos, dint, d0, d_max
real(rp) i_sum(7), rad_int(7)

logical do_int(7), complete

type ri_struct
  real(rp) h(0:jmax)
  real(rp) sum(0:jmax)
end type

type (ri_struct) ri(7)

!

ele => info%lat%ele(info%ix_ele)

eps_int = 1e-4
eps_sum = 1e-6

ri(:)%h(0) = 4.0
ri(:)%sum(0) = 0
rad_int = 0

ll = ele%value(l$)

    ix_ele = info%ix_ele
start = info%orbit(ix_ele-1)
end   = info%orbit(ix_ele)

! Go to the local element frame if there has been caching.
if (associated(info%cache_ele)) then
  call offset_particle (ele, info%lat%param, start, set$, &
       set_canonical = .false., set_multipoles = .false., set_hvkicks = .false.)
  call offset_particle (ele, info%lat%param, end, set$, &
       set_canonical = .false., set_multipoles = .false., set_hvkicks = .false., s_pos = ll)
endif

! Loop until integrals converge.
! ri(k) holds the info for the k^th integral.

do j = 1, jmax

  ri(:)%h(j) = ri(:)%h(j-1) / 4

!---------------
! This is trapzd from Numerical Recipes

  if (j == 1) then
    n_pts = 2
    del_z = ll
    l_ref = 0
  else
    n_pts = 2**(j-2)
    del_z = ll / n_pts
    l_ref = del_z / 2
  endif

  i_sum = 0

  do n = 1, n_pts
    z_pos = l_ref + (n-1) * del_z
    call propagate_part_way (start, end, pt, info, z_pos, j, n)
    i_sum(1) = i_sum(1) + info%g_x * (info%eta_a(1) + info%eta_b(1)) + &
                  info%g_y * (info%eta_a(3) + info%eta_b(3))
    i_sum(2) = i_sum(2) + info%g2
    i_sum(3) = i_sum(3) + info%g2 * info%g
    i_sum(4) = i_sum(4) + &
                  info%g2 * (info%g_x * info%eta_a(1) + info%g_y * info%eta_a(3)) + &
                  info%dg2_x * info%eta_a(1) + info%dg2_y * info%eta_a(3)
    i_sum(5) = i_sum(5) + &
                  info%g2 * (info%g_x * info%eta_b(1) + info%g_y * info%eta_b(3)) + &
                  info%dg2_x * info%eta_b(1) + info%dg2_y * info%eta_b(3)
    i_sum(6) = i_sum(6) + &
                  info%g2 * info%g * (info%a%gamma * info%a%eta**2 + &
                  2 * info%a%alpha * info%a%eta * info%a%etap + &
                  info%a%beta * info%a%etap**2)
    i_sum(7) = i_sum(7) + &
                  info%g2 * info%g * (info%b%gamma * info%b%eta**2 + &
                  2 * info%b%alpha * info%b%eta * info%b%etap + &
                  info%b%beta * info%b%etap**2)
  enddo

  ri(:)%sum(j) = (ri(:)%sum(j-1) + del_z * i_sum(:)) / 2

!--------------
! Back to qromb.
! For j >= 3 we test if the integral calculation has converged.
! Exception: Since wigglers have a periodic field, the calculation can 
! fool itself if we stop before j = 5.

  if (j < 3) cycle
  if (ele%key == wiggler$ .and. j < 5) cycle

  j0 = max(j-4, 1)

  complete = .true.
  d_max = 0

  do n = 1, 7
    if (.not. do_int(n)) cycle
    call polint (ri(n)%h(j0:j), ri(n)%sum(j0:j), 0.0_rp, rad_int(n), dint)
    d0 = eps_int * abs(rad_int(n)) + eps_sum * abs(int_tot(n))
    if (abs(dint) > d0)  complete = .false.
    if (d0 /= 0) d_max = abs(dint) / d0
  enddo

! If we have convergance or we are giving up (when j = jmax) then 
! stuff the results in the proper places.

  if (complete .or. j == jmax) then

    ric%n_steps(ix_ele) = j

    ! Note that ric%i... may already contain a contribution from edge
    ! affects (Eg bend face angles) so add it on to rad_int(i)

    ric%i1(ix_ele)  = ric%i1(ix_ele)  + rad_int(1)
    ric%i2(ix_ele)  = ric%i2(ix_ele)  + rad_int(2)
    ric%i3(ix_ele)  = ric%i3(ix_ele)  + rad_int(3)
    ric%i4a(ix_ele) = ric%i4a(ix_ele) + rad_int(4)
    ric%i4b(ix_ele) = ric%i4b(ix_ele) + rad_int(5)
    ric%i5a(ix_ele) = ric%i5a(ix_ele) + rad_int(6)
    ric%i5b(ix_ele) = ric%i5b(ix_ele) + rad_int(7)

    int_tot(1) = int_tot(1) + ric%i1(ix_ele)
    int_tot(2) = int_tot(2) + ric%i2(ix_ele)
    int_tot(3) = int_tot(3) + ric%i3(ix_ele)
    int_tot(4) = int_tot(4) + ric%i4a(ix_ele)
    int_tot(5) = int_tot(5) + ric%i4b(ix_ele)
    int_tot(6) = int_tot(6) + ric%i5a(ix_ele)
    int_tot(7) = int_tot(7) + ric%i5b(ix_ele)

  endif

  if (complete) return

end do

! We should not be here

print *, 'QROMB_RAD_INT: Note: Radiation Integral is not converging', d_max
print *, '     For element: ', ele%name

end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

subroutine propagate_part_way (start, end, pt, info, z_here, j_loop, n_pt)

implicit none

type (coord_struct) orb, start, end, orb0, orb1
type (ele_struct), pointer :: ele0, ele
type (ele_struct), save ::runt
type (twiss_struct) a0, b0, a1, b1
type (rad_int_info_struct) info
type (rad_int_track_point_struct) pt, pt0, pt1

real(rp) z_here, v(4,4), v_inv(4,4), s1, s2, error
real(rp) f0, f1, del_z, c, s, x, y
real(rp) eta_a0(4), eta_a1(4), eta_b0(4), eta_b1(4)
real(rp) dz, z1

integer i0, i1
integer i, ix, j_loop, n_pt, n, n1, n2
integer, save :: ix_ele = -1

! Init

ele0 => info%lat%ele(info%ix_ele-1)
ele  => info%lat%ele(info%ix_ele)

if (ix_ele /= info%ix_ele) then
  runt = ele
  ix_ele = info%ix_ele
endif

!--------------------------------------
! With caching
! Remember that here the start and stop coords, etc. are in the local ref frame.

if (associated(info%cache_ele)) then

  ! find cached point info near present z position

  del_z = info%cache_ele%del_z
  i0 = int(z_here/del_z)
  f1 = (z_here - del_z*i0) / del_z 
  f0 = 1 - f1
  if (ele%key == wiggler$ .and. ele%sub_key == periodic_type$) i0 = modulo (i0, 10)
  i1 = i0 + 1
  if (i1 > ubound(info%cache_ele%pt, 1)) i1 = i0  ! can happen with roundoff
  pt0 = info%cache_ele%pt(i0)
  pt1 = info%cache_ele%pt(i1)

  orb0%vec = matmul(pt0%mat6, start%vec) + pt0%vec0
  orb1%vec = matmul(pt1%mat6, start%vec) + pt1%vec0

  ! Interpolate information

  pt%dgx_dx = f0 * pt0%dgx_dx + f1 * pt1%dgx_dx
  pt%dgx_dy = f0 * pt0%dgx_dy + f1 * pt1%dgx_dy
  pt%dgy_dx = f0 * pt0%dgy_dx + f1 * pt1%dgy_dx
  pt%dgy_dy = f0 * pt0%dgy_dy + f1 * pt1%dgy_dy

  info%g_x   = f0 * (pt0%g_x0 + pt0%dgx_dx * (orb0%vec(1) - pt0%ref_orb_out%vec(1)) + &
                                pt0%dgx_dy * (orb0%vec(3) - pt0%ref_orb_out%vec(3))) + &
               f1 * (pt1%g_x0 + pt1%dgx_dx * (orb1%vec(1) - pt1%ref_orb_out%vec(1)) + &
                                pt1%dgx_dy * (orb1%vec(3) - pt1%ref_orb_out%vec(3)))
  info%g_y   = f0 * (pt0%g_y0 + pt0%dgy_dx * (orb0%vec(1) - pt0%ref_orb_out%vec(1)) + &
                                pt0%dgy_dy * (orb0%vec(3) - pt0%ref_orb_out%vec(3))) + &
               f1 * (pt1%g_y0 + pt1%dgy_dx * (orb1%vec(1) - pt1%ref_orb_out%vec(1)) + &
                                pt1%dgy_dy * (orb1%vec(3) - pt1%ref_orb_out%vec(3)))
               
  info%dg2_x = 2 * (info%g_x * pt%dgx_dx + info%g_y * pt%dgy_dx)
  info%dg2_y = 2 * (info%g_x * pt%dgx_dy + info%g_y * pt%dgy_dy) 
  info%g2 = info%g_x**2 + info%g_y**2
  info%g  = sqrt(info%g2)

  ! Now convert the g calc back to lab coords.
  
  if (ele%value(tilt_tot$) /= 0) then
    c = cos(ele%value(tilt_tot$)); s = sin(ele%value(tilt_tot$)) 
    x = info%g_x; y = info%g_y
    info%g_x = c * x - s * y
    info%g_y = s * x + c * y
    x = info%dg2_x; y = info%dg2_y
    info%dg2_x = c * x - s * y
    info%dg2_y = s * x + c * y
  endif

  ! Interpolate the dispersions

  runt%mat6 = pt0%mat6
  runt%vec0 = pt0%vec0
  runt%ref_orb_in  = pt0%ref_orb_in
  runt%ref_orb_out = pt0%ref_orb_out

  call mat6_add_offsets (runt)  ! back to lab coords
  call twiss_propagate1 (ele0, runt)
  a0 = runt%a; b0 = runt%b
  call make_v_mats (runt, v, v_inv)
  eta_a0 = matmul(v, (/ runt%a%eta, runt%a%etap, 0.0_rp,   0.0_rp    /))
  eta_b0 = matmul(v, (/ 0.0_rp,   0.0_rp,    runt%b%eta, runt%b%etap /))

  runt%mat6 = pt1%mat6
  runt%vec0 = pt1%vec0
  runt%ref_orb_in  = pt1%ref_orb_in
  runt%ref_orb_out = pt1%ref_orb_out

  call mat6_add_offsets (runt)  ! back to lab coords
  call twiss_propagate1 (ele0, runt)
  a1 = runt%a; b1 = runt%b
  call make_v_mats (runt, v, v_inv)
  eta_a1 = matmul(v, (/ runt%a%eta, runt%a%etap, 0.0_rp,   0.0_rp    /))
  eta_b1 = matmul(v, (/ 0.0_rp,   0.0_rp,    runt%b%eta, runt%b%etap /))

  info%a%beta  = a0%beta  * f0 + a1%beta  * f1
  info%a%alpha = a0%alpha * f0 + a1%alpha * f1
  info%a%gamma = a0%gamma * f0 + a1%gamma * f1
  info%a%eta   = a0%eta   * f0 + a1%eta   * f1
  info%a%etap  = a0%etap  * f0 + a1%etap  * f1

  info%b%beta  = b0%beta  * f0 + b1%beta  * f1
  info%b%alpha = b0%alpha * f0 + b1%alpha * f1
  info%b%gamma = b0%gamma * f0 + b1%gamma * f1
  info%b%eta   = b0%eta   * f0 + b1%eta   * f1
  info%b%etap  = b0%etap  * f0 + b1%etap  * f1

  info%eta_a = eta_a0 * f0 + eta_a1 * f1
  info%eta_b = eta_b0 * f0 + eta_b1 * f1

  return

endif

!--------------------------------------
! No caching

dz = 1e-3
z1 = z_here + dz
if (z1 > ele%value(l$)) z1 = max(0.0_rp, z_here - dz)

call twiss_and_track_partial (ele0, ele, info%lat%param, z_here, runt, start, orb0)
call twiss_and_track_partial (ele0, ele, info%lat%param, z1, start = start, end = orb1)
info%a = runt%a
info%b = runt%b

!

call make_v_mats (runt, v, v_inv)

info%eta_a = matmul(v, (/ info%a%eta, info%a%etap, 0.0_rp,   0.0_rp /))
info%eta_b = matmul(v, (/ 0.0_rp,   0.0_rp,    info%b%eta, info%b%etap /))

if (ele%key == wiggler$ .and. ele%sub_key == map_type$) then 
  call calc_wiggler_g_params (ele, z_here, orb, pt, info)
else
  info%g_x   = pt%g_x0 - (orb1%vec(2) - orb0%vec(2)) / (z1 - z_here)
  info%g_y   = pt%g_y0 - (orb1%vec(4) - orb0%vec(4)) / (z1 - z_here)
endif

info%dg2_x = 2 * (info%g_x * pt%dgx_dx + info%g_y * pt%dgy_dx)
info%dg2_y = 2 * (info%g_x * pt%dgx_dy + info%g_y * pt%dgy_dy) 

info%g2 = info%g_x**2 + info%g_y**2
info%g = sqrt(info%g2)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine calc_wiggler_g_params (ele, z, orb, pt, info)

implicit none

type (coord_struct) orb
type (rad_int_track_point_struct) pt
type (rad_int_info_struct) info
type (ele_struct) ele

real(rp) dk(3,3), z
real(rp) kick_0(6)

! 

call em_field_kick (ele, info%lat%param, z, orb%vec, .false., kick_0, dk)

pt%g_x0 = -kick_0(2)
pt%g_y0 = -kick_0(4)
pt%dgx_dx = -dk(1,1)
pt%dgx_dy = -dk(1,2)
pt%dgy_dx = -dk(2,1)
pt%dgy_dy = -dk(2,2)

info%g_x = -kick_0(2)
info%g_y = -kick_0(4)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine transfer_rad_int_struct (rad_int_in, rad_int_out)

implicit none

type (rad_int_common_struct) rad_int_in, rad_int_out
integer n

!

n = ubound(rad_int_in%i1, 1)

call re_allocate (rad_int_out%i1, 0, n)
call re_allocate (rad_int_out%i2, 0, n)
call re_allocate (rad_int_out%i3, 0, n)
call re_allocate (rad_int_out%i4a, 0, n)
call re_allocate (rad_int_out%i4b, 0, n)
call re_allocate (rad_int_out%i5a, 0, n)
call re_allocate (rad_int_out%i5b, 0, n)
call re_allocate (rad_int_out%n_steps, 0, n)
call re_allocate (rad_int_out%lin_i2_e4, 0, n)
call re_allocate (rad_int_out%lin_i3_e7, 0, n)
call re_allocate (rad_int_out%lin_i5a_e6, 0, n)
call re_allocate (rad_int_out%lin_i5b_e6, 0, n)

rad_int_out = rad_int_in

end subroutine

end module
