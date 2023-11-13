!+
! Module rad_int_common
!-

module rad_int_common               

use multipole_mod
use attribute_mod, only: attribute_index

! The "cache" is for saving values for g, etc through an element to speed
! up the calculation.

integer, parameter :: no_cache$ = 0, cache_no_misalign$ = 2

type rad_int_track_point_struct
  real(rp) :: s_body = 0
  real(rp) :: mat6(6,6) = 0
  real(rp) :: vec0(6) = 0
  type (coord_struct) :: ref_orb_in = coord_struct()
  type (coord_struct) :: ref_orb_out = coord_struct()
  real(rp) :: g_x0 = 0, g_y0 = 0      ! Additional g factors for bends.
  real(rp) :: dgx_dx = 0, dgx_dy = 0  ! bending strength gradient
  real(rp) :: dgy_dx = 0, dgy_dy = 0  ! bending strength gradient
end type

! Note: The points may not be evenly spaced.

type rad_int_cache1_struct
  type (rad_int_track_point_struct), allocatable :: pt(:)   ! pt(0:n_pt)
  integer :: n_pt = -1              ! Upper bound of pt(0:n_pt)
  integer :: cache_type = no_cache$
end type

type rad_int_cache_struct
  type (rad_int_cache1_struct), allocatable :: c_ele(:)
  logical :: in_use = .false.
end type

type (rad_int_cache_struct), allocatable, target, save :: rad_int_cache_common(:)

!

type rad_int_info_struct
  type (branch_struct), pointer :: branch
  type (ele_struct), pointer :: ele
  type (coord_struct), pointer :: orbit(:)
  type (twiss_struct)  a, b
  type (rad_int_cache1_struct), pointer :: cache_ele ! pointer to cache in use
  real(rp) eta_a(4), eta_b(4)
  real(rp) g, g2          ! bending strength (1/bending_radius)
  real(rp) g_x, g_y       ! components in x-y plane
  real(rp) dg2_x, dg2_y
end type

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! Subroutine qromb_rad_int(param, do_int, pt, info, int_tot, rad_int)
!
! Function to do integration using Romberg's method on the 7 radiation integrals.
! This is a modified version of QROMB from Num. Rec.
! See the Num. Rec. book for further details.
!
! This routine is only meant to be called by radiation_integrals and
! is not meant for general use.
!
! There are up to 9 integrals that are calculated:
!          I1, I2, I3, I4a, I4b, I5a, I5b, i0, i6b
! If do_int(1:7) is False for an integral that means that the integral has
! been calculated by the calling routine using a formula and therefore does
! not have to be done by this routine.
!-

subroutine qromb_rad_int (param, do_int, pt, info, int_tot, rad_int1)

use precision_def
use super_recipes_mod, only: super_polint

implicit none

integer, parameter :: num_int = 9

type ri_array_struct
  real(rp) h
  real(rp) sum(num_int)
end type

real(rp) d_err(num_int)

type (ri_array_struct) ri_array(0:4) ! ri_array(n) holds info for the integrals for a particular step.
type (ele_struct), pointer :: ele
type (ele_struct) runt
type (coord_struct) start, end
type (rad_int_track_point_struct) pt
type (rad_int_info_struct) info
type (rad_int1_struct) rad_int1, int_tot
type (lat_param_struct) param

integer n,  j, j1, j_min_test, j_max, n_pts

real(rp) :: eps_int, eps_sum, gamma
real(rp) :: ll, del_z, l_ref, z_pos, dint(num_int), d0(num_int), d_max
real(rp) i_sum(num_int), rad_int_vec(num_int), int_tot_vec(num_int)

logical do_int(num_int), converged

character(*), parameter :: r_name = 'qromb_rad_int'

!

ele => info%ele

eps_int = 1d-4
eps_sum = 1d-6

ri_array(0)%h = 4
ri_array(0)%sum = 0
rad_int_vec = 0
int_tot_vec = [int_tot%i1, int_tot%i2, int_tot%i3, int_tot%i4a, &
               int_tot%i4b, int_tot%i5a, int_tot%i5b, int_tot%i0,int_tot%i6b]

ll = ele%value(l$)

start = info%orbit(ele%ix_ele-1)
end   = info%orbit(ele%ix_ele)
gamma = ele%value(e_tot$) / mass_of(param%particle)
call transfer_ele (ele, runt)

! For j >= 3 we test if the integral calculation has converged.
! Exception: Since wigglers have a planar or helical field, the calculation can 
! fool itself if we stop before j = 5.

j_min_test = 3
j_max = 14

if (ele%key == wiggler$ .or. ele%key == undulator$) then
  if (ele%field_calc == planar_model$ .or. ele%field_calc == helical_model$) then
    j_min_test = 4 + log(max(1.0_rp, ele%value(n_period$))) / log(2.0_rp)
    j_max = j_min_test + 8
  else
    j_min_test = 5
    j_max = 16
  endif
endif

!---------------
! Loop until integrals converge.
! This is trapzd from Numerical Recipes.

do j = 1, j_max

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
    call propagate_part_way (start, param, pt, info, z_pos, runt)
    i_sum(1) = i_sum(1) + info%g_x * (info%eta_a(1) + info%eta_b(1)) + &
                          info%g_y * (info%eta_a(3) + info%eta_b(3))
    i_sum(2) = i_sum(2) + info%g2
    i_sum(3) = i_sum(3) + info%g2 * info%g
    i_sum(4) = i_sum(4) + info%g2 * (info%g_x * info%eta_a(1) + info%g_y * info%eta_a(3)) + &
                                   info%dg2_x * info%eta_a(1) + info%dg2_y * info%eta_a(3)
    i_sum(5) = i_sum(5) + info%g2 * (info%g_x * info%eta_b(1) + info%g_y * info%eta_b(3)) + &
                                   info%dg2_x * info%eta_b(1) + info%dg2_y * info%eta_b(3)
    i_sum(6) = i_sum(6) + info%g2 * info%g * (info%a%gamma * info%a%eta**2 + &
                      2 * info%a%alpha * info%a%eta * info%a%etap + info%a%beta * info%a%etap**2)
    i_sum(7) = i_sum(7) + info%g2 * info%g * (info%b%gamma * info%b%eta**2 + &
                      2 * info%b%alpha * info%b%eta * info%b%etap + info%b%beta * info%b%etap**2)
    i_sum(8) = i_sum(8) + gamma * info%g
    i_sum(9) = i_sum(9) + info%g2 * info%g * info%b%beta
  enddo

  if (j <= 4) then
    j1 = j
  else
    ri_array(0:3) = ri_array(1:4)
    j1 = 4
  endif

  ri_array(j1)%h = ri_array(j1-1)%h / 4
  ri_array(j1)%sum = (ri_array(j1-1)%sum + del_z * i_sum) / 2

  !--------------
  ! Back to qromb.

  if (j < j_min_test) cycle

  d_max = 0

  do n = 1, num_int
    if (.not. do_int(n)) cycle
    call super_polint (ri_array(1:j1)%h, ri_array(1:j1)%sum(n), 0.0_rp, rad_int_vec(n), dint(n))
    d0(n) = eps_int * abs(rad_int_vec(n)) + eps_sum * abs(int_tot_vec(n)) + 1d-30
    d_err(n) = abs(dint(n)) / d0(n)
    d_max = max(d_max, d_err(n))
  enddo

  ! If we have convergance or we are giving up (when j = j_max) then 
  ! stuff the results in the proper places.

  converged = (d_max <= 1)
  if (converged .or. j == j_max) then

    rad_int1%n_steps = j

    ! Note that rad_int%i... may already contain a contribution from edge
    ! affects (Eg bend face angles) so add it on to rad_int_vec(i)

    rad_int1%i1  = rad_int1%i1  + rad_int_vec(1)
    rad_int1%i2  = rad_int1%i2  + rad_int_vec(2)
    rad_int1%i3  = rad_int1%i3  + rad_int_vec(3)
    rad_int1%i4a = rad_int1%i4a + rad_int_vec(4)
    rad_int1%i4b = rad_int1%i4b + rad_int_vec(5)
    rad_int1%i5a = rad_int1%i5a + rad_int_vec(6)
    rad_int1%i5b = rad_int1%i5b + rad_int_vec(7)
    rad_int1%i0  = rad_int1%i0  + rad_int_vec(8)
    rad_int1%i6b = rad_int1%i6b + rad_int_vec(9)

    int_tot%i1  = int_tot_vec(1) + rad_int1%i1
    int_tot%i2  = int_tot_vec(2) + rad_int1%i2
    int_tot%i3  = int_tot_vec(3) + rad_int1%i3
    int_tot%i4a = int_tot_vec(4) + rad_int1%i4a
    int_tot%i4b = int_tot_vec(5) + rad_int1%i4b
    int_tot%i5a = int_tot_vec(6) + rad_int1%i5a
    int_tot%i5b = int_tot_vec(7) + rad_int1%i5b
    int_tot%i0  = int_tot_vec(8) + rad_int1%i0
    int_tot%i6b = int_tot_vec(9) + rad_int1%i6b

  endif

  if (converged) return

end do

! We should not be here

call out_io (s_warn$, r_name, 'Note: Radiation Integral is not converging \es12.3\ ', 'For element: ' // &
                              trim(ele%name) // ' ' // ele_loc_name(ele, parens = '()'), r_array = [d_max])

end subroutine qromb_rad_int

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

subroutine propagate_part_way (orb_start, param, pt, info, z_here, runt)

implicit none

type (coord_struct) orb, orb_start, orb0, orb1, orb_end, orb_end1
type (ele_struct), pointer :: ele0, ele, field_ele
type (ele_struct) :: runt, ele_end
type (twiss_struct) a0, b0, a1, b1, x0, y0, x1, y1
type (rad_int_info_struct) info
type (rad_int_track_point_struct) pt, pt0, pt1
type (lat_param_struct) param

real(rp) z_here, v(4,4), v_inv(4,4), s1, s2, error
real(rp) f0, f1, del_z, c, s, x, y
real(rp) dz, z1
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx), dk(2,2), kx, ky
real(rp) eta_a0(4), eta_a1(4), eta_b0(4), eta_b1(4)

integer i0, i1, tm_saved, m6cm_saved
integer i, ix, ip, j_loop, n_pt, n, n1, n2, ix_pole_max

logical is_special_wiggler

! Init

ele0 => info%branch%ele(info%ele%ix_ele-1)
ele  => info%ele

if (runt%ix_ele /= info%ele%ix_ele) then
  runt = ele
  runt%ix_ele = info%ele%ix_ele
endif

!--------------------------------------
! With caching
! Remember that here the start and stop coords, etc. are in the local ref frame.

if (associated(info%cache_ele)) then

  ! Find cached point info near present z position.
  ! Note: There may be a jump at the element end due to fringe fields. 
  ! We want to exclude any points outside the fringe.

  n_pt = info%cache_ele%n_pt
  i0 = bracket_index(z_here, info%cache_ele%pt(0:n_pt)%s_body, 0)
  i0 = min(i0, n_pt-1)
  ! Downstream fringe if s_body is same for n_pt-1 and n_pt.
  ! Note: bracket_index is such that upstream fringe does not have to be checked.
  if (i0 == n_pt-1 .and. info%cache_ele%pt(n_pt-1)%s_body == info%cache_ele%pt(n_pt)%s_body) i0 = n_pt - 2
  i1 = i0 + 1
  pt0 = info%cache_ele%pt(i0)
  pt1 = info%cache_ele%pt(i1)
  del_z = pt1%s_body - pt0%s_body 
  if (del_z == 0) then
    f1 = 1
    f0 = 0
  else
    f1 = (z_here - pt0%s_body) / del_z 
    f0 = 1 - f1
  endif

  orb0%vec = matmul(pt0%mat6, orb_start%vec) + pt0%vec0
  orb1%vec = matmul(pt1%mat6, orb_start%vec) + pt1%vec0

  ! Interpolate information

  pt%dgx_dx = f0 * pt0%dgx_dx + f1 * pt1%dgx_dx
  pt%dgx_dy = f0 * pt0%dgx_dy + f1 * pt1%dgx_dy
  pt%dgy_dx = f0 * pt0%dgy_dx + f1 * pt1%dgy_dx
  pt%dgy_dy = f0 * pt0%dgy_dy + f1 * pt1%dgy_dy

  info%g_x  = f0 * (pt0%g_x0 + pt0%dgx_dx * (orb0%vec(1) - pt0%ref_orb_out%vec(1)) + &
                               pt0%dgx_dy * (orb0%vec(3) - pt0%ref_orb_out%vec(3))) + &
              f1 * (pt1%g_x0 + pt1%dgx_dx * (orb1%vec(1) - pt1%ref_orb_out%vec(1)) + &
                               pt1%dgx_dy * (orb1%vec(3) - pt1%ref_orb_out%vec(3)))
  info%g_y  = f0 * (pt0%g_y0 + pt0%dgy_dx * (orb0%vec(1) - pt0%ref_orb_out%vec(1)) + &
                               pt0%dgy_dy * (orb0%vec(3) - pt0%ref_orb_out%vec(3))) + &
              f1 * (pt1%g_y0 + pt1%dgy_dx * (orb1%vec(1) - pt1%ref_orb_out%vec(1)) + &
                               pt1%dgy_dy * (orb1%vec(3) - pt1%ref_orb_out%vec(3)))
               
  info%dg2_x = 2 * (info%g_x * pt%dgx_dx + info%g_y * pt%dgy_dx)
  info%dg2_y = 2 * (info%g_x * pt%dgx_dy + info%g_y * pt%dgy_dy) 
  info%g2 = info%g_x**2 + info%g_y**2
  info%g  = sqrt(info%g2)

  ! Interpolate the dispersions

  runt%mat6 = pt0%mat6
  runt%vec0 = pt0%vec0
  runt%map_ref_orb_in  = pt0%ref_orb_in
  runt%map_ref_orb_out = pt0%ref_orb_out

  call twiss_propagate1 (ele0, runt)
  a0 = runt%a; b0 = runt%b

  call make_v_mats (runt, v, v_inv)
  eta_a0 = matmul(v, [runt%a%eta, runt%a%etap, 0.0_rp,   0.0_rp    ])
  eta_b0 = matmul(v, [0.0_rp,   0.0_rp,    runt%b%eta, runt%b%etap ])

  !

  runt%mat6 = pt1%mat6
  runt%vec0 = pt1%vec0
  runt%map_ref_orb_in  = pt1%ref_orb_in
  runt%map_ref_orb_out = pt1%ref_orb_out

  call twiss_propagate1 (ele0, runt)
  a1 = runt%a; b1 = runt%b

  call make_v_mats (runt, v, v_inv)
  eta_a1 = matmul(v, [runt%a%eta, runt%a%etap, 0.0_rp,   0.0_rp    ])
  eta_b1 = matmul(v, [0.0_rp,   0.0_rp,    runt%b%eta, runt%b%etap ])

  !

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
! No caching...
! Note: calc_wiggler_g_params assumes that orb is lab (not element) coords.

dz = 1d-3
z1 = z_here + dz
if (z1 > ele%value(l$)) then
  z_here = max(0.0_rp, z_here - dz)
  z1 = min(ele%value(l$), z_here + dz)
endif

! If wiggler/undulator with Taylor tracking then switch to symp_lie_bmad tracking 

field_ele => pointer_to_field_ele(ele, 1)
is_special_wiggler = ((ele%key == wiggler$ .or. ele%key == undulator$) .and. ele%tracking_method == taylor$ .and. &
                                  (field_ele%field_calc == planar_model$ .or. field_ele%field_calc == helical_model$))

if (is_special_wiggler) then
  tm_saved = ele%tracking_method  
  m6cm_saved = ele%mat6_calc_method  
  ele%tracking_method = symp_lie_bmad$
  ele%mat6_calc_method = symp_lie_bmad$
endif

call twiss_and_track_intra_ele (ele, info%branch%param, 0.0_rp, z_here, .true., .false., orb_start, orb_end, ele0, ele_end)

info%a = ele_end%a
info%b = ele_end%b

if (is_special_wiggler) then
  ele%tracking_method  = tm_saved 
  ele%mat6_calc_method = m6cm_saved 
endif

!

call make_v_mats (ele_end, v, v_inv)

info%eta_a = matmul(v, [info%a%eta, info%a%etap, 0.0_rp,   0.0_rp ])
info%eta_b = matmul(v, [0.0_rp,   0.0_rp,    info%b%eta, info%b%etap ])

is_special_wiggler = ((ele%key == wiggler$ .or. ele%key == undulator$) .and. ele%tracking_method /= custom$)

! bmad_standard will not properly do partial tracking through a planar or helical wiggler so
! use special calculation

if (is_special_wiggler) then
  call calc_wiggler_g_params (ele, info%branch%param, z_here, orb_end, pt, info)
else
  call twiss_and_track_intra_ele (ele, info%branch%param, z_here, z1, .false., .false., orb_end, orb_end1)
  info%g_x = pt%g_x0 - (orb_end1%vec(2) - orb_end%vec(2)) / (z1 - z_here)
  info%g_y = pt%g_y0 - (orb_end1%vec(4) - orb_end%vec(4)) / (z1 - z_here)
  info%dg2_x = 0
  info%dg2_y = 0
endif

info%g2 = info%g_x**2 + info%g_y**2
info%g = sqrt(info%g2)

! Add in multipole gradient

if (.not. is_special_wiggler) then
  call multipole_ele_to_ab (ele, .true., ix_pole_max, a_pole, b_pole, include_kicks = include_kicks$)
  do ip = 0, ix_pole_max
    if (a_pole(ip) == 0 .and. b_pole(ip) == 0) cycle
    call ab_multipole_kick (a_pole(ip), b_pole(ip), ip, param%particle, ele%orientation, orb_end, kx, ky, dk)
    info%dg2_x = info%dg2_x - 2 * (info%g_x * dk(1,1) + info%g_y * dk(2,1)) / ele%value(l$)
    info%dg2_y = info%dg2_y - 2 * (info%g_x * dk(1,2) + info%g_y * dk(2,2)) / ele%value(l$)
  enddo
endif

end subroutine propagate_part_way

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine calc_wiggler_g_params (ele, param, s_rel, orb, pt, info)

implicit none

type (coord_struct) orb
type (lat_param_struct) param
type (rad_int_track_point_struct) pt
type (rad_int_info_struct), optional :: info
type (ele_struct) ele

real(rp) s_rel, g(3), dg(3,3)
real(rp) f0, k_z, sinh_x, cosh_x, sinh_y, cosh_y, sin_z, cos_z

! Note: Using lab (not element) coords here.

call g_bending_strength_from_em_field (ele, param, s_rel, orb, .false., g, dg)

pt%g_x0 = g(1)
pt%g_y0 = g(2)
pt%dgx_dx = dg(1,1)
pt%dgx_dy = dg(1,2)
pt%dgy_dx = dg(2,1)
pt%dgy_dy = dg(2,2)

if (present(info)) then
  info%g_x = g(1)
  info%g_y = g(2)
endif

end subroutine calc_wiggler_g_params

end module
