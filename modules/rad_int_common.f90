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

  type track_point_cache_struct
    type (coord_struct) orb
    real(rp) mat6(6,6)
    real(rp) vec0(6)
    real(rp) g, g2          ! bending strength (1/bending_radius)
    real(rp) g_x0, g_y0     ! components on the reference orb.
    real(rp) g_x, g_y       ! components in x-y plane
    real(rp) dg2_x, dg2_y   ! bending strength gradiant
    real(rp) k1, s1         ! quad and skew quad components
  end type

  type ele_cache_struct
    type (track_point_cache_struct), allocatable :: pt(:)
    real(rp) del_z
    integer ix_ele
  end type

  type rad_int_cache_struct
    type (ele_cache_struct), allocatable :: ele(:)
    integer, allocatable :: ix_ele(:)
    logical :: set = .false.   ! is being used?
  end type

! This structure stores the radiation integrals for the individual elements
! eta_a(4) is the a-mode dispersion in the lab frame.

  type rad_int_common_struct
    type (ring_struct), pointer :: ring
    type (coord_struct), pointer :: orb0, orb1
    type (rad_int_cache_struct) cache(10)         ! up to 10 separate caches
    type (ele_cache_struct), pointer :: cache_ele ! pointer to cache in use
    type (track_point_cache_struct) pt
    real(rp) eta_a(4), eta_b(4), eta_a0(4), eta_a1(4), eta_b0(4), eta_b1(4)
    real(rp), allocatable :: i1_(:) 
    real(rp), allocatable :: i2_(:) 
    real(rp), allocatable :: i3_(:) 
    real(rp), allocatable :: i4a_(:)
    real(rp), allocatable :: i4b_(:)
    real(rp), allocatable :: i5a_(:) 
    real(rp), allocatable :: i5b_(:) 
    real(rp), allocatable :: n_steps(:)      ! number of qromb steps needed
    real(rp) :: int_tot(7)
    logical use_cache
  end type

  type (rad_int_common_struct), target, save :: ric

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! Subroutine qromb_rad_int(ele0, ele, do_int, ir)
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

subroutine qromb_rad_int (ele0, ele, do_int, ir)

  use precision_def
  use nrtype
  use nr, only: polint

  implicit none

  type (ele_struct) ele0, ele
  type (ele_struct), save :: runt

  integer, parameter :: jmax = 14
  integer j, j0, n, n_pts, ir

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

  eps_int = 1e-4
  eps_sum = 1e-6

  ri(:)%h(0) = 4.0
  ri(:)%sum(0) = 0
  rad_int = 0
  
  ll = ele%value(l$)

  runt = ele
  if (runt%tracking_method  == taylor$) runt%tracking_method  = bmad_standard$
  if (runt%mat6_calc_method == taylor$) runt%mat6_calc_method = bmad_standard$


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
      call propagate_part_way (ele0, ele, runt, z_pos, j, n)
      i_sum(1) = i_sum(1) + ric%pt%g_x * (ric%eta_a(1) + ric%eta_b(1)) + &
                            ric%pt%g_y * (ric%eta_a(3) + ric%eta_b(3))
      i_sum(2) = i_sum(2) + ric%pt%g2
      i_sum(3) = i_sum(3) + ric%pt%g2 * ric%pt%g
      i_sum(4) = i_sum(4) + &
                ric%pt%g2 * (ric%pt%g_x * ric%eta_a(1) + ric%pt%g_y * ric%eta_a(3)) + &
                         (ric%pt%dg2_x * ric%eta_a(1) + ric%pt%dg2_y * ric%eta_a(3)) 
      i_sum(5) = i_sum(5) + &
                ric%pt%g2 * (ric%pt%g_x * ric%eta_b(1) + ric%pt%g_y * ric%eta_b(3)) + &
                         (ric%pt%dg2_x * ric%eta_b(1) + ric%pt%dg2_y * ric%eta_b(3))
      i_sum(6) = i_sum(6) + &
                    ric%pt%g2 * ric%pt%g * (runt%x%gamma * runt%x%eta**2 + &
                    2 * runt%x%alpha * runt%x%eta * runt%x%etap + &
                    runt%x%beta * runt%x%etap**2)
      i_sum(7) = i_sum(7) + &
                    ric%pt%g2 * ric%pt%g * (runt%y%gamma * runt%y%eta**2 + &
                    2 * runt%y%alpha * runt%y%eta * runt%y%etap + &
                    runt%y%beta * runt%y%etap**2)
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
      d0 = eps_int * abs(rad_int(n)) + eps_sum * abs(ric%int_tot(n))
      if (abs(dint) > d0)  complete = .false.
      if (d0 /= 0) d_max = abs(dint) / d0
    enddo

! If we have convergance or we are giving up (when j = jmax) then 
! stuff the results in the proper places.

    if (complete .or. j == jmax) then

      ric%n_steps(ir) = j

      ! Note that ric%i... may already contain a contribution from edge
      ! affects (Eg bend face angles) so add it on to rad_int(i)

      ric%i1_(ir)  = ric%i1_(ir)  + rad_int(1)
      ric%i2_(ir)  = ric%i2_(ir)  + rad_int(2)
      ric%i3_(ir)  = ric%i3_(ir)  + rad_int(3)
      ric%i4a_(ir) = ric%i4a_(ir) + rad_int(4)
      ric%i4b_(ir) = ric%i4b_(ir) + rad_int(5)
      ric%i5a_(ir) = ric%i5a_(ir) + rad_int(6)
      ric%i5b_(ir) = ric%i5b_(ir) + rad_int(7)

      ric%int_tot(1) = ric%int_tot(1) + ric%i1_(ir)
      ric%int_tot(2) = ric%int_tot(2) + ric%i2_(ir)
      ric%int_tot(3) = ric%int_tot(3) + ric%i3_(ir)
      ric%int_tot(4) = ric%int_tot(4) + ric%i4a_(ir)
      ric%int_tot(5) = ric%int_tot(5) + ric%i4b_(ir)
      ric%int_tot(6) = ric%int_tot(6) + ric%i5a_(ir)
      ric%int_tot(7) = ric%int_tot(7) + ric%i5b_(ir)

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
!+
! Subroutine bracket_index (s_arr, i_min, i_max, s, ix)
!
! Subroutine to find the index ix so that s_arr(ix) <= s < s_arr(ix+1).
! If s <  s_arr(i_min) then ix = i_min - 1
! If s >= s_arr(i_max) then ix = i_max  [where i_max = ubound(s_arr)]
!
! This routine assumes that s_arr is in assending order.
!
! Input:
!   s_arr(i_min:) -- Real(rp): Sequence of real numbers.
!   i_min         -- Integer: lower bound of s_arr
!   i_max         -- Integer: upper bound of s_arr
!   s             -- Real(rp): Number to bracket.
!
! Output:
!   ix    -- Integer: Index so that s_arr(ix) <= s < s_arr(ix+1).
!-

subroutine bracket_index (s_arr, i_min, i_max, s, ix)

  implicit none

  integer i_min, i_max
  real(rp) s_arr(i_min:), s

  integer ix, n1, n2, n3

!

  if (s < s_arr(i_min)) then
    ix = i_min - 1
    return
  endif

  if (s >= s_arr(i_max)) then
    ix = i_max
    return
  endif

!

  n1 = i_min
  n3 = i_max

  do

    if (n3 == n1 + 1) then
      ix = n1
      return
    endif

    n2 = (n1 + n3) / 2

    if (s < s_arr(n2)) then
      n3 = n2
    else
      n1 = n2
    endif

  enddo

end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

subroutine propagate_part_way (ele0, ele, runt, z_here, j_loop, n_pt)

  implicit none

  type (coord_struct) orb, orb_0
  type (ele_struct) ele0, ele
  type (ele_struct) runt, e0, e1
  type (track_point_cache_struct) pt0, pt1

  real(rp) z_here, v(4,4), v_inv(4,4), s1, s2, error
  real(rp) f0, f1, del_z, c, s, x, y

  integer i0, i1
  integer i, ix, j_loop, n_pt, n, n1, n2

!--------------------------------------
! For caching or map_type wiggers

  if (ric%use_cache) then

    del_z = ric%cache_ele%del_z
    i0 = int(z_here/del_z)
    i1 = i0 + 1
    if (i1 > ubound(ric%cache_ele%pt, 1)) i1 = i0  ! can happen with roundoff
    f1 = (z_here - del_z*i0) / del_z 
    f0 = 1 - f1
    pt0 = ric%cache_ele%pt(i0)
    pt1 = ric%cache_ele%pt(i1)

    orb%vec = ric%orb0%vec * f0 + ric%orb1%vec * f1

    ! g factors have been already calculated for non-wiggler elements.
    if (ele%key == wiggler$ .and. ele%sub_key == map_type$) then
      ric%pt%g      = f0 * pt0%g     + f1 * pt1%g
      ric%pt%g2     = f0 * pt0%g2    + f1 * pt1%g2
      ric%pt%g_x    = f0 * pt0%g_x   + f1 * pt1%g_x
      ric%pt%g_y    = f0 * pt0%g_y   + f1 * pt1%g_y
      ric%pt%dg2_x  = f0 * pt0%dg2_x + f1 * pt1%dg2_x
      ric%pt%dg2_y  = f0 * pt0%dg2_y + f1 * pt1%dg2_y
    else
      ric%pt%g_x = ric%pt%g_x0 + orb%vec(1) * ric%pt%k1 + orb%vec(3) * ric%pt%s1
      ric%pt%g_y = ric%pt%g_y0 - orb%vec(3) * ric%pt%k1 + orb%vec(1) * ric%pt%s1
      ric%pt%dg2_x = 2 * (ric%pt%g_x * ric%pt%k1 + ric%pt%g_y * ric%pt%s1)
      ric%pt%dg2_y = 2 * (ric%pt%g_x * ric%pt%s1 - ric%pt%g_y * ric%pt%k1) 
      ric%pt%g2 = ric%pt%g_x**2 + ric%pt%g_y**2
      ric%pt%g = sqrt(ric%pt%g2)
    endif

    if (.not. ele%map_with_offsets) then
      c = cos(ele%value(tilt_tot$)); s = sin(ele%value(tilt_tot$)) 
      x = ric%pt%g_x; y = ric%pt%g_y
      ric%pt%g_x = c * x - s * y
      ric%pt%g_y = s * x + c * y
      x = ric%pt%dg2_x; y = ric%pt%dg2_y
      ric%pt%dg2_x = c * x - s * y
      ric%pt%dg2_y = s * x + c * y
    endif

    runt%mat6 = pt0%mat6
    runt%vec0 = pt0%vec0
    if (.not. ele%map_with_offsets) call mat6_add_offsets (runt)
    call twiss_propagate1 (ele0, runt)
    e0%x = runt%x; e0%y = runt%y
    call make_v_mats (runt, v, v_inv)
    ric%eta_a0 = &
      matmul(v, (/ runt%x%eta, runt%x%etap, 0.0_rp,   0.0_rp    /))
    ric%eta_b0 = &
      matmul(v, (/ 0.0_rp,   0.0_rp,    runt%y%eta, runt%y%etap /))

    runt%mat6 = pt1%mat6
    runt%vec0 = pt1%vec0
    if (.not. ele%map_with_offsets) call mat6_add_offsets (runt)
    call twiss_propagate1 (ele0, runt)
    e1%x = runt%x; e1%y = runt%y
    call make_v_mats (runt, v, v_inv)
    ric%eta_a1 = &
      matmul(v, (/ runt%x%eta, runt%x%etap, 0.0_rp,   0.0_rp    /))
    ric%eta_b1 = &
      matmul(v, (/ 0.0_rp,   0.0_rp,    runt%y%eta, runt%y%etap /))

    runt%x%beta  = e0%x%beta  * f0 + e1%x%beta  * f1
    runt%x%alpha = e0%x%alpha * f0 + e1%x%alpha * f1
    runt%x%gamma = e0%x%gamma * f0 + e1%x%gamma * f1
    runt%x%eta   = e0%x%eta   * f0 + e1%x%eta   * f1
    runt%x%etap  = e0%x%etap  * f0 + e1%x%etap  * f1

    runt%y%beta  = e0%y%beta  * f0 + e1%y%beta  * f1
    runt%y%alpha = e0%y%alpha * f0 + e1%y%alpha * f1
    runt%y%gamma = e0%y%gamma * f0 + e1%y%gamma * f1
    runt%y%eta   = e0%y%eta   * f0 + e1%y%eta   * f1
    runt%y%etap  = e0%y%etap  * f0 + e1%y%etap  * f1

    ric%eta_a = ric%eta_a0 * f0 + ric%eta_a1 * f1
    ric%eta_b = ric%eta_b0 * f0 + ric%eta_b1 * f1

    return

  ! map type wiggler with no caching case
  elseif (ele%key == wiggler$ .and. ele%sub_key == map_type$) then 

    f0 = (ele%value(l$) - z_here) / ele%value(l$)
    f1 = z_here / ele%value(l$)

    orb%vec = ric%orb0%vec * f0 + ric%orb1%vec * f1
    call calc_wiggler_g_params (ele, z_here, orb, ric%pt)

    runt%x%beta  = ele0%x%beta  * f0 + ele%x%beta  * f1
    runt%x%alpha = ele0%x%alpha * f0 + ele%x%alpha * f1
    runt%x%gamma = ele0%x%gamma * f0 + ele%x%gamma * f1
    runt%x%eta   = ele0%x%eta   * f0 + ele%x%eta   * f1
    runt%x%etap  = ele0%x%etap  * f0 + ele%x%etap  * f1

    runt%y%beta  = ele0%y%beta  * f0 + ele%y%beta  * f1
    runt%y%alpha = ele0%y%alpha * f0 + ele%y%alpha * f1
    runt%y%gamma = ele0%y%gamma * f0 + ele%y%gamma * f1
    runt%y%eta   = ele0%y%eta   * f0 + ele%y%eta   * f1
    runt%y%etap  = ele0%y%etap  * f0 + ele%y%etap  * f1

    ric%eta_a = ric%eta_a0 * f0 + ric%eta_a1 * f1
    ric%eta_b = ric%eta_b0 * f0 + ric%eta_b1 * f1

    return

  endif

!--------------------------------------
! Not a map type wiggler and No caching involved case.

  if (j_loop == 1 .and. n_pt == 1) then  ! z_here = 0
    runt%x       = ele0%x
    runt%y       = ele0%y
    runt%c_mat   = ele0%c_mat
    runt%gamma_c = ele0%gamma_c
    orb = ric%orb0

  elseif (j_loop == 1 .and. n_pt == 2) then  ! z_here = l$
    runt%x       = ele%x
    runt%y       = ele%y
    runt%c_mat   = ele%c_mat
    runt%gamma_c = ele%gamma_c
    orb = ric%orb1

  else
    runt%value(l$) = z_here
    if (ele%key == sbend$) runt%value(e2$) = 0
    call track1 (ric%orb0, runt, ric%ring%param, orb)
    call make_mat6 (runt, ric%ring%param, ric%orb0, orb, .true.)
    call twiss_propagate1 (ele0, runt)

  endif

!

  call make_v_mats (runt, v, v_inv)

  ric%eta_a = &
      matmul(v, (/ runt%x%eta, runt%x%etap, 0.0_rp,   0.0_rp    /))
  ric%eta_b = &
      matmul(v, (/ 0.0_rp,   0.0_rp,    runt%y%eta, runt%y%etap /))

  ric%pt%g_x = ric%pt%g_x0 + orb%vec(1) * ric%pt%k1 + orb%vec(3) * ric%pt%s1
  ric%pt%g_y = ric%pt%g_y0 - orb%vec(3) * ric%pt%k1 + orb%vec(1) * ric%pt%s1

  ric%pt%dg2_x = 2 * (ric%pt%g_x * ric%pt%k1 + ric%pt%g_y * ric%pt%s1)
  ric%pt%dg2_y = 2 * (ric%pt%g_x * ric%pt%s1 - ric%pt%g_y * ric%pt%k1) 

  ric%pt%g2 = ric%pt%g_x**2 + ric%pt%g_y**2
  ric%pt%g = sqrt(ric%pt%g2)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine calc_wiggler_g_params (ele, z, orb, pt)

  implicit none

  type (coord_struct) orb
  type (track_point_cache_struct) pt
  type (ele_struct) ele

  real(rp) dk(3,3), z
  real(rp) kick_0(6)

! Standard non-cache calc.

  call derivs_bmad (ele, ric%ring%param, z, orb%vec, kick_0, dk)

  pt%g_x = -kick_0(2)
  pt%g_y = -kick_0(4)
  pt%g_x0 = -kick_0(2)
  pt%g_y0 = -kick_0(4)
  pt%g2 = pt%g_x**2 + pt%g_y**2
  pt%g  = sqrt(pt%g2)
  pt%dg2_x = 2*kick_0(2)*dk(1,1) + 2*kick_0(4)*dk(2,1) 
  pt%dg2_y = 2*kick_0(2)*dk(1,2) + 2*kick_0(4)*dk(2,2) 

end subroutine

end module
