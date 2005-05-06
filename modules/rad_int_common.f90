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

! The "cache" is for saving values for g, etc through a wiggler to speed
! up the calculation

  type wig_cache_struct
    real(rp) g
    real(rp) g2
    real(rp) g_x, g_y
    real(rp) dg2_x, dg2_y
  end type

  type quad_bend_cache_struct
    real(rp) mat6(6,6)
    type (coord_struct) orb
  end type

  type ele_cache_struct
    type (wig_cache_struct), allocatable :: v(:)
    type (quad_bend_cache_struct) qb(3)
    real(rp) ds
    integer ix_ele
  end type

  type rad_int_cache_struct
    type (ele_cache_struct), allocatable :: ele(:)
    logical :: set = .false.
  end type

! This structure stores the radiation integrals for the individual elements

  type rad_int_common_struct
    real(rp) g_x0, g_y0, k1, s1
    real(rp) eta_a(4), eta_b(4), eta_a0(4), eta_a1(4), eta_b0(4), eta_b1(4)
    real(rp) g, g2, g_x, g_y, dg2_x, dg2_y 
    real(rp), allocatable :: i1_(:) 
    real(rp), allocatable :: i2_(:) 
    real(rp), allocatable :: i3_(:) 
    real(rp), allocatable :: i4a_(:)
    real(rp), allocatable :: i4b_(:)
    real(rp), allocatable :: i5a_(:) 
    real(rp), allocatable :: i5b_(:) 
    real(rp), allocatable :: n_steps(:)      ! number of qromb steps needed
    real(rp) :: int_tot(7)
    type (ring_struct), pointer :: ring
    type (coord_struct), pointer :: orb0, orb1
    type (track_struct) :: track(0:6)
    type (coord_struct) d_orb
    type (rad_int_cache_struct) cache(10)
    type (ele_cache_struct), pointer :: cache_ele
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
! See the Num. Rec. book for further details
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
  real(rp) :: ll, ds, s0, s_pos, dint, d0, d_max
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
! There are up to 7 integrals that are calculated:
!          I1, I2, I3, I4a, I4b, I5a, I5b
! ri(k) holds the info for the k^th integral.

  do j = 1, jmax

    ri(:)%h(j) = ri(:)%h(j-1) / 4

!---------------
! This is trapzd from Num. Rec.

    if (j == 1) then
      n_pts = 2
      ds = ll
      s0 = 0
    else
      n_pts = 2**(j-2)
      ds = ll / n_pts
      s0 = ds / 2
    endif

    i_sum = 0

    do n = 1, n_pts
      s_pos = s0 + (n-1) * ds
      call propagate_part_way (ele0, ele, runt, s_pos, j, n)
      i_sum(1) = i_sum(1) + ric%g_x * (ric%eta_a(1) + ric%eta_b(1)) + &
                            ric%g_y * (ric%eta_a(3) + ric%eta_b(3))
      i_sum(2) = i_sum(2) + ric%g2
      i_sum(3) = i_sum(3) + ric%g2 * ric%g
      i_sum(4) = i_sum(4) + &
                ric%g2 * (ric%g_x * ric%eta_a(1) + ric%g_y * ric%eta_a(3)) + &
                         (ric%dg2_x * ric%eta_a(1) + ric%dg2_y * ric%eta_a(3)) 
      i_sum(5) = i_sum(5) + &
                ric%g2 * (ric%g_x * ric%eta_b(1) + ric%g_y * ric%eta_b(3)) + &
                         (ric%dg2_x * ric%eta_b(1) + ric%dg2_y * ric%eta_b(3))
      i_sum(6) = i_sum(6) + &
                    ric%g2 * ric%g * (runt%x%gamma * runt%x%eta**2 + &
                    2 * runt%x%alpha * runt%x%eta * runt%x%etap + &
                    runt%x%beta * runt%x%etap**2)
      i_sum(7) = i_sum(7) + &
                    ric%g2 * ric%g * (runt%y%gamma * runt%y%eta**2 + &
                    2 * runt%y%alpha * runt%y%eta * runt%y%etap + &
                    runt%y%beta * runt%y%etap**2)
    enddo

    ri(:)%sum(j) = (ri(:)%sum(j-1) + ds * i_sum(:)) / 2

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

! If we have convergance or we are giving up (when j = jmax) then ...

    if (complete .or. j == jmax) then

      ric%n_steps(ir) = j

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
! Subroutine bracket_index (s_arr, i_min, s, ix)
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
!   s             -- Real(rp): Number to bracket.
!
! Output:
!   ix    -- Integer: Index so that s_arr(ix) <= s < s_arr(ix+1).
!-

subroutine bracket_index (s_arr, i_min, s, ix)

  implicit none

  integer i_min, i_max
  real(rp) s_arr(i_min:), s

  integer ix, n1, n2, n3

!

  i_max = ubound(s_arr, 1)

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

subroutine propagate_part_way (ele0, ele, runt, s, j_loop, n_pt)

  implicit none

  type (coord_struct) orb, orb_0
  type (ele_struct) ele0, ele
  type (ele_struct) runt

  real(rp) s, v(4,4), v_inv(4,4), s1, s2, error, f0, f1

  integer i, ix, j_loop, n_pt, n, n1, n2

! exact calc

  if (ele%exact_rad_int_calc) then

    do i = 0, 6
      n1 = lbound(ric%track(i)%pt, 1)
      n2 = ric%track(i)%n_pt
      call bracket_index (ric%track(i)%pt(:)%s, n1, s, ix)

      if (ix == n2) then
        orb = ric%track(i)%pt(n2)%orb
      else
        s1 = s - ric%track(i)%pt(ix)%s
        s2 = ric%track(i)%pt(ix+1)%s - s
        orb%vec = (s2 * ric%track(i)%pt(ix)%orb%vec + &
                s1 * ric%track(i)%pt(ix+1)%orb%vec) / (s1 + s2)
      endif

      if (i == 0) then
        orb_0 = orb
        call calc_g_params (ele, s, orb)
      else
        runt%mat6(1:6, i) = (orb%vec - orb_0%vec) / ric%d_orb%vec(i)
      endif
    enddo

    call mat_symp_check (runt%mat6, error)
    call mat_symplectify (runt%mat6, runt%mat6)

    call twiss_propagate1 (ele0, runt)

    call make_v_mats (runt, v, v_inv)

    ric%eta_a = &
          matmul(v, (/ runt%x%eta, runt%x%etap, 0.0_rp, 0.0_rp /))
    ric%eta_b = &
          matmul(v, (/ 0.0_rp, 0.0_rp, runt%y%eta, runt%y%etap /))

    return
  endif

! non-exact wiggler calc

  if (ele%key == wiggler$ .and. ele%sub_key == map_type$) then

    f0 = (ele%value(l$) - s) / ele%value(l$)
    f1 = s / ele%value(l$)

    orb%vec = ric%orb0%vec * f0 + ric%orb1%vec * f1
    call calc_g_params (ele, s, orb)

    ric%eta_a = ric%eta_a0 * f0 + ric%eta_a1 * f1
    ric%eta_b = ric%eta_b0 * f0 + ric%eta_b1 * f1

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

    return
  endif

! non-exact calc

  if (j_loop == 1 .and. n_pt == 1) then  ! s = 0
    runt%x       = ele0%x
    runt%y       = ele0%y
    runt%c_mat   = ele0%c_mat
    runt%gamma_c = ele0%gamma_c
    orb = ric%orb0
  elseif (j_loop == 1 .and. n_pt == 2) then  ! s = l$
    runt%x       = ele%x
    runt%y       = ele%y
    runt%c_mat   = ele%c_mat
    runt%gamma_c = ele%gamma_c
    orb = ric%orb1
  else
    if (ric%use_cache .and. (j_loop == 2 .or. j_loop == 3)) then
      n = j_loop + n_pt - 2
      runt%mat6 = ric%cache_ele%qb(n)%mat6
      orb = ric%cache_ele%qb(n)%orb
    else
      runt%value(l$) = s
      if (ele%key == sbend$) runt%value(e2$) = 0
      call track1 (ric%orb0, runt, ric%ring%param, orb)
      call make_mat6 (runt, ric%ring%param, ric%orb0, orb, .true.)
    endif
    call twiss_propagate1 (ele0, runt)
  endif

  call make_v_mats (runt, v, v_inv)

  ric%eta_a = &
      matmul(v, (/ runt%x%eta, runt%x%etap, 0.0_rp,   0.0_rp    /))
  ric%eta_b = &
      matmul(v, (/ 0.0_rp,   0.0_rp,    runt%y%eta, runt%y%etap /))

  ric%g_x = ric%g_x0 + orb%vec(1) * ric%k1 + orb%vec(3) * ric%s1
  ric%g_y = ric%g_y0 - orb%vec(3) * ric%k1 + orb%vec(1) * ric%s1
                   
  ric%dg2_x = 2 * (ric%g_x * ric%k1 + ric%g_y * ric%s1)
  ric%dg2_y = 2 * (ric%g_x * ric%s1 - ric%g_y * ric%k1) 

  ric%g2 = ric%g_x**2 + ric%g_y**2
  ric%g = sqrt(ric%g2)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine calc_g_params (ele, s, orb)

  implicit none

  type (coord_struct) orb
  type (wig_cache_struct) v0, v1
  type (ele_struct) ele

  real(rp) dk(3,3), s, ds
  real(rp) kick_0(6), f0, f1

  integer i0, i1

! Using the cache is faster if we have one.

  if (ric%use_cache) then
    ds = ric%cache_ele%ds
    i0 = int(s/ds)
    i1 = i0 + 1
    if (i1 > ubound(ric%cache_ele%v, 1)) i1 = i0  ! can happen with roundoff
    f1 = (s - ds*i0) / ds 
    f0 = 1 - f1
    v0 = ric%cache_ele%v(i0)
    v1 = ric%cache_ele%v(i1)
    ric%g      = f0 * v0%g     + f1 * v1%g
    ric%g2     = f0 * v0%g2    + f1 * v1%g2
    ric%g_x    = f0 * v0%g_x   + f1 * v1%g_x
    ric%g_y    = f0 * v0%g_y   + f1 * v1%g_y
    ric%dg2_x  = f0 * v0%dg2_x + f1 * v1%dg2_x
    ric%dg2_y  = f0 * v0%dg2_y + f1 * v1%dg2_y
    return
  endif

! Standard non-cache calc.

  call derivs_bmad (ele, ric%ring%param, s, orb%vec, kick_0, dk)

  ric%g_x = -kick_0(2)
  ric%g_y = -kick_0(4)
  ric%g2 = ric%g_x**2 + ric%g_y**2
  ric%g  = sqrt(ric%g2)

  ric%dg2_x = 2*kick_0(2)*dk(1,1) + 2*kick_0(4)*dk(2,1) 
  ric%dg2_y = 2*kick_0(2)*dk(1,2) + 2*kick_0(4)*dk(2,2) 

end subroutine

end module
