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

  type g_cache_struct
    real(rp) g
    real(rp) g2
    real(rp) g_x, g_y
    real(rp) dg2_x, dg2_y
  end type

  type ele_cache_struct
    type (g_cache_struct), allocatable :: v(:)
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
    real(rp), pointer :: i1_(:) => null()
    real(rp), pointer :: i2_(:) => null()
    real(rp), pointer :: i3_(:) => null()
    real(rp), pointer :: i4a_(:) => null()
    real(rp), pointer :: i4b_(:) => null()
    real(rp), pointer :: i5a_(:) => null()
    real(rp), pointer :: i5b_(:) => null()
    type (ring_struct), pointer :: ring
    type (ele_struct), pointer :: ele0, ele
    type (ele_struct) runt
    type (coord_struct), pointer :: orb0, orb1
    type (runge_kutta_com_struct) :: rk_track(0:6)
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
! Function qromb_rad_int (func, a, b, sum, which_int, &
!                                      eps_int, eps_sum) result (rad_int)
!
! Function to do integration using Romberg's method.
! This is a modified version of QROMB from Num. Rec.
! See the Num. Rec. book for further details
!-

function qromb_rad_int (func, a, b, sum, which_int, &
                                      eps_int, eps_sum) result (rad_int)

  use precision_def
  use nrtype; use nrutil, only : nrerror
  use nr, only : polint,trapzd

  implicit none

  interface
    function func(x)
    use precision_def
    real(rp), dimension(:), intent(in) :: x
    real(rp), dimension(size(x)) :: func
    end function func
  end interface

  integer(i4b), parameter :: jmax = 14, k = 5, km = k-1
  integer(i4b) :: j

  real(rp), intent(in), optional :: eps_int, eps_sum
  real(rp), intent(in) :: a, b, sum
  real(rp) :: rad_int
  real(rp) :: dqromb
  real(rp) :: eps_i, eps_s
  real(rp) :: h(0:jmax), s(0:jmax)

  character*(*) which_int

!

  eps_i = 1e-4
  if (present(eps_int)) eps_i = eps_int

  eps_s = 1e-6
  if (present(eps_sum)) eps_s = eps_sum

  h(0) = 4.0
  s(0) = 0.0

  do j = 1, jmax

    s(j) = s(j-1)
    h(j) = h(j-1) / 4
    call trapzd(func, a, b, s(j), j)

    if (ric%ele%key == wiggler$ .and. j < 4) cycle

    if (j >=  k) then
      call polint(h(j-km:j), s(j-km:j), 0.0_rp, rad_int, dqromb)
      if (abs(dqromb) <= eps_i * abs(rad_int) + eps_s * abs(sum)) return
    elseif (j >= 3) then
      call polint(h(1:j), s(1:j), 0.0_rp, rad_int, dqromb)
      if (abs(dqromb) <= eps_i * abs(rad_int) + eps_s * abs(sum)) return
    end if

  end do

  print *, &
      'Warning in QROMB_RAD_INT: Radiation integral is not converging: ', &
                                                            which_int
  print '(a, 1p3e11.2)', '     Error:', dqromb, rad_int, sum
  print *, '     For element: ', ric%ele%name

end function qromb_rad_int

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

function  eval_i1 (s_vec)

  implicit none

  real(rp), intent(in) :: s_vec(:)
  real(rp), dimension(size(s_vec)) :: eval_i1

  integer i

!                      
                                         
  do i = 1, size(s_vec)
    call propagate_part_way (s_vec(i))
    eval_i1(i) = ric%g_x * (ric%eta_a(1) + ric%eta_b(1)) + &
                 ric%g_y * (ric%eta_a(3) + ric%eta_b(3))
  enddo

end function
  
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

function  eval_i2 (s_vec)

  implicit none

  real(rp), intent(in) :: s_vec(:)
  real(rp), dimension(size(s_vec)) :: eval_i2

  integer i

!                      
                                         
  do i = 1, size(s_vec)
    call propagate_part_way (s_vec(i))
    eval_i2(i) = ric%g2
  enddo

end function
  
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

function  eval_i3 (s_vec)

  implicit none

  real(rp), intent(in) :: s_vec(:)
  real(rp), dimension(size(s_vec)) :: eval_i3

  integer i

!                      
                                         
  do i = 1, size(s_vec)
    call propagate_part_way (s_vec(i))
    eval_i3(i) = ric%g2 * ric%g
  enddo

end function
  
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

function  eval_i4a (s_vec)

  implicit none

  real(rp), intent(in) :: s_vec(:)
  real(rp), dimension(size(s_vec)) :: eval_i4a

  integer i

!

  do i = 1, size(s_vec)
    call propagate_part_way (s_vec(i))
    eval_i4a(i) = ric%g2 * (ric%g_x * ric%eta_a(1) + ric%g_y * ric%eta_a(3)) + &
                           (ric%dg2_x * ric%eta_a(1) + ric%dg2_y * ric%eta_a(3))
  enddo

end function
  
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

function  eval_i4b (s_vec)

  implicit none

  real(rp), intent(in) :: s_vec(:)
  real(rp), dimension(size(s_vec)) :: eval_i4b

  integer i

!

  do i = 1, size(s_vec)
    call propagate_part_way (s_vec(i))
    eval_i4b(i) = ric%g2 * (ric%g_x * ric%eta_b(1) + ric%g_y * ric%eta_b(3)) + &
                           (ric%dg2_x * ric%eta_b(1) + ric%dg2_y * ric%eta_b(3))
  enddo

end function
  
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

function  eval_i5a (s_vec)

  implicit none

  real(rp), intent(in) :: s_vec(:)
  real(rp), dimension(size(s_vec)) :: eval_i5a

!

  integer i

  do i = 1, size(s_vec)
    call propagate_part_way (s_vec(i))
    eval_i5a(i) = ric%g2 * ric%g * (ric%runt%x%gamma * ric%runt%x%eta**2 + &
                    2 * ric%runt%x%alpha * ric%runt%x%eta * ric%runt%x%etap + &
                    ric%runt%x%beta * ric%runt%x%etap**2)
  enddo             

end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

function  eval_i5b (s_vec)

  implicit none

  real(rp), intent(in) :: s_vec(:)
  real(rp), dimension(size(s_vec)) :: eval_i5b

  integer i

!

  do i = 1, size(s_vec)
    call propagate_part_way (s_vec(i))
    eval_i5b(i) = ric%g2 * ric%g * (ric%runt%y%gamma * ric%runt%y%eta**2 + &
                    2 * ric%runt%y%alpha * ric%runt%y%eta * ric%runt%y%etap + &
                    ric%runt%y%beta * ric%runt%y%etap**2)
  enddo

end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

subroutine transfer_rk_track (rk1, rk2)

  implicit none

  type (runge_kutta_com_struct) rk1, rk2

  integer n

!

  n = size(rk1%s)

  if (associated(rk2%s)) then
    if (size(rk2%s) < n) then
      deallocate (rk2%s, rk2%orb)
      allocate (rk2%s(n), rk2%orb(n))
    endif
  else
    allocate (rk2%s(n), rk2%orb(n))    
  endif

  n = rk1%n_pts

  rk2%n_pts    = rk1%n_pts
  rk2%s(1:n)   = rk1%s(1:n)
  rk2%orb(1:n) = rk1%orb(1:n)

end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! Subroutine bracket_index (s_, s, ix)
!
! Subroutine to find the index ix so that s_(ix) <= s < s_(ix+1).
! If s <  s_(1) then ix = 0
! If s >= s_(n) then ix = n  [where n = size(s_)]
!
! This routine assumes that s_ is in assending order.
!
! Input:
!   s_(:) -- Real(rp): Sequence of real numbers.
!   s     -- Real(rp): Number to bracket.
!
! Output:
!   ix    -- Integer: Index so that s_(ix) <= s < s_(ix+1).
!-

subroutine bracket_index (s_, s, ix)

  implicit none

  real(rp) s_(:), s

  integer i, ix, n, n1, n2, n3

!

  n = size(s_)

  if (s < s_(1)) then
    ix = 0
    return
  endif

  if (s >= s_(n)) then
    ix = n
    return
  endif

!

  n1 = 1
  n3 = n

  do

    if (n3 == n1 + 1) then
      ix = n1
      return
    endif

    n2 = (n1 + n3) / 2

    if (s < s_(n2)) then
      n3 = n2
    else
      n1 = n2
    endif

  enddo

end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

subroutine propagate_part_way (s)

  implicit none

  type (coord_struct) orb, orb_0

  real(rp) s, v(4,4), v_inv(4,4), s1, s2, error, f0, f1

  integer i, ix, n_pts

! exact calc

  if (ric%ele%exact_rad_int_calc) then

    do i = 0, 6
      n_pts = ric%rk_track(i)%n_pts
      call bracket_index (ric%rk_track(i)%s(1:n_pts), s, ix)

      if (ix == n_pts) then
        orb = ric%rk_track(i)%orb(n_pts)
      else
        s1 = s - ric%rk_track(i)%s(ix)
        s2 = ric%rk_track(i)%s(ix+1) - s
        orb%vec = (s2 * ric%rk_track(i)%orb(ix)%vec + &
                s1 * ric%rk_track(i)%orb(ix+1)%vec) / (s1 + s2)
      endif

      if (i == 0) then
        orb_0 = orb
        call calc_g_params (s, orb)
      else
        ric%runt%mat6(1:6, i) = (orb%vec - orb_0%vec) / ric%d_orb%vec(i)
      endif
    enddo

    call mat_symp_check (ric%runt%mat6, error)
    call mat_symplectify (ric%runt%mat6, ric%runt%mat6)

    call twiss_propagate1 (ric%ele0, ric%runt)

    call make_v_mats (ric%runt, v, v_inv)

    ric%eta_a = &
          matmul(v, (/ ric%runt%x%eta, ric%runt%x%etap, 0.0_rp, 0.0_rp /))
    ric%eta_b = &
          matmul(v, (/ 0.0_rp, 0.0_rp, ric%runt%y%eta, ric%runt%y%etap /))

    return
  endif

! non-exact wiggler calc

  if (ric%ele%key == wiggler$ .and. ric%ele%sub_key == map_type$) then

    f0 = (ric%ele%value(l$) - s) / ric%ele%value(l$)
    f1 = s / ric%ele%value(l$)

    orb%vec = ric%orb0%vec * f0 + ric%orb1%vec * f1
    call calc_g_params (s, orb)

    ric%eta_a = ric%eta_a0 * f0 + ric%eta_a1 * f1
    ric%eta_b = ric%eta_b0 * f0 + ric%eta_b1 * f1

    ric%runt%x%beta  = ric%ele0%x%beta  * f0 + ric%ele%x%beta  * f1
    ric%runt%x%alpha = ric%ele0%x%alpha * f0 + ric%ele%x%alpha * f1
    ric%runt%x%gamma = ric%ele0%x%gamma * f0 + ric%ele%x%gamma * f1
    ric%runt%x%eta   = ric%ele0%x%eta   * f0 + ric%ele%x%eta   * f1
    ric%runt%x%etap  = ric%ele0%x%etap  * f0 + ric%ele%x%etap  * f1

    ric%runt%y%beta  = ric%ele0%y%beta  * f0 + ric%ele%y%beta  * f1
    ric%runt%y%alpha = ric%ele0%y%alpha * f0 + ric%ele%y%alpha * f1
    ric%runt%y%gamma = ric%ele0%y%gamma * f0 + ric%ele%y%gamma * f1
    ric%runt%y%eta   = ric%ele0%y%eta   * f0 + ric%ele%y%eta   * f1
    ric%runt%y%etap  = ric%ele0%y%etap  * f0 + ric%ele%y%etap  * f1

    return
  endif

! non-exact calc

  if (s == 0) then
    ric%runt%x       = ric%ele0%x
    ric%runt%y       = ric%ele0%y
    ric%runt%c_mat   = ric%ele0%c_mat
    ric%runt%gamma_c = ric%ele0%gamma_c
    orb = ric%orb0
  elseif (s == ric%ele%value(l$)) then
    ric%runt = ric%ele
    orb = ric%orb1
  else
    ric%runt = ric%ele
    ric%runt%value(l$) = s
    if (ric%ele%key == sbend$) ric%runt%value(e2$) = 0
    call track1 (ric%orb0, ric%runt, ric%ring%param, orb)
    call make_mat6 (ric%runt, ric%ring%param, ric%orb0, orb)
    call twiss_propagate1 (ric%ele0, ric%runt)
  endif

  call make_v_mats (ric%runt, v, v_inv)

  ric%eta_a = &
      matmul(v, (/ ric%runt%x%eta, ric%runt%x%etap, 0.0_rp,   0.0_rp    /))
  ric%eta_b = &
      matmul(v, (/ 0.0_rp,   0.0_rp,    ric%runt%y%eta, ric%runt%y%etap /))

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

subroutine calc_g_params (s, orb)

  implicit none

  type (coord_struct) orb
  type (g_cache_struct) v0, v1

  real(rp) dk(3,3), s, ds
  real(rp) kick_0(6), f0, f1

  integer i0, i1

! Using the cache is faster if we have one.

  if (ric%use_cache) then
    ds = ric%cache_ele%ds
    i0 = int(s/ds)
    i1 = i0 + 1
    if (i1 > size(ric%cache_ele%v)) i1 = i0
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

  call derivs_bmad (ric%ele, ric%ring%param, s, orb%vec, kick_0, dk)

  ric%g_x = -kick_0(2)
  ric%g_y = -kick_0(4)
  ric%g2 = ric%g_x**2 + ric%g_y**2
  ric%g  = sqrt(ric%g2)

  ric%dg2_x = 2*kick_0(2)*dk(1,1) + 2*kick_0(4)*dk(2,1) 
  ric%dg2_y = 2*kick_0(2)*dk(1,2) + 2*kick_0(4)*dk(2,2) 

end subroutine

end module
