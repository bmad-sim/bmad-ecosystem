!+
! Module rad_int_common
!
! Module needed:
!   use rad_int_common
!-

!$Id$
!$Log$
!Revision 1.8  2003/05/02 15:44:35  dcs
!F90 standard conforming changes.
!
!Revision 1.7  2003/01/27 14:41:01  dcs
!bmad_version = 56
!
!Revision 1.5  2002/06/13 14:55:00  dcs
!Interfaced with FPP/PTC
!
!Revision 1.4  2002/02/23 20:32:32  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2001/10/02 18:50:26  rwh24
!Extended track_input_struct in bmad_struct.
!Fixed qromb_rad_int definition in rad_int_common.
!
!Revision 1.2  2001/09/27 18:32:13  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

module rad_int_common               

  use ptc_interface_mod
  use runge_kutta_mod

! i1_(:), etc. are for diagnostics

  type rad_int_common_struct
    real(rdef) g_x0, g_y0, k1, s1
    real(rdef) eta_a(4), eta_b(4), eta_a0(4), eta_a1(4), eta_b0(4), eta_b1(4)
    real(rdef) g, g2, g_x, g_y, dg2_x, dg2_y 
    real(rdef) i1_(n_ele_maxx), i2_(n_ele_maxx), i3_(n_ele_maxx)
    real(rdef) i4a_(n_ele_maxx), i4b_(n_ele_maxx)
    real(rdef) i5a_(n_ele_maxx), i5b_(n_ele_maxx)
    type (ring_struct), pointer :: ring
    type (ele_struct), pointer :: ele0, ele
    type (ele_struct) runt
    type (coord_struct), pointer :: orb0, orb1
    type (runge_kutta_com_struct) :: rk_track(0:6)
    type (coord_struct) d_orb
  end type

  type (rad_int_common_struct), save :: ric

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! This is a modified version of QROMB from Num. Rec.

contains

function qromb_rad_int (func, a, b, sum, which_int) result (rad_int)

  use precision_def
  use nrtype; use nrutil, only : nrerror
  use nr, only : polint,trapzd

  implicit none

  interface
    function func(x)
    use precision_def
    real(rdef), dimension(:), intent(in) :: x
    real(rdef), dimension(size(x)) :: func
    end function func
  end interface

  integer(i4b), parameter :: jmax = 14, k = 5, km = k-1
  integer(i4b) :: j

  real(rdef), intent(in) :: a, b, sum
  real(rdef) :: rad_int
  real(rdef) :: dqromb
  real(rdef), parameter :: eps = 1.0e-4_rdef, eps2 = 1.0e-6_rdef
  real(rdef) :: h(0:jmax), s(0:jmax)

  character*(*) which_int

!

  h(0) = 4.0
  s(0) = 0.0

  do j = 1, jmax

    s(j) = s(j-1)
    h(j) = h(j-1) / 4
    call trapzd(func, a, b, s(j), j)

    if (ric%ele%key == wiggler$ .and. j < 4) cycle

    if (j >=  k) then
      call polint(h(j-km:j), s(j-km:j), 0.0_rdef, rad_int, dqromb)
      if (abs(dqromb) <= eps * abs(rad_int) + eps2 * abs(sum)) return
    elseif (j >= 3) then
      call polint(h(1:j), s(1:j), 0.0_rdef, rad_int, dqromb)
      if (abs(dqromb) <= eps * abs(rad_int) + eps2 * abs(sum)) return
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

  real(rdef), intent(in) :: s_vec(:)
  real(rdef), dimension(size(s_vec)) :: eval_i1

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

  real(rdef), intent(in) :: s_vec(:)
  real(rdef), dimension(size(s_vec)) :: eval_i2

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

  real(rdef), intent(in) :: s_vec(:)
  real(rdef), dimension(size(s_vec)) :: eval_i3

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

  real(rdef), intent(in) :: s_vec(:)
  real(rdef), dimension(size(s_vec)) :: eval_i4a

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

  real(rdef), intent(in) :: s_vec(:)
  real(rdef), dimension(size(s_vec)) :: eval_i4b

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

  real(rdef), intent(in) :: s_vec(:)
  real(rdef), dimension(size(s_vec)) :: eval_i5a

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

  real(rdef), intent(in) :: s_vec(:)
  real(rdef), dimension(size(s_vec)) :: eval_i5b

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
!   s_(:) -- Real(rdef): Sequence of real numbers.
!   s     -- Real(rdef): Number to bracket.
!
! Output:
!   ix    -- Integer: Index so that s_(ix) <= s < s_(ix+1).
!-

subroutine bracket_index (s_, s, ix)

  implicit none

  real(rdef) s_(:), s

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

  real(rdef) s, v(4,4), v_inv(4,4), s1, s2, error, f0, f1
  real(rdef) here(6), kick_0(6), kick_x(6), kick_y(6)

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
        call calc_g_params (orb)
      else
        ric%runt%mat6(1:6, i) = (orb%vec - orb_0%vec) / ric%d_orb%vec(i)
      endif
    enddo

    call mat_symp_check (ric%runt%mat6, error)
    call mat_symplectify (ric%runt%mat6, ric%runt%mat6)

    call twiss_propagate1 (ric%ele0, ric%runt)

    call make_v_mats (ric%runt, v, v_inv)

    ric%eta_a = &
          matmul(v, (/ ric%runt%x%eta, ric%runt%x%etap, 0.0_rdef, 0.0_rdef /))
    ric%eta_b = &
          matmul(v, (/ 0.0_rdef, 0.0_rdef, ric%runt%y%eta, ric%runt%y%etap /))

    return
  endif

! non-exact wiggler calc

  if (ric%ele%key == wiggler$ .and. ric%ele%sub_key == map_type$) then
!    n_pts = ric%rk_track(i)%n_pts
!    call bracket_index (ric%rk_track(i)%s(1:n_pts), s, ix)

!    if (ix == n_pts) then
!      orb = ric%rk_track(i)%orb(n_pts)
!    else
!      s1 = s - ric%rk_track(i)%s(ix)
!      s2 = ric%rk_track(i)%s(ix+1) - s
!      orb%vec = (s2 * ric%rk_track(i)%orb(ix)%vec + &
!                s1 * ric%rk_track(i)%orb(ix+1)%vec) / (s1 + s2)
!    endif

    f0 = (ric%ele%value(l$) - s) / ric%ele%value(l$)
    f1 = s / ric%ele%value(l$)

    orb%vec = ric%orb0%vec * f0 + ric%orb1%vec * f1
    call calc_g_params (orb)

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
      matmul(v, (/ ric%runt%x%eta, ric%runt%x%etap, 0.0_rdef,   0.0_rdef    /))
  ric%eta_b = &
      matmul(v, (/ 0.0_rdef,   0.0_rdef,    ric%runt%y%eta, ric%runt%y%etap /))

  ric%g_x = ric%g_x0 + orb%x%pos * ric%k1 + orb%y%pos * ric%s1
  ric%g_y = ric%g_y0 - orb%y%pos * ric%k1 + orb%x%pos * ric%s1
                   
  ric%dg2_x = 2 * (ric%g_x * ric%k1 + ric%g_y * ric%s1)
  ric%dg2_y = 2 * (ric%g_x * ric%s1 - ric%g_y * ric%k1) 

  ric%g2 = ric%g_x**2 + ric%g_y**2
  ric%g = sqrt(ric%g2)

!----------------------------------------------------------------------------
contains

subroutine calc_g_params (orb)

  type (coord_struct) orb
  real(rdef) del

  del = 1e-6
  here = orb%vec
  call derivs_bmad (ric%ele, ric%ring%param, s, here, kick_0)
  here(1) = here(1) + del
  call derivs_bmad (ric%ele, ric%ring%param, s, here, kick_x)
  here = orb%vec
  here(3) = here(3) + del
  call derivs_bmad (ric%ele, ric%ring%param, s, here, kick_y)
  ric%g_x = -kick_0(2)
  ric%g_y = -kick_0(4)
  ric%g2 = ric%g_x**2 + ric%g_y**2
  ric%g  = sqrt(ric%g2)
  ric%dg2_x = ((kick_x(2)**2 + kick_x(4)**2) - ric%g2) / del
  ric%dg2_y = ((kick_y(2)**2 + kick_y(4)**2) - ric%g2) / del

end subroutine

end subroutine

end module
