!+
!-

!$Id$
!$Log$
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

  use bmad
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

  type (rad_int_common_struct) ric

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! This is a modified version of QROMB from Num. Rec.

contains

function qromb_rad_int (func, a, b, sum) result (rad_int)

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

  integer(i4b), parameter :: jmax = 6, jmaxp = jmax+1, k = 5, km = k-1
  integer(i4b) :: j

  real(rdef), intent(in) :: a, b, sum
  real(rdef) :: rad_int
  real(rdef) :: dqromb
  real(rdef), parameter :: eps = 1.0e-4_rdef, eps2 = 1.0e-6_rdef
  real(rdef), dimension(jmaxp) :: h, s

  logical :: debug = .false.

!

  h(1) = 1.0

  do j = 1, jmax
    call trapzd(func, a, b, s(j), j)
    if (j >=  k) then
      call polint(h(j-km:j), s(j-km:j), 0.0_rdef, rad_int, dqromb)
      if (abs(dqromb) <= eps * abs(rad_int) + eps2 * abs(sum)) return
    elseif (j >= 3) then
      call polint(h(1:j), s(1:j), 0.0_rdef, rad_int, dqromb)
      if (abs(dqromb) <= eps * abs(rad_int) + eps2 * abs(sum)) return
    end if
    s(j+1) = s(j)
    h(j+1) = 0.25_rdef * h(j)
  end do

  if (debug) then
    type *, 'Warning in RADIATION_INTEGRALS: Integral does not converge.'
  endif

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

  call propagate_part_way (s_vec)

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

end module
