module integration_timer_mod

  use ptc_interface_mod
  use mad_like, only: kill, alloc, fibre, real_8, universal_taylor

  private get_taylor, diff, term_diff  

!+ 
! Subroutine integration_timer
!
! Subroutine to set the number of integration steps for a fibre or element
! so that the error of the Taylor map is within tol of the exact map.
!
! This subroutine overloades:
!    integration_timer_ele (ele, param, start, orb_max, tol)
!    integration_timer_fibre (a_fibre, orbit, orbit_max, tol_dp)
! 
! Input:
!   ele      -- Ele_struct: Element to track through.
!   param    -- lat_param_struct:
!   start    -- Coord_struct: Starting orbit.
!   orb_max  -- Coord_struct: maximum orbit.
!   tol      -- Real(rp): Tolerance for the map coefs. 
!               (Good number is, say, 1d-10).
!
!   a_fibre      -- Fibre: Fibre to track through.
!   orbit(6)     -- Real(dp): Starting orbit.   
!   orbit_max(6) -- Real(dp): Maximum orbit
!   tol_dp       -- Real(8): Tolerance for the map coefs. 
!                   (Good number is, say, 1d-10).
!
! Output:
!-

interface integration_timer
  module procedure integration_timer_fibre
  module procedure integration_timer_ele
end interface

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine integration_timer_ele (ele, param, start, orb_max, tol)

implicit none

type (ele_struct) ele
type (fibre), pointer :: a_fibre
type (coord_struct), intent(in) :: start, orb_max
type (lat_param_struct) param

real(rp) tol
real(dp) tol_dp, orbit(6), orbit_max(6)
logical err_flag

!

call ele_to_fibre (ele, a_fibre, .true., err_flag)

orbit = start%vec
orbit_max = orb_max%vec
orbit_max = abs(orbit_max)
tol_dp = tol

call integration_timer_fibre (a_fibre, orbit, orbit_max, tol_dp)

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine integration_timer_fibre (a_fibre, orbit, orbit_max, tol_dp)

use s_status, only: assignment(=)

implicit none

type (fibre) a_fibre
type (real_8) y(6)
type (universal_taylor) :: ut1(6), ut2(6)

real(dp), intent(in) :: tol_dp, orbit(6), orbit_max(6)
real t0, t_del, secnds

integer i, n_steps0, order, order0
integer n1_step, n2_step

! Init

n_steps0 = a_fibre%mag%p%nst
order0 = a_fibre%mag%p%method

do i = 1, 6
  ut1(i) = 0
  ut2(i) = 0
enddo

print *
print '(a, 1p, e10.2)', ' Tolerance:', tol_dp
print '(a, 1p, 6d10.2)', ' Orbit    :', (orbit(i), i = 1, 6)
print '(a, 1p, 6d10.2)', ' Orbit_max:', (orbit_max(i), i = 1, 6)
print *, 'Order Steps      Diff      Time'

! loop over all orders

do order = 2, 6, 2

  a_fibre%mag%p%method = order    
  a_fibre%magp%p%method = order    

  n1_step = 4
  call get_taylor (n1_step, ut1, a_fibre, y, orbit)

! double the step number until the taylor series converges

  do
    n2_step = 2 * n1_step
    t0 = secnds(0.0)
    call get_taylor (n2_step, ut2, a_fibre, y, orbit)
    t_del = secnds(t0)
    print '(2i6, 1pe10.2, 0pf10.1)', order, &
        n2_step, diff(ut1, ut2, orbit_max), t_del
    if (diff(ut1, ut2, orbit_max) < tol_dp) exit
    n1_step = n2_step
    ut1 = ut2
    if (n1_step > 1000) then
      print *, 'ERROR IN INTEGRATION_TIMER: ', &
                              'NOT CONVGING FOR INTEGRATION ORDER:', order
      if (global_com%exit_on_error) call err_exit
    endif
  enddo

end do

a_fibre%mag%p%nst  = n_steps0
a_fibre%magp%p%nst = n_steps0

a_fibre%mag%p%method  = order0
a_fibre%magp%p%method = order0

do i = 1, 6
  ut1(i) = -1
  ut2(i) = -1
enddo

end subroutine

!---------------------------------------------------------------------------

subroutine get_taylor (ns, ut, a_fibre, y, orbit)

use s_status, only: assignment(=)
use mad_like, only: ptc_track => track

implicit none

type (universal_taylor) ut(6)
type (fibre) a_fibre
type (real_8) y(6)

real(dp), intent(in) :: orbit(6)

integer ns, iz

!

a_fibre%mag%p%nst  = ns
a_fibre%magp%p%nst = ns

call alloc (y)

y = orbit
call ptc_track (a_fibre, y, ptc_private%base_state)  ! "track" in PTC

do iz = 1, 6
  if (associated(ut(iz)%c)) ut(iz) = -1
  ut(iz) = y(iz)%t
enddo

call kill (y)

end subroutine

!---------------------------------------------------------------------------

function diff (ut1, ut2, max_orbit) result (dif)

use s_status, only: assignment(=)

implicit none

type (universal_taylor) ut1(6), ut2(6), ut11, ut22

real(dp) dif
real(dp) max_orbit(:)

integer i, n1, n2, ord1, ord2, expn(6)

!

dif = 0  
ut11 = 0  ! nullify
ut22 = 0  ! nullify

do i = 1, 6

  call sort_universal_terms(ut1(i), ut11)
  call sort_universal_terms(ut2(i), ut22)

  n1 = 1
  n2 = 1

  do

    if (n1 <= ut11%n) then
      expn = ut11%j(n1,:) 
      ord1 = sum(expn)*10**6 + expn(6)*10**5 + expn(5)*10**4 + &
               expn(4)*10**3 + expn(3)*10**2 + expn(2)*10**1 + expn(1)
    else
      ord1 = 10**8  ! something large
    endif

    if (n2 <= ut22%n) then
      expn = ut22%j(n2,:) 
      ord2 = sum(expn)*10**6 + expn(6)*10**5 + expn(5)*10**4 + &
               expn(4)*10**3 + expn(3)*10**2 + expn(2)*10**1 + expn(1)
    else
      ord2 = 10**8  ! something large
    endif

    if (ord1 < ord2) then
      call term_diff(dif, ut11%c(n1), ut11%j(n1,:), max_orbit)
      n1 = n1 + 1
    elseif (ord1 > ord2) then
      call term_diff(dif, ut22%c(n2), ut22%j(n2,:), max_orbit)
      n2 = n2 + 1
    else
      call term_diff(dif, ut11%c(n1)-ut22%c(n2), ut11%j(n1,:), max_orbit)
      n1 = n1 + 1
      n2 = n2 + 1

    endif

    if (n1 > ut11%n .and. n2 > ut22%n) exit

  enddo

enddo

ut11 = -1; ut22 = -1  ! deallocate

end function

!---------------------------------------------------------------------------

subroutine term_diff (dif, coef, expn, max_orbit)

implicit none

real(dp) dif, coef, max_orbit(:), del
integer expn(6), k

!

del = coef
do k = 1, 6
  if (expn(k) /= 0) del = del * max_orbit(k)**expn(k)
enddo
dif = max(dif, abs(del))

end subroutine

end module
