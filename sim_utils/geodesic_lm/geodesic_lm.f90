module geodesic_lm

use sim_utils

type geodesic_lm_param_struct
  integer :: mode = 0           ! LM damping matrix. 0->id, 1->dynamic jacob-based
  integer :: maxiter = 500      ! max # of routine iterations
  integer :: maxfev = 0         ! max # of func evals (0-> no limit)
  integer :: maxjev = 0         ! max number of jac evals (0->no limit)
  integer :: maxaev = 0         ! max number of dir second derivs (0->no limit)
  integer :: print_level = 3    ! how many details to be printed (0-5) 
  integer :: print_unit = 6     ! unit number details written to
  integer :: imethod = 10       ! method choice for updating LM parameter
  integer :: iaccel = 1         ! use geodesic acceleration or not
  integer :: ibold = 0          ! 'boldness' in accepting uphill (0-4) with 0->downhill
  integer :: ibroyden = 0       ! number of iterations using approximate jacobian

  !!!!real(rp) :: eps = 1.5E-6      ! function evaluation precision
  real(rp) :: h1=1.D-6,h2=1.D-1 ! controls step sizes for finite diff derivatives
                                ! h1 for jacobian, h2 for dir second deriv
  !! Stopping criterion
  real(rp) :: maxlam = 1E7      ! maximum limit on damping term lambda (if <0 no limit)
  real(rp) :: minlam = -1.0     ! minimum limit on damping term lambda (if <0 no limit)
  real(rp) :: artol = 1.E-3     ! cos of angle between residual and tangent plane
  real(rp) :: Cgoal  = 1.E-23   ! Cost lower limit (ends when falls below)
  real(rp) :: gtol  = 1.5E-8    ! gradient lower limit
  real(rp) :: xtol = 1.E-10     ! step size lower limit (ll)
  real(rp) :: xrtol = 1.5E-8    ! relative parameter change ll
  real(rp) :: ftol = 1.5E-8     ! consecutive cost difference ll
  real(rp) :: frtol = 1.5E-8    ! relative consecutive cost diff ll
  
  real(rp) :: initialfactor = 1. ! initial LM param or step size
  real(rp) :: factoraccept  = 5. ! (if imethod=0 or 10) adjusts initialfactor
  real(rp) :: factorreject  = 2. ! adjusts initialfactor for rejected step
  real(rp) :: avmax = 0.8        ! limits geo accel w.r.t. velocity

  logical :: analytic_jac = .true.
  logical :: analytic_avv = .false.
  logical :: center_diff = .true.

  logical :: geo_hit_limit= .false.  !flag for when limits are hit
end type 

type (geodesic_lm_param_struct), save, target ::  geodesic_lm_param

contains

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Subroutine type_geodesic_lm (lines, n_lines)
!
! Routine to print or put information into a string array of the geodesic_lm parameters.
! If "lines" is not present, the information will be printed to the screen.
!
! Input:
!   print_coords -- logical, optional: If True then print coordinate and  patch information.
!                     Default is True.
!   prefix       -- character(*), optional: prefix for output lines.
!                                  Default: 'geodesic_lm_param%'
! Output:
!   lines(:)  -- character(120), optional, allocatable: Character array to hold the output.
!   n_lines   -- integer, optional: Number of lines used in lines(:)

subroutine type_geodesic_lm (lines, n_lines, prefix)

integer, optional :: n_lines
integer i, nl

character(*), allocatable, optional :: lines(:)
character(160) :: li(40)
character(20) imt, rmt, lmt
character(*), optional :: prefix
character(20) :: pre
integer :: nchar
if (present(prefix)) then
  pre = prefix
  nchar=len(prefix)
else
  pre = 'geodesic_lm_param%'
  nchar = len(trim(pre))
endif

!

rmt  = '(2a, 9es16.8)'
imt  = '(2a, 9i8)'
lmt  = '(2a, 9(l3))'

nl = 0
nl=nl+1; write (li(nl), imt) pre(1:nchar), 'mode             =', geodesic_lm_param%mode
nl=nl+1; write (li(nl), imt) pre(1:nchar), 'maxiter          =', geodesic_lm_param%maxiter
nl=nl+1; write (li(nl), imt) pre(1:nchar), 'maxfev           =', geodesic_lm_param%maxfev
nl=nl+1; write (li(nl), imt) pre(1:nchar), 'maxjev           =', geodesic_lm_param%maxjev
nl=nl+1; write (li(nl), imt) pre(1:nchar), 'maxaev           =', geodesic_lm_param%maxaev
nl=nl+1; write (li(nl), imt) pre(1:nchar), 'print_level      =', geodesic_lm_param%print_level
nl=nl+1; write (li(nl), imt) pre(1:nchar), 'print_unit       =', geodesic_lm_param%print_unit
nl=nl+1; write (li(nl), imt) pre(1:nchar), 'imethod          =', geodesic_lm_param%imethod
nl=nl+1; write (li(nl), imt) pre(1:nchar), 'iaccel           =', geodesic_lm_param%iaccel
nl=nl+1; write (li(nl), imt) pre(1:nchar), 'ibold            =', geodesic_lm_param%ibold
nl=nl+1; write (li(nl), imt) pre(1:nchar), 'ibroyden         =', geodesic_lm_param%ibroyden

!!!nl=nl+1; write (li(nl), rmt) pre(1:nchar), 'eps              =', geodesic_lm_param%eps
nl=nl+1; write (li(nl), rmt) pre(1:nchar), 'h1               =', geodesic_lm_param%h1
nl=nl+1; write (li(nl), rmt) pre(1:nchar), 'h2               =', geodesic_lm_param%h2
nl=nl+1; write (li(nl), rmt) pre(1:nchar), 'maxlam           =', geodesic_lm_param%maxlam
nl=nl+1; write (li(nl), rmt) pre(1:nchar), 'minlam           =', geodesic_lm_param%minlam
nl=nl+1; write (li(nl), rmt) pre(1:nchar), 'artol            =', geodesic_lm_param%artol
nl=nl+1; write (li(nl), rmt) pre(1:nchar), 'Cgoal            =', geodesic_lm_param%Cgoal
nl=nl+1; write (li(nl), rmt) pre(1:nchar), 'gtol             =', geodesic_lm_param%gtol
nl=nl+1; write (li(nl), rmt) pre(1:nchar), 'xtol             =', geodesic_lm_param%xtol
nl=nl+1; write (li(nl), rmt) pre(1:nchar), 'xrtol            =', geodesic_lm_param%xrtol
nl=nl+1; write (li(nl), rmt) pre(1:nchar), 'ftol             =', geodesic_lm_param%ftol
nl=nl+1; write (li(nl), rmt) pre(1:nchar), 'frtol            =', geodesic_lm_param%frtol
nl=nl+1; write (li(nl), rmt) pre(1:nchar), 'initialfactor    =', geodesic_lm_param%initialfactor
nl=nl+1; write (li(nl), rmt) pre(1:nchar), 'factoraccept     =', geodesic_lm_param%factoraccept
nl=nl+1; write (li(nl), rmt) pre(1:nchar), 'factorreject     =', geodesic_lm_param%factorreject
nl=nl+1; write (li(nl), rmt) pre(1:nchar), 'avmax            =', geodesic_lm_param%avmax

nl=nl+1; write (li(nl), lmt) pre(1:nchar), 'analytic_jac     =', geodesic_lm_param%analytic_jac
nl=nl+1; write (li(nl), lmt) pre(1:nchar), 'analytic_avv     =', geodesic_lm_param%analytic_avv
nl=nl+1; write (li(nl), lmt) pre(1:nchar), 'center_diff      =', geodesic_lm_param%center_diff

if (present(lines)) then
  call re_allocate(lines, nl, .false.)
  n_lines = nl
  lines(1:nl) = li(1:nl)
else
  do i = 1, nl
    print *, trim(li(i))
  enddo
endif

end subroutine type_geodesic_lm

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!Perhaps remove geo_func, jacobian, Avv, callback from signature as they don't
!need to be there
!+
! Subroutine run_geodesic_lm (geo_func, jacobian, Avv, a, y_fit, fjac, callback, info, &
!                             dtd,  niters, nfev, njev, naev, converged)
!-

subroutine run_geodesic_lm (geo_func, jacobian, Avv, a, y_fit, fjac, callback, info, &
                             dtd,  niters, nfev, njev, naev, converged)

use geolevmar_module, only: geodesicLM ! contained in leastsq.f90

implicit none

type (geodesic_lm_param_struct), pointer :: g

real(rp) a(:), y_fit(:)
real(rp) dtd(:,:), fjac(:,:)

integer info, niters, nfev, njev, naev, converged
integer m, n

!User provides these subroutines. Can be trivial if parameters in
!geodesic_lm_param are set to not use jacobian, Avv. geo_func and callback are
!necessary for all cases.

interface
  subroutine geo_func (mm,nn, a, y_fit)
    import
    implicit none
    ! a is the array of variables, y_fit is the model data
    integer :: mm, nn
    real(rp) :: a(nn)
    real(rp) :: y_fit(mm)
    !returns an array of residuals
  end subroutine

  subroutine jacobian(m, n, x, fjac)
    implicit none
    integer :: m, n
    real(8) :: x(n), fjac(m,n)
    !calculates the jacobian of the residuals w.r.t. parameters
  end subroutine

  subroutine Avv(m,n,x,v,acc)
    implicit none
    integer :: m,n
    real(8) :: x(n),v(n),acc(n)
    !Calculates second directional derivative of geo_func along direction v
    !returning acceleration vector acc
  end subroutine

  subroutine callback(m,n,x,fvec,fjac,accepted,info)     
    implicit none
    integer m, n, info, accepted
    real(8) x(n), fvec(m)
    real(8) fjac(m,n)
    !Handles cases when algorithm oversteps bounds of calling routine
  end subroutine
end interface

!

g => geodesic_lm_param

m = size(fjac, 1) !number of parameters
n = size(fjac, 2) !number of functions

call geodesicLM(geo_func, jacobian, Avv, a, y_fit, fjac, n, m, callback, info, &
            g%analytic_jac, g%analytic_Avv, g%center_diff, g%h1, g%h2, &
            dtd, g%mode, niters, nfev, njev, naev, &
            g%maxiter, g%maxfev, g%maxjev, g%maxaev, g%maxlam, g%minlam, &
            g%artol, g%Cgoal, g%gtol, g%xtol, g%xrtol, g%ftol, g%frtol,  &
            converged, g%print_level, g%print_unit, &
            g%imethod, g%iaccel, g%ibold, g%ibroyden, &
            g%initialfactor, g%factoraccept, g%factorreject, g%avmax)


end subroutine run_geodesic_lm

end module
