MODULE longitudinal_profile_mod

USE bmad
USE fgsl
USE, INTRINSIC :: iso_c_binding

IMPLICIT none

INTEGER(fgsl_size_t), PARAMETER :: limit = 1000_fgsl_size_t

CONTAINS

!+
!-
FUNCTION psi_prime(t, p, dpdt, params) BIND(c)
  IMPLICIT none

  REAL(c_double), VALUE :: t
  REAL(c_double), DIMENSION(*), INTENT(IN) :: p
  REAL(c_double), DIMENSION(*) :: dpdt
  TYPE(c_ptr), VALUE :: params
  INTEGER(c_int) :: psi_prime

  REAL(c_double), POINTER :: args(:)
  REAL(c_double) A, Vrf, Q, omega, phi, R, L, U0

  CALL c_f_pointer(params,args,[8])
  A = args(1)
  Vrf = args(2)
  Q = args(3)
  omega = args(4)
  phi = args(5)
  R = args(6)
  L = args(7)
  U0 = args(8)
  dpdt(1) = -A*p(1)*((Vrf*COS(omega*t+phi) + Q*R*p(1) - U0)/(1 + A*Q*L*p(1)))

  psi_prime = fgsl_success
END FUNCTION psi_prime

!+
!-
SUBROUTINE psi_prime_sca(t, p, dpdt, args)
  IMPLICIT NONE

  REAL(rp) t, p, dpdt
  REAL(rp) args(1:8)
  REAL(c_double) pa(1), dpdta(1)
  TYPE(c_ptr) ptr

  INTEGER status

  ptr = c_loc(args)

  pa(1) = p
  status = psi_prime(t, pa, dpdta, ptr)
  dpdt = dpdta(1)
END SUBROUTINE psi_prime_sca

!+
!-
FUNCTION jac(t, p, dfdp, dfdt, params) BIND(c)
  IMPLICIT none
  REAL(c_double), VALUE :: t
  REAL(c_double), DIMENSION(*), INTENT(IN) :: p
  REAL(c_double), DIMENSION(*) :: dfdp
  REAL(c_double), DIMENSION(*) :: dfdt
  TYPE(c_ptr), VALUE :: params
  INTEGER(c_int) :: jac

  REAL(c_double), POINTER :: args(:)
  REAL(c_double) A, Vrf, Q, omega, phi, R, L, U0

  CALL c_f_pointer(params,args,[8])
  A = args(1)
  Vrf = args(2)
  Q = args(3)
  omega = args(4)
  phi = args(5)
  R = args(6)
  L = args(7)
  U0 = args(8)

  dfdp(1) = -A*((Vrf*COS(omega*t+phi) - U0 + Q*R*p(1)*(2.0_rp + A*L*Q*p(1))) / ((1.0_rp + A*L*Q*p(1))**2))

  dfdt(1) = A*Vrf*p(1)*omega*SIN(omega*t + phi) / (1.0_rp + A*L*Q*p(1))

  jac = fgsl_success
END FUNCTION jac

!+
!-
SUBROUTINE solve_psi_adaptive(t0,t1,p0,args,p1)
  ! t0    initial time
  ! t1    final time
  ! p0    boundary condition psi(t0)
  ! args  arguments for DEQ
  ! p1    psi(t1)
  ! 
  ! Step size is adaptive.
  !
  IMPLICIT none
  
  REAL(fgsl_double), VALUE :: t0  ! fgsl_odeiv2_driver_* fails without VALUE
  REAL(fgsl_double) t1
  REAL(rp) p0
  REAL(rp) args(1:8)
  REAL(rp) p1
  type(fgsl_odeiv2_system) :: ode_system
  TYPE(fgsl_odeiv2_driver) :: ode_drv
  TYPE(c_ptr) ptr

  INTEGER status
  REAL(fgsl_double) y(1)
  REAL(fgsl_double) default_step
  REAL(fgsl_double), PARAMETER :: abs_err_goal = 0.0_fgsl_double
  REAL(fgsl_double), PARAMETER :: rel_err_goal = 1.0d-6

  default_step = (t1-t0)/100.

  ptr = c_loc(args)
  ode_system = fgsl_odeiv2_system_init(psi_prime, 1_c_size_t, ptr, jac)
  ode_drv = fgsl_odeiv2_driver_alloc_y_new(ode_system, fgsl_odeiv2_step_bsimp, default_step, abs_err_goal, rel_err_goal)

  y(1) = p0

  status = fgsl_odeiv2_driver_apply(ode_drv, t0, t1, y)

  IF(status .ne. fgsl_success) THEN
    WRITE(*,'(A)') "ERROR: fgsl_odeiv2_driver_apply failed during bunch length calculation."
    WRITE(*,'(A,2I)') "fgsl_odeiv2_driver_apply returned (success is zero): ", status
    STOP
  ENDIF

  p1 = y(1)
  
  CALL fgsl_odeiv2_system_free(ode_system)
  CALL fgsl_odeiv2_driver_free(ode_drv)
END SUBROUTINE solve_psi_adaptive

!+
!-
SUBROUTINE solve_psi_fixed_steps(t0,t1,p0,args,t,p)
  ! t0    initial time
  ! t1    final time
  ! p0    boundary condition psi(t0)
  ! args  arguments for DEQ
  ! t(:)  array of times
  ! p(:)  array of psi
  ! 
  ! Step size is fixed and determined by (t1-t0)/size(p)
  !
  IMPLICIT none
  
  REAL(rp) t0
  REAL(rp) t1
  REAL(rp) p0
  REAL(rp) args(1:8)
  REAL(rp) t(:)
  REAL(rp) p(:)
  type(fgsl_odeiv2_system) :: ode_system
  TYPE(fgsl_odeiv2_driver) :: ode_drv
  TYPE(c_ptr) ptr

  INTEGER i
  INTEGER n
  INTEGER status
  REAL(fgsl_double) y(1)
  REAL(fgsl_double) tcur
  REAL(fgsl_double) step_size

  n = SIZE(p)
  step_size = (t1-t0)/(n-1)

  ptr = c_loc(args)
  ode_system = fgsl_odeiv2_system_init(psi_prime, 1_c_size_t, ptr, jac)
  ode_drv = fgsl_odeiv2_driver_alloc_y_new(ode_system, fgsl_odeiv2_step_bsimp, 1.0e-6_fgsl_double, 0.0_fgsl_double, 1.0E-6_fgsl_double)

  tcur = t0
  y(1) = p0

  t(1) = t0
  p(1) = p0

  DO i=2,n
    status = fgsl_odeiv2_driver_apply_fixed_step(ode_drv, tcur, step_size, 1, y)

    IF(status .ne. fgsl_success) THEN
      WRITE(*,'(A)') "ERROR: fgsl_odeiv2_driver_apply_fixed_step failed during bunch length calculation."
      WRITE(*,'(A,2I)') "fgsl_odeiv2_driver_apply_fixed_step returned ", status, fgsl_success
      STOP
    ENDIF

    t(i) = tcur
    p(i) = y(1)
  ENDDO
  
  CALL fgsl_odeiv2_system_free(ode_system)
  CALL fgsl_odeiv2_driver_free(ode_drv)
END SUBROUTINE solve_psi_fixed_steps

!+
!-
SUBROUTINE integrate_psi(bound,p0,args,result)
  REAL(rp) bound
  REAL(rp) p0
  INTEGER npts
  REAL(rp) args(1:8)
  REAL(rp) result

  INTEGER i

  REAL(rp), ALLOCATABLE :: t(:)
  REAL(rp), ALLOCATABLE :: tminus(:)
  REAL(rp), ALLOCATABLE :: tplus(:)
  REAL(rp), ALLOCATABLE :: p(:)
  REAL(rp), ALLOCATABLE :: pminus(:)
  REAL(rp), ALLOCATABLE :: pplus(:)

  npts = 100

  ALLOCATE(tplus(1:npts))
  ALLOCATE(tminus(1:npts))
  ALLOCATE(t(1:2*npts-1))
  ALLOCATE(pplus(1:npts))
  ALLOCATE(pminus(1:npts))
  ALLOCATE(p(1:2*npts-1))

  CALL solve_psi_fixed_steps(0.0_rp,-bound,p0,args,tminus,pminus)
  CALL solve_psi_fixed_steps(0.0_rp,bound,p0,args,tplus,pplus)

  DO i=1,npts
    t(i) = tminus(npts-i+1)
    p(i) = pminus(npts-i+1)
  ENDDO
  t(npts+1:2*npts-1) = tplus(2:npts)
  p(npts+1:2*npts-1) = pplus(2:npts)

  result = 0.0d0
  DO i=1,2*npts-2
    result = result + (t(i+1)-t(i))*(p(i+1)+p(i))/2.0d0
  ENDDO

  result = result - 1.0d0

  DEALLOCATE(tplus)
  DEALLOCATE(tminus)
  DEALLOCATE(t)
  DEALLOCATE(pplus)
  DEALLOCATE(pminus)
  DEALLOCATE(p)

END SUBROUTINE integrate_psi

!+
!-
SUBROUTINE find_normalization(bound,p0,args,pnrml)
  !Secant Method
  REAL(rp) bound
  REAL(rp) p0
  REAL(rp) args(1:8)
  REAL(rp) pnrml

  REAL(rp) f(1:3) !f(1) is f(n), f(2) is f(n-1), f(3) is f(n-2)
  REAL(rp) x(1:3) !x(1) is x(n), x(2) is x(n-1), x(3) is x(n-2)

  x(3) = p0*0.95d0
  CALL integrate_psi(bound,x(3),args,f(3))

  x(2) = p0*1.05d0
  DO WHILE(.true.)
    CALL integrate_psi(bound,x(2),args,f(2))
    x(1) = x(2) - f(2)*(x(2)-x(3))/(f(2)-f(3))

    IF(ABS((x(1)-x(2))/x(2)) .lt. 1.0d-6) EXIT

    x(3) = x(2)
    x(2) = x(1)
    f(3) = f(2)
  ENDDO

  pnrml = x(1)
END SUBROUTINE find_normalization

!+
!-
SUBROUTINE find_fwhm(bound,args,fwhm)
  IMPLICIT NONE

  REAL(rp) bound
  REAL(rp) args(1:8)
  REAL(rp) fwhm
  REAL(rp) p0
  REAL(rp) pnrml
  REAL(rp) ta(1:2)
  REAL(rp) pa(1:2)
  REAL(rp) p2
  REAL(rp) half_max_psi
  REAL(rp) lower_value
  REAL(rp) upper_value
  REAL(rp) xmax, xmin, xnew
  REAL(rp) fnew, dpdt
  REAL(rp) f(1:3)
  REAL(rp) x(1:3)

  REAL(rp) max_time, max_psi

  INTEGER status

  !-
  !- Step 1: Find Max
  !-
  !- Secant method not guaranteed to converge ... using bisection instead
  p0 = 1.0d9 !initial guess for normalization
  CALL find_normalization(bound,p0,args,pnrml)

  xmax = 0.0d0
  xmin = -bound

  DO WHILE(.true.)
    xnew = xmin+(xmax-xmin)/2.0d0
    CALL solve_psi_adaptive(0.0d0,xnew,pnrml,args,fnew)
    CALL psi_prime_sca(xnew, fnew, dpdt, args)

    IF( dpdt .gt. 0.0d0 ) THEN
      xmin = xnew
    ELSE
      xmax = xnew
    ENDIF

    IF(ABS((xmax-xmin)/xmax) .lt. 1.0d-6) EXIT
  ENDDO

  max_time = xnew
  max_psi = fnew
  half_max_psi = max_psi / 2.0d0

  !-
  !- Step 2: Find Lower Value
  !-
  x(3) = -bound
  CALL solve_psi_adaptive(0.0d0,x(3),pnrml,args,f(3))
  f(3) = f(3) - half_max_psi

  x(2) = max_time
 
  DO WHILE(.true.)
    CALL solve_psi_adaptive(0.0d0,x(2),pnrml,args,f(2))
    f(2) = f(2) - half_max_psi

    x(1) = x(2) - f(2) * (x(2)-x(3)) / (f(2)-f(3))
    x(1) = MIN(x(1), max_time)

    IF( ABS((x(1)-x(2))/x(1)) .lt. 1.0d-6 ) EXIT

    x(3) = x(2)
    f(3) = f(2)
    x(2) = x(1)
  ENDDO

  lower_value = x(1)

  !-
  !- Step 3: Find Upper Value
  !-
  x(3) = bound
  CALL solve_psi_adaptive(0.0d0,x(3),pnrml,args,f(3))
  f(3) = f(3) - half_max_psi

  x(2) = max_time
 
  DO WHILE(.true.)
    CALL solve_psi_adaptive(0.0d0,x(2),pnrml,args,f(2))
    f(2) = f(2) - half_max_psi

    x(1) = x(2) - f(2) * (x(2)-x(3)) / (f(2)-f(3))
    x(1) = MAX(x(1), max_time)

    IF( ABS((x(1)-x(2))/x(1)) .lt. 1.0d-6 ) EXIT

    x(3) = x(2)
    f(3) = f(2)
    x(2) = x(1)
  ENDDO

  upper_value = x(1)

  fwhm = upper_value - lower_value

END SUBROUTINE find_fwhm

!+
!-
SUBROUTINE get_bl_from_fwhm(bound,args,sigma)
  IMPLICIT none

  REAL(rp) bound
  REAL(rp) args(1:8)
  REAL(rp) sigma
  REAL(rp) fwhm
  REAL(rp), PARAMETER :: TwoRtTwoLnTwo = 2.354820045

  CALL find_fwhm(bound,args,fwhm)

  sigma = fwhm * c_light / TwoRtTwoLnTwo
END SUBROUTINE get_bl_from_fwhm


END MODULE longitudinal_profile_mod










