!+
! Subroutine derivs (x, y, dydx)
!
! Subroutine used by odeint. A Numerical Recipes routine.
! For more information see the odeint writeup in the Numerical Recipes book.
! This subroutine calls field_rk which must be supplied by the user.
!
! Note: odeint is called by track_runge_kutta.
!
! Input:
!   x       -- real(rdef): Dependent var. Typically: s
!   y(:)    -- real(rdef): Generally: (x, p_x, y, p_y, z, p_z)
!
! Output:
!   dydx(:) -- real(rdef): Derivatives.
!-

subroutine derivs (x, y, dydx)

  use bmad
  use nrtype

  implicit none
                                   
  real(sp), intent(in) :: x         ! s-position
  real(sp), intent(in) :: y(:)      ! (x, x', y, y', z, z')
  real(sp), intent(out) :: dydx(:)

  type (coord_struct) here

  real(rdef) field(3)
  real(rdef) vel_x, vel_y, vel_s, dvel_x, dvel_y, dvel_s, f

! calculate the field

  here%vec = y
  here%z%pos = x

  call field_rk (here, field)

! if this is a kick field then field gives us directly dydx

  if (bmad_common%func_type == 'KICK_FIELD') then
    dydx(1) = y(2)    ! dx/ds = 
    dydx(2) = field(1)
    dydx(3) = y(4)
    dydx(4) = field(2)
    dydx(5) = -(y(2)**2 + y(4)**2) / 2
    dydx(6) = field(3)
    return
  endif

! Here for func_type = B_FIELD
! The computation (up to some constant factors):
!     x' = dx/ds = v_x / v_s
!     dx'/ds = (dv_x/dt * v_s - v_s * dv_s/dt) / v_s^3  ! ds/dt == v_s
! where
!   dv_x/dt = v_y * B_s - v_s * B_y  
!   dv_y/dt = v_s * B_x - v_x * B_s  
!   dv_s/dt = v_x * B_y - v_y * B_x

  vel_x = y(2)                              ! proportional to x-velosity
  vel_y = y(4)                              ! proportional to y-velosity
  vel_s = 1/sqrt(1 + y(2)**2 + y(4)**2)     ! proportional to s-velosity

  dvel_x = vel_y * field(3) - vel_s * field(2)
  dvel_y = vel_s * field(1) - vel_x * field(3)
  dvel_s = vel_x * field(2) - vel_y * field(1)

  f = bmad_common%factor / (1 + y(6))

  dydx(1) = y(2)    
  dydx(2) = f * (dvel_x * vel_s - vel_x * dvel_s) / vel_s**3
  dydx(3) = y(4)
  dydx(4) = f * (dvel_y * vel_s - vel_y * dvel_s) / vel_s**3
  dydx(5) = -(y(2)**2 + y(4)**2) / 2
  dydx(6) = 0           ! dE/ds = 0

end subroutine
                                          
