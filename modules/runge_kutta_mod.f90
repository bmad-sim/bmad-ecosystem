#include "CESR_platform.inc"      

module runge_kutta_mod

  use bmad_struct

  type runge_kutta_com_struct
    logical :: save_steps = .true.         ! save orbit?
    integer :: n_ok, n_bad, n_pts          ! number of points
    real(rp) :: ds_save = 1e-3           ! min distance between points
    real(rp), pointer :: s(:) => null()  ! s-distance of a point
    type (coord_struct), pointer :: orb(:) ! position of a point
  end type

  type (runge_kutta_com_struct), save :: rk_com

  integer, parameter :: kick_field$ = 2   ! b_field$ is defined in bmad_struct

  interface 
    subroutine field_rk_custom (ele, param, s, orb, field, field_type)
      use bmad_struct
      implicit none
      type (ele_struct), intent(in) :: ele
      type (param_struct) param
      type (coord_struct), intent(in) :: orb
      real(rp), intent(in) :: s
      real(rp), intent(out) :: field(3)
      integer, intent(out) :: field_type
    end subroutine
  end interface

contains

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine odeint_bmad (start, ele, param, end, &
!                                 s1, s2, rel_tol, abs_tol, h1, hmin)
! 
! Subroutine to do Runge Kutta tracking. This routine is adapted from Numerical
! Recipes.  See the NR book for more details.
!
! Notice that this routine has an two tollerance arguments rel_tol and abs_tol.
! Odein only has 1. rel_tol (essentually equivalent to eps in odeint) 
! is scalled by the step size to to able to relate it to the final accuracy.
!
! Essentually (assuming random errors) one of these conditions holds:
!             %error in tracking < rel_tol
! or
!     absolute error in tracking < abs_tol
!
! Modules needed:
!   use bmad
!
! Input: 
!   start   -- Coord_struct: Starting coords.
!   ele     -- Ele_struct: Element to track through.
!     %tracking_method -- Determines which subroutine to use to calculate the 
!                         field. Note: BMAD does no supply field_rk_custom.
!                           == custom$ then use field_rk_custom
!                           /= custom$ then use field_rk_standard
!   param   -- Param_struct: Beam parameters.
!     %enegy       -- Energy in GeV
!     %particle    -- Particle type [positron$, or electron$]
!   s1      -- Real: Starting point.
!   s2      -- Real: Ending point.
!   rel_tol -- Real: Same as eps for odeint scalled by sqrt(h/(s2-s1))
!               where h is the step size for the current step. rel_tol
!               sets the %error of the result
!   abs_tol -- Real: Sets the absolute error of the result
!   h1      -- Real: Initial guess for a step size.
!   h_min   -- Real: Minimum step size (can be zero).
!
! Output:
!   end        -- Coord_struct: Ending coords.
!
! Common block:
!   rk_com -- common_block that holds the path.
!     %save_steps -- Set True if you want to save the path
!     %n_pts -- The number of data points
!     %s(:)   -- Real: S positions of the data points
!     %orb(:) -- Coord_struct: Coordinates.
!
!-

subroutine odeint_bmad (start, ele, param, end, &
                                  s1, s2, rel_tol, abs_tol, h1, h_min)

  implicit none

  type (coord_struct), intent(in) :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct) ele
  type (param_struct) param

  real(rp), intent(in) :: s1, s2, rel_tol, abs_tol, h1, h_min

  real(rp), parameter :: tiny = 1.0e-30_rp
  real(rp) :: h, h_did, h_next, s, s_sav, rel_tol_eff, abs_tol_eff, sqrt_N
  real(rp) :: dr_ds(6), r(6), r_scal(6)

  integer, parameter :: max_step = 10000
  integer :: n_step, n_pts_save_max

! init

  s = s1
  h = sign(h1, s2-s1)
  rk_com%n_ok = 0
  rk_com%n_bad = 0
  rk_com%n_pts = 0
  r = start%vec

! if we are saving the trajectory then allocate enough space in the arrays

  if (rk_com%save_steps) then
    s_sav = s - 2.0_rp * rk_com%ds_save
    n_pts_save_max = 2 + abs(s2-s1)/rk_com%ds_save
    if (associated(rk_com%s)) then
      if (size(rk_com%s) < n_pts_save_max) then
        deallocate(rk_com%s, rk_com%orb)
        allocate(rk_com%s(n_pts_save_max), rk_com%orb(n_pts_save_max))
      endif
    else
      allocate(rk_com%s(n_pts_save_max), rk_com%orb(n_pts_save_max))
    endif
  endif

! now track

  bmad_status%ok = .true.

  do n_step = 1, max_step

    call derivs_bmad (ele, param, s, r, dr_ds)

    sqrt_N = sqrt(abs((s2-s1)/h))  ! number of steps we would take with this h
    rel_tol_eff = rel_tol / sqrt_N
    abs_tol_eff = abs_tol / sqrt_N
    r_scal(:) = abs(r(:)) + abs(h*dr_ds(:)) + TINY

    if (rk_com%save_steps .and. (abs(s-s_sav) > rk_com%ds_save)) &
                                                       call save_a_step
    if ((s+h-s2)*(s+h-s1) > 0.0) h = s2-s

    call rkqs_bmad (ele, param, r, dr_ds, s, h, rel_tol_eff, abs_tol_eff, &
                                                       r_scal, h_did, h_next)
    if (.not. bmad_status%ok) return

    if (h_did == h) then
      rk_com%n_ok = rk_com%n_ok + 1
    else
      rk_com%n_bad = rk_com%n_bad+1
    end if

    if ((s-s2)*(s2-s1) >= 0.0) then
      end%vec = r
      if (rk_com%save_steps) call save_a_step
      return
    end if

    if (abs(h_next) < h_min) then
      bmad_status%ok = .false.
      if (bmad_status%type_out) print *, &
                'ERROR IN ODEINT_BMAD: STEPSIZE SMALLER THAN MINIMUM.' 
      if (bmad_status%exit_on_error) call err_exit
    endif

    h = h_next

  end do

  bmad_status%ok = .false.
  if (bmad_status%type_out) print *, 'ERROR IN ODEINT_BMAD: TOO MANY STEPS'
  if (bmad_status%exit_on_error) call err_exit

!-----------------------------------------------------------------

contains

subroutine save_a_step
  integer n_pts
  rk_com%n_pts = rk_com%n_pts + 1
  n_pts = rk_com%n_pts 
  if (n_pts > size(rk_com%s)) then
    print *, 'ERROR IN ODEINT_BMAD SAVE_A_STEP: ARRAY OVERFLOW!'
    call err_exit
  end if
  rk_com%s(n_pts) = s
  rk_com%orb(n_pts)%vec = r
  s_sav = s
end subroutine save_a_step

end subroutine odeint_bmad

!-----------------------------------------------------------------
!-----------------------------------------------------------------

subroutine rkqs_bmad (ele, param, r, dr_ds, s, h_try, rel_tol, abs_tol, &
                                             r_scal, h_did, h_next)

  implicit none

  type (ele_struct) ele
  type (param_struct) param

  real(rp), intent(inout) :: r(6)
  real(rp), intent(in)    :: dr_ds(6), r_scal(6)
  real(rp), intent(inout) :: s
  real(rp), intent(in)    :: h_try, rel_tol, abs_tol
  real(rp), intent(out)   :: h_did, h_next

  real(rp) :: err_max, h, h_temp, s_new
  real(rp) :: r_err(6), r_temp(6)
  real(rp), parameter :: safety = 0.9_rp, p_grow = -0.2_rp
  real(rp), parameter :: p_shrink = -0.25_rp, err_con = 1.89e-4

!

  h = h_try

  do

    call rkck_bmad (ele, param, r, dr_ds, s, h, r_temp, r_err)
    err_max = maxval(abs(r_err(:)/(r_scal(:)*rel_tol + abs_tol)))
    if (err_max <=  1.0) exit
    h_temp = safety * h * (err_max**p_shrink)
    h = sign(max(abs(h_temp), 0.1_rp*abs(h)), h)
    s_new = s + h

    if (s_new == s) then
      bmad_status%ok = .false.
      if (bmad_status%type_out) print *, &
                          'ERROR IN RKQS_BMAD: STEPSIZE UNDERFLOW'
      if (bmad_status%exit_on_error) call err_exit
      return
    endif

  end do

  if (err_max > err_con) then
    h_next = safety*h*(err_max**p_grow)
  else
    h_next = 5.0_rp * h
  end if

  h_did = h
  s = s+h
  r(:) = r_temp(:)

end subroutine rkqs_bmad

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine rkck_bmad (ele, param, r, dr_ds, s, h, r_out, r_err)

	implicit none

  type (ele_struct) ele
  type (param_struct) param

	real(rp), intent(in) :: r(6), dr_ds(6)
	real(rp), intent(in) :: s, h
	real(rp), intent(out) :: r_out(6), r_err(6)
	real(rp) :: ak2(6), ak3(6), ak4(6), ak5(6), ak6(6), r_temp(6)
	real(rp), parameter :: a2=0.2_rp, a3=0.3_rp, a4=0.6_rp, &
    a5=1.0_rp, a6=0.875_rp, b21=0.2_rp, b31=3.0_rp/40.0_rp, &
    b32=9.0_rp/40.0_rp, b41=0.3_rp, b42=-0.9_rp, b43=1.2_rp, &
    b51=-11.0_rp/54.0_rp, b52=2.5_rp, b53=-70.0_rp/27.0_rp, &
    b54=35.0_rp/27.0_rp, &
		b61=1631.0_rp/55296.0_rp, b62=175.0_rp/512.0_rp, &
		b63=575.0_rp/13824.0_rp, b64=44275.0_rp/110592.0_rp, &
		b65=253.0_rp/4096.0_rp, c1=37.0_rp/378.0_rp, &
		c3=250.0_rp/621.0_rp, c4=125.0_rp/594.0_rp, &
		c6=512.0_rp/1771.0_rp, dc1=c1-2825.0_rp/27648.0_rp, &
		dc3=c3-18575.0_rp/48384.0_rp, dc4=c4-13525.0_rp/55296.0_rp, &
		dc5=-277.0_rp/14336.0_rp, dc6=c6-0.25_rp

!

	r_temp=r+b21*h*dr_ds
	call derivs_bmad(ele, param, s+a2*h, r_temp, ak2)
	r_temp=r+h*(b31*dr_ds+b32*ak2)
	call derivs_bmad(ele, param, s+a3*h, r_temp, ak3)
	r_temp=r+h*(b41*dr_ds+b42*ak2+b43*ak3)
	call derivs_bmad(ele, param, s+a4*h, r_temp, ak4)
	r_temp=r+h*(b51*dr_ds+b52*ak2+b53*ak3+b54*ak4)
	call derivs_bmad(ele, param, s+a5*h, r_temp, ak5)
	r_temp=r+h*(b61*dr_ds+b62*ak2+b63*ak3+b64*ak4+b65*ak5)
	call derivs_bmad(ele, param, s+a6*h, r_temp, ak6)
	r_out=r+h*(c1*dr_ds+c3*ak3+c4*ak4+c6*ak6)
	r_err=h*(dc1*dr_ds+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6)

end subroutine rkck_bmad


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine derivs_bmad (ele, param, s, r, dr_ds)

  implicit none

  type (ele_struct) ele
  type (param_struct) param

  real(rp), intent(in) :: s         ! s-position
  real(rp), intent(in) :: r(6)      ! (x, x', y, y', z, z')
  real(rp), intent(out) :: dr_ds(6)

  type (coord_struct) here

  real(rp) field(3)
  real(rp) vel_x, vel_y, vel_s, dvel_x, dvel_y, dvel_s, f

  integer field_type

! calculate the field

  here%vec = r

  if (ele%tracking_method == custom$) then
    call field_rk_custom (ele, param, s, here, field, field_type)
  else 
    call field_rk_standard (ele, param, s, here, field, field_type)
  endif

! if this is a kick field then field gives us directly dr_ds

  if (field_type == KICK_FIELD$) then
    dr_ds(1) = r(2)    ! dx/ds =
    dr_ds(2) = field(1)
    dr_ds(3) = r(4)
    dr_ds(4) = field(2)
    dr_ds(5) = -(r(2)**2 + r(4)**2) / 2
    dr_ds(6) = field(3)
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

  if (field_type == B_FIELD$) then
    vel_x = r(2)                              ! proportional to x-velosity
    vel_y = r(4)                              ! proportional to y-velosity
    vel_s = 1/sqrt(1 + r(2)**2 + r(4)**2)     ! proportional to s-velosity

    dvel_x = vel_y * field(3) - vel_s * field(2)
    dvel_y = vel_s * field(1) - vel_x * field(3)
    dvel_s = vel_x * field(2) - vel_y * field(1)

    f = c_light / (param%beam_energy * (1 + r(6)))

    dr_ds(1) = r(2)
    dr_ds(2) = f * (dvel_x * vel_s - vel_x * dvel_s) / vel_s**3
    dr_ds(3) = r(4)
    dr_ds(4) = f * (dvel_y * vel_s - vel_y * dvel_s) / vel_s**3
    dr_ds(5) = -(r(2)**2 + r(4)**2) / 2
    dr_ds(6) = 0           ! dE/ds = 0
    return
  endif
  
  print *, 'ERROR IN DERIVS: UNKNOWN FIELD_TYPE: ', field_type
  call err_exit

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine field_rk (ele, param, s, here, field, field_type)

  implicit none

  type (ele_struct), intent(in) :: ele
  type (param_struct) param
  type (coord_struct), intent(in) :: here

  real(rp) :: s, field(3)

  integer, intent(out) :: field_type

!

  if (ele%tracking_method == custom$) then
    call field_rk_custom (ele, param, s, here, field, field_type)
  else
    call field_rk_standard (ele, param, s, here, field, field_type)
  endif

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine field_rk_standard (ele, param, s_pos, here, field, field_type)

  implicit none

  type (ele_struct), target, intent(in) :: ele
  type (param_struct) param
  type (coord_struct), intent(in) :: here
  type (wig_term_struct), pointer :: t

  real(rp) :: field(3), x, y, s, s_pos
  real(rp) :: c_x, s_x, c_y, s_y, c_z, s_z, coef

  integer, intent(out) :: field_type
  integer i

!

  x = here%vec(1)
  y = here%vec(3)
  s = s_pos

  field = 0

! Wiggler

  select case (ele%key)

  case(wiggler$)

    field_type = B_FIELD$

    do i = 1, size(ele%wig_term)
      t => ele%wig_term(i)

      if (t%type == hyper_y$) then
        c_x = cos(t%kx * x)
        s_x = sin(t%kx * x)
      else
        c_x =  cosh(t%kx * x)
        s_x = -sinh(t%kx * x)
      endif

      if (t%type == hyper_y$ .or. t%type == hyper_xy$) then
        c_y = cosh (t%ky * y)
        s_y = sinh (t%ky * y)
      else
        c_y = cos (t%ky * y)
        s_y = sin (t%ky * y)
      endif

      c_z = cos (t%kz * s + t%phi_z)
      s_z = sin (t%kz * s + t%phi_z)

      coef = t%coef * ele%value(polarity$)

      field(1) = field(1) - coef  * (t%kx / t%ky) * s_x * s_y * c_z
      field(2) = field(2) + coef  *                 c_x * c_y * c_z
      field(3) = field(3) - coef  * (t%kz / t%ky) * c_x * s_y * s_z
    enddo

! Quad

  case (quadrupole$)

    field_type = KICK_FIELD$
    field(1) = -ele%value(k1$) * x / (1 + here%vec(6))
    field(2) =  ele%value(k1$) * y / (1 + here%vec(6))

! Sextupole

  case (sextupole$)

    field_type = KICK_FIELD$
    field(1) = 0.5 * ele%value(k2$) * (y*y - x*x) / (1 + here%vec(6))
    field(2) =       ele%value(k2$) * x * y / (1 + here%vec(6))

! Error

  case default
    print *, 'ERROR IN FIELD_RK_STANDARD: ELEMENT NOT CODED: ', &
                                                         key_name(ele%key)
    print *, '      FOR: ', ele%name
    call err_exit
  end select

end subroutine

end module
