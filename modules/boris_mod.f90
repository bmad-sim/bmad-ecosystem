!+
! Module boris_mod
!
! Module to do Boris integration tracking.
!
! Reference: 
!  "Efficiency of a Boris-like Integration Scheme with Spatial Stepping", 
!   P. H. Stolz et al., Physical Review Special Topics.
!   5, 094001 (2002).
!
! The variables in Boris tracking are:
!     (x, p_x, y, p_y, s-c*t, dE/E)
! At high energy s-c*t = z which is the distance of the particle from the 
! reference particle.
!-


#include "CESR_platform.inc"      

module boris_mod

  use bmad_struct
  use bmad_interface

  type track_com_struct
    logical :: save_steps = .false.        ! save orbit?
    integer :: n_pts                       ! number of points
    real(rdef), pointer :: s(:) => null()  ! s-distance of a point
    type (coord_struct), pointer :: orb(:) ! position of a point
    real(rdef) :: ds_save = 1e-3           ! min distance between points
    real(rdef) :: step0 = 1e-3             ! Initial step size.
    real(rdef) :: step_min = 1e-8          ! min step size to step below which 
                                           !   track1_adaptive_boris gives up
    integer :: max_step = 10000            ! maximum number of steps allowed
  end type

  type (track_com_struct) track_com

  interface 
    subroutine em_field_custom (ele, param, s, orb, field)
      use bmad_struct
      implicit none
      type (ele_struct), intent(in) :: ele
      type (param_struct) param
      type (coord_struct), intent(in) :: orb
      real(rdef), intent(in) :: s
      type (em_field_struct), intent(out) :: field
    end subroutine
  end interface

contains

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine track1_adaptive_boris (start, ele, param, end, s_start, s_end)
! 
! Subroutine to do Boris tracking with adaptive step size control.
! This routine is adapted from odeint in Numerical Recipes. 
! See the NR book for more details.
!
! Note:
!   For each step the error in the orbit must be:
!     error < (|orbit|*rel_tol + abs_tol) / sqrt(N)
! Where N is the number of steps that would be needed at the 
! present step size.
!
! Modules needed:
!   use bmad
!
! Input: 
!   start    -- Coord_struct: Starting coords.
!   ele      -- Ele_struct: Element to track through.
!     %tracking_method -- Determines which subroutine to use to calculate the 
!                         field. Note: BMAD does no supply em_field_custom.
!                           == custom$ then use em_field_custom
!                           /= custom$ then use em_field_standard
!     %value(rel_tol$) -- Real: Relative error tollerance.
!                           Default if zero: 1e-6.
!     %value(abs_tol$) -- Real: Absolute error tollerance.
!                           Default if zero: 1e-7.
!   param    -- Param_struct: Beam parameters.
!     %enegy       -- Energy in GeV
!     %particle    -- Particle type [positron$, or electron$]
!   s_start  -- Real, optional: Starting point.
!   s_end    -- Real, optional: Ending point.
!
! Output:
!   end        -- Coord_struct: Ending coords.
!
! Common block:
!   track_com   -- Common_block that holds the path.
!     %save_steps -- Set True if you want to save the path
!     %ds_save    -- min distance between points.
!     %n_pts      -- The number of data points
!     %s(:)       -- Real: S positions of the data points
!     %orb(:)     -- Coord_struct: Coordinates.
!     %step_min   -- Real: min step size to step below which 
!                      track1_adaptive_boris gives up.
!     %step0      -- Real: Initial guess for a step size.
!-

subroutine track1_adaptive_boris (start, ele, param, end, s_start, s_end)

  implicit none

  type (coord_struct), intent(in) :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct) ele
  type (param_struct) param

  type (coord_struct) here, orb1, orb2

  real(rdef), optional, intent(in) :: s_start, s_end

  real(rdef) :: ds, ds_did, ds_next, s, s_sav, rel_tol, abs_tol, sqrt_N
  real(rdef), parameter :: err_5 = 0.0324, safety = 0.9
  real(rdef) :: s1, s2, scale_orb, err_max, ds_temp, rel_tol_N, abs_tol_N

  integer :: n_step

! init

  if (present(s_start)) then
    s1 = s_start
  else
    s1 = 0
  endif

  if (present(s_end)) then
    s2 = s_end
  else 
    s2 = ele%value(l$)
  endif

  if (ele%value(rel_tol$) == 0) then
    rel_tol = 1e-6 
  else
    rel_tol = ele%value(rel_tol$) 
  endif

  if (ele%value(abs_tol$) == 0) then
    abs_tol = 1e-7 
  else
    abs_tol = ele%value(abs_tol$) 
  endif

  s = s1
  ds = sign(track_com%step0, s2-s1)

  here = start
  call offset_particle (ele, param, here, set$, set_canonical = .false.)
  call track_solenoid_edge (ele, param, set$, here)

! if we are saving the trajectory then allocate enough space in the arrays

  if (track_com%save_steps) then
    s_sav = s - 2 * track_com%ds_save
    call allocate_saved_orbit (int(abs((s2-s1)/track_com%ds_save)))
  endif

! now track

  bmad_status%ok = .true.

  do n_step = 1, track_com%max_step

    sqrt_N = sqrt(abs((s2-s1)/ds))  ! N = estimated number of steps
    rel_tol_N = rel_tol / sqrt_N
    abs_tol_N = abs_tol / sqrt_N

! record a track if we went far enough.

    if (track_com%save_steps .and. (abs(s-s_sav) > track_com%ds_save)) &
                               call save_a_step (ele, param, s, here, s_sav)

    if ((s+ds-s2)*(s+ds-s1) > 0.0) ds = s2-s

! Make A step. Keep shrinking the step until the error is within bounds.
! The error in a step is estimated by the difference in making one whole step
! or two half steps.

    do

      call track1_boris_partial (here, ele, param, s, ds/2, orb2) 
      call track1_boris_partial (orb2, ele, param, s+ds/2, ds/2, orb2)
      call track1_boris_partial (here, ele, param, s, ds, orb1) 
      scale_orb = maxval((abs(orb1%vec) + abs(orb2%vec))) / 2

      err_max = maxval(abs(orb2%vec - orb1%vec) / &
                                       (scale_orb*rel_tol_N + abs_tol_N))
      if (err_max <= 1) exit

      ds_temp = safety * ds / sqrt(err_max)
      ds = sign(max(abs(ds_temp), 0.1_rdef*abs(ds)), ds)

      if (abs(ds) < track_com%step_min) then
        bmad_status%ok = .false.
        if (bmad_status%type_out) print *, &
            'ERROR IN TRACK1_ADAPTIVE_BORIS: STEPSIZE SMALLER THAN MINIMUM.' 
        if (bmad_status%exit_on_error) call err_exit
        return
      endif

    enddo

! now that we have a good step record the present position and 
! calculate a new step size.

    here = orb2
    s = s + ds

    if (err_max > err_5) then   ! limit increase to no more than a factor of 5
      ds = safety * ds / sqrt(err_max)
    else
      ds = 5 * ds   
    endif

! check if we are done

    if ((s-s2)*(s2-s1) >= 0.0) then
      if (track_com%save_steps) call save_a_step (ele, param, s, here, s_sav)
      call track_solenoid_edge (ele, param, unset$, here)
      call offset_particle (ele, param, here, unset$, set_canonical = .false.)
      end = here
      return
    end if

  end do

  bmad_status%ok = .false.
  if (bmad_status%type_out) &
               print *, 'ERROR IN TRACK1_ADAPTIVE_BORIS: TOO MANY STEPS'
  if (bmad_status%exit_on_error) call err_exit

end subroutine track1_adaptive_boris

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine track1_boris (start, ele, param, end, s_start, s_end)
!
! Subroutine to do Boris tracking.
! 
! Modules needed:
!   use bmad
!
! Input: 
!   start    -- Coord_struct: Starting coords.
!   ele      -- Ele_struct: Element to track through.
!     %tracking_method -- Determines which subroutine to use to calculate the 
!                         field. Note: BMAD does no supply em_field_custom.
!                           == custom$ then use em_field_custom
!                           /= custom$ then use em_field_standard
!     %num_steps      -- number of steps to take
!   param    -- Param_struct: Beam parameters.
!     %enegy       -- Energy in GeV
!     %particle    -- Particle type [positron$, or electron$]
!   s_start  -- Real, optional: Starting point.
!   s_end    -- Real, optional: Ending point.
!
! Output:
!   end        -- Coord_struct: Ending coords.
!
! Common block:
!   track_com   -- Common_block that holds the path.
!     %save_steps -- Set True if you want to save the path
!     %n_pts      -- The number of data points
!     %s(:)       -- Real: S positions of the data points
!     %orb(:)     -- Coord_struct: Coordinates.
!-

subroutine track1_boris (start, ele, param, end, s_start, s_end)

  implicit none

  type (coord_struct), intent(in) :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct) ele
  type (param_struct) param

  type (coord_struct) here

  real(rdef), optional, intent(in) :: s_start, s_end
  real(rdef) s1, s2, s_sav, ds, s

  integer i

! init

  if (present(s_start)) then
    s1 = s_start
  else
    s1 = 0
  endif

  if (present(s_end)) then
    s2 = s_end
  else 
    s2 = ele%value(l$)
  endif

  if (ele%num_steps == 0) then
    print *, 'ERROR IN TRACK1_BORIS: ELEMENT HAS ZERO "NUM_STEPS": ', ele%name
    call err_exit
  endif

  s = s1
  ds = (s2-s1) / ele%num_steps

! go to local coords

  here = start
  call offset_particle (ele, param, here, set$, set_canonical = .false.)
  call track_solenoid_edge (ele, param, set$, here)

! if we are saving the trajectory then allocate enough space in the arrays

  if (track_com%save_steps) then
    s_sav = s - 2.0_rdef * track_com%ds_save
    call allocate_saved_orbit (ele%num_steps)
    call save_a_step (ele, param, s, here, s_sav)
  endif

! track through the body

  do i = 1, ele%num_steps
    call track1_boris_partial (here, ele, param, s, ds, here)
    s = s + ds
    if (track_com%save_steps) call save_a_step (ele, param, s, here, s_sav)
  enddo

! back to lab coords

  call track_solenoid_edge (ele, param, unset$, here)
  call offset_particle (ele, param, here, unset$, set_canonical = .false.)

  end = here

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine track1_boris_partial (start, ele, param, s, ds, end)

  implicit none

  type (ele_struct), intent(in) :: ele
  type (param_struct) param
  type (coord_struct) :: start, end  
  type (em_field_struct) :: field

  real(rdef), intent(in) :: s, ds
  real(rdef) :: f, p_z, d2, alpha, dxv, dyv, ds2_f, charge
  real(rdef) :: r(3,3), w(3), ex, ey, ex2, ey2, exy, bz, bz2, m2c4

!

  charge = charge_of(param%particle)
  m2c4 = (mass_of(param%particle) / param%energy)**2

  end = start

! 1) Push the position 1/2 step

  p_z = sqrt((1 + end%z%vel)**2 - end%x%vel**2 - end%y%vel**2 - m2c4)
  ds2_f = ds / (2 * p_z)

  end%x%pos = end%x%pos + end%x%vel * ds2_f
  end%y%pos = end%y%pos + end%y%vel * ds2_f
  end%z%pos = end%z%pos - ds2_f + ds / (2 * (1 + end%z%vel)) 

! 2) Evaluate the fields .
! 3) Push the momenta a 1/2 step using only "b".

  call em_field (ele, param, s+ds/2, end, field)

  f = ds * charge * c_light / (2 * 1e9 * param%energy)

  end%x%vel = end%x%vel - field%b(2) * f
  end%y%vel = end%y%vel + field%b(1) * f

  if (field%e(3) /= 0) then
    end%z%vel = end%z%vel + field%e(3) * f / c_light
    p_z = sqrt((1 + end%z%vel)**2 - end%x%vel**2 - end%y%vel**2 - m2c4)
  endif

! 4) Push the momenta a full step using "R".

  d2 = ds * charge * c_light / (2 * 1e9 * p_z * param%energy) 

  if (field%e(1) == 0 .and. field%e(2) == 0) then
    if (field%b(3) /= 0) then
      d2 = d2 * field%b(3)
      alpha = 2 * d2 / (1 + d2**2)
      dxv = -d2 * end%x%vel + end%y%vel
      dyv = -end%x%vel - d2 * end%y%vel
      end%x%vel = end%x%vel + alpha * dxv
      end%y%vel = end%y%vel + alpha * dyv
    endif
  else
    ex = field%e(1) / c_light;     ex2 = ex**2
    ey = field%e(2) / c_light;     ey2 = ey**2
    bz = field%b(3);               bz2 = bz**2
    exy = ex * ey
    alpha = 2 * d2 / (1 + d2**2 * (bz2 - ex2 - ey2))
    r(1,1:3) = (/ d2 * (ex2 - bz2), bz2 + d2*exy,     bz*(ex + d2*ey)  /)
    r(2,1:3) = (/ -bz2 + d2*exy,    d2 * (ey2 - bz2), bz*(ey - d2*ex)  /)
    r(3,1:3) = (/ bz*(ex - d2*ey),  bz*(ey + d2*ex),  d2 * (ex2 + ey2) /)
    w = end%vec(2:6:2)
    w = w + matmul(r, w)
    end%vec(2:6:2) = w
  endif

! 5) Push the momenta a 1/2 step using only "b"

  end%x%vel = end%x%vel - field%b(2) * f
  end%y%vel = end%y%vel + field%b(1) * f
  end%z%vel = end%z%vel + field%e(3) * f / c_light

! 6) Push the position a 1/2 step.

  p_z = sqrt((1 + end%z%vel)**2 - end%x%vel**2 - end%y%vel**2 - m2c4)
  ds2_f = ds / (2 * p_z)

  end%x%pos = end%x%pos + end%x%vel * ds2_f
  end%y%pos = end%y%pos + end%y%vel * ds2_f
  end%z%pos = end%z%pos - ds2_f + ds / (2 * (1 + end%z%vel)) 

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine em_field (ele, param, s, here, field)

  implicit none

  type (ele_struct), intent(in) :: ele
  type (param_struct) param
  type (coord_struct), intent(in) :: here
  type (em_field_struct), intent(out) :: field

  real(rdef) :: s

!

  if (ele%tracking_method == custom$) then
    call em_field_custom (ele, param, s, here, field)
  else
    call em_field_standard (ele, param, s, here, field)
  endif

end subroutine

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------

subroutine allocate_saved_orbit (n_steps)

  implicit none

  integer n_steps, nn

!

  nn = 2 + n_steps

  if (associated(track_com%s)) then
    if (size(track_com%s) < nn) then
      deallocate(track_com%s, track_com%orb)
      allocate(track_com%s(nn), track_com%orb(nn))
    endif
  else
    allocate(track_com%s(nn), track_com%orb(nn))
  endif

  track_com%n_pts = 0

end subroutine

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------

subroutine save_a_step (ele, param, s, here, s_sav)

  implicit none

  type (ele_struct), intent(in) :: ele
  type (param_struct), intent(in) :: param
  type (coord_struct), intent(in) :: here
  type (coord_struct) orb
  integer n_pts
  real(rdef) s, s_sav

!

  track_com%n_pts = track_com%n_pts + 1
  n_pts = track_com%n_pts 

  if (n_pts > size(track_com%s)) then
    print *, 'ERROR IN SAVE_A_STEP: ARRAY OVERFLOW!'
    call err_exit
  end if

  orb = here
  call offset_particle (ele, param, orb, unset$, set_canonical = .false., s_pos = s)

  track_com%s(n_pts) = s
  track_com%orb(n_pts) = orb
  s_sav = s

end subroutine save_a_step

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine track_solenoid_edge (ele, param, set, orb)

  implicit none

  type (ele_struct), intent(in) :: ele
  type (param_struct), intent(in) :: param
  type (coord_struct) :: orb

  logical, intent(in) :: set

!

  if (ele%key /= solenoid$ .and. ele%key /= sol_quad$) return

  if (set == set$) then

    orb%x%vel = orb%x%vel + orb%y%pos * ele%value(ks$) / 2
    orb%y%vel = orb%y%vel - orb%x%pos * ele%value(ks$) / 2

  else

    orb%x%vel = orb%x%vel - orb%y%pos * ele%value(ks$) / 2
    orb%y%vel = orb%y%vel + orb%x%pos * ele%value(ks$) / 2

  endif

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine em_field_standard (ele, param, s_pos, here, field)
!
! Subroutine to calculate the E and B fields for an element.
!
! Note: The position and fields are calculated in the frame of referene of
!   the element. Not the Laboratory frame.
!
! The variables in Boris tracking are:
!     here%vec = (x, p_x, y, p_y, s_pos-c*t, dE/E)
! At high energy s-c*t = z which is the distance of the particle from the 
! reference particle.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element
!   param  -- Param_struct: Parameters
!   s_pos  -- Real(rdef): Longitudinal position.
!   here   -- Coord_struct: Transverse coordinates.
!
! Output:
!   field -- em_field_struct: E and B fields
!-


subroutine em_field_standard (ele, param, s_pos, here, field)

  implicit none

  type (ele_struct), target, intent(in) :: ele
  type (param_struct) param
  type (coord_struct), intent(in) :: here
  type (wig_term_struct), pointer :: t
  type (em_field_struct), intent(out) :: field

  real(rdef) :: x, y, s, s_pos, f, c
  real(rdef) :: c_x, s_x, c_y, s_y, c_z, s_z, coef

  integer i

!

  x = here%x%pos
  y = here%y%pos
  s = s_pos

  field%e = 0
  field%b = 0

! Wiggler

  select case (ele%key)

  case(wiggler$)

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

      field%b(1) = field%b(1) + -coef  * (t%kx / t%ky) * s_x * s_y * c_z
      field%b(2) = field%b(2) +  coef  *                 c_x * c_y * c_z
      field%b(3) = field%b(3) + -coef  * (t%kz / t%ky) * c_x * s_y * s_z
    enddo

! Elseparator

  case (elseparator$)

    f = 1e9 * param%energy / ele%value(l$)

    field%e(1) = ele%value(hkick$) * f
    field%e(2) = ele%value(vkick$) * f

    if (ele%value(tilt$) /= 0) then
      c = cos(ele%value(tilt$))
      s = sin(ele%value(tilt$))
      field%e(1:2) = &
          (/ c*field%e(1) + s*field%e(2), -s*field%e(1) + c*field%e(2) /)
    endif

! Quadrupole

  case (quadrupole$) 

    f = 1e9 * param%energy / c_light

    field%b(1) = y * ele%value(k1$) * f 
    field%b(2) = x * ele%value(k1$) * f 

! Sol_quad

  case (sol_quad$)

    f = 1e9 * param%energy / c_light

    field%b(1) = y * ele%value(k1$) * f 
    field%b(2) = x * ele%value(k1$) * f 
    field%b(3) = ele%value(ks$) * f

! Solenoid

  case (solenoid$)

    f = 1e9 * param%energy / c_light

    field%b(3) = ele%value(ks$) * f

! Error

  case default
    print *, 'ERROR IN EM_FIELD_STANDARD: ELEMENT NOT CODED: ', &
                                                         key_name(ele%key)
    print *, '      FOR: ', ele%name
    call err_exit
  end select

end subroutine

end module
