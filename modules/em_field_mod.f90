!+
! Module em_field_mod
!
! Module to define the electric and magnetic fields for an elemet.
!-

#include "CESR_platform.inc"      

module em_field_mod

  use bmad_struct
  use bmad_interface

! for tracking integration

  type track_com_struct
    real(rp), pointer :: s(:) => null()    ! s-distance of a point
    type (coord_struct), allocatable :: orb(:) ! position of a point
    real(rp) :: ds_save = 1e-3             ! min distance between points
    real(rp) :: step0 = 1e-3               ! Initial step size.
    real(rp) :: step_min = 1e-8            ! min step size to step below which
                                           !   track1_adaptive_boris gives up
    integer :: max_step = 10000            ! maximum number of steps allowed
    logical :: save_track = .false.        ! save orbit?
    integer :: n_pts                       ! number of points
    integer :: n_bad
    integer :: n_ok
  end type

  type (track_com_struct), save :: track_com

! Interface for custom field calc

  interface 
    subroutine em_field_custom (ele, param, s, orb, field, calc_dfield)
      use bmad_struct
      implicit none
      type (ele_struct), intent(in) :: ele
      type (param_struct) param
      type (coord_struct), intent(in) :: orb
      real(rp), intent(in) :: s
      type (em_field_struct), intent(out) :: field
      logical, optional :: calc_dfield
    end subroutine
  end interface

contains

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

  track_com%n_ok = 0
  track_com%n_bad = 0
  track_com%n_pts = 0

end subroutine

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------

subroutine save_a_step (ele, param, s, here, s_sav)

  implicit none

  type (ele_struct), intent(in) :: ele
  type (param_struct), intent(in) :: param
  type (coord_struct) orb
  integer n_pts
  real(rp) s, s_sav, here(:)

!

  track_com%n_pts = track_com%n_pts + 1
  n_pts = track_com%n_pts 

  if (n_pts > size(track_com%s)) then
    print *, 'ERROR IN SAVE_A_STEP: ARRAY OVERFLOW!'
    call err_exit
  end if

  orb%vec = here
  call offset_particle (ele, param, orb, unset$, set_canonical = .false., s_pos = s)

  track_com%s(n_pts) = s
  track_com%orb(n_pts) = orb
  s_sav = s

end subroutine save_a_step

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine em_field (ele, param, s_pos, here, field, calc_dfield)
!
! Subroutine to calculate the E and B fields for an element
! in the local frame of reference.
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
!   s_pos  -- Real(rp): Longitudinal position.
!   here   -- Coord_struct: Transverse coordinates.
!   calc_dfield -- Optional, logical: If present 
!        and True then calculate the field derivatives.
!
! Output:
!   field -- em_field_struct: E and B fields
!-

subroutine em_field (ele, param, s_pos, here, field, calc_dfield)

  implicit none

  type (ele_struct), target, intent(in) :: ele
  type (param_struct) param
  type (coord_struct) :: here
  type (wig_term_struct), pointer :: t
  type (em_field_struct), intent(out) :: field

  real(rp) :: x, y, s, s_pos, f, dk(3,3)
  real(rp) :: c_x, s_x, c_y, s_y, c_z, s_z, coef, fd(3)
  real(rp) :: cos_ang, sin_ang, s_rel, sgn_x, dc_x, dc_y

  integer i

  logical, optional :: calc_dfield
  logical offset, df_calc
  character(20) :: r_name = 'em_field'

! custom field_calc

  select case (ele%field_calc)
  case (custom$) 
    call em_field_custom (ele, param, s_pos, here, field, calc_dfield)
    return
  case (bmad_standard$)
  case default
    call out_io (s_fatal$, r_name, 'BAD FIELD_CALC METHOD FOR ELEMENT: ' // ele%name)
    call err_exit
  end select

!----------------------------------------------------------------------------
! convert to local coords

  x = here%vec(1)
  y = here%vec(3)
  s = s_pos

  offset = .false.
  if (ele%value(x_offset$) /= 0 .or. ele%value(y_offset$) /= 0 .or. &
       ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0) offset = .true.

  if (offset) then
    s_rel = s_pos - ele%value(l$) / 2  ! position relative to center.
    x = x - ele%value(x_offset$) - ele%value(x_pitch$) * s_rel
    y = y - ele%value(y_offset$) - ele%value(y_pitch$) * s_rel
  endif

  if (ele%value(tilt_tot$) /= 0) then
    cos_ang = cos(ele%value(tilt_tot$))
    sin_ang = sin(ele%value(tilt_tot$))
    x =  cos_ang * x + sin_ang * y
    y = -sin_ang * x + cos_ang * y
  endif

  field%e = 0
  field%b = 0
  field%type = em_field$

  df_calc = .false.
  if (present(calc_dfield)) df_calc = calc_dfield

  if (df_calc) then
    field%dB = 0
    field%dE = 0
  endif

!------------------------------------------

  select case (ele%key)

!------------------------------------------
! Wiggler

  case(wiggler$)

    if (ele%sub_key /= map_type$) then
      print *, 'ERROR IN EM_FIELD: PERIODIC WIGGLER NOT YET IMPLEMENTED!'
      call err_exit
    endif

    do i = 1, size(ele%wig_term)
      t => ele%wig_term(i)

      if (t%type == hyper_y$) then
        c_x = cos(t%kx * x)
        s_x = sin(t%kx * x)
        sgn_x = 1
        dc_x = -1
      else
        c_x = cosh(t%kx * x)
        s_x = sinh(t%kx * x)
        sgn_x = -1
        dc_x = 1
      endif

      if (t%type == hyper_y$ .or. t%type == hyper_xy$) then
        c_y = cosh (t%ky * y)
        s_y = sinh (t%ky * y)
        dc_y = 1
      else
        c_y = cos (t%ky * y)
        s_y = sin (t%ky * y)
        dc_y = -1
      endif

      c_z = cos (t%kz * s + t%phi_z)
      s_z = sin (t%kz * s + t%phi_z)

      coef = t%coef * ele%value(polarity$)

      field%B(1) = field%B(1) - coef  * (t%kx / t%ky) * s_x * s_y * c_z * sgn_x
      field%B(2) = field%B(2) + coef  *                 c_x * c_y * c_z
      field%B(3) = field%B(3) - coef  * (t%kz / t%ky) * c_x * s_y * s_z

      if (df_calc) then
        f = coef * t%kx
        field%dB(1,1) = field%dB(1,1) - f  * (t%kx / t%ky) * c_x * s_y * c_z * sgn_x
        field%dB(2,1) = field%dB(2,1) + f  *                 s_x * c_y * c_z * dc_x
        field%dB(3,1) = field%dB(3,1) - f  * (t%kz / t%ky) * s_x * s_y * s_z * dc_x
        f = coef * t%ky
        field%dB(1,2) = field%dB(1,2) - f  * (t%kx / t%ky) * s_x * c_y * c_z * sgn_x
        field%dB(2,2) = field%dB(2,2) + f  *                 c_x * s_y * c_z * dc_y
        field%dB(3,2) = field%dB(3,2) - f  * (t%kz / t%ky) * c_x * c_y * s_z 
      endif

    enddo

!------------------------------------------
! Drift, et. al.

  case (drift$, ecollimator$, rcollimator$, instrument$, monitor$)

!------------------------------------------
! Quadrupole

  case (quadrupole$) 

    f = ele%value(beam_energy$) / c_light
    field%b(1) = y * ele%value(k1$) * f 
    field%b(2) = x * ele%value(k1$) * f 

    if (df_calc) then
      field%dB = 0
      field%dB(1,1) = -ele%value(k1$) * f
      field%dB(2,2) =  ele%value(k1$) * f
    endif

!------------------------------------------
! Sol_quad

  case (sol_quad$)

    f = ele%value(beam_energy$) / c_light
    field%b(1) = y * ele%value(k1$) * f 
    field%b(2) = x * ele%value(k1$) * f 
    field%b(3) = ele%value(ks$) * f

    if (df_calc) then
      print *, 'ERROR IN EM_FIELD: dFIELD NOT YET IMPLEMENTED FOR SOL_QUAD!'
      call err_exit
    endif

!------------------------------------------
! Solenoid

  case (solenoid$)

    f = ele%value(beam_energy$) / c_light
    field%b(3) = ele%value(ks$) * f

    if (df_calc) then
      print *, 'ERROR IN EM_FIELD: dFIELD NOT YET IMPLEMENTED FOR SOLENOID!'
      call err_exit
    endif

!------------------------------------------
! Error

  case default
    print *, 'ERROR IN EM_FIELD: ELEMENT NOT YET CODED: ', key_name(ele%key)
    print *, '      FOR: ', ele%name
    call err_exit
  end select

!----------------------
! convert fields to lab coords

  if (ele%value(tilt_tot$) /= 0) then

    fd = field%B
    field%B(1) = cos_ang * fd(1) - sin_ang * fd(2)
    field%B(2) = sin_ang * fd(1) + cos_ang * fd(2)

    fd = field%E
    field%E(1) = cos_ang * fd(1) - sin_ang * fd(2)
    field%E(2) = sin_ang * fd(1) + cos_ang * fd(2)

    if (df_calc) then

      dk(1,:) = cos_ang * field%dB(1,:) - sin_ang * field%dB(2,:)
      dk(2,:) = sin_ang * field%dB(1,:) + cos_ang * field%dB(2,:)
      dk(3,:) = field%dB(3,:)

      field%dB(:,1) = dk(:,1) * cos_ang - dk(:,2) * sin_ang
      field%dB(:,2) = dk(:,1) * sin_ang - dk(:,2) * cos_ang
      field%dB(:,3) = dk(:,3) 

      dk(1,:) = cos_ang * field%dE(1,:) - sin_ang * field%dE(2,:)
      dk(2,:) = sin_ang * field%dE(1,:) + cos_ang * field%dE(2,:)
      dk(3,:) = field%dE(3,:)

      field%dE(:,1) = dk(:,1) * cos_ang - dk(:,2) * sin_ang
      field%dE(:,2) = dk(:,1) * sin_ang - dk(:,2) * cos_ang
      field%dE(:,3) = dk(:,3) 

    endif
  endif

  if (offset) then
    field%B(1) = field%B(1) + ele%value(x_pitch$) * field%B(3)
    field%B(2) = field%B(2) + ele%value(y_pitch$) * field%B(3)
    field%E(1) = field%E(1) + ele%value(x_pitch$) * field%E(3)
    field%E(2) = field%E(2) + ele%value(y_pitch$) * field%E(3)
  endif

end subroutine

end module
