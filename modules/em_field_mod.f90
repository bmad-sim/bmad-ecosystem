!+
! Module em_field_mod
!
! Module to define the electric and magnetic fields for an elemet.
!-

#include "CESR_platform.inc"      

module em_field_mod

  use bmad_struct
  use bmad_interface

! track_com is the common block track variable.

  type (track_struct), save :: track_com

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
!+
! Subroutine allocate_saved_orbit (track, n_pt)
!
! Subroutine to allocate space for the track structure.
!
! Input:
!   track -- Track_struct: Structure to initialize.
!   n_pt  -- Integer: Upper bound of track%pt(1:n).
!
! Output:
!   track -- Track_struct: structure for holding the track
!     %pt(1:n)  -- Track_point_struct: n will be at least n_pt
!     %n_bad    -- Reset to 0
!     %n_ok     -- Reset to 0
!     %n_pt     -- Reset to -1
!-

subroutine allocate_saved_orbit (track, n_pt)

  implicit none

  type (track_struct) track
  integer n_pt

!

  if (.not. associated (track%pt)) allocate(track%pt(0:n_pt))
  if (ubound(track%pt, 1) < n_pt) then
    deallocate(track%pt)
    allocate(track%pt(0:n_pt))
  endif

  track%n_ok = 0
  track%n_bad = 0
  track%n_pt = -1

end subroutine

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------

subroutine save_a_step (track, ele, param, s, here, s_sav)

  implicit none

  type (track_struct) track
  type (ele_struct), intent(in) :: ele
  type (param_struct), intent(in) :: param
  type (coord_struct) orb
  integer n_pt
  real(rp) s, s_sav, here(:)

!

  track%n_pt = track%n_pt + 1
  n_pt = track%n_pt 

  if (n_pt > ubound(track%pt, 1)) then
    print *, 'ERROR IN SAVE_A_STEP: ARRAY OVERFLOW!'
    call err_exit
  end if

  orb%vec = here
  call offset_particle (ele, param, orb, unset$, set_canonical = .false., s_pos = s)

  track%pt(n_pt)%s = s
  track%pt(n_pt)%orb = orb
  track%pt(n_pt)%mat6 = 0
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
  real(rp) phase, gradient, dEz_dz, theta, phi, r, E_r, B_phi

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
  if (ele%value(x_offset_tot$) /= 0 .or. ele%value(y_offset_tot$) /= 0 .or. &
       ele%value(x_pitch_tot$) /= 0 .or. ele%value(y_pitch_tot$) /= 0) offset = .true.

  if (offset) then
    s_rel = s_pos - ele%value(l$) / 2  ! position relative to center.
    x = x - ele%value(x_offset_tot$) - ele%value(x_pitch_tot$) * s_rel
    y = y - ele%value(y_offset_tot$) - ele%value(y_pitch_tot$) * s_rel
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

    f = ele%value(p0c$) / c_light
    field%b(1) = y * ele%value(k1$) * f 
    field%b(2) = x * ele%value(k1$) * f 

    if (df_calc) then
      field%dB = 0
      field%dB(1,1) = -ele%value(k1$) * f
      field%dB(2,2) =  ele%value(k1$) * f
    endif

!------------------------------------------
! Sextupole 

  case (sextupole$)

    f = ele%value(p0c$) / c_light
    field%b(1) = x * y * ele%value(k2$) * f
    field%b(2) = 1.0
    field%b(2) = (1/2.0) * ele%value(k2$) * f * (x**2 - y**2)

    if (df_calc) then
      print *, 'ERROR IN EM_FIELD: dFIELD NOT YET IMPLEMENTED FOR SEXTUPOLE!'
      call err_exit
    endif

!------------------------------------------
! Sol_quad

  case (sol_quad$)

    f = ele%value(p0c$) / c_light
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

    f = ele%value(p0c$) / c_light
    field%b(3) = ele%value(ks$) * f

    if (df_calc) then
      print *, 'ERROR IN EM_FIELD: dFIELD NOT YET IMPLEMENTED FOR SOLENOID!'
      call err_exit
    endif

!------------------------------------------
! Lcavity
!
! This is taken from the gradient as calculated in
!       J. Rosenzweig and L. Serafini
!       Phys Rev E, Vol. 49, p. 1599, (1994)
!
! Right now only works at relativistic energies

   case (lcavity$)

     !***
     ! This is taken from track1_bmad
     phase = twopi * (ele%value(phi0$) + ele%value(dphi0$) + ele%value(phi0_err$) - &
                        here%vec(5) * ele%value(rf_frequency$) / c_light)
     gradient = (ele%value(gradient$) + ele%value(gradient_err$)) * cos(phase)
     if (.not. ele%is_on) gradient = 0
     
     if (bmad_com%sr_wakes_on) then
       if (bmad_com%grad_loss_sr_wake /= 0) then  
         ! use grad_loss_sr_wake and ignore e_loss
         gradient = gradient - bmad_com%grad_loss_sr_wake
       else
         gradient = gradient - ele%value(e_loss$) * param%n_part * &
                                                     e_charge / ele%value(l$)
       endif
     endif
     !***

     gradient =  charge_of (param%particle) * gradient
     ! This gives the averave gradient, I need the maximum gradient  
     gradient = (pi/2.0) * gradient
                                                                     
     ! only use first pi phase of sine wave                                
     s = modulo (s, c_light/(2.0*ele%value(rf_frequency$)))          
                                                                     
                                                                     
     f = s * twopi * ele%value(rf_frequency$) / c_light              
     dEz_dz = (f/s) * gradient * cos (f + phase)                     
     theta = atan(y/x)                                               
     r = sqrt(x**2 + y**2)                                           
     E_r =  - (r/2.0) * dEz_dz                                       
     B_phi = (r/(2.0*c_light)) * dEz_dz                              
                                                                     
                                                                     
     field%E(1) = E_r * cos (theta)                                  
     field%E(2) = E_r * sin (theta)
     field%E(3) = gradient * sin (f + phase)
    
     phi = pi - theta
     field%B(1) =   B_phi * cos (phi)
     field%B(2) = - B_phi * sin (phi)

    if (df_calc) then
      print *, 'ERROR IN EM_FIELD: dFIELD NOT YET IMPLEMENTED FOR LCAVITY!'
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
      field%dB(:,2) = dk(:,1) * sin_ang + dk(:,2) * cos_ang
      field%dB(:,3) = dk(:,3) 

      dk(1,:) = cos_ang * field%dE(1,:) - sin_ang * field%dE(2,:)
      dk(2,:) = sin_ang * field%dE(1,:) + cos_ang * field%dE(2,:)
      dk(3,:) = field%dE(3,:)

      field%dE(:,1) = dk(:,1) * cos_ang - dk(:,2) * sin_ang
      field%dE(:,2) = dk(:,1) * sin_ang + dk(:,2) * cos_ang
      field%dE(:,3) = dk(:,3) 

    endif
  endif

  if (offset) then
    field%B(1) = field%B(1) + ele%value(x_pitch_tot$) * field%B(3)
    field%B(2) = field%B(2) + ele%value(y_pitch_tot$) * field%B(3)
    field%E(1) = field%E(1) + ele%value(x_pitch_tot$) * field%E(3)
    field%E(2) = field%E(2) + ele%value(y_pitch_tot$) * field%E(3)
  endif

end subroutine

end module
