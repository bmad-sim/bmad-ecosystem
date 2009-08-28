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
    subroutine em_field_custom (ele, param, s, orb, local_ref_frame, field, calc_dfield)
      use bmad_struct
      implicit none
      type (ele_struct) :: ele
      type (lat_param_struct) param
      type (coord_struct), intent(in) :: orb
      real(rp), intent(in) :: s
      logical local_ref_frame
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
  type (ele_struct) :: ele
  type (lat_param_struct), intent(in) :: param
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
! Subroutine em_field_calc (ele, param, s_pos, here, local_ref_frame, field, calc_dfield)
!
! Subroutine to calculate the E and B fields for an element.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element
!   param  -- lat_param_struct: Lattice parameters.
!   s_pos  -- Real(rp): Longitudinal position relative to the start of the element.
!   here   -- Coord_struct: Transverse coordinates.
!   local_ref_frame 
!          -- Logical, If True then take the input coordinates and output fields 
!                as being with respect to the frame of referene of the element. 
!   calc_dfield     
!         -- Logical, optional: If present and True 
!                then calculate the field derivatives.
!
! Output:
!   field -- em_field_struct: E and B fields
!-

subroutine em_field_calc (ele, param, s_pos, here, local_ref_frame, field, calc_dfield)

  implicit none

  type (ele_struct), target :: ele
  type (lat_param_struct) param
  type (coord_struct) :: here, local_here
  type (wig_term_struct), pointer :: t
  type (em_field_struct), intent(out) :: field

  real(rp) :: x, y, xx, yy, s, s_pos, f, dk(3,3), charge
  real(rp) :: c_x, s_x, c_y, s_y, c_z, s_z, coef, fd(3)
  real(rp) :: cos_ang, sin_ang, s_rel, sgn_x, dc_x, dc_y
  real(rp) phase, gradient, dEz_dz, theta, phi, r, E_r, B_phi

  integer i, sign_charge

  logical :: local_ref_frame
  logical, optional :: calc_dfield

  logical df_calc
  character(20) :: r_name = 'em_field_calc'

! custom field_calc

  select case (ele%field_calc)
  case (custom$) 
    call em_field_custom (ele, param, s_pos, here, local_ref_frame, field, calc_dfield)
    return
  case (bmad_standard$)
  case default
    call out_io (s_fatal$, r_name, 'BAD FIELD_CALC METHOD FOR ELEMENT: ' // ele%name)
    call err_exit
  end select

  charge = charge_of(param%particle)
  sign_charge = sign(1.0_rp, charge)
  
!----------------------------------------------------------------------------
! convert to local coords

  local_here = here
  if (.not. local_ref_frame) then
    call offset_particle (ele, param, local_here, set$, &
            set_canonical = .false., set_multipoles = .false., set_hvkicks = .false.)
  endif

  x = local_here%vec(1)
  y = local_here%vec(3)
  s = s_pos

! Init

  field%e = 0
  field%B = 0
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
      print *, 'ERROR IN EM_FIELD_CALC: PERIODIC WIGGLER NOT YET IMPLEMENTED!'
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
        f = coef * t%kz
        field%dB(1,3) = field%dB(1,3) + f  * (t%kx / t%ky) * s_x * s_y * s_z * sgn_x
        field%dB(2,3) = field%dB(2,3) - f  *                 c_x * c_y * s_z * dc_y
        field%dB(3,3) = field%dB(3,3) - f  * (t%kz / t%ky) * c_x * s_y * c_z 
      endif

    enddo

!------------------------------------------
! Drift, et. al.

  case (drift$, ecollimator$, rcollimator$, instrument$, monitor$, pipe$)

!------------------------------------------
! Quadrupole

  case (quadrupole$) 

    f = ele%value(p0c$) / c_light
    field%b(1) = y * sign_charge * ele%value(k1$) * f 
    field%b(2) = x * sign_charge * ele%value(k1$) * f 

    if (df_calc) then
      field%dB = 0
      field%dB(1,1) =  ele%value(k1$) * f
      field%dB(2,2) = -ele%value(k1$) * f
    endif

!------------------------------------------
! Sextupole 

  case (sextupole$)

    f = ele%value(p0c$) / c_light
    field%b(1) = x * y * sign_charge * ele%value(k2$) * f
    field%b(2) = 1.0
    field%b(2) = - (1/2.0) * sign_charge * ele%value(k2$) * f * (x**2 - y**2)

    if (df_calc) then
      print *, 'ERROR IN EM_FIELD_CALC: dFIELD NOT YET IMPLEMENTED FOR SEXTUPOLE!'
      call err_exit
    endif

!------------------------------------------
! Sol_quad

  case (sol_quad$)

    f = ele%value(p0c$) / c_light
    field%b(1) = y * sign_charge * ele%value(k1$) * f 
    field%b(2) = x * sign_charge * ele%value(k1$) * f 
    field%b(3) = sign_charge * ele%value(ks$) * f

    if (df_calc) then
      print *, 'ERROR IN EM_FIELD_CALC: dFIELD NOT YET IMPLEMENTED FOR SOL_QUAD!'
      call err_exit
    endif

!------------------------------------------
! Solenoid

  case (solenoid$)

    f = ele%value(p0c$) / c_light
    field%b(3) = sign_charge * ele%value(ks$) * f

    if (df_calc) then
      print *, 'ERROR IN EM_FIELD_CALC: dFIELD NOT YET IMPLEMENTED FOR SOLENOID!'
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
     if (x .eq. 0.0) then
       theta = 0.0
     else
       theta = atan(y/x)   
     endif
     r = sqrt(x**2 + y**2)                                           
     E_r =  - (r/2.0) * dEz_dz                                       
     B_phi = (r/(2.0*(c_light**2))) * dEz_dz                              
                                                                     
                                                                     
     field%E(1) = E_r * cos (theta)                                  
     field%E(2) = E_r * sin (theta)
     field%E(3) = gradient * sin (f + phase)
    
     phi = pi - theta
     field%B(1) =   B_phi * cos (phi)
     field%B(2) = - B_phi * sin (phi)

    if (df_calc) then
      print *, 'ERROR IN EM_FIELD_CALC: dFIELD NOT YET IMPLEMENTED FOR LCAVITY!'
      call err_exit
    endif

!------------------------------------------
! Error

  case default
    print *, 'ERROR IN EM_FIELD_CALC: ELEMENT NOT YET CODED: ', key_name(ele%key)
    print *, '      FOR: ', ele%name
    call err_exit
  end select

!----------------------
! convert fields to lab coords

  if (ele%value(tilt_tot$) /= 0) then

    sin_ang = sin(ele%value(tilt_tot$))
    cos_ang = cos(ele%value(tilt_tot$))

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

  !

  field%B(1) = field%B(1) + ele%value(x_pitch_tot$) * field%B(3)
  field%B(2) = field%B(2) + ele%value(y_pitch_tot$) * field%B(3)
  field%E(1) = field%E(1) + ele%value(x_pitch_tot$) * field%E(3)
  field%E(2) = field%E(2) + ele%value(y_pitch_tot$) * field%E(3)

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine em_field_kick (ele, param, s, r, local_ref_frame, dr_ds, dkick)
!
! Subroutine to essentially calculate the kick felt by a particle in a
! element. 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele   -- Ele_struct: Element being tracked thorugh.
!   param -- lat_param_struct: Lattice parameters.
!   s     -- Real(rp): Distance from the start of the element to the particle.
!   r(6)  -- Real(rp): Position vector: (x, x', y, y', z, Pz)
!   local_ref_frame 
!         -- Logical, If True then take the input coordinates and output fields 
!               as being with respect to the frame of referene of the element. 
!
! Output:
!   dr_ds(6) -- Real(rp):Kick vector: 
!                 (dx/ds, dx'/ds, dy/ds, dy'/ds, dz/ds, dPz/ds)
!   dkick(3,3) -- Real(rp), optional: dKick/dx
!-

subroutine em_field_kick (ele, param, s, r, local_ref_frame, dr_ds, dkick)

  implicit none

  type (ele_struct) ele
  type (lat_param_struct) param
  type (em_field_struct) field

  real(rp), intent(in) :: s         ! s-position
  real(rp), intent(in) :: r(6)      ! (x, x', y, y', z, z')
  real(rp), intent(out) :: dr_ds(6)
  real(rp), optional :: dkick(3,3)

  type (coord_struct) here

  real(rp) energy
  real(rp) vel_x, vel_y, vel_s, dvel_x, dvel_y, dvel_s, f

  logical :: local_ref_frame

! calculate the field

  call init_coord(here, r)
  if (present (dkick)) then
    call em_field_calc (ele, param, s, here, local_ref_frame, field, .true.)
  else
    call em_field_calc (ele, param, s, here, local_ref_frame, field, .false.)
  endif

! if this is a kick field then field gives us directly dr_ds

  if (field%type == kick_field$) then
    dr_ds(1) = r(2)    ! dx/ds =
    dr_ds(2) = field%kick(1)
    dr_ds(3) = r(4)
    dr_ds(4) = field%kick(2)
    dr_ds(5) = -(r(2)**2 + r(4)**2) / 2
    dr_ds(6) = field%kick(3)
    return
  endif

! Here for field%type = em_field
! The computation (up to some constant factors):
!     x' = dx/ds = v_x / v_s
!     dx'/ds = (dv_x/dt * v_s - v_s * dv_s/dt) / v_s^3  ! ds/dt == v_s
! where
!   dv_x/dt = v_y * B_s - v_s * B_y
!   dv_y/dt = v_s * B_x - v_x * B_s
!   dv_s/dt = v_x * B_y - v_y * B_x

  if (field%type == em_field$) then
    vel_x = r(2)                              ! proportional to x-velosity
    vel_y = r(4)                              ! proportional to y-velosity
    vel_s = 1/sqrt(1 + r(2)**2 + r(4)**2)     ! proportional to s-velosity

    f = param%particle / c_light
    dvel_x = vel_y * field%B(3) - vel_s * field%B(2) + field%E(1) * f
    dvel_y = vel_s * field%B(1) - vel_x * field%B(3) + field%E(2) * f
    dvel_s = vel_x * field%B(2) - vel_y * field%B(1)

    energy = (ele%value(E_TOT$) * (1 + r(6)))
    f = c_light / energy

    dr_ds(1) = r(2)
    dr_ds(2) = f * (dvel_x * vel_s - vel_x * dvel_s) / vel_s**3
    dr_ds(3) = r(4)
    dr_ds(4) = f * (dvel_y * vel_s - vel_y * dvel_s) / vel_s**3
    dr_ds(5) = -(r(2)**2 + r(4)**2) / 2
    dr_ds(6) = field%E(3) / energy

! We make the small angle approximation for the dkick calc

    if (present(dkick)) then
      dkick = 0
      f = c_light / (energy * vel_s**2)
      dkick(1,1) = f * (vel_y * field%dB(3,1) - vel_s * field%dB(2,1)) 
      dkick(1,2) = f * (vel_y * field%dB(3,2) - vel_s * field%dB(2,2)) 
      dkick(2,1) = f * (vel_s * field%dB(1,1) - vel_x * field%dB(3,1)) 
      dkick(2,2) = f * (vel_s * field%dB(1,2) - vel_x * field%dB(3,2)) 
    endif

    return
  endif
  
  print *, 'ERROR IN DERIVS: UNKNOWN FIELD_TYPE: ', field%type
  call err_exit

end subroutine

end module
