!+
! Module em_field_mod
!
! Module to define the electric and magnetic fields for an elemet.
!-

module em_field_mod

use bmad_struct
use bmad_interface

! Interface for custom field calc

interface 
  subroutine em_field_custom (ele, param, s_rel, t_rel, orb, local_ref_frame, field, calc_dfield)
    use bmad_struct
    implicit none
    type (ele_struct) :: ele
    type (lat_param_struct) param
    type (coord_struct), intent(in) :: orb
    real(rp), intent(in) :: s_rel, t_rel
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
! Subroutine init_saved_orbit (track, n_pt)
!
! Subroutine to initialize the track structure.
!
! Input:
!   track -- Track_struct: Structure to initialize.
!   n_pt  -- Integer: Upper bound of track%orb(0:n) and track%map(0:n).
!
! Output:
!   track -- Track_struct: structure for holding the track
!     %orb(1:n) -- Coord_struct: n will be at least n_pt
!     %map(1:n) -- Track_map_struct: n will be at least n_pt
!     %n_bad    -- Reset to 0
!     %n_ok     -- Reset to 0
!     %n_pt     -- Reset to -1
!-

subroutine init_saved_orbit (track, n_pt)

implicit none

type (track_struct) track
integer n_pt

!

if (.not. allocated (track%orb)) then
  allocate(track%orb(0:n_pt))
  allocate(track%map(0:n_pt))
endif

if (ubound(track%orb, 1) < n_pt) then
  deallocate(track%orb, track%map)
  allocate(track%orb(0:n_pt), track%map(0:n_pt))
endif

track%n_ok = 0
track%n_bad = 0
track%n_pt = -1

end subroutine

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! +
! Subroutine save_a_step (track, ele, param, local_ref_frame, s, here, s_sav)
!
! Routine used by Runge-Kutta and Boris tracking routines to save
! the trajectory through an element.
!
! Note: It is assumed by this routine that here(:) is the orbit in local 
! element coordinates. The actual track saved will be in laboratory coordinates.
!
! Input:
!   ele      -- ele_struct: Element being tracked through.
!   param    -- lat_param_struct: Lattice parameters.
!   local_ref_frame -- Logical: If True then coordinates are wrt the frame of ref of the element.
!   s        -- Real(rp): S-position with respect to start of element
!   orb      -- Coord_struct: trajectory at s with respect to element coordinates.
!
! Ouput:
!   track    -- track_struct: Trajectory structure to save to.
!   s_sav    -- Real(rp): Set equal to s.


subroutine save_a_step (track, ele, param, local_ref_frame, s, orb, s_sav)

implicit none

type (track_struct) track, track2
type (ele_struct) :: ele
type (lat_param_struct), intent(in) :: param
type (coord_struct) orb, orb2
integer n_pt, n, n_old
real(rp) s, s_sav
logical local_ref_frame

!

track%n_pt = track%n_pt + 1
n_pt = track%n_pt 

if (n_pt > ubound(track%orb, 1)) then
  n = 1.5 * n_pt
  n_old = ubound(track%orb, 1)
  allocate(track2%orb(0:n_old), track2%map(0:n_old))
  track2 = track
  deallocate(track%orb, track%map)
  allocate(track%orb(0:n), track%map(0:n))
  track%orb(:n_old) = track2%orb; track%map(:n_old) = track2%map
  deallocate(track2%orb, track2%map)
end if

! Notice that a translation due to a finite ele%value(s_offset$) is not wanted here.

orb2 = orb
if (local_ref_frame) call offset_particle (ele, param, orb2, unset$, &
                                    set_canonical = .false., set_s_offset = .false., s_pos = s)

track%orb(n_pt) = orb2
track%map(n_pt)%mat6 = 0
s_sav = s

end subroutine save_a_step

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine em_field_calc (ele, param, s_rel, t_rel, orbit, local_ref_frame, field, calc_dfield)
!
! Subroutine to calculate the E and B fields for an element.
!
! Note: Zero field will be returned if an element is turned off.
!
! Note: The fields due to any kicks will be present. It therefore important in tracking to make sure that 
! offset_particle does not add in kicks at the beginning and end which would result in double counting the kicks.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element
!   param  -- lat_param_struct: Lattice parameters.
!   s_rel  -- Real(rp): Longitudinal position relative to the start of the element.
!   t_rel  -- Real(rp): Time relative to the time the reference particle passed the element entrance end.
!   orbit  -- Coord_struct: Transverse coordinates.
!     %vec(1), %vec(3)  -- Transverse coords. These are the only components used in the calculation.
!   local_ref_frame 
!          -- Logical, If True then take the input coordinates and output fields 
!                as being with respect to the frame of referene of the element. 
!   calc_dfield     
!          -- Logical, optional: If present and True 
!                then calculate the field derivatives.
!
! Output:
!   field       -- em_field_struct: E and B fields and derivatives.
!-

subroutine em_field_calc (ele, param, s_rel, t_rel, orbit, local_ref_frame, field, calc_dfield)

implicit none

type (ele_struct), target :: ele
type (lat_param_struct) param
type (coord_struct) :: orbit, local_orb
type (wig_term_struct), pointer :: t
type (em_field_struct), intent(out) :: field
type (em_field_point_struct) :: local_field
type (em_field_mode_struct), pointer :: mode
type (em_field_map_term_struct), pointer :: term

real(rp) :: x, y, xx, yy, t_rel, s_rel, z,   f, dk(3,3), charge, f_p0c
real(rp) :: c_x, s_x, c_y, s_y, c_z, s_z, coef, fd(3)
real(rp) :: cos_ang, sin_ang, sgn_x, dc_x, dc_y, kx, ky, dkm(2,2)
real(rp) phase, gradient, dEz_dz, theta, r, E_r
real(rp) k_t, k_zn, kappa2_n, kap_rho, s_pos
real(rp) radius, phi, t_ref, tilt

complex(rp) E_rho, E_phi, E_z, Er, Ep, Ez, B_rho, B_phi, B_z, Br, Bp, Bz, expi, expt, dEp, dEr
complex(rp) Im_0, Im_plus, Im_minus, Im_0_R, kappa_n, Im_plus2, cm, sm

integer i, j, m, n, sign_charge

logical :: local_ref_frame
logical, optional :: calc_dfield

logical df_calc
character(20) :: r_name = 'em_field_calc'

! Initialize field
! If element is turned off then return zero
field%E = 0
field%B = 0

df_calc = logic_option (.false., calc_dfield)

if (df_calc) then
  field%dB = 0
  field%dE = 0
endif

if (.not. ele%is_on) return

!----------------------------------------------------------------------------
! Custom field calc 
!----------------------------------------------------------------------------
if (ele%field_calc == custom$) then 
  call em_field_custom (ele, param, s_rel, t_rel, orbit, local_ref_frame, field, calc_dfield)
  return
end if



!----------------------------------------------------------------------------
!Set up common variables for all (non-custom) methods

charge = charge_of(param%particle)
sign_charge = sign(1.0_rp, charge)

!----------------------------------------------------------------------------
! convert to local coords

local_orb = orbit
if (.not. local_ref_frame) then
  call offset_particle (ele, param, local_orb, set$, &
          set_canonical = .false., set_multipoles = .false., set_hvkicks = .false.)
endif

x = local_orb%vec(1)
y = local_orb%vec(3)

f_p0c = sign_charge * ele%value(p0c$) / c_light

!------------------------------------------

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! field_calc methods
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

select case (ele%field_calc)
  
!----------------------------------------------------------------------------
! Bmad_standard field calc 
!----------------------------------------------------------------------------
case (bmad_standard$)

select case (ele%key)

!------------------------------------------
! RFcavity and Lcavity

case(rfcavity$, lcavity$)

    ! This is taken from the gradient as calculated in
    !       J. Rosenzweig and L. Serafini
    !       Phys Rev E, Vol. 49, p. 1599, (1994)
    !
    ! Right now only works at relativistic energies

    ! This is taken from track1_bmad
    phase = twopi * (ele%value(phi0$) + ele%value(dphi0$) + ele%value(phi0_err$) - &
                        local_orb%vec(5) * ele%value(rf_frequency$) / c_light)
    gradient = (ele%value(gradient$) + ele%value(gradient_err$)) * cos(phase)
    if (.not. ele%is_on) gradient = 0
    gradient = gradient + gradient_shift_sr_wake(ele, param)
    
    dEz_dz = gradient * sign(1, charge_of(param%particle))

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
    field%E(3) = gradient * sin (phase)
    
    phi = pi - theta
    field%B(1) =  B_phi * cos (phi)
    field%B(2) = -B_phi * sin (phi)

    if (df_calc) then
      print *, 'ERROR IN EM_FIELD_CALC: dFIELD NOT YET IMPLEMENTED FOR LCAVITY!'
      call err_exit
    endif


!------------------------------------------
! Wiggler

case(wiggler$)

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

    c_z = cos (t%kz * s_rel + t%phi_z)
    s_z = sin (t%kz * s_rel + t%phi_z)

    coef = sign_charge * t%coef * ele%value(polarity$)

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
! Drift, et. al. Note that kicks get added at the end for all elements

case (drift$, ecollimator$, rcollimator$, instrument$, monitor$, pipe$, marker$)

!------------------------------------------
! Quadrupole

case (quadrupole$) 

  field%b(1) = y * ele%value(k1$) * f_p0c 
  field%b(2) = x * ele%value(k1$) * f_p0c 

  if (df_calc) then
    field%dB(1,1) =  ele%value(k1$) * f_p0c
    field%dB(2,2) = -ele%value(k1$) * f_p0c
  endif

!------------------------------------------
! Sextupole 

case (sextupole$)

  field%b(1) = x * y * ele%value(k2$) * f_p0c
  field%b(2) = -(x**2 - y**2) / 2 * ele%value(k2$) * f_p0c 

  if (df_calc) then
    field%dB(1,1) =  y * ele%value(k2$) * f_p0c
    field%dB(1,2) =  x * ele%value(k2$) * f_p0c
    field%dB(2,1) = -x * ele%value(k2$) * f_p0c
    field%dB(2,2) = -y * ele%value(k2$) * f_p0c
  endif

!------------------------------------------
! Sol_quad

case (sol_quad$)

  field%b(1) = y * ele%value(k1$) * f_p0c 
  field%b(2) = x * ele%value(k1$) * f_p0c 
  field%b(3) = ele%value(ks$) * f_p0c

  if (df_calc) then
    field%dB(1,1) =  ele%value(k1$) * f_p0c
    field%dB(2,2) = -ele%value(k1$) * f_p0c
  endif

!------------------------------------------
! Solenoid

case (solenoid$)

  field%b(3) = ele%value(ks$) * f_p0c

  if (df_calc) then
  endif

!------------------------------------------
! SBend

case (sbend$)

  field%b(1) = (y * ele%value(k1$) + x * y * ele%value(k2$)) * f_p0c 
  field%b(2) = (x * ele%value(k1$) - ele%value(k2$) * (x**2 - y**2) / 2 + ele%value(g$) + ele%value(g_err$)) * f_p0c 

  if (df_calc) then
    field%dB(1,1) =  ele%value(k1$) * f_p0c + y * ele%value(k2$) * f_p0c
    field%dB(1,2) =  x * ele%value(k2$) * f_p0c
    field%dB(2,1) = -x * ele%value(k2$) * f_p0c
    field%dB(2,2) = -ele%value(k1$) * f_p0c - y * ele%value(k2$) * f_p0c
  endif

!------------------------------------------
! HKicker

case (hkicker$)
  field%b(2) = -ele%value(kick$) * f_p0c 

!------------------------------------------
! VKicker

case (vkicker$)
  field%b(1) =  ele%value(kick$) * f_p0c 

!------------------------------------------
! Kicker  

case (kicker$)
  field%b(1) =  ele%value(vkick$) * f_p0c 
  field%b(2) = -ele%value(hkick$) * f_p0c 

!------------------------------------------
! Error

case default
  print *, 'ERROR IN EM_FIELD_CALC: ELEMENT NOT YET CODED: ', key_name(ele%key)
  print *, '      FOR: ', ele%name
  call err_exit
end select

!---------------------------------------------------------------------
! Add multipoles

if (associated(ele%a_pole)) then
  if (ele%value(l$) == 0) then
    print *, 'ERROR IN EM_FILED_CALC: dField NOT YET IMPLEMENTED FOR MULTIPOLES!'
    call err_exit
  endif

  do i = 0, ubound(ele%a_pole, 1)
    if (ele%a_pole(i) == 0 .and. ele%b_pole(i) == 0) cycle
    if (df_calc) then
      call ab_multipole_kick(ele%a_pole(i), ele%b_pole(i), i, local_orb, kx, ky, dkm)
    else
      call ab_multipole_kick(ele%a_pole(i), ele%b_pole(i), i, local_orb, kx, ky)
    endif
    field%B(1) = field%B(1) +  f_p0c * ky / ele%value(l$)
    field%B(2) = field%B(2) -  f_p0c * kx / ele%value(l$)
    if (df_calc) then
      field%dB(1,1) = field%dB(1,1) + f_p0c * dkm(2,1) / ele%value(l$)
      field%dB(1,2) = field%dB(1,2) + f_p0c * dkm(2,2) / ele%value(l$)
      field%dB(2,1) = field%dB(2,1) - f_p0c * dkm(1,1) / ele%value(l$)
      field%dB(2,2) = field%dB(2,2) - f_p0c * dkm(1,2) / ele%value(l$)
    endif
  enddo

endif

!----------------------
! convert fields to lab coords

if (.not. local_ref_frame) then

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

  if (ele%value(x_pitch_tot$) /= 0) then
    field%B(1) = field%B(1) + ele%value(x_pitch_tot$) * field%B(3)
    field%E(1) = field%E(1) + ele%value(x_pitch_tot$) * field%E(3)
  endif

  if (ele%value(y_pitch_tot$) /= 0) then
    field%B(2) = field%B(2) + ele%value(y_pitch_tot$) * field%B(3)
    field%E(2) = field%E(2) + ele%value(y_pitch_tot$) * field%E(3)
  endif

!-------------------------------
! Add kicks. Since the kick field is not rotated by a tilt then we have to unrotate if in the local_ref_frame

if (has_kick_attributes(ele%key) .and. (ele%value(hkick$) /= 0 .or. ele%value(vkick$) /= 0)) then
  select case (ele%key)
  case (kicker$, hkicker$, vkicker$, elseparator$)
  case default
    if (.not. local_ref_frame .or. ele%value(tilt$) == 0) then
      field%b(1) = field%b(1) + ele%value(Vkick$) * f_p0c 
      field%b(2) = field%b(2) - ele%value(Hkick$) * f_p0c 
    else
      tilt = -ele%value(tilt$)
      field%b(1) = field%b(1) + (ele%value(Vkick$) * cos(tilt) + ele%value(hkick$) * sin(tilt)) * f_p0c 
      field%b(2) = field%b(2) - (ele%value(Hkick$) * cos(tilt) - ele%value(vkick$) * sin(tilt)) * f_p0c 
    endif
  end select
endif

endif

!----------------------------------------------------------------------------
! Map field calc 
!----------------------------------------------------------------------------
case(map$)


!------------------------------------------

select case (ele%key)

!------------------------------------------
! RFcavity and Lcavity

case(rfcavity$, lcavity$)
  if (.not. associated(ele%em_field)) then
      print *, 'ERROR IN EM_FIELD_CALC: No accociated em_field for field calc = Map'
      call err_exit
  endif

    E_rho = 0; E_phi = 0; E_z = 0
    B_rho = 0; B_phi = 0; B_z = 0

    radius = sqrt(x**2 + y**2)
    phi = atan2(y, x)

    s_pos = s_rel + ele%value(ds_slave_offset$)

    ! 
    if (ele%key == rfcavity$) then
      t_ref = (ele%value(phi0$) + ele%value(dphi0$) + ele%value(phi0_err$)) 
    else
      t_ref = -(ele%value(phi0$) + ele%value(dphi0$))
    endif
    t_ref = t_ref / ele%em_field%mode(1)%freq

    do i = 1, size(ele%em_field%mode)
      mode => ele%em_field%mode(i)
      m = mode%m

      k_t = twopi * mode%freq / c_light

      Er = 0; Ep = 0; Ez = 0
      Br = 0; Bp = 0; Bz = 0

      do n = 1, size(mode%map%term)

        term => mode%map%term(n)
        k_zn = twopi * (n - 1) / (size(mode%map%term) * mode%map%dz)
        if (2 * n > size(mode%map%term)) k_zn = k_zn - twopi / mode%map%dz

        expi = cmplx(cos(k_zn * s_pos), sin(k_zn * s_pos))

        kappa2_n = k_zn**2 - k_t**2
        kappa_n = sqrt(abs(kappa2_n))
        kap_rho = kappa_n * radius
        if (kappa2_n < 0) then
          kappa_n = -i_imaginary * kappa_n
          kap_rho = -kap_rho
        endif

        if (m == 0) then
          Im_0    = I_bessel(0, kap_rho)
          Im_plus = I_bessel(1, kap_rho) / kappa_n

          Er = Er - term%e_coef * Im_plus * expi * I_imaginary * k_zn
          Ep = Ep + term%b_coef * Im_plus * expi
          Ez = Ez + term%e_coef * Im_0    * expi

          Br = Br - term%b_coef * Im_plus * expi * k_zn
          Bp = Bp - term%e_coef * Im_plus * expi * k_t**2 * I_imaginary
          Bz = Bz - term%b_coef * Im_0    * expi * I_imaginary

        else
          cm = expi * cos(m * phi - mode%phi_0)
          sm = expi * sin(m * phi - mode%phi_0)
          Im_plus  = I_bessel(m+1, kap_rho) / kappa_n**(m+1)
          Im_minus = I_bessel(m-1, kap_rho) / kappa_n**(m-1)

          ! Reason for computing Im_0_R like this is to avoid divide by zero when radius = 0.
          Im_0_R  = (Im_minus - Im_plus * kappa_n**2) / (2 * m) ! = Im_0 / radius
          Im_0    = radius * Im_0_R       

          Er = Er - i_imaginary * (k_zn * term%e_coef * Im_plus + term%b_coef * Im_0_R) * cm
          Ep = Ep - i_imaginary * (k_zn * term%e_coef * Im_plus + term%b_coef * (Im_0_R - Im_minus / m)) * sm
          Ez = Ez + term%e_coef * Im_0 * cm
   
          Br = Br + i_imaginary * sm * (term%e_coef * (m * Im_0_R + k_zn**2 * Im_plus) + &
                                        term%b_coef * k_zn * (m * Im_0_R - Im_minus / m))
          Bp = Bp + i_imaginary * cm * (term%e_coef * (Im_minus - (k_zn**2 + k_t**2) * Im_plus) / 2 - &
                                        term%b_coef * k_zn * Im_0_R)
          Bz = Bz +               sm * (-term%e_coef * k_zn * Im_0 + term%b_coef * kappa2_n * Im_0 / m)


       endif
        
      enddo

      expt = mode%field_scale * exp(-I_imaginary * twopi * (mode%freq * (t_rel + t_ref) + mode%theta_t0))
      E_rho = E_rho + Er * expt
      E_phi = E_phi + Ep * expt
      E_z   = E_z   + Ez * expt

      expt = expt / (twopi * mode%freq)
      B_rho = B_rho + Br * expt
      B_phi = B_phi + Bp * expt
      B_z   = B_z   + Bz * expt

    enddo

    field%E = [cos(phi) * real(E_rho) - sin(phi) * real(E_phi), sin(phi) * real(E_rho) + cos(phi) * real(E_phi), real(E_z)]
    field%B = [cos(phi) * real(B_rho) - sin(phi) * real(B_phi), sin(phi) * real(B_rho) + cos(phi) * real(B_phi), real(B_z)]

!------------------------------------------
! Error

case default
  print *, 'ERROR IN EM_FIELD_CALC: ELEMENT NOT YET CODED FOR MAP METHOD: ', key_name(ele%key)
  print *, '      FOR: ', ele%name
  call err_exit
end select

!----------------------------------------------------------------------------
! Grid field calc 
!----------------------------------------------------------------------------
case(grid$)

  !-----------------------------------------
  ! 

  if (.not. associated(ele%em_field)) then
    print *, 'ERROR IN EM_FIELD_CALC: No accociated em_field for field calc = Grid'
    call err_exit
  endif
  
  !radial coordinate
  r = sqrt(x**2 + y**2)

  !???
  s_pos = s_rel + ele%value(ds_slave_offset$)

  !reference time for oscillating elements
  select case (ele%key)
    case(rfcavity$) 
      t_ref = (ele%value(phi0$) + ele%value(dphi0$) + ele%value(phi0_err$)) / ele%em_field%mode(1)%freq
  
    case(lcavity$)
      t_ref = -(ele%value(phi0$) + ele%value(dphi0$))  / ele%em_field%mode(1)%freq

    case default
      t_ref = 0
  end select

  !Loop over modes
  do i = 1, size(ele%em_field%mode)
    mode => ele%em_field%mode(i)
    m = mode%m
    
    !DC modes should have mode%freq = 0
    expt = mode%field_scale * exp(-I_imaginary * twopi * &
					(mode%freq * (t_rel + t_ref) + mode%theta_t0))

    !Check for grid
    if (.not. associated(mode%grid)) then
        call out_io (s_fatal$, r_name, 'ERROR IN EM_FIELD_CALC: Missing grid for ele: ' // ele%name)
      call err_exit
    endif

    !calculate field based on grid type
    select case(mode%grid%type)
    
      case(rotationally_symmetric_rz$)
      
      !Format should be: pt (ir, iz) = ( Er, 0, Ez, 0, Bphi, 0 ) 
        
        !Interpolate 2D (r, z) grid
        !local_field is a em_field_pt_struct, which has complex E and B
        call em_grid_linear_interpolate(mode%grid, local_field, r, s_pos)
        !Transverse field is zero on axis. Otherwise:
	    if (r /= 0) then
		  !Get non-rotated field
		  E_rho = real(expt*local_field%E(1))
 		  E_phi = real(expt*local_field%E(2))
 		  B_rho = real(expt*local_field%B(1)) 
	 	  B_phi = real(expt*local_field%B(2))

		 !rotate field and output Ex, Ey, Bx, By
		 field%e(1) = field%e(1) +  (x*E_rho - y*E_phi)/r
	 	 field%e(2) = field%e(2) +  (y*E_rho + x*E_phi)/r
		 field%b(1) = field%b(1) +  (x*B_rho - y*B_phi)/r
		 field%b(2) = field%b(2) +  (y*B_rho + x*B_phi)/r
	   endif
	
	   !Ez, Bz 
	     field%e(3) = field%e(3) + real(expt*local_field%E(3))
	     field%b(3) = field%b(3) + real(expt*local_field%B(3)) 
  
      case default
        call out_io (s_fatal$, r_name, 'UNKOWN GRID TYPE FOR ELEMENT: ' // ele%name)
        call err_exit
    end select
  enddo


!----------------------------------------------------------------------------
! Unknown field calc
!----------------------------------------------------------------------------
case default
  call out_io (s_fatal$, r_name, 'BAD FIELD_CALC METHOD FOR ELEMENT: ' // ele%name)
  call err_exit
end select




end subroutine em_field_calc 

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine em_grid_linear_interpolate (grid, field, x1, x2, x3 )
!
! Subroutine to interpolate the E and B fields on a rectilinear grid
!
! Note: No error checking is done for providing x2 or x3 for 2D and 3D calculations!
!
! Modules needed:
!   use bmad_struct
!
! Input:
!   grid    --em_field_grid_struct
!   x1      -- real(rp) : dimension 1 interpolation point
!   x2      -- real(rp), optional : dimension 2 interpolation point
!   x3      -- real(rp), optional : dimension 3 interpolation point
!
! Output:
!   field  -- em_field_pt_struct: Interpolated field (complex)
!-

subroutine em_grid_linear_interpolate (grid, field, x1, x2, x3)

type (em_field_grid_struct), intent(in) :: grid
type (em_field_point_struct), intent(out)     :: field
real(rp) :: x1
real(rp), optional :: x2, x3
real(rp) rel_x1, rel_x2, rel_x3, approx_i1, approx_i2, approx_i3
integer i1, i2, i3


!Pick appropriate dimension 
select case(em_grid_dimension(grid%type))
	
  case (1)

	approx_i1 = (x1 - grid%r0(1)) / grid%dr(1); i1 = floor(approx_i1) ! index of lower x1 data point
    rel_x1 = approx_i1 - i1 !Relative distance from lower x1 grid point

  	!Check for bad indices
	if (    (i1 < lbound(grid%pt, 1)) .or. &
			(i1 > (ubound(grid%pt, 1) - 1)) ) then
	  	print *, 'Warning in  1D GRID interpolation: indicies out of bounds:'
	  	print *, '            i1 =', i1
	  	print *, 'Setting field to zero'
  		field%E = 0
  		field%B = 0
  		return
	end if			
  


    ! Do linear interpolation
    field%E(:) = (1-rel_x1) * grid%pt(i1,   1, 1)%E(:) &
               + (rel_x1)   * grid%pt(i1+1, 1, 1)%E(:) 
    field%B(:) = (1-rel_x1) * grid%pt(i1,   1, 1)%B(:) &
               + (rel_x1)   * grid%pt(i1+1, 1, 1)%B(:) 

  case (2)
  

	approx_i1 = (x1 - grid%r0(1)) / grid%dr(1); i1 = floor(approx_i1) ! index of lower x1 data point
	approx_i2 = (x2 - grid%r0(2)) / grid%dr(2); i2 = floor(approx_i2) ! index of lower x2 data point

    rel_x1 = approx_i1 - i1 !Relative distance from lower x1 grid point
    rel_x2 = approx_i2 - i2 !Relative distance from lower x2 grid point
  
  	!Check for bad indices
	if (    (i1 < lbound(grid%pt, 1)) .or. &
			(i1 > (ubound(grid%pt, 1) - 1)) .or. &
			(i2 < lbound(grid%pt, 2)) .or. &
			(i2 > (ubound(grid%pt, 2) - 1)) ) then
	  	print *, 'Warning in 2D GRID interpolation: indicies out of bounds:'
	  	print *, '            i1 =', i1
	  	print *, '            i2 =', i2
	  	print *, 'Setting field to zero'
  		field%E = 0
  		field%B = 0
  		return
	end if			
  

    
    ! Do bilinear interpolation
    field%E(:) = (1-rel_x1)*(1-rel_x2) * grid%pt(i1, i2,    1)%E(:) &
               + (1-rel_x1)*(rel_x2)   * grid%pt(i1, i2+1,  1)%E(:) &
               + (rel_x1)*(1-rel_x2)   * grid%pt(i1+1, i2,  1)%E(:) &
               + (rel_x1)*(rel_x2)     * grid%pt(i1+1, i2+1,1)%E(:) 
    field%B(:) = (1-rel_x1)*(1-rel_x2) * grid%pt(i1, i2,    1)%B(:) &
               + (1-rel_x1)*(rel_x2)   * grid%pt(i1, i2+1,  1)%B(:) &
               + (rel_x1)*(1-rel_x2)   * grid%pt(i1+1, i2,  1)%B(:) &
               + (rel_x1)*(rel_x2)     * grid%pt(i1+1, i2+1,1)%B(:)	
						
  case (3)

	approx_i1 = (x1 - grid%r0(1)) / grid%dr(1); i1 = floor(approx_i1) ! index of lower x1 data point
	approx_i2 = (x2 - grid%r0(2)) / grid%dr(2); i2 = floor(approx_i2) ! index of lower x2 data point
	approx_i3 = (x3 - grid%r0(3)) / grid%dr(3); i3 = floor(approx_i3) ! index of lower x3 data point

    rel_x1 = approx_i1 - i1 !Relative distance from lower x1 grid point
    rel_x2 = approx_i2 - i2 !Relative distance from lower x2 grid point
    rel_x3 = approx_i3 - i3 !Relative distance from lower x3 grid point

  	!Check for bad indices
	if (    (i1 < lbound(grid%pt, 1)) .or. &
			(i1 > (ubound(grid%pt, 1) - 1)) .or. &
			(i2 < lbound(grid%pt, 2)) .or. &
			(i2 > (ubound(grid%pt, 2) - 1)) .or. &
			(i3 < lbound(grid%pt, 3)) .or. &
			(i3 > (ubound(grid%pt, 3) - 1)) ) then
	  	print *, 'Warning in 3D GRID interpolation: indicies out of bounds:'
	  	print *, '            i1 =', i1
	  	print *, '            i2 =', i2
	  	print *, '            i3 =', i3
	  	print *, 'Setting field to zero'
  		field%E = 0
  		field%B = 0
  		return
	end if			

    
    ! Do trilinear interpolation
    field%E(:) = (1-rel_x1)*(1-rel_x2)*(1-rel_x3) * grid%pt(i1, i2,    i3  )%E(:) &
               + (1-rel_x1)*(rel_x2)  *(1-rel_x3) * grid%pt(i1, i2+1,  i3  )%E(:) &
               + (rel_x1)  *(1-rel_x2)*(1-rel_x3) * grid%pt(i1+1, i2,  i3  )%E(:) &
               + (rel_x1)  *(rel_x2)  *(1-rel_x3) * grid%pt(i1+1, i2+1,i3  )%E(:) &
               + (1-rel_x1)*(1-rel_x2)*(rel_x3)   * grid%pt(i1, i2,    i3+1)%E(:) &
               + (1-rel_x1)*(rel_x2)  *(rel_x3)   * grid%pt(i1, i2+1,  i3+1)%E(:) &
               + (rel_x1)  *(1-rel_x2)*(rel_x3)   * grid%pt(i1+1, i2,  i3+1)%E(:) &
               + (rel_x1)  *(rel_x2)  *(rel_x3)   * grid%pt(i1+1, i2+1,i3+1)%E(:)               
               
    ! Do bilinear interpolation
    field%B(:) = (1-rel_x1)*(1-rel_x2)*(1-rel_x3) * grid%pt(i1, i2,    i3  )%B(:) &
               + (1-rel_x1)*(rel_x2)  *(1-rel_x3) * grid%pt(i1, i2+1,  i3  )%B(:) &
               + (rel_x1)  *(1-rel_x2)*(1-rel_x3) * grid%pt(i1+1, i2,  i3  )%B(:) &
               + (rel_x1)  *(rel_x2)  *(1-rel_x3) * grid%pt(i1+1, i2+1,i3  )%B(:) &
               + (1-rel_x1)*(1-rel_x2)*(rel_x3)   * grid%pt(i1, i2,    i3+1)%B(:) &
               + (1-rel_x1)*(rel_x2)  *(rel_x3)   * grid%pt(i1, i2+1,  i3+1)%B(:) &
               + (rel_x1)  *(1-rel_x2)*(rel_x3)   * grid%pt(i1+1, i2,  i3+1)%B(:) &
               + (rel_x1)  *(rel_x2)  *(rel_x3)   * grid%pt(i1+1, i2+1,i3+1)%B(:) 


  case default
    print *, 'ERROR IN EM_GRID_INTERPOLATE, BAD DIMENSION:', em_grid_dimension(grid%type)
    call err_exit   
end select

end subroutine em_grid_linear_interpolate	

end module
