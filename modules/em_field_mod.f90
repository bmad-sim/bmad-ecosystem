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

orb2 = orb
if (local_ref_frame) call offset_particle (ele, param, orb2, unset$, set_canonical = .false., s_pos = s)

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
! Note: B field is appropriate for particle of positive charge.
! Note: Zero field will be returned if an element is turned off.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element
!   param  -- lat_param_struct: Lattice parameters.
!   s_rel  -- Real(rp): Longitudinal position relative to the start of the element.
!   t_rel  -- Real(rp): Time relative to the reference particle.
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
type (rf_field_mode_struct), pointer :: mode

real(rp) :: x, y, xx, yy, cm, sm, t_rel, s_rel, f, dk(3,3), charge, f_p0c
real(rp) :: c_x, s_x, c_y, s_y, c_z, s_z, coef, fd(3)
real(rp) :: cos_ang, sin_ang, sgn_x, dc_x, dc_y, kx, ky, dkm(2,2)
real(rp) phase, gradient, dEz_dz, theta, r, E_r
real(rp) k_t, k_zn, kappa2_n, kap_rho
real(rp) radius, phi, t_ref

complex(rp) E_rho, E_phi, E_z, Er, Ep, Ez, B_rho, B_phi, B_z, Br, Bp, Bz, expi, expt, dEp, dEr
complex(rp) Im0, Im_plus, Im_minus, Im_norm, kappa_n, Im_plus2

integer i, j, m, n, sign_charge

logical :: local_ref_frame
logical, optional :: calc_dfield

logical df_calc
character(20) :: r_name = 'em_field_calc'

! If element is turned off then return zero

field%e = 0
field%B = 0

df_calc = logic_option (.false., calc_dfield)

if (df_calc) then
  field%dB = 0
  field%dE = 0
endif

if (.not. ele%is_on) return

! custom field_calc

select case (ele%field_calc)
case (custom$) 
  call em_field_custom (ele, param, s_rel, t_rel, orbit, local_ref_frame, field, calc_dfield)
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

local_orb = orbit
if (.not. local_ref_frame) then
  call offset_particle (ele, param, local_orb, set$, &
          set_canonical = .false., set_multipoles = .false., set_hvkicks = .false.)
endif

x = local_orb%vec(1)
y = local_orb%vec(3)

f_p0c = sign_charge * ele%value(p0c$) / c_light

!------------------------------------------

select case (ele%key)

!------------------------------------------
! RFcavity and Lcavity

case(rfcavity$, lcavity$)

  if (.not. associated(ele%rf%field)) then
    print *, 'ERROR IN EM_FIELD_CALC: RF FIELD NOT DEFINED FOR: ' // ele%name
    call err_exit
  endif

  E_rho = 0; E_phi = 0; E_z = 0
  B_rho = 0; B_phi = 0; B_z = 0

  radius = sqrt(x**2 + y**2)
  phi = atan2(y, x)

  ! 
  if (ele%key == rfcavity$) then
    t_ref = (ele%value(phi0$) + ele%value(dphi0$) + ele%value(phi0_err$)) 
  else
    t_ref = -(ele%value(phi0$) + ele%value(dphi0$))
  endif
  t_ref = t_ref / ele%rf%field%mode(1)%freq

  do i = 1, size(ele%rf%field%mode)
    mode => ele%rf%field%mode(i)
    m = mode%m

    if (m /= 0) then
      print *, 'ERROR IN EM_FIELD_CALC: RF FIELD WITH M /= 0 NOT YET IMPLEMENTED FOR: ' // ele%name
      call err_exit
    endif

    if (any(mode%term%b /= 0)) then
      print *, 'ERROR IN EM_FIELD_CALC: RF FIELD WITH M == 0 NON-ACCELERATING NOT YET IMPLEMENTED FOR: ' // ele%name
      call err_exit
    endif

    k_t = twopi * mode%freq / c_light

    Er = 0; Ep = 0; Ez = 0
    Br = 0; Bp = 0; Bz = 0

    do n = 1, size(mode%term)

      k_zn = twopi * (n - 1) / (size(mode%term) * mode%dz)
      if (2 * n > size(mode%term)) k_zn = k_zn - twopi / mode%dz

      expi = cmplx(cos(k_zn * s_rel), sin(k_zn * s_rel))

      kappa2_n = k_zn**2 - k_t**2
      kappa_n = sqrt(abs(kappa2_n))
      kap_rho = kappa_n * radius
      if (kappa2_n < 0) then
        kappa_n = -i_imaginary * kappa_n
        kap_rho = -kap_rho
      endif

      if (m == 0) then
        Im0     = I_bessel(0, kap_rho)
        Im_plus = I_bessel(1, kap_rho) / kappa_n

        Er = Er - mode%term(n)%e * Im_plus * expi * I_imaginary * k_zn
        Ep = Ep + mode%term(n)%b * Im_plus * expi
        Ez = Ez + mode%term(n)%e * Im0     * expi

        Br = Br - mode%term(n)%b * Im_plus * expi * I_imaginary * k_zn
        Bp = Bp + mode%term(n)%e * Im_plus * expi * k_t**2
        Bz = Bz + mode%term(n)%b * Im0     * expi

      else
        cm = cos(m * phi - mode%phi_0)
        sm = sin(m * phi - mode%phi_0)
        Im_plus  = I_bessel(m+1, kap_rho) / kappa_n**(m+1)
        Im_minus = I_bessel(m-1, kap_rho) / kappa_n**(m-1)

        Im_plus2 = I_bessel(m+2, kap_rho) / kappa_n**(m+2)

        Im_norm  = (Im_minus - Im_plus * kappa_n**2) / (2 * m) ! = Im / radius
        Im0      = radius * Im_norm       

        dEr = -i_imaginary * (k_zn * mode%term(n)%e * Im_plus + mode%term(n)%b * Im_norm) * cm * expi
        dEp = -i_imaginary * (k_zn * mode%term(n)%e * Im_plus + mode%term(n)%b * (Im_norm - Im_minus / m)) * expi

        Er = Er + dEr
        Ep = Ep + dEp * sm
        Ez = Ez + mode%term(n)%e * Im0 * cm * expi
 
        Br = Br - m * mode%term(n)%e * Im_norm * sm * expi - i_imaginary * k_zn * dEp * sm
        Bp = Bp + i_imaginary * k_zn * dEr - &
                      mode%term(n)%e * cm * expi * (Im_minus - m * Im_norm)
        Bz = Bz - i_imaginary * sm * expi * (k_zn * mode%term(n)%e * (Im0 - Im_plus2 * kappa_n**2) / 2 + &
                    mode%term(n)%b * (Im_norm - Im_minus / m)) * expi

     endif
      
    enddo

    expt = mode%field_scale * exp(-I_imaginary * twopi * (mode%freq * (t_rel + t_ref) + mode%theta_t0))
    E_rho = E_rho + Er * expt
    E_phi = E_phi + Ep * expt
    E_z   = E_z   + Ez * expt

    expt = -I_imaginary * expt / (twopi * mode%freq)
    B_rho = B_rho + Br * expt
    B_phi = B_phi + Bp * expt
    B_z   = B_z   + Bz * expt

  enddo

  field%E = [cos(phi) * real(E_rho) - sin(phi) * real(E_phi), sin(phi) * real(E_rho) + cos(phi) * real(E_phi), real(E_z)]
  field%B = [cos(phi) * real(B_rho) - sin(phi) * real(B_phi), sin(phi) * real(B_rho) + cos(phi) * real(B_phi), real(B_z)]

  return

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
   
  if (bmad_com%sr_wakes_on) then
    if (bmad_com%grad_loss_sr_wake /= 0) then  
      ! use grad_loss_sr_wake and ignore e_loss
      gradient = gradient - bmad_com%grad_loss_sr_wake
    elseif (ele%value(e_loss$) /= 0) then
      gradient = gradient - e_loss_sr_wake(ele%value(e_loss$), param) / ele%value(l$)
    endif
  endif

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
  field%E(3) = gradient * sin (f + phase)
  
  phi = pi - theta
  field%B(1) =   B_phi * cos (phi)
  field%B(2) = - B_phi * sin (phi)

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
! Drift, et. al.

case (drift$, ecollimator$, rcollimator$, instrument$, monitor$, pipe$)

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

endif

end subroutine em_field_calc 

end module
