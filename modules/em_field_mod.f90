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

recursive subroutine em_field_calc (ele, param, s_rel, t_rel, orbit, local_ref_frame, field, calc_dfield)

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
type (lat_param_struct) param
type (coord_struct) :: orbit, local_orb
type (wig_term_struct), pointer :: wig
type (em_field_struct) :: field, field2
type (em_field_point_struct) :: local_field
type (em_field_mode_struct), pointer :: mode
type (em_field_map_term_struct), pointer :: term

real(rp) :: x, y, s, t, xx, yy, t_rel, s_rel, z,   f, dk(3,3), charge, f_p0c
real(rp) :: c_x, s_x, c_y, s_y, c_z, s_z, coef, fd(3)
real(rp) :: cos_ang, sin_ang, sgn_x, dc_x, dc_y, kx, ky, dkm(2,2)
real(rp) phase, gradient, theta, r, E_r, E_s, k_wave
real(rp) k_t, k_zn, kappa2_n, kap_rho
real(rp) radius, phi, t_ref, tilt, omega

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
! convert to local coords

local_orb = orbit
if (.not. local_ref_frame) then
  call offset_particle (ele, param, local_orb, set$, &
          set_canonical = .false., set_multipoles = .false., set_hvkicks = .false.)
endif

!----------------------------------------------------------------------------
! super_slave and multipass_slave elements have their field info stored in the associated lord elements.

if (ele%field_calc == refer_to_lords$) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    s = s_rel + (ele%s - ele%value(l$)) - (lord%s - lord%value(l$))
    t = t_rel + (lord%ref_time - lord%value(delta_ref_time$)) - (ele%ref_time - ele%value(delta_ref_time$))
    call em_field_calc (lord, param, s, t, local_orb, .false., field2, calc_dfield)
    field%E = field%E + field2%E
    field%B = field%B + field2%B
    if (df_calc) then
      field%dE = field%dE + field2%dE
      field%dB = field%dB + field2%dB
    endif
  enddo
  call convert_fields_to_lab_coords
  return
endif

!----------------------------------------------------------------------------
! Custom field calc 

if (ele%field_calc == custom$) then 
  call em_field_custom (ele, param, s_rel, t_rel, orbit, local_ref_frame, field, calc_dfield)
  return
end if

!----------------------------------------------------------------------------
!Set up common variables for all (non-custom) methods

charge = charge_of(param%particle)
sign_charge = sign(1.0_rp, charge)

x = local_orb%vec(1)
y = local_orb%vec(3)

f_p0c = sign_charge * ele%value(p0c$) / c_light

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! field_calc methods

select case (ele%field_calc)
  
  !----------------------------------------------------------------------------
  ! Bmad_standard field calc 

  case (bmad_standard$)

  select case (ele%key)

  !------------------------------------------
  ! RFcavity and Lcavity
  ! Note: This does not include the edge focusing which is calcuated in the 
  !  apply_element_edge_kick routine.

  case(rfcavity$, lcavity$)

    phase = twopi * (ele%value(theta_t0$) + ele%value(phi0$) + ele%value(dphi0$) + ele%value(phi0_err$))
    if (ele%key == rfcavity$) phase = phase + pi/2

    gradient = (ele%value(gradient$) + ele%value(gradient_err$)) * ele%value(field_scale$)
    if (.not. ele%is_on) gradient = 0
    gradient = gradient + gradient_shift_sr_wake(ele, param)

    ! Use pillbox formulas for standing wave TM_011 mode with infinite wall radius.
    ! See S.Y. Lee, "Accelerator Physics"
    !   E_s   = 2 * gradient * cos(k s) * cos(omega t + phase)
    !   E_r   = gradient * k * r * sin(k s) * cos(omega t + phase)
    !   B_phi = -gradient * k * r * cos(k s) * sin(omega t + phase) / c_light
    ! where
    !   k = pi / L
    !   omega = c * k

    theta = atan2(y, x)   
    r = sqrt(x**2 + y**2)
    k_wave = pi / ele%value(l$)
    omega = c_light * k_wave

    E_r = gradient * r * k_wave * sin(k_wave*s_rel) * cos(omega * t_rel + phase)

    field%E(1) = E_r * cos (theta)                                  
    field%E(2) = E_r * sin (theta)
    field%E(3) = 2 * gradient * cos(k_wave * s_rel) * cos(omega * t_rel + phase)
    
    B_phi = -gradient * r * k_wave * cos(k_wave*s_rel) * sin(omega * t_rel + phase) / c_light 
    field%B(1) = -B_phi * sin(theta)
    field%B(2) =  B_phi * cos(theta)

    if (df_calc) then
      call out_io (s_fatal$, r_name, 'dFIELD NOT YET IMPLEMENTED FOR LCAVITY!')
      if (bmad_status%exit_on_error) call err_exit
    endif

  !------------------------------------------
  ! Wiggler

  case(wiggler$)

    do i = 1, size(ele%wig%term)
      wig => ele%wig%term(i)

      if (wig%type == hyper_y$) then
        c_x = cos(wig%kx * x)
        s_x = sin(wig%kx * x)
        sgn_x = 1
        dc_x = -1
      else
        c_x = cosh(wig%kx * x)
        s_x = sinh(wig%kx * x)
        sgn_x = -1
        dc_x = 1
      endif

      if (wig%type == hyper_y$ .or. wig%type == hyper_xy$) then
        c_y = cosh (wig%ky * y)
        s_y = sinh (wig%ky * y)
        dc_y = 1
      else
        c_y = cos (wig%ky * y)
        s_y = sin (wig%ky * y)
        dc_y = -1
      endif

      c_z = cos (wig%kz * s_rel + wig%phi_z)
      s_z = sin (wig%kz * s_rel + wig%phi_z)

      coef = sign_charge * wig%coef * ele%value(polarity$)

      field%B(1) = field%B(1) - coef  * (wig%kx / wig%ky) * s_x * s_y * c_z * sgn_x
      field%B(2) = field%B(2) + coef  *                 c_x * c_y * c_z
      field%B(3) = field%B(3) - coef  * (wig%kz / wig%ky) * c_x * s_y * s_z

      if (df_calc) then
        f = coef * wig%kx
        field%dB(1,1) = field%dB(1,1) - f  * (wig%kx / wig%ky) * c_x * s_y * c_z * sgn_x
        field%dB(2,1) = field%dB(2,1) + f  *                 s_x * c_y * c_z * dc_x
        field%dB(3,1) = field%dB(3,1) - f  * (wig%kz / wig%ky) * s_x * s_y * s_z * dc_x
        f = coef * wig%ky
        field%dB(1,2) = field%dB(1,2) - f  * (wig%kx / wig%ky) * s_x * c_y * c_z * sgn_x
        field%dB(2,2) = field%dB(2,2) + f  *                 c_x * s_y * c_z * dc_y
        field%dB(3,2) = field%dB(3,2) - f  * (wig%kz / wig%ky) * c_x * c_y * s_z 
        f = coef * wig%kz
        field%dB(1,3) = field%dB(1,3) + f  * (wig%kx / wig%ky) * s_x * s_y * s_z * sgn_x
        field%dB(2,3) = field%dB(2,3) - f  *                 c_x * c_y * s_z * dc_y
        field%dB(3,3) = field%dB(3,3) - f  * (wig%kz / wig%ky) * c_x * s_y * c_z 
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
    field%b(2) = (x**2 - y**2) / 2 * ele%value(k2$) * f_p0c 

    if (df_calc) then
      field%dB(1,1) =  y * ele%value(k2$) * f_p0c
      field%dB(1,2) =  x * ele%value(k2$) * f_p0c
      field%dB(2,1) = -x * ele%value(k2$) * f_p0c
      field%dB(2,2) = -y * ele%value(k2$) * f_p0c
    endif

  !------------------------------------------
  ! Sextupole 

  case (octupole$)

    field%b(1) = -(y**3 - 3*y*x**2) / 6 * ele%value(k2$) * f_p0c 
    field%b(2) =  (x**3 - 3*x*y**2) / 6 * ele%value(k2$) * f_p0c 

    if (df_calc) then
      field%dB(1,1) =  y * ele%value(k2$) * f_p0c
      field%dB(1,2) =  x * ele%value(k2$) * f_p0c
      field%dB(2,1) = -x * ele%value(k2$) * f_p0c
      field%dB(2,2) = -y * ele%value(k2$) * f_p0c
    endif

  !------------------------------------------
  ! Sol_quad

  case (sol_quad$)

    if (bmad_status%type_out) &
      call out_io (s_fatal$, r_name, 'BMAD_STANDARD SOL_QUAD FIELD NOT IMPLEMENTED SINCE EDGE FIELD IS NOT DEFINED.')
    if (bmad_status%exit_on_error) call err_exit
    return

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

    if (bmad_status%type_out) &
      call out_io (s_fatal$, r_name, 'BMAD_STANDARD SOLENOID FIELD NOT IMPLEMENTED SINCE EDGE FIELD IS NOT DEFINED.')
    if (bmad_status%exit_on_error) call err_exit
    return

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
    call out_io (s_fatal$, r_name, 'ELEMENT NOT YET CODED: ' // key_name(ele%key), 'FOR: ' // ele%name)
    if (bmad_status%exit_on_error) call err_exit
    return
  end select

  !---------------------------------------------------------------------
  ! Add multipoles

  if (associated(ele%a_pole)) then
    if (ele%value(l$) == 0) then
      call out_io (s_fatal$, r_name, 'dField NOT YET IMPLEMENTED FOR MULTIPOLES!', 'FOR: ' // ele%name)
      if (bmad_status%exit_on_error) call err_exit
      return
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

  !-------------------------------
  ! Add kicks. Since the kick field is not rotated by a tilt then we have to unrotate if in the local_ref_frame

  if (has_kick_attributes(ele%key) .and. (ele%value(hkick$) /= 0 .or. ele%value(vkick$) /= 0)) then
    select case (ele%key)
    ! Handled above
    case (kicker$, hkicker$, vkicker$, elseparator$)  
    ! Everything else
    case default
      if (.not. local_ref_frame .or. ele%value(tilt$) == 0) then
        field%b(1) = field%b(1) + ele%value(Vkick$) * f_p0c 
        field%b(2) = field%b(2) - ele%value(Hkick$) * f_p0c 
      else
        ! Rotate from lab to local
        tilt = ele%value(tilt$)
        field%b(1) = field%b(1) + (ele%value(Vkick$) * cos(tilt) - ele%value(hkick$) * sin(tilt)) * f_p0c 
        field%b(2) = field%b(2) - (ele%value(Hkick$) * cos(tilt) - ele%value(vkick$) * sin(tilt)) * f_p0c 
      endif
    end select
  endif

!----------------------------------------------------------------------------
! Map field calc 

case(map$)

!------------------------------------------

  select case (ele%key)

  !------------------------------------------
  ! RFcavity and Lcavity

  case(rfcavity$, lcavity$)
    if (.not. associated(ele%em_field)) then
      call out_io (s_fatal$, r_name, 'No accociated em_field for field calc = Map', 'FOR: ' // ele%name) 
      if (bmad_status%exit_on_error) call err_exit
      return  
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

        expi = cmplx(cos(k_zn * s_rel), sin(k_zn * s_rel))

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
    call out_io (s_fatal$, r_name, 'ELEMENT NOT YET CODED FOR MAP METHOD: ' // key_name(ele%key), &
                                  'FOR: ' // ele%name)
    if (bmad_status%exit_on_error) call err_exit
  end select

!----------------------------------------------------------------------------
! Grid field calc 

case(grid$)

  !-----------------------------------------
  ! 

  if (.not. associated(ele%em_field)) then
    call out_io (s_fatal$, r_name, 'No accociated em_field for field calc = Grid', 'FOR: ' // ele%name)
    if (bmad_status%exit_on_error) call err_exit
    return
  endif
  
  ! radial coordinate
  r = sqrt(x**2 + y**2)

  ! reference time for oscillating elements
  select case (ele%key)
    case(rfcavity$) 
      t_ref = (ele%value(phi0$) + ele%value(dphi0$) + ele%value(phi0_err$)) / ele%em_field%mode(1)%freq
  
    case(lcavity$)
      t_ref = -(ele%value(phi0$) + ele%value(dphi0$))  / ele%em_field%mode(1)%freq

    case default
      t_ref = 0
  end select

  ! Loop over modes
  do i = 1, size(ele%em_field%mode)
    mode => ele%em_field%mode(i)
    m = mode%m
    
    ! DC modes should have mode%freq = 0
    expt = mode%field_scale * exp(-I_imaginary * twopi * &
          (mode%freq * (t_rel + t_ref) + mode%theta_t0))

    ! Check for grid
    if (.not. associated(mode%grid)) then
      call out_io (s_fatal$, r_name, 'MISSING GRID FOR ELE: ' // ele%name)
      if (bmad_status%exit_on_error) call err_exit
      return
    endif

    ! calculate field based on grid type
    select case(mode%grid%type)
    
    case(rotationally_symmetric_rz$)
      
      ! Format should be: pt (ir, iz) = ( Er, 0, Ez, 0, Bphi, 0 ) 
        
      ! Interpolate 2D (r, z) grid
      ! local_field is a em_field_pt_struct, which has complex E and B
      call em_grid_linear_interpolate(ele, mode%grid, local_field, r, s_rel)
      ! Transverse field is zero on axis. Otherwise:

      if (r /= 0) then
        ! Get non-rotated field
        E_rho = real(expt*local_field%E(1))
        E_phi = real(expt*local_field%E(2))
        B_rho = real(expt*local_field%B(1)) 
        B_phi = real(expt*local_field%B(2))

        ! rotate field and output Ex, Ey, Bx, By
        field%e(1) = field%e(1) +  (x*E_rho - y*E_phi)/r
        field%e(2) = field%e(2) +  (y*E_rho + x*E_phi)/r
        field%b(1) = field%b(1) +  (x*B_rho - y*B_phi)/r
        field%b(2) = field%b(2) +  (y*B_rho + x*B_phi)/r
      endif
  
      ! Ez, Bz 
      field%e(3) = field%e(3) + real(expt*local_field%E(3))
      field%b(3) = field%b(3) + real(expt*local_field%B(3)) 
  
    case default
      call out_io (s_fatal$, r_name, 'UNKNOWN GRID TYPE: \i0\ ', &
                                     'FOR ELEMENT: ' // ele%name, i_array = [mode%grid%type])
      if (bmad_status%exit_on_error) call err_exit
      return
    end select
  enddo

!----------------------------------------------------------------------------
! Unknown field calc

case default
  call out_io (s_fatal$, r_name, 'BAD FIELD_CALC METHOD FOR ELEMENT: ' // ele%name)
  if (bmad_status%exit_on_error) call err_exit
  return
end select

call convert_fields_to_lab_coords

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! convert fields to lab coords

contains

subroutine convert_fields_to_lab_coords

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

end subroutine convert_fields_to_lab_coords

end subroutine em_field_calc 

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine em_grid_linear_interpolate (ele, grid, field, x1, x2, x3 )
!
! Subroutine to interpolate the E and B fields on a rectilinear grid
!
! Note: No error checking is done for providing x2 or x3 for 2D and 3D calculations!
!
! Modules needed:
!   use bmad_struct
!
! Input:
!   ele     -- ele_struct: Element containing the grid
!   grid    -- em_field_grid_struct: Grid to interpolate
!   x1      -- real(rp) : dimension 1 interpolation point
!   x2      -- real(rp), optional : dimension 2 interpolation point
!   x3      -- real(rp), optional : dimension 3 interpolation point
!
! Output:
!   field  -- em_field_pt_struct: Interpolated field (complex)
!-

subroutine em_grid_linear_interpolate (ele, grid, field, x1, x2, x3)

type (ele_struct) ele
type (em_field_grid_struct) :: grid
type (em_field_point_struct), intent(out) :: field
real(rp) :: x1
real(rp), optional :: x2, x3
real(rp) rel_x1, rel_x2, rel_x3
integer i1, i2, i3, grid_dim
logical err1, err2, err3

character(32), parameter :: r_name = 'em_grid_linear_interpolate'

! Pick appropriate dimension 

grid_dim = em_grid_dimension(grid%type)
select case(grid_dim)

case (1)

  call get_this_index(x1, 1, i1, rel_x1, err1)
  if (err1) return

  field%E(:) = (1-rel_x1) * grid%pt(i1, 1, 1)%E(:) + (rel_x1) * grid%pt(i1+1, 1, 1)%E(:) 
  field%B(:) = (1-rel_x1) * grid%pt(i1, 1, 1)%B(:) + (rel_x1) * grid%pt(i1+1, 1, 1)%B(:) 

case (2)

  call get_this_index(x1, 1, i1, rel_x1, err1)
  call get_this_index(x2, 2, i2, rel_x2, err2)
  if (err1 .or. err2) return

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

  call get_this_index(x1, 1, i1, rel_x1, err1)
  call get_this_index(x2, 2, i2, rel_x2, err2)
  call get_this_index(x3, 3, i3, rel_x3, err3)
  if (err1 .or. err2 .or. err3) return
    
  ! Do trilinear interpolation
  field%E(:) = (1-rel_x1)*(1-rel_x2)*(1-rel_x3) * grid%pt(i1, i2,    i3  )%E(:) &
             + (1-rel_x1)*(rel_x2)  *(1-rel_x3) * grid%pt(i1, i2+1,  i3  )%E(:) &
             + (rel_x1)  *(1-rel_x2)*(1-rel_x3) * grid%pt(i1+1, i2,  i3  )%E(:) &
             + (rel_x1)  *(rel_x2)  *(1-rel_x3) * grid%pt(i1+1, i2+1,i3  )%E(:) &
             + (1-rel_x1)*(1-rel_x2)*(rel_x3)   * grid%pt(i1, i2,    i3+1)%E(:) &
             + (1-rel_x1)*(rel_x2)  *(rel_x3)   * grid%pt(i1, i2+1,  i3+1)%E(:) &
             + (rel_x1)  *(1-rel_x2)*(rel_x3)   * grid%pt(i1+1, i2,  i3+1)%E(:) &
             + (rel_x1)  *(rel_x2)  *(rel_x3)   * grid%pt(i1+1, i2+1,i3+1)%E(:)               
             
  field%B(:) = (1-rel_x1)*(1-rel_x2)*(1-rel_x3) * grid%pt(i1, i2,    i3  )%B(:) &
             + (1-rel_x1)*(rel_x2)  *(1-rel_x3) * grid%pt(i1, i2+1,  i3  )%B(:) &
             + (rel_x1)  *(1-rel_x2)*(1-rel_x3) * grid%pt(i1+1, i2,  i3  )%B(:) &
             + (rel_x1)  *(rel_x2)  *(1-rel_x3) * grid%pt(i1+1, i2+1,i3  )%B(:) &
             + (1-rel_x1)*(1-rel_x2)*(rel_x3)   * grid%pt(i1, i2,    i3+1)%B(:) &
             + (1-rel_x1)*(rel_x2)  *(rel_x3)   * grid%pt(i1, i2+1,  i3+1)%B(:) &
             + (rel_x1)  *(1-rel_x2)*(rel_x3)   * grid%pt(i1+1, i2,  i3+1)%B(:) &
             + (rel_x1)  *(rel_x2)  *(rel_x3)   * grid%pt(i1+1, i2+1,i3+1)%B(:) 


case default
  call out_io (s_fatal$, r_name, 'BAD DIMENSION: \i0\ ', em_grid_dimension(grid%type))
  if (bmad_status%exit_on_error) call err_exit
  return
end select

!-------------------------------------------------------------------------------------
contains

subroutine get_this_index (x, ix_x, i0, rel_x0, err_flag)

implicit none

real(rp) x, rel_x0, x_norm
integer ix_x, i0
logical err_flag

!

x_norm = (x - grid%r0(ix_x)) / grid%dr(ix_x)
i0 = floor(x_norm)     ! index of lower 1 data point
rel_x0 = x_norm - i0   ! Relative distance from lower x1 grid point

if (i0 == ubound(grid%pt, ix_x) .and. rel_x0 < bmad_com%significant_length) then
  i0 = ubound(grid%pt, ix_x) - 1
  rel_x0 = 1
endif

if (i0 == lbound(grid%pt, ix_x) - 1 .and. abs(rel_x0 - 1) < bmad_com%significant_length) then
  i0 = lbound(grid%pt, ix_x)
  rel_x0 = 0
endif

if (i0 < lbound(grid%pt, ix_x) .or. i0 >= ubound(grid%pt, ix_x)) then
  call out_io (s_warn$, r_name, '\i0\D GRID interpolation index out of bounds: i\i0\ = \i0\ ', &
                                'For element: ' // ele%name, &
                                'Setting field to zero', i_array = [grid_dim, ix_x, i0])
  
  field%E = 0
  field%B = 0
  err_flag = .true.
  return
endif

err_flag = .false.

end subroutine get_this_index 

end subroutine em_grid_linear_interpolate

end module
