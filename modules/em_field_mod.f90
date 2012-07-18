!+
! Module em_field_mod
!
! Module to define the electric and magnetic fields for an elemet.
!-

module em_field_mod

use bmad_struct
use bmad_interface

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
!     %field(1:n) -- em_field_stuct: n will be at least n_pt
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
  allocate(track%field(0:n_pt))
  allocate(track%map(0:n_pt))
endif

if (ubound(track%orb, 1) < n_pt) then
  deallocate(track%orb, track%field, track%map)
  allocate(track%orb(0:n_pt),  track%field(0:n_pt), track%map(0:n_pt))
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
if (local_ref_frame) call offset_particle (ele, orb2, param%particle, unset$, &
                                    set_canonical = .false., set_s_offset = .false., ds_pos = s)

track%orb(n_pt) = orb2
track%map(n_pt)%mat6 = 0
s_sav = s

end subroutine save_a_step

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine em_field_calc (ele, param, s_rel, time, orbit, local_ref_frame, field, calc_dfield, err_flag)
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
!   time   -- Real(rp): Particle time.
!                 For absolute time tracking this is the absolute time.
!                 For relative time tracking this is relative to the reference particle entering the element.
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
!   err_flag    -- logical, optional: Set True if there is an error. False otherwise.
!-

recursive subroutine em_field_calc (ele, param, s_rel, time, orbit, local_ref_frame, field, calc_dfield, err_flag)

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
type (lat_param_struct) param
type (coord_struct) :: orbit, local_orb
type (wig_term_struct), pointer :: wig
type (em_field_struct) :: field, field2
type (em_field_grid_pt_struct) :: local_field
type (em_field_mode_struct), pointer :: mode
type (em_field_map_term_struct), pointer :: term

real(rp) :: x, y, s, t, xx, yy, time, s_rel, z,   f, dk(3,3), charge, f_p0c
real(rp) :: c_x, s_x, c_y, s_y, c_z, s_z, coef, fd(3), s0
real(rp) :: cos_ang, sin_ang, sgn_x, dc_x, dc_y, kx, ky, dkm(2,2)
real(rp) phase, gradient, r, E_r_coef, E_s, k_wave, s_eff, t_eff
real(rp) k_t, k_zn, kappa2_n, kap_rho, s_hard_offset, beta_start
real(rp) radius, phi, t_ref, tilt, omega, freq, B_phi_coef
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)

complex(rp) E_rho, E_phi, E_z, BEr, Er, Ep, Ez, B_rho, B_phi, B_z, Br, Bp, Bz, expi, expt, dEp, dEr
complex(rp) Im_0, Im_plus, Im_minus, Im_0_R, kappa_n, Im_plus2, cm, sm

integer i, j, m, n, sign_charge

logical :: local_ref_frame, local_ref, has_nonzero_pole
logical, optional :: calc_dfield, err_flag
logical df_calc, err

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

if (present(err_flag)) err_flag = .false.
if (.not. ele%is_on) return

!----------------------------------------------------------------------------
! convert to local coords

local_orb = orbit
if (.not. local_ref_frame) then
  call offset_particle (ele, local_orb, param%particle, set$, &
          set_canonical = .false., set_multipoles = .false., set_hvkicks = .false.)
endif

!----------------------------------------------------------------------------
! super_slave, slice_slave, and  multipass_slave elements have their 
! field info stored in the associated lord elements.
! Note: The lord of an em_field element has independent misalignments.

if (ele%field_calc == refer_to_lords$) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (ele%slave_status /= slice_slave$ .and. lord%lord_status /= super_lord$) cycle

    local_ref = .true.
    if (ele%key == em_field$) local_ref = .false.

    s = s_rel + (ele%s - ele%value(l$)) - (lord%s - lord%value(l$))
    t = time 
    if (.not. absolute_time_tracking(ele)) then
      t = t + ele%value(ref_time_start$) - lord%value(ref_time_start$) 
    endif
    call em_field_calc (lord, param, s, t, local_orb, local_ref, field2, calc_dfield)

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
  call em_field_custom (ele, param, s_rel, time, orbit, local_ref_frame, field, calc_dfield)
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
  ! Use half-wave pillbox formulas for standing wave TM_011 mode with infinite wall radius.
  ! See S.Y. Lee, "Accelerator Physics"
  !   E_s   = 2 * gradient * cos(k s) * cos(omega t + phase)
  !   E_r   = gradient * k * r * sin(k s) * cos(omega t + phase)
  !   B_phi = -gradient * k * r * cos(k s) * sin(omega t + phase) / c_light
  ! where
  !   k = pi / L
  !   omega = c * k
  ! Field extends to +/- c_light * freq / 2 from centerline of element.
  ! Note: There is a discontinuity in the field at the edge. Edge focusing due to this 
  !  discontinuity can be handled in the apply_hard_edge_kick routine.

  case(rfcavity$, lcavity$)

    if (ele%value(rf_frequency$) == 0) return

    phase = twopi * (ele%value(phi0$) + ele%value(dphi0$) + ele%value(phi0_err$) + ele%value(dphi0_ref$))
    if (ele%key == rfcavity$) phase = pi/2 - phase

    if (ele%key == lcavity$) then
      gradient = (ele%value(gradient$) + ele%value(gradient_err$)) * ele%value(field_scale$)
    else
      gradient = ele%value(voltage$) * ele%value(field_scale$) / ele%value(l$)
    endif

    if (.not. ele%is_on) gradient = 0
    gradient = (gradient + gradient_shift_sr_wake(ele, param)) * sign_charge
    gradient = gradient * ele%value(l$) / ele%value(l_hard_edge$)

    omega = twopi * ele%value(rf_frequency$)
    k_wave = omega / c_light

    s_hard_offset = (ele%value(l$) - ele%value(l_hard_edge$)) / 2  ! Relative to entrance end of the cavity
    s_eff = s_rel - s_hard_offset
    if (s_eff < 0 .or. s_eff > ele%value(l_hard_edge$)) return  ! Zero field outside

    beta_start = ele%value(p0c_start$) / ele%value(e_tot_start$)
    t_eff = time - s_hard_offset / (c_light * beta_start)

    E_r_coef = gradient * k_wave * sin(k_wave*s_eff) * cos(omega * t_eff + phase)

    field%E(1) = E_r_coef * x
    field%E(2) = E_r_coef * y
    field%E(3) = 2 * gradient * cos(k_wave * s_eff) * cos(omega * t_eff + phase)
    
    B_phi_coef = -gradient * k_wave * cos(k_wave*s_eff) * sin(omega * t_eff + phase) / c_light 
    field%B(1) = -B_phi_coef * y
    field%B(2) =  B_phi_coef * x

    if (df_calc) then
      call out_io (s_fatal$, r_name, 'dFIELD NOT YET IMPLEMENTED FOR LCAVITY!')
      if (bmad_status%exit_on_error) call err_exit
    endif

  !------------------------------------------
  ! Wiggler

  case(wiggler$)

    n = 0
    if (associated(ele%wig)) n = size(ele%wig%term)

    do i = 1, n
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
  ! E_Gun

  case (e_gun$)
    field%e(3) = sign_charge * ele%value(gradient$) * ele%value(field_scale$)

  !------------------------------------------
  ! Error

  case default
    call out_io (s_fatal$, r_name, 'ELEMENT NOT YET CODED: ' // key_name(ele%key), 'FOR: ' // ele%name)
    if (bmad_status%exit_on_error) call err_exit
    if (present(err_flag)) err_flag = .true.
    return
  end select

  !---------------------------------------------------------------------
  ! Add multipoles

  call multipole_ele_to_ab(ele, param%particle, .not. local_ref_frame, has_nonzero_pole, a_pole, b_pole)
  if (has_nonzero_pole) then

    if (ele%value(l$) == 0) then
      call out_io (s_fatal$, r_name, 'dField NOT YET IMPLEMENTED FOR MULTIPOLES!', 'FOR: ' // ele%name)
      if (bmad_status%exit_on_error) call err_exit
      if (present(err_flag)) err_flag = .true.
      return
    endif

    do i = 0, n_pole_maxx
      if (a_pole(i) == 0 .and. b_pole(i) == 0) cycle
      if (df_calc) then
        call ab_multipole_kick(a_pole(i), b_pole(i), i, local_orb, kx, ky, dkm)
      else
        call ab_multipole_kick(a_pole(i), b_pole(i), i, local_orb, kx, ky)
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
      if (present(err_flag)) err_flag = .true.
      return  
    endif

    E_rho = 0; E_phi = 0; E_z = 0
    B_rho = 0; B_phi = 0; B_z = 0

    radius = sqrt(x**2 + y**2)
    phi = atan2(y, x)

    ! Notice that it is mode%dphi0_ref that is used below. Not ele%value(dphi0_ref$).

    freq = ele%value(rf_frequency$) * ele%em_field%mode(1)%harmonic
    if (freq == 0) then
      call out_io (s_fatal$, r_name, 'Frequency is zero for map in cavity: ' // ele%name)
      if (ele%em_field%mode(1)%harmonic == 0) &
            call out_io (s_fatal$, r_name, '   ... due to harmonic = 0')
      if (bmad_status%exit_on_error) call err_exit
      if (present(err_flag)) err_flag = .true.
      return  
    endif
    t_ref = (ele%value(phi0$) + ele%value(dphi0$) + ele%value(phi0_err$)) / freq
    if (ele%key == rfcavity$) t_ref = 0.25/freq - t_ref

    do i = 1, size(ele%em_field%mode)
      mode => ele%em_field%mode(i)
      m = mode%m

      k_t = twopi * ele%value(rf_frequency$) * mode%harmonic / c_light

      Er = 0; Ep = 0; Ez = 0
      Br = 0; Bp = 0; Bz = 0

      select case (mode%map%ele_anchor_pt)
      case (anchor_beginning$); s0 = 0
      case (anchor_center$);    s0 = ele%value(l$) / 2
      case (anchor_end$);       s0 = ele%value(l$)
      case default
        call out_io (s_fatal$, r_name, 'BAD ELE_ANCHOR_PT FOR FIELD MODE IN ELEMENT: ' // ele%name)
        if (bmad_status%exit_on_error) call err_exit
      end select

      do n = 1, size(mode%map%term)

        term => mode%map%term(n)
        k_zn = twopi * (n - 1) / (size(mode%map%term) * mode%map%dz)
        if (2 * n > size(mode%map%term)) k_zn = k_zn - twopi / mode%map%dz

        expi = cmplx(cos(k_zn * (s_rel-s0)), sin(k_zn * s_rel))

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
          cm = expi * cos(m * phi - mode%phi0_azimuth)
          sm = expi * sin(m * phi - mode%phi0_azimuth)
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

      ! Notice that phi0, dphi0, and phi0_err are folded into t_ref above.

      freq = ele%value(rf_frequency$) * mode%harmonic
      expt = mode%field_scale * exp(-I_imaginary * twopi * (freq * (time + t_ref) + mode%dphi0_ref))
      if (mode%master_scale > 0) expt = expt * ele%value(mode%master_scale)
      E_rho = E_rho + Er * expt
      E_phi = E_phi + Ep * expt
      E_z   = E_z   + Ez * expt

      expt = expt / (twopi * freq)
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
    if (present(err_flag)) err_flag = .true.
    return
  endif
  
  ! radial coordinate
  r = sqrt(x**2 + y**2)

  ! reference time for oscillating elements
  select case (ele%key)
  case(rfcavity$, lcavity$) 
    ! Notice that it is mode%dphi0_ref that is used below. Not ele%value(dphi0_ref$).
    freq = ele%value(rf_frequency$) * ele%em_field%mode(1)%harmonic
    if (freq == 0) then
      call out_io (s_fatal$, r_name, 'Frequency is zero for grid in cavity: ' // ele%name)
      if (ele%em_field%mode(1)%harmonic == 0) &
            call out_io (s_fatal$, r_name, '   ... due to harmonic = 0')
      if (bmad_status%exit_on_error) call err_exit
      if (present(err_flag)) err_flag = .true.
      return  
    endif
    t_ref = (ele%value(phi0$) + ele%value(dphi0$) + ele%value(phi0_err$)) / freq
    if (ele%key == rfcavity$) t_ref = 0.25/freq - t_ref

  case default
    t_ref = 0
  end select

  ! Loop over modes
  do i = 1, size(ele%em_field%mode)
    mode => ele%em_field%mode(i)
    m = mode%m
    
    select case (mode%grid%ele_anchor_pt)
    case (anchor_beginning$); s0 = 0
    case (anchor_center$);    s0 = ele%value(l$) / 2
    case (anchor_end$);       s0 = ele%value(l$)
    case default
      call out_io (s_fatal$, r_name, 'BAD ELE_ANCHOR_PT FOR FIELD GRID IN ELEMENT: ' // ele%name)
      if (bmad_status%exit_on_error) call err_exit
      if (present(err_flag)) err_flag = .true.
      return
    end select

    ! DC modes should have mode%harmonic = 0
    freq = ele%value(rf_frequency$) * mode%harmonic
    expt = mode%field_scale * exp(-I_imaginary * twopi * (freq * (time + t_ref) + mode%dphi0_ref))
    if (mode%master_scale > 0) expt = expt * ele%value(mode%master_scale)

    ! Check for grid
    if (.not. associated(mode%grid)) then
      call out_io (s_fatal$, r_name, 'MISSING GRID FOR ELE: ' // ele%name)
      if (bmad_status%exit_on_error) call err_exit
      if (present(err_flag)) err_flag = .true.
      return
    endif

    ! calculate field based on grid type
    select case(mode%grid%type)
    
    case(rotationally_symmetric_rz$)
      
      ! Format should be: pt (ir, iz) = ( Er, 0, Ez, 0, Bphi, 0 ) 
        
      ! Interpolate 2D (r, z) grid
      ! local_field is a em_field_pt_struct, which has complex E and B

      call em_grid_linear_interpolate(ele, mode%grid, local_field, err, r, s_rel-s0)
      if (err) then
        if (bmad_status%exit_on_error) call err_exit
        if (present(err_flag)) err_flag = .true.
        return
      endif

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
      if (present(err_flag)) err_flag = .true.
      return
    end select
  enddo

!----------------------------------------------------------------------------
! Unknown field calc

case default
  call out_io (s_fatal$, r_name, 'BAD FIELD_CALC METHOD FOR ELEMENT: ' // ele%name)
  if (bmad_status%exit_on_error) call err_exit
  if (present(err_flag)) err_flag = .true.
  return
end select

call convert_fields_to_lab_coords

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! convert fields to lab coords

contains

subroutine convert_fields_to_lab_coords

if (local_ref_frame) return

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

end subroutine convert_fields_to_lab_coords

end subroutine em_field_calc 

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine em_grid_linear_interpolate (ele, grid, field, err_flag, x1, x2, x3)
!
! Subroutine to interpolate the E and B fields on a rectilinear grid
!
! Note: No error checking is done for providing x2 or x3 for 2D and 3D calculations!
!
! Modules needed:
!   use bmad_struct
!
! Input:
!   ele      -- ele_struct: Element containing the grid
!   grid     -- em_field_grid_struct: Grid to interpolate
!   err_flag -- Logical: Set to true if there is an error. False otherwise.
!   x1       -- real(rp) : dimension 1 interpolation point
!   x2       -- real(rp), optional : dimension 2 interpolation point
!   x3       -- real(rp), optional : dimension 3 interpolation point
!
! Output:
!   field  -- em_field_pt_struct: Interpolated field (complex)
!-

subroutine em_grid_linear_interpolate (ele, grid, field, err_flag, x1, x2, x3)

type (ele_struct) ele
type (em_field_grid_struct) :: grid
type (em_field_grid_pt_struct), intent(out) :: field
real(rp) :: x1
real(rp), optional :: x2, x3
real(rp) rel_x1, rel_x2, rel_x3
integer i1, i2, i3, grid_dim
logical out_of_bounds1, out_of_bounds2, out_of_bounds3, err_flag

character(32), parameter :: r_name = 'em_grid_linear_interpolate'

! Pick appropriate dimension 

err_flag = .false.

grid_dim = em_grid_dimension(grid%type)
select case(grid_dim)

case (1)

  call get_this_index(x1, 1, i1, rel_x1, out_of_bounds1, err_flag); if (err_flag) return
  if (out_of_bounds1) return

  field%E(:) = (1-rel_x1) * grid%pt(i1, 1, 1)%E(:) + (rel_x1) * grid%pt(i1+1, 1, 1)%E(:) 
  field%B(:) = (1-rel_x1) * grid%pt(i1, 1, 1)%B(:) + (rel_x1) * grid%pt(i1+1, 1, 1)%B(:) 

case (2)

  call get_this_index(x1, 1, i1, rel_x1, out_of_bounds1, err_flag); if (err_flag) return
  call get_this_index(x2, 2, i2, rel_x2, out_of_bounds2, err_flag); if (err_flag) return
  if (out_of_bounds1 .or. out_of_bounds2) return

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

  call get_this_index(x1, 1, i1, rel_x1, out_of_bounds1, err_flag); if (err_flag) return
  call get_this_index(x2, 2, i2, rel_x2, out_of_bounds2, err_flag); if (err_flag) return
  call get_this_index(x3, 3, i3, rel_x3, out_of_bounds3, err_flag); if (err_flag) return
  if (out_of_bounds1 .or. out_of_bounds2 .or. out_of_bounds3) return
    
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
  err_flag = .true.
  return
end select

!-------------------------------------------------------------------------------------
contains

subroutine get_this_index (x, ix_x, i0, rel_x0, out_of_bounds, err_flag)

implicit none

real(rp) x, rel_x0, x_norm
integer ix_x, i0, ig0, ig1, idg
logical out_of_bounds, err_flag

!

ig0 = lbound(grid%pt, ix_x)
ig1 = ubound(grid%pt, ix_x)

x_norm = (x - grid%r0(ix_x)) / grid%dr(ix_x)
i0 = floor(x_norm)     ! index of lower 1 data point
rel_x0 = x_norm - i0   ! Relative distance from lower x1 grid point

if (i0 == ig1 .and. rel_x0 < bmad_com%significant_length) then
  i0 = ig1 - 1
  rel_x0 = 1
endif

if (i0 == ig0 - 1 .and. abs(rel_x0 - 1) < bmad_com%significant_length) then
  i0 = ig0
  rel_x0 = 0
endif

! Outside of the gird the field is considered to be zero.

! Only generate a warning message if the particle is grossly outside of the grid region.
! Here "gross" is defined as dOut > L_grid/2 where dOut is the distance between the
! particle and the grid edge and L_grid is the length of the grid.

out_of_bounds = .false.

if (i0 < ig0 .or. i0 >= ig1) then
  field%E = 0
  field%B = 0
  out_of_bounds = .true.
  idg = ig1 - ig0

  if (abs(i0 - idg/2) > idg) then
    err_flag = .true.
    call out_io (s_error$, r_name, '\i0\D GRID interpolation index out of bounds: i\i0\ = \i0\ ', &
                                'For element: ' // ele%name, &
                                'Setting field to zero', i_array = [grid_dim, ix_x, i0])
  endif
endif

end subroutine get_this_index 

end subroutine em_grid_linear_interpolate

end module
