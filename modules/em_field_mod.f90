!+
! Module em_field_mod
!
! Module to define the electric and magnetic fields for an elemet.
!-

module em_field_mod

use bmad_struct
use bmad_interface

implicit none

contains

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!+
! Function g_bend_from_em_field (B, E, orbit) result (g_bend)
!
! Routine to calculate the bending strength (1/bending_radius) for
! a given particle for a given field.
!
! Input:
!   B(3)  -- real(rp): Magnetic field.
!   E(3)  -- real(rp): Electric field
!   orbit -- coord_struct: particle orbit
!
! Output:
!   g_bend(3) -- real(rp): bending strength vector.
!-

function g_bend_from_em_field (B, E, orbit) result (g_bend)

type (coord_struct) orbit
real(rp) b(3), e(3), g_bend(3)
real(rp) vel(3), rel_pc, force(3)

! vel is normalized velocity

rel_pc = 1 + orbit%vec(6)
vel(1:2) = [orbit%vec(2), orbit%vec(4)] / rel_pc
vel(3) = sqrt(1 - vel(1)**2 - vel(2)**2) * orbit%direction

force = (E + cross_product(vel, B) * orbit%beta * c_light) * charge_of(orbit%species)
g_bend = -(force - vel * (dot_product(force, vel))) / (orbit%p0c * rel_pc)

end function g_bend_from_em_field

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!+
! Subroutine save_a_step (track, ele, param, local_ref_frame, s, here, s_sav, t_ref)
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
!   t_ref    -- real(rp), optional: Reference time. If present then track%field call be evaluated.
!
! Ouput:
!   track    -- track_struct: Trajectory structure to save to.
!   s_sav    -- Real(rp): Set equal to s.
!-

subroutine save_a_step (track, ele, param, local_ref_frame, s, orb, s_sav, t_ref)

type (track_struct) track, track2
type (ele_struct), target :: ele
type (lat_param_struct), intent(in) :: param
type (coord_struct) orb, orb2
integer n_pt, n, n_old
real(rp) s, s_sav
real(rp), optional :: t_ref
logical local_ref_frame

! Not allocated

if (.not. allocated (track%orb)) then
  allocate(track%orb(0:100))
  allocate(track%field(0:100))
  allocate(track%map(0:100))
  track%n_ok = 0
  track%n_bad = 0
  track%n_pt = -1
endif

!

track%n_pt = track%n_pt + 1
n_pt = track%n_pt
n_old = ubound(track%orb, 1)

if (n_pt > n_old) then
  n = 1.5 * n_pt
  call move_alloc (track%orb, track2%orb)
  call move_alloc (track%field, track2%field)
  call move_alloc (track%map, track2%map)
  allocate(track%orb(0:n), track%field(0:n), track%map(0:n))
  track%orb(:n_old) = track2%orb; track%field(:n_old) = track2%field; track%map(:n_old) = track2%map
  deallocate(track2%orb, track2%field, track2%map)
end if

! Notice that a translation due to a finite ele%value(z_offset$) is not wanted here.

orb2 = orb
if (local_ref_frame) call offset_particle (ele, param, unset$, orb2, set_z_offset = .false., ds_pos = s)

track%orb(n_pt) = orb2
track%orb(n_pt)%ix_ele = ele%ix_ele
track%map(n_pt)%mat6 = 0

s_sav = s

if (present(t_ref)) then
  call em_field_calc (ele, param, s, t_ref, orb, local_ref_frame, track%field(n_pt), .false.)
endif

end subroutine save_a_step

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine em_field_calc (ele, param, s_pos, time, orbit, local_ref_frame, field, calc_dfield, err_flag, potential)
!
! Subroutine to calculate the E and B fields for an element.
!
! s_pos is measured from the upstream edge of the element. For elements with ele%orientation = -1, s_pos
! is measured from the element's exit end. However, the field is always referenced to the element coordinates.
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
!   s_pos  -- Real(rp): Longitudinal position relative to the upstream edge of the element.
!   time   -- Real(rp): Particle time.
!                 For absolute time tracking this is the absolute time.
!                 For relative time tracking this is relative to the reference particle entering the element.
!   orbit  -- Coord_struct: Transverse coordinates.
!     %vec(1), %vec(3)  -- Transverse coords. These are the only components used in the calculation.
!   local_ref_frame 
!          -- Logical, If True then take the input coordinates and output fields 
!                as being with respect to the frame of referene of the element (ignore misalignments). 
!   calc_dfield     
!          -- Logical, optional: If present and True 
!                then calculate the field derivatives.
!
! Output:
!   field       -- em_field_struct: E and B fields and derivatives.
!   err_flag    -- logical, optional: Set True if there is an error. False otherwise.
!   potential   -- em_potential_struct, optional: The electric and magnetic potentials.
!                   This is experimental and only implemented for wigglers at present.
!-

recursive subroutine em_field_calc (ele, param, s_pos, time, orbit, local_ref_frame, field, &
                                                                           calc_dfield, err_flag, potential)

use geometry_mod

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
type (lat_param_struct) param
type (coord_struct) :: orbit, local_orb
type (em_potential_struct), optional :: potential
type (wig_term_struct), pointer :: wig
type (em_field_struct) :: field, field2
type (em_field_grid_pt_struct) :: local_field
type (em_field_mode_struct), pointer :: mode
type (em_field_map_term_struct), pointer :: term

real(rp) :: x, x_save, y, s, t, time, s_pos, s_rel, z, f, dk(3,3), ref_charge, f_p0c
real(rp) :: c_x, s_x, c_y, s_y, c_z, s_z, coef, fd(3), s0, Ex, Ey
real(rp) :: cos_ang, sin_ang, sgn_x, sgn_y, kx, ky, dkm(2,2), cos_ks, sin_ks
real(rp) phase, gradient, r, E_r_coef, E_s, k_wave, s_eff, t_eff
real(rp) k_t, k_zn, kappa2_n, kap_rho, s_hard_offset, beta_start
real(rp) radius, phi, t_ref, tilt, omega, freq0, freq, B_phi_coef
real(rp) Er_dc, Ep_dc, Ez_dc, Br_dc, Bp_dc, Bz_dc
real(rp) E_rho, E_phi, E_z, B_rho, B_phi, B_z, sx_over_kx, sy_over_ky
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)

complex(rp) Er, Ep, Ez, Br, Bp, Bz
complex(rp) exp_kz, exp_m, expt, dEp, dEr
complex(rp) Im_0, Im_plus, Im_minus, Im_0_R, kappa_n, Im_plus2, cm, sm, q

integer i, j, m, n, trig_x, trig_y

logical :: local_ref_frame, local_ref, has_nonzero_pole
logical, optional :: calc_dfield, err_flag
logical df_calc, err

character(20) :: r_name = 'em_field_calc'

! Initialize field
! If element is turned off then return zero

field = em_field_struct()
if (present(potential)) potential = em_potential_struct()

df_calc = logic_option (.false., calc_dfield)

if (present(err_flag)) err_flag = .false.
if (.not. ele%is_on) return

if (ele%orientation == 1) then
  s_rel = s_pos
else
  s_rel = ele%value(l$) - s_pos
endif

!----------------------------------------------------------------------------
! convert to local coords

local_orb = orbit
if (.not. local_ref_frame) then
  call offset_particle (ele, param, set$, local_orb, set_multipoles = .false., set_hvkicks = .false., ds_pos = s_rel)
endif

!----------------------------------------------------------------------------
! super_slave, and slice_slave, have their field info stored in the associated lord elements.
! Note: The lord of an em_field element has independent misalignments.
! Note: multipass_slave elements do store their own field info. This should be changed.

if (ele%field_calc == refer_to_lords$) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%field_calc == no_field$) cycle   ! Group, overlay and girder elements do not have fields.

    local_ref = .true.
    if (ele%key == em_field$) local_ref = .false.

    s = s_pos + (ele%s - ele%value(l$)) - (lord%s - lord%value(l$))
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
  call em_field_custom (ele, param, s_pos, time, orbit, local_ref_frame, field, calc_dfield)
  return
end if

!----------------------------------------------------------------------------
!Set up common variables for all (non-custom) methods

ref_charge = charge_of(param%particle)

x = local_orb%vec(1)
y = local_orb%vec(3)

if (ref_charge == 0) then
  f_p0c = 0
else
  f_p0c = ele%value(p0c$) / (c_light * ref_charge)
endif

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! field_calc methods

select case (ele%field_calc)
  
!----------------------------------------------------------------------------
! Bmad_standard field calc 

case (bmad_standard$)

  select case (ele%key)

  !------------------------------------------
  ! Drift, et. al. Note that kicks get added at the end for all elements

  case (drift$, ecollimator$, rcollimator$, instrument$, monitor$, pipe$, marker$, detector$)

  !------------------------------------------
  ! E_Gun

  case (e_gun$)
    field%e(3) = e_accel_field (ele, gradient$)

  !------------------------------------------
  ! Elseparator

  case (elseparator$)
    field%e(1) = ele%value(hkick$) * ele%value(p0c$) / ele%value(l$)
    field%e(2) = ele%value(vkick$) * ele%value(p0c$) / ele%value(l$)

  !------------------------------------------
  ! HKicker

  case (hkicker$)
    field%b(2) = -ele%value(kick$) * f_p0c / ele%value(l$)

  !------------------------------------------
  ! Kicker  

  case (kicker$)
    field%b(1) =  ele%value(vkick$) * f_p0c / ele%value(l$)
    field%b(2) = -ele%value(hkick$) * f_p0c / ele%value(l$)

  !------------------------------------------
  ! RFcavity and Lcavity  bmad_standard
  !
  ! For standing wave cavity:
  ! Use N_cell half-wave pillbox formulas for TM_011 mode with infinite wall radius.
  ! See S.Y. Lee, "Accelerator Physics"
  !   E_s   = 2 * gradient * cos(k s) * cos(omega t + phase)
  !   E_r   =     gradient * k * r * sin(k s) * cos(omega t + phase)
  !   B_phi =    -gradient * k * r * cos(k s) * sin(omega t + phase) / c_light
  ! where
  !   k = n_cell * pi / L
  !   omega = c * k
  ! Field extends to +/- c_light * freq / 2 from centerline of element.
  ! Note: There is a discontinuity in the field at the edge. Edge focusing due to this 
  !  discontinuity can be handled in the apply_element_edge_kick routine.
  !
  ! For traveling wave cavity:
  !   E_s   =  gradient * cos(omega t + phase - k s)
  !   E_r   = -gradient * k * r * sin(omega t + phase - k s) / 2
  !   B_phi = -gradient * k * r * sin(omega t + phase - k s) / c_light / 2

  case(rfcavity$, lcavity$)

    if (ele%value(rf_frequency$) == 0) return

    phase = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_err$) + ele%value(phi0_ref$))
    if (ele%key == rfcavity$) phase = pi/2 - phase
    orbit%phase(1) = phase  ! RF phase is needed by apply_element_edge_kick when calling rf_coupler_kick.

    gradient = e_accel_field (ele, gradient$)

    if (.not. ele%is_on) gradient = 0
    gradient = (gradient + gradient_shift_sr_wake(ele, param)) / ref_charge
    gradient = gradient * ele%value(l$) / hard_edge_model_length(ele)
    omega = twopi * ele%value(rf_frequency$)
    k_wave = omega / c_light

    s_hard_offset = (ele%value(l$) - hard_edge_model_length(ele)) / 2  ! Relative to entrance end of the cavity
    s_eff = s_rel - s_hard_offset
    if (s_eff < 0 .or. s_eff > hard_edge_model_length(ele)) return  ! Zero field outside

    beta_start = ele%value(p0c_start$) / ele%value(e_tot_start$)
    t_eff = time - s_hard_offset / (c_light * beta_start)

    if (nint(ele%value(cavity_type$)) == traveling_wave$) then
      phi = omega * t_eff + phase - k_wave * s_eff
      E_z        =  gradient * cos(phi)
      E_r_coef   = -gradient * k_wave * sin(phi) / 2
      B_phi_coef = -gradient * k_wave * sin(phi) / c_light / 2
    else
      E_z        = 2 * gradient * cos(k_wave * s_eff) * cos(omega * t_eff + phase)
      E_r_coef   =     gradient * k_wave * sin(k_wave*s_eff) * cos(omega * t_eff + phase)
      B_phi_coef =    -gradient * k_wave * cos(k_wave*s_eff) * sin(omega * t_eff + phase) / c_light 
    endif

    field%E(1) = E_r_coef * x
    field%E(2) = E_r_coef * y
    field%E(3) = E_z
    
    field%B(1) = -B_phi_coef * y
    field%B(2) =  B_phi_coef * x

    if (df_calc) then
      call out_io (s_fatal$, r_name, 'dFIELD NOT YET IMPLEMENTED FOR LCAVITY AND RFCAVITY!')
      if (global_com%exit_on_error) call err_exit
    endif

  !------------------------------------------
  ! Octupole 

  case (octupole$)

    field%b(1) = -(y**3 - 3*y*x**2) / 6 * ele%value(k3$) * f_p0c 
    field%b(2) =  (x**3 - 3*x*y**2) / 6 * ele%value(k3$) * f_p0c 

    if (df_calc) then
      field%dB(1,1) =  x*y * ele%value(k3$) * f_p0c
      field%dB(1,2) = (x**2 - y**2) / 2 * ele%value(k3$) * f_p0c
      field%dB(2,1) = (x**2 - y**2) / 2 * ele%value(k3$) * f_p0c
      field%dB(2,2) = -x*y * ele%value(k3$) * f_p0c
    endif

  !------------------------------------------
  ! Patch: There are no fields

  case (patch$)

    return

  !------------------------------------------
  ! Quadrupole

  case (quadrupole$) 

    field%b(1) = y * ele%value(k1$) * f_p0c 
    field%b(2) = x * ele%value(k1$) * f_p0c 

    if (df_calc) then
      field%dB(1,2) =  ele%value(k1$) * f_p0c
      field%dB(2,1) =  ele%value(k1$) * f_p0c
    endif

  !------------------------------------------
  ! Sextupole 

  case (sextupole$)

    field%b(1) = x * y * ele%value(k2$) * f_p0c
    field%b(2) = (x**2 - y**2) / 2 * ele%value(k2$) * f_p0c 

    if (df_calc) then
      field%dB(1,1) =  y * ele%value(k2$) * f_p0c
      field%dB(1,2) =  x * ele%value(k2$) * f_p0c
      field%dB(2,1) =  x * ele%value(k2$) * f_p0c
      field%dB(2,2) = -y * ele%value(k2$) * f_p0c
    endif

  !------------------------------------------
  ! VKicker

  case (vkicker$)
    field%b(1) =  ele%value(kick$) * f_p0c / ele%value(l$)

  !------------------------------------------
  ! SBend

  case (sbend$)

    field%b(1) = (y * ele%value(k1$) + x * y * ele%value(k2$)) * f_p0c 
    field%b(2) = (x * ele%value(k1$) + ele%value(k2$) * (x**2 - y**2) / 2 + ele%value(g$) + ele%value(g_err$)) * f_p0c 

    if (df_calc) then
      field%dB(1,1) =  y * ele%value(k2$) * f_p0c
      field%dB(1,2) =  (x * ele%value(k2$) + ele%value(k1$)) * f_p0c
      field%dB(2,1) =  (x * ele%value(k2$) + ele%value(k1$)) * f_p0c
      field%dB(2,2) = -y * ele%value(k2$) * f_p0c
    endif

  !------------------------------------------
  ! Sol_quad

  case (sol_quad$)

    field%b(1) = y * ele%value(k1$) * f_p0c 
    field%b(2) = x * ele%value(k1$) * f_p0c 
    field%b(3) = ele%value(ks$) * f_p0c

    if (df_calc) then
      field%dB(1,2) = ele%value(k1$) * f_p0c
      field%dB(2,1) = ele%value(k1$) * f_p0c
    endif

  !------------------------------------------
  ! Solenoid

  case (solenoid$)

    field%b(3) = ele%value(ks$) * f_p0c

    if (df_calc) then
    endif

  !------------------------------------------
  ! Wiggler

  case(wiggler$, undulator$)

    n = 0
    if (associated(ele%wig)) n = size(ele%wig%term)

    do i = 1, n
      wig => ele%wig%term(i)
      sgn_x = 1; sgn_y = 1

      select case (wig%type)
      case (hyper_y_plane_x$, hyper_y_plane_y$)
        coef = wig%coef * ele%value(polarity$) / ref_charge / wig%ky
        c_x = cos(wig%kx * (x + wig%x0))
        s_x = sin(wig%kx * (x + wig%x0))
        c_y = cosh (wig%ky * (y + wig%y0))
        s_y = sinh (wig%ky * (y + wig%y0))
        if (wig%type == hyper_y_plane_y$) sgn_x = -1
        trig_x = -1; trig_y = 1

      case (hyper_xy_plane_x$, hyper_xy_plane_y$)
        coef = wig%coef * ele%value(polarity$) / ref_charge / wig%kz
        c_x = cosh(wig%kx * (x + wig%x0))
        s_x = sinh(wig%kx * (x + wig%x0))
        c_y = cosh (wig%ky * (y + wig%y0))
        s_y = sinh (wig%ky * (y + wig%y0))
        trig_x = 1; trig_y = 1

      case (hyper_x_plane_x$, hyper_x_plane_y$)
        coef = wig%coef * ele%value(polarity$) / ref_charge / wig%kx
        c_x = cosh(wig%kx * (x + wig%x0))
        s_x = sinh(wig%kx * (x + wig%x0))
        c_y = cos (wig%ky * (y + wig%y0))
        s_y = sin (wig%ky * (y + wig%y0))
        if (wig%type == hyper_x_plane_x$) sgn_y = -1
        trig_x = 1; trig_y = -1
      end select

      c_z = cos (wig%kz * s_rel + wig%phi_z)
      s_z = sin (wig%kz * s_rel + wig%phi_z)

      select case (wig%type)
      case (hyper_y_plane_x$, hyper_xy_plane_x$, hyper_x_plane_x$)
        field%B(1) = field%B(1) + coef  * wig%kx * c_x * c_y * c_z
        field%B(2) = field%B(2) + coef  * wig%ky * s_x * s_y * c_z * sgn_y
        field%B(3) = field%B(3) - coef  * wig%kz * s_x * c_y * s_z
      case default
        field%B(1) = field%B(1) + coef  * wig%kx * s_x * s_y * c_z * sgn_x
        field%B(2) = field%B(2) + coef  * wig%ky * c_x * c_y * c_z
        field%B(3) = field%B(3) - coef  * wig%kz * c_x * s_y * s_z
      end select

      if (df_calc) then
        select case (wig%type)
        case (hyper_y_plane_x$, hyper_xy_plane_x$, hyper_x_plane_x$)
          f = coef * wig%kx
          field%dB(1,1) = field%dB(1,1) + f * wig%kx * s_x * c_y * c_z * trig_x
          field%dB(2,1) = field%dB(2,1) + f * wig%ky * c_x * s_y * c_z * sgn_y
          field%dB(3,1) = field%dB(3,1) - f * wig%kz * c_x * c_y * s_z 
          f = coef * wig%ky
          field%dB(1,2) = field%dB(1,2) + f * wig%kx * c_x * s_y * c_z * trig_y
          field%dB(2,2) = field%dB(2,2) + f * wig%ky * s_x * c_y * c_z * sgn_y
          field%dB(3,2) = field%dB(3,2) - f * wig%kz * s_x * s_y * s_z * trig_y
          f = coef * wig%kz
          field%dB(1,3) = field%dB(1,3) - f * wig%kx * c_x * c_y * s_z
          field%dB(2,3) = field%dB(2,3) - f * wig%ky * s_x * s_y * s_z * sgn_y 
          field%dB(3,3) = field%dB(3,3) - f * wig%kz * s_x * c_y * c_z 
        case default
          f = coef * wig%kx
          field%dB(1,1) = field%dB(1,1) + f * wig%kx * c_x * s_y * c_z * sgn_x
          field%dB(2,1) = field%dB(2,1) + f * wig%ky * s_x * c_y * c_z * trig_x
          field%dB(3,1) = field%dB(3,1) - f * wig%kz * s_x * s_y * s_z * trig_x
          f = coef * wig%ky
          field%dB(1,2) = field%dB(1,2) + f * wig%kx * s_x * c_y * c_z * sgn_x
          field%dB(2,2) = field%dB(2,2) + f * wig%ky * c_x * s_y * c_z * trig_y
          field%dB(3,2) = field%dB(3,2) - f * wig%kz * c_x * c_y * s_z 
          f = coef * wig%kz
          field%dB(1,3) = field%dB(1,3) - f * wig%kx * s_x * s_y * s_z * sgn_x
          field%dB(2,3) = field%dB(2,3) - f * wig%ky * c_x * c_y * s_z 
          field%dB(3,3) = field%dB(3,3) - f * wig%kz * c_x * s_y * c_z 
        end select
      endif

      if (present(potential)) then
        coef = wig%coef * ele%value(polarity$) / ref_charge
        select case (wig%type)
        case (hyper_y_plane_x$, hyper_xy_plane_x$, hyper_x_plane_x$)
          if (abs(wig%ky * (y + wig%y0)) < 1d-10) then
            sy_over_ky = y + wig%y0
          else
            sy_over_ky = s_y / wig%ky
          endif
        case (hyper_y_plane_y$, hyper_xy_plane_y$, hyper_x_plane_y$)
          if (abs(wig%kx * (x + wig%x0)) < 1d-10) then
            sx_over_kx = x + wig%x0
          else
            sx_over_kx = s_x / wig%kx
          endif
        end select

        select case (wig%type)
        case (hyper_y_plane_x$)
          potential%a(1) = potential%a(1) + coef * s_x * sy_over_ky * s_z * wig%kz / wig%ky
          potential%a(3) = potential%a(3) + coef * c_x * sy_over_ky * c_z * wig%kx / wig%ky
        case (hyper_xy_plane_x$)
          potential%a(1) = potential%a(1) + coef * s_x * sy_over_ky * s_z
          potential%a(3) = potential%a(3) + coef * c_x * sy_over_ky * c_z * wig%kx / wig%kz
        case (hyper_x_plane_x$)
          potential%a(1) = potential%a(1) + coef * s_x * sy_over_ky * s_z * wig%kz / wig%kx
          potential%a(3) = potential%a(3) + coef * c_x * sy_over_ky * c_z
        case (hyper_y_plane_y$)
          potential%a(2) = potential%a(2) - coef * sx_over_kx * s_y * s_z * wig%kz / wig%ky
          potential%a(3) = potential%a(3) - coef * sx_over_kx * c_y * c_z
        case (hyper_xy_plane_y$)
          potential%a(2) = potential%a(2) - coef * sx_over_kx * s_y * s_z
          potential%a(3) = potential%a(3) - coef * sx_over_kx * c_y * c_z * wig%ky / wig%kz
        case (hyper_x_plane_y$)
          potential%a(2) = potential%a(2) - coef * sx_over_kx * s_y * s_z * wig%kz / wig%kx
          potential%a(3) = potential%a(3) - coef * sx_over_kx * c_y * c_z * wig%ky / wig%kx
        end select
      endif

    enddo

  !------------------------------------------
  ! Error

  case default
    call out_io (s_fatal$, r_name, 'EM FIELD NOT YET CODED FOR ELEMENT OF TYPE: ' // key_name(ele%key), &
                                   'FOR ELEMENT: ' // ele%name, 'PERHAPS "FIELD_CALC" NEEDS TO BE SET FOR THIS ELEMENT?')
    if (global_com%exit_on_error) call err_exit
    if (present(err_flag)) err_flag = .true.
    return
  end select

  !---------------------------------------------------------------------
  ! Add multipoles

  call multipole_ele_to_ab(ele, .not. local_ref_frame, has_nonzero_pole, a_pole, b_pole)
  if (has_nonzero_pole) then

    if (ele%value(l$) == 0) then
      call out_io (s_fatal$, r_name, 'dField NOT YET IMPLEMENTED FOR MULTIPOLES!', 'FOR: ' // ele%name)
      if (global_com%exit_on_error) call err_exit
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
      field%B(1) = field%B(1) + f_p0c * ky / ele%value(l$)
      field%B(2) = field%B(2) - f_p0c * kx / ele%value(l$)
      if (df_calc) then
        field%dB(1,1) = field%dB(1,1) + f_p0c * dkm(2,1) / ele%value(l$)
        field%dB(1,2) = field%dB(1,2) + f_p0c * dkm(2,2) / ele%value(l$)
        field%dB(2,1) = field%dB(2,1) - f_p0c * dkm(1,1) / ele%value(l$)
        field%dB(2,2) = field%dB(2,2) - f_p0c * dkm(1,2) / ele%value(l$)
      endif
    enddo

  endif

  !---------------------------------------------------------------------
  ! Add electric multipoles

  call multipole_ele_to_ab(ele, .not. local_ref_frame, has_nonzero_pole, a_pole, b_pole, electric$)
  if (has_nonzero_pole) then
    do i = 0, n_pole_maxx
      if (a_pole(i) == 0 .and. b_pole(i) == 0) cycle
      call elec_multipole_field(a_pole(i), b_pole(i), i, local_orb, Ex, Ey, dkm, df_calc)
      field%E(1) = field%E(1) + Ex
      field%E(2) = field%E(2) + Ey
      if (df_calc) field%dE(1:2,1:2) = field%dE(1:2,1:2) + dkm
    enddo
  endif

  !-------------------------------
  ! Add kicks. Since the kick field is not rotated by a tilt then we have to unrotate if in the local_ref_frame

  if (has_hkick_attributes(ele%key) .and. (ele%value(hkick$) /= 0 .or. ele%value(vkick$) /= 0)) then
    select case (ele%key)
    ! Handled above
    case (kicker$, hkicker$, vkicker$, elseparator$)  
    ! Everything else
    case default
      if (.not. local_ref_frame .or. ele%value(tilt_tot$) == 0) then
        field%b(1) = field%b(1) + ele%value(Vkick$) * f_p0c / ele%value(l$)
        field%b(2) = field%b(2) - ele%value(Hkick$) * f_p0c / ele%value(l$)
      else
        ! Rotate from lab to local
        tilt = ele%value(tilt_tot$)
        if (ele%key == sbend$) tilt = ele%value(ref_tilt_tot$)
        field%b(1) = field%b(1) + (ele%value(Vkick$) * cos(tilt) - ele%value(hkick$) * sin(tilt)) * f_p0c / ele%value(l$)
        field%b(2) = field%b(2) - (ele%value(Hkick$) * cos(tilt) + ele%value(vkick$) * sin(tilt)) * f_p0c / ele%value(l$)
      endif
    end select
  endif

!----------------------------------------------------------------------------
! Map field calc...

case(map$)

  if (.not. associated(ele%em_field)) then
    call out_io (s_fatal$, r_name, 'No accociated em_field for field calc = Map', 'FOR: ' // ele%name) 
    if (global_com%exit_on_error) call err_exit
    if (present(err_flag)) err_flag = .true.
    return  
  endif

  radius = sqrt(x**2 + y**2)
  phi = atan2(y, x)

  ! Notice that it is mode%phi0_ref that is used below. Not ele%value(phi0_ref$).

  select case (ele%key)
  case (lcavity$, rfcavity$, em_field$)
    freq0 = ele%value(rf_frequency$) * ele%em_field%mode(1)%harmonic
    if (freq0 == 0 .and. ele%key /= em_field$) then
      call out_io (s_fatal$, r_name, 'Frequency is zero for map in cavity: ' // ele%name)
      if (ele%em_field%mode(1)%harmonic == 0) call out_io (s_fatal$, r_name, '   ... due to harmonic = 0')
      if (global_com%exit_on_error) call err_exit
      if (present(err_flag)) err_flag = .true.
      return  
    endif
    t_ref = (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_err$)) / freq0
    if (ele%key == rfcavity$) t_ref = 0.25/freq0 - t_ref

  case default ! Where rf_frequency is not defined
    if (any(ele%em_field%mode%harmonic /= 0)) then
      call out_io (s_fatal$, r_name, 'RF FIELD MAPS NOT IMPLEMENTED FOR ELEMENT: ' // ele%name)
      if (global_com%exit_on_error) call err_exit
      return
    endif
  end select

  !

  E_rho = 0; E_phi = 0; E_z = 0
  B_rho = 0; B_phi = 0; B_z = 0

  do i = 1, size(ele%em_field%mode)
    mode => ele%em_field%mode(i)
    m = mode%m

    Er = 0; Ep = 0; Ez = 0
    Br = 0; Bp = 0; Bz = 0

    Er_dc = 0; Ep_dc = 0; Ez_dc = 0
    Br_dc = 0; Bp_dc = 0; Bz_dc = 0

    if (mode%harmonic /= 0) k_t = twopi * ele%value(rf_frequency$) * mode%harmonic / c_light

    select case (mode%map%ele_anchor_pt)
    case (anchor_beginning$); s0 = 0
    case (anchor_center$);    s0 = ele%value(l$) / 2
    case (anchor_end$);       s0 = ele%value(l$)
    case default
      call out_io (s_fatal$, r_name, 'BAD ELE_ANCHOR_PT FOR FIELD MODE IN ELEMENT: ' // ele%name)
      if (global_com%exit_on_error) call err_exit
    end select

    do n = 1, size(mode%map%term)

      term => mode%map%term(n)
      k_zn = twopi * (n - 1) / (size(mode%map%term) * mode%map%dz)
      if (2 * n > size(mode%map%term)) k_zn = k_zn - twopi / mode%map%dz

      cos_ks = cos(k_zn * (s_rel-s0))
      sin_ks = sin(k_zn * (s_rel-s0))
      exp_kz = cmplx(cos_ks, sin_ks, rp)

      ! DC
      if (mode%harmonic == 0) then

        kap_rho = k_zn * radius
        if (m == 0) then
          Im_0    = I_bessel(0, kap_rho)
          Im_plus = I_bessel(1, kap_rho)
          Er_dc = Er_dc + real(term%e_coef * exp_kz * Im_plus)
          Br_dc = Br_dc + real(term%b_coef * exp_kz * Im_plus)
          Ez_dc = Ez_dc + real(term%e_coef * exp_kz * Im_0 * i_imaginary)
          Bz_dc = Bz_dc + real(term%b_coef * exp_kz * Im_0 * i_imaginary)
        else
          Im_plus  = I_bessel(m+1, kap_rho)
          Im_minus = I_bessel(m-1, kap_rho)
          Im_0     = kap_rho * (Im_minus - Im_plus) / (2 * m)
          exp_m = cmplx(cos(m * phi), sin(m * phi), rp)

          q = exp_kz * exp_m * (Im_minus + Im_plus) / 2
          Er_dc = Er_dc + real(term%e_coef * q)
          Br_dc = Br_dc + real(term%b_coef * q)

          q = i_imaginary * exp_kz * exp_m * (Im_minus - Im_plus) / 2
          Ep_dc = Ep_dc + real(term%e_coef * q)
          Bp_dc = Bp_dc + real(term%b_coef * q)

          q = i_imaginary * exp_kz * exp_m * Im_0
          Ez_dc = Ez_dc + real(term%e_coef * q)
          Bz_dc = Bz_dc + real(term%b_coef * q)
        endif

      ! RF mode 
      else
        kappa2_n = k_zn**2 - k_t**2
        kappa_n = sqrt(abs(kappa2_n))
        kap_rho = kappa_n * radius
        if (kappa2_n < 0) then
          kappa_n = -i_imaginary * kappa_n
          kap_rho = -kap_rho
        endif

        if (m == 0) then
          Im_0    = I_bessel_extended(0, kap_rho)
          Im_plus = I_bessel_extended(1, kap_rho) / kappa_n

          Er = Er - term%e_coef * Im_plus * exp_kz * I_imaginary * k_zn
          Ep = Ep + term%b_coef * Im_plus * exp_kz
          Ez = Ez + term%e_coef * Im_0    * exp_kz

          Br = Br - term%b_coef * Im_plus * exp_kz * k_zn
          Bp = Bp - term%e_coef * Im_plus * exp_kz * k_t**2 * I_imaginary
          Bz = Bz - term%b_coef * Im_0    * exp_kz * I_imaginary

        else
          cm = exp_kz * cos(m * phi - mode%phi0_azimuth)
          sm = exp_kz * sin(m * phi - mode%phi0_azimuth)
          Im_plus  = I_bessel_extended(m+1, kap_rho) / kappa_n**(m+1)
          Im_minus = I_bessel_extended(m-1, kap_rho) / kappa_n**(m-1)

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
      endif ! mode%harmonic /= 0
        
    enddo  ! mode%term

    ! Notice that phi0, phi0_multipass, and phi0_err are folded into t_ref above.

    if (mode%harmonic == 0) then
      E_rho = E_rho + Er_dc
      E_phi = E_phi + Ep_dc
      E_z   = E_z   + Ez_dc

      B_rho = B_rho + Br_dc
      B_phi = B_phi + Bp_dc
      B_z   = B_z   + Bz_dc
    else
      freq = ele%value(rf_frequency$) * mode%harmonic
      expt = mode%field_scale * exp(-I_imaginary * twopi * (freq * (time + t_ref) + mode%phi0_ref))
      if (mode%master_scale > 0) expt = expt * ele%value(mode%master_scale)
      E_rho = E_rho + real(Er * expt)
      E_phi = E_phi + real(Ep * expt)
      E_z   = E_z   + real(Ez * expt)

      expt = expt / (twopi * freq)
      B_rho = B_rho + real(Br * expt)
      B_phi = B_phi + real(Bp * expt)
      B_z   = B_z   + real(Bz * expt)
    endif

  enddo

  field%E = [cos(phi) * E_rho - sin(phi) * E_phi, sin(phi) * E_rho + cos(phi) * E_phi, E_z]
  field%B = [cos(phi) * B_rho - sin(phi) * B_phi, sin(phi) * B_rho + cos(phi) * B_phi, B_z]

!----------------------------------------------------------------------------
! Grid field calc 

case(grid$)

  if (.not. associated(ele%em_field)) then
    call out_io (s_fatal$, r_name, 'No accociated em_field for field calc = Grid', 'FOR: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    if (present(err_flag)) err_flag = .true.
    return
  endif
  
  !-------------------
  ! First calc reference time for oscillating elements

  select case (ele%key)
  case(rfcavity$, lcavity$)
    ! Notice that it is mode%phi0_ref that is used below. Not ele%value(phi0_ref$).
    freq = ele%value(rf_frequency$) * ele%em_field%mode(1)%harmonic
    if (freq == 0) then
      call out_io (s_fatal$, r_name, 'Frequency is zero for grid in cavity: ' // ele%name)
      if (ele%em_field%mode(1)%harmonic == 0) &
            call out_io (s_fatal$, r_name, '   ... due to harmonic = 0')
      if (global_com%exit_on_error) call err_exit
      if (present(err_flag)) err_flag = .true.
      return  
    endif
    t_ref = (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_err$)) / freq
    if (ele%key == rfcavity$) t_ref = 0.25/freq - t_ref

  case(e_gun$) 
    ! Same as above, but no error checking for zero frequency
    freq = ele%value(rf_frequency$) * ele%em_field%mode(1)%harmonic
    if (freq == 0) then
      t_ref = 0
    else
      t_ref = (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_err$)) / freq
    endif

  case default
    t_ref = 0
  end select

  !-------------------
  ! Now loop over grid modes

  do i = 1, size(ele%em_field%mode)
    mode => ele%em_field%mode(i)
    m = mode%m

    ! Check for grid
    if (.not. associated(mode%grid)) then
      call out_io (s_fatal$, r_name, 'MISSING GRID FOR ELE: ' // ele%name)
      if (global_com%exit_on_error) call err_exit
      if (present(err_flag)) err_flag = .true.
      return
    endif

    !

    select case (mode%grid%ele_anchor_pt)
    case (anchor_beginning$); s0 = 0
    case (anchor_center$);    s0 = ele%value(l$) / 2
    case (anchor_end$);       s0 = ele%value(l$)
    case default
      call out_io (s_fatal$, r_name, 'BAD ELE_ANCHOR_PT FOR FIELD GRID IN ELEMENT: ' // ele%name)
      if (global_com%exit_on_error) call err_exit
      if (present(err_flag)) err_flag = .true.
      return
    end select

    z = s_rel-s0
       
    ! Sbend grids are in cartesian coordinates
    if (ele%key == sbend$ .and. .not. mode%grid%curved_coords) call convert_curvilinear_to_cartesian()

    ! radial coordinate
    r = sqrt(x**2 + y**2)

    ! DC modes should have mode%harmonic = 0

    if (mode%harmonic == 0) then
      expt = 1
    else
      freq = ele%value(rf_frequency$) * mode%harmonic
      expt = mode%field_scale * exp(-I_imaginary * twopi * (freq * (time + t_ref) + mode%phi0_ref))
    endif

    if (mode%master_scale > 0) expt = expt * ele%value(mode%master_scale)

    ! calculate field based on grid type
    select case(mode%grid%type)

    case (xyz$)
    
      call em_grid_linear_interpolate(ele, mode%grid, local_field, err, x, y, z)
      if (err) then
        if (global_com%exit_on_error) call err_exit
        if (present(err_flag)) err_flag = .true.
        return
      endif

      field%e = field%e + real(expt * local_field%e)
      field%b = field%b + real(expt * local_field%B)


    case(rotationally_symmetric_rz$)
      
      ! Format should be: pt (ir, iz) = ( Er, 0, Ez, 0, Bphi, 0 ) 
        
      ! Interpolate 2D (r, z) grid
      ! local_field is a em_field_pt_struct, which has complex E and B

      call em_grid_linear_interpolate(ele, mode%grid, local_field, err, r, z)
      if (err) then
        if (global_com%exit_on_error) call err_exit
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
      if (global_com%exit_on_error) call err_exit
      if (present(err_flag)) err_flag = .true.
      return
    end select
    
    if (ele%key == sbend$ .and. .not. mode%grid%curved_coords) call restore_curvilinear()

  enddo

!----------------------------------------------------------------------------
! Unknown field_calc

case default
  call out_io (s_fatal$, r_name, 'BAD FIELD_CALC METHOD FOR ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
  if (present(err_flag)) err_flag = .true.
  return
end select

call convert_fields_to_lab_coords

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! Convert fields to lab coords

contains

subroutine convert_fields_to_lab_coords()

real(rp) w_mat(3,3), w_inv(3,3), w_s(3,3), w_rt(3,3), w_rt_inv(3,3)
real(rp) theta

!

if (local_ref_frame) return

!

if (ele%key == sbend$) then
  call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(roll$), w_mat)
  theta = ele%value(g$) * s_rel - ele%value(angle$)/2
  call w_mat_for_x_pitch (theta, w_s)
  if (ele%value(ref_tilt_tot$) == 0) then
    w_mat = matmul(matmul(w_s, w_mat), transpose(w_s))
  else
    call w_mat_for_tilt (ele%value(ref_tilt_tot$), w_rt, w_rt_inv)
    w_mat = matmul(matmul(matmul(matmul(matmul(w_rt, w_s), w_rt_inv), w_mat), w_rt), transpose(w_s))
  endif
else
  call floor_angles_to_w_mat (ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$), w_mat)
endif

field%B = matmul(w_mat, field%B)
field%E = matmul(w_mat, field%E)

if (present(potential)) then
  potential%A = matmul(w_mat, potential%A)
endif

if (df_calc) then
  w_inv = transpose(w_mat)
  field%dB = matmul(w_mat, matmul(field%dB, w_inv))
  field%dE = matmul(w_mat, matmul(field%dE, w_inv))
endif

end subroutine convert_fields_to_lab_coords

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! convert_curvilinear_to_cartesian()
!
! For sbend with Grid calculation.

subroutine convert_curvilinear_to_cartesian()
real(rp) :: temp

if (ele%value(g$) == 0) return

cos_ang = cos( (s_rel-s0)*ele%value(g$) )
sin_ang = sin( (s_rel-s0)*ele%value(g$) )

! Save values because modes may have different anchor points
x_save = x
x = (x_save + ele%value(rho$) )*cos_ang - ele%value(rho$)
z = (x_save + ele%value(rho$) )*sin_ang 

! Rotate current field into this frame
! Note that we are rotating the zx plane 
temp       = field%e(3)*cos_ang + field%e(1)*sin_ang
field%e(1) = field%e(3)*sin_ang - field%e(1)*cos_ang
field%e(3) = temp
temp       = field%b(3)*cos_ang + field%b(1)*sin_ang
field%b(1) = field%b(3)*sin_ang - field%b(1)*cos_ang
field%b(3) = temp 

end subroutine convert_curvilinear_to_cartesian

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! restore_curvilinear()
!
subroutine restore_curvilinear()
real(rp) :: temp
!For sbend with Grid calculation Restores x and s_rel, and rotates output fields.
if (ele%value(g$) == 0) return
x = x_save
temp       = field%e(3)*cos_ang - field%e(1)*sin_ang
field%e(1) = field%e(3)*sin_ang + field%e(1)*cos_ang
field%e(3) = temp
temp       = field%b(3)*cos_ang - field%b(1)*sin_ang
field%b(1) = field%b(3)*sin_ang + field%b(1)*cos_ang
field%b(3) = temp 
end subroutine restore_curvilinear


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
  if (global_com%exit_on_error) call err_exit
  err_flag = .true.
  return
end select

!-------------------------------------------------------------------------------------
contains

subroutine get_this_index (x, ix_x, i0, rel_x0, out_of_bounds, err_flag)

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

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Function field_interpolate_3d (position, field_mesh, deltas, position0) result (field)
!
! Function to interpolate a 3d field.
! The interpolation is such that the derivative is continuous.
!
! Note: For "interpolation" outside of the region covered by the field_mesh
! it is assumed that the field is constant, Equal to the field at the
! boundary.
!
! Module needed:
!   use em_field_mod
!
! Input:
!   position(3)       -- Real(rp): (x, y, z) position.
!   field_mesh(:,:,:) -- Real(rp): Grid of field points.
!   deltas(3)         -- Real(rp): (dx, dy, dz) distances between mesh points.
!   position0(3)      -- Real(rp), optional:  position at (ix0, iy0, iz0) where
!                            (ix0, iy0, iz0) is the lower bound of the
!                            filed_mesh(i, j, k) array. If not present then
!                            position0 is taken to be (0.0, 0.0, 0.0)
! Output:
!   field -- Real(rp): interpolated field.
!-

function field_interpolate_3d (position, field_mesh, deltas, position0) result (field)

real(rp), optional, intent(in) :: position0(3)
real(rp), intent(in) :: position(3), field_mesh(0:,0:,0:), deltas(3)
real(rp) field

real(rp) r(3), f(-1:2), g(-1:2), h(-1:2), r_frac(3)

integer i0(3), ix, iy, iz, iix, iiy, iiz

!

if (present(position0)) then
  r = (position - position0) / deltas
else
  r = position / deltas
endif

i0 = int(r)
r_frac = r - i0

do ix = -1, 2
 iix = min(max(ix + i0(1), 0), ubound(field_mesh, 1))
 do iy = -1, 2
    iiy = min(max(iy + i0(2), 0), ubound(field_mesh, 2))
    do iz = -1, 2
      iiz = min(max(iz + i0(3), 0), ubound(field_mesh, 3))
      f(iz) = field_mesh(iix, iiy, iiz)
    enddo
    g(iy) = interpolate_1d (r_frac(3), f)
  enddo
  h(ix) = interpolate_1d (r_frac(2), g)
enddo
field = interpolate_1d (r_frac(1), h)

!---------------------------------------------------------------

contains

! interpolation in 1 dimension using 4 equally spaced points: P1, P2, P3, P4.
!   x = interpolation point.
!           x = 0 -> point is at P2.
!           x = 1 -> point is at P3.
! Interpolation is done so that the derivative is continuous.
! The interpolation uses a cubic polynomial

function interpolate_1d (x, field1_in) result (field1)

real(rp) field1, x, field1_in(4), df_2, df_3
real(rp) c0, c1, c2, c3

!

df_2 = (field1_in(3) - field1_in(1)) / 2   ! derivative at P2
df_3 = (field1_in(4) - field1_in(2)) / 2   ! derivative at P3

c0 = field1_in(2)
c1 = df_2
c2 = 3 * field1_in(3) - df_3 - 3 * field1_in(2) - 2 * df_2
c3 = df_3 - 2 * field1_in(3) + 2 * field1_in(2) + df_2

field1 = c0 + c1 * x + c2 * x**2 + c3 * x**3

end function interpolate_1d

end function field_interpolate_3d 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function e_accel_field (ele, voltage_or_gradient) result (field)
!
! Routine to return the gradient or voltage through an e_gun, lcavity or rfcavity element.
!
! Moudle needed:
!   use track1_mod
!
! Input:
!   ele                 -- ele_struct: Lcavity or rfcavity element.
!   voltage_or_gradient -- integer: voltage$ or gradient$
!
! Output:
!   gradient -- real(rp): cavity gradient
!-

function e_accel_field (ele, voltage_or_gradient) result (field)

type (ele_struct) ele
real(rp) field
integer voltage_or_gradient 

!

if (.not. ele%is_on) then
  field = 0
  return
endif

select case (ele%key)
case (lcavity$)
  select case (voltage_or_gradient)
  case (voltage$)
    field = (ele%value(gradient$) + ele%value(gradient_err$)) * ele%value(field_factor$) * ele%value(l$)
  case (gradient$)
    field = (ele%value(gradient$) + ele%value(gradient_err$)) * ele%value(field_factor$)
  end select

case (rfcavity$)
  select case (voltage_or_gradient)
  case (voltage$)
    field = ele%value(voltage$) * ele%value(field_factor$)
  case (gradient$)
    field = ele%value(voltage$) * ele%value(field_factor$) / ele%value(l$)
  end select

case (e_gun$)
  select case (voltage_or_gradient)
  case (voltage$)
    field = ele%value(gradient$) * ele%value(field_factor$) * ele%value(l$)
  case (gradient$)
    field = ele%value(gradient$) * ele%value(field_factor$) 
  end select

end select

end function e_accel_field

end module
