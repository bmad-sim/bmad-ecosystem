!+
! Subroutine em_field_calc (ele, param, s_pos, orbit, local_ref_frame, field, calc_dfield, err_flag,
!               calc_potential, use_overlap, grid_allow_s_out_of_bounds, rf_time, used_eles, print_err)
!
! Routine to calculate the E and B fields at a particular place in an element.
!
! Note: Zero field will be returned if an element is turned off.
!
! Note: The fields due to any kicks will be present. It therefore important in tracking to make sure that 
! offset_particle does not add in kicks at the beginning and end which would result in double counting the kicks.
!
! Input:
!   ele             -- Ele_struct: Lattice element.
!   param           -- lat_param_struct: Lattice parameters.
!   s_pos           -- Real(rp): Longitudinal position.
!                        If local_ref_frame = T: In Body coords relative to the entrance edge of the element.
!                        If local_ref_frame = F: In Lab coords relative to the upstream edge of the element.
!   orbit           -- Coord_struct: Transverse coordinates.
!     %vec(1), %vec(3) -- Transverse coords.
!     %t               -- Used with absolute time tracking.
!     %vec(5)          -- Used with relative time tracking (except with time Runge-Kutta).
!   local_ref_frame  -- Logical, If True then take the input coordinates and output fields 
!                         as being with respect to the frame of referene of the element (ignore misalignments). 
!   calc_dfield      -- Logical, optional: If present and True then calculate the field derivatives.
!   use_overlap      -- logical, optional: Add in overlap fields from other elements? Default is True.
!   calc_potential   -- logical, optional: Calc electric and magnetic potentials? Default is false. 
!                         This is experimental and only implemented for wigglers at present.
!   grid_allow_s_out_of_bounds 
!                    -- logical, optional: For grids, allow s-coordinate to be grossly out of bounds 
!                         and return zero instead of an error? Default: False. Used internally for overlapping fields.
!   rf_time          -- real(rp), optional: Set the time relative to the RF clock. Normally this time is calculated using
!                         orbit%t or orbit%vec(5) but sometimes it is convenient to be able to override this.
!                         For example, time_runge_kutta uses this.
!   used_eles(:)     -- ele_pointer_struct, allocatable, optional: For internal use only when this routine is
!                         called recursively. Used to prevent double counting when there is field overlap.
!   print_err        -- logical, optional: Print an error message? Default is True.
!                         For example, if the particle is out of bounds when the field is defined on a grid.
!
! Output:
!   field       -- em_field_struct: E and B fields and derivatives.
!   err_flag    -- logical, optional: Set True if there is an error. False otherwise.
!-

recursive subroutine em_field_calc (ele, param, s_pos, orbit, local_ref_frame, field, calc_dfield, err_flag, &
             calc_potential, use_overlap, grid_allow_s_out_of_bounds, rf_time, used_eles, print_err)

use super_recipes_mod
use em_field_mod, dummy => em_field_calc
use multipole_mod, dummy2 => em_field_calc

implicit none

type (ele_struct), target :: ele, ele2
type (ele_pointer_struct), allocatable, optional :: used_eles(:)
type (ele_pointer_struct), allocatable :: used_list(:)
type (ele_struct), pointer :: lord
type (lat_param_struct) param
type (coord_struct) :: orbit, local_orb, lab_orb, lord_orb, this_orb
type (em_field_struct) :: field, field1, field2, lord_field, l1_field, mode_field
type (cartesian_map_struct), pointer :: ct_map
type (cartesian_map_term1_struct), pointer :: ct_term
type (cylindrical_map_struct), pointer :: cl_map
type (cylindrical_map_term1_struct), pointer :: cl_term
type (gen_grad_map_struct), pointer :: gg_map
type (grid_field_struct), pointer :: g_field, g_field_ptr
type (grid_field_pt1_struct) g_pt
type (floor_position_struct) lab_position, global_position, lord_position
type (spline_struct) spline
type (branch_struct), pointer :: branch

real(rp), optional :: rf_time
real(rp) :: x, y, j1, dj1, time, s_pos, s_body, s_lab, s_lab2, z, ff, dk(3,3), ref_charge, f_p0c
real(rp) :: c_x, s_x, c_y, s_y, c_z, s_z, ch_x, ch_y, sh_x, sh_y, coef, fd(3), Ex, Ey, amp
real(rp) :: cos_ang, sin_ang, sgn_x, sgn_y, sgn_z, dkm(2,2), cos_ks, sin_ks, length
real(rp) phase, gradient, r, E_r_coef, E_s, k_wave, s_eff, a_amp, inte
real(rp) k_t, k_zn, kappa2_n, kap_rho, s_active_offset, beta_start, f, f1, f2, f3, kx, ky, kz
real(rp) radius, phi, t_ref, tilt, omega, freq0, freq, B_phi_coef, z_center
real(rp) sx_over_kx, sy_over_ky, sz_over_kz, rot2(2,2)
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx), pot
real(rp) w_ele_mat(3,3), w_lord_mat(3,3), Er, Ep, Ez, Br, Bp, Bx, By, Bz
real(rp) :: fld(3), dfld(3,3), fld0(3), fld1(3), dfld0(3,3), dfld1(3,3)
real(rp) phi0_autoscale, field_autoscale, ds, beta_ref, ds_small, abs_tol
real(rp) rho, a, b, B0, gamma, Brho, voltage, k_rf
real(rp) rad_p, z_p, alpha_p, beta_p, k_p, rad_m, z_m, alpha_m, beta_m, k_m

complex(rp) exp_kz, dEp, dEr, E_rho, E_phi, E_z, B_rho, B_phi, B_z
complex(rp) Im_0, Im_plus, Im_minus, Im_0_R, kappa_n, Im_plus2, cm, sm, q
complex(rp), target :: expt
complex(rp), pointer :: expt_ptr

integer i, j, m, n, ix, trig_x, trig_y, status, im, iz0, iz1, izp, ix_pole_max

logical :: local_ref_frame
logical, optional :: calc_dfield, calc_potential, err_flag, use_overlap, grid_allow_s_out_of_bounds, print_err
logical do_df_calc, err, dfield_computed, add_kicks

character(*), parameter :: r_name = 'em_field_calc'

! Initialize field
! If element is turned off then return zero

field = em_field_struct()

do_df_calc = logic_option (.false., calc_dfield)
dfield_computed = .false.

if (present(err_flag)) err_flag = .false.
if (.not. ele%is_on) return

! Has this element been used before? If so nothing to be done.

if (present(used_eles)) then
  do j = 1, size(used_eles)
    if (.not. associated (used_eles(j)%ele)) exit
    if (associated(used_eles(j)%ele, ele)) return
  enddo
endif

!----------------------------------------------------------------------------
! super_slave, multipass_slave, and slice_slave, have their field info stored in the associated lord elements.

if (ele%field_calc == refer_to_lords$) then
  if (.not. present(used_eles)) allocate (used_list(ele%n_lord+5))

  ! The lord of an element may have independent misalignments.
  ! So use an orbit that is not in the slave's reference frame.

  lab_orb = orbit

  if (local_ref_frame) then
    call offset_particle (ele, unset$, lab_orb, set_hvkicks = .false., s_pos = s_pos, s_out = s_lab)
  else
    s_lab = s_pos
  endif

  !

  lord_loop: do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)

    if (lord%field_calc == no_field$) cycle   ! Group, overlay and girder elements do not have fields.

    ! Multipass_lords do not have a well defined global position so take the lord position equal to the slave position.
    ! This is justified if all the slaves have the same position as they should in a realistic lattice.
    if (lord%lord_status == multipass_lord$) then
      s_lab2 = s_lab
      lord%floor = ele%floor  ! Needed if there is field overlap.
    else
      ds = ele%s_start - lord%s_start
      if (lord%value(l$) > 0 .and. lord%s_start > lord%s) then ! Element wraps around zero
        branch => pointer_to_branch(lord)
        ds = modulo2(ds, 0.5_rp*branch%param%total_length)
      endif
      s_lab2 = s_lab + ds
    endif

    if (present(used_eles)) then
      do j = 1, size(used_eles)
        if (.not. associated(used_eles(j)%ele)) exit
        if (associated(used_eles(j)%ele, lord)) cycle lord_loop
      enddo
      call em_field_calc (lord, param, s_lab2, lab_orb, .false., field2, calc_dfield, err, calc_potential, &
                            use_overlap, grid_allow_s_out_of_bounds, rf_time, used_eles, print_err)
    else
      call em_field_calc (lord, param, s_lab2, lab_orb, .false., field2, calc_dfield, err, calc_potential, &
                            use_overlap, grid_allow_s_out_of_bounds, rf_time, used_list, print_err)
    endif

    if (err) then
      if (present(err_flag)) err_flag = .true.
      return
    endif

    field%E = field%E + field2%E
    field%B = field%B + field2%B
    if (do_df_calc) then
      field%dE = field%dE + field2%dE
      field%dB = field%dB + field2%dB
    endif

  enddo lord_loop

  if (local_ref_frame) call convert_field_ele_to_lab(ele, s_lab, .false., field, calc_dfield, calc_potential)
  return
endif

!----------------------------------------------------------------------------
! Custom field calc 

if (ele%field_calc == custom$) then
  if (.not. associated(em_field_custom_ptr)) then
    call out_io (s_error$, r_name, 'FIELD_CALC == CUSTOM INVALID SINCE EM_FIELD_CUSTOM_PTR HAS NOT BEEN SET IN THIS PROGRAM!')
    orbit%state = lost$
    return
  endif
  call em_field_custom_ptr (ele, param, s_pos, orbit, local_ref_frame, field, calc_dfield, err_flag, &
                                    calc_potential, use_overlap, grid_allow_s_out_of_bounds, rf_time, used_eles)
  return
end if

!-----
! If the used_eles list is present then put the present element in the list.

if (present(used_eles)) then
  do j = 1, size(used_eles)
    if (associated (used_eles(j)%ele)) cycle
    used_eles(j)%ele => ele
    exit
  enddo

  if (j == size(used_eles) + 1) then
    call move_alloc(used_eles, used_list)
    allocate(used_eles(2*j))
    used_eles(1:j-1) = used_list
    used_eles(j)%ele => ele
  endif
endif

!----------------------------------------------------------------------------
! convert to local coords

local_orb = orbit
if (local_ref_frame) then
  s_body = s_pos
else
  call offset_particle (ele, set$, local_orb, set_hvkicks = .false., s_pos = s_pos, s_out = s_body)
endif

!----------------------------------------------------------------------------
! Set up common variables for all (non-custom) methods

ref_charge = charge_of(ele%ref_species)

x = local_orb%vec(1)
y = local_orb%vec(3)

if (ref_charge == 0) then
  f_p0c = 0
else
  f_p0c = ele%value(p0c$) / (c_light * ref_charge)
endif

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
! field_calc methods

select case (ele%field_calc)

!----------------------------------------------------------------------------------------------
! Soft_edge

case (soft_edge$)

  select case (ele%key)

  ! See:
  !   "Cylindrical Magnets and Ideal Solenoids"
  !   Norman Derby & Stanislaw Olbert
  !   https://arxiv.org/pdf/0909.3880.pdf

  case (solenoid$)
    rho = norm2([x, y])
    b = ele%value(l_soft_edge$) / 2
    a = ele%value(r_solenoid$)
    B0 = ele%value(bs_field$) / pi

    if (a == 0) then
      call out_io (s_fatal$, r_name, 'R_SOLENOID NOT SET WHEN USING "SOFT_EDGE" SOLENOID FIELD.', &
                                     'FOR ELEMENT: ' // ele%name)
      orbit%state = lost$
      if (present(err_flag)) err_flag = .true.
      return
    endif

    z = s_body - ele%value(l$) / 2
    z_p = z + b;                                  z_m = z - b
    rad_p = sqrt(z_p**2 + (rho + a)**2);          rad_m = sqrt(z_m**2 + (rho + a)**2)
    alpha_p = a / rad_p;                          alpha_m = a / rad_m
    beta_p = z_p / rad_p;                         beta_m = z_m / rad_m
    k_p = sqrt(z_p**2 + (a - rho)**2) / rad_p;    k_m = sqrt(z_m**2 + (a - rho)**2) / rad_m
    gamma = (a - rho) / (a + rho)

    if (rho /= 0) then
      Brho = B0 * (alpha_p * gen_complete_elliptic(k_p, 1.0_rp, 1.0_rp, -1.0_rp) - &
                   alpha_m * gen_complete_elliptic(k_m, 1.0_rp, 1.0_rp, -1.0_rp))
      field%B(1) = x * Brho / rho
      field%B(2) = y * Brho / rho
    endif

    field%B(3) = B0 * a * (beta_p * gen_complete_elliptic(k_p, gamma**2, 1.0_rp, gamma) - &
                           beta_m * gen_complete_elliptic(k_m, gamma**2, 1.0_rp, gamma)) / (a + rho)

  case default
    call out_io (s_fatal$, r_name, '"SOFT_EDGE" FIELD NOT YET CODED FOR ELEMENT OF TYPE: ' // key_name(ele%key), &
                                   'FOR ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    orbit%state = lost$
    if (present(err_flag)) err_flag = .true.
    return
  end select

!----------------------------------------------------------------------------------------------
! Bmad_standard field calc 

case (bmad_standard$)

  ! Field outside of element is zero. 

  ds_small = bmad_com%significant_length / 10.0_rp
  if (s_body < -ds_small .or. s_body > ele%value(l$) + ds_small) goto 8000   ! Goto field overlap code.

  select case (ele%key)

  !------------------
  ! Crab cavity

  case (crab_cavity$)

    ! The crab cavity is modeled as a TM110 traveling wave mode
    if (ele%value(l$) /= 0) then
      voltage = e_accel_field(ele, voltage$) / ref_charge
      if (present(rf_time)) then
        time = rf_time
      else
        time = particle_rf_time(orbit, ele, .false., s_body)
      endif
      phase = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_autoscale$) - &
                      (time - rf_ref_time_offset(ele) - s_body/c_light) * ele%value(rf_frequency$))

      k_rf = twopi * ele%value(rf_frequency$) / c_light
      field%B(2) = -voltage * sin(phase) / (c_light * ele%value(l$))
      field%E(3) = voltage * k_rf * orbit%beta * orbit%vec(1) * cos(phase) / ele%value(l$)
    endif

  !------------------
  ! Drift, et. al. Note that kicks get added at the end for all elements

  case (drift$, ecollimator$, rcollimator$, instrument$, monitor$, pipe$, marker$, detector$, thick_multipole$)

  !------------------
  ! E_Gun

  case (e_gun$)
    if (ele%value(rf_frequency$) == 0) then
      field%e(3) = e_accel_field (ele, gradient$) / ref_charge
    else
      if (present(rf_time)) then
        time = rf_time
      else
        time = particle_rf_time(orbit, ele, .true., s_body)
      endif
      phase = (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_err$) + ele%value(phi0_autoscale$))
      field%e(3) = e_accel_field (ele, gradient$) * cos(twopi * (time * ele%value(rf_frequency$) + phase)) / ref_charge
    endif

  !------------------
  ! Elseparator

  case (elseparator$)
    field%e(1) = ele%value(hkick$) * ele%value(p0c$) / ele%value(l$)
    field%e(2) = ele%value(vkick$) * ele%value(p0c$) / ele%value(l$)

  !------------------
  ! HKicker

  case (hkicker$)
    field%b(2) = -ele%value(kick$) * f_p0c / ele%value(l$)

  !------------------
  ! Kicker  

  case (kicker$, ac_kicker$)
    field%b(1) =  ele%value(vkick$) * f_p0c / ele%value(l$)
    field%b(2) = -ele%value(hkick$) * f_p0c / ele%value(l$)

  !------------------
  ! RFcavity and Lcavity  bmad_standard
  !
  ! For standing wave cavity:
  ! Use N_cell half-wave pillbox formulas for TM_011 mode with infinite wall radius.
  ! See S.Y. Lee, "Accelerator Physics"
  !   E_s   = 2 * gradient *         cos(k s) * cos(omega t + phase)
  !   E_r   =     gradient * k * r * sin(k s) * cos(omega t + phase)
  !   B_phi =    -gradient * k * r * cos(k s) * sin(omega t + phase) / c_light
  ! For traveling wave cavity:
  !   E_s   =  gradient * cos(omega t + phase - k s)
  !   E_r   = -gradient * k * r * sin(omega t + phase - k s) / 2
  !   B_phi = -gradient * k * r * sin(omega t + phase - k s) / c_light / 2
  ! 
  ! Note: Length of pillbox is 1/2 wavelength. Not the length of the element.  
  ! That is, the field extends to +/- c_light * freq / 2 from centerline of element.
  !
  ! Since the active (pillbox) length is different from the element length, the gradient used
  ! is different from the element-gradient = voltage / element-length stored in the element struct so that
  !   gradient-used * pillbox-length = element-gradient * element-length = voltage
  !
  ! Note: There is a discontinuity in the field at the edge. Edge focusing due to this 
  !  discontinuity can be handled in the apply_element_edge_kick routine.

  case(rfcavity$, lcavity$)

    if (ele%value(rf_frequency$) == 0) return

    phase = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_err$) + ele%value(phi0_autoscale$))
    if (ele%key == rfcavity$) phase = pi/2 - phase
    orbit%phase(1) = phase  ! RF phase is needed by apply_element_edge_kick when calling rf_coupler_kick.

    gradient = e_accel_field (ele, gradient$)

    if (.not. ele%is_on) gradient = 0
    gradient = (gradient + gradient_shift_sr_wake(ele, param)) / ref_charge
    gradient = gradient * ele%value(l$) / ele%value(l_active$)
    omega = twopi * ele%value(rf_frequency$)
    k_wave = omega / c_light

    s_active_offset = (ele%value(l$) - ele%value(l_active$)) / 2  ! Relative to entrance end of the cavity
    s_eff = s_body - s_active_offset
    if (s_eff < 0 .or. s_eff > ele%value(l_active$)) then
      dfield_computed = .true.
      goto 8000  ! Zero field outside
    endif

    beta_start = ele%value(p0c_start$) / ele%value(e_tot_start$)

    if (present(rf_time)) then
      time = rf_time
    else
      time = particle_rf_time(orbit, ele, .true., s_body)
    endif
    
    if (nint(ele%value(cavity_type$)) == traveling_wave$) then
      phi = omega * time + phase - k_wave * s_eff
      E_z        =  gradient * cos(phi)
      E_r_coef   = -gradient * k_wave * sin(phi) / 2.0_rp
      B_phi_coef = -gradient * k_wave * sin(phi) / (2.0_rp * c_light)
    elseif (nint(ele%value(longitudinal_mode$)) == 0) then
      E_z        = 2.0_rp * gradient *          sin(k_wave*s_eff) * sin(omega * time + phase)
      E_r_coef   =         -gradient * k_wave * cos(k_wave*s_eff) * sin(omega * time + phase)
      B_phi_coef =          gradient * k_wave * sin(k_wave*s_eff) * cos(omega * time + phase) / c_light 
    elseif (nint(ele%value(longitudinal_mode$)) == 1) then
      E_z        = 2.0_rp * gradient *          cos(k_wave*s_eff) * cos(omega * time + phase)
      E_r_coef   =          gradient * k_wave * sin(k_wave*s_eff) * cos(omega * time + phase)
      B_phi_coef =         -gradient * k_wave * cos(k_wave*s_eff) * sin(omega * time + phase) / c_light 
    else
      call out_io (s_error$, r_name, 'LONGITUDINAL_MODE PARAMETER IS NOT 0 NOR 1 FOR ELEMENT: ' // ele%name)
      return
    endif

    field%E(1) = E_r_coef * x
    field%E(2) = E_r_coef * y
    field%E(3) = E_z
    
    field%B(1) = -B_phi_coef * y
    field%B(2) =  B_phi_coef * x

    if (do_df_calc) then
      dfield_computed = .true.
      field%dE(1,1) =  E_r_coef
      field%dE(2,2) =  E_r_coef
      field%dB(1,2) = -B_phi_coef
      field%dB(2,1) =  B_phi_coef
      if (nint(ele%value(cavity_type$)) == traveling_wave$) then
        f = gradient * k_wave**2 * cos(phi) / 2.0_rp
        field%dE(1,3) =  x * f
        field%dE(2,3) =  y * f
        field%dE(3,3) =  k_wave * gradient * sin(phi)
        field%dB(1,3) = -y * f / c_light
        field%dB(2,3) =  x * f / c_light
      else
        if (nint(ele%value(longitudinal_mode$)) == 0) then
          f1 = gradient * k_wave**2 * sin(k_wave*s_eff) * sin(omega * time + phase)
          f2 = gradient * k_wave**2 * cos(k_wave*s_eff) * cos(omega * time + phase)
          f3 = gradient * k_wave * cos(k_wave*s_eff) * sin(omega * time + phase)
        else
          f1 =  gradient * k_wave**2 * cos(k_wave*s_eff) * cos(omega * time + phase)
          f2 =  gradient * k_wave**2 * sin(k_wave*s_eff) * sin(omega * time + phase)
          f3 = -gradient * k_wave * sin(k_wave*s_eff) * cos(omega * time + phase)
        endif
        field%dE(1,3) =  x * f1
        field%dE(2,3) =  y * f1
        field%dE(3,3) =  2.0_rp * f3
        field%dB(1,3) = -y * f2 / c_light
        field%dB(2,3) =  x * f2 / c_light
      endif
    endif

  !------------------
  ! Octupole 

  case (octupole$)

    field%b(1) = -(y**3 - 3*y*x**2) / 6 * ele%value(k3$) * f_p0c 
    field%b(2) =  (x**3 - 3*x*y**2) / 6 * ele%value(k3$) * f_p0c 

    if (do_df_calc) then
      dfield_computed = .true.
      field%dB(1,1) =  x*y * ele%value(k3$) * f_p0c
      field%dB(1,2) = (x**2 - y**2) / 2 * ele%value(k3$) * f_p0c
      field%dB(2,1) = (x**2 - y**2) / 2 * ele%value(k3$) * f_p0c
      field%dB(2,2) = -x*y * ele%value(k3$) * f_p0c
    endif

  !------------------
  ! Patch: There are no fields

  case (patch$)

  !------------------
  ! Quadrupole

  case (quadrupole$) 

    field%b(1) = y * ele%value(k1$) * f_p0c 
    field%b(2) = x * ele%value(k1$) * f_p0c 

    if (do_df_calc) then
      dfield_computed = .true.
      field%dB(1,2) =  ele%value(k1$) * f_p0c
      field%dB(2,1) =  ele%value(k1$) * f_p0c
    endif

    if (logic_option(.false., calc_potential)) then
      field%A(3) = 0.5_rp * (y * field%b(1) - x * field%b(2)) 
    endif

  !------------------
  ! Sextupole 

  case (sextupole$)

    field%b(1) = x * y * ele%value(k2$) * f_p0c
    field%b(2) = (x**2 - y**2) / 2 * ele%value(k2$) * f_p0c 

    if (do_df_calc) then
      dfield_computed = .true.
      field%dB(1,1) =  y * ele%value(k2$) * f_p0c
      field%dB(1,2) =  x * ele%value(k2$) * f_p0c
      field%dB(2,1) =  x * ele%value(k2$) * f_p0c
      field%dB(2,2) = -y * ele%value(k2$) * f_p0c
    endif

  !------------------
  ! VKicker

  case (vkicker$)
    field%b(1) =  ele%value(kick$) * f_p0c / ele%value(l$)

  !------------------
  ! SBend

  case (sbend$)

    ! Finite dg, k1 or k2 is handled with rest of multipoles
    field%b(2) = ele%value(g$) * f_p0c 


  !------------------
  ! Sol_quad

  case (sol_quad$)

    field%b(1) = y * ele%value(k1$) * f_p0c 
    field%b(2) = x * ele%value(k1$) * f_p0c 
    field%b(3) = ele%value(ks$) * f_p0c

    if (do_df_calc) then
      dfield_computed = .true.
      field%dB(1,2) = ele%value(k1$) * f_p0c
      field%dB(2,1) = ele%value(k1$) * f_p0c
    endif

    if (logic_option(.false., calc_potential)) then
      field%A = (0.5_rp * field%b(3)) * [-y, x, 0.0_rp]
      field%A(3) = 0.5_rp * (y * field%b(1) - x * field%b(2)) 
    endif

  !------------------
  ! Solenoid / sad_mult

  case (solenoid$, sad_mult$)

    field%b(3) = ele%value(ks$) * f_p0c

    if (do_df_calc) then
      dfield_computed = .true.
    endif

    if (logic_option(.false., calc_potential)) then
      if (ele%key == sad_mult$) then
        field%A = (0.5_rp * field%b(3)) * [-y+ele%value(y_offset_mult$), x-ele%value(x_offset_mult$), 0.0_rp]      
      else
        field%A = (0.5_rp * field%b(3)) * [-y, x, 0.0_rp]      
      endif
    endif

  !------------------
  ! Wiggler

  case(wiggler$, undulator$)

    ! Should not be here. 
    call out_io (s_fatal$, r_name, 'BOOKKEEPING ERROR. PLEASE GET HELP. FOR: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    if (present(err_flag)) err_flag = .true.
    return

  !------------------
  ! Error

  case default
    if (logic_option(.true., print_err)) call out_io (s_fatal$, r_name, &
                          'BMAD_STANDARD FIELD NOT YET CODED FOR ELEMENT OF TYPE: ' // key_name(ele%key), &
                          'FOR ELEMENT: ' // ele%name, 'PERHAPS "FIELD_CALC" NEEDS TO BE SET FOR THIS ELEMENT?')
    if (global_com%exit_on_error) call err_exit
    if (present(err_flag)) err_flag = .true.
    return
  end select

  !---------------------------------------------
  ! Add multipoles

  if (ele%key == sbend$ .and. nint(ele%value(exact_multipoles$)) /= off$ .and. ele%value(g$) /= 0) then
    call bend_exact_multipole_field (ele, param, orbit, local_ref_frame, field2, do_df_calc, calc_potential)
    field%e = field%e + field2%e
    field%b = field%b + field2%b
    if (do_df_calc) then
      field%de = field2%de
      field%db = field2%db
    endif
    add_kicks = .false.  ! h/vkicks are accounted for in bend_exact_multipole_field

  ! Everything but exact bend multipoles
  else
    ! First magnetic
    ! This should be cleaned up so that include_kicks is always present.
    ! Do do this, the code above which puts in the kick, dg, k1, k2, k3 kicks needs to be removed.
    if (ele%key == sbend$) then
      call multipole_ele_to_ab(ele, .false., ix_pole_max, a_pole, b_pole, magnetic$, include_kicks$)
      add_kicks = .false.
    else
      call multipole_ele_to_ab(ele, .false., ix_pole_max, a_pole, b_pole, magnetic$)
      add_kicks = .true.
    endif

    if (ix_pole_max > -1) then

      if (ele%value(l$) == 0) then
        call out_io (s_fatal$, r_name, 'CANNOT COMPUTE FIELD OF ZERO LENGTH ELEMENT WITH MULTIPOLES. FOR: ' // ele%name)
        if (global_com%exit_on_error) call err_exit
        if (present(err_flag)) err_flag = .true.
        return
      endif

      ! Note: The spin_integration routine assumes that the field is zero outside of the active region.
      select case (ele%key)
      case (lcavity$, rfcavity$, crab_cavity$)
        length = ele%value(l_active$)
      case default
        length = ele%value(l$)
      end select

      do i = 0, ix_pole_max
        if (a_pole(i) == 0 .and. b_pole(i) == 0) cycle
        if (do_df_calc) then
          call ab_multipole_kick(a_pole(i), b_pole(i), i, ele%ref_species, 0, local_orb, kx, ky, dkm)
        else
          call ab_multipole_kick(a_pole(i), b_pole(i), i, ele%ref_species, 0, local_orb, kx, ky)
        endif
        field%B(1) = field%B(1) + f_p0c * ky / length
        field%B(2) = field%B(2) - f_p0c * kx / length
        if (do_df_calc) then
          field%dB(1,1) = field%dB(1,1) + f_p0c * dkm(2,1) / length
          field%dB(1,2) = field%dB(1,2) + f_p0c * dkm(2,2) / length
          field%dB(2,1) = field%dB(2,1) - f_p0c * dkm(1,1) / length
          field%dB(2,2) = field%dB(2,2) - f_p0c * dkm(1,2) / length
        endif
      enddo
    endif

    ! Add electric multipoles

    call multipole_ele_to_ab(ele, .not. local_ref_frame, ix_pole_max, a_pole, b_pole, electric$)
    do i = 0, ix_pole_max
      if (a_pole(i) == 0 .and. b_pole(i) == 0) cycle
      call elec_multipole_field(a_pole(i), b_pole(i), i, local_orb, Ex, Ey, dkm, do_df_calc)
      field%E(1) = field%E(1) + Ex
      field%E(2) = field%E(2) + Ey
      if (do_df_calc) field%dE(1:2,1:2) = field%dE(1:2,1:2) + dkm
    enddo

  endif

  !-------------------------------
  ! Add kicks. Since the kick field is not rotated by a tilt then we have to unrotate if in the local_ref_frame

  if (add_kicks .and. has_hkick_attributes(ele%key) .and. (ele%value(hkick$) /= 0 .or. ele%value(vkick$) /= 0)) then
    select case (ele%key)
    ! Kickers and elsep handled above
    case (ac_kicker$, kicker$, hkicker$, vkicker$, elseparator$)  
    ! Everything else...
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

!----------------------------------------------------------------------------------------------
! planar_model

case(planar_model$)

  kz = twopi / ele%value(l_period$)
  kx = ele%value(kx$)
  ky = sqrt(kx**2 + kz**2)
  c_x = cos(kx * x)
  s_x = sin(kx * x)
  ch_y = cosh(ky * y)
  sh_y = sinh(ky * y)
  c_z = cos(kz * (s_body - ele%value(l$)/2))
  s_z = sin(kz * (s_body - ele%value(l$)/2))
  bx = (kx/ky) * ele%value(b_max$) 
  by =           ele%value(b_max$) 
  bz = (kz/ky) * ele%value(b_max$)

  field%B(1) = -bx * s_x * sh_y * c_z
  field%B(2) =  by * c_x * ch_y * c_z
  field%B(3) = -bz * c_x * sh_y * s_z

  if (do_df_calc) then
    dfield_computed = .true.
    field%db(1,1) = -kx * bx * c_x * sh_y * c_z
    field%db(1,2) = -ky * bx * s_x * ch_y * c_z
    field%db(1,3) =  kz * bx * s_x * sh_y * s_z

    field%db(2,1) = -kx * by * s_x * ch_y * c_z
    field%db(2,2) =  ky * by * c_x * sh_y * c_z
    field%db(2,3) = -kz * by * c_x * ch_y * s_z

    field%db(3,1) =  kx * by * s_x * ch_y * s_z
    field%db(3,2) = -ky * by * c_x * sh_y * s_z
    field%db(3,3) = -kz * by * c_x * ch_y * c_z
  endif

!----------------------------------------------------------------------------------------------
! helical_model

case(helical_model$)

  kz = twopi * ele%value(n_period$) / ele%value(l$)
  ch_x = cosh(kz * x)
  sh_x = sinh(kz * x)
  ch_y = cosh(kz * y)
  sh_y = sinh(kz * y)
  c_z = cos(kz * (s_body - ele%value(l$)/2))
  s_z = sin(kz * (s_body - ele%value(l$)/2))

  field%B(1) = -ele%value(b_max$) * ch_x * s_z
  field%B(2) =  ele%value(b_max$) * ch_y * c_z
  field%B(3) = -ele%value(b_max$) * (sh_x * c_z + sh_y * s_z)

  if (do_df_calc) then
    dfield_computed = .true.
    field%db(1,1) = -kz * ele%value(b_max$) * sh_x * s_z
    field%db(1,3) = -kz * ele%value(b_max$) * ch_x * c_z
    field%db(2,2) =  kz * ele%value(b_max$) * sh_y * c_z
    field%db(2,3) = -kz * ele%value(b_max$) * ch_y * s_z
    field%db(3,1) = -kz * ele%value(b_max$) * ch_x * c_z
    field%db(3,2) = -kz * ele%value(b_max$) * ch_y * s_z
    field%db(3,3) =  kz * ele%value(b_max$) * (sh_x * s_z - sh_y * c_z)
  endif

!----------------------------------------------------------------------------------------------
! FieldMap

case(fieldmap$)

  if (present(rf_time)) then
    time = rf_time
  else
    time = particle_rf_time(orbit, ele, .false., s_body)
  endif

  if (.not. associated(ele%cylindrical_map) .and. .not. associated(ele%cartesian_map) .and. &
      .not. associated(ele%gen_grad_map) .and. .not. associated(ele%grid_field)) then
    call out_io (s_fatal$, r_name, 'No associated fieldmap (cartesican_map, grid_field, etc) FOR: ' // ele%name) 
    if (global_com%exit_on_error) call err_exit
    if (present(err_flag)) err_flag = .true.
    return  
  endif

  select case (ele%key)
  case (e_gun$, em_field$, lcavity$, rfcavity$)
    phi0_autoscale = ele%value(phi0_autoscale$)
    field_autoscale = ele%value(field_autoscale$)
  case default
    phi0_autoscale = 0
    field_autoscale = 1
  end select

  !------------------------------------
  ! Cartesian map field

  if (associated(ele%cartesian_map)) then
    do im = 1, size(ele%cartesian_map)
      ct_map => ele%cartesian_map(im)

      fld = 0; dfld = 0

      call to_fieldmap_coords (ele, local_orb, s_body, ct_map%ele_anchor_pt, ct_map%r0, .false., x, y, z, cos_ang, sin_ang, err)
      if (err) then
        if (present(err_flag)) err_flag = .true.
        return
      endif

      n = size(ct_map%ptr%term)
      do i = 1, n
        ct_term => ct_map%ptr%term(i)
        sgn_x = 1; sgn_y = 1; sgn_z = 1

        select case (ct_term%form)
        case (hyper_y$)
          coef = ct_term%coef / ct_term%ky
          c_x = cos(ct_term%kx * (x + ct_term%x0))
          s_x = sin(ct_term%kx * (x + ct_term%x0))
          c_y = cosh (ct_term%ky * (y + ct_term%y0))
          s_y = sinh (ct_term%ky * (y + ct_term%y0))
          select case (ct_term%family)
          case (family_y$);  sgn_x = -1
          case (family_sq$); sgn_x = -1; sgn_z = -1
          end select
          trig_x = -1; trig_y = 1

        case (hyper_xy$)
          coef = ct_term%coef / ct_term%kz
          c_x = cosh(ct_term%kx * (x + ct_term%x0))
          s_x = sinh(ct_term%kx * (x + ct_term%x0))
          c_y = cosh (ct_term%ky * (y + ct_term%y0))
          s_y = sinh (ct_term%ky * (y + ct_term%y0))
          select case (ct_term%family)
          case (family_sq$); sgn_z = -1
          end select
          trig_x = 1; trig_y = 1

        case (hyper_x$)
          coef = ct_term%coef / ct_term%kx
          c_x = cosh(ct_term%kx * (x + ct_term%x0))
          s_x = sinh(ct_term%kx * (x + ct_term%x0))
          c_y = cos (ct_term%ky * (y + ct_term%y0))
          s_y = sin (ct_term%ky * (y + ct_term%y0))
          select case (ct_term%family)
          case (family_x$);  sgn_y = -1
          case (family_sq$); sgn_x = -1
          end select
          trig_x = 1; trig_y = -1
        end select

        c_z = cos (ct_term%kz * z + ct_term%phi_z)
        s_z = sin (ct_term%kz * z + ct_term%phi_z)

        select case (ct_term%family)
        case (family_x$)
          fld(1) = fld(1) + coef  * ct_term%kx * c_x * c_y * c_z
          fld(2) = fld(2) + coef  * ct_term%ky * s_x * s_y * c_z * sgn_y
          fld(3) = fld(3) - coef  * ct_term%kz * s_x * c_y * s_z
        case (family_y$)
          fld(1) = fld(1) + coef  * ct_term%kx * s_x * s_y * c_z * sgn_x
          fld(2) = fld(2) + coef  * ct_term%ky * c_x * c_y * c_z
          fld(3) = fld(3) - coef  * ct_term%kz * c_x * s_y * s_z
        case (family_qu$)
          fld(1) = fld(1) + coef  * ct_term%kx * c_x * s_y * c_z
          fld(2) = fld(2) + coef  * ct_term%ky * s_x * c_y * c_z
          fld(3) = fld(3) - coef  * ct_term%kz * s_x * s_y * s_z
        case (family_sq$)
          fld(1) = fld(1) + coef  * ct_term%kx * s_x * c_y * c_z * sgn_x
          fld(2) = fld(2) + coef  * ct_term%ky * c_x * s_y * c_z
          fld(3) = fld(3) + coef  * ct_term%kz * c_x * c_y * s_z * sgn_z
        end select

        if (do_df_calc) then
          dfield_computed = .true.
          select case (ct_term%family)
          case (family_x$)
            ff = coef * ct_term%kx
            dfld(1,1) = dfld(1,1) + ff * ct_term%kx * s_x * c_y * c_z * trig_x
            dfld(2,1) = dfld(2,1) + ff * ct_term%ky * c_x * s_y * c_z * sgn_y
            dfld(3,1) = dfld(3,1) - ff * ct_term%kz * c_x * c_y * s_z 
            ff = coef * ct_term%ky
            dfld(1,2) = dfld(1,2) + ff * ct_term%kx * c_x * s_y * c_z * trig_y
            dfld(2,2) = dfld(2,2) + ff * ct_term%ky * s_x * c_y * c_z * sgn_y
            dfld(3,2) = dfld(3,2) - ff * ct_term%kz * s_x * s_y * s_z * trig_y
            ff = coef * ct_term%kz
            dfld(1,3) = dfld(1,3) - ff * ct_term%kx * c_x * c_y * s_z
            dfld(2,3) = dfld(2,3) - ff * ct_term%ky * s_x * s_y * s_z * sgn_y 
            dfld(3,3) = dfld(3,3) - ff * ct_term%kz * s_x * c_y * c_z
          case (family_y$)
            ff = coef * ct_term%kx
            dfld(1,1) = dfld(1,1) + ff * ct_term%kx * c_x * s_y * c_z * sgn_x
            dfld(2,1) = dfld(2,1) + ff * ct_term%ky * s_x * c_y * c_z * trig_x
            dfld(3,1) = dfld(3,1) - ff * ct_term%kz * s_x * s_y * s_z * trig_x
            ff = coef * ct_term%ky
            dfld(1,2) = dfld(1,2) + ff * ct_term%kx * s_x * c_y * c_z * sgn_x
            dfld(2,2) = dfld(2,2) + ff * ct_term%ky * c_x * s_y * c_z * trig_y
            dfld(3,2) = dfld(3,2) - ff * ct_term%kz * c_x * c_y * s_z 
            ff = coef * ct_term%kz
            dfld(1,3) = dfld(1,3) - ff * ct_term%kx * s_x * s_y * s_z * sgn_x
            dfld(2,3) = dfld(2,3) - ff * ct_term%ky * c_x * c_y * s_z 
            dfld(3,3) = dfld(3,3) - ff * ct_term%kz * c_x * s_y * c_z 
          case (family_qu$)
            ff = coef * ct_term%kx
            dfld(1,1) = dfld(1,1) + ff * ct_term%kx * s_x * s_y * c_z * trig_x
            dfld(2,1) = dfld(2,1) + ff * ct_term%ky * c_x * c_y * c_z 
            dfld(3,1) = dfld(3,1) - ff * ct_term%kz * c_x * s_y * s_z
            ff = coef * ct_term%ky
            dfld(1,2) = dfld(1,2) + ff * ct_term%kx * c_x * c_y * c_z 
            dfld(2,2) = dfld(2,2) + ff * ct_term%ky * s_x * s_y * c_z * trig_y
            dfld(3,2) = dfld(3,2) - ff * ct_term%kz * s_x * c_y * s_z 
            ff = coef * ct_term%kz
            dfld(1,3) = dfld(1,3) - ff * ct_term%kx * c_x * s_y * s_z
            dfld(2,3) = dfld(2,3) - ff * ct_term%ky * s_x * c_y * s_z 
            dfld(3,3) = dfld(3,3) - ff * ct_term%kz * s_x * s_y * c_z 
          case (family_sq$)
            ff = coef * ct_term%kx
            dfld(1,1) = dfld(1,1) + ff * ct_term%kx * c_x * c_y * c_z * sgn_x
            dfld(2,1) = dfld(2,1) + ff * ct_term%ky * s_x * s_y * c_z * trig_x
            dfld(3,1) = dfld(3,1) + ff * ct_term%kz * s_x * c_y * s_z * sgn_z * trig_x
            ff = coef * ct_term%ky
            dfld(1,2) = dfld(1,2) + ff * ct_term%kx * s_x * s_y * c_z * sgn_x * trig_y
            dfld(2,2) = dfld(2,2) + ff * ct_term%ky * c_x * c_y * c_z 
            dfld(3,2) = dfld(3,2) + ff * ct_term%kz * c_x * s_y * s_z * sgn_z * trig_y
            ff = coef * ct_term%kz
            dfld(1,3) = dfld(1,3) - ff * ct_term%kx * s_x * c_y * s_z * sgn_x
            dfld(2,3) = dfld(2,3) - ff * ct_term%ky * c_x * s_y * s_z 
            dfld(3,3) = dfld(3,3) + ff * ct_term%kz * c_x * c_y * c_z * sgn_z
          end select
        endif

        if (logic_option(.false., calc_potential)) then
          coef = ct_term%coef 
          select case (ct_term%family)
          case (family_x$)
            if (abs(ct_term%ky * (y + ct_term%y0)) < 1d-10) then
              sy_over_ky = y + ct_term%y0
            else
              sy_over_ky = s_y / ct_term%ky
            endif
          case (family_y$)
            if (abs(ct_term%kx * (x + ct_term%x0)) < 1d-10) then
              sx_over_kx = x + ct_term%x0
            else
              sx_over_kx = s_x / ct_term%kx
            endif
          case default
            if (abs(ct_term%kz * z + ct_term%phi_z) < 1d-10) then
              if (ct_term%kz == 0) then
                sz_over_kz = 0
              else
                sz_over_kz = z + ct_term%phi_z / ct_term%kz
              endif
            else
              sz_over_kz = s_z / ct_term%kz
            endif
          end select

          if (ct_map%field_type == magnetic$) then
            select case (ct_term%family)
            case (family_x$)
              select case (ct_term%form)
              case (hyper_y$)
                field%a(1) = field%a(1) + coef * s_x * sy_over_ky * s_z * ct_term%kz / ct_term%ky
                field%a(3) = field%a(3) + coef * c_x * sy_over_ky * c_z * ct_term%kx / ct_term%ky
              case (hyper_xy$)
                field%a(1) = field%a(1) + coef * s_x * sy_over_ky * s_z
                field%a(3) = field%a(3) + coef * c_x * sy_over_ky * c_z * ct_term%kx / ct_term%kz
              case (hyper_x$)
                field%a(1) = field%a(1) + coef * s_x * sy_over_ky * s_z * ct_term%kz / ct_term%kx
                field%a(3) = field%a(3) + coef * c_x * sy_over_ky * c_z
              end select

            case (family_y$)
              select case (ct_term%form)
              case (hyper_y$)
                field%a(2) = field%a(2) - coef * sx_over_kx * s_y * s_z * ct_term%kz / ct_term%ky
                field%a(3) = field%a(3) - coef * sx_over_kx * c_y * c_z
              case (hyper_xy$)
                field%a(2) = field%a(2) - coef * sx_over_kx * s_y * s_z
                field%a(3) = field%a(3) - coef * sx_over_kx * c_y * c_z * ct_term%ky / ct_term%kz
              case (hyper_x$)
                field%a(2) = field%a(2) - coef * sx_over_kx * s_y * s_z * ct_term%kz / ct_term%kx
                field%a(3) = field%a(3) - coef * sx_over_kx * c_y * c_z * ct_term%ky / ct_term%kx
              end select

            case (family_qu$)
              select case (ct_term%form)
              case (hyper_y$)
                field%a(1) = field%a(1) + coef * s_x * c_y * sz_over_kz
                field%a(2) = field%a(2) - coef * c_x * s_y * sz_over_kz * ct_term%kx / ct_term%ky
              case (hyper_xy$)
                field%a(1) = field%a(1) + coef * s_x * c_y * sz_over_kz * ct_term%ky / ct_term%kz
                field%a(2) = field%a(2) - coef * c_x * s_y * sz_over_kz * ct_term%kx / ct_term%kz
              case (hyper_x$)
                field%a(1) = field%a(1) + coef * s_x * c_y * sz_over_kz * ct_term%ky / ct_term%kx
                field%a(2) = field%a(2) - coef * c_x * s_y * sz_over_kz
              end select

            case (family_sq$)
              select case (ct_term%form)
              case (hyper_y$)
                field%a(1) = field%a(1) + coef * c_x * s_y * sz_over_kz
                field%a(2) = field%a(2) + coef * s_x * c_y * sz_over_kz * ct_term%kx / ct_term%ky
              case (hyper_xy$)
                field%a(1) = field%a(1) + coef * c_x * s_y * sz_over_kz * ct_term%ky / ct_term%kz
                field%a(2) = field%a(2) - coef * s_x * c_y * sz_over_kz * ct_term%kx / ct_term%kz
              case (hyper_x$)
                field%a(1) = field%a(1) + coef * c_x * s_y * sz_over_kz * ct_term%ky / ct_term%kx
                field%a(2) = field%a(2) + coef * s_x * c_y * sz_over_kz
              end select
            end select

          else  ! electric$
            select case (ct_term%family)
            case (family_x$)
              select case (ct_term%form)
              case (hyper_y$);   field%phi = field%phi + coef * s_x * c_y * c_z / ct_term%ky
              case (hyper_xy$);  field%phi = field%phi + coef * s_y * c_y * c_z / ct_term%kz
              case (hyper_x$);   field%phi = field%phi + coef * s_x * c_y * c_z / ct_term%kx
              end select

            case (family_y$)
              select case (ct_term%form)
              case (hyper_y$);   field%phi = field%phi + coef * c_x * s_y * c_z / ct_term%ky
              case (hyper_xy$);  field%phi = field%phi + coef * c_x * s_y * c_z / ct_term%kz
              case (hyper_x$);   field%phi = field%phi + coef * c_x * s_y * c_z / ct_term%kx
              end select

            case (family_qu$)
              select case (ct_term%form)
              case (hyper_y$);   field%phi = field%phi + coef * s_x * s_y * c_z / ct_term%ky
              case (hyper_xy$);  field%phi = field%phi + coef * s_x * s_y * c_z / ct_term%kz
              case (hyper_x$);   field%phi = field%phi + coef * s_x * s_y * c_z / ct_term%kx
              end select

            case (family_sq$)
              select case (ct_term%form)
              case (hyper_y$);   field%phi = field%phi + coef * c_x * c_y * c_z / ct_term%ky
              case (hyper_xy$);  field%phi = field%phi + coef * c_x * c_y * c_z / ct_term%kz
              case (hyper_x$);   field%phi = field%phi - coef * c_x * c_y * c_z / ct_term%kx
              end select
            end select
          endif

        endif
      enddo

      !

      fld = fld * ct_map%field_scale * master_parameter_value(ct_map%master_parameter, ele)
      if (ele%key == sbend$ .or. ele%key == rf_bend$) call restore_curvilinear_field(fld)

      select case (ct_map%field_type)
      case (electric$)
        field%E = field%E + fld
      case (magnetic$)
        field%B = field%B + fld
      case default
        if (global_com%exit_on_error) call err_exit
      end select

      if (do_df_calc) then
        dfld = dfld * ct_map%field_scale * master_parameter_value(ct_map%master_parameter, ele)
        if ((ele%key == sbend$ .or. ele%key == rf_bend$) .and. ele%value(g$) /= 0) then
          rot2(1,:) = [ cos_ang, sin_ang]
          rot2(2,:) = [-sin_ang, cos_ang]
          dfld(1:3:2,1:3:2) = matmul(dfld(1:3:2,1:3:2), rot2)
          rot2(1,2) = -sin_ang
          rot2(2,1) =  sin_ang
          dfld(1:3:2,1:3:2) = matmul(rot2, dfld(1:3:2,1:3:2))
        endif

        select case (ct_map%field_type)
        case (electric$)
          field%dE = field%dE + dfld
        case (magnetic$)
          field%dB = field%dB + dfld
        end select
      endif

    enddo
  endif

  !------------------------------------
  ! Cylindrical map field

  if (associated(ele%cylindrical_map)) then

    do i = 1, size(ele%cylindrical_map)
      cl_map => ele%cylindrical_map(i)

      if (cl_map%harmonic /= 0) then
        freq0 = ele%value(rf_frequency$)
        freq = ele%value(rf_frequency$) * cl_map%harmonic

        if (freq0 == 0) then
          call out_io (s_fatal$, r_name, 'Element frequency is zero but cylindrical_map harmonic is not in: ' // ele%name)
          if (global_com%exit_on_error) call err_exit
          if (present(err_flag)) err_flag = .true.
          return  
        endif
        t_ref = (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_err$) + &
                                             phi0_autoscale + cl_map%phi0_fieldmap) / freq0
        if (ele%key == rfcavity$) t_ref = 0.25/freq0 - t_ref
      endif

      coef = field_autoscale * cl_map%field_scale * master_parameter_value(cl_map%master_parameter, ele)

      !

      m = cl_map%m

      if (cl_map%harmonic /= 0) k_t = twopi * freq / c_light

      call to_fieldmap_coords (ele, local_orb, s_body, cl_map%ele_anchor_pt, cl_map%r0, .false., x, y, z, cos_ang, sin_ang, err)
      if (err) then
        if (present(err_flag)) err_flag = .true.
        return
      endif

      radius = sqrt(x**2 + y**2)
      phi = atan2(y, x)

      E_rho = 0; E_phi = 0; E_z = 0
      B_rho = 0; B_phi = 0; B_z = 0

      do n = 1, size(cl_map%ptr%term)

        cl_term => cl_map%ptr%term(n)
        if (cl_term%e_coef == 0 .and. cl_term%b_coef == 0) cycle

        k_zn = twopi * (n - 1) / (size(cl_map%ptr%term) * cl_map%dz)
        if (n > 1 .and. 2 * n > size(cl_map%ptr%term)) k_zn = k_zn - twopi / cl_map%dz

        cos_ks = cos(k_zn * z)
        sin_ks = sin(k_zn * z)
        exp_kz = cmplx(cos_ks, sin_ks, rp)

        ! DC
        if (cl_map%harmonic == 0) then

          kap_rho = k_zn * radius
          if (m == 0) then
            Im_0    = I_bessel(0, kap_rho)
            Im_plus = I_bessel(1, kap_rho)
            E_rho = E_rho + real(cl_term%e_coef * exp_kz * Im_plus)
            E_z   = E_z   + real(cl_term%e_coef * exp_kz * Im_0 * i_imag)
            B_rho = B_rho + real(cl_term%b_coef * exp_kz * Im_plus)
            B_z   = B_z   + real(cl_term%b_coef * exp_kz * Im_0 * i_imag)
          else
            cm = exp_kz * cos(m * phi - cl_map%theta0_azimuth)
            sm = exp_kz * sin(m * phi - cl_map%theta0_azimuth)
            Im_plus  = I_bessel(m+1, kap_rho)
            Im_minus = I_bessel(m-1, kap_rho)
            Im_0     = kap_rho * (Im_minus - Im_plus) / (2 * m)

            q = cm * (Im_minus + Im_plus) / 2
            E_rho = E_rho + real(cl_term%e_coef * q)
            B_rho = B_rho + real(cl_term%b_coef * q)

            q = -sm * (Im_minus - Im_plus) / 2
            E_phi = E_phi + real(cl_term%e_coef * q)
            B_phi = B_phi + real(cl_term%b_coef * q)

            q = i_imag * cm * Im_0
            E_z = E_z + real(cl_term%e_coef * q)
            B_z = B_z + real(cl_term%b_coef * q)
          endif

          if (logic_option(.false., calc_potential)) then
            if (k_zn == 0) then
              if (m == 0) then
                field%phi = field%phi - coef * real(cl_term%e_coef * z * i_imag)
              elseif (m == 1) then
                field%phi = field%phi - coef * real(cl_term%e_coef * cm * radius / 2)
              endif
            elseif (m == 0) then
              field%phi = field%phi - coef * real(cl_term%e_coef * exp_kz * Im_0 / k_zn)
            else
              field%phi = field%phi - coef * real(cl_term%e_coef * cm * Im_0 / k_zn)
            endif
          endif

        ! RF mode 
        else
          kappa2_n = k_zn**2 - k_t**2
          kappa_n = sqrt(abs(kappa2_n))
          kap_rho = kappa_n * radius
          if (kappa2_n < 0) then
            kappa_n = -i_imag * kappa_n
            kap_rho = -kap_rho
          endif

          if (m == 0) then
            Im_0    = I_bessel_extended(0, kap_rho)
            Im_plus = I_bessel_extended(1, kap_rho) / kappa_n

            E_rho = E_rho - cl_term%e_coef * Im_plus * exp_kz * I_imaginary * k_zn
            E_phi = E_phi + cl_term%b_coef * Im_plus * exp_kz
            E_z   = E_z   + cl_term%e_coef * Im_0    * exp_kz

            B_rho = B_rho - cl_term%b_coef * Im_plus * exp_kz * k_zn
            B_phi = B_phi - cl_term%e_coef * Im_plus * exp_kz * k_t**2 * I_imaginary
            B_z   = B_z   - cl_term%b_coef * Im_0    * exp_kz * I_imaginary

          else
            cm = exp_kz * cos(m * phi - cl_map%theta0_azimuth)
            sm = exp_kz * sin(m * phi - cl_map%theta0_azimuth)
            Im_plus  = I_bessel_extended(m+1, kap_rho) / kappa_n**(m+1)
            Im_minus = I_bessel_extended(m-1, kap_rho) / kappa_n**(m-1)

            ! Reason for computing Im_0_R like this is to avoid divide by zero when radius = 0.
            Im_0_R  = (Im_minus - Im_plus * kappa_n**2) / (2 * m) ! = Im_0 / radius
            Im_0    = radius * Im_0_R       

            E_rho = E_rho - i_imag * (k_zn * cl_term%e_coef * Im_plus + cl_term%b_coef * Im_0_R) * cm
            E_phi = E_phi - i_imag * (k_zn * cl_term%e_coef * Im_plus + cl_term%b_coef * (Im_0_R - Im_minus / m)) * sm
            E_z   = E_z +                         cl_term%e_coef * Im_0 * cm
     
            B_rho = B_rho + i_imag * sm * (cl_term%e_coef * (m * Im_0_R + k_zn**2 * Im_plus) + &
                                          cl_term%b_coef * k_zn * (m * Im_0_R - Im_minus / m))
            B_phi = B_phi + i_imag * cm * (cl_term%e_coef * (Im_minus - (k_zn**2 + k_t**2) * Im_plus) / 2 - &
                                          cl_term%b_coef * k_zn * Im_0_R)
            B_z   = B_z +                 sm * (-cl_term%e_coef * k_zn * Im_0 + cl_term%b_coef * kappa2_n * Im_0 / m)

         endif
        endif ! cl_map%harmonic /= 0
          
      enddo  ! cl_map%ptr%term

      ! Notice that phi0, phi0_multipass, and phi0_err are folded into t_ref above.

      if (cl_map%harmonic /= 0) then
        expt = exp(-I_imaginary * twopi * (freq * (time + t_ref)))
        E_rho = E_rho * expt
        E_phi = E_phi * expt
        E_z   = E_z * expt

        expt = expt / (twopi * freq)
        B_rho = B_rho * expt
        B_phi = B_phi * expt
        B_z   = B_z * expt
      endif

      Er = real(E_rho, rp); Ep = real(E_phi, rp); Ez = real(E_z, rp)
      Br = real(B_rho, rp); Bp = real(B_phi, rp); Bz = real(B_z, rp)

      mode_field%E = coef * [cos(phi) * Er - sin(phi) * Ep, sin(phi) * Er + cos(phi) * Ep, Ez]
      mode_field%B = coef * [cos(phi) * Br - sin(phi) * Bp, sin(phi) * Br + cos(phi) * Bp, Bz]

      if (ele%key == sbend$ .or. ele%key == rf_bend$) call restore_curvilinear_field(mode_field%E, mode_field%B)

      field%E = field%E + mode_field%E
      field%B = field%B + mode_field%B

      if (logic_option(.false., calc_potential)) then
      endif

    enddo
  endif

  !------------------------------------
  ! Gen grid map calc 

  if (associated(ele%gen_grad_map)) then

    do i = 1, size(ele%gen_grad_map)
      gg_map => ele%gen_grad_map(i)

      call to_fieldmap_coords (ele, local_orb, s_body, gg_map%ele_anchor_pt, gg_map%r0, gg_map%curved_ref_frame, x, y, z, cos_ang, sin_ang, err)
      if (err) then
        if (present(err_flag)) err_flag = .true.
        return
      endif

      call gen_grad_add_em_field (ele, gg_map, [x,y,z], s_body, field)
    enddo
  endif

  !------------------------------------
  ! Grid field calc 

  if (associated(ele%grid_field)) then
  
    expt_ptr => expt  ! To get around ifort bug where debug info for variables used in contained routines is missing.
    ! loop over grid modes

    do i = 1, size(ele%grid_field)
      g_field => ele%grid_field(i)
      g_field_ptr => ele%grid_field(i)  ! To get around ifort bug.

      if (g_field%harmonic /= 0) then
        freq0 = ele%value(rf_frequency$)
        freq = freq0 * g_field%harmonic
        if (freq0 == 0) then
          call out_io (s_fatal$, r_name, 'ELEMENT FREQUENCY IS ZERO BUT GRID_FIELD HARMONIC IS NOT FOR: ' // ele%name)
          if (global_com%exit_on_error) call err_exit
          if (present(err_flag)) err_flag = .true.
          return  
        endif

        t_ref = (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_err$) + &
                                                  phi0_autoscale + g_field%phi0_fieldmap) / freq0
        if (ele%key == rfcavity$) t_ref = 0.25/freq0 - t_ref
      endif

      call to_fieldmap_coords (ele, local_orb, s_body, g_field%ele_anchor_pt, g_field%r0, g_field%curved_ref_frame, x, y, z, cos_ang, sin_ang, err)
      if (err) then
        if (present(err_flag)) err_flag = .true.
        return
      endif

      ! DC modes should have g_field%harmonic = 0

      expt = field_autoscale * g_field%field_scale * master_parameter_value(g_field%master_parameter, ele)
      if (g_field%harmonic /= 0) expt = expt * exp(-I_imaginary * twopi * (freq * (time + t_ref)))

      ! calculate field based on grid type
      select case(g_field%geometry)

      case (xyz$)
      
        call grid_field_interpolate(ele, local_orb, g_field, g_pt, err, x, y, z, &
                    allow_s_out_of_bounds = grid_allow_s_out_of_bounds, print_err = print_err)
        if (err) then
          if (present(err_flag)) err_flag = .true.
          orbit%state = local_orb%state
          return
        endif

        mode_field%e = real(expt * g_pt%e)
        mode_field%b = real(expt * g_pt%B)

      case(rotationally_symmetric_rz$)
        
        ! Format should be: pt (ir, iz) = ( Er, 0, Ez, 0, Bphi, 0 ) 
          
        ! Interpolate 2D (r, z) grid
        ! g_pt is a grid_field_pt_struct, which has complex E and B

        r = sqrt(x**2 + y**2)

        call grid_field_interpolate(ele, local_orb, g_field, g_pt, err, r, z, &
                     allow_s_out_of_bounds = grid_allow_s_out_of_bounds, print_err = print_err)
        if (err) then
          if (present(err_flag)) err_flag = .true.
          orbit%state = local_orb%state
          return
        endif

        ! Transverse field is zero on axis. Otherwise:

        if (r /= 0) then
          ! Get non-rotated field
          E_rho = real(expt * g_pt%E(1))
          E_phi = real(expt * g_pt%E(2))
          B_rho = real(expt * g_pt%B(1)) 
          B_phi = real(expt * g_pt%B(2))

          ! rotate field and output Ex, Ey, Bx, By
          mode_field%e(1) = (x*E_rho - y*E_phi)/r
          mode_field%e(2) = (y*E_rho + x*E_phi)/r
          mode_field%b(1) = (x*B_rho - y*B_phi)/r
          mode_field%b(2) = (y*B_rho + x*B_phi)/r
        endif
    
        ! Ez, Bz 
        mode_field%e(3) = real(expt*g_pt%E(3))
        mode_field%b(3) = real(expt*g_pt%B(3)) 
    
        ! Vector potential.
        ! Right now only good for DC magnetic fields. This should be generalized.

        if (logic_option(.false., calc_potential)) then
          if (r /= 0) then
            abs_tol = abs(1e-10_rp * r * orbit%p0c * (1 + orbit%vec(6)) / (c_light * charge_of(ele%ref_species)))
            inte = super_qromb(rb_field, 0.0_rp, r, 1e-12_rp, abs_tol, 2, err) / r
            field%A(1:2) = field%A(1:2) + inte * [-y, x] / r
          endif
        endif

      case default
        call out_io (s_fatal$, r_name, 'UNKNOWN GRID GEOMETRY: \i0\ ', &
                                       'FOR ELEMENT: ' // ele%name, i_array = [g_field%geometry])
        if (global_com%exit_on_error) call err_exit
        if (present(err_flag)) err_flag = .true.
        return
      end select
      
      if ((ele%key == sbend$ .or. ele%key == rf_bend$) .and. .not. g_field%curved_ref_frame) call restore_curvilinear_field(mode_field%E, mode_field%B)

      field%E = field%E + mode_field%E
      field%B = field%B + mode_field%B

    enddo
  endif

! Beginning_ele, for example, has no field

case (no_field$)

  return

! Unknown field_calc

case default
  call out_io (s_fatal$, r_name, 'BAD FIELD_CALC METHOD FOR ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
  if (present(err_flag)) err_flag = .true.
  return
end select

! Scale ac_kicker element field

if (ele%key == ac_kicker$) then
  a_amp = ac_kicker_amp (ele, orbit, rf_time)

  field%E = a_amp * field%E
  field%B = a_amp * field%B
  field%dE = a_amp * field%dE
  field%dB = a_amp * field%dB
endif

!----------------------------------------------------------------------------------------------
! overlapping of fields from other elements

8000 continue

if (ele%n_lord_field /= 0 .and. logic_option(.true., use_overlap)) then
  lab_orb = orbit
  if (local_ref_frame) then
    call offset_particle (ele, unset$, lab_orb, set_hvkicks = .false., s_pos = s_body, s_out = s_lab)
  else
    s_lab = s_body
  endif

  lab_position%r = [lab_orb%vec(1), lab_orb%vec(3), s_lab]
  global_position = coords_local_curvilinear_to_floor (lab_position, ele, w_mat = w_ele_mat, calculate_angles = .false.)

  lord_orb = lab_orb
  do i = 1, ele%n_lord_field
    lord => pointer_to_lord(ele, i, lord_type = field_lord$)
    lord_position = coords_floor_to_local_curvilinear (global_position, lord, status, w_lord_mat)
    lord_orb%vec(1) = lord_position%r(1)
    lord_orb%vec(3) = lord_position%r(2)
    ! Set use_overlap = False to prevent recursion.
    call em_field_calc (lord, param, lord_position%r(3), lord_orb, .false., l1_field, calc_dfield, err, &
          use_overlap = .false., grid_allow_s_out_of_bounds = .true., used_eles = used_eles, &
          print_err = print_err)
    if (err) then
      if (present(err_flag)) err_flag = .true.
      return
    endif
    ! Field in lord lab coords to field in global coords
    call rotate_em_field (l1_field, transpose(w_lord_mat), w_lord_mat, calc_dfield)
    if (i == 1) then
      lord_field = l1_field
    else
      lord_field = lord_field + l1_field
    endif
  enddo

  ! Field in global coords to field in lab coords

  call rotate_em_field (lord_field, transpose(w_ele_mat), transpose(w_ele_mat))

  if (local_ref_frame) then
    call convert_field_ele_to_lab (ele, s_lab, .false., lord_field, calc_dfield, calc_potential)  ! lab -> ele
    field = field + lord_field
  else
    call convert_field_ele_to_lab (ele, s_body, .true., field, calc_dfield, calc_potential)
    field = field + lord_field
  endif

  return
endif

! Final

if (do_df_calc .and. .not. dfield_computed) then
  call em_field_derivatives (ele, param, s_pos, orbit, .true., field, grid_allow_s_out_of_bounds, rf_time)
endif

if (.not. local_ref_frame) call convert_field_ele_to_lab (ele, s_body, .true., field, calc_dfield, calc_potential)

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
contains

! Function for vector potential calc.

function rb_field(x)

real(rp), intent(in) :: x(:)
real(rp) :: rb_field(size(x))
integer i

!

do i = 1, size(x)
  call grid_field_interpolate(ele, local_orb, g_field_ptr, g_pt, err, x(i), z, &
              allow_s_out_of_bounds = .true., print_err = print_err)
  rb_field(i) = x(i) * expt_ptr * g_pt%b(3)
enddo

end function rb_field

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

! restore_curvilinear_field(field_a, field_b)
!
! For sbend with Grid calculation.

subroutine restore_curvilinear_field(field_a, field_b)

real(rp) temp, field_a(3)
real(rp), optional :: field_b(3)

! For sbend with Grid calculation Restores x and s_body, and rotates output fields.

if (ele%value(g$) == 0) return

temp       = field_a(3)*cos_ang - field_a(1)*sin_ang
field_a(1) = field_a(3)*sin_ang + field_a(1)*cos_ang
field_a(3) = temp

if (present(field_b)) then
  temp       = field_b(3)*cos_ang - field_b(1)*sin_ang
  field_b(1) = field_b(3)*sin_ang + field_b(1)*cos_ang
  field_b(3) = temp 
endif

end subroutine restore_curvilinear_field

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine gen_grad_add_em_field (ele, gg_map, r_pos, s_body, field)

type (ele_struct) ele
type (gen_grad_map_struct), target :: gg_map
type (gen_grad1_struct), pointer :: gg
type (em_field_struct) field

real(rp) r_pos(3), z_rel, theta, rho, s_body
real(rp), allocatable :: d0(:), d1(:), der(:)
integer iz0, j, id, nd

!

iz0 = floor(r_pos(3) / gg_map%dz)
if (iz0 < gg_map%iz0) iz0 = iz0 + 1 ! Allow one dz width out-of-bounds.
if (iz0 >= gg_map%iz1) iz0 = iz0 - 1 ! Allow one dz width out-of-bounds.

if (iz0 < gg_map%iz0 .or. iz0 >= gg_map%iz1) then
  if (.not. logic_option(.false., grid_allow_s_out_of_bounds)) then
    call out_io (s_error$, r_name, 'PARTICLE Z  \F10.3\ POSITION OUT OF BOUNDS.', &
                                   'FOR GEN_GRAD_MAP IN ELEMENT: ' // ele%name, r_array = [s_body])
  endif
  return
endif

!

z_rel = r_pos(3) - iz0 * gg_map%dz

theta = atan2(r_pos(2), r_pos(1))
rho = norm2(r_pos(1:2))

do j = 1, size(gg_map%gg)
  gg => gg_map%gg(j)
  nd = gg%n_deriv_max

  call re_allocate2(der, 0, nd, .false.)
  do id = 0, nd
    der(id) = poly_eval(gg%deriv(iz0, id:), z_rel, diff_coef=.true.)
  enddo

  if (gg_map%field_type == magnetic$) then
    field%B = field%B + gen_grad_field (der, gg, rho, theta) * (gg_map%field_scale * master_parameter_value(gg_map%master_parameter, ele))
  else
    field%E = field%E + gen_grad_field (der, gg, rho, theta) * (gg_map%field_scale * master_parameter_value(gg_map%master_parameter, ele))
  endif
enddo

end subroutine gen_grad_add_em_field

end subroutine em_field_calc 
