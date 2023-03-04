!+
! Function tao_param_value_at_s (data_type, ele, orbit, err_flag, why_invalid, print_err) result (value)
!
! Routine to evaluate a parameter at a lattice s-position.
!
! Input:
!   data_type     -- character(*): Parameter name.
!   ele           -- ele_struct: Lattice element whose exit end is at the evaluation s-position.
!   orbit         -- coord_struct: Orbit at the evaluation s-position.
!
! Output:
!   err_flag      -- logical: Set true if data_type does not have a corresponding Bmad parameter.
!   value         -- real(rp): Parameter value.
!   why_invalid   -- character(*), optional: Set if  err_flag = True to document why is there a problem.
!   print_err     -- logical, optional: Print error message on error? Default is True.
!-

function tao_param_value_at_s (data_type, ele, orbit, err_flag, why_invalid, print_err) result (value)

use tao_interface, except_dummy => tao_param_value_at_s
use measurement_mod
use em_field_mod, only: em_field_derivatives

implicit none

type (ele_struct) ele
type (coord_struct) orbit, orb
type (bpm_phase_coupling_struct) bpm_data
type (floor_position_struct) floor
type (branch_struct), pointer :: branch
type (em_field_struct) field, field0, field1

real(rp) value, cbar(2,2), f, amp_a, amp_b, amp_na, amp_nb, time, dt, amp, phase

character(*) data_type
character(*), optional :: why_invalid
character(40) name, prefix, d_type

integer i, j, ix

logical err_flag
logical, optional :: print_err

!

err_flag = .false.

if (data_type == 'state') then
  value = orbit%state
  return
endif

if (orbit%state /= alive$) then
  err_flag = .true.
  if (present(why_invalid)) why_invalid = 'PARTICLE DEAD AT ELEMENT: ' // ele%name
  return
endif

branch => pointer_to_branch(ele)

d_type = data_type
if (d_type(1:6) == 'orbit_') then;           d_type(1:6) = 'orbit.'
elseif (d_type(1:5) == 'spin_') then;        d_type(1:5) = 'spin.'
elseif (d_type(1:10) == 'intensity_') then;  d_type(1:10) = 'intensity.'
endif

ix = index(d_type, '.')
if (ix == 0) then
  prefix = d_type
else
  prefix = d_type(1:ix-1)
endif

!

select case (prefix)

case ('alpha')
  select case (d_type)
  case ('alpha.a');          value = ele%a%alpha
  case ('alpha.b');          value = ele%b%alpha
  case ('alpha.z');          value = ele%z%alpha
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('b_curl', 'b0_curl')
  orb = orbit
  if (prefix == 'b0_curl') then
    time = orbit%t
  else
    time = particle_rf_time(orb, ele, (ele%field_calc /= fieldmap$), orb%s-ele%s_start)
  endif

  call em_field_derivatives (ele, branch%param, orb%s-ele%s_start, orb, .false., field, rf_time = time)
  dt = bmad_com%d_orb(5) / c_light
  call em_field_calc (ele, branch%param, orb%s-ele%s_start, orb, .false., field0, rf_time = time-dt)
  call em_field_calc (ele, branch%param, orb%s-ele%s_start, orb, .false., field1, rf_time = time+dt)

  select case (d_type)
  case ('b_curl.x', 'b0_curl.x');  value = field%dB(3,2) - field%dB(2,3) - (field1%E(1) - field0%E(1)) / (2 * dt * c_light**2)
  case ('b_curl.y', 'b0_curl.y');  value = field%dB(1,3) - field%dB(3,1) - (field1%E(2) - field0%E(2)) / (2 * dt * c_light**2)
  case ('b_curl.z', 'b0_curl.z');  value = field%dB(2,1) - field%dB(1,2) - (field1%E(3) - field0%E(3)) / (2 * dt * c_light**2)
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('b_div', 'b0_div')
  orb = orbit
  if (prefix == 'b0_div') then
    time = orbit%t
  else
    time = particle_rf_time(orb, ele, (ele%field_calc /= fieldmap$), orb%s-ele%s_start)
  endif

  call em_field_derivatives (ele, branch%param, orb%s-ele%s_start, orb, .false., field, rf_time = time)
  value = field%dB(1,1) + field%dB(2,2) + field%dB(3,3)

case ('b_field', 'b0_field')
  orb = orbit
  if (prefix == 'b0_field') then
    time = orbit%t
  else
    time = particle_rf_time(orb, ele, (ele%field_calc /= fieldmap$), orb%s-ele%s_start)
  endif

  call em_field_calc (ele, branch%param, orb%s-ele%s_start, orb, .false., field, &
                                                  err_flag = err_flag, rf_time = time, print_err = print_err)
  if (err_flag) return
  select case (d_type)
  case ('b_field.x', 'b0_field.x');  value = field%b(1)
  case ('b_field.y', 'b0_field.y');  value = field%b(2)
  case ('b_field.z', 'b0_field.z');  value = field%b(3)
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('beta')
  select case (d_type)
  case ('beta');             value = orbit%beta
  case ('beta.a');           value = ele%a%beta
  case ('beta.b');           value = ele%b%beta
  case ('beta.z');           value = ele%z%beta
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('bpm_cbar')
  call to_phase_and_coupling_reading (ele, s%com%add_measurement_noise, bpm_data, err_flag)
  select case (d_type)
  case ('bpm_cbar.22a');     value = bpm_data%cbar22_a
  case ('bpm_cbar.12a');     value = bpm_data%cbar12_a
  case ('bpm_cbar.11b');     value = bpm_data%cbar11_b
  case ('bpm_cbar.12b');     value = bpm_data%cbar12_b
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('bpm_eta')
  select case (d_type)
  case ('bpm_eta.x');  call to_eta_reading ([ele%x%eta, ele%y%eta], ele, x_plane$, s%com%add_measurement_noise, value, err_flag)
  case ('bpm_eta.y');  call to_eta_reading ([ele%x%eta, ele%y%eta], ele, y_plane$, s%com%add_measurement_noise, value, err_flag)
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('bpm_orbit')
  select case (d_type)
  case ('bpm_orbit.x');      call to_orbit_reading (orbit, ele, x_plane$, s%com%add_measurement_noise, value, err_flag)
  case ('bpm_orbit.y');      call to_orbit_reading (orbit, ele, y_plane$, s%com%add_measurement_noise, value, err_flag)
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('bpm_phase')
  call to_phase_and_coupling_reading (ele, s%com%add_measurement_noise, bpm_data, err_flag)
  if (err_flag) return
  select case (d_type)
  case ('bpm_phase.a');      value = bpm_data%phi_a
  case ('bpm_phase.b');      value = bpm_data%phi_b
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('bpm_k')
  call to_phase_and_coupling_reading (ele, s%com%add_measurement_noise, bpm_data, err_flag)
  select case (d_type)
  case ('bpm_k.22a');        value = bpm_data%k_22a
  case ('bpm_k.12a');        value = bpm_data%k_12a
  case ('bpm_k.11b');        value = bpm_data%k_11b
  case ('bpm_k.12b');        value = bpm_data%k_12b
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('c_mat')
  select case (d_type)
  case ('cmat.11');         value = ele%c_mat(1,1)
  case ('cmat.12');         value = ele%c_mat(1,2)
  case ('cmat.21');         value = ele%c_mat(2,1)
  case ('cmat.22');         value = ele%c_mat(2,2)
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('cbar')
  call c_to_cbar (ele, cbar)
  select case (d_type)
  case ('cbar.11');          value = cbar(1,1)
  case ('cbar.12');          value = cbar(1,2)
  case ('cbar.21');          value = cbar(2,1)
  case ('cbar.22');          value = cbar(2,2)
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('coupling')
  call c_to_cbar (ele, cbar)  
  select case (d_type)
  case ('coupling.11b');  value = cbar(1,1) * sqrt(ele%a%beta/ele%b%beta) / ele%gamma_c
  case ('coupling.12a');  value = cbar(1,2) * sqrt(ele%b%beta/ele%a%beta) / ele%gamma_c
  case ('coupling.12b');  value = cbar(1,2) * sqrt(ele%a%beta/ele%b%beta) / ele%gamma_c
  case ('coupling.22a');  value = cbar(2,2) * sqrt(ele%b%beta/ele%a%beta) / ele%gamma_c
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('curly_h')
  select case (d_type)
  case ('curly_h.a');     value = ele%a%gamma * ele%a%eta**2 + 2 * ele%a%alpha * ele%a%eta * ele%a%etap + ele%a%beta * ele%a%etap**2
  case ('curly_h.b');     value = ele%b%gamma * ele%b%eta**2 + 2 * ele%b%alpha * ele%b%eta * ele%b%etap + ele%b%beta * ele%b%etap**2
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('e_curl', 'e0_curl')
  orb = orbit
  if (prefix == 'e0_curl') then
    time = orbit%t
  else
    time = particle_rf_time(orb, ele, (ele%field_calc /= fieldmap$), orb%s-ele%s_start)
  endif

  call em_field_derivatives (ele, branch%param, orb%s-ele%s_start, orb, .false., field, rf_time = time)
  time = particle_rf_time(orb, ele, (ele%field_calc /= fieldmap$), orb%s-ele%s_start)
  dt = bmad_com%d_orb(5) / c_light
  call em_field_calc (ele, branch%param, orb%s-ele%s_start, orb, .false., field0, rf_time = time-dt)
  call em_field_calc (ele, branch%param, orb%s-ele%s_start, orb, .false., field1, rf_time = time+dt)

  select case (d_type)
  case ('e_curl.x', 'e0_curl.x');  value = field%dE(3,2) - field%dE(2,3) + (field1%B(1) - field0%B(1)) / (2 * dt)
  case ('e_curl.y', 'e0_curl.y');  value = field%dE(1,3) - field%dE(3,1) + (field1%B(2) - field0%B(2)) / (2 * dt)
  case ('e_curl.z', 'e0_curl.z');  value = field%dE(2,1) - field%dE(1,2) + (field1%B(3) - field0%B(3)) / (2 * dt)
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('e_div', 'e0_div')
  orb = orbit
  if (prefix == 'e0_div') then
    time = orbit%t
  else
    time = particle_rf_time(orb, ele, (ele%field_calc /= fieldmap$), orb%s-ele%s_start)
  endif

  call em_field_derivatives (ele, branch%param, orb%s-ele%s_start, orb, .false., field, rf_time = time)
  value = field%dE(1,1) + field%dE(2,2) + field%dE(3,3)

case ('e_field', 'e0_field')
  orb = orbit
  if (prefix == 'e0_field') then
    time = orbit%t
  else
    time = particle_rf_time(orb, ele, (ele%field_calc /= fieldmap$), orb%s-ele%s_start)
  endif

  call em_field_calc (ele, branch%param, orb%s-ele%s_start, orb, .false., field, rf_time = time)
  select case (d_type)
  case ('e_field.x', 'e0_field.x');  value = field%e(1)
  case ('e_field.y', 'e0_field.y');  value = field%e(2)
  case ('e_field.z', 'e0_field.z');  value = field%e(3)
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('e_tot_ref');          value = ele%value(e_tot$)

case ('energy')
  if (orbit%species == photon$) then
    value = orbit%p0c
  else
    call convert_pc_to(orbit%p0c * (1 + orbit%vec(6)), orbit%species, e_tot = value)
  endif

case ('eta')
  select case (d_type)
  case ('eta.a');            value = ele%a%eta
  case ('eta.b');            value = ele%b%eta
  case ('eta.x');            value = ele%x%eta
  case ('eta.y');            value = ele%y%eta
  case ('eta.z');            value = ele%z%eta
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('etap')
  select case (d_type)
  case ('etap.a');           value = ele%a%etap
  case ('etap.b');           value = ele%b%etap
  case ('etap.x');           value = ele%x%etap
  case ('etap.y');           value = ele%y%etap
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('floor')
  select case (d_type)
  case ('floor.x');          value = ele%floor%r(1)
  case ('floor.y');          value = ele%floor%r(2)
  case ('floor.z');          value = ele%floor%r(3)
  case ('floor.theta');      value = ele%floor%theta
  case ('floor.phi');        value = ele%floor%phi
  case ('floor.psi');        value = ele%floor%psi
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('floor_actual')
  floor = ele_geometry_with_misalignments(ele)
  select case (d_type)
  case ('floor_actual.x');          value = floor%r(1)
  case ('floor_actual.y');          value = floor%r(2)
  case ('floor_actual.z');          value = floor%r(3)
  case ('floor_actual.theta');      value = floor%theta
  case ('floor_actual.phi');        value = floor%phi
  case ('floor_actual.psi');        value = floor%psi
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('floor_orbit')
  if (orbit%location == downstream_end$) then
    floor = orbit_to_local_curvilinear(orbit, ele, relative_to = orbit%location)
    floor = coords_local_curvilinear_to_floor (floor, ele, .false., relative_to = orbit%location)
  else
    floor = orbit_to_local_curvilinear(orbit, ele, orbit%direction)
    floor = coords_local_curvilinear_to_floor (floor, ele, .false.)
  endif
  select case (d_type)
  case ('floor_orbit.x');          value = floor%r(1)
  case ('floor_orbit.y');          value = floor%r(2)
  case ('floor_orbit.z');          value = floor%r(3)
  case ('floor_orbit.theta');      value = floor%theta
  case ('floor_orbit.phi');        value = floor%phi
  case ('floor_orbit.psi');        value = floor%psi
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('gamma')
  select case (d_type)
  case ('gamma.a');          value = ele%a%gamma
  case ('gamma.b');          value = ele%b%gamma
  case ('gamma.z');          value = ele%z%gamma
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('intensity')
  select case (d_type)
  case ('intensity');                     value = orbit%field(1)**2 + orbit%field(2)**2
  case ('intensity_x', 'intensity.x');    value = orbit%field(1)**2
  case ('intensity_y', 'intensity.y');    value = orbit%field(2)**2
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('k')
  call c_to_cbar (ele, cbar)
  f = sqrt(ele%a%beta/ele%b%beta) 
  select case (d_type)
  case ('k.11b');            value = cbar(1,1) * f / ele%gamma_c
  case ('k.12a');            value = cbar(1,2) / (f * ele%gamma_c)
  case ('k.12b');            value = cbar(1,2) * f / ele%gamma_c
  case ('k.22a');            value = cbar(2,2) / (f * ele%gamma_c)
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('momentum');            value = (1 + orbit%vec(6)) * orbit%p0c

case ('orbit')
  if (d_type(7:9) == 'amp' .or. d_type(7:9) == 'nor') &
          call orbit_amplitude_calc (ele, orbit, amp_a, amp_b, amp_na, amp_nb)
  select case (d_type)
  case ('orbit.x');           value = orbit%vec(1)
  case ('orbit.y');           value = orbit%vec(3)
  case ('orbit.z');           value = orbit%vec(5)
  case ('orbit.px');          value = orbit%vec(2)
  case ('orbit.py');          value = orbit%vec(4)
  case ('orbit.pz');          value = orbit%vec(6)
  case ('orbit.amp_a');       value = amp_a
  case ('orbit.amp_b');       value = amp_b
  case ('orbit.norm_amp_a');  value = amp_na
  case ('orbit.norm_amp_b');  value = amp_nb
  case ('orbit.e_tot', 'orbit.energy', 'orbit.kinetic')  ! orbit.e_tot is old style
    if (orbit%beta == 0) then
      value = mass_of(branch%param%particle)
    else
      value = orbit%p0c * (1 + orbit%vec(6)) / orbit%beta
    endif
    if (d_type == 'orbit.kinetic') value = value - mass_of(orbit%species)

  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('pc')
  if (orbit%species == photon$) then
    value = orbit%p0c
  else
    value = orbit%p0c * (1 + orbit%vec(6))
  endif

case ('phase', 'phase_frac')
  select case (d_type)
  case ('phase.a');           value = ele%a%phi
  case ('phase_frac.a');      value = modulo2 (ele%a%phi, pi)
  case ('phase.b');           value = ele%b%phi
  case ('phase_frac.b');      value = modulo2 (ele%b%phi, pi)
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('ping_a')
  call c_to_cbar (ele, cbar)
  select case (d_type)
  case ('ping_a.amp_x');          value = ele%gamma_c * sqrt(ele%a%beta)
  case ('ping_a.phase_x');        value = ele%a%phi
  case ('ping_a.amp_y');          value = sqrt(ele%b%beta * (cbar(1,2)**2 + cbar(2,2)**2))
  case ('ping_a.phase_y');        value = ele%a%phi + atan2(cbar(1,2), -cbar(2,2))
  case ('ping_a.amp_sin_rel_y');  value = -sqrt(ele%b%beta) * cbar(1,2)
  case ('ping_a.amp_cos_rel_y');  value = -sqrt(ele%b%beta) * cbar(2,2)
  case ('ping_a.amp_sin_y')
    amp = sqrt(ele%b%beta * (cbar(1,2)**2 + cbar(2,2)**2))
    phase = ele%a%phi + atan2(cbar(1,2), -cbar(2,2))
    value = amp * sin(phase)
  case ('ping_a.amp_cos_y')
    amp = sqrt(ele%b%beta * (cbar(1,2)**2 + cbar(2,2)**2))
    phase = ele%a%phi + atan2(cbar(1,2), -cbar(2,2))
    value = amp * cos(phase)
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('ping_b')
  call c_to_cbar (ele, cbar)
  select case (d_type)
  case ('ping_b.amp_y');          value = ele%gamma_c * sqrt(ele%b%beta)
  case ('ping_b.phase_y');        value = ele%b%phi
  case ('ping_b.amp_x');          value = sqrt(ele%a%beta * (cbar(1,2)**2 + cbar(1,1)**2))
  case ('ping_b.phase_x');        value = ele%b%phi + atan2(cbar(1,2), cbar(1,1))
  case ('ping_b.amp_sin_rel_x');  value = -sqrt(ele%a%beta) * cbar(1,2)
  case ('ping_b.amp_cos_rel_x');  value = sqrt(ele%a%beta) * cbar(1,1)
  case ('ping_b.amp_sin_x')
    amp = sqrt(ele%a%beta * (cbar(1,2)**2 + cbar(1,1)**2))
    phase = ele%b%phi + atan2(cbar(1,2), cbar(1,1))
    value = amp * sin(phase)
  case ('ping_b.amp_cos_x')
    amp = sqrt(ele%a%beta * (cbar(1,2)**2 + cbar(1,1)**2))
    phase = ele%b%phi + atan2(cbar(1,2), cbar(1,1))
    value = amp * cos(phase)
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('r')
  i = tao_read_phase_space_index (data_type, 3, .false.)
  j = tao_read_phase_space_index (data_type, 4, .false.)
  if (i == 0 .or. j == 0 .or. len_trim(data_type) /= 4) then
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'BAD DATA_TYPE = "' // trim(data_type)
    return
  endif

  value = ele%mat6(i,j)

case ('ref_p0c')
  value = orbit%p0c

case ('ref_time')
  value = ele%ref_time

case ('spin')
  select case (d_type)
  case ('spin.x');        value = orbit%spin(1)
  case ('spin.y');        value = orbit%spin(2)
  case ('spin.z');        value = orbit%spin(3)
  case ('spin.amp');      value = norm2(orbit%spin)
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

  if (.not. bmad_com%spin_tracking_on) then
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'SPIN TRACKING IS NOT ENABLED.'
  endif

case ('s_position');         value = ele%s

case ('time', 't');          value = orbit%t

case ('velocity')
  select case (d_type)
  case ('velocity');  value = orbit%beta
  case ('velocity.x');  value = orbit%vec(2) * (1 + orbit%vec(6)) * orbit%beta
  case ('velocity.y');  value = orbit%vec(4) * (1 + orbit%vec(6)) * orbit%beta
  case ('velocity.z');  value = sqrt(1 - (orbit%vec(2) * (1 + orbit%vec(6)))**2 - (orbit%vec(4) * (1 + orbit%vec(6)))**2) * orbit%beta
  case default
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case default
  err_flag = .true.
  if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)

end select

end function
