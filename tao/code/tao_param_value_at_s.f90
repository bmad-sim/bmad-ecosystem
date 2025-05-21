!+
! Function tao_param_value_at_s (data_type, ele_to_s, ele_here, orbit, err_flag, why_invalid, print_err, bad_datum) result (value)
!
! Routine to evaluate a parameter at a lattice s-position given by ele_to_s%s which is the same at orbit%s.
!
! Typically, ele_to_s is a slice_slave. However, sometimes ele_to_s could be a hybrid$ element representing the integrated Twiss 
! up to ele_to_s%s. In this case, if the electric or magnetic field is to be calculated, needed is the ele_here lattice element 
! that overlaps the position ele_to_s%s. If ele_to_s is a slice_slave, ele_here can be equal to ele_to_s.
!
! Input:
!   data_type     -- character(*): Parameter name.
!   ele_to_s      -- ele_struct: Element whose exit end is at the evaluation s-position.
!   ele_here      -- ele_struct: Lattice element that overlaps the s-position ele%s.
!   orbit         -- coord_struct: Orbit at the evaluation s-position.
!
! Output:
!   err_flag      -- logical: Set true if parameter cannot be evaluated.
!   value         -- real(rp): Parameter value.
!   why_invalid   -- character(*), optional: Set if  err_flag = True to document why is there a problem.
!   print_err     -- logical, optional: Print error message on error? Default is True.
!   bad_datum     -- logical, optional: Data_type is malformed.
!-

function tao_param_value_at_s (data_type, ele_to_s, ele_here, orbit, err_flag, why_invalid, print_err, bad_datum) result (value)

use tao_interface, except_dummy => tao_param_value_at_s
use measurement_mod
use em_field_mod, only: em_field_derivatives

implicit none

type (ele_struct) ele_to_s, ele_here
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
logical, optional :: print_err, bad_datum

!

err_flag = .false.
if (present(bad_datum)) bad_datum = .false.

if (data_type == 'state') then
  value = orbit%state
  return
endif

if (orbit%state /= alive$) then
  err_flag = .true.
  if (present(why_invalid)) why_invalid = 'PARTICLE DEAD AT ELEMENT: ' // ele_to_s%name
  return
endif

branch => pointer_to_branch(ele_here)

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
  case ('alpha.a');          value = ele_to_s%a%alpha
  case ('alpha.b');          value = ele_to_s%b%alpha
  case ('alpha.z');          value = ele_to_s%z%alpha
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('b_curl', 'b0_curl')
  orb = orbit
  if (prefix == 'b0_curl') then
    time = orbit%t
  else
    time = particle_rf_time(orb, ele_here, (ele_here%field_calc /= fieldmap$), orb%s-ele_here%s_start)
  endif

  call em_field_derivatives (ele_here, branch%param, orb%s-ele_here%s_start, orb, .false., field, rf_time = time)
  dt = bmad_com%d_orb(5) / c_light
  call em_field_calc (ele_here, branch%param, orb%s-ele_here%s_start, orb, .false., field0, rf_time = time-dt)
  call em_field_calc (ele_here, branch%param, orb%s-ele_here%s_start, orb, .false., field1, rf_time = time+dt)

  select case (d_type)
  case ('b_curl.x', 'b0_curl.x');  value = field%dB(3,2) - field%dB(2,3) - (field1%E(1) - field0%E(1)) / (2 * dt * c_light**2)
  case ('b_curl.y', 'b0_curl.y');  value = field%dB(1,3) - field%dB(3,1) - (field1%E(2) - field0%E(2)) / (2 * dt * c_light**2)
  case ('b_curl.z', 'b0_curl.z');  value = field%dB(2,1) - field%dB(1,2) - (field1%E(3) - field0%E(3)) / (2 * dt * c_light**2)
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('b_div', 'b0_div')
  orb = orbit
  if (prefix == 'b0_div') then
    time = orbit%t
  else
    time = particle_rf_time(orb, ele_here, (ele_here%field_calc /= fieldmap$), orb%s-ele_here%s_start)
  endif

  call em_field_derivatives (ele_here, branch%param, orb%s-ele_here%s_start, orb, .false., field, rf_time = time)
  value = field%dB(1,1) + field%dB(2,2) + field%dB(3,3)

case ('b_field', 'b0_field')
  orb = orbit
  if (prefix == 'b0_field') then
    time = orbit%t
  else
    time = particle_rf_time(orb, ele_here, (ele_here%field_calc /= fieldmap$), orb%s-ele_here%s_start)
  endif

  call em_field_calc (ele_here, branch%param, orb%s-ele_here%s_start, orb, .false., field, &
                                                  err_flag = err_flag, rf_time = time, print_err = print_err)
  if (err_flag) return
  select case (d_type)
  case ('b_field.x', 'b0_field.x');  value = field%b(1)
  case ('b_field.y', 'b0_field.y');  value = field%b(2)
  case ('b_field.z', 'b0_field.z');  value = field%b(3)
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('beta')
  select case (d_type)
  case ('beta');             value = orbit%beta
  case ('beta.a');           value = ele_to_s%a%beta
  case ('beta.b');           value = ele_to_s%b%beta
  case ('beta.z');           value = ele_to_s%z%beta
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('bpm_cbar')
  call to_phase_and_coupling_reading (ele_to_s, s%com%add_measurement_noise, bpm_data, err_flag)
  select case (d_type)
  case ('bpm_cbar.22a');     value = bpm_data%cbar22_a
  case ('bpm_cbar.12a');     value = bpm_data%cbar12_a
  case ('bpm_cbar.11b');     value = bpm_data%cbar11_b
  case ('bpm_cbar.12b');     value = bpm_data%cbar12_b
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('bpm_eta')
  select case (d_type)
  case ('bpm_eta.x');  call to_eta_reading ([ele_to_s%x%eta, ele_to_s%y%eta], ele_to_s, x_plane$, s%com%add_measurement_noise, value, err_flag)
  case ('bpm_eta.y');  call to_eta_reading ([ele_to_s%x%eta, ele_to_s%y%eta], ele_to_s, y_plane$, s%com%add_measurement_noise, value, err_flag)
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('bpm_orbit')
  select case (d_type)
  case ('bpm_orbit.x');      call to_orbit_reading (orbit, ele_to_s, x_plane$, s%com%add_measurement_noise, value, err_flag)
  case ('bpm_orbit.y');      call to_orbit_reading (orbit, ele_to_s, y_plane$, s%com%add_measurement_noise, value, err_flag)
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('bpm_phase')
  call to_phase_and_coupling_reading (ele_to_s, s%com%add_measurement_noise, bpm_data, err_flag)
  if (err_flag) return
  select case (d_type)
  case ('bpm_phase.a');      value = bpm_data%phi_a
  case ('bpm_phase.b');      value = bpm_data%phi_b
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('bpm_k')
  call to_phase_and_coupling_reading (ele_to_s, s%com%add_measurement_noise, bpm_data, err_flag)
  select case (d_type)
  case ('bpm_k.22a');        value = bpm_data%k_22a
  case ('bpm_k.12a');        value = bpm_data%k_12a
  case ('bpm_k.11b');        value = bpm_data%k_11b
  case ('bpm_k.12b');        value = bpm_data%k_12b
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('c_mat')
  select case (d_type)
  case ('cmat.11');         value = ele_to_s%c_mat(1,1)
  case ('cmat.12');         value = ele_to_s%c_mat(1,2)
  case ('cmat.21');         value = ele_to_s%c_mat(2,1)
  case ('cmat.22');         value = ele_to_s%c_mat(2,2)
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('cbar')
  call c_to_cbar (ele_to_s, cbar)
  select case (d_type)
  case ('cbar.11');          value = cbar(1,1)
  case ('cbar.12');          value = cbar(1,2)
  case ('cbar.21');          value = cbar(2,1)
  case ('cbar.22');          value = cbar(2,2)
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('coupling')
  call c_to_cbar (ele_to_s, cbar)  
  select case (d_type)
  case ('coupling.11b');  value = cbar(1,1) * sqrt(ele_to_s%a%beta/ele_to_s%b%beta) / ele_to_s%gamma_c
  case ('coupling.12a');  value = cbar(1,2) * sqrt(ele_to_s%b%beta/ele_to_s%a%beta) / ele_to_s%gamma_c
  case ('coupling.12b');  value = cbar(1,2) * sqrt(ele_to_s%a%beta/ele_to_s%b%beta) / ele_to_s%gamma_c
  case ('coupling.22a');  value = cbar(2,2) * sqrt(ele_to_s%b%beta/ele_to_s%a%beta) / ele_to_s%gamma_c
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('curly_h')
  select case (d_type)
  case ('curly_h.a');     value = ele_to_s%a%gamma * ele_to_s%a%eta**2 + 2 * ele_to_s%a%alpha * ele_to_s%a%eta * ele_to_s%a%etap + ele_to_s%a%beta * ele_to_s%a%etap**2
  case ('curly_h.b');     value = ele_to_s%b%gamma * ele_to_s%b%eta**2 + 2 * ele_to_s%b%alpha * ele_to_s%b%eta * ele_to_s%b%etap + ele_to_s%b%beta * ele_to_s%b%etap**2
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('deta_ds')
  select case (d_type)
  case ('deta_ds.a');           value = ele_to_s%a%deta_ds
  case ('deta_ds.b');           value = ele_to_s%b%deta_ds
  case ('deta_ds.x');           value = ele_to_s%x%deta_ds
  case ('deta_ds.y');           value = ele_to_s%y%deta_ds
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('e_curl', 'e0_curl')
  orb = orbit
  if (prefix == 'e0_curl') then
    time = orbit%t
  else
    time = particle_rf_time(orb, ele_here, (ele_here%field_calc /= fieldmap$), orb%s-ele_here%s_start)
  endif

  call em_field_derivatives (ele_here, branch%param, orb%s-ele_here%s_start, orb, .false., field, rf_time = time)
  time = particle_rf_time(orb, ele_here, (ele_here%field_calc /= fieldmap$), orb%s-ele_here%s_start)
  dt = bmad_com%d_orb(5) / c_light
  call em_field_calc (ele_here, branch%param, orb%s-ele_here%s_start, orb, .false., field0, rf_time = time-dt)
  call em_field_calc (ele_here, branch%param, orb%s-ele_here%s_start, orb, .false., field1, rf_time = time+dt)

  select case (d_type)
  case ('e_curl.x', 'e0_curl.x');  value = field%dE(3,2) - field%dE(2,3) + (field1%B(1) - field0%B(1)) / (2 * dt)
  case ('e_curl.y', 'e0_curl.y');  value = field%dE(1,3) - field%dE(3,1) + (field1%B(2) - field0%B(2)) / (2 * dt)
  case ('e_curl.z', 'e0_curl.z');  value = field%dE(2,1) - field%dE(1,2) + (field1%B(3) - field0%B(3)) / (2 * dt)
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('e_div', 'e0_div')
  orb = orbit
  if (prefix == 'e0_div') then
    time = orbit%t
  else
    time = particle_rf_time(orb, ele_here, (ele_here%field_calc /= fieldmap$), orb%s-ele_here%s_start)
  endif

  call em_field_derivatives (ele_here, branch%param, orb%s-ele_here%s_start, orb, .false., field, rf_time = time)
  value = field%dE(1,1) + field%dE(2,2) + field%dE(3,3)

case ('e_field', 'e0_field')
  orb = orbit
  if (prefix == 'e0_field') then
    time = orbit%t
  else
    time = particle_rf_time(orb, ele_here, (ele_here%field_calc /= fieldmap$), orb%s-ele_here%s_start)
  endif

  call em_field_calc (ele_here, branch%param, orb%s-ele_here%s_start, orb, .false., field, rf_time = time)
  select case (d_type)
  case ('e_field.x', 'e0_field.x');  value = field%e(1)
  case ('e_field.y', 'e0_field.y');  value = field%e(2)
  case ('e_field.z', 'e0_field.z');  value = field%e(3)
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('e_tot_ref');          value = ele_to_s%value(e_tot$)

case ('energy')
  if (orbit%species == photon$) then
    value = orbit%p0c
  else
    call convert_pc_to(orbit%p0c * (1 + orbit%vec(6)), orbit%species, e_tot = value)
  endif

case ('eta')
  select case (d_type)
  case ('eta.a');            value = ele_to_s%a%eta
  case ('eta.b');            value = ele_to_s%b%eta
  case ('eta.x');            value = ele_to_s%x%eta
  case ('eta.y');            value = ele_to_s%y%eta
  case ('eta.z');            value = ele_to_s%z%eta
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('etap')
  select case (d_type)
  case ('etap.a');           value = ele_to_s%a%etap
  case ('etap.b');           value = ele_to_s%b%etap
  case ('etap.x');           value = ele_to_s%x%etap
  case ('etap.y');           value = ele_to_s%y%etap
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('floor')
  select case (d_type)
  case ('floor.x');          value = ele_to_s%floor%r(1)
  case ('floor.y');          value = ele_to_s%floor%r(2)
  case ('floor.z');          value = ele_to_s%floor%r(3)
  case ('floor.theta');      value = ele_to_s%floor%theta
  case ('floor.phi');        value = ele_to_s%floor%phi
  case ('floor.psi');        value = ele_to_s%floor%psi
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('gamma')
  select case (d_type)
  case ('gamma.a');          value = ele_to_s%a%gamma
  case ('gamma.b');          value = ele_to_s%b%gamma
  case ('gamma.z');          value = ele_to_s%z%gamma
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('intensity')
  select case (d_type)
  case ('intensity');                     value = orbit%field(1)**2 + orbit%field(2)**2
  case ('intensity_x', 'intensity.x');    value = orbit%field(1)**2
  case ('intensity_y', 'intensity.y');    value = orbit%field(2)**2
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('k')
  call c_to_cbar (ele_to_s, cbar)
  f = sqrt(ele_to_s%a%beta/ele_to_s%b%beta) 
  select case (d_type)
  case ('k.11b');            value = cbar(1,1) * f / ele_to_s%gamma_c
  case ('k.12a');            value = cbar(1,2) / (f * ele_to_s%gamma_c)
  case ('k.12b');            value = cbar(1,2) * f / ele_to_s%gamma_c
  case ('k.22a');            value = cbar(2,2) / (f * ele_to_s%gamma_c)
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('momentum');            value = (1 + orbit%vec(6)) * orbit%p0c

case ('orbit')
  if (d_type(7:9) == 'amp' .or. d_type(7:9) == 'nor') &
          call orbit_amplitude_calc (ele_to_s, orbit, amp_a, amp_b, amp_na, amp_nb)
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
    if (present(bad_datum)) bad_datum = .true.
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
  case ('phase.a');           value = ele_to_s%a%phi
  case ('phase_frac.a');      value = modulo2 (ele_to_s%a%phi, pi)
  case ('phase.b');           value = ele_to_s%b%phi
  case ('phase_frac.b');      value = modulo2 (ele_to_s%b%phi, pi)
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('ping_a')
  call c_to_cbar (ele_to_s, cbar)
  select case (d_type)
  case ('ping_a.amp_x');          value = ele_to_s%gamma_c * sqrt(ele_to_s%a%beta)
  case ('ping_a.phase_x');        value = ele_to_s%a%phi
  case ('ping_a.amp_y');          value = sqrt(ele_to_s%b%beta * (cbar(1,2)**2 + cbar(2,2)**2))
  case ('ping_a.phase_y');        value = ele_to_s%a%phi + atan2(cbar(1,2), -cbar(2,2))
  case ('ping_a.amp_sin_rel_y');  value = -sqrt(ele_to_s%b%beta) * cbar(1,2)
  case ('ping_a.amp_cos_rel_y');  value = -sqrt(ele_to_s%b%beta) * cbar(2,2)
  case ('ping_a.amp_sin_y')
    amp = sqrt(ele_to_s%b%beta * (cbar(1,2)**2 + cbar(2,2)**2))
    phase = ele_to_s%a%phi + atan2(cbar(1,2), -cbar(2,2))
    value = amp * sin(phase)
  case ('ping_a.amp_cos_y')
    amp = sqrt(ele_to_s%b%beta * (cbar(1,2)**2 + cbar(2,2)**2))
    phase = ele_to_s%a%phi + atan2(cbar(1,2), -cbar(2,2))
    value = amp * cos(phase)
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('ping_b')
  call c_to_cbar (ele_to_s, cbar)
  select case (d_type)
  case ('ping_b.amp_y');          value = ele_to_s%gamma_c * sqrt(ele_to_s%b%beta)
  case ('ping_b.phase_y');        value = ele_to_s%b%phi
  case ('ping_b.amp_x');          value = sqrt(ele_to_s%a%beta * (cbar(1,2)**2 + cbar(1,1)**2))
  case ('ping_b.phase_x');        value = ele_to_s%b%phi + atan2(cbar(1,2), cbar(1,1))
  case ('ping_b.amp_sin_rel_x');  value = -sqrt(ele_to_s%a%beta) * cbar(1,2)
  case ('ping_b.amp_cos_rel_x');  value = sqrt(ele_to_s%a%beta) * cbar(1,1)
  case ('ping_b.amp_sin_x')
    amp = sqrt(ele_to_s%a%beta * (cbar(1,2)**2 + cbar(1,1)**2))
    phase = ele_to_s%b%phi + atan2(cbar(1,2), cbar(1,1))
    value = amp * sin(phase)
  case ('ping_b.amp_cos_x')
    amp = sqrt(ele_to_s%a%beta * (cbar(1,2)**2 + cbar(1,1)**2))
    phase = ele_to_s%b%phi + atan2(cbar(1,2), cbar(1,1))
    value = amp * cos(phase)
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case ('r')
  i = tao_read_phase_space_index (data_type, 3, .false.)
  j = tao_read_phase_space_index (data_type, 4, .false.)
  if (i == 0 .or. j == 0 .or. len_trim(data_type) /= 4) then
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'BAD DATA_TYPE = "' // trim(data_type)
    return
  endif

  value = ele_to_s%mat6(i,j)

case ('ref_p0c')
  value = orbit%p0c

case ('ref_time')
  value = ele_to_s%ref_time

case ('spin')
  select case (d_type)
  case ('spin.x');        value = orbit%spin(1)
  case ('spin.y');        value = orbit%spin(2)
  case ('spin.z');        value = orbit%spin(3)
  case ('spin.amp');      value = norm2(orbit%spin)
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

  if (.not. bmad_com%spin_tracking_on) then
    err_flag = .true.
    if (present(why_invalid)) why_invalid = 'SPIN TRACKING IS NOT ENABLED.'
  endif

case ('s_position');         value = ele_to_s%s

case ('time', 't');          value = orbit%t

case ('velocity')
  select case (d_type)
  case ('velocity');  value = orbit%beta
  case ('velocity.x');  value = orbit%vec(2) * (1 + orbit%vec(6)) * orbit%beta
  case ('velocity.y');  value = orbit%vec(4) * (1 + orbit%vec(6)) * orbit%beta
  case ('velocity.z');  value = sqrt(1 - (orbit%vec(2) * (1 + orbit%vec(6)))**2 - (orbit%vec(4) * (1 + orbit%vec(6)))**2) * orbit%beta
  case default
    err_flag = .true.
    if (present(bad_datum)) bad_datum = .true.
    if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
  end select

case default
  err_flag = .true.
  if (present(bad_datum)) bad_datum = .true.
  if (present(why_invalid)) why_invalid = 'INVALID DATA_TYPE: ' // quote(data_type)
end select

end function
