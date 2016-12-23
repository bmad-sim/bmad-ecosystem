!+
! Function tao_bmad_parameter_value (dat_name, ele, orbit, err_flag) result (value)
!
! Routine to take a Tao datam name ("beta.a", "orbit.x", etc.), translate to the corresponding Bmad parameter
! and then evaluate it at the given element or orbit position.
!
! Input:
!   dat_name       -- character(*): Data name.
!   ele           -- ele_struct: Lattice element to evaluate the parameter at.
!   orbit         -- coord_struct: Orbit to evaluate the parameter at.
!
! Output:
!   err_flag      -- logical: Set true if dat_name does not have a corresponding Bmad parameter.
!   value         -- real(rp): Parameter value.
!-

function tao_bmad_parameter_value (dat_name, ele, orbit, err_flag) result (value)

use tao_mod, except_dummy => tao_bmad_parameter_value
use measurement_mod

implicit none

type (ele_struct) ele
type (coord_struct) orbit
type (bpm_phase_coupling_struct) bpm_data
type (spin_polar_struct) polar_spin

real(rp) value, cbar(2,2), f, amp_a, amp_b, amp_na, amp_nb

character(*) dat_name
character(40) name, prefix

integer ix

logical err_flag

!

err_flag = .false.
ix = index(dat_name, '.')
if (ix == 0) then
  prefix = dat_name
else
  prefix = dat_name(1:ix)
endif

!

select case (prefix)

case ('alpha.')
  select case (dat_name)
  case ('alpha.a');          value = ele%a%alpha
  case ('alpha.b');          value = ele%b%alpha
  case ('alpha.z');          value = ele%z%alpha
  case default;              err_flag = .true.
  end select

case ('beta.')
  select case (dat_name)
  case ('beta.z');           value = ele%z%beta
  case ('beta.a');           value = ele%a%beta
  case ('beta.b');           value = ele%b%beta
  case default;              err_flag = .true.
  end select

case ('bpm_cbar.')
  call to_phase_and_coupling_reading (ele, bpm_data, err_flag)
  select case (dat_name)
  case ('bpm_cbar.22a');     value = bpm_data%cbar22_a
  case ('bpm_cbar.12a');     value = bpm_data%cbar12_a
  case ('bpm_cbar.11b');     value = bpm_data%cbar11_b
  case ('bpm_cbar.12b');     value = bpm_data%cbar12_b
  case default
  end select

case ('bpm_eta.')
  select case (dat_name)
  case ('bpm_eta.x');        call to_eta_reading ([ele%x%eta, ele%y%eta], ele, x_plane$, value, err_flag)
  case ('bpm_eta.y');        call to_eta_reading ([ele%x%eta, ele%y%eta], ele, y_plane$, value, err_flag)
  case default;              err_flag = .true.
  end select

case ('bpm_orbit.')
  select case (dat_name)
  case ('bpm_orbit.x');      call to_orbit_reading (orbit, ele, x_plane$, value, err_flag)
  case ('bpm_orbit.y');      call to_orbit_reading (orbit, ele, y_plane$, value, err_flag)
  case default;              err_flag = .true.
  end select

case ('bpm_phase.')
  call to_phase_and_coupling_reading (ele, bpm_data, err_flag)
  if (err_flag) return
  select case (dat_name)
  case ('bpm_phase.a');      value = bpm_data%phi_a
  case ('bpm_phase.b');      value = bpm_data%phi_b
  case default;              err_flag = .true.
  end select

case ('bpm_k.')
  call to_phase_and_coupling_reading (ele, bpm_data, err_flag)
  select case (dat_name)
  case ('bpm_k.22a');        value = bpm_data%k_22a
  case ('bpm_k.12a');        value = bpm_data%k_12a
  case ('bpm_k.11b');        value = bpm_data%k_11b
  case ('bpm_k.12b');        value = bpm_data%k_12b
  case default;              err_flag = .true.
  end select

case ('c_mat.')
  select case (dat_name)
  case ('c_mat.11');         value = ele%c_mat(1,1)
  case ('c_mat.12');         value = ele%c_mat(1,2)
  case ('c_mat.21');         value = ele%c_mat(2,1)
  case ('c_mat.22');         value = ele%c_mat(2,2)
  case default;              err_flag = .true.
  end select

case ('cbar.')
  call c_to_cbar (ele, cbar)
  select case (dat_name)
  case ('cbar.11');           value = cbar(1,1)
  case ('cbar.12');           value = cbar(1,2)
  case ('cbar.21');           value = cbar(2,1)
  case ('cbar.22');           value = cbar(2,2)
  case default
  end select

case ('e_tot');               value = (1 + orbit%vec(6)) * orbit%p0c / orbit%beta

case ('eta.')
  select case (dat_name)
  case ('eta.a');            value = ele%a%eta
  case ('eta.b');            value = ele%b%eta
  case ('eta.x');            value = ele%x%eta
  case ('eta.y');            value = ele%y%eta
  case ('eta.z');            value = ele%z%eta
  case default;              err_flag = .true.
  end select

case ('etap.')
  select case (dat_name)
  case ('etap.a');           value = ele%a%etap
  case ('etap.b');           value = ele%b%etap
  case ('etap.x');           value = ele%x%etap
  case ('etap.y');           value = ele%y%etap
  case default;              err_flag = .true.
  end select

case ('floor.')
  select case (dat_name)
  case ('floor.x');          value = ele%floor%r(1)
  case ('floor.y');          value = ele%floor%r(2)
  case ('floor.z');          value = ele%floor%r(3)
  case ('floor.theta');      value = ele%floor%theta
  case ('floor.phi');        value = ele%floor%phi
  case ('floor.psi');        value = ele%floor%psi
  case default;              err_flag = .true.
  end select

case ('gamma.')
  select case (dat_name)
  case ('gamma.a');          value = ele%a%gamma
  case ('gamma.b');          value = ele%b%gamma
  case ('gamma.z');          value = ele%z%gamma
  case default;              err_flag = .true.
  end select

case ('k.')
  call c_to_cbar (ele, cbar)
  f = sqrt(ele%a%beta/ele%b%beta) 
  select case (dat_name)
  case ('k.11b');            value = cbar(1,1) * f / ele%gamma_c
  case ('k.12a');            value = cbar(1,2) / (f * ele%gamma_c)
  case ('k.12b');            value = cbar(1,2) * f / ele%gamma_c
  case ('k.22a');            value = cbar(2,2) / (f * ele%gamma_c)
  case default;              err_flag = .true.
  end select

case ('momentum');            value = (1 + orbit%vec(6)) * orbit%p0c

case ('orbit.')
  if (dat_name(7:9) == 'amp' .or. dat_name(7:9) == 'nor') &
          call orbit_amplitude_calc (ele, orbit, amp_a, amp_b, amp_na, amp_nb)
  select case (dat_name)
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
  case default;               err_flag = .true.
  end select

case ('pc');                  value = (1 + orbit%vec(6)) * orbit%p0c

case ('phase.', 'phase_frac.')
  select case (dat_name)
  case ('phase.a');           value = ele%a%phi
  case ('phase_frac.a');      value = modulo2 (ele%a%phi, pi)
  case ('phase.b');           value = ele%a%phi
  case ('phase_frac.b');      value = modulo2 (ele%b%phi, pi)
  case default;               err_flag = .true.
  end select

case ('ref_time');            value = ele%ref_time

case ('spin.')
  select case (dat_name)
  case ('spin.x');        value = orbit%spin(1)
  case ('spin.y');        value = orbit%spin(2)
  case ('spin.z');        value = orbit%spin(3)
  case ('spin.theta', 'spin.phi', 'spin.amplitude')
    polar_spin = vec_to_polar(orbit%spin)
    select case (dat_name)
    case ('spin.theta');     value = polar_spin%theta
    case ('spin.phi');       value = polar_spin%phi
    case ('spin.amp');       value = polar_spin%polarization
    end select
  case default;              err_flag = .true.
  end select

case ('s_position');         value = ele%s

case ('time');               value = orbit%t

case default;                err_flag = .true.

end select


end function
