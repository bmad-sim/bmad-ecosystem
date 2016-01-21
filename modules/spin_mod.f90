module spin_mod

use bmad_struct
use bmad_interface

! This includes the phase of the spinor.
! Polarization is not 1 when the spin_polar struct represents an ensamble of spins.

type spin_polar_struct
  real(rp) :: polarization = 1
  real(rp) :: theta = 0
  real(rp) :: phi   = 0
  real(rp) :: xi    = 0
end type

! Pauli Matrices
type pauli_struct
  complex(rp) sigma(2,2)
end type

type (pauli_struct) pauli(0:3)

logical, private :: init_pauli_vector = .true. ! Does pauli vector needs to be set up?
private trapzd_omega, sbend_omega_func, quad_etc_omega_func

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine Initialize_pauli_vector ()
!
! This subroutine is not intended for public use.
!
! initialize pauli vector, if needed.
!
! If init_pauli_vector = T then pauli vector will be set up.
!-

subroutine initialize_pauli_vector ()

implicit none

!

if (.not. init_pauli_vector) return

pauli(0)%sigma(1,1) = ( 1.0,  0.0)
pauli(0)%sigma(2,1) = ( 0.0,  0.0)
pauli(0)%sigma(1,2) = ( 0.0,  0.0)
pauli(0)%sigma(2,2) = ( 1.0,  0.0)

pauli(1)%sigma(1,1) = ( 0.0,  0.0)
pauli(1)%sigma(2,1) = ( 1.0,  0.0)
pauli(1)%sigma(1,2) = ( 1.0,  0.0)
pauli(1)%sigma(2,2) = ( 0.0,  0.0)

pauli(2)%sigma(1,1) = ( 0.0,  0.0)
pauli(2)%sigma(2,1) = ( 0.0,  1.0)
pauli(2)%sigma(1,2) = ( 0.0, -1.0)
pauli(2)%sigma(2,2) = ( 0.0,  0.0)

pauli(3)%sigma(1,1) = ( 1.0,  0.0)
pauli(3)%sigma(2,1) = ( 0.0,  0.0)
pauli(3)%sigma(1,2) = ( 0.0,  0.0)
pauli(3)%sigma(2,2) = (-1.0,  0.0)

init_pauli_vector = .false.

end subroutine initialize_pauli_vector

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function spinor_to_polar (spinor) result (polar)
!
! Converts a spinor into a spin polar vector of unit length
!
! Modules needed:
!   use spin_mod
!
! Input:
!   spinor(2)  -- complex(rp): Spinor
!
! Output:
!   polar      -- Spin_polar_struct: The resultant Unitary Vector in polar coordinates
!-

function spinor_to_polar (spinor) result (polar)

implicit none

type (spin_polar_struct) ::  polar

complex(rp) spinor(2)
real(rp) temp(2)

character(20) :: r_name = "spinor_to_polar"

!

temp(1) = atan2 (imag(spinor(1)), real(spinor(1)))
temp(2) = atan2 (imag(spinor(2)), real(spinor(2)))
polar%xi = temp(1)
polar%phi = temp(2) - temp(1)

temp=abs(spinor)
polar%theta = 2 * atan2(temp(2), temp(1))
! no sqrt here! spinor scales with sqrt(r)
polar%polarization = dot_product(temp, temp)


end function spinor_to_polar

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function polar_to_vec (polar) result (vec)
!
! Comverts a spinor in polar coordinates to a spin vector. This will ignore the
! spinor phase.
!
! Modules needed:
!   use spin_mod
!
! Input:
!   polar         -- Spin_polar_struct
!
! Output:
!   vec(3)        -- Real(3)
!-

function polar_to_vec (polar) result (vec)

implicit none

type (spin_polar_struct) polar

real(rp) vec(3)

vec(1) = polar%polarization * sin(polar%theta) * cos(polar%phi)
vec(2) = polar%polarization * sin(polar%theta) * sin(polar%phi)
vec(3) = polar%polarization * cos(polar%theta)

end function polar_to_vec

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function polar_to_spinor (polar) result (spin)
!
! Converts a spin vector in polar coordinates to a spinor
!
! Modules needed:
!   use spin_mod
!
! Input:
!   polar     -- spin_polar_struct: includes polar phase
!
! Output:
!   spin(2)   -- complex(rp): the particle spin
!-

function polar_to_spinor (polar) result (spin)

implicit none

type (spin_polar_struct) polar
complex(rp) :: spin(2)

!

spin(1) = sqrt(polar%polarization) * Exp(i_imaginary * polar%xi) * cos(polar%theta / 2.0d0)
spin(2) = sqrt(polar%polarization) * Exp(i_imaginary * (polar%xi+polar%phi)) * sin(polar%theta / 2.0d0)

end function polar_to_spinor

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function vec_to_polar (vec, phase) result (polar)
!
! Converts a spin vector to a spin polar
!
! Modules needed:
!   use spin_mod
!
! Input:
!   vec(3)   -- real(rp): unitary spin vector
!   phase    -- real(rp), Optional: Phase of the spinor, if not given then
!                                   set to zero
!
! Output:
!   polar    -- spin_polar_struct:
!-

function vec_to_polar (vec, phase) result (polar)

implicit none

type (spin_polar_struct) :: polar

real(rp) vec(3)
real(rp), optional :: phase

!

polar%xi = real_option (0.0d0, phase)
polar%theta = atan2 (sqrt(vec(1)**2 + vec(2)**2), vec(3))
polar%phi = atan2(vec(2), vec(1))
polar%polarization = sqrt(dot_product(vec, vec))

end function vec_to_polar

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function spinor_to_vec (spinor) result (vec)
!
! Converts a spinor to a spin cartesian vector
!
! Modules needed:
!   use spin_mod
!
! Input:
!   spinor  -- complex(rp): Spinor
!
! Output
!   vec(3) -- Real(rp): spin vector in cartesian coordinates
!-

function spinor_to_vec (spinor) result (vec)

implicit none

complex(rp) spinor(2)
real(rp) vec(3)

!

! vec = conjg(spinor) * pauli(i)%sigma * spinor done explicitly
vec(1) = 2 * (real(spinor(1))*real(spinor(2)) + aimag(spinor(1))*aimag(spinor(2)))
vec(2) = 2 * (real(spinor(1))*aimag(spinor(2)) - aimag(spinor(1))*real(spinor(2)))
vec(3) = real(spinor(1))**2 + aimag(spinor(1))**2 - real(spinor(2))**2 - aimag(spinor(2))**2

end function spinor_to_vec

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function vec_to_spinor (vec, phase) result (spinor)
!
! Converts a spin cartesian vector to a spinor.
!
! Modules needed:
!   use spin_mod
!
! Input:
!   vec(3)   -- real(rp): Spin vector in cartesian coordinates
!   phase    -- real(rp)(Optional): Phase of the spinor, if not given then
!                                   set to zero
!
! Output:
!   spinor(2)-- complex(rp): Spinor.
!-

function vec_to_spinor (vec, phase) result (spinor)

implicit none

type (spin_polar_struct) :: polar

complex(rp) :: spinor(2)
real(rp) vec(3)
real(rp), optional :: phase

!

polar = vec_to_polar(vec, phase)
spinor = polar_to_spinor(polar)

end function vec_to_spinor

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! function angle_between_polars (polar1, polar2) result (angle)
!
! Finds the angle between two polar vectors.
! Note: This function is currently not being used by anything.
!
! Modules needed:
!   use spin_mod
!
! Input:
!   polar1    -- (spin_polar_struct)
!   polar2    -- (spin_polar_struct)
!
! Output:
!   angle     -- Real(rp): Angle between the polar vectors
!-

function angle_between_polars (polar1, polar2) result (angle)

implicit none

type (spin_polar_struct), intent(in) :: polar1, polar2

real(rp) :: angle, arg
real(rp) :: vec1(3), vec2(3)

! Round-off can make |arg| > 1 so need to check this.

vec1 = polar_to_vec (polar1) 
vec2 = polar_to_vec (polar2)

arg = dot_product(vec1,vec2) / (sqrt(dot_product(vec1, vec1) * dot_product(vec2,vec2)))
if (arg >= 1) then
  angle = 0
elseif (arg <= -1) then
  angle = pi
else
  angle = acos(arg)
endif

end function angle_between_polars

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine rotate_spinor_given_field (orbit, ele, BL, EL)
!
! Routine to rotate a spinor given the integrated magnetic and/or electric field strengths.
!
! Modules needed:
!   use spin_mod
!
! Input:
!   orbit       -- coord_struct: Initial orbit.
!   ele         -- ele_struct: Element being tracked through.
!   BL(3)       -- real(rp), optional: Integrated field strength. Assumed zero if not present.
!   EL(3)       -- real(rp), optional: Integrated field strength. Assumed zero if not present.
!
! Output:
!   orbit       -- coord_struct: Orbit with rotated spin
!-

subroutine rotate_spinor_given_field (orbit, ele, BL, EL)

implicit none

type (ele_struct) ele
type (coord_struct) orbit
type (em_field_struct) field

real(rp), optional :: BL(3), EL(3)
real(rp)  omega(3)

!

if (present(BL)) then
  field%B = BL
else
  field%B = 0
endif

if (present(EL)) then
  field%E = EL
else
  field%E = 0
endif

omega = spin_omega (field, orbit, ele)
call rotate_spinor (omega, orbit%spin)

end subroutine rotate_spinor_given_field

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine rotate_spinor (rot_vec, spin)
!
! Routine to rotate a spinor.
!
! Modules needed:
!   use spin_mod
!
! Input:
!   rot_vec(3)  -- real(rp): Rotation axis. Magnitude of rot_vec is the rotation angle.
!   spin(2)     -- complex(rp): Initial coords.
!
! Output:
!   spin(2)     -- complex(rp): Final coords.
!-

subroutine rotate_spinor (rot_vec, spin)

implicit none

complex(rp) :: spin(2), mat(2,2)
real(rp) :: rot_vec(3), angle, n_vec(3), c, s

!

angle = norm2(rot_vec)
if (angle == 0) return

c = cos(angle/2)
s = sin(angle/2)

n_vec = rot_vec * (s / angle)

mat(1,:) = [cmplx(c, -n_vec(3), rp),        cmplx(-n_vec(2), -n_vec(1), rp)]
mat(2,:) = [cmplx(n_vec(2), -n_vec(1), rp), cmplx(c, n_vec(3), rp)]
spin = matmul(mat, spin)

end subroutine rotate_spinor

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine quaternion_track (a_quat, spin)
!
! Transports a spinor using the quaternion a_quat
!
! Modules needed:
!   use spin_mod
!
! Input:
!   a_quat(0:3) -- real(rp): Euler four-vector (Quaternion)
!   spin(2)     -- complex(rp): Incoming spinor
!
! Output:
!   spin(2)    -- complex(rp): Resultant spinor
!-

subroutine quaternion_track (a_quat, spin)

implicit none

complex(rp), intent(inout) :: spin(2)

real(rp), intent(in) :: a_quat(0:3)

complex(rp) a_mat(2,2) ! The matrix associated with a_quat.

!

if (init_pauli_vector) call initialize_pauli_vector

a_mat = a_quat(0)*pauli(0)%sigma + i_imaginary * &
          (a_quat(1) * pauli(1)%sigma + a_quat(2) * pauli(2)%sigma + a_quat(3) * pauli(3)%sigma)

spin = matmul (a_mat, spin)

end subroutine quaternion_track

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function calc_rotation_quaternion (n_vec, angle) result (a)
!
! Calculates the quaternion for a rotation of a vector by an angle about (nx,ny,nz).
! (nx,ny,nz) has to be a unit vector, i.e. nx^2 + y^2 + nz^2 = 1.
!
! Modules needed:
!   use spin_mod
!
! Input:
!   n_vec(3)   -- Real(rp): Unit rotation axis vector.
!   angle      -- Real(rp): Rotation angle
!
! Output:
!   a(0:3)       -- Real(rp): Resultant quaternion
!-

function calc_rotation_quaternion (n_vec, angle) result (a)

real(rp) , intent(in) :: n_vec(3), angle
real(rp) :: a(0:3)

real(rp) half_angle, s

! 

half_angle = angle/2.
s = sin(half_angle)
a = [cos(half_angle), -n_vec(1)*s, -n_vec(2)*s, -n_vec(3)*s]

end function calc_rotation_quaternion

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! subroutine track1_spin (start_orb, ele, param, end_orb)
!
! Particle spin tracking through a single element.
!
! Typically this routine should not be directly called. 
! Instead, use track1 which calls this routine.
!
! Modules needed:
!   use spin_mod
!
! Input :
!   start_orb  -- Coord_struct: Starting coords.
!   ele        -- Ele_struct: Element to track through.
!   param      -- lat_param_struct: Beam parameters.
!   end_orb    -- Coord_struct: Ending coords.
!     %vec          -- Ending particle position needed for bmad_standard spin tracking.
!
! Output:
!   end_orb    -- Coord_struct: Ending coords.
!      %spin(2)   -- complex(rp): Ending spinor
!-

subroutine track1_spin (start_orb, ele, param, end_orb)

use ptc_spin, rename_dummy => dp, rename2_dummy => twopi
use ptc_interface_mod

implicit none

type (coord_struct) :: start_orb, end_orb, temp_orb
type (ele_struct) :: ele
type (lat_param_struct) :: param

integer method
character(*), parameter :: r_name = 'track1_spin'
logical err

! Use bmad_standard if spin_tracking_method = tracking$ and particle tracking is not using an integration method.

if (start_orb%species == photon$) return

method = ele%spin_tracking_method
if (method == tracking$) then
  select case (ele%tracking_method)
  case (boris$, runge_kutta$, symp_lie_ptc$)
    return ! Spin tracking is done at the same time orbital tracking is done
  case default
    method = bmad_standard$
  end select
endif

!

select case (method)
case (bmad_standard$)
  if (bmad_com%electric_dipole_moment /= 0) then
    call out_io (s_fatal$, r_name, &
          'TRACKING WITH AN ELECTRIC DIPOLE MOMENT NOT YET DEVELOPED FOR BMAD_STANDARD TRACKING')
    if (global_com%exit_on_error) call err_exit
  endif
  call track1_spin_bmad (start_orb, ele, param, end_orb)

case (custom$)
  call track1_spin_custom (start_orb, ele, param, end_orb, err)

! Notice that PTC spin tracking is only done here only when the (orbital) tracking_method is *not* symp_lie_ptc
case (symp_lie_ptc$)
  if (bmad_com%electric_dipole_moment /= 0) then
    call out_io (s_fatal$, r_name, &
          'TRACKING WITH AN ELECTRIC DIPOLE MOMENT NOT YET DEVELOPED FOR SYMP_LIE_PTC TRACKING')
    if (global_com%exit_on_error) call err_exit
  endif
  call track1_symp_lie_ptc (start_orb, ele, param, temp_orb)

case default
  call out_io (s_fatal$, r_name, 'BAD SPIN_TRACKING_METHOD: ' // spin_tracking_method_name(ele%spin_tracking_method), &
                                 'FOR ELEMENT: ', ele%name)
  if (global_com%exit_on_error) call err_exit
end select

end subroutine track1_spin

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! subroutine track1_spin_bmad (start_orb, ele, param, end_orb)
!
! Bmad_standard particle spin tracking through a single element.
!
! Note: spin tracking through a patch element is handled in track_a_patch since
! this is needed by runge_kutta tracking.
!
! Modules needed:
!   use spin_mod
!
! Input :
!   start_orb  -- Coord_struct: Starting coords.
!   ele        -- Ele_struct: Element to track through.
!   param      -- lat_param_struct: Beam parameters.
!   end_orb    -- Coord_struct: Ending coords.
!
! Output:
!   end_orb    -- Coord_struct:
!     %spin(2)   -- complex(rp): Ending spinor
!-

subroutine track1_spin_bmad (start_orb, ele, param, end_orb)

use ptc_spin, rename_dummy => dp, rename2_dummy => twopi
use ptc_interface_mod

implicit none

type (coord_struct) :: start_orb
type (coord_struct) :: temp_start, temp_middle, temp_end, end_orb
type (ele_struct) :: ele
type (lat_param_struct) :: param

real(rp) a_quat(0:3) ! quaternion four-vector
real(rp) omega1, xi, gamma0, gammaf, v, x, u, spline_x(0:3), spline_y(0:3)
real(rp) alpha, phase, cos_phi, gradient, pc_start, pc_end, omega(3)
real(rp) e_start, e_end, g_ratio, edge_length, beta_start, beta_end
real(rp) anomalous_moment, m_particle, sign_k, abs_a, coef
real(rp) vec0(6), ks, kss, length, kl, kl2, c, s, m21, m22, m23, m24

integer key

character(*), parameter :: r_name = 'track1_spin_bmad'

!

if (ele%key == patch$) return  ! Spin tracking handled by track_a_patch for patch elements.

m_particle = mass_of(start_orb%species)
anomalous_moment = anomalous_moment_of(start_orb%species)

key = ele%key
if (.not. ele%is_on .and. key /= lcavity$) key = drift$

! 

temp_start = start_orb
call offset_particle (ele, param, set$, temp_start, .true., .true., .true., set_spin = .true.)
call apply_element_edge_kick (temp_start, 0.0_rp, temp_start%t, ele, ele, param, first_track_edge$, .true.)

temp_end   = end_orb
call offset_particle (ele, param, set$, temp_end, .true., .false., .false., .true., ele%value(l$))

temp_end%spin = temp_start%spin

select case (key)

!-----------------------------------------------
! quadrupole, etc.

case (quadrupole$, sextupole$, octupole$, sol_quad$)

  call spline_fit_orbit (temp_start, temp_end, spline_x, spline_y)
  omega = trapzd_omega (ele, quad_etc_omega_func, temp_start, spline_x, spline_y)
  call rotate_spinor(omega, temp_end%spin)

!-----------------------------------------------
! sbend
! does not take k1, k2 (quadrupole and sextupole terms) into account

case (sbend$)

  call spline_fit_orbit (temp_start, temp_end, spline_x, spline_y)
  omega = trapzd_omega (ele, sbend_omega_func, temp_start, spline_x, spline_y) + [0.0_rp, ele%value(g$)*ele%value(l$), 0.0_rp]
  call rotate_spinor(omega, temp_end%spin)

!-----------------------------------------------
! solenoid

case (solenoid$)

  vec0 = temp_start%vec
  ks = ele%value(ks$) / (1+vec0(6))
  kss = ks / 2
  length = ele%value(l$)

  if (abs(length * kss) < 1d-10) then
    kl = kss * length 

    m21 = -kl*kss/2
    m22 = 1
    m23 = -kl**2 * kss/3
    m24 = kl/2

  else
    c = cos(ks * length)
    s = sin(ks * length)
    m21 = (c - 1) / (4 * length)
    m22 = (s / (ks * length) + 1) / 2
    m23 = kss * (m22 - 1)
    m24 = -m21 / kss
  endif

  temp_end%vec(2) =  m21 * vec0(1) + m22 * vec0(2) + m23 * vec0(3) + m24 * vec0(4)
  temp_end%vec(4) = -m23 * vec0(1) - m24 * vec0(2) + m21 * vec0(3) + m22 * vec0(4)
  call rotate_spinor_given_field (temp_end, ele, [0.0_rp, 0.0_rp, length*ele%value(bs_field$)])

!-----------------------------------------------
! LCavity
!
! Simulates the cavity edge field kicks as electrostatic quadrupoles
! since the quaternions for these have already been found.
!
! Uses the fringe field as calculated by Hartman and Rosenzweig

case (lcavity$)

  temp_middle%vec = (temp_start%vec + temp_end%vec) / 2 ! rough estimate of particle coordinates in the element
  a_quat = 0

  ! For now, just set to one
  g_ratio = 1

  gamma0 = ((1+temp_middle%vec(6)) * ele%value(E_TOT$)) / m_particle

  if (ele%value(E_TOT_START$) == 0) then
    call out_io (s_fatal$, r_name, 'E_TOT_START IS 0 FOR A LCAVITY!')
    if (global_com%exit_on_error) call err_exit
  endif

  phase = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_err$) - &
                      temp_end%vec(5) * ele%value(rf_frequency$) / c_light)
  cos_phi = cos(phase)
  gradient = (ele%value(gradient$) + ele%value(gradient_err$)) * cos_phi
  if (.not. ele%is_on) gradient = 0

  gradient = gradient + gradient_shift_sr_wake(ele, param)

  !

  if (gradient /= 0) then
    pc_start = ele%value(p0c_start$) * (1+temp_middle%vec(6))
    call convert_pc_to (pc_start, start_orb%species, E_tot = e_start, beta = beta_start)
    e_end = e_start + gradient * ele%value(l$)
    gammaf = gamma0 * (e_end / e_start)
    call convert_total_energy_to (e_end, start_orb%species, pc = pc_end, beta = beta_end)

    ! The edge field length of a cavity is about 1 quarter wavelength

    edge_length = (c_light * beta_start / ele%value(rf_frequency$)) / 4.0

    call lcav_edge_track (pc_start,  gradient, gamma0, anomalous_moment, edge_length, a_quat)
    call lcav_edge_track (pc_end,   -gradient, gammaf, anomalous_moment, edge_length, a_quat)

    ! Questionable equation for a_quat(0)
    a_quat(0) = sqrt(1.0 - (a_quat(1)**2 + a_quat(2)**2 + a_quat(3)**2))

  endif

  if (key /= lcavity$ .or. gradient /= 0) then
    ! Negative sign due to fact that original quaternion_track routine used a left handed rule for the sign of rotations.
    abs_a = sqrt(a_quat(1)**2 + a_quat(2)**2 + a_quat(3)**2 + a_quat(0)**2)
    a_quat = [a_quat(0), -a_quat(1), -a_quat(2), -a_quat(3)] / abs_a
    call quaternion_track (a_quat, temp_end%spin)
  endif

end select

! 

call apply_element_edge_kick (temp_end, ele%value(l$), temp_start%t, ele, ele, param, second_track_edge$, .true.)
call offset_particle (ele, param, unset$, temp_end, .true., .true., .true., set_spin = .true.)
end_orb%spin = temp_end%spin

!-------------------------------------------------------------------------
contains

subroutine add_to_quaternion (quat, expn, coef)

implicit none

type (taylor_term_struct), pointer :: map(:)

real(rp) quat, coef, qq

integer j, expn(6)

!

qq = coef

do j = 1, 6
  qq = qq * temp_middle%vec(j)**expn(j)
enddo

quat = quat + qq

end subroutine add_to_quaternion

!--------------------------------------------------------------------------
! contains

subroutine lcav_edge_track (pc, grad, gam, anomalous_moment, edge_length, a_quat)

real(rp) pc, grad, gam, anomalous_moment, edge_length, k_el, k_el_tilde, omega_el, a_quat(4)

! Is this correct? e_mass is in GeV and not eV!

k_el = abs(grad / (2 * pc))
omega_el = sqrt(k_el)
k_el_tilde = (e_charge * k_el * (1 + anomalous_moment + (anomalous_moment*gam))) / (omega_el * e_mass * c_light**2 * (1 + gam))

! Focusing kick

if (grad > 0) then
  call add_to_quaternion (a_quat(1), [0, 0, 1, 0, 0, 0], -(k_el_tilde/2.0) * sin (omega_el * edge_length))
  call add_to_quaternion (a_quat(1), [0, 0, 0, 1, 0, 0], -(k_el_tilde/omega_el) * (sin (omega_el * edge_length / 2.0))**2)
  call add_to_quaternion (a_quat(2), [0, 0, 1, 0, 0, 0], -(k_el_tilde/2.0) * sin (omega_el * edge_length))
  call add_to_quaternion (a_quat(2), [0, 0, 0, 1, 0, 0], -(k_el_tilde/omega_el) * (sin (omega_el * edge_length / 2.0))**2)

! Defocus kick

else
  call add_to_quaternion (a_quat(1), [0, 0, 1, 0, 0, 0], (k_el_tilde/2.0) * sinh (omega_el * edge_length))
  call add_to_quaternion (a_quat(1), [0, 0, 0, 1, 0, 0], (k_el_tilde/omega_el) * (sinh (omega_el * edge_length / 2.0))**2)
  call add_to_quaternion (a_quat(2), [0, 0, 1, 0, 0, 0], (k_el_tilde/2.0) * sinh (omega_el * edge_length))
  call add_to_quaternion (a_quat(2), [0, 0, 0, 1, 0, 0], (k_el_tilde/omega_el) * (sinh (omega_el * edge_length / 2.0))**2)
endif

end subroutine lcav_edge_track

end subroutine track1_spin_bmad

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine spline_fit_orbit (start_orb, end_orb, spline_x, spline_y)

implicit none

type (coord_struct) start_orb, end_orb
real(rp) spline_x(0:3), spline_y(0:3)
real(rp) ds, relp, alpha, beta

!

ds = end_orb%s - start_orb%s
relp = 1 + start_orb%vec(6)

! X

spline_x(0) = start_orb%vec(1)
spline_x(1) = start_orb%vec(2) / relp 

alpha = end_orb%vec(1) - spline_x(0) - spline_x(1) * ds
beta = end_orb%vec(2) / relp - spline_x(1)

spline_x(2) = 3 * alpha / ds**2 - beta / ds
spline_x(3) = beta / ds**2 - 2 * alpha / ds**3

! Y

spline_y(0) = start_orb%vec(3)
spline_y(1) = start_orb%vec(4) / relp 

alpha = end_orb%vec(3) - spline_y(0) - spline_y(1) * ds
beta = end_orb%vec(4) / relp - spline_y(1)

spline_y(2) = 3 * alpha / ds**2 - beta / ds
spline_y(3) = beta / ds**2 - 2 * alpha / ds**3


end subroutine spline_fit_orbit

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

function trapzd_omega (ele, omega_func, orbit, spline_x, spline_y) result (omega)

use nr, only: polint

implicit none

interface func
  function omega_func(s, spline_x, spline_y, orbit, ele) result (om)
    import
    implicit none
    type (coord_struct) orbit
    type (ele_struct) ele
    real(rp) s, spline_x(0:3), spline_y(0:3), om(3)
  end function
end interface

type q_array_struct
  real(rp) h
  real(rp) omega(3)
end type

integer, parameter ::  j_max = 5

type (q_array_struct) q_array(j_max)
type (ele_struct) ele
type (coord_struct) orbit, orb

real(rp) s1, del_s, s, spline_x(0:3), spline_y(0:3), om(3), omega(3)
real(rp) dint, eps
real(rp), parameter :: eps_rel = 1d-5, eps_abs = 1d-8

integer j, k, n, n_pts

!

s1 = ele%value(l$)

q_array(1)%h = 1
q_array(1)%omega = s1 * (omega_func(0.0_rp, spline_x, spline_y, orbit, ele) + omega_func(s1, spline_x, spline_y, orbit, ele)) / 2

do j = 2, j_max
  ! This is trapzd from NR
  n_pts = 2**(j-2)
  del_s = s1 / n_pts
  omega = 0

  do n = 1, n_pts
    s = del_s * (n - 0.5_rp)
    omega = omega + omega_func(s, spline_x, spline_y, orbit, ele)
  enddo

  q_array(j)%h = q_array(j-1)%h / 4
  q_array(j)%omega = (q_array(j-1)%omega + omega * del_s) / 2

  eps = eps_abs + eps_rel * sum(abs(q_array(j)%omega))

  do k = 1, 3
    call polint (q_array(1:j)%h, q_array(1:j)%omega(k), 0.0_rp, omega(k), dint)
    if (abs(dint) > eps .and. j < j_max) exit ! Failed test. Note: Last loop with j = j_max -> no test.
    if (k == 3) return                        ! Passed all tests or last loop
  enddo

enddo

end function trapzd_omega

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

function sbend_omega_func (s, spline_x, spline_y, orbit, ele) result (omega) 

implicit none

type (coord_struct) orbit, orb
type (ele_struct) ele
type (em_field_struct) field

real(rp) s, spline_x(0:3), spline_y(0:3), omega(3), B(3)
real(rp) x, y, s1

!

field%b = [0.0_rp, ele%value(b_field$) + ele%value(b_field_err$), 0.0_rp]
field%e = 0

if (ele%value(k1$) /= 0 .or. ele%value(k2$) /= 0) then
  x = spline_x(0) + spline_x(1) * s + spline_x(2) * s**2 + spline_x(3) * s**3
  y = spline_y(0) + spline_y(1) * s + spline_y(2) * s**2 + spline_y(3) * s**3
endif

if (ele%value(k1$) /= 0) field%b = field%b + ele%value(b1_gradient$) * [y, x, 0.0_rp]
if (ele%value(k2$) /= 0) field%b = field%b + ele%value(b2_gradient$) * [x*y, x*x-y*y, 0.0_rp]

!

orb = orbit
orb%vec(2) = (1 + orb%vec(6)) * (spline_x(1) + spline_x(2) * s + spline_x(3) * s**2)
orb%vec(4) = (1 + orb%vec(6)) * (spline_y(1) + spline_y(2) * s + spline_y(3) * s**2)

omega = spin_omega (field, orb, ele)

end function sbend_omega_func

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

function quad_etc_omega_func (s, spline_x, spline_y, orbit, ele) result (omega) 

implicit none

type (coord_struct) orbit, orb
type (ele_struct) ele
type (em_field_struct) field

real(rp) s, spline_x(0:3), spline_y(0:3), omega(3), B(3)
real(rp) x, y, s1

!

field%e = 0

x = spline_x(0) + spline_x(1) * s + spline_x(2) * s**2 + spline_x(3) * s**3
y = spline_y(0) + spline_y(1) * s + spline_y(2) * s**2 + spline_y(3) * s**3

select case (ele%key)
case (quadrupole$);   field%b = ele%value(b1_gradient$) * [y, x, 0.0_rp]
case (sextupole$);    field%b = ele%value(b2_gradient$) * [x*y, (x*x-y*y)/2, 0.0_rp]
case (octupole$);     field%b = ele%value(b3_gradient$) * [(3*x*x*y - y**3)/6, (x**3 - 3*x*y*y)/6, 0.0_rp]
case (sol_quad$);     field%b = ele%value(b1_gradient$) * [y, x, 0.0_rp] + [0.0_rp, 0.0_rp, ele%value(bs_field$)]
end select

!

orb = orbit
orb%vec(2) = (1 + orb%vec(6)) * (spline_x(1) + spline_x(2) * s + spline_x(3) * s**2)
orb%vec(4) = (1 + orb%vec(6)) * (spline_y(1) + spline_y(2) * s + spline_y(3) * s**2)

omega = spin_omega (field, orb, ele)

end function quad_etc_omega_func

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function spin_omega (field, coord, ele), result (omega)
!
! Return the modified T-BMT spin omega vector.
!
! Modules needed:
!   use spin_mod
!   use em_field_mod
!
! Input:
!   field      -- em_field_struct: E and B fields
!   coord      -- coord_struct: particle momentum
!   ele        -- ele_struct: element evauluated in
!      %value(E_TOT$) -- reaL(rp): needed to find momentum
!
! Output:
!   omega(3)   -- real(rp): Omega_TBMT/v_z in cartesian coordinates
!-

function spin_omega (field, coord, ele) result (omega)

implicit none

type (em_field_struct) :: field
type (coord_struct) :: coord
type (ele_struct) :: ele

real(rp) omega(3),  p_vec(3)
real(rp) anomalous_moment, charge, mc2, p_z, gamma0
real(rp) e_particle, pc

! Want everything in units of eV

if (init_pauli_vector) call initialize_pauli_vector

pc = ele%value(p0c$) * (1 + coord%vec(6))
call convert_pc_to (pc, coord%species, E_tot = e_particle)

anomalous_moment = anomalous_moment_of(coord%species)
charge = charge_of(coord%species)
mc2 = mass_of(coord%species)
gamma0 = e_particle / mc2
p_z = (ele%value(p0c$)/c_light) * sqrt((1 + coord%vec(6))**2 - coord%vec(2)**2 - coord%vec(4)**2)
p_vec(1:2) = (ele%value(p0c$)/c_light)* [coord%vec(2), coord%vec(4)]
p_vec(3) = p_z

omega = (1 + anomalous_moment*gamma0) * field%B

omega = omega - p_vec * (anomalous_moment*dot_product(p_vec,field%B) / ((gamma0+1)*(mc2**2/c_light**2)))

omega = omega - (1/mc2) * (anomalous_moment + 1/(1+gamma0)) * cross_product(p_vec,field%E)

if (bmad_com%electric_dipole_moment /= 0) then
  omega = omega - (gamma0 * bmad_com%electric_dipole_moment / (2 * c_light)) * &
            (field%E - gamma0 * dot_product(p_vec, field%E) * field%E / ((1 + gamma0) * mc2) + &
             c_light**2 * cross_product(p_vec, field%B))
endif

omega = -(charge/p_z) * omega

end function spin_omega

end module spin_mod
