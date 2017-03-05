module spin_mod

use em_field_mod
use rotation_3d_mod

! Pauli matrices

complex(rp), parameter :: pauli_1(2,2) = reshape([(1,0), (0,0), (0,0), (1,0)], [2,2])
complex(rp), parameter :: pauli_x(2,2) = reshape([(0,0), (1,0), (1,0), (0,0)], [2,2])
complex(rp), parameter :: pauli_y(2,2) = reshape([(0,0), (0,1), (0,-1), (0,0)], [2,2])
complex(rp), parameter :: pauli_z(2,2) = reshape([(1,0), (0,0), (0,0), (-1,0)], [2,2])

type pauli_struct
  complex(rp) sigma(2,2)
end type

type (pauli_struct), parameter :: pauli(0:3) = [pauli_struct(pauli_1), pauli_struct(pauli_x), &
                                                pauli_struct(pauli_y), pauli_struct(pauli_z)]

private trapzd_omega, sbend_omega_func, quad_etc_omega_func

contains

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
! Function polar_to_spinor (polar) result (spinor)
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
!   spinor(2)   -- complex(rp): Spinor
!-

function polar_to_spinor (polar) result (spinor)

implicit none

type (spin_polar_struct) polar
complex(rp) :: spinor(2)

!

spinor(1) = sqrt(polar%polarization) * Exp(i_imaginary * polar%xi) * cos(polar%theta / 2.0d0)
spinor(2) = sqrt(polar%polarization) * Exp(i_imaginary * (polar%xi+polar%phi)) * sin(polar%theta / 2.0d0)

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
!   spinor(2)  -- complex(rp): Spinor
!
! Output
!   vec(3)     -- Real(rp): spin vector in cartesian coordinates
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
! Subroutine rotate_spin_given_field (orbit, sign_z_vel, BL, EL, ds)
!
! Routine to rotate a spin given the integrated magnetic and/or electric field strengths.
!
! Integrated field is the field * length which is independent of particle direction of travel.
! Spin rotation will be done taking into account the particle direction of travel.
!
! Modules needed:
!   use spin_mod
!
! Input:
!   orbit   -- coord_struct: Initial orbit.
!   BL(3)   -- real(rp), optional: Integrated field strength. Assumed zero if not present.
!   EL(3)   -- real(rp), optional: Integrated field strength. Assumed zero if not present.
!
! Output:
!   orbit   -- coord_struct: Orbit with rotated spin
!-

subroutine rotate_spin_given_field (orbit, sign_z_vel, BL, EL)

implicit none

type (coord_struct) orbit
type (em_field_struct) field

real(rp), optional :: BL(3), EL(3)
real(rp)  omega(3)

integer sign_z_vel

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

omega = spin_omega (field, orbit, sign_z_vel)
call rotate_spin (omega, orbit%spin)

end subroutine rotate_spin_given_field

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine rotate_spin_a_step (orbit, field, ele, ds)
!
! Routine to rotate the spin through an integration step.
! Note: It is assumed that the orbit coords are in the element ref frame and not the lab frame.
!
! Modules needed:
!   use spin_mod
!
! Input:
!   orbit   -- coord_struct: Initial orbit.
!   field   -- em_field_struct: EM Field 
!   ele     -- ele_struct, Element being tracked through. 
!   ds      -- real(rp): Longitudinal step
!
! Output:
!   orbit   -- coord_struct: Orbit with rotated spin
!-

subroutine rotate_spin_a_step (orbit, field, ele, ds)

implicit none

type (ele_struct) ele
type (coord_struct) orbit
type (em_field_struct) field

real(rp) ds, omega(3)
integer sign_z_vel


!

sign_z_vel = ele%orientation * orbit%direction

if (ele%key == sbend$) then
  omega = (1 + ele%value(g$) * orbit%vec(1)) * spin_omega (field, orbit, sign_z_vel) + &
                      [0.0_rp, ele%value(g$)*sign_z_vel, 0.0_rp]
else
  omega = spin_omega (field, orbit, orbit%direction * ele%orientation)
endif

call rotate_spin (abs(ds) * omega, orbit%spin)

end subroutine rotate_spin_a_step

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine rotate_spin (rot_vec, spin)
!
! Routine to rotate a spin.
!
! Modules needed:
!   use spin_mod
!
! Input:
!   rot_vec(3)  -- real(rp): Rotation axis. Magnitude of rot_vec is the rotation angle.
!   spin(3)     -- real(rp): Initial coords.
!
! Output:
!   spin(3)     -- real(rp): Final coords.
!-

subroutine rotate_spin (rot_vec, spin)

implicit none

real(rp) :: spin(3), rot_vec(3), axis(3), angle

!

angle = norm2(rot_vec)
if (angle == 0) return

axis = rot_vec / angle

call rotate_vec_given_axis_angle (spin, axis, angle)

end subroutine rotate_spin

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function spin_omega (field, coord, sign_z_vel, phase_space_coords), result (omega)
!
! Return the modified T-BMT spin omega vector:
!   dOmega/d|s|   With phase space coords.
!   dOmega/dt     With time_RK coords.
!
! With phase_space_coords, d|s| is positive in the longitudinal direction of propagation 
! independent of the sign of sign_z_vel.
!
! In a bend, the omega returned should be modified:
!   true_omega = (1 + g*x) * omega + [0, g, 0]
!
! Modules needed:
!   use spin_mod
!   use em_field_mod
!
! Input:
!   field              -- em_field_struct: E and B fields.
!   coord              -- coord_struct: particle momentum.
!   sign_z_vel         -- integer: Direction of the z-velocity. 
!                           With phase space coords: sign_z_vel = coord%direciton * ele%orientation.
!                           With time RK coords:     sign_z_vel = ele%orientation
!   phase_space_coords -- logical, optional: Is coord in standard phase space coordinates or
!                           is it time Runge Kutta coords?
!
! Output:
!   omega(3)   -- real(rp): If phase_space_coords: Omega_TBMT/v_z
!                           If not: Omega_TBMT
!-

function spin_omega (field, coord, sign_z_vel, phase_space_coords) result (omega)

implicit none

type (em_field_struct) :: field
type (coord_struct) :: coord

real(rp) omega(3),  beta_vec(3), pc
real(rp) anomalous_moment, gamma, rel_p, e_particle, mc2, bz2

integer sign_z_vel

logical, optional :: phase_space_coords

! Want everything in units of eV

anomalous_moment = anomalous_moment_of(coord%species)

if (logic_option(.true., phase_space_coords)) then
  rel_p = 1 + coord%vec(6)
  e_particle = coord%p0c * rel_p / coord%beta
  mc2 = mass_of(coord%species)
  gamma = e_particle / mc2
  bz2 = rel_p**2 - coord%vec(2)**2 - coord%vec(4)**2
  if (bz2 < 0) then  ! Particle has unphysical velocity
    omega = 0
    return
  endif
  beta_vec = (coord%beta / rel_p) * [coord%vec(2), coord%vec(4), sign_z_vel * sqrt(bz2)]

else
  e_particle = sqrt(coord%vec(2)**2 + coord%vec(4)**2 + coord%vec(6)**2) / coord%beta
  beta_vec = [coord%vec(2), coord%vec(4), sign_z_vel * coord%vec(6)] / e_particle
  mc2 = mass_of(coord%species)
  gamma = e_particle / mc2
endif

omega = c_light * (1/gamma + anomalous_moment) * field%B
omega = omega - c_light * (gamma * anomalous_moment * dot_product(beta_vec, field%B) / (gamma + 1)) * beta_vec
omega = omega - (anomalous_moment + 1/(1+gamma)) * cross_product(beta_vec, field%E)

if (bmad_com%electric_dipole_moment /= 0) then
  omega = omega + (bmad_com%electric_dipole_moment / 2) * &
            (field%E - (gamma * dot_product(beta_vec, field%E)/ (1 + gamma)) * beta_vec + &
             c_light * cross_product(beta_vec, field%B))
endif

if (logic_option(.true., phase_space_coords)) then
  omega = -(charge_of(coord%species) / (mc2 * abs(beta_vec(3)))) * omega
else
  omega = -(charge_of(coord%species) * c_light / mc2) * omega
endif

end function spin_omega

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
!      %spin(2)   -- complex(rp): Ending spin
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
  case (boris$, runge_kutta$, time_runge_kutta$, symp_lie_ptc$)
    return ! Spin tracking is done at the same time orbital tracking is done
  case (taylor$)
    method = taylor$
  case default
    method = bmad_standard$
  end select
endif

!

select case (method)
case (bmad_standard$)
  call track1_spin_bmad (start_orb, ele, param, end_orb)

case (custom$)
  call track1_spin_custom (start_orb, ele, param, end_orb, err)

! Notice that PTC spin tracking is only done here only when the (orbital) tracking_method is *not* symp_lie_ptc
case (symp_lie_ptc$)
  call track1_symp_lie_ptc (start_orb, ele, param, temp_orb)
  end_orb%spin = temp_orb%spin

case (taylor$)
  call track1_spin_taylor (start_orb, ele, param, end_orb)

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
! subroutine track1_spin_taylor (start_orb, ele, param, end_orb)
!
! Particle spin tracking through a single element with a spin map.
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
!   end_orb     -- Coord_struct:
!     %spin(3)   -- Ending spin
!-

subroutine track1_spin_taylor (start_orb, ele, param, end_orb)

type (coord_struct) :: start_orb, end_orb
type (ele_struct) ele
type (lat_param_struct) param

real(rp) rot(3,3)
character(*), parameter :: r_name = 'track1_spin_taylor'

!

if (.not. associated(ele%spin_taylor(1,1)%term)) then
  call out_io (s_error$, r_name, 'NO SPIN TAYLOR MAP ASSOCIATED WITH ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
  end_orb%spin = start_orb%spin
endif

call track_taylor (start_orb%vec, ele%spin_taylor(1,:), rot(1,:), ele%taylor%ref)
call track_taylor (start_orb%vec, ele%spin_taylor(2,:), rot(2,:), ele%taylor%ref)
call track_taylor (start_orb%vec, ele%spin_taylor(3,:), rot(3,:), ele%taylor%ref)

end_orb%spin = matmul(rot, start_orb%spin)

end subroutine track1_spin_taylor

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
!     %spin(3)       -- Ending spin
!-

subroutine track1_spin_bmad (start_orb, ele, param, end_orb)

use ptc_spin, rename_dummy => dp, rename2_dummy => twopi
use ptc_interface_mod
use fringe_edge_track_mod

implicit none

type (coord_struct) :: start_orb
type (coord_struct) :: temp_start, temp_end, end_orb
type (ele_struct) :: ele
type (lat_param_struct) :: param
type (fringe_edge_info_struct) fringe_info

real(rp) a_quat(0:3) ! quaternion four-vector
real(rp) omega1, xi, gamma0, gammaf, v, x, u, spline_x(0:3), spline_y(0:3)
real(rp) alpha, phase, cos_phi, gradient, pc_start, pc_end, omega(3)
real(rp) e_start, e_end, g_ratio, edge_length, beta_start, beta_end
real(rp) anomalous_moment, m_particle, sign_k, abs_a, coef, s_edge_track
real(rp) vec0(6), ks, kss, length, kl, kl2, c, s, m21, m22, m23, m24

integer key, sign_z_vel

character(*), parameter :: r_name = 'track1_spin_bmad'

!

if (ele%key == patch$) return  ! Spin tracking handled by track_a_patch for patch elements.

m_particle = mass_of(start_orb%species)
anomalous_moment = anomalous_moment_of(start_orb%species)
sign_z_vel = start_orb%direction * ele%orientation

! A slice_slave may or may not span a fringe. calc_next_fringe_edge will figure this out.

temp_start = start_orb
call calc_next_fringe_edge (ele, s_edge_track, fringe_info, temp_start, .true.)

call offset_particle (ele, param, set$, temp_start, .true., .true., .true., set_spin = .true.)
if (fringe_info%particle_at == first_track_edge$) then
  if (fringe_info%ds_edge /= 0) call track_a_drift (temp_start, fringe_info%ds_edge)
  call apply_element_edge_kick (temp_start, fringe_info, ele, param, .true.)
  call calc_next_fringe_edge (ele, s_edge_track, fringe_info, end_orb, .false.)
endif

temp_end   = end_orb
call offset_particle (ele, param, set$, temp_end, .true., .false., .false., .true., ele%value(l$))
if (fringe_info%particle_at == second_track_edge$ .and. fringe_info%ds_edge /= 0) &
                                                  call track_a_drift (temp_end, fringe_info%ds_edge)
temp_end%spin = temp_start%spin

! Notice that the effects of multipoles on spin is handled in offset_particle.

key = ele%key
if (.not. ele%is_on .and. key /= lcavity$) key = null_ele$
if (ele%value(l$) == 0) key = null_ele$

select case (key)

!-----------------------------------------------
! quadrupole, etc.

case (quadrupole$, sextupole$, octupole$, sol_quad$)

  call spline_fit_orbit (temp_start, temp_end, spline_x, spline_y)
  omega = trapzd_omega (ele, quad_etc_omega_func, spline_x, spline_y, temp_start, temp_end, param)
  call rotate_spin(omega, temp_end%spin)

!-----------------------------------------------
! sbend
! does not take k1, k2 (quadrupole and sextupole terms) into account

case (sbend$)

  call spline_fit_orbit (temp_start, temp_end, spline_x, spline_y)
  omega = trapzd_omega (ele, sbend_omega_func, spline_x, spline_y, temp_start, temp_end, param) + &
                                            [0.0_rp, ele%value(g$)*ele%value(l$)*start_orb%direction*ele%orientation, 0.0_rp]
  call rotate_spin(omega, temp_end%spin)

!-----------------------------------------------
! solenoid
! Notice that ks is independent of the ele orientation

case (solenoid$)

  vec0 = temp_start%vec
  ks = ele%value(ks$) / (1+vec0(6)) * relative_tracking_charge(temp_start, param)
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
  call rotate_spin_given_field (temp_end, sign_z_vel, [0.0_rp, 0.0_rp, length * ele%value(bs_field$)])

!-----------------------------------------------
! LCavity, etc

case (lcavity$, rfcavity$, wiggler$)
  temp_start%s = temp_start%s + 1d-10  ! Just to make sure inside field region
  temp_end%s   = temp_end%s   - 1d-10  ! Just to make sure inside field region
  if (ele%key == lcavity$) call reference_energy_correction (ele, temp_start, first_track_edge$)
  call spline_fit_orbit (temp_start, temp_end, spline_x, spline_y)
  omega = trapzd_omega (ele, misc_omega_func, spline_x, spline_y, temp_start, temp_end, param)
  call rotate_spin(omega, temp_end%spin)
end select

!----------

if (fringe_info%particle_at == second_track_edge$) then
  call apply_element_edge_kick (temp_end, fringe_info, ele, param, .true.)
endif
call offset_particle (ele, param, unset$, temp_end, .true., .true., .true., set_spin = .true.)

end_orb%spin = temp_end%spin

end subroutine track1_spin_bmad

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine spline_fit_orbit (start_orb, end_orb, spline_x, spline_y)

implicit none

type (coord_struct) start_orb, end_orb
real(rp) spline_x(0:3), spline_y(0:3)
real(rp) ds, alpha, beta

!

ds = abs(end_orb%s - start_orb%s)

! X

spline_x(0) = start_orb%vec(1)
spline_x(1) = start_orb%vec(2) / (1 + start_orb%vec(6))

alpha = end_orb%vec(1) - spline_x(0) - spline_x(1) * ds
beta = end_orb%vec(2) / (1 + end_orb%vec(6)) - spline_x(1)

spline_x(2) = 3 * alpha / ds**2 - beta / ds
spline_x(3) = beta / ds**2 - 2 * alpha / ds**3

! Y

spline_y(0) = start_orb%vec(3)
spline_y(1) = start_orb%vec(4) / (1 + start_orb%vec(6))

alpha = end_orb%vec(3) - spline_y(0) - spline_y(1) * ds
beta = end_orb%vec(4) / (1 + end_orb%vec(6)) - spline_y(1)

spline_y(2) = 3 * alpha / ds**2 - beta / ds
spline_y(3) = beta / ds**2 - 2 * alpha / ds**3


end subroutine spline_fit_orbit

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

function trapzd_omega (ele, omega_func, spline_x, spline_y, start_orb, end_orb, param) result (omega)

use nr, only: polint

implicit none

interface func
  function omega_func(s, spline_x, spline_y, start_orb, end_orb, ele, param) result (om)
    import
    implicit none
    type (coord_struct) start_orb, end_orb
    type (ele_struct) ele
    type (lat_param_struct) param
    real(rp) s, spline_x(0:3), spline_y(0:3), om(3)
  end function
end interface

type q_array_struct
  real(rp) h
  real(rp) omega(3)
end type

integer, parameter ::  j_max = 10

type (q_array_struct) q_array(j_max), z(0:256)
type (ele_struct) ele
type (coord_struct) start_orb, end_orb, orb
type (lat_param_struct) param

real(rp) s0, s1, del_s, s, spline_x(0:3), spline_y(0:3), omega(3)
real(rp) dint, eps, quat(0:3)
real(rp), parameter :: eps_rel = 1d-5, eps_abs = 1d-8

integer j, k, n, n_pts

! Note: Integrate over start_orb%s to end_orb%s. This will be the whole element except for RF cavities.

s0 = start_orb%s - ele%s_start
s1 = end_orb%s - ele%s_start

q_array(1)%h = 1
z(0)%omega = omega_func(s0, spline_x, spline_y, start_orb, end_orb, ele, param)
z(1)%omega = omega_func(s1, spline_x, spline_y, start_orb, end_orb, ele, param)

del_s = abs(s1 - s0)
q_array(1)%omega = quat_to_omega(quat_mul(omega_to_quat(z(1)%omega * del_s / 2), omega_to_quat(z(0)%omega * del_s / 2)))

do j = 2, j_max
  ! This is trapzd from NR
  n_pts = 2**(j-2)
  del_s = (s1 - s0) / (2 * n_pts)
  quat = omega_to_quat(z(0)%omega * abs(del_s) / 2)

  z(2:2*n_pts:2) = z(1:n_pts) 

  do n = 1, n_pts
    s = s0 + del_s * (2*n - 1)
    z(2*n-1)%omega = omega_func(s, spline_x, spline_y, start_orb, end_orb, ele, param)
    quat = quat_mul(omega_to_quat(z(2*n-1)%omega * abs(del_s)), quat)
    if (n == n_pts) del_s = del_s / 2
    quat = quat_mul(omega_to_quat(z(2*n)%omega * abs(del_s)), quat)
  enddo

  q_array(j)%omega = quat_to_omega(quat)
  q_array(j)%h = q_array(j-1)%h / 4

  eps = eps_abs + eps_rel * sum(abs(q_array(j)%omega))

  if (ele%key == wiggler$ .and. j < 5) cycle  ! Cannot trust until have enough points

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

function sbend_omega_func (s_eval, spline_x, spline_y, start_orb, end_orb, ele, param) result (omega) 

implicit none

type (coord_struct) start_orb, end_orb, orb
type (ele_struct) ele
type (em_field_struct) field
type (lat_param_struct) param

real(rp) s_eval, spline_x(0:3), spline_y(0:3), omega(3), B(3)
real(rp) x, y, ds

!

field%b = [0.0_rp, ele%value(b_field$) + ele%value(b_field_err$), 0.0_rp]
field%e = 0

ds = abs(s_eval - start_orb%s + ele%s_start)
x = spline_x(0) + spline_x(1) * ds + spline_x(2) * ds**2 + spline_x(3) * ds**3
y = spline_y(0) + spline_y(1) * ds + spline_y(2) * ds**2 + spline_y(3) * ds**3

if (ele%value(k1$) /= 0) field%b = field%b + ele%value(b1_gradient$) * [y, x, 0.0_rp]
if (ele%value(k2$) /= 0) field%b = field%b + ele%value(b2_gradient$) * [x*y, (x*x-y*y)/2, 0.0_rp]

! 1 + g*x term comes from the curved coordinates.

orb = start_orb
orb%vec(2) = (1 + orb%vec(6)) * (spline_x(1) + 2 * spline_x(2) * ds + 3 * spline_x(3) * ds**2)
orb%vec(4) = (1 + orb%vec(6)) * (spline_y(1) + 2 * spline_y(2) * ds + 3 * spline_y(3) * ds**2)

omega = (1 + ele%value(g$) * x) * spin_omega (field, orb, start_orb%direction * ele%orientation)

end function sbend_omega_func

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

function quad_etc_omega_func (s_eval, spline_x, spline_y, start_orb, end_orb, ele, param) result (omega) 

implicit none

type (coord_struct) start_orb, end_orb, orb
type (ele_struct) ele
type (em_field_struct) field
type (lat_param_struct) param

real(rp) s_eval, spline_x(0:3), spline_y(0:3), omega(3), B(3)
real(rp) x, y, ds

!

field%e = 0

ds = abs(s_eval - start_orb%s + ele%s_start)
x = spline_x(0) + spline_x(1) * ds + spline_x(2) * ds**2 + spline_x(3) * ds**3
y = spline_y(0) + spline_y(1) * ds + spline_y(2) * ds**2 + spline_y(3) * ds**3

select case (ele%key)
case (quadrupole$);   field%b = ele%value(b1_gradient$) * [y, x, 0.0_rp]
case (sextupole$);    field%b = ele%value(b2_gradient$) * [x*y, (x*x-y*y)/2, 0.0_rp]
case (octupole$);     field%b = ele%value(b3_gradient$) * [(3*x*x*y - y**3)/6, (x**3 - 3*x*y*y)/6, 0.0_rp]
case (sol_quad$);     field%b = ele%value(b1_gradient$) * [y, x, 0.0_rp] + [0.0_rp, 0.0_rp, ele%value(bs_field$)]
end select

!

orb = start_orb
orb%vec(2) = (1 + orb%vec(6)) * (spline_x(1) + spline_x(2) * ds + spline_x(3) * ds**2)
orb%vec(4) = (1 + orb%vec(6)) * (spline_y(1) + spline_y(2) * ds + spline_y(3) * ds**2)

omega = spin_omega (field, orb, start_orb%direction * ele%orientation)

end function quad_etc_omega_func

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

function misc_omega_func (s_eval, spline_x, spline_y, start_orb, end_orb, ele, param) result (omega) 

implicit none

type (coord_struct) start_orb, end_orb, orb
type (ele_struct) ele
type (em_field_struct) field
type (lat_param_struct) param

real(rp) s_eval, spline_x(0:3), spline_y(0:3), omega(3), B(3)
real(rp) x, y, s_tot, ds

!

ds = abs(s_eval - start_orb%s + ele%s_start)
x = spline_x(0) + spline_x(1) * ds + spline_x(2) * ds**2 + spline_x(3) * ds**3
y = spline_y(0) + spline_y(1) * ds + spline_y(2) * ds**2 + spline_y(3) * ds**3

!

s_tot = abs(end_orb%s - start_orb%s)

orb = end_orb
orb%vec(6) = start_orb%vec(6) * (s_tot - ds) / s_tot + end_orb%vec(6) * ds / s_tot

orb%vec(2) = (1 + orb%vec(6)) * (spline_x(1) + spline_x(2) * ds + spline_x(3) * ds**2)
orb%vec(4) = (1 + orb%vec(6)) * (spline_y(1) + spline_y(2) * ds + spline_y(3) * ds**2)
orb%vec(5) = start_orb%vec(5) * (s_tot - ds) / s_tot + end_orb%vec(5) * ds / s_tot
orb%t      = start_orb%t      * (s_tot - ds) / s_tot + end_orb%t      * ds / s_tot
orb%beta   = start_orb%beta   * (s_tot - ds) / s_tot + end_orb%beta   * ds / s_tot

call em_field_calc (ele, param, s_eval, orb, .true., field)

omega = spin_omega (field, orb, start_orb%direction * ele%orientation)

end function misc_omega_func

end module spin_mod
