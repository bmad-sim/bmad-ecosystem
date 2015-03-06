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

! tracking maps are taylor series
type spin_map_struct
  type (taylor_term_struct), pointer :: gamma1(:) => null() ! quaternion four-vector (gamma1)
  type (taylor_term_struct), pointer :: gamma2(:) => null() ! quaternion four-vector (gamma2)
  type (taylor_term_struct), pointer :: gamma3(:) => null() ! quaternion four-vector (gamma3)
  type (taylor_term_struct), pointer :: kappa(:)  => null() ! quaternion four-vector (kappa)
end type

type (pauli_struct) pauli(3)

logical :: init_pauli_vector = .true. ! Does pauli vector needs to be set up?
logical :: do_print = .true.

! taylor maps for elements
! Keeping map allocationg between calls should speed things up
! So, a map for each element is required

type (spin_map_struct), save, target :: maps(n_key$)

private initialize_pauli_vector

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
! Subroutine spinor_to_polar (coord, polar)
!
! Converts a spinor into a spin polar vector of unit length
!
! Modules needed:
!   use spin_mod
!
! Input:
!   coord%spin(2) -- coord_struct: The particle
!
! Output:
!   polar         -- Spin_polar_struct: The resultant Unitary Vector in polar coordinates
!-

subroutine spinor_to_polar (coord, polar)

implicit none

type (coord_struct) :: coord
type (spin_polar_struct) ::  polar

real(rp) temp(2)

character(20) :: r_name = "spinor_to_polar"

!

temp(1) = atan2 (imag(coord%spin(1)), real(coord%spin(1)))
temp(2) = atan2 (imag(coord%spin(2)), real(coord%spin(2)))
polar%xi = temp(1)
polar%phi = temp(2) - temp(1)

temp=abs(coord%spin)
polar%theta = 2 * atan2(temp(2), temp(1))
! no sqrt here! spinor scales with sqrt(r)
polar%polarization = dot_product(temp, temp)


end subroutine spinor_to_polar

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine polar_to_vec (polar, vec)
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

subroutine polar_to_vec (polar, vec)

implicit none

type (spin_polar_struct) polar

real(rp) vec(3)

vec(1) = polar%polarization * sin(polar%theta) * cos(polar%phi)
vec(2) = polar%polarization * sin(polar%theta) * sin(polar%phi)
vec(3) = polar%polarization * cos(polar%theta)

end subroutine polar_to_vec

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine polar_to_spinor (polar, coord)
!
! Converts a spin vector in polar coordinates to a spinor
!
! Modules needed:
!   use spin_mod
!
! Input:
!   polar          -- spin_polar_struct: includes polar phase
!
! Output:
!   coord%spin(2)   -- coord_struct: the particle spin
!-

subroutine polar_to_spinor (polar, coord)

implicit none

type (spin_polar_struct) polar
type (coord_struct) coord

!

coord%spin(1) = sqrt(polar%polarization) * Exp(i_imaginary * polar%xi) * cos(polar%theta / 2.0d0)
coord%spin(2) = sqrt(polar%polarization) * Exp(i_imaginary * (polar%xi+polar%phi)) * sin(polar%theta / 2.0d0)

end subroutine polar_to_spinor

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine vec_to_polar (vec, polar, phase)
!
! Converts a spin vector to a spin polar
!
! Modules needed:
!   use spin_mod
!
! Input:
!   vec(3)   -- real(rp): unitary spin vector
!   phase    -- real(rp)(Optional): Phase of the spinor, if not given then
!                                   set to zero
!
! Output:
!   polar    -- spin_polar_struct:
!-

subroutine vec_to_polar (vec, polar, phase)

implicit none

type (spin_polar_struct) :: polar

real(rp) vec(3)
real(rp), optional :: phase

!

polar%xi = real_option (0.0d0, phase)
polar%theta = atan2 (sqrt(vec(1)**2 + vec(2)**2), vec(3))
polar%phi = atan2(vec(2), vec(1))
polar%polarization = sqrt(dot_product(vec, vec))

end subroutine vec_to_polar

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine spinor_to_vec (coord, vec)
!
! Converts a spinor to a spin cartesian vector
!
! Modules needed:
!   use spin_mod
!
! Input:
!   coord  -- coord_struct: the particle
!
! Output
!   vec(3) -- Real(rp): spin vector in cartesian coordinates
!-

subroutine spinor_to_vec (coord, vec)

implicit none

type (coord_struct) coord
! type (spin_polar_struct) polar

real(rp) vec(3)

!

! call spinor_to_polar (coord, polar)
! call polar_to_vec (polar, vec)

! vec = conjg(coord%spin) * pauli(i)%sigma * coord%spin done explicitly
vec(1) = 2.*( real(coord%spin(1))*real(coord%spin(2))+aimag(coord%spin(1))*aimag(coord%spin(2)) )
vec(2) = 2.*( real(coord%spin(1))*aimag(coord%spin(2))-aimag(coord%spin(1))*real(coord%spin(2)) )
vec(3) = real(coord%spin(1))**2+aimag(coord%spin(1))**2-real(coord%spin(2))**2-aimag(coord%spin(2))**2

end subroutine spinor_to_vec

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine vec_to_spinor (vec, coord, phase)
!
! Converts a spin cartesian vector to a spinor.
!
! Modules needed:
!   use spin_mod
!
! Input:
!   vec(3)   -- Real(rp): spin vector in cartesian coordinates
!   phase    -- real(rp)(Optional): Phase of the spinor, if not given then
!                                   set to zero
!
! Output:
!   spinor   -- Coord_struct: the particle
!-

subroutine vec_to_spinor (vec, coord, phase)

implicit none

type (coord_struct) coord
type (spin_polar_struct) :: polar

real(rp) vec(3)
real(rp), optional :: phase

real(rp) set_phase

!

set_phase = real_option (0.0d0, phase)

call vec_to_polar (vec, polar, set_phase)
call polar_to_spinor (polar, coord)

end subroutine vec_to_spinor

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

call polar_to_vec (polar1, vec1)
call polar_to_vec (polar2, vec2)

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
! Subroutine quaternion_track (a, spin)
!
! Transports a spinor through the quaternion a
!
! Modules needed:
!   use spin_mod
!
! Input:
!   a(4)       -- Real(rp): Euler four-vector (Quaternion)
!   spin(2)    -- complex(rp): Incoming spinor
!
! Output:
!   spin(2)    -- complex(rp): Resultant spinor
!-

subroutine quaternion_track (a, spin)

implicit none

complex(rp), intent(inout) :: spin(2)

real(rp), intent(in) :: a(4)

complex(rp) a_quat(2,2) ! The quaternion from the Euler parameters

!

call initialize_pauli_vector

a_quat(1,:) = [(1.0d0, 0.0d0), (0.0d0, 0.0d0)]
a_quat(2,:) = [(0.0d0, 0.0d0), (1.0d0, 0.0d0)]

a_quat = a(4) * a_quat

a_quat = a_quat - i_imaginary * &
          (a(1) * pauli(1)%sigma + a(2) * pauli(2)%sigma + a(3) * pauli(3)%sigma)

spin = matmul (a_quat, spin)

end subroutine quaternion_track

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine calc_rotation_quaternion (x, y, z, angle, a)
!
! Calculates the quaternion for a rotation of the spin vector by angle about (x,y,z)
! (x,y,z) has to be a unit vector, i.e. x*x+y*y+z*z=1
!
! Modules needed:
!   use spin_mod
!
! Input:
!   x, y, z    -- Real(rp): Unit rotation axis vector.
!   angle      -- Real(rp): Rotation angle
!
! Output:
!   a(4)       -- Real(rp): Resultant quaternion
!-

subroutine calc_rotation_quaternion (x, y, z, angle, a)

real(rp) , intent(in) :: x, y, z, angle
real(rp) , intent(out) :: a(4)

real(rp) half_angle, s

half_angle = angle/2.
s = -sin(half_angle)

a = [x*s, y*s, z*s, cos(half_angle)]

end subroutine calc_rotation_quaternion

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

type (coord_struct) :: start_orb, end_orb
type (ele_struct) :: ele
type (lat_param_struct) :: param

integer method
character(16), parameter :: r_name = 'track1_spin'
logical err

! Use bmad_standard if spin tracking = tracking$ but particle tracking is not using an integration method.

if (start_orb%species == photon$) return

method = ele%spin_tracking_method
if (.not. any(ele%tracking_method == [boris$, runge_kutta$, symp_lie_ptc$]) .and. &
                                            method == tracking$) method = bmad_standard$

select case (method)
case (bmad_standard$)
  call track1_spin_bmad (start_orb, ele, param, end_orb)

case (custom$)
  call track1_spin_custom (start_orb, ele, param, end_orb, err)

case (symp_lie_ptc$)
  call track1_spin_symp_lie_ptc (start_orb, ele, param, end_orb)

case (tracking$)
  ! Spin tracking here is handled by the particle coordinate tracker.

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
! subroutine track1_spin_symp_lie_ptc (start_orb, ele, param, end_orb)
!
! Particle spin tracking through a single element using Symplectic integration with
! Etienne Forset's PTC code.
!
! Modules needed:
!   use spin_mod
!
! Input :
!   start_orb  -- Coord_struct: Starting coords.
!   ele        -- Ele_struct: Element to track through.
!   param      -- lat_param_struct: Beam parameters.
!
! Output:
!   end_orb    -- Coord_struct: Ending coords.
!      %spin(2)   -- complex(rp): Ending spinor
!-

subroutine track1_spin_symp_lie_ptc (start_orb, ele, param, end_orb)

use ptc_spin, rename_dummy => dp, rename2_dummy => twopi
use ptc_interface_mod

implicit none

type (coord_struct) :: start_orb, end_orb
type (ele_struct) :: ele
type (lat_param_struct) :: param
type (spin_map_struct), pointer :: map
type (probe) spin_probe
type (fibre), pointer :: fibre_ele

real(rp) spin_vec(3), re(6), beta0, beta1

integer key

logical isTreatedHere, isKicker

character(16), parameter :: r_name = 'track1_spin_symp_lie_ptc'

! PTC spin tracking

beta0 = ele%value(p0c_start$) / ele%value(e_tot_start$)

call vec_bmad_to_ptc (start_orb%vec, beta0, re)
call spinor_to_vec (start_orb, spin_vec)

spin_probe = re
spin_probe%s(1)%x = real(spin_vec, dp)

call ele_to_fibre (ele, fibre_ele, param, .true.)
call track_probe (spin_probe, DEFAULT+SPIN0, fibre1 = fibre_ele)

spin_vec = spin_probe%s(1)%x
call vec_to_spinor (spin_vec, end_orb)

end subroutine track1_spin_symp_lie_ptc

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! subroutine track1_spin_bmad (start_orb, ele, param, end_orb)
!
! Particle spin tracking through a single element.
!
! Uses "Nonlinear Spin Transfer Maps" from C. Weissbaecker and G. H. Hoffstaetter
! proceedings of 1999 workshop on Polarized Protons at High Energies
!
! For now just does first order transport. The kappa term is determined from the
! unitarity condition.
!
! Modules needed:
!   use spin_mod
!
! Input :
!   start_orb  -- Coord_struct: Starting coords.
!   ele        -- Ele_struct: Element to track through.
!   param      -- lat_param_struct: Beam parameters.
!   end_orb    -- Coord_struct: Ending coords.
!     %vec          -- Ending particle position
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
type (spin_map_struct), pointer :: map

real(rp) a(4) ! quaternion four-vector
real(rp) omega1, xi, gamma0, gammaf, v, x, u
real(rp) alpha, phase, cos_phi, gradient, pc_start, pc_end
real(rp) e_start, e_end, g_ratio, edge_length, beta_start, beta_end
real(rp) anomalous_moment, m_particle, sign_k, abs_a

integer key

logical isTreatedHere, isKicker

character(16), parameter :: r_name = 'track1_spin_bmad'

!

m_particle = mass_of(start_orb%species)
anomalous_moment = anomalous_moment_of(start_orb%species)

end_orb%spin = start_orb%spin     ! transfer start to end

temp_start = start_orb
temp_end   = end_orb

key = ele%key
if (.not. ele%is_on .and. key /= lcavity$) key = drift$

select case (key)
case (quadrupole$, sbend$, solenoid$, lcavity$)
  isTreatedHere = .true.
  isKicker = .false.
case (kicker$, hkicker$, vkicker$) !elseparator$
  isTreatedHere = .false.
  isKicker = .true.
case default
  isTreatedHere = .false.
  isKicker = .false.
end select

! offset particle coordinates at entrance of element
call offset_particle (ele, param, set$, temp_start, .true., .false., .false.)
! offset particle coordinates at exit of element
call offset_particle (ele, param,   set$, temp_end, .true., .false., .false., .true., ele%value(l$))

call offset_spin (ele, param, temp_start, set$, (isTreatedHere .or. isKicker))

temp_middle%spin = temp_start%spin

if(isTreatedHere) then

  ! rough estimate of particle coordinates in the element
  temp_middle%vec = (temp_start%vec + temp_end%vec)/2.

  select case (key)

  !-----------------------------------------------
  ! drift: no change to spin

!   case (drift$, rcollimator$, ecollimator$, monitor$, instrument$, pipe$)
!
!     return

  !-----------------------------------------------
  ! kicker, separator
  ! note: these are taken into account in offset_spin

!     case (elseparator$, kicker$, hkicker$, vkicker$)
!
!     return

  !-----------------------------------------------
  ! sextupole, octupole, multipole
  ! note: these are taken into account in multipole_spin_precession,
  !       which is called in offset_spin

!     case (sextupole$, octupole$, multipole$)
!
!     return

  !-----------------------------------------------
  ! quadrupole

  case (quadrupole$)

    ! initial:
    omega1 = sqrt(abs(ele%value(k1$)))
    u = omega1*ele%value(l$)

    xi = 1 + anomalous_moment * &
          ((1+temp_middle%vec(6)) * ele%value(E_TOT$)) / m_particle

    map => maps(quadrupole$)

    call allocate_map (map, 2, 2, 0, 1)
    ! take into account sign of quadrupole (focusing or defocusing)
    sign_k = sign(1.0_rp, ele%value(k1$))

    map%gamma1(1)%expn(:) = [0, 0, 1, 0, 0, 0]
    map%gamma1(1)%coef   = -sign_k * xi * omega1 * sinh(u) / 2

    map%gamma1(2)%expn(:) = [0, 0, 0, 1, 0, 0]
    map%gamma1(2)%coef   = -xi * (sinh(u / 2.0))**2

    map%gamma2(1)%expn(:) = [1, 0, 0, 0, 0, 0]
    map%gamma2(1)%coef   = -sign_k * xi * omega1 * sin(u) / 2

    map%gamma2(2)%expn(:) = [0, 1, 0, 0, 0, 0]
    map%gamma2(2)%coef   = -xi * (sin(u / 2.0))**2

    ! no gamma3 terms

    map%kappa(1)%expn(:)  = [0, 0, 0, 0, 0, 0]
    map%kappa(1)%coef    = 1.0

  !-----------------------------------------------
  ! sbend
  ! does not take k1, k2 (quadrupole and sextupole terms) into account

  case (sbend$)

    gamma0 = ((1+temp_middle%vec(6)) * ele%value(E_TOT$)) / m_particle
    xi = 1 + anomalous_moment * gamma0
    v = ele%value(g$)*ele%value(l$)
    x = anomalous_moment*gamma0*v

    map => maps(sbend$)

    call allocate_map (map, 0, 4, 1, 4)

    ! No first order gamma1

    map%gamma2(1)%expn(:) = [0, 0, 0, 0, 0, 0]
    map%gamma2(1)%coef   = -sin(x / 2.0d0)
    map%gamma2(2)%expn(:) = [1, 0, 0, 0, 0, 0]
    map%gamma2(2)%coef   = -0.5d0 * xi * ele%value(g$) * sin(v) * cos(x / 2.0d0)
    map%gamma2(3)%expn(:) = [0, 1, 0, 0, 0, 0]
    map%gamma2(3)%coef   = -xi * cos(x / 2.0d0) * (sin(v / 2.0d0))**2
    map%gamma2(4)%expn(:) = [0, 0, 0, 0, 0, 1]
    map%gamma2(4)%coef = ((xi * gamma0 * sin(v) - anomalous_moment * (1+gamma0) * (gamma0-1) * v) / &
        (2.0d0 * (1+gamma0))) * cos(x / 2.0d0)

    map%gamma3(1)%expn(:) = [0, 0, 0, 1, 0, 0]
    map%gamma3(1)%coef   = (gamma0-1)/gamma0 * sin(x / 2.0d0)

    map%kappa(1)%expn(:) = [0, 0, 0, 0, 0, 0]
    map%kappa(1)%coef   = cos(x / 2.0d0)
    map%kappa(2)%expn(:) = [1, 0, 0, 0, 0, 0]
    map%kappa(2)%coef   = -0.5 * xi * ele%value(g$) * sin(v) *  sin(x / 2.0d0)
    map%kappa(3)%expn(:) = [0, 1, 0, 0, 0, 0]
    map%kappa(3)%coef   =  -xi * (sin(v / 2.0d0))**2 * sin( x / 2.0d0)
    map%kappa(4)%expn(:) = [0, 0, 0, 0, 0, 1]
    map%kappa(4)%coef   = ((xi * gamma0 * sin(v) - anomalous_moment * (1+gamma0) * (gamma0-1) * v) / &
         (2.0d0 * (1+gamma0))) * sin(x / 2.0d0)

  !-----------------------------------------------
  ! solenoid

  case (solenoid$)

    ! This is a simple zeroeth order transfer matrix

    ! rotation angle
    alpha = - (1-anomalous_moment)*ele%value(bs_field$)*ele%value(l$) / (ele%value(p0c$)/c_light)

    map => maps(solenoid$)

    call allocate_map (map, 0, 0, 1, 1)

    map%gamma3(1)%expn(:) = [0, 0, 0, 0, 0, 0]
    map%gamma3(1)%coef   = sin(alpha/2.0)

    map%kappa(1)%expn(:)  = [0, 0, 0, 0, 0, 0]
    map%kappa(1)%coef    = cos(alpha/2.0)

  !-----------------------------------------------
  ! LCavity
  !
  ! Simulates the cavity edge field kicks as electrostatic quadrupoles
  ! since the quaternions for these have already been found.
  !
  ! Uses the fringe field as calculated by Hartman and Rosenzweig

  case (lcavity$)

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

      map => maps(lcavity$)

      call allocate_map (map, 2, 2, 0, 0)

      map%gamma1(1)%expn(:) = [0, 0, 1, 0, 0, 0]
      map%gamma1(2)%expn(:) = [0, 0, 0, 1, 0, 0]
      map%gamma2(1)%expn(:) = [0, 0, 1, 0, 0, 0]
      map%gamma2(2)%expn(:) = [0, 0, 0, 1, 0, 0]

      map%gamma1%coef = 0
      map%gamma2%coef = 0
      call lcav_edge_track (pc_start,  gradient, gamma0, anomalous_moment, edge_length, map)
      call lcav_edge_track (pc_end,   -gradient, gammaf, anomalous_moment, edge_length, map)

    endif

  !-----------------------------------------------
  ! everything else, just use a drift
  ! This should be fixed!!!!

!  case default

  end select

  if ( (key /= lcavity$) .or. (gradient /= 0) ) then
    call compute_quaternion (map%gamma1, a(1))
    call compute_quaternion (map%gamma2, a(2))
    call compute_quaternion (map%gamma3, a(3))
    ! call compute_quaternion (map%kappa, a(4))

    ! Need to insert kappa terms for lcavity. For now, use questionable renormalization for lcavity
    if (key == lcavity$) then
      a(4) = sqrt(1.0 - (a(1)**2 + a(2)**2 + a(3)**2))
    else
      call compute_quaternion (map%kappa, a(4))
      abs_a = sqrt(a(1)**2 + a(2)**2 + a(3)**2 + a(4)**2)
      a = a / abs_a
    endif

    call quaternion_track (a, temp_middle%spin)
  endif
endif

temp_end%spin = temp_middle%spin

call offset_spin (ele, param, temp_end, unset$, (isTreatedHere .or. isKicker))

end_orb%spin = temp_end%spin

!-------------------------------------------------------------------------
contains

subroutine allocate_map (map, n_gamma1, n_gamma2, n_gamma3, n_kappa)

implicit none

type (spin_map_struct) map
integer n_gamma1, n_gamma2, n_gamma3, n_kappa

!

if (n_gamma1 == 0) then
  if (associated (map%gamma1)) deallocate (map%gamma1)
else
  if (.not. associated (map%gamma1)) then
    allocate(map%gamma1(n_gamma1))
  elseif (size(map%gamma1) .ne. n_gamma1) then
    deallocate(map%gamma1)
    allocate(map%gamma1(n_gamma1))
  endif
endif

if (n_gamma2 == 0) then
  if (associated (map%gamma2)) deallocate (map%gamma2)
else
  if (.not. associated (map%gamma2)) then
    allocate(map%gamma2(n_gamma2))
  elseif (size(map%gamma2) .ne. n_gamma2) then
    deallocate(map%gamma2)
    allocate(map%gamma2(n_gamma2))
  endif
endif

if (n_gamma3 == 0) then
  if (associated (map%gamma3)) deallocate (map%gamma3)
else
  if (.not. associated (map%gamma3)) then
    allocate(map%gamma3(n_gamma3))
  elseif (size(map%gamma3) .ne. n_gamma3) then
    deallocate(map%gamma3)
    allocate(map%gamma3(n_gamma3))
  endif
endif

if (n_kappa == 0) then
  if (associated (map%kappa)) deallocate (map%kappa)
else
  if (.not. associated (map%kappa)) then
    allocate(map%kappa(n_kappa))
  elseif (size(map%kappa) .ne. n_kappa) then
    deallocate(map%kappa)
    allocate(map%kappa(n_kappa))
  endif
endif

end subroutine allocate_map

!-------------------------------------------------------------------------
! contains

subroutine compute_quaternion (map, a)

implicit none

type (taylor_term_struct), pointer :: map(:)

real(rp) a, a_part

integer i, j

!

a = 0.0
if (.not. associated(map)) return
do i = 1, size(map)
  a_part = 1.0
  do j = 1, 6
    a_part = a_part * temp_middle%vec(j)**map(i)%expn(j)
  enddo
  a_part = map(i)%coef * a_part
  a = a + a_part
enddo

end subroutine compute_quaternion

end subroutine track1_spin_bmad

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine lcav_edge_track (pc, grad, gam, anomalous_moment, edge_length, map)

implicit none

real(rp) pc, grad, gam, anomalous_moment, edge_length, k_el, k_el_tilde, omega_el
type (spin_map_struct) map

! Is this correct? e_mass is in GeV and not eV!

k_el = abs(grad / (2 * pc))
omega_el = sqrt(k_el)
k_el_tilde = (e_charge * k_el * (1 + anomalous_moment + (anomalous_moment*gam))) / &
                (omega_el * e_mass * c_light**2 * (1 + gam))

! Focusing kick

if (grad > 0) then

  map%gamma1(1)%coef = map%gamma1(1)%coef - (k_el_tilde/2.0) * sin (omega_el * edge_length)
  map%gamma1(2)%coef = map%gamma1(2)%coef - (k_el_tilde/omega_el) * (sin (omega_el * edge_length / 2.0))**2

  map%gamma2(1)%coef = map%gamma2(1)%coef - (k_el_tilde/2.0) * sin (omega_el * edge_length)
  map%gamma2(2)%coef = map%gamma2(2)%coef - (k_el_tilde/omega_el) * (sin (omega_el * edge_length / 2.0))**2

! Defocus kick

else

  map%gamma1(1)%coef = map%gamma1(1)%coef + (k_el_tilde/2.0) * sinh (omega_el * edge_length)
  map%gamma1(2)%coef = map%gamma1(2)%coef + (k_el_tilde/omega_el) * (sinh (omega_el * edge_length / 2.0))**2

  map%gamma2(1)%coef = map%gamma2(1)%coef + (k_el_tilde/2.0) * sinh (omega_el * edge_length)
  map%gamma2(2)%coef = map%gamma2(2)%coef + (k_el_tilde/omega_el) * (sinh (omega_el * edge_length / 2.0))**2

endif

end subroutine lcav_edge_track

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function spin_omega_at (field, coord, ele), result (omega)
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

function spin_omega_at (field, coord, ele) result (omega)

implicit none

type (em_field_struct) :: field
type (coord_struct) :: coord
type (ele_struct) :: ele

real(rp) omega(3),  p_vec(3)
real(rp) anomalous_moment, charge, m_particle, p_z, gamma0
real(rp) e_particle, pc

! Want everything in units of eV

call initialize_pauli_vector

pc = ele%value(p0c$) * (1 + coord%vec(6))
call convert_pc_to (pc, coord%species, E_tot = e_particle)

anomalous_moment = anomalous_moment_of (coord%species)
charge = charge_of(coord%species)
m_particle = mass_of(coord%species)
gamma0 = e_particle / m_particle
p_z = (ele%value(p0c$)/c_light) * sqrt((1 + coord%vec(6))**2 - coord%vec(2)**2 - coord%vec(4)**2)
p_vec(1:2) = (ele%value(p0c$)/c_light)* [coord%vec(2), coord%vec(4)]
p_vec(3) = p_z

omega = (1 + anomalous_moment*gamma0) * field%B

omega = omega - ( anomalous_moment*dot_product(p_vec,field%B) / &
                  ((gamma0+1)*(m_particle**2/c_light**2)) )*p_vec

omega = omega - (1/m_particle) * (anomalous_moment + 1/(1+gamma0)) * cross_product(p_vec,field%E)

omega = -(charge/p_z)*omega

end function spin_omega_at

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function normalized_quaternion (quat)
!
! Returns the normalized quaternion (preserves the spin unitarity)
!
! Modules needed:
!   use spin_mod
!
! Output :
!   quat(2,2)       -- complex(rp): the quaternion to normalize
!-

function normalized_quaternion (quat) result (quat_norm)

implicit none

complex(rp) quat(2,2), q11, q21, q12, q22, quat_norm(2,2)

real(rp) a(0:4) ! Euler four-vector

!

q11 = quat(1,1)
q21 = quat(2,1)
q12 = quat(1,2)
q22 = quat(2,2)

a(0) = (0.0, 0.0)
a(1) = (i_imaginary/2) *  (q12 + q21)
a(2) = (1/2) * (q21 - q12)
a(3) = (i_imaginary/2) * (q11 - q22)

a(0) = sqrt(1.0 - (a(1)**2 + a(2)**2 + a(3)**2))

quat_norm(1,1) = a(0) - i_imaginary * a(3)
quat_norm(2,1) = - i_imaginary * a(1) + a(2)
quat_norm(1,2) = - i_imaginary * a(1) - a(2)
quat_norm(2,2) = a(0) + i_imaginary * a(3)

end function normalized_quaternion

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine offset_spin (ele, param, coord, set, set_tilt,
!                               set_multipoles, set_hvkicks)
! Subroutine to effectively offset an element by instead offsetting
! the spin vectors to correspond to the local element coordinates.
!
! set = set$ assumes the particle is at the entrance end of the element.
! set = unset$ assumes the particle is at the exit end of the element.
!
! Options:
!   Using the element tilt in the offset.
!   Using the HV kicks.
!   [Using the multipoles.]
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele       -- Ele_struct: Element
!     %value(x_pitch$)  -- Horizontal roll of element.
!     %value(y_pitch$)  -- Vertical roll of element.
!     %value(tilt$)     -- Tilt of element.
!     %value(roll$)     -- Roll of dipole.
!   coord     -- Coord_struct: Coordinates of the particle.
!     %spin(2)          -- Particle spin
!   param     -- lat_param_struct:
!   set       -- Logical:
!                   T (= set$)   -> Translate from lab coords to the local
!                                     element coords.
!                   F (= unset$) -> Translate back to lab coords.
!   set_tilt       -- Logical, optional: Default is True.
!                   T -> Rotate using ele%value(tilt$) and
!                            ele%value(roll$) for sbends.
!                   F -> Do not rotate
!   set_multipoles -- Logical, optional: Default is True.
!                   T -> 1/2 of the multipole is applied.
!   set_hvkicks    -- Logical, optional: Default is True.
!                   T -> Apply 1/2 any hkick or vkick.
!
! Output:
!     coord -- Coord_struct: Coordinates of particle.
!
! Currently not implemented: elseparators
!-

subroutine offset_spin (ele, param, coord, set, set_tilt, set_multipoles, set_hvkicks)

use bmad_interface

implicit none

type (ele_struct) :: ele
type (lat_param_struct) :: param
type (coord_struct), intent(inout) :: coord

real(rp), save :: old_angle = 0, old_roll = 0
real(rp), save :: del_x_vel = 0, del_y_vel = 0
real(rp) angle, a_gamma_plus, a(4)

logical, intent(in) :: set
logical, optional, intent(in) :: set_tilt, set_multipoles
logical, optional, intent(in) :: set_hvkicks
logical set_multi, set_hv, set_t, set_hv1, set_hv2

!---------------------------------------------------------------

set_multi = logic_option (.true., set_multipoles)
set_hv    = logic_option (.true., set_hvkicks) .and. ele%is_on .and. &
                   (has_kick_attributes(ele%key) .or. has_hkick_attributes(ele%key))
set_t     = logic_option (.true., set_tilt)  .and. has_orientation_attributes(ele)

if (set_hv) then
  select case (ele%key)
  case (elseparator$, kicker$, hkicker$, vkicker$)
    set_hv1 = .false.
    set_hv2 = .true.
  case default
    set_hv1 = .true.
    set_hv2 = .false.
  end select
else
  set_hv1 = .false.
  set_hv2 = .false.
endif

if (set_t .and. ele%key == sbend$) then
  angle = ele%value(l$) * ele%value(g$)
  if (angle /= old_angle .or. ele%value(roll_tot$) /= old_roll) then
    if (ele%value(roll$) == 0) then
      del_x_vel = 0
      del_y_vel = 0
    else if (abs(ele%value(roll$)) < 0.001) then
      del_x_vel = angle * ele%value(roll$)**2 / 4
      del_y_vel = -angle * sin(ele%value(roll_tot$)) / 2
    else
      del_x_vel = angle * (1 - cos(ele%value(roll_tot$))) / 2
      del_y_vel = -angle * sin(ele%value(roll_tot$)) / 2
    endif
    old_angle = angle
    old_roll = ele%value(roll_tot$)
  endif
endif

a_gamma_plus = anomalous_moment_of(coord%species) * ele%value(e_tot$) * (1 + coord%vec(6)) / mass_of(coord%species) + 1

!----------------------------------------------------------------
! Set...

if (set) then

  ! Setting z_offset done already in offset_particle

  ! Set: pitch
  ! contrary to offset_particle no dependence on E_rel
  if (ele%value(x_pitch_tot$) /= 0) then
    call calc_rotation_quaternion (0._rp, -1._rp, 0._rp, ele%value(x_pitch_tot$), a)
    call quaternion_track (a, coord%spin)
  endif
  if (ele%value(y_pitch_tot$) /= 0) then
    call calc_rotation_quaternion (1._rp, 0._rp, 0._rp, ele%value(y_pitch_tot$), a)
    call quaternion_track (a, coord%spin)
  endif

  ! Set: HV kicks for quads, etc. but not hkicker, vkicker, elsep and kicker elements.
  ! HV kicks must come after z_offset but before any tilts are applied.
  ! Note: Since this is applied before tilt_coords, kicks are independent of any tilt.

  if (set_hv1) then
      if (ele%value(hkick$) /= 0) then
        call calc_rotation_quaternion (0._rp, -1._rp, 0._rp, a_gamma_plus * ele%value(hkick$) / 2, a)
        call quaternion_track (a, coord%spin)
      endif
      if (ele%value(vkick$) /= 0) then
        call calc_rotation_quaternion (1._rp, 0._rp, 0._rp, a_gamma_plus * ele%value(vkick$) / 2, a)
        call quaternion_track (a, coord%spin)
      endif
  endif

  ! Set: Multipoles

  if (set_multi) then
    call multipole_spin_precession (ele, param, coord, .true., &
                               .true., (ele%key==multipole$ .or. ele%key==ab_multipole$))
  endif

  ! Set: Tilt
  ! A non-zero roll has a zeroth order effect that must be included

  if (set_t) then 

    if (ele%key == sbend$) then
      if (ele%value(roll_tot$) /= 0) then
        call calc_rotation_quaternion (0._rp, -1._rp, 0._rp, -a_gamma_plus * del_x_vel, a)
        call quaternion_track (a, coord%spin)
        call calc_rotation_quaternion (1._rp, 0._rp, 0._rp, -a_gamma_plus * del_y_vel, a)
        call quaternion_track (a, coord%spin)
      endif
      if (ele%value(ref_tilt_tot$)+ele%value(roll_tot$) /= 0) then
        call calc_rotation_quaternion (0._rp, 0._rp, 1._rp, ele%value(ref_tilt_tot$)+ele%value(roll_tot$), a)
        call quaternion_track (a, coord%spin)
      endif

    else
      if (ele%value(tilt_tot$) /= 0) then
        call calc_rotation_quaternion (0._rp, 0._rp, 1._rp, ele%value(tilt_tot$), a)
        call quaternion_track (a, coord%spin)
      endif
    endif

  endif

  ! Set: HV kicks for kickers and separators only.
  ! Note: Since this is applied after tilt_coords, kicks are dependent on any tilt.

  if (set_hv2) then
    if (ele%key == elseparator$) then
!     NOT IMPLEMENTED YET
!       if (coord%species < 0) then
!       else
!       endif
    elseif (ele%key == hkicker$) then
      if (ele%value(kick$) /= 0) then
        call calc_rotation_quaternion (0._rp, -1._rp, 0._rp, a_gamma_plus * ele%value(kick$) / 2, a)
        call quaternion_track (a, coord%spin)
      endif
    elseif (ele%key == vkicker$) then
      if (ele%value(kick$) /= 0) then
        call calc_rotation_quaternion (1._rp, 0._rp, 0._rp, a_gamma_plus * ele%value(kick$) / 2, a)
        call quaternion_track (a, coord%spin)
      endif
    else ! i.e. elseif (ele%key == kicker$) then
      if (ele%value(hkick$) /= 0) then
        call calc_rotation_quaternion (0._rp, -1._rp, 0._rp, a_gamma_plus * ele%value(hkick$) / 2, a)
        call quaternion_track (a, coord%spin)
      endif
      if (ele%value(vkick$) /= 0) then
        call calc_rotation_quaternion (1._rp, 0._rp, 0._rp, a_gamma_plus * ele%value(vkick$) / 2, a)
        call quaternion_track (a, coord%spin)
      endif
    endif
  endif

!----------------------------------------------------------------
! Unset...

else

  ! Unset: HV kicks for kickers and separators only.

  if (set_hv2) then
    if (ele%key == elseparator$) then
!     NOT IMPLEMENTED YET
!       if (coord%species < 0) then
!       else
!       endif
    elseif (ele%key == hkicker$) then
      if (ele%value(kick$) /= 0) then
        call calc_rotation_quaternion (0._rp, -1._rp, 0._rp, a_gamma_plus * ele%value(kick$) / 2, a)
        call quaternion_track (a, coord%spin)
      endif
    elseif (ele%key == vkicker$) then
      if (ele%value(kick$) /= 0) then
        call calc_rotation_quaternion (1._rp, 0._rp, 0._rp, a_gamma_plus * ele%value(kick$) / 2, a)
        call quaternion_track (a, coord%spin)
      endif
    else ! i.e. elseif (ele%key == kicker$) then
      if (ele%value(vkick$) /= 0) then
        call calc_rotation_quaternion (1._rp, 0._rp, 0._rp, a_gamma_plus * ele%value(vkick$) / 2, a)
        call quaternion_track (a, coord%spin)
      endif
      if (ele%value(hkick$) /= 0) then
        call calc_rotation_quaternion (0._rp, -1._rp, 0._rp, a_gamma_plus * ele%value(hkick$) / 2, a)
        call quaternion_track (a, coord%spin)
      endif
    endif
  endif


  ! Unset: Tilt

  if (set_t) then
    if (ele%key == sbend$) then
      if (ele%value(tilt_tot$)+ele%value(roll_tot$) /= 0) then
        call calc_rotation_quaternion (0._rp, 0._rp, 1._rp, -(ele%value(tilt_tot$)+ele%value(roll_tot$)), a)
        call quaternion_track (a, coord%spin)
      endif
      if (ele%value(roll_tot$) /= 0) then
        call calc_rotation_quaternion (1._rp, 0._rp, 0._rp, -a_gamma_plus * del_y_vel, a)
        call quaternion_track (a, coord%spin)
        call calc_rotation_quaternion (0._rp, -1._rp, 0._rp, -a_gamma_plus * del_x_vel, a)
        call quaternion_track (a, coord%spin)
      endif

    else
      if (ele%value(tilt_tot$) /= 0) then
        call calc_rotation_quaternion (0._rp, 0._rp, 1._rp, -ele%value(tilt_tot$), a)
        call quaternion_track (a, coord%spin)
      endif
    endif

  endif

  ! Unset: Multipoles
  if (set_multi) then
    call multipole_spin_precession (ele, param, coord, .true., &
                               .true., (ele%key==multipole$ .or. ele%key==ab_multipole$))
  endif

  ! UnSet: HV kicks for quads, etc. but not hkicker, vkicker, elsep and kicker elements.
  ! HV kicks must come after z_offset but before any tilts are applied.

  if (set_hv1) then
      if (ele%value(vkick$) /= 0) then
        call calc_rotation_quaternion (1._rp, 0._rp, 0._rp, a_gamma_plus * ele%value(vkick$) / 2, a)
        call quaternion_track (a, coord%spin)
      endif
      if (ele%value(hkick$) /= 0) then
        call calc_rotation_quaternion (0._rp, -1._rp, 0._rp, a_gamma_plus * ele%value(hkick$) / 2, a)
        call quaternion_track (a, coord%spin)
      endif
  endif

  ! Unset: (Offset and) pitch

  if (ele%value(y_pitch_tot$) /= 0) then
    call calc_rotation_quaternion (1._rp, 0._rp, 0._rp, -ele%value(y_pitch_tot$), a)
    call quaternion_track (a, coord%spin)
  endif
  if (ele%value(x_pitch_tot$) /= 0) then
    call calc_rotation_quaternion (0._rp, -1._rp, 0._rp, -ele%value(x_pitch_tot$), a)
    call quaternion_track (a, coord%spin)
  endif

endif

end subroutine offset_spin


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine multipole_spin_precession (ele, param, orbit, do_half_prec, &
!                                       include_sextupole_octupole, ref_orb_offset)
!
! Subroutine to track the spins in a multipole field
!
! Track1_spin uses quaternions which are calculated only up to second order,
! which does not take into account higher-order magnets (sextupoles etc.).
! This subroutine tracks spins through those higher-order magnets assuming simple
! T-BMT precession to get a rough estimate of their effects.
!
! Input:
!   ele              -- Ele_struct: Element
!     %value(x_pitch$)        -- Horizontal roll of element.
!     %value(y_pitch$)        -- Vertical roll of element.
!     %value(tilt$)           -- Tilt of element.
!     %value(roll$)           -- Roll of dipole.
!   param            -- Lat_param_struct
!   orbit            -- coord_struct: Coordinates of the particle.
!   do_half_prec     -- Logical, optional: Default is False.
!                          Apply half multipole effect only (for kick-drift-kick model)
!   include_sextupole_octupole  -- Logical, optional: Default is False.
!                          Include the effects of sextupoles and octupoles
!   ref_orb_offset   -- Logical, optional: Default is False.
!                          Rotate the local coordinate system according to the
!                          the dipole component of the multipole
!
! Output:
!   spin(2)          -- Complex(rp): Resultant spinor
!-

subroutine multipole_spin_precession (ele, param, orbit, do_half_prec, &
                                      include_sextupole_octupole, ref_orb_offset)

use multipole_mod, only: multipole_ele_to_ab

implicit none

type (ele_struct), intent(in) :: ele
type (lat_param_struct) param
type (coord_struct) orbit

complex(rp) kick, pos

real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), kick_angle, Bx, By, knl, a_coord(4), a_field(4)
real(rp) charge_dir
integer n

logical, optional, intent(in) :: do_half_prec, include_sextupole_octupole, ref_orb_offset
logical half_prec, sext_oct, ref_orb, has_nonzero_pole

!

half_prec = logic_option (.false., do_half_prec)
sext_oct  = logic_option (.false., include_sextupole_octupole)
ref_orb   = logic_option (.false., ref_orb_offset)

call multipole_ele_to_ab(ele, .true., has_nonzero_pole, an, bn)
if (.not. has_nonzero_pole) return
charge_dir = relative_tracking_charge(orbit, param) * ele%orientation
an = an * charge_dir
bn = bn * charge_dir

select case (ele%key)
  case (sextupole$)
    knl = ele%value(k2$)*ele%value(l$)
  case (octupole$)
    knl = ele%value(k3$)*ele%value(l$)
  case default
    knl = 0.
end select

if (half_prec) then
  an  = an/2.
  bn  = bn/2.
  knl = knl/2.
endif

if (sext_oct) then
  ! add half effect of element to take sextupoles/octupoles into account (kick-drift-kick model)
  select case (ele%key)
  case (sextupole$)
    bn(2) = bn(2) + knl*cos(3.*ele%value(tilt_tot$))/2.
    an(2) = an(2) - knl*sin(3.*ele%value(tilt_tot$))/2.
  case (octupole$)
    bn(3) = bn(3) + knl*cos(4.*ele%value(tilt_tot$))/6.
    an(3) = an(3) - knl*sin(4.*ele%value(tilt_tot$))/6.
  end select
endif

! calculate kick_angle (for particle) and unit vector (Bx, By) parallel to B-field
! according to bmad manual, chapter "physics", section "Magnetic Fields"
! kick = qL/P_0*(B_y+i*Bx) = \sum_n (b_n+i*a_n)*(x+i*y)^n
kick = bn(0)+i_imaginary*an(0)

if (ref_orb) then
  ! calculate rotation of local coordinate system due to dipole component
  kick_angle = abs(kick)
  if ( kick_angle == 0. ) then
    ref_orb = .false.
  else
    Bx = aimag(kick / kick_angle)
    By = real (kick / kick_angle)
    call calc_rotation_quaternion(Bx, By, 0._rp, -kick_angle, a_coord)
  endif
endif

pos = orbit%vec(1)+i_imaginary*orbit%vec(3)
if (pos /= 0.) then
  kick = kick + (bn(1)+i_imaginary*an(1))*pos
  do n = 2, n_pole_maxx
    pos = pos * (orbit%vec(1)+i_imaginary*orbit%vec(3))
    kick = kick + (bn(n)+i_imaginary*an(n))*pos
  enddo
endif

kick_angle = abs(kick)
if ( kick_angle /= 0. ) then
  kick = kick / kick_angle
  Bx = aimag(kick)
  By = real(kick)
  ! precession_angle = kick_angle*(a*gamma+1)
  kick_angle = kick_angle * (anomalous_moment_of(orbit%species) * ele%value(e_tot$) * &
                        (1 + orbit%vec(6)) / mass_of(orbit%species) + 1)
  call calc_rotation_quaternion(Bx, By, 0._rp, kick_angle, a_field)
  call quaternion_track (a_field, orbit%spin)
endif

if (ref_orb) then
  call quaternion_track (a_coord, orbit%spin)
endif

end subroutine multipole_spin_precession

end module spin_mod
