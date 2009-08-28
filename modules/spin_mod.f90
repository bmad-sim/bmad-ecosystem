module spin_mod

use bmad_struct
use bmad_interface

! right now, jst for electrons (and positrons)
real(rp), parameter :: g_factor = 0.001159657

! so no components are zero
real(rp), parameter :: fudge = 1e-30

! This includes the phase of the spinor
type spin_polar_struct
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

type (spin_map_struct), save, target :: maps(n_key)

private initialize_pauli_vector

real(rp), parameter :: g_factor_of(-2:2) = (/ g_factor_proton, g_factor_electron, 0.0_rp, &
                                            g_factor_electron, g_factor_proton /)

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
!   use bmad
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

real(rp) phi(2), val

character(20) :: r_name = "spinor_to_polar"


  phi(1) = atan (abs(imag(coord%spin(1))) / abs(real(coord%spin(1))))  
  ! get quadrant correct
  if (real(coord%spin(1)) .lt. 0.0 .and. imag(coord%spin(1)) .gt. 0.0) then
    phi(1) = pi - phi(1)
  elseif (real(coord%spin(1)) .lt. 0.0 .and. imag(coord%spin(1)) .lt. 0.0) then
    phi(1) = pi + phi(1)
  elseif (real(coord%spin(1)) .gt. 0.0 .and. imag(coord%spin(1)) .lt. 0.0) then
    phi(1) = -phi(1)
  endif

  phi(2) = atan (abs(imag(coord%spin(2))) / abs(real(coord%spin(2))))  
  ! get quadrant correct
  if (real(coord%spin(2)) .lt. 0.0 .and. imag(coord%spin(2)) .gt. 0.0) then
    phi(2) = pi - phi(2)
  elseif (real(coord%spin(2)) .lt. 0.0 .and. imag(coord%spin(2)) .lt. 0.0) then
    phi(2) = pi + phi(2)
  elseif (real(coord%spin(2)) .gt. 0.0 .and. imag(coord%spin(2)) .lt. 0.0) then
    phi(2) = -phi(2)
  endif

  polar%xi = phi(1)
! polar%phi = modulo(phi(2) - phi(1), 2.0*pi)
  polar%phi = phi(2) - phi(1)

  if (abs(coord%spin(1)) .gt. 1.0) then
    polar%theta = 0.0
  else
    val = abs(coord%spin(1))
    polar%theta  = modulo(2.0d0 * acos (val), pi)
  endif

end subroutine spinor_to_polar

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Suborutine polar_to_vec (polar, vec)
!
! Comverts a spinor in polar coordinates to a spin vector. This will ignore the
! spinor phase.
!
! Modules needed:
!   use bmad
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

  vec(1) = sin(polar%theta) * cos(polar%phi)
  vec(2) = sin(polar%theta) * sin(polar%phi)
  vec(3) = cos(polar%theta)

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
!   use bmad
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

  coord%spin(1) = Exp(i_imaginary * polar%xi) * cos(polar%theta / 2.0d0)
  coord%spin(2) = Exp(i_imaginary * polar%xi) * sin(polar%theta / 2.0d0) * &
                              Exp(i_imaginary * polar%phi)

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
!   use bmad
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

  polar%xi = real_option (0.0d0, phase)

  
  polar%theta = atan(sqrt(vec(1)**2 + vec(2)**2) / abs(vec(3)))
  ! get hemisphere correct
  if (vec(3) .lt. 0.0) polar%theta = pi - polar%theta

  polar%phi   = atan(abs(vec(2)) / abs(vec(1)))

  ! get quadrant right
  if (vec(1) .lt. 0.0 .and. vec(2) .gt. 0.0) then
    polar%phi = pi - polar%phi
  elseif (vec(1) .lt. 0.0 .and. vec(2) .lt. 0.0) then
    polar%phi = pi + polar%phi
  elseif (vec(1) .gt. 0.0 .and. vec(2) .lt. 0.0) then
    polar%phi = -polar%phi
  endif

  ! special case where component is zero
  if (vec(1) .eq. 0.0 .and. vec(2) .gt. 0.0) then
    polar%phi = pi/2.0 
  elseif (vec(1) .eq. 0.0 .and. vec(2) .lt. 0.0) then
    polar%phi = -pi/2.0
  elseif (vec(2) .eq. 0.0 .and. vec(1) .gt. 0.0) then
    polar%phi = 0.0
  elseif (vec(2) .eq. 0.0 .and. vec(1) .lt. 0.0) then
    polar%phi = pi
  endif

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
!   use bmad
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
type (spin_polar_struct) polar

real(rp) vec(3)

complex(rp) complex(3), complex2(2)

integer i

  call initialize_pauli_vector

! complex2 = conjg(coord%spin)
! 
! forall (i=1:3)
!   complex(i) = dot_product(conjg(coord%spin), matmul(pauli(i)%sigma, coord%spin))
! endforall

! vec = real(complex)
  
  call spinor_to_polar (coord, polar)
 
  call polar_to_vec (polar, vec)

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
!   use bmad
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
! Finds the angle between two polar vectors
!
! Modules needed:
!   use bmad
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

real(rp) :: angle
real(rp) :: vec1(3), vec2(3)

  call polar_to_vec (polar1, vec1)
  call polar_to_vec (polar2, vec2)
  
  angle = acos(dot_product(vec1,vec2) / &
               (sqrt(dot_product(vec1, vec1)) * sqrt(dot_product(vec2,vec2))))

end function angle_between_polars

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine quaternion_track (a, start, end)
!
! Transports a spinor through the quaternion a
!
! Modules needed:
!   use bmad
!
! Input:
!   a           -- Real(rp): Euler four-vector (Quaternion)
!   start       -- Coord_struct: Incoming spinor
!
! Output:
!   end      -- Coord_struct
!      %spin    -- complex(rp): Resultant spinor
!-

subroutine quaternion_track (a, start, end)

implicit none

type (coord_struct) start
type (coord_struct) end

real(rp) a(4)

complex(rp) a_quat(2,2) ! The quaternion from the Euler parameters

  call initialize_pauli_vector

  a_quat(1,:) = (/ (1.0d0, 0.0d0), (0.0d0, 0.0d0) /) 
  a_quat(2,:) = (/ (0.0d0, 0.0d0), (1.0d0, 0.0d0) /)

  a_quat = a(4) * a_quat

  a_quat = a_quat - i_imaginary * &
            (a(1) * pauli(1)%sigma + a(2) * pauli(2)%sigma + a(3) * pauli(3)%sigma)

 end%spin =  matmul (a_quat, start%spin)

end subroutine quaternion_track

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! subroutine track1_spin (start, ele, param, end)
!
! Particle spin tracking through a single element.
!
! Uses Nonlinear Spin Transfer Maps from C. Weissbaecker and G. H. Hoffstaetter
!
! For now just does first order transpport. The kappa term is determined fom the
! unitarity condition.
!
! Modules needed:
!   use bmad
!
! Input :
!   start      -- Coord_struct: Starting coords.
!   ele        -- Ele_struct: Element to track through.
!   param      -- lat_param_struct: Beam parameters.
!     %particle     -- Type of particle used
!
! Output:
!   end        -- Coord_struct: Ending coords.
!      %spin      -- complex(rp): Ending spinor
!-
 
subroutine track1_spin (start, ele, param, end)

implicit none

type (coord_struct), intent(in) :: start
type (coord_struct) :: temp
type (coord_struct) :: end
type (ele_struct) :: ele
type (lat_param_struct), intent(in) :: param
type (spin_map_struct), pointer :: map

real(rp) a(4) ! quaternion four-vector
real(rp) omega1, omega_el, xi, gamma0, gammaf, v, x, u
real(rp) alpha, phase, cos_phi, gradient, pc_start, pc_end, k_el, k_el_tilde
real(rp) e_start, e_end, g_ratio, edge_length, beta_start, beta_end
real(rp) g_factor, m_particle

integer key

! Boris tracking does it's own spin tracking
  if (ele%tracking_method .eq. boris$ .or. &
      ele%tracking_method .eq. adaptive_boris$) return

  m_particle = mass_of(param%particle)
  g_factor = g_factor_of(param%particle)

  end%spin = start%spin     ! transfer start to end
  
  key = ele%key
  if (.not. ele%is_on .and. key /= lcavity$) key = drift$

  select case (key)


!-----------------------------------------------
! drift: no change to spin

  case (drift$, rcollimator$, ecollimator$, monitor$, instrument$, pipe$) 

    return
    
!-----------------------------------------------
! kicker, separator

  case (elseparator$, kicker$, hkicker$, vkicker$) 

    return

!-----------------------------------------------
! quadrupole

  case (quadrupole$)

    ! initial:
    omega1 = sqrt(abs(ele%value(k1$)))
    u = omega1*ele%value(l$)

    xi = 1 + g_factor * &
          ((1+start%vec(6)) * ele%value(E_TOT$)) / m_particle
    
    map => maps(quadrupole$)
    
    call allocate_map (map, 2, 2, 0, 0)

    map%gamma1(1)%exp(:) = (/ 0, 0, 1, 0, 0, 0 /)
    map%gamma1(1)%coef   = -(1.0/2.0) * xi * omega1 * sinh(u)
    ! take into account sign of quadrupole (focusing or defocusing)
    map%gamma1(1)%coef   = sign(1.0_rp, ele%value(k1$)) * map%gamma1(1)%coef
    map%gamma1(2)%exp(:) = (/ 0, 0, 0, 1, 0, 0 /)
    map%gamma1(2)%coef   = -xi * (sinh (u / 2.0))**2

    map%gamma2(1)%exp(:) = (/ 1, 0, 0, 0, 0, 0 /)
    map%gamma2(1)%coef   = -(1.0/2.0) * xi * omega1 * sin(u)
    ! take into account sign of quadrupole (focusing or defocusing)
    map%gamma2(1)%coef   = sign(1.0_rp, ele%value(k1$)) * map%gamma2(1)%coef
    map%gamma2(2)%exp(:) = (/ 0, 1, 0, 0, 0, 0 /)
    map%gamma2(2)%coef   = -xi * (sin (u / 2.0))**2

    ! no gamma3 terms
    
!   map%kappa(1)%exp(:)  = (/ 0, 0, 0, 0, 0, 0 /)
!   map%kappa(1)%coef    = 1.0

!-----------------------------------------------
! sbend

  case (sbend$)

    gamma0 = ((1+start%vec(6)) * ele%value(E_TOT$)) / m_particle
    xi = 1 + g_factor * &
          ((1+start%vec(6)) * ele%value(E_TOT$)) / m_particle
    v = ele%value(g$)*ele%value(l$)
    x = g_factor*gamma0*v
    
    map => maps(sbend$)
    
    call allocate_map (map, 0, 4, 1, 0)

    ! No first order gamma1
    
    map%gamma2(1)%exp(:) = (/ 0, 0, 0, 0, 0, 0 /)
    map%gamma2(1)%coef   = -sin(x / 2.0d0)
    map%gamma2(2)%exp(:) = (/ 1, 0, 0, 0, 0, 0 /)
    map%gamma2(2)%coef   = -(1.0d0/2.0d0) * xi * ele%value(g$) * sin(v) * cos(x / 2.0d0)
    map%gamma2(3)%exp(:) = (/ 0, 1, 0, 0, 0, 0 /)
    map%gamma2(3)%coef   = -xi * cos(x / 2.0d0) * (sin(v / 2.0d0))**2
    map%gamma2(4)%exp(:) = (/ 0, 0, 0, 0, 0, 1 /)
    map%gamma2(4)%coef = ((xi * gamma0 * sin(v) - g_factor * (1+gamma0) * (gamma0-1) * v) / &
         (2.0d0 * (1+gamma0))) * cos(x / 2.0d0)

    map%gamma3(1)%exp(:) = (/ 0, 0, 0, 1, 0, 0 /)
    map%gamma3(1)%coef   = (gamma0-1)/gamma0 * sin(x / 2.0d0)

!   map%kappa(1)%exp(:) = (/ 0, 0, 0, 0, 0, 0 /)
!   map%kappa(1)%coef   = cos(x / 2.0d0)
!   map%kappa(2)%exp(:) = (/ 1, 0, 0, 0, 0, 0 /)
!   map%kappa(2)%coef   = -(1.0/2.0) * xi * ele%value(g$) * sin(v) *  sin(x / 2.0d0)
!   map%kappa(3)%exp(:) = (/ 0, 1, 0, 0, 0, 0 /)
!   map%kappa(3)%coef   =  -xi * (sin(v / 2.0d0))**2 * sin( x / 2.0d0)
!   map%kappa(4)%exp(:) = (/ 0, 0, 0, 0, 0, 1 /)
!   map%kappa(4)%coef   = ((xi * gamma0 * sin(v) - g_factor * (1+gamma0) * (gamma0-1) * v) / &
!        (2.0d0 * (1+gamma0))) * sin(x / 2.0d0)


!-----------------------------------------------
! solenoid

  case (solenoid$)
 
    ! This is a simple zeroeth order transfer matrix
    
    ! rotation angle
    alpha = - (1-g_factor)*ele%value(bs_field$)*ele%value(l$) / (ele%value(p0c$)/c_light)
    
    map => maps(solenoid$)
    
    call allocate_map (map, 0, 0, 1, 0)

    map%gamma3(1)%exp(:) = (/ 0, 0, 0, 0, 0, 0 /)
    map%gamma3(1)%coef   = sin(alpha/2.0)

!   map%kappa(1)%exp(:)  = (/ 0, 0, 0, 0, 0, 0 /)
!   map%kappa(1)%coef    = cos(alpha/2.0)
    
!-----------------------------------------------
! LCavity
!
! Simulates the cavity edge field kicks as electrostatic quadrupoles
! since the quaternions for these have already been found.
!
! Uses the fringe field as calulcated by Hartman and Rosenzweig

  case (lcavity$)

    ! For now, just set to one
    g_ratio = 1

    gamma0 = ((1+start%vec(6)) * ele%value(E_TOT$)) / m_particle
  
    if (ele%value(E_TOT_START$) == 0) then
      print *, 'ERROR IN TRACK1_BMAD: E_TOT_START IS 0 FOR A LCAVITY!'
      call err_exit
    endif

    phase = twopi * (ele%value(phi0$) + ele%value(dphi0$) + ele%value(phi0_err$) - &
                        end%vec(5) * ele%value(rf_frequency$) / c_light)
    cos_phi = cos(phase)
    gradient = (ele%value(gradient$) + ele%value(gradient_err$)) * cos_phi 
    if (.not. ele%is_on) gradient = 0

    if (bmad_com%sr_wakes_on) then
      if (bmad_com%grad_loss_sr_wake /= 0) then  
        ! use grad_loss_sr_wake and ignore e_loss
        gradient = gradient - bmad_com%grad_loss_sr_wake
      else
        gradient = gradient - ele%value(e_loss$) * param%n_part * &
                                                    e_charge / ele%value(l$)
      endif
    endif

    if (gradient == 0) then
      return
    endif

    pc_start = ele%value(p0c_start$) * (1+start%vec(6))
    call convert_pc_to (pc_start, param%particle, &
                                      E_tot = e_start, beta = beta_start)
    e_end = e_start + gradient * ele%value(l$)
    gammaf = gamma0 * (e_end / e_start)

    ! entrance kick is a focusing kick
    
    k_el = gradient / (2 * pc_start)
    omega_el = sqrt(k_el)
    
    k_el_tilde = (e_charge * k_el * (1 + g_factor + (g_factor*gamma0))) / &
                   (omega_el * e_mass * c_light**2 * (1 + gamma0))
    ! The edge field length of a cavity is about 1 quarter wavelength
    edge_length = (c_light * beta_start / ele%value(rf_frequency$)) / 4.0
    
    map => maps(lcavity$)
    
    call allocate_map (map, 2, 2, 0, 0)

    map%gamma1(1)%exp(:) = (/ 0, 0, 1, 0, 0, 0 /)
    map%gamma1(1)%coef   = - (k_el_tilde/2.0) * sin (omega_el * edge_length)
    map%gamma1(2)%exp(:) = (/ 0, 0, 0, 1, 0, 0 /)
    map%gamma1(2)%coef   = - (k_el_tilde/omega_el) * (sin (omega_el * edge_length / 2.0))**2

    map%gamma2(1)%exp(:) = (/ 0, 0, 1, 0, 0, 0 /)
    map%gamma2(1)%coef   = - (k_el_tilde/2.0) * sin (omega_el * edge_length)
    map%gamma2(2)%exp(:) = (/ 0, 0, 0, 1, 0, 0 /)
    map%gamma2(2)%coef   = - (k_el_tilde/omega_el) * (sin (omega_el * edge_length / 2.0))**2

    ! exit kick is a defocusing kick (just add to the entrance kick)
    
    call convert_total_energy_to (e_end, param%particle, &
                                             pc = pc_end, beta = beta_end)
    k_el = gradient / (2 * pc_end)
    omega_el = sqrt(k_el)
    k_el_tilde = (e_charge * k_el * (1 + g_factor + (g_factor*gammaf))) / &
                   (omega_el * e_mass * c_light**2 * (1 + gammaf))

 !   map%gamma1(1)%exp(:) = (/ 0, 0, 1, 0, 0, 0 /)
     map%gamma1(1)%coef   = map%gamma1(1)%coef + (k_el_tilde/2.0) * sinh (omega_el * edge_length)
 !   map%gamma1(2)%exp(:) = (/ 0, 0, 0, 1, 0, 0 /)
     map%gamma1(2)%coef   = map%gamma1(2)%coef + &
                                  (k_el_tilde/omega_el) * (sinh (omega_el * edge_length / 2.0))**2
 
 !   map%gamma2(1)%exp(:) = (/ 0, 0, 1, 0, 0, 0 /)
     map%gamma2(1)%coef   = map%gamma2(1)%coef + (k_el_tilde/2.0) * sinh (omega_el * edge_length)
 !   map%gamma2(2)%exp(:) = (/ 0, 0, 0, 1, 0, 0 /)
     map%gamma2(2)%coef   = map%gamma2(2)%coef + &
                                  (k_el_tilde/omega_el) * (sinh (omega_el * edge_length / 2.0))**2

!-----------------------------------------------
! everything else, just use a drift
! This should be fixed!!!!

  case default

    return

  end select

  temp = end

  call offset_particle (ele, param, temp, set$, set_canonical = .false., &
                        set_hvkicks = .false.)

  call compute_quaternion (map%gamma1, a(1))
  call compute_quaternion (map%gamma2, a(2))
  call compute_quaternion (map%gamma3, a(3))
! call compute_quaternion (map%kappa, a(4))
  
  a(4) = sqrt(1.0 - (a(1)**2 + a(2)**2 + a(3)**2))
  
  call quaternion_track (a, start, temp)
  
  call offset_particle (ele, param, temp, unset$, set_canonical = .false., & 
                        set_hvkicks = .false.)
    
  end%spin = temp%spin

contains

!-------------------------------------------------------------------------
subroutine allocate_map (map, n_gamma1, n_gamma2, n_gamma3, n_kappa)

implicit none

type (spin_map_struct) map
integer n_gamma1, n_gamma2, n_gamma3, n_kappa
    
  
  if (n_gamma1 .eq. 0) then
    if (associated (map%gamma1)) deallocate (map%gamma1)
  else
    if (.not. associated (map%gamma1)) then
      allocate(map%gamma1(n_gamma1))
    elseif (size(map%gamma1) .ne. n_gamma1) then
      deallocate(map%gamma1)
      allocate(map%gamma1(n_gamma1))
    endif
  endif

  if (n_gamma2 .eq. 0) then
    if (associated (map%gamma2)) deallocate (map%gamma2)
  else
    if (.not. associated (map%gamma2)) then
      allocate(map%gamma2(n_gamma2))
    elseif (size(map%gamma2) .ne. n_gamma2) then
      deallocate(map%gamma2)
      allocate(map%gamma2(n_gamma2))
    endif
  endif

  if (n_gamma3 .eq. 0) then
    if (associated (map%gamma3)) deallocate (map%gamma3)
  else
    if (.not. associated (map%gamma3)) then
      allocate(map%gamma3(n_gamma3))
    elseif (size(map%gamma3) .ne. n_gamma3) then
      deallocate(map%gamma3)
      allocate(map%gamma3(n_gamma3))
    endif
  endif

  if (n_kappa .eq. 0) then
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
subroutine compute_quaternion (map, a)

implicit none

type (taylor_term_struct), pointer :: map(:)

real(rp) a, a_part

integer i, j

  a = 0.0
  if (.not. associated(map)) return
  do i = 1, size(map)
    a_part = 1.0
    do j = 1, 6
      a_part = a_part * start%vec(j)**map(i)%exp(j)
    enddo
    a_part = map(i)%coef * a_part
    a = a + a_part
  enddo

end subroutine compute_quaternion

end subroutine track1_spin

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function spin_omega_at (field, coord, ele, param), result (omega)
!
! Return the modified T-BMT spin omega vector.
!
! Omega = - Omega_TBMT / v_z
!
! Modules needed:
!   use bmad
!   use em_field_mod
!
! Input :
!   field      -- em_field_struct: E and B fields
!   coord      -- coord_struct: particle momentum
!   ele        -- ele_struct: element evauluated in
!      %value(E_TOT$) -- reaL(rp): needed to find momentum
!   param      -- lat_param_struct: Beam parameters.
!     %particle     -- Type of particle used
!   omega(3)   -- Real(rp): Omega in cartesian coordinates
!   s          -- Real(rp): evaluate at position s in element
!-

function spin_omega_at (field, coord, ele, param, s) result (omega)

implicit none

type (em_field_struct) :: field
type (coord_struct) :: coord
type (ele_struct) :: ele
type (lat_param_struct) :: param

real(rp) omega(3),  p_vec(3)
real(rp) g_factor, charge, m_particle, p_z, gamma0
real(rp) s, e_particle, pc, phase, cos_phi, gradient

  call initialize_pauli_vector

  ! FIX_ME!!!
  ! get e_particle and pc at position in element
  if (ele%key .eq. lcavity$) then
    phase = twopi * (ele%value(phi0$) + ele%value(dphi0$) + ele%value(phi0_err$) - &
                        coord%vec(5) * ele%value(rf_frequency$) / c_light)
    cos_phi = cos(phase)
    gradient = (ele%value(gradient$) + ele%value(gradient_err$)) * cos_phi 
    if (.not. ele%is_on) gradient = 0
    if (bmad_com%sr_wakes_on) then
      if (bmad_com%grad_loss_sr_wake /= 0) then  
        ! use grad_loss_sr_wake and ignore e_loss
        gradient = gradient - bmad_com%grad_loss_sr_wake
      else
        gradient = gradient - ele%value(e_loss$) * param%n_part * &
                                                    e_charge / ele%value(l$)
      endif
    endif
    pc = ele%value(p0c_start$) * (1 + coord%vec(6))
    call convert_pc_to (pc, param%particle, E_tot = e_particle)
    e_particle = e_particle + gradient*s
  else
    pc = ele%value(p0c$) * (1 + coord%vec(6))
    call convert_pc_to (pc, param%particle, E_tot = e_particle)
  endif
  
  ! want everything in units of Ev
  g_factor = g_factor_of (param%particle)
  charge = charge_of(param%particle)
  m_particle = mass_of(param%particle)
  gamma0 = e_particle / m_particle
  p_z = (ele%value(p0c$)/c_light)*&
                     sqrt((1 + coord%vec(6))**2 - coord%vec(2)**2 - coord%vec(4)**2)
  p_vec(1:2) = (ele%value(p0c$)/c_light)*(/ coord%vec(2), coord%vec(4) /)
  p_vec(3) = p_z
  
  omega = (1 + g_factor*gamma0) * field%B

  omega = omega - ( g_factor*dot_product(p_vec,field%B)   /&
                    ((gamma0+1)*(m_particle**2/c_light**2))  )*p_vec
  
  omega = omega - (1/m_particle) * (g_factor + 1/(1+gamma0))*&
                     cross_product(p_vec,field%E)

  omega = (charge/p_z)*omega
                    
end function spin_omega_at

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function cross_product (vec1, vec2)
!
! Returns the cross product of vec1 x vec2
!
! Modules needed:
!   use bmad
!
! Input :
!   vec1(3)      -- Real(rp)
!   vec2(3)      -- Real(rp)
!-

function cross_product (vec1, vec2) result (vec)

implicit none

real(rp) :: vec1(3), vec2(3), vec(3)

  vec(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
  vec(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
  vec(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

end function cross_product  

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function normalized_quaternion (quat)
!
! Returns the normalized quaternion (preserves the spin unitarity)
!
! Modules needed:
!   use bmad
!
! Input :
!   quat(2,2)       -- complex(rp): the quaternion to normalize
!-

function normalized_quaternion (quat) result (quat_norm)

implicit none

complex(rp) quat(2,2), q11, q21, q12, q22, quat_norm(2,2)

real(rp) a(0:4) ! Euler four-vector

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

end module spin_mod
