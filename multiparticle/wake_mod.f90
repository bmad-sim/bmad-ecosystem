module wake_mod

use bmad_struct
use bmad_interface
use multipole_mod, only: ab_multipole_kick
use random_mod

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine randomize_lr_wake_frequencies (ele, set_done)
! 
! Routine to randomize the frequencies of the lr wake HOMs according to:
!   freq = freq_in * (1 + lr_freq_spread) * rr)
! where rr is a Gaussian distributed random number with unit variance.
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele -- ele_struct: Element with wake. If no wake then nothing is done.
!      %value(freq_in$)        -- Frequency.
!
! Output:
!   ele      -- ele_struct: Element with wake frequencies set.
!     %wake%lr(:)%freq -- Set frequency.
!   set_done -- Logical, optional: Set True if there where lr wakes to be set.
!                 False otherwise.
!-

subroutine randomize_lr_wake_frequencies (ele, set_done)

implicit none

type (ele_struct) ele
logical, optional :: set_done
integer n
real(rp) rr

!

if (present(set_done)) set_done = .false.
if (ele%wake%lr_freq_spread == 0 .or. .not. associated(ele%wake)) return

do n = 1, size(ele%wake%lr)
  call ran_gauss (rr)
  ele%wake%lr(n)%freq = ele%wake%lr(n)%freq_in * (1 + ele%wake%lr_freq_spread * rr)
  if (present(set_done)) set_done = .true.
enddo

end subroutine randomize_lr_wake_frequencies

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine zero_lr_wakes_in_lat (lat)
!
! Routine to zero the long range wake amplitudes for the elements that have
! long range wakes in a lattice.
!
! Modules needed:
!   use wake_mod
!
! Input:
!   lat -- Lat_struct: Lattice
!
! Output:
!   lat -- Lat_struct: Lattice
!     %ele(:) -- Lattice elements
!       %wake%lr(:)%b_sin -> Set to zero
!       %wake%lr(:)%b_cos -> Set to zero
!       %wake%lr(:)%a_sin -> Set to zero
!       %wake%lr(:)%a_cos -> Set to zero
!-       

subroutine zero_lr_wakes_in_lat (lat)

implicit none

type (lat_struct) lat
integer i

!

do i = 1, lat%n_ele_max
  if (.not. associated(lat%ele(i)%wake)) cycle
  lat%ele(i)%wake%lr%b_sin = 0; lat%ele(i)%wake%lr%b_cos = 0
  lat%ele(i)%wake%lr%a_sin = 0; lat%ele(i)%wake%lr%a_cos = 0
  lat%ele(i)%wake%lr%t_ref = 0
enddo

end subroutine zero_lr_wakes_in_lat

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine lr_wake_add_to (ele, t0_bunch, orbit)
!
! Subroutine to add to the existing long-range wake the contribution from
! a passing particle.
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele      -- Ele_struct: Element with wakes.
!   t0_bunch -- Real(rp): Time when the bench center was at the start of the lattice.
!   orbit    -- Coord_struct: Particle coords.
!
! Output:
!   ele      -- Ele_struct: Element with wakes.
!     %wake%lr(:)%b_sin -- Non-skew sin-like wake components.
!     %wake%lr(:)%b_cos -- Non-skew cos-like wake components.
!     %wake%lr(:)%a_sin -- Skew sin-like wake components.
!     %wake%lr(:)%a_cos -- Skew cos-like wake components.
!     %wake%lr(:)%t_ref -- Set to t0_bunch.
!+

subroutine lr_wake_add_to (ele, t0_bunch, orbit)

type (ele_struct), target :: ele
type (coord_struct) orbit
type (wake_lr_struct), pointer :: lr

integer i
real(rp) t0_bunch, dt, omega, f_exp, ff, c, s, kx, ky
real(rp) c_a, s_a, kxx, exp_shift, a_sin, b_sin

! Check if we have to do any calculations

if (.not. bmad_com%lr_wakes_on) return  
if (.not. associated(ele%wake)) return
  
! Loop over all modes
! Note: The spatial variation of the normal and skew
! components is the same as the spatial variation of a multipole kick.

! To prevent floating point overflow, the %a and %b factors are shifted 
! to be with respect to lr%t_ref which is the bunch reference time.

do i = 1, size(ele%wake%lr)

  lr => ele%wake%lr(i)

  omega = twopi * lr%freq
  if (lr%freq == 0) omega = twopi * ele%value(rf_frequency$)  ! fundamental mode wake.

  f_exp = omega / (2 * lr%Q)

  if (t0_bunch /= lr%t_ref) then
    dt = t0_bunch - lr%t_ref 
    exp_shift = exp(-dt * f_exp) 
    lr%t_ref = t0_bunch
    lr%b_sin = exp_shift * lr%b_sin
    lr%b_cos = exp_shift * lr%b_cos
    lr%a_sin = exp_shift * lr%a_sin
    lr%a_cos = exp_shift * lr%a_cos
    ! Need to shift a_sin, etc, since particle z is with respect to the bunch center.
    if (lr%freq /= 0) then  ! If not fundamental mode
      c = cos (dt * omega)
      s = sin (dt * omega)
      b_sin = lr%b_sin
      lr%b_sin =  c * b_sin + s * lr%b_cos
      lr%b_cos = -s * b_sin + c * lr%b_cos
      a_sin = lr%a_sin
      lr%a_sin =  c * a_sin + s * lr%a_cos
      lr%a_cos = -s * a_sin + c * lr%a_cos
    endif
  endif

  dt = -orbit%vec(5) * ele%value(p0c$) / (c_light * ele%value(e_tot$))
  if (lr%freq == 0) dt = dt + ele%value(phi0_multipass$) / omega

  ff = abs(orbit%charge) * lr%r_over_q * c_light * exp(dt * f_exp) 

  call ab_multipole_kick (0.0_rp, ff, lr%m, orbit, kx, ky)

  if (lr%polarized) then
    c_a = cos(twopi*lr%angle); s_a = sin(twopi*lr%angle)
    kxx = kx
    kx = kxx * c_a * c_a + ky * s_a * c_a
    ky = kxx * c_a * s_a + ky * s_a * s_a
  endif

  c = cos (-dt * omega)
  s = sin (-dt * omega)

  lr%b_sin = lr%b_sin - kx * c
  lr%b_cos = lr%b_cos + kx * s
  lr%a_sin = lr%a_sin - ky * c
  lr%a_cos = lr%a_cos + ky * s

enddo

end subroutine lr_wake_add_to

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine lr_wake_apply_kick (ele, t0_bunch, orbit)
!
! Subroutine to apply the long-range wake kick to a particle.
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele      -- Ele_struct: Element with wakes
!   t0_bunch -- Real(rp): Time when the bench center was at the start of the lattice.
!   orbit    -- Coord_struct: Starting coords of the particle.
!
! Output:
!   orbit    -- Coord_struct: coords after the kick.
!+

subroutine lr_wake_apply_kick (ele, t0_bunch, orbit)

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (wake_lr_struct), pointer :: lr

integer i
real(rp) t0_bunch, dt, dt_part, omega, f_exp, ff, c, s
real(rp) w_norm, w_skew, kx, ky, k_dum, ff_self, kx_self, ky_self, c_a, s_a

! Check if we have to do any calculations

if (.not. bmad_com%lr_wakes_on) return
if (.not. associated(ele%wake)) return

! Loop over all modes

do i = 1, size(ele%wake%lr)

  lr => ele%wake%lr(i)
  dt_part = -orbit%vec(5) * ele%value(p0c$) / (c_light * ele%value(e_tot$)) 
  dt = t0_bunch + dt_part - lr%t_ref

  omega = twopi * lr%freq
  if (lr%freq == 0) omega = twopi * ele%value(rf_frequency$)  ! fundamental mode wake.

  f_exp = omega / (2 * lr%Q)
  ff = exp(-dt * f_exp) / ele%value(p0c$) 

  if (lr%freq == 0) dt = dt_part + ele%value(phi0_multipass$) / omega

  c = cos (-dt * omega)
  s = sin (-dt * omega)

  ! Self wake component

  ff_self = abs(orbit%charge) * lr%r_over_q * omega / (2 * ele%value(p0c$))

  call ab_multipole_kick (0.0_rp, ff_self, lr%m, orbit, kx_self, ky_self)

  if (lr%polarized) then
    c_a = cos(twopi*lr%angle); s_a = sin(twopi*lr%angle)
    w_norm = -(kx_self * c_a * c_a + ky_self * s_a * c_a)
    w_skew = -(kx_self * c_a * s_a + ky_self * s_a * s_a)
  else
    w_norm = -kx_self
    w_skew = -ky_self
  endif

  ! longitudinal kick

  w_norm = w_norm + (lr%b_sin * ff * (f_exp * s + omega * c) + lr%b_cos * ff * (f_exp * c - omega * s)) / c_light
  w_skew = w_skew + (lr%a_sin * ff * (f_exp * s + omega * c) + lr%a_cos * ff * (f_exp * c - omega * s)) / c_light

  call ab_multipole_kick (0.0_rp, w_norm, lr%m, orbit, kx, k_dum)
  call ab_multipole_kick (0.0_rp, w_skew, lr%m, orbit, k_dum, ky)

  orbit%vec(6) = orbit%vec(6) + kx + ky

  ! transverse kick

  if (lr%m == 0) cycle

  w_norm = lr%b_sin * ff * s + lr%b_cos * ff * c
  w_skew = lr%a_sin * ff * s + lr%a_cos * ff * c

  call ab_multipole_kick (w_skew, w_norm, lr%m-1, orbit, kx, ky)

  orbit%vec(2) = orbit%vec(2) + lr%m * kx
  orbit%vec(4) = orbit%vec(4) + lr%m * ky

enddo

end subroutine lr_wake_apply_kick

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr_long_wake_particle (ele, orbit)
!
! Subroutine to apply the short-range wake kick to a particle and then add 
! to the existing short-range wake the contribution from the particle.
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele     -- Ele_struct: Element with wakes.
!   orbit   -- Coord_struct: Particle coords.
!
! Output:
!   ele     -- Ele_struct: Element with wakes.
!   orbit   -- Coord_struct: coords after the kick.
!+

subroutine sr_long_wake_particle (ele, orbit)

type (ele_struct), target :: ele
type (wake_sr_mode_struct), pointer :: mode
type (coord_struct) orbit

integer i
real(rp) arg, ff, c, s, dz, exp_factor, w_norm

!

dz = orbit%vec(5) - ele%wake%sr_long%z_ref ! Should be negative
ele%wake%sr_long%z_ref = orbit%vec(5)

! Check if we have to do any calculations

do i = 1, size(ele%wake%sr_long%mode)

  mode => ele%wake%sr_long%mode(i)

  ! Kick particle

  exp_factor = exp(dz * mode%damp)

  arg = orbit%vec(5) * mode%k 
  c = cos (arg)
  s = sin (arg)

  ff = abs(orbit%charge) * mode%amp * ele%value(l$) / ele%value(p0c$)

  w_norm = mode%b_sin * exp_factor * s + mode%b_cos * exp_factor * c

  select case (mode%transverse_dependence)
  case (none$, linear_leading$)
    orbit%vec(6) = orbit%vec(6) - w_norm
  case default  ! linear_trailing$
    if (mode%polarization == x_polarization$) then
      orbit%vec(6) = orbit%vec(6) - w_norm * orbit%vec(1)
    else  ! y_axis
      orbit%vec(6) = orbit%vec(6) - w_norm * orbit%vec(3)
    endif
  end select

  ! Self kick

  select case (mode%transverse_dependence)
  case (none$)
    orbit%vec(6) = orbit%vec(6) - ff * sin(mode%phi) / 2
  case default  ! linear_leading or linear_trailing
    if (mode%polarization == x_polarization$) then
      orbit%vec(6) = orbit%vec(6) - orbit%vec(1) * ff * sin(mode%phi) / 2
    else  ! y_axis
      orbit%vec(6) = orbit%vec(6) - orbit%vec(3) * ff * sin(mode%phi) / 2
    endif
  end select

  ! Add to wake

  arg = mode%phi - orbit%vec(5) * mode%k 
  c = cos (arg)
  s = sin (arg)

  ! The monipole wake does not have any skew components.

  select case (mode%transverse_dependence)
  case (none$, linear_trailing$)
    mode%b_sin = mode%b_sin * exp_factor + ff * c
    mode%b_cos = mode%b_cos * exp_factor + ff * s
  case default
    if (mode%polarization == x_polarization$) then
      mode%b_sin = mode%b_sin * exp_factor + orbit%vec(1) * ff * c
      mode%b_cos = mode%b_cos * exp_factor + orbit%vec(1) * ff * s
    else  ! y_axis
      mode%b_sin = mode%b_sin * exp_factor + orbit%vec(3) * ff * c
      mode%b_cos = mode%b_cos * exp_factor + orbit%vec(3) * ff * s
    endif
  end select

enddo

end subroutine sr_long_wake_particle

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr_trans_wake_particle (ele, orbit)
!
! Subroutine to apply the short-range wake kick to a particle and then add 
! to the existing short-range wake the contribution from the particle.
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele     -- Ele_struct: Element with wakes.
!   orbit   -- Coord_struct: Starting particle coords.
!
! Output:
!   ele     -- Ele_struct: Element with wakes.
!   orbit   -- Coord_struct: Ending particle coords.
!+

subroutine sr_trans_wake_particle (ele, orbit)

type (ele_struct), target :: ele
type (wake_sr_mode_struct), pointer :: mode
type (coord_struct) orbit

integer i
real(rp) arg, ff, c, s, dz, exp_factor, w_norm, w_skew

!

dz = orbit%vec(5) - ele%wake%sr_trans%z_ref ! Should be negative
ele%wake%sr_trans%z_ref = orbit%vec(5)

! Add to wake

do i = 1, size(ele%wake%sr_trans%mode)

  mode => ele%wake%sr_trans%mode(i)

  ! Kick particle...

  exp_factor = exp(dz * mode%damp)

  arg = orbit%vec(5) * mode%k 
  c = cos (arg)
  s = sin (arg)

  ! X-axis kick

  if (mode%polarization /= y_polarization$) then
    w_norm = mode%b_sin * exp_factor * s + mode%b_cos * exp_factor * c
    if (mode%transverse_dependence == linear_trailing$) then
      orbit%vec(2) = orbit%vec(2) - w_norm * orbit%vec(1)
    else
      orbit%vec(2) = orbit%vec(2) - w_norm
    endif
  endif

  ! Y-axis kick

  if (mode%polarization /= x_polarization$) then
    w_skew = mode%a_sin * exp_factor * s + mode%a_cos * exp_factor * c
    if (mode%transverse_dependence == linear_trailing$) then
      orbit%vec(4) = orbit%vec(4) - w_skew * orbit%vec(3)
    else
      orbit%vec(4) = orbit%vec(4) - w_skew
    endif
  endif

  ! Add to wake...

  ff = abs(orbit%charge) * mode%amp * ele%value(l$) / ele%value(p0c$)

  arg =  mode%phi - orbit%vec(5) * mode%k 
  c = cos (arg)
  s = sin (arg)

  ! Add to x-axis wake (b_sin, b_cos)

  if (mode%polarization /= y_polarization$) then
    if (mode%transverse_dependence == linear_leading$) then
      mode%b_sin = mode%b_sin * exp_factor + ff * c * orbit%vec(1)
      mode%b_cos = mode%b_cos * exp_factor + ff * s * orbit%vec(1)
    else
      mode%b_sin = mode%b_sin * exp_factor + ff * c
      mode%b_cos = mode%b_cos * exp_factor + ff * s
    endif
  endif

  ! Add to y-axis wake (a_sin, a_cos)

  if (mode%polarization /= x_polarization$) then
    if (mode%transverse_dependence == linear_leading$) then
      mode%a_sin = mode%a_sin * exp_factor + ff * c * orbit%vec(3)
      mode%a_cos = mode%a_cos * exp_factor + ff * s * orbit%vec(3)
    else
      mode%a_sin = mode%a_sin * exp_factor + ff * c
      mode%a_cos = mode%a_cos * exp_factor + ff * s
    endif
  endif

enddo

end subroutine sr_trans_wake_particle

end module
