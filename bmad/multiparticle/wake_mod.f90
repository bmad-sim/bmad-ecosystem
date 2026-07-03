module wake_mod

use multipole_mod
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
! Input:
!   ele -- ele_struct: Element with wake. If no wake then nothing is done.
!      %value(freq_in$)        -- Frequency.
!
! Output:
!   ele      -- ele_struct: Element with wake frequencies set.
!     %wake%lr%mode(:)%freq -- Set frequency.
!   set_done -- Logical, optional: Set True if there where lr wakes to be set.
!                 False otherwise.
!-

subroutine randomize_lr_wake_frequencies (ele, set_done)

implicit none

type (ele_struct), target :: ele
type (wake_lr_struct), pointer :: lr
logical, optional :: set_done
integer n
real(rp) rr

!

lr => ele%wake%lr
if (present(set_done)) set_done = .false.
if (lr%freq_spread == 0 .or. .not. associated(ele%wake)) return

do n = 1, size(lr%mode)
  call ran_gauss (rr)
  if (lr%mode(n)%freq_in < 0) cycle
  lr%mode(n)%freq = lr%mode(n)%freq_in * (1 + lr%freq_spread * rr)
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
! Input:
!   lat -- Lat_struct: Lattice
!
! Output:
!   lat -- Lat_struct: Lattice
!     %ele(:) -- Lattice elements
!       %wake%lr%mode(:)%b_sin -> Set to zero
!       %wake%lr%mode(:)%b_cos -> Set to zero
!       %wake%lr%mode(:)%a_sin -> Set to zero
!       %wake%lr%mode(:)%a_cos -> Set to zero
!-       

subroutine zero_lr_wakes_in_lat (lat)

implicit none

type (lat_struct) lat
integer i

!

do i = 1, lat%n_ele_max
  if (.not. associated(lat%ele(i)%wake)) cycle
  lat%ele(i)%wake%lr%mode%b_sin = 0; lat%ele(i)%wake%lr%mode%b_cos = 0
  lat%ele(i)%wake%lr%mode%a_sin = 0; lat%ele(i)%wake%lr%mode%a_cos = 0
  lat%ele(i)%wake%lr%t_ref = 0
enddo

end subroutine zero_lr_wakes_in_lat

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_lr_wake (bunch, ele)
!
! Subroutine to put in the long-range wakes for particle tracking.
!
! Input:
!   ele         -- ele_struct: Element with wakes.
!   bunch       -- bunch_struct: Bunch to track.
!
! Output:
!   ele         -- Ele_struct: Element with updated wake amplitudes.
!   bunch       -- bunch_struct: Kicked bunch.
!+

subroutine track1_lr_wake (bunch, ele)

implicit none

type (bunch_struct), target :: bunch
type (ele_struct) ele
type (coord_struct), pointer :: particle
type (wake_lr_mode_struct), pointer :: mode

real(rp) t0, dt, dt_phase, kx0, ky0, ff0, w_norm, w_skew
real(rp) omega, f_exp, ff, c_dt, s_dt, kx, ky, vec(6)
real(rp) c_a, s_a, kxx, exp_shift, a_sin, b_sin, t_cut
real(rp) da_sin, da_cos, db_sin, db_cos

integer n_mode, i, j, k, i0, n

! Check to see if we need to do any calc

if (.not. bmad_com%lr_wakes_on) return
if (.not. associated(ele%wake)) return
if (bunch%n_live == 0) return

if (.not. associated(ele%wake)) return
if (ele%wake%lr%amp_scale == 0) return

!

n_mode = size(ele%wake%lr%mode)
if (n_mode == 0) return  

! To prevent floating point overflow, the %a and %b factors are shifted 
! to be with respect to lr%t_ref which is the wake reference time.

i0 = bunch%ix_z(1)
t0 = bunch%particle(i0)%t   ! Time of particle at head of bunch.

do i = 1, size(ele%wake%lr%mode)
  mode => ele%wake%lr%mode(i)

  omega = twopi * mode%freq
  f_exp = mode%damp
  dt = ele%wake%lr%time_scale * (t0 - ele%wake%lr%t_ref)
  exp_shift = exp(-dt * f_exp)

  mode%b_sin = exp_shift * mode%b_sin
  mode%b_cos = exp_shift * mode%b_cos
  mode%a_sin = exp_shift * mode%a_sin
  mode%a_cos = exp_shift * mode%a_cos

  ! Need to shift a_sin, etc, since particle z is with respect to the bunch center.
  if (mode%freq_in >= 0) then  ! If not fundamental mode
    c_dt = cos (dt * omega)
    s_dt = sin (dt * omega)
    b_sin = mode%b_sin
    mode%b_sin =  c_dt * b_sin + s_dt * mode%b_cos
    mode%b_cos = -s_dt * b_sin + c_dt * mode%b_cos
    a_sin = mode%a_sin
    mode%a_sin =  c_dt * a_sin + s_dt * mode%a_cos
    mode%a_cos = -s_dt * a_sin + c_dt * mode%a_cos
  endif
enddo

ele%wake%lr%t_ref = t0

! Loop over all modes: kick particles and update wakes.

do i = 1, size(ele%wake%lr%mode)

  mode => ele%wake%lr%mode(i)

  omega = twopi * mode%freq
  f_exp = mode%damp

  if (mode%polarized) then
    c_a = cos(twopi*mode%angle)
    s_a = sin(twopi*mode%angle)
  endif

  !

  da_sin = 0; da_cos = 0; db_sin = 0; db_cos = 0

  do k = 1, size(bunch%particle)
    particle => bunch%particle(bunch%ix_z(k))
    if (particle%state /= alive$) cycle
    ff0 = ele%wake%lr%amp_scale * abs(particle%charge) * mode%r_over_q

    dt = ele%wake%lr%time_scale * (particle%t - ele%wake%lr%t_ref)
    dt_phase = dt
    if (mode%freq_in < 0 .and. .not. bmad_com%absolute_time_tracking) then
      dt_phase = dt_phase + ele%value(phi0_multipass$) / omega ! Fundamental mode phase shift
    endif

    ! The spatial variation of the normal and skew
    ! components is the same as the spatial variation of a multipole kick.

    call ab_multipole_kick (0.0_rp, 1.0_rp, mode%m, particle%species, +1, particle, kx0, ky0)

    ! longitudinal self-wake

    if (ele%wake%lr%self_wake_on) then
      ff = ff0 * omega / (2 * ele%value(p0c$))

      kx = ff * kx0
      ky = ff * ky0

      if (mode%polarized) then
        w_norm = -(kx * c_a * c_a + ky * s_a * c_a)
        w_skew = -(kx * c_a * s_a + ky * s_a * s_a)
      else
        w_norm = -kx
        w_skew = -ky
      endif

      particle%vec(6) = particle%vec(6) + (w_norm * kx0 + w_skew * ky0) * cos(twopi * mode%phi)
    endif

    ! Longitudinal non-self-wake kick

    ff = exp(-dt * f_exp) / ele%value(p0c$)
    c_dt = cos (dt_phase * omega + twopi * mode%phi)
    s_dt = sin (dt_phase * omega + twopi * mode%phi)

    w_norm = (mode%b_sin * ff * (-f_exp * s_dt + omega * c_dt) + mode%b_cos * ff * (f_exp * c_dt + omega * s_dt)) / c_light
    w_skew = (mode%a_sin * ff * (-f_exp * s_dt + omega * c_dt) + mode%a_cos * ff * (f_exp * c_dt + omega * s_dt)) / c_light

    particle%vec(6) = particle%vec(6) + w_norm * kx0 + w_skew * ky0

    ! Transverse wake kick (Note: Transverse has no self-wake kick)

    if (mode%m /= 0) then
      w_norm = ff * (-mode%b_sin * s_dt + mode%b_cos * c_dt)
      w_skew = ff * (-mode%a_sin * s_dt + mode%a_cos * c_dt)

      call ab_multipole_kick (w_skew, w_norm, mode%m-1, particle%species, +1, particle, kx, ky)

      particle%vec(2) = particle%vec(2) + mode%m * kx
      particle%vec(4) = particle%vec(4) + mode%m * ky
    endif

    ! Update wake amplitudes

    ff = ff0 * c_light * exp(dt * f_exp) 

    if (mode%polarized) then
      kx = ff * (kx0 * c_a * c_a + ky0 * s_a * c_a)
      ky = ff * (kx0 * c_a * s_a + ky0 * s_a * s_a)
    else
      kx = ff * kx0 
      ky = ff * ky0
    endif

    db_sin = db_sin - kx * cos(dt_phase * omega)
    db_cos = db_cos - kx * sin(dt_phase * omega)
    da_sin = da_sin - ky * cos(dt_phase * omega)
    da_cos = da_cos - ky * sin(dt_phase * omega)

  enddo  ! Particles

  ! Add to wake mode.

  mode%b_sin = mode%b_sin + db_sin
  mode%b_cos = mode%b_cos + db_cos
  mode%a_sin = mode%a_sin + da_sin
  mode%a_cos = mode%a_cos + da_cos

enddo  ! Wake modes

end subroutine track1_lr_wake

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr_longitudinal_wake_particle (ele, orbit)
!
! Routine to apply the short-range wake longitudinal component kick to a particle and then add 
! to the existing longitudinal wake the contribution from the particle.
!
! Input:
!   ele     -- Ele_struct: Element with wakes.
!   orbit   -- Coord_struct: Particle coords.
!
! Output:
!   ele     -- Ele_struct: Element with wakes.
!   orbit   -- Coord_struct: coords after the kick.
!+

subroutine sr_longitudinal_wake_particle (ele, orbit)

type (ele_struct), target :: ele
type (wake_sr_mode_struct), pointer :: mode
type (coord_struct) orbit

integer i
real(rp) arg, f0, ff, c, s, dz, exp_factor, w_norm

!

if (ele%wake%sr%amp_scale == 0) return
dz = ele%wake%sr%z_scale * (orbit%vec(5) - ele%wake%sr%z_ref_long) ! Should be negative
ele%wake%sr%z_ref_long = orbit%vec(5)

f0 = ele%wake%sr%amp_scale * abs(orbit%charge) / ele%value(p0c$)
if (ele%wake%sr%scale_with_length) f0 = f0 * ele%value(l$) 

! Loop over wakes

do i = 1, size(ele%wake%sr%long)

  mode => ele%wake%sr%long(i)
  ff = f0 * mode%amp

  ! Kick particle from existing wake

  exp_factor = exp(dz * mode%damp)

  arg = ele%wake%sr%z_scale * orbit%vec(5) * mode%k
  c = cos (arg)
  s = sin (arg)
  w_norm = mode%b_sin * exp_factor * s + mode%b_cos * exp_factor * c

  select case (mode%position_dependence)
  case (none$, x_leading$, y_leading$)
    orbit%vec(6) = orbit%vec(6) - w_norm
  case (x_trailing$)
    orbit%vec(6) = orbit%vec(6) - w_norm * orbit%vec(1)
  case (y_trailing$)
    orbit%vec(6) = orbit%vec(6) - w_norm * orbit%vec(3)
  end select

  ! Self kick

  select case (mode%position_dependence)
  case (none$)
    orbit%vec(6) = orbit%vec(6) - ff * sin(twopi * mode%phi) / 2
  case (x_leading$, x_trailing$)
    orbit%vec(6) = orbit%vec(6) - orbit%vec(1) * ff * sin(twopi * mode%phi) / 2
  case (y_leading$, y_trailing$)
    orbit%vec(6) = orbit%vec(6) - orbit%vec(3) * ff * sin(twopi * mode%phi) / 2
  end select

  ! Add to wake

  arg = twopi * mode%phi - ele%wake%sr%z_scale * orbit%vec(5) * mode%k
  c = cos (arg)
  s = sin (arg)

  ! The monopole wake does not have any skew components.

  select case (mode%position_dependence)
  case (none$, x_trailing$, y_trailing$)
    mode%b_sin = mode%b_sin * exp_factor + ff * c
    mode%b_cos = mode%b_cos * exp_factor + ff * s
  case (x_leading$)
    mode%b_sin = mode%b_sin * exp_factor + orbit%vec(1) * ff * c
    mode%b_cos = mode%b_cos * exp_factor + orbit%vec(1) * ff * s
  case (y_leading$)
    mode%b_sin = mode%b_sin * exp_factor + orbit%vec(3) * ff * c
    mode%b_cos = mode%b_cos * exp_factor + orbit%vec(3) * ff * s
  end select

enddo

end subroutine sr_longitudinal_wake_particle

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr_transverse_wake_particle (ele, orbit)
!
! Subroutine to apply the short-range wake transverse component of the kick to a particle and then add 
! to the existing transverse wake the contribution from the particle.
!
! Input:
!   ele     -- Ele_struct: Element with wakes.
!   orbit   -- Coord_struct: Starting particle coords.
!
! Output:
!   ele     -- Ele_struct: Element with wakes.
!   orbit   -- Coord_struct: Ending particle coords.
!+

subroutine sr_transverse_wake_particle (ele, orbit)

type (ele_struct), target :: ele
type (wake_sr_mode_struct), pointer :: mode
type (coord_struct) orbit

integer i
real(rp) arg, f0, ff, c, s, dz, exp_factor, w_norm, w_skew

!

if (ele%wake%sr%amp_scale == 0) return
dz = ele%wake%sr%z_scale * (orbit%vec(5) - ele%wake%sr%z_ref_trans) ! Should be negative
ele%wake%sr%z_ref_trans = orbit%vec(5)

f0 = ele%wake%sr%amp_scale * abs(orbit%charge) / ele%value(p0c$)
if (ele%wake%sr%scale_with_length) f0 = f0 * ele%value(l$) 

! Loop over all wakes

do i = 1, size(ele%wake%sr%trans)

  mode => ele%wake%sr%trans(i)
  ff = f0 * mode%amp

  ! Kick particle...

  exp_factor = exp(dz * mode%damp)

  arg = ele%wake%sr%z_scale * orbit%vec(5) * mode%k
  c = cos (arg)
  s = sin (arg)

  ! X-axis kick

  if (mode%polarization /= y_polarization$) then
    w_norm = mode%b_sin * exp_factor * s + mode%b_cos * exp_factor * c
    if (mode%position_dependence == trailing$) then
      orbit%vec(2) = orbit%vec(2) - w_norm * orbit%vec(1)
    else
      orbit%vec(2) = orbit%vec(2) - w_norm
    endif
  endif

  ! Y-axis kick

  if (mode%polarization /= x_polarization$) then
    w_skew = mode%a_sin * exp_factor * s + mode%a_cos * exp_factor * c
    if (mode%position_dependence == trailing$) then
      orbit%vec(4) = orbit%vec(4) - w_skew * orbit%vec(3)
    else
      orbit%vec(4) = orbit%vec(4) - w_skew
    endif
  endif

  ! Add to wake...

  arg = twopi * mode%phi - ele%wake%sr%z_scale * orbit%vec(5) * mode%k
  c = cos (arg)
  s = sin (arg)

  ! Add to x-axis wake (b_sin, b_cos)

  if (mode%polarization /= y_polarization$) then
    if (mode%position_dependence == leading$) then
      mode%b_sin = mode%b_sin * exp_factor + ff * c * orbit%vec(1)
      mode%b_cos = mode%b_cos * exp_factor + ff * s * orbit%vec(1)
    else
      mode%b_sin = mode%b_sin * exp_factor + ff * c
      mode%b_cos = mode%b_cos * exp_factor + ff * s
    endif
  endif

  ! Add to y-axis wake (a_sin, a_cos)

  if (mode%polarization /= x_polarization$) then
    if (mode%position_dependence == leading$) then
      mode%a_sin = mode%a_sin * exp_factor + ff * c * orbit%vec(3)
      mode%a_cos = mode%a_cos * exp_factor + ff * s * orbit%vec(3)
    else
      mode%a_sin = mode%a_sin * exp_factor + ff * c
      mode%a_cos = mode%a_cos * exp_factor + ff * s
    endif
  endif

enddo

end subroutine sr_transverse_wake_particle

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr_z_long_wake (ele, bunch, z_ave)
!
! Subroutine to apply the short-range z-wake kick to a particle.
!
! Input:
!   ele           -- ele_struct: Element with wake.
!   bunch         -- bunch_struct: Bunch before wake applied.
!   z_ave         -- real(rp): Average z-position of all live particles.
!
! Output:
!   orbit         -- coord_struct: Ending particle coords.
!   bunch         -- bunch_struct: Bunch before wake applied.
!+

subroutine sr_z_long_wake (ele, bunch, z_ave)

use spline_mod

type (ele_struct), target :: ele
type (bunch_struct), target :: bunch
type (wake_sr_struct), pointer :: sr
type (wake_sr_z_long_struct), pointer :: srz
type (coord_struct), pointer :: p, orbit

real(rp) x, f0, ff, f_add, kick, dz, rz_rel, r1, r2
real(rp) z_ave

integer i, j, ix1, ix2, n2, n_bad, nn
logical ok

character(*), parameter :: r_name = 'sr_z_long_wake'

!

sr => ele%wake%sr
if (sr%amp_scale == 0) return

srz => sr%z_long
if (srz%dz == 0) return

f0 = sr%amp_scale / ele%value(p0c$)
if (sr%scale_with_length) f0 = f0 * ele%value(l$) 

! Compute binned bunch distribution and wake

nn = size(srz%w_out)
n2 = (nn - 1) / 2
srz%w_out = 0
n_bad = 0

do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  if (p%state /= alive$) cycle

  rz_rel = sr%z_scale * (p%vec(5) - z_ave) / srz%dz + n2 + 1 
  ix1 = floor(rz_rel)
  ix2 = ix1 + 1
  if (ix1 < 1 .or. ix2 > nn) then
    n_bad = n_bad + 1
    cycle
  endif

  r1 = (ix2 - rz_rel) * p%charge
  r2 = (rz_rel - ix1) * p%charge

  select case (srz%position_dependence)
  case (none$, x_trailing$, y_trailing$)
    srz%w_out(ix1) = srz%w_out(ix1) + r1
    srz%w_out(ix2) = srz%w_out(ix2) + r2
  case (x_leading$)
    srz%w_out(ix1) = srz%w_out(ix1) + r1 * p%vec(1)
    srz%w_out(ix2) = srz%w_out(ix2) + r2 * p%vec(1)
  case (y_leading$)
    srz%w_out(ix1) = srz%w_out(ix2) + r1 * p%vec(3)
    srz%w_out(ix2) = srz%w_out(ix2) + r2 * p%vec(3)
  end select
enddo

if (n_bad > 0.01 * size(bunch%particle)) then
  call out_io (s_error$, r_name, &
      'The bunch is longer than the sr-z wake can handle for element: ' // ele%name)
  p%state = lost$
  return
endif

call fft_1d(srz%w_out, -1)
srz%w_out = srz%w_out * srz%fw * f0 / nn
call fft_1d(srz%w_out, 1)

! Apply wake
! Notice that p%charge does not appear here.

do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  if (p%state /= alive$) cycle

  rz_rel = sr%z_scale * (p%vec(5) - z_ave) / srz%dz + n2 + 1
  ix1 = MOD(floor(rz_rel) + n2 - 1, nn) + 1
  ix2 = MOD(ix1, nn) + 1

  r2 = MOD(rz_rel, 1.0_rp)
  r1 = 1 - r2

  select case (srz%position_dependence)
  case (none$, x_leading$, y_leading$)
    p%vec(6) = p%vec(6) - (r1 * srz%w_out(ix1) + r2 * srz%w_out(ix2))
  case (x_trailing$)
    p%vec(6) = p%vec(6) - (r1 * srz%w_out(ix1) + r2 * srz%w_out(ix2)) * p%vec(1)
  case (y_trailing$)
    p%vec(6) = p%vec(6) - (r1 * srz%w_out(ix1) + r2 * srz%w_out(ix2)) * p%vec(3)
  end select
enddo

end subroutine sr_z_long_wake

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr_z_taylor_wake (ele, bunch, z_ave)
!
! Subroutine to apply the Taylor-expanded 3D short-range wake (Zagorodnov style)
! to a bunch. The longitudinal point wake is expanded to second order in the
! transverse coordinates of the source and witness particles (13 tabulated
! terms h_ab, see the wake_sr_z_taylor_struct documentation). Longitudinal and
! transverse (Panofsky-Wenzel) kicks are applied.
!
! The algorithm bins the bunch charge, and the charge weighted by the source
! x, y, x*y, and x^2-y^2 moments, onto the wake grid ("generalized currents"),
! convolves with the appropriate term kernels via FFT, and interpolates the
! resulting wake potentials back to the particles.
!
! Input:
!   ele     -- ele_struct: Element with wake.
!   bunch   -- bunch_struct: Bunch of particles.
!   z_ave   -- real(rp): Average z-position of the bunch particles.
!
! Output:
!   bunch   -- bunch_struct: Bunch with wake kicks applied.
!-

subroutine sr_z_taylor_wake (ele, bunch, z_ave)

type (ele_struct), target :: ele
type (bunch_struct), target :: bunch
type (wake_sr_struct), pointer :: sr
type (wake_sr_z_taylor_struct), pointer :: srzt
type (coord_struct), pointer :: p

real(rp) z_ave, f0, rz_rel, r1, r2, x, y
real(rp) wt_z, wt_zy, wt_y, wt_zx, wt_x, wt_zxy, wt_xy, wt_zq, wt_q

integer i, j, k, ix1, ix2, n2, n_bad, nn
logical need_cur(5), need_curd(5), use_pot(9), present_term(13), has_deriv(13)

character(*), parameter :: r_name = 'sr_z_taylor_wake'

! Scratch column indices: generalized currents (fbunch, columns 6-10 hold the
! derivatives of columns 1-5) and wake potentials (w_out). Potentials with "z"
! in the name kick vec(6) (times the witness monomial); the others are
! Panofsky-Wenzel integrated potentials kicking vec(2)/vec(4).

integer, parameter :: ic00 = 1, ic10 = 2, ic01 = 3, ic11 = 4, ic20 = 5
integer, parameter :: iwz = 1, iwzy = 2, iwy = 3, iwzx = 4, iwx = 5, &
                      iwzxy = 6, iwxy = 7, iwzq = 8, iwq = 9

! Generalized current used by each term (indexed like sr_z_taylor_w00$, etc.)

integer, parameter :: term_cur(13) = [ic00, ic10, ic01, ic00, ic00, ic20, ic11, &
                                      ic10, ic10, ic01, ic01, ic00, ic00]

!

sr => ele%wake%sr
if (sr%amp_scale == 0) return

srzt => sr%z_taylor
if (srzt%dz == 0) return

do k = 1, 13
  present_term(k) = allocated(srzt%term(k)%fw) .or. allocated(srzt%term(k)%fw1) .or. &
          srzt%term(k)%r /= 0 .or. srzt%term(k)%l /= 0 .or. srzt%term(k)%c_inv /= 0
  has_deriv(k) = allocated(srzt%term(k)%fw1) .or. srzt%term(k)%l /= 0
enddo
if (.not. any(present_term)) return

f0 = sr%amp_scale / ele%value(p0c$)
if (sr%scale_with_length) f0 = f0 * ele%value(l$)

! Which generalized currents (and their derivatives) are needed?

do i = 1, 5
  need_cur(i) = any(present_term .and. term_cur == i)
  need_curd(i) = any(has_deriv .and. term_cur == i)
enddo
need_cur(ic00) = .true.

use_pot = .false.
use_pot(iwz)   = present_term(sr_z_taylor_w00$) .or. present_term(sr_z_taylor_w01$) .or. &
                 present_term(sr_z_taylor_w02$) .or. present_term(sr_z_taylor_w11$) .or. &
                 present_term(sr_z_taylor_w12$)
use_pot(iwzy)  = present_term(sr_z_taylor_w04$) .or. present_term(sr_z_taylor_w14$) .or. &
                 present_term(sr_z_taylor_w24$)
use_pot(iwy)   = use_pot(iwzy)
use_pot(iwzx)  = present_term(sr_z_taylor_w03$) .or. present_term(sr_z_taylor_w13$) .or. &
                 present_term(sr_z_taylor_w23$)
use_pot(iwx)   = use_pot(iwzx)
use_pot(iwzxy) = present_term(sr_z_taylor_w34$)
use_pot(iwxy)  = use_pot(iwzxy)
use_pot(iwzq)  = present_term(sr_z_taylor_w33$)
use_pot(iwq)   = use_pot(iwzq)

! Allocate scratch space if needed.

nn = size(srzt%f_step)
n2 = (nn - 1) / 2

if (allocated(srzt%fbunch)) then
  if (size(srzt%fbunch, 1) /= nn .or. size(srzt%fbunch, 2) /= 10) deallocate(srzt%fbunch, srzt%w_out)
endif
if (.not. allocated(srzt%fbunch)) allocate(srzt%fbunch(nn,10), srzt%w_out(nn,9))

! Bin the generalized currents.

srzt%fbunch = 0
n_bad = 0

do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  if (p%state /= alive$) cycle

  rz_rel = sr%z_scale * (p%vec(5) - z_ave) / srzt%dz + n2 + 1
  ix1 = floor(rz_rel)
  ix2 = ix1 + 1
  if (ix1 < 1 .or. ix2 > nn) then
    n_bad = n_bad + 1
    cycle
  endif

  r1 = (ix2 - rz_rel) * p%charge
  r2 = (rz_rel - ix1) * p%charge
  x = p%vec(1)
  y = p%vec(3)

  srzt%fbunch(ix1,ic00) = srzt%fbunch(ix1,ic00) + r1
  srzt%fbunch(ix2,ic00) = srzt%fbunch(ix2,ic00) + r2
  if (need_cur(ic10)) then
    srzt%fbunch(ix1,ic10) = srzt%fbunch(ix1,ic10) + r1 * x
    srzt%fbunch(ix2,ic10) = srzt%fbunch(ix2,ic10) + r2 * x
  endif
  if (need_cur(ic01)) then
    srzt%fbunch(ix1,ic01) = srzt%fbunch(ix1,ic01) + r1 * y
    srzt%fbunch(ix2,ic01) = srzt%fbunch(ix2,ic01) + r2 * y
  endif
  if (need_cur(ic11)) then
    srzt%fbunch(ix1,ic11) = srzt%fbunch(ix1,ic11) + r1 * x * y
    srzt%fbunch(ix2,ic11) = srzt%fbunch(ix2,ic11) + r2 * x * y
  endif
  if (need_cur(ic20)) then
    srzt%fbunch(ix1,ic20) = srzt%fbunch(ix1,ic20) + r1 * (x*x - y*y)
    srzt%fbunch(ix2,ic20) = srzt%fbunch(ix2,ic20) + r2 * (x*x - y*y)
  endif
enddo

if (n_bad > 0.01 * size(bunch%particle)) then
  call out_io (s_error$, r_name, &
      'The bunch is longer than the sr z_taylor wake can handle for element: ' // ele%name)
  bunch%particle%state = lost$
  return
endif

! Derivatives of the binned generalized currents (for the w1 and l terms).
! Central differences with one-sided ends, mirroring the ocelot reference.

do k = 1, 5
  if (.not. need_curd(k)) cycle
  srzt%fbunch(1,k+5)  = real(srzt%fbunch(2,k))  - real(srzt%fbunch(1,k))
  srzt%fbunch(nn,k+5) = real(srzt%fbunch(nn,k)) - real(srzt%fbunch(nn-1,k))
  do i = 2, nn-1
    srzt%fbunch(i,k+5) = 0.5_rp * (real(srzt%fbunch(i+1,k)) - real(srzt%fbunch(i-1,k)))
  enddo
  call fft_1d(srzt%fbunch(:,k+5), -1)
enddo

do k = 1, 5
  if (need_cur(k)) call fft_1d(srzt%fbunch(:,k), -1)
enddo

! Accumulate the wake potentials in frequency space. The kernel coefficients
! follow the second-order Taylor expansion: cross terms h_ab with both indices
! transverse appear twice in the double sum, hence the factors of 2. For the
! w33 term the transverse (Panofsky-Wenzel) coefficient is twice the
! longitudinal one since the witness monomial is x^2 - y^2.

srzt%w_out = 0

call zt_add (sr_z_taylor_w00$, ic00, 1.0_rp, 1.0_rp, iwz, 0)
call zt_add (sr_z_taylor_w01$, ic10, 1.0_rp, 1.0_rp, iwz, 0)
call zt_add (sr_z_taylor_w02$, ic01, 1.0_rp, 1.0_rp, iwz, 0)
call zt_add (sr_z_taylor_w11$, ic20, 1.0_rp, 1.0_rp, iwz, 0)
call zt_add (sr_z_taylor_w12$, ic11, 2.0_rp, 2.0_rp, iwz, 0)

call zt_add (sr_z_taylor_w04$, ic00, 1.0_rp, 1.0_rp, iwzy, iwy)
call zt_add (sr_z_taylor_w14$, ic10, 2.0_rp, 2.0_rp, iwzy, iwy)
call zt_add (sr_z_taylor_w24$, ic01, 2.0_rp, 2.0_rp, iwzy, iwy)

call zt_add (sr_z_taylor_w03$, ic00, 1.0_rp, 1.0_rp, iwzx, iwx)
call zt_add (sr_z_taylor_w13$, ic10, 2.0_rp, 2.0_rp, iwzx, iwx)
call zt_add (sr_z_taylor_w23$, ic01, 2.0_rp, 2.0_rp, iwzx, iwx)

call zt_add (sr_z_taylor_w34$, ic00, 2.0_rp, 2.0_rp, iwzxy, iwxy)
call zt_add (sr_z_taylor_w33$, ic00, 1.0_rp, 2.0_rp, iwzq, iwq)

do j = 1, 9
  if (.not. use_pot(j)) cycle
  srzt%w_out(:,j) = srzt%w_out(:,j) * f0 / nn
  call fft_1d(srzt%w_out(:,j), 1)
enddo

! Apply the kicks. Positive tabulated h_ab corresponds to energy loss for the
! longitudinal potentials; the integrated (Panofsky-Wenzel) potentials are
! built so that they add directly to the transverse momenta.

do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  if (p%state /= alive$) cycle

  rz_rel = sr%z_scale * (p%vec(5) - z_ave) / srzt%dz + n2 + 1
  ix1 = MOD(floor(rz_rel) + n2 - 1, nn) + 1
  ix2 = MOD(ix1, nn) + 1

  r2 = MOD(rz_rel, 1.0_rp)
  r1 = 1 - r2
  x = p%vec(1)
  y = p%vec(3)

  wt_z   = 0;  wt_zy = 0;  wt_y = 0;  wt_zx = 0;  wt_x = 0
  wt_zxy = 0;  wt_xy = 0;  wt_zq = 0; wt_q = 0

  if (use_pot(iwz))   wt_z   = r1 * real(srzt%w_out(ix1,iwz))   + r2 * real(srzt%w_out(ix2,iwz))
  if (use_pot(iwzy))  wt_zy  = r1 * real(srzt%w_out(ix1,iwzy))  + r2 * real(srzt%w_out(ix2,iwzy))
  if (use_pot(iwy))   wt_y   = r1 * real(srzt%w_out(ix1,iwy))   + r2 * real(srzt%w_out(ix2,iwy))
  if (use_pot(iwzx))  wt_zx  = r1 * real(srzt%w_out(ix1,iwzx))  + r2 * real(srzt%w_out(ix2,iwzx))
  if (use_pot(iwx))   wt_x   = r1 * real(srzt%w_out(ix1,iwx))   + r2 * real(srzt%w_out(ix2,iwx))
  if (use_pot(iwzxy)) wt_zxy = r1 * real(srzt%w_out(ix1,iwzxy)) + r2 * real(srzt%w_out(ix2,iwzxy))
  if (use_pot(iwxy))  wt_xy  = r1 * real(srzt%w_out(ix1,iwxy))  + r2 * real(srzt%w_out(ix2,iwxy))
  if (use_pot(iwzq))  wt_zq  = r1 * real(srzt%w_out(ix1,iwzq))  + r2 * real(srzt%w_out(ix2,iwzq))
  if (use_pot(iwq))   wt_q   = r1 * real(srzt%w_out(ix1,iwq))   + r2 * real(srzt%w_out(ix2,iwq))

  p%vec(6) = p%vec(6) - (wt_z + wt_zy * y + wt_zx * x + wt_zxy * x * y + wt_zq * (x*x - y*y))
  p%vec(2) = p%vec(2) + wt_x + wt_xy * y + wt_q * x
  p%vec(4) = p%vec(4) + wt_y + wt_xy * x - wt_q * y
enddo

!--------------------------------------------------------------------------
contains

! Accumulate, in frequency space, all pieces of one term (tabulated w, tabulated
! w1, and lumped r, l, c_inv) into the longitudinal potential column iwl and,
! if iwt > 0, the Panofsky-Wenzel integrated potential column iwt.
!
! Sign convention: The accumulated potentials are positive for energy loss (the
! kick applies vec(6) = vec(6) - potential). Relative to this, the w1 and l
! pieces enter with the OPPOSITE sign of the w, r, and c_inv pieces. This
! reproduces the ocelot reference implementation (which itself differs from
! Eq. (B1) of Zagorodnov, Dohlus & Tomin as printed) and is required for
! compatibility with existing wake table files. See the PR notes.

subroutine zt_add (k, ic, coef_l, coef_t, iwl, iwt)

type (wake_sr_z_taylor_term_struct), pointer :: tm
real(rp) coef_l, coef_t, cdz
integer k, ic, iwl, iwt

!

tm => srzt%term(k)
cdz = c_light / srzt%dz

if (allocated(tm%fw)) then
  srzt%w_out(:,iwl) = srzt%w_out(:,iwl) + coef_l * srzt%fbunch(:,ic) * tm%fw
  if (iwt > 0) srzt%w_out(:,iwt) = srzt%w_out(:,iwt) + coef_t * srzt%fbunch(:,ic) * tm%fw_int
endif

! Note: the derivative columns of fbunch are d(charge)/d(index) on the z grid.
! Ocelot differentiates the current with respect to its tau grid (tau = -z),
! which flips the sign; the '+' here reproduces ocelot (verified numerically
! against ocelot on smooth bunches).

if (allocated(tm%fw1)) then
  srzt%w_out(:,iwl) = srzt%w_out(:,iwl) + coef_l * cdz * srzt%fbunch(:,ic+5) * tm%fw1
  if (iwt > 0) srzt%w_out(:,iwt) = srzt%w_out(:,iwt) + coef_t * cdz * srzt%fbunch(:,ic+5) * tm%fw1_int
endif

if (tm%r /= 0) then
  srzt%w_out(:,iwl) = srzt%w_out(:,iwl) + coef_l * tm%r * cdz * srzt%fbunch(:,ic)
  if (iwt > 0) srzt%w_out(:,iwt) = srzt%w_out(:,iwt) + coef_t * tm%r * c_light * srzt%fbunch(:,ic) * srzt%f_step
endif

if (tm%l /= 0) then
  srzt%w_out(:,iwl) = srzt%w_out(:,iwl) + coef_l * tm%l * cdz**2 * srzt%fbunch(:,ic+5)
  if (iwt > 0) srzt%w_out(:,iwt) = srzt%w_out(:,iwt) + coef_t * tm%l * cdz * c_light * srzt%fbunch(:,ic+5) * srzt%f_step
endif

if (tm%c_inv /= 0) then
  srzt%w_out(:,iwl) = srzt%w_out(:,iwl) + coef_l * tm%c_inv * srzt%fbunch(:,ic) * srzt%f_step
  if (iwt > 0) srzt%w_out(:,iwt) = srzt%w_out(:,iwt) + coef_t * tm%c_inv * srzt%fbunch(:,ic) * srzt%f_step_int
endif

end subroutine zt_add

end subroutine sr_z_taylor_wake

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine order_particles_in_z (bunch)
!
! Routine to order the particles longitudinally in terms of decreasing %vec(5).
! That is from large z (head of bunch) to small z.
! Only live particles are ordered.
!
! Input:
!   bunch     -- Bunch_struct: collection of particles.
!     %particle(j)%vec(5) -- Longitudinal position of j^th particle.
!
! Output:
!   bunch     -- bunch_struct: collection of particles.
!     %ix_z(:)   -- Index for the ordering. 
!                   Order is from large z (head of bunch) to small z.
!                   That is: %bunch%ix_z(1) is the particle at the bunch head.
!                   Only live particles are ordered so if particle with index %bunch%ix_z(i)
!                     is dead, all particles with index %bunch%ix_z(j) with j > i are dead.
!-

Subroutine order_particles_in_z (bunch)

implicit none

type (bunch_struct), target :: bunch
type (coord_struct), pointer :: particle(:)
type (coord_struct) temp
integer ix, k, nm, i0, i1, n_max, kk
real(rp) z1, z2

! Init if needed. 

particle => bunch%particle
n_max = size(particle)
nm = n_max

! If first time through

if (bunch%ix_z(1) < 1) then
  call indexer (bunch%particle%vec(5), bunch%ix_z)
  bunch%ix_z(1:nm) = bunch%ix_z(nm:1:-1)
endif

! Order is from large z (head of bunch) to small z.
! This ordering calc is efficient when the particles are already more-or-less ordered to start with.  

ix = 1
do
  if (ix > nm) exit
  i0 = bunch%ix_z(ix)

  if (particle(i0)%state /= alive$) then
    bunch%ix_z(ix:nm) = [bunch%ix_z(ix+1:nm), i0]
    nm = nm - 1
    cycle
  endif

  if (ix >= nm) exit
  i1 = bunch%ix_z(ix+1)

  if (particle(i1)%state /= alive$) then
    bunch%ix_z(ix+1:nm) = [bunch%ix_z(ix+2:nm), i1]
    nm = nm - 1
    cycle
  endif

  if (particle(i0)%vec(5) < particle(i1)%vec(5)) then
    do k = ix-1, 1, -1
      kk = bunch%ix_z(k)
      if (particle(kk)%vec(5) >= particle(i1)%vec(5)) exit
    enddo
  
    bunch%ix_z(k+1:ix+1) = [bunch%ix_z(ix+1), bunch%ix_z(k+1:ix)]
  endif

  ix = ix + 1
enddo

bunch%n_live = nm

end subroutine order_particles_in_z

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_sr_wake (bunch, ele)
!
! Subroutine to apply the short range wake fields to a bunch. 
!
! Input:
!   bunch         -- bunch_struct: Bunch of particles.
!   ele           -- ele_struct: Element with wakefields.
!
! Output:
!   bunch         -- bunch_struct: Bunch with wakefields applied to the particles.
!-

subroutine track1_sr_wake (bunch, ele)

implicit none

type (bunch_struct), target :: bunch
type (ele_struct) ele
type (coord_struct), pointer :: particle
type (coord_struct), pointer :: p(:)

real(rp) sr02
integer i, j, k, i1, i2, n_sr_long, n_sr_trans, k_start, n_live

logical wake_here
character(16) :: r_name = 'track1_sr_wake'

!-----------------------------------

if (.not. bmad_com%sr_wakes_on) return
if (.not. associated(ele%wake)) return

n_live = bunch%n_live
if (n_live == 0) return    ! No one left alive.
p => bunch%particle

! error check and zero wake sums and order particles in z

i1 = bunch%ix_z(1)
i2 = bunch%ix_z(n_live)

if (ele%wake%sr%z_max > 0 .and. p(i1)%vec(5) - p(i2)%vec(5) > ele%wake%sr%z_max) then
  call out_io (s_error$, r_name, &
      'The bunch is longer than the sr wake can handle for element: ' // ele%name)
  p%state = lost$
  return
endif

! Mode wakes

ele%wake%sr%long%b_sin = 0
ele%wake%sr%long%b_cos = 0
ele%wake%sr%long%a_sin = 0
ele%wake%sr%long%a_cos = 0
ele%wake%sr%z_ref_long = p(i1)%vec(5)

ele%wake%sr%trans%b_sin = 0
ele%wake%sr%trans%b_cos = 0
ele%wake%sr%trans%a_sin = 0
ele%wake%sr%trans%a_cos = 0
ele%wake%sr%z_ref_trans = p(i1)%vec(5)

! Z-wake

call sr_z_long_wake(ele, bunch, p((i1+i2)/2)%vec(5))
call sr_z_taylor_wake(ele, bunch, p((i1+i2)/2)%vec(5))

! Loop over all particles in the bunch and apply the mode wakes

do j = 1, n_live
  particle => p(bunch%ix_z(j))  ! Particle to kick
  call sr_longitudinal_wake_particle (ele, particle)
  call sr_transverse_wake_particle (ele, particle)
enddo

end subroutine track1_sr_wake

end module
