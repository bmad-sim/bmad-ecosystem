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
!      %value(lr_freq_spread$) -- Fractional RMS frequency spread.
!
! Output:
!   ele      -- ele_struct: Element with wake frequencies set.
!     %rf_wake%lr(:)%freq -- Set frequency.
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
if (ele%value(lr_freq_spread$) == 0 .or. .not. associated(ele%rf_wake)) return

do n = 1, size(ele%rf_wake%lr)
  call ran_gauss (rr)
  ele%rf_wake%lr(n)%freq = ele%rf_wake%lr(n)%freq_in * (1 + ele%value(lr_freq_spread$) * rr)
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
!       %rf_wake%lr(:)%b_sin -> Set to zero
!       %rf_wake%lr(:)%b_cos -> Set to zero
!       %rf_wake%lr(:)%a_sin -> Set to zero
!       %rf_wake%lr(:)%a_cos -> Set to zero
!-       

subroutine zero_lr_wakes_in_lat (lat)

implicit none

type (lat_struct) lat
integer i

!

do i = 1, lat%n_ele_max
  if (.not. associated(lat%ele(i)%rf_wake)) cycle
  lat%ele(i)%rf_wake%lr%b_sin = 0; lat%ele(i)%rf_wake%lr%b_cos = 0
  lat%ele(i)%rf_wake%lr%a_sin = 0; lat%ele(i)%rf_wake%lr%a_cos = 0
  lat%ele(i)%rf_wake%lr%t_ref = 0
enddo

end subroutine zero_lr_wakes_in_lat

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine lr_wake_add_to (ele, t_ref, orbit, charge)
!
! Subroutine to add to the existing long-range wake the contribution from
! a passing (macro)particle.
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele     -- Ele_struct: Element with wakes.
!   t_ref   -- Real(rp): Time of the reference particle.
!   orbit   -- Coord_struct: Starting coords.
!   charge  -- Real(rp): Charge of passing (macro)particle.
!
! Output:
!   ele     -- Ele_struct: Element with wakes.
!     %rf_wake%lr(:)%b_sin -- Non-skew sin-like wake components.
!     %rf_wake%lr(:)%b_cos -- Non-skew cos-like wake components.
!     %rf_wake%lr(:)%a_sin -- Non-skew sin-like wake components.
!     %rf_wake%lr(:)%a_cos -- Non-skew cos-like wake components.
!     %rf_wake%lr(:)%t_ref -- Set to t_ref.
!+

subroutine lr_wake_add_to (ele, t_ref, orbit, charge)

type (ele_struct), target :: ele
type (coord_struct) orbit
type (rf_wake_lr_struct), pointer :: lr

integer i
real(rp) charge, t_ref, dt, omega, f_exp, ff, c, s, kx, ky
real(rp) c_a, s_a, kxx, exp_shift, a_sin, b_sin

! Check if we have to do any calculations

if (.not. bmad_com%lr_wakes_on) return  
if (.not. associated(ele%rf_wake)) return
  
! Loop over all modes
! We use the following trick: The spatial variation of the normal and skew
! components is the same as the spatial variation of a multipole kick.

! To prevent floating point overflow, the %a and %b factors are shifted 
! to be with respect to lr%t_ref.

do i = 1, size(ele%rf_wake%lr)

  lr => ele%rf_wake%lr(i)

  omega = twopi * lr%freq
  if (lr%freq == 0) omega = twopi * ele%value(rf_frequency$)  ! fundamental mode wake.

  f_exp = omega / (2 * lr%Q)

  if (t_ref /= lr%t_ref) then
    dt = t_ref - lr%t_ref 
    exp_shift = exp(-dt * f_exp) 
    lr%t_ref = t_ref
    lr%b_sin = exp_shift * lr%b_sin
    lr%b_cos = exp_shift * lr%b_cos
    lr%a_sin = exp_shift * lr%a_sin
    lr%a_cos = exp_shift * lr%a_cos
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
  if (lr%freq == 0) dt = dt + ele%value(dphi0$) / omega

  ff = abs(charge) * lr%r_over_q * c_light * exp(dt * f_exp) 

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
! Subroutine lr_wake_apply_kick (ele, t_ref, orbit, charge)
!
! Subroutine to apply the long-range wake kick to a particle.
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele     -- Ele_struct: Element with wakes
!   t_ref   -- Real(rp): Time of the reference particle.
!   orbit   -- Coord_struct: Starting coords of the particle.
!   charge  -- Real(rp): Charge of passing (macro)particle. Needed for self wake.
!
! Output:
!   orbit   -- Coord_struct: coords after the kick.
!+

subroutine lr_wake_apply_kick (ele, t_ref, orbit, charge)

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (rf_wake_lr_struct), pointer :: lr

integer i
real(rp) t_ref, charge, dt, dt_part, omega, f_exp, ff, c, s
real(rp) w_norm, w_skew, kx, ky, k_dum, ff_self, kx_self, ky_self, c_a, s_a

! Check if we have to do any calculations

if (.not. bmad_com%lr_wakes_on) return
if (.not. associated(ele%rf_wake)) return

! Loop over all modes

do i = 1, size(ele%rf_wake%lr)

  lr => ele%rf_wake%lr(i)
  dt_part = -orbit%vec(5) * ele%value(p0c$) / (c_light * ele%value(e_tot$)) 
  dt = t_ref + dt_part - lr%t_ref

  omega = twopi * lr%freq
  if (lr%freq == 0) omega = twopi * ele%value(rf_frequency$)  ! fundamental mode wake.

  f_exp = omega / (2 * lr%Q)
  ff = exp(-dt * f_exp) / ele%value(p0c$) 

  if (lr%freq == 0) dt = dt_part + ele%value(dphi0$) / omega

  c = cos (-dt * omega)
  s = sin (-dt * omega)

  ! Self wake component

  ff_self = abs(charge) * lr%r_over_q * omega / (2 * ele%value(p0c$))

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
! Subroutine sr_table_add_long_kick (ele, leader, charge, follower)
!
! Subroutine to add the component of the gradient loss from the leading particle
! on the following particle.
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele      -- Ele_struct: Element with wakes.
!   leader   -- Coord_struct: Coordinates of the leading particle.
!   charge   -- Real(rp): Charge of leader particle (in Coul).
!   follower -- Coord_struct: Starting coords of particle to kick.
!
! Output:
!   ele      -- Ele_struct: Element with wakes.
!     %value(grad_loss_sr_wake$) -- Adds the effects of the specified leader.
!-

subroutine sr_table_add_long_kick (ele, leader, charge, follower)

implicit none

type (ele_struct) ele
type (coord_struct) leader, follower

real(rp) z, dz, f1, f2, charge, fact

integer iw, n_sr_table

if (.not. bmad_com%sr_wakes_on) return
if (.not. associated(ele%rf_wake)) return

z = follower%vec(5) - leader%vec(5)
n_sr_table = size(ele%rf_wake%sr_table) - 1
dz = ele%rf_wake%sr_table(n_sr_table)%z / n_sr_table

iw = z / dz     ! integer part of z/dz
f2 = z/dz - iw  ! fractional part of z/dz
f1 = 1 - f2

if (iw .lt. 0 .or. iw .gt. ubound(ele%rf_wake%sr_table,1)) return

ele%value(grad_loss_sr_wake$) = ele%value(grad_loss_sr_wake$) + &
      (ele%rf_wake%sr_table(iw)%long*f1 + ele%rf_wake%sr_table(iw+1)%long*f2) * abs(charge) 

end subroutine sr_table_add_long_kick

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr_table_apply_trans_kick (ele, leader, charge, follower)
!
! Subroutine to put in the transverse kick for the short-range wakes.
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele      -- Ele_struct: Element with wakes.
!   leader   -- Coord_struct: Coordinates of the leading particle.
!   charge   -- Real(rp): Charge of leader particle (in Coul).
!   follower -- Coord_struct: Starting coords of particle to kick.
!
! Output:
!   follower -- Coord_struct: coords after the kick.
!+

subroutine sr_table_apply_trans_kick (ele, leader, charge, follower)

implicit none

type (ele_struct) ele
type (coord_struct) leader, follower

real(rp) z, dz, f1, f2, charge, fact
integer iw, n_sr_table

!

if (.not. bmad_com%sr_wakes_on) return
if (.not. associated(ele%rf_wake)) return

z = follower%vec(5) - leader%vec(5)
n_sr_table = size(ele%rf_wake%sr_table) - 1
dz = ele%rf_wake%sr_table(n_sr_table)%z / n_sr_table

iw = z / dz     ! integer part of z/dz
f2 = z/dz - iw  ! fractional part of z/dz
f1 = 1 - f2

fact = (ele%rf_wake%sr_table(iw)%trans*f1 + ele%rf_wake%sr_table(iw+1)%trans*f2) * &
                              abs(charge) * ele%value(l$) / ele%value(p0c$)
follower%vec(2) = follower%vec(2) - fact * leader%vec(1)
follower%vec(4) = follower%vec(4) - fact * leader%vec(3)

end subroutine sr_table_apply_trans_kick

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr_mode_long_wake_add_to (ele, orbit, charge)
!
! Subroutine to add to the existing short-range wake the contribution from
! a passing (macro)particle.
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele     -- Ele_struct: Element with wakes.
!   orbit   -- Coord_struct: Starting coords.
!   charge  -- Real(rp): Charge of passing (macro)particle.
!
! Output:
!   ele     -- Ele_struct: Element with wakes.
!+

subroutine sr_mode_long_wake_add_to (ele, orbit, charge)

type (ele_struct), target :: ele
type (rf_wake_sr_mode_struct), pointer :: sr_mode_long
type (coord_struct) orbit

integer i
real(rp) charge, arg, ff, c, s

! Check if we have to do any calculations

if (.not. bmad_com%sr_wakes_on) return
if (.not. associated(ele%rf_wake)) return

! Add to wake
! The monipole wake does not have any skew components.

do i = 1, size(ele%rf_wake%sr_mode_long)

  sr_mode_long => ele%rf_wake%sr_mode_long(i)

  ff = abs(charge) * sr_mode_long%amp * exp(-orbit%vec(5) * &
                   sr_mode_long%damp) * ele%value(l$) / ele%value(p0c$)

  arg = sr_mode_long%phi - orbit%vec(5) * sr_mode_long%k 
  c = cos (arg)
  s = sin (arg)

  sr_mode_long%b_sin = sr_mode_long%b_sin + ff * c
  sr_mode_long%b_cos = sr_mode_long%b_cos + ff * s

enddo

end subroutine sr_mode_long_wake_add_to

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr_mode_long_wake_apply_kick (ele, charge, orbit)
!
! Subroutine to put in the kick for the short-range wakes.
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele     -- Ele_struct: Element with wakes
!   charge  -- Real(rp): Charge of passing (macro)particle.
!   orbit   -- Coord_struct: Starting coords.
!
! Output:
!   orbit   -- Coord_struct: coords after the kick.
!+

subroutine sr_mode_long_wake_apply_kick (ele, charge, orbit)

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (rf_wake_sr_mode_struct), pointer :: sr_mode_long

integer i
real(rp) arg, ff, c, s, w_norm, charge

! Check if we have to do any calculations

if (.not. bmad_com%sr_wakes_on) return
if (.not. associated(ele%rf_wake)) return

! Loop over all modes

do i = 1, size(ele%rf_wake%sr_mode_long)

  sr_mode_long => ele%rf_wake%sr_mode_long(i)

  ff = exp(orbit%vec(5) * sr_mode_long%damp)

  arg = orbit%vec(5) * sr_mode_long%k 
  c = cos (arg)
  s = sin (arg)

  w_norm = sr_mode_long%b_sin * ff * s + sr_mode_long%b_cos * ff * c
  orbit%vec(6) = orbit%vec(6) - w_norm

  ! Self kick

  orbit%vec(6) = orbit%vec(6) - abs(charge) * sin(sr_mode_long%phi) * &
                         sr_mode_long%amp * ele%value(l$) / (2 * ele%value(p0c$))
enddo

end subroutine sr_mode_long_wake_apply_kick

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr_mode_trans_wake_add_to (ele, orbit, charge)
!
! Subroutine to add to the existing short-range wake the contribution from
! a passing (macro)particle.
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele     -- Ele_struct: Element with wakes.
!   orbit   -- Coord_struct: Starting coords.
!   charge  -- Real(rp): Charge of passing (macro)particle.
!
! Output:
!   ele     -- Ele_struct: Element with wakes.
!+

subroutine sr_mode_trans_wake_add_to (ele, orbit, charge)

type (ele_struct), target :: ele
type (rf_wake_sr_mode_struct), pointer :: sr_mode_trans
type (coord_struct) orbit

integer i
real(rp) charge, arg, ff, c, s

! Check if we have to do any calculations

if (.not. bmad_com%sr_wakes_on) return  
if (.not. associated(ele%rf_wake)) return

! Add to wake
! The monipole wake does not have any skew components.

do i = 1, size(ele%rf_wake%sr_mode_trans)

  sr_mode_trans => ele%rf_wake%sr_mode_trans(i)

  ff = abs(charge) * sr_mode_trans%amp * exp(-orbit%vec(5) * sr_mode_trans%damp) * &
                                           ele%value(l$) / ele%value(p0c$)

  arg =  sr_mode_trans%phi - orbit%vec(5) * sr_mode_trans%k 
  c = cos (arg)
  s = sin (arg)

  sr_mode_trans%b_sin = sr_mode_trans%b_sin + ff * orbit%vec(1) * c
  sr_mode_trans%b_cos = sr_mode_trans%b_cos + ff * orbit%vec(1) * s
  sr_mode_trans%a_sin = sr_mode_trans%a_sin + ff * orbit%vec(3) * c
  sr_mode_trans%a_cos = sr_mode_trans%a_cos + ff * orbit%vec(3) * s

enddo

end subroutine sr_mode_trans_wake_add_to

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr_mode_trans_wake_apply_kick (ele, orbit)
!
! Subroutine to put in the kick for the short-range wakes
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele     -- Ele_struct: Element with wakes
!   orbit   -- Coord_struct: Starting coords.
!
! Output:
!   orbit   -- Coord_struct: coords after the kick.
!+

subroutine sr_mode_trans_wake_apply_kick (ele, orbit)

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (rf_wake_sr_mode_struct), pointer :: sr_mode_trans

integer i
real(rp) arg, ff, c, s, w_norm, w_skew

! Check if we have to do any calculations

if (.not. bmad_com%sr_wakes_on) return
if (.not. associated(ele%rf_wake)) return

! Loop over all modes

do i = 1, size(ele%rf_wake%sr_mode_trans)

  sr_mode_trans => ele%rf_wake%sr_mode_trans(i)

  ff = exp(orbit%vec(5) * sr_mode_trans%damp)

  arg = orbit%vec(5) * sr_mode_trans%k 
  c = cos (arg)
  s = sin (arg)

  w_norm = sr_mode_trans%b_sin * ff * s + sr_mode_trans%b_cos * ff * c
  w_skew = sr_mode_trans%a_sin * ff * s + sr_mode_trans%a_cos * ff * c

  orbit%vec(2) = orbit%vec(2) - w_norm
  orbit%vec(4) = orbit%vec(4) - w_skew

enddo

end subroutine sr_mode_trans_wake_apply_kick

end module
