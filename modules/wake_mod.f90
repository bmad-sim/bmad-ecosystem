module wake_mod

use bmad_struct
use bmad_interface
use multipole_mod, only: ab_multipole_kick

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine lr_wake_apply_kick (ele, s_ref, orbit)
!
! Subroutine to put in the kick for the long-range wakes
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele     -- Ele_struct: Element with wakes
!   s_ref   -- Real(rp): S position of the reference particle.
!   orbit   -- Coord_struct: Starting coords.
!
! Output:
!   orbit   -- Coord_struct: coords after the kick.
!+

subroutine lr_wake_apply_kick (ele, s_ref, orbit)

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (lr_wake_struct), pointer :: lr

integer i, n_mode
real(rp) s_ref, ds, k, f, f_exp, ff, c, s, w_norm, w_skew, kx, ky

! Check if we have to do any calculations

if (.not. bmad_com%lr_wakes_on) return  
n_mode = lr_wake_array_size(ele)

! Loop over all modes

do i = 1, n_mode

  lr => ele%wake%lr(i)
  ds = (s_ref + orbit%vec(5)) - lr%s_ref  ! Note: ds < 0

  k = twopi * lr%freq / c_light
  f = ele%value(l$) * lr%r_over_q * (k**2 / 2)
  f_exp = k / (2 * lr%Q)
  ff = f * exp(ds * f_exp) / ele%value(beam_energy$) 

  c = cos (ds * k)
  s = sin (ds * k)

! longitudinal kick

  w_norm = lr%norm_sin * ff * (f_exp * s + k * c) + &
           lr%norm_cos * ff * (f_exp * c - k * s)

  w_skew = lr%skew_sin * ff * (f_exp * s + k * c) + &
           lr%skew_cos * ff * (f_exp * c - k * s)


  call ab_multipole_kick (w_skew, w_norm, lr%m, orbit, kx, ky)
  orbit%vec(6) = orbit%vec(6) + kx + ky

! transverse kick

  if (lr%m == 0) cycle

  w_norm = lr%norm_sin * ff * s + lr%norm_cos * ff * c
  w_skew = lr%skew_sin * ff * s + lr%skew_cos * ff * c

  call ab_multipole_kick (w_skew, w_norm, lr%m-1, orbit, kx, ky)

  orbit%vec(2) = orbit%vec(2) + kx
  orbit%vec(4) = orbit%vec(4) + ky

enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine lr_wake_add_to (ele, s_ref, orbit, charge)
!
! Subroutine to add to the existing long-range wake the contribution from
! a passing (macro)particle.
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele     -- Ele_struct: Element with wakes
!   s_ref   -- Real(rp): S position of the reference particle.
!   orbit   -- Coord_struct: Starting coords.
!   charge  -- Real(rp): Charge of passing (macro)particle.
!
! Output:
!   ele     -- Coord_struct: coords after the kick.
!+

subroutine lr_wake_add_to (ele, s_ref, orbit, charge)

type (ele_struct), target :: ele
type (coord_struct) orbit
type (lr_wake_struct), pointer :: lr

integer i, n_mode
real(rp) charge, s_ref, ds, k, f, f_exp, ff, c, s, kx, ky

! Check if we have to do any calculations

if (.not. bmad_com%lr_wakes_on) return  
n_mode = lr_wake_array_size(ele)

!

do i = 1, n_mode

  lr => ele%wake%lr(i)
  ds = (s_ref + orbit%vec(5)) - lr%s_ref  ! Note: ds < 0

  k = twopi * lr%freq / c_light
  f_exp = k / (2 * lr%Q)
  ff = exp(ds * f_exp)

  c = cos (ds * k)
  s = sin (ds * k)

  call ab_multipole_kick (0.0_rp, -1.0_rp, lr%m, orbit, kx, ky)

  lr%norm_sin = lr%norm_sin + charge * kx * c
  lr%norm_cos = lr%norm_cos - charge * kx * s
  lr%skew_sin = lr%skew_sin + charge * ky * c
  lr%skew_cos = lr%skew_cos - charge * ky * s

enddo

end subroutine

end module
