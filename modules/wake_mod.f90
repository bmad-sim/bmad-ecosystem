module wake_mod

use bmad_struct
use bmad_interface
use multipole_mod, only: ab_multipole_kick

contains

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
!   ele     -- Ele_struct: Element with wakes.
!   s_ref   -- Real(rp): S position of the reference particle.
!   orbit   -- Coord_struct: Starting coords.
!   charge  -- Real(rp): Charge of passing (macro)particle.
!
! Output:
!   ele     -- Ele_struct: Element with wakes.
!     %wake%lr(:)%norm_sin -- Non-skew sin-like wake components.
!     %wake%lr(:)%norm_cos -- Non-skew cos-like wake components.
!     %wake%lr(:)%skew_sin -- Non-skew sin-like wake components.
!     %wake%lr(:)%skew_cos -- Non-skew cos-like wake components.
!+

subroutine lr_wake_add_to (ele, s_ref, orbit, charge)

type (ele_struct), target :: ele
type (coord_struct) orbit
type (lr_wake_struct), pointer :: lr

integer i
real(rp) charge, s_ref, ds, k, f_exp, ff, c, s, kx, ky
real(rp) c_a, s_a, kxx

! Check if we have to do any calculations

if (.not. bmad_com%lr_wakes_on) return  
if (.not. associated(ele%wake)) return
  
! Loop over all modes
! We use the following trick: The spatial variation of the normal and skew
! components is the same as the spatial variation of a multipole kick.

do i = 1, size(ele%wake%lr)

  lr => ele%wake%lr(i)
  ds = (s_ref + orbit%vec(5)) ! Note: ds < 0

  k = twopi * lr%freq / c_light
  f_exp = k / (2 * lr%Q)
  ff = charge * lr%r_over_q * (c_light / 2) * exp(-ds * f_exp) / &
                                                    ele%value(beam_energy$) 

  c = cos (ds * k)
  s = sin (ds * k)

  call ab_multipole_kick (0.0_rp, ff, lr%m, orbit, kx, ky)

  if (lr%polarized) then
    c_a = cos(twopi*lr%angle); s_a = sin(twopi*lr%angle)
    kxx = kx
    kx = kxx * c_a * c_a + ky * s_a * c_a
    ky = kxx * c_a * s_a + ky * s_a * s_a
  endif

  lr%norm_sin = lr%norm_sin - kx * c
  lr%norm_cos = lr%norm_cos + kx * s
  lr%skew_sin = lr%skew_sin - ky * c
  lr%skew_cos = lr%skew_cos + ky * s

enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine lr_wake_apply_kick (ele, s_ref, orbit)
!
! Subroutine to apply the long-range wake kick to a particle.
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele     -- Ele_struct: Element with wakes
!   s_ref   -- Real(rp): S position of the reference particle.
!   orbit   -- Coord_struct: Starting coords of the particle.
!
! Output:
!   orbit   -- Coord_struct: coords after the kick.
!+

subroutine lr_wake_apply_kick (ele, s_ref, orbit)

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (lr_wake_struct), pointer :: lr

integer i
real(rp) s_ref, ds, k, f_exp, ff, c, s, w_norm, w_skew, kx, ky, k_dum

! Check if we have to do any calculations

if (.not. bmad_com%lr_wakes_on) return
if (.not. associated(ele%wake)) return

! Loop over all modes

do i = 1, size(ele%wake%lr)

  lr => ele%wake%lr(i)
  ds = s_ref + orbit%vec(5)  ! Note: ds < 0

  k = twopi * lr%freq / c_light
  f_exp = k / (2 * lr%Q)
  ff = exp(ds * f_exp)

  c = cos (ds * k)
  s = sin (ds * k)

! longitudinal kick

  w_norm = lr%norm_sin * ff * (f_exp * s + k * c) + &
           lr%norm_cos * ff * (f_exp * c - k * s)

  w_skew = lr%skew_sin * ff * (f_exp * s + k * c) + &
           lr%skew_cos * ff * (f_exp * c - k * s)

  call ab_multipole_kick (0.0_rp, w_norm, lr%m, orbit, kx, k_dum)
  call ab_multipole_kick (0.0_rp, w_skew, lr%m, orbit, k_dum, ky)

  orbit%vec(6) = orbit%vec(6) + kx + ky

! transverse kick

  if (lr%m == 0) cycle

  w_norm = lr%norm_sin * ff * s + lr%norm_cos * ff * c
  w_skew = lr%skew_sin * ff * s + lr%skew_cos * ff * c

  call ab_multipole_kick (w_skew, w_norm, lr%m-1, orbit, kx, ky)

  orbit%vec(2) = orbit%vec(2) + lr%m * kx
  orbit%vec(4) = orbit%vec(4) + lr%m * ky

enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr1_apply_kick (ele, leader, charge, follower)
!
! Subroutine to put in the kick for the short-range wakes.
!
! Modules needed:
!   use wake_mod
!
! Input:
!   ele      -- Ele_struct: Element with wakes.
!   leader   -- Coord_struct: Coordinates of the leading particle
!   charge   -- Real(rp): Charge of leader particle.
!   follower -- Coord_struct: Starting coords of particle to kick.
!
! Output:
!   follower -- Coord_struct: coords after the kick.
!+

subroutine sr1_apply_kick (ele, leader, charge, follower)

implicit none

type (ele_struct) ele
type (coord_struct) leader, follower

real(rp) z, dz, f1, f2, charge, fact
integer iw, n_sr1

!

z = follower%vec(6) - leader%vec(6)
n_sr1 = size(ele%wake%sr1) - 1
dz = ele%wake%sr1(n_sr1)%z / n_sr1

iw = z / dz     ! integer part of z/dz
f2 = z/dz - iw  ! fractional part of z/dz
f1 = 1 - f2

fact = (ele%wake%sr1(iw)%trans*f1 + ele%wake%sr1(iw+1)%trans*f2) * &
                              charge * ele%value(l$) / ele%value(p0c$)
follower%vec(2) = follower%vec(2) + fact * charge * leader%vec(1)
follower%vec(4) = follower%vec(4) + fact * charge * leader%vec(3)

fact = (ele%wake%sr1(iw)%long*f1 + ele%wake%sr1(iw+1)%long*f2) * &
                              charge * ele%value(l$) / ele%value(p0c$)
follower%vec(6) = follower%vec(6) + fact * charge / ele%value(p0c$)

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr2_long_wake_add_to (ele, orbit, charge)
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

subroutine sr2_long_wake_add_to (ele, orbit, charge)

type (ele_struct), target :: ele
type (sr2_wake_struct), pointer :: sr2_long
type (coord_struct) orbit

integer i
real(rp) charge, arg, ff, c, s, kx, ky

! Check if we have to do any calculations

if (.not. bmad_com%sr_wakes_on) return  
if (.not. associated(ele%wake)) return

! Add to wake
! The monipole wake does not have any skew components.

do i = 1, size(ele%wake%sr2_long)

  sr2_long => ele%wake%sr2_long(i)

  ff = sr2_long%amp * exp(-orbit%vec(5) * sr2_long%damp) / &
                                                  ele%value(beam_energy$) 

  arg = orbit%vec(5) * sr2_long%k + sr2_long%phi
  c = cos (arg)
  s = sin (arg)

  call ab_multipole_kick (0.0_rp, ff, 0, orbit, kx, ky)

  sr2_long%norm_sin = sr2_long%norm_sin - charge * ff * c
  sr2_long%norm_cos = sr2_long%norm_cos + charge * ff * s

enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr2_trans_wake_add_to (ele, orbit, charge)
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

subroutine sr2_trans_wake_add_to (ele, orbit, charge)

type (ele_struct), target :: ele
type (sr2_wake_struct), pointer :: sr2_trans
type (coord_struct) orbit

integer i
real(rp) charge, arg, ff, c, s, kx, ky

! Check if we have to do any calculations

if (.not. bmad_com%sr_wakes_on) return  
if (.not. associated(ele%wake)) return

! Add to wake
! The monipole wake does not have any skew components.

do i = 1, size(ele%wake%sr2_trans)

  sr2_trans => ele%wake%sr2_trans(i)

  ff = sr2_trans%amp * exp(-orbit%vec(5) * sr2_trans%damp) / &
                                                     ele%value(beam_energy$) 

  arg = orbit%vec(5) * sr2_trans%k + sr2_trans%phi
  c = cos (arg)
  s = sin (arg)

  call ab_multipole_kick (0.0_rp, ff, 1, orbit, kx, ky)

  sr2_trans%norm_sin = sr2_trans%norm_sin - charge * kx * c
  sr2_trans%norm_cos = sr2_trans%norm_cos + charge * kx * s
  sr2_trans%skew_sin = sr2_trans%skew_sin - charge * ky * c
  sr2_trans%skew_cos = sr2_trans%skew_cos + charge * ky * s

enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr2_long_wake_apply_kick (ele, orbit)
!
! Subroutine to put in the kick for the short-range wakes.
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

subroutine sr2_long_wake_apply_kick (ele, orbit)

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (sr2_wake_struct), pointer :: sr2_long

integer i
real(rp) arg, ff, c, s, w_norm

! Check if we have to do any calculations

if (.not. bmad_com%sr_wakes_on) return
if (.not. associated(ele%wake)) return

! Loop over all modes

do i = 1, size(ele%wake%sr2_long)

  sr2_long => ele%wake%sr2_long(i)

  ff = exp(orbit%vec(5) * sr2_long%damp)

  arg = orbit%vec(5) * sr2_long%k 
  c = cos (arg)
  s = sin (arg)

  w_norm = sr2_long%norm_sin * ff * s + sr2_long%norm_cos * ff * c
  orbit%vec(6) = orbit%vec(6) - w_norm

enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr2_long_self_wake_apply_kick (ele, charge, orbit)
!
! Subroutine to put in the kick for the short-range wakes
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

subroutine sr2_long_self_wake_apply_kick (ele, charge, orbit)

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (sr2_wake_struct), pointer :: sr2_long

integer i
real(rp) k, c, s, w_norm, charge

! Check if we have to do any calculations

if (.not. bmad_com%sr_wakes_on) return
if (.not. associated(ele%wake)) return

! Loop over all modes

do i = 1, size(ele%wake%sr2_long)
  sr2_long => ele%wake%sr2_long(i)
  orbit%vec(6) = orbit%vec(6) - charge * sr2_long%amp / 2
enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr2_trans_wake_apply_kick (ele, orbit)
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

subroutine sr2_trans_wake_apply_kick (ele, orbit)

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (sr2_wake_struct), pointer :: sr2_trans

integer i
real(rp) arg, ff, c, s, w_norm, w_skew, kx, ky

! Check if we have to do any calculations

if (.not. bmad_com%sr_wakes_on) return
if (.not. associated(ele%wake)) return

! Loop over all modes

do i = 1, size(ele%wake%sr2_trans)

  sr2_trans => ele%wake%sr2_trans(i)

  ff = exp(orbit%vec(5) * sr2_trans%damp)

  arg = orbit%vec(5) * sr2_trans%k 
  c = cos (arg)
  s = sin (arg)

  w_norm = sr2_trans%norm_sin * ff * s + sr2_trans%norm_cos * ff * c
  w_skew = sr2_trans%skew_sin * ff * s + sr2_trans%skew_cos * ff * c

  call ab_multipole_kick (w_skew, w_norm, 0, orbit, kx, ky)

  orbit%vec(2) = orbit%vec(2) + kx
  orbit%vec(4) = orbit%vec(4) + ky

enddo

end subroutine

end module
