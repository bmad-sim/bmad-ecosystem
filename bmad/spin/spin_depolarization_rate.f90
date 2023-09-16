!+
! Function spin_depolarization_rate (branch, match_info, rad_int_by_ele) result (depol_rate)
!
! Routine to calculate the depolarization rate of a beam in the linear model where the dependence
! of dn_isf/dDelta on transverse position is ignored. 
!
! See the writeup by Barber & Ripken in the Handbook of Accelerator Physics and Engineering.
!
! Input:
!   branch          -- branch_struct: Lattice branch the beam is going through.
!   match_info(0:)  -- spin_matching_struct:
!     %dn_dpz           -- ISR derivative.
!   rad_int_by_ele  -- rad_int_all_ele_struct: Element-by-element radiation integrals.
!     %i3               -- I3 radiation integral.
!
! Output:
!   depol_rate      -- real(rp): Depolarization rate (1/sec). Will be positive.
!-

function spin_depolarization_rate (branch, match_info, rad_int_by_ele) result (depol_rate)

use bmad_interface, dummy => spin_depolarization_rate

implicit none

type (branch_struct), target :: branch
type (ele_struct), pointer :: ele
type (spin_matching_struct) match_info(0:)
type (rad_int_all_ele_struct) rad_int_by_ele

real(rp) depol_rate
real(rp), parameter :: factor = 55.0_rp * sqrt(3.0_rp) * classical_radius_factor * h_bar_planck / 144.0_rp
integer ie, ib

!

depol_rate = 0
ib = branch%ix_branch

do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)
  depol_rate = depol_rate + 0.5_rp * rad_int_by_ele%branch(ib)%ele(ie)%i3 * ele%value(l$) * &
                  (norm2(match_info(ie-1)%dn_dpz) + norm2(match_info(ie)%dn_dpz))
                    
enddo

depol_rate = depol_rate * factor * gamma_ref(branch%ele(0))**5 / branch%param%total_length

end function
