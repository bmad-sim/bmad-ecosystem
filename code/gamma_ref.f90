!+
! Function gamma_ref(ele) result (gamma)
! 
! Routine to compute the referece relativisitic gamma factor for an element.
! This is the gamma at the exit end if the reference energy is changing (EG lcavity).
!
! Input:
!   ele       -- ele_struct: Element to evaluate at.
!
! Output:
!   gamma     -- real(rp): Relativistic gamma factor Energy/mass*c^2.
!-

function gamma_ref(ele) result (gamma)

use bmad_struct
implicit none
type (ele_struct) ele
real(rp) gamma

!

gamma = ele%value(e_tot$) / mass_of(ele%ref_species)

end function
