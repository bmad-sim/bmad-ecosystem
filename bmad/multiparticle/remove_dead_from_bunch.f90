!+
! Subroutine remove_dead_from_bunch(bunch_in, bunch_out)
!
! Routine to return a bunch with dead particles (not alive and not pre_born) removed
!
! Input:
!   bunch_in        -- bunch_struct: Input bunch with alive and dead particles.
!
! Output:
!   bunch_out       -- bunch_struct: Output bunch with only alive and pre_born particles.
!                       Note: bunch_out can be the same actual argument as bunch_in.
!-

subroutine remove_dead_from_bunch(bunch_in, bunch_out)

use bmad_struct

implicit none

type (bunch_struct) bunch_in, bunch_out
type (coord_struct), allocatable :: p(:)
integer ip, n

!

p = bunch_in%particle
n = count(p%state == alive$) + count(p%state == pre_born$)
if (allocated(bunch_out%particle)) deallocate(bunch_out%particle)
allocate(bunch_out%particle(n))

n = 0
do ip = 1, size(p)
  if (p(ip)%state /= alive$ .and. p(ip)%state /= pre_born$) cycle
  n = n + 1
  bunch_out%particle(n) = p(ip)
enddo

bunch_out%charge_live = bunch_in%charge_live
bunch_out%charge_tot  = bunch_in%charge_live
bunch_out%z_center    = bunch_in%z_center
bunch_out%t_center    = bunch_in%t_center
bunch_out%ix_ele      = bunch_in%ix_ele
bunch_out%ix_bunch    = bunch_in%ix_bunch
bunch_out%ix_turn     = bunch_in%ix_turn
bunch_out%n_live      = bunch_in%n_live
bunch_out%drift_between_t_and_s = bunch_in%drift_between_t_and_s

end subroutine
