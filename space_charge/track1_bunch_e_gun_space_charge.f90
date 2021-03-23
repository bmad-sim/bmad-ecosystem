!+
! Subroutine track1_bunch_e_gun_space_charge (bunch_start, ele, bunch_end, err)
!
! Subroutine to track a bunch of particles through an element.
!
! Input:
!   bunch_start  -- bunch_struct: Starting bunch position.
!   ele          -- Ele_struct: element to track through. Must be part of a lattice.
!
! Output:
!   bunch_end -- Bunch_struct: Ending bunch position.
!   err       -- Logical: Set true if there is an error. EG: Too many particles lost for a CSR calc.
!-

subroutine track1_bunch_e_gun_space_charge (bunch_start, ele, bunch_end, err)

use bmad_routine_interface, dummy => track1_bunch_e_gun_space_charge

implicit none

type (bunch_struct) bunch_start, bunch_end
type (ele_struct), target :: ele
type (lat_struct), pointer :: lat

logical err

character(*), parameter :: r_name = 'track1_bunch_e_gun_space_charge'

!

bunch_end = bunch_start
bunch_end%particle%state = pre_born$
lat => ele%branch%lat

! Question: When a particle is born, how to handle fact that the birth time is not comensurate with the time steps.
! Question: When to end tracking?
! Question: How to handle slices/super_slaves

do
  ! call track_bunch_time(lat, bunch, t_end, dt_step)

enddo

end subroutine track1_bunch_e_gun_space_charge

