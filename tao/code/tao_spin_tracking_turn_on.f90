!+
! Subroutine tao_spin_tracking_turn_on()
!
! Routine to turn on spin tracking.
!-

subroutine tao_spin_tracking_turn_on()

use tao_interface, dummy => tao_spin_tracking_turn_on
implicit none

logical ok
character(*), parameter :: r_name = 'tao_spin_tracking_turn_on'

!

if (bmad_com%spin_tracking_on) return

bmad_com%spin_tracking_on = .true.
s%u%calc%lattice = .true.
call tao_lattice_calc(ok)

call out_io (s_info$, r_name, 'Note: Setting bmad_com%spin_tracking_on to True for spin tracking.')

end subroutine
