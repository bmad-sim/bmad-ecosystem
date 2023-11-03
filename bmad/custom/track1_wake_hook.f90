!+
! Subroutine track1_wake_hook (bunch, ele, finished)
!
! Routine that can be customized for tracking through a wake.
! To use, see the Bmad manual.
!
! Input:
!   bunch       -- bunch_struct: Bunch of particles.
!   ele         -- ele_struct: Element being tracked through.
!
! Output:
!   bunch       -- bunch_struct: Bunch of particles with wake applied.
!   ele         -- ele_struct: Element with updated wake.
!   finished    -- logical: When set True, the standard wake code will not be called.
!-

subroutine track1_wake_hook (bunch, ele, finished)

use bmad_interface, dummy => track1_wake_hook

implicit none

type (bunch_struct) bunch
type (ele_struct) ele

logical finished

!

finished = .false.   ! Must set this

end subroutine track1_wake_hook
