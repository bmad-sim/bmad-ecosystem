!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine kill_ptc_genfield (ptc_genfield)
!
! Subroutine to kill a ptc_genfield.
!
! Input:
!   ptc_genfield -- Genfield, pointer: ptc_genfield to kill.
!
! Output:
!   ptc_genfield -- Genfield, pointer: Killed ptc_genfield.
!-

subroutine kill_ptc_genfield (ptc_genfield)

use bmad_routine_interface, dummy => kill_ptc_genfield

use tpsalie_analysis, only: kill 

implicit none

type (genfield), pointer :: ptc_genfield

!

if (associated(ptc_genfield)) then
  call kill (ptc_genfield)
  deallocate (ptc_genfield)
endif

end subroutine kill_ptc_genfield

