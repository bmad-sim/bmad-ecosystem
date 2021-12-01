!+
! Subroutine track1_bunch_hook (bunch_start, ele, bunch_end, err, centroid, direction, finished)
!
! Routine that can be customized for tracking a bunch through a single element.
!
! Input:
!   bunch_start   -- Bunch_struct: Starting bunch position.
!   ele          -- Ele_struct: Element to track through.
!   centroid(0:) -- coord_struct, optional: Approximate centroid orbit. Only needed if CSR is on.
!                     Hint: Calculate this before bunch tracking by tracking a single particle.
!   direction    -- integer, optional: +1 (default) -> Track forward, -1 -> Track backwards.
!
! Output:
!   bunch_end    -- bunch_struct: Ending bunch position.
!   err         -- Logical: Set true if there is an error. 
!                    EG: Too many particles lost for a CSR calc.
!   finished    -- logical: When set True, the standard track1_bunch code will not be called.
!-

subroutine track1_bunch_hook (bunch_start, ele, bunch_end, err, centroid, direction, finished)

use bmad, dummy => track1_bunch_hook

implicit none

type (bunch_struct), target :: bunch_start
type (bunch_struct), target :: bunch_end
type (ele_struct), target :: ele
type (coord_struct), optional :: centroid(0:)

integer, optional :: direction
logical err, finished

!

finished = .false.

end subroutine
