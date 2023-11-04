!+
! Subroutine track1_bunch_hook (bunch, ele, err, centroid, direction, finished, bunch_track)
!
! Prototype routine that can be customized for tracking a bunch through a single element.
!
! Input:
!   bunch          -- Bunch_struct: Starting bunch position.
!   ele           -- Ele_struct: Element to track through.
!   centroid(0:)  -- coord_struct, optional: Approximate centroid orbit. Only needed if CSR is on.
!                      Hint: Calculate this before bunch tracking by tracking a single particle.
!   direction     -- integer, optional: +1 (default) -> Track forward, -1 -> Track backwards.
!   bunch_track   -- bunch_track_struct, optional: Existing tracks. If bunch_track%n_pt = -1 then
!                        Overwrite any existing track.
!
! Output:
!   bunch       -- bunch_struct: Ending bunch position.
!   err         -- Logical: Set true if there is an error. 
!                    EG: Too many particles lost for a CSR calc.
!   finished    -- logical: When set True, the standard track1_bunch code will not be called.
!   bunch_track -- bunch_track_struct, optional: track information if the tracking method does
!                        tracking step-by-step. When tracking through multiple elements, the 
!                        trajectory in an element is appended to the existing trajectory. 
!-

subroutine track1_bunch_hook (bunch, ele, err, centroid, direction, finished, bunch_track)

use bmad

implicit none

type (bunch_struct), target :: bunch
type (ele_struct), target :: ele
type (coord_struct), optional :: centroid(0:)
type (bunch_track_struct), optional :: bunch_track

integer, optional :: direction
logical err, finished

!

finished = .false.

end subroutine
