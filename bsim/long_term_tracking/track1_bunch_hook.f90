!+
! Subroutine track1_bunch_hook (bunch, ele, err, centroid, direction, finished, bunch_track)
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

subroutine track1_bunch_hook (bunch, ele, err, centroid, direction, finished, bunch_track)

use lt_tracking_mod, dummy => track1_bunch_hook

implicit none

type (bunch_struct), target :: bunch
type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0
type (bunch_track_struct), optional :: bunch_track
type (coord_struct), optional :: centroid(0:)
type (coord_struct), pointer :: orb

real(rp) t, r

integer, optional :: direction
integer ip, ir, ie, n, iv

logical err, finished

! The reason to separate this routine from ltt_track1_bunch_hook is to allow for the possibility
! to add custom functionality by overriding track1_bunch_hook.

call ltt_track1_bunch_hook (bunch, ele, err, centroid, direction, finished, bunch_track)

end subroutine
