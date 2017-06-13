!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_beam_hook (beam_start, lat, ele, beam_end, err, centroid, direction)
!
! Routine that can be customized for tracking a beam through a single element.
!
! Input:
!   beam_start   -- Beam_struct: Starting beam position.
!   lat          -- lat_struct: Lattice containing element to be tracked through.
!   ele          -- Ele_struct: Element to track through.
!   centroid(0:) -- coord_struct, optional: Approximate centroid orbit. Only needed if CSR is on.
!                     Hint: Calculate this before beam tracking by tracking a single particle.
!   direction    -- integer, optional: +1 (default) -> Track forward, -1 -> Track backwards.
!
! Output:
!   beam_end    -- beam_struct: Ending beam position.
!   err         -- Logical: Set true if there is an error. 
!                    EG: Too many particles lost for a CSR calc.
!   finished    -- logical: When set True, the standard track1_beam code will not be called.
!-

subroutine track1_beam_hook (beam_start, lat, ele, beam_end, err, centroid, direction, finished)

use bmad, dummy => track1_beam_hook

implicit none

type (beam_struct) beam_start
type (beam_struct) :: beam_end
type (lat_struct) :: lat
type (ele_struct) ele
type (coord_struct), optional :: centroid(0:)

integer, optional :: direction
logical err, finished

!

finished = .false.

end subroutine
