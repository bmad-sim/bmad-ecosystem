!+
! Subroutine track1_custom (start_orb, ele, param, end_orb, track, err_flag, entry_pt, finished)
!
! Dummy routine for custom tracking. 
! This routine needs to be replaced for a custom calculation.
! If not replaced and this routine is called, this routine will generate an error message.
!
! This routine is potentially called twice by track1. 
! The entry_pt argument indicates from which point in track1 this routine is being called.
! 
! General rule: Your code may NOT modify any argument that is not listed as
! an output agument below.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start_orb  -- coord_struct: Starting position.
!   ele        -- ele_struct: Element.
!   param      -- lat_param_struct: Lattice parameters.
!   entry_pt   -- integer: Flag indicating from which point in track1 this routine is called.
!                   Possibilities are: entry_pt1$, entry_pt2$.
!
! Output:
!   end_orb     -- coord_struct: End position.
!   track       -- track_struct, optional: Structure holding the track information if the 
!                    tracking method does tracking step-by-step.
!   err_flag    -- logical: Set true if there is an error. False otherwise.
!   finished    -- logical: When set True, track1 will halt processing and return to its calling routine.
!-

subroutine track1_custom (start_orb, ele, param, end_orb, track, err_flag, entry_pt, finished)

use bmad_interface, except_dummy => track1_custom

implicit none

type (coord_struct) :: start_orb
type (coord_struct) :: end_orb
type (ele_struct) :: ele
type (lat_param_struct) :: param
type (track_struct), optional :: track

integer entry_pt
logical err_flag, finished

character(32) :: r_name = 'track1_custom'

!

finished = .false.

call out_io (s_fatal$, r_name, 'THIS DUMMY ROUTINE SHOULD NOT HAVE BEEN CALLED IN THE FIRST PLACE.')
err_flag = .true.

! Remember to also set end_orb%t

end_orb = start_orb
end_orb%s = ele%s 

end subroutine
