!+
! Subroutine track1_preprocess (start_orb, ele, param, err_flag, finished, radiation_included, track)
!
! Dummy routine for pre-processing at the start of the track1 routine.
!
! Also see:
!   track1_postprocess
!   track1_custom
!
! The radiation_included argument should be set to True if this routine (or a modified version of track1_custom)
! takes into account radiation damping and/or excitation. This will prevent track1 from calling track1_radiation.
! Note: If symp_lie_bmad is being called by this routine, symp_lie_bmad does take into account radiation effects.
! 
! General rule: Your code may NOT modify any argument that is not listed as an output agument below.
!
! Input:
!   start_orb  -- coord_struct: Starting position at the beginning of ele.
!   ele        -- ele_struct: Element.
!   param      -- lat_param_struct: Lattice parameters.
!
! Output:
!   start_orb   -- coord_struct: Modified starting position for track1 to use.
!   err_flag    -- logical: Set true if there is an error. False otherwise.
!   finished    -- logical: When set True, track1 will halt further processing and use the modified 
!                     start_orb as the final end position of the particle.
!   radiation_included
!               -- logical: Should be set True if radiation damping/excitation is included in the tracking.
!   track       -- track_struct, optional: Structure holding the track information if the 
!                    tracking method does tracking step-by-step.
!-

subroutine track1_preprocess (start_orb, ele, param, err_flag, finished, radiation_included, track)

use bmad, except_dummy => track1_preprocess

implicit none

type (coord_struct) :: start_orb
type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (track_struct), optional :: track

logical err_flag, finished, radiation_included

character(*), parameter :: r_name = 'track1_preprocess'

!

err_flag = .false.

end subroutine
