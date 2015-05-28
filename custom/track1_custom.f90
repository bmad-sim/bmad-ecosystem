!+
! Subroutine track1_custom (start_orb, ele, param, end_orb, err_flag, finished, track)
!
! Dummy routine for custom tracking. 
! This routine needs to be replaced for a custom calculation.
! If not replaced and this routine is called, this routine will generate an error message.
!
! Also see:
!   track1_preprocess
!   track1_postprocess
!
! If this routine takes into account radiation damping and/or excitation when bmad_com%radiation_damping_on 
! and/or bmad_com%radiation_fluctuations_on is True, a custom version of track1_preprocess should be 
! constructed to set its radiation_included argument to True.
! If not, the track1 routine will use track1_radiation to include the radiation effects.
! Note: If this routine calles symp_lie_bmad, the symp_lie_bmad routine does take into account radiation effects.
! 
! General rule: Your code may NOT modify any argument that is not listed as an output agument below.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start_orb  -- coord_struct: Starting position.
!   ele        -- ele_struct: Element.
!   param      -- lat_param_struct: Lattice parameters.
!
! Output:
!   end_orb     -- coord_struct: End position.
!   err_flag    -- logical: Set true if there is an error. False otherwise.
!   finished    -- logical: When set True, track1 will halt processing and return to its calling routine.
!   track       -- track_struct, optional: Structure holding the track information if the 
!                    tracking method does tracking step-by-step.
!-

subroutine track1_custom (start_orb, ele, param, end_orb, err_flag, finished, track)

use bmad, except_dummy => track1_custom

implicit none

type (coord_struct) :: start_orb
type (coord_struct) :: end_orb
type (ele_struct) :: ele
type (lat_param_struct) :: param
type (track_struct), optional :: track

logical err_flag, finished

character(*), parameter :: r_name = 'track1_custom'

!

call out_io (s_fatal$, r_name, 'THIS DUMMY ROUTINE SHOULD NOT HAVE BEEN CALLED IN THE FIRST PLACE.')
err_flag = .true.

! Remember to also set end_orb%t

end_orb = start_orb
end_orb%s = ele%s 

end subroutine
