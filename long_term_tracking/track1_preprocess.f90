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
use lt_tracking_mod

implicit none

type (coord_struct) :: start_orb
type (ele_struct) :: ele
type (lat_param_struct) :: param
type (track_struct), optional :: track

real(rp) r
integer i, ie
logical err_flag, finished, radiation_included

character(*), parameter :: r_name = 'track1_preprocess'

! If bunch tracking, ramper bookkeeping is handled by track1_bunch_hook.

err_flag = .false.
if (ltt_params_global%tracking_method /= 'SINGLE') return 

do i = 1, size(ltt_com_global%ix_ramper)
  ie = ltt_com_global%ix_ramper(i)
  if (ie == -1) exit
  ltt_com_global%lat%ele(ie)%control%var(1)%value = start_orb%t + 0.5_rp * ele%value(delta_ref_time$)
  call apply_ramper (ele, ltt_com_global%lat%ele(ie), err_flag)
enddo

! Adjust particle reference energy if needed.

if (start_orb%p0c == ele%value(p0c_start$)) return

r = start_orb%p0c / ele%value(p0c_start$)
start_orb%vec(2) = r * start_orb%vec(2)
start_orb%vec(4) = r * start_orb%vec(4)
start_orb%vec(6) = r * start_orb%vec(6) + (start_orb%p0c - ele%value(p0c_start$)) / ele%value(p0c_start$)
start_orb%p0c = ele%value(p0c_start$)

end subroutine
