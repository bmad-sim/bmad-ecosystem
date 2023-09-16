!+
! Subroutine save_a_beam_step (ele, beam, bunch_tracks, s_body, is_time_coords)
!
! Routine to save the bunch parameters when tracking a beam step-by-step.
!     
!
! Input:
!   ele             -- ele_struct: Element being tracked through.
!   beam            -- beam_struct: Bunches in the beam whose parameters are to be saved.
!   bunch_tracks(:) -- bunch_track_struct, optional: Track up to now. If bunch_tracks%n_pt < 0, the structure will be reinitialized.
!   s_body          -- real(rp), optional: Body s-position from beginning of element.
!   is_time_coords  -- logical, optional: Default is False. If True, input beam is using time coordinates in which
!                       case there will be a conversion to s-coords before bunch_params are computed.
!
! Ouput:
!   bunch_tracks(:) -- bunch_track_struct, optional: Track with current bunch info appended on. This routine does nothing
!                       if this argument is not present.
!-

subroutine save_a_beam_step (ele, beam, bunch_tracks, s_body, is_time_coords)

use beam_utils, dummy => save_a_beam_step

implicit none

type (ele_struct), target :: ele
type (beam_struct) beam
type (bunch_track_struct), optional, target :: bunch_tracks(:)

integer i
real(rp), optional :: s_body
logical, optional :: is_time_coords

!

if (.not. present(bunch_tracks)) return
do i = 1, size(beam%bunch)
  call save_a_bunch_step(ele, beam%bunch(i), bunch_tracks(i), s_body, is_time_coords)
enddo

end subroutine save_a_beam_step
