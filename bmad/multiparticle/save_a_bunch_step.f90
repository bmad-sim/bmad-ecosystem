!+
! Subroutine save_a_bunch_step (ele, bunch, bunch_track, s_body, is_time_coords)
!
! Routine to save the bunch parameters when tracking a bunch step-by-step.
!     
!
! Input:
!   ele             -- ele_struct: Element being tracked through.
!   bunch           -- bunch_struct: Bunch whose parameters are to be saved.
!   bunch_track     -- bunch_track_struct, optional: Track up to now. If bunch_track%n_pt < 0, the structure will be reinitialized.
!   s_body          -- real(rp), optional: Body s-position from beginning of element.
!   is_time_coords  -- logical, optional: Default is False. If True, input bunch is using time coordinates in which
!                       case there will be a conversion to s-coords before bunch_params are computed.
!
! Ouput:
!   bunch_track     -- bunch_track_struct, optional: Track with current bunch info appended on. This routine does nothing
!                       if this argument is not present.
!-

subroutine save_a_bunch_step (ele, bunch, bunch_track, s_body, is_time_coords)

use beam_utils, dummy => save_a_bunch_step

implicit none

type (ele_struct), target :: ele
type (bunch_struct) bunch
type (bunch_track_struct), optional, target :: bunch_track
type (bunch_track_struct), target :: bunch_track2

integer n_pt, n, n_old
real(rp), optional :: s_body
real(rp) s_last, s_here
logical, optional :: is_time_coords
logical err

! Init

if (.not. present(bunch_track)) return
if (bunch_track%ds_save < 0) return

if (.not. allocated (bunch_track%pt)) then
  allocate(bunch_track%pt(0:100))
  bunch_track%n_pt = -1
endif

if (bunch_track%n_pt < 0) then
  bunch_track%n_pt = -1
  s_last = -1e30  ! Something large and negative
else
  s_last = bunch_track%pt(bunch_track%n_pt)%s
endif

if (present(s_body)) then
  s_here = s_body + ele%s_start
else
  n = count(bunch%particle(:)%state==alive$)
  if (n == 0) return
  s_here = sum(bunch%particle(:)%s, bunch%particle(:)%state==alive$)/n
endif

if (s_here < s_last + bunch_track%ds_save) return

! Save a track point

bunch_track%n_pt = bunch_track%n_pt + 1
n_pt = bunch_track%n_pt
n_old = ubound(bunch_track%pt, 1)

if (n_pt > n_old) then
  n = 1.5 * n_pt
  call move_alloc (bunch_track%pt, bunch_track2%pt)
  allocate(bunch_track%pt(0:n))
  bunch_track%pt(:n_old) = bunch_track2%pt
end if

call calc_bunch_params (bunch, bunch_track%pt(n_pt), err, is_time_coords = is_time_coords, ele = ele)

end subroutine save_a_bunch_step
