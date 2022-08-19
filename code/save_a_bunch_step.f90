!+
! Subroutine save_a_bunch_step (bunch_track, ele, bunch)
!
! Routine to save the bunch parameters when tracking a bunch step-by-step.
!     
!
! Input:
!   bunch_track     -- bunch_track_struct: Track up to now. If bunch_track%n_pt < 0, the structure will be reinitialized.
!   ele             -- ele_struct: Element being tracked through.
!   bunch           -- bunch_struct: Bunch whose parameters are to be saved.
!
! Ouput:
!   bunch_track     -- bunch_track_struct: Track with current bunch info appended on.
!-

subroutine save_a_bunch_step (bunch_track, ele, bunch)

use beam_utils, dummy => save_a_bunch_step

implicit none

type (bunch_track_struct), target :: bunch_track, bunch_track2
type (bunch_params_struct), pointer :: bp
type (ele_struct), target :: ele
type (bunch_struct) bunch

integer n_pt, n, n_old
real(rp) s_rel
logical err

! Init

if (.not. allocated (bunch_track%pt)) then
  allocate(bunch_track%pt(0:100))
  bunch_track%n_pt = -1
endif

if (bunch_track%n_pt < 0) then
  bunch_track%n_pt = -1
endif

!

bunch_track%n_pt = bunch_track%n_pt + 1
n_pt = bunch_track%n_pt
n_old = ubound(bunch_track%pt, 1)

if (n_pt > n_old) then
  n = 1.5 * n_pt
  call move_alloc (bunch_track%pt, bunch_track2%pt)
  allocate(bunch_track%pt(0:n))
  bunch_track%pt(:n_old) = bunch_track2%pt
end if

!

bp => bunch_track%pt(n_pt)
call calc_bunch_params (bunch, bp, err)

end subroutine save_a_bunch_step
