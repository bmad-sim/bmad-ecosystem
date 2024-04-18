module beam_mod

use beam_utils

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track_beam (lat, beam, ele1, ele2, err, centroid, direction, bunch_tracks)
!
! Subroutine to track a beam of particles from the end of
! ele1 Through to the end of ele2. Both must be in the same lattice branch.
!
! Note: To zero wakes between runs, zero_lr_wakes_in_lat needs to be called.
!
! Input:
!   lat             -- lat_struct: Lattice to track through.
!   beam            -- beam_struct: Beam at end of element ix1.
!   ele1            -- ele_struct, optional: Starting element (this element 
!                        is NOT tracked through). Default is lat%ele(0).
!   ele2            -- ele_struct, optional: Ending element.
!                        Default is lat%ele(lat%n_ele_track).
!   centroid(0:)    -- coord_struct, optional: Approximate centroid orbit. Only needed if CSR is on.
!                        Hint: Calculate this before beam tracking by tracking a single particle.
!   direction       -- integer, optional: +1 (default) -> Track forward, -1 -> Track backwards.
!   bunch_tracks(:) -- bunch_track_struct, optional: Existing tracks. If bunch_track%n_pt = -1 then
!                        Overwrite any existing track.
!
! Output:
!   beam            -- beam_struct: Beam at end of element ix2.
!   err             -- logical: Set true if there is an error. 
!                        EG: Too many particles lost for a CSR calc.
!   bunch_tracks(:) -- bunch_track_struct, optional: track information if the tracking method does 
!                        tracking step-by-step. When tracking through multiple elements, the 
!                        trajectory in an element is appended to the existing trajectory. 
!-

subroutine track_beam (lat, beam, ele1, ele2, err, centroid, direction, bunch_tracks)

implicit none

type (lat_struct), target :: lat
type (beam_struct) :: beam
type (ele_struct), optional, target :: ele1, ele2
type (coord_struct), optional :: centroid(0:)
type (bunch_track_struct), optional :: bunch_tracks(:)

integer, optional :: direction
integer i

logical err

!

do i = 1, size(beam%bunch)
  if (present(bunch_tracks)) then
    call track_bunch(lat, beam%bunch(i), ele1, ele2, err, centroid, direction, bunch_tracks(i))
  else
    call track_bunch(lat, beam%bunch(i), ele1, ele2, err, centroid, direction)
  endif
  if (err) return
enddo

end subroutine track_beam

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track_bunch (lat, bunch, ele1, ele2, err, centroid, direction, bunch_track)
!
! Subroutine to track a particle bunch from the end of ele1 Through to the end of ele2.
! Both must be in the same lattice branch.
! With forward tracking, if ele2 is at or before ele1, the tracking will "wrap" around 
! the ends of the lattice.
!
! Note: To zero wakes between runs, zero_lr_wakes_in_lat needs to be called.
!
! Input:
!   lat           -- lat_struct: Lattice to track through.
!   bunch         -- Bunch_struct: Bunch at end of element ix1.
!   ele1          -- Ele_struct, optional: Starting element (this element 
!                      is NOT tracked through). Default is lat%ele(0).
!   ele2          -- Ele_struct, optional: Ending element.
!                      Default is lat%ele(lat%n_ele_track).
!   centroid(0:)  -- coord_struct, optional: Approximate centroid orbit. Only needed if CSR is on.
!                      Hint: Calculate this before bunch tracking by tracking a single particle.
!   direction     -- integer, optional: +1 (default) -> Track forward, -1 -> Track backwards.
!   bunch_track   -- bunch_track_struct, optional: Existing tracks. If bunch_track%n_pt = -1 then
!                        Overwrite any existing track.
!
! Output:
!   bunch       -- bunch_struct: Bunch at end of element ix2.
!   err         -- Logical: Set true if there is an error. 
!                    EG: Too many particles lost for a CSR calc.
!   bunch_track -- bunch_track_struct, optional: track information if the tracking method does 
!                        tracking step-by-step. When tracking through multiple elements, the 
!                        trajectory in an element is appended to the existing trajectory. 
!-

subroutine track_bunch (lat, bunch, ele1, ele2, err, centroid, direction, bunch_track)

implicit none

type (lat_struct), target :: lat
type (bunch_struct) :: bunch
type (branch_struct), pointer :: branch
type (ele_struct), optional, target :: ele1, ele2
type (ele_struct), pointer :: e1, e2
type (coord_struct), optional :: centroid(0:)
type (bunch_track_struct), optional :: bunch_track

integer, optional :: direction
integer i, j

logical err

! Init

e1 => lat%ele(0)
if (present(ele1)) e1 => ele1
e2 => lat%ele(lat%n_ele_track)
if (present(ele2)) e2 => ele2

branch => lat%branch(e1%ix_branch)

! 

if (integer_option(1, direction) == -1) then
  if (e1%ix_ele > e2%ix_ele) then
    do i = e1%ix_ele, e2%ix_ele+1, -1
      call track1_bunch (bunch, branch%ele(i), err, centroid, direction, bunch_track)
      if (err) return
    enddo

  else
    do i = e1%ix_ele, 1, -1
      call track1_bunch (bunch, branch%ele(i), err, centroid, direction, bunch_track)
      if (err) return
    enddo
    do i = branch%n_ele_track, e2%ix_ele+1, -1
      call track1_bunch (bunch, branch%ele(i), err, centroid, direction, bunch_track)
      if (err) return
    enddo
  endif

!

else
  if (e1%ix_ele < e2%ix_ele) then
    do i = e1%ix_ele+1, e2%ix_ele
      call track1_bunch (bunch, branch%ele(i), err, centroid, direction, bunch_track)
      if (err) return
    enddo

  else
    do i = e1%ix_ele+1, branch%n_ele_track
      call track1_bunch (bunch, branch%ele(i), err, centroid, direction, bunch_track)
      if (err) return
    enddo
    do i = 1, e2%ix_ele
      call track1_bunch (bunch, branch%ele(i), err, centroid, direction, bunch_track)
      if (err) return
    enddo
  endif
endif

end subroutine track_bunch

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_bunch (bunch, ele, err, centroid, direction, bunch_track)
!
! Subroutine to track a bunch of particles through an element.
!
! Input:
!   bunch         -- bunch_struct: Starting bunch position.
!   ele           -- Ele_struct: element to track through.
!   centroid(0:)  -- coord_struct, optional: Approximate centroid orbit. Only needed if CSR is on.
!                      Hint: Calculate this before beam tracking by tracking a single particle.
!   direction     -- integer, optional: +1 (default) -> Track forward, -1 -> Track backwards.
!   bunch_track   -- bunch_track_struct, optional: Existing tracks. If bunch_track%n_pt = -1 then
!                        Overwrite any existing track.
!
! Output:
!   bunch       -- Bunch_struct: Ending bunch position.
!   err         -- Logical: Set true if there is an error. 
!                    EG: Too many particles lost for a CSR calc.
!   bunch_track -- bunch_track_struct, optional: Track information appended to track.
!-

subroutine track1_bunch (bunch, ele, err, centroid, direction, bunch_track)

use csr_and_space_charge_mod, only: track1_bunch_csr, track1_bunch_csr3d
use beam_utils, only: track1_bunch_hom
use space_charge_mod, only: track_bunch_to_t, track_bunch_to_s

implicit none

type (bunch_struct) bunch
type (ele_struct) :: ele
type (ele_struct), pointer :: lord, slave, wake_ele
type (wake_lr_mode_struct), pointer :: lr, lr_chain
type (ele_pointer_struct), allocatable :: chain_ele(:)
type (coord_struct), optional :: centroid(0:)
type (bunch_track_struct), optional :: bunch_track

integer, optional :: direction
integer i, j, n, im, ix_pass, ixs, ix, n_links

logical csr_sc_on, err, finished, sc_fft_on, time_rk_tracking, track1_bunch_space_charge_called

character(*), parameter :: r_name = 'track1_bunch'

! Custom tracking

if (associated(track1_bunch_hook_ptr)) then
  call track1_bunch_hook_ptr (bunch, ele, err, centroid, direction, finished, bunch_track)
  if (finished) return
endif

!

if (integer_option(1, direction) == -1 .and. bmad_com%csr_and_space_charge_on) then
  call out_io (s_fatal$, r_name, 'BACKWARDS BUNCH TRACKING WITH CSR/SPACE_CHARGE NOT ALLOWED.')
  if (global_com%exit_on_error) call err_exit
  return
endif

!------------------------------------------------
! Tracking

track1_bunch_space_charge_called = .false.
csr_sc_on = (bmad_com%csr_and_space_charge_on .and. (ele%csr_method /= off$ .or. ele%space_charge_method /= off$))
sc_fft_on = (ele%space_charge_method == cathode_fft_3d$ .or. ele%space_charge_method == fft_3d$)
time_rk_tracking = (ele%tracking_method == time_runge_kutta$ .or. ele%tracking_method == fixed_step_time_runge_kutta$)

if (csr_sc_on .and. ele%space_charge_method == cathode_fft_3d$ .and. ele%csr_method /= off$) then
  call out_io (s_error$, r_name, 'WITH SPACE_CHARGE_METHOD SET TO CATHODE_FFT_3D, CSR EFFECTS CANNOT BE HANDLED SO', &
                                 'CSR_METHOD NEEDS TO BE SET TO OFF. FOR LATTICE ELEMENT: ' // ele%name, &
                                 'ALL PARTICLES OF THE BUNCH WILL BE MARKED AS LOST.')
  goto 9000  ! Mark all particles as lost and return
endif

if (csr_sc_on .and. ele%csr_method == off$ .and. sc_fft_on) then
  if (ele%tracking_method /= time_runge_kutta$ .and. ele%tracking_method /= fixed_step_time_runge_kutta$) then
    call out_io (s_error$, r_name, 'WITH SPACE_CHARGE_METHOD SET TO CATHODE_FFT_3D, THE TRACKING_METHOD SHOULD BE SET TO', &
                                   'TIME_RUNGE_KUTTA OR FIXED_STEP_TIME_RUNGE_KUTTA. FOR LATTICE ELEMENT: ' // ele%name, &
                                   'ALL PARTICLES OF THE BUNCH WILL BE MARKED AS LOST.')
    goto 9000  ! Mark all particles as lost and return
  endif
endif


if (ele%csr_method /= off$ .and. time_rk_tracking) then
  call out_io (s_error$, r_name, 'CSR_METHOD IS NOT OFF FOR LATTICE ELEMENT: ' // ele%name, &
                'THIS IS INCOMPATIBLE WITH TRACKING_METHOD SET TO TIME_RUNGE_KUTTA OR FIXED_STEP_TIME_RUNGE_KUTTA.', &
                'ALL PARTICLES OF THE BUNCH WILL BE MARKED AS LOST.')
  goto 9000  ! Mark all particles as lost and return
endif

! 

if (csr_sc_on .and. ele%key /= match$) then
  if (ele%csr_method == off$ .and. sc_fft_on .and. time_rk_tracking) then 
    call track1_bunch_space_charge (bunch, ele, err, bunch_track = bunch_track)
    track1_bunch_space_charge_called = .true.

  elseif (ele%csr_method == steady_state_3d$) then
    if (bunch%drift_between_t_and_s) call correct_s_t_tracking_conversion(bunch, ele)
    call track1_bunch_csr3d(bunch, ele, centroid, err, bunch_track = bunch_track)
     
  else
    if (.not. present(centroid)) then
      call out_io (s_fatal$, r_name, 'BUNCH CENTROID MUST BE SUPPLIED FOR CSR CALCULATION!')
      if (global_com%exit_on_error) call err_exit
      return
    endif
    if (bunch%drift_between_t_and_s) call correct_s_t_tracking_conversion(bunch, ele)
    call track1_bunch_csr (bunch, ele, centroid, err, bunch_track = bunch_track)

  endif
  bunch%ix_ele = ele%ix_ele

! Non csr / non space-charge tracking
else
  err = .false.
  if (bunch%drift_between_t_and_s) call correct_s_t_tracking_conversion(bunch, ele)
  call track1_bunch_hom (bunch, ele, direction, bunch_track = bunch_track)
  bunch%ix_ele = ele%ix_ele
endif

if (err) return

! Set bunch%t0 to real_garbage if there has been significant tracking outside of track1_bunch_space_charge

if (.not. track1_bunch_space_charge_called .and. ele%value(l$) /= 0) bunch%t0 = real_garbage$

! If there are wakes...

wake_ele => pointer_to_wake_ele(ele)
if (associated(wake_ele)) then

  ! If part of multipass: Transfer the lr wake to all the elements in the chain.
  ! A chain is a set of elements in the tracking lattice that all represent 
  ! the same physical element.

  call multipass_chain (wake_ele, ix_pass, n_links, chain_ele)
  do i = 1, n_links
    if (i == ix_pass) cycle
    chain_ele(i)%ele%wake%lr = wake_ele%wake%lr
  enddo

  lord => pointer_to_multipass_lord (wake_ele)
  if (associated(lord)) then
    lord%wake%lr = wake_ele%wake%lr
  endif

endif

return

!-----------------------------------------------------------------------------------
! Mark all particles as lost and return

9000 continue
bunch%particle%state = lost$
err = .true.
return

!-----------------------------------------------------------------------------------
contains

subroutine correct_s_t_tracking_conversion(bunch, ele)

type (bunch_struct) bunch
type (ele_struct) ele

!

call track_bunch_to_t(bunch, bunch%t0, ele%branch)
bunch%drift_between_t_and_s = .false.
call track_bunch_to_s(bunch, ele%s_start, ele%branch)

end subroutine correct_s_t_tracking_conversion

end subroutine track1_bunch

end module
