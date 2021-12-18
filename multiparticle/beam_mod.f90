module beam_mod

use beam_utils

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track_beam (lat, beam, ele1, ele2, err, centroid, direction)
!
! Subroutine to track a beam of particles from the end of
! ele1 Through to the end of ele2. Both must be in the same lattice branch.
!
! Note: To zero wakes between runs, zero_lr_wakes_in_lat needs to be called.
!
! Input:
!   lat          -- lat_struct: Lattice to track through.
!   beam         -- Beam_struct: Beam at end of element ix1.
!   ele1         -- Ele_struct, optional: Starting element (this element 
!                     is NOT tracked through). Default is lat%ele(0).
!   ele2         -- Ele_struct, optional: Ending element.
!                     Default is lat%ele(lat%n_ele_track).
!   centroid(0:) -- coord_struct, optional: Approximate centroid orbit. Only needed if CSR is on.
!                     Hint: Calculate this before beam tracking by tracking a single particle.
!   direction    -- integer, optional: +1 (default) -> Track forward, -1 -> Track backwards.
!
! Output:
!   beam   -- beam_struct: Beam at end of element ix2.
!   err    -- Logical: Set true if there is an error. 
!                  EG: Too many particles lost for a CSR calc.
!-

subroutine track_beam (lat, beam, ele1, ele2, err, centroid, direction)

implicit none

type (lat_struct), target :: lat
type (beam_struct) :: beam
type (ele_struct), optional, target :: ele1, ele2
type (coord_struct), optional :: centroid(0:)

integer, optional :: direction
integer i

logical err

!

do i = 1, size(beam%bunch)
  call track_bunch(lat, beam%bunch(i), ele1, ele2, err, centroid, direction)
  if (err) return
enddo

end subroutine track_beam

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track_bunch (lat, bunch, ele1, ele2, err, centroid, direction)
!
! Subroutine to track a particle bunch from the end of ele1 Through to the end of ele2.
! Both must be in the same lattice branch.
! With forward tracking, if ele2 is at or before ele1, the tracking will "wrap" around 
! the ends of the lattice.
!
! Note: To zero wakes between runs, zero_lr_wakes_in_lat needs to be called.
!
! Input:
!   lat          -- lat_struct: Lattice to track through.
!   bunch        -- Bunch_struct: Bunch at end of element ix1.
!   ele1         -- Ele_struct, optional: Starting element (this element 
!                     is NOT tracked through). Default is lat%ele(0).
!   ele2         -- Ele_struct, optional: Ending element.
!                     Default is lat%ele(lat%n_ele_track).
!   centroid(0:) -- coord_struct, optional: Approximate centroid orbit. Only needed if CSR is on.
!                     Hint: Calculate this before bunch tracking by tracking a single particle.
!   direction    -- integer, optional: +1 (default) -> Track forward, -1 -> Track backwards.
!
! Output:
!   bunch     -- bunch_struct: Bunch at end of element ix2.
!   err       -- Logical: Set true if there is an error. 
!                  EG: Too many particles lost for a CSR calc.
!-

subroutine track_bunch (lat, bunch, ele1, ele2, err, centroid, direction)

implicit none

type (lat_struct), target :: lat
type (bunch_struct) :: bunch
type (branch_struct), pointer :: branch
type (ele_struct), optional, target :: ele1, ele2
type (ele_struct), pointer :: e1, e2
type (coord_struct), optional :: centroid(0:)

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
      call track1_bunch (bunch, branch%ele(i), err, centroid, direction)
      if (err) return
    enddo

  else
    do i = e1%ix_ele, 1, -1
      call track1_bunch (bunch, branch%ele(i), err, centroid, direction)
      if (err) return
    enddo
    do i = branch%n_ele_track, e2%ix_ele+1, -1
      call track1_bunch (bunch, branch%ele(i), err, centroid, direction)
      if (err) return
    enddo
  endif

!

else
  if (e1%ix_ele < e2%ix_ele) then
    do i = e1%ix_ele+1, e2%ix_ele
      call track1_bunch (bunch, branch%ele(i), err, centroid, direction)
      if (err) return
    enddo

  else
    do i = e1%ix_ele+1, branch%n_ele_track
      call track1_bunch (bunch, branch%ele(i), err, centroid, direction)
      if (err) return
    enddo
    do i = 1, e2%ix_ele
      call track1_bunch (bunch, branch%ele(i), err, centroid, direction)
      if (err) return
    enddo
  endif
endif

end subroutine track_bunch

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_bunch (bunch, ele, err, centroid, direction)
!
! Subroutine to track a bunch of particles through an element.
!
! Input:
!   bunch        -- bunch_struct: Starting bunch position.
!   ele          -- Ele_struct: element to track through.
!   centroid(0:) -- coord_struct, optional: Approximate centroid orbit. Only needed if CSR is on.
!                     Hint: Calculate this before beam tracking by tracking a single particle.
!   direction    -- integer, optional: +1 (default) -> Track forward, -1 -> Track backwards.
!
! Output:
!   bunch     -- Bunch_struct: Ending bunch position.
!   err       -- Logical: Set true if there is an error. 
!                  EG: Too many particles lost for a CSR calc.
!-

subroutine track1_bunch (bunch, ele, err, centroid, direction)

use csr_old_mod, only: track1_bunch_csr_old
use csr_and_space_charge_mod, only: track1_bunch_csr, track1_bunch_csr3d
use beam_utils, only: track1_bunch_hom

implicit none

type (bunch_struct) bunch
type (ele_struct) :: ele
type (ele_struct), pointer :: lord, slave, wake_ele
type (wake_lr_mode_struct), pointer :: lr, lr_chain
type (ele_pointer_struct), allocatable :: chain_ele(:)
type (coord_struct), optional :: centroid(0:)

integer, optional :: direction
integer i, j, n, im, ix_pass, ixs, ix, n_links

logical csr_sc_on, err, finished

character(*), parameter :: r_name = 'track1_bunch'

! Custom tracking

call track1_bunch_hook (bunch, ele, err, centroid, direction, finished)
if (finished) return

!

if (integer_option(1, direction) == -1 .and. bmad_com%csr_and_space_charge_on) then
  call out_io (s_fatal$, r_name, 'BACKWARDS BUNCH TRACKING WITH CSR/SPACE_CHARGE NOT ALLOWED.')
  if (global_com%exit_on_error) call err_exit
  return
endif


!------------------------------------------------
! Tracking

csr_sc_on = bmad_com%csr_and_space_charge_on .and. (ele%csr_method /= off$ .or. ele%space_charge_method /= off$)

if (csr_sc_on .and. ele%key /= match$) then
  if (ele%key == e_gun$ .and. ele%value(l_cathode_region$) /= 0) then
    call track1_bunch_e_gun_space_charge (bunch, ele, err)
    
  elseif (ele%csr_method == steady_state_3d$) then
     call track1_bunch_csr3d(bunch, ele, centroid, err)
     
  else
    if (.not. present(centroid)) then
      call out_io (s_fatal$, r_name, 'BUNCH CENTROID MUST BE SUPPLIED FOR CSR CALCULATION!')
      if (global_com%exit_on_error) call err_exit
      return
    endif
    call track1_bunch_csr (bunch, ele, centroid, err)

  endif
  bunch%ix_ele = ele%ix_ele

! Non csr / non space-charge tracking
else
  err = .false.
  call track1_bunch_hom (bunch, ele, direction)
  bunch%ix_ele = ele%ix_ele
endif

if (err) return

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

end subroutine track1_bunch

end module
