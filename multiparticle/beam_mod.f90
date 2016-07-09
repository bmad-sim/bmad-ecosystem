module beam_mod

use multipass_mod
use beam_utils

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track_beam (lat, beam, ele1, ele2, err, centroid)
!
! Subroutine to track a beam of particles from the end of
! ele1 Through to the end of ele2. Both must be in the same lattice branch.
!
! Note: zero_lr_wakes_in_lat needs to be called to initial wakes before tracking.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   lat          -- lat_struct: Lattice to track through.
!   beam         -- Beam_struct: Beam at end of element ix1.
!   ele1         -- Ele_struct, optional: Starting element (this element 
!                     is NOT tracked through). Default is lat%ele(0).
!   ele2         -- Ele_struct, optional: Ending element.
!                     Default is lat%ele(lat%n_ele_track).
!   centroid(0:) -- coord_struct, optional: Centroid orbit. Only needed if CSR is on and the
!                     central orbit is far from the zero orbit.
!
! Output:
!   beam   -- beam_struct: Beam at end of element ix2.
!   err    -- Logical: Set true if there is an error. 
!                  EG: Too many particles lost for a CSR calc.
!-

subroutine track_beam (lat, beam, ele1, ele2, err, centroid)

implicit none

type (lat_struct), target :: lat
type (beam_struct) :: beam
type (branch_struct), pointer :: branch
type (ele_struct), optional, target :: ele1, ele2
type (ele_struct), pointer :: e1, e2
type (coord_struct), optional :: centroid(0:)

integer i, j

logical err

! Init

e1 => lat%ele(0)
if (present(ele1)) e1 => ele1
e2 => lat%ele(lat%n_ele_track)
if (present(ele2)) e2 => ele2

! Loop over all elements in the lattice

branch => lat%branch(e1%ix_branch)
do i = e1%ix_ele+1, e2%ix_ele
  call track1_beam (beam, lat, branch%ele(i), beam, err, centroid)
  if (err) return
enddo

end subroutine track_beam

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_beam (beam_start, lat, ele, beam_end, err, centroid)
!
! Routine to track a beam of particles through a single element.
! See also track1_beam_simple.
!
! Note: zero_lr_wakes_in_lat needs to be called to initial wakes before tracking.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   beam_start   -- Beam_struct: Starting beam position.
!   lat          -- lat_struct: Lattice containing element to be tracked through.
!   ele          -- Ele_struct: Element to track through.
!   centroid(0:) -- coord_struct, optional: Centroid orbit. Only needed if CSR is on and the
!                     central orbit is far from the zero orbit.
!
! Output:
!   beam_end    -- beam_struct: Ending beam position.
!   err         -- Logical: Set true if there is an error. 
!                    EG: Too many particles lost for a CSR calc.
!-

subroutine track1_beam (beam_start, lat, ele, beam_end, err, centroid)

implicit none

type (beam_struct) beam_start
type (beam_struct) :: beam_end
type (lat_struct) :: lat
type (ele_struct) ele
type (coord_struct), optional :: centroid(0:)

integer i, n_mode
logical err

! loop over all bunches in a beam

do i = 1, size(beam_start%bunch)
  call track1_bunch (beam_start%bunch(i), lat, ele, beam_end%bunch(i), err, centroid)
  if (err) return
enddo

end subroutine track1_beam

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_beam_simple (beam_start, ele, param, beam_end)
!
! Routine to track a beam of particles through a single element.
! This routine does *not* include multiparticle effects such as
! wakefields or CSR. This routine should only be used when the
! routine track1_beam cannot be used.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   beam_start  -- beam_struct: starting beam position
!   ele         -- Ele_struct: The element to track through.
!   param       -- lat_param_struct: General parameters.
!
! Output:
!   beam_end    -- beam_struct: ending beam position.
!-

subroutine track1_beam_simple (beam_start, ele, param, beam_end)

implicit none

type (beam_struct) beam_start
type (beam_struct), target :: beam_end
type (ele_struct) ele
type (lat_param_struct) param

integer i, j

! loop over all bunches in a beam

do i = 1, size(beam_start%bunch)
  do j = 1, size(beam_start%bunch(i)%particle)
    call track1 (beam_start%bunch(i)%particle(j), ele, param, beam_end%bunch(i)%particle(j))
  enddo
enddo

end subroutine track1_beam_simple

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_bunch (bunch_start, lat, ele, bunch_end, err, centroid)
!
! Subroutine to track a bunch of particles through an element.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   bunch_start  -- bunch_struct: Starting bunch position.
!   lat          -- lat_struct: Lattice containing element to be tracked through.
!   ele          -- Ele_struct: element to track through.
!   centroid(0:) -- coord_struct, optional: Centroid orbit. Only needed if CSR is on and the
!                     central orbit is far from the zero orbit.
!
! Output:
!   bunch_end -- Bunch_struct: Ending bunch position.
!   err       -- Logical: Set true if there is an error. 
!                  EG: Too many particles lost for a CSR calc.
!-

subroutine track1_bunch (bunch_start, lat, ele, bunch_end, err, centroid)

use csr_old_mod, only: track1_bunch_csr_old
use csr_mod, only: track1_bunch_csr
use beam_utils, only: track1_bunch_hom

implicit none

type (bunch_struct) bunch_start, bunch_end
type (lat_struct), target :: lat
type (ele_struct) :: ele
type (ele_struct), pointer :: lord, slave
type (wake_lr_struct), pointer :: lr, lr_chain
type (ele_pointer_struct), allocatable :: chain_ele(:)
type (coord_struct), optional :: centroid(0:)

integer i, j, n, im, ix_pass, ixs, ix, n_links

logical csr_on, err

character(*), parameter :: r_name = 'track1_bunch'

!------------------------------------------------
! Tracking

csr_on = bmad_com%coherent_synch_rad_on .and. ele%csr_calc_on
if (csr_param%ix1_ele_csr > -1) csr_on = csr_on .and. (ele%ix_ele > csr_param%ix1_ele_csr) 
if (csr_param%ix2_ele_csr > -1) csr_on = csr_on .and. (ele%ix_ele <= csr_param%ix2_ele_csr) 

if (csr_on) then
  if (csr_param%use_csr_old) then
    call track1_bunch_csr_old (bunch_start, lat, ele, bunch_end, err)
  else
    if (.not. present(centroid)) then
      call out_io (s_fatal$, r_name, 'BUNCH CENTROID MUST BE SUPPLIED FOR CSR CALCULATION!')
      if (global_com%exit_on_error) call err_exit
      return
    endif
    call track1_bunch_csr (bunch_start, ele, centroid, bunch_end, err)
  endif
  bunch_end%ix_ele = ele%ix_ele

! Non csr / non space-charge tracking
else
  err = .false.
  call track1_bunch_hom (bunch_start, ele, lat%param, bunch_end)
  bunch_end%ix_ele = ele%ix_ele
endif

if (err) return

! If there are wakes...

if (associated(ele%wake)) then

  ! If a super_slave, the lr wake in the lord is the sum of the slaves.

  if (ele%slave_status == super_slave$) then
    do i = 1, ele%n_lord
      lord => pointer_to_lord(ele, i)
      if (lord%lord_status /= super_lord$) cycle
      lord%wake%lr%b_sin = 0;  lord%wake%lr%b_cos = 0
      lord%wake%lr%a_sin = 0;  lord%wake%lr%a_cos = 0
      lord%wake%lr%t_ref = 0
      do j = 1, lord%n_slave
        slave => pointer_to_slave(lord, 1)
        lord%wake%lr%b_sin = lord%wake%lr%b_sin + slave%wake%lr%b_sin
        lord%wake%lr%b_cos = lord%wake%lr%b_cos + slave%wake%lr%b_cos
        lord%wake%lr%a_sin = lord%wake%lr%a_sin + slave%wake%lr%a_sin
        lord%wake%lr%a_cos = lord%wake%lr%a_cos + slave%wake%lr%a_cos
      enddo
    enddo
  endif

  ! If part of multipass: Transfer the lr wake to all the elements in the chain.
  ! A chain is a set of elements in the tracking lattice that all represent 
  ! the same physical element.

  call multipass_chain (ele, ix_pass, n_links, chain_ele)
  do i = 1, n_links
    if (i == ix_pass) cycle
    do j = 1, size(ele%wake%lr)
      chain_ele(i)%ele%wake%lr(j) = ele%wake%lr(j)
    enddo
  enddo

  lord => pointer_to_multipass_lord (ele)
  if (associated(lord)) then 
    do j = 1, size(ele%wake%lr)
      lord%wake%lr(j) = ele%wake%lr(j)
    enddo
  endif

endif

end subroutine track1_bunch


end module
