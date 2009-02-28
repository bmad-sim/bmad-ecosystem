module beam_mod

use csr_mod
use multipass_mod

interface assignment (=)
  module procedure bunch_equal_bunch
  module procedure beam_equal_beam
end interface

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track_beam (lat, beam, ix1, ix2, err)
!
! Subroutine to track a beam of particles from the end of
! lat%ele(ix1) Through to the end of lat%ele(ix2).
!
! Note: zero_lr_wakes_in_lat needs to be called to initial wakes before tracking.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   lat    -- lat_struct: Lattice to track through.
!   beam   -- Beam_struct: Beam at end of element ix1.
!   ix1    -- Integer, optional: Index of starting element (this element 
!               is NOT tracked through). Default is 0.
!   ix2    -- Integer, optional: Index of ending element.
!               Default is lat%n_ele_track.
!
! Output:
!   beam   -- beam_struct: Beam at end of element ix2.
!   err       -- Logical: Set true if there is an error. 
!                  EG: Too many particles lost for a CSR calc.
!-

subroutine track_beam (lat, beam, ix1, ix2, err)

implicit none

type (lat_struct) :: lat
type (beam_struct) :: beam

integer, optional, intent(in) :: ix1, ix2
integer i, i1, i2, j

logical err

! Init

i1 = 0
if (present(ix1)) i1 = ix1
i2 = lat%n_ele_track
if (present(ix2)) i2 = ix2

! Loop over all elements in the lattice

do i = i1+1, i2
  call track1_beam (beam, lat, i, beam, err)
  if (err) return
enddo

end subroutine track_beam

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_beam (beam_start, lat, ix_ele, beam_end, err)
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
!   beam_start  -- Beam_struct: Starting beam position.
!   lat         -- lat_struct: Lattice containing element to be tracked through.
!   ix_ele      -- Integer: Index of element to track through.
!
! Output:
!   beam_end    -- beam_struct: Ending beam position.
!   err         -- Logical: Set true if there is an error. 
!                    EG: Too many particles lost for a CSR calc.
!-

subroutine track1_beam (beam_start, lat, ix_ele, beam_end, err)

implicit none

type (beam_struct) beam_start
type (beam_struct) :: beam_end
type (lat_struct) :: lat

integer i, ix_ele, n_mode
logical err

! loop over all bunches in a beam

do i = 1, size(beam_start%bunch)
  call track1_bunch (beam_start%bunch(i), lat, ix_ele, beam_end%bunch(i), err)
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
    call track1_particle (beam_start%bunch(i)%particle(j), ele, param, beam_end%bunch(i)%particle(j))
  enddo
enddo

end subroutine track1_beam_simple

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_bunch (bunch_start, lat, ix_ele, bunch_end, err)
!
! Subroutine to track a bunch of particles through an element.
!
! Each particle experiences a different longitudinal short-range wakefield.
! bmad_com%grad_loss_sr_wake is used to tell track1_bmad the appropriate loss
! for each particle.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   bunch_start -- bunch_struct: Starting bunch position.
!   lat         -- lat_struct: Lattice containing element to be tracked through.
!   ix_ele      -- Integer: Index of element to track through.
!
! Output:
!   bunch_end -- Bunch_struct: Ending bunch position.
!   err       -- Logical: Set true if there is an error. 
!                  EG: Too many particles lost for a CSR calc.
!-

subroutine track1_bunch (bunch_start, lat, ix_ele, bunch_end, err)

implicit none

type (bunch_struct) bunch_start, bunch_end
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, lord, chain_ele
type (ele_struct), save :: rf_ele
type (lr_wake_struct), pointer :: lr

real(rp) charge, dt, c, s, k

integer i, j, n, ix_ele, ix_pass, ixs, ix
integer, save, allocatable :: ix_chain(:)

logical csr_on, err

!------------------------------------------------
! space charge tracking will also include wakes if they are on too.

bunch_end%ix_ele = ix_ele

if (lat%ele(ix_ele)%tracking_method == custom$) then
  call track1_bunch_custom (bunch_start, lat, ix_ele, bunch_end)
  return
endif

csr_on = bmad_com%coherent_synch_rad_on .and. lat%ele(ix_ele)%csr_calc_on
if (csr_param%ix1_ele_csr > -1) csr_on = csr_on .and. (ix_ele > csr_param%ix1_ele_csr) 
if (csr_param%ix2_ele_csr > -1) csr_on = csr_on .and. (ix_ele <= csr_param%ix2_ele_csr) 

if (csr_on) then
  call track1_bunch_csr (bunch_start, lat, ix_ele, bunch_end, err)
  return
endif

! Non-csr tracking

err = .false.
ele => lat%ele(ix_ele)
call track1_bunch_hom (bunch_start, ele, lat%param, bunch_end)

! If there are wakes...
! If a super_slave, the lr wake in the lord is the sum of the slaves.

if (associated(ele%wake)) then
  if (ele%slave_status == super_slave$) then
    do i = ele%ic1_lord, ele%ic2_lord
      ixs = lat%control(lat%ic(i))%ix_lord
      lord => lat%ele(ixs)
      lord%wake%lr%norm_sin = 0;  lord%wake%lr%norm_cos = 0
      lord%wake%lr%skew_sin = 0;  lord%wake%lr%skew_cos = 0
      do j = lord%ix1_slave, lord%ix2_slave
        ix = lat%control(j)%ix_slave
        lord%wake%lr%norm_sin = lord%wake%lr%norm_sin + lat%ele(ix)%wake%lr%norm_sin
        lord%wake%lr%norm_cos = lord%wake%lr%norm_cos + lat%ele(ix)%wake%lr%norm_cos
        lord%wake%lr%skew_sin = lord%wake%lr%skew_sin + lat%ele(ix)%wake%lr%skew_sin
        lord%wake%lr%skew_cos = lord%wake%lr%skew_cos + lat%ele(ix)%wake%lr%skew_cos
      enddo
    enddo
  endif

  ! If part of multipass: Transfer the lr wake to all the elements in the chain.
  ! A chain is a set of elements in the tracking lattice that all represent 
  ! the same physical element.

  call multipass_chain (ix_ele, lat, ix_pass, ix_chain)
  if (ix_pass > 0) then
    do i = 1, size(ix_chain)
      if (i == ix_pass) cycle
      lr => ele%wake%lr(i)
      dt = chain_ele%ref_time - ele%ref_time
      k = twopi * lr%freq
      c = cos (-dt * k)
      s = sin (-dt * k)

      chain_ele => lat%ele(ix_chain(i))
      chain_ele%wake%lr%norm_sin =  c * lr%norm_sin + s * lr%norm_cos
      chain_ele%wake%lr%norm_cos = -s * lr%norm_sin + c * lr%norm_cos
      chain_ele%wake%lr%skew_sin =  c * lr%skew_sin + s * lr%skew_cos
      chain_ele%wake%lr%skew_cos = -s * lr%skew_sin + c * lr%skew_cos
      chain_ele%wake%lr%t_ref    = lr%t_ref + dt
    enddo
  endif

endif

end subroutine track1_bunch


end module
