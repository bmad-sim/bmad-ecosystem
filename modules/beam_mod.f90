module beam_mod

use csr_mod

interface assignment (=)
  module procedure bunch_equal_bunch
  module procedure beam_equal_beam
end interface

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track_beam (lat, beam, ix1, ix2)
!
! Subroutine to track a beam of particles from the end of
! lat%ele(ix1) Through to the end of lat%ele(ix2).
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
!-

subroutine track_beam (lat, beam, ix1, ix2)

implicit none

type (lat_struct) :: lat
type (beam_struct) :: beam

integer, optional, intent(in) :: ix1, ix2
integer i, i1, i2, j

! Init

i1 = 0
if (present(ix1)) i1 = ix1
i2 = lat%n_ele_track
if (present(ix2)) i2 = ix2

! Loop over all elements in the lattice

do i = i1+1, i2
  call track1_beam_lat (beam, lat, i, beam)
enddo

end subroutine track_beam

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_beam_lat (beam_start, lat, ix_ele, beam_end)
!
! Subroutine to track a beam of particles through a single element.
!
! Note: This routine is overloaded by the routine track1_beam. See this
! routine for more details.
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
!-

subroutine track1_beam_lat (beam_start, lat, ix_ele, beam_end)

implicit none

type (beam_struct) beam_start
type (beam_struct) :: beam_end
type (lat_struct) :: lat

integer i, ix_ele, n_mode

! zero the long-range wakes if they exist.

if (associated(lat%ele(ix_ele)%wake)) then
  lat%ele(ix_ele)%wake%lr%norm_sin = 0; lat%ele(ix_ele)%wake%lr%norm_cos = 0
  lat%ele(ix_ele)%wake%lr%skew_sin = 0; lat%ele(ix_ele)%wake%lr%skew_cos = 0
endif

! loop over all bunches in a beam

do i = 1, size(beam_start%bunch)
  call track1_bunch_lat (beam_start%bunch(i), lat, ix_ele, beam_end%bunch(i))
enddo

end subroutine track1_beam_lat

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_beam_ele (beam_start, ele, param, beam_end)
!
! Subroutine to track a beam of particles through a single element.
!
! Note: This routine is overloaded by the routine track1_beam. See this
! routine for more details.
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

subroutine track1_beam_ele (beam_start, ele, param, beam_end)

implicit none

type (beam_struct) beam_start
type (beam_struct), target :: beam_end
type (ele_struct) ele
type (lat_param_struct) param

integer i, n_mode

! zero the long-range wakes if they exist.

if (associated(ele%wake)) then
  ele%wake%lr%norm_sin = 0; ele%wake%lr%norm_cos = 0
  ele%wake%lr%skew_sin = 0; ele%wake%lr%skew_cos = 0
endif

! loop over all bunches in a beam

do i = 1, size(beam_start%bunch)
  call track1_bunch_ele (beam_start%bunch(i), ele, param, beam_end%bunch(i))
enddo

end subroutine track1_beam_ele

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_bunch_lat (bunch_start, lat, ix_ele, bunch_end)
!
! Subroutine to track a bunch of particles through an element.
!
! Note: This routine is overloaded by the routine track1_bunch. See this
! routine for more details.
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
!-

subroutine track1_bunch_lat (bunch_start, lat, ix_ele, bunch_end)

implicit none

type (bunch_struct) bunch_start, bunch_end
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (ele_struct), save :: rf_ele

real(rp) charge
integer i, j, n, ix_ele

!------------------------------------------------
! space charge tracking will also include wakes if they are on too.

if (lat%ele(ix_ele)%tracking_method == custom$) then
  call track1_bunch_custom (bunch_start, lat, ix_ele, bunch_end)

elseif (bmad_com%coherent_synch_rad_on .and. lat%ele(ix_ele)%csr_calc_on) then
  call track1_bunch_csr (bunch_start, lat, ix_ele, bunch_end)

else
  call track1_bunch_ele (bunch_start, lat%ele(ix_ele), lat%param, bunch_end)

endif

end subroutine track1_bunch_lat

end module
