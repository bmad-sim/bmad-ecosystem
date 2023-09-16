!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_ele_to_ele_track (ix_universe, ix_branch, ix_ele, ix_ele_track)
!
! Subroutine to compute ix_ele_track:
!   = ix_ele                 if ix_ele <= lat%branch(ix)branch)%n_ele_track
!   = ix_slave_at_exit_end   if ix_ele is a super_lord  
!   = -1                     otherwise
!
! Input:
!   ix_universe -- Integer: Universe index.
!   ix_branch   -- Integer: Branch index.
!   ix_ele      -- Integer: Element index
!
! Output:
!   ix_ele_track -- Integer: Corresponding element in the tracking 
!                         part of the lattice.
!-

subroutine tao_ele_to_ele_track (ix_universe, ix_branch, ix_ele, ix_ele_track)

use tao_interface, dummy => tao_ele_to_ele_track

implicit none

type (lat_struct), pointer :: lat
type (ele_struct), pointer :: slave

integer ix_universe, ix_branch, ix_ele, ix_ele_track
integer i_uni, ix_c

!

i_uni = tao_universe_index(ix_universe)
lat => s%u(i_uni)%model%lat

if (ix_ele < 0) then
  ix_ele_track = -1

elseif (ix_ele <= lat%branch(ix_branch)%n_ele_track) then
  ix_ele_track = ix_ele

elseif (lat%ele(ix_ele)%lord_status == super_lord$) then
  slave => pointer_to_slave (lat%ele(ix_ele), lat%ele(ix_ele)%n_slave)
  ix_ele_track = slave%ix_ele ! element at exit end.

else  ! overlays, multipass_lords, etc.
  ix_ele_track = -1
endif

end subroutine tao_ele_to_ele_track
