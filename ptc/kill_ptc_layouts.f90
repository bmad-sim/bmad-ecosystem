!+
! Subroutine kill_ptc_layouts (lat)
!
! Routine to kill the layouts associated with a Bmad lattice.
!
! Input: 
!   lat  -- lat_struct: Bmad lattice with associated layouts.
!-

subroutine kill_ptc_layouts (lat)

use ptc_layout_mod, except_dummy => kill_ptc_layouts
use madx_ptc_module

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch

integer ib, il

!

if (.not. allocated(lat%branch)) return

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  if (associated(branch%ptc%m_t_layout)) then
    call kill_layout_in_universe(branch%ptc%m_t_layout)
    deallocate(branch%ptc%m_t_layout)
  endif

  if (allocated(branch%ptc%m_u_layout)) then
    do il = 1, size(branch%ptc%m_u_layout)
      call kill_layout_in_universe(branch%ptc%m_u_layout(il)%ptr)
    enddo
    deallocate(branch%ptc%m_u_layout)
  endif
enddo

end subroutine kill_ptc_layouts

