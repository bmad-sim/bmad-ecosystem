!+
! Function lord_edge_aligned (slave, slave_edge, lord) result (is_aligned)
!
! Routine to determine if the edge of a super_lord is aligned with a given edge 
! of a super_slave or slice_slave.
!
! Input:
!   slave       -- ele_struct: Slave element.
!   slave_edge  -- integer: End under consideration: entrance_end$, exit_end$, in_between$, etc.
!   lord        -- ele_struct: Lord element.
!
! Output:
!   is_aligned  -- integer: True if a lord edge is aligned with the slave edge.
!                   If slave_edge is not entrance_end$ nor exit_end$ then is_aligned is False.
!- 

function lord_edge_aligned (slave, slave_edge, lord) result (is_aligned)

use equal_mod, dummy => lord_edge_aligned

implicit none

type (ele_struct), target :: slave, lord
type (branch_struct), pointer :: branch
integer slave_edge, ix_slave
real(rp) s_lord
logical is_aligned
character(*), parameter :: r_name = 'lord_edge_aligned'

! 

select case (stream_ele_end(slave_edge, slave%orientation))
case (upstream_end$)
  s_lord = lord%s_start
  branch => slave%branch
  if (associated(branch)) then
    if (s_lord < branch%ele(0)%s) s_lord = s_lord + branch%param%total_length
  endif
  is_aligned = (abs(slave%s_start - s_lord) < bmad_com%significant_length)

case (downstream_end$)
  is_aligned = (abs(slave%s - lord%s) < bmad_com%significant_length)

case default
  is_aligned = .false.
end select

end function lord_edge_aligned
