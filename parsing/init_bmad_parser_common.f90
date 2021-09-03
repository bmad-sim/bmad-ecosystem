!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine init_bmad_parser_common()

use bmad_parser_struct

implicit none

integer nn, i

! 

if (allocated(bp_com2%const)) deallocate (bp_com2%const)

nn = size(physical_const_list)  ! number of standard (non-user defined) constants
allocate (bp_com2%const(nn))

do i = 1, nn
  bp_com2%const(i) = bp_const_struct(upcase(physical_const_list(i)%name), physical_const_list(i)%value, 0)
enddo

bp_com%i_const_init = nn
bp_com%i_const_tot  = nn

call indexer (bp_com2%const(1:nn)%name, bp_com2%const(1:nn)%index)

end subroutine init_bmad_parser_common

