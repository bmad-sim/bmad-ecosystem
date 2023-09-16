!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine init_bmad_parser_common(lat)

use bmad_parser_struct

implicit none

type (lat_struct), optional :: lat
integer n1, n2, i

! 

if (allocated(bp_com2%const)) deallocate (bp_com2%const)

n1 = size(physical_const_list)  ! number of standard (non-user defined) constants
n2 = n1
if (present(lat)) then
  if (allocated(lat%constant)) n2 = n2 + size(lat%constant)
endif

allocate (bp_com2%const(n2))

do i = 1, n1
  bp_com2%const(i) = bp_const_struct(upcase(physical_const_list(i)%name), physical_const_list(i)%value, 0)
enddo

do i = 1, n2-n1
  bp_com2%const(i+n1) = bp_const_struct(lat%constant(i)%name, lat%constant(i)%value, 0)
enddo

bp_com%i_const_init = n1
bp_com%i_const_tot  = n2

call indexer (bp_com2%const(1:n2)%name, bp_com2%const(1:n2)%index)

end subroutine init_bmad_parser_common

