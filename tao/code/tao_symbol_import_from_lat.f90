!+
! Subroutine tao_symbol_import_from_lat (lat)
!
! Routine to import symbols from a lattic

subroutine tao_symbol_import_from_lat (lat)

use tao_set_mod, dummy => tao_symbol_import_from_lat

implicit none
type (lat_struct) lat

integer i

!

if (.not. allocated(lat%constant)) return

do i = 1, size(lat%constant)
  call tao_set_symbolic_number_cmd (downcase(lat%constant(i)%name), val = lat%constant(i)%value)
enddo

end subroutine
