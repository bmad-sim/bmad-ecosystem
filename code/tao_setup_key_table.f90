!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Put the variables marked by key_bound in the key table.

subroutine tao_setup_key_table ()

use tao_struct

implicit none

integer i, j

! Key table setup

call re_allocate (s%key, count(s%var%key_bound))
s%key = -1

j = 0
do i = 1, s%n_var_used
  if (s%var(i)%key_bound .and. s%var(i)%exists) then
    j = j + 1
    s%key(j) = i
    s%var(i)%key_val0 = s%var(i)%model_value
    s%var(i)%ix_key_table = j
  else
    s%var(i)%key_bound = .false.
    s%var(i)%ix_key_table = -1
  endif
enddo

end subroutine tao_setup_key_table

