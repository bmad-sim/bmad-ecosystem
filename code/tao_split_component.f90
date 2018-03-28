!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_split_component (comp_str, comp, err)
!
! Routine to split a component string.
!
! Input:
!   comp_str -- Character(*): Components. EG: 'meas - design'
!
! Output:
!   comp(:) -- Tao_data_var_component_struct, allocatable: Array of individual components.
!   err     -- Logical: Set True if there is an error, False otherwise.
!-

subroutine tao_split_component (comp_str, comp, err)

use tao_struct

implicit none

type (tao_data_var_component_struct), allocatable :: comp(:)

integer i, n, ix, ix1, ix2

character(*) comp_str
character(60) str
character(40) :: r_name = 'tao_split_component'

logical err

! Count number of components.

err = .true.
call string_trim (comp_str, str, ix)
if (.not. (str(1:1) == '+' .or. str(1:1) == '-')) str = '+' // trim(str)

n = 0
do i = 1, len_trim(str)
  if (str(i:i) == '+' .or. str(i:i) == '-') n = n + 1
enddo

! Allocate space and transfer info.

if (allocated(comp)) then
  if (size(comp) /= n) deallocate (comp)
endif
if (.not. allocated(comp)) allocate (comp(n))

n = 0
do n = 1, size(comp)
  if (str(1:1) == '+') then
    comp(n)%sign = 1
  elseif (str(1:1) == '-') then
    comp(n)%sign = -1
  else
    call out_io (s_error$, r_name, 'BAD COMPONENT LIST: ' // comp_str)
    return
  endif

  call string_trim(str(2:), str, ix)

  if (n == size(comp)) then
    comp(n)%name = str
    exit
  endif

  ix1 = index(str, "+")
  ix2 = index(str, "-")
  if (ix1 /= 0 .and. ix2 /= 0) then
    ix = min(ix1, ix2) - 1
  else
    ix = max(ix1, ix2) - 1
  endif

  comp(n)%name = str(1:ix)
  call string_trim(str(ix+1:), str, ix)

enddo

err = .false.

end subroutine tao_split_component
