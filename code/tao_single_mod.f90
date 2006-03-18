module tao_single_mod

use tao_mod
use quick_plot
use tao_top10_mod

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_var_write (out_file, print_message)

implicit none


integer i, j, iu, ix, ios, ix_hash
character(*) out_file
character(80) file_name
character(20) :: r_name = 'tao_var_write'
logical print_message

!

ix_hash = index (out_file, '#')

do i = 1, size(s%u)

  file_name = out_file
  if (ix_hash /= 0) write (file_name(ix_hash:ix_hash), '(i0)') i

  iu = lunget()
  open (iu, file = file_name, carriagecontrol = 'list', iostat = ios)
    
  if (ios /= 0) then
    call out_io (s_error$, r_name, &
            'ERROR IN VAR_WRITE: CANNOT OPEN FILE: ' // file_name)
    return
  endif
    
  do j = 1, size(s%var)
    if (.not. s%var(j)%exists) cycle
    if (.not. any (s%var(j)%this(:)%ix_uni == i)) cycle
    write (iu, '(4a, es22.14)')  trim(s%var(j)%ele_name), &
              '[', trim(s%var(j)%attrib_name), '] = ', s%var(j)%model_value
  enddo
    
  write (iu, *)
  write (iu, *) 'end_file'
  write (iu, *)
    
  call tao_show_constraints (iu, 'ALL')
  call tao_show_constraints (iu, 'TOP10')
    
  close (iu)

!!  if (print_message) call out_io (s_blank$, r_name, 'Written: ' // file_name)
call out_io (s_blank$, r_name, 'Written: ' // file_name)

enddo

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

function tao_var_uni_string (var) result (str)

implicit none

type (tao_var_struct) var
character(10) str
integer i, iu, ct
logical uni(100)

!

iu = size(s%u)
uni = .false.

do i = 1, size (var%this)
  uni(var%this(i)%ix_uni) = .true.
enddo

ct = count(uni(1:iu))

if (ct == 1) then
  write (str, '(i2)') var%this(1)%ix_uni
elseif (ct == iu) then
  str = 'All'
else
  str = '?'
endif

end function

end module

