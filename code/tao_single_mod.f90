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
logical print_message

!

ix_hash = index (out_file, '#')

do i = 1, size(s%u)

  file_name = out_file
  if (ix_hash /= 0) write (file_name(ix_hash:ix_hash), '(i1)') i

  iu = lunget()
  open (iu, file = file_name, carriagecontrol = 'list', iostat = ios)
    
  if (ios /= 0) then
    print *, 'ERROR IN VAR_WRITE: CANNOT OPEN FILE: ', trim(file_name)
    return
  endif
    
  do j = 1, size(s%var)
    if (.not. s%var(j)%exists) cycle
    if (.not. any (s%var(j)%this(:)%ix_uni == iu)) cycle
    write (iu, '(4a, g14.8)')  trim(s%var(j)%ele_name), &
              '[', trim(s%var(j)%attrib_name), '] = ', s%var(j)%model_value
  enddo
    
  write (iu, *)
  write (iu, *) 'end_file'
  write (iu, *)
    
  call tao_show_constraints (iu, 'ALL')
  call tao_show_constraints (iu, 'TOP10')
    
  close (iu)

  if (print_message) print *, 'Written: ', trim(file_name)

enddo

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_var_print_auto ()

implicit none

integer i
character(20) var_name

!

do i = 1, size(s%var)
    if (.not. s%var(i)%exists) cycle
    var_name = trim(s%var(i)%ele_name) // tao_var_uni_string(s%var(i))
    print '(4a, g14.8)', trim(var_name), '[', &
              trim(s%var(i)%attrib_name), '] = ', s%var(i)%model_value
enddo

end subroutine


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_var_print_manual ()

implicit none

integer i, j
character(20) var_name

!

do j = 1, size(s%key)
    i = s%key(j)%ix_var 
    if (.not. s%var(i)%exists) cycle
    var_name = trim(s%var(i)%ele_name) // tao_var_uni_string(s%var(i))
    print '(4a, g14.8)', trim(var_name), '[', &
              trim(s%var(i)%attrib_name), '] = ', s%var(i)%model_value
enddo

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

function tao_var_uni_string (var) result (str)

implicit none

type (tao_var_struct) var
character(10) str
integer i

!

str = ':'
do i = 1, size (var%this)
  write (str, '(a, i1)') trim(str), var%this(i)%ix_uni
enddo

end function

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_read_single_input_file (file_name)

implicit none


character(*) file_name

!

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_show_constraints (iunit, form)

use nr

implicit none

type (tao_top10_struct) top_merit(10)

real(rp) value
real(rp), allocatable :: merit(:)

integer i, j, iunit, nc, ir, n_max
integer ir1, ir2
integer, allocatable :: ixm(:)

character(*) form
character(16) location, con_var, max_loc, loc1, loc2
character(80) fmt
character(1) plane

type constraint_struct
  character(16) name
  character(16) loc1, loc2, max_loc
  real(rp) target_value
  real(rp) actual_value
  real(rp) merit
  integer plane
end type

type (constraint_struct), allocatable :: con(:)

!

fmt  = '(i3, 1x, a9, 1x, a1, 2x, a10, 1x, a9, 1pe10.2, 1pe12.3, e10.2, 2x, a8)'
top_merit(:)%valid  = .false.; top_merit(:)%name  = ' '

nc = count (s%var(:)%useit_opt)
do i = 1, size(s%u)
  nc = nc + count (s%u(i)%data(:)%useit_opt)
enddo
allocate (con(nc), ixm(nc))

do i = 1, size(s%u)
  do j = 1, size(s%u(i)%data)
    if (.not. s%u(i)%data(j)%useit_opt) cycle
    nc = nc + 1
  enddo
enddo

do i = 1, size(s%var)
enddo

!

if (form == 'TOP10') then
  call indexx(con(:)%merit, ixm)
  n_max = min(nc, 10)
  ixm(1:n_max) = ixm(nc:nc-n_max+1:-1)
  if (iunit == 0) then
    print *
    print *, '! Top 10'
  else
    write (iunit, *)
    write (iunit, *) '! Top 10'
  endif
elseif (form == 'ALL') then
  n_max = nc
  forall (i = 1:n_max) ixm(i) = i
else
  print *, 'ERROR IN SHOW_CONSTRAINTS: UNKNOWN FORM: ', form
  call err_exit
endif

if (iunit == 0) then
  print *
  print *,         'Constrnt Plane  Where1' //  &
                     '     Where2     Target     Value      Merit     Max'
else
  write (iunit, *)
  write (iunit, *) 'Constrnt Plane  Where1' //  &
                     '     Where2     Target     Value      Merit     Max'
endif

!

do j = 1, n_max

  i = ixm(j)
  if (iunit == 0) then
    print fmt, i, con(i)%name, &
            con(i)%plane, con(i)%loc1, con(i)%loc2, con(i)%target_value, &
            con(i)%actual_value, con(i)%merit, con(i)%max_loc
  else
    write (iunit, fmt) i, con(i)%name, &
            con(i)%plane, con(i)%loc1, con(i)%loc2, con(i)%target_value, &
            con(i)%actual_value, con(i)%merit, con(i)%max_loc
  endif

end do

!

if (iunit == 0) then
  print *
  print '(1x, a, 1pe12.6)', 'figure of merit: ', tao_merit()
else
  write (iunit, *)
  write (iunit, '(1x, a, 1pe12.6)') 'figure of merit: ', tao_merit()
endif

end subroutine

end module

