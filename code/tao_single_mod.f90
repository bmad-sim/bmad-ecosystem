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
  if (ix_hash /= 0) write (file_name(ix_hash:ix_hash), '(i1)') i

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
    write (iu, '(4a, g14.8)')  trim(s%var(j)%ele_name), &
              '[', trim(s%var(j)%attrib_name), '] = ', s%var(j)%model_value
  enddo
    
  write (iu, *)
  write (iu, *) 'end_file'
  write (iu, *)
    
  call tao_show_constraints (iu, 'ALL')
  call tao_show_constraints (iu, 'TOP10')
    
  close (iu)

  if (print_message) call out_io (s_blank$, r_name, 'Written: ' // file_name)

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

str = ' '
do i = 1, size (var%this)
  write (str, '(a, i1)') trim(str), var%this(i)%ix_uni
enddo

end function

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
integer ir1, ir2, iu, ie, nl
integer, allocatable :: ixm(:)

character(*) form
character(16) location, con_var, max_loc, loc1, loc2
character(80) fmt
character(1) plane
character(24) :: r_name = 'tao_show_constraints'
character(100) line(100)

type constraint_struct
  character(16) name
  character(16) loc1, loc2, max_loc
  real(rp) target_value
  real(rp) actual_value
  real(rp) merit
end type

type (constraint_struct), allocatable :: con(:)

!

fmt  = '(i3, 1x, a9, 2x, a10, 1x, a9, 1pe10.2, 1pe12.3, e10.2, 2x, a8)'
top_merit(:)%valid  = .false.; top_merit(:)%name  = ' '

nc = count (s%var(:)%useit_opt)
do i = 1, size(s%u)
  nc = nc + count (s%u(i)%data(:)%useit_opt)
enddo
allocate (con(nc), ixm(nc))

nc = 0
do i = 1, size(s%u)
  do j = 1, size(s%u(i)%data)
    if (.not. s%u(i)%data(j)%useit_opt) cycle
    nc = nc + 1
    con(nc)%name = s%u(i)%data(j)%data_type
    if (size(s%u) > 1) write (con(nc)%name, '(2a, i1)') &
                                     trim(con(nc)%name), ';', i
    con(nc)%loc1 = s%u(i)%model%ele_(s%u(i)%data(j)%ix_ele)%name
    ie = s%u(i)%data(j)%ix_ele2
    if (ie < 0) then
      con(nc)%loc2 = ' '
    else
      con(nc)%loc2 = s%u(i)%model%ele_(ie)%name
    endif
    ie = s%u(i)%data(j)%ix_ele_merit
    if (ie < 0) then
      con(nc)%max_loc = ' '
    else
      con(nc)%max_loc = s%u(i)%model%ele_(ie)%name
    endif
    con(nc)%target_value = s%u(i)%data(j)%meas_value
    con(nc)%actual_value = s%u(i)%data(j)%model_value
    con(nc)%merit = s%u(i)%data(j)%merit
  enddo
enddo

do i = 1, size(s%var(:))
  if (.not. s%var(i)%useit_opt) cycle
  if (s%var(i)%merit == 0) cycle
  nc = nc + 1
  con(nc)%name = s%var(i)%name
  iu = s%var(i)%this(1)%ix_uni
  con(nc)%loc1 = s%u(iu)%model%ele_(s%var(i)%this(1)%ix_ele)%name
  con(nc)%loc2 = ' '
  if (abs(s%var(i)%model_value - s%var(i)%high_lim) < &
      abs(s%var(i)%model_value - s%var(i)%low_lim)) then
    con(nc)%target_value = s%var(i)%high_lim
  else
    con(nc)%target_value = s%var(i)%low_lim
  endif
  con(nc)%actual_value = s%var(i)%model_value
  con(nc)%merit = s%var(i)%merit
  con(nc)%max_loc = ' '
enddo

!

if (form == 'TOP10') then
  call indexx(con(:)%merit, ixm)
  n_max = min(nc, 10)
  ixm(1:n_max) = ixm(nc:nc-n_max+1:-1)
  line(1) = ' '
  line(2) = '! Top 10'
  nl = 2
elseif (form == 'ALL') then
  n_max = nc
  forall (i = 1:n_max) ixm(i) = i
  nl = 0
else
  call out_io (s_abort$, r_name, &
              'ERROR IN SHOW_CONSTRAINTS: UNKNOWN FORM: ' // form)
  call err_exit
endif

  nl=nl+1; line(nl) = ' '
  nl=nl+1; line(nl) = 'Constrnt      Where1' //  &
                     '     Where2     Target     Value      Merit     Max'

!

do j = 1, n_max
  i = ixm(j)
  nl = nl + 1
  write (line(nl), fmt) i, con(i)%name, &
            con(i)%loc1, con(i)%loc2, con(i)%target_value, &
            con(i)%actual_value, con(i)%merit, con(i)%max_loc
end do

!

  nl=nl+1; line(nl) = ' '
  nl=nl+1; write (line(nl), '(1x, a, 1pe12.6)') &
                                  'figure of merit: ', tao_merit()

!

do i = 1, nl
  if (iunit == 0) then
    call out_io (s_blank$, r_name, line(i))
  else
    write (iunit, *) line(i)
  endif
enddo

end subroutine

end module

