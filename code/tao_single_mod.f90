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

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_show_constraints (iunit, form)

use nr

implicit none

type (tao_top10_struct) top_merit(10)
type (tao_var_struct), pointer :: var
type (tao_data_struct), pointer :: data

real(rp) value, this_merit

integer i, j, n, iunit, nc, ir, n_max
integer ir1, ir2, iu, ie, nl
integer, allocatable :: ixm(:)
integer n_name, n_d2_d1_name, n_loc1, n_loc2

character(*) form
character(16) location, con_var, max_loc, loc1, loc2
character(80) fmt
character(1) plane
character(24) :: r_name = 'tao_show_constraints'
character(200), allocatable, save :: line(:)
character(200) l1

type constraint_struct
  character(32) d2_d1_name
  character(32) name
  character(16) loc1, loc2, max_loc
  real(rp) target_value
  real(rp) actual_value
  real(rp) merit
  integer ix
end type

type (constraint_struct), allocatable :: con(:)

!

call re_allocate (line, 200, 100)
this_merit = tao_merit()
top_merit(:)%valid  = .false.; top_merit(:)%name  = ' '

nc = count (s%var(:)%useit_opt)
do i = 1, size(s%u)
  nc = nc + count (s%u(i)%data(:)%useit_opt)
enddo
allocate (con(nc), ixm(nc))

nc = 0
do i = 1, size(s%u)
  do j = 1, size(s%u(i)%data)
    data => s%u(i)%data(j)
    if (.not. data%useit_opt) cycle
    nc = nc + 1
    con(nc)%name = data%name
    con(nc)%d2_d1_name = trim(data%d1%d2%name) // ':' // data%d1%name
    if (size(s%u) > 1) write (con(nc)%d2_d1_name, '(2a, i0)') &
                                     trim(con(nc)%d2_d1_name), ';', i
    if (data%ix_ele < 0) then
      con(nc)%loc1 = ' '
    else
      con(nc)%loc1 = s%u(i)%model%lat%ele_(data%ix_ele)%name
    endif
    ie = data%ix_ele2
    if (ie < 1) then
      con(nc)%loc2 = ' '
    else
      con(nc)%loc2 = s%u(i)%model%lat%ele_(ie)%name
    endif
    ie = data%ix_ele_merit
    if (ie < 0) then
      con(nc)%max_loc = ' '
    else
      con(nc)%max_loc = s%u(i)%model%lat%ele_(ie)%name
    endif
    con(nc)%target_value = data%meas_value
    con(nc)%actual_value = data%model_value
    con(nc)%merit = data%merit
    con(nc)%ix = data%ix_d1
  enddo
enddo

do i = 1, size(s%var(:))
  var => s%var(i)
  if (.not. var%useit_opt) cycle
!!  if (var%merit == 0) cycle
  nc = nc + 1
  con(nc)%d2_d1_name = var%v1%name
  con(nc)%name       = var%name
  iu = var%this(1)%ix_uni
  con(nc)%loc1 = s%u(iu)%model%lat%ele_(var%this(1)%ix_ele)%name
  con(nc)%loc2 = ' '
  if (var%merit_type == 'target') then
    con(nc)%target_value = var%meas_value
  elseif (var%merit_type == 'limit') then
    if (abs(var%model_value - var%high_lim) < &
                  abs(var%model_value - var%low_lim)) then
      con(nc)%target_value = var%high_lim
    else
      con(nc)%target_value = var%low_lim
    endif
  else
    call out_io (s_abort$, r_name, 'UNKNOWN VARIABLE MERIT_TYPE: ' // &
                                                            var%merit_type)
  endif
  con(nc)%actual_value = var%model_value
  con(nc)%merit = var%merit
  con(nc)%max_loc = ' '
  con(nc)%ix = var%ix_v1
enddo

!

if (form == 'TOP10') then
  call indexx(con(1:nc)%merit, ixm(1:nc))
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

! find string widths

n_d2_d1_name = 9;  n_name = 1; n_loc1 = 8;  n_loc2 = 8

do j = 1, n_max
  i = ixm(j)
  n_d2_d1_name = max(n_d2_d1_name, len_trim(con(i)%d2_d1_name))
  n_name = max(n_name, len_trim(con(i)%name))
  n_loc1 = max(n_loc1, len_trim(con(i)%loc1))
  n_loc2 = max(n_loc2, len_trim(con(i)%loc2))
enddo

!

l1 = 'Constraint'
n=8+n_d2_d1_name+2+n_name; l1(n:) = 'Where1'
n=len_trim(l1)+n_loc1-4;   l1(n:) = 'Where2'
n=len_trim(l1)+n_loc2-2;   l1(n:) = 'Target       Value     Merit   Max'

nl=nl+1; line(nl) = ' '
nl=nl+1; line(nl) = l1

!

fmt = '(a, i5, 2x, a, 2x, a, 1x, a, 1pe10.2, 1pe12.3, e10.2, 2x, a)'

call re_allocate (line, 200, nl+n_max+100)
do j = 1, n_max
  i = ixm(j)
  nl=nl+1; write (line(nl), fmt) con(i)%d2_d1_name(1:n_d2_d1_name), con(i)%ix, &
            con(i)%name(1:n_name), &
            con(i)%loc1(1:n_loc1), con(i)%loc2(1:n_loc2), con(i)%target_value, &
            con(i)%actual_value, con(i)%merit, con(i)%max_loc
end do
nl=nl+1; line(nl) = l1

!

nl=nl+1; line(nl) = ' '
nl=nl+1; write (line(nl), '(1x, a, 1pe12.6)') &
                                  'figure of merit: ', this_merit

!

do i = 1, nl
  if (iunit == 0) then
    call out_io (s_blank$, r_name, line(i))
  else
    write (iunit, *) trim(line(i))
  endif
enddo

end subroutine

end module

