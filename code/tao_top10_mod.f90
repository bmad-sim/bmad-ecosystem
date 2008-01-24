module tao_top10_mod

use tao_struct
use tao_interface
use cesr_utils
use tao_dmerit_mod

! structure for making lists of the biggest contributors to the merit function.

type tao_top10_struct
  character(40) name   ! name of contributor
  real(rp) value       ! contribution to the merit function
  integer index        ! index of contributor.
  logical valid        ! valid entry?
end type

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_top10_print ()
!
! Routine to print out the top10 contributors to the merit function.
!
! Input:
!-

subroutine tao_top10_print ()

implicit none

type (tao_top10_struct) top_merit(10)
type (tao_top10_struct) top_dmerit(10)
type (tao_top10_struct) top_delta(10)
type (tao_data_struct), pointer :: data

real(rp) delta, a_max, merit
integer i, j, n, nl, nu

character(40) name
character(100) fmt, lines(20)
character(20) :: r_name = 'tao_top10_print'

! tao_merit also calculates the contrribution of the individual
! variables and data to the merit function.

merit = tao_merit()
call tao_dmerit_calc ()

! top_merit stores the top contributors to the merit function.
! top_dmerit stores the top dmerit/dvar values
! top_delta stores the top |var_model - var_design| 

top_merit(:)%valid  = .false.; top_merit(:)%name  = ' '
top_merit(:)%value = 0; top_merit(:)%index = 0
top_dmerit(:)%valid = .false.; top_dmerit(:)%name = ' '
top_dmerit(:)%value = 0; top_dmerit(:)%index = 0
top_delta(:)%valid  = .false.; top_delta(:)%name  = ' '
top_delta(:)%value = 0; top_delta(:)%index = 0

nu = ubound(s%u, 1)
do i = 1, nu
  do j = 1, size(s%u(i)%data)
    data => s%u(i)%data(j)
    if (.not. data%useit_opt) cycle
    name = trim(data%d1%d2%name) // '.' // data%d1%name
    if (nu > 1) write (name, '(i0, 2a)') i, '@', trim(name)
    call tao_to_top10 (top_merit, data%merit, name, data%ix_d1, 'max')
  enddo
enddo


do j = 1, size(s%var)
  if (.not. s%var(j)%useit_opt) cycle
  name = s%var(j)%v1%name
  call tao_to_top10 (top_merit, s%var(j)%merit, name, s%var(j)%ix_v1, 'max')
  call tao_to_top10 (top_dmerit, s%var(j)%dmerit_dvar, name, &
                                                    s%var(j)%ix_v1, 'abs_max')
  delta = s%var(j)%model_value - s%var(j)%design_value
  call tao_to_top10 (top_delta, delta, name, s%var(j)%ix_v1, 'abs_max')
enddo

! write results

call tao_show_constraints (0, 'TOP10')

a_max = max(1.1, maxval(abs(top_delta(:)%value)))
n = max(0, 6 - int(log10(a_max)))

write (fmt, '(a, i1, a)') '((a10, i5, 1pe12.3, 3x), (a10, i5, 0pf11.', n, '))'

nl = 0
lines(nl+1) = ' '
lines(nl+2) = '      Top10 derivative       |       Top10 delta'
lines(nl+3) = ' Name         ix  Derivative | Name         ix     delta'
nl = nl + 3

do i = 1, 10
  nl = nl + 1
  write (lines(nl), fmt) &
      top_dmerit(i)%name, top_dmerit(i)%index, top_dmerit(i)%value,  &
      top_delta(i)%name,  top_delta(i)%index,  top_delta(i)%value
enddo

nl = nl + 1
write (lines(nl), *) 'Merit:  ', merit

call out_io (s_blank$, r_name, lines(1:nl))

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_to_top10 (top10, value, name, c_index, order)
!
! Routine to order the largest contributors to the merit function in
! a list. Call this routine for each contributor.
!
! Note: Before first calling this routine set:
!   top10(:)%valid = .false.
!
! Input:
!   value   -- Real(rp): value of the contributor.
!   name    -- Character(40): Name of the contributor..
!   c_index -- Integer: Index of the contributor.
!   order   -- Character(16): Ordering of the list. Possibilities are:
!                 'max'     -- #1 has the maximum value.
!                 'min'     -- #1 has the minimum value.
!                 'abs_max' -- #1 has the maximum aplitude.
!                 'abs_min' -- #1 has the maximum aplitude.
!
! Output:
!   top10(:) -- Tao_top10_struct: List of top contributors.
!                 Note that the list is not limited to 10 entries.
!-

subroutine tao_to_top10 (top10, value, name, c_index, order)

implicit none

type (tao_top10_struct) top10(:)

integer c_index, ix, n
real(rp) value

character(*) name, order
character(20) :: r_name = 'tao_to_top10'

! Find where in list the current contributor is.

n = size(top10)
do ix = n, 1, -1
  if (.not. top10(ix)%valid) cycle
  select case (order)
  case ('max')
    if (value < top10(ix)%value) exit
  case ('min')
    if (value > top10(ix)%value) exit
  case ('abs_max')  
    if (abs(value) < abs(top10(ix)%value)) exit
  case ('abs_min')  
    if (abs(value) > abs(top10(ix)%value)) exit
  case default
    call out_io (s_abort$, r_name, 'BAD "ORDER" ARGUMENT: ' // order)
  end select
enddo

ix = ix + 1          ! place to put current contributor.
if (ix > n) return   ! not big enough to be in list.

! Move the people below the current contributor down to make room and
! then put the contributor in.

top10(ix+1:n) = top10(ix:n-1) 
top10(ix) = tao_top10_struct(name, value, c_index, .true.)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_show_constraints (iunit, form)

use nr

implicit none

type (tao_top10_struct) top_merit(10)
type (tao_var_struct), pointer :: var
type (tao_data_struct), pointer :: data
type (tao_universe_struct), pointer :: u

real(rp) value, this_merit

integer i, j, n, iunit, nc, ir, n_max
integer ir1, ir2, iu, ie, nl
integer, allocatable, save :: ixm(:)
integer n_name, n_d2_d1_name, n_loc1, n_loc0

character(*) form
character(40) location, con_var, max_loc, loc1, loc0
character(80) fmt
character(1) plane
character(24) :: r_name = 'tao_show_constraints'
character(200), allocatable, save :: line(:)
character(200) l1

type constraint_struct
  character(40) d2_d1_name
  character(40) name
  character(40) loc1, loc0, max_loc
  real(rp) target_value
  real(rp) actual_value
  real(rp) merit
end type

type (constraint_struct), allocatable, save :: con(:)

! Init
 
call re_allocate (line, 200, 100)
this_merit = tao_merit()
top_merit(:)%valid  = .false.; top_merit(:)%name  = ' '

nc = count (s%var(:)%useit_opt)
do i = lbound(s%u, 1), ubound(s%u, 1)
  nc = nc + count (s%u(i)%data(:)%useit_opt)
enddo
allocate (con(nc), ixm(nc))

! Data constraints

nc = 0
do i = lbound(s%u, 1), ubound(s%u, 1)
  do j = 1, size(s%u(i)%data)
    data => s%u(i)%data(j)
    if (.not. data%useit_opt) cycle
    nc = nc + 1
    con(nc)%name = tao_datum_type_name(data)
    con(nc)%d2_d1_name = trim(tao_datum_name(data))

    if (data%ix_ele < 0) then
      con(nc)%loc1 = ' '
    else
      con(nc)%loc1 = data%ele_name  ! s%u(i)%model%lat%ele(data%ix_ele)%name
    endif

    con(nc)%loc0 = data%ele0_name

    ie = data%ix_ele_merit
    if (ie < 0) then
      con(nc)%max_loc = ' '
    else
      con(nc)%max_loc = s%u(i)%model%lat%ele(ie)%name
    endif

    con(nc)%target_value = data%meas_value
    con(nc)%actual_value = data%model_value
    con(nc)%merit = data%merit

  enddo
enddo

! Variable constraints

do i = 1, size(s%var(:))
  var => s%var(i)
  if (.not. var%useit_opt) cycle
  nc = nc + 1
  con(nc)%d2_d1_name = trim(tao_var1_name(var))
  con(nc)%name = trim(tao_var_attrib_name(var))
  u => s%u(var%this(1)%ix_uni)
  write (con(nc)%loc1, '(f8.2)') u%model%lat%ele(var%this(1)%ix_ele)%s
  con(nc)%loc0 = ''
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
  con(nc)%max_loc = ''
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

n_d2_d1_name = 9;  n_name = 1; n_loc1 = 8;  n_loc0 = 8

do j = 1, n_max
  i = ixm(j)
  n_d2_d1_name = max(n_d2_d1_name, len_trim(con(i)%d2_d1_name))
  n_name = max(n_name, len_trim(con(i)%name))
  n_loc1 = max(n_loc1, len_trim(con(i)%loc1))
  n_loc0 = max(n_loc0, len_trim(con(i)%loc0))
enddo

!

l1 = 'Constraint'
n=3+n_d2_d1_name+2+n_name; l1(n:) = 'Where0'
n=len_trim(l1)+n_loc0-3;   l1(n:) = 'Where'
n=len_trim(l1)+n_loc1-2;   l1(n:) = 'Target     Value      Merit     Max'

nl=nl+1; line(nl) = ' '
nl=nl+1; line(nl) = l1

!

fmt = '(a, 2x, a, 2x, a, 2x, a, 1pe10.2, 1pe12.3, e10.2, 2x, a)'

call re_allocate (line, 200, nl+n_max+100)
do j = 1, n_max
  i = ixm(j)
  nl=nl+1; write (line(nl), fmt) con(i)%d2_d1_name(1:n_d2_d1_name), &
            con(i)%name(1:n_name), &
            con(i)%loc0(1:n_loc0), con(i)%loc1(1:n_loc1), con(i)%target_value, &
            con(i)%actual_value, con(i)%merit, con(i)%max_loc
end do
nl=nl+1; line(nl) = l1

!

nl=nl+1; line(nl) = ' '
nl=nl+1; write (line(nl), '(1x, a, 1pe12.6)') &
                                  'figure of merit: ', this_merit

call tao_write_out (iunit, line, nl)

deallocate (con, ixm)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_var_write (out_file, good_opt_vars_only)

implicit none


integer i, j, iu, ix, ios, ix_hash
character(*) out_file
character(200) file_name, str(1)
character(20) :: r_name = 'tao_var_write'
logical, optional :: good_opt_vars_only
logical printed

!

ix_hash = index (out_file, '#')
printed = .false.

do i = lbound(s%u, 1), ubound(s%u, 1)

  if (.not. s%u(i)%is_on) cycle

  file_name = out_file
  if (ix_hash /= 0) write (file_name, '(a, i0, a)') &
                  file_name(1:ix_hash-1), i, trim(file_name(ix_hash+1:))

  if (file_name == ' ') then
    iu = 0
  else
    iu = lunget()
    open (iu, file = file_name, carriagecontrol = 'list', recl = 100, iostat = ios)
    if (ios /= 0) then
      call out_io (s_error$, r_name, &
            'ERROR IN VAR_WRITE: CANNOT OPEN FILE: ' // file_name)
      return
    endif
  endif

  ! If printing then only write vars once.

  if (iu == 0 .and. printed) return

  do j = 1, size(s%var)
    if (.not. s%var(j)%exists) cycle
    if (iu /= 0 .and. .not. any (s%var(j)%this(:)%ix_uni == i)) cycle
    if (logic_option(.false., good_opt_vars_only) .and. .not. s%var(j)%useit_opt) cycle
    write (str(1), '(4a, es22.14)')  trim(s%var(j)%ele_name), &
              '[', trim(s%var(j)%attrib_name), '] = ', s%var(j)%model_value
    call tao_write_out (iu, str, 1)
  enddo
  
  if (iu /= 0) then
    call tao_write_out (iu, (/ '        ', 'end_file', '        ' /), 3)
    call tao_show_constraints (iu, 'ALL')
    call tao_show_constraints (iu, 'TOP10')
    close (iu)
    call out_io (s_blank$, r_name, 'Written: ' // file_name)
  endif

  printed = .true.

enddo

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_write_out (iunit, line, nl)
!
! Subroutine to write out a series of lines to a file or to the terminal.
! It is assumed that any file has already been opened.
!
! Input:
!   iunit   -- Integer: File unit to write to. 0 => print to the terminal.
!   line(:) -- Character(*): A series of lines.
!   nl      -- Integer: The number of lines to write.
!-

subroutine tao_write_out (iunit, line, nl)

implicit none

integer iunit, i, nl
character(*) line(:)

!

do i = 1, nl
  if (iunit == 0) then
    call out_io (s_blank$, ' ', line(i))
  else
    write (iunit, *) trim(line(i))
  endif
enddo

end subroutine

end module
