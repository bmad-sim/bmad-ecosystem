module tao_top10_mod

use tao_interface
use tao_dmerit_mod
use sim_utils

! structure for making lists of the biggest contributors to the merit function.

type tao_top10_struct
  character(40) :: name = ''   ! name of contributor
  real(rp) :: value = 0        ! contribution to the merit function
  integer :: index = 0         ! index of contributor.
  logical :: valid = .false.   ! valid entry?
end type

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_top10_merit_categories_print (iunit)
!
! Routine to print the top data and variable categories that contribute to
! the merit function.
!
! Input:
!   iunit   -- Integer: File unit to write to. 0 => print to the terminal.
!-

subroutine tao_top10_merit_categories_print (iunit)

implicit none

type (tao_top10_struct), allocatable :: tao_merit(:), temp_merit(:)
type (tao_d1_data_struct), pointer :: d1

integer i, j, k, m, n, id, nl, num, iunit
integer, allocatable :: n_points(:)
real(rp), allocatable :: chi2(:)

character(100) fmt, lines(100)
character(40) name
character(*), parameter :: r_name = 'tao_top10_print'

! Need to deal with the fact that the merit value of a datum or variable may be negative.

n = 0
do i = lbound(s%u, 1), ubound(s%u, 1)
  do j = 1, s%u(i)%n_d2_data_used
    n = n + size(s%u(i)%d2_data(j)%d1(:))
  enddo
enddo

allocate (temp_merit(n), n_points(n), chi2(n))

temp_merit(:)%value = 0
chi2(:) = 0
n_points(:) = 0

num = 0
do i = lbound(s%u, 1), ubound(s%u, 1)
  do j = 1, s%u(i)%n_d2_data_used
    do k = 1, size(s%u(i)%d2_data(j)%d1)
      d1 => s%u(i)%d2_data(j)%d1(k)
      name = trim(s%u(i)%d2_data(j)%name) // '.' // trim(d1%name)

      do m = 1, num
        if (name == temp_merit(m)%name) exit
      enddo

      if (m == num + 1) then
        num = num + 1
        temp_merit(num)%name = name
      endif

      do id = lbound(d1%d, 1), ubound(d1%d, 1)
        if (d1%d(id)%weight == 0) cycle
        if (d1%d(id)%merit == 0) cycle
        temp_merit(m)%value = temp_merit(m)%value + d1%d(id)%merit
        chi2(m) = chi2(m) + d1%d(id)%delta_merit**2
        n_points(m) = n_points(m) + 1
      enddo
    enddo
  enddo
enddo

allocate (tao_merit(num+s%n_v1_var_used))

do i = 1, s%n_v1_var_used
  call tao_to_top10 (tao_merit, abs(sum(s%v1_var(i)%v(:)%merit)), s%v1_var(i)%name, -i, 'max')
enddo

do i = 1, num
  call tao_to_top10 (tao_merit, abs(temp_merit(i)%value), temp_merit(i)%name, i, 'max')
enddo

!

nl = 0
nl=nl+1; lines(nl) = ''
nl=nl+1; lines(nl) = 'List of non-zero contributors to the Merit Function:'
nl=nl+1; write (lines(nl), '(a, 42x, a)') 'Name', 'Merit           Sigma [= sqrt(Chi^2/N)]'
do i = 1, size(tao_merit)
  if (.not. tao_merit(i)%valid) exit
  if (tao_merit(i)%value == 0) exit
  m = tao_merit(i)%index
  if (m < 0) then
    nl=nl+1; write (lines(nl), '(a, 2es16.4)') tao_merit(i)%name, sum(s%v1_var(-m)%v(:)%merit)
  else
    nl=nl+1; write (lines(nl), '(a, 2es16.4)') tao_merit(i)%name, temp_merit(m)%value, sqrt(chi2(m)/n_points(m))
  endif
enddo

call tao_write_lines (iunit, lines(1:nl))

deallocate (tao_merit, temp_merit, n_points, chi2)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_top10_derivative_print ()
!
! Routine to print out the top10 contributors to the merit function.
!-

subroutine tao_top10_derivative_print ()

implicit none

type (tao_top10_struct) top_dmerit(s%global%n_top10_merit)
type (tao_top10_struct) top_delta(s%global%n_top10_merit)
type (tao_data_struct), pointer :: data

real(rp) delta, a_max, merit
integer i, j, n, nl, nu

character(40) name
character(100) fmt
character(100), allocatable :: lines(:)
character(20) :: r_name = 'tao_top10_print'

! tao_merit also calculates the contrribution of the individual
! variables and data to the merit function.

merit = tao_merit()
call tao_dmerit_calc ()

! top_dmerit stores the top dmerit/dvar values
! top_delta stores the top |var_model - var_design| 

top_dmerit(:)%valid = .false.; top_dmerit(:)%name = ' '
top_dmerit(:)%value = 0; top_dmerit(:)%index = 0
top_delta(:)%valid  = .false.; top_delta(:)%name  = ' '
top_delta(:)%value = 0; top_delta(:)%index = 0

do j = 1, s%n_var_used
  if (.not. s%var(j)%useit_opt) cycle
  name = s%var(j)%v1%name
  call tao_to_top10 (top_dmerit, s%var(j)%dmerit_dvar, name, &
                                                    s%var(j)%ix_v1, 'abs_max')
  delta = s%var(j)%model_value - s%var(j)%design_value
  call tao_to_top10 (top_delta, delta, name, s%var(j)%ix_v1, 'abs_max')
enddo

! write results

a_max = max(1.1, maxval(abs(top_delta(:)%value)))
n = max(0, 6 - int(log10(a_max)))

write (fmt, '(a, i1, a)') '((a15, es12.3, 3x), (a15, f11.', n, '))'
allocate (lines(10+s%global%n_top10_merit))

nl = 0
nl=nl+1; lines(nl) = ' '
nl=nl+1; lines(nl) = '      Top10 derivative       |       Top10 delta'
nl=nl+1; lines(nl) = ' Name            Derivative  | Name                Delta'

do i = 1, s%global%n_top10_merit
  if (.not. top_dmerit(i)%valid .and. .not. top_delta(i)%valid) exit
  nl=nl+1; write (lines(nl), fmt) &
        this_name(top_dmerit(i)%name, top_dmerit(i)%index), top_dmerit(i)%value,  &
        this_name(top_delta(i)%name,  top_delta(i)%index),  top_delta(i)%value
enddo

nl=nl+1; write (lines(nl), *) 'Merit:  ', merit

call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------------
contains

function this_name (name, ix) result (full_name)

character(*) name
character(15) full_name
integer ix

full_name = trim(name) // '[' // int_str(ix) // ']'

end function this_name

end subroutine tao_top10_derivative_print

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
    call out_io (s_error$, r_name, 'BAD "ORDER" ARGUMENT: ' // order)
  end select
enddo

ix = ix + 1          ! place to put current contributor.
if (ix > n) return   ! not big enough to be in list.

! Move the people below the current contributor down to make room and
! then put the contributor in.

top10(ix+1:n) = top10(ix:n-1) 
top10(ix) = tao_top10_struct(name, value, c_index, .true.)

end subroutine tao_to_top10

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_show_constraints (iunit, form)
!
! Routine to show a list of datums and variables and how they contribute to the merit function.
!
! Input:
!   iunit   -- integer: File unit to write to. 0 => print to the terminal.
!   form    -- character(*): What to output:
!                 'ALL'   -> All datums and variables.
!                 'TOP10' -> Top datums and variables that contribute to the merit function.
!-

subroutine tao_show_constraints (iunit, form)

implicit none

type (tao_top10_struct) top_merit(s%global%n_top10_merit)
type (tao_var_struct), pointer :: var
type (tao_data_struct), pointer :: data
type (tao_universe_struct), pointer :: u
type (branch_struct), pointer :: branch

real(rp) value, this_merit

integer i, j, n, ix, iunit, nc, ir, n_max
integer ir1, ir2, iu, ie, nl, ct
integer, allocatable :: ixm(:)
integer n_name, n_d2_d1_name, n_loc_ele, n_loc_start, n_loc_ref, n_tot

character(*) form
character(40) location, con_var, max_loc, loc_ele, loc_start, loc_ref
character(80) fmt, fmt2
character(1) plane
character(24) :: r_name = 'tao_show_constraints'
character(400) name
character(300), allocatable :: line(:)
character(300) l1

type constraint_struct
  character(40) d2_d1_name
  character(80) name
  character(40) loc_ele, loc_start, loc_ref, max_loc
  real(rp) target_value
  real(rp) actual_value
  real(rp) merit
  logical :: expression = .false.
end type

type (constraint_struct), allocatable, target :: con(:)
type (constraint_struct), pointer :: c

! Init
 
call re_allocate (line, 100)
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
    name = tao_constraint_type_name(data)
    n = len_trim(name)
    if (n > len(con(nc)%name)) then
      con(nc)%name = name(1:40) // ' ... ' // name(n-34:n)
    else
      con(nc)%name = name
    endif
    con(nc)%d2_d1_name = trim(tao_datum_name(data))
    if (substr(data%data_type,1,11) == 'expression:') con(nc)%expression = .true.

    if (data%ix_ele < 0) then
      con(nc)%loc_ele = ' '
    else
      con(nc)%loc_ele = data%ele_name 
    endif

    con(nc)%loc_start = data%ele_start_name
    con(nc)%loc_ref = data%ele_ref_name

    branch => s%u(i)%model%lat%branch(data%ix_branch)
    ie = data%ix_ele_merit

    if (ie < 0) then
      con(nc)%max_loc = ''
    else
      con(nc)%max_loc = branch%ele(ie)%name
    endif

    con(nc)%target_value = data%meas_value
    con(nc)%actual_value = data%model_value
    con(nc)%merit = data%merit

  enddo
enddo

! Variable constraints

do i = 1, s%n_var_used
  var => s%var(i)
  if (.not. var%useit_opt) cycle
  nc = nc + 1
  con(nc)%d2_d1_name = trim(tao_var1_name(var))
  con(nc)%name = trim(tao_var_attrib_name(var))
  u => s%u(var%slave(1)%ix_uni)

  con(nc)%loc_ref = ''
  con(nc)%loc_start = ''
  con(nc)%loc_ele = ''

  if (var%slave(1)%ix_ele < 0) then  ! EG particle_start
    con(nc)%loc_ele = '0'
  else
    branch => u%model%lat%branch(var%slave(1)%ix_branch)
    ct = branch%ele(var%slave(1)%ix_ele)%lord_status
    if (ct /= group_lord$ .and. ct /= overlay_lord$ .and. ct /= multipass_lord$) then
      write (con(nc)%loc_ele, '(f8.2)') branch%ele(var%slave(1)%ix_ele)%s
      call string_trim (con(nc)%loc_ele, con(nc)%loc_ele, ix)
    endif
  endif

  if (var%merit_type == 'target') then
    con(nc)%target_value = var%meas_value
  elseif (var%merit_type == 'limit') then
    if (abs(var%model_value - var%high_lim) < abs(var%model_value - var%low_lim)) then
      con(nc)%target_value = var%high_lim
    else
      con(nc)%target_value = var%low_lim
    endif
  else
    call out_io (s_error$, r_name, 'UNKNOWN VARIABLE MERIT_TYPE: ' // var%merit_type)
  endif
  con(nc)%actual_value = var%model_value
  con(nc)%merit = var%merit
  con(nc)%max_loc = ''
enddo

!

if (form == 'TOP10') then
  ! Merit may be negative when trying to maximize a datum instead of minimizing it.
  call indexer(abs(con(1:nc)%merit), ixm(1:nc))
  n_max = min(nc, s%global%n_top10_merit)
  ixm(1:n_max) = ixm(nc:nc-n_max+1:-1)
  line(1) = ' '
  line(2) = '! Top 10'
  nl = 2
elseif (form == 'ALL') then
  n_max = nc
  forall (i = 1:n_max) ixm(i) = i
  nl = 0
else
  call out_io (s_abort$, r_name, 'UNKNOWN FORM: ' // form)
  return
endif

! Find string widths
! Expressions generally have very long strings so we let this spill over to
! the where_start, where_ref and where fields

n_d2_d1_name = 9;  n_name = 1; n_loc_ele = 8;  n_loc_start = 10; n_loc_ref = 8; n_tot = 0

do j = 1, n_max
  i = ixm(j)
  n_d2_d1_name = max(n_d2_d1_name, len_trim(con(i)%d2_d1_name))
  if (con(i)%expression) then
    n_tot = max(n_tot,  len_trim(con(i)%name))
  else
    n_name = max(n_name, len_trim(con(i)%name))
  endif
  n_loc_ele   = max(n_loc_ele, len_trim(con(i)%loc_ele))
  n_loc_start = max(n_loc_start, len_trim(con(i)%loc_start))
  n_loc_ref   = max(n_loc_ref, len_trim(con(i)%loc_ref))
enddo

n_tot = max(n_tot,  n_name + n_loc_ref + n_loc_start + n_loc_ele + 6)
n_loc_ele = min(n_tot - n_name - n_loc_ref - n_loc_start - 6, len(con(1)%loc_ele))
n_tot = n_name + n_loc_ref + n_loc_start + n_loc_ele + 6

!

l1 = 'Constraints'
n=3+n_d2_d1_name+2+n_name;    l1(n:) = 'Ref_Ele'
n=len_trim(l1)+n_loc_ref-3;   l1(n:) = 'Start_Ele'
n=len_trim(l1)+n_loc_start-6; l1(n:) = 'Ele/S'
n=len_trim(l1)+n_loc_ele-3;   l1(n:) = 'Meas-or-Lim    Model       Merit      Max_Loc'

nl=nl+1; line(nl) = ' '
nl=nl+1; line(nl) = l1

!

fmt  = '(a, 4(2x, a), es13.4, es13.4, es11.3, 2x, a)'
fmt2 = '(a, 1(2x, a), es13.4, es13.4, es11.3, 2x, a)'

call re_allocate (line, nl+n_max+100)
do j = 1, n_max
  i = ixm(j)
  c => con(i)
  if (c%expression) then
    nl=nl+1; write (line(nl), fmt2) c%d2_d1_name(1:n_d2_d1_name), &
            c%name(1:n_tot), c%target_value, c%actual_value, c%merit, c%max_loc
  else
    nl=nl+1; write (line(nl), fmt) c%d2_d1_name(1:n_d2_d1_name), c%name(1:n_name), c%loc_ref(1:n_loc_ref), &
            c%loc_start(1:n_loc_start), c%loc_ele(1:n_loc_ele), c%target_value, c%actual_value, c%merit, c%max_loc
  endif
enddo
nl=nl+1; line(nl) = l1

!

nl=nl+1; line(nl) = ' '
nl=nl+1; write (line(nl), '(1x, a, es13.6)') 'figure of merit:', this_merit

call tao_write_lines (iunit, line(1:nl))

deallocate (con, ixm)

end subroutine tao_show_constraints

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_var_write (out_file, show_good_opt_only, tao_format)
!
! Routine to write the optimized variables. One file will be created for each universe.
! The created file will have three sections:
!   1) The variable values
!   2) The list of constraints.
!   3) A list of the top 10 constraints.
! If out_file = '' the information will be dumped to the terminal.
! In this case, only the variable values will be printed.
!
! When tao_format = True, the output is in the form "set variable <name> = <value>"
! so the file can be used as a Tao command file. If tao_format = False, the format
! is suitable for inclusion in a Bmad lattice file.
!
! Input:
!   out_file            -- character(*): Name of output file. If blank. Ouput to the terminal.
!   show_good_opt_only  -- logical, optional: Write only the variables used in the optimization?
!                           Default is False.
!   tao_format          -- logical, optional: Output format. Default False. See above.
!-

subroutine tao_var_write (out_file, show_good_opt_only, tao_format)

implicit none

integer i, j, iu, iu2, ix, ios, ix_hash

character(*) out_file
character(200) file_name, file_name2
character(20) :: r_name = 'tao_var_write'
character(200) str(1)

logical, optional :: show_good_opt_only, tao_format

! Output to terminal?

if (out_file == '') then
  call tao_print_vars (0, 0, show_good_opt_only, tao_format)
  return
endif

! Output to file

iu = lunget()
ix_hash = index (out_file, '#')

do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. s%u(i)%is_on .and. .not. logic_option(.false., tao_format)) cycle

  if (ix_hash /= 0) then
    write (file_name, '(a, i0, a)') out_file(1:ix_hash-1), i, trim(out_file(ix_hash+1:))
  else
    file_name = out_file
  endif

  open (iu, file = file_name, recl = 300, iostat = ios)
  if (ios /= 0) then
    call out_io (s_error$, r_name, 'ERROR IN VAR_WRITE: CANNOT OPEN FILE: ' // file_name)
    return
  endif

  call tao_print_vars (iu, i, show_good_opt_only, tao_format)
  call tao_write_lines (iu, ['        ', 'end_file', '        '])
  call tao_show_constraints (iu, 'TOP10')
  call tao_top10_merit_categories_print (iu)
  if (size(s%u) == 1) call tao_show_constraints (iu, 'ALL')

  close (iu)
  call out_io (s_blank$, r_name, 'Written: ' // file_name)

  if (logic_option(.false., tao_format)) return
enddo

! Write all constraints to a separate file when there are multiple universes.
! This can save a lot of memory when the number of universes is large.

if (size(s%u) > 1) then
  file_name = 'all_constraints.out'
  open (iu, file = file_name, recl = 300, iostat = ios)
  call tao_show_constraints (iu, 'ALL')
  close (iu)
  call out_io (s_blank$, r_name, 'Written constraints file: ' // file_name)
endif

end subroutine tao_var_write

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_print_vars (iu, ix_uni, show_good_opt_only, tao_format, v_array)
!
! Routine to print a list of set statements for the Bmad parameters controlled by the Tao variables.
! The set statements are Bmad lattice format compatible.
!
! When tao_format = True, the output is in the form "set variable <name> = <value>"
! so the file can be used as a Tao command file. If tao_format = False, the format
! is suitable for inclusion in a Bmad lattice file.
!
! Input:
!   iu                  -- integer: File unit number. 0 => print to the terminal.
!   ix_uni              -- integer: Universe index. If zero print slave parameters for all universes.
!                            If non-zero, only print set statements for slave parameters of this universe.
!                            Ignored if tao_format = True.
!   show_good_opt_only  -- logical: If True, only show slave parameters of variables used in optimization.
!   tao_format          -- logical, optional: Output format. Default False. See above.
!   v_array(:)          -- tao_var_array_struct, optional: Variable array. 
!                            If present, restrict printing to parameters of these variables.
!-

subroutine tao_print_vars (iu, ix_uni, show_good_opt_only, tao_format, v_array)

implicit none

type (tao_var_array_struct), optional :: v_array(:)
type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: lat

integer iu, ix_uni, j, n_line
character(40) useit_str
character(200), allocatable :: str(:)
logical, optional :: show_good_opt_only, tao_format

!

n_line = 0

if (present(v_array)) then
  allocate (str(size(v_array)))
  do j = 1, size(v_array)
    call print_this_var(v_array(j)%v, str, n_line, tao_format)
  enddo

else
  allocate (str(s%n_var_used))
  do j = 1, s%n_var_used
    call print_this_var(s%var(j), str, n_line, tao_format)
  enddo
endif

call tao_write_lines (iu, str(1:n_line))

!----------------------------------------------------------------------------------
contains
subroutine print_this_var (var, str, n_line, tao_format)

type (tao_var_struct), target :: var
type (ele_struct), pointer :: ele
type (tao_var_slave_struct), pointer :: slave(:)
integer n_line
integer ix, n_ele, i_uni
logical, optional :: tao_format
character(*), allocatable :: str(:)

!

slave => var%slave

if (.not. var%exists) return
if (ix_uni /= 0 .and. .not. any (slave(:)%ix_uni == ix_uni) .and. .not. logic_option(.false., tao_format)) return
if (logic_option(.false., show_good_opt_only) .and. .not. var%useit_opt) return
if (var%useit_opt) then
  useit_str = ''
else
  useit_str = '! Not used in optimizing'
endif

if (n_line+1 >= size(str)) call re_allocate (str, 2*(n_line+1))

! A potential problem is that the variable may not control all elements named var%ele_name in a universe.
! If this is the case, must use ele name order index to qualify the name.

! Tao "set" format.

if (logic_option(.false., tao_format)) then
  n_line=n_line+1; write (str(n_line), '(3a, i0, a, es22.15)') 'set var ', trim(var%v1%name), '[', var%ix_v1, ']|model =', var%model_value

! Controlling non-element parameter.

elseif (slave(1)%ix_ele == -1) then
  n_line=n_line+1; write (str(n_line), '(4a, es25.17e3, 3x, a)')  trim(var%ele_name), '[', trim(var%attrib_name), '] = ', var%model_value, useit_str

! Universe independent.

elseif (ix_uni == 0) then
  i_uni = slave(1)%ix_uni
  u => s%u(i_uni)
  lat => u%model%lat
  ele => lat%branch(slave(1)%ix_branch)%ele(slave(1)%ix_ele)
  n_line=n_line+1; write (str(n_line), '(4a, es25.17e3, 3x, a)')  trim(ele_unique_name(ele, u%ele_order)), &
                                                            '[', trim(var%attrib_name), '] = ', var%model_value, useit_str

! Universe is given.

else
  u => s%u(ix_uni)
  lat => u%model%lat

  ix = nametable_bracket_indexx(lat%nametable, var%ele_name, n_ele)
  if (n_ele == count(slave%ix_uni == ix_uni)) then
    n_line=n_line+1; write (str(n_line), '(4a, es25.17e3, 3x, a)')  trim(var%ele_name), '[', trim(var%attrib_name), '] = ', var%model_value, useit_str
  else
    if (n_line+size(slave)+1 >= size(str)) call re_allocate (str, 2*(n_line+size(slave)))

    ! Qualified names can only be used after an expand_lattice command.

    if (str(1) /= 'expand_lattice') then
      str(2:n_line+1) = str(1:n_line)
      n_line = n_line + 1
      str(1) = 'expand_lattice'
    endif

    do ix = 1, size(slave)
      if (slave(ix)%ix_uni /= ix_uni) cycle
      ele => lat%branch(slave(ix)%ix_branch)%ele(slave(ix)%ix_ele)
      n_line=n_line+1; write (str(n_line), '(4a, es25.17e3, 3x, a)')  trim(ele_unique_name(ele, u%ele_order)), &
                                                              '[', trim(var%attrib_name), '] = ', var%model_value, useit_str
    enddo
  endif
endif

end subroutine print_this_var

end subroutine tao_print_vars

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_write_lines (iunit, line)
!
! Subroutine to write out a series of lines to a file or to the terminal.
! It is assumed that any file has already been opened.
!
! Input:
!   iunit   -- Integer: File unit to write to. 0 => print to the terminal.
!   line(:) -- Character(*): A series of lines.
!-

subroutine tao_write_lines (iunit, line)

implicit none

integer iunit, i
character(*) line(:)

!

if (iunit == 0) then
  call out_io (s_blank$, ' ', line)
else
 do i = 1, size(line)
    write (iunit, *) trim(line(i))
  enddo
endif

end subroutine tao_write_lines

end module
