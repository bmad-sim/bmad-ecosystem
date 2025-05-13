module write_lattice_file_mod

use element_modeling_mod
use binary_parser_mod

type multipass_region_ele_struct
  integer ix_region
  logical region_start_pt
  logical region_stop_pt
end type

type multipass_region_branch_struct
  type (multipass_region_ele_struct), allocatable :: ele(:)
end type

type multipass_region_lat_struct
  type (multipass_region_branch_struct), allocatable :: branch(:)
end type

logical, save :: write_lat_debug_flag = .false.

interface re_str
  module procedure re_str_rp
  module procedure re_str_qp
end interface

contains

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
! Create the information on multipass regions.

subroutine multipass_region_info(lat, mult_lat, m_info)

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (multipass_region_lat_struct), target :: mult_lat
type (multipass_all_info_struct), target :: m_info

type (multipass_region_ele_struct), pointer :: mult_ele(:), m_ele
type (multipass_ele_info_struct), pointer :: e_info
type (ele_pointer_struct), pointer :: ss1(:), ss2(:)

integer ib, ix_r, ie, ix_pass, ix_lord, ix_super
logical in_multi_region, need_new_region

!

allocate (mult_lat%branch(0:ubound(lat%branch, 1)))
do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  allocate (mult_lat%branch(ib)%ele(0:branch%n_ele_max))
  mult_lat%branch(ib)%ele(:)%ix_region = 0
  mult_lat%branch(ib)%ele(:)%region_start_pt = .false.
  mult_lat%branch(ib)%ele(:)%region_stop_pt   = .false.
enddo

call multipass_all_info (lat, m_info)

if (size(m_info%lord) == 0) return

! Go through and mark all 1st pass regions
! In theory the original lattice file could have something like:
!   lat: line = (..., m1, m2, ..., m1, -m2, ...)
! where m1 and m2 are multipass lines. The first pass region (m1, m2) looks 
! like this is one big region but the later (m1, -m2) signals that this 
! is not so.
! We thus go through all the first pass regions and compare them to the
! corresponding higher pass regions. If we find two elements that are contiguous
! in the first pass region but not contiguous in some higher pass region, 
! we need to break the first pass region into two.

ix_r = 0
do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  mult_ele => mult_lat%branch(ib)%ele

  in_multi_region = .false.

  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    e_info => m_info%branch(ib)%ele(ie)
    ix_pass = e_info%ix_pass

    if (ix_pass /= 1) then  ! Not a first pass region
      if (in_multi_region) mult_ele(ie-1)%region_stop_pt = .true.
      in_multi_region = .false.
      cycle
    endif

    ! If start of a new region...
    if (.not. in_multi_region) then  
      ix_r = ix_r + 1
      mult_ele(ie)%ix_region = ix_r
      mult_ele(ie)%region_start_pt = .true.
      in_multi_region = .true.
      ix_lord = e_info%ix_lord(1)
      ix_super = e_info%ix_super(1)
      ss1 => m_info%lord(ix_lord)%slave(:,ix_super)
      cycle
    endif
    ix_lord = e_info%ix_lord(1)
    ix_super = e_info%ix_super(1)
    ss2 => m_info%lord(ix_lord)%slave(:, ix_super)

    need_new_region = .false.
    if (size(ss1) /= size(ss2)) then
      need_new_region = .true.
    else
      do ix_pass = 2, size(ss1)
        if (abs(ss1(ix_pass)%ele%ix_ele - ss2(ix_pass)%ele%ix_ele) == 1) cycle
        ! not contiguous then need a new region
        need_new_region = .true.
        exit
      enddo
    endif

    if (need_new_region) then
      ix_r = ix_r + 1
      mult_ele(ie-1)%region_stop_pt = .true.
      mult_ele(ie)%region_start_pt = .true.
    endif

    ss1 => ss2
    mult_ele(ie)%ix_region = ix_r
  enddo

enddo

if (in_multi_region) mult_ele(branch%n_ele_track)%region_stop_pt = .true.

end subroutine multipass_region_info

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

subroutine write_line_element (line, iu, ele, lat)

implicit none

type (lat_struct), target :: lat
type (ele_struct) :: ele
type (ele_struct), pointer :: lord, m_lord, slave

character(*) line
character(40) lord_name

integer iu, ix

!

if (ele%slave_status == super_slave$) then
  if (ele%orientation == 1) then
    write (line, '(a, 2(a, i0), a)') trim(line), ' slave_drift_', ele%ix_branch, '_', ele%ix_ele, ','
  else
    write (line, '(a, 2(a, i0), a)') trim(line), ' --slave_drift_', ele%ix_branch, '_', ele%ix_ele, ','
  endif

elseif (ele%slave_status == multipass_slave$) then
  lord => pointer_to_lord(ele, 1)
  write (line, '(4a)') trim(line), ' ', trim(lord%name), ','

else
  if (ele%orientation == 1) then
    write (line, '(4a)') trim(line), ' ', trim(ele%name), ','
  else
    write (line, '(4a)') trim(line), ' --', trim(ele%name), ','
  endif
endif

if (len_trim(line) > 100) call write_lat_line(line, iu, .false.)

end subroutine write_line_element

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

function re_str_rp(rel) result (str_out)

implicit none

real(rp) rel
character(25) :: str_out

!

if (write_lat_debug_flag) then
  write(str_out, '(es13.5)') rel
else
  write(str_out, '(es25.17e3)') rel
endif

end function re_str_rp

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

function re_str_qp(rel) result (str_out)

implicit none

real(qp) rel
character(25) :: str_out

!

if (write_lat_debug_flag) then
  write(str_out, '(es13.5)') rel
else
  write(str_out, '(es25.17e3)') rel
endif

end function re_str_qp

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

function array_re_str(arr, parens_in) result (str_out)

real(rp) arr(:)
integer i
character(120) str_out
character(*), optional :: parens_in
character(2) parens

!

parens = '()'
if (present(parens_in)) parens = parens_in

str_out = parens(1:1) // re_str(arr(1))
do i = 2, size(arr)
  str_out = trim(str_out) // ', ' // re_str(arr(i))
enddo
str_out = trim(str_out) // parens(2:2)

end function array_re_str

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

function cmplx_re_str(cmp) result (str_out)

complex(rp) cmp
character(40) str_out

!

if (imag(cmp) == 0) then
  str_out = re_str(real(cmp))
else
  str_out = '(' // re_str(real(cmp)) // ', ' // re_str(imag(cmp)) // ')'
endif

end function cmplx_re_str

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------

function rchomp (rel, plc) result (out)

implicit none

real(rp) rel
character(25) out
character(8) :: fmt = '(f24.xx)'
integer it, plc, ix

! The output when running with debug has less precision to prevent slight shifts in numbers (which happens
! when the translation code is run with different compilers) from changing the output. Output produced
! when running debug is used in regression testing.

if (write_lat_debug_flag) then
  write (fmt(6:7), '(i2.2)') 5-plc  ! 6 digits of accuracy
else
  write (fmt(6:7), '(i2.2)') 14-plc  ! 15 digits of accuracy
endif

write (out, fmt) rel
do it = len(out), 1, -1
  if (out(it:it) == ' ') cycle
  if (out(it:it) == '0') then
    out(it:it) = ' '
    cycle
  endif
  if (out(it:it) == '.') out(it:it) = ' '
  call string_trim(out, out, ix)
  return
enddo

end function rchomp

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
!+
! Subroutine write_lat_line (line, iu, end_is_neigh, do_split)
!
! Routine to write strings to a lattice file.
! This routine will break the string up into multiple lines
! if the string is too long and add a continuation character if needed.
!
! If the "line" arg does not represent a full "sentence" (end_is_neigh = False), 
! then only part of the line may be written and the part not written will be returned.
!
! Input:
!   line          -- character(*): String of text.
!   iu            -- integer: Unit number to write to.
!   end_is_neigh  -- logical: If true then write out everything.
!                      Otherwise wait for a full line of max_char characters or so.
!   do_split      -- logical, optional: Split line if overlength? Default is True.
!                      False is used when line has already been split for expressions since
!                      the expression splitting routine does a much better job of it.
!   julia         -- logical, optional: Default False. If True then do not include "&" line continuation
!
! Output:
!   line          -- character(*): part of the string not written. 
!                       If end_is_neigh = T then line will be blank.
!-

subroutine write_lat_line (line, iu, end_is_neigh, do_split, julia)

implicit none

character(*) line
integer i, iu, n
integer, parameter :: max_char = 105
logical end_is_neigh
logical, save :: init = .true.
logical, optional :: do_split, julia

!

if (.not. logic_option(.true., do_split)) then
  n = len_trim(line)
  if (end_is_neigh) then
    call write_this (line)
    init = .true.
  elseif (index(',[{(=', line(n:n)) /= 0) then
    call write_this (line)
  else
    if (logic_option(.false., julia)) then
      call write_this (trim(line))
    else
      call write_this (trim(line) // ' &')
    endif
  endif

  line = ''
  return
endif

!

outer_loop: do 

  if (len_trim(line) <= max_char) then
    if (end_is_neigh) then
      call write_this (line)
      line = ''
      init = .true.
    endif
    return
  endif

  i = index(line(1:max_char), ',', back = .true.)
  if (i /= 0) then
    call write_this (line(:i))
    line = line(i+1:)
    cycle outer_loop
  endif

  i = index(line, ',', back = .true.)
  if (i /= 0) then
    call write_this (line(:i))
    line = line(i+1:)
    cycle outer_loop
  endif

  if (end_is_neigh) then
    call write_this (line)
    init = .true.
    return
  endif

  if (logic_option(.false., julia)) then
    call write_this (trim(line))
  else
    call write_this (trim(line) // ' &')
  endif
  line = ''
  return

enddo outer_loop

!-----------------------------------

contains

subroutine write_this (line2)

character(*) line2
character(20) fmt

!

if (init) then
  fmt = '(a, 1x, a)'
  init = .false.
else
  fmt = '(2x, a, 1x, a)'
endif

write (iu, fmt) trim(line2)

end subroutine write_this

end subroutine write_lat_line

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

subroutine value_to_line (line, value, str, typ, ignore_if_zero, use_comma)

use precision_def

implicit none

character(*) line, str
character(40) fmt, val_str
character(*) typ

real(rp) value

integer ix

logical, optional :: ignore_if_zero, use_comma

!

if (value == 0 .and. logic_option(.true., ignore_if_zero)) return

if (logic_option(.true., use_comma)) then
  if (str == '') then
    line = trim(line) // ','
  else
    line = trim(line) // ', ' // trim(str) // ' ='
  endif
else
  if (str /= '') then
    line = trim(line) // ' ' // trim(str) // ' ='
  endif
endif

if (value == 0) then
  line = trim(line) // ' 0'
  return
endif

if (typ == 'R') then
  val_str = re_str(value)
elseif (typ == 'I') then
  write (val_str, '(i0)') nint(value)
else
  print *, 'ERROR IN VALUE_TO_LINE. BAD "TYP": ', typ 
  if (global_com%exit_on_error) call err_exit
endif

call string_trim(val_str, val_str, ix)
line = trim(line) // ' ' // trim(val_str)

end subroutine value_to_line

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

subroutine add_this_name_to_list (ele, names, an_indexx, n_names, ix_match, has_been_added, named_eles)

type (ele_struct), target :: ele
type (ele_pointer_struct), allocatable :: named_eles(:)  ! List of unique element names 

integer, allocatable :: an_indexx(:)
integer n_names, ix_match
logical has_been_added
character(40), allocatable :: names(:)

!

if (size(names) < n_names + 1) then
  call re_allocate(names, 2*size(names))
  call re_allocate(an_indexx, 2*size(names))
  call re_allocate_eles(named_eles, 2*size(names), .true.)
endif
call find_index (ele%name, names, an_indexx, n_names, ix_match, add_to_list = .true., has_been_added = has_been_added)
if (has_been_added) named_eles(n_names)%ele => ele

end subroutine add_this_name_to_list

end module

