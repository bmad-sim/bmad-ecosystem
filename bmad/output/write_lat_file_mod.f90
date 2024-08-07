module write_lat_file_mod

use element_modeling_mod
use binary_parser_mod
use multipole_mod

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

function re_str(rel) result (str_out)

implicit none

real(rp) rel
integer pl, n
character(:), allocatable :: str_out
character(24) str
character(16) fmt

!

if (rel == 0) then
  allocate(character(1):: str_out)
  str_out = '0'
  return
endif

pl = floor(log10(abs(rel)))

if (pl > 5) then
  fmt = '(2a, i0)'
  write (str, fmt) trim(rchomp(rel/10.0_rp**pl, 0)), 'E', pl

elseif (pl > -3) then
  str = rchomp(rel, pl)

else
  fmt = '(2a, i0)'
  write (str, fmt) trim(rchomp(rel*10.0_rp**(-pl), 0)), 'E', pl
endif

n = len_trim(str)
allocate(character(n):: str_out)
str_out = str(1:n)


end function re_str

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

!

write (fmt(6:7), '(i2.2)') 14-plc  ! 15 digits of accuracy
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

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+ 
! Subroutine write_lat_in_sad_format (out_file_name, lat, include_apertures, ix_start, ix_end, ix_branch, converted_lat, err)
!
! Private routine used by write_lat_in_sad_format and not for general use. 
! See write_lat_in_sad_format for details about the argument list.
!-

subroutine write_lat_in_sad_format (out_file_name, lat, include_apertures, ix_start, ix_end, ix_branch, converted_lat, err)

implicit none

type (lat_struct), target :: lat, lat_model, lat_out
type (lat_struct), optional, target :: converted_lat
type (ele_struct), pointer :: ele, ele1, ele2, lord, sol_ele, first_sol_edge
type (branch_struct), pointer :: branch, branch_out
type (ele_struct), save :: drift_ele, ab_ele, taylor_ele, col_ele, kicker_ele, null_ele, bend_ele, quad_ele

real(rp), pointer :: val(:)
real(rp) knl(0:n_pole_maxx), tilts(0:n_pole_maxx), a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
real(rp) bs_field, old_bs_field, a, b
real(rp) xp, yp, xo, yo, zo, sad_fshift
real(rp) hk, vk, tilt

! These will be used to point to slots in ele%old_value(:) array for extra needed information in converting to SAD.
integer, parameter :: sad_f1$          = 1
integer, parameter :: sad_geo$         = 2
integer, parameter :: sad_bound$       = 3
integer, parameter :: sad_bz$          = 4
integer, parameter :: sad_fshift$      = 5
integer, parameter :: sad_mark_offset$ = 6

integer, optional :: ix_start, ix_end, ix_branch
integer, allocatable :: n_repeat(:), an_indexx(:)
integer i, j, n, ib, iout, iu, ix, ix1, ix2, ios, ie1, ie2, ie2_orig, n_taylor_order_saved, ix_ele
integer s_count, j_count, n_elsep_warn, n_name_change_warn
integer ix_manch, n_names, aperture_at, ix_pole_max, ix_match
integer :: ix_line_min, ix_line_max, n_warn_max = 10

character(*), parameter :: r_name = "write_lat_in_sad_format"
character(*) out_file_name
character(40), allocatable :: names(:)
character(40) str, orig_name
character(300) line, knl_str, ksl_str
character(2000) line_out

logical, optional :: include_apertures, err
logical converted, init_needed, in_solenoid, ptc_exact_model

! Use ptc exact_model = True since this is needed to get the drift nonlinear terms

ptc_exact_model = ptc_com%exact_model
ptc_com%exact_model = .true.

! Init

ix = integer_option(0, ix_branch)
if (ix < 0 .or. ix > ubound(lat%branch, 1)) then
  call out_io (s_error$, r_name, 'BRANCH INDEX OUT OF RANGE: /i0/ ', i_array = [ix])
  return
endif

branch => lat%branch(ix)

ix_line_max = 120
ix_line_min = ix_line_max - 20

call init_ele (col_ele)
call init_ele (drift_ele, drift$)
call init_ele (ab_ele, ab_multipole$)
call init_ele (kicker_ele, kicker$) 
call init_ele (quad_ele, quadrupole$)
call init_ele (bend_ele, sbend$)
call multipole_init (ab_ele, magnetic$)
null_ele%key = null_ele$

ie1 = integer_option(1, ix_start)
ie2 = integer_option(branch%n_ele_track, ix_end)
ie2_orig = ie2

allocate (names(branch%n_ele_max), an_indexx(branch%n_ele_max)) ! list of element names

call out_io (s_info$, r_name, &
      'Note: Bmad lattice elements have attributes that cannot be translated. ', &
      '      For example, higher order terms in a Taylor element.', &
      '      Please use caution when using a translated lattice.')

! open file

if (present(err)) err = .true.
n_taylor_order_saved = ptc_private%taylor_order_ptc

iu = lunget()
call fullfilename (out_file_name, line)
open (iu, file = line, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(out_file_name))
  return
endif

!-----------------------------------------------------------------------------
! Translation is a two step process:
!   1) Create a new lattice called lat_out making substitutions for sol_quad and wiggler elements, etc..
!   2) Use lat_out to create the lattice file.

lat_out = lat
call allocate_lat_ele_array(lat_out, 2*branch%n_ele_max, branch%ix_branch)
branch_out => lat_out%branch(branch%ix_branch)

j_count = 0    ! drift around solenoid or sol_quad index. Also z shift count.
s_count = 0    ! SAD solenoid count
sad_fshift = 0
in_solenoid = .false.

! Remove all patch elements that are just time shifts.

do ix_ele = 1, branch_out%n_ele_track
  ele => branch_out%ele(ix_ele)
  ele%old_value(sad_f1$)     = value_of_attribute(ele, 'SAD_F1')
  ele%old_value(sad_geo$)    = value_of_attribute(ele, 'SAD_GEO')
  ele%old_value(sad_bound$)  = value_of_attribute(ele, 'SAD_BOUND')
  ele%old_value(sad_bz$)     = value_of_attribute(ele, 'SAD_BZ')
  ele%old_value(sad_fshift$) = value_of_attribute(ele, 'SAD_FSHIFT')

  if (ele%key == patch$ .and. ele%old_value(sad_fshift$) /= 0 .and. ix_ele <= ie2_orig) then
    ele%ix_ele = -1
    ie2 = ie2 - 1
  endif
enddo

call remove_eles_from_lat (lat_out)

! Loop over all input elements

branch_out%ele%old_value(sad_mark_offset$) = 0  ! SAD mark offset
nullify(first_sol_edge)
old_bs_field = 0
n_name_change_warn = 0
n_elsep_warn = 0
ix_ele = ie1 - 1

do
  ix_ele = ix_ele + 1
  if (ix_ele > ie2) exit
  ele => branch_out%ele(ix_ele)
  val => ele%value

  if (ele%sub_key == not_set$) cycle

  ! Must create SAD sol elements as needed if Bmad lattice was not derived from a SAD lattice.
  ! If there is a marker or patch element that can be used, convert it.
  ! Otherwise, insert a new element. 
  ! Bmad does not have a sol element so use null_ele to designate the sol element in the Bmad lattice.
  ! This works since there cannot be any actual null_eles in the lattice.

  if ((ele%key == patch$ .or. ele%key == marker$) .and. ele%old_value(sad_bound$) /= 0) then 
    ele%key = null_ele$  ! Convert to SOL
    if (in_solenoid) then
      ele%iyy = exit_end$
      bs_field = 0
      old_bs_field = 0
    else
      ele%iyy = entrance_end$
    endif
    in_solenoid = (.not. in_solenoid)

  elseif (ele%key == patch$ .or. ele%value(l$) /= 0 .or. ix_ele == branch_out%n_ele_max) then
    bs_field = 0
    if (has_attribute (ele, 'BS_FIELD')) bs_field = ele%value(bs_field$)

    if (bs_field /= old_bs_field) then
      if (ele%key == marker$ .or. ele%key == patch$) then
        sol_ele => ele
      else
        sol_ele => pointer_to_next_ele(ele, -1)
        ! Look to see if there is a marker or patch element that can be converted to a SAD SOL.
        do
          if (sol_ele%ix_ele == 1) exit
          if (sol_ele%key == marker$ .or. sol_ele%key == patch$ .or. sol_ele%key == null_ele$) exit
          if (sol_ele%value(l$) /= 0) exit  ! No suitable marker or patch found
          sol_ele => pointer_to_next_ele(sol_ele, -1)
        enddo
      endif

      if (sol_ele%key /= marker$ .and. sol_ele%key /= patch$ .and. sol_ele%key /= null_ele$) then
        s_count = s_count + 1
        write (sol_ele%name, '(a, i0)') 'SOL_', s_count  
        call insert_element (lat_out, sol_ele, ix_ele, branch_out%ix_branch)
        sol_ele => branch_out%ele(ix_ele)
        ie2 = ie2 + 1
        ix_ele = ix_ele + 1
        ele => branch_out%ele(ix_ele)
        val => ele%value
      endif

      if (sol_ele%key /= null_ele$) then
        if (old_bs_field == 0) then
          sol_ele%old_value(sad_geo$) = 1
          sol_ele%old_value(sad_bound$) = 1
          sol_ele%iyy = entrance_end$  ! Entering solenoid
          in_solenoid = .true.
          first_sol_edge => sol_ele
        elseif (bs_field == 0 .and. in_solenoid) then
          if (nint(ele%old_value(sad_geo$)) == 1 .and. nint(ele%old_value(sad_bound$)) == 1) then
            sol_ele%old_value(sad_geo$) = 1
            sol_ele%old_value(sad_bound$) = 1
            first_sol_edge%old_value(sad_bound$) = 1
            first_sol_edge%old_value(sad_geo$) = 0
          else
            sol_ele%old_value(sad_bound$) = 1
            sol_ele%old_value(sad_geo$) = 0
          endif
          sol_ele%iyy = exit_end$  ! Entering solenoid
          in_solenoid = .false.
        else
          sol_ele%old_value(sad_geo$) = 0
          sol_ele%old_value(sad_bound$) = 0
        endif
      endif

      sol_ele%key = null_ele$
      if (bs_field /= 0) sol_ele%old_value(sad_bz$) = bs_field
      old_bs_field = bs_field
    endif
  endif

  ! With an element superimposed with a marker create a whole element
  ! plus a marker with an offset

  if (ele%slave_status == super_slave$) then
    lord => pointer_to_lord(ele, 1, ix_slave_back = ix)
    ele1 => pointer_to_next_ele(ele)
    ele2 => pointer_to_slave(lord, lord%n_slave)
    if (lord%n_slave == 2 .and. ele1%key == marker$ .and. &
          num_lords(ele, super_lord$) == 1 .and. num_lords(ele2, super_lord$) == 1) then
      ele1%old_value(sad_mark_offset$) = -ele2%value(l$)/lord%value(l$) ! marker offset
      ele = lord                              ! Super_slave piece becomes the entire element.
      ele2%sub_key = not_set$                 ! Ignore this super_slave piece.
    endif
  endif

  ! Replace element name containing "/" or "#" with "_"

  orig_name = ele%name

  do
    j = max(index(ele%name, '#'), index(ele%name, '\'))       ! '
    if (j == 0) exit
    ele%name(j:j) = '_'
  enddo

  if (ele%name /= orig_name .and. n_name_change_warn <= n_warn_max) then
    call out_io (s_info$, r_name, 'Element name changed from: ' // trim(orig_name) // ' to: ' // ele%name)
    if (n_name_change_warn == n_warn_max) call out_io (s_info$, r_name, &
                           'Enough name change warnings. Will stop issuing them now.')
    n_name_change_warn = n_name_change_warn + 1
  endif

  ! SAD: If there is an aperture with an element that is not an ecoll or rcoll then need to make a separate
  ! element with the aperture info. 

  if ((val(x1_limit$) /= 0 .or. val(x2_limit$) /= 0 .or. val(y1_limit$) /= 0 .or. val(y2_limit$) /= 0) .and. &
      .not. (((ele%key == ecollimator$ .or. ele%key == rcollimator$) .and. ele%value(l$) == 0) .or. ele%key == sad_mult$) .and. &
      logic_option(.true., include_apertures)) then

    if (val(x1_limit$) /= val(x2_limit$)) then
      call out_io (s_warn$, r_name, 'Asymmetric x_limits cannot be converted for: ' // ele%name, &
                                    'Will use largest limit here.')
      val(x1_limit$) = max(val(x1_limit$), val(x2_limit$))
    endif

    if (val(y1_limit$) /= val(y2_limit$)) then
      call out_io (s_warn$, r_name, 'Asymmetric y_limits cannot be converted for: ' // ele%name, &
                                    'Will use largest limit here.')
      val(y1_limit$) = max(val(y1_limit$), val(y2_limit$))
    endif

    ! Create ecoll and rcoll elements.
    ! If original element is itself a collimator, turn it into a drift.

    if (ele%aperture_type == rectangular$) then
      col_ele%key = rcollimator$
    else
      col_ele%key = ecollimator$
    endif

    if (ele%key == ecollimator$ .or. ele%key == rcollimator$) then
      col_ele%name = ele%name
      ele%key = drift$
      ele%name = 'DRIFT_' // trim(ele%name)
    else
      col_ele%name = 'COLLIMATOR_' // trim(ele%name)
    endif

    col_ele%value = val
    col_ele%value(l$) = 0
    val(x1_limit$) = 0; val(x2_limit$) = 0; val(y1_limit$) = 0; val(y2_limit$) = 0; 
    aperture_at = ele%aperture_at  ! Save since ele pointer will be invalid after the insert
    if (aperture_at == both_ends$ .or. aperture_at == downstream_end$ .or. aperture_at == continuous$) then
      call insert_element (lat_out, col_ele, ix_ele+1, branch_out%ix_branch)
      ie2 = ie2 + 1
    endif
    if (aperture_at == both_ends$ .or. aperture_at == upstream_end$ .or. aperture_at == continuous$) then
      call insert_element (lat_out, col_ele, ix_ele, branch_out%ix_branch)
      ie2 = ie2 + 1
    endif
    ix_ele = ix_ele - 1 ! Want to process the element again on the next loop.

    cycle ! cycle since ele pointer is invalid
  endif

  ! If the bend has a roll then put kicker elements just before and just after

  if ((ele%key == sbend$) .and. val(roll$) /= 0) then
    j_count = j_count + 1
    write (kicker_ele%name,   '(a, i0)') 'ROLL_Z', j_count
    kicker_ele%value(hkick$) =  val(angle$) * (1 - cos(val(roll$))) / 2
    kicker_ele%value(vkick$) = -val(angle$) * sin(val(roll$)) / 2
    val(roll$) = 0   ! So on next iteration will not create extra kickers.
    call insert_element (lat_out, kicker_ele, ix_ele, branch_out%ix_branch)
    call insert_element (lat_out, kicker_ele, ix_ele+2, branch_out%ix_branch)
    ie2 = ie2 + 2
    cycle
  endif

  ! If there is a multipole component then put multipole elements at half strength 
  ! just before and just after the element.

  if (ele%key /= multipole$ .and. ele%key /= ab_multipole$ .and. ele%key /= null_ele$ .and. ele%key /= sad_mult$) then
    call multipole_ele_to_ab (ele, .true., ix_pole_max, ab_ele%a_pole, ab_ele%b_pole)
    if (ix_pole_max > -1) then
      ab_ele%a_pole = ab_ele%a_pole / 2
      ab_ele%b_pole = ab_ele%b_pole / 2
      if (associated(ele%a_pole)) deallocate (ele%a_pole, ele%b_pole)
      j_count = j_count + 1
      write (ab_ele%name, '(a1, a, i0)') key_name(ele%key), 'MULTIPOLE_', j_count
      call insert_element (lat_out, ab_ele, ix_ele, branch_out%ix_branch)
      call insert_element (lat_out, ab_ele, ix_ele+2, branch_out%ix_branch)
      ie2 = ie2 + 2
      cycle
    endif
  endif

  ! Convert wiggler elements to an "equivalent" set of elements.
  ! If the wiggler has been sliced due to superposition, throw 
  ! out the markers that caused the slicing.


  if (ele%key == wiggler$ .or. ele%key == undulator$) then
    call out_io (s_warn$, r_name, 'Converting element to drift-bend-drift model: ' // ele%name)
    if (ele%slave_status == super_slave$) then
      ! Create the wiggler model using the super_lord
      lord => pointer_to_lord(ele, 1)
      call create_planar_wiggler_model (lord, lat_model)
      ! Remove all the slave elements and markers in between.
      call out_io (s_warn$, r_name, 'Note: Not translating the markers within wiggler: ' // lord%name)
      lord%ix_ele = -1 ! mark for deletion
      call find_element_ends (lord, ele1, ele2)
      ix1 = ele1%ix_ele; ix2 = ele2%ix_ele
      ! If the wiggler wraps around the origin we are in trouble.
      if (ix2 < ix1) then 
        call out_io (s_fatal$, r_name, 'Wiggler wraps around origin. Cannot translate this!')
        if (global_com%exit_on_error) call err_exit
      endif
      do i = ix1+1, ix2
        branch_out%ele(i)%ix_ele = -1  ! mark for deletion
      enddo
      ie2 = ie2 - (ix2 - ix1 - 1)
    else
      call create_planar_wiggler_model (ele, lat_model)
    endif

    ele%ix_ele = -1 ! Mark for deletion
    call remove_eles_from_lat (lat_out)
    do j = 1, lat_model%n_ele_track
      call insert_element (lat_out, lat_model%ele(j), ix_ele+j-1, branch_out%ix_branch)
    enddo
    ie2 = ie2 + lat_model%n_ele_track - 1
    cycle
  endif

enddo

! If there is a finite bs_field then create a final null_ele element

if (bs_field /= 0) then
  ele => branch_out%ele(ie2)
  if (ele%key == marker$) then
    ele%key = null_ele$
    ele%old_value(sad_bz$) = bs_field
  else
    s_count = s_count + 1
    write (null_ele%name, '(a, i0)') 'SOL_', s_count  
    null_ele%old_value(sad_bz$) = bs_field
    ie2 = ie2 + 1
    call insert_element (lat_out, null_ele, ie2, branch_out%ix_branch)
    ele => branch_out%ele(ie2)
  endif

  ele%old_value(sad_bound$) = 1
  ele%old_value(sad_geo$) = 0
endif

! For a patch that is *not* associated with the edge of a solenoid: A z_offset must be split into a drift + patch

ix_ele = ie1 - 1

do
  ix_ele = ix_ele + 1
  if (ix_ele > ie2) exit
  ele => branch_out%ele(ix_ele)
  val => ele%value

  if (ele%key == patch$ .and. ele%value(z_offset$) /= 0) then
    drift_ele%name = 'DRIFT_' // ele%name
    drift_ele%value(l$) = val(z_offset$)
    call insert_element (lat_out, drift_ele, ix_ele, branch_out%ix_branch)
    ix_ele = ix_ele + 1
    ele => branch_out%ele(ix_ele)
    val => ele%value
    val(z_offset$) = 0
  endif
enddo

!-------------------------------------------------------------------------------------------------
! Now write info to the output file...
! lat lattice name

write (iu, '(3a)') '! File generated by Bmad from Bmad lattice file:'
write (iu, '(4x, 2a)') trim(lat%input_file_name), ';'
if (lat%lattice /= '') write (iu, '(4a)') '! Bmad lattice name: ', trim(lat%lattice), ';'
write (iu, '(a)')

write (iu, '(3a)') 'MOMENTUM = ',  re_str(ele%value(p0c$)), ';'
if (sad_fshift /= 0) write (iu, '(3a)') 'FSHIFT = ', re_str(sad_fshift), ';'

! write element parameters

n_names = 0                          ! number of names stored in the list
old_bs_field = 0

do ix_ele = ie1, ie2

  ele => branch_out%ele(ix_ele)
  val => ele%value

  if (ele%sub_key == not_set$) cycle 

  ! do not make duplicate specs

  call find_index (ele%name, names, an_indexx, n_names, ix_match)
  if (ix_match > 0) cycle

  ! Add to the list of elements

  if (size(names) < n_names + 1) then
    call re_allocate(names, 2*size(names))
    call re_allocate(an_indexx, 2*size(names))
  endif

  call find_index (ele%name, names, an_indexx, n_names, ix_match, add_to_list = .true.)

  converted = .false.
  if (.not. associated (ele%a_pole) .and. ele%value(hkick$) == 0 .and. ele%value(vkick$) == 0) then
    select case (ele%key)

    !
    case (octupole$)
      write (line_out, '(4a)') 'OCT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, val(k3$)*val(l$), 'K3', 'R', .true., .false.)
      converted = .true.

    !
    case (quadrupole$)
      write (line_out, '(4a)') 'QUAD ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, val(k1$)*val(l$), 'K1', 'R', .true., .false.)
      call value_to_line (line_out, -sign_of(val(fq1$)) * sqrt(24*abs(val(fq1$))), 'F1', 'R', .true., .false.)
      call value_to_line (line_out, val(fq2$), 'F2', 'R', .true., .false.)
      converted = .true.

    !
    case (sextupole$)
      write (line_out, '(4a)') 'SEXT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, val(k2$)*val(l$), 'K2', 'R', .true., .false.)
      converted = .true.
    end select
  endif

  ! If not yet converted due to the presence of multipoles or a kick

  if (.not. converted) then

    a_pole = 0; b_pole = 0
    if (ele%key /= null_ele$) call multipole_ele_to_ab (ele, .false., ix_pole_max, a_pole, b_pole)

    select case (ele%key)

    !
    case (drift$, pipe$)
      write (line_out, '(4a)') 'DRIFT ', trim(ele%name), ' = (L = ', re_str(val(l$))

    !
    case (instrument$, detector$, monitor$)
      if (ele%value(l$) == 0) then
        write (line_out, '(4a)') 'MONI ', trim(ele%name), ' = ('
      else
        write (line_out, '(4a)') 'DRIFT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      endif

    !
    case (ab_multipole$, multipole$)
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = ('
      do i = 0, ubound(a_pole, 1)
        write (str, '(a, i0)') 'K', i
        call value_to_line (line_out, a_pole(i) * factorial(i), 'S' // str, 'R', .true., .false.)
        call value_to_line (line_out, b_pole(i) * factorial(i), str, 'R', .true., .false.)
      enddo

    !
    case (ecollimator$)
      write (line_out, '(4a)') 'APERT ', trim(ele%name), ' = ('
      call value_to_line (line_out, val(x_offset$), 'DX', 'R', .true., .false.)
      call value_to_line (line_out, val(y_offset$), 'DY', 'R', .true., .false.)
      call value_to_line (line_out, val(x1_limit$), 'AX', 'R', .true., .false.)
      call value_to_line (line_out, val(y1_limit$), 'AY', 'R', .true., .false.)
      
    !
    case (rcollimator$)
      write (line_out, '(4a)') 'APERT ', trim(ele%name), ' = ('
      call value_to_line (line_out, -val(x1_limit$), 'DX1', 'R', .true., .false.)
      call value_to_line (line_out, -val(y1_limit$), 'DY1', 'R', .true., .false.)
      call value_to_line (line_out,  val(x2_limit$), 'DX2', 'R', .true., .false.)
      call value_to_line (line_out,  val(y2_limit$), 'DY2', 'R', .true., .false.)

    !
    case (elseparator$)
      call out_io (s_warn$, r_name, 'Elseparator will be converted into a mult: ' // ele%name)
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call multipole1_kt_to_ab (-val(hkick$), 0.0_rp, 0.0_rp, 0, a, b)
      a_pole = a_pole + a;  b_pole = b_pole + b
      call multipole1_kt_to_ab (-val(vkick$), pi/2, 0.0_rp, 0, a, b)
      a_pole = a_pole + a;  b_pole = b_pole + b

    !
    case (hkicker$)
      write (line_out, '(4a)') 'BEND ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, -val(kick$), 'K0', 'R', .true., .false.)
!      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
!      call multipole1_kt_to_ab (-val(kick$), 0.0_rp, 0.0_rp, 0, a, b)
!      a_pole = a_pole + a;  b_pole = b_pole + b

    !
    case (vkicker$)
      tilt = -val(tilt$) - pi/2
      write (line_out, '(4a)') 'BEND ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, -val(kick$), 'K0', 'R', .true., .false.)
      call value_to_line (line_out, tilt, 'ROTATE', 'R', .true., .false.)
!      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
!      call multipole1_kt_to_ab (-val(kick$), pi/2, 0.0_rp, 0, a, b)
!      a_pole = a_pole + a;  b_pole = b_pole + b

    !
    case (kicker$)
      write (line_out, '(4a)') 'BEND ', trim(ele%name), ' = (L = ', re_str(val(l$))
      hk = -val(hkick$)
      vk = -val(vkick$)
      tilt = atan2(hk, vk) + pi/2 - val(tilt$)
      if (hk /= 0 .or. vk /= 0) then
        call value_to_line (line_out, -sqrt(hk**2 + vk**2), 'K0', 'R', .true., .false.)
        call value_to_line (line_out, tilt, 'ROTATE', 'R', .true., .false.)
      endif
!      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
!      call multipole1_kt_to_ab (-val(hkick$), 0.0_rp, 0.0_rp, 0, a, b)
!      a_pole = a_pole + a;  b_pole = b_pole + b
!      call multipole1_kt_to_ab (-val(vkick$), pi/2, 0.0_rp, 0, a, b)
!      a_pole = a_pole + a;  b_pole = b_pole + b

    !
    case (lcavity$)
      write (line_out, '(4a)') 'CAVI ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, val(rf_frequency$), 'FREQ', 'R', .true., .false.)
      call value_to_line (line_out, val(voltage$), 'VOLT', 'R', .true., .false.)
      call value_to_line (line_out, 0.25 - val(phi0$), 'PHI', 'R', .true., .false.)
      call value_to_line (line_out, -val(phi0_err$), 'DPHI', 'R', .true., .false.)

    !
    case (marker$)
      write (line_out, '(4a)') 'MARK ', trim(ele%name), ' = ('
      call value_to_line (line_out, val(sad_mark_offset$), 'OFFSET', 'R', .true., .false.)
      call value_to_line (line_out, ele%old_value(sad_bound$), 'BOUND', 'R', .true., .false.)
      call value_to_line (line_out, ele%old_value(sad_geo$), 'GEO', 'R', .true., .false.)
      if (branch_out%param%geometry == open$ .and. ix_ele == 1) then
        call value_to_line (line_out, ele%a%beta, 'BX', 'R', .true., .false.)
        call value_to_line (line_out, ele%b%beta, 'BY', 'R', .true., .false.)
        call value_to_line (line_out, ele%a%alpha, 'AX', 'R', .true., .false.)
        call value_to_line (line_out, ele%b%alpha, 'AY', 'R', .true., .false.)
        call value_to_line (line_out, ele%x%eta, 'PEX', 'R', .true., .false.)
        call value_to_line (line_out, ele%y%eta, 'PEY', 'R', .true., .false.)
        call value_to_line (line_out, ele%x%etap, 'PEPX', 'R', .true., .false.)
        call value_to_line (line_out, ele%y%etap, 'PEPY', 'R', .true., .false.)
        call value_to_line (line_out, lat%a%emit, 'EMITX', 'R', .true., .false.)
        call value_to_line (line_out, lat%b%emit, 'EMITY', 'R', .true., .false.)
      endif

    !
    case (null_ele$)
      write (line_out, '(4a)') 'SOL ', trim(ele%name), ' = ('
      call value_to_line (line_out, ele%old_value(sad_bz$),    'BZ',    'R', .true., .false.)
      call value_to_line (line_out, ele%old_value(sad_f1$),    'F1',    'R', .true., .false.)
      call value_to_line (line_out, ele%old_value(sad_bound$), 'BOUND', 'R', .true., .false.)
      call value_to_line (line_out, ele%old_value(sad_geo$),   'GEO',   'R', .true., .false.)

    ! With nonzero multipoles or kick
    case (octupole$)
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call multipole1_kt_to_ab (val(k3$), 0.0_rp, 0.0_rp, 3, a, b)
      a_pole = a_pole + a;  b_pole = b_pole + b

    !
    case (patch$)
      write (line_out, '(4a)') 'COORD ', trim(ele%name), ' = ('

    ! With nonzero multipoles or kick
    case (quadrupole$) 
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, -sign_of(val(fq1$)) * sqrt(24*abs(val(fq1$))), 'F1', 'R', .true., .false.)
      call value_to_line (line_out, val(fq2$), 'F2', 'R', .true., .false.)

      call multipole1_kt_to_ab (val(k1$), 0.0_rp, 0.0_rp, 1, a, b)
      a_pole = a_pole + a;  b_pole = b_pole + b

    !
    case (rfcavity$)
      write (line_out, '(4a)') 'CAVI ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, val(rf_frequency$), 'FREQ', 'R', .true., .false.)
      call value_to_line (line_out, val(voltage$), 'VOLT', 'R', .true., .false.)
      call value_to_line (line_out, twopi * val(phi0$), 'DPHI', 'R', .true., .false.)

      select case (nint(val(fringe_at$)))
      case (entrance_end$)
        line_out = trim(line_out) // ' fringe = 1'
      case (exit_end$)
        line_out = trim(line_out) // ' fringe = 2'
      end select

    !
    case (sad_mult$)
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, -sign_of(val(fq1$)) * sqrt(24*abs(val(fq1$))), 'F1', 'R', .true., .false.)
      call value_to_line (line_out, val(fq2$), 'F2', 'R', .true., .false.)
      if (val(eps_step_scale$) /= 1) call value_to_line (line_out, val(eps_step_scale$), 'EPS', 'R', .true., .false.)
      call value_to_line (line_out, val(x_offset_mult$), 'DX', 'R', .true., .false.)
      call value_to_line (line_out, val(y_offset_mult$), 'DY', 'R', .true., .false.)
      if (val(x1_limit$) == val(y1_limit$)) then
        call value_to_line (line_out, val(x1_limit$), 'RADIUS', 'R', .true., .false.)
      else
        call out_io (s_warn$, r_name, 'Asymmetric x_limit vs y_limit cannot be converted for: ' // ele%name, &
                                  'Will use largest limit here.')
        if (val(x1_limit$) /= 0 .and. val(y1_limit$) /= 0) then
          call value_to_line (line_out, max(val(x1_limit$), val(y1_limit$)), 'RADIUS', 'R', .true., .false.)
        endif
      endif

    !
    case (sbend$)
      write (line_out, '(4a)') 'BEND ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call value_to_line (line_out, val(angle$), 'ANGLE', 'R', .true., .false.)
      call value_to_line (line_out, val(dg$)*val(l$), 'K0', 'R', .true., .false.)
      call value_to_line (line_out, val(k1$)*val(l$), 'K1', 'R', .true., .false.)
      call value_to_line (line_out, -val(ref_tilt$), 'ROTATE', 'R', .true., .false.)
      if (val(fintx$)*val(hgapx$) == val(fint$)*val(hgap$)) then
        call value_to_line (line_out, 12*val(fint$)*val(hgap$), 'F1', 'R', .true., .false.)
      else
        call value_to_line (line_out, 12*val(fint$)*val(hgap$), 'FB1', 'R', .true., .false.)
        call value_to_line (line_out, 12*val(fintx$)*val(hgapx$), 'FB2', 'R', .true., .false.)
      endif
      if (val(angle$) == 0) then
        call value_to_line (line_out, val(e1$), 'AE1', 'R', .true., .false.)
        call value_to_line (line_out, val(e2$), 'AE2', 'R', .true., .false.)
      else
        call value_to_line (line_out, val(e1$)/val(angle$), 'E1', 'R', .true., .false.)
        call value_to_line (line_out, val(e2$)/val(angle$), 'E2', 'R', .true., .false.)
      endif

      select case (nint(val(fringe_type$)))
      case (hard_edge_only$)
        ! Nothing to be done
      case (soft_edge_only$)
        line_out = trim(line_out) // ' FRINGE = 1 DISFRIN = 1'
      case (sad_full$)
        line_out = trim(line_out) // ' FRINGE = 1'
      end select

    case (beambeam$)
      write (line_out, '(4a)') 'BEAMBEAM ', trim(ele%name), ' = ('

    ! with nonzero multipoles or kick
    case (sextupole$)
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call multipole1_kt_to_ab (val(k3$), 0.0_rp, 0.0_rp, 1, a, b)
      a_pole = a_pole + a;  b_pole = b_pole + b

    !
    case (solenoid$)
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))

    !
    case (sol_quad$)
      write (line_out, '(4a)') 'MULT ', trim(ele%name), ' = (L = ', re_str(val(l$))
      call multipole1_kt_to_ab (val(k1$), 0.0_rp, 0.0_rp, 1, a, b)
      a_pole = a_pole + a;  b_pole = b_pole + b

    !
    case default
      call out_io (s_error$, r_name, 'I DO NOT KNOW HOW TO CONVERT AN ELEMENT OF TYPE: ' // key_name(ele%key), &
             'CONVERTING TO DRIFT')
      write (line_out, '(4a)') 'DRIFT ', trim(ele%name), ' = (L = ', re_str(val(l$))
    end select

    if (line_out(1:4) == 'MULT') then
      if (has_attribute (ele, 'HKICK') .and. ele%key /= kicker$) then
        call multipole1_kt_to_ab (-val(hkick$), -val(tilt_tot$), 0.0_rp, 0, a, b)
        a_pole = a_pole + a;  b_pole = b_pole + b
        call multipole1_kt_to_ab (-val(vkick$), pi/2-val(tilt_tot$), 0.0_rp, 0, a, b)
        a_pole = a_pole + a;  b_pole = b_pole + b
      endif

      do i = 0, 21
        write (str, '(i0)') i
        call value_to_line (line_out, b_pole(i)*factorial(i), 'K'//trim(str), 'R', .true., .false.)
        call value_to_line (line_out, a_pole(i)*factorial(i), 'SK'//trim(str), 'R', .true., .false.)
      enddo
    endif
  endif     ! Not converted

  ! fringe

  if (ele%key == quadrupole$ .or. ele%key == sad_mult$) then
    select case (nint(val(fringe_type$)))
    case (soft_edge_only$, full$)
      select case (nint(val(fringe_at$)))
      case (entrance_end$);         line_out = trim(line_out) // ' FRINGE = 1'
      case (exit_end$);             line_out = trim(line_out) // ' FRINGE = 2'
      case (both_ends$);            line_out = trim(line_out) // ' FRINGE = 3'
      end select
    end select

    select case (nint(val(fringe_type$)))
    case (none$, soft_edge_only$)
      line_out = trim(line_out) // ' DISFRIN = 1'
    case (hard_edge_only$)
      select case (nint(val(fringe_at$)))
      case (no_end$)
        line_out = trim(line_out) // ' DISFRIN = 1'
      case (entrance_end$, exit_end$)
        line_out = trim(line_out) // ' CANNOT TRANSLATE BMAD FRINGE_TYPE/FRINGE_AT!'
      end select
    case (full$)
      if (nint(val(fringe_at$)) == no_end$) line_out = trim(line_out) // ' DISFRIN = 1'
    end select
  endif

  ! misalignments
  ! Note: SAD applies pitches and offsets in reverse order to Bmad.

  xp = val(x_pitch$);  yp = val(y_pitch$)

  if (xp /= 0) then
    zo =  val(z_offset$) * cos(xp) + val(x_offset$) * sin(xp)
    xo = -val(z_offset$) * sin(xp) + val(x_offset$) * cos(xp)
    val(z_offset$) = zo
    val(x_offset$) = xo
  endif

  if (yp /= 0) then
    zo =  val(z_offset$) * cos(yp) + val(y_offset$) * sin(yp)
    yo = -val(z_offset$) * sin(yp) + val(y_offset$) * cos(yp)
    val(z_offset$) = zo
    val(y_offset$) = yo
  endif

  if (ele%key == null_ele$) then ! patch -> SOL
    if (ele%iyy == entrance_end$) then
      val(x_offset$) = -val(x_offset$)
      val(y_offset$) = -val(y_offset$)
      val(z_offset$) = -val(z_offset$)
    else
      val(z_offset$) = -val(z_offset$)
    endif
    val(z_offset$) = val(z_offset$) + ele%value(t_offset$) * c_light
  endif

  call value_to_line (line_out, val(x_offset$), 'DX', 'R', .true., .false.)
  call value_to_line (line_out, val(y_offset$), 'DY', 'R', .true., .false.)
  call value_to_line (line_out, val(z_offset$), 'DZ', 'R', .true., .false.)

  if (ele%key /= marker$) then
    call value_to_line (line_out, -val(x_pitch$), 'CHI1', 'R', .true., .false.)
    call value_to_line (line_out, -val(y_pitch$), 'CHI2', 'R', .true., .false.)
  endif

  if (ele%key == patch$ .or. ele%key == null_ele$) then   ! null_ele -> SOL
    call value_to_line (line_out, -val(tilt$),    'CHI3', 'R', .true., .false.)
  elseif (ele%key /= sbend$ .and. ele%key /= kicker$ .and. ele%key /= vkicker$) then
    call value_to_line (line_out, -val(tilt$),    'ROTATE', 'R', .true., .false.)
  endif

  if (ele%key == null_ele$ .and. nint(ele%old_value(sad_geo$)) == 1) then
    if (ele%iyy == entrance_end$) then
      call value_to_line (line_out, val(x_pitch$), 'DPX', 'R', .true., .false.)
      call value_to_line (line_out, val(y_pitch$), 'DPY', 'R', .true., .false.)
    else
      call value_to_line (line_out, -val(x_pitch$), 'DPX', 'R', .true., .false.)
      call value_to_line (line_out, -val(y_pitch$), 'DPY', 'R', .true., .false.)
    endif
  endif      

  !

  line_out = trim(line_out) // ')'
  call write_line(line_out)
  cycle

  ! write element spec to file

  call write_line(line_out)

enddo

!---------------------------------------------------------------------------------------
! Write the lattice line

write (iu, '(a)')
write (iu, '(3a)') '!---------------------------------;'
write (iu, '(a)')
write (line_out, '(2a)') 'LINE ASC = ('

do n = ie1, ie2

  ele => branch_out%ele(n)
  if (ele%sub_key == not_set$) cycle

  ix = len_trim(line_out) + len_trim(ele%name)

  if (ix > 75) then
    write (iu, '(3a)') trim(line_out)
    line_out = '   ' // ele%name
  else
    line_out = trim(line_out) // ' ' // ele%name
  endif

  ! Output line if long enough or at end

  if (n == ie2) then
    line_out = trim(line_out) // ')'
    write (iu, '(2a)') trim(line_out), ';'
    line_out = ' '
    init_needed = .true.
  endif

enddo

!------------------------------------------
! Use statement

write (iu, '(a)')
write (iu, '(3a)') '!---------------------------------;'
write (iu, '(a)')

write (iu, '(a)') 'FFS USE ASC;'

! End stuff

call out_io (s_info$, r_name, 'Written SAD lattice file: ' // trim(out_file_name))

deallocate (names)
if (present(err)) err = .false.

if (present(converted_lat)) then
  converted_lat = lat
  converted_lat%branch(branch%ix_branch) = branch_out
  converted_lat%n_ele_max = converted_lat%n_ele_track
  do ib = 0, ubound(converted_lat%branch, 1)
    branch => converted_lat%branch(ib)
    do i = 1, branch%n_ele_track
      branch%ele(i)%slave_status = free$
      branch%ele(i)%n_lord = 0
    enddo
  enddo
  converted_lat%n_control_max = 0
  converted_lat%n_ic_max = 0
endif

call deallocate_lat_pointers (lat_out)
call deallocate_lat_pointers (lat_model)

! Restore ptc settings

if (n_taylor_order_saved /= ptc_private%taylor_order_ptc) call set_ptc (taylor_order = n_taylor_order_saved) 
ptc_com%exact_model = ptc_exact_model

close (iu)

!------------------------------------------------------------------------
contains

subroutine write_line (line_out)

implicit none

character(*) line_out
integer ix, ix1, ix2, ix3

! Prefer to breakup a line after a comma

do
  if (len_trim(line_out) < ix_line_max) exit
  ix1 = index(line_out(ix_line_min+1:), ',')
  ix2 = index(line_out(ix_line_min+1:), '=')
  ix3 = index(line_out(ix_line_min+1:), ' ')

  if (ix1 /= 0 .and. ix1+ix_line_min < ix_line_max) then
    ix = ix1 + ix_line_min
  elseif (ix2 /= 0 .and. ix2+ix_line_min < ix_line_max) then
    ix = ix2 + ix_line_min
  elseif (ix3 /= 0 .and. ix3+ix_line_min < ix_line_max) then
    ix = ix3 + ix_line_min
  elseif (ix1 /= 0) then
    ix = ix1 + ix_line_min
  elseif (ix2 /= 0) then
    ix = ix2 + ix_line_min
  else
    ix = ix3 + ix_line_min
  endif

  write (iu, '(2a)') line_out(:ix)
  line_out = '    ' // line_out(ix+1:)
enddo

write (iu, '(2a)') trim(line_out), ';'

end subroutine write_line

end subroutine write_lat_in_sad_format

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

