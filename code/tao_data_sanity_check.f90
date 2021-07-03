!+
! Function tao_data_sanity_check (datum, print_err, default_data_type) result (is_valid)
!
! Routine to check if the data is internally consistent.
! Note: A datum whose data_type demands an associated lattice element will be invalid but will
! not generate an error message since d1 data arrays my have invalid datums that are just place holders.
!
! Input:
!   datum               -- tao_data_struct: Datum to check.
!   print_err           -- logical: Print error message if data is not valid?
!   default_data_type   -- character(*): Default data type associated with the datum's d2 structure.
!
! Output:
!   is_valid  -- logical: True if internally consistent.
!-

function tao_data_sanity_check (datum, print_err, default_data_type) result (is_valid)

use tao_interface, dummy => tao_data_sanity_check

implicit none

type (tao_data_struct) datum
type (branch_struct), pointer :: branch
type (tao_universe_struct), pointer :: u

integer has_associated_ele
logical print_err, is_valid
character(40) d_type
character(*) default_data_type
character(*), parameter :: r_name = 'tao_data_sanity_check'

!

u => s%u(datum%d1%d2%ix_universe)
branch => u%design%lat%branch(datum%ix_branch)

is_valid = .false.
d_type = datum%data_type

if (datum%ele_name /= '') then
  if (.not. check_ele_ok (datum%ele_name, datum%ix_ele, 'DATUM ELEMENT NOT FOUND: ' // datum%ele_name)) return
endif

if (datum%ele_ref_name /= '') then
  if (.not. check_ele_ok (datum%ele_ref_name, datum%ix_ele_ref, 'DATUM ELEMENT REFERENCE NOT FOUND: ' // datum%ele_ref_name)) return
endif

if (datum%ele_start_name /= '') then
  if (.not. check_ele_ok (datum%ele_start_name, datum%ix_ele_start, 'DATUM ELEMENT START NOT FOUND: ' // datum%ele_start_name)) return
endif

!

if (tao_chrom_calc_needed(d_type, datum%data_source)) then
  if (branch%param%geometry == open$) then
    if (print_err) call out_io (s_error$, r_name, 'CHROMATICITY DATUM NOT VALID FOR NON-CLOSED LATTICE!')
    return
  endif  
endif

!

if (datum%merit_type == '') then
  if (print_err) call out_io (s_error$, r_name, 'MERIT_TYPE NOT SET FOR DATUM: ' // tao_datum_name(datum))
  return
endif

if (d_type == '') then
  if (print_err) call out_io (s_error$, r_name, 'DATA_TYPE NOT SET FOR DATUM: ' // tao_datum_name(datum))
  return
endif

if (datum%data_source == '') then
  if (print_err) call out_io (s_error$, r_name, 'DATA_SOURCE NOT SET FOR DATUM: ' // tao_datum_name(datum))
  return
endif


!

has_associated_ele = tao_datum_has_associated_ele(d_type, branch%param%geometry)

if (has_associated_ele == maybe$) then
  ! Do nothing

elseif (has_associated_ele == no$) then
  if (datum%ele_name /= '') then
    call out_io (s_abort$, r_name, 'DATUM OF TYPE: ' // d_type, &
                                   'CANNOT HAVE AN ASSOCIATED ELEMENT: ' // datum%ele_name, &
                                   'FOR DATUM: ' // tao_datum_name(datum))
    return
  endif

else    ! has_associated_ele = yes$. Also invalid data_types here.
  if (datum%ele_name == '') then
    if (datum%data_type == default_data_type .or. datum%data_type == '') then
      ! Datum is invalid but this is do not generate an error since having "gaps" in the d1 array is a common situation.
      return
    endif
  endif
endif

is_valid = .true.

!-----------------------------------------------------------------------------------------------------
contains

function check_ele_ok(ele_name, ix_ele, err_str) result (is_ok)

character(*) ele_name, err_str
integer ix_ele
logical is_ok

!

is_ok = .false.

if (ix_ele < 0) then
  if (print_err) call out_io (s_error$, r_name, err_str, 'FOR DATUM: ' // tao_datum_name(datum))
  return
endif

! Problem with this test is that ele_name may be something like "1>>345".

!if (branch%ele(datum%ix_ele)%name /= datum%ele_name) then
!  if (print_err) call out_io (s_error$, r_name, err_str // ' LOCATION CONFUSION: ' // datum%ele_name, &
!                                                'FOR DATUM: ' // tao_datum_name(datum))
!  return
!endif

is_ok = .true.

end function check_ele_ok

end function tao_data_sanity_check

