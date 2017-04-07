!+
! Function tao_data_sanity_check (datum, print_err) result (is_valid)
!
! Routine to check if the data is internally consistent.
! Note: A datum whose data_type demands an associated lattice element will be invalid but will
! not generate an error message since d1 data arrays my have invalid datums that are just place holders.
!
! Input:
!   datum     -- tao_data_struct: Datum to check.
!   print_err -- logical: Print error message if data is not valid?
!
! Output:
!   is_valid  -- logical: True if internally consistent.
!-

function tao_data_sanity_check (datum, print_err) result (is_valid)

use tao_utils, dummy => tao_data_sanity_check

implicit none

type (tao_data_struct) datum
type (branch_struct), pointer :: branch
type (tao_universe_struct), pointer :: u

logical print_err, is_valid
character(*), parameter :: r_name = 'tao_data_sanity_check'

!

u => s%u(datum%d1%d2%ix_uni)
branch => u%design%lat%branch(datum%ix_branch)

is_valid = .false.

if (datum%ele_name /= '') then
  if (branch%ele(datum%ix_ele)%name /= datum%ele_name) then
    if (print_err) call out_io (s_error$, r_name, 'DATUM ELEMENT NOT LOCATED: ' // datum%ele_name, &
                                                  'FOR DATUM: ' // tao_datum_name(datum))
    return
  endif
endif

if (datum%ele_ref_name /= '') then
  if (branch%ele(datum%ix_ele_ref)%name /= datum%ele_ref_name) then
    if (print_err) call out_io (s_error$, r_name, 'DATUM ELEMENT REFERENCE NOT LOCATED: ' // datum%ele_ref_name, & 
                                                  'FOR DATUM: ' // tao_datum_name(datum))
    return
  endif
endif

if (datum%ele_start_name /= '') then
  if (branch%ele(datum%ix_ele_start)%name /= datum%ele_start_name) then
    if (print_err) call out_io (s_error$, r_name, 'DATUM ELEMENT START NOT LOCATED: ' // datum%ele_start_name, &
                                                  'FOR DATUM: ' // tao_datum_name(datum))
    return
  endif
endif

!

if (tao_chrom_calc_needed(datum%data_type, datum%data_source)) then
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

if (datum%data_type == '') then
  if (print_err) call out_io (s_error$, r_name, 'DATA_TYPE NOT SET FOR DATUM: ' // tao_datum_name(datum))
  return
endif

if (datum%data_source == '') then
  if (print_err) call out_io (s_error$, r_name, 'DATA_SOURCE NOT SET FOR DATUM: ' // tao_datum_name(datum))
  return
endif


!

if (datum%data_type == 'unstable.orbit') then
  if (datum%ele_name == '') then
    call out_io (s_abort$, r_name, 'DATUM OF TYPE: ' // datum%data_type, &
                                   'CANNOT HAVE AN ASSOCIATED ELEMENT: ' // datum%ele_name, &
                                   'FOR DATUM: ' // tao_datum_name(datum))
    return
  endif

elseif (branch%param%geometry == closed$ .and. &
            (datum%data_type(1:12) == 'chrom.dtune.' .or. datum%data_type(1:5) == 'damp.' .or. &
             datum%data_type(1:17) == 'multi_turn_orbit.' .or. datum%data_type(1:5) == 'tune.' .or. &
             datum%data_type(1:13) == 'unstable.ring' .or. index(datum%data_type, 'emit.') /= 0)) then
  if (datum%ele_name /= '') then
    if (print_err) call out_io (s_abort$, r_name, 'DATUM OF TYPE: ' // datum%data_type, &
                                                  'CANNOT HAVE AN ASSOCIATED ELEMENT IN A CIRCULAR LATTICE: ' // datum%ele_name, &
                                                  'FOR DATUM: ' // tao_datum_name(datum))
    return
  endif

else
  if (datum%ele_name == '') then
    ! Datum is invalid but there is no error.
    !if (print_err) call out_io (s_abort$, r_name, 'DATUM OF TYPE: ' // datum%data_type, &
    !                                              'MUST HAVE AN ASSOCIATED ELEMENT.', &
    !                                              'FOR DATUM: ' // tao_datum_name(datum))
    return
  endif
endif

is_valid = .true.

end function

