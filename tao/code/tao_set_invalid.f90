!+
! Subroutine tao_set_invalid (datum, message, why_invalid, exterminate, err_level, print_err)
!
! Routine to either print an error message to the terminal (if why_invalid
! is not present) or set the why_invalid string to the error message.
!
! Note: The exterminate argument should be set to False if the datum is invalid for
! reasons like the beam has been lost. In this case, the datum could possibly
! become valid in the future. Exterminate should be set to True if the datum could
! not possibly become valid. For example, the datum's reference element does not
! exist in the lattice.
!
! Input:
!   datum       -- tao_data_struct: Bad datum.
!   message     -- character(*): Error message.
!   exterminate -- logical, optional: Default is False. If True, set datum%exists
!                   to False so that Tao will ignore this datum from now on.
!   err_level   -- integer, optional: s_error$ (default), s_warn$, etc.
!   print_err   -- logical, optional: Default is True. If False, do not print an error message.
!
! Output:
!   why_invalid -- character(*), optional: Set to message if present.
!-

subroutine tao_set_invalid (datum, message, why_invalid, exterminate, err_level, print_err)

use tao_interface, dummy => tao_set_invalid

type (tao_data_struct) datum
type (tao_universe_struct), pointer :: u

logical, optional :: exterminate, print_err
logical identified_err_found, found

integer, optional :: err_level
integer i

character(*) message
character(*), optional :: why_invalid
character(*), parameter :: r_name = 'tao_set_invalid'
character(40) :: err_str(2) = [character(40):: 'NO BEAM TRACKING HAS BEEN DONE', &
                                               'CANNOT EVALUATE DUE TO PARTICLE LOSS']

! The idea with err_str is to limit the number of error messages generated of a given type.

datum%why_invalid = message

if (logic_option(.false., exterminate)) then
  datum%exists = .false. 
endif

if (present(why_invalid)) then
  why_invalid = message

elseif (.not. datum%err_message_printed) then
  if (.not. logic_option(.true., print_err)) return

  datum%err_message_printed = .true.
  identified_err_found = .false.
  do i = 1, size(err_str)
    if (index(message, trim(err_str(i))) == 0) cycle
    identified_err_found = .true.
    if (s%com%is_err_message_printed(i)) return
    s%com%is_err_message_printed(i) = .true.
    identified_err_found = .true.
  enddo

  s%com%n_err_messages_printed = s%com%n_err_messages_printed + 1
  if (s%com%n_err_messages_printed == s%global%datum_err_messages_max + 1) then
    call out_io (s_warn$, r_name, 'NUMBER OF ERROR MESSAGES EXCEEDS GLOBAL%DATUM_ERR_MESSAGES_MAX.', &
                                  'WILL NOT PRINT ANY MORE DATUM ERROR MESSAGES FOR THIS EVALUATION CYCLE.')
  endif
  if (s%com%n_err_messages_printed > s%global%datum_err_messages_max) return

  call out_io (integer_option(s_error$, err_level), r_name, message, &
                'FOR DATUM: ' // trim(tao_datum_name(datum)) // ' with data_type: ' // datum%data_type)

  if (index(upcase(message), 'UNSTABLE') /= 0) then
    u => tao_pointer_to_universe(datum%ix_uni)
    found = .false.
    do i = 1, size(u%data)
      if (substr(u%data(i)%data_type,1,8) == 'unstable') found = .true.
    enddo
    if (.not. found) call out_io(s_info$, r_name, &
              'NO unstable.orbit, unstable.ring, nor unstable.eigen FOR THE UNIVERSE WITH THE PROBLEM EXISTS.', &
              'YOU MIGHT WANT TO CONSIDER ADDING SUCH A DATUM.')
  endif

  if (identified_err_found) then
    call out_io (s_warn$, r_name, 'WILL NOT PRINT ANY MORE OF THIS KIND OF DATUM ERROR MESSAGE FOR THIS EVALUATION CYCLE.')
  endif
endif

end subroutine tao_set_invalid
