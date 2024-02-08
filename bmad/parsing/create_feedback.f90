!+
! Subroutine create_feedback (lord, modulator, kicker, err)
!
! Subroutine to add the lord/slave bookkeeping information for a feedback lord element.
!
! Input:
!   lord          -- ele_struct: Cooler element.
!   modulator     -- character(*): Name of modulator slave.
!   kicker        -- character(*): Name of kicker slave.
!   err           -- Logical: Set True if there is a problem.
!
! Output:
!   lord          -- ele_struct: Modified feedback elment.
!-

subroutine create_feedback (lord, modulator, kicker, err)

use bmad_parser_mod, except_dummy => create_feedback

implicit none

type (ele_struct), target :: lord
type (lat_struct), pointer :: lat

character(*) modulator, kicker

logical err

! Error check

lat => lord%branch%lat
err = .true.


end subroutine


