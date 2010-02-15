!+
! Subroutine tao_hook_evaluate_a_datum (found, datum, u, tao_lat, datum_value, valid_value)
!
! See the Programmer's manual for how to add custom data types here.
!
! Input:
!   datum        -- tao_data_struct: the current datum to evaluate
!   u            -- tao_universe_struct: universe this datum is in
!   tao_lat      -- Tao_lattice_struct: Lattice to use.
!
! Output:
!   datum_value  -- real(rp): which datum value to compute (model_value,
!                             design_value, etc...)
!   Found        -- Logical: TRUE if  this datum is evaluated in this subroutine.
!   valid_value  -- Logical: Set false when there is a problem. Set true otherwise.
!   why_invalid  -- Character(*), optional: Tells why datum value is invalid.
!-

subroutine tao_hook_evaluate_a_datum (found, datum, u, tao_lat, datum_value, valid_value, why_invalid)

use tao_mod

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct) datum
type (tao_lattice_struct), target :: tao_lat

real(rp) datum_value
logical found, valid_value

character(*), optional :: why_invalid
character(20) :: r_name = 'tao_hook_evaluate_a_datum'

!

found = .false.

end subroutine tao_hook_evaluate_a_datum
