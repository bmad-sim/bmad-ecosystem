!+
! Subroutine tao_hook_evaluate_a_datum (found, datum, u, tao_lat, datum_value, valid_value)
!
! See the Programmer's manual for how to add custom data types here.
!
! Input:
!   datum        -- tao_data_struct: The current datum to evaluate
!   u            -- tao_universe_struct: Universe this datum is in
!   tao_lat      -- tao_lattice_struct: Lattice to use.
!
! Output:
!   found        -- Logical: True if  this datum is evaluated in this subroutine.
!   datum_value  -- real(rp): Which datum value to compute (model_value, design_value, etc...)
!   valid_value  -- Logical: Set False when there is a problem. Set True otherwise.
!   why_invalid  -- Character(*), optional: Tells why datum value is invalid.
!-

subroutine tao_hook_evaluate_a_datum (found, datum, u, tao_lat, datum_value, valid_value, why_invalid)

use tao_data_and_eval_mod, dummy => tao_hook_evaluate_a_datum

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct) datum
type (tao_lattice_struct), target :: tao_lat

real(rp) datum_value
logical found, valid_value

character(*), optional :: why_invalid
character(*), parameter :: r_name = 'tao_hook_evaluate_a_datum'

!

found = .false.

end subroutine tao_hook_evaluate_a_datum
