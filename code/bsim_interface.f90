module bsim_interface
  
use bmad_interface

interface

function set_tune_3d (branch, target_tunes, mask, use_phase_trombone, z_tune_set, print_err) result (everything_ok)
  import
  implicit none
  type (branch_struct), target :: branch
  real(rp) target_tunes(3)
  logical everything_ok
  logical, optional :: use_phase_trombone, z_tune_set, print_err
  character(*), optional :: mask
end function

subroutine insert_phase_trombone(branch)
  import
  implicit none
  type (branch_struct), target :: branch
end subroutine
    
end interface
  
integer, private :: private_dummy ! This is to suppress the ranlib "has no symbols" message

end module bsim_interface


