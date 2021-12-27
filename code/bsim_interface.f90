module bsim_interface
  
use bmad_interface

interface

subroutine set_tune3(lat, target_tunes, use_phase_trombone, mask,  everything_ok)
  import
  implicit none
  type (lat_struct), target :: lat
  real(rp) target_tunes(3)
  logical everything_ok, use_phase_trombone
  character(*) :: mask
end subroutine set_tune3

subroutine insert_phase_trombone(branch)
  import
  implicit none
  type (branch_struct), target :: branch
end subroutine
    
end interface
  
integer, private :: private_dummy ! This is to suppress the ranlib "has no symbols" message

end module bsim_interface


