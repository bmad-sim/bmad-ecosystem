!+
! Function set_parameter (param_val, set_val) result (save_val)
!
! Routine to simultaneously set a parameter and return the old parameter value.
!
! Input:
!   set_val     -- real(rp): Value to set to.
!   param_val   -- real(rp): Parameter with some current value to save.
!
! Output:
!   param_val   -- real(rp): Parameter with value set to set_val.
!   save_val    -- real(rp): Old value of param_val.
!-

function set_parameter (param_val, set_val) result (save_val)

use precision_def
implicit none
real(rp) param_val, set_val, save_val

save_val = param_val
param_val = set_val

end function set_parameter
