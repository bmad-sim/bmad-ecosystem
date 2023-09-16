!+
! Function set_parameter (param_val, set_val) result (save_val)
!
! Routine to simultaneously set a parameter and return the old parameter value.
!
! This is the overloaded name for:
!   Function set_parameter_real  (param_val, set_val) result (save_val)
!   Function set_parameter_int   (param_val, set_val) result (save_val)
!   Function set_parameter_logic (param_val, set_val) result (save_val)
!
! Input:
!   set_val     -- real(rp)/integer/logical: Value to set to.
!   param_val   -- real(rp)/integer/logical: Parameter with some current value to save.
!
! Output:
!   param_val   -- real(rp)/integer/logical: Parameter with value set to set_val.
!   save_val    -- real(rp)/integer/logical: Old value of param_val.
!-

function set_parameter_real (param_val, set_val) result (save_val)

use precision_def
implicit none
real(rp) param_val, set_val, save_val

save_val = param_val
param_val = set_val

end function set_parameter_real

!--------------------------------------------------------------------------------

function set_parameter_int (param_val, set_val) result (save_val)

use precision_def
implicit none
integer param_val, set_val, save_val

save_val = param_val
param_val = set_val

end function set_parameter_int

!--------------------------------------------------------------------------------

function set_parameter_logic (param_val, set_val) result (save_val)

use precision_def
implicit none
logical param_val, set_val, save_val

save_val = param_val
param_val = set_val

end function set_parameter_logic
