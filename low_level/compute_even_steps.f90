!+
! Subroutine compute_even_steps (ds_in, length, ds_default, ds_out, n_step)
!
! Subroutine to compute a step size ds_out, close to ds_in, so that an 
! integer number of steps spans the length:
!   length = ds_out * n_step
!
! Modules needed:
!   use bmad
!
! Input:
!   ds_in      -- Real(rp): Input step size.
!   length     -- Real(rp): Total length.
!   ds_default -- Real(rp): Default to use if ds_in = 0.
!
! Output:
!   ds_out    -- Real(rp): Step size to use.
!   n_step    -- Integer: Number of steps needed.
!-

subroutine compute_even_steps (ds_in, length, ds_default, ds_out, n_step)

use sim_utils

implicit none

real(rp) ds_in, length, ds_default, ds_out
integer n_step

!

ds_out = ds_in
if (ds_out == 0) ds_out = ds_default
n_step = abs(nint(length / ds_out))
if (n_step == 0) n_step = 1
ds_out = length / n_step  

end subroutine compute_even_steps
