!+
! Subroutine tao_set_opt_vars (var_vec, print_limit_warning)
!
! Subrutine to set variable values from a vector of values. 
! This routine is used with optimization since optimimizers
! generally like their variables in the form of a vector.
!
! This routine assumes that the variables in each s%u(i) universe
! gets the same values. 
!
! Input:
!   var_vec(:) -- Real(rp): Vector of variables. 
!   print_limit_warning
!         -- Logical, optional: Print a warning if the value is past the variable's limits.
!             Default is True.
!
! Output:
!-

subroutine tao_set_opt_vars (var_vec, print_limit_warning)

use tao_interface, dummy => tao_set_opt_vars

implicit none

real(rp) var_vec(:)
integer i, j, k
logical, optional :: print_limit_warning

! Transfer the values from var_vec to the variables of each universe.

j = 0
do i = 1, s%n_var_used
  if (.not. s%var(i)%useit_opt) cycle
  j = j + 1
  call tao_set_var_model_value (s%var(i), var_vec(j), print_limit_warning)
enddo

end subroutine tao_set_opt_vars
