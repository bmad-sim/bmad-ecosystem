!+
! Subroutine tao_set_vars (s, var_vec)
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
!
! Output:
!   s  -- tao_super_universe_struct:
!-

subroutine tao_set_vars (s, var_vec)

use tao_mod

implicit none

type (tao_super_universe_struct), target :: s
real(rp) var_vec(:)

integer i, j, k

! Transfer the values from var_vec to the variables of each universe.

  j = 0
  do i = 1, size(s%var)
    if (.not. s%var(i)%useit_opt) cycle
    j = j + 1
    call tao_set_var_model_value (s, s%var(i), var_vec(j))
  enddo

end subroutine
