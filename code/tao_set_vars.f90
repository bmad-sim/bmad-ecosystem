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

! parallel_vars means that all variables s%u(:)%var(j) of different universes
! are in lock-step with each other.

if (s%global%parallel_vars) then

  j = 0
  do i = 1, size(s%u(1)%var)
    if (.not. s%u(1)%var(i)%useit_opt) cycle
    j = j + 1
    do k = 1, size(s%u)
      s%u(k)%var(i)%model_value = var_vec(j)
    enddo
  enddo

! Serial case where the variables of different universes are not correlated.

else

  j = 0
  do k = 1, size(s%u)
    do i = 1, size(s%u(k)%var)
      if (.not. s%u(k)%var(i)%useit_opt) cycle
      j = j + 1
      s%u(k)%var(i)%model_value = var_vec(j)
    enddo
  enddo

endif

end subroutine
