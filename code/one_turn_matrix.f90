!+         
! Subroutine one_turn_matrix (ring, rf_on, t1)
!
! Subroutine to calculate the full 6X6 1 turn matrix from the individual
! element tranfer matrices.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring     -- Ring_struct: Ring
!     %ele_(:)%mat6  -- Transfer matrices to use.
!   rf_on    -- Logical: False is what is needed to be able to calculate the
!                 Twiss parameters.
!
! Output:
!    t1(6,6) -- Real(rp): 1-turn matrix.
!-

#include "CESR_platform.inc"

subroutine one_turn_matrix (ring, rf_on, t1)

  use bmad_struct
  use bmad_interface, except => one_turn_matrix

  implicit none

  type (ring_struct)  ring

  real(rp), intent(out) :: t1(:,:)
  real(rp) rf_mat(6,6)

  logical, intent(in) :: rf_on
  integer i

!

  call mat_make_unit (t1)

  do i=1, ring%n_ele_use
    if (ring%ele_(i)%key == rfcavity$ .and. .not. rf_on) then
      rf_mat = ring%ele_(i)%mat6
      rf_mat(6,5) = 0  ! turn rf off
      t1 = matmul (rf_mat, t1)
    else
      t1 = matmul (ring%ele_(i)%mat6, t1)
    endif
  end do


end subroutine
                                          
