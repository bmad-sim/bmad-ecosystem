!+         
! Subroutine one_turn_matrix (ring, t1)
!
! Subroutine to calculate the full 6X6 1 turn matrix
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring     -- Ring_struct: Ring
!
! Output:
!    t1(n,n) -- Real(rp): 1-turn matrix.
!                   n = 4  -> Transverse matrix esentually with RF off.
!                   n = 6  -> Full 6x6 matrix.
!-

#include "CESR_platform.inc"

subroutine one_turn_matrix (ring, t1)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct)  ring

  real(rp) t1(:,:)
  integer i

!

  if (size(t1,1) == 4 .and. size(t1,2) == 4) then
    t1 = ring%ele_(1)%mat6(1:4,1:4)
    do i=2,ring%n_ele_ring
      t1 = matmul (ring%ele_(i)%mat6(1:4,1:4), t1)
    end do

  elseif (size(t1,1) == 6 .and. size(t1,2) == 6) then
    t1 = ring%ele_(1)%mat6
    do i=2,ring%n_ele_ring                
      t1 = matmul (ring%ele_(i)%mat6, t1)
    end do

  else
    print *, 'ERROR IN ONE_TURN_MATRIX: MATRIX HAS THE WRONG SIZE:', &
                                                    size(t1,1), size(t1,2) 
    call err_exit
  endif

end subroutine
                                          
