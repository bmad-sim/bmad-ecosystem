!+
! Subroutine track_taylor (start, bmad_taylor, end)
!
! Subroutine to track using a taylor series
!
! Modules needed:
!   use bmad
!
! Input:
!   bmad_taylor(6) -- Taylor_struct: Taylor series.
!   start          -- Coord_struct: Starting coords.
!
! Output:
!   end       -- Coord_struct: Ending coords.
!-

#include "CESR_platform.inc"   

subroutine track_taylor (start, bmad_taylor, end)

  use bmad
  
  implicit none
  
  type (taylor_struct), intent(in) :: bmad_taylor(6)
  type (coord_struct), intent(in) :: start
  type (coord_struct), intent(out) :: end
  
  real(rdef) vec_in(6), delta
  
  integer i, j, k, ie
  
!

  vec_in = start%vec
  
  end%vec = 0
  do i = 1, 6
    j_loop: do j = 1, size(bmad_taylor(i)%term)
      delta =  bmad_taylor(i)%term(j)%coef
      do k = 1, 6
        ie = bmad_taylor(i)%term(j)%exp(k) 
        if (ie == 0) cycle
        if (start%vec(k) == 0) cycle j_loop  ! delta = 0 in this case 
        delta = delta * start%vec(k) ** ie
      enddo
      end%vec(i) = end%vec(i) + delta
    enddo j_loop
  enddo

end subroutine
