!+
! Subroutine transfer_mat_from_tracking (ele, param, orb0, d_orb, error)
!
! Subroutine to compute the transfer map for an element from tracking.
! Because of nonlinearities, and the fact that d_orb is finite,
! the computed 6x6 matrix will not be exactly symplectic. To handle this
! the matrix is symplecitified.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele      -- Ele_struct: Element to track through.
!   param    -- param_struct: Needed for track1.
!   orb0     -- Coord_struct: Orbit at beginning of element that is used as
!                the origin for the analysis.
!   d_orb(6) -- Real(rdef): Vector of offsets to use. All 6 components must be
!                non-zero.
!
! Output:
!   ele     -- Ele_struct: 
!     %mat6   -- 6x6 transfer matrix.
!   error   -- Real(rdef): A measure of how symplectic the constructed matrix was
!              before symplecitification. See mat_symp_check for more details.
!-

!$Id$
!$Log$
!Revision 1.3  2002/02/23 20:32:27  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:59  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine transfer_mat_from_tracking (ele, param, orb0, d_orb, error)  

  use bmad

  implicit none

  type (ele_struct), intent(inout) :: ele
  type (param_struct), intent(in) :: param
  type (coord_struct), intent(in) :: orb0, d_orb
  real(rdef), intent(out) :: error

  type (coord_struct) :: start, end, end0

  integer i

! find 6x6 transfer matrix through tracking

  call track1 (orb0, ele, param, end0)

  do i = 1, 6
    start = orb0 
    start%vec(i) = start%vec(i) + d_orb%vec(i)
    call track1 (start, ele, param, end)
   ele%mat6(1:6, i) = (end%vec - end0%vec) / d_orb%vec(i)
  enddo

! symplectify

  call mat_symp_check (ele%mat6, error)
  call mat_symplectify (ele%mat6, ele%mat6)

end subroutine
