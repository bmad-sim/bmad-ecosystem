!+
! Subroutine transfer_mat_from_tracking (ele, param, start, d_orb, end, error)
!
! Subroutine to compute the transfer map for an element from tracking.
! Note: Because of nonlinearities, and the fact that d_orb is finite,
! the computed 6x6 matrix will not be exactly symplectic.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele      -- Ele_struct: Element to track through.
!     %symplectify -- If true the matrix will be symplecitified.
!   param    -- param_struct: Needed for track1.
!   start    -- Coord_struct: Orbit at beginning of element that is used as
!                the origin for the analysis.
!   d_orb(6) -- Real(rdef), optional: Vector of offsets to use. 
!               Default if d_orb = 0 is to set d_orb = 1e-5
!
! Output:
!   ele     -- Ele_struct: 
!     %mat6   -- 6x6 transfer matrix.
!   end     -- Coord_struct, optional: orbit at end of the element.
!   error   -- Real(rdef), optional: A measure of how symplectic the 
!               constructed matrix was before symplecitification. 
!               See mat_symp_check for more details.
!-

!$Id$
!$Log$
!Revision 1.4  2002/06/13 14:54:30  dcs
!Interfaced with FPP/PTC
!
!Revision 1.3  2002/02/23 20:32:27  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:59  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine transfer_mat_from_tracking (ele, param, start, d_orb, end, error)  

  use bmad

  implicit none

  type (ele_struct), intent(inout) :: ele
  type (param_struct), intent(inout) :: param
  type (coord_struct), intent(in) :: start
  type (coord_struct), optional, intent(out) :: end
  type (coord_struct), optional :: d_orb
  real(rdef), optional, intent(out) :: error

  type (coord_struct) :: start1, end1, end0, del_orb

  integer i

! find 6x6 transfer matrix through tracking

  call track1 (start, ele, param, end0)
  if (present(end)) end = end0

  if (present(d_orb)) del_orb = d_orb
  if (all(del_orb%vec == 0)) del_orb%vec = 1e-5

  do i = 1, 6
    start1 = start 
    start1%vec(i) = start1%vec(i) + del_orb%vec(i)
    call track1 (start1, ele, param, end1)
    ele%mat6(1:6, i) = (end1%vec - end0%vec) / del_orb%vec(i)
  enddo

! symplectify

  if (present(error)) call mat_symp_check (ele%mat6, error)
  if (ele%symplectify) call mat_symplectify (ele%mat6, ele%mat6)

end subroutine
