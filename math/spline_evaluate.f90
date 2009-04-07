!+
! Subroutine spline_evaluate (spline, x, ok, y, dy)
!
! Subroutine to evalueate a spline at a set of points.
!
! Modules used:
!   use sim_utils
!
! Input:
!   spline(:) -- Spline_struct: Spline structure.
!   x         -- Real(rp): point for evaluation.
!
! Output:
!   ok        -- Logical: Set .true. if everything ok
!   y         -- Real(rp), optional: Spline interpolation.
!   dy        -- Real(rp), optional: Spline derivative interpolation.
!
! Note:
!   The point x must lie between spline(1)%x and spline(max)%x
!-

#include "CESR_platform.inc"

subroutine spline_evaluate (spline, x, ok, y, dy)

  use sim_utils_struct
  use sim_utils_interface, except => spline_evaluate
  use sim_utils

  implicit none

  type (spline_struct), target :: spline(:)

  real(rp) :: x
  real(rp), optional :: y, dy
  real(rp) :: c(0:3)

  real(rp) dx       

  integer j1, j2, j3, smax
                    
  logical ok       
  character(16) :: r_name = 'spline_evaluate'

! Check if x value out of bounds.
            
  ok = .false.

  smax = ubound(spline(:), 1)
  dx = 1e-6 * (spline(smax)%x - spline(1)%x)   ! something small

  if (x < spline(1)%x - dx) then
    call out_io (s_error$, r_name, &
                'X EVALUATION POINT LESS THAN LOWER BOUND OF SPLINE INTERVAL')
    return
  endif
                                
  if (x > spline(smax)%x + dx) then
    call out_io (s_error$, r_name, &
             'X EVALUATION POINT GREATER THAN UPPER BOUND OF SPLINE INTERVAL')
    return
  endif

! find correct interval

  j1 = 1
  j2 = smax

  do while (.true.)
    j3 = (j1 + j2) / 2
    if (x > spline(j3)%x) then
      j1 = j3
    else
      j2 = j3
    endif
    if (j2 == j1 + 1) exit
  enddo

! evaluate

  dx = x - spline(j1)%x
  c = spline(j1)%coef

  if (present(y)) then
    y = (((c(3) * dx) + c(2)) * dx + c(1)) * dx + c(0)
  endif

  if (present(dy)) then
   dy = ((3*c(3) * dx) + 2*c(2)) * dx + c(1)
  endif

  ok = .true.
  return

end subroutine
