!+
! Subroutine twiss_from_mat2 (mat, det, twiss, stat, tol, type_out)
!
! Subroutine to extract the twiss parameters from one-turn
! 2x2 matrices
!
! Modules needed:
!   use bmad
!
! Input:
!   mat(:, :)   -- Real(rdef): Input matrix
!   type_out    -- Logical: If .true. then an error message is typed out
!                            for a non ok$ STAT
!   tol         -- Real(rdef): tollerence for nonsymplectiy
!
! Output:
!   stat -- Integer: status of results:
!                         OK$, UNSTABLE$, NON_SYMPLECTIC$
!   det    -- Real(rdef): Determinate of matrix. Should be = 1
!   twiss  -- Twiss_struct: Twiss parameters
!                            TWISS.PHI is in radians, 0 < TWISS.PHI < twopi
!-

!$Id$
!$Log$
!Revision 1.3  2002/02/23 20:32:28  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:32:00  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine twiss_from_mat2 (mat, det, twiss, stat, tol, type_out)

  use bmad
  implicit none
  type (twiss_struct)  twiss

  integer stat
  real(rdef) mat(:, :), t_cos, t_sin, det, tol, radical
  logical type_out

!

  stat = ok$

  t_cos = (mat(1,1) + mat(2,2)) / 2.0
  det = mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)

  if (abs(t_cos) >= 1.0) then
    if (type_out) type *, 'ERROR IN TWISS_FROM_MAT: UNSTABLE MATRIX'
    stat = unstable$
    return
  endif

  if (abs(det - 1) > tol) then
    if (type_out) type *,  &
                    'WARNING IN TWISS_FROM_MAT: MATRIX DETERMINANT >< 1:', det
    stat = non_symplectic$
  endif

  t_sin = sign(1.0, mat(1,2)) * sqrt(1.0 - t_cos**2)
  twiss%phi = atan2(t_sin, t_cos)
  if (twiss%phi < 0) twiss%phi = twiss%phi + twopi

  if (abs(t_sin) < 1.0e-7) then
    twiss%alpha = 0.0
    twiss%beta = 1.0
    twiss%gamma = 1.0
  else
    twiss%alpha = (mat(1,1) - mat(2,2)) / (2.0 * t_sin)
    radical = -(1 + twiss%alpha**2) * mat(1,2) / mat(2,1)
    if (radical <= 0) then
      if (type_out) type *, 'ERROR IN TWISS_FROM_MAT: NON-SYMPLECTIC MATRIX'
      stat = non_symplectic$
      return
    endif
    twiss%beta = sqrt(radical)
    twiss%gamma = (1 + twiss%alpha**2) / twiss%beta
  endif

end subroutine
