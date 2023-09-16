!+
! Subroutine twiss_from_mat2 (mat, twiss, stat, type_out)
!
! Subroutine to extract the twiss parameters from one-turn
! 2x2 matrices
!
! Input:
!   mat(:, :)   -- Real(rp): Input matrix.
!   type_out    -- Logical: If .true. then an error message is typed out
!                            for a non ok$ STAT
!
! Output:
!   stat -- Integer: status of results:
!                         OK$, UNSTABLE$, NON_SYMPLECTIC$
!   twiss  -- Twiss_struct: Twiss parameters
!                            TWISS.PHI is in radians, 0 < TWISS.PHI < twopi
!-

subroutine twiss_from_mat2 (mat_in, twiss, stat, type_out)

use bmad_struct

implicit none

type (twiss_struct)  twiss

integer stat
real(rp) mat_in(:, :), t_cos, t_sin, det, radical, mat(2,2)
logical type_out
character(16), parameter :: r_name = 'twiss_from_mat2'

!

stat = ok$

det = mat_in(1,1) * mat_in(2,2) - mat_in(1,2) * mat_in(2,1)
mat = mat_in / det

t_cos = (mat(1,1) + mat(2,2)) / 2.0
if (abs(t_cos) >= 1.0) then
  if (type_out) call out_io (s_warn$, r_name, 'UNSTABLE MATRIX')
  stat = unstable$
  return
endif

t_sin = sign(1.0_rp, mat(1,2)) * sqrt(1.0 - t_cos**2)
twiss%phi = atan2(t_sin, t_cos)
if (twiss%phi < 0) twiss%phi = twiss%phi + twopi

twiss%alpha = (mat(1,1) - mat(2,2)) / (2.0 * t_sin)
radical = -(1 + twiss%alpha**2) * mat(1,2) / mat(2,1)
if (radical <= 0) then
  if (type_out) call out_io (s_warn$, r_name, 'NON-SYMPLECTIC MATRIX')
  stat = non_symplectic$
  return
endif
twiss%beta = sqrt(radical)
twiss%gamma = (1 + twiss%alpha**2) / twiss%beta

end subroutine twiss_from_mat2

