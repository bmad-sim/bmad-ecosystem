!+
! Function mat_symp_conj (mat) result (mat_conj)
!
! Function to take the symplectic conjugate of a square matrix.
!
! Modules needed:
!   use sim_utils
!
! Input:
!   mat(:, :) -- Real(rp): Input matrix.
!
! Output:
!   mat_conj(:, :) -- Real(rp): Symplectic conjugate of mat.
!-

function mat_symp_conj(mat) result (mat_conj)

use output_mod, only: rp, out_io, s_fatal$, global_com

implicit none

integer i, j, nn

real(rp) mat(:,:)
real(rp) :: mat_conj(size(mat, 1), size(mat, 1))

character(*), parameter :: r_name = 'mat_symp_conj'

! Check bounds

nn = size(mat, 1)

if (mod(nn, 2) /= 0 .or. nn /= size(mat, 2)) then
  call out_io (s_fatal$, r_name, 'ARRAY SIZE IS NOT EVEN!')
  if (global_com%exit_on_error) call err_exit
endif

! Compute conjugate

do i = 1, nn, 2
  do j = 1, nn, 2
    mat_conj(i,   j)   =  mat(j+1, i+1)
    mat_conj(i,   j+1) = -mat(j,   i+1)
    mat_conj(i+1, j)   = -mat(j+1, i)
    mat_conj(i+1, j+1) =  mat(j,   i)
  enddo
enddo

end function
