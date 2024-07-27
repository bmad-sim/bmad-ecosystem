!+
! Function mat_symp_conj_i (mat) result (mat_conj)
!
! Function to take the symplectic conjugate of a square complex matrix.
!
! Input:
!   mat(:, :) -- complex(rp): Input matrix.
!
! Output:
!   mat_conj(:, :) -- complex(rp): Symplectic conjugate of mat.
!-

function mat_symp_conj_i(mat) result (mat_conj)

use output_mod, only: rp, out_io, s_error$, global_com

implicit none

integer i, j, nn, dims

complex(rp) mat(:,:)
complex(rp) :: mat_conj(size(mat, 1), size(mat, 1))
real(rp), allocatable :: S(:,:)

character(*), parameter :: r_name = 'mat_symp_conj'

! Check bounds

nn = size(mat, 1)

if (mod(nn, 2) /= 0 .or. nn /= size(mat, 2)) then
  call out_io (s_error$, r_name, 'ARRAY SIZE IS NOT EVEN!')
  if (global_com%exit_on_error) call err_exit
endif

dims = nn/2
allocate(S(nn,nn))
S=0.0d0
do i=1,nn,2
  S(i,i+1) =  1.0d0
  S(i+1,i) = -1.0d0
enddo

mat_conj = -1.0d0*matmul(S,matmul(transpose(mat),S))

deallocate(S)

end function
