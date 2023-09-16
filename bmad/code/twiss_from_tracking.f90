!+
! Subroutine twiss_from_tracking (lat, ref_orb0, symp_err, err_flag, d_orb)
!
! Subroutine to compute, from tracking, and for every element in the lat, 
! the Twiss parameters and the transfer matrices about some reference orbit.
!
! The is done by tracking 12 particle offset from the reference orbit 
! by +/- d_orb.
!
! Because of nonlinearities, and the fact that d_orb is finite,
! the computed 6x6 matrices will not be exactly symplectic. To handle this
! these matrices are symplecitified.
!     
! See twiss_at_start and twiss_propagate_all for alternative routines.
!
! Input:
!   lat            -- lat_struct: Lat to track through.
!   ref_orb0        -- Coord_struct: Reference orbit at lat%ele(0).
!   d_orb(6)        -- real(rp), optional: Vector of offsets to use. 
!                       If not present or zero bmad_com%d_orb(:) will be used.
!
! Output:
!   lat         -- lat_struct: Structure holding the Twiss parameters.
!     %mat6        -- transfer matrices.
!
!   symp_err    -- Real(rp): A measure of how symplectic the constructed 
!                   matrices were before symplecitification. 
!                  mat_symp_check for more details.
!   err_flag    -- Logical: Set True if there is an error. False otherwise.
!-

subroutine twiss_from_tracking (lat, ref_orb0, symp_err, err_flag, d_orb)

use bmad_interface, except_dummy => twiss_from_tracking
use bookkeeper_mod, only: set_on_off, restore_state$, off_and_save$

implicit none

type (lat_struct), intent(inout) :: lat
type (coord_struct), intent(in) :: ref_orb0

real(rp), optional, intent(in) :: d_orb(:)
real(rp), intent(out) :: symp_err
real(rp) delta(0:6), r(6), r0(6), r1(6)
real(rp), allocatable :: on_off_save(:)

type multi_orb_struct
  type (coord_struct), allocatable :: orb(:) 
end type
type (multi_orb_struct) :: mo(-6:6)

real(rp) mat(6,6), mat1(6,6), mat0(6,6), mat_inv(6,6), mat6_unit(6,6)
integer i, j, n_ele, status
logical err_flag, err

! Init

err_flag = .true.

do i = -6, 6
  call reallocate_coord (mo(i)%orb, lat%n_ele_max)
enddo

call mat_make_unit (mat6_unit)
n_ele = lat%n_ele_track

! 

call set_on_off (rfcavity$, lat, off_and_save$, saved_values = on_off_save)

delta(0) = 0
delta(1:6) = bmad_com%d_orb(1:6)
if (present(d_orb)) where (d_orb(1:6) /= 0) delta(1:6) = d_orb(1:6)

! track offset particles

mo(0)%orb(0) = ref_orb0 
call track_all (lat, mo(0)%orb, err_flag = err)
if (err) goto 9000
r = mo(0)%orb(n_ele)%vec - mo(0)%orb(0)%vec

do i = 1, 6
  mo(i)%orb(0) = ref_orb0 
  mo(i)%orb(0)%vec(i) = ref_orb0%vec(i) + delta(i)
  call track_all (lat, mo(i)%orb, err_flag = err)
  if (err) goto 9000
  mo(-i)%orb(0) = ref_orb0 
  mo(-i)%orb(0)%vec(i) = ref_orb0%vec(i) - delta(i)
  call track_all (lat, mo(-i)%orb, err_flag = err)
  if (err) goto 9000
  mat(:,i) = (mo(i)%orb(n_ele)%vec - mo(-i)%orb(n_ele)%vec) / (2*delta(i))
enddo

! compute the element transfer matrices.
! mat0 is the transfer matrix from the start to the beginning of element #j
! mat  is the transfer matrix from the start to the end of element #j
! mat1 is the transfer matrix through element #j

call mat_make_unit (mat0)  
r0 = 0

do j = 1, n_ele

  do i = 1, 6
    mat(:,i) = (mo(i)%orb(j)%vec - mo(-i)%orb(j)%vec) / (2*delta(i))
  enddo

  symp_err = mat_symp_error(mat)
  call mat_symplectify (mat, mat)

  r = mo(0)%orb(j)%vec - matmul(mat, mo(0)%orb(0)%vec)

  mat_inv = mat_symp_conj (mat0)   ! symp_conj is the inverse
  mat1 = matmul (mat, mat_inv)

  lat%ele(j)%mat6 = mat1
  lat%ele(j)%vec0 = r - matmul(mat1, r0)

  r0 = r
  mat0 = mat

enddo

! Now compute the twiss parameters.
! And turn the RF back on.

call twiss_at_start (lat, status)
if (status /= ok$) goto 9000

call twiss_propagate_all (lat, err_flag = err)
if (err) goto 9000

err_flag = .false.

! And reset

9000 continue

call set_on_off (rfcavity$, lat, restore_state$, saved_values = on_off_save)

end subroutine
