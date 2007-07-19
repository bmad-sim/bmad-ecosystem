!+         
! Subroutine transfer_matrix_calc (lat, rf_on, xfer_mat, xfer_vec, ix1, ix2)
!
! Subroutine to calculate the 6X6 transfer matrix between two elements.
! To calculate transfer maps see the routine: transfer_map_calc.
!
! The transfer matrix is from the end of element ix1 to the end of element ix2.
! If ix1 and ix2 are not present, the full 1-turn matrix is calculated.
!
! If ix2 < ix1 and lat%param%lattice_type is circular_lattice$ then the
! calculation will "wrap around" the lattice end.
! For example, if ix1 = 900 and ix2 = 10 then the xfer_mat is the matrix from
! element 900 to the lattice end plus from 0 through 10.
!
! If ix2 < ix1 and lat%param%lattice_type is linear_lattice$ then the backwards
! transfer matrix is computed.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat   -- lat_struct: Lattice used in the calculation.
!     %ele(:)%mat6  -- Transfer matrices used in the calculation.
!   rf_on -- Logical: Keep the RF on for the calculation? 
!              True does nothing. False forces the RF off.
!              False is what is needed for the Twiss parameter calculation.
!   ix1   -- Integer, optional: Element start index for the calculation.
!              Default is 0.
!   ix2   -- Integer, optional: Element end index for the calculation.
!              Default is lat%n_ele_track.
!
! Output:
!    xfr_mat(6,6) -- Real(rp): Transfer matrix.
!    xfr_vec(6)   -- Real(rp), optional: 0th order part of the transfer map.
!-

#include "CESR_platform.inc"

subroutine transfer_matrix_calc (lat, rf_on, xfer_mat, xfer_vec, ix1, ix2)

use bmad_struct
use bmad_interface, except_dummy => transfer_matrix_calc
use cesr_utils, only: integer_option

implicit none

type (lat_struct)  lat

real(rp), intent(out) :: xfer_mat(:,:)
real(rp), intent(out), optional :: xfer_vec(:)
real(rp) rf_mat(6,6)

integer, intent(in), optional :: ix1, ix2
integer i, i1, i2

logical, intent(in) :: rf_on
logical vec_present

!

vec_present = present(xfer_vec)
if (vec_present) xfer_vec = 0

call mat_make_unit (xfer_mat)
  
i1 = integer_option(0, ix1) 
i2 = integer_option(lat%n_ele_track, ix2) 

! Normal case

if (i1 <= i2) then  
  do i = i1+1, i2
    call add_on_to_xfer_mat
  enddo

! For a circular lattice push through the origin.

elseif (lat%param%lattice_type == circular_lattice$) then
  do i = i1+1, lat%n_ele_track
    call add_on_to_xfer_mat
  enddo
  do i = 1, i2
    call add_on_to_xfer_mat
  enddo

! For a linear lattice compute the backwards matrix

else
  do i = i2+1, i1 
    call add_on_to_xfer_mat
  enddo
  call mat_inverse (xfer_mat, xfer_mat)
  xfer_vec = -matmul(xfer_mat, xfer_vec)

endif

!--------------------------------------------------------
contains

subroutine add_on_to_xfer_mat

if (lat%ele(i)%key == rfcavity$ .and. .not. rf_on) then
  rf_mat = lat%ele(i)%mat6
  rf_mat(6,5) = 0  ! turn rf off
  xfer_mat = matmul (rf_mat, xfer_mat)
  if (vec_present) xfer_vec = matmul(rf_mat, xfer_vec) + lat%ele(i)%vec0
else
  xfer_mat = matmul (lat%ele(i)%mat6, xfer_mat)
  if (vec_present) xfer_vec = &
          matmul(lat%ele(i)%mat6, xfer_vec) + lat%ele(i)%vec0
endif

end subroutine

end subroutine
                                          
