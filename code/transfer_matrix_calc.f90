!+         
! Subroutine transfer_matrix_calc (lat, rf_on, t_mat, ix1, ix2)
!
! Subroutine to calculate the 6X6 transfer matrix between two elements.
!
! The transfer matrix is from the end of element ix1 to the end of element ix2.
! If ix1 and ix2 are not present, the full 1-turn matrix is calculated.
! If ix2 < ix1 then the calculation will "wrap around" the lattice end.
! For example if ix1 = 900 and ix2 = 10 then the t_mat is the matrix from
! element 900 to the lattice end plus from 0 through 10.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat   -- Ring_struct: Lattice used in the calculation.
!     %ele_(:)%mat6  -- Transfer matrices used in the calculation.
!   rf_on -- Logical: Keep the RF on for the calculation? 
!              True does nothing. False forces the RF off.
!              False is what is needed for the Twiss parameter calculation.
!   ix1   -- Integer, optional: Element start index for the calculation.
!              Default is 0.
!   ix2   -- Integer, optional: Element end index for the calculation.
!              Default is lat%n_ele_use.
!
! Output:
!    t_mat(6,6) -- Real(rp): Transfer matrix.
!-

#include "CESR_platform.inc"

subroutine transfer_matrix_calc (lat, rf_on, t_mat, ix1, ix2)

  use bmad_struct
  use bmad_interface, except => transfer_matrix_calc
  use cesr_utils, only: integer_option

  implicit none

  type (ring_struct)  lat

  real(rp), intent(out) :: t_mat(:,:)
  real(rp) rf_mat(6,6)

  integer, intent(in), optional :: ix1, ix2
  integer i, i1, i2

  logical, intent(in) :: rf_on

!

  call mat_make_unit (t_mat)

  i1 = integer_option(0, ix1) 
  i2 = integer_option(lat%n_ele_use, ix2) 

  if (i2 < i1) then
    do i = i1+1, lat%n_ele_use
      call add_on_to_t_mat
    enddo
    do i = 1, i2
      call add_on_to_t_mat
    enddo

  else
    do i = i1+1, i2
      call add_on_to_t_mat
    enddo
  endif

!--------------------------------------------------------
contains

subroutine add_on_to_t_mat

  if (lat%ele_(i)%key == rfcavity$ .and. .not. rf_on) then
    rf_mat = lat%ele_(i)%mat6
    rf_mat(6,5) = 0  ! turn rf off
    t_mat = matmul (rf_mat, t_mat)
  else
    t_mat = matmul (lat%ele_(i)%mat6, t_mat)
  endif

end subroutine

end subroutine
                                          
