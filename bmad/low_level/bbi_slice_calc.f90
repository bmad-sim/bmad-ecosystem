!+
! Subroutine bbi_slice_calc (ele, n_slice, z_slice)
!
! Routine to compute the longitudinal positions of the slices of a beambeam element.
!
! Input:
!   ele        -- ele_struct: beambeam element
!   n_slice    -- integer: Number of slices
!
! Output:
!   z_slice(:) -- real(rp): Array of slice positions 1:n_slice.
!                   zero padded for indexes greater than n_slice
!-

subroutine bbi_slice_calc (ele, n_slice, z_slice)

use bmad_struct

implicit none

type (ele_struct) ele

integer :: i, n_slice
real(rp) z_slice(:), y
real(rp) :: z_norm

!

z_slice = 0  

if (n_slice == 1) then
  z_slice(1) = ele%value(z_offset$)
elseif (n_slice > 1) then
  do i = 1, n_slice
    y = (i - 0.5) / n_slice - 0.5
    z_norm = inverse(probability_funct, y, -5.0_rp, 5.0_rp, 1.0e-5_rp)
    z_slice(i) = ele%value(sig_z$) * z_norm + ele%value(z_offset$)
  enddo
else
  print *, 'ERROR IN BBI_SLICE_CALC: N_SLICE IS NEGATIVE:', n_slice
  if (global_com%exit_on_error) call err_exit
endif

end subroutine bbi_slice_calc

