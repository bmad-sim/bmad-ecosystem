!+
! Subroutine tao_lattice_calc ()
!
! Routine to calculate the lattice functions.
! 
!
! Input:
!
!  Output:
!-

subroutine tao_lattice_calc ()

use tao_mod

implicit none


integer i

! Closed orbit and Twiss calculation.

do i = 1, size(s%u)
  call twiss_and_track (s%u(i)%model, s%u(i)%model_orb)
enddo

! Transfer info from %model to %data arrays.

call tao_load_data_array ()

end subroutine
