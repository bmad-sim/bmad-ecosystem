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
use tao_data_mod

implicit none


integer i

! make sure useit is up-to-date

call tao_set_var_useit_opt
call tao_set_data_useit_opt

! Closed orbit and Twiss calculation.

do i = 1, size(s%u)
  call twiss_and_track (s%u(i)%model, s%u(i)%model_orb)
enddo

! Transfer info from %model to %data arrays.

call tao_load_data_array ()

end subroutine
