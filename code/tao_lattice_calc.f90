!+
! Subroutine tao_lattice_calc (s)
!
! Routine to calculate the lattice functions.
! 
!
! Input:
!   s        -- tao_super_universe_struct
!
!  Output:
!   s        -- tao_super_universe_struct
!-

subroutine tao_lattice_calc (s)

use tao_mod

implicit none

type (tao_super_universe_struct) s

integer i

! Closed orbit and Twiss calculation.

do i = 1, size(s%u)
  call twiss_and_track (s%u(i)%model, s%u(I)%model_orb)
enddo

! Transfer info from %model to %data arrays.

call tao_load_data_array (s)

end subroutine
