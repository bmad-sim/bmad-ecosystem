module tao_cut_ring_mod

use tao_struct
use tao_interface
use tao_utils
use tao_lattice_calc_mod

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_cut_ring ()
!
! Routine to cut (make non-closed) a circular lattice.
!-

subroutine tao_cut_ring ()

type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: lat

logical ok

!

u => tao_pointer_to_universe(-1)
lat => u%model%lat

lat%param%lattice_type = linear_lattice$
u%calc%lattice = .true.
u%model%lat%beam_start%vec = 0
call tao_lattice_calc (ok)

end subroutine tao_cut_ring

end module
