!+
! Subroutine tao_lattice_calc ()
!
! Routine to calculate the lattice functions.
! 
!
! Input:
!
! Output:
!-

subroutine tao_lattice_calc ()

use tao_mod
use tao_data_mod

implicit none

integer i

character(20) :: r_name = "tao_lattice_calc"

! make sure useit is up-to-date

call tao_set_var_useit_opt
call tao_set_data_useit_opt

! if using custom tracking then s%global%lattice_recalc needs to be set to FALSE
! in tao_hook_lattice_calc
call tao_hook_lattice_calc ()

! Closed orbit and Twiss calculation.
! This can be slow for large lattices so only do it if the lattice changed.
if (s%global%lattice_recalc) then
  if (s%global%track_type .eq. singleparticle$) then
    do i = 1, size(s%u)
      call twiss_and_track (s%u(i)%model, s%u(i)%model_orb)
    enddo
    s%global%lattice_recalc = .false.
  elseif (s%global%track_type .eq. macroparticle$) then
    do i = 1, size(s%u)
      call macro_track (s%u(i)%model, s%u(i)%model_orb)
    enddo
  else
    call out_io (s_fatal$, r_name, &
                   "This tracking type has yet to be implemented!")
    call out_io (s_blank$, r_name, &
                   "No tracking or twiss calculations will be perfomred.")
  endif
endif

! Transfer info from %model to %data arrays.
call tao_load_data_array ()

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
contains

subroutine macro_track (lattice, orbit)

implicit none

type (ring_struct) :: lattice
type (coord_struct) :: orbit(:)

end subroutine macro_track

end subroutine tao_lattice_calc
