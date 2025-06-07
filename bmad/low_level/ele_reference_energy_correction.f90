!+
! Subroutine ele_reference_energy_correction (ele, orbit, particle_at, mat6, make_matrix)
!
! For elements where the reference energy is changing the reference energy in the 
! body is taken by convention to be the reference energy at the exit end.
! Elements where the reference energy can change:
!   lcavity
!   patch
!   custom
!
! This routine should be called at the start of any tracking integration.
!
! Input:
!   ele           -- ele_struct: Element being tracked through.
!   orbit         -- coord_struct: Coordinates to correct.
!   particle_at   -- integer: first_track_edge$ (that is, entering the element), or 
!                             second_track_edge$ (that is, leaving the element), or
!                             upstream_end$ (inherit ele%value(p0c_start$) ref), or
!                             downstream_end$ (inherit ele%value(p0c$)).
!                             inside$ (or anything else) -> Do nothing.
!   mat6(6,6)     -- real(rp), optional: Transfer matrix before correction.
!   make_matrix   -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit         -- coord_struct: Coordinates to correct.
!   mat6(6,6)     -- real(rp), optional: Transfer matrix transfer matrix including correction.
!-

subroutine ele_reference_energy_correction (ele, orbit, particle_at, mat6, make_matrix)

use bmad_struct

implicit none

type (ele_struct) :: ele
type (coord_struct) :: orbit

real(rp), optional :: mat6(6,6)

integer particle_at
logical, optional :: make_matrix

character(*), parameter :: r_name = 'ele_reference_energy_correction'

!

if (ele%value(p0c$) == ele%value(p0c_start$)) return

if ((orbit%direction*orbit%time_dir == 1 .and. particle_at == first_track_edge$) .or. particle_at == downstream_end$) then
  call orbit_reference_energy_correction(orbit, ele%value(p0c$), mat6, make_matrix)
elseif ((orbit%direction*orbit%time_dir == -1 .and. particle_at == second_track_edge$) .or. particle_at == upstream_end$) then
  call orbit_reference_energy_correction(orbit, ele%value(p0c_start$), mat6, make_matrix)
else
  return
endif

end subroutine ele_reference_energy_correction

