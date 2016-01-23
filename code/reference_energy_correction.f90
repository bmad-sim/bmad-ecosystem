!+
! Subroutine reference_energy_correction (ele, orbit, particle_at)
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
!   ele         -- Ele_struct: Element being tracked through.
!   orbit       -- Coord_struct: Coordinates to correct.
!   particle_at -- integer: first_track_edge$ (that is, entering the element), or 
!                           second_track_edge$ (that is, leaving the element).
!
! Output:
!   orbit     -- Coord_struct: Coordinates to correct.
!-

subroutine reference_energy_correction (ele, orbit, particle_at)

use bmad_struct

implicit none

type (ele_struct) :: ele
type (coord_struct) :: orbit

real(rp) p0, p1, e_start, p_rel
integer particle_at
character(*), parameter :: r_name = 'reference_energy_correction'

!

if (ele%value(p0c$) == ele%value(p0c_start$)) return

if (orbit%direction == 1 .and. particle_at == first_track_edge$) then
  p_rel = ele%value(p0c_start$) / ele%value(p0c$)
  orbit%p0c = ele%value(p0c$)
elseif (orbit%direction == -1 .and. particle_at == second_track_edge$) then
  p_rel = ele%value(p0c$) / ele%value(p0c_start$)
  orbit%p0c = ele%value(p0c_start$)
else
  return
endif

orbit%vec(2) = orbit%vec(2) * p_rel
orbit%vec(4) = orbit%vec(4) * p_rel
orbit%vec(6) = (1 + orbit%vec(6)) * p_rel - 1

end subroutine reference_energy_correction

