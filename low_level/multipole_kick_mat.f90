!+
! Subroutine multipole_kick_mat (knl, tilt, ref_species, ele, orbit, factor, mat6)
!
! Subroutine to return the multipole kick components needed to
! construct the transfer matrix.
! This routine is not meant for general use.
!
! Input:
!   knl(0:)     -- Real(rp): Strength of multipoles
!   tilt(0:)    -- Real(rp): Tilt of multipoles
!   ref_species -- integer: Reference species.
!   ele         -- ele_struct: Lattice element containing multipoles.
!   orbit       -- Coord_struct: coordinates of particle around which the
!                    multipole kick matrix is computed.
!   factor      -- real(rp): Factor to scale knl by.
!
! Output:
!   mat6(6,6)   -- Real(rp): matrix with kick values at mat6(2:4:2, 1:3:2).
!                   The rest of the matrix is untouched.
!-

subroutine multipole_kick_mat (knl, tilt, ref_species, ele, orbit, factor, mat6)

use equal_mod, dummy => multipole_kick_mat

implicit none

type (coord_struct) orbit
type (ele_struct) ele

real(rp) mat6(6,6), kmat1(4,4), factor, charge_dir
real(rp) knl(0:), tilt(0:)

integer ref_species, n

!                        

mat6(2:4:2, 1:3:2) = 0
if (orbit%vec(1) == 0 .and. orbit%vec(3) == 0 .and. knl(1) == 0) return
charge_dir = orbit%direction*orbit%time_dir * ele%orientation * charge_to_mass_of(orbit%species) / charge_to_mass_of(ref_species)

do n = 1, ubound(knl, 1)
  if (knl(n) /= 0) then
    call mat4_multipole (knl(n)*charge_dir, tilt(n), n, orbit, kmat1)
    mat6(2:4:2, 1:3:2) = mat6(2:4:2, 1:3:2) + factor * kmat1(2:4:2, 1:3:2)
  endif
enddo

end subroutine multipole_kick_mat

