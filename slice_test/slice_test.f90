program slice_test

use bmad
use transfer_map_mod

implicit none

type (lat_struct) lat
type (coord_struct), allocatable :: ref_orb(:)
type (coord_struct) orb1, orb2a, orb2b
type (ele_struct) ele1, ele2a, ele2b

real(rp) xmat12(6,6)
real(rp) vec12(6)

!

call bmad_parser ('slice_test.bmad', lat)
call reallocate_coord (ref_orb, lat)
ref_orb%vec(2) = 0.4

call track_all (lat, ref_orb)
call lat_make_mat6 (lat, -1, ref_orb)
call twiss_propagate_all (lat)

call mat6_from_s_to_s (lat, xmat12, vec12, 0.5_rp, 2.5_rp)
call twiss_and_track_at_s (lat, 0.5_rp, ele1, ref_orb, orb1)
call twiss_and_track_at_s (lat, 2.5_rp, ele2a, ref_orb, orb2a)

call twiss_and_track_from_s_to_s (lat, 0.5_rp, 2.5_rp, .true., .true.,  &
                                                   orb1, orb2b, ele1, ele2b)

open (1, file = 'correct.now')

write (1, '(a, es20.10)') '"vec(1)" REL  1E-10', orb2a%vec(1)
write (1, '(a, es20.10)') '"vec(2)" REL  1E-10', orb2a%vec(2)
write (1, '(a, es20.10)') '"vec(3)" REL  1E-10', orb2a%vec(3)
write (1, '(a, es20.10)') '"vec(4)" REL  1E-10', orb2a%vec(4)
write (1, '(a, es20.10)') '"vec(5)" REL  1E-10', orb2a%vec(5)
write (1, '(a, es20.10)') '"vec(6)" REL  1E-10', orb2a%vec(6)
write (1, '(a, es20.10)') '"vec(6)" REL  1E-10', orb2a%vec(6)
write (1, '(a, es20.10)') '"t"      REL  1E-10', orb2a%t
write (1, '(a, es20.10)') '"s"      REL  1E-10', orb2a%s

write (1, *)

write (1, '(a, es20.10)') '"D:vec(1)" ABS  1E-10', orb2b%vec(1) - orb2a%vec(1)
write (1, '(a, es20.10)') '"D:vec(2)" ABS  1E-10', orb2b%vec(2) - orb2a%vec(2)
write (1, '(a, es20.10)') '"D:vec(3)" ABS  1E-10', orb2b%vec(3) - orb2a%vec(3)
write (1, '(a, es20.10)') '"D:vec(4)" ABS  1E-10', orb2b%vec(4) - orb2a%vec(4)
write (1, '(a, es20.10)') '"D:vec(5)" ABS  1E-10', orb2b%vec(5) - orb2a%vec(5)
write (1, '(a, es20.10)') '"D:vec(6)" ABS  1E-10', orb2b%vec(6) - orb2a%vec(6)
write (1, '(a, es20.10)') '"D:vec(6)" ABS  1E-10', orb2b%vec(6) - orb2a%vec(6)
write (1, '(a, es20.10)') '"D:t"      ABS  1E-10', orb2b%t - orb2a%t
write (1, '(a, es20.10)') '"D:s"      ABS  1E-10', orb2b%s - orb2a%s

write (1, *)

write (1, '(a, es20.10)') '"a%beta " REL  1E-10', ele2a%a%beta
write (1, '(a, es20.10)') '"b%beta " REL  1E-10', ele2b%b%beta
write (1, '(a, es20.10)') '"a%alpha" REL  1E-10', ele2a%a%alpha
write (1, '(a, es20.10)') '"b%alpha" REL  1E-10', ele2b%b%alpha
write (1, '(a, es20.10)') '"a%eta  " REL  1E-10', ele2a%a%eta
write (1, '(a, es20.10)') '"b%eta  " REL  1E-10', ele2b%b%eta

write (1, *)

write (1, '(a, es20.10)') '"D:a%beta"  ABS  1E-10', ele2b%a%beta - ele2a%a%beta
write (1, '(a, es20.10)') '"D:b%beta"  ABS  1E-10', ele2b%b%beta - ele2b%b%beta
write (1, '(a, es20.10)') '"D:a%alpha" ABS  1E-10', ele2b%a%alpha - ele2b%a%alpha
write (1, '(a, es20.10)') '"D:b%alpha" ABS  1E-10', ele2b%b%alpha - ele2b%b%alpha
write (1, '(a, es20.10)') '"D:a%eta"   ABS  1E-10', ele2b%a%eta - ele2b%a%eta
write (1, '(a, es20.10)') '"D:b%eta"   ABS  1E-10', ele2b%b%eta - ele2b%b%eta


close (1)


end program 
