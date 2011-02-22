module synrad3d_output_mod

use synrad3d_utils

contains

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

subroutine print_hit_points (iu_hit_file, photon, wall_hit, fmt)

implicit none

type (sr3d_photon_track_struct), target :: photon
type (sr3d_photon_wall_hit_struct), pointer :: hit
type (sr3d_photon_wall_hit_struct), target :: wall_hit(0:)

integer iu, n, iu_hit_file

character(20) fm
character(*), optional :: fmt
!


fm = '(6f12.6)'
if (present(fmt)) fm = fmt

iu = iu_hit_file 
if (iu == 0) return

write (iu, *) '*********************************************'
write (iu, '(2i8, f10.1)') photon%ix_photon, 0, photon%start%energy
write (iu, fm) photon%start%vec

do n = 1, photon%n_wall_hit
  hit => wall_hit(n)
  write (iu, *) '*********************************************'
  write (iu, '(2i8, f10.1)') photon%ix_photon, n, hit%before_reflect%energy
  write (iu, fm) hit%before_reflect%vec
  write (iu, '(3(12x, f12.6))') hit%after_reflect%vec(2:6:2)
  write (iu, '(3f10.4, 10x, 2f12.6)') hit%dw_perp, hit%cos_perp, hit%reflectivity
enddo

end subroutine

end module
