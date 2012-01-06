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

!

if (iu == -1) then
  iu = lunget()
  open (iu, file = 'track_of_photon_that_generated_an_error')
endif

write (iu, *) '*********************************************'
write (iu, '(2i8, f10.1)') photon%ix_photon, 0, photon%start%energy
write (iu, fm) photon%start%vec

do n = 1, photon%n_wall_hit
  hit => wall_hit(n)
  write (iu, *) '*********************************************'
  write (iu, '(2i8, f10.1)') photon%ix_photon, n, hit%before_reflect%energy
  write (iu, fm) hit%before_reflect%vec
  write (iu, '(3(24x, f24.16))') hit%after_reflect%vec(2:6:2)
  write (iu, '(3f18.12, 10x, 3f16.10)') hit%dw_perp, hit%cos_perp_in, hit%cos_perp_out, hit%reflectivity
enddo

if (iu_hit_file == -1) then
  close (iu)
  print *, 'Written file: track_of_photon_that_generated_an_error'
endif

if (iu_hit_file == -1) then
  open (iu, file = 'error_photon_start')
  write (iu, *) '&start'
  write (iu, *) '  ran_state%ix            =', sr3d_params%ran_state%ix
  write (iu, *) '  ran_state%iy            =', sr3d_params%ran_state%iy
  write (iu, *) '  ran_state%number_stored =', sr3d_params%ran_state%number_stored
  write (iu, *) '  ran_state%engine        =', sr3d_params%ran_state%engine
  write (iu, *)
  write (iu, *) '  p%vec     =', photon%start%vec
  write (iu, *) '  p%energy  =', photon%start%energy
  write (iu, *) '/'
  close (1)
  print *, 'Written file: error_photon_start'
endif

end subroutine print_hit_points

end module
