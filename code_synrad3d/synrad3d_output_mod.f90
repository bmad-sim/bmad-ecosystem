module synrad3d_output_mod

use synrad3d_utils

contains

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

subroutine print_hit_points (iu_hit_file, photon, wall_hit, lots_of_digits)

implicit none

type (sr3d_photon_track_struct), target :: photon
type (sr3d_photon_wall_hit_struct), pointer :: hit
type (sr3d_photon_wall_hit_struct), target :: wall_hit(0:)

integer iu, n, iu_hit_file

logical, optional :: lots_of_digits
character(100) fm, fm2

!

if (logic_option(.false., lots_of_digits)) then
  fm  = '(6es25.15)'
  fm2 = '(2i8, f10.1, 5x, 3es25.15, 2(6x, 3es25.15), 10x, 3f18.12, 5x, 3f16.10)' 
else
  fm  = '(6f12.6)'
  fm2 = '(2i8, f10.1, 5x, 2f12.6, f14.6, 2(6x, 3f12.6), 10x, 3f12.6, 5x, 3f12.6)'
endif


iu = iu_hit_file 
if (iu == 0) return

!

if (iu == -1) then
  iu = lunget()
  open (iu, file = 'track_of_photon_that_generated_an_error')
endif

do n = 0, photon%n_wall_hit
  hit => wall_hit(n)
  write (iu, fm2) photon%ix_photon, n, hit%before_reflect%p0c, hit%after_reflect%vec(1:3:2), hit%after_reflect%s, &
          hit%before_reflect%vec(2:6:2), hit%after_reflect%vec(2:6:2), &
          hit%dw_perp, hit%cos_perp_in, hit%cos_perp_out, hit%reflectivity
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
  write (iu, *) '  p%vec  =', photon%start%orb%vec
  write (iu, *) '  p%p0c  =', photon%start%orb%p0c
  write (iu, *) '  p%direction =', photon%start%orb%direction
  write (iu, *) '/'
  close (iu)
  print *, 'Written file: error_photon_start'
endif

end subroutine print_hit_points

end module
