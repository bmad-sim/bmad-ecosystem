module synrad3d_output_mod

use synrad3d_utils

implicit none

contains

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

subroutine sr3d_print_hit_points (iu_hit_file, photon, wall_hit, branch, lots_of_digits)

type (sr3d_photon_track_struct), target :: photon
type (sr3d_photon_wall_hit_struct), pointer :: hit
type (sr3d_photon_wall_hit_struct), target :: wall_hit(0:)
type (branch_struct) branch

integer iu, n, iu_hit_file

logical, optional :: lots_of_digits
character(100) fm, fm2
character(40) wall_name

!

if (logic_option(.false., lots_of_digits)) then
  fm  = '(6es25.15)'
  fm2 = '(i7, i4, f10.2, 5x, 3es25.15, i4, 2(5x, 3es25.15), 10x, 3f18.12, 5x, 3f16.10, 3x, a)' 
else
  fm  = '(6f12.6)'
  fm2 = '(i7, i4, f10.2, 5x, 2f10.6, f14.6, i4, 2(5x, 3f10.6), 10x, 3f10.6, 5x, 3f10.6, 3x, a)'
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
  if (n == 0) then
    wall_name = '<inside chamber>'
  else
    wall_name = branch%wall3d(hit%ix_wall3d)%name
    if (wall_name == '') wall_name = '<default_subchamber>'
  endif
  write (iu, fm2) photon%ix_photon, n, hit%before_reflect%p0c, hit%after_reflect%vec(1:3:2), hit%after_reflect%s, &
          hit%ix_branch, hit%before_reflect%vec(2:6:2), hit%after_reflect%vec(2:6:2), &
          hit%dw_perp, hit%cos_perp_in, hit%cos_perp_out, hit%reflectivity, trim(wall_name)
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
  write (iu, *) '  ix_branch   =', photon%start%ix_branch
  write (iu, *) '/'
  close (iu)
  print *, 'Written file: error_photon_start'
endif

end subroutine sr3d_print_hit_points

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

subroutine sr3d_record_photon_position (action, photon)

type (sr3d_photon_track_struct), target, optional :: photon
type (coord_struct), pointer :: orb

integer ios
integer, save :: iu = 0
character(*) action
character(120) line

! Only tracks that pass the filter tests get into the official track file.

select case (action)

case ('START_RECORDING')
  ! Open a scratch file
  if (iu /= 0) call err_exit ! Should be zero
  iu = lunget()
  open (iu, status = 'scratch')

case ('ERASE_RECORDING')
  ! Close scratch file
  if (iu == 0) call err_exit ! Should be non-zero
  close (iu)
  iu = 0

case ('MOVE_TRACK_TO_FILE')
  ! Move info from scratch file to the official track file
  if (iu == 0) call err_exit ! Should be non-zero
  rewind (iu)
  do
    read (iu, '(a)', iostat = ios) line
    if (ios /= 0) exit
    write (sr3d_params%iu_photon_track, '(a)') trim(line)
  enddo
  close (iu)
  iu = 0

case ('RECORD_TRACK_POINT')
  ! Record a track point in the scratch file.
  if (iu == 0) call err_exit ! Should be non-zero
  orb => photon%now%orb
  write (iu, '(i8, i10, 2f11.6, f13.6, i4, 5x, 3f11.6)') &
     photon%ix_photon, photon%ix_photon_generated, orb%vec(1:3:2), orb%s, photon%now%ix_branch, orb%vec(2:6:2)

case default
  call err_exit   ! Should not be here
end select

end subroutine sr3d_record_photon_position

end module
