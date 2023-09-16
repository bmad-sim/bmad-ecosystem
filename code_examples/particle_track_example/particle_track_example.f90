!+ 
! Program particle_track_example
!
! Example program to track a particle a number of turns through a lattice.
!
! Command line syntax:
!   <path-to-bin-dir>/bin/particle_track_example <param-file-name>
! Default:
!   <param-file-name> = "particle_track.init"
!
! Output:
!   track.dat
!-

program particle_track_example

use bmad

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_pointer_struct), allocatable, target :: eles(:)
type (ele_struct), pointer :: ele1, ele2
type (coord_struct) start_orbit
type (coord_struct), allocatable :: orbit(:)

integer :: ran_seed = 0, n_turn = 1, ix_branch = 0
integer j, it, ie, track_state, ix_end, n_loc, n_track, n
logical err_flag, end_write, convert_from_prime_coords, output_prime_coords

character(40) write_track_at, pad
character(100) param_file, lat_file, dat_file, fmt

namelist / params / lat_file, dat_file, n_turn, ix_branch, start_orbit, ran_seed, bmad_com, write_track_at, &
                  convert_from_prime_coords, output_prime_coords

!------------------------------------------
! Get the input parameters

param_file = 'particle_track.init'    ! Default parameter file name
if (command_argument_count() > 1) then
  call get_command_argument(1, param_file)
endif

convert_from_prime_coords = .false.
output_prime_coords = .false.
write_track_at = 'beginning'

open (1, file = param_file)
read (1, nml = params)       ! Fortran namelist read
close (1)
call upcase_string(write_track_at)

call ran_seed_put(ran_seed)

!------------------------------------------
! Parse the lattice.
! Also reread the input parameters in case bmad_com components are set in the lattice.
! We want bmad_com components set in the param file to supercede.

print *,"Using lattice: ", quote(lat_file)
call bmad_parser (lat_file, lat)


open (1, file = param_file)
read (1, nml = params)       ! Fortran namelist read
close (1)

!------------------------------------------
! Tracking init.
! Convert_from_prime_coords = T means take the start_orbit%vec coords as being (x, x', y, y', z, pz) and
! convert to standard Bmad (x, px, y, py, z, pz) coords.

branch => lat%branch(ix_branch) ! Lattice branch to track through

if (write_track_at == 'ALL') then
  end_write = .true.
else
  call lat_ele_locator (write_track_at, branch%lat, eles, n_loc, err_flag, ix_dflt_branch = ix_branch)
  if (n_loc == 0) then
    print *, 'Lattice element not found: ' // trim(write_track_at)
    stop
  endif
  if (n_loc > 1) then
    print *, 'Note: Multiple lattice elements match: ' // trim(write_track_at) // '. Will output at all selected.'
    stop
  endif
  end_write = (eles(1)%ele%ix_ele == 0) 
endif

call ran_seed_put (ran_seed)                               ! Ran used if radiation is on.
call reallocate_coord(orbit, lat, ix_branch)               ! Allocate orbit(0:) array.

call init_coord (orbit(0), start_orbit, branch%ele(0), downstream_end$)  ! Init orbit(0).
if (convert_from_prime_coords) call angle_to_canonical_coords(orbit(0), 'ZGOUBI')

!------------------------------------------
! And track.

pad = ''
n_track = branch%n_ele_track
n = maxval(len_trim(branch%ele(0:n_track)%name))

open(1, file = dat_file)
write (1, '(2a, 4x, a, 6(10x,a,4x), 5x, 3(11x,a))') '#  Turn ix_ele  Name', pad(1:n-1), 's', 'x ', 'px', 'y ', 'py', 'z ', 'pz', 'Sx', 'Sy', 'Sz'
fmt = '(2i7, 2x, a' // int_str(n) // ', f14.6, 6es16.8, 4x, 3f13.8, 3x, a)'

ix_end = branch%n_ele_track - 1

do it = 0, n_turn-1
  call track_all (lat, orbit, ix_branch, track_state, err_flag)   ! Track 1-turn.

  if (track_state /= moving_forward$) then
    print *, 'Particle lost on turn: ' // int_str(it) // ' At element: ' // int_str(track_state)
    ix_end = track_state
  endif

  if (write_track_at == 'ALL') then
    do ie = 0, ix_end
      write (1, fmt) it, ie, branch%ele(ie)%name, branch%ele(ie)%s, vec(orbit(ie)), orbit(ie)%spin
    enddo

  else
    do ie = 1, n_loc
      call find_element_ends (eles(ie)%ele, ele1, ele2)
      if (ele2%ix_ele > ix_end) cycle
      j = ele2%ix_ele
      write (1, fmt) it, j, ele2%name, ele2%s, vec(orbit(j)), orbit(j)%spin
    enddo
  endif

  orbit(0) = orbit(branch%n_ele_track)                            ! Set initial orbit = final.
  if (track_state /= moving_forward$) exit
enddo

if (it == n_turn .and. end_write) then
  write (1, fmt) n_turn, 0, branch%ele(0)%name, branch%ele(0)%s, vec(orbit(0)), orbit(0)%spin
endif

close(1)
print *, "Written: ", trim(dat_file)

!-----------------------------------------------
contains

function vec(orbit) result (vec_out)

type (coord_struct) orbit, orb2
real(rp) vec_out(6)

if (output_prime_coords) then
  orb2 = orbit
  call canonical_to_angle_coords(orb2, 'ZGOUBI')
  vec_out = orb2%vec
else
  vec_out = orbit%vec
endif

end function vec

end program
