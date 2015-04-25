! +
! Program sad_to_bmad_postprocess
!
! Program to process a bmad lattice file after it has been created by the 
! sad_to_bmad.py script. 
!
! Syntax:
!   sad_to_bmad_postprocess <bmad-lattice-file>
!-

program sad_to_bmad_postprocess

use bmad

implicit none

type (lat_struct) lat
type (coord_struct), allocatable :: orbit(:)

real(rp) f_shift

integer ios, ix

logical ok, fshift_found

character(100) lat_file, temp_file
character(200) line

! Read in lattice file and write out temorary file

if (cesr_iargc() /= 1) then
  print *, 'Command line syntax:'
  print *, '  sad_to_bmad_postprocess <lattice-file-name>'
  stop
endif

call cesr_getarg (1, lat_file)
temp_file = 'temp.bmad'

open (1, file = lat_file)
open (2, file = temp_file)

write (2, '(a)') 'no_digested'
write (2, '(a)') 'fshift = 0'
do
  read (1, '(a)', iostat = ios) line
  if (ios /= 0) exit
  if (line(1:2) == '++') cycle   ! '++' => lines to ignore
  write (2, '(a)') trim(line)
enddo

close(1)
close(2)

! Parse temp file and calculate fshift

call bmad_parser (temp_file, lat)
call set_on_off (rfcavity$, lat, off$)
call reallocate_coord (orbit, lat)
call twiss_and_track (lat, orbit, ok)

f_shift = orbit(lat%n_ele_track)%vec(5) / lat%param%total_length

! Create file with f_shift value

open (1, file = lat_file)
open (2, file = temp_file)

fshift_found = .false.
do
  read (1, '(a)', iostat = ios) line
  if (ios /= 0) exit
  if (line(1:8) == '++fshift') then
    fshift_found = .true.
    ix = index(line, '!')
    write (2, '(a, es20.8, 2x, a)') 'fshift =', f_shift, trim(line(ix:))
  endif
  if (line(1:2) == '++') cycle
  write (2, '(a)') trim(line)
enddo

close(1)
close(2)

if (.not. fshift_found) then
  print *, 'fshift line not found! Something is wrong.'
  stop
endif

! And now move temp file to lattice file

call system_command ('mv ' // trim(temp_file) // ' ' // trim(lat_file))

end program
