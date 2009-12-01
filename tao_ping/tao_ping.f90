program tao_ping

use tao_ping_mod

implicit none

type (ping_universe_struct) pu
type (lat_struct) lat

! Init

call read_parameters(pu)

call read_data (pu)
print *, 'Number of bpms: ', size(pu%bpm)
print *, 'Number of turns:', size(pu%bpm(1)%x)

call bmad_parser (pu%param%lattice_file, lat)
call twiss_at_start (lat)
call twiss_propagate_all (lat)

! FFT data



! Fit data

! Plot data

! Print results


end program

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine read_parameters (pu)

use tao_ping_mod

implicit none

type (ping_universe_struct) :: pu
type (ping_param_struct) param

namelist / ping_params / param

!

open (1, file = 'tao_ping.init', status = 'old')
read (1, nml = ping_params)
close (1)

pu%param = param


end subroutine

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine read_data (pu)

use tao_ping_mod

implicit none

type (ping_universe_struct) :: pu

integer ix, ios, n_turn, n_turn_now, n_bpm, nt

character(200) line

! find number of BPM's and number of turns of data

open (1, file = pu%param%data_file, status = 'old')

n_turn = -1
n_bpm = 0

do

  ! skip header

  do 
    read (1, '(a)', iostat = ios) line
    if (ios /= 0) exit
    if (line == '') cycle
    if (line(1:1) /= '#') exit
  enddo

  if (ios /= 0) exit

  ! Count turn-by-turn data lines

  n_turn_now = 0
  do
    n_turn_now = n_turn_now + 1
    read (1, '(a)', iostat = ios) line
    if (ios /= 0) exit
    if (line == '') exit
    if (line(1:1) == '#') exit
  enddo

  if (n_turn == -1) then
    n_turn = n_turn_now
  elseif (n_turn /= n_turn_now) then
    print *, 'Error: inconsistant number of turns: ', n_turn, n_turn_now
    call err_exit
  endif

  n_bpm = n_bpm + 1
  if (ios /= 0) exit

enddo

! Now load the data

allocate (pu%bpm(n_bpm))

rewind(1)

n_bpm = 0
do

  n_bpm = n_bpm + 1
  allocate (pu%bpm(n_bpm)%x(n_turn), pu%bpm(n_bpm)%y(n_turn))

  ! Skip header

  do 
    read (1, '(a)', iostat = ios) line
    if (ios /= 0) exit
    if (line == '') cycle
    if (line(1:10) == '# Location') then
      ix = index(line, ':')
      call string_trim (line(ix+1:), pu%bpm(n_bpm)%name, ix)
    endif
    if (line(1:1) /= '#') exit
  enddo

  if (ios /= 0) exit

  ! Read turn-by-turn data

  nt = 0
  do
    nt = nt + 1
    call string_trim (line, line, ix)
    call string_trim (line(ix+1:), line, ix)  ! Trim button 1
    call string_trim (line(ix+1:), line, ix)  ! Trim button 2
    call string_trim (line(ix+1:), line, ix)  ! Trim button 3
    call string_trim (line(ix+1:), line, ix)  ! Trim button 4
  
    read (line(1:ix), *) pu%bpm(n_bpm)%x(nt)
    call string_trim (line(ix+1:), line, ix)  ! Trim x

    read (line(1:ix), *) pu%bpm(n_bpm)%y(nt)

    read (1, '(a)', iostat = ios) line
    if (ios /= 0) exit
    if (line == '') exit
    if (line(1:1) == '#') exit
  enddo

  if (ios /= 0) exit
  if (n_bpm == size(pu%bpm)) exit  ! So blank lines at end do not cause problems.

enddo

close(1)

end subroutine
