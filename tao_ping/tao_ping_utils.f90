module tao_ping_utils

use tao_ping_struct

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine ping_read_parameters ()

implicit none

type (ping_param_struct) param

namelist / ping_params / param

!

open (1, file = 'tao_ping.init', status = 'old')
read (1, nml = ping_params)
close (1)

ping_s%param = param


end subroutine

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine tao_hook_init_data ()

use tao_init_mod

implicit none

integer ix, ios, n_turn, n_turn_now, n_bpm, nt

character(200) line

! find number of BPM's and number of turns of data

open (1, file = ping_s%param%data_file, status = 'old')

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

rewind(1)

!---------------------------------------------------------------------
! Now load the data

call tao_init_data_in_universe (s%u(1), n_bpm)

n_bpm = 0
do

  n_bpm = n_bpm + 1

  ! Skip header

  do 
    read (1, '(a)', iostat = ios) line
    if (ios /= 0) exit
    if (line == '') cycle
    if (line(1:10) == '# Location') then
      ix = index(line, ':')
      call string_trim (line(ix+1:), ping_s%bpm(n_bpm)%name, ix)
    endif
    if (line(1:1) /= '#') exit
  enddo

  if (ios /= 0) exit

  call tao_d2_data_stuffit (s%u(1), 'orb_' // ping_s%bpm(n_bpm)%name, 2)
  do i = 1, 2
    n1 = u%n_data_used + 1
    n2 = u%n_data_used + n_turn
    u%n_data_used = n2
    if (n2 > size(u%data)) call tao_allocate_data_array (u, n2)



  ! Read turn-by-turn data

  nt = 0
  do
    nt = nt + 1
    call string_trim (line, line, ix)
    call string_trim (line(ix+1:), line, ix)  ! Trim button 1
    call string_trim (line(ix+1:), line, ix)  ! Trim button 2
    call string_trim (line(ix+1:), line, ix)  ! Trim button 3
    call string_trim (line(ix+1:), line, ix)  ! Trim button 4
  
    read (line(1:ix), *) ping_s%bpm(n_bpm)%x(nt)
    call string_trim (line(ix+1:), line, ix)  ! Trim x

    read (line(1:ix), *) ping_s%bpm(n_bpm)%y(nt)

    read (1, '(a)', iostat = ios) line
    if (ios /= 0) exit
    if (line == '') exit
    if (line(1:1) == '#') exit
  enddo

  if (ios /= 0) exit
  if (n_bpm == size(ping_s%bpm)) exit  ! So blank lines at end do not cause problems.

enddo

close(1)

end subroutine

end module
