!+
! Program tao_ping_program
!-

program tao_ping_program

use tao_ping_utils

implicit none

! Tell Tao to bypass the standard init sequence.

tao_com%init_beam           = .false.
tao_com%init_connected_uni  = .false.
tao_com%init_var            = .false.
tao_com%init_read_lat_info  = .false.
tao_com%init_data           = .false.
tao_com%parse_cmd_args      = .false.

! Do some init

call ping_read_parameters(ping_s)
call ping_read_data (ping_s)
print *, 'Number of bpms: ', size(ping_s%bpm)
print *, 'Number of turns:', size(ping_s%bpm(1)%x)

! And start Tao.

call tao_cl()

end program
