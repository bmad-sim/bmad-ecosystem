program beam_file_translate_format

use beam_file_io
use hdf5_interface

implicit none

type (beam_struct) beam
type (beam_init_struct) beam_init
logical err_flag
character(16) out_fmt
character(200) file_name, full_name

!

if (command_argument_count() < 2) then
  print *, 'Usage:'
  print *, '   beam_file_translate_format <beam_file> <out_format>'
  print *, ' where <out_format> is one of:'
  print *, '   ascii'
  print *, '   hdf5'
  print *, '   binary  (Old format! Do not use unless you know what you are doing!)'
  print *, 'NOTE: Time based coordinates for ASCII files not handled by this program.'
  stop
endif

call get_command_argument(1, file_name)
call fullfilename(file_name, full_name)

call get_command_argument(2, out_fmt)

call read_beam_file(full_name, beam, beam_init, err_flag)
if (err_flag) stop

select case (out_fmt)
case ('ascii')
  call file_suffixer(full_name, full_name, '.dat', .true.)
  call write_beam_file (full_name, beam, .true., ascii$)

case ('hdf5')
  call file_suffixer(full_name, full_name, '.hdf5', .true.)
  call write_beam_file (full_name, beam, .true., hdf5$)

case ('binary')
  call file_suffixer(full_name, full_name, '.bin', .true.)
  call write_beam_file (full_name, beam, .true., binary$)
end select

print *, 'Written: ', trim(full_name)

end program
