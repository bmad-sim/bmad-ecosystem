program beam_file_translate_format

use beam_file_io
use hdf5_interface

implicit none

type (beam_struct) beam
type (beam_init_struct) beam_init
type (lat_struct) lat

logical err_flag
character(16) out_fmt
character(200) beam_file, lat_name, full_name

!

call get_command_argument(1, beam_file)
call get_command_argument(2, lat_name)
call get_command_argument(3, out_fmt)

if (beam_file == '') then
  print '(a)', 'Usage:'
  print '(a)', '   beam_file_translate_format <beam_file> <lat_file> <out_format>'
  print '(a)', ' where:'
  print '(a)', '   <beam_file> is the name of the beamfile.'
  print '(a)', '   <lat_file> (optional) is the name of a lattice file which is possibly needed with hdf5 format input.'
  print '(a)', '   <out_format> (optional) Output format can be one of:'
  print '(a)', '     ascii      - ASCII format. Default if <beam_file> has ".h5" or ".hdf5" suffix.'
  print '(a)', '     hdf5       - HDF5 binary format. Default if <beam_file> does not have ".h5" or ".hdf5" suffix.'
  print '(a)', '     old_ascii  - Old ASCII format. Only to be used for testing.'
  print '(a)', '     old_binary - Old binary format. Only to be used for testing.'
  stop
endif

if (lat_name == '') then
  call read_beam_file(beam_file, beam, beam_init, err_flag)
else
  call bmad_parser(lat_name, lat)
  call read_beam_file(beam_file, beam, beam_init, err_flag, lat%ele(0))
endif  

if (err_flag) stop


if (out_fmt == '') then
  if (index(beam_file, '.h5') == 0 .and. index(beam_file, '.hdf5') == 0) then
    out_fmt = 'hdf5'
  else
    out_fmt = 'ascii'
  endif
endif

call fullfilename(beam_file, full_name)
select case (out_fmt)
case ('ascii')
  call file_suffixer(full_name, full_name, '.dat', .true.)
  call write_beam_file (full_name, beam, .true., ascii$)

case ('hdf5')
  call file_suffixer(full_name, full_name, '.hdf5', .true.)
  call write_beam_file (full_name, beam, .true., hdf5$)

case ('old_ascii')
  call file_suffixer(full_name, full_name, '.dat', .true.)
  call write_beam_file (full_name, beam, .true., old_ascii$)

case ('old_binary')
  call file_suffixer(full_name, full_name, '.bin', .true.)
  call write_beam_file (full_name, beam, .true., binary$)
end select

print '(a)', 'Written: ', trim(full_name)

end program
