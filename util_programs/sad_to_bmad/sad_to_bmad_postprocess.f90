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
use pointer_lattice, lat_ptc => lat, dpe => dp
use ptc_layout_mod

implicit none

type (lat_struct), target :: lat
type (coord_struct), allocatable :: orbit(:)
type (ele_struct), pointer :: ele

real(rp) f_shift, z_old
real(rp), allocatable :: t_shift(:)

integer nn, i, ios, ix

logical fshift_found

character(40), allocatable :: name(:)
character(40) pname, calc_fshift_for
character(100) lat_file
character(200) line

! Read in lattice file 

if (command_argument_count() /= 2) then
  print *, 'Command line syntax:'
  print *, '  sad_to_bmad_postprocess <lattice-file-name> <calc_fshift_for>'
  stop
endif

call get_command_argument (1, lat_file)
call get_command_argument (2, calc_fshift_for)
call bmad_parser (lat_file, lat)

calc_fshift_for = upcase(calc_fshift_for)

nn = 0
do i = 1, lat%n_ele_track
  ele => lat%ele(i)
  if (ele%key == patch$) nn = nn + 1
enddo

allocate (name(nn))
allocate (t_shift(nn))

!-------------
! PTC way
! Note: find_time_patch looks for rfcavities and assumes the patch is the element before the cavity.

if (calc_fshift_for == 'PTC') then
  call lat_to_ptc_layout (lat)
  call in_ptc_units() ! find_time_patch does not work with Bmad units
  call find_time_patch (lat%branch(0)%ptc%m_t_layout, DEFAULT, 1d40)

  nn = 0
  do i = 1, lat%n_ele_track
    ele => lat%ele(i)
    if (ele%key == patch$) then
!      print '(a, i6, 2x, a20, i2, 2es13.4)', 'Patch:', i, ele%name, &
!                              ele%ptc_fibre%patch%time, ele%ptc_fibre%patch%a_t, ele%ptc_fibre%patch%b_t
      if (ele%ptc_fibre%patch%time /= 0) then
        nn = nn + 1
        name(nn) = downcase(ele%name)
        t_shift(nn) = (ele%ptc_fibre%patch%a_t + ele%ptc_fibre%patch%b_t) / c_light
      endif

    elseif (ele%ptc_fibre%patch%time /= 0) then
      print '(a, i6, 2x, a20, i2, 2es13.4)', 'ERROR IN SAD_TO_BMAD_POSTPROCESS: ', &
                              i, ele%name, ele%ptc_fibre%patch%time, ele%ptc_fibre%patch%a_t, ele%ptc_fibre%patch%b_t
    endif
  enddo

!-------------
! Bmad way

else

  call set_on_off (rfcavity$, lat, off$)
  call reallocate_coord (orbit, lat)
  call twiss_and_track (lat, orbit)

  z_old = 0
  
  nn = 0
  do i = 1, lat%n_ele_track
    ele => lat%ele(i)
    if (ele%key /= patch$) cycle
    if (lat%ele(i-1)%key /= rfcavity$ .and. lat%ele(i+1)%key /= rfcavity$) cycle
    nn = nn + 1
    name(nn) = downcase(ele%name)
    t_shift(nn) = (z_old - orbit(i)%vec(5)) / c_light
    z_old = orbit(i)%vec(5)
  enddo

  ! Last patch must time offset to the end of the cavity
  ix = lat%n_ele_track
  t_shift(nn) = t_shift(nn) + (z_old - orbit(ix)%vec(5)) / c_light

endif

!-------------
! Create file with t_offset values

open (1, file = lat_file)
open (2, file = 'temp.temp')

do
  read (1, '(a)', iostat = ios) line
  if (ios /= 0) exit

  if (index(line, '! Will be replaced by sad_to_bmad_postprocess') /= 0) then
    ix = index(line, '=')
    pname = line(3:ix-1)
    do i = 1, nn
      if (name(i) == pname) exit
      if (i == nn) then
        print *, 'ERROR IN SAD_TO_BMAD_POSTPROCESS: CANNOT MATCH PATCH NAME: ', pname
        call err_exit
      endif
    enddo
    write (line, '(a, es16.8)') line(1:ix), t_shift(i)
    print '(i4, 2x, a, 2x, a20, es16.8)', i, 'Time Patch:', pname, t_shift(i)
  endif

  write (2, '(a)') trim(line)    
enddo

close(1)
close(2)

! And now move temp file to lattice file

call system_command ('mv temp.temp ' // trim(lat_file))
print *, 'Sad_to_bmad_postprocess: Modified file: ', trim(lat_file)

end program
