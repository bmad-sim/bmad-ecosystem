!+
! Program bmad_to_csrtrack
!
! Program to convert a Bmad lattice and bunch initialization to 
! CSRtrack format. Optionally only a portion of the lattice can be translated.
!
! Note: Markers before and after each element will be generated automatically.
! If the element name is unique then the marker names will be:
!   <ele_name>_m1
!   <ele_name>_m2
! Otherwise with the <n>th occurance of an element with a given name the
! marker names will be:
!   <ele_name>_m1_r<n>
!   <ele_name>_m2_r<n>
!
! Input:
!   A file named bmad_to_csrtrack.in containing the following:
!
!   1) A general_params namelist section which looks like:
!         &general_params
!           bmad_lattice      = <bmad lattice file name>
!           particle_out_file = <Particle output file name>
!                                 Default: "particles.fmt1"
!           ele_start         = <Starting Bmad element name>
!           ele_end           = <Ending Bmad element name>
!         /
!
!     Note: The lattace created for CSRtrack is from the beginning of
!       ele_start to the end of ele_end.
!
!   2) A beam_params namelist initialization section:
!         &beam_params
!           beam_init%a_norm_emitt = <Normalized emittance>
!           ... etc. ...
!         /
!
!     Note: beam_init is a variable of type beam_init_struct.
!     Use "getf beam_init_struct" to see the subelements of this structure.
!     beam_init%n_bunch = 1 by default and should not be changed. 
!
!   4) The CSRtrack parameters (See the CSRtrack manual).
!      There should be no lattice or particles sections. Instead there should
!      the line:
!         INSERT_LATTICE_AND_BUNCH_PARAMS_HERE
!      and the lattice and particle sections will be inserted at this place.
!
!      Note: Every place the string "$end_marker" is found the name of the last
!      marker in the lattice will be substituted.
!
! Output:
!  1) A "csrtrk.in" file that can be used as the input for CSRtrack.
!  2) A file of the initial particle distribution.
!-

program bmad_to_csrtrack

use bmad
use beam_mod

implicit none

type (lat_struct), target :: lat
type (beam_init_struct) beam_init
type (beam_struct), target :: beam
type (ele_struct), pointer :: ele
type (bunch_struct), pointer :: bunch
type (ele_pointer_struct), allocatable :: eles(:)

integer i, j, n, ix, ix_start, ix_end, n_arg, ios, n_list, n_loc
integer, allocatable :: repeat(:)

character(100) bmad_lattice, input_file, particle_out_file
character(100) line, fmt, csrtrk_file
character(40) ele_start, ele_end, ele_name, csr_type
character(20) prop_name(5), prop_fmt(5)
character(40), allocatable :: ele_name_list(:)
character(44) marker1, marker2

real(rp) prop_value(5), p0, c, r(6), s_old, ds
real(rp) delta_psi1, delta_psi2, deg, hk, vk

namelist / general_params / bmad_lattice, &
                particle_out_file, ele_start, ele_end
namelist / beam_params / beam_init

! Find input file name.

  n_arg = command_argument_count()
  input_file = "bmad_to_csrtrack.in"
  csrtrk_file = "csrtrk.in"

  if (n_arg == 1) call get_command_argument(1, input_file)

  if (n_arg > 1) then
    print *, "Usage: bmad_to_csrtrack {<input_file_name>}"
    print *, "       Default: <input_file_name> = 'bmad_to_csrtrack.in'"
    stop
  endif

! Read the input file.

print *, "Opening: ", trim(input_file)
open (1, file = input_file, status = "old")

ele_start = ""
ele_end = ""
particle_out_file = "particles.fmt1"

beam_init%n_bunch = 1
print *, "Now to read the general_params namelist..."
read (1, nml = general_params)
print *, "Now to read the beam_params namelist..."
read (1, nml = beam_params)

! Read the Bmad lattice and bunch init info.

call bmad_parser (bmad_lattice, lat)
call twiss_propagate_all (lat)

if (ele_start /= "") then
  call lat_ele_locator (ele_start, lat, eles, n_loc)
  if (n_loc == 0) then
    print *, 'Error: Unable to locate element: ', ele_start
    stop
  endif
else
  allocate (eles(1))
  eles(1)%ele => lat%ele(1)
endif
ix_start = eles(1)%ele%ix_ele

ix_end = lat%n_ele_track
if (ele_end /= "") then
  call lat_ele_locator (ele_end, lat, eles, n_loc)
  if (n_loc == 0) then
    print *, 'Error: Unable to locate element: ', ele_end
    stop
  endif
  ix_end = eles(1)%ele%ix_ele
endif

! Create the CSRtrack input file

open (2, file = csrtrk_file)

! Transfer from the input file to the output file until 
! the "INSERT_LATTICE_AND_BUNCH_PARAMS_HERE" line is found

do 
  read (1, '(a)', iostat = ios) line
  if (index(line, "INSERT_LATTICE_AND_BUNCH_PARAMS_HERE") /= 0) exit
  if (ios /= 0) then
    print *, 'Error: "INSERT_LATTICE_AND_BUNCH_PARAMS_HERE" line not found."'
    stop
  endif
  write (2, '(a)') trim(line)
enddo

! Now insert the lattice info.
! First mark which elements have unique names

n_list = 0
n = ix_end - ix_start + 1
allocate(repeat(n), ele_name_list(n))

ele_loop: do i = ix_start, ix_end
  ele_name = lat%ele(i)%name
  do j = 1, n_list
    if (ele_name == ele_name_list(j)) then
      repeat(j) = repeat(j) + 1
      cycle ele_loop
    endif
  enddo
  n_list = n_list + 1
  ele_name_list(n_list) = ele_name
  repeat(n_list) = 1
enddo ele_loop

do j = 1, n_list
  if (repeat(j) == 1) repeat(j) = 0
  if (repeat(j) > 1) repeat(j) = 1
enddo

! Write lattice info

write (2, *)
write (2, '(a)') "!------------------------------------------------"
write (2, '(a)') "! Lattice:"
write (2, *)
write (2, '(a)') "lattice {"

deg = 180 / pi

s_old = lat%ele(ix_start-1)%s
do i = ix_start, ix_end

  ele => lat%ele(i)
  ele_name = ele%name

  marker1 = trim(ele_name) // "_m1"
  marker2 = trim(ele_name) // "_m2"
  call downcase_string(marker1)
  call downcase_string(marker2)

  do j = 1, n_list
    if (ele_name == ele_name_list(j)) then
      if (repeat(j) > 0) then
        write(marker1, '(2a, i0)') trim(marker1), "_r", repeat(j)
        write(marker2, '(2a, i0)') trim(marker2), "_r", repeat(j)
        repeat(j) = repeat(j) + 1
        exit
      endif
    endif
  enddo

  delta_psi1 = 0
  delta_psi2 = 0
  prop_name = (/ "", "", "", "", "" /)

  select case (ele%key)
  case (drift$, marker$)
    cycle

  case (sbend$)
    csr_type = "dipole"
    delta_psi1 = -ele%value(e1$) * deg
    delta_psi2 =  ele%value(e2$) * deg
    prop_name(1) = "r"
    prop_value(1) = ele%value(rho$)
    prop_fmt(1) = "f12.4"

  case (quadrupole$)
    csr_type = "quadrupole"
    prop_name(1:4) = (/ "strength         ", "alpha            ", &
                        "horizontal_offset", "vertical_offset  " /)
    prop_value(1:4) = (/ ele%value(k1$), ele%value(tilt$)*deg, &
                          ele%value(x_offset$), ele%value(y_offset$) /)
    prop_fmt(1:4) = (/ "f0.6", "f0.4", "f0.6", "f0.6" /)

  case (sextupole$)
    csr_type = "multipole"
    prop_name(1:5) = (/ "strength         ", "alpha            ", &
                        "horizontal_offset", "vertical_offset  ", &
                        "poles            " /)
    prop_value(1:5) = (/ ele%value(k2$), ele%value(tilt$)*deg, &
                          ele%value(x_offset$), ele%value(y_offset$), 6.0_rp /)
    prop_fmt(1:5) = (/ "f0.6", "f0.4", "f0.6", "f0.6", "i3  " /)

  case (octupole$)
    csr_type = "multipole"
    prop_name(1:5) = (/ "strength         ", "alpha            ", &
                        "horizontal_offset", "vertical_offset  ", &
                        "poles            " /)
    prop_value(1:5) = (/ ele%value(k3$), ele%value(tilt$)*deg, &
                          ele%value(x_offset$), ele%value(y_offset$), 8.0_rp /)
    prop_fmt(1:5) = (/ "f0.6", "f0.4", "f0.6", "f0.6", "i3  " /)

  case (hkicker$)
    if (ele%value(kick$) == 0) cycle

    csr_type = "multipole"
    prop_name(1:5) = (/ "strength         ", "alpha            ", &
                        "horizontal_offset", "vertical_offset  ", &
                        "poles            " /)
    prop_value(1:5) = (/ ele%value(kick$), ele%value(tilt$)*deg, &
                          ele%value(x_offset$), ele%value(y_offset$), 2.0_rp /)
    prop_fmt(1:5) = (/ "f0.6", "f0.4", "f0.6", "f0.6", "i3  " /)

  case (vkicker$)
    if (ele%value(kick$) == 0) cycle

    csr_type = "multipole"
    prop_name(1:5) = (/ "strength         ", "alpha            ", &
                        "horizontal_offset", "vertical_offset  ", &
                        "poles            " /)
    prop_value(1:5) = (/ ele%value(kick$), ele%value(tilt$)*deg + 90, &
                          ele%value(x_offset$), ele%value(y_offset$), 2.0_rp /)
    prop_fmt(1:5) = (/ "f0.6", "f0.4", "f0.6", "f0.6", "i3  " /)


  case (kicker$)
    hk = ele%value(hkick$) 
    vk = ele%value(vkick$)
    if (hk == 0 .and. vk == 0) cycle

    csr_type = "multipole"
    prop_name(1:5) = (/ "strength         ", "alpha            ", &
                        "horizontal_offset", "vertical_offset  ", &
                        "poles            " /)
    prop_value(1:5) = (/ sqrt(hk**2 + vk**2), ele%value(tilt$)*deg + atan2(vk, hk)*180/pi, &
                          ele%value(x_offset$), ele%value(y_offset$), 2.0_rp /)
    prop_fmt(1:5) = (/ "f0.6", "f0.4", "f0.6", "f0.6", "i3  " /)

  case default
    print *, 'Error: Bmad lattice element: ', trim(ele%name)
    print *, '       Which is a: ', key_name(ele%key)
    print *, '       Does not have a corresponding CSRtrack element'
    stop
  end select

  ! write element info

  write (2, '(3x, 3a)') trim(csr_type), ' {   ! ', trim(ele%name)

  ds = lat%ele(i-1)%s - s_old
  if (ds == 0) ds = 1e-6  ! So CSRTrack does not complain.

  if (i == ix_start) then
    write (2, '(6x, 8a)') &
              'position {rho = ', trim(to_string(ds * cos(delta_psi1), 'f0.3')), &
              ', psi = ', trim(to_string(delta_psi1, 'f0.4')), &
              ', marker = ', trim(marker1), '}'
  else
    write (2, '(6x, 8a)') &
              'position {delta_s = ', trim(to_string(ds, 'f0.3')), &
              ', delta_psi = ', trim(to_string(delta_psi1, 'f0.4')), &
              ', marker = ', trim(marker1), '}'
  endif

  write (2, '(6x, a)') 'properties {'
  do j = 1, size(prop_name)
    if (prop_name(j) == "") exit
    write (2, '(9x, 3a)') trim(prop_name(j)), ' = ', &
                                  trim(to_string(prop_value(j), prop_fmt(j)))
  enddo
  write (2, '(6x, a)') '}'

  write (2, '(6x, 8a)') &
              'position {delta_s = ', trim(to_string(ele%value(l$), 'f0.3')), &
              ', delta_psi = ', trim(to_string(delta_psi2, 'f0.4')), &
              ', marker = ', trim(marker2), '}'

  write (2, '(3x, a)') '}' 

  s_old = ele%s

enddo

write (2, '(a)') "}"

! Now write the particle distribution info.

call convert_total_energy_to (lat%ele(ix_start-1)%value(E_tot$), &
                                              electron$, pc = p0)

write (2, *)
write (2, '(a)') "!---------------------------------------------"
write (2, '(a)') "! Particle distribution"
write (2, *)
write (2, '(a)')  "particles {"
write (2, '(2a)') "   reference_momentum  = reference_particle"  !, trim(to_string(p0, 'es15.6'))
write (2, '(a)')  "   reference_point_x   = 0.0"
write (2, '(a)')  "   reference_point_y   = 0.0"
write (2, '(a)')  "   reference_point_phi = 0.0"
write (2, '(a)')  "   format = fmt1"
write (2, '(3a)') "   array = #file{name = ", trim(particle_out_file), "}"
write (2, '(a)')  "}"


! Transfer the rest of file

do 
  read (1, '(a)', iostat = ios) line
  ix = index(line, '$end_marker')
  if (ix > 0) line = line(1:ix-1) // trim(marker2) // trim(line(ix+11:))
  if (ios /= 0) exit
  write (2, '(a)') trim(line)
enddo

close(2)
print *, 'Created: ', trim(csrtrk_file)

!----------------------------------------------
! Write the particle distribution file. 
! Conversion:
!   CSRtrack     Bmad
!     x =         s
!     y =         x
!     z =         y

call init_beam_distribution (lat%ele(ix_start-1), lat%param, beam_init, beam)
c = beam%bunch(1)%particle(1)%charge

open (2, file = particle_out_file)

fmt = '(7es15.6)'
write (2, '(7(14x, a))') "0", "0", "0", "0", "0", "0", "0"

r = beam_init%center
write (2, fmt) r(5), r(1), r(3), (1 + r(6))*p0, r(2)*p0, r(4)*p0, c

do i = 1, size(beam%bunch(1)%particle)
  r = beam%bunch(1)%particle(i)%vec - beam_init%center
  write (2, fmt) r(5), r(1), r(3), r(6)*p0, r(2)*p0, r(4)*p0, c
enddo

close(2)
print *, 'Created: ', trim(particle_out_file)

stop

!----------------------------------------------
9000 continue
print *, 'Error: "begin_csrtrack_params" line not found.'
stop

9010 continue
print *, 'Error: "io_path" line not found.'
stop

9020 continue
print *, 'Error: io_path ending "}" not found.'
stop

!-----------------------------------------------------------------------------
contains

function to_string(r_value, frmt)  result (str)

real(rp) r_value
character(*) frmt
character(20) str, fmt
integer ix

!

fmt = "(" // trim(frmt) // ")"

if (frmt(1:1) == "i") then
  write (str, fmt) nint(r_value)  
else
  write (str, fmt) r_value
endif

call string_trim(str, str, ix)

end function

end program
