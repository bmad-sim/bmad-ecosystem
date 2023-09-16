!+ 
! Program wall_generator
!
! Parses lattice file and outputs wall.out file with points that
! can be linearly interpolated to make wall geometry
! Note: has header for particle 1, as if it were particle track data
!
! Modules Needed:
!   use bmad
!
! Input (command line)
!   [lattice file]
!
! Output
!   wall.dat
!-

program wall_generator

use bmad
use wall3d_mod
use time_tracker_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (wall3d_struct), pointer :: wall3d
type (coord_struct) :: dummy_orb
type (branch_struct), pointer :: branch
type (wall3d_section_struct), pointer :: section
type (coord_struct) :: point
character(100) :: lat_name
character(10) ::  n_anglestring, ds_string
integer :: ix_ele, section_id
integer :: i, n_angles
integer, parameter :: wall_file = 1
real(rp) :: perp(3), r_wall, r0, d_radius, dr_dtheta, ds_offset, s_vertex, s_rel, theta, ds

character(32), parameter :: r_name = 'wall_generator'

logical err

!------------------------------------------

print *, ''
print *, '                         __   ___       ___  __       ___  __   __  '
print *, '|  |  /\  |    |        / _` |__  |\ | |__  |__)  /\   |  /  \ |__) '
print *, '|/\| /~~\ |___ |___ ___ \__> |___ | \| |___ |  \ /~~\  |  \__/ |  \ '
print *, ''





!Get data file name
!call getarg(1, lat_name)


if (command_argument_count() ==0) then
  print *, 'Please provide a lattice file. Format:'
  print *, 'wall_generator <lattice>'
  print *, '  Example: wall_generator lat.bmad'
  print *, 'wall_generator <lattice> <n_angles> <ds>'
  print *, '  Example: wall_generator lat.bmad 8'
  stop
endif 

if (command_argument_count() > 0) call get_command_argument(1, lat_name)
print *, "Creating wall for lattice file: ", lat_name


!Number of angles in a section
if (command_argument_count() > 1) then
  call get_command_argument(2, n_anglestring)
  read(n_anglestring, * ) n_angles 
  if (n_angles < 1) then
    print *, 'Not enough angles: ', n_angles
    stop
  endif
else
  n_angles = 2
endif
print *, 'Using number of angles: ', n_angles

!Number of angles in a section
if (command_argument_count() > 2) then
  call get_command_argument(3, ds_string)
  read(ds_string, * ) ds 
  print *, 'using ds: ', ds
else
  ds = 0
endif

!Parse lattice
call bmad_parser (lat_name, lat)

!-----------------------
!Open output file
open(wall_file, file = "wall.out")




!Make wall border points and write to file
!First set invariant coordinates
point%t = 0.0
point%vec(2) = 0.0
point%vec(4) = 0.0
point%vec(6) = 1.0
point%phase(2) = 0.0
point%vec(3) = 0.0




!Set rest of coordinates

branch => lat%branch(0)

if (.not. associated(branch%wall3d)) then
  print *, 'No branch%wall3d, exiting...'
  stop
endif

call write_header()

wall3d => branch%wall3d(1)
   
do section_id = 1, size(wall3d%section)
    section => wall3d%section(section_id)
    ele => branch%ele(section%ix_ele)
    
    s_rel = section%s - (ele%s -ele%value(L$))
    !Skip sections outside of our element
    !print *, 'ix_ele = ', ele%ix_ele
    !if (s_rel < -1*bmad_com%significant_length .or. s_rel > ele%value(L$) +bmad_com%significant_length) then
    !  !print *, 's_rel = ', s_rel, 'cycling, ', -1*bmad_com%significant_length
    !  cycle
    !end if
    
    do i = 0, n_angles -1
      
      theta = i*2*pi/n_angles
      if (ele%key == sbend$) theta = theta-ele%value(ref_tilt_tot$) ! Correct for ref tilt      
      call calc_wall_radius (section%v, cos(theta), sin(theta), r_wall, dr_dtheta)
    
      !Use a dummy orb to get the global wall position
      ! WARNING: NORMAL IS NOT COMPUTED!!!
      dummy_orb%vec = [r_wall*cos(theta), 0.0_rp, r_wall*sin(theta),  0.0_rp, s_rel, 0.0_rp]
      dummy_orb%ix_ele = ele%ix_ele
      dummy_orb = particle_in_global_frame(dummy_orb, branch, in_time_coordinates = .true., in_body_frame = .false.)
     
      !Write to file: wall border point
      write (wall_file, '(6es18.10, 2i18, es18.10)') dummy_orb%vec, ele%ix_ele, i, section%s
    end do

enddo

!Close data file
close(wall_file)
print *, "Written: wall.out"


if (ds ==0) stop

!-----------------------
!Open output file
open(wall_file, file = "wall_contour.out")

call write_header()

r0 = 1e-3 ! Test radius
point%s = 0.0
point%vec(5) = 0.0
point%vec(6) = 1.0
do
  if (point%s > branch%ele(branch%n_ele_track)%s) exit
  ix_ele = element_at_s (branch, point%s, .true., err)
  ele => branch%ele(ix_ele)
  
  point%vec(5) = point%s - (ele%s -ele%value(L$))
  do i = 0, n_angles -1
    theta = i*2*pi/n_angles
    if (ele%key == sbend$) theta = theta-ele%value(ref_tilt_tot$) ! Correct for ref tilt
    point%vec(1) = r0*cos(theta)
    point%vec(3) = r0*sin(theta)
    d_radius = wall3d_d_radius(point%vec, ele, 1, perp)
    r_wall = r0 - d_radius
  
    !Use a dummy orb to get the global wall position
    dummy_orb%vec = [r_wall*cos(theta), perp(1), r_wall*sin(theta), perp(2), point%vec(5), perp(3)]
    dummy_orb%ix_ele = ele%ix_ele
    dummy_orb = particle_in_global_frame(dummy_orb, branch, in_time_coordinates = .true., in_body_frame = .false.)
    !Write to file: wall border point
    write (wall_file, '(6es18.10, 2i18, es18.10)') dummy_orb%vec, ele%ix_ele, i, point%s
  enddo 
  
  !print *, point%s, point%vec(5), r_wall
  point%s = point%s + ds

enddo

!Close data file
close(wall_file)
print *, "Written: wall_contour.out"

contains

subroutine write_header()
!Write to file: header

! write (wall_file, '(2a)')   '# Wall for:', trim(lat%input_file_name)
write (wall_file, '(9a18)') '#          x', 'normal_x', 'y', 'normal_y', 'z', 'normal_z', 'ix_ele', 'angle_index', 's'
write (wall_file, '(9a18)') '#          m', '1',        'm', '1',        'm', '1',       '1',       '1',  'm'
end subroutine

end program
