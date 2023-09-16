module dark_current_mod

use bmad_struct
use wall3d_mod
use time_tracker_mod
use quick_plot



type dark_current_param_struct
  character(100) particle_file_name
  logical save_tracks                       ! track points will be saved
  logical save_field                        ! em field at track points will be saved
  logical global_frame                      ! tracks are written in the global frame 
  real(rp) dt_save                          ! time interval to save track point
  logical verbose                           ! print extra information to the screen
  logical plot_on                           ! tracks will be plotted
  integer id                                ! id for debugging parallel

end type dark_current_param_struct

type dark_current_tally_struct
  real(rp) :: charge = 0                 ! Total charge deposited
  real(rp) :: energy = 0                 ! Total energy deposited 
end type




contains


!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! subroutine plot_particle_track (track)
!
!
!
!-

subroutine plot_particle_track(lat, track, scale, in_global_frame)



implicit none

type (lat_struct) lat
type (track_struct) :: track
type (coord_struct) :: orb
real(rp), allocatable :: x(:), y(:)
real(rp) :: color, red, green, blue, sc
integer :: point_id
real(rp), optional :: scale 

logical :: g_frame
logical, optional :: in_global_frame

g_frame = logic_option(.true., in_global_frame)

!
if (present(scale)) then
 sc = scale
else 
 sc = .011_rp
endif

allocate(x(0:track%n_pt) )
allocate(y(0:track%n_pt) )

do point_id = 0, track%n_pt
  if (g_frame) then
    orb = particle_in_global_frame(track%pt(point_id)%orb, lat%branch(0))
    x(point_id) = orb%vec(5)
    y(point_id) = orb%vec(1)
  else
   x(point_id) = track%pt(point_id)%orb%vec(5)
   y(point_id) = track%pt(point_id)%orb%vec(1)
  endif

  !
  ! y(point_id) = hypot(orb%vec(1), orb%vec(3))
end do


! Color by s_origin
!color = sc*(track%pt(0)%orb%s/lat%ele(lat%n_ele_track)%s) 

! Color by initial radius
color = abs(sqrt(x(0)**2 + y(0)**2)/sc)

if (color > 1.0_rp) print *, 'Warning: color > 1: ', color
!red = color
!green = 0.0
!blue =  1.0 - color

!print *, 'color = ', color
!CALL PGSHLS(15, color*360, 0.0, 1.0)
!CALL PGSCR(15, red, green, blue)
!call qp_set_line_attrib ('PLOT', color =red$)
!print *, 'color ', color
!call plcol0(green$)
!call plcol1(color)
!global_com%exit_on_error = .false.
!call qp_set_line_attrib ('PLOT', color = 15)
!call plline (size(x), x, y)
!call qp_draw_data (x, y, symbol_every = 0)
call qp_set_color_basic (qp_continuous_color(color))

call qp_draw_polyline_no_set(x, y)

! Needed for macOS to finish drawing. TODO: debug this. 
call qp_draw_symbol (x(track%n_pt), y(track%n_pt))   

deallocate(x, y)

end subroutine plot_particle_track




!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! subroutine plot_wall (lat)
!
!
!
!-

subroutine plot_wall (lat)

type (lat_struct) :: lat
type (floor_position_struct), allocatable ::  wall_position(:)
real(rp), allocatable :: x(:), y(:)
real(rp) s0, s1, angle

integer n

!

s0 =  0.0_rp
s1 = lat%ele(lat%n_ele_track)%s
!call qp_set_line_attrib ('PLOT', color =black$)
angle = 0.0_rp
call wall_position_from_s_to_s ( s0, s1, angle, lat, wall_position, in_global_frame = .true.)
n = size(wall_position)
allocate(x(0:n-1))
allocate(y(0:n-1))

y = wall_position(:)%r(1)
x = wall_position(:)%r(3)
call qp_draw_polyline_no_set(x, y)

angle = pi
call wall_position_from_s_to_s ( s0, s1, angle, lat, wall_position, in_global_frame = .true.)
y = wall_position(:)%r(1)
x = wall_position(:)%r(3)
call qp_draw_polyline_no_set(x, y)


deallocate(x, y)

end subroutine plot_wall



!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! subroutine write_particle_track (outfile, track, write_field, write_header)
!
!
!
!-

subroutine write_particle_track (outfile, lat, track, global_frame, write_field, write_header)

implicit none

type (track_struct), optional :: track
type (lat_struct), optional, target :: lat
type (branch_struct), pointer :: branch
type (coord_struct) :: orb
integer :: outfile
logical, optional :: write_field, write_header, global_frame
integer :: point_id
logical :: w_track, w_header, w_field, g_frame
real(rp) :: w_mat(3,3)

character(30), parameter :: r_name = 'write_particle_track'
!

if (present(track)) then
  w_track = .true.
else
  w_track = .false.
endif

branch => null()
if (present(lat)) branch => lat%branch(0)

w_header = logic_option (.false., write_header)
w_field  = logic_option (.false., write_field)
g_frame = logic_option (.false., global_frame)

! Do some error checking
if (.not. present(lat) .and. w_track .and. g_frame ) then
  call out_io (s_error$, r_name, 'WRITE TRACKS IN GLOBAL FRAME CALLED WITHOUT A LAT. NOTHING WRITTEN.')
  return
endif 

if (g_frame) then
  if (w_field) then
    if (w_header) then
      write (outfile, '(15a19)') 'x', 'cp_x', 'y', 'cp_y', 'z', 'cp_z', 't', 'charge', 'hit_angle', 'E_x', 'E_y', 'E_z', 'B_x', 'B_y', 'B_z' 
      write (outfile, '(15a19)') 'm', 'eV', 	'm', 'eV',   'm', 'eV',  's', 'C',      'rad',       'V/m', 'V/m', 'V/m', 'T',   'T',   'T'      
    endif
    if (w_track) then 
      !Write lines with E and B fields. Use w_mat to rotate the field vectors into the global frame
      do point_id = 0, track%n_pt
        orb = particle_in_global_frame(track%pt(point_id)%orb, branch, w_mat_out = w_mat)
        write (outfile, '(15es19.10E3)') &
        orb%vec(1:6),&
        orb%t, &
        orb%charge, &
        orb%phase(2), &
        matmul(w_mat, track%pt(point_id)%field%E(1:3)), & 
        matmul(w_mat, track%pt(point_id)%field%B(1:3))
      end do
    endif
  else ! No field    
    if (w_header) then
      write (outfile, '(9a19)') 'x', 'cp_x', 'y', 'cp_y', 'z', 'cp_z', 't', 'charge', 'hit_angle'
      write (outfile, '(9a19)') 'm', 'eV', 	'm', 'eV',   'm', 'eV', 's', 'C',      'rad'   
    endif
    if (w_track) then 
      !Write lines with E and B fields. Use w_mat to rotate the field vectors into the global frame
      do point_id = 0, track%n_pt
        orb = particle_in_global_frame(track%pt(point_id)%orb, branch)
        write (outfile, '(9es19.10E3)') &
        orb%vec(1:6),&
        orb%t, &
        orb%charge, &
        orb%phase(2)
      end do
    endif
  endif

else
!s-coordinates. This is very similar to the above loops.
  if (w_field) then
    if (w_header) then
      write (outfile, '(15a19)') 'x', 'cp_x', 'y', 'cp_y', 's', 'cp_s', 't', 'charge', 'hit_angle', 'E_x', 'E_y', 'E_s', 'B_x', 'B_y', 'B_s' 
      write (outfile, '(15a19)') 'm', 'eV',   'm', 'eV',   'm', 'eV',   's', 'C',      'rad',       'V/m', 'V/m', 'V/m', 'T', 'T', 'T'      
    endif
    if (w_track) then 
      !Write raw orb data with E and B fields. 
      do point_id = 0, track%n_pt
        orb = track%pt(point_id)%orb
        write (outfile, '(15es19.10E3)') &
        orb%vec(1:4),&
        orb%s, &
        orb%vec(6), &
        orb%t, &
        orb%charge, &
        orb%phase(2), & 
        track%pt(point_id)%field%E(1:3), & 
        track%pt(point_id)%field%B(1:3)
      end do
    endif
  else ! No field    
    if (w_header) then
      write (outfile, '(9a19)') 'x', 'cp_x', 'y', 'cp_y', 's', 'cp_s',  't', 'charge', 'hit_angle' 
      write (outfile, '(9a19)') 'm', 'eV',   'm', 'eV',   'm', 'eV',    's', 'C',      'rad'  
    endif
    if (w_track) then 
      !Write raw orb data
      do point_id = 0, track%n_pt
        orb = track%pt(point_id)%orb
        write (outfile, '(9es19.10E3)') &
        orb%vec(1:4),&
        orb%s, &
        orb%vec(6), &
        orb%t, &
        orb%charge, &        
        orb%phase(2)
      end do
    endif
  endif
endif
   

end subroutine write_particle_track

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine   wall_position_from_s_to_s ( s0, s1, angle, lat, wall_position, in_global_frame)


subroutine wall_position_from_s_to_s ( s0, s1, angle, lat, wall_position, in_global_frame)

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (floor_position_struct), allocatable ::  wall_position(:)
type (ele_struct), pointer :: ele
type (wall3d_struct), pointer :: wall3d
type (coord_struct) :: dummy_orb
type (wall3d_section_struct), pointer :: section
type(floor_position_struct), allocatable :: array(:), array2(:)

real(rp) :: s0, s1, angle
real(rp) :: perp(3), r_wall, dr_dangle, ds_offset, s_vertex, s_rel

integer :: ix0, ix1, ix_ele, section_id, section_counter
integer :: n

logical :: err
logical, optional :: in_global_frame

character(30), parameter :: r_name = 'wall_position_from_s_to_s'

!

! work array 
allocate(array(1024))
section_counter = 0

ix0 = element_at_s(lat, s0, .false.)
ix1 = element_at_s(lat, s1, .true.)

branch => lat%branch(0)
wall3d => pointer_to_wall3d (branch%ele(1), 1, ds_offset)   

if (.not. associated(wall3d)) return

do section_id = 1, size(wall3d%section)
  section => wall3d%section(section_id)    

  !if (section%ix_ele < ix0 .or. section%ix_ele > ix1) cycle

  ele => branch%ele(section%ix_ele)
  s_rel = section%s - (ele%s - ele%value(l$))

  ! Main worker routine
  call calc_wall_radius (section%v, cos(angle), sin(angle), r_wall, dr_dangle)
  ! Use a dummy orb to get the global wall position
  dummy_orb%vec = [section%r0(1) + r_wall*cos(angle), 0.0_rp, section%r0(2) + r_wall*sin(angle), 0.0_rp, s_rel, 0.0_rp]
  dummy_orb%ix_ele = ele%ix_ele
  if (logic_option(.true., in_global_frame)) then
    dummy_orb = particle_in_global_frame(dummy_orb, branch, in_time_coordinates = .true., in_body_frame = .false.)
  endif

  ! Array control
  section_counter = section_counter +1
  ! Reallocate temporary structure if needed
  if (section_counter > size(array)) then
    n = size(array)
    allocate( array2(n))
    array2 = array
    deallocate(array)
    allocate(array(2*n))
    array(1:n) = array2
    deallocate(array2)
  end if
 
 array(section_counter)%r(1) = dummy_orb%vec(1) 
 array(section_counter)%r(2) = dummy_orb%vec(3) 
 array(section_counter)%r(3) = dummy_orb%vec(5) 
 array(section_counter)%theta = 0 
 array(section_counter)%phi   = 0 
 array(section_counter)%psi   = 0 
 
 !Write to file: wall border point
 !       write (wall_file, '(6es18.10, 2i18)') dummy_orb%vec, ele%ix_ele, i
 !     end do

end do

! Output array

if (allocated(wall_position)) then
  deallocate(wall_position)
endif

allocate(wall_position(0:section_counter-1) )
wall_position(0:section_counter-1) = array(1:section_counter)

deallocate(array)

end subroutine wall_position_from_s_to_s




!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine import_time_distribution(dist_file, mc2, orb)
!
! Simple routine to import time-based distribution particles from a file. 
!
! Modules Needed:
!   use bmad
!
! Input:
!   dist_file   -- character(*): File name for time distibution.
!                  The time format is: 
!                  nparticles(integer)
!                  x/m   m*c^2*\beta_x*\gamma_x/eV  y/m  m*c^2*\beta_y*\gamma_y z/m m*c^2*\beta_z*\gamma_z/eV time/s charge/C
!                  . 
!                  .
!
! Output:
!  orb(nparticles)   -- coord_struct: array containing all of these particles 
!                                     
!-

subroutine import_time_distribution(dist_file, particles)

use bmad_routine_interface

implicit none

character(*), intent(in) :: dist_file
type (coord_struct),  allocatable, intent(out)    :: particles(:)

real(rp) timevec(8)
integer nparticles, i
logical :: exist
integer :: dist = 999  
character(30), parameter :: r_name = 'import_time_distribution'


inquire(file=dist_file, exist=exist)
if (.not. exist) then
  write (*, *) 'ERROR: particle file does not exist: ', dist_file
  stop 
endif

dist = lunget()
open(dist, file = dist_file, action = 'read', status = 'old')

read(dist, '(i8)') nparticles

if (.not. allocated (particles)) then
  allocate(particles(1:nparticles))  
   !write (outfile, '(a)' )  "Start coordinates"
   !write (outfile, '(6es18.10)') particle_beg%vec
   !write (outfile, '(a)' )  "End coordinates"
   !write (outfile, '(6es18.10)') particle_end%vec
  
endif

do i = 1, nparticles
   read(dist, *) timevec
   !Bmad-T variables use c*px, c*py, c*ps in eV
   particles(i)%vec = [ timevec(1),  timevec(2),  timevec(3),   timevec(4),  timevec(5), timevec(6) ]
   particles(i)%t = timevec(7)
   particles(i)%s = timevec(5)
   particles(i)%charge = timevec(8)
end do

close(dist)

end subroutine import_time_distribution




!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine dark_current_tracker_parse_command_args (dc_param, error, cmd_words)
!
! Subroutine to parse the command line arguments.
!
! adapted from tao_parse_command_args
!
! Input:
!   dc_param     -- dark_current_param_struct 
!   cmd_words(:) -- Character(*), optional: If present then this is used
!                    in place of the command line.
! Output:
!   error -- Logical: Set True if there is an error. False otherwise.
!-

subroutine dark_current_tracker_parse_command_args (dc_param, error, cmd_words)

implicit none

type(dark_current_param_struct) :: dc_param
character(*), optional :: cmd_words(:)
character(80) arg0, base, switch
character(24) :: r_name = 'dark_current_tracker_parse_command_args'

integer n_arg, i_arg, ix
logical error

! Get command line input

error = .false.

if (present(cmd_words)) then
  n_arg = size(cmd_words)
  if (cmd_words(1) == '') return
else
  n_arg = command_argument_count()
  if (n_arg == 0) return
endif

! loop over all arguments

i_arg = 0

do 

  if (i_arg == n_arg) exit
  call get_next_arg (arg0)

  call match_word (arg0,            ['-help                    ', 'help                     ', &
        '-particle_file_name      '], &
      !  '-noplot                  ', '-lat                     ', '-log_startup             ', '-beam                    ', &
      !  '-var                     ', '-data                    ', '-building_wall           ', '-plot                    ', &
      !  '-startup                 ', 'help                     ', '-help                    ', '?                        ', &
      !  '-geometry                ', '-rf_on                   ', '-debug                   ', '-disable_smooth_line_calc'], &
              ix, .true., matched_name=switch)

  select case (switch)

  case ('-particle_file_name')
    call get_next_arg (dc_param%particle_file_name)
  
  case ('help', '-help')
     print *, 'HELP MESSAGE'
     stop
  case default
    call out_io (s_error$, r_name, 'BAD COMMAND LINE ARGUMENT: ' // arg0)
    error = .true.
    return
  end select

enddo

!-----------------------------
contains

subroutine get_next_arg(arg)

character(*) arg

!

if (i_arg == n_arg) then
  call out_io (s_error$, r_name, 'MISSING COMMAND LINE ARGUMENT FOR: ' // arg0)
  error = .true.
  return
endif

i_arg = i_arg + 1

if (present(cmd_words)) then
  arg = cmd_words(i_arg)
else
  call get_command_argument(i_arg, arg)
endif

end subroutine get_next_arg

end subroutine dark_current_tracker_parse_command_args

end module
