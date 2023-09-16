program synrad_aperture_to_wall

use bmad

implicit none

type (lat_struct) lat
type (coord_struct), allocatable :: orb(:)
type (coord_struct) orb_at_s
type (ele_struct), pointer :: ele
type (ele_struct), pointer :: slave

integer i, n
real(rp) s_position

logical first
real(rp) x1first
real(rp) x2first

character(100) lat_name
character(100) outfile_name
character(24)  :: fmt

!------------------

!Get lattice from command line
if (command_argument_count()==0) then
	print *,"ERROR: Please provide a lattice file"
	stop
endif

call get_command_argument(1, lat_name)
print *,"Using ", lat_name
call file_suffixer (lat_name, outfile_name, '.wall_dat', .true.)

call bmad_parser (lat_name, lat)
call twiss_and_track (lat, orb)

open(1, file = outfile_name)

!Format Statements
fmt  = '(3es18.10)'

!Header notes
   write (1, '(a)') '! Wall file for use with Synrad'
   write (1, '(a)') '!' 
   write (1, '(a)') '! Note: x_minus should be negative.'
   write (1, '(a)') '! Note: First s_position should be 0.' 
   write (1, '(a)') '! Note: Last s_position will be changed to the length of the ring' 
   write (1, '(a)') '!       so its value is ignored.' 
   write (1, '(a)') '!'
   write (1, '(a)') '! s_position(m)      x_minus(m)          x_plus(m) ' 

!---Go through lattice: only use branch 0

first = .true.

do i = 0, lat%branch(0)%n_ele_track
   ele => lat%branch(0)%ele(i)
   
   !Ignore zero (infinite in bmad) apertures
   if (ele%value(x1_limit$) < 1e-9) cycle
   if (ele%value(x2_limit$) < 1e-9) cycle

   if (first)then
      first = .false.
      x1first = -ele%value(x1_limit$)
      x2first =  ele%value(x2_limit$)
   endif
   
   !write to file
   if ( ele%aperture_at == both_ends$ ) then
      write (1, fmt) ele%s - ele%value(L$), -ele%value(x1_limit$), ele%value(x2_limit$)
      write (1, fmt) ele%s,                 -ele%value(x1_limit$), ele%value(x2_limit$)
   else if  ( ele%aperture_at == exit_end$ ) then
      write (1, fmt) ele%s,                 -ele%value(x1_limit$), ele%value(x2_limit$)
   else if  ( ele%aperture_at == entrance_end$ ) then
      write (1, fmt) ele%s - ele%value(L$), -ele%value(x1_limit$), ele%value(x2_limit$)
   endif
enddo

! Add final line with beginning aperture
ele => lat%branch(0)%ele(0)
write (1, fmt) 0.0_rp, x1first, x2first

close(1)

print *,'Created ', outfile_name

end program

