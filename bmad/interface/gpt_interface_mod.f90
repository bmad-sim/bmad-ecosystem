module gpt_interface_mod

use bmad_interface

implicit none

type gpt_lat_param_struct
  integer :: fieldmap_dimension = 3                   ! Dimensions for field map. 1 or 3
  logical :: only_write_autophase_parameters = .false. ! Option to only write phasing info
  character(100) :: gpt_filename = ''                  ! Blank => Append '.gpt' to Bmad lattice file name.
  character(100) :: header_file_name = ''              ! Header file to include in gpt file.
  character(40) :: tracking_end_element = ''           ! Bmad lattice element name or index.
end type

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine gpt_to_particle_bunch (gpt_file, ele, bunch, err_flag)
!
! Routine to initialize a bunch of particles from a GPT screen file.
! 
!
! Input:
!   gpt_file  -- character(*): Name of GPT data file.
!   ele       -- ele_struct: Lattice element whose downstream end coincident with the GPT screen.
!
! Output:
!   bunch     -- bunch_struct: Particle bunch
!   err_flag       -- logical: Set True if there is an error. False otherwise.
!-

subroutine gpt_to_particle_bunch (gpt_file, ele, bunch, err_flag)

type (ele_struct) ele
type (bunch_struct), target :: bunch
type (coord_struct), pointer :: p

real(rp) time, gamma, vec(6)
real(rp) :: zero_vec(6) = 0
integer i, ic, iu, ix, ios, i_col_max, col_id(20)
logical err_flag

character(*) gpt_file
character(350) line
character(*), parameter :: r_name = 'gpt_to_particle_bunch'

! Open file

err_flag = .true.
iu = lunget()
open (iu, file = gpt_file, status = 'old', iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN GPT FILE: ' // gpt_file)
  return
endif

! Read header to "time" line.

do 
  read (iu, *, iostat = ios) line 
  if (ios /= 0) goto 9000  ! And exit
  if (line(1:4) == 'time') then
    read (line(5:), *, iostat = ios) time
    if (ios /= 0) goto 9000  ! And exit
    exit
  endif
enddo

! Read columns

col_id = -1   
read (iu, *, iostat = ios) line 
if (ios /= 0) goto 9000  ! And exit

ix = 0
do i = 1, 100
  call string_trim(line(ix+1:), line, ix)
  if (ix == 0) exit
  select case (line(1:ix))
  case ('x');   i_col_max = i;  col_id(i) = 1  ! x
  case ('Bx');  i_col_max = i;  col_id(i) = 2  ! beta_x
  case ('y');   i_col_max = i;  col_id(i) = 3  ! y
  case ('By');  i_col_max = i;  col_id(i) = 4  ! beta_y
  case ('z');   i_col_max = i;  col_id(i) = 5  ! z
  case ('Bz');  i_col_max = i;  col_id(i) = 6  ! beta_z
  case ('G');   i_col_max = i;  col_id(i) = 7  ! Gamma factor
  case ('q');   i_col_max = i;  col_id(i) = 8  ! Charge
  end select
  
enddo

! Read data

call reallocate_coord(bunch%particle, 1000)
do i = 1, 100000
  read (iu, *, iostat = ios) line 
  if (ios /= 0) goto 9000  ! And exit
  if (line(1:4) == 'time') then
    call out_io (s_error$, r_name, 'CONFUSED! MULTIPLE "TIME" DATA BLOCKS IN GPT FILE: ' // gpt_file)
    return
  endif

  if (i > size(bunch%particle)) call reallocate_coord(bunch%particle, 2*i)
  p => bunch%particle(i)

  ix = 0
  do ic = 1, i_col_max
    call string_trim(line(ix+1:), line, ix)
    select case (col_id(ic))
    case (1, 2, 3, 4, 5, 6);    read (line(1:ix), *, err = 9000) vec(col_id(ic))
    case (7);                   read (line(1:ix), *, err = 9000) gamma
    case (8);                   read (line(1:ix), *, err = 9000) p%charge
    end select
  enddo

  call init_coord (p, zero_vec, ele, downstream_end$)

  p%beta = sqrt(vec(2)**2 + vec(4)**2 + vec(6)**2)
  p%vec(1) = vec(1) - vec(5) * vec(2) / vec(6)
  p%vec(2) = vec(2) / p%beta
  p%vec(3) = vec(3) - vec(5) * vec(4) / vec(6)
  p%vec(4) = vec(4) / p%beta
  p%vec(5) = p%beta * c_light * (ele%ref_time - (time - vec(5) / (c_light * vec(6))))
  p%vec(6) = (p%beta * gamma * mass_of(electron$) - ele%value(p0c$)) / ele%value(p0c$)
  p%t = time

enddo

! Cleanup

call reallocate_coord(bunch%particle, i-1)
err_flag = .false.
close (iu)
return

! Error reading file

9000 continue
call out_io (s_error$, r_name, 'ERROR READING GPT FILE: ' // gpt_file)
close (iu)

end subroutine gpt_to_particle_bunch

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine write_gpt_lattice_file (lat, gpt_lat_param, err)
!
! Subroutine to write a gpt lattice file using the information in a Bmad lattice.
!
! Input:
!   gpt_lat_param -- gpt_lat_param_struct: Parameters for constructing the lattice
!   lat           -- lat_struct: Holds the lattice information.
!
! Output:
!   err           -- Logical: Set True if there is an error
!-

subroutine write_gpt_lattice_file (lat, gpt_lat_param, err)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, end_ele
type (str_index_struct) :: fieldgrid_names
type (gpt_lat_param_struct), target :: gpt_lat_param
type (gpt_lat_param_struct), pointer :: param
type (ele_pointer_struct), allocatable :: eles(:)
type (branch_struct), pointer :: branch

integer :: n, ie, ix_start, ix_end, iu
integer :: gpt_file_unit, ios, n_loc

character(4000) :: line
character(40) :: name
character(*), parameter :: r_name = 'write_gpt_lat_file'

logical :: err, has_bend

! Open file

err = .true.
param => gpt_lat_param

gpt_file_unit = lunget()  
open (gpt_file_unit, file = param%gpt_filename, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(param%gpt_filename))
  return
endif

! Add header file

if (param%header_file_name /= '') then
  iu = lunget()
  open (iu, file = param%header_file_name, status = 'old', iostat = ios)
  if (ios /= 0) then
    call out_io (s_error$, r_name, 'CANNOT OPEN HEADER FILE: ' // trim(param%header_file_name))
    return
  endif
  do 
    read (iu, '(a)', iostat = ios) line
    if (ios /= 0) exit
    write (gpt_file_unit, '(a)') trim(line)
  enddo
endif

!

lat%ele%ix_pointer = 0
ix_start = 1

if (param%tracking_end_element == '') then
  end_ele => lat%ele(lat%n_ele_track)
else
  call lat_ele_locator (param%tracking_end_element, lat, eles, n_loc)
  if (n_loc == 0) then
    call out_io (s_error$, r_name, 'CANNOT FIND ELEMENT AT END OF GPT TRACKING: ' // param%tracking_end_element)
    return
  endif

  if (n_loc > 1) then
    call out_io (s_error$, r_name, 'MULTIPLE ELEMENTS MATCHING END ELEMENT NAME: ' // param%tracking_end_element)
    return
  endif

  end_ele => eles(1)%ele
  if (end_ele%lord_status == super_lord$) then
    end_ele => pointer_to_slave(end_ele, end_ele%n_slave)
  endif
endif

branch => lat%branch(end_ele%ix_branch)
branch%ele%ix_pointer = 0

! Initialize fieldgrid filename list

n = 0
do ie = ix_start, end_ele%ix_ele
  ele => branch%ele(ie)
  if (.not. skip_ele(ele)) n = n + 1
  if (ele%slave_status == super_slave$) n = n + ele%n_lord
enddo

fieldgrid_names%n_max = 0

! Write element info to GPT file.

write(gpt_file_unit, '(a)')  '# GPT lattice'
write(gpt_file_unit, '(2a)') '# ', trim(lat%lattice)

has_bend = .false.
do ie = ix_start, end_ele%ix_ele
  ele => lat%ele(ie)
  call write_to_file(ele)
  if (ele%key == sbend$) has_bend = .true.
enddo

!

if (has_bend) then
  call out_io (s_error$, r_name, 'SCREEN COMMAND MUST BE MODIFIED DUE TO BEND IN LATTICE!.')
endif

write (gpt_file_unit, '(a, f14.8, a)') 'screen("wcs", "I", ', ele%s, ');'

write(gpt_file_unit, '(a)') '# END GPT lattice'

err = .false.

!----------------------------------------------
! Function to skip non-physical elements

contains

function skip_ele(ele) result (skip)
type (ele_struct) ele
logical skip
! Skip these elements:
if (ele%slave_status == super_slave$ .or. &
      ele%slave_status == multipass_slave$ .or. &
      ele%key == drift$ .or. &
      ele%key == girder$ .or. &
      ele%key == overlay$ .or. &
      ele%key == group$ .or. &
      ele%ix_pointer == 1) then
  skip = .true.
else
  skip = .false.
endif
end function skip_ele

!----------------------------------------------
! contains

recursive subroutine write_to_file (ele)

type (ele_struct) ele
integer i

!

do i = 1, ele%n_lord
  call write_to_file(pointer_to_lord(ele, i))
enddo

if (skip_ele(ele)) return
ele%ix_pointer = 1

! Clean up symbols in element name

name = ele%name
call str_substitute (name, "#", "_part_")
call str_substitute (name, "\", "_and_")      ! "
call str_substitute (name, ".", "_")  

if (param%only_write_autophase_parameters) then
  if (ele%key == lcavity$ .or. ele%key == rfcavity$ .or. ele%key == e_gun$)  then
    call write_gpt_ele(gpt_file_unit, ele, name, fieldgrid_names, dimensions = param%fieldmap_dimension, only_phasing = .true.)
  endif  
else
  call write_gpt_ele(gpt_file_unit, ele, name, fieldgrid_names, dimensions = param%fieldmap_dimension)
endif

end subroutine write_to_file

end subroutine write_gpt_lattice_file

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 

subroutine write_gpt_ele(iu, ele, name, fieldgrid_names, dimensions, only_phasing)

use bmad_interface

type (ele_struct) :: ele
type (floor_position_struct) :: floor1, floor2
type (str_index_struct), optional :: fieldgrid_names
type (branch_struct), pointer :: branch

real(rp) :: w_mat(3,3)
real(rp) :: s, x(3), dx(3), d(3), d1(2), d2(2), d3(2), d4(2), w1, e_angle, ds_slice
real(rp) :: gap, strength, x0, y0, z0, x_vec(3), y_vec(3), x_center, y_center, z_center, theta_center
real(rp) :: max_field, freq, phase_lag, ref_time

integer :: i, iu, n_slice, i_slice, q_sign, dimensions

character(40)   :: fieldgrid_output_name, name
character(1) :: component
character(8), parameter :: rfmt = '(es15.7)'
character(*), parameter  :: r_name = 'write_gpt_ele'

logical, optional :: only_phasing

!

branch => pointer_to_branch(ele)
q_sign = sign(1,  charge_of(branch%param%particle) ) 

! Get global position and rotation of origin of the field map, at the center of the element

floor1%r = [0.0_rp, 0.0_rp, ele%value(L$)/2]
floor2 = coords_local_curvilinear_to_floor (floor1, ele, in_body_frame = .true., w_mat = w_mat)
x0 = floor2%r(1)
y0 = floor2%r(2)
z0 = floor2%r(3)
! Get x, y unit vectors from w_mat
x_vec = w_mat(:,1)
y_vec = w_mat(:,2)

write (iu, '(2a)')       '#  Element: ', trim(name)
write (iu, '(2a)')       '#  Bmad name: ', trim(ele%name)
write (iu, '(2a)')       '#  Bmad key: ', key_name(ele%key)
write (iu, '(a, i0)')    '#  ix_ele: ', ele%ix_ele
write (iu, '(a, f10.5)')    '#  theta: ', floor2%theta


select case (ele%key)

!---------------------------------------------------------------------
! Example: 1D cavity
!---------------------------------------------------------------------
! XCTB01 = 0.00;
! YCTB01 = 0.00;
! ZCTB01 = 467456456456;
! Map1D_TM("wcs", XCTB01,YCTB01,ZCTB01,  1,0,0, 0,1,0, "/home/colwyn/ued/fields/eindhoven_1D.gdf", 
!                                                           "Z", "Ez", ECTB01, phiCTB01, 2*pi*Master_RF);
!
! arguments:
!  (1     )  coordinate system - specify which coordinate system you are going to use to place the element in. 
!              If you don't want to make user defined systems, then use "wcs" for world coordinate system [x,y,z]
!  (2 -  4)  offset vector - the 3D offset of the origin of the field map [ (0,0,0) in the field map coordinate system ]
!  (5 -  7)  x-unit vector in the wcs basis for the element, in above example, there is no rotation so this is (1,0,0).  
!              This doesn't require the vector to be normalized.
!  (8 - 10)  y-unit vector in the wcs basis for the element, in above example, there is no rotation so this is (0,1,0).  
!              This doesn't require the vector to be normalized. 
!
! From these two unit vectors, GPT computes the z-unit vector = (x-unit) X (y-unit), which forms the rotation matrix 
! internally used to rotate the fields in the map into the wcs.
!
!  (11) "/home/colwyn/ued/fields/eindhoven_1D.gdf" - field map file 
!  (12-13) "Z", "Ez" are the names GPT will look for in the GDF file for the columns specifying the z, Ez data.  
!    You can call your columns whatever you want, but the order you see in the function call is important: it will 
!    assign any column named "Z" to the z data, which is the 12th arguement.
!  (14) scale factor - ECTB01 is the number that GPT will scale the field by.
!  (15) Phase in radians to apply to the map = phi_concrests + phi_relative_offset.
!  (16) Angular frequency of the cavity = 2*pi*f.  (f=1.3 GHz)

case (lcavity$, rfcavity$, e_gun$) 
  if (.not. present(fieldgrid_names)) call out_io (s_error$, r_name, 'fieldgrid_names required')
  
  ! Frequency
  if ( associated(ele%grid_field)) then
    freq = ele%value(rf_frequency$) * ele%grid_field(1)%harmonic
  else
    ! No grid present
    freq = ele%value(rf_frequency$) 
  endif
  
  ! Phase 
  !phase_lag = twopi*(ele%value(phi0$) +  ele%value(phi0_err$))
  !! adjust the lag for 'zero-crossing' 
  !if (ele%key == rfcavity$) phase_lag = twopi*( ele%grid_field%mode(1)%phi0_ref - ele%value(phi0_max$) ) - phase_lag
  !ref_time

  if (logic_option(.false., only_phasing)) then
    call gpt_field_grid_scaling(ele, dimensions, max_field, ref_time)
    call write_property('field_scale', max_field, 'Maximum on-axis Ez in V/m')
    call write_property('phase', ref_time*twopi*freq)
    return
  endif

  call get_gpt_fieldgrid_name_and_scaling(ele, fieldgrid_names, fieldgrid_output_name, max_field, ref_time, dimensions)
  call write_property('frequency', freq)
  call write_property('field_scale', max_field, 'Maximum on-axis Ez in V/m')
  call write_property('phase', ref_time*twopi*freq)
  
  select case (dimensions)
  case(1)
    write (iu, '(a)', advance='NO') 'Map1D_TM('
    call write_geometry()
    write (iu, '(3a)', advance='NO')  '"', trim(fieldgrid_output_name)//'.gdf', '", '
    write (iu, '(1a)', advance='NO') ' "z", "Ez", '
    call write_property('field_scale')
    call write_comma()
    call write_property('phase')
    call write_comma()
    call write_property('frequency*2*pi')
    write (iu, '(a)') ');'  
    
  case(2)
    write (iu, '(a)', advance='NO') 'Map25D_TM('
    call write_geometry()
    write (iu, '(3a)', advance='NO')  '"', trim(fieldgrid_output_name)//'.gdf', '", '
    write (iu, '(1a)', advance='NO') ' "r", "z", "Er", "Ez", "Bphi", '
    call write_property('field_scale')
    call write_comma()
    write (iu, '(i0)', advance='NO') 0 ! Longitudinal wavenumber k in [m-1]. Must be zero for the simulation of a standing wave cavity.
    call write_comma()
    call write_property('phase')
    call write_comma()
    call write_property('frequency*2*pi')
    write (iu, '(a)') ');'      
  
  case(3)
    ! Note: GPT fields oscillate as exp(+i wt)
    do i=1,2
    if (i==1) then
      component = 'E'
    else
      if (ele%key==e_gun$ .and. freq==0) exit ! No H field for dc guns!
      component = 'H'
    endif
    
    write (iu, '(a)', advance='NO') 'Map3D_'//component//'complex('
    call write_geometry()
    write (iu, '(3a)', advance='NO')  '"', trim(fieldgrid_output_name)//'_'//component//'.gdf', '", '
    if (i==1) then
      write (iu, '(1a)', advance='NO')  '"x", "y", "z", "ExRe", "EyRe", "EzRe", "ExIm", "EyIm", "EzIm", '
    else 
      write (iu, '(1a)', advance='NO')  '"x", "y", "z", "HxRe", "HyRe", "HzRe", "HxIm", "HyIm", "HzIm", '
    endif
    
    call write_property('field_scale')
    call write_comma()
    call write_property('phase')
    call write_comma()
    call write_property('frequency*2*pi')
    write (iu, '(a)') ');'    
    enddo
    
  end select

!--------------------------------

case (solenoid$) 
  if (.not. present(fieldgrid_names)) call out_io (s_error$, r_name, 'fieldgrid_names required')
  call get_gpt_fieldgrid_name_and_scaling( ele, fieldgrid_names, fieldgrid_output_name, max_field, ref_time, dimensions)
  
  call write_property('field_scale', max_field, 'signed abs max on-axis Bz in T')
  select case (dimensions)
  case(1)
    write (iu, '(a)', advance='NO') 'Map1D_B('
    call write_geometry()
    write (iu, '(3a)', advance='NO')  '"', trim(fieldgrid_output_name)//'.gdf', '", '
    write (iu, '(1a)', advance='NO') '"z", "Bz", '
    call write_property('field_scale')
    write (iu, '(a)') ');'    

  case(2)
    write (iu, '(a)', advance='NO') 'Map2D_B('
    call write_geometry()
    write (iu, '(3a)', advance='NO')  '"', trim(fieldgrid_output_name)//'.gdf', '", '
    write (iu, '(1a)', advance='NO') ' "r", "z", "Br", "Bz", '
    call write_property('field_scale')
    write (iu, '(a)') ');'    

  case(3)
    write (iu, '(a)', advance='NO') 'Map3D_B('
    call write_geometry()
    write (iu, '(3a)', advance='NO')  '"', trim(fieldgrid_output_name)//'_B.gdf', '", '
    write (iu, '(1a)', advance='NO') ' "x", "y", "z", "Bx", "By", "Bz", '
    call write_property('field_scale')
    write (iu, '(a)') ');'    
  
  end select

!--------------------------------

case (sbend$) 
  if (.not. present(fieldgrid_names)) call out_io (s_error$, r_name, 'fieldgrid_names required')
     !!!call get_gpt_fieldgrid_name_and_scaling( ele, fieldgrid_names, fieldgrid_output_name, &
     !!!                            max_field, ref_time, dimensions = dimensions)
    call out_io (s_error$, r_name, 'TRANSLATION TO GPT FOR BEND NOT YET IMPLEMENTED!', 'SIMULATION WILL BE INVALID')
    write(*, *) '#TEMP: NO FIELD MAP, GEOMETRY: '
    write(iu, '(a)', advance = 'NO') '#GEOMETRY: '
    call write_geometry()
    !call write_property('field_scale', max_field, 'signed abs max on-axis By in T')
    !select case (dimensions)
    !case(3)
    !  write (iu, '(a)', advance='NO') 'Map3D_B('
    !  call write_geometry()
    !  write (iu, '(3a)', advance='NO')  '"', trim(fieldgrid_output_name)//'_B.gdf', '", '
    !  write (iu, '(1a)', advance='NO') ' "x", "y", "z", "Bx", "By", "Bz", '
    !  call write_property('field_scale')
    !  write (iu, '(a)') ');'    
  

!--------------------------------

case (quadrupole$)   
   call write_property('gradient', ele%value(b1_Gradient$), 'T/m')
   call write_property('L', ele%value(L$), 'm') 
   write (iu, '(a)', advance='NO') 'quadrupole('
   call write_geometry()
   call write_property('L')
   call write_comma()
   call write_property('gradient')
   write (iu, '(a)') ' );'  

case (drift$, monitor$, instrument$, pipe$, girder$, overlay$, group$, marker$)
  ! Nothing to do

case default
  call out_io (s_error$, r_name, 'No translation to GPT for:' // ele%name, 'Of type: ' // key_name(ele%key))

end select

! Footer
write (iu, '(a)') '' 

!-------------------------------------------------------------------------------------

contains

subroutine write_geometry()
write (iu, '(a)', advance = 'NO') '"wcs", '
call write_real(x0); call write_comma()
call write_real(y0); call write_comma()
call write_real(z0); call write_comma()
call write_real(w_mat(1,1)); call write_comma()
call write_real(w_mat(2,1)); call write_comma()
call write_real(w_mat(3,1)); call write_comma()
call write_real(w_mat(1,2)); call write_comma()
call write_real(w_mat(2,2)); call write_comma()
call write_real(w_mat(3,2)); call write_comma()
end subroutine write_geometry

!-------------------------------------------------------------------------------------
! contains

subroutine write_property (property, value, comment)
character(*) :: property
real(rp), optional :: value
character(*), optional :: comment
if (present(value)) then
  if (value == nint(value)) then
    write (iu, '(a, i0, a)', advance = 'NO') trim(name)//'_'//trim(property)// ' = ', nint(value), ';'   
  else
    write (iu, '(a, '//rfmt//', a)', advance = 'NO') trim(name)//'_'//trim(property)// ' = ', value, ';'
  endif
  
  if (present(comment)) then
    write (iu, '(a)') ' # '//trim(comment)
  else
    write (iu, '(a)') ''
  endif
else
  write (iu, '(a)', advance = 'NO') trim(name)//'_'//trim(property)
endif
end subroutine write_property

!-------------------------------------------------------------------------------------
! contains

subroutine write_comma()
write (iu, '(a)', advance='NO') ', '
end subroutine write_comma

subroutine write_real(x)
real(rp) :: x
if (x == nint(x)) then
  write (iu, '(i0 )', advance = 'NO') nint(x)
else
  write (iu, rfmt, advance = 'NO') x
endif
end subroutine write_real

end subroutine write_gpt_ele

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine get_gpt_fieldgrid_name_and_scaling(ele, name_indexx, output_name, field_scale, ref_time, dimensions)
!
! Subroutine to get a field grid filename and its scaling. Calls write_gpt_field_grid_file.
!   If the field grid file does not exist, it is written.
!
!   Note: This is very similar to get_opal_fieldgrid_name_and_scaling
!
!
! Input:
!   ele              -- ele_struct: element to make map
!   name_indexx      -- str_index_struct: contains field grid filenames
!   dimensions       -- integer: 1, 2, or 3 dimensions.
!
! Output:   
!   name_indexx      -- char_indexx_struct: updated if new name is added
!   output_name      -- real(rp): output filename. 
!   field_scale      -- real(rp): the scaling of the field grid
!   ref_time         -- real(rp): time that the field was evaluated at
!-

subroutine get_gpt_fieldgrid_name_and_scaling(ele, name_indexx, output_name, field_scale, ref_time, dimensions)


type (ele_struct) :: ele
type (str_index_struct) :: name_indexx
character(*)  :: output_name
real(rp)      :: field_scale
real(rp)      :: ref_time
character(200) :: unique_grid_file

integer :: ix_match, iu_fieldgrid, ios
integer, optional :: dimensions
character(40), parameter:: r_name = 'get_gpt_fieldgrid_name_and_scaling'
character(10) :: suffix = '_ASCII.gpt'

!

output_name = ''

! Check for field grid. Force sbends to have unique fieldmaps
! because their geometry can change how the field is output.

if (associated (ele%grid_field) .and. ele%key /= sbend$) then
  unique_grid_file = ele%grid_field(1)%ptr%file
else
  ! Just us the element's index as an identifier
  write(unique_grid_file, '(a,i0)') 'ele_', ele%ix_ele
endif

! Check field map file. If file has not been written, create a new file. 

call find_index (unique_grid_file, name_indexx, ix_match)

! Check for match with existing grid

if (ix_match > 0) then
  ! File should exist  
  call get_output_name()
  call gpt_field_grid_scaling(ele, dimensions, field_scale, ref_time)
else
  ! File does not exist.
  ! Add name to list  
  call find_index (unique_grid_file, name_indexx, ix_match, add_to_list = .true.)
  ix_match = name_indexx%n_max
  call get_output_name()

  
  select case(dimensions)
  case(1)
    iu_fieldgrid = lunget()
    open(iu_fieldgrid, file = trim(output_name)//trim(suffix), iostat = ios)
    call write_gpt_field_grid_file_1D (iu_fieldgrid , ele, field_scale, ref_time)
    close(iu_fieldgrid)
  
  case(2)
    iu_fieldgrid = lunget()
    open(iu_fieldgrid, file = trim(output_name)//trim(suffix), iostat = ios)
    call write_gpt_field_grid_file_2D (iu_fieldgrid, ele, field_scale, ref_time)
    close(iu_fieldgrid)
  
  case(3)
    ! This will write several files
    call write_gpt_field_grid_file_3D (output_name, ele, field_scale, ref_time)
  
  case default
    call out_io (s_error$, r_name, 'Bad dimension request for field map ')
    if (global_com%exit_on_error) call err_exit
  end select  
  
end if

!------------------------------------------------------------------------
contains 

subroutine get_output_name()
  ! Names will be suffixed elsewhere
  write(output_name, '(i1, a, i0)') dimensions, 'D_fieldmap_', ix_match
end subroutine

end subroutine get_gpt_fieldgrid_name_and_scaling

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
! Get field_scale and ref_time depending on the type of grid dimension used

subroutine gpt_field_grid_scaling (ele, dimensions, field_scale, ref_time)

type(ele_struct) :: ele
integer :: dimensions
real(rp) :: field_scale, ref_time
character(40) :: r_name = 'gpt_field_grid_scaling'

! Call with iu=0 or '' to get field_scale and ref_time without writing a file

select case(dimensions)
  case(1)
    call write_gpt_field_grid_file_1D (0, ele, field_scale, ref_time)
  case(2)
    call write_gpt_field_grid_file_2D (0, ele, field_scale, ref_time)
  case(3)
    call write_gpt_field_grid_file_3D ('', ele, field_scale, ref_time)
  case default
    call out_io (s_error$, r_name, 'Bad dimension request for field map ')
    if (global_com%exit_on_error) call err_exit
end select

end subroutine gpt_field_grid_scaling

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine write_gpt_field_grid_file_1D (gpt_file_unit, ele, maxfield, ref_time, dz, err)
!
!   Write 1-D field map files for gpt. The format is:
!   z field
!   ...
!  
!   Note: Simplified from write_opal_field_grid_file
!
! Input:
!   gpt_file_unit -- Integer: unit number to write to, if > 0
!                        if < 0, nothing is written, and only maxfield is returned
!   ele            -- ele_struct: element to make map
!   dz             -- real(rp), optional: z step size in m. Default: 0.001 m
!
!
! Output:          
!   maxfield       -- Real(rp): absolute maximum found for element field scaling
!   ref_time       -- real(rp): time that the field was evaluated at
!   err            -- Logical, optional: Set True if, say a file could not be opened.
!-



subroutine write_gpt_field_grid_file_1D (gpt_file_unit, ele, maxfield, ref_time, dz, err)

type (coord_struct) :: orb
type(em_field_struct) :: field_re, field_im
type (grid_field_pt1_struct), allocatable :: pt(:)
type (grid_field_pt1_struct) :: ref_field
type (ele_struct) :: ele
type (lat_param_struct) :: param
type (branch_struct), pointer :: branch

real(rp) :: maxfield
real(rp) :: z_step, z_min, z_max, z0
real(rp) :: freq,  z, phase_ref
real(rp) :: gap, edge_range
real(rp) :: ref_time
real(rp) :: Ex_factor, Ez_factor, Bx_factor, By_factor, Bz_factor
real(rp), optional :: dz
complex ::  phasor_rotation

integer      :: gpt_file_unit
integer :: nx, nz, iz, ix

logical, optional :: err
logical loc_ref_frame

character(40)  :: r_name = 'write_gpt_field_grid_file'
character(10)   ::  rfmt 

!

if (present(err)) err = .true.

loc_ref_frame = .true. 

branch => pointer_to_branch(ele)
param = branch%param

! Format for numbers
rfmt = 'es13.5'

! Output will be relative to z0
z0 = ele%value(L$)/2

! Step size
z_step = real_option(0.001_rp, dz)

z_min = 0.0_rp
z_max = ele%value(L$)
 
nz = ceiling(z_max/z_step)

! These are 1D fields, so the field will be probed on-axis
orb%vec(1) = 0
orb%vec(3) = 0

Ez_factor = 1
Bz_factor = 1

! Reference time, will be adjusted by oscillating elements below  
ref_time = 0  
phase_ref = 0

select case (ele%key)

!-----------
! LCavity, RFCavity, E_GUN
!-----------
case (lcavity$, rfcavity$, e_gun$) 
  if ( associated(ele%grid_field)) then
    freq = ele%value(rf_frequency$) * ele%grid_field(1)%harmonic
  else
    ! No grid present
    freq = ele%value(rf_frequency$) 
  endif

  ! Allocate temporary pt array
  allocate (pt(0:nz))
  ! Write data points
  
  ! initialize maximum found field
  maxfield = 0
  do iz = 0, nz
    z = z_step * iz 
    ! Calculate field at \omegat*t=0 and \omega*t = \pi/2 to get real and imaginary parts
    call em_field_calc (ele, param, z, orb, loc_ref_frame, field_re, rf_time = 0.0_rp)
    ! if frequency is zero, zero out field_im
    if(freq == 0) then
      field_im%E=0
      field_im%B=0
    else 
      call em_field_calc (ele, param, z, orb, loc_ref_frame, field_im, rf_time = 0.25/freq)
    endif
    pt(iz)%E(:) = cmplx(field_re%E(:), field_im%E(:), rp)
    pt(iz)%B(:) = cmplx(field_re%B(:), field_im%B(:), rp)
    ! Update ref_field if larger Ez is found
    if(abs(pt(iz)%E(3)) > maxfield) then
      ref_field = pt(iz)
      maxfield = abs(ref_field%E(3))
    end if 
  end do
  
  ! Reference time
  if (freq /= 0) then
    ref_time = -phase_ref /(twopi*freq)
    phase_ref = atan2( aimag(ref_field%E(3) ), real(ref_field%E(3) ) )
  endif
  
  ! Write to file
  if (gpt_file_unit > 0 )  then
    ! Normalize
    if (maxfield > 0) Ez_factor = 1/maxfield

    ! Calculate complex rotation number to rotate Ez onto the real axis    
    phasor_rotation = cmplx(cos(phase_ref), -sin(phase_ref), rp)
    
    write (gpt_file_unit, '(2a13)') 'z', 'Ez'
    do iz = 0, nz   
      write (gpt_file_unit, '(2'//rfmt//')') z_step * iz - z0, Ez_factor * real ( pt(iz)%E(3) * phasor_rotation )
    enddo 
  end if
  deallocate(pt)

  !-----------
  ! Solenoid
  !-----------
  ! Note: This is similar to code for lcavity/rfcavity
  case (solenoid$) 
                                         
  ! Allocate temporary pt array
  allocate (pt(0:nz))
  ! Write data points
  
  ! initialize maximum found field
  maxfield = 0
  
  do iz = 0, nz
    z = z_step * iz 
    call em_field_calc (ele, param, z, orb, loc_ref_frame, field_re, rf_time = 0.0_rp)
    field_im%E = 0
    field_im%B = 0
    pt(iz)%E(:) = cmplx(field_re%E(:), field_im%E(:), rp)
    pt(iz)%B(:) = cmplx(field_re%B(:), field_im%B(:), rp)

    ! Update ref_field if larger Bz is found
    if (abs(pt(iz)%B(3)) > maxfield) then
         ref_field = pt(iz)
         maxfield = abs(ref_field%B(3))
    end if 
  end do
  
  ! Restore the sign
  maxfield = ref_field%B(3)
  
  ! Write to file
  if (gpt_file_unit > 0 )  then
    write (gpt_file_unit, '(2a13)') 'z', 'Bz'

    ! Scaling
    if (maxfield /= 0) Bz_factor = 1/maxfield
    
    do iz = 0, nz
      write (gpt_file_unit, '(2'//rfmt//')') z_step * iz - z0, Bz_factor * real ( pt(iz)%B(3))          
    enddo
  end if

  ! cleanup 
  deallocate(pt)

  !-----------
  ! Default (gives an error)
  !-----------
  case default
  call out_io (s_error$, r_name, 'MISSING gpt FIELD GRID CODE FOR ELEMENT TYPE: ' // key_name(ele%key), &
             '----')
  if (global_com%exit_on_error) call err_exit
  
end select 

if (maxfield == 0) then
  call out_io (s_error$, r_name, 'ZERO MAXIMUM FIELD IN ELEMENT: ' // key_name(ele%key), &
             '----')
  if (global_com%exit_on_error) call err_exit
end if

end subroutine write_gpt_field_grid_file_1D


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine write_gpt_field_grid_file_2D (gpt_file_unit, ele, maxfield, ref_time, dr, dz,  err)
!
! Subroutine to write an GPT lattice file using the information in
! a lat_struct. Optionally only part of the lattice can be generated.
!
!
! Input:
!   gpt_file_unit -- Integer: unit number to write to, if > 0
!                        if < 0, nothing is written, and only maxfield is returned
!   ele            -- ele_struct: element to make map
!   dr             -- real(rp), optional: r step size in m. Default: 0.001 m
!   dz             -- real(rp), optional: z step size in m. Default: 0.001 m
!   r_max          -- real(rp), optional: maximum radius in m. Default: 0.02 m
!
! Output:          
!   maxfield       -- Real(rp): absolute maximum found for element field scaling
!   ref_time       -- real(rp): time that the field was evaluated at
!   err            -- Logical, optional: Set True if, say a file could not be opened.
!-

subroutine write_gpt_field_grid_file_2D (gpt_file_unit, ele, maxfield, ref_time, dr, dz, r_max,  err)

type (coord_struct) :: orb
type(em_field_struct) :: field_re, field_im
type (grid_field_pt1_struct), allocatable :: pt(:,:,:)
type (grid_field_pt1_struct) :: ref_field
type (ele_struct) :: ele
type (lat_param_struct) :: param
type (branch_struct), pointer :: branch

real(rp) :: maxfield, ref_time
real(rp), optional :: dr, dz, r_max
real(rp) :: x_step, z_step, x_min, x_max, z_min, z_max, z0
real(rp) :: freq, x, z, phase_ref
real(rp) :: gap, edge_range
complex ::  phasor_rotation

integer :: gpt_file_unit
integer :: dimensions
integer :: nx, nz, iz, ix

logical, optional :: err
logical loc_ref_frame

character(40)  :: r_name = 'write_gpt_field_grid_file_2D '
character(10), parameter ::  rfmt  =  'es13.5'

real(rp) :: Ex_factor, Ez_factor, Bx_factor, By_factor, Bz_factor

!

if (present(err)) err = .true.

loc_ref_frame = .true. 

branch => pointer_to_branch(ele)
param = branch%param

! Output will be relative to z0
z0 = ele%value(L$)/2
  
! Step size
x_step = real_option(0.001_rp, dr)
z_step = real_option(0.001_rp, dz)

z_min = 0.0_rp
z_max = ele%value(L$)
 
nz = ceiling(z_max/z_step)

x_min = 0.0_rp
x_max = real_option(0.02_rp, r_max)

z_min = 0.0_rp
z_max = ele%value(L$)

nx = ceiling(x_max/x_step)  
nz = ceiling(z_max/z_step)


! Reference time, will be adjusted by oscillating elements below  
ref_time = 0  
phase_ref = 0

select case (ele%key)

!-----------
! LCavity, RFCavity, E_GUN

case (lcavity$, rfcavity$, e_gun$) 

  freq = ele%value(rf_frequency$) * ele%grid_field(1)%harmonic

  ! Example:
  ! r  z  Er  Ez  Bphi
  ! ...

  ! Allocate temporary pt array
  allocate (pt(0:nx, 0:nz, 1:1))
  ! Write data points
  
  ! initialize maximum found field
  maxfield = 0
  
  do ix = 0, nx
    do iz = 0, nz
      x = x_step * ix
      z = z_step * iz 
      orb%vec(1) = x
      orb%vec(3) = 0.0_rp
      
      ! Calculate field at \omegat*t=0 and \omega*t = \pi/2 to get real and imaginary parts
      call em_field_calc (ele, param, z, orb, loc_ref_frame, field_re, rf_time = 0.0_rp)
      ! if frequency is zero, zero out field_im
      if(freq == 0) then
        field_im%E=0
        field_im%B=0
      else 
        call em_field_calc (ele, param, z, orb, loc_ref_frame, field_im, rf_time = 0.25/freq)
      endif

      pt(ix, iz, 1)%E(:) = cmplx(field_re%E(:), field_im%E(:), rp)
      pt(ix, iz, 1)%B(:) = cmplx(field_re%B(:), field_im%B(:), rp)
      
      ! Update ref_field if larger Ez is found
      if(ix == 0 .and. abs(pt(ix, iz, 1)%E(3)) > maxfield) then
         ref_field = pt(ix, iz, 1)
         maxfield = abs(ref_field%E(3))
      end if 
    end do
  end do
  
  ! Reference time
  if (freq /= 0) then
    ref_time = -phase_ref /(twopi*freq)
    phase_ref = atan2( aimag(ref_field%E(3) ), real(ref_field%E(3) ) )
  endif
  
  ! Write to file
  if (gpt_file_unit > 0)  then
    write (gpt_file_unit, '(5a13)') 'r', 'z', 'Er', 'Ez', 'Bphi'
    ! Scaling
    Ex_factor =  (1/maxfield)
    Ez_factor =  (1/maxfield)
    By_factor = -(1/maxfield)
  
    ! Calculate complex rotation number to rotate Ez onto the real axis
    phasor_rotation = cmplx(cos(phase_ref), -sin(phase_ref), rp)
  
    do ix = 0, nx
      do iz = 0, nz
        write (gpt_file_unit, '(5'//rfmt//')')  x_step * ix, z_step * iz - z0 , &
                            Ex_factor * real ( pt(ix, iz, 1)%E(1) * phasor_rotation ), &
                                                Ez_factor * real ( pt(ix, iz, 1)%E(3) * phasor_rotation ), &
                                                By_factor * aimag (  pt(ix, iz, 1)%B(2)*phasor_rotation )
      enddo
    enddo
  
  end if
   
  deallocate(pt)

!-----------
! Solenoid

! Note: This is similar to code for lcavity/rfcavity
! Example:
! r z Br Bz

case (solenoid$) 
                                       
  ! Allocate temporary pt array
  allocate ( pt(0:nx, 0:nz, 1:1) )
  ! Write data points

  ! initialize maximum found field
  maxfield = 0

  do ix = 0, nx
    do iz = 0, nz
      x = x_step * ix
      z = z_step * iz 
      orb%vec(1) = x
      orb%vec(3) = 0.0_rp

      call em_field_calc (ele, param, z, orb, loc_ref_frame, field_re, rf_time = 0.0_rp)
      field_im%E = 0
      field_im%B = 0

      pt(ix, iz, 1)%E(:) = cmplx(field_re%E(:), field_im%E(:), rp)
      pt(ix, iz, 1)%B(:) = cmplx(field_re%B(:), field_im%B(:), rp)
    
      ! Update ref_field if larger Bz is found
      if (ix==0 .and. abs(pt(ix, iz, 1)%B(3)) > maxfield) then
         ref_field = pt(ix, iz, 1)
         maxfield = abs(ref_field%B(3))
      end if

    end do
  end do

  ! Restore the sign
  maxfield = ref_field%B(3)

  ! Write to file
  if (gpt_file_unit > 0 )  then

    ! Write header
    write (gpt_file_unit, '(4a13)') 'r', 'z', 'Br', 'Bz'

    ! Scaling
    Bx_factor = 1/maxfield
    Bz_factor = 1/maxfield
    
    ! XZ ordering: ix changes fastest (inner loop)
    do ix = 0, nx
      do iz = 0, nz
        write (gpt_file_unit, '(4'//rfmt//')')  x_step * ix, z_step * iz - z0, &
                 Bx_factor * real (  pt(ix, iz, 1)%B(1) ), &
                 Bz_factor * real (  pt(ix, iz, 1)%B(3) )            
      enddo
    enddo

  end if

  ! cleanup 
  deallocate(pt)

!-----------
! Default (gives an error)

case default
  call out_io (s_error$, r_name, 'MISSING GPT FIELD GRID CODE FOR: ' // key_name(ele%key), '----')
    if (global_com%exit_on_error) call err_exit
end select 

!-----------------------------------------------------------

if (maxfield == 0) then
  call out_io (s_error$, r_name, 'ZERO MAXIMUM FIELD IN ELEMENT: ' // key_name(ele%key), '----')
  if (global_com%exit_on_error) call err_exit
end if

end subroutine write_gpt_field_grid_file_2D


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine write_gpt_field_grid_file_3D (base_filename, ele, maxfield, ref_time, dz, err)
!
!   Writes 3-D field map files for gpt. The format is:
!
!   E-fields:
!   'x', 'y', 'z', 'ExRe', 'EyRe', 'EzRe', 'ExIm ', 'EyIm ', 'EzIm '
!   H-fields
!   'x', 'y', 'z', 'HxRe', 'HyRe', 'HzRe', 'HxIm ', 'HyIm ', 'HzIm '
!
!   where the fields oscillate as exp(+i \omega t)
!
!   Note: similar to write_gpt_field_grid_file
!
! Input:
!   base_filename  -- character(*): Base filename. Files will be written as:
!                                   base_filename_E_ASCII.gpt, _H_ASCII.gpt
!                                   If set to '', no files will be written
!   ele            -- ele_struct: element to make map
!   dz             -- real(rp), optional: z step size in m. Default: 0.001 m
!
! Output:          
!   maxfield       -- Real(rp): absolute maximum on-axis field found for element field scaling
!   ref_time       -- real(rp): time that the field was evaluated at
!   err            -- Logical, optional: Set True if, say a file could not be opened.
!-

subroutine write_gpt_field_grid_file_3D (base_filename, ele, maxfield, ref_time, dz, err)

implicit none

character(*) :: base_filename
integer      :: iu
type (ele_struct) :: ele
type (lat_param_struct) :: param

logical, optional :: err

character(40)  :: r_name = 'write_gpt_field_grid_file'
character(10)   ::  rfmt 

type (coord_struct) :: orb
type(em_field_struct) :: field_re, field_im
type (grid_field_pt1_struct), allocatable :: pt(:,:,:)
type (grid_field_pt1_struct) :: ref_field
type (branch_struct), pointer :: branch

real(rp)        :: maxfield, ref_time
real(rp) :: x_step, x_min, x_max
real(rp) :: y_step, y_min, y_max
real(rp) :: z_step, z_min, z_max, z0
real(rp) :: freq, x, y, z, phase_ref
real(rp) :: gap, edge_range
real(rp), optional :: dz
complex ::  phasor_rotation, test_field_value

integer :: nx, ny, nz, iz, ix, iy, ifield, ix_center, iy_center

logical loc_ref_frame, write_file

!

if (present(err)) err = .true.

if (base_filename == '') then
  write_file = .false.
else
  write_file = .true.
endif

loc_ref_frame = .true. 

branch => pointer_to_branch(ele)
param = branch%param

! Format for numbers
rfmt = 'es13.5'

! Output will be relative to z0
z0 = ele%value(L$)/2

z_step = real_option(0.001_rp, dz)

z_min = 0.0_rp
z_max = ele%value(L$)

! Special bounds for curvilinear coordinates
if (ele%key == sbend$) then
  if (ele%value(g$) /= 0) then
    z_max = ele%value(rho$)*sin(ele%value(angle$)/2)
    z_min = -z_max
  endif
endif


x_max =  .02_rp
x_step = z_step  ! Same as z TODO: generalize
x_min = -x_max
ix_center = nint(x_max/x_step) + 1
 
y_max =  .02_rp
y_step = z_step  ! Same as z 
y_min = -x_max
iy_center = nint(y_max/y_step) + 1
 
nz = ceiling((z_max-z_min)/z_step) + 1
nx = ceiling((x_max-x_min)/x_step) + 1
ny = ceiling((y_max-y_min)/y_step) + 1

! Reference time, will be adjusted by oscillating elements below  

ref_time = 0 

! Get frequency

freq = 0
select case (ele%key)
case (lcavity$, rfcavity$, e_gun$) 
  if ( associated(ele%grid_field)) then
    freq = ele%value(rf_frequency$) * ele%grid_field(1)%harmonic
  else
    ! No grid present
    freq = ele%value(rf_frequency$) 
  endif
end select

! Allocate all sample points.

allocate (pt(1:nx, 1:ny, 1:nz))

maxfield = 0

! 3D loop to load pt array, and to get the maximum field

do iz=1, nz
do iy=1, ny

! Get field along x line
do ix=1, nx
  orb%vec(1) = x_min + (ix-1)*x_step
  orb%vec(3) = y_min + (iy-1)*y_step
  z = z_min + (iz-1)*z_step
  if (ele%key == sbend$) then
    call sbend_field_at(orb%vec(1), z, field_re)
  else
    ! Calculate field at \omegat*t=0 and \omega*t = \pi/2 to get real and imaginary parts
    call em_field_calc (ele, param, z, orb, loc_ref_frame, field_re, rf_time = 0.0_rp)
  endif
  ! if frequency is zero, zero out field_im
  if(freq == 0) then
    field_im%E=0
    field_im%B=0
  else 
    call em_field_calc (ele, param, z, orb, loc_ref_frame, field_im, rf_time = 0.25/freq)
  endif
  pt(ix, iy, iz)%E(:) = cmplx(field_re%E(:), field_im%E(:), rp)
  pt(ix, iy, iz)%B(:) = cmplx(field_re%B(:), field_im%B(:), rp)

enddo
enddo
enddo

! Get reference field: The largest on-axis field

do iz = 1, nz
  test_field_value = gpt_max_field_reference(pt(ix_center, iy_center, iz), ele)
  if(abs(test_field_value) > maxfield) then
    ref_field = pt(ix_center, iy_center, iz)
    maxfield = abs(test_field_value)
  end if 
enddo

if (freq == 0) then
  ! Fields should be purely real
  ! restore the sign
  maxfield = sign(maxfield, real(gpt_max_field_reference(ref_field, ele), rp))
  phasor_rotation = cmplx(1.0_rp, 0.0_rp, rp)
else
  ! Calculate complex rotation number to rotate Ez onto the real axis
  phase_ref = atan2( aimag(ref_field%E(3) ), real(ref_field%E(3) ) )
  phasor_rotation = cmplx(cos(phase_ref), -sin(phase_ref), rp)
  ref_time = -phase_ref /(twopi*freq)
endif

! Put scaling in phasor rotation

if (maxfield /= 0) phasor_rotation = (1/maxfield)*phasor_rotation

! Write to file according to element

if (write_file) then
  select case(ele%key)
  
  case(e_gun$)
    call write_pts('E')
    if (freq /= 0) call write_pts('H')
  
  case(lcavity$, rfcavity$)
    call write_pts('E')
    call write_pts('H')
    
  case(sbend$)
    ! Curvilinear, remove z0 offset, because z_max and z_min already include it
    z0 = 0
    call write_pts('B')
    
  case(solenoid$)  
    call write_pts('B')
    
  case(quadrupole$)  
    call write_pts('B')    
    
  case default
    print *, 'not coded'
  end select
endif

! Cleanup

deallocate(pt)

!-----------------------------------------
contains

! write to file: according to component (1 for E, 2 for H)

subroutine write_pts(component)
implicit none
character(40) :: output_name
character(1) :: component
integer :: iu, ios

!

write(output_name, '(a)') trim(base_filename)//'_'//component

iu =lunget()
open(iu, file = trim(output_name)//'_ASCII.gpt', iostat = ios)

! Note: GPT fields oscillate as exp(+i wt), so the imaginary parts need a minus sign.

select case(component)
case('E')
! E field
  
  write(iu, '(9a13)') 'x', 'y', 'z', 'ExRe', 'EyRe', 'EzRe', 'ExIm ', 'EyIm ', 'EzIm '
  do iz=1, nz
  do iy=1, ny
  do ix=1, nx-1
    write(iu, '(9'//rfmt//')') x_min + (ix-1)*x_step, y_min + (iy-1)*y_step, z_min + (iz-1)*z_step - z0, &
                               real( pt(ix, iy, iz)%E(:) * phasor_rotation), &
                              -aimag( pt(ix, iy, iz)%E(:) * phasor_rotation)                        
  enddo
  enddo
  enddo
  
case('H')
  ! H field
  write(iu, '(9a13)') 'x', 'y', 'z', 'HxRe', 'HyRe', 'HzRe', 'HxIm ', 'HyIm ', 'HzIm '
  do iz=1, nz
  do iy=1, ny
  do ix=1, nx-1
    write(iu, '(9'//rfmt//')') x_min + (ix-1)*x_step, y_min + (iy-1)*y_step, z_min + (iz-1)*z_step - z0, &
                               real( pt(ix, iy, iz)%B(:) * phasor_rotation)/mu_0_vac, &
                               -aimag( pt(ix, iy, iz)%B(:) * phasor_rotation)/mu_0_vac         
  enddo
  enddo
  enddo
  
case('B')
  ! B field
  write(iu, '(6a13)') 'x', 'y', 'z', 'Bx', 'By', 'Bz'
  do iz=1, nz
  do iy=1, ny
  do ix=1, nx-1
    write(iu, '(6'//rfmt//')') x_min + (ix-1)*x_step, y_min + (iy-1)*y_step, z_min + (iz-1)*z_step - z0, &
                               real( pt(ix, iy, iz)%B(:) * phasor_rotation)     
  enddo
  enddo
  enddo
  
end select

close(iu)

end subroutine write_pts

!-----------------------------------------
! contains

subroutine sbend_field_at(x, z, field)
type(em_field_struct) :: field
real(rp) :: x, z, s

!s is relative to the origin of the Cartesian system. Add z0 to get s relative
!to the element
call convert_local_cartesian_to_local_curvilinear(x, z, ele%value(g$),  orb%vec(1), s)

if (s+z0>ele%value(L$) .or. s+z0 < 0) then
  field%e = 0
  field%b = 0
  return
endif  
  
call em_field_calc (ele, param, s + z0, orb, loc_ref_frame, field, rf_time = 0.0_rp)

call rotate_field_zx(field, -s*ele%value(g$))
end subroutine sbend_field_at

end subroutine write_gpt_field_grid_file_3D

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
!  Returns the relevant field value for the field scaling calculation
!-

function gpt_max_field_reference (pt0, ele) result (field_value)

implicit none

type (grid_field_pt1_struct) pt0 
type (ele_struct) :: ele
real(rp) :: field_value
character(40), parameter :: r_name = 'gpt_max_field_reference'

!

select case (ele%key)
case(e_gun$, lcavity$, rfcavity$)
  ! Ez
  field_value = pt0%E(3) ! Note: pt0 is complex
case(sbend$)
  ! By
  field_value = pt0%B(2)
case(solenoid$)
  ! Bz
  field_value = pt0%B(3)
  
case default
  call out_io(s_error$, r_name, 'Max field reference not specified for element type '// key_name(ele%key))

end select

end function gpt_max_field_reference 

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

subroutine rotate_field_zx (field, theta)

type (em_field_struct) :: field
real(rp) :: theta, temp, ca, sa

!

ca = cos(theta)
sa = sin(theta)

temp       = field%E(3)*ca - field%E(1)*sa
field%E(1) = field%E(3)*sa + field%E(1)*ca
field%E(3) = temp

temp       = field%B(3)*ca - field%B(1)*sa
field%B(1) = field%B(3)*sa + field%B(1)*ca
field%B(3) = temp

end subroutine rotate_field_zx

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

subroutine convert_local_curvilinear_to_local_cartesian(x, s, g, xout, zout)
real(rp) :: x, s, g, rho, xout, zout, theta
if (g == 0) then
  xout = x
  zout = s
  return
endif
rho = 1/g
theta = s/rho
xout = (x + rho)*cos(theta) + rho
zout = (x+rho)*sin(theta)
end subroutine convert_local_curvilinear_to_local_cartesian

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

subroutine convert_local_cartesian_to_local_curvilinear (x, z, g, xout, sout)
real(rp) :: x, z, g, rho, xout, sout, theta
if (g == 0) then
  xout = x
  sout = z
  return
endif
rho = 1/g
theta = atan2(z, abs(x + rho)) 
sout = theta*abs(rho)
if (abs(theta) < pi/2) then
  xout = (x + rho)/cos(theta) - rho
else
  xout = z/sin(theta) - rho
endif

end subroutine convert_local_cartesian_to_local_curvilinear

end module
