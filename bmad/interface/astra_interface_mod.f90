module astra_interface_mod

use bmad_interface

implicit none

type astra_lattice_param_struct
  integer :: fieldmap_dimension = 1    ! Dimensions for field map. 1 or 3
end type

contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine write_astra_lattice_file (astra_file_unit, lat, astra_lattice_param, err)
!
! Subroutine to write an Astra lattice file using the information in a lat_struct.
!
! Input:
!   astra_file_unit -- Integer: unit number to write to
!   lat             -- lat_struct: Holds the lattice information.
!
! Output:
!   err    -- Logical, optional: Set True if, say a file could not be opened.
!-

subroutine write_astra_lattice_file (astra_file_unit, lat, astra_lattice_param,  err)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (str_index_struct) :: fieldgrid_names
type (astra_lattice_param_struct) :: astra_lattice_param

integer :: n, ie, ix_start, ix_end, iu, id
integer :: astra_file_unit

character(4000)  :: line

logical :: err

character(40)  :: r_name = 'write_astra_lattice_file'

!

print *, 'This is '//r_name

ix_start = 1
ix_end = lat%n_ele_max ! Include lords


! Initialize fieldgrid filename list
n = ix_end - ix_start + 1
fieldgrid_names%n_max = 0

! Loop over bends
write(astra_file_unit, '(a)') '&DIPOLE'
id = 0
! No need to iterate id, the routine below will do so
do ie = ix_start, ix_end
  ele => lat%ele(ie)
  if (skip_ele()) cycle
  ! Point to first slave for multipass lords
  if (ele%lord_status == multipass_lord$) ele => pointer_to_slave (ele, 1)
  if (ele%key==sbend$) then
    if (id == 0) write(astra_file_unit, '(a)') '  LDIPOLE = T'
    call write_astra_ele(astra_file_unit, ele, id)
  endif
enddo
write(astra_file_unit, '(a)') '/'

! Loop over quadrupoles
write(astra_file_unit, '(a)') '&QUADRUPOLE'

id = 0
do ie = ix_start, ix_end
  ele => lat%ele(ie)
  if (skip_ele()) cycle
  if (ele%lord_status == multipass_lord$) ele => pointer_to_slave (ele, 1)
  if (ele%key==quadrupole$) then
    if (id == 0) write(astra_file_unit, '(a)') '  LQuad = T'
    call write_astra_ele(astra_file_unit, ele, id)
    
  endif
enddo
write(astra_file_unit, '(a)') '/'

! Loop over solenoids
write(astra_file_unit, '(a)') '&SOLENOID'

id = 0
do ie = ix_start, ix_end
  ele => lat%ele(ie)
  if (skip_ele()) cycle
  if (ele%lord_status == multipass_lord$) ele => pointer_to_slave (ele, 1)
  if (ele%key == solenoid$)  then
    if (id == 0) write(astra_file_unit, '(a)') '  LBfield = T'
    call write_astra_ele(astra_file_unit, ele, id, fieldgrid_names, dimensions = astra_lattice_param%fieldmap_dimension)
  endif
enddo
write(astra_file_unit, '(a)') '/'


! Loop over cavities
write(astra_file_unit, '(a)') '&CAVITY'

id = 0
do ie = ix_start, ix_end
  ele => lat%ele(ie)
  if (skip_ele()) cycle
  if (ele%lord_status == multipass_lord$) ele => pointer_to_slave (ele, 1)
  if (ele%key == lcavity$ .or. ele%key == rfcavity$ .or. ele%key == e_gun$)  then
    if (id == 0) write(astra_file_unit, '(a)') '  LEfield = T'
    call write_astra_ele(astra_file_unit, ele, id, fieldgrid_names, dimensions = astra_lattice_param%fieldmap_dimension)
  endif
enddo
write(astra_file_unit, '(a)') '/'

contains

!--- Function to skip non-physical elements
function skip_ele() result (skip)
implicit none
logical skip
  ! Skip these elements:
if (ele%slave_status == super_slave$ .or. &
      ele%slave_status == multipass_slave$ .or. &
      ele%key == girder$ .or. &
      ele%key == overlay$ .or. &
      ele%key == group$) then
  skip = .true.
else
  skip = .false.
endif
end function

end subroutine write_astra_lattice_file

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
subroutine write_astra_ele(iu, ele, id, fieldgrid_names, dimensions)

implicit none

type (ele_struct) :: ele
type (floor_position_struct) :: floor1, floor2
type (str_index_struct), optional :: fieldgrid_names
type (branch_struct), pointer :: branch

real(rp) :: s, x(3), dx(3), d(3), d1(2), d2(2), d3(2), d4(2), w1, e_angle, ds_slice
real(rp) :: gap, strength, x_center, y_center, z_center, theta_center
real(rp) :: absmax_Ez, absmax_Bz, freq, phase_lag

integer :: iu, id, n_slice, i_slice, q_sign, i_dim
integer, optional :: dimensions

character(40)   :: fieldgrid_output_name
character(24), parameter :: rfmt = 'es15.7'
character(40)  :: r_name = 'write_astra_ele'

!

i_dim = integer_option(1, dimensions)
branch => pointer_to_branch(ele)
q_sign = sign(1,  charge_of(branch%param%particle) ) 

! Get global position and rotation of the center of the element
floor1%r = [0.0_rp, 0.0_rp, ele%value(L$)/2]
floor2 = coords_local_curvilinear_to_floor (floor1, ele)
x_center = floor2%r(1)
y_center = floor2%r(2)
z_center = floor2%r(3)
theta_center = floor2%theta

write (iu, '(2a)')    '  ! Element: ', trim(ele%name)
write (iu, '(a, i0)') '  !  ix_ele: ', ele%ix_ele

select case (ele%key)

case (lcavity$, rfcavity$, e_gun$) 
  ! Example: 
  !&CAVITY
  !	 LEfield=.TRUE.
  !  C_pos(1)=0.00000000000000000
  !  File_Efield(1)='/fieldmaps/dcgun_GHV.dat'
  !	 MaxE(1)=-8.10316500000000062
  !	 Nue(1)=0 ! Frequency in Hz
  !	 Phi(1)=0.00000000000000000   ! degrees
  !	 C_smooth(1)=10
  !	 C_higher_order(1)=.TRUE.
  !  Com_grid(1) = 'all'     ! For 3D fields, if the grid is aligned
  id = id + 1

  
  if (.not. present(fieldgrid_names)) call out_io (s_error$, r_name, 'fieldgrid_names required')
  call get_astra_fieldgrid_name_and_scaling( ele, fieldgrid_names, fieldgrid_output_name, absmax_Ez, dimensions = i_dim)
  
  ! Frequency
  if ( associated(ele%grid_field)) then
    freq = ele%value(rf_frequency$) * ele%grid_field(1)%harmonic
  else
    ! No grid present
    freq = ele%value(rf_frequency$) 
  endif
  
  ! Phase 
  phase_lag = twopi*(ele%value(phi0$) +  ele%value(phi0_err$))
  ! adjust the lag for 'zero-crossing' 
  if (ele%key == rfcavity$) phase_lag = twopi*( ele%grid_field(1)%phi0_fieldmap - ele%value(phi0_max$) ) - phase_lag
  if (ele%key == e_gun$) then 
    absmax_Ez = q_sign*absmax_Ez
  endif

  write (iu, '(a, i0, 2a)') '  file_Efield(', id, ') = ', "'"//trim(fieldgrid_output_name)//"'"
  write (iu, '(a, i0, a, '//rfmt//')') '  C_pos(', id, ') = ', z_center
  if (abs(x_center) > 1e-15_rp) write (iu, '(a, i0, a, '//rfmt//')') '  C_xoff(', id, ') = ', x_center
  if (abs(y_center) > 1e-15_rp) write (iu, '(a, i0, a, '//rfmt//')') '  C_yoff(', id, ') = ', y_center
  if (abs(theta_center) > 1e-15_rp) write (iu, '(a, i0, a, '//rfmt//')') '  C_xrot(', id, ') = ', theta_center
  write (iu, '(a, i0, a, '//rfmt//', a)') '  maxE(', id, ') = ', 1e-6_rp*absmax_Ez, ' ! MV/m' ! the absolute maximum, on-axis, longitudinal electric (TM mode) in MV/m
  write (iu, '(a, i0, a, '//rfmt//', a)') '  nue(', id, ') = ', 1e-9_rp*freq, ' ! GHz'! frequency of the RF field in GHz
  write (iu, '(a, i0, a, '//rfmt//', a)') '  phi(', id, ') = ', phase_lag*180/pi, ' ! deg' 
  if (i_dim == 1) then 
    write (iu, '(a, i0, a)') '  C_smooth(', id, ') = 10'
    write (iu, '(a, i0, a)') '  C_higher_order(', id, ') = T'
  else
    write (iu, '(a, i0, a)') "  Com_grid(", id, ") = 'all'"
  endif


case (quadrupole$)
  id = id + 1
  write (iu, '(a, i0, a, '//rfmt//')') '  Q_pos(', id, ') = ', z_center
  if (abs(x_center) > 1e-15_rp) write (iu, '(a, i0, a, '//rfmt//')') '  Q_xoff(', id, ') = ', x_center
  if (abs(y_center) > 1e-15_rp) write (iu, '(a, i0, a, '//rfmt//')') '  Q_yoff(', id, ') = ', y_center
  if (abs(theta_center) > 1e-15_rp) write (iu, '(a, i0, a, '//rfmt//')') '  Q_xrot(', id, ') = ', theta_center
  write (iu, '(a, i0, a, es15.7)') '  Q_length(', id, ')  = ', ele%value(L$)
  write (iu, '(a, i0, a, es15.7, a)') '  Q_bore(', id, ')  = ', 0.005_rp, ' ! make a harder edge than the Astra default. Hard-coded at 5 mm'
  write (iu, '(a, i0, a, es15.7, a)') '  Q_grad(', id, ')    = ', q_sign* ele%value(b1_gradient$), '  ! T/m'

case (sbend$)
  
  w1 = 0.1_rp
  floor1%r = 0
  
  ! Slice if the angle is more than .2 rad (about 11 degrees)
  n_slice = floor(abs(ele%value(angle$))/.2_rp)
  ds_slice = ele%value(L$)/(n_slice+1)
  strength = q_sign*(ele%value(b_field$) + ele%value(db_field$))
  gap = 2*ele%value(HGAP$)

  if (n_slice > 0) write (iu, '(a, i0, a)') '  !  sliced into ', n_slice +1, ' parts'

  do i_slice = 0, n_slice
    id = id + 1
    
    ! D_type
    write (iu, '(a, i0, a)') '  D_type(', id, ") = 'horizontal' "
    
    ! D1, D2
    s = i_slice*ds_slice
    floor1%r = [0.0_rp, 0.0_rp, s]
    floor2 = coords_local_curvilinear_to_floor (floor1, ele)
    x = floor2%r ! Center
    floor1%r = [w1, 0.0_rp, s] ! Offset
    floor2 = coords_local_curvilinear_to_floor (floor1, ele)
    dx = floor2%r - x
    if (i_slice == 0) dx = rotate3(dx, -ele%value(e1$)) 
    D = x + dx; d1(1) = D(1); d1(2) = D(3);
    D = x - dx; d2(1) = D(1); d2(2) = D(3);
    ! D3, D4
    s = s + ds_slice
    floor1%r = [0.0_rp, 0.0_rp, s]
    floor2 = coords_local_curvilinear_to_floor (floor1, ele)
    x = floor2%r ! Center
    floor1%r = [w1, 0.0_rp, s] ! Offset
    floor2 = coords_local_curvilinear_to_floor (floor1, ele)
    dx = floor2%r - x
    if (i_slice == n_slice) dx = rotate3(dx, ele%value(e2$)) 
    D = x + dx; d3(1) = D(1); d3(2) = D(3);
    D = x - dx; d4(1) = D(1); d4(2) = D(3);
   
    call write_astra_bend(iu, strength, id, d1, d2, d3, d4)
 
    ! D_gap
    if (i_slice == 0) then
      write (iu, '(a, i0, a, '//rfmt//', a)') '  D_gap(1, ', id, ') = ', gap, ' ! Entrance edge'
    else 
      write (iu, '(a, i0, a)')  '  D_gap(1, ', id, ') = 0 ! interior'
    endif

    ! D_gap
    if (i_slice == n_slice) then
      write (iu, '(a, i0, a, '//rfmt//', a)') '  D_gap(2, ', id, ') = ', gap, ' ! Exit edge'
    else 
      write (iu, '(a, i0, a)')  '  D_gap(2, ', id, ') = 0 ! interior'
    endif
    

  enddo
  

case (solenoid$)
  ! Example
  !&SOLENOID
  !  LBField=.TRUE.'
  !  File_Bfield(1)='../fieldmaps/solenoid_SLA_L60.dat'
  !  S_pos(1)=0.303
  !  MaxB(1) = 0.054140711087  !From bmad em_field query on axis
  !  S_smooth(1)=10
  !  S_higher_order(1)=.TRUE.
  !/
  id = id + 1
  if (.not. present(fieldgrid_names)) call out_io (s_error$, r_name, 'fieldgrid_names required')
  
  ! TODO: Solenoids are 1D only. 3D maps need to use 'cavity' 
  
  call get_astra_fieldgrid_name_and_scaling( ele, fieldgrid_names, fieldgrid_output_name, absmax_bz, dimensions = 1)
  write (iu, '(a, i0, 2a)') '  file_Bfield(', id, ') = ', "'"//trim(fieldgrid_output_name)//"'"
  write (iu, '(a, i0, a, '//rfmt//')') '  S_pos(', id, ') = ', z_center
  if (abs(x_center) > 1e-15_rp) write (iu, '(a, i0, a, '//rfmt//')') '  S_xoff(', id, ') = ', x_center
  if (abs(y_center) > 1e-15_rp) write (iu, '(a, i0, a, '//rfmt//')') '  S_yoff(', id, ') = ', y_center
  if (abs(theta_center) > 1e-15_rp) write (iu, '(a, i0, a, '//rfmt//')') '  S_xrot(', id, ') = ', theta_center
  write (iu, '(a, i0, a, '//rfmt//', a)') '  maxB(', id, ') = ', absmax_bz, ' ! T'
  write (iu, '(a, i0, a)') '  s_smooth(', id, ') = 10'
  write (iu, '(a, i0, a)') '  s_higher_order(', id, ') = T'
  
case default

end select

! Footer
write (iu, '(a)') '' 

end subroutine

!--------
!
!--------
subroutine write_astra_bend(iu, strength, id, d1, d2, d3, d4)
implicit none
real(rp) :: strength, d1(2), d2(2), d3(2), d4(2)
integer :: id, iu
character(24) :: rfmt = 'es15.7'
!
write (iu, '(a, i0, a, '//rfmt//', a)') '  D_strength(', id, ') = ', strength,  '  ! T'
write (iu, '(a, i0, a, '//rfmt//', a, '//rfmt//', a)') '  D1(', id, ') = (', d1(1), ',',  d1(2), ')'
write (iu, '(a, i0, a, '//rfmt//', a, '//rfmt//', a)') '  D2(', id, ') = (', d2(1), ',',  d2(2), ')'
write (iu, '(a, i0, a, '//rfmt//', a, '//rfmt//', a)') '  D3(', id, ') = (', d3(1), ',',  d3(2), ')'
write (iu, '(a, i0, a, '//rfmt//', a, '//rfmt//', a)') '  D4(', id, ') = (', d4(1), ',',  d4(2), ')'

end subroutine


!--------
!
!--------
function rotate3(vec, angle) result (rvec)
implicit none
real(rp) :: vec(3), angle, rvec(3), ca, sa
integer :: ix1, ix2
ix1 = 3 ! Z 
ix2 = 1 ! X
ca = cos(angle)
sa = sin(angle)
rvec = vec
rvec(ix1) =  ca*vec(ix1) - sa*vec(ix2)
rvec(ix2) =  sa*vec(ix1) + ca*vec(ix2)

end function

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine  get_astra_fieldgrid_name_and_scaling(ele, name_indexx, output_name, field_scale, dimensions)
!
! Subroutine to get a field grid filename and its scaling. Calls write_astra_field_grid_file.
!   If the field grid file does not exist, it is written.
!
!   Note: This is very similar to get_opal_fieldgrid_name_and_scaling
!
!
! Input:
!   ele              -- ele_struct: element to make map
!   name_indexx      -- str_index_struct: contains field grid filenames
!   dimensions       -- integer, optional: 1 or 3 dimensions. Default: 1
!
! Output:   
!   name_indexx      -- str_index_struct: updated if new name is added
!   output_name      -- Real(rp): output filename. 
!   field_scale      -- Real(rp): the scaling of the field grid
!
!-

subroutine get_astra_fieldgrid_name_and_scaling(ele, name_indexx, output_name, field_scale, dimensions)

                                          
implicit none

type (ele_struct) :: ele
type (str_index_struct) :: name_indexx
character(*)  :: output_name
real(rp)      :: field_scale
character(200) :: unique_grid_file

integer :: ix_match, iu_fieldgrid, ios, i_dim
integer, optional :: dimensions
character(40), parameter:: r_name = 'get_astra_fieldgrid_name_and_scaling'

!

i_dim = integer_option(1, dimensions)

output_name = ''

! Check for field grid
if (associated (ele%grid_field)) then
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
  ! Call with iu=0 or '' to get field_scale
  if (i_dim ==1) then
    call write_astra_field_grid_file (0, ele, field_scale)
  else
    call write_astra_field_grid_file_3D ('', ele, field_scale)
  endif
else
  ! File does not exist.
  ! Add name to list  
  call find_index (unique_grid_file, name_indexx, ix_match, add_to_list = .true.)
  ix_match = name_indexx%n_max
  call get_output_name()
  ! Write new fieldgrid file
  if (i_dim ==1) then
    ! We need to open a file
    iu_fieldgrid = lunget()
    open(iu_fieldgrid, file = output_name, iostat = ios)
    call write_astra_field_grid_file (iu_fieldgrid , ele, field_scale)
    close(iu_fieldgrid)
  else
    ! This will write several files
    call write_astra_field_grid_file_3D (output_name, ele, field_scale)
  endif
end if

contains 

subroutine get_output_name()
implicit none
!
if (i_dim ==1) then
  write(output_name, '(a, i0, a)') '1D_fieldmap_', ix_match, '.astra'
else if (i_dim ==3) then
  write(output_name, '(a, i0)') '3D_fieldmap_', ix_match
else
  call out_io (s_error$, r_name, 'ONLY 1 and 3 dimensional field maps allowed')
  call err_exit
endif

end subroutine

end subroutine get_astra_fieldgrid_name_and_scaling

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine write_astra_field_grid_file (astra_file_unit, ele, maxfield, err)
!
!   Write 1-D field map files for Astra. The format is:
!   z field
!   ...
!  
!   Note: Simplified from write_opal_field_grid_file
!
! Input:
!   astra_file_unit -- Integer: unit number to write to, if > 0
!                        if < 0, nothing is written, and only maxfield is returned
!   ele            -- ele_struct: element to make map
!   dz             -- real(rp), optional: z step size in m. Default: 0.001 m
!
!
! Output:          
!   maxfield       -- Real(rp): absolute maximum found for element field scaling
!   err            -- Logical, optional: Set True if, say a file could not be opened.
!-



subroutine write_astra_field_grid_file (astra_file_unit, ele, maxfield, dz,  err)

implicit none


integer      :: astra_file_unit
type (ele_struct) :: ele
type (lat_param_struct) :: param
real(rp)        :: maxfield
logical, optional :: err

character(40)  :: r_name = 'write_astra_field_grid_file'
character(10)   ::  rfmt 

type (coord_struct) :: orb
type(em_field_struct) :: field_re, field_im
type (grid_field_pt1_struct), allocatable :: pt(:)
type (grid_field_pt1_struct) :: ref_field
type (branch_struct), pointer :: branch

real(rp) :: z_step, z_min, z_max
real(rp) :: freq,  z, phase_ref
real(rp) :: gap, edge_range
real(rp), optional :: dz
complex ::  phasor_rotation

integer :: nx, nz, iz, ix

real(rp) :: Ex_factor, Ez_factor, Bx_factor, By_factor, Bz_factor

logical loc_ref_frame

!
if (present(err)) err = .true.


loc_ref_frame = .true. 
branch => pointer_to_branch(ele)
param = branch%param

! Format for numbers
rfmt = 'es13.5'

! Step size
if (present(dz)) then
  z_step = dz
else
  z_step = 0.0001_rp
endif

z_min = 0.0_rp
z_max = ele%value(L$)
 
nz = ceiling(z_max/z_step)

! These are 1D fields, so the field will be probed on-axis
orb%vec(1) = 0
orb%vec(3) = 0

Ez_factor = 1
Bz_factor = 1
  
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
  ! if (freq == 0) freq = 1e-30_rp ! To prevent divide by zero

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
      print *, field_im%E(3)
    endif
    pt(iz)%E(:) = cmplx(field_re%E(:), field_im%E(:), rp)
    pt(iz)%B(:) = cmplx(field_re%B(:), field_im%B(:), rp)
    ! Update ref_field if larger Ez is found
    if(abs(pt(iz)%E(3)) > maxfield) then
      ref_field = pt(iz)
      maxfield = abs(ref_field%E(3))
    end if 
  end do
  
  ! Write to file
  if (astra_file_unit > 0 )  then
    ! Normalize
    if (maxfield > 0) Ez_factor = 1/maxfield

    ! Calculate complex rotation number to rotate Ez onto the real axis
    phase_ref = atan2( aimag(ref_field%E(3) ), real(ref_field%E(3) ) )
    phasor_rotation = cmplx(cos(phase_ref), -sin(phase_ref), rp)
    do iz = 0, nz   
      write (astra_file_unit, '(2'//rfmt//')') z_step * iz - (z_max - z_min)/2, Ez_factor * real ( pt(iz)%E(3) * phasor_rotation )
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
    if(abs(pt(iz)%B(3)) > maxfield) then
         ref_field = pt(iz)
         maxfield = abs(ref_field%B(3))
    end if 
  end do
  
  ! Restore the sign
  maxfield = ref_field%B(3)
  
  ! Write to file
  if (astra_file_unit > 0 )  then
    ! Scaling
    if (maxfield > 0) Bz_factor = 1/maxfield
    do iz = 0, nz
      write (astra_file_unit, '(2'//rfmt//')') z_step * iz - (z_max - z_min)/2, Bz_factor * real ( pt(iz)%B(3))          
    enddo
  end if

  ! cleanup 
  deallocate(pt)


  !-----------
  ! Default (gives an error)
  !-----------
  case default
  call out_io (s_error$, r_name, 'MISSING ASTRA FIELD GRID CODE FOR ELEMENT TYPE: ' // key_name(ele%key), &
             '----')
  if (global_com%exit_on_error) call err_exit
  
end select 


if (maxfield == 0) then
  call out_io (s_error$, r_name, 'ZERO MAXIMUM FIELD IN ELEMENT: ' // key_name(ele%key), &
             '----')
  if (global_com%exit_on_error) call err_exit
end if



end subroutine write_astra_field_grid_file

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine write_astra_field_grid_file_3D (base_filename, ele, maxfield, dz, err)
!
!   Writes 3-D field map files for Astra. The format is:
!   Nx x[1] x[2] ....... x[Nx-1] x[Nx] 
!   Ny y[1] y[2] ....... y[Ny-1] y[Ny] 
!   Nz z[1] z[2] ....... z[Nz-1] z[Nz]
!   <field values>
!   where field values are produced from a loop as in:
!   do iz = 1, Nz
!     do iy = 1, Ny
!       write single line: field(:, iy, iz)
!
!  
!   Note: similar to write_astra_field_grid_file
!
! Input:
!   base_filename  -- character(*): Base filename. Files will be written as:
!                                   base_filename.ex, .ey, .ez, .bx, .by, .bz
!                                   If set to '', no files will be written
!   ele            -- ele_struct: element to make map
!   dz             -- real(rp), optional: z step size in m. Default: 0.001 m
!
!
! Output:          
!   maxfield       -- Real(rp): absolute maximum on-axis field found for element field scaling
!   err            -- Logical, optional: Set True if, say a file could not be opened.
!-

subroutine write_astra_field_grid_file_3D (base_filename, ele, maxfield, dz, err)

implicit none


character(*) :: base_filename
integer      :: iu
type (ele_struct) :: ele
type (lat_param_struct) :: param

logical, optional :: err

character(40)  :: r_name = 'write_astra_field_grid_file'
character(10)   ::  rfmt 

type (coord_struct) :: orb
type(em_field_struct) :: field_re, field_im
type (grid_field_pt1_struct), allocatable :: pt(:,:,:)
type (grid_field_pt1_struct) :: ref_field
type (branch_struct), pointer :: branch

real(rp)        :: maxfield, test_field_value
real(rp) :: x_step, x_min, x_max
real(rp) :: y_step, y_min, y_max
real(rp) :: z_step, z_min, z_max
real(rp) :: freq, x, y, z, phase_ref
real(rp) :: gap, edge_range
real(rp), optional :: dz

complex ::  phasor_rotation

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

! Step size
if (present(dz)) then
  z_step = dz
else
  z_step = 0.001_rp
endif

z_min = 0.0_rp
z_max = ele%value(L$)


x_max =  .0127_rp 
x_step = z_step  ! Same as z TODO: generalize
x_min = -x_max
ix_center = nint(x_max/x_step) + 1
 
y_max =  .0127_rp
y_step = z_step  ! Same as z 
y_min = -x_max
iy_center = nint(y_max/y_step) + 1
 
nz = ceiling((z_max-z_min)/z_step) + 1
nx = ceiling((x_max-x_min)/x_step) + 1
ny = ceiling((y_max-y_min)/y_step) + 1

! Frequency

!-----------
! LCavity, RFCavity, E_GUN
!-----------
select case (ele%key)
case (lcavity$, rfcavity$, e_gun$) 
  if ( associated(ele%grid_field)) then
    freq = ele%value(rf_frequency$) * ele%grid_field(1)%harmonic
  else
    ! No grid present
    freq = ele%value(rf_frequency$) 
  endif
case default
  freq = 0
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
  ! Calculate field at \omegat*t=0 and \omega*t = \pi/2 to get real and imaginary parts
  call em_field_calc (ele, param, z, orb, loc_ref_frame, field_re, rf_time = 0.0_rp)
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
  
  
! Get reference field: The largest on-axis Ez
do iz = 1, nz
  test_field_value = astra_max_field_reference(pt(ix_center, iy_center, iz), ele)
  if(abs(test_field_value) > maxfield) then
    ref_field = pt(ix_center, iy_center, iz)
    maxfield = abs(test_field_value)
  end if 
enddo

! Calculate complex rotation number to rotate Ez onto the real axis
phase_ref = atan2( aimag(ref_field%E(3) ), real(ref_field%E(3) ) )
phasor_rotation = cmplx(cos(phase_ref), -sin(phase_ref), rp)

! Put scaling in phasor rotation
if (maxfield > 0) phasor_rotation = (1/maxfield)*phasor_rotation

! Write to file according to element

if (write_file) then
  select case(ele%key)
  case(e_gun$)
    do ifield = 1, 3
      call write_pt(ifield)
    enddo
    ! further write B-fields for nonzero frequency (RF gun)
    if (freq /= 0) then
      do ifield = 4, 6
        call write_pt(ifield)
      enddo
    endif
  case(lcavity$, rfcavity$)
    do ifield = 1, 6
      call write_pt(ifield)
    enddo
  case(solenoid$)
    do ifield = 4, 6
      call write_pt(ifield)
    enddo
  end select
endif

! Cleanup
deallocate(pt)

contains

!--- write to file: <base_filename>.ex, .ey, .ez, .bx, .by, .bz according to component
subroutine write_pt(component)
implicit none
character(40) :: output_name
integer :: iu, component, ios
!

select case(component)
  case(1)
    write(output_name, '(a)') trim(base_filename)//'.ex'
  case(2)
    write(output_name, '(a)') trim(base_filename)//'.ey'
  case(3)
    write(output_name, '(a)') trim(base_filename)//'.ez'
  case(4)
    write(output_name, '(a)') trim(base_filename)//'.bx'
  case(5)
    write(output_name, '(a)') trim(base_filename)//'.by'
  case(6)
    write(output_name, '(a)') trim(base_filename)//'.bz'
end select

iu =lunget()
open(iu, file = output_name, iostat = ios)

! Write first three lines: 
!Nx x[1] x[2] ....... x[Nx-1] x[Nx] 
!Ny y[1] y[2] ....... y[Ny-1] y[Ny] 
!Nz z[1] z[2] ....... z[Nz-1] z[Nz]

write(iu, '(i8)') nx
do ix=1, nx-1
  x = x_min + (ix-1)*x_step
  write(iu, '('//rfmt//')', advance = 'NO') x
enddo
x = x_min + (nx-1)*x_step
write(iu, '('//rfmt//')', advance = 'YES') x

write(iu, '(i8)') ny
do iy=1, ny-1
  y = y_min + (iy-1)*y_step
  write(iu, '('//rfmt//')', advance = 'NO') y
enddo
y = y_min + (ny-1)*y_step
write(iu, '('//rfmt//')', advance = 'YES') y

write(iu, '(i8)') nz
do iz=1, nz-1
  z = z_min + (iz-1)*z_step
  write(iu, '('//rfmt//')', advance = 'NO') z - ele%value(L$)/2
enddo
z = z_min + (nz-1)*z_step
write(iu, '('//rfmt//')', advance = 'YES') z - ele%value(L$)/2


if (component .le. 3) then
! E field
  do iz=1, nz
  do iy=1, ny
    do ix=1, nx-1
      write(iu, '('//rfmt//')',  advance = 'NO') real ( pt(ix, iy, iz)%E(component) * phasor_rotation)
    enddo
    write(iu, '('//rfmt//')',    advance = 'YES') real ( pt(ix, iy, nx)%E(component) * phasor_rotation)
  enddo
  enddo
  
else
  ! B field
  do iz=1, nz
  do iy=1, ny
    do ix=1, nx-1
      write(iu, '('//rfmt//')',  advance = 'NO')  real ( pt(ix, iy, iz)%B(component-3) * phasor_rotation * cmplx(0.0_rp, 1.0_rp, rp) ) ! further rotate by i
    enddo
    write(iu, '('//rfmt//')',    advance = 'YES') real ( pt(ix, iy, nz)%B(component-3) * phasor_rotation * cmplx(0.0_rp, 1.0_rp, rp) )
  enddo
  enddo
  
endif

close(iu)

end subroutine

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
!  Returns the relevant field value for the field scaling calculation
!-
function astra_max_field_reference(pt0, ele) result(field_value)
implicit none
type (grid_field_pt1_struct) pt0 
type (ele_struct) :: ele
real(rp) :: field_value
character(40), parameter :: r_name = 'astra_max_field_reference'
!
select case (ele%key)
case(e_gun$, lcavity$, rfcavity$)
  ! Ez
  field_value = abs(pt0%E(3)) ! Note: pt0 is complex
case(solenoid$)
  ! Bz
  field_value = abs(pt0%B(3))
case default
  call out_io(s_error$, r_name, 'Max field reference not specified for element type '// key_name(ele%key))

end select

end function


end module astra_interface_mod
