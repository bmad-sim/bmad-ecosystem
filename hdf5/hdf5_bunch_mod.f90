module hdf5_bunch_mod

use openpmd_interface
use bmad_interface

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine hdf5_write_beam (file_name, bunches, lat, error)

implicit none

type (bunch_struct), target :: bunches(:)
type (bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p(:), p_live
type (lat_struct) lat
type (pmd_unit_struct) unit

real(rp), allocatable :: rvec(:)
real(rp) factor

integer(HID_T) f_id, g_id, z_id, z2_id, r_id, b_id, b2_id
integer i, n, ib, ix, h5_err
integer, allocatable :: ivec(:)

character(*) file_name
character(20) date_time, root_path, bunch_path, particle_path, fmt
character(100) this_path

logical error

! Open a new file using default properties.

call h5open_f(h5_err)  ! Init Fortran interface
call h5fcreate_f (file_name, H5F_ACC_TRUNC_F, f_id, h5_err)

! Write some header stuff

call date_and_time_stamp (date_time, .true.)
root_path = '/bunch/'
bunch_path = '%T/'
particle_path = 'particles/'

call pmd_write_string_attrib(f_id, 'openPMD', '2.0.0')
call pmd_write_string_attrib(f_id, 'openPMDextension', 'BeamPhysics; SpeciesType')
call pmd_write_string_attrib(f_id, 'basePath', trim(root_path) // trim(bunch_path))
call pmd_write_string_attrib(f_id, 'particlesPath', trim(particle_path))
call pmd_write_string_attrib(f_id, 'software', 'Bmad')
call pmd_write_string_attrib(f_id, 'softwareVersion', '1.0')
call pmd_write_string_attrib(f_id, 'date', date_time)
call pmd_write_string_attrib(f_id, 'latticeFile', lat%input_file_name)
if (lat%lattice /= '') then
  call pmd_write_string_attrib(f_id, 'latticeName', lat%lattice)
endif

! Loop over bunches

call h5gcreate_f(f_id, trim(root_path), r_id, h5_err) ! Not actually needed if root_path = '/'

n = log10(1.000001 * size(bunches)) + 1
write (fmt, '(3(a, i0))') '(a, i', n, '.', n, ', 2a)'    ! EG get 'i2.2' if size(bunches) = 16

do ib = 1, size(bunches)
  bunch => bunches(ib)
  p => bunch%particle

  call re_allocate (rvec, size(p))
  call re_allocate (ivec, size(p))

  ! Write bunch particle info.
  ix = index(bunch_path, '%T')
  write (this_path, fmt) bunch_path(1:ix-1), ib, trim(bunch_path(ix+2:))
  call h5gcreate_f(r_id, trim(this_path), b_id, h5_err)
  call h5gcreate_f(b_id, particle_path, b2_id, h5_err)

  call pmd_write_string_attrib(b2_id, 'speciesType', openpmd_species_name(p(1)%species))
  call pmd_write_real_attrib(b2_id, 'totalCharge', bunch%charge_tot)
  call pmd_write_real_attrib(b2_id, 'chargeLive', bunch%charge_live)
  
  p_live => p(1)    ! In case everyone is dead.
  do i = 1, size(p)
    if (p(i)%state /= alive$) cycle
    p_live => p(i)
    exit
  enddo

  ! Position. The z-position (as opposed to z-cononical = %vec(5)) is always zero by construction.

  call h5gcreate_f(b2_id, 'position', z_id, h5_err)
  call pmd_write_real_vector_to_dataset(z_id, 'x', 'x', unit_m, p(:)%vec(1), error)
  call pmd_write_real_vector_to_dataset(z_id, 'y', 'y', unit_m, p(:)%vec(3), error)
  call pmd_write_real_to_pseudo_dataset(z_id, 'z', 'z', unit_m, 0.0_rp, size(p), error)
  call h5gclose_f(z_id, h5_err)

  ! Momentum

  factor = p_live%p0c * e_charge / c_light
  unit = pmd_unit_struct('', factor, dim_momentum)

  call h5gcreate_f(b2_id, 'momentum', z_id, h5_err)

  rvec = (p(:)%vec(2) * p(:)%p0c) / p_live%p0c
  call pmd_write_real_vector_to_dataset(z_id, 'x', 'px', unit, rvec, error)

  rvec = (p(:)%vec(4) * p(:)%p0c) / p_live%p0c
  call pmd_write_real_vector_to_dataset(z_id, 'y', 'py', unit, rvec, error)

  rvec = p(:)%direction * (sqrt((1 + p(:)%vec(6))**2 - p(:)%vec(2)**2 - p(:)%vec(4)**2) * p(:)%p0c) / p_live%p0c
  call pmd_write_real_vector_to_dataset(z_id, 'z', 'ps', unit, rvec, error)

  call h5gclose_f(z_id, h5_err)

  ! Spin.

  if (p(1)%species /= photon$) then
    call h5gcreate_f(b2_id, 'spin', z_id, h5_err)
    call pmd_write_real_vector_to_dataset(z_id, 'x', 'Sx', unit_1, p(:)%spin(1), error)
    call pmd_write_real_vector_to_dataset(z_id, 'y', 'Sy', unit_1, p(:)%spin(2), error)
    call pmd_write_real_vector_to_dataset(z_id, 'z', 'Sz', unit_1, p(:)%spin(3), error)
    call h5gclose_f(z_id, h5_err)
  endif

  ! Photon polarization

  if (p(1)%species == photon$) then
    call h5gcreate_f(b2_id, 'photonPolarization', z_id, h5_err)
    call h5gcreate_f(z_id, 'x', z2_id, h5_err)
    call pmd_write_real_vector_to_dataset(z2_id, 'amplitude', 'Field Amp_x', unit_1, p(:)%field(1), error)
    call pmd_write_real_vector_to_dataset(z2_id, 'phase', 'Field Phase_x', unit_1, p(:)%phase(1), error)
    call h5gclose_f(z2_id, h5_err)
    call h5gcreate_f(z_id, 'y', z2_id, h5_err)
    call pmd_write_real_vector_to_dataset(z2_id, 'amplitude', 'Field Amp_y', unit_1, p(:)%field(2), error)
    call pmd_write_real_vector_to_dataset(z2_id, 'phase', 'Field Phase_y', unit_1, p(:)%phase(2), error)
    call h5gclose_f(z2_id, h5_err)
    call h5gclose_f(z_id, h5_err)

    call pmd_write_real_vector_to_dataset(b2_id, 'pathLength', 'Path Length', unit_m, p(:)%field(2), error)
  endif

  !

  call pmd_write_real_vector_to_dataset(b2_id, 'sPosition', 's', unit_m, p(:)%s, error)

  rvec = -p(:)%vec(5) * p(:)%beta * c_light
  call pmd_write_real_vector_to_dataset(b2_id, 'time', 't - t_ref', unit_sec, rvec, error)
  call pmd_write_real_vector_to_dataset(b2_id, 'timeOffset', 't_ref', unit_sec, p(:)%t - rvec, error)
  call pmd_write_real_vector_to_dataset(b2_id, 'speed', 'beta', unit_c_light, p(:)%beta, error)
  call pmd_write_real_vector_to_dataset(b2_id, 'weight', 'macro-charge', unit_1, p(:)%charge, error)
  call pmd_write_real_vector_to_dataset(b2_id, 'momentumOffset', 'p0c', unit_eV, p(:)%p0c, error)

  call pmd_write_int_vector_to_dataset(b2_id, 'particleStatus', 'state', unit_1, p(:)%state, error)
  call pmd_write_int_vector_to_dataset(b2_id, 'branchIndex', 'ix_branch', unit_1, p(:)%state, error)
  call pmd_write_int_vector_to_dataset(b2_id, 'elementIndex', 'ix_ele', unit_1, p(:)%state, error)

  do i = 1, size(p)
    select case (p(i)%location)
    case (upstream_end$);   ivec(i) = -1
    case (downstream_end$); ivec(i) =  0
    case (inside$);         ivec(i) =  1
    end select
  enddo
  call pmd_write_int_vector_to_dataset(b2_id, 'locationInElement', 'location', unit_1, p(:)%location, error)

  call h5gclose_f(b2_id, h5_err)
  call h5gclose_f(b_id, h5_err)
enddo

call h5gclose_f(r_id, h5_err)
call h5fclose_f(f_id, h5_err)

end subroutine hdf5_write_beam 

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine hdf5_read_beam (file_name, beam, pmd_header, error)

use iso_c_binding

implicit none

type (beam_struct), target :: beam
type (pmd_header_struct) pmd_header
type (lat_struct) lat
type (pmd_unit_struct) unit
type(c_ptr) cv_ptr
type(c_funptr) c_func_ptr

real(rp), allocatable :: rvec(:)
real(rp) factor

integer(HID_T) f_id, g_id, z_id, z2_id, r_id, b_id
integer(HSIZE_T) idx
integer i, ib, ix, is, it, state, h5_err, n_bunch
integer, allocatable :: ivec(:)

character(*) file_name
character(*), parameter :: r_name = 'hdf5_read_beam'
logical error, err

! Init

error = .true.

call h5open_f(h5_err)  ! Init Fortran interface
call h5fopen_f(file_name, H5F_ACC_RDONLY_F, f_id, h5_err)

! Get header info

call pmd_read_attribute_string (f_id, '/', 'openPMD', pmd_header%openPMD, err, .true.);                    if (err) return
call pmd_read_attribute_string (f_id, '/', 'openPMDextension', pmd_header%openPMDextension, err, .true.);  if (err) return
call pmd_read_attribute_string (f_id, '/', 'basePath', pmd_header%basePath, err, .true.);                  if (err) return
call pmd_read_attribute_string (f_id, '/', 'particlesPath', pmd_header%particlesPath, err, .true.);        if (err) return
call pmd_read_attribute_string (f_id, '/', 'software', pmd_header%software, err, .false.)
call pmd_read_attribute_string (f_id, '/', 'softwareVersion', pmd_header%softwareVersion, err, .false.)
call pmd_read_attribute_string (f_id, '/', 'date', pmd_header%date, err, .false.)
call pmd_read_attribute_string (f_id, '/', 'latticeFile', pmd_header%latticeFile, err, .false.)
call pmd_read_attribute_string (f_id, '/', 'latticeName', pmd_header%latticeName, err, .false.)

! Find root group

it = index(pmd_header%basePath, '/%T/')
if (it /= len_trim(pmd_header%basePath) - 3) then
  call out_io (s_error$, r_name, 'PARSING OF BASE PATH FAILURE. PLEASE REPORT THIS. ' // pmd_header%basePath)
  return
endif

call pmd_find_group(f_id, pmd_header%basePath(1:it), z_id, err, .true.)

! Loop over all bunches

n_bunch = 0
idx = 0
call H5Literate_f (z_id, H5_INDEX_NAME_F, H5_ITER_INC_F, idx, c_funloc(hdf5_count_bunches), C_NULL_PTR, state, h5_err)

call reallocate_beam(beam, n_bunch)

cv_ptr = c_loc(beam)
c_func_ptr = c_funloc(hdf5_read_bunch)
idx = 0
n_bunch = 0
call H5Literate_f (z_id, H5_INDEX_NAME_F, H5_ITER_INC_F, idx, c_func_ptr, cv_ptr, state, h5_err)

!

call h5fclose_f(f_id, h5_err)
error = .false.
return

!------------------------------------------------------------
! Close and return

9000 continue

call h5fclose_f(f_id, h5_err)
return

!------------------------------------------------------------------------------------------
contains

function hdf5_count_bunches (g_id, g_name, info, c_point) result (stat) bind(C)

use iso_c_binding
use fortran_cpp_utils

implicit none

type(c_ptr) info
type(c_ptr) c_point
type(H5O_info_t) :: infobuf 

integer(hid_t), value :: g_id
integer stat, h5_stat

character(1) :: g_name(*)
character(100) group_name

!

stat = 0

call to_f_str(g_name, group_name)
call H5Oget_info_by_name_f(g_id, group_name, infobuf, h5_stat)

if (infobuf%type /= H5O_TYPE_GROUP_F) return
n_bunch = n_bunch + 1

end function  hdf5_count_bunches

!------------------------------------------------------------------------------------------
! contains

function hdf5_read_bunch (root_id, g_name_c, info, dummy_c_ptr) result (stat) bind(C)

use iso_c_binding
use fortran_cpp_utils

implicit none

type(c_ptr) info
type(c_ptr) dummy_c_ptr
type(H5O_info_t) :: infobuf 
type(bunch_struct), pointer :: bunch

integer(hid_t), value :: root_id
integer(hid_t) g_id, g2_id
integer stat, h5_stat

logical error

character(1) :: g_name_c(*)
character(100) g_name

! Return if not a group

stat = 0

call to_f_str(g_name_c, g_name)
call H5Oget_info_by_name_f(root_id, g_name, infobuf, h5_stat)

if (infobuf%type /= H5O_TYPE_GROUP_F) return

!

print *, trim(g_name)
n_bunch = n_bunch + 1
bunch => beam%bunch(n_bunch)

call pmd_find_group(root_id, g_name, g_id, error, .true.)
call pmd_find_group(g_id, pmd_header%particlesPath, g2_id, error, .true.);  if (error) return

! Get attributes.

! Loop over all sub-groups/datasets.



end function  hdf5_read_bunch

end subroutine hdf5_read_beam

end module hdf5_bunch_mod
