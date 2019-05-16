module hdf5_bunch_mod

use hdf5_interface
use bmad_interface
use iso_fortran_env

! Common units

type pmd_unit_struct
 character(8) :: unitSymbol = ''     ! Native units name. EG 'eV'
 real(rp) :: unitSI = 0              ! Conversion to SI
 real(rp) :: unitDimension(7) = 0    ! SI Base Exponents
end type

real(rp), parameter :: dim_1(7)              = [0,0,0,0,0,0,0]
real(rp), parameter :: dim_length(7)         = [1,0,0,0,0,0,0]
real(rp), parameter :: dim_mass(7)           = [0,1,0,0,0,0,0]
real(rp), parameter :: dim_time(7)           = [0,0,1,0,0,0,0]
real(rp), parameter :: dim_current(7)        = [0,0,0,1,0,0,0]
real(rp), parameter :: dim_temperture(7)     = [0,0,0,0,1,0,0]
real(rp), parameter :: dim_mol(7)            = [0,0,0,0,0,1,0]
real(rp), parameter :: dim_luminous(7)       = [0,0,0,0,0,0,1]
real(rp), parameter :: dim_charge(7)         = [0,0,1,1,0,0,0]
real(rp), parameter :: dim_elec_potential(7) = [1,1,-3,-1,0,0,0]
real(rp), parameter :: dim_velocity(7)       = [1,0,-1,0,0,0,0]
real(rp), parameter :: dim_energy(7)         = [2,1,-2,0,0,0,0]
real(rp), parameter :: dim_momentum(7)       = [1,1,-1,0,0,0,0]

type(pmd_unit_struct), parameter :: unit_1          = pmd_unit_struct('', 1.0_rp, dim_1)
type(pmd_unit_struct), parameter :: unit_m          = pmd_unit_struct('m', 1.0_rp, dim_length)
type(pmd_unit_struct), parameter :: unit_kg         = pmd_unit_struct('kg', 1.0_rp, dim_mass)
type(pmd_unit_struct), parameter :: unit_sec        = pmd_unit_struct('sec', 1.0_rp, dim_time)
type(pmd_unit_struct), parameter :: unit_amp        = pmd_unit_struct('Amp', 1.0_rp, dim_current)
type(pmd_unit_struct), parameter :: unit_K          = pmd_unit_struct('K', 1.0_rp, dim_temperture)
type(pmd_unit_struct), parameter :: unit_mol        = pmd_unit_struct('mol', 1.0_rp, dim_mol)
type(pmd_unit_struct), parameter :: unit_cd         = pmd_unit_struct('cd', 1.0_rp, dim_luminous)
type(pmd_unit_struct), parameter :: unit_Coulomb    = pmd_unit_struct('Coulomb', 1.0_rp, dim_charge)
type(pmd_unit_struct), parameter :: unit_charge_num = pmd_unit_struct('charge #', 1.0_rp, dim_charge)
type(pmd_unit_struct), parameter :: unit_V_per_m    = pmd_unit_struct('V/m', 1.0_rp, dim_elec_potential)
type(pmd_unit_struct), parameter :: unit_c_light    = pmd_unit_struct('vel/c', c_light, dim_velocity)
type(pmd_unit_struct), parameter :: unit_m_per_s    = pmd_unit_struct('m/s', 1.0_rp, dim_velocity)
type(pmd_unit_struct), parameter :: unit_eV         = pmd_unit_struct('eV', e_charge, dim_energy)
type(pmd_unit_struct), parameter :: unit_eV_per_c   = pmd_unit_struct('eV/c', e_charge/c_light, dim_momentum)

! Header information

type pmd_header_struct
  character(:), allocatable :: openPMD
  character(:), allocatable :: openPMDextension
  character(:), allocatable :: basePath
  character(:), allocatable :: particlesPath
  character(:), allocatable :: author          != 'anonymous'
  character(:), allocatable :: software        != 'Bmad'
  character(:), allocatable :: softwareVersion ! = 'Revision XXX'
  character(:), allocatable :: date            ! = ''
  character(:), allocatable :: latticeFile
  character(:), allocatable :: latticeName
end type

!

private pmd_write_real_vector_to_dataset, pmd_write_real_to_pseudo_dataset
private pmd_write_int_vector_to_dataset, pmd_write_int_to_pseudo_dataset
private pmd_write_units_to_dataset 
private pmd_read_int_dataset, pmd_read_real_dataset

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_write_beam (file_name, bunches, append, error)
!
! Routine to write particle positions of a beam to an HDF5 binary file.
! See also hdf5_read_beam.
!
! Input:
!   file_name     -- character(*): Name of the file to create.
!   bunches(:)    -- bunch_struct: Array of bunches. 
!                       Use "[bunch]" if you have a single bunch.
!                       use "beam%bunch" if you have beam_struct instance.
!   append        -- logical: If True then append if the file already exists.
!
! Output:
!   error         -- logical: Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_beam (file_name, bunches, append, error)

implicit none

type (bunch_struct), target :: bunches(:)
type (bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p(:), p_live

real(rp), allocatable :: rvec(:)

integer(HID_T) f_id, g_id, z_id, z2_id, r_id, b_id, b2_id
integer i, n, ib, ix, h5_err, h_err
integer, allocatable :: ivec(:)

character(*) file_name
character(20) date_time, root_path, bunch_path, particle_path, fmt
character(100) this_path
character(*), parameter :: r_name = 'hdf5_write_beam'

logical append, error, err

! Open a new file using default properties.

call hdf5_open_file (file_name, 'WRITE', f_id, err);  if (err) return

! Write some header stuff

call date_and_time_stamp (date_time, .true., .true.)
root_path = '/data/'
bunch_path = '%T/'
particle_path = 'particles/'

call hdf5_write_attribute_string(f_id, 'fileType', 'openPMD', err)
call hdf5_write_attribute_string(f_id, 'openPMD', '2.0.0', err)
call hdf5_write_attribute_string(f_id, 'openPMDextension', 'BeamPhysics;SpeciesType', err)
call hdf5_write_attribute_string(f_id, 'basePath', trim(root_path) // trim(bunch_path), err)
call hdf5_write_attribute_string(f_id, 'particlesPath', trim(particle_path), err)
call hdf5_write_attribute_string(f_id, 'software', 'Bmad', err)
call hdf5_write_attribute_string(f_id, 'softwareVersion', '1.0', err)
call hdf5_write_attribute_string(f_id, 'date', date_time, err)
!! call hdf5_write_attribute_string(f_id, 'latticeFile', lat%input_file_name, err)
!! if (lat%lattice /= '') then
!!   call hdf5_write_attribute_string(f_id, 'latticeName', lat%lattice, err)
!! endif

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

  call hdf5_write_attribute_string(b2_id, 'speciesType', openpmd_species_name(p(1)%species), err)
  call hdf5_write_attribute_real(b2_id, 'totalCharge', bunch%charge_tot, err)
  call hdf5_write_attribute_real(b2_id, 'chargeLive', bunch%charge_live, err)
  call hdf5_write_attribute_real(b2_id, 'chargeUnitSI', 1.0_rp, err)
  call hdf5_write_attribute_int(b2_id, 'numParticles', size(bunch%particle), err)
  
  p_live => p(1)    ! In case everyone is dead.
  do i = 1, size(p)
    if (p(i)%state /= alive$) cycle
    p_live => p(i)
    exit
  enddo

  !-----------------
  ! Photons...

  if (p(1)%species == photon$) then
    ! Photon polarization
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

    ! Spin.

    call h5gcreate_f(b2_id, 'spin', z_id, h5_err)
    call pmd_write_real_vector_to_dataset(z_id, 'x', 'Sx', unit_1, p(:)%spin(1), error)
    call pmd_write_real_vector_to_dataset(z_id, 'y', 'Sy', unit_1, p(:)%spin(2), error)
    call pmd_write_real_vector_to_dataset(z_id, 'z', 'Sz', unit_1, p(:)%spin(3), error)
    call h5gclose_f(z_id, h5_err)

    ! Velocity

    call h5gcreate_f(b2_id, 'velocity', z_id, h5_err)

    call pmd_write_real_vector_to_dataset(z_id, 'x', 'Vx', unit_m_per_s, p(:)%vec(2), error)
    call pmd_write_real_vector_to_dataset(z_id, 'y', 'Vy', unit_m_per_s, p(:)%vec(4), error)
    call pmd_write_real_vector_to_dataset(z_id, 'z', 'Vz', unit_m_per_s, p(:)%vec(6), error)

    call h5gclose_f(z_id, h5_err)

    call pmd_write_real_vector_to_dataset (b2_id, 'totalMomentum', 'p0c', unit_m_per_s, p(:)%p0c, error)

  ! Non-photons...

  else
    ! Momentum

    call h5gcreate_f(b2_id, 'momentum', z_id, h5_err)

    call pmd_write_real_vector_to_dataset(z_id, 'x', 'px * p0c', unit_ev_per_c, p(:)%vec(2) * p(:)%p0c, error)
    call pmd_write_real_vector_to_dataset(z_id, 'y', 'py * p0c', unit_ev_per_c, p(:)%vec(4) * p(:)%p0c, error)
    rvec = p(:)%direction * (sqrt((1 + p(:)%vec(6))**2 - p(:)%vec(2)**2 - p(:)%vec(4)**2) * p(:)%p0c)
    call pmd_write_real_vector_to_dataset(z_id, 'z', 'ps * p0c', unit_ev_per_c, rvec, error)

    call h5gclose_f(z_id, h5_err)

    call pmd_write_real_vector_to_dataset (b2_id, 'totalMomentumOffset', 'p0c', unit_eV_per_c, p(:)%p0c, error)
    call pmd_write_real_vector_to_dataset (b2_id, 'totalMomentum', 'pz * p0c', unit_eV_per_c, p(:)%vec(6)*p(:)%p0c, error)
  endif

  !-----------------
  ! Position. 

  call h5gcreate_f(b2_id, 'position', z_id, h5_err)

  call pmd_write_real_vector_to_dataset(z_id, 'x', 'x', unit_m, p(:)%vec(1), error)
  call pmd_write_real_vector_to_dataset(z_id, 'y', 'y', unit_m, p(:)%vec(3), error)

  ! For charged particles: The z-position (as opposed to z-cononical = %vec(5)) is always zero by construction.
  if (p(1)%species == photon$) then
    call pmd_write_real_to_pseudo_dataset(z_id, 'z', 'z', unit_m, 0.0_rp, size(p), error)
  else
  call pmd_write_real_vector_to_dataset(z_id, 'z', 'z', unit_m, p(:)%vec(5), error)
  endif

  call h5gclose_f(z_id, h5_err)

  !

  call pmd_write_real_vector_to_dataset(b2_id, 'sPosition', 's', unit_m, p(:)%s, error)

  rvec = -p(:)%vec(5) / (p(:)%beta * c_light)  ! t - t_ref
  call pmd_write_real_vector_to_dataset(b2_id, 'time', 't - t_ref', unit_sec, rvec, error)
  call pmd_write_real_vector_to_dataset(b2_id, 'timeOffset', 't_ref', unit_sec, p(:)%t - rvec, error)
  call pmd_write_real_vector_to_dataset(b2_id, 'speed', 'beta', unit_c_light, p(:)%beta, error)
  call pmd_write_real_vector_to_dataset(b2_id, 'weight', 'macro-charge', unit_1, p(:)%charge, error)

  call pmd_write_int_vector_to_dataset(b2_id, 'particleStatus', 'state', unit_1, p(:)%state, error)
  call pmd_write_int_vector_to_dataset(b2_id, 'branchIndex', 'ix_branch', unit_1, p(:)%ix_branch, error)
  call pmd_write_int_vector_to_dataset(b2_id, 'elementIndex', 'ix_ele', unit_1, p(:)%ix_ele, error)

  if (.not. is_fundamental_species(p(1)%species)) then
    do i = 1, size(p)
      ivec = charge_of(p(i)%species)
    enddo
    call pmd_write_int_vector_to_dataset(b2_id, 'chargeState', 'particle charge', unit_1, ivec, error)
  endif

  do i = 1, size(p)
    select case (p(i)%location)
    case (upstream_end$);   ivec(i) = -1
    case (inside$);         ivec(i) =  0
    case (downstream_end$); ivec(i) =  1
    end select
  enddo
  call pmd_write_int_vector_to_dataset(b2_id, 'locationInElement', 'location', unit_1, ivec, error)

  call h5gclose_f(b2_id, h5_err)
  call h5gclose_f(b_id, h5_err)
enddo

call h5gclose_f(r_id, h5_err)
call h5fclose_f(f_id, h5_err)

end subroutine hdf5_write_beam 

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_real_vector_to_dataset(root_id, dataset_name, bmad_name, unit, vector, error)
!
! Private routine used by hdf5_write_beam.
!-

subroutine pmd_write_real_vector_to_dataset(root_id, dataset_name, bmad_name, unit, vector, error)

type (pmd_unit_struct) unit
real(rp) vector(:), v_max, v_min
integer err
integer(HID_T) :: root_id, v_size
character(*) dataset_name, bmad_name
logical error

!

v_max = maxval(vector)
v_min = minval(vector)

if (v_max == v_min) then
  call pmd_write_real_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, v_max, size(vector), error) 
  return
endif

!

v_size = size(vector)
call H5LTmake_dataset_double_f(root_id, dataset_name, 1, [v_size], vector, err)    

call H5LTset_attribute_double_f(root_id, dataset_name, 'minValue', [v_min], 1_hz, err)
call H5LTset_attribute_double_f(root_id, dataset_name, 'maxValue', [v_max], 1_hz, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_real_vector_to_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_real_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, value, v_size, error)
!
! Private routine used by hdf5_write_beam.
!-

subroutine pmd_write_real_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, value, v_size, error)

type (pmd_unit_struct) unit
integer(HID_T) :: root_id, group_id
real(rp) value
integer v_size, err
character(*) dataset_name, bmad_name
logical error

!

call h5gcreate_f(root_id, dataset_name, group_id, err)

call H5LTset_attribute_double_f(root_id, dataset_name, 'value', [value], 1_hz, err)
call H5LTset_attribute_int_f(root_id, dataset_name, 'shape', [v_size], 1_hz, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

call h5gclose_f(group_id, err)

end subroutine pmd_write_real_to_pseudo_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_int_vector_to_dataset(root_id, dataset_name, bmad_name, unit, vector, error)
!
! Private routine used by hdf5_write_beam.
!-

subroutine pmd_write_int_vector_to_dataset(root_id, dataset_name, bmad_name, unit, vector, error)

type (pmd_unit_struct) unit
integer vector(:), v_max, v_min
integer err
integer(HID_T) :: root_id, v_size
character(*) dataset_name, bmad_name
logical error

!

v_max = maxval(vector)
v_min = minval(vector)

if (v_max == v_min) then
  call pmd_write_int_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, v_max, size(vector), error) 
  return
endif

!

v_size = size(vector)
call H5LTmake_dataset_int_f(root_id, dataset_name, 1, [v_size], vector, err)    

call H5LTset_attribute_int_f(root_id, dataset_name, 'minValue', [v_min], 1_hz, err)
call H5LTset_attribute_int_f(root_id, dataset_name, 'maxValue', [v_max], 1_hz, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_int_vector_to_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_int_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, value, v_size, error)
!
! Private routine used by hdf5_write_beam.
!-

subroutine pmd_write_int_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, value, v_size, error)

type (pmd_unit_struct) unit
integer(HID_T) :: root_id, group_id
integer value
integer err, v_size
character(*) dataset_name, bmad_name
logical error

!

call h5gcreate_f(root_id, dataset_name, group_id, err)

call H5LTset_attribute_int_f(root_id, dataset_name, 'value', [value], 1_hz, err)
call H5LTset_attribute_int_f(root_id, dataset_name, 'shape', [v_size], 1_hz, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

call h5gclose_f(group_id, err)

end subroutine pmd_write_int_to_pseudo_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)
!
! Private routine used by hdf5_write_beam.
!-

subroutine pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

type (pmd_unit_struct) unit
integer(HID_T) :: root_id
integer err
logical error
character(*) dataset_name, bmad_name

if (bmad_name /= '') call H5LTset_attribute_string_f(root_id, dataset_name, 'localName', bmad_name, err)
call H5LTset_attribute_double_f(root_id, dataset_name, 'unitSI', [unit%unitSI], 1_hz, err)
call H5LTset_attribute_double_f(root_id, dataset_name, 'unitDimension', [unit%unitDimension], 7_hz, err)
call H5LTset_attribute_string_f(root_id, dataset_name, 'unitSymbol', unit%unitSymbol, err)

end subroutine pmd_write_units_to_dataset 

!===============================================================================================
!===============================================================================================
!===============================================================================================
!+
! Subroutine hdf5_read_beam (file_name, beam, pmd_header, error)
!
! Routine to read a beam data file. 
! See also hdf5_write_beam
!
! Input:
!   file_name         -- character(*): Name of the beam data file.
!
! Output:
!   beam              -- beam_struct: Particle positions.
!   pmd_header        -- pmd_header_struct: Extra info like file creation date.
!   error             -- logical: Set True if there is a read error. False otherwise.
!-

subroutine hdf5_read_beam (file_name, beam, pmd_header, error)

use iso_c_binding
use fortran_cpp_utils

implicit none

type (beam_struct), target :: beam
type (pmd_header_struct) pmd_header
type (pmd_unit_struct) unit
type(H5O_info_t) :: infobuf 
type(c_ptr) cv_ptr
type(c_funptr) c_func_ptr

real(rp), allocatable :: rvec(:)
real(rp) factor

integer(HID_T) f_id, g_id, z_id, z2_id, r_id, b_id
integer(HSIZE_T) idx
integer(SIZE_T) g_size
integer i, ik, ib, ix, is, it, state, h5_err, n_bunch, storage_type, n_links, max_corder, h5_stat
integer h_err
integer, allocatable :: ivec(:)

character(*) file_name
character(*), parameter :: r_name = 'hdf5_read_beam'
character(100) c_name, name
logical error, err

! Init

error = .true.

call hdf5_open_file (file_name, 'READ', f_id, err);  if (err) return

! Get header info

call hdf5_read_attribute_alloc_string (f_id, 'openPMD', pmd_header%openPMD, err, .true.);                    if (err) return
call hdf5_read_attribute_alloc_string (f_id, 'openPMDextension', pmd_header%openPMDextension, err, .true.);  if (err) return
call hdf5_read_attribute_alloc_string (f_id, 'basePath', pmd_header%basePath, err, .true.);                  if (err) return
call hdf5_read_attribute_alloc_string (f_id, 'particlesPath', pmd_header%particlesPath, err, .true.);        if (err) return
call hdf5_read_attribute_alloc_string (f_id, 'software', pmd_header%software, err, .false.)
call hdf5_read_attribute_alloc_string (f_id, 'softwareVersion', pmd_header%softwareVersion, err, .false.)
call hdf5_read_attribute_alloc_string (f_id, 'date', pmd_header%date, err, .false.)
call hdf5_read_attribute_alloc_string (f_id, 'latticeFile', pmd_header%latticeFile, err, .false.)
call hdf5_read_attribute_alloc_string (f_id, 'latticeName', pmd_header%latticeName, err, .false.)

! Find root group
! Note: Right now it is assumed that the basepath uses "/%T/" to sort bunches.

it = index(pmd_header%basePath, '/%T/')
if (it /= len_trim(pmd_header%basePath) - 3) then
  call out_io (s_error$, r_name, 'PARSING OF BASE PATH FAILURE. PLEASE REPORT THIS. ' // pmd_header%basePath)
  return
endif

z_id = hdf5_open_group(f_id, pmd_header%basePath(1:it), err, .true.)

! Count bunches

n_bunch = 0
call H5Gget_info_f (z_id, storage_type, n_links, max_corder, h5_err)
do idx = 0, n_links-1
  call H5Lget_name_by_idx_f (z_id, '.', H5_INDEX_NAME_F, H5_ITER_INC_F, idx, c_name, h5_err, g_size)
  call to_f_str(c_name, name)
  call H5Oget_info_by_name_f(z_id, name, infobuf, h5_stat)
  if (infobuf%type /= H5O_TYPE_GROUP_F) cycle    ! Ignore non-group elements.
  if (.not. is_integer(name)) then    ! This assumes basepath uses "/%T/" to designate different bunches.
    call out_io (s_warn$, r_name, 'NAME OF DIRECTORY IN PATH IS NOT AN INTEGER: ' // quote(name))
    cycle
  endif
  n_bunch = n_bunch + 1
enddo

call reallocate_beam(beam, n_bunch)

! Loop over all bunches.
! Note: There is a GCC compiler bug where, if building shared, there is an error if 
! "c_funloc(pmd_read_bunch)" is used as the actual arg to H5Literate_f in place of c_fun_ptr.

cv_ptr = c_loc(beam)
c_func_ptr = c_funloc(pmd_read_bunch)
idx = 0
n_bunch = 0
call H5Literate_f (z_id, H5_INDEX_NAME_F, H5_ITER_INC_F, idx, c_func_ptr, cv_ptr, state, h5_err)

! And end

call H5Gclose_f(z_id, h5_err)
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

function pmd_read_bunch (root_id, g_name_c, info, dummy_c_ptr) result (stat) bind(C)

type(c_ptr) info
type(c_ptr) dummy_c_ptr
type(H5O_info_t) :: infobuf 
type(bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p
type(c_funptr) c_func_ptr

real(rp) charge_factor, f_ev
real(rp), allocatable :: dt(:)

integer(hid_t), value :: root_id
integer(hid_t) g_id, g2_id, g3_id, g4_id, a_id
integer(HSIZE_T) idx
integer n, stat, h5_stat, ia, ip, species
integer, allocatable :: charge_state(:)

logical error

character(1) :: g_name_c(*)
character(:), allocatable :: string
character(100) g_name, a_name, name, c_name
character(*), parameter :: r_name = 'pmd_read_bunch'

! Return if not a group with the proper name

stat = 0

call to_f_str(g_name_c, g_name)
call H5Oget_info_by_name_f(root_id, g_name, infobuf, h5_stat)

if (infobuf%type /= H5O_TYPE_GROUP_F) return
if (.not. is_integer(g_name)) return     ! This assumes basepath uses "/%T/" to designate different bunches.

!

n_bunch = n_bunch + 1
bunch => beam%bunch(n_bunch)
bunch%charge_tot = 0
bunch%charge_live = 0

g_id = hdf5_open_group(root_id, g_name, error, .true.);  if (error) return
g2_id = hdf5_open_group(g_id, pmd_header%particlesPath, error, .true.);  if (error) return

! Get number of particles.

call hdf5_read_attribute_int(g2_id, 'numParticles', n, error, .true.)
if (error) return
call reallocate_bunch(bunch, n)
bunch%particle = coord_struct()
allocate (dt(n), charge_state(n))

! Get attributes.

charge_factor = 0
species = int_garbage$  ! Garbage number
charge_state = 0

do ia = 1, hdf5_num_attributes (g2_id)
  call hdf5_get_attribute_by_index(g2_id, ia, a_id, a_name)
  select case (a_name)
  case ('speciesType')
    call hdf5_read_attribute_alloc_string(g2_id, a_name, string, error, .true.);  if (error) return
    species = species_id_from_openpmd(string, 0)
  case ('numParticles')
    ! Nothing to be done
  case ('totalCharge')
    call hdf5_read_attribute_real(g2_id, a_name, bunch%charge_tot, error, .true.);  if (error) return
  case ('chargeLive')
    call hdf5_read_attribute_real(g2_id, a_name, bunch%charge_live, error, .true.);  if (error) return
  case ('chargeUnitSI')
    call hdf5_read_attribute_real(g2_id, a_name, charge_factor, error, .true.);  if (error) return
  case default
    call out_io (s_warn$, r_name, 'Unknown bunch attribute: ', quote(a_name))
  end select
enddo

if (charge_factor == 0 .and. (bunch%charge_live /= 0 .or. bunch%charge_tot /= 0)) then
  call out_io (s_error$, r_name, 'chargeUnitSI is not set for bunch')
else
  bunch%charge_live = bunch%charge_live * charge_factor
  bunch%charge_tot = bunch%charge_tot * charge_factor
endif

! Loop over all datasets.

f_ev = e_charge / c_light

call H5Gget_info_f (g2_id, storage_type, n_links, max_corder, h5_err)
do idx = 0, n_links-1
  call H5Lget_name_by_idx_f (g2_id, '.', H5_INDEX_NAME_F, H5_ITER_INC_F, idx, name, h5_err, g_size)
  call H5Oget_info_by_name_f(g2_id, name, infobuf, h5_stat)
  select case (name)
  case ('spin')
    call pmd_read_real_dataset(g2_id, 'spin/x', 1.0_rp, bunch%particle%spin(1), error)
    call pmd_read_real_dataset(g2_id, 'spin/y', 1.0_rp, bunch%particle%spin(2), error)
    call pmd_read_real_dataset(g2_id, 'spin/z', 1.0_rp, bunch%particle%spin(3), error)
  case ('position')
    call pmd_read_real_dataset(g2_id, 'position/x', 1.0_rp, bunch%particle%vec(1), error)
    call pmd_read_real_dataset(g2_id, 'position/y', 1.0_rp, bunch%particle%vec(3), error)
    call pmd_read_real_dataset(g2_id, 'position/z', 1.0_rp, bunch%particle%vec(5), error)
  case ('momentum')
    call pmd_read_real_dataset(g2_id, 'momentum/x', f_ev, bunch%particle%vec(2), error)
    call pmd_read_real_dataset(g2_id, 'momentum/y', f_ev, bunch%particle%vec(4), error)
  case ('velocity')
    call pmd_read_real_dataset(g2_id, 'velocity/x', f_ev, bunch%particle%vec(2), error)
    call pmd_read_real_dataset(g2_id, 'velocity/y', f_ev, bunch%particle%vec(4), error)
    call pmd_read_real_dataset(g2_id, 'velocity/z', f_ev, bunch%particle%vec(6), error)
  case ('pathLength')
    call pmd_read_real_dataset(g2_id, name, 1.0_rp, bunch%particle%path_len, error)
  case ('totalMomentumOffset')
    call pmd_read_real_dataset(g2_id, name, f_ev, bunch%particle%p0c, error)
  case ('totalMomentum')
    call pmd_read_real_dataset(g2_id, name, f_ev, bunch%particle%vec(6), error)
  case ('photonPolarization')
    g3_id = hdf5_open_group(g2_id, 'photonPolarization', error, .true.)
    g4_id = hdf5_open_group(g3_id, 'x/amplitude', error, .false.)
    if (error) then
      call pmd_read_real_dataset(g3_id, 'x/real', 1.0_rp, bunch%particle%field(1), error)
      call pmd_read_real_dataset(g3_id, 'x/imaginary', 1.0_rp, bunch%particle%phase(1), error)
      call to_amp_phase(bunch%particle%field(1), bunch%particle%phase(1))
      call pmd_read_real_dataset(g3_id, 'y/real', 1.0_rp, bunch%particle%field(2), error)
      call pmd_read_real_dataset(g3_id, 'y/imaginary', 1.0_rp, bunch%particle%phase(2), error)
      call to_amp_phase(bunch%particle%field(2), bunch%particle%phase(2))
    else
      call pmd_read_real_dataset(g3_id, 'x/amplitude', 1.0_rp, bunch%particle%field(1), error)
      call pmd_read_real_dataset(g3_id, 'x/phase', 1.0_rp, bunch%particle%phase(1), error)
      call pmd_read_real_dataset(g3_id, 'x/amplitude', 1.0_rp, bunch%particle%field(2), error)
      call pmd_read_real_dataset(g3_id, 'x/phase', 1.0_rp, bunch%particle%phase(2), error)
    endif
    call H5Gclose_f(g4_id, h5_err)
    call H5Gclose_f(g3_id, h5_err)
  case ('sPosition')
    call pmd_read_real_dataset(g2_id, name, 1.0_rp, bunch%particle%s, error)
  case ('time')
    call pmd_read_real_dataset(g2_id, name, 1.0_rp, dt, error)
  case ('timeOffset')
    call pmd_read_real_dataset(g2_id, name, 1.0_rp, bunch%particle%t, error)
  case ('speed')
    call pmd_read_real_dataset(g2_id, name, c_light, bunch%particle%beta, error)
  case ('weight')
    call pmd_read_real_dataset(g2_id, name, 1.0_rp, bunch%particle%charge, error)
  case ('particleStatus')
    call pmd_read_int_dataset(g2_id, name, bunch%particle%state, error)
  case ('chargeState')
    call pmd_read_int_dataset(g2_id, name, charge_state, error)
  case ('branchIndex')
    call pmd_read_int_dataset(g2_id, name, bunch%particle%ix_branch, error)
  case ('elementIndex')
    call pmd_read_int_dataset(g2_id, name, bunch%particle%ix_ele, error)
  case ('locationInElement')
    call pmd_read_int_dataset(g2_id, name, bunch%particle%location, error)
    do ip = 1, size(bunch%particle)
      p => bunch%particle(ip)
      select case(p%location)
      case (-1);    p%location = upstream_end$
      case ( 0);    p%location = inside$
      case ( 1);    p%location = downstream_end$
      end select
    enddo
  end select
enddo

call H5Gclose_f(g2_id, h5_err)
call H5Gclose_f(g_id, h5_err)

do ip = 1, size(bunch%particle)
  p => bunch%particle(ip)
  if (species == photon$) then
  else
    p%vec(5) = -p%beta * c_light * dt(ip)
    p%t = p%t + dt(ip)
    p%vec(2) = p%vec(2) / p%p0c
    p%vec(4) = p%vec(4) / p%p0c
    p%vec(6) = p%vec(6) / p%p0c
    p%species = set_species_charge(species, charge_state(ip))
  endif
enddo

end function  pmd_read_bunch

!------------------------------------------------------------------------------------------
! contains

subroutine to_amp_phase(re_amp, im_phase)

real(rp) re_amp(:), im_phase(:)
real(rp) amp
integer i

!

do i = 1, size(re_amp)
  amp = norm2([re_amp(i), im_phase(i)])
  im_phase(i) = atan2(im_phase(i), re_amp(i))
  re_amp(i) = amp
enddo

end subroutine

end subroutine hdf5_read_beam

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_int_dataset(root_id, name, value, error)
!
! Private routine used by hdf5_read_beam.
!-

subroutine pmd_read_int_dataset(root_id, name, value, error)

type (hdf5_info_struct) info
type (pmd_unit_struct) unit

real(rp) unit_si

integer(HID_T) :: root_id, obj_id
integer h5_err, value(:)

logical error, err

character(*) name
character(*), parameter :: r_name = 'pmd_read_int_dataset'

!

info = hdf5_object_info(root_id, name, error, .true.)
obj_id = hdf5_open_object(root_id, name, info, error, .true.)

call hdf5_read_attribute_real(obj_id, 'unitSI', unit_si, error, .true.);  if (error) return
if (abs(unit_si - 1.0_rp) > 1d-6) then
  call out_io (s_error$, r_name, 'CONVERSION TO SI OF A VALUE OTHER THAN 1 DOES NOT MAKE SENSE.', &
                                  'FOR BUNCH PARAMETER: ' // name)
endif

!

if (info%element_type == H5O_TYPE_DATASET_F) then
  if (any(info%data_dim(2:) /= 0)) then
    call out_io (s_error$, r_name, 'DATA ARRAY IS NOT ONE-DIMENSIONAL! FOR DATA: ' // name)
    return
  endif

  call hdf5_read_dataset_int(root_id, name, value, err)

!

else  ! Must be a "constant record component" as defined by the openPMD standard
  call hdf5_read_attribute_int(obj_id, 'value', value(1), error, .true.)
  value = value(1)
endif

!

call hdf5_close_object(obj_id, info)

end subroutine pmd_read_int_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_real_dataset(root_id, name, conversion_factor, value, error)
!
! Private routine used by hdf5_read_beam.
!-

subroutine pmd_read_real_dataset(root_id, name, conversion_factor, value, error)

type (hdf5_info_struct) info

real(rp) conversion_factor, value(:)
real(rp) unit_si

integer(HID_T) :: root_id, obj_id
integer h5_err

logical error, err

character(*) name
character(*), parameter :: r_name = 'pmd_read_real_dataset'

!

info = hdf5_object_info(root_id, name, error, .true.)
obj_id = hdf5_open_object(root_id, name, info, error, .true.)

call hdf5_read_attribute_real(obj_id, 'unitSI', unit_si, error, .true.);  if (error) return

!

if (info%element_type == H5O_TYPE_DATASET_F) then
  if (any(info%data_dim(2:) /= 0)) then
    call out_io (s_error$, r_name, 'DATA ARRAY IS NOT ONE-DIMENSIONAL! FOR DATA: ' // name)
    return
  endif

  call hdf5_read_dataset_real(root_id, name, value, err)

!

else  ! Must be a "constant record component" as defined by the openPMD standard
  call hdf5_read_attribute_real(obj_id, 'value', value(1), error, .true.)
  value = value(1)
endif

!

if (abs(unit_si - conversion_factor) > 1e-15 * conversion_factor) value = value * (conversion_factor / unit_si)

call hdf5_close_object(obj_id, info)

end subroutine pmd_read_real_dataset

end module hdf5_bunch_mod
