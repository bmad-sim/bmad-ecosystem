!+
! Subroutine hdf5_read_beam (file_name, beam, error, pmd_header)
!
! Routine to read a beam data file. 
! See also hdf5_write_beam
!
! Input:
!   file_name         -- character(*): Name of the beam data file.
!
! Output:
!   beam              -- beam_struct: Particle positions.
!   error             -- logical: Set True if there is a read error. False otherwise.
!   pmd_header        -- pmd_header_struct, optional: Extra info like file creation date.
!-

subroutine hdf5_read_beam (file_name, beam, error, pmd_header)

use bmad_interface, dummy => hdf5_read_beam
use hdf5_openpmd_mod
use fortran_cpp_utils
use iso_fortran_env
use iso_c_binding

implicit none

type (beam_struct), target :: beam
type (pmd_header_struct), optional :: pmd_header
type (pmd_header_struct) :: pmd_head
type (pmd_unit_struct) unit
type(H5O_info_t) :: infobuf 
type(c_ptr) cv_ptr
type(c_funptr) c_func_ptr

real(rp), allocatable :: rvec(:)
real(rp) factor

integer(hid_t) f_id, g_id, z_id, z2_id, r_id, b_id
integer(hsize_t) idx
integer(size_t) g_size
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

call hdf5_read_attribute_alloc_string (f_id, 'openPMD', pmd_head%openPMD, err, .true.);                    if (err) return
call hdf5_read_attribute_alloc_string (f_id, 'openPMDextension', pmd_head%openPMDextension, err, .true.);  if (err) return
call hdf5_read_attribute_alloc_string (f_id, 'basePath', pmd_head%basePath, err, .true.);                  if (err) return
call hdf5_read_attribute_alloc_string (f_id, 'particlesPath', pmd_head%particlesPath, err, .true.);        if (err) return
! It is not an error if any of the following is missing
call hdf5_read_attribute_alloc_string (f_id, 'software', pmd_head%software, err, .false.)
call hdf5_read_attribute_alloc_string (f_id, 'softwareVersion', pmd_head%softwareVersion, err, .false.)
call hdf5_read_attribute_alloc_string (f_id, 'date', pmd_head%date, err, .false.)
call hdf5_read_attribute_alloc_string (f_id, 'latticeFile', pmd_head%latticeFile, err, .false.)
call hdf5_read_attribute_alloc_string (f_id, 'latticeName', pmd_head%latticeName, err, .false.)

if (present(pmd_header)) pmd_header = pmd_head

! Find root group
! Note: Right now it is assumed that the basepath uses "/%T/" to sort bunches.

it = index(pmd_head%basePath, '/%T/')
if (it /= len_trim(pmd_head%basePath) - 3) then
  call out_io (s_error$, r_name, 'PARSING OF BASE PATH FAILURE. PLEASE REPORT THIS. ' // pmd_head%basePath)
  return
endif

z_id = hdf5_open_group(f_id, pmd_head%basePath(1:it), err, .true.)

! Count bunches

n_bunch = 0
call H5Gget_info_f (z_id, storage_type, n_links, max_corder, h5_err)
do idx = 0, n_links-1
  call H5Lget_name_by_idx_f (z_id, '.', H5_INDEX_NAME_F, H5_ITER_INC_F, idx, c_name, h5_err, g_size)
  call to_f_str(c_name, name)
  call H5Oget_info_by_name_f(z_id, name, infobuf, h5_stat)
  if (infobuf%type /= H5O_TYPE_GROUP_F) cycle    ! Ignore non-group elements.
  if (.not. is_integer(name)) then               ! This assumes basepath uses "/%T/" to designate different bunches.
    call out_io (s_warn$, r_name, 'NAME OF DIRECTORY IN PATH IS NOT AN INTEGER: ' // quote(name))
    cycle
  endif
  n_bunch = n_bunch + 1
enddo

call reallocate_beam(beam, n_bunch)

! Loop over all bunches.
! Note: There is a GCC compiler bug where, if building shared, there is an error if 
! "c_funloc(hdf5_read_bunch)" is used as the actual arg to H5Literate_f in place of c_fun_ptr.

cv_ptr = c_loc(beam)
c_func_ptr = c_funloc(hdf5_read_bunch)
idx = 0
n_bunch = 0
call H5Literate_f (z_id, H5_INDEX_NAME_F, H5_ITER_INC_F, idx, c_func_ptr, cv_ptr, state, h5_err)
call H5Gclose_f(z_id, h5_err)

! And end

error = .false.

9000 continue
call h5fclose_f(f_id, h5_err)

!------------------------------------------------------------------------------------------
contains

function hdf5_read_bunch (root_id, g_name_c, info, dummy_c_ptr) result (stat) bind(C)

type(c_ptr) info
type(c_ptr) dummy_c_ptr
type(H5O_info_t) :: infobuf 
type(bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p
type(c_funptr) c_func_ptr

real(rp) charge_factor, f_ev
real(rp), allocatable :: dt(:), pz(:)

integer(hid_t), value :: root_id
integer(hid_t) g_id, g2_id, g3_id, g4_id, a_id
integer(hsize_t) idx
integer n, stat, h5_stat, ia, ip, species
integer, allocatable :: charge_state(:)

logical error

character(1) :: g_name_c(*)
character(:), allocatable :: string
character(100) g_name, a_name, name, c_name
character(*), parameter :: r_name = 'hdf5_read_bunch'

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
g2_id = hdf5_open_group(g_id, pmd_head%particlesPath, error, .true.);  if (error) return

! Get number of particles.

call hdf5_read_attribute_int(g2_id, 'numParticles', n, error, .true.)
if (error) return
call reallocate_bunch(bunch, n)
bunch%particle = coord_struct()
allocate (dt(n), charge_state(n), pz(n))

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
    call pmd_read_real_dataset(g2_id, 'momentum/z', f_ev, bunch%particle%vec(6), error)
  case ('velocity')
    call pmd_read_real_dataset(g2_id, 'velocity/x', c_light, bunch%particle%vec(2), error)
    call pmd_read_real_dataset(g2_id, 'velocity/y', c_light, bunch%particle%vec(4), error)
    call pmd_read_real_dataset(g2_id, 'velocity/z', c_light, bunch%particle%vec(6), error)
  case ('pathLength')
    call pmd_read_real_dataset(g2_id, name, 1.0_rp, bunch%particle%path_len, error)
  case ('totalMomentumOffset')
    call pmd_read_real_dataset(g2_id, name, f_ev, bunch%particle%p0c, error)
  case ('totalMomentum')
    call pmd_read_real_dataset(g2_id, name, f_ev, pz, error)
  case ('photonPolarizationAmplitude')
    call pmd_read_real_dataset(g2_id, 'photonPolarizationAmplitude/x', 1.0_rp, bunch%particle%field(1), error)
    call pmd_read_real_dataset(g2_id, 'photonPolarizationAmplitude/y', 1.0_rp, bunch%particle%field(2), error)
  case ('photonPolarizationPhase')
    call pmd_read_real_dataset(g2_id, 'photonPolarizationPhase/x', 1.0_rp, bunch%particle%phase(1), error)
    call pmd_read_real_dataset(g2_id, 'photonPolarizationPhase/y', 1.0_rp, bunch%particle%phase(2), error)
  case ('sPosition')
    call pmd_read_real_dataset(g2_id, name, 1.0_rp, bunch%particle%s, error)
  case ('time')
    call pmd_read_real_dataset(g2_id, name, 1.0_rp, dt, error)
  case ('timeOffset')
    call pmd_read_real_dataset(g2_id, name, 1.0_rp, bunch%particle%t, error)
  case ('weight')
    call pmd_read_real_dataset(g2_id, name, 1.0_rp, bunch%particle%charge, error)
  case ('particleStatus')
    call pmd_read_int_dataset(g2_id, name, 1.0_rp, bunch%particle%state, error)
  case ('chargeState')
    call pmd_read_int_dataset(g2_id, name, 1.0_rp, charge_state, error)
  case ('branchIndex')
    call pmd_read_int_dataset(g2_id, name, 1.0_rp, bunch%particle%ix_branch, error)
  case ('elementIndex')
    call pmd_read_int_dataset(g2_id, name, 1.0_rp, bunch%particle%ix_ele, error)
  case ('locationInElement')
    call pmd_read_int_dataset(g2_id, name, 1.0_rp, bunch%particle%location, error)
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
  p%t = p%t + dt(ip)
  if (species == photon$) then
    p%species = species
    p%beta = 1
  else
    p%species = set_species_charge(species, charge_state(ip))
    call convert_pc_to (sqrt(p%vec(2)**2 + p%vec(4)**2 + p%vec(6)**2), p%species, beta = p%beta)
    p%vec(5) = -p%beta * c_light * dt(ip)
    p%vec(2) = p%vec(2) / p%p0c
    p%vec(4) = p%vec(4) / p%p0c
    p%vec(6) = pz(ip) / p%p0c
  endif
enddo

end function  hdf5_read_bunch

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
