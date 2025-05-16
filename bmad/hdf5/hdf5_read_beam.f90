!+
! Subroutine hdf5_read_beam (file_name, beam, error, ele, pmd_header, print_mom_shift_warning, conserve_momentum)
!
! Routine to read a beam data file. 
! See also hdf5_write_beam.
!
! What is stored in the file is the momentum and ref momentum for each particle. If ele has a different ref p0c than
! the stored values, the particles that are created, since they inherit the ele p0c, can either have the same momentum
! as the particles there were used to create the beam file or they can have the same phase space px, py, pz values. But 
! not both. This is determined by the conserve_momentum arg.
!
! Input:
!   file_name         -- character(*): Name of the beam data file.
!                           Default is False.
!   ele               -- ele_struct, optional: Element where beam is to be started from.
!                           The element reference momentum will be used for the particle reference momentum.
!   print_mom_shift_warning   -- logical, optional: Default is True. Oly relavent if ele arg is present.
!                                 Print warning if element p0c reference momentum
!                                 is different from what is stored for the particles?
!   conserve_momentum         -- logical, optional: Default is False. Only relavent if ele arg is present and 
!                                 ele ref p0c is different from the ref p0c stored in the file.
!                                 If True and the ref p0c's differ, the output particles will have the same actual momentum 
!                                 as the particles used to create the beam file but phase space px, py, pz will differ. 
!                                 If False, phase space px, py, pz will be the same and the actual momentums can differ. See above.
!                                   
!
! Output:
!   beam              -- beam_struct: Particle positions.
!   error             -- logical: Set True if there is a read error. False otherwise.
!   pmd_header        -- pmd_header_struct, optional: Extra info like file creation date.
!-

subroutine hdf5_read_beam (file_name, beam, error, ele, pmd_header, print_mom_shift_warning, conserve_momentum)

use hdf5_openpmd_mod
use bmad_interface, dummy => hdf5_read_beam

implicit none

type (beam_struct), target :: beam
type (ele_struct), optional :: ele
type (pmd_header_struct), optional :: pmd_header
type (pmd_header_struct) :: pmd_head
type (pmd_unit_struct) unit
type(H5O_info_t) :: infobuf 

real(rp), allocatable :: rvec(:)
real(rp) factor

integer(hid_t) f_id, z_id, z2_id, r_id, b_id
integer(hsize_t) idx
integer(size_t) g_size
integer i, j, ik, ib, ix, is, it, nn, nt
integer state, h5_err, n_bunch, storage_type, n_links, max_corder, h5_stat
integer h_err
integer, allocatable :: ivec(:)

character(*) file_name
character(*), parameter :: r_name = 'hdf5_read_beam'
character(100) c_name, name, t_match, sub_dir

logical, optional :: print_mom_shift_warning, conserve_momentum
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

! Single bunch case with "%T" not present.

it = index(pmd_head%basePath, '%T')
if (it == 0) then  ! No "%T" means only one bunch present
  z_id = hdf5_open_group(f_id, pmd_head%basePath, err, .true.);  if (err) return
  call reallocate_beam (beam, 1)
  call hdf5_read_bunch(z_id, '.', beam%bunch(1), error, pmd_head, ele)
  call H5Gclose_f(z_id, h5_err)
  call H5Fclose_f(f_id, h5_err)
  return
endif

! Multiple bunch case with "%T" present.
! basePath could be something like "abc/zzz%Txxx/yyy"

is = index(pmd_head%basePath(1:it), '/', back = .true.)
z_id = hdf5_open_group(f_id, pmd_head%basePath(1:is), err, .true.);  if (err) return
t_match = pmd_head%basePath
t_match = t_match(is+1:)
i = len_trim(t_match)
if (t_match(i:i) == '/') t_match(i:i) = ' ' ! Erase trailing slash.
is = index(t_match, '/')
if (is == 0) then
  sub_dir = ''
else
  sub_dir = t_match(is:)
  t_match = t_match(1:is-1)  ! Something like "zzz%Txxx"
endif

! Loop over all bunches

n_bunch = 0
call H5Gget_info_f (z_id, storage_type, n_links, max_corder, h5_err)
do idx = 0, n_links-1
  call H5Lget_name_by_idx_f (z_id, '.', H5_INDEX_NAME_F, H5_ITER_INC_F, idx, c_name, h5_err, g_size)
  call to_f_str(c_name, name)
  call H5Oget_info_by_name_f(z_id, name, infobuf, h5_stat)
  if (infobuf%type /= H5O_TYPE_GROUP_F) cycle    ! Ignore non-group elements.
  nn = len_trim(name)
  nt = len_trim(t_match)
  it = index(t_match, '%T')
  if (nn < nt - 1) cycle
  if (it > 1) then
    if (name(1:it-1) /= t_match(1:it-1)) cycle
  endif
  j = nt - (it + 1)  ! Num char after "%T"
  if (j > 0) then
    if (name(nn-j+1:nn) /= t_match(nt-j+1:nt)) cycle
  endif
  if (.not. is_integer(name(it:nn-j))) cycle
  name = trim(name) // sub_dir
  n_bunch = n_bunch + 1
  call reallocate_beam (beam, n_bunch, 0, extend = .true.)
  call hdf5_read_bunch(z_id, name, beam%bunch(n_bunch), error, pmd_head, ele)
enddo

if (n_bunch == 0) then
  call out_io (s_error$, r_name, 'UNABLE TO LOCATE BEAM DATA IN FILE: ' // file_name)
  return
endif

call reallocate_beam(beam, n_bunch)

call H5Gclose_f(z_id, h5_err)
call H5Fclose_f(f_id, h5_err)

!------------------------------------------------------------------------------------------
contains

subroutine hdf5_read_bunch (root_id, bunch_obj_name, bunch, error, pmd_head, ele)

type (pmd_header_struct) pmd_head
type (ele_struct), optional :: ele
type(H5O_info_t) :: infobuf 
type(bunch_struct), target :: bunch
type (coord_struct), pointer :: p

real(rp) charge_factor, f_ev, p0c_initial, p0c_final
real(rp), allocatable :: dt(:), t0(:), tot_mom(:), mom_x_off(:), mom_y_off(:), mom_z_off(:)
real(rp), allocatable :: pos_x_off(:), pos_y_off(:), pos_z_off(:)

integer(hid_t), value :: root_id
integer(hid_t) g_id, g1_id, g2_id, g3_id, a_id
integer(hsize_t) idx
integer(size_t) g_size
integer n, stat, h5_stat, ia, ip, species, h5_err, storage_type, n_links, max_corder
integer, allocatable :: charge_state(:)

character(*) bunch_obj_name
character(:), allocatable :: string
character(100) g_name, a_name, name, c_name

logical error, momentum_warning_printed

! Get number of particles and init.
! Old style: Beam info tree is in g1_id directory.
! New style: Beam info tree is in subdirectory
! General note: There is no way to check names for misspellings since any name may just be an extension standard.

g_id = hdf5_open_group(root_id, bunch_obj_name, error, .true.);  if (error) return
g1_id = hdf5_open_group(g_id, pmd_head%particlesPath, error, .true.);  if (error) return
g2_id = g1_id    ! Assume old style

call hdf5_read_attribute_int(g1_id, 'numParticles', n, error, .false.)

if (error) then   ! Is new style.
  call H5Gget_info_f (g1_id, storage_type, n_links, max_corder, h5_err)
  do idx = 0, n_links-1
    call H5Lget_name_by_idx_f (g1_id, '.', H5_INDEX_NAME_F, H5_ITER_INC_F, idx, c_name, h5_err, g_size)
    call to_f_str(c_name, name)
    call H5Oget_info_by_name_f(g1_id, name, infobuf, h5_stat)
    if (infobuf%type /= H5O_TYPE_GROUP_F) cycle    ! Ignore non-group elements.
    g3_id = hdf5_open_group(g1_id, name, error, .false.)
    call hdf5_read_attribute_int(g3_id, 'numParticles', n, error, .false.)
    if (error) then
      call H5Gclose_f(g3_id, h5_err)
      cycle
    endif
    g2_id = g3_id
    exit
  enddo

  if (g2_id == g1_id) then ! Not found
    call hdf5_read_attribute_int(g2_id, 'numParticles', n, error, .true.)  ! Generate error message
    return
  endif
endif

!

allocate (dt(n), t0(n), charge_state(n), tot_mom(n), mom_x_off(n), mom_y_off(n), mom_z_off(n))
allocate (pos_x_off(n), pos_y_off(n), pos_z_off(n))

momentum_warning_printed = .false.
charge_factor = 0
species = int_garbage$  ! Garbage number
tot_mom = real_garbage$
charge_state = 0
mom_x_off = 0;  mom_y_off = 0;  mom_z_off = 0
pos_x_off = 0;  pos_y_off = 0;  pos_z_off = 0
t0 = 0; dt = 0

call reallocate_bunch(bunch, n)
bunch%particle = coord_struct()
bunch%particle%state = alive$  ! alive$ = 1 which is same as OpenPMD standard.
bunch%particle%charge = 1

bunch%charge_tot = 0
bunch%charge_live = 0

if (present(ele)) then
  bunch%particle%ix_ele = ele%ix_ele
  bunch%particle%ix_branch = ele%ix_branch
  bunch%particle%s = ele%s
  bunch%particle%location = exit_end$
  if (associated(ele%branch)) species = ele%branch%param%particle
endif

! Get attributes.

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
  end select
enddo

if (charge_factor == 0 .and. (bunch%charge_live /= 0 .or. bunch%charge_tot /= 0)) then
  call out_io (s_error$, r_name, 'chargeUnitSI is not set for bunch')
else
  bunch%charge_live = bunch%charge_live * charge_factor
  bunch%charge_tot = bunch%charge_tot * charge_factor
endif

! Loop over all datasets.

f_ev = unit_ev_per_c%unitSI

call H5Gget_info_f (g2_id, storage_type, n_links, max_corder, h5_err)
do idx = 0, n_links-1
  call H5Lget_name_by_idx_f (g2_id, '.', H5_INDEX_NAME_F, H5_ITER_INC_F, idx, name, h5_err, g_size)
  call H5Oget_info_by_name_f(g2_id, name, infobuf, h5_stat)
  select case (name)
  case ('spin')
    call pmd_read_real_dataset(g2_id, 'spin/x', unit_hbar%unitSI, bunch%particle%spin(1), error)
    call pmd_read_real_dataset(g2_id, 'spin/y', unit_hbar%unitSI, bunch%particle%spin(2), error)
    call pmd_read_real_dataset(g2_id, 'spin/z', unit_hbar%unitSI, bunch%particle%spin(3), error)
  case ('position')
    call pmd_read_real_dataset(g2_id, 'position/x', 1.0_rp, bunch%particle%vec(1), error)
    call pmd_read_real_dataset(g2_id, 'position/y', 1.0_rp, bunch%particle%vec(3), error)
    call pmd_read_real_dataset(g2_id, 'position/z', 1.0_rp, bunch%particle%vec(5), error)
  case ('positionOffset')
    call pmd_read_real_dataset(g2_id, 'positionOffset/x', 1.0_rp, pos_x_off, error)
    call pmd_read_real_dataset(g2_id, 'positionOffset/y', 1.0_rp, pos_y_off, error)
    call pmd_read_real_dataset(g2_id, 'positionOffset/z', 1.0_rp, pos_z_off, error)
  case ('momentum')
    call pmd_read_real_dataset(g2_id, 'momentum/x', f_ev, bunch%particle%vec(2), error)
    call pmd_read_real_dataset(g2_id, 'momentum/y', f_ev, bunch%particle%vec(4), error)
    call pmd_read_real_dataset(g2_id, 'momentum/z', f_ev, bunch%particle%vec(6), error)
  case ('momentumOffset')
    call pmd_read_real_dataset(g2_id, 'momentumOffset/x', f_ev, mom_x_off, error)
    call pmd_read_real_dataset(g2_id, 'momentumOffset/y', f_ev, mom_y_off, error)
    call pmd_read_real_dataset(g2_id, 'momentumOffset/z', f_ev, mom_z_off, error)
  case ('velocity')
    call pmd_read_real_dataset(g2_id, 'velocity/x', c_light, bunch%particle%vec(2), error)
    call pmd_read_real_dataset(g2_id, 'velocity/y', c_light, bunch%particle%vec(4), error)
    call pmd_read_real_dataset(g2_id, 'velocity/z', c_light, bunch%particle%vec(6), error)
  case ('velocityOffset')
    call pmd_read_real_dataset(g2_id, 'velocityOffset/x', c_light, mom_x_off, error)
    call pmd_read_real_dataset(g2_id, 'velocityOffset/y', c_light, mom_y_off, error)
    call pmd_read_real_dataset(g2_id, 'velocityOffset/z', c_light, mom_z_off, error)
  case ('pathLength')
    call pmd_read_real_dataset(g2_id, name, 1.0_rp, bunch%particle%dt_ref, error)
    bunch%particle%dt_ref = bunch%particle%dt_ref/c_light
  case ('photonPolarizationAmplitude')
    call pmd_read_real_dataset(g2_id, 'photonPolarizationAmplitude/x', 1.0_rp, bunch%particle%field(1), error)
    call pmd_read_real_dataset(g2_id, 'photonPolarizationAmplitude/y', 1.0_rp, bunch%particle%field(2), error)
  case ('photonPolarizationPhase')
    call pmd_read_real_dataset(g2_id, 'photonPolarizationPhase/x', 1.0_rp, bunch%particle%phase(1), error)
    call pmd_read_real_dataset(g2_id, 'photonPolarizationPhase/y', 1.0_rp, bunch%particle%phase(2), error)
  case ('sPosition')
    if (.not. present(ele)) call pmd_read_real_dataset(g2_id, name, 1.0_rp, bunch%particle%s, error)
  case ('time')
    call pmd_read_real_dataset(g2_id, name, 1.0_rp, dt, error)
  case ('timeOffset')
    call pmd_read_real_dataset(g2_id, name, 1.0_rp, t0, error)
  case ('totalMomentum')
    call pmd_read_real_dataset(g2_id, name, f_ev, tot_mom, error)
  case ('totalMomentumOffset')
    call pmd_read_real_dataset(g2_id, name, f_ev, bunch%particle%p0c, error)
  case ('weight')
    call pmd_read_real_dataset(g2_id, name, 1.0_rp, bunch%particle%charge, error)
  case ('particleStatus')
    call pmd_read_int_dataset(g2_id, name, 1.0_rp, bunch%particle%state, error)
  case ('chargeState')
    call pmd_read_int_dataset(g2_id, name, 1.0_rp, charge_state, error)
  case ('branchIndex')
    if (.not. present(ele)) call pmd_read_int_dataset(g2_id, name, 1.0_rp, bunch%particle%ix_branch, error)
  case ('elementIndex')
    if (.not. present(ele)) call pmd_read_int_dataset(g2_id, name, 1.0_rp, bunch%particle%ix_ele, error)
  case ('locationInElement')
    if (.not. present(ele)) then
      call pmd_read_int_dataset(g2_id, name, 1.0_rp, bunch%particle%location, error)
      do ip = 1, size(bunch%particle)
        p => bunch%particle(ip)
        select case(p%location)
        case (-1);    p%location = upstream_end$
        case ( 0);    p%location = inside$
        case ( 1);    p%location = downstream_end$
        end select
      enddo
    endif
  end select

  if (error) exit
enddo

! g_id = g2_id when pmd_head%particlesPath = "./". In this case both refer to the same group.

if (g2_id /= root_id)                     call H5Gclose_f(g2_id, h5_err)
if (g_id /= g2_id .and. g_id /= root_id)  call H5Gclose_f(g_id, h5_err)
if (g1_id /= g2_id)                       call H5Gclose_f(g1_id, h5_err)

if (error) then
  call out_io (s_error$, r_name, 'ERROR READING BEAM DATA FROM: ' // file_name, 'ABORTING READ...')
  return
endif

if (species == int_garbage$) then
  call out_io (s_error$, r_name, 'SPECIES OF PARTICLE NOT PRESENT WHILE READING BEAM DATA FROM: ' // file_name, 'ABORTING READ...')
  error = .true.
  return
endif

bunch%particle%vec(1) = bunch%particle%vec(1) + pos_x_off
bunch%particle%vec(3) = bunch%particle%vec(3) + pos_y_off
bunch%particle%vec(5) = bunch%particle%vec(5) + pos_z_off

bunch%particle%vec(2) = bunch%particle%vec(2) + mom_x_off
bunch%particle%vec(4) = bunch%particle%vec(4) + mom_y_off
bunch%particle%vec(6) = bunch%particle%vec(6) + mom_z_off

do ip = 1, size(bunch%particle)
  p => bunch%particle(ip)
  p%t = t0(ip) + dt(ip)

  if (species == photon$) then
    p%species = species
    p%beta = 1
    cycle
  endif

  ! Not a photon.

  p0c_initial = p%p0c
  p0c_final   = p%p0c
  if (present(ele)) then
    if (p0c_initial == 0) p0c_initial   = ele%value(p0c$)
    p0c_final   = ele%value(p0c$)
  endif

  if (p0c_final == 0) then
    call out_io (s_error$, r_name, 'REFERENCE MOMENTUM NOT PRESENT WHILE READING BEAM DATA FROM: ' // file_name, 'ABORTING READ...')
    error = .true.
    return
  endif

  select case (p%state)
  case (alive$)
    if (abs(p0c_initial - p0c_final) > 1e-12 * p0c_final .and. &
            logic_option(.true., print_mom_shift_warning) .and. .not. momentum_warning_printed) then
      if (logic_option(.false., conserve_momentum)) then
        call out_io (s_warn$, r_name, 'REFERENCE MOMENTUM OF PARTICLE IN BEAM FILE:  \es20.12\ ', &
                                      'FROM FILE: ' // file_name, &
                                      'DIFFERENT FROM REFERNECE MOMENTUM IN LATTICE: \es20.12\ ', &
                                      'THIS WILL CAUSE A SHIFT IN PARTICLE''S MOMENTUM (BUT NOT PZ)', &
                                      r_array = [p0c_initial, p0c_final])
      else
        call out_io (s_warn$, r_name, 'REFERENCE MOMENTUM OF PARTICLE IN BEAM FILE:  \es20.12\ ', &
                                      'FROM FILE: ' // file_name, &
                                      'DIFFERENT FROM REFERNECE MOMENTUM IN LATTICE: \es20.12\ ', &
                                      'THIS WILL CAUSE A SHIFT IN PARTICLE''S PHASE SPACE pz = (P - P0)/P (BUT NOT THE MOMENTUM)', &
                                      r_array = [p0c_initial, p0c_final])
      endif
    endif
    momentum_warning_printed = .true.
  case (pre_born$)
    ! Nothing to be done

  case default
    if (pmd_head%software /= 'Bmad') p%state = lost$
  end select

  p%species = set_species_charge(species, charge_state(ip))

  if (tot_mom(ip) == real_garbage$ .or. p0c_initial == 0) then 
    p%vec(6) = (sqrt(p%vec(2)**2 + p%vec(4)**2 + p%vec(6)**2) - p0c_final) / p0c_final
  elseif (logic_option(.false., conserve_momentum)) then
    p%vec(6) = (tot_mom(ip) + (p0c_initial - p0c_final)) / p0c_final
  else
    p%vec(6) = tot_mom(ip) / p0c_initial
  endif

  if (logic_option(.true., conserve_momentum)) then
    call convert_pc_to ((1 + p%vec(6)) * p0c_final, p%species, beta = p%beta)
    p%vec(2) = p%vec(2) / p0c_final
    p%vec(4) = p%vec(4) / p0c_final
    p%vec(5) = -p%beta * c_light * dt(ip)
  else
    call convert_pc_to ((1 + p%vec(6)) * p0c_initial, p%species, beta = p%beta)
    p%vec(2) = p%vec(2) / p0c_initial
    p%vec(4) = p%vec(4) / p0c_initial
    p%vec(5) = -p%beta * c_light * dt(ip)
  endif

  p%p0c = p0c_final
enddo

end subroutine  hdf5_read_bunch

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

end subroutine to_amp_phase

end subroutine hdf5_read_beam
