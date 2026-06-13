!+
! Subroutine write_lattice_pals_format(pals_file, lat, err_flag)
!
! Routine to create a Pals lattice file.
!
! Input:
!   lat           -- lat_struct: Lattice
!
! Output:
!   pals_file  -- character(*): Pals lattice file name.
!   err_flag      -- logical: Error flag
!-

subroutine write_lattice_pals_format(pals_file, lat, err_flag)

use write_lattice_file_mod, dummy => write_lattice_pals_format
use bmad_routine_interface, dummy2 => write_lattice_pals_format
use expression_mod
use taylor_mod, only: mat6_to_taylor

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, lord, slave, slave2, multi_lord
type (ele_struct) ele_dflt
type (coord_struct), pointer :: orb
type (multipass_region_lat_struct), target :: mult_lat
type (multipass_all_info_struct), target :: m_info
type (multipass_region_ele_struct), pointer :: mult_ele(:), m_ele
type (multipass_ele_info_struct), pointer :: e_info
type (ele_pointer_struct), allocatable :: named_eles_ptr(:)  ! List of unique element names 
type (lat_ele_order_struct) order
type (ele_attribute_struct) info
type (taylor_struct) taylor(6), spin_taylor(0:3)
type (nametable_struct) var_nametab
type (control_struct), pointer :: ctl
type (control_struct) control

real(rp) f, length, ang2
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)

integer n, i, j, k, ix, ib, ie, iu, is, n_names, ix_match, ix_pass, ix_r
integer :: ix_lord, ix_super, ie1, ib1, id, id1, id2, id_del = 2
integer, allocatable :: an_indexx(:), index_list(:)

logical header_needed, has_been_added, in_multi_region, have_expand_lattice_line, err, is_added
logical header_need1, header_need2
logical, optional :: err_flag

character(*) pals_file
character(1) prefix
character(3), parameter :: unit_spin_map(0:3) = ['1.0', '0.0', '0.0', '0.0']
character(100) :: name, look_for, ele_name, blank = ''
character(40), allocatable :: names(:)
character(240) fname
character(1000) line, bline
character(*), parameter :: r_name = 'write_lattice_pals_format'

character(20) :: pals_ele_type(n_key$)

!

pals_ele_type(ab_multipole$)         = 'Multipole'
pals_ele_type(ac_kicker$)            = 'ACKicker'
pals_ele_type(beambeam$)             = 'BeamBeam'
pals_ele_type(beginning_ele$)        = 'BeginningEle'
pals_ele_type(capillary$)            = 'Capillary??'
pals_ele_type(converter$)            = 'Converter'
pals_ele_type(crab_cavity$)          = 'CrabCavity'
pals_ele_type(crystal$)              = 'Crystal??'
pals_ele_type(custom$)               = 'Custom'
pals_ele_type(def_line$)             = '!Line'
pals_ele_type(def_ptc_com$)          = '!PTC_Com'
pals_ele_type(detector$)             = 'Instrument'
pals_ele_type(diffraction_plate$)    = 'DiffractionPlate'
pals_ele_type(drift$)                = 'Drift'
pals_ele_type(e_gun$)                = 'EGun'
pals_ele_type(ecollimator$)          = 'Drift'
pals_ele_type(elseparator$)          = 'ELSeparator??'
pals_ele_type(em_field$)             = 'EMField'
pals_ele_type(feedback$)             = 'Instrument'
pals_ele_type(fiducial$)             = 'Fiducial'
pals_ele_type(fixer$)                = 'Fixer'
pals_ele_type(floor_shift$)          = 'FloorShift'
pals_ele_type(foil$)                 = 'Foil'
pals_ele_type(fork$)                 = 'Fork'
pals_ele_type(girder$)               = 'Girder'
pals_ele_type(gkicker$)              = 'Kicker'
pals_ele_type(group$)                = '??Group'
pals_ele_type(hkicker$)              = 'Kicker'
pals_ele_type(hybrid$)               = 'Hybrid??'
pals_ele_type(instrument$)           = 'Instrument'
pals_ele_type(kicker$)               = 'Kicker'
pals_ele_type(lcavity$)              = 'RFCavity'
pals_ele_type(lens$)                 = 'Lens'
pals_ele_type(marker$)               = 'Marker'
pals_ele_type(mask$)                 = 'Mask'
pals_ele_type(match$)                = 'Match'
pals_ele_type(mirror$)               = 'Mirror'
pals_ele_type(monitor$)              = 'Instrument'
pals_ele_type(multilayer_mirror$)    = 'MultilayerMirror'
pals_ele_type(multipole$)            = 'Multipole'
pals_ele_type(null_ele$)             = 'Placeholder'
pals_ele_type(octupole$)             = 'Octupole'
pals_ele_type(overlay$)              = '??Overlay'
pals_ele_type(patch$)                = 'Patch'
pals_ele_type(photon_fork$)          = 'Fork'
pals_ele_type(photon_init$)          = 'PhotonInit'
pals_ele_type(pickup$)               = 'Instrument'
pals_ele_type(pipe$)                 = 'Drift'
pals_ele_type(quadrupole$)           = 'Quadrupole'
pals_ele_type(ramper$)               = 'Ramper'
pals_ele_type(rbend$)                = 'SBend'
pals_ele_type(rcollimator$)          = 'Drift'
pals_ele_type(rf_bend$)              = 'RFBend'
pals_ele_type(rfcavity$)             = 'RFCavity'
pals_ele_type(sad_mult$)             = 'SadMult'
pals_ele_type(sample$)               = 'Sample'
pals_ele_type(sbend$)                = 'SBend'
pals_ele_type(sextupole$)            = 'Sextupole'
pals_ele_type(sol_quad$)             = 'Solenoid'
pals_ele_type(solenoid$)             = 'Solenoid'
pals_ele_type(taylor$)               = 'Taylor'
pals_ele_type(thick_multipole$)      = 'ThickMultipole'
pals_ele_type(undulator$)            = 'Undulator'
pals_ele_type(vkicker$)              = 'Kicker'
pals_ele_type(wiggler$)              = 'Wiggler'

! Open file

call fullfilename(pals_file, fname)
iu = lunget()
open (iu, file = fname, status = 'unknown')

! Header

write (iu, '(a)')  'PALS:'
write (iu, '(a)')  '  notes:'
write (iu, '(4a)') '    - "Translated from Bmad lattice file: ', trim(lat%input_file_name), '"'
write (iu, '(a)')  ''
write (iu, '(a)')  '  extension_names:'
write (iu, '(a)')  '    names:'
write (iu, '(a)')  '      Bmad'
write (iu, '(a)')  '    prefixes:'
write (iu, '(a)')  '      Bmad_'
write (iu, '(a)')  ''

!

id = id_del

header_need1 = .true.
header_need2 = .true.
do ie = lat%n_ele_track+1, lat%n_ele_max
  ele => lat%ele(ie)
  select case (ele%lord_status)
  case (group_lord$, overlay_lord$)
    line = ''
    do j = 1, size(ele%control%var)
      line = trim(line) // trim(ele%control%var(j)%name) // ','
    enddo
    if (ele%lord_status == group_lord$) then
      call write_str_dict(id, 'Bmad_groups', ele%name, downcase(line), header_need1)
    else
      call write_str_dict(id, 'Bmad_overlays', ele%name, downcase(line), header_need2)
    endif
  end select
enddo

!

id = 2 * id_del
id1 = id + id_del
id2 = id + 2*id_del

call write_str(0, '')
call write_str(id_del, '#---------------------------------------------------------------------------------------')
call write_str(id_del, '# Constants')
call write_str(0, '')
write (iu, '(a)')  '  facility:'

!----------------------------
! Write element defs

! Note: Beamlines cannot currently handle multipass nor superimpose so ignore.
! Stuff that is commented out due to this is marked by "!!!"

n_names = 0
n = lat%n_ele_max
allocate (names(n), an_indexx(n), named_eles_ptr(n))

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  ele_loop: do ie = 0, branch%n_ele_track   !!! Note: Not n_ele_max since superimpose/multipass not handled
    ele => branch%ele(ie)
    length = ele%value(l$)
    if (ele%key == overlay$ .or. ele%key == group$ .or. ele%key == ramper$ .or. ele%key == girder$) cycle   ! Not currently handled
    if (ele%key == null_ele$) cycle

    call init_ele (ele_dflt, ele%key)

    ! Do not write anything for elements that have a duplicate name.

    call add_this_name_to_list (ele, names, an_indexx, n_names, ix_match, has_been_added, named_eles_ptr)
    if (.not. has_been_added) cycle
    ele_name = pals_ele_name(ele)
    if (ie == 0) ele_name = 'beginning_b' // int_str(ib+1)
    write (iu, '(a)')
    call write_list_header(id, ele_name)

    call write_str_dict(id1, '', 'kind', pals_ele_type(ele%key))
    call write_logic_dict(id1, '', 'is_on', ele%is_on, .false., .true.)

    call write_str_dict(id1, '', 'label', quote(ele%type), header_needed)
    call write_str_dict(id1, '', 'alias', quote(ele%alias), header_needed)
    if (associated(ele%descrip)) call write_str_dict(id1, '', 'description', quote(ele%descrip))

    if (has_attribute(ele, 'L') .and. length /= 0) call write_real_dict(id1, '', 'length', length)

    if (ie == 0) then
      header_needed = .true.
      call write_real_dict(id2, 'ReferenceP', 'pc_ref', ele%value(p0c$), header_needed)
      call write_str_dict(id2, 'ReferenceP', 'species_ref: ', quote(openpmd_species_name(ele%ref_species)), header_needed)

      header_needed = .true.
      call write_real_dict(id2, 'TwissP', 'beta_a',  ele%a%beta, header_needed)
      call write_real_dict(id2, 'TwissP', 'beta_b',  ele%b%beta, header_needed)
      call write_real_dict(id2, 'TwissP', 'alpha_a', ele%a%alpha, header_needed)
      call write_real_dict(id2, 'TwissP', 'alpha_b', ele%b%alpha, header_needed)
      call write_real_dict(id2, 'TwissP', 'eta_x',   ele%x%eta, header_needed)
      call write_real_dict(id2, 'TwissP', 'eta_y',   ele%y%eta, header_needed)
      call write_real_dict(id2, 'TwissP', 'etap_x',  ele%x%etap, header_needed)
      call write_real_dict(id2, 'TwissP', 'etap_y',  ele%y%etap, header_needed)
      call write_real_dict(id2, 'TwissP', 'c_mat11', ele%c_mat(1,1), header_needed)
      call write_real_dict(id2, 'TwissP', 'c_mat12', ele%c_mat(1,2), header_needed)
      call write_real_dict(id2, 'TwissP', 'c_mat21', ele%c_mat(2,1), header_needed)
      call write_real_dict(id2, 'TwissP', 'c_mat22', ele%c_mat(2,2), header_needed)

      orb => lat%particle_start
      header_needed = .true.
      call write_real_dict(id2, 'ParticleP', 'x',  orb%vec(1), header_needed)
      call write_real_dict(id2, 'ParticleP', 'px', orb%vec(2), header_needed)
      call write_real_dict(id2, 'ParticleP', 'y',  orb%vec(3), header_needed)
      call write_real_dict(id2, 'ParticleP', 'py', orb%vec(4), header_needed)
      call write_real_dict(id2, 'ParticleP', 'z',  orb%vec(5), header_needed)
      call write_real_dict(id2, 'ParticleP', 'pz', orb%vec(6), header_needed)
      call write_real_dict(id2, 'ParticleP', 'spin_x', orb%spin(1), header_needed)
      call write_real_dict(id2, 'ParticleP', 'spin_y', orb%spin(2), header_needed)
      call write_real_dict(id2, 'ParticleP', 'spin_z', orb%spin(3), header_needed)
    endif

    !

    if (ele%key == sbend$) then
      header_needed = .true.
      call write_real_dict(id2, 'BendP', 'g_ref',     ele%value(g$), header_needed)
      call write_real_dict(id2, 'BendP', 'tilt_ref',  ele%value(ref_tilt$), header_needed)
      call write_real_dict(id2, 'BendP', 'h1',        ele%value(h1$), header_needed)
      call write_real_dict(id2, 'BendP', 'h2',        ele%value(h2$), header_needed)
      call write_real_dict(id2, 'BendP', 'edge_int1', ele%value(fint$)*ele%value(hgap$), header_needed)
      call write_real_dict(id2, 'BendP', 'edge_int2', ele%value(fintx$)*ele%value(hgapx$), header_needed)
      if (ele%sub_key == sbend$) then
        call write_real_dict(id2, 'BendP', 'e1',        ele%value(e1$), header_needed)
        call write_real_dict(id2, 'BendP', 'e2',        ele%value(e2$), header_needed)
      else
        call write_real_dict(id2, 'BendP', 'e1_rect',   ele%value(e1$)-0.5_rp*ele%value(angle$), header_needed)
        call write_real_dict(id2, 'BendP', 'e2_rect',   ele%value(e2$)-0.5_rp*ele%value(angle$), header_needed)
      endif
    endif

    ! Magnetic multipoles

    call multipole_ele_to_ab(ele, .false., ix, a_pole, b_pole, magnetic$, include_kicks$)
    if (ele%key == sbend$) then 
      b_pole(0) = b_pole(0) + ele%value(angle$)
      ix = max(0, ix)
    endif

    if (ele%field_master) then
      f = ele%value(p0c$) / (charge_of(ele%ref_species) * c_light)
      prefix = 'B'
    else
      f = 1
      prefix = 'K'
    endif

    if (length /= 0) f = f / length
    header_needed = .true.

    do j = 0, ix
      if (length == 0) then
        call write_real_dict(id2, 'MagneticMultipoleP', prefix // 's' // int_str(j) // 'L', f * factorial(j) * a_pole(j), header_needed)
        call write_real_dict(id2, 'MagneticMultipoleP', prefix // 'n' // int_str(j) // 'L', f * factorial(j) * b_pole(j), header_needed)
      else
        call write_real_dict(id2, 'MagneticMultipoleP', prefix // 's' // int_str(j), f * factorial(j) * a_pole(j), header_needed)
        call write_real_dict(id2, 'MagneticMultipoleP', prefix // 'n' // int_str(j), f * factorial(j) * b_pole(j), header_needed)
      endif
    enddo

    ! Electric multipoles

    call multipole_ele_to_ab(ele, .false., ix, a_pole, b_pole, electric$, include_kicks$)
    header_needed = .true.

    do j = 0, ix
      call write_real_dict(id2, 'ElectricMultipoleP', 'Es' // int_str(j), factorial(j) * a_pole(j), header_needed)
      call write_real_dict(id2, 'ElectricMultipoleP', 'En' // int_str(j), factorial(j) * b_pole(j), header_needed)
    enddo

    !

    if (has_attribute(ele, 'X1_LIMIT')) then
      header_needed = .true.
      call write_real_dict(id2,   'ApertureP', 'x_min', -ele%value(x1_limit$), header_needed)
      call write_real_dict(id2,   'ApertureP', 'x_max',  ele%value(x2_limit$), header_needed)
      call write_real_dict(id2,   'ApertureP', 'y_min', -ele%value(y1_limit$), header_needed)
      call write_real_dict(id2,   'ApertureP', 'y_max',  ele%value(y2_limit$), header_needed)
      if (header_needed) then
        call write_switch_dict(id2, 'ApertureP', 'shape', 1, [character(16):: 'Auto', 'RECTANGULAR', 'ELLIPTICAL'], ele%aperture_type, header_needed)
        call write_switch_dict(id2, 'ApertureP', 'location', 0, [character(16):: 'ENTRANCE_END', 'EXIT_END', 'BOTH_ENDS', 'NOWHERE', 'EVERYWHERE'], ele%aperture_at, header_needed)
        call write_logic_dict(id2,  'ApertureP', 'aperture_shifts_with_body', ele%offset_moves_aperture, header_needed, .false.)
      endif
    endif

    !

    select case (ele%key)
    case (match$, taylor$)
      length = ele%value(l$)

      select case (ele%key)
      case (match$)
        call mat6_to_taylor(ele%vec0, ele%mat6, taylor)
        call write_this_taylor(id2, ele, taylor)
        cycle

      case (taylor$)
        call write_this_taylor(id2, ele, ele%taylor)
        cycle
      end select
    end select

    !

    if (ele%key == patch$) then
      header_needed = .true.
      call write_real_dict(id2,   'PatchP', 'x_offset', ele%value(x_offset$), header_needed)
      call write_real_dict(id2,   'PatchP', 'y_offset', ele%value(y_offset$), header_needed)
      call write_real_dict(id2,   'PatchP', 'z_offset', ele%value(z_offset$), header_needed)
      call write_real_dict(id2,   'PatchP', 'x_rot',   -ele%value(y_pitch$), header_needed)
      call write_real_dict(id2,   'PatchP', 'y_rot',    ele%value(x_pitch$), header_needed)
      call write_real_dict(id2,   'PatchP', 'z_rot',    ele%value(tilt$), header_needed)
      call write_logic_dict(id2,  'PatchP', 'flexible', is_true(ele%value(flexible$)), header_needed, .false.)
      call write_logic_dict(id2,  'PatchP', 'user_sets_length', is_true(ele%value(user_sets_length$)), header_needed, .false.)
      call write_switch_dict(id2, 'PatchP', 'ref_coords', 2, [character(16):: 'entrance_end ', 'exit_end'], nint(ele%value(ref_coords$)), header_needed)


      header_needed = .true.
      call write_real_dict(id2, 'ReferenceChangeP', 'dE_ref',    ele%value(E_tot_offset$), header_needed)
      call write_real_dict(id2, 'ReferenceChangeP', 'E_tot_ref', ele%value(E_tot_set$), header_needed)
      call write_real_dict(id2, 'ReferenceChangeP', 't_offset',  ele%value(t_offset$), header_needed)


    elseif (has_attribute(ele, 'X_PITCH')) then
      header_needed = .true.
      call write_real_dict(id2, 'BodyShiftP', 'x_offset', ele%value(x_offset$), header_needed)
      call write_real_dict(id2, 'BodyShiftP', 'y_offset', ele%value(y_offset$), header_needed)
      call write_real_dict(id2, 'BodyShiftP', 'z_offset', ele%value(z_offset$), header_needed)
      call write_real_dict(id2, 'BodyShiftP', 'x_rot',   -ele%value(y_pitch$), header_needed)
      call write_real_dict(id2, 'BodyShiftP', 'y_rot',    ele%value(x_pitch$), header_needed)
      call write_real_dict(id2, 'BodyShiftP', 'z_rot',    ele%value(tilt$), header_needed)
    endif

    !

    if (has_attribute(ele, 'KS')) then
      header_needed = .true.
      if (ele%field_master) then
        call write_real_dict(id2, 'SolenoidP', 'Bsol', ele%value(bs_field$), header_needed)
      else
        call write_real_dict(id2, 'SolenoidP', 'Ksol', ele%value(ks$), header_needed)
      endif
    endif

    !

    if (has_attribute(ele, 'RF_FREQUENCY')) then
      header_needed = .true.

      if (is_true(ele%value(harmon_master$))) then
        call write_real_dict(id2, 'RFP', 'rf_frequency', ele%value(rf_frequency$), header_needed)
      else
        call write_real_dict(id2, 'RFP', 'rf_frequency', ele%value(rf_frequency$), header_needed)
      endif

      if (ele%key == lcavity$) then
        if (ele%field_master) then
          call write_real_dict(id2, 'RFP', 'gradient', ele%value(gradient$) + ele%value(gradient_err$), header_needed)
        else
          call write_real_dict(id2, 'RFP', 'voltage', ele%value(voltage$) + ele%value(voltage_err$), header_needed)
        endif

        call write_real_dict(id2, 'RFP', 'phi0', ele%value(phi0$) + ele%value(phi0_err$), header_needed)
        call write_real_dict(id2, 'RFP', 'num_cells', ele%value(n_rf_steps$), header_needed)
        call write_real_dict(id2, 'RFP', 'L_active', ele%value(l_active$), header_needed)

      else
        if (ele%field_master) then
          call write_real_dict(id2, 'RFP', 'gradient', ele%value(gradient$), header_needed)
        else
          call write_real_dict(id2, 'RFP', 'voltage', ele%value(voltage$), header_needed)
        endif

        call write_real_dict(id2, 'RFP', 'phi0', ele%value(phi0$), header_needed)
        call write_str_dict(id2, 'RFP', 'zero_phase', 'ABOVE_TRANSITION', header_needed)
      endif


      if (has_attribute(ele, 'CAVITY_TYPE')) then
        call write_switch_dict(id2, 'RFP', 'cavity_type', 0, [character(20):: 'STANDING_WAVE', 'TRAVELING_WAVE', 'STANDING_WAVE'], nint(ele%value(cavity_type$)), header_needed)
      endif
    endif

    header_needed = .true.
    if (key_name(ele%key) /= pals_ele_type(ele%key)) call write_str_dict(id2, 'BmadP', 'Bmand_key', key_name(ele%key), header_needed)
    if (has_attribute(ele, 'TRACKING_METHOD') .and. (ele%tracking_method /= ele_dflt%tracking_method)) &
                                call write_str_dict(id2, 'BmadP', 'tracking_method', tracking_method_name(ele%key), header_needed)
    if (has_attribute(ele, 'MAT6_CALC_METHOD') .and. (ele%mat6_calc_method /= ele_dflt%mat6_calc_method)) &
                                call write_str_dict(id2, 'BmadP', 'mat6_calc_method', mat6_calc_method_name(ele%key), header_needed)

  enddo ele_loop
enddo

!------------------------------------------------------------------------------------------------------
! Write branch lines
! First write multipass lines

call multipass_region_info(lat, mult_lat, m_info)

bline = 'line ['
do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  mult_ele => mult_lat%branch(ib)%ele
  in_multi_region = .false.

  do ie = 0, branch%n_ele_track
    ele => branch%ele(ie)
    ix_pass = m_info%branch(ib)%ele(ie)%ix_pass
    if (ix_pass /= 1) cycle 

    if (mult_ele(ie)%region_start_pt) then
      if (in_multi_region) then
        call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #1! PLEASE REPORT THIS!')
      endif
      in_multi_region = .true.
      ix_r = mult_ele(ie)%ix_region
      write (iu, '(a)')
      call write_list_header(id, 'BeamLine')
      call write_str_list(id, 'BeamLine', 'name', 'multi_line_' // int_str(ix_r))
      call write_str_list(id, 'BeamLine', 'multipass', 'true')
    endif

    if (mult_ele(ie)%ix_region /= ix_r) then
      call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #2! PLEASE REPORT THIS!')
    endif

    call write_element_to_beamline (id2+id_del, ele, lat, bline)

    if (mult_ele(ie)%region_stop_pt) in_multi_region = .false.
  enddo

  if (in_multi_region) then
    call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #3! PLEASE REPORT THIS!')
  endif
enddo  ! ib branch loop

!------------------------------
! Lines for all the branches.
! If we get into a multipass region then name in the main_line list is "multi_line_nn".
! But only write this once.

id1 = id + id_del
id2 = id + 2*id_del

call write_str(0, '')
call write_str(id, '#---------------------------------------------------------------------------------------')
call write_str(id, '# BeamLines and Lattices')

bline = 'line: ['
do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  name = downcase(branch%name)
  write (iu, '(a)')
  call write_list_header(id, name)
  call write_str_dict(id, '-', 'kind', 'BeamLine')

  in_multi_region = .false.
  do ie = 0, branch%n_ele_track
    ele => branch%ele(ie)

    e_info => m_info%branch(ib)%ele(ie)

    if (.not. e_info%multipass) then
      call write_element_to_beamline (id2, ele, lat, bline)
      cycle
    endif

    ix_lord = e_info%ix_lord(1)
    ix_super = e_info%ix_super(1)
    ie1 = m_info%lord(ix_lord)%slave(1,ix_super)%ele%ix_ele
    ib1 = m_info%lord(ix_lord)%slave(1,ix_super)%ele%ix_branch
    m_ele => mult_lat%branch(ib1)%ele(ie1)
    ix_r = m_ele%ix_region

    ! If entering new multipass region
    if (.not. in_multi_region) then
      in_multi_region = .true.
        call write_list_header(id1, 'multi_line_' // int_str(ix_r))
      if (m_ele%region_start_pt) then
        look_for = 'stop'
      else
        call write_str_list(id1+id_del, '', 'direction', '-1')
        look_for = 'start'
      endif
    endif

    if (look_for == 'start' .and. m_ele%region_start_pt .or. &
        look_for == 'stop' .and. m_ele%region_stop_pt) then 
      in_multi_region = .false.
    endif
  enddo
enddo

bline = trim(bline) // ']'
if (bline(1:5) == 'line:') then
  call write_str(id2, bline)
else
  call write_str(id2+6, bline)
endif

! Define lat

if (.false.) then
  line = 'lat = expand(' // quote(downcase(lat%use_name)) // '- ['
  do ib = 0, ubound(lat%branch, 1)
    branch => lat%branch(ib)
    if (branch%ix_from_branch > -1) cycle
    name = downcase(branch%name)
    if (name == ')') name = 'lat_line'
  enddo
endif


! If there are multipass lines then expand the lattice and write out
! the post-expand info as needed.

have_expand_lattice_line = .false.
do ie = 0, lat%n_ele_max
  ele => lat%ele(ie)
  !!! if (ele%slave_status == super_slave$) cycle

  if (ele%key == lcavity$ .or. ele%key == rfcavity$) then
    if (ele%value(phi0_multipass$) == 0) cycle
    if (.not. have_expand_lattice_line) call write_expand_lat_header (iu, have_expand_lattice_line)
    !!write (iu, '(3a)') trim(pals_ele_name(ele)), '[phi0_multipass] = '- ele%value(phi0_multipass$))
  endif

enddo

! If there are lattice elements with duplicate names but differing parameters then
! Write the differences.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  do ie = 0, branch%n_ele_max
    ele => branch%ele(ie)
    if (ele%slave_status == super_slave$) cycle
    if (ele%slave_status == multipass_slave$) cycle
    !!! call eles_with_same_name_handler(ele, named_eles_ptr, an_indexx, names, n_names, order)
  enddo
enddo

!---------------------------
! Define lattice

write (iu, '(a)')
name = lat%machine
if (name == '') name = 'machine'
call write_list_header(id, name)
call write_str_dict(id, '-', 'kind', 'Lattice')

id1 = id + id_del
id2 = id + 2*id_del
call write_dict_header(id2, 'branches')

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  if (branch%ix_from_branch > -1) cycle
  name = branch%name
  if (name == '') name = 'beam_line' // int_str(ib)
  call write_list_header(id2+id_del, downcase(name))
  call write_logic_dict(id2+id_del, '-', 'periodic', branch%param%geometry == closed$)
enddo


!-----------------------------------
! Constants

call write_str(0, '')
call write_str(id, '#---------------------------------------------------------------------------------------')
call write_str(id, '# Constants')
call write_str(0, '')

! First print constants used in expressions.

call nametable_init(var_nametab)

header_needed = .true.
do ie = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(ie)
  if (lord%key == group$) then
    print *, 'GROUP ELEMENTS CANNOT YET BE TRANSLATED!'
    cycle
  endif

  if (lord%key == girder$) then
    print *, 'GIRDER ELEMENTS CANNOT YET BE TRANSLATED!'
    cycle
  endif

  if (lord%key == overlay$) then
    do is = 1, lord%n_slave
      slave => pointer_to_slave(lord, is, ctl)
      control = ctl
      if (.not. allocated(control%stack)) then
        print *, 'Overlay', trim(lord%name), ' uses knot points for the control curve. This cannot yet be translated!'
        exit
      endif

      do k = 1, size(control%stack)
        if (control%stack(k)%type /= variable$) cycle
        call find_index(control%stack(k)%name, var_nametab, ix_match, add_to_list = .true., has_been_added = is_added)
        if (is_added) call write_real_list (id, '- constants', control%stack(k)%name, control%stack(k)%value, header_needed)
      enddo
    enddo
  endif
enddo

!-------------------------------
! Overlays and Groups

call write_str(0, '')
call write_str(id, '#---------------------------------------------------------------------------------------')
call write_str(id, '# Bmad Overlay Translation.')
call write_str(id, '# Note: Translated overlay parameter names use a double underscore of the form <Overlay-name>__<param-name>')
call write_str(0, '')

lat%ele%select = .false.

do ie = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(ie)
  if (lord%key == overlay$) then
    call overlay_out(id, lord, lat)
  endif
enddo

! cleanup

close(iu)
deallocate (names, an_indexx)

!----------------------------------------------------------------------------------------------
contains

subroutine write_element_to_beamline (indnt, ele, lat, line)

type (lat_struct), target :: lat
type (ele_struct) :: ele
type (ele_struct), pointer :: lord, m_lord, slave

character(*) line
character(40) lord_name

integer indnt, iu, ix, ii, ix_slave

!

if (ele%slave_status == super_slave$) then
  do ii = 1, slave%n_lord
    lord => pointer_to_lord(ele, ii, ix_slave_back = ix_slave)
    if (lord%lord_status /= super_lord$) cycle
    if (ix_slave /= 1) cycle
    line = trim(line) // ' ' // trim(pals_ele_name(lord)) // ', '
  enddo

elseif (ele%slave_status == multipass_slave$) then
  lord => pointer_to_lord(ele, 1)
  line = trim(line) // ' ' // trim(pals_ele_name(lord)) // ', '

else
  line = trim(line) // ' ' // trim(pals_ele_name(ele)) // ', '
endif

if (len_trim(line) > 80) then
  if (line(1:5) == 'line:') then
    call write_str(indnt, line)
  else
    call write_str(indnt+6, line)
  endif
  line = ''
endif

end subroutine write_element_to_beamline

!--------------------------------------------------------------------------------
! contains

subroutine write_expand_lat_header (iu, have_expand_lattice_line)

integer iu
logical have_expand_lattice_line

write (iu, '(a)')
write (iu, '(a)') '!-------------------------------------------------------'
write (iu, '(a)')
write (iu, '(a)') 'expand_lattice'
write (iu, '(a)')
have_expand_lattice_line = .true.

end subroutine write_expand_lat_header

!--------------------------------------------------------------------------------
! contains

subroutine eles_with_same_name_handler(ele, named_eles_ptr, an_indexx, names, n_names, order)

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0
type (lat_struct), pointer :: lat
type (ele_pointer_struct) :: named_eles_ptr(:)
type (lat_ele_order_struct) order

real(rp), pointer :: a0(:), b0(:), ksl0(:), a(:), b(:), ksl(:)
real(rp), target :: az(0:n_pole_maxx) = 0, bz(0:n_pole_maxx) = 0
character(40), allocatable :: names(:)
integer, allocatable :: an_indexx(:)
integer n_names, ix_match
integer i, iv

!

lat => ele%branch%lat
if (ele%slave_status == multipass_slave$) return
call find_index (ele%name, names, an_indexx, n_names, ix_match)
ele0 => named_eles_ptr(ix_match)%ele   ! Element with this name whose attributes were written to the lattice file.
if (ele%ix_ele == ele0%ix_ele .and. ele%ix_branch == ele0%ix_branch) return

do iv = 1, num_ele_attrib$
  if (ele%value(iv) == ele0%value(iv)) cycle
  info = attribute_info(ele, iv)
  if (info%state /= is_free$ .and. info%state /= quasi_free$) cycle
  if (info%state == quasi_free$) then
    if (.not. attribute_free(ele, info%name, .false.)) cycle
  endif
  ! Have a differing attribute
  call write_this_differing_attrib(iu, ele, attribute_name(ele, iv), ele%value(iv), order)
enddo

if (associated(ele%a_pole) .or. associated(ele0%a_pole)) then
  call pointer_to_ele_multipole(ele0, a0, b0, ksl0, magnetic$)
  if (.not. associated(a0)) a0 => az
  if (.not. associated(b0)) b0 => bz
  call pointer_to_ele_multipole(ele, a, b, ksl, magnetic$)
  if (.not. associated(a)) a => az
  if (.not. associated(b)) b => bz

  if (ele%key == multipole$) then
    do i = 0, n_pole_maxx
      if (a(i) /= a0(i)) call write_this_differing_attrib(iu, ele, 'k' // int_str(i) // 'l', a(i), order)
      if (b(i) /= b0(i)) call write_this_differing_attrib(iu, ele, 't' // int_str(i), b(i), order)
      if (ksl(i) /= ksl0(i)) call write_this_differing_attrib(iu, ele, 'k' // int_str(i) // 'sl', ksl(i), order)
    enddo
  else
    do i = 0, n_pole_maxx
      if (a(i) /= a0(i)) call write_this_differing_attrib(iu, ele, 'a' // int_str(i), a(i), order)
      if (b(i) /= b0(i)) call write_this_differing_attrib(iu, ele, 'b' // int_str(i), b(i), order)
    enddo
  endif
endif

if (associated(ele%b_pole_elec) .or. associated(ele0%b_pole_elec)) then
  call pointer_to_ele_multipole(ele0, a0, b0, ksl0, electric$)
  if (.not. associated(a0)) a0 => az
  if (.not. associated(b0)) b0 => bz
  call pointer_to_ele_multipole(ele, a, b, ksl, electric$)
  if (.not. associated(a)) a => az
  if (.not. associated(b)) b => bz

  do i = 0, n_pole_maxx
    if (a(i) /= a0(i)) call write_this_differing_attrib(iu, ele, 'a' // int_str(i) // '_elec', a(i), order)
    if (b(i) /= b0(i)) call write_this_differing_attrib(iu, ele, 'b' // int_str(i) // '_elec', b(i), order)
  enddo
endif

end subroutine eles_with_same_name_handler

!--------------------------------------------------------------------------------
! contains

subroutine write_this_differing_attrib(iu, ele, attrib_name, value, order)

type (ele_struct) ele
type (lat_ele_order_struct) order

integer iu
real(rp) value
character(*) attrib_name

!

if (.not. have_expand_lattice_line) call write_expand_lat_header (iu, have_expand_lattice_line)

!!write (iu, '(5a)') trim(ele_unique_name(ele, order)), '['- trim(attrib_name), '] = '- value)

end subroutine write_this_differing_attrib

!--------------------------------------------------------------------------------
! contains

function pals_ele_name(ele) result (name_out)

type (ele_struct) ele
character(40) name_out
integer ix

name_out = downcase(ele%name)
if (name_out == 'end') name_out = 'end_b' // int_str(ele%ix_branch)

if (ele%ix_ele == 0) name_out = 'beginning_b' // int_str(ele%ix_branch)

ix = index(name_out, '#')
if (ix /= 0) name_out = name_out(1:ix-1) // '!s' // name_out(ix+1:)

ix = index(name_out, '\')     !'
if (ix /= 0) name_out = name_out(1:ix-1) // '!m' // name_out(ix+1:)

end function pals_ele_name

!--------------------------------------------------------------------------------
! contains

subroutine write_this_taylor(indnt, ele, taylor)

type (ele_struct) ele
type (taylor_struct), target :: taylor(6)
type (taylor_term_struct) term

integer indnt, i, j, k
integer e_max(6)
character(200) line

end subroutine write_this_taylor

!------------------------------------------------------
! contains

recursive subroutine overlay_out(indnt, overlay, lat)

type (lat_struct), target :: lat
type (ele_struct) overlay
type (ele_struct), pointer :: lord, slave
type (control_struct), pointer :: ctl
type (control_struct) control

integer indnt, ix, j, iv, it
logical header_needed
character(1000) c_str(40)

! Output is top down.
! Do not output if overlay is already outputted or has lords that have not yet been outputted.

if (overlay%select) return
do ix = 1, overlay%n_lord
  lord => pointer_to_lord(overlay, ix)
  if (.not. lord%select) return
enddo

! Output vars.

overlay%select = .true.
c_str = ''

do ix = 1, overlay%n_lord
  lord => pointer_to_lord(overlay, ix)
  do j = 1, lord%n_slave
    slave => pointer_to_slave(lord, j, ctl)
    control = ctl
    if (slave%ix_ele /= overlay%ix_ele) cycle
    it = control%ix_attrib - var_offset$
    if (c_str(it) == '') then
      c_str(it) = this_expression(control%stack, lord)
    else
      c_str(it) = trim(c_str(it)) // ' + ' // this_expression(control%stack, lord)
    endif
  enddo
enddo

do iv = 1, size(overlay%control%var)
  name = trim(overlay%name) // '__' // trim(downcase(overlay%control%var(iv)%name))
  header_needed = .true.
  if (c_str(iv) == '') then
    call write_real_dict (indnt, '- constants', name, overlay%control%var(iv)%value, header_needed)
  else
    call write_str_dict (indnt, '- constants', name, c_str(iv), header_needed)
  endif
enddo

! Now that this overlay has been outputted, check if any overlay slaves need outputting.

do ix = 1, overlay%n_slave
  slave => pointer_to_slave(overlay, ix)
  if (slave%key == overlay$) then
    call overlay_out(indnt, slave, lat)
  else
    call overlay_slave_out(indnt, slave, lat)
  endif
enddo

end subroutine overlay_out

!------------------------------------------------------
! contains

recursive subroutine overlay_slave_out(indnt, slave, lat)

type (lat_struct), target :: lat
type (ele_struct) slave
type (ele_struct), pointer :: lord
type (control_struct), pointer :: ctl
type (control_struct)  control

real(rp) f
integer indnt, ix, j, iv, it, n_contl, indx(40), ixm
logical header_needed
character(40) pals_name, attrib_names(40)
character(100) name
character(1000) :: c_str(40)

! Do not output if slave is already outputted or has overlay lords that have not yet been outputted.

if (slave%select) return
do ix = 1, slave%n_lord
  lord => pointer_to_lord(slave, ix)
  if (.not. lord%key == overlay$) cycle
  if (.not. lord%select) return
enddo

! Output slave dependentcies.

slave%select = .true.

c_str = ''
attrib_names = ''
n_contl = 0

do ix = 1, slave%n_lord
  lord => pointer_to_lord(slave, ix, ctl)
  control = ctl
  call find_index(control%attribute, attrib_names, indx, n_contl, ixm)
  if (ixm == 0) call find_index(control%attribute, attrib_names, indx, n_contl, ixm, add_to_list = .true.)

  if (c_str(ixm) == '') then
    c_str(ixm) = '(' // trim(this_expression(control%stack, lord)) // ')'
  else
    c_str(ixm) = trim(c_str(ixm)) // ' + (' // trim(this_expression(control%stack, lord)) // ')'
  endif
enddo

do iv = 1, n_contl
  pals_name = pals_attrib_name(attrib_names(iv), slave, f)
  name = trim(downcase(slave%name)) // '.' // trim(pals_name)
  if (f /= 1.0_rp) c_str(iv) = re_str(f) // ' * ' // trim(c_str(iv))
  header_needed = .true.
  call write_str_list(indnt, '- sets', name, c_str(iv), header_needed)
enddo

end subroutine overlay_slave_out

!------------------------------------------------------
! contains

recursive function this_expression(stack, lord) result(expr)

type (expression_atom_struct) :: stack(:)
type (ele_struct) lord

integer ix_match
character(1000) expr

!

do i = 1, size(stack)
  select case (stack(i)%type)
  case (constant$)      ! Something like "c_light"
    stack(i)%name = upcase(stack(i)%name)
  case (variable$)

  case default
    if (stack(i)%type > var_offset$ .and. stack(i)%type < var_offset$ + n_var_max$) then
      stack(i)%name = trim(lord%name) // '__' // downcase(stack(i)%name)
    endif
  end select
enddo

expr = expression_stack_to_string(stack)

end function this_expression

!------------------------------------------------------
! contains

! Return Pals attribute name given Bmad attribute name.

function pals_attrib_name(bmad_name, ele, factor) result (pals_name)

type (ele_struct) ele

character(*) bmad_name
character(40) pals_name
real(rp) factor

! factor is the conversion factor from Bmad parameter to PALS parameter.
! For example, accounting for a length normalization.

factor = 1.0_rp

if ((bmad_name(1:1) == 'A' .or. bmad_name(1:1) == 'B') .and. is_integer(bmad_name(2:), ix)) then
  if (ele%field_master) then
    pals_name = 'B'
    factor = ele%value(p0c$) / (charge_of(ele%ref_species) * c_light)
  else
    pals_name = 'K'
    factor = 1
  endif

  if (bmad_name(1:1) == 'A') then
    pals_name = pals_name(1:1) // 's'
  else
    pals_name = pals_name(1:1) // 'n'
  endif

  pals_name = trim(pals_name) // bmad_name(2:)
  factor = factor * factorial(ix)

  if (ele%value(l$) == 0) then
    pals_name = trim(pals_name) // 'L'
  else
    factor = factor / ele%value(l$)
  endif

  return
endif

select case (bmad_name)
case ('BL_KICK');       pals_name = 'Bn0L'
case ('BL_HKICK');      pals_name = 'Bn0L'
case ('BL_VKICK');      pals_name = 'Bs0L'
case ('B1_GRADIENT');   pals_name = 'Bn1'
case ('B2_GRADIENT');   pals_name = 'Bn2'
case ('B3_GRADIENT');   pals_name = 'Bn3'
case ('K1');            pals_name = 'Kn1'
case ('K2');            pals_name = 'Kn2'
case ('K3');            pals_name = 'Kn3'
case ('E1');            pals_name = 'e1'
case ('E2');            pals_name = 'e2'
case ('G');             pals_name = 'g_ref'
case ('ANGLE');         pals_name = 'g_ref'
case ('L');             pals_name = 'L'
case ('X_OFFSET', 'Y_OFFSET', 'Z_OFFSET', 'X_PITCH', 'Y_PITCH', 'TILT')
  if (ele%key == patch$) then
    select case (bmad_name)
    case ('X_OFFSET');      pals_name = 'dx'
    case ('Y_OFFSET');      pals_name = 'dy'
    case ('Z_OFFSET');      pals_name = 'dz'
    case ('X_PITCH');       pals_name = 'dy_rot'
    case ('Y_PITCH');       pals_name = 'dx_rot'; factor = -1
    case ('TILT');          pals_name = 'dz_rot'
    end select
  else
    select case (bmad_name)
    case ('X_OFFSET');      pals_name = 'x_offset'
    case ('Y_OFFSET');      pals_name = 'y_offset'
    case ('Z_OFFSET');      pals_name = 'z_offset'
    case ('X_PITCH');       pals_name = 'y_rot'
    case ('Y_PITCH');       pals_name = 'x_rot'; factor = -1
    case ('TILT');          pals_name = 'z_rot'
    end select
  endif

case ('T_OFFSET');      pals_name = 't_offset'
case ('KS');            pals_name = 'Ksol'
case ('BS_FIELD');      pals_name = 'Bsol'
case default
  print *, 'Attribute not yet coded for translation', trim(bmad_name)
  print *, 'Please report this.'
end select

end function pals_attrib_name

!------------------------------------------------------
! contains

subroutine write_list_header(indnt, name)

integer indnt
character(*) name
character(100) :: blank = ''

!

write (iu, '(4a)') blank(1:indnt), '- ', trim(name), ': '

end subroutine write_list_header

!------------------------------------------------------
! contains

subroutine write_list_item(indnt, name)

integer indnt
character(*) name
character(100) :: blank = ''

!

write (iu, '(4a)') blank(1:indnt+id_del), '- ', trim(name)

end subroutine write_list_item

!------------------------------------------------------
! contains

subroutine write_str(indnt, name)

integer indnt
character(*) name

write (iu, '(2a)') blank(1:indnt), trim(name)

end subroutine write_str

!------------------------------------------------------
! contains

subroutine write_str_list(indnt, group_name, name, value, header_needed)

integer indnt, ind
character(*) group_name, name, value
character(100) :: blank = ''
logical, optional :: header_needed

!

if (value == '' .or. value == '""') return
if (logic_option(.false., header_needed) .and. group_name /= '' .and. group_name /= '-') write (iu, '(4a)') blank(1:indnt), trim(group_name), ':'
ind = indent_tot(group_name, indnt+id_del)
write (iu, '(5a)') blank(1:ind), '- ', trim(name), ': ', trim(value)
if (present(header_needed)) header_needed = .false.

end subroutine write_str_list

!------------------------------------------------------
! contains

subroutine write_real_list(indnt, group_name, name, value, header_needed)

real(rp) value, ind
integer indnt
logical, optional :: header_needed
character(*) group_name, name
character(100) :: blank = ''

!

if (value == 0) return
if (logic_option(.false., header_needed) .and. group_name /= '' .and. group_name /= '-') write (iu, '(4a)') blank(1:indnt), trim(group_name), ':'  
ind = indent_tot(group_name, indnt+id_del)
write (iu, '(5a)') blank(1:ind), '- ', trim(name), ': ', re_str(value)
if (present(header_needed)) header_needed = .false.

end subroutine write_real_list

!------------------------------------------------------
! contains

subroutine write_logic_list(indnt, group_name, name, value, header_needed, default)

logical value
logical, optional :: header_needed, default
integer indnt, ind
character(*) group_name, name
character(100) :: blank = ''

if (present(default)) then
  if (value .eqv. default) return
endif

if (logic_option(.false., header_needed) .and. group_name /= '' .and. group_name /= '-') write (iu, '(4a)') blank(1:indnt), trim(group_name), ':'  
ind = indent_tot(group_name, indnt+id_del)
if (value) then
  write (iu, '(4a)') blank(1:ind), '- ', trim(name), ': true'
else
  write (iu, '(4a)') blank(1:ind), '- ', trim(name), ': false'
endif
if (present(header_needed)) header_needed = .false.

end subroutine write_logic_list

!------------------------------------------------------
! contains

subroutine write_switch_list(indnt, group_name, name, default, value_names, value, header_needed)

integer indnt, default, value, ind
logical, optional :: header_needed
character(*) group_name, name, value_names(:)
character(100) :: blank = ''

!

if (value == default) return
if (logic_option(.false., header_needed) .and. group_name /= '' .and. group_name /= '-') write (iu, '(4a)') blank(1:indnt), trim(group_name), ':'  
ind = indent_tot(group_name, indnt+id_del)
write (iu, '(5a)') blank(1:ind), '- ', trim(name), ': ', value_names(value)
if (present(header_needed)) header_needed = .false.

end subroutine write_switch_list

!------------------------------------------------------
! contains

function indent_tot(group_name, indnt) result (indnt_tot)

integer indnt, indnt_tot
character(*) group_name

if (group_name == '') then
  indnt_tot = indnt
elseif (group_name(1:1) == '-') then
  indnt_tot = indnt + id_del
else
  indnt_tot = indnt
endif

end function indent_tot

!------------------------------------------------------
! contains

subroutine write_dict_header(indnt, name)

integer indnt
character(*) name
character(100) :: blank = ''

!

write (iu, '(4a)') blank(1:indnt), trim(name), ': '

end subroutine write_dict_header

!------------------------------------------------------
! contains

subroutine write_str_dict(indnt, group_name, name, value, header_needed)

integer indnt, ind
character(*) group_name, name, value
character(100) :: blank = ''
logical, optional :: header_needed

!

if (value == '' .or. value == '""') return
if (logic_option(.false., header_needed) .and. group_name /= '' .and. group_name /= '-') write (iu, '(4a)') blank(1:indnt), trim(group_name), ':'  
ind = indent_tot(group_name, indnt+id_del)
write (iu, '(5a)') blank(1:ind), trim(name), ': ', trim(value)
if (present(header_needed)) header_needed = .false.

end subroutine write_str_dict

!------------------------------------------------------
! contains

subroutine write_real_dict(indnt, group_name, name, value, header_needed)

real(rp) value
integer indnt, ind
logical, optional :: header_needed
character(*) group_name, name
character(100) :: blank = ''

!

if (value == 0) return
if (logic_option(.false., header_needed) .and. group_name /= '' .and. group_name /= '-') write (iu, '(4a)') blank(1:indnt), trim(group_name), ':'  
ind = indent_tot(group_name, indnt+id_del)
write (iu, '(5a)') blank(1:ind), trim(name), ': ', re_str(value)
if (present(header_needed)) header_needed = .false.

end subroutine write_real_dict

!------------------------------------------------------
! contains

subroutine write_logic_dict(indnt, group_name, name, value, header_needed, default)

logical value
logical, optional :: header_needed, default
integer indnt, ind
character(*) group_name, name
character(100) :: blank = ''

if (present(default)) then
  if (value .eqv. default) return
endif

if (logic_option(.false., header_needed) .and. group_name /= '' .and. group_name /= '-') write (iu, '(4a)') blank(1:indnt), trim(group_name), ':'  
ind = indent_tot(group_name, indnt+id_del)
if (value) then
  write (iu, '(4a)') blank(1:ind), trim(name), ': true'
else
  write (iu, '(4a)') blank(1:ind), trim(name), ': false'
endif
if (present(header_needed)) header_needed = .false.

end subroutine write_logic_dict

!------------------------------------------------------
! contains

subroutine write_switch_dict(indnt, group_name, name, default, value_names, value, header_needed)

integer indnt, default, value, ind
logical, optional :: header_needed
character(*) group_name, name, value_names(:)
character(100) :: blank = ''

!

if (value == default) return
if (logic_option(.false., header_needed) .and. group_name /= '' .and. group_name /= '-') write (iu, '(4a)') blank(1:indnt), trim(group_name), ':'  
ind = indent_tot(group_name, indnt+id_del)
write (iu, '(5a)') blank(1:ind), trim(name), ': ', value_names(value)
if (present(header_needed)) header_needed = .false.

end subroutine write_switch_dict

end subroutine write_lattice_pals_format
