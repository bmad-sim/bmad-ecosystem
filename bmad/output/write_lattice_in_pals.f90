!+
! Subroutine write_lattice_in_pals(pals_file, lat, err_flag)
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

subroutine write_lattice_in_pals(pals_file, lat, err_flag)

use write_lattice_file_mod, dummy => write_lattice_in_pals
use bmad_routine_interface, dummy2 => write_lattice_in_pals
use expression_mod
use taylor_mod, only: mat6_to_taylor

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, lord, slave, slave2, multi_lord
type (coord_struct), pointer :: orb
type (multipass_region_lat_struct), target :: mult_lat
type (multipass_all_info_struct), target :: m_info
type (multipass_region_ele_struct), pointer :: mult_ele(:), m_ele
type (multipass_ele_info_struct), pointer :: e_info
type (ele_pointer_struct), allocatable :: named_eles_ptr(:)  ! List of unique element names 
type (lat_ele_order_struct) order
type (ele_attribute_struct) info
type (taylor_struct) taylor(6), spin_taylor(0:3)
type (nametable_struct) var_nametab, defexpr_nametab
type (control_struct), pointer :: ctl
type (control_struct) control

real(rp) f, length, ang2
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)

integer n, i, j, k, ix, ib, ie, iu, is, n_names, ix_match, ix_pass, ix_r
integer ix_lord, ix_super, ie1, ib1, indent, id
integer, allocatable :: an_indexx(:), index_list(:)

logical header_out, has_been_added, in_multi_region, have_expand_lattice_line, err, is_added
logical has_defexpr_var
logical, optional :: err_flag

character(*) pals_file
character(1) prefix
character(3), parameter :: unit_spin_map(0:3) = ['1.0', '0.0', '0.0', '0.0']
character(100) :: name, look_for, ele_name, blank = ''
character(40), allocatable :: names(:)
character(240) fname
character(1000) line
character(*), parameter :: r_name = 'write_lattice_in_pals'

character(20) :: pals_ele_type(n_key$)

!

pals_ele_type(drift$)                = 'Drift'
pals_ele_type(sbend$)                = 'SBend'
pals_ele_type(quadrupole$)           = 'Quadrupole'
pals_ele_type(group$)                = '??Group'
pals_ele_type(sextupole$)            = 'Sextupole'
pals_ele_type(overlay$)              = '??Overlay'
pals_ele_type(custom$)               = 'Custom'
pals_ele_type(taylor$)               = 'LineElement'
pals_ele_type(rfcavity$)             = 'RFCavity'
pals_ele_type(elseparator$)          = 'ELSeparator'
pals_ele_type(beambeam$)             = 'BeamBeam'
pals_ele_type(wiggler$)              = 'Wiggler'
pals_ele_type(sol_quad$)             = 'Solenoid'
pals_ele_type(marker$)               = 'Marker'
pals_ele_type(kicker$)               = 'Kicker'
pals_ele_type(hybrid$)               = 'Hybrid'
pals_ele_type(octupole$)             = 'Octupole'
pals_ele_type(rbend$)                = 'SBend'
pals_ele_type(multipole$)            = 'Multipole'
pals_ele_type(ab_multipole$)         = 'Multipole'
pals_ele_type(solenoid$)             = 'Solenoid'
pals_ele_type(patch$)                = 'Patch'
pals_ele_type(lcavity$)              = 'RFCavity'
pals_ele_type(null_ele$)             = 'NullEle'
pals_ele_type(beginning_ele$)        = 'BeginningEle'
pals_ele_type(def_line$)             = '!Line'
pals_ele_type(match$)                = 'LineElement'
pals_ele_type(monitor$)              = 'Drift'
pals_ele_type(instrument$)           = 'Drift'
pals_ele_type(hkicker$)              = 'Kicker'
pals_ele_type(vkicker$)              = 'Kicker'
pals_ele_type(rcollimator$)          = 'Drift'
pals_ele_type(ecollimator$)          = 'Drift'
pals_ele_type(girder$)               = 'Girder'
pals_ele_type(converter$)            = 'Converter'
pals_ele_type(photon_fork$)          = 'Fork'
pals_ele_type(fork$)                 = 'Fork'
pals_ele_type(mirror$)               = 'Mirror'
pals_ele_type(crystal$)              = 'Crystal'
pals_ele_type(pipe$)                 = 'Drift'
pals_ele_type(capillary$)            = 'Capillary'
pals_ele_type(multilayer_mirror$)    = 'MultilayerMirror'
pals_ele_type(e_gun$)                = 'EGun'
pals_ele_type(em_field$)             = 'EMField'
pals_ele_type(floor_shift$)          = 'FloorShift'
pals_ele_type(fiducial$)             = 'Fiducial'
pals_ele_type(undulator$)            = 'Undulator'
pals_ele_type(diffraction_plate$)    = 'DiffractionPlate'
pals_ele_type(photon_init$)          = 'PhotonInit'
pals_ele_type(sample$)               = 'Sample'
pals_ele_type(detector$)             = 'Detector'
pals_ele_type(sad_mult$)             = 'SadMult'
pals_ele_type(mask$)                 = 'Mask'
pals_ele_type(ac_kicker$)            = 'ACKicker'
pals_ele_type(lens$)                 = 'Lens'
pals_ele_type(crab_cavity$)          = 'CrabCavity'
pals_ele_type(ramper$)               = 'Ramper'
pals_ele_type(def_ptc_com$)          = '!PTC_Com'
pals_ele_type(rf_bend$)              = 'RFBend'
pals_ele_type(gkicker$)              = 'Kicker'
pals_ele_type(foil$)                 = 'Foil'
pals_ele_type(thick_multipole$)      = 'ThickMultipole'
pals_ele_type(pickup$)               = 'Drift'
pals_ele_type(feedback$)             = 'Drift'
pals_ele_type(fixer$)                = 'Fixer'

! Open file

call fullfilename(pals_file, fname)
iu = lunget()
open (iu, file = fname, status = 'unknown')

! Header

write (iu, '(a)')  'PALS:'
write (iu, '(a)')  '  notes:'
write (iu, '(a)')  '    - "File generated by: write_lattice_in_foreign_format."'
write (iu, '(4a)') '    - "From Bmad lattice file: ', trim(lat%input_file_name), '"'
write (iu, '(a)')  ''
write (iu, '(a)')  '  extension_names:'
write (iu, '(a)')  '    - names:'
write (iu, '(a)')  '        SciBmad'
write (iu, '(a)')  '    - Prefixes:'
write (iu, '(a)')  '        SciBmad_'
write (iu, '(a)')  ''

indent = 4

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

    !!! multi_lord => pointer_to_multipass_lord (ele, ix_pass) 
    !!! if (ele%lord_status == super_lord$ .and. ix_pass > 0) cycle
    !!! if (ele%slave_status == super_slave$ .and. ix_pass > 1) cycle

    !!!if (ele%slave_status == super_slave$) then
    !!!  lord => pointer_to_lord(ele, 1)
    !!!  slave => pointer_to_slave(lord, 1)
    !!!  slave2 => pointer_to_slave(lord, lord%n_slave)
    !!!  write (iu, '(2(a, i0), 2a)') '  slave_drift_', ib, '_', ele%ix_ele, ' = Drift(L = ', re_str(length) // ')'
    !!!  cycle
    !!!endif

    !!! if (ix_pass > 0) cycle

    ! Do not write anything for elements that have a duplicate name.

    call add_this_name_to_list (ele, names, an_indexx, n_names, ix_match, has_been_added, named_eles_ptr)
    if (.not. has_been_added) cycle
    id = indent + 2

    ele_name = pals_ele_name(ele)
    if (ie == 0) ele_name = 'begin' // int_str(ib+1)
    call write_param_group(id, ele_name)

    call write_param(id, '', 'kind', pals_ele_type(ele%key), .true.)
    call write_logic_param(id, '', 'is_on', .true., ele%is_on, .true.)


    call write_param(id, '', 'label', quote(ele%type), header_out)
    call write_param(id, '', 'alias', quote(ele%alias), header_out)
    if (associated(ele%descrip)) call write_param(id, '', 'description', quote(ele%descrip), .true.)

    if (has_attribute(ele, 'L') .and. length /= 0) call write_real_param(id, '', 'L', length, .true.)

    if (ie == 0) then
      header_out = .false.
      call write_real_param(id, 'ReferenceP', 'pc_ref', ele%value(p0c$), header_out)
      call write_param(id, 'ReferenceP', 'species_ref: ', quote(openpmd_species_name(ele%ref_species)), header_out)

      header_out = .false.
      write (iu, '(2a)') blank(1:id), '- TwissP:'
      call write_real_param(id, 'TwissP', 'beta_a',  ele%a%beta, header_out)
      call write_real_param(id, 'TwissP', 'beta_b',  ele%b%beta, header_out)
      call write_real_param(id, 'TwissP', 'alpha_a', ele%a%alpha, header_out)
      call write_real_param(id, 'TwissP', 'alpha_b', ele%b%alpha, header_out)
      call write_real_param(id, 'TwissP', 'eta_x',   ele%x%eta, header_out)
      call write_real_param(id, 'TwissP', 'eta_y',   ele%y%eta, header_out)
      call write_real_param(id, 'TwissP', 'etap_x',  ele%x%etap, header_out)
      call write_real_param(id, 'TwissP', 'etap_y',  ele%y%etap, header_out)
      call write_real_param(id, 'TwissP', 'c_mat11', ele%c_mat(1,1), header_out)
      call write_real_param(id, 'TwissP', 'c_mat12', ele%c_mat(1,2), header_out)
      call write_real_param(id, 'TwissP', 'c_mat21', ele%c_mat(2,1), header_out)
      call write_real_param(id, 'TwissP', 'c_mat22', ele%c_mat(2,2), header_out)

      orb => lat%particle_start
      header_out = .false.
      call write_real_param(id, 'ParticleP', 'x',  orb%vec(1), header_out)
      call write_real_param(id, 'ParticleP', 'px', orb%vec(2), header_out)
      call write_real_param(id, 'ParticleP', 'y',  orb%vec(3), header_out)
      call write_real_param(id, 'ParticleP', 'py', orb%vec(4), header_out)
      call write_real_param(id, 'ParticleP', 'z',  orb%vec(5), header_out)
      call write_real_param(id, 'ParticleP', 'pz', orb%vec(6), header_out)
      call write_real_param(id, 'ParticleP', 'spin_x', orb%spin(1), header_out)
      call write_real_param(id, 'ParticleP', 'spin_y', orb%spin(2), header_out)
      call write_real_param(id, 'ParticleP', 'spin_z', orb%spin(3), header_out)
    endif

    !

    if (ele%key == sbend$) then
      header_out = .false.
      call write_real_param(id, 'ParticleP', 'g_ref',     ele%value(g$), header_out)
      call write_real_param(id, 'ParticleP', 'tilt_ref',  ele%value(ref_tilt$), header_out)
      call write_real_param(id, 'ParticleP', 'e1',        ele%value(e1$), header_out)
      call write_real_param(id, 'ParticleP', 'e2',        ele%value(e2$), header_out)
      call write_real_param(id, 'ParticleP', 'h1',        ele%value(h1$), header_out)
      call write_real_param(id, 'ParticleP', 'h2',        ele%value(h2$), header_out)
      call write_real_param(id, 'ParticleP', 'edge_int1', ele%value(fint$)*ele%value(hgap$), header_out)
      call write_real_param(id, 'ParticleP', 'edge_int2', ele%value(fintx$)*ele%value(hgapx$), header_out)
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
    header_out = .false.

    do j = 0, ix
      if (length == 0) then
        call write_real_param(id, 'MagneticMultipoleP', prefix // 's' // int_str(j) // 'L', f * factorial(j) * a_pole(j), header_out)
        call write_real_param(id, 'MagneticMultipoleP', prefix // 'n' // int_str(j) // 'L', f * factorial(j) * b_pole(j), header_out)
      else
        call write_real_param(id, 'MagneticMultipoleP', prefix // 's' // int_str(j), f * factorial(j) * a_pole(j), header_out)
        call write_real_param(id, 'MagneticMultipoleP', prefix // 'n' // int_str(j), f * factorial(j) * b_pole(j), header_out)
      endif
    enddo

    ! Electric multipoles

    call multipole_ele_to_ab(ele, .false., ix, a_pole, b_pole, electric$, include_kicks$)
    header_out = .false.

    do j = 0, ix
      call write_real_param(id, 'ElectricMultipoleP', 'Es' // int_str(j), factorial(j) * a_pole(j), header_out)
      call write_real_param(id, 'ElectricMultipoleP', 'En' // int_str(j), factorial(j) * b_pole(j), header_out)
    enddo

    !

    if (has_attribute(ele, 'X1_LIMIT')) then
      header_out = .false.
      call write_real_param(id,   'ApertureP', 'x_min', -ele%value(x1_limit$), header_out)
      call write_real_param(id,   'ApertureP', 'x_max',  ele%value(x2_limit$), header_out)
      call write_real_param(id,   'ApertureP', 'y_min', -ele%value(y1_limit$), header_out)
      call write_real_param(id,   'ApertureP', 'y_max',  ele%value(y2_limit$), header_out)
      call write_switch_param(id, 'ApertureP', 'shape', 0, [character(16):: 'Auto', 'RECTANGULAR', 'ELLIPTICAL'], ele%aperture_type, header_out)
      call write_switch_param(id, 'ApertureP', 'location', 1, [character(16):: 'ENTRANCE_END', 'EXIT_END', 'BOTH_ENDS', 'NOWHERE', 'EVERYWHERE'], ele%aperture_at, header_out)
      call write_logic_param(id,  'ApertureP', 'aperture_shifts_with_body', .false., ele%offset_moves_aperture, header_out)
    endif

    !

    select case (ele%key)
    case (match$, taylor$)
      length = ele%value(l$)

      select case (ele%key)
      case (match$)
        call mat6_to_taylor(ele%vec0, ele%mat6, taylor)
        call write_this_taylor(indent, ele, taylor)
        cycle

      case (taylor$)
        call write_this_taylor(indent, ele, ele%taylor)
        cycle
      end select
    end select

    !

    if (ele%key == patch$) then
      header_out = .false.
      call write_real_param(id,   'PatchP', 'x_offset', ele%value(x_offset$), header_out)
      call write_real_param(id,   'PatchP', 'y_offset', ele%value(y_offset$), header_out)
      call write_real_param(id,   'PatchP', 'z_offset', ele%value(z_offset$), header_out)
      call write_real_param(id,   'PatchP', 'x_rot',   -ele%value(y_pitch$), header_out)
      call write_real_param(id,   'PatchP', 'y_rot',    ele%value(x_pitch$), header_out)
      call write_real_param(id,   'PatchP', 'z_rot',    ele%value(tilt$), header_out)
      call write_logic_param(id,  'PatchP', 'flexible', .false., is_true(ele%value(flexible$)), header_out)
      call write_logic_param(id,  'PatchP', 'user_sets_length', .false., is_true(ele%value(user_sets_length$)), header_out)
      call write_switch_param(id, 'PatchP', 'ref_coords', 2, [character(16):: 'entrance_end ', 'exit_end'], nint(ele%value(ref_coords$)), header_out)


      header_out = .false.
      call write_real_param(id, 'ReferenceChangeP', 'dE_ref',    ele%value(E_tot_offset$), header_out)
      call write_real_param(id, 'ReferenceChangeP', 'E_tot_ref', ele%value(E_tot_set$), header_out)
      call write_real_param(id, 'ReferenceChangeP', 't_offset',  ele%value(t_offset$), header_out)


    elseif (has_attribute(ele, 'X_PITCH')) then
      header_out = .false.
      call write_real_param(id, 'BodyShiftP', 'x_offset', ele%value(x_offset$), header_out)
      call write_real_param(id, 'BodyShiftP', 'y_offset', ele%value(y_offset$), header_out)
      call write_real_param(id, 'BodyShiftP', 'z_offset', ele%value(z_offset$), header_out)
      call write_real_param(id, 'BodyShiftP', 'x_rot',   -ele%value(y_pitch$), header_out)
      call write_real_param(id, 'BodyShiftP', 'y_rot',    ele%value(x_pitch$), header_out)
      call write_real_param(id, 'BodyShiftP', 'z_rot',    ele%value(tilt$), header_out)
    endif

    !

    if (has_attribute(ele, 'KS')) then
      header_out = .false.
      if (ele%field_master) then
        call write_real_param(id, 'SolenoidP', 'Bsol', ele%value(bs_field$), header_out)
      else
        call write_real_param(id, 'SolenoidP', 'Ksol', ele%value(ks$), header_out)
      endif
    endif

    !

    if (ele%key == lcavity$) then
      call write_param_group(id, 'TrackingP')
      call write_param_group(id+2, 'SciBmad')
      call write_param(id+2, 'SciBmad', 'tracking_method', 'SaganCavity', .true.)
    endif


    if (has_attribute(ele, 'RF_FREQUENCY')) then
      header_out = .false.

      if (is_true(ele%value(harmon_master$))) then
        call write_real_param(id, 'RFP', 'rf_frequency', ele%value(rf_frequency$), header_out)
      else
        call write_real_param(id, 'RFP', 'rf_frequency', ele%value(rf_frequency$), header_out)
      endif

      if (ele%key == lcavity$) then
        if (ele%field_master) then
          call write_real_param(id, 'RFP', 'gradient', ele%value(gradient$) + ele%value(gradient_err$), header_out)
        else
          call write_real_param(id, 'RFP', 'voltage', ele%value(voltage$) + ele%value(voltage_err$), header_out)
        endif

        call write_real_param(id, 'RFP', 'phi0', ele%value(phi0$) + ele%value(phi0_err$), header_out)
        call write_real_param(id, 'RFP', 'num_cells', ele%value(n_rf_steps$), header_out)
        call write_real_param(id, 'RFP', 'SciBmad_L_active', ele%value(l_active$), header_out)

      else
        if (ele%field_master) then
          call write_real_param(id, 'RFP', 'gradient', ele%value(gradient$), header_out)
        else
          call write_real_param(id, 'RFP', 'voltage', ele%value(voltage$), header_out)
        endif

        call write_real_param(id, 'RFP', 'phi0', ele%value(phi0$), header_out)
        call write_param(id, 'RFP', 'zero_phase', 'ABOVE_TRANSITION', header_out)
      endif


      if (has_attribute(ele, 'CAVITY_TYPE')) then
        call write_switch_param(id, 'RFP', 'cavity_type', 0, [character(20):: 'standing_wave', 'traveling_wave', 'standing_wave'], nint(ele%value(cavity_type$)), header_out)
      endif
    endif

  enddo ele_loop
enddo

!------------------------------------------------------------------------------------------------------
! Write branch lines
! First write multipass lines

call multipass_region_info(lat, mult_lat, m_info)

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
      call write_param_group(id, 'BeamLine')
      call write_param(id, 'BeamLine', 'name', 'multi_line_' // int_str(ix_r), .true.)
      call write_param(id, 'BeamLine', 'multipass', 'true', .true.)
    endif

    if (mult_ele(ie)%ix_region /= ix_r) then
      call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #2! PLEASE REPORT THIS!')
    endif

    call write_pals_element (line, iu, ele, lat)

    if (mult_ele(ie)%region_stop_pt) in_multi_region = .false.
  enddo

  if (in_multi_region) then
    call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #3! PLEASE REPORT THIS!')
  endif
enddo  ! ib branch loop


!------------------------------
! Overlay and group elements....

write (iu, '(a)') '#---------------------------------------------------------------------------------------'
write (iu, '(a)') '# Overlay and Group elements'
write (iu, '(a)')

! First print constants used in expressions.

call nametable_init(var_nametab)
call nametable_init(defexpr_nametab)

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
        if (is_added) write (iu, '(2a, es24.17)') trim(control%stack(k)%name), ' = '- control%stack(k)%value
      enddo
    enddo
  endif
enddo

! 

lat%ele%select = .false.

do ie = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(ie)
  if (lord%key == overlay$) then
    call overlay_out(lord, lat, defexpr_nametab)
  endif
enddo

!------------------------------
! Lines for all the branches.
! If we get into a multipass region then name in the main_line list is "multi_line_nn".
! But only write this once.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  write (iu, '(a)')
  name = downcase(branch%name)
  if (name == ')') name = 'lat_line'
  !!line = trim(name), ' = Beamline(['     !, quote(name), '- [', trim(pals_ele_name(branch%ele(0))), ','

  in_multi_region = .false.
  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)

    e_info => m_info%branch(ib)%ele(ie)

    if (.not. e_info%multipass) then
      call write_pals_element (line, iu, ele, lat)
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
        call write_param_group(id, 'multi_line_' // int_str(ix_r))
      if (m_ele%region_start_pt) then
        look_for = 'stop'
      else
        call write_param(id+2, '', 'direction', '-1', .true.)
        look_for = 'start'
      endif
    endif

    if (look_for == 'start' .and. m_ele%region_start_pt .or. &
        look_for == 'stop' .and. m_ele%region_stop_pt) then 
      in_multi_region = .false.
    endif
  enddo
enddo

! Define lat

if (.false.) then
  line = 'lat = expand(' // quote(downcase(lat%use_name)) // '- ['
  do ib = 0, ubound(lat%branch, 1)
    branch => lat%branch(ib)
    if (branch%ix_from_branch > -1) cycle
    name = downcase(branch%name)
    if (name == ')') name = 'lat_line'
  enddo

  ix = index(line, '[, ')
  line = line(:ix) // trim(line(ix+3:)) // '])'
  write (iu, '(a)')
  write (iu, '(a)') trim(line)
endif

! If there are multipass lines then expand the lattice and write out
! the post-expand info as needed.

have_expand_lattice_line = .false.
do ie = 1, lat%n_ele_max
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

  do ie = 1, branch%n_ele_max
    ele => branch%ele(ie)
    if (ele%slave_status == super_slave$) cycle
    if (ele%slave_status == multipass_slave$) cycle
    !!! call eles_with_same_name_handler(ele, named_eles_ptr, an_indexx, names, n_names, order)
  enddo
enddo

!---------------------------
! Define lattice

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  if (branch%ix_from_branch > -1) cycle
  name = branch%name
  if (name == '') name = 'lat_line'
  line = trim(line) // ', ' // name
enddo


ix = index(lat%use_name, ',')
if (ix == 0) ix = len_trim(lat%use_name)+1

call write_param_group(id, quote(lat%use_name(1:ix-1)))

if (ele%key == fork$ .or. ele%key == photon_fork$) then
  header_out = .false.
  !! call write_param(id, '', 'to_line', )

  n = nint(ele%value(ix_to_branch$))
  call write_param(id, '', 'to_line', trim(downcase(lat%branch(n)%name)), .true.)
  if (ele%value(ix_to_element$) > 0) then
    i = nint(ele%value(ix_to_element$))
    call write_param(id, '', 'to_element', trim(pals_ele_name(lat%branch(n)%ele(i))), .true.)
  endif
endif



! cleanup

close(iu)
deallocate (names, an_indexx)
!!! deallocate (mult_lat%branch)

!----------------------------------------------------------------------------------------------
contains

function jbool(logic) result (bool_str)
logical logic
character(5) bool_str

if (logic) then
  bool_str = 'true'
else
  bool_str = 'false'
endif

end function jbool

!----------------------------------------------------------------------------------------------
! contains

subroutine write_pals_element (line, iu, ele, lat)

type (lat_struct), target :: lat
type (ele_struct) :: ele
type (ele_struct), pointer :: lord, m_lord, slave

character(*) line
character(40) lord_name

integer iu, ix

!

!!!if (ele%slave_status == super_slave$) then
!!!  if (ele%orientation == 1) then
!!!    write (line, '(a, 2(a, i0), a)') ' slave_drift_'- ele%ix_branch, '_'- ele%ix_ele, ','
!!!  else
!!!    write (line, '(a, 2(a, i0), a)') ' reverse(slave_drift_'- ele%ix_branch, '_'- ele%ix_ele, '),'
!!!  endif
!!!
!!!elseif (ele%slave_status == multipass_slave$) then
!!!  lord => pointer_to_lord(ele, 1)
!!!  write (line, '(4a)') ' '- trim(downcase(lord%name)), ','
!!!
!!!else
  if (ele%orientation == 1) then
    !! write (line, '(4a)') ' '- trim(pals_ele_name(ele)), ','
  else
    !! write (line, '(4a)') ' reverse('- trim(pals_ele_name(ele)), '),'
  endif
!!!endif

end subroutine write_pals_element

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

recursive subroutine overlay_out(overlay, lat, defexpr_nametab)

type (lat_struct), target :: lat
type (ele_struct) overlay
type (ele_struct), pointer :: lord, slave
type (control_struct), pointer :: ctl
type (control_struct) control
type (nametable_struct) defexpr_nametab

integer ix, j, iv, it

character(1000) c_str(40)

! Output is top down.
! Do not output if overlay is already outputted or has lords that have not yet been outputted.

if (overlay%select) return
do ix = 1, overlay%n_lord
  lord => pointer_to_lord(overlay, ix)
  if (.not. lord%select) return
enddo

! Output vars.
! Controled vars are defined with a defered expression.

overlay%select = .true.
c_str = ''
has_defexpr_var = .false.

do ix = 1, overlay%n_lord
  lord => pointer_to_lord(overlay, ix)
  do j = 1, lord%n_slave
    slave => pointer_to_slave(lord, j, ctl)
    control = ctl
    if (slave%ix_ele /= overlay%ix_ele) cycle
    it = control%ix_attrib - var_offset$
    if (c_str(it) == '') then
      c_str(it) = this_expression(control%stack, lord, has_defexpr_var)
    else
      c_str(it) = trim(c_str(it)) // ' + ' // this_expression(control%stack, lord, has_defexpr_var)
    endif
  enddo
enddo

do iv = 1, size(overlay%control%var)
  name = trim(overlay%name) // '_' // trim(downcase(overlay%control%var(iv)%name))
  if (c_str(iv) == ')') then
    write (iu, '(2a, es24.16)') trim(name), ' = '- overlay%control%var(iv)%value
  elseif (has_defexpr_var) then
    write (iu, '(3a)') 'if !@isdefined(', trim(name), ')'
    write (iu, '(7a)') '  const ', trim(name), ' = ', trim(c_str(iv))
    write (iu, '(a)')  'end'
    call nametable_add(defexpr_nametab, name, 1)
  else
    write (iu, '(3a)') 'if !@isdefined(', trim(name), ')'
    write (iu, '(7a)') '  const ', trim(name), ' = DefExpr((3a) -> ', trim(c_str(iv)), ')'
    write (iu, '(a)')  'end'
    call nametable_add(defexpr_nametab, name, 1)
  endif
enddo

! Now that this overlay has been outputted, check if any overlay slaves need outputting.

do ix = 1, overlay%n_slave
  slave => pointer_to_slave(overlay, ix)
  if (slave%key == overlay$) then
    call overlay_out(slave, lat, defexpr_nametab)
  else
    call overlay_slave_out(slave, lat, defexpr_nametab)
  endif
enddo

end subroutine overlay_out

!------------------------------------------------------
! contains

recursive subroutine overlay_slave_out(slave, lat, defexpr_nametab)

type (lat_struct), target :: lat
type (ele_struct) slave
type (ele_struct), pointer :: lord
type (control_struct), pointer :: ctl
type (control_struct)  control
type (nametable_struct) defexpr_nametab

real(rp) f
integer ix, j, iv, it, n_contl, indx(40), ixm

character(40) pals_name, attrib_names(40)
character(100) name
character(1000) :: c_str(40)
logical has_defexpr_var(40)

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
has_defexpr_var = .false.

do ix = 1, slave%n_lord
  lord => pointer_to_lord(slave, ix, ctl)
  control = ctl
  call find_index(control%attribute, attrib_names, indx, n_contl, ixm)
  if (ixm == 0) call find_index(control%attribute, attrib_names, indx, n_contl, ixm, add_to_list = .true.)

  if (c_str(ixm) == ')') then
    c_str(ixm) = this_expression(control%stack, lord, has_defexpr_var(ixm))
  else
    c_str(ixm) = trim(c_str(ixm)) // ' + ' // this_expression(control%stack, lord, has_defexpr_var(ixm))
  endif
enddo

do iv = 1, n_contl
  pals_name = pals_attrib_name(attrib_names(iv), slave, f)
  name = trim(downcase(slave%name)) // '.' // trim(pals_name)
  if (f /= 1.0_rp) c_str(iv) = re_str(f) // ' * (' // trim(c_str(iv)) // ')'

  if (has_defexpr_var(iv)) then
    write (iu, '(4a)') trim(name), ' = ', trim(c_str(iv))
  else
    write (iu, '(6a)') trim(name), ' = DefExpr((3a) -> ', trim(c_str(iv)), ')'
  endif

enddo

end subroutine overlay_slave_out

!------------------------------------------------------
! contains

recursive function this_expression(stack, lord, has_defexpr_var) result(expr)

type (expression_atom_struct) :: stack(:)
type (ele_struct) lord

integer ix_match
character(1000) expr
logical has_defexpr_var

!

do i = 1, size(stack)
  select case (stack(i)%type)
  case (constant$)      ! Something like "c_light"
    stack(i)%name = upcase(stack(i)%name)
  case (variable$)

  case default
    if (stack(i)%type > var_offset$ .and. stack(i)%type < var_offset$ + n_var_max$) then
      stack(i)%name = trim(lord%name) // '_' // downcase(stack(i)%name)
      call find_index(stack(i)%name, defexpr_nametab, ix_match)
      if (ix_match >0) has_defexpr_var = .true.
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

!

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

subroutine write_param(indnt, group_name, name, value, header_out)

integer indnt
character(*) group_name, name, value
logical header_out

end subroutine write_param

!------------------------------------------------------
! contains

subroutine write_param_group(indnt, name)

integer indnt
character(*) name

end subroutine write_param_group

!------------------------------------------------------
! contains

subroutine write_real_param(indnt, group_name, name, value, header_out)

real(rp) value
integer indnt
logical header_out
character(*) group_name, name


end subroutine write_real_param

!------------------------------------------------------
! contains

subroutine write_logic_param(indnt, group_name, name, default, value, header_out)

logical default, value, header_out
integer indnt
character(*) group_name, name


end subroutine write_logic_param

!------------------------------------------------------
! contains

subroutine write_switch_param(indnt, group_name, name, default, value_names, value, header_out)


integer indnt, default, value
logical header_out
character(*) group_name, name, value_names(:)


end subroutine write_switch_param

end subroutine write_lattice_in_pals
