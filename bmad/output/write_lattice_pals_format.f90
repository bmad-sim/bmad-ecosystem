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
type (ele_pointer_struct), allocatable :: ctl_slave(:)       ! Slaves of a controller being written
type (lat_ele_order_struct) order
type (ele_attribute_struct) info
type (taylor_struct) taylor(6), spin_taylor(0:3)
type (nametable_struct) var_nametab
type (control_struct), pointer :: ctl
type (control_struct) control
type (control_var1_struct), pointer :: var

real(rp) f, length, ang2, resid, factor(2)
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)

integer n_pals, n_ctl, ixm
integer, allocatable :: ctl_indx(:)
integer n, i, j, k, ix, ib, ie, iu, is, n_names, ix_match, ix_pass, ix_r
integer :: ix_lord, ix_super, ie1, ib1, id, id1, id2, id_del = 2
integer, allocatable :: an_indexx(:), index_list(:)

logical header_needed, has_been_added, in_multi_region, have_expand_lattice_line, err, is_added
logical, optional :: err_flag

character(*) pals_file
character(1) prefix
character(3), parameter :: unit_spin_map(0:3) = ['1.0', '0.0', '0.0', '0.0']
character(40) pals_name(2)
character(40), allocatable :: ctl_attrib(:)
character(100), allocatable :: ctl_param(:)
character(1000), allocatable :: ctl_expr(:)
character(100) :: name, look_for, ele_name, attrib_name, blank = ''
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
pals_ele_type(group$)                = 'Controller'
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
pals_ele_type(overlay$)              = 'Controller'
pals_ele_type(patch$)                = 'Patch'
pals_ele_type(photon_fork$)          = 'Fork'
pals_ele_type(photon_init$)          = 'PhotonInit'
pals_ele_type(pickup$)               = 'Instrument'
pals_ele_type(pipe$)                 = 'Drift'
pals_ele_type(quadrupole$)           = 'Quadrupole'
pals_ele_type(ramper$)               = 'Ramper'
pals_ele_type(rbend$)                = 'Bend'
pals_ele_type(rcollimator$)          = 'Drift'
pals_ele_type(rf_bend$)              = 'RFBend'
pals_ele_type(rfcavity$)             = 'RFCavity'
pals_ele_type(sad_mult$)             = 'SadMult'
pals_ele_type(sample$)               = 'Sample'
pals_ele_type(sbend$)                = 'Bend'
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
write (iu, '(a)')  '  extension_labels:'
write (iu, '(a)')  '    names:'
write (iu, '(a)')  '      Bmad: Bmad specific data'
write (iu, '(a)')  '    prefixes:'
write (iu, '(a)')  '      Bmad_: Bmad specific data'
write (iu, '(a)')  ''

!

id = 2 * id_del
id1 = id + id_del
id2 = id + 2*id_del

call write_str(0, '')
call write_str(id_del, '#---------------------------------------------------------------------------------------')
write (iu, '(a)')  '  facility:'

!

call write_str(0, '')
call write_str(id, '#---------------------------------------------------------------------------------------')
call write_str(id, '# Constants')

if (allocated(lat%constant)) then
  header_needed = .true.
  call write_str(0, '')
  do i = 1, size(lat%constant)
    call write_real_list(id, '- constants', downcase(lat%constant(i)%name), lat%constant(i)%value, &
                                                                                header_needed, .false.)
  enddo
endif

!----------------------------
! Write element defs

! Note: Beamlines cannot currently handle multipass nor superimpose so ignore.
! Stuff that is commented out due to this is marked by "!!!"

call write_str(0, '')
call write_str(id, '#---------------------------------------------------------------------------------------')
call write_str(id, '# Lattice Elements')

n_names = 0
n = lat%n_ele_max
allocate (names(n), an_indexx(n), named_eles_ptr(n))

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  ele_loop: do ie = 0, branch%n_ele_track   !!! Note: Not n_ele_max since superimpose not handled
    ele => branch%ele(ie)
    length = ele%value(l$)
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

    header_needed = .true.
    call write_str_dict(id2, 'MetaP', 'label', quote(ele%type), header_needed)
    call write_str_dict(id2, 'MetaP', 'alias', quote(ele%alias), header_needed)
    if (associated(ele%descrip)) call write_str_dict(id2, 'MetaP', 'description', quote(ele%descrip), header_needed)

    if (has_attribute(ele, 'L') .and. length /= 0) call write_real_dict(id1, '', 'length', length)

    if (ie == 0) then
      header_needed = .true.
      call write_real_dict(id2, 'ReferenceP', 'pc_ref', ele%value(p0c$), header_needed)
      call write_str_dict(id2, 'ReferenceP', 'species_ref', quote(openpmd_species_name(ele%ref_species)), header_needed)

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
      call write_real_dict(id2, 'TwissP', 'dbeta_a_dpz',  ele%a%dbeta_dpz, header_needed)
      call write_real_dict(id2, 'TwissP', 'dbeta_b_dpz',  ele%b%dbeta_dpz, header_needed)
      call write_real_dict(id2, 'TwissP', 'dalpha_a_dpz', ele%a%dalpha_dpz, header_needed)
      call write_real_dict(id2, 'TwissP', 'dalpha_b_dpz', ele%b%dalpha_dpz, header_needed)
      call write_real_dict(id2, 'TwissP', 'deta_x_dpz',   ele%x%deta_dpz, header_needed)
      call write_real_dict(id2, 'TwissP', 'deta_y_dpz',   ele%y%deta_dpz, header_needed)
      call write_real_dict(id2, 'TwissP', 'detap_x_dpz',  ele%x%detap_dpz, header_needed)
      call write_real_dict(id2, 'TwissP', 'detap_y_dpz',  ele%y%detap_dpz, header_needed)

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
      call write_real_dict(id2, 'BendP', 'edge1_int', ele%value(fint$)*ele%value(hgap$), header_needed)
      call write_real_dict(id2, 'BendP', 'edge2_int', ele%value(fintx$)*ele%value(hgapx$), header_needed)
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
      if (.not. header_needed) then
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
    if (key_name(ele%key) /= pals_ele_type(ele%key)) call write_str_dict(id2, 'Bmad', 'Bmad_key', key_name(ele%key), header_needed)
    if (has_attribute(ele, 'TRACKING_METHOD') .and. (ele%tracking_method /= ele_dflt%tracking_method)) &
                                call write_str_dict(id2, 'Bmad', 'tracking_method', tracking_method_name(ele%key), header_needed)
    if (has_attribute(ele, 'MAT6_CALC_METHOD') .and. (ele%mat6_calc_method /= ele_dflt%mat6_calc_method)) &
                                call write_str_dict(id2, 'Bmad', 'mat6_calc_method', mat6_calc_method_name(ele%key), header_needed)

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
call write_str(id, '# Lord elements')
call write_str(0, '')

! First print constants used in expressions.

call nametable_init(var_nametab)

header_needed = .true.
do ie = lat%n_ele_track+1, lat%n_ele_max
  ele => lat%ele(ie)
  select case (ele%key)

  case (girder$)
    print *, 'GIRDER ELEMENTS CANNOT YET BE TRANSLATED!'
    cycle

  case (overlay$, group$)
    ele_name = pals_ele_name(ele)
    write (iu, '(a)')
    call write_list_header(id, ele_name)
    call write_str_dict(id1, '', 'kind', 'Controller')
    if (ele%key == overlay$) then
      call write_str_dict(id1, '', 'control_type', 'ABSOLUTE')
    else
      call write_str_dict(id1, '', 'control_type', 'RELATIVE')
    endif

    header_needed = .true.
    do i = 1, size(ele%control%var)
      var => ele%control%var(i)
      call write_real_dict(id2, 'variables', downcase(var%name), var%value, header_needed, .false.)
    enddo

    ! Note: A Bmad attribute may map to more than one PALS attribute (EG: HKICK of a tilted element
    ! maps to both the normal and skew n = 0 multipole components) and multiple Bmad attributes may
    ! map to a single PALS attribute. So the expressions are collected by PALS parameter name.

    n_ctl = 0
    if (allocated(ctl_param)) deallocate (ctl_param, ctl_attrib, ctl_expr, ctl_indx, ctl_slave)
    n = max(1, 2*ele%n_slave)
    allocate (ctl_param(n), ctl_attrib(n), ctl_expr(n), ctl_indx(n), ctl_slave(n))
    ctl_param = ''; ctl_expr = ''

    do is = 1, ele%n_slave
      slave => pointer_to_slave(ele, is, ctl)
      control = ctl
      if (.not. allocated(control%stack)) then
        print *, 'Overlay', trim(ele%name), ' uses knot points for the control curve. This cannot yet be translated!'
        exit
      endif

      call pals_attrib_name(ctl%attribute, slave, n_pals, pals_name, factor)
      if (n_pals == 0) cycle   ! Attribute cannot be translated.

      do i = 1, size(control%stack)
        select case (control%stack(i)%type)
        case (constant$)      ! Something like "c_light". Note: AtomicAndPhysicalConstants is uppercase by PALS is lower case
          control%stack(i)%name = downcase(control%stack(i)%name)
        case default
          control%stack(i)%name = downcase(control%stack(i)%name)
        end select
      enddo

      bline = expression_stack_to_string(control%stack)

      do k = 1, n_pals
        name = trim(pals_ele_name(slave)) // '>' // trim(pals_name(k))

        if (factor(k) == 1.0_rp) then
          line = bline
        else
          line = re_str(factor(k)) // ' * (' // trim(bline) // ')'
        endif

        call find_index(name, ctl_param, ctl_indx, n_ctl, ixm)
        if (ixm == 0) then
          call find_index(name, ctl_param, ctl_indx, n_ctl, ixm, add_to_list = .true.)
          ctl_attrib(ixm) = pals_name(k)
          ctl_slave(ixm)%ele => slave
        endif

        if (ctl_expr(ixm) == '') then
          ctl_expr(ixm) = line
        else
          ctl_expr(ixm) = trim(ctl_expr(ixm)) // ' + ' // trim(line)
        endif
      enddo
    enddo

    ! A PALS multipole component may have contributions from Bmad attributes that are not
    ! controlled (EG: The bend angle contribution to Kn0). Such contributions are constant so
    ! just add them in.

    header_needed = .true.
    do i = 1, n_ctl
      line = ctl_expr(i)
      if (ele%key == overlay$) then
        resid = pals_control_residual(ctl_slave(i)%ele, ele, ctl_attrib(i))
        if (resid /= 0) line = trim(line) // ' + ' // re_str(resid)
      endif

      call write_str_dict(id2, 'controls', '- parameter', ctl_param(i), header_needed)
      call write_str_dict(id2+id_del, 'controls', 'expression', line, header_needed)
    enddo
  end select
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

if (ele%ix_ele == 0) name_out = 'beginning_b' // int_str(ele%ix_branch+1)

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

! Return Pals attribute name(s) given Bmad attribute name.
! Since a Bmad attribute may map to more than one PALS attribute (EG: HKICK of a tilted element
! maps to both the normal and skew n = 0 multipole components), up to two names are returned.
! n_pals = 0 => Attribute cannot be translated.
! The value of the PALS attribute pals_name(i) is factor(i) times the value of the Bmad attribute.
! no_print = True suppresses the error message for an attribute that cannot be translated.

subroutine pals_attrib_name(bmad_name, ele, n_pals, pals_name, factor, no_print)

type (ele_struct) ele

integer n_pals
character(*) bmad_name
character(40) pals_name(2)
real(rp) factor(2)
integer n
logical, optional :: no_print

! factor is the conversion factor from Bmad parameter to PALS parameter.
! For example, accounting for a length normalization.

n_pals = 1
pals_name = ''
factor = 1.0_rp
n = len_trim(bmad_name)

if ((bmad_name(1:1) == 'A' .or. bmad_name(1:1) == 'B') .and. n >= 7 .and. bmad_name(max(1,n-4):n) == '_ELEC' .and. is_integer(bmad_name(2:max(2,n-5)), ix)) then
  pals_name(1) = 'ElectricMultipoleP.E'
  if (bmad_name(1:1) == 'A') then
    pals_name(1) = trim(pals_name(1)) // 's'
  else
    pals_name(1) = trim(pals_name(1)) // 'n'
  endif

  pals_name(1) = trim(pals_name(1)) // bmad_name(2:n-5)
  factor(1) = factorial(ix)
  return
endif

!

if ((bmad_name(1:1) == 'A' .or. bmad_name(1:1) == 'B') .and. is_integer(bmad_name(2:), ix)) then
  if (ele%field_master) then
    pals_name(1) = 'MagneticMultipoleP.B'
    factor(1) = ele%value(p0c$) / (charge_of(ele%ref_species) * c_light)
  else
    pals_name(1) = 'MagneticMultipoleP.K'
    factor(1) = 1
  endif

  if (bmad_name(1:1) == 'A') then
    pals_name(1) = trim(pals_name(1)) // 's'
  else
    pals_name(1) = trim(pals_name(1)) // 'n'
  endif

  pals_name(1) = trim(pals_name(1)) // bmad_name(2:)
  factor(1) = factor(1) * factorial(ix)

  ! Note: The element definition writes the integrated multipole (with an "L" suffix) only if the
  ! element has zero length. Otherwise the length normalized multipole is written.

  if (ele%value(l$) == 0) then
    pals_name(1) = trim(pals_name(1)) // 'L'
  else
    factor(1) = factor(1) / ele%value(l$)
  endif

  return
endif

! Kick attributes translate to n = 0 multipole components.

select case (bmad_name)
case ('KICK', 'HKICK', 'VKICK', 'BL_KICK', 'BL_HKICK', 'BL_VKICK')
  call pals_kick_attrib_name(bmad_name, ele, n_pals, pals_name, factor)
  return
end select

!

select case (bmad_name)
case ('B1_GRADIENT');   pals_name(1) = 'MagneticMultipoleP.Bn1'
case ('B2_GRADIENT');   pals_name(1) = 'MagneticMultipoleP.Bn2'
case ('B3_GRADIENT');   pals_name(1) = 'MagneticMultipoleP.Bn3'
case ('K1');            pals_name(1) = 'MagneticMultipoleP.Kn1'
case ('K2');            pals_name(1) = 'MagneticMultipoleP.Kn2'
case ('K3');            pals_name(1) = 'MagneticMultipoleP.Kn3'

case ('L');             pals_name(1) = 'length'

case ('E1');            pals_name(1) = 'BendP.e1'
case ('E2');            pals_name(1) = 'BendP.e2'
case ('G');             pals_name(1) = 'BendP.g_ref'
case ('ANGLE');         pals_name(1) = 'BendP.angle'
case ('TILT_REF');      pals_name(1) = 'BendP.tilt_ref'

case ('X1_LIMIT');      pals_name(1) = 'ApertureP.x_min'
case ('X2_LIMIT');      pals_name(1) = 'ApertureP.x_max'; factor(1) = -1
case ('Y1_LIMIT');      pals_name(1) = 'ApertureP.y_min'
case ('Y2_LIMIT');      pals_name(1) = 'ApertureP.y_max'; factor(1) = -1

case ('VOLTAGE');       pals_name(1) = 'RFP.voltage'
case ('GRADINET');      pals_name(1) = 'RFP.gradient'
case ('PHASE');         pals_name(1) = 'RFP.phase'
case ('RF_FREQUENCY');  pals_name(1) = 'RFP.frequency'

case ('T_OFFSET');      pals_name(1) = 't_offset'
case ('KS');            pals_name(1) = 'Ksol'
case ('BS_FIELD');      pals_name(1) = 'Bsol'

case ('X_OFFSET', 'Y_OFFSET', 'Z_OFFSET', 'X_PITCH', 'Y_PITCH', 'TILT')
  if (ele%key == patch$) then
    select case (bmad_name)
    case ('X_OFFSET');      pals_name(1) = 'PatchP.dx'
    case ('Y_OFFSET');      pals_name(1) = 'PatchP.dy'
    case ('Z_OFFSET');      pals_name(1) = 'PatchP.dz'
    case ('X_PITCH');       pals_name(1) = 'PatchP.dy_rot'
    case ('Y_PITCH');       pals_name(1) = 'PatchP.dx_rot'; factor(1) = -1
    case ('TILT');          pals_name(1) = 'PatchP.dz_rot'
    end select
  elseif (ele%key == floor_shift$ .or. ele%key == fiducial$) then
    select case (bmad_name)
    case ('X_OFFSET');      pals_name(1) = 'CoordinateSetP.dx'
    case ('Y_OFFSET');      pals_name(1) = 'CoordinateSetP.dy'
    case ('Z_OFFSET');      pals_name(1) = 'CoordinateSetP.dz'
    case ('X_PITCH');       pals_name(1) = 'CoordinateSetP.dy_rot'
    case ('Y_PITCH');       pals_name(1) = 'CoordinateSetP.dx_rot'; factor(1) = -1
    case ('TILT');          pals_name(1) = 'CoordinateSetP.dz_rot'
    end select
  else
    select case (bmad_name)
    case ('X_OFFSET');      pals_name(1) = 'BodyShiftP.x_offset'
    case ('Y_OFFSET');      pals_name(1) = 'BodyShiftP.y_offset'
    case ('Z_OFFSET');      pals_name(1) = 'BodyShiftP.z_offset'
    case ('X_PITCH');       pals_name(1) = 'BodyShiftP.y_rot'
    case ('Y_PITCH');       pals_name(1) = 'BodyShiftP.x_rot'; factor(1) = -1
    case ('TILT');          pals_name(1) = 'BodyShiftP.z_rot'
    end select
  endif

case default
  n_pals = 0
  if (.not. logic_option(.false., no_print)) then
    print *, 'Attribute not yet coded for translation: ' // trim(bmad_name)
    print *, 'Please report this.'
  endif
end select

end subroutine pals_attrib_name

!------------------------------------------------------
! contains

! Return PALS attribute name(s) for the Bmad kick attributes KICK, HKICK, VKICK and the
! corresponding integrated field attributes BL_KICK, BL_HKICK, BL_VKICK.
! In PALS a kick is represented by the n = 0 multipole components so a kick attribute of a
! tilted element must be distributed between the normal and skew components.
! The conversion mirrors what multipole_ele_to_ab does with the kick attributes (which is what is
! used when writing the element definitions) so that controlled values are consistent with the
! element definition values.

subroutine pals_kick_attrib_name(bmad_name, ele, n_pals, pals_name, factor)

type (ele_struct) ele

integer n_pals, key, i
real(rp) factor(2), f0, tilt, coef(2)
character(*) bmad_name
character(40) pals_name(2)
character(1) prefix
logical is_hkick

! is_hkick = True if the attribute gives a kick in the horizontal plane (in the element body frame).

key = ele%key
is_hkick = (index(bmad_name, 'VKICK') == 0)
if (key == vkicker$) is_hkick = .false.
if (key == hkicker$) is_hkick = .true.

! coef(1) is the normal (Kn0/Bn0) coefficient and coef(2) is the skew (Ks0/Bs0) coefficient.
! Note: For kicker type elements the kick is defined in the element body frame so there is no
! rotation by the element tilt.

select case (key)
case (hkicker$, vkicker$, kicker$, ac_kicker$)
  if (is_hkick) then
    coef = [-1.0_rp, 0.0_rp]
  else
    coef = [0.0_rp, 1.0_rp]
  endif

case (elseparator$)   ! Kick is electric
  if (ele%value(l$) == 0) then
    n_pals = 0
    return
  endif

  if (is_hkick) then
    pals_name(1) = 'ElectricMultipoleP.En0'
    factor(1) = -ele%value(p0c$) / ele%value(l$)
  else
    pals_name(1) = 'ElectricMultipoleP.Es0'
    factor(1) = ele%value(p0c$) / ele%value(l$)
  endif
  n_pals = 1
  return

case default
  if (key == sbend$ .or. key == rf_bend$) then
    tilt = ele%value(ref_tilt_tot$)
  else
    tilt = ele%value(tilt_tot$)
  endif

  if (is_hkick) then
    coef = [-cos(tilt), -sin(tilt)]
  else
    coef = [-sin(tilt), cos(tilt)]
  endif
end select

! BL_KICK, BL_HKICK and BL_VKICK are integrated field values so no scaling by the reference momentum.

if (bmad_name(1:3) == 'BL_') then
  prefix = 'B'
  f0 = 1
elseif (ele%field_master) then
  prefix = 'B'
  f0 = ele%value(p0c$) / (charge_of(ele%ref_species) * c_light)
else
  prefix = 'K'
  f0 = 1
endif

if (ele%value(l$) /= 0) f0 = f0 / ele%value(l$)

!

n_pals = 0

if (coef(1) /= 0) then
  n_pals = n_pals + 1
  pals_name(n_pals) = 'MagneticMultipoleP.' // prefix // 'n0'
  factor(n_pals) = f0 * coef(1)
endif

if (coef(2) /= 0) then
  n_pals = n_pals + 1
  pals_name(n_pals) = 'MagneticMultipoleP.' // prefix // 's0'
  factor(n_pals) = f0 * coef(2)
endif

if (ele%value(l$) == 0) then
  do i = 1, n_pals
    pals_name(i) = trim(pals_name(i)) // 'L'
  enddo
endif

end subroutine pals_kick_attrib_name

!------------------------------------------------------
! contains

! Return the constant (not controlled) part of the PALS multipole attribute pals_name of slave.
! Zero is returned if pals_name is not a multipole attribute or if there is nothing to add.
! Since several controllers may control a given component, and each controller writes its own
! expression for it, the constant part is only returned for the first controller of the component.
! Note: Group elements are ignored here since a group control is relative (a delta).

function pals_control_residual(slave, lord, pals_name) result (resid)

type (ele_struct) slave, lord
type (ele_struct), pointer :: this_lord
type (control_struct), pointer :: this_ctl

real(rp) resid, tot, sum_ctl, fact(2)
integer ix1, il, kk, np, n_done, ix_done(40)
character(*) pals_name
character(40) nam(2)
logical is_mult, is_new, contributes

!

resid = 0
tot = pals_multipole_value(slave, pals_name, is_mult)
if (.not. is_mult) return

! sum_ctl is the present value of the part of the component that is controlled.
! Multiple controllers may control a given attribute. In this case the attribute value is the sum
! of the contributions of all the controllers so only count the attribute value once.

sum_ctl = 0
n_done = 0
ix1 = 0

do il = 1, slave%n_lord
  this_lord => pointer_to_lord(slave, il, this_ctl)
  if (this_lord%key /= overlay$) cycle
  if (this_ctl%ix_attrib < 1 .or. this_ctl%ix_attrib > num_ele_attrib$) cycle
  call pals_attrib_name(this_ctl%attribute, slave, np, nam, fact, no_print = .true.)

  is_new = .true.
  do kk = 1, n_done
    if (ix_done(kk) == this_ctl%ix_attrib) is_new = .false.
  enddo

  contributes = .false.
  do kk = 1, np
    if (nam(kk) /= pals_name) cycle
    contributes = .true.
    if (is_new) sum_ctl = sum_ctl + fact(kk) * slave%value(this_ctl%ix_attrib)
  enddo

  if (.not. contributes) cycle
  if (ix1 == 0) ix1 = this_lord%ix_ele    ! First controller of this component.
  if (is_new .and. n_done < size(ix_done)) then
    n_done = n_done + 1
    ix_done(n_done) = this_ctl%ix_attrib
  endif
enddo

! Only the first controller of the component gets the constant part.

if (ix1 /= lord%ix_ele) return
if (abs(tot - sum_ctl) > 1e-14_rp * max(abs(tot), abs(sum_ctl))) resid = tot - sum_ctl

end function pals_control_residual

!------------------------------------------------------
! contains

! Return the value of the PALS multipole attribute pals_name as computed when writing the
! element definition. is_multipole is set False if pals_name is not a multipole attribute.
! This is needed since a given PALS multipole component may get contributions from several Bmad
! attributes (EG: Kn0 of a bend gets contributions from HKICK, VKICK, DG and the bend angle) and
! the part not controlled by a controller must be added in when writing a controlled value.

function pals_multipole_value(ele, pals_name, is_multipole) result (value)

type (ele_struct) ele

real(rp) value, ff, a_p(0:n_pole_maxx), b_p(0:n_pole_maxx)
integer nlen, nord, ixp
character(*) pals_name
character(40) nam
logical is_multipole

!

value = 0
is_multipole = .false.

ixp = index(pals_name, '.')
if (ixp == 0) return
nam = pals_name(ixp+1:)
if (nam(2:2) /= 'n' .and. nam(2:2) /= 's') return

! Electric multipoles are not scaled by the element length.

if (pals_name(1:ixp) == 'ElectricMultipoleP.') then
  if (.not. is_integer(nam(3:), nord)) return
  if (nord > n_pole_maxx) return
  call multipole_ele_to_ab(ele, .false., ixp, a_p, b_p, electric$, include_kicks$)
  ff = 1

elseif (pals_name(1:ixp) == 'MagneticMultipoleP.') then
  if ((nam(1:1) == 'B') .neqv. ele%field_master) return   ! Prefix is 'B' if and only if field_master.

  nlen = len_trim(nam)
  if (nam(nlen:nlen) == 'L') then
    if (ele%value(l$) /= 0) return
    nam = nam(1:nlen-1)
  else
    if (ele%value(l$) == 0) return
  endif

  if (.not. is_integer(nam(3:), nord)) return
  if (nord > n_pole_maxx) return

  call multipole_ele_to_ab(ele, .false., ixp, a_p, b_p, magnetic$, include_kicks$)
  if (ele%key == sbend$) b_p(0) = b_p(0) + ele%value(angle$)

  if (ele%field_master) then
    ff = ele%value(p0c$) / (charge_of(ele%ref_species) * c_light)
  else
    ff = 1
  endif

  if (ele%value(l$) /= 0) ff = ff / ele%value(l$)

else
  return
endif

!

if (nam(2:2) == 's') then
  value = ff * factorial(nord) * a_p(nord)
else
  value = ff * factorial(nord) * b_p(nord)
endif

is_multipole = .true.

end function pals_multipole_value

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

subroutine write_real_list(indnt, group_name, name, value, header_needed, ignore_zero)

real(rp) value
integer indnt, ind
logical, optional :: header_needed, ignore_zero
character(*) group_name, name
character(100) :: blank = ''

!

if (value == 0 .and. logic_option(.true., ignore_zero)) return
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

subroutine write_real_dict(indnt, group_name, name, value, header_needed, ignore_zero)

real(rp) value
integer indnt, ind
logical, optional :: header_needed, ignore_zero
character(*) group_name, name
character(100) :: blank = ''

!

if (value == 0 .and. logic_option(.true., ignore_zero)) return
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
