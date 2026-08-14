!+
! Subroutine write_lattice_scibmad_format(scibmad_file, lat, err_flag)
!
! Routine to create a SciBmad lattice file.
!
! Input:
!   scibmad_file  -- character(*): SciBmad lattice file name.
!   lat           -- lat_struct: Lattice
!
! Output:
!   err_flag      -- logical, optional: Set True if there is a problem. That is, if the file could
!                     not be opened or if the lattice contains something that could not be
!                     translated. Set False otherwise.
!-

subroutine write_lattice_scibmad_format(scibmad_file, lat, err_flag)

use write_lattice_file_mod, dummy => write_lattice_scibmad_format
use bmad_routine_interface, dummy2 => write_lattice_scibmad_format
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

integer n, i, j, k, ix, ib, ie, iu, is, n_names, ix_match, ix_pass, ix_r, ios
integer ix_lord, ix_super, ie1, ib1
integer, allocatable :: an_indexx(:), index_list(:)

logical has_been_added, in_multi_region, have_expand_lattice_line, err, is_added, has_defexpr_var
logical xlate_err    ! Set True if something in the lattice cannot be translated.
logical, optional :: err_flag

character(*) scibmad_file
character(1) prefix
character(3), parameter :: unit_spin_map(0:3) = ['1.0', '0.0', '0.0', '0.0']
character(100) name, look_for, ele_name
character(40), allocatable :: names(:)
character(240) fname
character(4000) line
character(*), parameter :: r_name = 'write_lattice_scibmad_format'

character(20) :: scibmad_ele_type(n_key$)

! err_flag is set False only after the file has been successfully written.

if (present(err_flag)) err_flag = .true.
xlate_err = .false.

scibmad_ele_type(drift$)                = 'Drift'
scibmad_ele_type(sbend$)                = 'SBend'
scibmad_ele_type(quadrupole$)           = 'Quadrupole'
scibmad_ele_type(group$)                = '??Group'
scibmad_ele_type(sextupole$)            = 'Sextupole'
scibmad_ele_type(overlay$)              = '??Overlay'
scibmad_ele_type(custom$)               = 'Custom'
scibmad_ele_type(taylor$)               = 'LineElement'
scibmad_ele_type(rfcavity$)             = 'RFCavity'
scibmad_ele_type(elseparator$)          = 'ELSeparator'
scibmad_ele_type(beambeam$)             = 'BeamBeam'
scibmad_ele_type(wiggler$)              = 'Wiggler'
scibmad_ele_type(sol_quad$)             = 'Solenoid'
scibmad_ele_type(marker$)               = 'Marker'
scibmad_ele_type(kicker$)               = 'Kicker'
scibmad_ele_type(hybrid$)               = 'Hybrid'
scibmad_ele_type(octupole$)             = 'Octupole'
scibmad_ele_type(rbend$)                = 'SBend'
scibmad_ele_type(multipole$)            = 'Multipole'
scibmad_ele_type(ab_multipole$)         = 'Multipole'
scibmad_ele_type(solenoid$)             = 'Solenoid'
scibmad_ele_type(patch$)                = 'Patch'
scibmad_ele_type(lcavity$)              = 'RFCavity'
scibmad_ele_type(null_ele$)             = 'NullEle'
scibmad_ele_type(beginning_ele$)        = 'BeginningEle'
scibmad_ele_type(def_line$)             = '!Line'
scibmad_ele_type(match$)                = 'LineElement'
scibmad_ele_type(monitor$)              = 'Drift'
scibmad_ele_type(instrument$)           = 'Drift'
scibmad_ele_type(hkicker$)              = 'Kicker'
scibmad_ele_type(vkicker$)              = 'Kicker'
scibmad_ele_type(rcollimator$)          = 'Drift'
scibmad_ele_type(ecollimator$)          = 'Drift'
scibmad_ele_type(girder$)               = 'Girder'
scibmad_ele_type(converter$)            = 'Converter'
scibmad_ele_type(photon_fork$)          = 'Fork'
scibmad_ele_type(fork$)                 = 'Fork'
scibmad_ele_type(mirror$)               = 'Mirror'
scibmad_ele_type(crystal$)              = 'Crystal'
scibmad_ele_type(pipe$)                 = 'Drift'
scibmad_ele_type(capillary$)            = 'Capillary'
scibmad_ele_type(multilayer_mirror$)    = 'MultilayerMirror'
scibmad_ele_type(e_gun$)                = 'EGun'
scibmad_ele_type(em_field$)             = 'EMField'
scibmad_ele_type(floor_shift$)          = 'FloorShift'
scibmad_ele_type(fiducial$)             = 'Fiducial'
scibmad_ele_type(undulator$)            = 'Undulator'
scibmad_ele_type(diffraction_plate$)    = 'DiffractionPlate'
scibmad_ele_type(photon_init$)          = 'PhotonInit'
scibmad_ele_type(sample$)               = 'Sample'
scibmad_ele_type(detector$)             = 'Detector'
scibmad_ele_type(sad_mult$)             = 'SadMult'
scibmad_ele_type(mask$)                 = 'Mask'
scibmad_ele_type(ac_kicker$)            = 'ACKicker'
scibmad_ele_type(lens$)                 = 'Lens'
scibmad_ele_type(crab_cavity$)          = 'CrabCavity'
scibmad_ele_type(ramper$)               = 'Ramper'
scibmad_ele_type(def_ptc_com$)          = '!PTC_Com'
scibmad_ele_type(rf_bend$)              = 'RFBend'
scibmad_ele_type(gkicker$)              = 'Kicker'
scibmad_ele_type(foil$)                 = 'Foil'
scibmad_ele_type(thick_multipole$)      = 'ThickMultipole'
scibmad_ele_type(pickup$)               = 'Drift'
scibmad_ele_type(feedback$)             = 'Drift'
scibmad_ele_type(fixer$)                = 'Fixer'

! Open file

call fullfilename(scibmad_file, fname)
iu = lunget()
open (iu, file = fname, status = 'unknown', iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE FOR WRITING: ' // trim(fname))
  return
endif

! Header

write (iu, '(4a)') '# Translated from Bmad lattice file: ', trim(lat%input_file_name)
write (iu, '(a)')
write (iu, '(a)')  'using Beamlines'

! Write functions for Taylor elements

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do ie = 1, branch%n_ele_max
    ele => branch%ele(ie)
    length = ele%value(l$)

    select case (ele%key)
    case (match$)
      call mat6_to_taylor(ele%vec0, ele%mat6, taylor)
      call write_this_taylor(iu, ele, taylor)
      cycle

    case (taylor$)
      call write_this_taylor(iu, ele, ele%taylor)
      cycle
    end select

  enddo
enddo

! Write element defs

! Note: Beamlines cannot currently handle multipass nor superimpose so ignore.
! Stuff that is commented out due to this is marked by "!!!"

n_names = 0
n = lat%n_ele_max
allocate (names(n), an_indexx(n), named_eles_ptr(n))

write (iu, '(a)')
write (iu, '(a)') '@elements begin'

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  ele_loop: do ie = 1, branch%n_ele_track   !!! Note: Not n_ele_max since superimpose/multipass not handled
    ele => branch%ele(ie)
    length = ele%value(l$)
    ele_name = scibmad_ele_name(ele)

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

    ! Write element def
    ! The beginning element for all branches has the same name so use a unique name here.

    if (ie == 0) ele_name = 'begin' // int_str(ib+1)
    line = '  ' // trim(ele_name) // ' = ' // trim(scibmad_ele_type(ele%key)) // '('

    if (ie == 0) then  ! Currently not used since ie starts at 1.
      line = trim(line) // ', pc_ref = ' // re_str(ele%value(p0c$))
      line = trim(line) // ', species_ref = species(' // quote(openpmd_species_name(ele%ref_species)) // ')'
      if (ele%a%beta /= 0) line = trim(line) // ', beta_a = ' // re_str(ele%a%beta)
      if (ele%b%beta /= 0) line = trim(line) // ', beta_b = ' // re_str(ele%b%beta)
      if (ele%a%alpha /= 0) line = trim(line) // ', alpha_a = ' // re_str(ele%a%alpha)
      if (ele%b%alpha /= 0) line = trim(line) // ', alpha_b = ' // re_str(ele%b%alpha)
      if (ele%x%eta /= 0) line = trim(line) // ', eta_x = ' // re_str(ele%x%eta)
      if (ele%y%eta /= 0) line = trim(line) // ', eta_y = ' // re_str(ele%y%eta)
      if (ele%x%etap /= 0) line = trim(line) // ', etap_x = ' // re_str(ele%x%etap)
      if (ele%y%etap /= 0) line = trim(line) // ', etap_y = ' // re_str(ele%y%etap)
      !! if (any(ele%c_mat /= 0)) line = trim(line) // ', c_mat = [' // re_str(ele%c_mat(1,1)) // ', ' // re_str(ele%c_mat(1,2)) // &
      !!                                                             '; ' // re_str(ele%c_mat(2,1)) // ', ' // re_str(ele%c_mat(2,2)) // ']'
      orb => lat%particle_start
      if (any(orb%vec /= 0)) line = trim(line) // ', particle.orbit = [' // re_str(orb%vec(1)) // ', ' // re_str(orb%vec(2)) // ', ' // &
                        re_str(orb%vec(3)) // ', ' // re_str(orb%vec(4)) // ', ' // re_str(orb%vec(5)) // ', ' // re_str(orb%vec(6)) // ']'
      if (any(orb%spin /= 0)) line = trim(line) // ', particle.spin = [' // &
                                             re_str(orb%spin(1)) // ', ' // re_str(orb%spin(2)) // ', ' //re_str(orb%spin(3)) // ']'

    endif

    if (.not. ele%is_on) write (line, '(3a)') trim(line), ', is_on = ', jbool(ele%is_on)

    !

    if (ele%key == sbend$) then
      line = trim(line) // ', L = ' // re_str(length)
      if (ele%value(e1$) /= 0) line = trim(line) // ', e1 = ' // re_str(ele%value(e1$))
      if (ele%value(e2$) /= 0) line = trim(line) // ', e2 = ' // re_str(ele%value(e2$))

      if (ele%value(g$) /= 0)  line = trim(line) // ', g_ref = ' // re_str(ele%value(g$))
      if (ele%value(ref_tilt$) /= 0)  line = trim(line) // ', tilt_ref = ' // re_str(ele%value(ref_tilt$))
      if (ele%value(roll$) /= 0)  line = trim(line) // ', roll = ' // re_str(ele%value(roll$))
      !!! if (ele%value(fint$)*ele%value(hgap$) /= 0)    line = trim(line) // ', edge_int1 = ' // re_str(ele%value(fint$)*ele%value(hgap$))
      !!! if (ele%value(fintx$)*ele%value(hgapx$) /= 0)  line = trim(line) // ', edge_int2 = ' // re_str(ele%value(fintx$)*ele%value(hgapx$))
      if (ele%value(fint$)*ele%value(hgap$) /= 0 .or. ele%value(fintx$)*ele%value(hgapx$) /= 0) then
        print *, 'BEND EDGE_INT PARAMETER CANNOT YET BE TRANSLATED!'
        xlate_err = .true.
      endif

    elseif (has_attribute(ele, 'L')) then
      if (length /= 0) line = trim(line) // ', L = ' // re_str(length)
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

    do j = 0, ix
      if (length == 0) then
        if (a_pole(j) /= 0) line = trim(line) // ', ' // prefix // 's' // int_str(j) // 'L = ' // re_str(f * factorial(j) * a_pole(j))
        if (b_pole(j) /= 0) line = trim(line) // ', ' // prefix // 'n' // int_str(j) // 'L = ' // re_str(f * factorial(j) * b_pole(j))
      else
        if (a_pole(j) /= 0) line = trim(line) // ', ' // prefix // 's' // int_str(j) // ' = ' // re_str(f * factorial(j) * a_pole(j))
        if (b_pole(j) /= 0) line = trim(line) // ', ' // prefix // 'n' // int_str(j) // ' = ' // re_str(f * factorial(j) * b_pole(j))
      endif
    enddo

    ! Electric multipoles

    call multipole_ele_to_ab(ele, .false., ix, a_pole, b_pole, electric$, include_kicks$)

    do j = 0, ix
      if (a_pole(j) /= 0) line = trim(line) // ', Es' // int_str(j) // ' = ' // re_str(factorial(j) * a_pole(j))
      if (b_pole(j) /= 0) line = trim(line) // ', En' // int_str(j) // ' = ' // re_str(factorial(j) * b_pole(j))
    enddo

    !

    if (has_attribute(ele, 'X1_LIMIT')) then
      if (ele%value(x1_limit$) /= 0) line = trim(line) // ', x1_limit = ' // trim(aper_str(-ele%value(x1_limit$)))
      if (ele%value(x2_limit$) /= 0) line = trim(line) // ', x2_limit = ' // trim(aper_str(ele%value(x2_limit$)))
      if (ele%value(y1_limit$) /= 0) line = trim(line) // ', y1_limit = ' // trim(aper_str(-ele%value(y1_limit$)))
      if (ele%value(y2_limit$) /= 0) line = trim(line) // ', y2_limit = ' // trim(aper_str(ele%value(y2_limit$)))

      if (ele%value(x1_limit$) /= 0 .or. ele%value(x2_limit$) /= 0 .or. &
          ele%value(y1_limit$) /= 0 .or. ele%value(y2_limit$) /= 0) then
        if (ele%aperture_type == elliptical$) then
          line = trim(line) // ', aperture_shape = ApertureShape.Elliptical'
        else
          line = trim(line) // ', aperture_shape = ApertureShape.Rectangular'
        endif
      endif
    endif

    !


    select case (ele%key)
    case (match$, taylor$)
      line = trim(line) // ', transport_map = map_' // trim(ele_name)
    end select

    !

    if (ele%key == patch$) then
      if (ele%value(t_offset$) /= 0)      line = trim(line) // ', dt = ' // re_str(ele%value(t_offset$))
      if (ele%value(x_offset$) /= 0)      line = trim(line) // ', dx = ' // re_str(ele%value(x_offset$))
      if (ele%value(y_offset$) /= 0)      line = trim(line) // ', dy = ' // re_str(ele%value(y_offset$))
      if (ele%value(z_offset$) /= 0)      line = trim(line) // ', dz = ' // re_str(ele%value(z_offset$))
      if (ele%value(y_pitch$) /= 0)       line = trim(line) // ', dx_rot = ' // re_str(-ele%value(y_pitch$))
      if (ele%value(x_pitch$) /= 0)       line = trim(line) // ', dy_rot = ' // re_str(ele%value(x_pitch$))
      if (ele%value(tilt$) /= 0)          line = trim(line) // ', dz_rot = ' // re_str(ele%value(tilt$))
      if (ele%value(E_tot_offset$) /= 0)  line = trim(line) // ', dE_ref = ' // re_str(ele%value(E_tot_offset$))
      if (ele%value(E_tot_set$) /= 0)     line = trim(line) // ', E_ref = ' // re_str(ele%value(E_tot_set$))

    else
      if (has_attribute(ele, 'X_PITCH')) then
        if (ele%value(x_offset$) /= 0)  line = trim(line) // ', x_offset = ' // re_str(ele%value(x_offset$))
        if (ele%value(y_offset$) /= 0)  line = trim(line) // ', y_offset = ' // re_str(ele%value(y_offset$))
        if (ele%value(z_offset$) /= 0)  line = trim(line) // ', z_offset = ' // re_str(ele%value(z_offset$))
        if (ele%value(y_pitch$) /= 0)  line = trim(line) // ', x_rot = ' // re_str(-ele%value(y_pitch$))
        if (ele%value(x_pitch$) /= 0)  line = trim(line) // ', y_rot = ' // re_str(ele%value(x_pitch$))
      endif

      if (has_attribute(ele, 'TILT')) then
        if (ele%value(tilt$) /= 0)  line = trim(line) // ', tilt = ' // re_str(ele%value(tilt$))
      endif
    endif

    !

    if (has_attribute(ele, 'KS')) then
      if (ele%field_master) then
        if (ele%value(bs_field$) /= 0)  line = trim(line) // ', bsol_field = ' // re_str(ele%value(bs_field$))
      else
        if (ele%value(ks$) /= 0)  line = trim(line) // ', Ksol = ' // re_str(ele%value(ks$))
      endif
    endif

    !

    if (has_attribute(ele, 'RF_FREQUENCY')) then
      if (is_true(ele%value(harmon_master$))) then
        if (ele%value(harmon$) /= 0)  line = trim(line) // ', harmon = ' // re_str(ele%value(harmon$))
      else
        if (ele%value(rf_frequency$) /= 0)  line = trim(line) // ', rf_frequency = ' // re_str(ele%value(rf_frequency$))
      endif
    endif

    if (ele%key == lcavity$) then
      if (ele%value(voltage$)+ele%value(voltage_err$) /= 0)  line = trim(line) // ', voltage = ' // re_str((ele%value(voltage$) + ele%value(voltage_err$)))
      if (ele%value(phi0$) /= 0)  line = trim(line) // ', phi0 = ' // re_str(ele%value(phi0$) + ele%value(phi0_err$))
      line = trim(line) // ', tracking_method = SaganCavity(num_cells = ' // re_str(ele%value(n_rf_steps$)) // ', L_active = ' // re_str(ele%value(L_active$)) // ')'

    elseif (has_attribute(ele, 'RF_FREQUENCY')) then
      if (ele%key == rfcavity$) line = trim(line) // ', zero_phase = PhaseRef.AboveTransition'
      if (ele%value(voltage$) /= 0)  line = trim(line) // ', voltage = ' // re_str(ele%value(voltage$)/abs(charge_of(branch%param%particle)))
      if (ele%value(phi0$) /= 0)  line = trim(line) // ', phi0 = ' // re_str(ele%value(phi0$))
    endif

    if (has_attribute(ele, 'CAVITY_TYPE')) then
      if (nint(ele%value(cavity_type$)) == standing_wave$) then
        line = trim(line) // ', traveling_wave = false'
      else
        line = trim(line) // ', traveling_wave = true'
      endif
    endif

    !

    if (ele%type /= ' ') line = trim(line) // ', label = ' // quote(ele%type)
    if (ele%alias /= ' ') line = trim(line) // ', alias = ' // quote(ele%alias)
    if (associated(ele%descrip)) line = trim(line) // ', description = ' // quote(ele%descrip)

    !

    if (ele%key == fork$ .or. ele%key == photon_fork$) then
      n = nint(ele%value(ix_to_branch$))
      line = trim(line) // ', to_line = ' // trim(downcase(lat%branch(n)%name))
      if (ele%value(ix_to_element$) > 0) then
        i = nint(ele%value(ix_to_element$))
        line = trim(line) // ', to_element = ' // trim(scibmad_ele_name(lat%branch(n)%ele(i)))
      endif
    endif

    !

    ix = index(line, '(, ')
    if (ix == 0) then
      line = trim(line) // ')'
    else
      line = line(1:ix) // trim(line(ix+3:)) // ')'
    endif

    call write_lat_line(line, iu, .true., ampersand_at_ends = .false.)

  enddo ele_loop
enddo

write (iu, '(a)') 'end    # @elements'

!------------------------------------------------------------------------------------------------------
! Write branch lines
! First write multipass lines

!!!!!!!!!!!!!!!!!!!!
if (.false.) then   !!!
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
          xlate_err = .true.
        endif
        in_multi_region = .true.
        ix_r = mult_ele(ie)%ix_region
        write (iu, '(a)')
        write (line, '(a, i2.2, a)') 'multi_line_', ix_r, ' = Beamline('
      endif

      if (mult_ele(ie)%ix_region /= ix_r) then
        call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #2! PLEASE REPORT THIS!')
        xlate_err = .true.
      endif

      call write_scibmad_element (line, iu, ele, lat)

      if (mult_ele(ie)%region_stop_pt) then
        line = line(:len_trim(line)-1) // ')'
        call write_lat_line (line, iu, .true., ampersand_at_ends = .false.)
        in_multi_region = .false.
      endif
    enddo

    if (in_multi_region) then
      call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #3! PLEASE REPORT THIS!')
      xlate_err = .true.
    endif
  enddo  ! ib branch loop
endif  !!!
!!!!!!!!!!!!!!!!!

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

  if (lord%key == girder$) then
    print *, 'GIRDER ELEMENTS CANNOT YET BE TRANSLATED!'
    xlate_err = .true.
    cycle
  endif

  if (lord%key == overlay$ .or. lord%key == group$) then
    do is = 1, lord%n_slave
      slave => pointer_to_slave(lord, is, ctl)
      control = ctl
      if (.not. allocated(control%stack)) then
        print *, trim(key_name(lord%key)) // ': ' // trim(lord%name) // &
                             ' uses knot points for the control curve. This cannot yet be translated!'
        xlate_err = .true.
        exit
      endif

      do k = 1, size(control%stack)
        if (control%stack(k)%type /= variable$) cycle
        call find_index(control%stack(k)%name, var_nametab, ix_match, add_to_list = .true., has_been_added = is_added)
        if (is_added) write (iu, '(2a, es24.17)') trim(control%stack(k)%name), ' = ', control%stack(k)%value
      enddo
    enddo
  endif
enddo

! 

lat%ele%select = .false.

do ie = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(ie)
  if (lord%key == overlay$ .or. lord%key == group$) then
    call controller_out(lord, lat, defexpr_nametab)
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
  if (name == '') name = 'lat_line'
  line = trim(name) // ' = Beamline(['     ! // quote(name) // ', [' // trim(scibmad_ele_name(branch%ele(0))) // ','

  in_multi_region = .false.
  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)

    !!! e_info => m_info%branch(ib)%ele(ie)

    !!!if (.not. e_info%multipass) then
      call write_scibmad_element (line, iu, ele, lat)
      cycle
    !!!endif

    ix_lord = e_info%ix_lord(1)
    ix_super = e_info%ix_super(1)
    ie1 = m_info%lord(ix_lord)%slave(1,ix_super)%ele%ix_ele
    ib1 = m_info%lord(ix_lord)%slave(1,ix_super)%ele%ix_branch
    m_ele => mult_lat%branch(ib1)%ele(ie1)
    ix_r = m_ele%ix_region

    ! If entering new multipass region
    if (.not. in_multi_region) then
      in_multi_region = .true.
      if (m_ele%region_start_pt) then
        write (line, '(2a, i2.2, a)') trim(line), ' multi_line_', ix_r, ','
        look_for = 'stop'
      else
        write (line, '(2a, i2.2, a)') trim(line), ' -multi_line_', ix_r, ','
        look_for = 'start'
      endif
    endif

    if (look_for == 'start' .and. m_ele%region_start_pt .or. &
        look_for == 'stop' .and. m_ele%region_stop_pt) then 
      in_multi_region = .false.
    endif
  enddo

  !!! line = line(:len_trim(line)-1) // '], geometry = ' // trim(downcase(geometry_name(branch%param%geometry))) // ')'
  !!! line = line(:len_trim(line)-1) // ']; R_ref = ' // &
  !!!                  trim(re_str(branch%ele(0)%value(p0c$)/charge_of(branch%param%particle))) // &
  !!!                  ', species_ref = Species(' // quote(openpmd_species_name(branch%param%particle)) // '))'
  line = line(:len_trim(line)-1) // ']; pc_ref = ' // trim(re_str(branch%ele(0)%value(p0c$))) // &
                    ', species_ref = Species(' // quote(openpmd_species_name(branch%param%particle)) // '))'

  call write_lat_line (line, iu, .true., ampersand_at_ends = .false.)
enddo

! Define lat

if (.false.) then
  line = 'lat = expand(' // quote(downcase(lat%use_name)) // ', ['
  do ib = 0, ubound(lat%branch, 1)
    branch => lat%branch(ib)
    if (branch%ix_from_branch > -1) cycle
    name = downcase(branch%name)
    if (name == '') name = 'lat_line'
    line = trim(line) // ', ' // name
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
    write (iu, '(3a)') trim(scibmad_ele_name(ele)), '[phi0_multipass] = ', re_str(ele%value(phi0_multipass$))
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

! cleanup

close(iu)
deallocate (names, an_indexx)
!!! deallocate (mult_lat%branch)

if (present(err_flag)) err_flag = xlate_err

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

function aper_str (limit) result (ap_str)

real(rp) limit
character(24) ap_str

!

if (limit == 0) then
  ap_str = 'NaN'
else
  ap_str = re_str(limit)
endif

end function aper_str

!----------------------------------------------------------------------------------------------
! contains

subroutine write_scibmad_element (line, iu, ele, lat)

type (lat_struct), target :: lat
type (ele_struct) :: ele
type (ele_struct), pointer :: lord, m_lord, slave

character(*) line
character(40) lord_name

integer iu, ix

!

!!!if (ele%slave_status == super_slave$) then
!!!  if (ele%orientation == 1) then
!!!    write (line, '(a, 2(a, i0), a)') trim(line), ' slave_drift_', ele%ix_branch, '_', ele%ix_ele, ','
!!!  else
!!!    write (line, '(a, 2(a, i0), a)') trim(line), ' reverse(slave_drift_', ele%ix_branch, '_', ele%ix_ele, '),'
!!!  endif
!!!
!!!elseif (ele%slave_status == multipass_slave$) then
!!!  lord => pointer_to_lord(ele, 1)
!!!  write (line, '(4a)') trim(line), ' ', trim(downcase(lord%name)), ','
!!!
!!!else
  if (ele%orientation == 1) then
    write (line, '(4a)') trim(line), ' ', trim(scibmad_ele_name(ele)), ','
  else
    write (line, '(4a)') trim(line), ' reverse(', trim(scibmad_ele_name(ele)), '),'
  endif
!!!endif

if (len_trim(line) > 100) call write_lat_line(line, iu, .false., ampersand_at_ends = .false.)

end subroutine write_scibmad_element

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

write (iu, '(5a)') trim(ele_unique_name(ele, order)), '[', trim(attrib_name), '] = ', re_str(value)

end subroutine write_this_differing_attrib

!--------------------------------------------------------------------------------
! contains

function scibmad_ele_name(ele) result (name_out)

type (ele_struct) ele
character(40) name_out
integer ix

name_out = downcase(ele%name)
if (name_out == 'end') name_out = 'end_b' // int_str(ele%ix_branch)

ix = index(name_out, '#')
if (ix /= 0) name_out = name_out(1:ix-1) // '!s' // name_out(ix+1:)

ix = index(name_out, '\')     !'
if (ix /= 0) name_out = name_out(1:ix-1) // '!m' // name_out(ix+1:)

end function scibmad_ele_name

!--------------------------------------------------------------------------------
! contains

subroutine write_this_taylor(iu, ele, taylor)

type (ele_struct) ele
type (taylor_struct), target :: taylor(6)
type (taylor_term_struct) term

integer iu, i, j, k
integer e_max(6)
character(200) line

!

write (iu, '(a)')
write (iu, '(9a)') 'function map_', trim(scibmad_ele_name(ele)), '(v, q)'

e_max = 0
do i = 1, 6
  do j = 1, 6
    e_max(j) = max(e_max(j), maxval(taylor(i)%term(:)%expn(j)))
  enddo
enddo

do i = 0, 3
  if (.not. associated(ele%spin_taylor(i)%term)) cycle
  if (size(ele%spin_taylor(i)%term) == 0) cycle
  do j = 1, 6
    e_max(j) = max(e_max(j), maxval(ele%spin_taylor(i)%term(:)%expn(j)))
  enddo
enddo

!

do i = 1, 6
  write (iu, '(2(a, i0))') '  v_out', i, '= '
  do j = 1, size(taylor(i)%term)
    term = taylor(i)%term(j)
    if (write_lat_debug_flag) then  ! Used for regression tests
      write (line, '(4x, es12.4)') term%coef
    else
      write (line, '(4x, es24.16)') term%coef
    endif

    do k = 1, 6
      if (term%expn(k) == 0) cycle
      if (term%expn(k) == 1) then
        write (line, '(a, 3(a, i0))') trim(line), '*v[', k, ']'
      else
        write (line, '(a, 3(a, i0))') trim(line), '*v[', k, ']^', term%expn(k)
      endif
    enddo
    if (j  == size(taylor(i)%term)) then
      write (iu, '(a)') line 
    else
      write (iu, '(a)') trim(line) // ' +' 
    endif
  enddo
enddo

!

write (iu, '(a)') 

do i = 0, 3
  if (.not. associated(ele%spin_taylor(i)%term)) then
    write (iu, '(a, i0, 2a)') ' q_out', i, ' = ', unit_spin_map(i)
    cycle
  elseif (size(ele%spin_taylor(i)%term) == 0) then
    write (iu, '(a, i0, 2a)') '  q_out', i, ' = ', unit_spin_map(i)
    cycle
  endif

  write (iu, '(2(a, i0))') '  q_out', i, ' = '
  do j = 1, size(ele%spin_taylor(i)%term)
    term = ele%spin_taylor(i)%term(j)
    if (write_lat_debug_flag) then  ! Used for regression tests
      write (line, '(4x, es13.5)') term%coef
    else
      write (line, '(4x, es24.16)') term%coef
    endif

    do k = 1, 6
      if (term%expn(k) == 0) cycle
      if (term%expn(k) == 1) then 
        write (line, '(a, 3(a, i0))') trim(line), '*q[', k, ']'
      else
        write (line, '(a, 3(a, i0))') trim(line), '*q[', k, ']^', term%expn(k)
      endif
    enddo
    if (j  == size(ele%spin_taylor(i)%term)) then
      write (iu, '(a)') line 
    else
      write (iu, '(a)') trim(line) // ' +' 
    endif
  enddo
enddo

write (iu, '(a)') 
write (iu, '(a)') '  return (v_out1, v_out2, v_out3, v_out4, v_out5, v_out6), (q_out1, q_out2, q_out3, q_out4)'
write (iu, '(a)') 'end'

end subroutine write_this_taylor

!------------------------------------------------------
! contains

! Output the control variables of an overlay or group element (called a "controller" here).

recursive subroutine controller_out(controller, lat, defexpr_nametab)

type (lat_struct), target :: lat
type (ele_struct) controller
type (ele_struct), pointer :: lord, slave
type (control_struct), pointer :: ctl
type (control_struct) control
type (nametable_struct) defexpr_nametab

integer ix, j, iv, it

character(1000) c_str(40), str
logical is_group, has_ovl(40), has_grp(40)

! Output is top down.
! Do not output if controller is already outputted or has lords that have not yet been outputted.

if (controller%select) return
do ix = 1, controller%n_lord
  lord => pointer_to_lord(controller, ix)
  if (lord%key /= overlay$ .and. lord%key /= group$) cycle
  if (.not. lord%select) return
enddo

! Output vars.
! Controled vars are defined with a defered expression.

controller%select = .true.
c_str = ''
has_ovl = .false.
has_grp = .false.
has_defexpr_var = .false.

do ix = 1, controller%n_lord
  lord => pointer_to_lord(controller, ix)
  if (lord%key /= overlay$ .and. lord%key /= group$) cycle
  is_group = (lord%key == group$)

  do j = 1, lord%n_slave
    slave => pointer_to_slave(lord, j, ctl)
    control = ctl
    if (slave%ix_ele /= controller%ix_ele) cycle
    if (.not. allocated(control%stack)) cycle   ! Knot point control. Message already given above.
    it = control%ix_attrib - var_offset$
    if (it < 1 .or. it > size(c_str)) cycle

    if (is_group) then
      has_grp(it) = .true.
      str = delta_expression(control, lord, has_defexpr_var)
    else
      has_ovl(it) = .true.
      str = this_expression(control%stack, lord, has_defexpr_var)
    endif

    if (c_str(it) == '') then
      c_str(it) = str
    else
      c_str(it) = trim(c_str(it)) // ' + ' // trim(str)
    endif
  enddo
enddo

do iv = 1, size(controller%control%var)
  name = trim(controller%name) // '_' // trim(downcase(controller%control%var(iv)%name))

  ! A group varies its slave incrementally so the present value of the var is the base value
  ! that the group deltas are added to.

  if (has_grp(iv) .and. .not. has_ovl(iv)) &
                c_str(iv) = re_str(controller%control%var(iv)%value) // ' + ' // trim(c_str(iv))

  if (c_str(iv) == '') then
    write (iu, '(2a, es24.16)') trim(name), ' = ', controller%control%var(iv)%value
  elseif (has_defexpr_var) then
    write (iu, '(3a)') 'if !@isdefined(', trim(name), ')'
    write (iu, '(7a)') '  const ', trim(name), ' = ', trim(c_str(iv))
    write (iu, '(a)')  'end'
    call nametable_add(defexpr_nametab, name, 1)
  else
    write (iu, '(3a)') 'if !@isdefined(', trim(name), ')'
    write (iu, '(7a)') '  const ', trim(name), ' = DefExpr(() -> ', trim(c_str(iv)), ')'
    write (iu, '(a)')  'end'
    call nametable_add(defexpr_nametab, name, 1)
  endif
enddo

! Now that this controller has been outputted, check if any slaves need outputting.

do ix = 1, controller%n_slave
  slave => pointer_to_slave(controller, ix)
  if (slave%key == overlay$ .or. slave%key == group$) then
    call controller_out(slave, lat, defexpr_nametab)
  else
    call controller_slave_out(slave, lat, defexpr_nametab)
  endif
enddo

end subroutine controller_out

!------------------------------------------------------
! contains

recursive subroutine controller_slave_out(slave, lat, defexpr_nametab)

type (lat_struct), target :: lat
type (ele_struct) slave
type (ele_struct), pointer :: lord
type (control_struct), pointer :: ctl
type (control_struct)  control
type (nametable_struct) defexpr_nametab

real(rp) factor(2), sum_ctl(40), grp_base(40), tot, resid
integer ix, j, k, iv, n_sci, n_contl, indx(40), ixm, n_done, ix_done(40)

character(40) sci_name(2), sci_names(40)
character(100) name
character(1000) :: c_str(40), str, term
logical has_defexpr_var(40), has_ovl(40), has_grp(40), hdv, is_new, is_mult, is_group

! Do not output if slave is already outputted or has controller lords that have not yet been outputted.

if (slave%select) return
do ix = 1, slave%n_lord
  lord => pointer_to_lord(slave, ix)
  if (lord%key /= overlay$ .and. lord%key /= group$) cycle
  if (.not. lord%select) return
enddo

! Output slave dependentcies.

slave%select = .true.

c_str = ''
sci_names = ''
n_contl = 0
n_done = 0
sum_ctl = 0
grp_base = 0
has_defexpr_var = .false.
has_ovl = .false.
has_grp = .false.

! Note: A single Bmad attribute (EG: HKICK of a tilted element) may map to multiple SciBmad
! attributes and multiple Bmad attributes may map to a single SciBmad attribute. So the
! expressions are collected by SciBmad attribute name.
! sum_ctl(i) is the present value of the part of the SciBmad attribute that is controlled.
! grp_base(i) is the present value of the part that is varied by a group. Unlike an overlay, a
! group varies an attribute incrementally so the present value is the base that deltas add to.

do ix = 1, slave%n_lord
  lord => pointer_to_lord(slave, ix, ctl)
  if (lord%key /= overlay$ .and. lord%key /= group$) cycle
  control = ctl
  if (.not. allocated(control%stack)) cycle   ! Knot point control. Message already given.
  is_group = (lord%key == group$)

  call scibmad_attrib_name(control%attribute, slave, n_sci, sci_name, factor)
  if (n_sci == 0) cycle    ! Attribute cannot be translated.

  ! Multiple lords may control a given attribute. In this case the attribute value is the sum of
  ! the contributions of all the lords so only count the attribute value once.

  is_new = (control%ix_attrib > 0 .and. control%ix_attrib <= num_ele_attrib$)
  do j = 1, n_done
    if (ix_done(j) == control%ix_attrib) is_new = .false.
  enddo
  if (is_new .and. n_done < size(ix_done)) then
    n_done = n_done + 1
    ix_done(n_done) = control%ix_attrib
  endif

  hdv = .false.
  if (is_group) then
    str = delta_expression(control, lord, hdv)
  else
    str = this_expression(control%stack, lord, hdv)
  endif

  do k = 1, n_sci
    call find_index(sci_name(k), sci_names, indx, n_contl, ixm)
    if (ixm == 0) call find_index(sci_name(k), sci_names, indx, n_contl, ixm, add_to_list = .true.)

    if (factor(k) == 1.0_rp) then
      term = str
    else
      term = re_str(factor(k)) // ' * (' // trim(str) // ')'
    endif

    if (c_str(ixm) == '') then
      c_str(ixm) = term
    else
      c_str(ixm) = trim(c_str(ixm)) // ' + ' // trim(term)
    endif

    if (is_group) then
      has_grp(ixm) = .true.
      if (is_new) grp_base(ixm) = grp_base(ixm) + factor(k) * slave%value(control%ix_attrib)
    else
      has_ovl(ixm) = .true.
      if (is_new) sum_ctl(ixm) = sum_ctl(ixm) + factor(k) * slave%value(control%ix_attrib)
    endif

    if (hdv) has_defexpr_var(ixm) = .true.
  enddo
enddo

! A SciBmad multipole component may have contributions from Bmad attributes that are not controlled
! (EG: The bend angle contribution to Kn0). Such contributions are constant so just add them in.
! Note: Since a group contribution is a delta with respect to the present attribute value, the
! group base value is part of the "not controlled" residual for a multipole component.

do iv = 1, n_contl
  name = trim(downcase(slave%name)) // '.' // trim(sci_names(iv))

  tot = scibmad_multipole_value(slave, sci_names(iv), is_mult)
  if (is_mult) then
    resid = tot - sum_ctl(iv)
    if (abs(resid) > 1e-14_rp * max(abs(tot), abs(sum_ctl(iv)))) &
                                            c_str(iv) = re_str(resid) // ' + ' // trim(c_str(iv))
  elseif (has_grp(iv) .and. .not. has_ovl(iv)) then
    c_str(iv) = re_str(grp_base(iv)) // ' + ' // trim(c_str(iv))
  endif

  if (has_defexpr_var(iv)) then
    write (iu, '(4a)') trim(name), ' = ', trim(c_str(iv))
  else
    write (iu, '(6a)') trim(name), ' = DefExpr(() -> ', trim(c_str(iv)), ')'
  endif

enddo

end subroutine controller_slave_out

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
  select case (downcase(stack(i)%name))
  case ('c_light', 'm_electron', 'm_proton', 'm_neutron', 'm_muon', 'm_pion_0', 'm_pion_charged', &
        'm_deuteron', 'm_helion', 'h_planck')
    stack(i)%name = upcase(stack(i)%name)
  case ('pi', 'sqrt', 'log', 'exp', 'sin', 'cos', 'tan', 'cot', 'asin', 'acos', 'atan', 'sinh', 'cosh', &
        'tanh', 'coth', 'asinh', 'acosh', 'atanh', 'acoth', 'abs', 'factorial', 'sign')
    stack(i)%name = downcase(stack(i)%name)
  case ('twopi')
    stack(i)%name = '2*pi'
  case ('fourpi')
    stack(i)%name = '4*pi'
  case ('e', 'e_log')
    stack(i)%name = 'exp(1.0)'
  case ('sqrt_2')
    stack(i)%name = 'sqrt(2.0)'
  case ('degrad')
    stack(i)%name = '(180 / pi)'
  case ('degrees', 'raddeg')
    stack(i)%name = '(pi / 180)'
  case ('r_e')
    stack(i)%name = 'R_ELECTRON'
  case ('r_p')
    stack(i)%name = 'R_PROTON'
  case ('h_bar_planck')
    stack(i)%name = 'H_BAR'
  case ('e_charge')
    stack(i)%name = 'E_CHARGE'
  case ('fine_struct_const')
    stack(i)%name = 'FINE_STRUCTURE'
  case ('emass')
    stack(i)%name = '(1e-9 * M_ELECTRON)'
  case ('pmass')
    stack(i)%name = '(1e-9 * M_PROTON)'
  case ('anom_moment_electron')
    stack(i)%name = 'ANOMALY_ELECTRON'
  case ('anom_moment_muon')
    stack(i)%name = 'ANOMALY_MUON'
  case ('anom_moment_proton')
    stack(i)%name = 'gyromagnetic_anomaly(Species("proton"))'
  case ('anom_moment_deuteron')
    stack(i)%name = 'gyromagnetic_anomaly(Species("deuteron"))'

  case ('atan2')
    stack(i)%name = 'atan'
  case ('modulo')
    stack(i)%name = 'mod'
  case ('sinc')
    stack(i)%name = 'sincu'
  case ('ran')
    stack(i)%name = 'rand'
  case ('ran_gauss')
    stack(i)%name = 'randn'
  case ('int')
    stack(i)%name = 'trunc'
  case ('nint')
    stack(i)%name = 'round'
  case ('floor')
    stack(i)%name = 'floor'
  case ('ceiling')
    stack(i)%name = 'ceil'
  case ('mass_of')
    stack(i)%name = 'massof'
  case ('charge_of')
    stack(i)%name = 'chargeof'
  case ('anomalous_modment_of')
    stack(i)%name = ''
  case ('species')
    stack(i)%name = 'Species'

  case default
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
  end select
enddo

expr = expression_stack_to_string(stack)

end function this_expression

!------------------------------------------------------
! contains

! A group element varies a controlled quantity Q incrementally:
!   Q -> Q + (E(v) - E(v0))
! where E is the control expression, v are the group variables and v0 are the variable values
! corresponding to the present value of Q. So the group contribution to Q is the returned
! delta expression and the present value of Q is the base value that the delta is added to.

function delta_expression(control, lord, has_defexpr_var) result (expr)

type (control_struct) control
type (ele_struct) lord

real(rp) val0
logical has_defexpr_var, err
character(1000) expr
character(100) err_str

! Note: E(v0) must be evaluated before this_expression is called since this_expression
! translates the atom names in the stack to their SciBmad equivalents.

val0 = 0
if (allocated(lord%control%var)) then
  val0 = expression_stack_value(control%stack, err, err_str, lord%control%var, .false.)
  if (err) then
    print *, 'Cannot evaluate control expression of group: ' // trim(lord%name)
    print *, err_str
    xlate_err = .true.
    val0 = 0
  endif
endif

expr = this_expression(control%stack, lord, has_defexpr_var)

if (val0 == 0) then
  expr = '(' // trim(expr) // ')'
else
  expr = '(' // trim(expr) // ' - (' // trim(re_str(val0)) // '))'
endif

end function delta_expression

!------------------------------------------------------
! contains

! Return SciBmad attribute name(s) given Bmad attribute name.
! Since a Bmad attribute may map to more than one SciBmad attribute (EG: HKICK of a tilted element
! maps to both the normal and skew n = 0 multipole components), up to two names are returned.
! n_sci = 0 => Attribute cannot be translated.
! The value of the SciBmad attribute sci_name(i) is factor(i) times the value of the Bmad attribute.

subroutine scibmad_attrib_name(bmad_name, ele, n_sci, sci_name, factor)

type (ele_struct) ele

integer n_sci
character(*) bmad_name
character(40) sci_name(2)
real(rp) factor(2)

!

n_sci = 1
sci_name = ''
factor = 1.0_rp

if ((bmad_name(1:1) == 'A' .or. bmad_name(1:1) == 'B') .and. is_integer(bmad_name(2:), ix)) then
  if (ele%field_master) then
    sci_name(1) = 'B'
    factor(1) = ele%value(p0c$) / (charge_of(ele%ref_species) * c_light)
  else
    sci_name(1) = 'K'
    factor(1) = 1
  endif

  if (bmad_name(1:1) == 'A') then
    sci_name(1) = sci_name(1)(1:1) // 's'
  else
    sci_name(1) = sci_name(1)(1:1) // 'n'
  endif

  sci_name(1) = trim(sci_name(1)) // bmad_name(2:)
  factor(1) = factor(1) * factorial(ix)

  if (ele%value(l$) == 0) then
    sci_name(1) = trim(sci_name(1)) // 'L'
  else
    factor(1) = factor(1) / ele%value(l$)
  endif

  return
endif

! Kick attributes translate to n = 0 multipole components.

select case (bmad_name)
case ('KICK', 'HKICK', 'VKICK', 'BL_KICK', 'BL_HKICK', 'BL_VKICK')
  call scibmad_kick_attrib_name(bmad_name, ele, n_sci, sci_name, factor)
  return
end select

select case (bmad_name)
case ('B1_GRADIENT');   sci_name(1) = 'Bn1'
case ('B2_GRADIENT');   sci_name(1) = 'Bn2'
case ('B3_GRADIENT');   sci_name(1) = 'Bn3'
case ('K1');            sci_name(1) = 'Kn1'
case ('K2');            sci_name(1) = 'Kn2'
case ('K3');            sci_name(1) = 'Kn3'
case ('E1');            sci_name(1) = 'e1'
case ('E2');            sci_name(1) = 'e2'
case ('G');             sci_name(1) = 'g_ref'
case ('ANGLE');         sci_name(1) = 'g_ref'
case ('L');             sci_name(1) = 'L'
case ('X_OFFSET', 'Y_OFFSET', 'Z_OFFSET', 'X_PITCH', 'Y_PITCH', 'TILT')
  if (ele%key == patch$) then
    select case (bmad_name)
    case ('X_OFFSET');      sci_name(1) = 'dx'
    case ('Y_OFFSET');      sci_name(1) = 'dy'
    case ('Z_OFFSET');      sci_name(1) = 'dz'
    case ('X_PITCH');       sci_name(1) = 'dy_rot'
    case ('Y_PITCH');       sci_name(1) = 'dx_rot'; factor(1) = -1
    case ('TILT');          sci_name(1) = 'dz_rot'
    end select
  else
    select case (bmad_name)
    case ('X_OFFSET');      sci_name(1) = 'x_offset'
    case ('Y_OFFSET');      sci_name(1) = 'y_offset'
    case ('Z_OFFSET');      sci_name(1) = 'z_offset'
    case ('X_PITCH');       sci_name(1) = 'y_rot'
    case ('Y_PITCH');       sci_name(1) = 'x_rot'; factor(1) = -1
    case ('TILT');          sci_name(1) = 'z_rot'
    end select
  endif

case ('T_OFFSET');      sci_name(1) = 't_offset'
case ('KS');            sci_name(1) = 'Ksol'
case ('BS_FIELD');      sci_name(1) = 'Bsol'

! These group specific attributes vary the lengths of neighboring elements. There is no
! SciBmad equivalent.

case ('START_EDGE', 'END_EDGE', 'ACCORDION_EDGE', 'S_POSITION', 'LORD_PAD1', 'LORD_PAD2')
  n_sci = 0
  xlate_err = .true.
  print *, 'Group control of the ' // trim(bmad_name) // ' attribute of element ' // &
                                                trim(ele%name) // ' cannot be translated.'

case default
  n_sci = 0
  xlate_err = .true.
  print *, 'Attribute not yet coded for translation: ' // trim(bmad_name)
  print *, 'Please report this.'
end select

end subroutine scibmad_attrib_name

!------------------------------------------------------
! contains

! Return SciBmad attribute name(s) for the Bmad kick attributes KICK, HKICK, VKICK and the
! corresponding integrated field attributes BL_KICK, BL_HKICK, BL_VKICK.
! In SciBmad a kick is represented by the n = 0 multipole components so a kick attribute of a
! tilted element must be distributed between the normal and skew components.
! The conversion mirrors what multipole_ele_to_ab does with the kick attributes (which is what is
! used when writing the element definitions) so that overlay/group controlled values are consistent
! with the element definition values.

subroutine scibmad_kick_attrib_name(bmad_name, ele, n_sci, sci_name, factor)

type (ele_struct) ele

integer n_sci, key, i
real(rp) factor(2), f0, tilt, coef(2)
character(*) bmad_name
character(40) sci_name(2)
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
    n_sci = 0
    return
  endif

  if (is_hkick) then
    sci_name(1) = 'En0'
    factor(1) = -ele%value(p0c$) / ele%value(l$)
  else
    sci_name(1) = 'Es0'
    factor(1) = ele%value(p0c$) / ele%value(l$)
  endif
  n_sci = 1
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

n_sci = 0

if (coef(1) /= 0) then
  n_sci = n_sci + 1
  sci_name(n_sci) = prefix // 'n0'
  factor(n_sci) = f0 * coef(1)
endif

if (coef(2) /= 0) then
  n_sci = n_sci + 1
  sci_name(n_sci) = prefix // 's0'
  factor(n_sci) = f0 * coef(2)
endif

if (ele%value(l$) == 0) then
  do i = 1, n_sci
    sci_name(i) = trim(sci_name(i)) // 'L'
  enddo
endif

end subroutine scibmad_kick_attrib_name

!------------------------------------------------------
! contains

! Return the value of the SciBmad multipole attribute sci_name as computed when writing the
! element definition. is_multipole is set False if sci_name is not a multipole attribute.
! This is needed since a given SciBmad multipole component may get contributions from several
! Bmad attributes (EG: Kn0 of a bend gets contributions from HKICK, VKICK, DG and the bend angle)
! and the part not controlled by an overlay must be added in when writing a controlled value.

function scibmad_multipole_value(ele, sci_name, is_multipole) result (value)

type (ele_struct) ele

real(rp) value, ff, a_p(0:n_pole_maxx), b_p(0:n_pole_maxx)
integer nlen, nord, ixp
character(*) sci_name
character(40) nam
logical is_multipole

!

value = 0
is_multipole = .false.

nam = sci_name
if (nam(2:2) /= 'n' .and. nam(2:2) /= 's') return

! Electric multipoles are not scaled by the element length.

if (nam(1:1) == 'E') then
  if (.not. is_integer(nam(3:), nord)) return
  if (nord > n_pole_maxx) return
  call multipole_ele_to_ab(ele, .false., ixp, a_p, b_p, electric$, include_kicks$)
  ff = 1

else
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
endif

!

if (nam(2:2) == 's') then
  value = ff * factorial(nord) * a_p(nord)
else
  value = ff * factorial(nord) * b_p(nord)
endif

is_multipole = .true.

end function scibmad_multipole_value

end subroutine write_lattice_scibmad_format
