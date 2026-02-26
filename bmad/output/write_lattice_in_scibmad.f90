!+
! Subroutine write_lattice_in_scibmad(scibmad_file, lat, err_flag)
!
! Routine to create a SciBmad lattice file.
!
! Input:
!   lat           -- lat_struct: Lattice
!
! Output:
!   scibmad_file  -- character(*): SciBmad lattice file name.
!   err_flag      -- logical: Error flag
!-

subroutine write_lattice_in_scibmad(scibmad_file, lat, err_flag)

use write_lattice_file_mod, dummy => write_lattice_in_scibmad
use bmad_routine_interface, dummy2 => write_lattice_in_scibmad
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
type (ele_pointer_struct), allocatable :: named_eles(:)  ! List of unique element names 
type (lat_ele_order_struct) order
type (ele_attribute_struct) info
type (taylor_struct) taylor(6), spin_taylor(0:3)
type (nametable_struct) nametab
type (control_struct), pointer :: ctl
type (control_struct) control

real(rp) f, length, ang2
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)

integer n, i, j, k, ix, ib, ie, iu, is, n_names, ix_match, ix_pass, ix_r
integer ix_lord, ix_super, ie1, ib1
integer, allocatable :: an_indexx(:), index_list(:)

logical has_been_added, in_multi_region, have_expand_lattice_line, err, is_added
logical, optional :: err_flag

character(*) scibmad_file
character(1) prefix
character(3), parameter :: unit_spin_map(0:3) = ['1.0', '0.0', '0.0', '0.0']
character(40) name, look_for, ele_name
character(40), allocatable :: names(:)
character(240) fname
character(1000) line
character(*), parameter :: r_name = 'write_lattice_in_scibmad'

character(20) :: scibmad_ele_type(n_key$)

!

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
open (iu, file = fname, status = 'unknown')

! Header

write (iu, '(a)')  '# File generated by: write_lattice_in_foreign_format'
write (iu, '(4a)') '# Bmad lattice file: ', trim(lat%input_file_name)
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
allocate (names(n), an_indexx(n), named_eles(n))

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

    call add_this_name_to_list (ele, names, an_indexx, n_names, ix_match, has_been_added, named_eles)
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
      if (ele%value(fint$)*ele%value(hgap$) /= 0 .or. ele%value(fintx$)*ele%value(hgapx$) /= 0) print *, 'BEND EDGE_INT PARAMETER CANNOT YET BE TRANSLATED!'

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

    if (ele%key == lcavity$) then
      if (ele%value(rf_frequency$) /= 0)  line = trim(line) // ', rf_frequency = ' // re_str(ele%value(rf_frequency$))
      if (ele%value(voltage$) /= 0)  line = trim(line) // ', voltage = ' // re_str(ele%value(voltage$) + ele%value(voltage_err$))
      if (ele%value(phi0$) /= 0)  line = trim(line) // ', phi0 = ' // re_str(ele%value(phi0$) + ele%value(phi0_err$))
      line = trim(line) // ', tracking_method = SaganCavity(num_cells = ' // re_str(ele%value(n_rf_steps$)) // ', L_active = ' // re_str(ele%value(L_active$)) // ')'

    elseif (has_attribute(ele, 'RF_FREQUENCY')) then
      if (ele%key == rfcavity$) line = trim(line) // ', zero_phase = PhaseReference.AboveTransition'
      if (ele%value(rf_frequency$) /= 0)  line = trim(line) // ', rf_frequency = ' // re_str(ele%value(rf_frequency$))
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

    call write_lat_line(line, iu, .true., scibmad = .true.)

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
        endif
        in_multi_region = .true.
        ix_r = mult_ele(ie)%ix_region
        write (iu, '(a)')
        write (line, '(a, i2.2, a)') 'multi_line_', ix_r, ' = Beamline('
      endif

      if (mult_ele(ie)%ix_region /= ix_r) then
        call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #2! PLEASE REPORT THIS!')
      endif

      call write_scibmad_element (line, iu, ele, lat)

      if (mult_ele(ie)%region_stop_pt) then
        line = line(:len_trim(line)-1) // ')'
        call write_lat_line (line, iu, .true., scibmad = .true.)
        in_multi_region = .false.
      endif
    enddo

    if (in_multi_region) then
      call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #3! PLEASE REPORT THIS!')
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

call nametable_init(nametab)

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
        print *, 'Overlay: ' // trim(lord%name) // ' uses knot points for the control curve. This cannot yet be translated!'
        exit
      endif

      do k = 1, size(control%stack)
        if (control%stack(k)%type /= variable$) cycle
        call find_index(control%stack(k)%name, nametab, ix_match, add_to_list = .true., has_been_added = is_added)
        if (is_added) write (iu, '(2a, es24.17)') trim(control%stack(k)%name), ' = ', control%stack(k)%value
      enddo
    enddo
  endif
enddo

! 

lat%ele%select = .false.

do ie = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(ie)
  if (lord%key == overlay$) then
    call overlay_out(lord, lat)
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

  call write_lat_line (line, iu, .true., scibmad = .true.)
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
    !!! call eles_with_same_name_handler(ele, named_eles, an_indexx, names, n_names, order)
  enddo
enddo

! cleanup

close(iu)
deallocate (names, an_indexx)
!!! deallocate (mult_lat%branch)

!----------------------------------------------------------------------------------------------
contains

function unique_name(ele, lat) result (name)

type (lat_struct) lat
type (ele_struct) ele
integer n_match, ine, imax
character(50) name

!

imax = nametable_bracket_indexx(lat%nametable, ele%name, n_match)
if (n_match == 1) then
  name = ele%name
  return
endif

ine = ele_nametable_index(ele)
name = trim(ele%name) // '_' // int_str(ine - (imax - n_match))

end function unique_name

!----------------------------------------------------------------------------------------------
! contains

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

if (len_trim(line) > 100) call write_lat_line(line, iu, .false., scibmad = .true.)

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

subroutine eles_with_same_name_handler(ele, named_eles, an_indexx, names, n_names, order)

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0
type (lat_struct), pointer :: lat
type (ele_pointer_struct) :: named_eles(:)
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
ele0 => named_eles(ix_match)%ele   ! Element with this name whose attributes were written to the lattice file.
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

recursive subroutine overlay_out(overlay, lat)

type (lat_struct), target :: lat
type (ele_struct) overlay
type (ele_struct), pointer :: lord, slave
type (control_struct), pointer :: ctl
type (control_struct) control

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
  if (c_str(iv) == '') then
    write (iu, '(4a, es24.16)') trim(overlay%name), '_', trim(downcase(overlay%control%var(iv)%name)), ' = ', overlay%control%var(iv)%value
  else
    write (iu, '(7a)') 'const ', trim(overlay%name), '_', trim(downcase(overlay%control%var(iv)%name)), ' = DefExpr(() -> ', trim(c_str(iv)), ')'
  endif
enddo

! Now that this overlay has been outputted, check if any overlay slaves need outputting.

do ix = 1, overlay%n_slave
  slave => pointer_to_slave(overlay, ix)
  if (slave%key == overlay$) then
    call overlay_out(slave, lat)
  else
    call overlay_slave_out(slave, lat)
  endif
enddo

end subroutine overlay_out

!------------------------------------------------------
! contains

recursive subroutine overlay_slave_out(slave, lat)

type (lat_struct), target :: lat
type (ele_struct) slave
type (ele_struct), pointer :: lord
type (control_struct), pointer :: ctl
type (control_struct)  control

real(rp) f
integer ix, j, iv, it, n_contl, indx(40), ixm

character(40) sci_name, attrib_names(40)
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
  lord => pointer_to_lord(slave, ix)
  do j = 1, lord%n_slave
    slave2 => pointer_to_slave(lord, j, ctl)
    control = ctl
    if (slave2%ix_ele /= slave%ix_ele) cycle
    call find_index(control%attribute, attrib_names, indx, n_contl, ixm)
    if (ixm == 0) call find_index(control%attribute, attrib_names, indx, n_contl, ixm, add_to_list = .true.)

    if (c_str(ixm) == '') then
      c_str(ixm) = this_expression(control%stack, lord)
    else
      c_str(ixm) = trim(c_str(ixm)) // ' + ' // this_expression(control%stack, lord)
    endif
  enddo
enddo

do iv = 1, n_contl
  sci_name = scibmad_attrib_name(attrib_names(iv), slave, f)
  if (f == 1.0_rp) then
    write (iu, '(7a)') trim(downcase(slave%name)), '.', trim(sci_name), ' = DefExpr(() -> ', trim(c_str(iv)), ')'
  else
    write (iu, '(4a, es24.16, 3a)') trim(downcase(slave%name)), '.', trim(sci_name), ' = DefExpr(() -> ', f, ' * (', trim(c_str(iv)), '))'
  endif
enddo

end subroutine overlay_slave_out

!------------------------------------------------------
! contains

recursive function this_expression(stack, lord) result(expr)

type (expression_atom_struct) :: stack(:)
type (ele_struct) lord

integer i

character(1000) expr

!

do i = 1, size(stack)
  select case (stack(i)%type)
  case (constant$)      ! Something like "c_light"
    stack(i)%name = upcase(stack(i)%name)
  case (variable$)

  case default
    if (stack(i)%type > var_offset$ .and. stack(i)%type < var_offset$ + n_var_max$) then
      stack(i)%name = trim(lord%name) // '_' // downcase(stack(i)%name)
    endif
  end select
enddo

expr = expression_stack_to_string(stack)

end function this_expression

!------------------------------------------------------
! contains

! Return SciBmad attribute name given Bmad attribute name.

function scibmad_attrib_name(bmad_name, ele, factor) result (sci_name)

type (ele_struct) ele

character(*) bmad_name
character(40) sci_name
real(rp) factor

!

factor = 1.0_rp

if ((bmad_name(1:1) == 'A' .or. bmad_name(1:1) == 'B') .and. is_integer(bmad_name(2:), ix)) then
  if (ele%field_master) then
    sci_name = 'B'
    factor = ele%value(p0c$) / (charge_of(ele%ref_species) * c_light)
  else
    sci_name = 'K'
    factor = 1
  endif

  if (bmad_name(1:1) == 'A') then
    sci_name = sci_name(1:1) // 's'
  else
    sci_name = sci_name(1:1) // 'n'
  endif

  sci_name = trim(sci_name) // bmad_name(2:)
  factor = factor * factorial(ix)

  if (ele%value(l$) == 0) then
    sci_name = trim(sci_name) // 'L'
  else
    factor = factor / ele%value(l$)
  endif

  return
endif

select case (bmad_name)
case ('BL_KICK');       sci_name = 'Bn0L'
case ('BL_HKICK');      sci_name = 'Bn0L'
case ('BL_VKICK');      sci_name = 'Bs0L'
case ('B1_GRADIENT');   sci_name = 'Bn1'
case ('B2_GRADIENT');   sci_name = 'Bn2'
case ('B3_GRADIENT');   sci_name = 'Bn3'
case ('K1');            sci_name = 'Kn1'
case ('K2');            sci_name = 'Kn2'
case ('K3');            sci_name = 'Kn3'
case ('E1');            sci_name = 'e1'
case ('E2');            sci_name = 'e2'
case ('G');             sci_name = 'g_ref'
case ('ANGLE');         sci_name = 'g_ref'
case ('L');             sci_name = 'L'
case ('X_OFFSET', 'Y_OFFSET', 'Z_OFFSET', 'X_PITCH', 'Y_PITCH', 'TILT')
  if (ele%key == patch$) then
    select case (bmad_name)
    case ('X_OFFSET');      sci_name = 'dx'
    case ('Y_OFFSET');      sci_name = 'dy'
    case ('Z_OFFSET');      sci_name = 'dz'
    case ('X_PITCH');       sci_name = 'dy_rot'
    case ('Y_PITCH');       sci_name = 'dx_rot'; factor = -1
    case ('TILT');          sci_name = 'dz_rot'
    end select
  else
    select case (bmad_name)
    case ('X_OFFSET');      sci_name = 'x_offset'
    case ('Y_OFFSET');      sci_name = 'y_offset'
    case ('Z_OFFSET');      sci_name = 'z_offset'
    case ('X_PITCH');       sci_name = 'y_rot'
    case ('Y_PITCH');       sci_name = 'x_rot'; factor = -1
    case ('TILT');          sci_name = 'z_rot'
    end select
  endif

case ('T_OFFSET');      sci_name = 't_offset'
case ('KS');            sci_name = 'Ksol'
case ('BS_FIELD');      sci_name = 'Bsol'
case default
  print *, 'Attribute not yet coded for translation: ' // trim(bmad_name)
  print *, 'Please report this.'
end select

end function scibmad_attrib_name

end subroutine write_lattice_in_scibmad
