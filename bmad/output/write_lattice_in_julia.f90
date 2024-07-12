!+
! Subroutine write_lattice_in_julia(bmad_file, lat, julia_file)
!
! Routine to create a Bmad-Julia lattice file.
!
! Input:
!   lat           -- lat_struct: Lattice
!   bmad_file     -- character(*): Input Bmad lattice file name.
!                     If the name does not have a .jl suffix then this suffix will be added.
! Output:
!   julia_file    -- character(*), optional: Bmad-Julia lattice file name.
!-

subroutine write_lattice_in_julia(bmad_file, lat, julia_file)

use write_lat_file_mod, dummy => write_lattice_in_julia

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

real(rp) f, length, ang2
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)

integer n, i, j, ix, ib, ie, iu, n_names, ix_match, ix_pass, ix_r
integer ix_lord, ix_super, ie1, ib1
integer, allocatable :: an_indexx(:), index_list(:)

logical has_been_added, in_multi_region, have_expand_lattice_line, err

character(*) bmad_file
character(*), optional :: julia_file
character(1) prefix
character(40) name, look_for
character(40), allocatable :: names(:)
character(240) fname
character(1000) line
character(*), parameter :: r_name = 'write_lattice_in_julia'


character(20), parameter :: julia_name(n_key$) = [character(20):: &
    'Drift             ', 'Bend              ', 'Quadrupole        ', 'Group             ', 'Sextupole         ', &
    'Overlay           ', 'Custom            ', 'Taylor            ', 'RFCavity          ', 'ELSeparator       ', &
    'BeamBeam          ', 'Wiggler           ', 'Solenoid          ', 'Marker            ', 'Kicker            ', &
    'Hybrid            ', 'Octupole          ', 'Bend              ', 'Multipole         ', '!Bmad_Com         ', &
    '!Mad_Beam         ', 'Multipole         ', 'Solenoid          ', 'Patch             ', 'LCavity           ', &
    '!Parameter        ', 'NullEle           ', 'BeginningEle      ', '!Line             ', 'Match             ', &
    'Instrument        ', 'Instrument        ', 'Kicker            ', 'Kicker            ', 'Collimator        ', &
    'Collimator        ', 'Girder            ', 'Converter         ', '!Particle_Start   ', 'Fork              ', &
    'Fork              ', 'Mirror            ', 'Crystal           ', 'Pipe              ', 'Capillary         ', &
    'MultilayerMirror  ', 'EGun              ', 'EMField           ', 'FloorShift        ', 'Fiducial          ', &
    'Undulator         ', 'DiffractionPlatee ', 'PhotonInit        ', 'Sample            ', 'Detector          ', &
    'SadMult           ', 'Mask              ', 'ACKicker          ', 'Lens              ', '!Space_Charge_Com ', &
    'CrabCavity        ', 'Ramper            ', '!PTC_Com          ', 'RFBend            ', 'Kicker            ', &
    'Foil              ', 'ThickMultipole    ', 'Instrument        ', 'Instrument        ']


! Open file

call fullfilename(bmad_file, fname)
call file_suffixer(fname, fname, '.jl', .true.)
iu = lunget()
open (iu, file = fname, status = 'unknown')
if (present(julia_file)) julia_file = fname

! Write element defs

write (iu, '(a)') '# Lattice file translated from Fortran Bmad.'
write (iu, '(a)')

n_names = 0
n = lat%n_ele_max
allocate (names(n), an_indexx(n), named_eles(n))

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  ele_loop: do ie = 0, branch%n_ele_max
    ele => branch%ele(ie)
    length = ele%value(l$)

    if (ele%key == overlay$ .or. ele%key == group$ .or. ele%key == ramper$ .or. ele%key == girder$) cycle   ! Not currently handled
    if (ele%key == null_ele$) cycle

    multi_lord => pointer_to_multipass_lord (ele, ix_pass) 
    if (ele%lord_status == super_lord$ .and. ix_pass > 0) cycle
    if (ele%slave_status == super_slave$ .and. ix_pass > 1) cycle

    if (ele%slave_status == super_slave$) then
      lord => pointer_to_lord(ele, 1)
      slave => pointer_to_slave(lord, 1)
      slave2 => pointer_to_slave(lord, lord%n_slave)
      write (iu, '(2(a, i0), 2a)') '@ele slave_drift_', ib, '_', ele%ix_ele, ' = Drift(L = ', re_str(length) // ')'
      cycle
    endif

    if (ix_pass > 0) cycle

    ! Do not write anything for elements that have a duplicate name.

    call add_this_name_to_list (ele, names, an_indexx, n_names, ix_match, has_been_added, named_eles)
    if (.not. has_been_added) cycle

    ! Write element def
    ! The beginning element for all branches has the same name so use a unique name here.

    if (ie == 0) ele%name = 'begin' // int_str(ib+1)
    if (ie == branch%n_ele_track .and. ele%name == 'END') ele%name = 'end' // int_str(ib+1)
    line = '@ele ' // trim(downcase(ele%name)) // ' = ' // trim(julia_name(ele%key)) // '('

    if (ie == 0) then
      line = trim(line) // ', pc_ref = ' // re_str(ele%value(p0c$))
      line = trim(line) // ', species_ref = ' // trim(species_name(ele%ref_species))
      if (ele%a%beta /= 0) line = trim(line) // ', twiss.a.beta = ' // re_str(ele%a%beta)
      if (ele%b%beta /= 0) line = trim(line) // ', twiss.b.beta = ' // re_str(ele%b%beta)
      if (ele%a%alpha /= 0) line = trim(line) // ', twiss.a.alpha = ' // re_str(ele%a%alpha)
      if (ele%b%alpha /= 0) line = trim(line) // ', twiss.b.alpha = ' // re_str(ele%b%alpha)
      if (ele%x%eta /= 0) line = trim(line) // ', twiss.x.eta = ' // re_str(ele%x%eta)
      if (ele%y%eta /= 0) line = trim(line) // ', twiss.y.eta = ' // re_str(ele%y%eta)
      if (ele%x%etap /= 0) line = trim(line) // ', twiss.x.etap = ' // re_str(ele%x%etap)
      if (ele%y%etap /= 0) line = trim(line) // ', twiss.y.etap = ' // re_str(ele%y%etap)
      if (any(ele%c_mat /= 0)) line = trim(line) // ', twiss.c_mat = [' // re_str(ele%c_mat(1,1)) // ', ' // re_str(ele%c_mat(1,2)) // &
                                                                   '; ' // re_str(ele%c_mat(2,1)) // ', ' // re_str(ele%c_mat(2,2)) // ']'
      orb => lat%particle_start
      if (any(orb%vec /= 0)) line = trim(line) // ', particle.orbit = [' // re_str(orb%vec(1)) // ', ' // re_str(orb%vec(2)) // ', ' // &
                        re_str(orb%vec(3)) // ', ' // re_str(orb%vec(4)) // ', ' // re_str(orb%vec(5)) // ', ' // re_str(orb%vec(6)) // ']'
      if (any(orb%spin /= 0)) line = trim(line) // ', particle.spin = [' // &
                                             re_str(orb%spin(1)) // ', ' // re_str(orb%spin(2)) // ', ' //re_str(orb%spin(3)) // ']'

    endif

    if (ele%field_master) write (line, '(3a)') trim(line), ', field_master = ', jbool(ele%field_master)
    if (.not. ele%is_on) write (line, '(3a)') trim(line), ', is_on = ', jbool(ele%is_on)

    !

    if (ele%key == sbend$) then
      if (ele%sub_key == rbend$) then
        ang2 = 0.5_rp * ele%value(angle$) 
        line = trim(line) // ', bend_type = rbend'
        line = trim(line) // ', L_chord = ' // re_str(ele%value(l_chord$))
        if (abs(ele%value(e1$) - ang2) > 1d-14) line = trim(line) // ', e1_rect = ' // re_str(ele%value(e1$) - ang2)
        if (abs(ele%value(e2$) - ang2) > 1d-14) line = trim(line) // ', e2_rect = ' // re_str(ele%value(e2$) - ang2)
      else
        line = trim(line) // ', L = ' // re_str(length)
        if (ele%value(e1$) /= 0) line = trim(line) // ', e1 = ' // re_str(ele%value(e1$))
        if (ele%value(e2$) /= 0) line = trim(line) // ', e2 = ' // re_str(ele%value(e2$))
      endif

      if (ele%value(g$) /= 0)  line = trim(line) // ', g = ' // re_str(ele%value(g$))
      if (ele%value(ref_tilt$) /= 0)  line = trim(line) // ', ref_tilt = ' // re_str(ele%value(ref_tilt$))
      if (ele%value(roll$) /= 0)  line = trim(line) // ', roll = ' // re_str(ele%value(roll$))
      if (ele%value(fint$) /= 0)  line = trim(line) // ', fint1 = ' // re_str(ele%value(fint$))
      if (ele%value(fintx$) /= 0) line = trim(line) // ', fint2 = ' // re_str(ele%value(fintx$))
      if (ele%value(hgap$) /= 0)  line = trim(line) // ', hgap1 = ' // re_str(ele%value(hgap$))
      if (ele%value(hgapx$) /= 0) line = trim(line) // ', hgap2 = ' // re_str(ele%value(hgapx$))

    elseif (has_attribute(ele, 'L')) then
      if (length /= 0) line = trim(line) // ', L = ' // re_str(length)
    endif

    ! Magnetic multipoles

    call multipole_ele_to_ab(ele, .false., ix, a_pole, b_pole, magnetic$, include_kicks$)

    if (ele%field_master) then
      f = ele%value(p0c$) / (charge_of(ele%ref_species) * c_light)
      prefix = 'B'
    else
      f = 1
      prefix = 'K'
    endif

    do j = 0, ix
      if (length == 0) then
        if (a_pole(j) /= 0) line = trim(line) // ', ' // prefix // 's' // int_str(j) // 'L = ' // re_str(f * factorial(j) * a_pole(j))
        if (b_pole(j) /= 0) line = trim(line) // ', ' // prefix // 'n' // int_str(j) // 'L = ' // re_str(f * factorial(j) * b_pole(j))
      else
        f = f / length
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
      if (ele%value(x1_limit$) /= 0 .or. ele%value(x2_limit$) /= 0) line = trim(line) // 'x_limit = [' // &
                trim(aper_str(-ele%value(x1_limit$))) // ', ' // trim(aper_str(ele%value(x2_limit$))) // ']'
      if (ele%value(y1_limit$) /= 0 .or. ele%value(y2_limit$) /= 0) line = trim(line) // 'y_limit = [' // &
                trim(aper_str(-ele%value(y1_limit$))) // ', ' // trim(aper_str(ele%value(y2_limit$))) // ']'
    endif

    !

    if (has_attribute(ele, 'X_PITCH')) then
      if (ele%value(x_offset$) /= 0)  line = trim(line) // ', x_offset = ' // re_str(ele%value(x_offset$))
      if (ele%value(y_offset$) /= 0)  line = trim(line) // ', y_offset = ' // re_str(ele%value(y_offset$))
      if (ele%value(z_offset$) /= 0)  line = trim(line) // ', z_offset = ' // re_str(ele%value(z_offset$))
      if (ele%value(x_pitch$) /= 0)  line = trim(line) // ', x_pitch = ' // re_str(ele%value(x_pitch$))
      if (ele%value(y_pitch$) /= 0)  line = trim(line) // ', y_pitch = ' // re_str(ele%value(y_pitch$))
    endif

    if (has_attribute(ele, 'TILT')) then
      if (ele%value(tilt$) /= 0)  line = trim(line) // ', tilt = ' // re_str(ele%value(tilt$))
    endif

    !

    if (has_attribute(ele, 'KS')) then
      if (ele%field_master) then
        if (ele%value(bs_field$) /= 0)  line = trim(line) // ', bsol_field = ' // re_str(ele%value(bs_field$))
      else
        if (ele%value(ks$) /= 0)  line = trim(line) // ', ksol = ' // re_str(ele%value(ks$))
      endif
    endif

    !

    if (ele%key == lcavity$) then
      if (ele%value(rf_frequency$) /= 0)  line = trim(line) // ', rf_frequency = ' // re_str(ele%value(rf_frequency$))
      if (ele%value(voltage$) /= 0)  line = trim(line) // ', voltage_ref = ' // re_str(ele%value(voltage$))
      if (ele%value(voltage_err$) /= 0)  line = trim(line) // ', voltage_err = ' // re_str(ele%value(voltage_err$))
      if (ele%value(phi0$) /= 0)  line = trim(line) // ', phase_ref = ' // re_str(ele%value(phi0$))
      if (ele%value(phi0_err$) /= 0)  line = trim(line) // ', phase_err = ' // re_str(ele%value(phi0_err$))

    elseif (has_attribute(ele, 'RF_FREQUENCY')) then
      if (ele%value(rf_frequency$) /= 0)  line = trim(line) // ', rf_frequency = ' // re_str(ele%value(rf_frequency$))
      if (ele%value(voltage$) /= 0)  line = trim(line) // ', voltage = ' // re_str(ele%value(voltage$))
      if (ele%value(phi0$) /= 0)  line = trim(line) // ', phase = ' // re_str(ele%value(phi0$))
    endif

    if (has_attribute(ele, 'N_CELL')) then
      if (ele%value(n_cell$) /= 0)  line = trim(line) // ', n_cell = ' // re_str(ele%value(n_cell$))
      if (nint(ele%value(cavity_type$)) == standing_wave$) then
        line = trim(line) // ', cavity_type = standing_wave'
      else
        line = trim(line) // ', cavity_type = traveling_wave'
      endif
    endif

    !

    if (ele%type /= ' ') line = trim(line) // ', type = ' // quote(ele%type)
    if (ele%alias /= ' ') line = trim(line) // ', alias = ' // quote(ele%alias)
    if (associated(ele%descrip)) line = trim(line) // ', descrip = ' // quote(ele%descrip)

    !

    if (ele%key == fork$ .or. ele%key == photon_fork$) then
      n = nint(ele%value(ix_to_branch$))
      line = trim(line) // ', to_line = ' // trim(lat%branch(n)%name)
      if (ele%value(ix_to_element$) > 0) then
        i = nint(ele%value(ix_to_element$))
        line = trim(line) // ', to_element = ' // trim(lat%branch(n)%ele(i)%name)
      endif
    endif

    !

    ix = index(line, '(, ')
    if (ix == 0) then
      line = trim(line) // ')'
    else
      line = line(1:ix) // trim(line(ix+3:)) // ')'
    endif

    call write_lat_line(line, iu, .true., julia = .true.)

  enddo ele_loop
enddo

write (iu, '(a)')

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
      write (line, '(a, i2.2, a)') 'multi_line_', ix_r, ' = beamline('
    endif

    if (mult_ele(ie)%ix_region /= ix_r) then
      call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #2! PLEASE REPORT THIS!')
    endif

    call write_julia_element (line, iu, ele, lat)

    if (mult_ele(ie)%region_stop_pt) then
      line = line(:len_trim(line)-1) // ')'
      call write_lat_line (line, iu, .true., julia = .true.)
      in_multi_region = .false.
    endif
  enddo

  if (in_multi_region) then
    call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #3! PLEASE REPORT THIS!')
  endif

enddo  ! ib branch loop


!------------------------------

! Lines for all the branches.
! If we get into a multipass region then name in the main_line list is "multi_line_nn".
! But only write this once.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  write (iu, '(a)')
  name = downcase(branch%name)
  if (name == '') name = 'lat_line'
  line = trim(name) // ' = beamline(' // quote(name) // ', [' // trim(branch%ele(0)%name) // ','

  in_multi_region = .false.
  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)

    e_info => m_info%branch(ib)%ele(ie)

    if (.not. e_info%multipass) then
      call write_julia_element (line, iu, ele, lat)
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

  line = line(:len_trim(line)-1) // '], geometry = ' // trim(downcase(geometry_name(branch%param%geometry))) // ')'
  call write_lat_line (line, iu, .true., julia = .true.)
enddo

! Define lat

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

! If there are multipass lines then expand the lattice and write out
! the post-expand info as needed.

have_expand_lattice_line = .false.
do ie = 1, lat%n_ele_max
  ele => lat%ele(ie)
  if (ele%slave_status == super_slave$) cycle

  if (ele%key == lcavity$ .or. ele%key == rfcavity$) then
    if (ele%value(phi0_multipass$) == 0) cycle
    if (.not. have_expand_lattice_line) call write_expand_lat_header (iu, have_expand_lattice_line)
    write (iu, '(3a)') trim(ele%name), '[phi0_multipass] = ', re_str(ele%value(phi0_multipass$))
  endif

enddo

! If there are lattice elements with duplicate names but differing parameters then
! Write the differences.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  do ie = 1, branch%n_ele_max
    ele => branch%ele(ie)
    if (ele%key == marker$ .and. ele%name == 'END') cycle
    if (ele%slave_status == super_slave$) cycle
    if (ele%slave_status == multipass_slave$) cycle
    !!! call eles_with_same_name_handler(ele, named_eles, an_indexx, names, n_names, order)
  enddo
enddo

! cleanup

close(iu)
deallocate (names, an_indexx)
deallocate (mult_lat%branch)

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

subroutine write_julia_element (line, iu, ele, lat)

type (lat_struct), target :: lat
type (ele_struct) :: ele
type (ele_struct), pointer :: lord, m_lord, slave

character(*) line
character(40) lord_name

integer iu, ix

!

if (ele%slave_status == super_slave$) then
  if (ele%orientation == 1) then
    write (line, '(a, 2(a, i0), a)') trim(line), ' slave_drift_', ele%ix_branch, '_', ele%ix_ele, ','
  else
    write (line, '(a, 2(a, i0), a)') trim(line), ' reverse(slave_drift_', ele%ix_branch, '_', ele%ix_ele, '),'
  endif

elseif (ele%slave_status == multipass_slave$) then
  lord => pointer_to_lord(ele, 1)
  write (line, '(4a)') trim(line), ' ', trim(downcase(lord%name)), ','

else
  if (ele%orientation == 1) then
    write (line, '(4a)') trim(line), ' ', trim(downcase(ele%name)), ','
  else
    write (line, '(4a)') trim(line), ' reverse(', trim(downcase(ele%name)), '),'
  endif
endif

if (len_trim(line) > 100) call write_lat_line(line, iu, .false., julia = .true.)

end subroutine write_julia_element

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

end subroutine write_lattice_in_julia
