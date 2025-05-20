!+
! Subroutine type_ele (ele, type_zero_attrib, type_mat6, type_taylor, twiss_out, type_control, 
!              type_wake, type_floor_coords, type_field, type_wall, type_rad_kick, lines, n_lines)
!
! Subroutine to print or put in a string array information on a lattice element.
! If the lines(:) argument is not present, the element information is printed to the terminal.
!
! Input:
!   ele               -- ele_struct: Element
!   type_zero_attrib  -- logical, optional: If False then surpress printing of
!                           real attributes whose value is 0 or switch attributes that have
!                           their default value. Default is False.
!   type_mat6         -- integer, optional:
!                            = 0   => Do not type ele%mat6
!                            = 4   => Type 4X4 xy submatrix
!                            = 6   => Type full 6x6 matrix (Default)
!   type_taylor       -- logical, optional: Print out taylor map terms?
!                          If ele%taylor is not allocated then this is ignored.
!                          Default is False.
!   twiss_out         -- integer, optional: Print the Twiss parameters at the element end?
!                          = 0         => Do not print the Twiss parameters
!                          = radians$  => Print Twiss, phi in radians (Default).
!                          = degrees$  => Print Twiss, phi in degrees.
!                          = cycles$   => Print Twiss, phi in radians/2pi.
!   type_control      -- integer, optional: Print control status? 
!                           If ele%branch%lat is not associated cannot print status info.
!                           = no$      => One line of info.
!                           = short$   => Almost all info except long knot point lists are truncated (default).
!                           = all$     => Everything.
!   type_wake         -- logical, optional: If True then print the long-range and 
!                          short-range wakes information. If False then just print
!                          how many terms the wake has. Default is True.
!                          If ele%wake is not allocated then this is ignored.
!   type_floor_coords -- logical, optional: Default is False. If True then print the global ("floor")
!                          coordinates at the exit end of the element.
!   type_field        -- integer, optional: Print field maps?
!                           = no$      => One line of info (default).
!                           = short$   => Header info. No tables.
!                           = all$     => Everything.
!   type_wall         -- logical, optional: Default is False. If True, print wall info. 
!   type_rad_kick     -- logical, optional: Default is False. If True, print synch rad kick info.
!
! Output       
!   lines(:)     -- character(200), allocatable, optional: Character array to hold the output. 
!                     If not present, the information is printed to the terminal.
!   n_lines      -- integer, optional: Number of lines in lines(:) that hold valid output.
!                     n_lines must be present if lines(:) is. 
!-

subroutine type_ele (ele, type_zero_attrib, type_mat6, type_taylor, twiss_out, type_control, &
             type_wake, type_floor_coords, type_field, type_wall, type_rad_kick, lines, n_lines)

use bmad_interface, except_dummy => type_ele
use expression_mod

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord, slave, ele0, lord2
type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (floor_position_struct) :: floor, f0, floor2
type (wake_lr_mode_struct), pointer :: lr
type (wake_sr_mode_struct), pointer :: mode
type (wake_sr_z_long_struct), pointer :: srz
type (cartesian_map_struct), pointer :: ct_map
type (cartesian_map_term1_struct), pointer :: ct_term
type (cylindrical_map_struct), pointer :: cl_map
type (cylindrical_map_term1_struct), pointer :: cl_term
type (grid_field_struct), pointer :: g_field
type (grid_field_pt1_struct), pointer :: g_pt
type (gen_grad_map_struct), pointer :: gg_map
type (gen_grad1_struct), pointer :: gg
type (wall3d_struct), pointer :: wall3d
type (wall3d_section_struct), pointer :: section
type (wall3d_vertex_struct), pointer :: v
type (photon_element_struct), pointer :: ph
type (ele_attribute_struct) attrib, attrib2
type (lat_param_struct) param
type (control_struct), pointer :: ctl
type (control_ramp1_struct), pointer :: rmp
type (all_pointer_struct) a_ptr
type (ac_kicker_struct), pointer :: ac
type (str_index_struct) str_index
type (rad_map_struct), pointer :: rm0, rm1
type (photon_reflect_table_struct), pointer :: rt
type (material_struct), pointer :: matter
type (surface_segmented_struct), pointer :: seg
type (surface_displacement_struct), pointer :: disp
type (surface_h_misalign_struct), pointer :: hmis
type (pixel_detec_struct), pointer :: pixel


integer, optional, intent(in) :: type_control, type_mat6, twiss_out, type_field
integer, optional, intent(out) :: n_lines
integer ia, ie, im, i1, ig, i, j, n, is, ix, iw, ix2_attrib, iv, ic, nl2, l_status, a_type, default_val
integer nl, nt, n_term, n_att, attrib_type, n_char, iy, particle, ix_pole_max, lb(2), ub(2)
integer id1, id2, id3, ne, na, nn

real(rp) coef, val, L_mis(3), S_mis(3,3), value
real(rp) a(0:n_pole_maxx), b(0:n_pole_maxx)
real(rp) a2(0:n_pole_maxx), b2(0:n_pole_maxx)
real(rp) knl(0:n_pole_maxx), tn(0:n_pole_maxx)
real(rp), pointer :: r_ptr

character(*), optional, allocatable :: lines(:)
character(200), allocatable, target :: li(:)
character(200), allocatable :: li2(:)
character(200) :: line
character(60) str1, str2
character(40) a_name, name, fmt_r, fmt_a, fmt_i, fmt_l, fmt
character(12) attrib_val_str, units, q_factor
character(8) angle, index_str

character(*), parameter :: r_name = 'type_ele'

logical, optional, intent(in) :: type_taylor, type_wake
logical, optional, intent(in) :: type_zero_attrib
logical, optional, intent(in) :: type_floor_coords, type_wall, type_rad_kick
logical type_zero, err_flag, print_it, is_default, has_it, has_been_added, is_zero1, is_zero2

! init

allocate (li(300))

type_zero = logic_option(.false., type_zero_attrib)
branch => pointer_to_branch(ele)

if (associated(branch)) then
  lat => branch%lat
  particle = branch%param%particle
  call lat_sanity_check(branch%lat, err_flag)
else
  nullify(lat)
  particle = electron$
endif

! Encode element name and type

nl = 0  

if (ele%ix_branch /= 0) then
  if (associated(lat)) then
    nl=nl+1; write (li(nl), '(a, i0, 2a)') 'Branch # ', ele%ix_branch, ': ', trim(branch%name)
  else
    nl=nl+1; write (li(nl), '(a, i0)') 'Branch # ', ele%ix_branch
  endif
endif
nl=nl+1; write (li(nl), '(a, i0)')     'Element # ', ele%ix_ele
nl=nl+1; write (li(nl), '(2a)')     'Element Name: ', trim(ele%name)

if (ele%type /= blank_name$) then
  nl=nl+1; write (li(nl), '(2a)') 'Element Type:  ', quote(ele%type)
endif

if (ele%alias /= blank_name$) then
  nl=nl+1; write (li(nl), '(2a)') 'Element Alias:  ', quote(ele%alias)
endif

if (associated(ele%descrip)) then
  nl=nl+1; write (li(nl), '(2a)') 'Descrip:  ', quote(ele%descrip)
endif

! Encode element key and attributes

if (ele%key <= 0) then
  nl=nl+1; write (li(nl), '(a, i0)') 'Key: BAD VALUE!', ele%key
else
  nl=nl+1; write (li(nl), '(2a)') 'Key: ', key_name(ele%key)
endif

if (ele%sub_key /= 0) then
  nl=nl+1; write (li(nl), '(2a)') 'Sub Key: ', sub_key_name(ele%sub_key)
endif

if (ele%key /= girder$ .and. ele%key /= ramper$ .and. ele%lord_status /= control_lord$) then
  nl=nl+1; write (li(nl), '(2(a, f14.6))')  'S_start, S:',  ele%s_start, ',', ele%s
  nl=nl+1; write (li(nl), '(2(a, es14.6))') 'Ref_time_start, Ref_time:', ele%value(ref_time_start$), ',', ele%ref_time
endif

nl=nl+1; li(nl) = ''
if (type_zero) then
  nl=nl+1; write (li(nl), '(a)') 'Attribute values:'
else
  nl=nl+1; write (li(nl), '(a)') 'Attribute values [Only non-zero values shown]:'
endif

n_att = n_attrib_string_max_len() + 1
write (fmt_a, '(a, i0, a)') '(7x, a, t', n_att+9, ', a, 2x, 3a)'
write (fmt_i, '(a, i0, a)') '(7x, a, t', n_att+9, ', a, i6)'
write (fmt_l, '(a, i0, a)') '(7x, a, t', n_att+9, ', a, 2x, l1)'
write (fmt_r, '(a, i0, a)') '(7x, a, t', n_att+9, ', a, 2x, es15.7)'


do ia = 1, num_ele_attrib$
  attrib = attribute_info(ele, ia)
  attrib%value = ele%value(ia)
  a_name = attrib%name
  if (a_name == null_name$) cycle
  if (a_name == 'MULTIPASS_REF_ENERGY' .and. (ele%lord_status /= multipass_lord$ .and. ele%slave_status /= multipass_slave$)) cycle
  if (a_name == 'LORD_PAD1' .and. ele%lord_status /= super_lord$) cycle
  if (attrib%state == private$) cycle
  if (is_2nd_column_attribute(ele, a_name, ix2_attrib)) cycle
  if (a_name == 'REF_SPECIES' .and. nint(ele%value(ia)) == not_set$) cycle

  attrib2 = ele_attribute_struct()

  select case (a_name)
  case ('SPECIES_STRONG'); attrib2 = ele_attribute_struct('BETA_STRONG', &
                              dependent$, is_real$, '', -1, ele%value(pc_strong$) / ele%value(E_tot_strong$))
  case ('LONGITUDINAL_MODE')
    if (ele%key == rfcavity$) attrib2 = ele_attribute_struct('RF_BUCKET_LENGTH', &
                                dependent$, is_real$, 'm', -1, ele%value(rf_wavelength$) * ele%value(p0c$) / ele%value(E_tot$))
  case ('P0C')
    if (particle == photon$) then
      attrib2 = ele_attribute_struct('REF_WAVELENGTH', dependent$, is_real$, 'm', -1, c_light * h_planck / ele%value(p0c$))
    else
      attrib2 = ele_attribute_struct('BETA', dependent$, is_real$, '', -1, ele%value(p0c$) / ele%value(e_tot$))
    endif
  case ('E_TOT'); if (particle /= photon$) attrib2 = ele_attribute_struct('GAMMA', &
                                             dependent$, is_real$, '', -1, ele%value(e_tot$) / mass_of(particle))
  case ('P0C_START'); attrib2 = ele_attribute_struct('BETA_START', &
                                  dependent$, is_real$, '', -1, ele%value(p0c_start$) / ele%value(e_tot_start$))
  case ('E_TOT_START'); attrib2 = ele_attribute_struct('DELTA_E', &
                                        dependent$, is_real$, 'eV', -1, ele%value(e_tot$) - ele%value(e_tot_start$))
  case ('DARWIN_WIDTH_SIGMA', 'DARWIN_WIDTH_PI')
    attrib2 = ele_attribute_struct(a_name, dependent$, is_real$, 'eV', -1, ele%value(ia) / ele%value(dbragg_angle_de$))
  case ('DBRAGG_ANGLE_DE'); attrib2 = ele_attribute_struct(a_name, dependent$, is_real$, 'deg/eV', -1, ele%value(ia) * 180 / pi)
  case ('DPHI_A', 'DPHI_B')
    attrib2 = ele_attribute_struct(a_name, dependent$, is_real$, 'rad/2pi', -1, ele%value(ia) * 1.0_rp / twopi)
  case default
    if (index(a_name, 'ANGLE') /= 0 .and. a_name /= 'CRITICAL_ANGLE_FACTOR') then
      attrib2 = ele_attribute_struct(a_name, dependent$, is_real$, 'deg', -1, ele%value(ia) * 180 / pi)
    else
      if (ix2_attrib > 0) then
        attrib2 = attribute_info(ele, ix2_attrib)
        attrib2%value = ele%value(ix2_attrib)
      endif
    endif
  end select

  is_zero1 = ((attrib%kind == is_real$ .or. attrib%kind == is_integer$) .and. attrib%value == 0)
  is_zero2 = (((attrib2%kind == is_real$ .or. attrib2%kind == is_integer$) .and. attrib2%value == 0) .or. &
                                                                                     attrib2%name == null_name$)
  if (is_zero1 .and. is_zero2 .and. .not. type_zero) cycle

  line = ''
  call write_this_attribute (attrib, ia, n_att, line(3:))
  call write_this_attribute (attrib2, attrib2%ix_attrib, 28, line(n_att+33:))
  nl=nl+1; li(nl) = line
enddo

! Custom attributes

if (ele%slave_status /= super_slave$) then
  do im = 1, custom_attribute_num$
    attrib = attribute_info(ele, im+custom_attribute0$)
    if (attrib%name(1:1) == '!') cycle
    nl=nl+1; write (li(nl), '(i5, 3x, 2a, es15.7, 3x, a)') &
                  im, attrib%name(1:n_att), '=', value_of_attribute(ele, attrib%name, err_flag), '! Custom attribute'
  enddo
endif

! Multipoles

if (associated(ele%a_pole) .or. associated(ele%a_pole_elec)) then
  nl=nl+1; write (li(nl), '(5x, a, l1)') 'MULTIPOLES_ON    = ', ele%multipoles_on 
endif

if (associated(ele%a_pole)) then
  if (attribute_index(ele, 'SCALE_MULTIPOLES') == scale_multipoles$) then
    nl=nl+1; write (li(nl), '(5x, a, l1, 2x, a)') 'SCALE_MULTIPOLES = ', ele%scale_multipoles, &
                                    '! Magnet strength scaling? Reference momentum scaling done if FIELD_MASTER = T.'
  endif

  if (associated(branch)) param = branch%param

  a = 0; b = 0; a2 = 0; b2 = 0; knl = 0; tn = 0
  if (ele%key == multipole$) then
    call multipole_ele_to_ab (ele, .false., ix_pole_max, a,  b)
    call multipole_ele_to_kt (ele, .true.,  ix_pole_max, knl, tn)
  else
    call multipole_ele_to_ab (ele, .false., ix_pole_max, a,  b)
    call multipole_ele_to_ab (ele, .true.,  ix_pole_max, a2, b2)
    call multipole_ele_to_kt (ele, .true.,  ix_pole_max, knl, tn)
  endif

  do im = 0, n_pole_maxx
    if (ele%key == multipole$) then
      if (a(im) == 0 .and. b(im) == 0 .and. tn(im) == 0) cycle

      nl=nl+1; write (li(nl), '(2x, 3(3x, a, i0, a, es11.3))') &
              'K', im, 'L       =', ele%a_pole(im), 'KS', im, '       =', ele%a_pole_elec(im), 'T', im, '        =', ele%b_pole(im)
      nl=nl+1; write (li(nl), '(2x, 3(3x, a, i0, a, es11.3))') &
              'B', im, '(equiv) =', b(im),          'A', im,  '(equiv) =', a(im),              'T', im, '(equiv) =', tn(im)

    elseif (ele%key == ab_multipole$) then
      if (a(im) == 0 .and. b(im) == 0) cycle

      nl=nl+1; write (li(nl), '(2x, 3(3x, a, i0, a, es11.3))') &
                 'A', im, ' =', ele%a_pole(im), 'A', im, '(w/Tilt) =', a2(im), 'K', im, 'L(equiv) =', knl(im)
      nl=nl+1; write (li(nl), '(2x, 3(3x, a, i0, a, es11.3))') &
                 'B', im, ' =', ele%b_pole(im), 'B', im, '(w/Tilt) =', b2(im), 'T', im, '(equiv)  =', tn(im)

    else
      if (ele%a_pole(im) == 0 .and. ele%b_pole(im) == 0 .and. a(im) == 0 .and. b(im) == 0) cycle

      nl=nl+1; write (li(nl), '(2x, 4(3x, a, i0, a, es11.3))') 'A', im, ' =', ele%a_pole(im), &
                 'A', im, '(Scaled) =', a(im), 'A', im, '(w/Tilt) =', a2(im), 'K', im, 'L(equiv) =', knl(im)
      nl=nl+1; write (li(nl), '(2x, 4(3x, a, i0, a, es11.3))') 'B', im, ' =', ele%b_pole(im), &
                 'B', im, '(Scaled) =', b(im), 'B', im, '(w/Tilt) =', b2(im), 'T', im, '(equiv)  =', tn(im)
    endif

  enddo
endif

! Electric Multipoles

if (associated(ele%a_pole_elec)) then
  call multipole_ele_to_ab (ele, .false., ix_pole_max, a, b, electric$)

  do im = 0, n_pole_maxx
    if (a(im) == 0 .and. b(im) == 0) cycle
    nl=nl+1; write (li(nl), '(2x, 4(3x, a, i0, a, es11.3))') 'A', im, '_elec =', ele%a_pole_elec(im), 'A', im, '_elec(Scaled) =', a(im)
    nl=nl+1; write (li(nl), '(2x, 4(3x, a, i0, a, es11.3))') 'B', im, '_elec =', ele%b_pole_elec(im), 'B', im, '_elec(Scaled) =', b(im)
  enddo
endif

! Encode on/off status etc.

if (.not. ele%is_on) then
  nl=nl+1; write (li(nl), *) '*** Note: Element is turned OFF ***'
endif

! Encode methods, etc.

nl=nl+1; write (li(nl), *) ''
nl2 = nl     ! For 2nd column parameters

if (attribute_name(ele, crystal_type$) == 'CRYSTAL_TYPE') then
  nl=nl+1; write (li(nl), fmt_a) 'CRYSTAL_TYPE', '=', ele%component_name
endif

if (attribute_name(ele, material_type$) == 'MATERIAL_TYPE') then
  nl=nl+1; write (li(nl), fmt_a) 'MATERIAL_TYPE', '=', ele%component_name
endif

if (attribute_name(ele, origin_ele$) == 'ORIGIN_ELE') then
  nl=nl+1; write (li(nl), fmt_a) 'ORIGIN_ELE', '=', quote(ele%component_name)
endif

if (attribute_name(ele, physical_source$) == 'PHYSICAL_SOURCE') then
  nl=nl+1; write (li(nl), fmt_a) 'PHYSICAL_SOURCE', '=', quote(ele%component_name)
endif

if (attribute_name(ele, tracking_method$) == 'TRACKING_METHOD') then
  nl=nl+1; write (li(nl), fmt_a) &
                  'TRACKING_METHOD', '=', tracking_method_name(ele%tracking_method)
endif

if (attribute_name(ele, mat6_calc_method$) == 'MAT6_CALC_METHOD') then
  nl=nl+1; write (li(nl), fmt_a) &
                  'MAT6_CALC_METHOD', '=', mat6_calc_method_name(ele%mat6_calc_method)
endif

if (attribute_name(ele, spin_tracking_method$) == 'SPIN_TRACKING_METHOD') then
  nl=nl+1; write (li(nl), fmt_a) &
                  'SPIN_TRACKING_METHOD', '=', spin_tracking_method_name(ele%spin_tracking_method)
endif

if (attribute_name(ele, ptc_integration_type$) == 'PTC_INTEGRATION_TYPE') then
  nl=nl+1; write (li(nl), fmt_a) &
                  'PTC_INTEGRATION_TYPE', '=', ptc_integration_type_name(ele%ptc_integration_type)
endif

! csr_method and space_charge_method not defined for multipass_lord elements.

if (ele%lord_status /= multipass_lord$ .and. attribute_name(ele, csr_method$) == 'CSR_METHOD') then
  nl=nl+1; write (li(nl), fmt_a) &
                  'CSR_METHOD', '=', csr_method_name(ele%csr_method)
endif

if (ele%lord_status /= multipass_lord$ .and. attribute_name(ele, space_charge_method$) == 'SPACE_CHARGE_METHOD') then
  nl=nl+1; write (li(nl), fmt_a) &
                  'SPACE_CHARGE_METHOD', '=', space_charge_method_name(ele%space_charge_method)
endif

if (attribute_name(ele, field_calc$) == 'FIELD_CALC') then
  nl=nl+1; write (li(nl), fmt_a) 'FIELD_CALC', '=', field_calc_name(ele%field_calc)
endif

! Write second column parameters

if (ele%key == beambeam$ .and. associated(branch)) then
  call encode_2nd_column_parameter (li, nl2, nl, 'PARAMETER[N_PART]', re_val = branch%param%n_part)
endif

if (attribute_name(ele, aperture_at$) == 'APERTURE_AT' .and. ele%aperture_at /= 0) then
  call encode_2nd_column_parameter (li, nl2, nl, 'APERTURE_AT', str_val = aperture_at_name(ele%aperture_at))
  call encode_2nd_column_parameter (li, nl2, nl, 'APERTURE_TYPE', str_val = aperture_type_name(ele%aperture_type))
  call encode_2nd_column_parameter (li, nl2, nl, 'OFFSET_MOVES_APERTURE', logic_val = ele%offset_moves_aperture)
endif

if (attribute_index(ele, 'SYMPLECTIFY') /= 0) then
  call encode_2nd_column_parameter (li, nl2, nl, 'SYMPLECTIFY', logic_val = ele%symplectify)
endif
  
if (attribute_index(ele, 'FIELD_MASTER') /= 0) then
  call encode_2nd_column_parameter (li, nl2, nl, 'FIELD_MASTER', logic_val = ele%field_master)
endif

if (ele%key /= overlay$ .and. ele%key /= group$ .and. &
          ele%key /= girder$ .and. ele%key /= ramper$ .and. ele%lord_status /= control_lord$) then
  call encode_2nd_column_parameter (li, nl2, nl, 'LONGITUDINAL ORIENTATION', int_val = ele%orientation)
endif

if (ele%ref_species /= not_set$) then
  call encode_2nd_column_parameter (li, nl2, nl, 'REF_SPECIES', str_val = species_name(ele%ref_species))
endif

! Foil

if (associated(ele%foil)) then
    nl=nl+1; li(nl) = ''
    nl=nl+1; li(nl) = 'Material_type: ' // ele%component_name

    do ix = 1, size(ele%foil%material)
      matter => ele%foil%material(ix)
      if (size(ele%foil%material) > 1) then
        nl=nl+1; li(nl) = ''
        nl=nl+1; write (li(nl), '(a, t24, a, t58, a, t97, a)') '  Component: ' // trim(species_name(matter%species)), &
                          attribute_units('DENSITY'), attribute_units('AREA_DENSITY'), attribute_units('RADIATION_LENGTH')
      else
        nl=nl+1; write (li(nl), '(t24, a, t58, a, t97, a)') &
                          attribute_units('DENSITY'), attribute_units('AREA_DENSITY'), attribute_units('RADIATION_LENGTH')
      endif

      nl=nl+1; write(li(nl), '(3(a, a14))') '  Density      =', this_real(matter%density, 'es14.6', '  Not_Set'), &
             '  Area_Density      =', this_real(matter%area_density, 'es14.6', '  Not_Set'), &
             '  Radiation_Length      =', this_real(matter%radiation_length, 'es14.6', '  Not_Set')
      nl=nl+1; write(li(nl), '((a, a14), 2(a, es14.6))') '  Density_Used =', this_real(matter%density_used, 'es14.6', '  Not_Used'), &
             '  Area_Density_Used =', matter%area_density_used, '  Radiation_Length_Used =', matter%radiation_length_used
    enddo
endif

! Converter

if (associated(ele%converter)) then
  do im = 1, size(ele%converter%dist)
  enddo
endif

! Cartesian map. The type_field logical is useful since field tables can be very large.

if (associated(ele%cartesian_map)) then
  if (integer_option(no$, type_field) == no$) then
    nl=nl+1; write (li(nl), '(a, i5)') 'Number of Cartesian_map modes:', size(ele%cartesian_map)
  else
    nl2 = 10; if (type_field == all$) nl2 = 999
    nl=nl+1; li(nl) = ''
    if (ele%field_calc == bmad_standard$) then
      nl=nl+1; li(nl) = 'Cartesian_map: [NOT USED SINCE FIELD_CALC = BMAD_STANDARD]'
    else
      nl=nl+1; li(nl) = 'Cartesian_map:'
    endif
    do im = 1, size(ele%cartesian_map)
      ct_map => ele%cartesian_map(im)
      if (ct_map%master_parameter == 0) then
        name = '<None>'
      else
        name = attribute_name(ele, ct_map%master_parameter)
      endif

      nl=nl+1; write (li(nl), '(a, i0)')      '  Cartesian_map mode #', im
      nl=nl+1; write (li(nl), '(2a)')         '    From file:        ', trim(ct_map%ptr%file)
      nl=nl+1; write (li(nl), '(2a)')         '    field_type        ', trim(em_field_type_name(ct_map%field_type))
      nl=nl+1; write (li(nl), '(2a)')         '    master_parameter: ', trim(name)
      nl=nl+1; write (li(nl), '(2a)')         '    ele_anchor_pt:    ', anchor_pt_name(ct_map%ele_anchor_pt)
      nl=nl+1; write (li(nl), '(a, es16.8)')  '    field_scale:      ', ct_map%field_scale
      nl=nl+1; write (li(nl), '(a, 3es16.8)') '    r0:               ', ct_map%r0
      nl=nl+1; write (li(nl), '(a, i0)')      '    n_link:           ', ct_map%ptr%n_link
      nl=nl+1; write (li(nl), '(3x, a, 13x, a, 3(9x, a), 2(12x, a), 9x, a, 3x, a)') 'Term#', &
                                    'A', 'K_x', 'K_y', 'K_z', 'x0', 'y0', 'phi_z', 'Family    Form'
      do j = 1, min(nl2, size(ct_map%ptr%term))
        if (nl+1 > size(li)) call re_allocate(li, 2 * nl, .false.)
        ct_term => ct_map%ptr%term(j)
        nl=nl+1; write (li(nl), '(i8, 1x, a13, 3f12.6, 3f14.6, 3x, a, 2x, a)') j, adjustr(real_to_string(ct_term%coef, 13, 7)), &
                        ct_term%kx, ct_term%ky, ct_term%kz, ct_term%x0, ct_term%y0, ct_term%phi_z, &
                        cartesian_map_family_name(ct_term%family), trim(cartesian_map_form_name(ct_term%form))
      enddo
      if (size(ct_map%ptr%term) > nl2) then
        nl=nl+1; write (li(nl), '(a, i0, a)') '     .... etc ... (#Terms = ', size(ct_map%ptr%term), ')' 
      endif
    enddo
  endif
endif

! Cylindrical_map

if (associated(ele%cylindrical_map)) then
  if (integer_option(no$, type_field) == no$) then
    nl=nl+1; write (li(nl), '(a, i5)') 'Number of Cylindrical_map modes:', size(ele%cylindrical_map)
  else
    nl2 = 10; if (type_field == all$) nl2 = 999
    nl=nl+1; li(nl) = ''
    if (ele%field_calc == bmad_standard$) then
      nl=nl+1; li(nl) = 'Cylindrical_map: [NOT USED SINCE FIELD_CALC = BMAD_STANDARD]'
    else
      nl=nl+1; li(nl) = 'Cylindrical_map:'
    endif
    do im = 1, size(ele%cylindrical_map)
      cl_map => ele%cylindrical_map(im)
      if (cl_map%master_parameter == 0) then
        name = '<None>'
      else
        name = attribute_name(ele, cl_map%master_parameter)
      endif

      nl=nl+1; write (li(nl), '(a, i0)')      '  Cylindrical_map mode #', im
      nl=nl+1; write (li(nl), '(2a)')         '    From file:        ', trim(cl_map%ptr%file)
      nl=nl+1; write (li(nl), '(2a)')         '    master_parameter: ', trim(name)
      nl=nl+1; write (li(nl), '(a, i0)')      '    harmonic:         ', cl_map%harmonic
      nl=nl+1; write (li(nl), '(a, i0)')      '    m:                ', cl_map%m
      nl=nl+1; write (li(nl), '(a, es16.8)')  '    field_scale:      ', cl_map%field_scale
      nl=nl+1; write (li(nl), '(a, es16.8)')  '    dz:               ', cl_map%dz
      nl=nl+1; write (li(nl), '(a, 3es16.8)') '    r0:               ', cl_map%r0
      nl=nl+1; write (li(nl), '(a, es16.8)')  '    phi0_fieldmap:    ', cl_map%phi0_fieldmap
      nl=nl+1; write (li(nl), '(a, es16.8)')  '    theta0_azimuth:   ', cl_map%theta0_azimuth
      nl=nl+1; write (li(nl), '(2a)')         '    ele_anchor_pt:    ', anchor_pt_name(cl_map%ele_anchor_pt)
      nl=nl+1; write (li(nl), '(a, i0)')      '    n_link:           ', cl_map%ptr%n_link
      nl=nl+1; write (li(nl), '(a)')          '    Term                E                           B'
      do j = 1, min(nl2, size(cl_map%ptr%term))
        if (nl+1 > size(li)) call re_allocate(li, 2 * nl, .false.)
        cl_term => cl_map%ptr%term(j)
        nl=nl+1; write (li(nl), '(i5, 3x, 2(a, 2es12.4), a)') j, '(', cl_term%e_coef, ')  (', cl_term%b_coef, ')'
      enddo
      if (size(cl_map%ptr%term) > nl2) then
        nl=nl+1; write (li(nl), '(a, i0, a)') '     .... etc ... (#Terms = ', size(cl_map%ptr%term), ')' 
      endif
    enddo
  endif
endif

! Gen_Grad_field

if (associated(ele%gen_grad_map)) then
  if (integer_option(no$, type_field) == no$) then
    nl=nl+1; write (li(nl), '(a, i5)') 'Number of Gen_Grad_field modes:', size(ele%gen_grad_map)
  else
    nl2 = 10; if (type_field == all$) nl2 = 999
    nl=nl+1; li(nl) = ''
    if (ele%field_calc == bmad_standard$) then
      nl=nl+1; li(nl) = 'Gen_Grad_field: [NOT USED SINCE FIELD_CALC = BMAD_STANDARD]'
    else
      nl=nl+1; li(nl) = 'Gen_Grad_field:'
    endif
    do im = 1, size(ele%gen_grad_map)
      gg_map => ele%gen_grad_map(im)
      if (gg_map%master_parameter == 0) then
        name = '<None>'
      else
        name = attribute_name(ele, gg_map%master_parameter)
      endif

      nl=nl+1; write (li(nl), '(2a)')         '  field_type:        ', em_field_type_name(gg_map%field_type)
      nl=nl+1; write (li(nl), '(a, es16.8)')  '  field_scale:       ', gg_map%field_scale
      nl=nl+1; write (li(nl), '(a, 3es16.8)') '  r0:                ', gg_map%r0
      nl=nl+1; write (li(nl), '(a, es16.8)')  '  dz:                ', gg_map%dz
      nl=nl+1; write (li(nl), '(a, 2i5)')     '  iz0, iz1:          ', gg_map%iz0, gg_map%iz1
      nl=nl+1; write (li(nl), '(2a)')         '  master_parameter:  ', trim(name)
      nl=nl+1; write (li(nl), '(2a)')         '  ele_anchor_pt:     ', anchor_pt_name(gg_map%ele_anchor_pt)
      nl=nl+1; write (li(nl), '(a, l1)')      '  curved_ref_frame   ', gg_map%curved_ref_frame
      do ig = 1, size(gg_map%gg)
        gg => gg_map%gg(ig)
        nl=nl+1; write (li(nl), '(4x, i0, a, i0, 2x, a)') ig, ': Curve m = ', gg%m, expression_op_name(gg%sincos)
        nl=nl+1; write (li(nl), '(9x, a)') 'Z  Derivs....'
        do n = lbound(gg%deriv,1), ubound(gg%deriv,1)
          if (n - lbound(gg%deriv,1) > 10) then
            nl=nl+1; li(nl) = '      ... etc. ...'
            exit
          endif
          nl=nl+1; write(li(nl), '(i10, f10.6, 99es14.6)') n, n*gg_map%dz, gg%deriv(n, 0:gg%n_deriv_max)
        enddo
      enddo
    enddo
  endif
endif

! Grid_field

if (associated(ele%grid_field)) then
  if (integer_option(no$, type_field) == no$) then
    nl=nl+1; write (li(nl), '(a, i5)') 'Number of Grid_field modes:', size(ele%grid_field)
  else
    nl2 = 10; if (type_field == all$) nl2 = 999
    nl=nl+1; li(nl) = ''
    if (ele%field_calc == bmad_standard$) then
      nl=nl+1; li(nl) = 'Grid_field: [NOT USED SINCE FIELD_CALC = BMAD_STANDARD]'
    else
      nl=nl+1; li(nl) = 'Grid_field:'
    endif
    do im = 1, size(ele%grid_field)
      g_field => ele%grid_field(im)
      if (g_field%master_parameter == 0) then
        name = '<None>'
      else
        name = attribute_name(ele, g_field%master_parameter)
      endif

      nl=nl+1; write (li(nl), '(a, i0)')      '  Grid_field mode #', im
      nl=nl+1; write (li(nl), '(2a)')         '    From file:          ', trim(g_field%ptr%file)
      nl=nl+1; write (li(nl), '(2a)')         '    field_type:         ', em_field_type_name(g_field%field_type)
      nl=nl+1; write (li(nl), '(2a)')         '    geometry:           ', grid_field_geometry_name(g_field%geometry)
      nl=nl+1; write (li(nl), '(2a)')         '    master_parameter:   ', trim(name)
      nl=nl+1; write (li(nl), '(2a)')         '    ele_anchor_pt:      ', anchor_pt_name(g_field%ele_anchor_pt)
      nl=nl+1; write (li(nl), '(a, i0)')      '    harmonic:           ', g_field%harmonic
      nl=nl+1; write (li(nl), '(a, i0)')      '    interpolation_order ', g_field%interpolation_order
      nl=nl+1; write (li(nl), '(a, es16.8)')  '    field_scale:        ', g_field%field_scale
      nl=nl+1; write (li(nl), '(a, es16.8)')  '    phi0_fieldmap:      ', g_field%phi0_fieldmap
      nl=nl+1; write (li(nl), '(a, l1)')      '    curved_ref_frame    ', g_field%curved_ref_frame
      nl=nl+1; write (li(nl), '(a, i0)')      '    n_link:             ', g_field%ptr%n_link
      nl=nl+1; write (li(nl), '(a, 3f14.6)')  '    r0:                 ', g_field%r0
      if (g_field%geometry == rotationally_symmetric_rz$) then
        lb = [lbound(g_field%ptr%pt,1), lbound(g_field%ptr%pt,2)]
        ub = [ubound(g_field%ptr%pt,1), ubound(g_field%ptr%pt,2)]
        nl=nl+1; write (li(nl), '(a, 3f14.6)')  '    dr:                 ', g_field%dr(1:2)
        nl=nl+1; write (li(nl), '(a, 3i14)')    '    Index_max:          ', ub
        nl=nl+1; write (li(nl), '(a, 3i14)')    '    Index_min:          ', lb
        nl=nl+1; write (li(nl), '(a, 3f14.6)')  '    r_max:              ', ub*g_field%dr(1:2) + g_field%r0(1:2)
        nl=nl+1; write (li(nl), '(a, 3f14.6)')  '    r_min:              ', lb*g_field%dr(1:2) + g_field%r0(1:2)
      else
        nl=nl+1; write (li(nl), '(a, 3f14.6)')  '    dr:                 ', g_field%dr
        nl=nl+1; write (li(nl), '(a, 3i14)')    '    Index_max:          ', ubound(g_field%ptr%pt)
        nl=nl+1; write (li(nl), '(a, 3i14)')    '    Index_min:          ', lbound(g_field%ptr%pt)
        nl=nl+1; write (li(nl), '(a, 3f14.6)')  '    r_max:              ', ubound(g_field%ptr%pt)*g_field%dr + g_field%r0
        nl=nl+1; write (li(nl), '(a, 3f14.6)')  '    r_min:              ', lbound(g_field%ptr%pt)*g_field%dr + g_field%r0
      endif
    enddo

    j = 0
    do id1 = lbound(g_field%ptr%pt, 1), ubound(g_field%ptr%pt, 1)          
    do id2 = lbound(g_field%ptr%pt, 2), ubound(g_field%ptr%pt, 2)
    do id3 = lbound(g_field%ptr%pt, 3), ubound(g_field%ptr%pt, 3)
      j = j + 1
      if (j > nl2) exit
      if (nl+1 > size(li)) call re_allocate(li, 2 * nl, .false.)
      nl=nl + 1

      select case (grid_field_dimension(g_field%geometry))
      case (1)
        write (line, '(2x, a, i0, 13a)') 'pt(', id1, ') = ('
      case (2)
        write (line, '(2x, 2(a, i0), 13a)') 'pt(', id1, ',', id2, ') = ('
      case (3)
        write (line, '(2x, 3(a, i0), 13a)') 'pt(', id1, ',', id2, ',', id3, ') = ('
      end select

      select case (g_field%field_type)
      case (mixed$)
        write (li(nl), '(2x, a, 13a)') trim(line), &
                                           trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(1))), ',', &
                                           trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(2))), ',', &
                                           trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(3))), ',', &
                                           trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(1))), ',', &
                                           trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(2))), ',', &
                                           trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(3))), ')'
      case (electric$)
        write (li(nl), '(2x, a, 13a)') trim(line), &
                                           trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(1))), ',', &
                                           trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(2))), ',', &
                                           trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%E(3))), ')'
      case (magnetic$)
        write (li(nl), '(2x, a, 13a)') trim(line), &
                                           trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(1))), ',', &
                                           trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(2))), ',', &
                                           trim(cmplx_re_str(g_field%ptr%pt(id1,id2,id3)%B(3))), ')'
      end select
    enddo
    enddo
    enddo

    if (size(g_field%ptr%pt) > nl2) then
      nl=nl+1; write (li(nl), '(a, i0, a)') '     .... etc ... (#Terms = ', size(g_field%ptr%pt), ')' 
    endif
  endif
endif

! ac_kick

if (associated(ele%ac_kick)) then
  ac => ele%ac_kick
  nl=nl+1; li(nl) = ''

  if (allocated(ac%amp_vs_time)) then
    nl=nl+1; li(nl) = 'AC_kicker (time, amplitude) knot points'
    nl=nl+1; li(nl) = '     Indx      Time       Amplitude'    
    do im = 1, size(ac%amp_vs_time)
      nl=nl+1; write (li(nl), '(i9, 2es14.6)') im, ac%amp_vs_time(im)%time, ac%amp_vs_time(im)%amp
    enddo
  endif

  if (allocated(ac%frequency)) then
    nl=nl+1; li(nl) = 'AC_kicker frequency components'
    nl=nl+1; li(nl) = '     Indx          Freq     Amplitude           Phi    Harmonic_Num'    
    do im = 1, size(ac%frequency)
      n = branch%n_ele_track
      nl=nl+1; write (li(nl), '(i9, 3es14.6, a16)') im, ac%frequency(im)%f, ac%frequency(im)%amp, &
            ac%frequency(im)%phi, real_to_string(ac%frequency(im)%f * branch%ele(n)%ref_time, 16, 6, 6)
    enddo
  endif
endif

! Radiation kick values

if (logic_option(.false., type_rad_kick) .and. associated(ele%rad_map)) then
  rm0 => ele%rad_map%rm0
  rm1 => ele%rad_map%rm1

  nl=nl+1; li(nl) = ''
  nl=nl+1; li(nl) = 'Matrices used for radiation stochastic and damping kicks:'
  nl=nl+1; write (li(nl), '(23x, a, 41x, a)') '1st half of element', '2nd half of element'
  nl=nl+1; li(nl) = 'Damping mat:'
  do i = 1, 6
    nl=nl+1; write (li(nl), '(2x, 6es10.2, 5x, 6es10.2)') rm0%damp_dmat(i,:), rm1%damp_dmat(i,:)
  enddo
  nl=nl+1; li(nl) = 'Damping vec:'
  nl=nl+1; write (li(nl), '(2x, 6es10.2, 5x, 6es10.2)') rm0%xfer_damp_vec(:), rm1%xfer_damp_vec(:)
  nl=nl+1; li(nl) = 'Damping ref_orb:'
  nl=nl+1; write (li(nl), '(2x, 6es10.2, 5x, 6es10.2)') rm0%ref_orb(:), rm1%ref_orb(:)
  nl=nl+1; li(nl) = 'Stochastic mat:'
  do i = 1, 6
    nl=nl+1; write (li(nl), '(2x, 6es10.2, 5x, 6es10.2)') rm0%stoc_mat(i,:), rm1%stoc_mat(i,:)
  enddo
endif

! wall3d cross-sections.
! Do not print more than 100 sections.

if (associated(ele%wall3d)) then
  do iw = 1, size(ele%wall3d)
    wall3d => ele%wall3d(iw)
    nl=nl+1; li(nl) = ''
    nl=nl+1; write (li(nl), '(2a)') 'Wall name: ', trim(wall3d%name)
    nl=nl+1; write (li(nl), '(a, i5)') 'Number of Wall Sections:', size(wall3d%section)
    nl=nl+1; write (li(nl), '(a, 2f11.5)') 'Wall region:',  wall3d%section(1)%s, wall3d%section(size(wall3d%section))%s
    if (logic_option(.false., type_wall)) then
      nl=nl+1; write (li(nl), '(2a)') 'Wall%ele_anchor_pt = ', anchor_pt_name(wall3d%ele_anchor_pt)
      select case (ele%key)
      case (capillary$)
      case (diffraction_plate$, mask$)
        nl=nl+1; write (li(nl), '(a, f10.6)') 'Wall%thickness       = ', wall3d%thickness
        nl=nl+1; write (li(nl), '(3a)')       'Wall%clear_material  = ', quote(wall3d%clear_material)
        nl=nl+1; write (li(nl), '(3a)')       'Wall%opaque_material = ', quote(wall3d%opaque_material)
      case default
        nl=nl+1; write (li(nl), '(a, l)') 'Wall%superimpose     = ', wall3d%superimpose
      end select
      n = min(size(wall3d%section), 100)
      do im = 1, n
        call re_allocate (li, nl+100, .false.) 
        section => wall3d%section(im)
        if (section%dr_ds == real_garbage$) then
          write (str1, '(a)')        ',  dr_ds = Not-set'
        else
          write (str1, '(a, f10.6)') ',  dr_ds =', section%dr_ds
        endif
        str2 = ''
        if (ele%key /= capillary$) then
          write (str2, '(2a)') ',   Type = ', trim(wall3d_section_type_name(section%type))
        endif

        nl=nl+1; write (li(nl), '(a, i0, a, f10.6, 2a, 2(f11.6, a), a)') &
                    'Wall%Section(', im, '):  S =', section%s, trim(str1),  ',   r0 = (', &
                    section%r0(1), ',', section%r0(2), ')', trim(str2)


        do j = 1, size(section%v)
          v => section%v(j)
          nl=nl+1; write (li(nl), '(4x, a, i0, a, 5f11.6)') &
                                'v(', j, ') =', v%x, v%y, v%radius_x, v%radius_y, v%tilt
        enddo
      enddo
    endif
  enddo
elseif (logic_option(.false., type_wall)) then
  nl=nl+1; write (li(nl), '(a)') 'No associated Wall.'
endif

! surface info

ph => ele%photon
if (associated(ph)) then
 
  nl=nl+1; li(nl) = ''
  if (ph%curvature%has_curvature) then
    nl=nl+1; li(nl) = 'Curviture:'
    nl=nl+1; write (li(nl), '(4x, a, f11.6)')  'Spherical        = ', ph%curvature%spherical
    nl=nl+1; write (li(nl), '(4x, a, 3f11.6)') 'Elliptical_X/Y/Z = ', ph%curvature%elliptical
    do ix = 0, ubound(ph%curvature%xy, 1)
    do iy = 0, ubound(ph%curvature%xy, 2)
      if (ph%curvature%xy(ix,iy) == 0) cycle
      nl=nl+1; write (li(nl), '(4x, 2(2x, 2(a, i0), a, es14.6))') 'X', ix, 'Y', iy, '   =', ph%curvature%xy(ix,iy)
    enddo
    enddo
  else
    nl=nl+1; li(nl) = 'No Curvature'
  endif

  seg => ph%segmented
  if (allocated(seg%pt)) then
    lb = lbound(seg%pt)
    ub = ubound(seg%pt)
    nl=nl+1; li(nl) = ''
    nl=nl+1; li(nl) = 'Surface Segmented'
    nl=nl+1; write (li(nl), '(4x, a, l1)')          'active:       ', seg%active
    nl=nl+1; write (li(nl), '(4x, a, 2f10.6)')      'dr:           ', seg%dr
    nl=nl+1; write (li(nl), '(4x, a, 2f10.6)')      'r0:           ', seg%r0
    nl=nl+1; write (li(nl), '(4x, a, 5(a, f10.6))') 'Bounds:       ', &
                        '(', seg%r0(1) + lb(1)*seg%dr(1), ',', seg%r0(1) + ub(1)*seg%dr(1), ')   (', & 
                             seg%r0(2) + lb(2)*seg%dr(2), ',', seg%r0(2) + ub(2)*seg%dr(2), ')'
    nl=nl+1; write (li(nl), '(4x, a, 5(a, i4))')    'Index Bounds: ', '(', lb(1), ',', ub(1), ')   (', lb(2), ',', ub(2), ')'
  endif

  disp => ph%displacement
  if (allocated(disp%pt)) then
    lb = lbound(disp%pt)
    ub = ubound(disp%pt)
    nl=nl+1; li(nl) = ''
    nl=nl+1; li(nl) = 'Surface Displacement'
    nl=nl+1; write (li(nl), '(4x, a, l1)')          'active:       ', disp%active
    nl=nl+1; write (li(nl), '(4x, a, 2f10.6)')      'dr:           ', disp%dr
    nl=nl+1; write (li(nl), '(4x, a, 2f10.6)')      'r0:           ', disp%r0
    nl=nl+1; write (li(nl), '(4x, a, 5(a, f10.6))') 'Bounds:       ', &
                        '(', disp%r0(1) + lb(1)*disp%dr(1), ',', disp%r0(1) + ub(1)*disp%dr(1), ')   (', & 
                             disp%r0(2) + lb(2)*disp%dr(2), ',', disp%r0(2) + ub(2)*disp%dr(2), ')'
    nl=nl+1; write (li(nl), '(4x, a, 5(a, i4))')    'Index Bounds: ', '(', lb(1), ',', ub(1), ')   (', lb(2), ',', ub(2), ')'
  endif


  hmis => ph%h_misalign
  if (allocated(hmis%pt)) then
    lb = lbound(hmis%pt)
    ub = ubound(hmis%pt)
    nl=nl+1; li(nl) = ''
    nl=nl+1; li(nl) = 'Surface H_misalign'
    nl=nl+1; write (li(nl), '(4x, a, l1)')          'active:       ', hmis%active
    nl=nl+1; write (li(nl), '(4x, a, 2f10.6)')      'dr:           ', hmis%dr
    nl=nl+1; write (li(nl), '(4x, a, 2f10.6)')      'r0:           ', hmis%r0
    nl=nl+1; write (li(nl), '(4x, a, 5(a, f10.6))') 'Bounds:       ', &
                        '(', hmis%r0(1) + lb(1)*hmis%dr(1), ',', hmis%r0(1) + ub(1)*hmis%dr(1), ')   (', & 
                             hmis%r0(2) + lb(2)*hmis%dr(2), ',', hmis%r0(2) + ub(2)*hmis%dr(2), ')'
    nl=nl+1; write (li(nl), '(4x, a, 5(a, i4))')    'Index Bounds: ', '(', lb(1), ',', ub(1), ')   (', lb(2), ',', ub(2), ')'
  endif

  pixel => ph%pixel
  if (allocated(pixel%pt)) then
    lb = lbound(pixel%pt)
    ub = ubound(pixel%pt)
    nl=nl+1; li(nl) = ''
    nl=nl+1; li(nl) = 'Surface Pixel'
    nl=nl+1; write (li(nl), '(4x, a, 2f10.6)')      'dr:           ', pixel%dr
    nl=nl+1; write (li(nl), '(4x, a, 2f10.6)')      'r0:           ', pixel%r0
    nl=nl+1; write (li(nl), '(4x, a, 5(a, f10.6))') 'Bounds:       ', &
                        '(', pixel%r0(1) + lb(1)*pixel%dr(1), ',', pixel%r0(1) + ub(1)*pixel%dr(1), ')   (', & 
                             pixel%r0(2) + lb(2)*pixel%dr(2), ',', pixel%r0(2) + ub(2)*pixel%dr(2), ')'
    nl=nl+1; write (li(nl), '(4x, a, 5(a, i4))')    'Index Bounds: ', '(', lb(1), ',', ub(1), ')   (', lb(2), ',', ub(2), ')'
  endif

  if (ph%material%f_h /= 0) then
    nl=nl+1; li(nl) = ''
    nl=nl+1; li(nl) = 'Structure Factors (Re, Im):'
    if (ele%key == multilayer_mirror$) then
      nl=nl+1; write (li(nl), '(4x, 2(a,f10.4))') 'F_0 (Material 1):', real(ph%material%f0_m1), ',', aimag(ph%material%f0_m1)
      nl=nl+1; write (li(nl), '(4x, 2(a,f10.4))') 'F_0 (Material 2):', real(ph%material%f0_m2), ', ', aimag(ph%material%f0_m2)
    else
      nl=nl+1; write (li(nl), '(4x, 2(a,f10.4))') 'F_0:             ', real(ph%material%f_0), ', ', aimag(ph%material%f_0)
    endif
    nl=nl+1; write (li(nl), '(4x, 2(a,f10.4))') 'F_H:             ', real(ph%material%f_h), ', ', aimag(ph%material%f_h)
    nl=nl+1; write (li(nl), '(4x, 2(a,f10.4))') 'F_Hbar:          ', real(ph%material%f_hbar), ', ', aimag(ph%material%f_hbar)
    nl=nl+1; write (li(nl), '(4x, 2(a,f10.4))') 'Sqrt(F_H*F_Hbar):', real(ph%material%f_hkl), ', ', aimag(ph%material%f_hkl)
  endif

  if (allocated(ph%init_energy_prob)) then
    nl=nl+1; li(nl) = ''
    nl=nl+1; li(nl) = 'Energy_probability_curve:'
    nl=nl+1; li(nl) = '   Energy(eV)      Prob_density/eV (normalized)'
    n = size(ph%init_energy_prob)
    do i = 1, min(5, n)
      nl=nl+1; write(li(nl), '(4x, 2es12.4)') ph%init_energy_prob(i)%x0, ph%init_energy_prob(i)%y0 / ph%integrated_init_energy_prob(n)
    enddo
    if (n > 10) then
      nl=nl+1; write(li(nl), '(3a)') '   ..... ', int_str(n), ' terms total ...'
    endif
    do i = 1, min(5, n-5)
      j = i + n - 5
      nl=nl+1; write(li(nl), '(4x, 2es12.4)') ph%init_energy_prob(j)%x0, ph%init_energy_prob(j)%y0 / ph%integrated_init_energy_prob(n)
    enddo
  endif

  rt => ph%reflectivity_table_sigma
  if (allocated(rt%p_reflect)) then
    na = size(rt%angle)
    ne = size(rt%energy)
    nl=nl+1; li(nl) = ''
    select case (ph%reflectivity_table_type)
    case (polarized$);  nl=nl+1; li(nl) = 'Reflectivity table for Sigma polarization:'
    case default;       nl=nl+1; li(nl) = 'Reflectivity table for both polarizations:'
    end select
    nl=nl+1; write(li(nl), '(a, 2f14.9, i6)') '  Angles Min, Max, #Points:           ', rt%angle(1),  rt%angle(na),  na
    nl=nl+1; write(li(nl), '(a, 2f14.1, i6)') '  Energy Min, Max, #Points:           ', rt%energy(1), rt%energy(ne), ne
    if (ele%key == crystal$) then
      nl=nl+1; write(li(nl), '(a, 2f14.9)')   '  Bragg angle (uncorrected) Min, Max: ', rt%bragg_angle(1), rt%bragg_angle(ne)
    endif
  endif

  rt => ph%reflectivity_table_pi
  if (allocated(rt%p_reflect)) then
    na = size(rt%angle)
    ne = size(rt%energy)
    nl=nl+1; li(nl) = ''
    nl=nl+1; li(nl) = 'Reflectivity table for Pi polarization:'
    nl=nl+1; write(li(nl), '(a, 2f14.9, i6)') '  Angles Min, Max, #Points:           ', rt%angle(1),  rt%angle(na),  na
    nl=nl+1; write(li(nl), '(a, 2f14.1, i6)') '  Energy Min, Max, #Points:           ', rt%energy(1), rt%energy(ne), ne
    if (ele%key == crystal$) then
      nl=nl+1; write(li(nl), '(a, 2f14.9)')   '  Bragg angle (uncorrected) Min, Max: ', rt%bragg_angle(1), rt%bragg_angle(ne)
    endif
  endif
endif



! Encode branch info

if (ele%key == fork$ .or. ele%key == photon_fork$) then
  if (li(nl) /= '') then
    nl=nl+1; li(nl) = ' '
  endif

  n = nint(ele%value(ix_to_branch$))
  im = nint(ele%value(ix_to_element$))
  if (associated(lat)) then
    nl=nl+1; write (li(nl), '(5a, 2(i0, a))') 'Branch to: ', trim(lat%branch(n)%name), '>>', &
                                  trim(lat%branch(n)%ele(im)%name), '  [', n, '>>', im, ']'
  else
    nl=nl+1; write (li(nl), '(a, i0, a, i0)') 'Branch to: ', n, '>>', im
  endif
endif

! Encode lord/slave info.
! For super_lords there is no attribute_name associated with a slave.
! For slaves who are overlay_lords then the attribute_name is obtained by
!   looking at the overlay_lord's 1st slave (slave of slave of the input ele).

if (associated(lat) .and. integer_option(short$, type_control) /= no$) then
  ! Print info on element's lords

  if (li(nl) /= '') then
    nl=nl+1; li(nl) = ' '
  endif

  if (ele%slave_status <= 0) then
    nl=nl+1; write (li(nl), '(a)') 'Slave_status: BAD! PLEASE SEEK HELP!', ele%slave_status
  else
    nl=nl+1; write (li(nl), '(2a)') 'Slave_status: ', control_name(ele%slave_status)
  endif

  if (nL + ele%n_lord + 100 > size(li)) call re_allocate(li, nl + ele%n_lord + 100)

  select case (ele%slave_status)
  case (multipass_slave$)
    lord => pointer_to_lord(ele, 1)
    nl=nl+1; write (li(nl), '(3a, i0, a)') 'Associated Multipass_Lord: ', trim(lord%name), '  (Index: ', lord%ix_ele, ')'
    nl=nl+1; li(nl) = 'Other slaves of this Lord:'
    nl=nl+1; li(nl) = '     Index   Name'
    do im = 1, lord%n_slave
      slave => pointer_to_slave(lord, im)
      if (slave%ix_ele == ele%ix_ele .and. slave%ix_branch == ele%ix_branch) cycle
      nl=nl+1; write (li(nl), '(a, 3x, a)') adjustr(ele_loc_name(slave)), trim(slave%name)
    enddo

  case (super_slave$)
    nl=nl+1; write (li(nl), '(3a, i0, a)') 'Associated Super_Lord(s):'
    nl=nl+1; li(nl) = '   Index   Name                             Type'
    do im = 1, ele%n_lord
      lord => pointer_to_lord(ele, im)
      if (lord%lord_status /= super_lord$) cycle
      if (lord%slave_status == multipass_slave$) then
        lord2 => pointer_to_lord(lord, 1)
        nl=nl+1; write (li(nl), '(i8, 3x, a, t45, 3a, 2x, a)') lord%ix_ele, trim(lord%name), trim(key_name(lord%key)), &
                      '   --> Multipass_slave of: ', trim(ele_loc_name(lord2)), lord2%name
      else
        nl=nl+1; write (li(nl), '(i8, 3x, a, t45, a)') lord%ix_ele, trim(lord%name), trim(key_name(lord%key))
      endif
    enddo
  end select

  ! Print controller lords

  has_it = .false.
  do im = 1, ele%n_lord
    lord => pointer_to_lord(ele, im)
    if (lord%lord_status == multipass_lord$ .or. lord%lord_status == super_lord$) cycle
    has_it = .true.
    exit
  enddo

  if (has_it) then

    nn = 10
    na = 12
    nt = 12

    do im = 1, ele%n_lord
      lord => pointer_to_lord (ele, im, ctl)
      nn = max(nn, len_trim(lord%name)+2)
      select case (lord%lord_status)
      case (super_lord$, multipass_lord$, girder_lord$)
      case default
        na = max(na, len_trim(ctl%attribute)+2)
        if (.not. lord%is_on) nt = 18
      end select
    enddo

    str1 = ''
    nl=nl+1; li(nl) = 'Controller Lord(s):'
    nl=nl+1; write (li(nl), '(9a)') '   Index   Name',  str1(:nn-4), 'Attribute', str1(:na-9), 'Lord_Type', str1(:nt-9), 'Expression'

    do im = 1, ele%n_lord
      lord => pointer_to_lord (ele, im, ctl)
      str1 = key_name(lord%key)
      select case (lord%lord_status)
      case (super_lord$, multipass_lord$)
        cycle
      case (control_lord$)
        call re_allocate (li2, 1)
        li2(1) = trim(ctl%attribute) // ' Slave'
        a_name = ''
      case (girder_lord$)
        call re_allocate (li2, 1)
        li2(1) = ''
        a_name = ''
      case default
        if (allocated(ctl%stack)) then
          call split_expression_string (expression_stack_to_string(ctl%stack), 70, 5, li2)
        else   ! Knots
          call split_expression_string (knots_to_string(lord%control%x_knot, ctl%y_knot), 70, 5, li2)
        endif
        a_name = ctl%attribute
        if (.not. lord%is_on) str1 = trim(str1) // ' [Is off]'
      end select

      if (nl+size(li2)+100 > size(li)) call re_allocate (li, nl+size(li2)+100)
      nl=nl+1; write (li(nl), '(i8, 3x, 4a)') lord%ix_ele, lord%name(:nn), a_name(:na), str1(:nt), trim(li2(1))
      do j = 2, size(li2)
        nl=nl+1; li(nl) = ''; li(nl)(17+nn+na+nt:) = trim(li2(j))
      enddo
    enddo
  endif

  !

  has_it = .false.
  nl=nl+1; li(nl) = ' '
  nl=nl+1; li(nl) = 'Elements whose fields overlap this one:'
  nl=nl+1; li(nl) = '   Index   Name                               Type'

  if (ele%slave_status == super_slave$ .or. ele%slave_status == multipass_slave$) then
    do im = 1, ele%n_lord
      lord => pointer_to_lord(ele, im)
      if (lord%slave_status == multipass_slave$) lord => pointer_to_lord(lord, 1)
      do j = 1, lord%n_lord_field
        has_it = .true.
        lord2 => pointer_to_lord(lord, lord%n_lord+im)
        nl=nl+1; write (li(nl), '(a8, t12, a35, a16, f10.3)') &
                      trim(ele_loc_name(lord2)), lord2%name, key_name(lord2%key)
      enddo
    enddo
  else
    do im = 1, ele%n_lord_field
      has_it = .true.
      lord => pointer_to_lord(ele, ele%n_lord+im)
      nl=nl+1; write (li(nl), '(a8, t12, a35, a16, f10.3)') &
                    trim(ele_loc_name(lord)), lord%name, key_name(lord%key)

    enddo
  endif

  if (.not. has_it) nl = nl - 3

  !

  if (ele%n_lord_ramper /= 0) then
    nl=nl+1; li(nl) = ' '
    nl=nl+1; li(nl) = 'Ramper Lords:'
    nl=nl+1; li(nl) = '   Ramper_Name               Attribute           Value                Expression/Knot Points'
    do ix = 1, ele%n_lord_ramper
      lord => pointer_to_lord(ele, ix, lord_type = ramper_lord$, ix_control = ic)
      rmp => lord%control%ramp(ic)
      if (allocated(rmp%stack)) then
        call split_expression_string (expression_stack_to_string(rmp%stack), 80, 5, li2)
      else  ! Spline
        call split_expression_string (knots_to_string(lord%control%x_knot, rmp%y_knot), 80, 5, li2, '),')
      endif
      value = ramper_value(lord, rmp, err_flag)

      nl=nl+1; write (li(nl), '(3x, a, t30, a18, es20.12, 4x, a)') ele_full_name(lord), rmp%attribute, value, trim(li2(1))
      if (nl+size(li2)+100 > size(li)) call re_allocate (li, nl+size(li2)+100)
      do im = 2, size(li2)
        n = 70
        if (im == 3 .and. integer_option(short$, type_control) /= all$) then
          nl=nl+1; li(nl) = ''; li(nl)(n:) = trim(li2(im)) // ' ... etc.'
          exit
        else
          nl=nl+1; li(nl) = ''; li(nl)(n:) = trim(li2(im))
        endif
      enddo
    enddo
  endif

  ! Print info on elements slaves.

  nl=nl+1; li(nl) = ' '

  if (ele%lord_status <= 0) then
    nl=nl+1; write (li(nl), '(a)') 'Lord_status: BAD!', ele%lord_status
  else
    nl=nl+1; write (li(nl), '(2a)') 'Lord_status:  ', control_name(ele%lord_status)
  endif

  if (associated(ele%control)) then
    nl=nl+1; li(nl) = 'Control Variables:'

    if (ele%lord_status == group_lord$) then
      n_att = maxval(len_trim(ele%control%var%name))
      do im = 1, size(ele%control%var)
        a_name = ele%control%var(im)%name
        nl=nl+1; write (li(nl), '(i5, 3x, 2a, es15.7, 11x, 3a, es15.7)')  im, &
                      a_name(1:n_att), '  =', ele%control%var(im)%value, &
                      'OLD_', a_name(1:n_att), '  =', ele%control%var(im)%old_value
      enddo

    elseif (allocated(ele%control%var)) then  ! overlay_lord or ramper_lord. A ramper slave does not have this allocated.
      do im = 1, size(ele%control%var)
        nl=nl+1; write (li(nl), '(i5, 3x, 2a, es15.7)')  im, ele%control%var(im)%name, '  =', ele%control%var(im)%value
      enddo
    endif

    ! Print named constants if present

    print_it = .true.
    do is = 1, ele%n_slave
      slave => pointer_to_slave (ele, is, ctl)
      if (.not. allocated(ctl%stack)) cycle
      call print_this_stack(nl, li, ctl%stack, print_it, str_index)
    enddo
  endif

  !

  if (ele%n_slave /= 0 .or. ele%key == ramper$) then
    if (nl + ele%n_slave + 100 > size(li)) call re_allocate(li, nl + ele%n_slave + 100)

    n_char = 10
    do im = 1, ele%n_slave
      slave => pointer_to_slave (ele, im)
      n_char = max(n_char, len_trim(slave%name))
    enddo

    select case (ele%lord_status)

    case (multipass_lord$, super_lord$, girder_lord$, control_lord$)
      if (ele%key == feedback$) then
        nl=nl+1; write (li(nl), '(a, i4)') 'Slaves:'
        nl=nl+1; li(nl) = '   Index   Name';  li(nl)(n_char+14:) = 'Ele_Type           Slave_Type      S'
        do im = 1, ele%n_slave
          slave => pointer_to_slave (ele, im, ctl)
          nl=nl+1; write (li(nl), '(a8, t12, a, 2x, a16, 3x, a6, f14.6)') trim(ele_loc_name(slave)), &
                                      slave%name(1:n_char), key_name(slave%key), ctl%attribute, slave%s
        enddo

      else
        nl=nl+1; write (li(nl), '(a, i4)') 'Slaves:'
        nl=nl+1; li(nl) = '   Index   Name';  li(nl)(n_char+14:) = 'Type                     S'
        do im = 1, ele%n_slave
          slave => pointer_to_slave (ele, im)
          nl=nl+1; write (li(nl), '(a8, t12, a, 2x, a16, 3x, f14.6)') &
                      trim(ele_loc_name(slave)), slave%name(1:n_char), key_name(slave%key), slave%s
        enddo
      endif

    case (ramper_lord$)
      print_it = .true.
      do ix = 1, size(ele%control%ramp)
        rmp => ele%control%ramp(ix)
        if (.not. allocated(rmp%stack)) cycle
        call print_this_stack(nl, li, rmp%stack, print_it, str_index)
      enddo

      nl=nl+1; write (li(nl), '(a, i4)') 'Slaves:'
      nl=nl+1; li(nl) = '   Ele_Name            Attribute           Value                  Expression/Knot Points'
      do ix = 1, size(ele%control%ramp)
        rmp => ele%control%ramp(ix)

        if (allocated(rmp%stack)) then
          call split_expression_string (expression_stack_to_string(rmp%stack), 80, 5, li2)
        else  ! Spline
          call split_expression_string (knots_to_string(ele%control%x_knot, rmp%y_knot), 80, 5, li2, '),')
        endif
        value = ramper_value(ele, rmp, err_flag)

        nl=nl+1; write (li(nl), '(3x, a20, a18, es20.12, 4x, a)') rmp%slave_name, rmp%attribute, value, trim(li2(1))
        if (nl+size(li2)+100 > size(li)) call re_allocate (li, nl+size(li2)+100)
        do im = 2, size(li2)
          n = 65
          nl=nl+1; li(nl) = ''; li(nl)(n:) = trim(li2(im))
        enddo
      enddo

    case default
      if (ele%key == overlay$) then
        nl=nl+1; write (li(nl), '(a, i4)') 'Slaves: [Attrib_Value = Expression_Val summed over all Overlays controlling the attribute.]'
      else ! Group
        nl=nl+1; write (li(nl), '(a, i4)') 'Slaves: [Attrib_Value = Value of the controlled attribute, Expression_Val = Value calculated by this Group element.]'
      endif

      nl=nl+1; li(nl) = ' Ele_Loc   Ele_Name';  li(nl)(n_char+14:) = 'Attribute         Attrib_Value  Expression_Val     Expression/Knot Points'
      do ix = 1, ele%n_slave
        slave => pointer_to_slave (ele, ix, ctl)

        select case (ctl%attribute)
        case ('ACCORDION_EDGE', 'START_EDGE', 'END_EDGE', 'S_POSITION')  ! Special group constructs
          attrib_val_str = ' ----'
        case default
          call pointer_to_attribute (slave, ctl%attribute, .false., a_ptr, err_flag)
          attrib_val_str = ' ----'
          if (associated(a_ptr%r)) write (attrib_val_str, '(es12.4)') a_ptr%r
        end select

        write (str1, '(es12.4)') ctl%value

        if (allocated(ctl%stack)) then
          call split_expression_string (expression_stack_to_string(ctl%stack), 70, 5, li2)
        else  ! Spline
          call split_expression_string (knots_to_string(ele%control%x_knot, ctl%y_knot), 70, 5, li2)
        endif

        nl=nl+1; write (li(nl), '(a8, t12, a, 2x, a18, 2a, 4x, a)') trim(ele_loc_name(slave)), &
                                  slave%name(1:n_char), ctl%attribute, attrib_val_str, str1(1:17), trim(li2(1))
        if (nl+size(li2)+100 > size(li)) call re_allocate (li, nl+size(li2)+100)
        do im = 2, size(li2)
          n = 50 + n_char + len(attrib_val_str)
          nl=nl+1; li(nl) = ''; li(nl)(n:) = trim(li2(im))
        enddo

      enddo
    end select
  endif

  ! Print elements that are field overlapped.

  has_it = .false.
  nl=nl+1; li(nl) = ' '
  nl=nl+1; li(nl) = "This element's field overlaps:"
  nl=nl+1; li(nl) = '   Index   Name                                      Type '

  if (ele%slave_status == super_slave$ .or. ele%slave_status == multipass_slave$) then
    do im = 1, ele%n_lord
      lord => pointer_to_lord(ele, im)
      if (lord%slave_status == multipass_slave$) lord => pointer_to_lord(lord, 1)
      do j = 1, lord%n_slave_field
        has_it = .true.
        slave => pointer_to_slave(ele, im, lord_type = field_lord$)
        nl=nl+1; write (li(nl), '(a8, t12, a30, a16, f10.3)') &
                      trim(ele_loc_name(slave)), slave%name, trim(key_name(slave%key))
      enddo
    enddo
  else
    do im = 1, ele%n_slave_field
      has_it = .true.
      slave => pointer_to_slave(ele, im, lord_type = field_lord$)
      nl=nl+1; write (li(nl), '(a8, t12, a30, a16, f10.3)') &
                    trim(ele_loc_name(slave)), slave%name, trim(key_name(slave%key))
    enddo
  endif

  if (.not. has_it) nl=nl - 3
endif

! Encode Twiss info

if (associated(branch) .and. ele%lord_status == not_a_lord$) then
  if (.not. branch%param%live_branch) then
    nl=nl+1; li(nl) = ''
    nl=nl+1; li(nl) = 'NOTE: Branch containing element is not live (branch''s live_branch parameter has been set to False).'
  endif
endif

if (ele%key /= ramper$ .and. ele%key /= girder$) then
  if (integer_option(radians$, twiss_out) /= 0 .and. ele%a%beta /= 0) then
    nl=nl+1; li(nl) = ''
    nl=nl+1; li(nl) = 'Twiss at end of element:'
    call type_twiss (ele, twiss_out, .false., li(nl+1:), nl2)
    nl=nl + nl2
  else
    nl=nl+1; li(nl) = ''
    nl=nl+1; li(nl) = 'Twiss: Not computed at element'  
  endif
endif

l_status = ele%lord_status
if (l_status /= overlay_lord$ .and. l_status /= multipass_lord$ .and. &
    l_status /= group_lord$ .and. l_status /= girder_lord$ .and. l_status /= ramper_lord$) then

  ! Encode mat6 info

  n = integer_option (6, type_mat6)
  if (n > 6) n = 6
  if (n < 0) n = 0

  if (n /= 0) then
    nl=nl+1; li(nl) = ' '
    nl=nl+1; write (li(nl), '(a, es12.3, a, t82, a, t95, a, t108, a)') 'Transfer Matrix : Kick  [Mat symplectic error:', &
                 mat_symp_error(ele%mat6), ']', 'Vec0'
  endif

  if (any(abs(ele%mat6(1:n,1:n)) >= 1d3)) then
    write (fmt, '(a, i0, a)') '(', n, 'es12.4, a, es13.5, 2es15.5)'
  else
    write (fmt, '(a, i0, a)') '(', n, 'f12.7, a, es13.5, 2es15.5)'
  endif

  do im = 1, n
    nl=nl+1; write (li(nl), fmt) (ele%mat6(im, j), j = 1, n), '   : ', ele%vec0(im)
  enddo

  ! Encode taylor series

  if (associated(ele%taylor(1)%term)) then
    nl=nl+1; li(nl) = ' '
    nl=nl+1; write (li(nl), '(a, l1)') 'taylor_map_includes_offsets: ', ele%taylor_map_includes_offsets
    if (logic_option(.false., type_taylor)) then
      call type_taylors (ele%taylor, lines = li2, n_lines = nt)
      call re_allocate (li, nl+nt+100, .false.)
      li(1+nl:nt+nl) = li2(1:nt)
      deallocate (li2)
      nl = nl + nt
    else
      n_term = 0
      do im = 1, size(ele%taylor)
        n_term = n_term + size(ele%taylor(im)%term)
      enddo
      nl=nl+1; write (li(nl), '(a, i6)') 'Taylor map total number of terms:', n_term
    endif
  endif

  if (associated(ele%spin_taylor(0)%term) .or. ele%key == taylor$) then
    if (logic_option(.false., type_taylor)) then
      nl=nl+1; li(nl) = ''
      call type_taylors (ele%spin_taylor, lines = li2, n_lines = nt)
      call re_allocate (li, nl+nt+100, .false.)
      li(1+nl:nt+nl) = li2(1:nt)
      deallocate (li2)
      nl = nl + nt
    else
      n_term = 0
      do im = 0, 3
        n_term = n_term + size(ele%spin_taylor(im)%term)
      enddo
      nl=nl+1; write (li(nl), '(a, i6)') 'Spin_Taylor map total number of terms:', n_term
    endif
  endif

endif

! Print wake info

if (associated(ele%wake)) then

  if (logic_option (.true., type_wake) .and. (size(ele%wake%sr%long) /= 0 .or. &
                                              size(ele%wake%sr%trans) /= 0 .or.size(ele%wake%sr%z_long%w) /= 0 )) then
    nl=nl+1; li(nl) = ''
    nl=nl+1; li(nl) = 'Short-Range Wake:'
    if (ele%wake%sr%file /= '') then
      nl=nl+1; li(nl) = '  SR_File: ' // trim(ele%wake%sr%file)
    endif
    nl=nl+1; write (li(nl), '(2x, a, l2)') 'scale_with_length =', ele%wake%sr%scale_with_length
    nl=nl+1; write (li(nl), '(2x, 2a)')    'amp_scale         = ', to_str(ele%wake%sr%amp_scale)
    nl=nl+1; write (li(nl), '(2x, 2a)')    'z_scale           = ', to_str(ele%wake%sr%z_scale)
    nl=nl+1; write (li(nl), '(2x, 2a)')    'z_max             = ', to_str(ele%wake%sr%z_max)
  endif

  if (size(ele%wake%sr%long) /= 0) then
    nl=nl+1; write (li(nl), *)
    if (logic_option (.true., type_wake)) then
      call re_allocate (li, nl+size(ele%wake%sr%long)+100, .false.)
      nl=nl+1; li(nl) = '  Short-Range Longitudinal Pseudo Modes:'
      nl=nl+1; li(nl) = &
            '   #        Amp        Damp           K         Phi   Transverse_Dependence'
      do im = 1, size(ele%wake%sr%long)
        mode => ele%wake%sr%long(im)
        nl=nl+1; write (li(nl), '(i4, 4es12.4, a15, a16)') im, mode%amp, mode%damp, mode%k, mode%phi, &
                  sr_longitudinal_position_dep_name(mode%position_dependence)
      enddo
    else
      nl=nl+1; li(nl) = '  No short-range longitudinal pseudo modes.'
    endif
  endif

  if (size(ele%wake%sr%trans) /= 0) then
    nl=nl+1; write (li(nl), *)
    if (logic_option (.true., type_wake)) then
      call re_allocate (li, nl+size(ele%wake%sr%trans)+100, .false.)
      nl=nl+1; li(nl) = '  Short-Range Transverse Pseudo Modes:'
      nl=nl+1; li(nl) = &
            '   #        Amp        Damp           K         Phi   Polarization  Transverse_Dependence'
      do im = 1, size(ele%wake%sr%trans)
        mode => ele%wake%sr%trans(im)
        nl=nl+1; write (li(nl), '(i4, 4es12.4, a15, a16)') im, mode%amp, mode%damp, mode%k, mode%phi, &
                  sr_transverse_polarization_name(mode%polarization), sr_transverse_position_dep_name(mode%position_dependence)
      enddo
    else
     nl=nl+1; li(nl) = '  No short-range transverse pseudo modes.'
    endif
  endif

  if (size(ele%wake%sr%z_long%w) /= 0) then
    nl=nl+1; write (li(nl), *)
    if (logic_option (.true., type_wake)) then
      call re_allocate (li, nl+size(ele%wake%sr%z_long%w)+100, .false.)
      nl=nl+1; li(nl) = '  Short-Range Z-dependent Longitudinal wake:'
      srz => ele%wake%sr%z_long
      nl=nl+1; li(nl) = '    smoothing_sigma     = ' // to_str(srz%smoothing_sigma)
      nl=nl+1; li(nl) = '    position_dependence = ' // trim(sr_transverse_position_dep_name(srz%position_dependence))
      nl=nl+1; li(nl) = '    dz (scaled) [m]  = ' // to_str(srz%dz / ele%wake%sr%z_scale)
      nl=nl+1; li(nl) = '    Wake range (scaled) [m]: +/-' // to_str(srz%z0 / ele%wake%sr%z_scale)
      nl=nl+1; write (li(nl), '(a, i0)') '  # wake points: ', size(srz%w)

    else
     nl=nl+1; li(nl) = '  No short-range z-dependent modes.'
    endif
  endif

  if (logic_option (.true., type_wake) .and. size(ele%wake%lr%mode) /= 0) then
    nl=nl+1; li(nl) = ''
    nl=nl+1; li(nl) = 'Long-Range Wake:'
    if (ele%wake%lr%file /= '') then
      nl=nl+1; li(nl) = '  LR_File: ' // trim(ele%wake%lr%file)
    endif
    nl=nl+1; write (li(nl), '(2x, 2a)')    'amp_scale    = ', to_str(ele%wake%lr%amp_scale)
    nl=nl+1; write (li(nl), '(2x, 2a)')    'time_scale   = ', to_str(ele%wake%lr%time_scale)
    nl=nl+1; write (li(nl), '(2x, 2a)')    'freq_spread  = ', to_str(ele%wake%lr%freq_spread)
    nl=nl+1; write (li(nl), '(2x, a, l1)') 'self_wake_on = ', ele%wake%lr%self_wake_on
    nl=nl+1; write (li(nl), '(2x, 2a)')    't_ref        = ', to_str(ele%wake%lr%t_ref)
  endif

  if (size(ele%wake%lr%mode) /= 0) then
    nl=nl+1; write (li(nl), *)
    if (logic_option (.true., type_wake)) then
      call re_allocate (li, nl+size(ele%wake%lr%mode)+100, .false.)
      nl=nl+1; li(nl) = '  Long-Range Wake Modes [Note: Freq will not be equal to Freq_in if there is a frequency spread]:'
      nl=nl+1; li(nl) = &
            '  #     Freq_in        Freq         R/Q        Damp           Q        Phi   m   Angle    b_sin     b_cos     a_sin     a_cos'
      do im = 1, size(ele%wake%lr%mode)
        lr => ele%wake%lr%mode(im)
        angle = ' unpolar'
        if (lr%polarized) write (angle, '(f8.3)') lr%angle
        if (lr%damp == 0 .or. lr%freq <= 0) then
          q_factor = '      ------'
        else
          write (q_factor, '(es12.4)') pi * lr%freq / lr%damp
        endif


        nl=nl+1; write (li(nl), '(i3, 4es12.4, a, es12.4, i3, a, 5es10.2)') im, &
                lr%freq_in, lr%freq, lr%R_over_Q, lr%damp, q_factor, lr%phi, lr%m, angle, &
                lr%b_sin, lr%b_cos, lr%a_sin, lr%a_cos
      enddo
    else
      nl=nl+1; li(nl) = '  No long-range HOM modes.'
    endif
  endif

endif

! Encode Floor coords. 
! Elements not associated with a lattice do not have floor coords.

if (logic_option(.false., type_floor_coords) .and. associated(ele%branch)) then
  if (ele%lord_status == multipass_lord$) then
    nl=nl+1; write (li(nl), '(27x, a)')         'X           Y           Z       Theta         Phi         Psi'
    do ie = 1, ele%n_slave
      slave => pointer_to_slave(ele, ie)
      floor = slave%floor
      nl=nl+1; write (li(nl), '(a, i2, 6f12.5, 3x, a)') 'Multipass Slave', ie, floor%r, floor%theta, floor%phi, floor%psi, '! Position without misalignments'
    enddo

  !
  else
    ele0 => pointer_to_next_ele(ele, -1)

    select case (ele%key)
    case (crystal$, mirror$, multilayer_mirror$)
      call ele_geometry (ele0%floor, ele, floor2, 0.5_rp)
      floor = ele_geometry_with_misalignments (ele, 0.5_rp)

      nl=nl+1; li(nl) = ''
      nl=nl+1; li(nl) = 'Global Floor Coords at Surface of Element:'
      nl=nl+1; write (li(nl), '(a)')         '                   X           Y           Z       Theta         Phi         Psi'
      nl=nl+1; write (li(nl), '(a, 6f12.5, 3x, a)') 'Reference', floor2%r, floor2%theta, floor2%phi, floor2%psi, '! Position without misalignments'
      nl=nl+1; write (li(nl), '(a, 6f12.5, 3x, a)') 'Actual   ', floor%r, floor%theta, floor%phi, floor%psi, '! Position with offset/pitch/tilt misalignments'

    case (girder$)
      floor = ele_geometry_with_misalignments (ele)
      nl=nl+1; li(nl) = ''
      nl=nl+1; li(nl) = 'Global Floor Coords at Reference Point:'
      nl=nl+1; write (li(nl), '(a)')         '                   X           Y           Z       Theta         Phi         Psi'
      nl=nl+1; write (li(nl), '(a, 6f12.5, 3x, a)') 'Reference', ele%floor%r, ele%floor%theta, ele%floor%phi, ele%floor%psi, '! Position without misalignments'
      nl=nl+1; write (li(nl), '(a, 6f12.5, 3x, a)') 'Actual   ', floor%r, floor%theta, floor%phi, floor%psi, '! Position with offset/pitch/tilt misalignments'
    end select

    !

    select case (ele%key)
    case (girder$)
    case default
      floor = ele_geometry_with_misalignments (ele)

      nl=nl+1; li(nl) = ''
      nl=nl+1; li(nl) = 'Global Floor Coords at End of Element:'
      nl=nl+1; write (li(nl), '(a)')         '                   X           Y           Z       Theta         Phi         Psi'
      nl=nl+1; write (li(nl), '(a, 6f12.5, 3x, a)') 'Reference', ele%floor%r, ele%floor%theta, ele%floor%phi, ele%floor%psi, '! Position without misalignments'
      nl=nl+1; write (li(nl), '(a, 6f12.5, 3x, a)') 'Actual   ', floor%r, floor%theta, floor%phi, floor%psi, '! Position with offset/pitch/tilt misalignments'
    end select

    !

    if (associated(ele0) .and. (ele%ix_ele /= 0 .or. branch%param%geometry == closed$)) then
      f0 = ele0%floor
      nl=nl+1; write (li(nl), '(a, 6f12.5, 3x, a)') 'delta Ref', floor%r-f0%r, floor%theta-f0%theta, floor%phi-f0%phi, floor%psi-f0%psi, &
                                                                                                         '! Delta of reference with respect to last element'  
    endif
  endif
endif

! finish

if (present(lines)) then
  call re_allocate(lines, nl, .false.)
  n_lines = nl
  lines(1:nl) = li(1:nl)
else
  do im = 1, nl
    print '(1x, a)', trim(li(im))
  enddo
endif

!----------------------------------------------------------------------------------------------------
contains

subroutine encode_2nd_column_parameter (li, nl2, nl, attrib_name, re_val, str_val, logic_val, int_val)

integer nl2, nl, ix0
integer, optional :: int_val
real(rp), optional :: re_val
logical, optional :: logic_val

character(*), target :: li(:), attrib_name
character(*), optional :: str_val
character(200), pointer :: line
character(40) value, name

!

nl2 = nl2 + 1
line => li(nl2)

! If the number of second column parameters exceeds number of first column parameters.
if (nl < nl2) then
  nl = nl2
  line = ''  ! In case there is garbage
endif

if (present(re_val)) then
  write (value, '(es15.7)') re_val
elseif (present(str_val)) then
  value = str_val
elseif (present(logic_val)) then
  write (value, '(l1)') logic_val
elseif (present(int_val)) then
  write (value, '(i6)') int_val
else
  call err_exit
endif

name = attrib_name
n = 8 + n_attrib_string_max_len() + 31
write (line(n:), '(a27, a, 2x, a)') name, '=', value

end subroutine encode_2nd_column_parameter 

!--------------------------------------------------------------------------
! contains
!+
! Function is_2nd_column_attribute (ele, attrib_name, ix2_attrib) result (is_2nd_col_attrib)
!
! Function to:
!     1) Return True if attribute is to be displayed in the 2nd column.
!     2) If attribute is a 1st column attribute with a corresponding 2nd column attribute: Set ix2_attrib accordingly.
! 
! Input:
!   ele           -- ele_struct: Element
!   attrib_name   -- character(*): Name of attribute     
!
! Output:
!   is_2nd_col_attrib -- logical: True if a second column attribute. False otherwise.
!   ix2_attrib        -- integer: If  > 0 --> Index of corresponding second column attribute.
!-

function is_2nd_column_attribute (ele, attrib_name, ix2_attrib) result (is_2nd_col_attrib)

type (ele_struct) ele
integer ix, ia, ix_attrib, ix2_attrib
character(*) attrib_name
character(40) a_name, a2_name
logical is_2nd_col_attrib

character(42), parameter :: att_name(102) = [character(42):: 'X_PITCH', 'Y_PITCH', 'X_OFFSET', &
                'Y_OFFSET', 'Z_OFFSET', 'REF_TILT', 'TILT', 'ROLL', 'X1_LIMIT', 'Y1_LIMIT', &
                'FB1', 'FQ1', 'LORD_PAD1', 'HKICK', 'VKICK', 'KICK', 'FRINGE_TYPE', 'DS_STEP', 'R0_MAG', &
                'KS', 'K1', 'K2', 'G', 'DG', 'G_TOT', 'H1', 'E1', 'FINT', 'HGAP', &
                'L_CHORD', 'PTC_FIELD_GEOMETRY', 'AUTOSCALE_AMPLITUDE', 'FIELD_AUTOSCALE', 'COUPLER_AT', &
                'VOLTAGE', 'VOLTAGE_TOT', 'PHI0', 'N_CELL', 'X_GAIN_ERR', 'X_GAIN_CALIB', 'X_OFFSET_CALIB', &
                'BETA_A', 'ALPHA_A', 'CRAB_X1', 'CRAB_X2', 'CRAB_TILT', 'PX_APERTURE_WIDTH2', 'PY_APERTURE_WIDTH2', &
                'PZ_APERTURE_WIDTH2', 'Z_APERTURE_WIDTH2', 'CMAT_11', 'CMAT_21', 'X_DISPERSION_ERR', &
                'X_DISPERSION_CALIB', 'K1X', 'RF_FREQUENCY', 'UPSTREAM_ELE_DIR', 'SIG_X', &
                'BETA_A0', 'BETA_B0', 'ALPHA_A0', 'ALPHA_B0', 'ETA_X0', 'ETAP_X0', 'X1_EDGE', 'Y1_EDGE', &
                'ETA_Y0', 'ETAP_Y0', 'KICK0', 'X0', 'PX0', 'Y0', 'PY0', 'Z0', 'PZ0', &
                'C11_MAT0', 'C12_MAT0', 'C21_MAT0', 'C22_MAT0', 'HARMON', 'FINAL_CHARGE', &
                'MODE_FLIP0', 'BETA_A_STRONG', 'BETA_B_STRONG', 'REF_TIME_START', 'THICKNESS', &
                'PX_KICK', 'PY_KICK', 'PZ_KICK', 'E_TOT_OFFSET', 'FLEXIBLE', 'CRUNCH', 'NOISE', &
                'F_FACTOR', 'EXACT_MULTIPOLES', 'CROSSING_TIME', 'SPIN_TRACKING_MODEL', &
                'SPIN_DN_DPZ_X', 'INHERIT_FROM_FORK', 'N_PERIOD', 'G_MAX', 'PC_STRONG']

character(42), parameter :: att2_name(102) = [character(42):: 'X_PITCH_TOT', 'Y_PITCH_TOT', 'X_OFFSET_TOT', &
                'Y_OFFSET_TOT', 'Z_OFFSET_TOT', 'REF_TILT_TOT', 'TILT_TOT', 'ROLL_TOT', 'X2_LIMIT', 'Y2_LIMIT', &
                'FB2', 'FQ2', 'LORD_PAD2', 'BL_HKICK', 'BL_VKICK', 'BL_KICK', 'FRINGE_AT', 'NUM_STEPS', 'R0_ELEC', &
                'BS_FIELD', 'B1_GRADIENT', 'B2_GRADIENT', 'B_FIELD', 'DB_FIELD', 'B_FIELD_TOT', 'H2', 'E2', 'FINTX', 'HGAPX', &
                'L_RECTANGLE', 'PTC_FRINGE_GEOMETRY', 'AUTOSCALE_PHASE', 'PHI0_AUTOSCALE', 'COUPLER_STRENGTH', &
                'GRADIENT', 'GRADIENT_TOT', 'PHI0_MULTIPASS', 'CAVITY_TYPE', 'Y_GAIN_ERR', 'Y_GAIN_CALIB', 'Y_OFFSET_CALIB', &
                'BETA_B', 'ALPHA_B', 'CRAB_X3', 'CRAB_X4', 'CRAB_X5', 'PX_APERTURE_CENTER', 'PY_APERTURE_CENTER', &
                'PZ_APERTURE_CENTER', 'Z_APERTURE_CENTER', 'CMAT_12', 'CMAT_22', 'Y_DISPERSION_ERR', &
                'Y_DISPERSION_CALIB', 'K1Y', 'RF_WAVELENGTH', 'DOWNSTREAM_ELE_DIR', 'SIG_Y', &
                'BETA_A1', 'BETA_B1', 'ALPHA_A1', 'ALPHA_B1', 'ETA_X1', 'ETAP_X1', 'X2_EDGE', 'Y2_EDGE', &
                'ETA_Y1', 'ETAP_Y1', 'MATRIX', 'X1', 'PX1', 'Y1', 'PY1', 'Z1', 'PZ1', &
                'C11_MAT1', 'C12_MAT1', 'C21_MAT1', 'C22_MAT1', 'HARMON_MASTER', 'SCATTER', &
                'MODE_FLIP1', 'ALPHA_A_STRONG', 'ALPHA_B_STRONG', 'DELTA_REF_TIME', 'DTHICKNESS_DX', &
                'X_KICK', 'Y_KICK', 'Z_KICK', 'E_TOT_START', 'REF_COORDS', 'CRUNCH_CALIB', 'N_SAMPLE', &
                'SCATTER_METHOD', 'FIDUCIAL_PT', 'S_BETA_MIN', 'RECALC', &
                'SPIN_DN_DPZ_Y', 'MODE_FLIP', 'L_PERIOD', 'B_MAX', 'E_TOT_STRONG']

! Exceptional cases

ix2_attrib = -1
is_2nd_col_attrib = .false.

select case (ele%key)
case (patch$)
  select case (attrib_name)
  case ('X_PITCH')
    ix2_attrib = y_pitch$
    return
  case ('Y_PITCH')
    is_2nd_col_attrib = .true.
    return
  case ('X_OFFSET')
    ix2_attrib = y_pitch$
    return
  case ('Y_OFFSET')
    is_2nd_col_attrib = .true.
    return
  end select
end select


select case (attrib_name)
case ('L')
  is_2nd_col_attrib = .false.

  if (ele%key == patch$) then
    ix2_attrib = user_sets_length$
  elseif (ele%key == sbend$) then
    ix2_attrib = l_sagitta$
  elseif (has_attribute(ele, 'L_ACTIVE')) then
    ix2_attrib = l_active$
  elseif (has_attribute(ele, 'L_SOFT_EDGE')) then
    ix2_attrib = l_soft_edge$
  endif

  return

case ('L_SOFT_EDGE', 'L_ACTIVE', 'L_SAGITTA', 'USER_SETS_LENGTH')
  is_2nd_col_attrib = .true.
  return
end select

! Is a 2nd column attribute if corresponding first column attribute exists
! Note: There are rare cases where the corresponding 1st column parameter does not exist for the particular element type.

call match_word (attrib_name, att2_name, ix, .true., .false.)
if (ix > 0) then
  if (.not. has_attribute(ele, att_name(ix))) return
  ia = attribute_index(ele, att_name(ix))
  is_2nd_col_attrib = (ia > 0)
  return
endif

! If the attribute has a corresponding 2nd column attribute, set ix2_attrib accordingly.

call match_word (attrib_name, att_name, ix, .true., .false.)
if (ix > 0) then
  if (.not. has_attribute(ele, att2_name(ix))) return
  ix2_attrib = attribute_index(ele, att2_name(ix))
endif

! Temp until bend fiducial_pt code finished

!! if (ix2_attrib == fiducial_pt$ .or. ix2_attrib == l_rectangle$) ix2_attrib = -1

end function is_2nd_column_attribute

!--------------------------------------------------------------------------
! contains

subroutine write_this_attribute (attrib, ia, n_name_width, line)

type (ele_attribute_struct) attrib
integer ia, n_name_width
character(*) line
character(40) name
character(3) str_ix

!

if (attrib%kind == does_not_exist$) return

if (attrib%ix_attrib > 0) then
  write (str_ix, '(i3)') attrib%ix_attrib
else
  str_ix = ''
endif

select case (attrib%kind)
case (is_logical$)
  if (ele%value(ia) /= 0) ele%value(ia) = 1
  write (line, '(a, 2x, 2a, l1, a, i0, a)')  str_ix, attrib%name(1:n_name_width), '=  ', is_true(attrib%value), ' (', nint(attrib%value), ')'
case (is_integer$)
  write (line, '(a, 2x, 2a, i0)')  str_ix, attrib%name(1:n_name_width), '= ', nint(attrib%value)
case (is_real$)
  write (line, '(a, 2x, 2a, es15.7, 1x, a8)')  str_ix, attrib%name(1:n_name_width), '=', attrib%value, attrib%units
case (is_switch$)
  if (attrib%name == 'SPECIES_STRONG' .and. attrib%value == real_garbage$) then
    name = 'Not_Set'
  else
    name = switch_attrib_value_name (attrib%name, attrib%value, ele, is_default)
  endif
  write (line, '(a, 2x, 4a, i0, a)')  str_ix, attrib%name(1:n_name_width), '=  ', trim(name), ' (', nint(attrib%value), ')'
end select

end subroutine write_this_attribute

!--------------------------------------------------------------------------
! contains

function cmplx_re_str(cmp) result (str_out)

complex(rp) cmp
character(40) str_out

!

if (imag(cmp) == 0) then
  str_out = real_str(real(cmp), 8)
else
  str_out = '(' // real_str(real(cmp), 8) // ', ' // real_str(imag(cmp), 8) // ')'
endif

end function cmplx_re_str

!--------------------------------------------------------------------------
! contains

subroutine print_this_stack(nl, li, stack, print_it, str_index)

type (expression_atom_struct) :: stack(:) 
type (str_index_struct) str_index

integer nl
logical print_it
character(200), allocatable, target :: li(:)

!

if (nl + size(stack) + 100 > size(li)) call re_allocate(li, nl + size(stack) + 100)

do im = 1, size(stack)
  if (stack(im)%type == end_stack$) return
  if (stack(im)%type /= variable$) cycle
  if (stack(im)%name == '') cycle
  if (any(stack(im)%name == physical_const_list%name)) cycle
  call find_index(stack(im)%name, str_index, ix, add_to_list = .true., has_been_added = has_been_added)
  if (.not. (has_been_added)) cycle  ! Avoid duuplicates
  if (print_it) then
    nl=nl+1; li(nl) = 'Named Constants:'
    print_it = .false.
  endif
  nl=nl+1; write (li(nl), '(8x, 2a, es15.7)') trim(stack(im)%name), ' = ', stack(im)%value
enddo

end subroutine print_this_stack

!--------------------------------------------------------------------------
! contains

function this_real(value, fmt, garbage_str) result (out_str)

real(rp) value
character(*) fmt, garbage_str
character(30) out_str

!

if (value == real_garbage$) then
  out_str = garbage_str
  return
endif

write (out_str, '(' // fmt // ')') value

end function this_real

end subroutine type_ele
