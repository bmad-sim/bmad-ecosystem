!+
! Subroutine type_ele (ele, type_zero_attrib, type_mat6, type_taylor, twiss_out, 
!        type_control, type_wake, type_floor_coords, type_field, type_wall, lines, n_lines)
!
! Subroutine to print or put in a string array information on a lattice element.
! If the lines(:) argument is not present, the element information is printed to the terminal.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele               -- Ele_struct: Element
!   type_zero_attrib  -- Logical, optional: If False then surpress printing of
!                           real attributes whose value is 0 or switch attributes that have
!                           their default value. Default is False.
!   type_mat6         -- Integer, optional:
!                            = 0   => Do not type ele%mat6
!                            = 4   => Type 4X4 xy submatrix
!                            = 6   => Type full 6x6 matrix (Default)
!   type_taylor       -- Logical, optional: Print out taylor map terms?
!                          If ele%taylor is not allocated then this is ignored.
!                          Default is False.
!   twiss_out         -- Integer, optional: Print the Twiss parameters at the element end?
!                          = 0         => Do not print the Twiss parameters
!                          = radians$  => Print Twiss, phi in radians (Default).
!                          = degrees$  => Print Twiss, phi in degrees.
!                          = cycles$   => Print Twiss, phi in radians/2pi.
!   type_control       -- Logical, optional: Print control status? Default is True.
!                           If ele%branch%lat is not associated cannot print status info.
!   type_wake         -- Logical, optional: If True then print the long-range and 
!                          short-range wakes information. If False then just print
!                          how many terms the wake has. Default is True.
!                          If ele%wake is not allocated then this is ignored.
!   type_floor_coords -- Logical, optional: If True then print the global ("floor")
!                          coordinates at the exit end of the element.
!                          Default is False.
!   type_field        -- Logical, optional: If True then print:
!                          Wiggler terms for a a map_type wiggler or
!                          RF field info for a lcavity or rfcavity, etc.
!                          Default is False.
!   type_wall         -- Logical, optional: If True then print wall info. Default is False.
!
! Output       
!   lines(:)     -- Character(200), allocatable, optional :: Character array to hold the output. 
!                     If not present, the information is printed to the terminal.
!   n_lines      -- Integer, optional: Number of lines in lines(:) that hold valid output.
!                     n_lines must be present if lines(:) is. 
!-

subroutine type_ele (ele, type_zero_attrib, type_mat6, type_taylor, twiss_out, &
             type_control, type_wake, type_floor_coords, type_field, type_wall, lines, n_lines)

use multipole_mod, except_dummy => type_ele
use geometry_mod, only: coords_relative_to_floor, ele_geometry
use expression_mod

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord, slave, ele0
type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (floor_position_struct) :: floor, f0
type (wake_lr_mode_struct), pointer :: lr
type (wake_sr_mode_struct), pointer :: mode
type (cartesian_map_struct), pointer :: ct_map
type (cartesian_map_term1_struct), pointer :: ct_term
type (cylindrical_map_struct), pointer :: cl_map
type (cylindrical_map_term1_struct), pointer :: cl_term
type (grid_field_struct), pointer :: g_field
type (grid_field_pt1_struct), pointer :: g_pt
type (taylor_field_struct), pointer :: t_field
type (taylor_field_plane1_struct), pointer :: t_term
type (wall3d_struct), pointer :: wall3d
type (wall3d_section_struct), pointer :: section
type (wall3d_vertex_struct), pointer :: v
type (photon_element_struct), pointer :: p
type (photon_surface_struct), pointer :: s
type (ele_attribute_struct) attrib
type (lat_param_struct) param
type (control_struct), pointer :: ctl
type (all_pointer_struct) a_ptr
type (ac_kicker_struct), pointer :: ac

integer, optional, intent(in) :: type_mat6, twiss_out
integer, optional, intent(out) :: n_lines
integer i, i1, j, n, ix, iw, ix_tot, iv, ic, nl2, l_status, a_type, default_val
integer nl, nt, n_term, n_att, attrib_type, n_char, iy, particle, ix_pole_max

real(rp) coef, val
real(rp) a(0:n_pole_maxx), b(0:n_pole_maxx)
real(rp) a2(0:n_pole_maxx), b2(0:n_pole_maxx)
real(rp) knl(0:n_pole_maxx), tn(0:n_pole_maxx)
real(rp), pointer :: r_ptr

character(*), optional, allocatable :: lines(:)
character(200), pointer :: li(:)
character(200), allocatable :: li2(:)
character(200) coef_str
character(60) str1, str2
character(40) a_name, name, fmt_r, fmt_a, fmt_i, fmt_l, fmt
character(12) val_str, units
character(8) angle, index_str

character(*), parameter :: r_name = 'type_ele'

logical, optional, intent(in) :: type_taylor, type_wake
logical, optional, intent(in) :: type_control, type_zero_attrib
logical, optional :: type_floor_coords, type_field, type_wall
logical type_zero, err_flag, print_it, is_default, has_it

! init

allocate (li(300))

type_zero = logic_option(.false., type_zero_attrib)

if (associated(ele%branch)) call lat_sanity_check(ele%branch%lat, err_flag)

if (associated(ele%branch)) then
  branch => ele%branch
  lat => branch%lat
  particle = branch%param%particle
else
  nullify(lat)
  particle = electron$
endif

! Encode element name and type

nl = 0  

if (ele%ix_branch /= 0) then
  if (associated(lat)) then
    nl=nl+1; write (li(nl), *) 'Branch #      ', ele%ix_branch, ': ', trim(branch%name)
  else
    nl=nl+1; write (li(nl), *) 'Branch #      ', ele%ix_branch
  endif
endif
nl=nl+1; write (li(nl), *)     'Element #     ', ele%ix_ele
nl=nl+1; write (li(nl), *)     'Element Name: ', trim(ele%name)

if (ele%type /= blank_name$) then
  nl=nl+1; write (li(nl), *) 'Element Type: "', trim(ele%type), '"'
endif

if (ele%alias /= blank_name$) then
  nl=nl+1; write (li(nl), *) 'Element Alias: "', trim(ele%alias), '"'
endif

if (associated(ele%descrip)) then
  nl=nl+1; write (li(nl), *) 'Descrip: "', trim(ele%descrip), '"'
endif

! Encode element key and attributes

if (ele%key <= 0) then
  nl=nl+1; write (li(nl), *) 'Key: BAD VALUE!', ele%key
else
  nl=nl+1; write (li(nl), *) 'Key: ', key_name(ele%key)
endif

if (ele%sub_key /= 0) then
  nl=nl+1; write (li(nl), *) 'Sub Key: ', sub_key_name(ele%sub_key)
endif

nl=nl+1; write (li(nl), '(1x, 3(a, f12.6))')  'S_start, S:',  ele%s_start, ',', ele%s
nl=nl+1; write (li(nl), '(1x, a, es14.6)') 'Ref_time:', ele%ref_time

nl=nl+1; li(nl) = ''
if (type_zero) then
  nl=nl+1; write (li(nl), *) 'Attribute values:'
else
  nl=nl+1; write (li(nl), *) 'Attribute values [Only non-zero/non-default values shown]:'
endif

n_att = n_attrib_string_max_len() + 2
write (fmt_a, '(a, i0, a)') '(9x, a, t', n_att+10, ', a, 2x, 3a)'
write (fmt_i, '(a, i0, a)') '(9x, a, t', n_att+10, ', a, i6)'
write (fmt_l, '(a, i0, a)') '(9x, a, t', n_att+10, ', a, 2x, l1)'
write (fmt_r, '(a, i0, a)') '(9x, a, t', n_att+10, ', a, 2x, es15.7)'


do i = 1, num_ele_attrib$
  attrib = attribute_info(ele, i)
  a_name = attrib%name
  if (a_name == null_name$) cycle
  if (attrib%type == private$) cycle
  ix_tot = corresponding_tot_attribute_index (ele, i)

  if (i >= custom_attribute1$ .and. i <= custom_attribute5$) then
    nl=nl+1; write (li(nl), '(i5, 3x, 2a, es15.7, 3x, a)') &
                      i, a_name(1:n_att), '=', ele%value(i), '! Custom attribute'

  elseif (ix_tot > 0) then
    if (ele%value(i) == 0 .and. ele%value(ix_tot) == 0 .and. .not. type_zero) cycle
    nl=nl+1; write (li(nl), '(i5, 3x, 2a, es15.7, 1x, a8, i3, 3x, a16, a, es15.7, 1x, a8)') &
                      i, a_name(1:n_att), '=', ele%value(i), attrib%units, &
                      ix_tot, attribute_name(ele, ix_tot), '=', ele%value(ix_tot), attrib%units

  elseif (a_name == 'RF_FREQUENCY' .and. ele%value(i) /= 0) then
    nl=nl+1; write (li(nl), '(i5, 3x, 2a, es15.7, 1x, a8, 6x, a, f13.9, 1x, a)') &
                      i, a_name(1:n_att), '=', ele%value(i), attrib%units, &
                      'RF_WAVELENGTH   =', c_light * ele%value(p0c$) / (ele%value(i) * ele%value(e_tot$)), 'm'

  elseif (a_name == 'P0C_START') then
    nl=nl+1; write (li(nl), '(i5, 3x, 2a, es15.7, 10x, a, f13.9)') &
                      i, a_name(1:n_att), '=', ele%value(i), &
                      'BETA_START      =', ele%value(p0c_start$) / ele%value(e_tot_start$)

  elseif (a_name == 'E_TOT_START') then
    nl=nl+1; write (li(nl), '(i5, 3x, 2a, es15.7, 1x, a8, 6x, a, es15.7, 1x, a8)') &
                      i, a_name(1:n_att), '=', ele%value(i), attrib%units, &
                      'DELTA_E         =', ele%value(e_tot$) - ele%value(e_tot_start$), attrib%units

  elseif (a_name == 'P0C') then
    if (particle == photon$) then
      nl=nl+1; write (li(nl), '(i5, 3x, 2a, es15.7, 1x, a8, 6x, a, es14.6, 1x, a)') &
                      i, a_name(1:n_att), '=', ele%value(i), attrib%units, &
                      'REF_WAVELENGTH  =', c_light * h_planck / ele%value(p0c$), 'm'
    else
      nl=nl+1; write (li(nl), '(i5, 3x, 2a, es15.7, 1x, a8, 6x, a, f13.9)') &
                      i, a_name(1:n_att), '=', ele%value(i), attrib%units, &
                      'BETA            =', ele%value(p0c$) / ele%value(e_tot$)
    endif

  elseif (a_name == 'E_TOT') then
    nl=nl+1; write (li(nl), '(i5, 3x, 2a, es15.7, 1x, a8, 6x, a, es15.7)') &
                      i, a_name(1:n_att), '=', ele%value(i), attrib%units, &
                      'GAMMA           =', ele%value(e_tot$) / mass_of(particle)

  elseif (index(a_name, 'ANGLE') /= 0 .and. a_name /= 'CRITICAL_ANGLE_FACTOR') then
    units = ' deg'
    if (a_name == 'DBRAGG_ANGLE_DE') units = ' deg/eV'
    if (.not. type_zero .and. ele%value(i) == 0) cycle
    nl=nl+1; write (li(nl), '(i5, 3x, 2a, es15.7, 1x, a8, f10.4, a)') &
                 i, a_name(1:n_att), '=', ele%value(i), attrib%units, ele%value(i) * 180 / pi, trim(units)
  else
    attrib_type = attribute_type(a_name)
    if (is_a_tot_attribute(ele, i)) cycle
    select case (attrib_type)
    case (is_logical$)
      if (ele%value(i) /= 0) ele%value(i) = 1
      nl=nl+1; write (li(nl), '(i5, 3x, 2a, l1, a, i0, a)')  i, a_name(1:n_att), '=  ', &
                                  is_true(ele%value(i)), ' (', nint(ele%value(i)), ')'
    case (is_integer$)
      if (ele%value(i) == 0 .and. .not. type_zero) cycle
      nl=nl+1; write (li(nl), '(i5, 3x, 2a, i0)')  i, a_name(1:n_att), '= ', nint(ele%value(i))
    case (is_real$)
      if (ele%value(i) == 0 .and. .not. type_zero) cycle
      nl=nl+1; write (li(nl), '(i5, 3x, 2a, es15.7, 1x, a8)')  i, a_name(1:n_att), '=', ele%value(i), attrib%units
    case (is_switch$)
      name = switch_attrib_value_name (a_name, ele%value(i), ele, is_default)
      if (.not. is_default .or. type_zero) then
        nl=nl+1; write (li(nl), '(i5, 3x, 4a, i0, a)')  i, a_name(1:n_att), '=  ', &
                                                      trim(name), ' (', nint(ele%value(i)), ')'
      endif
    end select
  endif
enddo

! Multipoles

if (associated(ele%a_pole) .or. associated(ele%a_pole_elec)) then
  nl=nl+1; write (li(nl), '(a, l1)')   'Multipoles_on:    ', ele%multipoles_on 
endif

if (associated(ele%a_pole)) then
  if (attribute_index(ele, 'SCALE_MULTIPOLES') == scale_multipoles$) then
    nl=nl+1; write (li(nl), '(a, l1)') 'Scale_Multipoles: ', ele%scale_multipoles
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

  do i = 0, n_pole_maxx
    if (ele%a_pole(i) == 0 .and. ele%b_pole(i) == 0) cycle

    if (ele%key == multipole$) then
      nl=nl+1; write (li(nl), '(2x, 3(3x, a, i0, a, es11.3))') &
             'K', i, 'L =', ele%a_pole(i), 'K', i, 'L(w/Tilt) =', knl(i), 'A', i, '(equiv) =', a(i)
      nl=nl+1; write (li(nl), '(2x, 3(3x, a, i0, a, es11.3))') &
             'T', i, '  =', ele%b_pole(i), 'T', i, '(w/Tilt)  =', tn(i),  'B', i, '(equiv) =', b(i)
    elseif (ele%key == ab_multipole$) then
      nl=nl+1; write (li(nl), '(2x, 3(3x, a, i0, a, es11.3))') &
             'A', i, ' =', ele%a_pole(i), 'A', i, '(w/Tilt) =', a2(i), 'K', i, 'L(equiv) =', knl(i)
      nl=nl+1; write (li(nl), '(2x, 3(3x, a, i0, a, es11.3))') &
             'B', i, ' =', ele%b_pole(i), 'B', i, '(w/Tilt) =', b2(i), 'T', i, '(equiv)  =', tn(i)
    else
      nl=nl+1; write (li(nl), '(2x, 4(3x, a, i0, a, es11.3))') &
             'A', i, ' =', ele%a_pole(i), &
             'A', i, '(Scaled) =', a(i), 'A', i, '(w/Tilt) =', a2(i), 'K', i, 'L(equiv) =', knl(i)
      nl=nl+1; write (li(nl), '(2x, 4(3x, a, i0, a, es11.3))') &
             'B', i, ' =', ele%b_pole(i), &
             'B', i, '(Scaled) =', b(i), 'B', i, '(w/Tilt) =', b2(i), 'T', i, '(equiv)  =', tn(i)
    endif

  enddo

endif

! Electric Multipoles

if (associated(ele%a_pole_elec)) then

  do i = 0, n_pole_maxx
    call multipole_ele_to_ab (ele, .false., ix_pole_max, a, b, electric$)
    if (ele%a_pole_elec(i) == 0 .and. ele%b_pole_elec(i) == 0) cycle
    nl=nl+1; write (li(nl), '(2x, 4(3x, a, i0, a, es11.3))') 'A', i, '_elec =', ele%a_pole_elec(i), 'A', i, '_elec(Scaled) =', a(i)
    nl=nl+1; write (li(nl), '(2x, 4(3x, a, i0, a, es11.3))') 'B', i, '_elec =', ele%b_pole_elec(i), 'B', i, '_elec(Scaled) =', b(i)
  enddo

endif

! Encode on/off status etc.

if (.not. ele%is_on) then
  nl=nl+1; write (li(nl), *) '*** Note: Element is turned OFF ***'
endif

! Encode methods, etc.

nl=nl+1; write (li(nl), *) ' '

if (ele%key == beambeam$ .and. associated(ele%branch)) then
  nl=nl+1; write (li(nl), fmt_r) 'PARAMETER[N_PART]', '=', ele%branch%param%n_part
endif

if (attribute_name(ele, crystal_type$) == 'CRYSTAL_TYPE') then
  nl=nl+1; write (li(nl), fmt_a) 'CRYSTAL_TYPE', '=', ele%component_name
endif

if (attribute_name(ele, material_type$) == 'MATERIAL_TYPE') then
  nl=nl+1; write (li(nl), fmt_a) 'MATERIAL_TYPE', '=', ele%component_name
endif

if (attribute_name(ele, origin_ele$) == 'ORIGIN_ELE') then
  nl=nl+1; write (li(nl), fmt_a) 'ORIGIN_ELE', '=', '"', trim(ele%component_name), '"'
endif

if (attribute_name(ele, physical_source$) == 'PHYSICAL_SOURCE') then
  nl=nl+1; write (li(nl), fmt_a) 'PHYSICAL_SOURCE', '=', '"', trim(ele%component_name), '"'
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

if (attribute_name(ele, field_calc$) == 'FIELD_CALC') then
  nl=nl+1; write (li(nl), fmt_a) 'FIELD_CALC', '=', field_calc_name(ele%field_calc)
endif

! Write aparture stuff if appropriate

if (attribute_name(ele, aperture_at$) == 'APERTURE_AT' .and. ele%aperture_at /= 0) then
  nl=nl+1; write (li(nl), fmt_a) 'APERTURE_AT', '=', aperture_at_name(ele%aperture_at)
  default_val = rectangular$
  if (ele%key == ecollimator$) default_val = elliptical$
  if (ele%aperture_type /= default_val .or. type_zero) then
    nl=nl+1; write (li(nl), fmt_a) 'APERTURE_TYPE', '=', aperture_type_name(ele%aperture_type)
  endif
  nl=nl+1; write (li(nl), fmt_l) 'OFFSET_MOVES_APERTURE', '=', ele%offset_moves_aperture
endif

if (attribute_index(ele, 'SYMPLECTIFY') /= 0) then
  nl=nl+1; write (li(nl), fmt_l) 'SYMPLECTIFY', '=', ele%symplectify
endif
  
if (attribute_index(ele, 'FIELD_MASTER') /= 0) then
  nl=nl+1; write (li(nl), fmt_l) 'FIELD_MASTER', '=', ele%field_master
endif

if (attribute_index(ele, 'HARMON_MASTER') /= 0) then
  nl=nl+1; write (li(nl), fmt_l) 'HARMON_MASTER', '=', ele%field_master
endif

if (attribute_index(ele, 'CSR_CALC_ON') /= 0) then
  nl=nl+1; write (li(nl), fmt_l) 'CSR_CALC_ON', '=', ele%csr_calc_on
endif

if (ele%key /= overlay$ .and. ele%key /= group$) then
  if (ele%orientation /= 1 .or. type_zero) then
    nl=nl+1; write (li(nl), fmt_i) 'LONGITUDINAL ORIENTATION', '=', ele%orientation
  endif
endif

! Cartesian map

if (associated(ele%cartesian_map)) then
  if (logic_option(.false., type_field)) then
    nl=nl+1; li(nl) = ''
    if (ele%field_calc == bmad_standard$) then
      nl=nl+1; li(nl) = 'Cartesian_map: [NOT USED SINCE FIELD_CALC = BMAD_STANDARD]'
    else
      nl=nl+1; li(nl) = 'Cartesian_map:'
    endif
    do i = 1, size(ele%cartesian_map)
      ct_map => ele%cartesian_map(i)
      if (ct_map%master_parameter == 0) then
        name = '<None>'
      else
        name = attribute_name(ele, ct_map%master_parameter)
      endif

      nl=nl+1; write (li(nl), '(a, i0)')      '  Mode #:', i
      nl=nl+1; write (li(nl), '(2a)')         '    From file:        ', trim(ct_map%ptr%file)
      nl=nl+1; write (li(nl), '(2a)')         '    field_type        ', trim(em_field_type_name(ct_map%field_type))
      nl=nl+1; write (li(nl), '(2a)')         '    master_parameter: ', trim(name)
      nl=nl+1; write (li(nl), '(2a)')         '    ele_anchor_pt:    ', anchor_pt_name(ct_map%ele_anchor_pt)
      nl=nl+1; write (li(nl), '(a, es16.8)')  '    field_scale:      ', ct_map%field_scale
      nl=nl+1; write (li(nl), '(a, 3es16.8)') '    r0:               ', ct_map%r0
      nl=nl+1; write (li(nl), '(a, i0)')      '    n_link:           ', ct_map%ptr%n_link
      nl=nl+1; write (li(nl), '(5x, a, 6x, a, 3(9x, a), 2(12x, a), 9x, a, 3x, a)') 'Term#', &
                                    'Coef', 'K_x', 'K_y', 'K_z', 'x0', 'y0', 'phi_z', 'Type'
      do j = 1, min(10, size(ct_map%ptr%term))
        if (nl+1 > size(li)) call re_associate(li, 2 * nl, .false.)
        ct_term => ct_map%ptr%term(j)
        nl=nl+1; write (li(nl+i), '(i8, 4f12.6, 3f14.6, 3x, a)') i, ct_term%coef, ct_term%kx, ct_term%ky, ct_term%kz, &
                               ct_term%x0, ct_term%y0, ct_term%phi_z, cartesian_map_term_name(ct_term%type)
      enddo
      if (size(ct_map%ptr%term) > 10) then
        nl=nl+1; write (li(nl), '(a, i0, a)') '     .... etc ... (#Terms = ', size(ct_map%ptr%term), ')' 
      endif
    enddo
  else
    nl=nl+1; write (li(nl), '(a, i5)') 'Number of Cartesian_map modes:', size(ele%cartesian_map)
  endif
endif

! Cylindrical_map

if (associated(ele%cylindrical_map)) then
  if (logic_option(.false., type_field)) then
    nl=nl+1; li(nl) = ''
    if (ele%field_calc == bmad_standard$) then
      nl=nl+1; li(nl) = 'Cylindrical_map: [NOT USED SINCE FIELD_CALC = BMAD_STANDARD]'
    else
      nl=nl+1; li(nl) = 'Cylindrical_map:'
    endif
    do i = 1, size(ele%cylindrical_map)
      cl_map => ele%cylindrical_map(i)
      if (cl_map%master_parameter == 0) then
        name = '<None>'
      else
        name = attribute_name(ele, cl_map%master_parameter)
      endif

      nl=nl+1; write (li(nl), '(a, i0)')      '  Mode #:', i
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
      do j = 1, min(10, size(cl_map%ptr%term))
        if (nl+1 > size(li)) call re_associate(li, 2 * nl, .false.)
        cl_term => cl_map%ptr%term(j)
        nl=nl+1; write (li(nl), '(i5, 3x, 2(a, 2es12.4), a)') j, '(', cl_term%e_coef, ')  (', cl_term%b_coef, ')'
      enddo
      if (size(cl_map%ptr%term) > 10) then
        nl=nl+1; write (li(nl), '(a, i0, a)') '     .... etc ... (#Terms = ', size(cl_map%ptr%term), ')' 
      endif
    enddo
  else
    nl=nl+1; write (li(nl), '(a, i5)') 'Number of Cylindrical_map modes:', size(ele%cylindrical_map)
  endif
endif

! Grid_field

if (associated(ele%grid_field)) then
  if (logic_option(.false., type_field)) then
    nl=nl+1; li(nl) = ''
    if (ele%field_calc == bmad_standard$) then
      nl=nl+1; li(nl) = 'Grid_field: [NOT USED SINCE FIELD_CALC = BMAD_STANDARD]'
    else
      nl=nl+1; li(nl) = 'Grid_field:'
    endif
    do i = 1, size(ele%grid_field)
      g_field => ele%grid_field(i)
      if (g_field%master_parameter == 0) then
        name = '<None>'
      else
        name = attribute_name(ele, g_field%master_parameter)
      endif

      nl=nl+1; write (li(nl), '(a, i0)')      '  Mode #:', i
      nl=nl+1; write (li(nl), '(2a)')         '    From file:        ', trim(g_field%ptr%file)
      nl=nl+1; write (li(nl), '(2a)')         '    field_type:       ', em_field_type_name(g_field%field_type)
      nl=nl+1; write (li(nl), '(2a)')         '    geometry:         ', grid_field_geometry_name(g_field%geometry)
      nl=nl+1; write (li(nl), '(2a)')         '    master_parameter: ', trim(name)
      nl=nl+1; write (li(nl), '(2a)')         '    ele_anchor_pt:    ', anchor_pt_name(g_field%ele_anchor_pt)
      nl=nl+1; write (li(nl), '(a, i0)')      '    harmonic:         ', g_field%harmonic
      nl=nl+1; write (li(nl), '(a, es16.8)')  '    field_scale:      ', g_field%field_scale
      nl=nl+1; write (li(nl), '(a, es16.8)')  '    phi0_fieldmap:    ', g_field%phi0_fieldmap
      nl=nl+1; write (li(nl), '(a, l1)')      '    curved_ref_frame  ', g_field%curved_ref_frame
      nl=nl+1; write (li(nl), '(a, i0)')      '    n_link:           ', g_field%ptr%n_link
      nl=nl+1; write (li(nl), '(a, 3f14.6)')  '    dr:               ', g_field%dr
      nl=nl+1; write (li(nl), '(a, 3f14.6)')  '    r0:               ', g_field%r0
      nl=nl+1; write (li(nl), '(a, 3i14)')    '    Index_max:        ', ubound(g_field%ptr%pt)
      nl=nl+1; write (li(nl), '(a, 3i14)')    '    Index_min:        ', lbound(g_field%ptr%pt)
      nl=nl+1; write (li(nl), '(a, 3f14.6)')  '    r_max:            ', ubound(g_field%ptr%pt)*g_field%dr + g_field%r0
      nl=nl+1; write (li(nl), '(a, 3f14.6)')  '    r_min:            ', lbound(g_field%ptr%pt)*g_field%dr + g_field%r0
    enddo
  else
    nl=nl+1; write (li(nl), '(a, i5)') 'Number of Grid_field modes:', size(ele%grid_field)
  endif
endif

! Taylor_field

if (associated(ele%taylor_field)) then
  if (logic_option(.false., type_field)) then
    nl=nl+1; li(nl) = ''
    if (ele%field_calc == bmad_standard$) then
      nl=nl+1; li(nl) = 'Taylor_field: [NOT USED SINCE FIELD_CALC = BMAD_STANDARD]'
    else
      nl=nl+1; li(nl) = 'Taylor_field:'
    endif
    do i = 1, size(ele%taylor_field)
      t_field => ele%taylor_field(i)
      if (t_field%master_parameter == 0) then
        name = '<None>'
      else
        name = attribute_name(ele, t_field%master_parameter)
      endif

      nl=nl+1; write (li(nl), '(a, i0)')      '  Mode #:', i
      nl=nl+1; write (li(nl), '(2a)')         '    From file:         ', trim(t_field%ptr%file)
      nl=nl+1; write (li(nl), '(2a)')         '    field_type:        ', em_field_type_name(t_field%field_type)
      nl=nl+1; write (li(nl), '(a, es16.8)')  '    field_scale:       ', t_field%field_scale
      nl=nl+1; write (li(nl), '(a, es16.8)')  '    dz:                ', t_field%dz
      nl=nl+1; write (li(nl), '(a, 3es16.8)') '    r0:                ', t_field%r0
      nl=nl+1; write (li(nl), '(2a)')         '    master_parameter:  ', trim(name)
      nl=nl+1; write (li(nl), '(2a)')         '    ele_anchor_pt:     ', anchor_pt_name(t_field%ele_anchor_pt)
      nl=nl+1; write (li(nl), '(a, l1)')      '    curved_ref_frame   ', t_field%curved_ref_frame
      nl=nl+1; write (li(nl), '(a, l1)')      '    canonical_tracking ', t_field%canonical_tracking
      nl=nl+1; write (li(nl), '(a, i0)')      '    n_link:            ', t_field%ptr%n_link
      nl=nl+1; write (li(nl), '(a, i0)')      '    n_plane:           ', size(t_field%ptr%plane)
    enddo
  else
    nl=nl+1; write (li(nl), '(a, i5)') 'Number of Taylor_field modes:', size(ele%taylor_field)
  endif
endif

! ac_kick

if (associated(ele%ac_kick) .and. logic_option(.false., type_field)) then
  ac => ele%ac_kick
  nl=nl+1; li(nl) = ''

  if (allocated(ac%amp_vs_time)) then
    nl=nl+1; li(nl) = '     Indx      Time     Amplitude'    
    do i = 1, size(ac%amp_vs_time)
      nl=nl+1; write (li(nl), '(i9, 2es14.6)') i, ac%amp_vs_time(i)%time, ac%amp_vs_time(i)%amp
    enddo
  endif

  if (allocated(ac%frequencies)) then
    nl=nl+1; li(nl) = '     Indx          Freq     Amplitude           Phi'    
    do i = 1, size(ac%frequencies)
      nl=nl+1; write (li(nl), '(i9, 3es14.6)') i, &
                                ac%frequencies(i)%f, ac%frequencies(i)%amp, ac%frequencies(i)%phi
    enddo
  endif
endif

! wall3d cross-sections.
! Do not print more than 100 sections.

if (associated(ele%wall3d)) then
  do iw = 1, size(ele%wall3d)
    wall3d => ele%wall3d(iw)
    nl=nl+1; write (li(nl), '(a, i5)') ''
    nl=nl+1; write (li(nl), '(2a)') 'Wall name: ', trim(wall3d%name)
    nl=nl+1; write (li(nl), '(a, i5)') 'Number of Wall Sections:', size(wall3d%section)
    nl=nl+1; write (li(nl), '(a, 2f11.5)') 'Wall region:',  wall3d%section(1)%s, wall3d%section(size(wall3d%section))%s
    if (logic_option(.false., type_wall)) then
      nl=nl+1; write (li(nl), '(2a)') 'Wall%ele_anchor_pt = ', anchor_pt_name(wall3d%ele_anchor_pt)
      select case (ele%key)
      case (capillary$)
      case (diffraction_plate$, mask$)
        nl=nl+1; write (li(nl), '(a, f10.6)') 'Wall%thickness       = ', wall3d%thickness
        nl=nl+1; write (li(nl), '(3a)') 'Wall%clear_material  = "', trim(wall3d%clear_material), '"'
        nl=nl+1; write (li(nl), '(3a)') 'Wall%opaque_material = "', trim(wall3d%opaque_material), '"'
      case default
        nl=nl+1; write (li(nl), '(a, l)') 'Wall%superimpose     = ', wall3d%superimpose
      end select
      n = min(size(wall3d%section), 100)
      do i = 1, n
        call re_associate (li, nl+100, .false.) 
        section => wall3d%section(i)
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
                    'Wall%Section(', i, '):  S =', section%s, trim(str1),  ',   r0 = (', &
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

p => ele%photon
if (associated(p)) then
  s => ele%photon%surface
  nl=nl+1; write (li(nl), *)
  nl=nl+1; write (li(nl), *) 'Surface:'
 
  if (s%has_curvature) then
    do ix = 0, ubound(s%curvature_xy, 1)
    do iy = 0, ubound(s%curvature_xy, 2)
      if (s%curvature_xy(ix,iy) == 0) cycle
      nl=nl+1; write (li(nl), '(2x, 2(2x, 2(a, i0), a, es14.6))') 'CURVATURE_X', ix, '_Y', iy, ' =', s%curvature_xy(ix,iy)
    enddo
    enddo
  else
    nl=nl+1; li(nl) = '    No Curvature'
  endif
  nl=nl+1; write (li(nl), '(4x, 2a)') 'Grid type:    ', surface_grid_type_name(s%grid%type)
  if (s%grid%type /= off$) then
    nl=nl+1; write (li(nl), '(4x, a, 2f10.6)')   'Grid dr:     ', s%grid%dr
    if (allocated(s%grid%pt)) then
      nl=nl+1; write (li(nl), '(4x, a, 2f10.6)') 'Grid r0:     ', s%grid%r0
      nl=nl+1; write (li(nl), '(4x, a, 2i10)')   'Num grid pts:', ubound(s%grid%pt) + 1
      nl=nl+1; write (li(nl), '(4x, a, 2(a, f10.6, a, f10.6, a, 4x))') &
                                                  'Grid bounds:', &
                        '(', -s%grid%r0(1), ',', -s%grid%r0(1) + ubound(s%grid%pt, 1) * s%grid%dr(1), ')', & 
                        '(', -s%grid%r0(2), ',', -s%grid%r0(2) + ubound(s%grid%pt, 2) * s%grid%dr(2), ')' 
    endif
  endif

  if (p%material%f_h /= 0) then
    if (ele%key == multilayer_mirror$) then
      nl = nl+1; write (li(nl), '(2(a,f10.3))') 'F_0 (Material 1):', real(p%material%f0_m1), ' + I ', aimag(p%material%f0_m1)
      nl = nl+1; write (li(nl), '(2(a,f10.3))') 'F_0 (Material 2):', real(p%material%f0_m2), ' + I ', aimag(p%material%f0_m2)
    else
      nl = nl+1; write (li(nl), '(2(a,f10.3))') 'F_0:             ', real(p%material%f_0), ' + I ', aimag(p%material%f_0)
    endif
    nl = nl+1; write (li(nl), '(2(a,f10.3))') 'F_H:             ', real(p%material%f_h), ' + I ', aimag(p%material%f_h)
    nl = nl+1; write (li(nl), '(2(a,f10.3))') 'F_Hbar:          ', real(p%material%f_hbar), ' + I ', aimag(p%material%f_hbar)
    nl = nl+1; write (li(nl), '(2(a,f10.3))') 'Sqrt(F_H*F_Hbar):', real(p%material%f_hkl), ' + I ', aimag(p%material%f_hkl)
  endif

endif

! Encode branch info

if (ele%key == fork$ .or. ele%key == photon_fork$) then
  
  if (li(nl) /= '') then
    nl=nl+1; li(nl) = ' '
  endif

  n = nint(ele%value(ix_to_branch$))
  i = nint(ele%value(ix_to_element$))
  if (associated(lat)) then
    nl=nl+1; write (li(nl), '(5a, 2(i0, a))') 'Branch to: ', trim(lat%branch(n)%name), '>>', &
                                  trim(lat%branch(n)%ele(i)%name), '  [', n, '>>', i, ']'
  else
    nl=nl+1; write (li(nl), '(a, i0, a, i0)') 'Branch to: ', n, '>>', i
  endif

endif

! Encode lord/slave info.
! For super_lords there is no attribute_name associated with a slave.
! For slaves who are overlay_lords then the attribute_name is obtained by
!   looking at the overlay_lord's 1st slave (slave of slave of the input ele).

if (associated(lat) .and. logic_option(.true., type_control)) then

  if (li(nl) /= '') then
    nl=nl+1; li(nl) = ' '
  endif

  ! Print info on element's lords

  if (ele%slave_status <= 0) then
    nl=nl+1; write (li(nl), '(a)') 'Slave_status: BAD! PLEASE SEEK HELP!', ele%slave_status
  else
    nl=nl+1; write (li(nl), '(2a)') 'Slave_status: ', control_name(ele%slave_status)
  endif

  select case (ele%slave_status)
  case (multipass_slave$)
    lord => pointer_to_lord(ele, 1)
    nl=nl+1; write (li(nl), '(3a, i0, a)') 'Associated Multipass_Lord: ', trim(ele%name), '  (Index: ', lord%ix_ele, ')'
    nl=nl+1; li(nl) = 'Other slaves of this Lord:'
    nl=nl+1; li(nl) = '   Index   Name'
    do i = 1, lord%n_slave
      slave => pointer_to_slave(lord, i)
      if (slave%ix_ele == ele%ix_ele) cycle
      nl=nl+1; write (li(nl), '(i8, 3x, a)') slave%ix_ele, trim(slave%name)
    enddo

  case (super_slave$)
    nl=nl+1; write (li(nl), '(3a, i0, a)') 'Associated Super_Lord(s):'
    nl=nl+1; li(nl) = '   Index   Name                             Type'
    do i = 1, ele%n_lord
      lord => pointer_to_lord(ele, i)
      if (lord%lord_status /= super_lord$) cycle
      nl=nl+1; write (li(nl), '(i8, 3x, a, t45, a)') lord%ix_ele, trim(lord%name), trim(key_name(lord%key))
    enddo
  end select

  ! Print controller lords

  has_it = .false.
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%lord_status == multipass_lord$ .or. lord%lord_status == super_lord$) cycle
    has_it = .true.
    exit
  enddo

  if (has_it) then
    nl=nl+1; li(nl) = 'Controller Lord(s):'
    nl=nl+1; li(nl) = '   Index   Name                            Attribute           Lord_Type           Expression'

    do i = 1, ele%n_lord
      lord => pointer_to_lord (ele, i, ctl)
      select case (lord%lord_status)
      case (super_lord$, multipass_lord$)
        cycle
      case (girder_lord$)
        coef_str = ''
        a_name = ''
        val_str = ''
      case default
        if (allocated(ctl%stack)) then
          coef_str = expression_stack_to_string (ctl%stack)
        else
          coef_str = ''
        endif
        iv = ctl%ix_attrib
        a_name = attribute_name(ele, iv)
      end select

      nl=nl+1; write (li(nl), '(i8, 3x, a32, a18, 2x, a20, a)') &
            lord%ix_ele, lord%name, a_name, key_name(lord%key), trim(coef_str)
    enddo
    nl=nl+1; li(nl) = ''
  endif

  !

  if (ele%n_lord_field /= 0) then
    nl=nl+1; li(nl) = 'Elements whose fields overlap this one:'
    nl=nl+1; li(nl) = '   Index   Name                                Type             S'
    do i = 1, ele%n_lord_field
      lord => pointer_to_lord(ele, ele%n_lord+i)
      nl=nl+1; write (li(nl), '(a8, t12, a30, a16, f10.3)') &
                    trim(ele_loc_to_string(lord)), lord%name, trim(key_name(lord%key)), lord%s
    enddo
  endif

  ! Print info on elements slaves.

  if (ele%lord_status <= 0) then
    nl=nl+1; write (li(nl), '(a)') 'Lord_status: BAD!', ele%lord_status
  else
    nl=nl+1; write (li(nl), '(2a)') 'Lord_status:  ', control_name(ele%lord_status)
  endif

  if (associated(ele%control_var)) then
    nl=nl+1; li(nl) = 'Control Variables:'
    n_att = maxval(len_trim(ele%control_var%name))
    if (ele%lord_status == group_lord$) then
      do i = 1, size(ele%control_var)
        a_name = ele%control_var(i)%name
        nl=nl+1; write (li(nl), '(i5, 3x, 2a, es15.7, 11x, 3a, es15.7)')  i, &
                      a_name(1:n_att), '  =', ele%control_var(i)%value, &
                      'OLD_', a_name(1:n_att), '  =', ele%control_var(i)%old_value
      enddo
    else  ! overlay_lord
      do i = 1, size(ele%control_var)
        nl=nl+1; write (li(nl), '(i5, 3x, 2a, es15.7)')  i, &
                      ele%control_var(i)%name, '  =', ele%control_var(i)%value
      enddo
    endif
  endif

  !

  if (ele%n_slave /= 0) then
    nl=nl+1; write (li(nl), '(a, i4)') 'Slaves:'

    n_char = 10
    do i = 1, ele%n_slave
      slave => pointer_to_slave (ele, i)
      n_char = max(n_char, len_trim(slave%name))
    enddo

    select case (ele%lord_status)

    case (multipass_lord$, super_lord$, girder_lord$)
      nl=nl+1; li(nl) = '   Index   Name';  li(nl)(n_char+14:) = 'Type'
      do i = 1, ele%n_slave
        slave => pointer_to_slave (ele, i)
        nl=nl+1; write (li(nl), '(a8, t12, a, 2x, a)') &
                    trim(ele_loc_to_string(slave)), slave%name(1:n_char), key_name(slave%key)
      enddo

    case default
      nl=nl+1; li(nl) = '   Index   Ele_Name';  li(nl)(n_char+14:) = 'Attribute         Attrib_Value  Expression_Val    Expression'
      do ix = 1, ele%n_slave
        slave => pointer_to_slave (ele, ix, ctl)
        if (allocated(ctl%stack)) then
          coef_str = expression_stack_to_string (ctl%stack)
          call evaluate_expression_stack(ctl%stack, val, err_flag, str1, ele%control_var)
          write (coef_str, '(es12.4, 4x, a)') val, trim(coef_str)
        else
          coef_str = ' ------ '
        endif

        iv = ctl%ix_attrib
        a_name = attribute_name(slave, iv)
        if (iv > num_ele_attrib$) then  ! Special group construct: accordion_edge$, start_edge$, end_edge$, or s_position$
          write (val_str, '(es12.4)') val
        else
          call pointer_to_attribute (slave, a_name, .false., a_ptr, err_flag)
          val_str = ' ----'
          if (associated(a_ptr%r)) write (val_str, '(es12.4)') a_ptr%r
        endif

        nl=nl+1; write (li(nl), '(a8, t12, a, 2x, a18, a, 4x, a)') trim(ele_loc_to_string(slave)), slave%name(1:n_char), a_name, val_str, trim(coef_str)
      enddo
    end select
  endif

  !

  if (ele%n_slave_field /= 0) then
    nl=nl+1; li(nl) = "This element's field overlaps:"
    nl=nl+1; li(nl) = '   Index   Name                                Type             S'
    do i = 1, ele%n_slave_field
      slave => pointer_to_slave(ele, ele%n_slave+i)
      nl=nl+1; write (li(nl), '(a8, t12, a30, a16, f10.3)') &
                    trim(ele_loc_to_string(slave)), slave%name, trim(key_name(slave%key)), slave%s
    enddo
  endif

endif

! Encode Twiss info

if (integer_option(radians$, twiss_out) /= 0 .and. ele%a%beta /= 0) then
  nl=nl+1; li(nl) = ' '
  nl=nl+1; li(nl) = 'Twiss at end of element:'
  call type_twiss (ele, twiss_out, .false., li(nl+1:), nl2)
  nl = nl + nl2
endif

l_status = ele%lord_status
if (l_status /= overlay_lord$ .and. l_status /= multipass_lord$ .and. &
    l_status /= group_lord$ .and. l_status /= girder_lord$) then

  ! Encode mat6 info

  n = integer_option (6, type_mat6)
  if (n > 6) n = 6
  if (n < 0) n = 0

  if (n /= 0) then
    nl=nl+1; li(nl) = ' '
    nl=nl+1; write (li(nl), '(a, es12.3, a)') 'Transfer Matrix : Kick  [Matrix symplectic error:', &
                  mat_symp_error(ele%mat6), ']'
  endif

  if (any(abs(ele%mat6(1:n,1:n)) >= 1d3)) then
    write (fmt, '(a, i0, a)') '(', n, 'es13.5, a, es13.5)'
  else
    write (fmt, '(a, i0, a)') '(', n, 'f12.7, a, es13.5)'
  endif

  do i = 1, n
    nl=nl+1; write (li(nl), fmt) (ele%mat6(i, j), j = 1, n), '   : ', ele%vec0(i)
  enddo

  ! Encode taylor series

  if (associated(ele%taylor(1)%term)) then
    nl=nl+1; li(nl) = ' '
    nl=nl+1; write (li(nl), '(a, l1)') 'taylor_map_includes_offsets: ', ele%taylor_map_includes_offsets
    if (logic_option(.false., type_taylor)) then
      call type_taylors (ele%taylor, lines = li2, n_lines = nt)
      call re_associate (li, nl+nt+100, .false.)
      li(1+nl:nt+nl) = li2(1:nt)
      deallocate (li2)
      nl = nl + nt
    else
      n_term = 0
      do i = 1, size(ele%taylor)
        n_term = n_term + size(ele%taylor(i)%term)
      enddo
      nl=nl+1; write (li(nl), '(a, i6)') 'Taylor map total number of terms:', n_term
    endif
  endif

  if (associated(ele%spin_taylor(1,1)%term) .or. ele%key == taylor$) then
    if (logic_option(.false., type_taylor)) then
      nl=nl+1; li(nl) = ' '
      call type_spin_taylors (ele%spin_taylor, lines = li2, n_lines = nt)
      call re_associate (li, nl+nt+100, .false.)
      li(1+nl:nt+nl) = li2(1:nt)
      deallocate (li2)
      nl = nl + nt
    else
      n_term = 0
      do i = 1, 3; do j = 1, 3
        n_term = n_term + size(ele%spin_taylor(i,j)%term)
      enddo; enddo
      nl=nl+1; write (li(nl), '(a, i6)') 'Spin_Taylor map total number of terms:', n_term
    endif
  endif

endif

! Print wake info

if (associated(ele%wake)) then

  if (logic_option (.true., type_wake) .and. (size(ele%wake%sr_long%mode) /= 0 .or. &
                                                       size(ele%wake%sr_trans%mode) /= 0)) then
    nl=nl+1; li(nl) = ''
    nl=nl+1; li(nl) = 'Short-Range Wake:'
    nl=nl+1; li(nl) = '  SR_File: ' // trim(ele%wake%sr_file)
    nl=nl+1; write (li(nl), '(2x, a, f10.6)') 'z_sr_max =', ele%wake%z_sr_max
  endif

  if (size(ele%wake%sr_long%mode) /= 0) then
    nl=nl+1; write (li(nl), *)
    if (logic_option (.true., type_wake)) then
      call re_associate (li, nl+size(ele%wake%sr_long%mode)+100, .false.)
      nl=nl+1; li(nl) = '  Short-Range Longitudinal Pseudo Modes:'
      nl=nl+1; li(nl) = &
            '   #        Amp        Damp           K         Phi   Polarization  Transverse_Dependence'
      do i = 1, size(ele%wake%sr_long%mode)
        mode => ele%wake%sr_long%mode(i)
        nl=nl+1; write (li(nl), '(i4, 4es12.4, a15, a16)') i, mode%amp, mode%damp, mode%k, mode%phi, &
                  sr_polarization_name(mode%polarization), sr_transverse_dependence_name(mode%transverse_dependence)
      enddo
    else
      nl=nl+1; li(nl) = '  No short-range longitudinal pseudo modes.'
    endif
  endif

  if (size(ele%wake%sr_trans%mode) /= 0) then
    nl=nl+1; write (li(nl), *)
    if (logic_option (.true., type_wake)) then
      call re_associate (li, nl+size(ele%wake%sr_trans%mode)+100, .false.)
      nl=nl+1; li(nl) = '  Short-Range Transverse Pseudo Modes:'
      nl=nl+1; li(nl) = &
            '   #        Amp        Damp           K         Phi   Polarization  Transverse_Dependence'
      do i = 1, size(ele%wake%sr_trans%mode)
        mode => ele%wake%sr_trans%mode(i)
        nl=nl+1; write (li(nl), '(i4, 4es12.4, a15, a16)') i, mode%amp, mode%damp, mode%k, mode%phi, &
                  sr_polarization_name(mode%polarization), sr_transverse_dependence_name(mode%transverse_dependence)
      enddo
    else
     nl=nl+1; li(nl) = '  No short-range transverse pseudo modes.'
    endif
  endif

  if (logic_option (.true., type_wake) .and. size(ele%wake%lr_mode) /= 0) then
    nl=nl+1; li(nl) = ''
    nl=nl+1; li(nl) = 'Long-Range Wake:'
    nl=nl+1; li(nl) = '  LR_File: ' // trim(ele%wake%lr_file)
    nl=nl+1; write (li(nl), '(2x, a, f10.6)') 'lr_freq_spread  =', ele%wake%lr_freq_spread
    nl=nl+1; write (li(nl), '(2x, a, l2)')    'lr_self_wake_on =', ele%wake%lr_self_wake_on
  endif

  if (size(ele%wake%lr_mode) /= 0) then
    nl=nl+1; write (li(nl), *)
    if (logic_option (.true., type_wake)) then
      call re_associate (li, nl+size(ele%wake%lr_mode)+100, .false.)
      nl=nl+1; li(nl) = '  Long-range HOM modes:'
      nl=nl+1; li(nl) = &
            '  #       Freq         R/Q           Q   m   Angle    b_sin     b_cos     a_sin     a_cos     t_ref'
      do i = 1, size(ele%wake%lr_mode)
        lr => ele%wake%lr_mode(i)
        angle = ' unpolar'
        if (lr%polarized) write (angle, '(f8.3)') lr%angle
        nl=nl+1; write (li(nl), '(i3, 3es12.4, i3, a, 5es10.2)') i, &
                lr%freq, lr%R_over_Q, lr%Q, lr%m, angle, &
                lr%b_sin, lr%b_cos, lr%a_sin, lr%a_cos, lr%t_ref
      enddo
    else
      nl=nl+1; li(nl) = '  No long-range HOM modes.'
    endif
  endif

endif

! Encode Floor coords

if (logic_option(.false., type_floor_coords)) then
  select case (ele%key)
  case (floor_shift$, group$, overlay$, hybrid$, beginning_ele$, match$, null_ele$, patch$)
    ! These elements do not have offsets
    floor = ele%floor

  case (crystal$, mirror$, multilayer_mirror$)
    ! Misalignments are referenced to beginning of element
    call ele_geometry (ele%floor, ele, floor, -1.0_rp)
    floor = coords_relative_to_floor (floor, [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)], &
                                        ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$))
    call ele_geometry (floor, ele, floor)

  case (girder$)
    floor = coords_relative_to_floor (ele%floor, [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)], &
                                        ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$)) 
  case default
    ! Misalignments referenced to center of element
    call ele_geometry (ele%floor, ele, floor, -0.5_rp)
    floor = coords_relative_to_floor (floor, [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)], &
                                        ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$)) 
    call ele_geometry (floor, ele, floor, 0.5_rp)
  end select

  nl=nl+1; li(nl) = ''
  nl=nl+1; li(nl) = 'Global Floor Coords at End of Element:'
  nl=nl+1; write (li(nl), '(a)')         '                   X           Y           Z       Theta         Phi         Psi'
  nl=nl+1; write (li(nl), '(a, 6f12.5, 3x, a)') 'Reference', ele%floor%r, ele%floor%theta, ele%floor%phi, ele%floor%psi, '! Position without misalignments'
  nl=nl+1; write (li(nl), '(a, 6f12.5, 3x, a)') 'Actual   ', floor%r, floor%theta, floor%phi, floor%psi, '! Position with offset/pitch/tilt misalignments'
  ele0 => pointer_to_next_ele(ele, -1)
  if (associated(ele0)) then
    if (ele%ix_ele /= 0 .or. ele%branch%param%geometry == closed$) then
      f0 = ele0%floor
      nl=nl+1; write (li(nl), '(a, 6f12.5, 3x, a)') 'delta Ref', floor%r-f0%r, floor%theta-f0%theta, floor%phi-f0%phi, floor%psi-f0%psi, &
                                                                                                         '! Delta with respect to last element'  
    endif
  endif
endif

! finish

if (present(lines)) then
  call re_allocate(lines, nl, .false.)
  n_lines = nl
  lines(1:nl) = li(1:nl)
else
  do i = 1, nl
    print '(1x, a)', trim(li(i))
  enddo
endif

end subroutine
