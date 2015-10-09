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
use geometry_mod, only: pointer_to_indexed_attribute, coords_relative_to_floor, ele_geometry
use expression_mod

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord, slave
type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (floor_position_struct) :: floor
type (wig_term_struct), pointer :: term
type (wake_lr_struct), pointer :: lr
type (wake_sr_mode_struct), pointer :: mode
type (em_field_mode_struct), pointer :: rfm
type (wall3d_struct), pointer :: wall3d
type (wall3d_section_struct), pointer :: section
type (wall3d_vertex_struct), pointer :: v
type (photon_element_struct), pointer :: p
type (photon_surface_struct), pointer :: s
type (ele_attribute_struct) attrib
type (lat_param_struct) param

integer, optional, intent(in) :: type_mat6, twiss_out
integer, optional, intent(out) :: n_lines
integer i, i1, j, n, ix, ix_tot, iv, ic, nl2, l_status, a_type, default_val
integer nl, nt, n_term, n_att, attrib_type, n_char, iy, particle

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
logical type_zero, err_flag, print_it, is_default, has_nonzero_pole

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

nl=nl+1; write (li(nl), '(1x, a, f15.6)')  'S:       ', ele%s
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
  if (ix_tot > 0) then
    if (ele%value(i) == 0 .and. ele%value(ix_tot) == 0 .and. .not. type_zero) cycle
    nl=nl+1; write (li(nl), '(i6, 3x, 2a, es15.7, a, i7, 3x, a16, a, es15.7)') &
                      i, a_name(1:n_att), '=', ele%value(i), ',', &
                      ix_tot, attribute_name(ele, ix_tot), '=', ele%value(ix_tot)

  elseif (a_name == 'P0C_START') then
    nl=nl+1; write (li(nl), '(i6, 3x, 2a, es15.7, a, 10x, a, f13.9)') &
                      i, a_name(1:n_att), '=', ele%value(i), ',', &
                      'BETA_START      =', ele%value(p0c_start$) / ele%value(e_tot_start$)

  elseif (a_name == 'E_TOT_START') then
    nl=nl+1; write (li(nl), '(i6, 3x, 2a, es15.7, a, 10x, a, es15.7)') &
                      i, a_name(1:n_att), '=', ele%value(i), ',', &
                      'DELTA_E         =', ele%value(e_tot$) - ele%value(e_tot_start$)

  elseif (a_name == 'P0C') then
    if (particle == photon$) then
      nl=nl+1; write (li(nl), '(i6, 3x, 2a, es15.7, a, 10x, a, es14.6)') &
                      i, a_name(1:n_att), '=', ele%value(i), ',', &
                      'REF_WAVELENGTH  =', c_light * h_planck / ele%value(p0c$)
    else
      nl=nl+1; write (li(nl), '(i6, 3x, 2a, es15.7, a, 10x, a, f13.9)') &
                      i, a_name(1:n_att), '=', ele%value(i), ',', &
                      'BETA            =', ele%value(p0c$) / ele%value(e_tot$)
    endif

  elseif (a_name == 'E_TOT') then
    nl=nl+1; write (li(nl), '(i6, 3x, 2a, es15.7, a, 10x, a, es15.7)') &
                      i, a_name(1:n_att), '=', ele%value(i), ',', &
                      'GAMMA           =', ele%value(e_tot$) / mass_of(particle)

  elseif (index(a_name, 'ANGLE') /= 0 .and. a_name /= 'CRITICAL_ANGLE_FACTOR') then
    units = ' deg]'
    if (a_name == 'DBRAGG_ANGLE_DE') units = ' deg/eV]'
    if (.not. type_zero .and. ele%value(i) == 0) cycle
    nl=nl+1; write (li(nl), '(i6, 3x, 2a, es15.7, 6x, a, f10.4, a)') &
                 i, a_name(1:n_att), '=', ele%value(i), '[', ele%value(i) * 180 / pi, trim(units)
  else
    attrib_type = attribute_type(a_name)
    if (is_a_tot_attribute(ele, i)) cycle
    select case (attrib_type)
    case (is_logical$)
      if (ele%value(i) /= 0) ele%value(i) = 1
      nl=nl+1; write (li(nl), '(i6, 3x, 2a, l1, a, i0, a)')  i, a_name(1:n_att), '=  ', &
                                  is_true(ele%value(i)), ' (', nint(ele%value(i)), ')'
    case (is_integer$)
      if (ele%value(i) == 0 .and. .not. type_zero) cycle
      nl=nl+1; write (li(nl), '(i6, 3x, 2a, i0)')  i, a_name(1:n_att), '= ', nint(ele%value(i))
    case (is_real$)
      if (ele%value(i) == 0 .and. .not. type_zero) cycle
      nl=nl+1; write (li(nl), '(i6, 3x, 2a, es15.7)')  i, a_name(1:n_att), '=', ele%value(i)
    case (is_switch$)
      name = switch_attrib_value_name (a_name, ele%value(i), ele, is_default)
      if (.not. is_default .or. type_zero) then
        nl=nl+1; write (li(nl), '(i6, 3x, 4a, i0, a)')  i, a_name(1:n_att), '=  ', &
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
    call multipole_ele_to_kt (ele, .false., has_nonzero_pole, a,  b)
    call multipole_ele_to_kt (ele, .true.,  has_nonzero_pole, knl, tn)
  else
    call multipole_ele_to_ab (ele, .false., has_nonzero_pole, a,  b)
    call multipole_ele_to_ab (ele, .true.,  has_nonzero_pole, a2, b2)
    call multipole_ele_to_kt (ele, .true.,  has_nonzero_pole, knl, tn)
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
    if (ele%a_pole_elec(i) == 0 .and. ele%b_pole_elec(i) == 0) cycle
    nl=nl+1; write (li(nl), '(2x, 4(3x, a, i0, a, es11.3))') 'A', i, '_elec =', ele%a_pole_elec(i)
    nl=nl+1; write (li(nl), '(2x, 4(3x, a, i0, a, es11.3))') 'B', i, '_elec =', ele%b_pole_elec(i)
  enddo

endif

! Encode on/off status etc.

if (.not. ele%is_on) then
  nl=nl+1; write (li(nl), *) '*** Note: Element is turned OFF ***'
endif

! wiggler terms

if (associated(ele%wig)) then
  if (logic_option(.false., type_field)) then
    nl=nl+1; write (li(nl), '(a, 6x, a, 3(9x, a), 3(9x, a))') ' Term#', &
                                'Coef', 'K_x', 'K_y', 'K_z', 'x0', 'y0', 'phi_z   Type'
    do i = 1, size(ele%wig%term)
      term => ele%wig%term(i)
      write (li(nl+i), '(i4, 4f12.6, 3f14.6, 3x, a)') i, term%coef, term%kx, term%ky, term%kz, &
                                 term%x0, term%y0, term%phi_z, wig_term_type_name(term%type)
    enddo
    nl = nl + size(ele%wig%term)
  else
    nl=nl+1; write (li(nl), '(a, i5)') 'Number of wiggler terms:', size(ele%wig%term)
  endif
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

if (ele%key /= overlay$ .and. ele%key /= group$) then
  if (ele%orientation /= 1 .or. type_zero) then
    nl=nl+1; write (li(nl), fmt_i) 'LONGITUDINAL ORIENTATION', '=', ele%orientation
  endif
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

! RF field coefs

if (associated(ele%em_field)) then
  if (logic_option(.false., type_field)) then
    do i = 1, size(ele%em_field%mode)
      rfm => ele%em_field%mode(i)
      if (rfm%master_scale == 0) then
        name = 'None'
      else
        name = attribute_name(ele, rfm%master_scale)
      endif
      nl=nl+1; li(nl) = ''
      nl=nl+1; write (li(nl), '(a, i0)')        'Mode #:', i
      nl=nl+1; write (li(nl), '(a, i0)')        '    m:             ', rfm%m
      nl=nl+1; write (li(nl), '(a, i0)')        '    harmonic:      ', rfm%harmonic
      nl=nl+1; write (li(nl), '(2a)')           '    master_scale:  ', trim(name)
      nl=nl+1; write (li(nl), '(a, es16.8)')    '    f_damp:        ', rfm%f_damp
      nl=nl+1; write (li(nl), '(a, es16.8)')    '    phi0_ref:     ', rfm%phi0_ref
      nl=nl+1; write (li(nl), '(a, es16.8)')    '    phi0_azimuth:  ', rfm%phi0_azimuth
      nl=nl+1; write (li(nl), '(a, es16.8)')    '    field_scale:   ', rfm%field_scale

      if (associated(rfm%map)) then
        nl=nl+1; write (li(nl), '(2a)')         '    File:          ', trim(rfm%map%file)
        nl=nl+1; write (li(nl), '(2a)')         '    ele_anchor_pt: ', anchor_pt_name(rfm%map%ele_anchor_pt)
        nl=nl+1; write (li(nl), '(a, i0)')      '    n_link:        ', rfm%map%n_link
        nl=nl+1; write (li(nl), '(a, es16.8)')  '    dz:            ', rfm%map%dz
        nl=nl+1; write (li(nl), '(a)')          '  Term                e                           b'
        do j = 1, min(10, size(rfm%map%term))
          if (nl+1 > size(li)) call re_associate(li, 2 * nl, .false.)
          nl=nl+1; write (li(nl), '(i5, 3x, 2(a, 2es12.4), a)') j, &
                                             '(', rfm%map%term(j)%e_coef, ')  (', rfm%map%term(j)%b_coef, ')'
        enddo
        if (size(rfm%map%term) > 10) then
          nl=nl+1; li(nl) = '     .... etc ...'
        endif
      endif

      if (associated(rfm%grid)) then
        nl=nl+1; write (li(nl), '(2a)')         '    File:          ', trim(rfm%grid%file)
        nl=nl+1; write (li(nl), '(2a)')         '    Type:          ', em_grid_type_name(rfm%grid%type)
        nl=nl+1; write (li(nl), '(2a)')         '    ele_anchor_pt: ', anchor_pt_name(rfm%grid%ele_anchor_pt)
        nl=nl+1; write (li(nl), '(a, i0)')      '    n_link:        ', rfm%grid%n_link
        nl=nl+1; write (li(nl), '(a, 3f14.6)')  '    dr:            ', rfm%grid%dr
        nl=nl+1; write (li(nl), '(a, 3f14.6)')  '    r0:            ', rfm%grid%r0
        nl=nl+1; write (li(nl), '(a, 3f14.6)')  '    r_max - r0:    ', ubound(rfm%grid%pt)*rfm%grid%dr - rfm%grid%r0
      endif

    enddo
  else
    nl=nl+1; write (li(nl), '(a, i5)') 'Number of Field modes:', size(ele%em_field%mode)
  endif
endif

! wall3d cross-sections.
! Do not print more than 100 sections.

wall3d => ele%wall3d
if (associated(wall3d)) then
  nl=nl+1; write (li(nl), '(a, i5)') ''
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
                  'Wall%Section(', i, '):  S =', section%s, trim(str1),  ',   (x0, y0) = (', &
                  section%x0, ',', section%y0, ')', trim(str2)


      do j = 1, size(section%v)
        v => section%v(j)
        nl=nl+1; write (li(nl), '(4x, a, i0, a, 5f11.6)') &
                              'v(', j, ') =', v%x, v%y, v%radius_x, v%radius_y, v%tilt
      enddo
    enddo
  endif
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
    nl=nl+1; write (li(nl), '(a)') 'Slave_status: BAD!', ele%slave_status
  else
    nl=nl+1; write (li(nl), '(2a)') 'Slave_status: ', control_name(ele%slave_status)
  endif

  if (ele%n_lord /= 0) then

    nl=nl+1; li(nl) = '   Index   Name                            Attribute           Lord_Type           Expression'

    do i = 1, ele%n_lord
      lord => pointer_to_lord (ele, i, ic)
      select case (lord%lord_status)
      case (super_lord$, multipass_lord$, girder_lord$)
        coef_str = ''
        a_name = ''
        val_str = ''
      case default
        if (allocated(lat%control(ic)%stack)) then
          coef_str = expression_stack_to_string (lat%control(ic)%stack)
        else
          coef_str = ''
        endif
        iv = lat%control(ic)%ix_attrib
        a_name = attribute_name(ele, iv)
      end select

      nl=nl+1; write (li(nl), '(i8, 3x, a32, a18, 2x, a20, a)') &
            lord%ix_ele, lord%name, a_name, key_name(lord%key), trim(coef_str)
    enddo
    nl=nl+1; li(nl) = ''

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
        nl=nl+1; write (li(nl), '(i6, 3x, 2a, es15.7, 11x, 3a, es15.7)')  i, &
                      a_name(1:n_att), '  =', ele%control_var(i)%value, &
                      'OLD_', a_name(1:n_att), '  =', ele%control_var(i)%old_value
      enddo
    else  ! overlay_lord
      do i = 1, size(ele%control_var)
        nl=nl+1; write (li(nl), '(i6, 3x, 2a, es15.7)')  i, &
                      ele%control_var(i)%name, '  =', ele%control_var(i)%value
      enddo
    endif
  endif

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
        nl=nl+1; write (li(nl), '(a8, t12, a, a)') &
                    trim(ele_loc_to_string(slave)), slave%name(1:n_char), key_name(slave%key)
      enddo

    case default
      nl=nl+1; li(nl) = '   Index   Name';  li(nl)(n_char+14:) = 'Attribute                Value    Expression'
      do ix = 1, ele%n_slave
        slave => pointer_to_slave (ele, ix, i)
        if (allocated(lat%control(i)%stack)) then
          coef_str = expression_stack_to_string (lat%control(i)%stack)
          if (ele%key == overlay$) then
            call evaluate_expression_stack(lat%control(i)%stack, val, err_flag, str1, ele%control_var)
            write (coef_str, '(es12.4, 4x, a)') val, trim(coef_str)
          endif
        else
          coef_str = ''
        endif
        iv = lat%control(i)%ix_attrib
        a_name = attribute_name(slave, iv)
        nl=nl+1; write (li(nl), '(a8, t12, a, 2x, a18, a)') trim(ele_loc_to_string(slave)), slave%name(1:n_char), a_name, trim(coef_str)
      enddo
    end select

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
      do i = 1, 6
        n_term = n_term + size(ele%taylor(i)%term)
      enddo
      nl=nl+1; write (li(nl), *) 'Taylor map total number of terms:', n_term
    endif
  endif

endif

! Encode HOM info

if (associated(ele%wake)) then

  if (size(ele%wake%sr_long%mode) /= 0) then
    nl=nl+1; write (li(nl), *)
    if (logic_option (.true., type_wake)) then
      call re_associate (li, nl+size(ele%wake%sr_long%mode)+100, .false.)
      nl=nl+1; li(nl) = 'Short-Range Longitudinal Pseudo Modes:'
      nl=nl+1; li(nl) = &
            '   #        Amp        Damp           K         Phi   Polarization  Transverse_Dependence'
      do i = 1, size(ele%wake%sr_long%mode)
        mode => ele%wake%sr_long%mode(i)
        nl=nl+1; write (li(nl), '(i4, 4es12.4, a15, a16)') i, mode%amp, mode%damp, mode%k, mode%phi, &
                  sr_polarization_name(mode%polarization), sr_transverse_dependence_name(mode%transverse_dependence)
      enddo
    else
      nl=nl+1; li(nl) = 'No short-range longitudinal pseudo modes.'
    endif
  endif

  if (size(ele%wake%sr_trans%mode) /= 0) then
    nl=nl+1; write (li(nl), *)
    if (logic_option (.true., type_wake)) then
      call re_associate (li, nl+size(ele%wake%sr_trans%mode)+100, .false.)
      nl=nl+1; li(nl) = 'Short-Range Transverse Pseudo Modes:'
      nl=nl+1; li(nl) = &
            '   #        Amp        Damp           K         Phi   Polarization  Transverse_Dependence'
      do i = 1, size(ele%wake%sr_trans%mode)
        mode => ele%wake%sr_trans%mode(i)
        nl=nl+1; write (li(nl), '(i4, 4es12.4, a15, a16)') i, mode%amp, mode%damp, mode%k, mode%phi, &
                  sr_polarization_name(mode%polarization), sr_transverse_dependence_name(mode%transverse_dependence)
      enddo
    else
     nl=nl+1; li(nl) = 'No short-range transverse pseudo modes.'
    endif
  endif

  if (size(ele%wake%lr) /= 0) then
    nl=nl+1; write (li(nl), *)
    if (logic_option (.true., type_wake)) then
      call re_associate (li, nl+size(ele%wake%lr)+100, .false.)
      nl=nl+1; li(nl) = 'Long-range HOM modes:'
      nl=nl+1; li(nl) = &
            '  #       Freq         R/Q           Q   m   Angle    b_sin     b_cos     a_sin     a_cos     t_ref'
      do i = 1, size(ele%wake%lr)
        lr => ele%wake%lr(i)
        angle = ' unpolar'
        if (lr%polarized) write (angle, '(f8.3)') lr%angle
        nl=nl+1; write (li(nl), '(i3, 3es12.4, i3, a, 5es10.2)') i, &
                lr%freq, lr%R_over_Q, lr%Q, lr%m, angle, &
                lr%b_sin, lr%b_cos, lr%a_sin, lr%a_cos, lr%t_ref
      enddo
    else
      nl=nl+1; li(nl) = 'No long-range HOM modes.'
    endif
  endif

endif

! Encode Floor coords

if (logic_option(.false., type_floor_coords)) then
  select case (ele%key)
  case (floor_shift$, group$, overlay$, hybrid$, beginning_ele$, match$, null_ele$, patch$)
    floor = ele%floor
  case (crystal$, mirror$, multilayer_mirror$)
    call ele_geometry (ele%floor, ele, floor, -1.0_rp)
    floor = coords_relative_to_floor (floor, [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)], &
                                        ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$)) 
    call ele_geometry (floor, ele, floor)
  case (girder$)
    floor = coords_relative_to_floor (ele%floor, [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)], &
                                        ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$)) 
  case default
    call ele_geometry (ele%floor, ele, floor, -0.5_rp)
    floor = coords_relative_to_floor (floor, [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)], &
                                        ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$)) 
    call ele_geometry (floor, ele, floor, 0.5_rp)
  end select
  nl=nl+1; li(nl) = ''
  nl=nl+1; li(nl) = 'Global Floor Coords:'
  nl=nl+1; write (li(nl), '(a)')         '                   X           Y           Z       Theta         Phi         Psi'
  nl=nl+1; write (li(nl), '(a, 6f12.5)') 'Reference', ele%floor%r, ele%floor%theta, ele%floor%phi, ele%floor%psi   
  nl=nl+1; write (li(nl), '(a, 6f12.5)') 'Actual   ', floor%r, floor%theta, floor%phi, floor%psi   
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
