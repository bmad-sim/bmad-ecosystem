module basic_attribute_mod

use equal_mod

implicit none

! attrib_array(key, ix_param)%type gives the type of the attribute.
! This may be one of:
!   does_not_exist$ -- Does not exist.
!   is_free$        -- Free to vary as long as attribute has not controlled by another element (overlay, super_lord, etc.)
!   quasi_free$     -- May be free or not. For example, k1 is only free if field_master = F.
!   dependent$      -- Value calculated by Bmad. Cannot be user varied as an independent parameter.
!   private$        -- Internal parameter used in calculations. Will not be displayed by type_ele.

integer, parameter :: does_not_exist$ = -1, is_free$ = 1, quasi_free$ = 2, dependent$ = 3, private$ = 4

type ele_attribute_struct
  character(40) :: name = null_name$
  integer :: type = does_not_exist$
end type

type (ele_attribute_struct), private, save :: attrib_array(n_key$, num_ele_attrib_extended$)

character(40), private, save :: short_attrib_array(n_key$, num_ele_attrib_extended$)
integer, private, save :: attrib_num(n_key$)
integer, private, save :: attrib_ix(n_key$, num_ele_attrib_extended$)
logical, private, save :: attribute_array_init_needed = .true.
logical, private, save :: has_orientation_attributes_key(n_key$)
private init_short_attrib_array

!+             
! Function attribute_index (...) result (attrib_index)
!
! Function to return the index of a attribute for a given BMAD element type
! and the name of the attribute. Abbreviations are by default permitted but must be at 
! least 3 characters. Exception: overlay and group varialbe names may not
! be abbreviated.
!
! This routine is an overloaded name for:
!   attribute_index1 (ele, name, full_name, can_abbreviate) result (attrib_index)
!   attribute_index2 (key, name, full_name, can_abbreviate) result (attrib_index)
!
! Note:
!   If ele%key or key = 0 -> Entire name table will be searched.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele     -- ele_struct: attribute_index will restrict the name search to 
!                valid attributes of the given element. 
!   key     -- integer: Equivalent to ele%key.
!   name    -- character(40): Attribute name. Must be uppercase.
!   can_abbreviate
!           -- logical, optional: Can abbreviate names? Default is True.
!
! Output:
!   full_name    -- character(40), optional: Non-abbreviated name.
!   attrib_index -- integer: Index of the attribute. If the attribute name
!                            is not appropriate then 0 will be returned.
!
! Example:
!     ele%key = sbend$
!     ix = attribute_index (ele, 'K1')
! Result:
!     ix -> k1$
!-

interface attribute_index
  module procedure attribute_index1
  module procedure attribute_index2
end interface

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_name (...) result (attrib_name)
!
! Function attribute_name1 (key, ix_att) result (attrib_name)
!
! Function to return the name of an attribute for a particular type of 
! Bmad element. 
!
! This routine is an overloaded name for:
!   attribute_name1 (ele, ix_att) result (attrib_name)
!   attribute_name2 (key, ix_att) result (attrib_name)
!
!
! Note: attribute_name (key, ix_att) is not able to handle overlay/group control variables.
! Use attributge_name (ele, ix_att) is this is needed.
!
! Input:
!   ele    -- Ele_struct: 
!     %key    -- Integer: Key name of element type (e.g. SBEND$, etc.)
!   key    -- Integer: Key name of element type (e.g. sbend$, etc.)
!   ix_att -- Integer: Index of attribute (e.g. k1$)
!
! Output:
!   attrib_name -- Character(40): Name of attribute. 
!      = "!BAD ELE KEY"                 %key is invalid
!      = "!BAD INDEX"                   ix_att is invalid (out of range).
!      = "!NULL" (null_name$)           ix_att does not correspond to an attribute or is private.
!
! Example:
!   ele%key = sbend$
!   name = attribute_name (ele, k1$)
! Result:
!   name -> "K1"
!-

interface attribute_name
  module procedure attribute_name1
  module procedure attribute_name2
end interface


contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+             
! Function attribute_index1 (ele, name, full_name, can_abbreviate) result (attrib_index)
!
! Overloaded by attribute_index. See attribute_index for more details.
!-

function attribute_index1 (ele, name, full_name, can_abbreviate) result (attrib_index)

type (ele_struct) ele
integer attrib_index, i, n
character(*) name
character(*), optional :: full_name
logical, optional :: can_abbreviate

! Note: ele%control_var may not be associated during parsing.

if ((ele%key == group$ .or. ele%key == overlay$) .and. associated(ele%control_var)) then
  n = min(4, len(name))
  if (name(1:n) == 'OLD_') then
    do i = 1, size(ele%control_var)
      if (full_name(5:) /= ele%control_var(i)%name) cycle
      attrib_index = i + old_control_var_offset$
      if (present(full_name)) full_name = name
      return
    enddo
  else  
    do i = 1, size(ele%control_var)
      if (name /= ele%control_var(i)%name) cycle
      attrib_index = i + var_offset$
      if (present(full_name)) full_name = name
      return
    enddo
  endif
endif

!

attrib_index = attribute_index2 (ele%key, name, full_name, can_abbreviate)

end function attribute_index1

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+             
! Function attribute_index2 (key, name, full_name, can_abbreviate) result (attrib_index)
!
! Overloaded by attribute_index. See attribute_index for more details.
!-

function attribute_index2 (key, name, full_name, can_abbreviate) result (attrib_index)

integer i, j, k, key, num, ilen, n_abbrev, ix_abbrev
integer attrib_index

character(*) name
character(*), optional :: full_name
character(40) name40
character(*), parameter :: r_name = 'attribute_index2'

logical, optional :: can_abbreviate

if (attribute_array_init_needed) call init_attribute_name_array

attrib_index = 0        ! match not found
if (present(full_name)) full_name = ''

ilen = len_trim(name)
if (ilen == 0) return
if (ilen < 3) ilen = 3  ! Need at least three characters.
n_abbrev = 0            ! number of abbreviation matches.

name40 = name           ! make sure we have 40 characters
if (ilen > 3 .and. name40(1:2) == 'TT') then
  do i = 3, ilen
    if (index('123456', name(i:i)) == 0) return
  enddo
  name40 = 'TT*'
endif

!-----------------------------------------------------------------------
! search for name

if (key == 0) then
  do k = 1, n_key$
    do i = 1, attrib_num(k)

      if (short_attrib_array(k, i) == name40) then
        attrib_index = attrib_ix(k, i)
        if (present(full_name)) then
          full_name = short_attrib_array(k, i)
          if (name40 == 'TT*') full_name = name
        endif
        return
      endif

      if (logic_option(.true., can_abbreviate)) then
        if (short_attrib_array(k, i)(1:ilen) == name40(1:ilen)) then
          n_abbrev = n_abbrev + 1
          ix_abbrev = attrib_ix(k, i)
          if (present(full_name)) full_name = short_attrib_array(k, i)
        endif 
      endif

    enddo
  enddo

! else only search this type of element

elseif (key > 0 .and. key <= n_key$) then
  do i = 1, attrib_num(key)

    if (short_attrib_array(key, i) == name40) then
      attrib_index = attrib_ix(key, i)
      if (present(full_name)) then
        full_name = short_attrib_array(key, i)
          if (name40 == 'TT*') full_name = name
        endif
      return
    endif

    if (logic_option(.true., can_abbreviate)) then
      if (short_attrib_array(key, i)(1:ilen) == name40(1:ilen)) then
        n_abbrev = n_abbrev + 1
        ix_abbrev = attrib_ix(key, i)
        if (present(full_name)) full_name = short_attrib_array(key, i)
      endif 
    endif

  enddo      

! error

else
  call out_io (s_fatal$, r_name, 'BAD KEY \i0\ ', i_array = [key])
  if (global_com%exit_on_error) call err_exit
endif

! If there is one unique abbreviation then use it.

if (n_abbrev == 1) attrib_index = ix_abbrev

end function attribute_index2

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_name1 (key, ix_att) result (attrib_name)
!
! Overloaded by attribute_name. See attribute_name for more details.
!-

function attribute_name1 (key, ix_att) result (attrib_name)

type (ele_struct) ele
integer i, key, ix_att, ix
character(40) attrib_name

!

ele%key = key
attrib_name = attribute_name2(ele, ix_att)

end function attribute_name1

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_name2 (ele, ix_att) result (attrib_name)
!
! Overloaded by attribute_name. See attribute_name for more details.
!-

function attribute_name2 (ele, ix_att) result (attrib_name)

type (ele_struct) ele
integer i, key, ix_att, ix
character(40) attrib_name

!

if (attribute_array_init_needed) call init_attribute_name_array()

key = ele%key

if (key <= 0 .or. key > n_key$) then
  attrib_name = '!BAD ELE KEY'

elseif ((key == group$ .or. key == overlay$) .and. is_attribute(ix_att, control_var$)) then
  ix = ix_att - var_offset$
  if (ix > size(ele%control_var)) then
    attrib_name = '!BAD INDEX'
  else
    attrib_name = ele%control_var(ix)%name
  endif

elseif (key == group$ .and. is_attribute(ix_att, old_control_var$)) then
  ix = ix_att - old_control_var_offset$
  if (ix > size(ele%control_var)) then
    attrib_name = '!BAD INDEX'
  else
    attrib_name = 'OLD_' // ele%control_var(ix)%name
  endif

elseif (ix_att <= 0 .or. ix_att > num_ele_attrib_extended$) then
  attrib_name = '!BAD INDEX'

else
  if (attrib_array(key, ix_att)%type == private$) then
    attrib_name = null_name$
  else
    attrib_name = attrib_array(key, ix_att)%name
  endif
endif

end function attribute_name2

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_info (ele, ix_att) result (attrib_info)
!
! Function to return the info structure associated with an attribute for 
! a particular type of BMAD element. 
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: 
!     %key    -- Integer: Key name of element type (e.g. SBEND$, etc.)
!   ix_att -- Integer: Index of attribute (e.g. k1$)
!
! Output:
!   attrib_info -- ele_attribute_struct: Info on this attribute.
!     %name -- Character(40): Name of attribute. 
!        = "!BAD ELE KEY"                 %key is invalid
!        = "!BAD INDEX"                   ix_att is invalid (out of range).
!        = "!INVALID INDEX"               ix_att is invalid for an overlay 
!        = "!NULL" (null_name$)           ix_att does not correspond to an attribute.
!
! Example:
!   ele%key = sbend$
!   name = attribute_name (ele, k1$)
! Result:
!   name -> "K1"
!-

function attribute_info (ele, ix_att) result (attrib_info)

type (ele_struct) ele
type (ele_attribute_struct) attrib_info
integer i, ix, key, ix_att

! Init

if (attribute_array_init_needed) call init_attribute_name_array()

attrib_info%type = does_not_exist$

! Overlay and group control vars.

if ((ele%key == group$ .or. ele%key == overlay$) .and. is_attribute(ix_att, control_var$)) then
  ix = ix_att - var_offset$
  if (ix > size(ele%control_var)) then
    attrib_info%name = '!BAD INDEX'
  else
    attrib_info%name = ele%control_var(ix)%name
    attrib_info%type = is_free$
  endif
  return
endif

if (ele%key == group$ .and. ix_att > old_control_var_offset$) then
  ix = ix_att - old_control_var_offset$
  if (ix > size(ele%control_var)) then
    attrib_info%name = '!BAD INDEX'
  else
    attrib_info%name = 'OLD_' // ele%control_var(ix)%name
    attrib_info%type = is_free$
  endif
  return
endif

! All else.

if (ele%key <= 0 .or. ele%key > n_key$) then
  attrib_info%name = '!BAD ELE KEY'
elseif (ix_att <= 0 .or. ix_att > num_ele_attrib_extended$) then
  attrib_info%name = '!BAD INDEX'
else
  attrib_info = attrib_array(ele%key, ix_att)
endif

end function attribute_info

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_attribute_name_array ()
!
! Private routine to initialize the attribute name array used by routines
! in attribute_mod. Not meant for general use.
!-

subroutine init_attribute_name_array ()

type (photon_surface_struct) surface
integer i, j, num, ix, iy
character(40) word

!

if (.not. attribute_array_init_needed) return

do i = 1, n_key$

  if (i == def_bmad_com$)   cycle
  if (i == def_mad_beam$)   cycle
  if (i == def_beam_start$) cycle
  if (i == line_ele$)       cycle
  if (i == def_parameter$)  cycle

  call init_attribute_name1 (i, custom_attribute1$,  'CUSTOM_ATTRIBUTE1', private$)
  call init_attribute_name1 (i, custom_attribute2$,  'CUSTOM_ATTRIBUTE2', private$)
  call init_attribute_name1 (i, custom_attribute3$,  'CUSTOM_ATTRIBUTE3', private$)
  call init_attribute_name1 (i, custom_attribute4$,  'CUSTOM_ATTRIBUTE4', private$)
  call init_attribute_name1 (i, custom_attribute5$,  'CUSTOM_ATTRIBUTE5', private$)
  call init_attribute_name1 (i, check_sum$, 'check_sum', private$)

  if (i == beginning_ele$)  cycle

  call init_attribute_name1 (i, type$,      'TYPE')
  call init_attribute_name1 (i, alias$,     'ALIAS')
  call init_attribute_name1 (i, descrip$,   'DESCRIP')

  if (i == group$)    cycle
  if (i == overlay$)  cycle
  if (i == girder$)   cycle

  call init_attribute_name1 (i, superimpose$,        'SUPERIMPOSE')
  call init_attribute_name1 (i, offset$,             'OFFSET')
  call init_attribute_name1 (i, reference$,          'REFERENCE')
  call init_attribute_name1 (i, ref_origin$,         'REF_ORIGIN')
  call init_attribute_name1 (i, ele_origin$,         'ELE_ORIGIN')

  if (i == null_ele$) cycle

  call init_attribute_name1 (i, mat6_calc_method$,       'MAT6_CALC_METHOD')
  call init_attribute_name1 (i, tracking_method$,        'TRACKING_METHOD')
  call init_attribute_name1 (i, spin_tracking_method$,   'SPIN_TRACKING_METHOD')
  call init_attribute_name1 (i, ptc_integration_type$,   'PTC_INTEGRATION_TYPE')

  call init_attribute_name1 (i, E_tot$,                  'E_TOT', quasi_free$) ! Free in multipass_lord
  call init_attribute_name1 (i, p0c$,                    'P0C', quasi_free$)   ! Free in multipass_lord
  call init_attribute_name1 (i, delta_ref_time$,         'DELTA_REF_TIME', dependent$)
  call init_attribute_name1 (i, ref_time_start$,         'ref_time_start', private$)

  if (i == fiducial$) cycle

  call init_attribute_name1 (i, create_jumbo_slave$,     'CREATE_JUMBO_SLAVE')

  call init_attribute_name1 (i, x_limit$,                'X_LIMIT')
  call init_attribute_name1 (i, x1_limit$,               'X1_LIMIT')
  call init_attribute_name1 (i, x2_limit$,               'X2_LIMIT')
  call init_attribute_name1 (i, y_limit$,                'Y_LIMIT')
  call init_attribute_name1 (i, y1_limit$,               'Y1_LIMIT')
  call init_attribute_name1 (i, y2_limit$,               'Y2_LIMIT')
  call init_attribute_name1 (i, aperture$,               'APERTURE')
  call init_attribute_name1 (i, aperture_at$,            'APERTURE_AT')
  call init_attribute_name1 (i, aperture_type$,          'APERTURE_TYPE')
  call init_attribute_name1 (i, offset_moves_aperture$,  'OFFSET_MOVES_APERTURE')

  if (i == hybrid$)   cycle
  if (i == match$) cycle
  if (i == photon_fork$)   cycle
  if (i == fork$)          cycle

  call init_attribute_name1 (i, x_offset$,      'X_OFFSET')
  call init_attribute_name1 (i, y_offset$,      'Y_OFFSET')
  call init_attribute_name1 (i, z_offset$,      'Z_OFFSET')
  call init_attribute_name1 (i, x_pitch$,       'X_PITCH')
  call init_attribute_name1 (i, y_pitch$,       'Y_PITCH')
  call init_attribute_name1 (i, tilt$,          'TILT' )
  call init_attribute_name1 (i, wall$,          'WALL' )

  if (i == floor_shift$)  cycle
  if (i == patch$)        cycle

  call init_attribute_name1 (i, tilt_tot$,      'TILT_TOT', dependent$)
  call init_attribute_name1 (i, x_offset_tot$,  'X_OFFSET_TOT', dependent$)
  call init_attribute_name1 (i, y_offset_tot$,  'Y_OFFSET_TOT', dependent$)
  call init_attribute_name1 (i, z_offset_tot$,  'Z_OFFSET_TOT', dependent$)
  call init_attribute_name1 (i, x_pitch_tot$,   'X_PITCH_TOT', dependent$)
  call init_attribute_name1 (i, y_pitch_tot$,   'Y_PITCH_TOT', dependent$)

  if (i == multilayer_mirror$) cycle
  if (i == mirror$)            cycle
  if (i == crystal$)           cycle
  if (i == sample$)            cycle
  if (i == capillary$)         cycle
  if (i == photon_init$)       cycle

  if (i /= drift$) call init_attribute_name1 (i, is_on$,        'IS_ON')

  if (i == diffraction_plate$) cycle
  if (i == mask$)              cycle
  if (i == marker$)            cycle
  if (i == detector$)          cycle
  if (i == beambeam$)          cycle
  if (i == multipole$)         cycle 
  if (i == ab_multipole$)      cycle

  call init_attribute_name1 (i, l$,                   'L')

  call init_attribute_name1 (i, symplectify$,         'SYMPLECTIFY')
  call init_attribute_name1 (i, taylor_map_includes_offsets$,    'TAYLOR_MAP_INCLUDES_OFFSETS')

  if (i == sad_mult$)      cycle

  call init_attribute_name1 (i, lord_pad1$,           'LORD_PAD1', quasi_free$)
  call init_attribute_name1 (i, lord_pad2$,           'LORD_PAD2', quasi_free$)

  if (i == taylor$)       cycle

  call init_attribute_name1 (i, start_edge$,          'START_EDGE', dependent$)
  call init_attribute_name1 (i, end_edge$,            'END_EDGE', dependent$)
  call init_attribute_name1 (i, accordion_edge$,      'ACCORDION_EDGE', dependent$)

  call init_attribute_name1 (i, integrator_order$,   'INTEGRATOR_ORDER')
  call init_attribute_name1 (i, num_steps$,          'NUM_STEPS', quasi_free$)
  call init_attribute_name1 (i, ds_step$,            'DS_STEP')
  call init_attribute_name1 (i, csr_calc_on$,        'CSR_CALC_ON')
  call init_attribute_name1 (i, n_ref_pass$,         'N_REF_PASS')

  if (i == drift$)        cycle

  call init_attribute_name1 (i, l_hard_edge$,        'L_HARD_EDGE', dependent$)
  call init_attribute_name1 (i, sr_wake_file$,       'SR_WAKE_FILE')

  call init_attribute_name1 (i, field_overlaps$,     'FIELD_OVERLAPS')

  if (i == pipe$)         cycle
  if (i == custom$)       cycle

  call init_attribute_name1 (i, fringe_type$,        'FRINGE_TYPE')
  call init_attribute_name1 (i, spin_fringe_on$,     'SPIN_FRINGE_ON')
  call init_attribute_name1 (i, fringe_at$,          'FRINGE_AT')

  if (i == hkicker$)      cycle
  if (i == vkicker$)      cycle
  if (i == e_gun$)        cycle
  if (i == em_field$)     cycle

  call init_attribute_name1 (i, hkick$,     'HKICK', quasi_free$)
  call init_attribute_name1 (i, vkick$,     'VKICK', quasi_free$)

  if (i == elseparator$) cycle

  call init_attribute_name1 (i, bl_hkick$,  'BL_HKICK', quasi_free$)
  call init_attribute_name1 (i, bl_vkick$,  'BL_VKICK', quasi_free$)
enddo

!

do i = 1, n_key$
  select case(i)
  case (monitor$, instrument$, marker$, detector$)
    call init_attribute_name1 (i, x_gain_err$,      'X_GAIN_ERR')
    call init_attribute_name1 (i, y_gain_err$,      'Y_GAIN_ERR')
    call init_attribute_name1 (i, crunch$,          'CRUNCH')
    call init_attribute_name1 (i, noise$,           'NOISE')
    call init_attribute_name1 (i, tilt_calib$,      'TILT_CALIB')
    call init_attribute_name1 (i, x_gain_calib$,    'X_GAIN_CALIB')
    call init_attribute_name1 (i, y_gain_calib$,    'Y_GAIN_CALIB')
    call init_attribute_name1 (i, crunch_calib$,    'CRUNCH_CALIB')
    call init_attribute_name1 (i, x_offset_calib$,  'X_OFFSET_CALIB')
    call init_attribute_name1 (i, y_offset_calib$,  'Y_OFFSET_CALIB')
    call init_attribute_name1 (i, n_sample$,        'N_SAMPLE')
    call init_attribute_name1 (i, de_eta_meas$,     'DE_ETA_MEAS')
    call init_attribute_name1 (i, osc_amplitude$,   'OSC_AMPLITUDE')
  end select
enddo

!

do i = 1, n_key$
  select case (i)
  case (crystal$, multilayer_mirror$, mirror$, sample$, diffraction_plate$, detector$)
    call init_attribute_name1 (i, surface_attrib$, 'SURFACE')
    num = a0$ - 1
    do ix = 0, ubound(surface%curvature_xy, 1)
    do iy = 0, ubound(surface%curvature_xy, 2)
      if (ix+iy < 2) cycle
      if (ix+iy > ubound(surface%curvature_xy, 1)) cycle
      write (word, '(a, i1, a, i1)') 'CURVATURE_X', ix, '_Y', iy
      num = num + 1
      call init_attribute_name1 (i, num, word) 
    enddo
    enddo
  end select
enddo

!

do i = 1, n_key$
  select case (i)
  case (elseparator$, kicker$, octupole$, quadrupole$, sbend$, rbend$, &
         sextupole$, solenoid$, sol_quad$, ab_multipole$, wiggler$, undulator$, bend_sol_quad$, &
         hkicker$, vkicker$, sad_mult$)
    attrib_array(i, a0$:a21$)%name = ['A0 ', &
                                   'A1 ', 'A2 ', 'A3 ', 'A4 ', 'A5 ', & 
                                   'A6 ', 'A7 ', 'A8 ', 'A9 ', 'A10', &
                                   'A11', 'A12', 'A13', 'A14', 'A15', &
                                   'A16', 'A17', 'A18', 'A19', 'A20', 'A21']
    attrib_array(i, b0$:b21$)%name = ['B0 ', &
                                   'B1 ', 'B2 ', 'B3 ', 'B4 ', 'B5 ', & 
                                   'B6 ', 'B7 ', 'B8 ', 'B9 ', 'B10', &
                                   'B11', 'B12', 'B13', 'B14', 'B15', &
                                   'B16', 'B17', 'B18', 'B19', 'B20', 'B21']
    attrib_array(i, a0$:b21$)%type = is_free$
    if (i == sad_mult$) cycle
    call init_attribute_name1 (i, multipoles_on$,     'MULTIPOLES_ON')

    if (i == multipole$) cycle
    if (i == ab_multipole$) cycle
    attrib_array(i, a0_elec$:a21_elec$)%name = ['A0_ELEC ', &
                                    'A1_ELEC ', 'A2_ELEC ', 'A3_ELEC ', 'A4_ELEC ', 'A5_ELEC ', & 
                                    'A6_ELEC ', 'A7_ELEC ', 'A8_ELEC ', 'A9_ELEC ', 'A10_ELEC', &
                                    'A11_ELEC', 'A12_ELEC', 'A13_ELEC', 'A14_ELEC', 'A15_ELEC', &
                                    'A16_ELEC', 'A17_ELEC', 'A18_ELEC', 'A19_ELEC', 'A20_ELEC', 'A21_ELEC']
    attrib_array(i, b0_elec$:b21_elec$)%name = ['B0_ELEC ', &
                                    'B1_ELEC ', 'B2_ELEC ', 'B3_ELEC ', 'B4_ELEC ', 'B5_ELEC ', & 
                                    'B6_ELEC ', 'B7_ELEC ', 'B8_ELEC ', 'B9_ELEC ', 'B10_ELEC', &
                                    'B11_ELEC', 'B12_ELEC', 'B13_ELEC', 'B14_ELEC', 'B15_ELEC', &
                                    'B16_ELEC', 'B17_ELEC', 'B18_ELEC', 'B19_ELEC', 'B20_ELEC', 'B21_ELEC']
    call init_attribute_name1 (i, scale_multipoles$,  'SCALE_MULTIPOLES')
  end select
enddo

attrib_array(hybrid$, a0$+1:a0$+6)%name = ['KICK1', 'KICK2', 'KICK3', 'KICK4', 'KICK5', 'KICK6']
attrib_array(hybrid$, a0_elec$+1:a0_elec$+36)%name = &
                                           ['MAT11', 'MAT12', 'MAT13', 'MAT14', 'MAT15', 'MAT16', &
                                            'MAT21', 'MAT22', 'MAT23', 'MAT24', 'MAT25', 'MAT26', &
                                            'MAT31', 'MAT32', 'MAT33', 'MAT34', 'MAT35', 'MAT36', &
                                            'MAT41', 'MAT42', 'MAT43', 'MAT44', 'MAT45', 'MAT46', &
                                            'MAT51', 'MAT52', 'MAT53', 'MAT54', 'MAT55', 'MAT56', &
                                            'MAT61', 'MAT62', 'MAT63', 'MAT64', 'MAT65', 'MAT66']

!----------------------------------------------------------
! Convention: Private attributes have names in lower case.

call init_attribute_name1 (photon_fork$, l$,                        'L', dependent$)
call init_attribute_name1 (photon_fork$, ix_to_branch$,             'IX_TO_BRANCH', dependent$)
call init_attribute_name1 (photon_fork$, ix_to_element$,            'IX_TO_ELEMENT', dependent$)
call init_attribute_name1 (photon_fork$, direction$,                'DIRECTION')
call init_attribute_name1 (photon_fork$, to_line$,                  'TO_LINE')
call init_attribute_name1 (photon_fork$, to_element$,               'TO_ELEMENT')
call init_attribute_name1 (photon_fork$, new_branch$,               'NEW_BRANCH')
call init_attribute_name1 (photon_fork$, is_on$,                    'IS_ON')

attrib_array(fork$, :) = attrib_array(photon_fork$, :)

call init_attribute_name1 (beambeam$, l$,                           'L', dependent$)
call init_attribute_name1 (beambeam$, sig_x$,                       'SIG_X')
call init_attribute_name1 (beambeam$, sig_y$,                       'SIG_Y')
call init_attribute_name1 (beambeam$, sig_z$,                       'SIG_Z')
call init_attribute_name1 (beambeam$, bbi_const$,                   'BBI_CONSTANT', dependent$)
call init_attribute_name1 (beambeam$, charge$,                      'CHARGE')
call init_attribute_name1 (beambeam$, n_slice$,                     'N_SLICE')
call init_attribute_name1 (beambeam$, symplectify$,                 'N_SLICE')
call init_attribute_name1 (beambeam$, field_calc$,                  'FIELD_CALC')
call init_attribute_name1 (beambeam$, field_master$,                'FIELD_MASTER')
call init_attribute_name1 (beambeam$, E_tot_start$,                 'E_tot_start', private$)
call init_attribute_name1 (beambeam$, beta_a$,                      'BETA_A')
call init_attribute_name1 (beambeam$, beta_b$,                      'BETA_B')
call init_attribute_name1 (beambeam$, alpha_a$,                     'ALPHA_A')
call init_attribute_name1 (beambeam$, alpha_b$,                     'ALPHA_B')
call init_attribute_name1 (beambeam$, cmat_11$,                     'CMAT_11')
call init_attribute_name1 (beambeam$, cmat_12$,                     'CMAT_12')
call init_attribute_name1 (beambeam$, cmat_21$,                     'CMAT_21')
call init_attribute_name1 (beambeam$, cmat_22$,                     'CMAT_22')

call init_attribute_name1 (beginning_ele$, l$,                           'l', dependent$)
call init_attribute_name1 (beginning_ele$, floor_set$,                   'floor_set', private$)
call init_attribute_name1 (beginning_ele$, delta_ref_time$,              'delta_ref_time', private$)
call init_attribute_name1 (beginning_ele$, ref_time_start$,              'ref_time_start', private$)
call init_attribute_name1 (beginning_ele$, e_tot_start$,                 'E_TOT_START', dependent$)
call init_attribute_name1 (beginning_ele$, p0c_start$,                   'P0C_START', dependent$)
call init_attribute_name1 (beginning_ele$, e_tot$,                       'E_TOT')
call init_attribute_name1 (beginning_ele$, p0c$,                         'P0C')
call init_attribute_name1 (beginning_ele$, x_position$,                  'X_POSITION')
call init_attribute_name1 (beginning_ele$, y_position$,                  'Y_POSITION')
call init_attribute_name1 (beginning_ele$, z_position$,                  'Z_POSITION')
call init_attribute_name1 (beginning_ele$, theta_position$,              'THETA_POSITION')
call init_attribute_name1 (beginning_ele$, phi_position$,                'PHI_POSITION')
call init_attribute_name1 (beginning_ele$, psi_position$,                'PSI_POSITION')
call init_attribute_name1 (beginning_ele$, beta_a_begin$,                'BETA_A')
call init_attribute_name1 (beginning_ele$, beta_b_begin$,                'BETA_B')
call init_attribute_name1 (beginning_ele$, alpha_a_begin$,               'ALPHA_A')
call init_attribute_name1 (beginning_ele$, alpha_b_begin$,               'ALPHA_B')
call init_attribute_name1 (beginning_ele$, eta_x$,                       'ETA_X')
call init_attribute_name1 (beginning_ele$, eta_y$,                       'ETA_Y')
call init_attribute_name1 (beginning_ele$, eta_z$,                       'ETA_Z')
call init_attribute_name1 (beginning_ele$, etap_x$,                      'ETAP_X')
call init_attribute_name1 (beginning_ele$, etap_y$,                      'ETAP_Y')
call init_attribute_name1 (beginning_ele$, phi_a$,                       'PHI_A')
call init_attribute_name1 (beginning_ele$, phi_b$,                       'PHI_B')
call init_attribute_name1 (beginning_ele$, cmat_11_begin$,               'CMAT_11')
call init_attribute_name1 (beginning_ele$, cmat_12_begin$,               'CMAT_12')
call init_attribute_name1 (beginning_ele$, cmat_21_begin$,               'CMAT_21')
call init_attribute_name1 (beginning_ele$, cmat_22_begin$,               'CMAT_22')
call init_attribute_name1 (beginning_ele$, s_long$,                      'S')
call init_attribute_name1 (beginning_ele$, ref_time$,                    'REF_TIME')

call init_attribute_name1 (beginning_ele$, spin_x$,                  'spin_x', private$)
call init_attribute_name1 (beginning_ele$, spin_y$,                  'spin_y', private$)
call init_attribute_name1 (beginning_ele$, spin_z$,                  'spin_z', private$)
call init_attribute_name1 (beginning_ele$, spinor_theta$,            'spinor_theta', private$)
call init_attribute_name1 (beginning_ele$, spinor_phi$,              'spinor_phi', private$)
call init_attribute_name1 (beginning_ele$, spinor_xi$,               'spinor_xi', private$)
call init_attribute_name1 (beginning_ele$, spinor_polarization$,     'spinor_polarization', private$)

attrib_array(line_ele$, :) = attrib_array(beginning_ele$, :)
call init_attribute_name1 (line_ele$, particle$,                    'PARTICLE')
call init_attribute_name1 (line_ele$, geometry$,                    'GEOMETRY')
call init_attribute_name1 (line_ele$, default_tracking_species$,    'DEFAULT_TRACKING_SPECIES')

call init_attribute_name1 (capillary$, l$,                          'L', dependent$)
call init_attribute_name1 (capillary$, s_spline$,                   'S_SPLINE')
call init_attribute_name1 (capillary$, n_slice_spline$,             'N_SLICE_SPLINE')
call init_attribute_name1 (capillary$, critical_angle_factor$,      'CRITICAL_ANGLE_FACTOR')
call init_attribute_name1 (capillary$, e_tot_start$,                'e_tot_start', private$)
call init_attribute_name1 (capillary$, p0c_start$,                  'p0c_start', private$)

call init_attribute_name1 (crystal$, l$,                            'L', dependent$)
call init_attribute_name1 (crystal$, bragg_angle_in$,               'BRAGG_ANGLE_IN', dependent$)
call init_attribute_name1 (crystal$, bragg_angle_out$,              'BRAGG_ANGLE_OUT', dependent$)
call init_attribute_name1 (crystal$, psi_angle$,                    'PSI_ANGLE')
call init_attribute_name1 (crystal$, alpha_angle$,                  'ALPHA_ANGLE', dependent$)
call init_attribute_name1 (crystal$, ref_tilt$,                     'REF_TILT')
call init_attribute_name1 (crystal$, ref_tilt_tot$,                 'REF_TILT_TOT', dependent$)
call init_attribute_name1 (crystal$, tilt_corr$,                    'TILT_CORR', dependent$)
call init_attribute_name1 (crystal$, d_spacing$,                    'D_SPACING')
call init_attribute_name1 (crystal$, v_unitcell$,                   'V_UNITCELL', dependent$)
call init_attribute_name1 (crystal$, b_param$,                      'B_PARAM')
call init_attribute_name1 (crystal$, bragg_angle$,                  'BRAGG_ANGLE' , dependent$)
call init_attribute_name1 (crystal$, ref_wavelength$,               'REF_WAVELENGTH', dependent$)
call init_attribute_name1 (crystal$, crystal_type$,                 'CRYSTAL_TYPE')
call init_attribute_name1 (crystal$, thickness$,                    'THICKNESS')
call init_attribute_name1 (crystal$, ref_orbit_follows$,            'REF_ORBIT_FOLLOWS')
call init_attribute_name1 (crystal$, e_tot_start$,                  'e_tot_start', private$)
call init_attribute_name1 (crystal$, p0c_start$,                    'p0c_start', private$)
call init_attribute_name1 (crystal$, ref_cap_gamma$,                'REF_CAP_GAMMA', dependent$)
call init_attribute_name1 (crystal$, darwin_width_sigma$,           'DARWIN_WIDTH_SIGMA', dependent$)
call init_attribute_name1 (crystal$, darwin_width_pi$,              'DARWIN_WIDTH_PI', dependent$)
call init_attribute_name1 (crystal$, pendellosung_period_sigma$,    'PENDELLOSUNG_PERIOD_SIGMA', dependent$)
call init_attribute_name1 (crystal$, pendellosung_period_pi$,       'PENDELLOSUNG_PERIOD_PI', dependent$)
call init_attribute_name1 (crystal$, dbragg_angle_de$,              'DBRAGG_ANGLE_DE', dependent$)
call init_attribute_name1 (crystal$, diffraction_limited$,          'DIFFRACTION_LIMITED')

call init_attribute_name1 (def_bmad_com$, max_aperture_limit$,         'MAX_APERTURE_LIMIT')
call init_attribute_name1 (def_bmad_com$, default_ds_step$,            'DEFAULT_DS_STEP')
call init_attribute_name1 (def_bmad_com$, significant_length$,         'SIGNIFICANT_LENGTH')
call init_attribute_name1 (def_bmad_com$, rel_tol_tracking$,           'REL_TOL_TRACKING')
call init_attribute_name1 (def_bmad_com$, abs_tol_tracking$,           'ABS_TOL_TRACKING')
call init_attribute_name1 (def_bmad_com$, rel_tol_adaptive_tracking$,  'REL_TOL_ADAPTIVE_TRACKING')
call init_attribute_name1 (def_bmad_com$, abs_tol_adaptive_tracking$,  'ABS_TOL_ADAPTIVE_TRACKING')
call init_attribute_name1 (def_bmad_com$, init_ds_adaptive_tracking$,  'INIT_DS_ADAPTIVE_TRACKING')
call init_attribute_name1 (def_bmad_com$, min_ds_adaptive_tracking$,   'MIN_DS_ADAPTIVE_TRACKING')
call init_attribute_name1 (def_bmad_com$, fatal_ds_adaptive_tracking$, 'FATAL_DS_ADAPTIVE_TRACKING')
call init_attribute_name1 (def_bmad_com$, aperture_limit_on$,          'APERTURE_LIMIT_ON')
call init_attribute_name1 (def_bmad_com$, max_num_runge_kutta_step$,   'MAX_NUM_RUNGE_KUTTA_STEP')

call init_attribute_name1 (def_mad_beam$, particle$,                      'PARTICLE')
call init_attribute_name1 (def_mad_beam$, e_tot$,                         'ENERGY')
call init_attribute_name1 (def_mad_beam$, p0c$,                           'PC')
call init_attribute_name1 (def_mad_beam$, n_part$,                        'N_PART')

call init_attribute_name1 (def_beam_start$, x$,                       'X')
call init_attribute_name1 (def_beam_start$, px$,                      'PX')
call init_attribute_name1 (def_beam_start$, y$,                       'Y')
call init_attribute_name1 (def_beam_start$, py$,                      'PY')
call init_attribute_name1 (def_beam_start$, z$,                       'Z')
call init_attribute_name1 (def_beam_start$, pz$,                      'PZ')
call init_attribute_name1 (def_beam_start$, field_x$,                 'FIELD_X')
call init_attribute_name1 (def_beam_start$, field_y$,                 'FIELD_Y')
call init_attribute_name1 (def_beam_start$, phase_x$,                 'PHASE_X')
call init_attribute_name1 (def_beam_start$, phase_y$,                 'PHASE_Y')
call init_attribute_name1 (def_beam_start$, t$,                       'T')
call init_attribute_name1 (def_beam_start$, e_photon$,                'E_PHOTON')
call init_attribute_name1 (def_beam_start$, spin_x$,                  'SPIN_X')
call init_attribute_name1 (def_beam_start$, spin_y$,                  'SPIN_Y')
call init_attribute_name1 (def_beam_start$, spin_z$,                  'SPIN_Z')
call init_attribute_name1 (def_beam_start$, spinor_theta$,            'SPINOR_THETA')
call init_attribute_name1 (def_beam_start$, spinor_phi$,              'SPINOR_PHI')
call init_attribute_name1 (def_beam_start$, spinor_xi$,               'SPINOR_XI')
call init_attribute_name1 (def_beam_start$, spinor_polarization$,     'SPINOR_POLARIZATION')
call init_attribute_name1 (def_beam_start$, emittance_a$,             'EMITTANCE_A')
call init_attribute_name1 (def_beam_start$, emittance_b$,             'EMITTANCE_B')
call init_attribute_name1 (def_beam_start$, emittance_z$,             'EMITTANCE_Z')
call init_attribute_name1 (def_beam_start$, sig_e$,                   'SIG_E')
call init_attribute_name1 (def_beam_start$, sig_z$,                   'SIG_Z')
call init_attribute_name1 (def_beam_start$, direction_beam_start$,    'DIRECTION')


call init_attribute_name1 (def_parameter$, custom_attribute1$,        'CUSTOM_ATTRIBUTE1')
call init_attribute_name1 (def_parameter$, custom_attribute2$,        'CUSTOM_ATTRIBUTE2')
call init_attribute_name1 (def_parameter$, custom_attribute3$,        'CUSTOM_ATTRIBUTE3')
call init_attribute_name1 (def_parameter$, custom_attribute4$,        'CUSTOM_ATTRIBUTE4')
call init_attribute_name1 (def_parameter$, custom_attribute5$,        'CUSTOM_ATTRIBUTE5')
call init_attribute_name1 (def_parameter$, e_tot$,                    'E_TOT')
call init_attribute_name1 (def_parameter$, p0c$,                      'P0C')
call init_attribute_name1 (def_parameter$, geometry$,                 'GEOMETRY')
call init_attribute_name1 (def_parameter$, lattice_type$,             'LATTICE_TYPE') ! For backwards compatibility
call init_attribute_name1 (def_parameter$, lattice$,                  'LATTICE')
call init_attribute_name1 (def_parameter$, taylor_order$,             'TAYLOR_ORDER')
call init_attribute_name1 (def_parameter$, ran_seed$,                 'RAN_SEED')
call init_attribute_name1 (def_parameter$, n_part$,                   'N_PART')
call init_attribute_name1 (def_parameter$, particle$,                 'PARTICLE')
call init_attribute_name1 (def_parameter$, photon_type$,              'PHOTON_TYPE')
call init_attribute_name1 (def_parameter$, no_end_marker$,            'NO_END_MARKER')
call init_attribute_name1 (def_parameter$, absolute_time_tracking$,   'ABSOLUTE_TIME_TRACKING')
call init_attribute_name1 (def_parameter$, ptc_exact_model$,          'PTC_EXACT_MODEL')
call init_attribute_name1 (def_parameter$, ptc_exact_misalign$,       'PTC_EXACT_MISALIGN')
call init_attribute_name1 (def_parameter$, default_tracking_species$, 'DEFAULT_TRACKING_SPECIES')
call init_attribute_name1 (def_parameter$, electric_dipole_moment$,   'ELECTRIC_DIPOLE_MOMENT')

call init_attribute_name1 (detector$, l$,                             'L', dependent$)
call init_attribute_name1 (detector$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (detector$, p0c_start$,                     'p0c_start', private$)

call init_attribute_name1 (diffraction_plate$, l$,                    'l', dependent$)
call init_attribute_name1 (diffraction_plate$, mode$,                 'MODE')
call init_attribute_name1 (diffraction_plate$, field_scale_factor$,   'FIELD_SCALE_FACTOR')
call init_attribute_name1 (diffraction_plate$, ref_wavelength$,       'REF_WAVELENGTH', dependent$)

call init_attribute_name1 (mask$, l$,                                 'l', dependent$)
call init_attribute_name1 (mask$, mode$,                              'MODE')
call init_attribute_name1 (mask$, field_scale_factor$,                'FIELD_SCALE_FACTOR')
call init_attribute_name1 (mask$, ref_wavelength$,                    'REF_WAVELENGTH', dependent$)

call init_attribute_name1 (drift$, E_tot_start$,                    'E_tot_start', private$)
call init_attribute_name1 (drift$, p0c_start$,                      'p0c_start', private$)
call init_attribute_name1 (drift$, spin_fringe_on$,                 'spin_fringe_on', private$)
call init_attribute_name1 (drift$, fringe_type$,                    'fringe_type', private$)
call init_attribute_name1 (drift$, fringe_at$,                      'fringe_at', private$)

call init_attribute_name1 (e_gun$, e_tot_ref_init$,                 'e_tot_ref_init', private$)
call init_attribute_name1 (e_gun$, p0c_ref_init$,                   'p0c_ref_init', private$)
call init_attribute_name1 (e_gun$, e_tot_start$,                    'e_tot_start', private$)
call init_attribute_name1 (e_gun$, p0c_start$,                      'p0c_start', private$)
call init_attribute_name1 (e_gun$, field_autoscale$,                'FIELD_AUTOSCALE', quasi_free$)
call init_attribute_name1 (e_gun$, autoscale_amplitude$,            'AUTOSCALE_AMPLITUDE')
call init_attribute_name1 (e_gun$, autoscale_phase$,                'AUTOSCALE_PHASE')
call init_attribute_name1 (e_gun$, voltage$,                        'VOLTAGE')
call init_attribute_name1 (e_gun$, voltage_err$,                    'VOLTAGE_ERR')
call init_attribute_name1 (e_gun$, gradient$,                       'GRADIENT')
call init_attribute_name1 (e_gun$, gradient_err$,                   'GRADIENT_ERR')
call init_attribute_name1 (e_gun$, field_calc$,                     'FIELD_CALC')
call init_attribute_name1 (e_gun$, cartesian_map$,                  'CARTESIAN_MAP')
call init_attribute_name1 (e_gun$, cylindrical_map$,                'CYLINDRICAL_MAP')
call init_attribute_name1 (e_gun$, grid_field$,                     'GRID_FIELD')
call init_attribute_name1 (e_gun$, taylor_field$,                   'TAYLOR_FIELD')
call init_attribute_name1 (e_gun$, rf_frequency$,                   'RF_FREQUENCY')
call init_attribute_name1 (e_gun$, phi0$,                           'PHI0')
call init_attribute_name1 (e_gun$, phi0_err$,                       'PHI0_ERR')
! e_gun attribute phi0_multipass should always be 0 and is used to make lcavity and e_gun equations similar
call init_attribute_name1 (e_gun$, phi0_multipass$,                 'phi0_multipass', private$) 
call init_attribute_name1 (e_gun$, phi0_autoscale$,                 'PHI0_AUTOSCALE', quasi_free$)

call init_attribute_name1 (elseparator$, gap$,                      'GAP')
call init_attribute_name1 (elseparator$, e_field$,                  'E_FIELD', quasi_free$)
call init_attribute_name1 (elseparator$, voltage$,                  'VOLTAGE', quasi_free$)
call init_attribute_name1 (elseparator$, r0_mag$,                   'R0_MAG')
call init_attribute_name1 (elseparator$, r0_elec$,                  'R0_ELEC')
call init_attribute_name1 (elseparator$, field_calc$,               'FIELD_CALC')
call init_attribute_name1 (elseparator$, field_master$,             'FIELD_MASTER')
call init_attribute_name1 (elseparator$, cartesian_map$,            'CARTESIAN_MAP')
call init_attribute_name1 (elseparator$, cylindrical_map$,          'CYLINDRICAL_MAP')
call init_attribute_name1 (elseparator$, grid_field$,               'GRID_FIELD')
call init_attribute_name1 (elseparator$, taylor_field$,             'TAYLOR_FIELD')
call init_attribute_name1 (elseparator$, E_tot_start$,              'E_tot_start', private$)
call init_attribute_name1 (elseparator$, p0c_start$,                'p0c_start', private$)

call init_attribute_name1 (em_field$, e_tot_start$,                 'E_TOT_START', dependent$)
call init_attribute_name1 (em_field$, p0c_start$,                   'P0C_START', dependent$)
call init_attribute_name1 (em_field$, cartesian_map$,               'CARTESIAN_MAP')
call init_attribute_name1 (em_field$, cylindrical_map$,             'CYLINDRICAL_MAP')
call init_attribute_name1 (em_field$, grid_field$,                  'GRID_FIELD')
call init_attribute_name1 (em_field$, taylor_field$,                'TAYLOR_FIELD')
call init_attribute_name1 (em_field$, field_calc$,                  'FIELD_CALC')
call init_attribute_name1 (em_field$, rf_frequency$,                'RF_FREQUENCY')
call init_attribute_name1 (em_field$, field_autoscale$,             'FIELD_AUTOSCALE', quasi_free$)
call init_attribute_name1 (em_field$, phi0_autoscale$,              'PHI0_AUTOSCALE', quasi_free$)
call init_attribute_name1 (em_field$, autoscale_amplitude$,         'AUTOSCALE_AMPLITUDE')
call init_attribute_name1 (em_field$, autoscale_phase$,             'AUTOSCALE_PHASE')
call init_attribute_name1 (em_field$, phi0$,                        'PHI0')
call init_attribute_name1 (em_field$, phi0_err$,                    'PHI0_ERR')

call init_attribute_name1 (girder$, l$,                             'L', dependent$)
call init_attribute_name1 (girder$, x_offset$,                      'X_OFFSET')
call init_attribute_name1 (girder$, y_offset$,                      'Y_OFFSET')
call init_attribute_name1 (girder$, z_offset$,                      'Z_OFFSET')
call init_attribute_name1 (girder$, x_pitch$,                       'X_PITCH')
call init_attribute_name1 (girder$, y_pitch$,                       'Y_PITCH')
call init_attribute_name1 (girder$, tilt$,                          'TILT')
call init_attribute_name1 (girder$, tilt_tot$,                      'TILT_TOT')
call init_attribute_name1 (girder$, ref_tilt$,                      'REF_TILT')
call init_attribute_name1 (girder$, ref_tilt_tot$,                  'REF_TILT_TOT')
call init_attribute_name1 (girder$, x_offset_tot$,                  'X_OFFSET_TOT')
call init_attribute_name1 (girder$, y_offset_tot$,                  'Y_OFFSET_TOT')
call init_attribute_name1 (girder$, z_offset_tot$,                  'Z_OFFSET_TOT')
call init_attribute_name1 (girder$, x_pitch_tot$,                   'X_PITCH_TOT')
call init_attribute_name1 (girder$, y_pitch_tot$,                   'Y_PITCH_TOT')
call init_attribute_name1 (girder$, origin_ele$,                    'ORIGIN_ELE')
call init_attribute_name1 (girder$, origin_ele_ref_pt$,             'ORIGIN_ELE_REF_PT')
call init_attribute_name1 (girder$, dx_origin$,                     'DX_ORIGIN')
call init_attribute_name1 (girder$, dy_origin$,                     'DY_ORIGIN')
call init_attribute_name1 (girder$, dz_origin$,                     'DZ_ORIGIN')
call init_attribute_name1 (girder$, dtheta_origin$,                 'DTHETA_ORIGIN')
call init_attribute_name1 (girder$, dphi_origin$,                   'DPHI_ORIGIN')
call init_attribute_name1 (girder$, dpsi_origin$,                   'DPSI_ORIGIN')
call init_attribute_name1 (girder$, is_on$,                         'IS_ON')

call init_attribute_name1 (overlay$, var$,                          'VAR')

call init_attribute_name1 (group$, var$,                            'VAR')
call init_attribute_name1 (group$, start_edge$,                     'START_EDGE')
call init_attribute_name1 (group$, end_edge$,                       'END_EDGE')
call init_attribute_name1 (group$, accordion_edge$,                 'ACCORDION_EDGE')
call init_attribute_name1 (group$, s_position$,                     'S_POSITION')

call init_attribute_name1 (lcavity$, field_autoscale$,              'FIELD_AUTOSCALE', quasi_free$)
call init_attribute_name1 (lcavity$, autoscale_amplitude$,          'AUTOSCALE_AMPLITUDE')
call init_attribute_name1 (lcavity$, autoscale_phase$,              'AUTOSCALE_PHASE')
call init_attribute_name1 (lcavity$, cavity_type$,                  'CAVITY_TYPE')
call init_attribute_name1 (lcavity$, p0c_start$,                    'P0C_START', dependent$)
call init_attribute_name1 (lcavity$, e_tot_start$,                  'E_TOT_START', dependent$)
call init_attribute_name1 (lcavity$, phi0_multipass$,               'PHI0_MULTIPASS')
call init_attribute_name1 (lcavity$, phi0$,                         'PHI0')
call init_attribute_name1 (lcavity$, gradient$,                     'GRADIENT')
call init_attribute_name1 (lcavity$, rf_frequency$,                 'RF_FREQUENCY')
call init_attribute_name1 (lcavity$, e_loss$,                       'E_LOSS')
call init_attribute_name1 (lcavity$, voltage$,                      'VOLTAGE', quasi_free$)
call init_attribute_name1 (lcavity$, field_calc$,                   'FIELD_CALC')
call init_attribute_name1 (lcavity$, field_master$,                 'FIELD_MASTER')
call init_attribute_name1 (lcavity$, lr_wake_file$,                 'LR_WAKE_FILE')
call init_attribute_name1 (lcavity$, lr_freq_spread$,               'LR_FREQ_SPREAD')
call init_attribute_name1 (lcavity$, coupler_strength$,             'COUPLER_STRENGTH')
call init_attribute_name1 (lcavity$, coupler_angle$,                'COUPLER_ANGLE')
call init_attribute_name1 (lcavity$, coupler_phase$,                'COUPLER_PHASE')
call init_attribute_name1 (lcavity$, coupler_at$,                   'COUPLER_AT')
call init_attribute_name1 (lcavity$, gradient_err$,                 'GRADIENT_ERR')
call init_attribute_name1 (lcavity$, voltage_err$,                  'VOLTAGE_ERR')
call init_attribute_name1 (lcavity$, phi0_err$,                     'PHI0_ERR')
call init_attribute_name1 (lcavity$, cartesian_map$,                'CARTESIAN_MAP')
call init_attribute_name1 (lcavity$, cylindrical_map$,              'CYLINDRICAL_MAP')
call init_attribute_name1 (lcavity$, grid_field$,                   'GRID_FIELD')
call init_attribute_name1 (lcavity$, taylor_field$,                 'TAYLOR_FIELD')
call init_attribute_name1 (lcavity$, phi0_autoscale$,               'PHI0_AUTOSCALE', quasi_free$)
call init_attribute_name1 (lcavity$, n_cell$,                       'N_CELL')

call init_attribute_name1 (marker$, l$,                             'L', dependent$)
call init_attribute_name1 (marker$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (marker$, p0c_start$,                     'p0c_start', private$)
call init_attribute_name1 (marker$, x_ray_line_len$,                'X_RAY_LINE_LEN')
call init_attribute_name1 (marker$, e_tot_ref_init$,                'e_tot_ref_init', private$)
call init_attribute_name1 (marker$, p0c_ref_init$,                  'p0c_ref_init', private$)

call init_attribute_name1 (match$, l$,                              'L')
call init_attribute_name1 (match$, beta_a0$,                        'BETA_A0')
call init_attribute_name1 (match$, alpha_a0$,                       'ALPHA_A0')
call init_attribute_name1 (match$, beta_b0$,                        'BETA_B0')
call init_attribute_name1 (match$, alpha_b0$,                       'ALPHA_B0')
call init_attribute_name1 (match$, beta_a1$,                        'BETA_A1')
call init_attribute_name1 (match$, alpha_a1$,                       'ALPHA_A1')
call init_attribute_name1 (match$, beta_b1$,                        'BETA_B1')
call init_attribute_name1 (match$, alpha_b1$,                       'ALPHA_B1')
call init_attribute_name1 (match$, dphi_a$,                         'DPHI_A')
call init_attribute_name1 (match$, dphi_b$,                         'DPHI_B')
call init_attribute_name1 (match$, eta_x0$,                         'ETA_X0')
call init_attribute_name1 (match$, etap_x0$,                        'ETAP_X0')
call init_attribute_name1 (match$, eta_y0$,                         'ETA_Y0')
call init_attribute_name1 (match$, etap_y0$,                        'ETAP_Y0')
call init_attribute_name1 (match$, eta_x1$,                         'ETA_X1')
call init_attribute_name1 (match$, etap_x1$,                        'ETAP_X1')
call init_attribute_name1 (match$, eta_y1$,                         'ETA_Y1')
call init_attribute_name1 (match$, etap_y1$,                        'ETAP_Y1')
call init_attribute_name1 (match$, match_end$,                      'MATCH_END')
call init_attribute_name1 (match$, x0$,                             'X0')
call init_attribute_name1 (match$, px0$,                            'PX0')
call init_attribute_name1 (match$, y0$,                             'Y0')
call init_attribute_name1 (match$, py0$,                            'PY0')
call init_attribute_name1 (match$, z0$,                             'Z0')
call init_attribute_name1 (match$, pz0$,                            'PZ0')
call init_attribute_name1 (match$, x1$,                             'X1')
call init_attribute_name1 (match$, px1$,                            'PX1')
call init_attribute_name1 (match$, y1$,                             'Y1')
call init_attribute_name1 (match$, py1$,                            'PY1')
call init_attribute_name1 (match$, z1$,                             'Z1')
call init_attribute_name1 (match$, pz1$,                            'PZ1')
call init_attribute_name1 (match$, match_end_orbit$,                'MATCH_END_ORBIT')
call init_attribute_name1 (match$, is_on$,                          'IS_ON')
call init_attribute_name1 (match$, E_tot_start$,                    'E_tot_start', private$)
call init_attribute_name1 (match$, p0c_start$,                      'p0c_start', private$)

call init_attribute_name1 (monitor$, field_calc$,                   'FIELD_CALC')
call init_attribute_name1 (monitor$, E_tot_start$,                  'E_tot_start', private$)
call init_attribute_name1 (monitor$, p0c_start$,                    'p0c_start', private$)

attrib_array(instrument$, :)                         = attrib_array(monitor$, :)
attrib_array(pipe$, :)                               = attrib_array(monitor$, :)

call init_attribute_name1 (hkicker$, kick$,                         'KICK', quasi_free$)
call init_attribute_name1 (hkicker$, field_calc$,                   'FIELD_CALC')
call init_attribute_name1 (hkicker$, field_master$,                 'FIELD_MASTER')
call init_attribute_name1 (hkicker$, bl_kick$,                      'BL_KICK', quasi_free$)
call init_attribute_name1 (hkicker$, E_tot_start$,                  'E_tot_start', private$)
call init_attribute_name1 (hkicker$, p0c_start$,                    'p0c_start', private$)
call init_attribute_name1 (hkicker$, cartesian_map$,                'CARTESIAN_MAP')
call init_attribute_name1 (hkicker$, cylindrical_map$,              'CYLINDRICAL_MAP')
call init_attribute_name1 (hkicker$, grid_field$,                   'GRID_FIELD')
call init_attribute_name1 (hkicker$, taylor_field$,                 'TAYLOR_FIELD')

attrib_array(vkicker$, :) = attrib_array(hkicker$, :)

call init_attribute_name1 (kicker$, h_displace$,                    'H_DISPLACE')
call init_attribute_name1 (kicker$, v_displace$,                    'V_DISPLACE')
call init_attribute_name1 (kicker$, r0_mag$,                        'R0_MAG')
call init_attribute_name1 (kicker$, r0_elec$,                       'R0_ELEC')
call init_attribute_name1 (kicker$, field_calc$,                    'FIELD_CALC')
call init_attribute_name1 (kicker$, field_master$,                  'FIELD_MASTER')
call init_attribute_name1 (kicker$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (kicker$, p0c_start$,                     'p0c_start', private$)
call init_attribute_name1 (kicker$, cartesian_map$,                 'CARTESIAN_MAP')
call init_attribute_name1 (kicker$, cylindrical_map$,               'CYLINDRICAL_MAP')
call init_attribute_name1 (kicker$, grid_field$,                    'GRID_FIELD')
call init_attribute_name1 (kicker$, taylor_field$,                  'TAYLOR_FIELD')

call init_attribute_name1 (bend_sol_quad$, angle$,                  'ANGLE')
call init_attribute_name1 (bend_sol_quad$, k1$,                     'K1')
call init_attribute_name1 (bend_sol_quad$, b1_gradient$,            'B1_GRADIENT', quasi_free$)
call init_attribute_name1 (bend_sol_quad$, g$,                      'G')
call init_attribute_name1 (bend_sol_quad$, b_field$,                'B_FIELD', quasi_free$)
call init_attribute_name1 (bend_sol_quad$, ks$,                     'KS')
call init_attribute_name1 (bend_sol_quad$, bs_field$,               'BS_FIELD', quasi_free$)
call init_attribute_name1 (bend_sol_quad$, dks_ds$,                 'DKS_DS')
call init_attribute_name1 (bend_sol_quad$, quad_tilt$,              'QUAD_TILT')
call init_attribute_name1 (bend_sol_quad$, bend_tilt$,              'BEND_TILT')
call init_attribute_name1 (bend_sol_quad$, x_quad$,                 'X_QUAD')
call init_attribute_name1 (bend_sol_quad$, y_quad$,                 'Y_QUAD')
call init_attribute_name1 (bend_sol_quad$, rho$,                    'RHO')
call init_attribute_name1 (bend_sol_quad$, r0_mag$,                 'R0_MAG')
call init_attribute_name1 (bend_sol_quad$, r0_elec$,                'R0_ELEC')
call init_attribute_name1 (bend_sol_quad$, field_calc$,             'FIELD_CALC')
call init_attribute_name1 (bend_sol_quad$, field_master$,           'FIELD_MASTER')
call init_attribute_name1 (bend_sol_quad$, cartesian_map$,          'CARTESIAN_MAP')
call init_attribute_name1 (bend_sol_quad$, cylindrical_map$,        'CYLINDRICAL_MAP')
call init_attribute_name1 (bend_sol_quad$, grid_field$,             'GRID_FIELD')
call init_attribute_name1 (bend_sol_quad$, taylor_field$,           'TAYLOR_FIELD')
call init_attribute_name1 (bend_sol_quad$, E_tot_start$,            'E_tot_start', private$)
call init_attribute_name1 (bend_sol_quad$, p0c_start$,              'p0c_start', private$)

call init_attribute_name1 (custom$, val1$,                          'VAL1')
call init_attribute_name1 (custom$, val2$,                          'VAL2')
call init_attribute_name1 (custom$, val3$,                          'VAL3')
call init_attribute_name1 (custom$, val4$,                          'VAL4')
call init_attribute_name1 (custom$, val5$,                          'VAL5')
call init_attribute_name1 (custom$, val6$,                          'VAL6')
call init_attribute_name1 (custom$, val7$,                          'VAL7')
call init_attribute_name1 (custom$, val8$,                          'VAL8')
call init_attribute_name1 (custom$, val9$,                          'VAL9')
call init_attribute_name1 (custom$, val10$,                         'VAL10')
call init_attribute_name1 (custom$, val11$,                         'VAL11')
call init_attribute_name1 (custom$, val12$,                         'VAL12')
call init_attribute_name1 (custom$, field_calc$,                    'FIELD_CALC')
call init_attribute_name1 (custom$, field_master$,                  'FIELD_MASTER')
call init_attribute_name1 (custom$, delta_e$,                       'DELTA_E')
call init_attribute_name1 (custom$, e_tot_start$,                   'E_TOT_START', dependent$)
call init_attribute_name1 (custom$, p0c_start$,                     'P0C_START', dependent$)

call init_attribute_name1 (floor_shift$, l$,                        'L')
call init_attribute_name1 (floor_shift$, origin_ele$,               'ORIGIN_ELE')
call init_attribute_name1 (floor_shift$, origin_ele_ref_pt$,        'ORIGIN_ELE_REF_PT')
call init_attribute_name1 (floor_shift$, upstream_ele_dir$,         'UPSTREAM_ELE_DIR', dependent$)
call init_attribute_name1 (floor_shift$, downstream_ele_dir$,       'DOWNSTREAM_ELE_DIR', dependent$)

call init_attribute_name1 (fiducial$, l$,                          'L', dependent$)
call init_attribute_name1 (fiducial$, origin_ele$,                 'ORIGIN_ELE')
call init_attribute_name1 (fiducial$, origin_ele_ref_pt$,          'ORIGIN_ELE_REF_PT')
call init_attribute_name1 (fiducial$, dx_origin$,                  'DX_ORIGIN')
call init_attribute_name1 (fiducial$, dy_origin$,                  'DY_ORIGIN')
call init_attribute_name1 (fiducial$, dz_origin$,                  'DZ_ORIGIN')
call init_attribute_name1 (fiducial$, dtheta_origin$,              'DTHETA_ORIGIN')
call init_attribute_name1 (fiducial$, dphi_origin$,                'DPHI_ORIGIN')
call init_attribute_name1 (fiducial$, dpsi_origin$,                'DPSI_ORIGIN')

call init_attribute_name1 (quadrupole$, k1$,                        'K1', quasi_free$)
call init_attribute_name1 (quadrupole$, B1_gradient$,               'B1_GRADIENT', quasi_free$)
call init_attribute_name1 (quadrupole$, r0_mag$,                    'R0_MAG')
call init_attribute_name1 (quadrupole$, r0_elec$,                   'R0_ELEC')
call init_attribute_name1 (quadrupole$, field_calc$,                'FIELD_CALC')
call init_attribute_name1 (quadrupole$, field_master$,              'FIELD_MASTER')
call init_attribute_name1 (quadrupole$, cartesian_map$,             'CARTESIAN_MAP')
call init_attribute_name1 (quadrupole$, cylindrical_map$,           'CYLINDRICAL_MAP')
call init_attribute_name1 (quadrupole$, grid_field$,                'GRID_FIELD')
call init_attribute_name1 (quadrupole$, taylor_field$,              'TAYLOR_FIELD')
call init_attribute_name1 (quadrupole$, E_tot_start$,               'E_tot_start', private$)
call init_attribute_name1 (quadrupole$, p0c_start$,                 'p0c_start', private$)
call init_attribute_name1 (quadrupole$, fq1$,                       'FQ1')
call init_attribute_name1 (quadrupole$, fq2$,                       'FQ2')

call init_attribute_name1 (sextupole$, k2$,                         'K2', quasi_free$)
call init_attribute_name1 (sextupole$, B2_gradient$,                'B2_GRADIENT', quasi_free$)
call init_attribute_name1 (sextupole$, r0_mag$,                     'R0_MAG')
call init_attribute_name1 (sextupole$, r0_elec$,                    'R0_ELEC')
call init_attribute_name1 (sextupole$, field_calc$,                 'FIELD_CALC')
call init_attribute_name1 (sextupole$, field_master$,               'FIELD_MASTER')
call init_attribute_name1 (sextupole$, cartesian_map$,              'CARTESIAN_MAP')
call init_attribute_name1 (sextupole$, cylindrical_map$,            'CYLINDRICAL_MAP')
call init_attribute_name1 (sextupole$, grid_field$,                 'GRID_FIELD')
call init_attribute_name1 (sextupole$, taylor_field$,               'TAYLOR_FIELD')
call init_attribute_name1 (sextupole$, E_tot_start$,                'E_tot_start', private$)
call init_attribute_name1 (sextupole$, p0c_start$,                  'p0c_start', private$)

call init_attribute_name1 (octupole$, k3$,                          'K3', quasi_free$)
call init_attribute_name1 (octupole$, B3_gradient$,                 'B3_GRADIENT', quasi_free$)
call init_attribute_name1 (octupole$, r0_mag$,                      'R0_MAG')
call init_attribute_name1 (octupole$, r0_elec$,                     'R0_ELEC')
call init_attribute_name1 (octupole$, field_calc$,                  'FIELD_CALC')
call init_attribute_name1 (octupole$, field_master$,                'FIELD_MASTER')
call init_attribute_name1 (octupole$, cartesian_map$,               'CARTESIAN_MAP')
call init_attribute_name1 (octupole$, cylindrical_map$,             'CYLINDRICAL_MAP')
call init_attribute_name1 (octupole$, grid_field$,                  'GRID_FIELD')
call init_attribute_name1 (octupole$, taylor_field$,                'TAYLOR_FIELD')
call init_attribute_name1 (octupole$, E_tot_start$,                 'E_tot_start', private$)
call init_attribute_name1 (octupole$, p0c_start$,                   'p0c_start', private$)

call init_attribute_name1 (patch$, l$,                              'L', dependent$)
call init_attribute_name1 (patch$, t_offset$,                       'T_OFFSET')
call init_attribute_name1 (patch$, p0c_start$,                      'P0C_START', dependent$)
call init_attribute_name1 (patch$, e_tot_start$,                    'E_TOT_START', dependent$)
call init_attribute_name1 (patch$, e_tot_offset$,                   'E_TOT_OFFSET')
call init_attribute_name1 (patch$, flexible$,                       'FLEXIBLE')
call init_attribute_name1 (patch$, field_calc$,                     'FIELD_CALC')
call init_attribute_name1 (patch$, upstream_ele_dir$,               'UPSTREAM_ELE_DIR', dependent$)
call init_attribute_name1 (patch$, downstream_ele_dir$,             'DOWNSTREAM_ELE_DIR', dependent$)
call init_attribute_name1 (patch$, ref_coordinates$,                'REF_COORDINATES')

call init_attribute_name1 (rfcavity$, field_autoscale$,             'FIELD_AUTOSCALE', quasi_free$)
call init_attribute_name1 (rfcavity$, autoscale_amplitude$,         'AUTOSCALE_AMPLITUDE')
call init_attribute_name1 (rfcavity$, autoscale_phase$,             'AUTOSCALE_PHASE')
call init_attribute_name1 (rfcavity$, cavity_type$,                 'CAVITY_TYPE')
call init_attribute_name1 (rfcavity$, voltage$,                     'VOLTAGE')
call init_attribute_name1 (rfcavity$, rf_frequency$,                'RF_FREQUENCY', quasi_free$)
call init_attribute_name1 (rfcavity$, phi0_multipass$,              'PHI0_MULTIPASS')
call init_attribute_name1 (rfcavity$, phi0$,                        'PHI0')
call init_attribute_name1 (rfcavity$, harmon$,                      'HARMON', quasi_free$)
call init_attribute_name1 (rfcavity$, field_calc$,                  'FIELD_CALC')
call init_attribute_name1 (rfcavity$, harmon_master$,               'HARMON_MASTER')
call init_attribute_name1 (rfcavity$, lr_wake_file$,                'LR_WAKE_FILE')
call init_attribute_name1 (rfcavity$, lr_freq_spread$,              'LR_FREQ_SPREAD')
call init_attribute_name1 (rfcavity$, coupler_strength$,            'COUPLER_STRENGTH')
call init_attribute_name1 (rfcavity$, coupler_angle$,               'COUPLER_ANGLE')
call init_attribute_name1 (rfcavity$, coupler_phase$,               'COUPLER_PHASE')
call init_attribute_name1 (rfcavity$, coupler_at$,                  'COUPLER_AT')
call init_attribute_name1 (rfcavity$, cartesian_map$,               'CARTESIAN_MAP')
call init_attribute_name1 (rfcavity$, cylindrical_map$,             'CYLINDRICAL_MAP')
call init_attribute_name1 (rfcavity$, grid_field$,                  'GRID_FIELD')
call init_attribute_name1 (rfcavity$, taylor_field$,                'TAYLOR_FIELD')
call init_attribute_name1 (rfcavity$, phi0_autoscale$,              'PHI0_AUTOSCALE', quasi_free$)
call init_attribute_name1 (rfcavity$, n_cell$,                      'N_CELL')
call init_attribute_name1 (rfcavity$, phi0_max$,                    'phi0_max', private$)
call init_attribute_name1 (rfcavity$, E_tot_start$,                 'E_tot_start', private$)
call init_attribute_name1 (rfcavity$, p0c_start$,                   'p0c_start', private$)
call init_attribute_name1 (rfcavity$, phi0_err$,                    'phi0_err', private$)
call init_attribute_name1 (rfcavity$, gradient$,                    'GRADIENT', dependent$)
call init_attribute_name1 (rfcavity$, gradient_err$,                'gradient_err', private$)
call init_attribute_name1 (rfcavity$, voltage_err$,                 'voltage_err', private$)

call init_attribute_name1 (sbend$, angle$,                          'ANGLE', quasi_free$)
call init_attribute_name1 (sbend$, ref_tilt$,                       'REF_TILT')
call init_attribute_name1 (sbend$, ref_tilt_tot$,                   'REF_TILT_TOT', dependent$)
call init_attribute_name1 (sbend$, e1$,                             'E1')
call init_attribute_name1 (sbend$, e2$,                             'E2')
call init_attribute_name1 (sbend$, h1$,                             'H1')
call init_attribute_name1 (sbend$, h2$,                             'H2')
call init_attribute_name1 (sbend$, k1$,                             'K1', quasi_free$)
call init_attribute_name1 (sbend$, k2$,                             'K2', quasi_free$)
call init_attribute_name1 (sbend$, g$,                              'G', quasi_free$)
call init_attribute_name1 (sbend$, g_err$,                          'G_ERR', quasi_free$)
call init_attribute_name1 (sbend$, roll$,                           'ROLL', override = .true.)
call init_attribute_name1 (sbend$, roll_tot$,                       'ROLL_TOT', dependent$, override = .true.)
call init_attribute_name1 (sbend$, hgap$,                           'HGAP')
call init_attribute_name1 (sbend$, hgapx$,                          'HGAPX')
call init_attribute_name1 (sbend$, fint$,                           'FINT')
call init_attribute_name1 (sbend$, fintx$,                          'FINTX')
call init_attribute_name1 (sbend$, rho$,                            'RHO', quasi_free$)
call init_attribute_name1 (sbend$, l_chord$,                        'L_CHORD', quasi_free$)
call init_attribute_name1 (sbend$, ptc_fringe_geometry$,            'PTC_FRINGE_GEOMETRY')
call init_attribute_name1 (sbend$, higher_order_fringe_type$,       'HIGHER_ORDER_FRINGE_TYPE')
call init_attribute_name1 (sbend$, b_field$,                        'B_FIELD', quasi_free$)
call init_attribute_name1 (sbend$, b_field_err$,                    'B_FIELD_ERR', quasi_free$)
call init_attribute_name1 (sbend$, b1_gradient$,                    'B1_GRADIENT', quasi_free$)
call init_attribute_name1 (sbend$, b2_gradient$,                    'B2_GRADIENT', quasi_free$)
call init_attribute_name1 (sbend$, r0_mag$,                         'R0_MAG')
call init_attribute_name1 (sbend$, r0_elec$,                        'R0_ELEC')
call init_attribute_name1 (sbend$, field_calc$,                     'FIELD_CALC')
call init_attribute_name1 (sbend$, field_master$,                   'FIELD_MASTER')
call init_attribute_name1 (sbend$, E_tot_start$,                    'E_tot_start', private$)
call init_attribute_name1 (sbend$, p0c_start$,                      'p0c_start', private$)
call init_attribute_name1 (sbend$, ptc_field_geometry$,             'PTC_FIELD_GEOMETRY')
call init_attribute_name1 (sbend$, cartesian_map$,                  'CARTESIAN_MAP')
call init_attribute_name1 (sbend$, cylindrical_map$,                'CYLINDRICAL_MAP')
call init_attribute_name1 (sbend$, grid_field$,                     'GRID_FIELD')
call init_attribute_name1 (sbend$, taylor_field$,                   'TAYLOR_FIELD')

attrib_array(rbend$, :) = attrib_array(sbend$, :)

call init_attribute_name1 (solenoid$, ks$,                          'KS', quasi_free$)
call init_attribute_name1 (solenoid$, bs_field$,                    'BS_FIELD', quasi_free$)
call init_attribute_name1 (solenoid$, r0_mag$,                      'R0_MAG')
call init_attribute_name1 (solenoid$, r0_elec$,                     'R0_ELEC')
call init_attribute_name1 (solenoid$, field_calc$,                  'FIELD_CALC')
call init_attribute_name1 (solenoid$, field_master$,                'FIELD_MASTER')
call init_attribute_name1 (solenoid$, cartesian_map$,               'CARTESIAN_MAP')
call init_attribute_name1 (solenoid$, cylindrical_map$,             'CYLINDRICAL_MAP')
call init_attribute_name1 (solenoid$, grid_field$,                  'GRID_FIELD')
call init_attribute_name1 (solenoid$, taylor_field$,                'TAYLOR_FIELD')
call init_attribute_name1 (solenoid$, E_tot_start$,                 'E_tot_start', private$)
call init_attribute_name1 (solenoid$, p0c_start$,                   'p0c_start', private$)

attrib_array(undulator$, :) = attrib_array(wiggler$, :)

call init_attribute_name1 (sample$, l$,                             'L')
call init_attribute_name1 (sample$, mode$,                          'MODE')
call init_attribute_name1 (sample$, material_type$,                 'MATERIAL_TYPE')

call init_attribute_name1 (sol_quad$, k1$,                          'K1', quasi_free$)
call init_attribute_name1 (sol_quad$, ks$,                          'KS', quasi_free$)
call init_attribute_name1 (sol_quad$, b1_gradient$,                 'B1_GRADIENT', quasi_free$)
call init_attribute_name1 (sol_quad$, bs_field$,                    'BS_FIELD', quasi_free$)
call init_attribute_name1 (sol_quad$, r0_mag$,                      'R0_MAG')
call init_attribute_name1 (sol_quad$, r0_elec$,                     'R0_ELEC') 
call init_attribute_name1 (sol_quad$, field_calc$,                  'FIELD_CALC')
call init_attribute_name1 (sol_quad$, field_master$,                'FIELD_MASTER')
call init_attribute_name1 (sol_quad$, cartesian_map$,               'CARTESIAN_MAP')
call init_attribute_name1 (sol_quad$, cylindrical_map$,             'CYLINDRICAL_MAP')
call init_attribute_name1 (sol_quad$, grid_field$,                  'GRID_FIELD')
call init_attribute_name1 (sol_quad$, taylor_field$,                'TAYLOR_FIELD')
call init_attribute_name1 (sol_quad$, E_tot_start$,                 'E_tot_start', private$)
call init_attribute_name1 (sol_quad$, p0c_start$,                   'p0c_start', private$)

attrib_array(multipole$, k0l$:k21l$)%name    = &
             ['K0L ', 'K1L ', 'K2L ', 'K3L ', 'K4L ', 'K5L ', 'K6L ', 'K7L ', 'K8L ', 'K9L ', 'K10L', &
                      'K11L', 'K12L', 'K13L', 'K14L', 'K15L', 'K16L', 'K17L', 'K18L', 'K19L', 'K20L', 'K21L']
attrib_array(multipole$, t0$:t21$)%name = ['T0 ', &
                               'T1 ', 'T2 ', 'T3 ', 'T4 ', 'T5 ', & 
                               'T6 ', 'T7 ', 'T8 ', 'T9 ', 'T10', &
                               'T11', 'T12', 'T13', 'T14', 'T15', &
                               'T16', 'T17', 'T18', 'T19', 'T20', 'T21']
attrib_array(multipole$, k0l$:t21$)%type = is_free$
call init_attribute_name1 (multipole$, l$,                          'L')
call init_attribute_name1 (multipole$, x_pitch$,          null_name$, does_not_exist$, .true.)
call init_attribute_name1 (multipole$, y_pitch$,          null_name$, does_not_exist$, .true.)
call init_attribute_name1 (multipole$, x_pitch_tot$,      null_name$, does_not_exist$, .true.)
call init_attribute_name1 (multipole$, y_pitch_tot$,      null_name$, does_not_exist$, .true.)

call init_attribute_name1 (ab_multipole$, l$,                       'L')
call init_attribute_name1 (ab_multipole$, x_pitch$,       null_name$, does_not_exist$, .true.)
call init_attribute_name1 (ab_multipole$, y_pitch$,       null_name$, does_not_exist$, .true.)
call init_attribute_name1 (ab_multipole$, x_pitch_tot$,   null_name$, does_not_exist$, .true.)
call init_attribute_name1 (ab_multipole$, y_pitch_tot$,   null_name$, does_not_exist$, .true.)

call init_attribute_name1 (sad_mult$, num_steps$,              'NUM_STEPS', dependent$)
call init_attribute_name1 (sad_mult$, ds_step$,                'DS_STEP', dependent$)
call init_attribute_name1 (sad_mult$, eps_step_scale$,         'EPS_STEP_SCALE')
call init_attribute_name1 (sad_mult$, angle$,                  'ANGLE')
call init_attribute_name1 (sad_mult$, e1$,                     'E1')
call init_attribute_name1 (sad_mult$, e2$,                     'E2')
!call init_attribute_name1 (sad_mult$, rf_frequency$,           'RF_FREQUENCY')  ! SAD: freq
!call init_attribute_name1 (sad_mult$, phi0$,                   'PHI0')         ! SAD: phi
!call init_attribute_name1 (sad_mult$, phi0_err$,               'PHI0_ERR')     ! SAD: dphi
!call init_attribute_name1 (sad_mult$, voltage$,                'VOLTAGE')      ! SAD: volt
!call init_attribute_name1 (sad_mult$, harmon$,                 'HARMON')       ! SAD: harm
call init_attribute_name1 (sad_mult$, fringe_at$,              'FRINGE_AT')    ! SAD: fringe
call init_attribute_name1 (sad_mult$, fringe_type$,            'FRINGE_TYPE')  ! SAD: disfrin
call init_attribute_name1 (sad_mult$, spin_fringe_on$,         'SPIN_FRINGE_ON')
call init_attribute_name1 (sad_mult$, fq1$,                    'FQ1')
call init_attribute_name1 (sad_mult$, fq2$,                    'FQ2')
call init_attribute_name1 (sad_mult$, fb1$,                    'FB1')
call init_attribute_name1 (sad_mult$, fb2$,                    'FB2')
call init_attribute_name1 (sad_mult$, bs_field$,               'BS_FIELD')
call init_attribute_name1 (sad_mult$, x_offset_mult$,          'X_OFFSET_MULT')
call init_attribute_name1 (sad_mult$, y_offset_mult$,          'Y_OFFSET_MULT')
call init_attribute_name1 (sad_mult$, x_pitch_mult$,           'X_PITCH_MULT')
call init_attribute_name1 (sad_mult$, y_pitch_mult$,           'Y_PITCH_MULT')
! sad_mult Attributes with no SAD equivalent
call init_attribute_name1 (sad_mult$, rho$,                    'RHO')   
call init_attribute_name1 (sad_mult$, g$,                      'G')
call init_attribute_name1 (sad_mult$, ks$,                     'KS')
call init_attribute_name1 (sad_mult$, b_field$,                'B_FIELD')

call init_attribute_name1 (hybrid$, l$,                             'L')
call init_attribute_name1 (hybrid$, delta_e$,                       'DELTA_E')
call init_attribute_name1 (hybrid$, e_tot_start$,                   'E_TOT_START', dependent$)
call init_attribute_name1 (hybrid$, p0c_start$,                     'P0C_START', dependent$)
call init_attribute_name1 (hybrid$, delta_ref_time$,                'DELTA_REF_TIME', override = .true.)

call init_attribute_name1 (mirror$, l$,                             'L', dependent$)
call init_attribute_name1 (mirror$, graze_angle$,                   'GRAZE_ANGLE')
call init_attribute_name1 (mirror$, critical_angle$,                'CRITICAL_ANGLE')
call init_attribute_name1 (mirror$, ref_tilt$,                      'REF_TILT')
call init_attribute_name1 (mirror$, ref_tilt_tot$,                  'REF_TILT_TOT', dependent$)
call init_attribute_name1 (mirror$, ref_wavelength$,                'REF_WAVELENGTH', dependent$)
call init_attribute_name1 (mirror$, diffraction_limited$,           'DIFFRACTION_LIMITED')

call init_attribute_name1 (multilayer_mirror$, l$,                    'L', dependent$)
call init_attribute_name1 (multilayer_mirror$, graze_angle$,          'GRAZE_ANGLE')
call init_attribute_name1 (multilayer_mirror$, ref_tilt$,             'REF_TILT')
call init_attribute_name1 (multilayer_mirror$, ref_tilt_tot$,         'REF_TILT_TOT', dependent$)
call init_attribute_name1 (multilayer_mirror$, n_cell$,               'N_CELL')
call init_attribute_name1 (multilayer_mirror$, d1_thickness$,         'D1_THICKNESS')
call init_attribute_name1 (multilayer_mirror$, d2_thickness$,         'D2_THICKNESS')
call init_attribute_name1 (multilayer_mirror$, v1_unitcell$,          'V1_UNITCELL')
call init_attribute_name1 (multilayer_mirror$, v2_unitcell$,          'V2_UNITCELL')
call init_attribute_name1 (multilayer_mirror$, ref_wavelength$,       'REF_WAVELENGTH', dependent$)
call init_attribute_name1 (multilayer_mirror$, material_type$,        'MATERIAL_TYPE')
call init_attribute_name1 (multilayer_mirror$, e_tot_start$,          'e_tot_start', private$)
call init_attribute_name1 (multilayer_mirror$, p0c_start$,            'p0c_start', private$)
call init_attribute_name1 (multilayer_mirror$, diffraction_limited$,  'DIFFRACTION_LIMITED')

call init_attribute_name1 (taylor$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (taylor$, p0c_start$,                     'p0c_start', private$)
call init_attribute_name1 (taylor$, ref_orbit$,                     'REF_ORBIT')
call init_attribute_name1 (taylor$, tt$,                            'TT*')

call init_attribute_name1 (wiggler$, k1$,                           'K1', dependent$)
call init_attribute_name1 (wiggler$, l_pole$,                       'L_POLE')
call init_attribute_name1 (wiggler$, b_max$,                        'B_MAX')
call init_attribute_name1 (wiggler$, rho$,                          'RHO', dependent$)
call init_attribute_name1 (wiggler$, n_pole$,                       'N_POLE')
call init_attribute_name1 (wiggler$, term$,                         'TERM')
call init_attribute_name1 (wiggler$, polarity$,                     'POLARITY')
call init_attribute_name1 (wiggler$, r0_mag$,                       'R0_MAG')
call init_attribute_name1 (wiggler$, r0_elec$,                      'R0_ELEC')
call init_attribute_name1 (wiggler$, field_calc$,                   'FIELD_CALC')
call init_attribute_name1 (wiggler$, field_master$,                 'FIELD_MASTER')
call init_attribute_name1 (wiggler$, x_ray_line_len$,               'X_RAY_LINE_LEN')
call init_attribute_name1 (wiggler$, cartesian_map$,                'CARTESIAN_MAP')
call init_attribute_name1 (wiggler$, cylindrical_map$,              'CYLINDRICAL_MAP')
call init_attribute_name1 (wiggler$, grid_field$,                   'GRID_FIELD')
call init_attribute_name1 (wiggler$, taylor_field$,                 'TAYLOR_FIELD')
call init_attribute_name1 (wiggler$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (wiggler$, p0c_start$,                     'p0c_start', private$)

attrib_array(undulator$, :) = attrib_array(wiggler$, :)

call init_attribute_name1 (photon_init$, l$,                         'L', dependent$)
call init_attribute_name1 (photon_init$, sig_x$,                     'SIG_X')
call init_attribute_name1 (photon_init$, sig_y$,                     'SIG_Y')
call init_attribute_name1 (photon_init$, sig_z$,                     'SIG_Z')
call init_attribute_name1 (photon_init$, sig_vx$,                    'SIG_VX')
call init_attribute_name1 (photon_init$, sig_vy$,                    'SIG_VY')
call init_attribute_name1 (photon_init$, sig_E$,                     'SIG_E')
call init_attribute_name1 (photon_init$, E_center$,                  'E_CENTER')
call init_attribute_name1 (photon_init$, E_center_relative_to_ref$,  'E_CENTER_RELATIVE_TO_REF')
call init_attribute_name1 (photon_init$, spatial_distribution$,      'SPATIAL_DISTRIBUTION')
call init_attribute_name1 (photon_init$, velocity_distribution$,     'VELOCITY_DISTRIBUTION')
call init_attribute_name1 (photon_init$, energy_distribution$,       'ENERGY_DISTRIBUTION')
call init_attribute_name1 (photon_init$, e_field_x$,                 'E_FIELD_X')
call init_attribute_name1 (photon_init$, e_field_y$,                 'E_FIELD_Y')
call init_attribute_name1 (photon_init$, scale_field_to_one$,        'SCALE_FIELD_TO_ONE')
call init_attribute_name1 (photon_init$, transverse_sigma_cut$,      'TRANSVERSE_SIGMA_CUT')
call init_attribute_name1 (photon_init$, ds_slice$,                  'DS_SLICE')
call init_attribute_name1 (photon_init$, physical_source$,           'PHYSICAL_SOURCE')
call init_attribute_name1 (photon_init$, ref_wavelength$,            'REF_WAVELENGTH', dependent$)

!-----------------------------------------------------------------------
! We make a short list to compare against to make things go faster.
! Note: A patch element has a z_offset, not an z_offset.

has_hkick_attributes = .false.  ! Defined in bmad_struct.f90
has_kick_attributes  = .false.  ! Defined in bmad_struct.f90
has_orientation_attributes_key = .false.  ! Defined in bmad_struct.f90

do i = 1, n_key$
  if (attrib_array(i, x_offset$)%name == 'X_OFFSET') has_orientation_attributes_key(i) = .true.
  if (attrib_array(i, kick$)%name     == 'KICK')  has_kick_attributes(i) = .true.
  if (attrib_array(i, hkick$)%name    == 'HKICK') has_hkick_attributes(i) = .true.

  call init_short_attrib_array(i)
enddo

attribute_array_init_needed = .false.

end subroutine init_attribute_name_array

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_short_attrib_array (ix_key)
!
! Internal routine to init the short_attrib_array array.
!-

subroutine init_short_attrib_array (ix_key)

integer ix_key, num, j

!

num = 0
do j = 1, num_ele_attrib_extended$
  if (attrib_array(ix_key, j)%name == null_name$) cycle
  if (attrib_array(ix_key, j)%type == private$) cycle
  num = num + 1
  short_attrib_array(ix_key, num) = attrib_array(ix_key, j)%name
  attrib_ix(ix_key, num) = j
enddo
attrib_num(ix_key) = num

end subroutine init_short_attrib_array 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_attribute_name1 (ix_key, ix_attrib, name, attrib_state, override)
!
! Routine to initialize a single name in the element attribute name table.
!
! Module needed:
!   use basic_attribute_mod
!
! Input:
!   ix_key       -- Integer: Key index.
!   ix_attrib    -- Integer: Attribute index.
!   name         -- Character(*): Attribute name.
!   attrib_state -- Integer, optional: Class of attribute: does_not_exist$, is_free$, etc.
!                     Defaults to is_free$.
!   override     -- Logical, optional: Normally this routine throws an error if 
!                     the [ix_key, ix_attrib] has been set previously. 
!                     If override = True then the set is done and no error is generated.
!-

subroutine init_attribute_name1 (ix_key, ix_attrib, name, attrib_state, override)

integer ix_key, ix_attrib
character(*) name
integer, optional :: attrib_state
logical, optional :: override

! Check that attrib_array(ix_key, ix_attrib)%name has not already been set.
! If so bomb program.

if (.not. logic_option(.false., override) .and. attrib_array(ix_key, ix_attrib)%name /= null_name$) then
  call out_io (s_fatal$, 'init_attribute_name1', 'ERROR IN INITIALIZING ATTRIB_ARRAY FOR: ' // key_name(ix_key), &
                  'IX_ATTRIB \i0\ ALREADY SET!', &
                  'OLD/NEW NAMES: ' // trim(attrib_array(ix_key, ix_attrib)%name) // ' : ' // name, &
                  i_array = [ix_attrib])
  if (global_com%exit_on_error) call err_exit
endif

attrib_array(ix_key, ix_attrib)%name    = name
attrib_array(ix_key, ix_attrib)%type = integer_option(is_free$, attrib_state)

! If things are done after the attribute array has been inited then the short table has
! to be reinited.

if (.not. attribute_array_init_needed) call init_short_attrib_array(ix_key)

end subroutine init_attribute_name1 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function has_orientation_attributes (ele) result (has_attribs)
!
! Routine to determine whether an element has orientation attributes like x_offset, etc.
!
! Module needed:
!   use bmad
!
! Input:
!   ele  -- ele_struct: Lattice element.
!
! Output:
!   has_attribs -- Logical: True if ele has orientation attributes. False otherwise.
!-

function has_orientation_attributes (ele) result (has_attribs)

type (ele_struct) ele
logical has_attribs

! Rule: The orientation attributes of an em_field control any fields that the element contains internally.
! If it is a super_slave or a mutipass_slave, the fields always reside in the lord elements and thus
! the em_field is not considered to have orientation attributes in this case.

has_attribs = has_orientation_attributes_key(ele%key)
if (ele%key == em_field$ .and. (ele%slave_status == super_slave$ .or. ele%slave_status == multipass_slave$)) &
      has_attribs = .false.

end function has_orientation_attributes

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_type (attrib_name) result (attrib_type)
!
! Routine to return the type of an attribute.
!
! A "switch" attribute is an attribute whose value corresponds to some string.
! For example, the "COUPLER_AT" attirbute with value 1 corresponds to "ENTRANCE_END", etc. 
!
! If attrib_name corresponds to a switch attribute, The routine switch_attrib_value_name can 
! be used to print the name corresponding to the attribute's value.
!
! Note: The "storage type" of an attribute is different from the "logical type" returned by
! this routine. For example, the logical type of attribute "n_slice" is integer. However, the 
! value of "n_slice" is stored as a real number in the ele_struct [in ele%value(n_slice$)]. 
!
! Modules needed:
!   use bmad
!
! Input:
!   attrib_name -- Character(*): Name of the attribute. Must be upper case.
!
! Output:
!   attrib_type  -- Integer: Attribute type: 
!                     is_string$, is_logical$, is_integer$, is_real$, is_switch$, or unknown$
!-

function attribute_type (attrib_name) result (attrib_type)

character(*) attrib_name
integer attrib_type

!

select case (attrib_name)
case ('MATCH_END', 'MATCH_END_ORBIT', 'NO_END_MARKER', 'SYMPLECTIFY', 'IS_ON', &
      'APERTURE_LIMIT_ON', 'ABSOLUTE_TIME_TRACKING', 'AUTOSCALE_PHASE', &
      'AUTOSCALE_AMPLITUDE', 'CSR_CALC_ON', 'PTC_EXACT_MODEL', 'PTC_EXACT_MISALIGN', &
      'TAYLOR_MAP_INCLUDES_OFFSETS', 'OFFSET_MOVES_APERTURE', 'FIELD_MASTER', 'SCALE_MULTIPOLES', &
      'FLEXIBLE', 'USE_HARD_EDGE_DRIFTS', 'NEW_BRANCH', 'HARMON_MASTER', 'SPIN_FRINGE_ON', &
      'BRANCHES_ARE_COHERENT', 'E_CENTER_RELATIVE_TO_REF', 'SCALE_FIELD_TO_ONE', 'DIFFRACTION_LIMITED', &
      'MULTIPOLES_ON')
  attrib_type = is_logical$

case ('TAYLOR_ORDER', 'N_SLICE', 'N_REF_PASS', 'DIRECTION', 'N_CELL', 'AD_N_DIV_MAX', &
      'IX_TO_BRANCH', 'IX_TO_ELEMENT', 'NUM_STEPS', 'INTEGRATOR_ORDER', &
      'PTC_MAX_FRINGE_ORDER', 'UPSTREAM_ELE_DIR', 'DOWNSTREAM_ELE_DIR')
  attrib_type = is_integer$

case ('APERTURE_AT', 'APERTURE_TYPE', 'COUPLER_AT', 'FIELD_CALC', &
      'FRINGE_TYPE', 'GEOMETRY', 'FRINGE_AT', 'MAT6_CALC_METHOD', 'HIGHER_ORDER_FRINGE_TYPE', &
      'ORIGIN_ELE_REF_PT', 'PARTICLE', 'PTC_FIELD_GEOMETRY', 'DEFAULT_TRACKING_SPECIES', &
      'PTC_INTEGRATION_TYPE', 'SPIN_TRACKING_METHOD', 'PTC_FRINGE_GEOMETRY', &
      'TRACKING_METHOD', 'REF_ORBIT_FOLLOWS', 'REF_COORDINATES', 'MODE', 'CAVITY_TYPE', &
      'SPATIAL_DISTRIBUTION', 'ENERGY_DISTRIBUTION', 'VELOCITY_DISTRIBUTION', 'KEY', 'SLAVE_STATUS', &
      'LORD_STATUS')
  attrib_type = is_switch$

case ('TYPE', 'ALIAS', 'DESCRIP', 'SR_WAKE_FILE', 'LR_WAKE_FILE', 'LATTICE', 'PHYSICAL_SOURCE', &
     'CRYSTAL_TYPE', 'MATERIAL_TYPE', 'REFERENCE', 'TO_LINE', 'TO_ELEMENT', 'ORIGIN_ELE', 'NAME')
  attrib_type = is_string$

case default
  if (attribute_index(0, attrib_name, can_abbreviate = .false.) == 0) then
    attrib_type = unknown$
  else
    attrib_type = is_real$
  endif
end select

end function attribute_type 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function corresponding_tot_attribute_index (ele, ix_attrib) result (ix_tot_attrib)
!
! Function to return the index of the corresponding "_tot" attribute.
! Corresponding attributes are:
!   Attribute         "_tot" attribute
!   ---------         -----------------
!   x_offset          x_offset_tot
!   x_pitch           x_pitch_tot
!   tilt              tilt_tot
!
! Module needed:
!   use bmad
!
! Input:
!   ele       -- Ele_struct: Element
!   ix_attrib -- Integer: Index of attribute
!
! Output:
!   ix_tot_attrib -- Integer: Index of corresponding _tot attribute.
!                      Set to -1 if no corresponding attribute.
!-

function corresponding_tot_attribute_index (ele, ix_attrib) result (ix_tot_attrib)

type (ele_struct) ele
integer ix_attrib, ix_tot_attrib
character(40) a_name

!

select case (attribute_name(ele, ix_attrib))
case ('X_PITCH');     ix_tot_attrib = x_pitch_tot$
case ('Y_PITCH');     ix_tot_attrib = y_pitch_tot$
case ('X_OFFSET');    ix_tot_attrib = x_offset_tot$
case ('Y_OFFSET');    ix_tot_attrib = y_offset_tot$
case ('Z_OFFSET');    ix_tot_attrib = z_offset_tot$
case ('REF_TILT');    ix_tot_attrib = ref_tilt_tot$
case ('TILT');        ix_tot_attrib = tilt_tot$
case ('ROLL');        ix_tot_attrib = roll_tot$
case default;         ix_tot_attrib = -1; return
end select

a_name = attribute_name(ele, ix_tot_attrib)
if (a_name(1:1) == '!') ix_tot_attrib = -1

end function corresponding_tot_attribute_index 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function is_a_tot_attribute (ele, ix_attrib) result (is_a_tot_attrib)
!
! Function returns True if ix_attrib corresponds to a "_tot" attribute.
! Example _tot attributes are:
!   -----------------
!   x_offset_tot
!   x_pitch_tot
!   tilt_tot
!
! Module needed:
!   use bmad
!
! Input:
!   ele       -- Ele_struct: Element
!   ix_attrib -- Integer: Index of attribute
!
! Output:
!   is_a_tot_attrib -- Logical: True if a _tot attribute. False otherwise.
!   
!-

function is_a_tot_attribute (ele, ix_attrib) result (is_a_tot_attrib)

type (ele_struct) ele
integer ix_attrib
logical is_a_tot_attrib

!

select case (attribute_name(ele, ix_attrib))
case ('X_PITCH_TOT', 'Y_PITCH_TOT', 'X_OFFSET_TOT', 'Y_OFFSET_TOT', &
      'REF_TILT_TOT', 'ROLL_TOT', 'Z_OFFSET_TOT', 'TILT_TOT')
  is_a_tot_attrib = .true.
case default
  is_a_tot_attrib = .false.
end select

end function is_a_tot_attribute

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine string_attrib (attrib_name, ele, attrib_value)
!
! Routine to return the value of a string attribute of a lattice element.
! This routine is useful when attrib_name is specified by the program user.
!
! For example:
!   call string_attrib ('NAME', ele, attrib_value)  ! Will return attrib_value = ele%name
!
! Input:
!   attrib_name   -- character(*): Name of the type of element attribute.
!   ele           -- ele_struct: Lattice element.
!
! Output:
!   attrib_value  -- character(*): The string associated with the attribute.
!-

subroutine string_attrib (attrib_name, ele, attrib_value)

type (ele_struct) ele
character(*) attrib_name, attrib_value
integer ib, ie

!

attrib_value = ''

select case (attrib_name)
case ('TYPE')
  attrib_value = ele%type
case ('ALIAS')
  attrib_value = ele%alias
case ('DESCRIP')
  if (associated(ele%descrip)) attrib_value = ele%descrip
case ('SR_WAKE_FILE')
  if (associated(ele%wake)) attrib_value = ele%wake%sr_file
case ('LR_WAKE_FILE')
  if (associated(ele%wake)) attrib_value = ele%wake%lr_file
case ('PHYSICAL_SOURCE')
  if (attribute_index(ele, attrib_name) /= 0) attrib_value = ele%component_name
case ('CRYSTAL_TYPE')
  if (attribute_index(ele, attrib_name) /= 0) attrib_value = ele%component_name
case ('MATERIAL_TYPE')
  if (attribute_index(ele, attrib_name) /= 0) attrib_value = ele%component_name
case ('TO_LINE')
  if (associated(ele%branch)) then
    ib = nint(ele%value(ix_to_branch$))
    attrib_value = ele%branch%lat%branch(ib)%name
  endif
case ('TO_ELEMENT')
  if (associated(ele%branch)) then
    ib = nint(ele%value(ix_to_branch$))
    ie = nint(ele%value(ix_to_element$))
    attrib_value = ele%branch%lat%branch(ib)%ele(ie)%name
  endif
case ('ORIGIN_ELE')
  if (attribute_index(ele, attrib_name) /= 0) attrib_value = ele%component_name
case ('NAME')
  attrib_value = ele%name
end select


end subroutine string_attrib

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function switch_attrib_value_name (attrib_name, attrib_value, ele, 
!                                     is_default) result (attrib_val_name)
!
! Routine to return the name corresponding to the value of a given switch attribute.
!
! This routine is for "switch" attributes. For example, the "aperture_type" attribute
! can have value names of "Entrance_End", "Exit_End", etc.
!
! Optionally, this routine can determine if the attribute value corresponds 
! to the default value. That is, the value that the attribute would have if 
! not specified in the lattice file.
!
! Use the routine attribute_type to first test if the type of the attribute
! corresponds to is_switch$. 
!
! Modules needed:
!   use bmad
!
! Input:
!   attrib_name  -- Character(*): Name of the attribute. Must be upper case.
!   attrib_value -- Real(rp): Value of the attribute.
!   ele          -- ele_struct" Lattice element that the attribute is contained in.
!
! Output:
!   attrib_val_name -- Character(40): Name corresponding to the value.
!   is_default      -- Logical, optional: If True then the value of the attiribute
!                        corresponds to the default value. If this argument is
!                        present, the ele argument must also be present.
!-

function switch_attrib_value_name (attrib_name, attrib_value, ele, is_default) result (attrib_val_name)

type (ele_struct) :: ele, ele2
character(*) attrib_name
real(rp) attrib_value
integer ix_attrib_val
character(40) attrib_val_name
character(24) :: r_name = 'switch_attrib_value_name'
logical, optional :: is_default

! 

ix_attrib_val = nint(attrib_value)

select case (attrib_name)

case ('APERTURE_AT')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, aperture_at_name, lbound(aperture_type_name, 1))
  if (present(is_default)) then
    is_default = (ix_attrib_val == exit_end$)
  endif

case ('APERTURE_TYPE')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, aperture_type_name, lbound(aperture_type_name, 1))
  if (present(is_default)) then
    if (ele%key == ecollimator$) then
      is_default = (ix_attrib_val == elliptical$)
    else
      is_default = (ix_attrib_val == rectangular$)
    endif
  endif

case ('CAVITY_TYPE')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, cavity_type_name, lbound(cavity_type_name, 1))
  if (present(is_default)) then
    is_default = (ix_attrib_val == default_value(ele, cavity_type$))
  endif

case ('COUPLER_AT')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, end_at_name, lbound(end_at_name, 1))
  if (present(is_default)) is_default = (ix_attrib_val == downstream_end$)

case ('DEFAULT_TRACKING_SPECIES')
  attrib_val_name = species_name(ix_attrib_val)
  if (present(is_default)) is_default = (ix_attrib_val == ref_particle$)

case ('ENERGY_DISTRIBUTION')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, distribution_name, lbound(distribution_name, 1))
  if (present(is_default)) is_default = (ix_attrib_val == gaussian$)

case ('PHOTON_TYPE')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, photon_type_name, lbound(photon_type_name, 1))
  if (present(is_default)) is_default = (ix_attrib_val == incoherent$)

case ('FIELD_CALC')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, field_calc_name, lbound(field_calc_name, 1))
  if (present(is_default)) is_default = (ix_attrib_val == bmad_standard$)

case ('FRINGE_TYPE')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, fringe_type_name, lbound(fringe_type_name, 1))
  if (present(is_default)) then
    select case (ele%key)
    case (sad_mult$)
      is_default = (ix_attrib_val == hard_edge_only$)      
    case (rbend$, sbend$)
      is_default = (ix_attrib_val == basic_bend$)
    case default
      is_default = (ix_attrib_val == none$)
    end select
  endif

case ('PTC_FRINGE_GEOMETRY')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, ptc_fringe_geometry_name, lbound(ptc_fringe_geometry_name, 1))
  if (present(is_default)) is_default = (ix_attrib_val == x_invariant$)

case ('GEOMETRY')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, geometry_name, lbound(geometry_name, 1))

case ('KEY')
    call get_this_attrib_name (attrib_val_name, ix_attrib_val, key_name, lbound(key_name, 1))
  if (present(is_default)) is_default = .false.

case ('LORD_STATUS')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, control_name, lbound(control_name, 1))  
  if (present(is_default)) is_default = .false.
  
case ('MODE')
  if (ele%key == diffraction_plate$ .or. ele%key == sample$ .or. ele%key == mask$) then
    call get_this_attrib_name (attrib_val_name, ix_attrib_val, mode_name, lbound(mode_name, 1))
    if (present(is_default)) then
      is_default = (ix_attrib_val == default_value(ele, mode$))
    endif
  else
    call get_this_attrib_name (attrib_val_name, ix_attrib_val, geometry_name, lbound(geometry_name, 1))
    if (present(is_default)) then
      is_default = (ix_attrib_val == open$)
    endif
  endif

case ('FRINGE_AT')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, end_at_name, lbound(end_at_name, 1))
  if (present(is_default)) then
    is_default = (ix_attrib_val == both_ends$)
  endif

case ('HIGHER_ORDER_FRINGE_TYPE')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, higher_order_fringe_type_name, lbound(higher_order_fringe_type_name, 1))
  if (present(is_default)) then
    is_default = (ix_attrib_val == none$)
  endif

case ('MAT6_CALC_METHOD')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, mat6_calc_method_name, lbound(mat6_calc_method_name, 1))
  if (present(is_default)) then
    call default_ele(ele, ele2)
    is_default = (ix_attrib_val == ele2%mat6_calc_method)
  endif

case ('ORIGIN_ELE_REF_PT')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, ref_pt_name, lbound(ref_pt_name, 1))
    if (present(is_default)) is_default = (ix_attrib_val == center_pt$)

case ('PARTICLE')
  attrib_val_name = species_name(ix_attrib_val)
  if (present(is_default)) then
    if (ele%key == photon_fork$) then
      is_default = (ix_attrib_val == photon$)
    else
      is_default = .false. ! Cannot tell so assume the worst.
    endif
  endif

case ('PTC_FIELD_GEOMETRY')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, ptc_field_geometry_name, &
                                                  lbound(ptc_field_geometry_name, 1))
  if (present(is_default)) is_default = (ix_attrib_val == sector$)

case ('PTC_INTEGRATION_TYPE')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, ptc_integration_type_name, &
                                                  lbound(ptc_integration_type_name, 1))
  if (present(is_default)) is_default = (ix_attrib_val == matrix_kick$)

case ('REF_COORDINATES')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, end_at_name(1:2), 1)
  if (present(is_default)) is_default = (ix_attrib_val == exit_end$)

case ('REF_ORBIT_FOLLOWS')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, ref_orbit_follows_name, lbound(ref_orbit_follows_name, 1))
  if (present(is_default)) then
    is_default = (ix_attrib_val == bragg_diffracted$)
  endif

case ('SLAVE_STATUS')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, control_name, lbound(control_name, 1))  
  if (present(is_default)) is_default = .false.
  
case ('SPATIAL_DISTRIBUTION')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, distribution_name, lbound(distribution_name, 1))
  if (present(is_default)) is_default = (ix_attrib_val == gaussian$)

case ('SPIN_TRACKING_METHOD')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, spin_tracking_method_name, lbound(spin_tracking_method_name, 1))
  if (present(is_default)) then
    call default_ele (ele, ele2)
    is_default = (ix_attrib_val == ele2%spin_tracking_method)
  endif

case ('TRACKING_METHOD')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, tracking_method_name, lbound(tracking_method_name, 1))
  if (present(is_default)) then
    call default_ele (ele, ele2)
    is_default = (ix_attrib_val == ele2%tracking_method)
  endif

case ('VELOCITY_DISTRIBUTION')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, distribution_name, lbound(distribution_name, 1))
  if (present(is_default)) is_default = (ix_attrib_val == gaussian$)

case default
  call out_io (s_fatal$, r_name, 'BAD ATTRIBUTE NAME: ' // attrib_name)
end select

!---------------------------------------
contains

subroutine get_this_attrib_name (val_name, ix_attrib_val, name_array, min_arr)

integer ix_attrib_val, min_arr
character(*) val_name
character(*) name_array(min_arr:)

!

if (ix_attrib_val < lbound(name_array, 1) .or. ix_attrib_val > ubound(name_array, 1)) then
  val_name = 'Garbage!'
else
  val_name = name_array(ix_attrib_val)
endif

end subroutine get_this_attrib_name

!---------------------------------------
! contains

! Note: This only works if switch is stored in ele%value(:) array.

function default_value(ele, ix_attrib) result (default_val)

type (ele_struct) ele, ele2
integer ix_attrib, default_val

!

ele2 = ele_struct()
ele2%key = ele%key
ele2%sub_key = ele%key
call set_ele_defaults (ele2, do_allocate = .false.)
default_val = nint(ele2%value(ix_attrib))

end function default_value

!---------------------------------------
! contains


subroutine default_ele(ele, ele_default)

type (ele_struct) ele, ele_default

!

ele_default = ele_struct()
ele_default%key = ele%key
ele_default%sub_key = ele%key
call set_ele_defaults (ele_default, do_allocate = .false.)

end subroutine default_ele

end function switch_attrib_value_name 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function n_attrib_string_max_len () result (max_len)
!
! Routine to return the the maximum number of characters in any attribute
! name known to bmad.
!
! Output:
!   max_len -- Integer: Maximum number of characters in any attribute name.
!-

function n_attrib_string_max_len () result (max_len)

integer max_len
integer, save :: max_length = 0

!

if (attribute_array_init_needed) call init_attribute_name_array
if (max_length == 0) max_length = maxval(len_trim(attrib_array(1:n_key$, 1:num_ele_attrib_extended$)%name))
max_len = max_length

end function

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Fuction ele_has (ele, attrib) result(has_it)
!
! Routine to determine if a given type of lattice element has a particular attribute.
!
! Input:
!   ele     -- ele_struct: Lattice element of a certain type
!   attrib  -- character(*): Name of the attribute. Must be upper case.
!
! Output:
!   has_it  -- logical: True if element has an attribute of that name.
!-

function ele_has (ele, attrib) result (has_it)

type (ele_struct) ele
character(*) attrib
logical has_it

!

has_it = (attribute_index(ele, attrib) > 0)

end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine set_attribute_alias (attrib_name, alias_name, err_flag, lat)
!
! Routine to setup an alias for element attributes like custom_attribute1$, etc. in
! the attribute name table.
!
! Input:
!   attrib_name -- Character(*): Name of attribute to define an alias for.
!   alias_name  -- Character(*): Name of alias. If prefixed by "<class>::" then
!                    the alias will be set only for that element class. Example:
!                    "quadrupole::error" will set the alias only for quadrupoles.
!
! Output:
!   err_flag    -- Logical: Set True if an error. False otherwise.
!   lat         -- lat_struct, optional: If present then the alias information will
!                     be added to the lat%attribute_alias(:) array. This is only 
!                     needed if the lattice is to be written to a file using 
!                     write_bmad_lattice_file
!-

subroutine set_attribute_alias (attrib_name, alias_name, err_flag, lat)

implicit none

type (lat_struct), optional :: lat
integer n
character(*) attrib_name, alias_name
character(40) attrib_str, alias_str 
character(20) :: r_name = 'set_attribute_alias'
logical err_flag

!

if (attribute_array_init_needed) call init_attribute_name_array

err_flag = .false.

attrib_str = upcase(attrib_name)
alias_str = upcase(alias_name)

select case(attrib_str)
case ('CUSTOM_ATTRIBUTE1'); call set_it (custom_attribute1$)
case ('CUSTOM_ATTRIBUTE2'); call set_it (custom_attribute2$)
case ('CUSTOM_ATTRIBUTE3'); call set_it (custom_attribute3$)
case ('CUSTOM_ATTRIBUTE4'); call set_it (custom_attribute4$)
case ('CUSTOM_ATTRIBUTE5'); call set_it (custom_attribute5$)
case default
  err_flag = .true.
  call out_io (s_error$, r_name, 'ATTRIBUTE NAME NOT VALID FOR ALIAS SETUP: ' // attrib_str)
end select

! Add to the lat%attribute_alias

if (present(lat)) then
  if (allocated(lat%attribute_alias)) then
    n = size(lat%attribute_alias) + 1
    call re_allocate(lat%attribute_alias, n)
  else
    n = 1
    allocate(lat%attribute_alias(1))
  endif
  lat%attribute_alias(n) = trim(attrib_str) // '=' // alias_str
endif

!---------------------------------------------------------------
contains

subroutine set_it (ix_attrib)

integer ix_attrib, ix, i, key
character(40) old_attrib, a_str
logical warning_given

! If alias_str is of the form "ele_class::alias_name" then need to split string

a_str = alias_str
key = 0
ix = index(alias_str, '::')
if (ix /= 0) then
  a_str = alias_str(ix+2:)
  key = key_name_to_key_index(alias_str(1:ix-1), .true.)
  if (key < 1) then
    call out_io (s_error$, r_name, 'ELEMENT CLASS NOT RECOGNIZED: ' // alias_str)
    err_flag = .true.
    return
  endif
endif

! Check alias name for invalid characters

if (str_find_first_not_in_set(trim(a_str), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_') > 0) then
  call out_io (s_error$, r_name, 'ALIAS NAME HAS INVALID CHARACTERS: ' // a_str)
  err_flag = .true.
  return
endif

! Set

warning_given = .false.

do i = 1, n_key$
  if (key /= 0 .and. key /= i) cycle
  if (i == def_parameter$ .or. i == def_mad_beam$) cycle
  old_attrib = attribute_name(i, ix_attrib)
  if (old_attrib /= null_name$ .and. old_attrib /= a_str .and. .not. warning_given) then
    if (key == 0) then
      call out_io (s_warn$, r_name, &
        'A CUSTOM_ATTRIBUTE IS BEING REDENFINED: ' // attrib_str, &
        'FROM: ' // old_attrib, &
        'TO:   ' // a_str)
    else
      call out_io (s_warn$, r_name, &
        'A CUSTOM_ATTRIBUTE IS BEING REDENFINED: ' // attrib_str, &
        'FOR ELEMENT CLASS: ' // key_name(i), &
        'FROM: ' // old_attrib, &
        'TO:   ' // a_str)
    endif
  endif

  call init_attribute_name1 (i, ix_attrib, a_str, override = .true.)
enddo

end subroutine set_it

end subroutine set_attribute_alias

end module
