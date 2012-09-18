module basic_attribute_mod

use bmad_struct

! The type of an attribute may be one of:
!   does_not_exist$ -- Does not exist.
!   is_free$        -- Free to vary as long as attribute has not controlled by, eg, overlay_lord.
!   quasi_free$     -- May be free or not. For example, k1 is only free if field_master = F.
!   dependent$      -- Cannot be varied as an independent parameter..
!   private$        -- Internal, parameter used in calculations. Is dependent by definition.

integer, parameter :: does_not_exist$ = -1, is_free$ = 1, quasi_free$ = 2, dependent$ = 3, private$ = 4

type ele_attribute_struct
  character(40) :: name = null_name$
  integer :: type = does_not_exist$
end type

type (ele_attribute_struct), private, save :: attrib_array(n_key, num_ele_attrib_extended$)

character(40), private, save :: short_attrib_array(n_key, num_ele_attrib_extended$)
integer, private, save :: attrib_num(n_key)
integer, private, save :: attrib_ix(n_key, num_ele_attrib_extended$)
logical, private, save :: init_needed = .true.
logical, private, save :: has_orientation_attributes_key(n_key)

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+             
! Function attribute_index (ele, name) result (attrib_index)
!
! Function to return the index of a attribute for a given BMAD element type
! and the name of the attribute. Abbreviations are permitted but must be at 
! least 3 characters.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele  -- Ele_struct: attribute_index will restrict the name search to 
!             valid attributes of the given element. Note: If 
!             ele%key = overlay$ then the entire name table will be searched.
!   name -- Character(40): Attribute name. Must be uppercase.
!
! Output:
!   attrib_index -- Integer: Index of the attribute. If the attribute name
!                            is not appropriate then 0 will be returned.
!
! Example:
!     ele%key = sbend$
!     ix = attribute_index (ele, 'K1')
! Result:
!     ix -> k1$
!-

function attribute_index (ele, name) result (attrib_index)

implicit none

type (ele_struct) ele

integer i, j, k, key, num, ilen, n_abbrev, ix_abbrev
integer attrib_index

character(*) name
character(40) name40

!-----------------------------------------------------------------------

if (init_needed) call init_attribute_name_array

name40 = name           ! make sure we have 40 characters
key = ele%key
attrib_index = 0        ! match not found

ilen = len_trim(name)
if (ilen == 0) return
if (ilen < 3) ilen = 3  ! Need at least three characters.
n_abbrev = 0            ! number of abbreviation matches.

!-----------------------------------------------------------------------
! search for name

! Overlays search all types of elements

if (key == overlay$) then
  do k = 1, n_key
    do i = 1, attrib_num(k)
      if (short_attrib_array(k, i) == name40) then
        attrib_index = attrib_ix(k, i)
        return
      endif
      if (short_attrib_array(k, i)(1:ilen) == name40(1:ilen)) then
        n_abbrev = n_abbrev + 1
        ix_abbrev = attrib_ix(k, i)
      endif 
    enddo
  enddo

  if (name40 == 'CURRENT') then
    attrib_index = current$
    return
  endif

! else only search this type of element

elseif (key > 0 .and. key <= n_key) then
  do i = 1, attrib_num(key)
    if (short_attrib_array(key, i) == name40) then
      attrib_index = attrib_ix(key, i)
      return
    endif
    if (short_attrib_array(key, i)(1:ilen) == name40(1:ilen)) then
      n_abbrev = n_abbrev + 1
      ix_abbrev = attrib_ix(key, i)
    endif 
  enddo      

  if (key == rfcavity$ .and. name40 == 'LAG') then
    attrib_index = phi0$
    return
  endif

! error

else
  print *, 'ERROR IN ATTRIBUTE_INDEX: BAD KEY', key
  if (bmad_status%exit_on_error) call err_exit
endif

! If there is one unique abbreviation then use it.

if (n_abbrev == 1) attrib_index = ix_abbrev

end function attribute_index 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_name (ele, ix_att) result (attrib_name)
!
! Function to return the name of an attribute for a particular type of 
! BMAD element. 
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
!   attrib_name -- Character(40): Name of attribute. 
!      = "!BAD ELE KEY"                 %key is invalid
!      = "!BAD INDEX"                   ix_att is invalid (out of range).
!      = "!INVALID INDEX"               ix_att is invalid for an overlay 
!      = "!NULL" (null_name$)           ix_att does not correspond to an attribute or is private.
!
! Example:
!   ele%key = sbend$
!   name = attribute_name (ele, k1$)
! Result:
!   name -> "K1"
!-

function attribute_name (ele, ix_att) result (attrib_name)

implicit none

type (ele_struct) ele
integer i, key, ix_att
character(40) attrib_name

!

if (init_needed) call init_attribute_name_array()

key = ele%key

if (key <= 0 .or. key > n_key) then
  attrib_name = '!BAD ELE KEY'
elseif (ix_att <= 0 .or. ix_att > num_ele_attrib_extended$) then
  attrib_name = '!BAD INDEX'
elseif (ele%lord_status == overlay_lord$) then
  if (ix_att == ele%ix_value) then
    attrib_name = ele%component_name
  else
    attrib_name = '!INVALID INDEX'
  endif
else
  if (attrib_array(key, ix_att)%type == private$) then
    attrib_name = null_name$
  else
    attrib_name = attrib_array(key, ix_att)%name
  endif
endif

end function attribute_name 

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

implicit none

type (ele_struct) ele
type (ele_attribute_struct) attrib_info
integer i, key, ix_att

!

if (init_needed) call init_attribute_name_array()

attrib_info%type = does_not_exist$

if (ele%key <= 0 .or. ele%key > n_key) then
  attrib_info%name = '!BAD ELE KEY'
elseif (ix_att <= 0 .or. ix_att > num_ele_attrib_extended$) then
  attrib_info%name = '!BAD INDEX'
elseif (ele%lord_status == overlay_lord$) then
  if (ix_att == ele%ix_value) then
    attrib_info%name = ele%component_name
    attrib_info%type = free$
  else
    attrib_info%name = '!INVALID INDEX'
  endif
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

implicit none

integer i, j, num

!

if (.not. init_needed) return

do i = 1, n_key

  call init_attribute_name1 (i, general1$,  'general1', private$)
  call init_attribute_name1 (i, general2$,  'general2', private$)
  call init_attribute_name1 (i, general3$,  'general3', private$)
  call init_attribute_name1 (i, check_sum$, 'check_sum', private$)
  call init_attribute_name1 (i, scratch$,   'scratch', private$)

  if (i == crystal$ .or. i == multilayer_mirror$ .or. i == mirror$) then
    call init_attribute_name1 (i, d_source$,             'D_SOURCE')
    call init_attribute_name1 (i, d_detec$,              'D_DETEC')
    call init_attribute_name1 (i, c2_curve$,             'C2_CURVE')
    call init_attribute_name1 (i, c3_curve$,             'C3_CURVE')
    call init_attribute_name1 (i, c4_curve$,             'C4_CURVE')
    call init_attribute_name1 (i, c2_curve_tot$,         'C2_CURVE_TOT')
    call init_attribute_name1 (i, c3_curve_tot$,         'C3_CURVE_TOT')
    call init_attribute_name1 (i, c4_curve_tot$,         'C4_CURVE_TOT')
  endif

  select case(i)
  case (monitor$, instrument$, marker$, pipe$)
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
    call init_attribute_name1 (i, x_length$,        'X_LENGTH')
    call init_attribute_name1 (i, y_length$,        'Y_LENGTH')
  end select

  if (i == hybrid$)         cycle
  if (i == def_beam$)       cycle
  if (i == def_parameter$)  cycle
  if (i == def_beam_start$) cycle
  if (i == init_ele$) cycle

  call init_attribute_name1 (i, type$,      'TYPE')
  call init_attribute_name1 (i, alias$,     'ALIAS')
  call init_attribute_name1 (i, descrip$,   'DESCRIP')

  if (i == group$)    cycle
  if (i == overlay$)  cycle
  if (i == girder$)   cycle

  call init_attribute_name1 (i, superimpose$,        'SUPERIMPOSE')
  call init_attribute_name1 (i, offset$,             'OFFSET')
  call init_attribute_name1 (i, reference$,          'REFERENCE')
  call init_attribute_name1 (i, ele_beginning$,      'ELE_BEGINNING')
  call init_attribute_name1 (i, ele_center$,         'ELE_CENTER')
  call init_attribute_name1 (i, ele_end$,            'ELE_END')
  call init_attribute_name1 (i, ref_beginning$,      'REF_BEGINNING')
  call init_attribute_name1 (i, ref_center$,         'REF_CENTER')
  call init_attribute_name1 (i, ref_end$,            'REF_END')
  call init_attribute_name1 (i, create_em_field_slave$, 'CREATE_EM_FIELD_SLAVE')

  if (i == null_ele$) cycle

  call init_attribute_name1 (i, wall_attribute$,         'WALL')
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
  call init_attribute_name1 (i, mat6_calc_method$,       'MAT6_CALC_METHOD')
  call init_attribute_name1 (i, tracking_method$,        'TRACKING_METHOD')
  call init_attribute_name1 (i, spin_tracking_method$,   'SPIN_TRACKING_METHOD')

  call init_attribute_name1 (i, E_tot$,                  'E_TOT', quasi_free$) ! Free in multipass_lord
  call init_attribute_name1 (i, p0c$,                    'P0C', quasi_free$)   ! Free in multipass_lord
  call init_attribute_name1 (i, delta_ref_time$,         'DELTA_REF_TIME', dependent$)
  call init_attribute_name1 (i, ref_time_start$,         'ref_time_start', private$)

  if (i == match$) cycle
  if (i == floor_position$) cycle

  call init_attribute_name1 (i, tilt$,          'TILT' )
  call init_attribute_name1 (i, x_offset$,      'X_OFFSET')
  call init_attribute_name1 (i, y_offset$,      'Y_OFFSET')
  call init_attribute_name1 (i, s_offset$,      'S_OFFSET')
  call init_attribute_name1 (i, x_pitch$,       'X_PITCH')
  call init_attribute_name1 (i, y_pitch$,       'Y_PITCH')
  call init_attribute_name1 (i, tilt_tot$,      'TILT_TOT', dependent$)
  call init_attribute_name1 (i, x_offset_tot$,  'X_OFFSET_TOT', dependent$)
  call init_attribute_name1 (i, y_offset_tot$,  'Y_OFFSET_TOT', dependent$)
  call init_attribute_name1 (i, s_offset_tot$,  'S_OFFSET_TOT', dependent$)
  call init_attribute_name1 (i, x_pitch_tot$,   'X_PITCH_TOT', dependent$)
  call init_attribute_name1 (i, y_pitch_tot$,   'Y_PITCH_TOT', dependent$)

  if (i == mirror$)     cycle
  if (i == crystal$)    cycle
  if (i == multilayer_mirror$) cycle
  if (i == capillary$)    cycle

  if (i /= drift$) call init_attribute_name1 (i, is_on$,        'IS_ON')

  if (i == photon_branch$) cycle
  if (i == branch$)       cycle
  if (i == marker$)       cycle
  if (i == patch$)        cycle
  if (i == beambeam$)     cycle
  if (i == hom$)          cycle
  if (i == multipole$)    cycle 
  if (i == ab_multipole$) cycle

  call init_attribute_name1 (i, symplectify$,         'SYMPLECTIFY')
  call init_attribute_name1 (i, map_with_offsets$,    'MAP_WITH_OFFSETS')

  if (i == taylor$)       cycle

  call init_attribute_name1 (i, l$,            'L')
  call init_attribute_name1 (i, l_hard_edge$,  'L_HARD_EDGE', dependent$)

  call init_attribute_name1 (i, integrator_order$,   'INTEGRATOR_ORDER')
  call init_attribute_name1 (i, num_steps$,          'NUM_STEPS', quasi_free$)
  call init_attribute_name1 (i, ds_step$,            'DS_STEP')
  call init_attribute_name1 (i, csr_calc_on$,        'CSR_CALC_ON')
  call init_attribute_name1 (i, n_ref_pass$,         'N_REF_PASS')
  call init_attribute_name1 (i, field_scale$,        'FIELD_SCALE', quasi_free$)

  if (i == hkicker$)      cycle
  if (i == vkicker$)      cycle
  if (i == custom$)       cycle
  if (i == e_gun$)        cycle
  if (i == em_field$)     cycle
  if (i == drift$)        cycle

  call init_attribute_name1 (i, hkick$,     'HKICK', quasi_free$)
  call init_attribute_name1 (i, vkick$,     'VKICK', quasi_free$)
  call init_attribute_name1 (i, bl_hkick$,  'BL_HKICK', quasi_free$)
  call init_attribute_name1 (i, bl_vkick$,  'BL_VKICK', quasi_free$)

enddo

!

do i = 1, n_key
  select case (i)
  case (elseparator$, kicker$, octupole$, quadrupole$, sbend$, rbend$, &
         sextupole$, solenoid$, sol_quad$, ab_multipole$, wiggler$, bend_sol_quad$, &
         hkicker$, vkicker$)
    attrib_array(i, a0$:a20$)%name = [ 'A0 ', &
                                   'A1 ', 'A2 ', 'A3 ', 'A4 ', 'A5 ', & 
                                   'A6 ', 'A7 ', 'A8 ', 'A9 ', 'A10', &
                                   'A11', 'A12', 'A13', 'A14', 'A15', &
                                   'A16', 'A17', 'A18', 'A19', 'A20' ]
    attrib_array(i, b0$:b20$)%name = [ 'B0 ', &
                                   'B1 ', 'B2 ', 'B3 ', 'B4 ', 'B5 ', & 
                                   'B6 ', 'B7 ', 'B8 ', 'B9 ', 'B10', &
                                   'B11', 'B12', 'B13', 'B14', 'B15', &
                                   'B16', 'B17', 'B18', 'B19', 'B20' ]
    attrib_array(i, a0$:b20$)%type = is_free$
    if (i == ab_multipole$) cycle
    call init_attribute_name1 (i, scale_multipoles$,  'SCALE_MULTIPOLES')
  end select
enddo



!

call init_attribute_name1 (photon_branch$, ix_branch_to$,         'ix_branch_to', private$)
call init_attribute_name1 (photon_branch$, direction$,            'DIRECTION')
call init_attribute_name1 (photon_branch$, to$,                   'TO')
call init_attribute_name1 (photon_branch$, particle$,             'PARTICLE')
call init_attribute_name1 (photon_branch$, lattice_type$,         'LATTICE_TYPE')
call init_attribute_name1 (photon_branch$, E_tot_start$,          'E_TOT_START')
call init_attribute_name1 (photon_branch$, p0c_start$,            'P0C_START')

attrib_array(branch$, :) = attrib_array(photon_branch$, :)

call init_attribute_name1 (init_ele$, delta_ref_time$,              'delta_ref_time', private$)
call init_attribute_name1 (init_ele$, ref_time_start$,              'ref_time_start', private$)
call init_attribute_name1 (init_ele$, e_tot_start$,                 'E_TOT_START')
call init_attribute_name1 (init_ele$, p0c_start$,                   'P0C_START')
call init_attribute_name1 (init_ele$, e_tot$,                       'E_TOT')
call init_attribute_name1 (init_ele$, p0c$,                         'P0C')
call init_attribute_name1 (init_ele$, x_position$,                  'X_POSITION')
call init_attribute_name1 (init_ele$, y_position$,                  'Y_POSITION')
call init_attribute_name1 (init_ele$, z_position$,                  'Z_POSITION')
call init_attribute_name1 (init_ele$, theta_position$,              'THETA_POSITION')
call init_attribute_name1 (init_ele$, phi_position$,                'PHI_POSITION')
call init_attribute_name1 (init_ele$, psi_position$,                'PSI_POSITION')
call init_attribute_name1 (init_ele$, beta_a$,                      'BETA_A')
call init_attribute_name1 (init_ele$, beta_b$,                      'BETA_B')
call init_attribute_name1 (init_ele$, alpha_a$,                     'ALPHA_A')
call init_attribute_name1 (init_ele$, alpha_b$,                     'ALPHA_B')
call init_attribute_name1 (init_ele$, eta_x$,                       'ETA_X')
call init_attribute_name1 (init_ele$, eta_y$,                       'ETA_Y')
call init_attribute_name1 (init_ele$, eta_z$,                       'ETA_Z')
call init_attribute_name1 (init_ele$, etap_x$,                      'ETAP_X')
call init_attribute_name1 (init_ele$, etap_y$,                      'ETAP_Y')
call init_attribute_name1 (init_ele$, phi_a$,                       'PHI_A')
call init_attribute_name1 (init_ele$, phi_b$,                       'PHI_B')
call init_attribute_name1 (init_ele$, cmat_11$,                     'CMAT_11')
call init_attribute_name1 (init_ele$, cmat_12$,                     'CMAT_12')
call init_attribute_name1 (init_ele$, cmat_21$,                     'CMAT_21')
call init_attribute_name1 (init_ele$, cmat_22$,                     'CMAT_22')
call init_attribute_name1 (init_ele$, s_long$,                      'S')
call init_attribute_name1 (init_ele$, ref_time$,                    'REF_TIME')
call init_attribute_name1 (init_ele$, e_field_x$,                   'E_FIELD_X')
call init_attribute_name1 (init_ele$, e_field_y$,                   'E_FIELD_Y')
call init_attribute_name1 (init_ele$, phase_x$,                     'PHASE_X')
call init_attribute_name1 (init_ele$, phase_y$,                     'PHASE_Y')
call init_attribute_name1 (init_ele$, wall_attribute$,              'WALL')

call init_attribute_name1 (def_parameter$, general1$,               'GENERAL1', override = .true.)
call init_attribute_name1 (def_parameter$, general2$,               'GENERAL2', override = .true.)
call init_attribute_name1 (def_parameter$, general3$,               'GENERAL3', override = .true.)
call init_attribute_name1 (def_parameter$, e_tot$,                  'E_TOT')
call init_attribute_name1 (def_parameter$, p0c$,                    'P0C')
call init_attribute_name1 (def_parameter$, lattice_type$,           'LATTICE_TYPE')
call init_attribute_name1 (def_parameter$, lattice$,                'LATTICE')
call init_attribute_name1 (def_parameter$, taylor_order$,           'TAYLOR_ORDER')
call init_attribute_name1 (def_parameter$, ran_seed$,               'RAN_SEED')
call init_attribute_name1 (def_parameter$, n_part$,                 'N_PART')
call init_attribute_name1 (def_parameter$, particle$,               'PARTICLE')
call init_attribute_name1 (def_parameter$, aperture_limit_on$,      'APERTURE_LIMIT_ON')
call init_attribute_name1 (def_parameter$, root_branch_name$,       'ROOT_BRANCH_NAME') 
call init_attribute_name1 (def_parameter$, no_end_marker$,          'NO_END_MARKER')
call init_attribute_name1 (def_parameter$, absolute_time_tracking$, 'ABSOLUTE_TIME_TRACKING')
call init_attribute_name1 (def_parameter$, use_ptc_layout$,         'USE_PTC_LAYOUT')
call init_attribute_name1 (def_parameter$, rf_auto_scale_phase$,    'RF_AUTO_SCALE_PHASE')
call init_attribute_name1 (def_parameter$, rf_auto_scale_amp$,      'RF_AUTO_SCALE_AMP')
call init_attribute_name1 (def_parameter$, ptc_exact_model$,        'PTC_EXACT_MODEL')
call init_attribute_name1 (def_parameter$, ptc_exact_misalign$,     'PTC_EXACT_MISALIGN')

call init_attribute_name1 (def_beam$, particle$,                    'PARTICLE')
call init_attribute_name1 (def_beam$, e_tot$,                       'ENERGY')
call init_attribute_name1 (def_beam$, p0c$,                         'PC')
call init_attribute_name1 (def_beam$, n_part$,                      'N_PART')

call init_attribute_name1 (def_beam_start$, x$,                     'X')
call init_attribute_name1 (def_beam_start$, px$,                    'PX')
call init_attribute_name1 (def_beam_start$, y$,                     'Y')
call init_attribute_name1 (def_beam_start$, py$,                    'PY')
call init_attribute_name1 (def_beam_start$, z$,                     'Z')
call init_attribute_name1 (def_beam_start$, pz$,                    'PZ')
call init_attribute_name1 (def_beam_start$, e_field_x$,             'E_FIELD_X')
call init_attribute_name1 (def_beam_start$, e_field_y$,             'E_FIELD_Y')
call init_attribute_name1 (def_beam_start$, phase_x$,               'PHASE_X')
call init_attribute_name1 (def_beam_start$, phase_y$,               'PHASE_Y')
call init_attribute_name1 (def_beam_start$, t$,                     'T')

call init_attribute_name1 (e_gun$, e_tot_start$,                    'e_tot_start', private$)
call init_attribute_name1 (e_gun$, p0c_start$,                      'p0c_start', private$)
call init_attribute_name1 (e_gun$, voltage$,                        'VOLTAGE')
call init_attribute_name1 (e_gun$, voltage_err$,                    'VOLTAGE_ERR')
call init_attribute_name1 (e_gun$, gradient$,                       'GRADIENT')
call init_attribute_name1 (e_gun$, gradient_err$,                   'GRADIENT_ERR')
call init_attribute_name1 (e_gun$, field_calc$,                     'FIELD_CALC')
call init_attribute_name1 (e_gun$, field$,                          'FIELD')

call init_attribute_name1 (em_field$, e_tot_start$,                 'E_TOT_START')
call init_attribute_name1 (em_field$, p0c_start$,                   'P0C_START')
call init_attribute_name1 (em_field$, field$,                       'FIELD')

call init_attribute_name1 (taylor$, l$,                             'L')
call init_attribute_name1 (taylor$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (taylor$, p0c_start$,                     'p0c_start', private$)

call init_attribute_name1 (marker$, l$,                             'l', private$)
call init_attribute_name1 (marker$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (marker$, p0c_start$,                     'p0c_start', private$)

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
call init_attribute_name1 (match$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (match$, p0c_start$,                     'p0c_start', private$)

call init_attribute_name1 (girder$, l$,                             'L')
call init_attribute_name1 (girder$, x_offset$,                      'X_OFFSET')
call init_attribute_name1 (girder$, y_offset$,                      'Y_OFFSET')
call init_attribute_name1 (girder$, s_offset$,                      'S_OFFSET')
call init_attribute_name1 (girder$, x_pitch$,                       'X_PITCH')
call init_attribute_name1 (girder$, y_pitch$,                       'Y_PITCH')
call init_attribute_name1 (girder$, s_min$,                         'S_MIN')
call init_attribute_name1 (girder$, s_max$,                         'S_MAX')
call init_attribute_name1 (girder$, s_center$,                      'S_CENTER')
call init_attribute_name1 (girder$, tilt$,                          'TILT')

call init_attribute_name1 (lcavity$, p0c_start$,                    'P0C_START')
call init_attribute_name1 (lcavity$, e_tot_start$,                  'E_TOT_START')
call init_attribute_name1 (lcavity$, dphi0$,                        'DPHI0')
call init_attribute_name1 (lcavity$, phi0$,                         'PHI0')
call init_attribute_name1 (lcavity$, gradient$,                     'GRADIENT')
call init_attribute_name1 (lcavity$, rf_frequency$,                 'RF_FREQUENCY')
call init_attribute_name1 (lcavity$, e_loss$,                       'E_LOSS')
call init_attribute_name1 (lcavity$, voltage$,                      'VOLTAGE', quasi_free$)
call init_attribute_name1 (lcavity$, sr_wake_file$,                 'SR_WAKE_FILE')
call init_attribute_name1 (lcavity$, lr_wake_file$,                 'LR_WAKE_FILE')
call init_attribute_name1 (lcavity$, field_calc$,                   'FIELD_CALC')
call init_attribute_name1 (lcavity$, field_master$,                 'FIELD_MASTER')
call init_attribute_name1 (lcavity$, lr_freq_spread$,               'LR_FREQ_SPREAD')
call init_attribute_name1 (lcavity$, coupler_strength$,             'COUPLER_STRENGTH')
call init_attribute_name1 (lcavity$, coupler_angle$,                'COUPLER_ANGLE')
call init_attribute_name1 (lcavity$, coupler_phase$,                'COUPLER_PHASE')
call init_attribute_name1 (lcavity$, coupler_at$,                   'COUPLER_AT')
call init_attribute_name1 (lcavity$, gradient_err$,                 'GRADIENT_ERR')
call init_attribute_name1 (lcavity$, phi0_err$,                     'PHI0_ERR')
call init_attribute_name1 (lcavity$, field$,                        'FIELD')
call init_attribute_name1 (lcavity$, dphi0_ref$,                    'DPHI0_REF', quasi_free$)
call init_attribute_name1 (lcavity$, n_cell$,                       'N_CELL')
call init_attribute_name1 (lcavity$, grad_loss_sr_wake$,            'grad_loss_sr_wake', private$)


call init_attribute_name1 (group$, command$,                        'COMMAND')
call init_attribute_name1 (group$, old_command$,                    'OLD_COMMAND')
call init_attribute_name1 (group$, coef$,                           'COEF')
call init_attribute_name1 (group$, start_edge$,                     'START_EDGE')
call init_attribute_name1 (group$, end_edge$,                       'END_EDGE')
call init_attribute_name1 (group$, accordion_edge$,                 'ACCORDION_EDGE')
call init_attribute_name1 (group$, symmetric_edge$,                 'SYMMETRIC_EDGE')

call init_attribute_name1 (drift$, field_calc$,                     'FIELD_CALC')
call init_attribute_name1 (drift$, field_master$,                   'FIELD_MASTER')
call init_attribute_name1 (drift$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (drift$, p0c_start$,                     'p0c_start', private$)

call init_attribute_name1 (monitor$, field_master$,                 'FIELD_MASTER')
call init_attribute_name1 (monitor$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (monitor$, p0c_start$,                     'p0c_start', private$)

attrib_array(instrument$, :)                         = attrib_array(monitor$, :)
attrib_array(pipe$, :)                               = attrib_array(monitor$, :)

call init_attribute_name1 (hkicker$, kick$,                         'KICK', quasi_free$)
call init_attribute_name1 (hkicker$, field_calc$,                   'FIELD_CALC')
call init_attribute_name1 (hkicker$, field_master$,                 'FIELD_MASTER')
call init_attribute_name1 (hkicker$, bl_kick$,                      'BL_KICK', quasi_free$)
call init_attribute_name1 (hkicker$, pole_radius$,                  'POLE_RADIUS')
call init_attribute_name1 (hkicker$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (hkicker$, p0c_start$,                     'p0c_start', private$)

attrib_array(vkicker$, :) = attrib_array(hkicker$, :)

call init_attribute_name1 (kicker$, h_displace$,                    'H_DISPLACE')
call init_attribute_name1 (kicker$, v_displace$,                    'V_DISPLACE')
call init_attribute_name1 (kicker$, radius$,                        'RADIUS')
call init_attribute_name1 (kicker$, field_calc$,                    'FIELD_CALC')
call init_attribute_name1 (kicker$, field_master$,                  'FIELD_MASTER')
call init_attribute_name1 (kicker$, pole_radius$,                   'POLE_RADIUS')
call init_attribute_name1 (kicker$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (kicker$, p0c_start$,                     'p0c_start', private$)

call init_attribute_name1 (sbend$, angle$,                          'ANGLE', quasi_free$)
call init_attribute_name1 (sbend$, exact_fringe$,                   'EXACT_FRINGE')
call init_attribute_name1 (sbend$, e1$,                             'E1')
call init_attribute_name1 (sbend$, e2$,                             'E2')
call init_attribute_name1 (sbend$, h1$,                             'H1')
call init_attribute_name1 (sbend$, h2$,                             'H2')
call init_attribute_name1 (sbend$, k1$,                             'K1', quasi_free$)
call init_attribute_name1 (sbend$, k2$,                             'K2', quasi_free$)
call init_attribute_name1 (sbend$, g$,                              'G', quasi_free$)
call init_attribute_name1 (sbend$, g_err$,                          'G_ERR', quasi_free$)
call init_attribute_name1 (sbend$, roll$,                           'ROLL')
call init_attribute_name1 (sbend$, hgap$,                           'HGAP')
call init_attribute_name1 (sbend$, hgapx$,                          'HGAPX')
call init_attribute_name1 (sbend$, fint$,                           'FINT')
call init_attribute_name1 (sbend$, fintx$,                          'FINTX')
call init_attribute_name1 (sbend$, rho$,                            'RHO', quasi_free$)
call init_attribute_name1 (sbend$, l_chord$,                        'L_CHORD', quasi_free$)
call init_attribute_name1 (sbend$, b_field$,                        'B_FIELD', quasi_free$)
call init_attribute_name1 (sbend$, b_field_err$,                    'B_FIELD_ERR', quasi_free$)
call init_attribute_name1 (sbend$, b1_gradient$,                    'B1_GRADIENT', quasi_free$)
call init_attribute_name1 (sbend$, b2_gradient$,                    'B2_GRADIENT', quasi_free$)
call init_attribute_name1 (sbend$, radius$,                         'RADIUS')
call init_attribute_name1 (sbend$, field_calc$,                     'FIELD_CALC')
call init_attribute_name1 (sbend$, field_master$,                   'FIELD_MASTER')
call init_attribute_name1 (sbend$, ref_orbit$,                      'REF_ORBIT')
call init_attribute_name1 (sbend$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (sbend$, p0c_start$,                     'p0c_start', private$)

attrib_array(rbend$, :) = attrib_array(sbend$, :)

call init_attribute_name1 (bend_sol_quad$, angle$,                  'ANGLE')
call init_attribute_name1 (bend_sol_quad$, k1$,                     'K1')
call init_attribute_name1 (bend_sol_quad$, g$,                      'G')
call init_attribute_name1 (bend_sol_quad$, ks$,                     'KS')
call init_attribute_name1 (bend_sol_quad$, dks_ds$,                 'DKS_DS')
call init_attribute_name1 (bend_sol_quad$, quad_tilt$,              'QUAD_TILT')
call init_attribute_name1 (bend_sol_quad$, bend_tilt$,              'BEND_TILT')
call init_attribute_name1 (bend_sol_quad$, x_quad$,                 'X_QUAD')
call init_attribute_name1 (bend_sol_quad$, y_quad$,                 'Y_QUAD')
call init_attribute_name1 (bend_sol_quad$, rho$,                    'RHO')
call init_attribute_name1 (bend_sol_quad$, radius$,                 'RADIUS')
call init_attribute_name1 (bend_sol_quad$, field_calc$,             'FIELD_CALC')
call init_attribute_name1 (bend_sol_quad$, field_master$,           'FIELD_MASTER')
call init_attribute_name1 (bend_sol_quad$, field$,                  'FIELD')
call init_attribute_name1 (bend_sol_quad$, E_tot_start$,            'E_tot_start', private$)
call init_attribute_name1 (bend_sol_quad$, p0c_start$,              'p0c_start', private$)

call init_attribute_name1 (patch$, l$,                              'L')
call init_attribute_name1 (patch$, t_offset$,                       'T_OFFSET')
call init_attribute_name1 (patch$, p0c_start$,                      'P0C_START')
call init_attribute_name1 (patch$, e_tot_start$,                    'E_TOT_START')
call init_attribute_name1 (patch$, e_tot_offset$,                   'E_TOT_OFFSET')
call init_attribute_name1 (patch$, ref_orbit$,                      'REF_ORBIT')
call init_attribute_name1 (patch$, ref_patch$,                      'REF_PATCH')
call init_attribute_name1 (patch$, n_ref_pass$,                     'N_REF_PASS')
call init_attribute_name1 (patch$, translate_after$,                'TRANSLATE_AFTER')
call init_attribute_name1 (patch$, z_offset$,                       'Z_OFFSET', override = .true.)

call init_attribute_name1 (floor_position$, l$,                           'l', private$)
call init_attribute_name1 (floor_position$, x_position$,                  'X_POSITION')
call init_attribute_name1 (floor_position$, y_position$,                  'Y_POSITION')
call init_attribute_name1 (floor_position$, z_position$,                  'Z_POSITION')
call init_attribute_name1 (floor_position$, theta_position$,              'THETA_POSITION')
call init_attribute_name1 (floor_position$, phi_position$,                'PHI_POSITION')
call init_attribute_name1 (floor_position$, psi_position$,                'PSI_POSITION')


call init_attribute_name1 (quadrupole$, k1$,                        'K1', quasi_free$)
call init_attribute_name1 (quadrupole$, B1_gradient$,               'B1_GRADIENT', quasi_free$)
call init_attribute_name1 (quadrupole$, include_fringe$,            'INCLUDE_FRINGE')
call init_attribute_name1 (quadrupole$, radius$,                    'RADIUS')
call init_attribute_name1 (quadrupole$, field_calc$,                'FIELD_CALC')
call init_attribute_name1 (quadrupole$, field_master$,              'FIELD_MASTER')
call init_attribute_name1 (quadrupole$, pole_radius$,               'POLE_RADIUS')
call init_attribute_name1 (quadrupole$, field$,                     'FIELD')
call init_attribute_name1 (quadrupole$, E_tot_start$,               'E_tot_start', private$)
call init_attribute_name1 (quadrupole$, p0c_start$,                 'p0c_start', private$)

call init_attribute_name1 (sextupole$, k2$,                         'K2', quasi_free$)
call init_attribute_name1 (sextupole$, B2_gradient$,                'B2_GRADIENT', quasi_free$)
call init_attribute_name1 (sextupole$, radius$,                     'RADIUS')
call init_attribute_name1 (sextupole$, field_calc$,                 'FIELD_CALC')
call init_attribute_name1 (sextupole$, field_master$,               'FIELD_MASTER')
call init_attribute_name1 (sextupole$, pole_radius$,                'POLE_RADIUS')
call init_attribute_name1 (sextupole$, field$,                      'FIELD')
call init_attribute_name1 (sextupole$, E_tot_start$,                'E_tot_start', private$)
call init_attribute_name1 (sextupole$, p0c_start$,                  'p0c_start', private$)

call init_attribute_name1 (octupole$, k3$,                          'K3', quasi_free$)
call init_attribute_name1 (octupole$, B3_gradient$,                 'B3_GRADIENT', quasi_free$)
call init_attribute_name1 (octupole$, radius$,                      'RADIUS')
call init_attribute_name1 (octupole$, field_calc$,                  'FIELD_CALC')
call init_attribute_name1 (octupole$, field_master$,                'FIELD_MASTER')
call init_attribute_name1 (octupole$, pole_radius$,                 'POLE_RADIUS')
call init_attribute_name1 (octupole$, field$,                       'FIELD')
call init_attribute_name1 (octupole$, E_tot_start$,                 'E_tot_start', private$)
call init_attribute_name1 (octupole$, p0c_start$,                   'p0c_start', private$)

call init_attribute_name1 (solenoid$, ks$,                          'KS', quasi_free$)
call init_attribute_name1 (solenoid$, bs_field$,                    'BS_FIELD', quasi_free$)
call init_attribute_name1 (solenoid$, radius$,                      'RADIUS')
call init_attribute_name1 (solenoid$, field_calc$,                  'FIELD_CALC')
call init_attribute_name1 (solenoid$, field_master$,                'FIELD_MASTER')
call init_attribute_name1 (solenoid$, pole_radius$,                 'POLE_RADIUS')
call init_attribute_name1 (solenoid$, field$,                       'FIELD')
call init_attribute_name1 (solenoid$, E_tot_start$,                 'E_tot_start', private$)
call init_attribute_name1 (solenoid$, p0c_start$,                   'p0c_start', private$)

call init_attribute_name1 (rfcavity$, dphi0$,                       'DPHI0')
call init_attribute_name1 (rfcavity$, voltage$,                     'VOLTAGE')
call init_attribute_name1 (rfcavity$, rf_frequency$,                'RF_FREQUENCY')
call init_attribute_name1 (rfcavity$, phi0$,                        'PHI0')
call init_attribute_name1 (rfcavity$, phi0_err$,                    'phi0_err', private$)
call init_attribute_name1 (rfcavity$, harmon$,                      'HARMON')
call init_attribute_name1 (rfcavity$, field_calc$,                  'FIELD_CALC')
call init_attribute_name1 (rfcavity$, field_master$,                'FIELD_MASTER')
call init_attribute_name1 (rfcavity$, sr_wake_file$,                'SR_WAKE_FILE')
call init_attribute_name1 (rfcavity$, lr_wake_file$,                'LR_WAKE_FILE')
call init_attribute_name1 (rfcavity$, lr_freq_spread$,              'LR_FREQ_SPREAD')
call init_attribute_name1 (rfcavity$, coupler_strength$,            'COUPLER_STRENGTH')
call init_attribute_name1 (rfcavity$, coupler_angle$,               'COUPLER_ANGLE')
call init_attribute_name1 (rfcavity$, coupler_phase$,               'COUPLER_PHASE')
call init_attribute_name1 (rfcavity$, coupler_at$,                  'COUPLER_AT')
call init_attribute_name1 (rfcavity$, field$,                       'FIELD')
call init_attribute_name1 (rfcavity$, dphi0_ref$,                   'DPHI0_REF')
call init_attribute_name1 (rfcavity$, n_cell$,                      'N_CELL')
call init_attribute_name1 (rfcavity$, grad_loss_sr_wake$,           'grad_loss_sr_wake', private$)
call init_attribute_name1 (rfcavity$, dphi0_max$,                   'dphi0_max', private$)
call init_attribute_name1 (rfcavity$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (rfcavity$, p0c_start$,                     'p0c_start', private$)

call init_attribute_name1 (elseparator$, gap$,                      'GAP')
call init_attribute_name1 (elseparator$, e_field$,                  'E_FIELD', quasi_free$)
call init_attribute_name1 (elseparator$, voltage$,                  'VOLTAGE', quasi_free$)
call init_attribute_name1 (elseparator$, radius$,                   'RADIUS')
call init_attribute_name1 (elseparator$, field_calc$,               'FIELD_CALC')
call init_attribute_name1 (elseparator$, field_master$,             'FIELD_MASTER')
call init_attribute_name1 (elseparator$, field$,                    'FIELD')
call init_attribute_name1 (elseparator$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (elseparator$, p0c_start$,                     'p0c_start', private$)

call init_attribute_name1 (beambeam$, l$,                           'l', private$)
call init_attribute_name1 (beambeam$, sig_x$,                       'SIG_X')
call init_attribute_name1 (beambeam$, sig_y$,                       'SIG_Y')
call init_attribute_name1 (beambeam$, sig_z$,                       'SIG_Z')
call init_attribute_name1 (beambeam$, bbi_const$,                   'BBI_CONSTANT', dependent$)
call init_attribute_name1 (beambeam$, charge$,                      'CHARGE')
call init_attribute_name1 (beambeam$, n_slice$,                     'N_SLICE')
call init_attribute_name1 (beambeam$, symplectify$,                 'N_SLICE')
call init_attribute_name1 (beambeam$, field_calc$,                  'FIELD_CALC')
call init_attribute_name1 (beambeam$, field_master$,                'FIELD_MASTER')
call init_attribute_name1 (beambeam$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (beambeam$, p0c_start$,                     'p0c_start', private$)

call init_attribute_name1 (wiggler$, k1$,                           'K1', dependent$)
call init_attribute_name1 (wiggler$, l_pole$,                       'L_POLE')
call init_attribute_name1 (wiggler$, b_max$,                        'B_MAX')
call init_attribute_name1 (wiggler$, rho$,                          'RHO', dependent$)
call init_attribute_name1 (wiggler$, n_pole$,                       'N_POLE')
call init_attribute_name1 (wiggler$, term$,                         'TERM')
call init_attribute_name1 (wiggler$, polarity$,                     'POLARITY')
call init_attribute_name1 (wiggler$, radius$,                       'RADIUS')
call init_attribute_name1 (wiggler$, field_calc$,                   'FIELD_CALC')
call init_attribute_name1 (wiggler$, field_master$,                 'FIELD_MASTER')
call init_attribute_name1 (wiggler$, x_ray_line_len$,               'X_RAY_LINE_LEN')
call init_attribute_name1 (wiggler$, field$,                        'FIELD')
call init_attribute_name1 (wiggler$, E_tot_start$,                   'E_tot_start', private$)
call init_attribute_name1 (wiggler$, p0c_start$,                     'p0c_start', private$)

call init_attribute_name1 (sol_quad$, k1$,                          'K1', quasi_free$)
call init_attribute_name1 (sol_quad$, ks$,                          'KS', quasi_free$)
call init_attribute_name1 (sol_quad$, b1_gradient$,                 'B1_GRADIENT', quasi_free$)
call init_attribute_name1 (sol_quad$, bs_field$,                    'BS_FIELD', quasi_free$)
call init_attribute_name1 (sol_quad$, radius$,                      'RADIUS')
call init_attribute_name1 (sol_quad$, field_calc$,                  'FIELD_CALC')
call init_attribute_name1 (sol_quad$, field_master$,                'FIELD_MASTER')
call init_attribute_name1 (sol_quad$, field$,                       'FIELD')
call init_attribute_name1 (sol_quad$, E_tot_start$,                 'E_tot_start', private$)
call init_attribute_name1 (sol_quad$, p0c_start$,                   'p0c_start', private$)

call init_attribute_name1 (multipole$, l$,                          'L')
attrib_array(multipole$, k0l$:k20l$)%name    = &
             ['K0L ', 'K1L ', 'K2L ', 'K3L ', 'K4L ', 'K5L ', 'K6L ', 'K7L ', 'K8L ', 'K9L ', 'K10L', &
                      'K11L', 'K12L', 'K13L', 'K14L', 'K15L', 'K16L', 'K17L', 'K18L', 'K19L', 'K20L']
attrib_array(multipole$, t0$:t20$)%name = ['T0 ', &
                               'T1 ', 'T2 ', 'T3 ', 'T4 ', 'T5 ', & 
                               'T6 ', 'T7 ', 'T8 ', 'T9 ', 'T10', &
                               'T11', 'T12', 'T13', 'T14', 'T15', &
                               'T16', 'T17', 'T18', 'T19', 'T20' ]
attrib_array(multipole$, k0l$:t20$)%type = is_free$

call init_attribute_name1 (ab_multipole$, l$,                       'L')

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
call init_attribute_name1 (custom$, e_tot_start$,                   'E_TOT_START')
call init_attribute_name1 (custom$, p0c_start$,                     'P0C_START')

call init_attribute_name1 (hybrid$, l$,                             'L')
call init_attribute_name1 (hybrid$, delta_e$,                       'DELTA_E')
call init_attribute_name1 (hybrid$, delta_ref_time$,                'DELTA_REF_TIME')
call init_attribute_name1 (hybrid$, ref_time_start$,                'ref_time_start', private$)
call init_attribute_name1 (hybrid$, e_tot_start$,                   'E_TOT_START')
call init_attribute_name1 (hybrid$, p0c_start$,                     'P0C_START')

call init_attribute_name1 (mirror$, l$,                             'l', private$)
call init_attribute_name1 (mirror$, graze_angle$,                   'GRAZE_ANGLE')
call init_attribute_name1 (mirror$, graze_angle_err$,               'GRAZE_ANGLE_ERR')
call init_attribute_name1 (mirror$, critical_angle$,                'CRITICAL_ANGLE')
call init_attribute_name1 (mirror$, tilt_err$,                      'TILT_ERR')
call init_attribute_name1 (mirror$, g_trans$,                       'G_TRANS')
call init_attribute_name1 (mirror$, ref_wavelength$,                'REF_WAVELENGTH')

call init_attribute_name1 (multilayer_mirror$, l$,                    'l', private$)
call init_attribute_name1 (multilayer_mirror$, graze_angle$,          'GRAZE_ANGLE')
call init_attribute_name1 (multilayer_mirror$, graze_angle_err$,      'GRAZE_ANGLE_ERR')
call init_attribute_name1 (multilayer_mirror$, tilt_err$,             'TILT_ERR')
call init_attribute_name1 (multilayer_mirror$, n_cell$,               'N_CELL')
call init_attribute_name1 (multilayer_mirror$, d1_thickness$,         'D1_THICKNESS')
call init_attribute_name1 (multilayer_mirror$, d2_thickness$,         'D2_THICKNESS')
call init_attribute_name1 (multilayer_mirror$, v1_unitcell$,          'V1_UNITCELL')
call init_attribute_name1 (multilayer_mirror$, v2_unitcell$,          'V2_UNITCELL')
call init_attribute_name1 (multilayer_mirror$, f0_re1$,               'F0_RE1')
call init_attribute_name1 (multilayer_mirror$, f0_im1$,               'F0_IM1')
call init_attribute_name1 (multilayer_mirror$, f0_re2$,               'F0_RE2')
call init_attribute_name1 (multilayer_mirror$, f0_im2$,               'F0_IM2')
call init_attribute_name1 (multilayer_mirror$, ref_wavelength$,       'REF_WAVELENGTH')
call init_attribute_name1 (multilayer_mirror$, crystal_type$,         'CRYSTAL_TYPE')
call init_attribute_name1 (multilayer_mirror$, ref_polarization$,     'REF_POLARIZATION')
call init_attribute_name1 (multilayer_mirror$, negative_graze_angle$, 'NEGATIVE_GRAZE_ANGLE')

call init_attribute_name1 (crystal$, l$,                            'l', private$)
call init_attribute_name1 (crystal$, graze_angle_in$,               'GRAZE_ANGLE_IN')
call init_attribute_name1 (crystal$, graze_angle_out$,              'GRAZE_ANGLE_OUT')
call init_attribute_name1 (crystal$, graze_angle_err$,              'GRAZE_ANGLE_ERR')
call init_attribute_name1 (crystal$, psi_angle$,                    'PSI_ANGLE')
call init_attribute_name1 (crystal$, alpha_angle$,                  'ALPHA_ANGLE')
call init_attribute_name1 (crystal$, tilt_err$,                     'TILT_ERR')
call init_attribute_name1 (crystal$, tilt_corr$,                    'TILT_CORR')
call init_attribute_name1 (crystal$, d_spacing$,                    'D_SPACING')
call init_attribute_name1 (crystal$, v_unitcell$,                   'V_UNITCELL')
call init_attribute_name1 (crystal$, f0_re$,                        'F0_RE')
call init_attribute_name1 (crystal$, f0_im$,                        'F0_IM')
call init_attribute_name1 (crystal$, fh_re$,                        'FH_RE')
call init_attribute_name1 (crystal$, fh_im$,                        'FH_IM')
call init_attribute_name1 (crystal$, b_param$,                      'B_PARAM')
call init_attribute_name1 (crystal$, bragg_angle$,                  'BRAGG_ANGLE' )
call init_attribute_name1 (crystal$, g_trans$,                      'G_TRANS')
call init_attribute_name1 (crystal$, ref_wavelength$,               'REF_WAVELENGTH')
call init_attribute_name1 (crystal$, diffraction_type$,             'DIFFRACTION_TYPE')  ! Note: Not currently used.
call init_attribute_name1 (crystal$, crystal_type$,                 'CRYSTAL_TYPE')
call init_attribute_name1 (crystal$, thickness$,                    'THICKNESS')
call init_attribute_name1 (crystal$, follow_diffracted_beam$,       'FOLLOW_DIFFRACTED_BEAM')
call init_attribute_name1 (crystal$, kh_x_norm$,                    'kh_x_norm', private$)
call init_attribute_name1 (crystal$, kh_y_norm$,                    'kh_y_norm', private$)
call init_attribute_name1 (crystal$, kh_z_norm$,                    'kh_z_norm', private$)
call init_attribute_name1 (crystal$, l_x$,                          'l_x', private$)
call init_attribute_name1 (crystal$, l_y$,                          'l_y', private$)
call init_attribute_name1 (crystal$, l_z$,                          'l_z', private$)
call init_attribute_name1 (crystal$, ref_polarization$,             'REF_POLARIZATION')
call init_attribute_name1 (crystal$, ref_cap_gamma$,                'REF_CAP_GAMMA')
call init_attribute_name1 (crystal$, negative_graze_angle$,         'NEGATIVE_GRAZE_ANGLE')

call init_attribute_name1 (capillary$, s_spline$,                   'S_SPLINE')
call init_attribute_name1 (capillary$, n_slice_spline$,             'N_SLICE_SPLINE')
call init_attribute_name1 (capillary$, critical_angle_factor$,      'CRITICAL_ANGLE_FACTOR')

!-----------------------------------------------------------------------
! We make a short list to compare against to make things go faster.
! For has_orientation_attributes_key check both tilt and x_offset attributes
! since, for example, a solenoid does not have a tilt.
! Also note: A patch element has a z_offset, not an s_offset.

has_hkick_attributes = .false.  ! Defined in bmad_struct.f90
has_kick_attributes  = .false.  ! Defined in bmad_struct.f90
has_orientation_attributes_key = .false.  ! Defined in bmad_struct.f90

do i = 1, n_key
  if (attrib_array(i, tilt$)%name     == 'TILT')  has_orientation_attributes_key(i) = .true.
  if (attrib_array(i, x_offset$)%name == 'X_OFFSET') has_orientation_attributes_key(i) = .true.
  if (attrib_array(i, kick$)%name     == 'KICK')  has_kick_attributes(i) = .true.
  if (attrib_array(i, hkick$)%name    == 'HKICK') has_hkick_attributes(i) = .true.

  num = 0
  do j = 1, num_ele_attrib_extended$
    if (attrib_array(i, j)%name == null_name$) cycle
    if (attrib_array(i, j)%type == private$) cycle
    num = num + 1
    short_attrib_array(i, num) = attrib_array(i, j)%name
    attrib_ix(i, num) = j
  enddo
  attrib_num(i) = num
enddo

init_needed = .false.

end subroutine init_attribute_name_array

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
  if (bmad_status%exit_on_error) call err_exit
endif

attrib_array(ix_key, ix_attrib)%name    = name
attrib_array(ix_key, ix_attrib)%type = integer_option(is_free$, attrib_state)

end subroutine init_attribute_name1 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function has_orientation_attributes (ele) result (has_attribs)
!
! Routine to determine whether an element has orientation attributes like x_offset, etc.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele  -- ele_struct: Lattice element.
!
! Output:
!   has_attribs -- Logical: True if ele has orientation attributes. False otherwise.
!-

function has_orientation_attributes (ele) result (has_attribs)

implicit none

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
! Routine to return the type (logical, integer, real, or named) of an attribute.
!
! A "Named" attribute is an attribute whose integer value corresponds to some string.
! For example, the "COUPLER_AT" attirbute with value 1 corresponds to "ENTRANCE_END", etc. 
!
! If this is a named attribute, attribute_value_name can be used to print the 
! name corresponding to the attribute's value.
!
! Modules needed:
!   use bmad
!
! Input:
!   attrib_name -- Character(*): Name of the attribute. Must be upper case.
!
! Output:
!   attrib_type  -- Integer: Attribute type: 
!                     is_logical$, is_integer$, is_real$, or is_name$
!-

function attribute_type (attrib_name) result (attrib_type)

implicit none

character(*) attrib_name
integer attrib_type

!

select case (attrib_name)
case ('MATCH_END', 'MATCH_END_ORBIT', 'TRANSLATE_AFTER', 'FOLLOW_DIFFRACTED_BEAM', &
      'NEGATIVE_GRAZE_ANGLE', 'SCALE_MULTIPOLES', 'NO_END_MARKER', 'EXACT_FRINGE', 'INCLUDE_FRINGE')
  attrib_type = is_logical$
case ('TAYLOR_ORDER', 'N_SLICE', 'N_REF_PASS', 'DIRECTION', 'N_CELL', &
      'IX_BRANCH_TO', 'NUM_STEPS', 'INTEGRATOR_ORDER', 'N_LAYERS')
  attrib_type = is_integer$
case ('PARTICLE', 'COUPLER_AT', 'ATTRIBUTE_TYPE', 'REF_POLARAIZATION', 'LATTICE_TYPE')
  attrib_type = is_name$
case default
  attrib_type = is_real$
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
!   c2_curve          c2_curve_tot
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

implicit none

type (ele_struct) ele
integer ix_attrib, ix_tot_attrib

!

select case (attribute_name(ele, ix_attrib))
case ('X_PITCH');     ix_tot_attrib = x_pitch_tot$
case ('Y_PITCH');     ix_tot_attrib = y_pitch_tot$
case ('X_OFFSET');    ix_tot_attrib = x_offset_tot$
case ('Y_OFFSET');    ix_tot_attrib = y_offset_tot$
case ('S_OFFSET');    ix_tot_attrib = s_offset_tot$
case ('TILT');        ix_tot_attrib = tilt_tot$
case ('C2_CURVE');    ix_tot_attrib = c2_curve_tot$
case ('C3_CURVE');    ix_tot_attrib = c3_curve_tot$
case ('C4_CURVE');    ix_tot_attrib = c4_curve_tot$
case default;         ix_tot_attrib = -1
end select

end function corresponding_tot_attribute_index 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function is_a_tot_attribute (ele, ix_attrib) result (is_a_tot_attrib)
!
! Function returns True if ix_attrib corresponds to a "_tot" attribute.
! Example _tot attributes are:
!   "_tot" attribute
!   -----------------
!   x_offset_tot
!   x_pitch_tot
!   tilt_tot
!   c2_curve_tot
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

implicit none

type (ele_struct) ele
integer ix_attrib
logical is_a_tot_attrib

!

select case (attribute_name(ele, ix_attrib))
case ('X_PITCH_TOT', 'Y_PITCH_TOT', 'X_OFFSET_TOT', 'Y_OFFSET_TOT', &
      'S_OFFSET_TOT', 'TILT_TOT', 'C2_CURVE_TOT', 'C3_CURVE_TOT', 'C4_CURVE_TOT')
  is_a_tot_attrib = .true.
case default
  is_a_tot_attrib = .false.
end select

end function is_a_tot_attribute

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_value_name (attrib_name, attrib_value, ele, 
!                                     is_default) result (attrib_val_name)
!
! Routine to return the name corresponding to the value of a given attribute.
!
! This routine is for "switch" attributes. For example, the "aperture_type" attribute
! can have value names of "Entrance_End", "Exit_End", etc.
!
! Optionally, this routine can determine if the attribute value corresponds 
! to the default value. That is, the value that the attribute would have if 
! not specified in the lattice file.
!
! Use the routine attribute_type to first test if the type of the attribute
! corresponds to is_name$. 
!
! Modules needed:
!   use bmad
!
! Input:
!   attrib_name  -- Character(*): Name of the attribute. Must be upper case.
!   attrib_value -- Real(rp): Value of the attribute.
!   ele          -- ele_struct, optional: Lattice element that the attribute is contained in.
!                     This argument is only needed if is_default is present.
!
! Output:
!   attrib_val_name -- Character(40): Name corresponding to the value.
!   is_default      -- Logical, optional: If True then the value of the attiribute
!                        corresponds to the default value. If this argument is
!                        present, the ele argument must also be present.
!-

function attribute_value_name (attrib_name, attrib_value, ele, &
                                     is_default) result (attrib_val_name)

implicit none

type (ele_struct), optional :: ele
character(*) attrib_name
real(rp) attrib_value
character(40) attrib_val_name
character(24) :: r_name = 'attribute_value_name'
logical, optional :: is_default

!

select case (attrib_name)

case ('PARTICLE')
  attrib_val_name = particle_name(nint(attrib_value))
  if (present(is_default)) then
    if (ele%key == photon_branch$) then
      is_default = (nint(attrib_value) == photon$)
    else
      is_default = .false. ! Cannot tell so assume the worst.
    endif
  endif

case ('LATTICE_TYPE')
  attrib_val_name = lattice_type(nint(attrib_value))
  if (present(is_default)) then
    is_default = (nint(attrib_value) == linear_lattice$)
  endif

case ('COUPLER_AT')
  attrib_val_name = end_at_name(nint(attrib_value))
  if (present(is_default)) then
    is_default = (nint(attrib_value) == exit_end$)
  endif

case ('APERTURE_TYPE')
  attrib_val_name = aperture_type_name(nint(attrib_value))
  if (present(is_default)) then
    if (ele%key == ecollimator$) then
      is_default = (nint(attrib_value) == elliptical$)
    else
      is_default = (nint(attrib_value) == rectangular$)
    endif
  endif

case ('REF_POLARIZATION')
  attrib_val_name = polarization_name(nint(attrib_value))
  if (present(is_default)) then
    is_default = (nint(attrib_value) == sigma_polarization$)
  endif


case default
  call out_io (s_fatal$, r_name, 'BAD ATTRIBUTE NAME: ' // attrib_name)
end select

end function attribute_value_name 

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

implicit none
integer max_len
integer, save :: max_length = 0

!

if (init_needed) call init_attribute_name_array
if (max_length == 0) max_length = maxval(len_trim(attrib_array(1:n_key, 1:num_ele_attrib_extended$)%name))
max_len = max_length

end function

end module
