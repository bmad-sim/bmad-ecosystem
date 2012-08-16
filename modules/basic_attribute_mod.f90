module basic_attribute_mod

use bmad_struct

type ele_attribute_struct
  character(40) :: name = null_name$
  logical :: private = .false.
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
  if (attrib_array(key, ix_att)%private) then
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
! Function attribute_record (ele, ix_att) result (attrib_record)
!
! Function to return the record structure associated with an attribute for 
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
!   attrib_record -- ele_attribute_struct: Info on this attribute.
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

function attribute_record (ele, ix_att) result (attrib_record)

implicit none

type (ele_struct) ele
type (ele_attribute_struct) attrib_record
integer i, key, ix_att

!

if (init_needed) call init_attribute_name_array()

attrib_record%private = .false.

if (ele%key <= 0 .or. ele%key > n_key) then
  attrib_record%name = '!BAD ELE KEY'
elseif (ix_att <= 0 .or. ix_att > num_ele_attrib_extended$) then
  attrib_record%name = '!BAD INDEX'
elseif (ele%lord_status == overlay_lord$) then
  if (ix_att == ele%ix_value) then
    attrib_record%name = ele%component_name
  else
    attrib_record%name = '!INVALID INDEX'
  endif
else
  attrib_record = attrib_array(ele%key, ix_att)
endif

end function attribute_record

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

  call init_attrib (i, general1$,  'general1', .true.)
  call init_attrib (i, general2$,  'general2', .true.)
  call init_attrib (i, general3$,  'general3', .true.)
  call init_attrib (i, check_sum$, 'check_sum', .true.)
  call init_attrib (i, scratch$,   'scratch', .true.)

  if (i == crystal$ .or. i == multilayer_mirror$ .or. i == mirror$) then
    call init_attrib (i, d_source$,             'D_SOURCE')
    call init_attrib (i, d_detec$,              'D_DETEC')
    call init_attrib (i, c2_curve$,             'C2_CURVE')
    call init_attrib (i, c3_curve$,             'C3_CURVE')
    call init_attrib (i, c4_curve$,             'C4_CURVE')
    call init_attrib (i, c2_curve_tot$,         'C2_CURVE_TOT')
    call init_attrib (i, c3_curve_tot$,         'C3_CURVE_TOT')
    call init_attrib (i, c4_curve_tot$,         'C4_CURVE_TOT')
  endif

  select case(i)
  case (monitor$, instrument$, marker$, pipe$)
    call init_attrib (i, x_gain_err$,      'X_GAIN_ERR')
    call init_attrib (i, y_gain_err$,      'Y_GAIN_ERR')
    call init_attrib (i, crunch$,          'CRUNCH')
    call init_attrib (i, noise$,           'NOISE')
    call init_attrib (i, tilt_calib$,      'TILT_CALIB')
    call init_attrib (i, x_gain_calib$,    'X_GAIN_CALIB')
    call init_attrib (i, y_gain_calib$,    'Y_GAIN_CALIB')
    call init_attrib (i, crunch_calib$,    'CRUNCH_CALIB')
    call init_attrib (i, x_offset_calib$,  'X_OFFSET_CALIB')
    call init_attrib (i, y_offset_calib$,  'Y_OFFSET_CALIB')
    call init_attrib (i, n_sample$,        'N_SAMPLE')
    call init_attrib (i, de_eta_meas$,     'DE_ETA_MEAS')
    call init_attrib (i, osc_amplitude$,   'OSC_AMPLITUDE')
    call init_attrib (i, x_length$,        'X_LENGTH')
    call init_attrib (i, y_length$,        'Y_LENGTH')
  end select

  if (i == hybrid$)         cycle
  if (i == def_beam$)       cycle
  if (i == def_parameter$)  cycle
  if (i == def_beam_start$) cycle
  if (i == init_ele$) cycle

  call init_attrib (i, type$,      'TYPE')
  call init_attrib (i, alias$,     'ALIAS')
  call init_attrib (i, descrip$,   'DESCRIP')

  if (i == group$)    cycle
  if (i == overlay$)  cycle
  if (i == girder$)   cycle

  call init_attrib (i, superimpose$,        'SUPERIMPOSE')
  call init_attrib (i, offset$,             'OFFSET')
  call init_attrib (i, reference$,          'REFERENCE')
  call init_attrib (i, ele_beginning$,      'ELE_BEGINNING')
  call init_attrib (i, ele_center$,         'ELE_CENTER')
  call init_attrib (i, ele_end$,            'ELE_END')
  call init_attrib (i, ref_beginning$,      'REF_BEGINNING')
  call init_attrib (i, ref_center$,         'REF_CENTER')
  call init_attrib (i, ref_end$,            'REF_END')
  call init_attrib (i, create_em_field_slave$, 'CREATE_EM_FIELD_SLAVE')

  if (i == null_ele$) cycle

  call init_attrib (i, wall_attribute$,         'WALL')
  call init_attrib (i, x_limit$,                'X_LIMIT')
  call init_attrib (i, x1_limit$,               'X1_LIMIT')
  call init_attrib (i, x2_limit$,               'X2_LIMIT')
  call init_attrib (i, y_limit$,                'Y_LIMIT')
  call init_attrib (i, y1_limit$,               'Y1_LIMIT')
  call init_attrib (i, y2_limit$,               'Y2_LIMIT')
  call init_attrib (i, aperture$,               'APERTURE')
  call init_attrib (i, aperture_at$,            'APERTURE_AT')
  call init_attrib (i, aperture_type$,          'APERTURE_TYPE')
  call init_attrib (i, offset_moves_aperture$,  'OFFSET_MOVES_APERTURE')
  call init_attrib (i, mat6_calc_method$,       'MAT6_CALC_METHOD')
  call init_attrib (i, tracking_method$,        'TRACKING_METHOD')
  call init_attrib (i, spin_tracking_method$,   'SPIN_TRACKING_METHOD')

  call init_attrib (i, E_tot$,                  'E_TOT')
  call init_attrib (i, p0c$,                    'P0C')
  call init_attrib (i, delta_ref_time$,         'DELTA_REF_TIME')
  call init_attrib (i, ref_time_start$,         'ref_time_start', .true.)

  if (i == match$) cycle
  if (i == floor_position$) cycle

  call init_attrib (i, tilt$,          'TILT' )
  call init_attrib (i, x_offset$,      'X_OFFSET')
  call init_attrib (i, y_offset$,      'Y_OFFSET')
  call init_attrib (i, s_offset$,      'S_OFFSET')
  call init_attrib (i, x_pitch$,       'X_PITCH')
  call init_attrib (i, y_pitch$,       'Y_PITCH')
  call init_attrib (i, tilt_tot$,      'TILT_TOT')
  call init_attrib (i, x_offset_tot$,  'X_OFFSET_TOT')
  call init_attrib (i, y_offset_tot$,  'Y_OFFSET_TOT')
  call init_attrib (i, s_offset_tot$,  'S_OFFSET_TOT')
  call init_attrib (i, x_pitch_tot$,   'X_PITCH_TOT')
  call init_attrib (i, y_pitch_tot$,   'Y_PITCH_TOT')

  if (i == mirror$)     cycle
  if (i == crystal$)    cycle
  if (i == multilayer_mirror$) cycle
  if (i == capillary$)    cycle

  if (i /= drift$) call init_attrib (i, is_on$,        'IS_ON')

  if (i == photon_branch$) cycle
  if (i == branch$)       cycle
  if (i == marker$)       cycle
  if (i == patch$)        cycle
  if (i == beambeam$)     cycle
  if (i == hom$)          cycle
  if (i == multipole$)    cycle 
  if (i == ab_multipole$) cycle

  call init_attrib (i, symplectify$,         'SYMPLECTIFY')
  call init_attrib (i, map_with_offsets$,    'MAP_WITH_OFFSETS')

  if (i == taylor$)       cycle

  call init_attrib (i, l$,            'L')
  call init_attrib (i, l_hard_edge$,  'L_HARD_EDGE')

  call init_attrib (i, integrator_order$,   'INTEGRATOR_ORDER')
  call init_attrib (i, num_steps$,          'NUM_STEPS')
  call init_attrib (i, ds_step$,            'DS_STEP')
  call init_attrib (i, csr_calc_on$,        'CSR_CALC_ON')
  call init_attrib (i, n_ref_pass$,         'N_REF_PASS')
  call init_attrib (i, field_scale$,        'FIELD_SCALE')

  if (i == hkicker$)      cycle
  if (i == vkicker$)      cycle
  if (i == custom$)       cycle
  if (i == e_gun$)        cycle
  if (i == em_field$)     cycle
  if (i == drift$)        cycle

  call init_attrib (i, hkick$,     'HKICK')
  call init_attrib (i, vkick$,     'VKICK')
  call init_attrib (i, bl_hkick$,  'BL_HKICK')
  call init_attrib (i, bl_vkick$,  'BL_VKICK')

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
    if (i == ab_multipole$) cycle
    call init_attrib (i, scale_multipoles$,  'SCALE_MULTIPOLES')
  end select
enddo

!

call init_attrib (photon_branch$, ix_branch_to$,         'ix_branch_to', .true.)
call init_attrib (photon_branch$, direction$,            'DIRECTION')
call init_attrib (photon_branch$, to$,                   'TO')
call init_attrib (photon_branch$, particle$,             'PARTICLE')
call init_attrib (photon_branch$, lattice_type$,         'LATTICE_TYPE')
call init_attrib (photon_branch$, E_tot_start$,          'E_TOT_START')
call init_attrib (photon_branch$, p0c_start$,            'P0C_START')

attrib_array(branch$, :) = attrib_array(photon_branch$, :)

call init_attrib (init_ele$, delta_ref_time$,              'delta_ref_time', .true.)
call init_attrib (init_ele$, ref_time_start$,              'ref_time_start', .true.)
call init_attrib (init_ele$, e_tot_start$,                 'E_TOT_START')
call init_attrib (init_ele$, p0c_start$,                   'P0C_START')
call init_attrib (init_ele$, e_tot$,                       'E_TOT')
call init_attrib (init_ele$, p0c$,                         'P0C')
call init_attrib (init_ele$, x_position$,                  'X_POSITION')
call init_attrib (init_ele$, y_position$,                  'Y_POSITION')
call init_attrib (init_ele$, z_position$,                  'Z_POSITION')
call init_attrib (init_ele$, theta_position$,              'THETA_POSITION')
call init_attrib (init_ele$, phi_position$,                'PHI_POSITION')
call init_attrib (init_ele$, psi_position$,                'PSI_POSITION')
call init_attrib (init_ele$, beta_a$,                      'BETA_A')
call init_attrib (init_ele$, beta_b$,                      'BETA_B')
call init_attrib (init_ele$, alpha_a$,                     'ALPHA_A')
call init_attrib (init_ele$, alpha_b$,                     'ALPHA_B')
call init_attrib (init_ele$, eta_x$,                       'ETA_X')
call init_attrib (init_ele$, eta_y$,                       'ETA_Y')
call init_attrib (init_ele$, eta_z$,                       'ETA_Z')
call init_attrib (init_ele$, etap_x$,                      'ETAP_X')
call init_attrib (init_ele$, etap_y$,                      'ETAP_Y')
call init_attrib (init_ele$, phi_a$,                       'PHI_A')
call init_attrib (init_ele$, phi_b$,                       'PHI_B')
call init_attrib (init_ele$, cmat_11$,                     'CMAT_11')
call init_attrib (init_ele$, cmat_12$,                     'CMAT_12')
call init_attrib (init_ele$, cmat_21$,                     'CMAT_21')
call init_attrib (init_ele$, cmat_22$,                     'CMAT_22')
call init_attrib (init_ele$, s_long$,                      'S')
call init_attrib (init_ele$, ref_time$,                    'REF_TIME')
call init_attrib (init_ele$, e_field_x$,                   'E_FIELD_X')
call init_attrib (init_ele$, e_field_y$,                   'E_FIELD_Y')
call init_attrib (init_ele$, phase_x$,                     'PHASE_X')
call init_attrib (init_ele$, phase_y$,                     'PHASE_Y')
call init_attrib (init_ele$, wall_attribute$,              'WALL')

call init_attrib (def_parameter$, e_tot$,                  'E_TOT')
call init_attrib (def_parameter$, p0c$,                    'P0C')
call init_attrib (def_parameter$, lattice_type$,           'LATTICE_TYPE')
call init_attrib (def_parameter$, lattice$,                'LATTICE')
call init_attrib (def_parameter$, taylor_order$,           'TAYLOR_ORDER')
call init_attrib (def_parameter$, ran_seed$,               'RAN_SEED')
call init_attrib (def_parameter$, n_part$,                 'N_PART')
call init_attrib (def_parameter$, particle$,               'PARTICLE')
call init_attrib (def_parameter$, aperture_limit_on$,      'APERTURE_LIMIT_ON')
call init_attrib (def_parameter$, root_branch_name$,       'ROOT_BRANCH_NAME') 
call init_attrib (def_parameter$, no_end_marker$,          'NO_END_MARKER')
call init_attrib (def_parameter$, absolute_time_tracking$, 'ABSOLUTE_TIME_TRACKING')
call init_attrib (def_parameter$, use_ptc_layout$,         'USE_PTC_LAYOUT')
call init_attrib (def_parameter$, rf_auto_scale_phase$,    'RF_AUTO_SCALE_PHASE')
call init_attrib (def_parameter$, rf_auto_scale_amp$,      'RF_AUTO_SCALE_AMP')
call init_attrib (def_parameter$, ptc_exact_model$,        'PTC_EXACT_MODEL')
call init_attrib (def_parameter$, ptc_exact_misalign$,     'PTC_EXACT_MISALIGN')

call init_attrib (def_beam$, particle$,                    'PARTICLE')
call init_attrib (def_beam$, e_tot$,                       'ENERGY')
call init_attrib (def_beam$, p0c$,                         'PC')
call init_attrib (def_beam$, n_part$,                      'N_PART')

call init_attrib (def_beam_start$, x$,                     'X')
call init_attrib (def_beam_start$, px$,                    'PX')
call init_attrib (def_beam_start$, y$,                     'Y')
call init_attrib (def_beam_start$, py$,                    'PY')
call init_attrib (def_beam_start$, z$,                     'Z')
call init_attrib (def_beam_start$, pz$,                    'PZ')
call init_attrib (def_beam_start$, e_field_x$,             'E_FIELD_X')
call init_attrib (def_beam_start$, e_field_y$,             'E_FIELD_Y')
call init_attrib (def_beam_start$, phase_x$,               'PHASE_X')
call init_attrib (def_beam_start$, phase_y$,               'PHASE_Y')
call init_attrib (def_beam_start$, t$,                     'T')

call init_attrib (e_gun$, e_tot_start$,                    'e_tot_start', .true.)
call init_attrib (e_gun$, p0c_start$,                      'p0c_start', .true.)
call init_attrib (e_gun$, voltage$,                        'VOLTAGE')
call init_attrib (e_gun$, voltage_err$,                    'VOLTAGE_ERR')
call init_attrib (e_gun$, gradient$,                       'GRADIENT')
call init_attrib (e_gun$, gradient_err$,                   'GRADIENT_ERR')
call init_attrib (e_gun$, field_calc$,                     'FIELD_CALC')
call init_attrib (e_gun$, field$,                          'FIELD')

call init_attrib (em_field$, e_tot_start$,                 'E_TOT_START')
call init_attrib (em_field$, p0c_start$,                   'P0C_START')
call init_attrib (em_field$, field$,                       'FIELD')

call init_attrib (taylor$, l$,                             'L')
call init_attrib (taylor$, E_tot_start$,                   'E_tot_start', .true.)
call init_attrib (taylor$, p0c_start$,                     'p0c_start', .true.)

call init_attrib (marker$, l$,                             'l', .true.)
call init_attrib (marker$, E_tot_start$,                   'E_tot_start', .true.)
call init_attrib (marker$, p0c_start$,                     'p0c_start', .true.)

call init_attrib (match$, l$,                              'L')
call init_attrib (match$, beta_a0$,                        'BETA_A0')
call init_attrib (match$, alpha_a0$,                       'ALPHA_A0')
call init_attrib (match$, beta_b0$,                        'BETA_B0')
call init_attrib (match$, alpha_b0$,                       'ALPHA_B0')
call init_attrib (match$, beta_a1$,                        'BETA_A1')
call init_attrib (match$, alpha_a1$,                       'ALPHA_A1')
call init_attrib (match$, beta_b1$,                        'BETA_B1')
call init_attrib (match$, alpha_b1$,                       'ALPHA_B1')
call init_attrib (match$, dphi_a$,                         'DPHI_A')
call init_attrib (match$, dphi_b$,                         'DPHI_B')
call init_attrib (match$, eta_x0$,                         'ETA_X0')
call init_attrib (match$, etap_x0$,                        'ETAP_X0')
call init_attrib (match$, eta_y0$,                         'ETA_Y0')
call init_attrib (match$, etap_y0$,                        'ETAP_Y0')
call init_attrib (match$, eta_x1$,                         'ETA_X1')
call init_attrib (match$, etap_x1$,                        'ETAP_X1')
call init_attrib (match$, eta_y1$,                         'ETA_Y1')
call init_attrib (match$, etap_y1$,                        'ETAP_Y1')
call init_attrib (match$, match_end$,                      'MATCH_END')
call init_attrib (match$, x0$,                             'X0')
call init_attrib (match$, px0$,                            'PX0')
call init_attrib (match$, y0$,                             'Y0')
call init_attrib (match$, py0$,                            'PY0')
call init_attrib (match$, z0$,                             'Z0')
call init_attrib (match$, pz0$,                            'PZ0')
call init_attrib (match$, x1$,                             'X1')
call init_attrib (match$, px1$,                            'PX1')
call init_attrib (match$, y1$,                             'Y1')
call init_attrib (match$, py1$,                            'PY1')
call init_attrib (match$, z1$,                             'Z1')
call init_attrib (match$, pz1$,                            'PZ1')
call init_attrib (match$, match_end_orbit$,                'MATCH_END_ORBIT')
call init_attrib (match$, is_on$,                          'IS_ON')
call init_attrib (match$, E_tot_start$,                   'E_tot_start', .true.)
call init_attrib (match$, p0c_start$,                     'p0c_start', .true.)

call init_attrib (girder$, l$,                             'L')
call init_attrib (girder$, x_offset$,                      'X_OFFSET')
call init_attrib (girder$, y_offset$,                      'Y_OFFSET')
call init_attrib (girder$, s_offset$,                      'S_OFFSET')
call init_attrib (girder$, x_pitch$,                       'X_PITCH')
call init_attrib (girder$, y_pitch$,                       'Y_PITCH')
call init_attrib (girder$, s_min$,                         'S_MIN')
call init_attrib (girder$, s_max$,                         'S_MAX')
call init_attrib (girder$, s_center$,                      'S_CENTER')
call init_attrib (girder$, tilt$,                          'TILT')

call init_attrib (lcavity$, p0c_start$,                    'P0C_START')
call init_attrib (lcavity$, e_tot_start$,                  'E_TOT_START')
call init_attrib (lcavity$, dphi0$,                        'DPHI0')
call init_attrib (lcavity$, phi0$,                         'PHI0')
call init_attrib (lcavity$, gradient$,                     'GRADIENT')
call init_attrib (lcavity$, rf_frequency$,                 'RF_FREQUENCY')
call init_attrib (lcavity$, e_loss$,                       'E_LOSS')
call init_attrib (lcavity$, voltage$,                      'VOLTAGE')
call init_attrib (lcavity$, sr_wake_file$,                 'SR_WAKE_FILE')
call init_attrib (lcavity$, lr_wake_file$,                 'LR_WAKE_FILE')
call init_attrib (lcavity$, field_calc$,                   'FIELD_CALC')
call init_attrib (lcavity$, field_master$,                 'FIELD_MASTER')
call init_attrib (lcavity$, lr_freq_spread$,               'LR_FREQ_SPREAD')
call init_attrib (lcavity$, coupler_strength$,             'COUPLER_STRENGTH')
call init_attrib (lcavity$, coupler_angle$,                'COUPLER_ANGLE')
call init_attrib (lcavity$, coupler_phase$,                'COUPLER_PHASE')
call init_attrib (lcavity$, coupler_at$,                   'COUPLER_AT')
call init_attrib (lcavity$, gradient_err$,                 'GRADIENT_ERR')
call init_attrib (lcavity$, phi0_err$,                     'PHI0_ERR')
call init_attrib (lcavity$, field$,                        'FIELD')
call init_attrib (lcavity$, dphi0_ref$,                    'DPHI0_REF')
call init_attrib (lcavity$, n_cell$,                       'N_CELL')
call init_attrib (lcavity$, grad_loss_sr_wake$,            'grad_loss_sr_wake', .true.)


call init_attrib (group$, command$,                        'COMMAND')
call init_attrib (group$, old_command$,                    'OLD_COMMAND')
call init_attrib (group$, coef$,                           'COEF')
call init_attrib (group$, start_edge$,                     'START_EDGE')
call init_attrib (group$, end_edge$,                       'END_EDGE')
call init_attrib (group$, accordion_edge$,                 'ACCORDION_EDGE')
call init_attrib (group$, symmetric_edge$,                 'SYMMETRIC_EDGE')

call init_attrib (drift$, field_calc$,                     'FIELD_CALC')
call init_attrib (drift$, field_master$,                   'FIELD_MASTER')
call init_attrib (drift$, E_tot_start$,                   'E_tot_start', .true.)
call init_attrib (drift$, p0c_start$,                     'p0c_start', .true.)

call init_attrib (monitor$, field_master$,                 'FIELD_MASTER')
call init_attrib (monitor$, E_tot_start$,                   'E_tot_start', .true.)
call init_attrib (monitor$, p0c_start$,                     'p0c_start', .true.)

attrib_array(instrument$, :)                         = attrib_array(monitor$, :)
attrib_array(pipe$, :)                               = attrib_array(monitor$, :)

call init_attrib (hkicker$, kick$,                         'KICK')
call init_attrib (hkicker$, field_calc$,                   'FIELD_CALC')
call init_attrib (hkicker$, field_master$,                 'FIELD_MASTER')
call init_attrib (hkicker$, bl_kick$,                      'BL_KICK')
call init_attrib (hkicker$, pole_radius$,                  'POLE_RADIUS')
call init_attrib (hkicker$, E_tot_start$,                   'E_tot_start', .true.)
call init_attrib (hkicker$, p0c_start$,                     'p0c_start', .true.)

attrib_array(vkicker$, :) = attrib_array(hkicker$, :)

call init_attrib (kicker$, h_displace$,                    'H_DISPLACE')
call init_attrib (kicker$, v_displace$,                    'V_DISPLACE')
call init_attrib (kicker$, radius$,                        'RADIUS')
call init_attrib (kicker$, field_calc$,                    'FIELD_CALC')
call init_attrib (kicker$, field_master$,                  'FIELD_MASTER')
call init_attrib (kicker$, pole_radius$,                   'POLE_RADIUS')
call init_attrib (kicker$, E_tot_start$,                   'E_tot_start', .true.)
call init_attrib (kicker$, p0c_start$,                     'p0c_start', .true.)

call init_attrib (sbend$, angle$,                          'ANGLE')
call init_attrib (sbend$, exact_fringe$,                   'EXACT_FRINGE')
call init_attrib (sbend$, e1$,                             'E1')
call init_attrib (sbend$, e2$,                             'E2')
call init_attrib (sbend$, h1$,                             'H1')
call init_attrib (sbend$, h2$,                             'H2')
call init_attrib (sbend$, k1$,                             'K1')
call init_attrib (sbend$, k2$,                             'K2')
call init_attrib (sbend$, g$,                              'G')
call init_attrib (sbend$, g_err$,                          'G_ERR')
call init_attrib (sbend$, roll$,                           'ROLL')
call init_attrib (sbend$, hgap$,                           'HGAP')
call init_attrib (sbend$, hgapx$,                          'HGAPX')
call init_attrib (sbend$, fint$,                           'FINT')
call init_attrib (sbend$, fintx$,                          'FINTX')
call init_attrib (sbend$, rho$,                            'RHO')
call init_attrib (sbend$, l_chord$,                        'L_CHORD')
call init_attrib (sbend$, b_field$,                        'B_FIELD')
call init_attrib (sbend$, b_field_err$,                    'B_FIELD_ERR')
call init_attrib (sbend$, b1_gradient$,                    'B1_GRADIENT')
call init_attrib (sbend$, b2_gradient$,                    'B2_GRADIENT')
call init_attrib (sbend$, radius$,                         'RADIUS')
call init_attrib (sbend$, field_calc$,                     'FIELD_CALC')
call init_attrib (sbend$, field_master$,                   'FIELD_MASTER')
call init_attrib (sbend$, ref_orbit$,                      'REF_ORBIT')
call init_attrib (sbend$, E_tot_start$,                   'E_tot_start', .true.)
call init_attrib (sbend$, p0c_start$,                     'p0c_start', .true.)

attrib_array(rbend$, :) = attrib_array(sbend$, :)

call init_attrib (bend_sol_quad$, angle$,                  'ANGLE')
call init_attrib (bend_sol_quad$, k1$,                     'K1')
call init_attrib (bend_sol_quad$, g$,                      'G')
call init_attrib (bend_sol_quad$, ks$,                     'KS')
call init_attrib (bend_sol_quad$, dks_ds$,                 'DKS_DS')
call init_attrib (bend_sol_quad$, quad_tilt$,              'QUAD_TILT')
call init_attrib (bend_sol_quad$, bend_tilt$,              'BEND_TILT')
call init_attrib (bend_sol_quad$, x_quad$,                 'X_QUAD')
call init_attrib (bend_sol_quad$, y_quad$,                 'Y_QUAD')
call init_attrib (bend_sol_quad$, rho$,                    'RHO')
call init_attrib (bend_sol_quad$, radius$,                 'RADIUS')
call init_attrib (bend_sol_quad$, field_calc$,             'FIELD_CALC')
call init_attrib (bend_sol_quad$, field_master$,           'FIELD_MASTER')
call init_attrib (bend_sol_quad$, field$,                  'FIELD')
call init_attrib (bend_sol_quad$, E_tot_start$,            'E_tot_start', .true.)
call init_attrib (bend_sol_quad$, p0c_start$,              'p0c_start', .true.)

call init_attrib (patch$, l$,                              'L')
call init_attrib (patch$, t_offset$,                       'T_OFFSET')
call init_attrib (patch$, p0c_start$,                      'P0C_START')
call init_attrib (patch$, e_tot_start$,                    'E_TOT_START')
call init_attrib (patch$, e_tot_offset$,                   'E_TOT_OFFSET')
call init_attrib (patch$, ref_orbit$,                      'REF_ORBIT')
call init_attrib (patch$, ref_patch$,                      'REF_PATCH')
call init_attrib (patch$, n_ref_pass$,                     'N_REF_PASS')
call init_attrib (patch$, translate_after$,                'TRANSLATE_AFTER')
call init_attrib (patch$, z_offset$,                       'Z_OFFSET', override = .true.)

call init_attrib (floor_position$, l$,                           'l', .true.)
call init_attrib (floor_position$, x_position$,                  'X_POSITION')
call init_attrib (floor_position$, y_position$,                  'Y_POSITION')
call init_attrib (floor_position$, z_position$,                  'Z_POSITION')
call init_attrib (floor_position$, theta_position$,              'THETA_POSITION')
call init_attrib (floor_position$, phi_position$,                'PHI_POSITION')
call init_attrib (floor_position$, psi_position$,                'PSI_POSITION')


call init_attrib (quadrupole$, k1$,                        'K1')
call init_attrib (quadrupole$, B1_gradient$,               'B1_GRADIENT')
call init_attrib (quadrupole$, radius$,                    'RADIUS')
call init_attrib (quadrupole$, field_calc$,                'FIELD_CALC')
call init_attrib (quadrupole$, field_master$,              'FIELD_MASTER')
call init_attrib (quadrupole$, pole_radius$,               'POLE_RADIUS')
call init_attrib (quadrupole$, field$,                     'FIELD')
call init_attrib (quadrupole$, E_tot_start$,               'E_tot_start', .true.)
call init_attrib (quadrupole$, p0c_start$,                 'p0c_start', .true.)

call init_attrib (sextupole$, k2$,                         'K2')
call init_attrib (sextupole$, B2_gradient$,                'B2_GRADIENT')
call init_attrib (sextupole$, radius$,                     'RADIUS')
call init_attrib (sextupole$, field_calc$,                 'FIELD_CALC')
call init_attrib (sextupole$, field_master$,               'FIELD_MASTER')
call init_attrib (sextupole$, pole_radius$,                'POLE_RADIUS')
call init_attrib (sextupole$, field$,                      'FIELD')
call init_attrib (sextupole$, E_tot_start$,                'E_tot_start', .true.)
call init_attrib (sextupole$, p0c_start$,                  'p0c_start', .true.)

call init_attrib (octupole$, k3$,                          'K3')
call init_attrib (octupole$, B3_gradient$,                 'B3_GRADIENT')
call init_attrib (octupole$, radius$,                      'RADIUS')
call init_attrib (octupole$, field_calc$,                  'FIELD_CALC')
call init_attrib (octupole$, field_master$,                'FIELD_MASTER')
call init_attrib (octupole$, pole_radius$,                 'POLE_RADIUS')
call init_attrib (octupole$, field$,                       'FIELD')
call init_attrib (octupole$, E_tot_start$,                 'E_tot_start', .true.)
call init_attrib (octupole$, p0c_start$,                   'p0c_start', .true.)

call init_attrib (solenoid$, ks$,                          'KS')
call init_attrib (solenoid$, bs_field$,                    'BS_FIELD')
call init_attrib (solenoid$, radius$,                      'RADIUS')
call init_attrib (solenoid$, field_calc$,                  'FIELD_CALC')
call init_attrib (solenoid$, field_master$,                'FIELD_MASTER')
call init_attrib (solenoid$, pole_radius$,                 'POLE_RADIUS')
call init_attrib (solenoid$, field$,                       'FIELD')
call init_attrib (solenoid$, E_tot_start$,                 'E_tot_start', .true.)
call init_attrib (solenoid$, p0c_start$,                   'p0c_start', .true.)

call init_attrib (rfcavity$, dphi0$,                       'DPHI0')
call init_attrib (rfcavity$, voltage$,                     'VOLTAGE')
call init_attrib (rfcavity$, rf_frequency$,                'RF_FREQUENCY')
call init_attrib (rfcavity$, phi0$,                        'PHI0')
call init_attrib (rfcavity$, phi0_err$,                    'phi0_err', .true.)
call init_attrib (rfcavity$, harmon$,                      'HARMON')
call init_attrib (rfcavity$, field_calc$,                  'FIELD_CALC')
call init_attrib (rfcavity$, field_master$,                'FIELD_MASTER')
call init_attrib (rfcavity$, sr_wake_file$,                'SR_WAKE_FILE')
call init_attrib (rfcavity$, lr_wake_file$,                'LR_WAKE_FILE')
call init_attrib (rfcavity$, lr_freq_spread$,              'LR_FREQ_SPREAD')
call init_attrib (rfcavity$, coupler_strength$,            'COUPLER_STRENGTH')
call init_attrib (rfcavity$, coupler_angle$,               'COUPLER_ANGLE')
call init_attrib (rfcavity$, coupler_phase$,               'COUPLER_PHASE')
call init_attrib (rfcavity$, coupler_at$,                  'COUPLER_AT')
call init_attrib (rfcavity$, field$,                       'FIELD')
call init_attrib (rfcavity$, dphi0_ref$,                   'DPHI0_REF')
call init_attrib (rfcavity$, n_cell$,                      'N_CELL')
call init_attrib (rfcavity$, grad_loss_sr_wake$,           'grad_loss_sr_wake', .true.)
call init_attrib (rfcavity$, dphi0_max$,                   'dphi0_max', .true.)
call init_attrib (rfcavity$, E_tot_start$,                   'E_tot_start', .true.)
call init_attrib (rfcavity$, p0c_start$,                     'p0c_start', .true.)

call init_attrib (elseparator$, gap$,                      'GAP')
call init_attrib (elseparator$, e_field$,                  'E_FIELD')
call init_attrib (elseparator$, voltage$,                  'VOLTAGE')
call init_attrib (elseparator$, radius$,                   'RADIUS')
call init_attrib (elseparator$, field_calc$,               'FIELD_CALC')
call init_attrib (elseparator$, field_master$,             'FIELD_MASTER')
call init_attrib (elseparator$, field$,                    'FIELD')
call init_attrib (elseparator$, E_tot_start$,                   'E_tot_start', .true.)
call init_attrib (elseparator$, p0c_start$,                     'p0c_start', .true.)

call init_attrib (beambeam$, l$,                           'l', .true.)
call init_attrib (beambeam$, sig_x$,                       'SIG_X')
call init_attrib (beambeam$, sig_y$,                       'SIG_Y')
call init_attrib (beambeam$, sig_z$,                       'SIG_Z')
call init_attrib (beambeam$, bbi_const$,                   'BBI_CONSTANT')
call init_attrib (beambeam$, charge$,                      'CHARGE')
call init_attrib (beambeam$, n_slice$,                     'N_SLICE')
call init_attrib (beambeam$, symplectify$,                 'N_SLICE')
call init_attrib (beambeam$, field_calc$,                  'FIELD_CALC')
call init_attrib (beambeam$, field_master$,                'FIELD_MASTER')
call init_attrib (beambeam$, E_tot_start$,                   'E_tot_start', .true.)
call init_attrib (beambeam$, p0c_start$,                     'p0c_start', .true.)

call init_attrib (wiggler$, k1$,                           'K1')
call init_attrib (wiggler$, l_pole$,                       'L_POLE')
call init_attrib (wiggler$, b_max$,                        'B_MAX')
call init_attrib (wiggler$, rho$,                          'RHO')
call init_attrib (wiggler$, n_pole$,                       'N_POLE')
call init_attrib (wiggler$, term$,                         'TERM')
call init_attrib (wiggler$, polarity$,                     'POLARITY')
call init_attrib (wiggler$, radius$,                       'RADIUS')
call init_attrib (wiggler$, field_calc$,                   'FIELD_CALC')
call init_attrib (wiggler$, field_master$,                 'FIELD_MASTER')
call init_attrib (wiggler$, x_ray_line_len$,               'X_RAY_LINE_LEN')
call init_attrib (wiggler$, field$,                        'FIELD')
call init_attrib (wiggler$, E_tot_start$,                   'E_tot_start', .true.)
call init_attrib (wiggler$, p0c_start$,                     'p0c_start', .true.)

call init_attrib (sol_quad$, k1$,                          'K1')
call init_attrib (sol_quad$, ks$,                          'KS')
call init_attrib (sol_quad$, radius$,                      'RADIUS')
call init_attrib (sol_quad$, field_calc$,                  'FIELD_CALC')
call init_attrib (sol_quad$, field_master$,                'FIELD_MASTER')
call init_attrib (sol_quad$, b1_gradient$,                 'B1_GRADIENT')
call init_attrib (sol_quad$, bs_field$,                    'BS_FIELD')
call init_attrib (sol_quad$, field$,                       'FIELD')
call init_attrib (sol_quad$, E_tot_start$,                 'E_tot_start', .true.)
call init_attrib (sol_quad$, p0c_start$,                   'p0c_start', .true.)

call init_attrib (multipole$, l$,                          'L')
attrib_array(multipole$, k0l$:k20l$)%name    = &
             ['K0L ', 'K1L ', 'K2L ', 'K3L ', 'K4L ', 'K5L ', 'K6L ', 'K7L ', 'K8L ', 'K9L ', 'K10L', &
                      'K11L', 'K12L', 'K13L', 'K14L', 'K15L', 'K16L', 'K17L', 'K18L', 'K19L', 'K20L']
attrib_array(multipole$, t0$:t20$)%name = ['T0 ', &
                               'T1 ', 'T2 ', 'T3 ', 'T4 ', 'T5 ', & 
                               'T6 ', 'T7 ', 'T8 ', 'T9 ', 'T10', &
                               'T11', 'T12', 'T13', 'T14', 'T15', &
                               'T16', 'T17', 'T18', 'T19', 'T20' ]

call init_attrib (ab_multipole$, l$,                       'L')

call init_attrib (custom$, val1$,                          'VAL1')
call init_attrib (custom$, val2$,                          'VAL2')
call init_attrib (custom$, val3$,                          'VAL3')
call init_attrib (custom$, val4$,                          'VAL4')
call init_attrib (custom$, val5$,                          'VAL5')
call init_attrib (custom$, val6$,                          'VAL6')
call init_attrib (custom$, val7$,                          'VAL7')
call init_attrib (custom$, val8$,                          'VAL8')
call init_attrib (custom$, val9$,                          'VAL9')
call init_attrib (custom$, val10$,                         'VAL10')
call init_attrib (custom$, val11$,                         'VAL11')
call init_attrib (custom$, val12$,                         'VAL12')
call init_attrib (custom$, field_calc$,                    'FIELD_CALC')
call init_attrib (custom$, field_master$,                  'FIELD_MASTER')
call init_attrib (custom$, delta_e$,                       'DELTA_E')
call init_attrib (custom$, e_tot_start$,                   'E_TOT_START')
call init_attrib (custom$, p0c_start$,                     'P0C_START')

call init_attrib (hybrid$, l$,                             'L')
call init_attrib (hybrid$, delta_e$,                       'DELTA_E')
call init_attrib (hybrid$, delta_ref_time$,                'DELTA_REF_TIME')
call init_attrib (hybrid$, ref_time_start$,                'ref_time_start', .true.)
call init_attrib (hybrid$, e_tot_start$,                   'E_TOT_START')
call init_attrib (hybrid$, p0c_start$,                     'P0C_START')

call init_attrib (mirror$, l$,                             'l', .true.)
call init_attrib (mirror$, graze_angle$,                   'GRAZE_ANGLE')
call init_attrib (mirror$, graze_angle_err$,               'GRAZE_ANGLE_ERR')
call init_attrib (mirror$, critical_angle$,                'CRITICAL_ANGLE')
call init_attrib (mirror$, tilt_err$,                      'TILT_ERR')
call init_attrib (mirror$, g_trans$,                       'G_TRANS')
call init_attrib (mirror$, ref_wavelength$,                'REF_WAVELENGTH')

call init_attrib (multilayer_mirror$, l$,                    'l', .true.)
call init_attrib (multilayer_mirror$, graze_angle$,          'GRAZE_ANGLE')
call init_attrib (multilayer_mirror$, graze_angle_err$,      'GRAZE_ANGLE_ERR')
call init_attrib (multilayer_mirror$, tilt_err$,             'TILT_ERR')
call init_attrib (multilayer_mirror$, n_cell$,               'N_CELL')
call init_attrib (multilayer_mirror$, d1_thickness$,         'D1_THICKNESS')
call init_attrib (multilayer_mirror$, d2_thickness$,         'D2_THICKNESS')
call init_attrib (multilayer_mirror$, v1_unitcell$,          'V1_UNITCELL')
call init_attrib (multilayer_mirror$, v2_unitcell$,          'V2_UNITCELL')
call init_attrib (multilayer_mirror$, f0_re1$,               'F0_RE1')
call init_attrib (multilayer_mirror$, f0_im1$,               'F0_IM1')
call init_attrib (multilayer_mirror$, f0_re2$,               'F0_RE2')
call init_attrib (multilayer_mirror$, f0_im2$,               'F0_IM2')
call init_attrib (multilayer_mirror$, ref_wavelength$,       'REF_WAVELENGTH')
call init_attrib (multilayer_mirror$, crystal_type$,         'CRYSTAL_TYPE')
call init_attrib (multilayer_mirror$, ref_polarization$,     'REF_POLARIZATION')
call init_attrib (multilayer_mirror$, negative_graze_angle$, 'NEGATIVE_GRAZE_ANGLE')

call init_attrib (crystal$, l$,                            'l', .true.)
call init_attrib (crystal$, graze_angle_in$,               'GRAZE_ANGLE_IN')
call init_attrib (crystal$, graze_angle_out$,              'GRAZE_ANGLE_OUT')
call init_attrib (crystal$, graze_angle_err$,              'GRAZE_ANGLE_ERR')
call init_attrib (crystal$, psi_angle$,                    'PSI_ANGLE')
call init_attrib (crystal$, alpha_angle$,                  'ALPHA_ANGLE')
call init_attrib (crystal$, tilt_err$,                     'TILT_ERR')
call init_attrib (crystal$, tilt_corr$,                    'TILT_CORR')
call init_attrib (crystal$, d_spacing$,                    'D_SPACING')
call init_attrib (crystal$, v_unitcell$,                   'V_UNITCELL')
call init_attrib (crystal$, f0_re$,                        'F0_RE')
call init_attrib (crystal$, f0_im$,                        'F0_IM')
call init_attrib (crystal$, fh_re$,                        'FH_RE')
call init_attrib (crystal$, fh_im$,                        'FH_IM')
call init_attrib (crystal$, b_param$,                      'B_PARAM')
call init_attrib (crystal$, bragg_angle$,                  'BRAGG_ANGLE' )
call init_attrib (crystal$, g_trans$,                      'G_TRANS')
call init_attrib (crystal$, ref_wavelength$,               'REF_WAVELENGTH')
call init_attrib (crystal$, diffraction_type$,             'DIFFRACTION_TYPE')  ! Note: Not currently used.
call init_attrib (crystal$, crystal_type$,                 'CRYSTAL_TYPE')
call init_attrib (crystal$, thickness$,                    'THICKNESS')
call init_attrib (crystal$, follow_diffracted_beam$,       'FOLLOW_DIFFRACTED_BEAM')
call init_attrib (crystal$, kh_x_norm$,                    'kh_x_norm', .true.)
call init_attrib (crystal$, kh_y_norm$,                    'kh_y_norm', .true.)
call init_attrib (crystal$, kh_z_norm$,                    'kh_z_norm', .true.)
call init_attrib (crystal$, l_x$,                          'l_x', .true.)
call init_attrib (crystal$, l_y$,                          'l_y', .true.)
call init_attrib (crystal$, l_z$,                          'l_z', .true.)
call init_attrib (crystal$, ref_polarization$,             'REF_POLARIZATION')
call init_attrib (crystal$, ref_cap_gamma$,                'REF_CAP_GAMMA')
call init_attrib (crystal$, negative_graze_angle$,         'NEGATIVE_GRAZE_ANGLE')

call init_attrib (capillary$, s_spline$,                   'S_SPLINE')
call init_attrib (capillary$, n_slice_spline$,             'N_SLICE_SPLINE')
call init_attrib (capillary$, critical_angle_factor$,      'CRITICAL_ANGLE_FACTOR')

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
    if (attrib_array(i, j)%private) cycle
    num = num + 1
    short_attrib_array(i, num) = attrib_array(i, j)%name
    attrib_ix(i, num) = j
  enddo
  attrib_num(i) = num
enddo

init_needed = .false.

!-------------------------------------------------------------------------
contains

subroutine init_attrib (ix_key, ix_attrib, name, private, override)

integer ix_key, ix_attrib
character(*) name
logical, optional :: private, override

! Check that attrib_array(ix_key, ix_attrib)%name has not already been set.
! If so bomb program.

if (.not. logic_option(.false., override) .and. attrib_array(ix_key, ix_attrib)%name /= null_name$) then
  call out_io (s_fatal$, 'init_attrib', 'ERROR IN INITIALIZING ATTRIB_ARRAY FOR: ' // key_name(ix_key), &
                  'IX_ATTRIB \i0\ ALREADY SET!', &
                  'OLD/NEW NAMES: ' // trim(attrib_array(ix_key, ix_attrib)%name) // ' : ' // name, &
                  i_array = [ix_attrib])
  if (bmad_status%exit_on_error) call err_exit
endif

attrib_array(ix_key, ix_attrib)%name    = name
attrib_array(ix_key, ix_attrib)%private = logic_option(.false., private)

end subroutine init_attrib 

end subroutine init_attribute_name_array

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
      'NEGATIVE_GRAZE_ANGLE', 'SCALE_MULTIPOLES', 'NO_END_MARKER', 'EXACT_FRINGE')
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
