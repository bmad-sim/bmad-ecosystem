module attribute_mod

use bmad_struct
use basic_bmad_interface
use multipole_mod


character(40), private, save :: attrib_array(n_key, n_attrib_special_maxx)
character(40), private, save :: short_attrib_array(n_key, n_attrib_special_maxx)
integer, private, save :: attrib_num(n_key)
integer, private, save :: attrib_ix(n_key, n_attrib_special_maxx)
logical, private, save :: init_needed = .true.

type(ele_struct), private, pointer, save :: ele0 ! For Error message purposes
character(40), private, save :: attrib_name0     ! For Error message purposes

private init_attribute_name_array, check_this_attribute_free, print_error


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free
!
! Overloaded function for:
!   Function attribute_free1 (ix_ele, attrib_name, lat,
!                                err_print_flag, except_overlay) result (free)
!   Function attribute_free2 (ix_ele, ix_branch, attrib_name, lat, 
!                                err_print_flag, except_overlay) result (free)
!   Function attribute_free3 (loc, attrib_name, lat, 
!                                err_print_flag, except_overlay) result (free)
!
! Subroutine to check if an attribute is free to vary.
!
! Attributes that cannot be changed directly include super_slave attributes (since
! these attributes are controlled by their super_lords) and attributes that
! are controlled by an overlay_lord.
!
! Also dependent variables such as the angle of a bend cannot be 
!   freely variable.
!
! Modules needed:
!   use bmad
!
! Input:
!   loc             -- Lat_ele_loc_struct: Element location.
!   ix_ele          -- Integer: Index of element in lat%ele(:) array.
!   ix_branch       -- Integer: Branch index. Default is 0.
!   attrib_name     -- Character(*): Name of the attribute. Assumed upper case.
!   lat             -- lat_struct: Lattice structure.
!   err_print_flag  -- Logical, optional: If present and False then supress
!                       printing of an error message if attribute is not free.
!   except_overlay  -- Logical, optional: If present and True then an attribute that
!                       is an overlay_slave will be treated as free. This is used by,
!                       for example, the create_overlay routine.
!
! Output:
!   free   -- Logical: Set True if attribtute not found or attriubte
!                     cannot be changed directly.
!-

interface attribute_free
  module procedure attribute_free1
  module procedure attribute_free2
  module procedure attribute_free3
end interface

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine pointer_to_indexed_attribute (ele, ix_attrib, do_allocation,
!                                      ptr_attrib, err_flag, err_print_flag)
!
! Returns a pointer to an attribute of an element ele with attribute index ix_attrib.
! 
! Use of this routine is restricted to attributes that have an index. That is,
! attributes in the ele%value(:) array and ele%a_pole(:), and ele%b_pole(:) values.
! A more general routine is pointer_to_attribute.
! Alternatively, consider the routine pointers_to_attribute.
! Note: Use attribute_free to see if the attribute may be varied independently.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele             -- Ele_struct: After this routine finishes Ptr_attrib 
!                        will point to a variable within this element.
!   ix_attrib       -- Integer, Attribute index.
!   do_allocation   -- Logical: If True then do an allocation if needed.
!                       EG: The multipole An and Bn arrays need to be allocated
!                       before their use.
!   err_print_flag  -- Logical, optional: If present and False then supress
!                       printing of an error message on error.
!
! Output:
!   ptr_attrib -- Real(rp), pointer: Pointer to the attribute.
!                     Pointer will be deassociated if there is a problem.
!   err_flag   -- Logical: Set True if attribtute not found. False otherwise.
!-

subroutine pointer_to_indexed_attribute (ele, ix_attrib, do_allocation, &
                                        ptr_attrib, err_flag, err_print_flag)

implicit none

type (ele_struct), target :: ele

real(rp), pointer :: ptr_attrib

integer :: ix_attrib

character(40) :: r_name = 'pointer_to_indexed_attribute'

logical err_flag, do_allocation, do_print
logical, optional :: err_print_flag

! Init

err_flag = .true.
nullify (ptr_attrib)
do_print = logic_option (.true., err_print_flag)

! multipole?

if (ix_attrib >= a0$ .and. ix_attrib <= b20$) then   ! multipole attribute

  if (.not. associated(ele%a_pole)) then
    if (do_allocation) then
      call multipole_init (ele)
    else
      if (do_print) call out_io (s_error$, r_name, &
                      'MULTIPOLE NOT ALLOCATED FOR ELEMENT: ' // ele%name)
      return
    endif
  endif

  if (ix_attrib >= b0$) then
    ptr_attrib => ele%b_pole(ix_attrib-b0$)
  else
    ptr_attrib => ele%a_pole(ix_attrib-a0$)
  endif

elseif (ix_attrib < 1 .or. ix_attrib > n_attrib_maxx) then
  if (do_print) call out_io (s_error$, r_name, &
          'INVALID ATTRIBUTE INDEX: \i0\ ', 'FOR THIS ELEMENT: ' // ele%name, &
          i_array = (/ ix_attrib /))
  return

! otherwise must be in ele%value(:) array

else
  ptr_attrib => ele%value(ix_attrib)
endif


err_flag = .false.
return

end subroutine pointer_to_indexed_attribute 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+             
! Function attribute_index (ele, name) result (attrib_index)
!
! Function to return the index of a attribute for a given BMAD element type
! and the name of the attribute.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele  -- Ele_struct: attribute_index will restrict the name search to 
!             valid attributes of the given element. Note: If 
!             ele%key = overlay then the entire name table will be searched.
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

name40 = name          ! make sure we have 40 characters
key = ele%key
attrib_index = 0           ! match not found

ilen = len_trim(name)
if (ilen == 0) return
n_abbrev = 0            ! number of abbreviation matches

!-----------------------------------------------------------------------
! Support old "B_GRADIENT" notation

if (name == 'B_GRADIENT') then
  select case (key)
  case (quadrupole$, sextupole$, octupole$, sol_quad$, overlay$, group$)
    attrib_index = b_gradient$
  end select
  return
endif

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
  call err_exit
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
!              If %key is invalid then                   attribute_name = "!BAD ELE KEY"
!              If ix_att is invalid then                 attribute_name = "!BAD INDEX"
!              If ix_att is invalid for an overlay then  attribute_name = "!INVALID INDEX"
!              If ix_att does not correspond to an attribute for the given key then
!                                                        attribute_name = "!NULL"
!
! Example:
!   ele%key = sbend$
!   name = attribute_name (ele, k1$)
! Result:
!   name -> 'K1'
!-

function attribute_name (ele, ix_att) result (at_name)

implicit none

type (ele_struct) ele

integer i, key, ix_att

character(40) at_name

!--------------------------------------------------------------------

if (init_needed) call init_attribute_name_array()

key = ele%key

if (key <= 0 .or. key > n_key) then
  at_name = '!BAD ELE KEY'
elseif (ix_att <= 0 .or. ix_att > n_attrib_special_maxx) then
  at_name = '!BAD INDEX'
elseif (ele%lord_status == overlay_lord$) then
  if (ix_att == ele%ix_value) then
    at_name = ele%attribute_name
  else
    at_name = '!INVALID INDEX'
  endif
else
  at_name = attrib_array(key, ix_att)
endif

end function attribute_name 

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

attrib_array = null_name$

do i = 1, n_key
                                
  if (i == monitor$ .or. i == instrument$ .or. i == marker$) then 
    attrib_array(i, x_gain_err$)     = 'X_GAIN_ERR'
    attrib_array(i, y_gain_err$)     = 'Y_GAIN_ERR'
    attrib_array(i, crunch$)         = 'CRUNCH'
    attrib_array(i, noise$)          = 'NOISE'
    attrib_array(i, tilt_calib$)     = 'TILT_CALIB'
    attrib_array(i, x_gain_calib$)   = 'X_GAIN_CALIB'
    attrib_array(i, y_gain_calib$)   = 'Y_GAIN_CALIB'
    attrib_array(i, crunch_calib$)   = 'CRUNCH_CALIB'
    attrib_array(i, x_offset_calib$) = 'X_OFFSET_CALIB'
    attrib_array(i, y_offset_calib$) = 'Y_OFFSET_CALIB'
    attrib_array(i, y_offset_calib$) = 'Y_OFFSET_CALIB'
    attrib_array(i, n_sample$)       = 'N_SAMPLE'
    attrib_array(i, de_eta_meas$)    = 'DE_ETA_MEAS'
    attrib_array(i, osc_amplitude$)  = 'OSC_AMPLITUDE'
  endif

  if (i == hybrid$)         cycle
  if (i == def_beam$)       cycle
  if (i == def_parameter$)  cycle
  if (i == def_beam_start$) cycle
  if (i == init_ele$) cycle

  attrib_array(i, type$)     = 'TYPE'
  attrib_array(i, alias$)    = 'ALIAS'
  attrib_array(i, descrip$)  = 'DESCRIP'

  if (i == group$)    cycle
  if (i == overlay$)  cycle
  if (i == girder$)   cycle
  if (i == mirror$)   cycle
  if (i == crystal$)  cycle

  attrib_array(i, superimpose$)       = 'SUPERIMPOSE'
  attrib_array(i, offset$)            = 'OFFSET'
  attrib_array(i, reference$)         = 'REFERENCE'
  attrib_array(i, ele_beginning$)     = 'ELE_BEGINNING'
  attrib_array(i, ele_center$)        = 'ELE_CENTER'
  attrib_array(i, ele_end$)           = 'ELE_END'
  attrib_array(i, ref_beginning$)     = 'REF_BEGINNING'
  attrib_array(i, ref_center$)        = 'REF_CENTER'
  attrib_array(i, ref_end$)           = 'REF_END'
  attrib_array(i, common_lord$)       = 'COMMON_LORD'

  if (i == photon_branch$) cycle
  if (i == branch$) cycle

  attrib_array(i, E_tot$)                 = 'E_TOT'
  attrib_array(i, p0c$)                   = 'P0C'
  attrib_array(i, x_limit$)               = 'X_LIMIT'
  attrib_array(i, x1_limit$)              = 'X1_LIMIT'
  attrib_array(i, x2_limit$)              = 'X2_LIMIT'
  attrib_array(i, y_limit$)               = 'Y_LIMIT'
  attrib_array(i, y1_limit$)              = 'Y1_LIMIT'
  attrib_array(i, y2_limit$)              = 'Y2_LIMIT'
  attrib_array(i, aperture$)              = 'APERTURE'
  attrib_array(i, aperture_at$)           = 'APERTURE_AT'
  attrib_array(i, offset_moves_aperture$) = 'OFFSET_MOVES_APERTURE'
  attrib_array(i, mat6_calc_method$)      = 'MAT6_CALC_METHOD'
  attrib_array(i, tracking_method$)       = 'TRACKING_METHOD'

  attrib_array(i, is_on$)       = 'IS_ON'

  if (i == marker$)       cycle
  if (i == match$)        cycle
  if (i == patch$)        cycle
  if (i == beambeam$)     cycle
  if (i == hom$)          cycle
  if (i == multipole$)    cycle 
  if (i == ab_multipole$) cycle

  attrib_array(i, symplectify$)        = 'SYMPLECTIFY'
  attrib_array(i, map_with_offsets$)   = 'MAP_WITH_OFFSETS'

  if (i == taylor$)       cycle

  attrib_array(i, integrator_order$)  = 'INTEGRATOR_ORDER'
  attrib_array(i, num_steps$)         = 'NUM_STEPS'
  attrib_array(i, ds_step$)           = 'DS_STEP'
  attrib_array(i, csr_calc_on$)       = 'CSR_CALC_ON'
  attrib_array(i, n_ref_pass$)        = 'N_REF_PASS'

  if (i == hkicker$)      cycle
  if (i == vkicker$)      cycle
  if (i == custom$)       cycle

  attrib_array(i, hkick$)    = 'HKICK'
  attrib_array(i, vkick$)    = 'VKICK'
  attrib_array(i, bl_hkick$) = 'BL_HKICK'
  attrib_array(i, bl_vkick$) = 'BL_VKICK'

  if (i == drift$)        cycle
  if (i == kicker$)       cycle
  if (i == monitor$)      cycle
  if (i == instrument$)   cycle

  attrib_array(i, x_offset$) = 'X_OFFSET'
  attrib_array(i, y_offset$) = 'Y_OFFSET'
  attrib_array(i, s_offset$) = 'S_OFFSET'

  if (i == rcollimator$) cycle
  if (i == ecollimator$) cycle

  attrib_array(i, x_pitch$)   = 'X_PITCH'
  attrib_array(i, y_pitch$)   = 'Y_PITCH'

enddo

!

do i = 1, n_key
  select case (i)
  case (elseparator$, kicker$, octupole$, quadrupole$, sbend$, rbend$, &
         sextupole$, solenoid$, sol_quad$, ab_multipole$, wiggler$, bend_sol_quad$, &
         hkicker$, vkicker$)
    attrib_array(i, a0$:a20$) = (/ 'A0 ', &
                                   'A1 ', 'A2 ', 'A3 ', 'A4 ', 'A5 ', & 
                                   'A6 ', 'A7 ', 'A8 ', 'A9 ', 'A10', &
                                   'A11', 'A12', 'A13', 'A14', 'A15', &
                                   'A16', 'A17', 'A18', 'A19', 'A20' /)
    attrib_array(i, b0$:b20$) = (/ 'B0 ', &
                                   'B1 ', 'B2 ', 'B3 ', 'B4 ', 'B5 ', & 
                                   'B6 ', 'B7 ', 'B8 ', 'B9 ', 'B10', &
                                   'B11', 'B12', 'B13', 'B14', 'B15', &
                                   'B16', 'B17', 'B18', 'B19', 'B20' /)
  end select
enddo

!

attrib_array(photon_branch$, direction$) = 'DIRECTION'
attrib_array(photon_branch$, to$)        = 'TO'

attrib_array(branch$, :) = attrib_array(photon_branch$, :)

attrib_array(init_ele$, E_tot$)       = 'E_TOT'
attrib_array(init_ele$, p0c$)         = 'P0C'

attrib_array(def_parameter$, E_TOT$)              = 'E_TOT'
attrib_array(def_parameter$, lattice_type$)       = 'LATTICE_TYPE'
attrib_array(def_parameter$, lattice$)            = 'LATTICE'
attrib_array(def_parameter$, taylor_order$)       = 'TAYLOR_ORDER'
attrib_array(def_parameter$, ran_seed$)           = 'RAN_SEED'
attrib_array(def_parameter$, n_part$)             = 'N_PART'
attrib_array(def_parameter$, particle$)           = 'PARTICLE'
attrib_array(def_parameter$, aperture_limit_on$)  = 'APERTURE_LIMIT_ON'

attrib_array(def_beam$, particle$)   = 'PARTICLE'
attrib_array(def_beam$, e_tot$)      = 'ENERGY'
attrib_array(def_beam$, p0c$)        = 'PC'
attrib_array(def_beam$, n_part$)     = 'N_PART'

attrib_array(def_beam_start$, x$)     = 'X'
attrib_array(def_beam_start$, px$)    = 'PX'
attrib_array(def_beam_start$, y$)     = 'Y'
attrib_array(def_beam_start$, py$)    = 'PY'
attrib_array(def_beam_start$, z$)     = 'Z'
attrib_array(def_beam_start$, pz$)    = 'PZ'

attrib_array(taylor$, l$)           = 'L'
attrib_array(taylor$, x_offset$)    = 'X_OFFSET'   
attrib_array(taylor$, y_offset$)    = 'Y_OFFSET'   
attrib_array(taylor$, s_offset$)    = 'S_OFFSET'   
attrib_array(taylor$, x_pitch$)     = 'X_PITCH'   
attrib_array(taylor$, y_pitch$)     = 'Y_PITCH'   
attrib_array(taylor$, tilt$)        = 'TILT' 

attrib_array(match$, l$)               = 'L'
attrib_array(match$, beta_a0$)         = 'BETA_A0'
attrib_array(match$, alpha_a0$)        = 'ALPHA_A0'
attrib_array(match$, beta_b0$)         = 'BETA_B0'
attrib_array(match$, alpha_b0$)        = 'ALPHA_B0'
attrib_array(match$, beta_a1$)         = 'BETA_A1'
attrib_array(match$, alpha_a1$)        = 'ALPHA_A1'
attrib_array(match$, beta_b1$)         = 'BETA_B1'
attrib_array(match$, alpha_b1$)        = 'ALPHA_B1'
attrib_array(match$, dphi_a$)          = 'DPHI_A'
attrib_array(match$, dphi_b$)          = 'DPHI_B'
attrib_array(match$, eta_x0$)          = 'ETA_X0'
attrib_array(match$, etap_x0$)         = 'ETAP_X0'
attrib_array(match$, eta_y0$)          = 'ETA_Y0'
attrib_array(match$, etap_y0$)         = 'ETAP_Y0'
attrib_array(match$, eta_x1$)          = 'ETA_X1'
attrib_array(match$, etap_x1$)         = 'ETAP_X1'
attrib_array(match$, eta_y1$)          = 'ETA_Y1'
attrib_array(match$, etap_y1$)         = 'ETAP_Y1'
attrib_array(match$, match_end$)       = 'MATCH_END'
attrib_array(match$, x0$)              = 'X0'
attrib_array(match$, px0$)             = 'PX0'
attrib_array(match$, y0$)              = 'Y0'
attrib_array(match$, py0$)             = 'PY0'
attrib_array(match$, z0$)              = 'Z0'
attrib_array(match$, pz0$)             = 'PZ0'
attrib_array(match$, x1$)              = 'X1'
attrib_array(match$, px1$)             = 'PX1'
attrib_array(match$, y1$)              = 'Y1'
attrib_array(match$, py1$)             = 'PY1'
attrib_array(match$, z1$)              = 'Z1'
attrib_array(match$, pz1$)             = 'PZ1'
attrib_array(match$, match_end_orbit$) = 'MATCH_END_ORBIT'

attrib_array(girder$, x_offset$)     = 'X_OFFSET'
attrib_array(girder$, y_offset$)     = 'Y_OFFSET'
attrib_array(girder$, s_offset$)     = 'S_OFFSET'
attrib_array(girder$, x_pitch$)      = 'X_PITCH'
attrib_array(girder$, y_pitch$)      = 'Y_PITCH'
attrib_array(girder$, s_center$)     = 'S_CENTER'
attrib_array(girder$, tilt$)         = 'TILT'

attrib_array(lcavity$, l$)                = 'L'
attrib_array(lcavity$, tilt$)             = 'TILT'
attrib_array(lcavity$, lrad$)             = 'LRAD'   ! This is for felv testing.
attrib_array(lcavity$, p0c_start$)        = 'P0C_START'
attrib_array(lcavity$, e_tot_start$)      = 'E_TOT_START'
attrib_array(lcavity$, dphi0$)            = 'DPHI0'
attrib_array(lcavity$, phi0$)             = 'PHI0'
attrib_array(lcavity$, gradient$)         = 'GRADIENT'
attrib_array(lcavity$, rf_frequency$)     = 'RF_FREQUENCY'
attrib_array(lcavity$, e_loss$)           = 'E_LOSS'
attrib_array(lcavity$, delta_e$)          = 'DELTA_E'
attrib_array(lcavity$, sr_wake_file$)     = 'SR_WAKE_FILE'
attrib_array(lcavity$, lr_wake_file$)     = 'LR_WAKE_FILE'
attrib_array(lcavity$, field_calc$)       = 'FIELD_CALC'
attrib_array(lcavity$, field_master$)     = 'FIELD_MASTER'
attrib_array(lcavity$, lr_freq_spread$)   = 'LR_FREQ_SPREAD'
attrib_array(lcavity$, coupler_at$)       = 'COUPLER_AT'
attrib_array(lcavity$, coupler_strength$) = 'COUPLER_STRENGTH'
attrib_array(lcavity$, coupler_angle$)    = 'COUPLER_ANGLE'
attrib_array(lcavity$, coupler_phase$)    = 'COUPLER_PHASE'
attrib_array(lcavity$, gradient_err$)     = 'GRADIENT_ERR'
attrib_array(lcavity$, phi0_err$)         = 'PHI0_ERR'


attrib_array(group$, command$)        = 'COMMAND'
attrib_array(group$, old_command$)    = 'OLD_COMMAND'
attrib_array(group$, coef$)           = 'COEF'
attrib_array(group$, start_edge$)     = 'START_EDGE'
attrib_array(group$, end_edge$)       = 'END_EDGE'
attrib_array(group$, accordion_edge$) = 'ACCORDION_EDGE'
attrib_array(group$, symmetric_edge$) = 'SYMMETRIC_EDGE'

attrib_array(drift$, l$)             = 'L'
attrib_array(drift$, is_on$)         =  null_name$    
attrib_array(drift$, field_calc$)    = 'FIELD_CALC'
attrib_array(drift$, field_master$)  = 'FIELD_MASTER'

attrib_array(monitor$, l$)           = 'L'
attrib_array(monitor$, x_offset$)    = 'X_OFFSET'
attrib_array(monitor$, y_offset$)    = 'Y_OFFSET'
attrib_array(monitor$, x_pitch$)     = 'X_PITCH'
attrib_array(monitor$, y_pitch$)     = 'Y_PITCH'
attrib_array(monitor$, tilt$)        = 'TILT'

attrib_array(instrument$, :) = attrib_array(monitor$, :)

attrib_array(marker$, x_offset$) = 'X_OFFSET'
attrib_array(marker$, y_offset$) = 'Y_OFFSET'
attrib_array(marker$, tilt$)     = 'TILT'

attrib_array(rcollimator$, l$)     = 'L'

attrib_array(ecollimator$, l$)     = 'L'

attrib_array(hkicker$, l$)            = 'L'
attrib_array(hkicker$, tilt$)         = 'TILT'
attrib_array(hkicker$, kick$)         = 'KICK'
attrib_array(hkicker$, field_calc$)   = 'FIELD_CALC'
attrib_array(hkicker$, field_master$) = 'FIELD_MASTER'
attrib_array(hkicker$, bl_kick$)      = 'BL_KICK'
attrib_array(hkicker$, s_offset$)     = 'S_OFFSET'

attrib_array(vkicker$, :) = attrib_array(hkicker$, :)

attrib_array(kicker$, l$)            = 'L'
attrib_array(kicker$, tilt$)         = 'TILT'
attrib_array(kicker$, h_displace$)   = 'H_DISPLACE'
attrib_array(kicker$, v_displace$)   = 'V_DISPLACE'
attrib_array(kicker$, radius$)       = 'RADIUS'
attrib_array(kicker$, field_calc$)   = 'FIELD_CALC'
attrib_array(kicker$, field_master$) = 'FIELD_MASTER'
attrib_array(kicker$, s_offset$)     = 'S_OFFSET'

attrib_array(sbend$, l$)                  = 'L'
attrib_array(sbend$, angle$)              = 'ANGLE'
attrib_array(sbend$, e1$)                 = 'E1'
attrib_array(sbend$, e2$)                 = 'E2'
attrib_array(sbend$, h1$)                 = 'H1'
attrib_array(sbend$, h2$)                 = 'H2'
attrib_array(sbend$, k1$)                 = 'K1'
attrib_array(sbend$, k2$)                 = 'K2'
attrib_array(sbend$, g$)                  = 'G'
attrib_array(sbend$, g_err$)              = 'G_ERR'
attrib_array(sbend$, tilt$)               = 'TILT'
attrib_array(sbend$, roll$)               = 'ROLL'
attrib_array(sbend$, hgap$)               = 'HGAP'
attrib_array(sbend$, hgapx$)              = 'HGAPX'
attrib_array(sbend$, fint$)               = 'FINT'
attrib_array(sbend$, fintx$)              = 'FINTX'
attrib_array(sbend$, rho$)                = 'RHO'
attrib_array(sbend$, l_chord$)            = 'L_CHORD'
attrib_array(sbend$, b_field$)            = 'B_FIELD'
attrib_array(sbend$, b_field_err$)        = 'B_FIELD_ERR'
attrib_array(sbend$, b1_gradient$)        = 'B1_GRADIENT'
attrib_array(sbend$, b2_gradient$)        = 'B2_GRADIENT'
attrib_array(sbend$, radius$)             = 'RADIUS'
attrib_array(sbend$, field_calc$)         = 'FIELD_CALC'
attrib_array(sbend$, field_master$)       = 'FIELD_MASTER'
attrib_array(sbend$, ref_orbit$)          = 'REF_ORBIT'

attrib_array(rbend$, :) = attrib_array(sbend$, :)

attrib_array(bend_sol_quad$, l$)            = 'L'
attrib_array(bend_sol_quad$, angle$)        = 'ANGLE'
attrib_array(bend_sol_quad$, k1$)           = 'K1'
attrib_array(bend_sol_quad$, g$)            = 'G'
attrib_array(bend_sol_quad$, ks$)           = 'KS'
attrib_array(bend_sol_quad$, dks_ds$)       = 'DKS_DS'
attrib_array(bend_sol_quad$, quad_tilt$)    = 'QUAD_TILT'
attrib_array(bend_sol_quad$, bend_tilt$)    = 'BEND_TILT'
attrib_array(bend_sol_quad$, x_quad$)       = 'X_QUAD'
attrib_array(bend_sol_quad$, y_quad$)       = 'Y_QUAD'
attrib_array(bend_sol_quad$, tilt$)         = 'TILT'
attrib_array(bend_sol_quad$, rho$)          = 'RHO'
attrib_array(bend_sol_quad$, radius$)       = 'RADIUS'
attrib_array(bend_sol_quad$, field_calc$)   = 'FIELD_CALC'
attrib_array(bend_sol_quad$, field_master$) = 'FIELD_MASTER'

attrib_array(patch$, l$)          = 'L'
attrib_array(patch$, x_pitch$)    = 'X_PITCH'
attrib_array(patch$, y_pitch$)    = 'Y_PITCH'
attrib_array(patch$, x_offset$)   = 'X_OFFSET'
attrib_array(patch$, y_offset$)   = 'Y_OFFSET'
attrib_array(patch$, z_offset$)   = 'Z_OFFSET'
attrib_array(patch$, pz_offset$)  = 'PZ_OFFSET'
attrib_array(patch$, tilt$)       = 'TILT'

attrib_array(quadrupole$, l$)             = 'L'
attrib_array(quadrupole$, tilt$)          = 'TILT'
attrib_array(quadrupole$, k1$)            = 'K1'
attrib_array(quadrupole$, B1_gradient$)   = 'B1_GRADIENT'
attrib_array(quadrupole$, radius$)        = 'RADIUS'
attrib_array(quadrupole$, field_calc$)    = 'FIELD_CALC'
attrib_array(quadrupole$, field_master$)  = 'FIELD_MASTER'

attrib_array(sextupole$, l$)             = 'L'
attrib_array(sextupole$, tilt$)          = 'TILT'
attrib_array(sextupole$, k2$)            = 'K2'
attrib_array(sextupole$, B2_gradient$)   = 'B2_GRADIENT'
attrib_array(sextupole$, radius$)        = 'RADIUS'
attrib_array(sextupole$, field_calc$)    = 'FIELD_CALC'
attrib_array(sextupole$, field_master$)  = 'FIELD_MASTER'

attrib_array(octupole$, l$)             = 'L'
attrib_array(octupole$, tilt$)          = 'TILT'
attrib_array(octupole$, k3$)            = 'K3'
attrib_array(octupole$, B3_gradient$)   = 'B3_GRADIENT'
attrib_array(octupole$, radius$)        = 'RADIUS'
attrib_array(octupole$, field_calc$)    = 'FIELD_CALC'
attrib_array(octupole$, field_master$)  = 'FIELD_MASTER'

attrib_array(solenoid$, l$)            = 'L'
attrib_array(solenoid$, ks$)           = 'KS'
attrib_array(solenoid$, bs_field$)     = 'BS_FIELD'
attrib_array(solenoid$, radius$)       = 'RADIUS'
attrib_array(solenoid$, field_calc$)   = 'FIELD_CALC'
attrib_array(solenoid$, field_master$) = 'FIELD_MASTER'

attrib_array(rfcavity$, l$)             = 'L'
attrib_array(rfcavity$, dphi0$)         = 'DPHI0'
attrib_array(rfcavity$, voltage$)       = 'VOLTAGE'
attrib_array(rfcavity$, rf_frequency$)  = 'RF_FREQUENCY'
attrib_array(rfcavity$, phi0$)          = 'PHI0'
attrib_array(rfcavity$, harmon$)        = 'HARMON'
attrib_array(rfcavity$, field_calc$)    = 'FIELD_CALC'
attrib_array(rfcavity$, field_master$)  = 'FIELD_MASTER'

attrib_array(elseparator$, l$)            = 'L'
attrib_array(elseparator$, gap$)          = 'GAP'
attrib_array(elseparator$, e_field$)      = 'E_FIELD'
attrib_array(elseparator$, voltage$)      = 'VOLTAGE'
attrib_array(elseparator$, tilt$)         = 'TILT'
attrib_array(elseparator$, radius$)       = 'RADIUS'
attrib_array(elseparator$, field_calc$)   = 'FIELD_CALC'
attrib_array(elseparator$, field_master$) = 'FIELD_MASTER'

attrib_array(beambeam$, sig_x$)         = 'SIG_X'
attrib_array(beambeam$, sig_y$)         = 'SIG_Y'
attrib_array(beambeam$, sig_z$)         = 'SIG_Z'
attrib_array(beambeam$, bbi_const$)     = 'BBI_CONSTANT'
attrib_array(beambeam$, charge$)        = 'CHARGE'
attrib_array(beambeam$, n_slice$)       = 'N_SLICE'
attrib_array(beambeam$, symplectify$)   = 'N_SLICE'
attrib_array(beambeam$, x_offset$)      = 'X_OFFSET'
attrib_array(beambeam$, y_offset$)      = 'Y_OFFSET'
attrib_array(beambeam$, s_offset$)      = 'S_OFFSET'
attrib_array(beambeam$, x_pitch$)       = 'X_PITCH'
attrib_array(beambeam$, y_pitch$)       = 'Y_PITCH'
attrib_array(beambeam$, tilt$)          = 'TILT'
attrib_array(beambeam$, field_calc$)    = 'FIELD_CALC'
attrib_array(beambeam$, field_master$)  = 'FIELD_MASTER'

attrib_array(wiggler$, l$)              = 'L'
attrib_array(wiggler$, k1$)             = 'K1'
attrib_array(wiggler$, l_pole$)         = 'L_POLE'
attrib_array(wiggler$, b_max$)          = 'B_MAX'
attrib_array(wiggler$, rho$)            = 'RHO'
attrib_array(wiggler$, n_pole$)         = 'N_POLE'
attrib_array(wiggler$, tilt$)           = 'TILT'
attrib_array(wiggler$, radius$)         = 'RADIUS'
attrib_array(wiggler$, term$)           = 'TERM'
attrib_array(wiggler$, polarity$)       = 'POLARITY'
attrib_array(wiggler$, z_patch$)        = 'Z_PATCH'
attrib_array(wiggler$, radius$)         = 'RADIUS'
attrib_array(wiggler$, field_calc$)     = 'FIELD_CALC'
attrib_array(wiggler$, field_master$)   = 'FIELD_MASTER'
attrib_array(wiggler$, x_ray_line_len$) = 'X_RAY_LINE_LEN'
attrib_array(wiggler$, l_start$)        = 'L_START'
attrib_array(wiggler$, l_end$)          = 'L_END'

attrib_array(sol_quad$, l$)             = 'L'
attrib_array(sol_quad$, k1$)            = 'K1'
attrib_array(sol_quad$, ks$)            = 'KS'
attrib_array(sol_quad$, tilt$)          = 'TILT'
attrib_array(sol_quad$, radius$)        = 'RADIUS'
attrib_array(sol_quad$, field_calc$)    = 'FIELD_CALC'
attrib_array(sol_quad$, field_master$)  = 'FIELD_MASTER'
attrib_array(sol_quad$, b1_gradient$)   = 'B1_GRADIENT'
attrib_array(sol_quad$, bs_field$)      = 'BS_FIELD'


attrib_array(multipole$, l$)         = 'L'
attrib_array(multipole$, tilt$)      = 'TILT'
attrib_array(multipole$, k0l$:k20l$) = (/ 'K0L ', &
                               'K1L ', 'K2L ', 'K3L ', 'K4L ', 'K5L ', & 
                               'K6L ', 'K7L ', 'K8L ', 'K9L ', 'K10L', &
                               'K11L', 'K12L', 'K13L', 'K14L', 'K15L', &
                               'K16L', 'K17L', 'K18L', 'K19L', 'K20L' /)
attrib_array(multipole$, t0$:t20$) = (/ 'T0 ', &
                               'T1 ', 'T2 ', 'T3 ', 'T4 ', 'T5 ', & 
                               'T6 ', 'T7 ', 'T8 ', 'T9 ', 'T10', &
                               'T11', 'T12', 'T13', 'T14', 'T15', &
                               'T16', 'T17', 'T18', 'T19', 'T20' /)
attrib_array(multipole$, x_offset$) = 'X_OFFSET'
attrib_array(multipole$, y_offset$) = 'Y_OFFSET'
attrib_array(multipole$, s_offset$) = 'S_OFFSET'

attrib_array(ab_multipole$, l$)        = 'L'
attrib_array(ab_multipole$, tilt$)     = 'TILT'
attrib_array(ab_multipole$, x_offset$) = 'X_OFFSET'
attrib_array(ab_multipole$, y_offset$) = 'Y_OFFSET'
attrib_array(ab_multipole$, s_offset$) = 'S_OFFSET'

attrib_array(accel_sol$, l$)             = 'L'
attrib_array(accel_sol$, voltage$)       = 'VOLTAGE'
attrib_array(accel_sol$, phi0$)          = 'PHI0'
attrib_array(accel_sol$, rf_wavelength$) = 'RF_WAVELENGTH'
attrib_array(accel_sol$, b_z$)           = 'B_Z'
attrib_array(accel_sol$, b_x1$)          = 'B_X1'
attrib_array(accel_sol$, b_y1$)          = 'B_Y1'
attrib_array(accel_sol$, s_st1$)         = 'S_ST1'
attrib_array(accel_sol$, l_st1$)         = 'L_ST1'
attrib_array(accel_sol$, b_x2$)          = 'B_X2'
attrib_array(accel_sol$, b_y2$)          = 'B_Y2'
attrib_array(accel_sol$, s_st2$)         = 'S_ST2'
attrib_array(accel_sol$, l_st2$)         = 'L_ST2'
attrib_array(accel_sol$, x_beg_limit$)   = 'X_BEG_LIMIT'
attrib_array(accel_sol$, y_beg_limit$)   = 'Y_BEG_LIMIT'
attrib_array(accel_sol$, field_calc$)    = 'FIELD_CALC'
attrib_array(accel_sol$, field_master$)  = 'FIELD_MASTER'

attrib_array(custom$, l$)            = 'L'
attrib_array(custom$, tilt$)         = 'TILT'
attrib_array(custom$,  val1$)        = 'VAL1'
attrib_array(custom$,  val2$)        = 'VAL2'
attrib_array(custom$,  val3$)        = 'VAL3'
attrib_array(custom$,  val4$)        = 'VAL4'
attrib_array(custom$,  val5$)        = 'VAL5'
attrib_array(custom$,  val6$)        = 'VAL6'
attrib_array(custom$,  val7$)        = 'VAL7'
attrib_array(custom$,  val8$)        = 'VAL8'
attrib_array(custom$,  val9$)        = 'VAL9'
attrib_array(custom$, val10$)        = 'VAL10'
attrib_array(custom$, val11$)        = 'VAL11'
attrib_array(custom$, val12$)        = 'VAL12'
attrib_array(custom$, x_offset$)     = 'X_OFFSET'
attrib_array(custom$, y_offset$)     = 'Y_OFFSET'
attrib_array(custom$, s_offset$)     = 'S_OFFSET'
attrib_array(custom$, x_pitch$)      = 'X_PITCH'
attrib_array(custom$, y_pitch$)      = 'Y_PITCH'
attrib_array(custom$, field_calc$)   = 'FIELD_CALC'
attrib_array(custom$, field_master$) = 'FIELD_MASTER'
attrib_array(custom$, delta_e$)      = 'DELTA_E'

attrib_array(hybrid$, l$)              = 'L'
attrib_array(hybrid$, delta_e$)        = 'DELTA_E'
attrib_array(hybrid$, delta_ref_time$) = 'DELTA_REF_TIME'
attrib_array(hybrid$, e_tot_start$)    = 'E_TOT_START'

attrib_array(mirror$, graze_angle$)     = 'GRAZE_ANGLE'
attrib_array(mirror$, graze_angle_err$) = 'GRAZE_ANGLE_ERR'
attrib_array(mirror$, critical_angle$)  = 'CRITICAL_ANGLE'
attrib_array(mirror$, tilt$)            = 'TILT'
attrib_array(mirror$, tilt_err$)        = 'TILT_ERR'
attrib_array(mirror$, x_offset$)        = 'X_OFFSET'   
attrib_array(mirror$, y_offset$)        = 'Y_OFFSET'   
attrib_array(mirror$, s_offset$)        = 'S_OFFSET'   
attrib_array(mirror$, x_pitch$)         = 'X_PITCH'
attrib_array(mirror$, y_pitch$)         = 'Y_PITCH'
attrib_array(mirror$, g_graze$)         = 'G_GRAZE'
attrib_array(mirror$, g_trans$)         = 'G_TRANS'

attrib_array(init_ele$, x_position$)      = 'X_POSITION'
attrib_array(init_ele$, y_position$)      = 'Y_POSITION'
attrib_array(init_ele$, z_position$)      = 'Z_POSITION'
attrib_array(init_ele$, theta_position$)  = 'THETA_POSITION'
attrib_array(init_ele$, phi_position$)    = 'PHI_POSITION'
attrib_array(init_ele$, psi_position$)    = 'PSI_POSITION'
attrib_array(init_ele$, beta_a$)          = 'BETA_A'
attrib_array(init_ele$, beta_b$)          = 'BETA_B'
attrib_array(init_ele$, alpha_a$)         = 'ALPHA_A'
attrib_array(init_ele$, alpha_b$)         = 'ALPHA_B'
attrib_array(init_ele$, eta_x$)           = 'ETA_X'
attrib_array(init_ele$, eta_y$)           = 'ETA_Y'
attrib_array(init_ele$, etap_x$)          = 'ETAP_X'
attrib_array(init_ele$, etap_y$)          = 'ETAP_Y'
attrib_array(init_ele$, phi_a$)           = 'PHI_A'
attrib_array(init_ele$, phi_b$)           = 'PHI_B'
attrib_array(init_ele$, cmat_11$)         = 'CMAT_11'
attrib_array(init_ele$, cmat_12$)         = 'CMAT_12'
attrib_array(init_ele$, cmat_21$)         = 'CMAT_21'
attrib_array(init_ele$, cmat_22$)         = 'CMAT_22'
attrib_array(init_ele$, s_long$)          = 'S'
attrib_array(init_ele$, ref_time$)        = 'REF_TIME'

!-----------------------------------------------------------------------
! We make a short list to compare against to make things go faster

do i = 1, n_key
  num = 0
  do j = 1, n_attrib_special_maxx
    if (attrib_array(i, j) /= null_name$) then
      num = num + 1
      short_attrib_array(i, num) = attrib_array(i, j)
      attrib_ix(i, num) = j
    endif
  enddo
  attrib_num(i) = num
enddo

init_needed = .false.

end subroutine init_attribute_name_array

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_type (attrib_name) result (attrib_type)
!
! Routine to return the type (logical, integer, or real) of an attribute.
!
! Modules needed:
!   use bmad
!
! Input:
!   attrib_name -- Character(*): Name of the attribute. Must be upper case.
!
! Output:
!   attrib_type  -- Integer: is_logical$, is_integer$, or is_real$
!-

function attribute_type (attrib_name) result (attrib_type)

implicit none

character(*) attrib_name
integer attrib_type

!

select case (attrib_name)
case ('MATCH_END', 'MATCH_END_ORBIT')
  attrib_type = is_logical$
case ('PARTICLE', 'TAYLOR_ORDER', 'N_SLICE', 'N_REF_PASS', 'N_POLE', 'DIRECTION', 'IX_BRANCH_TO')
  attrib_type = is_integer$
case default
  attrib_type = is_real$
end select

end function attribute_type 

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

!

if (init_needed) call init_attribute_name_array
max_len = maxval(len_trim(attrib_array(1:n_key, 1:n_attrib_special_maxx)))

end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free1 (ix_ele, attrib_name, 
!                                 lat, err_print_flag, except_overlay) result (free)
!
! This function overloaded by attribute_free. See attribute_free for more details.
!-

function attribute_free1 (ix_ele, attrib_name, lat, err_print_flag, except_overlay) result (free)

implicit none

type (lat_struct) :: lat

integer ix_ele

character(*) attrib_name

logical free, do_print, do_except_overlay
logical, optional :: err_print_flag, except_overlay

!

do_print = logic_option (.true., err_print_flag)
do_except_overlay = logic_option(.false., except_overlay)

ele0 => lat%ele(ix_ele)
attrib_name0 = attrib_name

call check_this_attribute_free (ele0, attrib_name, lat, do_print, do_except_overlay, free)

end function attribute_free1

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free2 (ix_ele, ix_branch, attrib_name, 
!                                 lat, err_print_flag, except_overlay) result (free)
!
! This function overloaded by attribute_free. See attribute_free for more details.
!-

function attribute_free2 (ix_ele, ix_branch, attrib_name, lat, &
                                        err_print_flag, except_overlay) result (free)

implicit none

type (lat_struct), target :: lat

integer ix_branch, ix_ele, ix_ele0

character(*) attrib_name

logical free, do_print, do_except_overlay
logical, optional :: err_print_flag, except_overlay

character(16) :: r_name = 'attribute_free'

! init & check

do_print = logic_option (.true., err_print_flag)
do_except_overlay = logic_option(.false., except_overlay)

ele0 => lat%branch(ix_branch)%ele(ix_ele)
attrib_name0 = attrib_name

call check_this_attribute_free (ele0, attrib_name, lat, do_print, do_except_overlay, free)

end function attribute_free2

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free3 (loc, lat, err_print_flag, except_overlay) result (free)
!
! This function overloaded by attribute_free. See attribute_free for more details.
!-

function attribute_free3 (loc, attrib_name, lat, err_print_flag, except_overlay) result (free)

implicit none

type (lat_struct) :: lat
type (lat_ele_loc_struct) loc

character(*) attrib_name

logical free, do_print, do_except_overlay
logical, optional :: err_print_flag, except_overlay

!

do_print = logic_option (.true., err_print_flag)
do_except_overlay = logic_option(.false., except_overlay)

ele0 => lat%branch(loc%ix_branch)%ele(loc%ix_ele)
attrib_name0 = attrib_name

call check_this_attribute_free (ele0, attrib_name, lat, do_print, do_except_overlay, free)

end function attribute_free3

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

recursive subroutine check_this_attribute_free (ele, attrib_name, lat, &
                                          do_print, do_except_overlay, free, ix_lord)

type (ele_struct) :: ele
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele_p

integer ix_branch, i, ir, ix_attrib, ix
integer, optional :: ix_lord

character(*) attrib_name

logical free, do_print, do_except_overlay

! super_slaves attributes cannot be varied

free = .false.

ix_attrib = attribute_index(ele, attrib_name)
if (ix_attrib < 1) then
  if (do_print) call print_error (ele, ix_attrib, &
          'THIS ATTRIBUTE INDEX DOES NOT CORRESPOND TO A VALID ATTRIBUTE.')
  return
endif



! If the attribute is controled by an overlay lord then it cannot be varied.
! Exception: Multiple overlays can control the same attribute.

if (.not. do_except_overlay) then
  do i = ele%ic1_lord, ele%ic2_lord
    ix = lat%ic(i)
    ir = lat%control(ix)%ix_lord
    if (present(ix_lord)) then
      if (ix_lord == ir) cycle
      if (lat%ele(ix_lord)%lord_status == overlay_lord$) cycle
    endif
    if (lat%ele(ir)%lord_status == overlay_lord$) then
      if (lat%control(ix)%ix_attrib == ix_attrib) then 
        if (do_print) call print_error (ele, ix_attrib, & 
           'IT IS CONTROLLED BY THE OVERLAY_LORD: ' // lat%ele(ir)%name)
        return
      endif
    endif
  enddo
endif

! Check for a super_slave

if (ele%slave_status == super_slave$) then
  if (do_print) call print_error (ele, ix_attrib, 'THIS ELEMENT IS A SUPER_SLAVE.')
  return
endif

! Check for a multipass_slave.
! Exception: dphi0 can be varied for lcavity and rfcavity slaves.

if (ele%slave_status == multipass_slave$) then
  if ((ele%key /= lcavity$ .and. ele%key /= rfcavity$) .or. ix_attrib /= dphi0$) then
    if (do_print) call print_error (ele, ix_attrib, 'THIS ELEMENT IS A SUPER_SLAVE.')
    return
  endif
endif

! only one particular attribute of an overlay lord is allowed to be adjusted

if (ele%lord_status == overlay_lord$) then
  if (ix_attrib /= ele%ix_value) then
    if (do_print) call print_error (ele, ix_attrib, &
           'FOR THIS OVERLAY ELEMENT THE ATTRIBUTE TO VARY IS: ' // ele%attribute_name)
    return
  endif
endif

if (ele%lord_status == group_lord$) then
  if (ix_attrib /= command$ .and. ix_attrib /= old_command$) then
    if (do_print) call print_error (ele, ix_attrib, &
          'FOR THIS GROUP ELEMENT THE ATTRIBUTE TO VARY IS: "COMMAND" OR "OLD_COMMAND"')
    return
  endif
endif

! check if it is a dependent variable.

free = .true.

select case (ele%key)
case (sbend$)
  if (any(ix_attrib == (/ angle$, l_chord$, rho$ /))) free = .false.
case (rfcavity$)
  if (ix_attrib == rf_frequency$ .and. ele%value(harmon$) /= 0) free = .false.
case (beambeam$)
  if (ix_attrib == bbi_const$) free = .false.
case (wiggler$)
  if (ix_attrib == k1$ .or. ix_attrib == rho$) free = .false. 
case (lcavity$)
  if (any(ix_attrib == (/ delta_e$, p0c_start$, e_tot_start$ /))) free = .false.
case (elseparator$)
  if (ix_attrib == e_field$ .or. ix_attrib == voltage$) free = .false.
end select

if (ix_attrib == tilt_tot$) free = .false.
if (ix_attrib == x_pitch_tot$) free = .false.
if (ix_attrib == y_pitch_tot$) free = .false.
if (ix_attrib == x_offset_tot$) free = .false.
if (ix_attrib == y_offset_tot$) free = .false.
if (ix_attrib == s_offset_tot$) free = .false.
if (ix_attrib == e_tot$) free = .false.
if (ix_attrib == p0c$) free = .false.

if (ele%key == sbend$ .and. ele%lord_status == multipass_lord$ .and. &
    ele%value(n_ref_pass$) == 0 .and. ix_attrib == p0c$) free = .true.

if (.not.free) then
  if (do_print) call print_error (ele, ix_attrib, 'THE ATTRIBUTE IS A DEPENDENT VARIABLE.')
  return
endif

! field_master on means that the b_field and bn_gradient values control
! the strength.

if (ele%field_master) then
  select case (ele%key)
  case (quadrupole$)
    if (ix_attrib == k1$) free = .false.
  case (sextupole$)
    if (ix_attrib == k2$) free = .false.
  case (octupole$)
    if (ix_attrib == k3$) free = .false.
  case (solenoid$)
    if (ix_attrib == ks$) free = .false.
  case (sol_quad$)
    if (ix_attrib == ks$) free = .false.
    if (ix_attrib == k1$) free = .false.
    if (ix_attrib == k2$) free = .false.
  case (sbend$)
    if (ix_attrib == g$) free = .false.
    if (ix_attrib == g_err$) free = .false.
  case (hkicker$, vkicker$)
    if (ix_attrib == kick$) free = .false.
  end select

  if (ix_attrib == hkick$) free = .false.
  if (ix_attrib == vkick$) free = .false.

else
  select case (ele%key)
  case (quadrupole$)
    if (ix_attrib == b1_gradient$) free = .false.
  case (sextupole$)
    if (ix_attrib == b2_gradient$) free = .false.
  case (octupole$)
    if (ix_attrib == b3_gradient$) free = .false.
  case (solenoid$)
    if (ix_attrib == bs_field$) free = .false.
  case (sol_quad$)
    if (ix_attrib == bs_field$) free = .false.
    if (ix_attrib == b1_gradient$) free = .false.
  case (sbend$)
    if (ix_attrib == b_field$) free = .false.
    if (ix_attrib == b_field_err$) free = .false.
    if (ix_attrib == b1_gradient$) free = .false.
    if (ix_attrib == b2_gradient$) free = .false.
  case (hkicker$, vkicker$)
    if (ix_attrib == bl_kick$) free = .false.
  end select

  if (ix_attrib == bl_hkick$) free = .false.
  if (ix_attrib == bl_vkick$) free = .false.

endif

if (.not. free) then
  if (do_print) call print_error (ele, ix_attrib, &
       "THE ATTRIBUTE IS A DEPENDENT VARIABLE SINCE", &
       "THE ELEMENT'S FIELD_MASTER IS " // on_off_logic (ele%field_master))
  return
endif

! check slaves

if (ele%lord_status == group_lord$ .or. ele%lord_status == overlay_lord$) then
  do i = ele%ix1_slave, ele%ix2_slave
    ele_p => lat%ele(lat%control(i)%ix_slave)
    call check_this_attribute_free (ele_p, attribute_name(ele_p, lat%control(i)%ix_attrib), &
                                    lat, do_print, do_except_overlay, free, ele%ix_ele)
    if (.not. free) return
  enddo
endif

end subroutine

!-------------------------------------------------------

subroutine print_error (ele, ix_attrib, l1, l2)

type (ele_struct) ele

integer ix_attrib, nl

character(*) l1
character(*), optional :: l2
character(100) li(8)
character (20) :: r_name = 'attribute_free'

!

nl = 0

nl=nl+1; li(nl) =   'THE ATTRIBUTE: ' // attrib_name0
nl=nl+1; li(nl) =   'OF THE ELEMENT: ' // ele0%name

if (ele%ix_branch == 0 .and. ele%ix_ele == ele0%ix_ele) then
  nl=nl+1; li(nl) = 'IS NOT FREE TO VARY SINCE:'
else 
  nl=nl+1; li(nl) = 'IS NOT FREE TO VARY SINCE IT IS TRYING TO CONTROL:'
  nl=nl+1; li(nl) = 'THE ATTRIBUTE: ' // attrib_name0
  nl=nl+1; li(nl) = 'OF THE ELEMENT: ' // ele%name
  nl=nl+1; li(nl) = 'AND THIS IS NOT FREE TO VARY SINCE'
endif

nl=nl+1; li(nl) = l1
if (present(l2)) then
  nl=nl+1; li(nl) = l2
endif

call out_io (s_error$, r_name, li(1:nl))   

end subroutine

end module
