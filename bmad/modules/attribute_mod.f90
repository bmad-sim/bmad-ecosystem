module attribute_mod

use bmad_routine_interface

implicit none

! Define parameters for attrib_array()%state and why_not_free argument in attribute_free(...) 
! Note: Slot 2 is reserved for super_slave$ and slot 9 is reserved for multipass_slave$

integer, parameter :: does_not_exist$ = -1, is_free$ = 0, quasi_free$ = 1, dependent$ = 3, private$ = 4 
integer, parameter :: overlay_slave$ = 5, field_master_dependent$ = 6, super_lord_align$ = 7

! attrib_array(key, ix_param)%state gives the state of the attribute.
! This may be one of:
!   does_not_exist$   -- Does not exist.
!   is_free$          -- Free to vary as long as attribute has not controlled by another 
!                         element (overlay, super_lord, etc.)
!   quasi_free$       -- May be free or not. For example, k1 is only free if field_master = F.
!   dependent$        -- Value calculated by Bmad. Cannot be user varied as an independent parameter.
!   private$          -- Internal parameter used in calculations. Will not be displayed by type_ele.
!   super_lord_align$ -- A super_lord alignment attribute, except for tilt, may not be changed if any super_slave 
!                         is not an em_field element and that super_slave has a second lord that is not a pipe.
!                         Exception: Solenoid and Quadrupole overlapping.

type ele_attribute_struct
  character(40) :: name = null_name$
  integer :: state = does_not_exist$  ! See above.
  integer :: kind = unknown$          ! Is_switch$, is_real$, etc. See attribute_type routine.
  character(16) :: units = ''         ! EG: 'T*m'.
  integer :: ix_attrib = -1           ! Attribute index. Frequently will be where in the 
                                      !   ele%value(:) array the attribute is.
  real(rp) :: value = real_garbage$   ! Used by type_ele.
end type

type (ele_attribute_struct), private, save :: attrib_array(n_key$, num_ele_attrib_extended$)

character(40), private, save :: short_attrib_array(n_key$, num_ele_attrib_extended$)
integer, private, save :: attrib_num(n_key$)
integer, private, save :: attrib_ix(n_key$, num_ele_attrib_extended$)
logical, private, save :: attribute_array_init_needed = .true.
logical, private, save :: has_orientation_attributes_key(n_key$)
private init_short_attrib_array

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free
!
! Overloaded function for:
!   Function attribute_free1 (ix_ele, attrib_name, lat, err_print_flag, 
!                                except_overlay, dependent_attribs_free, why_not_free) result (free)
!   Function attribute_free2 (ele, attrib_name, err_print_flag, 
!                                except_overlay, dependent_attribs_free, why_not_free) result (free)
!   Function attribute_free3 (ix_ele, ix_branch, attrib_name, lat, err_print_flag, 
!                                except_overlay, why_not_free) result (free)
!
! Routine to check if an attribute is free to vary.
!
! Attributes that cannot be changed directly include super_slave attributes (since
! these attributes are controlled by their super_lords) and attributes that
! are controlled by an overlay.
!
! Also dependent variables such as the angle of a bend cannot be 
!   freely variable.
!
! Input:
!   ix_ele                  -- integer: Index of element in element array.
!   ix_branch               -- integer: Branch index of element. 
!   ele                     -- ele_struct: Element containing the attribute
!   attrib_name             -- character(*): Name of the attribute. Assumed upper case.
!   lat                     -- lat_struct: Lattice structure.
!   err_print_flag          -- logical, optional: If present and False then suppress
!                               printing of an error message if attribute is not free.
!   except_overlay          -- logical, optional: If present and True then an attribute that
!                               is controlled by an overlay will be treated as free. 
!                               This is used by, for example, the create_overlay routine.
!   dependent_attribs_free  -- logical, optional: If present and True then mark as free 
!                               attributes that are dependent. For example, if ele%field_master = F,
!                               b1_field is dependent upon k1. Default is False. Use True when
!                               using intelligent bookkeeping.
!
! Output:
!   free          -- logical: Set True if attribtute not found or attriubte
!                       cannot be changed directly.
!   why_not_free  -- integer, optional: Possibilities are:
!                         field_master_dependent$  -> Dependent due to setting of ele%field_master.
!                         dependent$               -> Not field_master_dependent$ but value is 
!                                                       dependent upon the value of other attributes.
!                         does_not_exist$          -> Attribute name is unrecognized or does not exist
!                                                       for the type of element.
!                         overlay_slave$           -> Attribute is controlled by an overlay lord.
!                         super_slave$             -> Attribute is controlled by element's super_lord.
!                         multipass_slave$         -> Attribute is controlled by element's multipass_lord.
!-

interface attribute_free
  module procedure attribute_free1
  module procedure attribute_free2
  module procedure attribute_free3
end interface

private check_this_attribute_free

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
! See also:
!   has_attribute
!   attribute_info
!   attribute_name
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
! Function to return the name of an attribute for a particular type of 
! Bmad element. 
!
! This routine is an overloaded name for:
!   attribute_name1 (ele, ix_att, show_private) result (attrib_name)
!   attribute_name2 (key, ix_att, show_private) result (attrib_name)
!
!
! Note: attribute_name (key, ix_att) is not able to handle overlay/group control variables.
! Use attributge_name (ele, ix_att) is this is needed.
!
! Input:
!   ele             -- ele_struct: 
!     %key             -- Integer: Key name of element type (e.g. SBEND$, etc.)
!   key             -- integer: Key name of element type (e.g. sbend$, etc.)
!   ix_att          -- integer: Index of attribute (e.g. k1$)
!   show_private    -- logical, optional: If False (default) return null_name$ for private attributes.
!
! Output:
!   attrib_name     -- Character(40): Name of attribute. First character is a "!" if there is a problem.
!                       Will always be upper case (even with private attributes).
!                         = "!BAD ELE KEY"           %key is invalid
!                         = "!BAD INDEX"             ix_att is invalid (out of range).
!                         = "!NULL" (null_name$)     ix_att does not correspond to an attribute or is private.
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

! Note: ele%control or ele%control%var may not be allocated during parsing.

if ((ele%key == group$ .or. ele%key == overlay$ .or. ele%key == ramper$) .and. associated(ele%control)) then
  if (allocated(ele%control%var)) then
    n = min(4, len(name))
    if (name(1:n) == 'OLD_') then
      do i = 1, size(ele%control%var)
        if (name(5:) /= ele%control%var(i)%name) cycle
        attrib_index = i + old_control_var_offset$
        if (present(full_name)) full_name = name
        return
      enddo
    else  
      do i = 1, size(ele%control%var)
        if (name /= ele%control%var(i)%name) cycle
        attrib_index = i + var_offset$
        if (present(full_name)) full_name = name
        return
      enddo
    endif
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

integer i, j, k, key, num, ilen, n_abbrev, ix_abbrev, i0, ixs
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

select case (name)
case ('G_ERR');       name40 = 'DG'
case ('B_FIELD_ERR'); name40 = 'DB_FIELD'
case ('REF');         name40 = 'REFERENCE'
case default;         name40 = name
end select

if (key == rbend$) then  ! Note: Rbends only exist while parsing
  if (name == 'L')     name40 = 'L_CHORD'
endif

!

if (ilen > 2 .and. name40(1:2) == 'TT' .and. (key == 0 .or. key == taylor$)) then
  if (name40(3:3) == 'S') then
    ixs = index('1XYZ', name40(4:4))
    if (ixs == 0) return
    i0 = 5
  else
    ixs = 0
    i0 = 3
  endif

  do i = i0, ilen
    if (index('123456', name(i:i)) == 0) return
  enddo
  if (name40(i0:) == '') then
    attrib_index = 0
  else
    if (.not. is_integer(name40(i0:), attrib_index)) return
  endif
  attrib_index = attrib_index + taylor_offset$ + ixs * (taylor_offset$/10)
  full_name = name
  return
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

  if (key == rbend$ .and. name40 == 'L_ARC') then
    attrib_index = l$
    if (present(full_name)) full_name = 'L_ARC'
  endif 

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
! Function attribute_name1 (key, ix_att, show_private) result (attrib_name)
!
! Overloaded by attribute_name. See attribute_name for more details.
!-

function attribute_name1 (key, ix_att, show_private) result (attrib_name)

type (ele_struct) ele
integer i, key, ix_att, ix
character(40) attrib_name
logical, optional :: show_private

!

ele%key = key
attrib_name = attribute_name2(ele, ix_att, show_private)

end function attribute_name1

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_name2 (ele, ix_att, show_private) result (attrib_name)
!
! Overloaded by attribute_name. See attribute_name for more details.
!-

function attribute_name2 (ele, ix_att, show_private) result (attrib_name)

type (ele_struct) ele
integer i, key, ix_att, ix, ixs
character(40) attrib_name
character(10) str
character(1), parameter :: spn(4) = ['1', 'X', 'Y', 'Z']
logical, optional :: show_private

!

if (attribute_array_init_needed) call init_attribute_name_array()

key = ele%key

if (key <= 0 .or. key > n_key$) then
  attrib_name = '!BAD ELE KEY'

elseif (key == taylor$ .and. ix_att > taylor_offset$) then
  ixs = (10 * (ix_att - taylor_offset$)) / taylor_offset$

  write (str, '(i0)') ix_att - taylor_offset$ - ixs * (taylor_offset$/10)
  if (str(1:1) == ' ') return
  do i = 1, len(str)
    if (str(i:i) == ' ') exit
    if (index('123456', str(i:i)) == 0) return
  enddo

  if (ixs == 0) then
    attrib_name = 'TT' // str(1:i-1)
  else
    attrib_name = 'TTS' // spn(ixs) // str(1:i-1)
  endif

elseif ((key == group$ .or. key == overlay$) .and. is_attribute(ix_att, control_var$)) then
  ix = ix_att - var_offset$
  if (ix > size(ele%control%var)) then
    attrib_name = '!BAD INDEX'
  else
    attrib_name = ele%control%var(ix)%name
  endif

elseif (key == group$ .and. is_attribute(ix_att, old_control_var$)) then
  ix = ix_att - old_control_var_offset$
  if (ix > size(ele%control%var)) then
    attrib_name = '!BAD INDEX'
  else
    attrib_name = 'OLD_' // ele%control%var(ix)%name
  endif

elseif (ix_att <= 0 .or. ix_att > num_ele_attrib_extended$) then
  attrib_name = '!BAD INDEX'

else
  if (attrib_array(key, ix_att)%state == private$) then
    if (logic_option(.false., show_private)) then
      attrib_name = upcase(attrib_array(key, ix_att)%name)
    else
      attrib_name = null_name$
    endif
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
! Input:
!   ele    -- Ele_struct: 
!     %key    -- Integer: Key name of element type (e.g. SBEND$, etc.)
!   ix_att -- Integer: Index of attribute (e.g. k1$)
!
! Output:
!   attrib_info -- ele_attribute_struct: Info on this attribute.
!     Note: %value is not set since this info is contained in the ele argument.
!           Use pointer_to_attribute to access the attribute.
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

attrib_info%ix_attrib = ix_att
attrib_info%state = does_not_exist$
attrib_info%name = '!BAD INDEX'

! Taylor term

if (ele%key == taylor$ .and. ix_att > taylor_offset$) then
  attrib_info%name = attribute_name2(ele, ix_att)
  if (attrib_info%name /= '!BAD INDEX') attrib_info%state = is_free$
  return
endif

! Overlay and group control vars.

if ((ele%key == group$ .or. ele%key == overlay$) .and. is_attribute(ix_att, control_var$)) then
  ix = ix_att - var_offset$
  if (ix > size(ele%control%var)) return
  attrib_info%name = ele%control%var(ix)%name
  attrib_info%state = is_free$
  return
endif

if (ele%key == group$ .and. ix_att > old_control_var_offset$) then
  ix = ix_att - old_control_var_offset$
  if (ix > size(ele%control%var)) return
  attrib_info%name = 'OLD_' // ele%control%var(ix)%name
  attrib_info%state = is_free$
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

type (surface_curvature_struct) curv
type (ele_struct) ele
integer i, j, num, ix, iy
character(40) word

! 

if (.not. attribute_array_init_needed) return
attribute_array_init_needed = .false.

do i = 1, n_key$

  if (i == def_bmad_com$)         cycle
  if (i == def_space_charge_com$) cycle
  if (i == def_mad_beam$)         cycle
  if (i == def_particle_start$)   cycle
  if (i == def_line$)             cycle
  if (i == def_parameter$)        cycle
  if (i == def_ptc_com$)          cycle

  call init_attribute_name1 (i, check_sum$, 'check_sum', private$)

  if (i == beginning_ele$)  cycle

  call init_attribute_name1 (i, type$,      'TYPE')
  call init_attribute_name1 (i, alias$,     'ALIAS')
  call init_attribute_name1 (i, descrip$,   'DESCRIP')

  if (i == group$)    cycle
  if (i == overlay$)  cycle
  if (i == girder$)   cycle
  if (i == ramper$)   cycle
  if (i == feedback$)   cycle

  call init_attribute_name1 (i, superimpose$,        'SUPERIMPOSE')
  call init_attribute_name1 (i, super_offset$,       'OFFSET')
  call init_attribute_name1 (i, reference$,          'REFERENCE')
  call init_attribute_name1 (i, ref_origin$,         'REF_ORIGIN')
  call init_attribute_name1 (i, ele_origin$,         'ELE_ORIGIN')
  call init_attribute_name1 (i, wrap_superimpose$,   'WRAP_SUPERIMPOSE')

  if (i == null_ele$) cycle

  call init_attribute_name1 (i, mat6_calc_method$,       'MAT6_CALC_METHOD')
  call init_attribute_name1 (i, tracking_method$,        'TRACKING_METHOD')
  call init_attribute_name1 (i, spin_tracking_method$,   'SPIN_TRACKING_METHOD')
  call init_attribute_name1 (i, ptc_integration_type$,   'PTC_INTEGRATION_TYPE')

  select case (i)   ! Does not include beginning_ele
  case (lcavity$, em_field$, custom$)
    call init_attribute_name1 (i, E_tot$,                  'E_TOT', dependent$)
    call init_attribute_name1 (i, p0c$,                    'P0C', dependent$)  
    call init_attribute_name1 (i, p0c_start$,              'P0C_START', quasi_free$)
    call init_attribute_name1 (i, e_tot_start$,            'E_TOT_START', quasi_free$)
  case (converter$)
    call init_attribute_name1 (i, E_tot$,                  'E_TOT')
    call init_attribute_name1 (i, p0c$,                    'P0C')
    call init_attribute_name1 (i, p0c_start$,              'P0C_START', dependent$)
    call init_attribute_name1 (i, e_tot_start$,            'E_TOT_START', dependent$)
  case (patch$)
    call init_attribute_name1 (i, E_tot$,                  'E_TOT', quasi_free$) ! Free in multipass_lord
    call init_attribute_name1 (i, p0c$,                    'P0C', quasi_free$)   ! Free in multipass_lord
    call init_attribute_name1 (i, p0c_start$,              'P0C_START', dependent$)
    call init_attribute_name1 (i, e_tot_start$,            'E_TOT_START', dependent$)
  case default
    call init_attribute_name1 (i, E_tot$,                  'E_TOT', quasi_free$) ! Free in multipass_lord
    call init_attribute_name1 (i, p0c$,                    'P0C', quasi_free$)   ! Free in multipass_lord
    call init_attribute_name1 (i, e_tot_start$,            'e_tot_start', private$)
    call init_attribute_name1 (i, p0c_start$,              'p0c_start', private$)
  end select

  call init_attribute_name1 (i, delta_ref_time$,         'DELTA_REF_TIME', dependent$)
  call init_attribute_name1 (i, ref_time_start$,         'REF_TIME_START', dependent$)

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

  if (i == hybrid$)        cycle
  if (i == match$)         cycle
  if (i == photon_fork$)   cycle
  if (i == fork$)          cycle
  if (i == gkicker$)       cycle

  call init_attribute_name1 (i, x_offset$,      'X_OFFSET')
  call init_attribute_name1 (i, y_offset$,      'Y_OFFSET')
  call init_attribute_name1 (i, z_offset$,      'Z_OFFSET')
  call init_attribute_name1 (i, x_pitch$,       'X_PITCH')
  call init_attribute_name1 (i, y_pitch$,       'Y_PITCH')
  call init_attribute_name1 (i, tilt$,          'TILT' )
  call init_attribute_name1 (i, wall$,          'WALL' )

  if (i == floor_shift$)  cycle
  if (i == patch$)        cycle

  call init_attribute_name1 (i, tilt_tot$,        'TILT_TOT', dependent$)
  call init_attribute_name1 (i, x_offset_tot$,    'X_OFFSET_TOT', dependent$)
  call init_attribute_name1 (i, y_offset_tot$,    'Y_OFFSET_TOT', dependent$)
  call init_attribute_name1 (i, z_offset_tot$,    'Z_OFFSET_TOT', dependent$)
  call init_attribute_name1 (i, x_pitch_tot$,     'X_PITCH_TOT', dependent$)
  call init_attribute_name1 (i, y_pitch_tot$,     'Y_PITCH_TOT', dependent$)
  call init_attribute_name1 (i, dispatch$,        'dispatch', private$)

  if (i == multilayer_mirror$) cycle
  if (i == mirror$)            cycle
  if (i == crystal$)           cycle
  if (i == sample$)            cycle
  if (i == capillary$)         cycle
  if (i == lens$)              cycle
  if (i == photon_init$)       cycle

  if (i /= drift$) call init_attribute_name1 (i, is_on$,        'IS_ON')

  if (i == diffraction_plate$) cycle
  if (i == mask$)              cycle
  if (i == detector$)          cycle
  if (i == beambeam$)          cycle
  if (i == multipole$)         cycle 
  if (i == ab_multipole$)      cycle
  if (i == marker$)            cycle

  call init_attribute_name1 (i, l$,                    'L')

  if (i == converter$)         cycle
  if (i == foil$)              cycle
  if (i == pickup$)            cycle

  call init_attribute_name1 (i, symplectify$,          'SYMPLECTIFY')
  call init_attribute_name1 (i, taylor_map_includes_offsets$,    'TAYLOR_MAP_INCLUDES_OFFSETS')

  call init_attribute_name1 (i, lord_pad1$,            'LORD_PAD1', quasi_free$)
  call init_attribute_name1 (i, lord_pad2$,            'LORD_PAD2', quasi_free$)

  if (i == taylor$)            cycle

  call init_attribute_name1 (i, integrator_order$,     'INTEGRATOR_ORDER')
  call init_attribute_name1 (i, field_calc$,           'FIELD_CALC')
  call init_attribute_name1 (i, csr_method$,           'CSR_METHOD')
  call init_attribute_name1 (i, csr_ds_step$,          'CSR_DS_STEP')
  call init_attribute_name1 (i, space_charge_method$,  'SPACE_CHARGE_METHOD')
  call init_attribute_name1 (i, multipass_ref_energy$, 'MULTIPASS_REF_ENERGY', private$)
  call init_attribute_name1 (i, static_linear_map$,    'STATIC_LINEAR_MAP')

  if (i == sad_mult$)          cycle

  call init_attribute_name1 (i, num_steps$,            'NUM_STEPS', quasi_free$)
  call init_attribute_name1 (i, ds_step$,              'DS_STEP')

  if (i == drift$)             cycle

  call init_attribute_name1 (i, field_overlaps$,       'FIELD_OVERLAPS')

  ! Markers will also have these wake attributes. See below.
  call init_attribute_name1 (i, sr_wake$,                   'SR_WAKE')
  call init_attribute_name1 (i, lr_wake$,                   'LR_WAKE')
  call init_attribute_name1 (i, sr_wake_file$,              'SR_WAKE_FILE')
  call init_attribute_name1 (i, lr_wake_file$,              'LR_WAKE_FILE')
  call init_attribute_name1 (i, lr_freq_spread$,            'LR_FREQ_SPREAD')
  call init_attribute_name1 (i, lr_self_wake_on$,           'LR_SELF_WAKE_ON')

  if (i == pipe$)         cycle
  if (i == custom$)       cycle

  if (i /= crab_cavity$) then
    call init_attribute_name1 (i, fringe_type$,        'FRINGE_TYPE')
    call init_attribute_name1 (i, spin_fringe_on$,     'SPIN_FRINGE_ON')
    call init_attribute_name1 (i, fringe_at$,          'FRINGE_AT')
  endif

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
    call init_attribute_name1 (i, x_gain_err$,          'X_GAIN_ERR')
    call init_attribute_name1 (i, y_gain_err$,          'Y_GAIN_ERR')
    call init_attribute_name1 (i, crunch$,              'CRUNCH')
    call init_attribute_name1 (i, noise$,               'NOISE')
    call init_attribute_name1 (i, tilt_calib$,          'TILT_CALIB')
    call init_attribute_name1 (i, x_gain_calib$,        'X_GAIN_CALIB')
    call init_attribute_name1 (i, y_gain_calib$,        'Y_GAIN_CALIB')
    call init_attribute_name1 (i, crunch_calib$,        'CRUNCH_CALIB')
    call init_attribute_name1 (i, x_offset_calib$,      'X_OFFSET_CALIB')
    call init_attribute_name1 (i, y_offset_calib$,      'Y_OFFSET_CALIB')
    call init_attribute_name1 (i, n_sample$,            'N_SAMPLE')
    call init_attribute_name1 (i, de_eta_meas$,         'DE_ETA_MEAS')
    call init_attribute_name1 (i, osc_amplitude$,       'OSC_AMPLITUDE')
    call init_attribute_name1 (i, x_dispersion_err$,    'X_DISPERSION_ERR')
    call init_attribute_name1 (i, y_dispersion_err$,    'Y_DISPERSION_ERR')
    call init_attribute_name1 (i, x_dispersion_calib$,  'X_DISPERSION_CALIB')
    call init_attribute_name1 (i, y_dispersion_calib$,  'Y_DISPERSION_CALIB')
    call init_attribute_name1 (i, split_id$,            'split_id', private$)
  end select
enddo

!

do i = 1, n_key$
  select case (i)
  case (crystal$, multilayer_mirror$, mirror$, sample$, diffraction_plate$, detector$)
    call init_attribute_name1 (i, p89$, 'DISPLACEMENT')
    call init_attribute_name1 (i, p90$, 'SEGMENTED')
    num = a0$ - 1
    do ix = 0, ubound(curv%xy, 1)
    do iy = 0, ubound(curv%xy, 2)
      if (ix+iy < 2) cycle
      if (ix+iy > ubound(curv%xy, 1)) cycle
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
  case (ac_kicker$, elseparator$, kicker$, octupole$, quadrupole$, sbend$, rbend$, &
         sextupole$, solenoid$, sol_quad$, ab_multipole$, wiggler$, undulator$, &
         hkicker$, vkicker$, sad_mult$, thick_multipole$, rfcavity$, lcavity$, crab_cavity$)
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
    attrib_array(i, a0$:b21$)%state = is_free$
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
call init_attribute_name1 (photon_fork$, ref_species$,              'REF_SPECIES', dependent$)

attrib_array(fork$, :) = attrib_array(photon_fork$, :)

call init_attribute_name1 (beambeam$, l$,                           'L', dependent$)
call init_attribute_name1 (beambeam$, ks$,                          'KS', quasi_free$)
call init_attribute_name1 (beambeam$, bs_field$,                    'BS_FIELD', quasi_free$)
call init_attribute_name1 (beambeam$, sig_x$,                       'SIG_X')
call init_attribute_name1 (beambeam$, sig_y$,                       'SIG_Y')
call init_attribute_name1 (beambeam$, sig_z$,                       'SIG_Z')
call init_attribute_name1 (beambeam$, bbi_const$,                   'BBI_CONSTANT', dependent$)
call init_attribute_name1 (beambeam$, charge$,                      'CHARGE')
call init_attribute_name1 (beambeam$, n_particle$,                  'N_PARTICLE')
call init_attribute_name1 (beambeam$, n_slice$,                     'N_SLICE')
call init_attribute_name1 (beambeam$, symplectify$,                 'SYMPLECTIFY')
call init_attribute_name1 (beambeam$, field_calc$,                  'FIELD_CALC')
call init_attribute_name1 (beambeam$, beta_a_strong$,               'BETA_A_STRONG')
call init_attribute_name1 (beambeam$, beta_b_strong$,               'BETA_B_STRONG')
call init_attribute_name1 (beambeam$, alpha_a_strong$,              'ALPHA_A_STRONG')
call init_attribute_name1 (beambeam$, alpha_b_strong$,              'ALPHA_B_STRONG')
call init_attribute_name1 (beambeam$, cmat_11$,                     'CMAT_11')
call init_attribute_name1 (beambeam$, cmat_12$,                     'CMAT_12')
call init_attribute_name1 (beambeam$, cmat_21$,                     'CMAT_21')
call init_attribute_name1 (beambeam$, cmat_22$,                     'CMAT_22')
call init_attribute_name1 (beambeam$, crab_x1$,                     'CRAB_X1')
call init_attribute_name1 (beambeam$, crab_x2$,                     'CRAB_X2')
call init_attribute_name1 (beambeam$, crab_x3$,                     'CRAB_X3')
call init_attribute_name1 (beambeam$, crab_x4$,                     'CRAB_X4')
call init_attribute_name1 (beambeam$, crab_x5$,                     'CRAB_X5')
call init_attribute_name1 (beambeam$, crab_tilt$,                   'CRAB_TILT')
call init_attribute_name1 (beambeam$, crossing_time$,               'CROSSING_TIME')
call init_attribute_name1 (beambeam$, s_twiss_ref$,                 'S_TWISS_REF')
call init_attribute_name1 (beambeam$, repetition_frequency$,        'REPETITION_FREQUENCY')
call init_attribute_name1 (beambeam$, rf_clock_harmonic$,           'rf_clock_harminic', private$)
call init_attribute_name1 (beambeam$, species_strong$,              'SPECIES_STRONG')
call init_attribute_name1 (beambeam$, e_tot_strong$,                'E_TOT_STRONG')
call init_attribute_name1 (beambeam$, pc_strong$,                   'PC_STRONG')

call init_attribute_name1 (beginning_ele$, l$,                           'l', private$)
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
call init_attribute_name1 (beginning_ele$, beta_a$,                      'BETA_A')
call init_attribute_name1 (beginning_ele$, beta_b$,                      'BETA_B')
call init_attribute_name1 (beginning_ele$, alpha_a$,                     'ALPHA_A')
call init_attribute_name1 (beginning_ele$, alpha_b$,                     'ALPHA_B')
call init_attribute_name1 (beginning_ele$, dbeta_dpz_a$,                 'DBETA_DPZ_A')
call init_attribute_name1 (beginning_ele$, dbeta_dpz_b$,                 'DBETA_DPZ_B')
call init_attribute_name1 (beginning_ele$, dalpha_dpz_a$,                'DALPHA_DPZ_A')
call init_attribute_name1 (beginning_ele$, dalpha_dpz_b$,                'DALPHA_DPZ_B')
call init_attribute_name1 (beginning_ele$, deta_dpz_x$,                  'DETA_DPZ_X')
call init_attribute_name1 (beginning_ele$, deta_dpz_y$,                  'DETA_DPZ_Y')
call init_attribute_name1 (beginning_ele$, detap_dpz_x$,                 'DETAP_DPZ_X')
call init_attribute_name1 (beginning_ele$, detap_dpz_y$,                 'DETAP_DPZ_Y')
call init_attribute_name1 (beginning_ele$, eta_x$,                       'ETA_X')
call init_attribute_name1 (beginning_ele$, eta_y$,                       'ETA_Y')
call init_attribute_name1 (beginning_ele$, eta_z$,                       'ETA_Z')
call init_attribute_name1 (beginning_ele$, etap_x$,                      'ETAP_X')
call init_attribute_name1 (beginning_ele$, etap_y$,                      'ETAP_Y')
call init_attribute_name1 (beginning_ele$, phi_a$,                       'PHI_A')
call init_attribute_name1 (beginning_ele$, phi_b$,                       'PHI_B')
call init_attribute_name1 (beginning_ele$, cmat_11$,                     'CMAT_11')
call init_attribute_name1 (beginning_ele$, cmat_12$,                     'CMAT_12')
call init_attribute_name1 (beginning_ele$, cmat_21$,                     'CMAT_21')
call init_attribute_name1 (beginning_ele$, cmat_22$,                     'CMAT_22')
call init_attribute_name1 (beginning_ele$, mode_flip$,                   'MODE_FLIP')
call init_attribute_name1 (beginning_ele$, spin_dn_dpz_x$,               'SPIN_DN_DPZ_X')
call init_attribute_name1 (beginning_ele$, spin_dn_dpz_y$,               'SPIN_DN_DPZ_Y')
call init_attribute_name1 (beginning_ele$, spin_dn_dpz_z$,               'SPIN_DN_DPZ_Z')
call init_attribute_name1 (beginning_ele$, s_long$,                      'S')
call init_attribute_name1 (beginning_ele$, ref_time$,                    'REF_TIME')
call init_attribute_name1 (beginning_ele$, inherit_from_fork$,           'INHERIT_FROM_FORK')
call init_attribute_name1 (beginning_ele$, deta_ds_master$,              'deta_ds_master', private$)

call init_attribute_name1 (feedback$, input_ele$,                        'INPUT_ELE')
call init_attribute_name1 (feedback$, output_ele$,                       'OUTPUT_ELE')

attrib_array(def_line$, :) = attrib_array(beginning_ele$, :)
call init_attribute_name1 (def_line$, particle$,                    'PARTICLE')
call init_attribute_name1 (def_line$, live_branch$,                 'LIVE_BRANCH')
call init_attribute_name1 (def_line$, geometry$,                    'GEOMETRY')
call init_attribute_name1 (def_line$, default_tracking_species$,    'DEFAULT_TRACKING_SPECIES')
call init_attribute_name1 (def_line$, ix_branch$,                   'ix_branch', private$)
call init_attribute_name1 (def_line$, high_energy_space_charge_on$, 'HIGH_ENERGY_SPACE_CHARGE_ON')

call init_attribute_name1 (def_mad_beam$, particle$,                      'PARTICLE')
call init_attribute_name1 (def_mad_beam$, e_tot$,                         'ENERGY')
call init_attribute_name1 (def_mad_beam$, p0c$,                           'PC')
call init_attribute_name1 (def_mad_beam$, n_part$,                        'N_PART')

call init_attribute_name1 (def_particle_start$, x$,                         'X')
call init_attribute_name1 (def_particle_start$, px$,                        'PX')
call init_attribute_name1 (def_particle_start$, y$,                         'Y')
call init_attribute_name1 (def_particle_start$, py$,                        'PY')
call init_attribute_name1 (def_particle_start$, z$,                         'Z')
call init_attribute_name1 (def_particle_start$, pz$,                        'PZ')
call init_attribute_name1 (def_particle_start$, field_x$,                   'FIELD_X')
call init_attribute_name1 (def_particle_start$, field_y$,                   'FIELD_Y')
call init_attribute_name1 (def_particle_start$, phase_x$,                   'PHASE_X')
call init_attribute_name1 (def_particle_start$, phase_y$,                   'PHASE_Y')
call init_attribute_name1 (def_particle_start$, t$,                         'T')
call init_attribute_name1 (def_particle_start$, e_photon$,                  'E_PHOTON')
call init_attribute_name1 (def_particle_start$, spin_x$,                    'SPIN_X')
call init_attribute_name1 (def_particle_start$, spin_y$,                    'SPIN_Y')
call init_attribute_name1 (def_particle_start$, spin_z$,                    'SPIN_Z')
call init_attribute_name1 (def_particle_start$, emittance_a$,               'EMITTANCE_A')
call init_attribute_name1 (def_particle_start$, emittance_b$,               'EMITTANCE_B')
call init_attribute_name1 (def_particle_start$, emittance_z$,               'EMITTANCE_Z')
call init_attribute_name1 (def_particle_start$, sig_pz$,                    'SIG_PZ')
call init_attribute_name1 (def_particle_start$, sig_z$,                     'SIG_Z')

call init_attribute_name1 (def_parameter$, ix_branch$,                    'ix_branch', private$)
call init_attribute_name1 (def_parameter$, e_tot$,                        'E_TOT')
call init_attribute_name1 (def_parameter$, p0c$,                          'P0C')
call init_attribute_name1 (def_parameter$, live_branch$,                  'LIVE_BRANCH')
call init_attribute_name1 (def_parameter$, geometry$,                     'GEOMETRY')
call init_attribute_name1 (def_parameter$, lattice_type$,                 'LATTICE_TYPE') ! For backwards compatibility
call init_attribute_name1 (def_parameter$, lattice$,                      'LATTICE')
call init_attribute_name1 (def_parameter$, machine$,                      'MACHINE')
call init_attribute_name1 (def_parameter$, taylor_order$,                 'TAYLOR_ORDER')
call init_attribute_name1 (def_parameter$, ran_seed$,                     'RAN_SEED')
call init_attribute_name1 (def_parameter$, n_part$,                       'N_PART')
call init_attribute_name1 (def_parameter$, particle$,                     'PARTICLE')
call init_attribute_name1 (def_parameter$, photon_type$,                  'PHOTON_TYPE')
call init_attribute_name1 (def_parameter$, no_end_marker$,                'NO_END_MARKER')
call init_attribute_name1 (def_parameter$, absolute_time_tracking$,       'ABSOLUTE_TIME_TRACKING')   ! Deprecated
call init_attribute_name1 (def_parameter$, exact_model$,                  'PTC_EXACT_MODEL')    ! Deprecated
call init_attribute_name1 (def_parameter$, exact_misalign$,               'PTC_EXACT_MISALIGN') ! Deprecated
call init_attribute_name1 (def_parameter$, default_tracking_species$,     'DEFAULT_TRACKING_SPECIES')
call init_attribute_name1 (def_parameter$, electric_dipole_moment$,       'ELECTRIC_DIPOLE_MOMENT')
call init_attribute_name1 (def_parameter$, high_energy_space_charge_on$,  'HIGH_ENERGY_SPACE_CHARGE_ON')

call init_attribute_name1 (def_ptc_com$, exact_model$,                    'EXACT_MODEL')
call init_attribute_name1 (def_ptc_com$, exact_misalign$,                 'EXACT_MISALIGN')
call init_attribute_name1 (def_ptc_com$, old_integrator$,                 'OLD_INTEGRATOR')
call init_attribute_name1 (def_ptc_com$, max_fringe_order$,               'MAX_FRINGE_ORDER')
call init_attribute_name1 (def_ptc_com$, vertical_kick$,                  'VERTICAL_KICK')

call init_attribute_name1 (capillary$, l$,                          'L', dependent$)
call init_attribute_name1 (capillary$, n_slice_spline$,             'N_SLICE_SPLINE')
call init_attribute_name1 (capillary$, critical_angle_factor$,      'CRITICAL_ANGLE_FACTOR')

call init_attribute_name1 (converter$, distribution$,               'DISTRIBUTION')
call init_attribute_name1 (converter$, pc_out_min$,                 'PC_OUT_MIN')
call init_attribute_name1 (converter$, pc_out_max$,                 'PC_OUT_MAX')
call init_attribute_name1 (converter$, angle_out_max$,              'ANGLE_OUT_MAX')
call init_attribute_name1 (converter$, species_out$,                'SPECIES_OUT')

call init_attribute_name1 (crystal$, l$,                            'L', dependent$)
call init_attribute_name1 (crystal$, bragg_angle_in$,               'BRAGG_ANGLE_IN', dependent$)
call init_attribute_name1 (crystal$, bragg_angle_out$,              'BRAGG_ANGLE_OUT', dependent$)
call init_attribute_name1 (crystal$, graze_angle_in$,               'GRAZE_ANGLE_IN')
call init_attribute_name1 (crystal$, graze_angle_out$,              'GRAZE_ANGLE_OUT')
call init_attribute_name1 (crystal$, psi_angle$,                    'PSI_ANGLE')
call init_attribute_name1 (crystal$, alpha_angle$,                  'ALPHA_ANGLE', dependent$)
call init_attribute_name1 (crystal$, ref_tilt$,                     'REF_TILT')
call init_attribute_name1 (crystal$, ref_tilt_tot$,                 'REF_TILT_TOT', dependent$)
call init_attribute_name1 (crystal$, tilt_corr$,                    'TILT_CORR', dependent$)
call init_attribute_name1 (crystal$, d_spacing$,                    'D_SPACING', dependent$)
call init_attribute_name1 (crystal$, v_unitcell$,                   'V_UNITCELL', dependent$)
call init_attribute_name1 (crystal$, b_param$,                      'B_PARAM')
call init_attribute_name1 (crystal$, bragg_angle$,                  'BRAGG_ANGLE' , dependent$)
call init_attribute_name1 (crystal$, ref_wavelength$,               'REF_WAVELENGTH', dependent$)
call init_attribute_name1 (crystal$, crystal_type$,                 'CRYSTAL_TYPE')
call init_attribute_name1 (crystal$, thickness$,                    'THICKNESS')
call init_attribute_name1 (crystal$, ref_orbit_follows$,            'REF_ORBIT_FOLLOWS')
call init_attribute_name1 (crystal$, ref_cap_gamma$,                'REF_CAP_GAMMA', dependent$)
call init_attribute_name1 (crystal$, darwin_width_sigma$,           'DARWIN_WIDTH_SIGMA', dependent$)
call init_attribute_name1 (crystal$, darwin_width_pi$,              'DARWIN_WIDTH_PI', dependent$)
call init_attribute_name1 (crystal$, pendellosung_period_sigma$,    'PENDELLOSUNG_PERIOD_SIGMA', dependent$)
call init_attribute_name1 (crystal$, pendellosung_period_pi$,       'PENDELLOSUNG_PERIOD_PI', dependent$)
call init_attribute_name1 (crystal$, dbragg_angle_de$,              'DBRAGG_ANGLE_DE', dependent$)
call init_attribute_name1 (crystal$, is_mosaic$,                    'IS_MOSAIC')
call init_attribute_name1 (crystal$, mosaic_thickness$,             'MOSAIC_THICKNESS')
call init_attribute_name1 (crystal$, mosaic_angle_rms_in_plane$,    'MOSAIC_ANGLE_RMS_IN_PLANE')
call init_attribute_name1 (crystal$, mosaic_angle_rms_out_plane$,   'MOSAIC_ANGLE_RMS_OUT_PLANE')
call init_attribute_name1 (crystal$, mosaic_diffraction_num$,       'MOSAIC_DIFFRACTION_NUM')
call init_attribute_name1 (crystal$, p88$,                          'H_MISALIGN')   ! Defined so H_misalign will show up with type_ele
call init_attribute_name1 (crystal$, curvature$,                    'CURVATURE')
call init_attribute_name1 (crystal$, use_reflectivity_table$,       'USE_REFLECTIVITY_TABLE')
call init_attribute_name1 (crystal$, reflectivity_table$,           'REFLECTIVITY_TABLE')

call init_attribute_name1 (foil$, scatter_method$,                  'SCATTER_METHOD')
call init_attribute_name1 (foil$, scatter_test$,                    'SCATTER_TEST')
call init_attribute_name1 (foil$, thickness$,                       'THICKNESS')
call init_attribute_name1 (foil$, material_type$,                   'MATERIAL_TYPE')
call init_attribute_name1 (foil$, final_charge$,                    'FINAL_CHARGE')
call init_attribute_name1 (foil$, radiation_length$,                'RADIATION_LENGTH')
call init_attribute_name1 (foil$, radiation_length_used$,           'RADIATION_LENGTH_USED', dependent$)
call init_attribute_name1 (foil$, density$,                         'DENSITY')
call init_attribute_name1 (foil$, density_used$,                    'DENSITY_USED', dependent$)
call init_attribute_name1 (foil$, area_density$,                    'AREA_DENSITY', quasi_free$)
call init_attribute_name1 (foil$, area_density_used$,               'AREA_DENSITY_USED', dependent$)
call init_attribute_name1 (foil$, x1_edge$,                         'X1_EDGE')
call init_attribute_name1 (foil$, x2_edge$,                         'X2_EDGE')
call init_attribute_name1 (foil$, y1_edge$,                         'Y1_EDGE')
call init_attribute_name1 (foil$, y2_edge$,                         'Y2_EDGE')
call init_attribute_name1 (foil$, dthickness_dx$,                   'DTHICKNESS_DX')
call init_attribute_name1 (foil$, f_factor$,                        'F_FACTOR')
call init_attribute_name1 (foil$, num_steps$,                       'NUM_STEPS')

!! call init_attribute_name1 (foil$, probability_final_charge$         'PROBABILITY_FINAL_CHARGE')

call init_attribute_name1 (lens$, l$,                               'L')
call init_attribute_name1 (lens$, radius$,                          'RADIUS')
call init_attribute_name1 (lens$, focal_strength$,                  'FOCAL_STRENGTH')

call init_attribute_name1 (detector$, l$,                               'L', dependent$)
call init_attribute_name1 (detector$, pixel$,                           'PIXEL')
call init_attribute_name1 (detector$, curvature$,                       'CURVATURE')

call init_attribute_name1 (diffraction_plate$, l$,                      'l', private$)
call init_attribute_name1 (diffraction_plate$, mode$,                   'MODE')
call init_attribute_name1 (diffraction_plate$, field_scale_factor$,     'FIELD_SCALE_FACTOR')
call init_attribute_name1 (diffraction_plate$, ref_wavelength$,         'REF_WAVELENGTH', dependent$)
call init_attribute_name1 (diffraction_plate$, curvature$,              'CURVATURE')

call init_attribute_name1 (mask$, l$,                                   'l', private$)
call init_attribute_name1 (mask$, mode$,                                'MODE')
call init_attribute_name1 (mask$, field_scale_factor$,                  'FIELD_SCALE_FACTOR')
call init_attribute_name1 (mask$, ref_wavelength$,                      'REF_WAVELENGTH', dependent$)

call init_attribute_name1 (drift$, spin_fringe_on$,                 'spin_fringe_on', private$)
call init_attribute_name1 (drift$, fringe_type$,                    'fringe_type', private$)
call init_attribute_name1 (drift$, fringe_at$,                      'fringe_at', private$)
call init_attribute_name1 (drift$, split_id$,                       'split_id', private$)

call init_attribute_name1 (e_gun$, dt_max$,                         'DT_MAX')
call init_attribute_name1 (e_gun$, emit_fraction$,                  'EMIT_FRACTION')
call init_attribute_name1 (e_gun$, e_tot_ref_init$,                 'e_tot_ref_init', private$)
call init_attribute_name1 (e_gun$, p0c_ref_init$,                   'p0c_ref_init', private$)
call init_attribute_name1 (e_gun$, field_autoscale$,                'FIELD_AUTOSCALE', quasi_free$)
call init_attribute_name1 (e_gun$, autoscale_amplitude$,            'AUTOSCALE_AMPLITUDE')
call init_attribute_name1 (e_gun$, autoscale_phase$,                'AUTOSCALE_PHASE')
call init_attribute_name1 (e_gun$, voltage$,                        'VOLTAGE')
call init_attribute_name1 (e_gun$, voltage_err$,                    'VOLTAGE_ERR')
call init_attribute_name1 (e_gun$, gradient$,                       'GRADIENT')
call init_attribute_name1 (e_gun$, gradient_err$,                   'GRADIENT_ERR')
call init_attribute_name1 (e_gun$, cartesian_map$,                  'CARTESIAN_MAP')
call init_attribute_name1 (e_gun$, cylindrical_map$,                'CYLINDRICAL_MAP')
call init_attribute_name1 (e_gun$, gen_grad_map$,                   'GEN_GRAD_MAP')
call init_attribute_name1 (e_gun$, grid_field$,                     'GRID_FIELD')
call init_attribute_name1 (e_gun$, rf_frequency$,                   'RF_FREQUENCY')
call init_attribute_name1 (e_gun$, rf_wavelength$,                  'RF_WAVELENGTH', dependent$)
call init_attribute_name1 (e_gun$, rf_clock_harmonic$,              'rf_clock_harminic', private$)
call init_attribute_name1 (e_gun$, phi0$,                           'PHI0')
call init_attribute_name1 (e_gun$, phi0_err$,                       'PHI0_ERR')
! e_gun attribute phi0_multipass should always be 0 and is used to make lcavity and e_gun equations similar
call init_attribute_name1 (e_gun$, phi0_multipass$,                 'phi0_multipass', private$) 
call init_attribute_name1 (e_gun$, phi0_autoscale$,                 'PHI0_AUTOSCALE', quasi_free$)
call init_attribute_name1 (e_gun$, voltage_tot$,                    'VOLTAGE_TOT', dependent$)
call init_attribute_name1 (e_gun$, gradient_tot$,                   'GRADIENT_TOT', dependent$)

call init_attribute_name1 (ecollimator$, px_aperture_width2$,       'PX_APERTURE_WIDTH2')
call init_attribute_name1 (ecollimator$, px_aperture_center$,       'PX_APERTURE_CENTER')
call init_attribute_name1 (ecollimator$, py_aperture_width2$,       'PY_APERTURE_WIDTH2')
call init_attribute_name1 (ecollimator$, py_aperture_center$,       'PY_APERTURE_CENTER')
call init_attribute_name1 (ecollimator$, pz_aperture_width2$,       'PZ_APERTURE_WIDTH2')
call init_attribute_name1 (ecollimator$, pz_aperture_center$,       'PZ_APERTURE_CENTER')
call init_attribute_name1 (ecollimator$, z_aperture_width2$,        'Z_APERTURE_WIDTH2')
call init_attribute_name1 (ecollimator$, z_aperture_center$,        'Z_APERTURE_CENTER')

call init_attribute_name1 (elseparator$, gap$,                      'GAP')
call init_attribute_name1 (elseparator$, e_field$,                  'E_FIELD', quasi_free$)
call init_attribute_name1 (elseparator$, voltage$,                  'VOLTAGE', quasi_free$)
call init_attribute_name1 (elseparator$, voltage_err$,              'voltage_err', private$) ! To fool master_parameter_value
call init_attribute_name1 (elseparator$, r0_mag$,                   'R0_MAG')
call init_attribute_name1 (elseparator$, r0_elec$,                  'R0_ELEC')
call init_attribute_name1 (elseparator$, field_master$,             'FIELD_MASTER')
call init_attribute_name1 (elseparator$, cartesian_map$,            'CARTESIAN_MAP')
call init_attribute_name1 (elseparator$, cylindrical_map$,          'CYLINDRICAL_MAP')
call init_attribute_name1 (elseparator$, gen_grad_map$,             'GEN_GRAD_MAP')
call init_attribute_name1 (elseparator$, grid_field$,               'GRID_FIELD')
call init_attribute_name1 (elseparator$, ptc_canonical_coords$,     'PTC_CANONICAL_COORDS')

call init_attribute_name1 (rcollimator$, px_aperture_width2$,       'PX_APERTURE_WIDTH2')
call init_attribute_name1 (rcollimator$, px_aperture_center$,       'PX_APERTURE_CENTER')
call init_attribute_name1 (rcollimator$, py_aperture_width2$,       'PY_APERTURE_WIDTH2')
call init_attribute_name1 (rcollimator$, py_aperture_center$,       'PY_APERTURE_CENTER')
call init_attribute_name1 (rcollimator$, pz_aperture_width2$,       'PZ_APERTURE_WIDTH2')
call init_attribute_name1 (rcollimator$, pz_aperture_center$,       'PZ_APERTURE_CENTER')
call init_attribute_name1 (rcollimator$, z_aperture_width2$,        'Z_APERTURE_WIDTH2')
call init_attribute_name1 (rcollimator$, z_aperture_center$,        'Z_APERTURE_CENTER')

call init_attribute_name1 (em_field$, cartesian_map$,               'CARTESIAN_MAP')
call init_attribute_name1 (em_field$, cylindrical_map$,             'CYLINDRICAL_MAP')
call init_attribute_name1 (em_field$, gen_grad_map$,                'GEN_GRAD_MAP')
call init_attribute_name1 (em_field$, grid_field$,                  'GRID_FIELD')
call init_attribute_name1 (em_field$, ptc_canonical_coords$,        'PTC_CANONICAL_COORDS')
call init_attribute_name1 (em_field$, rf_frequency$,                'RF_FREQUENCY')
call init_attribute_name1 (em_field$, rf_wavelength$,               'RF_WAVELENGTH', dependent$)
call init_attribute_name1 (em_field$, rf_clock_harmonic$,           'rf_clock_harminic', private$)
call init_attribute_name1 (em_field$, field_autoscale$,             'FIELD_AUTOSCALE', quasi_free$)
call init_attribute_name1 (em_field$, phi0_autoscale$,              'PHI0_AUTOSCALE', quasi_free$)
call init_attribute_name1 (em_field$, autoscale_amplitude$,         'AUTOSCALE_AMPLITUDE')
call init_attribute_name1 (em_field$, autoscale_phase$,             'AUTOSCALE_PHASE')
call init_attribute_name1 (em_field$, phi0$,                        'PHI0')
call init_attribute_name1 (em_field$, phi0_err$,                    'PHI0_ERR')
call init_attribute_name1 (em_field$, constant_ref_energy$,         'CONSTANT_REF_ENERGY')
call init_attribute_name1 (em_field$, polarity$,                    'POLARITY')


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
call init_attribute_name1 (overlay$, gang$,                         'GANG')
call init_attribute_name1 (overlay$, x_knot$,                       'X_KNOT')
call init_attribute_name1 (overlay$, y_knot$,                       'Y_KNOT')
call init_attribute_name1 (overlay$, slave$,                        'SLAVE')
call init_attribute_name1 (overlay$, is_on$,                        'IS_ON')
call init_attribute_name1 (overlay$, interpolation$,                'INTERPOLATION')

call init_attribute_name1 (ramper$, var$,                          'VAR')
call init_attribute_name1 (ramper$, x_knot$,                       'X_KNOT')
call init_attribute_name1 (ramper$, y_knot$,                       'Y_KNOT')
call init_attribute_name1 (ramper$, slave$,                        'SLAVE')
call init_attribute_name1 (ramper$, is_on$,                        'IS_ON')
call init_attribute_name1 (ramper$, interpolation$,                'INTERPOLATION')

call init_attribute_name1 (group$, var$,                            'VAR')
call init_attribute_name1 (group$, gang$,                           'GANG')
call init_attribute_name1 (group$, start_edge$,                     'START_EDGE')
call init_attribute_name1 (group$, end_edge$,                       'END_EDGE')
call init_attribute_name1 (group$, accordion_edge$,                 'ACCORDION_EDGE')
call init_attribute_name1 (group$, s_position$,                     'S_POSITION')
call init_attribute_name1 (group$, x_knot$,                         'X_KNOT')
call init_attribute_name1 (group$, y_knot$,                         'Y_KNOT')
call init_attribute_name1 (group$, slave$,                          'SLAVE')
call init_attribute_name1 (group$, is_on$,                          'IS_ON')
call init_attribute_name1 (group$, interpolation$,                  'INTERPOLATION')

call init_attribute_name1 (lcavity$, longitudinal_mode$,            'LONGITUDINAL_MODE')
call init_attribute_name1 (lcavity$, field_autoscale$,              'FIELD_AUTOSCALE', quasi_free$)
call init_attribute_name1 (lcavity$, autoscale_amplitude$,          'AUTOSCALE_AMPLITUDE')
call init_attribute_name1 (lcavity$, autoscale_phase$,              'AUTOSCALE_PHASE')
call init_attribute_name1 (lcavity$, cavity_type$,                  'CAVITY_TYPE')
call init_attribute_name1 (lcavity$, phi0_multipass$,               'PHI0_MULTIPASS')
call init_attribute_name1 (lcavity$, phi0$,                         'PHI0')
call init_attribute_name1 (lcavity$, gradient$,                     'GRADIENT')
call init_attribute_name1 (lcavity$, rf_frequency$,                 'RF_FREQUENCY')
call init_attribute_name1 (lcavity$, rf_wavelength$,                'RF_WAVELENGTH', dependent$)
call init_attribute_name1 (lcavity$, rf_clock_harmonic$,            'rf_clock_harminic', private$)
call init_attribute_name1 (lcavity$, e_loss$,                       'E_LOSS')
call init_attribute_name1 (lcavity$, voltage$,                      'VOLTAGE', quasi_free$)
call init_attribute_name1 (lcavity$, field_master$,                 'FIELD_MASTER')
call init_attribute_name1 (lcavity$, coupler_strength$,             'COUPLER_STRENGTH')
call init_attribute_name1 (lcavity$, coupler_angle$,                'COUPLER_ANGLE')
call init_attribute_name1 (lcavity$, coupler_phase$,                'COUPLER_PHASE')
call init_attribute_name1 (lcavity$, coupler_at$,                   'COUPLER_AT')
call init_attribute_name1 (lcavity$, gradient_err$,                 'GRADIENT_ERR')
call init_attribute_name1 (lcavity$, voltage_err$,                  'VOLTAGE_ERR', quasi_free$)
call init_attribute_name1 (lcavity$, phi0_err$,                     'PHI0_ERR')
call init_attribute_name1 (lcavity$, cartesian_map$,                'CARTESIAN_MAP')
call init_attribute_name1 (lcavity$, cylindrical_map$,              'CYLINDRICAL_MAP')
call init_attribute_name1 (lcavity$, gen_grad_map$,                 'GEN_GRAD_MAP')
call init_attribute_name1 (lcavity$, grid_field$,                   'GRID_FIELD')
call init_attribute_name1 (lcavity$, phi0_autoscale$,               'PHI0_AUTOSCALE', quasi_free$)
call init_attribute_name1 (lcavity$, n_cell$,                       'N_CELL')
call init_attribute_name1 (lcavity$, l_active$,                     'L_ACTIVE', dependent$)
call init_attribute_name1 (lcavity$, voltage_tot$,                  'VOLTAGE_TOT', dependent$)
call init_attribute_name1 (lcavity$, gradient_tot$,                 'GRADIENT_TOT', dependent$)

call init_attribute_name1 (marker$, l$,                             'L', dependent$)
call init_attribute_name1 (marker$, e_tot_ref_init$,                'e_tot_ref_init', private$)
call init_attribute_name1 (marker$, p0c_ref_init$,                  'p0c_ref_init', private$)
call init_attribute_name1 (marker$, sr_wake$,                       'SR_WAKE')
call init_attribute_name1 (marker$, lr_wake$,                       'LR_WAKE')
call init_attribute_name1 (marker$, sr_wake_file$,                  'SR_WAKE_FILE')
call init_attribute_name1 (marker$, lr_wake_file$,                  'LR_WAKE_FILE')
call init_attribute_name1 (marker$, lr_freq_spread$,                'LR_FREQ_SPREAD')
call init_attribute_name1 (marker$, lr_self_wake_on$,               'LR_SELF_WAKE_ON')
call init_attribute_name1 (marker$, ref_species$,                   'REF_SPECIES', dependent$)

call init_attribute_name1 (match$, l$,                              'L')
call init_attribute_name1 (match$, delta_time$,                     'DELTA_TIME')
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
call init_attribute_name1 (match$, matrix$,                         'MATRIX')
call init_attribute_name1 (match$, kick0$,                          'KICK0')
call init_attribute_name1 (match$, recalc$,                         'RECALC')
call init_attribute_name1 (match$, is_on$,                          'IS_ON')
call init_attribute_name1 (match$, c11_mat0$,                       'C11_MAT0')
call init_attribute_name1 (match$, c12_mat0$,                       'C12_MAT0')
call init_attribute_name1 (match$, c21_mat0$,                       'C21_MAT0')
call init_attribute_name1 (match$, c22_mat0$,                       'C22_MAT0')
call init_attribute_name1 (match$, mode_flip0$,                     'MODE_FLIP0')
call init_attribute_name1 (match$, c11_mat1$,                       'C11_MAT1')
call init_attribute_name1 (match$, c12_mat1$,                       'C12_MAT1')
call init_attribute_name1 (match$, c21_mat1$,                       'C21_MAT1')
call init_attribute_name1 (match$, c22_mat1$,                       'C22_MAT1')
call init_attribute_name1 (match$, mode_flip1$,                     'MODE_FLIP1')

call init_attribute_name1 (pickup$, num_steps$,                     'NUM_STEPS')
call init_attribute_name1 (pickup$, ds_step$,                       'DS_STEP')

attrib_array(instrument$, :)                         = attrib_array(monitor$, :)
attrib_array(pipe$, :)                               = attrib_array(monitor$, :)

call init_attribute_name1 (gkicker$, x_kick$,                       'X_KICK')
call init_attribute_name1 (gkicker$, px_kick$,                      'PX_KICK')
call init_attribute_name1 (gkicker$, y_kick$,                       'Y_KICK')
call init_attribute_name1 (gkicker$, py_kick$,                      'PY_KICK')
call init_attribute_name1 (gkicker$, z_kick$,                       'Z_KICK')
call init_attribute_name1 (gkicker$, pz_kick$,                      'PZ_KICK')

call init_attribute_name1 (hkicker$, kick$,                         'KICK', quasi_free$)
call init_attribute_name1 (hkicker$, field_master$,                 'FIELD_MASTER')
call init_attribute_name1 (hkicker$, bl_kick$,                      'BL_KICK', quasi_free$)
call init_attribute_name1 (hkicker$, cartesian_map$,                'CARTESIAN_MAP')
call init_attribute_name1 (hkicker$, cylindrical_map$,              'CYLINDRICAL_MAP')
call init_attribute_name1 (hkicker$, gen_grad_map$,                 'GEN_GRAD_MAP')
call init_attribute_name1 (hkicker$, grid_field$,                   'GRID_FIELD')
call init_attribute_name1 (hkicker$, ptc_canonical_coords$,         'PTC_CANONICAL_COORDS')

attrib_array(vkicker$, :) = attrib_array(hkicker$, :)

call init_attribute_name1 (kicker$, h_displace$,                    'H_DISPLACE')
call init_attribute_name1 (kicker$, v_displace$,                    'V_DISPLACE')
call init_attribute_name1 (kicker$, r0_mag$,                        'R0_MAG')
call init_attribute_name1 (kicker$, r0_elec$,                       'R0_ELEC')
call init_attribute_name1 (kicker$, field_master$,                  'FIELD_MASTER')
call init_attribute_name1 (kicker$, cartesian_map$,                 'CARTESIAN_MAP')
call init_attribute_name1 (kicker$, cylindrical_map$,               'CYLINDRICAL_MAP')
call init_attribute_name1 (kicker$, gen_grad_map$,                  'GEN_GRAD_MAP')
call init_attribute_name1 (kicker$, grid_field$,                    'GRID_FIELD')
call init_attribute_name1 (kicker$, ptc_canonical_coords$,          'PTC_CANONICAL_COORDS')

call init_attribute_name1 (ac_kicker$, interpolation$,                 'INTERPOLATION')
call init_attribute_name1 (ac_kicker$, r0_mag$,                        'R0_MAG')
call init_attribute_name1 (ac_kicker$, r0_elec$,                       'R0_ELEC')
call init_attribute_name1 (ac_kicker$, field_master$,                  'FIELD_MASTER')
call init_attribute_name1 (ac_kicker$, cartesian_map$,                 'CARTESIAN_MAP')
call init_attribute_name1 (ac_kicker$, cylindrical_map$,               'CYLINDRICAL_MAP')
call init_attribute_name1 (ac_kicker$, gen_grad_map$,                  'GEN_GRAD_MAP')
call init_attribute_name1 (ac_kicker$, grid_field$,                    'GRID_FIELD')
call init_attribute_name1 (ac_kicker$, amp_vs_time$,                   'AMP_VS_TIME')
call init_attribute_name1 (ac_kicker$, frequencies$,                   'FREQUENCIES')
call init_attribute_name1 (ac_kicker$, t_offset$,                      'T_OFFSET')
call init_attribute_name1 (ac_kicker$, phi0_multipass$,                'PHI0_MULTIPASS')

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
call init_attribute_name1 (custom$, field_master$,                  'FIELD_MASTER')
call init_attribute_name1 (custom$, delta_e_ref$,                   'DELTA_E_REF')

call init_attribute_name1 (floor_shift$, l$,                        'L')
call init_attribute_name1 (floor_shift$, origin_ele$,               'ORIGIN_ELE')
call init_attribute_name1 (floor_shift$, origin_ele_ref_pt$,        'ORIGIN_ELE_REF_PT')
call init_attribute_name1 (floor_shift$, upstream_coord_dir$,         'UPSTREAM_ELE_DIR', dependent$)
call init_attribute_name1 (floor_shift$, downstream_coord_dir$,       'DOWNSTREAM_ELE_DIR', dependent$)

call init_attribute_name1 (fiducial$, l$,                           'L', dependent$)
call init_attribute_name1 (fiducial$, origin_ele$,                  'ORIGIN_ELE')
call init_attribute_name1 (fiducial$, origin_ele_ref_pt$,           'ORIGIN_ELE_REF_PT')
call init_attribute_name1 (fiducial$, dx_origin$,                   'DX_ORIGIN')
call init_attribute_name1 (fiducial$, dy_origin$,                   'DY_ORIGIN')
call init_attribute_name1 (fiducial$, dz_origin$,                   'DZ_ORIGIN')
call init_attribute_name1 (fiducial$, dtheta_origin$,               'DTHETA_ORIGIN')
call init_attribute_name1 (fiducial$, dphi_origin$,                 'DPHI_ORIGIN')
call init_attribute_name1 (fiducial$, dpsi_origin$,                 'DPSI_ORIGIN')

call init_attribute_name1 (null_ele$, ix_branch$,                   'ix_branch', private$)

call init_attribute_name1 (quadrupole$, k1$,                        'K1', quasi_free$)
call init_attribute_name1 (quadrupole$, B1_gradient$,               'B1_GRADIENT', quasi_free$)
call init_attribute_name1 (quadrupole$, r0_mag$,                    'R0_MAG')
call init_attribute_name1 (quadrupole$, r0_elec$,                   'R0_ELEC')
call init_attribute_name1 (quadrupole$, field_master$,              'FIELD_MASTER')
call init_attribute_name1 (quadrupole$, cartesian_map$,             'CARTESIAN_MAP')
call init_attribute_name1 (quadrupole$, cylindrical_map$,           'CYLINDRICAL_MAP')
call init_attribute_name1 (quadrupole$, gen_grad_map$,              'GEN_GRAD_MAP')
call init_attribute_name1 (quadrupole$, grid_field$,                'GRID_FIELD')
call init_attribute_name1 (quadrupole$, ptc_canonical_coords$,      'PTC_CANONICAL_COORDS')
call init_attribute_name1 (quadrupole$, fq1$,                       'FQ1')
call init_attribute_name1 (quadrupole$, fq2$,                       'FQ2')

call init_attribute_name1 (sextupole$, k2$,                         'K2', quasi_free$)
call init_attribute_name1 (sextupole$, B2_gradient$,                'B2_GRADIENT', quasi_free$)
call init_attribute_name1 (sextupole$, r0_mag$,                     'R0_MAG')
call init_attribute_name1 (sextupole$, r0_elec$,                    'R0_ELEC')
call init_attribute_name1 (sextupole$, field_master$,               'FIELD_MASTER')
call init_attribute_name1 (sextupole$, cartesian_map$,              'CARTESIAN_MAP')
call init_attribute_name1 (sextupole$, cylindrical_map$,            'CYLINDRICAL_MAP')
call init_attribute_name1 (sextupole$, gen_grad_map$,               'GEN_GRAD_MAP')
call init_attribute_name1 (sextupole$, grid_field$,                 'GRID_FIELD')
call init_attribute_name1 (sextupole$, ptc_canonical_coords$,       'PTC_CANONICAL_COORDS')

call init_attribute_name1 (octupole$, k3$,                          'K3', quasi_free$)
call init_attribute_name1 (octupole$, B3_gradient$,                 'B3_GRADIENT', quasi_free$)
call init_attribute_name1 (octupole$, r0_mag$,                      'R0_MAG')
call init_attribute_name1 (octupole$, r0_elec$,                     'R0_ELEC')
call init_attribute_name1 (octupole$, field_master$,                'FIELD_MASTER')
call init_attribute_name1 (octupole$, cartesian_map$,               'CARTESIAN_MAP')
call init_attribute_name1 (octupole$, cylindrical_map$,             'CYLINDRICAL_MAP')
call init_attribute_name1 (octupole$, gen_grad_map$,                'GEN_GRAD_MAP')
call init_attribute_name1 (octupole$, grid_field$,                  'GRID_FIELD')
call init_attribute_name1 (octupole$, ptc_canonical_coords$,        'PTC_CANONICAL_COORDS')

call init_attribute_name1 (thick_multipole$, field_master$,         'FIELD_MASTER')
call init_attribute_name1 (thick_multipole$, cartesian_map$,        'CARTESIAN_MAP')
call init_attribute_name1 (thick_multipole$, cylindrical_map$,      'CYLINDRICAL_MAP')
call init_attribute_name1 (thick_multipole$, gen_grad_map$,         'GEN_GRAD_MAP')
call init_attribute_name1 (thick_multipole$, grid_field$,           'GRID_FIELD')
call init_attribute_name1 (thick_multipole$, ptc_canonical_coords$, 'PTC_CANONICAL_COORDS')

call init_attribute_name1 (patch$, l$,                              'L', quasi_free$)
call init_attribute_name1 (patch$, user_sets_length$,               'USER_SETS_LENGTH')
call init_attribute_name1 (patch$, t_offset$,                       'T_OFFSET')
call init_attribute_name1 (patch$, e_tot_set$,                      'E_TOT_SET')
call init_attribute_name1 (patch$, p0c_set$,                        'P0C_SET')
call init_attribute_name1 (patch$, e_tot_offset$,                   'E_TOT_OFFSET')
call init_attribute_name1 (patch$, flexible$,                       'FLEXIBLE')
call init_attribute_name1 (patch$, upstream_coord_dir$,             'UPSTREAM_ELE_DIR', dependent$)
call init_attribute_name1 (patch$, downstream_coord_dir$,           'DOWNSTREAM_ELE_DIR', dependent$)
call init_attribute_name1 (patch$, ref_coords$,                     'REF_COORDS')
call init_attribute_name1 (patch$, field_calc$,                     'FIELD_CALC')
call init_attribute_name1 (patch$, csr_method$,                     'CSR_METHOD')
call init_attribute_name1 (patch$, csr_ds_step$,                    'CSR_DS_STEP')
call init_attribute_name1 (patch$, space_charge_method$,            'SPACE_CHARGE_METHOD')

call init_attribute_name1 (crab_cavity$, voltage$,                  'VOLTAGE')
call init_attribute_name1 (crab_cavity$, phi0$,                     'PHI0')
call init_attribute_name1 (crab_cavity$, phi0_multipass$,           'PHI0_MULTIPASS')
call init_attribute_name1 (crab_cavity$, harmon$,                   'HARMON')
call init_attribute_name1 (crab_cavity$, harmon_master$,            'HARMON_MASTER')
call init_attribute_name1 (crab_cavity$, cartesian_map$,            'CARTESIAN_MAP')
call init_attribute_name1 (crab_cavity$, cylindrical_map$,          'CYLINDRICAL_MAP')
call init_attribute_name1 (crab_cavity$, gen_grad_map$,             'GEN_GRAD_MAP')
call init_attribute_name1 (crab_cavity$, grid_field$,               'GRID_FIELD')
call init_attribute_name1 (crab_cavity$, gradient$,                 'GRADIENT', dependent$)
call init_attribute_name1 (crab_cavity$, rf_frequency$,             'RF_FREQUENCY')
call init_attribute_name1 (crab_cavity$, rf_wavelength$,            'RF_WAVELENGTH', dependent$)
call init_attribute_name1 (crab_cavity$, rf_clock_harmonic$,        'rf_clock_harminic', private$)
call init_attribute_name1 (crab_cavity$, field_autoscale$,          'FIELD_AUTOSCALE', private$)      ! Not yet used
call init_attribute_name1 (crab_cavity$, phi0_autoscale$,           'PHI0_AUTOSCALE', private$)       ! Not yet used
call init_attribute_name1 (crab_cavity$, autoscale_amplitude$,      'AUTOSCALE_AMPLITUDE', private$)  ! Not yet used
call init_attribute_name1 (crab_cavity$, autoscale_phase$,          'AUTOSCALE_PHASE', private$)      ! Not yet used
call init_attribute_name1 (crab_cavity$, field_master$,             'FIELD_MASTER')

call init_attribute_name1 (rfcavity$, longitudinal_mode$,           'LONGITUDINAL_MODE')
call init_attribute_name1 (rfcavity$, field_autoscale$,             'FIELD_AUTOSCALE', quasi_free$)
call init_attribute_name1 (rfcavity$, phi0_autoscale$,              'PHI0_AUTOSCALE', quasi_free$)
call init_attribute_name1 (rfcavity$, autoscale_amplitude$,         'AUTOSCALE_AMPLITUDE')
call init_attribute_name1 (rfcavity$, autoscale_phase$,             'AUTOSCALE_PHASE')
call init_attribute_name1 (rfcavity$, cavity_type$,                 'CAVITY_TYPE')
call init_attribute_name1 (rfcavity$, voltage$,                     'VOLTAGE')
call init_attribute_name1 (rfcavity$, rf_frequency$,                'RF_FREQUENCY', quasi_free$)
call init_attribute_name1 (rfcavity$, rf_wavelength$,               'RF_WAVELENGTH', dependent$)
call init_attribute_name1 (rfcavity$, rf_clock_harmonic$,           'rf_clock_harminic', private$)
call init_attribute_name1 (rfcavity$, phi0_multipass$,              'PHI0_MULTIPASS')
call init_attribute_name1 (rfcavity$, phi0$,                        'PHI0')
call init_attribute_name1 (rfcavity$, harmon$,                      'HARMON', quasi_free$)
call init_attribute_name1 (rfcavity$, harmon_master$,               'HARMON_MASTER')
call init_attribute_name1 (rfcavity$, coupler_strength$,            'COUPLER_STRENGTH')
call init_attribute_name1 (rfcavity$, coupler_angle$,               'COUPLER_ANGLE')
call init_attribute_name1 (rfcavity$, coupler_phase$,               'COUPLER_PHASE')
call init_attribute_name1 (rfcavity$, coupler_at$,                  'COUPLER_AT')
call init_attribute_name1 (rfcavity$, cartesian_map$,               'CARTESIAN_MAP')
call init_attribute_name1 (rfcavity$, cylindrical_map$,             'CYLINDRICAL_MAP')
call init_attribute_name1 (rfcavity$, gen_grad_map$,                'GEN_GRAD_MAP')
call init_attribute_name1 (rfcavity$, grid_field$,                  'GRID_FIELD')
call init_attribute_name1 (rfcavity$, n_cell$,                      'N_CELL')
call init_attribute_name1 (rfcavity$, phi0_max$,                    'phi0_max', private$)
call init_attribute_name1 (rfcavity$, phi0_err$,                    'phi0_err', private$)
call init_attribute_name1 (rfcavity$, gradient$,                    'GRADIENT', dependent$)
call init_attribute_name1 (rfcavity$, gradient_err$,                'gradient_err', private$)
call init_attribute_name1 (rfcavity$, voltage_err$,                 'voltage_err', private$)
call init_attribute_name1 (rfcavity$, l_active$,                    'L_ACTIVE', dependent$)

call init_attribute_name1 (rf_bend$, angle$,                        'ANGLE', quasi_free$)
call init_attribute_name1 (rf_bend$, ref_tilt$,                     'REF_TILT')
call init_attribute_name1 (rf_bend$, ref_tilt_tot$,                 'REF_TILT_TOT', dependent$)
call init_attribute_name1 (rf_bend$, g$,                            'G', quasi_free$)
call init_attribute_name1 (rf_bend$, roll$,                         'ROLL', override = .true.)
call init_attribute_name1 (rf_bend$, roll_tot$,                     'ROLL_TOT', dependent$, override = .true.)
call init_attribute_name1 (rf_bend$, rho$,                          'RHO', quasi_free$)
call init_attribute_name1 (rf_bend$, l_chord$,                      'L_CHORD', quasi_free$)
call init_attribute_name1 (rf_bend$, l_sagitta$,                    'L_SAGITTA', dependent$)
call init_attribute_name1 (rf_bend$, l_rectangle$,                  'L_RECTANGLE', quasi_free$)
call init_attribute_name1 (rf_bend$, fiducial_pt$,                  'FIDUCIAL_PT')
call init_attribute_name1 (rf_bend$, b_field$,                      'B_FIELD', quasi_free$)
call init_attribute_name1 (rf_bend$, init_needed$,                  'init_needed', private$)
call init_attribute_name1 (rf_bend$, field_master$,                 'FIELD_MASTER')
call init_attribute_name1 (rf_bend$, grid_field$,                   'GRID_FIELD')
call init_attribute_name1 (rf_bend$, rf_frequency$,                 'RF_FREQUENCY', quasi_free$)
call init_attribute_name1 (rf_bend$, rf_wavelength$,                'RF_WAVELENGTH', dependent$)
call init_attribute_name1 (rf_bend$, rf_clock_harmonic$,            'rf_clock_harminic', private$)
call init_attribute_name1 (rf_bend$, phi0_multipass$,               'PHI0_MULTIPASS')
call init_attribute_name1 (rf_bend$, phi0$,                         'PHI0')
call init_attribute_name1 (rf_bend$, harmon$,                       'HARMON', quasi_free$)
call init_attribute_name1 (rf_bend$, harmon_master$,                'HARMON_MASTER')
call init_attribute_name1 (rf_bend$, dg$,                           'DG', private$)
call init_attribute_name1 (rf_bend$, db_field$,                     'DB_FIELD', private$)

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
call init_attribute_name1 (sbend$, dg$,                             'DG', quasi_free$)
call init_attribute_name1 (sbend$, roll$,                           'ROLL', override = .true.)
call init_attribute_name1 (sbend$, roll_tot$,                       'ROLL_TOT', dependent$, override = .true.)
call init_attribute_name1 (sbend$, hgap$,                           'HGAP')
call init_attribute_name1 (sbend$, hgapx$,                          'HGAPX')
call init_attribute_name1 (sbend$, fint$,                           'FINT')
call init_attribute_name1 (sbend$, fintx$,                          'FINTX')
call init_attribute_name1 (sbend$, rho$,                            'RHO', quasi_free$)
call init_attribute_name1 (sbend$, init_needed$,                    'init_needed', private$)
call init_attribute_name1 (sbend$, l_chord$,                        'L_CHORD', quasi_free$)
call init_attribute_name1 (sbend$, l_sagitta$,                      'L_SAGITTA', dependent$)
call init_attribute_name1 (sbend$, l_rectangle$,                    'L_RECTANGLE', quasi_free$)
call init_attribute_name1 (sbend$, fiducial_pt$,                    'FIDUCIAL_PT')
call init_attribute_name1 (sbend$, ptc_fringe_geometry$,            'PTC_FRINGE_GEOMETRY')
call init_attribute_name1 (sbend$, b_field$,                        'B_FIELD', quasi_free$)
call init_attribute_name1 (sbend$, db_field$,                       'DB_FIELD', quasi_free$)
call init_attribute_name1 (sbend$, b1_gradient$,                    'B1_GRADIENT', quasi_free$)
call init_attribute_name1 (sbend$, b2_gradient$,                    'B2_GRADIENT', quasi_free$)
call init_attribute_name1 (sbend$, r0_mag$,                         'R0_MAG')
call init_attribute_name1 (sbend$, r0_elec$,                        'R0_ELEC')
call init_attribute_name1 (sbend$, field_master$,                   'FIELD_MASTER')
call init_attribute_name1 (sbend$, cartesian_map$,                  'CARTESIAN_MAP')
call init_attribute_name1 (sbend$, cylindrical_map$,                'CYLINDRICAL_MAP')
call init_attribute_name1 (sbend$, gen_grad_map$,                   'GEN_GRAD_MAP')
call init_attribute_name1 (sbend$, grid_field$,                     'GRID_FIELD')
call init_attribute_name1 (sbend$, ptc_canonical_coords$,           'PTC_CANONICAL_COORDS')
call init_attribute_name1 (sbend$, exact_multipoles$,               'EXACT_MULTIPOLES')
call init_attribute_name1 (sbend$, ptc_field_geometry$,             'PTC_FIELD_GEOMETRY')
call init_attribute_name1 (sbend$, g_tot$,                          'G_TOT', dependent$)
call init_attribute_name1 (sbend$, b_field_tot$,                    'B_FIELD_TOT', dependent$)

attrib_array(rbend$, :) = attrib_array(sbend$, :)

call init_attribute_name1 (solenoid$, ks$,                          'KS', quasi_free$)
call init_attribute_name1 (solenoid$, bs_field$,                    'BS_FIELD', quasi_free$)
call init_attribute_name1 (solenoid$, r0_mag$,                      'R0_MAG')
call init_attribute_name1 (solenoid$, r0_elec$,                     'R0_ELEC')
call init_attribute_name1 (solenoid$, field_master$,                'FIELD_MASTER')
call init_attribute_name1 (solenoid$, cartesian_map$,               'CARTESIAN_MAP')
call init_attribute_name1 (solenoid$, cylindrical_map$,             'CYLINDRICAL_MAP')
call init_attribute_name1 (solenoid$, gen_grad_map$,                'GEN_GRAD_MAP')
call init_attribute_name1 (solenoid$, grid_field$,                  'GRID_FIELD')
call init_attribute_name1 (solenoid$, ptc_canonical_coords$,        'PTC_CANONICAL_COORDS')
call init_attribute_name1 (solenoid$, l_soft_edge$,                 'L_SOFT_EDGE')
call init_attribute_name1 (solenoid$, r_solenoid$,                  'R_SOLENOID')

call init_attribute_name1 (sample$, l$,                             'L')
call init_attribute_name1 (sample$, mode$,                          'MODE')
call init_attribute_name1 (sample$, material_type$,                 'MATERIAL_TYPE')
call init_attribute_name1 (sample$, curvature$,                     'CURVATURE')

call init_attribute_name1 (sol_quad$, k1$,                          'K1', quasi_free$)
call init_attribute_name1 (sol_quad$, ks$,                          'KS', quasi_free$)
call init_attribute_name1 (sol_quad$, b1_gradient$,                 'B1_GRADIENT', quasi_free$)
call init_attribute_name1 (sol_quad$, bs_field$,                    'BS_FIELD', quasi_free$)
call init_attribute_name1 (sol_quad$, r0_mag$,                      'R0_MAG')
call init_attribute_name1 (sol_quad$, r0_elec$,                     'R0_ELEC') 
call init_attribute_name1 (sol_quad$, field_master$,                'FIELD_MASTER')
call init_attribute_name1 (sol_quad$, cartesian_map$,               'CARTESIAN_MAP')
call init_attribute_name1 (sol_quad$, cylindrical_map$,             'CYLINDRICAL_MAP')
call init_attribute_name1 (sol_quad$, gen_grad_map$,                'GEN_GRAD_MAP')
call init_attribute_name1 (sol_quad$, grid_field$,                  'GRID_FIELD')
call init_attribute_name1 (sol_quad$, ptc_canonical_coords$,        'PTC_CANONICAL_COORDS')

attrib_array(multipole$, k0l$:k21l$)%name = ['K0L ', &
             'K1L ', 'K2L ', 'K3L ', 'K4L ', 'K5L ', 'K6L ', 'K7L ', 'K8L ', 'K9L ', 'K10L', &
             'K11L', 'K12L', 'K13L', 'K14L', 'K15L', 'K16L', 'K17L', 'K18L', 'K19L', 'K20L', 'K21L']
attrib_array(multipole$, k0sl$:k21sl$)%name  = ['K0SL ', &
             'K1SL ', 'K2SL ', 'K3SL ', 'K4SL ', 'K5SL ', 'K6SL ', 'K7SL ', 'K8SL ', 'K9SL ', 'K10SL', &
             'K11SL', 'K12SL', 'K13SL', 'K14SL', 'K15SL', 'K16SL', 'K17SL', 'K18SL', 'K19SL', 'K20SL', 'K21SL']
attrib_array(multipole$, t0$:t21$)%name = ['T0 ', &
             'T1 ', 'T2 ', 'T3 ', 'T4 ', 'T5 ', 'T6 ', 'T7 ', 'T8 ', 'T9 ', 'T10', &
             'T11', 'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T18', 'T19', 'T20', 'T21']
attrib_array(multipole$, k0l$:t21$)%state = is_free$
call init_attribute_name1 (multipole$, l$,                          'L')
call init_attribute_name1 (multipole$, field_master$,               'FIELD_MASTER')
call init_attribute_name1 (multipole$, x_pitch$,          null_name$, does_not_exist$, .true.)
call init_attribute_name1 (multipole$, y_pitch$,          null_name$, does_not_exist$, .true.)
call init_attribute_name1 (multipole$, x_pitch_tot$,      null_name$, does_not_exist$, .true.)
call init_attribute_name1 (multipole$, y_pitch_tot$,      null_name$, does_not_exist$, .true.)

call init_attribute_name1 (ab_multipole$, l$,                       'L')
call init_attribute_name1 (ab_multipole$, field_master$,            'FIELD_MASTER')
call init_attribute_name1 (ab_multipole$, x_pitch$,       null_name$, does_not_exist$, .true.)
call init_attribute_name1 (ab_multipole$, y_pitch$,       null_name$, does_not_exist$, .true.)
call init_attribute_name1 (ab_multipole$, x_pitch_tot$,   null_name$, does_not_exist$, .true.)
call init_attribute_name1 (ab_multipole$, y_pitch_tot$,   null_name$, does_not_exist$, .true.)

call init_attribute_name1 (sad_mult$, eps_step_scale$,         'EPS_STEP_SCALE')
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
call init_attribute_name1 (sad_mult$, e1$,                     'E1')    ! SAD: ae1
call init_attribute_name1 (sad_mult$, e2$,                     'E2')    ! SAD: ae2
call init_attribute_name1 (sad_mult$, ds_step$,                'DS_STEP', dependent$)
call init_attribute_name1 (sad_mult$, num_steps$,              'NUM_STEPS', dependent$)
! sad_mult Attributes with no SAD equivalent
call init_attribute_name1 (sad_mult$, rho$,                    'RHO')   
call init_attribute_name1 (sad_mult$, ks$,                     'KS')
! Not implemented
!call init_attribute_name1 (sad_mult$, g$,                      'G')  ! No SAD equiv
!call init_attribute_name1 (sad_mult$, b_field$,                'B_FIELD')  ! No SAD equiv
!call init_attribute_name1 (sad_mult$, angle$,                  'ANGLE')
!call init_attribute_name1 (sad_mult$, rf_frequency$,           'RF_FREQUENCY') ! SAD: freq
!call init_attribute_name1 (sad_mult$, phi0$,                   'PHI0')         ! SAD: phi
!call init_attribute_name1 (sad_mult$, phi0_err$,               'PHI0_ERR')     ! SAD: dphi
!call init_attribute_name1 (sad_mult$, voltage$,                'VOLTAGE')      ! SAD: volt
!call init_attribute_name1 (sad_mult$, harmon$,                 'HARMON')       ! SAD: harm

call init_attribute_name1 (hybrid$, l$,                             'L')
call init_attribute_name1 (hybrid$, delta_e_ref$,                   'DELTA_E_REF')
call init_attribute_name1 (hybrid$, delta_ref_time$,                'DELTA_REF_TIME', override = .true.) ! Here not dependent

call init_attribute_name1 (mirror$, l$,                             'L', dependent$)
call init_attribute_name1 (mirror$, graze_angle$,                   'GRAZE_ANGLE')
call init_attribute_name1 (mirror$, critical_angle$,                'CRITICAL_ANGLE')
call init_attribute_name1 (mirror$, ref_tilt$,                      'REF_TILT')
call init_attribute_name1 (mirror$, ref_tilt_tot$,                  'REF_TILT_TOT', dependent$)
call init_attribute_name1 (mirror$, ref_wavelength$,                'REF_WAVELENGTH', dependent$)
call init_attribute_name1 (mirror$, curvature$,                     'CURVATURE')
call init_attribute_name1 (mirror$, use_reflectivity_table$,        'USE_REFLECTIVITY_TABLE')
call init_attribute_name1 (mirror$, reflectivity_table$,            'REFLECTIVITY_TABLE')

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
call init_attribute_name1 (multilayer_mirror$, curvature$,            'CURVATURE')
call init_attribute_name1 (multilayer_mirror$, material_type$,        'MATERIAL_TYPE')

call init_attribute_name1 (taylor$, ref_orbit$,                     'REF_ORBIT')
call init_attribute_name1 (taylor$, tt$,                            'TT<out><n1><n2>...')
call init_attribute_name1 (taylor$, x_ref$,                         'X_REF')
call init_attribute_name1 (taylor$, px_ref$,                        'PX_REF')
call init_attribute_name1 (taylor$, y_ref$,                         'Y_REF')
call init_attribute_name1 (taylor$, py_ref$,                        'PY_REF')
call init_attribute_name1 (taylor$, z_ref$,                         'Z_REF')
call init_attribute_name1 (taylor$, pz_ref$,                        'PZ_REF')
call init_attribute_name1 (taylor$, delta_e_ref$,                   'DELTA_E_REF')
call init_attribute_name1 (taylor$, delta_ref_time$,                'DELTA_REF_TIME', override = .true.) ! Here not dependent

call init_attribute_name1 (wiggler$, k1x$,                          'K1X', dependent$)
call init_attribute_name1 (wiggler$, k1y$,                          'K1Y', dependent$)
call init_attribute_name1 (wiggler$, kx$,                           'KX')
call init_attribute_name1 (wiggler$, l_period$,                     'L_PERIOD')
call init_attribute_name1 (wiggler$, n_period$,                     'N_PERIOD')
call init_attribute_name1 (wiggler$, b_max$,                        'B_MAX')
call init_attribute_name1 (wiggler$, g_max$,                        'G_MAX', dependent$)
call init_attribute_name1 (wiggler$, term$,                         'TERM')
call init_attribute_name1 (wiggler$, polarity$,                     'POLARITY')
call init_attribute_name1 (wiggler$, r0_mag$,                       'R0_MAG')
call init_attribute_name1 (wiggler$, r0_elec$,                      'R0_ELEC')
call init_attribute_name1 (wiggler$, field_master$,                 'FIELD_MASTER')
call init_attribute_name1 (wiggler$, cartesian_map$,                'CARTESIAN_MAP')
call init_attribute_name1 (wiggler$, cylindrical_map$,              'CYLINDRICAL_MAP')
call init_attribute_name1 (wiggler$, gen_grad_map$,                 'GEN_GRAD_MAP')
call init_attribute_name1 (wiggler$, grid_field$,                   'GRID_FIELD')
call init_attribute_name1 (wiggler$, ptc_canonical_coords$,         'PTC_CANONICAL_COORDS')
call init_attribute_name1 (wiggler$, osc_amplitude$,                'OSC_AMPLITUDE', dependent$)
call init_attribute_name1 (wiggler$, delta_ref_time_user_set$,      'DELTA_REF_TIME_USER_SET')
call init_attribute_name1 (wiggler$, delta_ref_time$,               'DELTA_REF_TIME', override = .true.)

attrib_array(undulator$, :) = attrib_array(wiggler$, :)

call init_attribute_name1 (photon_init$, l$,                         'L', dependent$)
call init_attribute_name1 (photon_init$, sig_x$,                     'SIG_X')
call init_attribute_name1 (photon_init$, sig_y$,                     'SIG_Y')
call init_attribute_name1 (photon_init$, sig_z$,                     'SIG_Z')
call init_attribute_name1 (photon_init$, sig_vx$,                    'SIG_VX')
call init_attribute_name1 (photon_init$, sig_vy$,                    'SIG_VY')
call init_attribute_name1 (photon_init$, sig_E$,                     'SIG_E')
call init_attribute_name1 (photon_init$, sig_E2$,                    'SIG_E2')
call init_attribute_name1 (photon_init$, E_center$,                  'E_CENTER')
call init_attribute_name1 (photon_init$, E2_center$,                 'E2_CENTER')
call init_attribute_name1 (photon_init$, E2_probability$,            'E2_PROBABILITY')
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
call init_attribute_name1 (photon_init$, energy_probability_curve$,  'ENERGY_PROBABILITY_CURVE')

do i = 1, n_key$
  if (attrib_array(i, l$)%name /= 'L') cycle
  if (attrib_array(i, l$)%state /= is_free$) cycle
  call init_attribute_name1 (i, accordion_edge$, 'Accordion_Edge', private$)
  call init_attribute_name1 (i, start_edge$,     'Start_Edge', private$)
  call init_attribute_name1 (i, end_edge$,       'End_Edge', private$)
  call init_attribute_name1 (i, s_position$,     'S_Position', private$)
enddo

!-----------------------------------------------------------------------
! Fill in the %units and %kind parameters.

do i = 1, n_key$
  ele%key = i
  do j = 1, ubound(attrib_array, 2)
    attrib_array(i,j)%units = attribute_units(attrib_array(i,j)%name)
    attrib_array(i,j)%kind = attribute_type(attrib_array(i,j)%name, ele)
  enddo
enddo

! We make a short list to compare against to make things go faster.

has_hkick_attributes = .false.  ! Defined in bmad_struct.f90
has_kick_attributes  = .false.  ! Defined in bmad_struct.f90
has_orientation_attributes_key = .false.  ! Defined in bmad_struct.f90

do i = 1, n_key$
  if (attrib_array(i, x_offset$)%name == 'X_OFFSET') has_orientation_attributes_key(i) = .true.
  if (attrib_array(i, kick$)%name     == 'KICK')  has_kick_attributes(i) = .true.
  if (attrib_array(i, hkick$)%name    == 'HKICK') has_hkick_attributes(i) = .true.

  call init_short_attrib_array(i)
enddo

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
  if (attrib_array(ix_key, j)%state == private$) cycle
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
! Input:
!   ix_key       -- Integer: Key index.
!   ix_attrib    -- Integer: Attribute index. 
!   name         -- Character(*): Attribute name.
!                     Should be uppercase if attrib_state = is_free$.
!                     Should contain non-uppercase characters if attrib_state = private$.
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

attrib_array(ix_key, ix_attrib)%name      = name
attrib_array(ix_key, ix_attrib)%state     = integer_option(is_free$, attrib_state)
attrib_array(ix_key, ix_attrib)%ix_attrib = ix_attrib

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
! Also see: has_attribute function.
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

if (attribute_array_init_needed) call init_attribute_name_array

has_attribs = has_orientation_attributes_key(ele%key)
if (ele%key == em_field$ .and. (ele%slave_status == super_slave$ .or. ele%slave_status == multipass_slave$)) &
      has_attribs = .false.

end function has_orientation_attributes

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_type (attrib_name, ele) result (attrib_type)
!
! Routine to return the logical type of an attribute.
!
! A "switch" attribute is an attribute whose value corresponds to some string.
! For example, the "COUPLER_AT" attirbute with value 1 corresponds to "ENTRANCE_END", etc. 
!
! A "struct" attribute is an attribute that is the name for a "structure". For example,
! CARTESIAN_MAP is the name of the structure hoding a Cartesian map.
!
! If attrib_name corresponds to a switch attribute, The routine switch_attrib_value_name can 
! be used to print the name corresponding to the attribute's value.
!
! Note: The "storage type" of an attribute is different from the "logical type" returned by
! this routine. For example, the logical type of attribute "n_slice" is integer. However, the 
! value of "n_slice" is stored as a real number in the ele_struct [in ele%value(n_slice$)]. 
!
! Input:
!   attrib_name -- Character(*): Name of the attribute. Must be upper case.
!   ele         -- ele_struct, optional: Element associated with the attribute. Needed if attrib_name
!                   can correspond to an overlay or group variable.
!
! Output:
!   attrib_type  -- Integer: Attribute type: 
!                     is_string$, is_logical$, is_integer$, is_real$, is_switch$, is_struct$ or invalid_name$
!                     Note: An overlay or group variable will be marked invalid_name$ if ele is missing.
!-

function attribute_type (attrib_name, ele) result (attrib_type)

type (ele_struct), optional :: ele
character(*) attrib_name
integer attrib_type, i, key

! Check if an overlay or group variable

key = int_garbage$
if (present(ele)) then
  key = ele%key
  if (associated(ele%control)) then
    do i = 1, size(ele%control%var)
      if (attrib_name(1:4) == 'OLD_') then
        if (attrib_name(5:) /= ele%control%var(i)%name) cycle
      else
        if (attrib_name /= ele%control%var(i)%name) cycle
      endif
      attrib_type = is_real$
      return
    enddo
  endif
endif

!

select case (attrib_name)
case ('NO_END_MARKER', 'SYMPLECTIFY', 'IS_ON', 'LIVE_BRANCH', 'HARMON_MASTER', &
      'APERTURE_LIMIT_ON', 'ABSOLUTE_TIME_TRACKING', 'AUTOSCALE_PHASE', 'GANG', &
      'AUTOSCALE_AMPLITUDE', 'PTC_EXACT_MODEL', 'PTC_EXACT_MISALIGN', 'HIGH_ENERGY_SPACE_CHARGE_ON', &
      'TAYLOR_MAP_INCLUDES_OFFSETS', 'OFFSET_MOVES_APERTURE', 'FIELD_MASTER', 'SCALE_MULTIPOLES', &
      'FLEXIBLE', 'NEW_BRANCH', 'SPIN_FRINGE_ON', 'REF_TIME_OFFSET', 'WRAP_SUPERIMPOSE', &
      'BRANCHES_ARE_COHERENT', 'E_CENTER_RELATIVE_TO_REF', 'SCALE_FIELD_TO_ONE', &
      'MULTIPOLES_ON', 'LR_SELF_WAKE_ON', 'GEO', 'SCATTER', 'SCATTER_TEST', 'DELTA_REF_TIME_USER_SET', &
      'CONSTANT_REF_ENERGY', 'CREATE_JUMBO_SLAVE', 'PTC_CANONICAL_COORDS', 'LR_WAKE%SELF_WAKE_ON', &
      'SR_WAKE%SCALE_WITH_LENGTH', 'IS_MOSAIC', 'INHERIT_FROM_FORK', 'MODE_FLIP', &
      'EXACT_MODEL', 'EXACT_MISALIGN', 'OLD_INTEGRATOR', 'RECALC', 'DETA_DS_MASTER', &
      'MODE_FLIP0', 'MODE_FLIP1', 'STATIC_LINEAR_MAP', 'USER_SETS_LENGTH', 'USE_REFLECTIVITY_TABLE')
  attrib_type = is_logical$

case ('TAYLOR_ORDER', 'N_SLICE', 'DIRECTION', 'TIME_DIR', 'VERTICAL_KICK', 'N_CELL', &
      'IX_TO_BRANCH', 'IX_TO_ELEMENT', 'NUM_STEPS', 'INTEGRATOR_ORDER', 'N_SLAVE', 'N_LORD', &
      'MAX_FRINGE_ORDER', 'UPSTREAM_ELE_DIR', 'DOWNSTREAM_ELE_DIR', 'RUNGE_KUTTA_ORDER', &
      'SAD_N_DIV_MAX', 'LONGITUDINAL_MODE', 'MOSAIC_DIFFRACTION_NUM', 'FINAL_CHARGE')
  attrib_type = is_integer$

case ('APERTURE_AT', 'APERTURE_TYPE', 'COUPLER_AT', 'FIELD_CALC', 'EXACT_MULTIPOLES', &
      'FRINGE_TYPE', 'GEOMETRY', 'FRINGE_AT', 'MAT6_CALC_METHOD', 'MATRIX', 'KICK0', &
      'ORIGIN_ELE_REF_PT', 'PARTICLE', 'PTC_FIELD_GEOMETRY', 'DEFAULT_TRACKING_SPECIES', &
      'PTC_INTEGRATION_TYPE', 'SPIN_TRACKING_METHOD', 'PTC_FRINGE_GEOMETRY', 'INTERPOLATION', &
      'TRACKING_METHOD', 'REF_ORBIT_FOLLOWS', 'REF_COORDS', 'MODE', 'CAVITY_TYPE', 'FIELD_TYPE', &
      'SPATIAL_DISTRIBUTION', 'ENERGY_DISTRIBUTION', 'VELOCITY_DISTRIBUTION', 'KEY', 'SLAVE_STATUS', &
      'LORD_STATUS', 'PHOTON_TYPE', 'ELE_ORIGIN', 'REF_ORIGIN', 'CSR_METHOD', 'SPACE_CHARGE_METHOD', &
      'MULTIPASS_REF_ENERGY', 'REF_SPECIES', 'SPECIES_OUT', 'DISTRIBUTION', 'LATTICE_TYPE', &
      'SPECIES_STRONG', 'SCATTER_METHOD', 'FIDUCIAL_PT')
  attrib_type = is_switch$

case ('TYPE', 'ALIAS', 'DESCRIP', 'SR_WAKE_FILE', 'LR_WAKE_FILE', 'LATTICE', 'PHYSICAL_SOURCE', &
     'CRYSTAL_TYPE', 'MATERIAL_TYPE', 'REFERENCE', 'TO_LINE', 'TO_ELEMENT', 'ORIGIN_ELE', 'NAME', &
     'MACHINE', 'START_EDGE', 'INPUT_ELE', 'OUTPUT_ELE')
  attrib_type = is_string$

case ('CARTESIAN_MAP', 'CYLINDRICAL_MAP', 'FIELD_OVERLAPS', 'GEN_GRAD_MAP', 'GRID_FIELD', 'REF_ORBIT', &
      'SUPERIMPOSE', 'H_MISALIGN', 'DISPLACEMENT', 'SEGMENTED', 'PIXEL', 'TERM', &
      'VAR', 'WALL', 'AMP_VS_TIME', 'FREQUENCIES', 'X_KNOT', 'SR_WAKE', 'LR_WAKE', 'CURVATURE', &
      'ENERGY_PROBABILITY_CURVE', 'REFLECTIVITY_TABLE')
  attrib_type = is_struct$

case default
  attrib_type = is_real$

  if (attrib_name(1:1) == '!') then
    attrib_type = invalid_name$
  elseif (key == int_garbage$ .or. key == overlay$ .or. key == group$) then
    if (attribute_index(0, attrib_name, can_abbreviate = .false.) == 0) attrib_type = invalid_name$
  else
    if (attribute_index(key, attrib_name, can_abbreviate = .false.) == 0) attrib_type = invalid_name$
  endif
end select

end function attribute_type 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_units (attrib_name, unrecognized_units) result (attrib_units)
!
! Routine to return the units associated with an attribute.
! Example: attrib_units('P0C') -> 'eV'
!
! Input:
!   attrib_name        -- character(*): Name of the attribute. Must be upper case.
!   unrecognized_units -- character(*), optional: String to use if the attribute name is unrecognized.
!                           Note: Non-real attributes (EG: 'TRACKING_METHOD') are not recognized.
!                           Default is ""
! Output:
!   attrib_units -- character(16): Units associated with the attribute.
!-

function attribute_units (attrib_name, unrecognized_units) result (attrib_units)

character(*) attrib_name
character(*), optional :: unrecognized_units
character(16) attrib_units

!

select case (attrib_name)

case ('DENSITY', 'DENSITY_USED')
  attrib_units = 'kg/m^3'

case ('AREA_DENSITY', 'AREA_DENSITY_USED', 'RADIATION_LENGTH', 'RADIATION_LENGTH_USED')
  attrib_units = 'kg/m^2'

case ('ALPHA_A', 'ALPHA_A0', 'ALPHA_A1', 'ALPHA_ANGLE', 'ALPHA_B', 'ALPHA_B0', 'ALPHA_B1', &
      'BBI_CONSTANT', 'B_PARAM', 'ALPHA_A_STRONG', 'ALPHA_B_STRONG', 'N_PARTICLE', 'DTHICKNESS_DX', &
      'CHARGE', 'CMAT_11', 'CMAT_12', 'CMAT_21', 'CMAT_22', 'COUPLER_STRENGTH', 'DE_ETA_MEAS', 'F_FACTOR', &
      'ELECTRIC_DIPOLE_MOMENT', 'ETAP_X', 'ETAP_X0', 'ETAP_X1', 'ETAP_Y', 'ETAP_Y0', 'ETAP_Y1', &
      'ETAP_X_OUT', 'ETAP_Y_OUT', 'EMIT_FRACTION', 'Y_KNOT', 'SLAVE', 'DALPHA_DPZ_A', 'DALPHA_DPZ_B', &
      'FIELD_AUTOSCALE', 'FIELD_SCALE_FACTOR', 'FIELD_X', 'FIELD_Y', 'FINT', 'FINTX', 'GAP', 'HARMON', 'HKICK', &
      'KICK', 'MAX_NUM_RUNGE_KUTTA_STEP', 'NOISE', 'N_PART', 'N_PERIOD', 'N_SAMPLE', 'N_SLICE_SPLINE', &
      'POLARITY', 'PX', 'PX0', 'PX1', 'PX_REF', 'PY', 'PY0', 'PY1', 'PY_REF', 'PZ', 'PZ0', 'PZ1', 'PZ_REF', &
      'RAN_SEED', 'REF_CAP_GAMMA', 'REL_TOL_ADAPTIVE_TRACKING', 'REL_TOL_TRACKING', 'SIG_PZ', 'DETAP_DPZ_X', 'DETAP_DPZ_Y', &
      'SPIN_X', 'SPIN_Y', 'SPIN_Z', 'TRANSVERSE_SIGMA_CUT', 'VKICK', 'LONGITUDINAL_MODE', 'MOSAIC_DIFFRACTION_NUM', &
      'AUTOSCALE_AMP_REL_TOL', 'PX_KICK', 'PY_KICK', 'PZ_KICK', 'SPIN_DN_DPZ_X', 'SPIN_DN_DPZ_Y', 'SPIN_DN_DPZ_Z', &
      'VAL1', 'VAL2', 'VAL3', 'VAL4', 'VAL5', 'VAL6', 'VAL7', 'VAL8', 'VAL9', 'VAL10', 'VAL11', 'VAL12', &
      'C11_MAT0', 'C11_MAT1', 'C22_MAT0', 'C22_MAT1', 'E2_PROBABILITY', 'CRAB_X1', 'PZ_APERTURE_CENTER', &
      'PX_APERTURE_WIDTH2', 'PX_APERTURE_CENTER', 'PY_APERTURE_WIDTH2', 'PY_APERTURE_CENTER', 'PZ_APERTURE_WIDTH2')
  attrib_units = ''

case ('SIG_VX', 'SIG_VY')
  attrib_units = 'm/s'

case ('ABS_TOL_ADAPTIVE_TRACKING', 'ABS_TOL_TRACKING', 'ACCORDION_EDGE', 'APERTURE', &
      'BETA_A', 'BETA_A0', 'BETA_A1', 'BETA_B', 'BETA_B0', 'BETA_B1', 'BETA_A_STRONG', 'BETA_B_STRONG', &
      'D1_THICKNESS', 'D2_THICKNESS', 'DEFAULT_DS_STEP', 'OSC_AMPLITUDE', 'R_SOLENOID', &
      'DS_SLICE', 'DS_STEP', 'DX_ORIGIN', 'DY_ORIGIN', 'DZ_ORIGIN', 'D_SPACING', 'END_EDGE', 'EPS_STEP_SCALE', &
      'ETA_X_OUT', 'ETA_Y_OUT', 'CSR_DS_STEP', 'X_KICK', 'Y_KICK', 'Z_KICK', 'DBETA_DPZ_A', 'DBETA_DPZ_B', &
      'ETA_X', 'ETA_X0', 'ETA_X1', 'ETA_Y', 'ETA_Y0', 'ETA_Y1', 'ETA_Z', 'FATAL_DS_ADAPTIVE_TRACKING', &
      'FB1', 'FB2', 'FQ1', 'FQ2', 'HGAP', 'HGAPX', 'H_DISPLACE', 'INIT_DS_ADAPTIVE_TRACKING', 'L', 'DETA_DPZ_X', 'DETA_DPZ_Y', &
      'LORD_PAD1', 'LORD_PAD2', 'L_CHORD', 'L_ACTIVE', 'L_SOFT_EDGE', 'L_PERIOD', 'L_SAGITTA', 'MAX_APERTURE_LIMIT', &
      'MIN_DS_ADAPTIVE_TRACKING', 'OFFSET', 'PENDELLOSUNG_PERIOD_PI', 'PENDELLOSUNG_PERIOD_SIGMA', 'R0_ELEC', 'R0_MAG', &
      'REF_WAVELENGTH', 'RHO', 'S', 'SIGNIFICANT_LENGTH', 'SIG_X', 'SIG_Y', 'SIG_Z', 'S_POSITION', 'THICKNESS', &
      'X', 'X0', 'X1', 'Y', 'Y0', 'Y1', 'X1_LIMIT', 'X2_LIMIT', 'Y1_LIMIT', 'Y2_LIMIT', 'X_LIMIT', 'Y_LIMIT', &
      'X_OFFSET', 'Y_OFFSET', 'X_OFFSET_CALIB', 'Y_OFFSET_CALIB', 'X_OFFSET_MULT', 'Y_OFFSET_MULT', &
      'X_OFFSET_TOT', 'Y_OFFSET_TOT', 'X_POSITION', 'Y_POSITION', 'X_QUAD', 'Y_QUAD', &
      'X_REF', 'Y_REF', 'Z', 'Z0', 'Z1', 'Z_OFFSET', 'Z_OFFSET_TOT', 'Z_POSITION', 'Z_REF', 'MOSAIC_THICKNESS', &
      'C12_MAT0', 'C12_MAT1', 'X_GAIN_CALIB', 'Y_GAIN_CALIB', 'X_GAIN_ERR', 'Y_GAIN_ERR', 'RADIUS', &
      'Z_APERTURE_WIDTH2', 'Z_APERTURE_CENTER', 'RF_WAVELENGTH', &
      'X1_EDGE', 'X2_EDGE', 'Y1_EDGE', 'Y2_EDGE', 'L_RECTANGLE', 'S_TWISS_REF', &
      'X_DISPERSION_ERR', 'Y_DISPERSION_ERR', 'X_DISPERSION_CALIB', 'Y_DISPERSION_CALIB')
  attrib_units = 'm'


case ('V1_UNITCELL', 'V2_UNITCELL', 'V_DISPLACE', 'V_UNITCELL')
  attrib_units = 'm^3'

case ('ANGLE', 'BEND_TILT', 'BRAGG_ANGLE', 'BRAGG_ANGLE_IN', 'BRAGG_ANGLE_OUT', &
      'COUPLER_ANGLE', 'CRITICAL_ANGLE', 'CRUNCH', 'CRUNCH_CALIB', 'DARWIN_WIDTH_PI', 'DARWIN_WIDTH_SIGMA', &
      'DPHI_A', 'DPHI_B', 'DPHI_ORIGIN', 'DPSI_ORIGIN', 'DTHETA_ORIGIN', 'E1', 'E2', 'GRAZE_ANGLE', & 
      'PHASE_X', 'PHASE_Y', 'PHI_A', 'PHI_B', 'PHI_POSITION', 'PSI_ANGLE', 'PSI_POSITION', 'QUAD_TILT', &
      'REF_TILT', 'REF_TILT_TOT', 'ROLL', 'ROLL_TOT', 'THETA_POSITION', &
      'X_PITCH', 'Y_PITCH', 'X_PITCH_TOT', 'Y_PITCH_TOT', &
      'TILT', 'TILT_CALIB', 'TILT_CORR', 'TILT_TOT', 'GRAZE_ANGLE_IN', 'GRAZE_ANGLE_OUT', &
      'T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', &
      'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T18', 'T19', 'T20', 'T21', &
      'MOSAIC_ANGLE_RMS_IN_PLANE', 'MOSAIC_ANGLE_RMS_OUT_PLANE', 'ANGLE_OUT_MAX', 'CRAB_TILT')
  attrib_units = 'rad'

case ('COUPLER_PHASE', 'PHI0', 'PHI0_AUTOSCALE', 'PHI0_ERR', 'PHI0_MULTIPASS', 'AUTOSCALE_PHASE_TOL')
  attrib_units = 'rad/2pi'

case ('CRITICAL_ANGLE_FACTOR')
  attrib_units = 'rad*eV'

case ('C21_MAT0', 'C21_MAT1', 'CURVATURE_X0_Y2', 'CURVATURE_X1_Y1', 'CURVATURE_X2_Y0', 'G', 'DG', 'G_MAX', &
      'G_TOT', 'H1', 'H2', 'CRAB_X2', 'SPHERICAL_CURVATURE', 'ELLIPTICAL_CURVATURE_X', 'ELLIPTICAL_CURVATURE_Y', &
      'ELLIPTICAL_CURVATURE_Z', 'FOCAL_STRENGTH')
  attrib_units = '1/m'

case ('CURVATURE_X0_Y3', 'CURVATURE_X1_Y2', 'CURVATURE_X2_Y1', 'CURVATURE_X3_Y0', 'DKS_DS', 'CRAB_X3')
  attrib_units = '1/m^2'

case ('CURVATURE_X0_Y4', 'CURVATURE_X1_Y3', 'CURVATURE_X2_Y2', 'CURVATURE_X3_Y1', 'CURVATURE_X4_Y0', 'CRAB_X4')
  attrib_units = '1/m^3'

case ('CURVATURE_X0_Y5', 'CURVATURE_X1_Y4', 'CURVATURE_X2_Y3', 'CURVATURE_X3_Y2', 'CURVATURE_X4_Y1', &
      'CURVATURE_X5_Y0', 'CRAB_X5')
  attrib_units = '1/m^4'

case ('CURVATURE_X0_Y6', 'CURVATURE_X1_Y5', 'CURVATURE_X2_Y4', 'CURVATURE_X3_Y3', 'CURVATURE_X4_Y2', &
      'CURVATURE_X5_Y1', 'CURVATURE_X6_Y0')
  attrib_units = '1/m^5'

case ('DBRAGG_ANGLE_DE')
  attrib_units = 'rad/eV'

case ('DELTA_E', 'ENERGY', 'E_CENTER', 'E2_CENTER', 'E_LOSS', 'E_PHOTON', 'E_TOT', 'E_TOT_OFFSET', 'E_TOT_START', &
      'P0C', 'P0C_START', 'PC', 'P0C_SET', 'E_TOT_SET', 'AUTOSCALE_AMP_ABS_TOL', 'DELTA_E_REF', 'SIG_E', 'SIG_E2', &
      'PC_OUT_MIN', 'PC_OUT_MAX', 'E_TOT_STRONG', 'PC_STRONG')
  attrib_units = 'eV'

case ('DELTA_REF_TIME', 'REF_TIME', 'REF_TIME_START', 'T', 'T_OFFSET', 'DELTA_TIME', 'DT_MAX', 'CROSSING_TIME')
  attrib_units = 'sec'

case ('EMITTANCE_A', 'EMITTANCE_B', 'EMITTANCE_Z')
  attrib_units = 'm*rad'

case ('E_FIELD', 'E_FIELD_X', 'E_FIELD_Y')
  attrib_units = 'V/m'

case ('VOLTAGE', 'VOLTAGE_ERR', 'VOLTAGE_TOT')
  attrib_units = 'Volt'

case ('GRADIENT', 'GRADIENT_ERR', 'GRADIENT_TOT')
  attrib_units = 'eV/m'

case ('LR_FREQ_SPREAD', 'RF_FREQUENCY', 'REPETITION_FREQUENCY')
  attrib_units = 'Hz'

case ('BS_FIELD', 'B_FIELD', 'B_FIELD_TOT', 'DB_FIELD', 'B_MAX')
  attrib_units = 'T'

case ('B1_GRADIENT');                                   attrib_units = 'T/m'
case ('B2_GRADIENT');                                   attrib_units = 'T/m^2'
case ('B3_GRADIENT');                                   attrib_units = 'T/m^3'
case ('BL_HKICK', 'BL_KICK', 'BL_VKICK');               attrib_units = 'T*m'
case ('A0', 'B0', 'K0L', 'K0SL');                       attrib_units = ''
case ('A1', 'B1', 'K1L', 'K1SL', 'KS', 'KX');           attrib_units = '1/m'
case ('A2', 'B2', 'K2L', 'K2SL', 'K1', 'K1X', 'K1Y');   attrib_units = '1/m^2'
case ('A3', 'B3', 'K3L', 'K3SL', 'K2');                 attrib_units = '1/m^3'
case ('A4', 'B4', 'K4L', 'K4SL', 'K3');                 attrib_units = '1/m^4'
case ('A5', 'B5', 'K5L', 'K5SL');                       attrib_units = '1/m^5'
case ('A6', 'B6', 'K6L', 'K6SL');                       attrib_units = '1/m^6'
case ('A7', 'B7', 'K7L', 'K7SL');                       attrib_units = '1/m^7'
case ('A8', 'B8', 'K8L', 'K8SL');                       attrib_units = '1/m^8'
case ('A9', 'B9', 'K9L', 'K9SL');                       attrib_units = '1/m^9'
case ('A10', 'B10', 'K10L', 'K10SL');                   attrib_units = '1/m^10'
case ('A11', 'B11', 'K11L', 'K11SL');                   attrib_units = '1/m^11'
case ('A12', 'B12', 'K12L', 'K12SL');                   attrib_units = '1/m^12'
case ('A13', 'B13', 'K13L', 'K13SL');                   attrib_units = '1/m^13'
case ('A14', 'B14', 'K14L', 'K14SL');                   attrib_units = '1/m^14'
case ('A15', 'B15', 'K15L', 'K15SL');                   attrib_units = '1/m^15'
case ('A16', 'B16', 'K16L', 'K16SL');                   attrib_units = '1/m^16'
case ('A17', 'B17', 'K17L', 'K17SL');                   attrib_units = '1/m^17'
case ('A18', 'B18', 'K18L', 'K18SL');                   attrib_units = '1/m^18'
case ('A19', 'B19', 'K19L', 'K19SL');                   attrib_units = '1/m^19'
case ('A20', 'B20', 'K20L', 'K20SL');                   attrib_units = '1/m^20'
case ('A21', 'B21', 'K21L', 'K21SL');                   attrib_units = '1/m^21'
case ('A0_ELEC', 'B0_ELEC');                            attrib_units = 'V/m'
case ('A1_ELEC', 'B1_ELEC');                            attrib_units = 'V/m^2'
case ('A2_ELEC', 'B2_ELEC');                            attrib_units = 'V/m^3'
case ('A3_ELEC', 'B3_ELEC');                            attrib_units = 'V/m^4'
case ('A4_ELEC', 'B4_ELEC');                            attrib_units = 'V/m^5'
case ('A5_ELEC', 'B5_ELEC');                            attrib_units = 'V/m^6'
case ('A6_ELEC', 'B6_ELEC');                            attrib_units = 'V/m^7'
case ('A7_ELEC', 'B7_ELEC');                            attrib_units = 'V/m^8'
case ('A8_ELEC', 'B8_ELEC');                            attrib_units = 'V/m^9'
case ('A9_ELEC', 'B9_ELEC');                            attrib_units = 'V/m^10'
case ('A10_ELEC', 'B10_ELEC');                          attrib_units = 'V/m^11'
case ('A11_ELEC', 'B11_ELEC');                          attrib_units = 'V/m^12'
case ('A12_ELEC', 'B12_ELEC');                          attrib_units = 'V/m^13'
case ('A13_ELEC', 'B13_ELEC');                          attrib_units = 'V/m^14'
case ('A14_ELEC', 'B14_ELEC');                          attrib_units = 'V/m^15'
case ('A15_ELEC', 'B15_ELEC');                          attrib_units = 'V/m^16'
case ('A16_ELEC', 'B16_ELEC');                          attrib_units = 'V/m^17'
case ('A17_ELEC', 'B17_ELEC');                          attrib_units = 'V/m^18'
case ('A18_ELEC', 'B18_ELEC');                          attrib_units = 'V/m^19'
case ('A19_ELEC', 'B19_ELEC');                          attrib_units = 'V/m^20'
case ('A20_ELEC', 'B20_ELEC');                          attrib_units = 'V/m^21'
case ('A21_ELEC', 'B21_ELEC');                          attrib_units = 'V/m^22'

case default
  if (present(unrecognized_units)) then
    attrib_units = unrecognized_units
  else
    attrib_units = ''
  endif

end select

end function attribute_units

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
case ('NAME')
  attrib_value = ele%name
case ('TYPE')
  attrib_value = ele%type
case ('ALIAS')
  attrib_value = ele%alias
case ('DESCRIP')
  if (associated(ele%descrip)) attrib_value = ele%descrip
case ('SR_WAKE_FILE')
  if (associated(ele%wake)) attrib_value = ele%wake%sr%file
case ('LR_WAKE_FILE')
  if (associated(ele%wake)) attrib_value = ele%wake%lr%file
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
end select


end subroutine string_attrib

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function switch_attrib_value_name (attrib_name, attrib_value, ele, 
!                                     is_default, name_list) result (attrib_val_name)
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
! Input:
!   attrib_name    -- character(*): Name of the attribute. Must be upper case.
!   attrib_value   -- real(rp): Value of the attribute.
!   ele            -- ele_struct: Lattice element that the attribute is contained in.
!                       Generally only needed to determine the default value.
!
! Output:
!   attrib_val_name -- character(40): Name corresponding to the value. Set to null_name$ if there is a problem.
!   is_default      -- logical, optional: If True then the value of the attiribute
!                        corresponds to the default value. If this argument is
!                        present, the ele argument must also be present.
!   name_list(:)    -- character(40), allocatable, optional: List of names the switch can take.
!                         Deallocated if there is an error.
!-

function switch_attrib_value_name (attrib_name, attrib_value, ele, is_default, name_list) result (attrib_val_name)

type (ele_struct) :: ele, ele2
character(*) attrib_name
real(rp) attrib_value
integer ix_attrib_val
logical, optional :: is_default
character(*), optional, allocatable :: name_list(:)
character(40) attrib_val_name
character(*), parameter :: r_name = 'switch_attrib_value_name'

! 

ix_attrib_val = nint(attrib_value)

if (present(name_list)) then
  if (allocated(name_list)) deallocate(name_list)
endif

!

select case (attrib_name)

case ('APERTURE_AT')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, aperture_at_name, lbound(aperture_type_name, 1), name_list)
  if (present(is_default)) then
    is_default = (ix_attrib_val == exit_end$)
  endif

case ('APERTURE_TYPE')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, aperture_type_name, lbound(aperture_type_name, 1), name_list)
  if (present(is_default)) then
    if (ele%key == ecollimator$) then
      is_default = (ix_attrib_val == elliptical$)
    else
      is_default = (ix_attrib_val == rectangular$)
    endif
  endif

case ('CAVITY_TYPE')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, cavity_type_name, lbound(cavity_type_name, 1), name_list)
  if (present(is_default)) then
    is_default = (ix_attrib_val == default_value(ele, cavity_type$))
  endif

case ('CSR_METHOD')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, csr_method_name, lbound(csr_method_name, 1), name_list)
  if (present(is_default)) then
    is_default = (ix_attrib_val == off$)
  endif

case ('COUPLER_AT')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, end_at_name, lbound(end_at_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == downstream_end$)

case ('DEFAULT_TRACKING_SPECIES')
  attrib_val_name = species_name(ix_attrib_val)
  if (present(is_default)) is_default = (ix_attrib_val == ref_particle$)

case ('ELE_ORIGIN', 'REF_ORIGIN')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, anchor_pt_name, lbound(anchor_pt_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == anchor_center$)

case ('ENERGY_DISTRIBUTION')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, distribution_name, lbound(distribution_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == gaussian$)

case ('EXACT_MULTIPOLES')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, exact_multipoles_name, lbound(exact_multipoles_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == off$)

case ('FIDUCIAL_PT')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, fiducial_pt_name, lbound(fiducial_pt_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == none_pt$)

case ('FIELD_CALC')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, field_calc_name, lbound(field_calc_name, 1), name_list)
  if (present(is_default)) then
    select case (ele%key)
    case (group$, overlay$, girder$, ramper$); is_default = (ix_attrib_val == no_field$)
    case default;                              is_default = (ix_attrib_val == bmad_standard$)
    end select
  endif

case ('FIELD_TYPE')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, em_field_type_name, lbound(em_field_type_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == bmad_standard$)

case ('MULTIPASS_REF_ENERGY')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, multipass_ref_energy_name, lbound(multipass_ref_energy_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == first_pass$)

case ('NONGRID^FIELD_TYPE')      ! This is for the Tao "python" command
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, em_field_type_name(1:2), lbound(em_field_type_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == bmad_standard$)

case ('FRINGE_AT')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, end_at_name, lbound(end_at_name, 1), name_list)
  if (present(is_default)) then
    is_default = (ix_attrib_val == both_ends$)
  endif

case ('FRINGE_TYPE')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, fringe_type_name, lbound(fringe_type_name, 1), name_list)
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

case ('GEOMETRY')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, geometry_name, lbound(geometry_name, 1), name_list)

case ('GRID_FIELD^GEOMETRY')      ! This is for the Tao "python" command
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, grid_field_geometry_name, lbound(grid_field_geometry_name, 1), name_list)

case ('INTERPOLATION')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, interpolation_name, lbound(interpolation_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == cubic$)

case ('KEY')
    call get_this_attrib_name (attrib_val_name, ix_attrib_val, key_name, lbound(key_name, 1), name_list)
  if (present(is_default)) is_default = .false.

case ('KICK0')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, kick0_name, lbound(kick0_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == standard$)

case ('LORD_STATUS')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, control_name, lbound(control_name, 1), name_list)
  if (present(is_default)) is_default = .false.
  
case ('MAT6_CALC_METHOD')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, mat6_calc_method_name, lbound(mat6_calc_method_name, 1), name_list)
  if (present(is_default)) then
    call default_ele(ele, ele2)
    is_default = (ix_attrib_val == ele2%mat6_calc_method)
  endif

case ('MATRIX')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, matrix_name, lbound(matrix_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == standard$)

case ('MODE')
  if (ele%key == diffraction_plate$ .or. ele%key == sample$ .or. ele%key == mask$) then
    call get_this_attrib_name (attrib_val_name, ix_attrib_val, mode_name, lbound(mode_name, 1), name_list)
    if (present(is_default)) then
      is_default = (ix_attrib_val == default_value(ele, mode$))
    endif
  else
    call get_this_attrib_name (attrib_val_name, ix_attrib_val, geometry_name, lbound(geometry_name, 1), name_list)
    if (present(is_default)) then
      is_default = (ix_attrib_val == open$)
    endif
  endif

case ('ORIGIN_ELE_REF_PT')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, ref_pt_name, lbound(ref_pt_name, 1), name_list)
    if (present(is_default)) is_default = (ix_attrib_val == center_pt$)

case ('PTC_FRINGE_GEOMETRY')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, ptc_fringe_geometry_name, lbound(ptc_fringe_geometry_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == x_invariant$)

case ('PHOTON_TYPE')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, photon_type_name, lbound(photon_type_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == incoherent$)

case ('PARTICLE', 'REF_SPECIES', 'SPECIES_STRONG')
  attrib_val_name = species_name(ix_attrib_val)
  if (present(is_default)) then
    if (ele%key == photon_fork$) then
      is_default = (ix_attrib_val == photon$)
    else
      is_default = .false. ! Cannot tell so assume the worst.
    endif
  endif

case ('PHASE_UNITS')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, angle_units_name, lbound(angle_units_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == radians$)

case ('PTC_FIELD_GEOMETRY')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, ptc_field_geometry_name, lbound(ptc_field_geometry_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == sector$)

case ('PTC_INTEGRATION_TYPE')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, ptc_integration_type_name, lbound(ptc_integration_type_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == matrix_kick$)

case ('REF_COORDS')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, ref_coords_name(1:4), 1, name_list)
  if (present(is_default)) is_default = (ix_attrib_val == exit_end$)

case ('REF_ORBIT_FOLLOWS')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, ref_orbit_follows_name, lbound(ref_orbit_follows_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == bragg_diffracted$)

case ('SCATTER_METHOD')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, scatter_method_name, lbound(scatter_method_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == highland$)

case ('SECTION^TYPE')    ! This is for the Tao "python" command
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, wall3d_section_type_name, lbound(wall3d_section_type_name, 1), name_list)
  if (present(is_default)) is_default = .false.

case ('SLAVE_STATUS')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, control_name, lbound(control_name, 1), name_list)
  if (present(is_default)) is_default = .false.
  
case ('SPACE_CHARGE_METHOD')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, space_charge_method_name, lbound(space_charge_method_name, 1), name_list)
  if (present(is_default)) then
    is_default = (ix_attrib_val == off$)
  endif

case ('SPATIAL_DISTRIBUTION')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, distribution_name, lbound(distribution_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == gaussian$)

case ('SPIN_TRACKING_METHOD')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, spin_tracking_method_name, lbound(spin_tracking_method_name, 1), name_list)
  if (present(is_default)) then
    call default_ele (ele, ele2)
    is_default = (ix_attrib_val == ele2%spin_tracking_method)
  endif

case ('TRACKING_METHOD')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, tracking_method_name, lbound(tracking_method_name, 1), name_list)
  if (present(is_default)) then
    call default_ele (ele, ele2)
    is_default = (ix_attrib_val == ele2%tracking_method)
  endif

case ('VELOCITY_DISTRIBUTION')
  call get_this_attrib_name (attrib_val_name, ix_attrib_val, distribution_name, lbound(distribution_name, 1), name_list)
  if (present(is_default)) is_default = (ix_attrib_val == gaussian$)

case default
  call out_io (s_fatal$, r_name, 'BAD ATTRIBUTE NAME: ' // attrib_name)
  attrib_val_name = null_name$
end select

!---------------------------------------
contains

subroutine get_this_attrib_name (val_name, ix_attrib_val, name_array, min_arr, name_list, exceptions)

integer ix_attrib_val, min_arr, i, j, n
integer, optional :: exceptions(:)
character(*) val_name
character(*) name_array(min_arr:)
character(*), optional, allocatable :: name_list(:)

!

if (ix_attrib_val < lbound(name_array, 1) .or. ix_attrib_val > ubound(name_array, 1)) then
  val_name = null_name$
else
  val_name = name_array(ix_attrib_val)
  if (present(exceptions)) then
    if (any (ix_attrib_val == exceptions)) val_name = null_name$
  endif
endif

n = 0
if (present(exceptions)) n = size(exceptions)

if (present(name_list)) then
  allocate (name_list(lbound(name_array,1):ubound(name_array,1)-n))
  j = lbound(name_array, 1) - 1
  do i = lbound(name_array, 1), ubound(name_array, 1)
    if (present(exceptions)) then
      if (any (i == exceptions)) cycle
    endif
    j = j + 1
    name_list(j) = name_array(i)
  enddo
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

end function n_attrib_string_max_len

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Fuction has_attribute (ele, attrib) result(has_it)
!
! Routine to determine if a given type of lattice element has a particular attribute.
! Also see: has_orientation_attributes function.
!
! Input:
!   ele     -- ele_struct: Lattice element of a certain type
!   attrib  -- character(*): Name of the attribute. Must be upper case.
!
! Output:
!   has_it  -- logical: True if element has an attribute of that name.
!-

function has_attribute (ele, attrib) result (has_it)

type (ele_struct) ele
character(*) attrib
logical has_it

!

has_it = (attribute_index(ele, attrib) > 0)

end function has_attribute

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function custom_attribute_ubound_index(ele_class) result (ix_ubound)
! 
! Routine to return, for a given element class, the upper bound index for the ele%custom(:) 
! array which is needed to accomodate the registered custom attributes for that class.
!
! Input:
!   ele_class   -- integer: Element class (key). EG: quadrupole$, etc.
!
! Output:
!   ix_ubound   -- integer: Maximum index needed.
!-

function custom_attribute_ubound_index(ele_class) result (ix_ubound)

integer ele_class, ix_ubound

!

do ix_ubound = custom_attribute_num$, 1, -1
  if (attribute_name(ele_class, ix_ubound+custom_attribute0$) /= null_name$) exit
enddo

end function custom_attribute_ubound_index

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine set_custom_attribute_name (custom_name, err_flag, custom_index)
!
! Routine to add custom element attributes to the element attribute name table.
!
! Input:
!   custom_name  -- Character(*): Name of the custom attribute. If prefixed by "<class>::" then
!                    the custom name will be set only for that element class. Example:
!                    "quadrupole::error" will set the alias custom namefor quadrupoles.
!   custom_index -- Integer, optional: Index used in assigning where in the ele_struct the custom
!                     attribute is put. If not present or 0 then the next unused slot is used.
!
! Output:
!   err_flag    -- Logical: Set True if an error. False otherwise.
!-

subroutine set_custom_attribute_name (custom_name, err_flag, custom_index)

implicit none

integer, optional :: custom_index
integer i, ix, ix_custom, key, ix_key
character(*) custom_name
character(40) custom_str, a_str
character(*), parameter :: r_name = 'set_custom_attribute_name'
logical :: err_flag

! If custom_str is of the form "ele_class::custom_name" then need to split string

if (attribute_array_init_needed) call init_attribute_name_array
err_flag = .true.

custom_str = upcase(custom_name)
a_str = custom_str
key = 0
ix = index(custom_str, '::')
if (ix /= 0) then
  a_str = custom_str(ix+2:)
  select case (custom_str(1:ix-1))
  case ('PARAMETER')
    key = def_parameter$
  case default
    key = key_name_to_key_index(custom_str(1:ix-1), .true.)
    if (key < 1) then
      call out_io (s_error$, r_name, 'ELEMENT CLASS NOT RECOGNIZED: ' // custom_name)
      return
    endif
  end select
endif

! Check custom name for invalid characters

if (str_first_not_in_set(trim(a_str), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_') > 0) then
  call out_io (s_error$, r_name, 'CUSTOM NAME HAS INVALID CHARACTERS: ' // a_str)
  return
endif

! Note: The null_ele$ slot is where the "common" attribute name for a given index is stored.

ix_custom = integer_option(0, custom_index)
if (ix_custom < 0 .or. ix_custom > custom_attribute_num$) then
  call out_io (s_error$, r_name, 'INDEX NUMBER FOR CUSTOM_ATTRIBUTE OUT OF RANGE: \i0\ ', ix_custom)
  return
endif

if (ix_custom == 0) then
  ix_key = key
  if (ix_key == 0) ix_key = null_ele$
  do i = 1, custom_attribute_num$
    if (attribute_name(ix_key, i+custom_attribute0$) /= null_name$) cycle
    call set_it (i)
    return
  enddo
  ! Should not be here
  call out_io (s_error$, r_name, 'NO BLANK SLOT FOUND FOR PUTTING THIS CUSTOM ATTRIBUTE: ' // custom_name)
  return

else
  call set_it (ix_custom)
endif

!---------------------------------------------------------------
contains

subroutine set_it (ix_attrib)

integer ix_attrib, i
character(40) old_attrib

! Set common custom attribute
! Note: A common custom attribute will not override a non-common custom attribute that has already been set.

if (key == 0) then
  old_attrib = attribute_name(null_ele$, ix_attrib+custom_attribute0$)
  if (old_attrib /= null_name$ .and. old_attrib /= a_str) then
    call out_io (s_warn$, r_name, &
      'A CUSTOM_ATTRIBUTE IS BEING REDEFINED: ' // custom_name, &
      'FROM: ' // attribute_name(null_ele$, ix_attrib), &
      'TO:   ' // a_str)
    return
  endif

  do i = 1, n_key$
    select case (i)
    case (def_parameter$, def_mad_beam$, def_bmad_com$, def_particle_start$, def_line$, def_space_charge_com$); cycle
    end select
    old_attrib = attribute_name(i, ix_attrib+custom_attribute0$)
    if (old_attrib /= null_name$) cycle
    call init_attribute_name1 (i, ix_attrib+custom_attribute0$, a_str, override = .true.)
  enddo

! Set a element type specific (that is, non-common) custom attribute.
! Note: If the custom attribute slot is already being used this is an error except if the new
! name matches the existing name.

else

  old_attrib = attribute_name(key, ix_attrib+custom_attribute0$)
  if (old_attrib /= null_name$ .and. old_attrib /= attribute_name(null_ele$, ix_attrib) .and. old_attrib /= a_str) then
    call out_io (s_warn$, r_name, &
      'A CUSTOM_ATTRIBUTE IS BEING REDEFINED: ' // custom_name, &
      'FOR ELEMENT CLASS: ' // key_name(key), &
      'FROM: ' // old_attrib, &
      'TO:   ' // a_str)
    return
  endif

  call init_attribute_name1 (key, ix_attrib+custom_attribute0$, a_str, override = .true.)
endif

err_flag = .false.

end subroutine set_it

end subroutine set_custom_attribute_name

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine custom_ele_attrib_name_list (index_list, name_list)
!
! Routine to create an array (index_list(i), name_list(i)) of custom element attribute names and indexes.
! Each name in the name_list is of the form:
!   "{<class>::}<attribute_name>"
! where:
!   <class>:: is an optional class prefix.
!   <attribute_name> is the name of the attribute.
!
! Output:
!   index_list(:) -- integer, allocatable: Index of the custom attribute.
!   name_list(:)  -- character(*), allocatable: List of custom attributes.
!-

subroutine custom_ele_attrib_name_list (index_list, name_list)

implicit none

integer, allocatable :: index_list(:)
integer n, j, ic, icc, is, n_name(n_key$)

character(*), allocatable :: name_list(:)
character(40) common_name, k_name

!

call re_allocate(name_list, 10)
call re_allocate(index_list, 10)

n = 0  ! Number of elements in the list

!

do ic = 1, custom_attribute_num$
  icc = ic + custom_attribute0$

  ! Common attribute is stored in null_ele's slot

  common_name = ''
  if (attribute_name(null_ele$, icc) /= null_name$) then
    common_name = attribute_name(null_ele$, icc)
    n = n + 1
    if (n > size(name_list)) call re_allocate(name_list, n+10)
    if (n > size(index_list)) call re_allocate(index_list, n+10)
    index_list(n) = ic
    name_list(n) = common_name
  endif

  ! For not common attributes

  do is = 1, n_key$
    select case (is)
    case (def_mad_beam$, def_bmad_com$, def_space_charge_com$, def_particle_start$, def_line$); cycle
    end select
    if (attribute_name(is, icc) == null_name$) cycle
    if (attribute_name(is, icc) == common_name) cycle
    n = n + 1
    if (n > size(name_list)) call re_allocate(name_list, n+10)
    if (n > size(index_list)) call re_allocate(index_list, n+10)
    k_name = key_name(is)
    if (is == def_parameter$) k_name = 'Parameter'   ! Not "def_parameter"
    index_list(n) = ic
    name_list(n) = trim(k_name) // '::' // attribute_name(is, icc)
  enddo
enddo

!

call re_allocate(index_list, n)
call re_allocate(name_list, n)

end subroutine custom_ele_attrib_name_list

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free1 (ix_ele, attrib_name, lat, err_print_flag, except_overlay, dependent_attribs_free, why_not_free) result (free)
!
! This function overloaded by attribute_free. See attribute_free for more details.
!-

function attribute_free1 (ix_ele, attrib_name, lat, err_print_flag, except_overlay, dependent_attribs_free, why_not_free) result (free)

implicit none

type (lat_struct) :: lat

integer ix_ele
integer, optional :: why_not_free

character(*) attrib_name

logical free
logical, optional :: err_print_flag, except_overlay, dependent_attribs_free

!

free = check_this_attribute_free (lat%ele(ix_ele), attrib_name, lat, err_print_flag, except_overlay, dependent_attribs_free, why_not_free)

end function attribute_free1

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free2 (ele, attrib_name, err_print_flag, except_overlay, dependent_attribs_free, why_not_free) result (free)
!
! This function overloaded by attribute_free. See attribute_free for more details.
!-

function attribute_free2 (ele, attrib_name, err_print_flag, except_overlay, dependent_attribs_free, why_not_free) result (free)

implicit none

type (lat_struct), target :: lat
type (ele_struct) ele

character(*) attrib_name

integer, optional :: why_not_free

logical free
logical, optional :: err_print_flag, except_overlay, dependent_attribs_free

character(16) :: r_name = 'attribute_free'

! Elements not assocaited with a lattice are considered free.

if (.not. associated(ele%branch)) then
  free = .true.
  return
endif

! init & check

free = check_this_attribute_free (ele, attrib_name, ele%branch%lat, err_print_flag, except_overlay, dependent_attribs_free, why_not_free)

end function attribute_free2

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free3 (ix_ele, ix_branch, attrib_name, lat, 
!                                   err_print_flag, except_overlay, dependent_attribs_free, why_not_free) result (free)
!
! This function overloaded by attribute_free. See attribute_free for more details.
!-

function attribute_free3 (ix_ele, ix_branch, attrib_name, lat, &
                                  err_print_flag, except_overlay, dependent_attribs_free, why_not_free) result (free)

implicit none

type (lat_struct) :: lat

integer ix_ele, ix_branch
integer, optional :: why_not_free

character(*) attrib_name

logical free
logical, optional :: err_print_flag, except_overlay, dependent_attribs_free

!

free = check_this_attribute_free (lat%branch(ix_branch)%ele(ix_ele), attrib_name, &
                                  lat, err_print_flag, except_overlay, dependent_attribs_free, why_not_free)

end function attribute_free3

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

function check_this_attribute_free (ele, attrib_name, lat, err_print_flag, except_overlay, &
                                                           dependent_attribs_free, why_not_free) result (free)

implicit none

type (ele_struct), target :: ele
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele_p, lord, slave
type (branch_struct), pointer :: branch
type (ele_attribute_struct) attrib_info
type (control_struct), pointer :: control
type (all_pointer_struct) a_ptr

integer, optional :: why_not_free
integer ix_branch, i, ir, ix_attrib, ix, ic, is, il

character(*) attrib_name
character(40) a_name

logical, optional :: err_print_flag, except_overlay, dependent_attribs_free
logical free, do_print, do_except_overlay, dep_attribs_free, err_flag

!

free = .true.

do_print = logic_option (.true., err_print_flag)
do_except_overlay = logic_option(.false., except_overlay)
dep_attribs_free = logic_option(.false., dependent_attribs_free)

branch => lat%branch(ele%ix_branch)

! Init & check that the name corresponds to an attribute

ix_attrib = attribute_index(ele, attrib_name)
attrib_info = attribute_info(ele, ix_attrib)

a_name = attribute_name (ele, ix_attrib)

if (attrib_info%state == private$) then
  call it_is_not_free (free, ele, ix_attrib, dependent$, 'THIS ATTRIBUTE IS PRIVATE.')
  return
endif

if (attrib_info%state == does_not_exist$) then
  ! Calculated quantities like the Twiss function do not have an entry in the attribute table but
  ! pointer_to_attribute will return a pointer.
  ! Note: Something like beginning element beta_a do have an entry in the attribute table. 
  select case (attrib_name)
  case ('ALPHA_A', 'ALPHA_B', 'BETA_A', 'BETA_B', 'PHI_A', 'PHI_B', 'DPHI_A', 'DPHI_B', &
        'ETA_A', 'ETAP_A', 'ETA_B', 'ETAP_B')
    call it_is_not_free (free, ele, ix_attrib, dependent$, 'THIS ATTRIBUTE IS A COMPUTED PARAMETER.')
  case default
    ! Something like 'cartesian_map(1)%field_scale' does not have an attribute index
    call pointer_to_attribute (ele, attrib_name, .true., a_ptr, err_flag, .false.)
    if (.not. err_flag) return
    call it_is_not_free (free, ele, ix_attrib, does_not_exist$, 'DOES NOT CORRESPOND TO A VALID ATTRIBUTE.', skip = .true.)
  end select
  return
endif

! Overlay or group lord check

if (ele%key == overlay$ .or. ele%key == group$) then
  if (attrib_name == 'IS_ON') return
  if (all(attrib_name /= ele%control%var%name)) then
    call it_is_not_free (free, ele, ix_attrib, does_not_exist$, 'IS NOT A VALID CONTROL VARIABLE', skip = .true.)
  endif
  return
endif

! Here if checking something that is not an overlay or group lord... 

if (attrib_info%state == dependent$) then
  call it_is_not_free (free, ele, ix_attrib, dependent$, 'THIS IS A DEPENDENT ATTRIBUTE.')
  return
endif

! If the attribute is controled by an overlay lord then it cannot be varied.
! Exception: Multiple overlays can control the same attribute.

if (.not. do_except_overlay) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i, control)
    if (lord%key /= overlay$ .or. .not. lord%is_on) cycle
    if (control%ix_attrib /= ix_attrib) cycle

    call it_is_not_free (free, ele, ix_attrib, overlay_slave$, 'IT IS CONTROLLED BY THE OVERLAY: ' // lord%name)
    return
  enddo
endif

! With a few exceptions, super_slaves attributes cannot be varied.

if (ele%slave_status == super_slave$) then
  select case (a_name)
  case ('L', 'CSR_METHOD', 'SPACE_CHARGE_METHOD')
  case default
    call it_is_not_free (free, ele, ix_attrib, super_slave$, 'THE ELEMENT IS A SUPER_SLAVE.', &
                                         '[ATTRIBUTES OF SUPER_SLAVE ELEMENTS ARE GENERALLY NOT FREE TO VARY.]')
  end select
  return
endif

! Check for a multipass_slave.
! Exception: phi0_multipass can be varied for lcavity and rfcavity slaves, etc.

if (ele%slave_status == multipass_slave$) then
  select case (a_name)
  case ('CSR_METHOD', 'SPACE_CHARGE_METHOD', 'DESCRIP', 'ALIAS', 'TYPE', 'TRACKING_METHOD', &
        'MAT6_CALC_METHOD', 'SPIN_TRACKING_METHOD', 'FIELD_CALC', 'PHI0_MULTIPASS', 'PTC_INTEGRATION_TYPE', &
        'INTEGRATOR_ORDER', 'DS_STEP', 'CSR_DS_STEP', 'NUM_STEPS'); return
  end select

  select case (ele%key)
  case (patch$)
    lord => pointer_to_lord(ele, 1)
  end select

  call it_is_not_free (free, ele, ix_attrib, multipass_slave$, 'THE ELEMENT IS A MULTIPASS_SLAVE.', &
                                       '[ATTRIBUTES OF MULTIPASS_SLAVE ELEMENTS ARE GENERALLY NOT FREE TO VARY.]')

  return
endif

! 

select case (a_name)
case ('NUM_STEPS')
  if (.not. dep_attribs_free) call it_is_not_free (free, ele, ix_attrib, dependent$, 'THIS  IS A DEPENDENT ATTRIBUTE.')
  return

case ('FIELD_SCALE', 'PHI0_FIELDMAP', 'CSR_METHOD', 'SPACE_CHARGE_METHOD')
  return

case ('E_TOT', 'P0C')
  if (ele%key == beginning_ele$) return

  if (.not. dep_attribs_free .and. ele%lord_status == multipass_lord$ .and. &
                      .not. ele%field_master .and. nint(ele%value(multipass_ref_energy$)) == user_set$) then
    select case (ele%key)
    case (quadrupole$, sextupole$, octupole$, thick_multipole$, solenoid$, sol_quad$, sbend$, rf_bend$, &
          hkicker$, vkicker$, kicker$, ac_kicker$, elseparator$, lcavity$, rfcavity$)
      return  ! Is free
    end select
  endif

  call it_is_not_free (free, ele, ix_attrib, dependent$, 'THIS IS A DEPENDENT ATTRIBUTE.')
  return
end select

! A super_lord alignment attribute may not be changed if any super_slave is not an em_field element and
! that super_slave has a second lord that is not a pipe.

if (ele%lord_status == super_lord$) then
  select case (a_name)
  case ('ROLL', 'REF_TILT', 'X_OFFSET', 'Y_OFFSET', 'Z_OFFSET', 'X_PITCH', 'Y_PITCH')
    do is = 1, ele%n_slave
      slave => pointer_to_slave(ele, is)
      do il = 1, slave%n_lord
        lord => pointer_to_lord(slave, il)
        if (lord%lord_status /= super_lord$ .or. lord%ix_ele == ele%ix_ele) cycle
        if (lord%key == pipe$ .or. lord%key == em_field$) cycle
        if ((lord%key == quadrupole$ .or. lord%key == solenoid$) .and. &
             (ele%key == quadrupole$ .or. ele%key == solenoid$)) cycle
        call it_is_not_free(free, ele, ix_attrib, super_lord_align$, &
                     'BMAD CANNOT HANDLE OVERLAPPING SUPER_LORD ELEMENTS THAT HAVE MISALIGNMENTS OTHER THAN TILT.')
      enddo
    enddo
  end select
endif


! Check if it is a dependent variable.

if (attrib_info%state == is_free$) return

select case (ele%key)
case (sbend$, rf_bend$)
  if (any(ix_attrib == [angle$, l_chord$, rho$])) free = .false.
case (rfcavity$)
  if (.not. dep_attribs_free) then
    if (ix_attrib == harmon$) free = .false.
    if (ix_attrib == gradient$) free = .false.
  endif
case (lcavity$)
  if (.not. dep_attribs_free) then
    if (ix_attrib == voltage$) free = .false.
    if (ix_attrib == voltage_err$) free = .false.
  endif
  if (ix_attrib == gradient$ .and. ele%value(l$) == 0) free = .false.
  if (ix_attrib == gradient_err$ .and. ele%value(l$) == 0) free = .false.
case (elseparator$)
  if (ix_attrib == voltage$) free = .false.
end select

if ((ele%key == sbend$ .or. ele%key == rf_bend$) .and. ele%lord_status == multipass_lord$ .and. &
    nint(ele%value(multipass_ref_energy$)) == user_set$ .and. ix_attrib == p0c$) free = .true.

if (.not. free) then
  call it_is_not_free (free, ele, ix_attrib, dependent$, 'THIS IS A DEPENDENT ATTRIBUTE.')
  return
endif

! field_master on means that the b_field and bn_gradient values control the strength.

if (.not. dep_attribs_free) then
  free = field_attribute_free (ele, a_name)
  if (.not. free) then
    call it_is_not_free (free, ele, ix_attrib, field_master_dependent$, &
         "THIS IS A DEPENDENT ATTRIBUTE SINCE", &
         "THE ELEMENT'S FIELD_MASTER IS SET TO: " // on_off_logic (ele%field_master, 'True', 'False'))
    return
  endif
endif

!-------------------------------------------------------
contains

subroutine it_is_not_free (free, ele, ix_attrib, why, l1, l2, skip)

type (ele_struct) ele

integer ix_attrib, why, nl

character(*) l1
character(*), optional :: l2
character(100) li(8)
character(*), parameter :: r_name = 'attribute_free'

logical free
logical, optional :: skip

!

free = .false.
if (present(why_not_free)) why_not_free = why

if (.not. do_print) return

nl = 0

nl=nl+1; li(nl) =   'THE ATTRIBUTE: ' // attrib_name
nl=nl+1; li(nl) =   'OF THE ELEMENT: ' // ele%name
if (.not. logic_option(.false., skip)) then
  nl=nl+1; li(nl) = 'IS NOT FREE TO VARY SINCE:'
endif

nl=nl+1; li(nl) = l1
if (present(l2)) then
  nl=nl+1; li(nl) = l2
endif

call out_io (s_error$, r_name, li(1:nl))   

end subroutine it_is_not_free

end function check_this_attribute_free

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function field_attribute_free (ele, attrib_name) result (free)
!
! Routine to check if a field attribute is free to vary.
!
! Field attributes are either normalized (EG K2 of a sextupole) or unnormalized (EG B2_GRADIENT of a sextupole).
! Whether normalized or unnormalized attributes are free to vary will depend on the setting  of ele%field_master.
!
! Generally, this routine should not be called directly. Use the routine attribute_free instead.
!
! Input:
!   ele           -- ele_struct: Element containing the attribute
!   attrib_name   -- character(*): Name of the field attribute. Assumed upper case.
!
! Output:
!   free          -- logical: Is the attribute free to vary? 
!                     If the attribute is not recognized, free = True will be returned.
!-

function field_attribute_free (ele, attrib_name) result (free)

implicit none

type (ele_struct) ele
character(*) attrib_name
logical free

!

free = .true.

if (ele%field_master) then
  select case (ele%key)
  case (quadrupole$)
    if (attrib_name == 'K1') free = .false.
  case (sextupole$)
    if (attrib_name == 'K2') free = .false.
  case (octupole$)
    if (attrib_name == 'K3') free = .false.
  case (solenoid$)
    if (attrib_name == 'KS') free = .false.
  case (sol_quad$)
    if (attrib_name == 'KS') free = .false.
    if (attrib_name == 'K1') free = .false.
    if (attrib_name == 'K2') free = .false.
  case (sbend$)
    if (attrib_name == 'G') free = .false.
    if (attrib_name == 'DG') free = .false.
  case (rf_bend$)
    if (attrib_name == 'G') free = .false.
  case (hkicker$, vkicker$)
    if (attrib_name == 'KICK') free = .false.
  end select

  if (has_hkick_attributes(ele%key)) then
    if (attrib_name == 'HKICK') free = .false.
    if (attrib_name == 'VKICK') free = .false.
  endif

else
  select case (ele%key)
  case (elseparator$)
    if (attrib_name == 'E_FIELD') free = .false.
  case (quadrupole$)
    if (attrib_name == 'B1_GRADIENT') free = .false.
  case (sextupole$)
    if (attrib_name == 'B2_GRADIENT') free = .false.
  case (octupole$)
    if (attrib_name == 'B3_GRADIENT') free = .false.
  case (solenoid$)
    if (attrib_name == 'BS_FIELD') free = .false.
  case (sol_quad$)
    if (attrib_name == 'BS_FIELD') free = .false.
    if (attrib_name == 'B1_GRADIENT') free = .false.
  case (sbend$)
    if (attrib_name == 'B_FIELD') free = .false.
    if (attrib_name == 'DB_FIELD') free = .false.
    if (attrib_name == 'B1_GRADIENT') free = .false.
    if (attrib_name == 'B2_GRADIENT') free = .false.
  case (rf_bend$)
    if (attrib_name == 'B_FIELD') free = .false.
  case (hkicker$, vkicker$)
    if (attrib_name == 'BL_KICK') free = .false.
  end select

  if (has_hkick_attributes(ele%key)) then
    if (attrib_name == 'BL_HKICK') free = .false.
    if (attrib_name == 'BL_VKICK') free = .false.
  endif
endif

end function field_attribute_free

end module
