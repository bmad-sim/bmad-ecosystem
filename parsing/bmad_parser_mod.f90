!+
! Module bmad_parser_mod
!
! This module is a collection of helper routines used by bmad_parser and bmad_parser2.
! The routines in this module are specifically taylored for bmad_parser and
! bmad_parser2 and cannot, in general, be used otherwise.
!-

module bmad_parser_mod

use ptc_interface_mod, only: set_ptc
use bookkeeper_mod
use wake_mod
use attribute_mod
use superimpose_mod
use track1_mod

implicit none

private parse_grid_field, parse_cylindrical_map, parser_get_integer, parser_get_logical

! A "sequence" is a line or a list.
! The information about a sequence is stored in a seq_struct.

! A seq_struct has an array of seq_ele_struct structures.
! Each seq_ele_struct represents an individual element in a sequence and, 
! since sequences can be nested, can itself be a line or a list.

type seq_ele_struct
  character(40) name                     ! name of element, subline, or sublist
  character(40), pointer :: actual_arg(:) => null()
  character(40) :: tag = ''              ! tag name.
  integer :: type = 0                    ! LINE$, REPLACEMENT_LINE$, LIST$, ELEMENT$
  integer :: ix_ele = 0                  ! if an element: pointer to ELE array
                                         ! if a list: pointer to SEQ array
  integer :: ix_arg  = 0                 ! index in arg list (for replacement lines)
  integer ::rep_count = 1                ! how many copies of an element
  logical :: ele_order_reflect = .false. ! Travel through ele sequence in reverse order
  integer :: ele_orientation = 1         ! Travel through elements in reverse.
end type

type seq_struct
  character(40) name              ! name of sequence
  type (seq_ele_struct), pointer :: ele(:) => null()
  character(40), pointer :: dummy_arg(:) => null()
  character(40), pointer :: corresponding_actual_arg(:) => null()
  integer type                    ! LINE$, REPLACEMENT_LINE$ or LIST$
  integer ix                      ! current index of element in %ele
  integer indexx                  ! alphabetical order sorted index
  character(200) file_name        ! file where sequence is defined
  integer ix_line                 ! line number in filewhere sequence is defined
  logical multipass
  logical ptc_layout              ! Put in separate PTC layout
end type

type used_seq_struct
  character(40) :: name = ''            ! name of sequence or element
  character(40) :: tag = ''             ! tag name.
  integer :: ix_multi = 0               ! Multipass indentifier
  integer :: orientation = 1            ! Element reversed?
  integer :: ix_ele_in_in_lat = -1
end type

! A LIFO stack structure is used in the final evaluation of the line that is
! used to form a lattice

type seq_stack_struct
  integer ix_seq                ! index to seq(:) array
  integer ix_ele                ! index to seq%ele(:) array
  integer rep_count             ! repetition count
  integer ele_order_direction   ! +1 => forwad, -1 => back reflection.
  integer orientation_direction ! +1 => forwad, -1 => back reflection.
  character(40) :: tag = ''
  logical multipass
end type

! A LIFO stack structure is used to hold the list of input lattice files
! that are currently open.

type stack_file_struct
  character(200) :: full_name = ''
  character(200) :: dir = './'
  character(200) parse_line_saved
  integer i_line
  integer f_unit
  logical inline_call_active
end type

! structure for holding the control names and pointers for superimpose and overlay elements

type parser_controller_struct ! For overlays and groups
  character(40) :: name     
  character(40) :: attrib_name
  type (expression_atom_struct), allocatable :: stack(:) ! Arithmetic expression stack
  integer n_stk 
end type

type parser_ele_struct
  type (parser_controller_struct), allocatable :: control(:)
  character(40), allocatable :: field_overlaps(:)
  character(40) :: ref_name = ''
  integer :: ix_ref_multipass = 0              ! multipass index for reference element.
  character(40) :: ele_name = ''               ! For patch.
  character(200) :: lat_file = ''                     ! File where element was defined.
  real(rp) :: offset = 0
  integer ix_line_in_file    ! Line in file where element was defined.
  integer ix_count
  integer ele_pt, ref_pt
  integer indexx
  logical :: superposition_has_been_set = .false.
  logical :: create_jumbo_slave = .false.
  logical :: is_range = .false.               ! For girders
  character(40) :: default_attrib = ''        ! For group/overlay elements: slave attribute 
end type

type parser_lat_struct
  type (parser_ele_struct), allocatable :: ele(:) 
end type

!

integer, parameter :: line$ = 1001, list$ = 1002, element$ = 1003
integer, parameter :: replacement_line$ = 1004
integer, parameter :: def$ = 1, redef$ = 2

!------------------------------------------------
! common stuff

integer, parameter :: n_parse_line = 280

type bp_var_struct
  character(40) name      ! variable name
  real(rp) value          ! variable value
  integer :: indexx = 0   ! variable sort index
end type

type bp_common_struct
  type (stack_file_struct), pointer :: current_file
  type (stack_file_struct), pointer :: calling_file
  type (lat_struct), pointer :: old_lat
  type (bp_var_struct), allocatable :: var(:)   ! variable name
  type (extra_parsing_info_struct) extra
  integer num_lat_files               ! Number of files opened
  integer ivar_tot, ivar_init
  character(200), allocatable :: lat_file_names(:) ! List of all files used to create lat
  character(200) line1_file_name               ! Name of file from which input_line1 was read
  character(200) line2_file_name               ! Name of file from which input_line1 was read
  character(n_parse_line) parse_line
  character(n_parse_line) input_line1          ! For debug messages
  character(n_parse_line) input_line2          ! For debug messages
  character(40) parser_name
  logical :: bmad_parser_calling = .false.              ! used for expand_lattice
  logical error_flag                           ! Set True on error
  logical fatal_error_flag                     ! Set True on fatal (must abort now) error 
  logical input_line_meaningful
  logical do_superimpose
  logical write_digested      ! For bmad_parser
  logical write_digested2     ! For bmad_parser2
  logical :: always_parse = .false. ! For debugging to force parsing
  logical input_from_file     ! Input is from a lattice file?
  logical inline_call_active
  logical :: print_err = .true.  ! Print error messages?
  logical :: use_local_lat_file = .false.
  logical :: used_line_set_by_calling_routine = .false.
end type

!

type (bp_common_struct), save, target :: bp_com

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_set_attribute (how, ele, lat, delim, delim_found, err_flag, pele, check_free, wild_and_key0)
!
! Subroutine used by bmad_parser and bmad_parser2 to get the value of
! an attribute from the input file and set the appropriate value in an element.
!
! This subroutine is not intended for general use.
!
! Input:
!   how           -- Integer: Either def$ if the element is being construct from scratch or
!                      redef$ if the element has already been formed and this is part of a
!                      "ele_name[attrib_name] = value" construct.
!   lat           -- lat_struct: Lattice. Needed if the attribute value is an expression
!                      that uses values of other elements.
!   check_free    -- Logical, optional: If present and True then an error will be generated
!                       if the attribute is not free to vary. Used by bmad_parser2.
!   wild_and_key0 -- Logical, optional: If True (default = False), calling routine is working on
!                       something like "*[tracking_method] = runge_kutta". In this case, 
!                       runge_kutta may not be valid for ele but this is not an error.
!
! Output
!   ele          -- ele_struct: Element whos attribute this is.
!   delim        -- Character(1): Delimiter found where the parsing of the input line stops.
!   delim_found  -- Logical: Delimiter found? False if end of input command.
!   err_flag     -- Logical: Set True if there is a problem parsing the input.
!   pele         -- parser_ele_struct, optional: Structure to hold additional 
!                     information that cannot be stored in the ele argument.
!-

subroutine parser_set_attribute (how, ele, lat, delim, delim_found, err_flag, pele, check_free, wild_and_key0)

use random_mod
use wall3d_mod
       
implicit none

type (lat_struct), target :: lat
type (parser_ele_struct), optional :: pele
type (ele_struct), target ::  ele
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), target, save ::  ele0
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: bele
type (all_pointer_struct), allocatable :: a_ptrs(:)
type (all_pointer_struct) a_ptr
type (wall3d_struct), pointer :: wall3d_arr(:), wall3d
type (wall3d_section_struct), pointer :: section
type (wall3d_vertex_struct), pointer :: v_ptr
type (photon_surface_struct), pointer :: surf
type (cylindrical_map_struct), pointer :: cl_map
type (cartesian_map_term1_struct), pointer :: ct_term
type (cartesian_map_term1_struct), allocatable :: ct_terms(:)
type (grid_field_struct), pointer :: g_field
type (taylor_field_struct), pointer :: t_field
type (cartesian_map_struct), pointer :: ct_map
type (wake_lr_spline_struct), allocatable :: lr_pa_temp(:)
type (wake_lr_spline_struct), pointer :: lr_pa

real(rp) kx, ky, kz, tol, value, coef, r_vec(10), r0(2)
real(rp), pointer :: r_ptr

integer i, i2, j, k, n, nn, ix_word, how, ix_word1, ix_word2, ios, ix, i_out, ix_coef, switch
integer expn(6), ix_attrib, i_section, ix_v, ix_sec, i_ptr, i_term, ib, ie, im
integer ix_bounds(2), iy_bounds(2), i_vec(2), family, n_sec, key

character(40) :: word, str_ix, attrib_word, word2, name
character(40), allocatable :: name_list(:)
character(1) delim, delim1, delim2
character(80) str, err_str, line

logical delim_found, err_flag, logic, set_done, end_of_file, do_evaluate, wild_key0
logical is_attrib, err_flag2, old_style_input
logical, optional :: check_free, wild_and_key0

! Get next WORD.
! If an overlay or group element then word is just an attribute to control
! [except for a "GROUP[COMMAND] = 0.343" redef construct]

err_flag = .true.  ! assume the worst
call get_next_word (word, ix_word, ':, =(){', delim, delim_found, call_check = .true.)

! Taylor

wild_key0 = logic_option(.false., wild_and_key0)

if ((ele%key == taylor$ .or. ele%key == hybrid$) .and. delim == '{' .and. word == '') then

  call get_next_word (word, ix_word, ':, =(){', delim, delim_found, call_check = .true.)
  call match_word (word, ['XX', 'XY', 'XZ', 'YX', 'YY', 'YZ', 'ZX', 'ZY', 'ZZ'], i_out, .true., .false.)
  if (i_out > 0) then
    i_out = i_out + 100 ! Make i_out not in range [1:6]
  else
    read (word, *, iostat = ios) i_out
    if (delim /= ':' .or. ix_word == 0 .or. ios /= 0) then
      call parser_error ('BAD "OUT" COMPONENT: ' // word, 'IN TERM FOR TAYLOR ELEMENT: ' // ele%name)
      return
    endif
  endif

  call evaluate_value (ele%name, coef, lat, delim, delim_found, err_flag, ',|');  if (err_flag) return
  delim2 = delim   ! Save
  expn = 0

  ! Need to check for "{1: xxx |}" case where there are no exponents.
  if (delim2 == '|') then
    call get_next_word (word, ix_word, '}', delim, delim_found)
    ! If there are exponents then rewind the get_next_word call.
    if (ix_word /= 0 .or. delim /= '}') then
      bp_com%parse_line = trim(word) // delim // bp_com%parse_line
      delim = delim2
    endif
  endif

  ! Parse exponents

  do j = 1, 100 
    if (delim == '}') exit
    call parser_get_integer (n, word, ix_word, delim, delim_found, err_flag, 'BAD EXPONENT');   if (err_flag) return
    if (.not. expect_one_of ('} ', .true., ele, delim, delim_found)) return
    if (delim2 == ',') then
      select case (j)
      case (6);      if( .not. expect_one_of ('}', .true., ele, delim, delim_found)) return
      case default;  if (.not. expect_one_of (' ', .true., ele, delim, delim_found)) return
      end select
      expn(j) = n
    else
      ! Where, for example, n = 34, must separate into 3 and 4.
      do
        nn = modulo(n, 10)
        if (nn < 1 .or. nn > 6) then
          call parser_error ('BAD EXPONENT VALUE FOR TAYLOR ELEMENT: ' // ele%name, 'CANNOT PARSE: ' // str)
          return
        endif
        expn(nn) = expn(nn) + 1
        n = (n - nn) / 10
        if (n == 0) exit
      enddo
    endif
  enddo

  call add_this_taylor_term (ele, i_out, coef, expn)

  call get_next_word (word, ix_word, '},', delim, delim_found)
  if (ix_word /= 0 .or. (delim_found .and. delim /= ',')) then
    call parser_error ('BAD TERM ENDING FOR TAYLOR ELEMENT: ' // ele%name, 'CANNOT PARSE: ' // str)
    return
  endif

  return
endif

! overlay or group

if (ele%key == overlay$ .or. ele%key == group$) then
  i = attribute_index(ele, word)       ! general attribute search

  select case (i)
  case (type$, alias$, descrip$)
    call bmad_parser_type_get (ele, word, delim, delim_found)
    err_flag = .false.
    return

  case (var$)
    if (how == redef$ .or. associated(ele%control_var)) then
      call parser_error ('RESETTING VAR = {...} IS NOT PERMITTED', 'FOR: ' // ele%name)
      return
    endif
    call get_overlay_group_names(ele, lat, pele, delim, delim_found, .true.)
    pele%default_attrib = ele%control_var(1)%name
    err_flag = .false.
    return

  case (gang$)
    call get_logical_real ('GANG', ele%value(gang$), err_flag)
    return
  end select

  ! Parse old style control var syntax: "i > num_ele_attrib$" handles accordion_edge for example.

  is_attrib = (attribute_index(0, word) > 0 .or. (ele%key == group$ .and. word == 'COMMAND'))
  if (how == def$ .and. .not. associated(ele%control_var) .and. (i < 1 .or. i > num_ele_attrib$) .and. is_attrib) then 
    allocate (ele%control_var(1))
    if (ele%key == group$) then
      ele%control_var(1)%name = 'COMMAND'
    else
      ele%control_var(1)%name = word
    endif
    pele%default_attrib = word
    i = 1 + var_offset$
  endif

  !

  if (i < 1) then
    if (wild_key0) then
      err_flag = .false.
      return
    endif
    call parser_error ('BAD OVERLAY/GROUP ATTRIBUTE: ' // word, 'FOR: ' // ele%name)
    return
  endif

  value = 0
  if (delim == '=') then  ! value
    call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag)
    if (err_flag) return
  endif

  call pointer_to_indexed_attribute (ele, i, .true., a_ptr, err_flag, .true.)
  a_ptr%r = value

  if (attrib_free_problem(word)) return

  err_flag = .false.
  return
endif

! beam_start and bmad_com element can have attributes that are not part of the element so
! Need to use pointers_to_attribute.

key = ele%key
if (ele%key == def_parameter$ .and. word == 'APERTURE_LIMIT_ON') key = def_bmad_com$
if (ele%key == def_parameter$ .and. word == 'ELECTRIC_DIPOLE_MOMENT') key = def_bmad_com$
if (ele%key == def_parameter$ .and. word == 'PTC_CUT_FACTOR') key = def_bmad_com$
if (ele%key == def_parameter$ .and. word == 'USE_HARD_EDGE_DRIFTS') key = def_bmad_com$

if (key == def_beam_start$ .or. key == def_bmad_com$) then
  name = ele%name
  if (ele%name == 'PARAMETER') name = 'BMAD_COM'

  if (delim == '(') then
    ix = index(bp_com%parse_line, '=')
    if (ix == 0) then
      call parser_error ('MALFORMED BMAD_COM OR PARAMETER SET')
      return
    endif
    word = trim(word) // '(' // bp_com%parse_line(:ix-1)
    delim = '='
    call string_trim(bp_com%parse_line(ix+1:), bp_com%parse_line, ix)
  endif

  call pointers_to_attribute (lat, name, word, .false., a_ptrs, err_flag, .false.)
  if (err_flag .or. size(a_ptrs) /= 1) then
    call parser_error ('BAD ATTRIBUTE: ' // word, 'FOR ELEMENT: ' // ele%name)
    return
  endif

  if (ele%key == def_parameter$ .and. word == 'APERTURE_LIMIT_ON') then
    call parser_error ('SYNTAX HAS CHANGED: PARAMETER[APERTURE_LIMIT_ON] = ... NEEDS TO BE REPLACED BY BMAD_COM[APERTURE_LIMIT_ON] = ...', &
                       'THIS IS A WARNING ONLY. THE PROGRAM WILL RUN NORMALLY.', level = s_warn$)
  endif

  if (associated(a_ptrs(1)%r)) then
    call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag) 
    if (err_flag) return
    a_ptrs(1)%r = value
    if (associated(a_ptrs(1)%r, bmad_com%max_aperture_limit))             bp_com%extra%max_aperture_limit_set              = .true.
    if (associated(a_ptrs(1)%r, bmad_com%default_ds_step))                bp_com%extra%default_ds_step_set                 = .true.
    if (associated(a_ptrs(1)%r, bmad_com%significant_length))             bp_com%extra%significant_length_set              = .true.
    if (associated(a_ptrs(1)%r, bmad_com%rel_tol_tracking))               bp_com%extra%rel_tol_tracking_set                = .true.
    if (associated(a_ptrs(1)%r, bmad_com%abs_tol_tracking))               bp_com%extra%abs_tol_tracking_set                = .true.
    if (associated(a_ptrs(1)%r, bmad_com%rel_tol_adaptive_tracking))      bp_com%extra%rel_tol_adaptive_tracking_set       = .true.
    if (associated(a_ptrs(1)%r, bmad_com%abs_tol_adaptive_tracking))      bp_com%extra%abs_tol_adaptive_tracking_set       = .true.
    if (associated(a_ptrs(1)%r, bmad_com%init_ds_adaptive_tracking))      bp_com%extra%init_ds_adaptive_tracking_set       = .true.
    if (associated(a_ptrs(1)%r, bmad_com%min_ds_adaptive_tracking))       bp_com%extra%min_ds_adaptive_tracking_set        = .true.
    if (associated(a_ptrs(1)%r, bmad_com%fatal_ds_adaptive_tracking))     bp_com%extra%fatal_ds_adaptive_tracking_set      = .true.
    if (associated(a_ptrs(1)%r, bmad_com%electric_dipole_moment))         bp_com%extra%electric_dipole_moment_set          = .true.
    if (associated(a_ptrs(1)%r, bmad_com%ptc_cut_factor))                 bp_com%extra%ptc_cut_factor_set                  = .true.
    if (associated(a_ptrs(1)%r, bmad_com%sad_eps_scale))                  bp_com%extra%sad_eps_scale_set                   = .true.
    if (associated(a_ptrs(1)%r, bmad_com%sad_amp_max))                    bp_com%extra%sad_amp_max_set                     = .true.
    if (name(1:5) == 'D_ORB')                                             bp_com%extra%d_orb_set                           = .true.

  elseif (associated(a_ptrs(1)%i)) then
    call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag) 
    if (err_flag) return
    a_ptrs(1)%i = nint(value)
    if (associated(a_ptrs(1)%i, bmad_com%taylor_order))                   bp_com%extra%taylor_order_set                    = .true.
    if (associated(a_ptrs(1)%i, bmad_com%default_integ_order))            bp_com%extra%default_integ_order_set             = .true.
    if (associated(a_ptrs(1)%i, bmad_com%ptc_max_fringe_order))           bp_com%extra%ptc_max_fringe_order_set            = .true.
    if (associated(a_ptrs(1)%i, bmad_com%runge_kutta_order))              bp_com%extra%runge_kutta_order_set               = .true.
    if (associated(a_ptrs(1)%i, bmad_com%sad_n_div_max))                  bp_com%extra%sad_n_div_max_set                   = .true.
    if (associated(a_ptrs(1)%i, bmad_com%max_num_runge_kutta_step))       bp_com%extra%max_num_runge_kutta_step_set        = .true.

  elseif (associated(a_ptrs(1)%l)) then
    call get_logical (trim(ele%name) // ' ' // word, a_ptrs(1)%l, err_flag)
    if (err_flag) return
    if (associated(a_ptrs(1)%l, bmad_com%use_hard_edge_drifts))           bp_com%extra%use_hard_edge_drifts_set            = .true.
    if (associated(a_ptrs(1)%l, bmad_com%sr_wakes_on))                    bp_com%extra%sr_wakes_on_set                     = .true.
    if (associated(a_ptrs(1)%l, bmad_com%lr_wakes_on))                    bp_com%extra%lr_wakes_on_set                     = .true.
    if (associated(a_ptrs(1)%l, bmad_com%mat6_track_symmetric))           bp_com%extra%mat6_track_symmetric_set            = .true.
    if (associated(a_ptrs(1)%l, bmad_com%auto_bookkeeper))                bp_com%extra%auto_bookkeeper_set                 = .true.
    if (associated(a_ptrs(1)%l, bmad_com%space_charge_on))                bp_com%extra%space_charge_on_set                 = .true.
    if (associated(a_ptrs(1)%l, bmad_com%coherent_synch_rad_on))          bp_com%extra%coherent_synch_rad_on_set           = .true.
    if (associated(a_ptrs(1)%l, bmad_com%spin_tracking_on))               bp_com%extra%spin_tracking_on_set                = .true.
    if (associated(a_ptrs(1)%l, bmad_com%radiation_damping_on))           bp_com%extra%radiation_damping_on_set            = .true.
    if (associated(a_ptrs(1)%l, bmad_com%radiation_fluctuations_on))      bp_com%extra%radiation_fluctuations_on_set       = .true.
    if (associated(a_ptrs(1)%l, bmad_com%conserve_taylor_maps))           bp_com%extra%conserve_taylor_maps_set            = .true.
    if (associated(a_ptrs(1)%l, bmad_com%absolute_time_tracking_default)) bp_com%extra%absolute_time_tracking_default_set  = .true.
    if (associated(a_ptrs(1)%l, bmad_com%convert_to_kinetic_momentum))    bp_com%extra%convert_to_kinetic_momentum_set     = .true.
    if (associated(a_ptrs(1)%l, bmad_com%aperture_limit_on))              bp_com%extra%aperture_limit_on_set               = .true.
    if (associated(a_ptrs(1)%l, bmad_com%debug))                          bp_com%extra%debug_set                           = .true.

  else
    call parser_error ('BOOKKEEPING ERROR. PLEASE CONTACT A BMAD MAINTAINER!')
  endif

  return
endif

! Redefinition of Long-range "wake()", "r_custom()", etc.
! Old style wiggler "term()" handled below.

if (delim == '(' .and. .not. (word == 'TERM' .and. how == def$)) then
  word2 = trim(word) // '('
  call get_next_word (word, ix_word, '=', delim, delim_found)

  if (.not. delim_found) then
    call parser_error ('NO "=" SIGN FOUND', 'FOR ELEMENT: ' // ele%name)
    return
  endif

  call pointer_to_attribute (ele, trim(word2) // word, how == def$, a_ptr, err_flag, .false.)

  if (err_flag .or. .not. associated(a_ptr%r)) then
    call parser_error ('BAD ATTRIBUTE: ' // word, 'FOR ELEMENT: ' // ele%name)
    return
  endif

  call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag)
  a_ptr%r = value
  return
endif

! "WALL%" redef

if (word(1:5) == 'WALL%') then

  select case (word(6:))

  ! Section redef

  case ('SECTION')

    if (delim /= '(') then
      call parser_error ('MALFORMED WALL COMPONENT REDEF IN ELEMENT: ' // ele%name)
      return
    endif

    ix_sec = evaluate_array_index (err_flag, ')', word2, '(=', delim)
    n_sec = -1
    if (associated(ele%wall3d)) n_sec = size(ele%wall3d(1)%section)
    if (err_flag .or. ix_sec < 0 .or. ix_sec > n_sec) then
      call parser_error('BAD ' // trim(word) // ' INDEX', 'FOR ELEMENT: ' // ele%name)
      return
    endif
    section => ele%wall3d(1)%section(ix_sec)

    if (word2 == '%S' .and. delim == '=') then
      r_ptr => section%s

    elseif (word2 == '%DR_DS' .and. delim == '=') then
      r_ptr => section%dr_ds

    elseif (word2 == '%V' .and. delim == '(') then
      ix_v = evaluate_array_index (err_flag, ')', word, '=', delim)
      if (err_flag .or. ix_v < 0 .or. ix_v > size(section%v)) then
        call parser_error('BAD VERTEX INDEX',  'FOR ELEMENT: ' // ele%name)
        return
      endif
      v_ptr => section%v(ix_v)

      select case (word)
      case ('%X');        r_ptr => v_ptr%x
      case ('%Y');        r_ptr => v_ptr%y
      case ('%RADIUS_X'); r_ptr => v_ptr%radius_x
      case ('%RADIUS_Y'); r_ptr => v_ptr%radius_y
      case ('%TILT');     r_ptr => v_ptr%tilt
      case default
        call parser_error('BAD WALL SECTION VERTEX COMPONENT: ' // word, 'FOR ELEMENT: ' // ele%name)
        return
      end select
    else
      call parser_error('BAD WALL SECTION COMPONENT: ' // word2, 'FOR ELEMENT: ' // ele%name)
      return
    endif

    call evaluate_value (ele%name, r_ptr, lat, delim, delim_found, err_flag)

  ! Not recognized

  case default
    call parser_error ('BAD WALL COMPONENT REDEF: ' // word, 'IN ELEMENT: ' // ele%name)
  end select

  return
endif

! if not an overlay/group then see if it is an ordinary attribute.
! if not an ordinary attribute then might be a superimpose switch

if (ix_word == 0) then  ! no word
  call parser_error  ('"," NOT FOLLOWED BY ATTRIBUTE NAME FOR: ' // ele%name)
  return
endif


select case (word)
case ('ELE_BEGINNING', 'ELE_CENTER', 'END_END', 'REF_BEGINNING', 'REF_CENTER', 'REF_END')
  call parser_error ('OLD SUPERPOSITION SYNTAX: ' // word, &
              'PLEASE CONVERT (SEE THE BMAD MANUAL)', 'WARNING ONLY, PROGRAM WILL RUN NORMALLY...', level = s_warn$)
end select

select case (word)
case ('TILT')
  if (ele%key == sbend$ .or. ele%key == rbend$) then
    call parser_error ('BENDS HAVE A "REF_TILT" ATTRIBUTE BUT NOT A "TILT" ATTRIBUTE.')
  endif

case ('DPHI0')
  call parser_error ('THE ATTRIBUTE NAME "DPHI0" HAS BEEN CHANGED TO "PHI0_MULTIPASS"', &
                     'PLEASE MAKE THE CHANGE IN THE LATTICE FILE.', &
                     '[THIS IS A WARNING ONLY. THIS PROGRAM WILL RUN NORMALLY]', level = s_warn$)
  word = 'PHI0_MULTIPASS'

case ('REL_TRACKING_CHARGE') 
  call parser_error ('THE ATTRIBUTE NAME "REL_TRACKING_CHARGE" HAS BEEN CHANGED TO "DEFAULT_TRACKING_SPECIES"', &
                     'PLEASE MAKE THE CHANGE IN THE LATTICE FILE.')

case ('RADIUS')
  call parser_error ('THE ATTRIBUTE "RADIUS" HAS BEEN CHANGED TO "R0_MAG"', &
                     'PLEASE MAKE THE CHANGE IN THE LATTICE FILE.', &
                     '[THIS IS A WARNING ONLY. THIS PROGRAM WILL RUN NORMALLY]', level = s_warn$)
  word = 'R0_MAG'

case ('FIELD')
  call parser_error ('THE "FIELD = {...}" SYNTAX HAS BEEN CHANGED TO "GRID_FIELD = {...} AND/OR CYLINDRICAL_MAP = {...}.', &
                     'NOTE: THIS INVOLVES MORE THAN CHANGING "FIELD" TO "GRID_FIELD" OR "CYLINDRICAL_MAP".', &
                     'PLEASE SEE THE BMAD MANUAL FOR MORE DETAILS')

case ('REF_BEGINNING')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  pele%ref_pt = anchor_beginning$
  err_flag = .false.
  return

case ('REF_CENTER')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  pele%ref_pt = anchor_center$
  err_flag = .false.
  return

case ('REF_END')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  pele%ref_pt = anchor_end$
  err_flag = .false.
  return

case ('ELE_BEGINNING')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  pele%ele_pt = anchor_beginning$
  err_flag = .false.
  return

case ('ELE_CENTER')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  pele%ele_pt = anchor_center$
  err_flag = .false.
  return

case ('ELE_END')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  pele%ele_pt = anchor_end$
  err_flag = .false.
  return
end select

word = parser_translate_attribute_name (ele%key, word)

ix_attrib = attribute_index(ele, word, attrib_word)
if (attrib_free_problem(attrib_word)) return

if (ix_attrib < 1) then
  if (ele%key == drift$ .and. (word == 'HKICK' .or. word == 'VKICK' .or. &
        word == 'BL_HKICK' .or. word == 'BL_VKICK')) then
    call parser_error ('BAD ATTRIBUTE: ' // word, 'FOR ELEMENT: ' // ele%name, &
                      'ONE SOLUTION IS TO MAKE THIS DRIFT A "PIPE" ELEMENT.')
  else
    call parser_error ('BAD ATTRIBUTE NAME: ' // word, 'FOR ELEMENT: ' // ele%name)
  endif
  return
endif

! wall cross-section definition

if (attrib_word == 'WALL') then

  i_section = 0
  if (.not. expect_this ('=', .true., .true., 'AFTER "WALL"', ele, delim, delim_found)) return
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)

  ! "ele1[wall] = ele2[wall]" construct

  if (delim == '[') then
    call parser_find_ele_for_attrib_transfer ('WALL')
    if (.not. associated(eles(1)%ele%wall3d)) then
      call parser_error ('NO WALL ASSOCIATED WITH LATTICE ELEMENT: ' // word)
      return
    endif
    call transfer_wall3d(eles(1)%ele%wall3d, ele%wall3d)
    return
  endif

  !

  if (.not. expect_this ('{', .true., .true., 'AFTER "WALL"', ele, delim, delim_found)) return

  ! Loop wall3d_struct components.

  if (associated (ele%wall3d)) then
    n = size(ele%wall3d)
    wall3d_arr => ele%wall3d
    allocate (ele%wall3d(n+1))
    do i = 1, n
      wall3d => ele%wall3d(i)
      n_sec = size(wall3d_arr(i)%section)
      allocate(wall3d%section(n_sec))
      do ix_sec = 1, n_sec
        nn = size(wall3d%section(ix_sec)%v)
        allocate(wall3d%section(ix_sec)%v(nn))
      enddo
      wall3d = wall3d_arr(i)
      wall3d%n_link = 1
    enddo
    call unlink_wall3d (wall3d_arr)
    wall3d => ele%wall3d(n+1)
  else
    allocate (ele%wall3d(1))
    wall3d => ele%wall3d(1)
  endif

  wall3d_loop: do

    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

    ! Possible "}" is end of wall 
    if (delim /= '}' .and. word == '') exit
    if (.not. expect_this ('=', .true., .false., 'AFTER ' // trim(word) // ' IN WALL CONSTRUCT', ele, delim, delim_found)) return

    select case (word)

    case ('NAME')
      call bmad_parser_type_get (ele, word, delim, delim_found, str_out = wall3d%name)

    case ('OPAQUE_MATERIAL') 
      call bmad_parser_type_get (ele, word, delim, delim_found, str_out = wall3d%opaque_material)

    case ('CLEAR_MATERIAL') 
      call bmad_parser_type_get (ele, word, delim, delim_found, str_out = wall3d%clear_material)

    case ('THICKNESS') 
      call evaluate_value (ele%name, wall3d%thickness, lat, delim, delim_found, err_flag, ',}')
      if (err_flag) return

    case ('ELE_ANCHOR_PT')
      call get_switch ('WALL ELE_ANCHOR_PT', anchor_pt_name(1:), wall3d%ele_anchor_pt, err_flag2, ele, delim, delim_found)
      if (err_flag2) return

    case ('SUPERIMPOSE')
      call get_logical ('WALL SUPERIMPOSE', wall3d%superimpose, err_flag2); if (err_flag2) return

    ! Must be "section = {"

    case ('SECTION')

      ! Read in section

      if (.not. expect_this ('{', .false., .true., 'AFTER "SECTION =" IN WALL CONSTRUCT', ele, delim, delim_found)) return

      i_section = i_section + 1
      ix_v = 0
      call re_allocate (wall3d%section, i_section)
      section => wall3d%section(i_section)

      wall3d_section_loop: do

        call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

        ! Possible "}" is end of wall 
        if (delim /= '}' .and. word == '') exit
        if (word == 'V') then
          if (.not. expect_this ('(', .true., .false., 'AFTER ' // trim(word) // ' IN WALL CONSTRUCT', ele, delim, delim_found)) return
        else
          if (.not. expect_this ('=', .true., .false., 'AFTER ' // trim(word) // ' IN WALL CONSTRUCT', ele, delim, delim_found)) return
        endif

        select case (word)

        case ('TYPE') 
          call get_switch ('WALL TYPE', wall3d_section_type_name(1:), section%type, err_flag2, ele, delim, delim_found)
          if (err_flag2) return

        case ('MATERIAL') 
          call bmad_parser_type_get (ele, word, delim, delim_found, str_out = section%material)

        case ('THICKNESS')
          call evaluate_value (trim(ele%name) // ' ' // word, section%thickness, lat, delim, delim_found, err_flag, ',}')
          if (err_flag) return
          if (ele%key == capillary$) ele%value(l$) = section%s

        case ('S')
          call evaluate_value (trim(ele%name) // ' ' // word, section%s, lat, delim, delim_found, err_flag, ',}')
          if (err_flag) return
          if (ele%key == capillary$) ele%value(l$) = section%s

        case ('DR_DS') 
          call evaluate_value (trim(ele%name) // ' ' // word, section%dr_ds, lat, delim, delim_found, err_flag, ',}')
          if (err_flag) return
                  
        case ('ABSOLUTE_VERTICES') 
          call get_logical (trim(ele%name) // ' ' // word, section%absolute_vertices_input, err_flag)
          if (err_flag) return

        case ('X0') 
          call evaluate_value (trim(ele%name) // ' ' // word, section%r0(1), lat, delim, delim_found, err_flag, ',}')
          if (err_flag) return

        case ('Y0') 
          call evaluate_value (trim(ele%name) // ' ' // word, section%r0(2), lat, delim, delim_found, err_flag, ',}')
          if (err_flag) return

        case ('R0')
          if (.not. parse_real_list (lat, trim(ele%name) // ' GRID R0', section%r0, .true.)) return

        ! Parse "V() = ..." constructs.

        case ('V')

          ix_v = ix_v + 1
          section%n_vertex_input = ix_v
          call re_allocate (section%v, ix_v)

          call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
          read (word, *, iostat = ios) j 
          if (ios /= 0 .or. ix_v /= j) then
            call parser_error ('BAD OR OUT OF ORDER WALL SECTION VERTEX INDEX NUMBER FOR: ' // ele%name)
            return
          endif

          if (.not. expect_this (')={', .true., .false., 'AFTER "V(n)" IN WALL CONSTRUCT', ele, delim, delim_found)) return

          call evaluate_value (trim(ele%name), section%v(ix_v)%x, lat, delim, delim_found, err_flag, ',')
          if (err_flag) return

          call evaluate_value (trim(ele%name), section%v(ix_v)%y, lat, delim, delim_found, err_flag, ',}')
          if (err_flag) return

          if (delim == ',') then
            call evaluate_value (trim(ele%name), section%v(ix_v)%radius_x, lat, delim, delim_found, err_flag, ',}')
            if (err_flag) return
          endif

          if (delim == ',') then
            call evaluate_value (trim(ele%name), section%v(ix_v)%radius_y, lat, delim, delim_found, err_flag, ',}')
            if (err_flag) return
          endif

          if (delim == ',') then
            call evaluate_value (trim(ele%name), section%v(ix_v)%tilt, lat, delim, delim_found, err_flag, '}')
            if (err_flag) return
          endif

          call get_next_word (word, ix_word, '{},()=', delim, delim_found)
          if (word /= '' .or. (delim /= '}' .and. delim /= ',')) then
            call parser_error ('BAD SYNTAX IN WALL SECTION DEFINITION FOR ELEMENT: ' // ele%name)
            return
          endif

        case default
          call parser_error ('WALL SECTION COMPONENT NOT RECOGNIZED: ' // word, 'FOR ELEMENT: ' // ele%name)
          return
        end select   ! section components

        if (.not. expect_one_of (',}', .true., ele, delim, delim_found)) return
        if (delim == '}') then
          if (.not. expect_one_of(',}', .false., ele, delim, delim_found)) return
          exit
        endif
      enddo wall3d_section_loop

    case default
      call parser_error ('WALL COMPONENT NOT RECOGNIZED: ' // word, 'FOR ELEMENT: ' // ele%name)
      return
    end select   ! wall components

    if (.not. expect_one_of (',}', .true., ele, delim, delim_found)) return
    if (delim == '}') exit

  enddo wall3d_loop

  ! Next thing on line should be either a "," or end-of-line

  logic = expect_one_of(', ', .false., ele, delim, delim_found)
  return

endif

!-------------------------------
! Reflecting Surface

if (attrib_word == 'SURFACE') then
  surf => ele%photon%surface

  if (.not. expect_this ('=', .true., .true., 'AFTER "SURFACE"', ele, delim, delim_found)) return
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)

  ! "ele1[surface] = ele2[surface]" construct

  if (delim == '[') then
    call parser_find_ele_for_attrib_transfer ('SURFACE')
    if (.not. associated(eles(1)%ele%photon)) then
      call parser_error ('NO SURFACE ASSOCIATED WITH LATTICE ELEMENT: ' // word)
      return
    endif
    ele%photon%surface = eles(1)%ele%photon%surface
    return
  endif

  !

  if (.not. expect_this ('{', .true., .true., 'AFTER "SURFACE"', ele, delim, delim_found)) return

  surface_loop: do

    ! Expect "GRID ={" or "TYPE ="

    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

    select case (word)

    case ('GRID')

      if (.not. expect_this ('={', .true., .true., 'AFTER "GRID"', ele, delim, delim_found)) return
      ix_bounds = int_garbage$; iy_bounds = int_garbage$

      do
        call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
        if (word /= 'PT') then
          if (.not. expect_this ('=', .true., .false., 'AFTER ' // trim(word) // ' IN SURFACE CONSTRUCT', ele, delim, delim_found)) return
        endif

        select case (word)
        case ('DR')
          if (.not. parse_real_list (lat, trim(ele%name) // ' GRID DR', surf%grid%dr, .true.)) return

        case ('R0')
          if (.not. parse_real_list (lat, trim(ele%name) // ' GRID R0', surf%grid%r0, .true.)) return

        case ('IX_BOUNDS', 'IY_BOUNDS')
          if (.not. parse_integer_list (trim(ele%name) // ' GRID ' // trim(word), lat, i_vec, .true.)) return
          if (word == 'IX_BOUNDS') ix_bounds = i_vec
          if (word == 'IY_BOUNDS') iy_bounds = i_vec

          if (any(ix_bounds /= int_garbage$) .and. any(iy_bounds /= int_garbage$)) then
            if (any(ix_bounds == int_garbage$) .or. any(iy_bounds == int_garbage$) .or. &
                ix_bounds(1) > ix_bounds(2) .or. iy_bounds(1) > iy_bounds(2)) then
              call parser_error ('SURFACE GRID X/IY_BOUNDS NOT PROPERLY SET', trim(ele%name))
              return
            endif
            if (allocated (surf%grid%pt)) deallocate (surf%grid%pt)
            allocate (surf%grid%pt(ix_bounds(1):ix_bounds(2), iy_bounds(1):iy_bounds(2)))
          endif

        case ('PT')
          bp_com%parse_line = delim // bp_com%parse_line
          if (.not. parse_integer_list (trim(ele%name) // ' GRID PT', lat, i_vec, .true.)) return
          if (.not. allocated(surf%grid%pt)) then
            call parser_error ('SURFACE PT_MAX MISSING', 'FOR: ' // ele%name)
            return
          endif
          if (any(i_vec < 0) .or. any(i_vec > ubound(surf%grid%pt))) then
            call parser_error ('SURFACE PT(I,J) INDEX OUT OF BOUNDS', 'FOR: ' // ele%name)
            return
          endif
          if (.not. expect_this ('=', .false., .false., 'IN GRID PT', ele, delim, delim_found)) return
          if (.not. parse_real_list (lat, trim(ele%name) // 'IN GRID PT', r_vec(1:4), .true.)) return
          surf%grid%pt(i_vec(1), i_vec(2))%orientation = surface_orientation_struct(r_vec(1), r_vec(2), r_vec(3), r_vec(4))

        case ('TYPE')
          call get_switch ('SURFACE GRID TYPE', surface_grid_type_name(1:), surf%grid%type, err_flag2, ele, delim, delim_found)
          if (err_flag2) return
          bp_com%parse_line = delim // bp_com%parse_line

        case default
          call parser_error ('GRID COMPONENT NOT RECOGNIZED: ' // word, 'FOR ELEMENT: ' // ele%name)
          return
        end select

        if (.not. expect_one_of (',}', .false., ele, delim, delim_found)) return
        if (delim == '}') then
          if (.not. expect_one_of (',}', .false., ele, delim, delim_found)) return
          exit
        endif

      enddo

    case default
      call parser_error ('UNKNOWN SURFACE COMPONENT: ' // word2, 'FOR: ' // ele%name)
      return
    end select

    if (delim == '}') exit

  enddo surface_loop

  if (.not. expect_one_of(', ', .false., ele, delim, delim_found)) return
  err_flag = .false.
  return

endif

!-------------------------------
! Cartesian_map field

if (attrib_word == 'CARTESIAN_MAP') then

  if (.not. expect_this ('=', .true., .true., 'AFTER "CARTESIAN_MAP"', ele, delim, delim_found)) return
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)

  ! "ele1[cartesian_map] = ele2[cartesian_map]" construct

  if (delim == '[') then
    call parser_find_ele_for_attrib_transfer ('CARTESIAN_MAP')
    if (err_flag) return
    if (.not. associated(eles(1)%ele%cartesian_map)) then
      call parser_error ('NO CARTESIAN_MAP ASSOCIATED WITH LATTICE ELEMENT: ' // word)
      return
    endif
    call transfer_fieldmap(eles(1)%ele, ele, cartesian_map$)
    return
  endif

  if (.not. expect_this ('{', .true., .true., 'AFTER "CARTESIAN_MAP"', ele, delim, delim_found)) return

  if (associated(ele%cartesian_map)) then
    i_ptr = size(ele%cartesian_map) + 1
    ele0%cartesian_map => ele%cartesian_map
    allocate(ele%cartesian_map(i_ptr))
    do i = 1, i_ptr-1
     ele%cartesian_map(i)%ptr => ele0%cartesian_map(i)%ptr
    enddo
  else
    allocate(ele%cartesian_map(1))
    i_ptr = 1
  endif

  allocate (ele%cartesian_map(i_ptr)%ptr)
  call parse_cartesian_map(ele%cartesian_map(i_ptr), ele, lat, delim, delim_found, err_flag)
  return
endif

!-------------------------------
! Cylindrical_map field

if (attrib_word == 'CYLINDRICAL_MAP') then

  if (.not. expect_this ('=', .true., .true., 'AFTER "CYLINDRICAL_MAP"', ele, delim, delim_found)) return
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)

  ! "ele1[cylindrical_map] = ele2[cylindrical_map]" construct

  if (delim == '[') then
    call parser_find_ele_for_attrib_transfer ('CYLINDRICAL_MAP')
    if (err_flag) return
    if (.not. associated(eles(1)%ele%cylindrical_map)) then
      call parser_error ('NO CYLINDRICAL_MAP ASSOCIATED WITH LATTICE ELEMENT: ' // word)
      return
    endif
    call transfer_fieldmap(eles(1)%ele, ele, cylindrical_map$)
    return
  endif

  if (.not. expect_this ('{', .true., .true., 'AFTER "CYLINDRICAL_MAP"', ele, delim, delim_found)) return

  if (associated(ele%cylindrical_map)) then
    i_ptr = size(ele%cylindrical_map) + 1
    ele0%cylindrical_map => ele%cylindrical_map
    allocate(ele%cylindrical_map(i_ptr))
    do i = 1, i_ptr-1
     ele%cylindrical_map(i)%ptr => ele0%cylindrical_map(i)%ptr
    enddo
  else
    allocate(ele%cylindrical_map(1))
    i_ptr = 1
  endif

  allocate (ele%cylindrical_map(i_ptr)%ptr)
  cl_map => ele%cylindrical_map(i_ptr)
  if (ele%key == lcavity$ .or. ele%key == rfcavity$) cl_map%harmonic = 1 ! Default

  call parse_cylindrical_map(cl_map, ele, lat, delim, delim_found, err_flag)
  return
endif

!-------------------------------
! grid_field field

if (attrib_word == 'GRID_FIELD') then

  if (.not. expect_this ('=', .true., .true., 'AFTER "GRID_FIELD"', ele, delim, delim_found)) return
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)

  ! "ele1[grid_field] = ele2[grid_field]" construct

  if (delim == '[') then
    call parser_find_ele_for_attrib_transfer ('GRID_FIELD')
    if (err_flag) return
    if (.not. associated(eles(1)%ele%grid_field)) then
      call parser_error ('NO GRID_FIELD ASSOCIATED WITH LATTICE ELEMENT: ' // word)
      return
    endif
    call transfer_fieldmap(eles(1)%ele, ele, grid_field$)
    return
  endif

  if (.not. expect_this ('{', .true., .true., 'AFTER "GRID_FIELD"', ele, delim, delim_found)) return

  if (associated(ele%grid_field)) then
    i_ptr = size(ele%grid_field) + 1
    ele0%grid_field => ele%grid_field
    allocate(ele%grid_field(i_ptr))
    do i = 1, i_ptr-1
     ele%grid_field(i) = ele0%grid_field(i)
    enddo
    deallocate (ele0%grid_field)
  else
    allocate(ele%grid_field(1))
    i_ptr = 1
  endif

  allocate (ele%grid_field(i_ptr)%ptr)
  g_field => ele%grid_field(i_ptr)
  if (ele%key == lcavity$ .or. ele%key == rfcavity$) g_field%harmonic = 1 ! Default

  call parse_grid_field(g_field, ele, lat, delim, delim_found, err_flag)
  return
endif

!-------------------------------
! Taylor_field field

if (attrib_word == 'TAYLOR_FIELD') then

  if (.not. expect_this ('=', .true., .true., 'AFTER "TAYLOR_FIELD"', ele, delim, delim_found)) return
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)

  ! "ele1[taylor_field] = ele2[taylor_field]" construct

  if (delim == '[') then
    call parser_find_ele_for_attrib_transfer ('TAYLOR_FIELD')
    if (err_flag) return
    if (.not. associated(eles(1)%ele%taylor_field)) then
      call parser_error ('NO TAYLOR_FIELD ASSOCIATED WITH LATTICE ELEMENT: ' // word)
      return
    endif
    call transfer_fieldmap(eles(1)%ele, ele, taylor_field$)
    return
  endif

  if (.not. expect_this ('{', .true., .true., 'AFTER "TAYLOR_FIELD"', ele, delim, delim_found)) return

  if (associated(ele%taylor_field)) then
    i_ptr = size(ele%taylor_field) + 1
    ele0%taylor_field => ele%taylor_field
    allocate(ele%taylor_field(i_ptr))
    do i = 1, i_ptr-1
     ele%taylor_field(i)%ptr => ele0%taylor_field(i)%ptr
    enddo
  else
    allocate(ele%taylor_field(1))
    i_ptr = 1
  endif

  allocate (ele%taylor_field(i_ptr)%ptr)
  t_field => ele%taylor_field(i_ptr)

  call parse_taylor_field(t_field, ele, lat, delim, delim_found, err_flag)
  return
endif

!------------------------------
! wiggler term attribute

if (ix_attrib == term$ .and. (ele%key == wiggler$ .or. ele%key == undulator$)) then

  err_flag = .true. ! assume the worst

  if (delim /= '(') then   ! ) then
    call parser_error ('"TERM" FOR A WIGGLER NOT FOLLOWED BY A "(" FOR: ' // ele%name)  ! )
    return
  endif

  call parser_get_integer (ix, word, ix_word, delim, delim_found, err_flag, 'BAD WIGGLER "TERM(IX)" CONSTRUCT'); if (err_flag) return

  if (delim /= ')') then
    call parser_error ('CANNOT FIND CLOSING ")" for a "TERM(i)" FOR A WIGGLER"', 'FOR: ' // ele%name)
    return
  endif

  write (str_ix, '(a, i3, a)') 'TERM(', ix, ')'

  if (.not. associated(ele%cartesian_map)) then
    allocate(ele%cartesian_map(1))
    ct_map => ele%cartesian_map(1)
    allocate(ct_map%ptr)
    allocate(ct_map%ptr%term(ix))
    ct_map%ptr%file = bp_com%line1_file_name
    ct_map%master_parameter = polarity$
  else
    ct_map => ele%cartesian_map(1)
    if (ix > size(ct_map%ptr%term)) then
      call move_alloc (ct_map%ptr%term, ct_terms)
      allocate (ct_map%ptr%term(ix))
      ct_map%ptr%term(1:size(ct_terms)) = ct_terms
      deallocate (ct_terms)
    endif
  endif

  ! 1) chop "=", 2) chop to "{", 3) chop to "}", 4) chop to "," if it exists

  call get_next_word (word, ix_word1, ':,={}', delim1, delim_found, .true.) 
  call get_next_word (word, ix_word2, ':,={}', delim2, delim_found, .true., call_check = .true.)  

  if (delim1 /= '=' .or. delim2 /= '{' .or. ix_word1 /= 0 .or. ix_word2 /= 0) then
    call parser_error ('CONFUSED SYNTAX FOR TERM IN WIGGLER: ' // ele%name, str_ix)
    return
  endif

  err_str = trim(ele%name) // ' ' // str_ix
  ct_term => ct_map%ptr%term(ix)

  call evaluate_value (err_str, ct_term%coef, lat, delim, delim_found, err_flag, ',');   if (err_flag) return
  call evaluate_value (err_str, ct_term%kx, lat, delim, delim_found, err_flag, ',');     if (err_flag) return
  call evaluate_value (err_str, ct_term%ky, lat, delim, delim_found, err_flag, ',');     if (err_flag) return
  call evaluate_value (err_str, ct_term%kz, lat, delim, delim_found, err_flag, ',');     if (err_flag) return
  call evaluate_value (err_str, ct_term%phi_z, lat, delim, delim_found, err_flag, ',}'); if (err_flag) return

  old_style_input = .true.
  family = y_family$

  if (delim == ',') then
    ct_term%x0 = ct_term%phi_z
    call evaluate_value (err_str, ct_term%y0, lat, delim, delim_found, err_flag, ','); if (err_flag) return
    call evaluate_value (err_str, ct_term%phi_z, lat, delim, delim_found, err_flag, ','); if (err_flag) return
    call get_switch ('FAMILY', ['X ', 'Y ', 'QU', 'SQ'], family, err_flag, ele, delim, delim_found); if (err_flag) return
    if (.not. expect_this ('}', .true., .false., 'AFTER "FAMILY" SWITCH', ele, delim, delim_found)) return
    old_style_input = .false.
    call parser_error ('"HYBRID" STYLE WIGGLER TERMS DEPRECATED. PLEASE CONVERT TO CARTESIAN_MAP FORM.', level = s_warn$)
  endif

  kx = ct_term%kx
  ky = ct_term%ky
  kz = ct_term%kz
  tol = 1d-5 * (kx**2 + ky**2 + kz**2)

  if (abs(ky**2 - kx**2 - kz**2) < tol) then
    select case (family)
    case (x_family$);   ct_term%type = hyper_y_family_x$
    case (y_family$);   ct_term%type = hyper_y_family_y$
    case (qu_family$);  ct_term%type = hyper_y_family_qu$
    case (sq_family$);  ct_term%type = hyper_y_family_sq$
    end select

    if (old_style_input) then
      if (ct_term%kx == 0) ct_term%kx = 1d-30  ! Something small to prevent divide by zero problems.
    endif

  elseif (abs(ky**2 + kx**2 - kz**2) < tol) then
    select case (family)
    case (x_family$);   ct_term%type = hyper_xy_family_x$
    case (y_family$);   ct_term%type = hyper_xy_family_y$
    case (qu_family$);  ct_term%type = hyper_xy_family_qu$
    case (sq_family$);  ct_term%type = hyper_xy_family_sq$
    end select

    if (old_style_input) then
      ct_term%coef = ct_term%coef * ct_term%kz / ct_term%ky
      if (ct_term%kx == 0) ct_term%kx = 1d-30  ! Something small to prevent divide by zero problems.
      if (ct_term%ky == 0) ct_term%ky = 1d-30  ! Something small to prevent divide by zero problems.
    endif

  elseif (abs(ky**2 - kx**2 + kz**2) < tol) then
    select case (family)
    case (x_family$);   ct_term%type = hyper_x_family_x$
    case (y_family$);   ct_term%type = hyper_x_family_y$
    case (qu_family$);  ct_term%type = hyper_x_family_qu$
    case (sq_family$);  ct_term%type = hyper_x_family_sq$
    end select

    if (old_style_input) then
      ct_term%coef = ct_term%coef * ct_term%kx / ct_term%ky
      if (ct_term%ky == 0) ct_term%ky = 1d-30  ! Something small to prevent divide by zero problems.
    endif

  else
    call parser_error ('WIGGLER TERM DOES NOT HAVE CONSISTANT Kx, Ky, and Kz', &
                  'FOR WIGGLER: ' // ele%name // '  ' // str_ix)
    err_flag = .true.
    return
  endif

  call get_next_word (word, ix_word,  ':,=()', delim,  delim_found, .true.)  
  if (ix_word /= 0) then
    call parser_error ('BAD SYNTAX FOR WIGGLER: ' // ele%name, str_ix)
    err_flag = .true.
    return
  endif

  ele%sub_key = map_type$
  ele%field_calc = fieldmap$
  return

endif

! Check that next delim is a "=". 
! If not, it might be a flag attribute or an attribute that has a default value.

if (delim /= '=')  then
  err_flag = .false.

  if (ele%key == multipole$ .and. ix_attrib >= t0$) then
    if (.not. associated(ele%a_pole)) call multipole_init (ele, magnetic$)
    ele%b_pole(ix_attrib-t0$) = pi / (2*(ix_attrib-t0$) + 2)
    return
  endif

  if (attrib_word == 'TILT') then
    select case (ele%key)
    case (quadrupole$, sol_quad$) 
      ele%value(tilt$) = pi / 4
      return
    case (sextupole$) 
      ele%value(tilt$) = pi / 6
      return
    case (octupole$) 
      ele%value(tilt$) = pi / 8
      return
    case default
      call parser_error ('SORRY I''M NOT PROGRAMMED TO USE A "TILT" DEFAULT' // &
                         'FOR A: ' // key_name(ele%key), 'FOR: ' // ele%name)
      err_flag = .true.
      return
    end select
  endif

  if (ele%key == sbend$ .or. ele%key == rbend$) then
    select case (ix_attrib)
    case (fint$)
      ele%value(fint$) = 0.5
      return
    case (fintx$)
      ele%value(fintx$) = 0.5
      return
    end select
  endif

  select case (attrib_word)

  case ('SUPERIMPOSE')
    ele%lord_status = super_lord$
    pele%superposition_has_been_set = .true.

  case default
    call parser_error ('EXPECTING "=" AFTER ATTRIBUTE: ' // word,  'FOR ELEMENT: ' // ele%name)
    err_flag = .true.
  end select

  return
endif

! lr_wake_spline

if (attrib_word == 'LR_WAKE_SPLINE') then
  if (associated(ele%wake)) then
    n = size(ele%wake%lr_spline)
    call move_alloc(ele%wake%lr_spline, lr_pa_temp)
    allocate (ele%wake%lr_spline(n+1))
    do i = 1, n
      allocate (ele%wake%lr_spline(i)%spline(0))
      allocate (ele%wake%lr_spline(i)%bunch(0))
    enddo
    ele%wake%lr_spline(1:n) = lr_pa_temp
    deallocate(lr_pa_temp)
    lr_pa => ele%wake%lr_spline(n+1)
  else
    call init_wake (ele%wake, 0, 0, 0, 1)
    lr_pa => ele%wake%lr_spline(1)
  endif

  if (.not. expect_this ('{', .false., .true., 'AFTER "LR_WAKE_SPLINE"', ele, delim, delim_found)) return

  call parse_wake_lr_spline(lr_pa, ele, lat, delim, delim_found, err_flag)
  return
endif

! get the value of the attribute.
! Stuff like TYPE, ALIAS, and DESCRIP attributes are special because their "values"
! are character strings

if (attrib_word(1:16) == 'CUSTOM_ATTRIBUTE') then
  call bmad_parser_type_get (ele, attrib_word, delim, delim_found, str)
  call set_attribute_alias(attrib_word, str, err_flag, lat)
  if (err_flag) call parser_error ('CANNOT SET PARAMETER[' // trim(attrib_word) // ']')
  return
endif

select case (attrib_word)

case ('REFERENCE')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  call get_next_word(pele%ref_name, ix_word,  '=,', delim, delim_found, .true.)

case ('OFFSET')
  call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag)
  if (err_flag) return
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  pele%offset = value

case ('FIELD_OVERLAPS')

  ! If pele is not present then bmad_parser2 is the parser and this is an element in the lattice.
  ! In this case, simple call create_field_overlap directly.

  call get_list_of_names (ele, 'FIELD_OVERLAPS', name_list, delim, delim_found, err_flag)
  nn = size(name_list)

  if (present(pele)) then
    n = 0
    if (allocated(pele%field_overlaps)) n = size(pele%field_overlaps)
    call re_allocate (pele%field_overlaps, n+nn)
    pele%field_overlaps(n+1:n+nn) = name_list

  else
    do i = 1, n
      call create_field_overlap (ele%branch%lat, ele%name, name_list(i), err_flag)
      if (err_flag) then
        call parser_error ('OVERLAP ELEMENT: ' // name_list(i), 'NOT FOUND FOR OVERLAPPING ELEMENT: ' // ele%name)
      endif
    enddo
  endif

case('TYPE', 'ALIAS', 'DESCRIP', 'SR_WAKE_FILE', 'LR_WAKE_FILE', 'LATTICE', 'TO', &
     'TO_LINE', 'TO_ELEMENT', 'CRYSTAL_TYPE', 'MATERIAL_TYPE', 'ORIGIN_ELE', 'PHYSICAL_SOURCE')
  call bmad_parser_type_get (ele, attrib_word, delim, delim_found, pele = pele)

case ('REF_ORBIT')
  if (.not. parse_real_list (lat, ele%name // ' REF_ORBIT', ele%taylor%ref, .true.)) return
  if (.not. expect_one_of (', ', .false., ele, delim, delim_found)) return

case ('PTC_MAX_FRINGE_ORDER')
  call parser_error ('PLEASE CONVERT "PARAMETER[PTC_MAX_FRINGE_ORDER]" TO "BMAD_COM[PTC_MAX_FRINGE_ORDER]"', level = s_warn$)
  call parser_get_integer (bmad_com%ptc_max_fringe_order, word, ix_word, delim, delim_found, err_flag)
  bp_com%extra%ptc_max_fringe_order_set = .true.

case ('TAYLOR_ORDER')
  call parser_get_integer (ix, word, ix_word, delim, delim_found, err_flag)
  if (ix <= 0) then
    call parser_error ('TAYLOR_ORDER IS LESS THAN 1')
    return
  endif
  ptc_com%taylor_order_saved = ix
  lat%input_taylor_order = ix

case ('RUNGE_KUTTA_ORDER')
  call parser_get_integer (ix, word, ix_word, delim, delim_found, err_flag)
  if (ix /= 2 .and. ix /= 4) then
    call parser_error ('RUNGE_KUTTA_ORDER NOT EQUAL TO 2 OR 4')
    return
  endif
  bmad_com%runge_kutta_order = ix
  bp_com%extra%runge_kutta_order_set = .true.

case ('SYMPLECTIFY') 
  if (how == def$ .and. (delim == ',' .or. .not. delim_found)) then
    ele%symplectify = .true.
  else
    call get_logical (attrib_word, ele%symplectify, err_flag)
  endif
  
case ('IS_ON')
  call get_logical (attrib_word, ele%is_on, err_flag)

case ('USE_HARD_EDGE_DRIFTS')
  call parser_error ('PLEASE CONVERT "PARAMETER[USE_HARD_EDGE_DRIFTS]" TO "BMAD_COM[USE_HARD_EDGE_DRIFTS]"', level = s_warn$)
  call get_logical (attrib_word, bmad_com%use_hard_edge_drifts, err_flag)
  bp_com%extra%use_hard_edge_drifts_set = .true.

case ('SUPERIMPOSE')  ! ele[superimpose] = False case
  call get_logical (attrib_word, logic, err_flag)
  if (logic) then
    ele%lord_status = super_lord$
  else
    ele%lord_status = not_a_lord$
  endif
  pele%superposition_has_been_set = .true.

case ('APERTURE_AT')
  call get_switch (attrib_word, aperture_at_name(1:), ele%aperture_at, err_flag, ele, delim, delim_found)

case ('APERTURE_TYPE')
  call get_switch (attrib_word, aperture_type_name(1:), ele%aperture_type, err_flag, ele, delim, delim_found)

case ('ABSOLUTE_TIME_TRACKING')
  call get_logical (attrib_word, lat%absolute_time_tracking, err_flag)

case ('CAVITY_TYPE')
  call get_switch (attrib_word, cavity_type_name(1:), ix, err_flag, ele, delim, delim_found)
  ele%value(cavity_type$) = ix

case ('COUPLER_AT')
  call get_switch (attrib_word, end_at_name(1:), ix, err_flag, ele, delim, delim_found)
  ele%value(coupler_at$) = ix

case ('CREATE_JUMBO_SLAVE')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  call get_logical (attrib_word, pele%create_jumbo_slave, err_flag)

case ('CSR_CALC_ON')
  call get_logical (attrib_word, ele%csr_calc_on, err_flag)

case ('DEFAULT_TRACKING_SPECIES')
  call get_next_word (word, ix_word, ':,=(){}', delim, delim_found, .false.)
  ix = species_id(word)
  if (ix == invalid$) then
    call parser_error ('INVALID PARTICLE SPECIES: ' // word)
    return
  endif
  ele%value(default_tracking_species$) = ix

case ('ENERGY_DISTRIBUTION')
  call get_switch (attrib_word, distribution_name(1:), ix, err_flag, ele, delim, delim_found)
  ele%value(energy_distribution$) = ix

case ('EXACT_MULTIPOLES')
  call get_switch (attrib_word, exact_multipoles_name(1:), ix, err_flag, ele, delim, delim_found)
  ele%value(exact_multipoles$) = ix

case ('PTC_EXACT_MODEL')
  call get_logical (attrib_word, logic, err_flag)
  if (.not. err_flag) call set_ptc (exact_modeling = logic)

case ('PTC_EXACT_MISALIGN')
  call get_logical (attrib_word, logic, err_flag)
  if (err_flag) return
  call set_ptc (exact_misalign = logic)

case ('OFFSET_MOVES_APERTURE')
  call get_logical (attrib_word, ele%offset_moves_aperture, err_flag)

case ('FIELD_MASTER', 'HARMON_MASTER')
  call get_logical (attrib_word, ele%field_master, err_flag)

case ('SCALE_MULTIPOLES')
  call get_logical (attrib_word, ele%scale_multipoles, err_flag)

case ('FIELD_CALC')
  call get_switch (attrib_word, field_calc_name(1:), ele%field_calc, err_flag, ele, delim, delim_found)

case ('REF_ORIGIN')
  call get_switch (attrib_word, anchor_pt_name(1:), pele%ref_pt, err_flag, ele, delim, delim_found)

case ('ELE_ORIGIN')
  call get_switch (attrib_word, anchor_pt_name(1:), pele%ele_pt, err_flag, ele, delim, delim_found)

case ('PTC_FRINGE_GEOMETRY')
  call get_switch (attrib_word, ptc_fringe_geometry_name(1:), ix, err_flag, ele, delim, delim_found)
  ele%value(ptc_fringe_geometry$) = ix

case ('FRINGE_TYPE')
  if (ele%key == rbend$ .or. ele%key == sbend$) then
    call get_switch (attrib_word, fringe_type_name(1:), ix, err_flag, ele, delim, delim_found)
  else
    call get_switch (attrib_word, fringe_type_name(1:n_non_bend_fringe_type$), ix, err_flag, ele, delim, delim_found)
  endif
  ele%value(fringe_type$) = ix

case ('HIGHER_ORDER_FRINGE_TYPE')
  call get_switch (attrib_word, higher_order_fringe_type_name(1:), ix, err_flag, ele, delim, delim_found)
  ele%value(higher_order_fringe_type$) = ix

case ('MAT6_CALC_METHOD')
  call get_switch (attrib_word, mat6_calc_method_name(1:), switch, err_flag, ele, delim, delim_found)
  if (err_flag) return
  if (.not. valid_mat6_calc_method (ele, not_set$, switch)) then
    if (wild_key0) then
      err_flag = .false.
    else
      err_flag = .true.
      call parser_error ('NOT A VALID MAT6_CALC_METHOD: ' // mat6_calc_method_name(switch), &
                         'FOR: ' // trim(ele%name), 'WHICH IS A: ' // key_name(ele%key))
    endif
    return
  endif
  ele%mat6_calc_method = switch

case ('REF_COORDINATES')
  call get_switch (attrib_word, end_at_name(1:2), ix, err_flag, ele, delim, delim_found)
  ele%value(ref_coordinates$) = ix

case ('REF_ORBIT_FOLLOWS')
  call get_switch (attrib_word, ref_orbit_follows_name(1:), ix, err_flag, ele, delim, delim_found)
  ele%value(ref_orbit_follows$) = ix

case ('MODE')
  call get_switch (attrib_word, mode_name(1:), ix, err_flag, ele, delim, delim_found)
  ele%value(mode$) = ix

case ('PTC_INTEGRATION_TYPE')
  call get_switch (attrib_word, ptc_integration_type_name(1:), ele%ptc_integration_type, err_flag, ele, delim, delim_found)

case ('PARTICLE')
  call get_next_word (word, ix_word, ':,=(){}', delim, delim_found, .false.)
  ix = species_id(word)
  if (ix == invalid$) then
    call parser_error ('INVALID PARTICLE SPECIES: ' // word)
    return
  endif
  branch => pointer_to_branch(ele%name, lat, .true.)
  if (associated(branch)) branch%param%particle = ix 
  ele%value(particle$) = ix

case ('PTC_FIELD_GEOMETRY')
  call get_switch (attrib_word, ptc_field_geometry_name(1:), ix, err_flag, ele, delim, delim_found)
  ele%value(ptc_field_geometry$) = ix

  if (ele%key == sbend$ .and. ix == true_rbend$) then
    call parser_error ('TRUE_RBEND IS NOT A VALID PTC_FIELD_GEOMETRY VALUE FOR AN SBEND')
    return
  endif

case ('GEOMETRY')
  call get_switch (attrib_word, geometry_name(1:), ix, err_flag, ele, delim, delim_found)
  branch => pointer_to_branch(ele%name, lat, .true.)
  if (associated(branch)) branch%param%geometry = ix
  ele%value(geometry$) = ix

case ('LIVE_BRANCH')
  call get_logical_real (attrib_word, ele%value(live_branch$), err_flag)
  if (err_flag) return
  branch => pointer_to_branch(ele%name, lat, .true.)
  if (associated(branch)) branch%param%live_branch = is_true(ele%value(live_branch$))

case ('PHOTON_TYPE')
  call get_switch (attrib_word, photon_type_name(1:), ix, err_flag, ele, delim, delim_found)
  lat%photon_type = ix   ! photon_type has been set.

case ('LATTICE_TYPE')   ! Old style
  call parser_error ('PARAMETER[LATTICE_TYPE] IS OLD SYNTAX.', &
                     'PLEASE REPLACE WITH PARAMETER[GEOMETRY] = OPEN/CLOSED', &
                     'THIS PROGRAM WILL RUN NORMALLY...', level = s_warn$)
  call get_switch (attrib_word, lattice_type_name(1:), ix, err_flag, ele, delim, delim_found)
  ele%value(geometry$) = ix

case ('FRINGE_AT')
  call get_switch (attrib_word, end_at_name(1:), ix, err_flag, ele, delim, delim_found)
  ele%value(fringe_at$) = ix

case ('ORIGIN_ELE_REF_PT')
  call get_switch (attrib_word, ref_pt_name(1:), ix, err_flag, ele, delim, delim_found)
  ele%value(origin_ele_ref_pt$) = ix

case ('SPATIAL_DISTRIBUTION')
  call get_switch (attrib_word, distribution_name(1:), ix, err_flag, ele, delim, delim_found)
  ele%value(spatial_distribution$) = ix

case ('SPIN_TRACKING_METHOD')
  if (attrib_word == 'BMAD_STANDARD') then
    call parser_error ('SPIN_TRACKING_METHOD = BMAD_STANDARD NOW NO LONGER VALID.', &
                     'PLEASE REPLACE WITH SPIN_TRACKING_METHOD = TRACKING.', &
                     'THIS PROGRAM WILL RUN NORMALLY...', level = s_warn$)
    attrib_word = 'TRACKING'
  endif
  call get_switch (attrib_word, spin_tracking_method_name(1:), switch, err_flag, ele, delim, delim_found)
  if (err_flag) return
  if (.not. valid_spin_tracking_method (ele, switch)) then
    if (wild_key0) then
      err_flag = .false.
    else
      call parser_error ('NOT A VALID SPIN_TRACKING_METHOD: ' // word, &
                         'FOR: ' // trim(ele%name), 'WHICH IS A: ' // key_name(ele%key))
    endif
    return
  endif
  ele%spin_tracking_method = switch

case ('TAYLOR_MAP_INCLUDES_OFFSETS')
  call get_logical (attrib_word, ele%taylor_map_includes_offsets, err_flag)

case ('TRACKING_METHOD')
  call get_switch (attrib_word, tracking_method_name(1:), switch, err_flag, ele, delim, delim_found)
  if (err_flag) return
  if (.not. valid_tracking_method (ele, not_set$, switch)) then
    if (wild_key0) then
      err_flag = .false.
    else
      call parser_error ('NOT A VALID TRACKING_METHOD: ' // word, &
                         'FOR: ' // trim(ele%name), 'WHICH IS A: ' // key_name(ele%key))
    endif
    return
  endif
  ele%tracking_method = switch

case ('VELOCITY_DISTRIBUTION')
  call get_switch (attrib_word, distribution_name(1:), ix, err_flag, ele, delim, delim_found)
  ele%value(velocity_distribution$) = ix

case default   ! normal attribute

  call pointer_to_attribute (ele, attrib_word, .true., a_ptr, err_flag2, .false.)

  if (attribute_type(attrib_word) == is_logical$) then
    if (associated (a_ptr%l)) then
      call get_logical (trim(ele%name) // ' ' // attrib_word, a_ptr%l, err_flag)
    else
      call get_logical_real (attrib_word, ele%value(ix_attrib), err_flag)
    endif
    if (err_flag) return

  else
    call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag)
    if (err_flag) return

    ! multipole attribute?
    if (ele%key == hybrid$ .and. is_attribute(ix_attrib, multipole$)) then
      ele%vec0(ix_attrib-a0$) = value
    elseif (ele%key == hybrid$ .and. is_attribute(ix_attrib, elec_multipole$)) then
      i = 1 + (ix_attrib - a0_elec$ - 1) / 6
      j = ix_attrib - a0_elec$ - 6 * (i - 1)
      ele%mat6(i,j) = value
    elseif (is_attribute(ix_attrib, multipole$) .and. attrib_word(1:4) /= 'CURV') then  
        if (.not. associated(ele%a_pole)) call multipole_init (ele, magnetic$)
        if (ix_attrib >= b0$) then
          ele%b_pole(ix_attrib-b0$) = value
        else
          ele%a_pole(ix_attrib-a0$) = value
        endif
    ! Electric multipole attribute
    elseif (is_attribute(ix_attrib, elec_multipole$)) then
        if (.not. associated(ele%a_pole_elec)) call multipole_init (ele, electric$)
        if (ix_attrib >= b0_elec$) then
          ele%b_pole_elec(ix_attrib-b0_elec$) = value
        else
          ele%a_pole_elec(ix_attrib-a0_elec$) = value
        endif
    !
    elseif (attrib_word == 'RAN_SEED') then
      call ran_seed_put (nint(value))  ! init random number generator
      if (nint(value) == 0) then  ! Using system clock -> Not determinisitc.
        bp_com%extra%deterministic = 0
      else
        bp_com%extra%deterministic = 2
      endif
    elseif (attrib_word == 'APERTURE') then
      ele%value(x1_limit$) = value
      ele%value(x2_limit$) = value
      ele%value(y1_limit$) = value
      ele%value(y2_limit$) = value
    elseif (attrib_word == 'X_LIMIT') then
      ele%value(x1_limit$) = value
      ele%value(x2_limit$) = value
    elseif (attrib_word == 'Y_LIMIT') then
      ele%value(y1_limit$) = value
      ele%value(y2_limit$) = value
    elseif (ix_attrib > num_ele_attrib$) then
      ! TT* Taylor Terms need do_allocation = True.
      call pointer_to_attribute (ele, attrib_word, .true., a_ptr, err_flag, .false.)
      if (err_flag .or. .not. associated(a_ptr%r)) then
        call parser_error ('BAD ATTRIBUTE: ' // attrib_word, 'FOR ELEMENT: ' // ele%name)
        return
      endif
      a_ptr%r = value
      if (attrib_word == 'X_POSITION' .or. attrib_word == 'X_POSITION' .or. &
          attrib_word == 'X_POSITION' .or. attrib_word == 'THETA_POSITION' .or. &
          attrib_word == 'PHI_POSITION' .or. attrib_word == 'PSI_POSITION') ele%value(floor_set$) = 1
    else
      ele%value(ix_attrib) = value
      ix = len_trim(attrib_word)
      if (ix > 9 .and. index(attrib_word, '_GRADIENT') == ix-8) ele%field_master = .true.
      if (ix > 6 .and. index(attrib_word, '_FIELD') == ix-5) ele%field_master = .true.
      if (ix > 10 .and. index(attrib_word, '_FIELD_ERR') == ix-9) ele%field_master = .true.
      if (attrib_word(1:3) == 'BL_') ele%field_master = .true.
      if (ele%key == elseparator$ .and. attrib_word == 'VOLTAGE') ele%field_master = .true.
      if (ele%key == elseparator$ .and. attrib_word == 'E_FIELD') ele%field_master = .true.

      select case (attrib_word)

      case ('E_TOT')
        if (ele%key == def_parameter$) then
          lat%ele(0)%value(e_tot$) = value
          lat%ele(0)%value(p0c$) = -1
        else
          ele%value(p0c$) = -1
        endif

        branch => pointer_to_branch(ele%name, lat, .true.)
        if (associated(branch)) then
          branch%ele(0)%value(e_tot$) = value
          call set_flags_for_changed_attribute (branch%ele(0), branch%ele(0)%value(e_tot$))
        endif

      case ('ENERGY')    ! Only in def_mad_beam
        lat%ele(0)%value(e_tot$) = 1d9 * value
        lat%ele(0)%value(p0c$) = -1

      case ('P0C')
        if (ele%key == def_parameter$) then
          lat%ele(0)%value(p0c$) = value
          lat%ele(0)%value(e_tot$) = -1
        else
          ele%value(e_tot$) = -1
        endif

        branch => pointer_to_branch(ele%name, lat, .true.)
        if (associated(branch)) then
          branch%ele(0)%value(p0c$) = value
          call set_flags_for_changed_attribute (branch%ele(0), branch%ele(0)%value(p0c$))
        endif

      case ('PC')    ! Only in def_mad_beam
        lat%ele(0)%value(p0c$) = 1d9 * value
        ele%value(e_tot$) = -1

      case ('LR_FREQ_SPREAD')
        call randomize_lr_wake_frequencies (ele, set_done)
        if (set_done) call bp_set_ran_status

      case ('N_PART')
        branch => pointer_to_branch(ele%name, lat, .true.)
        if (associated(branch)) branch%param%n_part = value

      end select

    endif

  endif

end select

err_flag = .false.

!--------------------------------------------------------
contains

subroutine parser_find_ele_for_attrib_transfer (attribute)

character(*) attribute

call get_next_word (word2, ix_word, '[],(){}', delim2, delim_found, call_check = .true.)
if (delim2 /= ']' .or. word2 /= attribute) then
  call parser_error ('BAD ' // attribute // ' CONSTRUCT')
  return
endif
if (.not. expect_this (' ', .false., .false., '', ele, delim, delim_found)) return
call lat_ele_locator (word, lat, eles, n, err_flag)
if (err_flag .or. n /= 1) then
  call parser_error ('LATTICE ELEMENT NOT FOUND: ' // word)
  return
endif



end subroutine parser_find_ele_for_attrib_transfer

!--------------------------------------------------------
! contains

function attrib_free_problem (attrib_name) result (is_problem)

type (ele_attribute_struct) attrib_info
character(*) attrib_name
logical is_problem, is_free

! Attributes may be definitely free, definitely dependent, or may be free or
! dependent depending upon the state of other element parameters.

! If not check_free then at least check if it is a dependent attribute.

is_problem = .false.

if (logic_option(.false., check_free)) then
  is_free = attribute_free (ele, attrib_name, bp_com%print_err)
  if (.not. is_free) then
    call parser_error ('ATTRIBUTE NOT FREE TO BE SET: ' // attrib_name, &
                                      'FOR: ' // ele%name)
    err_flag = .true.
    is_problem = .true.
  endif
else
  attrib_info = attribute_info(ele, attribute_index(ele, attrib_name))
  if (attrib_info%type == dependent$) then
    call parser_error ('DEPENDENT ATTRIBUTE NOT FREE TO BE SET: ' // attrib_name, &
                                      'FOR: ' // ele%name)
    is_problem = .true.
  endif
endif

end function attrib_free_problem

!--------------------------------------------------------
! contains

subroutine get_logical (attrib_name, this_logic, err)

character(*) attrib_name
logical this_logic, err

!

call get_next_word (word, ix_word, ':,=()', delim, delim_found, .true.)
this_logic = evaluate_logical (word, ios)
if (ios /= 0 .or. ix_word == 0) then
  call parser_error ('BAD "' // trim(attrib_name) // '" SWITCH FOR: ' // ele%name, 'I DO NOT UNDERSTAND: ' // word)
  err = .true.
else
  err = .false.
endif

end subroutine get_logical

!--------------------------------------------------------
! contains

subroutine get_logical_real (name, logic_real, err)

character(*) name
real(rp) logic_real
logical this_logical, err

!

call get_logical (name, this_logical, err)
if (err) return

if (this_logical) then
  logic_real = 1
else
  logic_real = 0
endif

err = .false.

end subroutine get_logical_real

end subroutine parser_set_attribute 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

function parser_translate_attribute_name (key, word) result (trans_word)

implicit none

character(*) word
character(40) trans_word
integer key

!

trans_word = word

if (word == 'REF') then
  trans_word = 'REFERENCE' ! allowed abbrev

elseif (key == rbend$) then
  if (word == 'L') then
    trans_word = 'L_CHORD'
  elseif (word == 'L_ARC') then
    trans_word = 'L'
  endif

elseif (key == rfcavity$ .and. word == 'LAG') then   ! For MAD compatibility
  trans_word = 'PHI0'
endif

end function parser_translate_attribute_name 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine get_called_file (delim, call_file, xsif_called, err)

implicit none

character(1) delim
character(*) call_file

integer ix_word, ix, n
logical delim_found, finished, xsif_called, err

!

err = .true.

if (delim /= ',')  call parser_error ('"CALL" NOT FOLLOWED BY COMMA', stop_here = .true.)
call get_next_word(call_file, ix_word, ':=,', delim, delim_found, .true.)

if (ix_word == 0) then
  call parser_error ('NOTHING AFTER "CALL"', stop_here = .true.)
  return
elseif (index('FILENAME', call_file(:ix_word)) /= 1) then
  call parser_error ('INVALID "CALL" COMMAND', stop_here = .true.)
  return
elseif (delim /= '=') then
  call parser_error ('NO "=" AFTER "FILENAME"', stop_here = .true.)
  return
endif

call get_next_word(call_file, ix_word, ',', delim, delim_found, .false.)
if (ix_word == 0) then
  call parser_error ('NO FILE NAME SPECIFIED', stop_here = .true.)
  return
endif

if (call_file(1:1) == '"') then
  call_file = call_file(2:)
  ix = index(call_file, '"')
  if (ix == 0 .or. ix /= len_trim(call_file)) then
    call parser_error ('MISSING DOUBLE QUOTE MARK (") FOR CALL STATEMENT', stop_here = .true.)
    return
  endif
  call_file(ix:ix) = ' '
endif

if (call_file(1:1) == "'") then
  call_file = call_file(2:)
  ix = index(call_file, "'")
  if (ix == 0 .or. ix /= len_trim(call_file)) then
    call parser_error ("MISSING SINGLE QUOTE MARK (') FOR CALL STATEMENT", stop_here = .true.)
    return
  endif
  call_file(ix:ix) = ' '
endif

if (call_file(1:6) == 'xsif::') then
  call_file = call_file(7:)
  n = size(bp_com%lat_file_names)
  if (n < bp_com%num_lat_files + 1) &
              call re_allocate (bp_com%lat_file_names, n + 100)
  bp_com%num_lat_files = bp_com%num_lat_files + 1 
  inquire (file = call_file, name = bp_com%lat_file_names(bp_com%num_lat_files))
  xsif_called = .true.
  err = .false.
  return
endif

xsif_called = .false.
call parser_file_stack ('push', call_file, finished, err) ! err gets set here

end subroutine get_called_file

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine add_this_taylor_term (ele, i_out, coef, expn) 
!
! Subroutine used by bmad_parser and bmad_parser2 to parse the input file.
! This subroutine is not intended for general use.
!-

subroutine add_this_taylor_term (ele, i_out, coef, expn)

implicit none

type (ele_struct), target :: ele
type (taylor_struct), pointer :: taylor

real(rp) coef
integer i, j, i_out, expn(6)

!

if (i_out > 100) then
  i_out = i_out - 101
  i = i_out / 3
  j = i_out - 3 * i
  taylor => ele%spin_taylor(i+1,j+1) 
else
  if (i_out < 1 .or. i_out > 6) then
    call parser_error ('"OUT" VALUE IN TAYLOR TERM NOT IN RANGE (1 - 6)', &
                  'FOR TAYLOR ELEMENT: ' // ele%name)
    return
  endif
  taylor => ele%taylor(i_out)
endif

call add_taylor_term (taylor, coef, expn, .true.)

end subroutine add_this_taylor_term

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_next_word (word, ix_word, delim_list, delim, delim_found, upper_case_word, call_check)
!
! Subroutine to get the next word from the input stream.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
! Input:
!   word            -- Character(*): Word returned
!   delim_list      -- Character(*): List of valid delimiters
!   upper_case_word -- Logical, optional: if True then convert word to 
!                       upper case. Default is True.
!   call_check      -- Logical, optional: If present and True then check for 'call::<filename>' construct.
!
! Output
!   ix_word     -- Integer: length of word argument
!   delim       -- Character(1): Actual delimiter found
!   delim_found -- Logical: Set true if a delimiter found. A delimiter
!                    may not be found if the end of the line is reached first.
!-

subroutine get_next_word (word, ix_word, delim_list, delim, delim_found, upper_case_word, call_check)

implicit none

integer ix_a, ix_word

character(*) word, delim_list, delim

integer n, ix

logical delim_found, end_of_file
logical, optional :: upper_case_word, call_check

character(n_parse_line) line
character(6) str

! Possible inline call...

if (logic_option(.false., call_check)) then
  call string_trim(bp_com%parse_line, bp_com%parse_line, ix)
  call str_upcase (str, bp_com%parse_line(1:6))
  if (str == 'CALL::') then
    bp_com%parse_line = bp_com%parse_line(7:)
    call word_read (bp_com%parse_line, ',} ',  line, ix_word, delim, delim_found, bp_com%parse_line)
    bp_com%parse_line = delim // bp_com%parse_line  ! put delim back on parse line.
    call parser_file_stack ('push_inline', line); if (bp_com%fatal_error_flag) return
  endif
endif

! Check for continuation character and, if found, then load more characters
! into the parse line from the lattice file. 
! If the input is not from a file then skip this.

if (bp_com%input_from_file) then 
  do
    n = len_trim(bp_com%parse_line)
    if (n == 0 .or. n > 90) exit

    select case (bp_com%parse_line(n:n))
    case (',', '(', '{', '[', '=')
      call load_parse_line('continue', n+2, end_of_file)
      if (end_of_file) exit

    case ('&')
      call load_parse_line('continue', n, end_of_file)
      if (end_of_file) exit

    case default
      ! If in an inline called file then make sure the rest of the file is blank and
      ! return to the calling file
      if (bp_com%inline_call_active) then
        call load_parse_line('continue', n+2, end_of_file)
        if (bp_com%parse_line(n+1:) /= '') then
          call string_trim (bp_com%parse_line(n+1:), line, ix)
          call str_upcase (line(1:10), line(1:10))
          if (line /= 'END_FILE') THEN
            call parser_error ('EXTRA STUFF IN FILE')
          endif
        endif
        bp_com%parse_line(n+1:) = ''
        call parser_file_stack ('pop')
      else
        exit
      endif

    end select

  enddo
endif

! Get the first word in bp_com%parse_line

call word_read (bp_com%parse_line, delim_list, word, ix_word, delim, delim_found, bp_com%parse_line)

if (len(word) < ix_word) then
  call parser_error ('BAD WORD: ' // bp_com%parse_line)
  ix_word = len(word)
endif

if (present(upper_case_word)) then
  if (upper_case_word) call str_upcase (word, word)
else
  call str_upcase (word, word)
endif

! Note: "var := num" is old-style variable definition syntax.
! If delim is ":" and next char is "=" then use "=" as the delim

if (delim == ':' .and. index(delim_list, '=') /= 0 .and. bp_com%parse_line(1:1) == '=') then
  delim = '='
  bp_com%parse_line = bp_com%parse_line(2:)
endif

end subroutine get_next_word

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_file_stack (how, file_name_in, finished, err)
!
! Subroutine to keep track of the files that are opened for reading.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_file_stack (how, file_name_in, finished, err)

implicit none

integer, parameter :: f_maxx = 20
type (stack_file_struct), save, target :: file(0:f_maxx)

integer, save :: i_level
integer i, ix, ios, n, n_file

character(*) how
character(*), optional :: file_name_in
character(200) file_name, basename, file_name2

logical, optional :: finished, err
logical found_it, is_relative, valid, err_flag

! "Init" means init

if (present(err)) err = .true.

if (how == 'init') then
  i_level = 0
  call fullfilename ('./', file(0)%dir)
  if (present(err)) err = .false.
  if (.not. allocated(bp_com%lat_file_names)) allocate(bp_com%lat_file_names(100))
  bp_com%inline_call_active = .false.
  return
endif

! "push" means open a file and put its name on the stack.

if (present(finished)) finished = .false.

select case (how)
case ('push', 'push_inline')

  i_level = i_level + 1    ! number of files currently open
  if (i_level > f_maxx) then
    call parser_error ('CALL NESTING GREATER THAN 20 LEVELS')
    if (global_com%exit_on_error) call err_exit
  endif

  if (how == 'push_inline') then
    file(i_level)%parse_line_saved = bp_com%parse_line
    file(i_level)%inline_call_active = .true.    
    bp_com%parse_line = '&'
    bp_com%inline_call_active = .true.
  endif

  bp_com%current_file => file(i_level)
  bp_com%calling_file => file(i_level-1)

  if (i_level == 1) then   ! if we are just starting out then init some vars.
    bp_com%num_lat_files = 0           ! total number of files opened
    bp_com%error_flag = .false.  ! set to true on an error
    bp_com%current_file%full_name = ' '
    bp_com%input_line_meaningful = .false.
    call init_bmad_parser_common
  endif

  call fullfilename (file_name_in, file_name2, valid)
  if (.not. valid) then
    call parser_error ('MALFORMED FILE NAME: ' // file_name_in, stop_here = .true.)
    if (global_com%exit_on_error) call err_exit
    do i = 1, i_level-1
      close (file(i_level)%f_unit)
    enddo
    return
  endif

  ix = splitfilename (file_name2, file(i_level)%dir, basename, is_relative)

  if (bp_com%use_local_lat_file) then
    inquire (file = basename, exist = found_it, name = file_name2)
    if (found_it) file(i_level)%dir = file(0)%dir
  else
    found_it = .false.
  endif

  if (is_relative .and. .not. found_it) then
    call append_subdirectory (trim(file(i_level-1)%dir), file(i_level)%dir, file(i_level)%dir, err_flag)
    if (err_flag) call parser_error ('BAD DIRECTORY SYNTAX FOR: ' // file_name, stop_here = .true.)
    call append_subdirectory (file(i_level-1)%dir, file_name2, file_name2, err_flag)
  endif

  inquire (file = file_name2, exist = found_it, name = file_name)

  file(i_level)%full_name = file_name
  file(i_level)%f_unit = lunget()

  open (file(i_level)%f_unit, file = file_name, status = 'OLD', action = 'READ', iostat = ios)
  if (ios /= 0 .or. .not. found_it) then
    bp_com%current_file => file(i_level-1)  ! For warning
    if (file_name2 == file_name)  then
      call parser_error ('UNABLE TO OPEN FILE: ' // file_name, stop_here = .true.)
    else
      call parser_error ('UNABLE TO OPEN FILE: ' // file_name, &
                         'THIS FROM THE LOGICAL FILE NAME: ' // file_name_in, stop_here = .true.)
    endif
    do i = 1, i_level-1
      close (file(i_level)%f_unit)
    enddo
    return
  endif

  bp_com%current_file%i_line = 0

  n = size(bp_com%lat_file_names)
  n_file = bp_com%num_lat_files + 1
  if (n < n_file) call re_allocate (bp_com%lat_file_names, n + 100)
  bp_com%num_lat_files = n_file
  inquire (file = file_name, name = bp_com%lat_file_names(n_file))

  ! Note: The same file may be validly called multiple times if it is an inline file.
  ! EG: A wall file called inline.
  ! Therefore the warning is disabled.

!  if (how == 'push') then
!    do i = 1, n_file - 1
!      if (bp_com%lat_file_names(i) /= bp_com%lat_file_names(n_file)) cycle
!      call parser_error ('Same lattice file called multiple times: ' // trim(bp_com%lat_file_names(n_file)), &
!                         level = s_warn$)
!      exit
!    enddo
!  endif

! "pop" means close the current file and pop its name off the stack

case ('pop')
  close (unit = bp_com%current_file%f_unit)
  i_level = i_level - 1
  if (i_level < 0) then
    call parser_error ('BAD "RETURN"')
    return
  elseif (i_level > 0) then
    bp_com%current_file => file(i_level)
    bp_com%calling_file => file(i_level-1)
  else    ! i_level == 0
    if (present(finished)) finished = .true.
  endif

  if (bp_com%inline_call_active) then
    bp_com%parse_line = trim(bp_com%parse_line) // ' ' // file(i_level+1)%parse_line_saved
    bp_com%inline_call_active = file(i_level+1)%inline_call_active
  endif

  bp_com%inline_call_active = .false.

! Programming error

case default
  call parser_error ('INTERNAL ERROR IN PARSER_FILE_STACK SUBROUTINE!')
  if (global_com%exit_on_error) call err_exit
end select

if (present(err)) err = .false.

end subroutine parser_file_stack

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine load_parse_line (action, ix_start, end_of_file) 
!
! Subroutine to load characters from the input file.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
! Input:
!   action -- Character(*): 'continue', 'new_command', or 'init'
!   ix_start  -- Integer: index in bp_com%parse_line string where to append stuff.
!
! Output:
!   end_of_file       -- Logical: 
!   bp_com%parse_line -- String to append to.
!-

recursive subroutine load_parse_line (action, ix_start, end_of_file)

implicit none

integer ix_start, ix, n

character(*) action
character(n_parse_line+20) :: line
character(n_parse_line+20), save :: saved_line
character(1), parameter :: tab = achar(9)

logical :: have_saved_line = .false., end_of_file, flush_this

! action = 'init'

if (action == 'init') then
  have_saved_line = .false.
  bp_com%parse_line = ''
  return
endif

!

end_of_file = .false.
flush_this = .false.

! If 'new_command' then will need to must lines that are part of the rest of the current command. 
! This will happen when there has been an error and the entire command was not parsed.

if (action == 'new_command' .and. .not. have_saved_line) then
  n = len_trim(bp_com%parse_line)
  if (n /= 0) then
    if (index (',+-*/({[=&', bp_com%parse_line(n:n)) /= 0) flush_this = .true.
  endif
endif

!

do

  ! Read a line or use saved_line if it exists

  if (have_saved_line) then
    line = saved_line
    have_saved_line = .false.
  else
    read (bp_com%current_file%f_unit, '(a)', end = 9000) line
    bp_com%current_file%i_line = bp_com%current_file%i_line + 1
    if (line(n_parse_line-ix_start-20:) /= ' ') &
      call parser_error ('INPUT LINE HAS TOO MANY CHARACTERS:', line)
  endif

  ! %input_line1 and %input_line2 are for error messages if needed.
  ! Only the input string being parsed is saved in these lines.
  ! 'new_command' action means we are loading a new input string so start from scratch.
  ! 'continue' action means keep the existing input string.

  if (action == 'continue') then
    bp_com%input_line1 = bp_com%input_line2
    bp_com%input_line2 = line
    bp_com%line1_file_name = bp_com%line2_file_name
    write (bp_com%line2_file_name, '(2a, i0)') trim(bp_com%current_file%full_name), ':', bp_com%current_file%i_line
  elseif (action == 'new_command') then
    bp_com%input_line1 = ''
    bp_com%input_line2 = line
    bp_com%line1_file_name = ''
    write (bp_com%line2_file_name, '(2a, i0)') trim(bp_com%current_file%full_name), ':', bp_com%current_file%i_line
  else
    call parser_error ('INTERNAL ERROR #4: CALL HELP')
    if (global_com%exit_on_error) call err_exit
  endif

  ! strip off comments

  ix = index(line, '!')
  if (ix == 1) then
    line = ' '
  elseif (ix > 1) then
    line = line(:ix-1)
  endif

  ! semi-colon delimiter means that we need to split the line
  ! and save the 2nd piece for the next time around.

  ix = index(line, ';')
  if (ix == 1) then
    have_saved_line = .true.
    saved_line = line(ix+1:)
    line = ' '
  elseif (ix > 1) then
    have_saved_line = .true.
    saved_line = line(ix+1:)
    line = line(:ix-1)
  else
    have_saved_line = .false.
  endif

  ! if the command line is blank then go back for more input

  call string_trim (line, line, ix)
  if (ix /= 0 .or. have_saved_line) exit
enddo

! now simply append the line to %parse_line starting at ix_start

call str_substitute (line)
bp_com%parse_line(ix_start:) = line

! Flush this line if needed

if (flush_this) call load_parse_line (action, ix_start, end_of_file)

return

!

9000  continue
end_of_file = .true.

if (bp_com%parse_line /= '' .and. .not. bp_com%inline_call_active) then
  call parser_error ('FILE ENDED BEFORE PARSING FINISHED', stop_here = .true.)
endif

end subroutine load_parse_line

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function evaluate_array_index (err_flag, delim_list1, word2, delim_list2, delim2) result (this_index)
!
! Function of evaluate the index of an array. Typically the text being parsed looks like:
!      "5) = ..."         or
!      "6).COMP = ..."
!
! Input:
!   delim_list1 -- Character(1): Delimitor after the integer. Normally ')'.
!   delim_list2 -- Character(*): Delimitor list to mark the end of word2. Normally '='.
!
! Output:
!   err_flag    -- Logical: Set True if there is an error. False otherwise.
!   word2       -- Character(*): Word found after delim1. Normally this should be blank.
!   delim2      -- Character(1): Actual delimitor found after word2.
!   this_index  -- Integer: Integer value
!-

function evaluate_array_index (err_flag, delim_list1, word2, delim_list2, delim2) result (this_index)

implicit none

character(*) delim_list1, word2, delim_list2, delim2
integer this_index, ix_word
character(1) delim
character(20) word
logical err_flag, delim_found


! Init

err_flag = .true.
this_index = -1

! Get integer

call get_next_word (word, ix_word, delim_list1, delim, delim_found)
if (.not. delim_found) return

if (.not. is_integer(word)) return
read (word, *) this_index

! Get word after integer

call get_next_word (word2, ix_word, delim_list2, delim2, delim_found)
if (delim_found) err_flag = .false.

end function evaluate_array_index

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function evaluate_logical (word, iostat) result (this_logic)
!
! Function of convert a string into a logical value.
! Accepted possibilities are:
!   .TRUE.  .FALSE. 
!    TRUE    FALSE
!    T       F
!
! Input:
!   word   -- Character(*): Input string.
!
! Output:
!   this_logic -- Logical: Result.
!   iostat     -- Integer: Status: Returns 0 if conversion successful. 
!-

function evaluate_logical (word, iostat) result (this_logic)

implicit none

character(*), intent(in) :: word
character(len(word)+8) wd
logical this_logic
integer, intent(out) :: iostat
integer i

!

iostat = -1
this_logic = .false.  ! To avoid uninit compiler warnings.
if (word == '') return

call str_upcase(wd, word)

do i = 1, len(word)
  if (wd(i:i) == '') cycle

  if (wd(i:i+6) == '.TRUE. ' .or. wd(i:i+4) == 'TRUE ' .or. wd(i:i+1) == 'T ') then
    this_logic = .true.
    iostat = 0
  elseif (wd(i:i+7) == '.FALSE. ' .or. wd(i:i+5) == 'FALSE ' .or. wd(i:i+1) == 'F ') then
    this_logic = .false.
    iostat = 0
  endif
  return
enddo

end function evaluate_logical

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------      
!+
! Subroutine evaluate_value (err_str, value, lat, delim, delim_found, err_flag, end_delims, string_out)
!
! This routine evaluates as a real number the characters at the beginning of bp_com%parse_line.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
! Input:
!   err_str     -- character(*): String to print as part of error message if there is an error.
!   lat         -- lat_struct: 
!   end_delims  -- character(*): List of delimiters that should be present after section of line used for evaluation.
!
! Output:
!   value       -- real(rp):
!   delim       -- character(1): Actual delimiter found. Set to blank is no delim found
!   delim_found -- logical: Set False if end-of-line found instead of a delimiter.
!   err_flag    -- logical:
!   string_out  -- character(*), optional: If present then parsed expression is returned but not evaluated.
!                   Useful for group and overlay control expressions.
!-

subroutine evaluate_value (err_str, value, lat, delim, delim_found, err_flag, end_delims, string_out)

use expression_mod

implicit none

type (lat_struct)  lat
type (expression_atom_struct), allocatable :: stk(:)

real(rp) value

integer i, ix_word, ix_str, n_parens, n_stk

character(*) err_str
character(*), optional :: end_delims
character(*), optional :: string_out
character(1) delim
character(200) str
character(100) err_str2
character(200) word

logical delim_found, ran_function_pending
logical err_flag, call_check

! Get string

call_check = .true.
str = ''
ix_str = 0
n_parens = 0
err_flag = .true.

do
  ! Include "+-" as delims to avoid error with sub-expression exceeding 90 characters and with ending "&" continuation char.
  call get_next_word (word, ix_word, '(),:}+-|', delim, delim_found, upper_case_word = .false., call_check = call_check)
  call_check = .false.
  str = str(1:ix_str) // word
  ix_str = ix_str + ix_word
  if (.not. delim_found) exit

  select case (delim)
  case (',', ')', '}') 
    if (n_parens == 0) exit
  case (':', '|')
    exit
  case default
    ! Nothing to do
  end select

  ix_str = ix_str + 1
  str(ix_str:ix_str) = delim
  if (delim == '(') n_parens = n_parens + 1
  if (delim == ')') n_parens = n_parens - 1
enddo

! Check that final delim matches.

if (present(end_delims)) then
  if (.not. delim_found .or. index(end_delims, delim) == 0) then
    call parser_error ('BAD DELIMITOR AFTER VALUE', 'FOR: ' // err_str)
    return
  endif
endif

! If the string_out argument is present then just return the string for later evaluation

if (present(string_out)) then
  string_out = str
  err_flag = .false.
  return
endif

! If expression is just a number then evaluate and return

if (is_real(str)) then
  read (str, *) value
  err_flag = .false.
  return
endif

! Make a stack

call expression_string_to_stack(str, stk, n_stk,  err_flag, err_str2)
if (err_flag) then
  call parser_error (err_str2, 'FOR: ' // err_str)
  if (err_str2 == 'MALFORMED EXPRESSION') bp_com%parse_line = ''
  return
endif

do i = 1, n_stk
  select case (stk(i)%type)
  case (ran$, ran_gauss$)
    call bp_set_ran_status
  case (variable$)
    call word_to_value (stk(i)%name, lat, stk(i)%value)
  case (species_var$)
    stk(i)%value = species_id(stk(i)%name)
    if (stk(i)%value == invalid$) then
      call parser_error ('INVALID PARTICLE SPECIES: ' // word)
      return
    endif
  end select
enddo

! Evaluate

call evaluate_expression_stack (stk, value, err_flag, err_str2)
if (err_flag) then
  call parser_error (err_str2, 'FOR: ' // err_str)
endif

end subroutine evaluate_value

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine word_to_value (word, lat, value)

implicit none

type (lat_struct), target ::  lat
type (ele_struct), pointer :: ele
type (all_pointer_struct), allocatable :: ptr(:)
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_attribute_struct) attrib_info

integer i, ix1, ix2, ix_word, ios, ix, n_loc, ix_attrib
integer ivar, ixm, ixm2
real(rp) value
real(rp), pointer :: v(:)
character(*) word
character(40) attrib_name, ele_name
character(80) var_name
logical err_flag

! see if this is numeric

if (index('-+.0123456789', word(1:1)) /= 0) then
  read (word, *, iostat = ios) value
  if (ios /= 0) call parser_error ('BAD VARIABLE: ' // word)
  return
endif

! If not numeric...

ix_word = len_trim(word)
var_name = upcase(word)
if (.not. verify_valid_name (var_name, ix_word)) return

! If word does not have a "[...]" then it must be a variable

ix1 = index(var_name, '[')
if (ix1 == 0) then   
  call find_indexx (var_name, bp_com%var%name, bp_com%var%indexx, bp_com%ivar_tot, i)
  if (i == 0) then
    call parser_error ('VARIABLE USED BUT NOT YET DEFINED: ' // word, 'WILL TREAT AS ZERO.', level = s_warn$)
    value = 0
    ! To prevent multiple error messages define this variable.
    bp_com%ivar_tot = bp_com%ivar_tot + 1
    if (bp_com%ivar_tot > size(bp_com%var%name)) call reallocate_bp_com_var()
    ivar = bp_com%ivar_tot
    bp_com%var(ivar)%name = var_name
    bp_com%var(ivar)%value = 0
    call find_indexx (var_name, bp_com%var%name, bp_com%var%indexx, ivar-1, ixm, ixm2)
    bp_com%var(ixm2+1:ivar)%indexx = bp_com%var(ixm2:ivar-1)%indexx
    bp_com%var(ixm2)%indexx = ivar

  else
    value = bp_com%var(i)%value
  endif
  return
endif

! Here if word does have a "[...]" then is a element attribute

ele_name = var_name(:ix1-1)    ! name of attribute

ix2 = index(var_name, ']')
attrib_name = var_name(ix1+1:ix2-1)

if (attrib_name == 'S' .and. bp_com%parser_name /= 'bmad_parser2') then
  call parser_error ('"S" ATTRIBUTE CAN ONLY BE USED WITH BMAD_PARSER2')
endif

! Apertures, etc.

select case (attrib_name)
case ('APERTURE', 'X_LIMIT', 'Y_LIMIT')
  call lat_ele_locator (ele_name, lat, eles, n_loc, err_flag)
  if (.not. err_flag .and. size(eles) > 0) then
    v => eles(1)%ele%value
    if (attrib_name == 'APERTURE') then
      if (v(x1_limit$) /= v(x2_limit$) .or. v(x1_limit$) /= v(y1_limit$) .or. &
          v(x1_limit$) /= v(y2_limit$)) then
        err_flag = .true.
      else
        value = v(x1_limit$)
      endif
    elseif (attrib_name == 'X_LIMIT') then
      if (v(x1_limit$) /= v(x2_limit$)) then
        err_flag = .true.
      else
        value = v(x1_limit$)
      endif
    elseif (attrib_name == 'Y_LIMIT') then
      if (v(y1_limit$) /= v(y2_limit$)) then
        err_flag = .true.
      else
        value = v(y1_limit$)
      endif
    endif
  endif

  ! If size(eles) > 1 then there must be more than one element of the same name.

  if (size(eles) > 1) then
    do i = 2, size(eles)
      if (eles(i)%ele%value(x1_limit$) == eles(1)%ele%value(x1_limit$) .and. &
          eles(i)%ele%value(y1_limit$) == eles(1)%ele%value(y1_limit$)) cycle
      call parser_error (&
            'MULTIPLE ELEMENTS OF THE SAME NAME BUT WITH DIFFERENT ATTRIBUTE VALUES REFERENCED IN: ' // word)
      exit
    enddo
  endif

! Everything else

case default
  call pointers_to_attribute (lat, ele_name, attrib_name, .false., ptr, err_flag, .false., eles, ix_attrib)
  if (err_flag .or. size(ptr) == 0) then
    call parser_error('BAD ATTRIBUTE: ' // word)
    return
  elseif (associated(ptr(1)%r)) then
    value = ptr(1)%r
  elseif (associated(ptr(1)%i)) then
    value = ptr(1)%i
  else  ! Must
    call parser_error('ATTRIBUTE IS NOT REAL OR INTEGER: ' // word)
    return
  endif

  ! If this is bmad_parser, and not bmad_parser2, then dependent attributes have not been set and cannot
  ! be used.

  if (ix_attrib > 0 .and. size(eles) > 0 .and. bp_com%parser_name == 'bmad_parser') then
    attrib_info = attribute_info (eles(1)%ele, ix_attrib)
    if (attrib_info%type == dependent$) then
      call parser_error ('DEPENDENT ATTRIBUTE IS NOT CALCULATED BEFORE LATTICE EXPANSION AND', &
                         'THEREFORE CANNOT BE USED BEFORE ANY EXPAND_LATTICE COMMAND: ' // word)
    endif
  endif

  ! If size(ptr) > 1 then there must be more than one element of the same name.

  if (size(ptr) > 1) then
    do i = 2, size(ptr)
      if (ptr(i)%r /= ptr(1)%r) then
        call parser_error (&
              'MULTIPLE ELEMENTS OF THE SAME NAME BUT WITH DIFFERENT ATTRIBUTE VALUES REFERENCED IN: ' // word)
        exit
      endif
    enddo
  endif

end select

end subroutine word_to_value

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_add_variable (word, lat)

implicit none

type (lat_struct) lat
character(*) word
character(1) delim
real(rp) old_val
integer i, ivar, ixm, ixm2
logical delim_found, err_flag

!

call find_indexx (word, bp_com%var%name, bp_com%var%indexx, bp_com%ivar_tot, i)
if (i /= 0) then
  old_val = bp_com%var(i)%value
  call evaluate_value (word, bp_com%var(i)%value, lat, delim, delim_found, err_flag)

  if (bp_com%var(i)%value == old_val) then
    call parser_error ('VARIABLES ARE NOT ALLOWED TO BE REDEFINED: ' // word, 'BUT SINCE OLD_VALUE = NEW_VALUE THIS IS ONLY A WARNING...', level = s_warn$)
  else
    call parser_error ('VARIABLES ARE NOT ALLOWED TO BE REDEFINED: ' // word)
  endif

else
  bp_com%ivar_tot = bp_com%ivar_tot + 1
  if (bp_com%ivar_tot > size(bp_com%var%name)) call reallocate_bp_com_var()
  ivar = bp_com%ivar_tot
  bp_com%var(ivar)%name = word
  ! Reindex.
  call find_indexx (word, bp_com%var%name, bp_com%var%indexx, ivar-1, ixm, ixm2)
  bp_com%var(ixm2+1:ivar)%indexx = bp_com%var(ixm2:ivar-1)%indexx
  bp_com%var(ixm2)%indexx = ivar
  ! Evaluate
  call evaluate_value (bp_com%var(ivar)%name, bp_com%var(ivar)%value, lat, delim, delim_found, err_flag)
endif

if (delim_found .and. .not. err_flag) call parser_error  &
                  ('EXTRA CHARACTERS ON RHS: ' // bp_com%parse_line,  &
                   'FOR VARIABLE: ' // bp_com%var(ivar)%name)

end subroutine parser_add_variable

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine bmad_parser_type_get (ele, attrib_name, delim, delim_found, name, pele, str_out)

implicit none

type (ele_struct)  ele
type (parser_ele_struct), optional :: pele

integer ix, ix_word, n

character(*) attrib_name
character(*), optional :: name, str_out
character(40)  word
character(1)   delim, str_end
character(200) type_name

logical delim_found, err

!

call string_trim(bp_com%parse_line, bp_com%parse_line, ix)

str_end = bp_com%parse_line(1:1)
if (str_end == '"' .or. str_end == "'") then
  bp_com%parse_line = bp_com%parse_line(2:)
  ix = index(bp_com%parse_line, str_end)
  if (ix == 0) then
    call parser_error ('MISSING ENDING QUOTE MARK FOR: ' // attrib_name,  &
                        'FOR ELEMENT: ' // ele%name)
    type_name = ' '
  else
    type_name = bp_com%parse_line(1:ix-1)
    bp_com%parse_line = bp_com%parse_line(ix+1:)
    call get_next_word (word, ix_word, ',=', delim, delim_found, .true.)
    if (ix_word /= 0) call parser_error ( &
              'EXTRA CHARACTERS FOUND AFTER VALUE OF: ' // attrib_name, &
              'I DO NOT UNDERSTAND: ' // word,  &
              'FOR ELEMENT: ' // ele%name)
  endif
else
  call get_next_word (type_name, ix_word, ',= ', delim, delim_found, .false.)
endif

select case (attrib_name)
case ('ALIAS')
  ele%alias = type_name
case ('CRYSTAL_TYPE', 'MATERIAL_TYPE', 'PHYSICAL_SOURCE')
  ele%component_name = type_name
case ('DESCRIP', 'LATTICE')
  if (.not. associated(ele%descrip)) allocate (ele%descrip) 
  ele%descrip = type_name
case ('LR_WAKE_FILE') 
  call read_lr_wake (ele, type_name)
case ('MATERIAL', 'CLEAR_MATERIAL', 'OPAQUE_MATERIAL')
  str_out = type_name
case ('TO', 'TO_LINE', 'ORIGIN_ELE')
  ele%component_name = type_name
  call upcase_string (ele%component_name)
case ('TO_ELEMENT')
  pele%ele_name = type_name
  call upcase_string (pele%ele_name)
case ('SR_WAKE_FILE') 
  call read_sr_wake (ele, type_name)
case ('TYPE')
  ele%type = type_name
case ('CUSTOM_ATTRIBUTE1', 'CUSTOM_ATTRIBUTE2', 'CUSTOM_ATTRIBUTE3', &
      'CUSTOM_ATTRIBUTE4', 'CUSTOM_ATTRIBUTE5')
  name = type_name
  call upcase_string (name)

case default
  call parser_error ('INTERNAL ERROR IN BMAD_PARSER_TYPE_GET: I NEED HELP!')
  if (global_com%exit_on_error) call err_exit
end select

end subroutine bmad_parser_type_get

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine read_lr_wake (ele, lr_file_name)
!
! Subroutine to read in a long-range wake field from an external file.
! This subroutine is used by bmad_parser and bmad_parser2.
!
! Input:
!   ele          -- Ele_struct: Element containing wake structure.
!   lr_file_name -- Character(*):  Name of long-range wake field file.
!
! Output:
!   ele -- Ele_struct: Element with wake information.
!     %wake%lr_mode(:)       -- Long-range wake potential.
!-
      
subroutine read_lr_wake (ele, lr_file_name)

implicit none

type lr_wake_input_struct
  real(rp) freq         ! Actual Frequency in Hz.
  real(rp) R_over_Q     ! Strength in V/C/m^2.
  real(rp) Q            ! Quality factor.
  integer m             ! Order (1 = dipole, 2 = quad, etc.)
  character(16) angle   ! polarization angle (radians/2pi).
  real(rp) b_sin, b_cos, a_sin, a_cos, t_ref
end type

type (ele_struct) ele
type (lr_wake_input_struct) lr(500)
integer n_row, iu, i, j, ios
character(*) lr_file_name
character(200) full_file_name
logical set_done, finished, err

namelist / long_range_modes / lr

! Init

if (.not. associated(ele%wake)) allocate (ele%wake)
if (.not. allocated(ele%wake%sr_long%mode))  allocate (ele%wake%sr_long%mode(0))
if (.not. allocated(ele%wake%sr_trans%mode)) allocate (ele%wake%sr_trans%mode(0))
if (.not. allocated(ele%wake%lr_spline)) allocate (ele%wake%lr_spline(0))
if (allocated(ele%wake%lr_mode)) deallocate (ele%wake%lr_mode)

! get data

call parser_file_stack ('push', lr_file_name, finished, err)
if (err) return
iu = bp_com%current_file%f_unit

ele%wake%lr_file = lr_file_name

lr%freq = -1
lr%angle = ''
lr%b_sin = 0
lr%b_cos = 0
lr%a_sin = 0
lr%a_cos = 0
lr%t_ref = 0

read (iu, nml = long_range_modes, iostat = ios)
call parser_file_stack ('pop')

if (ios > 0 .or. lr(1)%freq == -1) then
  call parser_error ('CANNOT READ LONG_RANGE_MODES NAMELIST FOR ELEMENT: ' // ele%name, & 
                'FROM FILE: '// full_file_name)
  return
endif

! Transfer info to ele structure.

n_row = count(lr%freq /= -1)
allocate (ele%wake%lr_mode(n_row))
j = 0
do i = 1, size(lr)
  if (lr(i)%freq == -1) cycle

  j = j + 1
  ele%wake%lr_mode(j)%freq_in   = lr(i)%freq
  ele%wake%lr_mode(j)%freq      = lr(i)%freq
  ele%wake%lr_mode(j)%r_over_q  = lr(i)%r_over_q
  ele%wake%lr_mode(j)%q         = lr(i)%q
  ele%wake%lr_mode(j)%m         = lr(i)%m
  ele%wake%lr_mode(j)%b_sin     = lr(i)%b_sin
  ele%wake%lr_mode(j)%b_cos     = lr(i)%b_cos
  ele%wake%lr_mode(j)%a_sin     = lr(i)%a_sin
  ele%wake%lr_mode(j)%a_cos     = lr(i)%a_cos
  ele%wake%lr_mode(j)%t_ref     = lr(i)%t_ref

  call downcase_string(lr(i)%angle)
  if (lr(i)%angle == '') then
    call parser_error ('LONG_RANGE_MODE ANGLE IS MISSING. MUST BE NUMBER OR "UNPOLARIZED"', & 
                  'FOR ELEMENT: ' // ele%name, &
                  'IN FILE: ' // full_file_name)
    cycle
  endif

  if (index('unpolarized', trim(lr(j)%angle)) == 1) then
    ele%wake%lr_mode(j)%polarized = .false.
    ele%wake%lr_mode(j)%angle     = 0
  else
    ele%wake%lr_mode(j)%polarized = .true.
    read (lr(j)%angle, *, iostat = ios) ele%wake%lr_mode(j)%angle
    if (ios /= 0) then
      call parser_error ('BAD LONG_RANGE_MODE ANGLE.', &
                    'FOR ELEMENT: ' // ele%name, &
                    'IN FILE: ' // full_file_name)
      cycle
    endif
  endif
enddo

call randomize_lr_wake_frequencies (ele, set_done)
if (set_done) call bp_set_ran_status

end subroutine read_lr_wake

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine read_sr_wake (ele, sr_file_name)
!
! Subroutine to read in a short-range wake field from an external file.
! This subroutine is used by bmad_parser and bmad_parser2.
!
! Input:
!   ele          -- Ele_struct: Element containing wake structure.
!   sr_file_name -- Character(*):  Name of short-range wake field file.
!
! Output:
!   ele -- Ele_struct: Element with wake information.
!     %wake%sr_table(:)       -- Short-range wake potential.
!     %wake%sr_long(:)  -- Short-range wake potential.
!     %wake%sr_trans(:) -- Short-range wake potential.
!-

subroutine read_sr_wake (ele, sr_file_name)

implicit none

type (ele_struct) ele
type (wake_sr_mode_struct), target :: longitudinal(100), transverse(100)
type (wake_sr_mode_struct), pointer :: sr(:), sr1

real(rp) z_max
integer n, j, iu, ios, ix, i, ixx

character(*) sr_file_name
character(140) line, line_in
character(200) full_file_name

logical in_namelist, finished, err

! init

if (.not. associated(ele%wake))   allocate (ele%wake)
if (.not. allocated(ele%wake%lr_mode)) allocate (ele%wake%lr_mode(0))
if (.not. allocated(ele%wake%lr_spline)) allocate (ele%wake%lr_spline(0))
if (allocated(ele%wake%sr_long%mode))  deallocate (ele%wake%sr_long%mode)
if (allocated(ele%wake%sr_trans%mode)) deallocate (ele%wake%sr_trans%mode)

! Open file

ele%wake%sr_file = sr_file_name
call parser_file_stack ('push', sr_file_name, finished, err)
if (err) return
iu = bp_com%current_file%f_unit

! Get data

longitudinal = wake_sr_mode_struct()
longitudinal%phi = real_garbage$

transverse = wake_sr_mode_struct()
transverse%phi = real_garbage$

z_max = real_garbage$
in_namelist = .false.

do
  read (iu, '(a)', iostat = ios) line_in
  if (ios /= 0) then
    call parser_error ('END OF FILE REACHED BEFORE SHORT_RANGE_MODES NAMELIST PARSED.', &
                       'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
    return
  endif

  line = line_in
  ix = index(line, '!')
  if (ix /= 0) line = line(1:ix-1)
  call string_trim(line, line, ixx)
  if (ixx == 0) cycle
  call downcase_string (line) 

  if (in_namelist) then
    if (line(1:1) == '/') exit
  else
    if (line == '&short_range_modes') then
      in_namelist = .true.
      cycle
    endif
  endif

  if (line(1:12) == 'longitudinal') then
    call string_trim(line(13:), line, ixx)
    sr => longitudinal
    if (.not. get_this_sr1 (sr, sr1)) return
    sr1%transverse_dependence = none$    ! Default
    sr1%polarization = none$             ! Default

    if (.not. expect_equal_sign()) return
    if (.not. get_this_param (sr1%amp)) return
    if (.not. get_this_param (sr1%damp)) return
    if (.not. get_this_param (sr1%k)) return
    if (.not. get_this_param (sr1%phi)) return
    if (.not. get_this_switch (sr1%polarization, sr_polarization_name)) return
    if (.not. get_this_switch (sr1%transverse_dependence, sr_transverse_dependence_name)) return
    if (.not. expect_nothing ()) return

  elseif (line(1:10) == 'transverse') then
    call string_trim(line(11:), line, ixx)
    sr => transverse
    if (.not. get_this_sr1 (sr, sr1)) return
    sr1%transverse_dependence = linear_leading$     ! Default

    if (.not. expect_equal_sign()) return
    if (.not. get_this_param (sr1%amp)) return
    if (.not. get_this_param (sr1%damp)) return
    if (.not. get_this_param (sr1%k)) return
    if (.not. get_this_param (sr1%phi)) return
    if (.not. get_this_switch (sr1%polarization, sr_polarization_name)) return
    if (.not. get_this_switch (sr1%transverse_dependence, sr_transverse_dependence_name)) return
    if (.not. expect_nothing ()) return

  elseif (line(1:ixx) == 'z_max') then
    call string_trim(line(ixx+1:), line, ixx)
    if (.not. expect_equal_sign()) return
    if (.not. get_this_param (z_max)) return
    if (.not. expect_nothing ()) return

  else
    call parser_error ('BAD PARAMETER NAME IN SHORT_RANGE_MODES NAMELIST: ' // line_in, &
                       '(MUST BE "LONGITUDINAL", "TRANSVERSE", OR "Z_MAX")', &
                       'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
    return
  endif

enddo

call parser_file_stack ('pop')

! Put data in element

n = count(longitudinal%phi /= real_garbage$)
allocate (ele%wake%sr_long%mode(n))
ele%wake%sr_long%mode = longitudinal(1:n)
if (any(longitudinal(1:n)%phi == real_garbage$)) call parser_error ( &
    'JUMBLED INDEX FOR LONGITUDINAL SHORT_RANGE_MODES FROM FILE: ' // full_file_name, &
    'FOR ELEMENT: ' // ele%name)

n = count(transverse%phi /= real_garbage$)
allocate (ele%wake%sr_trans%mode(n))
ele%wake%sr_trans%mode = transverse(1:n)
if (any(transverse(1:n)%phi == real_garbage$)) call parser_error ( &
    'JUMBLED INDEX FOR TRANSVERSE SHORT_RANGE_MODES FROM FILE: ' // full_file_name, &
    'FOR ELEMENT: ' // ele%name)


ele%wake%z_sr_max = z_max
if (z_max == real_garbage$) call parser_error ( &
    'Z_MAX NOT SET FOR SHORT_RANGE_MODES FROM FILE: ' // full_file_name, &
    'FOR ELEMENT: ' // ele%name)

!-------------------------------------------------------------------------
contains

function expect_equal_sign () result (is_ok)

logical is_ok

is_ok = .false.

if (line(1:1) /= '=') then
  call parser_error ('MISING "=" SIGN IN SHORT_RANGE_MODES NAMELIST: ' // line_in, &
                     'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

call string_trim(line(2:), line, ixx)
is_ok = .true.

end function expect_equal_sign

!-------------------------------------------------------------------------
! contains

function get_this_sr1 (sr, sr1) result (is_ok)

type (wake_sr_mode_struct), pointer :: sr(:), sr1
integer ixp, ios, n
logical is_ok


!
is_ok = .false.

if (line(1:1) /= '(') then
  call parser_error ('MISING "(" IN SHORT_RANGE_MODES NAMELIST: ' // line_in, &
                     'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

ixp = index(line, ')') 
if (ixp == 0) then
  call parser_error ('MISING "(" IN SHORT_RANGE_MODES NAMELIST: ' // line_in, &
                     'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

read(line(2:ixp-1), *, iostat = ios) n
if (ios /= 0 .or. n < 1 .or. n > size(sr)) then
  call parser_error ('BAD WAKE MODE INDEX IN SHORT_RANGE_MODES NAMELIST: ' // line_in, &
                     'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

call string_trim (line(ixp+1:), line, ixx)
sr1 => sr(n)

is_ok = .true.

end function get_this_sr1

!-------------------------------------------------------------------------
! contains

function get_this_param (param) result (is_ok)

real(rp) param
integer ios
logical is_ok

!

is_ok = .false.

if (ixx == 0) then
  call parser_error ('MISING NUMBER IN SHORT_RANGE_MODES NAMELIST: ' // line_in, &
                     'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

read (line(1:ixx), *, iostat = ios) param

if (ios /= 0) then
  call parser_error ('MALFORMED NUMBER IN SHORT_RANGE_MODES NAMELIST: ' // line_in, &
                     'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

call string_trim(line(ixx+1:), line, ixx)

is_ok = .true.

end function get_this_param

!-------------------------------------------------------------------------
! contains

function get_this_switch (switch, names) result (is_ok)

integer switch, ix
character(*) :: names(:)
logical is_ok

is_ok = .false.

if (ixx == 0) then
  is_ok = .true.
  return
endif

call match_word (line(1:ixx), names, switch)

if (switch < 1) then
  call parser_error ('BAD SWITCH NAME IN SHORT_RANGE_MODES NAMELIST: ' // line_in, &
                     'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

call string_trim(line(ixx+1:), line, ixx)

is_ok = .true.

end function get_this_switch

!-------------------------------------------------------------------------
! contains

function expect_nothing () result (is_ok)

logical is_ok

is_ok = .false.

if (line /= '') then    
  call parser_error ('EXTRA STUFF ON LINE: ' // line_in, &
                     'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif


is_ok = .true.

end function expect_nothing

end subroutine read_sr_wake

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_list_of_names (ele, err_str, delim, delim_found)
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-
      
subroutine get_list_of_names (ele, err_str, name_list, delim, delim_found, err_flag)

implicit none

type (ele_struct)  ele

character(*) err_str
character(*), allocatable :: name_list(:)
character(1) delim
character(40) word

integer ix_word, n_name 

logical delim_found, curly_parens_found, err_flag

! Opening "}" is optional if there is only one word

err_flag = .true.

call get_next_word (word, ix_word, '{,}()', delim, delim_found, .true.)
curly_parens_found = (delim == '{')
if (curly_parens_found .and. ix_word /= 0) then
  call parser_error ('ERROR PARSING: ' // err_str, 'FOUND STUFF BEFORE OPENING "{"')
  return
endif

if (curly_parens_found) call get_next_word (word, ix_word, '{,}()', delim, delim_found, .true.)

n_name = 0
call re_allocate (name_list, 10)

do
  n_name = n_name + 1
  if (n_name > size(name_list)) call re_allocate(name_list, n_name+10)
  name_list(n_name) = word

  if ((delim == '}' .and. .not. curly_parens_found) .or. (.not. delim_found .and. curly_parens_found)) then
    call parser_error ('ERROR PARSING: ' // err_str, 'MISMATCHED {} BRACES')
    return
  endif

  if (delim_found .and. delim /= '}' .and. delim /= ',') then
    call parser_error ('ERROR PARSING: ' // err_str, 'MALFORMED STATEMENT')
    return
  endif

  if (delim == '}' .or. .not. curly_parens_found) then 
    if (delim == '}') call get_next_word (word, ix_word, ',=:', delim, delim_found, .true.)
    exit
  endif

  call get_next_word (word, ix_word, '{,}()', delim, delim_found, .true.)
enddo

call re_allocate(name_list, n_name)
err_flag = .false.

end subroutine get_list_of_names

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_overlay_group_names (ele, lat, pele, delim, delim_found, is_control_var_list)
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-
      
subroutine get_overlay_group_names (ele, lat, pele, delim, delim_found, is_control_var_list)

implicit none

type (ele_struct)  ele
type (parser_ele_struct), target :: pele
type (lat_struct)  lat
type (parser_controller_struct), pointer :: pc

real(rp) value

integer ix_word, n_slave, i, j, k
                           
character(1) delim
character(40) word_in, word
character(40), allocatable :: name(:), attrib_name(:)
character(200), allocatable :: expression(:)
character(200) err_str

logical :: is_control_var_list
logical delim_found, err_flag, end_of_file, ele_names_only

!

allocate (name(40), attrib_name(40), expression(40))

call get_next_word (word_in, ix_word, '{,}', delim, delim_found, .true.)
if (delim /= '{' .or. ix_word /= 0) call parser_error  &
        ('BAD ' // control_name(ele%lord_status) // 'SPEC: ' // word_in,  &
        'FOR ELEMENT: ' // ele%name)

ele_names_only = (is_control_var_list .or. ele%key == girder$)

! loop over all names in "{...}" list

n_slave = 0
do 

  call get_next_word (word_in, ix_word, '{,}/:', delim, delim_found, .true.)
  if (delim == ':' .and. ele%key == girder$) pele%is_range = .true.

  ! If "{}" with no slaves... 
  if (delim == '}' .and. ix_word == 0 .and. n_slave == 0) then
    call get_next_word (word, ix_word, ',=:', delim, delim_found, .true.)
    exit
  endif

  n_slave = n_slave + 1
  word = word_in

  if (n_slave > size(name)) then
    call re_allocate(name, 2*n_slave)
    call re_allocate(attrib_name, 2*n_slave)
    call re_allocate(expression, 2*n_slave)
  endif

  j = index(word, '[')
  if (j > 1) then
    k = index(word, ']')
    if (k <= j+1) then
      call parser_error ('BAD ATTRIBUTE SPEC: ' // word_in, 'FOR: ' // ele%name)
      word = word(:k-1) // word(j+1:)
    else
      attrib_name(n_slave) = word(j+1:k-1)
      word = word(:j-1) // word(k+1:)
    endif
  else
    attrib_name(n_slave) = blank_name$
  endif

  name(n_slave) = word
  if (word == '') call parser_error ('SLAVE ELEMENT NAME MISSING WHEN PARSING LORD: ' // ele%name)

  ! If ele_names_only = True then evaluating "var = {...}" construct or is a girder.
  ! In this case, there are no expressions
  
  expression(n_slave) = '1'

  if (delim == '/' .or. (delim == ':' .and. ele%key /= girder$)) then
    if (ele_names_only) then
      call parser_error ('BAD VAR = {...} CONSTRUCT.', 'FOR: ' // ele%name)
      return
    else
      call evaluate_value (trim(ele%name), value, lat, delim, delim_found, err_flag, ',}', expression(n_slave))
      if (err_flag) then
        call parser_error ('BAD EXPRESSION: ' // word_in,  'FOR ELEMENT: ' // ele%name)
        call load_parse_line ('new_command', 1, end_of_file)         ! next line
        return
      endif
    endif
  endif

  if (delim == '}') then
    call get_next_word (word, ix_word, ',=:', delim, delim_found, .true.)
    exit
  elseif (delim /= ',' .and. .not. (delim == ':' .and. ele%key == girder$)) then
    call parser_error ('BAD ' // control_name(ele%lord_status) //  &
            'SPEC: ' // word_in, 'FOR: ' // ele%name)
    exit
  endif
                        
enddo

! if (n_slave == 0) call parser_error ( &
!        'NO SLAVE ELEMENTS ASSOCIATED WITH GROUP/OVERLAY ELEMENT: ' // ele%name)

if (is_control_var_list) then
  allocate(ele%control_var(n_slave))
  ele%control_var%name = name(1:n_slave)
else
  allocate (pele%control(n_slave))
  pele%control%name = name(1:n_slave)
  pele%control%attrib_name = attrib_name(1:n_slave)
  do i = 1, n_slave
    ! Expression string to stack
    pc => pele%control(i)
    call expression_string_to_stack (expression(i), pc%stack, pc%n_stk, err_flag, err_str)
    if (err_flag) then
      call parser_error (err_str, 'FOR ELEMENT: ' // ele%name, 'EXPRESSION: ' // trim(expression(i)))
    endif
  enddo
endif

end subroutine get_overlay_group_names

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-


function verify_valid_name (name, ix_name) result (is_valid)

implicit none

integer i, ix_name, ix1, ix2

character(*) name
character(27), parameter :: letters = '\ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
character(44), parameter :: valid_chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ\0123456789_[]().#'
character(1), parameter :: tab = achar(9)

logical OK, is_valid

! check for blank spaces

is_valid = .false.

do i = 1, min(ix_name, len(name))
  if (name(i:i) == ' ' .or. name(i:i) == tab) then
    call parser_error  ('NO DELIMITER BETWEEN NAMES: ' // name)
    return
  endif
enddo

! check for name too long

if (ix_name > len(name)) then
   call parser_error ('NAME TOO LONG: ' // name)
   ix_name = len(name)      ! chop name
  return
endif

! check for name too short

if (ix_name == 0) then
  call parser_error ('BLANK NAME')
  return
endif

! check for invalid characters in name

OK = .true.
if (index(letters, name(1:1)) == 0) OK = .false.

do i = 1, ix_name
  if (index(valid_chars, name(i:i)) == 0) OK = .false.
enddo

if (.not. OK) then
  call parser_error ('INVALID NAME: UNRECOGNIZED CHARACTERS IN: ' // name)
  return
endif

! Check for non-matched "(" ")" pairs

ix1 = index(name, '(')
ix2 = index(name, ')')
if (ix1 /= 0 .or. ix2 /= 0) then
  if (ix1 == 0) then
    call parser_error ('UNMATCHED PARENTHESIS: ' // name)
    return
  endif
  if (ix2 <= ix1+1) then
    call parser_error  ('INVALID: REVERSED PARENTHESES: ' // name)
    return
  endif
  if (index(name(ix1+1:), '(') /= 0 .or. index(name(ix2+1:), ')') /=  0) then
    call parser_error ('INVALID: BAD PARENTHESES: ' // name)
    return
  endif
endif

! Check for non matched "[" "]" pairs

ix1 = index(name, '[')
ix2 = index(name, ']')

if (ix1 /= 0 .or. ix2 /= 0) then
  if (ix1 == 0) then 
    call parser_error ('UNMATCHED BRACKET: ' // name)
    return
  endif

  if (ix2 <= ix1+1) then
    call parser_error  ('INVALID: REVERSED BRACKETS: ' // name)
    return
  endif

  if (index(name(ix1+1:), '[') /= 0 .or. index(name(ix2+1:), ']') /=  0) then
    call parser_error ('INVALID: BAD BRACKETS: ' // name)
    return
  endif

  if (ix2 /= len(name)) then
    if (name(ix2+1:ix2+1) /= ' ') then
      call parser_error  ('INVALID: SOMETHING AFTER CLOSING "]" BRACKET: ' // name)
      return
    endif
  endif

endif

! check for more than 40 characters

if ((ix1 == 0 .and. ix_name > 40) .or. (ix1 > 41 .or. ix2 - ix1 > 41)) then
  call parser_error ('NAME HAS > 40 CHARACTERS: ' // name)
  return
endif

is_valid = .true.

end function verify_valid_name

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_error (what1, what2, what3, what4, seq, pele, stop_here, level, r_array, i_array)

implicit none

type (seq_struct), optional :: seq
type (parser_ele_struct), optional :: pele

character(*) what1
character(*), optional :: what2, what3, what4
character(160) lines(12)
character(16), parameter :: r_name = 'parser_error'
real(rp), optional :: r_array(:)
integer, optional :: level
integer nl, err_level
integer, optional :: i_array(:)
logical, optional :: stop_here

! bp_com%error_flag is a common logical used so program will stop at end of parsing

err_level = integer_option(s_error$, level)

if (bp_com%print_err) then

  nl = 0

  select case (err_level)
  case (s_info$)
    nl=nl+1; lines(nl) = 'Note from: ' // trim(bp_com%parser_name) // ': ' // trim(what1)
  case (s_warn$)
    nl=nl+1; lines(nl) = 'WARNING IN ' // trim(bp_com%parser_name) // ': ' // trim(what1)
  case (s_error$)
    nl=nl+1; lines(nl) = 'ERROR IN ' // trim(bp_com%parser_name) // ': ' // trim(what1)
  end select

  if (present(what2)) then
    nl=nl+1; lines(nl) = '     ' // trim(what2)
  endif

  if (present(what3)) then
    nl=nl+1; lines(nl) = '     ' // trim(what3)
  endif

  if (present(what4)) then
    nl=nl+1; lines(nl) = '     ' // trim(what4)
  endif

  if (present(seq)) then
    nl=nl+1; lines(nl) = '      IN FILE: ' // trim(seq%file_name)
    nl=nl+1; write (lines(nl), '(a, i0)') '      AT LINE: ', seq%ix_line
  elseif (bp_com%current_file%full_name /= ' ') then
    if (bp_com%input_line_meaningful) then
      nl=nl+1; lines(nl) = '      IN FILE: ' // trim(bp_com%current_file%full_name)
      nl=nl+1; write (lines(nl), '(a, i0)') '      AT OR BEFORE LINE: ', bp_com%current_file%i_line
    else
      nl=nl+1; lines(nl) = '      ROOT FILE: ' // trim(bp_com%current_file%full_name)
    endif
  endif

  if (bp_com%input_line_meaningful) then
    if (len_trim(bp_com%input_line1) /= 0) then
      nl=nl+1; lines(nl) = '     ' // trim(bp_com%input_line1)
    endif
    if (len_trim(bp_com%input_line2) /= 0) then
      nl=nl+1; lines(nl) = '     ' // trim(bp_com%input_line2)
    endif
  endif

  if (present(pele)) then
    nl=nl+1; lines(nl) = '      ELEMENT DEFINED IN FILE: ' // trim(pele%lat_file)
    nl=nl+1; write (lines(nl), '(a, i0)') '      AT LINE: ', pele%ix_line_in_file
  endif

  nl=nl+1; lines(nl) = ''

  call out_io (err_level, r_name, lines(1:nl), r_array = r_array, i_array = i_array, insert_tag_line = .false.)

endif

! Warnings do not result in bp_com%error_flag being set. Just no digested file is generated.

if (err_level == s_warn$) then
  bp_com%write_digested = .false.
elseif (err_level == s_error$) then
  bp_com%error_flag = .true.
  if (logic_option(.false., stop_here)) then
    if (global_com%exit_on_error) stop
    bp_com%fatal_error_flag = .true.
  endif
endif

end subroutine parser_error

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine init_bmad_parser_common

implicit none

integer nn, i

! 

if (allocated(bp_com%var)) deallocate (bp_com%var)

nn = 32  ! number of standard (non-user defined) constants
allocate (bp_com%var(nn))

bp_com%var( 1) = bp_var_struct('PI', pi, 0)
bp_com%var( 2) = bp_var_struct('TWOPI', twopi, 0)
bp_com%var( 3) = bp_var_struct('FOURPI', fourpi, 0)
bp_com%var( 4) = bp_var_struct('E_LOG', 2.718281828459_rp, 0)
bp_com%var( 5) = bp_var_struct('SQRT_2', sqrt_2, 0)
bp_com%var( 6) = bp_var_struct('DEGRAD', 180 / pi, 0)
bp_com%var( 7) = bp_var_struct('DEGREES', pi / 180, 0) ! From degrees to radians.
bp_com%var( 8) = bp_var_struct('RADDEG', pi / 180, 0)
bp_com%var( 9) = bp_var_struct('M_ELECTRON', m_electron, 0)
bp_com%var(10) = bp_var_struct('M_MUON', m_muon, 0)
bp_com%var(11) = bp_var_struct('M_PION_0', m_pion_0)
bp_com%var(12) = bp_var_struct('M_PION_CHARGED', m_pion_charged, 0)
bp_com%var(13) = bp_var_struct('M_PROTON', m_proton, 0)
bp_com%var(14) = bp_var_struct('M_DEUTERON', m_deuteron, 0)
bp_com%var(15) = bp_var_struct('C_LIGHT', c_light, 0)
bp_com%var(16) = bp_var_struct('R_E', r_e, 0)
bp_com%var(17) = bp_var_struct('R_P', r_p, 0)
bp_com%var(18) = bp_var_struct('E_CHARGE', e_charge, 0)
bp_com%var(19) = bp_var_struct('H_PLANCK', h_planck, 0)
bp_com%var(20) = bp_var_struct('H_BAR_PLANCK', h_bar_planck, 0)
bp_com%var(21) = bp_var_struct('PMASS', p_mass, 0)
bp_com%var(22) = bp_var_struct('EMASS', e_mass, 0)
bp_com%var(23) = bp_var_struct('CLIGHT', c_light, 0)
bp_com%var(24) = bp_var_struct('FINE_STRUCT_CONST', fine_structure_constant, 0)
bp_com%var(25) = bp_var_struct('ANOM_MAG_ELECTRON', anomalous_mag_moment_electron, 0)  ! Old style. Deprecated.
bp_com%var(26) = bp_var_struct('ANOM_MAG_PROTON', anomalous_mag_moment_proton, 0)      ! Old style. Deprecated.
bp_com%var(27) = bp_var_struct('ANOM_MAG_MUON', anomalous_mag_moment_muon, 0)          ! Old style. Deprecated.
bp_com%var(28) = bp_var_struct('ANOM_MAG_DEUTERON', anomalous_mag_moment_deuteron, 0)  ! Old style. Deprecated.
bp_com%var(29) = bp_var_struct('ANOM_MOMENT_ELECTRON', anomalous_mag_moment_electron, 0)
bp_com%var(20) = bp_var_struct('ANOM_MOMENT_PROTON', anomalous_mag_moment_proton, 0)
bp_com%var(31) = bp_var_struct('ANOM_MOMENT_MUON', anomalous_mag_moment_muon, 0)
bp_com%var(32) = bp_var_struct('ANOM_MOMENT_DEUTERON', anomalous_mag_moment_deuteron, 0)

bp_com%ivar_init = nn
bp_com%ivar_tot  = nn

call indexx (bp_com%var(1:nn)%name, bp_com%var(1:nn)%indexx)

end subroutine init_bmad_parser_common

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine add_this_multipass (lat, m_slaves, lord_in)

use multipass_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: slave, lord, slave2, lord2, ele
type (ele_struct), optional :: lord_in
type (branch_struct), pointer :: branch
type (lat_ele_loc_struct) m_slaves(:)

integer i, j, k, n, i1, ix, ixc, ixic, ix_lord, ixb, ixb2, ix_n
integer n_multipass, ic, ix_l1, ix_l0, ix_pass, n_links, lmax

character(40) base_name
character(100) slave2_name

! Count slaves.
! If i > lat%n_ele_track we are looking at cloning a super_lord which should
! not happen.

n_multipass = size(m_slaves)

! setup multipass_lord

call new_control (lat, ix_lord)
lord => lat%ele(ix_lord)

if (present(lord_in)) then
  lord = lord_in   ! Use lord_in as template
else
  lord = pointer_to_ele(lat, m_slaves(1))  ! Set attributes equal to first slave.
endif

lord%logic = .false.  ! So parser_add_superimpose will not try to use as ref ele.
lord%lord_status = multipass_lord$
lord%n_slave = 0
lord%ix1_slave = 0
call add_lattice_control_structs (lord, n_add_slave = n_multipass)

! Multipass_lord does not have reference energy or s_position or map bookkeeping. 

lord%bookkeeping_state%ref_energy = ok$   
lord%bookkeeping_state%s_position = ok$   
lord%bookkeeping_state%mat6       = ok$   

! A multipass lord defaults to n_ref_pass = 1 if neither n_ref_pass, p0c and e_tot are not set.

if (nint(lord%value(n_ref_pass$)) == 0 .and. lord%value(p0c$) == 0 .and. lord%value(e_tot$) == 0) then
  lord%value(n_ref_pass$) = 1
endif

! Setup bookkeeping between lord and slaves

do i = 1, n_multipass
  slave => pointer_to_ele (lat, m_slaves(i))
  ixc = i + lord%ix1_slave - 1
  lat%control(ixc)%lord%ix_ele = ix_lord
  lat%control(ixc)%slave = lat_ele_loc_struct(slave%ix_ele, slave%ix_branch)
  if (slave%n_lord /= 0) then
    call parser_error ('INTERNAL ERROR: CONFUSED MULTIPASS SETUP.', &
                  'PLEASE GET EXPERT HELP!')
    if (global_com%exit_on_error) call err_exit
  endif
  write (slave%name, '(2a, i0)') trim(slave%name), '\', i   ! '
  call add_lattice_control_structs (slave, n_add_lord = 1)
  slave%slave_status = multipass_slave$
  ixic = slave%ic1_lord
  lat%ic(ixic) = ixc
  ! If slave is a super_lord then create the appropriate super_slave names
  i1 = 0
  do j = 1, slave%n_slave
    slave2 => pointer_to_slave(slave, j)
    if (slave2%n_lord == 1) then
      i1 = i1 + 1
      write (slave2%name, '(2a, i0, a, i0)') trim(lord%name), '\', i, '#', i1      ! '
    else
      slave2_name = ''
      lmax = len(slave2%name) - 2
      do k = 1, slave2%n_lord
        lord2 => pointer_to_lord(slave2, k)
        if (lord2%n_lord > 0) lord2 => pointer_to_lord(lord2, 1)
        slave2_name = trim(slave2_name) // trim(lord2%name) // '\'     ! '
        if (len_trim(slave2_name) > lmax) exit
      enddo
      if (len_trim(slave2_name) > lmax) slave2_name = slave2_name(1:lmax) // '\'   ! '
      write (slave2%name, '(a, i0)') trim(slave2_name), i  
    endif
  enddo
enddo

! If this is a drift multipass whose multipass_slave elements are the result
! of splitting a drift with superposition then make sure that all split drift elements 
! of the lattice with the same base name have a name of the form "<base_name>#<n>" where
! <n> is an index from 1 for the first split drift.

ixb = index(lord%name, '#') - 1
if (lord%key == drift$ .and. ixb > 0) then
  ix_n = 0
  base_name = lord%name(1:ixb) 

  do i = 0, ubound(lat%branch, 1)
    branch => lat%branch(i)
    do j = 1, branch%n_ele_track
      ele => branch%ele(j)
      if (ele%key /= drift$) cycle
      ixb2 = index(ele%name, '#') - 1
      if (ixb2 /= ixb) cycle
      if (base_name(1:ixb) /= ele%name(1:ixb)) cycle
      ! super_slave drifts are temporary constructs that need to be ignored.
      ! This routine will be called later to correct the name of such elements.
      if (ele%slave_status == super_slave$) cycle 
      if (ele%slave_status == multipass_slave$) then
        call multipass_chain (ele, ix_pass, n_links)
        if (ix_pass /= 1) cycle  ! Only do renaming once
        lord2 => pointer_to_lord(ele, 1)
        ix_n = ix_n + 1
        write (lord2%name, '(2a, i0)') base_name(1:ixb), '#', ix_n
        do k = 1, lord2%n_slave
          slave => pointer_to_slave(lord2, k)
          write (slave%name, '(2a, i0, a, i0)') base_name(1:ixb), '#', ix_n, '\', k      !'
        enddo
      else
        ix_n = ix_n + 1
        write (ele%name, '(2a, i0)') base_name(1:ixb), '#', ix_n
      endif
    enddo
  enddo

endif

!

call control_bookkeeper (lat, lord)

end subroutine add_this_multipass

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine reallocate_bp_com_var()

implicit none

type (bp_var_struct) :: var_temp(size(bp_com%var))

integer n

!

var_temp = bp_com%var

deallocate (bp_com%var)

n = bp_com%ivar_tot+200
allocate (bp_com%var(n))

n = size(var_temp)
bp_com%var(1:n) = var_temp


end subroutine reallocate_bp_com_var

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_add_superimpose (branch, super_ele_in, pele, in_lat, plat)

use multipass_mod

implicit none

type (ele_struct) super_ele_in
type (ele_struct), save :: super_ele_saved, super_ele
type (ele_struct), pointer :: ref_ele, ele, slave, lord, super_ele_out, ele_at_s
type (ele_pointer_struct), allocatable :: eles(:)
type (parser_lat_struct), optional, target :: plat
type (parser_ele_struct) :: pele
type (parser_ele_struct), pointer :: ref_pele
type (multipass_all_info_struct) m_info
type (lat_struct), optional :: in_lat
type (lat_ele_loc_struct), allocatable :: m_slaves(:)
type (branch_struct), target :: branch
type (branch_struct), pointer :: ref_branch
type (lat_struct), pointer :: lat

integer ix, i, j, k, it, nic, nn, i_ele, ib, il
integer n_con, ix_branch, n_loc, ix_insert

character(40) name, ref_name
character(40), allocatable :: multi_name(:)
character(80) line

logical have_inserted, found, err_flag, err

! init

if (.not. bp_com%do_superimpose) return

call settable_dep_var_bookkeeping (super_ele_in)

call init_ele(super_ele_saved)
call init_ele(super_ele)

super_ele = super_ele_in
super_ele%logic = .false.
super_ele_saved = super_ele     ! in case super_ele_in changes
lat => branch%lat
pele%ix_ref_multipass = 0

! If no refrence element then superposition is simple.
! Remember, no reference element implies superposition on branch 0.

if (pele%ref_name == blank_name$) then
  if (branch%ix_branch /= 0) return
  if (bp_com%used_line_set_by_calling_routine) return
  call compute_super_lord_s (branch%ele(0), super_ele, pele, ix_insert)
  ele_at_s => pointer_to_element_at_s (branch, super_ele%s, .true., err_flag)
  if (err_flag) then
    call parser_error ('S-POSITION OUT OF BOUNDS FOR SUPERPOSITION OF: ' // super_ele%name)
    return
  endif

  if (ele_at_s%iyy == 0) then  ! If not in multipass region proceed as normal.
    call check_for_multipass_superimpose_problem (branch, super_ele, err_flag); if (err_flag) return
    call add_superimpose (lat, super_ele, 0, err_flag, save_null_drift = .true., &
                                        create_jumbo_slave = pele%create_jumbo_slave)
    if (err_flag) bp_com%error_flag = .true.
    return
  endif

  ! Must be in multipass region
  if (ele_at_s%slave_status == super_slave$) ele_at_s => pointer_to_lord(ele_at_s, 1)
  pele%ref_name = ele_at_s%name
  pele%ref_pt = anchor_end$
  pele%ele_pt = anchor_end$
  pele%offset = super_ele%s - ele_at_s%s
  pele%ix_ref_multipass = ele_at_s%iyy
endif

! Find all matches

call lat_ele_locator (pele%ref_name, lat, eles, n_loc, err)
if (err) then
  call parser_error ('MALFORMED SUPERIMPOSE REFERENCE ELEMENT NAME: ' // pele%ref_name, &
                     'FOR SUPERPOSITION OF: ' // super_ele_saved%name, pele = pele)
  return
endif

if (pele%ix_ref_multipass /= 0) then ! throw out elements that are the same physical element
  i = 1
  do
    if (i > n_loc) exit
    if (eles(i)%ele%iyy == pele%ix_ref_multipass) then
      i = i + 1
    else
      eles(i:n_loc-1) = eles(i+1:n_loc)  ! Remove
      n_loc = n_loc - 1
    endif
  enddo
endif

! If the ref element has key = null_ele$ and sub_key = drift$, this element was originally a drift that 
! has superimposed on. In this case, the information on where the original drift was
! has been lost so this is an error.

if (n_loc > 0) then
  ref_ele => eles(1)%ele
  if (ref_ele%key == null_ele$ .and. ref_ele%sub_key == drift$) then
    call parser_error ('DRIFT: ' // trim(ref_ele%name) // ' NO LONGER EXISTS FOR SUPERPOSITION OF ' // super_ele_saved%name, &
                       '[IT HAS DISAPPEARED DUE TO A PREVIOUS SUPERPOSITION]')
  endif
endif


! If the reference element is a group or overlay then this is fine as long as there
! is only one slave element.

if (n_loc == 0) then
  call lat_ele_locator (pele%ref_name, in_lat, eles, n_loc, err)
  if (n_loc == 0) return
  if (n_loc > 1) then
    call parser_error ('CONFUSED SUPERPOSITION. PLEASE CONTACT A BMAD MAINTAINER.')
    return
  endif

  ref_ele => eles(1)%ele
  if (ref_ele%key /= group$ .and. ref_ele%key /= overlay$) return
  ref_pele => plat%ele(ref_ele%ixx)
  ref_name = ref_pele%control(1)%name
  do i = 2, size(ref_pele%control)
    if (ref_pele%control(i)%name /= ref_name) then
      call parser_error ('SUPERPOSITION USING A GROUP OR OVERLAY AS A REFERENCE ELEMENT IS ONLY ALLOWED', &
                         'WHEN THE GROUP OR OVERLAY CONTROLS A SINGLE ELEMENT.', &
                         'FOR SUPERPOSITION OF: ' // super_ele_saved%name, pele = pele)
      return
    endif
  enddo 
 
  call lat_ele_locator (ref_name, lat, eles, n_loc, err)
  if (err .or. n_loc == 0) return
endif

! If there is a single ref element, and the superposition offset puts the super_ele into a multipass region,
! then shift the ref element to the multipass region.
! If there are multiple ref elements then the situation is too complicated and is considered an error.

if (n_loc == 1) then
  ref_ele => eles(1)%ele
  ref_branch => pointer_to_branch(ref_ele)
  if (ref_ele%iyy == 0 .and. ref_branch%ix_branch == branch%ix_branch) then
    call compute_super_lord_s (eles(1)%ele, super_ele, pele, ix_insert)
    ele_at_s => pointer_to_element_at_s (branch, super_ele%s_start, .true., err_flag)
    if (ele_at_s%slave_status == super_slave$) ele_at_s => pointer_to_lord(ele_at_s, 1)
    if (.not. err_flag) then
      if (ele_at_s%iyy /= 0) then  ! If in multipass region...
        pele%ref_name = ele_at_s%name
        pele%ref_pt = anchor_end$
        pele%ele_pt = anchor_end$
        pele%offset = super_ele%s - ele_at_s%s
        pele%ix_ref_multipass = ele_at_s%iyy
        call lat_ele_locator (ele_at_s%name, lat, eles, n_loc, err)
        i = 1 ! throw out elements that are the same physical element
        do
          if (i > n_loc) exit
          if (eles(i)%ele%iyy == pele%ix_ref_multipass) then
            i = i + 1
          else
            eles(i:n_loc-1) = eles(i+1:n_loc)  ! Remove
            n_loc = n_loc - 1
          endif
        enddo
      endif
    endif
  endif
endif

! Tag reference elements using %logic flag which is not otherwise used during parsing.

branch%ele(:)%logic = .false.  
do i = 1, n_loc
  eles(i)%ele%logic = .true. ! Tag reference element.
enddo

! Note: branch%n_ele_max will vary when an element is superimposed.

do 

  have_inserted = .false.

  i_ele = -1
  ele_loop: do 
    i_ele = i_ele + 1
    if (i_ele > branch%n_ele_max) exit

    ref_ele => branch%ele(i_ele)
     
    if (ref_ele%slave_status == super_slave$ .or. ref_ele%slave_status == multipass_slave$) then
      do il = 1, ref_ele%n_lord
        lord => pointer_to_lord(ref_ele, il)
        if (lord%slave_status == multipass_slave$) lord => pointer_to_lord(lord, 1)
        if (.not. lord%logic) cycle
        ref_ele => lord
        exit
      enddo
    endif

    if (ref_ele%key == group$) cycle
    if (ref_ele%key == girder$) cycle
    if (ref_ele%slave_status == super_slave$) cycle
    if (.not. ref_ele%logic) cycle

    ref_ele%logic = .false.  ! So only use this reference once

    ! If superimposing on a multipass_lord (only happens with bmad_parser2) 
    ! then the superposition must be done at all multipass locations.

    if (ref_ele%lord_status == multipass_lord$) then
      allocate (m_slaves(ref_ele%n_slave), multi_name(ref_ele%n_slave))
      allocate (ele_loc_com%branch(1))
      allocate (ele_loc_com%branch(1)%ele(ref_ele%n_slave))
      do i = 1, ref_ele%n_slave
        slave => pointer_to_slave(ref_ele, i)
        ele_loc_com%branch(1)%ele(i)%ix_ele    = slave%ix_ele
        ele_loc_com%branch(1)%ele(i)%ix_branch = slave%ix_branch
      enddo
      call string_trim(super_ele_saved%name, super_ele_saved%name, ix)
      super_ele%name = super_ele_saved%name(:ix)

      ! Put in the superposition at the multipass locations.
      ! Since elements get shuffled around, tag the superimposed elements 
      !     with "temp_name!" to identify them later.

      do i = 1, ref_ele%n_slave
        ele => pointer_to_ele (lat, ele_loc_com%branch(1)%ele(i))
        call compute_super_lord_s (ele, super_ele, pele, ix_insert)
        super_ele%iyy = ele%iyy   ! Multipass info
        call check_for_multipass_superimpose_problem (branch, super_ele, err_flag, ele); if (err_flag) return
        ! Don't need to save drifts since a multipass_lord drift already exists.
        call add_superimpose (lat, super_ele, ix_branch, err_flag, super_ele_out, &
               save_null_drift = .false., create_jumbo_slave = pele%create_jumbo_slave, ix_insert = ix_insert)
        if (err_flag) bp_com%error_flag = .true.
        super_ele_out%name = 'temp_name!'
      enddo

      ! Mark any super_lord drifts (along with any lords that control the super_lord) for 
      ! future deletion at the end of lattice parsing.

      call multipass_all_info (lat, m_info) ! Save multipass info for later.

      do i = lat%n_ele_track+1, lat%n_ele_max 
        ele => lat%ele(i)
        if (ele%key /= drift$) cycle
        if (ele%lord_status /= super_lord$) cycle 
        ele%key = -1 ! mark for deletion

        do j = 1, ele%n_lord
          lord => pointer_to_lord(ele, j)
          lord%key = -1  ! Mark lord for deletion
        enddo

        ! Need to remove super_lord/super_slave links otherwise the code below gets confused
        ! when it tries to connect the former super_slave drifts.

        do while (ele%n_slave /= 0)
          call remove_lord_slave_link (ele, pointer_to_slave(ele, 1))
        enddo
      enddo

      ! Reconnect drifts that were part of the multipass region.

      do i = 1, size(m_info%lord) ! Loop over multipass lords
        do j = 1, size(m_info%lord(i)%slave, 2)   ! loop over super_slaves
          slave => m_info%lord(i)%slave(1, j)%ele
          if (slave%key /= drift$) cycle
          if (slave%slave_status == multipass_slave$) cycle
          do k = 1, size(m_info%lord(i)%slave(:, j))   ! Loop over all passes
            ele => m_info%lord(i)%slave(k, j)%ele
            m_slaves(k) = ele_to_lat_loc (ele)  ! Make a list slave elements
            ib = index(ele%name, '\') ! '
            if (ib /= 0) ele%name = ele%name(1:ib-1) // ele%name(ib+2:)
          enddo
          call add_this_multipass (lat, m_slaves) ! And create a multipass lord
        enddo
      enddo

      ! If the super_lords have a single super_slave and the super_slave
      ! has only a single super_lord, the super_lords
      ! can be eliminated and the created multipass_lord can control the
      ! super_slaves directly.

      j = 0

      do i = 1, branch%n_ele_track
        if (branch%ele(i)%name == 'temp_name!') then
          branch%ele(i)%name = super_ele_saved%name
          j = j + 1
          m_slaves(j) = ele_to_lat_loc (branch%ele(i))
        endif
      enddo

      do i = lat%n_ele_track+1, lat%n_ele_max
        if (lat%ele(i)%name == 'temp_name!') then
          lat%ele(i)%name = super_ele_saved%name
          j = j + 1
          m_slaves(j) = ele_to_lat_loc (lat%ele(i))
        endif
      enddo

      ele => pointer_to_ele (lat, m_slaves(1))
      if (ele%lord_status == super_lord$ .and. ele%n_slave == 1) then
        slave => pointer_to_slave(ele, 1)
        if (slave%n_lord == 1) then
          do i = 1, size(m_slaves)
            ele => pointer_to_ele (lat, m_slaves(i))
            ele%key = -1 ! Mark for deletion
            ele => pointer_to_slave(ele, 1)
            ele%name = super_ele_saved%name
            m_slaves(i) = ele_to_lat_loc (ele)
          enddo
        endif
      endif

      ! Remove eles marked for deletion. But first shift m_slaves list

      do i = 1, size(m_slaves)
        ele => pointer_to_ele (lat, m_slaves(i))
        m_slaves(i)%ix_ele = m_slaves(i)%ix_ele - count(lat%branch(ele%ix_branch)%ele(1:ele%ix_ele)%key == -1)
      enddo

      call remove_eles_from_lat (lat, .false.)

      ! Add a multipass_lord to control the created super_lords.

      call add_this_multipass (lat, m_slaves, super_ele_saved) 

      call deallocate_multipass_all_info_struct (m_info)
      deallocate (m_slaves, multi_name)
      deallocate (ele_loc_com%branch)

    !-----------------------
    ! Else not superimposing on a multipass_lord ...
    ! [Note: Only will superimpose on a multipass_lord in bmad_parser2.]

    else
      call compute_super_lord_s (ref_ele, super_ele, pele, ix_insert)
      super_ele%iyy = ref_ele%iyy   ! Multipass info
      call check_for_multipass_superimpose_problem (branch, super_ele, err_flag, ref_ele); if (err_flag) return
      call string_trim(super_ele_saved%name, super_ele_saved%name, ix)
      super_ele%name = super_ele_saved%name(:ix)            
      call add_superimpose (lat, super_ele, branch%ix_branch, err_flag, super_ele_out, &
         save_null_drift = .true., create_jumbo_slave = pele%create_jumbo_slave, ix_insert = ix_insert)
      if (err_flag) bp_com%error_flag = .true.
      call control_bookkeeper (lat, super_ele_out)
    endif

    !---------------------

    call s_calc (lat)

    have_inserted = .true.   

  enddo ele_loop

  if (.not. have_inserted) exit

enddo

end subroutine parser_add_superimpose

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_check_superimpose_valid_ref (ele, lat, pele, in_lat)

implicit none

type (ele_struct) ele
type (lat_struct) lat
type (parser_ele_struct) :: pele
type (lat_struct), optional :: in_lat
type (ele_pointer_struct), allocatable :: eles(:)

integer n_loc
logical found, err

! If a superimpose reference element does not appear in the list of elements in the lattice file and
! does not appear in the final lattice, this is considered a mispelling and is considered an error.
! Note: Something like "branch_name>>ele_name" may show up in the finished lattice but not the list of
! elements in the lattice file.

if (.not. bp_com%do_superimpose) return
if (pele%ref_name == blank_name$) return

call lat_ele_locator (pele%ref_name, lat, eles, n_loc, err)
if (err) return    ! this error already handled by parser_add_superimpose
if (n_loc /= 0) return

if (present(in_lat)) then
  call lat_ele_locator (pele%ref_name, in_lat, eles, n_loc, err)
  if (n_loc /= 0) return
endif

call parser_error ('NO MATCH FOR REFERENCE ELEMENT: ' //  pele%ref_name, &      
                   'FOR SUPERPOSITION OF: ' // ele%name, pele = pele)

end subroutine parser_check_superimpose_valid_ref

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine compute_super_lord_s (ref_ele, super_ele, pele, ix_insert)

implicit none

type (lat_struct), target :: lat
type (ele_struct) ref_ele, super_ele
type (ele_struct), pointer :: slave
type (parser_ele_struct) pele
type (branch_struct), pointer :: branch

integer i, ix, ix_insert

real(rp) s_ref_begin, s_ref_end

! Find the reference point on the element being superimposed.

ix_insert = -1
super_ele%s = pele%offset

if (pele%ele_pt == anchor_beginning$) then
  super_ele%s = super_ele%s + super_ele%value(l$)
elseif (pele%ele_pt == anchor_center$ .or. pele%ele_pt == not_set$) then
  super_ele%s = super_ele%s + super_ele%value(l$) / 2
elseif (pele%ele_pt /= anchor_end$) then
  call parser_error ('ERROR IN COMPUTE_SUPER_LORD_S: CONTROL #1 INTERNAL ERROR!')
  if (global_com%exit_on_error) call err_exit
endif

! Find the refernce point in the lattice.

select case (ref_ele%key)
case (overlay$, girder$)
  s_ref_begin = 1d10
  s_ref_end = 0
  do i = 1, ref_ele%n_slave
    slave => pointer_to_slave(ref_ele, i)
    s_ref_begin = min(s_ref_begin, slave%s_start)
    s_ref_end = max(s_ref_end, slave%s)
  enddo
case (group$)
  call parser_error ('SUPERPOSING: ' // super_ele%name, 'UPON GROUP' // pele%ref_name)
  return
case default
  s_ref_begin = ref_ele%s_start
  s_ref_end = ref_ele%s
end select

! Now compute the s position at the end of the element and put it in ele%s.

if (pele%ref_pt == anchor_beginning$) then
  super_ele%s = super_ele%s + s_ref_begin
elseif (pele%ref_pt == anchor_center$ .or. pele%ref_pt == not_set$) then
  super_ele%s = super_ele%s + (s_ref_begin + s_ref_end) / 2
elseif (pele%ref_pt == anchor_end$) then
  super_ele%s = super_ele%s + s_ref_end
else
  call parser_error ('ERROR IN COMPUTE_SUPER_LORD_S: CONTROL #2 INTERNAL ERROR!')
  if (global_com%exit_on_error) call err_exit
endif

! For circular lattices a superimpose can wrap around the beginning or 
! the end of the lattice.

branch => ref_ele%branch

if (branch%param%geometry == closed$) then
  if (super_ele%s > branch%ele(branch%n_ele_track)%s) then
    super_ele%s = super_ele%s - branch%param%total_length
  elseif (super_ele%s < 0) then
    super_ele%s = super_ele%s + branch%param%total_length
  endif
endif

super_ele%s_start = super_ele%s - super_ele%value(l$)
if (branch%param%geometry == closed$) then
  if (super_ele%s_start < 0) then
    super_ele%s_start = super_ele%s_start + branch%param%total_length
  endif
endif

! Compute ix_insert which is used for positioning zero length super_lords in case
! They are next to an element which is also zero length.

if (ref_ele%value(l$) == 0 .and. super_ele%value(l$) == 0 .and. pele%offset == 0) then
  if (ref_ele%ix_ele == 0) then  ! Automatically must be at downstream end.
    ix_insert = 1
  elseif (pele%ref_pt == anchor_beginning$) then
    ix_insert = ref_ele%ix_ele  
  elseif (pele%ref_pt == anchor_end$) then
    ix_insert = ref_ele%ix_ele + 1 
  endif
endif

end subroutine compute_super_lord_s

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine check_for_multipass_superimpose_problem (branch, super_ele, err_flag, ref_ele)
!
! Subroutine to check if there is a problem superimposing an element when there is multipass.
! In particular will check that:
!   1) If the ref_ele is part of a multipass region then super_ele must be superimposed
!      within the region.
! Or:
!   2) If the ref_ele is not part of a multipass region then super_ele must also not
!      be part of a multipass region.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine check_for_multipass_superimpose_problem (branch, super_ele, err_flag, ref_ele)

implicit none

type (ele_struct) super_ele
type (ele_struct), optional :: ref_ele
type (ele_struct), pointer :: ele1, ele2
type (branch_struct) :: branch
real(rp) eps
logical err_flag
integer ix1, ix2


!

eps = bmad_com%significant_length

ele1 => pointer_to_element_at_s (branch, super_ele%s_start + eps, .true., err_flag)
if (err_flag) then
  call parser_error ('BAD SUPERIMPOSE OF: ' // super_ele%name, 'UPSTREAM ELEMENT EDGE OUT OF BOUNDS.')
  return
endif
if (ele1%slave_status == super_slave$) ele1 => pointer_to_lord(ele1, 1)

ele2 => pointer_to_element_at_s (branch, super_ele%s - eps, .false., err_flag)
if (err_flag) then
  call parser_error ('BAD SUPERIMPOSE OF: ' // super_ele%name, 'DOWNSTREAM ELEMENT EDGE OUT OF BOUNDS.')
  return
endif
if (ele2%slave_status == super_slave$) ele2 => pointer_to_lord(ele2, 1)

if (present(ref_ele)) then
  if (ref_ele%iyy /= 0) then     ! Ref element in multipass region
    if (abs(super_ele%value(l$)) < eps .and. (ele1%iyy /= 0 .or. ele2%iyy /= 0)) return ! At multipass edge is OK
    if (ele1%iyy == 0 .or. ele2%iyy == 0) then
      call parser_error ('SUPERIMPOSE OF: ' // super_ele%name, &
           'USES MULTIPASS REFERENCE ELEMENT BUT OFFSET PLACES IT OUT OF THE MULTIPASS REGION!')
      return
    endif

  else
    if (abs(super_ele%value(l$)) < eps .and. (ele1%iyy == 0 .or. ele2%iyy == 0)) return  ! At multipass edge is OK
    if (ele1%iyy /= 0 .or. ele2%iyy /= 0) then
      call parser_error ('SUPERIMPOSE OF: ' // super_ele%name, &
                         'USES NON-MULTIPASS REFERENCE ELEMENT BUT OFFSET PLACES IT IN A MULTIPASS REGION!')
      return
    endif
  endif
endif

err_flag = .false.

end subroutine check_for_multipass_superimpose_problem 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_sequence_args (seq_name, arg_list, delim, err_flag)
!
! Subroutine to get the argument list for a replacement_line or a list.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine get_sequence_args (seq_name, arg_list, delim, err_flag)

implicit none

integer n_arg, ix_word

character(*), pointer :: arg_list(:)
character(*) seq_name
character(1) delim
character(40) name(20), word

logical delim_found, err_flag

!

n_arg = 0
err_flag = .true.

do
  call get_next_word (word, ix_word, '(): =,', delim, delim_found, .true.)
  if (ix_word == 0 .or. delim == '( :=') then
    call parser_error ('BAD ARGUMENT LIST FOR: ', seq_name)
    return
  endif
  n_arg = n_arg + 1
  name(n_arg) = word
  if (delim == ')') exit
enddo

err_flag = .false.
allocate (arg_list(n_arg))
arg_list = name(1:n_arg)

end subroutine get_sequence_args

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine seq_expand1 (sequence, iseq_tot, lat, top_level)
!
! Subroutine to expand a sequence.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

recursive subroutine seq_expand1 (sequence, iseq_tot, lat, top_level)

implicit none

type (seq_struct), target :: sequence(:)
type (seq_struct), pointer :: seq, sub_seq
type (seq_ele_struct), pointer :: s_ele(:)
type (seq_ele_struct), pointer :: s_ele2(:)
type (seq_ele_struct), pointer :: this_ele
type (lat_struct) lat

integer ix_ele, iseq_tot, ix_word, ix, i, n, ios, rcount
integer, save :: ix_internal = 0

character(40) word
character(1) delim, c_delim
character(40) str
character(n_parse_line) parse_line_saved

logical delim_found, replacement_line_here, c_delim_found
logical err_flag, top_level

! init

allocate (s_ele(ubound(lat%ele, 1)))

! save info on what file we are parsing for error messages.

seq => sequence(iseq_tot)
seq%file_name = bp_com%current_file%full_name 
seq%ix_line = bp_com%current_file%i_line

! first thing should be a "("

call get_next_word(word, ix_word, ':=(),', delim, delim_found, .true.)

if (delim /= '(') call parser_error  &
      ('EXPECTING "(", GOT: ' // delim, 'FOR LINE: ' // seq%name)
if (ix_word /= 0)  call parser_error  &
      ('EXTRANEOUS STUFF BEFORE "(", GOT: ' // word,  &
      'FOR LINE: ' // seq%name)

! now parse list proper

ix_ele = 1
this_ele => s_ele(ix_ele)

do 

  call get_next_word (word, ix_word, ':=(,)[]*', delim, delim_found, .true.)

  this_ele%rep_count = 1

  if (word(1:1) == '-') then
    this_ele%ele_order_reflect = .true.
    word = word(2:)
    ix_word = ix_word - 1
  endif

  if (word(1:1) == '-') then
    this_ele%ele_orientation = -1
    word = word(2:)
    ix_word = ix_word - 1
  endif

  if (delim == '*') then    ! E.g. '-3*(A,B)'
    ! Evaluate the rep count.
    read (word, *, iostat = ios) rcount
    if (ix_word == 0 .or. ios /= 0) then
      call parser_error ('MALFORMED REPETION COUNT FOUND IN SEQUENCE: ' // seq%name)
      return
    endif
    this_ele%rep_count = rcount
    call get_next_word (word, ix_word, ':=(,)[]*', delim, delim_found, .true.)
  endif

  this_ele%name = word
  if (word /= ' ') then
    if (.not. verify_valid_name (word, ix_word)) return
  endif

  ! Check for line tag

  if (delim == '[') then
    call get_next_word (this_ele%tag, ix_word, '[]:=(,)', delim, delim_found, .true.)
    if (delim /= ']') then
      call parser_error ('NO MATCHING "]" FOUND FOR OPENING "[" IN SEQUENCE: ' // seq%name)
      return
    endif
    call get_next_word (word, ix_word, '[]:=(,)', delim, delim_found, .true.)
    if (ix_word > 0) then
      call parser_error ('ILLEGAL CHARACTERS AFTER CLOSING "]" FOUND IN SEQUENCE: ' // seq%name)
      return
    endif   
  endif

  ! Check for a subline or replacement line.
  ! If there is one then save as an internal sequence.

  replacement_line_here = .false.

  if (delim == '(') then ! subline or replacement line

    ! if a subline...
    if (this_ele%name == '') then  
      ix_internal = ix_internal + 1
      write (str, '(a, i3.3)') '#Internal', ix_internal   ! unique name 
      this_ele%name = str
      iseq_tot = iseq_tot + 1
      sub_seq => sequence(iseq_tot) 
      sub_seq%name = str
      sub_seq%type = seq%type
      sub_seq%multipass = seq%multipass
      if (sub_seq%type == replacement_line$) then
        ix = size (seq%dummy_arg)
        allocate (sub_seq%dummy_arg(ix), &
              sub_seq%corresponding_actual_arg(ix), this_ele%actual_arg(ix))
        sub_seq%dummy_arg = seq%dummy_arg
        this_ele%actual_arg = seq%dummy_arg
      endif
      bp_com%parse_line = '(' // bp_com%parse_line
      call seq_expand1 (sequence, iseq_tot, lat, .false.)

    ! else this is a replacement line
    else    
      replacement_line_here = .true.
      call get_sequence_args (this_ele%name, this_ele%actual_arg, delim, err_flag)
      if (err_flag) return
    endif

    call get_next_word(word, ix_word, ':=(),', delim, delim_found, .true.)
    if (word /= ' ') call parser_error &
              ('NO COMMA AFTER SUBLINE OR REPLACEMENT LINE. FOUND: ' // &
               word, 'IN THE SEQUENCE: ' // seq%name)
  endif

  if (this_ele%name == ' ') call parser_error &
            ('SUB-ELEMENT NAME IS BLANK FOR LINE/LIST: ' // seq%name)

  ! if a replacement line then look for element in argument list

  this_ele%ix_arg = 0
  if (seq%type == replacement_line$) then
    do i = 1, size(seq%dummy_arg)
      if (seq%dummy_arg(i) == this_ele%name) then
        this_ele%ix_arg = i
        exit
      endif
    enddo
  endif

! 

  n = size(s_ele)
  ix_ele = ix_ele + 1

  if (ix_ele > n) then
    allocate (s_ele2(n))      
    s_ele2 = s_ele(1:n)
    deallocate (s_ele)
    allocate (s_ele(n+1000))
    s_ele(1:n) = s_ele2
    deallocate(s_ele2)
  endif

  this_ele => s_ele(ix_ele)

  if (delim == ')') exit

  if (delim /= ',') then
    call parser_error ('EXPECTING "," GOT: ' // delim, 'FOR LINE: ' // seq%name)
    exit
  endif
         
enddo

! make sure there is nothing else if at top level

if (top_level) then
  call get_next_word(word, ix_word, ':=() ', delim, delim_found, .true.)
  if (delim_found .or. ix_word /= 0) call parser_error  &
        ('EXTRA CHARACTERS AFTER CLOSING ")"',  'FOR LINE: ' // seq%name)
endif

! transfer

ix_ele = ix_ele - 1
allocate (seq%ele(ix_ele))

do i = 1, ix_ele
  seq%ele(i) = s_ele(i)
enddo

deallocate (s_ele)

end subroutine seq_expand1

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine allocate_plat (plat, n_ele_max) 
!
! Subroutine to allocate allocatable array sizes.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

Subroutine allocate_plat (plat, n_ele_max)

implicit none

type (parser_lat_struct) plat, temp_plat

integer i, n_now, n_ele_max

! assume all the arrays have the same size

if (allocated(plat%ele)) then
  n_now = ubound(plat%ele, 1)
  call move_alloc (plat%ele, temp_plat%ele)
  allocate (plat%ele(0:n_ele_max))
  plat%ele(0:n_now) = temp_plat%ele
  deallocate (temp_plat%ele)

else
  allocate (plat%ele(0:n_ele_max))
  n_now = -1
endif

! %ixx is used as a pointer from the in_lat%ele array to the plat%ele array

do i = n_now+1, ubound(plat%ele, 1)
  plat%ele(i)%ele_name = ''
  plat%ele(i)%ref_name = blank_name$
  plat%ele(i)%ref_pt  = not_set$
  plat%ele(i)%ele_pt  = not_set$
  plat%ele(i)%offset  = 0
  plat%ele(i)%create_jumbo_slave = .false.
enddo

end subroutine allocate_plat

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_add_lord (in_lat, n_ele_max, plat, lat)
!
! Subroutine to add overlay, group, and girder lords to the lattice.
! For overlays and groups: If multiple elements have the same name then 
! use all of them.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_add_lord (in_lat, n_ele_max, plat, lat)

implicit none

type multi_ele_pointer_struct
  type (ele_pointer_struct), allocatable :: eles(:)
  integer n_loc
end type  

type (lat_struct), target :: in_lat, lat
type (ele_struct), pointer :: lord, lord2, slave, ele, g_lord, g_slave0, g_slave1
type (parser_lat_struct), target :: plat
type (parser_ele_struct), pointer :: pele
type (control_struct), allocatable, target :: cs(:)
type (branch_struct), pointer :: branch
type (multi_ele_pointer_struct), allocatable :: m_eles(:)
type (ele_pointer_struct), allocatable :: in_eles(:)
type (parser_controller_struct), pointer :: pc
type (control_struct), pointer :: con

integer i, n_in, ic, ig, k, ix, ib, ie_start, n_list, ns, ixs, ii, ix_end, n_ele_max
integer n_slave, nn, n_loc, n_names, n_in_loc, ix1_slave
integer ix_lord, k_slave, ix_ele_now, ix_super_lord_end

character(60) err_str
character(40) input_slave_name, attrib_name, missing_slave_name

logical err, slave_not_in_lat, created_girder_lord, err_flag, matched_to_drift, have_ignored_a_drift

! loop over lord elements

main_loop: do n_in = 1, n_ele_max

  lord => in_lat%ele(n_in)  ! next lord to add
  pele => plat%ele(lord%ixx)
  
  !-----------------------------------------------------
  ! overlay and groups

  select case (lord%key)
  case (overlay$, group$)
 
    ! If all slave elements are defined in in_lat, but are not present in lat, then
    ! this lord can be ignored.

    slave_not_in_lat = .false.
    n_in_loc = 1
    n_slave = 0

    if (allocated(m_eles)) deallocate (m_eles)
    allocate (m_eles(size(pele%control)))

    do i = 1, size(pele%control)
      pc => pele%control(i)
      call lat_ele_locator (pc%name, lat, m_eles(i)%eles, m_eles(i)%n_loc, err)
      n_loc = m_eles(i)%n_loc
      n_slave = n_slave + n_loc

      if (n_loc == 0) then
        slave_not_in_lat = .true.
        missing_slave_name = pc%name
        call lat_ele_locator (pc%name, in_lat, in_eles, n_in_loc, err)
      endif

      if ((n_loc == 0 .and. n_slave > 0) .or. (n_loc > 0 .and. slave_not_in_lat) .or. &
          (n_loc == 0 .and. n_in_loc == 0)) then
        call parser_error ('CANNOT FIND SLAVE FOR: ' // lord%name, &
                           'CANNOT FIND: '// missing_slave_name, pele = pele)
        cycle main_loop
      endif
    enddo

    if (n_slave == 0) cycle main_loop

    ! Create the lord(s)

    if (is_true(lord%value(gang$))) then
      call make_this_overlay_group_lord(0, err_flag)
    else
      do i = 2, size(pele%control)
        if (m_eles(1)%n_loc /= m_eles(i)%n_loc) then
          call parser_error ('IN OVERLAY OR GROUP ELEMENT: ' // lord%name, &
                    'WITH GANG = FALSE NEED ALL SLAVES WITH A GIVEN NAME TO HAVE THE SAME NUMBER OF', &
                    'ELEMENTS IN THE LATTICE. BUT ' // trim(pele%control(1)%name) // ' HAS \i0\ ELEMENTS', &
                    'WHILE ' // trim(pele%control(I)%name) // ' HAS \i0\ ELEMENTS', &
                    i_array = [m_eles(1)%n_loc, m_eles(I)%n_loc])
          cycle main_loop
        endif
      enddo

      do nn = 1, m_eles(1)%n_loc
        call make_this_overlay_group_lord(nn, err_flag)
        if (err_flag) exit
      enddo
    endif

  !-----------------------------------------------------
  ! girder
  ! Create an girder element for each lattice section that matches the slave list names.
  ! If no lattice element names match any names on the girder slave list then assume 
  ! this girder is for a different lattice and ignore the girder. If some do match and 
  ! some don't then flag this as an error.

  case (girder$)

    if (pele%is_range .and. size(pele%control) /= 2) then
      call parser_error ('GIRDER HAS BAD "ELE_START:ELE_END" RANGE CONSTRUCT. ' // ele%name)
      cycle
    endif

    ! Loop over all elements in the lattice.

    if (allocated(cs)) deallocate(cs)
    allocate (cs(size(pele%control)))

    created_girder_lord = .false.

    branch_loop: do ib = 0, ubound(lat%branch, 1)
      branch => lat%branch(ib)
      ie_start = 1
      ix1_slave = -1

      ! Loop over all possible first elements

      ele_loop: do ie_start = 1, branch%n_ele_track


        ! Loop over girder slave list and see if this section matches.

        ix_ele_now = ie_start
        ix_super_lord_end = -1   ! Not in a super_lord
        ixs = 1                  ! Index of girder slave element we are looking for.
        n_slave = 0              ! Number of actual slaves found.
        matched_to_drift = .false.
        have_ignored_a_drift = .false.

        if (size(pele%control) == 0) then
          call parser_error ('GIRDER DOES NOT HAVE ANY ELEMENTS TO SUPPORT: ' // lord%name)
          cycle main_loop
        endif

        ! loop over all girder slaves and see if lattice eles match slaves.

        slave_loop: do
          if (n_slave > 0 .and. cs(1)%slave%ix_ele == ix1_slave) cycle ele_loop  ! Can happen with superposition
          if (ixs > size(pele%control)) exit

          ! Wrap around the origin if needed.
          if (ix_ele_now > branch%n_ele_track) then
            ix_ele_now = ix_ele_now - branch%n_ele_track
          endif

          input_slave_name = pele%control(ixs)%name

          ele => pointer_to_ele (lat, ix_ele_now, ib)

          if (girder_match_slave_element(ele, ele, n_slave, ixs, ix_ele_now)) cycle

          ! Here if no match. 
          ! If not previously in a super lord then there is no overall match to the slave list.

          if (ele%slave_status == super_slave$ .and. ix_super_lord_end > -1) then
            if (ix_ele_now == ix_super_lord_end) ix_super_lord_end = -1
            ix_ele_now = ix_ele_now + 1
            cycle
          endif

          ! No match to the slave list. If a marker or drift then ignore except 
          ! if this is the first slave in which case there is no match.

          if ((ele%key == marker$ .or. ele%key == drift$) .and. ixs > 1) then
            if (ix_ele_now == ix_super_lord_end) ix_super_lord_end = -1
            ix_ele_now = ix_ele_now + 1
            if (ele%key == drift$) have_ignored_a_drift = .true.
            cycle
          endif

          ! Match failed. Start again 

          cycle ele_loop

        enddo slave_loop

        if (matched_to_drift .and. have_ignored_a_drift .and. .not. pele%is_range) cycle  ! matching rules violated.

        ! create the girder element

        if (n_slave == 0) then
          call parser_error ('LIST OF GIRDER SLAVES IN LATTICE FILE DOES NOT INCLUDE A NON-DRIFT ELEMENT: ' // &
                              lord%name, level = s_warn$)
          cycle main_loop
        endif

        call new_control (lat, ix_lord)
        call create_girder (lat, ix_lord, cs(1:n_slave), lord)
        created_girder_lord = .true.
        ix1_slave = cs(1)%slave%ix_ele

      enddo ele_loop

    enddo branch_loop

    if (.not. created_girder_lord) then
      call parser_error ('FAILED TO FIND REGION IN LATTICE FOR CREATING GIRDER: ' // &
                          lord%name, level = s_warn$)
    endif

  end select

enddo main_loop

!-------------------------------------------------------------------------
contains

recursive function girder_match_slave_element (ele, slave, n_slave, ixs, ix_ele_now) result (is_matched)

type (ele_struct), target :: ele, slave
type (ele_struct), pointer :: lord, slave1, slave2
type (control_struct), allocatable :: cs_temp(:)
integer n_slave, ixs, ix_ele_now
integer ii, ls, n, ie
logical is_matched, add_slave

! Try to match

is_matched = match_wild(ele%name, input_slave_name)

if (is_matched .or. (pele%is_range .and. ixs == 2 .and. ele%slave_status == not_a_child$)) then

  if (is_matched) ixs = ixs + 1  ! Next element in list

  if (ele%key == drift$) then
    if (is_matched) matched_to_drift = .true.

  else ! If not a drift then ele will be a girder_slave
    ! First check for duplicates. For example, element with superimposed 
    ! marker will be duplicated unless a check is made.
    add_slave = .true.
    do ie = 1, n_slave
      if (ele%ix_ele == cs(ie)%slave%ix_ele .and. ele%ix_branch == cs(ie)%slave%ix_branch) add_slave = .false.
    enddo

    if (add_slave) then
      n_slave = n_slave + 1
      if (size(cs) < n_slave) then  ! Can happen if there is a range.
        n = size(cs)
        call move_alloc(cs, cs_temp)
        allocate (cs(2*n))
        cs(1:n) = cs_temp
      endif
      cs(n_slave)%slave = lat_ele_loc_struct(ele%ix_ele, ele%ix_branch)
    endif
  endif

  is_matched = .true.

  ! If a super_lord the logic here is complicated by the fact that 
  ! elements can, for example, be completely contained within another element.
  ! Note: slave can be a super_lord if ele is a multipass_lord

  if (ele%lord_status == super_lord$) then
    slave2 => pointer_to_slave(ele, ele%n_slave)
    ix_super_lord_end = ix_far_index(ix_ele_now-1, ix_super_lord_end, slave2%ix_ele)
  endif

  if (slave%lord_status == super_lord$) then
    slave2 => pointer_to_slave(slave, slave%n_slave)
    ix_super_lord_end = ix_far_index(ix_ele_now-1, ix_super_lord_end, slave2%ix_ele)
  endif

  ! If in a super_slave region then need to recheck the current slave against the next girder slave name.

  if (ix_super_lord_end == -1 .or. pele%is_range) ix_ele_now = ix_ele_now + 1

  ! If match to a girder then move pointers to element after last girder slave

  if (ele%key == girder$) then
    call find_element_ends (ele, slave1, slave2)
    ix_ele_now = slave2%ix_ele + 1
    ix_super_lord_end = -1
  endif

  return
endif

! Since ele does not match, look for match at a lord of this element

do ii = 1, ele%n_lord
  lord => pointer_to_lord (ele, ii)
  ls = lord%lord_status
  if (ls /= super_lord$ .and. ls /= multipass_lord$ .and. ls /= girder_lord$) cycle
  is_matched = girder_match_slave_element(lord, ele, n_slave, ixs, ix_ele_now)
  if (is_matched) return
enddo

end function girder_match_slave_element

!-------------------------------------------------------------------------
! contains

! Function to return the index that is farthest (reached last) when moving
! from ix_now in a positive direction. Tricky part is if there is wrap around.
! An index of -1 means that the corresponding point does not exist.

function ix_far_index (ix_now, ix1, ix2) result (ix_far)

implicit none

integer ix_now, ix1, ix2, ix_far

!

if (ix1 < 0) then      ! Point 1 does not exist so point 2 wins by default
  ix_far = ix2
elseif (ix2 < 0) then  ! Point 2 does not exist so point 1 wins by default
  ix_far = ix1
elseif (ix1 < ix_now .and. ix2 < ix_now) then  ! both wrapped case
  ix_far = max(ix1, ix2)
elseif (ix1 > ix_now .and. ix2 > ix_now) then ! No wrap ("normal") case
  ix_far = max(ix1, ix2)
else                      ! One is wrapped but not the other case
  ix_far = min(ix1, ix2)
endif

end function ix_far_index

!-------------------------------------------------------------------------
! contains

subroutine make_this_overlay_group_lord (ix_pick, err_flag)

integer ix_pick
logical err_flag

!

err_flag = .true.

call new_control (lat, ix_lord)  ! get index in lat where lord goes
lat%ele(ix_lord) = lord

if (allocated(cs)) then
  if (size(cs) < n_slave) deallocate(cs)
endif
if (.not. allocated(cs)) allocate (cs(n_slave))

! Slave setup

n_slave = 0 ! number of slaves found
do i = 1, size(pele%control)

  pc => pele%control(i)

  ! There might be more than 1 element with same name. 
  ! Loop over all elements whose name matches name.
  ! Put the info into the cs structure.

  do k = 1, m_eles(i)%n_loc
    if (ix_pick /= 0 .and. k /= ix_pick) cycle
    slave => pointer_to_ele (lat, m_eles(i)%eles(k)%loc)
    n_slave = n_slave + 1
    call reallocate_expression_stack (cs(n_slave)%stack, pc%n_stk)
    cs(n_slave)%stack = pc%stack(1:pc%n_stk)
    cs(n_slave)%slave = lat_ele_loc_struct(slave%ix_ele, slave%ix_branch)
    cs(n_slave)%lord%ix_ele = -1             ! dummy value
    attrib_name = pc%attrib_name
    if (attrib_name == blank_name$) attrib_name = pele%default_attrib
    ix = attribute_index(slave, attrib_name)
    ! If attribute not found it may be a special attribute like accordion_edge$.
    ! A special attribute will have ix > num_ele_attrib$
    if (ix < 1 .and. lord%key == group$) then
      ix = attribute_index(lord, attrib_name)
      if (ix <= num_ele_attrib$) ix = 0  ! Mark as not valid
    endif
    cs(n_slave)%ix_attrib = ix
    if (ix < 1) then
      call parser_error ('IN OVERLAY OR GROUP ELEMENT: ' // lord%name, &
                    'ATTRIBUTE: ' // attrib_name, &
                    'IS NOT A VALID ATTRIBUTE OF: ' // slave%name, &
                    pele = pele)
      return
    endif
  enddo

enddo

! If the lord has no slaves then discard it

if (n_slave == 0) then
  lat%n_ele_max = lat%n_ele_max - 1 ! Undo new_control call
  return
endif

! create the lord

select case (lord%key)
case (overlay$)
  call create_overlay (lat%ele(ix_lord), cs(1:n_slave), err)
case (group$)
  call create_group (lat%ele(ix_lord), cs(1:n_slave), err)
end select
if (err) then
  call parser_error ('MALFORMED OVERLAY OR GROUP: ' // lord%name, &
                     'IS TRYING TO CONTROL AN ATTRIBUTE THAT IS NOT FREE TO VARY!', &
                       pele = pele)
  return
endif

lat%ele(ix_lord)%value(gang$) = lord%value(gang$)

! Evaluate any variable values.

lord2 => lat%ele(ix_lord)
do k = lord2%ix1_slave, lord2%ix1_slave+lord2%n_slave-1
  con => lat%control(k)
  do ic = 1, size(con%stack)
    select case (con%stack(ic)%type)
    case (ran$, ran_gauss$)
      call parser_error ('RANDOM NUMBER FUNCITON MAY NOT BE USED WITH AN OVERLAY OR GROUP', &
                         'FOR ELEMENT: ' // lord%name)
    case (variable$)
      call word_to_value (con%stack(ic)%name, lat, con%stack(ic)%value)
      ! Variables in the arithmetic expression are immediately evaluated and never reevaluated.
      ! If the variable is an element attribute (looks like: "ele_name[attrib_name]") then this may
      ! be confusing if the attribute value changes later. To avoid some (but not all) confusion, 
      ! turn the variable into a numeric$ so the output from the type_ele routine looks "sane".
      if (index(con%stack(ic)%name, '[') /= 0) then
        con%stack(ic)%type = numeric$
        con%stack(ic)%name = ''
      endif
    end select
  enddo
enddo

err_flag = .false.

end subroutine make_this_overlay_group_lord

end subroutine parser_add_lord

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine settable_dep_var_bookkeeping (ele)
!
! Subroutine to initialize dependent variables in an element.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine settable_dep_var_bookkeeping (ele)

use random_mod

implicit none

type (ele_struct),  target :: ele
type (branch_struct), pointer :: branch

real(rp) angle, rr, v_inv_mat(4,4), eta_vec(4)

integer n
logical kick_set, length_set, set_done, err_flag
logical b_field_set, g_set

! Wall3d init.

if (associated(ele%wall3d)) then
  do n = 1, size(ele%wall3d)
    call wall3d_initializer (ele%wall3d(n), err_flag)
    if (err_flag) then
      call parser_error ('WALL INIT ERROR FOR ELEMENT: ' // ele%name)
      return
    endif
  enddo
endif

! Aperture init

if (ele%aperture_type == auto_aperture$) then
  call aperture_bookkeeper (ele)
endif

!

kick_set = (ele%value(hkick$) /= 0) .or. (ele%value(vkick$) /= 0)

select case (ele%key)

case (beginning_ele$)

  if (ele%a%beta /= 0) ele%a%gamma = (1 + ele%a%alpha**2) / ele%a%beta
  if (ele%b%beta /= 0) ele%b%gamma = (1 + ele%b%alpha**2) / ele%b%beta

  call make_v_mats (ele, v_inv_mat = v_inv_mat)
  eta_vec = matmul (v_inv_mat, [ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap])

  ele%a%eta  = eta_vec(1)
  ele%a%etap = eta_vec(2)
  ele%b%eta  = eta_vec(3)
  ele%b%etap = eta_vec(4)

! Convert rbends to sbends and evaluate G if needed.
! Needed is the length and either: angle, G, or rho.

case (sbend$, rbend$, sad_mult$) 

  b_field_set = (ele%value(b_field$) /= 0 .or. ele%value(b_field_err$) /= 0)
  g_set = (ele%value(g$) /= 0 .or. ele%value(g_err$) /= 0)

  if (ele%key /= sad_mult$) ele%sub_key = ele%key  ! Save sbend/rbend input type.
  angle = ele%value(angle$) 

  ! Only one of b_field, g, or rho may be set.
  ! B_field may not be set for an rbend since, in this case, L is not computable (we don't know the ref energy).

  if (b_field_set .and. ele%key == rbend$) call parser_error &
          ("B_FIELD NOT SETTABLE FOR AN RBEND (USE AN SBEND INSTEAD): " // ele%name)

  if (b_field_set .and. g_set) call parser_error &
          ('BOTH G (OR G_ERR) AND B_FIELD (OR B_FIELD_ERR) SET FOR A BEND: ' // ele%name)

  if (b_field_set .and. ele%value(rho$) /= 0) call parser_error &
          ('BOTH RHO AND B_FIELD (OR B_FIELD_ERR) SET FOR A BEND: ' // ele%name)

  if (ele%value(g$) /= 0 .and. ele%value(rho$) /= 0) &
            call parser_error ('BOTH G AND RHO SPECIFIED FOR BEND: ' // ele%name)

  ! if rho is set then this gives g

  if (ele%value(l$) /= 0 .and. angle /= 0 .and. ele%value(g$) /= 0) &
                      call parser_error ('ANGLE, G, AND L ARE ALL SPECIFIED FOR BEND: ' // ele%name)
  if (ele%value(l$) /= 0 .and. angle /= 0 .and. ele%value(rho$) /= 0) &
                      call parser_error ('ANGLE, RHO, AND L ARE ALL SPECIFIED FOR BEND: ' // ele%name)

  if (ele%value(rho$) /= 0) ele%value(g$) = 1 / ele%value(rho$)

  ! If g and angle are set then this determines l

  if (ele%value(g$) /= 0 .and. angle /= 0) ele%value(l$) = angle / ele%value(g$)

  ! Convert an rbend to an sbend

  if (ele%key == rbend$) then
    ! Note: L may not be zero if g and angle have both been specified and are non-zero.
    if (ele%value(l$) == 0) then
      if (ele%value(l_chord$) == 0) then
        angle = 0
      elseif (angle /= 0) then
        ele%value(l$) = ele%value(l_chord$) * angle / (2 * sin(angle/2))
      elseif (ele%value(g$) /= 0) then
        angle = 2 * asin(ele%value(l_chord$) * ele%value(g$) / 2)
        ele%value(l$) = ele%value(l_chord$) * angle / (2 * sin(angle/2))
      else  ! g and angle are zero.
        ele%value(l$) = ele%value(l_chord$)
      endif
    endif
    
    ele%value(e1$) = ele%value(e1$) + angle / 2
    ele%value(e2$) = ele%value(e2$) + angle / 2
    ele%key = sbend$

  endif

  ! 

  if (ele%value(angle$) /= 0 .and. ele%value(l$) == 0) then
    call parser_error ('THE BENDING ANGLE IS NONZERO IN A ZERO LENGTH BEND! ' // ele%name)
  elseif (ele%value(angle$) /= 0) then
    ele%value(g$) = ele%value(angle$) / ele%value(l$) 
  endif

  ! If fintx or hgapx are real_garbage then they have not been set.
  ! If so, set their valuse to fint and hgap.

  if (ele%value(hgapx$) == real_garbage$) ele%value(hgapx$) = ele%value(hgap$)
  if (ele%value(fintx$) == real_garbage$) ele%value(fintx$) = ele%value(fint$)

! Accept Use of Voltage for lcavities and vary the mode frequencies.

case (lcavity$) 

  if (ele%value(voltage$) /= 0) then
    if (ele%value(gradient$) /= 0) call parser_error &
                ('BOTH VOLTAGE AND GRADIENT NON-ZERO FOR A LCAVITY:', ele%name)
    if (ele%value(l$) == 0) then
      ele%value(gradient$) = 0
    else
      ele%value(gradient$) = ele%value(voltage$) / ele%value(l$)
    endif
  endif

! 

case (rfcavity$) 

  if (ele%value(rf_frequency$) /= 0 .and. ele%value(harmon$) /= 0) call parser_error &
              ('BOTH RF_FREQUENCY AND HARMON SET FOR RFCAVITY: ' // ele%name, &
               'SETTING OF HARMON WILL BE IGNORED!', level = s_warn$)

! for a periodic_type wiggler n_pole is a dependent attribute

case (wiggler$, undulator$)
  if (ele%sub_key == periodic_type$) then

    if (ele%value(l_pole$) == 0 .and. ele%value(n_pole$) /= 0) then
      ele%value(l_pole$) = ele%value(l$) / ele%value(n_pole$) 
    endif

  endif

!

case (match$)
  ele%value(match_end_input$) = ele%value(match_end$)
  ele%value(match_end_orbit_input$) = ele%value(match_end_orbit$)

! check for inconsistancies

case (solenoid$)
  if (ele%field_master .and. (ele%value(ks$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM FOR ELEMENT: ' // ele%name, &
       'BOTH STRENGTH (KS, HKICK, ETC.) AND FIELD SET FOR A SOLENOID.')

case (sol_quad$)
  if (ele%field_master .and. (ele%value(ks$) /= 0 .or. &
                            ele%value(k1$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM FOR ELEMENT: ' // ele%name, &
       'BOTH STRENGTH (K1, HKICK, ETC.) AND FIELD SET FOR A SOL_QUAD.')

case (quadrupole$)
  if (ele%field_master .and. (ele%value(k1$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM FOR ELEMENT: ' // ele%name, &
       'BOTH STRENGTH (K1, HKICK, ETC.) AND FIELD SET FOR A QUAD.')

case (sextupole$)
  if (ele%field_master .and. (ele%value(k2$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM FOR ELEMENT: ' // ele%name, &
       'BOTH STRENGTH (K2, HKICK, ETC.) AND FIELD SET FOR A SEXTUPOLE.')

case (octupole$)
  if (ele%field_master .and. (ele%value(k3$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM FOR ELEMENT: ' // ele%name, &
       'BOTH STRENGTH (K3, HKICK, ETC.) AND FIELD SET FOR A OCTUPOLE.')

case (hkicker$, vkicker$)
  if (ele%field_master .and. (ele%value(kick$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM FOR ELEMENT: ' // ele%name, &
       'BOTH STRENGTH AND BL_KICK SET FOR A H/VKICKER.')

case (elseparator$)
  if (ele%field_master .and. kick_set) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM FOR ELEMENT: ' // ele%name, &
       'BOTH KICK (HKICK OR VKICK) AND E_FIELD OR VOLTAGE SET FOR A ELSEPARATOR.')

  if (ele%field_master) then
    if (ele%value(voltage$) /= 0 .and. ele%value(e_field$) /= 0) call parser_error &
              ('INDEPENDENT VARIABLE PROBLEM FOR ELSEPARATOR: ' // ele%name, &
               'BOTH VOLTAGE AND E_FIELD SET FOR THIS ELEMENT.')

    if (ele%value(voltage$) /= 0) then
      if (ele%value(gap$) == 0) then
        call parser_error ('FOR ELSEPARATOR: ' // ele%name, 'VOLTAGE IS SET BUT GAP IS NOT!')
      else
        ele%value(e_field$) = ele%value(voltage$) / ele%value(gap$)
      endif
    endif
  endif

case (e_gun$)
  if (ele%value(gradient$) == 0 .and. ele%value(l$) /= 0) ele%value(gradient$) = ele%value(voltage$) / ele%value(l$)

end select

end subroutine settable_dep_var_bookkeeping 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine form_digested_bmad_file_name (lat_file, digested_file, full_lat_file, use_line)
!
! Subroutine to form the standard name of the Bmad digested file. 
! The standard digested file name has the suffix added to the file name:
!     suffix = '.digested' + bmad_inc_version$ 
! Exception: If the use_line argument is present and not blank, the suffix will be:
!     suffix = '.' + use_line + '.digested' + bmad_inc_version$ 
!   
!
! Modules needed:
!   use bmad_parser_mod
!
! Input:
!   lat_file -- Character(*): Input lattice file name.
!   use_line -- Character(*), optional: Line used for lattice expansion. If not present
!                 or blank, the line used is the one that was specified in the lattice file.
!
! Output:
!   digested_file -- Character(200): Name of the digested file.
!   full_lat_file -- Character(200), optional: Input lattice file name with full directory.
!                       Can be used for error messages.
!-

subroutine form_digested_bmad_file_name (lat_file, digested_file, full_lat_file, use_line)

character(*) lat_file, digested_file
character(*), optional :: full_lat_file, use_line
character(200) full_name

integer ix

! Get the full_lat_file name

call fullfilename (lat_file, full_name)
inquire (file = full_name, name = full_name)  ! full input file_name
if (present (full_lat_file)) full_lat_file = full_name

! Construct the digested_file name

if (present(use_line)) then
  if (use_line /= '') then
    write (digested_file, '(4a, i0)') trim(full_name), '.', trim(use_line), '.digested', bmad_inc_version$ 
    return
  endif
endif

write (digested_file, '(2a, i0)') trim(full_name), '.digested', bmad_inc_version$ 

end subroutine form_digested_bmad_file_name

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_add_branch (fork_ele, lat, sequence, in_name, in_indexx, &
!                                                        seq_name, seq_indexx, no_end_marker, in_lat, plat, created_new_branch)
!
! Subroutine to do line expansion.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_add_branch (fork_ele, lat, sequence, in_name, in_indexx, &
                                    seq_name, seq_indexx, no_end_marker, in_lat, plat, created_new_branch, new_branch_name)

implicit none

type (lat_struct), target :: lat, in_lat
type (parser_lat_struct) plat
type (ele_struct) fork_ele
type (ele_struct), pointer :: target_ele
type (seq_struct), target :: sequence(:)
type (branch_struct), pointer :: branch

integer, allocatable :: seq_indexx(:), in_indexx(:)
integer i, j, nb, n_ele_use, n, ix, key

character(*), optional :: new_branch_name
character(*), allocatable ::  in_name(:), seq_name(:)

logical created_new_branch, no_end_marker

!

created_new_branch = .true.

if (fork_ele%value(new_branch$) == 0) then ! Branch back if
  do i = 0, ubound(lat%branch, 1) - 1
    branch => lat%branch(i)
    if (branch%name /= fork_ele%component_name) cycle
    fork_ele%value(ix_to_branch$) = i
    created_new_branch = .false.
    if (present(new_branch_name)) new_branch_name = ''
  enddo
endif

if (created_new_branch) then
  call parser_expand_line (lat, fork_ele%component_name, sequence, in_name, &
                                in_indexx, seq_name, seq_indexx, in_lat, n_ele_use, no_end_marker)

  nb = ubound(lat%branch, 1)
  fork_ele%value(ix_to_branch$) = nb
  branch => lat%branch(nb)

  branch%ix_from_branch     = fork_ele%ix_branch
  branch%ix_from_ele        = fork_ele%ix_ele
endif

if (present(new_branch_name)) new_branch_name = fork_ele%component_name
fork_ele%component_name = plat%ele(fork_ele%ixx)%ele_name  ! Substitute element name for line name.

end subroutine parser_add_branch

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_identify_fork_to_element (lat)
!
! Routine to identify the elements the forks in a lattice are branching to.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_identify_fork_to_element (lat)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: target_ele
type (ele_struct), pointer :: fork_ele
type (branch_struct), pointer :: branch

integer ib, ie, j

character(40) name

!

do ib = 0, ubound(lat%branch, 1)
  do ie = 1, lat%branch(ib)%n_ele_max

    fork_ele => lat%branch(ib)%ele(ie)
    if (fork_ele%key /= fork$ .and. fork_ele%key /= photon_fork$) cycle

    branch => lat%branch(nint(fork_ele%value(ix_to_branch$)))
    nullify(target_ele)
    name = fork_ele%component_name

    if (name == '') then
      if (nint(fork_ele%value(direction$)) == 1) then
        target_ele => branch%ele(0)
      else
        target_ele => branch%ele(branch%n_ele_track)
      endif
    else
      do j = 0, branch%n_ele_max
        if (branch%ele(j)%name /= name) cycle
        if (associated(target_ele)) then
          call parser_error('DUPLICATE TO_ELEMENT: ' // name, 'FOR FORK ELEMENT: ' // fork_ele%name)
          exit
        endif
        target_ele => branch%ele(j)
      enddo
    endif

    if (.not. associated(target_ele)) then
      call parser_error('TO_ELEMENT NOT FOUND: ' // name, 'FOR FORK ELEMENT: ' // fork_ele%name)
      return
    endif

    fork_ele%value(ix_to_element$) = target_ele%ix_ele

    select case (target_ele%key)
    case (marker$, fork$, photon_fork$, fiducial$, beginning_ele$)
    case default
      call parser_error('TO_ELEMENT: ' // name, 'FOR FORK ELEMENT: ' // fork_ele%name, &
                        'IS NOT A ZERO-LENGTH MARKER-LIKE ELEMENT')
    end select

  enddo
enddo

end subroutine parser_identify_fork_to_element

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_expand_line (lat, use_name, sequence, in_name, &
!                in_indexx, seq_name, seq_indexx, in_lat, n_ele_use, no_end_marker, expanded_line)
!
! Subroutine to do line expansion.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
! Input:
!   use_name      -- character(*): Root line to expand.
!   sequence(:)   -- seq_struct: Array of sequencies.
!   in_name(:)    -- character(*): Array of element names.
!   in_indexx(:)  -- integer: Index array of for the element names.
!   seq_name(:)   -- character(*): Array of sequence names.
!   seq_indexx(:) -- integer: Index array for the sequence names.
!   in_lat        -- lat_struct: Lattice with array of elements defined in the lattice file.
!   no_end_marker -- logical: Put a marker named "end" at the end of the branch?
!
! Output:
!   lat           -- lat_struct: Lattice with new line. Except if expanded_line is present.
!   n_ele_use     -- integer: Number of elements in the finished line.
!   expanded_line(:) -- character(*), allocatable, optional: If present, lat argument will be
!                         ignored and expanded_line will have the expanded line.
!                         This arg is used for girder lords.
!-

subroutine parser_expand_line (lat, use_name, sequence, in_name, &
               in_indexx, seq_name, seq_indexx, in_lat, n_ele_use, no_end_marker, expanded_line)

implicit none

type (lat_struct), target :: lat, in_lat
type (ele_struct), pointer :: ele_line(:), ele, ele2
type (seq_struct), target :: sequence(:)
type (seq_ele_struct), pointer :: s_ele, this_seq_ele
type (seq_stack_struct) stack(40)
type (seq_struct), pointer :: seq, seq2
type (used_seq_struct), allocatable ::  used2(:)
type (seq_ele_struct), target :: dummy_seq_ele
type (branch_struct), pointer :: branch
type (used_seq_struct), allocatable ::  used_line(:)

integer, allocatable :: seq_indexx(:), in_indexx(:)
integer iseq_tot, i_lev, i_use, n_ele_use, n_max
integer i, j, k, n, ix, ix_multipass, ix_branch, flip

character(*), allocatable ::  in_name(:), seq_name(:)
character(*), allocatable, optional :: expanded_line(:)
character(*) use_name
character(40) name

logical no_end_marker

! find line corresponding to the "use" statement.

iseq_tot = size(seq_indexx)
n_max = in_lat%n_ele_max
allocate (used_line(n_max))

call find_indexx (use_name, seq_name, seq_indexx, iseq_tot, i_use)
if (i_use == 0) then
  call parser_error ('CANNOT FIND DEFINITION OF LINE IN "USE" STATEMENT: ' // use_name, stop_here = .true.)
  return
endif

if (sequence(i_use)%type /= line$) then
  call parser_error ('NAME IN "USE" STATEMENT IS NOT A LINE!', stop_here = .true.)
  return
endif

! Now to expand the lines and lists to find the elements to use.
! First go through the lines and lists and index everything.

do k = 1, iseq_tot
  do i = 1, size(sequence(k)%ele(:))

    s_ele => sequence(k)%ele(i)
    name = s_ele%name

    if (s_ele%ix_arg > 0) then   ! dummy arg
      s_ele%type = element$
      cycle
    endif

    ! Remember: sequence names also appear in the element list so search the sequence list first.

    call find_indexx (name, seq_name, seq_indexx, iseq_tot, j)
    if (j == 0) then  ! if not an sequence, it must be an element
      call find_indexx (name, in_name, 0, in_indexx, n_max, j)
      if (j < 0) then  ! if not an element, I don't know what it is
        s_ele%ix_ele = -1       ! Struggle on for now...
        s_ele%type = element$
      else
        s_ele%ix_ele = j
        s_ele%type = element$
      endif
    else
      s_ele%ix_ele = j
      s_ele%type = sequence(j)%type
      if (s_ele%type == list$ .and. s_ele%ele_order_reflect) call parser_error ( &
                          'A REFLECTION WITH A LIST IS NOT ALLOWED IN: '  &
                          // sequence(k)%name, 'FOR LIST: ' // s_ele%name, &
                          seq = sequence(k))
      if (sequence(k)%type == list$) &
                call parser_error ('A REPLACEMENT LIST: ' // sequence(k)%name, &
                'HAS A NON-ELEMENT MEMBER: ' // s_ele%name)
    endif

  enddo
enddo

! to expand the "used" line we use a stack for nested sublines.
! used_line(:) is the expanded array of elements in the lat.
! init stack

i_lev = 1                          ! level on the stack
seq => sequence(i_use)

stack(1)%ix_seq    = i_use           ! which sequence to use for the lat
stack(1)%ix_ele    =  1              ! we start at the beginning
stack(1)%ele_order_direction = +1              ! element order is forward
stack(1)%orientation_direction = +1            ! and propagate forward through elements
stack(1)%rep_count = seq%ele(1)%rep_count
stack(1)%multipass = seq%multipass
stack(1)%tag = ''

n_ele_use = 0
         
sequence(:)%ix = 1  ! Init. Used for replacement list index

! Note: if present(expanded_line) => expansion is for getting a girder slave list.

if (stack(1)%multipass) then
  if (present(expanded_line)) then
    ix_multipass = 1
  else
    call parser_error ('"USE"D LINE FOR LATTICE EXPANSION IS MARKED MULTIPASS!')
    if (global_com%exit_on_error) call err_exit
  endif
endif

!-------------------------------------------------------------------------
! Expand "used" line...

line_expansion: do

  ! If rep_count is zero then we are finished with the current element.

  if (stack(i_lev)%rep_count == 0) then      ! goto next element in the sequence
    ! goto the next element by changing %ix_ele index by +/- 1 
    stack(i_lev)%ix_ele = stack(i_lev)%ix_ele + stack(i_lev)%ele_order_direction 
    ix = stack(i_lev)%ix_ele

    ! Check if off the end of the current line...
    if (ix > 0 .and. ix <= size(seq%ele)) then  ! Nope. Still have more element to process
      stack(i_lev)%rep_count = seq%ele(ix)%rep_count  ! set the rep_count for the next ele.
      if (stack(i_lev)%rep_count == 0) cycle          ! For "0*sub_line" construct.

    ! If we have got to the end of the current line then pop the stack back to
    ! the next lower level.
    else  
      i_lev = i_lev - 1
      if (i_lev == 0) exit line_expansion    ! level 0 -> we are done.
      seq => sequence(stack(i_lev)%ix_seq)
      cycle
    endif

  endif

  stack(i_lev)%rep_count = stack(i_lev)%rep_count - 1

  ! if s_ele is a dummy arg then get corresponding actual arg.

  s_ele => seq%ele(stack(i_lev)%ix_ele)  ! next element, line, or list

  ix = s_ele%ix_arg
  if (ix /= 0) then  ! it is a dummy argument.
    name = seq%corresponding_actual_arg(ix)
    s_ele => dummy_seq_ele
    s_ele%name = name
    s_ele%ele_orientation = seq%ele(stack(i_lev)%ix_ele)%ele_orientation
    call find_indexx (name, in_name, 0, in_indexx, n_max, j)
    if (j < 0) then  ! if not an element it must be a sequence
      call find_indexx (name, seq_name, seq_indexx, iseq_tot, j)
      if (j == 0) then  ! if not a sequence then I don't know what it is
        call parser_error ('CANNOT FIND DEFINITION FOR: ' // name, &
                          'IN LINE: ' // seq%name, seq = seq)
        if (global_com%exit_on_error) call err_exit
        return
      endif
      s_ele%ix_ele = j
      s_ele%type = sequence(j)%type
    else
      s_ele%ix_ele = j 
      s_ele%type = element$
    endif
    
  endif

  ! if an element

  select case (s_ele%type)

  case (element$, list$) 

    if (s_ele%type == list$) then
      seq2 => sequence(s_ele%ix_ele)
      j = seq2%ix
      this_seq_ele => seq2%ele(j)
      seq2%ix = seq2%ix + 1
      if (seq2%ix > size(seq2%ele(:))) seq2%ix = 1
    else
      if (s_ele%tag /= '') then
        call parser_error ('ELEMENTS IN A LINE OR LIST ARE NOT ALLOWED TO HAVE A TAG.', &
                      'FOUND ILLEGAL TAG FOR ELEMENT: ' // s_ele%name, &
                      'IN THE LINE/LIST: ' // seq%name, seq = seq)
      endif
      this_seq_ele => s_ele
    endif

    if (this_seq_ele%ix_ele < 1) call parser_error('NOT A DEFINED ELEMENT: ' // &
                          s_ele%name, 'IN THE LINE/LIST: ' // seq%name, seq = seq)

    if (n_ele_use+2 > size(used_line)) then
      n = 1.5*n_ele_use
      ix = size(used_line) 
      if (ix < n) then
        call move_alloc (used_line, used2)
        allocate (used_line(1:n))
        used_line(1:ix) = used2(1:ix)
        deallocate (used2)
      endif
    endif

    n_ele_use = n_ele_use + 1
    used_line(n_ele_use)%ix_ele_in_in_lat = this_seq_ele%ix_ele

    used_line(n_ele_use)%name = this_seq_ele%name
    used_line(n_ele_use)%orientation = stack(i_lev)%orientation_direction * this_seq_ele%ele_orientation

    if (stack(i_lev)%tag /= '' .and. s_ele%tag /= '') then
      used_line(n_ele_use)%tag =  trim(stack(i_lev)%tag) // '.' // s_ele%tag
    elseif (s_ele%tag /= '') then
      used_line(n_ele_use)%tag = s_ele%tag
    else
      used_line(n_ele_use)%tag =  stack(i_lev)%tag
    endif

    if (stack(i_lev)%multipass) then
      ix_multipass = ix_multipass + 1
      used_line(n_ele_use)%ix_multi = ix_multipass + 1000000 * stack(i_lev)%ix_seq
    else
      used_line(n_ele_use)%ix_multi = 0
    endif

  ! if a line:
  !     a) move pointer on current level past line element
  !     b) go to the next higher level
  !     c) initialize pointers for the higher level to use the line

  case (line$, replacement_line$)
    i_lev = i_lev + 1
    if (i_lev > size(stack)) then
      call parser_error ('NESTED LINES EXCEED STACK DEPTH! SUSPECT INFINITE LOOP!')
      if (global_com%exit_on_error) call err_exit
    endif
    if (s_ele%type == replacement_line$) then
      seq2 => sequence(s_ele%ix_ele)
      if (size(seq2%dummy_arg) /= size(s_ele%actual_arg)) then
        call parser_error ('WRONG NUMBER OF ARGUMENTS FORREPLACEMENT LINE: ' // &
            s_ele%name, 'WHEN USED IN LINE: ' // seq%name, seq = seq)
        return
      endif
      arg_loop: do i = 1, size(seq2%dummy_arg)
        seq2%corresponding_actual_arg(i) = s_ele%actual_arg(i)
        if (associated(seq%dummy_arg)) then
          do j = 1, size(seq%dummy_arg)
            if (seq2%corresponding_actual_arg(i) == seq%dummy_arg(j)) then
              seq2%corresponding_actual_arg(i) = seq%corresponding_actual_arg(j)
              cycle arg_loop
            endif
          enddo
        endif
        name = seq2%corresponding_actual_arg(i)
      enddo arg_loop
    endif

    seq => sequence(s_ele%ix_ele)
    stack(i_lev)%ix_seq = s_ele%ix_ele
    stack(i_lev)%multipass = (stack(i_lev-1)%multipass .or. seq%multipass)

    if (stack(i_lev-1)%tag /= '' .and. s_ele%tag /= '') then
       stack(i_lev)%tag = trim(stack(i_lev-1)%tag) // '.' // s_ele%tag
    elseif (stack(i_lev-1)%tag /= '') then
       stack(i_lev)%tag = trim(stack(i_lev-1)%tag)
    else
       stack(i_lev)%tag = s_ele%tag
    endif

    stack(i_lev)%orientation_direction = stack(i_lev-1)%orientation_direction * s_ele%ele_orientation

    stack(i_lev)%ele_order_direction = stack(i_lev-1)%ele_order_direction
    if (s_ele%ele_order_reflect) stack(i_lev)%ele_order_direction = -stack(i_lev)%ele_order_direction

    if (stack(i_lev)%ele_order_direction == 1) then
      ix = 1
    else
      ix = size(seq%ele(:))
    endif

    stack(i_lev)%ix_ele = ix
    stack(i_lev)%rep_count = seq%ele(ix)%rep_count

    if (seq%multipass) then
      ix_multipass = 1
    endif

  case default
    call parser_error ('INTERNAL SEQUENCE ERROR!')

  end select

enddo line_expansion

! Stop here if there has been an error

if (bp_com%error_flag) return

! expanded_line present?

if (present(expanded_line)) then
  if (allocated(expanded_line)) deallocate(expanded_line)
  allocate(expanded_line(n_ele_use))
  do i = 1, n_ele_use
    expanded_line(i) = used_line(i)%name
  enddo
  return
endif

! Transfer the ele information from the in_lat to lat and
! do the bookkeeping for settable dependent variables.

if (lat%n_ele_max < 1) then
  ix_branch = 0
else
  ix_branch = ubound(lat%branch, 1) + 1
endif

if (ix_branch == 0) then  ! Main branch
  call allocate_lat_ele_array(lat, n_ele_use+1)
  lat%n_ele_track = n_ele_use
  lat%n_ele_max   = n_ele_use
  ele_line => lat%ele
else                    ! branch line
  call allocate_branch_array (lat, ix_branch)
  call allocate_lat_ele_array(lat, n_ele_use+1, ix_branch)
  ele_line => lat%branch(ix_branch)%ele
endif

branch => lat%branch(ix_branch)

do i = 1, n_ele_use
  ele_line(i) = in_lat%ele(used_line(i)%ix_ele_in_in_lat) 
  ele_line(i)%name        = used_line(i)%name
  ele_line(i)%iyy         = used_line(i)%ix_multi
  ele_line(i)%orientation = used_line(i)%orientation
  if (used_line(i)%tag /= '') ele_line(i)%name = trim(used_line(i)%tag) // '.' // ele_line(i)%name
  call settable_dep_var_bookkeeping (ele_line(i))
  if (ele_line(i)%lord_status == super_lord$) then
    call parser_error ('SUPERPOSITION ELEMENTS CANNOT BE USED IN A LINE: ' // ele_line(i)%name)
  endif
enddo

ele_line(0)%ix_branch = ix_branch
ele_line(0)%orientation = ele_line(1)%orientation

deallocate(used_line)

! Add End marker and make sure it's orientation is consistant

if (.not. no_end_marker) then
  n_ele_use = n_ele_use + 1
  ele => ele_line(n_ele_use)
  ele%name = 'END'
  ele%key = marker$
  call set_ele_defaults (ele)
  flip = 1
  do j = n_ele_use-1, 0, -1
    ele2 => ele_line(j)
    if (ele2%key == patch$ .or. ele2%key == floor_shift$) then
      if (patch_flips_propagation_direction (ele2%value(x_pitch$), ele2%value(y_pitch$))) flip = -flip
      cycle
    endif
    exit
  enddo
  ele%orientation = ele2%orientation * flip
endif

! Branch info

branch%n_ele_track = n_ele_use
branch%n_ele_max   = n_ele_use
branch%ix_branch   = ix_branch
branch%name        = use_name

end subroutine parser_expand_line

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine bp_set_ran_status
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine bp_set_ran_status

if (bp_com%extra%deterministic == 0) then
  bp_com%extra%ran_function_was_called = .true.
elseif (bp_com%extra%deterministic == 1) then
  bp_com%extra%deterministic_ran_function_was_called = .true.
endif

end subroutine bp_set_ran_status
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_debug_print_info (debug_line)
!
! Subroutine to remove all null_ele elements.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_debug_print_info (lat, debug_line)

type (lat_struct) lat
character(*) debug_line
integer i, ix

!

call str_upcase (debug_line, debug_line)

if (index(debug_line, 'VAR') /= 0) then
  print *
  print *, '----------------------------------------'
  print *, 'Number of Defined Variables:', bp_com%ivar_tot - bp_com%ivar_init
  do i = bp_com%ivar_init+1, bp_com%ivar_tot
    print '(i6, 2x, a, es18.10)', i, bp_com%var(i)%name, bp_com%var(i)%value
  enddo
endif

if (index(debug_line, 'SLAVE') /= 0) then
  print *
  print *, '----------------------------------------'
  print *, 'Number of Elements in Tracking Lattice:', lat%n_ele_track
  do i = 1, lat%n_ele_track
    print *, '-------------'
    print *, 'Ele #', i
    call type_ele (lat%ele(i), .false., 0, .false., 0, .true., .true., .false., .true., .true.)
  enddo
endif

if (index(debug_line, 'LORD') /= 0) then
  print *
  print *, '----------------------------------------'
  print *, 'LORD elements: ', lat%n_ele_max - lat%n_ele_track
  do i = lat%n_ele_track+1, lat%n_ele_max
    print *, '-------------'
    print *, 'Ele #', i
    call type_ele (lat%ele(i), .false., 0, .false., 0, .true., .true., .false., .true., .true.)
  enddo
endif

if (index(debug_line, 'LATTICE') /= 0) then  
  print *
  print *, '----------------------------------------'
  print *, 'Lattice Used: ', lat%use_name
  print *, 'Number of lattice elements:', lat%n_ele_track
  print *, 'List:                                 Key                 Length         S'
  do i = 1, lat%n_ele_track
    print '(i4, 2a, 3x, a, 2f10.2)', i, ') ', lat%ele(i)%name(1:30),  &
      key_name(lat%ele(i)%key), lat%ele(i)%value(l$), lat%ele(i)%s
  enddo
  print *, '---- Lord Elements ----'
  do i = lat%n_ele_track+1, lat%n_ele_max
    print '(2x, i4, 2a, 3x, a, 2f10.2)', i, ') ', lat%ele(i)%name(1:30),  &
           key_name(lat%ele(i)%key), lat%ele(i)%value(l$), lat%ele(i)%s
  enddo
endif

ix = index(debug_line, 'ELE')
if (ix /= 0) then
  print *
  print *, '----------------------------------------'
  call string_trim (debug_line(ix+3:), debug_line, ix)
  do
    if (ix == 0) exit
    read (debug_line, *) i
    print *
    print *, '----------------------------------------'
    print *, 'Element #', i
    call type_ele (lat%ele(i), .false., 0, .true., 0, .true., .true., .true., .true., .true.)
    call string_trim (debug_line(ix+1:), debug_line, ix)
  enddo
endif

if (index(debug_line, 'BEAM_START') /= 0) then
  print *
  print *, '----------------------------------------'
  print *, 'beam_start:'
  print '(3x, 6es13.4)', lat%beam_start%vec      
endif

end subroutine parser_debug_print_info

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parse_wake_lr_spline (lr_pa, ele, lat, delim, delim_found, err_flag)
!
! Subroutine to parse a "lr_spline = {}" construct
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is private to bmad_parser_mod.
!-

subroutine parse_wake_lr_spline (lr_pa, ele, lat, delim, delim_found, err_flag)

implicit none

type (wake_lr_spline_struct) lr_pa
type (ele_struct), target :: ele
type (lat_struct), target :: lat

integer ix_word

character(40) attrib_name
character(1) delim

logical err_flag, delim_found

!

err_flag = .true.

do

  ! Read attriubute
  call get_next_word (attrib_name, ix_word, '{}=,()', delim, delim_found, call_check = .true.)
  if (.not. expect_this ('=', .true., .false., 'IN WAKE_LR_SPLINE DEFINITION', ele, delim, delim_found)) return

  select case (attrib_name)

  case ('T_MAX')
    call evaluate_value (ele%name, lr_pa%t_max, lat, delim, delim_found, err_flag, ',}')

  case default
    if (attrib_name == '') then
      call parser_error ('MANGLED WAKE_LR_SPLINE DEFINITION FOR ELEMENT: ' // ele%name)
    else
      call parser_error ('UNKNOWN WAKE_LR_SPLINE COMPONENT: ' // attrib_name, 'FOR ELEMENT: ' // ele%name)
    endif
    return

  end select

  ! Possible "}" is end of mode
  if (delim == '}') exit

enddo

if (.not. expect_one_of (', ', .false., ele, delim, delim_found)) return
err_flag = .false.

end subroutine parse_wake_lr_spline

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parse_cartesian_map (ct_map, ele, lat, delim, delim_found, err_flag)
!
! Subroutine to parse a "cartesian_map = {}" construct
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is private to bmad_parser_mod.
! This must read in:
! {type = ,
!    dr = , 
!    r0 = , 
!    pt(i,j,k) = ( (ex_re, ex_im), .... (bz_re, bz_im) ) 
!    .
!    .
!    . ) },
!-

subroutine parse_cartesian_map (ct_map, ele, lat, delim, delim_found, err_flag)

implicit none

type (cartesian_map_struct) :: ct_map
type (ele_struct), target :: ele
type (ele_struct), pointer :: match_ele
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (cartesian_map_term1_struct), allocatable :: term(:)
type (cartesian_map_term1_struct), pointer :: tm

real(rp) kx, ky, kz, tol

complex(rp), pointer :: c_ptr(:)

integer n, ix_word, i_term, ib, ie, im, ix, family

character(80) err_str
character(40) word, word2, name, attrib_name
character(1) delim, delim2

logical err_flag, delim_found

!

name = 'xxx'
err_flag = .true.

if (ele%key == wiggler$) then
  ct_map%master_parameter = polarity$
  ele%sub_key = map_type$
endif

!

do

  ! Read attriubute
  call get_next_word (attrib_name, ix_word, '{}=,()', delim, delim_found, call_check = .true.)
  if (.not. expect_this ('=', .true., .false., 'IN CARTESIAN_MAP DEFINITION', ele, delim, delim_found)) return

  select case (attrib_name)

  case ('FIELD_SCALE')
    call evaluate_value (ele%name, ct_map%field_scale, lat, delim, delim_found, err_flag, ',}')

  case ('R0')
    if (.not. equal_sign_here(ele, delim)) return
    if (.not. parse_real_list (lat, trim(ele%name) // ' GRID_FIELD', ct_map%r0, .true.)) return
    if (.not. expect_one_of (',}', .false., ele, delim, delim_found)) return


  case ('ELE_ANCHOR_PT', 'FIELD_TYPE', 'MASTER_PARAMETER')
    ! Expect "<component> = "
    if (delim /= '=') then
      call parser_error ('NO "=" SIGN FOUND AFTER ' // attrib_name,  &
                         'IN ELEMENT: ' // ele%name)
      return
    endif
    call get_next_word (word2, ix_word, ',}', delim, delim_found)

    !

    select case (attrib_name)
 
    case ('MASTER_PARAMETER')
      if (word2 == 'NONE') then
        ix = 0
      else
        ix = attribute_index(ele, word2)
        if (ix < 1 .or. ix > num_ele_attrib$) then
          call parser_error ('BAD NAME FOR "MASTER_PARAMETER = <NAME>" CONSTRUCT', &
                               'FOUND IN ELEMENT: ' // ele%name)
          return
        endif
      endif
      ct_map%master_parameter = ix

    case ('ELE_ANCHOR_PT')
      call match_word(word2, anchor_pt_name(1:), ct_map%ele_anchor_pt, can_abbreviate = .false., matched_name = name)
  
    case ('FIELD_TYPE')
      call match_word(word2, em_field_type_name(1:), ct_map%field_type, can_abbreviate = .false., matched_name = name)
  
    end select

    !

    if (name == '') then
      call parser_error ('UNKNKOWN ' // trim(attrib_name) // ' VALUE:' // word2, &
                         'IN ELEMENT: ' // ele%name)
      return        
    endif      
      
  case ('TERM')

    if (.not. allocated(ct_map%ptr%term)) then
      allocate (ct_map%ptr%term(1))
      tm => ct_map%ptr%term(1)
      ! Set %file to be the last called file:<line_number>. 
      ct_map%ptr%file = bp_com%line1_file_name
    else
      n = size(ct_map%ptr%term) + 1
      call move_alloc(ct_map%ptr%term, term)
      allocate (ct_map%ptr%term(n))
      ct_map%ptr%term(1:n-1) = term
      deallocate (term)
      tm => ct_map%ptr%term(n)
    endif

    err_str = trim(ele%name) // ' CARTESIAN_MAP TERM'

    if (.not. expect_this ('{', .false., .false., 'AFTER "TERM =" IN CARTESIAN_MAP DEFINITION', ele, delim, delim_found)) return
    call evaluate_value (err_str, tm%coef, lat, delim, delim_found, err_flag, ',');  if (err_flag) return
    call evaluate_value (err_str, tm%kx, lat, delim, delim_found, err_flag, ',');    if (err_flag) return
    call evaluate_value (err_str, tm%ky, lat, delim, delim_found, err_flag, ',');    if (err_flag) return
    call evaluate_value (err_str, tm%kz, lat, delim, delim_found, err_flag, ',');    if (err_flag) return
    call evaluate_value (err_str, tm%x0, lat, delim, delim_found, err_flag, ',');    if (err_flag) return
    call evaluate_value (err_str, tm%y0, lat, delim, delim_found, err_flag, ',');    if (err_flag) return
    call evaluate_value (err_str, tm%phi_z, lat, delim, delim_found, err_flag, ','); if (err_flag) return
    call get_switch ('FAMILY', ['X ', 'Y ', 'QU', 'SQ'], family, err_flag, ele, delim, delim_found); if (err_flag) return
    if (.not. expect_this ('}', .true., .false., 'AFTER "FAMILY" SWITCH', ele, delim, delim_found)) return
    if (.not. expect_one_of(',}', .false., ele, delim, delim_found)) return

    kx = tm%kx
    ky = tm%ky
    kz = tm%kz
    tol = 1d-5 * (kx**2 + ky**2 + kz**2)

    if (abs(ky**2 - kx**2 - kz**2) < tol) then
      select case (family)
      case (x_family$);   tm%type = hyper_y_family_x$
      case (y_family$);   tm%type = hyper_y_family_y$
      case (qu_family$);  tm%type = hyper_y_family_qu$
      case (sq_family$);  tm%type = hyper_y_family_sq$
      end select

    elseif (abs(ky**2 + kx**2 - kz**2) < tol) then
      select case (family)
      case (x_family$);   tm%type = hyper_xy_family_x$
      case (y_family$);   tm%type = hyper_xy_family_y$
      case (qu_family$);  tm%type = hyper_xy_family_qu$
      case (sq_family$);  tm%type = hyper_xy_family_sq$
      end select

    elseif (abs(ky**2 - kx**2 + kz**2) < tol) then
      select case (family)
      case (x_family$);   tm%type = hyper_x_family_x$
      case (y_family$);   tm%type = hyper_x_family_y$
      case (qu_family$);  tm%type = hyper_x_family_qu$
      case (sq_family$);  tm%type = hyper_x_family_sq$
      end select

    else
      call parser_error ('CARTESIAN_MAP TERM DOES NOT HAVE CONSISTANT Kx, Ky, and Kz', &
                    'FOR ELEMENT: ' // ele%name)
      err_flag = .true.
      return
    endif

  case default
    if (attrib_name == '') then
      call parser_error ('MANGLED CARTESIAN_MAP DEFINITION FOR ELEMENT: ' // ele%name)
    else
      call parser_error ('UNKNOWN CARTESIAN_MAP COMPONENT: ' // attrib_name, 'FOR ELEMENT: ' // ele%name)
    endif
    return

  end select

  ! Possible "}" is end of mode
  if (delim == '}') exit

enddo

!

if (.not. expect_one_of (', ', .false., ele, delim, delim_found)) return
err_flag = .false.

! Check if data has already been read in for another element.
! If so, save space by pointing to the data.

call find_matching_fieldmap(ct_map%ptr%file, ele, cartesian_map$, match_ele, im)
if (im > 0) then
  deallocate(ct_map%ptr)
  ct_map%ptr => match_ele%cartesian_map(im)%ptr
  ct_map%ptr%n_link = ct_map%ptr%n_link + 1        
endif

end subroutine parse_cartesian_map

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! parse_cylindrical_map (cl_map, ele, lat, delim, delim_found, err_flag)
!
! Subroutine to parse a "cylindrical_map = {}" construct
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is private to bmad_parser_mod.
! This must read in:
! {type = ,
!    dr = , 
!    r0 = , 
!    pt(i,j,k) = ( (ex_re, ex_im), .... (bz_re, bz_im) ) 
!    .
!    .
!    . ) },
!-

subroutine parse_cylindrical_map (cl_map, ele, lat, delim, delim_found, err_flag)

implicit none

type (cylindrical_map_struct), pointer :: cl_map
type (ele_struct), target :: ele
type (ele_struct), pointer :: match_ele
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch

real(rp), allocatable :: array(:)

complex(rp), pointer :: c_ptr(:)

integer ix_word, i_term, ib, ie, im, ix

character(1) delim, delim2
character(40) word, word2, name, attrib_name

logical err_flag, delim_found, file_name_set

! Init

err_flag = .true.
allocate (array(1024))
file_name_set = .false.

!

do

  ! Read attriubute
  call get_next_word (attrib_name, ix_word, '{}=,()', delim, delim_found, call_check = .true.)

  select case (attrib_name)

  case ('DPHI0_REF')
    call parser_error ('THE ATTRIBUTE NAME "DPHI0_REF" HAS BEEN CHANGED TO "PHI0_FIELDMAP"', &
                       'PLEASE MAKE THE CHANGE IN THE LATTICE FILE.')

  case ('PHI0_FIELDMAP')
    call evaluate_value (ele%name, cl_map%phi0_fieldmap, lat, delim, delim_found, err_flag, ',}')

  case ('THETA0_AZIMUTH')
    call evaluate_value (ele%name, cl_map%theta0_azimuth, lat, delim, delim_found, err_flag, ',}')

  case ('FIELD_SCALE')
    call evaluate_value (ele%name, cl_map%field_scale, lat, delim, delim_found, err_flag, ',}')

  case ('DZ')            
    call evaluate_value (trim(ele%name), cl_map%dz, lat, delim, delim_found, err_flag, ',}')

  case ('R0')
    if (.not. equal_sign_here(ele, delim)) return
    if (.not. parse_real_list (lat, trim(ele%name) // ' GRID_FIELD', cl_map%r0, .true.)) return
    if (.not. expect_one_of (',}', .false., ele, delim, delim_found)) return

  case ('M')
   call parser_get_integer (cl_map%m, word, ix_word, delim, delim_found, err_flag, 'BAD CYLINDRICAL_MAP M CONSTRUCT', 'IN ELEMENT: ' // ele%name)

  case ('HARMONIC')
    call parser_get_integer (cl_map%harmonic, word, ix_word, delim, delim_found, err_flag, 'BAD CYLINDRICAL_MAP HARMONIC CONSTRUCT', 'IN ELEMENT: ' // ele%name)

  case ('MASTER_PARAMETER')
    call get_next_word (word, ix_word, ',}', delim, delim_found)
    if (word == 'NONE') then
      ix = 0
    else
      ix = attribute_index(ele, word)
      if (ix < 1 .or. ix > num_ele_attrib$) then
        call parser_error ('BAD NAME FOR "MASTER_PARAMETER = <NAME>" CONSTRUCT', &
                             'FOUND IN ELEMENT: ' // ele%name)
        return
      endif
    endif
    cl_map%master_parameter = ix


  case ('ELE_ANCHOR_PT')
    ! Expect "<component> = "
    if (delim /= '=') then
      call parser_error ('NO "=" SIGN FOUND AFTER ELE_ANCHOR_PT ' // attrib_name,  &
                         'IN ELEMENT: ' // ele%name)
      return
    endif
    call get_next_word (word2, ix_word, ',}', delim, delim_found)

    ! Evaluate string into integer.

    call match_word(word2, anchor_pt_name(1:), cl_map%ele_anchor_pt, can_abbreviate = .false., matched_name = name)
  
    if (name == '') then
      call parser_error ('UNKNKOWN ELE_ANCHOR_PT ' // trim(word) // ': ' // word2, &
                         'FOUND IN ELEMENT: ' // ele%name)
      return        
    endif      
      
  case ('E_COEF_RE', 'E_COEF_IM', 'B_COEF_RE', 'B_COEF_IM')

    if (.not. file_name_set) then
      ! Set %file to be the last called file:<line_number>. 
      cl_map%ptr%file = bp_com%line1_file_name
      file_name_set = .true.
    endif

    ! Expect "("
    call get_next_word (word, ix_word, ',({', delim, delim_found)
    if (word /= '' .or. delim /= '(') then
      call parser_error ('NO "(" FOUND AFTER "' // trim(attrib_name) // ' =" ', &
                           'IN ELEMENT: ' // ele%name)
      return
    endif

    ! Read list of values.
    call re_allocate(array, 1024, .false.)
    do i_term = 1, 100000
      call get_next_word (word, ix_word, '{},()', delim, delim_found)
      if ((delim /= ',' .and. delim /= ')') .or. .not. is_real(word)) then
        call parser_error ('ERROR PARSING CYLINDRICAL_MAP ARRAY: ' // word2, &
                             'IN ELEMENT: ' // ele%name)
        return
      endif
      if (i_term > size(array)) call re_allocate(array, 2*size(array))
      read (word, *) array(i_term)
      if (delim == ')') exit
    enddo

    if (allocated(cl_map%ptr%term)) then
      if (size(cl_map%ptr%term) /= i_term) then
        call parser_error ('ARRAY SIZE MISMATCH FOR: ' // word2, &
                           'IN CYLINDRICAL_MAP DEFINITION IN ELEMENT: ' // ele%name)
        return
      endif
    else
      allocate(cl_map%ptr%term(i_term))
    endif

    select case (attrib_name)
    case ('E_COEF_RE', 'E_COEF_IM'); c_ptr => cl_map%ptr%term%e_coef 
    case ('B_COEF_RE', 'B_COEF_IM'); c_ptr => cl_map%ptr%term%b_coef
    end select

    if (attrib_name(8:9) == 'RE') then
      if (any(real(c_ptr) /= 0)) then
        call parser_error ('DUPLICATE ARRAY FOR: ' // attrib_name, &
                           'IN CYLINDRICAL_MAP IN ELEMENT: ' // ele%name)
        return
      endif
      c_ptr = c_ptr + array(1:i_term)

    else
      if (any(aimag(c_ptr) /= 0)) then
        call parser_error ('DUPLICATE ARRAY FOR: ' // attrib_name, &
                           'IN CYLINDRICAL_MAP IN ELEMENT: ' // ele%name)
        return
      endif
      c_ptr = c_ptr + i_imaginary * array(1:i_term)
    endif

    ! Expect "," or "}"
    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
    if (word /= '' .or. (delim /= ',' .and. delim /= '}')) then
      call parser_error ('BAD ' // trim(attrib_name) // ' = (...) CONSTRUCT', &
                           'FOUND IN CYLINDRICAL_MAP DEFINITION IN ELEMENT: ' // ele%name)
      return
    endif

  case default
    if (attrib_name == '') then
      call parser_error ('MANGLED CYLINDRICAL_MAP DEFINITION FOR ELEMENT: ' // ele%name)
    else
      call parser_error ('UNKNOWN CYLINDRICAL_MAP COMPONENT: ' // attrib_name, &
                         'FOR ELEMENT: ' // ele%name)
    endif
    return
    
  end select

  ! Possible "}" is end of mode
  if (delim == '}') exit

enddo

! Get final separator after grid construct.
 
if (.not. expect_one_of (', ', .false., ele, delim, delim_found)) return

deallocate(array)
err_flag = .false.

! Check if data has already been read in for another element.
! If so, save space by pointing to the data.

call find_matching_fieldmap(cl_map%ptr%file, ele, cylindrical_map$, match_ele, im)
if (im > 0) then
  deallocate(cl_map%ptr)
  cl_map%ptr => match_ele%cylindrical_map(im)%ptr
  cl_map%ptr%n_link = cl_map%ptr%n_link + 1        
endif

end subroutine parse_cylindrical_map

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! parse_grid_field (grid, ele, lat, delim, delim_found, err_flag)
!
! Subroutine to parse a "grid_field = {}" construct
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is private to bmad_parser_mod.
! This must read in:
! {type = ,
!    dr = , 
!    r0 = , 
!    pt(i,j,k) = ( (ex_re, ex_im), .... (bz_re, bz_im) ) 
!    .
!    .
!    . ) },
!-

subroutine parse_grid_field(g_field, ele, lat, delim, delim_found, err_flag)

type grid_pt_struct
  integer :: ix(3) = [1, 1, 1]
  complex(rp) :: field(6) = 0
end type

type (grid_field_struct), pointer :: g_field
type (ele_struct) :: ele
type (ele_struct), pointer :: match_ele
type (lat_struct),  target :: lat
type (branch_struct), pointer :: branch
type (grid_pt_struct), allocatable :: array(:), array2(:)

character(1) delim, delim2
character(40) :: word, word2, name

integer ix_word, ix_word2, ix
integer pt_counter, n, i, ib, ie, im, ix0, ix1, iy0, iy1, iz0, iz1
integer grid_dim,  num_dr, num_r0, ios

logical delim_found, delim_found2, err_flag, err_flag2

! Init. Last thing read in was initial "{"

allocate(array(1024))
pt_counter = 0
err_flag = .true.

do    

  ! Read attriubute
  call get_next_word (word, ix_word, '{}=,()', delim, delim_found, call_check = .true.)

  select case (word)

  case ('PHI0_FIELDMAP')
    call evaluate_value (ele%name, g_field%phi0_fieldmap, lat, delim, delim_found, err_flag, ',}')

  case ('FIELD_SCALE')
    call evaluate_value (ele%name, g_field%field_scale, lat, delim, delim_found, err_flag, ',}')

  case ('HARMONIC')
    call parser_get_integer (g_field%harmonic, word, ix_word, delim, delim_found, err_flag, 'BAD GRID_FIELD HARMONIC CONSTRUCT', 'IN ELEMENT: ' // ele%name)

  case ('MASTER_PARAMETER')
    call get_next_word (word, ix_word, ',}', delim, delim_found)
    if (word == 'NONE') then
      ix = 0
    else
      ix = attribute_index(ele, word)
      if (ix < 1 .or. ix > num_ele_attrib$) then
        call parser_error ('BAD NAME FOR GRID_FIELD MASTER_PARAMETER', 'FOUND IN ELEMENT: ' // ele%name)
        return
      endif
    endif
    g_field%master_parameter = ix

  case ('FIELD_TYPE', 'ELE_ANCHOR_PT', 'GEOMETRY')
    if (.not. equal_sign_here(ele, delim)) return

    call get_next_word (word2, ix_word, ',}', delim, delim_found)
    ! Check to see if this is a valid type by checking against grid_field_geometry_name(:)

    if (word == 'FIELD_TYPE') then
      call match_word(word2, em_field_type_name, g_field%field_type, can_abbreviate = .false., matched_name = name)
    elseif (word == 'GEOMETRY') then
      call match_word(word2, grid_field_geometry_name(1:), g_field%geometry, can_abbreviate = .false., matched_name = name)
    else
      call match_word(word2, anchor_pt_name(1:), g_field%ele_anchor_pt, can_abbreviate = .false., matched_name = name)
    endif
  
    if (name == '') then
      call parser_error ('UNKNKOWN GRID_FIELD ' // trim(word) // ': ' // word2, &
                         'FOUND IN GRID_FIELD DEFINITION FOR ELEMENT: ' // ele%name)
      return        
    endif      
      
  case ('CURVED_COORDS', 'CURVED_REF_FRAME')  ! 'curved_coords' is old style.
    if (.not. equal_sign_here(ele, delim)) return
    call get_next_word (word2, ix_word, ':,=()', delim, delim_found, .true.)
    g_field%curved_ref_frame = evaluate_logical (word2, ios)
    if (ios /= 0 .or. ix_word == 0) then
      call parser_error ('BAD GRID_FIELD CURVED_REF_FRAME SETTING ' // word2, 'FOR: ' // ele%name)
    endif

  case ('R0')
    if (.not. equal_sign_here(ele, delim)) return
    if (.not. parse_real_list (lat, trim(ele%name) // ' GRID_FIELD', g_field%r0, .false.)) return
    if (.not. expect_one_of (',}', .false., ele, delim, delim_found)) return

    case ('DR')
    if (.not. equal_sign_here(ele, delim)) return
    ! expect ( 1.) or (1. , 2.) or (1., 2., 3.)
    if (.not. parse_real_list (lat, trim(ele%name) // ' GRID', g_field%dr, .false.)) return
    call get_next_word (word, ix_word, ',}', delim, delim_found)     
    if (word /= '') then
      call parser_error ('BAD INPUT AFTER DR DEFINITION: ' // word , &
                                 'FOUND IN GRID_FIELD DEFINITION FOR ELEMENT: ' // ele%name)
      return
    end if

  case ('PT')

    if (pt_counter == 0) then
      ! Set %file to be the last called file:<line_number>.
      g_field%ptr%file = bp_com%line1_file_name
    endif

    ! Increment 
    pt_counter = pt_counter + 1
    ! Reallocate temporary structure if needed
    n = size(array)
    if (pt_counter > n) then
      call move_alloc(array, array2)
      allocate(array(2*n))
      array(1:n) = array2
      deallocate(array2)
    end if

    ! Get indices
    bp_com%parse_line = delim // bp_com%parse_line
    if (.not. parse_integer_list (trim(ele%name) // ' GRID_FIELD PT', lat, array(pt_counter)%ix, .false.)) return
      
    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
    call get_next_word (word2, ix_word2, '{}=,()', delim2, delim_found2)
    if ((word /= '') .or. (word2 /= '') .or. (delim /= '=') .or. (delim2 /= '(')) then
      call parser_error ('BAD GRID_FIELD PT CONSTRUCT, NO  = "(" ', &
                 'FOUND IN GRID_FIELD DEFINITION FOR ELEMENT: ' // ele%name)
      return
    end if
    ! Get as many field components as listed
    do i = 1, 6
      call parse_complex_component(array(pt_counter)%field(i), delim, err_flag2)
      if (err_flag2) return
      if (delim == ')') exit
      if (delim /= ',') then
        call parser_error ('BAD GRID_FIELD PT CONSTRUCT, NO "," BETWEEN FIELD COMPONENTS', &
            'FOUND IN GRID_FIELD DEFINITION FOR ELEMENT: ' // ele%name)
        return
      end if
    end do

    select case (i)
    case (3)
      if (g_field%field_type == mixed$) then
        call parser_error ('FIELD_GRID WITH FIELD_TYPE = MIXED IN ELEMENT: ' // ele%name, 'DOES NOT SPECIFY BOTH E AND B FIELDS.')
        return
      endif        
    case (6)
      if (g_field%field_type /= mixed$) then
        call parser_error ('FIELD_GRID WITH FIELD_TYPE = ELECTRIC OR MAGNETIC IN ELEMENT: ' // ele%name, 'CANNOT SPECIFY BOTH E AND B FIELDS.')
        return
      endif
    case default
      call parser_error ('FIELD_GRID IN ELEMENT: ' // ele%name, 'DOES NOT HAVE THE CORRECT NUMBER OF FIELD COMPOENTS IN PT = (...) CONSTRUCT')
      return
    end select

    ! Expect , or }
    if (.not. expect_one_of(',}', .false., ele, delim, delim_found)) return
  
  case default
    if (word == '') then
      call parser_error ('MANGLED GRID_FIELD DEFINITION FOR ELEMENT: ' // ele%name)
    else
      call parser_error ('UNKNOWN GRID_FIELD COMPONENT: ' // word, &
                         'FOR ELEMENT: ' // ele%name)
    endif
    return
    
  end select 

  if (delim == '}') exit   

enddo

! Get final separator after grid_field construct.
 
if (.not. expect_one_of (', ', .false., ele, delim, delim_found)) return

! Clear pts

if (allocated(g_field%ptr%pt)) deallocate(g_field%ptr%pt)

! Allocate grid_field for different dimensions

grid_dim = grid_field_dimension(g_field%geometry)

if (grid_dim < 1 .or. grid_dim > 3) then
  call parser_error ('BAD GRID_FIELD DIMENSION', &
             'FOUND IN GRID_FIELD DEFINITION FOR ELEMENT: ' // ele%name)
  return
endif

ix0 = minval(array(1:pt_counter)%ix(1))
ix1 = maxval(array(1:pt_counter)%ix(1))
iy0 = minval(array(1:pt_counter)%ix(2))
iy1 = maxval(array(1:pt_counter)%ix(2))
iz0 = minval(array(1:pt_counter)%ix(3))
iz1 = maxval(array(1:pt_counter)%ix(3))

allocate(g_field%ptr%pt(ix0:ix1, iy0:iy1, iz0:iz1))

n = (ix1+1-ix0) * (iy1+1-iy0) * (iz1+1-iz0)
if (n /= pt_counter) then
  call out_io (s_warn$, bp_com%parser_name, &
                 'Note: Number of grid_field points (\i0\) in the file not equal to grid_field array size (\i0\ x \i0\ x \i0\).', &
                 'for element: ' // ele%name, &
                 i_array = [pt_counter, (ix1+1-ix0), (iy1+1-iy0), (iz1+1-iz0)])
endif

! Assign grid_field values
do i = 1, pt_counter
  ix1 = array(i)%ix(1)
  iy1 = array(i)%ix(2)
  iz1 = array(i)%ix(3)
  if (g_field%field_type == magnetic$) then
    g_field%ptr%pt(ix1, iy1, iz1)%B = array(i)%field(1:3)
  else
    g_field%ptr%pt(ix1, iy1, iz1)%E = array(i)%field(1:3)
    g_field%ptr%pt(ix1, iy1, iz1)%B = array(i)%field(4:6)
  endif
end do

! Clear temporary array

deallocate(array)

! Check if grid_field data has already been read in for another element.
! If so, save space by pointing to the existing grid.

call find_matching_fieldmap(g_field%ptr%file, ele, grid_field$, match_ele, im)
if (im > 0) then
  deallocate(g_field%ptr)
  g_field%ptr => match_ele%grid_field(im)%ptr
  g_field%ptr%n_link = g_field%ptr%n_link + 1        
endif

err_flag = .false.

!-----------------------------------------------------------
contains 

! subroutine parse_complex_component(complex_component, delim, err_flag2)
! looks for (x, y) or x followed by , or ) 
! returns complex field_component and next delim, which should be , or )

subroutine parse_complex_component(complex_component, delim, err_flag2)

character(1) delim, delim2
character(40) :: word, word2
integer  ix_word, ix_word2
logical delim_found, delim_found2, err_flag2
real(rp) x, y
complex(rp) complex_component

!

err_flag2 = .true.

! Expect "(" for complex, "," for real in the middle of the list, and ")" at the end of the list

call get_next_word (word, ix_word, ',()', delim, delim_found)
if (is_real(word)) then
  read (word, *) x
  complex_component = cmplx(x, 0.0_rp, rp)
else if (delim == '(') then
  call get_next_word (word, ix_word, ',', delim, delim_found)
  call get_next_word (word2, ix_word2, ')', delim2, delim_found2)
  ! Expect: real, real ) 
  if ((.not. is_real(word)) .or. (.not. is_real(word2)) &
    .or. (delim /= ',') .or. (delim2 /= ')') &
    .or. (.not. delim_found) .or. (.not. delim_found2)) then
    call parser_error ('BAD COMPLEX COMPONENT CONSTRUCT', &
                         'FOUND IN GRID_FIELD DEFINITION FOR ELEMENT: ' // ele%name)
    return
   end if
   !
   read (word, *) x
      read (word2, *) y
      complex_component = cmplx(x, y, rp)
      ! Look for "," or end of list ")"
   call get_next_word (word, ix_word, ',)', delim, delim_found)
end if

err_flag2 = .false.

end subroutine parse_complex_component

end subroutine parse_grid_field

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parse_taylor_field (grid, ele, lat, delim, delim_found, err_flag)
!
! Subroutine to parse a "taylor_field = {}" construct
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is private to bmad_parser_mod.
! This must read in:
! {type = ,
!    dr = , 
!    r0 = , 
!    pt(i,j,k) = ( (ex_re, ex_im), .... (bz_re, bz_im) ) 
!    .
!    .
!    . ) },
!-

subroutine parse_taylor_field (t_field, ele, lat, delim, delim_found, err_flag)

implicit none

type (taylor_field_struct), pointer :: t_field
type (ele_struct), target :: ele
type (ele_struct), pointer :: match_ele
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (em_taylor_term_struct), allocatable :: term(:)
type (em_taylor_term_struct), pointer :: tm
type (taylor_field_plane1_struct), allocatable :: plane(:)
type (taylor_field_plane1_struct), pointer :: t_plane

real(rp) coef

complex(rp), pointer :: c_ptr(:)

integer i, j, expn(2), nn, n, i_term, ib, ie, im, ix, family, ios
integer lb, i_out, ix_word

character(80) err_str
character(40) word, word2, name, attrib_name
character(1) delim, delim2

logical err_flag, delim_found

! Init

name = 'xxx'
err_flag = .true.

!

do

  ! Read attriubute
  call get_next_word (attrib_name, ix_word, '{}=,()', delim, delim_found, call_check = .true.)

  select case (attrib_name)

  case ('FIELD_SCALE')
    call evaluate_value (ele%name, t_field%field_scale, lat, delim, delim_found, err_flag, ',}')

  case ('R0')
    if (.not. equal_sign_here(ele, delim)) return
    if (.not. parse_real_list (lat, trim(ele%name) // ' GRID_FIELD', t_field%r0, .true.)) return
    if (.not. expect_one_of (',}', .false., ele, delim, delim_found)) return

  case ('DZ')
    call evaluate_value (ele%name, t_field%dz, lat, delim, delim_found, err_flag, ',}')

  case ('ELE_ANCHOR_PT', 'FIELD_TYPE', 'MASTER_PARAMETER')
    ! Expect "<component> = "
    if (delim /= '=') then
      call parser_error ('NO "=" SIGN FOUND AFTER ' // attrib_name,  &
                         'IN ELEMENT: ' // ele%name)
      return
    endif
    call get_next_word (word2, ix_word, ',}', delim, delim_found)

    !

    select case (attrib_name)
 
    case ('MASTER_PARAMETER')
      if (word2 == 'NONE') then
        ix = 0
      else
        ix = attribute_index(ele, word2)
        if (ix < 1 .or. ix > num_ele_attrib$) then
          call parser_error ('BAD NAME FOR "MASTER_PARAMETER = <NAME>" CONSTRUCT', &
                               'FOUND IN ELEMENT: ' // ele%name)
          return
        endif
      endif
      t_field%master_parameter = ix

    case ('ELE_ANCHOR_PT')
      call match_word(word2, anchor_pt_name(1:), t_field%ele_anchor_pt, can_abbreviate = .false., matched_name = name)
  
    case ('FIELD_TYPE')
      call match_word(word2, em_field_type_name(1:), t_field%field_type, can_abbreviate = .false., matched_name = name)
  
    end select

    !

    if (name == '') then
      call parser_error ('UNKNKOWN ' // trim(attrib_name) // ' VALUE:' // word2, &
                         'IN ELEMENT: ' // ele%name)
      return        
    endif      
      
  case ('CURVED_COORDS', 'CURVED_REF_FRAME')    ! 'curved_coords' is old style.
    if (.not. equal_sign_here(ele, delim)) return
    call get_next_word (word2, ix_word, ':,=()', delim, delim_found, .true.)
    t_field%curved_ref_frame = evaluate_logical (word2, ios)
    if (ios /= 0 .or. ix_word == 0) then
      call parser_error ('BAD TAYLOR_FIELD CURVED_REF_FRAME SETTING ' // word2, 'FOR: ' // ele%name)
    endif

  case ('CANONICAL_TRACKING')
    if (.not. equal_sign_here(ele, delim)) return
    call get_next_word (word2, ix_word, ':,=()', delim, delim_found, .true.)
    T_field%canonical_tracking = evaluate_logical (word2, ios)
    if (ios /= 0 .or. ix_word == 0) then
      call parser_error ('BAD TAYLOR_FIELD CANONICAL_TRACKING SETTING ' // word2, 'FOR: ' // ele%name)
    endif

  case ('PLANE')

    ! Expect "("
    if (.not.  expect_this ('(', .true., .false., 'AFTER "PLANE" IN TAYLOR_FIELD DEFINITION', ele, delim, delim_found)) return
    call parser_get_integer (ix, word, ix_word, delim, delim_found, err_flag, 'BAD PLANE INDEX IN TAYLOR_FIELD')
    if (.not.  expect_this (')={', .true., .false., 'AFTER "PLANE(IX" IN TAYLOR_FIELD DEFINITION', ele, delim, delim_found)) return

    if (allocated(t_field%ptr%plane)) then
      if (ix /= ubound(t_field%ptr%plane, 1) + 1) then
        call parser_error ('PLANE INDEX OUT OF ORDER IN TAYLOR_FIELD DEFININITION', 'FOR ELEMENT: ' // ele%name)
        return
      endif
      call move_alloc(t_field%ptr%plane, plane)
      lb = lbound(plane, 1)
      allocate (t_field%ptr%plane(lb:ix))
      t_field%ptr%plane(lb:ix-1) = plane
      deallocate(plane)
    else
      ! Set %file to be the last called file:<line_number>. 
      t_field%ptr%file = bp_com%line1_file_name
      allocate(t_field%ptr%plane(ix:ix))
    endif

    t_plane => t_field%ptr%plane(ix)
    allocate (t_plane%field(1)%term(0), t_plane%field(2)%term(0), t_plane%field(3)%term(0))

    do
      if (.not.  expect_this ('{', .false., .false., 'FOR PLANE TAYLOR TERM IN TAYLOR_FIELD DEFINITION', ele, delim, delim_found)) return
      call get_next_word (word, ix_word, ':{}=,()', delim, delim_found)
      if (.not.  expect_this (':', .true., .false., 'FOR PLANE TAYLOR TERM IN TAYLOR_FIELD DEFINITION', ele, delim, delim_found)) return
      call match_word (word, ['BX', 'BY', 'BZ'], i_out, .true., .false.)
      if (i_out < 1) then
        call parser_error ('BAD "OUT" COMPONENT: ' // word, 'IN TERM FOR TAYLOR_FIELD IN ELEMENT: ' // ele%name)
        return
      endif

      call evaluate_value (ele%name, coef, lat, delim, delim_found, err_flag, ',|');  if (err_flag) return
      delim2 = delim   ! Save
      expn = 0

      ! Need to check for "{x: 3.4 |}" case where there are no exponents.
      if (delim2 == '|') then
        call get_next_word (word, ix_word, '}', delim, delim_found)
        ! If there are exponents then rewind the get_next_word call.
        if (ix_word /= 0 .or. delim /= '}') then
          bp_com%parse_line = trim(word) // delim // bp_com%parse_line
          delim = delim2
        endif
      endif

      ! Parse exponents

      do j = 1, 100 
        if (delim == '}') exit
        call parser_get_integer (n, word, ix_word, delim, delim_found, err_flag, 'BAD EXPONENT');   if (err_flag) return
        if (.not. expect_one_of ('} ', .true., ele, delim, delim_found)) return
        if (delim2 == ',') then
          select case (j)
          case (2);      if( .not. expect_one_of ('}', .true., ele, delim, delim_found)) return
          case default;  if (.not. expect_one_of (' ', .true., ele, delim, delim_found)) return
          end select
          expn(j) = n
        else
          ! Where, for example, n = 34, must separate into 3 and 4.
          do
            if (word(1:1) == '') exit
            read (word(1:1), *) nn
            if (nn < 1 .or. nn > 2) then
              call parser_error ('BAD EXPONENT VALUE FOR TAYLOR ELEMENT IN TAYLOR_FIELD DEPFINITION IN: ' // ele%name)
              return
            endif
            expn(nn) = expn(nn) + 1
            word = word(2:)
          enddo
        endif
      enddo

      call add_em_taylor_term (t_plane%field(i_out), coef, expn)

      if (.not. expect_one_of ('},', .false., ele, delim, delim_found)) return

      if (delim == '}') then
        if (.not. expect_one_of ('},', .false., ele, delim, delim_found)) return
        exit
      endif

    enddo

  case default
    if (attrib_name == '') then
      call parser_error ('MANGLED GRID_FIELD DEFINITION FOR ELEMENT: ' // ele%name)
    else
      call parser_error ('UNKNOWN GRID_FIELD COMPONENT: ' // attrib_name, &
                         'FOR ELEMENT: ' // ele%name)
    endif
    return

  end select

  ! Possible "}" is end of mode
  if (delim == '}') exit

enddo

! Get final separator after grid construct.
 
if (.not. expect_one_of (', ', .false., ele, delim, delim_found)) return
err_flag = .false.

! Check if data has already been read in for another element.
! If so, save space by pointing to the data.

call find_matching_fieldmap(t_field%ptr%file, ele, taylor_field$, match_ele, im)
if (im > 0) then
  deallocate(t_field%ptr)
  t_field%ptr => match_ele%taylor_field(im)%ptr
  t_field%ptr%n_link = t_field%ptr%n_link + 1        
endif

end subroutine parse_taylor_field

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function parse_integer_list (err_str, lat, int_array, exact_size, open_delim, 
!           separator, close_delim, default_value) result (is_ok)
!
! Routine to parse a list of integers of the form:
!    open_delim integer_1 separator integer_2 . . . close_delim
! Example:   "(1.2, 2.3, 4.4, 8.5)"
! 
! Similar to parse_integer_list2 except does not use allocatable array.
! See parse_integer_list2 for more details
!-

function parse_integer_list (err_str, lat, int_array, exact_size, open_delim, &
                                    separator, close_delim, default_value) result (is_ok)

implicit none

type (lat_struct) lat

integer int_array(:)
integer, optional :: default_value
integer, allocatable :: vec(:)

integer num_found

character(*) err_str
character(*), optional :: open_delim, separator, close_delim

logical is_ok, exact_size

!

is_ok = .false.
if (.not. parse_integer_list2 (err_str, lat, vec, num_found, size(int_array), &
                          open_delim, separator, close_delim, default_value)) return

if (num_found > size(int_array) .or. (exact_size .and. num_found < size(int_array))) then
  call parser_error (err_str)
  return
endif

int_array(1:num_found) = vec(1:num_found)

is_ok = .true.

end function parse_integer_list

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function parse_integer_list2 (err_str, lat, int_array, num_found, num_expected, open_delim, 
!                               separator, close_delim, default_value) result (is_ok)
!
! Routine to parse a list of integers of the form
!    open_delim integer_1 separator integer_2 . . . close_delim
! Example:   (1, 2, 4, 8) 
!
! Input:
!  err_str       -- character(*): Error string to print if there is an error. 
!  lat           -- lat_struct: lattice
!  int_array(:)  -- Integer, allocatable: the array to be read in 
!
!   Optional: 
!   num_expected = 1     -- integer : number of expected arguments
!                             Used to initialize int_array
!   open_delim   = '('   -- character(1) : opening delimeter
!   separator    = ','   -- character(1) : separating character
!   close_delim  = ')'   -- character(1) : closing delimeter
!   default_value = 0    -- real(rp) : inital assignment of int_array elements
!
! Output:
!   is_ok                   -- logical: Set True if everything is ok
!   int_array(1:num_found) --integer(rp) : Array of values
!   num_found                  -- integer : number of elements
!-

function parse_integer_list2 (err_str, lat, int_array, num_found, num_expected, open_delim, &
                              separator, close_delim, default_value) result (is_ok)


type (lat_struct) lat

integer, allocatable :: int_array(:)
integer :: num_found
integer, optional :: num_expected, default_value
character(*) err_str
character(*), optional :: open_delim, close_delim, separator
logical is_ok

! Local
integer num_expect
character(1) delim, op_delim, cl_delim, sep
character(40) :: word
real(rp) rval
integer  ix_word
logical delim_found, err_flag

! Optional arguments

is_ok = .false.
num_expect = integer_option (1, num_expected)
op_delim = '('
cl_delim = ')'
sep      = ','
if (present(open_delim)) op_delim = open_delim
if (present(close_delim)) cl_delim = close_delim
if (present(separator)) sep = separator

! Expect op_delim
if (op_delim /= '') then
  call get_next_word (word, ix_word, op_delim, delim, delim_found)
  if ((word /= '') .or. (delim /= op_delim)) then
    call parser_error (err_str)
    return
  end if
endif

! Initial allocation
call re_allocate(int_array, num_expected, .false.)
int_array = integer_option(0, default_value)

! counter
num_found = 0

! Get integers
do 

  call evaluate_value ('BAD NUMBER IN: ' // err_str, rval, lat, delim, delim_found, err_flag, sep // cl_delim)
  if (err_flag) return
  if (abs(rval - nint(rval)) > 1d-10) then
    call parser_error ('BAD INTEGER NUMBER IN: ' // err_str)
    return
   end if    

  num_found = num_found + 1
  if (size(int_array) < num_found) then
    call re_allocate (int_array, 2*num_found, .false.)
    int_array(num_found:2*num_found) = integer_option(0, default_value)
  endif
  
  int_array(num_found) = nint(rval)
  
  ! Exit if cl_delim is found
  if (delim == cl_delim) exit
  
  ! Check separator
  if (delim /= sep) then
    call parser_error ('BAD SEPARATOR IN: ' // err_str)
    return  
  end if
  
end do

is_ok = .true.

end function parse_integer_list2

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function parse_real_list (lat, err_str, real_array, exact_size, open_delim, 
!           separator, close_delim, default_value) result (is_ok)
!
! Routine to parse a list of reals of the form:
!    open_delim real_1 separator real_2 . . . close_delim
! Example:   "(1.2, 2.3, 4.4, 8.5)"
! 
! Similar to parse_real_list2 except does not use allocatable array.
! See parse_real_list2 for more details
!-

function parse_real_list (lat, err_str, real_array, exact_size, open_delim, &
                          separator, close_delim, default_value) result (is_ok)

implicit none

type (lat_struct) lat

real(rp) real_array(:)
real(rp), optional :: default_value
real(rp), allocatable :: vec(:)

integer num_found

character(*) err_str
character(*), optional :: open_delim, separator, close_delim

logical is_ok, exact_size

!

is_ok = .false.
if (.not. parse_real_list2 (lat, err_str, vec, num_found, size(real_array), &
                          open_delim, separator, close_delim, default_value)) return

if (num_found > size(real_array) .or. (exact_size .and. num_found < size(real_array))) then
  call parser_error (err_str)
  return
endif

real_array(1:num_found) = vec(1:num_found)

is_ok = .true.

end function parse_real_list

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function parse_real_list2 (lat, err_str, real_array, num_found, num_expected, open_delim, 
!                            separator, close_delim, default_value) result (is_ok)
!
! Routine to parse a list of reals of the form:
!    open_delim real_1 separator real_2 . . . close_delim
! Example:   "(1.2, 2.3, 4.4, 8.5)"
!
! Input:
!  lat        -- lat_struct: lattice
!  err_str    -- character(*): Error string to print if there is an error. 
!  real_array -- real(rp), allocatable: the array to be read in 
!
!   Optional: 
!   num_expected = 1     -- integer : number of expected arguments
!                             Used to initialize real_array
!   open_delim   = '('   -- character(1) : opening delimeter
!   separator    = ','   -- character(1) : separating character
!   close_delim  = ')'   -- character(1) : closing delimeter
!   default_value = 0.0_rp -- real(rp) : inital assignment of real_array elements
!
! Output:
!   is_ok                   -- logical: Set True if everything is ok
!   real_array(1:num_found) -- real(rp) : Array of values
!   num_found               -- integer : number of elements
!-

function parse_real_list2 (lat, err_str, real_array, num_found, num_expected, open_delim, &
          separator, close_delim, default_value) result (is_ok)

! Arguments

type (lat_struct) lat

real(rp), allocatable :: real_array(:)
integer :: num_found
integer, optional :: num_expected
character(*) err_str
character(1), optional :: open_delim, close_delim, separator
logical is_ok
real(rp), optional :: default_value

! Local
integer num_expect
character(1) delim, op_delim, cl_delim, sep
character(40) :: word
integer  ix_word
logical delim_found, err_flag
real(rp) :: default_val, value

! Optional arguments

is_ok = .false.
num_expect = integer_option(1, num_expected)
default_val = real_option(0.0_rp, default_value)

op_delim = '('
cl_delim = ')'
sep      = ','
if (present(open_delim)) op_delim = open_delim
if (present(close_delim)) cl_delim = close_delim
if (present(separator)) sep = separator

! Expect op_delim
call get_next_word (word, ix_word, op_delim, delim, delim_found)
if ((word /= '') .or. (delim /= op_delim)) then
  call parser_error ('BAD OPENING DELIMITER IN: ' // err_str)
  return
end if

! Initial allocation
call re_allocate(real_array, num_expected, .false.)
real_array = default_val

! Get reals

num_found = 0

do 
  call evaluate_value ('BAD REAL NUMBER IN: ' // err_str, value, lat, delim, delim_found, err_flag, sep // cl_delim)
  if (err_flag) return
  ! real is found
  num_found = num_found + 1
  ! reallocate if needed  
  if (size(real_array) < num_found) then
    call re_allocate (real_array, 2*num_found, .false.)
    real_array(num_found:2*num_found) = default_val
  endif

  ! Set value
   real_array(num_found) = value
  
  ! Exit if cl_delim is found
  if (delim == cl_delim) exit
  
  ! Check separator
  if (delim /= sep) then
    call parser_error ('BAD SEPARATOR IN: ' // err_str)
    return  
  end if

end do

is_ok = .true.

end function parse_real_list2

!----------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------
!+
! Subroutine parser_set_spin (bs_ele, orbit)
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!- 

subroutine parser_set_spin (bs_ele, orbit)

use spin_mod

type (ele_struct) bs_ele    ! beam_start element
type (coord_struct) orbit   
type (spin_polar_struct) :: polar

real(rp) vec(3)

!

vec = [bs_ele%value(spin_x$), bs_ele%value(spin_y$), bs_ele%value(spin_z$)]
polar = spin_polar_struct(bs_ele%value(spinor_polarization$), bs_ele%value(spinor_theta$), &
                          bs_ele%value(spinor_phi$), bs_ele%value(spinor_xi$))

if (any(vec /= [0, 0, 1])) then
  if (polar%polarization /= 1 .or. polar%theta /= 0 .or. polar%phi /= 0 .or. polar%xi /= 0) &
              call parser_error ('ERROR SETTING BEAM_START. BOTH SPIN_X/Y/Z AND SPINOR_XXX QUANTITIES SET!')
  polar = vec_to_polar (vec)
endif

orbit%spin = polar_to_vec (polar)

end subroutine parser_set_spin

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

subroutine parser_get_integer (int_val, word, ix_word, delim, delim_found, err, str1, str2)

integer int_val, ix_word
logical err
character(*) word, delim
character(*), optional :: str1, str2
logical delim_found

!

call get_next_word (word, ix_word, ':,=(){} ', delim, delim_found, .true.)
if (.not. is_integer(word) ) then
  if (present(str1)) then
    call parser_error (str1, str2)
  else
    call parser_error ('INTEGER EXPECTED, I DO NOT UNDERSTAND: ' // word)
  endif
  err = .true.

else
  read (word, *) int_val
  err = .false.
endif

end subroutine parser_get_integer

!--------------------------------------------------------------------------------------

subroutine parser_get_logical (ele, attrib_name, this_logic, word, ix_word, delim, delim_found, err)

type (ele_struct) ele
character(*) attrib_name, word, delim
integer ix_word, ios
logical this_logic, err
logical delim_found

!

call get_next_word (word, ix_word, ':,=()', delim, delim_found, .true.)
this_logic = evaluate_logical (word, ios)
if (ios /= 0 .or. ix_word == 0) then
  call parser_error ('BAD "' // trim(attrib_name) // '" SWITCH FOR: ' // ele%name, 'I DO NOT UNDERSTAND: ' // word)
  err = .true.
else
  err = .false.
endif

end subroutine parser_get_logical

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
! delim_list    -- character(*): String of tokens expected. 
! check_delim   -- logical: If True then use "delim" character as first token to check.
!                    A blank character indicates end of command is expected.
! call_check    -- Logical: If True then check for 'call::<filename>' construct.

function expect_this (expecting, check_delim, call_check, err_str, ele, delim, delim_found) result (is_ok)

implicit none

type (ele_struct) ele
character(*) expecting, err_str
character(1) delim
character(40) word
logical is_ok, delim_found, check_delim, call_check
integer ix, ix_word

!

is_ok = .false.

do ix = 1, len(expecting)
  if (ix == 1 .and. check_delim) then
    word = ''
  else
    call get_next_word (word, ix_word, expecting(ix:ix), delim, delim_found, call_check = call_check)
  endif

  if (expecting(ix:ix) == ' ') then
    if (delim_found .or. word /= '') then
      call parser_error ('EXTRA STUFF ON LINE.', err_str)
      return
    endif
  elseif (delim /= expecting(ix:ix) .or. word /= '') then
    call parser_error ('NO "' // expecting // '" FOUND ' // err_str, 'FOR ELEMENT: ' // ele%name)
    return
  endif
enddo

is_ok = .true.

end function expect_this

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

subroutine get_switch (name, name_list, switch, err, ele, delim, delim_found)

type (ele_struct) :: ele
type (all_pointer_struct), allocatable :: ptr(:)

character(*) name, name_list(:)
character(1) delim
character(60) word
character(40) ele_name, attrib_name 
character(:), allocatable :: line
integer i, this_switch, switch, ix_word, ixp, ixp2
logical err, delim_found

!

err = .true.
call get_next_word (word, ix_word, ':,=(){}', delim, delim_found, .true.)

!

if (name == 'FIELD_CALC') then
  if (word == 'GRID' .or. word == 'MAP') then
    
    call parser_error ('FIELD_CALC = ' // word, 'HAS BEEN CHANGED TO FIELD_CALC = FIELDMAP', &
                       'Program will execute as normal...', &
                       '[But eventually this warning will be converted to an error. You have been warned!]', level = s_warn$)
    word = 'FIELDMAP'
  endif
endif

! If word is something like "q1[tracking_method]" then need retrieve this value.

ixp = index(word, '[')
if (ixp == 0) then
  call match_word (word, name_list, this_switch, can_abbreviate = .false.)
  if (this_switch < 1) then
    do  i = 1, size(name_list)
      if (name_list(i) == str_garbage$) cycle
      line = line // ' ' // trim(name_list(i)) // ','
    enddo
    line = line(1:len_trim(line)-1)
    call parser_error ('BAD "' // trim(name) // '" SWITCH FOR: ' // ele%name, 'I DO NOT UNDERSTAND: ' // word, &
                       'POSSIBILITIES ARE: ' // line)
    return
  else
    switch = this_switch
  endif

else
  ixp2 = len_trim(word)
  if (word(ixp2:ixp2) /= ']' .or. word(ixp+1:ixp2-1) /= name) then
    call parser_error ('BAD "' // trim(name) // '" SWITCH FOR: ' // ele%name, 'I DO NOT UNDERSTAND: ' // word)
    return
  endif

  ele_name = word(:ixp-1)
  attrib_name = word(ixp+1:ixp2-1)
  call pointers_to_attribute (ele%branch%lat, ele_name, attrib_name, .true., ptr, err, .false.)
  if (size(ptr) == 0) then
    call parser_error ('NO ELEMENT FOUND TO EVALUATE: ' // word, 'EVALUATING SWITCH IN ELEMENT: ' // ele%name)
    return
  endif
  if (size(ptr) > 1) then
    call parser_error ('MULTIPLE ELEMENTS FOUND FOR EVALUATING: ' // word, 'EVALUATING SWITCH IN ELEMENT: ' // ele%name)
    return
  endif
  if (associated(ptr(1)%i)) switch = ptr(1)%i
  if (associated(ptr(1)%r)) switch = nint(ptr(1)%r)
endif

err = .false.

end subroutine get_switch

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

! delim_list -- character(*): List of valid tokens. If list contains a space character
!                 then no token (indicating the end of the command) is a valid possibility.
!

function expect_one_of (delim_list, check_delim, ele, delim, delim_found) result (is_ok)

type (ele_struct) ele
integer ix_word
character(*) delim_list
character(1) delim
character(40) word
logical check_delim, delim_found, is_ok, must_have_delim

!

is_ok = .false.
must_have_delim = (index(delim_list, ' ') == 0)

if (check_delim) then
  if ((must_have_delim .and. .not. delim_found) .or. &
      (delim /= '' .and. index(delim_list, delim) == 0)) then
    call parser_error ('BAD DELIMITOR', 'FOR ELEMENT: ' // ele%name)
    return
  endif

else
  call get_next_word (word, ix_word, '{}=,()[]', delim, delim_found)
  if (word /= '' .or. (must_have_delim .and. .not. delim_found) .or. &
      (delim /= '' .and. index(delim_list, delim) == 0)) then
    call parser_error ('BAD DELIMITOR', 'FOR ELEMENT: ' // ele%name)
    return
  endif
endif

is_ok = .true.

end function expect_one_of

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

function equal_sign_here(ele, delim) result (is_here)

type (ele_struct) ele
character(40) word
character(1) delim
logical is_here

! Expect "<component> = "

is_here = .false.

if (delim /= '=') then
call parser_error ('NO "=" SIGN FOUND AFTER: ' // word,  &
                   'IN GRID_FIELD STRUCTURE IN ELEMENT: ' // ele%name)
  return
endif

is_here = .true.

end function equal_sign_here

end module
