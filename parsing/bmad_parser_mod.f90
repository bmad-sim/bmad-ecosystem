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
use add_superimpose_mod
use track1_mod

private parse_rf_grid, parse_rf_map

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
  integer ix                      ! current index of element in %ELE
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
  character(200) logical_name
  character(200) :: full_name = ''
  character(200) :: dir = './'
  character(200) parse_line_saved
  integer i_line
  integer f_unit
  logical inline_call_active
end type

!-----------------------------------------------------------

integer, parameter, private :: plus$ = 1, minus$ = 2, times$ = 3, divide$ = 4
integer, parameter, private :: l_parens$ = 5, r_parens$ = 6, power$ = 7
integer, parameter, private :: unary_minus$ = 8, unary_plus$ = 9, no_delim$ = 10
integer, parameter, private :: sin$ = 11, cos$ = 12, tan$ = 13
integer, parameter, private :: asin$ = 14, acos$ = 15, atan$ = 16, abs$ = 17, sqrt$ = 18
integer, parameter, private :: log$ = 19, exp$ = 20, ran$ = 21, ran_gauss$ = 22
integer, parameter, private :: atan2$ = 23, factorial$ = 24, numeric$ = 100

integer, parameter, private :: eval_level(24) = [1, 1, 2, 2, 0, 0, 4, 3, 3, -1, &
                            9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9]

type eval_stack_struct
  integer type
  real(rp) value
end type

! structure for holding the control names and pointers for superimpose and overlay elements

type parser_ele_struct
  character(40) ref_name
  character(40) :: ele_name = ''               ! For patch.
  character(40), allocatable :: name(:)        ! For overlays and groups
  character(40), allocatable :: attrib_name(:) ! For overlays and groups
  character(200) lat_file                      ! File where element was defined.
  real(rp), allocatable :: coef(:)             ! For overlays and groups
  real(rp) s
  integer ix_line_in_file    ! Line in file where element was defined.
  integer ix_count
  integer ele_pt, ref_pt
  integer indexx
  logical create_jumbo_slave
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
  type (ele_struct), pointer :: param_ele, beam_ele, beam_start_ele, root_branch_ele
  type (stack_file_struct), pointer :: current_file
  type (stack_file_struct), pointer :: calling_file
  type (lat_struct), pointer :: old_lat
  type (bp_var_struct), allocatable :: var(:)   ! variable name
  type (extra_parsing_info_struct) extra
  integer num_lat_files               ! Number of files opened
  integer ivar_tot, ivar_init
  character(200), allocatable :: lat_file_names(:) ! List of all files used to create lat
  character(n_parse_line) parse_line
  character(n_parse_line) input_line1          ! For debug messages
  character(n_parse_line) input_line2          ! For debug messages
  character(40) parser_name
  character(200) :: dirs(2) 
  logical :: bmad_parser_calling = .false.              ! used for expand_lattice
  logical error_flag                           ! Set True on error
  logical fatal_error_flag                     ! Set True on fatal (must abort now) error 
  logical input_line_meaningful
  logical do_superimpose
  logical write_digested      ! For bmad_parser
  logical write_digested2     ! For bmad_parser2
  logical input_from_file     ! Input is from a lattice file?
  logical inline_call_active
  logical :: always_parse = .false. ! For debugging to force parsing
  logical :: print_err = .true.  ! Print error messages?
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
type (ele_struct), target, save ::  ele0
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: bele
type (wig_term_struct), pointer :: wig_term(:)
type (real_pointer_struct), allocatable :: r_ptrs(:)
type (wall3d_section_struct), pointer :: section
type (wall3d_vertex_struct), pointer :: v_ptr
type (em_field_mode_struct), pointer :: em_modes(:)
type (em_field_mode_struct), pointer :: em_mode
type (photon_surface_struct), pointer :: surf

real(rp) kx, ky, kz, tol, value, coef, r_vec(10)
real(rp), pointer :: r_ptr

integer i, j, n, ix_word, how, ix_word1, ix_word2, ios, ix, i_out, ix_coef, switch
integer expn(6), ix_attrib, i_section, ix_v, ix_sec, i_mode, i_term, ib, ie, im
integer x_bounds(2), y_bounds(2), i_vec(2)

character(40) :: word, str_ix, attrib_word, word2
character(1) delim, delim1, delim2
character(80) str, err_str, line

logical delim_found, err_flag, logic, set_done, end_of_file, do_evaluate, wild_key0, err_flag2
logical, optional :: check_free, wild_and_key0

! Get next WORD.
! If an overlay or group element then word is just an attribute to control
! [except for a "GROUP[COMMAND] = 0.343" redef construct]

err_flag = .true.  ! assume the worst
call get_next_word (word, ix_word, ':, =()', delim, delim_found, call_check = .true.)

! taylor

wild_key0 = logic_option(.false., wild_and_key0)

if (ele%key == taylor$ .and. word(1:1) == '{') then

  word = word(2:)             ! strip off '{'
  read (word, *, iostat = ios) i_out
  if (delim /= ':' .or. ix_word == 0 .or. ios /= 0) then
    call parser_error ('BAD "OUT" IN TERM FOR TAYLOR ELEMENT: ' // ele%name, 'CANNOT PARSE: ' // str)
    return
  endif

  call evaluate_value (str, coef, lat, delim, delim_found, err_flag)
  if (err_flag) return

  call get_next_word (line, ix_word, '},', delim, delim_found)
  read (line, *, iostat = ios) expn
  if (delim /= '}' .or. ix_word == 0 .or. ios /= 0) then
    call parser_error ('BAD "EXPONENT" IN TERM FOR TAYLOR ELEMENT: ' // ele%name, 'CANNOT PARSE: ' // str)
    return
  endif

  call add_this_taylor_term (ele, i_out, coef, expn)
  call get_next_word (word, ix_word, '},', delim, delim_found)

  if (ix_word /= 0 .or. (delim_found .and. delim /= ',')) then
    call parser_error ('BAD TERM ENDING FOR TAYLOR ELEMENT: ' // ele%name, 'CANNOT PARSE: ' // str)
    return
  endif

  return
endif

! overlay

if (ele%key == overlay$) then
  i = attribute_index(ele, word)       ! general attribute search

  if (i == type$ .or. i == alias$ .or. i == descrip$) then
    call bmad_parser_type_get (ele, word, delim, delim_found)

  else

    if (i < 1) then
      if (wild_key0) then
        err_flag = .false.
        return
      endif
      call parser_error ('BAD OVERLAY ATTRIBUTE: ' // word, 'FOR: ' // ele%name)
      return
    endif

    if (how == def$) then
      ele%ix_value = i
      ele%component_name = word
    endif

    value = 0
    if (delim == '=') then  ! value
      call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag)
      if (err_flag) return
    endif

    if (ele%ix_value /= i) then
      if (wild_key0) then
        err_flag = .false.
        return
      endif
      call parser_error ('BAD OVERLAY ATTRIBUTE SET FOR: ' // ele%name, &
            'YOU ARE TRYING TO SET: ' // word, &
            'BUT YOU SHOULD BE SETTING: ' // ele%component_name)
      return
    endif

    call pointer_to_indexed_attribute (ele, i, .true., r_ptr, err_flag, .true.)
    r_ptr = value

    if (attrib_free_problem(word)) return

  endif

  err_flag = .false.
  return
endif

! group...

if (ele%key == group$) then

  if (how == def$) then
    ele0%key = overlay$
    i = attribute_index(ele0, word)       ! general attribute search
  else   ! how == redef$
    i = attribute_index(ele, word)
  endif

  if (i < 1) then
    call parser_error ('BAD GROUP ATTRIBUTE: ' // word, 'FOR: ' // ele%name)
    return
  endif

  if (i == type$ .or. i == alias$) then
    call bmad_parser_type_get (ele, word, delim, delim_found)
  else
    if (how == def$) then
      ele%ix_value = i
      ele%component_name = word
    endif
    if (delim == '=') then  ! value
      call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag)
      if (err_flag) return
      if (how == def$) then
        ele%value(command$) = value
      elseif (i > num_ele_attrib$) then
        call parser_error ('BAD GROUP ATTRIBUTE: ' // word, 'FOR: ' // ele%name)
        return
      else
        ele%value(i) = value
      endif
    elseif (how == redef$) then
      call parser_error ('NO VALUE GIVEN FOR ATTRIBUTE FOR: ' // ele%name)
    endif
    if (attrib_free_problem(word)) return
  endif

  err_flag = .false.
  return

endif

! beam_start element can have attributes that are not part of the element so
! Need to use pointers_to_attribute.

if (ele%key == def_beam_start$) then
  call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag) 
  if (err_flag) return
  call pointers_to_attribute (lat, ele%name, word, .false., r_ptrs, err_flag, .false.)
  if (err_flag .or. size(r_ptrs) /= 1) then
    call parser_error ('BAD ATTRIBUTE: ' // word, 'FOR ELEMENT: ' // ele%name)
    return
  endif

  r_ptrs(1)%r = value

  return
endif

! Long-range wake

if (word == 'LR' .and. delim == '(') then

  call get_next_word (word, ix_word, '=', delim, delim_found)
  if (.not. delim_found) then
    call parser_error ('NO "=" SIGN FOUND', 'FOR ELEMENT: ' // ele%name)
    return
  endif
  call pointer_to_attribute (ele, 'LR(' // word, .false., r_ptr, err_flag, .false.)
  if (err_flag) then
    call parser_error ('BAD ATTRIBUTE: ' // word, 'FOR ELEMENT: ' // ele%name)
    return
  endif
  call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag)
  r_ptr = value
  return
endif

! "WALL." redef

if (word(1:5) == 'WALL.') then

  select case (word(6:))

  ! Section redef

  case ('SECTION')

    if (delim /= '(') then
      call parser_error ('MALFORMED WALL COMPONENT REDEF IN ELEMENT: ' // ele%name)
      return
    endif

    ix_sec = evaluate_array_index (err_flag, ')', word2, '(=', delim)
    if (err_flag .or. .not. associated(ele%wall3d) .or. ix_sec < 0 .or. &
                                                  ix_sec > size(ele%wall3d%section)) then
      call parser_error('BAD ' // trim(word) // ' INDEX', 'FOR ELEMENT: ' // ele%name)
      return
    endif
    section => ele%wall3d%section(ix_sec)

    if (word2 == '.S' .and. delim == '=') then
      r_ptr => section%s

    elseif (word2 == '.DR_DS' .and. delim == '=') then
      r_ptr => section%dr_ds

    elseif (word2 == '.V' .and. delim == '(') then
      ix_v = evaluate_array_index (err_flag, ')', word, '=', delim)
      if (err_flag .or. ix_v < 0 .or. ix_v > size(section%v)) then
        call parser_error('BAD VERTEX INDEX',  'FOR ELEMENT: ' // ele%name)
        return
      endif
      v_ptr => section%v(ix_v)

      select case (word)
      case ('.X');        r_ptr => v_ptr%x
      case ('.Y');        r_ptr => v_ptr%y
      case ('.RADIUS_X'); r_ptr => v_ptr%radius_x
      case ('.RADIUS_Y'); r_ptr => v_ptr%radius_y
      case ('.TILT');     r_ptr => v_ptr%tilt
      case default
        call parser_error('BAD WALL SECTION VERTEX COMPONENT: ' // word, 'FOR ELEMENT: ' // ele%name)
        return
      end select
    else
      call parser_error('BAD WALL SECTION COMPONENT: ' // word2, 'FOR ELEMENT: ' // ele%name)
      return
    endif

    call evaluate_value (trim(ele%name), r_ptr, lat, delim, delim_found, err_flag)

  ! Not recognized

  case default
    call parser_error ('BAD WALL COMPONENT REDEF: ' // word, 'IN ELEMENT: ' // ele%name)
  end select

  return
endif

! if not an overlay then see if it is an ordinary attribute.
! if not an ordinary attribute then might be a superimpose switch

if (ix_word == 0) then  ! no word
  call parser_error  ('"," NOT FOLLOWED BY ATTRIBUTE NAME FOR: ' // ele%name)
  return
endif

if (word == 'TILT' .and. (ele%key == sbend$ .or. ele%key == rbend$)) then
  call parser_error ('BENDS HAVE A "REF_TILT" ATTRIBUTE BUT NOT A "TILT" ATTRIBUTE.')
endif

if (word == 'DPHI0') then
  call parser_error ('THE ATTRIBUTE NAME "DPHI0" HAS BEEN CHANGED TO "PHI0_MULTIPASS"', &
                     'PLEASE MAKE THE CHANGE IN THE LATTICE FILE.', &
                     '[THIS IS A WARNING ONLY. THIS PROGRAM WILL RUN NORMALLY]', warn_only = .true.)
  word = 'PHI0_MULTIPASS'
endif

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
  if (associated (ele%wall3d)) then
    call parser_error ('MULTIPLE WALL DEFINITIONS FOR ELEMENT: ' // ele%name)
    return
  endif

  i_section = 0
  allocate (ele%wall3d)
  if (.not. expect_this ('={', .true., .true., 'AFTER "WALL"')) return

  ! Loop wall3d_struct components.

  wall3d_loop: do    

    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

    ! Possible "}" is end of wall 
    if (delim /= '}' .and. word == '') exit
    if (.not. expect_this ('=', .true., .false., 'AFTER ' // trim(word) // ' IN WALL CONSTRUCT')) return

    select case (word)

    case ('OPAQUE_MATERIAL') 
      call bmad_parser_type_get (ele, word, delim, delim_found, str_out = ele%wall3d%opaque_material)

    case ('CLEAR_MATERIAL') 
      call bmad_parser_type_get (ele, word, delim, delim_found, str_out = ele%wall3d%clear_material)

    case ('THICKNESS') 
      call evaluate_value (err_str, ele%wall3d%thickness, lat, delim, delim_found, err_flag, ',}')
      if (err_flag) return

    case ('ELE_ANCHOR_PT')
      call get_switch ('WALL ELE_ANCHOR_PT', anchor_pt_name(1:), ele%wall3d%ele_anchor_pt, err_flag2)
      if (err_flag2) return

    case ('SUPERIMPOSE')
      call get_logical ('WALL SUPERIMPOSE', ele%wall3d%superimpose, err_flag2); if (err_flag2) return

    ! Must be "section = {"

    case ('SECTION')

      ! Read in section

      if (.not. expect_this ('{', .false., .true., 'AFTER "SECTION =" IN WALL CONSTRUCT')) return

      i_section = i_section + 1
      ix_v = 0
      call re_allocate (ele%wall3d%section, i_section)
      section => ele%wall3d%section(i_section)

      wall3d_section_loop: do

        call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

        ! Possible "}" is end of wall 
        if (delim /= '}' .and. word == '') exit
        if (word == 'V') then
          if (.not. expect_this ('(', .true., .false., 'AFTER ' // trim(word) // ' IN WALL CONSTRUCT')) return
        else
          if (.not. expect_this ('=', .true., .false., 'AFTER ' // trim(word) // ' IN WALL CONSTRUCT')) return
        endif

        select case (word)

        case ('TYPE') 
          call get_switch ('WALL TYPE', wall3d_section_type_name(1:), section%type, err_flag2)
          if (err_flag2) return

        case ('MATERIAL') 
          call bmad_parser_type_get (ele, word, delim, delim_found, str_out = section%material)

        case ('THICKNESS')
          call evaluate_value (trim(ele%name), section%thickness, lat, delim, delim_found, err_flag, ',}')
          if (err_flag) return
          if (ele%key == capillary$) ele%value(l$) = section%s

        case ('S')
          call evaluate_value (trim(ele%name), section%s, lat, delim, delim_found, err_flag, ',}')
          if (err_flag) return
          if (ele%key == capillary$) ele%value(l$) = section%s

        case ('DR_DS') 
          call evaluate_value (trim(ele%name), section%dr_ds, lat, delim, delim_found, err_flag, ',}')
          if (err_flag) return
                  
        case ('X0') 
          call evaluate_value (trim(ele%name), section%x0, lat, delim, delim_found, err_flag, ',}')
          if (err_flag) return
                  
        case ('Y0') 
          call evaluate_value (trim(ele%name), section%y0, lat, delim, delim_found, err_flag, ',}')
          if (err_flag) return

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

          if (.not. expect_this (')={', .true., .false., 'AFTER "V(n)" IN WALL CONSTRUCT')) return

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

        if (.not. expect_either (',}', .true.)) return
        if (delim == '}') then
          if (.not. expect_either(',}', .false.)) return
          exit
        endif
      enddo wall3d_section_loop

    case default
      call parser_error ('WALL COMPONENT NOT RECOGNIZED: ' // word, 'FOR ELEMENT: ' // ele%name)
      return
    end select   ! wall components

    if (.not. expect_either (',}', .true.)) return
    if (delim == '}') exit

  enddo wall3d_loop

  ! Next thing on line should be either a "," or end-of-line

  logic = expect_either(', ', .false.)
  return

endif

!-------------------------------
! Reflecting Surface

if (attrib_word == 'SURFACE') then
  surf => ele%photon%surface

  if (.not. expect_this ('={', .true., .true., 'AFTER "SURFACE"')) return

  surface_loop: do

    ! Expect "GRID ={" or "TYPE ="

    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

    select case (word)

    case ('GRID')

      if (.not. expect_this ('={', .true., .true., 'AFTER "GRID"')) return
      x_bounds = int_garbage$; y_bounds = int_garbage$

      do
        call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
        if (word /= 'PT') then
          if (.not. expect_this ('=', .true., .false., 'AFTER ' // trim(word) // ' IN SURFACE CONSTRUCT')) return
        endif

        select case (word)
        case ('DR')
          if (.not. parse_real_list (lat, trim(ele%name) // ' GRID DR', surf%grid%dr, .true.)) return

        case ('R0')
          if (.not. parse_real_list (lat, trim(ele%name) // ' GRID R0', surf%grid%r0, .true.)) return

        case ('X_BOUNDS', 'Y_BOUNDS')
          if (.not. parse_integer_list (trim(ele%name) // ' GRID ' // trim(word), i_vec, .true.)) return
          if (word == 'X_BOUNDS') x_bounds = i_vec
          if (word == 'Y_BOUNDS') y_bounds = i_vec

          if (any(x_bounds /= int_garbage$) .and. any(y_bounds /= int_garbage$)) then
            if (any(x_bounds == int_garbage$) .or. any(y_bounds == int_garbage$) .or. &
                x_bounds(1) > x_bounds(2) .or. y_bounds(1) > y_bounds(2)) then
              call parser_error ('SURFACE GRID X/Y_BOUNDS NOT PROPERLY SET', trim(ele%name))
              return
            endif
            if (allocated (surf%grid%pt)) deallocate (surf%grid%pt)
            allocate (surf%grid%pt(x_bounds(1):x_bounds(2), y_bounds(1):y_bounds(2)))
          endif

        case ('PT')
          bp_com%parse_line = delim // bp_com%parse_line
          if (.not. parse_integer_list (trim(ele%name) // ' GRID PT', i_vec, .true.)) return
          if (.not. allocated(surf%grid%pt)) then
            call parser_error ('SURFACE PT_MAX MISSING', 'FOR: ' // ele%name)
            return
          endif
          if (any(i_vec < 0) .or. any(i_vec > ubound(surf%grid%pt))) then
            call parser_error ('SURFACE PT(I,J) INDEX OUT OF BOUNDS', 'FOR: ' // ele%name)
            return
          endif
          if (.not. expect_this ('=', .false., .false., 'GRID PT')) return
          if (.not. parse_real_list (lat, trim(ele%name) // ' GRID PT', r_vec(1:4), .true.)) return
          surf%grid%pt(i_vec(1), i_vec(2)) = surface_grid_pt_struct(r_vec(1), r_vec(2), r_vec(3), r_vec(4))

        case ('TYPE')
          call get_switch ('SURFACE GRID TYPE', surface_grid_type_name(1:), surf%grid%type, err_flag2)
          if (err_flag2) return
          bp_com%parse_line = delim // bp_com%parse_line

        case default
          call parser_error ('GRID COMPONENT NOT RECOGNIZED: ' // word, 'FOR ELEMENT: ' // ele%name)
          return
        end select

        if (.not. expect_either (',}', .false.)) return
        if (delim == '}') then
          if (.not. expect_either (',}', .false.)) return
          exit
        endif

      enddo

    case default
      call parser_error ('UNKNOWN SURFACE COMPONENT: ' // word2, 'FOR: ' // ele%name)
      return
    end select

    if (delim == '}') exit

  enddo surface_loop

  if (.not. expect_either(', ', .false.)) return
  err_flag = .false.
  return

endif

!-------------------------------
! RF field

if (attrib_word == 'FIELD') then

  if (.not. expect_this ('={', .true., .true., 'AFTER "FIELD"')) return

  ! Loop over all modes

  do

    if (associated(ele%em_field)) then
      i_mode = size(ele%em_field%mode) + 1
      em_modes => ele%em_field%mode
      nullify(ele%em_field)
      call init_em_field (ele%em_field, i_mode)
      ele%em_field%mode(1:i_mode-1) = em_modes 
      deallocate(em_modes)
    else
      call init_em_field (ele%em_field, 1)
      i_mode = 1
    endif

    em_mode => ele%em_field%mode(i_mode)
    if (ele%key == lcavity$ .or. ele%key == rfcavity$) em_mode%harmonic = 1 ! Default

    ! Expect "MODE = {" or "MODE_TO_AUTOSCALE ="

    call get_next_word (word2, ix_word, '{}=,()', delim, delim_found)
    if (word2 /= 'MODE' .and. word2 /= 'MODE_TO_AUTOSCALE') then
      call parser_error ('EXPECTED "MODE" OR , "MODE_TO_AUTOSCALE" IN FIELD DEFINITION BUT FOUND: ' // word2, &
                         'FOR ELEMENT: ' // ele%name)
      return
    endif 

    ! MODE_TO_AUTOSCALE 

    if (word2 == 'MODE_TO_AUTOSCALE') then
      if (.not. expect_this ('=', .true., .true., 'AFTER "MODE_TO_AUTOSCALE"')) return
      call get_integer (ele%em_field%mode_to_autoscale, .false., err_flag, 'BAD MODE_TO_AUTOSCALE VALUE')
      if (err_flag) return
    endif

    ! MODE = { 

    if (.not. expect_this ('={', .true., .true., 'AFTER "MODE"')) return
    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

    do

      ! Expect "<component> = "

      if (delim /= '=') then
        call parser_error ('NO "=" SIGN FOUND IN MODE DEFINITION', 'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
        return
      endif

      do_evaluate = .true.
      word2 = word

      select case (word2)

      case ('F_DAMP');         r_ptr => em_mode%f_damp
      case ('PHI0_REF');       r_ptr => em_mode%phi0_ref
      case ('STORED_ENERGY');  r_ptr => em_mode%stored_energy
      case ('PHI0_AZIMUTH');   r_ptr => em_mode%phi0_azimuth
      case ('FIELD_SCALE');    r_ptr => em_mode%field_scale
      case ('DPHI0_REF')
        r_ptr => em_mode%phi0_ref  ! For backwards compatibility
        call parser_error ('THE ATTRIBUTE NAME "DPHI0_REF" HAS BEEN CHANGED TO "PHI0_REF"', &
                           'PLEASE MAKE THE CHANGE IN THE LATTICE FILE.', &
                           '[THIS IS A WARNING ONLY. THIS PROGRAM WILL RUN NORMALLY]', warn_only = .true.)

      case ('GRID') 
        call parse_rf_grid(em_mode%grid, ele, lat, delim, delim_found, err_flag)
        if (err_flag) return
        do_evaluate = .false.

      case ('MAP') 
        call parse_rf_map(em_mode%map, ele, lat, delim, delim_found, err_flag)
        if (err_flag) return
        do_evaluate = .false.

      case ('M', 'HARMONIC')

        select case (word2)
        case ('M');        call get_integer (em_mode%m, .false., err_flag, &
                                      'BAD "FIELD MODE M = <INT>" CONSTRUCT', 'IN ELEMENT: ' // ele%name)
        case ('HARMONIC'); call get_integer (em_mode%harmonic, .false., err_flag, &
                                      'BAD "FIELD MODE HARMONIC = <INT>" CONSTRUCT', 'IN ELEMENT: ' // ele%name)
        end select

        if (err_flag) return
        do_evaluate = .false.

      case ('MASTER_SCALE')
        call get_next_word (word, ix_word, ',}', delim, delim_found)
        if (word == 'NONE') then
          ix = 0
        else
          ix = attribute_index(ele, word)
          if (ix < 1 .or. ix > num_ele_attrib$) then
            call parser_error ('BAD NAME FOR "MASTER_SCALE = <NAME>" CONSTRUCT', &
                                 'FOUND IN MODE DEFINITION IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
            return
          endif
        endif
        em_mode%master_scale = ix
        do_evaluate = .false.

      case default
        call parser_error ('UNKNOWN MODE COMPONENT: ' // word, &
                             'FOUND IN MODE DEFINITION IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
        return
      end select

      if (do_evaluate) call evaluate_value (trim(ele%name), r_ptr, lat, delim, delim_found, err_flag, ',}')

      ! Possible "}" is end of mode
      if (delim == '}') exit

      call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

    enddo

    ! Check if map data has already been read in for another element.
    ! If so, save space by pointing to the existing map.

    branch_loop: do ib = 0, ubound(lat%branch, 1)
      branch => lat%branch(ib)
      do ie = 1, branch%n_ele_max
        bele => branch%ele(ie)
        if (.not. associated(bele%em_field)) cycle    
        if (bele%ix_ele == ele%ix_ele .and. bele%ix_branch == ele%ix_branch) cycle
        do im = 1, size(bele%em_field%mode)
          if (.not. associated(bele%em_field%mode(im)%map)) cycle
          if (bele%em_field%mode(im)%map%file /= em_mode%map%file) cycle
          deallocate(em_mode%map)
          em_mode%map => bele%em_field%mode(im)%map
          em_mode%map%n_link = em_mode%map%n_link + 1        
          exit branch_loop
        end do
      enddo 
    enddo branch_loop

    ! Expect "," or "}"
    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

    if (word /= '' .or. (delim /= '}' .and. delim /= ',')) then
      call parser_error ('NO "," OR "}" FOUND AFTER MODE DEFINITION', &
                           'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
      return
    endif

    if (delim == '}') exit

  enddo

  call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
  err_flag = .false.
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

  call get_integer (ix, .false., err_flag, 'BAD WIGGLER "TERM(IX)" CONSTRUCT'); if (err_flag) return

  if (delim /= ')') then
    call parser_error ('CANNOT FIND CLOSING ")" for a "TERM(i)" FOR A WIGGLER"', 'FOR: ' // ele%name)
    return
  endif

  write (str_ix, '(a, i3, a)') 'TERM(', ix, ')'

  if (.not. associated(ele%wig)) then
    allocate(ele%wig)
    allocate(ele%wig%term(ix))
  elseif (ix > size(ele%wig%term)) then
    allocate (wig_term(size(ele%wig%term)))
    wig_term = ele%wig%term
    deallocate (ele%wig%term)
    allocate (ele%wig%term(ix))
    ele%wig%term(1:size(wig_term)) = wig_term
    deallocate (wig_term)
  endif

! 1) chop "=", 2) chop to "{", 3) chop to "}", 4) chop to "," if it exists

  call get_next_word (word, ix_word1, ':,={}', delim1, delim_found, .true.) 
  call get_next_word (word, ix_word2, ':,={}', delim2, delim_found, .true., call_check = .true.)  

  if (delim1 /= '=' .or. delim2 /= '{' .or. ix_word1 /= 0 .or. ix_word2 /= 0) then
    call parser_error ('CONFUSED SYNTAX FOR TERM IN WIGGLER: ' // ele%name, str_ix)
    return
  endif

  err_str = trim(ele%name) // ' ' // str_ix

  call evaluate_value (err_str, ele%wig%term(ix)%coef, lat, delim, delim_found, err_flag, ',')
  if (err_flag) return
 
  call evaluate_value (err_str, ele%wig%term(ix)%kx, lat, delim, delim_found, err_flag, ',')
  if (err_flag) return

  call evaluate_value (err_str, ele%wig%term(ix)%ky, lat, delim, delim_found, err_flag, ',')
  if (err_flag) return

  call evaluate_value (err_str, ele%wig%term(ix)%kz, lat, delim, delim_found, err_flag, ',')
  if (err_flag) return

  call evaluate_value (err_str, ele%wig%term(ix)%phi_z, lat, delim, delim_found, err_flag, '}')
  if (err_flag) return

  kx = ele%wig%term(ix)%kx
  ky = ele%wig%term(ix)%ky
  kz = ele%wig%term(ix)%kz
  tol = 1e-5 * (kx**2 + ky**2 + kz**2)

  if (abs(ky**2 - kx**2 - kz**2) < tol) then
    ele%wig%term(ix)%type = hyper_y$
  elseif (abs(ky**2 + kx**2 - kz**2) < tol) then
    ele%wig%term(ix)%type = hyper_xy$
  elseif (abs(ky**2 - kx**2 + kz**2) < tol) then
    ele%wig%term(ix)%type = hyper_x$
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
  return

endif

! Check that next delim is a "=". 
! If not, it might be a flag attribute or an attribute that has a default value.

if (delim /= '=')  then
  err_flag = .false.

  if (ele%key == multipole$ .and. ix_attrib >= t0$) then
    if (.not. associated(ele%a_pole)) call multipole_init (ele)
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

  case ('REF_BEGINNING')
    if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
    pele%ref_pt = anchor_beginning$

  case ('REF_CENTER')
    if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
    pele%ref_pt = anchor_center$

  case ('REF_END')
    if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
    pele%ref_pt = anchor_end$

  case ('ELE_BEGINNING')
    if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
    pele%ele_pt = anchor_beginning$

  case ('ELE_CENTER')
    if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
    pele%ele_pt = anchor_center$

  case ('ELE_END')
    if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
    pele%ele_pt = anchor_end$

  case default
    call parser_error ('EXPECTING "=" AFTER ATTRIBUTE: ' // word,  'FOR ELEMENT: ' // ele%name)
    err_flag = .true.
  end select

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
  pele%s = value

case('TYPE', 'ALIAS', 'DESCRIP', 'SR_WAKE_FILE', 'LR_WAKE_FILE', 'LATTICE', 'TO', &
     'TO_LINE', 'TO_ELEMENT', 'CRYSTAL_TYPE', 'MATERIAL_TYPE', 'ORIGIN_ELE')
  call bmad_parser_type_get (ele, attrib_word, delim, delim_found, pele = pele)

case ('PTC_MAX_FRINGE_ORDER')
  call get_integer (bmad_com%ptc_max_fringe_order, .false., err_flag)
  bp_com%extra%ptc_max_fringe_order_set = .true.
  bp_com%extra%ptc_max_fringe_order = bmad_com%ptc_max_fringe_order

case ('TAYLOR_ORDER')
  call get_integer (ix, .false., err_flag)
  if (ix <= 0) then
    call parser_error ('TAYLOR_ORDER IS LESS THAN 1')
    return
  endif
  ptc_com%taylor_order_saved = ix

case ('SYMPLECTIFY') 
  if (how == def$ .and. (delim == ',' .or. .not. delim_found)) then
    ele%symplectify = .true.
  else
    call get_logical (attrib_word, ele%symplectify, err_flag)
  endif
  
case ('IS_ON')
  call get_logical (attrib_word, ele%is_on, err_flag)

case ('USE_HARD_EDGE_DRIFTS')
  call get_logical (attrib_word, bmad_com%use_hard_edge_drifts, err_flag)
  bp_com%extra%use_hard_edge_drifts_set = .true.
  bp_com%extra%use_hard_edge_drifts = bmad_com%use_hard_edge_drifts

case ('APERTURE_LIMIT_ON') 
  call get_logical (attrib_word, lat%param%aperture_limit_on, err_flag)

case ('ABSOLUTE_TIME_TRACKING')
  call get_logical (attrib_word, lat%absolute_time_tracking, err_flag)

case ('CREATE_JUMBO_SLAVE')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  call get_logical (attrib_word, pele%create_jumbo_slave, err_flag)

case ('USE_PTC_LAYOUT')
  call get_logical (attrib_word, lat%use_ptc_layout, err_flag)

case ('AUTO_SCALE_FIELD_PHASE')
  call get_logical (attrib_word, lat%auto_scale_field_phase, err_flag)

case ('AUTO_SCALE_FIELD_AMP')
  call get_logical (attrib_word, lat%auto_scale_field_amp, err_flag)

case ('CSR_CALC_ON')
  call get_logical (attrib_word, ele%csr_calc_on, err_flag)

case ('PTC_EXACT_MODEL')
  call get_logical (attrib_word, logic, err_flag)
  if (.not. err_flag) call set_ptc (exact_modeling = logic)

case ('PTC_EXACT_MISALIGN')
  call get_logical (attrib_word, logic, err_flag)
  if (err_flag) return
  call set_ptc (exact_misalign = logic)

case ('TAYLOR_MAP_INCLUDES_OFFSETS')
  call get_logical (attrib_word, ele%taylor_map_includes_offsets, err_flag)

case ('OFFSET_MOVES_APERTURE')
  call get_logical (attrib_word, ele%offset_moves_aperture, err_flag)

case ('FIELD_MASTER', 'HARMON_MASTER')
  call get_logical (attrib_word, ele%field_master, err_flag)

case ('SCALE_MULTIPOLES')
  call get_logical (attrib_word, ele%scale_multipoles, err_flag)

case ('FIELD_CALC')
  call get_switch (attrib_word, field_calc_name(1:), ele%field_calc, err_flag)

case ('APERTURE_AT')
  call get_switch (attrib_word, aperture_at_name(1:), ele%aperture_at, err_flag)

case ('REF_ORIGIN')
  call get_switch (attrib_word, anchor_pt_name(1:), pele%ref_pt, err_flag)

case ('ELE_ORIGIN')
  call get_switch (attrib_word, anchor_pt_name(1:), pele%ele_pt, err_flag)

case ('APERTURE_TYPE')
  call get_switch (attrib_word, aperture_type_name(1:), ele%aperture_type, err_flag)

case ('COUPLER_AT')
  call get_switch (attrib_word, end_at_name(1:), ix, err_flag)
  ele%value(coupler_at$) = ix

case ('FRINGE_TYPE')
  call get_switch (attrib_word, fringe_type_name(1:), ix, err_flag)
  ele%value(fringe_type$) = ix

case ('FRINGE_AT')
  call get_switch (attrib_word, end_at_name(1:), ix, err_flag)
  ele%value(fringe_at$) = ix

case ('ORIGIN_ELE_REF_PT')
  call get_switch (attrib_word, ref_pt_name(1:), ix, err_flag)
  ele%value(origin_ele_ref_pt$) = ix

case ('TRACKING_METHOD')
  call get_switch (attrib_word, tracking_method_name(1:), switch, err_flag)
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

case ('SPIN_TRACKING_METHOD')
  call get_switch (attrib_word, spin_tracking_method_name(1:), switch, err_flag)
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

case ('MAT6_CALC_METHOD')
  call get_switch (attrib_word, mat6_calc_method_name(1:), switch, err_flag)
  if (err_flag) return
  if (.not. valid_mat6_calc_method (ele, not_set$, switch)) then
    if (wild_key0) then
      err_flag = .false.
    else
      call parser_error ('NOT A VALID MAT6_CALC_METHOD: ' // word, &
                         'FOR: ' // trim(ele%name), 'WHICH IS A: ' // key_name(ele%key))
    endif
    return
  endif
  ele%mat6_calc_method = switch

case ('REF_COORDINATES')
  call get_switch (attrib_word, end_at_name(1:2), ix, err_flag)
  ele%value(ref_coordinates$) = ix

case ('REF_ORBIT_FOLLOWS')
  call get_switch (attrib_word, ref_orbit_follows_name(1:), ix, err_flag)
  ele%value(ref_orbit_follows$) = ix

case ('MODE')
  call get_switch (attrib_word, mode_name(1:), ix, err_flag)
  ele%value(geometry$) = ix

case ('PTC_INTEGRATION_TYPE')
  call get_switch (attrib_word, ptc_integration_type_name(1:), ele%ptc_integration_type, err_flag)

case ('PARTICLE')
  call get_switch (attrib_word, particle_name(:), ix, err_flag)
  ele%value(particle$) = ix + lbound(particle_name, 1) - 1 

case ('PTC_FIELD_GEOMETRY')
  call get_switch (attrib_word, ptc_field_geometry_name(1:), ix, err_flag)
  ele%value(ptc_field_geometry$) = ix

  if (ele%key == sbend$ .and. ix == true_rbend$) then
    call parser_error ('TRUE_RBEND IS NOT A VALID PTC_FIELD_GEOMETRY VALUE FOR AN SBEND')
    return
  endif

case ('GEOMETRY')
  call get_switch (attrib_word, geometry_name(1:), ix, err_flag)
  ele%value(geometry$) = ix

case ('PHOTON_TYPE')
  call get_switch (attrib_word, photon_type_name(1:), ix, err_flag)
  ele%value(photon_type$) = ix

case ('LATTICE_TYPE')   ! Old style
  call get_switch (attrib_word, lattice_type_name(1:), ix, err_flag)
  ele%value(geometry$) = ix

case default   ! normal attribute

  if (attribute_type(attrib_word) == is_logical$) then
    call get_logical_real (attrib_word, ele%value(ix_attrib), err_flag)

  else
    call evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag)
    if (err_flag) return

    ! multipole attribute?
    if (ix_attrib >= a0$ .and. ix_attrib <= b21$ .and. attrib_word(1:4) /= 'CURV') then  
        if (.not. associated(ele%a_pole)) call multipole_init (ele)
        if (ix_attrib >= b0$) then
          ele%b_pole(ix_attrib-b0$) = value
        else
          ele%a_pole(ix_attrib-a0$) = value
        endif
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
      call pointer_to_attribute (ele, attrib_word, .false., r_ptr, err_flag, .false.)
      if (err_flag) then
        call parser_error ('BAD ATTRIBUTE: ' // attrib_word, 'FOR ELEMENT: ' // ele%name)
        return
      endif
      r_ptr = value
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

      select case (attrib_word)

      case ('NUM_STEPS')
        ele%value(ds_step$) = abs(ele%value(l$) * nint(ele%value(num_steps$)))

      case ('E_TOT')
        if (ele%key == def_parameter$) then
          lat%ele(0)%value(e_tot$) = value
          lat%ele(0)%value(p0c$) = -1
        else
          ele%value(p0c$) = -1
        endif

      case ('ENERGY')    ! Only in def_beam
        lat%ele(0)%value(e_tot$) = 1d9 * value
        lat%ele(0)%value(p0c$) = -1

      case ('P0C')
        if (ele%key == def_parameter$) then
          lat%ele(0)%value(p0c$) = value
          lat%ele(0)%value(e_tot$) = -1
        else
          ele%value(e_tot$) = -1
        endif

      case ('PC')    ! Only in def_beam
        lat%ele(0)%value(p0c$) = 1d9 * value
        ele%value(e_tot$) = -1

      case ('LR_FREQ_SPREAD')
        call randomize_lr_wake_frequencies (ele, set_done)
        if (set_done) call bp_set_ran_status

      end select

    endif

  endif

end select

err_flag = .false.

!--------------------------------------------------------
contains

function expect_this (expecting, check_delim, call_check, err_str) result (is_ok)

implicit none

character(*) expecting, err_str
logical is_ok, check_delim, call_check
integer ix

!

is_ok = .false.

ix = 1
if (check_delim) then
  if (delim /= expecting(1:1)) then
    call parser_error ('NO "' // expecting // '" FOUND ' // err_str, 'FOR ELEMENT: ' // ele%name)
    return
  endif
  ix = 2
endif

do
  if (ix > len(expecting)) exit
  call get_next_word (word, ix_word, expecting(ix:ix), delim, delim_found, call_check = call_check)
  if (delim /= expecting(ix:ix) .or. word /= '') then
    call parser_error ('NO "' // expecting // '" FOUND ' // err_str, 'FOR ELEMENT: ' // ele%name)
    return
  endif
  ix = ix + 1
enddo

is_ok = .true.

end function expect_this

!--------------------------------------------------------
! contains

! delim_list -- character(*): List of valid tokens. If list contains a space character
!                 then no token (indicating the end of the command) is a valid possibility.

function expect_either (delim_list, check_delim) result (is_ok)

character(*) delim_list
logical check_delim, is_ok, must_have_delim

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
  call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
  if (word /= '' .or. (must_have_delim .and. .not. delim_found) .or. &
      (delim /= '' .and. index(delim_list, delim) == 0)) then
    call parser_error ('BAD DELIMITOR', 'FOR ELEMENT: ' // ele%name)
    return
  endif
endif

is_ok = .true.

end function expect_either

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
  is_free = attribute_free (ele, attrib_name, lat, bp_com%print_err)
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

subroutine get_integer (i, look_for_equal_sign, err, str1, str2)

integer i
logical look_for_equal_sign, err
character(*), optional :: str1, str2

!

if (look_for_equal_sign) then

endif

!

call get_next_word (word, ix_word, ':,=(){}', delim, delim_found, .true.)
if (.not. is_integer(word) ) then
  if (present(str1)) then
    call parser_error (str1, str2)
  else
    call parser_error ('INTEGER EXPECTED, I DO NOT UNDERSTAND: ' // word)
  endif
  err = .true.

else
  read (word, *) i 
  err = .false.
endif

end subroutine get_integer

!--------------------------------------------------------
! contains

subroutine get_switch (name, name_list, switch, err)

character(*) name, name_list(:)
integer this_switch, switch
logical err

!

call get_next_word (word, ix_word, ':,=(){}', delim, delim_found, .true.)
call match_word (word, name_list, this_switch, can_abbreviate = .false.)
if (this_switch < 1) then
  call parser_error ('BAD "' // trim(name) // '" SWITCH FOR: ' // ele%name, 'I DO NOT UNDERSTAND: ' // word)
  err = .true.
else
  err = .false.
  switch = this_switch
endif

end subroutine get_switch

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

type (ele_struct) ele
real(rp) coef
integer i_out, expn(6)

!

if (.not. associated(ele%taylor(1)%term)) then
  allocate (ele%taylor(1)%term(0), ele%taylor(2)%term(0))
  allocate (ele%taylor(3)%term(0), ele%taylor(4)%term(0))
  allocate (ele%taylor(5)%term(0), ele%taylor(6)%term(0))
endif

if (i_out < 1 .or. i_out > 6) then
  call parser_error ('"OUT" VALUE IN TAYLOR TERM NOT IN RANGE (1 - 6)', &
                'FOR TAYLOR ELEMENT: ' // ele%name)
  return
endif

call add_taylor_term (ele%taylor(i_out), coef, expn, .true.)

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

character(100) line
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
    if (n == 0 .or. n > 60) exit

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
  call fullfilename('./', file_name)
  bp_com%dirs(2) = file_name
  file(:)%dir = file_name
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
  if (is_relative) then
    call append_subdirectory (trim(file(i_level-1)%dir), file(i_level)%dir, file(i_level)%dir, err_flag)
    if (err_flag) call parser_error ('BAD DIRECTORY SYNTAX FOR: ' // file_name, stop_here = .true.)
  endif
  bp_com%dirs(1) = file(i_level-1)%dir
  call find_file (file_name2, found_it, file_name, bp_com%dirs)
  file(i_level)%logical_name = file_name2
  file(i_level)%full_name = file_name
  file(i_level)%f_unit = lunget()

  open (file(i_level)%f_unit, file = file_name, status = 'OLD', action = 'READ', iostat = ios)
  if (ios /= 0 .or. .not. found_it) then
    bp_com%current_file => file(i_level-1)  ! For warning
    if (file_name2 == file_name)  then
      call parser_error ('UNABLE TO OPEN FILE: ' // file_name, stop_here = .true.)
    else
      call parser_error ('UNABLE TO OPEN FILE: ' // file_name, &
                    'THIS FROM THE LOGICAL FILE NAME: ' // file_name2, &
                    stop_here = .true.)
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

  ! The same file may be validly called multiple times if it is an inline file.
  ! EG: A wall file called inline.

  if (how == 'push') then
    do i = 1, n_file - 1
      if (bp_com%lat_file_names(i) /= bp_com%lat_file_names(n_file)) cycle
      call parser_error ('Same lattice file called multiple times: ' // trim(bp_com%lat_file_names(n_file)), &
                         warn_only = .true.)
      exit
    enddo
  endif

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
  elseif (action == 'new_command') then
    bp_com%input_line1 = ' '
    bp_com%input_line2 = line
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
! Function of evaluate the index of an array. Typically the text peing parsed looks like:
!      "5) = ..."  or
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
character(len(word)+8) :: wd
logical this_logic
integer, intent(out) :: iostat
integer i

!

iostat = -1
this_logic = .false.  ! To avoid uninit compiler warnings.

call str_upcase(wd, word)
do i = 1, len(word)
  if (wd(i:i) /= ' ') then
    if (wd(i:i+6) == '.TRUE. ' .or. wd(i:i+4) == 'TRUE ' .or. &
                                          wd(i:i+1) == 'T ') then
      this_logic = .true.
      iostat = 0
    elseif (wd(i:i+7) == '.FALSE. ' .or. wd(i:i+5) == 'FALSE ' .or. &
                                          wd(i:i+1) == 'F ') then
      this_logic = .false.
      iostat = 0
    endif
    return
  endif
enddo

end function evaluate_logical

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine evaluate_value (err_str, value, lat, delim, delim_found, err_flag, end_delims)
!
! This routine creates an "evaluation stack" structure which can be used 
! to evaluate an arithmethic expression.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine evaluate_value (err_str, value, lat, delim, delim_found, err_flag, end_delims)

use random_mod

implicit none

type (lat_struct)  lat
type (eval_stack_struct) stk(200)

integer i_lev, i_op, i

integer op(200), ix_word, i_delim, i2, ix_word2

real(rp) value

character(*) err_str
character(*), optional :: end_delims
character(1) delim
character(80) word, word2

logical delim_found, split, ran_function_pending, first_get_next_word_call
logical err_flag

! The general idea is to rewrite the expression on a stack in reverse polish.
! Reverse polish means that the operand goes last so that 2 * 3 is writen 
! on the stack as: [2, 3, *]

! The stack is called: stk
! Since operations move towards the end of the stack we need a separate
! stack called op which keeps track of what operations have not yet
! been put on stk.

! init

err_flag = .true.
i_lev = 0
i_op = 0
ran_function_pending = .false.
first_get_next_word_call = .true.

! parsing loop to build up the stack.

parsing_loop: do

  ! get a word

  if (first_get_next_word_call) then
    call get_next_word (word, ix_word, '+-*/()^,:} ', delim, delim_found, call_check = .true.)
    first_get_next_word_call = .false.
  else
    call get_next_word (word, ix_word, '+-*/()^,:} ', delim, delim_found)
  endif

  if (delim == '*' .and. word(1:1) == '*') then
    call parser_error ('EXPONENTIATION SYMBOL IS "^" AS OPPOSED TO "**"!',  &
                  'for: ' // err_str)
    return
  endif

  if (ran_function_pending .and. (ix_word /= 0 .or. delim /= ')')) then
    call parser_error ('RAN AND RAN_GAUSS DO NOT TAKE AN ARGUMENT', 'FOR: ' // err_str)
    return
  endif

  !--------------------------
  ! Preliminary: If we have split up something that should have not been split
  ! then put it back together again...

  ! just make sure we are not chopping a number in two, e.g. "3.5e-7" should not
  ! get split at the "-" even though "-" is a delimiter

  split = .true.         ! assume initially that we have a split number
  if (ix_word == 0) then
    split = .false.
  elseif (word(ix_word:ix_word) /= 'E' .and. word(ix_word:ix_word) /= 'D') then
    split = .false.
  endif
  if (delim(1:1) /= '-' .and. delim(1:1) /= '+') split = .false.
  do i = 1, ix_word-1
    if (index('.0123456789', word(i:i)) == 0) split = .false.
  enddo

  ! If still SPLIT = .TRUE. then we need to unsplit

  if (split) then
    word = word(:ix_word) // delim
    call get_next_word (word2, ix_word2, '+-*/()^,:}', delim, delim_found)
    word = word(:ix_word+1) // word2
    ix_word = ix_word + ix_word2
  endif

  ! Something like "lcav[lr(2).freq]" will get split on the "("

  if (delim == '(' .and. index(word, '[LR') /= 0) then
    call get_next_word (word2, ix_word2, '+-*/(^,:}', delim, delim_found)
    word = word(:ix_word) // '(' // word2
    ix_word = ix_word + ix_word2 + 1
  endif

  !---------------------------
  ! Now see what we got...

  ! For a "(" delim we must have a function

  if (delim == '(') then

    ran_function_pending = .false.
    if (ix_word /= 0) then
      select case (word)
      case ('SIN') 
        call pushit (op, i_op, sin$)
      case ('COS') 
        call pushit (op, i_op, cos$)
      case ('TAN') 
        call pushit (op, i_op, tan$)
      case ('ASIN') 
        call pushit (op, i_op, asin$)
      case ('ACOS') 
        call pushit (op, i_op, acos$)
      case ('ATAN') 
        call pushit (op, i_op, atan$)
      case ('ATAN2') 
        call pushit (op, i_op, atan2$)
      case ('ABS') 
        call pushit (op, i_op, abs$)
      case ('SQRT') 
        call pushit (op, i_op, sqrt$)
      case ('LOG') 
        call pushit (op, i_op, log$)
      case ('EXP') 
        call pushit (op, i_op, exp$)
      case ('FACTORIAL') 
        call pushit (op, i_op, factorial$)
      case ('RAN') 
        call pushit (op, i_op, ran$)
        ran_function_pending = .true.
        call bp_set_ran_status
      case ('RAN_GAUSS') 
        call pushit (op, i_op, ran_gauss$)
        ran_function_pending = .true.
        call bp_set_ran_status
      case default
        call parser_error ('UNEXPECTED CHARACTERS ON RHS BEFORE "(": ' // word,  &
                                                'FOR: ' // err_str)
        return
      end select
    endif

    call pushit (op, i_op, l_parens$)
    cycle parsing_loop

  ! for a unary "-"

  elseif (delim == '-' .and. ix_word == 0) then
    call pushit (op, i_op, unary_minus$)
    cycle parsing_loop

  ! for a unary "+"

  elseif (delim == '+' .and. ix_word == 0) then
    call pushit (op, i_op, unary_plus$)
    cycle parsing_loop

  ! for a ")" delim

  elseif (delim == ')') then
    if (ix_word == 0) then
      if (.not. ran_function_pending) then
        call parser_error  ('CONSTANT OR VARIABLE MISSING BEFORE ")"', 'FOR: ' // err_str)
        return
      endif
      ran_function_pending = .false.
    else
      call word_to_value (word, lat, value)
      call pushit (stk%type, i_lev, numeric$)
      stk(i_lev)%value = value
    endif

    do
      do i = i_op, 1, -1     ! release pending ops
        if (op(i) == l_parens$) exit          ! break do loop
        call pushit (stk%type, i_lev, op(i))
      enddo

      if (i == 0) then
        if (index(end_delims, ')') /= 0) then
          i_op = 0
          exit   ! End of expression
        endif
        call parser_error ('UNMATCHED ")" ON RHS', 'FOR: ' // err_str)
        return
      endif

      i_op = i - 1

      call get_next_word (word, ix_word, '+-*/()^,:}', delim, delim_found)
      if (ix_word /= 0) then
        call parser_error ('UNEXPECTED CHARACTERS ON RHS AFTER ")"',  &
                                                  'FOR: ' // err_str)
        return
      endif

      if (delim /= ')') exit  ! if no more ')' then no need to release more
    enddo


    if (delim == '(') then
      call parser_error ('")(" CONSTRUCT DOES NOT MAKE SENSE FOR: ' // err_str)
      return
    endif

  ! For binary "+-/*^" delims

  else
    if (ix_word == 0) then
      call parser_error ('CONSTANT OR VARIABLE MISSING IN EVALUATING: ' // err_str)
      return
    endif
    call word_to_value (word, lat, value)
    call pushit (stk%type, i_lev, numeric$)
    stk(i_lev)%value = value
  endif

  ! If we are here then we have an operation that is waiting to be identified

  if (.not. delim_found) delim = ':'

  select case (delim)
  case ('+')
    i_delim = plus$
  case ('-')
    i_delim = minus$
  case ('*')
    i_delim = times$
  case ('/')
    i_delim = divide$
  case ('^')
    i_delim = power$
  case (',', '}', ':', ')')   ! End of expression delims
    i_delim = no_delim$
  case default
    call parser_error ('MALFORMED EXPRESSION')
    bp_com%parse_line = ' '
    return
  end select

  ! now see if there are operations on the OP stack that need to be transferred
  ! to the STK stack

  do i = i_op, 1, -1
    if (eval_level(op(i)) >= eval_level(i_delim)) then
      if (op(i) == l_parens$) then
        if (i > 1 .and. op(max(1,i-1)) == atan2$ .and. delim == ',') cycle parsing_loop
        call parser_error ('UNMATCHED "(" IN EVALUATING: ' // err_str)
        return
      endif
      call pushit (stk%type, i_lev, op(i))
    else
      exit
    endif
  enddo

  ! put the pending operation on the OP stack

  i_op = i
  if (i_delim == no_delim$) then
    exit parsing_loop
  else
    call pushit (op, i_op, i_delim)
  endif

enddo parsing_loop

!------------------------------------------------------------------
! now go through the stack and perform the operations

if (i_op /= 0) then
  call parser_error ('UNMATCHED "(" IN EVALUATING: ' // err_str)
  return
endif

if (i_lev == 0) call parser_error ('NO VALUE FOUND FOR: ' // err_str)

i2 = 0
do i = 1, i_lev
  if (stk(i)%type == numeric$) then
    i2 = i2 + 1
    stk(i2)%value = stk(i)%value

  elseif (stk(i)%type == unary_minus$) then
    stk(i2)%value = -stk(i2)%value

  elseif (stk(i)%type == unary_plus$) then
    stk(i2)%value = stk(i2)%value

  elseif (stk(i)%type == plus$) then
    stk(i2-1)%value = stk(i2-1)%value + stk(i2)%value
    i2 = i2 - 1

  elseif (stk(i)%type == minus$) then
    stk(i2-1)%value = stk(i2-1)%value - stk(i2)%value
    i2 = i2 - 1

  elseif (stk(i)%type == times$) then
    stk(i2-1)%value = stk(i2-1)%value * stk(i2)%value
    i2 = i2 - 1

  elseif (stk(i)%type == divide$) then
    if (stk(i2)%value == 0) then
      call parser_error ('DIVIDE BY 0 ON RHS', 'FOR: ' // err_str)
      return
    endif
    stk(i2-1)%value= stk(i2-1)%value / stk(i2)%value
    i2 = i2 - 1

  elseif (stk(i)%type == power$) then
    stk(i2-1)%value = stk(i2-1)%value**stk(i2)%value
    i2 = i2 - 1

  elseif (stk(i)%type == sin$) then
    stk(i2)%value = sin(stk(i2)%value)

  elseif (stk(i)%type == cos$) then
    stk(i2)%value = cos(stk(i2)%value)

  elseif (stk(i)%type == tan$) then
    stk(i2)%value = tan(stk(i2)%value)

  elseif (stk(i)%type == asin$) then
    if (stk(i2)%value < -1 .or. stk(i2)%value > 1) then
      call parser_error ('ASIN ARGUMENT HAS MAGNITUDE GREATER THAN 1', 'FOR: ' // err_str)
      return
    endif
    stk(i2)%value = asin(stk(i2)%value)

  elseif (stk(i)%type == acos$) then
    if (stk(i2)%value < -1 .or. stk(i2)%value > 1) then
      call parser_error ('ACOS ARGUMENT HAS MAGNITUDE GREATER THAN 1', 'FOR: ' // err_str)
      return
    endif
    stk(i2)%value = acos(stk(i2)%value)

  elseif (stk(i)%type == factorial$) then
    stk(i2)%value = factorial(nint(stk(i2)%value))
    if (stk(i2)%value < 0) then
      call parser_error ('FACTORIAL PROBLEM FOR: ' // err_str)
      return
    endif

  elseif (stk(i)%type == atan$) then
    stk(i2)%value = atan(stk(i2)%value)

  elseif (stk(i)%type == atan2$) then
    stk(i2-1)%value = atan2(stk(i2-1)%value, stk(i2)%value)
    i2 = i2 - 1

  elseif (stk(i)%type == abs$) then
    stk(i2)%value = abs(stk(i2)%value)

  elseif (stk(i)%type == sqrt$) then
    if (stk(i2)%value < 0) then
      call parser_error ('SQRT ARGUMENT IS NEGATIVE ', 'FOR: ' // err_str)
      return
    endif
    stk(i2)%value = sqrt(stk(i2)%value)

  elseif (stk(i)%type == log$) then
    if (stk(i2)%value < 0) then
      call parser_error ('LOG ARGUMENT IS NEGATIVE ', 'FOR: ' // err_str)
      return
    endif
    stk(i2)%value = log(stk(i2)%value)

  elseif (stk(i)%type == exp$) then
    stk(i2)%value = exp(stk(i2)%value)

  elseif (stk(i)%type == ran$) then
    i2 = i2 + 1
    call ran_uniform(stk(i2)%value)

  elseif (stk(i)%type == ran_gauss$) then
    i2 = i2 + 1
    call ran_gauss(stk(i2)%value)

  else
    call parser_error ('INTERNAL ERROR #02: GET HELP')
    if (global_com%exit_on_error) call err_exit
  endif
enddo

if (i2 /= 1) then
  call parser_error ('INTERNAL ERROR #03: GET HELP')
  if (global_com%exit_on_error) call err_exit
endif

value = stk(1)%value
err_flag = .false.

! Check that final delim matches.

if (present(end_delims)) then
  if (.not. delim_found .or. index(end_delims, delim) == 0) then
    call parser_error ('BAD DELIMITOR AFTER VALUE FOR: ' // err_str)
    err_flag = .true.
  endif
endif


end subroutine evaluate_value

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine pushit (stack, i_lev, value)

implicit none

integer stack(:), i_lev, value

!

i_lev = i_lev + 1

if (i_lev > size(stack)) then
  call parser_error ('STACK OVERFLOW.', 'EXPERT HELP IS NEEDED!')
  if (global_com%exit_on_error) call err_exit
endif

stack(i_lev) = value

end subroutine pushit
                     
!-------------------------------------------------------------------------
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
type (real_pointer_struct), allocatable :: ptr(:)
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_attribute_struct) attrib_info

integer i, ix1, ix2, ix_word, ios, ix, n_loc, ix_attrib
integer ivar, ixm, ixm2
real(rp) value
real(rp), pointer :: v(:)
character(*) word
character(40) attrib_name, ele_name
logical err_flag

! see if this is numeric

if (index('-+.0123456789', word(1:1)) /= 0) then
  read (word, *, iostat = ios) value
  if (ios /= 0) call parser_error ('BAD VARIABLE: ' // word)
  return
endif

! If not numeric...

ix_word = len_trim(word)
if (.not. verify_valid_name (word, ix_word)) return

! If word does not have a "[...]" then it must be a variable

ix1 = index(word, '[')
if (ix1 == 0) then   
  call find_indexx (word, bp_com%var%name, bp_com%var%indexx, bp_com%ivar_tot, i)
  if (i == 0) then
    call parser_error ('VARIABLE USED BUT NOT YET DEFINED: ' // word)
    value = 0
    ! To prevent multiple error messages define this variable.
    bp_com%ivar_tot = bp_com%ivar_tot + 1
    if (bp_com%ivar_tot > size(bp_com%var%name)) call reallocate_bp_com_var()
    ivar = bp_com%ivar_tot
    bp_com%var(ivar)%name = word
    bp_com%var(ivar)%value = 0
    call find_indexx (word, bp_com%var%name, bp_com%var%indexx, ivar-1, ixm, ixm2)
    bp_com%var(ixm2+1:ivar)%indexx = bp_com%var(ixm2:ivar-1)%indexx
    bp_com%var(ixm2)%indexx = ivar

  else
    value = bp_com%var(i)%value
  endif
  return
endif

! Here if word does have a "[...]" then is a element attribute

ele_name = word(:ix1-1)    ! name of attribute

ix2 = index(word, ']')
attrib_name = word(ix1+1:ix2-1)

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
  else
    value = ptr(1)%r
  endif

  ! If this is bmad_parser, and not bmad_parser2, then dependent attributes have not been set and cannot
  ! be used.

  if (ix_attrib > 0 .and. bp_com%parser_name == 'bmad_parser') then
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
integer i, ivar, ixm, ixm2
logical delim_found, err_flag

!

call find_indexx (word, bp_com%var%name, bp_com%var%indexx, bp_com%ivar_tot, i)
if (i /= 0) then
  call parser_error ('VARIABLES ARE NOT ALLOWED TO BE REDEFINED: ' // word)
  call evaluate_value (word, bp_com%var(i)%value, lat, delim, delim_found, err_flag)
  return
endif

bp_com%ivar_tot = bp_com%ivar_tot + 1
if (bp_com%ivar_tot > size(bp_com%var%name)) call reallocate_bp_com_var()
ivar = bp_com%ivar_tot
bp_com%var(ivar)%name = word
call evaluate_value (bp_com%var(ivar)%name, bp_com%var(ivar)%value, &
                                     lat, delim, delim_found, err_flag)
if (delim_found .and. .not. err_flag) call parser_error  &
                  ('EXTRA CHARACTERS ON RHS: ' // bp_com%parse_line,  &
                   'FOR VARIABLE: ' // bp_com%var(ivar)%name)

! Reindex.

call find_indexx (word, bp_com%var%name, bp_com%var%indexx, ivar-1, ixm, ixm2)

bp_com%var(ixm2+1:ivar)%indexx = bp_com%var(ixm2:ivar-1)%indexx
bp_com%var(ixm2)%indexx = ivar

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

integer ix, ix_word

character(*) attrib_name
character(*), optional :: name, str_out
character(40)  word
character(1)   delim, str_end
character(200) type_name

logical delim_found

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
  ele%alias = type_name
case ('CRYSTAL_TYPE', 'MATERIAL_TYPE')
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
!     %wake%lr(:)       -- Long-range wake potential.
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
logical set_done
namelist / long_range_modes / lr

! Init

if (.not. associated(ele%wake)) allocate (ele%wake)
if (.not. allocated(ele%wake%sr_long%mode))  allocate (ele%wake%sr_long%mode(0))
if (.not. allocated(ele%wake%sr_trans%mode)) allocate (ele%wake%sr_trans%mode(0))
if (allocated(ele%wake%lr)) deallocate (ele%wake%lr)

! get data

call find_this_file (iu, lr_file_name, full_file_name)
if (iu < 0) return

ele%wake%lr_file = lr_file_name

lr%freq = -1
lr%angle = ''
lr%b_sin = 0
lr%b_cos = 0
lr%a_sin = 0
lr%a_cos = 0
lr%t_ref = 0

read (iu, nml = long_range_modes, iostat = ios)
close (iu)
if (ios > 0 .or. lr(1)%freq == -1) then
  call parser_error ('CANNOT READ LONG_RANGE_MODES NAMELIST FOR ELEMENT: ' // ele%name, & 
                'FROM FILE: '// full_file_name)
  return
endif

n_row = count(lr%freq /= -1)
allocate (ele%wake%lr(n_row))
j = 0
do i = 1, size(lr)
  if (lr(i)%freq == -1) cycle

  j = j + 1
  ele%wake%lr(j)%freq_in   = lr(i)%freq
  ele%wake%lr(j)%freq      = lr(i)%freq
  ele%wake%lr(j)%r_over_q  = lr(i)%r_over_q
  ele%wake%lr(j)%q         = lr(i)%q
  ele%wake%lr(j)%m         = lr(i)%m
  ele%wake%lr(j)%b_sin     = lr(i)%b_sin
  ele%wake%lr(j)%b_cos     = lr(i)%b_cos
  ele%wake%lr(j)%a_sin     = lr(i)%a_sin
  ele%wake%lr(j)%a_cos     = lr(i)%a_cos
  ele%wake%lr(j)%t_ref     = lr(i)%t_ref

  call downcase_string(lr(i)%angle)
  if (lr(i)%angle == '') then
    call parser_error ('LONG_RANGE_MODE ANGLE IS MISSING. MUST BE NUMBER OR "UNPOLARIZED"', & 
                  'FOR ELEMENT: ' // ele%name, &
                  'IN FILE: ' // full_file_name)
    cycle
  endif

  if (index('unpolarized', trim(lr(j)%angle)) == 1) then
    ele%wake%lr(j)%polarized = .false.
    ele%wake%lr(j)%angle     = 0
  else
    ele%wake%lr(j)%polarized = .true.
    read (lr(j)%angle, *, iostat = ios) ele%wake%lr(j)%angle
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
type (wake_sr_mode_struct) longitudinal(100), transverse(100)

real(rp) z_max
integer n, j, iu, ios, ix, i

character(*) sr_file_name
character(80) line
character(200) full_file_name

logical found_it

namelist / short_range_modes / z_max, longitudinal, transverse

! init

if (.not. associated(ele%wake))   allocate (ele%wake)
if (.not. allocated(ele%wake%lr)) allocate (ele%wake%lr(0))
if (allocated(ele%wake%sr_long%mode))  deallocate (ele%wake%sr_long%mode)
if (allocated(ele%wake%sr_trans%mode)) deallocate (ele%wake%sr_trans%mode)

! Open file

iu = 0
ele%wake%sr_file = sr_file_name
call find_this_file (iu, sr_file_name, full_file_name)
if (iu < 0) return

! Get data

longitudinal(:)%phi = real_garbage$
transverse(:)%phi = real_garbage$
z_max = real_garbage$

read (iu, nml = short_range_modes, iostat = ios)
close (iu)
if (ios > 0) then
  call parser_error ('CANNOT READ SHORT_RANGE_MODES NAMELIST FROM FILE: ' & 
                    // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

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

end subroutine read_sr_wake

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine find_this_file (iu, file, full_file_name)
!
! Subroutine to open a file.
! This subroutine is used by bmad_parser and bmad_parser2.
!
! Input:
!   file -- Character(*): Name of wake field file
!
! Output:
!   iu             -- Integer: Open file unit. 
!                       If negative then there was an error. 
!   full_file_name -- Character(*): Full name with directory spec.
!-

subroutine find_this_file (iu, file, full_file_name)

implicit none

character(*) file, full_file_name
integer iu, ios, n
logical found_it

! open file

iu = lunget()
bp_com%dirs(2) = bp_com%calling_file%dir
call find_file (file, found_it, full_file_name, bp_com%dirs)
open (iu, file = full_file_name, status = 'OLD', action = 'READ', iostat = ios)
if (ios /= 0) then
call parser_error ('CANNOT OPEN WAKE FILE: ' // file)
iu = -1
return
endif

! If we have not read in this file before then add this to the list of files
! that are used to create the lattice.

n = size(bp_com%lat_file_names)
if (n < bp_com%num_lat_files + 1) call re_allocate (bp_com%lat_file_names, n + 100)
n = bp_com%num_lat_files
inquire (file = full_file_name, name = bp_com%lat_file_names(n+1))
if (all(bp_com%lat_file_names(n+1) /= bp_com%lat_file_names(1:n))) &
                                              bp_com%num_lat_files = n + 1
end subroutine find_this_file

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_overlay_group_names (ele, lat, pele, delim, delim_found)
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-
      
subroutine get_overlay_group_names (ele, lat, pele, delim, delim_found)

implicit none

type (ele_struct)  ele
type (parser_ele_struct) pele
type (lat_struct)  lat

real(rp) coef(200)
real(rp) value

integer ix_word, ixs, j, k
                           
character(1) delim
character(40) word_in, word
character(40) name(200), attrib_name(200)

logical delim_found, err_flag, end_of_file
                    
!

call get_next_word (word_in, ix_word, '{,}', delim, delim_found, .true.)
if (delim /= '{' .or. ix_word /= 0) call parser_error  &
        ('BAD ' // control_name(ele%lord_status) // 'SPEC: ' // word_in,  &
        'FOR ELEMENT: ' // ele%name)

! loop over all names in "{...}" list

do 

  call get_next_word (word_in, ix_word, '{,}/:', delim, delim_found, .true.)

  ! If "{}" with no slaves... 
  if (delim == '}' .and. ix_word == 0 .and. ele%n_slave == 0) then
    call get_next_word (word, ix_word, ',=:', delim, delim_found, .true.)
    exit
  endif

  ele%n_slave = ele%n_slave + 1
  ixs = ele%n_slave
  word = word_in

  j = index(word, '[')
  if (j > 1) then
    k = index(word, ']')
    if (k <= j+1) then
      call parser_error ('BAD ATTRIBUTE SPEC: ' // word_in, 'FOR: ' // ele%name)
      word = word(:k-1) // word(j+1:)
    else
      attrib_name(ixs) = word(j+1:k-1)
      word = word(:j-1) // word(k+1:)
    endif
  else
    attrib_name(ixs) = blank_name$
  endif

  name(ixs) = word

  if (delim == '/' .or. delim == ':') then
    call evaluate_value (trim(ele%name), value, lat, delim, delim_found, err_flag)
    if (err_flag) then
      call parser_error ('BAD COEFFICIENT: ' // word_in,  &
                                        'FOR ELEMENT: ' // ele%name)
      call load_parse_line ('new_command', 1, end_of_file)         ! next line
      return
    endif
    coef(ixs) = value
  else
    coef(ixs) = 1.0
  endif

  if (delim == '}') then
    call get_next_word (word, ix_word, ',=:', delim, delim_found, .true.)
    exit
  elseif (delim /= ',') then
    call parser_error ('BAD ' // control_name(ele%lord_status) //  &
            'SPEC: ' // word_in, 'FOR: ' // ele%name)
    exit
  endif
                        
enddo

!

ixs = ele%n_slave

! if (ixs == 0) call parser_error ( &
!        'NO SLAVE ELEMENTS ASSOCIATED WITH GROUP/OVERLAY ELEMENT: ' // ele%name)

allocate (pele%coef(ixs), pele%name(ixs), pele%attrib_name(ixs))
pele%coef = coef(1:ixs)
pele%name = name(1:ixs)
pele%attrib_name = attrib_name(1:ixs)

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

subroutine parser_error (what1, what2, what3, seq, pele, stop_here, warn_only)

implicit none

type (seq_struct), optional :: seq
type (parser_ele_struct), optional :: pele

character(*) what1
character(*), optional :: what2, what3
character(160) lines(12)
character(16), parameter :: r_name = 'parser_error'
integer nl, err_level
logical, optional :: stop_here, warn_only

! bp_com%error_flag is a common logical used so program will stop at end of parsing

if (global_com%type_out .and. bp_com%print_err) then

  nl = 0

  if (logic_option(.false., warn_only)) then
    nl=nl+1; lines(nl) = 'WARNING IN ' // trim(bp_com%parser_name) // ': ' // trim(what1)
    err_level = s_warn$
  else
    nl=nl+1; lines(nl) = 'ERROR IN ' // trim(bp_com%parser_name) // ': ' // trim(what1)
    err_level = s_error$
  endif

  if (present(what2)) then
    nl=nl+1; lines(nl) = '     ' // trim(what2)
  endif

  if (present(what3)) then
    nl=nl+1; lines(nl) = '     ' // trim(what3)
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

  call out_io (err_level, r_name, lines(1:nl))

endif

! Warnings do not result in bp_com%error_flag being set

if (.not. logic_option(.false., warn_only)) then
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

nn = 15  ! number of standard (non-user defined) constants
bp_com%ivar_init = nn
bp_com%ivar_tot  = nn

allocate (bp_com%var(nn))

bp_com%var( 1) = bp_var_struct('PI', pi, 0)
bp_com%var( 2) = bp_var_struct('TWOPI', twopi, 0)
bp_com%var( 3) = bp_var_struct('DEGRAD', 180 / pi, 0)
bp_com%var( 4) = bp_var_struct('RADDEG', pi / 180, 0)
bp_com%var( 5) = bp_var_struct('E_LOG', 2.718281828459_rp, 0)
bp_com%var( 6) = bp_var_struct('E_MASS', e_mass, 0)
bp_com%var( 7) = bp_var_struct('C_LIGHT', c_light, 0)
bp_com%var( 8) = bp_var_struct('M_ELECTRON', m_electron, 0)
bp_com%var( 9) = bp_var_struct('M_PROTON', m_proton, 0)
bp_com%var(10) = bp_var_struct('R_P', r_p, 0)
bp_com%var(11) = bp_var_struct('E_CHARGE', e_charge, 0)
bp_com%var(12) = bp_var_struct('EMASS', e_mass, 0)
bp_com%var(13) = bp_var_struct('CLIGHT', c_light, 0)
bp_com%var(14) = bp_var_struct('R_E', r_e, 0)
bp_com%var(15) = bp_var_struct('DEGREES', pi / 180, 0) ! From degrees to radians.

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

lord%old_is_on = .false.  ! So add_all_superimpose will not try to use as ref ele.
lord%lord_status = multipass_lord$
lord%n_slave = n_multipass
lord%ix1_slave = 0
lord%ix2_slave = -1
call add_lattice_control_structs (lat, lord)

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
  lat%control(ixc)%ix_lord = ix_lord
  lat%control(ixc)%ix_slave = slave%ix_ele
  lat%control(ixc)%ix_branch = slave%ix_branch
  if (slave%n_lord /= 0) then
    call parser_error ('INTERNAL ERROR: CONFUSED MULTIPASS SETUP.', &
                  'PLEASE GET EXPERT HELP!')
    if (global_com%exit_on_error) call err_exit
  endif
  slave%n_lord = 1
  write (slave%name, '(2a, i1)') trim(slave%name), '\', i   ! '
  call add_lattice_control_structs (lat, slave)
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

subroutine add_all_superimpose (branch, super_ele_in, pele, in_lat)

use multipass_mod

implicit none

type (ele_struct) super_ele_in
type (ele_struct), save :: super_ele_saved, super_ele
type (ele_struct), pointer :: ref_ele, ele, slave, lord, super_ele_out
type (ele_pointer_struct), allocatable :: eles(:)
type (parser_ele_struct) pele
type (multipass_all_info_struct) m_info
type (lat_struct), optional :: in_lat
type (lat_ele_loc_struct), allocatable :: m_slaves(:)
type (branch_struct), target :: branch
type (lat_struct), pointer :: lat

integer ix, i, j, k, it, nic, nn, i_ele, ib
integer n_con, ix_branch, n_loc

character(40) name
character(40), allocatable :: multi_name(:)
character(80) line

logical have_inserted, found, err_flag, err

! init

if (.not. bp_com%do_superimpose) return

call settable_dep_var_bookkeeping (super_ele_in)

call init_ele(super_ele_saved)
call init_ele(super_ele)

super_ele_saved = super_ele_in      ! in case super_ele_in changes
super_ele = super_ele_saved        ! 
lat => branch%lat

! If no refrence point then superposition is simple

if (pele%ref_name == blank_name$) then
  call compute_super_lord_s (branch%ele(0), super_ele, pele)
  call check_for_multipass_superimpose_problem (branch%ele(0), super_ele, err_flag); if (err_flag) return
  call add_superimpose (lat, super_ele, 0, err_flag, save_null_drift = .true., &
                                 create_jumbo_slave = pele%create_jumbo_slave)
  if (err_flag) bp_com%error_flag = .true.
  return
endif

! Find all matches

call lat_ele_locator (pele%ref_name, lat, eles, n_loc, err)
if (err) then
  call parser_error ('MALFORMED SUPERIMPOSE REFERENCE ELEMENT NAME: ' // pele%ref_name, &
                     'FOR SUPERPOSITION OF: ' // super_ele_saved%name, pele = pele)
  return
endif

! If no match and the reference element has been defined but not used in the lattice
! then this is not an error

if (n_loc == 0) then
  do i = 1, in_lat%n_ele_max
    found = match_wild(in_lat%ele(i)%name, pele%ref_name)
    if (found) exit
  enddo
  if (.not. found) call parser_error ('NO MATCH FOR REFERENCE ELEMENT: ' //  pele%ref_name, &
                                     'FOR SUPERPOSITION OF: ' // super_ele_saved%name, pele = pele)
endif

! Tag reference elements using %old_is_on flag which is not 

branch%ele(:)%old_is_on = .false.  
do i = 1, n_loc
  eles(i)%ele%old_is_on = .true. ! Tag reference element.
enddo

!

do 

  have_inserted = .false.

  ele_loop: do i_ele = 0, branch%n_ele_max

    ref_ele => branch%ele(i_ele)
     
    if (ref_ele%key == group$ .or. ref_ele%slave_status == super_slave$) cycle
    if (ref_ele%key == girder$) cycle
    if (.not. ref_ele%old_is_on) cycle

    ref_ele%old_is_on = .false.  ! So only use this reference once

    ! If superimposing on a multipass_lord then the superposition
    ! must be done at all multipass locations.

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
        call compute_super_lord_s (ele, super_ele, pele)
        super_ele%iyy = ele%iyy   ! Multipass info
        call check_for_multipass_superimpose_problem (ele, super_ele, err_flag); if (err_flag) return
        ! Don't need to save drifts since a multipass_lord drift already exists.
        call add_superimpose (lat, super_ele, ix_branch, err_flag, super_ele_out, &
                      save_null_drift = .false., create_jumbo_slave = pele%create_jumbo_slave)
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

      do i = 1, size(m_info%top) ! Loop over multipass lords
        do j = 1, size(m_info%top(i)%slave, 2)   ! loop over super_slaves
          slave => m_info%top(i)%slave(1, j)%ele
          if (slave%key /= drift$) cycle
          if (slave%slave_status == multipass_slave$) cycle
          do k = 1, size(m_info%top(i)%slave(:, j))   ! Loop over all passes
            ele => m_info%top(i)%slave(k, j)%ele
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
        m_slaves(i)%ix_ele = m_slaves(i)%ix_ele - &
                                count(lat%branch(ele%ix_branch)%ele(1:ele%ix_ele)%key == -1)
      enddo

      call remove_eles_from_lat (lat, .false.)

      ! Add a multipass_lord to control the created super_lords.

      call add_this_multipass (lat, m_slaves, super_ele_saved) 

      call deallocate_multipass_all_info_struct (m_info)
      deallocate (m_slaves, multi_name)
      deallocate (ele_loc_com%branch)

    !-----------------------
    ! Else not superimposing on a multipass_lord ...

    else
      call compute_super_lord_s (branch%ele(i_ele), super_ele, pele)
      super_ele%iyy = branch%ele(i_ele)%iyy   ! Multipass info
      call check_for_multipass_superimpose_problem (branch%ele(i_ele), super_ele, err_flag); if (err_flag) return
      call string_trim(super_ele_saved%name, super_ele_saved%name, ix)
      super_ele%name = super_ele_saved%name(:ix)            
      call add_superimpose (lat, super_ele, branch%ix_branch, err_flag, super_ele_out, &
                  save_null_drift = .true., create_jumbo_slave = pele%create_jumbo_slave)
      if (err_flag) bp_com%error_flag = .true.
      call control_bookkeeper (lat, super_ele_out)
    endif

    !---------------------

    call s_calc (lat)

    have_inserted = .true.   

  enddo ele_loop

  if (.not. have_inserted) exit

enddo

end subroutine add_all_superimpose

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine compute_super_lord_s (ref_ele, super_ele, pele)

implicit none

type (lat_struct), target :: lat
type (ele_struct) ref_ele, super_ele
type (ele_struct), pointer :: slave
type (parser_ele_struct) pele
type (branch_struct), pointer :: branch

integer i, ix

real(rp) s_ref_begin, s_ref_end

! Find the reference point on the element being superimposed.

super_ele%s = pele%s

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
  s_ref_begin = 1e10
  s_ref_end = 0
  do i = 1, ref_ele%n_slave
    slave => pointer_to_slave(ref_ele, i)
    s_ref_begin = min(s_ref_begin, slave%s - slave%value(l$))
    s_ref_end = max(s_ref_end, slave%s)
  enddo
case (group$)
  call parser_error ('SUPERPOSING: ' // super_ele%name, 'UPON GROUP' // pele%ref_name)
  return
case default
  s_ref_begin = ref_ele%s - ref_ele%value(l$)
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

end subroutine compute_super_lord_s

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine check_for_multipass_superimpose_problem (ref_ele, super_ele, err_flag)
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

subroutine check_for_multipass_superimpose_problem (ref_ele, super_ele, err_flag)

implicit none

type (ele_struct) ref_ele, super_ele
type (ele_struct), pointer :: ele1, ele2
type (branch_struct), pointer :: branch
real(rp) eps
logical err_flag
integer ix1, ix2


!

branch => ref_ele%branch
err_flag = .true.
eps = bmad_com%significant_length

ix1 = element_at_s (branch%lat, super_ele%s - super_ele%value(l$) + eps, .true., ref_ele%ix_branch, err_flag)
ele1 => branch%ele(ix1)
if (ele1%slave_status == super_slave$) ele1 => pointer_to_lord(ele1, 1)

ix2 = element_at_s (branch%lat, super_ele%s - eps, .false., ref_ele%ix_branch, err_flag)
ele2 => branch%ele(ix2)
if (ele2%slave_status == super_slave$) ele2 => pointer_to_lord(ele2, 1)

if (ref_ele%slave_status == multipass_slave$) then
  if (ele1%slave_status /= multipass_slave$ .or. ele2%slave_status /= multipass_slave$) then
    call parser_error ('SUPERIMPOSE OF: ' // super_ele%name, &
         'USES MULTIPASS REFERENCE ELEMENT BUT OFFSET PLACES IT OUT OF THE MULTIPASS REGION!')
    return
  endif

else
  if (ele1%slave_status == multipass_slave$ .or. ele2%slave_status == multipass_slave$) then
    call parser_error ('SUPERIMPOSE OF: ' // super_ele%name, &
         'USES NON-MULTIPASS REFERENCE ELEMENT BUT OFFSET PLACES IT IN A MULTIPASS REGION!')
    return
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
character(40) str, name
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
    if (name == ' ') then  
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
      call get_sequence_args (name, this_ele%actual_arg, delim, err_flag)
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
  plat%ele(i)%s       = 0
  plat%ele(i)%create_jumbo_slave = .false.
enddo

end subroutine allocate_plat

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_add_lord (in_lat, n2, plat, lat)
!
! Subroutine to add overlay, group, and girder lords to the lattice.
! For overlays and groups: If multiple elements have the same name then 
! use all of them.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_add_lord (in_lat, n2, plat, lat)

implicit none

type (lat_struct), target :: in_lat, lat
type (ele_struct), pointer :: lord, lord2, slave, ele, g_lord, g_slave0, g_slave1
type (parser_lat_struct), target :: plat
type (parser_ele_struct), pointer :: pele
type (control_struct), allocatable :: cs(:)
type (branch_struct), pointer :: branch

integer i, ic, ig, n, n2, k, k2, ix, j, ib, ie, n_list, ix2, ns, ixs, ii, ix_end
integer n_match, n_slave, nn
integer ix_lord, k_slave, ix_ele_now, ix_girder_end, ix_super_lord_end
integer, allocatable :: r_indexx(:), ix_ele(:), ix_branch(:)

character(40), allocatable :: name_list(:)
character(40) name, input_slave_name, attrib_name, missing_slave_name

logical err, slave_not_in_lat, created_girder_lord

! Setup...
! in_lat has the lords that are to be added to lat.
! We add an extra 1000 places to the arrays to give us some overhead.

n = n2 + 1000
do i = 0, ubound(lat%branch, 1)
  n = n + lat%branch(i)%n_ele_max
enddo

allocate (name_list(n), ix_ele(n), ix_branch(n), r_indexx(n))

n_list = 0
do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  n = branch%n_ele_max
  ix2 = n_list + n
  name_list(n_list+1:ix2) = branch%ele(1:n)%name
  ix_ele(n_list+1:ix2)    = branch%ele(1:n)%ix_ele
  ix_branch(n_list+1:ix2) = branch%ele(1:n)%ix_branch
  n_list = ix2
enddo

call indexx (name_list(1:n_list), r_indexx(1:n_list)) ! get sorted list

! loop over elements

main_loop: do n = 1, n2

  lord => in_lat%ele(n)  ! next lord to add
  pele => plat%ele(lord%ixx)
  
  !-----------------------------------------------------
  ! overlay and groups

  select case (lord%key)
  case (overlay$, group$)
 
    call new_control (lat, ix_lord)  ! get index in lat where lord goes
    lat%ele(ix_lord) = lord

    ! Find where the slave elements are. 
    ! If a slave element is not in lat but is in in_lat then the slave has 
    ! not yet been used in the lattice list. 
    ! In this case, do not add the lord to the lattice.

    ! First count the number of slaves
    n_match = 0
    do i = 1, lord%n_slave
      call find_indexx (pele%name(i), name_list, r_indexx, n_list, k, k2, n_match = nn)
      n_match = n_match + nn
    enddo

    if (allocated(cs)) then
      if (size(cs) < n_match) deallocate(cs)
    endif
    if (.not. allocated(cs)) allocate (cs(n_match))

    j = 0 ! number of slaves found
    slave_not_in_lat = .false.  ! Is there a slave that is not in the lattice?

    do i = 1, lord%n_slave

      name = pele%name(i)
      call find_indexx (name, name_list, r_indexx, n_list, k, k2)

      if (k == 0) then
        slave_not_in_lat = .true.
        missing_slave_name = name
      endif

      if ((k == 0 .and. j > 0) .or. (k > 0 .and. slave_not_in_lat) .or. &
          (k == 0 .and. all(in_lat%ele(1:n2)%name /= name))) then
        call parser_error ('CANNOT FIND SLAVE FOR: ' // lord%name, &
                      'CANNOT FIND: '// missing_slave_name, pele = pele)
        lat%n_ele_max = lat%n_ele_max - 1 ! Undo new_control call
        cycle main_loop
      endif

      if (k == 0) cycle

      ! There might be more than 1 element with %name = name. 
      ! Loop over all elements whose name matches name.
      ! Put the info into the cs structure.

      do 
        j = j + 1
        k = r_indexx(k2)
        cs(j)%coef = pele%coef(i)
        cs(j)%ix_slave = ix_ele(k)
        cs(j)%ix_branch = ix_branch(k)
        cs(j)%ix_lord = -1             ! dummy value
        attrib_name = pele%attrib_name(i)
        if (attrib_name == blank_name$) attrib_name = lord%component_name
        slave => pointer_to_ele (lat, ix_ele(k), ix_branch(k))
        ix = attribute_index(slave, attrib_name)
        ! If attribute not found it may be a special attribute like accordian_edge$.
        ! A special attribute will have ix > num_ele_attrib$
        if (ix < 1 .and. lord%key == group$) then
          ix = attribute_index(lord, attrib_name)
          if (ix <= num_ele_attrib$) ix = 0  ! Mark as not valid
        endif
        cs(j)%ix_attrib = ix
        if (ix < 1) then
          call parser_error ('IN OVERLAY OR GROUP ELEMENT: ' // lord%name, &
                        'ATTRIBUTE: ' // attrib_name, &
                        'IS NOT A VALID ATTRIBUTE OF: ' // slave%name, &
                        pele = pele)
          cycle main_loop
        endif
        k2 = k2 + 1
        if (k2 > n_list) exit
        k = r_indexx(k2)
        slave => pointer_to_ele (lat, ix_ele(k), ix_branch(k))
        ! exit loop if no more matches
        if (slave%name /= name) exit 
      enddo

    enddo

    lord%n_slave = j

    ! If the lord has no slaves then discard it

    if (j == 0) then
      lat%n_ele_max = lat%n_ele_max - 1 ! Undo new_control call
      cycle main_loop
    endif

    ! put the element name in the list r_indexx list

    call find_indexx (lord%name, name_list, r_indexx, n_list, k, add_to_list = .true.)
    n_list = n_list + 1
    ix_ele(k) = ix_lord
    ix_branch(k) = 0

    do ii = 1, n_list-1
      if (name_list(r_indexx(ii)) > name_list(r_indexx(ii+1))) then
        call parser_error ('PARER_ADD_LORD INTERNAL ERROR!')
        if (global_com%exit_on_error) call err_exit
      endif
    enddo

    ! create the lord

    ns = lord%n_slave

    select case (lord%key)
    case (overlay$)
      call create_overlay (lat, ix_lord, lord%component_name, cs(1:ns), err)
    case (group$)
      call create_group (lat, ix_lord, cs(1:ns), err)
    end select
    if (err) call parser_error ('ELEMENT OR GROUP: ' // lord%name, &
                           'IS TRYING TO CONTROL AN ATTRIBUTE THAT IS NOT FREE TO VARY!', &
                           pele = pele)

  !-----------------------------------------------------
  ! girder
  ! Create an girder element for each lattice section that matches the slave list names.
  ! If no lattice element names match any names on the girder slave list then assume 
  ! this girder is for a different lattice and ignore the girder. If some do match and 
  ! some don't then flag this as an error.

  case (girder$)

    ! Loop over all elements in the lattice.

    if (allocated(cs)) deallocate(cs)
    allocate (cs(lord%n_slave))

    created_girder_lord = .false.

    branch_loop: do ib = 0, ubound(lat%branch, 1)
      branch => lat%branch(ib)
      ix_girder_end = -1

      ! Loop over all possible first elements

      ele_loop: do ie = 1, branch%n_ele_track

        if (ie <= ix_girder_end) cycle

        ! Loop over girder slave list and see if this section matches.

        ix_ele_now = ie
        ix_super_lord_end = -1   ! Not in a super_lord
        ixs = 1       ! Index of girder slave element we are looking for
        n_slave = 0   ! Number of actual slaves found.

        ! loop over all girder slaves and see if lattice eles match slaves.

        slave_loop: do

          if (ix_ele_now > branch%n_ele_track) then
            ix_ele_now = ix_ele_now - branch%n_ele_track
          endif

          if (ixs > lord%n_slave) exit
          input_slave_name = pele%name(ixs)

          ele => pointer_to_ele (lat, ix_ele_now, ib)

          if (girder_match_this(ele)) cycle

          ! Here if no match. 
          ! If not previously in a super lord then there is no overall match to the slave list.

          if (ele%slave_status == super_slave$ .and. ix_super_lord_end > -1) then
            ix_ele_now = ix_ele_now + 1
            cycle
          endif

          ix_super_lord_end = -1  ! Not in a super_lord

          ! No match to the slave list. If a marker or drift then ignore except 
          ! if this is the first slave in which case there is no match.

          if ((ele%key == marker$ .or. ele%key == drift$) .and. ixs > 1) then
            ix_ele_now = ix_ele_now + 1
            cycle
          endif

          ! Match failed. Start again 

          cycle ele_loop

        enddo slave_loop

        ! create the girder element

        if (n_slave == 0) then
          call parser_error ('LIST OF GIRDER SLAVES IN LATTICE FILE DOES NOT INCLUDE A NON-DRIFT ELEMENT: ' // &
                              lord%name, warn_only = .true.)
          cycle main_loop
        endif

        call new_control (lat, ix_lord)
        call create_girder (lat, ix_lord, cs(1:n_slave), lord)
        created_girder_lord = .true.

        ix_girder_end = cs(lord%n_slave)%ix_slave
        if (ix_girder_end < ie) exit ele_loop  ! And on to next branch

      enddo ele_loop

    enddo branch_loop

    if (.not. created_girder_lord) then
      call parser_error ('FAILED TO FIND REGION IN LATTICE FOR CRATING GIRDER: ' // &
                          lord%name, warn_only = .true.)
    endif

  end select

enddo main_loop

! cleanup

deallocate (r_indexx)
deallocate (name_list)

!-------------------------------------------------------------------------
contains

recursive function girder_match_this (ele) result (is_matched)

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord, slave1, slave2

integer ii
logical is_matched

! Try to match

is_matched = .false.

if (match_wild(ele%name, input_slave_name)) then

  is_matched = .true.
  if (ele%key /= drift$) then
    n_slave = n_slave + 1
    cs(n_slave)%ix_slave  = ele%ix_ele
    cs(n_slave)%ix_branch = ele%ix_branch
  endif
  ixs = ixs + 1  ! Next girder slave
  ix_ele_now = ix_ele_now + 1

  ! If a super_lord the logic here is complicated by the fact that 
  ! elements can, for example, be completely contained within another element.

  if (ele%lord_status == super_lord$) then
    slave2 => pointer_to_slave(ele, ele%n_slave)
    ix_super_lord_end = ix_far_index(ix_ele_now-1, ix_super_lord_end, slave2%ix_ele)
  endif

  ! If match to a girder then move pointers to element after latst girder slave

  if (ele%key == girder$) then
    call find_element_ends (ele, slave1, slave2)
    ix_ele_now = slave2%ix_ele + 1
  endif

  return
endif

! Since ele does not match, look for match at a lord of this element

do ii = 1, ele%n_lord
  lord => pointer_to_lord (ele, ii)
  if (lord%key == group$ .or. lord%key == overlay$) cycle
  is_matched = girder_match_this(lord)
  if (is_matched) return
enddo

end function girder_match_this

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

real(rp) angle, rr
integer n
logical kick_set, length_set, set_done, err_flag

! Wall3d init.

if (associated(ele%wall3d)) then
  call wall3d_initializer (ele%wall3d, err_flag)
  if (err_flag) then
    call parser_error ('WALL INIT ERROR FOR ELEMENT: ' // ele%name)
    return
  endif
endif

! Aperture init

if (ele%aperture_type == auto_aperture$) then
  call aperture_bookkeeper (ele)
endif

!

kick_set = (ele%value(hkick$) /= 0) .or. (ele%value(vkick$) /= 0)

select case (ele%key)

! Convert rbends to sbends and evaluate G if needed.
! Needed is the length and either: angle, G, or rho.

case (sbend$, rbend$, sad_mult$) 

  if (ele%key /= sad_mult$) ele%sub_key = ele%key  ! Save sbend/rbend input type.
  angle = ele%value(angle$) 

  ! Only one of b_field, g, or rho may be set.

  if (ele%value(b_field$) /= 0 .and. ele%key == rbend$) call parser_error &
          ("B_FIELD NOT SETTABLE FOR AN RBEND (USE AN SBEND INSTEAD): " // ele%name)

  if (ele%value(b_field$) /= 0 .and. ele%value(g$) /= 0) call parser_error &
          ('BOTH G AND B_FIELD SET FOR A BEND: ' // ele%name)

  if (ele%value(b_field$) /= 0 .and. ele%value(rho$) /= 0) call parser_error &
          ('BOTH RHO AND B_FIELD SET FOR A BEND: ' // ele%name)

  if (ele%value(g$) /= 0 .and. ele%value(rho$) /= 0) &
            call parser_error ('BOTH G AND RHO SPECIFIED FOR BEND: ' // ele%name)

  ! if rho is set then this gives g

  if (ele%value(rho$) /= 0) ele%value(g$) = 1 / ele%value(rho$)

  ! If g and angle are set then this determines l

  if (ele%value(g$) /= 0 .and. angle /= 0) then
    if (ele%value(l$) /= 0) call parser_error ('ANGLE, G/RHO, AND L SPECIFIED FOR BEND: ' // ele%name)
    ele%value(l$) = angle / ele%value(g$)
  endif

  ! Convert an rbend to an sbend

  if (ele%key == rbend$) then

    if (ele%value(l$) == 0) then
      if (ele%value(l_chord$) == 0) then
        angle = 0
      elseif (angle /= 0) then
        ele%value(l$) = ele%value(l_chord$) * angle / (2 * sin(angle/2))
      elseif (ele%value(g$) /= 0) then
        angle = 2 * asin(ele%value(l_chord$) * ele%value(g$) / 2)
        ele%value(l$) = ele%value(l_chord$) * angle / (2 * sin(angle/2))
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
    ele%value(gradient$) = ele%value(voltage$) / ele%value(l$)
  endif

! 

case (rfcavity$) 

  if (ele%value(rf_frequency$) /= 0 .and. ele%value(harmon$) /= 0) call parser_error &
              ('BOTH RF_FREQUENCY AND HARMON SET FOR RFCAVITY: ' // ele%name, &
               'SETTING OF HARMON WILL BE IGNORED!', warn_only = .true.)

! for a periodic_type wiggler n_pole is a dependent attribute

case (wiggler$, undulator$)
  if (ele%sub_key == periodic_type$) then

    if (ele%value(l_pole$) == 0 .and. ele%value(n_pole$) /= 0) then
      ele%value(l_pole$) = ele%value(l$) / ele%value(n_pole$) 
    endif

  endif

! check for inconsistancies

case (solenoid$)
  if (ele%field_master .and. (ele%value(ks$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH (KS, HKICK, ETC.) AND FIELD SET FOR A SOLENOID.')

case (sol_quad$)
  if (ele%field_master .and. (ele%value(ks$) /= 0 .or. &
                            ele%value(k1$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH (K1, HKICK, ETC.) AND FIELD SET FOR A SOL_QUAD.')

case (quadrupole$)
  if (ele%field_master .and. (ele%value(k1$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH (K1, HKICK, ETC.) AND FIELD SET FOR A QUAD.')

case (sextupole$)
  if (ele%field_master .and. (ele%value(k2$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH (K2, HKICK, ETC.) AND FIELD SET FOR A SEXTUPOLE.')

case (octupole$)
  if (ele%field_master .and. (ele%value(k3$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH (K3, HKICK, ETC.) AND FIELD SET FOR A OCTUPOLE.')

case (hkicker$, vkicker$)
  if (ele%field_master .and. (ele%value(kick$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM: ' // ele%name, &
       'BOTH STRENGTH AND BL_KICK SET FOR A H/VKICKER.')

case (e_gun$)
  if (ele%value(gradient$) == 0 .and. ele%value(l$) /= 0) ele%value(gradient$) = ele%value(voltage$) / ele%value(l$)

end select

! set ds_step if not already set.

if (attribute_index(ele, 'DS_STEP') > 0) then  ! If this is an attribute for this element...
  if (ele%value(num_steps$) > 0) then
    ele%value(ds_step$) = abs(ele%value(l$) / ele%value(num_steps$))
  elseif (ele%value(ds_step$) == 0) then
    if ((ele%key == wiggler$ .or. ele%key == undulator$) .and. ele%value(l_pole$) /= 0) then
      ele%value(ds_step$) = ele%value(l_pole$) / 10
    else
      ele%value(ds_step$) = bmad_com%default_ds_step
    endif
  endif
endif

end subroutine settable_dep_var_bookkeeping 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine form_digested_bmad_file_name (lat_file, digested_file, full_lat_file)
!
! Subroutine to form the standard name of the Bmad digested file. 
! The standard digested file name has the suffix '_digested' added to the file name.
!
! Modules needed:
!   use bmad_parser_mod
!
! Input:
!   lat_file -- Character(200): Input lattice file name.
!
! Output:
!   digested_file -- Character(200): Name of the digested file.
!   full_lat_file  -- Character(200), optional: Input lattice file name with full directory
!
!-

subroutine form_digested_bmad_file_name (lat_file, digested_file, full_lat_file)

character(*) lat_file, digested_file
character(*), optional :: full_lat_file
character(200) full_name

integer ix

!

call fullfilename (lat_file, full_name)
inquire (file = full_name, name = full_name)  ! full input file_name
if (present (full_lat_file)) full_lat_file = full_name
write (digested_file, '(2a, i0)') trim(full_name), '.digested', bmad_inc_version$ 

end subroutine form_digested_bmad_file_name

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_add_branch (fork_ele, lat, sequence, in_name, in_indexx, &
!                                                        seq_name, seq_indexx, in_lat, plat)
!
! Subroutine to do line expansion.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

recursive subroutine parser_add_branch (fork_ele, lat, sequence, in_name, in_indexx, &
                                                   seq_name, seq_indexx, in_lat, plat, created_new_branch)

implicit none

type (lat_struct), target :: lat, in_lat
type (parser_lat_struct) plat
type (ele_struct) fork_ele
type (ele_struct), pointer :: target_ele
type (seq_struct), target :: sequence(:)
type (branch_struct), pointer :: branch

integer, allocatable :: seq_indexx(:), in_indexx(:)
integer i, j, nb, n_ele_use, n, ix, key

character(*), allocatable ::  in_name(:), seq_name(:)
character(40) name

logical created_new_branch

!

created_new_branch = .true.

if (fork_ele%value(new_branch$) == 0) then ! Branch back if
  do i = 0, ubound(lat%branch, 1) - 1
    branch => lat%branch(i)
    if (branch%name /= fork_ele%component_name) cycle
    fork_ele%value(ix_to_branch$) = i
    created_new_branch = .false.
  enddo
endif

if (created_new_branch) then
  call parser_expand_line (lat, fork_ele%component_name, sequence, in_name, &
                                in_indexx, seq_name, seq_indexx, in_lat, n_ele_use)

  nb = ubound(lat%branch, 1)
  fork_ele%value(ix_to_branch$) = nb
  branch => lat%branch(nb)

  branch%ix_from_branch     = fork_ele%ix_branch
  branch%ix_from_ele        = fork_ele%ix_ele
endif

name = plat%ele(fork_ele%ixx)%ele_name
nullify(target_ele)

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

key = target_ele%key
if (key /= marker$ .and. key /= fork$ .and. key /= photon_fork$ .and. &
    key /= fiducial$ .and. key /= beginning_ele$) then
  call parser_error('TO_ELEMENT: ' // name, 'FOR FORK ELEMENT: ' // fork_ele%name, &
                    'IS NOT A ZERO-LENGTH MARKER-LIKE ELEMENT')
endif

end subroutine parser_add_branch

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_expand_line (lat, use_name, sequence, in_name, &
!                       in_indexx, seq_name, seq_indexx, in_lat, n_ele_use, allow_end_marker)
!
! Subroutine to do line expansion.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_expand_line (lat, use_name, sequence, in_name, &
                       in_indexx, seq_name, seq_indexx, in_lat, n_ele_use, allow_end_marker)

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
integer iseq_tot, i_lev, i_use, n0_multi, n_ele_use, n_max
integer i, j, k, n, ix, ix_multipass, ix_branch, flip

character(*), allocatable ::  in_name(:), seq_name(:)
character(*) use_name
character(40) name

logical, optional :: allow_end_marker

! find line corresponding to the "use" statement.

iseq_tot = size(seq_indexx)
n_max = in_lat%n_ele_max
allocate (used_line(n_max))

call find_indexx (use_name, seq_name, seq_indexx, iseq_tot, i_use)
if (i_use == 0) then
  call parser_error ('CANNOT FIND DEFINITION OF LINE IN "USE" STATEMENT: ' // use_name, ' ')
  return
endif

if (sequence(i_use)%type /= line$) then
  call parser_error ('NAME IN "USE" STATEMENT IS NOT A LINE!', ' ')
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
      call find_indexx2 (name, in_name, in_indexx, 0, n_max, j)
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

if (stack(1)%multipass) then
  call parser_error ('"USE"D LINE FOR LATTICE EXPANSION IS MARKED MULTIPASS!')
  if (global_com%exit_on_error) call err_exit
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
    call find_indexx2 (name, in_name, in_indexx, 0, n_max, j)
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
                      'IN THE LINE/LIST: ' // seq%name, seq)
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

    if (stack(i_lev)%multipass .and. .not. stack(i_lev-1)%multipass) then
      ix_multipass = 1
      n0_multi = n_ele_use + 1
    endif

  case default
    call parser_error ('INTERNAL SEQUENCE ERROR!')

  end select

enddo line_expansion

! Stop here if there has been an error

if (bp_com%error_flag) return

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

if (nint(bp_com%param_ele%value(no_end_marker$)) == 0 .and. logic_option(.true., allow_end_marker)) then
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
    print *
    print *, 'Var #', i
    print *, 'Name: ', bp_com%var(i)%name
    print *, 'Value:', bp_com%var(i)%value
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
! parse_rf_map (grid, ele, lat, delim, delim_found, err_flag)
!
! Subroutine to parse a "map = {}" construct
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

subroutine parse_rf_map (map, ele, lat, delim, delim_found, err_flag)

implicit none

type (em_field_map_struct), pointer :: map
type (ele_struct), target :: ele
type (lat_struct) lat

real(rp), allocatable :: array(:)

complex(rp), pointer :: c_ptr(:)

integer ix_word, i_term

character(1) delim, delim2
character(40) word, word2, name, attrib_name

logical err_flag, delim_found

!

err_flag = .true.
if (.not. associated(map)) allocate (map)
allocate (array(1024))

! Expect {

call get_next_word (word, ix_word, '{', delim, delim_found, call_check = .true.)
if ((word /= '') .or. (delim /= '{')) then
  call parser_error ('NO { SIGN FOUND IN MAP DEFINITION',  &
                    'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
  return
endif

!

do

  ! Read attriubute
  call get_next_word (attrib_name, ix_word, '{}=,()', delim, delim_found)

  select case (attrib_name)

  case ('ELE_ANCHOR_PT')
    ! Expect "<component> = "
    if (delim /= '=') then
      call parser_error ('NO "=" SIGN FOUND AFTER MAP ' // attrib_name,  &
                         'IN MAP STRUCTURE IN ELEMENT: ' // ele%name)
      return
    endif
    call get_next_word (word2, ix_word, ',}', delim, delim_found)

    ! Evaluate string into integer.

    call match_word(word2, anchor_pt_name(1:), map%ele_anchor_pt, can_abbreviate = .false., matched_name = name)
  
    if (name == '') then
      call parser_error ('UNKNKOWN MAP ' // trim(word) // ': ' // word2, &
                         'FOUND IN MODE MAP DEFINITION FOR ELEMENT: ' // ele%name)
      return        
    endif      
      
  case ('E_COEF_RE', 'E_COEF_IM', 'B_COEF_RE', 'B_COEF_IM')

    ! Expect "("
    call get_next_word (word, ix_word, ',({', delim, delim_found)
    if (word /= '' .or. delim /= '(') then
      call parser_error ('NO "(" FOUND AFTER "' // trim(attrib_name) // ' =" ', &
                           'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
      return
    endif

    ! Read list of values.
    call re_allocate(array, 1024, .false.)
    do i_term = 1, 100000
      call get_next_word (word, ix_word, '{},()', delim, delim_found)
      if ((delim /= ',' .and. delim /= ')') .or. .not. is_real(word)) then
        call parser_error ('ERROR PARSING ARRAY FOR: ' // word2, &
                             'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
        return
      endif
      if (i_term > size(array)) call re_allocate(array, 2*size(array))
      read (word, *) array(i_term)
      if (delim == ')') exit
    enddo

    if (allocated(map%term)) then
      if (size(map%term) /= i_term) then
        call parser_error ('ARRAY SIZE MISMATCH FOR: ' // word2, &
                           'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
        return
      endif
    else
      allocate(map%term(i_term))
    endif

    select case (attrib_name)
    case ('E_COEF_RE', 'E_COEF_IM'); c_ptr => map%term%e_coef 
    case ('B_COEF_RE', 'B_COEF_IM'); c_ptr => map%term%b_coef
    end select

    if (attrib_name(8:9) == 'RE') then
      if (any(real(c_ptr) /= 0)) then
        call parser_error ('DUPLICATE ARRAY FOR: ' // attrib_name, &
                           'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
        return
      endif
      c_ptr = c_ptr + array(1:i_term)

    else
      if (any(aimag(c_ptr) /= 0)) then
        call parser_error ('DUPLICATE ARRAY FOR: ' // attrib_name, &
                           'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
        return
      endif
      c_ptr = c_ptr + i_imaginary * array(1:i_term)
    endif

    ! Expect "," or "}"
    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
    if (word /= '' .or. (delim /= ',' .and. delim /= '}')) then
      call parser_error ('BAD ' // trim(attrib_name) // ' = (...) CONSTRUCT', &
                           'FOUND IN MODE DEFINITION IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
      return
    endif

  case ('DZ');            
    if (.not. associated(map)) allocate (map)
    call evaluate_value (trim(ele%name), map%dz, lat, delim, delim_found, err_flag, ',}')

  end select

  ! Possible "}" is end of mode
  if (delim == '}') exit

enddo

! Get final separator after grid construct.
 
call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

if (ix_word /= 0) then
    call parser_error ('UNKNOWN TEXT FOUND AFTER MODE GRID DEFINITION: ' // word, &
                         'FOR ELEMENT: ' // ele%name)
    return 
endif

! Deallocate

deallocate(array)
err_flag = .false.

end subroutine parse_rf_map

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! parse_rf_grid (grid, ele, lat, delim, delim_found, err_flag)
!
! Subroutine to parse a "grid = {}" construct
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

subroutine parse_rf_grid(grid, ele, lat, delim, delim_found, err_flag)

type grid_pt_struct
  integer :: ix(3) = [1, 1, 1]
  complex(rp) :: field(6) = 0
end type


type (em_field_grid_struct), pointer :: grid
type (ele_struct) :: ele
type (ele_struct), pointer :: bele
type (lat_struct),  target :: lat
type (branch_struct), pointer :: branch
type (grid_pt_struct), allocatable :: array(:), array2(:)

character(1) delim, delim2
character(40) :: word, word2, name

integer ix_word, ix_word2
integer pt_counter, n, i, ib, ie, im, ix0, ix1, iy0, iy1, iz0, iz1
integer grid_dim,  num_dr, num_r0

logical delim_found, delim_found2, err_flag, err_flag2

!

if (.not. associated(grid)) allocate (grid)

! Expect {

call get_next_word (word, ix_word, '{', delim, delim_found, call_check = .true.)
if ((word /= '') .or. (delim /= '{')) then
  call parser_error ('NO { SIGN FOUND IN GRID DEFINITION',  &
                    'IN FIELD STRUCTURE IN ELEMENT: ' // ele%name)
  return
endif

! Set %file to be the last called file with 

write(grid%file, '(2a, i0)') trim(bp_com%current_file%full_name),  ':', bp_com%current_file%i_line

! Init

allocate(array(1024))
pt_counter = 0
err_flag = .true.

do    

  ! Read attriubute
  call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

  select case (word)

  case ('TYPE', 'ELE_ANCHOR_PT')
    ! Expect "<component> = "
    if (delim /= '=') then
      call parser_error ('NO "=" SIGN FOUND AFTER GRID ' // word,  &
                         'IN GRID STRUCTURE IN ELEMENT: ' // ele%name)
      return
    endif
    call get_next_word (word2, ix_word, ',}', delim, delim_found)
    ! Check to see if this is a valid type by checking against em_grid_type_name(:)

    if (word == 'TYPE') then
      call match_word(word2, em_grid_type_name, grid%type, can_abbreviate = .false., matched_name = name)
    else
      call match_word(word2, anchor_pt_name(1:), grid%ele_anchor_pt, can_abbreviate = .false., matched_name = name)
    endif
  
    if (name == '') then
      call parser_error ('UNKNKOWN GRID ' // trim(word) // ': ' // word2, &
                         'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
      return        
    endif      
      
  case ('R0')       
    ! Expect "<component> = "
    if (delim /= '=') then
      call parser_error ('NO "=" SIGN FOUND AFTER GRID R0',  &
                      'IN GRID STRUCTURE IN ELEMENT: ' // ele%name)
      return
    endif
    ! expect ( 1. ) or (1. , 2.) or (1., 2., 3.)
    if (.not. parse_real_list (lat, trim(ele%name) // ' GRID', grid%r0, .false.)) return
    ! Expect , or }
    call get_next_word (word, ix_word, ',}', delim, delim_found)     
    if (word /= '') then
      call parser_error ('BAD INPUT AFTER R0 DEFINITION: ' // word , &
                                 'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
      return
    end if

    case ('DR')      
    ! Expect "<component> = "
    if (delim /= '=') then
      call parser_error ('NO "=" SIGN FOUND AFTER GRID DR',  &
                      'IN GRID STRUCTURE IN ELEMENT: ' // ele%name)
      return
    endif      
    ! expect ( 1.) or (1. , 2.) or (1., 2., 3.)
    if (.not. parse_real_list (lat, trim(ele%name) // ' GRID', grid%dr, .false.)) return
    call get_next_word (word, ix_word, ',}', delim, delim_found)     
    if (word /= '') then
      call parser_error ('BAD INPUT AFTER DR DEFINITION: ' // word , &
                                 'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
      return
    end if

  case ('PT')
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
    if (.not. parse_integer_list (trim(ele%name) // ' GRID PT', array(pt_counter)%ix, .false.)) return
      
    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
    call get_next_word (word2, ix_word2, '{}=,()', delim2, delim_found2)
    if ((word /= '') .or. (word2 /= '') .or. (delim /= '=') .or. (delim2 /= '(')) then
      call parser_error ('BAD GRID PT CONSTRUCT, NO  = "(" ', &
                 'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
      return
    end if
    ! Get as many field components as listed
    do i = 1, 6
      call parse_complex_component(array(pt_counter)%field(i), delim, err_flag2)
      if (err_flag2) return
      if (delim == ')') exit
      if (delim /= ',') then
        call parser_error ('BAD GRID PT CONSTRUCT, NO "," BETWEEN FIELD COMPONENTS', &
            'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
        return
      end if    
    end do
    ! Expect , or }
    call get_next_word (word, ix_word, ',}', delim, delim_found) 
    if (word /= '') then
      call parser_error ('BAD INPUT AFTER PT DEFINITION: ' // word , &
                           'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
      return
    end if
  
  case default
    call parser_error ('UNKNOWN GRID ATTRIBUTE: ' // word, &
                         'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
    return 
        
  end select 

  if (delim == '}') exit   

enddo

! Get final separator after grid construct.
 
call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

if (ix_word /= 0) then
    call parser_error ('UNKNOWN TEXT FOUND AFTER MODE GRID DEFINITION: ' // word, &
                         'FOR ELEMENT: ' // ele%name)
    return 
endif

! Clear pts

if (allocated(grid%pt)) deallocate(grid%pt)

! Allocate grid for different dimensions

grid_dim = em_grid_dimension(grid%type)

if (grid_dim < 1 .or. grid_dim > 3) then
  call parser_error ('BAD GRID DIMENSION', &
             'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
  return
endif

ix0 = minval(array%ix(1))
ix1 = maxval(array%ix(1))
iy0 = minval(array%ix(2))
iy1 = maxval(array%ix(2))
iz0 = minval(array%ix(3))
iz1 = maxval(array%ix(3))

allocate(grid%pt(ix0:ix1, iy0:iy1, iz0:iz1))

! Assign grid values
do i = 1, pt_counter
  ix1 = array(i)%ix(1)
  iy1 = array(i)%ix(2)
  iz1 = array(i)%ix(3)
  grid%pt(ix1, iy1, iz1)%E(1:3) = array(i)%field(1:3)
  grid%pt(ix1, iy1, iz1)%B(1:3) = array(i)%field(4:6)
end do

! Clear temporary array

deallocate(array)

! Check if grid data has already been read in for another element.
! If so, save space by pointing to the existing grid.

branch_loop: do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do ie = 1, branch%n_ele_max
    bele => branch%ele(ie)
    if (.not. associated(bele%em_field)) cycle    
    if (bele%ix_ele == ele%ix_ele .and. bele%ix_branch == ele%ix_branch) cycle
    do im = 1, size(bele%em_field%mode)
      if (.not. associated(bele%em_field%mode(im)%grid)) cycle
      if (bele%em_field%mode(im)%grid%file /= grid%file) cycle
      deallocate(grid)
      grid => bele%em_field%mode(im)%grid
      grid%n_link = grid%n_link + 1        
      exit branch_loop
    end do
  enddo 
enddo branch_loop

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
  complex_component = cmplx(x, 0.0_rp)
else if (delim == '(') then
  call get_next_word (word, ix_word, ',', delim, delim_found)
  call get_next_word (word2, ix_word2, ')', delim2, delim_found2)
  ! Expect: real, real ) 
  if ((.not. is_real(word)) .or. (.not. is_real(word2)) &
    .or. (delim /= ',') .or. (delim2 /= ')') &
    .or. (.not. delim_found) .or. (.not. delim_found2)) then
    call parser_error ('BAD COMPLEX COMPONENT CONSTRUCT', &
                         'FOUND IN MODE GRID DEFINITION FOR ELEMENT: ' // ele%name)
    return
   end if
   !
   read (word, *) x
      read (word2, *) y
      complex_component = cmplx(x,y)
      ! Look for "," or end of list ")"
   call get_next_word (word, ix_word, ',)', delim, delim_found)
end if

err_flag2 = .false.

end subroutine parse_complex_component

end subroutine parse_rf_grid


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function parse_integer_list (err_str, int_array, exact_size, open_delim, 
!           separator, close_delim, default_value) result (is_ok)
!
! Routine to parse a list of integers of the form:
!    open_delim integer_1 separator integer_2 . . . close_delim
! Example:   "(1.2, 2.3, 4.4, 8.5)"
! 
! Similar to parse_integer_list2 except does not use allocatable array.
! See parse_integer_list2 for more details
!-

function parse_integer_list (err_str, int_array, exact_size, open_delim, &
                          separator, close_delim, default_value) result (is_ok)

implicit none

integer int_array(:)
integer, optional :: default_value
integer, allocatable :: vec(:)

integer num_found

character(*) err_str
character(*), optional :: open_delim, separator, close_delim

logical is_ok, exact_size

!

is_ok = .false.
if (.not. parse_integer_list2 (err_str, vec, num_found, size(int_array), &
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
! Function parse_integer_list2 (err_str, int_array, num_found, num_expected, open_delim, 
!                               separator, close_delim, default_value) result (is_ok)
!
! Routine to parse a list of integers of the form
!    open_delim integer_1 separator integer_2 . . . close_delim
! Example:   (1, 2, 4, 8) 
!
! Input:
!  err_str    -- character(*): Error string to print if there is an error. 
!  int_array -- Integer, allocatable: the array to be read in 
!
!   Optional: 
!   num_expected = 1     -- integer : number of expected arguments
!                             Used to initialize int_array
!   open_delim   = '('   -- character(1) : opening delimeter
!   separator    = ','   -- character(1) : separating character
!   close_delim  = ')'   -- character(1) : closing delimeter
!   default_value = 0    -- real(rp) : inital assignment of real_array elements
!
! Output:
!   is_ok                   -- logical: Set True if everything is ok
!   int_array(1:num_found) --integer(rp) : Array of values
!   num_found                  -- integer : number of elements
!-

function parse_integer_list2 (err_str, int_array, num_found, num_expected, open_delim, &
                              separator, close_delim, default_value) result (is_ok)

! Arguments
integer, allocatable :: int_array(:)
integer :: num_found
integer, optional :: num_expected, default_value
character(*) err_str
character(1), optional :: open_delim, close_delim, separator
logical is_ok

! Local
integer num_expect
character(1) delim, op_delim, cl_delim, sep
character(40) :: word
integer  ix_word
logical delim_found

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
call get_next_word (word, ix_word, op_delim, delim, delim_found)
if ((word /= '') .or. (delim /= op_delim)) then
  call parser_error (err_str)
  return
end if

! Initial allocation
call re_allocate(int_array, num_expected, .false.)
int_array = integer_option(0, default_value)

! counter
num_found = 0

! Get integers
do 
  call get_next_word (word, ix_word, sep // cl_delim, delim, delim_found)
  if (.not. is_integer(word) ) then 
    call parser_error ('BAD REAL NUMBER IN: ' // err_str)
    return
   end if    
  ! integer is found
  num_found = num_found + 1
  ! reallocate if needed  
  ! reallocate if needed  
  if (size(int_array) < num_found) then
    call re_allocate (int_array, 2*num_found, .false.)
    int_array(num_found:2*num_found) = integer_option(0, default_value)
  endif
  
  ! Read value
   read (word, *)  int_array(num_found) 
  
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

end module
