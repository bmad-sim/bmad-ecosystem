module bmad_parser_struct

use bmad_struct

integer, parameter :: n_parse_line = 500, n_parse_line_extended = n_parse_line + 100

! A "sequence" is a line or a list.
! The information about a sequence is stored in a seq_struct.

! A seq_struct has an array of seq_ele_struct structures.
! Each seq_ele_struct represents an individual element in a sequence and, 
! since sequences can be nested, can itself be a line or a list.

type seq_ele_struct
  character(40) name                     ! name of element, subline, or sublist
  character(40), allocatable :: actual_arg(:)
  character(40) :: tag = ''              ! tag name.
  character(40) :: slice_start = ''      ! For "my_line[start:end]" slice constructs.
  character(40) :: slice_end = ''        ! For "my_line[start:end]" slice constructs.
  integer :: type = 0                    ! LINE$, REPLACEMENT_LINE$, LIST$, ELEMENT$
  integer :: ix_ele = 0                  ! if an element: pointer to ELE array
                                         ! if a line or list: pointer to SEQ array
  integer :: ix_arg  = 0                 ! index in arg list (for replacement lines)
  integer ::rep_count = 1                ! how many copies of an element
  logical :: ele_order_reflect = .false. ! Travel through ele sequence in reverse order
  integer :: ele_orientation = 1         ! element has reverse orientation.
end type

type base_line_ele_struct
  character(40) :: name = ''             ! Name of sequence or element
  character(40) :: tag = ''              ! Tag name.
  integer :: ix_multi = 0                ! Multipass indentifier
  integer :: orientation = 1             ! Element reversed?
  integer :: ix_ele_in_in_lat = -1
  logical :: ele_order_reflect = .false. ! Part of reflection or reversed line?
end type

type seq_struct
  character(40) name                ! name of sequence
  type (seq_ele_struct), allocatable :: ele(:)   ! Elements in the sequence
  character(40), allocatable :: dummy_arg(:)
  character(40), allocatable :: corresponding_actual_arg(:)
  integer type                      ! LINE$, REPLACEMENT_LINE$ or LIST$
  integer ix_list                   ! Current index for lists
  integer :: list_upcount = 0
  integer index                     ! Alphabetical order sorted index
  character(200) :: file_name = ''  ! File where sequence is defined
  integer ix_file_line              ! Line number in file where sequence is defined
  logical multipass
  logical ptc_layout                ! Put in separate PTC layout
  logical :: active = .false.       ! Used to prevent infinite loops.
end type

! A LIFO stack structure is used to hold the list of input lattice files that are currently open.
! %parse_line_saved is used with inline calls to save the rest of the line of the current file
! while the called file is being parsed.

integer, parameter :: f_maxx = 20

type stack_file_struct
  character(200) :: full_name = ''
  character(200) :: dir = './'
  character(n_parse_line+20) :: input_line1_saved = ''
  character(n_parse_line+20) :: input_line2_saved = ''
  character(n_parse_line) :: rest_of_line_saved = ''
  character(n_parse_line) :: parse_line_saved = ''
  character(n_parse_line) :: next_chunk_saved = ''
  character(1) :: last_char_in_parse_line_saved = ''
  integer :: ios_next_chunk_saved = 0
  integer :: ios_this_chunk_saved = 0
  integer :: i_line = 0
  integer :: f_unit = 0
  logical :: inline_call_active = .false.
end type

! structure for holding the control names and pointers for superimpose and overlay elements

type parser_controller_struct ! For overlays and groups
  character(40) :: name     
  character(40) :: attrib_name
  type (expression_atom_struct), allocatable :: stack(:) ! Arithmetic expression stack
  real(rp), allocatable :: y_knot(:)
  integer n_stk 
end type

type parser_ele_struct
  type (parser_controller_struct), allocatable :: control(:)
  character(40), allocatable :: field_overlaps(:)
  character(40) :: ref_name = ''
  integer :: ix_super_ref_multipass = 0        ! Multipass index for superimpose reference element.
  character(40) :: ele_name = ''               ! For fork element or superimpose statement.
  character(40), allocatable :: names1(:)      ! Currently just used by feedback element.
  character(40), allocatable :: names2(:)      ! Currently just used by feedback element.
  character(200) :: lat_file = ''              ! File where element was defined.
  real(rp) :: offset = 0
  integer ix_line_in_file    ! Line in file where element was defined.
  integer ix_count
  integer ele_pt, ref_pt
  integer index
  logical :: superposition_command_here = .false.
  logical :: superposition_has_been_set = .false.
  logical :: wrap_superimpose = .true.
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

type bp_const_struct
  character(40) :: name = ''   ! Constant name
  real(rp) :: value = 0        ! Constant value
  integer :: index = 0         ! Constant sort index
end type

type bp_common_struct
  type (stack_file_struct) :: file(0:f_maxx) = stack_file_struct()
  type (stack_file_struct), pointer :: current_file => null()
  type (lat_struct), pointer :: old_lat => null()
  type (extra_parsing_info_struct) :: extra = extra_parsing_info_struct()
  integer :: i_file_level = 0
  integer :: num_lat_files = 0                        ! Number of files opened
  integer :: i_const_tot = 0, i_const_init = 0
  integer :: ios_next_chunk = 0
  integer :: ios_this_chunk = 0
  character(200), allocatable :: lat_file_names(:)     ! List of all files used to create lat
  ! Note: use %line2_file_name to ID line. %line1_file_name may be blank!
  character(200) :: line1_file_name = ''               ! Name of file from which %input_line1 was read
  character(200) :: line2_file_name = ''               ! Name of file from which %input_line2 was read
  character(n_parse_line_extended) :: parse_line = ''  ! Current string to be parsed.
  character(n_parse_line_extended) :: input_line1 = '' ! Line before current line. For debug messages.
  character(n_parse_line_extended) :: input_line2 = '' ! Current line. For debug messages.
  character(n_parse_line) :: rest_of_line = ''         ! Line after semicolon saved until current statement is completely parsed.
  character(n_parse_line) :: next_chunk = ''      ! Line waiting to be appended to the parse_line.
  character(1) :: last_char_in_parse_line =  ''   ! Needed for long lines read in pieces.
  ! parser_name is used by routines to tell if parsing is being done or not.
  character(40) :: parser_name = ''               ! Blank means not in bmad_parser nor bmad_parser2.
  character(100) :: last_word = ''                ! Last word to be parsed
  logical :: bmad_parser_calling = .false.        ! used for expand_lattice
  logical :: fatal_error_flag = .false.           ! Set True on fatal (must abort now) error 
  logical :: error_flag = .false.                 ! Set True on error
  logical :: input_line_meaningful = .false.   
  logical :: do_superimpose = .true.
  logical :: write_digested = .true.              ! For bmad_parser
  logical :: write_digested2 = .true.             ! For bmad_parser2
  logical :: always_parse = .false.               ! For debugging to force parsing
  logical :: input_from_file = .true.             ! Input is from a lattice file?
  logical :: inline_call_active = .false.
  logical :: print_err = .true.                   ! Print error messages?
  ! For compatibility with translated MAD files, treat undefined vars as having zero value.
  ! Note: When using the parser code for local evaluations (done by Tao), do not wnat this.
  logical :: undefined_vars_evaluate_to_zero = .true.
  logical :: use_local_lat_file = .false.
  logical :: used_line_set_by_calling_routine = .false.
  logical :: calc_reference_orbit = .false.
  logical :: detected_expand_lattice_cmd = .false.
  real :: time0 = 0, time1 = 0, time2 = 0, time3 = 0   ! For timing parsing
end type

! The const component is separated into a separate structure from bp_common_struct to get
! around a bug in gfortran where an error is generated with initialization in the form
!   bp_com = bp_common_struct()
! See GCC Bugzilla Bug report #87568.

type bp_common2_struct
  type (bp_const_struct), allocatable :: const(:)   ! Constant name
end type

!

type (bp_common_struct), save, target :: bp_com
type (bp_common2_struct), save, target :: bp_com2

end module
