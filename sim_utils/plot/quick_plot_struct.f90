module quick_plot_struct

use utilities_mod
use sim_utils_interface

integer, parameter :: white$ = 0, black$ = 1, red$ = 2, green$ = 3
integer, parameter :: blue$ = 4, cyan$ = 5, magenta$ = 6, yellow$ = 7
integer, parameter :: orange$ = 8, yellow_green$ = 9, light_green$ = 10
integer, parameter :: navy_blue$ = 11, purple$ = 12, reddish_purple$ = 13
integer, parameter :: dark_grey$ = 14, light_grey$ = 15, transparent$ = 16

character(16), parameter :: qp_color_name(-1:16) =   [ &
  'Not_Set       ', 'White         ', &
  'Black         ', 'Red           ', 'Green         ', 'Blue          ', &
  'Cyan          ', 'Magenta       ', 'Yellow        ', 'Orange        ', &
  'Yellow_Green  ', 'Light_Green   ', 'Navy_Blue     ', 'Purple        ', &
  'Reddish_Purple', 'Dark_Grey     ', 'Light_Grey    ', 'Transparent   ']

integer, parameter :: solid$ = 1, dashed$ = 2, dash_dot$ = 3
integer, parameter :: dotted$ = 4, dash_dot3$ = 5

character(16), parameter :: qp_line_pattern_name(5) = ['solid    ', &
    'dashed   ', 'dash_dot ', 'dotted   ', 'dash_dot3' ]

integer, parameter :: solid_fill$ = 1, no_fill$ = 2
integer, parameter :: hatched$ = 3, cross_hatched$ = 4

character(16), parameter :: qp_symbol_fill_pattern_name(4) = [character(16):: 'solid_fill', &
                                                          'no_fill', 'hatched', 'cross_hatched' ]

integer, parameter :: square_sym$ = 0, dot_sym$ = 1, plus_sym$ = 2, times_sym$ = 3
integer, parameter :: circle_sym$ = 4, x_symbol_sym$ = 5, triangle_sym$ = 7
integer, parameter :: circle_plus_sym$ = 8, circle_dot_sym$ = 9
integer, parameter :: square_concave_sym$ = 10, diamond_sym$ = 11
integer, parameter :: star5_sym$ = 12, triangle_filled_sym$ = 13, red_cross_sym$ = 14
integer, parameter :: star_of_david_sym$ = 15, square_filled_sym$ = 16
integer, parameter :: circle_filled_sym$ = 17, star5_filled_sym$ = 18

character(16), parameter :: qp_symbol_type_name(-1:18) = [ &
    'do_not_draw    ', 'square         ', &
    'dot            ', 'plus           ', 'times          ', 'circle         ', &
    'x_symbol       ', '---------------', 'triangle       ', 'circle_plus    ', &
    'circle_dot     ', 'square_concave ', 'diamond        ', 'star5          ', &
    'triangle_filled', 'red_cross      ', 'star_of_david  ', 'square_filled  ', &
    'circle_filled  ', 'star5_filled   ' ]

integer, parameter :: dflt_draw$ = 1, dflt_set$ = 2

real(rp), parameter :: print_page_long_len = 10.5
real(rp), parameter :: print_page_short_len = 7.8

integer, parameter :: filled_arrow_head$ = 1, outline_arrow_head$ = 2
character(16), parameter :: qp_arrow_head_type_name(2) = [character(16):: 'Filled', 'Outline']

!------------------------------------

type qp_axis_struct
  character(80) :: label = ' '
  real(rp) :: min = 0, max = 10          ! Axis min/max in data units.
  real(rp) :: tick_min = 0               ! Min tick location along axis in data units.
  real(rp) :: tick_max = 10              ! Max tick location along axis in data units.
  real(rp) :: eval_min = 0, eval_max = 0 ! For general use. Not set by quick_plot.
  real(rp) :: dtick = 2                  ! Distance between ticks. In data units. Ticks will be drawn between %min and %max.
  real(rp) :: number_offset = 0.05       ! Offset from axis line in inches.
  real(rp) :: label_offset = 0.05        ! Offset from numbers in inches.
  real(rp) :: major_tick_len = 0.10      ! In inches.
  real(rp) :: minor_tick_len = 0.06      ! In inches.
  character(16) :: label_color = 'black' ! Color of the label.
  integer :: major_div = 5               ! Actual numbrer of major divisions
  integer :: major_div_nominal = 5       ! Nominal value.
  integer :: minor_div = 0               ! 0 = auto choose.
  integer :: minor_div_max = 5           ! Max number for auto choose.
  integer :: places = 0                  ! Number of places after the decimal point to print.
  character(16) :: type = 'LINEAR'       ! Or 'LOG', or 'CUSTOM'
  character(16) :: bounds = 'GENERAL'    ! Or 'ZERO_AT_END' or 'ZERO_SYMMETRIC'
  integer :: tick_side = +1              ! +1 = Draw on the side inside the graph, 0 = both (longer tick), -1 = outside.
  integer :: number_side = -1            ! +1 = Draw to the side inside the graph, -1 = outside.
  logical :: draw_label = .true.
  logical :: draw_numbers  = .true.
end type

! x2_mirrors_x and y2_mirrors_y force the x2 or y2 axis to have the same tick markings.
! What can be still different is whether a label and/or numbering is drawn.

type qp_plot_struct
  character(80) :: title = ' '
  type (qp_axis_struct) x, y, x2, y2
  type (qp_axis_struct), pointer :: xx  ! Pointer to x or x2 active axes currently used for plotting.
  type (qp_axis_struct), pointer :: yy  ! Pointer to y or y2 active axes currently used for plotting.
  logical :: draw_box    = .true.
  logical :: draw_title  = .true.
  logical :: draw_grid   = .true.
  logical :: x2_mirrors_x = .true.   ! x2 axis same as x?
  logical :: y2_mirrors_y = .true.   ! y2 axis same as y?
  logical :: x_is_active_axis        ! True => x is the active axis. False => x2 is active.
  logical :: y_is_active_axis        ! True => y is the active axis. False => y2 is active.
end type          

type qp_point_struct     ! A point on the page.
  real(rp) :: x = 0, y = 0
  character(16) :: units = ' '
end type

type qp_rect_struct     ! Rectangular: A structure with 4 numbers
  real(rp) :: x1 = 0, x2 = 0, y1 = 0, y2 = 0
  character(16) :: units = ' '
end type

type qp_text_struct
  real(rp) :: height         = 12   ! in points
  character(16) :: color     = 'black'
  logical :: uniform_spacing = .false.
end type

type qp_line_struct
  integer :: width         = 3
  character(16) :: color   = 'black'
  character(16) :: pattern = 'solid'
end type

type qp_symbol_struct
  character(16) :: type         = 'circle_dot'
  real(rp) :: height            = 10d0  ! in points (same as text height)
  character(16) :: color        = 'black'
  character(16) :: fill_pattern = 'solid_fill'
  integer :: line_width         = 3
end type

type qp_arrow_struct
  real(rp) :: head_angle = 30      ! Acute angle of the arrow point in degrees.
  real(rp) :: head_barb  = 0.4     ! Fraction of triangular arrow head that is cut away from the back.
  real(rp) :: head_size  = 1.0
  character(16) :: head_type   = 'filled_arrow_head'    ! Or 'outline_arrow_head'
  character(16) :: color       = 'black'
end type

type qp_legend_struct
  type (qp_point_struct) :: origin = qp_point_struct(0.2_rp, 0.0_rp, '%/GRAPH/LB')
  real(rp) :: row_spacing = 1.0              ! Spacing between rows.
  real(rp) :: line_length = 72               ! Length of the line in points.
  real(rp) :: text_offset = 10               ! Horizontal offset in points between the line and the text.
  logical :: draw_line = .true.              ! Draw lines? 
  logical :: draw_symbol = .true.            ! Draw symbols?
  logical :: draw_text = .true.              ! Draw text? 
end type

type qp_state_struct
  type (qp_plot_struct) plot
  type (qp_rect_struct) :: page   = qp_rect_struct (0.0, 0.0, 0.0, 0.0, ' ')
  type (qp_rect_struct) :: box    = qp_rect_struct (1.0, 2.0, 1.0, 2.0, ' ')
  type (qp_rect_struct) :: graph  = qp_rect_struct (1.0, 2.0, 1.0, 2.0, ' ')
  type (qp_rect_struct) :: margin = qp_rect_struct (0.0, 0.0, 0.0, 0.0, ' ')
  type (qp_rect_struct) :: border = qp_rect_struct (0.0, 0.0, 0.0, 0.0, ' ')
  type (qp_text_struct) :: main_title = qp_text_struct(18.0, 'black', .false.)
  type (qp_text_struct) :: graph_title= qp_text_struct(20.0, 'black', .false.)
  type (qp_text_struct) :: legend     = qp_text_struct(13.0, 'black', .false.)
  type (qp_text_struct) :: text       = qp_text_struct(18.0, 'black', .false.)
  type (qp_text_struct) :: axis_number= qp_text_struct(12.0, 'black', .false.)
  type (qp_text_struct) :: axis_label = qp_text_struct(15.0, 'black', .false.)
  type (qp_text_struct) :: this_text  ! current settings.
  type (qp_symbol_struct) :: symbol = qp_symbol_struct()
  type (qp_arrow_struct) :: arrow = qp_arrow_struct()
  type (qp_line_struct) :: std_line  = qp_line_struct()
  type (qp_line_struct) :: plot_line = qp_line_struct()
  type (qp_line_struct) :: axis_line = qp_line_struct()
  type (qp_line_struct) :: legend_line = qp_line_struct()
  type (qp_line_struct) :: grid_line = qp_line_struct(1, 'light_grey', 'solid')
  real(rp) :: text_scale = 1
  real(rp) :: text_spacing_factor = 0.6
  real(rp) :: dflt_axis_slop_factor = 1d-3
  character(16) :: text_background = 'not_set'
  integer :: max_axis_zero_digits = 3
  integer :: dflt_units = dflt_draw$
  integer :: max_digits = 8
  character(200) plot_file
  character(16) page_type       ! 'PS', 'X', etc.
  character(8) :: dflt_draw_units(3) = ['DATA ', 'GRAPH', 'LB   ' ]
  character(8) :: dflt_set_units(3)  = ['INCH ', 'PAGE ', 'LB   ' ]
  logical :: subgraph_on = .false.
  logical :: clip = .false.
  logical :: buffer = .false.   ! to be used by qp_save_state only
  logical :: uniform_symbol_size = .true.
end type

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function qp_string_to_enum (enum_str, enum_type, ix_dflt) result (ix_enum)
! 
! Routine to convert a string to the corresponding integer for enumerated parameter.
!
! Input:
!   enum_str      -- character(*): String representation (case insensitive) for the enumerated parameter. 
!   enum_type     -- character(*): Type of enum. Possibilities are: (case sensitive)
!                       "color", "line_pattern", "fill_pattern", "symbol_type", "arrow_head_type"
!   ix_dflt       -- integer, optional: Default number to set ix_enum to if enum_str is not valid.
!                     If not present, the default appropriate to enum_type is used.
!   error         -- logical, optional: Set True if there is a problem. False otherwise.
!
! Output:
!   ix_enum       -- integer: Index corresponding to enum_str.
!-

function qp_string_to_enum (enum_str, enum_type, ix_dflt, error) result (ix_enum)

implicit none

real(rp) cnum
integer ix_enum, ix, ixd, lb
integer, optional :: ix_dflt
character(*) enum_str, enum_type
logical, optional :: error

!

if (present(error)) error = .true.

select case (enum_type)
case ('color')
  if (enum_str(1:1) == 'Z') then    ! Continuous color mapping.
    if (.not. is_real(enum_str(2:), real_num = cnum)) then
      ix_enum = integer_option(-1 + lbound(qp_color_name, 1), ix_dflt)
    else
      ix_enum  = nint(17 + cnum * (huge(ix) - 17))
      if (present(error)) error = .false.
    endif
    return
  else
    call match_word (enum_str, qp_color_name, ix, .false., .false.)
    ixd = integer_option(black$, ix_dflt)
    lb = lbound(qp_color_name, 1)
  endif
case ('fill_pattern')
  call match_word (enum_str, qp_symbol_fill_pattern_name, ix, .false., .false.)
  ixd = integer_option(solid_fill$, ix_dflt)
  lb = lbound(qp_symbol_fill_pattern_name, 1)
case ('line_pattern')
  call match_word (enum_str, qp_line_pattern_name, ix, .false., .false.)
  ixd = integer_option(solid$, ix_dflt)
  lb = lbound(qp_line_pattern_name, 1)
case ('symbol_type')
  call match_word (enum_str, qp_symbol_type_name, ix, .false., .false.)
  ixd = integer_option(circle_dot_sym$, ix_dflt)
  lb = lbound(qp_symbol_type_name, 1)
case ('arrow_head_type')
  call match_word (enum_str, qp_arrow_head_type_name, ix, .false., .false.)
  ixd = integer_option(filled_arrow_head$, ix_dflt)
  lb = lbound(qp_arrow_head_type_name, 1)
case default
  ix_enum = integer_option(-1, ix_dflt)
  return
end select

if (ix == 0) then
  ix_enum = ixd
else
  ix_enum = ix + (lb - 1)
  if (present(error)) error = .false.
endif

end function qp_string_to_enum

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function qp_enum_to_string (ix_enum, enum_type, str_dflt, error) result (enum_str)
! 
! Routine to convert a string to the corresponding integer for enumerated parameter.
!
! Input:
!   ix_enum       -- integer: Index of enumerated number.
!   enum_type     -- character(*): Type of enum. Possibilities are:
!                       "color", "line_pattern", "fill_pattern", "symbol_type", "arrow_head_type"
!   str_dflt      -- character(*), optional: Default string to set enum_str to if ix_enum is not valid.
!                     If not present then the default appropriate to enum_type is used.
!
! Output:
!   enum_str      -- character(16): String representation for the enumerated parameter.
!   error         -- logical, optional: Set True if there is a problem. False otherwise.
!-

function qp_enum_to_string (ix_enum, enum_type, str_dflt, error) result (enum_str)

implicit none

integer ix_enum, ix
character(*) enum_type
character(*), optional :: str_dflt
character(16) enum_str
logical, optional :: error

!

if (present(error)) error = .true.

select case (enum_type)
case ('color')
  if (ix_enum > 16) then
    enum_str(1:1) = 'Z'
    write (enum_str(2:), '(f15.13)') (ix_enum - 17) / (1.0_rp * (huge(ix) - 17))
    if (present(error)) error = .false.
  elseif (ix_enum < lbound(qp_color_name, 1) .or. ix_enum > ubound(qp_color_name, 1)) then
    enum_str = string_option('black', str_dflt)
  else
    enum_str = downcase(qp_color_name(ix_enum))
    if (present(error)) error = .false.
  endif

case ('fill_pattern')
  if (ix_enum < lbound(qp_symbol_fill_pattern_name, 1) .or. ix_enum > ubound(qp_symbol_fill_pattern_name, 1)) then
    enum_str = string_option('solid', str_dflt)
  else
    enum_str = downcase(qp_symbol_fill_pattern_name(ix_enum))
    if (present(error)) error = .false.
  endif

case ('line_pattern')
  if (ix_enum < lbound(qp_line_pattern_name, 1) .or. ix_enum > ubound(qp_line_pattern_name, 1)) then
    enum_str = string_option('solid', str_dflt)
  else
    enum_str = downcase(qp_line_pattern_name(ix_enum))
    if (present(error)) error = .false.
  endif

case ('symbol_type')
  if (ix_enum < lbound(qp_symbol_type_name, 1) .or. ix_enum > ubound(qp_symbol_type_name, 1)) then
    enum_str = string_option('circle_dot', str_dflt)
  else
    enum_str = downcase(qp_symbol_type_name(ix_enum))
    if (present(error)) error = .false.
  endif

case ('arrow_head_type')
  if (ix_enum < lbound(qp_arrow_head_type_name, 1) .or. ix_enum > ubound(qp_arrow_head_type_name, 1)) then
    enum_str = string_option('filled', str_dflt)
  else
    enum_str = downcase(qp_arrow_head_type_name(ix_enum))
    if (present(error)) error = .false.
  endif

case default
  enum_str = string_option ('????', str_dflt)
end select

end function qp_enum_to_string

end module
