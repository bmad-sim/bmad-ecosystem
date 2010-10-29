module quick_plot_struct

use precision_def

integer, parameter :: white$ = 0, black$ = 1, red$ = 2, green$ = 3
integer, parameter :: blue$ = 4, cyan$ = 5, magenta$ = 6, yellow$ = 7
integer, parameter :: orange$ = 8, yellow_green$ = 9, light_green$ = 10
integer, parameter :: navy_blue$ = 11, purple$ = 12, redish_purple$ = 13
integer, parameter :: dark_grey$ = 14, light_grey$ = 15, transparent$ = 16

character(16), parameter :: qp_color_name(0:16) =   ['White        ', &
  'Black        ', 'Red          ', 'Green        ', 'Blue         ', &
  'Cyan         ', 'Magenta      ', 'Yellow       ', 'Orange       ', &
  'Yellow_Green ', 'Light_Green  ', 'Navy_Blue    ', 'Purple       ', &
  'Redish_Purple', 'Dark_Grey    ', 'Light_Grey   ', 'Transparent  ']

integer, parameter :: solid$ = 1, dashed$ = 2, dash_dot$ = 3
integer, parameter :: dotted$ = 4, dash_dot3$ = 5

character(16) :: qp_line_pattern_name(5) = (/ 'solid    ', &
    'dashed   ', 'dash_dot ', 'dotted   ', 'dash_dot3' /)

integer, parameter :: solid_fill$ = 1, no_fill$ = 2
integer, parameter :: hatched$ = 3, cross_hatched$ = 4

character(16) :: qp_fill_name(4) = (/ 'solid_fill   ', 'no_fill      ', &
                                      'hatched      ', 'cross_hatched' /)

integer, parameter :: square_sym$ = 0, dot_sym$ = 1, plus_sym$ = 2, times_sym$ = 3
integer, parameter :: circle_sym$ = 4, x_symbol_sym$ = 5, triangle_sym$ = 7
integer, parameter :: circle_plus_sym$ = 8, circle_dot_sym$ = 9
integer, parameter :: square_concave_sym$ = 10, diamond_sym$ = 11
integer, parameter :: star5_sym$ = 12, triangle_filled_sym$ = 13, red_cross_sym$ = 14
integer, parameter :: star_of_david_sym$ = 15, square_filled_sym$ = 16
integer, parameter :: circle_filled_sym$ = 17, star5_filled_sym$ = 18

character(16) :: qp_symbol_type_name(0:18) = (/ 'square         ', &
    'dot            ', 'plus           ', 'times          ', 'circle         ', &
    'x_symbol       ', '---------------', 'triangle       ', 'circle_plus    ', &
    'circle_dot     ', 'square_concave ', 'diamond        ', 'star5          ', &
    'triangle_filled', 'red_cross      ', 'star_of_david  ', 'square_filled  ', &
    'circle_filled  ', 'star5_filled   ' /)

integer, parameter :: dflt_draw$ = 1, dflt_set$ = 2

integer, parameter :: print_page_long_len = 10.5
integer, parameter :: print_page_short_len = 7.8

!------------------------------------

type qp_axis_struct
  character(80) :: label = ' '
  real(rp) :: min = 0, max = 10         ! min is actually left or bottom axis number.
  real(rp) :: number_offset = 0.05      ! offset from axis line in inches.
  real(rp) :: label_offset = 0.05       ! offset from numbers in inches.
  real(rp) :: major_tick_len = 0.10     ! in inches.
  real(rp) :: minor_tick_len = 0.06     ! in inches.
  integer :: label_color = black$       
  integer :: major_div = 5
  integer :: major_div_nominal = 5      ! Nominal value.
  integer :: minor_div = 0              ! 0 = auto choose.
  integer :: minor_div_max = 5          ! max number for auto choose.
  integer :: places = 0
  character(16) :: type = 'LINEAR'      ! or 'LOG', or 'CUSTOM'
  character(16) :: bounds = 'GENERAL'   ! or 'ZERO_AT_END' or 'ZERO_SYMMETRIC'
  integer :: tick_side = +1    ! +1 = draw to the inside, 0 = both, -1 = outside.
  integer :: number_side = -1  ! +1 = draw to the inside, -1 = outside.
  logical :: draw_label = .true.
  logical :: draw_numbers  = .true.
end type

! x2_mirrors_x and y2_mirrors_y force the x2 or y2 axis to have the same tick markings.
! What can be still different is whether a label and/or numbering is drawn.

type qp_plot_struct
  character(80) :: title = ' '
  type (qp_axis_struct) x, y, x2, y2
  type (qp_axis_struct), pointer :: xx, yy  ! Pointer to axes used for data plotting, etc. 
  logical :: draw_box    = .true.
  logical :: draw_title  = .true.
  logical :: draw_grid   = .true.
  logical limit_plot
  logical :: x2_mirrors_x = .true.
  logical :: y2_mirrors_y = .true.
  logical :: xx_points_to_x
  logical :: yy_points_to_y
end type          

type qp_point_struct     ! A point on the page.
  real(rp) x, y
  character(16) :: units = ' '
end type

type qp_rect_struct     ! Rectangular: A structure with 4 numbers
  real(rp) x1, x2, y1, y2   
  character(16) :: units = ' '
end type

type qp_text_struct
  real(rp) height   ! in points
  integer :: color = black$
  logical :: uniform_spacing = .false.
end type

type qp_line_struct
  integer :: width = 1
  integer :: color = black$
  integer :: pattern = solid$
end type

type qp_symbol_struct
  integer :: type = circle_dot_sym$
  real(rp) :: height      = 10d0  ! in points (same as text height)
  integer :: color        = black$
  integer :: fill_pattern = solid_fill$
  integer :: line_width   = 1
end type

type qp_state_struct
  type (qp_plot_struct) plot
  type (qp_rect_struct) :: page   = qp_rect_struct (0.0, 0.0, 0.0, 0.0, ' ')
  type (qp_rect_struct) :: box    = qp_rect_struct (1.0, 2.0, 1.0, 2.0, ' ')
  type (qp_rect_struct) :: graph  = qp_rect_struct (1.0, 2.0, 1.0, 2.0, ' ')
  type (qp_rect_struct) :: margin = qp_rect_struct (0.0, 0.0, 0.0, 0.0, ' ')
  type (qp_rect_struct) :: border = qp_rect_struct (0.0, 0.0, 0.0, 0.0, ' ')
  type (qp_text_struct) :: main_title = qp_text_struct(18.0, black$, .false.)
  type (qp_text_struct) :: graph_title= qp_text_struct(20.0, black$, .false.)
  type (qp_text_struct) :: legend     = qp_text_struct(13.0, black$, .false.)
  type (qp_text_struct) :: text       = qp_text_struct(18.0, black$, .false.)
  type (qp_text_struct) :: axis_number= qp_text_struct(10.0, black$, .false.)
  type (qp_text_struct) :: axis_label = qp_text_struct(15.0, black$, .false.)
  type (qp_text_struct) :: this_text  ! current settings.
  type (qp_symbol_struct) :: symbol 
  type (qp_line_struct) :: std_line  = qp_line_struct (2, black$, solid$)
  type (qp_line_struct) :: plot_line = qp_line_struct (2, black$, solid$)
  type (qp_line_struct) :: axis_line = qp_line_struct (2, black$, solid$)
  type (qp_line_struct) :: legend_line = qp_line_struct (2, black$, solid$)
  type (qp_line_struct) :: grid_line = qp_line_struct(1, light_grey$, solid$)
  real(rp) :: text_scale = 1
  real(rp) :: text_spacing_factor = 0.6
  real(rp) :: dflt_axis_slop_factor = 1e-3
  integer :: text_background = -1
  integer :: max_axis_zero_digits = 3
  integer :: dflt_units = dflt_draw$
  integer :: max_digits = 8
  character(80) plot_file
  character(16) page_type       ! 'PS', 'X', etc.
  character(8) :: dflt_draw_units(3) = (/ 'DATA ', 'GRAPH', 'LB   ' /)
  character(8) :: dflt_set_units(3)  = (/ 'INCH ', 'PAGE ', 'LB   ' /)
  logical :: subgraph_on = .false.
  logical :: clip = .false.
  logical :: buffer = .false.   ! to be used by qp_save_state only
  logical :: uniform_symbol_size = .true.
end type

end module
