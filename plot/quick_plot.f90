!+
! Module quick_plot
!
! Module that defines the QUICK_PLOT graphics plotting routines. 
!
! To fix:
!   * Generalize sub/super scripts and greek letters.
!   * Sub/super scripts and greek letters with uniform spacing.
!
!--------------------------------------------------------------------
!
! QUICK_PLOT uses the following concepts:
!   PAGE  -- The entire drawing surface. 
!   BOX   -- The area that the graph with axes, titles, etc. is placed into.
!   GRAPH -- The actual plotting area within the bounds of the axes.
!
!--------------------------------------------------------------------
!
! To set the defaults used for optional arguments use the routines:
!   qp_set_line_attrib
!   qp_set_line
!   gp set_symbol_attrib
!   gp set_symbol
!   qp_set_text_attrib
!
!--------------------------------------------------------------------
!
! To get the defaults used for optional arguments use the routines:
!   qp_get_line
!   gp get_symbol
!
!--------------------------------------------------------------------
!
! A "justify" argument is a character string with 2 characters.
! The first character gives the horizontal justification:
!   'L' -- Left justify
!   'C' -- Center justify
!   'R' -- Right justify
! The second character gives the vertical justification
!   'B' -- Bottom justify
!   'C' -- Center justify
!   'T' -- Top justify
!
!--------------------------------------------------------------------
!
! The "units" argument in GP routines is a character string which is divided
! up into three parts. The syntax of the string is:
!     'unit_type/ref_object/corner'
! Where
!   unit_type -- Type of units:
!           '%'       - Percent.
!           'DATA'    - Data units. (Draw default)
!           'MM'      - millimeters.
!           'INCH'    - Inches. (Set default)
!           'POINTS'  - Printers points (72 points = 1 inch, 1pt ~ 1pixel).
!
!   ref_object -- Reference object (optional except if unit_type = "%").
!           'PAGE'  -- Relative to the page (Set default).
!           'BOX'   -- Relative to the box.
!           'GRAPH' -- Relative to the graph (Draw default).
!
!   corner    -- Origin location (optional).
!           'LB' -- Left Bottom (Set and Draw default).
!           'LT' -- Left Top.
!           'RB' -- Right Bottom.
!           'RT' -- Right Top.
!
! Notes:
!  1) The DATA unit type, by definition, always uses the 
!     lower left corner of the GRAPH as a reference point.
!  2) For the '%' unit_type the '/' between unit_type and ref_object 
!     can be omitted.
!  3) If the corner is specified then the ref_object must appear also.
!  4) Everything must be in upper case.
!  5) For some routines (qp_set_margin, etc.) only a relative distance is 
!     needed. In this case the ref_object/corner part, if present, is ignored.
!  6) There are two defaults: One for drawing and another for setting 
!     margins, etc. Initially the draw default is 'DATA/GRAPH/LB' and the
!     set default is 'INCH/PAGE/LB'. use qp_set_parameters to change this.
!    
!
! Examples:
!     'DATA'          -- This is the draw default. 
!     'DATA/GRAPH/LB' -- Same as above.
!     'DATA/BOX/RT'   -- ILLEGAL: DATA must always go with GRAPH/LB.
!     '%PAGE/LT'      -- Percentage of page so (0.0, 1.0) = RT of page.
!     '%BOX'          -- Percentage of box so (1.0, 1.0) = RT of box.
!     'INCH/PAGE'     -- Inches from LB of page.
!
!--------------------------------------------------------------------
!
! The line styles:
!     1 - solid$                  Solid
!     2 - dashed$                 Dashed
!     3 - dash_dot$               Dash dot 
!     4 - dotted$                 Dotted
!     5 - dash_dot3$              Dash dot dot dot        
!
!--------------------------------------------------------------------
!
! The colors are:
!     0 - White$   (actually the background color)
!     1 - Black$   (actually the foreground color)
!     2 - Red$
!     3 - Green$
!     4 - Blue$
!     5 - Cyan$
!     6 - Magenta$
!     7 - Yellow$ 
!     8 - Orange$
!     9 - Yellow_Green$
!    10 - Light_Green$
!    11 - Navy_Blue$
!    12 - Purple$
!    13 - Redish_Purple$
!    14 - Dark_Grey$
!    15 - Light_Grey$
! 
!--------------------------------------------------------------------
!
! The fill styles are:
!     1 - solid_fill$        
!     2 - no_fill$           
!     3 - hatched$           
!     4 - cross_hatched$     
!
!--------------------------------------------------------------------
!
! The symbols are:
!     0 - square$
!     1 - dot$
!     2 - plus$
!     3 - times$
!     4 - circle$
!     5 - x_symbol$
!     7 - triangle$
!     8 - circle_plus$
!     9 - circle_dot$
!    10 - square_concave$
!    11 - diamond$
!    12 - star5$
!    13 - triangle_filled$
!    14 - red_cross$
!    15 - star_of_david$
!    16 - square_filled$
!    17 - circle_filled$
!    18 - star5_filled$
!
!--------------------------------------------------------------------
!
! The following are the symbol types:
!   -3 ... -31 - a regular polygon with abs(type) edges.
!           -2 - Same as -1.
!           -1 - Dot with diameter = current line width.
!    0 ...  31 - Standard marker symbols.
!   32 ... 127 - ASCII characters (in current font).
!                   E.G. to use letter F as a marker, set type = ICHAR('F'). 
!        > 127 - A Hershey symbol number.
!
!--------------------------------------------------------------------
!
! The text background index is:
!          -1 - Transparent background.
!           0 - Erase underlying graphics before drawing text.
!    1 to 255 - Opaque with the number specifying the color index.
!-

module quick_plot

use physical_constants
use sim_utils_interface
use utilities_mod
use output_mod
use plplot_interface

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
  integer :: style = solid$
end type

type qp_symbol_struct
  integer :: type = circle_dot$
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

!---------------------------------------------------------------------------
! common block
! General NOTE: qp_com is made private so that you cannot change it directly.
! This was done since the layout of qp_com can change.

type (qp_state_struct), pointer, save :: qp_com
type (qp_state_struct), target, save :: qp_save_com(0:20)

integer, save :: ix_qp_com = 0

private qp_save_com, ix_qp_com, qp_com

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_layout (x_axis, y_axis, x2_axis, y2_axis, &
!                       x2_mirrors_x, y2_mirrors_y, box, margin, page_border)
!
! Subroutine to set varies attributes. This routine can be used
! in place of other qp_set_* routines.
!
! Input:
!   x_axis       -- Qp_axis_struct, optional: Bottom axis.
!   y_axis       -- Qp_axis_struct, optional: Left axis.
!   x2_axis      -- Qp_axis_struct, optional: Right axis.
!   y2_axis      -- Qp_axis_struct, optional: Top axis.
!   x2_mirrors_x -- Logical, optional: If True then the x2 axis will mirror the 
!                       x axis. Mirroring means that the major and minor ticks of 
!                       the x2 axis will be the same as the x axis. 
!   y2_mirrors_y -- Logical, optional: If True then the y2 axis will mirror the 
!                       y axis. Mirroring means that the major and minor ticks of 
!                       the y2 axis will be the same as the y axis. 
!   box(4)       -- Integer, optional: Box to plot in. See qp_set_box.
!   margin       -- Qp_rect_struct, optional: Margin around the graph.
!   page_border  -- Qp_rect_struct, optional: Space around the page edge.
!-

subroutine qp_set_layout (x_axis, y_axis, x2_axis, y2_axis, &
                      x2_mirrors_x, y2_mirrors_y, box, margin, page_border)

implicit none

type (qp_axis_struct), optional :: x_axis, y_axis, x2_axis, y2_axis
type (qp_rect_struct), optional :: margin, page_border

integer, optional :: box(4)
logical, optional :: x2_mirrors_x, y2_mirrors_y

! axes set

if (present(x_axis))  qp_com%plot%x = x_axis
if (present(y_axis))  qp_com%plot%y = y_axis
if (present(x2_axis)) qp_com%plot%x2 = x2_axis
if (present(y2_axis)) qp_com%plot%y2 = y2_axis
if (present(x2_mirrors_x)) qp_com%plot%x2_mirrors_x = x2_mirrors_x
if (present(y2_mirrors_y)) qp_com%plot%y2_mirrors_y = y2_mirrors_y

! set world coords if needed

if (present(box)) call qp_set_box (box(1), box(2), box(3), box(4))
if (present(margin)) call qp_set_margin (margin%x1, margin%x2, &
                                    margin%y1, margin%y2, margin%units)
if (present(page_border)) call qp_set_page_border (page_border%x1, &
           page_border%x2, page_border%y1, page_border%y2, page_border%units)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_get_layout_attrib (who, x1, x2, y1, y2, units)
!
! Subroutine to get the attributes of the layout.
!
! Input:
!   who -- Character(*): 
!                 "PAGE"   = Page size. In this case x1 and y1 are 0.
!                 "BOX"    = Placement of the box on the page.
!                 "GRAPH"  = Placement of the Graph on the page.
!                 "MARGIN" = Distances between the edges of the graph and 
!                              edges of the box.
!                 "BORDER" = Distances between the edges of the page and the
!                              area where boxes are placed.
!  units -- Character(*), optional: Units of returned numbers.
!                  'MM'     - milimeters.
!                  'INCH'   - Inches (default unless default_set_units is changed).
!                  'POINTS' - Point.
!
! Output:
!   x1  -- Real(rp): Left distance.
!   x2  -- Real(rp): Right distance or page width.
!   y1  -- Real(rp): Bottom distance.
!   y2  -- Real(rp): Top distance or page height.
!-

subroutine qp_get_layout_attrib (who, x1, x2, y1, y2, units)

implicit none

type (qp_rect_struct) rect

real(rp) x1, x2, y1, y2

character(*) who
character(*), optional :: units
character(40) :: r_name = "qp_get_layout_attrib"
!

if (who == 'PAGE') then
  rect = qp_com%page
elseif (who == 'BOX') then
  rect = qp_com%box
elseif (who == 'GRAPH') then
  rect = qp_com%graph
elseif (who == 'MARGIN') then
  rect = qp_com%margin
elseif (who == 'BORDER') then
  rect = qp_com%border
else
  call out_io (s_fatal$, r_name, 'BAD "WHO": ' // who)
  call err_exit
endif

qp_com%dflt_units = dflt_set$

call qp_from_inch_abs (rect%x1, rect%y1, x1, y1, units)
call qp_from_inch_abs (rect%x2, rect%y2, x2, y2, units)

qp_com%dflt_units = dflt_draw$

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_init_com_struct 
! 
! Subroutine to initialize the common block qp_state_struct.
! This subroutine is not for general use.
!-
             
subroutine qp_init_com_struct

implicit none

!

qp_com => qp_save_com(0)

qp_com%plot%x2%draw_numbers = .false.
qp_com%plot%y2 = qp_com%plot%x2
qp_com%plot%xx_points_to_x = .true.
qp_com%plot%yy_points_to_y = .true.

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_save_state (buffer_basic)
!
! Subroutine to save the current attributes. 
! Use qp_restore_state to restore the saved state.
!
! There are two buffers: One for quick_plot and the other for
! the underlying plot library (PGPLOT at present). The underlying buffer
! is only used to buffer the graphical output (nothing will be seen on
! the plot window until the buffering is ended) and does not buffer
! the underlying state.
!
! qp_save_state always buffers the quick_plot state and the buffer_basic
!   argument determines whether the base state is buffered as well.
! The base state is the state of the underlying plot library. 
!
! Note: Buffering can make things go faster but with buffering the display 
!   will not be updated until you call qp_restore_state.
!
! Input:
!   buffer_basic -- Logical: If True then buffer the base state.
!-

subroutine qp_save_state (buffer_basic)

implicit none

logical buffer_basic

character(16) :: r_name = 'qp_save_state'

!

if (buffer_basic) call qp_save_state_basic

if (ix_qp_com == size(qp_save_com)) then
  call out_io (s_fatal$, r_name, 'TRYING TO SAVE TOO MANY STATES!')
  call err_exit
endif

if (ix_qp_com == 0) call qp_init_com_struct
ix_qp_com = ix_qp_com + 1

qp_save_com(ix_qp_com) = qp_com
qp_com => qp_save_com(ix_qp_com)
qp_com%buffer = buffer_basic

if (qp_com%plot%xx_points_to_x) then
  qp_com%plot%xx => qp_com%plot%x
else
  qp_com%plot%xx => qp_com%plot%x2
endif

if (qp_com%plot%yy_points_to_y) then
  qp_com%plot%yy => qp_com%plot%y
else
  qp_com%plot%yy => qp_com%plot%y2
endif

end subroutine  

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_restore_state
!
! Subroutine to restore saved attributes. 
! Use qp_save_state to restore the saved state.
!-

subroutine qp_restore_state

implicit none

character(16) :: r_name = 'qp_restore_state'

!

if (qp_com%buffer) call qp_restore_state_basic 

if (ix_qp_com == 0) then
  call out_io (s_fatal$, r_name, 'NO STATE TO RESTORE!')
  call err_exit
endif

ix_qp_com = ix_qp_com - 1
qp_com => qp_save_com(ix_qp_com)

end subroutine  

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_pointer_to_axis (axis_str, axis_ptr)
!   
! Subroutine to return a pointer to an common block axis.
!
! Input:
!   axis_str -- Character(*): 
!                 'X' to set the Left x-axis
!                 'Y' to set the Bottom y-axis.
!                 'X2' to set the Right x-axis
!                 'Y2' to set the Top y-axis.
!
! Output:
!   axis_ptr -- Qp_axis_struct, pointer: Pointer to the common block axis.
!-

subroutine qp_pointer_to_axis (axis_str, axis_ptr)

implicit none

type (qp_axis_struct), pointer :: axis_ptr

character(*) axis_str
character(20) :: r_name = 'qp_pointer_to_axis'

!

if (axis_str == 'X') then
  axis_ptr => qp_com%plot%x
elseif (axis_str == 'Y') then
  axis_ptr => qp_com%plot%y
elseif (axis_str == 'X2') then
  axis_ptr => qp_com%plot%x2
elseif (axis_str == 'Y2') then
  axis_ptr => qp_com%plot%y2
else
  call out_io (s_fatal$, r_name, 'INVALID AXIS: ' // axis_str)
  call err_exit
endif

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_use_axis (x, y)
!
! Subroutine to set what axis to use: X or X2, Y or Y2.
!
! Input:
!   x -- Character(*), optional: 'X', or 'X2'
!   y -- Character(*), optional: 'Y', or 'Y2'
!-

subroutine qp_use_axis (x, y)

implicit none

character(*), optional :: x, y
character(10) :: r_name = 'qp_use_axis'

!

if (present(x)) then
  select case (x)
  case ('X')
    qp_com%plot%xx => qp_com%plot%x
    qp_com%plot%xx_points_to_x = .true.
  case ('X2')
    qp_com%plot%xx => qp_com%plot%x2
    qp_com%plot%xx_points_to_x = .false.
  case default
    call out_io (s_error$, r_name, 'BAD "X": ' // x)
  end select
endif

if (present(y)) then
  select case (y)
  case ('Y')
    qp_com%plot%yy => qp_com%plot%y
    qp_com%plot%yy_points_to_y = .true.
  case ('Y2')
    qp_com%plot%yy => qp_com%plot%y2
    qp_com%plot%yy_points_to_y = .false.
  case default
    call out_io (s_error$, r_name, 'BAD "Y": ' // y)
  end select
endif

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_axis (axis_str, a_min, a_max, div, places, label, draw_label, 
!               draw_numbers, minor_div, minor_div_max, mirror,
!               number_offset, label_offset, major_tick_len, minor_tick_len, ax_type)
!   
! Subroutine to set (but not plot) the min, max, divisions etc. for the
! X and Y axes. 
!
! Input:
!   axis_str  -- Character(*): 
!                 'X' to set the Left x-axis
!                 'Y' to set the Bottom y-axis.
!                 'X2' to set the Right x-axis
!                 'Y2' to set the Top y-axis.
!   a_min     -- Real(rp), optional: Axis minimum.
!   a_max     -- Real(rp), optional: Axis maximum.
!   div       -- Integer, optional: Number of major divisions.
!   places    -- Integer, optional: Number of decmal places after the decimal 
!                   point. A Negative number surpresses that number of zeros.
!                   E.g. For the number 30 then 
!                        places =  2 -> Output is: "30.00"
!                        places =  0 -> Output is: "30"
!                        places = -1 -> Output can be scalled to: "3"
!   label         -- Character(*), optional: Axis label.
!   draw_label    -- Logical, optional: Draw axis label.
!   draw_numbers  -- Logical, optional: Draw axis numbers
!   minor_div     -- Integer, optional: Number of minor divisions.
!   minor_div_max -- Integer, optional: Maximum number of minor divisions.
!                     This is used when you want Quick_Plot to pick the
!                     actual number of minor divisions
!   mirror        -- Logical, optional: If True and axis = "x" or axis = "x2" then
!                     the x2 axis will mirror the x axis. Mirroring means that
!                     the major and minor ticks of the x2 axis will be the 
!                     same as the x axis. A similar situation holds if axis = "y"
!                     or axis = "y2".
!   number_offset  -- Real(rp), optional: Offset from axis line in inches.
!   label_offset   -- Real(rp), optional: Offset form numbers in inches.
!   major_tick_len -- Real(rp), optional: Major tick length in inches.
!   minor_tick_len -- Real(rp), optional: Minor tick length in inches.
!   ax_type           -- Character(16): Axis type. 'LINEAR', or 'LOG'.
!-

subroutine qp_set_axis (axis_str, a_min, a_max, div, places, label, draw_label, &
                  draw_numbers, minor_div, minor_div_max, mirror, &
                  number_offset, label_offset, major_tick_len, minor_tick_len, ax_type)

implicit none

type (qp_axis_struct), pointer :: this_axis
real(rp), optional :: a_min, a_max, number_offset
real(rp), optional :: label_offset, major_tick_len, minor_tick_len

integer, optional :: div, places, minor_div, minor_div_max
logical, optional :: draw_label, draw_numbers, mirror

character(*), optional :: label, ax_type
character(*) axis_str

!

call qp_pointer_to_axis (axis_str, this_axis)

if (present(a_min))  this_axis%min = a_min
if (present(a_max))  this_axis%max = a_max
if (present(div))    this_axis%major_div = div
if (present(places)) this_axis%places = places

if (present(label)) this_axis%label = label
if (present(draw_label)) this_axis%draw_label     = draw_label
if (present(draw_numbers)) this_axis%draw_numbers = draw_numbers
if (present(minor_div)) this_axis%minor_div       = minor_div
if (present(minor_div_max)) then
  this_axis%minor_div_max = minor_div_max
  this_axis%minor_div = 0
endif

if (axis_str == 'X2') then
  if (present(a_min) .or. present(a_max)) qp_com%plot%x2_mirrors_x = .false.
endif

if (axis_str == 'Y2') then
  if (present(a_min) .or. present(a_max)) qp_com%plot%y2_mirrors_y = .false.
endif

if (present(mirror)) then
  if (axis_str(1:1) == 'X') qp_com%plot%x2_mirrors_x = mirror
  if (axis_str(1:1) == 'Y') qp_com%plot%y2_mirrors_y = mirror
endif


if (present(number_offset))  this_axis%number_offset  = number_offset
if (present(label_offset))   this_axis%label_offset   = label_offset
if (present(major_tick_len)) this_axis%major_tick_len = major_tick_len
if (present(minor_tick_len)) this_axis%minor_tick_len = minor_tick_len
if (present(ax_type))        this_axis%type = ax_type

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_get_axis (axis_str, a_min, a_max, div, places, label, 
!               draw_label, draw_numbers, minor_div, mirror, number_offset, 
!               label_offset, major_tick_len, minor_tick_len, ax_type)
!   
! Subroutine to get the min, max, divisions etc. for the X and Y axes.  
!
! Input:
!   axis_str  -- Character(*): 
!                 'X' to set the Left x-axis
!                 'Y' to set the Bottom y-axis.
!                 'X2' to set the Right x-axis
!                 'Y2' to set the Top y-axis.
!   a_min     -- Real(rp), optional: Axis minimum.
!   a_max     -- Real(rp), optional: Axis maximum.
!   div       -- Integer, optional: Number of major divisions.
!   places    -- Integer, optional: Number of decmal places after the decimal 
!                   point. A Negative number surpresses that number of zeros.
!                   E.g. For the number 30 then 
!                        places =  2 -> Output is: "30.00"
!                        places =  0 -> Output is: "30"
!                        places = -1 -> Output can be scalled to: "3"
!   label          -- Character(*), optional: Axis label.
!   draw_label     -- Logical, optional: Draw axis label.
!   draw_numbers   -- Logical, optional: Draw axis numbers
!   minor_div      -- Integer, optional: Number of minor divisions.
!   mirror         -- Logical, optional: Mirror x2 or y2 axis?
!   number_offset  -- Real(rp), optional: Offset from axis line in inches.
!   label_offset   -- Real(rp), optional: Offset form numbers in inches.
!   major_tick_len -- Real(rp), optional: Major tick length in inches.
!   minor_tick_len -- Real(rp), optional: Minor tick length in inches.
!   ax_type        -- Character(16): Axis type. 'LINEAR', or 'LOG'.
!-

subroutine qp_get_axis (axis_str, a_min, a_max, div, places, label, draw_label, &
                draw_numbers, minor_div, mirror, number_offset, &
                label_offset, major_tick_len, minor_tick_len, ax_type)

implicit none

type (qp_axis_struct), pointer :: this_axis
real(rp), optional :: a_min, a_max, number_offset
real(rp), optional :: label_offset, major_tick_len, minor_tick_len

integer, optional :: div, places, minor_div
logical, optional :: draw_label, draw_numbers, mirror

character(*), optional :: label, ax_type
character(*) axis_str

!

call qp_pointer_to_axis (axis_str, this_axis)

if (present(a_min))  a_min  = this_axis%min  
if (present(a_max))  a_max  = this_axis%max  
if (present(div))    div    = this_axis%major_div  
if (present(places)) places = this_axis%places  

if (present(label))        label        = this_axis%label
if (present(draw_label))   draw_label   = this_axis%draw_label   
if (present(draw_numbers)) draw_numbers = this_axis%draw_numbers 
if (present(minor_div))    minor_div    = this_axis%minor_div  

if (present(mirror)) then
  if (axis_str(1:1) == 'X') mirror = qp_com%plot%x2_mirrors_x
  if (axis_str(1:1) == 'Y') mirror = qp_com%plot%y2_mirrors_y
endif

if (present(number_offset))  number_offset  = this_axis%number_offset   
if (present(label_offset))   label_offset   = this_axis%label_offset    
if (present(major_tick_len)) major_tick_len = this_axis%major_tick_len 
if (present(minor_tick_len)) minor_tick_len = this_axis%minor_tick_len  
if (present(ax_type))        ax_type        = this_axis%type

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_calc_and_set_axis (axis_str, data_min, data_max, 
!                         div_min, div_max, bounds, axis_type, slop_factor)
!
! Subroutine to calculate a "nice" plot scale given the minimum and maximum
! of the data. 
!
! Note: If data_min > data_max then on output axis%min > axis%max
!
! Input:
!   axis_str    -- Character(*): 
!                   'X' to set the Left x-axis
!                   'Y' to set the Bottom y-axis.
!                   'X2' to set the Right x-axis
!                   'Y2' to set the Top y-axis.
!   data_min    -- Real(rp): Minimum of the data
!   data_max    -- Real(rp): Maximum of the data
!   div_min     -- Integer: Minimum number of divisions.
!   div_max     -- Integer: Maximum number of divisions.
!   bounds      -- Character(*):
!                        'ZERO_AT_END'    -- Make AXIS_MIN or AXIS_MAX zero.
!                        'ZERO_SYMMETRIC' -- Make AXIS_MIN = -AXIS_MAX
!                        'GENERAL'        -- No restriction on min or max.
!   axis_type   -- Character(*), optional: Type of axis. 'LINEAR' or 'LOG'.
!                      Default is 'LINEAR'
!   slop_factor -- Real(rp): See qp_calc_axis_scale for info.
!
! Example:
!   call qp_calc_and_set_axis ('X', 352, 378, 4, 6, 'ZERO_AT_END')
!
! Gives for the x-axis:
!   min       = 0
!   max       = 400
!   places    = 0      ! places after the decimal point needed
!   divisions = 4
!-

subroutine qp_calc_and_set_axis (axis_str, data_min, data_max, &
                      div_min, div_max, bounds, axis_type, slop_factor)

implicit none

type (qp_axis_struct), pointer :: ax

integer divisions, places
integer div_min, div_max

real(rp) data_max, data_min, axis_max, axis_min
real(rp), optional :: slop_factor

character(*) axis_str, bounds
character(*), optional :: axis_type

!

call qp_pointer_to_axis (axis_str, ax)
ax%bounds = bounds
ax%type = 'LINEAR'
if (present(axis_type)) ax%type = axis_type
call qp_calc_axis_params (data_min, data_max, div_min, div_max, ax, slop_factor)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_calc_axis_params (data_min, data_max, div_min, div_max, axis, slop_factor)
!
! Subroutine to calculate a "nice" plot scale given the minimum and maximum
! of the data. This is similar to CALC_AXIS_SCALE but here the subroutine will
! pick the number of divisions.
!
! Note: If data_min > data_max then on output axis_min > axis_max
!
! Input:
!   data_min    -- Real(rp): Minimum of the data
!   data_max    -- Real(rp): Maximum of the data
!   div_min     -- Integer: Minimum number of divisions.
!   div_max     -- Integer: Maximum number of divisions.
!   axis        -- qp_axis_struct: Structure holding the axis parameters.
!     %axis_type  -- Character(*): Type of axis. 'LINEAR' or 'LOG'.
!     %bounds     -- Character(*): Only used for 'LINEAR' axis types.
!                       'ZERO_AT_END'    -- Make AXIS_MIN or AXIS_MAX zero.
!                       'ZERO_SYMMETRIC' -- Make AXIS_MIN = -AXIS_MAX
!                       'GENERAL'        -- No restriction on min or max.
!   slop_factor -- Real(rp): See qp_calc_axis_scale for info.
!
! Output:
!   axis       -- qp_axis_struct: Structure holding the axis parameters.
!     %places     -- Integer: Number of places after the decimal point needed
!                     to display the axis numbers.
!     %min        -- Real(rp): Axis minimum.
!     %max        -- Real(rp): Axis maximum.
!     %major_div  -- Integer: How many divisions the axis is divided up into
!
!
! Example:
!   axis%bounds = 'ZERO_AT_END'
!   call qp_calc_axis_params (352, 378, 4, 6, axis)
!
! Gives:
!   axis%min       = 0
!   axis%max       = 400
!   axis%places    = 0
!   axis%major_div = 4
!-

subroutine qp_calc_axis_params (data_min, data_max, div_min, div_max, axis, slop_factor)

implicit none

type (qp_axis_struct) axis

integer i, div_min, div_max, div_best
real(rp) data_max, data_min, d_max, d_min, score, score_max
real(rp), optional :: slop_factor

! 

score_max = -1e20
axis%major_div = 0
d_min = data_min
d_max = data_max

do i = div_min, div_max
  axis%major_div = i
  call qp_calc_axis_scale (d_min, d_max, axis, score, slop_factor)
  if (score_max < score) then
    score_max = score
    div_best = i
  endif
enddo

axis%major_div = div_best
call qp_calc_axis_scale (d_min, d_max, axis, score, slop_factor)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_calc_axis_places (axis)
!
! Subroutine to calculate the number of decmal places needed to display the
! axis numbers. Note: Reversed axes with axis%min > axis%max is OK.
!
! Note: axis%places is ignored for axis%type = 'LOG'.
!
! Input:
!   axis       -- qp_axis_struct: Structure holding the axis parameters.
!     %min        -- Real(rp): Axis bound.
!     %max        -- Real(rp): Axis bound.
!     %major_div  -- Integer: How many divisions the axis is divided up into
!
! Output:
!   axis       -- qp_axis_struct: Structure holding the axis parameters.
!     %places    -- Integer: Number of places after the decimal point needed
!                      to display the axis numbers. places can be negative
!                      if axis numbers are divisible by factors of 10. For example,
!                      axis numbers of 100, 200, 300 will have places = -2.
!-

subroutine qp_calc_axis_places (axis)
      
implicit none

type (qp_axis_struct) axis

integer i
real(rp) width, num, num2, a_min, a_max, max_d, effective_zero

! LOG scale does not use places

if (axis%type == 'LOG') return

! sort true min and max

a_min = min (axis%min, axis%max)
a_max = max (axis%min, axis%max)
max_d = 10.0_rp**qp_com%max_digits
effective_zero = max(abs(a_max), abs(a_min)) / max_d

! If min = max then cannot do a calculation so don't do anything.

if (a_min == a_max) return

! First calculation: Take each axis number and find how many digits it has.
! The number of places is the maximum number of digits needed to represent
! all the numbers on the axis with a limit of qp_com%max_digits digits maximum for any one number.

axis%places = -1000
do i = 0, axis%major_div
  num = a_min + i * (a_max - a_min) / axis%major_div
  if (abs(num) < effective_zero) cycle  ! Ignore zero
  axis%places = max(axis%places, floor(-log10(abs(num))))
  do
    num2 = num * 10d0**axis%places
    if (abs(num2 - nint(num2)) < 0.01 .or. abs(num2) > max_d) exit
    axis%places = axis%places + 1
  enddo
enddo

! Second calculation: Places based upon the width of the plot.
! The number of places returned by the subroutine is the maximum of the
! two calculations

width = abs(a_max - a_min) / axis%major_div
axis%places = max(axis%places, floor(-log10(width)+0.9))

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_calc_axis_scale (data_min, data_max, axis, niceness_score, slop_factor)
!
! Subroutine to calculate a "nice" plot scale given the minimum and maximum
! of the data. If data_max - data_min < 1e-30 then axis%max - axis%min will
! be at least axis%major_div * 1e-30.
!
! The slop_factor argument is used to increase the effective data_min and
! decrease the effective data_max. For example, if data_min = 0 and 
! data_max = 1.0000001 then an axis with min = 0 and max = 1 may be acceptable.
! The effective data min/max is computed via:
!    data_min (eff) = data_min + (data_max - data_min) * slop_factor
!    data_max (eff) = data_max - (data_max - data_min) * slop_factor
! For a log scale substitute log(data_min) for data_min and log(data_max) 
! for data_max in the above equations.
!
! Note: If you are not sure how many divisions the axis should have then
! consider the routine qp_calc_axis_params.
!
! Note: If data_min > data_max then on output axis%min > axis%max
!
! Input:
!   data_min    -- Real(rp): Minimum of the data
!   data_max    -- Real(rp): Maximum of the data
!   axis        -- qp_axis_struct: Structure holding the axis parameters.
!     %axis_type  -- Character(*): Type of axis. 'LINEAR' or 'LOG'.
!     %bounds     -- Character(*): Only used for 'LINEAR' axis types.
!                      'ZERO_AT_END'    -- Make AXIS%MIN or AXIS%MAX zero.
!                      'ZERO_SYMMETRIC' -- Make AXIS%MIN = -AXIS%MAX
!                      'GENERAL'        -- No restriction on min or max.
!     %major_div  -- Integer: How many divisions the axis is divided up into
!   slop_factor -- Real(rp), optional: If not present then the a default is used.
!                      This default is set by qp_set_parameters.
!
! Output:
!   axis       -- qp_axis_struct: Structure holding the axis parameters.
!     %places    -- Integer: Number of axis%places after the decimal point needed
!                     to display the axis numbers. 
!     %min       -- Real(rp): Axis minimum.
!     %max       -- Real(rp): Axis maximum.
!   niceness_score -- Real(rp), optional: Score as to how "nice" 
!                 axis%min and axis%max are. The larger the number the nicer.
!
! Example:
!
!   call qp_calc_axis_scale (3.52, 3.78, 4, 'GENERAL', axis)
! Gives:
!   axis%min = 3.4
!   axis%max = 3.8
!   axis%places = 1
!
!   call qp_calc_axis_scale (352, 378, 4, 'ZERO_AT_END', axis)
! Gives:
!   axis%min = 0
!   axis%max = 400
!   axis%places = 0
!-

subroutine qp_calc_axis_scale (data_min, data_max, axis, niceness_score, slop_factor)

implicit none
           
type (qp_axis_struct) axis

integer div_eff, m_min, m_max, m
integer min1, min2, max1, max2, imin, imax, j, ave

real(rp), optional :: niceness_score
real(rp) data_max, data_min, r, a_min, a_max, slop
real(dp) data_width, data_width10, min_width, max_score, score, log_width, aa
real(rp), optional :: slop_factor

character(20) :: r_name = 'qp_calc_axis_scale'

! Error check

if (axis%major_div < 1) then
  call out_io (s_abort$, r_name, &
                    '"AXIS%MAJOR_DIV" NUMBER IS LESS THAN 1! \i\ ', axis%major_div)
  call err_exit
endif

! 'LOG' axis

if (axis%type == 'LOG') then

  if (data_min <= 0 .or. data_max <= 0) then
    call out_io (s_abort$, r_name, &
                    'DATA IS NEGATIVE FOR LOG AXIS CALC: \i\ ', min(data_min, data_max))
    call err_exit
  endif

  a_min = log10(data_min)
  a_max = log10(data_max)
  slop = (a_max - a_min) * real_option(qp_com%dflt_axis_slop_factor, slop_factor)
  a_min = a_min + slop
  a_max = a_max - slop
  m_min = floor(a_min)
  m_max = ceiling (a_max)

  r = real(m_max - m_min) / axis%major_div

  if (axis%places >= m_max-m_min) then 
    if (present(niceness_score)) niceness_score = 0.0  ! perfect
    m = 0
  else
    if (present(niceness_score)) niceness_score = -(r - floor(r)) 
    m = ceiling(r) * axis%major_div  - (m_max - m_min)
  endif

  if (a_min/m_min < m_max/a_max) then
    axis%max = 10d0 ** (m_max + m/2)
    axis%min = 10d0 ** (m_min - (m+1)/2)
  else
    axis%max = 10d0 ** (m_max + (m+1)/2)
    axis%min = 10d0 ** (m_min - m/2)
  endif

  return

endif  ! log scale

!-----------------------------------------------------------------------
! Below is for a non-log axis
! Error check           

if (axis%bounds == 'ZERO_AT_END' .and. data_max*data_min < 0) then
  call out_io (s_error$, r_name, 'DATA ABOVE AND BELOW ZERO!', &
                                 'WITH "ZERO_AT_END"')
  niceness_score = -1
  return
endif

! find width of data

if (axis%bounds == 'ZERO_AT_END') then
  a_max = max(abs(data_max), abs(data_min))
  a_min = 0
elseif (axis%bounds == 'ZERO_SYMMETRIC') then
  a_max = max(abs(data_max), abs(data_min))
  a_min = 0
elseif (axis%bounds == 'GENERAL') then
  a_min = min (data_min, data_max)
  a_max = max (data_min, data_max)
else
  call out_io (s_fatal$, r_name, 'I DO NOT UNDERSTAND "AXIS%BOUNDS": ' // axis%bounds)
  call err_exit
endif

! add slop factor

data_width = a_max - a_min
slop = real_option(qp_com%dflt_axis_slop_factor, slop_factor)
a_min = a_min + data_width * slop
a_max = a_max - data_width * slop

! find possible candidates
            
min_width = axis%major_div * max(abs(a_max)*1e-5, abs(a_min)*1e-5, 1e-29_rp)
data_width = max(data_width, min_width)
log_width = log10(data_width)
data_width10 = 10d0**(floor(log_width)-1)
              
if (axis%bounds == 'ZERO_SYMMETRIC') then
  div_eff = axis%major_div / 2
else
  div_eff = axis%major_div
endif
min_width = data_width10 * axis%major_div

if (axis%bounds == 'ZERO_AT_END') then
  min1 = 0
  min2 = 0
else
  min2 = floor(a_min / data_width10)
  min1 = floor(a_min / (100 * data_width10)) * 100
  min1 = min(min1, min2 - 10*div_eff)
endif

aa = a_max/data_width10 + 1
max1 = floor(aa)
aa = a_max / (100 * data_width10) + 1
max2 = max(floor(aa) * 100, max1 + 10*div_eff)

ave = nint((a_max + a_min)/data_width10)
min1 = max(min1, 2*min2 - max1)
min1 = min(min1, (ave - axis%major_div)/2 - 1)

max2 = min(max2, 2*max1 - min2)
max2 = max(max2, min1+axis%major_div+1)

! go through and rate all the possibilities and choose the one
! with the highest score

max_score = -1000
do imin = min1, min2
  do imax = max1, max2
    score = qp_axis_niceness (imin, imax, div_eff)
    if (score > max_score) then
      max_score = score
      axis%min = imin * data_width10
      axis%max = imax * data_width10
    endif
  enddo
enddo

if (present(niceness_score)) niceness_score = max_score 

! adjust the scale if necessary

if (axis%bounds == 'ZERO_AT_END' .and. (data_min < 0 .or. data_max < 0)) then
  axis%min = -axis%max
  axis%max = 0
elseif (axis%bounds == 'ZERO_SYMMETRIC') then
  axis%min = -axis%max
endif

! find number of places needed

call qp_calc_axis_places (axis)

! reverse max/min if data max/min is reversed

if (data_min > data_max) then
  aa = axis%min 
  axis%min = axis%max
  axis%max = aa
endif

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function qp_axis_niceness (imin, imax, divisions) result (score)
!
! Routine to calculate how "nicely" an axis will look.
! The higher the score the nicer.
!
! Input:
!   imin      -- Integer: Axis minimum
!   imax      -- Integer: Axis maximum
!   divisions -- Integer: Number of divisions
! 
! Output:
!   score -- Real(rp): Niceness score.
!-

function qp_axis_niceness (imin, imax, divisions) result (score)

implicit none

integer imin, imax, divisions
integer del, i, j, im
real(rp) score

!

del = imax - imin

if (del == 0) then
  score = -1000
  return
endif

im = mod(del, divisions) 
if (im /= 0) then
  do i = 1, 5
    im = mod(10*im , divisions)
    if (im == 0) then
      score = -i
      return
    endif
  enddo
  score = -500
  return
endif

!

score = 0
if (imin == 0) score = score + 1
if (imax == 0) score = score + 1

do j = imin, imax, (imax - imin) / divisions
  if (mod(j, 5) == 0) score = score + 1.0 / (divisions + 1)
  if (mod(j, 10) == 0) score = score + 1.5 / (divisions + 1)
  if (mod(j, 100) == 0) score = score + 2.0 / (divisions + 1)
  if (mod(j, 2) == 0) score = score + 1.0 / (divisions + 1)
  if (j == 0) score = score + 1
enddo

score = score + 20 * (2 - log10(float(imax - imin)))

end function

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_calc_axis_divisions (axis_min, axis_max, 
!                                       div_min, div_max, divisions)
!
! Routine to calculate the best (gives the nicest looking drawing) number 
! of major divisions for fixed axis minimum and maximum.
! 
! Input:
!   axis_min  -- Real(rp): Axis minimum.
!   axis_max  -- Real(rp): Axis maximum.
!   div_min   -- Integer: Smallest number of divisions possible.
!   div_max   -- Integer: Greatest number of divisions possible.
!
! Output:
!   divisions -- Integer: Best number of divisions.
!                 This will be between div_min and div_max.
!-

subroutine qp_calc_axis_divisions (axis_min, axis_max, &
                                      div_min, div_max, divisions)

implicit none

real(rp) axis_min, axis_max
real(rp) data_width, min_width, log_width, data_width10

integer div_min, div_max, divisions
integer i, score, this_score, imin, imax

!

data_width = abs(axis_max - axis_min)
min_width = div_min * max(abs(axis_max)*1e-5, abs(axis_min)*1e-5, 1e-29_rp)
data_width = max(data_width, min_width)
log_width = log10(data_width)
data_width10 = 10d0**(floor(log_width)-1)
                
imin = nint(axis_min / data_width10)
imax = nint(axis_max / data_width10)

score = -10000
do i = div_min, div_max
  this_score = qp_axis_niceness(imin, imax, i)
  if (this_score > score) then
    score = this_score
    divisions = i
  endif
enddo

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_graph_limits
!
! Subroutine to calculate the offsets for the graph.
! Note: This subroutine is for internal use only.
!-


subroutine qp_set_graph_limits 

implicit none

type (qp_rect_struct), pointer :: graph

! qp_com%graph%z1 is the lower left corner of plot
! qp_com%graph%z2 is the upper right corner of plot

graph => qp_com%graph

if (.not. qp_com%subgraph_on) then      ! if no subgraph
  graph%x1 = qp_com%box%x1 + qp_com%margin%x1  
  graph%y1 = qp_com%box%y1 + qp_com%margin%y1
  graph%x2 = qp_com%box%x2 - qp_com%margin%x2  
  graph%y2 = qp_com%box%y2 - qp_com%margin%y2
endif

! only set graph position if the coords look reasonable

if (graph%x1 < graph%x2 .and. graph%y1 < graph%y2) then
  call qp_set_graph_position_basic (graph%x1, graph%x2, graph%y1, graph%y2)
endif

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_clear_page
!
! Subroutine to clear all drawing from the page.
!-

subroutine qp_clear_page

type (qp_rect_struct) border  

! Clear all graphics.

call qp_clear_page_basic

! For gif must paint a white background.

if (qp_com%page_type(1:3) == 'GIF') then
  call qp_paint_rectangle (qp_com%page%x1, qp_com%page%x2, &
                           qp_com%page%y1, qp_com%page%y2, &
                           color = white$, fill_pattern = solid_fill$)
endif

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_clear_box
!
! Subroutine to clear the current box on the page.
!-

subroutine qp_clear_box

! Clear the box

call qp_paint_rectangle (qp_com%box%x1, qp_com%box%x2, &
                         qp_com%box%y1, qp_com%box%y2, &
                           color = white$, fill_pattern = solid_fill$)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_paint_rectangle (x1, x2, y1, y2, units, color, fill_pattern)
!
! Subroutine to paint a rectangular region a specified color.
! The default color is the background color (white$).
!
! Input:
!   x1  -- Real(rp): Left edge
!   x2  -- Real(rp): Right edge
!   y1  -- Real(rp): Bottom edge
!   y2  -- Real(rp): Top edge.
!   units        -- Character(*), optional: Units of returned numbers.
!                     Default = 'INCH'   (unless default_set_units is changed).
!   color        -- Integer, optional: Color to paint the rectangle.
!                     Default is to use the symbol color.
!   fill_pattern -- Integer, optional: Fill pattern. 
!                   Default is to use the symbol color.
!-

subroutine qp_paint_rectangle (x1, x2, y1, y2, units, color, fill_pattern)
              
implicit none

real(rp) x1, x2, y1, y2
real(rp) x1_inch, x2_inch, y1_inch, y2_inch
character(*), optional :: units
integer, optional :: color, fill_pattern

!

if (x1 == x2 .or. y1 == y2) return

qp_com%dflt_units = dflt_set$
call qp_to_inch_abs (x1, y1, x1_inch, y1_inch, units)
call qp_to_inch_abs (x2, y2, x2_inch, y2_inch, units)
qp_com%dflt_units = dflt_draw$
 
call qp_paint_rectangle_basic (x1_inch, x2_inch, y1_inch, y2_inch, &
              integer_option(qp_com%symbol%color, color), &
              integer_option(qp_com%symbol%fill_pattern, fill_pattern), &
              qp_com%page_type)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_box (ix, iy, ix_tot, iy_tot)
!                                              
! Subroutine to set the box on the physical page.
! This routine divides the page into a grid of boxes. 
! There are ix_tot boxes horizontally and iy_tot boxes vertically.
! the (ix, iy) = (1, 1) box is the lower left box and the
! (ix, iy) = (ix_tot, iy_tot) box is at the upper right.
!
! Use this routine with qp_set_margin and qp_set_page_border.
!
! Input:
!   ix_tot, iy_tot -- Integer: X and Y box divisions.
!   ix, iy         -- Integer: Index for box to be used.
!
!-

subroutine qp_set_box (ix, iy, ix_tot, iy_tot)

implicit none

integer ix, iy, ix_tot, iy_tot
real(rp) w_x, w_y

! calculate the placement of the box

w_x = qp_com%page%x2 - qp_com%page%x1 - qp_com%border%x1 - qp_com%border%x2
w_y = qp_com%page%y2 - qp_com%page%y1 - qp_com%border%y1 - qp_com%border%y2

qp_com%box%x1 = qp_com%page%x1 + qp_com%border%x1 + w_x * (ix - 1) / ix_tot
qp_com%box%y1 = qp_com%page%y1 + qp_com%border%y1 + w_y * (iy - 1) / iy_tot

qp_com%box%x2 = qp_com%box%x1 + w_x / ix_tot
qp_com%box%y2 = qp_com%box%y1 + w_y / iy_tot

! finally calculate placement of the graph within the box
! set qp_com%subgraph_on so qp_set_graph_limits will calculate the 
! graph boundry

qp_com%subgraph_on = .false.   
call qp_set_graph_limits 

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_page_border_to_box ()
!
! Subroutine to set the page border to correspond to the region of the
! current box. This allows qp_set_box to subdivide the current box.
!-

subroutine qp_set_page_border_to_box ()

implicit none

qp_com%border%x1 = qp_com%box%x1 - qp_com%page%x1
qp_com%border%y1 = qp_com%box%y1 - qp_com%page%y1

qp_com%border%x2 = qp_com%page%x2 - qp_com%box%x2
qp_com%border%y2 = qp_com%page%y2 - qp_com%box%y2
  
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_page_border (x1_b, x2_b, y1_b, y2_b, units)
!
! Subroutine to set the border around the physical page.
!               
! Input:
!   x1_b  -- Real(rp): Left border.
!   y1_b  -- Real(rp): Bottom border.
!   x2_b  -- Real(rp): Right border. 
!   y2_b  -- Real(rp): Top border.
!   units -- Character(*), optional: border units:
!               '%PAGE'  - Percent of page.
!               'MM'     - milimeters.
!               'INCH'   - Inches (default).
!               'POINTS' - Points.
!-

subroutine qp_set_page_border (x1_b, x2_b, y1_b, y2_b, units)

implicit none

real(rp) x1_b, y1_b, x2_b, y2_b
character(*), optional :: units

!

qp_com%dflt_units = dflt_set$

call qp_to_inch_rel (x1_b, y1_b, qp_com%border%x1, qp_com%border%y1, units)
call qp_to_inch_rel (x2_b, y2_b, qp_com%border%x2, qp_com%border%y2, units)

qp_com%dflt_units = dflt_draw$

call qp_set_graph_limits

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_from_inch_rel (x_inch, y_inch, x, y, units)
!
! Subroutine to convert from a relative position (an offset) in inches
! to other units.
!
! Input:
!   x_inch, y_inch -- Real(rp): lengths in inches.
!   units  -- Character(*), optional: Units of x_inch and y_inch
!
! Output:
!   x, y   -- Real(rp): Lengths to convert.
!-

Subroutine qp_from_inch_rel (x_inch, y_inch, x, y, units)

implicit none

real(rp) x, y
real(rp) x_inch, y_inch, dx, dy

character(*), optional :: units
character(8) u_type, region, corner

!  

call qp_split_units_string (u_type, region, corner, units)

if (u_type == 'MM') then
  x = x_inch * 25.4
  y = y_inch * 25.4
elseif (u_type == 'INCH') then
  x = x_inch 
  y = y_inch 
elseif (u_type == 'POINTS') then
  x = x_inch * 72
  y = y_inch * 72
elseif (u_type == '%' .and. region == 'PAGE') then
  x = x_inch / (qp_com%page%x2 - qp_com%page%x1)
  y = y_inch / (qp_com%page%y2 - qp_com%page%y1)
elseif (u_type == '%' .and. region == 'GRAPH') then
  x = x_inch / (qp_com%graph%x2 - qp_com%graph%x1)
  y = y_inch / (qp_com%graph%y2 - qp_com%graph%y1)
elseif (u_type == '%' .and. region == 'BOX') then
  x = x_inch / (qp_com%box%x2 - qp_com%box%x1)
  y = y_inch / (qp_com%box%y2 - qp_com%box%y1)
elseif (u_type == 'DATA') then
  dx = x_inch / (qp_com%graph%x2 - qp_com%graph%x1)
  if (qp_com%plot%xx%type == 'LOG') then
    x = (qp_com%plot%xx%max / qp_com%plot%xx%min) ** dx
  else
    x = dx * (qp_com%plot%xx%max - qp_com%plot%xx%min) 
  endif
  dy = y_inch / (qp_com%graph%y2 - qp_com%graph%y1)
  if (qp_com%plot%yy%type == 'LOG') then
    y = (qp_com%plot%yy%max / qp_com%plot%yy%min) ** dy
  else
    y = dy * (qp_com%plot%yy%max - qp_com%plot%yy%min) 
  endif
endif

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_from_inch_abs (x_inch, y_inch, x, y, units)
!
! Subroutine to convert to absolute position (x, y) from inches referenced
! to the Left Bottom corner of the page
!
! Input:
!   x_inch, y_inch -- Real(rp): Position in inches from LB corner of the page.
!   units  -- Character(*), optional: Units of x and y
!
! Output:
!   x, y   -- Real(rp): Position in new coords
!-

subroutine qp_from_inch_abs (x_inch, y_inch, x, y, units)

implicit none

type (qp_rect_struct) ref

real(rp) x, y, x0, y0
real(rp) x_inch, y_inch

character(*), optional :: units
character(8) u_type, region, corner

! Init

call qp_split_units_string (u_type, region, corner, units)

! Data units case.

if (u_type == 'DATA') then
  x0 = x_inch - qp_com%graph%x1
  y0 = y_inch - qp_com%graph%y1
  call qp_from_inch_rel (x0, y0, x, y, units)
  if (qp_com%plot%xx%type == 'LOG') then
    x = x * qp_com%plot%xx%min
  else
    x = x + qp_com%plot%xx%min
  endif
  if (qp_com%plot%yy%type == 'LOG') then
    y = y * qp_com%plot%yy%min
  else
    y = y + qp_com%plot%yy%min
  endif
  return
endif

! All other cases

if (region == 'PAGE') then
  ref = qp_com%page
elseif (region == 'BOX') then
  ref = qp_com%box
else
  ref = qp_com%graph
endif

if (corner(1:1) == 'L') then
  x0 = x_inch - ref%x1
else
  x0 = x_inch - ref%x2
endif

if (corner(2:2) == 'B') then
  y0 = y_inch - ref%y1 
else
  y0 = y_inch - ref%y2
endif

call qp_from_inch_rel (x0, y0, x, y, units)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_to_inch_rel (x, y, x_inch, y_inch, units)
!
! Subroutine to convert a relative (x, y) into inches.
!
! Input:
!   x, y   -- Real(rp): Lengths to convert.
!   units  -- Character(*), optional: Units of x and y
!
! Output:
!   x_inch, y_inch -- Real(rp): lengths in inches.
!-

subroutine qp_to_inch_rel (x, y, x_inch, y_inch, units)

implicit none

real(rp) x, y
real(rp) x_inch, y_inch, dgx, dgy

character(*), optional :: units
character(8) u_type, region, corner
character(16) :: r_name = 'qp_to_inch_rel'

!  

call qp_split_units_string (u_type, region, corner, units)

if (u_type == 'MM') then
  x_inch = x / 25.4
  y_inch = y / 25.4
elseif (u_type == 'INCH') then
  x_inch = x 
  y_inch = y 
elseif (u_type == 'POINTS') then
  x_inch = x / 72
  y_inch = y / 72
elseif (u_type == '%' .and. region == 'PAGE') then
  x_inch = x * (qp_com%page%x2 - qp_com%page%x1)
  y_inch = y * (qp_com%page%y2 - qp_com%page%y1)
elseif (u_type == '%' .and. region == 'GRAPH') then
  x_inch = x * (qp_com%graph%x2 - qp_com%graph%x1)
  y_inch = y * (qp_com%graph%y2 - qp_com%graph%y1)
elseif (u_type == '%' .and. region == 'BOX') then
  x_inch = x * (qp_com%box%x2 - qp_com%box%x1)
  y_inch = y * (qp_com%box%y2 - qp_com%box%y1)
elseif (u_type == 'DATA') then
  dgx = qp_com%graph%x2 - qp_com%graph%x1
  if (qp_com%plot%xx%type == 'LOG') then
    x_inch = log(x) * dgx / (log(qp_com%plot%xx%max) - log(qp_com%plot%xx%min))
  else
    x_inch = x * dgx / (qp_com%plot%xx%max - qp_com%plot%xx%min)
  endif
  dgy = qp_com%graph%y2 - qp_com%graph%y1
  if (qp_com%plot%yy%type == 'LOG') then
    y_inch = log(y) * dgy / (log(qp_com%plot%yy%max) - log(qp_com%plot%yy%min))
  else
    y_inch = y * dgy / (qp_com%plot%yy%max - qp_com%plot%yy%min)
  endif
endif

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_to_inches_rel (x, y, x_inch, y_inch, units)
!
! Subroutine to convert a relative (x, y) into inches.
!
! Input:
!   x(:), y(:) -- Real(rp): Lengths to convert.
!   units      -- Character(*), optional: Units of x and y
!
! Output:
!   x_inch(:), y_inch(:) -- Real(rp): lengths in inches.
!-

subroutine qp_to_inches_rel (x, y, x_inch, y_inch, units)

implicit none

real(rp) x(:), y(:)
real(rp) x_inch(:), y_inch(:)

integer i
character(*), optional :: units

!

do i = 1, size(x)
  call qp_to_inch_rel (x(i), y(i), x_inch(i), y_inch(i), units)
enddo

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_to_inch_abs (x, y, x_inch, y_inch, units)
!
! Subroutine to convert an absolute position (x, y) into inches referenced
! to the Left Bottom corner of the page.
!
! Input:
!   x, y   -- Real(rp): Position to convert.
!   units  -- Character(*), optional: Units of x and y
!
! Output:
!   x_inch, y_inch -- Real(rp): Position in inches referenced to the page.
!-

subroutine qp_to_inch_abs (x, y, x_inch, y_inch, units)

implicit none

type (qp_rect_struct) ref

real(rp) x, y, x0, y0
real(rp) x_inch, y_inch

character(*), optional :: units
character(8) u_type, region, corner

! Init

call qp_split_units_string (u_type, region, corner, units)

! Data units case

if (u_type == 'DATA') then
  if (qp_com%plot%xx%type == 'LOG') then
    x0 = x / qp_com%plot%xx%min
  else
    x0 = x - qp_com%plot%xx%min
  endif
  if (qp_com%plot%yy%type == 'LOG') then
    y0 = y / qp_com%plot%yy%min
  else
    y0 = y - qp_com%plot%yy%min
  endif
  call qp_to_inch_rel (x0, y0, x_inch, y_inch, units)
  x_inch = x_inch + qp_com%graph%x1
  y_inch = y_inch + qp_com%graph%y1
  return
endif

! Other cases.

call qp_to_inch_rel (x, y, x_inch, y_inch, units)

if (region == 'PAGE') then
  ref = qp_com%page
elseif (region == 'BOX') then
  ref = qp_com%box
else
  ref = qp_com%graph
endif

if (corner(1:1) == 'L') then
  x_inch = x_inch + ref%x1 
else
  x_inch = x_inch + ref%x2 
endif

if (corner(2:2) == 'B') then
  y_inch = y_inch + ref%y1 
else
  y_inch = y_inch + ref%y2 
endif

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_to_inches_abs (x, y, x_inch, y_inch, units)
!
! Subroutine to convert an absolute position (x, y) into inches referenced
! to the Left Bottom corner of the page.
!
! Input:
!   x(:), y(:) -- Real(rp): Position to convert.
!   units      -- Character(*), optional: Units of x and y
!
! Output:
!   x_inch(:), y_inch(:) -- Real(rp): Position in inches referenced to the page.
!-

subroutine qp_to_inches_abs (x, y, x_inch, y_inch, units)

implicit none

real(rp) x(:), y(:)
real(rp) x_inch(:), y_inch(:)

integer i
character(*), optional :: units

!

do i = 1, size(x)
  call qp_to_inch_abs (x(i), y(i), x_inch(i), y_inch(i), units)
enddo

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_convert_point_rel (x_in, y_in, units_in, x_out, y_out, units_out)
!
! Subroutine to convert a (x, y) point from from
! one set of relative units to another.
!
! Input:
!   x_in, y_in  -- Real(rp): Input coordinates.
!   units_in    -- Character(*): Units of the input coordinates.
!   units_out   -- Character(*): Units of the output coordinates.
!
! Output:
!   x_out, y_out  -- Real(rp): output coordinates.
!-

subroutine qp_convert_point_rel (x_in, y_in, units_in, x_out, y_out, units_out)

implicit none

real(rp) x_in, y_in, x_out, y_out
character(*) units_in, units_out

!

call qp_to_inch_rel   (x_in, y_in, x_out, y_out, units_in)
call qp_from_inch_rel (x_out, y_out, x_out, y_out, units_out)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_convert_point_abs (x_in, y_in, units_in, x_out, y_out, units_out)
!
! Subroutine to convert a (x, y) point from from
! one set of absolute units to another.
!
! Input:
!   x_in, y_in  -- Real(rp): Input coordinates.
!   units_in    -- Character(*): Units of the input coordinates.
!   units_out   -- Character(*): Units of the output coordinates.
!
! Output:
!   x_out, y_out  -- Real(rp): output coordinates.
!-

subroutine qp_convert_point_abs (x_in, y_in, units_in, x_out, y_out, units_out)

implicit none

real(rp) x_in, y_in, x_out, y_out
character(*) units_in, units_out

!

call qp_to_inch_abs   (x_in, y_in, x_out, y_out, units_in)
call qp_from_inch_abs (x_out, y_out, x_out, y_out, units_out)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_convert_rectangle_rel (rect1, rect2)
!
! Subroutine to convert a "rectangle" (structure of 4 points) from
! one set of relative units to another.
!
! Input:
!   rect1 -- Qp_rect_struct: Input rectangle.
!     %x1, %x2, %y1, %y2 -- Real(rp): The 4 points.
!     %units             -- Units of the input rectangle.
!   rect2 -- Qp_rect_struct: 
!     %units             -- Units of the output rectangle. 
!
! Output:
!   rect2 -- Qp_rect_struct: Output rectangle.
!     %x1, %x2, %y1, %y2 -- Real(rp): The 4 points.
!-

subroutine qp_convert_rectangle_rel (rect1, rect2)

implicit none

type (qp_rect_struct) rect1, rect2

!

call qp_convert_point_rel (rect1%x1, rect1%y1, rect1%units, &
                           rect2%x1, rect2%y1, rect2%units)

call qp_convert_point_rel (rect1%x2, rect1%y2, rect1%units, &
                           rect2%x2, rect2%y2, rect2%units)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_join_units_string (u_type, region, corner, units)
!
! Subroutine to form a units string from its components.
!
! Input:
!   u_type -- Character(*): Unit type.
!   region -- Character(*): Region.
!   corner -- Character(*): Origin Corner.
!
! Output:
!   units -- Character(*) Joined units string. Eg: 'DATA/GRAPH/LB'.
!-

subroutine qp_join_units_string (u_type, region, corner, units) 

implicit none

character(*) u_type, region, corner
character(*) units

!

units = u_type // '/' // region // '/' // corner

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_split_units_string (u_type, region, corner, units)
!
! Subroutine to split a units string into its components.
!
! Input:
!   units -- Character(*), optional: Units string to split.
!               If not present then 'DATA/GRAPH/LB' is assumed.
!
! Output:
!   u_type -- Character(*): Unit type.
!   region -- Character(*): Region.
!   corner -- Character(*): Origin Corner.
!-

subroutine qp_split_units_string (u_type, region, corner, units) 

implicit none

integer i, ix

character(*) u_type, region, corner
character(*), optional :: units
character(20) u
character(8) dflt_units(3)
character(28) :: r_name = 'qp_split_units_string'

! Default

if (qp_com%dflt_units == dflt_draw$) then
  dflt_units = qp_com%dflt_draw_units
else
  dflt_units = qp_com%dflt_set_units
endif

u_type = dflt_units(1)
region = dflt_units(2)
corner = dflt_units(3)

!

if (.not. present(units)) return
u = units

! strip '/'

do i = 1, 2
  ix = index(u, '/')
  if (ix /= 0) u(ix:ix) = ' '
enddo

! get u_type

call string_trim (u, u, ix)
if (ix == 0) return
if (u(1:1) == '%') ix = 1

u_type = u(:ix)

if (all(u_type /= (/ 'DATA  ', 'MM    ', 'INCH  ', 'POINTS', '%     ' /))) then
  call out_io (s_fatal$, r_name, 'BAD UNITS TYPE: "' // trim(units) // '"')
  call err_exit
endif

! get region

call string_trim (u(ix+1:), u, ix)
if (ix == 0) return
region = u(:ix)

if (all(region /= (/ 'PAGE ', 'BOX  ', 'GRAPH' /))) then
  call out_io (s_fatal$, r_name, 'BAD REGION: "' // trim(units) // '"')
  call err_exit
endif

! get corner

call string_trim (u(ix+1:), u, ix)
if (ix == 0) return
corner = u(:ix)

if (all(corner /= (/ 'LB', 'LT', 'RB', 'RT' /))) then
  call out_io (s_fatal$, r_name, 'BAD CORNER: "' // trim(units) // '"')
  call err_exit
endif

call string_trim (u(ix+1:), u, ix)
if (ix /= 0) then
  call out_io (s_fatal$, r_name, 'EXTRA CHARACTERS IN UNITS: "' // trim(units) // '"')
  call err_exit
endif

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_graph_placement (x1_marg, x_graph_len, y1_marg, 
!                                                        y_graph_len, units)
!
! Subroutine to set the placement of the current graph inside the box. 
! This routine can be used in place of qp_set_margin.
!
! Note: Subsequent to a call to this routine, if the size of the box is
! changed, then the graph size will be modified to keep the margins constant.
! Thus, qp_set_box and qp_set_page_border should be called first.
!
! Input:
!   x1_marg, y1_marg -- Real(rp): offset from the lower left corner 
!                  of the box to the lower left corner of the plotting region.
!   x_graph_len, y_graph_len -- Real(rp): Size of graph.
!   units    -- Character(*), optional: Units of the margins.
!                 Default is: 'DATA/GRAPH'
!                 See quick_plot writeup for more details.
!-

subroutine qp_set_graph_placement (x1_marg, x_graph_len, &
                                                y1_marg, y_graph_len, units)

implicit none

real(rp) x1_marg, y1_marg, x_graph_len, y_graph_len
real(rp) x_len, y_len
character(*), optional :: units

!

qp_com%dflt_units = dflt_set$

call qp_to_inch_rel (x1_marg, y1_marg, qp_com%margin%x1, &
                                        qp_com%margin%y1, units)

call qp_to_inch_rel (x_graph_len, y_graph_len, x_len, y_len, units)
qp_com%margin%x2 = (qp_com%box%x2 - qp_com%box%x1) - x_len - qp_com%margin%x1
qp_com%margin%y2 = (qp_com%box%y2 - qp_com%box%y1) - y_len - qp_com%margin%y1

qp_com%dflt_units = dflt_draw$

call qp_set_graph_limits

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_margin (x1_marg, x2_marg, y1_marg, y2_marg, units)
!
! Subroutine to set up the margins from the sides of the box (see qp_set_box)
! to the edges of the actual graph. Alternatively, qp_set_graph_placement
! can be used instead of this routine.
!
! Input:
!   x1_marg, y1_marg -- Real(rp): offset from the lower left corner 
!                  of the box to the lower left corner of the plotting region.
!   x2_marg, y2_marg -- Real(rp): offset from the upper right corner
!                  of the box to the upper right corner of the plotting region.
!   units    -- Character(*), optional: Units of the margins.
!                 Default is: 'DATA/GRAPH'
!                 See quick_plot writeup for more details.
!
! Use this routine with qp_set_box and qp_set_page_border.
!-

subroutine qp_set_margin (x1_marg, x2_marg, y1_marg, y2_marg, units)

implicit none

real(rp) x1_marg, y1_marg, x2_marg, y2_marg
character(*), optional :: units

!

qp_com%dflt_units = dflt_set$

call qp_to_inch_rel (x1_marg, y1_marg, qp_com%margin%x1, &
                                        qp_com%margin%y1, units)
call qp_to_inch_rel (x2_marg, y2_marg, qp_com%margin%x2, &
                                        qp_com%margin%y2, units)

qp_com%dflt_units = dflt_draw$

call qp_set_graph_limits

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_rectangle (x1, x2, y1, y2, units, color, width, 
!                                                          style, clip)
!
! Subroutine to draw a rectangular box.
!
! Input:
!   x1, y1 -- Real(rp): (x, y) corner of box.
!   x2, y2 -- Real(rp): (x, y) opposite corner of box.
!   units  -- Character(*), optional: Units of x and y.
!                 Default is: 'DATA/GRAPH/LB'
!                 See quick_plot writeup for more details.
!   color   -- Integer, optional: Color index for the box
!   width   -- Integer, optional: Width of the line. Default = 1
!   style   -- Integer, optional: Line style. 
!   clip    -- Logical, optional: Clip at the graph boundary?
!-

subroutine qp_draw_rectangle (x1, x2, y1, y2, units, color, width, style, clip)

implicit none

real(rp) x1, y1, x2, y2

integer, optional :: color, width, style

character(*), optional :: units

logical, optional :: clip

! If there is zero thickness nothing has to be drawn.

if (x1 == x2 .or. y1 == y2) return

call qp_draw_polyline ((/ x1, x1, x2, x2, x1 /), (/ y1, y2, y2, y1, y1 /), &
                                            units, width, color, style, clip)

end subroutine
                                     
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_symbol (x, y, units, type, height, color, 
!                                            fill_pattern, line_width, clip)
!
! Draws a symbol at (x, y) 
!
! Input:
!   x, y         -- Real(rp): Symbol coordinates in data units.
!                     x and y may be vectors.
!   units        -- Character(*), optional: Units of (x, y)
!   type         -- Integer, optional: Symbol type. 
!   height       -- Real(rp), optional: Size of the symbol.
!   color        -- Integer, optional: Symbol color.
!   fill_pattern -- Integer, optional: fill pattern.
!   line_width   -- Integer, optional: Line width.
!   clip         -- Logical, optional: Clip at the graph boundary?
!-

subroutine qp_draw_symbol (x, y, units, type, height, color, &
                                                 fill_pattern, line_width, clip)

implicit none
              
integer, optional :: type, color, fill_pattern, line_width

real(rp) x, y, x_inch, y_inch
real(rp), optional :: height

character(*), optional :: units

logical, optional :: clip

!

call qp_save_state (.true.)

call qp_set_symbol_attrib (type, height, color, fill_pattern, line_width, clip)
call qp_to_inch_abs (x, y, x_inch, y_inch, units)
call qp_draw_symbol_basic (x_inch, y_inch, qp_com%symbol%type)

call qp_restore_state

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_symbols (x, y, units, type, height, color, 
!                                 fill_pattern, line_width, clip, symbol_every)
!
! Draws a symbol at the (x, y) points. 
! Data units are assumed.
!
! Input:
!   x, y         -- Real(rp): Symbol coordinates in data units.
!                     x and y may be vectors.
!   type         -- Integer, optional: Symbol type. 
!   height       -- Real(rp), optional: Size of the symbol.
!   color        -- Integer, optional: Symbol color.
!   fill_pattern -- Integer, optional: fill pattern.
!   line_width   -- Integer, optional: Line width.
!   clip         -- Logical, optional: Clip at the graph boundary?
!   symbol_every -- Integer, optional: 
!                               0  --> Do not draw symbols.
!                               1  --> Draw symbols (Default).
!                               2  --> Draw every 2nd symbol.
!                               etc.
!-


subroutine qp_draw_symbols (x, y, units, type, height, color, &
                                      fill_pattern, line_width, clip, symbol_every)

implicit none

integer, optional :: type, color, fill_pattern, line_width
integer, optional :: symbol_every
integer i, i_skip

real(rp) x(:), y(:)
real(rp), optional :: height

character(*), optional :: units 
character(16) :: r_name = 'qp_draw_symbols'

logical, optional :: clip

!

if (size(x) /= size(y)) then
  call out_io (s_fatal$, r_name, 'X, Y COORD VECTORS HAVE UNEQUAL LENGTH!')
  call err_exit
endif

i_skip = 1
if (present(symbol_every)) then
  i_skip = symbol_every
endif

if (i_skip < 1) return

do i = 1, size(x)
  if (mod(i-1, i_skip) /= 0) cycle
  call qp_draw_symbol (x(i), y(i), units, type, height, color, fill_pattern, &
                                                          line_width, clip)
enddo


end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_graph (title)
!
! Subroutine to set certain graph attributes
!
! Input:
!   title  -- Character(*), optional: Graph title.
!-

subroutine qp_set_graph (title)

character(*), optional :: title

!

if (present(title)) qp_com%plot%title = title

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_graph (x, y, x_lab, y_lab, title, &
!                                        draw_line, symbol_every, clip)
!
! Subroutine to plot data, axes with labels, a grid, and a title.
! Use the routine qp_draw_data just to draw the data.
! Use the routine qp_draw_axes just to draw the axes.
!
! Input:
!   x(:), y(:)    -- Real(rp): data arrays.
!   x_lab, y_lab  -- Character(*), optional: x and y axes labels.
!   title         -- Character(*), optional: Graph Title.
!   draw_line     -- Logical, optional: Default = T.
!   symbol_every  -- Integer, optional: 
!                               0  --> Do not draw symbols.
!                               1  --> Draw symbols (Default).
!                               2  --> Draw every 2nd symbol.
!                               etc.
!   clip          -- Logical, optional: Clip at the graph boundary?
!
! See:
!   qp_open_page          for setting up the plot page.
!   qp_set_box            for setting up the box within the page.
!   qp_set_margin         for setting up the graph margins within the box.
!   qp_set_symbol_attrib  for setting the symbol attributes.
!   qp_set_line_attrib    for setting the line attributes.
!   qp_set_axis           for setting up the axis scales.
!-

subroutine qp_draw_graph (x_dat, y_dat, x_lab, y_lab, title, &
                                            draw_line, symbol_every, clip)

implicit none      

real(rp) x_dat(:), y_dat(:)                        

integer, optional :: symbol_every
logical, optional :: draw_line, clip

character(*), optional :: x_lab, y_lab, title  
character(16) :: r_name = 'qp_draw_graph'

! Error check

if (qp_com%plot%xx%max == qp_com%plot%xx%min) then
  call out_io (s_fatal$, r_name, 'X_MAX = X_MIN')
  call err_exit
endif

if (qp_com%plot%yy%max == qp_com%plot%yy%min) then
  call out_io (s_fatal$, r_name, 'Y_MAX = Y_MIN')
  call err_exit
endif
           
! Draw the axes and the data

call qp_draw_axes (x_lab, y_lab, title)
call qp_draw_data (x_dat, y_dat, draw_line, symbol_every, clip)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_data (x, y, draw_line, symbol_every, clip)
!
! Subroutine to plot data.
! See qp_draw_graph for a routine that will draw the axes as well.
! See qp_draw_axes for a routine to only draw the axes.
!
! Input:
!   x(:), y(:)    -- Real(rp): data arrays.
!   draw_line     -- Logical, optional: Default = T.
!   symbol_every  -- Integer, optional: 
!                               0  --> Do not draw symbols.
!                               1  --> Draw symbols (Default).
!                               2  --> Draw every 2nd symbol.
!                               etc.
!   clip          -- Logical, optional: Clip at the graph boundary?
!                      Default is set by qp_set_line_attrib.
!-

subroutine qp_draw_data (x_dat, y_dat, draw_line, symbol_every, clip)

implicit none      

real(rp) x_dat(:), y_dat(:)                        

integer, optional :: symbol_every
logical, optional :: draw_line, clip

! init

call qp_save_state (.true.)

! plot a polyline

if (logic_option (.true., draw_line)) then
  call qp_set_line_attrib ('PLOT', clip = clip)
  call qp_draw_polyline_no_set (x_dat, y_dat)
endif

! plot symbols

call qp_draw_symbols (x_dat, y_dat, symbol_every = symbol_every, clip = clip)

! 

call qp_restore_state

end subroutine      

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+       
! Subroutine qp_draw_axes (x_lab, y_lab, title, draw_grid)
!
! Subroutine to plot the axes, title, etc. of a plot.
!
! Input:
!   x_lab, y_lab  -- Character(*), optional: x and y axes labels.
!   title         -- Character(*), optional: Graph Title.
!   draw_grid     -- Logical, optional: Draw a grid? 
!                      Default is set by qp_set_graph_attrib.
!-

subroutine qp_draw_axes (x_lab, y_lab, title, draw_grid)

implicit none

character(*), optional :: x_lab, y_lab, title  
logical, optional :: draw_grid

!

call qp_save_state (.true.)

if (present(x_lab)) qp_com%plot%xx%label = x_lab
if (present(y_lab)) qp_com%plot%yy%label = y_lab
if (present(title)) qp_com%plot%title = title

if (logic_option (qp_com%plot%draw_grid, draw_grid)) call qp_draw_grid

call qp_draw_x_axis ('X',  0.0_rp)
call qp_draw_x_axis ('X2', 1.0_rp)
call qp_draw_y_axis ('Y',  0.0_rp)
call qp_draw_y_axis ('Y2', 1.0_rp)

if (qp_com%plot%draw_title) call qp_draw_graph_title (qp_com%plot%title)

call qp_restore_state

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_graph_title (title)
!
! Subroutine to draw the title for a graph.
!
! Input:
!   title -- Character(*): Title
!-

subroutine qp_draw_graph_title (title)

implicit none

character(*) title

real(rp) xt, yt                 

!

if (len_trim(title) == 0) return

call qp_set_text_attrib ('GRAPH_TITLE')  
call qp_to_inch_rel (0.5_rp, 1.0_rp, xt, yt, '%GRAPH')
call qp_draw_text_no_set (title, xt, yt+0.05, 'INCH', 'CB')

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_histogram (x_dat, y_dat, fill_color, fill_pattern, draw_line, clip)
!
! Subroutine to draw a histogram.
! The histogram is drawn as a series of rectangles, one for each
! (x_dat(i), y_dat(i)) data point. 
! Horizontally, the i^th bin extends between 
!    (x_dat(i-1)+x_dat(i)) / 2   to   (x_dat(i)+x_dat(i+1))/2
! For the first data point, The rectangle is drawn centered on
! x_dat(1). Similarly, for the last data point.
!
! Note: Use qp_draw_axes to draw the axes.
!
! Input:
!   x_dat(:), y_dat(:) -- Real(rp): Rectangle Data.
!   fill_color   -- Integer, optional: Color to fill the rectangles.
!                     If fill_color is 0 (white$) no color is added.
!                     Default = white$
!   fill_pattern -- Integer, optional: Default is set by symbol fill pattern.
!   draw_line    -- Logical, optional: Draw an outline of the rectangles?
!                     Default is True.
!   clip         -- Logical, optional: Clip at the graph boundary?
!                     Default is set by qp_set_line_attrib.
!-

subroutine qp_draw_histogram (x_dat, y_dat, fill_color, fill_pattern, draw_line, clip) 

implicit none

integer i, n, n_min, n_max
integer, optional :: fill_color, fill_pattern

real(rp) x_dat(:), y_dat(:)
real(rp) :: xh(2*size(x_dat)+2), yh(2*size(x_dat)+2)

logical, optional :: draw_line, clip

character(16) :: r_name = 'qp_draw_histogram'

! error check

if (qp_com%plot%xx%max == qp_com%plot%xx%min) then
  call out_io (s_fatal$, r_name, 'X_MAX = X_MIN')
  call err_exit
endif

if (qp_com%plot%yy%max == qp_com%plot%yy%min) then
  call out_io (s_fatal$, r_name, 'Y_MAX = Y_MIN')
  call err_exit
endif

! Find which points are inside the plot horizontally.

n_min = 1; n_max = 0
do i = 1, size(x_dat)
  if (i > 1) then
    if (x_dat(i-1) > x_dat(i)) then
      call out_io (s_error$, r_name, 'X_DAT POINTS NOT ORDERED IN INCREASING VALUE.')
      return
    endif
  endif
  if (x_dat(i) < qp_com%plot%xx%min) n_min = i + 1
  if (x_dat(i) > qp_com%plot%xx%max .and. n_max == 0) n_max = i - 1
enddo
if (n_max == 0) n_max = size(x_dat)

n = n_max - n_min + 1 ! number of rectangles to be plotted
if (n < 2) return

! Compute the rectangle boundary points

xh(1)     = max ((3*x_dat(n_min) - x_dat(n_min+1))/2, qp_com%plot%xx%min)
xh(2)     = xh(1)

xh(3:2*n-1:2) = (x_dat(n_min:n_max-1) + x_dat(n_min+1:n_max)) / 2
xh(4:2*n:2)   = (x_dat(n_min:n_max-1) + x_dat(n_min+1:n_max)) / 2

xh(2*n+1) = min ((3*x_dat(n_max) - x_dat(n_max-1))/2, qp_com%plot%xx%max)
xh(2*n+2) = xh(2*n+1)


yh(1) = 0

yh(2:2*n:2)   = y_dat(n_min:n_max)
yh(3:2*n+1:2) = y_dat(1:n)

yh(2*n+2) = 0

if (logic_option (qp_com%clip, clip)) then
  where (yh < qp_com%plot%yy%min) yh = qp_com%plot%yy%min
  where (yh > qp_com%plot%yy%max) yh = qp_com%plot%yy%max
endif

! Draw

call qp_save_state (.true.)

call qp_set_line_attrib ('PLOT', clip = .false.)

if (logic_option(.true., draw_line)) call qp_draw_polyline_no_set (xh, yh)

if (integer_option(white$, fill_color) /= white$) then
  do i = 1, n
    call qp_paint_rectangle (xh(2*i), xh(2*i+1), yh(1), yh(2*i), 'DATA', &
                                                       fill_color, fill_pattern)
  enddo
endif

call qp_restore_state

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_text_legend (text, x_origin, y_origin, units)
!
! Subroutine to draw a legend of lines of text.
! Note: If (x_origin, y_origin) is not given then by default the legend 
!   is drawn starting just to the right of the graph.
! Note: If the units are not specified then data units are assumed.
!
! Input:
!   text(:)  -- Character(*): Array of text lines to print.
!   x_origin -- Real(rp), optional: x-postion of start of the first line.
!   y_origin -- Real(rp), optional: y-postion of start of the first line.
!   units    -- Character(*), optional: Units of x_origin, y_origin.
!                 Default is: 'DATA/GRAPH/LB'
!                 See quick_plot writeup for more details.
!-

subroutine qp_draw_text_legend (text, x_origin, y_origin, units)

implicit none

real(rp), optional :: x_origin, y_origin
real(rp) xc, yc, height

character(*) text(:)
character(*), optional :: units

integer i 

!

call qp_save_state (.true.)

call qp_set_text_attrib ('LEGEND')
height = qp_text_height_to_inches(qp_com%this_text%height) 

if (present(x_origin)) then
  call qp_to_inch_abs (x_origin, y_origin, xc, yc, units)
else
  call qp_to_inch_abs (0.2_rp, 0.0_rp, xc, yc, 'INCH/GRAPH/RT')
endif

yc = yc - 1.1 * height
do i = 1, size(text)
  call qp_draw_text_no_set (text(i), xc, yc, 'INCH/PAGE/LB')
  yc = yc - 1.5 * height
enddo

call qp_restore_state

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_curve_legend (origin, text_offset, line_length, 
!                  line, symbol, text, draw_line, draw_symbol, draw_text)
!
! Subroutine to draw a legend with each line in the legend having
!   a line, a symbol, some text.
!
! Note: If origin is not given then by default the legend is drawn in the
!   upper left hand corner of the graph.
!
! Input:
!   origin      -- qp_point_struct: Position of upper left edge of the legend.
!   text_offset -- Real(rp), optional: Horizontal offset between the line and the text.
!   line_length -- Real(rp), optional: Length of the line in points.
!   line(:)     -- qp_line_struct, optional: Array of lines.
!                    Set line(i)%width < 0 to suppress drawing of the i^th line
!   symbol(:)   -- qp_symbol_struct, optional: Array of symbols.
!                    Set symbol(i)%type < 0 to suppress drawing of the i^th symbol.
!   test(:)     -- Character(*), optional: Array of text lines.
!   draw_line   -- Logical, optional: Draw lines? Default is True.
!   draw_symbol -- Logical, optional: Draw symbols? Default is True.
!   draw_text   -- Logical, optional: Draw text? Default is True.
!-

subroutine qp_draw_curve_legend (origin, text_offset, line_length, &
                    line, symbol, text, draw_line, draw_symbol, draw_text)

implicit none

type (qp_line_struct), optional :: line(:)
type (qp_symbol_struct), optional :: symbol(:)
type (qp_point_struct) origin

real(rp), optional :: line_length, text_offset
real(rp) height, xc, yc, yc2, line_len, dummy, text_off

integer i, n_rows

character(*), optional :: text(:)
character(20) :: r_name = 'qp_draw_curve_legend'

logical, optional :: draw_line, draw_symbol, draw_text
logical has_line, has_symbol, has_text

!

call qp_save_state (.true.)

call qp_set_text_attrib ('LEGEND')
call qp_set_line_attrib ('LEGEND')
height = qp_text_height_to_inches(qp_com%this_text%height) 

call qp_to_inch_abs (origin%x, origin%y, xc, yc, origin%units)

! Find out how many rows to draw

n_rows = 0
has_text = .false.
if (present(text) .and. logic_option(.true., draw_text)) then
  has_text = .true.
  n_rows = size(text)
  call qp_to_inch_rel (real_option(10.0_rp, text_offset), 0.0_rp, &
                                                       text_off, dummy, 'POINTS')
endif

has_line = .false.
line_len = 0
if (present(line) .and. logic_option(.true., draw_line)) then
  if (n_rows /= 0 .and. size(line) /= n_rows) then
    call out_io (s_error$, r_name, 'LINE ARRAY SIZE NOT THE SAME AS THE OTHERS.')
    return
  endif
  has_line = .true.
  n_rows = size(line)
  call qp_to_inch_rel (real_option(72.0_rp, line_length), 0.0_rp, &
                                                       line_len, dummy, 'POINTS')
endif

has_symbol = .false.
if (present(symbol) .and. logic_option(.true., draw_symbol)) then
  if (n_rows /= 0 .and. size(symbol) /= n_rows) then
    call out_io (s_error$, r_name, 'SYMBOL ARRAY SIZE NOT THE SAME AS THE OTHERS.')
    return
  endif
  has_symbol = .true.
  n_rows = size(symbol)
endif

! Draw the rows

yc = yc - 1.1 * height

do i = 1, n_rows
  yc2 = yc + 0.5 * height

  if (has_line .and. line(i)%width > -1) then
    call qp_set_line ('STD', line(i))
    call qp_draw_line (xc, xc + line_len, yc2, yc2, 'INCH/PAGE/LB')
  endif

  if (has_symbol .and. symbol(i)%type > -1) then
    call qp_set_symbol (symbol(i))
    call qp_draw_symbol (xc + line_len/2, yc2, 'INCH/PAGE/LB')
  endif

  if (has_text) then
    call qp_set_text_attrib ('LEGEND')
    call qp_draw_text_no_set (text(i), xc + line_len + text_off, yc, 'INCH/PAGE/LB')
  endif

  yc = yc - 1.5 * height
enddo

call qp_restore_state

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_ellipse (x0, y0, r_x, r_y, theta_xy, angle0, 
!                              del_angle, units, width, color, style, clip)
!
! Subroutine to plot a section of an ellipse.
! Drawn is:
!     (x, y) = (x0, y0) + M * (r_x * cos(theta), r_y * sin(theta))
! Where M is a rotation matrix with rotation angle theta_xy, 
! and theta goes from angle0 to angle0 + del_angle.
!
! Note: Currently this routine can only draw solid lines.
!
! Input:
!   x0, y0    -- Real(rp): Center of arc.
!   r_x       -- Real(rp): Horizontal radius.
!   r_y       -- Real(rp): Vertical radius.
!   theta_xy  -- Real(rp), optional: Ellipse tilt. Default = 0.
!   angle0    -- Real(rp), optional: Starting angle of arc to draw in radians.
!                 Default is 0
!   del_angle -- Real(rp), optional: Angle of arc to draw. 
!                 Default is 2pi.
!   units     -- Character(*), optional: Units of x, y.
!                 Default is: 'DATA/GRAPH/LB'
!                 See quick_plot writeup for more details.
!   width     -- Integer, optional: Width of line
!   color     -- Integer, optional: Line color.
!   style     -- Integer, optional: Line style. 
!                 Currently can only be 1 (solid line).
!   clip      -- Logical, optional: Clip at graph boundary?
!-

subroutine qp_draw_ellipse (x0, y0, r_x, r_y, theta_xy, angle0, &
                          del_angle, units, width, color, style, clip)

implicit none

real(rp) x0, y0, r_x, r_y
real(rp), optional :: theta_xy, angle0, del_angle
real(rp) x(1000), y(1000), ang22, del, ang, cos_xy, sin_xy
real(rp) xx0, yy0, rr_x, rr_y, ang0, del_ang, t_xy, dx, dy

integer, optional :: width, color, style
integer i

character(*), optional :: units

logical, optional :: clip

!

if (present(style)) i = style   ! so compiler will not complain

call qp_save_state (.true.)

call qp_to_inch_abs (x0, y0, xx0, yy0, units)
call qp_to_inch_rel (r_x, r_y, rr_x, rr_y, units)

! adjust angle if ang2 < ang0

t_xy = real_option(0.0_rp, theta_xy)
sin_xy = sin(t_xy)
cos_xy = cos(t_xy)

ang0 = real_option(0.0_rp, angle0)
del_ang = real_option(twopi, del_angle)

! This gives about a 0.03" line segment length

del = max(twopi/1000, twopi / (100 * (rr_x + rr_y)))
if (del_ang < 0) del = -del

! draw

ang = 0

do i = 1, size(x)
  dx = rr_x * cos(ang0+ang)
  dy = rr_y * sin(ang0+ang)
  x(i) = xx0 + dx * cos_xy - dy * sin_xy
  y(i) = yy0 + dx * sin_xy + dy * cos_xy
  if (abs(ang) >= abs(del_ang)) exit
  ang = ang + del
  if (abs(ang) >= abs(del_ang)) ang = del_ang
enddo

call qp_set_line_attrib ('STD', width, color, 1, clip)  ! solid line
call qp_draw_polyline_no_set (x(1:i), y(1:i), 'INCH/PAGE/LB')

call qp_restore_state

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_circle (x0, y0, r, angle0, del_angle, 
!                                    units, width, color, style, clip)
!
! Subroutine to plot a section of an ellipse.
! Drawn is:
!     (x, y) = (x0, y0) + M * (r * cos(theta), r * sin(theta))
! Where M is a rotation matrix with rotation angle theta_xy, 
! and theta goes from angle0 to angle0 + del_angle.
!
! Note: Currently this routine can only draw solid lines.
!
! Input:
!   x0, y0    -- Real(rp): Center of arc.
!   r         -- Real(rp): Radius.
!   angle0    -- Real(rp), optional: Starting angle of arc to draw in radians.
!                 Default is 0
!   del_angle -- Real(rp), optional: Angle of arc to draw. 
!                 Default is 2pi.
!   units     -- Character(*), optional: Units of x, y.
!                 Default is: 'DATA/GRAPH/LB'
!                 See quick_plot writeup for more details.
!   width     -- Integer, optional: Width of line
!   color     -- Integer, optional: Line color.
!   style     -- Integer, optional: Line style. 
!                 Currently can only be 1 (solid line).
!   clip      -- Logical, optional: Clip at graph boundary?
!-

subroutine qp_draw_circle (x0, y0, r, angle0, &
                        del_angle, units, width, color, style, clip)

implicit none

real(rp) x0, y0, r
real(rp), optional :: angle0, del_angle

integer, optional :: width, color, style

character(*), optional :: units

logical, optional :: clip

!

call qp_draw_ellipse (x0, y0, r, r, 0.0_rp, angle0, del_angle, &
                                        units, width, color, style, clip)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_polyline (x, y, units, width, color, style, clip)
!
! Subroutine to draw a polyline.
!
! Input:
!   x(:), y(:) -- Real(rp): (x, y) points for the polyline.
!   units      -- Character(*), optional: Units of x, y.
!                   Default is: 'DATA/GRAPH/LB'
!                   See quick_plot writeup for more details.
!   width      -- Integer, optional: Width of line
!   color      -- Integer, optional: Line color.
!   style      -- Integer, optional: Line style. 
!   clip       -- Logical, optional: Clip at graph boundary?
!-

subroutine qp_draw_polyline (x, y, units, width, color, style, clip)

implicit none

real(rp) :: x(:), y(:)
real(rp) :: xd(size(x)), yd(size(y))

integer, optional :: width, color, style

character(*), optional :: units

logical, optional :: clip

!

call qp_save_state (.true.)
call qp_set_line_attrib ('STD', width, color, style, clip)
call qp_to_inches_abs (x, y, xd, yd, units)
call qp_draw_polyline_basic (xd, yd)
call qp_restore_state

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_polyline_no_set (x, y, units)
!
! Subroutine to draw a polyline.
! This is similar to qp_draw_polyline except qp_set_line_attrib is not called.
!
! Input:
!   x(:), y(:) -- Real(rp): (x, y) points for the polyline.
!   units      -- Character(*), optional: Units of x, y.
!                   Default is: 'DATA/GRAPH/LB'
!                   See quick_plot writeup for more details.
!-

subroutine qp_draw_polyline_no_set (x, y, units)

implicit none

real(rp) :: x(:), y(:)
real(rp) :: xd(size(x)), yd(size(y))
character(*), optional :: units

!

call qp_to_inches_abs (x, y, xd, yd, units)
call qp_draw_polyline_basic (xd, yd)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_line (x1, x2, y1, y2, units, width, color, style, clip)
!
! Subroutine to draw a line.
!
! Input:
!   x1, x2     -- Real(rp): X-coords of line endpoints
!   y1, y2     -- Real(rp): Y-coords of line endpoints
!   units      -- Character(*), optional: Units of x, y.
!                   Default is: 'DATA/GRAPH/LB'
!                   See quick_plot writeup for more details.
!   width      -- Integer, optional: Width of line
!   color      -- Integer, optional: Line color.
!   style      -- Integer, optional: Line style. 
!   clip       -- Logical, optional: Clip at graph boundary?
!-

subroutine qp_draw_line (x1, x2, y1, y2, units, width, color, style, clip)

implicit none

real(rp) :: x1, x2, y1, y2

integer, optional :: width, color, style

character(*), optional :: units

logical, optional :: clip

!

call qp_save_state (.true.)
call qp_set_line_attrib ('STD', width, color, style, clip)
call qp_draw_polyline_no_set ((/ x1, x2 /), (/ y1, y2 /), units)
call qp_restore_state

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+                    
! Subroutine qp_open_page (page_type, i_chan, x_len, y_len, units, plot_file, scale)
!
! Subroutine to Initialize a page (window) for plotting.
!
! Use qp_close_page to finish a page.
! Use qp_select_page to select a page for input.
!
! if x_len and y_len are not supplied, x_len and y_len are set to the values used on a 
! previous call to qp_open_page. 
! If there has been no previous call, a size of 8-1/2 x 11 is used.
!
! The scale argument can be used to expand or shrink the drawing. 
! A value of 0.0, cannot be used with a page_type of 'X'.
! will auto-scale so that the graph of size (x_len, y_len) will fit the 8-1/2 x 11
! printed page.
!
! Note: pgplot does not do a good job with gif files. Consider making
! a postscript file and using the unix command pstogif to convert it.
!
! Input:
!   page_type -- Character(*). Device name for the type of plot.
!                 TYPE is passed to GG_SETUP. E.g.
!                 TYPE = 'X'     --> Open an X-window.
!                 TYPE = 'GIF'   --> To create a gif file.
!                 TYPE = 'GIF-L' --> Gif w/ landscape page orientation.
!                 TYPE = 'PS'    --> To create a Color PostScript file.
!                 TYPE = 'PS-L'  --> PostScript w/ landscape page orientation.
!   x_len     -- Real(rp), optional: Page horizontal width before scaling. Default: See above.
!   y_len     -- Real(rp), optional: Page vertical width before scaling. Default: See above.
!   units     -- Character(*), optional: units for x_len and y_len.
!                    'MM'     - milimeters.
!                    'INCH'   - Inches (default unless default_set_units is changed).
!                    'POINTS' - Point
!   plot_file -- Character(*), optional: Name for the plot file.
!                    Default is: 'quick_plot.ps' or 'quick_plot.gif'
!   scale     -- real(rp), optional: If positive then the entire plot will be scaled.
!                    Default is 1.0. A value of 0.0 will autoscale to fit the page.
!
! Output:
!   i_chan    -- Inteter, optional: Plot channel. 
!                 Like a unit number for a fortran OPEN.
!                 To be used with qp_select_page.
!-

subroutine qp_open_page (page_type, i_chan, x_len, y_len, units, plot_file, scale)

implicit none

real(rp), optional :: x_len, y_len, scale
real(rp) x_inch, y_inch, x_page, y_page, page_scale, x, y, l_long, l_short
real(rp), save :: x_inch_old, y_inch_old

integer, optional :: i_chan
integer ix

character(*) page_type
character(*), optional :: units, plot_file
character(16) :: r_name = 'qp_open_page'

logical saved_state

!

call qp_save_state (.false.)

qp_com%dflt_units = dflt_set$
qp_com%page_type = page_type

! set the name for the output plot file.

if (qp_com%page_type(1:3) == 'GIF') then
  qp_com%plot_file = 'quick_plot.gif'
else
  qp_com%plot_file = 'quick_plot.ps'
endif

if (present(plot_file)) then
  if (plot_file /= ' ') then
    qp_com%plot_file = plot_file
    call string_trim (qp_com%plot_file, qp_com%plot_file, ix)
  endif
endif

! set page size

l_long = 10.5
l_short = 7.8

if (present(y_len)) then
  call qp_to_inch_rel (x_len, y_len, x_inch, y_inch, units)
else
  call qp_to_inch_rel (1.0_rp, 1.0_rp, x_inch, y_inch, '%PAGE') ! Get present size.
endif

if (x_inch == 0) then
  if (qp_com%page_type == 'PS-L' .or. qp_com%page_type == 'GIF-L') then
    x_inch = l_long
    y_inch = l_short
  else
    x_inch = l_short
    y_inch = l_long
  endif
endif

qp_com%page%x1 = 0
qp_com%page%y1 = 0
qp_com%page%x2 = x_inch
qp_com%page%y2 = y_inch

! Calculate the page_scale

page_scale = real_option(1.0_rp, scale)

if (page_scale == 0) then ! Auto scale to fit standard paper size.

  if (qp_com%page_type == 'PS' .or. qp_com%page_type == 'GIF') then
    if (x_inch/y_inch > l_short/l_long) then  ! x_inch is larger
      page_scale = l_short / x_inch
    else
      page_scale = l_long / y_inch
    endif
  elseif (qp_com%page_type == 'PS-L' .or. qp_com%page_type == 'GIF-L') then
    if (x_inch/y_inch > l_long/l_short) then  ! x_inch is larger
      page_scale = l_long / x_inch
    else
      page_scale = l_short / y_inch
    endif
  else
    call out_io (s_abort$, r_name, &
                  'SCALE = 0 CAN NOT BE USE WITH PAGE_TYPE = ' // qp_com%page_type)
    call err_exit
  endif

  call qp_open_page_basic (qp_com%page_type, 0.0_rp, 0.0_rp, &
                        qp_com%plot_file, x_page, y_page, i_chan, page_scale)

else
    call qp_open_page_basic (qp_com%page_type, x_inch * page_scale, y_inch * page_scale, &
                        qp_com%plot_file, x_page, y_page, i_chan, page_scale)
endif

! set the graph parameters

call qp_set_box (1, 1, 1, 1)
qp_com%dflt_units = dflt_draw$

! set a white background for gif

if (qp_com%page_type(1:3) == 'GIF') call qp_clear_page

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_select_page (iw)
!
! Subroutine to switch to a particular page for drawing graphics.
!
! Input:
!   iw -- Integer: ID of page obtained from qp_open_page
!-

subroutine qp_select_page (iw)

implicit none

integer iw

!

call qp_select_page_basic (iw)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_close_page
!
! Subroutine to finish plotting on a page.
! For X this closes the window.
! You will need to call qp_open_page to do further graphics.
!-

subroutine qp_close_page

implicit none

character(16) :: r_name = 'qp_close_page'
integer i

!

call qp_restore_state

if (ix_qp_com > 1) then
  call out_io (s_error$, r_name, &
    'THERE HAVE BEEN MORE CALLS TO QP_SAVE_STATE THAN CALLS TO QP_RESTORE_STATE!')
  do i = 1, ix_qp_com 
    call qp_restore_state
  enddo
endif

call qp_close_page_basic

if (qp_com%page_type(1:2) == 'PS' .or. qp_com%page_type(1:3) == 'GIF') then

  call out_io (s_info$, r_name, 'Written: ' // qp_com%plot_file)

  if (qp_com%page_type(3:3) == '/') then
    call out_io (s_error$, r_name, 'TO_PRINTER SPAWN COMMAND DISABLED.')
!       call lib$spawn ('@com:ccwplot ' // trim(qp_com%plot_file) // '/noflag')
  endif

endif

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_text (text, x, y, units, justify, height 
!         color, angle, background, uniform_spacing, spacing_factor)
!
! Subroutine to draw text.
!
! Input:
!   text       -- Character(*): Character: Text to be displayed.
!   x, y       -- Real(rp): position of the text.
!   units      -- Character(*), optional: Units of x and y.
!                 Default is: 'DATA/GRAPH/LB'
!                 See quick_plot writeup for more details.
!   justify    -- Character(*), optional: Horizontal/vertical justification.
!                   Default is 'LB' (Left Bottom).
!   height     -- Real(rp), optional: height in points.
!   color      -- Integer, optional: Color index for the box
!   angle      -- Real(rp), optional: Angle to the horizontal (in degrees). 
!                   Positive angle is CCW.
!   background -- Integer, optional: Background color.
!   uniform_spacing -- Logical, optional: If T then the distance between 
!                      characters is uniform.
!   spacing_factor  -- Real(rp), optional: Spacing factor if uniform_spacing
!                       is used. 
!-

subroutine qp_draw_text (text, x, y, units, justify, height, color, &
                 angle, background, uniform_spacing, spacing_factor)

implicit none

real(rp) x, y
real(rp), optional :: angle, height, spacing_factor

integer, optional :: color, background

character(*) text
character(*), optional :: units, justify

logical, optional :: uniform_spacing

!

call qp_save_state (.true.)

call qp_set_text_attrib ('TEXT', height, color, background, uniform_spacing)
if (present(spacing_factor)) qp_com%text_spacing_factor = spacing_factor

call qp_draw_text_no_set (text, x, y, units, justify, angle)

call qp_restore_state

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_text_no_set (text, x, y, units, justify, angle)
!
! Subroutine to display on a plot a character string.
! See also: qp_draw_text.
!
! Input:
!   text       -- Character(*): Character: Text to be displayed.
!   x, y       -- Real(rp): position of the text.
!   units      -- Character(*), optional: Units of x and y.
!                 Default is: 'DATA/GRAPH/LB'
!                 See quick_plot writeup for more details.
!   justify    -- Character(*), optional: Horizontal/Vertical justification.
!                   Default is 'LB' (Left Bottom).
!   angle      -- Real(rp), optional: Angle to the horizontal (in degrees). 
!                   Positive angle is CCW.
!   uniform_spacing -- Logical, optional: If T then the distance between 
!                      characters is uniform.
!-

subroutine qp_draw_text_no_set (text, x, y, units, justify, angle)

implicit none

real(rp) x, y
real(rp), optional :: angle
real(rp) ang, ang_rad, x_inch, y_inch, dx, dy, x1, y1, h

integer ixx, i

character(*) text
character(*), optional :: units, justify
character(16) :: r_name = 'qp_draw_text_no_set'

!

if (len_trim(text) == 0) return

call qp_to_inch_abs (x, y, x1, y1, units)

ang = 0
if (present(angle)) ang = angle
ang_rad = pi * ang / 180

if (present(justify)) then
  h = qp_text_height_to_inches(qp_com%this_text%height) 
  dx = -h * sin(ang_rad) 
  dy =  h * cos(ang_rad)
  if (justify(2:2) == 'C') then
    x1 = x1 - dx / 2
    y1 = y1 - dy / 2
  elseif (justify(2:2) == 'T') then
    x1 = x1 - dx 
    y1 = y1 - dy 
  elseif (justify(2:2) /= 'B') then
    call out_io (s_fatal$, r_name, 'UNKNOWN "JUSTIFY": ' // justify)
    call err_exit
  endif
endif

if (qp_com%this_text%uniform_spacing) then
  ixx = len(trim(text)) + 1
  h = qp_text_height_to_inches(qp_com%this_text%height)
  dx = h * cos(ang_rad) 
  dy = h * sin(ang_rad) 
  x1 = x1 - dx * (ixx * qp_justify(justify) - 0.5)
  y1 = y1 - dy * (ixx * qp_justify(justify) - 0.5)
  do i = 1, len(trim(text))
    call qp_to_inch_abs (x1+i*dx, y1+i*dy, x_inch, y_inch, 'INCH/PAGE/LB')
    call qp_draw_text_basic (text(i:i), 1, x_inch, y_inch, ang, 0.5_rp)
  enddo
else
  call qp_to_inch_abs (x1, y1, x_inch, y_inch, 'INCH/PAGE/LB')
  call qp_draw_text_basic (text, len_trim(text), &
                                      x_inch, y_inch, ang, qp_justify(justify))
endif

end subroutine         

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function qp_justify (justify)
!
! Function to convert a justify character string to a real value
! representing the horizontal justification. 
!
! Input:
!   justify(1:1) -- Character(*), optional: Possibilities are:
!                   'L'  Left (Default)
!                   'C'  Center
!                   'R'  Right
!
! Output:
!   qp_justify -- Real(rp): Between 0.0, 0.5, or 1.0.
!-

function qp_justify (justify) result (horiz_justy)

implicit none

real(rp) horiz_justy
character(*), optional :: justify
character(16) :: r_name = 'qp_justify'

!

horiz_justy = 0.0

if (present(justify)) then
  if (justify(1:1) == 'C') then
    horiz_justy = 0.5
  elseif (justify(1:1) == 'R') then
    horiz_justy = 1.0
  elseif (justify(1:1) /= 'L') then
    call out_io (s_fatal$, r_name, 'BAD "JUSTIFY": ' // justify)
    call err_exit
  endif
endif

end function

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_main_title (lines, justify)
!
! Subroutine to plot the main title at the top of the page.
!
! Input:
!   lines(:) -- Character(*): Array of lines to print
!   justify  -- Character(*), optional: Horizontal justification 
!                    Vertical justification is ignored.
!                   'L' is default.
!-

subroutine qp_draw_main_title (lines, justify)

implicit none

real(rp) rx
real(rp) xt, yt, dy, x_inch, y_inch

integer i

character(*) lines(:)
character(*), optional :: justify

!

call qp_save_state (.true.)

call qp_set_text_attrib ('MAIN_TITLE')

rx = 0.1 + 0.8 * qp_justify(justify)
call qp_to_inch_abs (rx, 1.0_rp, xt, yt, '%PAGE')

dy = qp_text_height_to_inches(qp_com%this_text%height)

do i = 1, size(lines)
  yt = yt - 1.5 * dy
  call qp_to_inch_abs (xt, yt, x_inch, y_inch, 'INCH/PAGE/LB')
  call qp_draw_text_basic (lines(i), len_trim(lines(i)), &
                             x_inch, y_inch, 0.0_rp, qp_justify(justify))
enddo

call qp_restore_state

end subroutine         

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+  
! Subroutine qp_set_symbol (symbol)
!
! Subroutine to set the type and size of the symbols used in plotting data.
! See the quick_plot documentation for more details.
!
! Input:
!   symbol -- qp_symbol_struct:
!     %type          -- Integer: Symbol type. 
!     %height        -- Real(rp): Size of the symbol.
!     %color         -- Integer: Symbol color.
!     %fill_pattern  -- Integer: fill pattern.
!     %line_width    -- Integer: Line width.
!-

subroutine qp_set_symbol (symbol)

implicit none

type (qp_symbol_struct) symbol

qp_com%symbol = symbol

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+  
! Subroutine qp_set_symbol_attrib (type, height, color, fill_pattern, line_width, clip)
!
! Subroutine to set the type and size of the symbols used in plotting data.
! See the quick_plot documentation for more details.
!
! Input:
!   type         -- Integer, optional: Symbol type. 
!   height       -- Real(rp), optional: Size of the symbol.
!   color        -- Integer, optional: Symbol color.
!   fill_pattern -- Integer, optional: fill pattern.
!   line_width   -- Integer, optional: Line width.
!   clip         -- Logical, optional: Clip at graph boundary?
!                     Note: This sets the line clip also.
!-

subroutine qp_set_symbol_attrib (type, height, color, fill_pattern, line_width, clip)

implicit none

integer, optional :: type, line_width, fill_pattern, color
real(rp), optional :: height
logical, optional :: clip

!

if (present(type)) then
  qp_com%symbol%type = type
endif

if (present(height)) then
  qp_com%symbol%height = height
endif
 
if (present(color)) then
  qp_com%symbol%color = color
endif

if (present(fill_pattern)) then
  qp_com%symbol%fill_pattern = fill_pattern
endif

if (present(line_width)) then
  qp_com%symbol%line_width = line_width
endif

call qp_set_clip (clip)
call qp_set_symbol_size_basic (qp_com%symbol%height, &
         qp_com%symbol%type, qp_com%page_type, qp_com%uniform_symbol_size)
call qp_set_color_basic (qp_com%symbol%color, qp_com%page_type)
call qp_set_symbol_fill_basic (qp_com%symbol%fill_pattern)       
call qp_set_line_width_basic (qp_com%symbol%line_width)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+  
! Subroutine qp_get_symbol (symbol)
!
! Subroutine to get the symbol parameters used in plotting data.
! Use qp_set_symbol or qp_set_symbol_attrib to set symbol attributes.
! See the quick_plot documentation for more details.
!
! Output:
!   symbol -- qp_symbol_struct:
!     %type          -- Integer: Symbol type. 
!     %height        -- Real(rp): Size of the symbol.
!     %color         -- Integer: Symbol color.
!     %fill_pattern  -- Integer: fill pattern.
!     %line_width    -- Integer: Line width.
!-

subroutine qp_get_symbol (symbol)

implicit none

type (qp_symbol_struct) symbol

symbol = qp_com%symbol

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_line (who, line)
!
! Subroutine to set the default line attributes.
! See the quick_plot documentation for more details.
!
! Input:
!   who     -- Character(*): 
!                 'AXIS'     -- Graph axis.
!                 'GRID'     -- Graph grid.
!                 'LEGEND'   -- Line legend.
!                 'PLOT'     -- Plot data lines.
!                 'STD'      -- Everything else.
!   line    -- qp_line_struct: Attributes of a line
!     %style   -- Integer: Line style.
!     %width   -- Integer: Size of the line.
!     %color   -- Integer: Line color.
!-

subroutine qp_set_line (who, line)

implicit none

type (qp_line_struct) line
character(*) who
character(16) :: r_name = 'qp_set_line'

!

if (who == 'STD') then
  qp_com%std_line = line
elseif (who == 'GRID') then
  qp_com%grid_line = line
elseif (who == 'PLOT') then
  qp_com%plot_line = line
elseif (who == 'AXIS') then
  qp_com%axis_line = line
elseif (who == 'LEGEND') then
  qp_com%legend_line = line
else
  call out_io (s_fatal$, r_name, 'UNKNOWN LINE "WHO": ' // who)
  call err_exit
endif

call qp_set_line_attrib (who)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_get_line (who, line)
!
! Subroutine to get the default line attributes.
! See the quick_plot documentation for more details.
!
! Input:
!   who     -- Character(*): 
!                 'AXIS'     -- Graph axis.
!                 'GRID'     -- Graph grid.
!                 'LEGEND'   -- Line legend.
!                 'PLOT'     -- Plot data lines.
!                 'STD'      -- Everything else.
!
! Output:
!   line    -- qp_line_struct: Attributes of a line
!     %style   -- Integer: Line style.
!     %width   -- Integer: Size of the line.
!     %color   -- Integer: Line color.
!-

subroutine qp_get_line (who, line)

implicit none

type (qp_line_struct) line
character(*) who
character(16) :: r_name = 'qp_get_line'

!

if (who == 'STD') then
  line = qp_com%std_line
elseif (who == 'GRID') then
  line = qp_com%grid_line
elseif (who == 'PLOT') then
  line = qp_com%plot_line
elseif (who == 'AXIS') then
  line = qp_com%axis_line
elseif (who == 'LEGEND') then
  line = qp_com%legend_line
else
  call out_io (s_fatal$, r_name, 'UNKNOWN LINE "WHO": ' // who)
  call err_exit
endif

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_line_attrib (who, width, color, style, clip)
!        
! Subroutine to set the default line attributes.
! See the quick_plot documentation for more details.
!
! Input:
!   who     -- Character(*): 
!                 'PLOT'     -- Plot data lines.
!                 'GRID'     -- Graph grid.
!                 'AXIS'     -- Graph axis.
!                 'LEGEND    -- Line legend.
!                 'STD'      -- Everything else.
!   style   -- Integer, optional: Line style.
!   width   -- Integer, optional: Size of the line.
!   color   -- Integer, optional: Line color.
!   clip    -- Logical, optional: Clip at graph boundary?
!               Note: This sets the symbol clip also.
!-

subroutine qp_set_line_attrib (who, width, color, style, clip)

implicit none

type (qp_line_struct), pointer :: this

integer, optional :: style, width, color
logical, optional :: clip

character(*) who
character(16) :: r_name = 'qp_set_line_attrib'


!

if (who == 'STD') then
  this => qp_com%std_line
elseif (who == 'GRID') then
  this => qp_com%grid_line
elseif (who == 'PLOT') then
  this => qp_com%plot_line
elseif (who == 'AXIS') then
  this => qp_com%axis_line
elseif (who == 'LEGEND') then
  this => qp_com%legend_line
else
  call out_io (s_fatal$, r_name, 'UNKNOWN LINE "WHO": ' // who)
  call err_exit
endif

if (present(width)) this%width = width
if (present(color)) this%color = color
if (present(style)) this%style = style

call qp_set_clip (clip)
call qp_set_color_basic (this%color, qp_com%page_type)
call qp_set_line_width_basic (this%width)
call qp_set_line_style_basic (this%style)       

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_clip (clip)
!
! Subroutine to set the default clipping state.
! Note: This affects both lines and symbols.
!
! Input:
!   clip    -- Logical, optional: Clip at graph boundary?
!-

subroutine qp_set_clip (clip) 

implicit none

logical, optional :: clip

if (present(clip)) qp_com%clip = clip
call qp_set_clip_basic (qp_com%clip)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_subset_box (ix, iy, ix_tot, iy_tot, x_marg, y_marg)
!
! Subroutine to set the box for a graph. This is the same as
! QP_SET_BOX but the boundries of the page are taken to be the box boundries
! set by QP_SET_BOX. The margins between the graphs and the edges of the "page"
! are what is set by QP_SET_MARGIN. The distance between graphs is given by the
! arguments X_MARG, Y_MARG. SUBGRAPH is used to cluster graphs together.
!
! Input:
!     ix_tot, iy_tot -- Integer: X and Y box divisions of the "page".
!     ix, iy         -- Integer: Index for box to be used.
!     x_marg, y_marg -- Real(rp): Margins between boxes
!-

subroutine qp_subset_box (ix, iy, ix_tot, iy_tot, x_marg, y_marg)

implicit none

integer ix, iy, ix_tot, iy_tot

real(rp) x_marg, y_marg
real(rp) x_width, y_width

! x_width and y_width are the size of the graphs

x_width = (qp_com%box%x2 - qp_com%box%x1 - qp_com%margin%x1 - &
                qp_com%margin%x2 - (ix_tot - 1) * x_marg) / ix_tot
y_width = (qp_com%box%y2 - qp_com%box%y1 - qp_com%margin%y1 - &
                qp_com%margin%y2 - (iy_tot - 1) * y_marg) / iy_tot

qp_com%graph%x1 = qp_com%box%x1 + qp_com%margin%x1 + &
                                               (ix - 1) * (x_marg + x_width)
qp_com%graph%x2 = qp_com%graph%x1 + x_width
qp_com%graph%y1 = qp_com%box%y1 + qp_com%margin%y1 + &
                                               (iy - 1) * (y_marg + y_width)
qp_com%graph%y2 = qp_com%graph%y1 + y_width

! set QP_COM%SUBGRAPH_ON so SET_GRAPH_LIMITS will not 
! recalculate the graph boundry

qp_com%subgraph_on = .true.   
call qp_set_graph_limits 

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_text_attrib (who, height, color, background, &
!                                       uniform_spacing, spacing_factor)
!
! Routine to set the default text attributes for titles, legends, etc.
! This routine also sets the current text attributes used in drawing.
!
! See the quick_plot documentation for more details.
! Note: The background color index is common to all types of text.
!   If you set it you set it for all text.
!
! Input:                                        
!   who -- Character(*):     Used by:             Comment
!            "TEXT"          qp_draw_text         General text.
!            "MAIN_TITLE"    qp_draw_main_title   Title at top of page.
!            "GRAPH_TITLE"   qp_draw_graph_title  Title above a graph.
!            "LEGEND"        qp_draw_text_legend  Legend.
!            "LEGEND"        qp_draw_label_legend Legend.
!            "AXIS_NUMBERS"  qp_draw_graph        Axes Numbers.
!            "AXIS_LABEL"    qp_draw_graph        Axis label.
!   height      -- Real(rp), optional: Character height.
!   color       -- Integer, optional: Color index.
!   background  -- Integer, optional: Background color index.
!   uniform_spacing -- Logical, optional: If T then the distance between 
!                      characters is uniform.
!   spacing_factor  -- Real(rp), optional: Spacing factor for the 
!             uniform_spacing option. This is set globally for all "who".
!-

subroutine qp_set_text_attrib (who, height, color, background, &
                                             uniform_spacing, spacing_factor)

implicit none

integer, optional :: color, background
real(rp), optional :: height, spacing_factor
logical, optional :: uniform_spacing
character(*) who
character(16) :: r_name = 'qp_set_text_attrib'

!

if (who == "MAIN_TITLE") then
  call qp_set_this_char_size (qp_com%main_title)
elseif (who == "GRAPH_TITLE") then
  call qp_set_this_char_size (qp_com%graph_title)
elseif (who == "LEGEND") then
  call qp_set_this_char_size (qp_com%legend)
elseif (who == "TEXT") then
  call qp_set_this_char_size (qp_com%text)
elseif (who == "AXIS_NUMBERS") then
  call qp_set_this_char_size (qp_com%axis_number)
elseif (who == "AXIS_LABEL") then
  call qp_set_this_char_size (qp_com%axis_label)
else
  call out_io (s_fatal$, r_name, 'BAD "WHO": "' // trim(who) // '"' )
  call err_exit

endif

!----------------------------------------------------------------
contains

subroutine qp_set_this_char_size (this_text)

type (qp_text_struct) this_text
real(rp) text_height

!

if (present(spacing_factor)) qp_com%text_spacing_factor = spacing_factor

if (present(height)) then
  this_text%height = height
endif

if (present(color)) then
  this_text%color = color
endif

if (present(background)) then
  qp_com%text_background = background
endif

if (present(uniform_spacing)) then
  this_text%uniform_spacing = uniform_spacing
endif

text_height = this_text%height * qp_com%text_scale

call qp_set_char_size_basic (text_height)
call qp_set_color_basic (this_text%color, qp_com%page_type)
call qp_set_text_background_color_basic (qp_com%text_background)
call qp_set_line_width_basic (1)

qp_com%this_text = this_text

end subroutine

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_x_axis (who, y_pos)
!
! Subroutine to draw a horizontal axis.
!
! Input:
!   who    -- Character(*): Which axis:
!               'X'  -- Normal x-axis.
!               'X2' -- Secondary x-axis.
!   y_pos  -- Real: Vertical position (%GRAPH): 
!                0.0 = Bottom
!                1.0 = Top
!-
                                           
subroutine qp_draw_x_axis (who, y_pos) 

implicit none
                                                                      
type (qp_axis_struct) ax1, ax2

real(rp) del, dx0, dum, x1, y1, dy, x0, y0, y_pos, dy1, dy2
real(rp) dy11, dy22, x11, dg, xl0

integer i, j, m_div, who_sign, divisions, n_draw

character(*) who
character(16) justify, str
character(16) :: r_name = 'qp_draw_x_axis'

! save state

call qp_save_state (.true.)
call qp_set_clip (.false.)     ! no clipping of axis

! ax1 and ax2 are the same except when x2_mirrors_x is True.

if (who == 'X') then
  ax2 = qp_com%plot%x
  ax1 = ax2
  who_sign = +1
elseif (who == 'X2') then
  ax2 = qp_com%plot%x2
  ax1 = ax2
  if (qp_com%plot%x2_mirrors_x) ax1 = qp_com%plot%x 
  who_sign = -1
else
  call out_io (s_error$, r_name, 'BAD AXIS NAME: ' // who)
endif

! the axis line itself

call qp_set_line_attrib ('AXIS')
call qp_set_text_attrib ('AXIS_NUMBERS')

call qp_draw_polyline_no_set ((/ 0.0_rp, 1.0_rp /), &
                                (/y_pos, y_pos /), '%GRAPH')  

del = (ax1%max - ax1%min) / ax1%major_div

if (ax1%minor_div == 0) then
  call qp_calc_minor_div (del, ax1%minor_div_max, m_div)
else
  m_div = ax1%minor_div
endif

if (ax2%tick_side == 0) then
  dy1  =  ax1%major_tick_len
  dy2  = -ax1%major_tick_len
  dy11 =  ax1%minor_tick_len
  dy22 = -ax1%minor_tick_len
else
  dy1  = 0
  dy2  = ax1%major_tick_len * ax2%tick_side * who_sign
  dy11 = 0
  dy22 = ax1%minor_tick_len * ax2%tick_side * who_sign
endif

call qp_to_inch_rel (0.0_rp, y_pos, x0, y0, '%GRAPH')
call qp_to_inch_rel (1.0_rp, 0.0_rp, dg, dum, '%GRAPH')

if (ax1%type == 'LOG') then
  divisions = nint(log10(ax1%max)) - nint(log10(ax1%min))
  xl0 = (nint(log10(ax1%min)) - log10(ax1%min)) * dg / divisions
  n_draw = max(nint(real(divisions)/ax1%major_div), 1)
else
  divisions = ax1%major_div
  xl0 = 0
  n_draw = 1
endif

dx0 = dg / divisions

if (ax1%number_side*who_sign == +1) then
  justify = 'CB'
else
  justify = 'CT'
endif

!--------------------------------------------------------------------
! axis ticks and numbers


do i = 0, divisions

  x1 = i*dx0 + xl0
  y1 = y0 + ax1%number_side * who_sign * ax1%number_offset
 
  ! major ticks
  call qp_draw_polyline_no_set ((/x1, x1 /), (/ y0+dy1, y0+dy2 /), 'INCH')

  if (.not. ax2%draw_numbers) cycle
  if (n_draw * (i/n_draw) /= i) cycle

  ! numbers
  call qp_to_axis_number_text (ax1, i, str)

  if (i == 0) then
    call qp_draw_text_no_set (str, x1, y1, 'INCH', justify)
  elseif (i == ax1%major_div) then
    call qp_draw_text_no_set (str, x1, y1, 'INCH', justify)
  else
    call qp_draw_text_no_set (str, x1, y1, 'INCH', justify)
  endif

enddo

! minor ticks

do i = 0, divisions - 1
  if (ax1%type == 'LOG') then
    do j = 2, 9
      x11 = (i + log10(real(j))) * dx0 + xl0
      call qp_draw_polyline_no_set ((/ x11, x11 /), (/ y0+dy11, y0+dy22 /), 'INCH')
    enddo
  else
    do j = 1, m_div - 1
      x11 = (i + real(j) / m_div) * dx0 + xl0
      call qp_draw_polyline_no_set ((/ x11, x11 /), (/ y0+dy11, y0+dy22 /), 'INCH')
    enddo
  endif
enddo

!--------------------------------------------------------------------
! axis label

if (ax2%draw_label) then
  call qp_to_inch_rel (0.5_rp, y_pos, x0, y0, '%GRAPH')
  dy = ax1%number_side * who_sign * (ax1%number_offset + &
      ax1%label_offset + qp_text_height_to_inches(qp_com%axis_number%height))
  call qp_set_text_attrib ('AXIS_LABEL', color = ax2%label_color)
  call qp_draw_text_no_set (ax2%label, x0, y0+dy, 'INCH', justify)
endif

!

call qp_restore_state

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_y_axis (who, x_pos)
!
! Subroutine to draw a horizontal axis.
!
! Input:
!   who    -- Character(*): Which axis:
!               'Y'  -- Normal y-axis.
!               'Y2' -- Secondary y-axis.
!   x_pos  -- Real: Horizontal position (%GRAPH): 
!               0.0 = Left
!               1.0 = Right
!-

subroutine qp_draw_y_axis (who, x_pos) 

implicit none
                                                                      
type (qp_axis_struct) ax1, ax2

real(rp) del, dum, x1, y1, dx, x0, y0, x_pos, dx1, dx2, dx11, dx22
real(rp) number_len, dy0, y11, dg, yl0

integer i, j, m_div, who_sign, number_side, divisions, n_draw

character(*) who
character(2) justify
character(16) str
character(16) :: r_name = 'qp_draw_y_axis'

! save state
  
call qp_set_clip (.false.)     ! no clipping of axis
call qp_save_state (.true.)

! ax1 and ax2 are the same except when x2_mirrors_x is True.

if (who == 'Y') then
  ax2 = qp_com%plot%y
  ax1 = ax2
  who_sign = +1
elseif (who == 'Y2') then
  ax2 = qp_com%plot%y2
  ax1 = ax2
  if (qp_com%plot%y2_mirrors_y) ax1 = qp_com%plot%y 
  who_sign = -1
else
  call out_io (s_error$, r_name, 'BAD AXIS NAME: ' // who)
endif

number_side = ax1%number_side * who_sign

! draw axis line itself

call qp_set_line_attrib ('AXIS')
call qp_set_text_attrib ('AXIS_NUMBERS')

call qp_draw_polyline_no_set ((/x_pos, x_pos/), (/0.0_rp, 1.0_rp/), '%GRAPH')

! major and minor divisions calc

del = (ax1%max - ax1%min) / ax1%major_div

if (ax1%minor_div == 0) then
  call qp_calc_minor_div (del, ax1%minor_div_max, m_div)
else
  m_div = ax1%minor_div
endif

if (ax2%tick_side == 0) then
  dx1 =  ax1%major_tick_len
  dx2 = -ax1%major_tick_len
  dx11 =  ax1%minor_tick_len
  dx22 = -ax1%minor_tick_len
else
  dx1 = 0
  dx2 = ax1%major_tick_len * ax2%tick_side * who_sign
  dx11 = 0
  dx22 = ax1%minor_tick_len * ax2%tick_side * who_sign
endif

call qp_to_inch_rel (x_pos, 0.0_rp, x0, y0, '%GRAPH')
call qp_to_inch_rel (0.0_rp, 1.0_rp, dum, dg, '%GRAPH')

if (ax1%type == 'LOG') then
  divisions = nint(log10(ax1%max)) - nint(log10(ax1%min))
  yl0 = (nint(log10(ax1%min)) - log10(ax1%min)) * dg / divisions
  n_draw = max(nint(real(divisions/ax1%major_div)), 1)
else
  divisions = ax1%major_div
  yl0 = 0
  n_draw = 1
endif

dy0 = dg / divisions


! draw axis ticks and numbers

number_len = 0   ! length in inches

do i = 0, divisions

  x1 = x0 + number_side * ax1%number_offset
  y1 = i*dy0 + yl0

  ! major ticks
  call qp_draw_polyline_no_set ((/ x0+dx1, x0+dx2 /), (/ y1, y1 /), 'INCH')

  ! numbers

  if (.not. ax2%draw_numbers) cycle
  if (n_draw * (i/n_draw) /= i) cycle

  call qp_to_axis_number_text (ax1, i, str)

  if (number_side == +1) then
    justify = 'L'
  else
    justify = 'R'
  endif

  if (i == 0) then
    justify(2:2) = 'B'
  elseif (i == ax1%major_div) then
    justify(2:2) = 'T'
  else
    justify(2:2) = 'C'
  endif

  call qp_draw_text_no_set (str, x1, y1, 'INCH', justify)

  number_len = max (number_len, qp_text_len(str))

enddo

! minor ticks

do i = 0, divisions-1
  if (ax1%type == 'LOG') then
    do j = 2, 9
      y11 = (i + log10(real(j))) * dy0 + yl0
      call qp_draw_polyline_no_set ((/ x0+dx11, x0+dx22 /), (/ y11, y11 /), 'INCH')
    enddo
  else
    do j = 1, m_div - 1
      y11 = (i + real(j) / m_div) * dy0 + yl0
      call qp_draw_polyline_no_set ((/ x0+dx11, x0+dx22 /), (/ y11, y11 /), 'INCH')
    enddo
  endif
enddo

! draw label

if (ax2%draw_label) then
  call qp_to_inch_rel (x_pos, 0.5_rp, x0, y0, '%GRAPH')
  dx = number_side * (ax1%number_offset + ax1%label_offset + number_len)
  call qp_set_text_attrib ('AXIS_LABEL', color = ax2%label_color)
  call qp_draw_text_no_set (ax2%label, x0+dx, y0, 'INCH', 'CB', &
                                       angle = -90.0_rp * number_side)
endif

!

call qp_restore_state
                   
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_to_axis_number_text (axis, ix_n, text)
!
! Subroutine to form the text string for an axis number.
!
! Input:
!   axis -- qp_axis_struct:
!   ix_n -- Integer: Index of particular number.
!
! Output:
!   text -- Character(*): Character string.
!-

subroutine qp_to_axis_number_text (axis, ix_n, text)

implicit none

type (qp_axis_struct) axis

integer i, ix_n, ie, n_char, n_log_min, n_log_max, n_log_delta
integer n_zero_crit, ix, n_log, p

real(rp) val, delta, v, effective_zero

character(*) text
character(20) fmt

logical too_small, too_large

! LOG scale

if (axis%type == 'LOG') then
  ie =  (nint(log10(axis%min)) + ix_n)
  select case (ie)
  case (-2); text = '0.01'
  case (-1); text = '0.1'
  case  (0); text = '1'
  case  (1); text = '10'
  case  (2); text = '100'
  case  (3); text = '1000'
  case default; write (text, '(a, i0)') '1E', ie
  end select
  return
endif

! Calculate output number

delta = (axis%max - axis%min) / axis%major_div
val = axis%min + ix_n * delta
n_log_delta = floor(log10(abs(delta)) + 0.0001)
effective_zero = max(abs(axis%max), abs(axis%min)) / 10.0_rp**qp_com%max_digits

! If the number is essentially zero then life is simple

if (abs(val) < effective_zero) then
  text = '0'
  return
endif

! n_char is the number of characters we need to draw

n_log_min = 1000
n_log_max = -1000
do i = 0, axis%major_div
  v = axis%min + i * delta
  if (abs(v) < effective_zero) cycle  ! Ignore zero.
  n_log = floor(log10(abs(v)) + 0.0001)
  n_log_min   = min (n_log_min, n_log)
  n_log_max   = max (n_log_max, n_log)
enddo

n_char = max(0, n_log_max) + 1
n_char = n_char + max(axis%places, 0)
if (axis%max < 0 .or. axis%min < 0) n_char = n_char + 1  ! for negative sign
if (axis%places > 0) n_char = n_char + 1  ! add 1 character for the decimal place itself

! Special case: Switch to scientific notation if the number is too big or too small.

too_large = .false.
too_small = .false.
n_zero_crit = qp_com%max_axis_zero_digits 

if (axis%places <= 0) then  ! check for too big
  if (n_log_delta > n_zero_crit) then 
    too_large = .true.
    n_char = n_char - n_log_delta
  endif

else                  ! check for too small
  if (n_log_min < -n_zero_crit) then 
    too_small = .true.
    n_char = n_char + n_log_delta
  endif
endif

! Convert output number to a string.
! Special case: switch to scientific notation if the number is too big or string
! is not long enough.

if (n_char+1 > len(text)) then
  text = "*****"

elseif (too_large) then
  p = axis%places
  v = val*10d0**p
  ! Make sure that there is no integer overflow with the write statement
  if (abs(v) > huge(1)) p = p - ceiling(log10(abs(v)/huge(1)))
  write (text, '(i0, a, i0)') nint(val*10d0**p), 'E', -p

elseif (too_small) then 
  p = axis%places + n_log_min
  if (p <= 0) then
    write (text, '(i0, a, i0)') nint(val*10d0**(axis%places)), 'E', -(axis%places)
  else
    write (fmt, '(a, i0, a, i0, a)') '(f', p+4, '.', p, ', a, i0)' 
    write (text, fmt) val*10d0**(-n_log_min), 'E', n_log_min
  endif

elseif (axis%places <= 0) then
  write (fmt, '(a, i0, a)') '(f', n_char+1, '.0)'    
  write (text, fmt) val
  call string_trim (text, text, ix)
  ix = index(text, '.')
  text = text(1:ix-1)

else
  write (fmt, '(a, i0, a, i0, a)') '(f', n_char, '.', axis%places, ')'
  write (text, fmt) val
endif

call string_trim (text, text, i)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function qp_text_len (text)
!
! Function to find the length of a text string.
!
! input:
!   text -- Character(*): Text string.
!
! Output:
!   qp_text_len -- Real(rp): Length of text in inches.
!-

function qp_text_len (text) result (t_len)

implicit none

real(rp) t_len
real tl, dum

character(*) text

!

t_len = qp_text_len_basic (text, len_trim(text))

end function

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_graph_attrib (draw_grid, draw_title)
!
! Subroutine to set attributes of the current graph
!
! Input:
!   draw_grid  -- Logical, optional: Draw a grid?
!   draw_title -- Logical, optional: Draw the title?
!-

Subroutine qp_set_graph_attrib (draw_grid, draw_title)

implicit none

logical, optional :: draw_grid, draw_title

!

if (present(draw_grid))  qp_com%plot%draw_grid  = draw_grid
if (present(draw_title)) qp_com%plot%draw_title = draw_title

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_grid
!
! Subroutine to draw a grid on the current graph.
!-

subroutine qp_draw_grid 

implicit none

type (qp_axis_struct), pointer :: axis
real(rp) r0, z(2), r01(2)
integer i, divisions

!

call qp_save_state (.true.)
call qp_set_line_attrib ('GRID')
r01 = (/ 0.0, 1.0 /)

! horizontal lines

axis => qp_com%plot%xx

if (axis%type == 'LOG') then
  divisions = nint(log10(axis%max)) - nint(log10(axis%min))
  r0 = (nint(log10(axis%min)) - log10(axis%min)) / divisions
else
  divisions = axis%major_div
  r0 = 0
endif

do i = 1, divisions - 1
  z = real(i) / divisions + r0
  call qp_draw_polyline_no_set (z, r01, '%GRAPH')
enddo

! vertical lines

axis => qp_com%plot%yy

if (axis%type == 'LOG') then
  divisions = nint(log10(axis%max)) - nint(log10(axis%min))
  r0 = (nint(log10(axis%min)) - log10(axis%min)) / divisions
else
  divisions = axis%major_div
  r0 = 0
endif

do i = 1, divisions - 1
  z = real(i) / divisions + r0
  call qp_draw_polyline_no_set (r01, z, '%GRAPH')
enddo

!

call qp_restore_state 

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_calc_minor_div (delta, div_max, divisions)
!
! Subroutine to calculate the number of minor divisions an axis should have.
! This routine picks the maximum number of minor divisions consistant with
! the restrictions that the number of divisions is "nice" and that
!               divisions <= div_max
!
! Input:
!   delta    -- Real(rp): Axis width between main divisions.
!   div_max  -- Integer: Maximum number divisions can be.
!   
! Output:
!   divisions -- Integer: Number of minor divisions.
!                   Minimum number this can be is 1.
!-

subroutine qp_calc_minor_div (delta, div_max, divisions)

implicit none

real(rp) delta
real(dp) log_del

integer div_max, divisions, idel
integer id5, id10

! scale delta so it is in the range of [10, 100)

log_del = log10 (delta * 1.000000001_dp)
idel = nint(delta / 10d0**(floor(log_del)-1))

! First look for a division that gives a width that is a multiple of 5 or 10.
! Choose id10 over id5 except when id5 >= 2 * id10
! A division of 1 is not acceptable in this first step.

do id10 = div_max, 1, -1
  if (mod(idel, 10 * id10) == 0) exit
enddo

do id5 = div_max, 1, -1
  if (mod(idel, 5 * id5) == 0) exit
enddo

if (id10 > 1 .and. id5 > 1) then
  if (id5 >= 2 * id10) then
    divisions = id5
  else
    divisions = id10
  endif
  return
endif

if (id10 > 1) then
  divisions = id10
  return
endif

if (id5 > 1) then
  divisions = id5
  return
endif


! Now look for anything that divides evenly into idel.

do divisions = div_max, 1, -1
  if (mod(idel, divisions) == 0) return
enddo

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_translate_to_color_index (name, index)
!
! Subroutine to translate from a string to a color index.
! Translation is case insensitive.
!
! Input:
!   name -- Character(*): Name of the color.
!
! Output:
!   index -- Integer: Color index. -1 => Unknown name.
!-

subroutine qp_translate_to_color_index (name, index)

implicit none

integer index

character(*) name
character(16) color, this
character(16) :: r_name = 'qp_translate_to_color_index'

!

call str_upcase (this, name)

do index = lbound(qp_color_name, 1), ubound (qp_color_name, 1)
  call str_upcase (color, qp_color_name(index))
  if (color == this) return
enddo

call out_io (s_error$, r_name, 'UNKNOWN COLOR NAME: ' // name)
index = -1

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_get_parameters (text_scaledefault_draw_units, default_set_units, &
!                           default_axis_slop_factor)
!
! Routine to get various quick_plot parameters.
! The routine qp_set_parameters can be used to set these parameters.
!
! Output:
!   text_scale         -- Real(rp), optional: Overall text scale.
!   default_draw_units -- Character(*), optional: 
!   default_set_units  -- Character(*), optional: 
!   default_axis_slop_factor -- Real(rp), optional: 
!                           See qp_calc_axis_scale for more info. 
!-

subroutine qp_get_parameters (text_scale, default_draw_units, default_set_units, &
                         default_axis_slop_factor)

implicit none

real(rp), optional :: text_scale
real(rp), optional :: default_axis_slop_factor
character(*), optional :: default_draw_units, default_set_units

!

if (present(text_scale)) text_scale = qp_com%text_scale

if (present(default_set_units)) &
        call qp_join_units_string (qp_com%dflt_set_units(1), &
            qp_com%dflt_set_units(2), qp_com%dflt_set_units(3), default_set_units)

if (present(default_draw_units)) &
        call qp_join_units_string (qp_com%dflt_draw_units(1), &
            qp_com%dflt_draw_units(2), qp_com%dflt_draw_units(3), default_draw_units)

if (present(default_axis_slop_factor)) &
         default_axis_slop_factor = qp_com%dflt_axis_slop_factor 



end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_parameters (text_scale, default_draw_units, default_set_units, 
!                            default_axis_slop_factor)
!
! Routine to set various quick_plot parameters.
! The routine qp_get_parameters can be used to obtain these parameters.
!
! Input:
!   text_scale         -- Real(rp), optional: Overall text scale (1.0 = default).
!   default_draw_units -- Character(*), optional: 
!                           Initial default is: 'DATA/GRAPH/LB'
!   default_set_units  -- Character(*), optional: 
!                           Initial default is: 'INCH/PAGE/LB'
!   default_axis_slop_factor -- Real(rp), optional: 
!                           See qp_calc_axis_scale for more info. 
!                           Initial default is: 0.001
!-

subroutine qp_set_parameters (text_scale, default_draw_units, default_set_units, &
                                                           default_axis_slop_factor)

implicit none

real(rp), optional :: text_scale
real(rp), optional :: default_axis_slop_factor
character(*), optional :: default_draw_units, default_set_units

!

if (present(text_scale)) qp_com%text_scale = text_scale

if (present(default_set_units)) then
  qp_com%dflt_units = dflt_set$
  call qp_split_units_string (qp_com%dflt_set_units(1), &
        qp_com%dflt_set_units(2), qp_com%dflt_set_units(3), default_set_units)
endif

if (present(default_draw_units)) then
  qp_com%dflt_units = dflt_draw$
  call qp_split_units_string (qp_com%dflt_draw_units(1), &
         qp_com%dflt_draw_units(2), qp_com%dflt_draw_units(3), default_draw_units)
endif

if (present(default_axis_slop_factor)) &
         qp_com%dflt_axis_slop_factor = default_axis_slop_factor

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function qp_text_height_to_inches(height_pt) result (height_inch)
!
! Function to convert from a text height in points to a text height in
! inches taking into account the text_scale.
!
! Input:
!   height_pt -- Real(rp): Height in points.
!
! Output:
!   height_inch -- Real(rp): Height in inches.
!-

function qp_text_height_to_inches(height_pt) result (height_inch)

implicit none

real(rp) height_pt, height_inch, fudge_factor

fudge_factor = 0.68  ! pgplot fudge factor!!!
height_inch = fudge_factor * height_pt * qp_com%text_scale / 72.0

end function

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_read_data (iu, err_flag, x, ix_col, y, iy_col, z, iz_col, 
!                                                                t, it_col)
!
! Subroutine to read columns of data.
! Note: This routine assumes that the file has already been opened.
! The routine will read from the current record down until data is 
! encountered and then will read until an end-of-file or a non-data line
! is encountered. A non-data line is a line whose first non-blank character
! is not a digit or a ".", "-", or "+".
!
! By repeated calls to qp_read_data you can read in multiple data sets.
! A data set can be skipped by using only the first 2 arguments:
!     call qp_read_data (iu, err_flag)
!
! Input:
!   iu       -- Integer: File unit number.
!   err_flag -- Logical: Set True if there is an error reading data.
!   ix_col   -- Integer, optional: Column number of the x data.
!   iy_col   -- Integer, optional: Column number of the y data.
!   iz_col   -- Integer, optional: Column number of the z data.
!   it_col   -- Integer, optional: Column number of the t data.
!
! Output:
!   x        -- Real(rp), allocatable, optional: Array to put data. The array 
!                  size will be changed to match the number of data points.
!   y        -- Real(rp), allocatable, optional: Another array for data.
!   z        -- Real(rp), allocatable, optional: Another array for data.
!   t        -- Real(rp), allocatable, optional: Another array for data.
!-

subroutine qp_read_data (iu, err_flag, x, ix_col, y, iy_col, z, iz_col, t, it_col)

implicit none

integer, parameter :: i_del = 100

real(rp), allocatable, optional :: x(:), y(:), z(:), t(:)
real(rp), allocatable, save :: xyz(:)
real(rp) :: x1(i_del), y1(i_del), z1(i_del), t1(i_del)

integer i, j, ix, ios
integer iu, i_size
integer, optional :: ix_col, iy_col, iz_col, it_col

logical err_flag, good_x, good_y, good_z, good_t

character(140) line, line_in
character(16) :: r_name = 'qp_read_data'

!

call skip_header (iu, err_flag)
if (err_flag) return

allocate (xyz(0))

! loop over lines

i_size = 0  ! size of x, y, and z arrays

if (present(x)) then
  if (allocated(x)) deallocate(x)
  allocate (x(0))
endif

if (present(y)) then
  if (allocated(y)) deallocate(y)
  allocate (y(0))
endif

if (present(z)) then
  if (allocated(z)) deallocate(z)
  allocate (z(0))
endif

if (present(t)) then
  if (allocated(t)) deallocate(t)
  allocate (t(0))
endif

main_loop: do
  do i = 1, i_del  
    read (iu, '(a)', iostat = ios) line_in
    if (ios /= 0) exit main_loop

    call string_trim (line_in, line, ix)
    if (index('1234567890-+.', line(1:1)) == 0) exit main_loop

    ! loop over columns

    good_x = .false.;   good_y = .false.;   
    good_z = .false.;   good_t = .false.

    do j = 1, 200

      if (present(x)) then
        if (ix_col == j) then
          read (line, *, err = 9000) x1(i)
          good_x = .true.
        endif
      endif

      if (present(y)) then
        if (iy_col == j) then
          read (line, *, err = 9000) y1(i)
          good_y = .true.
        endif
      endif

      if (present(z)) then
        if (iz_col == j) then
          read (line, *, err = 9000) z1(i)
          good_z = .true.
        endif
      endif

      if (present(t)) then
        if (it_col == j) then
          read (line, *, err = 9000) t1(i)
          good_t = .true.
        endif
      endif

      call string_trim (line(ix+1:), line, ix)
      if (ix == 0) exit

    enddo

    if (.not. good_x .and. present(x)) goto 9000
    if (.not. good_y .and. present(y)) goto 9000
    if (.not. good_z .and. present(z)) goto 9000
    if (.not. good_t .and. present(t)) goto 9000

  enddo

  ! transfer x1, y1, and z1 arrays to x, y, and z.

  if (present(x)) call load_data(x, x1)
  if (present(y)) call load_data(y, y1)
  if (present(z)) call load_data(z, z1)
  if (present(t)) call load_data(t, t1)
  i_size = i_size + i_del

enddo main_loop

! last transfer

if (present(x)) call load_data(x, x1(1:i-1))
if (present(y)) call load_data(y, y1(1:i-1))
if (present(z)) call load_data(z, z1(1:i-1))
if (present(t)) call load_data(t, t1(1:i-1))
deallocate(xyz)

return

!

9000 continue
call out_io (s_error$, r_name, 'ERROR READING DATA LINE: ' // line_in)
err_flag = .true.
deallocate(xyz)

!--------------------------------------------------------------------------
contains

subroutine load_data (t, t1)

real(rp), allocatable :: t(:)
real(rp) :: t1(:)
integer id

id = size(t1)

if (size(xyz) < i_size+id) then
  deallocate (xyz)
  allocate (xyz(i_size+id))
endif

xyz(1:i_size) = t
xyz(i_size+1:i_size+id) = t1
deallocate (t)
allocate (t(i_size+id))
t = xyz(1:i_size+id)

end subroutine

end subroutine


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_eliminate_xy_distortion
!
! This subroutine will increase the x or y margins so that the conversion
! between data units and page units is the same for the x and y axes.
!
! In other words, This routine will make sure that the following:
!   call qp_draw_line (0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 'DATA')
!   call qp_draw_line (0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 'DATA')
!   call qp_draw_line (1.0_rp, 1.0_rp, 0.0_rp, 1.0_rp, 'DATA')
!   call qp_draw_line (0.0_rp, 1.0_rp, 1.0_rp, 1.0_rp, 'DATA')
! will look like a square (and not a rectangle) when drawn.
! This routine is useful in drawing maps.
!-

subroutine qp_eliminate_xy_distortion

implicit none

real(rp) x_scale, y_scale, r, d

!

call qp_to_inch_rel (1.0_rp, 1.0_rp, x_scale, y_scale)

if (x_scale > y_scale) then  ! shrink in x
  r = (x_scale - y_scale) / x_scale 
  d = qp_com%graph%x2 - qp_com%graph%x1
  qp_com%margin%x1 = qp_com%margin%x1 + r * d / 2
  qp_com%margin%x2 = qp_com%margin%x2 + r * d / 2

else  ! shrink in y
  r = (y_scale - x_scale) / y_scale 
  d = qp_com%graph%y2 - qp_com%graph%y1
  qp_com%margin%y1 = qp_com%margin%y1 + r * d / 2
  qp_com%margin%y2 = qp_com%margin%y2 + r * d / 2

endif

call qp_set_graph_limits

end subroutine

end module
