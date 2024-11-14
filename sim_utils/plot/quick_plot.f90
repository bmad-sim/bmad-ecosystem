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
!   qp_get_axis_attrib
!   qp_get_text_atttrib
!   qp_get_layout_attrib
!   qp_get_parameters
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
! A line's appearance is specified by three quantities:
!      pattern     [solid, dashed, etc.]
!      width   
!      color
!
! These three quantities are collectively called a "style".
! Different styles can be associated with different parts of a drawing using
! the qp_set_line_attrib routine. The parts of the graph are:
!    'AXIS'     -- Graph axis.
!    'GRID'     -- Graph grid.
!    'LEGEND'   -- Line legend.
!    'PLOT'     -- Plot data lines.
!    'STD'      -- Everything else.
!
!--------------------------------------------------------------------
!
! The appearance of text is specified by the factors:
!     height
!     color
!     background color
!     spacing between characters
!
! These three quantities are collectively called a "style".
! Different styles can be associated with different parts of a drawing using
! the qp_set_text_attrib routine. The parts of the graph are:
!            Style:          used by:             Comment:
!            "MAIN_TITLE"    qp_draw_main_title   Title at top of page.
!            "GRAPH_TITLE"   qp_draw_graph_title  Title above a graph.
!            "LEGEND"        qp_draw_text_legend  Legend.
!            "LEGEND"        qp_draw_label_legend Legend.
!            "AXIS_NUMBERS"  qp_draw_graph        Axes Numbers.
!            "AXIS_LABEL"    qp_draw_graph        Axis label.
!            "TEXT"          qp_draw_text         Everything else.
!
!--------------------------------------------------------------------
!
! The line patterns are:
!     1 - solid$                  Solid
!     2 - dashed$                 Dashed
!     3 - dash_dot$               Dash-dot 
!     4 - dotted$                 Dotted
!     5 - dash_dot3$              Dash-dot-dot-dot        
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
!    13 - Reddish_Purple$
!    14 - Dark_Grey$
!    15 - Light_Grey$
!    16 - Transparent$
!
!--------------------------------------------------------------------
!
! The fill patterns are:
!     1 - solid_fill$        
!     2 - no_fill$           
!     3 - hatched$           
!     4 - cross_hatched$     
!
!--------------------------------------------------------------------
!
! The symbols are:
!     0 - square_sym$
!     1 - dot_sym$
!     2 - plus_sym$
!     3 - times_sym$
!     4 - circle_sym$
!     5 - x_sym$
!     7 - triangle_sym$
!     8 - circle_plus_sym$
!     9 - circle_dot_sym$
!    10 - square_concave_sym$
!    11 - diamond_sym$
!    12 - star5_sym$
!    13 - triangle_filled_sym$
!    14 - red_cross_sym$
!    15 - star_of_david_sym$
!    16 - square_filled_sym$
!    17 - circle_filled_sym$
!    18 - star5_filled_sym$
!
! The correspoinding name strings are are constructed by removing the "_sym$" from the parameter.
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

#if defined (CESR_PLPLOT)
  use plplot_interface
  character(*), parameter :: qp_base_library = 'PLPLOT'
#elif defined (CESR_NOPLOT)
  use noplot_interface
  character(*), parameter :: qp_base_library = 'NOPLOT'
#else
  use pgplot_interface
  character(*), parameter :: qp_base_library = 'PGPLOT'
#endif

!---------------------------------------------------------------------------
! common block
! General NOTE: qp_com is made private so that you cannot change it directly.
! This was done since the layout of qp_com can change.

type (qp_state_struct), pointer, save :: qp_com
type (qp_state_struct), target, save :: qp_save_com(0:20)

integer, save, private :: ix_qp_com = 0

private qp_save_com, qp_com

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
if (present(margin)) call qp_set_margin (margin%x1, margin%x2, margin%y1, margin%y2, margin%units)
if (present(page_border)) call qp_set_page_border (page_border%x1, &
                                  page_border%x2, page_border%y1, page_border%y2, page_border%units)

end subroutine qp_set_layout

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
character(*), parameter :: r_name = "qp_get_layout_attrib"
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
  call out_io (s_error$, r_name, 'BAD "WHO": ' // who)
  if (global_com%exit_on_error) call err_exit
endif

qp_com%dflt_units = dflt_set$

call qp_from_inch_abs (rect%x1, rect%y1, x1, y1, units)
call qp_from_inch_abs (rect%x2, rect%y2, x2, y2, units)

qp_com%dflt_units = dflt_draw$

end subroutine qp_get_layout_attrib

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_wait_to_flush (wait)
!
! Routine to signal to quick_plot whether to wait for flushing the plot buffer.
!
! If in a wait state, the plot buffer will not be flushed until this routine is called
! with wait = False. Note: By default, quick_plot is not in a wait state.
!
! Input:
!   wait    -- logical: If True, go into a wait state for flushing. 
!                       If False, flush the buffer and go into a non-wait state.
!-

subroutine qp_wait_to_flush(wait)

implicit none
logical wait

!

call qp_wait_to_flush_basic(wait)

end subroutine qp_wait_to_flush

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
qp_com%plot%x_is_active_axis = .true.
qp_com%plot%y_is_active_axis = .true.

end subroutine qp_init_com_struct

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
! the underlying plot library (pgplot or plplot). 
!
! qp_save_state always buffers the quick_plot state and the buffer_basic
!   argument determines whether the base state is buffered as well.
! The base state is the state of the underlying plot library. 
!
! Note: Buffering can make things go faster but with buffering the display 
!   will not be updated until you call qp_restore_state.
!
! Also see qp_wait_to_flush which can be used to prevent flushing of the plot buffer.
!
! Input:
!   buffer_basic -- Logical: If True then buffer the base state.
!-

subroutine qp_save_state (buffer_basic)

implicit none

logical buffer_basic

character(*), parameter :: r_name = 'qp_save_state'

!

if (buffer_basic) call qp_save_state_basic

if (ix_qp_com == size(qp_save_com)) then
  call out_io (s_error$, r_name, 'TRYING TO SAVE TOO MANY STATES!')
  if (global_com%exit_on_error) call err_exit
endif

if (ix_qp_com == 0) call qp_init_com_struct
ix_qp_com = ix_qp_com + 1

qp_save_com(ix_qp_com) = qp_com
qp_com => qp_save_com(ix_qp_com)
qp_com%buffer = buffer_basic

if (qp_com%plot%x_is_active_axis) then
  qp_com%plot%xx => qp_com%plot%x
else
  qp_com%plot%xx => qp_com%plot%x2
endif

if (qp_com%plot%y_is_active_axis) then
  qp_com%plot%yy => qp_com%plot%y
else
  qp_com%plot%yy => qp_com%plot%y2
endif

end subroutine qp_save_state

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

character(*), parameter :: r_name = 'qp_restore_state'

!

if (qp_com%buffer) call qp_restore_state_basic()

if (ix_qp_com == 0) then
  call out_io (s_error$, r_name, 'NO STATE TO RESTORE!')
  if (global_com%exit_on_error) call err_exit
endif

ix_qp_com = ix_qp_com - 1
qp_com => qp_save_com(ix_qp_com)

end subroutine qp_restore_state

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
!                 'X' Left x-axis
!                 'X2' Right x-axis
!                 'XX' Active x-axis.
!                 'Y' Bottom y-axis.
!                 'Y2' Top y-axis.
!                 'YY' Active y-axis.
!
! Output:
!   axis_ptr -- Qp_axis_struct, pointer: Pointer to the common block axis.
!-

subroutine qp_pointer_to_axis (axis_str, axis_ptr)

implicit none

type (qp_axis_struct), pointer :: axis_ptr

character(*) axis_str
character(*), parameter :: r_name = 'qp_pointer_to_axis'

!

if (axis_str == 'X') then
  axis_ptr => qp_com%plot%x
elseif (axis_str == 'X2') then
  axis_ptr => qp_com%plot%x2
elseif (axis_str == 'XX') then
  axis_ptr => qp_com%plot%xx
elseif (axis_str == 'Y') then
  axis_ptr => qp_com%plot%y
elseif (axis_str == 'Y2') then
  axis_ptr => qp_com%plot%y2
elseif (axis_str == 'YY') then
  axis_ptr => qp_com%plot%yy
else
  call out_io (s_error$, r_name, 'INVALID AXIS: ' // axis_str)
  if (global_com%exit_on_error) call err_exit
endif

end subroutine qp_pointer_to_axis

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_use_axis (x, y)
!
! Subroutine to set what axes are active (being used) : X or X2, Y or Y2.
!
! Input:
!   x -- Character(*), optional: 'X', or 'X2'
!   y -- Character(*), optional: 'Y', or 'Y2'
!-

subroutine qp_use_axis (x, y)

implicit none

character(*), optional :: x, y
character(*), parameter :: r_name = 'qp_use_axis'

!

if (present(x)) then
  select case (x)
  case ('X')
    qp_com%plot%xx => qp_com%plot%x
    qp_com%plot%x_is_active_axis = .true.
  case ('X2')
    qp_com%plot%xx => qp_com%plot%x2
    qp_com%plot%x_is_active_axis = .false.
  case default
    call out_io (s_error$, r_name, 'BAD "X": ' // x)
  end select
endif

if (present(y)) then
  select case (y)
  case ('Y')
    qp_com%plot%yy => qp_com%plot%y
    qp_com%plot%y_is_active_axis = .true.
  case ('Y2')
    qp_com%plot%yy => qp_com%plot%y2
    qp_com%plot%y_is_active_axis = .false.
  case default
    call out_io (s_error$, r_name, 'BAD "Y": ' // y)
  end select
endif

end subroutine qp_use_axis

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_axis (axis_str, a_min, a_max, div, places, label, draw_label, 
!               draw_numbers, minor_div, minor_div_max, mirror, number_offset, label_offset, 
!               major_tick_len, minor_tick_len, ax_type, tick_min, tick_max, dtick, set_ticks, axis)
!   
! Routine to set (but not plot) the min, max, divisions etc. for the X and Y axes. 
!
! Input:
!   axis_str      -- character(*): 
!                     'X' to set the Left x-axis
!                     'Y' to set the Bottom y-axis.
!                     'X2' to set the Right x-axis
!                     'Y2' to set the Top y-axis.
!   a_min         -- real(rp), optional: Axis minimum in data units.
!   a_max         -- real(rp), optional: Axis maximum in data units.
!   div           -- integer, optional: Number of major divisions.
!   places        -- integer, optional: Number of decmal places after the decimal 
!                       point. A Negative number surpresses that number of zeros.
!                       E.g. For the number 30 then 
!                            places =  2 -> Output is: "30.00"
!                            places =  0 -> Output is: "30"
!                            places = -1 -> Output can be scalled to: "3"
!   label         -- character(*), optional: Axis label.
!   draw_label    -- logical, optional: Draw axis label.
!   draw_numbers  -- logical, optional: Draw axis numbers
!   minor_div     -- integer, optional: Number of minor divisions.
!   minor_div_max -- integer, optional: Maximum number of minor divisions.
!                     This is used when you want Quick_Plot to pick the
!                     actual number of minor divisions
!   mirror        -- logical, optional: If True and axis = "x" or axis = "x2" then
!                     the x2 axis will mirror the x axis. Mirroring means that
!                     the major and minor ticks of the x2 axis will be the 
!                     same as the x axis. A similar situation holds if axis = "y"
!                     or axis = "y2".
!   number_offset  -- real(rp), optional: Offset from axis line in inches.
!   label_offset   -- real(rp), optional: Offset form numbers in inches.
!   major_tick_len -- real(rp), optional: Major tick length in inches.
!   minor_tick_len -- real(rp), optional: Minor tick length in inches.
!   ax_type        -- character(16), optional: Axis type. 'LINEAR', or 'LOG'.
!   tick_min       -- real(rp), optional: Min tick location in data units.
!   tick_max       -- real(rp), optional: Max tick location in data units.
!   dtick          -- real(rp), optional: Distance between ticks in data units.
!   set_ticks      -- logical, optional: If True, set %tick_min, %tick_max, and %dtick using values of %max, %min, and %major_div.
!                       Default is True if a_min or a_max arguments are present. Otherwise the default is False.
!   axis           -- qp_axis_struct, optional: Axis. If present with other arguments then the other arguments
!                       will override components of this argument.
!-

subroutine qp_set_axis (axis_str, a_min, a_max, div, places, label, draw_label, &
                  draw_numbers, minor_div, minor_div_max, mirror, number_offset, label_offset, &
                  major_tick_len, minor_tick_len, ax_type, tick_min, tick_max, dtick, set_ticks, axis)

implicit none

type (qp_axis_struct), pointer :: this_axis
type (qp_axis_struct), optional :: axis

real(rp), optional :: a_min, a_max, number_offset, tick_min, tick_max, dtick
real(rp), optional :: label_offset, major_tick_len, minor_tick_len

integer, optional :: div, places, minor_div, minor_div_max
logical, optional :: draw_label, draw_numbers, mirror, set_ticks
logical tick_set

character(*), optional :: label, ax_type
character(*) axis_str

!

call qp_pointer_to_axis (axis_str, this_axis)

if (present(axis))   this_axis = axis

if (present(tick_min)) this_axis%tick_min = tick_min
if (present(tick_max)) this_axis%tick_max = tick_max
if (present(dtick))    this_axis%dtick = dtick
if (present(a_min))    this_axis%min = a_min
if (present(a_max))    this_axis%max = a_max
if (present(div))      this_axis%major_div = div
if (present(places))   this_axis%places = places

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

if (present(a_min) .or. present(a_max)) then
  tick_set = logic_option(.true., set_ticks)
else
  tick_set = logic_option(.false., set_ticks)
endif

if (tick_set) THEN
  this_axis%tick_max = this_axis%max
  this_axis%tick_min = this_axis%min
  if (this_axis%type == 'LOG') then
    this_axis%dtick = 10
  else
    this_axis%dtick = (this_axis%max - this_axis%min) / max(1, this_axis%major_div)
  endif
endif

end subroutine qp_set_axis

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_get_axis_attrib (axis_str, a_min, a_max, div, places, label, 
!               draw_label, draw_numbers, minor_div, mirror, number_offset, 
!               label_offset, major_tick_len, minor_tick_len, ax_type, tick_min, tick_max, dtick)
!   
! Subroutine to get the min, max, divisions etc. for the X and Y axes.  
!
! Input:
!   axis_str       -- Character(*): 
!                      'X' to set the Left x-axis
!                      'Y' to set the Bottom y-axis.
!                      'X2' to set the Right x-axis
!                      'Y2' to set the Top y-axis.
!   a_min          -- Real(rp), optional: Axis minimum.
!   a_max          -- Real(rp), optional: Axis maximum.
!   div            -- Integer, optional: Number of major divisions.
!   places         -- Integer, optional: Number of decmal places after the decimal 
!                        point. A Negative number surpresses that number of zeros.
!                        E.g. For the number 30 then 
!                             places =  2 -> Output is: "30.00"
!                             places =  0 -> Output is: "30"
!                             places = -1 -> Output can be scalled to: "3"
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
!   tick_min       -- real(rp), optional: Min tick location in data units.
!   tick_max       -- real(rp), optional: Max tick location in data units.
!   dtick          -- real(rp), optional: Distance between ticks in data units.
!-

subroutine qp_get_axis_attrib (axis_str, a_min, a_max, div, places, label, draw_label, &
                draw_numbers, minor_div, mirror, number_offset, label_offset, &
                major_tick_len, minor_tick_len, ax_type, tick_min, tick_max, dtick)

implicit none

type (qp_axis_struct), pointer :: this_axis
real(rp), optional :: a_min, a_max, number_offset, tick_min, tick_max, dtick
real(rp), optional :: label_offset, major_tick_len, minor_tick_len

integer, optional :: div, places, minor_div
logical, optional :: draw_label, draw_numbers, mirror

character(*), optional :: label, ax_type
character(*) axis_str

!

call qp_pointer_to_axis (axis_str, this_axis)

if (present(a_min))     a_min     = this_axis%min
if (present(a_max))     a_max     = this_axis%max
if (present(div))       div       = this_axis%major_div
if (present(places))    places    = this_axis%places
if (present(tick_min))  tick_min  = this_axis%tick_min
if (present(tick_max))  tick_max  = this_axis%tick_max
if (present(dtick))     dtick     = this_axis%dtick

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

end subroutine qp_get_axis_attrib

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_calc_and_set_axis (axis_str, data_min, data_max, div_min, div_max,
!                                                      bounds, axis_type, slop_factor)
!
! Subroutine to calculate a "nice" plot scale given the minimum and maximum of the data.
!
! Input:
!   axis_str    -- character(*): 
!                   'X' to set the Left x-axis
!                   'Y' to set the Bottom y-axis.
!                   'X2' to set the Right x-axis
!                   'Y2' to set the Top y-axis.
!   data_min    -- real(rp): Minimum of the data
!   data_max    -- real(rp): Maximum of the data
!   div_min     -- integer: Minimum number of divisions.
!   div_max     -- integer: Maximum number of divisions.
!   bounds      -- character(*):
!                    'ZERO_AT_END'    -- Make AXIS_MIN or AXIS_MAX zero.
!                    'ZERO_SYMMETRIC' -- Make AXIS_MIN = -AXIS_MAX
!                    'GENERAL'        -- No restriction on min or max.
!                    'EXACT'          -- Use the data_min and data_max as the plot bounds.
!                                          Is ignored if data_min is very close to data_max.
!   axis_type   -- character(*), optional: Type of axis. 'LINEAR' or 'LOG'.
!                      Default is 'LINEAR'
!   slop_factor -- real(rp), optional: See qp_calc_axis_scale for info.
!
! Example:
!   call qp_calc_and_set_axis ('X', 352.0_rp, 378.0_rp, 4, 6, 'ZERO_AT_END')
!
! Gives for the x-axis:
!   %min         = 0.0_rp
!   %max         = 400.0_rp
!   %tick_min    = 0.0_rp
!   %tick_max    = 400.0_rp
!   %dtick       = 100.0_rp
!   %major_div   = 4
!   %places      = 0         ! places after the decimal point needed
!-

subroutine qp_calc_and_set_axis (axis_str, data_min, data_max, div_min, div_max, &
                                                        bounds, axis_type, slop_factor)

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
ax%type = string_option('LINEAR', axis_type)
call qp_calc_axis_params (data_min, data_max, div_min, div_max, ax)

end subroutine qp_calc_and_set_axis

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_calc_axis_params (data_min, data_max, div_min, div_max, axis)
!
! Subroutine to calculate a "nice" plot scale given the minimum and maximum
! of the data. This is similar to CALC_AXIS_SCALE but here the subroutine will
! pick the number of divisions.
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
!                       'EXACT'          -- Use the data_min and data_max as the plot bounds.
!                                             Is ignored if data_min is very close to data_max.
!
! Output:
!   axis       -- qp_axis_struct: Structure holding the axis parameters.
!     %places
!     %min, %max
!     %major_div
!     %tick_min, tick_max      
!     %dtick
!-

subroutine qp_calc_axis_params (data_min, data_max, div_min, div_max, axis)

implicit none

type (qp_axis_struct) axis

integer i, div_min, div_max, div_best
real(rp) data_max, data_min, d_max, d_min, score, score_max

! For log axis the ticks are always spaced a factor of 10 apart

d_min = min(data_min, data_max)
d_max = max(data_min, data_max)

if (axis%type == 'LOG') then
  axis%major_div = (div_min + div_max) / 2
  call qp_calc_axis_scale (d_min, d_max, axis, score)

else
  score_max = -1d20
  axis%major_div = 0

  do i = div_min, div_max
    axis%major_div = i
    call qp_calc_axis_scale (d_min, d_max, axis, score)
    if (score_max < score) then
      score_max = score
      div_best = i
    endif
  enddo

  axis%major_div = div_best
  call qp_calc_axis_scale (d_min, d_max, axis, score)
endif

end subroutine qp_calc_axis_params

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_calc_axis_scale (data_min, data_max, axis, niceness_score)
!
! Subroutine to calculate a "nice" plot scale given the minimum and maximum
! of the data. If data_max - data_min < 1d-30 then axis%max - axis%min will
! be at least axis%major_div * 1d-30.
!
! Note: If you are not sure how many divisions the axis should have,
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
!                       'ZERO_AT_END'    -- Make AXIS%MIN or AXIS%MAX zero.
!                       'ZERO_SYMMETRIC' -- Make AXIS%MIN = -AXIS%MAX
!                       'GENERAL'        -- No restriction on min or max.
!                       'EXACT'          -- Use the data_min and data_max as the plot bounds.
!                                             Is ignored if data_min is very close to data_max.
!     %major_div  -- Integer: How many divisions the axis is divided up into
!
! Output:
!   axis           -- qp_axis_struct: Structure holding the axis parameters.
!     %places        -- Integer: Number of axis%places after the decimal point needed
!                         to display the axis numbers. 
!     %min           -- Real(rp): Axis minimum.
!     %max           -- Real(rp): Axis maximum.
!   niceness_score -- Real(rp), optional: Score as to how "nice" 
!                      axis%min and axis%max are. The larger the number the nicer.
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

subroutine qp_calc_axis_scale (data_min, data_max, axis, niceness_score)

implicit none
           
type (qp_axis_struct) axis

integer div_eff, m_min, m_max, m
integer i, i1_min, i2_min, i1_max, i2_max, imin, imax, j, ave, div_max, div

real(rp), optional :: niceness_score
real(rp) data_max, data_min, r, p, a_min, a_max
real(dp) data_width, data_width10, min_width, max_score, score, aa

character(*), parameter :: r_name = 'qp_calc_axis_scale'

! Error check

if (axis%major_div < 1) then
  call out_io (s_abort$, r_name, '"AXIS%MAJOR_DIV" NUMBER IS LESS THAN 1! \i\ ', axis%major_div)
  axis%major_div = 1
endif

a_min = min(data_min, data_max)
a_max = max(data_min, data_max)

! 'LOG' axis

if (axis%type == 'LOG') then

  div_max = 1.4 * axis%major_div
  if (a_min <= 0) then
    call out_io (s_abort$, r_name, 'DATA IS NEGATIVE FOR LOG AXIS CALC: \es12.2\ ', a_min)
    if (global_com%exit_on_error) call err_exit
    axis%tick_min = 1
    axis%tick_max = 10
    axis%dtick = 10
    axis%min = 1
    axis%max = 10
    return
  endif

  a_min = log10(data_min)
  a_max = log10(data_max)
  m_min = floor(a_min+0.3)
  m_max = ceiling(a_max-0.3)

  div = 1 + (m_max - m_min - 1) / div_max

  if (div == 1) then
    if (m_min == m_max) then
      if (a_min+a_max < m_min+m_max) then
        m_min = m_min - 1
      else
        m_max = m_max + 1
      endif
    endif
  else
    m_min = m_min - modulo(m_min, div)
    m_max = m_max + modulo(-m_max, div)
  endif

  axis%tick_min = 10.0_rp**m_min
  axis%tick_max = 10.0_rp**m_max
  axis%dtick = 10.0_rp**div

  axis%min = min(axis%tick_min, data_min, data_max)
  axis%max = max(axis%tick_max, data_min, data_max)

  return

endif  ! log scale

!-----------------------------------------------------------------------
! If a non-log axis...

if (axis%bounds == 'EXACT' .and. &
          abs(data_max-data_min)>axis%major_div * max(abs(data_max)*1e-5, abs(data_min)*1e-5, 1e-29_rp)) then
  axis%tick_min = a_min
  axis%tick_max = a_max
  axis%min = a_min
  axis%max = a_max
  axis%dtick = (axis%tick_max - axis%tick_min) / axis%major_div
  call qp_calc_axis_places (axis)
  if (present(niceness_score)) then
    r = axis%dtick / 10d0**floor(1.0000001_rp*log10(axis%dtick))
    do i = 0, 6
      p = r * 10**i
      if (abs(nint(p) - p) < 1d-5) exit
    enddo
    niceness_score = -i
  endif
  return
endif

! Find width of data

select case (axis%bounds)
case ('ZERO_AT_END')
  a_max = max(abs(data_max), abs(data_min))
  a_min = 0
case ('ZERO_SYMMETRIC')
  a_max = max(abs(data_max), abs(data_min))
  a_min = -a_max
case ('GENERAL', 'EXACT')
  ! Nothing to do
case default
  call out_io (s_error$, r_name, 'I DO NOT UNDERSTAND "AXIS%BOUNDS": ' // axis%bounds)
  if (global_com%exit_on_error) call err_exit
end select

! 

data_width = a_max - a_min
a_min = a_min + 0.3_rp * data_width / max(4, axis%major_div)
a_max = a_max - 0.3_rp * data_width / max(4, axis%major_div)

! Find possible candidates
            
min_width = axis%major_div * max(abs(a_max)*1e-5, abs(a_min)*1e-5, 1e-29_rp)
data_width = max(data_width, min_width)
data_width10 = 10d0**(floor(log10(data_width))-1)
div_eff = axis%major_div

if (axis%bounds == 'ZERO_AT_END') then
  i1_min = 0
  i2_min = 0
else
  i2_min = floor(a_min / data_width10)
  i1_min = floor(a_min / (100 * data_width10)) * 100
  i1_min = min(i1_min, i2_min - 10*div_eff)
endif

aa = a_max/data_width10 + 1
i1_max = floor(aa)
aa = a_max / (100 * data_width10) + 1
i2_max = max(floor(aa) * 100, i1_max + 10*div_eff)

ave = nint((a_max + a_min)/data_width10)
i1_min = max(i1_min, 2*i2_min - i1_max)
i1_min = min(i1_min, (ave - axis%major_div)/2 - 1)

i2_max = min(i2_max, 2*i1_max - i2_min)
i2_max = max(i2_max, i1_min+axis%major_div+1)

! Go through and rate all the possibilities and choose the one with the highest score.

max_score = -1000
do imin = i1_min, i2_min
  do imax = i1_max, i2_max
    score = qp_axis_niceness (imin, imax, div_eff)
    if (score > max_score) then
      max_score = score
      axis%tick_min = imin * data_width10
      axis%tick_max = imax * data_width10
    endif
  enddo
enddo

if (present(niceness_score)) niceness_score = max_score 

! adjust the scale if necessary

if (axis%bounds == 'ZERO_AT_END' .and. data_min + data_max < 0) then
  axis%tick_min = -axis%tick_max
  axis%tick_max = 0
elseif (axis%bounds == 'ZERO_SYMMETRIC') then
  axis%tick_min = -axis%tick_max
endif

! find number of places needed

axis%min = min(axis%tick_min, data_min, data_max)
axis%max = max(axis%tick_max, data_min, data_max)

if (axis%tick_min == axis%tick_max) then
  axis%dtick = 1
else
  axis%dtick = (axis%tick_max - axis%tick_min) / axis%major_div
endif

call qp_calc_axis_places (axis)

end subroutine qp_calc_axis_scale

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
integer del, i, j, im, d0
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

d0 = (imax - imin) / divisions
do j = imin, imax, d0
  if (mod(j, 5) == 0) score = score + 4.0 / (divisions + 1)
  if (mod(j, 20) == 0) score = score + 1.0 / (divisions + 1)
  if (mod(j, 2) == 0) score = score + 4.0 / (divisions + 1)
  if (j > 6 .and. mod(j, 2) == 1) score = score - 8.0 / (divisions + 1)
  if (j == 0) score = score + 4
enddo

score = score + 100 - 50 * log10(float(imax - imin))

select case (d0)
case (2, 3, 4, 5, 8, 10, 20, 25, 30, 40, 50, 80, 100)
  score = score + 3
case (6, 15, 60)
  score = score + 2
end select

end function qp_axis_niceness 

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_calc_axis_places (axis)
!
! Subroutine to calculate the number of decmal places needed to display the
! axis numbers. Note: If axis%min > axis%max on input they will be reversed.
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
real(rp) num, num2, max_d, effective_zero

! LOG scale does not use places.
! If min = max or major_div < 1 then cannot do a calculation so don't do anything.

if (axis%type == 'LOG') return
if (axis%min == axis%max) return
if (axis%major_div < 1) return
if (axis%dtick == 0) return

if (axis%min > axis%max) then
  num = axis%min
  axis%min = axis%max
  axis%max = num
endif

! Sort true min and max.
! Huge(0) is used here so that nint(num2) will not overflow.

max_d = min(10.0_rp**qp_com%max_digits, 0.999_rp * huge(0))
effective_zero = max(abs(axis%tick_max), abs(axis%tick_min)) / max_d

! First calculation: Take each axis number and find how many digits it has.
! The number of places is the maximum number of digits needed to represent
! all the numbers on the axis with a limit of qp_com%max_digits digits maximum for any one number.

axis%places = -1000
do i = 0, axis%major_div
  num = axis%tick_min + i * axis%dtick
  if (abs(num) < effective_zero) cycle  ! Ignore zero
  axis%places = max(axis%places, floor(-log10(abs(num))))
  do
    num2 = abs(num * 10d0**axis%places)
    if (num2 > max_d) exit
    if (abs(num2 - nint(num2)) < 0.01) exit
    axis%places = axis%places + 1
  enddo
enddo

! Second calculation: Places based upon the distance between ticks.
! The number of places returned by the subroutine is the maximum of the two calculations

axis%places = max(axis%places, floor(-log10(axis%dtick)+0.9))

end subroutine qp_calc_axis_places

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_graph_limits ()
!
! Subroutine to calculate the offsets for the graph.
! Note: This subroutine is for internal use only.
!-


subroutine qp_set_graph_limits ()

implicit none

type (qp_rect_struct), pointer :: graph

real(rp) y1, y2

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
! Also: The graph position sets the clip boundries. We enlarge the graph position
! slightly to prevent clipping of lines that are slightly outside the graph.

if (graph%x1 < graph%x2 .and. graph%y1 < graph%y2) then
  y1 = graph%y1 - 2 * qp_com%dflt_axis_slop_factor * (graph%y2 - graph%y1)
  y2 = graph%y2 + 2 * qp_com%dflt_axis_slop_factor * (graph%y2 - graph%y1)
  call qp_set_graph_position_basic (graph%x1, graph%x2, y1, y2)
endif

end subroutine qp_set_graph_limits

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_clear_page ()
!
! Subroutine to clear all drawing from the displayed page.
! If outputting to a file: End current page and start a new one.
!-

subroutine qp_clear_page ()

type (qp_rect_struct) border  

! Clear all graphics.

call qp_clear_page_basic

! For gif must paint a white background.

if (qp_com%page_type(1:3) == 'GIF') then
  call qp_paint_rectangle (qp_com%page%x1, qp_com%page%x2, &
                           qp_com%page%y1, qp_com%page%y2, &
                           'INCH', color = 'white', fill_pattern = 'solid_fill')
endif

end subroutine qp_clear_page

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
                         'INCH/PAGE', color = 'white', fill_pattern = 'solid_fill')

end subroutine qp_clear_box

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_paint_rectangle (x1, x2, y1, y2, units, color, fill_pattern)
!
! Subroutine to paint a rectangular region a specified color.
! The default color is the background color (white$).
! 
! Also see qp_draw_rectangle
!
! Input:
!   x1  -- Real(rp): Left edge
!   x2  -- Real(rp): Right edge
!   y1  -- Real(rp): Bottom edge
!   y2  -- Real(rp): Top edge.
!   units        -- Character(*), optional: Units of returned numbers.
!                     Default = 'DATA/GRAPH/LB'
!   color        -- character(*), optional: Color to paint the rectangle.
!                     Default is to use the symbol color.
!   fill_pattern -- character(*), optional: Fill pattern. 
!                   Default is to use the symbol color.
!-

subroutine qp_paint_rectangle (x1, x2, y1, y2, units, color, fill_pattern)
              
implicit none

real(rp) x1, x2, y1, y2
real(rp) x1_inch, x2_inch, y1_inch, y2_inch
character(*), optional :: color, units, fill_pattern

!

if (x1 == x2 .or. y1 == y2) return

call qp_to_inch_abs (x1, y1, x1_inch, y1_inch, units)
call qp_to_inch_abs (x2, y2, x2_inch, y2_inch, units)
 
call qp_paint_rectangle_basic (x1_inch, x2_inch, y1_inch, y2_inch, &
              qp_string_to_enum(string_option(qp_com%symbol%color, color), 'color'), &
              qp_string_to_enum(string_option(qp_com%symbol%fill_pattern, fill_pattern), 'fill_pattern'))

end subroutine qp_paint_rectangle

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

end subroutine qp_set_box

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
  
end subroutine qp_set_page_border_to_box

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

end subroutine qp_set_page_border

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
!   units  -- Character(*), optional: Units of x and y
!                   Default is: 'DATA/GRAPH/LB'
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
    if (qp_com%plot%xx%min == 0) then
      x = 1
    else
      x = (qp_com%plot%xx%max / qp_com%plot%xx%min) ** dx
    endif
  else
    x = dx * (qp_com%plot%xx%max - qp_com%plot%xx%min) 
  endif

  dy = y_inch / (qp_com%graph%y2 - qp_com%graph%y1)
  if (qp_com%plot%yy%type == 'LOG') then
    if (qp_com%plot%yy%min == 0) then
      y = 1
    else
      y = (qp_com%plot%yy%max / qp_com%plot%yy%min) ** dy
    endif
  else
    y = dy * (qp_com%plot%yy%max - qp_com%plot%yy%min) 
  endif
endif

end subroutine qp_from_inch_rel

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
!                   Default is: 'DATA/GRAPH/LB'
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

end subroutine qp_from_inch_abs

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
!                   Default is: 'DATA/GRAPH/LB'
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
character(*), parameter :: r_name = 'qp_to_inch_rel'

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
    if (x <= 0 .or. qp_com%plot%xx%max == qp_com%plot%xx%min) then ! x could be a dummy number so just return a dummy number
      x_inch = 1  ! Dummy number
    else
      x_inch = log(x) * dgx / (log(qp_com%plot%xx%max) - log(qp_com%plot%xx%min))
    endif
  elseif (qp_com%plot%xx%max == qp_com%plot%xx%min) then
    x_inch = 0
  else
    x_inch = x * dgx / (qp_com%plot%xx%max - qp_com%plot%xx%min)
  endif
  dgy = qp_com%graph%y2 - qp_com%graph%y1
  if (qp_com%plot%yy%type == 'LOG') then
    if (y <= 0 .or. qp_com%plot%yy%max == qp_com%plot%yy%min) then ! y could be a dummy number so just return a dummy number
      y_inch = 1 ! Dummy number
    else
      y_inch = log(y) * dgy / (log(qp_com%plot%yy%max) - log(qp_com%plot%yy%min))
    endif
  elseif (qp_com%plot%yy%max == qp_com%plot%yy%min) then
    y_inch = 0
  else
    y_inch = y * dgy / (qp_com%plot%yy%max - qp_com%plot%yy%min)
  endif
endif

end subroutine qp_to_inch_rel

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
!                   Default is: 'DATA/GRAPH/LB'
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

end subroutine qp_to_inches_rel

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
    if (qp_com%plot%xx%min == 0) then  ! Can happen when modifying plot parameters
      x0 = 0
    else
      x0 = x / qp_com%plot%xx%min
    endif
  else
    x0 = x - qp_com%plot%xx%min
  endif

  if (qp_com%plot%yy%type == 'LOG') then
    if (qp_com%plot%yy%min == 0) then  ! Can happen when modifying plot parameters
      y0 = 0
    else
      y0 = y / qp_com%plot%yy%min
    endif
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

end subroutine qp_to_inch_abs

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

end subroutine qp_to_inches_abs

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

end subroutine qp_convert_point_rel

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

end subroutine qp_convert_point_abs

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

end subroutine qp_convert_rectangle_rel

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

end subroutine qp_join_units_string

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
character(*), parameter :: r_name = 'qp_split_units_string'

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

if (all(u_type /= ['DATA  ', 'MM    ', 'INCH  ', 'POINTS', '%     '])) then
  call out_io (s_error$, r_name, 'BAD UNITS TYPE: "' // trim(units) // '"')
  if (global_com%exit_on_error) call err_exit
endif

! get region

call string_trim (u(ix+1:), u, ix)
if (ix == 0) return
region = u(:ix)

if (all(region /= ['PAGE ', 'BOX  ', 'GRAPH'])) then
  call out_io (s_error$, r_name, 'BAD REGION: "' // trim(units) // '"')
  if (global_com%exit_on_error) call err_exit
endif

! get corner

call string_trim (u(ix+1:), u, ix)
if (ix == 0) return
corner = u(:ix)

if (all(corner /= ['LB', 'LT', 'RB', 'RT'])) then
  call out_io (s_error$, r_name, 'BAD CORNER: "' // trim(units) // '"')
  if (global_com%exit_on_error) call err_exit
endif

call string_trim (u(ix+1:), u, ix)
if (ix /= 0) then
  call out_io (s_error$, r_name, 'EXTRA CHARACTERS IN UNITS: "' // trim(units) // '"')
  if (global_com%exit_on_error) call err_exit
endif

end subroutine qp_split_units_string

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

end subroutine qp_set_graph_placement

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

end subroutine qp_set_margin

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_rectangle (x1, x2, y1, y2, units, color, width, line_pattern, clip, style)
!
! Subroutine to draw a rectangular box.
!
! Also see qp_paint_rectangle
!
! Input:
!   x1, y1       -- Real(rp): (x, y) corner of box.
!   x2, y2       -- Real(rp): (x, y) opposite corner of box.
!   units        -- Character(*), optional: Units of x and y.
!                     Default is: 'DATA/GRAPH/LB'
!                     See quick_plot writeup for more details.
!   color        -- Character(*), optional: Color index for the box
!   width        -- Integer, optional: Width of the line. Default = 1
!   line_pattern -- Character(*), optional: Line type (dashed$, etc). 
!   clip         -- Logical, optional: Clip at the graph boundary?
!   style        -- Character(*): Default line style to use if not specified by the other arguments.
!                     Default is 'STD'. See qp_set_line_attrib for more details.
!-

subroutine qp_draw_rectangle (x1, x2, y1, y2, units, color, width, line_pattern, clip, style)

implicit none

real(rp) x1, y1, x2, y2

integer, optional :: width

character(*), optional :: units, style, color, line_pattern


logical, optional :: clip

!

call qp_draw_polyline ([x1, x1, x2, x2, x1], [y1, y2, y2, y1, y1], &
                                            units, width, color, line_pattern, clip, style)

end subroutine qp_draw_rectangle

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_arrow (r1, r2, units, color, head_size, head_type, head_angle, head_barb)
!
! Draws an arrow 
!
! Input:
!   r1(2), r2(2)  -- real(rp): Arrow tail and head (x, y) coordinates.
!   units         -- character(*), optional: Units of r1 and r2. Default is: 'DATA/GRAPH/LB'
!   color         -- character(*), optional: Arrow color.
!   head_size     -- real(rp), optional: Size of the arrow.
!   head_type     -- character(*), optional: Arrow head type: filled_arrow_head$ or outline_arrow_head$ 
!   head_angle    -- real(rp), optional: Acute angle of the arrow point in degrees.
!   head_barb     -- real(rp), optional: Fraction of triangular arrow head that is cut away from the back.
!-

subroutine qp_draw_arrow (r1, r2, units, color, head_size, head_type, head_angle, head_barb)

implicit none
              
real(rp) r1(2), r2(2), rinch1(2), rinch2(2)
real(rp), optional :: head_size, head_angle, head_barb

character(*), optional :: color, head_type

character(*), optional :: units

!

call qp_save_state (.true.)

call qp_set_arrow_attrib (color, head_size, head_type, head_angle, head_barb)
call qp_to_inch_abs (r1(1), r1(2), rinch1(1), rinch1(2), units)
call qp_to_inch_abs (r2(1), r2(2), rinch2(1), rinch2(2), units)
call qp_draw_arrow_basic (r1, r2, qp_com%arrow)

call qp_restore_state

end subroutine qp_draw_arrow

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_symbol (x, y, units, type, height, color, fill_pattern, line_width, clip)
!
! Draws a symbol at (x, y) 
!
! Also see: qp_draw_symbols.
!
! Input:
!   x, y         -- real(rp): Symbol coordinates.
!   units        -- character(*), optional: Units of (x, y). Default is: 'DATA/GRAPH/LB'
!   type         -- character(*), optional: Symbol type. 
!   height       -- real(rp), optional: Size of the symbol.
!   color        -- character(*), optional: Symbol color.
!   fill_pattern -- character(*), optional: fill pattern.
!   line_width   -- integer, optional: Line width.
!   clip         -- logical, optional: Clip at the graph boundary?
!-

subroutine qp_draw_symbol (x, y, units, type, height, color, fill_pattern, line_width, clip)

implicit none
              
integer, optional :: line_width

real(rp) x, y, x_inch, y_inch
real(rp), optional :: height

character(*), optional :: type, color, fill_pattern
character(*), optional :: units

logical, optional :: clip

!

call qp_save_state (.true.)

call qp_set_symbol_attrib (type, height, color, fill_pattern, line_width, clip)
call qp_to_inch_abs (x, y, x_inch, y_inch, units)
call qp_draw_symbol_basic (x_inch, y_inch, qp_string_to_enum(qp_com%symbol%type, 'symbol_type'))

call qp_restore_state

end subroutine qp_draw_symbol

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_symbols (x, y, units, type, height, color, fill_pattern, line_width, clip, symbol_every)
!
! Draws symbols at a set of (x, y) points. 
! Data units are assumed.
!
! Also see: qp_draw_symbol.
!
! Input:
!   x(:), y(:)   -- Real(rp): Symbol coordinates.
!   units        -- character(*), optional: Units of (x, y). Default is: 'DATA/GRAPH/LB'
!   type         -- Character(*), optional: Symbol type. 
!   height       -- Real(rp), optional: Size of the symbol.
!   color        -- Character(*), optional: Symbol color.
!   fill_pattern -- Character(*), optional: fill pattern.
!   line_width   -- Integer, optional: Line width.
!   clip         -- Logical, optional: Clip at the graph boundary?
!   symbol_every -- Integer, optional: 
!                               0  --> Do not draw symbols.
!                               1  --> Draw symbols (Default).
!                               2  --> Draw every 2nd symbol starting by drawing the first symbol.
!                               etc.
!-


subroutine qp_draw_symbols (x, y, units, type, height, color, fill_pattern, line_width, clip, symbol_every)

implicit none

integer, optional :: symbol_every, line_width
integer i, i_skip

real(rp) x(:), y(:)
real(rp), optional :: height

character(*), optional :: units, type, color, fill_pattern
character(*), parameter :: r_name = 'qp_draw_symbols'

logical, optional :: clip

!

if (size(x) /= size(y)) then
  call out_io (s_error$, r_name, 'X, Y COORD VECTORS HAVE UNEQUAL LENGTH!')
  if (global_com%exit_on_error) call err_exit
endif

i_skip = integer_option(1, symbol_every)
if (i_skip < 1) return

do i = 1, size(x)
  if (mod(i-1, i_skip) /= 0) cycle
  call qp_draw_symbol (x(i), y(i), units, type, height, color, fill_pattern, line_width, clip)
enddo


end subroutine qp_draw_symbols

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

end subroutine qp_set_graph

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_graph (x, y, x_lab, y_lab, title, draw_line, symbol_every, clip)
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

subroutine qp_draw_graph (x_dat, y_dat, x_lab, y_lab, title, draw_line, symbol_every, clip)

implicit none      

real(rp) x_dat(:), y_dat(:)                        

integer, optional :: symbol_every
logical, optional :: draw_line, clip

character(*), optional :: x_lab, y_lab, title  
character(*), parameter :: r_name = 'qp_draw_graph'

! Error check

if (qp_com%plot%xx%max == qp_com%plot%xx%min) then
  call out_io (s_error$, r_name, 'X_MAX = X_MIN')
  if (global_com%exit_on_error) call err_exit
endif

if (qp_com%plot%yy%max == qp_com%plot%yy%min) then
  call out_io (s_error$, r_name, 'Y_MAX = Y_MIN')
  if (global_com%exit_on_error) call err_exit
endif
           
! Draw the axes and the data

call qp_draw_axes (x_lab, y_lab, title)
call qp_draw_data (x_dat, y_dat, draw_line, symbol_every, clip)

end subroutine qp_draw_graph

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_data (x, y, draw_line, symbol_every, clip)
!
! Subroutine to draw a symbol at the data points and optionaly connect
! the symbols with line segments.
! The line segments use the PLOT line style attributes (see qp_set_line_attrib).
!
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

integer i_skip, n 
integer, optional :: symbol_every
logical, optional :: draw_line, clip

! init

call qp_save_state (.true.)

! plot a polyline

if (logic_option (.true., draw_line)) then
  call qp_set_line_attrib ('PLOT', clip = clip)
  i_skip = integer_option(1, symbol_every)
  if (i_skip > 1) then
    n = size(x_dat)
    call qp_draw_polyline_no_set (x_dat(1:n:i_skip), y_dat(1:n:i_skip))
  else
    call qp_draw_polyline_no_set (x_dat, y_dat)
  endif
endif

! plot symbols

call qp_draw_symbols (x_dat, y_dat, symbol_every = symbol_every, clip = clip)

! 

call qp_restore_state

end subroutine qp_draw_data

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

end subroutine qp_draw_axes

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

end subroutine qp_draw_graph_title

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_histogram (x_dat, y_dat, fill_color, fill_pattern, line_color, clip)
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
!   fill_color   -- Character(*), optional: Color of fill pattern and outline.
!                     If fill_color is transparent$ no color is added.
!                     Default = black$
!   fill_pattern -- Character(*), optional: Default is set by symbol fill pattern.
!   line_color   -- Character(*), optional: Outline of the rectangles color
!                     Default is black$. 
!   clip         -- Logical, optional: Clip at the graph boundary?
!                     Default is set by qp_set_line_attrib.
!-

subroutine qp_draw_histogram (x_dat, y_dat, fill_color, fill_pattern, line_color, clip) 

implicit none

integer i, n, n_min, n_max

real(rp) x_dat(:), y_dat(:)
real(rp) :: xh(2*size(x_dat)+2), yh(2*size(x_dat)+2)

logical, optional :: clip

character(*), optional :: line_color, fill_color, fill_pattern
character(*), parameter :: r_name = 'qp_draw_histogram'

! error check

if (qp_com%plot%xx%max == qp_com%plot%xx%min) then
  call out_io (s_error$, r_name, 'X_MAX = X_MIN')
  if (global_com%exit_on_error) call err_exit
endif

if (qp_com%plot%yy%max == qp_com%plot%yy%min) then
  call out_io (s_error$, r_name, 'Y_MAX = Y_MIN')
  if (global_com%exit_on_error) call err_exit
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
yh(3:2*n+1:2) = y_dat(n_min:n_max)

yh(2*n+2) = 0

if (logic_option (qp_com%clip, clip)) then
  where (yh < qp_com%plot%yy%min) yh = qp_com%plot%yy%min
  where (yh > qp_com%plot%yy%max) yh = qp_com%plot%yy%max
endif

! Draw

call qp_save_state (.true.)

call qp_set_line_attrib ('PLOT', clip = .false.)

if (string_option('black', fill_color) /= 'transparent') then
  do i = 1, n
    call qp_paint_rectangle (xh(2*i), xh(2*i+1), yh(1), yh(2*i), 'DATA', &
                                      string_option('black', fill_color), fill_pattern)
  enddo
endif

call qp_draw_polyline (xh(1:2*n+2), yh(1:2*n+2), color = string_option('black', line_color))

call qp_restore_state

end subroutine qp_draw_histogram

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
  ! This test should not be needed but there is some strange memory bug in cesrv linking shared
  ! and this test gets around it.
  if (text(i) == '') cycle
  call qp_draw_text_no_set (text(i), xc, yc, 'INCH/PAGE/LB')
  yc = yc - 1.5 * height
enddo

call qp_restore_state

end subroutine qp_draw_text_legend

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_curve_legend (x_origin, y_origin, units, line, line_length, 
!                  symbol, text, text_offset, row_gap_scale, draw_line, draw_symbol, draw_text)
!
! Subroutine to draw a legend with each line in the legend having
!   a line, a symbol, some text.
!
! Input:
!   x_origin      -- Real(rp), optional: x-postion of start of the first line.
!   y_origin      -- Real(rp), optional: y-postion of start of the first line.
!   units         -- Character(*), optional: Units of x_origin, y_origin.
!                      Default is: 'DATA/GRAPH/LB'
!                      See quick_plot writeup for more details.
!   line(:)       -- qp_line_struct, optional: Array of lines.
!                      Set line(i)%width < 0 to suppress drawing of the i^th line
!   line_length   -- Real(rp), optional: Length of the line in points. Default is 72 pts (~ 1 inch).
!   symbol(:)     -- qp_symbol_struct, optional: Array of symbols.
!                      Set symbol(i)%type < 0 to suppress drawing of the i^th symbol.
!   text(:)       -- Character(*), optional: Array of text lines.
!   text_offset   -- Real(rp), optional: Horizontal offset in points between the line and the text.
!                      Default is 10 pt.
!   row_gap_scale -- Real(rp), optional: Scale factor for gap size between entries in a legend. Default is 1.
!   draw_line     -- Logical, optional: Draw lines? Default is True if line arg is present.
!                      Line style set by the LEGEND line style. See qp_set_line_attrib.
!   draw_symbol   -- Logical, optional: Draw symbols? Default is True if symbol arg is present.
!   draw_text     -- Logical, optional: Draw text? Default is True if text arg is present.
!-

subroutine qp_draw_curve_legend (x_origin, y_origin, units, line, line_length, &
                      symbol, text, text_offset, row_gap_scale, draw_line, draw_symbol, draw_text)

implicit none

type (qp_line_struct), optional :: line(:)
type (qp_symbol_struct), optional :: symbol(:)

real(rp) x_origin, y_origin
real(rp), optional :: line_length, text_offset, row_gap_scale
real(rp) height, xc, yc, yc2, line_len, dummy, text_off

integer i, n_rows

character(*), optional :: units, text(:)
character(20) :: r_name = 'qp_draw_curve_legend'

logical, optional :: draw_line, draw_symbol, draw_text
logical has_line, has_symbol, has_text

!

call qp_save_state (.true.)

call qp_set_text_attrib ('LEGEND')
call qp_set_line_attrib ('LEGEND')
height = real_option(1.0_rp, row_gap_scale)*qp_text_height_to_inches(qp_com%this_text%height) 

call qp_to_inch_abs (x_origin, y_origin, xc, yc, units)

! Find out how many rows to draw

n_rows = 0
has_text = .false.
if (present(text) .and. logic_option(.true., draw_text)) then
  has_text = .true.
  n_rows = size(text)
  call qp_to_inch_rel (real_option(10.0_rp, text_offset), 0.0_rp, text_off, dummy, 'POINTS')
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
  call qp_to_inch_rel (real_option(72.0_rp, line_length), 0.0_rp, line_len, dummy, 'POINTS')
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

  if (has_line) then
    if (line(i)%width > -1) then
      call qp_set_line ('STD', line(i))
      call qp_draw_line (xc, xc + line_len, yc2, yc2, 'INCH/PAGE/LB')
    endif
  endif

  if (has_symbol) then 
    if (symbol(i)%type /= '') then
      call qp_set_symbol (symbol(i))
      call qp_draw_symbol (xc + line_len/2, yc2, 'INCH/PAGE/LB')
    endif
  endif

  if (has_text) then
    call qp_set_text_attrib ('LEGEND')
    call qp_draw_text_no_set (text(i), xc + line_len + text_off, yc, 'INCH/PAGE/LB')
  endif

  yc = yc - 1.5 * height
enddo

call qp_restore_state

end subroutine qp_draw_curve_legend

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_ellipse (x0, y0, r_x, r_y, theta_xy, angle0, 
!                              del_angle, units, width, color, line_pattern, clip)
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
!   x0, y0       -- Real(rp): Center of arc.
!   r_x          -- Real(rp): Horizontal radius.
!   r_y          -- Real(rp): Vertical radius.
!   theta_xy     -- Real(rp), optional: Ellipse tilt. Default = 0.
!   angle0       -- Real(rp), optional: Starting angle of arc to draw in radians.
!                    Default is 0
!   del_angle    -- Real(rp), optional: Angle of arc to draw. 
!                    Default is 2pi.
!   units        -- Character(*), optional: Units of x, y.
!                    Default is: 'DATA/GRAPH/LB'
!                    See quick_plot writeup for more details.
!   width        -- Integer, optional: Width of line
!   color        -- Character(*), optional: Line color.
!   line_pattern -- Character(*), optional: Line type. Currently ignored.
!   clip         -- Logical, optional: Clip at graph boundary?
!-

subroutine qp_draw_ellipse (x0, y0, r_x, r_y, theta_xy, angle0, &
                          del_angle, units, width, color, line_pattern, clip)

implicit none

real(rp) x0, y0, r_x, r_y
real(rp), optional :: theta_xy, angle0, del_angle
real(rp) x(1000), y(1000), ang22, del, ang, cos_xy, sin_xy
real(rp) xx0, yy0, rr_x, rr_y, ang0, del_ang, t_xy, dx, dy

integer, optional :: width
integer i

character(*), optional :: units, color, line_pattern

logical, optional :: clip

!

if (present(line_pattern)) i = len(line_pattern)   ! so compiler will not complain

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

call qp_set_line_attrib ('STD', width, color, 'solid', clip)  ! solid line
call qp_draw_polyline_no_set (x(1:i), y(1:i), 'INCH/PAGE/LB')

call qp_restore_state

end subroutine qp_draw_ellipse

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_circle (x0, y0, r, angle0, del_angle, units, width, color, line_pattern, clip)
!
! Subroutine to plot a section of a circle.
! Drawn is:
!     (x, y) = (x0, y0) + r * (cos(theta), sin(theta))
! where theta goes from angle0 to angle0 + del_angle.
!
! Note: Currently this routine can only draw solid lines.
!
! Input:
!   x0, y0       -- Real(rp): Center of arc.
!   r            -- Real(rp): Radius.
!   angle0       -- Real(rp), optional: Starting angle of arc to draw in radians.
!                    Default is 0
!   del_angle    -- Real(rp), optional: Angle of arc to draw. 
!                    Default is 2pi.
!   units        -- Character(*), optional: Units of x, y.
!                    Default is: 'DATA/GRAPH/LB'
!                    See quick_plot writeup for more details.
!   width        -- Integer, optional: Width of line
!   color        -- Character(*), optional: Line color.
!   line_pattern -- Character(*), optional: Line type. 
!                    Currently can only be 1 (solid line).
!   clip         -- Logical, optional: Clip at graph boundary?
!-

subroutine qp_draw_circle (x0, y0, r, angle0, del_angle, units, width, color, line_pattern, clip)

implicit none

real(rp) x0, y0, r
real(rp), optional :: angle0, del_angle

integer, optional :: width

character(*), optional :: units, color, line_pattern

logical, optional :: clip

!

call qp_draw_ellipse (x0, y0, r, r, 0.0_rp, angle0, del_angle, units, width, color, line_pattern, clip)

end subroutine qp_draw_circle

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_polyline (x, y, units, width, color, line_pattern, clip, style)
!
! Subroutine to draw a polyline.
!
! Input:
!   x(:), y(:)    -- Real(rp): (x, y) points for the polyline.
!   units         -- Character(*), optional: Units of x, y.
!                      Default is: 'DATA/GRAPH/LB'
!                      See quick_plot writeup for more details.
!   width         -- Integer, optional: Width of line
!   color         -- Character(*), optional: Line color.
!   line_pattern  -- Character(*), optional: Line type. 
!   clip          -- Logical, optional: Clip at graph boundary?
!   style         -- Character(*), optional: Part of drawing this is for.
!                   'PLOT'     -- Plot data lines.
!                   'GRID'     -- Graph grid.
!                   'AXIS'     -- Graph axis.
!                   'LEGEND    -- Line legend.
!                   'STD'      -- Everything else. Default.
!-

subroutine qp_draw_polyline (x, y, units, width, color, line_pattern, clip, style)

implicit none

real(rp) :: x(:), y(:)
real(rp) :: xd(size(x)), yd(size(y))

integer, optional :: width

character(*), optional :: units, style, color, line_pattern

logical, optional :: clip

!

call qp_save_state (.true.)
call qp_set_line_attrib (style, width, color, line_pattern, clip)
call qp_to_inches_abs (x, y, xd, yd, units)
call qp_draw_polyline_basic (xd, yd)
call qp_restore_state

end subroutine qp_draw_polyline

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

end subroutine qp_draw_polyline_no_set

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_line (x1, x2, y1, y2, units, width, color, line_pattern, clip, style)
!
! Subroutine to draw a line.
!
! Input:
!   x1, x2        -- Real(rp): X-coords of line endpoints
!   y1, y2        -- Real(rp): Y-coords of line endpoints
!   units         -- Character(*), optional: Units of x, y.
!                      Default is: 'DATA/GRAPH/LB'
!                      See quick_plot writeup for more details.
!   width         -- Integer, optional: Width of line
!   color         -- Character(*), optional: Line color.
!   line_pattern  -- Character(*), optional: Line type. 
!   clip          -- Logical, optional: Clip at graph boundary?
!   style         -- Character(*): Default line style to use if not specified by the other arguments.
!                      Default is 'STD'. See qp_set_line_attrib for more details.
!-

subroutine qp_draw_line (x1, x2, y1, y2, units, width, color, line_pattern, clip, style)

implicit none

real(rp) :: x1, x2, y1, y2

integer, optional :: width

character(*), optional :: units, style, color, line_pattern

logical, optional :: clip

!

if (x1 == x2 .and. y1 == y2) return

call qp_save_state (.true.)
call qp_set_line_attrib (style, width, color, line_pattern, clip)
call qp_draw_polyline_no_set ([x1, x2], [y1, y2], units)
call qp_restore_state

end subroutine qp_draw_line

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
! If there has been no previous call, this is an error.
!
! The scale argument can be used to expand or shrink the drawing.
!
! Note: pgplot does not do a good job with gif files. Consider making
! a postscript file and using the unix command pstogif to convert it.
!
! Input:
!   page_type -- Character(*). Device name for the type of plot.
!                 TYPE is passed to GG_SETUP. E.g.
!                 TYPE = 'GIF'   --> To create a gif file.
!                 TYPE = 'PDF'   -- > pdf file. [plplot only.]
!                 TYPE = 'PS'    --> To create a Color PostScript file.
!                 TYPE = 'X'     --> Open an X-window.
!   x_len     -- Real(rp), optional: Page horizontal width before scaling. Default: See above.
!   y_len     -- Real(rp), optional: Page vertical width before scaling. Default: See above.
!   units     -- Character(*), optional: units for x_len and y_len.
!                    'MM'     - milimeters.
!                    'INCH'   - Inches (default unless default_set_units is changed).
!                    'POINTS' - Point
!   plot_file -- Character(*), optional: Name for the plot file.
!                    Default is: 'quick_plot.ps' or 'quick_plot.gif'
!   scale     -- real(rp), optional: If positive then the entire page will be scaled.
!                    Default is 1.0.
!
! Output:
!   i_chan    -- Inteter, optional: Plot channel. Like a unit number for a fortran OPEN.
!                 To be used with qp_select_page.
!-

subroutine qp_open_page (page_type, i_chan, x_len, y_len, units, plot_file, scale)

implicit none

real(rp), optional :: x_len, y_len, scale
real(rp) x_inch, y_inch, x_page, y_page, page_scale, x, y

integer, optional :: i_chan
integer ix

character(*) page_type
character(*), optional :: units, plot_file
character(*), parameter :: r_name = 'qp_open_page'

logical saved_state

! Init

call qp_save_state (.false.)

qp_com%dflt_units = dflt_set$
qp_com%page_type = page_type

! set the name for the output plot file.

if (qp_com%page_type(1:3) == 'GIF') then
  qp_com%plot_file = 'quick_plot.gif'
elseif (qp_com%page_type(1:3) == 'PDF') then
  qp_com%plot_file = 'quick_plot.pdf'
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

if (present(y_len)) then
  call qp_to_inch_rel (x_len, y_len, x_inch, y_inch, units)
else
  call qp_to_inch_rel (1.0_rp, 1.0_rp, x_inch, y_inch, '%PAGE') ! Get present size.
endif

if (x_inch == 0) then
  x_inch = print_page_short_len
  y_inch = print_page_long_len
endif

qp_com%page%x1 = 0
qp_com%page%y1 = 0
qp_com%page%x2 = x_inch
qp_com%page%y2 = y_inch

! Calculate the page_scale

page_scale = real_option(1.0_rp, scale)

call qp_open_page_basic (qp_com%page_type, x_inch * page_scale, y_inch * page_scale, &
                        qp_com%plot_file, x_page, y_page, i_chan, page_scale)

! set the graph parameters

call qp_set_box (1, 1, 1, 1)
qp_com%dflt_units = dflt_draw$

! set a white background for gif

if (qp_com%page_type(1:3) == 'GIF') call qp_clear_page

end subroutine qp_open_page

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

end subroutine qp_select_page

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_close_page ()
!
! Subroutine to finish plotting on a page.
! For X this closes the window.
! You will need to call qp_open_page to do further graphics.
!-

subroutine qp_close_page ()

implicit none

character(*), parameter :: r_name = 'qp_close_page'
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

end subroutine qp_close_page

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
!   color      -- Character(*), optional: Color index for the box
!   angle      -- Real(rp), optional: Angle to the horizontal (in degrees). 
!                   Positive angle is CCW.
!   background -- Character(*), optional: Background color.
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

character(*) text
character(*), optional :: units, justify, color, background

logical, optional :: uniform_spacing

!

call qp_save_state (.true.)

call qp_set_text_attrib ('TEXT', height, color, background, uniform_spacing)
if (present(spacing_factor)) qp_com%text_spacing_factor = spacing_factor

call qp_draw_text_no_set (text, x, y, units, justify, angle)

call qp_restore_state

end subroutine qp_draw_text

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
character(*), parameter :: r_name = 'qp_draw_text_no_set'

!

if (text == '') return

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
    call out_io (s_error$, r_name, 'UNKNOWN "JUSTIFY": ' // justify)
    if (global_com%exit_on_error) call err_exit
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
  call qp_draw_text_basic (text, len_trim(text), x_inch, y_inch, ang, qp_justify(justify))
endif

end subroutine qp_draw_text_no_set

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
character(*), parameter :: r_name = 'qp_justify'

!

horiz_justy = 0.0

if (present(justify)) then
  if (justify(1:1) == 'C') then
    horiz_justy = 0.5
  elseif (justify(1:1) == 'R') then
    horiz_justy = 1.0
  elseif (justify(1:1) /= 'L') then
    call out_io (s_error$, r_name, 'BAD "JUSTIFY": ' // justify)
    if (global_com%exit_on_error) call err_exit
  endif
endif

end function qp_justify

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

end subroutine qp_draw_main_title

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_arrow (arrow)
!
! Subroutine to set the type of the arrows used in plotting data.
! See the quick_plot documentation for more details.
!
! Input:
!   arrow -- qp_arrow_struct: Arrow parameters.
!-

subroutine qp_set_arrow (arrow)

implicit none

type (qp_arrow_struct) arrow

qp_com%arrow = arrow

end subroutine qp_set_arrow

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
!     %type          -- Character(*): Symbol type. 
!     %height        -- Real(rp): Size of the symbol.
!     %color         -- Character(*): Symbol color.
!     %fill_pattern  -- Character(*): fill pattern.
!     %line_width    -- Integer: Line width.
!-

subroutine qp_set_symbol (symbol)

implicit none

type (qp_symbol_struct) symbol

qp_com%symbol = symbol

end subroutine qp_set_symbol

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+  
! Subroutine qp_set_arrow_attrib (color, head_size, head_type, head_angle, head_barb)
!
! Subroutine to set the arrow shape used in drawing arrows.
! See the quick_plot documentation for more details.
!
! Input:
!   color         -- character(*), optional: Arrow color.
!   head_size     -- real(rp), optional: Size of the arrow.
!   head_type     -- character(*), optional: Arrow head type: filled_arrow_head$ or outline_arrow_head$ 
!   head_angle    -- real(rp), optional: Acute angle of the arrow point in degrees.
!   head_barb     -- real(rp), optional: Fraction of triangular arrow head that is cut away from the back.
!-

subroutine qp_set_arrow_attrib (color, head_size, head_type, head_angle, head_barb)

implicit none

real(rp), optional :: head_size, head_angle, head_barb
character(*), optional :: color, head_type

!

if (present(color))      qp_com%arrow%color       = color
if (present(head_size))  qp_com%arrow%head_size   = head_size
if (present(head_type))  qp_com%arrow%head_type   = head_type
if (present(head_angle)) qp_com%arrow%head_angle  = head_angle
if (present(head_barb))  qp_com%arrow%head_barb   = head_barb

end subroutine qp_set_arrow_attrib

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
!   type         -- Character(*), optional: Symbol type. 
!   height       -- Real(rp), optional: Size of the symbol.
!   color        -- Character(*), optional: Symbol color.
!   fill_pattern -- Character(*), optional: fill pattern.
!   line_width   -- Integer, optional: Line width.
!   clip         -- Logical, optional: Clip at graph boundary?
!                     Note: This sets the line clip also.
!-

subroutine qp_set_symbol_attrib (type, height, color, fill_pattern, line_width, clip)

implicit none

integer, optional :: line_width
real(rp), optional :: height
logical, optional :: clip
character(*), optional :: type, fill_pattern, color

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
                    qp_string_to_enum(qp_com%symbol%type, 'symbol_type'), qp_com%uniform_symbol_size)
call qp_set_color_basic (qp_string_to_enum(qp_com%symbol%color, 'color'))
call qp_set_symbol_fill_basic (qp_string_to_enum(qp_com%symbol%fill_pattern, 'fill_pattern'))
call qp_set_line_width_basic (qp_com%symbol%line_width)

end subroutine qp_set_symbol_attrib

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+  
! Subroutine qp_get_arrow_attrib (arrow)
!
! Subroutine to get the arrow parameters used in plotting data.
! Use qp_set_arrow or qp_set_arrow_attrib to set arrow attributes.
! See the quick_plot documentation for more details.
!
! Output:
!   arrow -- qp_arrow_struct: Arrow parameters.
!-

subroutine qp_get_arrow_attrib (arrow)

implicit none

type (qp_arrow_struct) arrow

arrow = qp_com%arrow

end subroutine qp_get_arrow_attrib

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+  
! Subroutine qp_get_symbol_attrib (symbol)
!
! Subroutine to get the symbol parameters used in plotting data.
! Use qp_set_symbol or qp_set_symbol_attrib to set symbol attributes.
! See the quick_plot documentation for more details.
!
! Output:
!   symbol -- qp_symbol_struct:
!     %type          -- Character(*): Symbol type. 
!     %height        -- Real(rp): Size of the symbol.
!     %color         -- Character(*): Symbol color.
!     %fill_pattern  -- Character(*): fill pattern.
!     %line_width    -- Integer: Line width.
!-

subroutine qp_get_symbol_attrib (symbol)

implicit none

type (qp_symbol_struct) symbol

symbol = qp_com%symbol

end subroutine qp_get_symbol_attrib

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
!     %pattern   -- Character(*): 'solid', etc.
!     %width     -- Integer: Size of the line.
!     %color     -- Character(*): Line color.
!-

subroutine qp_set_line (who, line)

implicit none

type (qp_line_struct) line
character(*) who
character(*), parameter :: r_name = 'qp_set_line'

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
  call out_io (s_error$, r_name, 'UNKNOWN LINE "WHO": ' // who)
  if (global_com%exit_on_error) call err_exit
endif

call qp_set_line_attrib (who)

end subroutine qp_set_line

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_get_line_attrib (style, line)
!
! Subroutine to get the default line attributes.
! See the quick_plot documentation for more details.
!
! Input:
!   style   -- Character(*): 
!                 'AXIS'     -- Graph axis.
!                 'GRID'     -- Graph grid.
!                 'LEGEND'   -- Line legend.
!                 'PLOT'     -- Plot data lines.
!                 'STD'      -- Everything else.
!
! Output:
!   line    -- qp_line_struct: Attributes of a line
!     %pattern   -- Character(*): Line type.
!     %width     -- Integer: Size of the line.
!     %color     -- Character(*): Line color.
!-

subroutine qp_get_line_attrib (style, line)

implicit none

type (qp_line_struct) line
character(*) style
character(24) :: r_name = 'qp_get_line_attrib'

!

if (style == 'STD') then
  line = qp_com%std_line
elseif (style == 'GRID') then
  line = qp_com%grid_line
elseif (style == 'PLOT') then
  line = qp_com%plot_line
elseif (style == 'AXIS') then
  line = qp_com%axis_line
elseif (style == 'LEGEND') then
  line = qp_com%legend_line
else
  call out_io (s_error$, r_name, 'UNKNOWN LINE "STYLE": ' // style)
  if (global_com%exit_on_error) call err_exit
endif

end subroutine qp_get_line_attrib

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_line_attrib (style, width, color, pattern, clip)
!        
! Subroutine to set the default line attributes for different parts of the drawing.
! Also: Sets the current line attributes to the style indicated by "style".
! See the quick_plot documentation for more details.
!
! Input:
!   style       -- Character(*), optional: Part of drawing default is set for.
!                 'PLOT'     -- Plot data lines.
!                 'GRID'     -- Graph grid.
!                 'AXIS'     -- Graph axis.
!                 'LEGEND    -- Line legend.
!                 'STD'      -- Everything else. Default.
!   pattern   -- Character(*), optional: Line type.
!   width     -- Integer, optional: Size of the line.
!   color     -- Character(*), optional: Line color.
!   clip      -- Logical, optional: Clip at graph boundary?
!                  Note: This sets the symbol clip also.
!-

subroutine qp_set_line_attrib (style, width, color, pattern, clip)

implicit none

type (qp_line_struct), pointer :: this

integer, optional :: width
logical, optional :: clip

character(*), optional :: style, pattern, color
character(*), parameter :: r_name = 'qp_set_line_attrib'


!

if (.not. present(style)) then
  this => qp_com%std_line
elseif (style == 'STD') then
  this => qp_com%std_line
elseif (style == 'GRID') then
  this => qp_com%grid_line
elseif (style == 'PLOT') then
  this => qp_com%plot_line
elseif (style == 'AXIS') then
  this => qp_com%axis_line
elseif (style == 'LEGEND') then
  this => qp_com%legend_line
else
  call out_io (s_error$, r_name, 'UNKNOWN LINE "STYLE": ' // style)
  if (global_com%exit_on_error) call err_exit
endif

if (present(width)) this%width = width
if (present(color)) this%color = color
if (present(pattern)) this%pattern = pattern

call qp_set_clip (clip)
call qp_set_color_basic (qp_string_to_enum(this%color, 'color'))
call qp_set_line_width_basic (this%width)
call qp_set_line_pattern_basic (qp_string_to_enum(this%pattern, 'line_pattern'))

end subroutine qp_set_line_attrib

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

end subroutine qp_set_clip

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

end subroutine qp_subset_box

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_get_text_attrib (who, height, color, background, &
!                                       uniform_spacing, spacing_factor)
!
! Routine to get the default text attributes for titles, legends, etc.
! See the quick_plot documentation for more details.
! Note: The background color index is common to all types of text.
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
!   color       -- Character(*), optional: Color index.
!   background  -- Character(*), optional: Background color index.
!   uniform_spacing -- Logical, optional: If T then the distance between 
!                      characters is uniform.
!   spacing_factor  -- Real(rp), optional: Spacing factor for the 
!             uniform_spacing option. This is set globally for all "who".
!-

subroutine qp_get_text_attrib (who, height, color, background, &
                                             uniform_spacing, spacing_factor)

implicit none

real(rp), optional :: height, spacing_factor
logical, optional :: uniform_spacing
character(*) who
character(*), optional :: color, background
character(24) :: r_name = 'qp_set_text_attrib'

!

if (who == "MAIN_TITLE") then
  call qp_get_this_text_attrib (qp_com%main_title)
elseif (who == "GRAPH_TITLE") then
  call qp_get_this_text_attrib (qp_com%graph_title)
elseif (who == "LEGEND") then
  call qp_get_this_text_attrib (qp_com%legend)
elseif (who == "TEXT") then
  call qp_get_this_text_attrib (qp_com%text)
elseif (who == "AXIS_NUMBERS") then
  call qp_get_this_text_attrib (qp_com%axis_number)
elseif (who == "AXIS_LABEL") then
  call qp_get_this_text_attrib (qp_com%axis_label)
else
  call out_io (s_error$, r_name, 'BAD "WHO": "' // trim(who) // '"' )
  if (global_com%exit_on_error) call err_exit
endif

!----------------------------------------------------------------
contains

subroutine qp_get_this_text_attrib (this_text)

type (qp_text_struct) this_text
real(rp) text_height

!

if (present(spacing_factor)) spacing_factor = qp_com%text_spacing_factor 

if (present(height)) then
  height = this_text%height
endif

if (present(color)) then
  color = this_text%color
endif

if (present(background)) then
  background = qp_com%text_background
endif

if (present(uniform_spacing)) then
  uniform_spacing = this_text%uniform_spacing
endif

end subroutine qp_get_this_text_attrib

end subroutine qp_get_text_attrib

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
!            "MAIN_TITLE"    qp_draw_main_title   Title at top of page.
!            "GRAPH_TITLE"   qp_draw_graph_title  Title above a graph.
!            "LEGEND"        qp_draw_text_legend  Legend.
!            "LEGEND"        qp_draw_label_legend Legend.
!            "AXIS_NUMBERS"  qp_draw_graph        Axes Numbers.
!            "AXIS_LABEL"    qp_draw_graph        Axis label.
!            "TEXT"          qp_draw_text         Everything else.
!   height      -- Real(rp), optional: Character height.
!   color       -- Character(*), optional: Color index.
!   background  -- Character(*), optional: Background color index.
!   uniform_spacing -- Logical, optional: If T then the distance between 
!                      characters is uniform.
!   spacing_factor  -- Real(rp), optional: Spacing factor for the 
!             uniform_spacing option. This is set globally for all "who".
!-

subroutine qp_set_text_attrib (who, height, color, background, &
                                             uniform_spacing, spacing_factor)

implicit none

real(rp), optional :: height, spacing_factor
logical, optional :: uniform_spacing
character(*) who
character(*), optional :: color, background
character(*), parameter :: r_name = 'qp_set_text_attrib'

!

if (who == "MAIN_TITLE") then
  call qp_set_this_text_attrib (qp_com%main_title)
elseif (who == "GRAPH_TITLE") then
  call qp_set_this_text_attrib (qp_com%graph_title)
elseif (who == "LEGEND") then
  call qp_set_this_text_attrib (qp_com%legend)
elseif (who == "TEXT") then
  call qp_set_this_text_attrib (qp_com%text)
elseif (who == "AXIS_NUMBERS") then
  call qp_set_this_text_attrib (qp_com%axis_number)
elseif (who == "AXIS_LABEL") then
  call qp_set_this_text_attrib (qp_com%axis_label)
else
  call out_io (s_error$, r_name, 'BAD "WHO": "' // trim(who) // '"' )
  if (global_com%exit_on_error) call err_exit
endif

!----------------------------------------------------------------
contains

subroutine qp_set_this_text_attrib (this_text)

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
call qp_set_color_basic (qp_string_to_enum(this_text%color, 'color'))
call qp_set_text_background_color_basic (qp_string_to_enum(qp_com%text_background, 'color'))
call qp_set_line_width_basic (1)

qp_com%this_text = this_text

end subroutine qp_set_this_text_attrib 

end subroutine qp_set_text_attrib

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
                                                                      
type (qp_axis_struct), pointer :: ax1, ax2

real(rp) dx0, dum, x1, y1, dy, x0, y0, y_pos, dy1, dy2
real(rp) dy11, dy22, x11, tick_width, x0_tick, r, d, x1_inch

integer i, j, m_div, who_sign, divisions, m_min, m_max, d_div

character(*) who
character(16) justify, str
character(*), parameter :: r_name = 'qp_draw_x_axis'

! save state

call qp_save_state (.true.)
call qp_set_clip (.false.)     ! no clipping of axis

! ax1 and ax2 are the same except when x2_mirrors_x is True.

call qp_use_axis (x = who)
ax2 => qp_com%plot%xx

if (who == 'X') then
  who_sign = +1
elseif (who == 'X2') then
  if (qp_com%plot%x2_mirrors_x) call qp_use_axis (x = 'X')
  who_sign = -1
else
  call out_io (s_error$, r_name, 'BAD AXIS NAME: ' // who)
endif

ax1 => qp_com%plot%xx

! the axis line itself

call qp_set_line_attrib ('AXIS')
call qp_draw_polyline_no_set ([0.0_rp, 1.0_rp], [y_pos, y_pos], '%GRAPH')  

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

if (ax1%number_side*who_sign == +1) then
  justify = 'CB'
else
  justify = 'CT'
endif

call qp_to_inch_rel (1.0_rp, y_pos, x1_inch, y0, '%GRAPH')

!

if (ax1%type == 'LOG') then
  if (ax1%tick_min <= 0) then
    call out_io (s_abort$, r_name, 'NEGATIVE VALUES ENCOUNTERED WHEN DRAWING LOG X-AXIS!')
    if (global_com%exit_on_error) call err_exit
    call qp_restore_state
    return
  endif
  call qp_to_inch_rel (ax1%tick_max/ax1%tick_min, 0.0_rp, tick_width, dum, 'DATA')
  call qp_to_inch_rel (ax1%tick_min/ax1%min, 0.0_rp, x0_tick, dum, 'DATA')
  d_div = nint(log10(ax1%dtick))
  m_max = nint(log10(ax1%tick_max))
  m_min = nint(log10(ax1%tick_min))
  divisions = (m_max - m_min) / d_div
  dx0 = tick_width / divisions
  ! minor ticks
  ! If d_div > 1 then not enough room to subdivide a decade.
  if (d_div == 1) then
    do i = -1, divisions
      do j = 2, 9
        x11 = (i + log10(real(j))) * dx0 + x0_tick
        if (x11 < 0 .or. x11 > x1_inch) cycle
        call qp_draw_polyline_no_set ([x11, x11], [y0+dy11, y0+dy22], 'INCH')
      enddo
    enddo
  elseif (d_div <= max(ax1%minor_div, ax1%minor_div_max)) then
    do i = m_min, m_max-1
      do j = 0, d_div-1
        x11 = (i + real(j) / d_div) * dx0
        call qp_draw_polyline_no_set ([x11, x11], [y0+dy11, y0+dy22], 'INCH')
      enddo
    enddo
  endif

else ! Linear scale
  call qp_to_inch_rel (ax1%tick_max-ax1%tick_min, 0.0_rp, tick_width, dum, 'DATA')
  call qp_to_inch_rel (ax1%tick_min-ax1%min, 0.0_rp, x0_tick, dum, 'DATA')
  divisions = ax1%major_div
  if (divisions == 0) then
    dx0 = 0
  else
    dx0 = tick_width / divisions
  endif

  ! minor ticks
  if (ax1%minor_div == 0) then
    call qp_calc_minor_div (ax1%dtick, ax1%minor_div_max, m_div)
  else
    m_div = ax1%minor_div
  endif
  
  do i = -1, divisions
    do j = 1, m_div - 1
      r = i + real(j) / m_div
      d = ax1%tick_min + r * ax1%dtick
      if (d < ax1%min .or. d > ax1%max) cycle
      x11 = r * dx0 + x0_tick
      call qp_draw_polyline_no_set ([x11, x11], [y0+dy11, y0+dy22], 'INCH')
    enddo
  enddo
endif


! Axis numbers

call qp_set_text_attrib ('AXIS_NUMBERS')
do i = 0, divisions

  if (ax1%type == 'LOG') then
    x1 = i*dx0 + x0_tick
  else
    x1 = i*dx0 + x0_tick
  endif
  y1 = y0 + ax1%number_side * who_sign * ax1%number_offset
 
  ! major tick
  call qp_draw_polyline_no_set ([x1, x1], [y0+dy1, y0+dy2], 'INCH')

  ! numbers
  if (.not. ax2%draw_numbers) cycle

  call qp_to_axis_number_text (ax1, i, str)

  if (i == 0) then
    call qp_draw_text_no_set (str, x1, y1, 'INCH', justify)
  elseif (i == ax1%major_div) then
    call qp_draw_text_no_set (str, x1, y1, 'INCH', justify)
  else
    call qp_draw_text_no_set (str, x1, y1, 'INCH', justify)
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

end subroutine qp_draw_x_axis

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

type (qp_axis_struct), pointer :: ax1, ax2

real(rp) x1, y1, dx, x0, y0, x_pos, dx1, dx2, dx11, dx22, y1_inch
real(rp) number_len, dy0, y11, tick_width, ax_len, dum, y0_tick, r, d

integer i, j, m_div, who_sign, number_side, divisions, d_div, m_max, m_min

character(*) who
character(2) justify
character(16) str
character(*), parameter :: r_name = 'qp_draw_y_axis'

! save state
  
call qp_set_clip (.false.)     ! no clipping of axis
call qp_save_state (.true.)

! ax1 and ax2 are the same except when x2_mirrors_x is True.

call qp_use_axis (y = who)
ax2 => qp_com%plot%yy

if (who == 'Y') then
  who_sign = +1
elseif (who == 'Y2') then
  if (qp_com%plot%y2_mirrors_y) call qp_use_axis (y = 'Y')
  who_sign = -1
else
  call out_io (s_error$, r_name, 'BAD AXIS NAME: ' // who)
endif

ax1 => qp_com%plot%yy   ! Not necessarily the same as ax2
number_side = ax1%number_side * who_sign

! draw axis line itself

call qp_set_line_attrib ('AXIS')
call qp_draw_polyline_no_set ([x_pos, x_pos], [0.0_rp, 1.0_rp], '%GRAPH')

! major and minor divisions calc

if (ax2%tick_side == 0) then
  dx1  =  ax1%major_tick_len
  dx2  = -ax1%major_tick_len
  dx11 =  ax1%minor_tick_len
  dx22 = -ax1%minor_tick_len
else
  dx1  = 0
  dx2  = ax1%major_tick_len * ax2%tick_side * who_sign
  dx11 = 0
  dx22 = ax1%minor_tick_len * ax2%tick_side * who_sign
endif

call qp_to_inch_rel (x_pos, 1.0_rp, x0, y1_inch, '%GRAPH')

if (ax1%type == 'LOG') then
  if (ax1%tick_min <= 0 .or. ax1%min <= 0) then
    call out_io (s_abort$, r_name, 'NEGATIVE VALUES ENCOUNTERED WHEN DRAWING LOG Y-AXIS!')
    if (global_com%exit_on_error) call err_exit
    call qp_restore_state
    return
  endif
  call qp_to_inch_rel (0.0_rp, ax1%tick_max/ax1%tick_min, dum, tick_width, 'DATA')
  call qp_to_inch_rel (0.0_rp, ax1%tick_min/ax1%min, dum, y0_tick, 'DATA')
  divisions = max(nint(log10(ax1%tick_max)) - nint(log10(ax1%tick_min)), 1)
  d_div = nint(log10(ax1%dtick))
  m_max = nint(log10(ax1%tick_max))
  m_min = nint(log10(ax1%tick_min))
  divisions = (m_max - m_min) / d_div
  dy0 = tick_width / divisions
  ! If d_div > 1 then not enough room to subdivide a decade.
  if (d_div == 1) then
    do i = -1, divisions
      do j = 2, 9
        y11 = (i + log10(real(j))) * dy0 + y0_tick
        if (y11 < 0 .or. y11 > y1_inch) cycle
        call qp_draw_polyline_no_set ([x0+dx11, x0+dx22], [y11, y11], 'INCH')
      enddo
    enddo
  elseif (d_div <= max(ax1%minor_div, ax1%minor_div_max)) then
    do i = m_min, m_max-1
      do j = 0, d_div-1
        y11 = (i + real(j) / d_div) * dy0
        call qp_draw_polyline_no_set ([x0+dx11, x0+dx22], [y11, y11], 'INCH')
      enddo
    enddo
  endif

else ! Linear
  call qp_to_inch_rel (0.0_rp, ax1%tick_max-ax1%tick_min, dum, tick_width, 'DATA')
  call qp_to_inch_rel (0.0_rp, ax1%tick_min-ax1%min, dum, y0_tick, 'DATA')
  divisions = ax1%major_div
  if (divisions == 0) then
    dy0 = 0
  else
    dy0 = tick_width / divisions
  endif

  ! Minor ticks
  if (ax1%minor_div == 0) then
    call qp_calc_minor_div (ax1%dtick, ax1%minor_div_max, m_div)
  else
    m_div = ax1%minor_div
  endif

  do i = -1, divisions
    do j = 1, m_div - 1
      r = i + real(j) / m_div
      d = ax1%tick_min + r * ax1%dtick
      if (d < ax1%min .or. d > ax1%max) cycle
      y11 = r * dy0 + y0_tick
      call qp_draw_polyline_no_set ([x0+dx11, x0+dx22], [y11, y11], 'INCH')
    enddo
  enddo
endif

! draw axis numbers and major tick

call qp_set_text_attrib ('AXIS_NUMBERS')

number_len = 0   ! length in inches

do i = 0, divisions
  x1 = x0 + number_side * ax1%number_offset
  y1 = i*dy0 + y0_tick

  ! major tick
  call qp_draw_polyline_no_set ([x0+dx1, x0+dx2], [y1, y1], 'INCH')

  ! Number
  if (.not. ax2%draw_numbers) cycle

  call qp_to_axis_number_text (ax1, i, str)

  if (number_side == +1) then
    justify = 'L'
  else
    justify = 'R'
  endif

  justify(2:2) = 'C'

  call qp_draw_text_no_set (str, x1, y1, 'INCH', justify)

  number_len = max (number_len, qp_text_len(str))
enddo

! draw label

if (ax2%draw_label .and. ax2%label /= '') then
  if (.false.) print *, 'Here'   ! To get around ifort/production/openmp bomb when running Tao!
  call qp_to_inch_rel (x_pos, 0.5_rp, x0, y0, '%GRAPH')
  dx = number_side * (ax1%number_offset + ax1%label_offset + number_len)
  call qp_set_text_attrib ('AXIS_LABEL', color = ax2%label_color)
  call qp_draw_text_no_set (ax2%label, x0+dx, y0, 'INCH', 'CB', angle = -90.0_rp * number_side)
endif

!

call qp_restore_state
                   
end subroutine qp_draw_y_axis

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_to_axis_number_text (axis, ix_n, text)
!
! Subroutine to form the text string for an axis number.
!
! Input:
!   axis      -- qp_axis_struct:
!   ix_n      -- Integer: Index of particular number.
!
! Output:
!   text      -- Character(*): Character string.
!-

subroutine qp_to_axis_number_text (axis, ix_n, text)

implicit none

type (qp_axis_struct) axis

integer i, ix_n, ie, n_char, n_log_min, n_log_max, n_log_delta
integer n_zero_crit, ix, n_log, p

real(rp) val, v, effective_zero

character(*) text
character(20) fmt

logical too_small, too_large

! LOG scale

if (axis%type == 'LOG') then
  ie =  nint(log10(axis%tick_min) + ix_n * log10(axis%dtick))
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

val = axis%tick_min + ix_n * axis%dtick
effective_zero = max(abs(axis%tick_min), abs(axis%tick_max)) / 10.0_rp**qp_com%max_digits

! If the number is essentially zero then life is simple

if (abs(val) < effective_zero) then
  text = '0'
  return
endif

! n_char is the number of characters we need to draw

n_log_min = 1000
n_log_max = -1000
do i = 0, axis%major_div
  v = axis%tick_min + i * axis%dtick
  if (abs(v) < effective_zero) cycle  ! Ignore zero.
  n_log = floor(log10(abs(v)+1d-30) + 0.0001)
  n_log_min   = min (n_log_min, n_log)
  n_log_max   = max (n_log_max, n_log)
enddo

n_char = max(0, n_log_max) + 1
n_char = n_char + max(axis%places, 0)
if (axis%tick_max < 0 .or. axis%tick_min < 0) n_char = n_char + 1  ! for negative sign
if (axis%places > 0) n_char = n_char + 1  ! add 1 character for the decimal place itself

! Special case: Switch to scientific notation if the number is too big or too small.

too_large = .false.
too_small = .false.
n_zero_crit = qp_com%max_axis_zero_digits 

n_log_delta = floor(log10(axis%dtick+1d-30) + 0.0001)
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

end subroutine qp_to_axis_number_text

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

t_len = qp_text_len_basic (text)

end function qp_text_len

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

end subroutine qp_set_graph_attrib

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_grid ()
!
! Subroutine to draw a grid on the current graph.
!-

subroutine qp_draw_grid ()

implicit none

type (qp_axis_struct), pointer :: ax, ay
real(rp) r0, z(2), width_x, width_y
integer i, divisions

!

call qp_save_state (.true.)
call qp_set_line_attrib ('GRID')

ax => qp_com%plot%xx
ay => qp_com%plot%yy

width_x = ax%max - ax%min
width_y = ay%max - ay%min

! horizontal lines

if (ay%type == 'LOG') then
  if (ay%tick_min > 0) then
    divisions = nint((log10(ay%tick_max) - log10(ay%tick_min)) / log10(ay%dtick))
    do i = 0, divisions
      z = ay%tick_min * ay%dtick**i
      call qp_draw_polyline_no_set ([ax%min, ax%max], z, 'DATA')
    enddo
  endif

else
  do i = 0, ay%major_div
    z = ay%tick_min + i * ay%dtick
    call qp_draw_polyline_no_set ([ax%min, ax%max], z, 'DATA')
  enddo
endif

! vertical lines

if (ax%type == 'LOG') then
  if (ax%tick_min > 0) then
    divisions = nint((log10(ax%tick_max) - log10(ax%tick_min)) / log10(ax%dtick))
    do i = 0, divisions
      z = ax%tick_min * ax%dtick**i
      call qp_draw_polyline_no_set (z, [ay%min, ay%max], 'DATA')
    enddo
  endif

else
  do i = 0, ax%major_div
    z = ax%tick_min + i * ax%dtick
    call qp_draw_polyline_no_set (z, [ay%min, ay%max], 'DATA')
  enddo
endif

!

call qp_restore_state 

end subroutine qp_draw_grid 

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

! Sinple case

if (delta == 0) then
  divisions = 1
  return
endif

! scale delta so it is in the range of [10, 100)

log_del = log10 (abs(delta) * 1.000000001_dp)
idel = nint(abs(delta) / 10d0**(floor(log_del)-1))

! First look for a division that gives a width that is a multiple of 2, 5 or 10.
! A division of 1 is not acceptable in this first step.

do divisions = div_max, 2, -1
  if (idel/divisions == 5) return
  if (mod(idel/divisions + 5, 10) == 0) cycle  ! Reject 15, 25, 35, 45, etc.
  if (mod(idel, 5 * divisions) == 0) return
enddo

do divisions = div_max, 2, -1
  if (mod(idel, 10 * divisions) == 0) return
enddo

do divisions = div_max, 1, -1
  if (mod(idel, 2 * divisions) == 0) return
enddo

! Now look for anything that divides evenly into idel.

do divisions = div_max, 1, -1
  if (mod(idel, divisions) == 0) return
enddo

end subroutine qp_calc_minor_div

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_get_parameters (text_scale, default_draw_units, &
!                           default_set_units, default_axis_slop_factor)
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



end subroutine qp_get_parameters

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

end subroutine qp_set_parameters

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

end function qp_text_height_to_inches

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
real(rp), allocatable :: xyz(:)
real(rp) :: x1(i_del), y1(i_del), z1(i_del), t1(i_del)

integer i, j, ix, ios
integer iu, i_size
integer, optional :: ix_col, iy_col, iz_col, it_col

logical err_flag, good_x, good_y, good_z, good_t

character(140) line, line_in
character(*), parameter :: r_name = 'qp_read_data'

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

end subroutine load_data

end subroutine qp_read_data

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_eliminate_xy_distortion (axis_to_scale)
!
! This subroutine will vary the X or Y active axes scale so that the distance between ticks and the
! conversion between data units and inches on the page is the same for the x and y axes.
! In other words, this routine will make sure that a square in data units looks like a square when drawn.
! This routine is useful in drawing such things as maps.
!
! In varying an axis scale, the [min, max] range is never decreased.
!
! Note: Whether X or X2 axis is varied or whether the Y or Y2 axis is varied depends upon which axes are
! active as set by the qp_use_axis routine.
!
! Input:
!   axis_to_scale -- Logical, optional:
!                         'XY'   -> Vary the active x or active y axis (default).
!                         'X'    -> Vary the active x axis.
!                         'Y'    -> Vary the active y axis.
!-

subroutine qp_eliminate_xy_distortion (axis_to_scale)

implicit none

type (qp_axis_struct), pointer :: x_axis, y_axis
real(rp) x_unit_len, y_unit_len, x_graph_len, y_graph_len, x_data_len, y_data_len

character(*), optional :: axis_to_scale
character(8) ax_to_scale

! 

x_axis => qp_com%plot%xx
y_axis => qp_com%plot%yy

ax_to_scale = 'XY'
if (present(axis_to_scale)) ax_to_scale = axis_to_scale

call qp_to_inch_rel (1.0_rp, 1.0_rp, x_unit_len, y_unit_len)  ! Length in inches for a data delta of 1.
x_data_len = x_axis%max - x_axis%min
y_data_len = y_axis%max - y_axis%min

! Adjust one axis scale

x_graph_len = qp_com%graph%x2 - qp_com%graph%x1
y_graph_len = qp_com%graph%y2 - qp_com%graph%y1

if (ax_to_scale == 'X' .or. (ax_to_scale == 'XY' .and. x_data_len/x_graph_len < y_data_len/y_graph_len)) then  
  call axis_scale (x_axis, y_axis, y_unit_len/x_unit_len, y_graph_len/x_graph_len) ! Vary x-axis scale
else                                                                
  call axis_scale (y_axis, x_axis, x_unit_len/y_unit_len, x_graph_len/y_graph_len) ! Vary y-axis scale
endif

! Now adjust the margins

call qp_set_graph_limits
call qp_to_inch_rel (1.0_rp, 1.0_rp, x_unit_len, y_unit_len)

if (ax_to_scale == 'X' .or. (ax_to_scale == 'XY' .and. x_unit_len > y_unit_len)) then  
  ! Shrink x-axis margins
  call margin_scale (y_unit_len / x_unit_len, qp_com%graph%x2 - qp_com%graph%x1, qp_com%margin%x1, qp_com%margin%x2) 
else                                                                             
  ! Shrink y-axis margins
  call margin_scale (x_unit_len / y_unit_len, qp_com%graph%y2 - qp_com%graph%y1, qp_com%margin%y1, qp_com%margin%y2) 
endif

call qp_set_graph_limits

!----------------------------------------------------
contains

subroutine axis_scale (axis_z, axis_t, tz_scale_ratio, tz_graph_ratio)

type (qp_axis_struct) axis_z, axis_t

real(rp) tz_scale_ratio, tz_graph_ratio, dz

! axis_z is the axis to vary.
! axis_t is not changed.

axis_z%dtick = axis_t%dtick
axis_z%minor_div = axis_t%minor_div
axis_z%places = axis_t%places

if (axis_z%max - axis_z%min < axis_z%dtick) then
  axis_z%tick_max = axis_t%dtick * ceiling(axis_z%max / axis_t%dtick)
  axis_z%tick_min = axis_t%dtick * floor(axis_z%min / axis_t%dtick)
  if (axis_z%tick_max == axis_z%tick_min) then
    axis_z%tick_max = axis_z%tick_max + axis_t%dtick
    axis_z%tick_min = axis_z%tick_min - axis_t%dtick
  endif
else
  axis_z%tick_max = axis_t%dtick * ceiling(axis_z%max / axis_t%dtick - 0.3)
  axis_z%tick_min = axis_t%dtick * floor(axis_z%min / axis_t%dtick + 0.3)
endif

axis_z%major_div = nint((axis_z%tick_max - axis_z%tick_min) / axis_z%dtick)
axis_z%max = max(axis_z%max, axis_z%tick_max)
axis_z%min = min(axis_z%min, axis_z%tick_min)

end subroutine axis_scale

!----------------------------------------------------
! contains

subroutine margin_scale (tz_scale_ratio, graph_width, margin_z1, margin_z2)

real(rp) rd, tz_scale_ratio, graph_width, margin_z1, margin_z2

! Scale margins

rd = (1 - tz_scale_ratio) * graph_width / 2
margin_z1 = margin_z1 + rd
margin_z2 = margin_z2 + rd

end subroutine margin_scale

end subroutine qp_eliminate_xy_distortion

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function qp_continuous_color (real_color) result (integer_color)
! 
! Maps real colors from 0.0_rp -- 1.0_rp to integers between[17, huge(integer)]
!
! Input:
!   real_color    -- real(rp): between 0.0_rp -- 1.0_rp
!                              Values outside this range are mapped to 
!                              integer colors 1 and 0
!
! Output:
!   integer_color -- integer: between[17, huge(integer)]
!-
function qp_continuous_color (real_color) result (integer_color)

implicit none

real(rp) :: real_color
integer :: integer_color
if (real_color > 1.0_rp) then
  integer_color = 1
else if (real_color < 0.0_rp) then
  integer_color = 0
else
  integer_color = nint(real_color*(huge(integer_color)-17) + 17)
endif

end function qp_continuous_color

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function qp_continuous_color_inverse (integer_color) result (real_color)
! 
! Inverse of function qp_continuous_color
!
!-
function qp_continuous_color_inverse (integer_color) result (real_color)

implicit none

real(rp) :: real_color
integer :: integer_color

if (integer_color < 17) then
  real_color = 0.0_rp
else
  real_color = (integer_color - 17)/ (1.0_rp*(huge(integer_color) - 17) )
endif

end function qp_continuous_color_inverse

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_end()
!
! Cleanup routine at the end of plotting.
!-

subroutine qp_end()

call qp_end_basic()

end subroutine qp_end

end module
