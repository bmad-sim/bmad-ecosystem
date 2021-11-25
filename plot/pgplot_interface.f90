                                    
!+
! Module pgplot_interface
!
! Module to interface between pgplot and quick_plot routines.
! This module isolates quick_plot from pgplot.
!
! To tell the difference between quick_plot and pgplot: 
!   QUICK_PLOT subroutines start with "QP_"
!   PGPLOT routines start with "PG" and do not have any "_".
!
! Note: PGPLOT uses real(4) units, QUICK_PLOT tries to interface between
!       real(4) and real(rp). Watch out for mixed units in subroutines.
!
! The correspondence between PAGE, BOX, and GRAPH and PGPLOT is:
!   QUICK_PLOT    PGPLOT
!   ----------    ------
!   PAGE          VIEW SURFACE 
!   BOX           No corresponding entity.
!   GRAPH         VIEWPORT and WINDOW
!
! Essentially the VIEWPORT is the region outside of which lines and symbols
! will be clipped (if clipping is turned on) and the WINDOW defines the 
! plot area. I'm not sure why PGPLOT makes a distinction, but VIEWPORT and 
! WINDOW always are the same region.
!-

module pgplot_interface

use output_mod
use quick_plot_struct

! This #if def wraps the entire module

#if defined (CESR_PGPLOT)

type pg_interface_struct
character(16) page_type
character(100) plot_file
integer :: i_chan = -1
real qp_to_pg_text_height_factor
real page_scale   ! scaling for entire page
end type

type (pg_interface_struct), pointer, save, private :: pg_com
type (pg_interface_struct), save, target, private :: pg_interface_save_com(0:20)
integer, save, private :: i_save = 0
logical, save, private :: pg_allow_flush = .true.

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_graph_position_basic (x1, x2, y1, y2)
!
! Subroutine to set the position of a graph.
! Units are inches from lower left of page.
!
! Input:
!   x1, y1 -- Real(rp): Bottom left corner of graph.
!   x2, y2 -- Real(rp): Upper right corner of graph.
!-

subroutine qp_set_graph_position_basic (x1, x2, y1, y2)
implicit none
real(rp) x1, x2, y1, y2, f
real xx1, xx2, yy1, yy2
! Single precision roundoff can make things equal which PGPLOT does not like
f = pg_com%page_scale
xx1 = f * x1
xx2 = f * x2
if (xx2 == xx1) xx2 = xx2 + (1 + 10.0**(-7))
yy1 = f * y1
yy2 = f * y2
if (yy2 == yy1) yy2 = yy2 * (1 + 10.0**(-7))
call pgvsiz (xx1, xx2, yy1, yy2)
call pgswin (xx1, xx2, yy1, yy2)
end subroutine qp_set_graph_position_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_symbol_size_basic (height, symbol_type, uniform_size)
!
! Subroutine to set the symbol_size
!
! Input:
!   height       -- Real(rp): Symbol height.
!   symbol_type  -- Integer: Symbol type.
!   uniform_size -- Logical: Make all symbols the save size for constant height.
!-

subroutine qp_set_symbol_size_basic (height, symbol_type, uniform_size)

implicit none

real(rp) height
real h, f

integer symbol_type

logical uniform_size

! The PGPLOT symbol set does not have a constant symbol size.
! This generally does not look nice so renormalize to get a consistant size.
! This excludes the set of circles with different sizes.

f = pg_com%page_scale
h = height * pg_com%qp_to_pg_text_height_factor

if (uniform_size) then

  select case (symbol_type)
  case (dot_sym$)
    if (pg_com%page_type(1:3) /= 'GIF') h = h * 2.0       ! I like bigger dots
  case (square_filled_sym$)
    h = h * 1.56
  case (circle_filled_sym$)
    h = h * 1.60
  case (star5_filled_sym$)
    h = h * 1.30
  case (square_concave_sym$)
    h = h * 0.73
  end select

  if (pg_com%page_type == 'X' .or. pg_com%page_type == 'TK' .or. pg_com%page_type == 'QT') then
    select case (symbol_type)
    case (circle_sym$)
      h = h * 0.89
    case (circle_plus_sym$)
      h = h * 0.55
    case (circle_dot_sym$)
      h = h * 0.59
    case (triangle_filled_sym$)
      h = h * 1.22
    end select
  else
    if (symbol_type == triangle_filled_sym$) h = h * 1.03 
  endif

endif

call pgsch (f * h)   ! set symbol size

end subroutine qp_set_symbol_size_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_symbol_fill_basic (fill)
!
! Subroutine to set the symbol fill pattern.
!
! Input:
!   fill -- Integer: fill pattern.
!-

subroutine qp_set_symbol_fill_basic (fill)
implicit none
integer fill
call pgsfs (fill)       ! set fill
end subroutine qp_set_symbol_fill_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_line_width_basic (line_width)
!
! Subroutine to set the line width.
!
! Input:
!   line_width -- Integer: Line width.
!-

subroutine qp_set_line_width_basic (line_width)
implicit none
integer line_width
call pgslw (line_width) ! set line width
end subroutine qp_set_line_width_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_line_pattern_basic (line_pattern)
!
! Subroutine to set the line type.
!
! Input:
!   line_pattern -- Integer: Line type.
!-

subroutine qp_set_line_pattern_basic (line_pattern)
implicit none
integer line_pattern
call pgsls  (line_pattern)       ! Set line type
end subroutine qp_set_line_pattern_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_clip_basic (clip)
!
! Subroutine to set the clipping state.
! Note: This affects both lines and symbols.
!
! Input:
!   clip    -- Logical: Clip at graph boundary?
!-

subroutine qp_set_clip_basic (clip) 
implicit none
logical :: clip
if (clip) then
  call pgsclp (1)     ! Clip on
else
  call pgsclp (0)     ! Clip off
endif
end subroutine qp_set_clip_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_char_size_basic (height)
!
! Subroutine to set the character size.
!
! Input:
!   height -- Integer: Character size.
!-

subroutine qp_set_char_size_basic (height)
implicit none
real(rp) height, f
f = pg_com%page_scale
call pgsch(real(f * height * pg_com%qp_to_pg_text_height_factor))
end subroutine qp_set_char_size_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_text_background_color_basic (color)
!
! Subroutine to set the character text background color.
!
! Input:
!   color -- Integer: Color index.
!-

subroutine qp_set_text_background_color_basic (color)
implicit none
integer color
if (color < 0)  then
  call pgstbg(color)
else
  call qp_set_color_basic (color, set_background=.true.)
endif
end subroutine qp_set_text_background_color_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function qp_text_len_basic (text) result (t_len)
!
! Function to find the length of a text string.
!
! input:
!   text -- Character(*): Text string.
!
! Output:
!   t_len -- Real(rp): Length of text in inches.
!-

function qp_text_len_basic (text) result (t_len)

implicit none

real(rp) t_len
real tl, dum, f
character(*) text

!

f = pg_com%page_scale
call pglen (1, trim(text), tl, dum)
t_len = tl / f

end function

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_text_basic (text, len_text, x0, y0, angle, justify)
!
! Subroutine to draw text.
!
! Input:
!   text     -- Character(*): Text to draw.
!   len_text -- Integer: Length of text.
!   x0, y0   -- Real(rp): Position of text in inches from bottom left of page.
!   angle    -- Real(rp): Rotation angle of text.
!   justify  -- Real(rp): Left/Right justification.
!-

subroutine qp_draw_text_basic (text, len_text, x0, y0, angle, justify)
implicit none
character(*) text
integer len_text
real(rp) x0, y0, angle, justify, f
f = pg_com%page_scale
call pgptxt (real(f*x0), real(f*y0), real(angle), real(justify), trim(text))
end subroutine qp_draw_text_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_arrow_basic (r1, r2, arrow)
!
! Subroutine to draw an arrow.
!
! Input:
!   r1(2)  -- real(rp): tail (x, y) position in inches from bottom, left page edge.
!   r2(2)  -- real(rp): head (x, y) position in inches from bottom left page edge.
!   arrow  -- qp_arrow_struct: Arrow parameters.
!-

subroutine qp_draw_arrow_basic (r1, r2, arrow)

implicit none
real(rp) r1(2), r2(2), f
real char_height
integer color
type (qp_arrow_struct) arrow

call pgqch (char_height)   ! Get current character height
call pgqci (color)         ! Get current color

call pgsah (arrow%head_type, real(arrow%head_angle), real(arrow%head_barb))  ! Set arrow parameters
call pgsch (real(arrow%head_size))                         ! Set Arrow head size
call pgsci (qp_string_to_enum(arrow%color, 'color'))       ! Set Arrow color

f = pg_com%page_scale
call pgarro (real(f*r1(1)), real(f*r1(2)), real(f*r2(1)), real(f*r2(2)))  ! Draw arrow.

call pgsch (char_height)  ! Restore old character height.
call pgsci (color)        ! Restore old color.

end subroutine qp_draw_arrow_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_symbol_basic (x, y, symbol)
!
! Subroutine to draw a symbol.
!
! Input:
!   x  -- Real(rp): X-position in inches from left page edge.
!   y  -- Real(rp): Y-position in inches from bottom page edge.
!   symbol -- Integer: symbol index.
!-

subroutine qp_draw_symbol_basic (x, y, symbol)
implicit none
real(rp) x, y, f
integer symbol
f = pg_com%page_scale
call pgpt1 (real(f*x), real(f*y), symbol)
end subroutine qp_draw_symbol_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_save_state_basic()
!
! Subroutine to save the print state.
!-

subroutine qp_save_state_basic()
if (pg_allow_flush) call pgbbuf     ! buffer commands
end subroutine qp_save_state_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_restore_state_basic()
!
! Subroutine to restore the print state.
!
! Input:
!   flush   -- logical: Flush the plot buffer?
!-

subroutine qp_restore_state_basic()
if (pg_allow_flush) call pgebuf        ! flush buffer
end subroutine qp_restore_state_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_wait_to_flush_basic (wait)
!
! Routine to signal to quick_plot whether to wait for flushing the plot buffer.
! Note: By default, quick_plot is not in a wait state.
!
! Input:
!   wait    -- logical: If True, go into a wait state for flushing.
!                       If False, flush the buffer and go into a non-wait state.
!-

subroutine qp_wait_to_flush_basic(wait)
implicit none
logical wait
!
pg_allow_flush = .not. wait

if (pg_allow_flush) then
  call pgebuf        ! flush buffer
else
  call pgbbuf        ! buffer commands
endif

end subroutine qp_wait_to_flush_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_color_basic (ix_color)  
!
! Subroutine to set the color taking into accout that GIF
! inverts the black for white.
!
! Input:
!   ix_color -- Integer: Color index (0 - 15).
!-

subroutine qp_set_color_basic (ix_color, set_background)

implicit none
real(rp) :: real_color
integer ix_color
!integer, parameter :: inverse_color(0:15) = &
!        [1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 ]
integer, parameter :: clist1(0:15) = &
        [5, 10, 3, 9, 7, 8, 6, 13, 2, 11,  4, 12, 12, 15, 14, 1]
logical, optional :: set_background
logical :: set_bg

!
set_bg = .false.
if (present(set_background)) set_bg = set_background


! Error check

if (ix_color < 0) then
  print *, 'ERROR IN QP_SET_PGPLOT: IX_COLOR ARGUMENT OUT OF RANGE:', &
                                                                    ix_color
  if (global_com%exit_on_error) call err_exit
endif  

if (ix_color > 15) then 
  real_color = (ix_color - 17)/ (1.0_rp*(huge(ix_color) - 17) )
  ix_color=  floor( 12*real_color)
  ix_color = clist1(ix_color)
endif

if (set_bg) then
  call pgstbg(ix_color)
else
  call pgsci (ix_color)
endif

end subroutine qp_set_color_basic
  
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_clear_page_basic
!
! Subroutine to clear all drawing from the page.
!-

subroutine qp_clear_page_basic
implicit none
call pgpage
end subroutine qp_clear_page_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_paint_rectangle_basic (x1, x2, y1, y2, color, fill_pattern)
!
! Subroutine to fill a rectangle with a given color. 
! A color of white essentially eraces the rectangle.
! Units are inches from lower left of page.
!
! Input:
!   x1, y1       -- Real(rp): Bottom left corner of box.
!   x2, y2       -- Real(rp): Upper right corner of box.
!   color        -- Integer: Color of rectangle.
!   fill_pattern -- Integer: fill pattern.
!   page_type    -- Character(*): Type of page ('GIF', 'X', etc).
!-

subroutine qp_paint_rectangle_basic (x1, x2, y1, y2, color, fill_pattern)

implicit none

real(rp) x1, x2, y1, y2
integer ci, fs
real  xv1, xv2, yv1, yv2, xw1, xw2, yw1, yw2, f
integer color, fill_pattern

!

if (color == transparent$) return
if (fill_pattern == no_fill$) return

call qp_save_state_basic              ! Buffer the following calls

f = pg_com%page_scale

call qp_set_color_basic(color)        ! Set color index to background
call pgsfs(fill_pattern)              ! Set fill-area pattern to solid

call pgqwin (xw1, xw2, yw1, yw2)      ! get graph data min/max
call pgqvp (0, xv1, xv2, yv1, yv2)    ! get viewport coords

! set the viewport to the box

call pgvsiz (real(f*x1), real(f*x2), real(f*y1), real(f*y2))  

call pgrect (xw1, xw2, yw1, yw2)      ! color the box
call pgsvp (xv1, xv2, yv1, yv2)       ! reset the viewport coords

call qp_restore_state_basic()         ! Flush the buffer.

end subroutine qp_paint_rectangle_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_polyline_basic (x, y) 
!
! Subroutine to draw a polyline.
!
! Input:
!   x(:), y(:) -- Real(rp): (x, y) points for the polyline.
!                     in inches from lower left of page.
!-

subroutine qp_draw_polyline_basic (x, y)

implicit none

real(rp) :: x(:), y(:), f

!

f = pg_com%page_scale

if (size(x) /= size(y)) then
  print *, 'ERROR IN QP_DRAW_POLYLINE_BASIC: X, Y COORD VECTORS HAVE'
  print *, '      UNEQUAL LENGTH!', size(x), size(y)
  if (global_com%exit_on_error) call err_exit
endif

if (size(x) < 2) return
call pgline (size(x), real(f*x), real(f*y))

end subroutine qp_draw_polyline_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+                    
! Subroutine qp_open_page_basic (page_type, x_len, y_len, plot_file, &
!                                          x_page, y_page, i_chan, page_scale)
!
! Subroutine to Initialize a page (window) for plotting.
!
! Note: pgplot does not do a good job with gif files. Consider making
! a postscript file and using the unix command pstogif to convert it.
!
! Input:
!   page_type  -- Character(*). Device name for the type of plot.
!                  TYPE is passed to GG_SETUP. E.g.
!                  TYPE = 'X'     --> Open an X-window.
!                  TYPE = 'GIF'   --> To create a gif file.
!                  TYPE = 'PS'    --> To create a Color PostScript file.
!   x_len      -- Real(rp), optional: Horizontal width in inches.
!   y_len      -- Real(rp), optional: Vertical width in inches.
!   plot_file  -- Character(*), optional: Name for the plot file.
!                    Default is: 'quick_plot.ps' or 'quick_plot.gif'
!   page_scale -- Real(rp), optional: Scale to expand or shrink the drawing.
!                    Default is 1.0.  
!
! Output:
!   x_page    -- Real(rp): Horizontal page size in inches.
!   y_page    -- Real(rp): Vertical page size in inches.
!   i_chan    -- Integer, optional: Plot channel. 
!                 Like a unit number for a fortran OPEN.
!                 To be used with qp_select_page.
!-

subroutine qp_open_page_basic (page_type, x_len, y_len, plot_file, x_page, y_page, i_chan, page_scale)

implicit none

real(rp) x_len, y_len, x_page, y_page
real(rp), optional :: page_scale

real x, y
real x1i, x2i, y1i, y2i, h, xi, yi

integer, optional :: i_chan
integer pgopen, iw, ix

character(*) page_type, plot_file
character(16) :: r_name = 'qp_open_page_basic'

! Set plot type
! GIF does not work well so create a ps file and convert to gif on closing.

select case (page_type)
case ('X')
#if defined (CESR_WINCVF)
  iw = pgopen ('/WV')
#else
  iw = pgopen ('/XWINDOW')
#endif
  if (iw <= 0) stop
  call pgscr (0, 1.0, 1.0, 1.0)    ! white background
  call pgscr (1, 0.0, 0.0, 0.0)    ! black foreground

case ('PS', 'PS-L')
  iw = pgopen (trim(plot_file) // '/CPS')

case ('GIF', 'GIF-L')
  iw = pgopen (trim(plot_file) // '/GIF')
  call pgscr (1, 0.0, 0.0, 0.0)    ! black foreground
  call pgscr (0, 1.0, 1.0, 1.0)    ! white background

case default
  call out_io (s_abort$, r_name, 'ERROR: UNKNOWN PAGE_TYPE: ' // page_type)
  if (global_com%exit_on_error) call err_exit
end select

if (present(i_chan)) i_chan = iw

if (iw <= 0) then
  print *, 'ERROR IN QP_OPEN_PAGE: CANNONT OPEN OUTPUT DEVICE!'
  stop
endif

! set page size. 
! pgplot assumes 72 pixels / inch for X-windows so there needs to be a correction

if (x_len > 0 .and. y_len > 0) then
  call pgpap (real(x_len), real(y_len/x_len))
endif

! do not pause when clearing the screen.

call pgask (.false.)  

! Get page size info.

call pgqvsz (1, x1i, x2i, y1i, y2i)  ! page in inches
x_page = x2i
y_page = y2i

! clear page and set graph min/max

call pgpage
call pgsvp (0.0, 1.0, 0.0, 1.0)  ! viewport to entire page
call pgswin (0.0, x2i, 0.0, y2i) ! set min/max

! get the conversion factor for character height.

call pgqch (h)               ! text height in pgplot units.
call pgqcs (1, xi, yi)       ! size in inches

i_save = i_save + 1
pg_com => pg_interface_save_com(i_save)

pg_com%i_chan = iw
pg_com%qp_to_pg_text_height_factor = h / (xi * 72)
pg_com%page_scale = real_option(1.0_rp, page_scale)
pg_com%page_type = page_type
pg_com%plot_file = plot_file

end subroutine qp_open_page_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_select_page_basic (iw)
!
! Subroutine to switch to a particular page for drawing graphics.
!
! Input:
!   iw -- Integer: ID of page obtained from qp_open_page
!-

subroutine qp_select_page_basic (iw)

implicit none
integer i, iw
!
call out_io (s_abort$, 'qp_select_page_basic', 'NOT YET IMPLEMENTED!')
if (global_com%exit_on_error) call err_exit
call pgslct(iw)

end subroutine qp_select_page_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_close_page_basic ()
!
! Subroutine to finish plotting on a page.
! For X this closes the window.
!-

subroutine qp_close_page_basic ()

implicit none
!
call pgclos
i_save = i_save - 1
pg_com => pg_interface_save_com(i_save)
if (i_save /= 0) then
  call pgslct(pg_com%i_chan)
endif

end subroutine qp_close_page_basic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_end_basic ()
!
! Cleanup routine at the end of plotting.
!-

subroutine qp_end_basic ()

call pgend()

end subroutine qp_end_basic

#endif

end module
