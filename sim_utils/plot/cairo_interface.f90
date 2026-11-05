!+
! Module cairo_interface
!
! Module to interface between the Cairo plotting backend and quick_plot routines.
! This module isolates quick_plot from Cairo.
!
! The drawing is done in C (qp_cairo.c and qp_cairo_render.c). Hardcopy (PS, PDF, SVG, PNG) is
! rendered directly by Cairo. On-screen plots ("X" page type) are sent as a display list to a
! separate viewer program (qp_cairo_viewer) which draws the plot in its own window using the same
! renderer. Running the window in its own process keeps it responsive while the calling program
! (EG Tao) is waiting for input, and works the same on Linux, macOS, and Windows.
!
! Units passed to the C side are inches from the lower left corner of the page,
! heights are in points, and line widths are in quick_plot (PGPLOT) units.
!-

module cairo_interface

use output_mod
use quick_plot_struct

! This #if def wraps the entire module

#if defined (CESR_CAIRO)

use, intrinsic :: iso_c_binding

type cairo_interface_struct
  character(16) :: page_type = ''
  integer :: i_chan = -1
  real(rp) :: page_scale = 1
end type

type (cairo_interface_struct), pointer, save, private :: ca_com => null()
type (cairo_interface_struct), save, target, private :: ca_interface_save_com(0:20)
integer, save, private :: i_page = 0
logical, save, private :: ca_allow_flush = .true.

real(rp), parameter, private :: png_dpi = 100   ! Resolution of PNG ("GIF") output.

! Page kinds known to qp_cairo.c

integer(c_int), parameter, private :: qpc_screen$ = 0, qpc_ps$ = 1, qpc_png$ = 2, qpc_svg$ = 3, qpc_pdf$ = 4

private ca_flush, ca_color_rgb, ca_scale

interface
  function qpc_open (kind, file, nfile, w_in, h_in, dpi) result (ix) bind(c, name = 'qpc_open')
    import
    integer(c_int), value :: kind, nfile
    character(kind=c_char) :: file(*)
    real(c_double), value :: w_in, h_in, dpi
    integer(c_int) ix
  end function

  subroutine qpc_close () bind(c, name = 'qpc_close')
  end subroutine

  subroutine qpc_end () bind(c, name = 'qpc_end')
  end subroutine

  subroutine qpc_flush () bind(c, name = 'qpc_flush')
  end subroutine

  subroutine qpc_clear () bind(c, name = 'qpc_clear')
  end subroutine

  subroutine qpc_push () bind(c, name = 'qpc_push')
  end subroutine

  subroutine qpc_pop () bind(c, name = 'qpc_pop')
  end subroutine

  subroutine qpc_color (r, g, b) bind(c, name = 'qpc_color')
    import
    real(c_double), value :: r, g, b
  end subroutine

  subroutine qpc_line_width (w) bind(c, name = 'qpc_line_width')
    import
    real(c_double), value :: w
  end subroutine

  subroutine qpc_line_pattern (p) bind(c, name = 'qpc_line_pattern')
    import
    integer(c_int), value :: p
  end subroutine

  subroutine qpc_clip (on) bind(c, name = 'qpc_clip')
    import
    integer(c_int), value :: on
  end subroutine

  subroutine qpc_char_size (h) bind(c, name = 'qpc_char_size')
    import
    real(c_double), value :: h
  end subroutine

  subroutine qpc_sym_size (h) bind(c, name = 'qpc_sym_size')
    import
    real(c_double), value :: h
  end subroutine

  subroutine qpc_sym_fill (f) bind(c, name = 'qpc_sym_fill')
    import
    integer(c_int), value :: f
  end subroutine

  subroutine qpc_graph_pos (x1, x2, y1, y2) bind(c, name = 'qpc_graph_pos')
    import
    real(c_double), value :: x1, x2, y1, y2
  end subroutine

  subroutine qpc_text_bg (opaque, r, g, b) bind(c, name = 'qpc_text_bg')
    import
    integer(c_int), value :: opaque
    real(c_double), value :: r, g, b
  end subroutine

  subroutine qpc_polyline (n, x, y) bind(c, name = 'qpc_polyline')
    import
    integer(c_int), value :: n
    real(c_double) :: x(*), y(*)
  end subroutine

  subroutine qpc_rect (x1, x2, y1, y2, fill) bind(c, name = 'qpc_rect')
    import
    real(c_double), value :: x1, x2, y1, y2
    integer(c_int), value :: fill
  end subroutine

  subroutine qpc_symbol (x, y, type) bind(c, name = 'qpc_symbol')
    import
    real(c_double), value :: x, y
    integer(c_int), value :: type
  end subroutine

  subroutine qpc_text (x, y, angle, justify, text, n) bind(c, name = 'qpc_text')
    import
    real(c_double), value :: x, y, angle, justify
    character(kind=c_char) :: text(*)
    integer(c_int), value :: n
  end subroutine

  subroutine qpc_arrow (x1, y1, x2, y2, head_angle, head_barb, head_size, filled) bind(c, name = 'qpc_arrow')
    import
    real(c_double), value :: x1, y1, x2, y2, head_angle, head_barb, head_size
    integer(c_int), value :: filled
  end subroutine

  function qpc_text_len (text, n) result (t_len) bind(c, name = 'qpc_text_len')
    import
    character(kind=c_char) :: text(*)
    integer(c_int), value :: n
    real(c_double) t_len
  end function
end interface

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Page scale of the current page (1 if no page is open).

function ca_scale() result (f)
real(rp) f
f = 1
if (associated(ca_com)) f = ca_com%page_scale
end function

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Send pending screen updates to the viewer unless quick_plot has asked us to wait.

subroutine ca_flush()
if (ca_allow_flush) call qpc_flush()
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Convert a quick_plot color index to RGB in [0, 1].
! Continuous colors use the same HLS map as the PLPLOT and PGPLOT backends.

subroutine ca_color_rgb (ix_color, rgb)

implicit none

integer ix_color
real(rp) rgb(3), real_color, hue, chroma, x, m
integer h_sector

!

if (ix_color == qp_custom_color$) then
  rgb = qp_custom_rgb / 255.0_rp

elseif (ix_color >= qp_continuous_color_start$) then
  real_color = (ix_color - qp_continuous_color_start$) / (1.0_rp*(huge(ix_color) - qp_continuous_color_start$))
  hue = 240.0_rp * (1.0_rp - real_color)
  chroma = 0.64_rp
  m = 0.28_rp
  x = chroma * (1.0_rp - abs(mod(hue / 60.0_rp, 2.0_rp) - 1.0_rp))
  h_sector = int(hue / 60.0_rp)
  select case (h_sector)
  case (0);     rgb = [chroma, x, 0.0_rp]
  case (1);     rgb = [x, chroma, 0.0_rp]
  case (2);     rgb = [0.0_rp, chroma, x]
  case (3);     rgb = [0.0_rp, x, chroma]
  case default; rgb = [0.0_rp, 0.0_rp, chroma]
  end select
  rgb = rgb + m

elseif (ix_color >= 0 .and. ix_color <= qp_max_color_index$) then
  rgb = [qp_color_red(ix_color), qp_color_green(ix_color), qp_color_blue(ix_color)] / 255.0_rp

else
  rgb = 0
endif

end subroutine ca_color_rgb

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
!   x1, y1 -- Real(rp): Bottom left corner of graph
!   x2, y2 -- Real(rp): Upper right corner of graph
!-

subroutine qp_set_graph_position_basic (x1, x2, y1, y2)
implicit none
real(rp) x1, x2, y1, y2, f
if (x1 == x2 .or. y1 == y2) return
f = ca_scale()
call qpc_graph_pos (f*x1, f*x2, f*y1, f*y2)
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_symbol_size_basic (height, symbol_type, uniform_size)
!
! Subroutine to set the symbol_size.
! The Cairo symbols are drawn as geometry and so are already uniform in size.
!
! Input:
!   height       -- Real(rp): Symbol height in points.
!   symbol_type  -- Integer: Symbol type.
!   uniform_size -- Logical: Make all symbols the save size for constant height.
!-

subroutine qp_set_symbol_size_basic (height, symbol_type, uniform_size)
implicit none
real(rp) height
integer symbol_type
logical uniform_size
call qpc_sym_size (ca_scale() * height)
end subroutine

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
!   fill_pattern -- Integer: Fill pattern.
!-

subroutine qp_paint_rectangle_basic (x1, x2, y1, y2, color, fill_pattern)
implicit none
real(rp) x1, x2, y1, y2, f
integer color, fill_pattern

if (x1 == x2 .or. y1 == y2) return
if (color == transparent$) return

f = ca_scale()
call qpc_push()
call qp_set_color_basic (color)
call qpc_rect (f*x1, f*x2, f*y1, f*y2, int(fill_pattern, c_int))
call qpc_pop()
call ca_flush()
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_symbol_fill_basic (fill)
!
! Subroutine to set the symbol fill pattern.
!
! Input:
!   fill -- Integer: Fill pattern.
!-

subroutine qp_set_symbol_fill_basic (fill)
implicit none
integer fill
call qpc_sym_fill (int(fill, c_int))
end subroutine

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
call qpc_line_width (real(line_width, c_double))
end subroutine

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
call qpc_line_pattern (int(line_pattern, c_int))
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_clip_basic (clip)
!
! Subroutine to set the clipping state.
! Note: This affects lines, symbols, text, and arrows but not filled rectangles.
!
! Input:
!   clip    -- Logical: Clip at graph boundary?
!-

subroutine qp_set_clip_basic (clip)
implicit none
logical :: clip
if (clip) then
  call qpc_clip(1_c_int)
else
  call qpc_clip(0_c_int)
endif
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_char_size_basic (height)
!
! Subroutine to set the character size.
!
! Input:
!   height -- Real(rp): Character size in points.
!-

subroutine qp_set_char_size_basic (height)
implicit none
real(rp) height
call qpc_char_size (ca_scale() * height)
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_text_background_color_basic (color)
!
! Subroutine to set the character text background color.
!
! Input:
!   color -- Integer: Color index. Negative => transparent.
!-

subroutine qp_set_text_background_color_basic (color)
implicit none
integer color
real(rp) rgb(3)

if (color < 0 .and. color /= qp_custom_color$) then
  call qpc_text_bg (0_c_int, 1.0_rp, 1.0_rp, 1.0_rp)
else
  call ca_color_rgb (color, rgb)
  call qpc_text_bg (1_c_int, rgb(1), rgb(2), rgb(3))
endif
end subroutine

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
character(*) text
real(rp) t_len
t_len = qpc_text_len (text, int(len_trim(text), c_int)) / ca_scale()
end function

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_draw_text_basic (text, len_text, x0, y0, angle, justify)
!
! Subroutine to draw text.
! Text may contain PGPLOT escape sequences (\u, \d, \gX, \fX, etc.).
!
! Input:
!   text     -- Character(*): Text to draw.
!   len_text -- Integer: Length of text.
!   x0, y0   -- Real(rp): Position of text baseline in inches from bottom left of page.
!   angle    -- Real(rp): Rotation angle of text in degrees.
!   justify  -- Real(rp): Left/Right justification. 0 = left, 0.5 = center, 1 = right.
!-

subroutine qp_draw_text_basic (text, len_text, x0, y0, angle, justify)
implicit none
character(*) text
integer len_text, n
real(rp) x0, y0, angle, justify, f

f = ca_scale()
n = min(len_text, len(text))
call qpc_text (f*x0, f*y0, angle, justify, text, int(n, c_int))
call ca_flush()
end subroutine

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
type (qp_arrow_struct) arrow
integer(c_int) filled

f = ca_scale()
filled = 1
if (index(downcase(arrow%head_type), 'outline') /= 0) filled = 0

call qpc_push()
call qp_set_color_basic (qp_string_to_enum(arrow%color, 'color'))
call qpc_arrow (f*r1(1), f*r1(2), f*r2(1), f*r2(2), arrow%head_angle, arrow%head_barb, &
                f*arrow%head_size, filled)
call qpc_pop()
call ca_flush()
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
f = ca_scale()
call qpc_symbol (f*x, f*y, int(symbol, c_int))
call ca_flush()
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_save_state_basic
!
! Subroutine to save the print state.
!-

subroutine qp_save_state_basic
implicit none
call qpc_push()
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_restore_state_basic
!
! Subroutine to restore the print state.
!-

subroutine qp_restore_state_basic
implicit none
call qpc_pop()
call ca_flush()
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_wait_to_flush_basic(wait)
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
ca_allow_flush = .not. wait
call ca_flush()
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_color_basic (ix_color)
!
! Subroutine to set the color.
!
! Input:
!   ix_color -- Integer: Color index.
!-

subroutine qp_set_color_basic (ix_color)
implicit none
integer ix_color
real(rp) rgb(3)

if (ix_color < 0 .and. ix_color /= qp_custom_color$) then
  print *, 'ERROR IN QP_SET_COLOR_BASIC: IX_COLOR ARGUMENT OUT OF RANGE:', ix_color
  if (global_com%exit_on_error) call err_exit
endif

call ca_color_rgb (ix_color, rgb)
call qpc_color (rgb(1), rgb(2), rgb(3))
end subroutine

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
call qpc_clear()
call ca_flush()
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_clear_box_basic (x1, x2, y1, y2)
!
! Subroutine to clear all drawing from a box.
! That is, white out the box region.
! Units are inches from lower left of page.
!
! Input:
!   x1, y1 -- Real(rp): Bottom left corner of box in inches.
!   x2, y2 -- Real(rp): Upper right corner of box in inches.
!-

subroutine qp_clear_box_basic (x1, x2, y1, y2)
implicit none
real(rp) x1, x2, y1, y2, f
f = ca_scale()
call qpc_push()
call qpc_color (1.0_rp, 1.0_rp, 1.0_rp)
call qpc_rect (f*x1, f*x2, f*y1, f*y2, int(solid_fill$, c_int))
call qpc_pop()
call ca_flush()
end subroutine qp_clear_box_basic

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
real(rp) :: x(:), y(:)
real(c_double) :: xs(size(x)), ys(size(y))

if (size(x) /= size(y)) then
  print *, 'ERROR IN QP_DRAW_POLYLINE_BASIC: X, Y COORD VECTORS HAVE'
  print *, '      UNEQUAL LENGTH!', size(x), size(y)
  if (global_com%exit_on_error) call err_exit
  return
endif

if (size(x) < 2) return
xs = ca_scale() * x
ys = ca_scale() * y
call qpc_polyline (int(size(x), c_int), xs, ys)
call ca_flush()
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_open_page_basic (page_type, x_len, y_len, plot_file, x_page, y_page, i_chan, page_scale)
!
! Subroutine to Initialize a page (window) for plotting.
!
! Input:
!   page_type  -- Character(*). Device name for the type of plot:
!                   'X'               -> On screen window (via the qp_cairo_viewer program).
!                   'PS', 'PS-L'      -> Encapsulated PostScript.
!                   'GIF', 'GIF-L'    -> PNG image (PLPLOT compatible: GIF is not supported by Cairo).
!                   'PDF', 'PDF-L'    -> PDF.
!                   'SVG'             -> SVG.
!   x_len      -- Real(rp): Horizontal width in inches (already scaled by page_scale).
!   y_len      -- Real(rp): Vertical width in inches (already scaled by page_scale).
!   plot_file  -- Character(*): Name for the plot file.
!   page_scale -- Real(rp), optional: Scale to expand or shrink the drawing. Default is 1.0.
!
! Output:
!   x_page    -- Real(rp): Horizontal page size in inches.
!   y_page    -- Real(rp): Vertical page size in inches.
!   i_chan    -- Integer, optional: Plot channel.
!-

subroutine qp_open_page_basic (page_type, x_len, y_len, plot_file, x_page, y_page, i_chan, page_scale)

implicit none

real(rp) x_len, y_len, x_page, y_page
real(rp), optional :: page_scale
integer, optional :: i_chan
integer(c_int) kind, ix
character(*) page_type, plot_file
character(*), parameter :: r_name = 'qp_open_page_basic'

!

select case (page_type)
case ('X');            kind = qpc_screen$
case ('PS', 'PS-L');   kind = qpc_ps$
case ('GIF', 'GIF-L'); kind = qpc_png$
case ('SVG');          kind = qpc_svg$
case ('PDF', 'PDF-L'); kind = qpc_pdf$
case default
  call out_io (s_abort$, r_name, 'ERROR: UNKNOWN PAGE_TYPE: ' // page_type)
  if (global_com%exit_on_error) call err_exit
  return
end select

ix = qpc_open (kind, plot_file, int(len_trim(plot_file), c_int), x_len, y_len, png_dpi)
if (ix < 0) then
  call out_io (s_abort$, r_name, 'ERROR: CANNOT OPEN PAGE OF TYPE: ' // trim(page_type) // &
                                 '  FILE: ' // plot_file)
  if (global_com%exit_on_error) call err_exit
  return
endif

i_page = i_page + 1
ca_com => ca_interface_save_com(i_page)
ca_com%page_type = page_type
ca_com%i_chan = ix
ca_com%page_scale = real_option(1.0_rp, page_scale)

if (present(i_chan)) i_chan = ix
x_page = x_len
y_page = y_len

call qp_set_clip_basic(.false.)

end subroutine

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
integer iw
call out_io (s_abort$, 'qp_select_page_basic', 'NOT YET IMPLEMENTED!')
if (global_com%exit_on_error) call err_exit
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_close_page_basic
!
! Subroutine to finish plotting on a page.
! For X this closes the window.
!-

subroutine qp_close_page_basic
implicit none
call qpc_close()
i_page = max(0, i_page - 1)
ca_com => ca_interface_save_com(i_page)
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_end_basic ()
!
! Cleanup routine at the end of plotting.
!-

subroutine qp_end_basic ()
call qpc_end()
i_page = 0
ca_com => ca_interface_save_com(0)
end subroutine qp_end_basic

#endif

end module
