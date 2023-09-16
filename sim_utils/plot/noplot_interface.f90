!+
! Module noplot_interface
!
! Module to define dummy routines for quick_plot.
!
! Noplot_interface is used when quick_plot is not going to be used and it is desired that neither 
! plplot nor pgplot is to be compiled. For example, this is useful when the GUI version of Tao is used 
! on Windows where compiling plplot and pgplot can be problematical.
!-

module noplot_interface

use output_mod
use quick_plot_struct


! This #if def wraps the entire module

#if defined (CESR_NOPLOT)

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
!   x1, y1 -- Real(rp): Bottom left corner of graph
!   x2, y2 -- Real(rp): Upper right corner of graph
!+

subroutine qp_set_graph_position_basic (x1, x2, y1, y2)

implicit none
real(rp) x1, x2, y1, y2, x1m, x2m, y1m, y2m, xp1, xp2, yp1, yp2, fx, fy

!

end subroutine

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
!+

subroutine qp_set_symbol_size_basic (height, symbol_type, uniform_size)

implicit none

real(rp) height, h, d, dum

integer symbol_type

logical uniform_size

!

end subroutine

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

real(rp) x1, x2, y1, y2
real(rp) xv1, xv2, yv1, yv2, xw1, xw2, yw1, yw2, fx, fy
integer color, fill_pattern
integer ci, fs

!

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
!+

subroutine qp_set_symbol_fill_basic (fill)
implicit none
integer fill

!

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
!+

subroutine qp_set_line_width_basic (line_width)
implicit none  
integer line_width
 
!

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
!+

subroutine qp_set_line_pattern_basic (line_pattern)
implicit none
integer line_pattern

!

end subroutine

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
real(rp) xp1, xp2, yp1, yp2

!

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
!   height -- Integer: Character size.
!+

subroutine qp_set_char_size_basic (height)
implicit none
real(rp) height, d, h

!

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
!   color -- Integer: Color index.
!+

subroutine qp_set_text_background_color_basic (color)
implicit none
integer color

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

real(rp) t_len, d, h
integer i, n_text
character(*) text

! 

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
!+

subroutine qp_draw_text_basic (text, len_text, x0, y0, angle, justify)
implicit none
character(*) text
character(len(text)+20) text2
integer len_text, i, ix
real(rp) x0, y0, angle, justify, dx, dy, x0m, y0m, d, h, t_len
real(rp), parameter :: pi=3.141592

!

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

!

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
!+

subroutine qp_draw_symbol_basic (x, y, symbol)
implicit none
real(rp) x, y, xm, ym
integer symbol

!

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_save_state_basic 
!
! Subroutine to save the print state.
!+

subroutine qp_save_state_basic 
implicit none

!

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_restore_state_basic ()
!
! Subroutine to restore the print state.
!+

subroutine qp_restore_state_basic ()
implicit none
real(rp) def,dum,xp1,xp2,yp1,yp2

!

end subroutine

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
!   ix_color -- Integer: Color index (0 - 16).
!-

subroutine qp_set_color_basic (ix_color)

implicit none

integer ix_color

!

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

!

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

real(rp) x1, x2, y1, y2, x1m, x2m, y1m, y2m
real(rp) :: x_vec(0:4)
real(rp) :: y_vec(0:4)

!

end subroutine

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
real(rp) :: xm(size(x))
real(rp) :: ym(size(y))

!

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+                    
! Subroutine qp_open_page_basic (page_type, x_len, y_len, plot_file, &
!                                            x_page, y_page, i_chan, page_scale)
!
! Subroutine to Initialize a page (window) for plotting.
!
! Input:
!   page_type -- Character(*). Device name for the type of plot.
!                  See qp_open_page for more details.
!   x_len     -- Real(rp), optional: Horizontal width in inches, Not used with PS.
!   y_len     -- Real(rp), optional: Vertical width in inches, Not used with PS.
!   plot_file -- Character(*), optional: Name for the plot file.
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

real(rp) x_len, y_len, x_page, y_page, x1i, x2i, y1i, y2i, d, h
real(rp), optional :: page_scale

integer, optional :: i_chan
integer ix, xp, yp, ix_len, iy_len, i_ch, stat

character(*) page_type, plot_file
character(40) geom
character(16) :: r_name = 'qp_open_page_basic'

!

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

!

end subroutine

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

end subroutine qp_end_basic

#endif

end module
