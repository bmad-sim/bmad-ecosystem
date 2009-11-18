#include "CESR_platform.inc"
                                    
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

type pg_interface_struct
  character(16) page_type
  integer :: i_chan = -1
  real qp_to_pg_text_height_factor
  real page_scale   ! scaling for entire page
end type

type (pg_interface_struct), target, save, private :: pg_com
type (pg_interface_struct), save, private :: pg_interface_save_com(10)
integer, save, private :: i_save = 0

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
!+

subroutine qp_set_graph_position_basic (x1, x2, y1, y2)
  implicit none
  real(rp) x1, x2, y1, y2, f
  f = pg_com%page_scale
  call pgvsiz (real(f*x1), real(f*x2), real(f*y1), real(f*y2))
  call pgswin (real(f*x1), real(f*x2), real(f*y1), real(f*y2))
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
    case (1)            ! dot$
      h = h * 2.0       ! I like bigger dots
    case (16)           ! square_filled$
      h = h * 1.56
    case (17)           ! circle_filled$
      h = h * 1.60
    case (18)           ! star5_filled$
      h = h * 1.30
    case (10)           ! square_concave$
      h = h * 0.73
    end select

    if (pg_com%page_type == 'X' .or. pg_com%page_type == 'TK') then
      select case (symbol_type)
      case (8)           ! circle_plus$
        h = h * 0.55
      case (9)           ! circle_dot$
        h = h * 0.59
      case (4)           ! circle$
        h = h * 0.89
      case (13)          ! triangle_filled$
        h = h * 1.22
      end select
    else
      if (symbol_type == 13) h = h * 1.03   ! triangle_filled$
    endif

  endif

  call pgsch (f * h)   ! set symbol size

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_symbol_fill_basic (fill)
!
! Subroutine to set the symbol fill style.
!
! Input:
!   fill -- Integer: Fill style.
!+

subroutine qp_set_symbol_fill_basic (fill)
  implicit none
  integer fill
  call pgsfs (fill)       ! set fill
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
  call pgslw (line_width) ! set line width
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_line_style_basic (style)
!
! Subroutine to set the line style.
!
! Input:
!   style -- Integer: Line style.
!+

subroutine qp_set_line_style_basic (style)
  implicit none
  integer style
  call pgsls  (style)       ! Set style
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
  if (clip) then
    call pgsclp (1)     ! Clip on
  else
    call pgsclp (0)     ! Clip off
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
!   height -- Integer: Character size.
!+

subroutine qp_set_char_size_basic (height)
  implicit none
  real(rp) height, f
  f = pg_com%page_scale
  call pgsch(real(f * height * pg_com%qp_to_pg_text_height_factor))
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
  call pgstbg (color)            ! set text background color
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function qp_text_len_basic (text, len_text)
!
! Function to find the length of a text string.
!
! input:
!   text -- Character(*): Text string.
!
! Output:
!   qp_text_len -- Real(rp): Length of text in inches.
!-

function qp_text_len_basic (text, len_text) result (t_len)

  implicit none

  real(rp) t_len
  real tl, dum, f
  integer len_text
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
!+

subroutine qp_draw_text_basic (text, len_text, x0, y0, angle, justify)
  implicit none
  character(*) text
  integer len_text
  real(rp) x0, y0, angle, justify, f
  f = pg_com%page_scale
  call pgptxt (real(f*x0), real(f*y0), real(angle), real(justify), trim(text))
end subroutine

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
  real(rp) x, y, f
  integer symbol
  f = pg_com%page_scale
  call pgpt1 (real(f*x), real(f*y), symbol)
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
  call pgbbuf     ! buffer commands
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
  call pgebuf        ! flush buffer
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
!   ix_color -- Integer: Color index (0 - 15).
!-

subroutine qp_set_color_basic (ix_color)

  implicit none

  integer ix_color
  integer, parameter :: inverse_color(0:15) = &
          (/ 1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 /)
!            1, 0, 5, 6, 7, 2, 3, 4, 11, 12, 13,  8,  9, 10, 15, 14 
!            0  1  2  3  4  5  6  7   8   9  10  11  12  13  14  15 

! Error check

  if (ix_color < 0 .or. ix_color > 15) then
    print *, 'ERROR IN QP_SET_PGPLOT: IX_COLOR ARGUMENT OUT OF RANGE:', &
                                                                      ix_color
    call err_exit
  endif

! Set pgplot color

  if (pg_com%page_type == 'GIF') then
    call pgsci (inverse_color(ix_color))
  else
    call pgsci (ix_color)
  endif

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
  call pgpage
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
!   fill_pattern -- Integer: Fill style.
!   page_type    -- Character(*): Type of page ('GIF', 'X', etc).
!-

subroutine qp_paint_rectangle_basic (x1, x2, y1, y2, color, fill_pattern)

  implicit none

  real(rp) x1, x2, y1, y2
  integer ci, fs
  real  xv1, xv2, yv1, yv2, xw1, xw2, yw1, yw2, f
  integer color, fill_pattern

!

  call qp_save_state_basic              ! Buffer the following calls

  f = pg_com%page_scale

  call qp_set_color_basic(color)        ! Set color index to background
  call pgsfs(fill_pattern)              ! Set fill-area style to solid

  call pgqwin (xw1, xw2, yw1, yw2)      ! get graph data min/max
  call pgqvp (0, xv1, xv2, yv1, yv2)    ! get viewport coords

! set the viewport to the box

  call pgvsiz (real(f*x1), real(f*x2), real(f*y1), real(f*y2))  

  call pgrect (xw1, xw2, yw1, yw2)      ! color the box
  call pgsvp (xv1, xv2, yv1, yv2)       ! reset the viewport coords
  
  call qp_restore_state_basic           ! Flush the buffer.

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

  real(rp) :: x(:), y(:), f

!

  f = pg_com%page_scale

  if (size(x) /= size(y)) then
    print *, 'ERROR IN QP_DRAW_POLYLINE_BASIC: X, Y COORD VECTORS HAVE'
    print *, '      UNEQUAL LENGTH!', size(x), size(y)
    call err_exit
  endif

  if (size(x) < 2) return
  call pgline (size(f*x), real(f*x), real(f*y))

end subroutine

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
!                  TYPE = 'GIF-L' --> Gif w/ landscape page orientation.
!                  TYPE = 'PS'    --> To create a Color PostScript file.
!                  TYPE = 'PS-L'  --> PostScript w/ landscape page orientation.
!   x_len      -- Real(rp), optional: Horizontal width in inches. 0 => Ignore.
!   y_len      -- Real(rp), optional: Vertical width in inches. 0 => Ignore.
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

subroutine qp_open_page_basic (page_type, x_len, y_len, plot_file, &
                                              x_page, y_page, i_chan, page_scale)

  implicit none

  real(rp) x_len, y_len, x_page, y_page
  real(rp), optional :: page_scale

  real x, y
  real x1i, x2i, y1i, y2i, h, xi, yi

  integer, optional :: i_chan
  integer pgopen, iw, ix

  character(*) page_type, plot_file
  character(16) :: r_name = 'qp_open_page_basic'

! set plot type

  if (page_type == 'X' .or. page_type == 'TK') then
#if defined (CESR_WINCVF)
    iw = pgopen ('/WV')
#else
    iw = pgopen ('/XWINDOW')
#endif
    if (iw <= 0) stop
    call pgscr (0, 1.0, 1.0, 1.0)    ! white background
    call pgscr (1, 0.0, 0.0, 0.0)    ! black foreground

  elseif (page_type == 'PS') then
    iw = pgopen (trim(plot_file) // '/VCPS')

  elseif (page_type == 'PS-L') then
    iw = pgopen (trim(plot_file) // '/CPS')

  elseif (page_type == 'GIF') then
    iw = pgopen (trim(plot_file) // '/GIF')
    call pgscr (0, 1.0, 1.0, 1.0)    ! white background
    call pgscr (1, 0.0, 0.0, 0.0)    ! black foreground

  elseif (page_type == 'GIF-L') then
    iw = pgopen (trim(plot_file) // '/VGIF')
    call pgscr (0, 1.0, 1.0, 1.0)    ! white background
    call pgscr (1, 0.0, 0.0, 0.0)    ! black foreground

  else
    call out_io (s_abort$, r_name, 'ERROR: UNKNOWN PAGE_TYPE: ' // page_type)
    call err_exit
  endif

  if (present(i_chan)) i_chan = iw

  if (iw <= 0) then
    print *, 'ERROR IN QP_OPEN_PAGE: CANNONT OPEN OUTPUT DEVICE!'
    stop
  endif

! set page size

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
  pg_interface_save_com(i_save)%i_chan = iw
  pg_interface_save_com(i_save)%qp_to_pg_text_height_factor = h / (xi * 72)
  pg_interface_save_com(i_save)%page_scale = real_option(1.0_rp, page_scale)
  pg_interface_save_com(i_save)%page_type = page_type

  pg_com = pg_interface_save_com(i_save)

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
  integer i, iw
  call pgslct(iw)
  do i = 1, size(pg_interface_save_com)
    if (pg_interface_save_com(i)%i_chan == iw) then
      pg_interface_save_com(i_save) = pg_interface_save_com(i)
      pg_interface_save_com(i) = pg_com
      pg_com = pg_interface_save_com(i_save)
      return
    endif
  enddo
  call out_io (s_abort$, 'qp_select_page_basic', 'BAD PAGE ID: \i\ ', iw)
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
  call pgclos
  i_save = i_save - 1
  if (i_save /= 0) pg_com = pg_interface_save_com(i_save)
  call pgslct(pg_com%i_chan)
end subroutine

end module
