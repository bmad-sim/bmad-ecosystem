!+
! Module plplot_interface
!
! Module to interface between plplot and quick_plot routines.
! This module isolates quick_plot from plplot.
!
! To tell the difference between quick_plot and plplot: 
!   QUICK_PLOT subroutines start with "QP_"
!   PLPLOT routines start with "PL" and do not have any "_".
!
! The correspondence between PAGE, BOX, and GRAPH and PLPLOT is:
!   QUICK_PLOT    PLPLOT
!   ----------    ------
!   PAGE          VIEW SURFACE 
!   BOX           No corresponding entity.
!   GRAPH         VIEWPORT and WINDOW
!
! Essentially the VIEWPORT is the region outside of which lines and symbols
! will be clipped (if clipping is turned on) and the WINDOW defines the 
! plot area. I'm not sure why PLPLOT makes a distinction, but VIEWPORT and 
! WINDOW always are the same region.
!
! Note: plwind is called to reset the plplot internal level to 3.
!
! This is a list of all logicals which need to be set:
!   'uniform_size' in 'qp_set_symbol_size_basic'
!   'clip' in 'qp_set_clip_basic'
!-

module plplot_interface

use output_mod
use quick_plot_struct

type viewport_size
  real(rp) :: x1 !in mm
  real(rp) :: x2
  real(rp) :: y1
  real(rp) :: y2
end type

type pl_interface_struct
  character(16) page_type
  type (viewport_size) :: graph_pos
  integer :: i_chan = -1
  integer :: fill_style
  integer :: line_width
  integer :: line_style
  real(rp) :: char_size
  real(rp) :: sym_size
  integer :: fg_color
  logical :: clip
  real(rp) page_scale
  real(rp) x_inch_page
  real(rp) y_inch_page
  real(rp) x_inch_to_mm
  real(rp) y_inch_to_mm
end type

type (pl_interface_struct), pointer, save, private :: pl_com
type (pl_interface_struct), save, target, private :: pl_interface_save_com(0:20)
integer, save, private :: i_save = 0
real(rp), parameter, private :: point_to_mm_conv = .25   !approximate

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

  if (x1 == x2 .or. y1 == y2) return

  fx = pl_com%x_inch_to_mm 
  fy = pl_com%y_inch_to_mm 

  x1m = x1 * fx
  x2m = x2 * fx
  y1m = y1 * fy
  y2m = y2 * fy

  pl_com%graph_pos%x1 = x1m
  pl_com%graph_pos%x2 = x2m
  pl_com%graph_pos%y1 = y1m
  pl_com%graph_pos%y2 = y2m

  if (.not. pl_com%clip) then
    call plgspa(xp1, xp2, yp1, yp2)               ! Get current subpage in mm
    call plvpor(0.0_rp, 1.0_rp, 0.0_rp, 1.0_rp)   ! Set viewport in normalized coords.
    call plwind(0.0_rp, xp2-xp1, 0.0_rp, yp2-yp1) ! Set world coords
  else
    call plsvpa (x1m, x2m, y1m, y2m)  ! Set viewport in abs coords
    call plwind (x1m, x2m, y1m, y2m)  ! Set world coords
  endif

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

! The PLPLOT symbol set does not have a constant symbol size.
! This generally does not look nice so renormalize to get a consistant size.
! This excludes the set of circles with different sizes.
        
  h = height * pl_com%page_scale / 3  ! 3 => conversion from points to mm.

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

    if (pl_com%page_type == 'X' .or. pl_com%page_type == 'TK') then
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

  ! Set symbol size and save this state.

  if (pl_com%page_type == 'X' .or. pl_com%page_type(1:2) == 'PS') then
    h = h * 0.5
  endif

  call plssym(h, 1.0_rp)  ! Set symbol scale factor in mm.
  pl_com%sym_size = h

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
!   fill_pattern -- Integer: Fill style.
!-

subroutine qp_paint_rectangle_basic (x1, x2, y1, y2, color, fill_pattern)

  implicit none

  real(rp) x1, x2, y1, y2
  real(rp) xv1, xv2, yv1, yv2, xw1, xw2, yw1, yw2, fx, fy
  integer color, fill_pattern
  integer ci, fs

!

  if (x1 == x2 .or. y1 == y2) return

  call qp_save_state_basic              ! Buffer the following calls

  fx = pl_com%x_inch_to_mm
  fy = pl_com%y_inch_to_mm

  call qp_set_color_basic(color)    ! Set color index to background

  select case (fill_pattern)
  case (hatched$) 
    call plpsty(3)
  case (cross_hatched$)
    call plpsty(5)
  case (solid_fill$)
    call plpsty(0)
  end select

  if (fill_pattern == no_fill$) then
    call plline (4, fx * (/ x1, x2, x2, x1 /), fy * (/ y1, y1, y2, y2 /)) ! No fill
  else
    call plfill (4, fx * (/ x1, x2, x2, x1 /), fy * (/ y1, y1, y2, y2 /)) ! color the box
  endif

  call qp_restore_state_basic                 ! Flush the buffer.

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

  call plpsty (fill)       ! set fill

  !Save this state
  pl_com%fill_style = fill

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

  call plwid (line_width) ! set line width

  !Save this state
  pl_com%line_width = line_width

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

  call pllsty (style)       ! Set style

  !Save this state
  pl_com%line_style = style

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
  type(viewport_size), pointer :: gp

  !

  if (.not. clip) then
     call plgspa(xp1, xp2, yp1, yp2)               ! Get current subpage in mm
     call plvpor(0.0_rp, 1.0_rp, 0.0_rp, 1.0_rp)   ! Set viewport in normalized coords.
     call plwind(0.0_rp, xp2-xp1, 0.0_rp, yp2-yp1) ! Set world coords
  else
    gp => pl_com%graph_pos
    call plsvpa (gp%x1, gp%x2, gp%y1, gp%y2)  ! Set viewport in abs coords
    call plwind (gp%x1, gp%x2, gp%y1, gp%y2)  ! Set wold coords
  endif

   pl_com%clip = clip

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

  ! call plschr(pl_com%default_cs, height)
  h = height
  if (pl_com%page_type(1:2) == 'X') then
    h = h * 0.7
  elseif (pl_com%page_type == 'TK') then
    h = h * 1.3    
  elseif (pl_com%page_type(1:2) == 'PS') then
    h = h * 0.85
  endif

  call plschr(0.0_rp, h)  ! Set
  call plgchr(d,h)        ! Get

  !Save this state
  pl_com%char_size = h

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
  !call pgstbg (color)            ! set text background color
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function qp_text_len_basic (text, len_text) result (t_len)
!
! Function to find the length of a text string.
!
! input:
!   text -- Character(*): Text string.
!
! Output:
!   t_len -- Real(rp): Length of text in inches.
!-

function qp_text_len_basic (text, len_text) result (t_len)

  implicit none

  real(rp) t_len, d, h
  integer len_text
  character(*) text

  ! This is kind-of a 1st order approx since there is no subroutine for this action

  call plgchr(d,h)

  t_len = 0.8 * 0.03937 * len_trim(text) * h  

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
  character(len(text)) text2
  integer len_text, i, ix
  real(rp) x0, y0, angle, justify, dx, dy, x0m, y0m, d, h, t_len
  real(rp), parameter :: pi=3.141592

  ! plplot uses '#' for the meta character instead of pgplot's '\'

  text2 = text
  do
    ix = index(text2, '\') ! '
    if (ix == 0) exit
    text2(ix:ix) = '#'
  enddo

  !

  call plgchr(d,h)
  t_len = len_trim(text2)*h
  dx = cos(angle*pi/180)
  dy = sin(angle*pi/180)
  x0m = x0 * pl_com%x_inch_to_mm - 0.5*h*dy  ! x0, y0 specify coordinates of text baseline
  y0m = y0 * pl_com%y_inch_to_mm + 0.5*h*dx  ! but plptex needs the coordinates of the midline

  call plptex (x0m, y0m, dx, dy, justify, trim(text2))

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
  real(rp) x, y, xm, ym
  integer symbol
  
  xm = x * pl_com%x_inch_to_mm
  ym = y * pl_com%y_inch_to_mm

  call plpoin (1, xm, ym, symbol)
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

  i_save = i_save + 1
  pl_interface_save_com(i_save) = pl_com
  pl_com => pl_interface_save_com(i_save) 

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
  type(viewport_size), pointer :: gp
  real(rp) def,dum,xp1,xp2,yp1,yp2

  !

  i_save = i_save - 1
  pl_com => pl_interface_save_com(i_save)

  call plgchr(def, dum)
  call plschr(def, pl_com%char_size)

  call plssym(pl_com%sym_size, 1.0_rp)  ! Set symbol scale factor

  if(pl_com%fill_style /= 0) then
     call plpsty(pl_com%fill_style)
  endif

  call plwid(pl_com%line_width)
  
  if(pl_com%line_style /= 0) then
     call pllsty(pl_com%line_style)
  endif

  call plcol0(pl_com%fg_color)

  if (pl_com%clip) then
    gp => pl_com%graph_pos
    if (gp%x1 /= gp%x2 .and. gp%y1 /= gp%y2) then
      call plsvpa (gp%x1, gp%x2, gp%y1, gp%y2)  ! Set viewport in abs coords
      call plwind (gp%x1, gp%x2, gp%y1, gp%y2)  ! Set world coords
    endif
  else
    call plgspa(xp1, xp2, yp1, yp2)               ! Get current subpage in mm
    call plvpor(0.0_rp, 1.0_rp, 0.0_rp, 1.0_rp)   ! Set viewport in normalized coords.
    call plwind(0.0_rp, xp2-xp1, 0.0_rp, yp2-yp1) ! Set world coords
  endif
  
  call plflush()

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
    print *, 'ERROR IN QP_SET_COLOR_BASIC: IX_COLOR ARGUMENT OUT OF RANGE:', &
                                                                      ix_color
    call err_exit
  endif

! Set plplot color

  !if (page_type == 'GIF') then
    !call plcol0 (inverse_color(ix_color))
    !Save this state
     !pl_com%fg_color = inverse_color(ix_color)
  !else
    call plcol0 (ix_color)
    !Save this state
    pl_com%fg_color = ix_color
  !endif

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
  ! plclear should work but does not.
  ! So also call qp_clear_box_basic which does the job.
  call plclear()  
  call qp_clear_box_basic (0.0_rp, pl_com%x_inch_page, 0.0_rp, pl_com%y_inch_page)
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

  x1m = pl_com%x_inch_to_mm * x1
  x2m = pl_com%x_inch_to_mm * x2
  y1m = pl_com%y_inch_to_mm * y1
  y2m = pl_com%y_inch_to_mm * y2
!
  call qp_save_state_basic              ! Buffer the following calls

  call qp_set_color_basic(0)            ! Set color index to background
  call plpsty(0)                        ! Set fill-area style to solid
  
  x_vec = (/x1m, x2m, x2m, x1m, x1m/)
  y_vec = (/y1m, y1m, y2m, y2m, y1m/)
  call plfill (5, x_vec, y_vec)         ! Fills a polygon with 4 vertices
  
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

  real(rp) :: x(:), y(:)
  real(rp) :: xm(size(x))
  real(rp) :: ym(size(y))

!

  if (size(x) /= size(y)) then
    print *, 'ERROR IN QP_DRAW_POLYLINE_BASIC: X, Y COORD VECTORS HAVE'
    print *, '      UNEQUAL LENGTH!', size(x), size(y)
    call err_exit
  endif
  
  xm = pl_com%x_inch_to_mm * x
  ym = pl_com%y_inch_to_mm * y

  if (size(x) < 2) return
  call plline (size(x), xm, ym)

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
!                 TYPE is passed to GG_SETUP. E.g.
!                 TYPE = 'X'     --> Open an X-window.
!                 TYPE = 'GIF'   --> To create a gif file.
!                 TYPE = 'GIF-L' --> Gif w/ landscape page orientation.
!                 TYPE = 'PS'    --> To create a Color PostScript file.
!                 TYPE = 'PS-L'  --> PostScript w/ landscape page orientation.
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

subroutine qp_open_page_basic (page_type, x_len, y_len, plot_file, &
                                             x_page, y_page, i_chan, page_scale)

  implicit none

  type (pl_interface_struct), pointer :: pl_ptr

  real(rp) x_len, y_len, x_page, y_page, x1i, x2i, y1i, y2i, d, h
  real(rp), optional :: page_scale

  integer, optional :: i_chan
  integer ix, xp, yp, ix_len, iy_len, i_ch

  character(*) page_type, plot_file
  character(40) geom
  character(16) :: r_name = 'qp_open_page_basic'

  integer, parameter :: red(0:15) = (/255, 0, 255, 0, 0, 0, 255, 255, 255, 127&
                             , 0, 0, 127, 255, 85, 170/)
  integer, parameter :: green(0:15) = (/255, 0, 0, 255, 0, 255, 0, 255, 127,&
                             255, 255, 127, 0, 0, 85, 170/)
  integer, parameter :: blue(0:15) = (/255, 0, 0, 0, 255, 255, 255, 0, 0, 0,&
                             127, 255, 255, 127, 85, 170/)

! set plot type

  if (i_save == 0) then
    i_ch = 0
  else
    i_ch = pl_com%i_chan + 1
    call plsstrm(i_ch)
  endif


  if (page_type == 'X') then
    call plsdev ('xwin')

  elseif (page_type == 'TK') then
    call plsdev ('tk')

  elseif (page_type == 'PS') then
    call plsori(1)   ! portrait mode
    call plsdev ('psc')

  elseif (page_type == 'PS-L') then
    call plsdev ('psc')

  elseif (page_type == 'GIF') then
    call plsori(1)   ! portrait mode
    call plsdev ('png')

  elseif (page_type == 'GIF-L') then
    call plsdev ('png')

  else
    call out_io (s_abort$, r_name, 'ERROR: UNKNOWN PAGE_TYPE: ' // page_type)
    call err_exit
  endif

! Set output file name  

  if (page_type /= 'X' .and. page_type /= 'TK') then
     call plsfnam (trim(plot_file))
  endif

! Set color map
  call plscmap0(red, green, blue, 16)

! Set size of x-window.
! Work around for bug in plplot-5.9.5 is to set the geometry

  if (page_type == 'X' .or. page_type == 'TK') then
    ix_len = nint(85*x_len)
    iy_len = nint(85*y_len)
    call plspage (0.0_rp, 0.0_rp, ix_len, iy_len, 0, 0)
    write (geom, '(i0, a, i0, a)') ix_len, 'x', iy_len, '+10+10'
    call plsetopt ("geometry", trim(geom))
  endif

! Initialize plplot
  call plstar(1,1)
  call pladv(0)

! set viewport/window size

  call plvpor (0.0_rp, 1.0_rp, 0.0_rp, 1.0_rp)   ! viewport to entire page
  call plgspa(x1i, x2i, y1i, y2i)                ! Get current subpage in mm.
  call plwind (0.0_rp, x2i-x1i, 0.0_rp, y2i-y1i) ! set min/max for window

! Get page size info.
  call plgvpw(x1i,x2i,y1i,y2i)      ! Get viewport size in mm

  if (page_type == 'X' .or. page_type == 'TK') then
    !! call plschr(.7*point_to_mm_conv, 1.0_rp)
    call plschr(point_to_mm_conv, 1.0_rp)
  else
    call plschr(point_to_mm_conv, 1.0_rp)
  endif

  call plgchr(d, h)

! Remember plot area parameters

  if (present(i_chan)) i_chan = i_ch

  i_save = i_save + 1
  pl_com => pl_interface_save_com(i_save)

  pl_com%graph_pos%x1 = 0
  pl_com%graph_pos%x2 = x2i
  pl_com%graph_pos%y1 = 0
  pl_com%graph_pos%y2 = y2i
  pl_com%i_chan = i_ch
  pl_com%char_size = d*h
  pl_com%sym_size = 10
  pl_com%page_scale = 1    ! real_option(1.0_rp, page_scale)
  pl_com%page_type = page_type

  if (page_type == 'X' .or. page_type == 'TK') then
    pl_com%x_inch_to_mm = pl_com%page_scale * x2i / x_len
    pl_com%y_inch_to_mm = pl_com%page_scale * y2i / y_len
    x_page = x_len
    y_page = y_len
  else
    pl_com%x_inch_to_mm = pl_com%page_scale * 25.4
    pl_com%y_inch_to_mm = pl_com%page_scale * 25.4
    x_page = (x2i-x1i) / 25.4         ! convert to inches
    y_page = (y2i-y1i) / 25.4
  endif

  pl_com%x_inch_page = x_page
  pl_com%y_inch_page = y_page
  pl_com%page_type   = page_type

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
  !
  call out_io (s_abort$, 'qp_select_page_basic', 'NOT YET IMPLEMENTED!')
  call err_exit
  call plsstrm(iw)
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
  call plspause(.false.)  ! Disable end of page pause
  call plend1()
  i_save = i_save - 1
  pl_com => pl_interface_save_com(i_save)
  if (i_save /= 0) then
    call plsstrm(pl_com%i_chan)
  endif
end subroutine

end module
