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
! Note: PLPLOT uses real(4) units, QUICK_PLOT tries to interface between
!       real(4) and real(rp). Watch out for mixed units in subroutines.
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
! This is a list of all logicals which need to be set:
!   'uniform_size' in 'qp_set_symbol_size_basic'
!   'clip' in 'qp_set_clip_basic'
!-

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note: Needs scaling for PS !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module plplot_interface

use output_mod

type viewport_size
   real(rp) :: x1 !in mm
   real(rp) :: x2
   real(rp) :: y1
   real(rp) :: y2
end type

type pl_interface_struct
  type (viewport_size) :: graph_pos
  integer :: i_chan = -1
  integer :: fill_style
  integer :: line_width
  integer :: line_style
  real(rp) :: char_size
  integer :: fg_color
  logical :: clip
end type

type (pl_interface_struct), target, save :: pl_interface_com
type (pl_interface_struct), save :: pl_interface_save_com(10)
integer, save :: i_save = 0
real(rp), parameter :: point_to_mm_conv = .25   !approximate
private pl_interface_com, pl_interface_save_com, i_save

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
  real(rp) x1, x2, y1, y2, x1m, x2m, y1m, y2m, xp1, xp2, yp1, yp2

  x1m = x1 * 25.4
  x2m = x2 * 25.4
  y1m = y1 * 25.4
  y2m = y2 * 25.4

  call plsvpa (x1m, x2m, y1m, y2m)
  call plwind (x1m, x2m, y1m, y2m)

  pl_interface_com%graph_pos%x1 = x1m
  pl_interface_com%graph_pos%x2 = x2m
  pl_interface_com%graph_pos%y1 = y1m
  pl_interface_com%graph_pos%y2 = y2m

  if (.not. pl_interface_com%clip) then
     call plgspa(xp1, xp2, yp1, yp2)
     call plvpor(0.0_rp, 1.0_rp, 0.0_rp, 1.0_rp)
     call plwind(0.0_rp, xp2-xp1, 0.0_rp, yp2-yp1)
  endif

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_symbol_size_basic (height, symbol_type, page_type, uniform_size)
!
! Subroutine to set the symbol_size
!
! Input:
!   height       -- Real(rp): Symbol height.
!   symbol_type  -- Integer: Symbol type.
!   page_type    -- Character(*): Page type.
!   uniform_size -- Logical: Make all symbols the save size for constant height.
!+

subroutine qp_set_symbol_size_basic (height, symbol_type, page_type, uniform_size)

  implicit none

  real(rp) height, h, d, dum

  integer symbol_type

  logical uniform_size

  character(*) page_type

! The PLPLOT symbol set does not have a constant symbol size.
! This generally does not look nice so renormalize to get a consistant size.
! This excludes the set of circles with different sizes.
        
  h = height !* pl_interface_com%qp_to_pl_text_height_factor

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

    if (page_type == 'X') then
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

  call plschr(0.0_rp,h)
  call plgchr(d,dum)

  !Save this state
  pl_interface_com%char_size = h*d

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
  pl_interface_com%fill_style = fill

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
  pl_interface_com%line_width = line_width

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
  pl_interface_com%line_style = style

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

  if (.not. clip) then
     call plgspa(xp1, xp2, yp1, yp2)
     call plvpor(0.0_rp, 1.0_rp, 0.0_rp, 1.0_rp)
     call plwind(0.0_rp, xp2-xp1, 0.0_rp, yp2-yp1)

     pl_interface_com%clip = clip
  else
     call plsvpa (pl_interface_com%graph_pos%x1, pl_interface_com%graph_pos%x2, &
             pl_interface_com%graph_pos%y1, pl_interface_com%graph_pos%y2)
     call plwind (0.0_rp, pl_interface_com%graph_pos%x2-pl_interface_com%graph_pos%x1, &
             0.0_rp, pl_interface_com%graph_pos%y2-pl_interface_com%graph_pos%y1)

     pl_interface_com%clip = clip
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
  real(rp) height, d, h

  !call plschr(pl_interface_com%default_cs, height)
  call plschr(0.0_rp, height)
  call plgchr(d,h)

  !Save this state
  pl_interface_com%char_size = height*d

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

  real(rp) t_len, d, h
  integer len_text
  character(*) text

  call plgchr(d,h)

  t_len = .8*.03937*len_trim(text)*h  !This is kind-of a 1st order approx.
                            !since there is no subroutine for this action

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
  integer len_text, i
  real(rp) x0, y0, angle, justify, dx, dy, x0m, y0m, d, h, x1,x2,y1,y2, t_len
  real(rp), parameter :: pi=3.141592

  call plgvpw(x1,x2,y1,y2)
  call plgchr(d,h)
  t_len = len_trim(text)*h
  dx = cos(angle*pi/180)
  dy = sin(angle*pi/180)
  x0m = x0 * 25.4 - .5*h*dy  !x0, y0 specify coordinates of text baseline
  y0m = y0 * 25.4 + .5*h*dx  !but plptex needs the coordinates of the midline

  call plptex (x0m, y0m, dx, dy, justify, trim(text))

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
  
  xm = x * 25.4
  ym = y * 25.4

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
  real(rp) a,b,c,d
  integer style, fill, width, color

  i_save = i_save + 1

  pl_interface_save_com(i_save) = pl_interface_com
 
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
  real(rp) a,b,c,d,def,dum,xp1,xp2,yp1,yp2
  integer style, fill, width, color, state

  state = i_save
  
  call plgchr(def,dum)
  call plschr(def,(pl_interface_save_com(i_save)%char_size)/def)

  fill = pl_interface_save_com(i_save)%fill_style
  if(fill /= 0) then
     call plpsty(pl_interface_save_com(i_save)%fill_style)
  endif

  width = pl_interface_save_com(i_save)%line_width
  call plwid(pl_interface_save_com(i_save)%line_width)
  
  style = pl_interface_save_com(i_save)%line_style
  if(style /= 0) then
     call pllsty(pl_interface_save_com(i_save)%line_style)
  endif

  color = pl_interface_save_com(i_save)%fg_color
  call plcol0(pl_interface_save_com(i_save)%fg_color)

  a = pl_interface_save_com(i_save)%graph_pos%x1
  b = pl_interface_save_com(i_save)%graph_pos%x2
  c = pl_interface_save_com(i_save)%graph_pos%y1
  d = pl_interface_save_com(i_save)%graph_pos%y2

  if (pl_interface_save_com(i_save)%clip) then
     call plsvpa (a,b,c,d)
     call plwind (0.0_rp, b-a, 0.0_rp, d-c)
  else
     call plgspa(xp1, xp2, yp1, yp2)
     call plvpor(0.0_rp, 1.0_rp, 0.0_rp, 1.0_rp)
     call plwind(0.0_rp, xp2-xp1, 0.0_rp, yp2-yp1)
  endif
  
  pl_interface_com = pl_interface_save_com(i_save)
  i_save = i_save - 1
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_set_color_basic (ix_color, page_type)  
!
! Subroutine to set the color taking into accout that GIF
! inverts the black for white.
!
! Input:
!   ix_color -- Integer: Color index (0 - 15).
!   page_type -- Character(*): Type of page.
!-

subroutine qp_set_color_basic (ix_color, page_type)  

  implicit none

  integer ix_color
  integer, parameter :: inverse_color(0:15) = &
          (/ 1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 /)
!            1, 0, 5, 6, 7, 2, 3, 4, 11, 12, 13,  8,  9, 10, 15, 14 
!            0  1  2  3  4  5  6  7   8   9  10  11  12  13  14  15 

  character(*) page_type

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
     !pl_interface_com%fg_color = inverse_color(ix_color)
  !else
    call plcol0 (ix_color)
    !Save this state
     pl_interface_com%fg_color = ix_color
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
  call plclear()
end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine qp_clear_box_basic (x1, x2, y1, y2, page_type)
!
! Subroutine to clear all drawing from a box.
! That is, white out the box region.
! Units are inches from lower left of page.
!
! Input:
!   x1, y1 -- Real(rp): Bottom left corner of box.
!   x2, y2 -- Real(rp): Upper right corner of box.
!   page_type -- Character(*): Type of page.
!-

subroutine qp_clear_box_basic (x1, x2, y1, y2, page_type)
                
  implicit none

  real(rp) x1, x2, y1, y2, x1m, x2m, y1m, y2m
  real(rp) :: x_vec(0:3)
  real(rp) :: y_vec(0:3)
  character(*) page_type

  x1m = 25.4 * x1
  x2m = 25.4 * x2
  y1m = 25.4 * y1
  y2m = 25.4 * y2
!
  call qp_save_state_basic              ! Buffer the following calls

  call qp_set_color_basic(0, page_type) ! Set color index to background
  call plpsty(8)                        ! Set fill-area style to solid
  
  x_vec = (/x1m, x2m, x1m, x2m/)
  y_vec = (/y1m, y1m, y2m, y2m/)
  call plfill (4, x_vec, y_vec)         ! Fills a polygon with 4 vertices
  
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
  
  xm = 25.4 * x
  ym = 25.4 * y

  if (size(x) < 2) return
  call plline (size(x), xm, ym)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+                    
! Subroutine qp_open_page_basic (page_type, x_len, y_len, plot_file, &
!                                                x_page, y_page, i_chan)
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
!
! Output:
!   x_page    -- Real(rp): Horizontal page size in inches.
!   y_page    -- Real(rp): Vertical page size in inches.
!   i_chan    -- Integer, optional: Plot channel. 
!                 Like a unit number for a fortran OPEN.
!                 To be used with qp_select_page.
!-

subroutine qp_open_page_basic (page_type, x_len, y_len, plot_file, &
                                                       x_page, y_page, i_chan)

  implicit none

  real(rp) x_len, y_len, x_page, y_page, x1i, x2i, y1i, y2i, d, h
  integer, optional :: i_chan
  integer plsdev, iw, ix, xp, yp

  character(*) page_type, plot_file
  character(16) :: r_name = 'qp_open_page_basic'

  integer, parameter :: red(0:15) = (/255, 0, 255, 0, 0, 0, 255, 255, 255, 127&
                             , 0, 0, 127, 255, 85, 170/)
  integer, parameter :: green(0:15) = (/255, 0, 0, 255, 0, 255, 0, 255, 127,&
                             255, 255, 127, 0, 0, 85, 170/)
  integer, parameter :: blue(0:15) = (/255, 0, 0, 0, 255, 255, 255, 0, 0, 0,&
                             127, 255, 255, 127, 85, 170/)

! set plot type
  if (page_type == 'X') then
    iw = plsdev ('xwin')
  elseif (page_type == 'PS') then
    iw = plsdev ('psc')
  elseif (page_type == 'PS-L') then
    iw = plsdev ('psc')
  elseif (page_type == 'GIF') then
    iw = plsdev ('png')
  elseif (page_type == 'GIF-L') then
    iw = plsdev ('png')
  else
    call out_io (s_abort$, r_name, 'ERROR: UNKNOWN PAGE_TYPE: ' // page_type)
    call err_exit
  endif

! Set output file name  
  if (page_type /= 'X') then
     call plsetopt('-o', trim(plot_file))
  endif

  if (present(i_chan)) i_chan = iw

  if (iw <= 0) then
    print *, 'ERROR IN QP_OPEN_PAGE: CANNONT OPEN OUTPUT DEVICE!'
    stop
  endif

! Set color map
  call plscmap0(red, green, blue, 16)

! Set size of x-window
  if (page_type == 'X') then
     call plspage (0.0_rp, 0.0_rp, nint(100*x_len), nint(100*y_len), 0, 0)
  endif

! Initialize plplot
  call plstar(1,1)
  call pladv(0)

! set viewport/window size

! Adjusts page for portrait
  if (page_type == 'PS' .or. page_type == 'GIF') then
     call plsdiori(1.0_rp)
  endif

  call plvpor (0.0_rp, 1.0_rp, 0.0_rp, 1.0_rp)  ! viewport to entire page
  call plgspa(x1i,x2i,y1i,y2i)
  call plwind (0.0_rp, x2i-x1i, 0.0_rp, y2i-y1i) ! set min/max

! Get page size info.
  call plgvpw(x1i,x2i,y1i,y2i)      !in mm
  x_page = (x2i-x1i)*.03937               !convert to inches
  y_page = (y2i-y1i)*.03937

  if (page_type == 'X') then
     call plschr(.7*point_to_mm_conv, 1.0_rp)
  else
     call plschr(point_to_mm_conv, 1.0_rp)
  endif

  call plgchr(d,h)

! Remember plot area parameters
  pl_interface_com%graph_pos%x1 = 0
  pl_interface_com%graph_pos%x2 = x2i
  pl_interface_com%graph_pos%y1 = 0
  pl_interface_com%graph_pos%y2 = y2i

  call qp_set_clip_basic(.false.)

  pl_interface_com%i_chan = iw
  pl_interface_com%char_size = d*h

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
  call pladv(iw)
  do i = 1, size(pl_interface_save_com)
    if (pl_interface_save_com(i)%i_chan == iw) then
      pl_interface_save_com(i_save) = pl_interface_save_com(i)
      pl_interface_save_com(i) = pl_interface_com
      pl_interface_com = pl_interface_save_com(i_save)
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
  call plend()
  if (i_save /= 0) pl_interface_com = pl_interface_save_com(i_save)
end subroutine

end module
