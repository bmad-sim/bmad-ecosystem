!+
! Subroutine tao_init_plotting (s, plot_file)
!
! Subroutine to initialize the tao plotting structures.
! If plot_file is not in the current directory then it will be searched
! for in the directory:
!   TAO_INIT_DIR
!
! Input:
!   plot_file -- Character(*): Plot initialization file.
!
! Output:
!   s -- Tao_super_universe_struct:
!-

subroutine tao_init_plotting (s, plot_file)

use tao_mod
use tao_input_struct
use quick_plot

implicit none

type (tao_super_universe_struct), target :: s
type (tao_plot_page_struct), pointer :: page
type (tao_plot_struct), pointer :: plt
type (tao_graph_struct), pointer :: grph
type (tao_curve_struct), pointer :: crv
type (tao_plot_input) plot
type (tao_graph_input) graph
type (tao_plot_page_input) plot_page
type (tao_plot_region_struct) region(n_region_maxx)
type (tao_curve_input) curve(n_curve_maxx)
type (tao_place_input) place(10)

integer iu, i, j, ip, n, ng, ios
integer graph_index

character(200) namelist_file, plot_file
character(20) :: r_name = 'tao_init_plotting'

namelist / tao_plot_page / plot_page, region, place
namelist / tao_template_plot / plot
namelist / tao_template_graph / graph, graph_index, curve

! Read in the plot page parameters

call tao_open_file ('TAO_INIT_DIR', plot_file, iu, namelist_file)

place%region = ' '
region%name = ' '       ! a region exists only if its name is not blank 
read (iu, nml = tao_plot_page, err = 9000)

page => s%plot_page
page%size = plot_page%size
page%border = plot_page%border

! Open a plot window

call qp_open_page ('X', page%id_window, page%size(1), page%size(2), 'POINTS')
call qp_set_layout (page_border = page%border)

!set default font size

call qp_set_text_attrib ("MAIN_TITLE",  height = plot_page%text_height) 
call qp_set_text_attrib ("GRAPH_TITLE",  height = plot_page%text_height) 
call qp_set_text_attrib ("AXIS_LABEL",  height = plot_page%text_height) 

! allocate a s%plot_page%plot structure for each region defined and
! transfer the info from the input region structure.

n = count(region%name /= ' ')
allocate (page%plot(n))

do i = 1, n
  page%plot(i)%region = region(i)
enddo

! Read in the plot templates and transfer the info to the 
! s%tamplate_plot structures

ip = 0   ! number of template plots
do
  plot%who%name  = ' '                               ! set default
  plot%who(1) = tao_plot_who_struct('model', +1)     ! set default
  plot%who(2) = tao_plot_who_struct('design', -1)    ! set default
  plot%convert = .false.                             ! set default
  read (iu, nml = tao_template_plot, iostat = ios, err = 9100)  
  if (ios /= 0) exit                                 ! exit on end of file.
  ip = ip + 1
  plt => s%template_plot(ip)
  plt%class      = plot%class
  plt%box_layout = plot%box_layout
  plt%x          = plot%x
  plt%who        = plot%who
  plt%convert    = plot%convert
  ng = plot%n_graph
  allocate (plt%graph(ng))
  do i = 1, ng
    graph_index = 0
    curve(:)%symbol_every = 1
    curve(:)%ix_universe = 0
    curve(:)%draw_line = .true.
    read (iu, nml = tao_template_graph, err = 9200)
    if (graph_index /= i) then
      call out_io (s_error$, r_name, &
                                  'BAD "GRAPH_INDEX" FOR: ' // graph%name)
      call err_exit
    endif
    grph => plt%graph(i)
    grph%name       = graph%name
    grph%this_box   = graph%this_box
    grph%title      = graph%title
    grph%margin     = graph%margin
    grph%y          = graph%y
    allocate (grph%curve(graph%n_curve))
    do j = 1, graph%n_curve
      crv => grph%curve(j)
      crv%sub_class          = curve(j)%sub_class
      crv%units_factor      = curve(j)%units_factor
      crv%symbol_every      = curve(j)%symbol_every
      crv%ix_universe       = curve(j)%ix_universe
      crv%draw_line         = curve(j)%draw_line
    enddo
  enddo
enddo

close (1)

! initial placement of plots

do i = 1, size(place)
  if (place(i)%region == ' ') cycle
  call tao_place_cmd (s, place(i)%region, place(i)%plot)
enddo

return

!-----------------------------------------
! Error handling

9000 continue
call out_io (s_error$, r_name, &
        'TAO_PLOT_PAGE NAMELIST READ ERROR.', 'IN FILE: ' // namelist_file)
rewind (iu)
do
  read (iu, nml = tao_plot_page)  ! force printing of error message
enddo

!-----------------------------------------

9100 continue
call out_io (s_error$, r_name, &
        'TAO_TEMPLATE_PLOT NAMELIST READ ERROR.', 'IN FILE: ' // namelist_file)
rewind (iu)
do
  read (iu, nml = tao_template_plot)  ! force printing of error message
enddo

!-----------------------------------------

9200 continue
call out_io (s_error$, r_name, &
       'TAO_TEMPLATE_GRAPH NAMELIST READ ERROR.', 'IN FILE: ' // namelist_file)
rewind (iu)
do
  read (iu, nml = tao_template_graph)  ! force printing of error message
enddo


end subroutine tao_init_plotting
