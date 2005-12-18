!+
! Subroutine tao_init_plotting (plot_file)
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
!-

subroutine tao_init_plotting (plot_file)

use tao_mod
use tao_input_struct
use quick_plot
use tao_plot_window_mod

implicit none

type (tao_plot_page_struct), pointer :: page
type (tao_plot_struct), pointer :: plt
type (tao_graph_struct), pointer :: grph
type (tao_curve_struct), pointer :: crv
type (tao_plot_input) plot
type (tao_graph_input) graph
type (tao_plot_page_input) plot_page
type (tao_region_input) region(n_region_maxx)
type (tao_curve_input) curve(n_curve_maxx)
type (tao_place_input) place(10)
type (tao_ele_shape_struct) shape(20)
type (qp_symbol_struct) default_symbol
type (qp_line_struct) default_line
type (qp_axis_struct) init_axis

integer iu, i, j, ip, n, ng, ios, i_uni
integer graph_index
integer, allocatable :: ix_ele(:)

character(200) file_name, plot_file
character(20) :: r_name = 'tao_init_plotting'

logical lat_layout_here

namelist / tao_plot_page / plot_page, region, place
namelist / tao_template_plot / plot
namelist / tao_template_graph / graph, graph_index, curve
namelist / element_shapes / shape

! See if this routine has been called before

if (.not. s%global%init_plot_needed) return
s%global%init_plot_needed = .false.

! Read in the plot page parameters

call tao_open_file ('TAO_INIT_DIR', plot_file, iu, file_name)
call out_io (s_blank$, r_name, '*Init: Opening Plotting File: ' // file_name)

place%region = ' '
region%name = ' '       ! a region exists only if its name is not blank 
plot_page%title(:)%draw_it = .false.
plot_page%title(:)%string = ' '
plot_page%title(:)%justify = 'CC'
plot_page%title(:)%x = 0.50
plot_page%title(:)%y = 0.990
plot_page%title(1)%y = 0.996
plot_page%title(2)%y = 0.97
plot_page%title(:)%units = '%PAGE'
plot_page%text_scale = 1
call out_io (s_blank$, r_name, 'Init: Reading tao_plot_page namelist')
read (iu, nml = tao_plot_page, err = 9000)

page => s%plot_page
page%size = plot_page%size
page%border = plot_page%border
page%text_scale = plot_page%text_scale

! title

page%text_height = plot_page%text_height
page%title = plot_page%title
forall (i = 1:size(page%title), (page%title(i)%string .ne. ' ')) &
            page%title(i)%draw_it = .true.

! allocate a s%plot_page%plot structure for each region defined and
! transfer the info from the input region structure.

n = count(region%name /= ' ')
allocate (page%region(n))

do i = 1, n
  page%region(i)%name     = region(i)%name
  page%region(i)%location = region(i)%location
enddo

! Read in the plot templates and transfer the info to the 
! s%tamplate_plot structures

ip = 0   ! number of template plots
lat_layout_here = .false.

do
  plot%name = ' '
  plot%x_axis_type = 'index'
  plot%x = init_axis
  plot%independent_graphs = .false.
  read (iu, nml = tao_template_plot, iostat = ios, err = 9100)  
  if (ios /= 0) exit                                 ! exit on end of file.
  call out_io (s_blank$, r_name, &
                  'Init: Read tao_template_plot namelist: ' // plot%name)
  ip = ip + 1

  if (ip .gt. n_template_maxx) then
    call out_io (s_warn$, r_name, &
            "Number of plot templates exceeds maximum of \I2\ ", n_template_maxx)
    call out_io (s_blank$, r_name, &
                "Only first \I2\ will be used", n_template_maxx)
    exit
  endif
  
  plt => s%template_plot(ip)
  plt%name        = plot%name
  plt%x           = plot%x
  plt%x_divisions = plot%x%major_div
  plt%x_axis_type = plot%x_axis_type
  plt%independent_graphs = plot%independent_graphs

  ng = plot%n_graph
  if (ng == 0) then
    nullify (plt%graph)
  else
    allocate (plt%graph(ng))
  endif

  do i = 1, ng
    write (graph%name, '(i0)') i
    graph_index = 0                 ! setup defaults
    graph%type  = 'data'
    graph%legend_origin = qp_point_struct(1.0_rp, 1.0_rp, '%GRAPH')
    graph%y  = init_axis
    graph%y2 = init_axis
    graph%y2%draw_numbers = .false.
    graph%ix_universe = 0
    graph%clip = .true.
    graph%who%name  = ' '                               ! set default
    graph%who(1) = tao_plot_who_struct('model', +1)     ! set default
    do j = 1, size(curve)
      write (curve(j)%name, '(i0)') j
    enddo
    curve(:)%units_factor = 1
    curve(:)%convert = .false.                             ! set default
    curve(:)%symbol_every = 1
    curve(:)%ix_universe = 0
    curve(:)%draw_line = .true.
    curve(:)%draw_symbols = .true.
    curve(:)%use_y2 = .false.
    curve(:)%symbol = default_symbol
    curve(:)%line   = default_line
    curve(:)%ele2_name   = ' '
    curve(:)%ix_ele2 = -1
    curve(2:7)%symbol%type = (/ times$, square$, plus$, triangle$, &
                                  x_symbol$, diamond$ /)
    curve(2:7)%symbol%color = (/ blue$, red$, green$, cyan$, magenta$, yellow$ /)
    curve(2:7)%line%color = curve(2:7)%symbol%color
!    curve(2:7)%line%style = (/ dashed$, dotted$, dash_dot$, &
!                                               dash_dot3$, solid$, dotted$ /)
    read (iu, nml = tao_template_graph, err = 9200)
    call out_io (s_blank$, r_name, &
                 'Init: Read tao_template_graph namelist: ' // graph%name)
    if (graph_index /= i) then
      call out_io (s_error$, r_name, &
                                  'BAD "GRAPH_INDEX" FOR: ' // graph%name)
      call err_exit
    endif
    grph => plt%graph(i)
    grph%name          = graph%name
    grph%type          = graph%type
    grph%who           = graph%who
    grph%legend_origin = graph%legend_origin
    grph%box           = graph%box
    grph%title         = graph%title
    grph%margin        = graph%margin
    grph%y             = graph%y
    grph%y2            = graph%y2
    grph%ix_universe   = graph%ix_universe
    grph%clip          = graph%clip
    grph%title_suffix = ' '
    grph%legend = ' '

    if (grph%ix_universe < 0 .or. grph%ix_universe > size(s%u)) then
      call out_io (s_error$, r_name, 'UNIVERSE INDEX: \i4\ ', grph%ix_universe)
      call out_io (s_blank$, r_name, &
       'OUT OF RANGE FOR PLOT:GRAPH: ' // trim(plot%name) // ':' // graph%name)
      call err_exit
    endif

    if (grph%type == 'lat_layout') then
      lat_layout_here = .true.
      if (plt%x_axis_type /= 's') call out_io (s_error$, r_name, &
                            'A lat_layout must have x_axis_type = "s" for a visible plot!')
    endif

    if (graph%n_curve == 0) then
      nullify (grph%curve)
    else
      allocate (grph%curve(graph%n_curve))
    endif

    do j = 1, graph%n_curve
      crv => grph%curve(j)
      crv%data_source       = curve(j)%data_source
      crv%data_type         = curve(j)%data_type
      crv%units_factor      = curve(j)%units_factor
      crv%symbol_every      = curve(j)%symbol_every
      crv%ix_universe       = curve(j)%ix_universe
      crv%draw_line         = curve(j)%draw_line
      crv%draw_symbols      = curve(j)%draw_symbols
      crv%use_y2            = curve(j)%use_y2
      crv%symbol            = curve(j)%symbol
      crv%line              = curve(j)%line
      crv%convert           = curve(j)%convert
      crv%name              = curve(j)%name
      crv%ele2_name         = curve(j)%ele2_name
      call str_upcase (crv%ele2_name, crv%ele2_name)
      crv%ix_ele2           = curve(j)%ix_ele2

      i_uni = crv%ix_universe
      if (i_uni == 0) i_uni = s%global%u_view
 
      ! Find the ele2 info if either ele2_name or ix_ele2 has been set.
      ! If plotting something like the phase then the default is for ele2 to be the beginning element.

      if (crv%ele2_name == ' ' .and. crv%ix_ele2 >= 0) then ! if ix_ele2 has been set ...
        crv%ele2_name = s%u(i_uni)%design%ele_(crv%ix_ele2)%name ! then find the name
      elseif (crv%ele2_name /= ' ') then                    ! if ele2_name has been set ...
        call tao_locate_element (crv%ele2_name, i_uni, ix_ele, .true.) ! then find the index
        crv%ix_ele2 = ix_ele(1)
      elseif (crv%data_type(1:5) == 'phase' .or. crv%data_type(1:2) == 'r:' .or. &
              crv%data_type(1:2) == 't:' .or. crv%data_type(1:3) == 'tt:') then
        crv%ix_ele2 = 0
        crv%ele2_name = s%u(i_uni)%design%ele_(0)%name
      endif

    enddo
  enddo
enddo

! read in shapes

s%plot_page%ele_shape%key = 0

if (lat_layout_here) then

  rewind (iu)
  shape(:)%key_name = ' '
  shape(:)%key = 0
  read (iu, nml = element_shapes, iostat = ios)

  if (ios /= 0) then
    call out_io (s_error$, r_name, 'ERROR READING ELE_SHAPE NAMELIST IN FILE.')
    call err_exit
  endif

  do i = 1, size(shape)
    call str_upcase (shape(i)%key_name, shape(i)%key_name)
    call str_upcase (shape(i)%ele_name, shape(i)%ele_name)
    call str_upcase (shape(i)%shape,    shape(i)%shape)
    call str_upcase (shape(i)%color,    shape(i)%color)

    if (shape(i)%key_name == ' ') cycle

    do j = 1, n_key
      if (shape(i)%key_name == key_name(j)) then
        shape(i)%key = j
        exit
      endif
    enddo          

    if (shape(i)%key == 0) then
      print *, 'ERROR: CANNOT FIND KEY FOR: ', shape(i)%key_name
      call err_exit
    endif

  enddo
  s%plot_page%ele_shape = shape

endif

close (1)

! initial placement of plots

do i = 1, size(place)
  if (place(i)%region == ' ') cycle
  call tao_place_cmd (place(i)%region, place(i)%plot)
enddo

call tao_create_plot_window
if (s%global%plot_on) call qp_set_qp_parameters (text_scale = page%text_scale)

return

!-----------------------------------------
! Error handling

9000 continue
call out_io (s_error$, r_name, &
        'TAO_PLOT_PAGE NAMELIST READ ERROR.', 'IN FILE: ' // file_name)
rewind (iu)
do
  read (iu, nml = tao_plot_page)  ! force printing of error message
enddo

!-----------------------------------------

9100 continue
call out_io (s_error$, r_name, &
        'TAO_TEMPLATE_PLOT NAMELIST READ ERROR.', 'IN FILE: ' // file_name)
rewind (iu)
do
  read (iu, nml = tao_template_plot)  ! force printing of error message
enddo

!-----------------------------------------

9200 continue
call out_io (s_error$, r_name, &
       'TAO_TEMPLATE_GRAPH NAMELIST READ ERROR.', 'IN FILE: ' // file_name)
rewind (iu)
do
  read (iu, nml = tao_template_graph)  ! force printing of error message
enddo

end subroutine tao_init_plotting
