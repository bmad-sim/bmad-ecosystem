!+
! Subroutine tao_init_plotting (plot_file_in)
!
! Subroutine to initialize the tao plotting structures.
! If plot_file is not in the current directory then it will be searched
! for in the directory:
!   TAO_INIT_DIR
!
! Input:
!   plot_file_in -- Character(*): Plot initialization file.
!
! Output:
!-

subroutine tao_init_plotting (plot_file_in)

use tao_mod
use tao_input_struct
use quick_plot
use tao_plot_window_mod

implicit none

type old_tao_ele_shape_struct    ! for the element layout plot
  character(40) key_name     ! Element key name
  character(40) ele_name     ! element name
  character(16) shape        ! plot shape
  character(16) color        ! plot color
  real(rp) dy_pix            ! plot vertical height 
  Logical :: draw_name = .true.
  integer key                ! Element key index to match to
end type

type (tao_plot_page_struct), pointer :: page
type (tao_plot_page_struct) plot_page, plot_page_default
type (tao_plot_struct), pointer :: plt
type (tao_graph_struct), pointer :: grph
type (tao_curve_struct), pointer :: crv
type (tao_plot_input) plot, default_plot
type (tao_graph_input) graph, default_graph, master_default_graph
type (tao_region_input) region(n_region_maxx)
type (tao_curve_input) curve(n_curve_maxx)
type (tao_place_input) place(10)
type (old_tao_ele_shape_struct) shape(20)
type (tao_ele_shape_struct) ele_shape(20)
type (qp_symbol_struct) default_symbol
type (qp_line_struct) default_line
type (qp_axis_struct) init_axis
type (ele_pointer_struct), allocatable, save :: eles(:)

real(rp) shape_height_max, y1, y2

integer iu, i, j, k, ix, ip, n, ng, ios, i_uni
integer graph_index, color, i_graph

character(*) plot_file_in
character(len(plot_file_in)) plot_file_array
character(100) plot_file, graph_name, file_name
character(80) label
character(20) :: r_name = 'tao_init_plotting'

namelist / tao_plot_page / plot_page, default_plot, default_graph, region, place
namelist / tao_template_plot / plot, default_graph
namelist / tao_template_graph / graph, graph_index, curve
namelist / element_shapes / shape
namelist / element_shapes_floor_plan / ele_shape
namelist / element_shapes_lat_layout / ele_shape


logical err

! See if this routine has been called before

call qp_init_com_struct()  ! Init quick_plot
if (.not. s%global%init_plot_needed) return
s%global%init_plot_needed = .false.

! Init

init_axis%min = 0
init_axis%max = 0

place%region = ' '
region%name = ' '       ! a region exists only if its name is not blank 
plot_page = plot_page_default
plot_page%title(:)%draw_it = .false.
plot_page%title(:)%string = ' '
plot_page%title(:)%justify = 'CC'
plot_page%title(:)%x = 0.50
plot_page%title(:)%y = 0.990
plot_page%title(1)%y = 0.996
plot_page%title(2)%y = 0.97
plot_page%title(:)%units = '%PAGE'
plot_page%size = [600, 800]

default_plot%name = ' '
default_plot%x_axis_type = 'index'
default_plot%x = init_axis
default_plot%x%minor_div_max = 6
default_plot%x%major_div_nominal = -1
default_plot%x%major_div = -1
default_plot%independent_graphs = .false.
default_plot%autoscale_gang_x = .true.
default_plot%autoscale_gang_y = .true.
default_plot%autoscale_x = .false.
default_plot%autoscale_y = .false.
default_plot%n_graph = 0

default_graph%title                 = ''
default_graph%type                  = 'data'
default_graph%text_legend_origin    = qp_point_struct(5.0_rp, 0.0_rp, 'POINTS/GRAPH/RT')
default_graph%curve_legend_origin   = qp_point_struct(5.0_rp, -2.0_rp, 'POINTS/GRAPH/LT')
default_graph%y                     = init_axis
default_graph%y%major_div           = -1
default_graph%y%major_div_nominal   = -1
default_graph%y2                    = init_axis
default_graph%y2%major_div          = -1
default_graph%y2%major_div_nominal  = -1
default_graph%y2%label_color        = blue$
default_graph%y2%draw_numbers       = .false.
default_graph%ix_universe           = -1
default_graph%clip                  = .true.
default_graph%draw_axes             = .true.
default_graph%correct_xy_distortion = .false.
default_graph%draw_curve_legend     = .true.
default_graph%component             = 'model'
default_graph%who%name              = ''
default_graph%who%sign              = 1
default_graph%box     = [1, 1, 1, 1]
default_graph%margin  = qp_rect_struct(0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, '%GRAPH')
default_graph%n_curve = 0

! If there is no plot file then use the built-in defaults.

if (plot_file_in == '') then
  call tao_setup_default_plotting()
  return
endif

! Read in the plot page parameters
! plot_file_in may contain multiple file names separated by spaces.

plot_file_array = plot_file_in
call string_trim(plot_file_array, plot_file_array, ix)
plot_file = plot_file_array(1:ix)

call out_io (s_blank$, r_name, '*Init: Opening Plotting File: ' // plot_file)
call tao_open_file ('TAO_INIT_DIR', plot_file, iu, file_name)
if (iu == 0) then
  call out_io (s_fatal$, r_name, 'ERROR OPENING PLOTTING FILE. WILL EXIT HERE...')
  call err_exit
endif

call out_io (s_blank$, r_name, 'Init: Reading tao_plot_page namelist')
read (iu, nml = tao_plot_page, iostat = ios)
if (ios > 0) then
  call out_io (s_error$, r_name, 'ERROR READING TAO_PLOT_PAGE NAMELIST.')
  rewind (iu)
  read (iu, nml = tao_plot_page)  ! To give error message
endif
if (ios < 0) call out_io (s_blank$, r_name, 'Note: No tao_plot_page namelist found')

master_default_graph = default_graph

s%plot_page = plot_page
page => s%plot_page

! title

forall (i = 1:size(page%title), (page%title(i)%string .ne. ' ')) &
            page%title(i)%draw_it = .true.

! allocate a s%plot_page%plot structure for each region defined and
! transfer the info from the input region structure.

n = count(region%name /= ' ')
allocate (s%plot_region(n))

do i = 1, n
  s%plot_region(i)%name     = region(i)%name
  s%plot_region(i)%location = region(i)%location
enddo

!-----------------------------------------------------------------------------------
! Read in element shapes

rewind (iu)
shape(:)%key_name = ''
shape(:)%ele_name = ''
read (iu, nml = element_shapes, iostat = ios)

if (ios > 0) then
  call out_io (s_error$, r_name, 'ERROR READING ELEMENT_SHAPES NAMELIST IN FILE.')
  rewind (iu)
  read (iu, nml = element_shapes)  ! To generate error message
endif

if (ios == 0) then
!    call out_io (s_error$, r_name, 'DEPRECATED ELEMENT_SHAPES NAMELIST DETECTED.', &
!                                   'PLEASE CHANGE TO THE NEW SYNTAX (SEE THE MANUAL)!')
  do i = 1, size(shape)
    ele_shape(i)%ele_name = shape(i)%ele_name
    if (shape(i)%key_name /= '') &
            ele_shape(i)%ele_name = trim(shape(i)%key_name) // '::' // ele_shape(i)%ele_name
    ele_shape(i)%shape     = shape(i)%shape
    ele_shape(i)%color     = shape(i)%color
    ele_shape(i)%dy_pix    = shape(i)%dy_pix
    if (shape(i)%draw_name) then
      ele_shape(i)%label_type = 'name'
    else
      ele_shape(i)%label_type = 'none'
    endif
  enddo

  call tao_uppercase_shapes (n)
  allocate (tao_com%ele_shape_floor_plan(n), tao_com%ele_shape_lat_layout(n))
  tao_com%ele_shape_floor_plan = ele_shape(1:n)
  tao_com%ele_shape_lat_layout = ele_shape(1:n)
endif


if (ios < 0) then

  rewind (iu)
  ele_shape(:)%ele_name = ' '
  ele_shape(:)%label_type = 'name'
  read (iu, nml = element_shapes_floor_plan, iostat = ios)

  if (ios > 0) then
    call out_io (s_error$, r_name, 'ERROR READING ELEMENT_SHAPES_FLOOR_PLAN NAMELIST')
    rewind (iu)
    read (iu, nml = element_shapes_floor_plan)  ! To generate error message
  endif

  call tao_uppercase_shapes (n)
  allocate (tao_com%ele_shape_floor_plan(n))
  tao_com%ele_shape_floor_plan = ele_shape(1:n)

  rewind (iu)
  ele_shape(:)%ele_name = ' '
  read (iu, nml = element_shapes_lat_layout, iostat = ios)

  if (ios > 0) then
    call out_io (s_error$, r_name, 'ERROR READING ELEMENT_SHAPES_LAT_LAYOUT NAMELIST')
    rewind (iu)
    read (iu, nml = element_shapes_lat_layout)  ! To generate error message
  endif

  call tao_uppercase_shapes (n)
  allocate (tao_com%ele_shape_lat_layout(n))
  tao_com%ele_shape_lat_layout = ele_shape(1:n)

endif

! Error check

if (allocated(tao_com%ele_shape_lat_layout)) then
  do i = 1, size(tao_com%ele_shape_lat_layout)
    select case (tao_com%ele_shape_lat_layout(i)%shape)
    case ('BOX', 'VAR_BOX', 'ASYM_VAR_BOX', 'XBOX', 'DIAMOND', 'BOW_TIE', 'CIRCLE', 'X')
    case default
      print *, 'ERROR: UNKNOWN ELE_SHAPE: ', tao_com%ele_shape_lat_layout(i)%shape
      call err_exit
    end select
  enddo
endif

if (allocated(tao_com%ele_shape_floor_plan)) then
  do i = 1, size(tao_com%ele_shape_floor_plan)
    select case (tao_com%ele_shape_floor_plan(i)%shape)
    case ('BOX', 'VAR_BOX', 'ASYM_VAR_BOX', 'XBOX', 'DIAMOND', 'BOW_TIE', 'CIRCLE', 'X')
    case default
      print *, 'ERROR: UNKNOWN ELE_SHAPE: ', tao_com%ele_shape_floor_plan(i)%shape
      call err_exit
    end select
  enddo
endif

close (iu)

!------------------------------------------------------------------------------------
! Read in the plot templates and transfer the info to the 
! s%tamplate_plot structures

! First count the number of plots needed

ip = 0   ! number of template plots
plot_file_array = plot_file_in

do   ! Loop over plot files

  call string_trim(plot_file_array, plot_file_array, ix)
  if (ix == 0) exit
  plot_file = plot_file_array(1:ix)
  plot_file_array = plot_file_array(ix+1:)
  call tao_open_file ('TAO_INIT_DIR', plot_file, iu, file_name)
  if (iu == 0) then
    call out_io (s_fatal$, r_name, 'ERROR OPENING PLOTTING FILE: ' // plot_file, &
                                   'WILL EXIT HERE...')
    call err_exit
  endif

  do   ! Loop over templates in a file
    read (iu, nml = tao_template_plot, iostat = ios, err = 9100)  
    call out_io (s_blank$, r_name, 'Init: Read tao_template_plot ' // plot%name)
    if (ios /= 0) exit
    ip = ip + 1
  enddo

  close (iu)
enddo

allocate (s%template_plot(ip))

! Now read in the plots

ip = 0   ! number of template plots
plot_file_array = plot_file_in

do  ! Loop over plot files

  call string_trim(plot_file_array, plot_file_array, ix)
  if (ix == 0) exit
  plot_file = plot_file_array(1:ix)
  plot_file_array = plot_file_array(ix+1:)
  call out_io (s_blank$, r_name, '*Init: Opening Plotting File: ' // plot_file)
  call tao_open_file ('TAO_INIT_DIR', plot_file, iu, file_name)

  do   ! Loop over templates in a file

    plot = default_plot
    default_graph = master_default_graph

    read (iu, nml = tao_template_plot, iostat = ios, err = 9100)  
    if (ios /= 0) exit

    call out_io (s_blank$, r_name, 'Init: Read tao_template_plot namelist: ' // plot%name)
    do i = 1, ip
      if (plot%name == s%template_plot(ip)%name) then
        call out_io (s_error$, r_name, 'DUPLICATE PLOT NAME: ' // plot%name)
        exit
      endif
    enddo

    ip = ip + 1

    plt => s%template_plot(ip)
    nullify(plt%r)
    plt%name                 = plot%name
    plt%x_axis_type          = plot%x_axis_type
    plt%x                    = plot%x
    plt%autoscale_gang_x     = plot%autoscale_gang_x 
    plt%autoscale_gang_y     = plot%autoscale_gang_y 
    plt%autoscale_x          = plot%autoscale_x 
    plt%autoscale_y          = plot%autoscale_y 

    if (plt%x%major_div < 0 .and. plt%x%major_div_nominal < 0) plt%x%major_div_nominal = 6

    if (plot%independent_graphs) then  ! Old style
      call out_io (s_error$, r_name, [&
            '**********************************************************', &
            '**********************************************************', &
            '**********************************************************', &
            '***** SYNTAX CHANGE:                                 *****', &
            '*****     PLOT%INDEPENDENT_GRAPHS = True             *****', &
            '***** NEEDS TO BE CHANGED TO:                        *****', &
            '*****     PLOT%AUTOSCALE_GANG_Y = False              *****', &
            '***** TAO WILL RUN NORMALLY FOR NOW...               *****', &
            '***** SEE THE TAO MANUAL FOR MORE DETAILS!           *****', &
            '**********************************************************', &
            '**********************************************************', &
            '**********************************************************'] )
      plt%autoscale_gang_y = .false.
    endif

    call qp_calc_axis_places (plt%x)

    do
      ix = index(plt%name, '.')
      if (ix == 0) exit
      call out_io (s_error$, r_name, 'PLOT NAME HAS ".": ' // plt%name, &
                   'SUBSTITUTING "-"')
      plt%name(ix:ix) = '-'
    enddo

    ng = plot%n_graph
    if (allocated(plt%graph)) deallocate (plt%graph)
    if (ng /= 0) allocate (plt%graph(ng))

    do i_graph = 1, ng
      graph_index = 0         ! setup defaults
      graph = default_graph
      graph%x = plot%x
      graph%x_axis_scale_factor = 1
      write (graph%name, '(a, i0)') 'g', i_graph
      do j = 1, size(curve)
        write (curve(j)%name, '(a, i0)') 'c', j
      enddo
      curve(:)%data_source = 'lat'
      curve(:)%data_index  = ''
      curve(:)%data_type_x = ''
      curve(:)%data_type   = ''
      curve(:)%y_axis_scale_factor = 1
      curve(:)%symbol_every = 1
      curve(:)%ix_universe = -1
      curve(:)%ix_branch = 0
      curve(:)%draw_line = .true.
      if (plt%x_axis_type == 's' .or. plt%x_axis_type == 'lat' .or. plt%x_axis_type == 'var') then
        curve(:)%draw_symbols = .false.
      else
        curve(:)%draw_symbols = .true.
      endif
      curve(:)%draw_symbol_index = .false.
      curve(:)%use_y2 = .false.
      curve(:)%symbol = default_symbol
      curve(:)%line   = default_line
      curve(:)%ele_ref_name   = ' '
      curve(:)%ix_ele_ref = -1
      curve(:)%smooth_line_calc = .true.
      curve(:)%draw_interpolated_curve = .true.
      curve(:)%line%width = -1
      curve(:)%legend_text = ''
      curve(:)%x_axis_scale_factor = 0  ! This is old syntax. Not used.
      curve(2:7)%symbol%type = [times_sym$, square_sym$, &
                  plus_sym$, triangle_sym$, x_symbol_sym$, diamond_sym$]
      curve(2:7)%symbol%color = &
                  [blue$, red$, green$, cyan$, magenta$, yellow$]
      curve(2:7)%line%color = curve(2:7)%symbol%color
      read (iu, nml = tao_template_graph, err = 9200)
      call out_io (s_blank$, r_name, 'Init: Read tao_template_graph ' // graph%name)
      graph_name = trim(plot%name) // '.' // graph%name
      call out_io (s_blank$, r_name, &
              'Init: Read tao_template_graph namelist: ' // graph_name)
      if (graph_index /= i_graph) then
        call out_io (s_error$, r_name, &
              'BAD "GRAPH_INDEX" FOR PLOT: ' // plot%name, &
              'LOOKING FOR GRAPH_INDEX: \I0\ ', &
              'BUT TAO_TEMPLACE_GRAPH HAD GRAPH_INDEX: \I0\ ', &
              i_array = [i_graph, graph_index] )
        call err_exit
      endif
      grph => plt%graph(i_graph)
      grph%p             => plt
      grph%name          = graph%name
      grph%type          = graph%type
      grph%component     = graph%component
      grph%x_axis_scale_factor = graph%x_axis_scale_factor 
      if (graph%who(1)%name /= '') then  ! Old style
        call out_io (s_error$, r_name, [&
            '**********************************************************', &
            '**********************************************************', &
            '**********************************************************', &
            '***** SYNTAX CHANGE:          GRAPH%WHO              *****', &
            '***** NEEDS TO BE CHANGED TO: GRAPH%COMPONENT        *****', &
            '***** EXAMPLE:                                       *****', &
            '*****   GRAPH%WHO(1) = "MODEL", +1                   *****', &
            '*****   GRAPH%WHO(2) = "DESIGN", -1                  *****', &
            '***** GETS CHANGED TO:                               *****', &
            '*****   GRAPH%COMPONENT = "MODEL - DESIGN"           *****', &
            '***** TAO WILL RUN NORMALLY FOR NOW...               *****', &
            '***** SEE THE TAO MANUAL FOR MORE DETAILS!           *****', &
            '**********************************************************', &
            '**********************************************************', &
            '**********************************************************'] )
        grph%component = graph%who(1)%name
        do i = 2, size(graph%who)
          if (graph%who(i)%name == '') exit
          if (nint(graph%who(i)%sign) == 1) then
            grph%component = trim(grph%component) // ' + ' // graph%who(i)%name
          elseif (nint(graph%who(i)%sign) == -1) then
            grph%component = trim(grph%component) // ' - ' // graph%who(i)%name
          else
            call out_io (s_fatal$, r_name, 'BAD "WHO" IN PLOT TEMPLATE: ' // plot%name)
            call err_exit
          endif
        enddo
      endif
      grph%text_legend_origin    = graph%text_legend_origin
      grph%curve_legend_origin   = graph%curve_legend_origin
      grph%box                   = graph%box
      grph%title                 = graph%title
      grph%margin                = graph%margin
      grph%x                     = graph%x
      grph%y                     = graph%y
      grph%y2                    = graph%y2
      grph%ix_universe           = graph%ix_universe
      grph%ix_branch             = graph%ix_branch
      grph%clip                  = graph%clip
      grph%draw_axes             = graph%draw_axes
      grph%correct_xy_distortion = graph%correct_xy_distortion
      grph%draw_curve_legend     = graph%draw_curve_legend
      grph%title_suffix          = ''
      grph%text_legend           = ''
      grph%y2_mirrors_y          = .true.
      if (grph%x%major_div < 0 .and. grph%x%major_div_nominal < 0) grph%x%major_div_nominal = 6
      if (grph%y%major_div < 0 .and. grph%y%major_div_nominal < 0) grph%y%major_div_nominal = 4
      if (grph%y2%major_div < 0 .and. grph%y2%major_div_nominal < 0) grph%y2%major_div_nominal = 4

      call qp_calc_axis_places (grph%x)

      do
        ix = index(grph%name, '.')
        if (ix == 0) exit
        call out_io (s_error$, r_name, 'GRAPH NAME HAS ".": ' // grph%name, &
                     'SUBSTITUTING "-"')
        grph%name(ix:ix) = '-'
      enddo

      call qp_calc_axis_places (grph%y)

      if (.not. tao_com%common_lattice .and. grph%ix_universe == 0) then
        call out_io (s_error$, r_name, [&
            '**********************************************************', &
            '**********************************************************', &
            '**********************************************************', &
            '***** SYNTAX CHANGE: GRAPH%IX_UNIVERSE = 0           *****', &
            '***** NEEDS TO BE CHANGED TO: GRAPH%IX_UNIVERSE = -1 *****', &
            '**********************************************************', &
            '**********************************************************', &
            '**********************************************************'] )
        grph%ix_universe = -1
      endif

      if (grph%ix_universe < -1 .or. grph%ix_universe > ubound(s%u, 1)) then
        call out_io (s_error$, r_name, 'UNIVERSE INDEX: \i4\ ', & 
                                       'OUT OF RANGE FOR PLOT:GRAPH: ' // graph_name, &
                                       i_array = [grph%ix_universe] )
        call err_exit
      endif

      if (grph%type == 'floor_plan' .and. .not. allocated (tao_com%ele_shape_floor_plan)) &
                call out_io (s_error$, r_name, 'NO ELEMENT SHAPES DEFINED FOR FLOOR_PLAN PLOT.')
   
      if (grph%type == 'lat_layout') then
        if (.not. allocated (tao_com%ele_shape_lat_layout)) call out_io (s_error$, r_name, &
                              'NO ELEMENT SHAPES DEFINED FOR LAT_LAYOUT PLOT.')
        if (plt%x_axis_type /= 's') call out_io (s_error$, r_name, &
                              'A LAT_LAYOUT MUST HAVE X_AXIS_TYPE = "s" FOR A VISIBLE PLOT!')
      endif

      if (graph%n_curve == 0) then
        if (allocated(grph%curve)) deallocate (grph%curve)
      else
        allocate (grph%curve(graph%n_curve))
      endif

      do j = 1, graph%n_curve
        crv => grph%curve(j)
        crv%g                    => grph
        crv%data_source          = curve(j)%data_source
        crv%data_index           = curve(j)%data_index
        crv%data_type_x          = curve(j)%data_type_x
        crv%data_type            = curve(j)%data_type
        crv%y_axis_scale_factor  = curve(j)%y_axis_scale_factor
        crv%symbol_every         = curve(j)%symbol_every
        crv%ix_universe          = curve(j)%ix_universe
        crv%draw_line            = curve(j)%draw_line
        crv%draw_symbols         = curve(j)%draw_symbols
        crv%draw_symbol_index    = curve(j)%draw_symbol_index
        crv%use_y2               = curve(j)%use_y2
        crv%symbol               = curve(j)%symbol
        crv%line                 = curve(j)%line
        crv%smooth_line_calc     = curve(j)%smooth_line_calc
        crv%name                 = curve(j)%name
        crv%ele_ref_name         = curve(j)%ele_ref_name
        call str_upcase (crv%ele_ref_name, crv%ele_ref_name)
        crv%ix_ele_ref           = curve(j)%ix_ele_ref
        crv%ix_bunch             = curve(j)%ix_bunch
        crv%ix_branch            = curve(j)%ix_branch
        crv%legend_text          = curve(j)%legend_text

        ! Convert old syntax to new

        if (crv%data_source == 'beam_tracking') crv%data_source = 'beam'
        if (crv%data_source == 'lattice')       crv%data_source = 'lat'
        if (crv%data_source == 'data_array')    crv%data_source = 'dat'
        if (crv%data_source == 'var_array')     crv%data_source = 'var'

        if (curve(j)%x_axis_scale_factor /= 0) then
          call out_io (s_error$, r_name, [&
            '**********************************************************', &
            '**********************************************************', &
            '**********************************************************', &
            '***** SYNTAX CHANGE:                                 *****', &
            '*****         CURVE%X_AXIS_SCALE_FACTOR              *****', &
            '***** NEEDS TO BE CHANGED TO:                        *****', &
            '*****         GRAPH%X_AXIS_SCALE_FACTOR              *****', &
            '***** TAO WILL RUN NORMALLY FOR NOW...               *****', &
            '**********************************************************', &
            '**********************************************************', &
            '**********************************************************'] )
          crv%smooth_line_calc = .false.
        endif

        if (.not. curve(j)%draw_interpolated_curve) then
          call out_io (s_error$, r_name, [&
            '**********************************************************', &
            '**********************************************************', &
            '**********************************************************', &
            '***** SYNTAX CHANGE:                                 *****', &
            '*****         CURVE%DRAW_INTERPOLATED_CURVE          *****', &
            '***** NEEDS TO BE CHANGED TO:                        *****', &
            '*****         CURVE%SMOOTH_LINE_CALC                 *****', &
            '***** TAO WILL RUN NORMALLY FOR NOW...               *****', &
            '**********************************************************', &
            '**********************************************************', &
            '**********************************************************'] )
          crv%smooth_line_calc = .false.
        endif

        ! Convert old syntax to new

        ix = index(crv%data_type, 'emittance.')
        if (ix /= 0) crv%data_type = crv%data_type(1:ix-1) // 'emit.' // crv%data_type(ix+10:)

        ! Default data type

        if (crv%data_type == '') crv%data_type = trim(plt%name) // '.' // trim(grph%name)

        ! A dot in the name is verboten.
        do
          ix = index(crv%name, '.')
          if (ix == 0) exit
          call out_io (s_error$, r_name, 'CURVE NAME HAS ".": ' // crv%name, &
                       'SUBSTITUTING "-"')
          crv%name(ix:ix) = '-'
        enddo

        ! Convert old style phase_space data_type to new style

        if (grph%type == 'phase_space') then
          ix = index(crv%data_type, '-')
          if (ix /= 0 .and. crv%data_type_x == '') then
            crv%data_type_x = crv%data_type(1:ix-1)
            crv%data_type   = crv%data_type(ix+1:)
          endif
        endif  

        ! Turn on the y2 axis numbering if needed.

        if (crv%use_y2) then
          grph%y2%draw_numbers = .true.
          grph%y2_mirrors_y = .false.
          grph%y2%label_color = crv%symbol%color
        endif

        ! Set curve line width

        if (crv%line%width == -1) then
          if (plt%x_axis_type == 's') then
            crv%line%width = 2
          else
            crv%line%width = 1
          endif
        endif

        ! Enable the radiation integrals calculation if needed.

        if (.not. tao_com%common_lattice .and. crv%ix_universe == 0) then
          call out_io (s_error$, r_name, [&
            '**********************************************************', &
            '**********************************************************', &
            '**********************************************************', &
            '***** SYNTAX CHANGE: CURVE%IX_UNIVERSE = 0           *****', &
            '***** NEEDS TO BE CHANGED TO: CURVE%IX_UNIVERSE = -1 *****', &
            '**********************************************************', &
            '**********************************************************', &
            '**********************************************************'] )
          crv%ix_universe = -1
        endif


        i_uni = tao_universe_number (crv%ix_universe)
        if (i_uni > ubound(s%u, 1)) then
          call out_io (s_error$, r_name, &
                          'CURVE OF PLOT: ' // plot%name, &
                          'HAS UNIVERSE INDEX OUT OF RANGE: \I0\ ', &
                          i_array = [i_uni] )
          call err_exit
        endif

        ! Find the ele_ref info if either ele_ref_name or ix_ele_ref has been set.
        ! If plotting something like the phase then the default is for ele_ref 
        ! to be the beginning element.

        ! if ix_ele_ref has been set ...
        if (crv%ele_ref_name == ' ' .and. crv%ix_ele_ref >= 0) then 
          crv%ele_ref_name = s%u(i_uni)%design%lat%ele(crv%ix_ele_ref)%name ! find the name
        ! if ele_ref_name has been set ...
        elseif (crv%ele_ref_name /= ' ') then
          call tao_locate_elements (crv%ele_ref_name, i_uni, eles, err, .true.) ! find the index
          crv%ix_ele_ref = eles(1)%ele%ix_ele
          crv%ix_branch  = eles(1)%ele%ix_branch
        elseif (crv%data_type(1:5) == 'phase' .or. crv%data_type(1:2) == 'r.' .or. &
                crv%data_type(1:2) == 't.' .or. crv%data_type(1:3) == 'tt.') then
          crv%ix_ele_ref = 0
          crv%ele_ref_name = s%u(i_uni)%design%lat%ele(0)%name
        elseif (graph%type == 'phase_space') then
          plt%x_axis_type = 'phase_space'
          crv%ix_ele_ref = 0
          crv%ele_ref_name = s%u(i_uni)%design%lat%ele(0)%name
        elseif (graph%type == 'key_table') then
          plt%x_axis_type = 'none'
        elseif (graph%type == 'floor_plan') then
          plt%x_axis_type = 'floor'
        endif

        call tao_ele_to_ele_track (i_uni, crv%ix_branch, crv%ix_ele_ref, crv%ix_ele_ref_track)

      enddo  ! curve

      call qp_calc_axis_places (grph%y2)
      if (grph%y2%min == grph%y2%max .and. .not. grph%y2_mirrors_y) then
        label = grph%y2%label
        color = grph%y2%label_color
        grph%y2 = grph%y
        grph%y2%label = label
        grph%y2%label_color = color
      endif
    enddo  ! graph
  enddo  ! plot

  close(iu)

enddo  ! file

close (iu)

! initial placement of plots

do i = 1, size(place)
  if (place(i)%region == ' ') cycle
  call tao_place_cmd (place(i)%region, place(i)%plot)
enddo

call tao_create_plot_window

return

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

!----------------------------------------------------------------------------------------
contains

subroutine tao_uppercase_shapes (n_shape)

integer n, n_shape

!

n_shape = 0
do n = 1, size(ele_shape)
  ! Bmad wants ele names upper case but Tao data is case sensitive.
  if (ele_shape(n)%ele_name(1:5) /= 'dat::' .and. ele_shape(n)%ele_name(1:5) /= 'var::') &
                     call str_upcase (ele_shape(n)%ele_name, ele_shape(n)%ele_name)
  call str_upcase (ele_shape(n)%shape,    ele_shape(n)%shape)
  call str_upcase (ele_shape(n)%color,    ele_shape(n)%color)
  call downcase_string (ele_shape(n)%label_type)
  if (ele_shape(n)%label_type == '') ele_shape(n)%label_type = 'name'
  if (index('false', trim(ele_shape(n)%label_type)) == 1) ele_shape(n)%label_type = 'none'
  if (index('true', trim(ele_shape(n)%label_type)) == 1) ele_shape(n)%label_type = 'name'
  ! Convert old class:name format to new class::name format
  ix = index(ele_shape(n)%ele_name, ":")
  if (ix /= 0 .and. ele_shape(n)%ele_name(ix+1:ix+1) /= ':') &
     ele_shape(n)%ele_name = ele_shape(n)%ele_name(1:ix) // ':' // ele_shape(n)%ele_name(ix+1:)

  if (ele_shape(n)%ele_name /= '') n_shape = n
enddo

end subroutine

!----------------------------------------------------------------------------------------
! contains

subroutine tao_setup_default_plotting()

real(rp) y_top

!

s%plot_page = plot_page

allocate (tao_com%ele_shape_floor_plan(10), tao_com%ele_shape_lat_layout(10))

tao_com%ele_shape_floor_plan(:)%ele_name = ''
tao_com%ele_shape_floor_plan(1:5) = [&
          tao_ele_shape_struct('SBEND::*',      'BOX',  'BLUE',    08.0_rp, 'none'), &
          tao_ele_shape_struct('QUADRUPOLE::*', 'XBOX', 'MAGENTA', 15.0_rp, 'name'), &
          tao_ele_shape_struct('SEXTUPOLE::*',  'XBOX', 'GREEN',   15.0_rp, 'none'), &
          tao_ele_shape_struct('LCAVITY::*',    'XBOX', 'RED',     20.0_rp, 'none'), &
          tao_ele_shape_struct('RFCAVITY::*',   'XBOX', 'RED',     20.0_rp, 'none')]

tao_com%ele_shape_lat_layout = tao_com%ele_shape_floor_plan

! beta

allocate (s%template_plot(5))
plt => s%template_plot(1)

nullify(plt%r)
if (allocated(plt%graph)) deallocate (plt%graph)
allocate (plt%graph(2))
plt%graph(1)%p => plt
plt%graph(2)%p => plt
allocate (plt%graph(1)%curve(1))
allocate (plt%graph(2)%curve(1))

plt%name                 = 'beta'
plt%x_axis_type          = 's'
plt%x                    = init_axis
plt%x%major_div_nominal  = 8
plt%x%minor_div_max = 6
plt%autoscale_gang_x = .true.
plt%autoscale_gang_y = .true.


grph => plt%graph(1)
grph%name          = 'a'
grph%title         = 'Horizontal Beta'
grph%type          = 'data'
grph%margin        =  qp_rect_struct(0.15, 0.06, 0.12, 0.12, '%BOX')
grph%box           = [1, 2, 1, 2]
grph%y             = init_axis
grph%y%label       = '\gb\dA\u'
grph%y%major_div_nominal   = 4
grph%y2%draw_numbers = .false.
grph%component     = 'model'
crv => grph%curve(1)
crv%data_source = 'lat'
crv%draw_symbols = .false.
crv%data_type = 'beta.a'

grph => plt%graph(2)
grph               = plt%graph(1)
grph%name          = 'b'
grph%title         = 'Vertical Beta'
grph%y%label       = '\gb\dB\u'
grph%box           = [1, 1, 1, 2]
crv => grph%curve(1)
crv%data_type = 'beta.b'

! eta

plt => s%template_plot(2)

nullify(plt%r)
if (allocated(plt%graph)) deallocate (plt%graph)
allocate (plt%graph(2))
plt%graph(1)%p => plt
plt%graph(2)%p => plt
allocate (plt%graph(1)%curve(1))
allocate (plt%graph(2)%curve(1))

plt = s%template_plot(1)
plt%name           = 'eta'

grph => plt%graph(1)
grph%name          = 'x'
grph%title         = 'Horizontal Eta'
grph%y%label       = '\gy\dX\u'
grph%y%major_div_nominal   = 4
crv => grph%curve(1)
crv%data_type = 'eta.x'

grph => plt%graph(2)
grph%name          = 'y'
grph%title         = 'Vertical Eta'
grph%y%label       = '\gy\dY\u'
crv => grph%curve(1)
crv%data_type = 'eta.y'

! Orbit

plt => s%template_plot(3)

nullify(plt%r)
if (allocated(plt%graph)) deallocate (plt%graph)
allocate (plt%graph(2))
plt%graph(1)%p => plt
plt%graph(2)%p => plt
allocate (plt%graph(1)%curve(1))
allocate (plt%graph(2)%curve(1))

plt = s%template_plot(1)
plt%name           = 'orbit'

grph => plt%graph(1)
grph%name          = 'x'
grph%title         = 'Horizontal Orbit'
grph%y%label       = 'X'
grph%y%major_div_nominal   = 4
crv => grph%curve(1)
crv%data_type = 'orbit.x'

grph => plt%graph(2)
grph%name          = 'y'
grph%title         = 'Vertical Orbit'
grph%y%label       = 'Y'
crv => grph%curve(1)
crv%data_type = 'orbit.y'

! Lat Layout

plt => s%template_plot(4)

nullify(plt%r)
if (allocated(plt%graph)) deallocate (plt%graph)
allocate (plt%graph(1))
plt%graph(1)%p => plt

plt%name           = 'lat_layout'
plt%x_axis_type    = 's'
plt%x              = init_axis

grph => plt%graph(1)
grph%name          = 'g1'
grph%box           = [1, 1, 1, 1]
grph%type          = 'lat_layout'
grph%margin        =  qp_rect_struct(0.15, 0.06, 0.12, 0.12, '%BOX')
grph%x             = init_axis

! Floor Plan

plt => s%template_plot(5)

nullify(plt%r)
if (allocated(plt%graph)) deallocate (plt%graph)
allocate (plt%graph(1))
plt%graph(1)%p => plt

plt%name           = 'floor_plan'
plt%x_axis_type          = 'floor'
plt%x                    = init_axis

grph => plt%graph(1)
grph%name          = 'g1'
grph%box           = [1, 1, 1, 1]
grph%type          = 'floor_plan'
grph%margin        =  qp_rect_struct(0.15, 0.06, 0.12, 0.12, '%BOX')
grph%correct_xy_distortion = .true.
grph%x             = init_axis
grph%y             = init_axis

! Regions

allocate (s%plot_region(20))

y_top = 0.85
s%plot_region(1)%name = 'r_top'
s%plot_region(1)%location = [0.0_rp, 1.0_rp, y_top, 1.00_rp]

k = 1
do i = 1, 4
  do j = 1, i
    k = k + 1
    write (s%plot_region(k)%name, '(a, 2i0)') 'r', j, i
    y1 = y_top * real(j-1)/ i
    y2 = y_top * real(j) / i
    s%plot_region(k)%location = [0.0_rp, 1.0_rp, y1, y2]
  enddo
enddo

call tao_place_cmd ('r_top', 'lat_layout')
call tao_place_cmd ('r12', 'beta')
call tao_place_cmd ('r22', 'eta')

end subroutine

end subroutine tao_init_plotting
