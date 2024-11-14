!+
! Subroutine tao_init_plotting (plot_file_in)
!
! Subroutine to initialize the tao plotting structures.
!
! Input:
!   plot_file_in -- Character(*): Plot initialization file.
!
! Output:
!-

subroutine tao_init_plotting (plot_file_in)

use tao_input_struct
use tao_plot_window_mod, dummy => tao_init_plotting
use quick_plot

implicit none

interface
  function old_style_title_syntax(iu) result (is_old_style)
  integer iu
  logical is_old_style
  end function
end interface

integer, parameter :: n_region_maxx   = 100     ! number of plotting regions.
integer, parameter :: n_curve_maxx    = 40      ! number of curves per graph

type old_tao_ele_shape_struct    ! for the element layout plot
  character(40) key_name     ! Element key name
  character(40) ele_name     ! element name
  character(16) shape        ! plot shape
  character(16) color        ! plot color
  real(rp) size              ! plot vertical height 
  Logical :: draw_name  = .true.
  integer key                ! Element key index to match to
end type

type (tao_plot_page_input) plot_page, plot_page_default
type (tao_plot_struct), pointer :: plt
type (tao_graph_struct), pointer :: grph
type (tao_curve_struct), pointer :: crv
type (tao_plot_input) plot, default_plot
type (tao_graph_input) :: graph, default_graph, master_default_graph
type (tao_region_input) region(n_region_maxx)
type (tao_curve_input) curve(n_curve_maxx), default_curve
type (tao_place_input) place(30)
type (old_tao_ele_shape_struct) shape(30)
type (tao_ele_shape_input) ele_shape(60)
type (tao_ele_shape_struct), pointer :: e_shape
type (ele_pointer_struct), allocatable :: eles(:)
type (qp_axis_struct) init_axis

real(rp) y1, y2
real(rp), parameter :: re_g$ = real_garbage$

integer iu, i, j, k, k1, k2, ix, ip, n, ng, ios, ios1, ios2, i_uni
integer graph_index, i_graph, ic, n_curve, nc
integer, parameter :: int_g$ = int_garbage$

character(*) plot_file_in
character(len(plot_file_in)) plot_file_array
character(200) plot_file, full_file_name
character(100) graph_name
character(80) label
character(40) str, old_name
character(16) color
character(16), parameter :: symbol_type(8) = [character(16):: 'circle_dot', 'times', 'square', 'plus', &
                                                  'triangle', 'x_symbol', 'diamond', 'diamond']
character(16), parameter :: draw_color(8) = [character(16):: 'black', 'blue', 'red', 'green', 'cyan', &
                                                   'magenta', 'yellow', 'black']
character(*), parameter :: r_name = 'tao_init_plotting'

logical err, include_default_plots, include_default_shapes, include_dflt_lat_layout, include_dflt_floor_plan

namelist / tao_plot_page / plot_page, default_plot, default_graph, region, place, include_default_plots
namelist / tao_template_plot / plot, default_graph, default_curve
namelist / tao_template_graph / graph, graph_index, curve

namelist / floor_plan_drawing / ele_shape, include_default_shapes
namelist / lat_layout_drawing / ele_shape, include_default_shapes

! These are old style

namelist / element_shapes / shape
namelist / element_shapes_floor_plan / ele_shape
namelist / element_shapes_lat_layout / ele_shape

! See if this routine has been called before

if (.not. allocated(s%plot_page%lat_layout%ele_shape)) allocate(s%plot_page%lat_layout%ele_shape(0))
if (.not. allocated(s%plot_page%floor_plan%ele_shape)) allocate(s%plot_page%floor_plan%ele_shape(0))

call qp_init_com_struct()  ! Init quick_plot
if (.not. s%com%init_plot_needed) return
s%com%init_plot_needed = .false.

! Init

init_axis%min = 0
init_axis%max = 0

place%region = ''
region%name  = ''       ! a region exists only if its name is not blank 
include_default_plots = .true.

plot_page = plot_page_default
plot_page%title%y = 0.99
plot_page%subtitle%y = 0.97
plot_page%size = [500, 600]
plot_page%border = qp_rect_struct(0.001_rp, 0.001_rp, 0.001_rp, 0.001_rp, '%PAGE')

default_plot%name = ''
default_plot%description = ''
default_plot%x_axis_type = 'index'
default_plot%x = init_axis
default_plot%x%minor_div_max = 6
default_plot%x%major_div_nominal = -1
default_plot%x%major_div = -1
default_plot%autoscale_gang_x = .true.
default_plot%autoscale_gang_y = .true.
default_plot%autoscale_x = .false.
default_plot%autoscale_y = .false.
default_plot%n_graph = 0
default_plot%n_curve_pts = -1

default_graph = tao_graph_input()
default_graph%x                     = init_axis
default_graph%x%major_div           = -1
default_graph%x%major_div_nominal   = -1
default_graph%y                     = init_axis
default_graph%y%major_div           = -1
default_graph%y%major_div_nominal   = -1
default_graph%x2                    = init_axis
default_graph%x2%major_div          = -1
default_graph%x2%major_div_nominal  = -1
default_graph%x2%draw_label         = .false.
default_graph%y2                    = init_axis
default_graph%y2%major_div          = -1
default_graph%y2%major_div_nominal  = -1
default_graph%y2%label_color        = null_name$
default_graph%y2%draw_numbers       = .false.
default_graph%margin                = qp_rect_struct(0.15, 0.06, 0.12, 0.12, '%BOX')

! If there is no plot file then use the built-in defaults.

if (allocated(s%plot_page%pattern)) deallocate(s%plot_page%pattern)
allocate (s%plot_page%pattern(0))

if (plot_file_in == '') then
  call tao_setup_default_plotting(.false., .false., default_graph)
  call number_template_plots()
  return
endif

! Read in the plot page parameters
! plot_file_in may contain multiple file names separated by spaces.

plot_file_array = plot_file_in
call string_trim(plot_file_array, plot_file_array, ix)
plot_file = plot_file_array(1:ix)

call out_io (s_blank$, r_name, '*Init: Opening Plotting File: ' // plot_file)
call tao_open_file (plot_file, iu, full_file_name, s_error$)
if (iu == 0) then
  call out_io (s_fatal$, r_name, 'ERROR OPENING PLOTTING FILE: ' // plot_file)
  call tao_setup_default_plotting(.false., .false., default_graph)
  call number_template_plots()
  return
endif

call out_io (s_blank$, r_name, 'Init: Reading tao_plot_page namelist')
read (iu, nml = tao_plot_page, iostat = ios)
if (ios > 0) then
  ! First check to see if this is due to old style "title(i)" syntax.
  if (old_style_title_syntax(iu)) then
    call out_io (s_error$, r_name, 'ERROR READING TAO_PLOT_PAGE NAMELIST IN FILE:' // plot_file, &
                                   '"PLOT_PAGE%TITLE(1) = ..." IS REPLACED BY "PLOT_PAGE%TITLE = ..." and', &
                                   '"PLOT_PAGE%TITLE(2) = ..." IS REPLACED BY "PLOT_PAGE%SUBTITLE = ...".')
  else
    call out_io (s_error$, r_name, 'ERROR READING TAO_PLOT_PAGE NAMELIST IN FILE:' // plot_file)
  endif
  rewind (iu)
  read (iu, nml = tao_plot_page)  ! To give error message
endif
if (ios < 0) call out_io (s_blank$, r_name, 'Note: No tao_plot_page namelist found')

master_default_graph = default_graph

call tao_set_plotting (plot_page, s%plot_page, .true.)

! Plot window geometry specified on cmd line?

if (s%init%geometry_arg /= '') then
  str = s%init%geometry_arg
  ix = index(str, 'x')
  if (ix == 0) then
    call out_io (s_error$, r_name, 'MALFORMED -geometry ARGUMENT. NO "x" PRESENT: ' // str, 'IN FILE: ' // plot_file)
  else
    if (.not. is_integer(str(1:ix-1)) .or. .not. is_integer(str(ix+1:))) then
      call out_io (s_error$, r_name, 'MALFORMED -geometry ARGUMENT: ' // str, 'IN FILE: ' // plot_file)
    else
      read (str(:ix-1), *) plot_page%size(1)
      read (str(ix+1:), *) plot_page%size(2)
    endif
  endif
endif
 
! allocate a s%plot_page%plot structure for each region defined and
! transfer the info from the input region structure.

n = count(region%name /= '')
allocate (s%plot_page%region(n))

do i = 1, n
  s%plot_page%region(i)%name      = region(i)%name
  s%plot_page%region(i)%location  = region(i)%location
  s%plot_page%region(i)%plot%name = ''
enddo

!------------------------------------------------------------------------------------
! Read in shape patterns
! Reason this is in a separate routine is due to conflict with "curve" variable in namelist.

call tao_read_in_patterns(iu, plot_file)

!-----------------------------------------------------------------------------------
! Read in element shapes...
! First look for old style namelist 

rewind (iu)
ele_shape(:) = tao_ele_shape_input()

shape(:)%key_name = ''
shape(:)%ele_name = ''

read (iu, nml = element_shapes, iostat = ios)

if (ios > 0) then
  call out_io (s_error$, r_name, 'ERROR READING ELEMENT_SHAPES NAMELIST.', 'IN FILE: ' // plot_file)
  rewind (iu)
  read (iu, nml = element_shapes)  ! To generate error message
endif

if (ios == 0) then
  call out_io (s_error$, r_name, 'ELEMENT_SHAPES NAMELIST IS DEPRECATED.', &
                                 'PLEASE CONVERT TO FLOOR_PLAN_DRAWING AND/OR LAT_LAYOUT_DRAWING NAMELISTS.', 'IN FILE: ' // plot_file)
  do i = 1, size(shape)
    ele_shape(i)%ele_id = shape(i)%ele_name
    if (shape(i)%key_name /= '') shape(i)%ele_name = trim(shape(i)%key_name) // '::' // ele_shape(i)%ele_id
    ele_shape(i)%shape      = shape(i)%shape
    ele_shape(i)%color      = shape(i)%color
    ele_shape(i)%size       = shape(i)%size
    ele_shape(i)%draw       = .true.
    if (shape(i)%draw_name) then
      ele_shape(i)%label = 'name'
    else
      ele_shape(i)%label = 'none'
    endif
  enddo

  call tao_transfer_shape (ele_shape, s%plot_page%lat_layout%ele_shape, 'ELEMENT_SHAPES')
  call tao_transfer_shape (ele_shape, s%plot_page%floor_plan%ele_shape, 'ELEMENT_SHAPES')
endif

! Look for new style shape namelist if could not find old style

if (ios < 0) then

  ! Read floor_plan_drawing namelist

  include_default_shapes = .false.
  rewind (iu)
  read (iu, nml = element_shapes_floor_plan, iostat = ios1)  ! Deprecated name
  rewind (iu)
  read (iu, nml = floor_plan_drawing, iostat = ios2)
  include_dflt_floor_plan = include_default_shapes

  if (ios1 >= 0) then
    call out_io (s_error$, r_name, &
            'Note: The "element_shapes_floor_plan" namelist has been renamed to', &
            '      "floor_plan_drawing" to reflect the fact that this namelist ', &
            '      now is used to specify more than element shapes. Please     ', &
            '      make the appropriate change in your input file...           ')
  endif

  if (ios1 > 0) then 
    rewind (iu)
    call out_io (s_error$, r_name, 'ERROR READING ELEMENT_SHAPES_FLOOR_PLAN NAMELIST', 'IN FILE: ' // plot_file)
    read (iu, nml = element_shapes_floor_plan)  ! To generate error message
  endif

  if (ios2 > 0) then 
    rewind (iu)
    call out_io (s_error$, r_name, 'ERROR READING FLOOR_PLAN_DRAWING NAMELIST', 'IN FILE: ' // plot_file)
    read (iu, nml = floor_plan_drawing)
  endif

  call tao_transfer_shape (ele_shape, s%plot_page%floor_plan%ele_shape, 'FLOOR_PLAN_DRAWING')

  ! Read element_shapes_lat_layout namelist

  ele_shape(:) = tao_ele_shape_input()

  include_default_shapes = .false.
  rewind (iu)
  read (iu, nml = element_shapes_lat_layout, iostat = ios1)
  rewind (iu)
  read (iu, nml = lat_layout_drawing, iostat = ios2)
  include_dflt_lat_layout = include_default_shapes

  if (ios1 == 0) then
    call out_io (s_error$, r_name, &
            'Note: The "element_shapes_lattice_list" namelist has been renamed to', &
            '      "lat_layout_drawing" to reflect the fact that this namelist   ', &
            '      now is used to specify more than element shapes. Please       ', &
            '      make the appropriate change in your input file...             ')
  endif

  if (ios1 == 0) then
    ele_shape(:)%size = ele_shape(:)%size * 1.0 / 40.0 ! scale to current def.
  endif 

  if (ios1 > 0) then 
    rewind (iu)
    call out_io (s_error$, r_name, 'ERROR READING ELEMENT_SHAPES_LAT_LAYOUT NAMELIST', 'IN FILE: ' // plot_file)
    read (iu, nml = element_shapes_lat_layout)  ! To generate error message
  endif

  if (ios2 > 0) then 
    rewind (iu)
    call out_io (s_error$, r_name, 'ERROR READING LAT_LAYOUT_DRAWING NAMELIST', 'IN FILE: ' // plot_file)
    read (iu, nml = lat_layout_drawing)
  endif

  call tao_transfer_shape (ele_shape, s%plot_page%lat_layout%ele_shape, 'LAT_LAYOUT_DRAWING')
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
  call tao_open_file (plot_file, iu, full_file_name, s_fatal$)

  do   ! Loop over templates in a file
    read (iu, nml = tao_template_plot, iostat = ios)
    if (ios > 0) then
      call out_io (s_error$, r_name, &
                'TAO_TEMPLATE_PLOT NAMELIST READ ERROR.', 'IN FILE: ' // full_file_name, &
                'THIS IS THE ' // ordinal_str(ip+1) // ' TAO_TEMPLATE_PLOT NAMELIST IN THE FILE.')

      rewind (iu)
      do
        read (iu, nml = tao_template_plot)  ! force printing of error message
      enddo
      return
    endif

    call out_io (s_blank$, r_name, 'Init: Read tao_template_plot ' // quote(plot%name))
    if (ios /= 0) exit
    ip = ip + 1
  enddo

  close (iu)
enddo

!---------------
! Allocate the template plot and define a scratch plot

allocate (s%plot_page%template(ip+1))

plt => s%plot_page%template(ip+1)

nullify(plt%r)
if (allocated(plt%graph)) deallocate (plt%graph)
allocate (plt%graph(1))
plt%graph(1)%p => plt
plt%name = 'scratch'
plt%graph(1)%name = 'g'
plt%graph(1)%p => plt

! Now read in the plots

ip = 0   ! template plot index
plot_file_array = plot_file_in

do  ! Loop over plot files

  call string_trim(plot_file_array, plot_file_array, ix)
  if (ix == 0) exit
  plot_file = plot_file_array(1:ix)
  plot_file_array = plot_file_array(ix+1:)
  call out_io (s_blank$, r_name, 'Init: Opening another plotting file: ' // plot_file)
  call tao_open_file (plot_file, iu, full_file_name, s_fatal$)

  do   ! Loop over templates in a file

    plot = default_plot
    default_graph = master_default_graph
    default_curve = tao_curve_input()
    default_curve%symbol%type = ''
    default_curve%symbol%color = ''
    default_curve%line%color = ''

    plot%x = qp_axis_struct(null_name$, re_g$, re_g$, re_g$, re_g$, re_g$, re_g$, re_g$, re_g$, re_g$, re_g$, re_g$, &
                null_name$, int_g$, int_g$, int_g$, int_g$, int_g$, null_name$, null_name$, &
                int_g$, int_g$, .true., .true.)

    read (iu, nml = tao_template_plot, iostat = ios)  
    if (ios /= 0) exit

    call out_io (s_blank$, r_name, 'Init: Read tao_template_plot namelist: ' // plot%name)
    do i = 1, ip
      if (plot%name == s%plot_page%template(ip)%name) then
        call out_io (s_error$, r_name, 'DUPLICATE PLOT NAME: ' // plot%name, 'IN FILE: ' // plot_file)
        exit
      endif
    enddo

    ip = ip + 1

    plt => s%plot_page%template(ip)
    nullify(plt%r)
    plt%name                 = plot%name
    plt%description          = plot%description
    plt%x_axis_type          = plot%x_axis_type
    plt%n_curve_pts          = plot%n_curve_pts
    plt%autoscale_gang_x     = plot%autoscale_gang_x 
    plt%autoscale_gang_y     = plot%autoscale_gang_y 
    plt%autoscale_x          = plot%autoscale_x 
    plt%autoscale_y          = plot%autoscale_y 

    call transfer_this_axis(plot%x, default_graph%x) ! Remember: plot%x is deprecated.

    if (default_graph%x%major_div < 0 .and. default_graph%x%major_div_nominal < 0) default_graph%x%major_div_nominal = 7

    call qp_calc_axis_places (default_graph%x)

    do
      ix = index(plt%name, '.')
      if (ix == 0) exit
      call out_io (s_error$, r_name, 'PLOT NAME HAS ".": ' // plt%name, &
                   'SUBSTITUTING "-"', 'IN FILE: ' // plot_file)
      plt%name(ix:ix) = '-'
    enddo

    ng = plot%n_graph
    if (allocated(plt%graph)) deallocate (plt%graph)
    if (ng /= 0) allocate (plt%graph(ng))

    do i_graph = 1, ng
      graph_index = 0         ! setup defaults
      graph = default_graph
      write (graph%name, '(a, i0)') 'g', i_graph
      curve(:) = default_curve
      do j = 1, size(curve)
        write (curve(j)%name, '(a, i0)') 'c', j
      enddo
      if (plt%x_axis_type == 's' .or. plt%x_axis_type == 'lat' .or. &
                  plt%x_axis_type == 'var' .or. plt%x_axis_type == 'curve') then
        curve(:)%draw_symbols = .false.
      else
        curve(:)%draw_symbols = .true.
      endif

      read (iu, nml = tao_template_graph, iostat = ios)
      if (ios > 0) then
        call out_io (s_error$, r_name, &
               'TAO_TEMPLATE_GRAPH NAMELIST READ ERROR.', 'IN FILE: ' // full_file_name, &
               'AFTER TAO_TEMPLATE_PLOT WITH PLOT%NAME = ' // quote(plot%name), &
               'READING TAO_TEMPLATE_GRAPH WITH INDEX: ' // int_str(i_graph))
        rewind (iu)
        do
          read (iu, nml = tao_template_graph)  ! force printing of error message
        enddo
      endif

      call out_io (s_blank$, r_name, 'Init: Read tao_template_graph ' // graph%name)
      graph_name = trim(plot%name) // '.' // graph%name
      call out_io (s_blank$, r_name, 'Init: Read tao_template_graph namelist: ' // graph_name)
      if (graph_index /= i_graph) then
        call out_io (s_error$, r_name, &
              'BAD "GRAPH_INDEX" FOR PLOT: ' // quote(plot%name), &
              'LOOKING FOR GRAPH_INDEX: \I0\ ', &
              'BUT TAO_TEMPLATE_GRAPH HAD GRAPH_INDEX: \I0\ ', 'IN FILE: ' // plot_file, &
              i_array = [i_graph, graph_index] )
        cycle
      endif
      grph   => plt%graph(i_graph)
      grph%p => plt
      grph%name                             = graph%name
      grph%type                             = graph%type
      grph%x_axis_scale_factor              = graph%x_axis_scale_factor 
      grph%symbol_size_scale                = graph%symbol_size_scale   
      grph%curve_legend_row_gap_scale       = graph%curve_legend_row_gap_scale
      grph%text_legend_origin               = graph%text_legend_origin
      grph%curve_legend_origin              = graph%curve_legend_origin
      grph%box                              = graph%box
      grph%title                            = graph%title
      grph%margin                           = graph%margin
      grph%scale_margin                     = graph%scale_margin
      grph%x                                = graph%x
      grph%y                                = graph%y
      grph%x2                               = graph%x2
      grph%y2                               = graph%y2
      grph%ix_universe                      = graph%ix_universe
      grph%ix_branch                        = graph%ix_branch
      grph%clip                             = graph%clip
      grph%allow_wrap_around                = graph%allow_wrap_around
      grph%draw_title                       = graph%draw_title
      grph%draw_axes                        = graph%draw_axes
      grph%draw_grid                        = graph%draw_grid
      grph%draw_only_good_user_data_or_vars = graph%draw_only_good_user_data_or_vars
      grph%draw_curve_legend                = graph%draw_curve_legend
      grph%floor_plan                       = graph%floor_plan      
      grph%title_suffix                     = ''
      grph%text_legend                      = ''
      grph%y2_mirrors_y                     = .true.

      if (grph%x%major_div < 0 .and. grph%x%major_div_nominal < 0) grph%x%major_div_nominal = 7
      if (grph%y%major_div < 0 .and. grph%y%major_div_nominal < 0) grph%y%major_div_nominal = 4
      if (grph%y2%major_div < 0 .and. grph%y2%major_div_nominal < 0) grph%y2%major_div_nominal = 4
      if (graph%floor_plan_orbit_color /= '')           grph%floor_plan%orbit_color = graph%floor_plan_orbit_color ! Old style
      if (graph%floor_plan_orbit_scale /= -1)           grph%floor_plan%orbit_scale = graph%floor_plan_orbit_scale ! Old style
      if (graph%floor_plan_view /=  '')                 grph%floor_plan%view = graph%floor_plan_view
      if (graph%floor_plan_rotation /=  real_garbage$)  grph%floor_plan%rotation = graph%floor_plan_rotation
      if (graph%floor_plan_flip_label_side)             grph%floor_plan%flip_label_side      = graph%floor_plan_flip_label_side
      if (graph%floor_plan_size_is_absolute)            grph%floor_plan%size_is_absolute     = graph%floor_plan_size_is_absolute
      if (graph%floor_plan_draw_only_first_pass)        grph%floor_plan%draw_only_first_pass = graph%floor_plan_draw_only_first_pass
      if (.not. graph%correct_xy_distortion)            grph%floor_plan%correct_distortion   = graph%correct_xy_distortion

      call qp_calc_axis_places (grph%x)

      do
        ix = index(grph%name, '.')
        if (ix == 0) exit
        call out_io (s_error$, r_name, 'GRAPH NAME HAS ".": ' // grph%name, &
                                       'SUBSTITUTING "-"', 'IN FILE: ' // plot_file)
        grph%name(ix:ix) = '-'
      enddo

      call qp_calc_axis_places (grph%y)

      if (grph%ix_universe == 0) then
        call out_io (s_error$, r_name, &
            '**********************************************************', &
            '***** SYNTAX CHANGE: GRAPH%IX_UNIVERSE = 0           *****', &
            '***** NEEDS TO BE CHANGED TO: GRAPH%IX_UNIVERSE = -1 *****', &
            '**********************************************************')
        grph%ix_universe = -1
        stop
      endif

      if (grph%ix_universe < -2 .or. grph%ix_universe > ubound(s%u, 1)) then
        call out_io (s_error$, r_name, 'UNIVERSE INDEX: \i4\ ', & 
                                       'OUT OF RANGE FOR PLOT:GRAPH: ' // graph_name, 'IN FILE: ' // plot_file, &
                                       i_array = [grph%ix_universe] )
        cycle
      endif

      if (grph%type == 'lat_layout') then
        if (plt%x_axis_type /= 's') call out_io (s_error$, r_name, &
                              'A LAT_LAYOUT MUST HAVE X_AXIS_TYPE = "s" FOR A VISIBLE PLOT!', 'IN FILE: ' // plot_file)
        plt%autoscale_gang_y = .false.  ! True does not make sense.
      endif

      ! Count curves. Any curve must have %data_type set or, for %data_source = "data",
      ! %data_type can be constructed from "graph%name.curve%name".

      if (graph%n_curve > 0) then   ! Old style before n_curve was eliminated.
        n_curve = graph%n_curve
      else
        n_curve = 0
        do j = 1, size(curve)
          if (curve(j)%data_type == '' .and. (curve(j)%data_source /= 'data' .or. curve(j)%name == '')) cycle
          n_curve = n_curve + 1
        enddo
      endif

      if (n_curve == 0) then
        if (allocated(grph%curve)) deallocate (grph%curve)
      else
        allocate (grph%curve(n_curve))
      endif

      nc = 0
      do j = 1, size(curve)
        if ((graph%n_curve < 1 .or. j > graph%n_curve) .and. &
              (curve(j)%data_type == '' .and. (curve(j)%data_source /= 'data' .or. curve(j)%name == ''))) cycle
        nc = nc + 1
        if (nc > n_curve) exit
        crv => grph%curve(nc)

        crv%g                    => grph
        crv%data_source          = curve(j)%data_source
        crv%component            = curve(j)%component
        crv%data_index           = curve(j)%data_index
        crv%data_type_x          = curve(j)%data_type_x
        crv%data_type            = curve(j)%data_type
        crv%y_axis_scale_factor  = curve(j)%y_axis_scale_factor
        crv%symbol_every         = curve(j)%symbol_every
        crv%ix_universe          = curve(j)%ix_universe
        crv%n_turn               = curve(j)%n_turn
        crv%draw_line            = curve(j)%draw_line
        crv%draw_symbols         = curve(j)%draw_symbols
        crv%draw_symbol_index    = curve(j)%draw_symbol_index
        crv%draw_error_bars      = curve(j)%draw_error_bars
        crv%use_y2               = curve(j)%use_y2
        crv%symbol               = curve(j)%symbol
        crv%line                 = curve(j)%line
        crv%smooth_line_calc     = curve(j)%smooth_line_calc
        crv%name                 = curve(j)%name
        crv%ele_ref_name         = curve(j)%ele_ref_name
        call str_upcase (crv%ele_ref_name, crv%ele_ref_name)
        crv%ix_bunch             = curve(j)%ix_bunch
        crv%ix_branch            = curve(j)%ix_branch
        crv%legend_text          = curve(j)%legend_text
        crv%hist                 = curve(j)%hist
        crv%z_color              = curve(j)%z_color

        k = min(j, size(symbol_type))
        if (crv%symbol%type == '')  crv%symbol%type  = symbol_type(k)
        if (crv%symbol%color == '') crv%symbol%color = draw_color(k)
        if (crv%line%color == '')   crv%line%color   = draw_color(k)

        if (curve(j)%z_color0 /= invalid$) crv%z_color%min = curve(j)%z_color0
        if (curve(j)%z_color1 /= invalid$) crv%z_color%max = curve(j)%z_color1
        if (curve(j)%autoscale_z_color)    crv%z_color%autoscale = .true.
        if (curve(j)%use_z_color)          crv%z_color%is_on = .true.
        if (curve(j)%data_type_z /= '')    crv%z_color%data_type = curve(j)%data_type_z

        if (crv%data_source == '') crv%data_source = 'lat'

        if (is_integer(crv%symbol%type, ix))         crv%symbol%type          = qp_enum_to_string(ix, 'symbol_type')
        if (is_integer(crv%symbol%color, ix))        crv%symbol%color         = qp_enum_to_string(ix, 'color')
        if (is_integer(crv%symbol%fill_pattern, ix)) crv%symbol%fill_pattern  = qp_enum_to_string(ix, 'fill_pattern')

        if (is_integer(crv%line%color, ix))     crv%line%color     = qp_enum_to_string(ix, 'color')
        if (is_integer(crv%line%pattern, ix))   crv%line%pattern   = qp_enum_to_string(ix, 'line_pattern')

        ! Convert old syntax to new

        if (crv%data_source == 'beam_tracking') crv%data_source = 'beam'
        if (crv%data_source == 'lattice')       crv%data_source = 'lat'
        if (crv%data_source == 'data_array')    crv%data_source = 'data'
        if (crv%data_source == 'var_array')     crv%data_source = 'var'

        ! A floor_plan does not have any associated curves.

        if (grph%type == 'floor_plan') then
          call out_io (s_warn$, r_name, 'FLOOR_PLAN GRAPHS DO NOT HAVE ASSOCIATED CURVES.')
          crv%why_invalid = 'FLOOR_PLAN GRAPHS DO NOT HAVE ASSOCIATED CURVES.'
        endif

        ! multi_turn

        if (crv%n_turn < 0 .and. &
              (crv%data_source == 'multi_turn_orbit' .or. crv%data_source == 'rel_multi_turn_orbit')) then
          call out_io (s_warn$, r_name, 'Curves with %data_source set to "multi_turn_orbit" must now explicity set', &
                                        'the curve''s %n_turn component. In plot: ' // plot%name)
        endif

        if (crv%n_turn > 0 .and. &
              .not. (crv%data_source == 'multi_turn_orbit' .or. crv%data_source == 'rel_multi_turn_orbit')) then
          call out_io (s_warn$, r_name, 'Curve has %n_turn component set positive but %data_source is not set to', &
                                        '"multi_turn_orbit" or "rel_multi_turn_orbit". Multi-turn tracking will', &
                                        '*not* be done for this curve in plot: ' // plot%name)
        endif

        ! Default data type

        if (crv%data_type == '') then
          crv%data_type = trim(plt%name) // '.' // trim(grph%name)
          if (grph%type /= 'dynamic_aperture' .and. size(grph%curve) > 1) call out_io (s_warn$, r_name, &
                        'CURVE(' // int_str(j) // ') OF GRAPH ' // trim(crv%data_type) // ' HAS %data_type NOT SET.', &
                        ' WILL DEFAULT TO: ' // crv%data_type)
        endif

        ! A dot in the name is verboten.
        do
          ix = index(crv%name, '.')
          if (ix == 0) exit
          call out_io (s_error$, r_name, 'CURVE NAME HAS ".": ' // crv%name, 'SUBSTITUTING "-"', 'IN FILE: ' // plot_file)
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
        endif

        if (grph%y2%label_color == null_name$) then
          if (crv%use_y2) then
            grph%y2%label_color = crv%symbol%color
          else
            grph%y2%label_color = 'blue'
          endif
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

        if (crv%ix_universe == 0) then
          call out_io (s_error$, r_name, &
            '**********************************************************', &
            '***** SYNTAX CHANGE: CURVE%IX_UNIVERSE = 0           *****', &
            '***** NEEDS TO BE CHANGED TO: CURVE%IX_UNIVERSE = -1 *****', &
            '**********************************************************')
          crv%ix_universe = -1
          stop
        endif

        i_uni = tao_universe_index (crv%ix_universe)
        if (i_uni > ubound(s%u, 1)) then
          call out_io (s_warn$, r_name, &
                          'CURVE OF PLOT: ' // plot%name, &
                          'HAS UNIVERSE INDEX OUT OF RANGE: \I0\ ', 'IN FILE: ' // plot_file, &
                          i_array = [i_uni] )
          cycle
        endif

        ! Find the ele_ref info.
        ! If plotting something like the phase then the default is for ele_ref 
        ! to be the beginning element.

        if (crv%ele_ref_name /= '') then
          call tao_locate_elements (crv%ele_ref_name, i_uni, eles, err, model$, .true., s_warn$) ! find the index
          if (err) then
            call out_io (s_warn$, r_name, 'CANNOT LOCATE ELEMENT FOR PLOT CURVE: ' // tao_curve_name(crv))
            cycle  ! Check
          endif
          crv%ix_ele_ref = eles(1)%ele%ix_ele
          crv%ix_branch  = eles(1)%ele%ix_branch
        elseif (substr(crv%data_type,1,5) == 'phase' .or. substr(crv%data_type,1,2) == 'r.' .or. &
                substr(crv%data_type,1,2) == 't.' .or. substr(crv%data_type,1,3) == 'tt.') then
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

      ! 

      if (allocated(grph%curve)) then
        do ic = 1, size(grph%curve)
          if (grph%curve(ic)%component == '') grph%curve(ic)%component = graph%component  ! Deprecated syntax.
          if (grph%curve(ic)%component == '') grph%curve(ic)%component = 'model'
        enddo
      endif
    enddo  ! graph
  enddo  ! plot

  close(iu)

enddo  ! file

close (iu)

! If no plots have been defined or default plots wanted then use default

if (ip == 0 .or. include_default_plots) then
  if (size(s%plot_page%region) == 0) deallocate (s%plot_page%region)
  call tao_setup_default_plotting(include_dflt_lat_layout, include_dflt_floor_plan, master_default_graph)
endif

! Initial placement of plots

if (s%global%plot_on .or. s%global%external_plotting) then
  do i = 1, size(place)
    if (place(i)%region == '') cycle
    call tao_place_cmd (place(i)%region, place(i)%plot)
  enddo
endif

! Hook

call number_template_plots()
  if (associated(tao_hook_init_plotting_ptr)) call tao_hook_init_plotting_ptr()

! And finish

call tao_create_plot_window

!----------------------------------------------------------------------------------------
contains

subroutine tao_transfer_shape (shape_input, shape_array, namelist_name)

type (tao_ele_shape_input), target :: shape_input(:)
type (tao_ele_shape_struct), allocatable :: shape_array(:)
type (tao_ele_shape_struct), pointer :: es
integer n, n_max
character(*) namelist_name
character(40) shape, prefix
logical err

!

do n_max = size(shape_input), 1, -1
  if (shape_input(n_max)%ele_id /= '') exit
enddo

if (allocated(shape_array)) deallocate (shape_array)
allocate (shape_array(n_max))

do n = 1, n_max
  shape_array(n) = tao_ele_shape_input_to_struct(shape_input(n), namelist_name)
  call tao_shape_init (shape_array(n), err, .true.)
enddo

end subroutine tao_transfer_shape

!----------------------------------------------------------------------------------------
! contains

subroutine number_template_plots()
integer i
do i = 1, size(s%plot_page%template)
  s%plot_page%template(i)%ix_plot = i
enddo
end subroutine number_template_plots

!----------------------------------------------------------------------------------------
! contains

subroutine tao_setup_default_plotting(include_dflt_lat_layout, include_dflt_floor_plan, default_graph)

type (tao_plot_struct), pointer :: plt
type (tao_graph_struct), pointer :: grph
type (tao_graph_input) default_graph
type (tao_curve_struct), pointer :: crv
type (branch_struct), pointer :: branch
type (tao_plot_struct), target :: default_plot_g1c1, default_plot_g1c2, default_plot_g2c1, default_plot_g1c3, default_plot_g1c4
type (tao_plot_struct), allocatable :: temp_template(:)
type (tao_plot_region_struct), allocatable :: temp_region(:)
type (tao_plot_struct), pointer :: plot
type (tao_ele_shape_struct), allocatable :: temp_shape(:)
type (tao_ele_shape_struct) :: dflt_shapes(30) = [&
      tao_ele_shape_struct('fork::*',              'circle',      'red',     0.15_rp, 'name', .true.,   .false., 1, 0.0_rp, fork$, '*', null()), &
      tao_ele_shape_struct('crystal::*',           'circle',      'red',     0.15_rp, 'name', .true.,   .false., 1, 0.0_rp, crystal$, '*', null()), &
      tao_ele_shape_struct('detector::*',          'box',         'black',   0.30_rp, 'name', .true.,   .false., 1, 0.0_rp, detector$, '*', null()), &
      tao_ele_shape_struct('diffraction_plate::*', 'box',         'cyan',    0.30_rp, 'name', .true.,   .false., 1, 0.0_rp, diffraction_plate$, '*', null()), &
      tao_ele_shape_struct('e_gun::*',             'xbox',        'red',     0.40_rp, 'name', .true.,   .false., 1, 0.0_rp, e_gun$, '*', null()), &
      tao_ele_shape_struct('em_field::*',          'xbox',        'blue',    0.40_rp, 'name', .true.,   .false., 1, 0.0_rp, em_field$, '*', null()), &
      tao_ele_shape_struct('ecollimator::*',       'xbox',        'blue',    0.20_rp, 'name', .false.,  .false., 1, 0.0_rp, ecollimator$, '*', null()), &
      tao_ele_shape_struct('instrument::*',        'box',         'blue',    0.30_rp, 'name', .false.,  .false., 1, 0.0_rp, instrument$, '*', null()), &
      tao_ele_shape_struct('kicker::*',            'u_triangle',  'red',     0.40_rp, 'name', .true.,   .false., 1, 0.0_rp, kicker$, '*', null()), &
      tao_ele_shape_struct('hkicker::*',           'd_triangle',  'red',     0.40_rp, 'name', .true.,   .false., 1, 0.0_rp, hkicker$, '*', null()), &
      tao_ele_shape_struct('vkicker::*',           'u_triangle',  'yellow',  0.40_rp, 'name', .true.,   .false., 1, 0.0_rp, vkicker$, '*', null()), &
      tao_ele_shape_struct('lcavity::*',           'xbox',        'red',     0.50_rp, 'none', .true.,   .false., 1, 0.0_rp, lcavity$, '*', null()), &
      tao_ele_shape_struct('marker::*',            'box',         'blue',    0.30_rp, 'name', .false.,  .false., 1, 0.0_rp, marker$, '*', null()), &
      tao_ele_shape_struct('mirror::*',            'circle',      'red',     0.15_rp, 'name', .true.,   .false., 1, 0.0_rp, mirror$, '*', null()), &
      tao_ele_shape_struct('monitor::*',           'box',         'black',   0.30_rp, 'name', .false.,  .false., 1, 0.0_rp, monitor$, '*', null()), &
      tao_ele_shape_struct('multilayer_mirror::*', 'circle',      'red',     0.15_rp, 'name', .true.,   .false., 1, 0.0_rp, multilayer_mirror$, '*', null()), &
      tao_ele_shape_struct('octupole::*',          'box',         'black',   0.30_rp, 'name', .false.,  .false., 1, 0.0_rp, octupole$, '*', null()), &
      tao_ele_shape_struct('patch::*',             'box',         'yellow',  0.25_rp, 'none', .false.,  .false., 1, 0.0_rp, patch$, '*', null()), &
      tao_ele_shape_struct('photon_fork::*',       'circle',      'red',     0.15_rp, 'name', .true.,   .false., 1, 0.0_rp, photon_fork$, '*', null()), &
      tao_ele_shape_struct('quadrupole::*',        'xbox',        'magenta', 0.37_rp, 'name', .true.,   .false., 1, 0.0_rp, quadrupole$, '*', null()), &
      tao_ele_shape_struct('rcollimator::*',       'xbox',        'blue',    0.20_rp, 'name', .false.,  .false., 1, 0.0_rp, rcollimator$, '*', null()), &
      tao_ele_shape_struct('rfcavity::*',          'xbox',        'red',     0.50_rp, 'name', .true.,   .false., 1, 0.0_rp, rfcavity$, '*', null()), &
      tao_ele_shape_struct('sample::*',            'box',         'black',   0.30_rp, 'name', .true.,   .false., 1, 0.0_rp, sample$, '*', null()), &
      tao_ele_shape_struct('sbend::*',             'box',         'black',   0.20_rp, 'none', .true.,   .false., 1, 0.0_rp, sbend$, '*', null()), &
      tao_ele_shape_struct('sextupole::*',         'xbox',        'green',   0.37_rp, 'none', .true.,   .false., 1, 0.0_rp, sextupole$, '*', null()), &
      tao_ele_shape_struct('sol_quad::*',          'box',         'black',   0.40_rp, 'name', .false.,  .false., 1, 0.0_rp, sol_quad$, '*', null()), &
      tao_ele_shape_struct('solenoid::*',          'box',         'blue',    0.30_rp, 'name', .true.,   .false., 1, 0.0_rp, solenoid$, '*', null()), &
      tao_ele_shape_struct('wiggler::*',           'xbox',        'cyan',    0.50_rp, 'name', .true.,   .false., 1, 0.0_rp, wiggler$, '*', null()), &
      tao_ele_shape_struct('photon_init::*',       'box',         'black',   0.30_rp, 'name', .true.,   .false., 1, 0.0_rp, photon_init$, '*', null()), &
      tao_ele_shape_struct('building_wall::*',     'solid_line',  'black',   0.30_rp, 'name', .true.,   .false., 3, 0.0_rp, 999, '*', null())]

real(rp) y_layout, dx, dy, dz, x1, x2, y1, y2
integer np, n, nr, n_plots
integer i, j, k, ie, ic, n_old
character(20) :: color(19)
logical include_dflt_lat_layout, include_dflt_floor_plan
character(40) name

!

call tao_set_plotting (plot_page, s%plot_page, .true.)

n = size(dflt_shapes)

if (include_dflt_lat_layout .or. size(s%plot_page%lat_layout%ele_shape) == 0) then
  n_old = 0
  if (allocated(s%plot_page%lat_layout%ele_shape)) then
    n_old = size(s%plot_page%lat_layout%ele_shape)
    call move_alloc (s%plot_page%lat_layout%ele_shape, temp_shape)
  endif
  allocate (s%plot_page%lat_layout%ele_shape(40+n_old))
  s%plot_page%lat_layout%ele_shape(:)%ele_id = ''
  if (n_old /= 0) s%plot_page%lat_layout%ele_shape(1:n_old) = temp_shape
  s%plot_page%lat_layout%ele_shape(n_old+1:n_old+n) = dflt_shapes
  if (n_old /= 0) deallocate(temp_shape)
endif

if (include_dflt_floor_plan .or. size(s%plot_page%floor_plan%ele_shape) == 0) then
  if (allocated(s%plot_page%floor_plan%ele_shape)) then
    n_old = size(s%plot_page%floor_plan%ele_shape)
    call move_alloc (s%plot_page%floor_plan%ele_shape, temp_shape)
  endif
  allocate (s%plot_page%floor_plan%ele_shape(40+n_old))
  s%plot_page%floor_plan%ele_shape(:)%ele_id = ''
  if (n_old /= 0) s%plot_page%floor_plan%ele_shape(1:n_old) = temp_shape
  s%plot_page%floor_plan%ele_shape(n_old+1:n_old+n) = dflt_shapes
  s%plot_page%floor_plan%ele_shape(n_old+1:n_old+n)%size = 20 * s%plot_page%floor_plan%ele_shape(n_old+1:n_old+n)%size
  if (n_old /= 0) deallocate(temp_shape)
endif

!---------------------------------

n_plots = 100

if (allocated(s%plot_page%template)) then
  n = size(s%plot_page%template)
  call move_alloc(s%plot_page%template, temp_template)
  allocate (s%plot_page%template(n + n_plots))
  s%plot_page%template(1:n) = temp_template
  deallocate (temp_template)
  np = n
  if (s%plot_page%template(np)%name == 'scratch') np = np - 1

  do j = 1, n
    plot => s%plot_page%template(j)
    if (.not. allocated(plot%graph)) cycle
    do k = 1, size(plot%graph)
      plot%graph(k)%p => plot
      if (.not. allocated(plot%graph(k)%curve)) cycle
      do ic = 1, size(plot%graph(k)%curve)
        plot%graph(k)%curve(ic)%g => plot%graph(k)
      enddo
    enddo
  enddo

else
  allocate (s%plot_page%template(n_plots))
  np = 0
endif

!---------------
! This plot defines the default 2-graph, 1-curve/graph plot

plt => default_plot_g2c1

nullify(plt%r)
if (allocated(plt%graph)) deallocate (plt%graph)
allocate (plt%graph(2))
allocate (plt%graph(1)%curve(1))
allocate (plt%graph(2)%curve(1))

plt%x_axis_type          = 's'

grph => plt%graph(1)
grph%name                 = 'g1'
grph%type                 = 'data'
grph%margin               = qp_rect_struct(0.15, 0.06, 0.12, 0.12, '%BOX')
grph%scale_margin         = qp_rect_struct(0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, '%GRAPH')
grph%box                  = [1, 1, 1, 1]
grph%y                    = init_axis
grph%y%label_offset       = 0.15
grph%y%major_div_nominal  = 4
grph%y2                   = init_axis
grph%y2%major_div_nominal = 4
grph%y2%draw_numbers      = .false.
grph%draw_axes            = .true.
grph%draw_grid            = .true.
grph%x                    = init_axis
grph%x%major_div_nominal  = 7
grph%x%minor_div_max      = 6
grph%x%label = 's [m]'

crv => grph%curve(1)
crv%name         = 'c'
crv%component    = 'model'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = 'blue'
crv%line%width   = 2
crv%symbol%color = crv%line%color

grph => plt%graph(2)
grph%name                 = 'g2'
grph%type                 = 'data'
grph%margin               = qp_rect_struct(0.15, 0.06, 0.12, 0.12, '%BOX')
grph%scale_margin         = qp_rect_struct(0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, '%GRAPH')
grph%box                  = [1, 2, 1, 2]
grph%y                    = init_axis
grph%y%label_offset       = 0.15
grph%y%major_div_nominal  = 4
grph%y2                   = init_axis
grph%y2%major_div_nominal = 4
grph%y2%draw_numbers      = .false.
grph%draw_axes            = .true.
grph%draw_grid            = .true.
grph%x                    = init_axis

crv => grph%curve(1)
crv%name         = 'c'
crv%component    = 'model'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = 'blue'
crv%line%width   = 2
crv%symbol%color = crv%line%color

!---------------
! This plot defines the default 1-graph, 3-curve/graph plot

plt => default_plot_g1c3

nullify(plt%r)
if (allocated(plt%graph)) deallocate (plt%graph)
allocate (plt%graph(1))
allocate (plt%graph(1)%curve(3))

plt%x_axis_type          = 's'

grph => plt%graph(1)
grph%name                 = 'g'
grph%type                 = 'data'
grph%margin               = qp_rect_struct(0.15, 0.06, 0.12, 0.12, '%BOX')
grph%scale_margin         = qp_rect_struct(0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, '%GRAPH')
grph%box                  = [1, 1, 1, 1]
grph%y                    = init_axis
grph%y%label_offset       = 0.15
grph%y%major_div_nominal  = 4
grph%y2                   = init_axis
grph%y2%major_div_nominal = 4
grph%y2%draw_numbers      = .false.
grph%draw_curve_legend    = .true.
grph%draw_axes            = .true.
grph%draw_grid            = .true.
grph%text_legend_origin   = default_graph%text_legend_origin
grph%curve_legend_origin  = default_graph%curve_legend_origin
grph%x                    = init_axis
grph%x%major_div_nominal  = 7
grph%x%minor_div_max      = 6
grph%x%label              = 's [m]'

crv => grph%curve(1)
crv%name         = 'c1'
crv%component    = 'model'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = 'blue'
crv%line%width   = 2
crv%symbol%color = crv%line%color

crv => grph%curve(2)
crv%name         = 'c2'
crv%component    = 'model'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = 'orange'
crv%line%width   = 2
crv%symbol%color = crv%line%color

crv => grph%curve(3)
crv%name         = 'c3'
crv%component    = 'model'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = 'green'
crv%line%width   = 2
crv%symbol%color = crv%line%color

!---------------
! This plot defines the default 1-graph, 4-curve/graph plot

plt => default_plot_g1c4

nullify(plt%r)
if (allocated(plt%graph)) deallocate (plt%graph)
allocate (plt%graph(1))
allocate (plt%graph(1)%curve(4))

plt%x_axis_type          = 's'

grph => plt%graph(1)
grph%name                 = 'g'
grph%type                 = 'data'
grph%margin               = qp_rect_struct(0.15, 0.06, 0.12, 0.12, '%BOX')
grph%scale_margin         = qp_rect_struct(0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, '%GRAPH')
grph%box                  = [1, 1, 1, 1]
grph%y                    = init_axis
grph%y%label_offset       = 0.15
grph%y%major_div_nominal  = 4
grph%y2                   = init_axis
grph%y2%major_div_nominal = 4
grph%y2%draw_numbers      = .false.
grph%draw_curve_legend    = .true.
grph%draw_axes            = .true.
grph%draw_grid            = .true.
grph%text_legend_origin   = default_graph%text_legend_origin
grph%curve_legend_origin  = default_graph%curve_legend_origin
grph%x                    = init_axis
grph%x%major_div_nominal  = 7
grph%x%minor_div_max      = 6
grph%x%label              = 's [m]'

crv => grph%curve(1)
crv%name         = 'c1'
crv%component    = 'model'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = 'blue'
crv%line%width   = 2
crv%symbol%color = crv%line%color

crv => grph%curve(2)
crv%name         = 'c2'
crv%component    = 'model'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = 'orange'
crv%line%width   = 2
crv%symbol%color = crv%line%color

crv => grph%curve(3)
crv%name         = 'c3'
crv%component    = 'model'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = 'green'
crv%line%width   = 2
crv%symbol%color = crv%line%color

crv => grph%curve(4)
crv%name         = 'c4'
crv%component    = 'model'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = 'magenta'
crv%line%width   = 2
crv%symbol%color = crv%line%color

!---------------
! This plot defines the default 1-graph, 2-curve/graph plot

plt => default_plot_g1c2

nullify(plt%r)
if (allocated(plt%graph)) deallocate (plt%graph)
allocate (plt%graph(1))
allocate (plt%graph(1)%curve(2))

plt%x_axis_type          = 's'

grph => plt%graph(1)
grph%name                 = 'g'
grph%type                 = 'data'
grph%margin               = qp_rect_struct(0.15, 0.06, 0.12, 0.12, '%BOX')
grph%scale_margin         = qp_rect_struct(0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, '%GRAPH')
grph%box                  = [1, 1, 1, 1]
grph%y                    = init_axis
grph%y%label_offset       = 0.15
grph%y%major_div_nominal  = 4
grph%y2                   = init_axis
grph%y2%major_div_nominal = 4
grph%y2%draw_numbers      = .false.
grph%draw_curve_legend    = .true.
grph%draw_axes            = .true.
grph%draw_grid            = .true.
grph%text_legend_origin   = default_graph%text_legend_origin
grph%curve_legend_origin  = default_graph%curve_legend_origin
grph%x                    = init_axis
grph%x%major_div_nominal  = 7
grph%x%minor_div_max      = 6
grph%x%label              = 's [m]'

crv => grph%curve(1)
crv%name         = 'c1'
crv%component    = 'model'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = 'blue'
crv%line%width   = 2
crv%symbol%color = crv%line%color

crv => grph%curve(2)
crv%name         = 'c2'
crv%component    = 'model'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = 'orange'
crv%line%width   = 2
crv%symbol%color = crv%line%color

!---------------
! This plot defines the default 1-graph, 1-curve/graph plot

plt => default_plot_g1c1 

nullify(plt%r)
if (allocated(plt%graph)) deallocate (plt%graph)
allocate (plt%graph(1))
allocate (plt%graph(1)%curve(1))

plt%x_axis_type          = 's'

grph => plt%graph(1)
grph%name                 = 'g'
grph%type                 = 'data'
grph%margin               = qp_rect_struct(0.15, 0.06, 0.12, 0.12, '%BOX')
grph%scale_margin         = qp_rect_struct(0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, '%GRAPH')
grph%box                  = [1, 1, 1, 1]
grph%y                    = init_axis
grph%y%label_offset       = 0.15
grph%y%major_div_nominal  = 4
grph%y2                   = init_axis
grph%y2%major_div_nominal = 4
grph%y2%draw_numbers      = .false.
grph%draw_axes            = .true.
grph%draw_grid            = .true.
grph%x                    = init_axis
grph%x%major_div_nominal  = 7
grph%x%minor_div_max      = 6
grph%x%label              = 's [m]'

crv => grph%curve(1)
crv%name         = 'c'
crv%component    = 'model'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = 'blue'
crv%line%width   = 2
crv%symbol%color = crv%line%color

!---------------
! alpha plot

if (all(s%plot_page%template%name /= 'alpha')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'alpha'
  plt%description          = 'Twiss alpha function'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Alpha Function'
  grph%y%label             = '\ga\fn\dA\u, \ga\fn\dB\u [m]'

  crv => grph%curve(1)
  crv%name         = 'a'
  crv%g => grph
  crv%data_type    = 'alpha.a'
  crv%legend_text  = '\ga\fn\dA\u'

  crv => grph%curve(2)
  crv%name         = 'b'
  crv%g => grph
  crv%data_type    = 'alpha.b'
  crv%legend_text  = '\ga\fn\dB\u'
endif

!---------------
! b_div_curl plot

if (all(s%plot_page%template%name /= 'b_div_curl')) then
  call default_plot_init (np, plt, default_plot_g1c4)
  plt%name                 = 'b_div_curl'
  plt%description          = 'Mag Field Divergence and (Curl B - (dE/dt)/c\u2\d) Along Orbit'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'B-Field Div and (Curl B - (dE/dt)/c\u2\d) Along Orbit'
  grph%y%label             = 'Mag Div, Curl (T/m)'
  grph%y%label_offset= .15

  crv => grph%curve(1)
  crv%name         = 'div'
  crv%g => grph
  crv%data_type    = 'b_div'
  crv%legend_text  = 'B Divergence at particle position'

  crv => grph%curve(2)
  crv%name         = 'cx'
  crv%g => grph
  crv%data_type    = 'b_curl.x'
  crv%legend_text  = '(curl B)x - (dEx/dt)/c\u2\d at particle position'

  crv => grph%curve(3)
  crv%name         = 'cy'
  crv%g => grph
  crv%data_type    = 'b_curl.y'
  crv%legend_text  = '(curl B)y - (dEy/dt)/c\u2\d at particle position'

  crv => grph%curve(4)
  crv%name         = 'cz'
  crv%g => grph
  crv%data_type    = 'b_curl.z'
  crv%legend_text  = '(curl B)z - (dEy/dt)/c\u2\d at particle position'
endif
 
!---------------
! B_field plot

if (all(s%plot_page%template%name /= 'b_field')) then
  call default_plot_init (np, plt, default_plot_g1c3)
  plt%name                 = 'b_field'
  plt%description          = 'Magnetic Field Along Orbit'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Magnetic Field Along Orbit'
  grph%y%label             = 'B-Field (T)'
  grph%y%label_offset= .15

  crv => grph%curve(1)
  crv%name         = 'bx'
  crv%g => grph
  crv%data_type    = 'b_field.x'
  crv%legend_text  = 'B_field.x at particle position'

  crv => grph%curve(2)
  crv%name         = 'by'
  crv%g => grph
  crv%data_type    = 'b_field.y'
  crv%legend_text  = 'B_field.y at particle position'

  crv => grph%curve(3)
  crv%name         = 'bz'
  crv%g => grph
  crv%data_type    = 'b_field.z'
  crv%legend_text  = 'B_field.z at particle position'
endif
 
!---------------
! b0_div_curl plot

if (all(s%plot_page%template%name /= 'b0_div_curl')) then
  call default_plot_init (np, plt, default_plot_g1c4)
  plt%name                 = 'b0_div_curl'
  plt%description          = 'B-Field Divergence and (Curl B - (dE/dt)/c\u2\d) at const x,y,t'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'B-Field Div and (Curl B - (dE/dt)/c\u2\d) at const x,y,t'
  grph%y%label             = 'B-Field Div, Curl (T/m)'
  grph%y%label_offset= .15

  crv => grph%curve(1)
  crv%name         = 'div'
  crv%g => grph
  crv%data_type    = 'b0_div'
  crv%legend_text  = 'B Divergence'

  crv => grph%curve(2)
  crv%name         = 'cx'
  crv%g => grph
  crv%data_type    = 'b0_curl.x'
  crv%legend_text  = '(curl B)x - (dEx/dt)/c\u2\d'

  crv => grph%curve(3)
  crv%name         = 'cy'
  crv%g => grph
  crv%data_type    = 'b0_curl.y'
  crv%legend_text  = '(curl B)y - (dEy/dt)/c\u2\d'

  crv => grph%curve(4)
  crv%name         = 'cz'
  crv%g => grph
  crv%data_type    = 'b0_curl.z'
  crv%legend_text  = '(curl B)z - (dEy/dt)/c\u2\d'
endif
 
!---------------
! B0_field plot

if (all(s%plot_page%template%name /= 'b0_field')) then
  call default_plot_init (np, plt, default_plot_g1c3)
  plt%name                 = 'b0_field'
  plt%description          = 'Magnetic Field at const x,y,t'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Magnetic Field at const x,y,t'
  grph%y%label             = 'B-Field (T)'
  grph%y%label_offset= .15

  crv => grph%curve(1)
  crv%name         = 'bx'
  crv%g => grph
  crv%data_type    = 'b0_field.x'
  crv%legend_text  = 'B_field.x'

  crv => grph%curve(2)
  crv%name         = 'by'
  crv%g => grph
  crv%data_type    = 'b0_field.y'
  crv%legend_text  = 'B_field.y'

  crv => grph%curve(3)
  crv%name         = 'bz'
  crv%g => grph
  crv%data_type    = 'b0_field.z'
  crv%legend_text  = 'B_field.z'
endif
 
!---------------
! beta plot

if (all(s%plot_page%template%name /= 'beta')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'beta'
  plt%description          = 'Twiss beta function'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Beta Function'
  grph%y%label             = '\gb\fn\dA\u, \gb\fn\dB\u [m]'


  crv => grph%curve(1)
  crv%name         = 'a'
  crv%g => grph
  crv%data_type    = 'beta.a'
  crv%legend_text  = '\gb\fn\dA\u'

  crv => grph%curve(2)
  crv%name         = 'b'
  crv%g => grph
  crv%data_type    = 'beta.b'
  crv%legend_text  = '\gb\fn\dB\u'
endif

!---------------
! bunch_sigma_xy plot

if (all(s%plot_page%template%name /= 'bunch_sigma_xy')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'bunch_sigma_xy'
  plt%description          = 'Bunch transverse sigmas'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Bunch transverse sigmas'
  grph%y%label             = '\gs\fn\dx\u, \gs\fn\dy\u'

  crv => grph%curve(1)
  crv%name         = 'x'
  crv%g => grph
  crv%data_type    = 'sigma.x'
  crv%legend_text  = 'Sigma_x'
  crv%data_source  = 'beam' 
  crv%smooth_line_calc = .false.

  crv => grph%curve(2)
  crv%name         = 'y'
  crv%g => grph
  crv%data_type    = 'sigma.y'
  crv%legend_text  = 'Sigma_y'
  crv%data_source  = 'beam' 
  crv%smooth_line_calc = .false.
endif

!---------------
! bunch_x_px, etc. plots

np = np + 1
plt => s%plot_page%template(np)
plt%phantom = .true.
plt%name = 'bunch_<R1>_<R2>'
plt%description = 'Bunch phase space plot. <R1>, <R2> -> x, px, y, py, z, or pz. EG: bunch_z_pz'

do i = 1, 6
do j = 1, 6

  if (i == j) cycle
  name = 'bunch_' // trim(coord_name(i)) // '_' // trim(coord_name(j))
  if (any(s%plot_page%template%name == name)) cycle

  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = name
  plt%description          = 'Bunch phase space'
  plt%list_with_show_plot_command = .false.
  plt%x_axis_type = 'phase_space'

  grph => plt%graph(1)
  grph%p => plt
  grph%type                = 'phase_space'
  grph%title               = trim(coord_name_cap(i)) // ' Vs ' // coord_name_cap(j)
  grph%x%label             = coord_name_cap(i)
  grph%y%label             = coord_name_cap(j)
  grph%x%major_div_nominal = 4

  crv => grph%curve(1)
  crv%name         = 'c'
  crv%g => grph
  crv%data_source = 'beam'
  crv%data_type_x  = coord_name(i)
  crv%data_type    = coord_name(j)
  crv%draw_symbols = .true.
  crv%draw_line    = .false.
  ie = s%u(1)%design%lat%n_ele_track
  crv%ix_ele_ref = ie
  crv%ix_ele_ref_track = ie
  crv%ele_ref_name = s%u(1)%design%lat%ele(ie)%name
  if (modulo(i,2) == 0) then
  else
  endif

enddo
enddo

!---------------
! bunch_density_x, etc. plots

np = np + 1
plt => s%plot_page%template(np)
plt%phantom = .true.
plt%name = 'bunch_density_<R>'
plt%description = 'Bunch Charge Density Histogram. <R> -> x, px, y, py, z, or pz. EG: bunch_density_z'

do i = 1, 6
  name = 'bunch_density_' // trim(coord_name(i))
  if (any(s%plot_page%template%name == name)) cycle

  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = name
  plt%description          = 'Bunch density histogram'
  plt%list_with_show_plot_command = .false.
  plt%x_axis_type = 'histogram'

  grph => plt%graph(1)
  grph%p => plt
  graph%name               = coord_name(i)
  grph%type                = 'histogram'
  grph%title               = 'Bunch Density Histogram in ' // trim(coord_name_cap(i)) 
  grph%x%major_div_nominal = 4
  grph%x%label = coord_name_cap(i)
  grph%x_axis_scale_factor = 1d3  ! mm

  crv => grph%curve(1)
  crv%g => grph
  crv%name         = 'c'
  crv%data_source  = 'beam'
  crv%data_type    = coord_name(i)

  crv%hist%density_normalized = .true.
  crv%hist%weight_by_charge = .true.
  crv%hist%number = 100

  crv%draw_symbols = .true.
  crv%symbol%type = 'dot'

  crv%draw_line    = .true.
  crv%line%color = 'blue'
  crv%line%pattern = 'dashed'

  crv%y_axis_scale_factor = 1e9   ! nC

  ie = s%u(1)%design%lat%n_ele_track
  crv%ix_ele_ref = ie
  crv%ix_ele_ref_track = ie
  crv%ele_ref_name = s%u(1)%design%lat%ele(ie)%name

  if (modulo(i,2) == 0) then
    grph%x%label     = trim(coord_name_cap(i)) // ' [* 10^3]'
    grph%y%label     = 'Charge Density [nC/' // trim(graph%x%label) // ']'
  else
    grph%x%label     = trim(coord_name_cap(i)) // ' [mm]'
    grph%y%label     = 'Charge Density [nC/mm]'
  endif

enddo

!---------------
! cbar plot

if (all(s%plot_page%template%name /= 'cbar')) then
  call default_plot_init (np, plt, default_plot_g1c4)
  plt%name                 = 'cbar'
  plt%description          = 'Cbar coupling matrix'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Cbar coupling matrix'
  grph%y%label             = 'Cbar'

  crv => grph%curve(1)
  crv%name         = '11'
  crv%g => grph
  crv%data_type    = 'cbar.11'
  crv%legend_text  = 'Cbar11'

  crv => grph%curve(2)
  crv%name         = '12'
  crv%g => grph
  crv%data_type    = 'cbar.12'
  crv%legend_text  = 'Cbar12'

  crv => grph%curve(3)
  crv%name         = '21'
  crv%g => grph
  crv%data_type    = 'cbar.21'
  crv%legend_text  = 'Cbar21'

  crv => grph%curve(4)
  crv%name         = '22'
  crv%g => grph
  crv%data_type    = 'cbar.22'
  crv%legend_text  = 'Cbar22'
endif

!---------------
! dbeta (chrom.dbeta) plot

if (all(s%plot_page%template%name /= 'dbeta')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'dbeta'
  plt%description          = 'Chromatic normalized beta beat'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Chromatic normalized beta beat'
  grph%y%label             = '\gb\fn\u-1\d \(2265)\gb\fn/\(2265)\gd\fn (1)'
  grph%y%label_offset= .15

  crv => grph%curve(1)
  crv%name         = 'a'
  crv%g => grph
  crv%data_type    = 'chrom.dbeta.a'
  crv%legend_text  = 'a'
  crv%smooth_line_calc = .false.

  crv => grph%curve(2)
  crv%name         = 'b'
  crv%g => grph
  crv%data_type    = 'chrom.dbeta.b'
  crv%legend_text  = 'b'
  crv%smooth_line_calc = .false.
endif

!---------------
! deta (chrom.deta) plot

if (all(s%plot_page%template%name /= 'deta')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'deta'
  plt%description          = 'Second order dispersion'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Second order dispersion'
  grph%y%label             = '\(2265)\gy\fn/\(2265)\gd\fn (m)'
  grph%y%label_offset= .15

  crv => grph%curve(1)
  crv%name         = 'x'
  crv%g => grph
  crv%data_type    = 'chrom.deta.x'
  crv%legend_text  = 'x'
  crv%smooth_line_calc = .false.

  crv => grph%curve(2)
  crv%name         = 'y'
  crv%g => grph
  crv%data_type    = 'chrom.deta.y'
  crv%legend_text  = 'y'
  crv%smooth_line_calc = .false.
endif

!---------------
! detap (chrom.detap) plot

if (all(s%plot_page%template%name /= 'detap')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'detap'
  plt%description          = 'Second order dispersion slope'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Second order dispersion slope'
  grph%y%label             = "\(2265)\gy\fn'/\(2265)\gd\fn (1)"
  grph%y%label_offset= .15

  crv => grph%curve(1)
  crv%name         = 'x'
  crv%g => grph
  crv%data_type    = 'chrom.detap.x'
  crv%legend_text  = 'x'
  crv%smooth_line_calc = .false.

  crv => grph%curve(2)
  crv%name         = 'y'
  crv%g => grph
  crv%data_type    = 'chrom.detap.y'
  crv%legend_text  = 'y'
  crv%smooth_line_calc = .false.
endif

!---------------
! Dispersion plot

if (all(s%plot_page%template%name /= 'dispersion')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name           = 'dispersion'
  plt%description    = 'X & Y Dispersion'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Dispersion'
  grph%y%label       = '\gy\fn\dX\u, \gy\fn\dY\u [m]'

  crv => grph%curve(1)
  crv%name         = 'x'
  crv%g => grph
  crv%data_type    = 'eta.x'
  crv%legend_text  = '\gy\fn\dX\u'

  crv => grph%curve(2)
  crv%name         = 'y'
  crv%g => grph
  crv%data_type = 'eta.y'
  crv%legend_text  = '\gy\fn\dY\u'
endif

!---------------
! Dispersion derivative plot

if (all(s%plot_page%template%name /= 'ddispersion')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name           = 'ddispersion'
  plt%description    = 'X & Y Dispersion Derivative'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Dispersion Derivative'
  grph%y%label       = 'd\gy\fn\dX\u, d\gy\fn\dY\u'

  crv => grph%curve(1)
  crv%name         = 'x'
  crv%g => grph
  crv%data_type    = 'etap.x'
  crv%legend_text  = 'd\gy\fn\dX\u'

  crv => grph%curve(2)
  crv%name         = 'y'
  crv%g => grph
  crv%data_type = 'etap.y'
  crv%legend_text  = 'd\gy\fn\dY\u'
endif

!---------------
! Normal mode Dispersion plot

if (all(s%plot_page%template%name /= 'mode_dispersion')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name           = 'mode_dispersion'
  plt%description    = 'A & B Normal Mode Dispersion'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'A & B Normal Mode Dispersion'
  grph%y%label       = '\gy\fn\dA\u, \gy\fn\dB\u [m]'

  crv => grph%curve(1)
  crv%name         = 'a'
  crv%g => grph
  crv%data_type    = 'eta.a'
  crv%legend_text  = '\gy\fn\dA\u'

  crv => grph%curve(2)
  crv%name         = 'b'
  crv%g => grph
  crv%data_type = 'eta.b'
  crv%legend_text  = '\gy\fn\dB\u'
endif

!---------------
! dphi (chrom.dphi) plot

if (all(s%plot_page%template%name /= 'dphi')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'dphi'
  plt%description          = 'Chromatic phase deviation'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Chromatic phase deviation'
  grph%y%label             = '\(2265)\gf\fn/\(2265)\gd\fn (1)'
  grph%y%label_offset= .15

  crv => grph%curve(1)
  crv%name         = 'a'
  crv%g => grph
  crv%data_type    = 'chrom.dphi.a'
  crv%legend_text  = 'a'
  crv%smooth_line_calc = .false.

  crv => grph%curve(2)
  crv%name         = 'b'
  crv%g => grph
  crv%data_type    = 'chrom.dphi.b'
  crv%legend_text  = 'b'
  crv%smooth_line_calc = .false.
endif
 
!---------------
! dynamic_aperture plot

if (all(s%plot_page%template%name /= 'dynamic_aperture')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  deallocate(plt%graph(1)%curve)
  allocate (plt%graph(1)%curve(19))

  plt%name                 = 'dynamic_aperture'
  plt%description          = 'Dynamic aperture using universe calc'
  plt%x_axis_type = 'phase_space'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'dynamic aperture'
  grph%type                = 'dynamic_aperture'
  grph%x_axis_scale_factor = 1000
  grph%x%label             = 'x (mm)'
  grph%y%label             = 'y (mm)'
  grph%y%label_offset= .15

  crv => grph%curve(1)
  crv%g => grph
  crv%name = 'e'
  crv%line%color = "red"
  crv%data_type = "beam_ellipse"
  crv%line%width = 4
  crv%legend_text = "10 sigma beam ellipse"
  crv%draw_symbols = .false.
  crv%y_axis_scale_factor = 1000

  color(1:9) = [character(20):: 'blue', 'green', 'cyan', 'magenta', 'yellow', 'black', 'orange', 'purple', 'light_grey']
  color(10:18) = color(1:9)

  do i = 2, size(grph%curve)
    crv => grph%curve(i)
    crv%name         = 'c' // int_str(i-1)
    crv%g => grph
    crv%smooth_line_calc = .false.
    crv%y_axis_scale_factor = 1000
    crv%line%width = 4
    crv%draw_symbols = .true.
    crv%data_type = 'dynamic_aperture'
    crv%data_index = int_str(i-1)
    crv%symbol%color = color(i-1)
    crv%line%color = color(i-1)
  enddo
endif

!---------------
! e_div_curl plot

if (all(s%plot_page%template%name /= 'e_div_curl')) then
  call default_plot_init (np, plt, default_plot_g1c4)
  plt%name                 = 'e_div_curl'
  plt%description          = 'Electric Field Divergence and (Curl E - dB/dt) Along Orbit'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'E-Field Div and (Curl E - dB/dt) Along Orbit'
  grph%y%label             = 'Elec Div, Curl (V/m^2)'
  grph%y%label_offset= .15

  crv => grph%curve(1)
  crv%name         = 'div'
  crv%g => grph
  crv%data_type    = 'e_div'
  crv%legend_text  = 'divergence'

  crv => grph%curve(2)
  crv%name         = 'cx'
  crv%g => grph
  crv%data_type    = 'e_curl.x'
  crv%legend_text  = '(curl E)x - dBx/dt at particle position'

  crv => grph%curve(3)
  crv%name         = 'cy'
  crv%g => grph
  crv%data_type    = 'e_curl.y'
  crv%legend_text  = '(curl E)y - dBy/dt at particle position'

  crv => grph%curve(4)
  crv%name         = 'cz'
  crv%g => grph
  crv%data_type    = 'e_curl.z'
  crv%legend_text  = '(curl E)z - dBz/dt at particle position'
endif
 
!---------------
! E_field plot

if (all(s%plot_page%template%name /= 'e_field')) then
  call default_plot_init (np, plt, default_plot_g1c3)
  plt%name                 = 'e_field'
  plt%description          = 'Electric Field Along Orbit'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Electric Field Along Orbit'
  grph%y%label             = 'E-Field (V/m)'
  grph%y%label_offset= .15

  crv => grph%curve(1)
  crv%name         = 'ex'
  crv%g => grph
  crv%data_type    = 'e_field.x'
  crv%legend_text  = 'E_field.x at particle position'

  crv => grph%curve(2)
  crv%name         = 'ey'
  crv%g => grph
  crv%data_type    = 'e_field.y'
  crv%legend_text  = 'E_field.y at particle position'

  crv => grph%curve(3)
  crv%name         = 'ez'
  crv%g => grph
  crv%data_type    = 'e_field.z'
  crv%legend_text  = 'E_field.z at particle position'

endif

!---------------
! e0_div_curl plot

if (all(s%plot_page%template%name /= 'e0_div_curl')) then
  call default_plot_init (np, plt, default_plot_g1c4)
  plt%name                 = 'e0_div_curl'
  plt%description          = 'Electric Field Divergence and (Curl E - dB/dt) at const x,y,t'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'E-Field Div and (Curl E - dB/dt) at const x,y,t'
  grph%y%label             = 'Elec Div, Curl (V/m^2)'
  grph%y%label_offset= .15

  crv => grph%curve(1)
  crv%name         = 'div'
  crv%g => grph
  crv%data_type    = 'e0_div'
  crv%legend_text  = 'divergence'

  crv => grph%curve(2)
  crv%name         = 'cx'
  crv%g => grph
  crv%data_type    = 'e0_curl.x'
  crv%legend_text  = '(curl E)x - dBx/dt'

  crv => grph%curve(3)
  crv%name         = 'cy'
  crv%g => grph
  crv%data_type    = 'e0_curl.y'
  crv%legend_text  = '(curl E)y - dBy/dt'

  crv => grph%curve(4)
  crv%name         = 'cz'
  crv%g => grph
  crv%data_type    = 'e0_curl.z'
  crv%legend_text  = '(curl E)z - dBz/dt'
endif
 
!---------------
! E0_field plot

if (all(s%plot_page%template%name /= 'e0_field')) then
  call default_plot_init (np, plt, default_plot_g1c3)
  plt%name                 = 'e0_field'
  plt%description          = 'Electric Field at const x,y,t'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'E-Field at const x,y,t'
  grph%y%label             = 'E-Field (V/m)'
  grph%y%label_offset= .15

  crv => grph%curve(1)
  crv%name         = 'ex'
  crv%g => grph
  crv%data_type    = 'e0_field.x'
  crv%legend_text  = 'E_field.x'

  crv => grph%curve(2)
  crv%name         = 'ey'
  crv%g => grph
  crv%data_type    = 'e0_field.y'
  crv%legend_text  = 'E_field.y'

  crv => grph%curve(3)
  crv%name         = 'ez'
  crv%g => grph
  crv%data_type    = 'e0_field.z'
  crv%legend_text  = 'E_field.z'

endif

!---------------
! emittance growth

if (all(s%plot_page%template%name /= 'emittance')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'emittance'
  plt%description          = 'Linac emittance'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Linac Emittance Growth [Rad_Integrals: (I5a_E6, I5b_E6) * 2 * C_q * r_e / 3]'
  grph%y%label       = 'Emit Growth [mm-mrad]'

  crv => grph%curve(1)
  crv%name         = 'a'
  crv%g => grph
  crv%data_type    = 'rad_int.i5a_e6'
  crv%legend_text  = 'a-mode emit'
  crv%y_axis_scale_factor = 7.213927194325027E-22 !for mm-mrad
  crv%smooth_line_calc = .false.

  crv => grph%curve(2)
  crv%name         = 'b'
  crv%g => grph
  crv%data_type    = 'rad_int.i5b_e6'
  crv%legend_text  = 'b-mode emit'
  crv%y_axis_scale_factor = 7.213927194325027E-22 !for mm-mrad
  crv%smooth_line_calc = .false.
endif

!---------------
! energy

if (all(s%plot_page%template%name /= 'energy')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'energy'
  plt%description          = 'Particle energy'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Total Energy'
  grph%y%label       = 'E\dTot\u [eV]'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'orbit.energy'
endif

!---------------
! Floor Plan plot

if (all(s%plot_page%template%name /= 'floor_plan')) then
  call default_plot_init (np, plt)
  allocate (plt%graph(1))

  plt%name               = 'floor_plan'
  plt%description        = 'Floor plan drawing of lattice elements.'
  plt%x_axis_type        = 'floor'

  grph => plt%graph(1)
  grph%p => plt
  grph%name                  = 'g'
  grph%box                   = [1, 1, 1, 1]
  grph%type                  = 'floor_plan'
  grph%margin                = qp_rect_struct(0.15, 0.06, 0.05, 0.05, '%BOX')
  grph%scale_margin          = qp_rect_struct(0.02_rp, 0.02_rp, 0.02_rp, 0.02_rp, '%GRAPH')
  grph%draw_only_good_user_data_or_vars = .true.
  grph%x                     = init_axis
  grph%y                     = init_axis
  grph%draw_axes             = .true.
  grph%draw_grid             = .true.
  grph%x%label               = 'SMART LABEL'
  grph%x%major_div_nominal   = 4
  grph%y%label               = 'SMART LABEL'
  grph%y%major_div_nominal   = 4
  grph%y2%major_div_nominal  = 4
  grph%y2_mirrors_y          = .true.

  ! X-ray lines can mainly lie in the vertical plane. 
  ! If so choose as the default the "best" view.

  branch => s%u(1)%model%lat%branch(0)
  if (branch%param%particle == photon$) then
    n = branch%n_ele_track
    dx = maxval(branch%ele(1:n)%floor%r(1)) - minval(branch%ele(1:n)%floor%r(1))
    dy = maxval(branch%ele(1:n)%floor%r(2)) - minval(branch%ele(1:n)%floor%r(2))
    dz = maxval(branch%ele(1:n)%floor%r(3)) - minval(branch%ele(1:n)%floor%r(3))
    if (dx < min(dy, dz)) then
      grph%floor_plan%view = 'yz'
    elseif (dy < min(dx, dz)) then
      grph%floor_plan%view = 'zx'
    else
      grph%floor_plan%view = 'yx'
    endif
  endif

endif

!---------------
! gamma plot

if (all(s%plot_page%template%name /= 'gamma')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'gamma'
  plt%description          = 'Twiss gamma function'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Gamma Function'
  grph%y%label             = '\gg\fn\dA\u, \gg\fn\dB\u [m]'

  crv => grph%curve(1)
  crv%name         = 'a'
  crv%g => grph
  crv%data_type    = 'gamma.a'
  crv%legend_text  = '\gg\fn\dA\u'

  crv => grph%curve(2)
  crv%name         = 'b'
  crv%g => grph
  crv%data_type    = 'gamma.b'
  crv%legend_text  = '\gg\fn\dB\u'
endif

!---------------
! i1

if (all(s%plot_page%template%name /= 'i1')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'i1'
  plt%description          = 'Integrated I1 Radiation integral'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Integrated I1 Radiation Integral'
  grph%y%label       = 'Integrated I1'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'rad_int.i1'
  crv%smooth_line_calc = .false.
endif

!---------------
! i2

if (all(s%plot_page%template%name /= 'i2')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'i2'
  plt%description          = 'Integrated I2 Radiation integral'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Integrated I2 Radiation Integral'
  grph%y%label       = 'Integrated I2'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'rad_int.i2'
  crv%smooth_line_calc = .false.
endif

  !---------------
  ! i3

if (all(s%plot_page%template%name /= 'i3')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'i3'
  plt%description          = 'Integrated I3 Radiation integral'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Integrated I3 Radiation Integral'
  grph%y%label       = 'Integrated I3'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'rad_int.i3'
  crv%smooth_line_calc = .false.
endif

!---------------
! i4a

if (all(s%plot_page%template%name /= 'i4a')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'i4a'
  plt%description          = 'Integrated I4A Radiation integral'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Integrated I4A Radiation Integral'
  grph%y%label       = 'Integrated I4A'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type      = 'rad_int.i4a'
  crv%smooth_line_calc = .false.
endif

!---------------
! i4b

if (all(s%plot_page%template%name /= 'i4b')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'i4b'
  plt%description          = 'Integrated I4B Radiation integral'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Integrated I4B Radiation Integral'
  grph%y%label       = 'Integrated I4B'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'rad_int.i4b'
  crv%smooth_line_calc = .false.
endif

!---------------
! i5a

if (all(s%plot_page%template%name /= 'i5a')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'i5a'
  plt%description          = 'Integrated I5A Radiation integral'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Integrated I5A Radiation Integral'
  grph%y%label       = 'Integrated I5A'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'rad_int.i5a'
  crv%smooth_line_calc = .false.
endif

!---------------
! i5b

if (all(s%plot_page%template%name /= 'i5b')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'i5b'
  plt%description          = 'Integrated I5B Radiation integral'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Integrated I5B Radiation Integral'
  grph%y%label       = 'Integrated I5B'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'rad_int.i5b'
  crv%smooth_line_calc = .false.
endif

!---------------
! Key table

if (all(s%plot_page%template%name /= 'key_table')) then
  call default_plot_init (np, plt)
  allocate (plt%graph(1))

  plt%name           = 'key_table'
  plt%description    = 'Table of keyboard keys bound to variable values'
  plt%x_axis_type    = 'none'

  grph => plt%graph(1)
  grph%p => plt
  grph%name          = 'g'
  grph%box           = [1, 1, 1, 1]
  grph%type          = 'key_table'
  grph%margin        =  qp_rect_struct(0.00, 0.00, 0.03, 0.12, '%BOX')
  grph%x             = init_axis
  grph%y%min         = -1
  grph%y%max         =  1
endif

!---------------
! Lat Layout plot

if (all(s%plot_page%template%name /= 'lat_layout')) then
  call default_plot_init (np, plt)
  allocate (plt%graph(1))

  plt%name           = 'lat_layout'
  plt%description    = 'Lattice elements drawn as a function of S'
  plt%x_axis_type    = 's'

  grph => plt%graph(1)
  grph%p => plt
  grph%name          = 'g'
  grph%box           = [1, 1, 1, 1]
  grph%type          = 'lat_layout'
  grph%margin        =  qp_rect_struct(0.15, 0.06, 0.12, 0.03, '%BOX')
  grph%x             = init_axis
  grph%x%major_div_nominal  = 7
  grph%x%minor_div_max      = 6
  grph%y%min         = -1
  grph%y%max         =  1
endif

!---------------
! kinetic energy

if (all(s%plot_page%template%name /= 'kinetic')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'kinetic'
  plt%description          = 'Kinetic Energy'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Kinetic Energy'
  grph%y%label       = 'K [eV]'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'orbit.kinetic'
endif

!---------------
! Orbit plot

if (all(s%plot_page%template%name /= 'orbit')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name           = 'orbit'
  plt%description    = 'x-y particle orbit'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Orbit'
  grph%y%label       = 'Orbit [mm]'

  crv => grph%curve(1)
  crv%name                = 'x'
  crv%g => grph
  crv%data_type           = 'orbit.x'
  crv%legend_text         = 'X'
  crv%y_axis_scale_factor = 1000

  crv => grph%curve(2)
  crv%name                = 'y'
  crv%g => grph
  crv%data_type           = 'orbit.y'
  crv%legend_text         = 'Y'
  crv%y_axis_scale_factor = 1000
endif

!---------------
! photon_intensity

if (all(s%plot_page%template%name /= 'photon_intensity')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt = default_plot_g1c1 
  plt%name                 = 'photon_intensity'
  plt%description          = 'Photon intensity'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Photon Intensity'
  grph%y%label       = 'Intens'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'photon.intensity'
  crv%smooth_line_calc = .false.
  crv%y_axis_scale_factor = 1
endif

!---------------
! Momentum

if (all(s%plot_page%template%name /= 'momentum')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name           = 'momentum'
  plt%description    = 'Particle momentum'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Particle Momentum'
  grph%y%label       = 'PC [eV]'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type    = 'momentum'
endif

!---------------
! momentum_compaction

if (all(s%plot_page%template%name /= 'momentum_compaction')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'momentum_compaction'
  plt%description          = 'Momentum compaction'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Momentum Compaction'
  grph%y%label       = 'Momentum Compaction [m]'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'momentum_compaction'
endif

!---------------
! phase

if (all(s%plot_page%template%name /= 'phase')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'phase'
  plt%description          = 'Betatron phase'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Betatron Phase'
  grph%y%label       = '\gf\fn\dA\u, \gf\fn\dB\u (deg)'

  crv => grph%curve(1)
  crv%name         = 'a'
  crv%component     = 'model - design'
  crv%g => grph
  crv%data_type    = 'phase.a'
  crv%legend_text  = '\gf\fn\dA\u'
  crv%y_axis_scale_factor = 180/pi

  crv => grph%curve(2)
  crv%name         = 'b'
  crv%component     = 'model - design'
  crv%g => grph
  crv%data_type    = 'phase.b'
  crv%legend_text  = '\gf\fn\dB\u'
  crv%y_axis_scale_factor = 180/pi
endif

!---------------
! ping_a_skew

if (all(s%plot_page%template%name /= 'ping_a_skew')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'ping_a_skew'
  plt%description          = 'Pinged a-mode out-of-plane oscillations'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Pinged A-mode Skew'
  grph%y%label       = 'ping_a: amp_sin_y, amp_cos_y'

  crv => grph%curve(1)
  crv%name         = 'sin_y'
  crv%g => grph
  crv%data_type    = 'ping_a.amp_sin_y'
  crv%legend_text  = 'ping_a.amp_sin_y'

  crv => grph%curve(2)
  crv%name         = 'cos_y'
  crv%g => grph
  crv%data_type    = 'ping_a.amp_cos_y'
  crv%legend_text  = 'ping_a.amp_cos_y'
endif

!---------------
! ping_a_rel_skew

if (all(s%plot_page%template%name /= 'ping_a_rel_skew')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'ping_a_rel_skew'
  plt%description          = 'Pinged a-mode relative out-of-plane oscillations'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Pinged A-mode Rel Skew'
  grph%y%label       = 'ping_a: amp_sin_rel_y, amp_cos_rel_y'

  crv => grph%curve(1)
  crv%name         = 'sin_rel_y'
  crv%g => grph
  crv%data_type    = 'ping_a.amp_sin_rel_y'
  crv%legend_text  = 'ping_a.amp_sin_rel_y'

  crv => grph%curve(2)
  crv%name         = 'cos_rel_y'
  crv%g => grph
  crv%data_type    = 'ping_a.amp_cos_rel_y'
  crv%legend_text  = 'ping_a.amp_cos_rel_y'
endif

!---------------
! ping_a_y_amp_phase

if (all(s%plot_page%template%name /= 'ping_a_y_amp_phase')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'ping_a_y_amp_phase'
  plt%description          = 'Pinged a-mode out-of-plane phase and amplitude'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Pinged A-mode Y-plane Amp & Phase'
  grph%y%label       = 'ping_a: amp_y, phase_y'

  crv => grph%curve(1)
  crv%name         = 'amp_y'
  crv%g => grph
  crv%data_type    = 'ping_a.amp_y'
  crv%legend_text  = 'ping_a.amp_y'

  crv => grph%curve(2)
  crv%name         = 'phase_y'
  crv%g => grph
  crv%data_type    = 'ping_a.phase_y'
  crv%legend_text  = 'ping_a.phase_y'
endif

!---------------
! ping_b_skew

if (all(s%plot_page%template%name /= 'ping_b_skew')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'ping_b_skew'
  plt%description          = 'Pinged b-mode out-of-plane oscillations'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Ping B-mode Skew'
  grph%y%label       = 'ping_b: sin_x, cos_x'

  crv => grph%curve(1)
  crv%name         = 'sin_x'
  crv%g => grph
  crv%data_type    = 'ping_b.amp_sin_x'
  crv%legend_text  = 'ping_b.amp_sin_x'

  crv => grph%curve(2)
  crv%name         = 'cos_x'
  crv%g => grph
  crv%data_type    = 'ping_b.amp_cos_x'
  crv%legend_text  = 'ping_b.amp_cos_x'
endif

!---------------
! ping_b_rel_skew

if (all(s%plot_page%template%name /= 'ping_b_rel_skew')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'ping_b_rel_skew'
  plt%description          = 'Pinged b-mode relative out-of-plane oscillations'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Ping B-mode Rel Skew'
  grph%y%label       = 'ping_b: sin_rel_x, cos_rel_x'

  crv => grph%curve(1)
  crv%name         = 'sin_rel_x'
  crv%g => grph
  crv%data_type    = 'ping_b.amp_sin_rel_x'
  crv%legend_text  = 'ping_b.amp_sin_rel_x'

  crv => grph%curve(2)
  crv%name         = 'cos_rel_x'
  crv%g => grph
  crv%data_type    = 'ping_b.amp_cos_rel_x'
  crv%legend_text  = 'ping_b.amp_cos_rel_x'
endif

!---------------
! ping_b_x_amp_phase

if (all(s%plot_page%template%name /= 'ping_b_x_amp_phase')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'ping_b_x_amp_phase'
  plt%description          = 'Pinged b-mode out-of-plane phase and amplitude'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Pinged B-mode X-plane Amp & Phase'
  grph%y%label       = 'ping_b: amp_x, phase_x'

  crv => grph%curve(1)
  crv%name         = 'amp_x'
  crv%g => grph
  crv%data_type    = 'ping_b.amp_x'
  crv%legend_text  = 'ping_b.amp_x'

  crv => grph%curve(2)
  crv%name         = 'phase_x'
  crv%g => grph
  crv%data_type    = 'ping_b.phase_x'
  crv%legend_text  = 'ping_b.phase_x'
endif

!---------------
! ping_amp

if (all(s%plot_page%template%name /= 'ping_amp')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'ping_amp'
  plt%description          = 'Pinged beam in-plane oscillation amplitudes'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Pinged Beam Osc Amplitudes'
  grph%y%label       = 'A\dax\u, A\dby\u'

  crv => grph%curve(1)
  crv%name         = 'amp_x'
  crv%g => grph
  crv%data_type    = 'ping_a.amp_x'
  crv%legend_text  = 'ping_a.amp_x'

  crv => grph%curve(2)
  crv%name         = 'amp_y'
  crv%g => grph
  crv%data_type    = 'ping_b.amp_y'
  crv%legend_text  = 'ping_b.amp_y'
endif

!---------------
! ping_phase

if (all(s%plot_page%template%name /= 'ping_phase')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'ping_phase'
  plt%description          = 'Pinged beam in-plane oscillation phase'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Pinged Beam Osc Phase'
  grph%y%label       = '\gf\fn\dax\u, \gf\fn\dby\u'

  crv => grph%curve(1)
  crv%name         = 'phase_x'
  crv%component     = 'model - design'
  crv%g => grph
  crv%data_type    = 'ping_a.phase_x'
  crv%legend_text  = 'ping_a.phase_x'
  crv%y_axis_scale_factor = 180/pi

  crv => grph%curve(2)
  crv%name         = 'phase_y'
  crv%component     = 'model - design'
  crv%g => grph
  crv%data_type    = 'ping_b.phase_y'
  crv%legend_text  = 'ping_b.phase_y'
  crv%y_axis_scale_factor = 180/pi
endif

!---------------
! px_py

if (all(s%plot_page%template%name /= 'px_py')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'px_py'
  plt%description          = 'Particle px and py phase space coordinates'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Particle px, py Phase Space Coords'
  grph%y%label       = 'px, py'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'orbit.px'
  crv%legend_text         = 'px'

  crv => grph%curve(2)
  crv%g => grph
  crv%data_type     = 'orbit.py'
  crv%legend_text         = 'py'
endif

!---------------
! pz

if (all(s%plot_page%template%name /= 'pz')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'pz'
  plt%description          = 'Particle pz phase space coordinate'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Particle pz Phase Space Coordinate'
  grph%y%label       = 'pz'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'orbit.pz'
  crv%legend_text         = 'pz'
endif

!---------------
! Reference momentom

if (all(s%plot_page%template%name /= 'ref_p0c')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'ref_p0c'
  plt%description          = 'Reference momentum'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Reference Momentum'
  grph%y%label       = 'Ref p0c'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'ref_p0c'

endif

!---------------
! Reference time

if (all(s%plot_page%template%name /= 'ref_time')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'ref_time'
  plt%description          = 'Reference time'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Reference Time'
  grph%y%label       = 'Ref Time'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'ref_time'

endif

!---------------
! spin dn/dpz

if (all(s%plot_page%template%name /= 'spin_dn_dpz')) then
  call default_plot_init (np, plt, default_plot_g1c4)
  plt%name                 = 'spin_dn_dpz'
  plt%description          = 'Spin dn/dpz (x,y,z) components & amplitude'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Spin dn/dpz (x,y,z) & Amplitude'
  grph%y%label       = 'dn/dpz X, Y, Z, Amp'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'spin_dn_dpz.x'
  crv%smooth_line_calc = .false.

  crv => grph%curve(2)
  crv%g => grph
  crv%data_type     = 'spin_dn_dpz.y'
  crv%smooth_line_calc = .false.

  crv => grph%curve(3)
  crv%g => grph
  crv%data_type     = 'spin_dn_dpz.z'
  crv%smooth_line_calc = .false.

  crv => grph%curve(4)
  crv%g => grph
  crv%data_type     = 'spin_dn_dpz.amp'
  crv%smooth_line_calc = .false.
endif

!---------------
! spin xyz

if (all(s%plot_page%template%name /= 'spin_xyz')) then
  call default_plot_init (np, plt, default_plot_g1c3)
  plt%name                 = 'spin_xyz'
  plt%description          = 'Spin x, y, z components'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Spin x, y, z'
  grph%y%label       = 'X, Y, Z'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'spin.x'

  crv => grph%curve(2)
  crv%g => grph
  crv%data_type     = 'spin.y'

  crv => grph%curve(3)
  crv%g => grph
  crv%data_type     = 'spin.z'
endif

!---------------
! sr energy loss

if (all(s%plot_page%template%name /= 'sr_energy_loss')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'sr_energy_loss'
  plt%description          = 'Synch Radiation energy loss'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Synch Rad Energy Loss [Rad_Integral:I2_E4 * r_e * mc\u2\d * 2 / 3]'
  grph%y%label       = 'E\dLoss\u [MeV * m]'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type        = 'rad_int.i2_e4'
  crv%y_axis_scale_factor = 9.59976e-16 ! (2/3) * r_e *mec2 in MeV*m
  crv%smooth_line_calc = .false.
endif

!---------------
! time

if (all(s%plot_page%template%name /= 'time')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'time'
  plt%description          = 'particle time'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Particle Time (sec)'
  grph%y%label       = 'Time [sec]'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'time'
endif

!---------------
! velocity

if (all(s%plot_page%template%name /= 'velocity')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'velocity'
  plt%description          = 'Velocity/c_light'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Particle Velocity/c'
  grph%y%label       = 'v/c [mm]'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'velocity'
  crv%y_axis_scale_factor = 1
endif

!---------------
! w-function

if (all(s%plot_page%template%name /= 'w_function')) then
  call default_plot_init (np, plt, default_plot_g1c2)
  plt%name                 = 'w_function'
  plt%description          = 'Chromatic W-function'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'W-function'
  grph%y%label       = 'W\da\u, W\db\u'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'chrom.w.a'

  crv => grph%curve(2)
  crv%g => grph
  crv%data_type     = 'chrom.w.b'
endif

!---------------
! z

if (all(s%plot_page%template%name /= 'z')) then
  call default_plot_init (np, plt, default_plot_g1c1)
  plt%name                 = 'z'
  plt%description          = 'Particle z position'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Particle Z-position'
  grph%y%label       = 'Z [mm]'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'orbit.z'
  crv%y_axis_scale_factor = 1000
endif

!-------------------------------------------------
! Scratch plot

plt => s%plot_page%template(size(s%plot_page%template))

nullify(plt%r)
if (allocated(plt%graph)) deallocate (plt%graph)
allocate (plt%graph(1))
plt%graph(1)%p => plt
plt%graph(1)%name = 'g'
plt%name = 'scratch'

!-------------------------------------------------
! Regions

if (.not. allocated(s%plot_page%region)) then
  allocate (s%plot_page%region(130))
  nr = 0
else
  nr = size(s%plot_page%region)
  call move_alloc (s%plot_page%region, temp_region)
  allocate (s%plot_page%region(nr + 130))
  s%plot_page%region(1:nr) = temp_region
  deallocate(temp_region)
endif

y_layout = 0.15

if (all(s%plot_page%region(:)%name /= 'layout')) then
  nr = nr + 1
  s%plot_page%region(nr)%name = 'layout'
  s%plot_page%region(nr)%location = [0.0_rp, 1.0_rp, 0.0_rp, y_layout]
endif

do i = 1, 5
  do j = 1, i
    write (name, '(a, 2i0)') 'r', j, i
    if (any(s%plot_page%region(:)%name == name)) cycle
    nr = nr + 1
    s%plot_page%region(nr)%name = name
    y1 = y_layout + (1 - y_layout) * real(i-j)/ i
    y2 = y_layout + (1 - y_layout) * real(i-j+1) / i
    s%plot_page%region(nr)%location = [0.0_rp, 1.0_rp, y1, y2]
    if (i > 4) s%plot_page%region(nr)%list_with_show_plot_command = .false.
  enddo
enddo

!

do k1 = 2, 4

  do i = 1, k1
    write (name, '(a, 2i0)') 'layout', i, k1
    if (all(s%plot_page%region(:)%name /= name)) then
      nr = nr + 1
      if (k1 > 2) s%plot_page%region(nr)%list_with_show_plot_command = .false.
      s%plot_page%region(nr)%name = name
      s%plot_page%region(nr)%location = [real(i-1, rp)/k1, real(i, rp)/k1, 0.0_rp, y_layout]
    endif
  enddo

  do k2 = 1, 4
    do i = 1, k1
    do j = 1, k2
      write (name, '(a, 4i0)') 'r', i, j, k1, k2
      if (any(s%plot_page%region(:)%name == name)) cycle
      nr = nr + 1
      s%plot_page%region(nr)%name = name
      if (k1 > 2) s%plot_page%region(nr)%list_with_show_plot_command = .false.
      x1 = real(i-1)/ k1
      x2 = real(i) / k1
      y1 = y_layout + (1 - y_layout) * real(k2-j)/ k2
      y2 = y_layout + (1 - y_layout) * real(k2-j+1) / k2
      s%plot_page%region(nr)%location = [x1, x2, y1, y2]
    enddo
    enddo
  enddo

enddo

! For parametric plots

do i = 0, 9
  nr = size(s%plot_page%region) - i
  s%plot_page%region(nr)%name = 'scratch' // int_str(i)
  s%plot_page%region(nr)%list_with_show_plot_command = .false.
  s%plot_page%region(nr)%location = [0, 1, 0, 1]
enddo

!

if (all (place(:)%region == '')) then
  branch => s%u(1)%model%lat%branch(0)
  if (branch%param%particle == photon$) then
    call tao_place_cmd ('r13', 'floor_plan')
    call tao_place_cmd ('r23', 'photon_intensity')
    call tao_place_cmd ('r33', 'orbit')
    call tao_place_cmd ('layout', 'lat_layout')
  else  ! Charged particle
    call tao_place_cmd ('r13', 'beta')
    call tao_place_cmd ('r23', 'dispersion')
    call tao_place_cmd ('r33', 'orbit')
    call tao_place_cmd ('layout', 'lat_layout')
  endif
endif

end subroutine tao_setup_default_plotting

!----------------------------------------------------------------------------------------
! contains

subroutine default_plot_init(np, plt, default_plot)

type (tao_plot_struct), pointer :: plt
type (tao_plot_struct), optional :: default_plot
integer np

!

np = np + 1
plt => s%plot_page%template(np)
plt%default_plot = .true.

nullify(plt%r)
if (allocated(plt%graph)) deallocate (plt%graph)

if (present(default_plot)) plt = default_plot

end subroutine default_plot_init

!----------------------------------------------------------------------------------------
! contains

subroutine transfer_this_axis (ax_in, ax_out)

type (qp_axis_struct) ax_in, ax_out

if (ax_in%label       /= null_name$) ax_out%label       = ax_in%label
if (ax_in%label_color /= null_name$) ax_out%label_color = ax_in%label_color
if (ax_in%type        /= null_name$) ax_out%type        = ax_in%type
if (ax_in%bounds      /= null_name$) ax_out%bounds      = ax_in%bounds

if (ax_in%min            /= real_garbage$) ax_out%min            = ax_in%min
if (ax_in%max            /= real_garbage$) ax_out%eval_max       = ax_in%max
if (ax_in%min            /= real_garbage$) ax_out%eval_min       = ax_in%min
if (ax_in%max            /= real_garbage$) ax_out%max            = ax_in%max
if (ax_in%tick_min       /= real_garbage$) ax_out%tick_min       = ax_in%tick_min
if (ax_in%tick_max       /= real_garbage$) ax_out%tick_max       = ax_in%tick_max
if (ax_in%dtick          /= real_garbage$) ax_out%dtick          = ax_in%dtick
if (ax_in%number_offset  /= real_garbage$) ax_out%number_offset  = ax_in%number_offset
if (ax_in%label_offset   /= real_garbage$) ax_out%label_offset   = ax_in%label_offset
if (ax_in%major_tick_len /= real_garbage$) ax_out%major_tick_len = ax_in%major_tick_len
if (ax_in%minor_tick_len /= real_garbage$) ax_out%minor_tick_len = ax_in%minor_tick_len

if (ax_in%major_div         /= int_garbage$) ax_out%major_div         = ax_in%major_div
if (ax_in%major_div_nominal /= int_garbage$) ax_out%major_div_nominal = ax_in%major_div_nominal
if (ax_in%minor_div         /= int_garbage$) ax_out%minor_div         = ax_in%minor_div
if (ax_in%minor_div_max     /= int_garbage$) ax_out%minor_div_max     = ax_in%minor_div_max
if (ax_in%places            /= int_garbage$) ax_out%places            = ax_in%places
if (ax_in%tick_side         /= int_garbage$) ax_out%tick_side         = ax_in%tick_side
if (ax_in%number_side       /= int_garbage$) ax_out%number_side       = ax_in%number_side

if (.not. ax_in%draw_label)   ax_out%draw_label   = ax_in%draw_label
if (.not. ax_in%draw_numbers) ax_out%draw_numbers = ax_in%draw_numbers

end subroutine transfer_this_axis

end subroutine tao_init_plotting

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------

function old_style_title_syntax(iu) result (is_old_style)

use tao_input_struct

type (tao_plot_page_test_input) plot_page
type (tao_plot_input) default_plot
type (tao_graph_input) :: default_graph
type (tao_region_input) region(100)
type (tao_place_input) place(30)

integer iu, ios
logical is_old_style, include_default_plots

namelist / tao_plot_page / plot_page, default_plot, default_graph, region, place, include_default_plots

!

rewind (iu)
read (iu, nml = tao_plot_page, iostat = ios) 

is_old_style = (ios == 0)

end function old_style_title_syntax

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine tao_read_in_patterns(iu, plot_file)

use tao_struct

implicit none

type (tao_shape_pattern_struct), allocatable :: temp_pat(:)
type (tao_shape_pattern_struct), pointer :: pat
type (qp_line_struct) :: line 
type (tao_shape_pattern_point_struct) :: pt(30)

integer iu, ios, nn, j, jc, jpt, nc, npt

character(40) name
character(8) :: scale   ! This is no longer used.
character(*) plot_file
character(*), parameter :: r_name = 'tao_read_in_patterns'

namelist / shape_pattern / name, line, pt, scale

!

rewind (iu)
if (allocated(s%plot_page%pattern)) deallocate(s%plot_page%pattern)
allocate (s%plot_page%pattern(0))

do  ! Loop over all patterns
  line  = qp_line_struct(1, '', 'solid')
  pt    = tao_shape_pattern_point_struct()
  name = ''

  read (iu, nml = shape_pattern, iostat = ios) 
  if (ios < 0) exit
  if (ios > 0) then
    call out_io (s_error$, r_name, 'ERROR READING SHAPE_PATTERN NAMELIST.', 'IN FILE: ' // plot_file)
    rewind (iu)
    do
      read (iu, nml = shape_pattern)  ! To generate error message
    enddo
  endif

  !

  nn = size(s%plot_page%pattern)
  call move_alloc(s%plot_page%pattern, temp_pat)
  allocate (s%plot_page%pattern(nn+1))

  do j = 1, nn
    pat => s%plot_page%pattern(j)
    npt = size(temp_pat(j)%pt)
    allocate (pat%pt(npt))
    pat = temp_pat(j)
  enddo

  deallocate (temp_pat)

  !

  pat => s%plot_page%pattern(nn+1)
  pat%name = name

  do jpt = size(pt), 1, -1
    if (pt(jpt)%s /= real_garbage$) exit
  enddo
  allocate (pat%pt(jpt))
  pat%line = line
  pat%pt   = pt(1:jpt)
enddo

end subroutine tao_read_in_patterns

