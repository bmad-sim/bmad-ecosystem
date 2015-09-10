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

use tao_input_struct
use tao_plot_window_mod, dummy => tao_init_plotting
use quick_plot

implicit none

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
type (tao_curve_input) curve(n_curve_maxx), curve1, curve2, curve3, curve4
type (tao_place_input) place(10)
type (old_tao_ele_shape_struct) shape(20)
type (tao_ele_shape_struct) ele_shape(20)
type (tao_ele_shape_struct), pointer :: e_shape
type (qp_symbol_struct) default_symbol
type (qp_line_struct) default_line
type (qp_axis_struct) init_axis
type (ele_pointer_struct), allocatable, save :: eles(:)

real(rp) y1, y2

integer iu, i, j, k, ix, ip, n, ng, ios, ios1, ios2, i_uni
integer graph_index, color, i_graph

character(*) plot_file_in
character(len(plot_file_in)) plot_file_array
character(100) plot_file, graph_name, full_file_name
character(80) label
character(40) str
character(20) :: r_name = 'tao_init_plotting'

logical err, include_default_plots

namelist / tao_plot_page / plot_page, default_plot, default_graph, region, place, &
                include_default_plots
namelist / tao_template_plot / plot, default_graph
namelist / tao_template_graph / graph, graph_index, curve, curve1, curve2, curve3, curve4

namelist / floor_plan_drawing / ele_shape
namelist / lat_layout_drawing / ele_shape

! These are old style

namelist / element_shapes / shape
namelist / element_shapes_floor_plan / ele_shape
namelist / element_shapes_lat_layout / ele_shape

! See if this routine has been called before

call qp_init_com_struct()  ! Init quick_plot
if (.not. s%global%init_plot_needed) return
s%global%init_plot_needed = .false.

! Init

init_axis%min = 0
init_axis%max = 0

place%region = ' '
region%name  = ' '       ! a region exists only if its name is not blank 
include_default_plots = .false.

plot_page = plot_page_default
plot_page%title(:)%draw_it = .false.
plot_page%title(:)%string = ' '
plot_page%title(:)%justify = 'CC'
plot_page%title(:)%x = 0.50
plot_page%title(:)%y = 0.990
plot_page%title(1)%y = 0.996
plot_page%title(2)%y = 0.97
plot_page%title(:)%units = '%PAGE'
plot_page%size = [500, 600]

default_plot%name = ' '
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

default_graph = tao_graph_input()
default_graph%x                     = init_axis
default_graph%y                     = init_axis
default_graph%y%major_div           = -1
default_graph%y%major_div_nominal   = -1
default_graph%y2                    = init_axis
default_graph%y2%major_div          = -1
default_graph%y2%major_div_nominal  = -1
default_graph%y2%label_color        = blue$
default_graph%y2%draw_numbers       = .false.
default_graph%margin                = qp_rect_struct(0.15, 0.06, 0.12, 0.12, '%BOX')

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
call tao_open_file (plot_file, iu, full_file_name, s_fatal$)
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

call tao_set_plotting (plot_page, s%plot_page, .true.)

! title

forall (i = 1:size(s%plot_page%title), (s%plot_page%title(i)%string .ne. ' ')) &
            s%plot_page%title(i)%draw_it = .true.

! Plot window geometry specified on cmd line?

if (s%com%plot_geometry /= '') then
   str = s%com%plot_geometry
   ix = index(str, 'x')
   if (ix == 0) then
     call out_io (s_error$, r_name, 'Malformed -geometry argument. No "x" present: ' // str)
   else
     if (.not. is_integer(str(1:ix-1)) .or. .not. is_integer(str(ix+1:))) then
       call out_io (s_error$, r_name, 'Malformed -geometry argument: ' // str)
     else
       read (str(:ix-1), *) plot_page%size(1)
       read (str(ix+1:), *) plot_page%size(2)
     endif
   endif
 endif
 
! allocate a s%plot_page%plot structure for each region defined and
! transfer the info from the input region structure.

n = count(region%name /= ' ')
allocate (s%plot_page%region(n))

do i = 1, n
  s%plot_page%region(i)%name     = region(i)%name
  s%plot_page%region(i)%location = region(i)%location
enddo

!-----------------------------------------------------------------------------------
! Read in element shapes
! Look for old style namelist 

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
  call out_io (s_warn$, r_name, 'ELEMNET_SHAPES NAMELIST IS DEPRECATED.', &
                                'PLEASE CONVERT TO FLOOR_PLAN_DRAWING AND/OR LAT_LAYOUT_DRAWING NAMELISTS.')
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
    call tao_string_to_element_id (ele_shape(i)%ele_id, ele_shape(i)%ix_ele_key, ele_shape(i)%name_ele, err, .true.)
  enddo

  call tao_uppercase_shapes (ele_shape, n, 'f')
  allocate (s%plot_page%floor_plan%ele_shape(n), s%plot_page%lat_layout%ele_shape(n))
  s%plot_page%floor_plan%ele_shape = ele_shape(1:n)
  s%plot_page%lat_layout%ele_shape = ele_shape(1:n)
endif

! Look for new style shape namelist if could not find old style

if (ios < 0) then

  ! Read floor_plan_drawing namelist

  ele_shape(:)%ele_id     = ''
  ele_shape(:)%label      = 'name'
  ele_shape(:)%draw       = .true.

  rewind (iu)
  read (iu, nml = element_shapes_floor_plan, iostat = ios1)  ! Deprecated name
  rewind (iu)
  read (iu, nml = floor_plan_drawing, iostat = ios2)

!  if (ios1 >= 0) then
!    call out_io (s_warn$, r_name, &
!            'Note: The "element_shapes_floor_plan" namelist has been renamed to', &
!            '      "floor_plan_drawing to reflect the fact that this namelist  ', &
!            '      now is used to specify more than element shapes. Please     ', &
!            '      make the appropriate change in your input file.             ', &
!            'For now, Tao will accept the old namelist name...                 ')
!  endif

  if (ios1 > 0) then 
    rewind (iu)
    call out_io (s_error$, r_name, 'ERROR READING ELEMENT_SHAPES_FLOOR_PLAN NAMELIST')
    read (iu, nml = element_shapes_floor_plan)  ! To generate error message
  endif

  if (ios2 > 0) then 
    rewind (iu)
    call out_io (s_error$, r_name, 'ERROR READING FLOOR_PLAN_DRAWING NAMELIST')
    read (iu, nml = floor_plan_drawing)
  endif

  call tao_uppercase_shapes (ele_shape, n, 'f')
  allocate (s%plot_page%floor_plan%ele_shape(n))
  s%plot_page%floor_plan%ele_shape = ele_shape(1:n)

  ! Read element_shapes_lat_layout namelist

  ele_shape(:)%ele_id     = ''
  ele_shape(:)%label      = 'name'
  ele_shape(:)%draw       = .true.

  rewind (iu)
  read (iu, nml = element_shapes_lat_layout, iostat = ios1)
  rewind (iu)
  read (iu, nml = lat_layout_drawing, iostat = ios2)

!  if (ios1 == 0) then
!    call out_io (s_warn$, r_name, &
!            'Note: The "element_shapes_lattice_list" namelist has been renamed to', &
!            '      "lattice_list_drawing to reflect the fact that this namelist  ', &
!            '      now is used to specify more than element shapes. Please       ', &
!            '      make the appropriate change in your input file.               ', &
!            'For now, Tao will accept the old namelist name...                   ')
!  endif

  if (ios1 == 0) then
    ele_shape(:)%size = ele_shape(:)%size * 1.0 / 40.0 ! scale to current def.
  endif 

  if (ios1 > 0) then 
    rewind (iu)
    call out_io (s_error$, r_name, 'ERROR READING ELEMENT_SHAPES_LAT_LAYOUT NAMELIST')
    read (iu, nml = element_shapes_lat_layout)  ! To generate error message
  endif

  if (ios2 > 0) then 
    rewind (iu)
    call out_io (s_error$, r_name, 'ERROR READING LAT_LAYOUT_DRAWING NAMELIST')
    read (iu, nml = lat_layout_drawing)
  endif

  call tao_uppercase_shapes (ele_shape, n, 'l')
  allocate (s%plot_page%lat_layout%ele_shape(n))
  s%plot_page%lat_layout%ele_shape  = ele_shape(1:n)

endif

! Error check

if (allocated(s%plot_page%lat_layout%ele_shape)) then
  do i = 1, size(s%plot_page%lat_layout%ele_shape)
    e_shape => s%plot_page%lat_layout%ele_shape(i)
    if (e_shape%ele_id(1:6) == 'wall::') cycle
    select case (e_shape%shape)
    case ('BOX', 'VAR_BOX', 'ASYM_VAR_BOX', 'XBOX', 'DIAMOND', 'BOW_TIE', 'CIRCLE', 'X', 'NONE')
    case default
      call out_io (s_fatal$, r_name, 'ERROR: UNKNOWN ELE_SHAPE: ' // e_shape%shape)
      call err_exit
    end select
  enddo
endif

if (allocated(s%plot_page%floor_plan%ele_shape)) then
  do i = 1, size(s%plot_page%floor_plan%ele_shape)
    e_shape => s%plot_page%floor_plan%ele_shape(i)
    select case (e_shape%shape)
    case ('BOX', 'VAR_BOX', 'ASYM_VAR_BOX', 'XBOX', 'DIAMOND', 'BOW_TIE', 'CIRCLE', 'X', 'NONE')
    case default
      call out_io (s_fatal$, r_name, 'ERROR: UNKNOWN ELE_SHAPE: ' // e_shape%shape)
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
  call tao_open_file (plot_file, iu, full_file_name, s_fatal$)

  do   ! Loop over templates in a file
    read (iu, nml = tao_template_plot, iostat = ios, err = 9100)  
    call out_io (s_blank$, r_name, 'Init: Read tao_template_plot ' // plot%name)
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
plt%graph(1)%name = 'g1'

! Now read in the plots

ip = 0   ! template plot index
plot_file_array = plot_file_in

do  ! Loop over plot files

  call string_trim(plot_file_array, plot_file_array, ix)
  if (ix == 0) exit
  plot_file = plot_file_array(1:ix)
  plot_file_array = plot_file_array(ix+1:)
  call out_io (s_blank$, r_name, '*Init: Opening Plotting File: ' // plot_file)
  call tao_open_file (plot_file, iu, full_file_name, s_fatal$)

  do   ! Loop over templates in a file

    plot = default_plot
    default_graph = master_default_graph

    read (iu, nml = tao_template_plot, iostat = ios, err = 9100)  
    if (ios /= 0) exit

    call out_io (s_blank$, r_name, 'Init: Read tao_template_plot namelist: ' // plot%name)
    do i = 1, ip
      if (plot%name == s%plot_page%template(ip)%name) then
        call out_io (s_error$, r_name, 'DUPLICATE PLOT NAME: ' // plot%name)
        exit
      endif
    enddo

    ip = ip + 1

    plt => s%plot_page%template(ip)
    nullify(plt%r)
    plt%name                 = plot%name
    plt%description          = plot%description
    plt%x_axis_type          = plot%x_axis_type
    plt%x                    = plot%x
    plt%autoscale_gang_x     = plot%autoscale_gang_x 
    plt%autoscale_gang_y     = plot%autoscale_gang_y 
    plt%autoscale_x          = plot%autoscale_x 
    plt%autoscale_y          = plot%autoscale_y 

    if (plt%x%major_div < 0 .and. plt%x%major_div_nominal < 0) plt%x%major_div_nominal = 6

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
      write (graph%name, '(a, i0)') 'g', i_graph
      curve(:) = tao_curve_input()
      do j = 1, size(curve)
        write (curve(j)%name, '(a, i0)') 'c', j
      enddo
      if (plt%x_axis_type == 's' .or. plt%x_axis_type == 'lat' .or. plt%x_axis_type == 'var') then
        curve(:)%draw_symbols = .false.
      else
        curve(:)%draw_symbols = .true.
      endif
      curve(2:7)%symbol%type = [times_sym$, square_sym$, plus_sym$, triangle_sym$, x_symbol_sym$, diamond_sym$]
      curve(2:7)%symbol%color = [blue$, red$, green$, cyan$, magenta$, yellow$]
      curve(2:7)%line%color = curve(2:7)%symbol%color
      ! to get around gfortran compiler bug.
      curve1 = curve(1); curve2 = curve(2); curve3 = curve(3); curve4 = curve(4)

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
      grph%symbol_size_scale   = graph%symbol_size_scale   
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
      grph%scale_margin          = graph%scale_margin
      grph%x                     = graph%x
      grph%y                     = graph%y
      grph%y2                    = graph%y2
      grph%ix_universe           = graph%ix_universe
      grph%ix_branch             = graph%ix_branch
      grph%clip                  = graph%clip
      grph%draw_axes             = graph%draw_axes
      grph%draw_grid             = graph%draw_grid
      grph%correct_xy_distortion = graph%correct_xy_distortion
      grph%draw_only_good_user_data_or_vars = &
                                   graph%draw_only_good_user_data_or_vars
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

      if (.not. s%com%common_lattice .and. grph%ix_universe == 0) then
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

      if (grph%type == 'floor_plan' .and. .not. allocated (s%plot_page%floor_plan%ele_shape)) &
                call out_io (s_error$, r_name, 'NO ELEMENT SHAPES DEFINED FOR FLOOR_PLAN PLOT.')
   
      if (grph%type == 'lat_layout') then
        if (.not. allocated (s%plot_page%lat_layout%ele_shape)) call out_io (s_error$, r_name, &
                              'NO ELEMENT SHAPES DEFINED FOR LAT_LAYOUT PLOT.')
        if (plt%x_axis_type /= 's') call out_io (s_error$, r_name, &
                              'A LAT_LAYOUT MUST HAVE X_AXIS_TYPE = "s" FOR A VISIBLE PLOT!')
        plt%autoscale_gang_y = .false.  ! True does not make sense.
      endif

      if (graph%n_curve == 0) then
        if (allocated(grph%curve)) deallocate (grph%curve)
      else
        allocate (grph%curve(graph%n_curve))
      endif

      do j = 1, graph%n_curve
        crv => grph%curve(j)

        select case (j)
        case (1); if (curve1%data_type /= '') curve(1) = curve1
        case (2); if (curve2%data_type /= '') curve(2) = curve2
        case (3); if (curve3%data_type /= '') curve(3) = curve3
        case (4); if (curve4%data_type /= '') curve(4) = curve4
        end select

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
        crv%hist                 = curve(j)%hist

        ! Convert old syntax to new

        if (crv%data_source == 'beam_tracking') crv%data_source = 'beam'
        if (crv%data_source == 'lattice')       crv%data_source = 'lat'
        if (crv%data_source == 'data_array')    crv%data_source = 'data'
        if (crv%data_source == 'dat')           crv%data_source = 'data'
        if (crv%data_source == 'var_array')     crv%data_source = 'var'

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

        if (.not. s%com%common_lattice .and. crv%ix_universe == 0) then
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
          call tao_locate_elements (crv%ele_ref_name, i_uni, eles, err, ignore_blank = .true.) ! find the index
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

! If no plots have been defined or default plots wanted then use default

if (ip == 0 .or. include_default_plots) then
  if (size(s%plot_page%lat_layout%ele_shape) == 0) deallocate(s%plot_page%lat_layout%ele_shape)
  if (size(s%plot_page%floor_plan%ele_shape) == 0) deallocate(s%plot_page%floor_plan%ele_shape)
  if (size(s%plot_page%region) == 0) deallocate (s%plot_page%region)
  call tao_setup_default_plotting()
endif

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
        'TAO_TEMPLATE_PLOT NAMELIST READ ERROR.', 'IN FILE: ' // full_file_name)
rewind (iu)
do
  read (iu, nml = tao_template_plot)  ! force printing of error message
enddo

!-----------------------------------------

9200 continue
call out_io (s_error$, r_name, &
       'TAO_TEMPLATE_GRAPH NAMELIST READ ERROR.', 'IN FILE: ' // full_file_name)
rewind (iu)
do
  read (iu, nml = tao_template_graph)  ! force printing of error message
enddo

!----------------------------------------------------------------------------------------
contains

subroutine tao_uppercase_shapes (ele_shape, n_shape, prefix)

type (tao_ele_shape_struct), target :: ele_shape(:)
type (tao_ele_shape_struct), pointer :: s
integer n, n_shape
character(1) prefix

!

n_shape = 0
do n = 1, size(ele_shape)
  s => ele_shape(n)
  ! Bmad wants ele names upper case but Tao data is case sensitive.
  if (s%ele_id(1:5) /= 'dat::' .and. s%ele_id(1:5) /= 'var::' .and. &
      s%ele_id(1:5) /= 'lat::' .and. s%ele_id(1:6) /= 'wall::') call str_upcase (s%ele_id, s%ele_id)
  call str_upcase (s%shape,    s%shape)
  call str_upcase (s%color,    s%color)
  call downcase_string (s%label)
  if (s%label == '') s%label = 'name'
  if (index('false', trim(s%label)) == 1) s%label = 'none'
  if (index('true', trim(s%label)) == 1) s%label = 'name'
  ! Convert old class:name format to new class::name format
  ix = index(s%ele_id, ":")
  if (ix /= 0 .and. s%ele_id(ix+1:ix+1) /= ':') &
     s%ele_id = s%ele_id(1:ix) // ':' // s%ele_id(ix+1:)
  call tao_string_to_element_id (s%ele_id, s%ix_ele_key, s%name_ele, err, .true.)
  if (s%ele_id /= '') n_shape = n
enddo

end subroutine tao_uppercase_shapes

!----------------------------------------------------------------------------------------
! contains

subroutine tao_setup_default_plotting()

type (branch_struct), pointer :: branch
type (tao_plot_struct), target :: default_plot_g1c1, default_plot_g1c2, default_plot_g2c1, default_plot_g1_c4
type (tao_plot_struct), allocatable :: temp_template(:)
type (tao_plot_region_struct), allocatable :: temp_region(:)
type (tao_ele_shape_struct) :: dflt_lat_layout(25) = [&
      tao_ele_shape_struct('FORK::*',              'CIRCLE', 'RED',     0.15_rp, 'name', .true.,   fork$, '*'), &
      tao_ele_shape_struct('CRYSTAL::*',           'CIRCLE', 'RED',     0.15_rp, 'name', .true.,   crystal$, '*'), &
      tao_ele_shape_struct('DETECTOR::*',          'BOX',    'BLACK',   0.30_rp, 'name', .true.,   detector$, '*'), &
      tao_ele_shape_struct('DIFFRACTION_PLATE::*', 'BOX',    'CYAN',    0.30_rp, 'name', .true.,   diffraction_plate$, '*'), &
      tao_ele_shape_struct('E_GUN::*',             'XBOX',   'RED',     0.40_rp, 'name', .true.,   e_gun$, '*'), &
      tao_ele_shape_struct('EM_FIELD::*',          'XBOX',   'BLUE',    0.40_rp, 'name', .true.,   em_field$, '*'), &
      tao_ele_shape_struct('ECOLLIMATOR::*',       'XBOX',   'BLUE',    0.20_rp, 'name', .false.,  ecollimator$, '*'), &
      tao_ele_shape_struct('INSTRUMENT::*',        'BOX',    'BLUE',    0.30_rp, 'name', .false.,  instrument$, '*'), &
      tao_ele_shape_struct('LCAVITY::*',           'XBOX',   'RED',     0.50_rp, 'none', .true.,   lcavity$, '*'), &
      tao_ele_shape_struct('MARKER::*',            'BOX',    'BLUE',    0.30_rp, 'name', .false.,  marker$, '*'), &
      tao_ele_shape_struct('MIRROR::*',            'CIRCLE', 'RED',     0.15_rp, 'name', .true.,   mirror$, '*'), &
      tao_ele_shape_struct('MONITOR::*',           'BOX',    'BLACK',   0.30_rp, 'name', .false.,  monitor$, '*'), &
      tao_ele_shape_struct('MULTILAYER_MIRROR::*', 'CIRCLE', 'RED',     0.15_rp, 'name', .true.,   multilayer_mirror$, '*'), &
      tao_ele_shape_struct('OCTUPOLE::*',          'BOX',    'BLACK',   0.30_rp, 'name', .false.,  octupole$, '*'), &
      tao_ele_shape_struct('PHOTON_FORK::*',       'CIRCLE', 'RED',     0.15_rp, 'name', .true.,   photon_fork$, '*'), &
      tao_ele_shape_struct('QUADRUPOLE::*',        'XBOX',   'MAGENTA', 0.37_rp, 'name', .true.,   quadrupole$, '*'), &
      tao_ele_shape_struct('RCOLLIMATOR::*',       'XBOX',   'BLUE',    0.20_rp, 'name', .false.,  rcollimator$, '*'), &
      tao_ele_shape_struct('RFCAVITY::*',          'XBOX',   'RED',     0.50_rp, 'name', .true.,   rfcavity$, '*'), &
      tao_ele_shape_struct('SAMPLE::*',            'BOX',    'BLACK',   0.30_rp, 'name', .true.,   sample$, '*'), &
      tao_ele_shape_struct('SBEND::*',             'BOX',    'BLACK',   0.20_rp, 'none', .true.,   sbend$, '*'), &
      tao_ele_shape_struct('SEXTUPOLE::*',         'XBOX',   'GREEN',   0.37_rp, 'none', .true.,   sextupole$, '*'), &
      tao_ele_shape_struct('SOL_QUAD::*',          'BOX',    'BLACK',   0.40_rp, 'name', .false.,  sol_quad$, '*'), &
      tao_ele_shape_struct('SOLENOID::*',          'BOX',    'BLUE',    0.30_rp, 'name', .true.,   solenoid$, '*'), &
      tao_ele_shape_struct('WIGGLER::*',           'XBOX',   'CYAN',    0.50_rp, 'name', .true.,   wiggler$, '*'), &
      tao_ele_shape_struct('PHOTON_INIT::*',       'BOX',    'BLACK',   0.30_rp, 'name', .true.,   photon_init$, '*') ]

real(rp) y_layout, dx, dy, dz
integer np, n, nr
character(40) name

!

call tao_set_plotting (plot_page, s%plot_page, .true.)

n = size(dflt_lat_layout)

if (.not. allocated(s%plot_page%lat_layout%ele_shape)) then
  allocate (s%plot_page%lat_layout%ele_shape(30))
  s%plot_page%lat_layout%ele_shape(:)%ele_id = ''
  s%plot_page%lat_layout%ele_shape(1:n) = dflt_lat_layout
endif

if (.not. allocated(s%plot_page%floor_plan%ele_shape)) then
  allocate (s%plot_page%floor_plan%ele_shape(30))
  s%plot_page%floor_plan%ele_shape(:)%ele_id = ''
  s%plot_page%floor_plan%ele_shape(1:n) = dflt_lat_layout
  s%plot_page%floor_plan%ele_shape%size = 20 * s%plot_page%floor_plan%ele_shape%size
endif

!---------------------------------

if (allocated(s%plot_page%template)) then
  n = size(s%plot_page%template)
  call move_alloc(s%plot_page%template, temp_template)
  allocate (s%plot_page%template(n + 29))
  s%plot_page%template(1:n) = temp_template
  deallocate (temp_template)
  np = n
  if (s%plot_page%template(np)%name == 'scratch') np = np - 1
else
  allocate (s%plot_page%template(29))
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
plt%x                    = init_axis
plt%x%major_div_nominal  = 8
plt%x%minor_div_max      = 6

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
grph%component            = 'model'
grph%x                    = init_axis
grph%x%label = 's [m]'

crv => grph%curve(1)
crv%name         = 'c'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = blue$
crv%line%width   = 2
crv%symbol%color = crv%line%color

grph => plt%graph(2)
grph%name                 = 'g2'
grph%type                 = 'data'
grph%margin               = qp_rect_struct(0.15, 0.06, 0.12, 0.12, '%BOX')
grph%scale_margin         = qp_rect_struct(0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, '%GRAPH')
grph%box                  = [1, 2, 1, 2]
grph%x                    = init_axis
grph%y                    = init_axis
grph%y%label_offset       = 0.15
grph%y%major_div_nominal  = 4
grph%y2                   = init_axis
grph%y2%major_div_nominal = 4
grph%y2%draw_numbers      = .false.
grph%draw_axes            = .true.
grph%draw_grid            = .true.
grph%component            = 'model'

crv => grph%curve(1)
crv%name         = 'c'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = blue$
crv%line%width   = 2
crv%symbol%color = crv%line%color

!---------------
! This plot defines the default 1-graph, 4-curve/graph plot

plt => default_plot_g1_c4

nullify(plt%r)
if (allocated(plt%graph)) deallocate (plt%graph)
allocate (plt%graph(1))
allocate (plt%graph(1)%curve(4))

plt%x_axis_type          = 's'
plt%x                    = init_axis
plt%x%major_div_nominal  = 8
plt%x%minor_div_max = 6

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
grph%component            = 'model'
grph%draw_curve_legend    = .true.
grph%draw_axes            = .true.
grph%draw_grid            = .true.
grph%text_legend_origin   = default_graph%text_legend_origin
grph%curve_legend_origin  = default_graph%curve_legend_origin
grph%x                    = init_axis
grph%x%label = 's [m]'

crv => grph%curve(1)
crv%name         = 'c1'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = blue$
crv%line%width   = 2
crv%symbol%color = crv%line%color

crv => grph%curve(2)
crv%name         = 'c2'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = orange$
crv%line%width   = 2
crv%symbol%color = crv%line%color

crv => grph%curve(3)
crv%name         = 'c3'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = green$
crv%line%width   = 2
crv%symbol%color = crv%line%color

crv => grph%curve(4)
crv%name         = 'c4'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = magenta$
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
plt%x                    = init_axis
plt%x%major_div_nominal  = 8
plt%x%minor_div_max = 6

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
grph%component            = 'model'
grph%draw_curve_legend    = .true.
grph%draw_axes            = .true.
grph%draw_grid            = .true.
grph%text_legend_origin   = default_graph%text_legend_origin
grph%curve_legend_origin  = default_graph%curve_legend_origin
grph%x                    = init_axis
grph%x%label = 's [m]'

crv => grph%curve(1)
crv%name         = 'c1'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = blue$
crv%line%width   = 2
crv%symbol%color = crv%line%color

crv => grph%curve(2)
crv%name         = 'c2'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = orange$
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
plt%x                    = init_axis

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
grph%component            = 'model'
grph%x                    = init_axis
grph%x%label              = 's [m]'

crv => grph%curve(1)
crv%name         = 'c'
crv%data_source  = 'lat'
crv%draw_symbols = .false.
crv%line%color   = blue$
crv%line%width   = 2
crv%symbol%color = crv%line%color

!---------------
! beta plot

if (all(s%plot_page%template%name /= 'beta')) then
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(2))

  plt = default_plot_g1c2
  plt%name                 = 'beta'
  plt%description          = 'Twiss beta function'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Beta Function'
  grph%y%label             = '\gb\dA\u, \gb\dB\u [m]'


  crv => grph%curve(1)
  crv%name         = 'a'
  crv%g => grph
  crv%data_type    = 'beta.a'
  crv%legend_text  = '\gb\dA\u'

  crv => grph%curve(2)
  crv%name         = 'b'
  crv%g => grph
  crv%data_type    = 'beta.b'
  crv%legend_text  = '\gb\dB\u'
endif

!---------------
! cbar plot

if (all(s%plot_page%template%name /= 'cbar')) then
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(4))

  plt = default_plot_g1_c4
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
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(2))

  plt = default_plot_g1c2
  plt%name                 = 'dbeta'
  plt%description          = 'Chromatic normalized beta beat'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Chromatic normalized beta beat'
  grph%y%label             = '\gb\u-1\d \(2265)\gb/\(2265)\gd (1)'
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
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(2))

  plt = default_plot_g1c2
  plt%name                 = 'deta'
  plt%description          = 'Second order dispersion'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Second order dispersion'
  grph%y%label             = '\(2265)\gy/\(2265)\gd (m)'
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
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(2))

  plt = default_plot_g1c2
  plt%name                 = 'detap'
  plt%description          = 'Second order dispersion slope'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Second order dispersion slope'
  grph%y%label             = "\(2265)\gy'/\(2265)\gd (1)"
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
! dphi (chrom.dphi) plot

if (all(s%plot_page%template%name /= 'dphi')) then
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(2))

  plt = default_plot_g1c2
  plt%name                 = 'dphi'
  plt%description          = 'Chromatic phase deviation'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'Chromatic phase deviation'
  grph%y%label             = '\(2265)\gf/\(2265)\gd (1)'
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
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(1))

  plt = default_plot_g1c2
  plt%name                 = 'dynamic_aperture'
  plt%description          = 'Dynamic aperture using universe calc'
  plt%x%label = 'x (mm)'
  plt%x_axis_type = 'phase_space'

  grph => plt%graph(1)
  grph%p => plt
  grph%title               = 'dynamic aperture'
  grph%type                = 'dynamic_aperture'
  grph%y%label             = 'y (mm)'
  grph%x_axis_scale_factor = 1000
  grph%y%label_offset= .15

  crv => grph%curve(1)
  crv%name         = 'c1'
  crv%g => grph
  !crv%legend_text  = 'a'  ! Legend text is automatically generated
  crv%smooth_line_calc = .false.
  crv%y_axis_scale_factor = 1000

  crv => grph%curve(2)
  crv%name         = 'c2'
  crv%g => grph
  crv%smooth_line_calc = .false.
  crv%y_axis_scale_factor = 1000
endif

!---------------
! emittance growth

if (all(s%plot_page%template%name /= 'emittance')) then
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(2))

  plt = default_plot_g1c2
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

  crv => grph%curve(2)
  crv%name         = 'b'
  crv%g => grph
  crv%data_type    = 'rad_int.i5b_e6'
  crv%legend_text  = 'b-mode emit'
  crv%y_axis_scale_factor = 7.213927194325027E-22 !for mm-mrad
endif

!---------------
! energy

if (all(s%plot_page%template%name /= 'energy')) then
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(1))

  plt = default_plot_g1c1 
  plt%name                 = 'energy'
  plt%description          = 'Particle energy'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Total Energy'
  grph%y%label       = 'E\dTot\u [eV]'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'e_tot'
endif

!---------------
! eta plot

if (all(s%plot_page%template%name /= 'eta')) then
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(2))

  plt = default_plot_g1c2
  plt%name           = 'eta'
  plt%description    = 'Dispersion'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Dispersion'
  grph%y%label       = '\gy\dX\u, \gy\dY\u [m]'

  crv => grph%curve(1)
  crv%name         = 'x'
  crv%g => grph
  crv%data_type    = 'eta.x'
  crv%legend_text  = '\gy\dX\u'

  crv => grph%curve(2)
  crv%name         = 'y'
  crv%g => grph
  crv%data_type = 'eta.y'
  crv%legend_text  = '\gy\dY\u'
endif

!---------------
! Floor Plan plot

if (all(s%plot_page%template%name /= 'floor_plan')) then
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))

  plt%name               = 'floor_plan'
  plt%description        = 'Floor plan drawing of lattice elements.'
  plt%x_axis_type        = 'floor'
  plt%x                  = init_axis

  grph => plt%graph(1)
  grph%p => plt
  grph%name                  = 'g1'
  grph%box                   = [1, 1, 1, 1]
  grph%type                  = 'floor_plan'
  grph%margin                = qp_rect_struct(0.15, 0.06, 0.05, 0.05, '%BOX')
  grph%scale_margin          = qp_rect_struct(0.02_rp, 0.02_rp, 0.02_rp, 0.02_rp, '%GRAPH')
  grph%correct_xy_distortion = .true.
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
      grph%floor_plan_view = 'yz'
    elseif (dy < min(dx, dz)) then
      grph%floor_plan_view = 'zx'
    else
      grph%floor_plan_view = 'yx'
    endif
  endif

endif

!---------------
! i1

if (all(s%plot_page%template%name /= 'i1')) then
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(1))

  plt = default_plot_g1c1 
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
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(1))

  plt = default_plot_g1c1 
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

  !---------------
  ! i3

  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(1))

  plt = default_plot_g1c1 
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
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(1))

  plt = default_plot_g1c1 
  plt%name                 = 'i4a'
  plt%description          = 'Integrated I4A Radiation integral'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Integrated I4A Radiation Integral'
  grph%y%label       = 'Integrated I4A'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'rad_int.i4a'
endif

!---------------
! i4b

if (all(s%plot_page%template%name /= 'i4b')) then
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(1))

  plt = default_plot_g1c1 
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
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(1))

  plt = default_plot_g1c1 
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
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(1))

  plt = default_plot_g1c1 
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
! Lat Layout plot

if (all(s%plot_page%template%name /= 'lat_layout')) then
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))

  plt%name           = 'lat_layout'
  plt%description    = 'Lattice elements drawn as a function of S'
  plt%x_axis_type    = 's'
  plt%x              = init_axis

  grph => plt%graph(1)
  grph%p => plt
  grph%name          = 'g1'
  grph%box           = [1, 1, 1, 1]
  grph%type          = 'lat_layout'
  grph%margin        =  qp_rect_struct(0.15, 0.06, 0.12, 0.03, '%BOX')
  grph%x             = init_axis
  grph%y%min         = -1
  grph%y%max         =  1
endif

!---------------
! Orbit plot

if (all(s%plot_page%template%name /= 'orbit')) then
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(2))

  plt = default_plot_g1c2
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
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(1))

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
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(1))

  plt = default_plot_g1c1
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
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(1))

  plt = default_plot_g1c1 
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
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(2))

  plt = default_plot_g1c2
  plt%name                 = 'phase'
  plt%description          = 'Betatron phase'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Betatron Phase'
  grph%y%label       = '\gf\dA\u, \gf\dB\u (deg)'
  grph%component     = 'model - design'

  crv => grph%curve(1)
  crv%name         = 'a'
  crv%g => grph
  crv%data_type    = 'phase.a'
  crv%legend_text  = '\gf\dA\u'

  crv => grph%curve(2)
  crv%name         = 'b'
  crv%g => grph
  crv%data_type    = 'phase.b'
  crv%legend_text  = '\gf\dB\u'
endif

!---------------
! pz

if (all(s%plot_page%template%name /= 'pz')) then
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(1))

  plt = default_plot_g1c1 
  plt%name                 = 'pz'
  plt%description          = 'Particle Pz momentum deviation'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Particle Momentum deviation Delta_P / P0'
  grph%y%label       = 'Pz'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'orbit.pz'
endif

!---------------
! sr energy loss

if (all(s%plot_page%template%name /= 'sr_energy_loss')) then
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(1))

  plt = default_plot_g1c1 
  plt%name                 = 'sr_energy_loss'
  plt%description          = 'Synch Radiation energy loss'

  grph => plt%graph(1)
  grph%p => plt
  grph%title         = 'Synch Rad Energy Loss [Rad_Integral:I2_E4 * r_e * mc^2 * 2 / 3]'
  grph%y%label       = 'E\dLoss\u [MeV * m]'

  crv => grph%curve(1)
  crv%g => grph
  crv%data_type     = 'rad_int.i2_e4'
  crv%y_axis_scale_factor = 9.59976e-16 ! (2/3) * r_e *mec2 in MeV*m
endif

!---------------
! time

if (all(s%plot_page%template%name /= 'time')) then
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(1))

  plt = default_plot_g1c1 
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
! z

if (all(s%plot_page%template%name /= 'z')) then
  np = np + 1
  plt => s%plot_page%template(np)

  nullify(plt%r)
  if (allocated(plt%graph)) deallocate (plt%graph)
  allocate (plt%graph(1))
  allocate (plt%graph(1)%curve(1))

  plt = default_plot_g1c1 
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
plt%name = 'scratch'

!-------------------------------------------------
! Regions

if (.not. allocated(s%plot_page%region)) then
  allocate (s%plot_page%region(20))
  nr = 0
else
  nr = size(s%plot_page%region)
  call move_alloc (s%plot_page%region, temp_region)
  allocate (s%plot_page%region(nr + 20))
  s%plot_page%region(1:nr) = temp_region
  deallocate(temp_region)
endif

y_layout = 0.15

if (all(s%plot_page%region(:)%name /= 'layout')) then
  nr = nr + 1
  s%plot_page%region(nr)%name = 'layout'
  s%plot_page%region(nr)%location = [0.0_rp, 1.0_rp, 0.0_rp, y_layout]
endif

do i = 1, 4
  do j = 1, i
    write (name, '(a, 2i0)') 'r', j, i
    if (any(s%plot_page%region(:)%name == name)) cycle
    nr = nr + 1
    s%plot_page%region(nr)%name = name
    y1 = y_layout + (1 - y_layout) * real(i-j)/ i
    y2 = y_layout + (1 - y_layout) * real(i-j+1) / i
    s%plot_page%region(nr)%location = [0.0_rp, 1.0_rp, y1, y2]
  enddo
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
    call tao_place_cmd ('r23', 'eta')
    call tao_place_cmd ('r33', 'orbit')
    call tao_place_cmd ('layout', 'lat_layout')
  endif
endif

end subroutine tao_setup_default_plotting

end subroutine tao_init_plotting
