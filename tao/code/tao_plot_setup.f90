!+
! Subroutine tao_plot_setup ()
!
! Subroutine to set the data for plotting.
!-

subroutine tao_plot_setup ()

use quick_plot
use tao_interface, dummy => tao_plot_setup
use tao_graph_setup_mod, only: tao_graph_setup
use tao_scale_mod, only: tao_scale_graph
use tao_x_scale_mod, only: tao_x_scale_graph
use tao_wave_mod, only: tao_wave_analysis

implicit none

type (tao_plot_region_struct), pointer :: region
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve
type (tao_plot_array_struct), allocatable :: plots(:)
integer ir, ig, ic, ip
logical err
character(20) :: r_name = 'tao_plot_setup'

!

if ((.not. s%global%plot_on  .and. .not. s%global%external_plotting) .or. .not. s%global%lattice_calc_on) return
s%plot_page%region%setup_done = .false.

! First do plots that are needed by parametric plots.

do ir = 1, size(s%plot_page%region)
  region => s%plot_page%region(ir)
  if (.not. region%visible) cycle  
  if (.not. allocated(region%plot%graph)) cycle
  do ig = 1, size(region%plot%graph)
    graph => region%plot%graph(ig)
    if (.not. allocated(graph%curve)) cycle
    ic_loop: do ic = 1, size(graph%curve)
      curve => graph%curve(ic)
      if (curve%data_source /= 'curve') cycle
      call tao_find_plots(err, curve%data_type, 'REGION', plots, print_flag = .false., only_visible = .false.)
      if (size(plots) == 0) then
        call setup_this_dependent_plot(curve)
      else
        do ip = 1, size(plots)
          if (plots(ip)%p%r%setup_done) cycle  ! Already used by some other parametric plot
          call setup_this_dependent_plot(curve, plots(ip)%p%r)
          cycle ic_loop 
        enddo
        call setup_this_dependent_plot (curve)
      endif
    enddo ic_loop
  enddo
enddo

! Normal plots

do ir = 1, size(s%plot_page%region)
  region => s%plot_page%region(ir)
  ! Don't worry about invisible graphs
  if (.not. region%visible) cycle  
  call setup_this_plot (region)
enddo

!---------------------------------------------------------------------
contains

subroutine setup_this_dependent_plot(curve, region)

type (tao_curve_struct) curve
type (tao_plot_region_struct), optional, target :: region
type (tao_plot_region_struct), pointer :: r
type (tao_plot_array_struct), allocatable :: plots(:)
integer ix, ir
logical err
character(40) plot_name, graph_name

!

ix = index(curve%data_type, '.')
if (ix == 0) then
  plot_name = curve%data_type
  graph_name = ''
else
  plot_name = curve%data_type(:ix-1)
  graph_name = curve%data_type(ix+1:)
endif

!

if (.not. present(region)) then 
  do ir = size(s%plot_page%region), 1, -1
    r => s%plot_page%region(ir)
    if (r%plot%name /= '') cycle
    call tao_find_plots (err, plot_name, 'TEMPLATE', plots, print_flag = .false.)
    if (size(plots) == 0) return
    call tao_plot_struct_transfer(plots(1)%p, r%plot)
    r%plot%r => r
    exit
  enddo
else
  r => region
endif

!

if (r%name /= plot_name) then
  r%name = trim(plot_name) // '_' // trim(curve%g%p%r%name)
  curve%data_type = trim(r%name)
  if (graph_name /= '') curve%data_type = trim(curve%data_type) // '.' // graph_name
endif

call setup_this_plot (r)

end subroutine setup_this_dependent_plot

!---------------------------------------------------------------------
! contains

subroutine setup_this_plot (region)

type (tao_plot_region_struct), target :: region
type (tao_plot_struct), pointer :: plot
type (lat_struct), pointer :: lat
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve
type (ele_struct), pointer :: ele
type (tao_universe_struct), pointer :: u

integer i, ii, k, m, i_uni, ie, jj, ic
integer ix, jg

real(rp) ax_min, ax_max, slop

integer, allocatable :: ix_ele(:)
logical err

!

if (region%setup_done) return
plot => region%plot

plot%graph%is_valid = .true.
plot%graph%why_invalid = ''
region%setup_done = .true.

select case (plot%x_axis_type)
case ('index', 's', 'z', 'ele_index', 'phase_space', 'data', 'none', 'floor', 'lat',  'var', 'histogram', 'curve')
case default
  call out_io (s_error$, r_name, &
                  'BAD X_AXIS_TYPE: ' // plot%x_axis_type, &
                  'FOR PLOT: ' // plot%name)
  plot%graph%is_valid = .false.
  plot%graph%why_invalid = 'BAD X_AXIS_TYPE: ' // plot%x_axis_type
  return
endselect

! Scale the x-axis if needed

if (plot%type == 'wave') then
  call tao_wave_analysis(plot)
  return
endif

! Loop over all graphs

do jg = 1, size(plot%graph)

  graph => plot%graph(jg)

  ! May need to x_scale after calling tao_graph_setup if scaling uses data.

  call tao_x_scale_graph (graph, graph%x%eval_min, graph%x%eval_max, include_wall = .true.)
  call tao_graph_setup (plot, graph)
  if (.not. graph%is_valid) cycle
  call tao_x_scale_graph (graph, graph%x%eval_min, graph%x%eval_max, include_wall = .true.)

  ! Scale the y-axis and determine if any points are out-of-bounds.

  call tao_scale_graph (graph, graph%y%min, graph%y%max, 'y', include_wall = .true.)
  call tao_scale_graph (graph, graph%y2%min, graph%y2%max, 'y2', include_wall = .true.)
  graph%limited = .false.
  if (allocated(graph%curve)) then
    do ic = 1, size(graph%curve)
      curve => graph%curve(ic)
      if (.not. curve%valid) cycle

      call qp_get_parameters (default_axis_slop_factor = slop)
      if (curve%use_y2) then
        ax_min = graph%y2%min + slop * (graph%y2%min - graph%y2%max)
        ax_max = graph%y2%max + slop * (graph%y2%max - graph%y2%min)
      else
        ax_min = graph%y%min + slop * (graph%y%min - graph%y%max)
        ax_max = graph%y%max + slop * (graph%y%max - graph%y%min)
      endif

      if (allocated(curve%y_symb)) then
        if (any(curve%y_symb < ax_min)) graph%limited = .true.
        if (any(curve%y_symb > ax_max)) graph%limited = .true.
      endif
      if (allocated(curve%y_line)) then
        if (any(curve%y_line < ax_min)) graph%limited = .true.
        if (any(curve%y_line > ax_max)) graph%limited = .true.
      endif
    enddo
  endif
enddo

end subroutine setup_this_plot

end subroutine
