!+
! Subroutine tao_plot_setup ()
!
! Subroutine to set the data for plotting.
! Essentially transfer info from the s%u(:)%data arrays
! to the s%plot_page%region(:)%plot%graph(:)%curve(:) arrays.
!
! Input/Output:
!-

subroutine tao_plot_setup ()

use quick_plot
use tao_interface
use tao_graph_setup_mod
use tao_scale_mod
use tao_x_scale_mod
use tao_wave_mod

implicit none

type (lat_struct), pointer :: lat
type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve
type (ele_struct), pointer :: ele
type (tao_universe_struct), pointer :: u

integer i, ii, k, m, i_uni, ie, jj, ic
integer ix, ir, jg

real(rp) ax_min, ax_max, slop

integer, allocatable, save :: ix_ele(:)
logical err

character(20) :: r_name = 'tao_plot_setup'

! setup the plots

if (.not. s%global%plot_on  .and. .not. s%global%external_plotting) return

plot_loop: do ir = 1, size(s%plot_page%region)

  plot => s%plot_page%region(ir)%plot

  ! Don't worry about invisible graphs
  if (.not. s%plot_page%region(ir)%visible) cycle  

  plot%graph%is_valid = .true.
  plot%graph%why_invalid = ''

  select case (plot%x_axis_type)
  case ('index', 's', 'ele_index', 'phase_space', 'data', 'none', 'floor', 'lat',  'var', 'histogram')
  case default
    call out_io (s_error$, r_name, &
                    'BAD X_AXIS_TYPE: ' // plot%x_axis_type, &
                    'FOR PLOT: ' // plot%name)
    plot%graph%is_valid = .false.
    plot%graph%why_invalid = 'BAD X_AXIS_TYPE: ' // plot%x_axis_type
    cycle
  endselect

  ! Scale the x-axis if needed

  if ((plot%x%major_div < 0 .or. plot%x%min == plot%x%max) .and. plot%autoscale_gang_x .and. &
            size(plot%graph) > 1) call tao_x_scale_plot (plot, plot%x%min, plot%x%max, include_wall = .true.)

  if (plot%type == 'wave') then
    call tao_wave_analysis(plot)
    cycle
  endif

  ! Loop over all graphs

  do jg = 1, size(plot%graph)

    graph => plot%graph(jg)

    ! May need to x_scale after calling tao_graph_setup if scaling uses data.

    call tao_x_scale_graph (graph, graph%x%min, graph%x%max)
    call tao_graph_setup (plot, graph)
    call tao_x_scale_graph (graph, graph%x%min, graph%x%max)

    ! Scale the y-axis and determine if any points are out-of-bounds.

    call tao_scale_graph (graph, graph%y%min, graph%y%max, 'y')
    call tao_scale_graph (graph, graph%y2%min, graph%y2%max, 'y2')
    graph%limited = .false.
    if (allocated(graph%curve)) then
      do ic = 1, size(graph%curve)
        curve => graph%curve(ic)
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

enddo plot_loop

end subroutine
