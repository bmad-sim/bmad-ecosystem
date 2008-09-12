!+
! Subroutine tao_plot_setup ()
!
! Subroutine to set the data for plotting.
! Essentially transfer info from the s%u(:)%data arrays
! to the s%plot_region(:)%plot%graph(:)%curve(:) arrays.
!
! Input/Output:
!-

subroutine tao_plot_setup ()

use tao_x_scale_mod
use tao_scale_mod
use tao_graph_setup_mod

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

if (.not. s%global%plot_on) return

plot_loop: do ir = 1, size(s%plot_region)

  plot => s%plot_region(ir)%plot

  ! Don't worry about invisable graphs
  if (.not. s%plot_region(ir)%visible) cycle  

  select case (plot%x_axis_type)
  case ('index', 's', 'ele_index', 'phase_space', 'data', 'none')
  case default
    call out_io (s_abort$, r_name, &
                    'BAD X_AXIS_TYPE: ' // plot%x_axis_type, &
                    'FOR PLOT: ' // plot%name)
    plot%graph%valid = .false.
    cycle
  endselect

! loop over all graphs and curves

  ! If x%min = x%max we need to autoscale but first we need to call tao_graph_setup to make sure
  ! the curves contain all the data.

  if (plot%x%min == plot%x%max .and. plot%autoscale_gang_x) then
    call tao_x_scale_plot (plot, -1.0e30_rp, 1.0e30_rp) ! To include all the data
    do jg = 1, size(plot%graph)
      call tao_graph_setup (plot, plot%graph(jg))  ! And populate the curves
    enddo
    call tao_x_scale_plot (plot, 0.0_rp, 0.0_rp)  ! Now we can autoscale

  else
    do jg = 1, size(plot%graph)
      graph => plot%graph(jg)
      if (graph%x%min == graph%x%max) then
        call tao_x_scale_graph (graph, -1.0e30_rp, 1.0e30_rp)
        call tao_graph_setup (plot, graph)
        call tao_x_scale_graph (graph, 0.0_rp, 0.0_rp)
      else
        call tao_graph_setup (plot, graph)
      endif
    enddo
  endif

  ! Now scale the y axis if needed and determine if any points are out-of-bounds.

  do jg = 1, size(plot%graph)
    graph => plot%graph(jg)
    if (graph%y%min == graph%y%max) call tao_scale_graph (graph, 0.0_rp, 0.0_rp)
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
