!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_plot_data_setup ()
!
! Subroutine to set the data for plotting.
! Essentially transfer info from the s%u(:)%data arrays
! to the s%plot_region(:)%plot%graph(:)%curve(:) arrays.
!
! Input/Output:
!-

subroutine tao_plot_data_setup ()

use tao_x_scale_mod
use tao_scale_mod
use tao_plot_data_mod

implicit none

type (lat_struct), pointer :: lat
type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve
type (ele_struct), pointer :: ele

integer i, ii, k, m, i_uni, ie, jj, ic
integer ix, ir, jg

real(rp) ax_min, ax_max, slop

character(20) :: r_name = 'tao_plot_data_setup'

! Find which elements are to be drawn for a lattice layout and floor plan plots.
! base%lat%ele%ix_pointer will be used for storing the shape index.
! model%lat%ele%ix_pointer will be used to store the element end points.

do i_uni = 1, size(s%u)
  lat => s%u(i_uni)%model%lat
  s%u(i_uni)%base%lat%ele(:)%ix_pointer = 0
  lat%ele(:)%ix_pointer = 0
  do ie = 1, lat%n_ele_max
    ele => lat%ele(ie)
    if (ele%control_type == group_lord$) cycle
    if (ele%control_type == multipass_lord$) cycle
    if (ele%control_type == overlay_lord$) cycle
    if (ele%control_type == super_slave$) cycle
    if (ie > lat%n_ele_track .and. ele%control_type == free$) cycle
    do k = 1, size(s%plot_page%ele_shape(:))
      if (s%plot_page%ele_shape(k)%key == 0) cycle
      if (ele%key == s%plot_page%ele_shape(k)%key .and. &
               match_wild(ele%name, s%plot_page%ele_shape(k)%ele_name)) then
        ele%ix_pointer = k
        s%u(i_uni)%base%lat%ele(ie)%ix_pointer = ie
        s%u(i_uni)%base%lat%ele(ie-1)%ix_pointer = ie-1
        exit
      endif
    enddo
  enddo
enddo

! setup the plots

plot_loop: do ir = 1, size(s%plot_region)

  plot => s%plot_region(ir)%plot

  ! Don't worry about invisable graphs
  if (.not. s%plot_region(ir)%visible) cycle  

  select case (plot%x_axis_type)
  case ('index', 's', 'ele_index', 'phase_space', 'none')
  case default
    call out_io (s_abort$, r_name, 'BAD X_AXIS_TYPE: ' // plot%x_axis_type)
    plot%graph%valid = .false.
    cycle
  endselect

! loop over all graphs and curves

  do jg = 1, size(plot%graph)
    graph => plot%graph(jg)
    call tao_graph_data_setup (plot, graph)
    if (graph%y%min == graph%y%max) call tao_scale_graph (graph, '', 0.0_rp, 0.0_rp)
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

  if (plot%x%min == plot%x%max) call tao_x_scale_plot (plot, 0.0_rp, 0.0_rp)

enddo plot_loop

end subroutine
