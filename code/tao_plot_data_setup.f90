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
type (tao_universe_struct), pointer :: u
type (tao_ele_shape_struct), pointer :: shape(:)

integer i, ii, k, m, i_uni, ie, jj, ic
integer ix, ir, jg

real(rp) ax_min, ax_max, slop

integer, allocatable, save :: ix_ele(:)
logical err

character(20) :: r_name = 'tao_plot_data_setup'

! Find which elements are to be drawn for a lattice layout and floor plan plots.

if (.not. s%global%plot_on) return

do i_uni = lbound(s%u, 1), ubound(s%u, 1)

  u => s%u(i_uni)
  if (.not. associated(u%ele)) cycle

  lat => u%model%lat

  u%ele%ix_shape_lat_layout = 0
  u%ele%ix_ele_end_lat_layout = 0
  u%ele%ix_shape_floor_plan = 0
  u%ele%ix_ele_end_floor_plan = 0

  ! floor_plan

  shape => tao_com%ele_shape_floor_plan
  do k = 1, size(shape)
    if (shape(k)%ele_name == '') cycle
    call tao_ele_locations_given_name (u, shape(k)%ele_name, ix_ele, err, .false.)
    if (err) then
      call out_io (s_error$, r_name, 'BAD ELEMENT KEY IN SHAPE: ' // shape(k)%ele_name)
      cycle
    endif

    do i = 1, size(ix_ele)
      ie = ix_ele(i)
      ele => lat%ele(ie)
      if (ele%control_type == group_lord$) cycle
      if (ele%control_type == multipass_lord$) cycle
      if (ele%control_type == overlay_lord$) cycle
      if (ele%control_type == super_slave$) cycle
      if (ie > lat%n_ele_track .and. ele%control_type == free$) cycle

      u%ele(ie)%ix_shape_floor_plan = k
      u%ele(ie)%ix_ele_end_floor_plan = ie
      u%ele(ie-1)%ix_ele_end_floor_plan = ie-1
    enddo
  enddo

! lat_layout

  shape => tao_com%ele_shape_lat_layout
  do k = 1, size(shape)
    if (shape(k)%ele_name == '') cycle
    call tao_ele_locations_given_name (u, shape(k)%ele_name, ix_ele, err, .false.)
    if (err) then
      call out_io (s_error$, r_name, 'BAD ELEMENT KEY IN SHAPE: ' // shape(k)%ele_name)
      cycle
    endif

    do i = 1, size(ix_ele)
      ie = ix_ele(i)
      ele => lat%ele(ie)
      if (ele%control_type == group_lord$) cycle
      if (ele%control_type == multipass_lord$) cycle
      if (ele%control_type == overlay_lord$) cycle
      if (ele%control_type == super_slave$) cycle
      if (ie > lat%n_ele_track .and. ele%control_type == free$) cycle

      u%ele(ie)%ix_shape_lat_layout = k
      u%ele(ie)%ix_ele_end_lat_layout = ie
      u%ele(ie-1)%ix_ele_end_lat_layout = ie-1
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

  if (plot%x%min == plot%x%max) call tao_x_scale_plot (plot, 0.0_rp, 0.0_rp)

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

enddo plot_loop

end subroutine
