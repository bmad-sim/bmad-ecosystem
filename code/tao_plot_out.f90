!+
! Subroutine tao_plot_out (plot_page)
!
! Subroutine to draw the plots on the plot window.
!
! Input:
!   plot_page  -- Tao_plot_page_struct:
!-

subroutine tao_plot_out (plot_page)

  use tao_mod
  use quick_plot

  implicit none

  type (tao_plot_page_struct), target :: plot_page
  type (tao_plot_struct), pointer :: plot
  type (tao_graph_struct), pointer :: graph
  type (tao_curve_struct), pointer :: curve
  type (qp_rect_struct) border1, border2

  real(rp) region(4), dx, dy
  integer i, j, k, graph_box(4)

! inits

  call qp_clear_page

! loop over all plots

  do i = 1, size(plot_page%plot)

    plot => plot_page%plot(i)
    if (.not. plot%visible) cycle

! set the plot_page border for this particular region

    region = plot_page%plot(i)%region%location
    border1%units = '%PAGE'
    call qp_convert_rectangle_rel (plot_page%border, border1)
    dx = 1 - (border1%x2 - border1%x1)
    dy = 1 - (border1%y2 - border1%y1)
    border2%x1 = border1%x1 + dx * region(1)
    border2%x2 = border1%x2 + dx * (1 - region(2))
    border2%y1 = border1%y1 + dy * region(3)
    border2%y2 = border1%y2 + dy * (1 - region(4))
    border2%units = '%PAGE'
    call qp_set_layout (page_border = border2)

! loop over all the graphs of the plot and draw the axes.

    do j = 1, size(plot%graph)
      graph => plot%graph(j)
      graph_box = (/ graph%this_box, plot%box_layout /)
      call qp_set_layout (box = graph_box, margin = graph%margin)
      call qp_set_layout (x_axis = plot%x, y_axis = graph%y)
      call qp_set_graph (title = trim(graph%title) // ' ' // graph%title_suffix)
      call qp_draw_axes

! loop over all the curves of the graph and draw them

      do k = 1, size(graph%curve)
        curve => graph%curve(k)
        call qp_set_symbol (curve%symbol)
        call qp_set_line ('PLOT', curve%line) 
        call qp_draw_data (curve%x_dat, curve%y_dat, curve%draw_line, &
                                               curve%symbol_every, graph%clip)
      enddo
    enddo
  enddo

end subroutine
