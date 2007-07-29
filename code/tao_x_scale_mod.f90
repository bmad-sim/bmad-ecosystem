module tao_x_scale_mod

use tao_mod
use quick_plot

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_x_scale_cmd (where, x_min, x_max, err, force)
!
! Routine to scale a plot. If x_min = x_max
! Then the scales will be chosen to show all the data.
! 
! Input:
!   where -- Character(*): Region to scale. Eg: "top"
!   x_min -- Real(rp): Plot x-axis min value.
!   x_max -- Real(rp): Plot x-axis max value.
!   force -- Logical: (Optional) If true then do not use nice x-scales but
!                     exactly what is specified
!
!  Output:
!   err -- Logical: Set to True if the plot cannot be found. False otherwise.
!-

subroutine tao_x_scale_cmd (where, x_min, x_max, err, force)

implicit none

type (tao_plot_struct), pointer :: p, p2
type (tao_plot_array_struct), allocatable, save :: plot(:)

real(rp) x_min, x_max
integer i, j, n, ix, places, im
character(*) where
logical err
logical, optional :: force

! find plots to scale

if (len_trim(where) == 0 .or. where == 'all' .or. where == 's') then
  if (allocated(plot)) deallocate(plot)
  n = 0
  do j = 1, size(s%plot_region)
    if (.not. s%plot_region(j)%visible) cycle
    if (where == 's' .and. s%plot_region(j)%plot%x_axis_type /= 's') cycle
    n = n + 1
  enddo
  allocate (plot(n))
  n = 0
  do j = 1, size(s%plot_region)
    if (.not. s%plot_region(j)%visible) cycle
    if (where == 's' .and. s%plot_region(j)%plot%x_axis_type /= 's') cycle
    n = n + 1
    plot(n)%p => s%plot_region(j)%plot
  enddo
else
  call tao_find_plots (err, where, 'REGION', plot)
  if (err) return
endif

! If autoscaling then figure out the limits of data.

if (x_min == x_max) then
  do i = 1, size(plot)
    plot(i)%p%x%min = -1e30  ! So no cliping of points
    plot(i)%p%x%max = 1e30   ! So no cliping of points
  enddo
  call tao_plot_data_setup()
endif

! Set max and min

  do i = 1, size(plot)
    if (present(force)) then
      call tao_x_scale_plot (plot(i)%p, x_min, x_max, force)
    else
      call tao_x_scale_plot (plot(i)%p, x_min, x_max)
    endif
  enddo

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_x_scale_plot (plot, x_min, x_max, force)
!
! Routine to scale a plot. If x_min = x_max
! Then the scales will be chosen to show all the data.
! 
! Input:
!   plot  -- Tao_plot_struct: Plot to scale. Eg: "top"
!   x_min -- Real(rp): Plot x-axis min value.
!   x_max -- Real(rp): Plot x-axis max value.
!   force -- Logical: (Optional) forces the use of the specified min and max
!-

subroutine tao_x_scale_plot (plot, x_min, x_max, force)

type (tao_plot_struct), target :: plot
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve
integer j, k, n, p1, p2, iu, ix
real(rp) x_min, x_max
real(rp) x1, x2
logical curve_here
logical, optional :: force

! Check if the thing exists

if (.not. allocated (plot%graph)) return

! auto scale

if (x_max == x_min) then

  plot%x%min = -1e20     ! to setup all data
  plot%x%max = 1e20

  x1 =  1e20
  x2 = -1e20
  curve_here = .false.
  do j = 1, size(plot%graph)
    graph => plot%graph(j)
    if (graph%type == 'key_table') then
      cycle
    else if (graph%type == 'floor_plan') then
      x1 = min(x1, graph%x_min)
      x2 = max(x2, graph%x_max)
      curve_here = .true.
    else if (plot%x_axis_type == 's') then
      iu = graph%ix_universe
      if (iu == 0) iu = s%global%u_view
      x1 = min (x1, s%u(iu)%model%lat%ele(0)%s)
      ix = s%u(iu)%model%lat%n_ele_track
      x2 = max (x2, s%u(iu)%model%lat%ele(ix)%s)
      curve_here = .true.
    else
      if (.not. allocated(graph%curve)) cycle
      do k = 1, size(graph%curve)
        curve => graph%curve(k)
        if (allocated (curve%x_symb)) then
          curve_here = .true.
          x1 = min (x1, minval(curve%x_symb(:)))
          x2 = max (x2, maxval(curve%x_symb(:)))
        endif
        if (allocated (curve%x_line)) then
          curve_here = .true.
          x1 = min (x1, minval(curve%x_line(:)))
          x2 = max (x2, maxval(curve%x_line(:)))
        endif
      enddo
    endif
  enddo

  if (.not. curve_here) return

else
  x1 = x_min
  x2 = x_max
endif

! calculate divisions and places

p1 = nint(0.7 * plot%x_divisions)  ! Used to be 8
p2 = nint(1.3 * plot%x_divisions)  ! Used to be 15
call qp_calc_and_set_axis ('X', x1, x2, p1, p2, 'GENERAL', plot%x%type)
call qp_get_axis ('X', plot%x%min, plot%x%max, plot%x%major_div, plot%x%places)
if (present(force)) then
  if (force) then
    plot%x%min = x_min
    plot%x%max = x_max
  endif
endif
  
end subroutine 

end module
