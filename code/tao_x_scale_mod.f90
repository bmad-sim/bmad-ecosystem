module tao_x_scale_mod

use tao_mod
use quick_plot

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_x_scale_cmd (where, x_min, x_max, err)
!
! Routine to scale a plot. If x_min = x_max
! Then the scales will be chosen to show all the data.
! 
! Input:
!   where -- Character(*): Region to scale. Eg: "top"
!   x_min -- Real(rp): Plot x-axis min value.
!   x_max -- Real(rp): Plot x-axis max value.
!
!  Output:
!   err -- Logical: Set to True if the plot cannot be found. False otherwise.
!-

subroutine tao_x_scale_cmd (where, x_min, x_max, err)

implicit none

type (tao_plot_struct), pointer :: p, p2
type (tao_plot_array_struct), allocatable, save :: plot(:)

real(rp) x_min, x_max

integer i, j, ix, places, im

character(*) where

real(rp) this_min, this_max

logical err

! If the where argument is blank then scale all graphs.
! Make sure that all plots using an 's' scale have the same scale

if (len_trim(where) == 0 .or. where == 'all') then

  im = 0
  do j = 1, size(s%plot_page%region)
    p => s%plot_page%region(j)%plot
    if (.not. s%plot_page%region(j)%visible) cycle
    call tao_x_scale_plot (p, x_min, x_max)
    if (p%x_axis_type == 's' .and. &
              any (p%graph(:)%type /= 'lat_layout')) then
      if (im == 0) im = j
      if (p%x%max < s%plot_page%region(im)%plot%x%max) im = j
    endif
  enddo

  if (im /= 0) then
    do j = 1, size(s%plot_page%region)
      p => s%plot_page%region(j)%plot
      if (.not. s%plot_page%region(j)%visible) cycle
      if (p%x_axis_type /= 's') cycle
      p2 => s%plot_page%region(im)%plot
      p%x%max = p2%x%max
      p%x%min = p2%x%min
      p%x%major_div = p2%x%major_div
    enddo
  endif

  return
endif

! locate the plot by the region name given by the where argument.
! If where has a '.' then we are dealing with just one graph of the plot.
! Otherwise we scale all the graphs of the plot.

call tao_find_plots (err, where, 'REGION', plot)
if (err) return
do i = 1, size(plot)
  call tao_x_scale_plot (plot(i)%p, x_min, x_max)
enddo

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_x_scale_plot (plot, x_min, x_max)
!
! Routine to scale a plot. If x_min = x_max
! Then the scales will be chosen to show all the data.
! 
! Input:
!   plot  -- Tao_plot_struct: Plot to scale. Eg: "top"
!   x_min -- Real(rp): Plot x-axis min value.
!   x_max -- Real(rp): Plot x-axis max value.
!-

subroutine tao_x_scale_plot (plot, x_min, x_max)

type (tao_plot_struct), target :: plot
type (tao_curve_struct), pointer :: curve
integer j, k, n, p1, p2, iu, ix
real(rp) x_min, x_max
real(rp) x1, x2
logical curve_here

! Check if the thing exists

if (.not. associated (plot%graph)) return

! auto scale

if (x_max == x_min) then

  plot%x%min = -1e20     ! to setup all data
  plot%x%max = 1e20

  x1 =  1e20
  x2 = -1e20
  curve_here = .false.
  do j = 1, size(plot%graph)
    if (plot%graph(j)%type == 'key_table') cycle
    if (plot%graph(j)%type == 'lat_layout') then
      iu = plot%graph(j)%ix_universe
      if (iu == 0) iu = s%global%u_view
      x1 = min (x1, s%u(iu)%model%lat%ele_(0)%s)
      ix = s%u(iu)%model%lat%n_ele_use
      x2 = max (x2, s%u(iu)%model%lat%ele_(ix)%s)
      curve_here = .true.
    else
      if (.not. associated(plot%graph(j)%curve)) cycle
      do k = 1, size(plot%graph(j)%curve)
        curve => plot%graph(j)%curve(k)
        if (associated (curve%x_symb)) then
          curve_here = .true.
          x1 = min (x1, minval(curve%x_symb(:)))
          x2 = max (x2, maxval(curve%x_symb(:)))
        endif
        if (associated (curve%x_line)) then
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

end subroutine 

end module
