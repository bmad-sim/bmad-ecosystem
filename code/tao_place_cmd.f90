!+
! Subroutine tao_place_cmd (where, who)
!
! Subroutine to determine the placement of a plot in the plot window.
! The appropriate s%tamplate_plot(i) determined by the who argument is
! transfered to the appropriate s%plot_page%plot(j) determined by the where
! argument.
!
! Input:
!    %template_plot(i) -- template matched to who.
!   where -- Character(*): Region where the plot goes. Eg: 'top'.
!   who   -- Character(*): Type of plot. Eg: 'orbit'.
!
! Output
!    %plot_page%plot(j) -- Plot matched to where.
!-

subroutine tao_place_cmd (where, who)

use tao_mod
use tao_x_scale_mod
use quick_plot
use beam_mod

implicit none

type (tao_plot_array_struct), allocatable, save :: template(:)
type (qp_axis_struct), pointer :: ax
type (tao_universe_struct), pointer :: u
type (tao_plot_region_struct), pointer :: region
type (tao_curve_struct), pointer :: curve

integer i, j, k, i_uni
logical err

character(*) who, where
character(20) :: r_name = 'tao_place_cmd'

! Find the region where the plot is to be placed.
! The plot pointer will point to the plot associated with the region.

call tao_find_plot_region (err, where, region)
if (err) return

! If who = 'none' then no plot is wanted here so just turn off
! plotting in the region

if (who == 'none') then
  region%visible = .false.
  return
endif

! Find the template for the type of plot.

call tao_find_plots (err, who, 'TEMPLATE', template)
if (err) return

! transfer the plotting information from the template to the plot 
! representing the region

call tao_plot_struct_transfer (template(1)%p, region%plot)
region%visible = .true.
region%plot%r => region

! If the plot has a phase_space curve then recalculate the lattice

do i = 1, size(region%plot%graph)
  if (region%plot%graph(i)%type /= 'phase_space') cycle
  do j = 1, size(region%plot%graph(i)%curve)
    curve => region%plot%graph(i)%curve(j)
    if (curve%ix_ele_ref_track < 0) then
      call out_io (s_error$, r_name, &
                'BAD REFERENCE ELEMENT: ' // curve%ele_ref_name, &
                'CANNOT PLOT PHASE SPACE FOR: ' // tao_curve_name(curve))
      return
    endif
    u => s%u(tao_universe_number(curve%ix_universe))
    if (.not. allocated(u%ele(curve%ix_ele_ref_track)%beam%bunch)) then
      tao_com%lattice_recalc = .true.
      u%ele(curve%ix_ele_ref_track)%save_beam = .true.
    endif
  enddo
enddo

! If the plot uses 's' for the x-axis and min/max have not been set then use the scale for 
! an existing 's' plot.

if (region%plot%x_axis_type == 's' .and. region%plot%x%min == region%plot%x%max) then
  do i = 1, size(s%plot_region)
    if (.not. s%plot_region(i)%visible) cycle
    if (s%plot_region(i)%plot%x_axis_type /= 's') cycle
    region%plot%x = s%plot_region(i)%plot%x
    exit
  enddo
endif

end subroutine
