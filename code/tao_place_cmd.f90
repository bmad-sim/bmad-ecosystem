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

implicit none

type (tao_plot_struct), pointer :: plot
type (tao_plot_struct), pointer :: template
type (qp_axis_struct), pointer :: ax

integer i
logical err

character(*) who, where
character(20) :: r_name = 'tao_place_cmd'

! Find the region where the plot is to be placed.
! The plot pointer will point to the plot associated with the region.

call tao_find_plot (err, s%plot_page%plot, 'BY_REGION', where, plot)
if (err) return

! If who = 'non' then no plot is wanted here so just turn off
! plotting in the region

if (who == 'none') then
  plot%visible = .false.
  return
endif

! Find the template for the type of plot.

call tao_find_plot (err, s%template_plot, 'BY_TYPE', who, template)
if (err) return

! transfer the plotting information from the template to the plot 
! representing the region

call tao_plot_struct_transfer (template, plot, .true.)
plot%visible = .true.

! auto scale and calculate places

if (plot%x%min == plot%x%max) then
  call tao_x_scale_cmd (where, 0.0_rp, 0.0_rp, err)
else
  ax => plot%x
  call qp_calc_axis_places (ax%min, ax%max, ax%major_div, ax%places)
endif

if (associated (plot%graph)) then
  do i = 1, size (plot%graph)
    ax => plot%graph(i)%y
    call qp_calc_axis_places (ax%min, ax%max, ax%major_div, ax%places)
    ax => plot%graph(i)%y2
    call qp_calc_axis_places (ax%min, ax%max, ax%major_div, ax%places)
  enddo
endif

end subroutine
