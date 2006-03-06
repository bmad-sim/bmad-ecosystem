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

type (tao_plot_struct), pointer :: template
type (qp_axis_struct), pointer :: ax
type (tao_plot_region_struct), pointer :: region

integer i
logical err

character(*) who, where
character(20) :: r_name = 'tao_place_cmd'

! Find the region where the plot is to be placed.
! The plot pointer will point to the plot associated with the region.

call tao_find_plot_by_region (err, where, region = region)
if (err) return

! If who = 'non' then no plot is wanted here so just turn off
! plotting in the region

if (who == 'none') then
  region%visible = .false.
  return
endif

! Find the template for the type of plot.

call tao_find_template_plot (err, who, template)
if (err) return

! transfer the plotting information from the template to the plot 
! representing the region

call tao_plot_struct_transfer (template, region%plot)
region%visible = .true.

! If the plot has a phase_space curve then recalculate the lattice

do i = 1, size(template%graph)
  if (template%graph(i)%type == 'phase_space') s%global%lattice_recalc = .true.
enddo

! auto scale and calculate places

if (region%plot%x%min == region%plot%x%max) then
  call tao_x_scale_cmd (where, 0.0_rp, 0.0_rp, err)
else
  ax => region%plot%x
  call qp_calc_axis_places (ax)
endif

if (associated (region%plot%graph)) then
  do i = 1, size (region%plot%graph)
    ax => region%plot%graph(i)%y
    call qp_calc_axis_places (ax)
    ax => region%plot%graph(i)%y2
    call qp_calc_axis_places (ax)
  enddo
endif

end subroutine
