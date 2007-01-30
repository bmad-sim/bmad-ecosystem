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

type (tao_plot_array_struct), allocatable, save :: template(:)
type (qp_axis_struct), pointer :: ax
type (tao_plot_region_struct), pointer :: region

integer i, j, k
logical err

character(*) who, where
character(20) :: r_name = 'tao_place_cmd'

! Find the region where the plot is to be placed.
! The plot pointer will point to the plot associated with the region.

call tao_find_plot_region (err, where, region)
if (err) return

! free up any s%beam_save structures.

do i = 1, size(region%plot%graph)
  if (region%plot%graph(i)%type /= 'phase_space') cycle
  do j = 1, size(region%plot%graph(i)%curve)
    region%plot%graph(i)%curve(j)%beam_save%ix_universe = -1  ! free-up
  enddo
enddo

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
  s%global%lattice_recalc = .true.
  curve_loop: do j = 1, size(region%plot%graph(i)%curve)
    do k = 1, size(s%beam_save)
      if (s%beam_save(k)%ix_universe > -1) cycle ! look for unused
      region%plot%graph(i)%curve(j)%beam_save => s%beam_save(k)
      s%beam_save(k)%ix_universe = region%plot%graph(i)%curve(j)%ix_universe
      s%beam_save(k)%ix_ele = region%plot%graph(i)%curve(j)%ix_ele_ref
      cycle curve_loop
    enddo
    print *, 'ERROR IN TAO_PLACE_CMD: NO AVAILABLE BEAM_SAVE SLOTS!'
    call err_exit
  enddo curve_loop
enddo

end subroutine
