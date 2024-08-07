!+
! Subroutine tao_place_cmd (where, who, no_buffer)
!
! Subroutine to determine the placement of a plot in the plot window.
! The appropriate s%tamplate_plot(i) determined by the who argument is
! transfered to the appropriate s%plot_page%plot(j) determined by the where
! argument.
!
! Input:
!   s%plot_page%template(i) -- template matched to who.
!   where       -- Character(*): Region where the plot goes. Eg: 'top'.
!   who         -- Character(*): Type of plot. Eg: 'orbit'.
!   no_buffer   -- logical, optional: If present and True then prevents buffering 
!                   in the case when s%global%external_plotting = T
!
! Output
!   s%plot_page%plot(j) -- Plot matched to where.
!-

subroutine tao_place_cmd (where, who, no_buffer)

use tao_x_scale_mod, dummy => tao_place_cmd
use beam_mod

implicit none

type (tao_plot_array_struct), allocatable :: template(:)
type (qp_axis_struct), pointer :: ax
type (tao_universe_struct), pointer :: u
type (tao_plot_region_struct), pointer :: region, r2
type (tao_curve_struct), pointer :: curve
type (tao_plot_region_struct), allocatable :: temp_buf(:)

real(rp) x1, x2, y1, y2
integer i, j, k, i_uni, ib, n

logical err
logical, optional :: no_buffer

character(*) who, where
character(*), parameter :: r_name = 'tao_place_cmd'

! If external plotting is being done then just save the arguments and the GUI
! can query for the info via a pipe command.

if (s%global%external_plotting .and. .not. logic_option(.false., no_buffer)) then
  if (allocated (s%com%plot_place_buffer)) then
    n = size(s%com%plot_place_buffer) + 1
    call move_alloc(s%com%plot_place_buffer, temp_buf)
    allocate(s%com%plot_place_buffer(n))
    s%com%plot_place_buffer(1:n-1) = temp_buf
  else
    allocate(s%com%plot_place_buffer(1))
    n = 1
  endif

  s%com%plot_place_buffer(n)%name = where
  s%com%plot_place_buffer(n)%plot%name = who

  return
endif


! Find the region where the plot is to be placed.
! The plot pointer will point to the plot associated with the region.

if (where == '*' .and. who == 'none') then
  do i = 1, size(s%plot_page%region)
    s%plot_page%region(i)%visible = .false.
    s%plot_page%region(i)%plot%name = ''
    s%plot_page%region(i)%plot%r => null()
  enddo
  return
endif

call tao_find_plot_region (err, where, region)
if (err) return

! If who = 'none' then no plot is wanted here so just turn off plotting in the region

if (who == 'none') then
  region%visible = .false.
  region%plot%name = ''
  region%plot%r => null()
  return
endif

! Find the template for the type of plot.

call tao_find_plots (err, who, 'TEMPLATE', template)
if (err) return
if (size(template) > 1) then
  call out_io (s_error$, r_name, 'MULTIPLE TEMPLATE PLOTS MATCH: ' // who, &
                                 'USE THE COMMAND "show plot -templates" TO SEE A LIST OF PLOTS.')
  return
endif

! Transfer the plotting information from the template to the plot 
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
    u => s%u(tao_universe_index(tao_curve_ix_uni(curve)))
    ib = tao_branch_index(curve%ix_branch)
    if (.not. allocated(u%model_branch(ib)%ele(curve%ix_ele_ref_track)%beam%bunch)) then
      if (s%global%plot_on) u%calc%lattice = .true.
      u%model_branch(ib)%ele(curve%ix_ele_ref_track)%save_beam_internally = .true.
    endif
  enddo
enddo

! Turn off overlapping plots.
! Do this even if plot_on = F since in a command file might have plot_on = F.

if (s%plot_page%delete_overlapping_plots .and. .not. s%global%external_plotting) then
  do i = 1, size(s%plot_page%region)
    r2 => s%plot_page%region(i)
    if (.not. r2%visible) cycle
    if (r2%name == region%name) cycle
    if (r2%location(1) > region%location(2) - 0.02) cycle
    if (r2%location(2) < region%location(1) + 0.02) cycle
    if (r2%location(3) > region%location(4) - 0.02) cycle
    if (r2%location(4) < region%location(3) + 0.02) cycle
    r2%visible = .false. 
    r2%plot%name = ''
    r2%plot%r => null()
  enddo
endif

! Check to see if radiation integrals need be computed

call tao_turn_on_special_calcs_if_needed_for_plotting()

end subroutine
