!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_turn_on_special_calcs_if_needed_for_plotting ()
! 
! Routine to set u%dynch_rad_int_clac, etc to True if needed for a plot.
!-

subroutine tao_turn_on_special_calcs_if_needed_for_plotting ()

use tao_interface, dummy => tao_turn_on_special_calcs_if_needed_for_plotting

implicit none

type (tao_universe_struct), pointer :: u
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve

integer i, j, k

character(*), parameter :: r_name = 'tao_turn_on_special_calcs_if_needed_for_plotting'

! Go through all the plots and find which universes need special_calcs.
! u%picked_uni = True => need calc.

s%u(:)%picked_uni  = s%u(:)%calc%rad_int_for_plotting

s%u(:)%calc%rad_int_for_plotting    = .false.
s%u(:)%calc%chrom_for_plotting      = .false.
s%u(:)%calc%lat_sigma_for_plotting  = .false.

do i = 1, size(s%plot_page%region)
  if (.not. s%plot_page%region(i)%visible) cycle
  if (.not. allocated(s%plot_page%region(i)%plot%graph)) cycle

  do j = 1, size(s%plot_page%region(i)%plot%graph)
    graph => s%plot_page%region(i)%plot%graph(j)
    if (.not. allocated(graph%curve)) cycle

    do k = 1, size(graph%curve)
      curve => graph%curve(k)
      u => tao_pointer_to_universe(tao_curve_ix_uni(curve))
      if (.not. associated(u)) cycle

      if (.not. u%picked_uni .and. tao_rad_int_calc_needed(curve%data_type, curve%data_source)) then
        u%calc%rad_int_for_plotting = .true.
        u%calc%lattice = .true.
        u%picked_uni = .true.
      endif

      if (tao_chrom_calc_needed(curve%data_type, curve%data_source)) then
        u%calc%chrom_for_plotting = .true.
        u%calc%lattice = .true.
      endif

      if (tao_lat_sigma_calc_needed(curve%data_type, curve%data_source)) then
        u%calc%lat_sigma_for_plotting = .true.
        u%calc%lattice = .true.
      endif

    enddo     ! graph%curve
  enddo     ! plot%graph
enddo     ! plotting%region

end subroutine tao_turn_on_special_calcs_if_needed_for_plotting
