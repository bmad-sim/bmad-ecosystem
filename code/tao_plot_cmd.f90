!+
! Subroutine tao_plot_cmd (where, who)
!
! Routine to set what is plotted Eg: model - design, etc.
!
! Input:
!   where  -- Character(*): Region name to identify the plot to set.
!   who(:) -- character(*): Who to plot. First character must be a "+"
!               or a "-" to designate which goes in the baseline.
!
!  Output:
!-

subroutine tao_plot_cmd (where, who)

use tao_mod

implicit none

type (tao_plot_array_struct), allocatable, save :: plot(:)
type (tao_graph_array_struct), allocatable, save :: graph(:)
type (tao_plot_region_struct), pointer :: region
type (tao_plot_who_struct) :: p_who(n_who_maxx)
integer i, j
integer ix, ix_line, ix_cmd, which, n_word

character(*) :: where
character(*) :: who(:)
character(20) :: r_name = 'tao_plot_cmd'

logical err

! Check the who.

err = .true.
p_who%name = ' '
p_who%sign = 0

do i = 1, size(who)

  call string_trim (who(i)(2:), p_who(i)%name, j)

  if (.not. any(p_who(i)%name == s%global%valid_plot_who)) then
    call out_io (s_error$, r_name, 'BAD "WHO": ' // who(i))
    return
  endif

  select case (who(i)(1:1))
  case ('+')
    p_who(i)%sign = +1
  case ('-')
    p_who(i)%sign = -1
  case default
    call out_io (s_error$, r_name, 'NO +/- SIGN: ' // who(i))
    return
  end select

enddo

! Find plot for the region given by "where"

if (where == 'all') then
  do i = 1, size(s%plot_page%region)
    do j = 1, size(s%plot_page%region(i)%plot%graph)
      s%plot_page%region(i)%plot%graph(j)%who = p_who
    enddo
  enddo

else
  call tao_find_plots (err, where, 'REGION', plot, graph)
  if (err) return
  if (allocated(graph)) then
    do i = 1, size(graph)
      graph(i)%g%who = p_who
    enddo
  else
    do i = 1, size(plot)
      do j = 1, size(plot(i)%p%graph)
        plot(i)%p%graph(j)%who = p_who
      enddo
    enddo
  endif
endif

err = .false.

end subroutine 




