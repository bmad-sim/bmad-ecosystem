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

type (tao_plot_struct), pointer :: plot

integer i, j
integer ix, ix_line, ix_cmd, which, n_word

character(*) :: where
character(*) :: who(:)
character(20) :: r_name = 'tao_plot_cmd'

logical err

! Check the who.

err = .true.

do i = 1, size(who)
  if (who(i)(1:1) /= '+' .and. who(i)(1:1) /= '-') then
    call out_io (s_error$, r_name, 'NO +/- SIGN: ' // who(i))
    return
  endif

  if (.not. any(who(i)(2:) == s%global%valid_plot_who)) then
    call out_io (s_error$, r_name, 'BAD "WHO": ' // who(i)(2:))
    return
  endif

enddo

! Find plot for the region given by "where"

call tao_find_plot (err, s%plot_page%plot, 'BY_REGION', where, plot)
if (err) return

! set the who

plot%who(:)%name = ' '  ! default

do i = 1, size(who)

  if (who(i)(1:1) == '+') then
    plot%who(i)%sign = +1
  elseif (who(i)(1:1) == '-') then
    plot%who(i)%sign = -1
  endif

  call string_trim (who(i)(2:), plot%who(i)%name, j)
 
enddo

end subroutine 




