!+
! Program tune_plane_res_plot 
!
! Program to find the resonances in a bounded region of the tune plane.
! See the tune_plane_res_plot.in file for more details.
!
! resonances are of the form:
!                 p * Q_x + q * Q_y + r*Q_s = n
!
! with the restrictions that:
!                 |p| + |q|  + |r| =< pqr_max
!                 |p| + |q|        =< pq_max
!                 |p|              =< p_max
!                 |q|              =< q_max
!                 |r|              =< r_max
!
! and we search the tune plane in the rectangle:
!                 x_min <= Q_x <= x_max
!                 y_min <= Q_y <= y_max
!
! UNITS = 'kHz' => kHz units 
!       = '1'  => tune units (1 = 390 kHz)
!
! ANGLE is the orientation so that the label is parallel to the resonance line
!-

program tune_plane_res_plot

  use tune_plane_res_mod

  implicit none

  type (res_params_struct) param
  type (res_struct), target :: res
  type (res_line_struct), pointer :: line

  integer i

  character(1) ans

  namelist / res_params / param

! get the parameters

  param%y_min = 0; param%y_max = 0

  open (unit = 1, file = 'tune_plane_res_plot.in', status = 'old')
  read (1, nml = res_params)
  close (1)

! calculate resonance lines and plot

  call res_line_calc (param, res)

  
  open (1, file = 'tune_plane_res_plot.dat')

  write (1, '(a)')        '&res_params'
  write (1, '(a, f12.6)') '  param%x_min = ', param%x_min
  write (1, '(a, f12.6)') '  param%x_max = ', param%x_max
  write (1, '(a, f12.6)') '  param%y_min = ', param%y_min
  write (1, '(a, f12.6)') '  param%y_max = ', param%y_max
  write (1, '(a, i0)')    '  param%x_div = ', param%x_div
  write (1, '(a, i0)')    '  param%y_div = ', param%y_div
  write (1, '(a, f12.6)') '  param%Q_s   = ', param%Q_s
  write (1, '(a, a)')     '  param%units = ', param%units
  write (1, '(a, i0)')    '  param%pqr_max = ', param%pqr_max
  write (1, '(a, i0)')    '  param%pq_max  = ', param%pq_max
  write (1, '(a, i0)')    '  param%p_max   = ', param%p_max  
  write (1, '(a, i0)')    '  param%q_max   = ', param%q_max  
  write (1, '(a, i0)')    '  param%r_max   = ', param%r_max  
  write (1, '(a, i0)')    '  param%p_restrict = ', param%p_restrict
  write (1, '(a, i0)')    '  param%q_restrict = ', param%q_restrict
  write (1, '(a, i0)')    '  param%sum_diff   = ', param%sum_diff
  write (1, '(a, l1)')    '  param%plot_all_1st_and_2nd = ', param%plot_all_1st_and_2nd
  write (1, '(a, l1)')    '  param%show_labels          = ', param%show_labels
  write (1, '(a)')        '/'
  write (1, *)

  write (1, *) '    P     Q     R     N'
  do i = 1, res%num_line
    line => res%line(i)
    write (1, '(4i6)') line%p, line%q, line%r, line%n
  enddo
  print *, 'Resonance Line list in: tune_plane_res_plot.dat'
  close (1)

! plot results

  call res_line_plot ('PS-L', param, res)
  call qp_close_page

  call res_line_plot ('X', param, res)

  write (*, '(a)', advance = 'NO') ' Hit any key to end program: '
  read (*, '(a)') ans

end program
