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
  type (res_struct) res

  character(1) ans

  namelist / res_params / param

! get the parameters

  param%y_min = 0; param%y_max = 0

  open (unit = 1, file = 'tune_plane_res_plot.in', status = 'old')
  read (1, nml = res_params)
  close (1)

! calculate resonance lines and plot

  call res_line_calc (param, res)

! Write results to a file.
  
  open (1, file = 'tune_plane_res_plot.dat')
  call res_line_write_list (1, param, res)
  close (1)
  print *, 'Resonance Line list in: tune_plane_res_plot.dat'

! plot results

  call res_line_plot ('PS-L', param, res)
  call qp_close_page

  call res_line_plot ('X', param, res)

  write (*, '(a)', advance = 'NO') ' Hit any key to end program: '
  read (*, '(a)') ans

end program
