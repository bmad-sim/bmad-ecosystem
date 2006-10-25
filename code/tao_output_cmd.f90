!+
! Subroutine tao_output_cmd (what, arg1, arg2)
!
! 
! Input:
!
!  Output:
!-

subroutine tao_output_cmd (what, arg1, arg2)

use tao_mod
use tao_top10_mod
use quick_plot
use tao_plot_mod
use io_mod

implicit none

type (tao_curve_array_struct), allocatable, save :: curve(:)
type (tao_curve_struct), pointer :: c

character(*) what
character(*) arg1, arg2
character(20) action
character(20) :: r_name = 'tao_output_cmd'
character(100) file_name

character(20) :: names(8) = (/ &
      'hard             ', 'gif              ', 'ps               ', 'variable         ', &
      'lattice          ', 'derivative_matrix', 'digested         ', 'curve            ' /)

integer i, j, ix, iu, nd, ii
logical err

!

call string_trim (what, action, ix)
call match_word (action, names, ix)
if (ix == 0) then
  call out_io (s_error$, r_name, 'UNRECOGNIZED "WHAT": ' // action)
  return
elseif (ix < 0) then
  call out_io (s_error$, r_name, 'AMBIGUOUS "WHAT": ' // action)
  return
endif
action = names(ix)

select case (action)

! hard

case ('hard')
  call qp_open_page ('PS')
  call tao_plot_out ()   ! Update the plotting window
  call qp_close_page
  call qp_select_page (s%plot_page%id_window)  ! Back to X-windows
  call tao_plot_out ()   ! Update the plotting window

  if (s%global%print_command == ' ') then
    call out_io (s_fatal$, r_name, &
        'P%PRINT_COMMAND NEEDS TO BE SET TO SEND THE PS FILE TO THE PRINTER!')
    return
  endif

  call system (trim(s%global%print_command) // ' quick_plot.ps')
  call out_io (s_blank$, r_name, 'Printing with command: ' // &
                                              s%global%print_command)

case ('gif')
  file_name = arg1
  if (file_name == "") file_name = "tao.gif"
  call qp_open_page ('GIF', x_len = s%plot_page%size(1), &
           y_len = s%plot_page%size(2), units = 'POINTS', plot_file = file_name)
  call tao_plot_out ()   ! Update the plotting window
  call qp_close_page
  call qp_select_page (s%plot_page%id_window)  ! Back to X-windows
  call tao_plot_out ()   ! Update the plotting window
  call out_io (s_info$, r_name, "Created GIF file: " // file_name)

! ps

case ('ps')
  file_name = arg1
  if (file_name == "") file_name = "tao.ps"
  call qp_open_page ('PS', plot_file = file_name)
  call tao_plot_out ()   ! Update the plotting window
  call qp_close_page
  call qp_select_page (s%plot_page%id_window)  ! Back to X-windows
  call tao_plot_out ()   ! Update the plotting window
  call out_io (s_blank$, r_name, "Created PS file: " // file_name)

! variables

case ('variable')
  if (arg1 == ' ') then
    call tao_var_write (s%global%var_out_file)
  else
    call tao_var_write (arg1)
  endif

case ('lattice')
  file_name = arg1
  if (file_name == ' ') file_name = 'lat_universe_#.bmad'
  do i = 1, size(s%u)
    ix = index(file_name, '#')
    if (ix /= 0) write (file_name, '(a, i0, a)') file_name(1:ix-1), i, trim(file_name(ix+1:))
    call write_bmad_lattice_file (file_name, s%u(i)%model%lat)
    call out_io (s_info$, r_name, 'Writen: ' // file_name)
  enddo

case ('derivative_matrix')

  nd = 0
  do i = 1, size(s%u)  
    if (.not. s%u(i)%is_on) cycle
    nd = nd + count(s%u(i)%data%useit_opt)
    if (.not. associated(s%u(i)%dmodel_dvar)) then
      call out_io (s_error$, r_name, 'DERIVATIVE MATRIX NOT YET CALCULATED!')
      return
    endif
  enddo

  file_name = arg1
  if (file_name == ' ') file_name = 'derivative_matrix.dat'
  iu = lunget()
  open (iu, file = file_name)

  write (iu, *) count(s%var%useit_opt), '  ! n_var'
  write (iu, *) nd, '  ! n_data'

  write (iu, *)
  write (iu, *) '! Index   Variable'

  do i = 1, size(s%var)
    if (.not. s%var(i)%useit_opt) cycle
    write (iu, '(i7, 3x, a)') s%var(i)%ix_dvar, tao_var1_name(s%var(i))
  enddo

  write (iu, *)
  write (iu, *) '! Index   Data'

  do i = 1, size(s%u)
    if (.not. s%u(i)%is_on) cycle
    do j = 1, size(s%u(i)%data)
      if (.not. s%u(i)%data(j)%useit_opt) cycle
      write (iu, '(i7, 3x, a)') s%u(i)%data(j)%ix_dModel, tao_datum_name(s%u(i)%data(j))
    enddo
  enddo

  write (iu, *)
  write (iu, *) ' ix_dat ix_var  dModel_dVar'
  nd = 0
  do i = 1, size(s%u)
    if (.not. s%u(i)%is_on) cycle
    do ii = 1, size(s%u(i)%dmodel_dvar, 1)
      do j = 1, size(s%u(i)%dmodel_dvar, 2)
        write (iu, '(2i7, es15.5)') nd + ii, j, s%u(i)%dmodel_dvar(ii, j)
      enddo
    enddo
    nd = nd + count(s%u(i)%data%useit_opt)
  enddo


  call out_io (s_info$, r_name, 'Writen: ' // file_name)
  close(iu)

case ('digested')
  file_name = arg1
  if (file_name == ' ') file_name = 'digested_lat_universe_#.bmad'
  do i = 1, size(s%u)
    ix = index(file_name, '#')
    if (ix /= 0) write (file_name, '(a, i0, a)') file_name(1:ix-1), i, trim(file_name(ix+1:))
    call write_digested_bmad_file (file_name, s%u(i)%model%lat)
    call out_io (s_info$, r_name, 'Writen: ' // file_name)
  enddo

! curve

case ('curve')

  call tao_find_plots (err, arg1, 'BOTH', curve = curve)
  if (.not. allocated(curve)) then
    call out_io (s_error$, r_name, 'NO CURVE SPECIFIED.')
    return
  endif

  iu = lunget()

  file_name = 'curve'
  if (arg2 /= ' ') file_name = arg2
  c => curve(1)%c

  if (allocated(c%beam%bunch)) then
    call file_suffixer (file_name, file_name, 'particle_dat', .true.)
    open (iu, file = file_name)
    write (iu, '(a, 6(12x, a))') '  Ix', '  x', 'p_x', '  y', 'p_y', '  z', 'p_z'
    do i = 1, size(c%beam%bunch(1)%particle)
      write (iu, '(i6, 6es15.7)') i, (c%beam%bunch(1)%particle(i)%r%vec(j), j = 1, 6)
    enddo
    call out_io (s_info$, r_name, 'Writen: ' // file_name)
    close(iu)
  endif

  call file_suffixer (file_name, file_name, 'symbol_dat', .true.)
  open (iu, file = file_name)
  write (iu, '(a, 6(12x, a))') '  Ix', '  x', '  y'
  do i = 1, size(c%x_symb)
    write (iu, '(i6, 2es15.7)') i, c%x_symb(i), c%y_symb(i)
  enddo
  call out_io (s_info$, r_name, 'Writen: ' // file_name)
  close(iu)

  call file_suffixer (file_name, file_name, 'line_dat', .true.)
  open (iu, file = file_name)
  write (iu, '(a, 6(12x, a))') '  Ix', '  x', '  y'
  do i = 1, size(c%x_line)
    write (iu, '(i6, 2es15.7)') i, c%x_line(i), c%y_line(i)
  enddo
  call out_io (s_info$, r_name, 'Writen: ' // file_name)
  close(iu)


! error

case default

  call out_io (s_error$, r_name, 'UNKNOWN "WHAT": ' // what)

end select

end subroutine 
