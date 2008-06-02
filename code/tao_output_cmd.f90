!+
! Subroutine tao_output_cmd (what)
!
! 
! Input:
!
!  Output:
!-

subroutine tao_output_cmd (what)

use tao_mod
use tao_top10_mod
use quick_plot
use tao_plot_mod
use io_mod
use tao_command_mod

implicit none

type (tao_curve_array_struct), allocatable, save :: curve(:)
type (tao_curve_struct), pointer :: c
type (beam_struct), pointer :: beam
type (bunch_struct), pointer :: bunch
type (tao_universe_struct), pointer :: u
type (tao_arg_struct) arg(10)

character(*) what
character(20) action
character(20) :: r_name = 'tao_output_cmd'
character(100) file_name0, file_name

character(20) :: names(12) = (/ &
      'hard             ', 'gif              ', 'ps               ', 'variable         ', &
      'bmad_lattice     ', 'derivative_matrix', 'digested         ', 'curve            ', &
      'mad_lattice      ', 'beam             ', 'ps-l             ', 'hard-l           ' /)
integer :: n_arg_max(12) = (/ &
      1, 2, 2, 2, &
      2, 2, 2, 3, &
      2, 4, 2, 2 /)      

character(20) :: arg_names(2) = (/ '-ascii', '-at   ' /)
integer :: n_arg_values(2) = (/ 0, 1 /)

integer i, j, ix, iu, nd, ii, i_uni, ib, ip, ios, loc
integer n_arg, i_chan
integer, allocatable, save :: ix_ele_at(:)

logical is_open, ascii, ok, err

!

call tao_arg_split (what, arg_names, n_arg_values, arg, n_arg, err)
if (err) return

action = arg(1)%name
call match_word (action, names, ix, .true.)

if (ix == 0) then
  call out_io (s_error$, r_name, 'UNRECOGNIZED "WHAT": ' // action)
  return
elseif (ix < 0) then
  call out_io (s_error$, r_name, 'AMBIGUOUS "WHAT": ' // action)
  return
endif
action = names(ix)

! Make sure number of arguments is ok

if (n_arg_max(ix) < n_arg) then
  call out_io (s_error$, r_name, 'THERE ARE TOO MANY ARGUMENTS HERE.')
  return
endif

select case (action)

!---------------------------------------------------
! hard

case ('hard', 'hard-l')

  if (action == 'hard') then
    call qp_open_page ('PS', scale = 0.0_rp)
  else
    call qp_open_page ('PS-L', scale = 0.0_rp)
  endif
  call tao_draw_plots ()   ! Update the plotting window
  call qp_close_page
  call qp_select_page (s%plot_page%id_window)  ! Back to X-windows
  call tao_draw_plots ()   ! Update the plotting window

  if (s%global%print_command == ' ') then
    call out_io (s_fatal$, r_name, &
        'P%PRINT_COMMAND NEEDS TO BE SET TO SEND THE PS FILE TO THE PRINTER!')
    return
  endif

  call system (trim(s%global%print_command) // ' quick_plot.ps')
  call out_io (s_blank$, r_name, 'Printing with command: ' // &
                                              s%global%print_command)
!---------------------------------------------------
! gif

case ('gif')

  file_name = "tao.gif"
  if (n_arg == 2) file_name = arg(2)%name

  call qp_open_page ('GIF', i_chan, s%plot_page%size(1), &
                       s%plot_page%size(2), 'POINTS', file_name)
  call qp_select_page (i_chan)
  call tao_draw_plots ()   ! Update the plotting window
  call qp_close_page
  call qp_select_page (s%plot_page%id_window)  ! Back to X-windows
  call tao_draw_plots ()   ! Update the plotting window
  call out_io (s_info$, r_name, "Created GIF file: " // file_name)

!---------------------------------------------------
! ps

case ('ps', 'ps-l')

  file_name = "tao.ps"
  if (n_arg == 2) file_name = arg(2)%name

  if (action == 'ps') then
    call qp_open_page ('PS', plot_file = file_name, scale = 0.0_rp)
  else
    call qp_open_page ('PS-L', plot_file = file_name, scale = 0.0_rp)
  endif
  call tao_draw_plots ()   ! Update the plotting window
  call qp_close_page
  call qp_select_page (s%plot_page%id_window)  ! Back to X-windows
  call tao_draw_plots ()   ! Update the plotting window
  call out_io (s_blank$, r_name, "Created PS file: " // file_name)

!---------------------------------------------------
! variables

case ('variable')

  if (arg(2)%name == ' ') then
    call tao_var_write (s%global%var_out_file)
  else
    call tao_var_write (arg(2)%name)
  endif

!---------------------------------------------------
! bmad_lattice

case ('bmad_lattice')

  file_name0 = arg(2)%name
  if (file_name0 == ' ') file_name0 = 'lat_#.bmad'

  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. subin_uni_number (file_name0, i, file_name)) return
    call write_bmad_lattice_file (file_name, s%u(i)%model%lat, err)
    if (err) return
    call out_io (s_info$, r_name, 'Writen: ' // file_name)
  enddo

!---------------------------------------------------
! mad_lattice

case ('mad_lattice')

  file_name0 = arg(2)%name
  if (file_name0 == ' ') file_name0 = 'lat_#.mad'

  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. subin_uni_number (file_name0, i, file_name)) return
    call bmad_to_mad (file_name, s%u(i)%model%lat, err = err)
    if (err) return
    call out_io (s_info$, r_name, 'Writen: ' // file_name)
  enddo

!---------------------------------------------------
! derivative_matrix

case ('derivative_matrix')

  nd = 0
  do i = lbound(s%u, 1), ubound(s%u, 1)  
    if (.not. s%u(i)%is_on) cycle
    nd = nd + count(s%u(i)%data%useit_opt)
    if (.not. allocated(s%u(i)%dmodel_dvar)) then
      call out_io (s_error$, r_name, 'DERIVATIVE MATRIX NOT YET CALCULATED!')
      return
    endif
  enddo

  file_name = arg(2)%name
  if (file_name == ' ') file_name = 'derivative_matrix.dat'
  call fullfilename (file_name, file_name)

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

  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. s%u(i)%is_on) cycle
    do j = 1, size(s%u(i)%data)
      if (.not. s%u(i)%data(j)%useit_opt) cycle
      write (iu, '(i7, 3x, a)') s%u(i)%data(j)%ix_dModel, tao_datum_name(s%u(i)%data(j))
    enddo
  enddo

  write (iu, *)
  write (iu, *) ' ix_dat ix_var  dModel_dVar'
  nd = 0
  do i = lbound(s%u, 1), ubound(s%u, 1)
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

!---------------------------------------------------
! digested

case ('digested')

  file_name0 = arg(2)%name
  if (file_name0 == ' ') file_name0 = 'digested_lat_universe_#.bmad'

  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. subin_uni_number (file_name0, i, file_name)) return
    call write_digested_bmad_file (file_name, s%u(i)%model%lat)
    call out_io (s_info$, r_name, 'Writen: ' // file_name)
  enddo

!---------------------------------------------------
! curve

case ('curve')

  call tao_find_plots (err, arg(2)%name, 'BOTH', curve = curve)
  if (.not. allocated(curve)) then
    call out_io (s_error$, r_name, 'NO CURVE SPECIFIED.')
    return
  endif

  file_name = 'curve'
  if (arg(3)%name /= ' ') file_name = arg(3)%name
  call fullfilename (file_name, file_name)

  c => curve(1)%c
  iu = lunget()

  if (c%g%type == "phase_space") then
    i_uni = c%ix_universe
    if (i_uni == 0) i_uni = s%global%u_view
    beam => s%u(i_uni)%ele(c%ix_ele_ref_track)%beam
    call file_suffixer (file_name, file_name, 'particle_dat', .true.)
    open (iu, file = file_name)
    write (iu, '(a, 6(12x, a))') '  Ix', '  x', 'p_x', '  y', 'p_y', '  z', 'p_z'
    do i = 1, size(beam%bunch(1)%particle)
      write (iu, '(i6, 6es15.7)') i, (beam%bunch(1)%particle(i)%r%vec(j), j = 1, 6)
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

!---------------------------------------------------
! beam

case ('beam')

  ascii = .false.
  file_name0 = 'beam_#.dat'
  loc = -1
  is_open = .false.

  do i = 2, n_arg
    select case (arg(i)%name)
    case ('-ascii') 
      ascii = .true.
    case ('-at')
      call tao_locate_elements (arg(i)%value(1), s%global%u_view, ix_ele_at)
      loc = ix_ele_at(1)
      if (loc < 0) return
    case default
      file_name0 = arg(i)%name
    end select
  enddo

  iu = lunget()

  do i = lbound(s%u, 1), ubound(s%u, 1)

    u => s%u(i)
    if (.not. subin_uni_number (file_name0, i, file_name)) return
    call fullfilename (file_name, file_name)

  ! Write file

    do j = lbound(u%ele, 1), ubound(u%ele, 1)
      if (loc > -1 .and. loc /= j) cycle
      beam => u%ele(j)%beam
      if (.not. allocated(beam%bunch)) cycle

      if (.not. is_open) then
        if (ascii) then
          open (iu, file = file_name)
        else
          open (iu, file = file_name, form = 'unformatted')
          write (iu) '!BINARY'
        endif
        is_open = .true.
      endif

      if (ascii) then
        write (iu, *) j, u%beam_init%n_bunch, &
                              u%beam_init%n_particle
        do ib = 1, size(beam%bunch)
          bunch => beam%bunch(ib)
          write (iu, *) bunch%charge
          write (iu, *) bunch%z_center
          write (iu, *) bunch%t_center
          do ip = 1, size(bunch%particle)
            write (iu, '(6es19.10, i6, 4es19.10)') &
                          bunch%particle(ip)%r%vec, bunch%particle(ip)%ix_lost, &
                          bunch%particle(ip)%r%spin 
          enddo
        enddo
      else
        write (iu) j, u%beam_init%n_bunch, &
                              u%beam_init%n_particle, u%beam_init%bunch_charge
        do ib = 1, size(beam%bunch)
          bunch => beam%bunch(ib)
          write (iu) bunch%charge
          write (iu) bunch%z_center
          write (iu) bunch%t_center
          do ip = 1, size(bunch%particle)
            write (iu) bunch%particle(ip)%r%vec, bunch%particle(ip)%ix_lost, &
                       bunch%particle(ip)%r%spin
          enddo
        enddo
      endif

    enddo 

    if (is_open) then
      close (iu)
      call out_io (s_info$, r_name, 'Writen: ' // file_name)
    else
      call out_io (s_error$, r_name, 'NO ALLOCATED BEAM FOUND!')
    endif

  enddo

!---------------------------------------------------
! error

case default

  call out_io (s_error$, r_name, 'UNKNOWN "WHAT": ' // what)

end select

!---------------------------------------------------
contains

function subin_uni_number (name_in, ix_uni, name_out) result (ok)

  character(*) name_in, name_out
  integer ix, ix_uni
  logical ok

!

  ok = .true.
  name_out = name_in

  ix = index(name_out, '#')
  if (size(s%u) > 1 .and. ix == 0) then
    call out_io (s_info$, r_name, 'FILE NAME DOES NOT HAVE A "#" CHARACTER!', &
      ' YOU NEED THIS TO GENERATE A UNIQUE FILE NAME FOR EACH UNIVERSE!')
    ok = .false.
  endif

  if (ix /= 0) write (name_out, '(a, i0, a)') &
                        name_out(1:ix-1), ix_uni, trim(name_out(ix+1:))

end function

end subroutine 
