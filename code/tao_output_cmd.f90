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
use write_lat_file_mod
use tao_command_mod

implicit none

type (tao_curve_array_struct), allocatable, save :: curve(:)
type (tao_curve_struct), pointer :: c
type (beam_struct), pointer :: beam
type (bunch_struct), pointer :: bunch
type (tao_universe_struct), pointer :: u

real(rp) scale

character(*) what
character(20) action, name
character(40) switch
character(20) :: r_name = 'tao_output_cmd'
character(100) file_name0, file_name, what2
character(80) :: word(10)

character(20) :: names(14) = (/ &
      'hard             ', 'gif              ', 'ps               ', 'variable         ', &
      'bmad_lattice     ', 'derivative_matrix', 'digested         ', 'curve            ', &
      'mad_lattice      ', 'beam             ', 'ps-l             ', 'hard-l           ', &
      'covariance_matrix', 'orbit            ' /)

integer i, j, n, ix, iu, nd, ii, i_uni, ib, ip, ios, loc
integer i_chan, ix_beam
integer, allocatable, save :: ix_ele_at(:)

logical is_open, ascii, ok, err, good_opt_only

!

call string_trim (what, what2, ix)
action = what2(1:ix)
what2 = what2(ix+1:)

call tao_cmd_split (what2, 10, word, .true., err)
if (err) return

call match_word (action, names, ix, .true.)

if (ix == 0) then
  call out_io (s_error$, r_name, 'UNRECOGNIZED "WHAT": ' // action)
  return
elseif (ix < 0) then
  call out_io (s_error$, r_name, 'AMBIGUOUS "WHAT": ' // action)
  return
endif
action = names(ix)

select case (action)

!---------------------------------------------------
! beam

case ('beam')

  ascii = .false.
  file_name0 = 'beam_#.dat'
  loc = -1
  is_open = .false.

  do
    call tao_next_switch (what2, (/ '-ascii', '-at   ' /), switch, err, ix)
    if (err) return
    if (switch == '') exit
    if (switch == '-ascii') ascii = .true.
    if (switch == '-at') then
      call tao_locate_elements (what2(1:ix), s%global%u_view, ix_ele_at)
      what2 = what2(ix+1:)
      loc = ix_ele_at(1)
      if (loc < 0) return
    endif
  enddo

  if (what2 /= '') then
    file_name0 = what2(1:ix)
    if (what2(ix+1:) /= '') then
      call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
      return
    endif
  endif

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
        write (iu, *) 'BEGIN_BUNCH'
        write (iu, *) j, '  ! ix_ele' 
        write (iu, *) size(beam%bunch), '  ! n_bunch'
        write (iu, *) size(beam%bunch(1)%particle), '  ! n_particle'
        do ib = 1, size(beam%bunch)
          bunch => beam%bunch(ib)
          write (iu, *) bunch%charge, '  ! bunch_charge'
          write (iu, *) bunch%z_center, '  ! z_center'
          write (iu, *) bunch%t_center, '  ! t_center'
          do ip = 1, size(bunch%particle)
            write (iu, '(6es19.10, es14.5, i6, 4es19.10)') &
                          bunch%particle(ip)%r%vec, bunch%particle(ip)%charge, &
                          bunch%particle(ip)%ix_lost, bunch%particle(ip)%r%spin 
          enddo
          write (iu, *) 'END_BUNCH'
        enddo
      else
        write (iu) j, size(beam%bunch), size(beam%bunch(1)%particle)
        do ib = 1, size(beam%bunch)
          bunch => beam%bunch(ib)
          write (iu) bunch%charge, bunch%z_center, bunch%t_center, size(bunch%particle)
          do ip = 1, size(bunch%particle)
            write (iu) bunch%particle(ip)%r%vec, bunch%particle(ip)%charge, &
                               bunch%particle(ip)%ix_lost, bunch%particle(ip)%r%spin
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
! bmad_lattice

case ('bmad_lattice')

  file_name0 = 'lat_#.bmad'
  if (word(1) /= '') file_name0 = word(1) 

  if (word(2) /= '') then
    call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
    return
  endif

  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. subin_uni_number (file_name0, i, file_name)) return
    call write_bmad_lattice_file (file_name, s%u(i)%model%lat, err)
    if (err) return
    call out_io (s_info$, r_name, 'Writen: ' // file_name)
  enddo

!---------------------------------------------------
case ('covariance_matrix')

  if (.not. allocated (tao_com%covar)) then
    call out_io (s_error$, r_name, 'COVARIANCE MATRIX NOT YET CALCULATED!')
    return
  endif

  file_name = 'lat_#.bmad'
  if (word(1) /= '') file_name = word(1) 

  if (word(2) /= '') then
    call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
    return
  endif

  call fullfilename (file_name, file_name)

  iu = lunget()
  open (iu, file = file_name)

  write (iu, '(i7, 2x, a)') count(s%var%useit_opt), '! n_var'

  write (iu, *)
  write (iu, *) '! Index   Variable'

  do i = 1, size(s%var)
    if (.not. s%var(i)%useit_opt) cycle
    write (iu, '(i7, 3x, a)') s%var(i)%ix_dvar, tao_var1_name(s%var(i))
  enddo

  write (iu, *)
  write (iu, *) '!   i     j    Covar_Mat    Alpha_Mat'

  do i = 1, ubound(tao_com%covar, 1)
    do j = 1, ubound(tao_com%covar, 2)
      write (iu, '(2i6, 2es13.4)') i, j, tao_com%covar(i,j), tao_com%alpha(i,j)
    enddo
  enddo

  call out_io (s_info$, r_name, 'Writen: ' // file_name)
  close(iu)

!---------------------------------------------------
! curve

case ('curve')

  call tao_find_plots (err, word(1), 'BOTH', curve = curve, always_allocate = .true.)
  if (err .or. size(curve) == 0) then
    call out_io (s_error$, r_name, 'CANNOT FIND CURVE')
    return
  endif

  if (size(curve) > 1) then
    call out_io (s_error$, r_name, 'MULTIPLE CURVES FIT NAME')
    return
  endif

  file_name = 'curve'
  if (word(2) /= ' ') file_name = word(2)
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

  file_name = word(1)
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

  file_name0 = word(1)
  if (file_name0 == ' ') file_name0 = 'digested_lat_universe_#.bmad'

  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. subin_uni_number (file_name0, i, file_name)) return
    call write_digested_bmad_file (file_name, s%u(i)%model%lat)
    call out_io (s_info$, r_name, 'Writen: ' // file_name)
  enddo

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
! mad_lattice

case ('mad_lattice')

  file_name0 = word(1)
  if (file_name0 == ' ') file_name0 = 'lat_#.mad'

  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. subin_uni_number (file_name0, i, file_name)) return
    call bmad_to_mad (file_name, s%u(i)%model%lat, err = err)
    if (err) return
    call out_io (s_info$, r_name, 'Writen: ' // file_name)
  enddo

!---------------------------------------------------
! orbit

case ('orbit')

  file_name0 = 'orbit.dat'
  i = 0
  do while (i <= size(word))
    i = i + 1
    if (word(i) == '') exit
    call match_word (word(i), &
                      (/ '-beam_index', '-design    ', '-base      ' /), n, .true., name)
    if (n < 0 .or. (n == 0 .and. word(i)(1:1) == '-')) then
      call out_io (s_error$, r_name, 'AMBIGUOUS SWITCH: ' // word(i))
      return
    endif
    select case (name)
    case ('-beam_index') 
      i=i+1; read (word(i), *, iostat = ios) ix_beam
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'CANNOT READ BEAM INDEX.')
        return
      endif
      action = name
    case ('-design')
      action = name
    case ('-base')
      action = name
    case default
      i=i+1; file_name0 = word(i)
      if (word(i+1) /= '') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
        return
      endif
    end select
  enddo

  if (i < size(word)) then
    call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
    return
  endif

  !

  u => tao_pointer_to_universe (-1) 

  iu = lunget()
  open (iu, file = file_name0)

  write (iu, '(a)') '&particle_orbit'

  do i = 0, u%model%lat%n_ele_track
    select case (action)
    case ('-beam_index') 
      
    case ('-design')
    case ('-base')
    end select
  enddo

  write (iu, '(a)') '/'
  close (iu)

!---------------------------------------------------
! ps

case ('ps', 'ps-l', 'gif', 'gif-l')

  file_name = "tao.ps"
  if (action(1:3) == 'gif') file_name = 'tao.gif'
  scale = 0
  call str_upcase (action, action)

  do
    call tao_next_switch (what2, (/ '-scale' /), switch, err, ix)
    if (err) return
    if (switch == '') exit
    if (switch == '-scale') then
      read (what2(1:ix), *, iostat = ios) scale
      if (ios /= 0 .or. what2 == '') then
        call out_io (s_error$, r_name, 'BAD SCALE NUMBER.')
        return
      endif
      what2 = what2(ix+1:)
    endif
  enddo

  if (what2 /= '') then
    file_name = what2(1:ix)
    if (what2(ix+1:) /= '') then
      call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
      return
    endif
  endif

  if (action(1:3) == 'GIF') then
    call qp_open_page ('PS' // trim(action(4:)), plot_file = 'tao_out.ps', scale = scale)
  else
    call qp_open_page (action, plot_file = file_name, scale = scale)
  endif

  call tao_draw_plots ()   ! Update the plotting window
  call qp_close_page
  call qp_select_page (s%plot_page%id_window)  ! Back to X-windows
  call tao_draw_plots ()   ! Update the plotting window

  if (action(1:3) == 'GIF') then
    call ps2gif ('tao_out.ps', file_name, .true.)
    call out_io (s_blank$, r_name, "Created GIF file: " // file_name)
  else
    call out_io (s_blank$, r_name, "Created PS file: " // file_name)
  endif

!---------------------------------------------------
! variables

case ('variable')

  good_opt_only = .false.
  do 
    call tao_next_switch (what2, (/ '-good_opt_only' /), switch, err, ix)
    if (err) return
    if (switch == '') exit
    if (switch == '-good_opt_only') good_opt_only = .true.
  enddo  

  if (what2 == ' ') then
    call tao_var_write (s%global%var_out_file, good_opt_only)
  else
    call tao_var_write (what2, good_opt_only)
  endif

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
