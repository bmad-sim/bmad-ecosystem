!+
! Subroutine tao_write_cmd (what)
!
! Routine to write output to a file or files or send the output to the printer.
! 
! Input:
!   what -- Character(*): What to output. See the code for more details.
!-

subroutine tao_write_cmd (what)

use quick_plot
use write_lat_file_mod
use blender_interface_mod
use tao_mod, dummy => tao_write_cmd
use tao_command_mod, dummy2 => tao_write_cmd
use tao_plot_mod, dummy3 => tao_write_cmd
use tao_top10_mod, dummy4 => tao_write_cmd
use madx_ptc_module, only: m_u, m_t, print_universe_pointed, print_complex_single_structure, print_new_flat, print_universe

implicit none

type (tao_curve_array_struct), allocatable, save :: curve(:)
type (tao_curve_struct), pointer :: c
type (beam_struct), pointer :: beam
type (bunch_struct), pointer :: bunch
type (branch_struct), pointer :: branch
type (tao_universe_struct), pointer :: u
type (ele_pointer_struct), allocatable, save :: eles(:)
type (ele_struct), pointer :: ele
type (coord_struct), pointer :: p

real(rp) scale

character(*) what
character(20) action, name, lat_type, which
character(40) switch
character(20) :: r_name = 'tao_write_cmd'
character(200) file_name0, file_name, what2
character(200) :: word(10)

character(20) :: names(24) = [ &
      'hard             ', 'gif              ', 'ps               ', 'variable         ', &
      'bmad_lattice     ', 'derivative_matrix', 'digested         ', 'curve            ', &
      'mad_lattice      ', 'beam             ', 'ps-l             ', 'hard-l           ', &
      'covariance_matrix', 'orbit            ', 'mad8_lattice     ', 'madx_lattice     ', &
      'pdf              ', 'pdf-l            ', 'opal_lattice     ', '3d_model         ', &
      'gif-l            ', 'ptc              ', 'sad_lattice      ', 'blender          ']

integer i, j, n, ie, ix, iu, nd, ii, i_uni, ib, ip, ios, loc
integer i_chan, ix_beam, ix_word, ix_w2

logical is_open, ascii, ok, err, good_opt_only, at_switch

!

call string_trim (what, what2, ix)
action = what2(1:ix)
call string_trim(what2(ix+1:), what2, ix_w2)

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
  is_open = .false.
  at_switch = .false.
  ix_word = 0

  do 
    ix_word = ix_word + 1
    if (ix_word == size(word)-1) exit

    call tao_next_switch (word(ix_word), ['-ascii', '-at   '], switch, err, ix)
    if (err) return

    select case (switch)
    case ('');       exit
    case ('-ascii'); ascii = .true.
    case ('-at')
      ix_word = ix_word + 1
      call tao_locate_elements (word(ix_word), s%com%default_universe, eles, err)
      if (err .or. size(eles) == 0) return
      at_switch = .true.
    end select
  enddo

  if (word(ix_word) /= '') then
    file_name0 = word(ix_word)
    if (word(ix_word+1) /= '' .or. file_name0(1:1) == '-') then
      call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
      return
    endif
  endif

  if (.not. at_switch) then
    call out_io (s_error$, r_name, 'YOU NEED TO SPECIFY "-at".')
    return
  endif 

  iu = lunget()

  uni_loop: do i = lbound(s%u, 1), ubound(s%u, 1)

    u => s%u(i)

    if (.not. tao_subin_uni_number (file_name0, i, file_name)) return
    call fullfilename (file_name, file_name)

    do ie = 1, size(eles)

      ele => eles(ie)%ele
      ! Write file

      beam => u%uni_branch(ele%ix_branch)%ele(ele%ix_ele)%beam
      if (.not. allocated(beam%bunch)) cycle

      if (.not. is_open) then
        if (ascii) then
          open (iu, file = file_name)
        else
          open (iu, file = file_name, form = 'unformatted')
          write (iu) '!BIN::2'
        endif
        is_open = .true.
      endif

      if (ascii) then
        write (iu, *) ele%ix_ele, '  ! ix_ele' 
        write (iu, *) size(beam%bunch), '  ! n_bunch'
        write (iu, *) size(beam%bunch(1)%particle), '  ! n_particle'
        do ib = 1, size(beam%bunch)
          bunch => beam%bunch(ib)
          write (iu, *) 'BEGIN_BUNCH'
          write (iu, *) '  ', trim(particle_name(bunch%particle(1)%species))
          write (iu, *) bunch%charge_tot, '  ! bunch_charge_tot'
          write (iu, *) bunch%z_center,   '  ! z_center'
          write (iu, *) bunch%t_center,   '  ! t_center'
          do ip = 1, size(bunch%particle)
            p => bunch%particle(ip)
            write (iu, '(6es19.10, es14.5, i6, 2(a, es19.10, a, es19.10, a), 3i3)') &
                  p%vec, p%charge, p%state, ('  (', real(p%spin(j)), ',', aimag(p%spin(j)), ')', j = 1, 2), &
                  p%ix_ele, p%location
          enddo
          write (iu, *) 'END_BUNCH'
        enddo
      else
        write (iu) ie, size(beam%bunch), size(beam%bunch(1)%particle)
        do ib = 1, size(beam%bunch)
          bunch => beam%bunch(ib)
          write (iu) bunch%particle(1)%species, bunch%charge_tot, bunch%z_center, bunch%t_center, size(bunch%particle)
          do ip = 1, size(bunch%particle)
            p => bunch%particle(ip)
            write (iu) p%vec, p%charge, p%state, p%spin, p%ix_ele, p%location
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

  enddo uni_loop


!---------------------------------------------------
! 3D model script for Blender

case ('3d_model', 'blender')

  file_name0 = 'blender_lat_#.py'
  if (word(1) /= '') file_name0 = word(1) 

  if (word(2) /= '') then
    call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
    return
  endif

  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. tao_subin_uni_number (file_name0, i, file_name)) return
    call write_blender_lat_layout (file_name, s%u(i)%model%lat)
    call out_io (s_info$, r_name, 'Written: ' // file_name)
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
    if (.not. tao_subin_uni_number (file_name0, i, file_name)) return
    call write_bmad_lattice_file (file_name, s%u(i)%model%lat, err)
    if (err) return
    call out_io (s_info$, r_name, 'Writen: ' // file_name)
  enddo

!---------------------------------------------------
case ('covariance_matrix')

  if (.not. allocated (s%com%covar)) then
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

  do i = 1, s%n_var_used
    if (.not. s%var(i)%useit_opt) cycle
    write (iu, '(i7, 3x, a)') s%var(i)%ix_dvar, tao_var1_name(s%var(i))
  enddo

  write (iu, *)
  write (iu, *) '!   i     j    Covar_Mat    Alpha_Mat'

  do i = 1, ubound(s%com%covar, 1)
    do j = 1, ubound(s%com%covar, 2)
      write (iu, '(2i6, 2es13.4)') i, j, s%com%covar(i,j), s%com%alpha(i,j)
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
    if (i_uni == 0) i_uni = s%com%default_universe
    beam => s%u(i_uni)%uni_branch(c%ix_branch)%ele(c%ix_ele_ref_track)%beam
    call file_suffixer (file_name, file_name, 'particle_dat', .true.)
    open (iu, file = file_name)
    write (iu, '(a, 6(12x, a))') '  Ix', '  x', 'px', '  y', 'py', '  z', 'pz'
    do i = 1, size(beam%bunch(1)%particle)
      write (iu, '(i6, 6es15.7)') i, (beam%bunch(1)%particle(i)%vec(j), j = 1, 6)
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

  do i = 1, s%n_var_used
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
    if (.not. tao_subin_uni_number (file_name0, i, file_name)) return
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
  call tao_draw_plots ()   ! PS out
  call qp_close_page
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
! Foreign lattice format

case ('mad_lattice', 'mad8_lattice', 'madx_lattice', 'opal_latice', 'sad_lattice')

  if (word(1) == '') then
    select case (action)
    case ('mad_lattice');   file_name0 = 'lat_#.mad8'; lat_type = 'MAD-8'
    case ('mad8_lattice');  file_name0 = 'lat_#.mad8'; lat_type = 'MAD-8'
    case ('madx_lattice');  file_name0 = 'lat_#.madX'; lat_type = 'MAD-X'
    case ('opal_latice');   file_name0 = 'lat_#.opal'; lat_type = 'OPAL-T'
    case ('sad_lattice');   file_name0 = 'lat_#.sad';  lat_type = 'SAD'
    end select
  else
    file_name0 = word(1)
  endif

  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. tao_subin_uni_number (file_name0, i, file_name)) return
    call write_lattice_in_foreign_format (lat_type, file_name, s%u(i)%model%lat, &
                                             s%u(i)%model%lat_branch(0)%orbit, err = err)
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
        ['-beam_index', '-design    ', '-base      '], n, .true., .true., name)
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

case ('ps', 'ps-l', 'gif', 'gif-l', 'pdf', 'pdf-l')

  if (qp_base_library == 'PGPLOT' .and. action(1:3) == 'pdf') then
    call out_io (s_error$, r_name, 'PGPLOT DOES NOT SUPPORT PDF!')
    return
  endif

  file_name = "tao.ps"
  if (action(1:3) == 'gif') file_name = 'tao.gif'
  if (action(1:3) == 'pdf') file_name = 'tao.pdf'
  call str_upcase (action, action)

  ix_word = 0
  scale = 0
  do
    ix_word = ix_word + 1
    if (ix_word == size(word)-1) exit

    call tao_next_switch (word(ix_word), ['-scale'], switch, err, ix)
    if (err) return

    select case (switch)
    case ('');  exit
    case ('-scale')
      ix_word = ix_word + 1
      read (word(ix_word), *, iostat = ios) scale
      if (ios /= 0 .or. word(ix_word) == '') then
        call out_io (s_error$, r_name, 'BAD SCALE NUMBER.')
        return
      endif
    end select
  enddo

  if (word(ix_word) /= '') then
    file_name = word(ix_word)
    if (word(ix_word+1) /= '' .or. file_name(1:1) == '-') then
      call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
      return
    endif
  endif

  if (action(1:3) == 'gif') then
    call qp_open_page (action, plot_file = file_name, x_len = s%plot_page%size(1), y_len = s%plot_page%size(2), scale = scale)
  else
    call qp_open_page (action, plot_file = file_name, scale = scale)
  endif
  call tao_draw_plots (.false.)   ! GIF plot
  call qp_close_page

  call tao_draw_plots ()   ! Update the plotting window

  call out_io (s_blank$, r_name, "Created " // trim(action) // " file: " // file_name)

!---------------------------------------------------
! ptc

case ('ptc')

which = '-old'
u => tao_pointer_to_universe(-1)
branch => u%model%lat%branch(0)

do 
  call tao_next_switch (what2, ['-old   ', '-new   ', '-branch', '-all   '], switch, err, ix_w2)
  if (err) return
  if (switch == '') exit

  select case (switch)
  case ('-old', '-new', '-all')
    which = switch
  case ('-branch')
    branch => pointer_to_branch (what2(1:ix_w2), u%model%lat)
    if (.not. associated(branch)) then
      call out_io (s_fatal$, r_name, 'Bad branch name or index: ' // what2(:ix_w2))
      return
    endif
  end select
enddo

file_name = what2(:ix_w2)
if (file_name == '') file_name = 'ptc.flatfile'

select case (which)
case ('-old', '-new')
  if (.not. associated(branch%ptc%m_t_layout)) then
    call out_io (s_fatal$, r_name, 'No associated PTC layout exists.', &
                                  'You must use the command "ptc init" before creating a flat file.')
    return
  endif

  if (which == '-old') then
    call print_complex_single_structure (branch%ptc%m_t_layout, file_name)
  else
    call print_new_flat (branch%ptc%m_t_layout, file_name)
  endif

  call out_io (s_info$, r_name, 'Writen: ' // file_name)

case ('-all')
  call print_universe (M_u, trim(file_name) // '.m_u')
  call print_universe_pointed (M_u, M_t, trim(file_name) // '.m_t')
  call out_io (s_info$, r_name, 'Writen: ' // trim(file_name) // '.m_u')
  call out_io (s_info$, r_name, 'Writen: ' // trim(file_name) // '.m_t')
end select

!---------------------------------------------------
! variables

case ('variable')

  good_opt_only = .false.
  ix_word = 0

  do 
    ix_word = ix_word + 1
    if (ix_word >= size(word)-1) exit
    call tao_next_switch (word(ix_word), ['-good_opt_only'], switch, err, ix)
    if (err) return
    select case (switch)
    case (''); exit
    case ('-good_opt_only'); good_opt_only = .true.
    end select
  enddo  

  if (word(ix_word) == ' ') then
    call tao_var_write (s%global%var_out_file, good_opt_only)
  else
    call tao_var_write (word(ix_word), good_opt_only)
  endif

!---------------------------------------------------
! error

case default

  call out_io (s_error$, r_name, 'UNKNOWN "WHAT": ' // what)

end select

end subroutine 
