!+
! Subroutine tao_write_cmd (what)
!
! Routine to write output to a file or files or send the output to the printer.
! 
! Input:
!   what -- Character(*): What to output. See the code for more details.
!-

subroutine tao_write_cmd (what)

use tao_interface, dummy => tao_write_cmd
use tao_command_mod, only: tao_cmd_split, tao_next_switch, tao_next_word
use tao_plot_mod, only: tao_draw_plots
use tao_top10_mod, only: tao_var_write

use quick_plot, only: qp_open_page, qp_base_library, qp_close_page
use blender_interface_mod, only: write_blender_lat_layout
use madx_ptc_module, only: m_u, m_t, print_universe_pointed, &
                           print_new_flat, print_universe
use beam_file_io, only: write_beam_file
use ptc_layout_mod, only: ptc_emit_calc, lat_to_ptc_layout

implicit none

type (tao_curve_array_struct), allocatable :: curve(:)
type (tao_curve_struct), pointer :: c
type (tao_plot_struct), pointer :: tp
type (tao_plot_region_struct), pointer :: r
type (tao_universe_struct), pointer :: u
type (beam_struct), pointer :: beam
type (bunch_struct), pointer :: bunch
type (branch_struct), pointer :: branch
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), pointer :: ele
type (lat_param_struct), pointer :: param
type (coord_struct), pointer :: p
type (coord_struct)  :: orbit
type (tao_d2_data_struct), pointer :: d2
type (tao_d1_data_struct), pointer :: d1
type (tao_data_struct), pointer :: dat
type (tao_data_struct), target :: datum
type (tao_v1_var_struct), pointer :: v1
type (tao_spin_map_struct), pointer :: sm
type (bunch_track_struct), pointer :: bunch_params_comb(:)
type (bunch_track_struct), pointer :: comb1
type (bunch_params_struct), pointer :: bpt
type (em_field_struct) field

real(rp) scale, mat6(6,6)
real(rp), allocatable :: values(:)

character(*) what
character(1) delim
character(20) action, name, lat_type, which, last_col, b_name
character(40), allocatable :: z(:)
character(100) str, ele_name, c1, c2
character(200) line, switch, header1, header2, aname
character(200) file_name0, file_name, what2
character(200) :: word(40)
character(1000) :: aline
character(*), parameter :: r_name = 'tao_write_cmd'

real(rp) dr(3), r_max(3), r_min(3), rr(3)

integer i, j, k, m, n, ie, ic, id, ix, iu, nd, ii, i_uni, ib, ip, ios, loc, iy, iz
integer n_max(3), n_min(3), i_chan, ix_beam, ix_word, ix_w2, file_format
integer n_type, n_ref, n_start, n_ele, n_merit, n_meas, n_weight, n_good, n_bunch, n_eval, n_s
integer i_min, i_max, n_len, len_d_type, ix_branch, ix_bunch, n_loc

logical is_open, ok, err, good_opt_only, at_switch, new_file, append, write_floor
logical write_data_source, write_data_type, write_merit_type, write_weight, write_attribute, write_step
logical write_high_lim, write_low_lim, tao_format, eq_d_type, delim_found, found_plot_command

!

call string_trim (what, what2, ix)
action = what2(1:ix)
call string_trim(what2(ix+1:), what2, ix_w2)

call tao_cmd_split (what2, size(word), word, .true., err, ',')
if (err) return

call match_word (action, [character(20):: &
              '3d_model', 'beam', 'bmad', 'blender', 'bunch_comb', 'covariance_matrix', 'curve', &
              'derivative_matrix', 'digested', 'elegant', 'field', &
              'gif', 'gif-l', 'hard', 'hard-l', 'mad', 'mad8', 'madx', 'matrix', &
              'namelist', 'opal', 'pdf', 'pdf-l', 'plot_commands', 'ps', 'ps-l', 'ptc', &
              'sad', 'spin_mat8', 'tao', 'variable', 'xsif'], &
              ix, .true., matched_name = action)

if (ix == 0) then
  call out_io (s_error$, r_name, 'UNRECOGNIZED "WHAT": ' // action)
  return
elseif (ix < 0) then
  call out_io (s_error$, r_name, 'AMBIGUOUS "WHAT": ' // action)
  return
endif

iu = lunget()

select case (action)

!---------------------------------------------------
! beam

case ('beam')

  file_format = hdf5$
  is_open = .false.
  at_switch = .false.
  ix_word = 0
  file_name0 = ''
  write_floor = .false.

  do
    ix_word = ix_word + 1
    if (ix_word == size(word)-1) exit

    call tao_next_switch (word(ix_word), [character(16):: '-ascii', '-at', '-hdf5', '-floor_position'], .true., switch, err)
    if (err) return

    select case (switch)
    case ('');                 exit
    case ('-ascii');           file_format = ascii$
    case ('-floor_position');  write_floor = .true.
    case ('-hdf5');            file_format = hdf5$
    case ('-at')
      ix_word = ix_word + 1
      call tao_locate_elements (word(ix_word), s%global%default_universe, eles, err)
      if (err .or. size(eles) == 0) return
      at_switch = .true.
    case default
      if (file_name0 /= '') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
        return
      endif
      file_name0 = switch
    end select
  enddo

  if (write_floor) then
    file_name0 = 'beam_floor_#.dat'
  elseif (file_format == hdf5$) then
    if (file_name0 == '') then
      file_name0 = 'beam_#.hdf5'
    else
      n = len_trim(file_name0)
      if (file_name0(n-2:n) /= '.h5' .and. file_name0(n-4:n) /= '.hdf5') then
        file_name0 = trim(file_name0) // '.hdf5'
      endif
    endif

  elseif (file_name0 == '') then
    file_name0 = 'beam_#.dat'
  endif

  if (.not. at_switch) then
    call out_io (s_error$, r_name, 'YOU NEED TO SPECIFY "-at".')
    return
  endif 

  uni_loop: do i = lbound(s%u, 1), ubound(s%u, 1)
    u => s%u(i)

    if (.not. tao_subin_uni_number (file_name0, i, file_name)) return
    call fullfilename (file_name, file_name)
    new_file = .true.

    do ie = 1, size(eles)
      ele => eles(ie)%ele
      ! Write file

      beam => u%model_branch(ele%ix_branch)%ele(ele%ix_ele)%beam
      if (.not. allocated(beam%bunch)) cycle

      if (write_floor) then
        call write_beam_floor_positions(file_name, beam, ele, new_file)
      else
        call write_beam_file (file_name, beam, new_file, file_format, u%model%lat)
      endif
      new_file = .false.
    enddo 

    if (new_file) then
      call out_io (s_error$, r_name, 'BEAM NOT SAVED AT THIS ELEMENT.', &
                    'CHECK THE SETTING OF THE SAVED_AT COMPONENT OF THE TAO_BEAM_INIT NAMELIST.', &
                    'ANOTHER POSSIBILITY IS THAT GLOBAL%TRACK_TYPE = "single" SO NO BEAM TRACKING HAS BEEN DONE.')
    else
      call out_io (s_info$, r_name, 'Written: ' // file_name)
    endif

  enddo uni_loop


!---------------------------------------------------
! bmad

case ('bmad')

  file_format = binary$
  file_name0 = 'lat_#.bmad'
  ix_word = 0

  do 
    ix_word = ix_word + 1
    if (ix_word == size(word)-1) exit

    call tao_next_switch (word(ix_word), [character(16):: '-one_file', '-format'], .true., switch, err)
    if (err) return

    select case (switch)
    case ('');       exit
    case ('-one_file'); file_format = one_file$
    case ('-format')
      ix_word = ix_word + 1
      call tao_next_switch(word(ix_word), [character(16):: 'one_file', 'binary', 'ascii'], .true., switch, err)
      if (err) return
      select case (switch)
      case ('one_file');   file_format = one_file$
      case ('binary');     file_format = binary$
      case ('ascii');      file_format = ascii$
      case default
        call out_io (s_error$, r_name, 'UNKNOWN -format SWITCH: ' // word(ix_word))
        return
      end select

    case default
      if (file_name0 /= 'lat_#.bmad') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
        return
      endif
      file_name0 = switch
    end select
  enddo

  do i = lbound(s%u, 1), ubound(s%u, 1)
    u => s%u(i)
    if (.not. tao_subin_uni_number (file_name0, i, file_name)) return
    call write_bmad_lattice_file (file_name, u%model%lat, err, file_format, u%model%tao_branch(0)%orbit(0))
    if (err) return
    call out_io (s_info$, r_name, 'Written: ' // file_name)
  enddo

!---------------------------------------------------
! 3D model script for Blender
! Note: Old cubit interface code was in tao_write_3d_floor_plan.f90 which was deleted 9/2015.

case ('blender', '3d_model')

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
case ('bunch_comb')

  file_name = ''
  switch = '-sigma'
  i_uni = -1
  b_name = ''
  ix_bunch = 1

  do 
    ix_word = ix_word + 1
    call tao_next_switch (word(ix_word), [character(16):: '-sigma', '-min_max', &
             '-universe', '-centroid', '-ix_bunch', '-branch'], .true., switch, err)

    if (err) return

    select case (switch)
    case ('')
      exit
    case ('-sigma', '-min_max', '-centroid')
      which = switch
    case ('-universe')
      ix_word = ix_word + 1
      if (.not. is_integer(word(ix_word), i_uni)) then
        call out_io (s_error$, r_name, 'BAD UNIVERSE INDEX: ' // word(ix_word))
        return
      endif
    case ('ix_bunch')
      ix_word = ix_word + 1
      if (.not. is_integer(word(ix_word), ix_bunch)) then
        call out_io (s_error$, r_name, 'BAD BUNCH INDEX: ' // word(ix_word))
        return
      endif

    case ('-branch')
      ix_word = ix_word + 1
      b_name = word(ix_word)
    case default
      if (file_name /= '') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
        return
      endif
      file_name = switch
    end select
  enddo

  u => tao_pointer_to_universe(i_uni)
  if (.not. associated(u)) then
    call out_io (s_error$, r_name, 'BAD UNIVERSE INDEX: ' // word(i_uni))
    return
  endif

  branch => pointer_to_branch(b_name, u%model%lat, blank_branch = s%global%default_branch)
  if (.not. associated(branch)) then
    call out_io (s_error$, r_name, 'BAD LATTICE BRANCH NAME OR INDEX: ' // b_name)
    return
  endif

  if (.not. allocated(u%model%tao_branch(branch%ix_branch)%bunch_params_comb)) then
    call out_io (s_error$, r_name, 'COMB ARRAY NOT ALLOCATED. PROBABLY CAUSED BY NO BUNCH TRACKING.')
    return
  endif
  bunch_params_comb => u%model%tao_branch(branch%ix_branch)%bunch_params_comb

  if (ix_bunch > size(bunch_params_comb)) then
    call out_io (s_error$, r_name, 'IX_BUNCH INDEX OUT OF RANGE.')
    return
  endif
  comb1 => bunch_params_comb(ix_bunch)


  if (file_name == '') file_name = 'bunch_comb.' // which(2:)
  open (iu, file = file_name, recl = 500)

  select case (which)
  case ('-sigma')
    write (iu, '(a3, a7, 2a12, 2x, 21a14)') '# 1', '2', '3', '4', '5', '6', '7', '8', '9', '10', &
                    '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25'

    write (iu, '(a4, a7, 2a12, 2x, 21a14)') '# Ix', 'N-Live', 'S-Pos', 'Time', &
      '<x.x>', '<x.px>', '<x.y>', '<x.py>', '<x.z>', '<x.pz>', '<px.px>', '<px.y>', '<px.py>', '<px.z>', '<px.pz>', &
      '<y.y>', '<y.py>', '<y.z>', '<y.pz>', '<py.py>', '<py.z>', '<py.pz>', '<z.z>', '<z.pz>', '<pz.pz>'

    do ic = 0, comb1%n_pt
      bpt => comb1%pt(ic)
      write (iu, '(i4, i7, f12.6, es14.6, 2x, 21es14.6)') ic, bpt%n_particle_live, bpt%centroid%s, bpt%centroid%t, &
                            ((bpt%sigma(i,j), j = i,6), i = 1,6)
    enddo

  case ('-centroid')
    write (iu, '(a3, a7, 2a12, 4a14, 2x, 7a14, 2x, 6a14, 2x, 3a14, 2x, 3a14)') '# 1', '2', '3', '4', &
                    '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', &
                    '20', '21', '22', '23', '24', '25', '26', '27'

    write (iu, '(a4, a7, 2a12, a16, a12, 2a14, 2x, 7a14, 2x, 6a14, 2x, 3a14, 3x, 3a14)') '# Ix', &
                   'N-Live', 'S-Pos', 'Time', &
                   'Polarization', '<Sx>', '<Sy>', '<Sz>', &
                   '<x>', '<px>', '<y>', '<py>', '<z>', '<pz>', '<p0c>', &
                   'Sig_x', 'Sig_px', 'Sig_y', 'Sig_py', 'Sig_z', 'Sig_pz', &
                   'Emit_a', 'Emit_b', 'Emit_c', 'Norm_Emit_a', 'Norm_Emit_b', 'Norm_Emit_c'

    do ic = 0, comb1%n_pt
      bpt => comb1%pt(ic)
      write (iu, '(i4, i7, f12.6, es14.6, 4f14.9, 2x, 7es14.6, 2x, 6es14.6, 2(2x, 3es14.6))') ic, &
              bpt%n_particle_live, bpt%centroid%s, bpt%centroid%t, &
              norm2(bpt%centroid%spin), bpt%centroid%spin, bpt%centroid%vec, bpt%centroid%p0c, &
              bpt%a%sigma, bpt%a%sigma_p, bpt%b%sigma, bpt%b%sigma_p, bpt%c%sigma, bpt%c%sigma_p, &
              bpt%a%emit, bpt%b%emit, bpt%c%emit, bpt%a%norm_emit, bpt%b%norm_emit, bpt%c%norm_emit
    enddo

  case ('-min_max')
    write (iu, '(a3, a7, 2a12, 6(2x, 2a14))') '# 1', '2', '3', '4', '5', '6', '7', '8', '9', '10', &
                    '11', '12', '13', '14', '15', '16'

    write (iu, '(a4, a7, 2a12, 6(2x, 2a14))') '# Ix', 'N-Live', 'S-Pos', 'Time', &
                  'x_min', 'x_max', 'px_min', 'px_max', 'y_min', 'y_max', &
                  'py_min', 'py_max', 'z_min', 'z_max', 'pz_min', 'pz_max'

    do ic = 0, comb1%n_pt
      bpt => comb1%pt(ic)
      write (iu, '(i4, i7, f12.6, es14.6, 6(2x, 2es14.6))') ic, bpt%n_particle_live, bpt%centroid%s, bpt%centroid%t, &
                            (bpt%rel_min(i), bpt%rel_max(i), i = 1,6)
    enddo
  end select

  call out_io (s_info$, r_name, 'Written: ' // file_name)
  close (iu)

!---------------------------------------------------
case ('covariance_matrix')

  if (.not. allocated (s%com%covar)) then
    call out_io (s_error$, r_name, 'COVARIANCE MATRIX NOT YET CALCULATED!')
    return
  endif

  file_name = 'covar.matrix'
  if (word(1) /= '') file_name = word(1) 

  if (word(2) /= '') then
    call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
    return
  endif

  call fullfilename (file_name, file_name)
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

  call out_io (s_info$, r_name, 'Written: ' // file_name)
  close(iu)

!---------------------------------------------------
! curve

case ('curve')

  call out_io (s_info$, r_name, &
      '"show curve" command superseded by the more versatile "show -write <file> curve ..." command.')

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


  call out_io (s_info$, r_name, 'Written: ' // file_name)
  close(iu)

!---------------------------------------------------
! digested

case ('digested')

  file_name0 = word(1)
  if (file_name0 == ' ') file_name0 = 'lat_#.digested'

  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. tao_subin_uni_number (file_name0, i, file_name)) return
    call write_digested_bmad_file (file_name, s%u(i)%model%lat)
    call out_io (s_info$, r_name, 'Written: ' // file_name)
  enddo

!---------------------------------------------------
! field

case ('field')

  dr = real_garbage$
  r_min = real_garbage$
  r_max = real_garbage$
  n_min = int_garbage$
  n_max = int_garbage$
  file_name = ''
  ix_word = 0

  do
    ix_word = ix_word + 1
    if (ix_word == size(word)-1) exit    
    call tao_next_switch (word(ix_word), [character(16):: '-dr', '-nmax', '-nmin', &
                                                  '-rmax', '-rmin', '-ele'], .true., switch, err)
    if (err) return

    select case (switch)
    case ('');       exit
    case ('-dr');    if (.not. read_real3(word, ix_word, dr, '-DR')) return
    case ('-rmin');  if (.not. read_real3(word, ix_word, r_min, '-rmin')) return
    case ('-rmax');  if (.not. read_real3(word, ix_word, r_max, '-rmax')) return
    case ('-nmin');  if (.not. read_int3(word, ix_word, n_min, '-nmin')) return
    case ('-nmax');  if (.not. read_int3(word, ix_word, n_max, '-nmax')) return
    case ('-ele')
      ix_word = ix_word + 1
      ele_name = word(ix_word)
    case default
      if (file_name /= '') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
        return
      endif
      file_name = switch
    end select
  enddo

  !

  if (file_name == '') file_name = 'field.dat'

  n = count([dr(1), r_max(1)] /= real_garbage$) + count([n_max(1) /= int_garbage$])
  if (n /= 2) then
    call out_io(s_error$, r_name, 'EXACTLY TWO OF -dr, -rmax, -nmax MUST BE SPECIFIED. NO MORE AND NO LESS.')
    return
  endif

  if (r_max(1) /= real_garbage$ .and. r_min(1) == real_garbage$) r_min = [-r_max(1), -r_max(2), 0.0_rp]
  if (n_max(1) /= int_garbage$ .and. n_min(1) == int_garbage$) n_min = [-n_max(1), -n_max(2), 0]

  if (r_max(1) == real_garbage$) then
    r_max = n_max * dr
    r_min = n_min * dr
  elseif (n_max(1) == int_garbage$) then
    n_max = nint(r_max/dr)
    n_min = nint(r_min/dr)
  else
    dr = (r_max - r_min) / (n_max - n_min)
  endif

  call tao_locate_elements (ele_name, s%global%default_universe, eles, err)
  if (err) return

  if (size(eles) > 1) then
    call out_io(s_warn$, r_name, 'Element name matches to multiple elements: ' // ele_name, 'Will use first one.')
  elseif (size(eles) == 0) then
    call out_io(s_error$, r_name, 'ELEMENT NAME DOES NOT MATCH TO ANY ELEMENTS: ' // ele_name)
    return
  endif

  ele => eles(1)%ele

  !

  open (iu, file = file_name, recl = 200)
  write (iu, '(9a)') '# ele = ', quote(ele%name) 
  write (iu, '(9a)') '# nx = [', int_str(n_min(1)), ', ', int_str(n_max(1)), ']' 
  write (iu, '(9a)') '# ny = [', int_str(n_min(2)), ', ', int_str(n_max(2)), ']' 
  write (iu, '(9a)') '# nz = [', int_str(n_min(3)), ', ', int_str(n_max(3)), ']' 
  write (iu, '(9a)') '# rx = [', real_str(r_min(1)), ', ', real_str(r_max(1)), ']' 
  write (iu, '(9a)') '# ry = [', real_str(r_min(2)), ', ', real_str(r_max(2)), ']' 
  write (iu, '(9a)') '# rz = [', real_str(r_min(3)), ', ', real_str(r_max(3)), ']' 
  write (iu, '(9a)') '# dr = [', real_str(dr(1)), ', ', real_str(dr(2)), ', ', real_str(dr(3)), ']'
  write (iu, '(9a)') '##' 
  write (iu, '(9a)') '##        x           y           z                    Bx                    By                    Bz                    Ex                    Ey                    Ez' 

  branch => pointer_to_branch(ele)

  do ix = n_min(1), n_max(1)
  do iy = n_min(2), n_max(2)
  do iz = n_min(3), n_max(3)
    rr = [ix, iy, iz] * dr
    orbit%vec(1:3:2) = rr(1:2)
    call em_field_calc(ele, branch%param, rr(3), orbit, .false., field)
    write (iu, '(3f12.6, 6es22.13)') rr, field%b, field%e
  enddo
  enddo
  enddo

  close(iu)
  call out_io (s_info$, r_name, 'Written: ' // file_name)


!---------------------------------------------------
! hard

case ('hard', 'hard-l')

  call qp_open_page ('PS', ix, s%plot_page%size(1), s%plot_page%size(2), 'POINTS')
  call tao_draw_plots ()   ! PS out
  call qp_close_page
  call tao_draw_plots ()   ! Update the plotting window

  if (s%global%print_command == ' ') then
    call out_io (s_fatal$, r_name, 'P%PRINT_COMMAND NEEDS TO BE SET TO SEND THE PS FILE TO THE PRINTER!')
    return
  endif

  call system (trim(s%global%print_command) // ' quick_plot.ps')
  call out_io (s_blank$, r_name, 'Printing with command: ' // s%global%print_command)

!---------------------------------------------------
! Foreign lattice format

case ('mad', 'mad8', 'madx', 'opal_latice', 'sad', 'xsif', 'elegant')

  select case (action)
  case ('mad');     file_name0 = 'lat_#.mad8'; lat_type = 'MAD-8'
  case ('mad8');    file_name0 = 'lat_#.mad8'; lat_type = 'MAD-8'
  case ('madx');    file_name0 = 'lat_#.madx'; lat_type = 'MAD-X'
  case ('opal');    file_name0 = 'lat_#.opal'; lat_type = 'OPAL-T'
  case ('xsif');    file_name0 = 'lat_#.xsif'; lat_type = 'XSIF'
  case ('sad');     file_name0 = 'lat_#.sad';  lat_type = 'SAD'
  case ('elegant'); file_name0 = 'lat_#.lte';  lat_type = 'ELEGANT'
  end select

  if (word(1) /= '') file_name0 = word(1) 

  if (word(2) /= '') then
    call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
    return
  endif

  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. tao_subin_uni_number (file_name0, i, file_name)) return
    call write_lattice_in_foreign_format (lat_type, file_name, s%u(i)%model%lat, &
                                             s%u(i)%model%tao_branch(0)%orbit, err = err)
    if (err) return
    call out_io (s_info$, r_name, 'Written: ' // file_name)
  enddo

!---------------------------------------------------
! matrix

case ('matrix')

  ix_word = 0
  file_name = ''
  which = '-single'
  append = .false.
  i_uni = -1
  b_name = ''

  do
    ix_word = ix_word + 1
    if (ix_word == size(word)-1) exit    
    call tao_next_switch (word(ix_word), [character(16):: '-single', '-from_start', '-combined', &
                      '-universe', '-branch'], .true., switch, err)
    if (err) return

    select case (switch)
    case ('')
        exit
    case ('-single', '-from_start', '-combined')
      which = switch
    case ('-universe')
      ix_word = ix_word + 1
      if (.not. is_integer(word(ix_word), i_uni)) then
        call out_io (s_error$, r_name, 'BAD UNIVERSE INDEX: ' // word(ix_word))
        return
      endif
    case ('-branch')
      ix_word = ix_word + 1
      b_name = word(ix_word)
    case default
      if (file_name /= '') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
        return
      endif
      file_name = switch
    end select
  enddo

  !

  u => tao_pointer_to_universe(i_uni)
  if (.not. associated(u)) then
    call out_io (s_error$, r_name, 'BAD UNIVERSE INDEX: ' // word(i_uni))
    return
  endif

  branch => pointer_to_branch(b_name, u%model%lat, blank_branch = s%global%default_branch)
  if (.not. associated(branch)) then
    call out_io (s_error$, r_name, 'BAD LATTICE BRANCH NAME OR INDEX: ' // b_name)
    return
  endif

  if (file_name == '') file_name = 'matrix.dat'
  open (iu, file = file_name)

  call mat_make_unit(mat6)

  do i = 1, branch%n_ele_track
    ele => branch%ele(i)
    mat6 = matmul(ele%mat6, mat6)

    if (which == '-single' .or. which == '-combined') then
      write (iu, *)
      write (iu, '(i6, 2x, a, a16, f16.9)') i, ele%name, key_name(ele%key), ele%s
      call mat_type (ele%mat6, iu, num_form = '(4x, 6f14.8)')
    endif

    if (which == '-from_start' .or. which == '-combined') then
      write (iu, *)
      write (iu, '(a, i6, 2x, a, a16, f16.9)') 'From start to:', i, ele%name, key_name(ele%key), ele%s
      call mat_type (mat6, iu, num_form = '(4x, 6f14.8)')
    endif
  enddo

  close (iu)

!---------------------------------------------------
! namelist

case ('namelist')

  ix_word = 0
  file_name = ''
  which = ''
  append = .false.

  do
    ix_word = ix_word + 1
    if (ix_word == size(word)-1) exit
    call tao_next_switch (word(ix_word), [character(16):: '-data', '-plot', '-variable', '-append'], .true., switch, err)
    if (err) return

    select case (switch)
    case ('');                             exit
    case ('-data', '-plot', '-variable');  which = switch
    case ('-append');                      append = .true.
    case default
      if (file_name /= '') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
        return
      endif
      file_name = switch
    end select
  enddo

  !

  if (which == '') then
    call out_io (s_error$, r_name, 'WHICH NAMELIST (-data, -variable) NOT SET.')
    return
  endif

  if (file_name == '') file_name = 'tao.namelist'

  if (append) then
    open (iu, file = file_name, access = 'append')
  else
    open (iu, file = file_name)
  endif

  !--------------
  ! namelist -data

  select case (which)
  case ('-data')
    do i = 1, size(s%u)
      u => s%u(i)

      do j = 1, u%n_d2_data_used
        d2 => u%d2_data(j)
        write (iu, *)
        write (iu, '(a)') '!---------------------------------------'
        write (iu, *)
        write (iu, '(a)')     '&tao_d2_data'
        write (iu, '(2a)')    '  d2_data%name = ', quote(d2%name)
        write (iu, '(a, i0)') '  universe = ', i
        write (iu, '(a, i0)') '  n_d1_data = ', size(d2%d1)
        write (iu, '(a)')     '/'

        do k = 1, size(d2%d1)
          d1 => d2%d1(k)
          write (iu, *)
          write (iu, '(a)')      '&tao_d1_data'
          write (iu, '(2a)')     '  d1_data%name   = ', quote(d1%name)
          write (iu, '(a, i0)')  '  ix_d1_data     = ', k
          i_min = lbound(d1%d, 1);   i_max = ubound(d1%d, 1)
          write (iu, '(a, i0)')  '  ix_min_data    = ', i_min
          write (iu, '(a, i0)')  '  ix_max_data    = ', i_max

          ! Data output parameter-by-parameter

          len_d_type = 0
          eq_d_type = .true.
          do id = i_min, i_max
            len_d_type = max(len_d_type, len_trim(d1%d(id)%data_type))
            eq_d_type = eq_d_type .and. d1%d(id)%data_type == d1%d(i_min)%data_type
          enddo

          if ((eq_d_type .and. (size(d1%d) > 10)) .or. len_d_type > 30) then
            write_data_source = .true.
            if (all(d1%d%data_source == d1%d(i_min)%data_source)) then
              if (d1%d(i_min)%data_source /= tao_d2_d1_name(d1, .false.)) write (iu, '(2a)') '  default_data_source = ', quote(d1%d(i_min)%data_source)
              write_data_source = .false.
            endif

            write_data_type = .true.
            if (eq_d_type) then
              write (iu, '(2a)') '  default_data_type = ', quote(d1%d(i_min)%data_type)
              write_data_type = .false.
            endif

            write_merit_type = .true.
            if (all(d1%d%merit_type == 'target')) then
              write_merit_type = .false.
            endif

            write_weight = .true.
            if (all(d1%d%weight == d1%d(i_min)%weight)) then
              write (iu, '(2a)') '  default_weight = ', real_to_string(d1%d(i_min)%weight, 12, 5)
              write_weight = .false.
            endif

            if (write_data_source) call namelist_param_out ('d', 'data_source', i_min, i_max, d1%d%data_source)
            if (write_data_type)   call namelist_param_out ('d', 'data_type', i_min, i_max, data_type_arr = d1%d)
            call namelist_param_out ('d', 'ele_name', i_min, i_max, d1%d%ele_name)
            call namelist_param_out ('d', 'ele_start_name', i_min, i_max, d1%d%ele_start_name, '')
            call namelist_param_out ('d', 'ele_ref_name', i_min, i_max, d1%d%ele_ref_name, '')
            if (write_merit_type)  call namelist_param_out ('d', 'merit_type', i_min, i_max, d1%d%merit_type, '')

            if (any(d1%d%good_meas)) then
              call namelist_param_out ('d', 'good_meas', i_min, i_max, logic_arr = d1%d%good_meas, logic_dflt = .false.)
              call namelist_param_out ('d', 'meas', i_min, i_max, re_arr = d1%d%meas_value)
            endif

            if (any(d1%d%good_ref)) then
              call namelist_param_out ('d', 'good_ref', i_min, i_max, logic_arr = d1%d%good_ref, logic_dflt = .false.)
              call namelist_param_out ('d', 'ref', i_min, i_max, re_arr = d1%d%ref_value)
            endif

            if (write_weight)      call namelist_param_out ('d', 'weight', i_min, i_max, re_arr = d1%d%weight)
            call namelist_param_out ('d', 'good_user', i_min, i_max, logic_arr = d1%d%good_user, logic_dflt = .true.)
            call namelist_param_out ('d', 'eval_point', i_min, i_max, anchor_pt_name(d1%d%eval_point), anchor_pt_name(anchor_end$))
            call namelist_param_out ('d', 's_offset', i_min, i_max, re_arr = d1%d%s_offset, re_dflt = 0.0_rp)
            call namelist_param_out ('d', 'ix_bunch', i_min, i_max, int_arr = d1%d%ix_bunch, int_dflt = 0)

          ! Data output datum-by-datum
          else
            n_type   = max(11, len_d_type)
            n_ref    = max(11, maxval(len_trim(d1%d%ele_ref_name)))
            n_start  = max(11, maxval(len_trim(d1%d%ele_start_name)))
            n_ele    = max(11, maxval(len_trim(d1%d%ele_name)))
            n_merit  = max(10, maxval(len_trim(d1%d%merit_type)))
            n_meas   = 14
            n_weight = 12
            n_good   = 6
            n_bunch  = 6
            n_eval   = max(8, maxval(len_trim(anchor_pt_name(d1%d%eval_point))))
            n_s      = 12

            last_col = 'merit'
            if (any(d1%d%meas_value /= 0)) last_col = 'meas'
            if (any(d1%d%weight /= 0)) last_col = 'weight'
            if (any(d1%d%good_user .neqv. .true.)) last_col = 'good'
            if (any(d1%d%ix_bunch /= 0)) last_col = 'bunch'
            if (any(d1%d%eval_point /= anchor_end$)) last_col = 'eval'
            if (any(d1%d%s_offset /= 0)) last_col = 's'

            do m = i_min, i_max
              dat => d1%d(m)
              header1 =                  '  !'
              header2 =                  '  !'
              write (line, '(a, i3, a)') '  datum(', m, ') ='
              n_len = len_trim(line) + 1
              call namelist_item_out (header1, header2, line, n_len, n_type,    'data_', 'type', dat%data_type)
              call namelist_item_out (header1, header2, line, n_len, n_ref,     'ele_ref', 'name', dat%ele_ref_name)
              call namelist_item_out (header1, header2, line, n_len, n_start,   'ele_start', 'name', dat%ele_start_name)
              call namelist_item_out (header1, header2, line, n_len, n_ele,     'ele', 'name', dat%ele_name)
              call namelist_item_out (header1, header2, line, n_len, n_merit,   'merit', 'type', dat%merit_type)
              call namelist_item_out (header1, header2, line, n_len, n_meas,    'meas', 'value', re_val = dat%meas_value)
              call namelist_item_out (header1, header2, line, n_len, n_weight,  'weight', '', re_val = dat%weight)
              call namelist_item_out (header1, header2, line, n_len, n_good ,   'good', 'user', logic_val = dat%good_user)
              call namelist_item_out (header1, header2, line, n_len, n_bunch,   'ix', 'bunch', int_val = dat%ix_bunch)
              call namelist_item_out (header1, header2, line, n_len, n_eval,    'eval', 'point', anchor_pt_name(dat%eval_point))
              call namelist_item_out (header1, header2, line, n_len, n_s,       's', 'offset', re_val = dat%s_offset)
            enddo
          endif

          ! spin out

          do m = i_min, i_max
            if (any(d1%d(m)%spin_map%axis0%n0 /= 0)) &
                    write (iu, '(a, i0, a, 3f12.6)') 'datum(', m, ')%spin_axis%n0 = ', d1%d(m)%spin_map%axis0%n0
          enddo

          write (iu, '(a)') '/'
        enddo

      enddo
    enddo

  !--------------------------------------------
  ! namelist -plot

  case ('-plot')

    write (iu, '(a)') '&tao_plot_page'
    j = 0
    do i = 1, size(s%plot_page%region)
      r => s%plot_page%region(i)
      if (r%plot%name == '' .or. .not. r%visible) cycle
      j = j + 1
      write (iu, '(a, i0, a, i0, 2a)') '  place(', j, ') = @R', r%plot%ix_plot, ', ',  r%plot%name 
    enddo
    write (iu, '(a, /)') '/'

  !--------------------------------------------
  ! namelist -variable

  case ('-variable')

    do i = 1, s%n_v1_var_used
      v1 => s%v1_var(i)
      write (iu, *)
      write (iu, '(a)') '!---------------------------------------'
      write (iu, *)
      write (iu, '(a)')    '&tao_var'
      write (iu, '(2a)')   '  v1_var%name   = ', quote(v1%name)
      i_min = lbound(v1%v, 1);   i_max = ubound(v1%v, 1)
      write (iu, '(a, i0)')  '  ix_min_var    = ', i_min
      write (iu, '(a, i0)')  '  ix_max_var    = ', i_max
      
      if (size(s%u) > 1) then
        call re_allocate2(z, i_min, i_max)
        do j = i_min, i_max
          z(j) = ''
          if (.not. v1%v(j)%exists) cycle 
          s%u%picked_uni = .false.
          do k = 1, size(v1%v(j)%slave)
            s%u(v1%v(j)%slave(k)%ix_uni)%picked_uni = .true.
          enddo
          if (all(s%u%picked_uni)) then
            z(j) = '*'
          else
            z(j) = ''
            do k = lbound(s%u, 1), ubound(s%u, 1)
              if (.not. s%u(k)%picked_uni) cycle
              if (z(j) == '') then
                z(j) = int_str(k)
              else 
                z(j) = trim(z(j)) // ', ' // int_str(k)
              endif
            enddo
          endif
        enddo

        if (all(z == z(i_min))) then
          write (iu, '(2a)') '  default_universe = ', quote(z(i_min))
        else
          call namelist_param_out ('v', 'universe', i_min, i_max, z)
        endif
      endif

      call namelist_param_out ('v', 'ele_name', i_min, i_max, v1%v%ele_name)

      if (all(v1%v%attrib_name == v1%v(i_min)%attrib_name)) then
        write (iu, '(2a)') '  default_attribute = ', quote(v1%v(i_min)%attrib_name)
      else
        call namelist_param_out ('v', 'attribute', i_min, i_max, v1%v%attrib_name)
      endif

      if (all(v1%v%step == v1%v(i_min)%step)) then
        write (iu, '(2a)') '  default_step = ', real_to_string(v1%v(i_min)%step, 12, 5)
      else
        call namelist_param_out ('v', 'step', i_min, i_max, re_arr = v1%v%step)
      endif

      if (all(v1%v%weight == v1%v(i_min)%weight)) then
        write (iu, '(2a)') '  default_weight = ', real_to_string(v1%v(i_min)%weight, 12, 5)
      else
        call namelist_param_out ('v', 'weight', i_min, i_max, re_arr = v1%v%weight)
      endif

      if (all(v1%v%merit_type == v1%v(i_min)%merit_type)) then
        write (iu, '(2a)') '  default_merit_type = ', v1%v(i_min)%merit_type
      else
        call namelist_param_out ('v', 'merit_type', i_min, i_max, v1%v%merit_type)
      endif

      if (all(v1%v%low_lim == v1%v(i_min)%low_lim)) then
        write (iu, '(2a)') '  default_low_lim = ', real_to_string(v1%v(i_min)%low_lim, 12, 5)
      else
        call namelist_param_out ('v', 'low_lim', i_min, i_max, re_arr = v1%v%low_lim, re_dflt = 0.0_rp)
      endif

      if (all(v1%v%high_lim == v1%v(i_min)%high_lim)) then
        write (iu, '(2a)') '  default_high_lim = ', real_to_string(v1%v(i_min)%high_lim, 12, 5)
      else
        call namelist_param_out ('v', 'high_lim', i_min, i_max, re_arr = v1%v%high_lim, re_dflt = 0.0_rp)
      endif

      call namelist_param_out ('v', 'good_user', i_min, i_max, logic_arr = v1%v%good_user, logic_dflt = .true.)
      call namelist_param_out ('v', 'key_bound', i_min, i_max, logic_arr = v1%v%key_bound, logic_dflt = .false.)
      call namelist_param_out ('v', 'key_delta', i_min, i_max, re_arr = v1%v%key_delta, re_dflt = 0.0_rp)
      write (iu, '(a)') '/'
    enddo
  end select

  close (iu)

!---------------------------------------------------
! plot_commands

case ('plot_commands')

  file_name = ''

  do 
    call tao_next_switch (what2, [character(16):: '-XXX'], .true., switch, err)
    if (err) return
    if (switch == '') exit

    select case (switch)
    case ('-XXX')

    case default
      if (file_name /= '') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
        return
      endif
      file_name = switch
    end select
  enddo

  if (file_name == '') file_name = 'plot_commands.tao'
  call fullfilename (file_name, file_name)
  open (iu, file = file_name)

  !

  found_plot_command = .false.

  do i = 1, size(s%history)
    n = modulo(i - 1 + s%com%ix_history, size(s%history)) + 1
    if (n > s%com%ix_history .and. size(s%history) - n + s%com%ix_history >= s%com%n_history) cycle

    call string_trim(s%history(n)%cmd, aline, ix)
    if (aline(1:1) == '!') then ! Is command read from command file.
      if (found_plot_command) write (iu, '(a)') s%history(n)%cmd
      cycle
    endif

    c1 = aline(1:ix)
    call string_trim(aline(ix+1:), aline, ix)
    c2 = aline(1:ix)

    if (c1 == 'set' .and. (index('curve', trim(c2)) == 1 .or. index('floor_plan', trim(c2)) == 1 .or. &
        index('graph', trim(c2)) == 1 .or. index('key', trim(c2)) == 1 .or. index('lat_layout', trim(c2)) == 1 .or. &
        index('plot', trim(c2)) == 1 .or. index('plot_page', trim(c2)) == 1 .or. index('region', trim(c2)) == 1) .or. &
        index('scale', trim(c1)) == 1 .or. index('x_axis', trim(c1)) == 1 .or. index('x_scale', trim(c1)) == 1 .or. &
        index('xy_scale', trim(c1)) == 1) then
      write (iu, '(a)') s%history(n)%cmd
      found_plot_command = .true.
    else
      found_plot_command = .false.
    endif
  enddo

  call out_io (s_info$, r_name, 'Written: ' // file_name)

!---------------------------------------------------
! ps

case ('ps', 'ps-l', 'gif', 'gif-l', 'pdf', 'pdf-l')

  if (qp_base_library == 'PGPLOT' .and. action(1:3) == 'pdf') then
    call out_io (s_error$, r_name, 'PGPLOT DOES NOT SUPPORT PDF!')
    return
  endif

  ix_word = 0
  scale = 1
  file_name = ''

  do
    ix_word = ix_word + 1
    if (ix_word == size(word)-1) exit

    call tao_next_switch (word(ix_word), ['-scale'], .true., switch, err)
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
    case default
      if (file_name /= '') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
        return
      endif
      file_name = switch
    end select
  enddo

  if (word(ix_word) /= '') then
    file_name = word(ix_word)
    if (word(ix_word+1) /= '' .or. file_name(1:1) == '-') then
      call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
      return
    endif
  endif

  if (file_name == '') then
    file_name = "tao.ps"
    if (action(1:3) == 'gif') file_name = 'tao.gif'
    if (action(1:3) == 'pdf') file_name = 'tao.pdf'
  endif

  call str_upcase (action, action)

  call qp_open_page (action, ix, s%plot_page%size(1), s%plot_page%size(2), 'POINTS', file_name, scale)
  call tao_draw_plots (.false.)   ! GIF plot
  call qp_close_page

  call tao_draw_plots ()   ! Update the plotting window

  call out_io (s_blank$, r_name, "Created " // trim(action) // " file: " // file_name)

!---------------------------------------------------
! ptc

case ('ptc')

  which = '-new'
  u => tao_pointer_to_universe(-1)
  branch => u%model%lat%branch(0)
  file_name = ''

  do 
    call tao_next_switch (what2, [character(16):: '-branch', '-all'], .true., switch, err)
    if (err) return
    if (switch == '') exit

    select case (switch)
    case ('-all')
      which = switch
    case ('-branch')
      call tao_next_word(what2, aname)
      branch => pointer_to_branch (aname, u%model%lat, blank_branch = s%global%default_branch)
      if (.not. associated(branch)) then
        call out_io (s_fatal$, r_name, 'Bad branch name or index: ' // aname)
        return
      endif
    case default
      if (file_name /= '') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
        return
      endif
      file_name = switch
    end select
  enddo

  if (file_name == '') file_name = 'ptc.flatfile'

  if (.not. associated(branch%ptc%m_t_layout)) then
    call out_io (s_info$, r_name, 'Note: Creating PTC layout (equivalent to "ptc init").')
    call lat_to_ptc_layout (branch%lat)
  endif

  select case (which)
  case ('-new')
    call print_new_flat (branch%ptc%m_t_layout, file_name)
    call out_io (s_info$, r_name, 'Written: ' // file_name)

  case ('-all')
    call print_universe (M_u, trim(file_name) // '.m_u')
    call print_universe_pointed (M_u, M_t, trim(file_name) // '.m_t')
    call out_io (s_info$, r_name, 'Written: ' // trim(file_name) // '.m_u')
    call out_io (s_info$, r_name, 'Written: ' // trim(file_name) // '.m_t')
  end select

!---------------------------------------------------
! spin

case ('spin_mat8')

  u => tao_pointer_to_universe(-1)
  branch => u%model%lat%branch(0)
  file_name = ''
  sm => datum%spin_map

  do 
    call tao_next_switch (what2, [character(16):: '-l_axis'], .true., switch, err)
    if (err) return
    if (switch == '') exit

    select case (switch)
    case ('-l_axis')
      read (what2, *, iostat = ios) sm%axis_input%l
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'CANNOT PARSE L-AXIS: ' // what2)
        return
      endif
      call word_read(what2, ' ,', str, ix, delim, delim_found, what2) ! Strip Axis from what2
      call word_read(what2, ' ,', str, ix, delim, delim_found, what2)
      call word_read(what2, ' ,', str, ix, delim, delim_found, what2)

    case default
      if (file_name /= '') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
        return
      endif
      file_name = switch
    end select
  enddo

  if (file_name == '') file_name = 'spin_mat8.dat'
  call fullfilename (file_name, file_name)
  open (iu, file = file_name)

  !

  sm%axis_input%n0 = u%model%tao_branch(branch%ix_branch)%orbit(0)%spin
  datum%ix_branch = branch%ix_branch

  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    call tao_spin_matrix_calc (datum, u, pointer_to_next_ele(ele,-1), ele)

    write (iu, *)
    write (iu, '(i6, 2x, a, a16, f16.9)') ie, ele%name, key_name(ele%key), ele%s
    write (iu, '(2(a, 3f14.8))') 'l_start: ', sm%axis0%l,  '  l_end: ', sm%axis1%l
    write (iu, '(2(a, 3f14.8))') 'n0_start:', sm%axis0%n0, '  n0_end:', sm%axis1%n0
    write (iu, '(2(a, 3f14.8))') 'm_start: ', sm%axis0%m,  '  m_end: ', sm%axis1%m
    do i = 1, 8
      write(iu, '(5x, a)') reals_to_table_row(sm%mat8(i,:), 13, 7)
    enddo

    sm%axis_input%n0 = sm%axis1%n0
    sm%axis_input%l  = sm%axis1%l
    sm%axis_input%m  = sm%axis1%m
  enddo

  call out_io (s_info$, r_name, 'Written: ' // file_name)

!---------------------------------------------------
! tao

case ('tao')

  file_name = ''

  do 
    call tao_next_switch (what2, [character(16):: '-XXX'], .true., switch, err)
    if (err) return
    if (switch == '') exit

    select case (switch)
    case ('-XXX')

    case default
      if (file_name /= '') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
        return
      endif
      file_name = switch
    end select
  enddo

  if (file_name == '') file_name = 'tao_new.init'
  call fullfilename (file_name, file_name)
  open (iu, file = file_name)

  ! tao_start namelist

  write (iu, '(a)') '!----------------------------------------------------------'
  write (iu, '(a)') ''
  write (iu, '(a)') '&tao_start'


!---------------------------------------------------
! variable

case ('variable')
  good_opt_only = .false.
  ix_word = 0
  file_name = ''
  tao_format = .false.

  do 
    ix_word = ix_word + 1
    if (ix_word >= size(word)-1) exit
    call tao_next_switch (word(ix_word), [character(20):: '-good_opt_only', '-tao_format'], .true., switch, err)
    if (err) return
    select case (switch)
    case (''); exit
    case ('-tao_format'); tao_format = .true.
    case ('-good_opt_only'); good_opt_only = .true.
    case default
      if (file_name /= '') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON THE COMMAND LINE. NOTHING DONE.')
        return
      endif
      file_name = switch
    end select
  enddo  

  if (file_name == '') then
    call tao_var_write (s%global%var_out_file, good_opt_only, tao_format)
  else
    call tao_var_write (file_name, good_opt_only, tao_format)
  endif

!---------------------------------------------------
! error

case default

  call out_io (s_error$, r_name, 'UNKNOWN "WHAT": ' // what)

end select

!-----------------------------------------------------------------------------
contains

subroutine namelist_item_out (header1, header2, line, n_len, n_add, h1, h2, str_val, re_val, logic_val, int_val)

real(rp), optional :: re_val
integer, optional :: int_val
integer n_len, n_add
logical, optional :: logic_val
character(*) header1, header2, line, h1, h2
character(*), optional :: str_val
character(n_add) add_str

!

header1 = header1(1:n_len) // h1
header2 = header2(1:n_len) // h1

!

if (present(str_val)) then
  add_str = quote(str_val)
elseif (present(re_val)) then
  add_str = real_to_string(re_val, n_add-1, n_add-9)
elseif (present(logic_val)) then
  write (add_str, '(l2)') logic_val
elseif (present(int_val)) then
  write (add_str, '(i4)') int_val
endif

!

line = line(1:n_len) // add_str
n_len = n_len + n_add

end subroutine namelist_item_out

!-----------------------------------------------------------------------------
! contains

subroutine namelist_param_out (who, name, i_min, i_max, str_arr, str_dflt, data_type_arr, re_arr, re_dflt, logic_arr, logic_dflt, int_arr, int_dflt)

integer i_min, i_max

type (tao_data_struct), optional :: data_type_arr(i_min:)
type (var_length_string_struct) :: out_str(i_min:i_max)

real(rp), optional :: re_arr(i_min:), re_dflt

integer i
integer, optional :: int_arr(i_min:), int_dflt
logical, optional :: logic_arr(i_min:), logic_dflt

character(*) who, name
character(*), optional :: str_arr(i_min:), str_dflt
character(600) line


! Encode values

if (present(data_type_arr)) then
  do i = i_min, i_max
    out_str(i)%str = quote(data_type_arr(i)%data_type)
  enddo

elseif (present(str_arr)) then
  if (present(str_dflt)) then
    if (all(str_arr == str_dflt)) return
  endif

  do i = i_min, i_max
    out_str(i)%str = quote(str_arr(i))
  enddo

elseif (present(re_arr)) then
  if (present(re_dflt)) then
    if (all(re_arr == re_dflt)) return
  endif

  do i = i_min, i_max
    out_str(i)%str = real_to_string(re_arr(i), 15, 8)
  enddo

elseif (present(logic_arr)) then
  if (present(logic_dflt)) then
    if (all(logic_arr .eqv. logic_dflt)) return
  endif

  do i = i_min, i_max
    write (out_str(i)%str, '(l1)') logic_arr(i)
  enddo

elseif (present(int_arr)) then
  if (present(int_dflt)) then
    if (all(int_arr == int_dflt)) return
  endif

  do i = i_min, i_max
    write (out_str(i)%str, '(i0)') int_arr(i) 
  enddo
endif

! Write to output
! Note: Using an array multiplyer is not valid for strings.

if (who == 'd') then
  write (line, '(2x, 2(a, i0), 4a)') 'datum(', i_min, ':', i_max, ')%', trim(name), ' = '
else
  write (line, '(2x, 2(a, i0), 4a)') 'var(', i_min, ':', i_max, ')%', trim(name), ' = '
endif

if (all_equal_var_str(out_str, out_str(i_min)%str)) then
  if (present(str_arr)) then
    write (iu, '(a, i0, 2a)') trim(line), i_max-i_min+1, '*', quote(out_str(i_min)%str)
  else
    write (iu, '(a, i0, 2a)') trim(line), i_max-i_min+1, '*', trim(out_str(i_min)%str)
  endif
  return
endif

write (iu, '(a)') trim(line)
line = ''

do i = i_min, i_max
  if (line == '') then
    line = out_str(i)%str
  else
    line = trim(line) // ', ' // out_str(i)%str
  endif

  if (i == i_max) then
    write (iu, '(6x, a)') trim(line)
    exit
  elseif (len_trim(line) +len_trim(out_str(i+1)%str) > 100) then
    write (iu, '(6x, a)') trim(line)
    line = ''
  endif
enddo

end subroutine namelist_param_out

!-----------------------------------------------------------------------------
! contains

function read_real3(word, ix_word, rvec, err_str) result (ok)

character(*) word(:), err_str
real(rp) rvec(:)
integer ix_word
logical ok

!

ok = .false.

do i = 1, 3
  ix_word = ix_word + 1
  if (word(ix_word) == ',') ix_word = ix_word + 1
  read (word(ix_word), *, iostat = ios) rvec(i)
  if (ios /= 0) then
    call out_io(s_error$, r_name, 'ERROR READING ' // err_str // ' VALUE: ' // word(ix_word))
    return
  endif
enddo

ok = .true.

end function read_real3

!-----------------------------------------------------------------------------
! contains

function read_int3(word, ix_word, ivec, err_str) result (ok)

character(*) word(:), err_str
integer ivec(:)
integer ix_word
logical ok

!

ok = .false.

do i = 1, 3
  ix_word = ix_word + 1
  if (word(ix_word) == ',') ix_word = ix_word + 1
  read (word(ix_word), *, iostat = ios) ivec(i)
  if (ios /= 0) then
    call out_io(s_error$, r_name, 'ERROR READING ' // err_str // ' VALUE: ' // word(ix_word))
    return
  endif
enddo

ok = .true.

end function read_int3

end subroutine tao_write_cmd
