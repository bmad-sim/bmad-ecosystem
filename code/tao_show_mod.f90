module tao_show_mod

use tao_mod
use tao_evaluate_mod
use tao_top10_mod
use tao_command_mod, only: tao_cmd_split
use random_mod
use csr_mod, only: csr_param
use tao_lattice_calc_mod

type show_common_struct
  type (ele_struct), pointer :: ele 
  type (coord_struct), pointer :: orbit 
  type (bunch_params_struct), pointer :: bunch_params
  type (tao_universe_struct), pointer :: u
  integer ix_ele
end type

type (show_common_struct), private, save :: show_common
integer, parameter, private :: n_char = 500


contains

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine tao_show_cmd (what, stuff)
!
! Show information on variable, parameters, elements, etc.
!
! Input:
!   what  -- Character(*): What to show.
!   stuff -- Character(*): Particular stuff to show.
!-

subroutine tao_show_cmd (what, stuff)

implicit none


integer iu, ix, n
integer :: n_write_file = 0            ! used for indexing 'show write' files

character(*) what, stuff
character(100) file_name
character(len(stuff)) stuff2
character(16) :: r_name = 'tao_show_cmd'

! See if the results need to be written to a file.

n = len_trim(what)
if (n > 1 .and. (index('-append', trim(what)) == 1 .or. &
                                index('-write', trim(what)) == 1)) then

  call string_trim(stuff, stuff2, ix)
  file_name = stuff2(:ix)
  call string_trim(stuff2(ix+1:), stuff2, ix)

  ix = index(file_name, '*')
  if (ix /= 0) then
    n_write_file = n_write_file + 1
    write (file_name, '(a, i3.3, a)') file_name(1:ix-1), n_write_file, trim(file_name(ix+1:))
  endif

  iu = lunget()
  if (what(1:2) == '-a') then
    open (iu, file = file_name, position = 'APPEND', status = 'UNKNOWN', recl = n_char)
  else
    open (iu, file = file_name, status = 'REPLACE', recl = n_char)
  endif

  call output_direct (iu)  ! tell out_io to write to a file

  call string_trim (stuff2, stuff2, ix)
  call tao_show_this (stuff2(1:ix), stuff2(ix+1:))  

  call output_direct (0)  ! reset to not write to a file
  close (iu)
  call out_io (s_blank$, r_name, 'Written to file: ' // file_name)

else
  call tao_show_this (what, stuff)
endif

end subroutine

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine tao_show_this (what, stuff)

implicit none

type (tao_universe_struct), pointer :: u
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_ptr
type (tao_data_struct), pointer :: d_ptr
type (tao_v1_var_array_struct), allocatable, save, target :: v1_array(:)
type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_var_struct), pointer :: v_ptr
type (tao_var_array_struct), allocatable, save, target :: v_array(:)
type (tao_plot_array_struct), allocatable, save :: plot(:)
type (tao_graph_array_struct), allocatable, save :: graph(:)
type (tao_curve_array_struct), allocatable, save :: curve(:)
type (tao_plot_struct), pointer :: p
type (tao_graph_struct), pointer :: g
type (tao_curve_struct), pointer :: c
type (tao_plot_region_struct), pointer :: region
type (tao_d1_data_array_struct), allocatable, save :: d1_array(:)
type (tao_data_array_struct), allocatable, save :: d_array(:)
type (tao_ele_shape_struct), pointer :: shape
type (branch_struct), pointer :: branch
type (beam_struct), pointer :: beam
type (lat_struct), pointer :: lat
type (bunch_struct), pointer :: bunch
type (lr_wake_struct), pointer :: lr
type (ele_struct), pointer :: ele
type (coord_struct), target :: orb
type (ele_struct), target :: ele3
type (bunch_params_struct) bunch_params
type (bunch_params_struct), pointer :: bunch_p

type show_lat_column_struct
  character(80) name
  character(16) format
  integer field_width
  character(32) label
end type

type (show_lat_column_struct) column(40)

real(rp) f_phi, s_pos, l_lat, gam, s_ele, s1, s2
real(rp) :: delta_e = 0
real(rp), allocatable, save :: value(:)

character(*) :: what, stuff
character(24) :: var_name
character(24)  :: plane, imt, lmt, amt, iamt, f3mt, rmt, irmt, iimt
character(80) :: word(2), fmt, fmt2, fmt3
character(20) :: r_name = "tao_show_cmd"
character(24) show_name, show2_name
character(100), pointer :: ptr_lines(:)
character(100) file_name
character(40) ele_name, name, sub_name
character(60) nam

character(16) :: show_what, show_names(22) = (/ &
   'data        ', 'variable    ', 'global      ', 'alias       ', 'top10       ', &
   'optimizer   ', 'element     ', 'lattice     ', 'constraints ', 'plot        ', &
   'beam        ', 'tune        ', 'graph       ', 'curve       ', 'particle    ', &
   'hom         ', 'opt_vars    ', 'universe    ', 'orbit       ', 'derivative  ', &
   'branches    ', 'use         ' /)

character(n_char), allocatable, save :: lines(:)
character(n_char) line, line1, line2, line3
character(n_char) stuff2
character(9) angle

integer :: data_number, ix_plane, ix_class, n_live, n_tot
integer nl, loc, ixl, iu, nc, n_size, ix_u, ios, ie, nb, id, iv, jd, jv
integer ix, ix1, ix2, ix_s2, i, j, k, n, show_index, ju, ios1, ios2, i_uni
integer num_locations, ix_ele, n_name, n_e0, n_e1
integer, allocatable, save :: ix_eles(:)

logical err, found, at_ends, first_time, by_s, print_header_lines
logical show_sym, show_line, show_shape, print_data, ok
logical show_all, name_found, print_taylor, print_wig_terms, print_all
logical, allocatable, save :: picked_uni(:)
logical, allocatable, save :: picked_ele(:)
logical, allocatable, save :: good(:)

namelist / custom_show_list / column

!

call re_allocate (ix_eles, 1)
call re_allocate (lines, n_char, 200)

err = .false.

lines = " "
nl = 0

rmt  = '(a, 9es16.8)'
f3mt  = '(a, 9f0.3)'
irmt = '(a, i0, a, es16.8)'
imt  = '(a, 9i8)'
iimt = '(a, i0, a, i8)'
lmt  = '(a, 9l)'
amt  = '(9a)'
iamt = '(a, i0, 9a)'

u => tao_pointer_to_universe(-1)
lat => u%model%lat

if (s%global%phase_units == radians$) f_phi = 1
if (s%global%phase_units == degrees$) f_phi = 180 / pi
if (s%global%phase_units == cycles$)  f_phi = 1 / twopi

! find what to show

if (what == ' ') then
  call out_io (s_error$, r_name, 'SHOW WHAT?')
  return
endif

call match_word (what, show_names, ix, matched_name = show_what)
if (ix == 0) then
  call out_io (s_error$, r_name, 'SHOW WHAT? WORD NOT RECOGNIZED: ' // what)
  return
endif

if (ix < 0) then
  call out_io (s_error$, r_name, 'SHOW WHAT? AMBIGUOUS: ' // what)
  return
endif

call tao_cmd_split (stuff, 2, word, .false., err)

select case (show_what)

!----------------------------------------------------------------------
! alias

case ('alias')

  call re_allocate (lines, len(lines(1)), tao_com%n_alias+10, .false.)
  lines(1) = 'Aliases:'
  nl = 1
  do i = 1, tao_com%n_alias
    nl=nl+1; lines(nl) = trim(tao_com%alias(i)%name) // ' = "' // &
                                    trim(tao_com%alias(i)%string) // '"'
  enddo
  
  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! beam

case ('beam')

  ! no element index

  if (word(1) == '') then

    nl=nl+1; write(lines(nl), '(a, i3)') 'Universe: ', u%ix_uni
    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), rmt) 'beam_init%a_norm_emitt      = ', u%beam_init%a_norm_emitt
    nl=nl+1; write(lines(nl), rmt) 'beam_init%b_norm_emitt      = ', u%beam_init%b_norm_emitt
    nl=nl+1; write(lines(nl), rmt) 'beam_init%dPz_dz            = ', u%beam_init%dPz_dz
    nl=nl+1; write(lines(nl), rmt) 'beam_init%center            = ', u%beam_init%center
    nl=nl+1; write(lines(nl), rmt) 'beam_init%ds_bunch          = ', u%beam_init%ds_bunch
    nl=nl+1; write(lines(nl), rmt) 'beam_init%sig_z             = ', u%beam_init%sig_z
    nl=nl+1; write(lines(nl), rmt) 'beam_init%sig_e             = ', u%beam_init%sig_e
    nl=nl+1; write(lines(nl), rmt) 'beam_init%bunch_charge      = ', u%beam_init%bunch_charge
    nl=nl+1; write(lines(nl), rmt) 'beam_init%center_jitter     = ', u%beam_init%center_jitter
    nl=nl+1; write(lines(nl), rmt) 'beam_init%emitt_jitter      = ', u%beam_init%emitt_jitter
    nl=nl+1; write(lines(nl), rmt) 'beam_init%sig_z_jitter      = ', u%beam_init%sig_z_jitter
    nl=nl+1; write(lines(nl), rmt) 'beam_init%sig_e_jitter      = ', u%beam_init%sig_e_jitter
    nl=nl+1; write(lines(nl), rmt) 'beam_init%spin%polarization = ', u%beam_init%spin%polarization
    nl=nl+1; write(lines(nl), rmt) 'beam_init%spin%theta        = ', u%beam_init%spin%theta
    nl=nl+1; write(lines(nl), rmt) 'beam_init%spin%phi          = ', u%beam_init%spin%phi
    nl=nl+1; write(lines(nl), imt) 'beam_init%n_particle        = ', u%beam_init%n_particle
    nl=nl+1; write(lines(nl), imt) 'beam_init%n_bunch           = ', u%beam_init%n_bunch
    nl=nl+1; write(lines(nl), lmt) 'beam_init%renorm_center     = ', u%beam_init%renorm_center
    nl=nl+1; write(lines(nl), lmt) 'beam_init%renorm_sigma      = ', u%beam_init%renorm_sigma
    nl=nl+1; write(lines(nl), lmt) 'beam_init%preserve_dist     = ', u%beam_init%preserve_dist
    nl=nl+1; write(lines(nl), lmt) 'beam_init%init_spin         = ', u%beam_init%init_spin
    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), lmt) 'bmad_com%sr_wakes_on               = ', bmad_com%sr_wakes_on
    nl=nl+1; write(lines(nl), lmt) 'bmad_com%lr_wakes_on               = ', bmad_com%lr_wakes_on
    nl=nl+1; write(lines(nl), lmt) 'bmad_com%trans_space_charge_on     = ', bmad_com%trans_space_charge_on
    nl=nl+1; write(lines(nl), lmt) 'bmad_com%coherent_synch_rad_on     = ', bmad_com%coherent_synch_rad_on
    nl=nl+1; write(lines(nl), lmt) 'bmad_com%spin_tracking_on          = ', bmad_com%spin_tracking_on
    nl=nl+1; write(lines(nl), lmt) 'bmad_com%radiation_damping_on      = ', bmad_com%radiation_damping_on
    nl=nl+1; write(lines(nl), lmt) 'bmad_com%radiation_fluctuations_on = ', bmad_com%radiation_fluctuations_on
    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), rmt) 'csr_param%ds_track_step        = ', csr_param%ds_track_step
    nl=nl+1; write(lines(nl), imt) 'csr_param%n_bin                = ', csr_param%n_bin
    nl=nl+1; write(lines(nl), imt) 'csr_param%particle_bin_span    = ', csr_param%particle_bin_span
    nl=nl+1; write(lines(nl), lmt) 'csr_param%lcsr_component_on    = ', csr_param%lcsr_component_on
    nl=nl+1; write(lines(nl), lmt) 'csr_param%lsc_component_on     = ', csr_param%lsc_component_on
    nl=nl+1; write(lines(nl), lmt) 'csr_param%tsc_component_on     = ', csr_param%tsc_component_on
    nl=nl+1; write(lines(nl), lmt) 'csr_param%ix1_ele_csr          = ', csr_param%ix1_ele_csr
    nl=nl+1; write(lines(nl), lmt) 'csr_param%ix2_ele_csr          = ', csr_param%ix2_ele_csr
    nl=nl+1; lines(nl) = ''
    call convert_total_energy_to (lat%ele(0)%value(e_tot$), lat%param%particle, gamma = gam)
    nl=nl+1; write(lines(nl), rmt) 'model%lat%a%emit               = ', lat%a%emit
    nl=nl+1; write(lines(nl), rmt) '          a%emit (normalized)  = ', lat%a%emit * gam
    nl=nl+1; write(lines(nl), rmt) 'model%lat%b%emit               = ', lat%b%emit
    nl=nl+1; write(lines(nl), rmt) '          b%emit (normalized)  = ', lat%b%emit * gam
    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), amt)  'global%track_type           = ', s%global%track_type
    nl=nl+1; write (lines(nl), imt) 'u%ix_track_start           = ', u%ix_track_start 
    nl=nl+1; write (lines(nl), imt) 'u%ix_track_end             = ', u%ix_track_end
    nl=nl+1; write(lines(nl), amt)  'u%save_beam_at:'
    do i = lbound(u%save_beam_at, 1), ubound(u%save_beam_at, 1)
      nl=nl+1; write (lines(nl), '(a, i0, 2a)') '           (', i, ') = ', u%save_beam_at(i)
    enddo

  ! have element index

  else
    call tao_to_int (word(1), ix_ele, err)
    if (err) return
    n = s%global%bunch_to_plot
    bunch_p => u%model%bunch_params(ix_ele)
    nl=nl+1; lines(nl) = 'Cashed bunch parameters:'
    nl=nl+1; write (lines(nl), rmt) '  Centroid:', bunch_p%centroid%vec
    nl=nl+1; write (lines(nl), rmt) '  RMS:     ', &
                              sqrt(bunch_p%sigma((/s11$, s22$, s33$, s44$, s55$, s66$/)))
    nl=nl+1; write (lines(nl), rmt) '             norm_emitt           beta'
    nl=nl+1; write (lines(nl), rmt) '  a:       ', bunch_p%a%norm_emitt, bunch_p%a%beta
    nl=nl+1; write (lines(nl), rmt) '  b:       ', bunch_p%b%norm_emitt, bunch_p%b%beta
    nl=nl+1; write (lines(nl), rmt) '  x:       ', bunch_p%x%norm_emitt, bunch_p%x%beta
    nl=nl+1; write (lines(nl), rmt) '  y:       ', bunch_p%y%norm_emitt, bunch_p%y%beta
    nl=nl+1; write (lines(nl), rmt) '  z:       ', bunch_p%z%norm_emitt, bunch_p%z%beta

    beam => u%ele(ix_ele)%beam
    if (allocated(beam%bunch)) then
      bunch => beam%bunch(n)
      call calc_bunch_params (bunch, lat%ele(ix_ele), bunch_params, err)
      n_live = bunch_params%n_live_particle
      n_tot = size(bunch%particle)

      nl=nl+1; lines(nl) = 'Parameters from saved beam at element:'
      nl=nl+1; write (lines(nl), imt)  '  Parameters for bunch:       ', n
      nl=nl+1; write (lines(nl), imt)  '  Particles surviving:        ', n_live
      nl=nl+1; write (lines(nl), imt)  '  Particles lost:             ', n_tot - n_live
      nl=nl+1; write (lines(nl), f3mt) '  Particles lost (%):         ', &
                                                  100 * real(n_tot - n_live) / n_tot
      nl=nl+1; write (lines(nl), rmt) '  Centroid:', bunch_params%centroid%vec
      nl=nl+1; write (lines(nl), rmt) '  RMS:     ', &
                         sqrt(bunch_params%sigma((/s11$, s22$, s33$, s44$, s55$, s66$/)))
      nl=nl+1; write (lines(nl), rmt) '             norm_emitt           beta'
      nl=nl+1; write (lines(nl), rmt) '  a:       ', &
                                bunch_params%a%norm_emitt, bunch_params%a%beta
      nl=nl+1; write (lines(nl), rmt) '  b:       ', &
                                bunch_params%b%norm_emitt, bunch_params%b%beta
      nl=nl+1; write (lines(nl), rmt) '  x:       ', &
                                bunch_params%x%norm_emitt, bunch_params%x%beta
      nl=nl+1; write (lines(nl), rmt) '  y:       ', &
                                bunch_params%y%norm_emitt, bunch_params%y%beta
      nl=nl+1; write (lines(nl), rmt) '  z:       ', &
                                bunch_params%z%norm_emitt, bunch_params%z%beta
    else
      nl=nl+1; lines(nl) = 'No Allocated Beam At Element.'
    endif
  
  endif

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! constraints

case ('branches')

  if (ubound(lat%branch, 1) == 0) then
    call out_io (s_blank$, r_name, 'No branches')
    return
  endif

  nl=nl+1; lines(nl) = '                        N_ele  N_ele    From  From  From Ele'
  nl=nl+1; lines(nl) = 'Ix  Name                Track    Max  Branch   Ele  Name'
  do i = 1, ubound(lat%branch, 1)
    branch => lat%branch(i)
    nl=nl+1; write (lines(nl), '(i2, 2x, a21, 1x, i3, i7, i8, i6, 2x, a)') &
                i, branch%name, branch%n_ele_track, &
                branch%n_ele_max, branch%ix_from_branch, branch%ix_from_ele, &
                lat%branch(branch%ix_from_branch)%ele(branch%ix_from_ele)%name
  enddo

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! constraints

case ('constraints')

  call tao_show_constraints (0, 'ALL')
  call tao_show_constraints (0, 'TOP10')
  return

!----------------------------------------------------------------------
! curve

case ('curve')

  ! Look for switches

  show_sym = .false.
  show_line = .false.
  call string_trim(stuff, stuff2, ix)

  do
    if (ix == 0) exit
    if (stuff2(1:1) /= '-') exit
    call match_word (stuff2(:ix), (/ '-symbol', '-line  ' /), n, .true., name)
    if (n < 1) then
      call out_io (s_error$, r_name, 'UNKNOWN SWITCH: ' // stuff2(:ix))
      return
    endif
    if (name == '-symbol') show_sym = .true.
    if (name == '-line')   show_line = .true.
    call string_trim(stuff2(ix+1:), stuff2, ix)
  enddo

  ! Find particular plot

  call tao_find_plots (err, stuff2, 'BOTH', curve = curve, &
                                    print_flag = .false., always_allocate = .true.)
  if (err) return

  ! print info on particular plot, graph, or curve

  if (allocated(curve)) then
    c => curve(1)%c
    nl=nl+1; lines(nl) = 'Region.Graph.Curve: ' // trim(tao_curve_name(c, .true.))
    do i = 2, size(curve)
      nl=nl+1; lines(nl) = '                    ' // trim(tao_curve_name(curve(i)%c, .true.))
    enddo
    nl=nl+1; lines(nl) = 'Plot.Graph.Curve:   ' // trim(tao_curve_name(c))
    do i = 2, size(curve)
      nl=nl+1; lines(nl) = '                    ' // trim(tao_curve_name(curve(i)%c))
    enddo
    nl=nl+1; write (lines(nl), amt) 'data_source          = ', c%data_source
    nl=nl+1; write (lines(nl), amt) 'data_index           = ', c%data_index
    nl=nl+1; write (lines(nl), amt) 'data_type_x          = ', c%data_type_x
    nl=nl+1; write (lines(nl), amt) 'data_type            = ', c%data_type
    nl=nl+1; write (lines(nl), amt) 'ele_ref_name         = ', c%ele_ref_name
    nl=nl+1; write (lines(nl), imt) 'ix_ele_ref           = ', c%ix_ele_ref
    nl=nl+1; write (lines(nl), imt) 'ix_ele_ref_track     = ', c%ix_ele_ref_track
    nl=nl+1; write (lines(nl), imt) 'ix_bunch             = ', c%ix_bunch
    nl=nl+1; write (lines(nl), imt) 'ix_universe          = ', c%ix_universe
    nl=nl+1; write (lines(nl), imt) 'symbol_every         = ', c%symbol_every
    nl=nl+1; write (lines(nl), rmt) 'x_axis_scale_factor  = ', c%x_axis_scale_factor
    nl=nl+1; write (lines(nl), rmt) 'y_axis_scale_factor  = ', c%y_axis_scale_factor
    nl=nl+1; write (lines(nl), lmt) 'use_y2               = ', c%use_y2
    nl=nl+1; write (lines(nl), lmt) 'draw_line            = ', c%draw_line
    nl=nl+1; write (lines(nl), lmt) 'draw_symbols         = ', c%draw_symbols
    nl=nl+1; write (lines(nl), lmt) 'draw_symbol_index    = ', c%draw_symbol_index
    nl=nl+1; write (lines(nl), lmt) 'smooth_line_calc     = ', c%smooth_line_calc
    
    if (show_sym) then
      n = nl + size(c%x_symb) + 10
      if (n > size(lines)) call re_allocate(lines, len(lines(1)), n, .false.)
      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = 'Symbol points:'
      nl=nl+1; lines(nl) = '     i index             x             y'
      err = .false.
      do j = 2, size(curve)
        if (size(curve(j)%c%y_symb) /= size(c%y_symb)) then
          nl=nl+1; lines(nl) = 'NUMBER OF SYMBOL POINTS NOT THE SAME IN ALL CURVES!'
          err = .true.
          exit
        endif
      enddo
      if (.not. err) then
        do i = 1, size(c%x_symb)
          nl=nl+1; write (lines(nl), '(i6, 10es14.6)') i, c%ix_symb(i), c%x_symb(i), &
                                        (/ (curve(j)%c%y_symb(i), j = 1, size(curve)) /)
        enddo
      endif
    endif

    if (show_line) then
      n = nl + size(c%x_line) + 10
      if (n > size(lines)) call re_allocate(lines, len(lines(1)), n, .false.)
      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = 'Smooth line points:'
      nl=nl+1; lines(nl) = '             x             y'
      do j = 2, size(curve)
        if (size(curve(j)%c%y_line) /= size(c%y_line)) then
          nl=nl+1; lines(nl) = 'NUMBER OF LINE POINTS NOT THE SAME IN ALL CURVES!'
          err = .true.
          exit
        endif
      enddo
      if (.not. err) then
        do i = 1, size(c%x_line)
          nl=nl+1; write (lines(nl), '(2es14.6)') c%x_line(i), &
                                        (/ (curve(j)%c%y_line(i), j = 1, size(curve)) /)
        enddo
      endif
    endif

  else
    call out_io (s_error$, r_name, 'This is not a curve name')
    return
  endif

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! data

case ('data')


  ! If just "show data" then show all names

  call tao_pick_universe (word(1), line1, picked_uni, err)
  if (err) return

  if (line1 == ' ') then  ! just specified a universe

    do iu = lbound(s%u, 1), ubound(s%u, 1)

      if (.not. picked_uni(iu)) cycle

      u => s%u(iu)

      nl=nl+1; lines(nl) = ''
      nl=nl+1; write(lines(nl), '(t40, a)') 'Using'

      do i = 1, size(u%d2_data)
        d2_ptr => u%d2_data(i)
        if (d2_ptr%name == ' ') cycle
        nl=nl+1; lines(nl) = ' '
        do j = lbound(d2_ptr%d1, 1), ubound(d2_ptr%d1, 1)
          d1_ptr => d2_ptr%d1(j)
          call location_encode(line, d1_ptr%d%useit_opt, &
                                  d1_ptr%d%exists, lbound(d1_ptr%d, 1))
          nl=nl+1; write (lines(nl), '(2a, i0, a, i0, a, t40, a)') &
                  trim(tao_d2_d1_name(d1_ptr)), '[', &
                  lbound(d1_ptr%d, 1), ':', ubound(d1_ptr%d, 1), ']', trim(line)
        enddo
      enddo
    enddo

    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif

  ! get pointers to the data

  call tao_find_data (err, word(1), d2_ptr, d1_array, d_array)
  if (err) return

  n_size = 0
  if (allocated(d_array)) n_size = size(d_array)

  ! If d_ptr points to something then show the datum info.

  if (n_size == 1) then
    d_ptr => d_array(1)%d
    nl=nl+1; lines(nl) = ''
    if (size(s%u) > 1) then
      nl=nl+1; write(lines(nl), '(a, i4)') 'Universe:', d_ptr%d1%d2%ix_uni
    endif
    nl=nl+1; write(lines(nl), amt)  '%ele0_name         = ', d_ptr%ele0_name
    nl=nl+1; write(lines(nl), amt)  '%ele_name          = ', d_ptr%ele_name
    nl=nl+1; write(lines(nl), amt)  '%data_type         = ', d_ptr%data_type
    nl=nl+1; write(lines(nl), amt)  '%data_source       = ', d_ptr%data_source
    nl=nl+1; write(lines(nl), imt)  '%ix_ele0           = ', d_ptr%ix_ele0
    nl=nl+1; write(lines(nl), imt)  '%ix_ele            = ', d_ptr%ix_ele
    nl=nl+1; write(lines(nl), imt)  '%ix_ele_merit      = ', d_ptr%ix_ele_merit
    nl=nl+1; write(lines(nl), imt)  '%ix_dmodel         = ', d_ptr%ix_dModel
    nl=nl+1; write(lines(nl), imt)  '%ix_d1             = ', d_ptr%ix_d1
    nl=nl+1; write(lines(nl), imt)  '%ix_data           = ', d_ptr%ix_data
    nl=nl+1; write(lines(nl), imt)  '%ix_bunch          = ', d_ptr%ix_bunch
    nl=nl+1; write(lines(nl), rmt)  '%model             = ', d_ptr%model_value
    nl=nl+1; write(lines(nl), rmt)  '%design            = ', d_ptr%design_value
    nl=nl+1; write(lines(nl), rmt)  '%meas              = ', d_ptr%meas_value
    nl=nl+1; write(lines(nl), rmt)  '%ref               = ', d_ptr%ref_value
    nl=nl+1; write(lines(nl), rmt)  '%base              = ', d_ptr%base_value
    nl=nl+1; write(lines(nl), rmt)  '%old               = ', d_ptr%old_value   
    nl=nl+1; write(lines(nl), rmt)  '%fit               = ', d_ptr%fit_value
    nl=nl+1; write(lines(nl), rmt)  '%s                 = ', d_ptr%s
    nl=nl+1; write(lines(nl), amt)  '%merit_type        = ', d_ptr%merit_type
    nl=nl+1; write(lines(nl), rmt)  '%merit             = ', d_ptr%merit
    nl=nl+1; write(lines(nl), rmt)  '%delta_merit       = ', d_ptr%delta_merit
    nl=nl+1; write(lines(nl), rmt)  '%weight            = ', d_ptr%weight
    nl=nl+1; write(lines(nl), lmt)  '%relative          = ', d_ptr%relative
    nl=nl+1; write(lines(nl), lmt)  '%exists            = ', d_ptr%exists
    nl=nl+1; write(lines(nl), lmt)  '%good_model        = ', d_ptr%good_model
    nl=nl+1; write(lines(nl), lmt)  '%good_meas         = ', d_ptr%good_meas
    nl=nl+1; write(lines(nl), lmt)  '%good_ref          = ', d_ptr%good_ref
    nl=nl+1; write(lines(nl), lmt)  '%good_user         = ', d_ptr%good_user
    nl=nl+1; write(lines(nl), lmt)  '%good_opt          = ', d_ptr%good_opt
    nl=nl+1; write(lines(nl), lmt)  '%good_plot         = ', d_ptr%good_plot
    nl=nl+1; write(lines(nl), lmt)  '%useit_plot        = ', d_ptr%useit_plot
    nl=nl+1; write(lines(nl), lmt)  '%useit_opt         = ', d_ptr%useit_opt

  ! Else show the d1_data info.

  elseif (size(d1_array) == 1) then

    d1_ptr => d1_array(1)%d1
    if (size(s%u) > 1) then
      nl=nl+1; write(lines(nl), '(a, i4)') 'Universe:', d1_ptr%d2%ix_uni
    endif
    
    nl=nl+1; write(lines(nl), '(2a)') 'Data name: ', trim(d2_ptr%name) // '.' // d1_ptr%name

    ! find string widths

    n_name = 9
    n_e0 = 8
    n_e1 = 8
    do i = 1, size(d_array)
      d_ptr => d_array(i)%d
      if (.not. d_ptr%exists) cycle
      n_name = max(n_name, len_trim(tao_datum_type_name(d_ptr)))
      n_e0 = max(n_e0, len_trim(d_ptr%ele0_name))
      n_e1 = max(n_e1, len_trim(d_ptr%ele_name))
    enddo

    ! Write header

    line1 = ''; line2 = ''
    n=9+n_name;               line2(n:) = 'Where0'
    n=len_trim(line2)+n_e0-3; line2(n:) = 'Where'
    n=len_trim(line2)+n_e1+4; line2(n:) = 'Meas         Model        Design    | Opt  Plot'
                              line1(n:) = '                                    |   Useit'

    nl=nl+1; lines(nl) = line1
    nl=nl+1; lines(nl) = line2

    ! if a range is specified, show the data range   

    call re_allocate (lines, len(lines(1)), nl+100+size(d_array), .false.)

    fmt = '(i4, 2x, a, 2x, a, 2x, a, 3es14.4, 2l6)'

    do i = 1, size(d_array)
      d_ptr => d_array(i)%d
      if (.not. d_ptr%exists) cycle
      name = tao_datum_type_name(d_ptr)
      nl=nl+1; write(lines(nl), fmt) d_ptr%ix_d1, name(1:n_name), & 
                     d_ptr%ele0_name(1:n_e0), d_ptr%ele_name(1:n_e1), &
                     d_ptr%meas_value, d_ptr%model_value, &
                     d_ptr%design_value, d_ptr%useit_opt, d_ptr%useit_plot
    enddo

    nl=nl+1; lines(nl) = line2
    nl=nl+1; lines(nl) = line1

  ! else we must have a valid d2_ptr.

  elseif (associated(d2_ptr)) then

    call re_allocate (lines, len(lines(1)), nl+100+size(d2_ptr%d1), .false.)

    nl=nl+1; write(lines(nl), '(t40, a)')     'Using' 

    do i = 1, size(d2_ptr%d1)
      if (size(lines) > nl + 50) call re_allocate (lines, len(lines(1)), nl+100, .false.)
      call location_encode(line, d2_ptr%d1(i)%d%useit_opt, &
                      d2_ptr%d1(i)%d%exists, lbound(d2_ptr%d1(i)%d, 1))
      nl=nl+1; write(lines(nl), '(2x, 2a, i0, a, i0, a, t40, a)') &
                  trim(tao_d2_d1_name(d2_ptr%d1(i))), '[', lbound(d2_ptr%d1(i)%d, 1), &
                  ':', ubound(d2_ptr%d1(i)%d, 1), ']', trim(line)
    enddo

    if (any(d2_ptr%descrip /= ' ')) then
      call re_allocate (lines, len(lines(1)), nl+100+size(d2_ptr%descrip), .false.)
      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = 'Descrip:'
      do i = 1, size(d2_ptr%descrip)
        if (d2_ptr%descrip(i) /= ' ') then
          nl=nl+1; write (lines(nl), '(i4, 2a)') i, ': ', d2_ptr%descrip(i)
        endif
      enddo
    endif

  ! error

  else
    lines(1) = 'TRY BEING MORE SPECIFIC.'
    nl = 1
  endif

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! derivative

case ('derivative')

  call tao_find_data (err, word(1), d_array = d_array)
  if (err) return

  call tao_find_var(err, word(2), v_array = v_array) 
  if (err) return

  do id = 1, size(d_array)
    do iv = 1, size(v_array)
      d_ptr => d_array(id)%d
      v_ptr => v_array(iv)%v
      u => s%u(d_ptr%d1%d2%ix_uni)
      jd = d_ptr%ix_dmodel
      jv = v_ptr%ix_dvar
      if (jd == 0 .or. jv == 0) then
        nl=nl+1; write (lines(nl), '(2a20, a14, 2i5)') tao_datum_name(d_ptr), &
                          tao_var1_name(v_ptr), ' No Derivative', jd, jv
      else
        nl=nl+1; write (lines(nl), '(2a20, es14.5, 2i5)') tao_datum_name(d_ptr), &
                                  tao_var1_name(v_ptr), u%dModel_dVar(jd, jv), jd, jv
      endif
    enddo
  enddo

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! ele

case ('element')

  print_taylor = .false.
  print_wig_terms = .false.
  print_all = .false.
  print_data = .false.

  call string_trim(stuff, stuff2, ix)
  do
    if (ix == 0) exit
    if (stuff2(1:1) /= '-') exit
    call match_word (stuff2(:ix), (/ '-taylor        ', '-wig_terms     ', &
                      '-all_attributes', '-data          ' /), n, .true., name)
    if (n < 1) then
      call out_io (s_error$, r_name, 'UNKNOWN SWITCH: ' // stuff2(:ix))
      return
    endif
    if (name == '-taylor') print_taylor = .true.
    if (name == '-wig_term') print_wig_terms = .true.
    if (name == '-all_attributes') print_all = .true.
    if (name == '-data') print_data = .true.
    call string_trim(stuff2(ix+1:), stuff2, ix)
  enddo

  call str_upcase (ele_name, stuff2)
  call tao_pick_universe (ele_name, ele_name, picked_uni, err, ix_u)
  if (err) return

  ! Wildcard: show all elements.

  if (index(ele_name, '*') /= 0 .or. index(ele_name, '%') /= 0 .or. &
                (ele_name(1:2) /= 'S:' .and. index(ele_name, ':') /= 0) .or. &
                count(picked_uni) > 1) then

    nl=1; write (lines(1), '(a, i0)') 'Matches: ', size(ix_eles)

    do i_uni = lbound(s%u, 1), ubound(s%u, 1)
      if (.not. picked_uni(i_uni)) cycle
      call tao_ele_locations_given_name (u, ele_name, ix_eles, err, .true.)
      if (err) return
      do i = 1, size(ix_eles)
        loc = ix_eles(i)
        if (size(lines) < nl+100) call re_allocate (lines, len(lines(1)), nl+200, .false.)
        if (count(picked_uni) > 1) then
          nl=nl+1; write (lines(nl), '(i0, a, i8, 2x, a)') &
                                                i_uni, '@', loc, lat%ele(loc)%name
        else
          nl=nl+1; write (lines(nl), '(i8, 2x, a)') loc, lat%ele(loc)%name
        endif
      enddo
    enddo

    deallocate(ix_eles)

    if (nl == 1) then
      call out_io (s_blank$, r_name, '*** No Matches to Name Found ***')
      return
    endif

    call out_io (s_blank$, r_name, lines(1:nl))
    return

  endif

  ! No wildcard case...
  ! Normal: Show the element info

  call tao_locate_elements (ele_name, ix_u, ix_eles)
  loc = ix_eles(1)
  if (loc < 0) return
  ele => lat%ele(loc)

  ! Show data associated with this element
  if (print_data) then
    call show_ele_data (u, loc, lines, nl)
    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif

  if (tao_com%common_lattice) then
    call tao_lattice_calc (ok, ix_u, model$)
  endif
  call type2_ele (ele, ptr_lines, n, print_all, 6, print_taylor, s%global%phase_units, &
            .true., s%u(ix_u)%model%lat, .true., .true., print_wig_terms)

  if (size(lines) < nl+n+100) call re_allocate (lines, len(lines(1)), nl+n+100, .false.)
  lines(nl+1:nl+n) = ptr_lines(1:n)
  nl = nl + n
  deallocate (ptr_lines)

  if (show_what == 'element') then
    nl=nl+1; lines(nl) = '[Conversion from Global to Screen: (Z, X) -> (-X, -Y)]'
  endif

  orb = u%model%orb(loc)
  fmt = '(2x, a, 3p2f15.8)'
  lines(nl+1) = ' '
  lines(nl+2) = 'Orbit: [mm, mrad]'
  write (lines(nl+3), fmt) "X  X':", orb%vec(1:2)
  write (lines(nl+4), fmt) "Y  Y':", orb%vec(3:4)
  write (lines(nl+5), fmt) "Z  Z':", orb%vec(5:6)
  nl = nl + 5

  found = .false.
  do i = loc + 1, lat%n_ele_max
    if (lat%ele(i)%name /= ele_name) cycle
    if (size(lines) < nl+2) call re_allocate (lines, len(lines(1)), nl+10, .false.)
    if (found) then
      nl=nl+1; lines(nl) = ''
      found = .true.
    endif 
    nl=nl+1;  write (lines(nl), '(a, i0)') 'Note: Found another element with same name at: ', i
  enddo

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! global

case ('global')

  nl=nl+1; lines(nl) = 'Global:'
  nl=nl+1; write (lines(nl), lmt) '%auto_scale                 = ', s%global%auto_scale
  nl=nl+1; write (lines(nl), imt) '%bunch_to_plot              = ', s%global%bunch_to_plot
  nl=nl+1; write (lines(nl), rmt) '%de_lm_step_ratio           = ', s%global%de_lm_step_ratio
  nl=nl+1; write (lines(nl), lmt) '%derivative_recalc          = ', s%global%derivative_recalc
  nl=nl+1; write (lines(nl), lmt) '%label_lattice_elements     = ', s%global%label_lattice_elements
  nl=nl+1; write (lines(nl), lmt) '%label_keys                 = ', s%global%label_keys
  nl=nl+1; write (lines(nl), rmt) '%lm_opt_deriv_reinit        = ', s%global%lm_opt_deriv_reinit
  nl=nl+1; write (lines(nl), rmt) '%lmdif_eps                  = ', s%global%lmdif_eps
  nl=nl+1; write (lines(nl), rmt) '%merit_finish               = ', s%global%merit_finish
  nl=nl+1; write (lines(nl), imt) '%n_top10                    = ', s%global%n_top10
  nl=nl+1; write (lines(nl), imt) '%n_curve_pts                = ', s%global%n_curve_pts
  nl=nl+1; write (lines(nl), imt) '%n_opti_loops               = ', s%global%n_opti_loops
  nl=nl+1; write (lines(nl), imt) '%n_opti_cycles              = ', s%global%n_opti_cycles
  nl=nl+1; write (lines(nl), lmt) '%opt_with_ref               = ', s%global%opt_with_ref 
  nl=nl+1; write (lines(nl), lmt) '%opt_with_base              = ', s%global%opt_with_base
  nl=nl+1; write (lines(nl), amt) '%optimizer                  = ', s%global%optimizer
  nl=nl+1; write (lines(nl), amt) '%phase_units                = ', &
                                                  frequency_units_name(s%global%phase_units)
  nl=nl+1; write (lines(nl), lmt) '%plot_on                    = ', s%global%plot_on
  nl=nl+1; write (lines(nl), lmt) '%lattice_calc_on            = ', s%global%lattice_calc_on
  nl=nl+1; write (lines(nl), lmt) '%command_file_print_on      = ', s%global%command_file_print_on
  nl=nl+1; write (lines(nl), amt) '%prompt_string              = ', s%global%prompt_string
  nl=nl+1; write (lines(nl), amt) '%print_command              = ', s%global%print_command
  nl=nl+1; write (lines(nl), amt) '%random_engine              = ', s%global%random_engine
  nl=nl+1; write (lines(nl), amt) '%random_gauss_converter     = ', s%global%random_gauss_converter
  nl=nl+1; write (lines(nl), rmt) '%random_sigma_cutoff        = ', s%global%random_sigma_cutoff
  nl=nl+1; write (lines(nl), imt) '%random_seed                = ', s%global%random_seed
  if (s%global%random_seed == 0) then
    call ran_seed_get(ix)
    nl=nl+1; write (lines(nl), imt) ' random_seed (generated)    = ', ix
  endif
  nl=nl+1; write (lines(nl), amt) '%track_type                 = ', s%global%track_type
  nl=nl+1; write (lines(nl), imt) '%u_view                     = ', s%global%u_view
  nl=nl+1; write (lines(nl), lmt) '%var_limits_on              = ', s%global%var_limits_on
  nl=nl+1; write (lines(nl), amt) '%var_out_file               = ', s%global%var_out_file
  nl=nl+1; write (lines(nl), rmt) '%y_axis_plot_dmin           = ', s%global%y_axis_plot_dmin

  nl=nl+1; lines(nl) = ''
  nl=nl+1; lines(nl) = 'Internal Variables:'
  nl=nl+1; write (lines(nl), imt) 'Universe index range:        = ', lbound(s%u, 1), ubound(s%u, 1)
  nl=nl+1; write (lines(nl), lmt) 'common_lattice             = ', tao_com%common_lattice
  nl=nl+1; write (lines(nl), amt) 'tao_com%beam_all_file        = ', tao_com%beam_all_file
  nl=nl+1; write (lines(nl), amt) 'tao_com%beam0_file           = ', tao_com%beam0_file
  nl=nl+1; write (lines(nl), lmt) 'tao_com%combine_consecutive_elements_of_like_name = ', &
                                              tao_com%combine_consecutive_elements_of_like_name
  nl=nl+1; write (lines(nl), amt) 'tao_com%init_lat_file        = ', tao_com%init_lat_file
  nl=nl+1; write (lines(nl), amt) 'tao_com%init_tao_file        = ', tao_com%init_tao_file
  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! graph

case ('graph')

  ! Look for switches

  call string_trim(stuff, stuff2, ix)

  ! Find particular graph

  call tao_find_plots (err, stuff2, 'BOTH', graph = graph, &
                                      print_flag = .false., always_allocate = .true.)
  if (err) return

  if (allocated(graph)) then
    g => graph(1)%g
    if (associated(g%p%r)) then
      nl=nl+1; lines(nl) = 'Region.Graph: ' // trim(g%p%r%name) // '.' // trim(g%name)
    endif
    nl=nl+1; lines(nl) = 'Plot.Graph:   ' // trim(g%p%name) // '.' // trim(g%name)
    nl=nl+1; write (lines(nl), amt) 'type                  = ', g%type
    nl=nl+1; write (lines(nl), amt) 'title                 = ', g%title
    nl=nl+1; write (lines(nl), amt) 'title_suffix          = ', g%title_suffix
    nl=nl+1; write (lines(nl), amt) 'component             = ', g%component
    nl=nl+1; write (lines(nl), '(a, 4f10.2, 2x, a)') &
                                    'margin                = ', g%margin
    nl=nl+1; write (lines(nl), imt) 'box                   = ', g%box
    nl=nl+1; write (lines(nl), imt) 'ix_universe           = ', g%ix_universe
    nl=nl+1; write (lines(nl), lmt) 'valid                 = ', g%valid
    nl=nl+1; write (lines(nl), lmt) 'y2_mirrors_y          = ', g%y2_mirrors_y
    nl=nl+1; write (lines(nl), rmt) 'y%max                 = ', g%y%max
    nl=nl+1; write (lines(nl), rmt) 'y%min                 = ', g%y%min
    nl=nl+1; write (lines(nl), imt) 'y%major_div           = ', g%y%major_div
    nl=nl+1; write (lines(nl), imt) 'y%places              = ', g%y%places
    nl=nl+1; write (lines(nl), lmt) 'y%draw_label          = ', g%y%draw_label
    nl=nl+1; write (lines(nl), lmt) 'y%draw_numbers        = ', g%y%draw_numbers
    nl=nl+1; write (lines(nl), rmt) 'y2%max                = ', g%y2%max
    nl=nl+1; write (lines(nl), rmt) 'y2%min                = ', g%y2%min
    nl=nl+1; write (lines(nl), imt) 'y2%major_div          = ', g%y2%major_div
    nl=nl+1; write (lines(nl), imt) 'y2%places             = ', g%y2%places
    nl=nl+1; write (lines(nl), lmt) 'y2%draw_label         = ', g%y2%draw_label
    nl=nl+1; write (lines(nl), lmt) 'y2%draw_numbers       = ', g%y2%draw_numbers
    nl=nl+1; write (lines(nl), lmt) 'limited               = ', g%limited
    nl=nl+1; write (lines(nl), lmt) 'clip                  = ', g%clip
    nl=nl+1; write (lines(nl), lmt) 'draw_axes             = ', g%draw_axes
    nl=nl+1; write (lines(nl), lmt) 'correct_xy_distortion = ', g%correct_xy_distortion
    nl=nl+1; lines(nl) = 'Curves:'
    do i = 1, size(g%curve)
      nl=nl+1; write (lines(nl), amt) '   ', g%curve(i)%name
    enddo

  else
    call out_io (s_error$, r_name, 'This is not a graph')
    return
  endif

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! hom

case ('hom')

  nl=nl+1; lines(nl) = &
        '       #        Freq         R/Q           Q   m  Polarization_Angle'
  do i = 1, size(lat%ele)
    ele => lat%ele(i)
    if (ele%key /= lcavity$) cycle
    if (ele%control_type == multipass_slave$) cycle
    nl=nl+1; write (lines(nl), '(a, i6)') ele%name, i
    do j = 1, size(ele%wake%lr)
      lr => ele%wake%lr(j)
      angle = '-'
      if (lr%polarized) write (angle, '(f9.4)') lr%angle
      nl=nl+1; write (lines(nl), '(i8, 3es12.4, i4, a)') j, &
                  lr%freq, lr%R_over_Q, lr%Q, lr%m, angle
    enddo
    nl=nl+1; lines(nl) = ' '
  enddo
  nl=nl+1; lines(nl) = '       #        Freq         R/Q           Q   m  Polarization_Angle'

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! lattice

case ('lattice')
  
  column(:)%name = ""
  column(:)%label = ""
  column(1)  = show_lat_column_struct("index",   "i6",       6, "")
  column(2)  = show_lat_column_struct("x",       "x",        2, "")
  column(3)  = show_lat_column_struct("name",    "a",        0, "")
  column(4)  = show_lat_column_struct("key",     "a16",     16, "")
  column(5)  = show_lat_column_struct("s",       "f10.3",   10, "")
  column(6)  = show_lat_column_struct("beta_a",  "f7.2",     7, "Beta| A")
  column(7)  = show_lat_column_struct("phi_a",   "f8.3",     8, "")
  column(8)  = show_lat_column_struct("eta_a",   "f5.1",     5, "")
  column(9)  = show_lat_column_struct("orbit_x", "3p, f8.3", 8, "")
  column(10) = show_lat_column_struct("beta_b",  "f7.2",     7, "")
  column(11) = show_lat_column_struct("phi_b",   "f8.3",     8, "")
  column(12) = show_lat_column_struct("eta_b",   "f5.1",     5, "")
  column(13) = show_lat_column_struct("orbit_y", "3p, f8.3", 8, "")

  at_ends = .true.
  by_s = .false.
  print_header_lines = .true.
  ele_name = ''
  if (allocated (picked_ele)) deallocate (picked_ele)
  allocate (picked_ele(0:lat%n_ele_max))

  ! get command line switches

  call string_trim(stuff, stuff2, ix)
  do
    if (ix <= 1) exit

    if (index('-middle', stuff2(1:ix)) == 1) then
      at_ends = .false.

    elseif (index('-no_header', stuff2(1:ix)) == 1) then
      print_header_lines = .false.

    elseif (index('-custom', stuff2(1:ix)) == 1) then
      call string_trim(stuff2(ix+1:), stuff2, ix)
      file_name = stuff2(1:ix)
      iu = lunget()
      open (iu, file = file_name, status = 'old', iostat = ios)
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // file_name)
        return
      endif
      column(:)%name = ""
      column(:)%label = ""
      read (iu, nml = custom_show_list, iostat = ios)
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'CANNOT READ "CUSTOM_SHOW_LIST" NAMELIST IN FILE: ' // file_name)
        return
      endif

    elseif (index('-elements', stuff2(1:ix)) == 1) then
      call string_trim(stuff2(ix+1:), stuff2, ix)
      ele_name = stuff2(1:ix)

    elseif (index('-s', stuff2(1:ix)) == 1) then
      by_s = .true.

    else
      exit
    endif

    call string_trim(stuff2(ix+1:), stuff2, ix)
  enddo
  
  ! Find elements to use

  call re_allocate (picked_ele, 0, lat%n_ele_max)

  if (ele_name /= '') then
    call tao_ele_locations_given_name (u, ele_name, ix_eles, err, .true.)
    if (err) return
    picked_ele = .false.
    picked_ele(ix_eles) = .true.

  elseif (ix == 0 .or. stuff2(1:ix) == 'all') then
    picked_ele = .true.

  elseif (by_s) then
    ix = index(stuff2, ':')
    if (ix == 0) then
      call out_io (s_error$, r_name, 'NO ":" FOUND FOR RANGE SELECTION')
      return
    endif
    read (stuff2(1:ix-1), *, iostat = ios1) s1
    read (stuff2(ix+1:), *, iostat = ios2) s2
    if (ios1 /= 0 .or. ios2 /= 0) then
      call out_io (s_error$, r_name, 'ERROR READING RANGE SELECTION: ' // stuff2)
      return
    endif

    picked_ele = .false.
    do ie = 1, lat%n_ele_track
      if (at_ends) then;
        s_ele = lat%ele(ie)%s
      else
        s_ele = (lat%ele(ie-1)%s + lat%ele(ie)%s) / 2
      endif
      if (s_ele >= s1 .and. s_ele <= s2) picked_ele(ie) = .true.
    enddo

  else
    call location_decode (stuff2, picked_ele, 0, num_locations)
    if (num_locations .eq. -1) then
      call out_io (s_error$, r_name, "SYNTAX ERROR IN RANGE LIST:" // stuff2)
      deallocate(picked_ele)
      return
    endif
  endif

  if (at_ends) then
    write (line1, '(6x, a)') 'Model values at End of Element:'
  else
    write (line1, '(6x, a)') 'Model values at Center of Element:'
  endif

  ix1 = 1
  line2 = ""
  line3 = ""
  do i = 1, size(column)
    if (column(i)%name == "") cycle
    if (column(i)%field_width == 0) then
      do ie = 0, lat%n_ele_track
        if (.not. picked_ele(ie)) cycle
        column(i)%field_width = &
                  max(column(i)%field_width, len_trim(lat%ele(ie)%name)+1)
      enddo
    endif
    ix2 = ix1 + column(i)%field_width

    if (column(i)%label == '') then
      name = column(i)%name
      select case (name)
      case ("name")
        line2(ix1:) = "Name" 
      case ("key")
        line2(ix1:) = "Key" 
      case ("index")
        line2(ix2-5:) = "Index"
      case ("x")
      case default
        if (name(1:5) == "beta_" .or. name(1:6) == "alpha_" .or. name(1:4) == "phi_" .or. &
            name(1:4) == "eta_" .or. name(1:5) == "etap_" .or. name(1:6) == "orbit_") then
          ix = index(name, '_')
          call upcase_string(name(1:1))
          call upcase_string(name(ix+1:ix+1))
          line2(ix2-ix+1:) = name(1:ix-1)
          line3(ix2-ix+2:) = name(ix+1:)
        else
          ix = min (len_trim(name), ix2-ix1+2)
          line2(ix2-ix-1:) = name(1:ix)
        endif
      end select
    else
      name = column(i)%label
      ix = index(name, '|')
      if (ix == 0) then
        j = len_trim(name)
        line2(ix2-j:) = name(1:j)
      else
        j = max(ix-1, len_trim(name(ix+1:)))
        line2(ix2-j:) = name(1:ix-1)
        line3(ix2-j:) = trim(name(ix+1:))
      endif
    endif

    column(i)%format = "(" // trim(column(i)%format) // ")"
    ix1 = ix2
  enddo

  if (print_header_lines) then
    nl=nl+1; lines(nl) = line1
    nl=nl+1; lines(nl) = line2
    if (line3 /= '') then
      nl=nl+1; lines(nl) = line3
    endif
  endif

  do ie = 0, lat%n_ele_track
    if (.not. picked_ele(ie)) cycle
    if (size(lines) < nl+100) call re_allocate (lines, len(lines(1)), nl+200, .false.)
    ele => lat%ele(ie)
    if (ie == 0 .or. at_ends) then
      ele3 = ele
      orb = u%model%orb(ie)
    else
      call twiss_and_track_partial (lat%ele(ie-1), ele, &
                lat%param, ele%value(l$)/2, ele3, u%model%orb(ie-1), orb)
      ele3%s = ele%s-ele%value(l$)/2
    endif

    line = ""
    ix = 1
    do i = 1, size(column)
      if (column(i)%name == "") cycle
      select case (column(i)%name)
      case ("index")
        write (line(ix:), column(i)%format, iostat = ios) ie
      case ("name")
        write (line(ix:), column(i)%format, iostat = ios) ele%name
      case ("key")
        write (line(ix:), column(i)%format, iostat = ios) key_name(ele%key)
      case ("x")
        ios = 0
      case default
        show_common%ele => ele3
        show_common%orbit => orb
        show_common%ix_ele = ie
        show_common%u => u
        call tao_evaluate_expression (column(i)%name, 1, value, good, &
                                             .false., err, tao_ele_value_routine)
        if (err) then
          j = ix + column(i)%field_width - 5
          line(j:) = '-----'
        else
          write (line(ix:), column(i)%format, iostat = ios) value(1)
        endif
      end select

      if (ios /= 0) then
        call out_io (s_error$, r_name, &
            'TOO MANY CHARACTERS ON A LINE OR BAD FORMAT: ' // column(i)%format, &
            'FOR DISPLAYING: ' // column(i)%name)
        return
      endif

      ix  = ix + column(i)%field_width

    enddo

    nl=nl+1; lines(nl) = line

  enddo

  if (print_header_lines) then
    nl=nl+1; lines(nl) = line2
    if (line3 /= '') then
      nl=nl+1; lines(nl) = line3
    endif
    nl=nl+1; lines(nl) = line1
  endif

  if (print_header_lines) then
    first_time = .true.  
    do ie = lat%n_ele_track+1, lat%n_ele_max
      if (.not. picked_ele(ie)) cycle
      if (size(lines) < nl+100) call re_allocate (lines, len(lines(1)), nl+200, .false.)
      ele => lat%ele(ie)
      if (first_time) then
        nl=nl+1; lines(nl) = ' '
        nl=nl+1; lines(nl) = 'Lord Elements:'
        first_time = .false.
      endif
      nl=nl+1
      write (lines(nl), '(i6, 1x, a24, 1x, a16, f10.3, 2(f7.2, f8.3, f5.1, f8.3))') &
            ie, ele%name, key_name(ele%key)
    enddo
  endif

  call out_io (s_blank$, r_name, lines(1:nl))

  deallocate(picked_ele)

!----------------------------------------------------------------------
! optimizer

case ('optimizer')

  do i = lbound(s%u, 1), ubound(s%u, 1)
    u => s%u(i)
    call out_io (s_blank$, r_name, ' ', 'Data Used:')
    write (lines(1), '(a, i4)') 'Universe: ', i
    if (size(s%u) > 1) call out_io (s_blank$, r_name, lines(1))
    do j = 1, size(u%d2_data)
      if (u%d2_data(j)%name == ' ') cycle
      call tao_data_show_use (u%d2_data(j))
    enddo
  enddo

  call out_io (s_blank$, r_name, ' ', 'Variables Used:')
  do j = 1, size(s%v1_var)
    if (s%v1_var(j)%name == ' ') cycle
    call tao_var_show_use (s%v1_var(j))
  enddo

  nl=nl+1; lines(nl) = ' '
  nl=nl+1; write (lines(nl), amt) 'optimizer:        ', s%global%optimizer
  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! optimized_vars

case ('opt_vars')

  call tao_var_write (' ')

!----------------------------------------------------------------------
! particle

case ('orbit')

  read (word(1), *, iostat = ios) ix
  nl=nl+1; write (lines(nl), imt) '  Orbit at Element:', ix
  do i = 1, 6
    nl=nl+1; write (lines(nl), rmt) '     ', u%model%orb(ix)%vec(i)
  enddo
  call out_io (s_blank$, r_name, lines(1:nl))


!----------------------------------------------------------------------
! particle

case ('particle')

  nb = s%global%bunch_to_plot

  if (index('lost', word(1)) == 1) then
    bunch => u%ele(lat%n_ele_track)%beam%bunch(nb)
    nl=nl+1; lines(nl) = 'Particles lost at:'
    nl=nl+1; lines(nl) = '    Ix Ix_Ele  Ele_Name '
    do i = 1, size(bunch%particle)
      if (bunch%particle(i)%ix_lost == not_lost$) cycle
      if (nl == size(lines)) call re_allocate (lines, len(lines(1)), nl+100, .false.)
      ie = bunch%particle(i)%ix_lost
      nl=nl+1; write (lines(nl), '(i6, i7, 2x, a)') i, ie, lat%ele(ie)%name
    enddo
    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif

  bunch => u%ele(0)%beam%bunch(nb)
  read (word(1), *, iostat = ios) ix

  nl=nl+1; write (lines(nl), imt) '  Starting Coords for Particle: ', ix
  do i = 1, 6
    nl=nl+1; write (lines(nl), rmt) '     ', bunch%particle(ix)%r%vec(i)
  enddo
  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! plot

case ('plot')

  ! Look for switches

  show_shape = .false.
  call string_trim(stuff, stuff2, ix)

  do
    if (ix == 0) exit
    if (stuff2(1:1) /= '-') exit
    call match_word (stuff2(:ix), (/ '-shapes' /), n, .true., name)
    if (n < 1) then
      call out_io (s_error$, r_name, 'UNKNOWN SWITCH: ' // stuff2(:ix))
      return
    endif
    if (name == '-shapes') show_shape = .true.
    call string_trim(stuff2(ix+1:), stuff2, ix)
  enddo

  ! Element shapes

  if (show_shape) then

    nl=nl+1; lines(nl) = ' '
    nl=nl+1; lines(nl) = 'Floor_plan Element Shapes:'
    nl=nl+1; lines(nl) = &
          'Ele_Name                        Shape         Color        dy_pix   Draw_Name?'
    nl=nl+1; lines(nl) = &
          '----------------------------    --------      -----        -------  ---------'

    do i = 1, size(tao_com%ele_shape_floor_plan)
      shape => tao_com%ele_shape_floor_plan(i)
      if (shape%ele_name == '') cycle
      nl=nl+1; write (lines(nl), '(3a, f10.4, l3)') &
                shape%ele_name(1:32), shape%shape(1:14), shape%color(1:10), &
                shape%dy_pix, shape%draw_name
    enddo

    nl=nl+1; lines(nl) = ' '
    nl=nl+1; lines(nl) = 'Lat_layout Element Shapes:'
    nl=nl+1; lines(nl) = &
          'Ele_Name                        Shape         Color        dy_pix   Draw_Name?'
    nl=nl+1; lines(nl) = &
          '----------------------------    --------      -----        -------  ---------'

    do i = 1, size(tao_com%ele_shape_lat_layout)
      shape => tao_com%ele_shape_lat_layout(i)
      if (shape%ele_name == '') cycle
      nl=nl+1; write (lines(nl), '(3a, f10.4, l3)') &
                shape%ele_name(1:32), shape%shape(1:14), shape%color(1:10), &
                shape%dy_pix, shape%draw_name
    enddo

    call out_io (s_blank$, r_name, lines(1:nl))
    return 
   
  endif

  ! stuff2 is blank => print overall info

  if (stuff2 == ' ') then

    nl=nl+1; write (lines(nl), f3mt) 'plot_page%text_height            = ', s%plot_page%text_height 
    nl=nl+1; write (lines(nl), f3mt) 'plot_page%main_title_text_scale  = ', s%plot_page%main_title_text_scale 
    nl=nl+1; write (lines(nl), f3mt) 'plot_page%graph_title_text_scale = ', s%plot_page%graph_title_text_scale 
    nl=nl+1; write (lines(nl), f3mt) 'plot_page%axis_number_text_scale = ', s%plot_page%axis_number_text_scale 
    nl=nl+1; write (lines(nl), f3mt) 'plot_page%axis_label_text_scale  = ', s%plot_page%axis_label_text_scale 
    nl=nl+1; write (lines(nl), f3mt) 'plot_page%key_table_text_scale   = ', s%plot_page%key_table_text_scale 
    nl=nl+1; write (lines(nl), f3mt) 'plot_page%shape_height_max       = ', s%plot_page%shape_height_max  

    nl=nl+1; lines(nl) = ' '
    nl=nl+1; lines(nl) = 'Templates:        Plot.Graph'
    nl=nl+1; lines(nl) = '             --------- ----------'
    do i = 1, size(s%template_plot)
      p => s%template_plot(i)
      if (p%name == ' ') cycle
      ix = 21 - len_trim(p%name)
      name = ' '
      name(ix:) = trim(p%name)
      if (allocated(p%graph)) then
        do j = 1, size(p%graph)
          nl=nl+1; write (lines(nl), '(2x, 3a)') name(1:20), '.', p%graph(j)%name
          name = ' '
        enddo
      else
        nl=nl+1; write (lines(nl), '(2x, 2a)') name(1:20), '.'
      endif
      nl=nl+1; lines(nl) = ' '
    enddo

    nl=nl+1; lines(nl) = ' '
    nl=nl+1; lines(nl) = '[Visible]     Plot Region         <-->  Template' 
    nl=nl+1; lines(nl) = '---------     -----------               ------------'
    do i = 1, size(s%plot_region)
      region => s%plot_region(i)
      nl=nl+1; write (lines(nl), '(3x l1, 10x, a20, 2a)') region%visible, &
                                    region%name, '<-->  ', region%plot%name
    enddo

    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif

! Find particular plot

  call tao_find_plots (err, stuff2, 'BOTH', plot, print_flag = .false.)
  if (err) return

  if (allocated(plot)) then
    p => plot(1)%p
    if (associated(p%r)) then
      nl=nl+1; lines(nl) = 'Region:  ' // trim(p%r%name)
    endif
    nl=nl+1; lines(nl) = 'Plot:  ' // p%name
    nl=nl+1; write (lines(nl), amt) 'x_axis_type          = ', p%x_axis_type
    nl=nl+1; write (lines(nl), rmt) 'x%max                = ', p%x%max
    nl=nl+1; write (lines(nl), rmt) 'x%min                = ', p%x%min
    nl=nl+1; write (lines(nl), imt) 'x%major_div_nominal  = ', p%x%major_div_nominal
    nl=nl+1; write (lines(nl), imt) 'x%major_div          = ', p%x%major_div
    nl=nl+1; write (lines(nl), imt) 'x%places             = ', p%x%places
    nl=nl+1; write (lines(nl), lmt) 'x%draw_label         = ', p%x%draw_label
    nl=nl+1; write (lines(nl), lmt) 'x%draw_numbers       = ', p%x%draw_numbers
    nl=nl+1; write (lines(nl), lmt) 'independent_graphs   = ', p%independent_graphs
    
    nl=nl+1; lines(nl) = 'Graphs:'
    do i = 1, size(p%graph)
      nl=nl+1; write (lines(nl), amt) '   ', p%graph(i)%name
    enddo

  else
    call out_io (s_error$, r_name, 'This is not a name of a plot')
    return
  endif

  call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------
! top10

case ('top10')

  call string_trim(stuff, stuff2, ix)
  if (ix == 0) then
    call tao_show_constraints (0, 'TOP10')
    call tao_top10_merit_categories_print (0)
  elseif (index('-derivative', trim(stuff2)) == 1) then 
    call tao_top10_derivative_print ()
  else
    call out_io (s_error$, r_name, 'UNKNOWN SWITCH: ' // stuff2)
    return
  endif

!----------------------------------------------------------------------
! tune

case ('tune')



!----------------------------------------------------------------------
! universe
    
case ('universe')

  if (word(1) == ' ') then
    ix_u = s%global%u_view
  else
    read (word(1), *, iostat = ios) ix_u
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'BAD UNIVERSE NUMBER')
      return
    endif
    if (ix_u < lbound(s%u, 1) .or. ix_u > ubound(s%u, 1)) then
      call out_io (s_error$, r_name, 'UNIVERSE NUMBER OUT OF RANGE')
      return
    endif
  endif

  u => s%u(ix_u)

  nl = 0
  nl=nl+1; write (lines(nl), imt) 'Universe: ', ix_u
  nl=nl+1; write (lines(nl), imt) '%n_d2_data_used        = ', u%n_d2_data_used
  nl=nl+1; write (lines(nl), imt) '%n_data_used           = ', u%n_data_used
  nl=nl+1; write (lines(nl), lmt) '%do_synch_rad_int_calc = ', u%do_synch_rad_int_calc
  nl=nl+1; write (lines(nl), lmt) '%do_chrom_calc         = ', u%do_chrom_calc
  nl=nl+1; write (lines(nl), lmt) '%calc_beam_emittance   = ', u%calc_beam_emittance
  nl=nl+1; write (lines(nl), lmt) '%is_on                 = ', u%is_on
  nl=nl+1; write (lines(nl), amt) '%beam0_file            = ', trim(u%beam0_file)
  nl=nl+1; write (lines(nl), amt) '%beam_all_file         = ', trim(u%beam_all_file)
  nl=nl+1; write (lines(nl), amt) '%save_beam_at:'
  do i = lbound(u%save_beam_at, 1), ubound(u%save_beam_at, 1)
    nl=nl+1; write (lines(nl), '(a, i0, 2a)') '           (', i, ') = ', u%save_beam_at(i)
  enddo
  nl=nl+1; lines(nl) = ''
  nl=nl+1; write(lines(nl), amt) 'Lattice name:    ', lat%lattice
  nl=nl+1; write(lines(nl), amt) 'Input_file_name: ', lat%input_file_name
  nl=nl+1; lines(nl) =           'Lattice Type:    ' // lattice_type(u%model%lat%param%lattice_type)
  nl=nl+1; write (lines(nl), imt) &
                'Elements used in tracking: From 1 through ', lat%n_ele_track
  if (lat%n_ele_max .gt. lat%n_ele_track) then
    nl=nl+1; write (lines(nl), '(a, i0, a, i0)') 'Lord elements:   ', &
                      lat%n_ele_track+1, '  through ', lat%n_ele_max
  else
    nl=nl+1; write (lines(nl), '(a)') "There are NO Lord elements"
  endif

  nl=nl+1; write (lines(nl), '(a, f0.3)') "Lattice length: ", lat%param%total_length
  nl=nl+1; write (lines(nl), '(a)') 'This universe is turned: on_off_logic(u%is_on)'  
  nl=nl+1; write (lines(nl), lmt) 'Aperture limits on?: ', lat%param%aperture_limit_on

  if (.not. lat%param%stable .or. .not. lat%param%stable) then
    nl=nl+1; write (lines(nl), '(a, l)') 'Model lattice stability: ', &
                                                          lat%param%stable
    nl=nl+1; write (lines(nl), '(a, l)') 'Design lattice stability:', &
                                                          u%design%lat%param%stable
    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif
 
  call radiation_integrals (lat, &
                                u%model%orb, u%model%modes, u%ix_rad_int_cache)
  call radiation_integrals (u%design%lat, &
                                u%design%orb, u%design%modes, u%ix_rad_int_cache)
  if (lat%param%lattice_type == circular_lattice$) then
    call chrom_calc (lat, delta_e, &
                        u%model%a%chrom, u%model%b%chrom, exit_on_error = .false.)
    call chrom_calc (u%design%lat, delta_e, &
                        u%design%a%chrom, u%design%b%chrom, exit_on_error = .false.)
  endif

  nl=nl+1; lines(nl) = ''
  nl=nl+1; write (lines(nl), '(17x, a)') '       X          |            Y'
  nl=nl+1; write (lines(nl), '(17x, a)') 'Model     Design  |     Model     Design'
  fmt = '(1x, a10, 1p 2e11.3, 2x, 2e11.3, 2x, a)'
  fmt2 = '(1x, a10, 2f11.3, 2x, 2f11.3, 2x, a)'
  fmt3 = '(1x, a10, 2f11.4, 2x, 2f11.4, 2x, a)'
  f_phi = 1 / twopi
  l_lat = lat%param%total_length
  n = lat%n_ele_track
  if (lat%param%lattice_type == circular_lattice$) then
    nl=nl+1; write (lines(nl), fmt2) 'Q', f_phi*lat%ele(n)%a%phi, &
          f_phi*u%design%lat%ele(n)%a%phi, f_phi*lat%ele(n)%b%phi, &
          f_phi*u%design%lat%ele(n)%b%phi,  '! Tune'
    nl=nl+1; write (lines(nl), fmt2) 'Chrom', u%model%a%chrom, & 
          u%design%a%chrom, u%model%b%chrom, u%design%b%chrom, '! dQ/(dE/E)'
    nl=nl+1; write (lines(nl), fmt2) 'J_damp', u%model%modes%a%j_damp, &
        u%design%modes%a%j_damp, u%model%modes%b%j_damp, &
        u%design%modes%b%j_damp, '! Damping Partition #'
    nl=nl+1; write (lines(nl), fmt) 'Emittance', u%model%modes%a%emittance, &
        u%design%modes%a%emittance, u%model%modes%b%emittance, &
        u%design%modes%b%emittance, '! Meters'
  endif
  nl=nl+1; write (lines(nl), fmt) 'Alpha_damp', u%model%modes%a%alpha_damp, &
        u%design%modes%a%alpha_damp, u%model%modes%b%alpha_damp, &
        u%design%modes%b%alpha_damp, '! Damping per turn'
  nl=nl+1; write (lines(nl), fmt) 'I4', u%model%modes%a%synch_int(4), &
        u%design%modes%a%synch_int(4), u%model%modes%b%synch_int(4), &
        u%design%modes%b%synch_int(4), '! Radiation Integral'
  nl=nl+1; write (lines(nl), fmt) 'I5', u%model%modes%a%synch_int(5), &
        u%design%modes%a%synch_int(5), u%model%modes%b%synch_int(5), &
        u%design%modes%b%synch_int(5), '! Radiation Integral'

  nl=nl+1; lines(nl) = ''
  nl=nl+1; write (lines(nl), '(19x, a)') 'Model     Design'
  fmt = '(1x, a12, 1p2e11.3, 3x, a)'
  nl=nl+1; write (lines(nl), fmt) 'Sig_E/E:', u%model%modes%sigE_E, &
            u%design%modes%sigE_E
  nl=nl+1; write (lines(nl), fmt) 'Energy Loss:', u%model%modes%e_loss, &
            u%design%modes%e_loss, '! Energy_Loss (eV / Turn)'
  nl=nl+1; write (lines(nl), fmt) 'J_damp:', u%model%modes%z%j_damp, &
        u%design%modes%z%j_damp, '! Longitudinal Damping Partition #'
  nl=nl+1; write (lines(nl), fmt) 'Alpha_damp:', u%model%modes%z%alpha_damp, &
        u%design%modes%z%alpha_damp, '! Longitudinal Damping per turn'
  nl=nl+1; write (lines(nl), fmt) 'Alpha_p:', u%model%modes%synch_int(1)/l_lat, &
               u%design%modes%synch_int(1)/l_lat, '! Momentum Compaction'
  nl=nl+1; write (lines(nl), fmt) 'I1:', u%model%modes%synch_int(1), &
               u%design%modes%synch_int(1), '! Radiation Integral'
  nl=nl+1; write (lines(nl), fmt) 'I2:', u%model%modes%synch_int(2), &
               u%design%modes%synch_int(2), '! Radiation Integral'
  nl=nl+1; write (lines(nl), fmt) 'I3:', u%model%modes%synch_int(3), &
               u%design%modes%synch_int(3), '! Radiation Integral'

  call out_io (s_blank$, r_name, lines(1:nl)) 

!----------------------------------------------------------------------
! variable
    
case ('use')  

  nl=nl+1; lines(nl) = 'veto data *@*'
  nl=nl+1; lines(nl) = ''

  do i = lbound(s%u, 1), ubound(s%u, 1)
    do j = 1, s%u(i)%n_d2_data_used
      d2_ptr => s%u(i)%d2_data(j)
      call re_allocate (lines, n_char, nl+size(d2_ptr%d1)+10, .false.)
      do k = lbound(d2_ptr%d1, 1), ubound(d2_ptr%d1, 1)
        d1_ptr => d2_ptr%d1(k)
        call location_encode(line, d1_ptr%d%useit_opt, &
                                  d1_ptr%d%exists, lbound(d1_ptr%d, 1))
        if (line == '') cycle
        nl=nl+1; write (lines(nl), '(5a)') 'use data ', &
                      trim(tao_d2_d1_name(d1_ptr)), '[', trim(line), ']'
      enddo
    enddo
  enddo
  nl=nl+1; lines(nl) = ''

  call re_allocate (lines, n_char, nl+size(s%v1_var)+10, .false.)
  do i = 1, size(s%v1_var)
    v1_ptr => s%v1_var(i)
    if (v1_ptr%name == ' ') cycle
    call re_allocate (lines, len(lines(1)), nl+200, .false.)
    call location_encode (line, v1_ptr%v%useit_opt, v1_ptr%v%exists, lbound(v1_ptr%v, 1))
    nl=nl+1; write (lines(nl), '(5a)') 'use var ', trim(v1_ptr%name), '[', trim(line), ']'
  enddo

  call out_io (s_blank$, r_name, lines(1:nl))


!----------------------------------------------------------------------
! variable
    
case ('variable')

  if (.not. allocated (s%v1_var)) then
    call out_io (s_error$, r_name, 'NO VARIABLES HAVE BEEN DEFINED IN THE INPUT FILES!')
    return 
  endif

! If 'n@' is present then write out stuff for universe n

  ix = index(word(1), '@')
  if (ix /= 0) then
    if (ix == 1) then
      ix_u = s%global%u_view
    else
      read (word(1)(:ix-1), *, iostat = ios) ix_u
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'BAD UNIVERSE NUMBER')
        return
      endif
      if (ix_u == 0) ix_u = s%global%u_view
      if (ix_u < 1 .or. ix_u > ubound(s%u, 1)) then
        call out_io (s_error$, r_name, 'UNIVERSE NUMBER OUT OF RANGE')
        return
      endif
    endif
    write (lines(1), '(a, i4)') 'Variables controlling universe:', ix_u
    write (lines(2), '(5x, a)') '                    '
    write (lines(3), '(5x, a)') 'Name                '
    nl = 3
    do i = 1, size(s%var)
      if (.not. s%var(i)%exists) cycle
      found = .false.
      do j = 1, size(s%var(i)%this)
        if (s%var(i)%this(j)%ix_uni == ix_u) found = .true.
      enddo
      if (.not. found) cycle
      nl=nl+1; write(lines(nl), '(5x, a25, a40)') tao_var1_name(s%var(i)), &
                                                      tao_var_attrib_name(s%var(i))
    enddo
    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif

! If just "show var" then show all namees

  if (word(1) == '*') then
    call tao_var_write (' ')
    return
  endif

  if (word(1) == ' ') then
    nl=nl+1; write (lines(nl), '(7x, a, t50, a)') 'Name', 'Using'
    do i = 1, size(s%v1_var)
      v1_ptr => s%v1_var(i)
      if (v1_ptr%name == ' ') cycle
      call re_allocate (lines, len(lines(1)), nl+200, .false.)
      call location_encode (line, v1_ptr%v%useit_opt, v1_ptr%v%exists, lbound(v1_ptr%v, 1))
      nl=nl+1; write(lines(nl), '(i5, 2x, 2a, i0, a, i0, a, t50, a)') v1_ptr%ix_v1, &
                      trim(v1_ptr%name), '[', lbound(v1_ptr%v, 1), ':', &
                      ubound(v1_ptr%v, 1), ']', trim(line)
      
    enddo
    call out_io (s_blank$, r_name, lines(1:nl))
    return
  endif

! get pointers to the variables

  call string_trim (word(2), word(2), ix)

! are we looking at a range of locations?

  call tao_find_var(err, word(1), v1_array, v_array) 
  if (err) return
  n_size = 0
  if (allocated(v_array)) n_size = size(v_array)

! v_ptr is valid then show the variable info.

  if (n_size == 1) then

    v_ptr => v_array(1)%v

    nl=nl+1; write(lines(nl), amt)  'Ele_name      = ', v_ptr%ele_name
    nl=nl+1; write(lines(nl), amt)  'Attrib_name   = ', v_ptr%attrib_name 
    nl=nl+1; write(lines(nl), imt)  'Ix_attrib     = ', v_ptr%ix_attrib 
    nl=nl+1; write(lines(nl), imt)  'Ix_var        = ', v_ptr%ix_var
    nl=nl+1; write(lines(nl), imt)  'Ix_dvar       = ', v_ptr%ix_dvar           
    nl=nl+1; write(lines(nl), imt)  'Ix_v1         = ', v_ptr%ix_v1
    nl=nl+1; write(lines(nl), rmt)  'Model_value   = ', v_ptr%model_value
    nl=nl+1; write(lines(nl), rmt)  'Base_value    = ', v_ptr%base_value

    if (.not. allocated (v_ptr%this)) then
      nl=nl+1; write(lines(nl), imt)  'this(:) -- Not associated!'
    else
      do i = 1, size(v_ptr%this)
        nl=nl+1; write(lines(nl), iimt)  '%this(', i, ')%ix_uni:        ', &
                                                            v_ptr%this(i)%ix_uni
        nl=nl+1; write(lines(nl), iimt)  '%this(', i, ')%ix_ele:        ', v_ptr%this(i)%ix_ele
        if (associated (v_ptr%this(i)%model_value)) then
          nl=nl+1; write(lines(nl), irmt)  '%this(', i, ')%Model_value: ', &
                                                            v_ptr%this(i)%model_value
        else
          nl=nl+1; write(lines(nl), irmt)  '%this(', i, ')%Model_value: <not associated>'
        endif
        if (associated (v_ptr%this(i)%base_value)) then
          nl=nl+1; write(lines(nl), irmt)  '%this(', i, ')%Base_value:  ', &
                                                            v_ptr%this(i)%base_value
        else
          nl=nl+1; write(lines(nl), irmt)  '%this(', i, ')%Base_value:  <not associated>'
        endif
      enddo
    endif

    if (associated (v_ptr%common%model_value)) then
      nl=nl+1; write(lines(nl), imt)  '%common%ix_uni:      ', v_ptr%common%ix_uni
      nl=nl+1; write(lines(nl), imt)  '%common%ix_ele:      ', v_ptr%common%ix_ele
      nl=nl+1; write(lines(nl), rmt)  '%common%Model_value: ', v_ptr%common%model_value
      nl=nl+1; write(lines(nl), rmt)  '%common%Base_value:  ', v_ptr%common%base_value
    endif

    nl=nl+1; write(lines(nl), rmt)  '%Design_value     = ', v_ptr%design_value
    nl=nl+1; write(lines(nl), rmt)  '%Old_value        = ', v_ptr%old_value
    nl=nl+1; write(lines(nl), rmt)  '%Meas_value       = ', v_ptr%meas_value
    nl=nl+1; write(lines(nl), rmt)  '%Ref_value        = ', v_ptr%ref_value
    nl=nl+1; write(lines(nl), rmt)  '%Correction_value = ', v_ptr%correction_value
    nl=nl+1; write(lines(nl), rmt)  '%High_lim         = ', v_ptr%high_lim
    nl=nl+1; write(lines(nl), rmt)  '%Low_lim          = ', v_ptr%low_lim
    nl=nl+1; write(lines(nl), rmt)  '%Step             = ', v_ptr%step
    nl=nl+1; write(lines(nl), rmt)  '%Weight           = ', v_ptr%weight
    nl=nl+1; write(lines(nl), rmt)  '%delta_merit      = ', v_ptr%delta_merit
    nl=nl+1; write(lines(nl), amt)  '%Merit_type       = ', v_ptr%merit_type
    nl=nl+1; write(lines(nl), rmt)  '%Merit            = ', v_ptr%merit
    nl=nl+1; write(lines(nl), rmt)  '%dMerit_dVar      = ', v_ptr%dMerit_dVar
    nl=nl+1; write(lines(nl), lmt)  '%Key_bound        = ', v_ptr%key_bound
    if (v_ptr%key_bound) then
      nl=nl+1; write(lines(nl), rmt)  '%Key_val0         = ', v_ptr%key_val0
      nl=nl+1; write(lines(nl), rmt)  '%Key_delta        = ', v_ptr%key_delta
    endif
    nl=nl+1; write(lines(nl), lmt)  '%Exists           = ', v_ptr%exists
    nl=nl+1; write(lines(nl), lmt)  '%Good_var         = ', v_ptr%good_var
    nl=nl+1; write(lines(nl), lmt)  '%Good_user        = ', v_ptr%good_user
    nl=nl+1; write(lines(nl), lmt)  '%Good_opt         = ', v_ptr%good_opt
    nl=nl+1; write(lines(nl), lmt)  '%Useit_opt        = ', v_ptr%useit_opt
    nl=nl+1; write(lines(nl), lmt)  '%Useit_plot       = ', v_ptr%useit_plot

    call tao_hook_show_variable (v_ptr, lines, nl)

! check if there is a variable number
! if no variable number requested, show a range

  elseif (size(v1_array) == 1) then

    nc = 0
    do i = 1, size(v_array)
      v_ptr => v_array(i)%v
      if (.not. v_ptr%exists) cycle
      nc = max(nc, len_trim(tao_var_attrib_name(v_ptr)))
    enddo

    write(lines(1), '(2a)') 'Variable name:  ', v1_array(1)%v1%name
    lines(2) = ' '
    line1 = '       Name'
    line1(nc+17:) = 'Meas         Model        Design  Useit_opt'
    lines(3) = line1
    nl = 3
    ! if a range is specified, show the variable range   
    do i = 1, size(v_array)
      v_ptr => v_array(i)%v
      if (.not. v_ptr%exists) cycle
      call re_allocate (lines, len(lines(1)), nl+200, .false.)
      nl=nl+1
      write(lines(nl), '(i6, 2x, a)') v_ptr%ix_v1, tao_var_attrib_name(v_ptr)
      write(lines(nl)(nc+9:), '(3es14.4, 7x, l)') v_ptr%meas_value, &
                 v_ptr%model_value, v_ptr%design_value, v_ptr%useit_opt
    enddo
    nl=nl+1; lines(nl) = line1

  else
    lines(1) = '???'
    nl = 1
  endif

! print out results

  call out_io (s_blank$, r_name, lines(1:nl))


!----------------------------------------------------------------------

case default

  call out_io (s_error$, r_name, "INTERNAL ERROR, SHOULDN'T BE HERE!")
  return

end select

!----------------------------------------------------------------------
!----------------------------------------------------------------------
contains

subroutine show_ele_data (u, i_ele, lines, nl)

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct), pointer :: datum
character(*) :: lines(:)
character(100) l1
integer i_ele, nl, i

logical found_one

!

nl=nl+1; write (lines(nl), '(a)') "  "
write (l1, '(a, 20x, a)') "Data Name", &
          "Data Type             |  Model Value  |  Design Value |  Base Value"
nl=nl+1; lines(nl) = l1

found_one = .false.
do i = 1, size(u%data)
  if (u%data(i)%ix_ele .eq. i_ele) then
    found_one = .true.
    datum => u%data(i)
    nl=nl+1; write (lines(nl), "(a, t30, a20, 3(1x, es15.5))") &
                trim(tao_datum_name(datum)),  datum%data_type, datum%model_value, &
                datum%design_value, datum%base_value 
  endif
enddo

if (found_one) then
  nl=nl+1; lines(nl) = l1
else
  write (lines(nl), '(a)') "No data associated with this element."
endif

end subroutine show_ele_data

end subroutine tao_show_this

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine tao_ele_value_routine (str, value, good, err_flag)

implicit none

real(rp), allocatable :: value(:)
logical, allocatable :: good(:)
type (tao_universe_struct), pointer :: u
integer ios, i, n, ie

character(*) str
character(40) attribute
character(24) :: r_name = 'tao_ele_value_routine'

real(rp), pointer :: real_ptr

logical err_flag

!

call re_allocate (value, 1, .true.)
value = 0  ! Default
err_flag = .true.

ie = show_common%ix_ele
u => show_common%u

attribute = str
select case (attribute)
case ("orbit_x")
  value(1) = show_common%orbit%vec(1)
case ("orbit_px")
  value(1) = show_common%orbit%vec(2)
case ("orbit_y")
  value(1) = show_common%orbit%vec(3)
case ("orbit_py")
  value(1) = show_common%orbit%vec(4)
case ("orbit_z")
  value(1) = show_common%orbit%vec(5)
case ("orbit_pz")
  value(1) = show_common%orbit%vec(6)
case ("sigma_x", "sigma_y", "sigma_z", "sigma_px", "sigma_py", "sigma_pz")
  if (attribute == "sigma_x")  value(1) = sqrt(u%model%bunch_params(ie)%sigma(s11$))
  if (attribute == "sigma_px") value(1) = sqrt(u%model%bunch_params(ie)%sigma(s22$))
  if (attribute == "sigma_y")  value(1) = sqrt(u%model%bunch_params(ie)%sigma(s33$))
  if (attribute == "sigma_py") value(1) = sqrt(u%model%bunch_params(ie)%sigma(s44$))
  if (attribute == "sigma_z")  value(1) = sqrt(u%model%bunch_params(ie)%sigma(s55$))
  if (attribute == "sigma_pz") value(1) = sqrt(u%model%bunch_params(ie)%sigma(s66$))
case ("i5a_e6")
  if (.not. allocated (u%model%rad_int%lin_i5a_e6)) return
  value(1) = sum(u%model%rad_int%lin_i5a_e6(1:ie))
case ("i5b_e6")
  if (.not. allocated (u%model%rad_int%lin_i5b_e6)) return
  value(1) = sum(u%model%rad_int%lin_i5b_e6(1:ie))

! Must be an element attribute
case default
  call upcase_string(attribute)
  call pointer_to_attribute (show_common%ele, attribute, .true., real_ptr, err_flag, .false.)
  if (.not. err_flag) value(1) = real_ptr
end select

err_flag = .false.

end subroutine

end module
