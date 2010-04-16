module tao_show_mod

use tao_mod
use tao_data_and_eval_mod
use tao_top10_mod
use tao_command_mod, only: tao_next_switch
use tao_lattice_calc_mod

type show_common_struct
  type (ele_struct), pointer :: ele 
  type (coord_struct), pointer :: orbit 
  type (bunch_params_struct), pointer :: bunch_params
  type (tao_universe_struct), pointer :: u
  integer ix_ele
end type

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


integer iu, ix, n, nl
integer :: n_write_file = 0            ! used for indexing 'show write' files

character(*) what, stuff
character(100) file_name, result_id
character(len(what)) what2
character(len(stuff)) stuff2
character(n_char), allocatable, save :: lines(:)
character(16) :: r_name = 'tao_show_cmd'
character(20) switch

logical opened, err

! Init

what2 = what
stuff2 = stuff
opened = .false.

! See if the results need to be written to a file.

call tao_next_switch (what2, ['-append', '-write '], switch, err, ix)
if (err) return

if (switch /= '') then
  call string_trim(stuff, stuff2, ix)
  file_name = stuff2(:ix)
  call string_trim(stuff2(ix+1:), stuff2, ix)

  ix = index(file_name, '*')
  if (ix /= 0) then
    n_write_file = n_write_file + 1
    write (file_name, '(a, i3.3, a)') file_name(1:ix-1), n_write_file, trim(file_name(ix+1:))
  endif

  iu = lunget()
  if (switch == '-append') then
    open (iu, file = file_name, position = 'APPEND', status = 'UNKNOWN', recl = n_char)
  else
    open (iu, file = file_name, status = 'REPLACE', recl = n_char)
  endif

  call output_direct (iu)  ! tell out_io to write to a file

  call string_trim (stuff2, stuff2, ix)
  what2 = stuff2(1:ix)
  stuff2 = stuff2(ix+1:)
  opened = .true.
endif

! Get restults.
! Result_id is for tao_show_this to show exactly what it did.
! This info can be helpful in tao_hook_show_cmd.

call tao_show_this (what2, stuff2, result_id, lines, nl)  
call tao_hook_show_cmd (what2, stuff2, result_id, lines, nl)

if (nl > 0) then
  if (result_id == 'ERROR') then
    call out_io (s_error$, r_name, lines(1:nl))
  else
    call out_io (s_blank$, r_name, lines(1:nl))
  endif
endif

! Finish

if (opened) then
  call output_direct (0)  ! reset to not write to a file
  close (iu)
  call out_io (s_blank$, r_name, 'Written to file: ' // file_name)

  call tao_show_this (what, stuff, result_id, lines, nl)
endif

end subroutine

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine tao_show_this (what, stuff, result_id, lines, nl)

use random_mod
use csr_mod, only: csr_param
use location_encode_mod
use transfer_map_mod
use opti_de_mod

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
type (tao_curve_struct), pointer :: c1
type (tao_plot_region_struct), pointer :: region
type (tao_d1_data_array_struct), allocatable, save :: d1_array(:)
type (tao_data_array_struct), allocatable, save :: d_array(:)
type (tao_ele_shape_struct), pointer :: shape
type (beam_struct), pointer :: beam
type (lat_struct), pointer :: lat
type (bunch_struct), pointer :: bunch
type (lr_wake_struct), pointer :: lr
type (ele_struct), pointer :: ele
type (coord_struct), target :: orb
type (ele_struct), target :: ele3
type (bunch_params_struct) bunch_params
type (bunch_params_struct), pointer :: bunch_p
type (taylor_struct) taylor(6)
type (ele_pointer_struct), allocatable, save :: eles(:)
type (branch_struct), pointer :: branch
type (tao_universe_branch_struct), pointer :: uni_branch

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
character(24) :: var_name, blank_str = ''
character(24)  :: plane, imt, lmt, amt, iamt, ramt, f3mt, rmt, irmt, iimt
character(80) :: word1, fmt, fmt2, fmt3
character(20) :: r_name = "tao_show_cmd"
character(24) show_name, show2_name
character(200), pointer :: ptr_lines(:)
character(100), allocatable :: alloc_lines(:)
character(100) file_name, name
character(120) header, str
character(40) ele_name, sub_name, ele1, ele2, switch
character(60) nam
character(3) undef_str
character(40) replacement_for_blank

character(16) :: show_what, show_names(25) = [ &
   'data        ', 'variable    ', 'global      ', 'alias       ', 'top10       ', &
   'optimizer   ', 'element     ', 'lattice     ', 'constraints ', 'plot        ', &
   'beam        ', 'tune        ', 'graph       ', 'curve       ', 'particle    ', &
   'hom         ', 'key_bindings', 'universe    ', 'orbit       ', 'derivative  ', &
   'branches    ', 'use         ', 'taylor_map  ', 'value       ', 'wave        ' ]

character(*), allocatable :: lines(:)
character(*) result_id
character(n_char) line, line1, line2, line3
character(n_char) stuff2
character(9) angle

integer :: data_number, ix_plane, ix_class, n_live, n_order, i1, i2, ix_branch
integer nl, loc, ixl, iu, nc, n_size, ix_u, ios, ie, nb, id, iv, jd, jv
integer ix, ix1, ix2, ix_s2, i, j, k, n, show_index, ju, ios1, ios2, i_uni
integer num_locations, ix_ele, n_name, n_start, n_ele, n_ref, n_tot, ix_p, ix_word

logical bmad_format, good_opt_only, show_lords, show_custom
logical err, found, at_ends, first_time, by_s, print_header_lines, all_lat, limited
logical show_sym, show_line, show_shape, print_data, ok, print_tail_lines
logical show_all, name_found, print_taylor, print_wig_terms, print_all
logical, allocatable, save :: picked_uni(:)
logical, allocatable, save :: picked_ele(:)
logical, allocatable, save :: good(:)
logical print_global, print_optimization, print_bmad_com, print_csr_param 

namelist / custom_show_list / column

!

call re_allocate (lines, 200)

err = .false.

lines = " "
nl = 0

rmt  = '(a, 9es16.8)'
f3mt  = '(a, 9f0.3)'
irmt = '(a, i0, a, es16.8)'
imt  = '(a, 9i8)'
iimt = '(a, i0, a, i8)'
lmt  = '(a, 9(l1, 2x))'
amt  = '(9a)'
iamt = '(a, i0, 2x, 9a)'
ramt = '(a, f0.3, 2x, 9a)'

u => tao_pointer_to_universe(-1)
lat => u%model%lat

if (s%global%phase_units == radians$) f_phi = 1
if (s%global%phase_units == degrees$) f_phi = 180 / pi
if (s%global%phase_units == cycles$)  f_phi = 1 / twopi

! find what to show

result_id = 'ERROR'

if (what == ' ') then
  nl=1; lines(1) = 'SHOW WHAT?' 
  return
endif

call match_word (what, show_names, ix, matched_name = show_what)
if (ix == 0) then
  nl=1; lines(1) = 'SHOW WHAT? WORD NOT RECOGNIZED: ' // what
  return
endif

if (ix < 0) then
  nl=1; lines(1) = 'SHOW WHAT? AMBIGUOUS: ' // what
  return
endif

call string_trim (stuff, stuff2, ix_word)
word1 = stuff2(:ix_word)

select case (show_what)

!----------------------------------------------------------------------
! alias

case ('alias')

  call re_allocate (lines, tao_com%n_alias+10, .false.)
  lines(1) = 'Aliases:'
  nl = 1
  do i = 1, tao_com%n_alias
    nl=nl+1; lines(nl) = trim(tao_com%alias(i)%name) // ' = "' // &
                                    trim(tao_com%alias(i)%string) // '"'
  enddo
  
  result_id = show_what

!----------------------------------------------------------------------
! beam

case ('beam')

  ! no element index

  if (word1 == '') then

    ix_branch = 0
    uni_branch => u%uni_branch(ix_branch)
    nl=nl+1; write(lines(nl), '(a, i3)') 'Universe: ', u%ix_uni
    nl=nl+1; write(lines(nl), '(a, i3)') 'Branch:   ', ix_branch
    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), amt) 'beam0_file                  = ', u%beam_info%beam0_file
    nl=nl+1; write(lines(nl), amt) 'beam_all_file               = ', u%beam_info%beam_all_file
    beam => uni_branch%ele(0)%beam
    if (allocated(beam%bunch)) then
      nl=nl+1; write(lines(nl), imt) 'n_particle                  = ', size(beam%bunch(1)%particle)
      nl=nl+1; write(lines(nl), imt) 'n_bunch                     = ', size(beam%bunch)
      nl=nl+1; write(lines(nl), rmt) 'bunch_charge                = ', beam%bunch(1)%charge
    endif
    if (u%beam_info%beam_all_file == '' .and. u%beam_info%beam0_file == '') then
      nl=nl+1; write(lines(nl), rmt) 'beam_init%center            = ', u%beam_info%beam_init%center
      nl=nl+1; write(lines(nl), rmt) 'beam_init%a_norm_emitt      = ', u%beam_info%beam_init%a_norm_emitt
      nl=nl+1; write(lines(nl), rmt) 'beam_init%b_norm_emitt      = ', u%beam_info%beam_init%b_norm_emitt
      nl=nl+1; write(lines(nl), rmt) 'beam_init%dPz_dz            = ', u%beam_info%beam_init%dPz_dz
      nl=nl+1; write(lines(nl), rmt) 'beam_init%dt_bunch          = ', u%beam_info%beam_init%dt_bunch
      nl=nl+1; write(lines(nl), rmt) 'beam_init%sig_z             = ', u%beam_info%beam_init%sig_z
      nl=nl+1; write(lines(nl), rmt) 'beam_init%sig_e             = ', u%beam_info%beam_init%sig_e
      nl=nl+1; write(lines(nl), rmt) 'beam_init%center_jitter     = ', u%beam_info%beam_init%center_jitter
      nl=nl+1; write(lines(nl), rmt) 'beam_init%emitt_jitter      = ', u%beam_info%beam_init%emitt_jitter
      nl=nl+1; write(lines(nl), rmt) 'beam_init%sig_z_jitter      = ', u%beam_info%beam_init%sig_z_jitter
      nl=nl+1; write(lines(nl), rmt) 'beam_init%sig_e_jitter      = ', u%beam_info%beam_init%sig_e_jitter
      nl=nl+1; write(lines(nl), rmt) 'beam_init%spin%polarization = ', u%beam_info%beam_init%spin%polarization
      nl=nl+1; write(lines(nl), rmt) 'beam_init%spin%theta        = ', u%beam_info%beam_init%spin%theta
      nl=nl+1; write(lines(nl), rmt) 'beam_init%spin%phi          = ', u%beam_info%beam_init%spin%phi
      nl=nl+1; write(lines(nl), lmt) 'beam_init%renorm_center     = ', u%beam_info%beam_init%renorm_center
      nl=nl+1; write(lines(nl), lmt) 'beam_init%renorm_sigma      = ', u%beam_info%beam_init%renorm_sigma
      nl=nl+1; write(lines(nl), lmt) 'beam_init%init_spin         = ', u%beam_info%beam_init%init_spin
    endif
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
    nl=nl+1; write(lines(nl), amt)  'global%track_type          = ', s%global%track_type
    nl=nl+1; write(lines(nl), lmt)  'global%beam_timer_on       = ', s%global%beam_timer_on
    nl=nl+1; write (lines(nl), imt) 'ix_track_start             = ', uni_branch%ix_track_start 
    nl=nl+1; write (lines(nl), imt) 'ix_track_end               = ', uni_branch%ix_track_end
    nl=nl+1; write(lines(nl), amt)  'u%beam_saved_at:           = ', trim(u%beam_saved_at)

  ! have element index

  else
    call tao_pick_universe (word1, word1, picked_uni, err, ix_u)
    if (err) return
    u => s%u(ix_u)
    call tao_locate_elements (word1, ix_u, eles, err)
    if (err .or. size(eles) == 0) return
    ix_ele = eles(1)%ele%ix_ele
    n = s%global%bunch_to_plot
    bunch_p => u%model%lat_branch(eles(1)%ele%ix_branch)%bunch_params(ix_ele)
    nl=nl+1; lines(nl) = 'Cashed bunch parameters:'
    nl=nl+1; write (lines(nl), rmt) '  Centroid:', bunch_p%centroid%vec
    nl=nl+1; write (lines(nl), rmt) '  RMS:     ', &
                              sqrt(bunch_p%sigma([s11$, s22$, s33$, s44$, s55$, s66$]))
    nl=nl+1; write (lines(nl), rmt) '             norm_emitt           beta'
    nl=nl+1; write (lines(nl), rmt) '  a:       ', bunch_p%a%norm_emit, bunch_p%a%beta
    nl=nl+1; write (lines(nl), rmt) '  b:       ', bunch_p%b%norm_emit, bunch_p%b%beta
    nl=nl+1; write (lines(nl), rmt) '  x:       ', bunch_p%x%norm_emit, bunch_p%x%beta
    nl=nl+1; write (lines(nl), rmt) '  y:       ', bunch_p%y%norm_emit, bunch_p%y%beta
    nl=nl+1; write (lines(nl), rmt) '  z:       ', bunch_p%z%norm_emit, bunch_p%z%beta

    beam => u%uni_branch(eles(1)%ele%ix_branch)%ele(ix_ele)%beam
    if (allocated(beam%bunch)) then
      bunch => beam%bunch(n)
      call calc_bunch_params (bunch, lat%ele(ix_ele), lat%param, bunch_params, err)
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
                         sqrt(bunch_params%sigma([s11$, s22$, s33$, s44$, s55$, s66$]))
      nl=nl+1; write (lines(nl), rmt) '             norm_emitt           beta'
      nl=nl+1; write (lines(nl), rmt) '  a:       ', &
                                bunch_params%a%norm_emit, bunch_params%a%beta
      nl=nl+1; write (lines(nl), rmt) '  b:       ', &
                                bunch_params%b%norm_emit, bunch_params%b%beta
      nl=nl+1; write (lines(nl), rmt) '  x:       ', &
                                bunch_params%x%norm_emit, bunch_params%x%beta
      nl=nl+1; write (lines(nl), rmt) '  y:       ', &
                                bunch_params%y%norm_emit, bunch_params%y%beta
      nl=nl+1; write (lines(nl), rmt) '  z:       ', &
                                bunch_params%z%norm_emit, bunch_params%z%beta
    else
      nl=nl+1; lines(nl) = 'No allocated beam at element.'
    endif
  
  endif

  result_id = show_what

!----------------------------------------------------------------------
! constraints

case ('branches')

  if (ubound(lat%branch, 1) == 0) then
    nl=1; lines(1) = 'NO BRANCHES'
    result_id = show_what
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

  result_id = show_what

!----------------------------------------------------------------------
! constraints

case ('constraints')

  call tao_show_constraints (0, '*')
  call tao_show_constraints (0, 'TOP10')
  result_id = show_what

!----------------------------------------------------------------------
! curve

case ('curve')

  ! Look for switches

  show_sym = .false.
  show_line = .false.

  do
    call tao_next_switch (stuff2, ['-symbol', '-line  ' ], switch, err, ix)
    if (err) return
    if (switch == '') exit
    if (switch == '-symbol') show_sym = .true.
    if (switch == '-line')   show_line = .true.
  enddo

  ! Find particular plot

  call tao_find_plots (err, stuff2, 'BOTH', curve = curve, always_allocate = .true.)
  if (err) return

  ! print info on particular plot, graph, or curve

  if (allocated(curve)) then
    c1 => curve(1)%c
    nl=nl+1; lines(nl) = 'Region.Graph.Curve: ' // trim(tao_curve_name(c1, .true.))
    do i = 2, size(curve)
      nl=nl+1; lines(nl) = '                    ' // trim(tao_curve_name(curve(i)%c, .true.))
    enddo
    nl=nl+1; lines(nl) = 'Plot.Graph.Curve:   ' // trim(tao_curve_name(c1))
    do i = 2, size(curve)
      nl=nl+1; lines(nl) = '                    ' // trim(tao_curve_name(curve(i)%c))
    enddo
    nl=nl+1; write (lines(nl), amt)  'data_source          = ', c1%data_source
    nl=nl+1; write (lines(nl), amt)  'data_index           = ', c1%data_index
    nl=nl+1; write (lines(nl), amt)  'data_type_x          = ', c1%data_type_x
    nl=nl+1; write (lines(nl), amt)  'data_type            = ', c1%data_type
    nl=nl+1; write (lines(nl), amt)  'legend_text          = ', c1%legend_text
    nl=nl+1; write (lines(nl), amt)  'ele_ref_name         = ', c1%ele_ref_name
    nl=nl+1; write (lines(nl), imt)  'ix_branch            = ', c1%ix_branch
    nl=nl+1; write (lines(nl), imt)  'ix_ele_ref           = ', c1%ix_ele_ref
    nl=nl+1; write (lines(nl), imt)  'ix_ele_ref_track     = ', c1%ix_ele_ref_track
    nl=nl+1; write (lines(nl), imt)  'ix_bunch             = ', c1%ix_bunch
    nl=nl+1; write (lines(nl), imt)  'ix_universe          = ', c1%ix_universe
    nl=nl+1; write (lines(nl), imt)  'symbol_every         = ', c1%symbol_every
    nl=nl+1; write (lines(nl), rmt)  'y_axis_scale_factor  = ', c1%y_axis_scale_factor
    nl=nl+1; write (lines(nl), lmt)  'use_y2               = ', c1%use_y2
    nl=nl+1; write (lines(nl), lmt)  'draw_line            = ', c1%draw_line
    nl=nl+1; write (lines(nl), lmt)  'draw_symbols         = ', c1%draw_symbols
    nl=nl+1; write (lines(nl), lmt)  'draw_symbol_index    = ', c1%draw_symbol_index
    nl=nl+1; write (lines(nl), lmt)  'smooth_line_calc     = ', c1%smooth_line_calc
    nl=nl+1; write (lines(nl), iamt) 'line%width           = ', c1%line%width
    nl=nl+1; write (lines(nl), iamt) 'line%color           = ', c1%line%color, qp_color_name(c1%line%color)
    nl=nl+1; write (lines(nl), iamt) 'line%style           = ', c1%line%style, qp_line_style_name(c1%line%style)
    nl=nl+1; write (lines(nl), iamt) 'symbol%type          = ', c1%symbol%type, qp_symbol_type_name(c1%symbol%type)
    nl=nl+1; write (lines(nl), f3mt) 'symbol%height        = ', c1%symbol%height
    nl=nl+1; write (lines(nl), iamt) 'symbol%fill_pattern  = ', c1%symbol%fill_pattern, qp_fill_name(c1%symbol%fill_pattern)
    nl=nl+1; write (lines(nl), iamt) 'symbol%line_width    = ', c1%symbol%line_width
    
    if (show_sym) then
      n = nl + size(c1%x_symb) + 10
      if (n > size(lines)) call re_allocate(lines, n, .false.)
      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = 'Symbol points:'
      nl=nl+1; lines(nl) = '      i  index             x             y'
      err = .false.
      do j = 2, size(curve)
        if (size(curve(j)%c%y_symb) /= size(c1%y_symb)) then
          nl=nl+1; lines(nl) = 'NUMBER OF SYMBOL POINTS NOT THE SAME IN ALL CURVES!'
          err = .true.
          exit
        endif
      enddo
      if (.not. err) then
        do i = 1, size(c1%x_symb)
          nl=nl+1; write (lines(nl), '(2i7, 10es14.6)') i, c1%ix_symb(i), &
                      c1%x_symb(i), [ (curve(j)%c%y_symb(i), j = 1, size(curve)) ]
        enddo
      endif
    endif

    if (show_line) then
      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = 'Smooth line points:'
      nl=nl+1; lines(nl) = '             x             y'
      do j = 2, size(curve)
        if (size(curve(j)%c%y_line) /= size(c1%y_line)) then
          nl=nl+1; lines(nl) = 'NUMBER OF LINE POINTS NOT THE SAME IN ALL CURVES!'
          err = .true.
          exit
        endif
      enddo
      if (.not. err) then
        call re_allocate (lines, nl+size(c1%x_line)+100, .false.)
        do i = 1, size(c1%x_line)
          nl=nl+1; write (lines(nl), '(10es14.6)') c1%x_line(i), &
                                        [ (curve(j)%c%y_line(i), j = 1, size(curve)) ]
        enddo
      endif
    endif

  else
    nl=1; lines(1) = 'THIS IS NOT A CURVE NAME'
    return
  endif

  result_id = show_what

!----------------------------------------------------------------------
! data

case ('data')


  ! If just "show data" then show all names

  call tao_pick_universe (word1, line1, picked_uni, err)
  if (err) return

  if (line1 == ' ') then  ! just specified a universe

    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), '(a, t40, a)') 'Name', 'Using for Optimization'

    do iu = lbound(s%u, 1), ubound(s%u, 1)

      if (.not. picked_uni(iu)) cycle

      u => s%u(iu)

      do i = 1, size(u%d2_data)
        d2_ptr => u%d2_data(i)
        if (d2_ptr%name == ' ') cycle
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

    result_id = 'data:'
    return
  endif

  ! get pointers to the data

  call tao_find_data (err, word1, d2_ptr, d1_array, d_array)
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
    nl=nl+1; write(lines(nl), amt)  '%ele_name          = ', d_ptr%ele_name
    nl=nl+1; write(lines(nl), amt)  '%ele_start_name    = ', d_ptr%ele_start_name
    nl=nl+1; write(lines(nl), amt)  '%ele_ref_name      = ', d_ptr%ele_ref_name
    nl=nl+1; write(lines(nl), amt)  '%data_type         = ', d_ptr%data_type
    nl=nl+1; write(lines(nl), amt)  '%data_source       = ', d_ptr%data_source
    nl=nl+1; write(lines(nl), imt)  '%ix_branch         = ', d_ptr%ix_branch
    nl=nl+1; write(lines(nl), imt)  '%ix_ele            = ', d_ptr%ix_ele
    nl=nl+1; write(lines(nl), imt)  '%ix_ele_start      = ', d_ptr%ix_ele_start
    nl=nl+1; write(lines(nl), imt)  '%ix_ele_ref        = ', d_ptr%ix_ele_ref
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
    nl=nl+1; write(lines(nl), rmt)  '%invalid           = ', d_ptr%invalid_value
    nl=nl+1; write(lines(nl), rmt)  '%s                 = ', d_ptr%s
    nl=nl+1; write(lines(nl), amt)  '%merit_type        = ', d_ptr%merit_type
    nl=nl+1; write(lines(nl), rmt)  '%merit             = ', d_ptr%merit
    nl=nl+1; write(lines(nl), rmt)  '%delta_merit       = ', d_ptr%delta_merit
    nl=nl+1; write(lines(nl), rmt)  '%weight            = ', d_ptr%weight
    nl=nl+1; write(lines(nl), lmt)  '%exists            = ', d_ptr%exists
    nl=nl+1; write(lines(nl), lmt)  '%good_model        = ', d_ptr%good_model
    nl=nl+1; write(lines(nl), lmt)  '%good_design       = ', d_ptr%good_design
    nl=nl+1; write(lines(nl), lmt)  '%good_base         = ', d_ptr%good_base 
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
    ! Expressions generally have very long strings so we let this spill over to
    ! the where0 and where fields

    n_name  = 9      ! Set mimimum field widths
    n_start = 10
    n_ref   = 8
    n_ele   = 8

    do i = 1, size(d_array)
      d_ptr => d_array(i)%d
      if (.not. d_ptr%exists) cycle
      name = tao_datum_type_name(d_ptr)
      if (d_ptr%data_type(1:11) /= 'expression:') then
        n_name  = max(n_name,  len_trim(name))
        n_start = max(n_start, len_trim(d_ptr%ele_start_name))
        n_ref   = max(n_ref,   len_trim(d_ptr%ele_ref_name))
        n_ele   = max(n_ele,   len_trim(d_ptr%ele_name))
      endif
    enddo

    ! Write header
    ! Element names are left justified and real quantities are right justified

    line1 = ''; line2 = ''
    n=9+n_name;    line2(n:) = 'Ref_Ele'  ! n = i4 + 2x + n_name + 2x + 1
    n=n+n_ref+2;   line2(n:) = 'Start_Ele'
    n=n+n_start+2; line2(n:) = 'Ele'
    n=n+n_ele+10;  line2(n:) = 'Meas         Model        Design | Opt  Plot'
                   line1(n:) = '                                 |   Useit'

    nl=nl+1; lines(nl) = line1
    nl=nl+1; lines(nl) = line2

    ! if a range is specified, show the data range   

    call re_allocate (lines, nl+100+size(d_array), .false.)

    fmt  = '(i4, 4(2x, a), 3es14.4, 2l6)'
    fmt2 = '(4x, 4(2x, a), 3es14.4, 2l6)'

    do i = 1, size(d_array)
      d_ptr => d_array(i)%d
      if (.not. d_ptr%exists) cycle
      name = tao_datum_type_name(d_ptr)
      if (d_ptr%data_type(1:11) == 'expression:') then
        nl=nl+1; write(lines(nl), fmt) d_ptr%ix_d1, trim(name)
        nl=nl+1; write(lines(nl), fmt2) blank_str(1:n_name), &
                     d_ptr%ele_ref_name(1:n_ref), d_ptr%ele_start_name(1:n_start), &
                     d_ptr%ele_name(1:n_ele), d_ptr%meas_value, d_ptr%model_value, &
                     d_ptr%design_value, d_ptr%useit_opt, d_ptr%useit_plot
      else
        nl=nl+1; write(lines(nl), fmt) d_ptr%ix_d1, name(1:n_name), &
                     d_ptr%ele_ref_name(1:n_ref), d_ptr%ele_start_name(1:n_start), &
                     d_ptr%ele_name(1:n_ele), d_ptr%meas_value, d_ptr%model_value, &
                     d_ptr%design_value, d_ptr%useit_opt, d_ptr%useit_plot
      endif
    enddo

    nl=nl+1; lines(nl) = line2
    nl=nl+1; lines(nl) = line1

  ! else we must have a valid d2_ptr.

  elseif (associated(d2_ptr)) then

    call re_allocate (lines, nl+100+size(d2_ptr%d1), .false.)

    nl=nl+1; write(lines(nl), '(t40, a)')     'Using' 

    do i = 1, size(d2_ptr%d1)
      if (size(lines) < nl + 50) call re_allocate (lines, nl+100, .false.)
      call location_encode(line, d2_ptr%d1(i)%d%useit_opt, &
                      d2_ptr%d1(i)%d%exists, lbound(d2_ptr%d1(i)%d, 1))
      nl=nl+1; write(lines(nl), '(2x, 2a, i0, a, i0, a, t40, a)') &
                  trim(tao_d2_d1_name(d2_ptr%d1(i))), '[', lbound(d2_ptr%d1(i)%d, 1), &
                  ':', ubound(d2_ptr%d1(i)%d, 1), ']', trim(line)
    enddo

    if (any(d2_ptr%descrip /= ' ')) then
      call re_allocate (lines, nl+100+size(d2_ptr%descrip), .false.)
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
    nl=1; lines(1) = 'TRY BEING MORE SPECIFIC.'
    return
  endif

  result_id = show_what

!----------------------------------------------------------------------
! derivative

case ('derivative')

  call tao_find_data (err, word1, d_array = d_array)
  if (err) return

  call string_trim(stuff2(ix_word+1:), stuff2, ix_word)
  call tao_find_var(err, stuff2(:ix_word), v_array = v_array) 
  if (err) return

  found = .false.
  do id = 1, size(d_array)
    do iv = 1, size(v_array)
      d_ptr => d_array(id)%d
      v_ptr => v_array(iv)%v
      u => s%u(d_ptr%d1%d2%ix_uni)
      jd = d_ptr%ix_dmodel
      jv = v_ptr%ix_dvar
      if (jd /= 0 .and. jv /= 0) then
        nl=nl+1; write (lines(nl), '(2a20, es14.5, 2i5)') tao_datum_name(d_ptr), &
                                  tao_var1_name(v_ptr), u%dModel_dVar(jd, jv), jd, jv
        if (size(lines) < nl + 1) call re_allocate (lines, nl+200, .false.)
        found = .true.
      endif
    enddo
  enddo

  if (.not. found) then
    nl=nl+1; write (lines(nl), '(a)') 'No Derivatives'
  endif
  
  result_id = show_what

!----------------------------------------------------------------------
! ele

case ('element')

  print_taylor = .false.
  print_wig_terms = .false.
  print_all = .false.
  print_data = .false.

  do
    call tao_next_switch (stuff2, ['-taylor        ', '-wig_terms     ', &
                                   '-all_attributes', '-data          '], switch, err, ix)
    if (err) return
    if (switch == '') exit
    if (switch == '-taylor') print_taylor = .true.
    if (switch == '-wig_terms') print_wig_terms = .true.
    if (switch == '-all_attributes') print_all = .true.
    if (switch == '-data') print_data = .true.
  enddo

  call str_upcase (ele_name, stuff2)
  call tao_pick_universe (ele_name, ele_name, picked_uni, err, ix_u)
  if (err) return

  ! Wildcard: show all elements.

  if (index(ele_name, '*') /= 0 .or. index(ele_name, '%') /= 0 .or. &
                (ele_name(1:2) /= 'S:' .and. index(ele_name, ':') /= 0) .or. &
                count(picked_uni) > 1) then

    n_tot = 0
    do i_uni = lbound(s%u, 1), ubound(s%u, 1)
      if (.not. picked_uni(i_uni)) cycle
      call tao_locate_elements (ele_name, i_uni,eles, err, .true.)
      if (err) return
      lat => s%u(i_uni)%model%lat
      n_tot = n_tot + size(eles)
      do i = 1, size(eles)
        ele => eles(i)%ele
        if (size(lines) < nl+100) call re_allocate (lines, nl+200, .false.)
        if (count(picked_uni) > 1) then
          nl=nl+1; write (lines(nl), '(i8, 2x, i0, 2a)') ele%ix_ele, i_uni, '@', ele%name
        else
          nl=nl+1; write (lines(nl), '(i8, 2x, a)') ele%ix_ele, ele%name
        endif
      enddo
    enddo

    deallocate(eles)
    nl=nl+1; write (lines(nl), '(a, i0)') 'Number of Matches: ', n_tot

    if (nl == 0) then
      lines(1) = '*** No Matches to Name Found ***'
      return
    endif

    result_id = 'element:*'
    return

  endif

  ! No wildcard case...
  ! Normal: Show the element info

  call tao_locate_elements (ele_name, ix_u, eles, err)
  if (err) return
  ele => eles(1)%ele

  ! Show data associated with this element

  if (print_data) then
    call show_ele_data (u, ele, lines, nl)
    result_id = 'element:data'
    return
  endif

  if (tao_com%common_lattice) then
    s%u(ix_u)%lattice_recalc = .true.
    call tao_lattice_calc (ok)
  endif
  call type2_ele (ele, ptr_lines, n, print_all, 6, print_taylor, s%global%phase_units, &
            .true., s%u(ix_u)%model%lat, .true., .true., print_wig_terms)

  if (size(lines) < nl+n+100) call re_allocate (lines, nl+n+100, .false.)
  lines(nl+1:nl+n) = ptr_lines(1:n)
  nl = nl + n
  deallocate (ptr_lines)

  if (show_what == 'element') then
    nl=nl+1; lines(nl) = '[Conversion from Global to Screen: (Z, X) -> (-X, -Y)]'
  endif

  orb = u%model%lat_branch(eles(1)%ele%ix_branch)%orbit(eles(1)%ele%ix_ele)
  fmt = '(2x, a, 3p2f15.8)'
  lines(nl+1) = ' '
  lines(nl+2) = 'Orbit: [mm, mrad]'
  write (lines(nl+3), fmt) "X  X':", orb%vec(1:2)
  write (lines(nl+4), fmt) "Y  Y':", orb%vec(3:4)
  write (lines(nl+5), fmt) "Z  Z':", orb%vec(5:6)
  nl = nl + 5

  found = .false.
  do i = 2, size(eles)
    if (size(lines) < nl+2) call re_allocate (lines, nl+10, .false.)
    if (found) then
      nl=nl+1; lines(nl) = ''
      found = .true.
    endif
    nl=nl+1
    if (eles(i)%ele%ix_branch == 0) then
      write (lines(nl), '(a, i0)') &
              'Note: Found another element with same name at: ', eles(i)%ele%ix_ele
    else
      write (lines(nl), '(a, i0, a,i0)') &
              'Note: Found another element with same name at: ', &
              eles(i)%ele%ix_branch, '.', eles(i)%ele%ix_ele
    endif
  enddo

  result_id = show_what

!----------------------------------------------------------------------
! global

case ('global')

  print_global = .true.
  print_optimization = .false.
  print_bmad_com = .false.
  print_csr_param = .false.

  do
    call tao_next_switch (stuff2, ['-optimization', '-bmad_com    ', &
                                   '-csr_param   '], switch, err, ix)
    if (err) return
    if (switch == '') exit
    print_global = .false.
    if (switch == '-optimization') print_optimization = .true.
    if (switch == '-bmad_com') print_bmad_com = .true.
    if (switch == '-csr_param') print_csr_param = .true.
  enddo

  if (print_global) then
    nl=nl+1; lines(nl) = 'Global parameters:'
    nl=nl+1; write (lines(nl), imt) '  %bunch_to_plot                 = ', s%global%bunch_to_plot
    nl=nl+1; write (lines(nl), lmt) '  %label_lattice_elements        = ', s%global%label_lattice_elements
    nl=nl+1; write (lines(nl), lmt) '  %label_keys                    = ', s%global%label_keys
    nl=nl+1; write (lines(nl), amt) '  %phase_units                   = ', &
                                                    frequency_units_name(s%global%phase_units)
    nl=nl+1; write (lines(nl), lmt) '  %plot_on                       = ', s%global%plot_on
    nl=nl+1; write (lines(nl), lmt) '  %lattice_calc_on               = ', s%global%lattice_calc_on
    nl=nl+1; write (lines(nl), lmt) '  %command_file_print_on         = ', s%global%command_file_print_on
    nl=nl+1; write (lines(nl), lmt) '  %beam_timer_on                 = ', s%global%beam_timer_on
    nl=nl+1; write (lines(nl), lmt) '  %init_lats_with_rf_off         = ', s%global%init_lats_with_rf_off
    nl=nl+1; write (lines(nl), amt) '  %prompt_string                 = ', s%global%prompt_string
    nl=nl+1; write (lines(nl), amt) '  %print_command                 = ', s%global%print_command
    nl=nl+1; write (lines(nl), amt) '  %random_engine                 = ', s%global%random_engine
    nl=nl+1; write (lines(nl), amt) '  %random_gauss_converter        = ', s%global%random_gauss_converter
    nl=nl+1; write (lines(nl), rmt) '  %random_sigma_cutoff           = ', s%global%random_sigma_cutoff
    nl=nl+1; write (lines(nl), imt) '  %random_seed                   = ', s%global%random_seed
    if (s%global%random_seed == 0) then
      call ran_seed_get(ix)
      nl=nl+1; write (lines(nl), imt) '   random_seed (generated)      = ', ix
    endif
    nl=nl+1; write (lines(nl), amt) '  %track_type                    = ', s%global%track_type
    nl=nl+1; write (lines(nl), imt) '  %u_view                        = ', s%global%u_view
    nl=nl+1; write (lines(nl), lmt) '  %var_limits_on                 = ', s%global%var_limits_on
    nl=nl+1; write (lines(nl), amt) '  %var_out_file                  = ', s%global%var_out_file
    nl=nl+1; write (lines(nl), rmt) '  %y_axis_plot_dmin              = ', s%global%y_axis_plot_dmin

    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Internal Tao Parameters:'
    nl=nl+1; write (lines(nl), imt) 'Universe index range:          = ', lbound(s%u, 1), ubound(s%u, 1)
    nl=nl+1; write (lines(nl), lmt) 'common_lattice                 = ', tao_com%common_lattice
    nl=nl+1; write (lines(nl), amt) 'tao_com%beam_all_file          = ', tao_com%beam_all_file
    nl=nl+1; write (lines(nl), amt) 'tao_com%beam0_file             = ', tao_com%beam0_file
    nl=nl+1; write (lines(nl), lmt) 'tao_com%combine_consecutive_elements_of_like_name = ', &
                                                tao_com%combine_consecutive_elements_of_like_name
    nl=nl+1; write (lines(nl), amt) 'tao_com%init_lat_file          = ', tao_com%init_lat_file
    nl=nl+1; write (lines(nl), amt) 'tao_com%init_tao_file          = ', tao_com%init_tao_file
    nl=nl+1; write (lines(nl), imt) 'Number paused command files    = ', count(tao_com%cmd_file%paused)
  endif

  if (print_optimization) then
    call show_opt()
  endif

  if (print_bmad_com) then
    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Bmad_com Parameters:'
    nl=nl+1; write (lines(nl), imt) '  %taylor_order              = ', bmad_com%taylor_order
    nl=nl+1; write (lines(nl), lmt) '  %auto_bookkeeper           = ', bmad_com%auto_bookkeeper
    nl=nl+1; write (lines(nl), lmt) '  %trans_space_charge_on     = ', bmad_com%trans_space_charge_on
    nl=nl+1; write (lines(nl), lmt) '  %coherent_synch_rad_on     = ', bmad_com%coherent_synch_rad_on
    nl=nl+1; write (lines(nl), lmt) '  %spin_tracking_on          = ', bmad_com%spin_tracking_on
    nl=nl+1; write (lines(nl), lmt) '  %radiation_damping_on      = ', bmad_com%radiation_damping_on
    nl=nl+1; write (lines(nl), lmt) '  %radiation_fluctuations_on = ', bmad_com%radiation_fluctuations_on
    nl=nl+1; write (lines(nl), lmt) '  %spin_tracking_on          = ', bmad_com%spin_tracking_on
  endif

  if (print_csr_param) then
    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'CSR_param Parameters:'
    nl=nl+1; write(lines(nl), rmt) '  %ds_track_step        = ', csr_param%ds_track_step
    nl=nl+1; write(lines(nl), imt) '  %n_bin                = ', csr_param%n_bin
    nl=nl+1; write(lines(nl), imt) '  %particle_bin_span    = ', csr_param%particle_bin_span
    nl=nl+1; write(lines(nl), lmt) '  %lcsr_component_on    = ', csr_param%lcsr_component_on
    nl=nl+1; write(lines(nl), lmt) '  %lsc_component_on     = ', csr_param%lsc_component_on
    nl=nl+1; write(lines(nl), lmt) '  %tsc_component_on     = ', csr_param%tsc_component_on
    nl=nl+1; write(lines(nl), lmt) '  %ix1_ele_csr          = ', csr_param%ix1_ele_csr
    nl=nl+1; write(lines(nl), lmt) '  %ix2_ele_csr          = ', csr_param%ix2_ele_csr
  endif

  result_id = show_what

!----------------------------------------------------------------------
! graph

case ('graph')

  ! Look for switches

  call string_trim(stuff, stuff2, ix)

  ! Find particular graph

  call tao_find_plots (err, stuff2, 'BOTH', graph = graph, always_allocate = .true.)
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

    nl=nl+1; write (lines(nl), rmt) 'x_axis_scale_factor   = ', g%x_axis_scale_factor
    nl=nl+1; write (lines(nl), rmt) 'x%max                 = ', g%x%max
    nl=nl+1; write (lines(nl), rmt) 'x%min                 = ', g%x%min
    nl=nl+1; write (lines(nl), imt) 'x%major_div           = ', g%x%major_div
    nl=nl+1; write (lines(nl), imt) 'x%major_div_nominal   = ', g%x%major_div_nominal
    nl=nl+1; write (lines(nl), imt) 'x%places              = ', g%x%places
    nl=nl+1; write (lines(nl), lmt) 'x%draw_label          = ', g%x%draw_label
    nl=nl+1; write (lines(nl), lmt) 'x%draw_numbers        = ', g%x%draw_numbers

    nl=nl+1; write (lines(nl), lmt) 'y2_mirrors_y          = ', g%y2_mirrors_y
    nl=nl+1; write (lines(nl), rmt) 'y%max                 = ', g%y%max
    nl=nl+1; write (lines(nl), rmt) 'y%min                 = ', g%y%min
    nl=nl+1; write (lines(nl), imt) 'y%major_div           = ', g%y%major_div
    nl=nl+1; write (lines(nl), imt) 'y%major_div_nominal   = ', g%y%major_div_nominal
    nl=nl+1; write (lines(nl), imt) 'y%places              = ', g%y%places
    nl=nl+1; write (lines(nl), lmt) 'y%draw_label          = ', g%y%draw_label
    nl=nl+1; write (lines(nl), lmt) 'y%draw_numbers        = ', g%y%draw_numbers

    nl=nl+1; write (lines(nl), rmt) 'y2%max                = ', g%y2%max
    nl=nl+1; write (lines(nl), rmt) 'y2%min                = ', g%y2%min
    nl=nl+1; write (lines(nl), imt) 'y2%major_div          = ', g%y2%major_div
    nl=nl+1; write (lines(nl), imt) 'y2%major_div_nominal  = ', g%y2%major_div_nominal
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
    nl=1; lines(1) = 'This is not a graph'
    return
  endif

  result_id = show_what

!----------------------------------------------------------------------
! hom

case ('hom')

  nl=nl+1; lines(nl) = &
        '       #        Freq         R/Q           Q   m  Polarization_Angle'
  do i = 1, ubound(lat%ele, 1)
    ele => lat%ele(i)
    if (ele%key /= lcavity$) cycle
    if (ele%slave_status == multipass_slave$) cycle
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

  result_id = show_what

!----------------------------------------------------------------------
! keys

case ('key_bindings')

  call tao_key_info_to_str (1, 1, size(s%key), str, header)
  nl=nl+1; lines(nl) = ' Ix  ' // header

  do i = 1, size(s%key)
    call tao_key_info_to_str (i, 1, size(s%key), str, header)
    nl=nl+1; write (lines(nl), '(i3, 2x, a)') i, str
  enddo

  result_id = show_what

!----------------------------------------------------------------------
! lattice

case ('lattice')
  
  limited = .false.
  all_lat = .false.
  at_ends = .true.
  by_s = .false.
  print_header_lines = .true.
  print_tail_lines = .true.
  replacement_for_blank = ''
  ix_branch = 0
  undef_str = '---'
  show_lords = .false.
  show_custom = .false.
  column(:)%name = ""
  column(:)%label = ""

  ! get command line switches

  do
    call tao_next_switch (stuff2, [ &
        '-branch           ', '-blank_replacement', '-lords            ', '-middle           ', &
        '-all_tracking     ', '-0undef           ', '-no_label_lines   ', '-no_tail_lines    ', &
        '-custom           ', '-s                '], switch, err, ix)
    if (err) return
    if (switch == '') exit
    select case (switch)

    case ('-branch')
      call string_trim(stuff2(ix+1:), stuff2, ix)
      read (stuff2(1:ix), *, iostat = ios) ix_branch
      if (ios /= 0 .or. ix_branch < 0 .or. ix_branch > ubound(u%model%lat%branch, 1)) then
        nl=1; write (lines(1), *) 'Branch index out of bounds:', ix_branch
        return
      endif

    case ('-blank_replacement')
      call string_trim(stuff2(ix+1:), stuff2, ix)
      replacement_for_blank = stuff2(1:ix)

    case ('-lords')
      show_lords = .true.

    case ('-middle')
      at_ends = .false.

    case ('-all_tracking')
      all_lat = .true. 

    case ('-0undef')
      undef_str = '  0'

    case ('-no_label_lines')
      print_header_lines = .false.
      print_tail_lines = .false.

    case ('-no_tail_lines')
      print_tail_lines = .false.

    case ('-custom')
      show_custom = .true.
      call string_trim(stuff2(ix+1:), stuff2, ix)
      file_name = stuff2(1:ix)
      iu = lunget()
      open (iu, file = file_name, status = 'old', iostat = ios)
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT OPEN FILE: ' // file_name
        return
      endif
      column(:)%name = ""
      column(:)%label = ""
      read (iu, nml = custom_show_list, iostat = ios)
      close (iu)
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ "CUSTOM_SHOW_LIST" NAMELIST IN FILE: ' // file_name
        return
      endif

    case ('-s')
      by_s = .true.
    end select

  enddo
  
  branch => lat%branch(ix_branch)

  ! Construct columns if needed.

  if (.not. show_custom) then
    if (show_lords) then
      column(1)  = show_lat_column_struct('#',                   'i6',        6, '')
      column(2)  = show_lat_column_struct('x',                   'x',         2, '')
      column(3)  = show_lat_column_struct('ele::#[name]',        'a',         0, '')
      column(4)  = show_lat_column_struct('ele::#[key]',         'a16',      16, '')
      column(5)  = show_lat_column_struct('ele::#[s]',           'f10.3',    10, '')
      column(6)  = show_lat_column_struct('x',                   'x',         2, '')
      column(7)  = show_lat_column_struct("ele::#[lord_status]", 'a16',      16, '') 
    else
      column(1)  = show_lat_column_struct('#',                 'i6',        6, '')
      column(2)  = show_lat_column_struct('x',                 'x',         2, '')
      column(3)  = show_lat_column_struct('ele::#[name]',      'a',         0, '')
      column(4)  = show_lat_column_struct('ele::#[key]',       'a16',      16, '')
      column(5)  = show_lat_column_struct('ele::#[s]',         'f10.3',    10, '')
      column(6)  = show_lat_column_struct('ele::#[l]',         'f8.3',      8, '')
      column(7)  = show_lat_column_struct('ele::#[beta_a]',    'f7.2',      7, 'beta|  a')
      column(8)  = show_lat_column_struct('ele::#[phi_a]',     'f8.3',      8, '')
      column(9)  = show_lat_column_struct('ele::#[eta_a]',     'f5.1',      5, '')
      column(10) = show_lat_column_struct('ele::#[orbit_x]',   '3p, f8.3',  8, '')
      column(11) = show_lat_column_struct('ele::#[beta_b]',    'f7.2',      7, '')
      column(12) = show_lat_column_struct('ele::#[phi_b]',     'f8.3',      8, '')
      column(13) = show_lat_column_struct('ele::#[eta_b]',     'f5.1',      5, '')
      column(14) = show_lat_column_struct('ele::#[orbit_y]',   '3p, f8.3',  8, '')
    endif
  endif

  ! Find elements to use

  if (allocated (picked_ele)) deallocate (picked_ele)
  allocate (picked_ele(0:branch%n_ele_max))
  if (show_lords) then
    picked_ele(0:branch%n_ele_track) = .false.
    picked_ele(branch%n_ele_track+1:branch%n_ele_max) = .true.
  else
    picked_ele(0:branch%n_ele_track) = .true.
    picked_ele(branch%n_ele_track+1:branch%n_ele_max) = .false.
  endif

  if (by_s) then
    ix = index(stuff2, ':')
    if (ix == 0) then
      nl=1; lines(1) = 'NO ":" FOUND FOR RANGE SELECTION'
      return
    endif
    read (stuff2(1:ix-1), *, iostat = ios1) s1
    read (stuff2(ix+1:), *, iostat = ios2) s2
    if (ios1 /= 0 .or. ios2 /= 0) then
      nl=1; lines(1) = 'ERROR READING RANGE SELECTION: ' // stuff2
      return
    endif

    picked_ele = .false.
    do ie = 1, branch%n_ele_track
      if (at_ends) then;
        s_ele = branch%ele(ie)%s
      else
        s_ele = (branch%ele(ie-1)%s + branch%ele(ie)%s) / 2
      endif
      if (s_ele >= s1 .and. s_ele <= s2) picked_ele(ie) = .true.
    enddo

  elseif (stuff2(1:ix) == '*' .or. all_lat) then
    ! picked_ele already set

  elseif (ix /= 0) then
    call tao_locate_elements (stuff2, u%ix_uni, eles, err, .true.)
    if (err) return
    picked_ele = .false.
    do i = 1, size(eles)
      picked_ele(eles(i)%ele%ix_ele) = .true.
    enddo

  elseif (.not. show_lords) then
    if (size(picked_ele) > 200) then
      picked_ele(201:) = .false.
      limited = .true.
    endif

  endif

  !

  if (at_ends) then
    write (line1, '(6x, a)') 'Model values at End of Element:'
  else
    write (line1, '(6x, a)') 'Model values at Center of Element:'
  endif

  ! Setup columns

  ix1 = 1
  line2 = ""
  line3 = ""
  do i = 1, size(column)
    if (column(i)%name == "") cycle

    ! Convert from old 'dat::' format to 'lat::' format.
    ix = index(column(i)%name, 'dat::')
    if (ix /= 0) column(i)%name = column(i)%name(1:ix-1) // 'lat' // column(i)%name(ix+3:)

    column(i)%format = '(' // trim(column(i)%format) // ')'

    if (column(i)%field_width == 0) then
      if (column(i)%name /= 'ele::#[name]') then
        call out_io (s_error$, r_name, &
            'FIELD_WIDTH = 0 CAN ONLY BE USED WITH "ele::#[name]" TYPE COLUMNS')
        return
      endif
      column(i)%field_width = 5
      do ie = 0, branch%n_ele_max
        if (.not. picked_ele(ie)) cycle
        column(i)%field_width = max(column(i)%field_width, len_trim(branch%ele(ie)%name)+1)
      enddo
    endif

    ix2 = ix1 + column(i)%field_width

    if (column(i)%label == '') then
      name = column(i)%name
      if (index(name, 'ele::') /= 0) then
        i1 = index(name, '[')
        i2 = index(name, ']')
        name = name(i1+1:i2-1)
      elseif (index(name, 'lat::') /= 0) then
        ix = index(name, 'lat::')
        name = name(ix+5:)
        i2 = index(name, '[')
        name = name(1:i2-1)
      elseif  (name == '#') then
        line2(ix2-5:) = 'Index'
        ix1 = ix2
        cycle    
      elseif (name == 'x') then
        ix1 = ix2
        cycle    
      else
        name = ''
      endif

      ix = index(name, '_')
      n = len_trim(name)
      if (column(i)%format(2:2) == 'a') then
        line2(ix1:) = name
      elseif (ix == 0) then
        line2(ix2-n:) = name
      else
        line2(ix2-ix+1:) = name(1:ix-1)
        line3(ix2-n+ix:) = name(ix+1:)
      endif

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

    ix1 = ix2
  enddo

  ! Collect lines

  if (print_header_lines) then
    nl=nl+1; lines(nl) = line1
    nl=nl+1; lines(nl) = line2
    if (line3 /= '') then
      nl=nl+1; lines(nl) = line3
    endif
  endif

  do ie = 0, branch%n_ele_max
    if (.not. picked_ele(ie)) cycle
    if (size(lines) < nl+100) call re_allocate (lines, nl+200, .false.)
    line = ''
    nc = 1
    ele => branch%ele(ie)
    do i = 1, size(column)
      name = column(i)%name
      if (name == '') cycle

      if (name == '#') then
        write (line(nc:), column(i)%format, iostat = ios) ie

      elseif (name == 'ele::#[name]') then
        write (line(nc:), column(i)%format, iostat = ios) ele%name

      elseif (name == 'ele::#[key]') then
        write (line(nc:), column(i)%format, iostat = ios) key_name(ele%key)

      elseif (name == 'ele::#[slave_status]') then
        write (line(nc:), column(i)%format, iostat = ios) control_name(ele%slave_status)

      elseif (name == 'ele::#[lord_status]') then
        write (line(nc:), column(i)%format, iostat = ios) control_name(ele%lord_status)

      elseif (name == 'ele::#[type]') then
        if (ele%type == '') then
          write (line(nc:), column(i)%format, iostat = ios) replacement_for_blank
        else
          write (line(nc:), column(i)%format, iostat = ios) ele%type
        endif

      elseif (name == 'x') then
        ios = 0

      else
        write (nam, '(i0, a, i0)') ix_branch, '>>', ie
        call str_substitute (name, '#', trim(nam))
        ix = index(name, 'ele::')
        if (.not. at_ends .and. ix /= 0) then
          name = name(:ix+2) // '_mid' // trim(name(ix+3:))
        endif
        call tao_evaluate_expression (name, 1, .false., value, good, err, .false.)
        if (err .or. .not. allocated(value) .or. size(value) /= 1) then
          n = len(undef_str)
          k = min(n, column(i)%field_width - 1)
          j = nc + column(i)%field_width - k
          line(j:) = undef_str(n-k+1:n)
        else
          if (index(column(i)%format, 'i') /= 0 .or. index(column(i)%format, 'I') /= 0) then
            write (line(nc:), column(i)%format, iostat = ios) nint(value(1))
          else
            write (line(nc:), column(i)%format, iostat = ios) value(1)
          endif
        endif
      endif

      if (ios /= 0) then
        lines(1) = 'WIDTH TOO SMALL FOR NUMBER OR BAD FORMAT: ' // column(i)%format
        lines(2) = 'FOR DISPLAYING: ' // column(i)%name
        nl = 2
        return
      endif

      nc  = nc + column(i)%field_width

    enddo

    nl=nl+1; lines(nl) = line

  enddo

  if (print_tail_lines) then
    nl=nl+1; lines(nl) = line2
    if (line3 /= '') then
      nl=nl+1; lines(nl) = line3
    endif
    nl=nl+1; lines(nl) = line1
  endif

  if (limited) then
    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Note: Since no range is given, the number of elements shown is limited to 200.'
  endif

  deallocate(picked_ele)

  result_id = show_what

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
  call show_opt
  call out_io (s_blank$, r_name, lines(1:nl))
  nl = 0

  result_id = show_what

!----------------------------------------------------------------------
! particle

case ('orbit')

  call tao_locate_elements (word1, u%ix_uni, eles, err)
  if (err) return
  do i = 1, 6
    nl=nl+1; write (lines(nl), rmt) '     ', &
                u%model%lat_branch(eles(1)%ele%ix_branch)%orbit%vec(eles(1)%ele%ix_ele)
  enddo

  result_id = show_what

!----------------------------------------------------------------------
! particle

case ('particle')

  nb = s%global%bunch_to_plot
  ix_branch = 0

  call tao_next_switch (stuff2, ['-lost'], switch, err, ix)
  if (err) return

  if (switch == '-lost') then
    if (ix /= 0) then
      read (stuff2(1:ix), *, iostat = ios) nb
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'CANNOT READ BUNCH INDEX.')
        return
      endif
    endif
    bunch => u%uni_branch(ix_branch)%ele(lat%n_ele_track)%beam%bunch(nb)
    nl=nl+1; write (lines(nl), *) 'Bunch:', nb
    nl=nl+1; lines(nl) = 'Particles lost at:'
    nl=nl+1; lines(nl) = '    Ix Ix_Ele  Ele_Name '
    do i = 1, size(bunch%particle)
      if (bunch%particle(i)%ix_lost == not_lost$) cycle
      if (nl == size(lines)) call re_allocate (lines, nl+100, .false.)
      ie = bunch%particle(i)%ix_lost
      nl=nl+1; write (lines(nl), '(i6, i7, 2x, a)') i, ie, lat%ele(ie)%name
    enddo
    result_id = 'particle:lost'
    return
  endif

  ix = index(stuff2, '.')
  if (ix > 0 .and. ix < ix_word) then
    read (stuff2(1:ix-1), *, iostat = ios) nb
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'CANNOT READ BUNCH INDEX.')
      return
    endif
    call string_trim (stuff2(ix+1:), stuff2, ix_word)
  endif

  read (stuff2, *, iostat = ios) ix_p
  if (ios /= 0) then
    call out_io (s_error$, r_name, 'CANNOT READ PARTICLE INDEX')
    return
  endif

  ix_ele = 0
  call string_trim(stuff2(ix_word+1:), stuff2, ix_word)
  ele_name = stuff2(:ix_word)
  if (ele_name /= '') then
    call tao_pick_universe (ele_name, ele_name, picked_uni, err, ix_u)
    if (err) return
    call tao_locate_elements (ele_name, ix_u, eles, err)
    if (err) return
    ix_ele = eles(1)%ele%ix_ele
    ix_branch = eles(2)%ele%ix_ele
  endif

  uni_branch => u%uni_branch(ix_branch)

  if (.not. allocated(uni_branch%ele(ix_ele)%beam%bunch)) then
    call out_io (s_error$, r_name, 'BUNCH NOT ASSOCIATED WITH THIS ELEMENT.')
    return
  endif

  if (nb < 1 .or. nb > size(uni_branch%ele(ix_ele)%beam%bunch)) then
    call out_io (s_error$, r_name, 'BUNCH INDEX OUT OF RANGE: \i0\ ', i_array = [ nb ])
    return
  endif

  bunch => uni_branch%ele(ix_ele)%beam%bunch(nb)

  if (ix_p < 1 .or. ix_p > size(bunch%particle)) then
    call out_io (s_error$, r_name, 'PARTICLE INDEX OUT OF RANGE: \i0\ ', i_array = [ ix_p ])
    return
  endif

  nl=nl+1; write (lines(nl), imt) 'At lattice element:', ix_ele
  nl=nl+1; write (lines(nl), imt) 'Bunch:    ', nb
  nl=nl+1; write (lines(nl), imt) 'Particle: ', ix_p
  nl=nl+1; write (lines(nl), lmt) 'Is Alive? ', bunch%particle(ix_p)%ix_lost == not_lost$
  nl=nl+1; write (lines(nl), lmt) 'Coords: '
  nl=nl+1; write (lines(nl), '(a, 6es13.5)') '  ', bunch%particle(ix_p)%r%vec

  result_id = show_what

!----------------------------------------------------------------------
! plot

case ('plot')

  ! Look for switches

  show_shape = .false.

  do
    call tao_next_switch (stuff2, ['-shapes'], switch, err, ix)
    if (err) return
    if (switch == '') exit
    if (switch == '-shapes') show_shape = .true.
  enddo

  ! Element shapes

  if (show_shape) then

    nl=nl+1; lines(nl) = ' '
    nl=nl+1; lines(nl) = 'Floor_plan Element Shapes:'
    nl=nl+1; lines(nl) = &
          'Ele_Name                        Shape         Color        dy_pix   Label    '
    nl=nl+1; lines(nl) = &
          '----------------------------    --------      -----        -------  ---------'

    do i = 1, size(tao_com%ele_shape_floor_plan)
      shape => tao_com%ele_shape_floor_plan(i)
      if (shape%ele_name == '') cycle
      nl=nl+1; write (lines(nl), '(3a, f10.4, 2x, a8)') &
                shape%ele_name(1:32), shape%shape(1:14), shape%color(1:10), &
                shape%dy_pix, shape%label_type
    enddo

    nl=nl+1; lines(nl) = ' '
    nl=nl+1; lines(nl) = 'Lat_layout Element Shapes:'
    nl=nl+1; lines(nl) = &
          'Ele_Name                        Shape         Color        dy_pix   Label    '
    nl=nl+1; lines(nl) = &
          '----------------------------    --------      -----        -------  ---------'

    do i = 1, size(tao_com%ele_shape_lat_layout)
      shape => tao_com%ele_shape_lat_layout(i)
      if (shape%ele_name == '') cycle
      nl=nl+1; write (lines(nl), '(3a, f10.4, 2x, a8)') &
                shape%ele_name(1:32), shape%shape(1:14), shape%color(1:10), &
                shape%dy_pix, shape%label_type
    enddo

    result_id = 'plot:shape'
    return 
   
  endif

  ! stuff2 is blank => print overall info

  if (stuff2 == ' ') then

    nl=nl+1; lines(nl) = 'plot_page parameters:'
    nl=nl+1; write (lines(nl), rmt)  '%size                       = ', s%plot_page%size       
    nl=nl+1; write (lines(nl), imt)  '%n_curve_pts                = ', s%plot_page%n_curve_pts
    nl=nl+1; write (lines(nl), f3mt) '%text_height                = ', s%plot_page%text_height 
    nl=nl+1; write (lines(nl), f3mt) '%main_title_text_scale      = ', s%plot_page%main_title_text_scale 
    nl=nl+1; write (lines(nl), f3mt) '%graph_title_text_scale     = ', s%plot_page%graph_title_text_scale 
    nl=nl+1; write (lines(nl), f3mt) '%axis_number_text_scale     = ', s%plot_page%axis_number_text_scale 
    nl=nl+1; write (lines(nl), f3mt) '%axis_label_text_scale      = ', s%plot_page%axis_label_text_scale 
    nl=nl+1; write (lines(nl), f3mt) '%key_table_text_scale       = ', s%plot_page%key_table_text_scale 
    nl=nl+1; write (lines(nl), f3mt) '%legend_text_scale          = ', s%plot_page%legend_text_scale 
    nl=nl+1; write (lines(nl), f3mt) '%shape_height_max           = ', s%plot_page%shape_height_max  

    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Templates:'
    nl=nl+1; lines(nl) = '   Plot                .Graph'
    nl=nl+1; lines(nl) = '   ------------------  ----------'
    do i = 1, size(s%template_plot)
      p => s%template_plot(i)
      if (p%name == '') cycle
      if (allocated(p%graph)) then
          nl=nl+1; write (lines(nl), '(3x, a, 100(2a, 2x))') &
                            p%name(1:20), ('.', trim(p%graph(j)%name), j = 1, size(p%graph))
      else
        nl=nl+1; write (lines(nl), '(3x, a)') p%name 
      endif
    enddo

    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = '                    ' // &
                         '                                    Location on Page'
    nl=nl+1; lines(nl) = 'Visible  Plot Region' // &
                         '         <-->  Template             x1    x2    y1    y2'  
    nl=nl+1; lines(nl) = '-------  -----------' // &
                         '               -----------------------------------------'
    do i = 1, size(s%plot_region)
      region => s%plot_region(i)
      if (region%name == '') cycle
      nl=nl+1; write (lines(nl), '(3x l1, 5x, a20, a, a18, 4f6.2)') region%visible, &
                                    region%name, '<-->  ', region%plot%name, region%location
    enddo

    result_id = 'plot:'
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
    nl=nl+1; write (lines(nl), imt) 'x%major_div          = ', p%x%major_div
    nl=nl+1; write (lines(nl), imt) 'x%major_div_nominal  = ', p%x%major_div_nominal
    nl=nl+1; write (lines(nl), imt) 'x%places             = ', p%x%places
    nl=nl+1; write (lines(nl), lmt) 'x%draw_label         = ', p%x%draw_label
    nl=nl+1; write (lines(nl), lmt) 'x%draw_numbers       = ', p%x%draw_numbers
    nl=nl+1; write (lines(nl), lmt) 'autoscale_x          = ', p%autoscale_x
    nl=nl+1; write (lines(nl), lmt) 'autoscale_y          = ', p%autoscale_y
    nl=nl+1; write (lines(nl), lmt) 'autoscale_gang_x     = ', p%autoscale_gang_x
    nl=nl+1; write (lines(nl), lmt) 'autoscale_gang_y     = ', p%autoscale_gang_y
    
    nl=nl+1; lines(nl) = 'Graphs:'
    do i = 1, size(p%graph)
      nl=nl+1; write (lines(nl), amt) '   ', p%graph(i)%name
    enddo

  else
    nl=1; lines(1) = 'This is not a name of a plot'
    return
  endif

  result_id = show_what

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
    nl=1; lines(1) = 'UNKNOWN SWITCH: ' // stuff2
    return
  endif

  result_id = show_what

!----------------------------------------------------------------------
! taylor_map

case ('taylor_map')

  by_s = .false.
  n_order = -1

  do
    call tao_next_switch (stuff2, [ '-order', '-s' ], switch, err, ix)
    if (err) return
    if (switch == '') exit
    select case (switch)
    case ('-s')
      by_s = .true.
    case ('-order')
      read (stuff2(:ix), *, iostat = ios) n_order
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ ORDER NUMBER!'
        return
      endif
      call string_trim (stuff2(ix+1:), stuff2, ix)
    end select
  enddo

  ele1 = stuff2(:ix)
  call string_trim(stuff2(ix+1:), stuff2, ix)
  ele2 = stuff2(:ix)
  if (stuff2(ix+1:) /= '') then
    nl=1; lines(1) = 'EXTRA STUFF ON LINE!'
    return
  endif

  ! By s

  if (by_s) then
    if (ele1 == '') then
      s2 = lat%ele(lat%n_ele_track)%s
      s1 = 0
    else
      read (ele1, *, iostat = ios) s1
      if (ios /= 0) then
        nl=1; lines(1) = 'BAD S1 VALUE:' // ele1
        return
      endif
    endif

    if (ele2 == '') then
      s2 = s1
      s1 = 0
    else 
      read (ele2, *, iostat = ios) s2
      if (ios /= 0) then
        nl=1; lines(1) = 'BAD S2 VALUE:' // ele2
        return
      endif
    endif

    call transfer_map_calc_at_s (lat, taylor, s1, s2)

  ! By element

  else

    if (ele1 == '') then
      ix2 = lat%n_ele_track
      ix1 = 0
    else
      call tao_locate_elements (ele1, u%ix_uni, eles, err)
      if (size(eles) > 1) then
        nl=1; lines(1) = 'MULTIPLE ELEMENTS BY THIS NAME: ' // ele1
        return
      endif
      if (err .or. size(eles) == 0) return
      ix1 = eles(1)%ele%ix_ele
    endif

    if (ele2 == '') then
      ix2 = ix1
      ix1 = 0
    else
      call tao_locate_elements (ele2, u%ix_uni, eles, err)
      if (size(eles) > 1) then
        nl=1; lines(1) = 'MULTIPLE ELEMENTS BY THIS NAME: ' // ele2
        return
      endif
      if (err .or. size(eles) == 0) return
      ix2 = eles(1)%ele%ix_ele
    endif

    call transfer_map_calc (lat, taylor, ix1, ix2)

  endif

  ! Print results

  if (n_order > -1) call truncate_taylor_to_order (taylor, n_order, taylor)
  call type2_taylors (taylor, lines, nl)

  result_id = show_what

!----------------------------------------------------------------------
! tune

case ('tune')

  nl=nl+1; lines(nl) = 'Use "show universe" instead.'

  result_id = show_what

!----------------------------------------------------------------------
! universe
    
case ('universe')

  if (len_trim(word1) > 1 .and. index('-connections', trim(word1)) == 1) then
    do i = lbound(s%u, 1), ubound(s%u, 1)
      u => s%u(i)
      nl=nl+1; write (lines(nl), imt) 'For universe:', i
      if (u%connect%connected) then
        n = u%connect%from_uni_ix_ele
        ix_u = u%connect%from_uni
        nl=nl+1; write (lines(nl), imt) '  Injection from universe:', ix_u
        nl=nl+1; write (lines(nl), '(3a, i0, a)') '  From element:    ', &
                                      trim(s%u(ix_u)%model%lat%ele(n)%name), ' (# ', n, ')'
        nl=nl+1; write (lines(nl), '(a, f10.2)')  '  From s:          ', u%connect%from_uni_s
        nl=nl+1; write (lines(nl), '(a, l1)')     '  Match_to_design: ', u%connect%match_to_design
      else
        nl=nl+1; write (lines(nl), '(a)') '  No injection into this universe.'
      endif
      if (u%connect%to_uni > 0) then
        nl=nl+1; write (lines(nl), '(a, i0)') '  Injection from this universe to universe: ', &
                                                                           u%connect%to_uni 
      else
        nl=nl+1; write (lines(nl), '(a)') '  No injection from this universe'
      endif
      nl=nl+1; lines(nl) = ''
    enddo
    if (nl > 1) nl = nl - 1  ! Erase last blank line.
    result_id = show_what
    return
  endif



  if (word1 == ' ') then
    ix_u = s%global%u_view
  else
    read (word1, *, iostat = ios) ix_u
    if (ios /= 0) then
      nl=1; lines(1) = 'BAD UNIVERSE NUMBER'
      return
    endif
    if (ix_u < lbound(s%u, 1) .or. ix_u > ubound(s%u, 1)) then
      nl=1; lines(1) = 'UNIVERSE NUMBER OUT OF RANGE'
      return
    endif
  endif

  u => s%u(ix_u)
  ix_branch = 0
  uni_branch => u%uni_branch(ix_branch)
  branch => lat%branch(ix_branch)

  nl = 0
  nl=nl+1; write (lines(nl), imt) 'Universe: ', ix_u
  nl=nl+1; write (lines(nl), imt) 'Branch:   ', ix_branch
  nl=nl+1; write (lines(nl), imt) '%n_d2_data_used        = ', u%n_d2_data_used
  nl=nl+1; write (lines(nl), imt) '%n_data_used           = ', u%n_data_used
  nl=nl+1; write (lines(nl), lmt) '%do_synch_rad_int_calc = ', u%do_synch_rad_int_calc
  nl=nl+1; write (lines(nl), lmt) '%do_chrom_calc         = ', u%do_chrom_calc
  nl=nl+1; write (lines(nl), lmt) '%calc_beam_emittance   = ', u%calc_beam_emittance
  nl=nl+1; write (lines(nl), lmt) '%mat6_recalc_on        = ', u%mat6_recalc_on
  nl=nl+1; write (lines(nl), lmt) '%is_on                 = ', u%is_on
  nl=nl+1; write (lines(nl), amt) '%beam0_file            = ', trim(u%beam_info%beam0_file)
  nl=nl+1; write (lines(nl), amt) '%beam_all_file         = ', trim(u%beam_info%beam_all_file)
  nl=nl+1; write (lines(nl), amt) '%beam_saved_at:        = ', trim(u%beam_saved_at)
  nl=nl+1; lines(nl) = ''
  nl=nl+1; write(lines(nl), amt) 'Lattice name:    ', lat%lattice
  nl=nl+1; write(lines(nl), amt) 'Input_file_name: ', lat%input_file_name
  nl=nl+1; lines(nl) =           'Lattice Type:    ' // lattice_type(branch%param%lattice_type)
  nl=nl+1; write (lines(nl), imt) &
                'Elements used in tracking: From 1 through ', branch%n_ele_track
  if (branch%n_ele_max .gt. branch%n_ele_track) then
    nl=nl+1; write (lines(nl), '(a, i0, a, i0)') 'Lord elements:   ', &
                      branch%n_ele_track+1, '  through ', branch%n_ele_max
  else
    nl=nl+1; write (lines(nl), '(a)') 'There are NO Lord elements'
  endif

  nl=nl+1; write (lines(nl), '(a, f0.3)')   'Lattice length:             ', branch%param%total_length
  nl=nl+1; write (lines(nl), lmt)           'Aperture limits on?:        ', branch%param%aperture_limit_on

  if (branch%param%lattice_type == linear_lattice$ .and. branch%param%ix_lost /= not_lost$) then
    if (s%global%track_type == 'beam') then
      nl=nl+1; write (lines(nl), '(a, i0)') 'Tracking: Lost beam at:     ', branch%param%ix_lost
    else
      nl=nl+1; write (lines(nl), '(a, i0)') 'Tracking: Lost particle at: ', branch%param%ix_lost
    endif
  endif

  if (.not. branch%param%stable) then
    nl=nl+1; write (lines(nl), '(a, l)') 'Model lattice stability: ', branch%param%stable
    nl=nl+1; write (lines(nl), '(a, l)') 'Design lattice stability:', u%design%lat%param%stable
    result_id = 'universe:unstable'
    return
  endif
 
  call radiation_integrals (lat, &
                     u%model%lat_branch(0)%orbit, u%model%modes, u%ix_rad_int_cache)
  call radiation_integrals (u%design%lat, &
                     u%design%lat_branch(0)%orbit, u%design%modes, u%ix_rad_int_cache)
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
  if (lat%param%lattice_type == circular_lattice$) then
    call calc_z_tune(u%model%lat)
    nl=nl+1; write (lines(nl), '(1x, a12, 2f11.4, 3x, a)') 'Z_tune:', &
         -u%model%lat%z%tune/twopi, -u%design%lat%z%tune/twopi, '! The design value is calculated with RF on'
  endif
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

  result_id = show_what

!----------------------------------------------------------------------
! variable
    
case ('use')  

  nl=nl+1; lines(nl) = 'veto data *@*'
  nl=nl+1; lines(nl) = ''

  do i = lbound(s%u, 1), ubound(s%u, 1)
    do j = 1, s%u(i)%n_d2_data_used
      d2_ptr => s%u(i)%d2_data(j)
      call re_allocate (lines, nl+size(d2_ptr%d1)+10, .false.)
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

  call re_allocate (lines, nl+size(s%v1_var)+10, .false.)
  do i = 1, size(s%v1_var)
    v1_ptr => s%v1_var(i)
    if (v1_ptr%name == ' ') cycle
    call re_allocate (lines, nl+200, .false.)
    call location_encode (line, v1_ptr%v%useit_opt, v1_ptr%v%exists, lbound(v1_ptr%v, 1))
    nl=nl+1; write (lines(nl), '(5a)') 'use var ', trim(v1_ptr%name), '[', trim(line), ']'
  enddo

  result_id = show_what

!----------------------------------------------------------------------
! variable

case ('value')

  call tao_evaluate_expression (stuff2, 0, .false., value, good, err)
  if (err) return

  if (size(value) == 1) then
    nl=nl+1; write (lines(nl), '(3x, es17.8)') value(1)
  else
    call re_allocate (lines, size(value)+100, .false.)
    do i = 1, size(value)
      nl=nl+1; write (lines(nl), '(i4, a, es17.8)') i, ':  ', value(i)
    enddo
  endif

  result_id = show_what

!----------------------------------------------------------------------
! variable

case ('variable')

  good_opt_only = .false.
  bmad_format = .false.
  do
    call tao_next_switch (stuff2, [ '-bmad_format  ', '-good_opt_only' ], switch, err, ix)
    if (err) return
    if (switch == '') exit
    if (switch == '-bmad_format') bmad_format = .true.
    if (switch == '-good_opt_only') good_opt_only = .true.
  enddo
  
  if (.not. allocated (s%v1_var)) then
    nl=1; lines(1) = 'NO VARIABLES HAVE BEEN DEFINED IN THE INPUT FILES!'
    return 
  endif

  ! Bmad format

  if (bmad_format) then
    call tao_var_write (' ', good_opt_only)
    result_id = 'variable:bmad'
    return
  endif

  ! If 'n@' is present then write out stuff for universe n

  ix = index(word1, '@')
  if (ix /= 0) then
    if (ix == 1) then
      ix_u = s%global%u_view
    else
      read (word1(:ix-1), *, iostat = ios) ix_u
      if (ios /= 0) then
        nl=1; lines(1) = 'BAD UNIVERSE NUMBER'
        return
      endif
      if (ix_u == 0) ix_u = s%global%u_view
      if (ix_u < 1 .or. ix_u > ubound(s%u, 1)) then
        nl=1; lines(1) = 'UNIVERSE NUMBER OUT OF RANGE'
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
    result_id = 'variable:@'
    return
  endif

! If just "show var" then show all names

  if (word1 == ' ') then
    nl=nl+1; write (lines(nl), '(7x, a, t50, a)') 'Name', 'Using for Optimization'
    do i = 1, size(s%v1_var)
      v1_ptr => s%v1_var(i)
      if (v1_ptr%name == ' ') cycle
      call re_allocate (lines, nl+200, .false.)
      call location_encode (line, v1_ptr%v%useit_opt, v1_ptr%v%exists, lbound(v1_ptr%v, 1))
      nl=nl+1; write(lines(nl), '(i5, 2x, 2a, i0, a, i0, a, t50, a)') v1_ptr%ix_v1, &
                      trim(v1_ptr%name), '[', lbound(v1_ptr%v, 1), ':', &
                      ubound(v1_ptr%v, 1), ']', trim(line)
      
    enddo

    result_id = 'variable:'
    return
  endif

! are we looking at a range of locations?

  call tao_find_var(err, word1, v1_array, v_array) 
  if (err) return
  n_size = 0
  if (allocated(v_array)) n_size = size(v_array)

! v_ptr is valid then show the variable info.

  if (n_size == 1) then

    v_ptr => v_array(1)%v

    nl=nl+1; write(lines(nl), amt)  '%ele_name         = ', v_ptr%ele_name
    nl=nl+1; write(lines(nl), amt)  '%attrib_name      = ', v_ptr%attrib_name 
    nl=nl+1; write(lines(nl), imt)  '%ix_attrib        = ', v_ptr%ix_attrib 
    nl=nl+1; write(lines(nl), imt)  '%ix_var           = ', v_ptr%ix_var
    nl=nl+1; write(lines(nl), imt)  '%ix_dvar          = ', v_ptr%ix_dvar           
    nl=nl+1; write(lines(nl), imt)  '%ix_v1            = ', v_ptr%ix_v1
    nl=nl+1; write(lines(nl), rmt)  '%model            = ', v_ptr%model_value
    nl=nl+1; write(lines(nl), rmt)  '%base             = ', v_ptr%base_value

    if (.not. allocated (v_ptr%this)) then
      nl=nl+1; write(lines(nl), imt)  'this(:) -- Not associated!'
    else
      do i = 1, size(v_ptr%this)
        nl=nl+1; write(lines(nl), iimt)  '%this(', i, ')%ix_uni:        ', &
                                                            v_ptr%this(i)%ix_uni
        nl=nl+1; write(lines(nl), iimt)  '%this(', i, ')%ix_branch:     ', v_ptr%this(i)%ix_branch
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

    nl=nl+1; write(lines(nl), rmt)  '%design           = ', v_ptr%design_value
    nl=nl+1; write(lines(nl), rmt)  '%old              = ', v_ptr%old_value
    nl=nl+1; write(lines(nl), rmt)  '%meas             = ', v_ptr%meas_value
    nl=nl+1; write(lines(nl), rmt)  '%ref              = ', v_ptr%ref_value
    nl=nl+1; write(lines(nl), rmt)  '%correction       = ', v_ptr%correction_value
    nl=nl+1; write(lines(nl), rmt)  '%high_lim         = ', v_ptr%high_lim
    nl=nl+1; write(lines(nl), rmt)  '%low_lim          = ', v_ptr%low_lim
    nl=nl+1; write(lines(nl), rmt)  '%step             = ', v_ptr%step
    nl=nl+1; write(lines(nl), rmt)  '%weight           = ', v_ptr%weight
    nl=nl+1; write(lines(nl), rmt)  '%delta_merit      = ', v_ptr%delta_merit
    nl=nl+1; write(lines(nl), amt)  '%merit_type       = ', v_ptr%merit_type
    nl=nl+1; write(lines(nl), rmt)  '%merit            = ', v_ptr%merit
    nl=nl+1; write(lines(nl), rmt)  '%dmerit_dvar      = ', v_ptr%dMerit_dVar
    nl=nl+1; write(lines(nl), imt)  '%ix_key_table     = ', v_ptr%ix_key_table
    if (v_ptr%ix_key_table > 0 ) then
      nl=nl+1; write(lines(nl), rmt)  '%key_val0         = ', v_ptr%key_val0
      nl=nl+1; write(lines(nl), rmt)  '%key_delta        = ', v_ptr%key_delta
    endif
    nl=nl+1; write(lines(nl), lmt)  '%exists           = ', v_ptr%exists
    nl=nl+1; write(lines(nl), lmt)  '%good_var         = ', v_ptr%good_var
    nl=nl+1; write(lines(nl), lmt)  '%good_user        = ', v_ptr%good_user
    nl=nl+1; write(lines(nl), lmt)  '%good_opt         = ', v_ptr%good_opt
    nl=nl+1; write(lines(nl), lmt)  '%useit_opt        = ', v_ptr%useit_opt
    nl=nl+1; write(lines(nl), lmt)  '%useit_plot       = ', v_ptr%useit_plot

    result_id = 'variable:1:' // word1

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
      call re_allocate (lines, nl+200, .false.)
      nl=nl+1
      write(lines(nl), '(i6, 2x, a)') v_ptr%ix_v1, tao_var_attrib_name(v_ptr)
      write(lines(nl)(nc+9:), '(3es14.4, 7x, l)') v_ptr%meas_value, &
                 v_ptr%model_value, v_ptr%design_value, v_ptr%useit_opt
    enddo
    nl=nl+1; lines(nl) = line1

  else
    nl=1; lines(1) = '???'
    result_id = 'variable:?'
  endif

  result_id = show_what

!----------------------------------------------------------------------
! wave

case ('wave')

  nl=nl+1; write(lines(nl), '(a, 2i4)') 'ix_a:', s%wave%ix_a1, s%wave%ix_a2
  nl=nl+1; write(lines(nl), '(a, 2i4)') 'ix_b:', s%wave%ix_b1, s%wave%ix_b2

  select case (s%wave%data_type)
  case ('orbit.x', 'orbit.y', 'eta.x', 'eta.y', 'beta.a', 'beta.b')
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'A Region Sigma_Fit/Amp_Fit:  ', s%wave%rms_rel_a
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'B Region Sigma_Fit/Amp_Fit:  ', s%wave%rms_rel_b
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_Kick/Kick: ', s%wave%rms_rel_k
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_phi:       ', s%wave%rms_phi
    nl=nl+1; lines(nl) = ' '
    if (s%wave%data_type(1:4) == 'beta') then
      nl=nl+1; lines(nl) = 'Normalized Kick = kick * beta  [urad * meter]'
    else
      nl=nl+1; lines(nl) = 'Normalized Kick = kick * sqrt(beta)  [urad * sqrt(meter)]'
    endif
    nl=nl+1; lines(nl) = 'After Dat#    Norm_K       phi'
    do i = 1, min(s%wave%n_kick, 10)
      nl=nl+1; write (lines(nl), '(i9, f12.2, 1f10.3)') s%wave%kick(i)%ix_dat, &
                  1e6*s%wave%kick(i)%amp, s%wave%kick(i)%phi
    enddo

  case ('phase.a', 'phase.b')
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'A Region Sigma_Fit/Amp_Fit:  ', s%wave%rms_rel_a
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'B Region Sigma_Fit/Amp_Fit:  ', s%wave%rms_rel_b
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_Kick/Kick: ', s%wave%rms_rel_k
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_phi:       ', s%wave%rms_phi
    nl=nl+1; write(lines(nl), '(a, f8.3, a)') &
                                    'Chi_C:           ', s%wave%chi_c, ' [Figure of Merit]'
    nl=nl+1; lines(nl) = ' '
    nl=nl+1; lines(nl) = 'Normalized Kick = k * l * beta [dimensionless]'
    nl=nl+1; lines(nl) = '   where k = quadrupole gradient [rad/m^2].'
    nl=nl+1; lines(nl) = 'After Dat#     Norm_K       phi'
    do i = 1, min(s%wave%n_kick, 10)
      nl=nl+1; write (lines(nl), '(i9, f12.4, f10.3)') s%wave%kick(i)%ix_dat, &
                  s%wave%kick(i)%amp, s%wave%kick(i)%phi
    enddo

  case ('cbar.11', 'cbar.12', 'cbar.22')
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'A Region Sigma_+/Amp_+:  ', s%wave%rms_rel_as
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'A Region Sigma_-/Amp_-:  ', s%wave%rms_rel_ar
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'B Region Sigma_+/Amp_+:  ', s%wave%rms_rel_bs
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'B Region Sigma_-/Amp_-:  ', s%wave%rms_rel_br
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Kick |K+|   = ', 2*s%wave%amp_ba_s
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_K+/K+ = ', 2*s%wave%rms_rel_ks
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Kick |K-|   = ', 2*s%wave%amp_ba_r
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_K-/K- = ', 2*s%wave%rms_rel_kr
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_phi+:      ', s%wave%rms_phi_s
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_phi-:      ', s%wave%rms_phi_r
    nl=nl+1; write(lines(nl), '(a, f8.3, a)') &
                                    'Chi_a:           ', s%wave%chi_a, ' [Figure of Merit]'

    nl=nl+1; lines(nl) = 'After Dat#     Norm_K    phi+    phi-   phi_a   phi_b'
    do i = 1, min(s%wave%n_kick, 10)
      nl=nl+1; write (lines(nl), '(i11, f10.4, 4f8.3, 2f10.3)') &
            s%wave%kick(i)%ix_dat, &
            s%wave%kick(i)%amp, s%wave%kick(i)%phi_s, s%wave%kick(i)%phi_r, &
            (s%wave%kick(i)%phi_s+s%wave%kick(i)%phi_r)/2, &
            (s%wave%kick(i)%phi_s-s%wave%kick(i)%phi_r)/2
    enddo

  end select

  if (s%wave%n_kick > 10) then
    nl=nl+1; lines(nl) = ' etc...'
  endif


  result_id = show_what

!----------------------------------------------------------------------

case default

  nl=1; lines(1) = "INTERNAL ERROR, SHOULDN'T BE HERE!"

end select

!----------------------------------------------------------------------
!----------------------------------------------------------------------
contains

subroutine show_ele_data (u, ele, lines, nl)

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct), pointer :: datum
type (ele_struct) ele

character(*) :: lines(:)
character(100) l1
integer nl, i

logical found_one

!

nl=nl+1; write (lines(nl), '(a)') "  "
write (l1, '(a, 20x, a)') "Data Name", &
          "Data Type             |  Model Value  |  Design Value |  Base Value"
nl=nl+1; lines(nl) = l1

found_one = .false.
do i = 1, size(u%data)
  if (u%data(i)%ix_ele == ele%ix_ele .and. u%data(i)%ix_branch == ele%ix_branch) then
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

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! contains

subroutine show_opt ()

implicit none

nl=nl+1; lines(nl) = 'Global optimization parameters:'
nl=nl+1; write (lines(nl), rmt) '  %de_lm_step_ratio              = ', s%global%de_lm_step_ratio
nl=nl+1; write (lines(nl), rmt) '  %de_var_to_population_factor   = ', s%global%de_var_to_population_factor
nl=nl+1; write (lines(nl), rmt) '  %lm_opt_deriv_reinit           = ', s%global%lm_opt_deriv_reinit
nl=nl+1; write (lines(nl), rmt) '  %lmdif_eps                     = ', s%global%lmdif_eps
nl=nl+1; write (lines(nl), rmt) '  %merit_finish                  = ', s%global%merit_finish
nl=nl+1; write (lines(nl), imt) '  %n_top10                       = ', s%global%n_top10
nl=nl+1; write (lines(nl), imt) '  %n_opti_loops                  = ', s%global%n_opti_loops
nl=nl+1; write (lines(nl), imt) '  %n_opti_cycles                 = ', s%global%n_opti_cycles
nl=nl+1; write (lines(nl), lmt) '  %derivative_recalc             = ', s%global%derivative_recalc
nl=nl+1; write (lines(nl), lmt) '  %svd_retreat_on_merit_increase = ', s%global%svd_retreat_on_merit_increase 
nl=nl+1; write (lines(nl), lmt) '  %derivative_uses_design        = ', s%global%derivative_uses_design
nl=nl+1; write (lines(nl), lmt) '  %opt_with_ref                  = ', s%global%opt_with_ref 
nl=nl+1; write (lines(nl), lmt) '  %opt_with_base                 = ', s%global%opt_with_base
nl=nl+1; write (lines(nl), amt) '  %optimizer                     = ', s%global%optimizer
nl=nl+1; lines(nl) = ''
nl=nl+1; lines(nl) = 'opti_de_param Parameters:'
nl=nl+1; write (lines(nl), rmt) '  %CR                   = ', opti_de_param%CR
nl=nl+1; write (lines(nl), rmt) '  %F                    = ', opti_de_param%F
nl=nl+1; write (lines(nl), rmt) '  %l_best               = ', opti_de_param%l_best
nl=nl+1; write (lines(nl), lmt) '  %binomial_cross       = ', opti_de_param%binomial_cross
nl=nl+1; write (lines(nl), lmt) '  %use_2nd_diff         = ', opti_de_param%use_2nd_diff
nl=nl+1; write (lines(nl), lmt) '  %randomize_F          = ', opti_de_param%randomize_F
nl=nl+1; write (lines(nl), lmt) '  %minimize_merit       = ', opti_de_param%minimize_merit

end subroutine show_opt

end subroutine tao_show_this

end module
