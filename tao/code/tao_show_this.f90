!+
! Subroutine tao_show_this (what, result_id, lines, nl)
!
! Routine called by tao_show_cmd and tao_pipe_cmd. 
! This routine is not meant for general use.
!
! Input:
!   what      -- character(*): What to show.
!
! Output:
!   result_id -- character(*): ID string idendifying what was shown. Used by the Tao GUI.
!   lines(:)  -- character(*), allocatable: Output.
!   nl        -- integer: Number of lines(:).
!-

subroutine tao_show_this (what, result_id, lines, nl)

use tao_top10_mod, dummy => tao_show_this
use tao_c_interface_mod, only: tao_c_interface_com
use tao_command_mod, only: tao_next_switch, tao_next_word
use tao_version_mod
use location_encode_mod, only: location_encode
use transfer_map_mod, only: transfer_map_from_s_to_s, mat6_from_s_to_s
use opti_de_mod, only: opti_de_param
use wall3d_mod, only: calc_wall_radius
use twiss_and_track_mod, only: twiss_and_track_at_s
use ptc_spin, only: c_linear_map, assignment(=)
use pointer_lattice, only: operator(.sub.), operator(**), operator(*), alloc, kill, print, ci_phasor, assignment(=)
use ptc_layout_mod, only: ptc_emit_calc, lat_to_ptc_layout, type_ptc_fibre, assignment(=), taylor_inverse
use ptc_map_with_radiation_mod, only: ptc_rad_map_struct, ptc_setup_map_with_radiation, tree_element_zhe
use photon_target_mod, only: to_surface_coords
use expression_mod, only: expression_stack_to_string, split_expression_string
use beam_utils, only: calc_bunch_params_z_slice

implicit none

type (tao_universe_struct), pointer :: u
type (tao_lattice_branch_struct), pointer :: tao_branch, design_tao_branch
type (tao_lattice_struct), pointer :: tao_lat
type (tao_d2_data_array_struct), allocatable :: d2_array(:)
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_ptr
type (tao_data_struct), pointer :: d_ptr
type (tao_data_struct), target :: datum
type (tao_dynamic_aperture_struct), pointer :: da
type (tao_building_wall_section_struct), pointer :: section
type (tao_building_wall_point_struct), pointer :: pt
type (tao_v1_var_array_struct), allocatable, target :: v1_array(:)
type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_var_struct), pointer :: v_ptr
type (tao_var_array_struct), allocatable, target :: v_array(:)
type (tao_plot_array_struct), allocatable :: plot(:)
type (tao_graph_array_struct), allocatable :: graph(:)
type (tao_curve_array_struct), allocatable :: curve(:)
type (tao_plot_struct), pointer :: p
type (tao_graph_struct), pointer :: g
type (tao_curve_struct), pointer :: cc, c1
type (tao_plot_region_struct), pointer :: region
type (tao_d1_data_array_struct), allocatable :: d1_array(:)
type (tao_data_array_struct), allocatable :: d_array(:)
type (tao_ele_shape_struct), pointer :: shapes(:)
type (tao_ele_shape_struct), pointer :: shape
type (tao_shape_pattern_struct), pointer :: pattern
type (all_pointer_struct) a_ptr
type (beam_struct), pointer :: beam
type (beam_init_struct), pointer :: beam_init
type (tao_beam_branch_struct), pointer :: bb
type (lat_struct), pointer :: lat, design_lat
type (ele_struct), pointer :: ele, ele1, ele2, slave, lord, ele_ref
type (ele_struct), target :: ele3, ele0
type (em_field_struct) field, field0, field1
type (control_struct), pointer :: contl
type (controller_struct), pointer :: cntr
type (control_ramp1_struct), pointer :: rmp
type (ramper_lord_struct), pointer :: rml
type (bunch_struct), pointer :: bunch
type (wake_struct), pointer :: wake
type (wake_lr_mode_struct), pointer :: lr_mode
type (coord_struct), target :: orb, orb0, orb2, orbit
type (bunch_params_struct), target :: bunch_params
type (bunch_params_struct), pointer :: bunch_p
type (taylor_struct) taylor(6)
type (ele_pointer_struct), allocatable :: eles(:)
type (branch_struct), pointer :: branch, branch2, design_branch
type (tao_model_branch_struct), pointer :: model_branch
type (wall3d_struct), pointer :: wall
type (wall3d_section_struct), pointer :: wall_sec
type (wall3d_vertex_struct), pointer :: v
type (random_state_struct) ran_state
type (bmad_normal_form_struct), pointer :: bmad_nf
type (ptc_normal_form_struct), pointer :: ptc_nf
type (aperture_scan_struct), pointer :: da_scan
type (aperture_point_struct), pointer :: da_pt
type (tao_wave_kick_pt_struct), pointer :: wk
type (tao_spin_map_struct), pointer :: sm
type (normal_modes_struct) mode_ptc, mode_ptc_no_vert, mode_6d, mode_6d_no_vert
type (normal_modes_struct), pointer :: mode_m, mode_d
type (ptc_rad_map_struct), target :: rad_map
type (tree_element_zhe), pointer :: rmap(:)
type (tao_lat_sigma_struct), pointer :: lat_sig
type (track_struct), target :: track
type (track_point_struct), pointer :: tp
type (strong_beam_struct), pointer :: sb
type (c_taylor) ptc_ctaylor
type (complex_taylor_struct) bmad_ctaylor, ctaylor(3), ctaylor1
type (rad_map_ele_struct), pointer :: ri
type (grid_field_pt1_struct), pointer :: g_pt
type (tao_expression_info_struct), allocatable :: info(:)

type old_show_lat_column_struct
  character(80) :: name = ''
  character(40) :: format = ''
  integer :: width = 0
  character(40) :: label = ''
  logical :: remove_line_if_zero = .false.
  real(rp) :: scale_factor = 1
end type

type show_lat_column_struct
  character(80) :: name = ''
  character(40) :: format = ''
  character(40) :: label = ''
  logical :: remove_line_if_zero = .false.
  real(rp) :: scale_factor = 1
  integer :: width = 0
end type

type show_lat_column_info_struct
  integer :: attrib_type = -1
  character(40) :: attrib_name = '' ! Is Upper case
  real(rp) :: sum = 0
  real(rp) :: sum2 = 0
  real(rp) :: int_sum = 0
  real(rp) :: int_sum2 = 0
  real(rp) :: int_s = 0
  real(rp) :: val_last = real_garbage$
  integer :: n_sum = 0
  integer :: indent = 0   ! Spacing from start of line to beginning of column.
end type

type (old_show_lat_column_struct) column(60)
type (show_lat_column_struct) col(60)
type (show_lat_column_info_struct) col_info(60) 

real(rp) phase_units, gam, s_ele, s0, s1, s2, s3, val, z, z1, z2, z_in, s_pos, dt, angle, r
real(rp) sig_mat(6,6), sig0_mat(6,6), mat6(6,6), vec0(6), vec_in(6), vec3(3), l_lat
real(rp) pc, e_tot, value_min, value_here, pz1, phase, zb(2)
real(rp) g_vec(3), dr(3), v0(3), v2(3), g_bend, c_const, mc2, del, time1, ds, ref_vec(6), beta
real(rp) gamma, E_crit, E_ave, c_gamma, P_gam, N_gam, N_E2, H_a, H_b, rms, mean, s_last, s_now, n0(3)
real(rp) pz2, qs, q, x, xi_sum, xi_diff, dn_dpz(3), dn_partial(3,3), dn_partial2(3,3)
real(rp), allocatable :: value(:)

complex(rp) eval(6), evec(6,6), n_eigen(6,3)

character(*) :: what
character(*), parameter :: r_name = "tao_show_cmd"
character(*), allocatable :: lines(:)
character(*) result_id
character(n_char_show) line, line1, line2, line3
character(n_char_show) what2

character(1) delim
character(3) undef_str 
character(6) :: mode(4) = [character(6):: 'a-mode', 'b-mode', 'c-mode', 'spin']
character(9) angle_str
character(16) velocity_fmt, momentum_fmt, e_field_fmt, b_field_fmt, position_fmt, energy_fmt
character(16) spin_fmt, t_fmt, twiss_fmt, disp_fmt, str1, str2, where, output_type
character(24) show_name, show2_name, what_to_show
character(24) :: blank_str = '', phase_units_str, val_str
character(24) :: plane, imt, imt2, lmt, lmt2, amt, iamt, ramt, f3mt, rmt, rmt2, rmt3, irmt, iimt
character(40) ele_name, ele1_name, ele2_name, ele_ref_name, b_name, param_name, uni_str
character(40) replacement_for_blank, component, s_fmt
character(60) var_name, aname, myname, style
character(100) :: word1, word2, fmt, fmt2, fmt3, switch, why_invalid
character(200) header, str, attrib0, file_name, name, excite_zero(3), veto
character(200), allocatable :: alloc_lines(:)
character(400) attrib_list(20), attrib, sub_name

character(2), parameter :: q_name(0:3) = ['q0', 'qx', 'qy', 'qz']
character(1), parameter :: abc_name(1:3) = ['A', 'B', 'C']

character(20) :: show_what

integer data_number, ix_plane, ix_class, n_live, n_order, i0, i1, i2, ix_branch, width, expo(6)
integer nl, nl0, loc, ixl, iu, nc, n_size, ix_u, ios, ie, ig, nb, id, iv, jd, jv, stat, lat_type, print_control
integer ix, ix0, ix1, ix2, i, j, k, n, n_print, show_index, ju, ios1, ios2, i_uni, i_con, i_ic, ix_comb
integer num_locations, ix_ele, n_name, n_start, n_ele, n_ref, n_tot, ix_p, print_lords, ix_word, species
integer xfer_mat_print, twiss_out, ix_sec, n_attrib, ie0, a_type, ib, ix_min, n_remove, n_zeros_found
integer eval_pt, n_count, print_field, nt
integer, allocatable :: ix_c(:), ix_remove(:)

logical bmad_format, good_opt_only, print_wall, show_lost, logic, aligned, undef_uses_column_format, print_debug
logical err, found, first_time, by_s, print_header_lines, all_lat, limited, show_labels, do_calc, flip, show_energy
logical show_sym, show_line, show_shape, print_data, ok, print_tail_lines, print_slaves, print_super_slaves
logical show_all, name_found, print_taylor, print_rad, print_attributes, err_flag, angle_units, map_calc, clean
logical print_ptc, force_use_ptc, called_from_pipe_cmd, print_eigen, show_mat, show_q, print_rms, do_inverse
logical valid_value, print_floor, show_section, is_complex, print_header, print_by_uni, do_field, delim_found
logical, allocatable :: picked_uni(:), valid(:), picked2(:)
logical, allocatable :: picked_ele(:)

namelist / custom_show_list / column, col

!

call re_allocate (lines, 200)

err = .false.

lines = " "
nl = 0

rmt  = '(a, 9es16.8)'
rmt2 = '(a, es16.8, t40, es16.8)'
rmt3 = '(a, es16.8, 10a)'
f3mt = '(a, 9(f0.3, 2x))'
irmt = '(a, i0, a, es16.8)'
imt  = '(a, 9(i0, 2x))'
imt2 = '(a, i0, t40, a, i0)'
iimt = '(a, i0, a, i8)'
lmt  = '(a, 9(l1, 2x))'
lmt2 = '(a, l1, t40, a, l1)'
amt  = '(9a)'
iamt = '(a, i0, 2x, 9a)'
ramt = '(a, f0.3, 2x, 9a)'

ix_branch = s%global%default_branch
call point_to_this_uni(-1, ix_branch, u, lat, branch, model_branch, tao_branch, design_tao_branch)

phase_units = radians_to_angle_units(s%global%phase_units)
phase_units_str = short_angle_units_name(s%global%phase_units)

! find what to show

result_id = 'ERROR'

if (what == '') then
  nl=1; lines(1) = 'SHOW WHAT? (Use "help show" to see list)' 
  return
endif

call string_trim (what, what2, ix_word)
if (ix_word > 4 .and. index('plot_page', what2(:ix_word)) == 1)        what2 = 'plot -page ' // what2(ix_word+1:)
if (ix_word > 3 .and. index('bmad_com', what2(:ix_word)) == 1)         what2 = 'global -bmad_com ' // what2(ix_word+1:)
if (ix_word > 3 .and. index('ptc_com', what2(:ix_word)) == 1)          what2 = 'global -ptc_com ' // what2(ix_word+1:)
if (ix_word > 2 .and. index('csr_param', what2(:ix_word)) == 1)        what2 = 'global -csr_param ' // what2(ix_word+1:)
if (ix_word > 3 .and. index('space_charge_com', what2(:ix_word)) == 1) what2 = 'global -space_charge_com ' // what2(ix_word+1:)
if (ix_word > 2 .and. index('floor_plan', what2(:ix_word)) == 1)       what2 = 'plot -floor_plan ' // what2(ix_word+1:)

call match_word (what2, [character(20):: 'alias', 'beam', 'branch', 'building_wall', &
        'chromaticity', 'constraints', 'control', 'curve', 'data', 'debug', &
        'derivative', 'dynamic_aperture', 'element', 'emittance', 'field', 'global', 'graph', &
        'history', 'hom', 'internal', 'key_bindings', 'lattice', 'matrix', 'merit', 'normal_form', &
        'optimizer', 'orbit', 'particle', 'plot', 'ptc', 'radiation_integrals', 'rampers', 'spin', 'string', &
        'symbolic_numbers', 'taylor_map',  'top10', &
        'track', 'tune', 'twiss_and_orbit', 'universe', 'use', 'value', 'variables', 'version', &
        'wake_elements', 'wall', 'wave'], ix, matched_name = show_what)
if (ix == 0) then
  nl=1; lines(1) = 'SHOW WHAT? WORD NOT RECOGNIZED: ' // what
  return
endif

if (ix < 0) then
  nl=1; lines(1) = 'SHOW WHAT? AMBIGUOUS: ' // what
  return
endif

call string_trim (what2, what2, ix_word)
call string_trim (what2(ix_word+1:), what2, ix_word)
word1 = what2(:ix_word)
result_id = show_what  ! Default

select case (show_what)

!----------------------------------------------------------------------
! alias

case ('alias')

  call re_allocate (lines, s%com%n_alias+10, .false.)
  lines(1) = 'Aliases:'
  nl = 1
  do i = 1, s%com%n_alias
    nl=nl+1; lines(nl) = trim(s%com%alias(i)%name) // ' = "' // trim(s%com%alias(i)%expanded_str) // '"'
  enddo

!----------------------------------------------------------------------
! beam

case ('beam')

  ele_name = ''
  what_to_show = ''
  zb = -1

  do 
    call tao_next_switch (what2, [character(16):: '-universe', '-lattice', '-comb', '-z'], .true., switch, err)
    if (err) return
    if (switch == '') exit

    select case (switch)

    case ('-universe')
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios) ix
      u => tao_pointer_to_universe(ix)
      if (aname == '' .or. ios /= 0 .or. .not. associated(u)) then
        nl=1; lines(1) = 'CANNOT READ OR OUT-OF RANGE "-universe" ARGUMENT'
        return
      endif

    case ('-lattice', '-comb')
      what_to_show = switch

    case ('-z')
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios1) zb(1)
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios2) zb(2)

      if (ios1 /= 0 .or. ios2 /= 0) then
        nl=1; lines(1) = 'CANNOT READ Z VALUES.'
        return
      endif
      if (any(zb < 0) .or. any(zb > 1)) then
        nl=1; lines(1) = 'Z VALUES MUST BE IN RANGE [0, 1].'
        return
      endif
      if (zb(1) >= zb(2)) then
        nl=1; lines(1) = 'Z(2) MUST BE GREATER THAN Z(1).'
        return
      endif
  
    case default
      if (ele_name /= '') then
        nl=1; lines(1) = 'EXTRA STUFF ON THE COMMAND LINE: ' // switch
        return
      endif

      ele_name = switch

    end select
  enddo

  lat => u%model%lat
  branch => lat%branch(0)
  bb => u%model_branch(0)%beam

  !

  if (what_to_show == '-comb' .and. zb(1) >= 0) then
    nl=nl+1; lines(nl) = '"-comb" and "-z" cannot both be used at the same time'
    return
  endif

  if (zb(1) >= 0 .and. ele_name == '') then
    nl=nl+1; lines(nl) = 'Lattice element must be specified when using "-z" option.'
    return
  endif

  if (what_to_show == '-comb' .and. ele_name == '') then   ! -comb used without comb index.
    if (.not. allocated(tao_branch%bunch_params_comb)) then
      nl=nl+1; lines(nl) = 'Beam parameter comb not calculated (check comb_ds_save)' 
      return
    endif

    nl=nl+1; lines(nl) = ' Index           S  N_Live  N_lost'
    do i = 0, tao_branch%bunch_params_comb(1)%n_pt
      bunch_p => tao_branch%bunch_params_comb(1)%pt(i)
      nl=nl+1; write (lines(nl), '(i6, f12.3, 2i8)') i, bunch_p%s, bunch_p%n_particle_live, &
                                                                 bunch_p%n_particle_tot-bunch_p%n_particle_live
    enddo

    return
  endif

  !

  if (ele_name == '') then
    ele => branch%ele(bb%ix_track_start)
  elseif (what_to_show == '-comb') then
    read (ele_name, *, iostat = ios) ix_comb
    if (ios /= 0) then
      nl=1; lines(1) = 'CANNOT DECODE COMB INDEX: ' // ele_name
      return
    endif
  else
    call tao_pick_universe (ele_name, ele_name, picked_uni, err, ix_u)
    if (err) return
    u => s%u(ix_u)
    call tao_locate_elements (ele_name, ix_u, eles, err)
    if (err .or. size(eles) == 0) return
    ele => eles(1)%ele
    branch => ele%branch
  endif

  ix_branch = branch%ix_branch
  tao_branch => u%model%tao_branch(ix_branch)

  ! Sigma calc from lattice Twiss.

  if (what_to_show == '-lattice') then
    lat_sig => tao_branch%lat_sigma(ele%ix_ele)

    nl=nl+1; lines(nl) = 'Sigma calc from lattice Twiss at: ' // trim(ele%name) // ' ' //  ele_loc_name(ele, .false., '()')
    nl=nl+1; write(lines(nl), rmt)  '  S-position:                 ', ele%s
    nl=nl+1; write(lines(nl), imt)  '  In branch:                  ', ix_branch
    nl=nl+1; write(lines(nl), rmt) '  RMS:     ', &
                      sqrt(lat_sig%mat(1,1)), sqrt(lat_sig%mat(2,2)), sqrt(lat_sig%mat(3,3)), &
                      sqrt(lat_sig%mat(4,4)), sqrt(lat_sig%mat(5,5)), sqrt(lat_sig%mat(6,6))
    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Sigma Mat       x              px               y              py              z             pz'
    do i = 1, 6
      nl=nl+1; write (lines(nl), '(a2, 2x, 6es16.8)') coord_name(i), lat_sig%mat(i,:)
    enddo

    if (all(lat_sig%mat == 0)) then
      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = 'Note: Emittances are set in the beam_init structure.'
    endif

    result_id = 'beam:lat'
    return
  endif

  ! No element index

  if (ele_name == '') then
    if (.not. u%beam%track_beam_in_universe) then
      nl=nl+1; lines(nl) = 'Beam tracking not done in universe: ' // int_str(u%ix_uni)
      nl=nl+1; lines(nl) = 'Create a tao_beam_init namelist for this universe in the appropriate init file if beam tracking wanted.'
      return
    endif

    nl=nl+1; write(lines(nl), '(2(a, i0))') 'Universe: ', u%ix_uni, '  of: ', ubound(s%u, 1)
    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), amt) 'global%track_type           = ', quote(s%global%track_type)
    nl=nl+1; write(lines(nl), lmt) 'global%beam_timer_on        = ', s%global%beam_timer_on
    nl=nl+1; write(lines(nl), rmt) 'global%beam_dead_cutoff     = ', s%global%beam_dead_cutoff

    fmt = '(3a, i0, a)'
    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'General beam components (set by "set beam ..."):'
    nl=nl+1; write(lines(nl), lmt)  'always_reinit     = ', u%beam%always_reinit
    nl=nl+1; write(lines(nl), amt)  'saved_at          = ', quote(u%beam%saved_at)
    nl=nl+1; write(lines(nl), amt)  'dump_at           = ', quote(u%beam%dump_at)
    nl=nl+1; write(lines(nl), amt)  'dump_file         = ', quote(u%beam%dump_file)
    nl=nl+1; write(lines(nl), rmt3) 'comb_ds_save      = ', tao_branch%comb_ds_save, '  ! Note: -1 => Use (latice branch length)/plot_page%n_curve_pts'
    if (allocated(tao_branch%bunch_params_comb)) then
      nl=nl+1; write(lines(nl), amt)  'comb index range  = [0, ', int_str(tao_branch%bunch_params_comb(1)%n_pt), ']'
    endif
    nl=nl+1; write(lines(nl), amt)  'track_start       = ', quote(bb%track_start), '  ! ', ele_full_name(branch%ele(bb%ix_track_start))
    nl=nl+1; write(lines(nl), amt)  'track_end         = ', quote(bb%track_end),   '  ! ', ele_full_name(branch%ele(bb%ix_track_end))

    beam => u%model_branch(0)%ele(bb%ix_track_start)%beam
    if (allocated(beam%bunch)) then
      nl=nl+1; write(lines(nl), imt) 'n_particle (actual)       = ', size(beam%bunch(1)%particle)
      nl=nl+1; write(lines(nl), imt) 'n_bunch                   = ', size(beam%bunch)
      nl=nl+1; write(lines(nl), rmt) 'bunch_charge_tot          = ', beam%bunch(1)%charge_tot
      nl=nl+1; write(lines(nl), amt) 'bunch_species             = ', species_name(beam%bunch(1)%particle(1)%species)
    endif

    beam_init => bb%beam_init
    species = species_id(beam_init%species, default_tracking_species(u%model%lat%branch(0)%param))

    nl=nl+1; lines(nl) = 'beam_init components (set by "set beam_init ..."):'
    nl=nl+1; write(lines(nl), amt) '  %position_file          = ', quote(beam_init%position_file)
    nl=nl+1; write(lines(nl), amt) '  %distribution_type      = ', quoten(beam_init%distribution_type)
    nl=nl+1; write(lines(nl), lmt) '  %full_6D_coupling_calc  = ', beam_init%full_6d_coupling_calc
    nl=nl+1; write(lines(nl), lmt) '  %use_particle_start     = ', beam_init%use_particle_start
    nl=nl+1; write(lines(nl), lmt) '  %use_t_coords           = ', beam_init%use_t_coords
    nl=nl+1; write(lines(nl), lmt) '  %use_z_as_t             = ', beam_init%use_z_as_t
    if (beam_init%use_particle_start) then
      nl=nl+1; write(lines(nl), '(a, 6es16.8, 3x, a)') '  %center                 = ', beam_init%center, '! Note: Will use particle_start instead'
      nl=nl+1; write(lines(nl), '(a, 3es16.8, 3x, a)') '  %spin                   = ', beam_init%spin,   '! Note: Will use particle_start instead'
      nl=nl+1; write(lines(nl), '(a, 6es16.8, 3x, a)') '  particle_start[x]       = ', u%model%lat%particle_start%vec(1)
      nl=nl+1; write(lines(nl), '(a, 6es16.8, 3x, a)') '  particle_start[px]      = ', u%model%lat%particle_start%vec(2)
      nl=nl+1; write(lines(nl), '(a, 6es16.8, 3x, a)') '  particle_start[y]       = ', u%model%lat%particle_start%vec(3)
      nl=nl+1; write(lines(nl), '(a, 6es16.8, 3x, a)') '  particle_start[py]      = ', u%model%lat%particle_start%vec(4)
      nl=nl+1; write(lines(nl), '(a, 6es16.8, 3x, a)') '  particle_start[z]       = ', u%model%lat%particle_start%vec(5)
      nl=nl+1; write(lines(nl), '(a, 6es16.8, 3x, a)') '  particle_start[pz]      = ', u%model%lat%particle_start%vec(6)
    else
      nl=nl+1; write(lines(nl), rmt) '  %spin                   = ', beam_init%spin
      nl=nl+1; write(lines(nl), rmt) '  %center                 = ', beam_init%center
      nl=nl+1; write(lines(nl), rmt) '  %center_jitter          = ', beam_init%center_jitter
    endif
    nl=nl+1; write(lines(nl), rmt)  '  %t_offset               = ', beam_init%t_offset
    nl=nl+1; write(lines(nl), imt)  '  %n_particle             = ', beam_init%n_particle
    nl=nl+1; write(lines(nl), imt)  '  %ix_turn                = ', beam_init%ix_turn
    nl=nl+1; write(lines(nl), imt2) '  %n_bunch                = ', beam_init%n_bunch, '! Note: 0 => Create one bunch if not reading from a file'
    nl=nl+1; write(lines(nl), rmt)  '  %bunch_charge           = ', beam_init%bunch_charge
    if (u%model%lat%branch(0)%param%particle == photon$) then
      nl=nl+1; write(lines(nl), '(2(a, es16.8))') '  %a_emit                 = ', beam_init%a_emit
      nl=nl+1; write(lines(nl), '(2(a, es16.8))') '  %b_emit                 = ', beam_init%b_emit
    elseif (bb%beam_init_used%a_emit == real_garbage$) then
      nl=nl+1; write(lines(nl), '(2(a, es16.8))') '  %a_norm_emit            = ', beam_init%a_norm_emit
      nl=nl+1; write(lines(nl), '(2(a, es16.8))') '  %b_norm_emit            = ', beam_init%b_norm_emit
      nl=nl+1; write(lines(nl), '(2(a, es16.8))') '  %a_emit                 = ', beam_init%a_emit
      nl=nl+1; write(lines(nl), '(2(a, es16.8))') '  %b_emit                 = ', beam_init%b_emit
      nl=nl+1; write(lines(nl), '(2(a, es16.8))') '  %sig_z                  = ', beam_init%sig_z
      nl=nl+1; write(lines(nl), '(2(a, es16.8))') '  %sig_pz                 = ', beam_init%sig_pz
    else
      nl=nl+1; write(lines(nl), '(2(a, es16.8))') '  %a_norm_emit            = ', beam_init%a_norm_emit, '  ! Value used:', bb%beam_init_used%a_norm_emit
      nl=nl+1; write(lines(nl), '(2(a, es16.8))') '  %b_norm_emit            = ', beam_init%b_norm_emit, '  ! Value used:', bb%beam_init_used%b_norm_emit
      nl=nl+1; write(lines(nl), '(2(a, es16.8))') '  %a_emit                 = ', beam_init%a_emit,      '  ! Value used:', bb%beam_init_used%a_emit
      nl=nl+1; write(lines(nl), '(2(a, es16.8))') '  %b_emit                 = ', beam_init%b_emit,      '  ! Value used:', bb%beam_init_used%b_emit
      nl=nl+1; write(lines(nl), '(2(a, es16.8))') '  %sig_z                  = ', beam_init%sig_z,       '  ! Value used:', bb%beam_init_used%sig_z
      nl=nl+1; write(lines(nl), '(2(a, es16.8))') '  %sig_pz                 = ', beam_init%sig_pz,      '  ! Value used:', bb%beam_init_used%sig_pz
    endif
    nl=nl+1; write(lines(nl), rmt) '  %dPz_dz                 = ', beam_init%dPz_dz
    nl=nl+1; write(lines(nl), rmt) '  %dt_bunch               = ', beam_init%dt_bunch
    nl=nl+1; write(lines(nl), rmt) '  %emit_jitter            = ', beam_init%emit_jitter
    nl=nl+1; write(lines(nl), rmt) '  %sig_z_jitter           = ', beam_init%sig_z_jitter
    nl=nl+1; write(lines(nl), rmt) '  %sig_pz_jitter          = ', beam_init%sig_pz_jitter
    nl=nl+1; write(lines(nl), amt) '  %species                = ', quote(beam_init%species)
    nl=nl+1; write(lines(nl), lmt) '  %renorm_center          = ', beam_init%renorm_center
    nl=nl+1; write(lines(nl), lmt) '  %renorm_sigma           = ', beam_init%renorm_sigma
    nl=nl+1; write(lines(nl), amt) '  %random_engine          = ', quote(beam_init%random_engine)
    nl=nl+1; write(lines(nl), amt) '  %random_gauss_converter = ', quote(beam_init%random_gauss_converter)
    nl=nl+1; write(lines(nl), f3mt)'  %random_sigma_cutoff    = ', beam_init%random_sigma_cutoff
    fmt = '(a, i1, a, es16.8)'
    do i = 1, 3
      if (upcase(beam_init%distribution_type(i)) == 'ELLIPSE') then
        nl=nl+1; write(lines(nl), iimt) '  %ellipse(', i, ')%part_per_ellipse  = ', beam_init%ellipse(i)%part_per_ellipse
        nl=nl+1; write(lines(nl), iimt) '  %ellipse(', i, ')%n_ellipse         = ', beam_init%ellipse(i)%n_ellipse
        nl=nl+1; write(lines(nl), irmt) '  %ellipse(', i, ')%sigma_cutoff      = ', beam_init%ellipse(i)%sigma_cutoff
      elseif (upcase(beam_init%distribution_type(i)) == 'GRID') then
        nl=nl+1; write(lines(nl), iimt) '  %grid(', i, ')%n_x            = ', beam_init%grid(i)%n_x
        nl=nl+1; write(lines(nl), iimt) '  %grid(', i, ')%n_px           = ', beam_init%grid(i)%n_px
        nl=nl+1; write(lines(nl), irmt) '  %grid(', i, ')%x_min          = ', beam_init%grid(i)%x_min
        nl=nl+1; write(lines(nl), irmt) '  %grid(', i, ')%x_max          = ', beam_init%grid(i)%x_max
        nl=nl+1; write(lines(nl), irmt) '  %grid(', i, ')%px_min         = ', beam_init%grid(i)%px_min
        nl=nl+1; write(lines(nl), irmt) '  %grid(', i, ')%px_max         = ', beam_init%grid(i)%px_max
      endif
    enddo
    if (any(upcase(beam_init%distribution_type) == 'KV')) then
      nl=nl+1; write(lines(nl), imt) '  %kv%part_per_phi(1:2) = ', beam_init%kv%part_per_phi
      nl=nl+1; write(lines(nl), imt) '  %kv%n_i2              = ', beam_init%kv%n_i2
      nl=nl+1; write(lines(nl), rmt) '  %kv%a                 = ', beam_init%kv%a
    endif

    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'bmad_com components (set by "set bmad_com ..."):'
    nl=nl+1; write(lines(nl), lmt) '  %sr_wakes_on                     = ', bmad_com%sr_wakes_on
    nl=nl+1; write(lines(nl), lmt) '  %lr_wakes_on                     = ', bmad_com%lr_wakes_on
    nl=nl+1; write(lines(nl), lmt) '  %csr_and_space_charge_on         = ', bmad_com%csr_and_space_charge_on
    nl=nl+1; write(lines(nl), lmt) '  %spin_tracking_on                = ', bmad_com%spin_tracking_on
    nl=nl+1; write(lines(nl), lmt) '  %spin_sokolov_ternov_flipping_on = ', bmad_com%spin_sokolov_ternov_flipping_on
    nl=nl+1; write(lines(nl), lmt) '  %radiation_damping_on            = ', bmad_com%radiation_damping_on
    nl=nl+1; write(lines(nl), lmt) '  %radiation_zero_average          = ', bmad_com%radiation_zero_average
    nl=nl+1; write(lines(nl), lmt) '  %radiation_fluctuations_on       = ', bmad_com%radiation_fluctuations_on

    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'space_charge_com components (set by "set space_charge_com ..."):'
    nl=nl+1; write(lines(nl), rmt) '  %ds_track_step                  = ', space_charge_com%ds_track_step
    nl=nl+1; write(lines(nl), rmt) '  %dt_track_step                  = ', space_charge_com%dt_track_step
    nl=nl+1; write(lines(nl), rmt) '  %cathode_strength_cutoff        = ', space_charge_com%cathode_strength_cutoff
    nl=nl+1; write(lines(nl), rmt) '  %rel_tol_tracking               = ', space_charge_com%rel_tol_tracking
    nl=nl+1; write(lines(nl), rmt) '  %abs_tol_tracking               = ', space_charge_com%abs_tol_tracking
    nl=nl+1; write(lines(nl), rmt) '  %beam_chamber_height            = ', space_charge_com%beam_chamber_height
    nl=nl+1; write(lines(nl), rmt) '  %lsc_sigma_cutoff               = ', space_charge_com%lsc_sigma_cutoff
    nl=nl+1; write(lines(nl), rmt) '  %particle_sigma_cutoff          = ', space_charge_com%particle_sigma_cutoff
    nl=nl+1; write(lines(nl), imt) '  %space_charge_mesh_size         = ', space_charge_com%space_charge_mesh_size
    nl=nl+1; write(lines(nl), imt) '  %csr3d_mesh_size                = ', space_charge_com%csr3d_mesh_size
    nl=nl+1; write(lines(nl), imt) '  %n_bin                          = ', space_charge_com%n_bin
    nl=nl+1; write(lines(nl), imt) '  %particle_bin_span              = ', space_charge_com%particle_bin_span
    nl=nl+1; write(lines(nl), imt) '  %n_shield_images                = ', space_charge_com%n_shield_images
    nl=nl+1; write(lines(nl), imt) '  %sc_min_in_bin                  = ', space_charge_com%sc_min_in_bin
    nl=nl+1; write(lines(nl), lmt) '  %lsc_kick_transverse_dependence = ', space_charge_com%lsc_kick_transverse_dependence
    nl=nl+1; write(lines(nl), amt) '  %diagnostic_output_file         = ', quote(space_charge_com%diagnostic_output_file)

    if (size(lat%branch) > 1) then
      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = 'Branch  Start_Track                            End_Track'
      do i = 0, ubound(lat%branch, 1)
        branch => lat%branch(i)
        bb => u%model_branch(i)%beam

        if (bb%ix_track_start == not_set$) then
          nl=nl+1; write (lines(nl), '(i5, 6x, a)') i, '.... No Tracking Done ...'
        else
          nl=nl+1; write (lines(nl), '(i5, 4x, a, t48, a)') i, trim(bb%track_start), trim(bb%track_end)
        endif
      enddo
    endif

  !-------------------------------
  ! have element index

  else
    if (what_to_show == '-comb') then
      if (.not. allocated(tao_branch%bunch_params_comb)) then
        nl=nl+1; lines(nl) = 'Beam parameter comb not calculated (check comb_ds_step)' 
        return
      endif

      n = tao_branch%bunch_params_comb(1)%n_pt

      if (n < 0 .and. tao_branch%bunch_params_comb(1)%ds_save < 0) then
        nl=nl+1; lines(nl) = 'Comb not calculated since comb_ds_save not set positive.'
        return
      endif

      if (ix_comb < 0 .or. ix_comb > n) then
        nl=nl+1; lines(nl) = 'Comb index out of range: [0, ' // int_str(n) // ']'
        return
      endif
      bunch_p => tao_branch%bunch_params_comb(1)%pt(ix_comb)
      nl=nl+1; lines(nl) = 'Bunch parameters at comb index: ' // int_str(ix_comb)

    else
      if (zb(1) >= 0) then
        beam => u%model_branch(0)%ele(ele%ix_ele)%beam
        if (.not. allocated(beam%bunch)) then
          nl=nl+1; lines(nl) = 'Bunch parameters at: ' // trim(ele%name) // ' ' //  ele_loc_name(ele, .false., '()')
          return
        endif
        call calc_bunch_params_z_slice(beam%bunch(s%global%bunch_to_plot), bunch_params, zb, err, .true., ele = ele)
        if (err) return
        bunch_p => bunch_params
      else
        bunch_p => tao_branch%bunch_params(ele%ix_ele)
      endif

      nl=nl+1; lines(nl) = 'Bunch parameters at: ' // trim(ele%name) // ' ' //  ele_loc_name(ele, .false., '()')
      found = (allocated(u%model_branch(ele%ix_branch)%ele(ele%ix_ele)%beam%bunch))
    endif

    n_live = bunch_p%n_particle_live
    n_tot = bunch_p%n_particle_tot
    n = s%global%bunch_to_plot

    if (n_tot == 0) then
      nl=nl+1; lines(nl) = 'Beam not tracked through this element!'
      result_id = 'beam:no_particles'
      return
    endif

    if (n_live == 0) then
      nl=nl+1; lines(nl) = 'No live particles!'
      result_id = 'beam:no_live'
      return
    endif

    nl=nl+1; write(lines(nl), imt)  '  Parameters for bunch:       ', n
    nl=nl+1; write(lines(nl), rmt)  '  S-position:                 ', bunch_p%s
    nl=nl+1; write(lines(nl), imt)  '  In branch:                  ', ix_branch
    nl=nl+1; write(lines(nl), imt)  '  Particles surviving:        ', n_live
    nl=nl+1; write(lines(nl), imt)  '  Particles lost:             ', n_tot - n_live
    nl=nl+1; write(lines(nl), f3mt) '  Particles lost (%):         ', 100 * real(n_tot - n_live) / n_tot
    nl=nl+1; write(lines(nl), iimt) '  Good steps:                 ', bunch_p%n_good_steps, '    ! Set when tracking with space charge.'
    nl=nl+1; write(lines(nl), iimt) '  Bad steps:                  ', bunch_p%n_bad_steps,  '    ! Set when tracking with space charge.'

    if (branch%param%particle == photon$) then
      nl=nl+1; write(lines(nl), rmt)  '  Intensity:                  ', &
                        bunch_p%centroid%field(1)**2 + bunch_p%centroid%field(2)**2
    else
      nl=nl+1; write(lines(nl), rmt) '  Charge live (C):            ', bunch_p%charge_live
    endif

    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), '(24x,a)') 'x               px               y              py                z               pz          time'
    nl=nl+1; write(lines(nl), rmt) '  Centroid:    ', bunch_p%centroid%vec, bunch_p%t
    nl=nl+1; write(lines(nl), rmt) '  RMS:         ', sqrt(bunch_p%sigma(1,1)), sqrt(bunch_p%sigma(2,2)), sqrt(bunch_p%sigma(3,3)), &
                                                      sqrt(bunch_p%sigma(4,4)), sqrt(bunch_p%sigma(5,5)), sqrt(bunch_p%sigma(6,6)), bunch_p%sigma_t
    nl=nl+1; write(lines(nl), rmt) '  Max-Centroid:', bunch_p%rel_max
    nl=nl+1; write(lines(nl), rmt) '  Min-Centroid:', bunch_p%rel_min

    if (u%model%lat%branch(ix_branch)%param%particle /= photon$) then
      nl=nl+1; lines(nl) = ''
      nl=nl+1; write(lines(nl), rmt) '               norm_emitt            emit            beta           alpha'
      nl=nl+1; write(lines(nl), rmt) '  a:       ', bunch_p%a%norm_emit, bunch_p%a%emit, bunch_p%a%beta, bunch_p%a%alpha
      nl=nl+1; write(lines(nl), rmt) '  b:       ', bunch_p%b%norm_emit, bunch_p%b%emit, bunch_p%b%beta, bunch_p%b%alpha
      nl=nl+1; write(lines(nl), rmt) '  x:       ', bunch_p%x%norm_emit, bunch_p%x%emit, bunch_p%x%beta, bunch_p%x%alpha
      nl=nl+1; write(lines(nl), rmt) '  y:       ', bunch_p%y%norm_emit, bunch_p%y%emit, bunch_p%y%beta, bunch_p%y%alpha
      nl=nl+1; write(lines(nl), rmt) '  z:       ', bunch_p%z%norm_emit, bunch_p%z%emit, bunch_p%z%beta
      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = 'Sigma Mat       x              px               y              py              z             pz'
      do i = 1, 6
        nl=nl+1; write (lines(nl), '(a2, 2x, 6es16.8)') coord_name(i), bunch_p%sigma(i,:)
      enddo

      if (abs(bunch_p%a%norm_emit - bunch_p%b%norm_emit) < 0.2 * (bunch_p%a%norm_emit + bunch_p%b%norm_emit)) then
        nl=nl+1; lines(nl) = 'Note! The calculation of the emittances from the beam sigma matrix is singular when the emittances are equal!'
        nl=nl+1; lines(nl) = 'This is due to modes being degenerate. Emittance and Twiss values are untrustworthy in this limit.'
      endif
    endif

    if (what_to_show == '') then
      nl=nl+1; lines(nl) = ''
      if (found) then
        nl=nl+1; lines(nl) = 'Note: Individual particle positions are saved at this element.'
      else
        nl=nl+1; lines(nl) = 'Note: Individual particle positions are not saved at this element.'
      endif
    endif
  
  endif

!----------------------------------------------------------------------
! constraints

case ('branch')

  do 
    call tao_next_switch (what2, ['-universe'], .false., switch, err);  if (err) return
    if (switch == '') exit

    select case (switch)

    case ('-universe')
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios) ix
      u => tao_pointer_to_universe(ix)
      if (aname == '' .or. ios /= 0 .or. .not. associated(u)) then
        nl=1; lines(1) = 'CANNOT READ OR OUT-OF RANGE "-universe" ARGUMENT'
        return
      endif

    case default
      call out_io (s_error$, r_name, 'EXTRA STUFF ON LINE: ' // what2)
      return

    end select
  enddo

  lat => u%model%lat

  if (size(s%u) > 1) then
    nl=nl+1; write(lines(nl), '(a, i0)') 'For the lattice of universe: ', ix_u
  endif

  nl=nl+1; lines(nl) = '                          N_ele  N_ele   Reference      Default_                      Live'  
  nl=nl+1; lines(nl) = '  Branch                  Track    Max   Particle       Tracking_Species    Geometry  Branch  From_Fork'


  fmt = '((i3, 2a), t26, i6, i7, t42, a, t57, a, t77, a, t87, l2, 6x, a)'
  do i = 0, ubound(lat%branch, 1)
    branch => lat%branch(i)
    ele_name = ''
    if (branch%ix_from_ele > 0) write (ele_name, '(i0, a, i0)') branch%ix_from_branch, '>>', branch%ix_from_ele

    nl=nl+1; write(lines(nl), fmt) i, ': ', branch%name, branch%n_ele_track, branch%n_ele_max, &
              trim(species_name(branch%param%particle)), trim(species_name(branch%param%default_tracking_species)), &
              trim(geometry_name(branch%param%geometry)), branch%param%live_branch, ele_name
  enddo

  nl=nl+1; lines(nl) = ''
  nl=nl+1; lines(nl) = '                                                                                                      Defines'
  nl=nl+1; lines(nl) = '  Fork_Element                              Forking_To                                   Direction    To_Branch?'
  nl0 = nl

  fmt = '((i3, a, i0, 4a), t45, (2(i0, a), 3a), t90, i2, l14)'
  do i = 0, ubound(lat%branch, 1)
    branch => lat%branch(i)
    do j = 1, branch%n_ele_max
      ele => branch%ele(j)
      if (ele%key /= fork$ .and. ele%key /= photon_fork$) cycle
      branch2 => lat%branch(nint(ele%value(ix_to_branch$)))
      logic = (ele%ix_ele == branch2%ix_from_ele .and. ele%ix_branch == branch2%ix_from_branch)
      ele2 => branch2%ele(nint(ele%value(ix_to_element$)))
      nl=nl+1; write(lines(nl), fmt) i, '>>', j, ': ', trim(branch%name), '>>', trim(ele%name), &
                branch2%ix_branch, '>>', ele2%ix_ele, ': ', trim(branch2%name), '>>', trim(ele2%name), &
                nint(ele%value(direction$)), logic
    enddo
  enddo

  if (nl == nl0) nl = nl - 3 ! Erase header if no info.

!----------------------------------------------------------------------
! building_wall

case ('building_wall')

  if (size(s%building_wall%section) == 0) then
    nl=nl+1; lines(nl) = 'No building wall defined.'
    result_id = 'building_wall:none'
    return
  endif

  nl=nl+1; lines(nl) = 'Orientation:'
  nl=nl+1; write (lines(nl), '(a, f10.4)') '  theta    =', s%building_wall%orientation%theta
  nl=nl+1; write (lines(nl), '(a, f10.4)') '  x_offset =', s%building_wall%orientation%x_offset
  nl=nl+1; write (lines(nl), '(a, f10.4)') '  z_offset =', s%building_wall%orientation%z_offset

  do i = 1, size(s%building_wall%section)
    section => s%building_wall%section(i)
    n = nl + size(section%point)
    if (n + 10 > size(lines)) call re_allocate (lines, n, .false.)
    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), '(a, i0, 4a)') 'Section(', i, ')    Name: ', quote(section%name), '  Constraint: ', trim(section%constraint)

    nl=nl+1; lines(nl) = '                        Z           X      Radius    Z_center    X_center'
    do j = 1, size(section%point)
      pt => section%point(j)
      if (pt%radius == 0) then
        nl=nl+1; write(lines(nl), '(a, i3, a, 3f12.3)') &
              '  point(', j, '):', pt%z, pt%x, pt%radius
      else
        nl=nl+1; write(lines(nl), '(a, i3, a, 5f12.3)') &
              '  point(', j, '):', pt%z, pt%x, pt%radius, pt%z_center, pt%x_center
      endif
    enddo
  enddo

!----------------------------------------------------------------------
! chromaticity

case ('chromaticity')

  what_to_show = ''

  do 
    call tao_next_switch (what2, [character(16):: '-universe', '-taylor'], .false., switch, err);  if (err) return
    if (switch == '') exit

    select case (switch)

    case ('-universe')
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios) ix
      u => tao_pointer_to_universe(ix)
      if (aname == '' .or. ios /= 0 .or. .not. associated(u)) then
        nl=1; lines(1) = 'CANNOT READ OR OUT-OF RANGE "-universe" ARGUMENT'
        return
      endif

    case ('-taylor')
      what_to_show = switch

    case default
      call out_io (s_error$, r_name, 'EXTRA STUFF ON LINE: ' // what2)
      return

    end select
  enddo

  !

  if (branch%param%geometry /= closed$) then
    nl=1;    lines(1) = 'BRANCH GEOMETRY NOT CLOSED'
    nl=nl+1; lines(2) = 'USE THE "set branch ' // int_str(branch%ix_branch) // ' geometry = closed" TO CHANGE THIS.'
    return
  endif

  tao_lat => tao_pointer_to_tao_lat (u, model$)
  if (.not. u%calc%one_turn_map) call tao_ptc_normal_form (.true., tao_lat, ix_branch, rf_on = no$)

  bmad_nf => tao_branch%bmad_normal_form
  ptc_nf  => tao_branch%ptc_normal_form

  nl=nl+1; lines(nl) = '  Note: Calculation is done with RF off.'
  nl=nl+1; lines(nl) = '  N     chrom_ptc.a.N     chrom_ptc.b.N'

  do i = 0, ptc_private%taylor_order_ptc-1
    expo = [0, 0, 0, 0, 0, i]
    z1 =  real(ptc_nf%phase(1) .sub. expo)
    z2 =  real(ptc_nf%phase(2) .sub. expo)
    if (i == 0) then
      nl=nl+1; write (lines(nl), '(i3, 2es18.7, a)') i, z1, z2, '  ! 0th order are the tunes'
    else
      nl=nl+1; write (lines(nl), '(i3, 2es18.7)') i, z1, z2
    endif
  enddo

  nl=nl+1; lines(nl) = ''
  nl=nl+1; lines(nl) = '  N   slip_factor_ptc.N   momentum_compaction_ptc.N'

  do i = 1, ptc_private%taylor_order_ptc
    expo = [0, 0, 0, 0, 0, i]
    z1 = -real(ptc_nf%phase(3) .sub. expo) / branch%param%total_length
    z2 =  real(ptc_nf%path_length .sub. expo) / branch%param%total_length
    nl=nl+1; write (lines(nl), '(i3, 2x, 2es18.7)') i, z1, z2
  enddo

  !

  if (what_to_show == '-taylor') then
    call alloc(ptc_ctaylor)

    do i = 1, 3
        ! ptc_ctaylor = ptc_nf%spin_tune * ci_phasor() * ptc_nf%normal_form%atot**(-1)
      ptc_ctaylor = ptc_nf%phase(i) * ci_phasor() * ptc_nf%normal_form%atot**(-1)
      bmad_ctaylor = ptc_ctaylor
      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = 'Taylor series: ' // trim(mode(i)) // ' tune'
      call type_complex_taylors ([bmad_ctaylor], lines = alloc_lines, n_lines = n, out_type = 'NONE')
      if (size(lines) < nl+n+100) call re_allocate (lines, nl+n+100, .false.)
      lines(nl+1:nl+n) = alloc_lines(1:n)
      nl = nl + n

      !if (i == 4 .and. .not. associated(bmad_ctaylor%term) .and. .not. bmad_com%spin_tracking_on) then
      !  nl=nl+1; lines(nl) = 'Spin tracking is turned on with: "set bmad_com spin_tracking_on = T".'
      !endif
    enddo

    call kill(ptc_ctaylor)
  endif

!----------------------------------------------------------------------
! constraints

case ('constraints')

  call tao_show_constraints (0, 'ALL')
  call tao_show_constraints (0, 'TOP10')

!----------------------------------------------------------------------
! control

case ('control')

  ele_name = what2

  call tao_pick_universe (ele_name, ele_name, picked_uni, err, ix_u)
  if (err) return
  u => s%u(ix_u)
  call tao_locate_elements (ele_name, ix_u, eles, err, multiple_eles_is_err = .true.)
  if (err .or. size(eles) == 0) return

  ele => eles(1)%ele
  call tao_control_tree_list(ele, eles)

  do i = size(eles), 1, -1  ! Show lords first
    ele => eles(i)%ele
    call re_allocate(lines, nl+10+ele%n_lord+ele%n_slave+ele%n_slave_field)

    do j = 1, ele%n_lord
      lord => pointer_to_lord(ele, j)
      nl=nl+1; write (lines(nl), '(5a)') 'Lord: ', ele_loc_name(lord, .true.), lord%name, key_name(lord%key), control_name(lord%lord_status)
    enddo

    nl=nl+1; write (lines(nl), '(5a)') 'Element: ', ele_loc_name(ele, .true.), ele%name, key_name(ele%key), control_name(ele%lord_status)

    do j = 1, ele%n_slave+ele%n_slave_field
      slave => pointer_to_slave(ele, j)
      nl=nl+1; write (lines(nl), '(5a)') 'Slave: ', ele_loc_name(slave, .true.), slave%name, key_name(slave%key), control_name(slave%slave_status)
    enddo

    if (i /= 1) then
      nl=nl+1; lines(nl) = ''
    endif
  enddo

!----------------------------------------------------------------------
! curve

case ('curve')

  ! Look for switches

  show_sym = .false.
  show_line = .false.
  print_header = .true.
  attrib0 = ''

  do
    call tao_next_switch (what2, [character(20):: '-symbol', '-line', '-no_header'], .true., switch, err)
    if (err) return
    select case (switch)
    case ('');           exit
    case ('-symbol');    show_sym = .true.
    case ('-line');      show_line = .true.
    case ('-no_header'); print_header = .false.
    case default
      if (attrib0 /= '') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON LINE: ' // attrib0)
        return
      endif
      attrib0 = switch
    end select
  enddo

  ! Find particular plot

  call tao_find_plots (err, attrib0, 'BOTH', curve = curve, blank_means_all = .true., only_visible = .false.)
  if (err .or. size(curve) == 0) then
    nl=1; lines(1) = 'THIS IS NOT A CURVE NAME'
    return
  endif

  ! print info on particular plot, graph, or curve

  if (size(curve) > 0 .and. .not. show_line .and. .not. show_sym) then
    c1 => curve(1)%c

    if (print_header) then
      if (associated(c1%g%p%r)) then
        nl=nl+1; lines(nl) = '                                   |          Line           |          Symbols'
        nl=nl+1; lines(nl) = 'Curves:                            |     Mean         RMS    |      Mean         RMS'
        do i = 1, size(curve)
          cc => curve(i)%c
          nl=nl+1; lines(nl) = '   ' // trim(tao_curve_name(cc, .true.))
          if (allocated(cc%x_line)) then
            call tao_curve_rms_calc(cc, 'LINE', rms, mean)
            write(lines(nl)(36:), '(2es12.4)') mean, rms
          endif
          if (allocated(cc%x_symb)) then
            call tao_curve_rms_calc(cc, 'SYMBOL', rms, mean)
            write(lines(nl)(63:), '(2es12.4)') mean, rms
          endif
        enddo

      else
        nl=nl+1; lines(nl) = 'Curves:'
        do i = 1, size(curve)
          nl=nl+1; lines(nl) = '   ' // trim(tao_curve_name(curve(i)%c))
        enddo
      endif

      ele => tao_curve_ele_ref(c1, .false.)

      nl=nl+1; write(lines(nl), amt)  'data_source          = ', quote(c1%data_source)
      nl=nl+1; write(lines(nl), amt)  'data_index           = ', quote(c1%data_index)
      nl=nl+1; write(lines(nl), amt)  'data_type_x          = ', quote(c1%data_type_x)
      nl=nl+1; write(lines(nl), amt)  'data_type            = ', quote(c1%data_type)
      nl=nl+1; write(lines(nl), amt)  'legend_text          = ', quote(c1%legend_text)
      if (associated(ele)) then
        nl=nl+1; write(lines(nl), amt)  'ele_ref_name         = ', quote(c1%ele_ref_name), ele_full_name(ele, '  (&#)')
      else
        nl=nl+1; write(lines(nl), amt)  'ele_ref_name         = ', quote(c1%ele_ref_name), '  (--)'
      endif
      nl=nl+1; write(lines(nl), amt)  'component            = ', quote(c1%component)
      nl=nl+1; write(lines(nl), imt)  'ix_branch            = ', c1%ix_branch
      nl=nl+1; write(lines(nl), imt)  'ix_bunch             = ', c1%ix_bunch
      nl=nl+1; write(lines(nl), imt)  'ix_universe          = ', c1%ix_universe
      nl=nl+1; write(lines(nl), imt)  'n_turn               = ', c1%n_turn
      nl=nl+1; write(lines(nl), imt)  'symbol_every         = ', c1%symbol_every
      nl=nl+1; write(lines(nl), rmt)  'y_axis_scale_factor  = ', c1%y_axis_scale_factor
      nl=nl+1; write(lines(nl), lmt)  'use_y2               = ', c1%use_y2
      nl=nl+1; write(lines(nl), lmt)  'draw_line            = ', c1%draw_line
      nl=nl+1; write(lines(nl), lmt)  'draw_symbols         = ', c1%draw_symbols
      nl=nl+1; write(lines(nl), lmt)  'draw_symbol_index    = ', c1%draw_symbol_index
      nl=nl+1; write(lines(nl), lmt)  'draw_error_bars      = ', c1%draw_error_bars
      nl=nl+1; write(lines(nl), lmt)  'smooth_line_calc     = ', c1%smooth_line_calc
      nl=nl+1; write(lines(nl), lmt)  'valid                = ', c1%valid
      nl=nl+1; write(lines(nl), iamt) 'line%width           = ', c1%line%width
      nl=nl+1; write(lines(nl), amt)  'line%color           = ', c1%line%color
      nl=nl+1; write(lines(nl), amt)  'line%pattern         = ', c1%line%pattern
      nl=nl+1; write(lines(nl), amt)  'symbol%type          = ', c1%symbol%type
      nl=nl+1; write(lines(nl), f3mt) 'symbol%height        = ', c1%symbol%height
      nl=nl+1; write(lines(nl), amt)  'symbol%fill_pattern  = ', c1%symbol%fill_pattern
      nl=nl+1; write(lines(nl), imt)  'symbol%line_width    = ', c1%symbol%line_width
      nl=nl+1; write(lines(nl), amt)  'symbol%color         = ', c1%symbol%color
      nl=nl+1; write(lines(nl), amt)  'z_color%data_type    = ', quote(c1%z_color%data_type)
      nl=nl+1; write(lines(nl), lmt)  'z_color%is_on        = ', c1%z_color%is_on
      nl=nl+1; write(lines(nl), lmt)  'z_color%autoscale    = ', c1%z_color%autoscale
      nl=nl+1; write(lines(nl), rmt)  'z_color_min          = ', c1%z_color%min
      nl=nl+1; write(lines(nl), rmt)  'z_color_max          = ', c1%z_color%max
      
      ! Histogram specific parameters
      if (c1%g%type == 'histogram') then
        nl=nl+1; write(lines(nl), lmt)  'hist%density_normalized = ', c1%hist%density_normalized 
        nl=nl+1; write(lines(nl), lmt)  'hist%weight_by_charge   = ', c1%hist%weight_by_charge
        nl=nl+1; write(lines(nl), rmt)  'hist%minimum            = ', c1%hist%minimum
        nl=nl+1; write(lines(nl), rmt)  'hist%maximum            = ', c1%hist%maximum
        nl=nl+1; write(lines(nl), rmt)  'hist%width              = ', c1%hist%width
        nl=nl+1; write(lines(nl), rmt)  'hist%center             = ', c1%hist%center
        nl=nl+1; write(lines(nl), imt)  'hist%number             = ', c1%hist%number
      endif

      if (.not. c1%valid) then
        nl=nl+1; write(lines(nl), amt)  'Curve is not valid since: ', c1%why_invalid
      endif
    endif

  ! Field at constant offset parameters
  if (substr(c1%data_type, 1, 3) == 'b0_' .or. substr(c1%data_type, 1, 3) == 'e0_') then
    nl=nl+1; write(lines(nl), rmt)      'orbit%x                 = ', c1%orbit%x
    nl=nl+1; write(lines(nl), rmt)      'orbit%y                 = ', c1%orbit%y
    nl=nl+1; write(lines(nl), rmt)      'orbit%t                 = ', c1%orbit%t
  endif


    return
  endif

  ! Show symbol points
  
  if (show_sym) then
    nc = 0
    do j = 1, size(curve)
      if (.not. allocated(curve(j)%c%x_symb)) cycle
      nc = max(nc, size(curve(j)%c%x_symb))
    enddo

    if (print_header) then
      nl=nl+1; lines(nl)   = ''
      nl=nl+1; lines(nl)   = '# Symbol points:'
      nl=nl+1; lines(nl)   = '#     i  ix_sym    x-axis '
      do j = 1, size(curve)
        str = curve(j)%c%name
        lines(nl) = lines(nl)(1:28+(j-1)*14) // adjustr(str(1:14))
      enddo

      if (nc == 0) then
        nl=nl+1; lines(nl) = '#     No Symbol Points'
      endif
    endif

    n = size(curve)
    allocate (ix_c(n), value(n), valid(n))
    ix_c = 1
    
    id = 0
    do
      value_min = 1e30
      valid = .false.
      do i = 1, n
        if (ix_c(i) > size(curve(i)%c%x_symb)) cycle
        value(i) = curve(i)%c%x_symb(ix_c(i))
        valid(i) = .true.
        value_min = min(value_min, value(i))
      enddo

      if (all(.not. valid)) exit

      ix_min = 100000
      do i = 1, n
        if (.not. valid(i) .or. value(i) /= value_min) cycle
        ix_min = min(ix_min, curve(i)%c%ix_symb(ix_c(i)))
      enddo

      id = id + 1
      nl=nl+1; write (lines(nl), '(2i7, 10es14.6)') id, ix_min, value_min
      do i = 1, n
        if (valid(i) .and. value(i) == value_min .and. curve(i)%c%ix_symb(ix_c(i)) == ix_min) then
          write (lines(nl)(29+(i-1)*14:), '(es14.6)') curve(i)%c%y_symb(ix_c(i))
          ix_c(i) = ix_c(i) + 1
        endif
      enddo

      if (nl+10 > size(lines)) call re_allocate(lines, nl+100, .false.)
    enddo
  endif

  ! Show line points

  if (show_line) then

    nc = 0
    ix0 = 0
    aligned = .true.    ! True => can have one x column for all curves.
    do j = 1, size(curve)
      if (.not. allocated(curve(j)%c%x_line)) cycle
      ! Physical aperture curves (part of dynamic_aperture plots) are not interesting and generate too many lines.
      if (size(curve) > 1 .and. curve(j)%c%data_type == 'physical_aperture') cycle
      nc = max(nc, size(curve(j)%c%x_line))
      if (ix0 == 0) ix0 = j
      if (size(curve(j)%c%x_line) /= size(curve(ix0)%c%x_line)) then
        aligned = .false.
      else
        if (any(curve(j)%c%x_line /= curve(ix0)%c%x_line)) aligned = .false.
      endif
    enddo

    n = nl + nc + 10
    if (n > size(lines)) call re_allocate(lines, n, .false.)

    if (print_header) then
      nl=nl+1; lines(nl)   = ''
      nl=nl+1; lines(nl)   = '# Smooth line points:'
      if (aligned) then
        nl=nl+1; lines(nl)   = '# index        x-axis'
      else
        nl=nl+1; lines(nl)   = '# index'
      endif

      nl0 = nl

      if (nc == 0) then
        nl=nl+1; lines(nl) = '#     No Line Points'
      endif

      k = 0
      do j = 1, size(curve)
        if (.not. allocated(curve(j)%c%x_line)) cycle
        if (size(curve) > 1 .and. curve(j)%c%data_type == 'physical_aperture') cycle
        k = k + 1
        str = curve(j)%c%name
        if (aligned) then
          lines(nl0) = lines(nl0)(1:21+(k-1)*14) // adjustr(str(1:14))
        else
          lines(nl0) = lines(nl0)(1:7+(k-1)*28) // '        x-axis' // adjustr(str(1:14))
        endif
      enddo
    endif

    do i = 1, nc
      if (aligned) then
        nl=nl+1; write(lines(nl), '(i7, es14.6)') i, curve(ix0)%c%x_line(i)
      else
        nl=nl+1; write(lines(nl), '(i7, es14.6)') i
      endif

      k = 0
      do j = 1, size(curve)
        if (.not. allocated(curve(j)%c%x_line)) cycle
        if (size(curve(j)%c%x_line) < i) cycle
        ! Physical aperture curves (part of dynamic_aperture plots) are not interesting and generate too many lines.
        if (size(curve) > 1 .and. curve(j)%c%data_type == 'physical_aperture') cycle
        k = k + 1
        if (aligned) then
          write(lines(nl)(22+(k-1)*14:), '(10es14.6)') curve(j)%c%y_line(i)
        else
          write(lines(nl)(8+(k-1)*28:), '(10es14.6)') curve(j)%c%x_line(i), curve(j)%c%y_line(i)
        endif
      enddo

    enddo
  endif

!----------------------------------------------------------------------
! data

case ('data')

  ! If just "show data" then show all names

  call tao_pick_universe (word1, line1, picked_uni, err)
  if (err) return

  ! get pointers to the data

  if (word1 == '') word1 = '*@*'
  call tao_find_data (err, word1, d2_array, d1_array, d_array)
  if (err) return

  ! If d_ptr points to something then show the datum info.

  if (size(d_array) == 1 .and. word1 /= '*@*') then
    d_ptr => d_array(1)%d
    nl=nl+1; lines(nl) = ''
    if (size(s%u) > 1) then
      nl=nl+1; write(lines(nl), '(2(a, i0))') 'Universe: ', d_ptr%d1%d2%ix_universe, '  of: ', ubound(s%u, 1)
    endif
    nl=nl+1; write(lines(nl), amt)    '%ele_name          = ', quote(d_ptr%ele_name)
    nl=nl+1; write(lines(nl), amt)    '%ele_start_name    = ', quote(d_ptr%ele_start_name)
    nl=nl+1; write(lines(nl), amt)    '%ele_ref_name      = ', quote(d_ptr%ele_ref_name)
    nl=nl+1; write(lines(nl), amt)    '%data_type         = ', quote(d_ptr%data_type)
    nl=nl+1; write(lines(nl), amt)    '%data_source       = ', quote(d_ptr%data_source)
    if (d_ptr%id /= '') then
      nl=nl+1; write(lines(nl), amt)  '%id                = ', quote(d_ptr%id)
    endif
    nl=nl+1; write(lines(nl), imt)    '%ix_branch         = ', d_ptr%ix_branch
    nl=nl+1; write(lines(nl), imt)    '%ix_ele            = ', d_ptr%ix_ele
    nl=nl+1; write(lines(nl), imt)    '%ix_ele_start      = ', d_ptr%ix_ele_start
    nl=nl+1; write(lines(nl), imt)    '%ix_ele_ref        = ', d_ptr%ix_ele_ref
    nl=nl+1; write(lines(nl), imt)    '%ix_ele_merit      = ', d_ptr%ix_ele_merit
    nl=nl+1; write(lines(nl), imt)    '%ix_dmodel         = ', d_ptr%ix_dModel
    nl=nl+1; write(lines(nl), imt)    '%ix_d1             = ', d_ptr%ix_d1
    nl=nl+1; write(lines(nl), imt)    '%ix_data           = ', d_ptr%ix_data
    nl=nl+1; write(lines(nl), imt)    '%ix_bunch          = ', d_ptr%ix_bunch
    nl=nl+1; write(lines(nl), rmt)    '%model             = ', d_ptr%model_value
    nl=nl+1; write(lines(nl), rmt)    '%design            = ', d_ptr%design_value
    nl=nl+1; write(lines(nl), rmt)    '%meas              = ', d_ptr%meas_value
    nl=nl+1; write(lines(nl), rmt)    '%ref               = ', d_ptr%ref_value
    nl=nl+1; write(lines(nl), rmt)    '%base              = ', d_ptr%base_value
    nl=nl+1; write(lines(nl), rmt)    '%error_rms         = ', d_ptr%error_rms
    nl=nl+1; write(lines(nl), rmt)    '%old               = ', d_ptr%old_value   
    nl=nl+1; write(lines(nl), rmt)    '%invalid           = ', d_ptr%invalid_value
    nl=nl+1; write(lines(nl), amt)    '%eval_point        = ', anchor_pt_name(d_ptr%eval_point)
    nl=nl+1; write(lines(nl), rmt)    '%s_offset          = ', d_ptr%s_offset
    if (d_ptr%s == real_garbage$) then  ! Happens with expressions, etc.
      nl=nl+1; write(lines(nl), rmt)    '%s                 = UNDEFINED S-POSITION'
    else
      nl=nl+1; write(lines(nl), rmt)    '%s                 = ', d_ptr%s
    endif
    nl=nl+1; write(lines(nl), amt)    '%merit_type        = ', quote(d_ptr%merit_type)
    nl=nl+1; write(lines(nl), rmt)    '%merit             = ', d_ptr%merit
    nl=nl+1; write(lines(nl), rmt)    '%delta_merit       = ', d_ptr%delta_merit
    nl=nl+1; write(lines(nl), rmt)    '%weight            = ', d_ptr%weight
    if (substr(d_ptr%data_type,1,14) == 'spin_g_matrix.') then
      nl=nl+1; write(lines(nl), rmt)  '%spin_axis%l       = ', d_ptr%spin_map%axis_input%l
      nl=nl+1; write(lines(nl), rmt)  '%spin_axis%n0      = ', d_ptr%spin_map%axis_input%n0
      nl=nl+1; write(lines(nl), rmt)  '%spin_axis%m       = ', d_ptr%spin_map%axis_input%m
      if (d_ptr%spin_map%valid) then
        nl=nl+1; write (lines(nl), '(2x, a, 16x, a, 34x, a)') 'Axes:', 'Initial', 'Final'
        nl=nl+1; write (lines(nl), '(a, 3f12.8, 5x, 3f12.8)') '  L-axis:', d_ptr%spin_map%axis0%l, d_ptr%spin_map%axis1%l
        nl=nl+1; write (lines(nl), '(a, 3f12.8, 5x, 3f12.8)') '  N-axis:', d_ptr%spin_map%axis0%n0, d_ptr%spin_map%axis1%n0
        nl=nl+1; write (lines(nl), '(a, 3f12.8, 5x, 3f12.8)') '  M-axis:', d_ptr%spin_map%axis0%m, d_ptr%spin_map%axis1%m
      else
        nl=nl+1; lines(nl) = 'Spin calculation not valid since: ' // d_ptr%why_invalid
      endif
    endif
    nl=nl+1; write(lines(nl), lmt)    '%exists            = ', d_ptr%exists
    nl=nl+1; write(lines(nl), lmt)    '%good_model        = ', d_ptr%good_model
    nl=nl+1; write(lines(nl), lmt)    '%good_design       = ', d_ptr%good_design
    nl=nl+1; write(lines(nl), lmt)    '%good_base         = ', d_ptr%good_base 
    nl=nl+1; write(lines(nl), lmt)    '%good_meas         = ', d_ptr%good_meas
    nl=nl+1; write(lines(nl), lmt)    '%good_ref          = ', d_ptr%good_ref
    nl=nl+1; write(lines(nl), lmt)    '%good_user         = ', d_ptr%good_user
    nl=nl+1; write(lines(nl), lmt)    '%good_opt          = ', d_ptr%good_opt
    nl=nl+1; write(lines(nl), lmt)    '%good_plot         = ', d_ptr%good_plot
    nl=nl+1; write(lines(nl), lmt)    '%useit_plot        = ', d_ptr%useit_plot
    nl=nl+1; write(lines(nl), '(a, l1, 3x, a)')    '%useit_opt         = ', d_ptr%useit_opt, tao_optimization_status(d_ptr)

    if (d_ptr%why_invalid /= '') then
      nl=nl+1; lines(nl) = 'Model value is invalid since: ' // d_ptr%why_invalid
    endif

    if (d_ptr%d1%d2%name(1:4) == 'ping') call show_ping(d_ptr%d1%d2%ix_universe)

  ! Else show the d1_data info.

  elseif (size(d1_array) == 1 .and. word1 /= '*@*') then

    d1_ptr => d1_array(1)%d1
    if (size(s%u) > 1) then
      nl=nl+1; write(lines(nl), '(2(a, i0))') 'Universe: ', d1_ptr%d2%ix_universe, '  of: ', ubound(s%u, 1)
    endif
    
    nl=nl+1; write(lines(nl), '(2a)') 'Data name: ', trim(d1_ptr%d2%name) // '.' // d1_ptr%name

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
      name = tao_constraint_type_name(d_ptr)
      if (substr(d_ptr%data_type,1,11) /= 'expression:') then
        n_name  = max(n_name,  len_trim(name))
        n_start = max(n_start, len_trim(d_ptr%ele_start_name))
        n_ref   = max(n_ref,   len_trim(d_ptr%ele_ref_name))
        n_ele   = max(n_ele,   len_trim(d_ptr%ele_name))
      endif
    enddo

    ! Write header
    ! Element names are left justified and real quantities are right justified

    line1 = ''; line2 = ''
    n=10+n_name;    line2(n:) = 'Ref_Ele'  ! n = i4 + 2x + n_name + 2x + 1
    n=n+n_ref+2;   line2(n:) = 'Start_Ele'
    n=n+n_start+2; line2(n:) = 'Ele'
    n=n+n_ele+12;  line2(n:) = 'Meas           Model          Design | Opt  Plot'
                   line1(n:) = '                                     |   Useit'

    nl=nl+1; lines(nl) = line1
    nl=nl+1; lines(nl) = line2

    ! if a range is specified, show the data range   

    call re_allocate (lines, nl+100+2*size(d_array), .false.)

    fmt  = '(i5, 4(2x, a), 3es16.7, 2l6)'
    fmt2 = '(5x, 4(2x, a), 3es16.7, 2l6)'

    do i = 1, size(d_array)
      d_ptr => d_array(i)%d
      if (.not. d_ptr%exists) cycle
      name = tao_constraint_type_name(d_ptr)
      if (substr(d_ptr%data_type,1,11) == 'expression:') then
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

    if (d1_ptr%d2%name(1:4) == 'ping') call show_ping(d1_ptr%d2%ix_universe)

  ! else if a single d2 structure

  elseif (size(d2_array) == 1 .and. word1 /= '*@*') then

    d2_ptr => d2_array(1)%d2

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

    if (d2_ptr%data_read_in) then
      nl=nl+1; write(lines(nl), amt)  '%data_file_name    = ', quote(d2_ptr%data_file_name)
      nl=nl+1; write(lines(nl), amt)  '%data_date         = ', quote(d2_ptr%data_date)
    endif

    if (d2_ptr%ref_read_in) then
      nl=nl+1; write(lines(nl), amt)  '%ref_file_name    = ', quote(d2_ptr%data_file_name)
      nl=nl+1; write(lines(nl), amt)  '%ref_date         = ', quote(d2_ptr%data_date)
    endif

    if (any(d2_ptr%descrip /= ' ')) then
      call re_allocate (lines, nl+100+size(d2_ptr%descrip), .false.)
      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = 'Descrip:'
      do i = 1, size(d2_ptr%descrip)
        if (d2_ptr%descrip(i) /= ' ') then
          nl=nl+1; write(lines(nl), '(i5, 2a)') i, ': ', d2_ptr%descrip(i)
        endif
      enddo
    endif

    if (d2_ptr%name(1:4) == 'ping') call show_ping(d2_ptr%ix_universe)

  ! Else several d2 structures

  elseif (size(d2_array) > 0) then

    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), '(a, t40, a)') '  Name', 'Using for Optimization'

    found = .false.
    do i = 1, size(d2_array)
      d2_ptr => d2_array(i)%d2
      if (d2_ptr%name == ' ') cycle
      if (d2_ptr%name(1:4) == 'ping') found = .true.
      call tao_data_show_use (d2_ptr, lines, nl)
    enddo

    if (found) call show_ping (d2_ptr%ix_universe)

  ! error

  else
    nl=1; lines(1) = 'NO MATCHING DATA FOUND.'
    return
  endif

!----------------------------------------------------------------------
! debug
! Temp code used for debugging. Also see 'internal'

case ('debug')

  print *, 'Absolute_time_tracking: ', bmad_com%absolute_time_tracking
  print *, 'E_tot: ', lat%ele(1)%value(E_tot$)

!----------------------------------------------------------------------
! derivative

case ('derivative')

  do_calc = .false.
  word1 = ''    ! data
  word2 = ''    ! variables

  do
    call tao_next_switch (what2, ['-derivative_recalc'], .true., switch, err)

    if (err) return
    if (switch == '') exit

    select case (switch)
    case ('-derivative_recalc')
      do_calc = .true.  

    case default
      if (word1 == '') then
        word1 = switch
      elseif (word2 == '') then
        word2 = switch
      else
        call out_io (s_error$, r_name, 'EXTRA STUFF ON LINE: ' // attrib0)
        return
      endif
    end select
  enddo

  if (word1 == '') word1 = '*@*'
  if (word2 == '') word2 = '*'

  call tao_find_data (err, word1, d_array = d_array);  if (err) return
  call tao_find_var (err, word2, v_array = v_array);  if (err) return

  call tao_dmodel_dvar_calc(do_calc, err);  if (err) return

  i1 = 0
  do id = 1, size(d_array)
    if (d_array(id)%d%ix_dmodel == -1) cycle
    i1 = i1 + 1
  enddo

  i2 = 0
  do iv = 1, size(v_array)
    if (v_array(iv)%v%ix_dvar == 0) cycle
    i2 = i2 + 1
  enddo

  call re_allocate (lines, nl+i1*i2+10, .false.)

  found = .false.
  nl=nl+1; lines(nl) = 'Data                          Variable                         Derivative   ix_dat  ix_var'
  do id = 1, size(d_array)
    do iv = 1, size(v_array)
      d_ptr => d_array(id)%d
      v_ptr => v_array(iv)%v
      u => s%u(d_ptr%d1%d2%ix_universe)
      jd = d_ptr%ix_dmodel
      jv = v_ptr%ix_dvar
      if (jd > 0 .and. jv > 0) then
        nl=nl+1; write(lines(nl), '(2a30, es14.5, 2i6)') tao_datum_name(d_ptr), &
                                  tao_var1_name(v_ptr), u%dModel_dVar(jd, jv), jd, jv
        found = .true.
      endif
    enddo
  enddo

  if (.not. found) then
    nl=nl+1; lines(nl) = 'No Derivative(s) Found. Note: Derivatives are only calculated by Tao for'
    nl=nl+1; lines(nl) = 'Data and variables that are used in an optimization.'
  endif
  
!----------------------------------------------------------------------
! dynamic_aperture

case ('dynamic_aperture')
  da => u%dynamic_aperture

  if (.not. allocated(da%scan)) then
    nl=nl+1; lines(nl) ='No dynamic aperture specified for this universe'
    return
  endif
  
  ! Count lines needed
  i1 = size(da%scan) * (da%param%n_angle + 20)
  call re_allocate (lines, nl+i1, .false.)

  nl=nl+1; write(lines(nl), '(a, a)')       'da_param%start_ele:     ', da%param%start_ele
  nl=nl+1; write(lines(nl), '(a, i10)')     'da_param%n_angle:       ', da%param%n_angle
  nl=nl+1; write(lines(nl), '(a, f10.6)')   'da_param%min_angle:     ', da%param%min_angle
  nl=nl+1; write(lines(nl), '(a, f10.6)')   'da_param%max_angle:     ', da%param%max_angle
  nl=nl+1; write(lines(nl), '(a, i10)')     'da_param%n_turn:        ', da%param%n_turn
  nl=nl+1; write(lines(nl), '(a, f10.6)')   'da_param%x_init:        ', da%param%x_init
  nl=nl+1; write(lines(nl), '(a, f10.6)')   'da_param%y_init:        ', da%param%y_init
  nl=nl+1; write(lines(nl), '(a, f10.6)')   'da_param%rel_accuracy:  ', da%param%rel_accuracy
  nl=nl+1; write(lines(nl), '(a, f10.6)')   'da_param%abs_accuracy:  ', da%param%abs_accuracy
  nl=nl+1; write(lines(nl), rmt)            'a_emit:                 ', da%a_emit
  nl=nl+1; write(lines(nl), rmt)            'b_emit:                 ', da%b_emit
  nl=nl+1; write(lines(nl), rmt)            'ellipse_scale:          ', da%ellipse_scale

  do k = 1, size(da%scan)
    da_scan => da%scan(k) 
    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), '(a, 99f11.6)') 'pz:        ', da%pz(k)
    if (.not. allocated(da_scan%point)) then
      nl=nl+1; write(lines(nl), '(a)') 'aperture not calculated for this universe'
    else
      branch => pointer_to_branch(pointer_to_ele(u%design%lat, da_scan%ref_orb%ix_ele, da_scan%ref_orb%ix_branch))
      nl=nl+1; write(lines(nl), '(a, 6es14.5)') 'ref_orb%vec:   ', da_scan%ref_orb%vec
      nl=nl+1; write(lines(nl), '(2a15)') 'aperture.x', 'aperture.y' 
      do j = 1, size(da_scan%point)
        da_pt => da_scan%point(j)
        nl=nl+1; write(lines(nl), '(2f11.6, 6x, a, 3x, a)')  da_pt%x, da_pt%y, &
                          coord_state_name(da_pt%plane), trim(branch%ele(da_pt%ix_ele)%name)
      enddo
    endif
  enddo

!----------------------------------------------------------------------
! element

case ('element')

  print_floor = .false.
  print_taylor = .false.
  print_field = no$
  print_attributes = .false.
  print_control = short$
  print_data = .false.
  print_wall = .false.
  print_rad  = .false.
  xfer_mat_print = 0
  print_slaves = .true.
  print_super_slaves = .true.
  lat_type = model$
  print_ptc = .false.
  attrib0 = ''
  name = ''

  do
    call tao_next_switch (what2, [character(16):: '-taylor', '-em_field', &
                '-all', '-data', '-design', '-no_slaves', '-wall', '-base', &
                '-field', '-floor_coords', '-xfer_mat', '-ptc', '-everything', &
                '-attributes', '-no_super_slaves', '-radiation_kick'], .true., switch, err)
    if (err) return
    select case (switch)
    case ('');                  exit
    case ('-xfer_mat');         xfer_mat_print = 6
    case ('-floor_coords');     print_floor = .true.
    case ('-taylor');           print_taylor = .true.
    case ('-design');           lat_type = design$
    case ('-base');             lat_type = base$
    case ('-em_field');         print_field = all$  ! Old style. Use "-field".
    case ('-field');            print_field = all$
    case ('-attributes');       print_attributes = .true.
    case ('-data');             print_data = .true.
    case ('-no_slaves');        print_slaves = .false.
    case ('-no_super_slaves');  print_super_slaves = .false.
    case ('-wall');             print_wall = .true.
    case ('-radiation_kick');   print_rad = .true.
    case ('-ptc');              print_ptc = .true.
    case ('-everything', '-all')
      print_attributes = .true.
      xfer_mat_print = 6
      print_taylor = .true.
      print_floor = .true.
      if (print_field == no$) print_field = short$
      print_wall = .true.
      print_rad  = .true.
      print_control = all$
    case default
      if (attrib0 /= '') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON LINE: ' // switch)
        return
      endif
      name = upcase(switch)
    end select
  enddo

  call tao_pick_universe (name, name, picked_uni, err, ix_u)
  if (err) return
  u => s%u(ix_u)

  ! Wildcard: show all elements.

  if (index(name, '*') /= 0 .or. index(name, '%') /= 0 .or. &
                (name(1:2) /= 'S:' .and. index(name, ':') /= 0) .or. &
                count(picked_uni) > 1) then

    n_tot = 0
    do i_uni = lbound(s%u, 1), ubound(s%u, 1)
      if (.not. picked_uni(i_uni)) cycle
      call tao_locate_elements (name, i_uni, eles, err, ignore_blank = .true.)
      if (err) return
      lat => s%u(i_uni)%model%lat
      do i = 1, size(eles)
        ele => eles(i)%ele
        if (.not. print_slaves .and. (ele%slave_status == super_slave$ .or. ele%slave_status == multipass_slave$)) cycle
        if (.not. print_super_slaves .and. ele%slave_status == super_slave$) cycle
        if (size(lines) < nl+100) call re_allocate (lines, nl+200, .false.)
        n_tot = n_tot + 1
        if (size(ele%branch%lat%branch) == 1) then
          write (str, '(i10)') ele%ix_ele
        else
          write (str, '(i0, a, i0)') ele%ix_branch, '>>', ele%ix_ele
          str(1:10) = adjustr(str(1:10))
        endif
        if (count(picked_uni) > 1) then
          nl=nl+1; write(lines(nl), '(a10, 2x, i0, 2a, f14.3)') str, i_uni, '@', ele%name, ele%s
        else
          nl=nl+1; write(lines(nl), '(a10, 2x, a, f14.3)') str, ele%name, ele%s
        endif
      enddo
    enddo

    deallocate(eles)
    nl=nl+1; write(lines(nl), '(a, i0)') 'Number of Matches: ', n_tot

    if (nl == 0) then
      lines(1) = '*** No Matches to Name Found ***'
      return
    endif

    result_id = 'element:*'
    return
  endif

  ! No wildcard case...
  ! Normal: Show the element info

  call tao_locate_elements (name, ix_u, eles, err, lat_type)
  if (err) return
  ele => eles(1)%ele

  tao_lat => tao_pointer_to_tao_lat (u, lat_type)
  branch => tao_lat%lat%branch(ele%ix_branch)
  tao_branch => tao_lat%tao_branch(ele%ix_branch)

  ! Show data associated with this element

  if (print_data) then
    call show_ele_data (u, ele, lines, nl)
    result_id = 'element:data'
    return
  endif

  if (print_ptc) then
    if (.not. associated (ele%ptc_fibre)) then
      nl=nl+1; lines(nl) = 'No Fibre associated  with element.'
      return
    endif
    call type_ptc_fibre (ele%ptc_fibre, .true., alloc_lines, n)
    if (size(lines) < nl+n+100) call re_allocate (lines, nl+n+100, .false.)
    lines(nl+1:nl+n) = alloc_lines(1:n)
    nl = nl + n
    result_id = 'element:ptc'
    return
  endif

  if (mass_of(branch%param%particle) /= 0) then
    ele%a%sigma = sqrt(ele%a%beta * tao_branch%modes_ri%a%emittance)
    ele%b%sigma = sqrt(ele%b%beta * tao_branch%modes_ri%b%emittance)
    ele%x%sigma = sqrt(tao_branch%lat_sigma(ele%ix_ele)%mat(1,1))
    ele%y%sigma = sqrt(tao_branch%lat_sigma(ele%ix_ele)%mat(3,3))
  endif

  twiss_out = s%global%phase_units
  if (lat%branch(ele%ix_branch)%param%particle == photon$) twiss_out = 0
  call type_ele (ele, print_attributes, xfer_mat_print, print_taylor, twiss_out, print_control, .true., &
            print_floor, print_field, print_wall, print_rad, lines = alloc_lines, n_lines = n)
  if (size(s%u) > 1) alloc_lines(1) = trim(alloc_lines(1)) // ',   Universe: ' // int_str(ix_u)

  if (size(lines) < nl+n+100) call re_allocate (lines, nl+n+100, .false.)
  lines(nl+1:nl+n) = alloc_lines(1:n)
  nl = nl + n

  stat = ele%lord_status
  if (stat /= ramper_lord$ .and. stat /= girder_lord$) then
    ele2 => ele
    if (ele2%lord_status == super_lord$) ele2 => pointer_to_slave(ele2, ele2%n_slave)
    orb = tao_lat%tao_branch(ele2%ix_branch)%orbit(ele2%ix_ele)
    if (orb%state == not_set$) then
      nl=nl+1; lines(nl) = 'Tracking: ' // trim(species_name(orb%species)) // ',   State: Orbit not computed.'
    else
      nl=nl+1; lines(nl) = ' '
      nl=nl+1; write(lines(nl), '(4a)') 'Tracking: ', trim(species_name(orb%species)), ',   State: ', trim(coord_state_name(orb%state))
      if (lat%branch(ele%ix_branch)%param%particle == photon$) then
        fmt  = '(2x, a, 2f15.8, f15.6, f11.6, 7x, a, f11.3)'
        fmt2 = '(2x, a, 2f15.8, a, es16.8)'
        nl=nl+1; lines(nl) = '         Position[mm]            V/C      Intensity      Phase  '
        nl=nl+1; write(lines(nl), fmt)  'X:  ', 1000*orb%vec(1), orb%vec(2), orb%field(1)**2, orb%phase(1), 'E: ', orb%p0c
        nl=nl+1; write(lines(nl), fmt)  'Y:  ', 1000*orb%vec(3), orb%vec(4), orb%field(2)**2, orb%phase(2), 'dE:', orb%p0c - ele%value(p0c$)
        nl=nl+1; write(lines(nl), fmt2) 'Z:  ', 1000*orb%vec(5), orb%vec(6)

        if (associated(ele%photon) .and. orb%state == alive$) then
          call to_surface_coords(orb, ele, orb2)
          nl=nl+1; lines(nl) = ''
          nl=nl+1; write(lines(nl), '(2x, a, 3f15.8)') 'Surface (x,y,z) [mm]:', orb2%vec(1:5:2)
        endif

      else
        z = (ele%ref_time - orb%t) * orb%beta * c_light
        dt = orb%t - ele%ref_time
        pc = orb%p0c * (1 + orb%vec(6))
        call convert_pc_to (pc, orb%species, e_tot = e_tot) 
        nl=nl+1; lines(nl) = '         Position[mm] Momentum[1E-3]        Spin   |'
        if (bmad_com%spin_tracking_on) then
          fmt  = '(2x, a, 2f15.8, x, a, a, es16.8, 2x, a, es12.5)'
          fmt2 = '(2x, a, 2f15.8, x, a, a, es16.8, 2x, a, f12.9)'
          nl=nl+1; write(lines(nl), fmt)  'X:  ', 1000*orb%vec(1:2), real_to_string(orb%spin(1), 12, 4, 8), '  | t_particle [sec]:      ', orb%t, 'E_tot:', e_tot
          nl=nl+1; write(lines(nl), fmt)  'Y:  ', 1000*orb%vec(3:4), real_to_string(orb%spin(2), 12, 4, 8), '  | t_part-t_ref [sec]:    ', dt,    'PC:   ', pc
          nl=nl+1; write(lines(nl), fmt2) 'Z:  ', 1000*orb%vec(5:6), real_to_string(orb%spin(3), 12, 4, 8), '  | (t_ref-t_part)*Vel [m]:', z,     'Beta: ', orb%beta
        else
          fmt  = '(2x, a, 2f15.8, 13x, a, es16.8, 2x, a, es12.5)'
          fmt2 = '(2x, a, 2f15.8, 13x, a, es16.8, 2x, a, f12.9)'
          nl=nl+1; write(lines(nl), fmt)  'X:  ', 1000*orb%vec(1:2), '  | t_particle [sec]:      ', orb%t, 'E_tot:', e_tot
          nl=nl+1; write(lines(nl), fmt)  'Y:  ', 1000*orb%vec(3:4), '  | t_part-t_ref [sec]:    ', dt,    'PC:   ', pc
          nl=nl+1; write(lines(nl), fmt2) 'Z:  ', 1000*orb%vec(5:6), '  | (t_ref-t_part)*Vel [m]:', z,     'Beta: ', orb%beta
        endif
        if (ele%key == rfcavity$ .or. ele%key == lcavity$) then
          if (ele%key == rfcavity$) then
            phase = ele%value(phi0$) + ele%value(phi0_multipass$) - &
                  (particle_rf_time (orb, ele, .false.) - rf_ref_time_offset(ele)) * ele%value(rf_frequency$)
          else
            phase = ele%value(phi0_err$) + ele%value(phi0$) + ele%value(phi0_multipass$) + &
                  (particle_rf_time (orb, ele, .false.) - rf_ref_time_offset(ele)) * ele%value(rf_frequency$)
          endif
          if (ele%tracking_method /= bmad_standard$) phase = phase + ele%value(phi0_autoscale$)
          nl=nl+1; write(lines(nl), '(2x, 2a)') 'Particle Phase relative to RF Phase (rad/2pi): ', real_str(phase, 9, 6)
        endif
      endif
    endif
  endif

  if (size(eles) > 1) then
    nl=nl+1; lines(nl) = ''
  endif

  if (size(eles) > 10) then
      nl=nl+1; write(lines(nl), '(a, i0)') 'NOTE: The number of other elements in the lattice with the same name is: ', size(eles) - 1
      nl=nl+1; write(lines(nl), '(a, i0)') '      To see a  list of these elements use a wild card character in the element name.'
  else
    call re_allocate (lines, nl+size(eles), .false.)
    do i = 2, size(eles)
      nl=nl+1; write(lines(nl), '(2a)') &
                'NOTE: There is another element with the same name at: ', trim(ele_loc_name(eles(i)%ele))
    enddo
  endif

!----------------------------------------------------------------------
! field

case ('emittance')

  ele_name = ''
  what_to_show = ''

  do 
    call tao_next_switch (what2, [character(16):: '-universe', '-element', '-xmatrix', '-sigma_matrix'], .true., switch, err)
    if (err) return
    if (switch == '') exit

    select case (switch)

    case ('-universe')
      read (aname, *, iostat = ios) ix
      u => tao_pointer_to_universe(ix)
      if (aname == '' .or. ios /= 0 .or. .not. associated(u)) then
        nl=1; lines(1) = 'CANNOT READ OR OUT-OF RANGE "-universe" ARGUMENT'
        return
      endif

    case ('-element')
      call tao_next_word(what2, ele_name)

    case ('-xmatrix', '-sigma_matrix')
      what_to_show = switch

    case default
      call out_io (s_error$, r_name, 'EXTRA STUFF ON LINE: ' // attrib0)
      return
    end select

  enddo

  if (ele_name == '') then
    ele => u%model%lat%branch(s%global%default_branch)%ele(0)
  else
    call tao_locate_elements (ele_name, u%ix_uni, eles, err)
    if (err) return
    ele => eles(1)%ele
  endif

  !

  if (what_to_show == '-xmatrix') then
    if (.not. associated(ele%rad_map)) then
      nl=nl+1; lines(nl) = 'No radiation matrices associated with element.'
      return
    endif

    ri => ele%rad_map

    nl=nl+1; write (lines(nl), '(10x, 30x, a, 64x, a)') 'Upstream half', 'Downstream half'

    nl=nl+1; lines(nl) = ''
    nl=nl+1; write (lines(nl), '(a10, 6es12.4, 4x, 6es12.4)') 'Ref Orb:  ', ri%rm0%ref_orb, ri%rm1%ref_orb
    nl=nl+1; write (lines(nl), '(a10, 6es12.4, 4x, 6es12.4)') 'Damp Vec: ', ri%rm0%xfer_damp_vec, ri%rm1%xfer_damp_vec


    nl=nl+1; lines(nl) = ''
    nl=nl+1; write (lines(nl), '(a10, 6es12.4, 4x, 6es12.4)') 'Damp Mat: ', ri%rm0%xfer_damp_mat(1,:), ri%rm1%xfer_damp_mat(1,:)
    do i = 2, 6
      nl=nl+1; write (lines(nl), '(10x, 6es12.4, 4x, 6es12.4)') ri%rm0%xfer_damp_mat(i,:), ri%rm1%xfer_damp_mat(i,:)
    enddo

    nl=nl+1; lines(nl) = ''
    nl=nl+1; write (lines(nl), '(a10, 6es12.4, 4x, 6es12.4)') 'Stoc Mat: ', ri%rm0%stoc_mat(1,:), ri%rm1%stoc_mat(1,:)
    do i = 2, 6
      nl=nl+1; write (lines(nl), '(10x, 6es12.4, 4x, 6es12.4)') ri%rm0%stoc_mat(i,:), ri%rm1%stoc_mat(i,:)
    enddo

    return
  endif

  !

  tao_branch => u%model%tao_branch(ele%ix_branch)
  mode_m => tao_branch%modes_ri

  u%model%high_e_lat = u%model%lat
  ele2 => u%model%high_e_lat%branch(ele%ix_branch)%ele(ele%ix_ele)
  call emit_6d (ele2, .false., mode_6d_no_vert, sig0_mat)
  call emit_6d (ele2, .true., mode_6d, sig0_mat)
  call radiation_integrals (u%model%lat, tao_branch%orbit, tao_branch%modes_ri, tao_branch%ix_rad_int_cache, ele%ix_branch)

  if (.not. associated(branch%ptc%m_t_layout)) then
    call out_io (s_info$, r_name, 'Note: Creating PTC layout (equivalent to "ptc init").')
    call lat_to_ptc_layout (lat)
  endif

  r = ptc_com%vertical_kick
  ptc_com%vertical_kick = 0
  call ptc_emit_calc (branch%ele(ele%ix_ele), mode_ptc_no_vert, sig_mat, orb)
  ptc_com%vertical_kick = 1
  call ptc_emit_calc (branch%ele(ele%ix_ele), mode_ptc, sig_mat, orb)
  ptc_com%vertical_kick = r

  nl=nl+1; lines(nl) = '      |       Photon Opening Angle Included       |      Photon Opening Angle Ignored         |'
  nl=nl+1; lines(nl) = ' Mode |   PTC_Emit    Bmad_6D_Emit   Rad_Int_Emit | PTC_Emit      Bmad_6D_Emit   Rad_Int_Emit |'
  nl=nl+1; write(lines(nl), '(1x, a, 2x, 6es15.7)') 'A', mode_ptc%a%emittance, mode_6d%a%emittance, mode_m%a%emittance, &
                                             mode_ptc_no_vert%a%emittance, mode_6d_no_vert%a%emittance_no_vert, mode_m%a%emittance_no_vert
  nl=nl+1; write(lines(nl), '(1x, a, 2x, 6es15.7)') 'B', mode_ptc%b%emittance, mode_6d%b%emittance, mode_m%b%emittance, &
                                             mode_ptc_no_vert%b%emittance, mode_6d_no_vert%b%emittance_no_vert, mode_m%b%emittance_no_vert
  nl=nl+1; write(lines(nl), '(1x, a, 2x, 6es15.7)') 'C', mode_ptc%z%emittance, mode_6d%z%emittance, mode_m%z%emittance, &
                                             mode_ptc_no_vert%z%emittance, mode_6d_no_vert%z%emittance_no_vert, mode_m%z%emittance

  nl=nl+1; lines(nl) = ''
  nl=nl+1; write(lines(nl), '(a, 3es12.4)') 'J_damp:    ', mode_6d%a%j_damp, mode_6d%b%j_damp, mode_6d%z%j_damp
  nl=nl+1; write(lines(nl), '(a, 3es12.4)') 'Alpha_damp:', mode_6d%a%alpha_damp, mode_6d%b%alpha_damp, mode_6d%z%alpha_damp
  nl=nl+1; write(lines(nl), '(a, 3es12.4)') 'SigE/E:    ', mode_6d%sige_e
  nl=nl+1; write(lines(nl), '(a, 3es12.4)') 'Sig_z:     ', mode_6d%sig_z
  nl=nl+1; write(lines(nl), '(a, 3es12.4)') 'Alpha_p:   ', momentum_compaction(branch)


  if (what_to_show == '-sigma_matrix') then
    nl=nl+1; lines(nl) = ''
    nl=nl+1; write (lines(nl), '(a)') 'PTC Sigma Mat (w/vert angle):'
    do i = 1, 6
      nl=nl+1; write (lines(nl), '(5x, 6es12.4, 4x, 6es12.4)') sig_mat(i,:)
    enddo

    nl=nl+1; lines(nl) = ''
    nl=nl+1; write (lines(nl), '(a)') 'Bmad Sigma Mat (w/vert angle):'
    do i = 1, 6
      nl=nl+1; write (lines(nl), '(5x, 6es12.4, 4x, 6es12.4)') sig0_mat(i,:)
    enddo
  endif

!----------------------------------------------------------------------
! field

case ('field')

  lat_type = model$
  orb%vec = 0
  orb%t = 0
  s_pos = 0
  show_all = .false.   ! Show derivatives?
  what_to_show = 'standard'
  s_fmt = 'relative'
  n_count = 0          ! Counter for non-switch args. 

  do
    call tao_next_switch (what2, [character(16):: '-derivatives', '-grid_pt', '-percent_len', '-absolute_s'], &
                                                                       .true., switch, err, .true.)
    if (err) return
    select case (switch)
    case ('');               exit
    case ('-derivatives');                show_all = .true.
    case ('-percent_len', '-absolute_s'); s_fmt = switch
    case ('-grid_pt');                    what_to_show = switch 
    case default
      n_count = n_count + 1

      select case (n_count)
      case (1)
        call str_upcase(ele_name, switch)
        call tao_pick_universe (ele_name, ele_name, picked_uni, err, ix_u)
        if (err) return
        u => s%u(ix_u)
        call tao_locate_elements (ele_name, ix_u, eles, err, lat_type)
        if (err) return
        if (size(eles) == 0) then
          nl=1; lines(nl) = 'No lattice elements found that match to: ' // quote(ele_name)
          return
        endif

      case (2)
        call tao_evaluate_expression(switch, 1, .false., value, err)
        if (err) then
          nl=1; lines(nl) = 'Bad X value';  result_id = 'field:bad-x'
          return
        endif
        orb%vec(1) = value(1)

      case (3)
        call tao_evaluate_expression(switch, 1, .false., value, err)
        if (err) then
          nl = 1; lines(1) = 'Bad Y value';  result_id = 'field:bad-y'
          return
        endif
        orb%vec(3) = value(1)

      case (4)
        call tao_evaluate_expression(switch, 1, .false., value, err)
        if (err) then
          nl = 1; lines(1) = 'Bad S value';  result_id = 'field:bad-s'
          return
        endif
        s_pos = value(1)

      case (5)
        call tao_evaluate_expression(switch, 1, .false., value, err)
        if (err) then
          nl = 1; lines(1) = 'Bad Z or T value';  result_id = 'field:bad-t'
          return
        endif

        if (bmad_com%absolute_time_tracking) then
          orb%t = value(1)
        else
          orb%vec(5) = value(1)
        endif

      case (6)
        nl = 1; lines(1) = 'Extra stuff on line. Remember: expressions with blanks need to be quoted'
        result_id = 'field:extra'
        return
      end select

    end select
  enddo

  if (n_count == 0) then
    nl=1; lines(nl) = 'No lattice element name given.'
    return
  endif

  !

  do i = 1, size(eles)
    ele => eles(i)%ele

    if (what_to_show == '-grid_pt') then
      if (.not. associated (ele%grid_field)) then
        call out_io (s_error$, r_name, 'NO GRID FIELD ASSOCIATED WITH ELEMENT: ' // ele%name)
        return
      endif
      ix0 = nint(orb%vec(1))
      ix1 = nint(orb%vec(3))
      ix2 = nint(s_pos)
      g_pt => ele%grid_field(1)%ptr%pt(ix0, ix1, ix2)
      nl=nl+1; 
      cycle
    endif

    if (size(lines) < nl + 20) call re_allocate (lines, nl + 100) 

    !

    call init_coord (orb, orb, ele, downstream_end$)

    select case (s_fmt)
    case ('-percent_len');   s_pos = s_pos * ele%value(l$)
    case ('-absolute_s');    s_pos = s_pos - ele%s_start
    case ('relative');       s_pos = s_pos
    end select

    time1 = particle_rf_time(orb, ele, .true., s_pos)
    dt = 1d-12

    call em_field_calc (ele, ele%branch%param, s_pos, orb, .false., field, calc_dfield = show_all, &
                                                          err_flag = err, rf_time = time1);  if (err) return
    call em_field_calc (ele, ele%branch%param, s_pos, orb, .false., field0, rf_time = time1-dt)
    call em_field_calc (ele, ele%branch%param, s_pos, orb, .false., field1, rf_time = time1+dt)
    field1%E = (field1%E - field0%E) / (2 * dt)
    field1%B = (field1%B - field0%B) / (2 * dt)

    if (i > 1) then
      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = '=================================================='
    endif
    nl=nl+1; lines(nl) = trim(ele_loc_name(ele)) // ': ' // trim(ele%name) 
    nl=nl+1; write (lines(nl), '(a24, 2x, a)') '  B [T]:', reals_to_string(field%B, 14, 2, 6, 6)
    nl=nl+1; write (lines(nl), '(a24, 2x, a)') 'E [V/m]:', reals_to_string(field%E, 14, 2, 6, 2)

    if (show_all) then
      nl=nl+1; lines(nl) = ''
      nl=nl+1; write (lines(nl), '(a24, 3es16.8)') 'dB/dx [T/m]:', field%dB(1,:)
      nl=nl+1; write (lines(nl), '(a24, 3es16.8)') 'dB/dy [T/m]:', field%dB(2,:)
      nl=nl+1; write (lines(nl), '(a24, 3es16.8)') 'dB/dz [T/m]:', field%dB(3,:)

      nl=nl+1; lines(nl) = ''
      nl=nl+1; write (lines(nl), '(a24, 3es16.8)') 'dB/dt [T/sec]:', field1%B
      nl=nl+1; write (lines(nl), '(a24, 3es16.8, a, es16.8)') 'Curl_B, Div_B [T/m]:', field%dB(3,2) - field%dB(2,3), &
                                            field%dB(1,3) - field%dB(3,1), field%dB(2,1) - field%dB(1,2), ',', &
                                            field%dB(1,1) + field%dB(2,2) + field%dB(3,3)
      nl=nl+1; lines(nl) = ''
      nl=nl+1; write (lines(nl), '(a24, 3es16.8)') 'dE/dx [V/m^2]:', field%dE(1,:)
      nl=nl+1; write (lines(nl), '(a24, 3es16.8)') 'dE/dy [V/m^2]:', field%dE(2,:)
      nl=nl+1; write (lines(nl), '(a24, 3es16.8)') 'dE/dz [V/m^2]:', field%dE(3,:)

      nl=nl+1; lines(nl) = ''
      nl=nl+1; write (lines(nl), '(a24, 3es16.8)') 'dE/dt [V/m/sec]:', field1%E
      nl=nl+1; write (lines(nl), '(a24, 3es16.8, a, es16.8)') 'Curl_E, Div_E [V/m^2]:', field%dE(3,2) - field%dE(2,3), &
                                            field%dE(1,3) - field%dE(3,1), field%dE(2,1) - field%dE(1,2), ',', &
                                            field%dE(1,1) + field%dE(2,2) + field%dE(3,3)
    endif
  enddo

!----------------------------------------------------------------------
! global

case ('global')

  call ran_default_state (get_state = ran_state)
  what_to_show = 'global'

  do
    call tao_next_switch (what2, [character(20):: '-optimization', '-bmad_com', &
                    '-csr_param', '-space_charge_com', '-ran_state', '-ptc_com', '-internal'], .true., switch, err)
    if (err) return

    select case (switch)
    case ('')
      exit
    case ('-optimization', '-bmad_com', '-csr_param', '-space_charge_com', '-ran_state', '-ptc_com', '-internal')
      what_to_show = switch
    case default
      call out_io (s_error$, r_name, 'EXTRA STUFF ON LINE: ' // switch)
      return
    end select
  enddo

  select case (what_to_show)
  case ('global')
    nl=nl+1; lines(nl) = 'Tao Global parameters [Note: To print optimizer globals use: "show optimizer"]'
    nl=nl+1; write(lines(nl), lmt) '  %beam_timer_on                 = ', s%global%beam_timer_on
    nl=nl+1; write(lines(nl), imt) '  %bunch_to_plot                 = ', s%global%bunch_to_plot
    nl=nl+1; write(lines(nl), imt) '  %datum_err_messages_max        = ', s%global%datum_err_messages_max
    nl=nl+1; write(lines(nl), lmt) '  %cmd_file_abort_on_error       = ', s%global%cmd_file_abort_on_error
    nl=nl+1; write(lines(nl), lmt) '  %concatenate_maps              = ', s%global%concatenate_maps
    nl=nl+1; write(lines(nl), rmt) '  %beam_dead_cutoff              = ', s%global%beam_dead_cutoff
    nl=nl+1; write(lines(nl), rmt) '  %delta_e_chrom                 = ', s%global%delta_e_chrom
    nl=nl+1; write(lines(nl), lmt) '  %disable_smooth_line_calc      = ', s%global%disable_smooth_line_calc
    nl=nl+1; write(lines(nl), lmt) '  %draw_curve_off_scale_warn     = ', s%global%draw_curve_off_scale_warn
    nl=nl+1; write(lines(nl), lmt) '  %init_lat_sigma_from_beam      = ', s%global%init_lat_sigma_from_beam
    nl=nl+1; write(lines(nl), lmt) '  %label_lattice_elements        = ', s%global%label_lattice_elements
    nl=nl+1; write(lines(nl), lmt) '  %label_keys                    = ', s%global%label_keys
    nl=nl+1; write(lines(nl), lmt) '  %lattice_calc_on               = ', s%global%lattice_calc_on
    nl=nl+1; write(lines(nl), rmt) '  %max_plot_time                 = ', s%global%max_plot_time
    nl=nl+1; write(lines(nl), imt) '  %n_threads                     = ', s%global%n_threads
    nl=nl+1; write(lines(nl), lmt) '  %only_limit_opt_vars           = ', s%global%only_limit_opt_vars
    nl=nl+1; write(lines(nl), lmt) '  %opt_match_auto_recalc         = ', s%global%opt_match_auto_recalc
    nl=nl+1; write(lines(nl), lmt) '  %opti_write_var_file           = ', s%global%opti_write_var_file
    nl=nl+1; write(lines(nl), lmt) '  %optimizer_var_limit_warn      = ', s%global%optimizer_var_limit_warn
    nl=nl+1; write(lines(nl), amt) '  %phase_units                   = ', angle_units_name(s%global%phase_units)
    nl=nl+1; write(lines(nl), lmt) '  %rad_int_calc_on               = ', s%global%rad_int_user_calc_on
    nl=nl+1; write(lines(nl), amt) '  %history_file                  = ', s%global%history_file
    nl=nl+1; write(lines(nl), lmt) '  %plot_on                       = ', s%global%plot_on
    nl=nl+1; write(lines(nl), lmt) '  %external_plotting             = ', s%global%external_plotting
    nl=nl+1; write(lines(nl), amt) '  %print_command                 = ', quote(s%global%print_command)
    nl=nl+1; write(lines(nl), amt) '  %prompt_string                 = ', quote(s%global%prompt_string)
    nl=nl+1; write(lines(nl), amt) '  %prompt_color                  = ', quote(s%global%prompt_color)
    nl=nl+1; write(lines(nl), amt) '  %random_engine                 = ', quote(s%global%random_engine)
    nl=nl+1; write(lines(nl), amt) '  %random_gauss_converter        = ', quote(s%global%random_gauss_converter)
    nl=nl+1; write(lines(nl), amt) '  %quiet                         = ', quote(s%global%quiet)

    nl=nl+1; write(lines(nl), amt) '  %random_engine (input)         = ', quote(s%global%random_engine)
    if (s%global%random_engine == '') then
      nl=nl+1; write(lines(nl), amt) '   random engine used            = ', quote(ran_engine_name(ran_state%engine))
    endif

    nl=nl+1; write(lines(nl), amt) '  %random_gauss_converter (input)= ', quote(s%global%random_gauss_converter)
    if (s%global%random_gauss_converter == '') then
      nl=nl+1; write(lines(nl), amt) '   random gauss converter used   = ', quote(ran_gauss_converter_name(ran_state%gauss_converter))
    endif

    nl=nl+1; write(lines(nl), imt) '  %random_seed (input)           = ', s%global%random_seed
    nl=nl+1; write(lines(nl), imt) '     random seed used              = ', ran_state%seed

    nl=nl+1; write(lines(nl), rmt) '  %random_sigma_cutoff           = ', s%global%random_sigma_cutoff

    nl=nl+1; write(lines(nl), lmt) '  %rf_on                         = ', s%global%rf_on
    nl=nl+1; write(lines(nl), amt) '  %track_type                    = ', quote(s%global%track_type)
    nl=nl+1; write(lines(nl), lmt) '  %var_limits_on                 = ', s%global%var_limits_on
    nl=nl+1; write(lines(nl), amt) '  %var_out_file                  = ', quote(s%global%var_out_file)
    nl=nl+1; write(lines(nl), lmt) '  %wait_for_CR_in_single_mode    = ', s%global%wait_for_CR_in_single_mode
    nl=nl+1; write(lines(nl), lmt) '  %blank_line_between_commands   = ', s%global%blank_line_between_commands

    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Tao Parameters:'
    nl=nl+1; write(lines(nl), imt) '  Universe index range:        = ', lbound(s%u, 1), ubound(s%u, 1)
    nl=nl+1; write(lines(nl), amt) '  default_universe:            = ', int_str(s%global%default_universe), '  ! Set using: "set default universe = ..."'
    nl=nl+1; write(lines(nl), amt) '  default_branch:              = ', int_str(s%global%default_branch),   '  ! Set using: "set default branch = ..."'
    nl=nl+1; write(lines(nl), imt) '  Number paused command files  = ', count(s%com%cmd_file%paused)
    nl=nl+1; write(lines(nl), amt) '  unique_name_suffix           = ', quote(s%init%unique_name_suffix)
    nl=nl+1; write(lines(nl), lmt) '  Combine_consecutive_elements_of_like_name = ', &
                                                s%com%combine_consecutive_elements_of_like_name

    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Tao command line startup arguments:'
    call write_this_arg (nl, lines, '  -beam_file', s%init%beam_file_arg)
    call write_this_arg (nl, lines, '  -beam_init_position_file', s%init%beam_init_position_file_arg)
    call write_this_arg (nl, lines, '  -command', s%init%command_arg)
    call write_this_arg (nl, lines, '  -building_wall_file', s%init%building_wall_file_arg)
    call write_this_arg (nl, lines, '  -prompt_color', s%init%prompt_color_arg)
    call write_this_arg (nl, lines, '  -data_file', s%init%data_file_arg)
    call write_this_arg (nl, lines, '  -disable_smooth_line_calc', s%init%disable_smooth_line_calc_arg)
    call write_this_arg (nl, lines, '  -debug', s%init%debug_arg)
    call write_this_arg (nl, lines, '  -geometry', s%init%geometry_arg)
    call write_this_arg (nl, lines, '  -hook_init_file', s%init%hook_init_file_arg)
    call write_this_arg (nl, lines, '  -init_file', s%init%init_file_arg)
    call write_this_arg (nl, lines, '  -lattice_file', s%init%lattice_file_arg)
    call write_this_arg (nl, lines, '  -log_startup', s%init%log_startup_arg)
    call write_this_arg (nl, lines, '  -no_stopping', s%init%no_stopping_arg)
    call write_this_arg (nl, lines, '  -noinit', s%init%noinit_arg)
    call write_this_arg (nl, lines, '  -noplot', s%init%noplot_arg)
    call write_this_arg (nl, lines, '  -nostartup', s%init%nostartup_arg)
    call write_this_arg (nl, lines, '  -no_rad_int', s%init%no_rad_int_arg)
    call write_this_arg (nl, lines, '  -plot_file', s%init%plot_file_arg)
    call write_this_arg (nl, lines, '  -reverse', s%init%reverse_arg)
    call write_this_arg (nl, lines, '  -rf_on', s%init%rf_on_arg)
    call write_this_arg (nl, lines, '  -quiet', s%init%quiet_arg)
    call write_this_arg (nl, lines, '  -slice_lattice', s%init%slice_lattice_arg)
    call write_this_arg (nl, lines, '  -start_branch_at', s%init%start_branch_at_arg)
    call write_this_arg (nl, lines, '  -startup_file', s%init%startup_file_arg)
    call write_this_arg (nl, lines, '  -symbol_import', s%init%symbol_import_arg)
    call write_this_arg (nl, lines, '  -var_file', s%init%var_file_arg)

  case ('-ran_state')
    nl=nl+1; write(lines(nl), imt) '  %ix              = ', ran_state%ix
    nl=nl+1; write(lines(nl), imt) '  %iy              = ', ran_state%iy
    nl=nl+1; write(lines(nl), lmt) '  %number_stored   = ', ran_state%number_stored
    nl=nl+1; write(lines(nl), lmt) '  %h_saved         = ', ran_state%h_saved
    nl=nl+1; write(lines(nl), imt) '  %engine          = ', ran_state%engine
    nl=nl+1; write(lines(nl), imt) '  %seed            = ', ran_state%seed
    nl=nl+1; write(lines(nl), rmt) '  %am              = ', ran_state%am
    nl=nl+1; write(lines(nl), imt) '  %gauss_converter = ', ran_state%gauss_converter
    nl=nl+1; write(lines(nl), rmt) '  %gauss_sigma_cut = ', ran_state%gauss_sigma_cut
    nl=nl+1; write(lines(nl), imt) '  %in_sobseq       = ', ran_state%in_sobseq

  case ('-optimization')
    call show_opt()

  case ('-bmad_com')
    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Bmad_com Parameters (use "set bmad_com" to change):'
    nl=nl+1; write(lines(nl), rmt) '  %max_aperture_limit              = ', bmad_com%max_aperture_limit
    nl=nl+1; write(lines(nl), rmt) '  %d_orb                           = ', bmad_com%d_orb
    nl=nl+1; write(lines(nl), rmt) '  %default_ds_step                 = ', bmad_com%default_ds_step
    nl=nl+1; write(lines(nl), rmt) '  %significant_length              = ', bmad_com%significant_length
    nl=nl+1; write(lines(nl), rmt) '  %rel_tol_tracking                = ', bmad_com%rel_tol_tracking
    nl=nl+1; write(lines(nl), rmt) '  %abs_tol_tracking                = ', bmad_com%abs_tol_tracking
    nl=nl+1; write(lines(nl), rmt) '  %rel_tol_adaptive_tracking       = ', bmad_com%rel_tol_adaptive_tracking
    nl=nl+1; write(lines(nl), rmt) '  %abs_tol_adaptive_tracking       = ', bmad_com%abs_tol_adaptive_tracking
    nl=nl+1; write(lines(nl), rmt) '  %autoscale_amp_rel_tol           = ', bmad_com%autoscale_amp_rel_tol
    nl=nl+1; write(lines(nl), rmt) '  %autoscale_amp_abs_tol           = ', bmad_com%autoscale_amp_abs_tol
    nl=nl+1; write(lines(nl), rmt) '  %autoscale_phase_tol             = ', bmad_com%autoscale_phase_tol
    nl=nl+1; write(lines(nl), rmt) '  %init_ds_adaptive_tracking       = ', bmad_com%init_ds_adaptive_tracking
    nl=nl+1; write(lines(nl), rmt) '  %min_ds_adaptive_tracking        = ', bmad_com%min_ds_adaptive_tracking
    nl=nl+1; write(lines(nl), rmt) '  %electric_dipole_moment          = ', bmad_com%electric_dipole_moment
    nl=nl+1; write(lines(nl), rmt) '  %synch_rad_scale                 = ', bmad_com%synch_rad_scale
    nl=nl+1; write(lines(nl), rmt) '  %sad_eps_scale                   = ', bmad_com%sad_eps_scale
    nl=nl+1; write(lines(nl), rmt) '  %sad_amp_max                     = ', bmad_com%sad_amp_max

    nl=nl+1; write(lines(nl), imt) '  %sad_n_div_max                   = ', bmad_com%sad_n_div_max
    nl=nl+1; write(lines(nl), iimt)'  %taylor_order                    = ', bmad_com%taylor_order, ' ! Input order. 0 => Use default'
    nl=nl+1; write(lines(nl), imt)'     Actual Taylor order in PTC: ', ptc_private%taylor_order_ptc
    nl=nl+1; write(lines(nl), imt) '  %default_integ_order             = ', bmad_com%default_integ_order

    nl=nl+1; write(lines(nl), lmt) '  %rf_phase_below_transition_ref   = ', bmad_com%rf_phase_below_transition_ref
    nl=nl+1; write(lines(nl), lmt) '  %sr_wakes_on                     = ', bmad_com%sr_wakes_on
    nl=nl+1; write(lines(nl), lmt) '  %lr_wakes_on                     = ', bmad_com%lr_wakes_on
    nl=nl+1; write(lines(nl), lmt) '  %auto_bookkeeper                 = ', bmad_com%auto_bookkeeper
    nl=nl+1; write(lines(nl), lmt) '  %high_energy_space_charge_on     = ', bmad_com%high_energy_space_charge_on
    nl=nl+1; write(lines(nl), lmt) '  %csr_and_space_charge_on         = ', bmad_com%csr_and_space_charge_on
    nl=nl+1; write(lines(nl), lmt) '  %spin_tracking_on                = ', bmad_com%spin_tracking_on
    nl=nl+1; write(lines(nl), lmt) '  %spin_sokolov_ternov_flipping_on = ', bmad_com%spin_sokolov_ternov_flipping_on
    nl=nl+1; write(lines(nl), lmt) '  %radiation_damping_on            = ', bmad_com%radiation_damping_on
    nl=nl+1; write(lines(nl), lmt) '  %radiation_zero_average          = ', bmad_com%radiation_zero_average
    nl=nl+1; write(lines(nl), lmt) '  %radiation_fluctuations_on       = ', bmad_com%radiation_fluctuations_on
    nl=nl+1; write(lines(nl), lmt) '  %conserve_taylor_maps            = ', bmad_com%conserve_taylor_maps
    nl=nl+1; write(lines(nl), lmt) '  %absolute_time_tracking          = ', bmad_com%absolute_time_tracking
    nl=nl+1; write(lines(nl), lmt) '  %absolute_time_ref_shift         = ', bmad_com%absolute_time_ref_shift
    nl=nl+1; write(lines(nl), lmt) '  %convert_to_kinetic_momentum     = ', bmad_com%convert_to_kinetic_momentum
    nl=nl+1; write(lines(nl), lmt) '  %aperture_limit_on               = ', bmad_com%aperture_limit_on

    if (allocated(lat%custom)) then
      nl=nl+1; lines(nl) = 'Custom lattice parameters defined in lattice file:'
      do i = 1, size(lat%custom)
        myname = attribute_name(def_parameter$, i+custom_attribute0$)
        if (myname(1:1) == '!') cycle
        nl= nl+1; write (lines(nl), rmt) '  parameter[' // trim(myname) // ']: ', lat%custom(i)
      enddo
    endif

  case ('-csr_param')
    nl=nl+1; lines(nl) = '"-csr_param" is now "-space_charge_com".'

  case ('-space_charge_com')
    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'space_charge_com Parameters (set by "set space_charge_com ..."):'
    nl=nl+1; write(lines(nl), rmt) '  %ds_track_step                  = ', space_charge_com%ds_track_step
    nl=nl+1; write(lines(nl), rmt) '  %dt_track_step                  = ', space_charge_com%dt_track_step
    nl=nl+1; write(lines(nl), rmt) '  %cathode_strength_cutoff        = ', space_charge_com%cathode_strength_cutoff
    nl=nl+1; write(lines(nl), rmt) '  %rel_tol_tracking               = ', space_charge_com%rel_tol_tracking
    nl=nl+1; write(lines(nl), rmt) '  %abs_tol_tracking               = ', space_charge_com%abs_tol_tracking
    nl=nl+1; write(lines(nl), rmt) '  %beam_chamber_height            = ', space_charge_com%beam_chamber_height
    nl=nl+1; write(lines(nl), rmt) '  %lsc_sigma_cutoff               = ', space_charge_com%lsc_sigma_cutoff
    nl=nl+1; write(lines(nl), rmt) '  %particle_sigma_cutoff          = ', space_charge_com%particle_sigma_cutoff
    nl=nl+1; write(lines(nl), imt) '  %space_charge_mesh_size         = ', space_charge_com%space_charge_mesh_size
    nl=nl+1; write(lines(nl), imt) '  %csr3d_mesh_size                = ', space_charge_com%csr3d_mesh_size
    nl=nl+1; write(lines(nl), imt) '  %n_bin                          = ', space_charge_com%n_bin
    nl=nl+1; write(lines(nl), imt) '  %particle_bin_span              = ', space_charge_com%particle_bin_span
    nl=nl+1; write(lines(nl), imt) '  %n_shield_images                = ', space_charge_com%n_shield_images
    nl=nl+1; write(lines(nl), imt) '  %sc_min_in_bin                  = ', space_charge_com%sc_min_in_bin
    nl=nl+1; write(lines(nl), lmt) '  %lsc_kick_transverse_dependence = ', space_charge_com%lsc_kick_transverse_dependence
    nl=nl+1; write(lines(nl), amt) '  %diagnostic_output_file         = ', quote(space_charge_com%diagnostic_output_file)

  case ('-ptc_com')
    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'PTC_com Parameters (set using "set ptc_com ..."):'
    nl=nl+1; write(lines(nl), rmt) '  %vertical_kick              = ', ptc_com%vertical_kick
    nl=nl+1; write(lines(nl), imt) '  %max_fringe_order           = ', ptc_com%max_fringe_order
    nl=nl+1; write(lines(nl), imt) '  %old_integrator             = ', ptc_com%old_integrator
    nl=nl+1; write(lines(nl), lmt) '  %exact_model                = ', ptc_com%exact_model
    nl=nl+1; write(lines(nl), lmt) '  %exact_misalign             = ', ptc_com%exact_misalign
    nl=nl+1; write(lines(nl), lmt) '  %translate_patch_drift_time = ', ptc_com%translate_patch_drift_time
    nl=nl+1; write(lines(nl), lmt) '  %print_info_messages        = ', ptc_com%print_info_messages

  ! Internal parameters are not of general interest.
  case ('-internal')
    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), imt) '  u%calc%srdt_for_data              = ', u%calc%srdt_for_data
    nl=nl+1; write(lines(nl), lmt) '  u%calc%rad_int_for_data           = ', u%calc%rad_int_for_data
    nl=nl+1; write(lines(nl), lmt) '  u%calc%rad_int_for_plotting       = ', u%calc%rad_int_for_plotting
    nl=nl+1; write(lines(nl), lmt) '  u%calc%chrom_for_data             = ', u%calc%chrom_for_data
    nl=nl+1; write(lines(nl), lmt) '  u%calc%chrom_for_plotting         = ', u%calc%chrom_for_plotting
    nl=nl+1; write(lines(nl), lmt) '  u%calc%lat_sigma_for_data         = ', u%calc%lat_sigma_for_data
    nl=nl+1; write(lines(nl), lmt) '  u%calc%lat_sigma_for_plotting     = ', u%calc%lat_sigma_for_plotting
    nl=nl+1; write(lines(nl), lmt) '  u%calc%dynamic_aperture           = ', u%calc%dynamic_aperture
    nl=nl+1; write(lines(nl), lmt) '  u%calc%one_turn_map               = ', u%calc%one_turn_map
    nl=nl+1; write(lines(nl), lmt) '  u%calc%lattice                    = ', u%calc%lattice
    nl=nl+1; write(lines(nl), lmt) '  u%calc%twiss                      = ', u%calc%twiss
    nl=nl+1; write(lines(nl), lmt) '  u%calc%track                      = ', u%calc%track
    nl=nl+1; write(lines(nl), lmt) '  u%calc%spin_matrices              = ', u%calc%spin_matrices
    nl=nl+1; write(lines(nl), lmt) '  bmad_com%debug                    = ', bmad_com%debug
  end select

!----------------------------------------------------------------------
! graph

case ('graph')

  print_debug = .false.
  print_rms = .false.
  if (allocated(graph)) deallocate(graph)

  do
    call tao_next_switch (what2, [character(8):: '-debug', '-rms'], .true., switch, err)
    if (err) return
    select case (switch)
    case ('')
      if (allocated(graph)) exit
      nl=1; lines(1) = 'Graph name is blank!'
      return
    case ('-debug')
      print_debug = .true.
    case default
      call tao_find_plots (err, switch, 'BOTH', graph = graph, blank_means_all = .true., only_visible = .false.)
      if (err) return
      if (size(graph) == 0) then
        nl=1; lines(1) = 'This is not a graph'
        return
      endif
    case ('-rms')
      print_rms = .true.
    end select
  enddo    

  ! RMS

  nl=nl+1; lines(nl) = '                         |          Line           |          Symbols'
  nl=nl+1; lines(nl) = 'Curves:                  |     Mean         RMS    |      Mean         RMS'

  do ig = 1, size(graph)
    g => graph(ig)%g
    if (.not. allocated(g%curve)) cycle
    do i = 1, size(g%curve)
      nl=nl+1; write(lines(nl), '(3x, a)') tao_curve_name(g%curve(i), .true.)
      if (allocated(g%curve(i)%x_line)) then
        call tao_curve_rms_calc(g%curve(i), 'LINE', rms, mean)
        write(lines(nl)(26:), '(2es12.4)') mean, rms
      endif
      if (allocated(g%curve(i)%x_symb)) then
        call tao_curve_rms_calc(g%curve(i), 'SYMBOL', rms, mean)
        write(lines(nl)(53:), '(2es12.4)') mean, rms
      endif
    enddo
  enddo

  ! Standard

  g => graph(1)%g
  fmt = '(a, f6.3)'

  if (associated(g%p%r)) then
    nl=nl+1; lines(nl) = 'Region.Graph: ' // trim(g%p%r%name) // '.' // trim(g%name)
  endif
  nl=nl+1; lines(nl) = 'Plot.Graph:   ' // trim(g%p%name) // '.' // trim(g%name)
  nl=nl+1; write(lines(nl), amt)  'type                             = ', quote(g%type)
  nl=nl+1; write(lines(nl), amt)  'title                            = ', quote(g%title)
  nl=nl+1; write(lines(nl), amt)  'title_suffix                     = ', quote(g%title_suffix)
  nl=nl+1; write(lines(nl), '(a, 4f10.2, 2x, a)') &
                                  'margin                           = ', g%margin
  nl=nl+1; write(lines(nl), '(a, 4f10.2, 2x, a)') &
                                  'scale_margin                     = ', g%scale_margin
  nl=nl+1; write(lines(nl), imt)  'box                              = ', g%box
  nl=nl+1; write(lines(nl), imt)  'ix_universe                      = ', g%ix_universe
  if (.not. allocated(g%curve)) then   ! %ix_branch is only used if there are no curves.
    nl=nl+1; write(lines(nl), imt)  'ix_branch                        = ', g%ix_branch
  endif
  nl=nl+1; write(lines(nl), lmt)  'is_valid                         = ', g%is_valid

  nl=nl+1; write(lines(nl), rmt)  'x_axis_scale_factor              = ', g%x_axis_scale_factor
  nl=nl+1; write(lines(nl), rmt)  'symbol_size_scale                = ', g%symbol_size_scale
  nl=nl+1; write(lines(nl), amt)  'text_legend_origin%x,y,units     = ', real_str(g%text_legend_origin%x, 3), ', ', &
                                                       real_str(g%text_legend_origin%y, 3), ', ', quote(g%text_legend_origin%units)
  nl=nl+1; write(lines(nl), amt)  'curve_legend_origin%x,y,units    = ', real_str(g%curve_legend_origin%x, 3), ', ', &
                                                       real_str(g%curve_legend_origin%y, 3), ', ', quote(g%curve_legend_origin%units)
  nl=nl+1; write(lines(nl), f3mt) 'curve_legend%line_length         = ', g%curve_legend%line_length
  nl=nl+1; write(lines(nl), f3mt) 'curve_legend%text_offset         = ', g%curve_legend%text_offset
  nl=nl+1; write(lines(nl), f3mt) 'curve_legend%row_spacing         = ', g%curve_legend%row_spacing
  nl=nl+1; write(lines(nl), lmt)  'curve_legend%draw_line           = ', g%curve_legend%draw_line
  nl=nl+1; write(lines(nl), lmt)  'curve_legend%draw_symbol         = ', g%curve_legend%draw_symbol
  nl=nl+1; write(lines(nl), lmt)  'curve_legend%draw_text           = ', g%curve_legend%draw_text

  if (g%type == 'floor_plan') then
    nl=nl+1; write(lines(nl), amt)  'floor_plan%view                  = ', quote(g%floor_plan%view)
    nl=nl+1; write(lines(nl), fmt)  'floor_plan%rotation              = ', g%floor_plan%rotation
    nl=nl+1; write(lines(nl), lmt)  'floor_plan%correct_distortion    = ', g%floor_plan%correct_distortion
    nl=nl+1; write(lines(nl), lmt)  'floor_plan%flip_label_side       = ', g%floor_plan%flip_label_side
    nl=nl+1; write(lines(nl), lmt)  'floor_plan%size_is_absolute      = ', g%floor_plan%size_is_absolute
    nl=nl+1; write(lines(nl), lmt)  'floor_plan%draw_building_wall    = ', g%floor_plan%draw_building_wall
    nl=nl+1; write(lines(nl), lmt)  'floor_plan%draw_only_first_pass  = ', g%floor_plan%draw_only_first_pass
    nl=nl+1; write(lines(nl), fmt)  'floor_plan%orbit_scale           = ', g%floor_plan%orbit_scale
    nl=nl+1; write(lines(nl), amt)  'floor_plan%orbit_color           = ', quote(g%floor_plan%orbit_color)
    nl=nl+1; write(lines(nl), amt)  'floor_plan%orbit_lattice         = ', quote(g%floor_plan%orbit_lattice)
    nl=nl+1; write(lines(nl), amt)  'floor_plan%orbit_pattern         = ', quote(g%floor_plan%orbit_pattern)
    nl=nl+1; write(lines(nl), imt)  'floor_plan%orbit_width           = ', g%floor_plan%orbit_width
  endif

  do i = 1, size(g%text_legend)
    if (g%text_legend(i) == '') cycle
    nl=nl+1; write(lines(nl), '(a, i0, 2a)') 'text_legend(', i, ')                = ', quote(g%text_legend(i))
  enddo

  nl=nl+1; write(lines(nl), amt)  'x%label                          = ', quote(g%x%label)
  nl=nl+1; write(lines(nl), rmt)  'x%max                            = ', g%x%max
  nl=nl+1; write(lines(nl), rmt)  'x%min                            = ', g%x%min
  nl=nl+1; write(lines(nl), imt)  'x%major_div                      = ', g%x%major_div
  nl=nl+1; write(lines(nl), imt)  'x%major_div_nominal              = ', g%x%major_div_nominal
  nl=nl+1; write(lines(nl), imt)  'x%places                         = ', g%x%places
  nl=nl+1; write(lines(nl), lmt)  'x%draw_label                     = ', g%x%draw_label
  nl=nl+1; write(lines(nl), lmt)  'x%draw_numbers                   = ', g%x%draw_numbers
  nl=nl+1; write(lines(nl), rmt)  'x%tick_max                       = ', g%x%tick_max
  nl=nl+1; write(lines(nl), rmt)  'x%tick_min                       = ', g%x%tick_min
  nl=nl+1; write(lines(nl), rmt)  'x%dtick                          = ', g%x%dtick
  nl=nl+1; write(lines(nl), rmt)  'x%eval_max                       = ', g%x%eval_max
  nl=nl+1; write(lines(nl), rmt)  'x%eval_min                       = ', g%x%eval_min

  nl=nl+1; write(lines(nl), lmt)  'y2_mirrors_y                     = ', g%y2_mirrors_y
  nl=nl+1; write(lines(nl), amt)  'y%label                          = ', quote(g%y%label)
  nl=nl+1; write(lines(nl), rmt)  'y%label_offset                   = ', g%y%label_offset
  nl=nl+1; write(lines(nl), rmt)  'y%max                            = ', g%y%max
  nl=nl+1; write(lines(nl), rmt)  'y%min                            = ', g%y%min
  nl=nl+1; write(lines(nl), imt)  'y%major_div                      = ', g%y%major_div
  nl=nl+1; write(lines(nl), imt)  'y%major_div_nominal              = ', g%y%major_div_nominal
  nl=nl+1; write(lines(nl), imt)  'y%places                         = ', g%y%places
  nl=nl+1; write(lines(nl), lmt)  'y%draw_label                     = ', g%y%draw_label
  nl=nl+1; write(lines(nl), lmt)  'y%draw_numbers                   = ', g%y%draw_numbers
  nl=nl+1; write(lines(nl), rmt)  'y%tick_max                       = ', g%y%tick_max
  nl=nl+1; write(lines(nl), rmt)  'y%tick_min                       = ', g%y%tick_min
  nl=nl+1; write(lines(nl), rmt)  'y%dtick                          = ', g%y%dtick

  nl=nl+1; write(lines(nl), amt)  'y2%label                         = ', quote(g%y2%label)
  nl=nl+1; write(lines(nl), rmt)  'y2%label_offset                  = ', g%y2%label_offset
  nl=nl+1; write(lines(nl), rmt)  'y2%max                           = ', g%y2%max
  nl=nl+1; write(lines(nl), rmt)  'y2%min                           = ', g%y2%min
  nl=nl+1; write(lines(nl), imt)  'y2%major_div                     = ', g%y2%major_div
  nl=nl+1; write(lines(nl), imt)  'y2%major_div_nominal             = ', g%y2%major_div_nominal
  nl=nl+1; write(lines(nl), imt)  'y2%places                        = ', g%y2%places
  nl=nl+1; write(lines(nl), lmt)  'y2%draw_label                    = ', g%y2%draw_label
  nl=nl+1; write(lines(nl), lmt)  'y2%draw_numbers                  = ', g%y2%draw_numbers
  nl=nl+1; write(lines(nl), lmt)  'limited                          = ', g%limited
  nl=nl+1; write(lines(nl), lmt)  'clip                             = ', g%clip
  nl=nl+1; write(lines(nl), lmt)  'draw_axes                        = ', g%draw_axes
  nl=nl+1; write(lines(nl), lmt)  'draw_curve_legend                = ', g%draw_curve_legend
  nl=nl+1; write(lines(nl), lmt)  'draw_grid                        = ', g%draw_grid
  nl=nl+1; write(lines(nl), lmt)  'draw_title                       = ', g%draw_title
  nl=nl+1; write(lines(nl), lmt)  'draw_only_good_user_data_or_vars = ', g%draw_only_good_user_data_or_vars
  nl=nl+1; write(lines(nl), lmt)  'allow_wrap_around                = ', g%allow_wrap_around


  if (allocated(g%curve)) then
    nl=nl+1; lines(nl) = '                         |          Line           |          Symbols'
    nl=nl+1; lines(nl) = 'Curves:                  |     Mean         RMS    |      Mean         RMS'
    do i = 1, size(g%curve)
      nl=nl+1; write(lines(nl), '(3x, a)') quote(g%curve(i)%name)
      if (allocated(g%curve(i)%x_line)) then
        call tao_curve_rms_calc(g%curve(i), 'LINE', rms, mean)
        write(lines(nl)(26:), '(2es12.4)') mean, rms
      endif
      if (allocated(g%curve(i)%x_symb)) then
        call tao_curve_rms_calc(g%curve(i), 'SYMBOL', rms, mean)
        write(lines(nl)(53:), '(2es12.4)') mean, rms
      endif
    enddo
  else
    nl=nl+1; lines(nl) = 'Curves: None associated'
  endif

!----------------------------------------------------------------------
! history

case ('history')

  show_labels = .true.
  n_print = 50
  show_all = .false.

  do 
    call tao_next_switch (what2, [character(8):: '-no_num', '-all', '-filed'], .true., switch, err)

    if (err) return
    if (switch == '') exit

    select case (switch)
    case ('-no_num');   show_labels = .false.
    case ('-all');      n_print = 9999
    case ('-filed');    show_all = .true.
    case default
      read (switch, *, iostat = ios) n_print
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'ERROR READING HISTORY NUMBER')
        return
      endif
    end select
  enddo

  !

  if (n_print < 1) return
  n_count = n_print - s%com%ix_history

  ! If commands from previous sessions wanted...
  if (n_count > 0 .and. show_all) then  
    iu = lunget()
    call fullfilename(s%global%history_file, file_name)
    open (iu, file = file_name, status = 'old', iostat = ios)
    ix1 = 0
    do
      if (ix1+1 >= size(lines)) call re_allocate (lines, 2*size(lines))
      read (iu, '(a)', iostat = ios) lines(ix1+1)
      if (ios /= 0) exit
      ix1 = ix1 + 1
    enddo
    close(1)

    ix0 = max(1, ix1-n_count+1)
    do i = ix0, ix1
      if (show_labels) then
        nl=nl+1; write (lines(nl), '(i5, 2a)') i-ix1, ': ', trim(lines(i))
      else
        nl=nl+1; write (lines(nl), '(a)') lines(i)
      endif
    enddo

  ! Here if commands from previous sessions not printed. 
  else 
    ix = max(1, s%com%ix_history - n_print + 1)

    do
      if (ix > s%com%ix_history) exit
      if (nl >= size(lines)) call re_allocate (lines, 2*size(lines))

      if (s%history(ix)%ix /= 0) then
        if (show_labels) then
          nl=nl+1; write (lines(nl), '(i5, 2a)') s%history(ix)%ix, ': ', trim(s%history(ix)%cmd)
        else
          nl=nl+1; write (lines(nl), '(a)') s%history(ix)%cmd
        endif
      endif

      ix = ix + 1
    enddo

    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Note: Use "-filed" switch to see commands from previous sessions stored in: ' // s%global%history_file
  endif

!----------------------------------------------------------------------
! hom

case ('hom')

  nl=nl+1; lines(nl) = &
        '       #        Freq         R/Q           Q   m  Polarization_Angle'
  do i = 1, ubound(lat%ele, 1)
    ele => lat%ele(i)
    if (ele%key /= lcavity$) cycle
    if (ele%slave_status == multipass_slave$) cycle
    nl=nl+1; write(lines(nl), '(a, i6)') ele%name, i
    do j = 1, size(ele%wake%lr%mode)
      lr_mode => ele%wake%lr%mode(j)
      angle_str = '-'
      if (lr_mode%polarized) write (angle_str, '(f9.4)') lr_mode%angle
      nl=nl+1; write(lines(nl), '(i8, 3es12.4, i4, a)') j, &
                  lr_mode%freq, lr_mode%R_over_Q, lr_mode%Q, lr_mode%m, angle_str
    enddo
    nl=nl+1; lines(nl) = ' '
  enddo
  nl=nl+1; lines(nl) = '       #        Freq         R/Q           Q   m  Polarization_Angle'

!----------------------------------------------------------------------
! Internal
! Used for debugging purposes. See also 'debug' which is for "temparary" code.

case ('internal')

  call tao_next_switch (what2, [character(16):: '-pipe', '-control'], .true., switch, err)
  select case (switch)

  ! Format: show -pipe_buffer
  ! This is useful for debugging the real and integer array passing which is used in the pipe interface.
  case ('-pipe')

    nl=nl+1; write (lines(nl), imt) 'N_real: ', tao_c_interface_com%n_real

    do i = 1, min(tao_c_interface_com%n_real, 5)
      nl=nl+1; write (lines(nl), '(a, i0, es12.4)') 'Real:  ', i, tao_c_interface_com%c_real(i)
    enddo

    if (tao_c_interface_com%n_real > 10) then
      nl=nl+1; lines(nl) = ' ... etc...'
    endif

    do i = max(tao_c_interface_com%n_real-5, 6), tao_c_interface_com%n_real
      nl=nl+1; write (lines(nl), '(a, i0, es12.4)') 'Real:  ', i, tao_c_interface_com%c_real(i)
    enddo

    !

    nl=nl+1; lines(nl) = ''
    nl=nl+1; write (lines(nl), imt) 'N_int:  ', tao_c_interface_com%n_int

    do i = 1, min(tao_c_interface_com%n_int, 5)
      nl=nl+1; write (lines(nl), '(a, i0, i12)') 'Int:  ', i, tao_c_interface_com%c_integer(i)
    enddo

    if (tao_c_interface_com%n_int > 10) then
      nl=nl+1; lines(nl) = ' ... etc...'
    endif

    do i = max(tao_c_interface_com%n_int-5, 6), tao_c_interface_com%n_int
      nl=nl+1; write (lines(nl), '(a, i0, i12)') 'Int:  ', i, tao_c_interface_com%c_integer(i)
    enddo

  ! Format: show -control <element-name>
  ! Lattice lord/slave control info.
  case ('-control')
    call tao_locate_elements (what2, -1, eles, err)
    if (err .or. size(eles) == 0) then
      nl=nl+1; lines(nl) = 'Cannot find lattice element: ' // what2
      return
    endif

    ele => eles(1)%ele
    nl=nl+1; lines(nl) = 'For element: (' // trim(ele_loc_name(ele)) // ')  ' // ele%name
    fmt = '(2x, a14, i6, i8, i6, 3x, a, t70, a15, a, es12.4, 3x, a)'

    if (ele%n_slave + ele%n_slave_field /= 0) then 
      nl=nl+1; lines(nl) = 'Slaves: Type       %ic  %cntrl  back   Slave                         Param          Expression'
      do i = 1, ele%n_slave + ele%n_slave_field
        slave => pointer_to_slave (ele, i, contl, all$, j, i_con, i_ic)
        if (i <= ele%n_slave) then
          nl=nl+1; write (lines(nl), fmt) control_name(ele%lord_status), i_ic, i_con, j, &
                      trim(slave%name) // ' ' // trim(ele_loc_name(slave, .true., '()')), &
                      contl%attribute, expression_stack_to_string(contl%stack)
        else
          nl=nl+1; write (lines(nl), fmt)  'Field_Overlap', i_ic, i_con, j, &
                                              ele_loc_name(slave, .true.), slave%name
        endif
      enddo
    endif

    if (ele%n_lord + ele%n_lord_field /= 0) then 
      nl=nl+1; lines(nl) = 'Lords:  Type       %ic  %cntrl  back   Lord                          Param          Attrib_Value  Expssn_val  Expression'
      do i = 1, ele%n_lord + ele%n_lord_field
        lord => pointer_to_lord (ele, i, contl, j, all$, i_con, i_ic)
        if (i <= ele%n_lord) then
          call pointer_to_attribute(ele, contl%attribute, .false., a_ptr, err_flag)
          if (associated(a_ptr%r)) then
            write (val_str, '(es12.4)') a_ptr%r
          else
            val_str = ' ----'
          endif
          nl=nl+1; write (lines(nl), fmt) control_name(lord%lord_status), i_ic, i_con, j, &
                      trim(lord%name) // ' ' // trim(ele_loc_name(lord, .true., '()')), &
                      contl%attribute, val_str(1:12), contl%value, expression_stack_to_string(contl%stack)

        else
          nl=nl+1; write (lines(nl), fmt)  &
                                     'Field_Overlap', i_ic, i_con, j, ele_loc_name(lord, .true.), lord%name
        endif
      enddo
    endif

  end select

!----------------------------------------------------------------------
! keys

case ('key_bindings')

  call tao_key_info_to_str (1, 1, size(s%key), str, header)
  nl=nl+1; lines(nl) = ' Ix  ' // header

  do i = 1, size(s%key)
    call tao_key_info_to_str (i, 1, size(s%key), str, header)
    nl=nl+1; write(lines(nl), '(i3, 2x, a)') i, str
  enddo

  ! Custom keys
  do i = 1, size(s%com%key)
    if (s%com%key(i)%name /= '') then
      nl=nl+1; write(lines(nl), '(a, 2x, a)') trim(s%com%key(i)%name), trim(s%com%key(i)%expanded_str)
    endif
  enddo

!----------------------------------------------------------------------
! lattice

case ('lattice')

  print_slaves = .true.
  print_super_slaves = .true.
  limited = .false.
  all_lat = .false.
  where = 'exit'
  by_s = .false.
  print_header_lines = .true.
  print_tail_lines = .true.
  replacement_for_blank = ''
  ix_branch = s%global%default_branch
  undef_str = '---'
  print_lords = maybe$
  what_to_show = 'standard'
  allocate (ix_remove(size(col)))
  n_remove = 0
  lat_type = model$
  n_attrib = 0
  attrib0 = ''
  undef_uses_column_format = .false.
  called_from_pipe_cmd = .false.
  print_rms = .false.

  column(:) = old_show_lat_column_struct()
  col(:)    = show_lat_column_struct()

  ! get command line switches

  do
    call tao_next_switch (what2, [character(32):: &
        '-branch', '-blank_replacement', '-lords', '-center', '-middle', &
        '-tracking_elements', '-0undef', '-beginning', '-pipe', &
        '-no_label_lines', '-no_tail_lines', '-custom', '-s', '-radiation_integrals', '-remove_line_if_zero', &
        '-base', '-design', '-floor_coords', '-orbit', '-attribute', '-all', '-no_slaves', '-energy', &
        '-spin', '-undef0', '-no_super_slaves', '-sum_radiation_integrals', '-python', '-universe', '-rms', &
        '-6d_radiation_integrals', '-ri_radiation_integrals'], .true., switch, err)
    if (err) return
    if (switch == '') exit
    select case (switch)

    case ('-0undef')
      undef_str = '  0'

    case ('-undef0')
      undef_str = '  0'
      undef_uses_column_format = .true.

    case ('-tracking_elements')
      print_lords = no$

    case ('-all')
      all_lat = .true. 

    case ('-attribute')
      what_to_show = 'attributes'
      n_attrib = n_attrib + 1
      call tao_next_word(what2, attrib_list(n_attrib))

    case ('-base')
      lat_type = base$

    case ('-blank_replacement')
      call tao_next_word(what2, replacement_for_blank)

    case ('-branch')
      call tao_next_word(what2, aname)
      branch => pointer_to_branch(aname, u%model%lat)
      if (.not. associated(branch)) then
        nl=1; write(lines(1), *) 'Bad branch name or index: ', aname
        return
      endif
      ix_branch = branch%ix_branch

    case ('-custom')
      what_to_show = 'custom'
      call tao_next_word(what2, file_name)
      iu = lunget()
      open (iu, file = file_name, status = 'old', iostat = ios)
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT OPEN FILE: ' // file_name
        return
      endif
      column(:) = old_show_lat_column_struct()
      col(:)    = show_lat_column_struct()
      read (iu, nml = custom_show_list, iostat = ios)
      close (iu)
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ "CUSTOM_SHOW_LIST" NAMELIST IN FILE: ' // file_name
        return
      endif

      if (column(1)%name /= '') then ! Old style
        do i = 1, size(column)
          col(i) = show_lat_column_struct(column(i)%name, column(i)%format, column(i)%label, column(i)%remove_line_if_zero, &
                                                          column(i)%scale_factor, column(i)%width)
        enddo
      else
        do i = 1, size(column)
          col(i) = setup_lat_column(col(i)%name, col(i)%format, col(i)%label, col(i)%remove_line_if_zero, col(i)%scale_factor)
        enddo
      endif

    case ('-design')
      lat_type = design$

    case ('-energy')
      what_to_show = 'energy'

    case ('-floor_coords')
      select case (what_to_show)
      case ('orbit', 'orbit:spin');   what_to_show = 'floor_orbit'
      case default;                   what_to_show = 'floor_coords'
      end select

    case ('-lords')
      print_lords = yes$

    case ('-beginning')
      where = 'beginning'

    case ('-center', '-middle')
      where = 'middle'

    case ('-no_label_lines')
      print_header_lines = .false.
      print_tail_lines = .false.

    case ('-no_slaves')
      print_slaves = .false.

    case ('-no_super_slaves')
      print_super_slaves = .false.

    case ('-no_tail_lines')
      print_tail_lines = .false.

    case ('-orbit')
      select case (what_to_show)
      case ('spin');          what_to_show = 'orbit:spin'
      case ('floor_coords');  what_to_show = 'floor_orbit'
      case default;           what_to_show = 'orbit'
      end select

    case ('-pipe', '-python')
      called_from_pipe_cmd = .true.
      print_tail_lines = .false.

    case ('-radiation_integrals')
      what_to_show = 'rad_int'
      all_lat = .true.  ! Will only print where radiation integrals is non-zero

    case ('-6d_radiation_integrals')
      what_to_show = '6d_rad_int'
      all_lat = .true.  ! Will only print where radiation integrals is non-zero

    case ('-ri_radiation_integrals')
      what_to_show = 'ri_rad_int'
      all_lat = .true.  ! Will only print where radiation integrals is non-zero

    case ('-rms')
      print_rms = .true.

    case ('-sum_radiation_integrals')
      what_to_show = 'sum_rad_int'
      all_lat = .true. ! Will only print where radiation integrals is non-zero

    case ('-remove_line_if_zero')
      call tao_next_word(what2, aname)
      n_remove = n_remove + 1
      read (aname, *, iostat = ios) ix_remove(n_remove)
      if (ios /= 0 .or. ix_remove(n_remove) < 1 .or. ix_remove(n_remove) > size(col)) then
        nl=1; lines(1) = 'CANNOT READ OR OUT-OF RANGE "-remove_line_if_zero" argument'
        return
      endif

    case ('-s')
      by_s = .true.

    case ('-spin')
      if (what_to_show == 'orbit') then
        what_to_show = 'orbit:spin'
      else
        what_to_show = 'spin'
      endif

    case ('-universe')
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios) ix
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ "-universe" ARGUMENT'
        return
      endif

      call point_to_this_uni(ix, ix_branch, u, lat, branch, model_branch, tao_branch, design_tao_branch)
      if (.not. associated(u)) return
      
      if (aname == '' .or. ios /= 0 .or. .not. associated(u)) then
        nl=1; lines(1) = 'CANNOT READ OR OUT-OF RANGE "-universe" ARGUMENT'
        return
      endif

    case default
      attrib0 = trim(attrib0) // ' ' // trim(switch)
    end select

  enddo

  ! Find elements to use

  branch => lat%branch(ix_branch)
  if (allocated (picked_ele)) deallocate (picked_ele)
  allocate (picked_ele(0:branch%n_ele_max))

  if (by_s) then
    ix2 = index(attrib0, ':')
    if (ix2 == 0) then
      nl=1; lines(1) = 'NO ":" FOUND FOR RANGE SELECTION'
      return
    endif
    read (attrib0(1:ix2-1), *, iostat = ios1) s1
    read (attrib0(ix2+1:), *, iostat = ios2) s2
    if (ios1 /= 0 .or. ios2 /= 0) then
      nl=1; lines(1) = 'ERROR READING RANGE SELECTION: ' // attrib0
      return
    endif

    picked_ele = .false.
    do ie = 1, branch%n_ele_track
      select case (where)
      case ('exit');      s_ele = branch%ele(ie)%s
      case ('middle');    s_ele = (branch%ele(ie)%s_start + branch%ele(ie)%s) / 2
      case ('beginning'); s_ele = branch%ele(ie)%s_start
      end select
      if (s_ele >= s1 .and. s_ele <= s2) picked_ele(ie) = .true.
    enddo

  elseif (attrib0 == '*' .or. all_lat) then
    picked_ele = .true.

  elseif (attrib0 /= '') then
    call tao_locate_elements (attrib0, u%ix_uni, eles, err, lat_type, &
                  ignore_blank = .true., above_ubound_is_err = .false., ix_branch = ix_branch)
    if (err .or. size(eles) == 0) return

    if (ix_branch /= eles(1)%ele%ix_branch) then
      ix_branch = eles(1)%ele%ix_branch
      branch => lat%branch(ix_branch)
      call re_allocate2 (picked_ele, 0, branch%n_ele_max)
    endif

    picked_ele = .false.
    do i = 1, size(eles)
      if (print_lords == yes$ .and. eles(i)%ele%lord_status == not_a_lord$) cycle
      picked_ele(eles(i)%ele%ix_ele) = .true.
    enddo

  elseif (index(what_to_show, 'rad_int') /= 0) then
    picked_ele = .true.

  else
    picked_ele = .true.
    if (count(picked_ele) > 300 .and. print_lords == maybe$) then
      picked_ele(201:) = .false.
      limited = .true.
    endif

  endif

  !

  tao_lat => tao_pointer_to_tao_lat(u, lat_type)
  lat => tao_lat%lat
  model_branch => u%model_branch(ix_branch)
  tao_branch => u%model%tao_branch(ix_branch)
  design_tao_branch => u%design%tao_branch(ix_branch)

  ! Construct columns if needed.

  select case (what_to_show)
  case ('attributes')
    col( 1)  = setup_lat_column('#',             'i7',    '', .false., 1.0_rp)
    col( 2)  = setup_lat_column('x',             '2x',    '', .false., 1.0_rp)
    col( 3)  = setup_lat_column('ele::#[name]',  'a0',    '', .false., 1.0_rp)
    col( 4)  = setup_lat_column('ele::#[key]',   'a16',   '', .false., 1.0_rp)
    col( 5)  = setup_lat_column('ele::#[s]',     'f10.3', '', .false., 1.0_rp)
    col( 6)  = setup_lat_column('ele::#[l]',     'f8.3',  '', .false., 1.0_rp)
    i0 = 6

    do i = 1, n_attrib
      attrib = attrib_list(i)
      sub_name = upcase(attrib)
      a_type = attribute_type(sub_name)

      select case (a_type)
      case (is_logical$)
        fmt = 'l12'
      case (is_real$)
        fmt = 'es12.4'
      case (is_integer$)
        fmt = 'i12'
      case (is_switch$)
        col(i0+i) = setup_lat_column('x', '2x', '', .false., 1.0_rp)
        i0 = i0 + 1
        fmt = 'a20'
      case (is_string$)
        col(i0+i) = setup_lat_column('x', '2x', '', .false., 1.0_rp)
        i0 = i0 + 1
        fmt = 'a20'
      case default
        fmt = 'es12.4'
      end select

      ix = index(attrib, '@')
      if (ix /= 0) then
        fmt = attrib(ix+1:)
        attrib = attrib(1:ix-1)
        nc = str_first_in_set(fmt, '0123456789')
      endif

      col(i0+i) = setup_lat_column(attrib, fmt, attrib, .false., 1.0_rp)
    enddo

  case ('energy')
    col( 1)  = setup_lat_column('#',                      'i7',       '', .false., 1.0_rp)
    col( 2)  = setup_lat_column('x',                      '2x',       '', .false., 1.0_rp)
    col( 3)  = setup_lat_column('ele::#[name]',           'a0',       '', .false., 1.0_rp)
    col( 4)  = setup_lat_column('ele::#[key]',            'a17',      '', .false., 1.0_rp)
    col( 5)  = setup_lat_column('ele::#[s]',              'f10.3',    '', .false., 1.0_rp)
    col( 6)  = setup_lat_column('ele::#[orbit_x]',        'es14.6',   '', .false., 1.0_rp)
    col( 7)  = setup_lat_column('ele::#[orbit_y]',        'es14.6',   '', .false., 1.0_rp)
    col( 8)  = setup_lat_column('ele::#[orbit_z]',        'es14.6',   '', .false., 1.0_rp)
    col( 9)  = setup_lat_column('ele::#[orbit_pz]',       'es14.6',   '', .false., 1.0_rp)
    col(10)  = setup_lat_column('ele::#[e_tot]',          'es14.6',   '', .false., 1.0_rp)
    col(11)  = setup_lat_column('ele::#[pc]',             'es14.6',   '', .false., 1.0_rp)

  case ('floor_coords')
    col( 1)  = setup_lat_column('#',                      'i7',       '', .false., 1.0_rp)
    col( 2)  = setup_lat_column('x',                      '2x',       '', .false., 1.0_rp)
    col( 3)  = setup_lat_column('ele::#[name]',           'a0',       '', .false., 1.0_rp)
    col( 4)  = setup_lat_column('ele::#[key]',            'a17',      '', .false., 1.0_rp)
    col( 5)  = setup_lat_column('ele::#[s]',              'f12.5',    '', .false., 1.0_rp)
    col( 6)  = setup_lat_column('ele::#[l]',              'f10.5',    '', .false., 1.0_rp)
    col( 7)  = setup_lat_column('ele::#[x_position]',     'f12.5',    'X',     .false., 1.0_rp)
    col( 8)  = setup_lat_column('ele::#[y_position]',     'f12.5',    'Y',     .false., 1.0_rp)
    col( 9)  = setup_lat_column('ele::#[z_position]',     'f12.5',    'Z',     .false., 1.0_rp)
    col(10)  = setup_lat_column('ele::#[theta_position]', 'f12.5',    'Theta', .false., 1.0_rp)
    col(11)  = setup_lat_column('ele::#[phi_position]',   'f12.5',    'Phi',   .false., 1.0_rp)
    col(12)  = setup_lat_column('ele::#[psi_position]',   'f12.5',    'Psi',   .false., 1.0_rp)

  case ('floor_orbit')
    col( 1)  = setup_lat_column('#',                         'i7',       '', .false., 1.0_rp)
    col( 2)  = setup_lat_column('x',                         '2x',       '', .false., 1.0_rp)
    col( 3)  = setup_lat_column('ele::#[name]',              'a0',       '', .false., 1.0_rp)
    col( 4)  = setup_lat_column('ele::#[key]',               'a17',      '', .false., 1.0_rp)
    col( 5)  = setup_lat_column('ele::#[s]',                 'f12.5',    '', .false., 1.0_rp)
    col( 6)  = setup_lat_column('ele::#[l]',                 'f10.5',    '', .false., 1.0_rp)
    col( 7)  = setup_lat_column('ele::#[floor_orbit.x]',     'f12.5',    'X_Orb',     .false., 1.0_rp)
    col( 8)  = setup_lat_column('ele::#[floor_orbit.y]',     'f12.5',    'Y_Orb',     .false., 1.0_rp)
    col( 9)  = setup_lat_column('ele::#[floor_orbit.z]',     'f12.5',    'Z_Orb',     .false., 1.0_rp)
    col(10)  = setup_lat_column('ele::#[floor_orbit.theta]', 'f12.5',    'Theta_Orb', .false., 1.0_rp)
    col(11)  = setup_lat_column('ele::#[floor_orbit.phi]',   'f12.5',    'Phi_Orb',   .false., 1.0_rp)
    col(12)  = setup_lat_column('ele::#[floor_orbit.psi]',   'f12.5',    'Psi_Orb',   .false., 1.0_rp)

  case ('orbit')
    col( 1)  = setup_lat_column('#',                      'i7',       '', .false., 1.0_rp)
    col( 2)  = setup_lat_column('x',                      '2x',       '', .false., 1.0_rp)
    col( 3)  = setup_lat_column('ele::#[name]',           'a0',       '', .false., 1.0_rp)
    col( 4)  = setup_lat_column('ele::#[key]',            'a17',      '', .false., 1.0_rp)
    col( 5)  = setup_lat_column('ele::#[s]',              'f13.6',    '', .false., 1.0_rp)
    col( 6)  = setup_lat_column('ele::#[orbit_x]',        'es14.6',   '', .false., 1.0_rp)
    col( 7)  = setup_lat_column('ele::#[orbit_px]',       'es14.6',   '', .false., 1.0_rp)
    col( 8)  = setup_lat_column('ele::#[orbit_y]',        'es14.6',   '', .false., 1.0_rp)
    col( 9)  = setup_lat_column('ele::#[orbit_py]',       'es14.6',   '', .false., 1.0_rp)
    col(10)  = setup_lat_column('ele::#[orbit_z]',        'es14.6',   '', .false., 1.0_rp)
    col(11)  = setup_lat_column('ele::#[orbit_pz]',       'es14.6',   '', .false., 1.0_rp)

  case ('orbit:spin')
    col( 1)  = setup_lat_column('#',                      'i7',       '', .false., 1.0_rp)
    col( 2)  = setup_lat_column('x',                      '2x',       '', .false., 1.0_rp)
    col( 3)  = setup_lat_column('ele::#[name]',           'a0',       '', .false., 1.0_rp)
    col( 4)  = setup_lat_column('ele::#[key]',            'a17',      '', .false., 1.0_rp)
    col( 5)  = setup_lat_column('ele::#[s]',              'f13.6',    '', .false., 1.0_rp)
    col( 6)  = setup_lat_column('ele::#[orbit_x]',        'es14.6',   '', .false., 1.0_rp)
    col( 7)  = setup_lat_column('ele::#[orbit_px]',       'es14.6',   '', .false., 1.0_rp)
    col( 8)  = setup_lat_column('ele::#[orbit_y]',        'es14.6',   '', .false., 1.0_rp)
    col( 9)  = setup_lat_column('ele::#[orbit_py]',       'es14.6',   '', .false., 1.0_rp)
    col(10)  = setup_lat_column('ele::#[orbit_z]',        'es14.6',   '', .false., 1.0_rp)
    col(11)  = setup_lat_column('ele::#[orbit_pz]',       'es14.6',   '', .false., 1.0_rp)
    col(12)  = setup_lat_column('ele::#[spin_x]',         'es14.6',   '', .false., 1.0_rp)
    col(13)  = setup_lat_column('ele::#[spin_y]',         'es14.6',   '', .false., 1.0_rp)
    col(14)  = setup_lat_column('ele::#[spin_z]',         'es14.6',   '', .false., 1.0_rp)

  case ('spin')
    col( 1)  = setup_lat_column('#',                        'i7',       '', .false., 1.0_rp)
    col( 2)  = setup_lat_column('x',                        '2x',       '', .false., 1.0_rp)
    col( 3)  = setup_lat_column('ele::#[name]',             'a0',       '', .false., 1.0_rp)
    col( 4)  = setup_lat_column('ele::#[key]',              'a17',      '', .false., 1.0_rp)
    col( 5)  = setup_lat_column('ele::#[s]',                'f13.6',    '', .false., 1.0_rp)
    col( 6)  = setup_lat_column('ele::#[spin_x]',           'es14.6',   '', .false., 1.0_rp)
    col( 7)  = setup_lat_column('ele::#[spin_y]',           'es14.6',   '', .false., 1.0_rp)
    col( 8)  = setup_lat_column('ele::#[spin_z]',           'es14.6',   '', .false., 1.0_rp)
    if (branch%param%geometry == closed$) then
      col( 9)  = setup_lat_column('ele::#[spin_dn_dpz.x]',    'es14.6', '', .false., 1.0_rp)
      col(10)  = setup_lat_column('ele::#[spin_dn_dpz.y]',    'es14.6', '', .false., 1.0_rp)
      col(11)  = setup_lat_column('ele::#[spin_dn_dpz.z]',    'es14.6', '', .false., 1.0_rp)
      col(12)  = setup_lat_column('ele::#[spin_dn_dpz.amp]',  'es14.6', '', .false., 1.0_rp)
    endif

  case ('6d_rad_int')
    col(1)  = setup_lat_column('#',                     'i7',       '', .false., 1.0_rp)
    col(2)  = setup_lat_column('x',                     '2x',       '', .false., 1.0_rp)
    col(3)  = setup_lat_column('ele::#[name]',          'a0',       '', .false., 1.0_rp)
    col(4)  = setup_lat_column('ele::#[key]',           'a17',      '', .false., 1.0_rp)
    col(5)  = setup_lat_column('ele::#[s]',             'f10.3',    '', .false., 1.0_rp)
    if (branch%param%geometry == open$) then
      col(6)  = setup_lat_column('lat::rad_int1_6d.i0[#]',     'es10.2',  '', .true., 1.0_rp)
      col(7)  = setup_lat_column('lat::rad_int1_6d.i1[#]',     'es10.2',  '', .true., 1.0_rp)
      col(8)  = setup_lat_column('lat::rad_int1_6d.i2_e4[#]',  'es10.2',  '', .false., 1.0_rp)
      col(9)  = setup_lat_column('lat::rad_int1_6d.i3_e7[#]',  'es10.2',  '', .false., 1.0_rp)
      col(10) = setup_lat_column('lat::rad_int1_6d.i5a_e6[#]', 'es10.2',  '', .false., 1.0_rp)
      col(11) = setup_lat_column('lat::rad_int1_6d.i5b_e6[#]', 'es10.2',  '', .false., 1.0_rp)
      col(12) = setup_lat_column('lat::emit.a[#]',         'es10.2',  '', .false., 1.0_rp)
      col(13) = setup_lat_column('lat::emit.b[#]',         'es10.2',  '', .false., 1.0_rp)
      col(14) = setup_lat_column('lat::sigma.pz[#]',       'es10.2',  '', .false., 1.0_rp)
    else
      col(6)  = setup_lat_column('lat::rad_int1_6d.i0[#]',     'es10.2',  '', .true., 1.0_rp)
      col(7)  = setup_lat_column('lat::rad_int1_6d.i1[#]',     'es10.2',  '', .true., 1.0_rp)
      col(8)  = setup_lat_column('lat::rad_int1_6d.i2[#]',     'es10.2',  '', .false., 1.0_rp)
      col(9)  = setup_lat_column('lat::rad_int1_6d.i3[#]',     'es10.2',  '', .false., 1.0_rp)
      col(10) = setup_lat_column('lat::rad_int1_6d.i4a[#]',    'es10.2',  '', .false., 1.0_rp)
      col(11) = setup_lat_column('lat::rad_int1_6d.i5a[#]',    'es10.2',  '', .false., 1.0_rp)
      col(12) = setup_lat_column('lat::rad_int1_6d.i4b[#]',    'es10.2',  '', .false., 1.0_rp)
      col(13) = setup_lat_column('lat::rad_int1_6d.i5b[#]',    'es10.2',  '', .false., 1.0_rp)
      col(14) = setup_lat_column('lat::rad_int1_6d.i6b[#]',    'es10.2',  '', .false., 1.0_rp)
    endif

  case ('ri_rad_int')
    col(1)  = setup_lat_column('#',                     'i7',       '', .false., 1.0_rp)
    col(2)  = setup_lat_column('x',                     '2x',       '', .false., 1.0_rp)
    col(3)  = setup_lat_column('ele::#[name]',          'a0',       '', .false., 1.0_rp)
    col(4)  = setup_lat_column('ele::#[key]',           'a17',      '', .false., 1.0_rp)
    col(5)  = setup_lat_column('ele::#[s]',             'f10.3',    '', .false., 1.0_rp)
    if (branch%param%geometry == open$) then
      col(6)  = setup_lat_column('lat::rad_int1_ri.i0[#]',     'es10.2',  '', .true., 1.0_rp)
      col(7)  = setup_lat_column('lat::rad_int1_ri.i1[#]',     'es10.2',  '', .true., 1.0_rp)
      col(8)  = setup_lat_column('lat::rad_int1_ri.i2_e4[#]',  'es10.2',  '', .false., 1.0_rp)
      col(9)  = setup_lat_column('lat::rad_int1_ri.i3_e7[#]',  'es10.2',  '', .false., 1.0_rp)
      col(10) = setup_lat_column('lat::rad_int1_ri.i5a_e6[#]', 'es10.2',  '', .false., 1.0_rp)
      col(11) = setup_lat_column('lat::rad_int1_ri.i5b_e6[#]', 'es10.2',  '', .false., 1.0_rp)
      col(12) = setup_lat_column('lat::emit.a[#]',         'es10.2',  '', .false., 1.0_rp)
      col(13) = setup_lat_column('lat::emit.b[#]',         'es10.2',  '', .false., 1.0_rp)
      col(14) = setup_lat_column('lat::sigma.pz[#]',       'es10.2',  '', .false., 1.0_rp)
    else
      col(6)  = setup_lat_column('lat::rad_int1_ri.i0[#]',     'es10.2',  '', .true., 1.0_rp)
      col(7)  = setup_lat_column('lat::rad_int1_ri.i1[#]',     'es10.2',  '', .true., 1.0_rp)
      col(8)  = setup_lat_column('lat::rad_int1_ri.i2[#]',     'es10.2',  '', .false., 1.0_rp)
      col(9)  = setup_lat_column('lat::rad_int1_ri.i3[#]',     'es10.2',  '', .false., 1.0_rp)
      col(10) = setup_lat_column('lat::rad_int1_ri.i4a[#]',    'es10.2',  '', .false., 1.0_rp)
      col(11) = setup_lat_column('lat::rad_int1_ri.i5a[#]',    'es10.2',  '', .false., 1.0_rp)
      col(12) = setup_lat_column('lat::rad_int1_ri.i4b[#]',    'es10.2',  '', .false., 1.0_rp)
      col(13) = setup_lat_column('lat::rad_int1_ri.i5b[#]',    'es10.2',  '', .false., 1.0_rp)
      col(14) = setup_lat_column('lat::rad_int1_ri.i6b[#]',    'es10.2',  '', .false., 1.0_rp)
    endif

  case ('rad_int', 'sum_rad_int')
    col(1)  = setup_lat_column('#',                     'i7',       '', .false., 1.0_rp)
    col(2)  = setup_lat_column('x',                     '2x',       '', .false., 1.0_rp)
    col(3)  = setup_lat_column('ele::#[name]',          'a0',       '', .false., 1.0_rp)
    col(4)  = setup_lat_column('ele::#[key]',           'a17',      '', .false., 1.0_rp)
    col(5)  = setup_lat_column('ele::#[s]',             'f10.3',    '', .false., 1.0_rp)
    if (branch%param%geometry == open$) then
      col(6)  = setup_lat_column('lat::rad_int.i0[#]',     'es10.2',  '', .true., 1.0_rp)
      col(7)  = setup_lat_column('lat::rad_int.i1[#]',     'es10.2',  '', .true., 1.0_rp)
      col(8)  = setup_lat_column('lat::rad_int.i2_e4[#]',  'es10.2',  '', .false., 1.0_rp)
      col(9)  = setup_lat_column('lat::rad_int.i3_e7[#]',  'es10.2',  '', .false., 1.0_rp)
      col(10) = setup_lat_column('lat::rad_int.i5a_e6[#]', 'es10.2',  '', .false., 1.0_rp)
      col(11) = setup_lat_column('lat::rad_int.i5b_e6[#]', 'es10.2',  '', .false., 1.0_rp)
      col(12) = setup_lat_column('lat::emit.a[#]',         'es10.2',  '', .false., 1.0_rp)
      col(13) = setup_lat_column('lat::emit.b[#]',         'es10.2',  '', .false., 1.0_rp)
      col(14) = setup_lat_column('lat::sigma.pz[#]',       'es10.2',  '', .false., 1.0_rp)
    else
      col(6)  = setup_lat_column('lat::rad_int.i0[#]',     'es10.2',  '', .true., 1.0_rp)
      col(7)  = setup_lat_column('lat::rad_int.i1[#]',     'es10.2',  '', .true., 1.0_rp)
      col(8)  = setup_lat_column('lat::rad_int.i2[#]',     'es10.2',  '', .false., 1.0_rp)
      col(9)  = setup_lat_column('lat::rad_int.i3[#]',     'es10.2',  '', .false., 1.0_rp)
      col(10) = setup_lat_column('lat::rad_int.i4a[#]',    'es10.2',  '', .false., 1.0_rp)
      col(11) = setup_lat_column('lat::rad_int.i5a[#]',    'es10.2',  '', .false., 1.0_rp)
      col(12) = setup_lat_column('lat::rad_int.i4b[#]',    'es10.2',  '', .false., 1.0_rp)
      col(13) = setup_lat_column('lat::rad_int.i5b[#]',    'es10.2',  '', .false., 1.0_rp)
      col(14) = setup_lat_column('lat::rad_int.i6b[#]',    'es10.2',  '', .false., 1.0_rp)
    endif

  case ('standard')
    if (branch%param%particle == photon$) then
      col( 1) = setup_lat_column('#',                   'i7',         '', .false., 1.0_rp)
      col( 2) = setup_lat_column('x',                   '2x',         '', .false., 1.0_rp)
      col( 3) = setup_lat_column('ele::#[name]',        'a0',         '', .false., 1.0_rp)
      col( 4) = setup_lat_column('ele::#[key]',         'a17',        '', .false., 1.0_rp)
      col( 5) = setup_lat_column('ele::#[s]',           'f10.3',      '', .false., 1.0_rp)
      col( 6) = setup_lat_column('ele::#[l]',           'f8.3',       '', .false., 1.0_rp)
      col( 7) = setup_lat_column('ele::#[orbit_x]',     '3pf10.5',    'x [mm]', .false., 1.0_rp)
      col( 8) = setup_lat_column('ele::#[orbit_px]',    '3pf10.5',    'px [mr]', .false., 1.0_rp)
      col( 9) = setup_lat_column('ele::#[orbit_y]',     '3pf10.5',    'y [mm]', .false., 1.0_rp)
      col(10) = setup_lat_column('ele::#[orbit_py]',    '3pf10.5',    'py [mr]', .false., 1.0_rp)
      col(11) = setup_lat_column('ele::#[energy] - ele::#[E_tot]', 'f10.4', 'dE [eV]', .false., 1.0_rp)
      col(12) = setup_lat_column('ele::#[intensity_x]', 'f8.4',       'I_x', .false., 1.0_rp)
      col(13) = setup_lat_column('ele::#[intensity_y]', 'f8.4',       'I_y', .false., 1.0_rp)
      col(14) = setup_lat_column('ele::#[phase_x]',     'f10.4',      'phase_x', .false., 1.0_rp)
      col(15) = setup_lat_column('ele::#[phase_y]',     'f10.4',      'phase_y', .false., 1.0_rp)
      col(16) = setup_lat_column('x',                   '2x',         '', .false., 1.0_rp)
      col(17) = setup_lat_column('ele::#[state]',       'a11',        'Track|State', .false., 1.0_rp)
    else
      col(1)  = setup_lat_column('#',                   'i7',         '', .false., 1.0_rp)
      col(2)  = setup_lat_column('x',                   '2x',         '', .false., 1.0_rp)
      col(3)  = setup_lat_column('ele::#[name]',        'a0',         '', .false., 1.0_rp)
      col(4)  = setup_lat_column('ele::#[key]',         'a17',        '', .false., 1.0_rp)
      col(5)  = setup_lat_column('ele::#[s]',           'f10.3',      '', .false., 1.0_rp)
      col(6)  = setup_lat_column('ele::#[l]',           'f8.3',       '', .false., 1.0_rp)
      col(7)  = setup_lat_column('ele::#[beta_a]',      'f8.2',       '', .false., 1.0_rp)
      col(8)  = setup_lat_column('ele::#[phi_a]',       'f8.3',       'phi_a|[2pi]', .false., 1/twopi)
      col(9)  = setup_lat_column('ele::#[eta_x]',       'f7.2',       '', .false., 1.0_rp)
      col(10) = setup_lat_column('ele::#[orbit_x]',     '3pf8.3',     'orbit|x [mm]', .false., 1.0_rp)
      col(11) = setup_lat_column('ele::#[beta_b]',      'f8.2',       '', .false., 1.0_rp)
      col(12) = setup_lat_column('ele::#[phi_b]',       'f8.3',       'phi_b|[2pi]', .false., 1/twopi)
      col(13) = setup_lat_column('ele::#[eta_y]',       'f7.2',       '', .false., 1.0_rp)
      col(14) = setup_lat_column('ele::#[orbit_y]',     '3pf8.3',     'orbit|y [mm]', .false., 1.0_rp)
      col(15) = setup_lat_column('x',                   '2x',         '', .false., 1.0_rp)
      col(16) = setup_lat_column('ele::#[state]',       'a11',        'Track|State', .false., 1.0_rp)
    endif

  end select

  ! remove_line_if_zero bookkeeping. Ignore space lines (name = 'x')

  do ix = 1, n_remove
    j = 0
    do i = 1, size(column)
      if (col(i)%name == 'x') cycle
      j = j + 1
      if (j == ix_remove(ix)) then
        col(i)%remove_line_if_zero = .true.
        exit
      endif
      if (i == size(column)) then
        nl=1; lines(1) = 'ARGUMENT FOR "-remove_line_if_zero" OUT OF RANGE!'
        return
      endif
    enddo
  enddo

  ! Need to compute radiation integrals?

  do i = 1, size(column)
    if (index(col(i)%name, 'rad_int') /= 0) then
      ix = ix_branch
      call radiation_integrals (u%model%lat, tao_branch%orbit, tao_branch%modes_ri, tao_branch%ix_rad_int_cache, ix, u%model%rad_int_by_ele_ri)
      call radiation_integrals (u%design%lat, design_tao_branch%orbit, design_tao_branch%modes_ri, &
                                                            design_tao_branch%ix_rad_int_cache, ix, u%design%rad_int_by_ele_ri)
      call radiation_integrals (u%base%lat, u%base%tao_branch(ix)%orbit, u%base%tao_branch(ix)%modes_ri, &
                                                  u%base%tao_branch(ix)%ix_rad_int_cache, ix, u%base%rad_int_by_ele_ri)
      exit
    endif
  enddo

  ! Compute some column info

  col_info = show_lat_column_info_struct()

  do i = 1, size(column)
    name = col(i)%name

    if (name(1:7) == 'ele::#[' .and. index(name, ']') /= 0) then
      ix = index(name, ']')-1
      col_info(i)%attrib_name = upcase(name(8:ix))
      col_info(i)%attrib_type = attribute_type(col_info(i)%attrib_name)
    else
      col_info(i)%attrib_name = upcase(name)
      col_info(i)%attrib_type = attribute_type(col_info(i)%attrib_name)
    endif
  enddo

  !

  if (print_lords == yes$) then
    picked_ele(0:branch%n_ele_track) = .false.
  elseif (print_lords == no$) then
    picked_ele(branch%n_ele_track+1:branch%n_ele_max) = .false.
  endif

  !

  if (.not. print_super_slaves) then
    do ie = 0, branch%n_ele_max
      ele => branch%ele(ie)
      if (ele%slave_status == super_slave$) picked_ele(ie) = .false.
    enddo
  endif

  if (.not. print_slaves) then
    do ie = 0, branch%n_ele_max
      ele => branch%ele(ie)
      if (ele%slave_status == super_slave$ .or. ele%slave_status == multipass_slave$) picked_ele(ie) = .false.
    enddo
  endif

  !

  if (called_from_pipe_cmd) then
    line1 = ''
    line2 = ''
    line3 = ''
    ix1 = 0

  else
    select case (where)
    case ('exit');      line1 = '# Values shown are for the Downstream End of each Element (Girder elements shown at ref point):'
    case ('middle');    line1 = '# Values shown are for the Center of each Element (Girder elements shown at ref point):'
    case ('beginning'); line1 = '# Values shown are for the Upstream of each Element (Girder elements shown at ref point):'
    end select

    if (size(lat%branch) > 1) line1 = '# Branch ' // int_str(branch%ix_branch) // '.' // line1(2:)
    if (size(s%u) > 1) line1 = '# Universe ' // int_str(u%ix_uni) // '.' // line1(2:)


    ix1 = 1
    line2 = "#"
    line3 = "#"
  endif

  ! Setup columns

  do i = 1, size(column)
    if (col(i)%name == '') cycle

    ! Use finer scale for s and l if needed.

    if (branch%param%total_length < 1.0_rp .and. what_to_show /= 'custom') then
      select case (col(i)%name)
      case ('ele::#[s]')
        col(i)%format = 'F13.6'
        col(i)%width = 13

      case ('ele::#[l]')
        col(i)%format = 'F11.6'
        col(i)%width = 11
      end select
    endif

    !

    col(i)%format = '(' // trim(upcase(col(i)%format)) // ')'

    if (col(i)%width == 0) then
      if (col(i)%name /= 'ele::#[name]') then
        call out_io (s_error$, r_name, &
            'WIDTH = 0 CAN ONLY BE USED WITH "ele::#[name]" TYPE COLUMNS')
        return
      endif
      col(i)%width = 5
      do ie = 0, branch%n_ele_max
        if (.not. picked_ele(ie)) cycle
        col(i)%width = max(col(i)%width, len_trim(branch%ele(ie)%name)+1)
      enddo
    endif

    ix2 = ix1 + col(i)%width

    if (col(i)%label == '') then
      name = col(i)%name
      if (index(name, 'ele::') /= 0) then
        i1 = index(name, '[')
        i2 = index(name, ']')
        name = name(i1+1:i2-1)
      elseif (index(name, 'beam::') /= 0) then
        ix = index(name, 'beam::')
        name = name(ix+6:)
        i2 = index(name, '[')
        name = name(1:i2-1)
      elseif (index(name, 'lat::') /= 0) then
        ix = index(name, 'lat::')
        name = name(ix+5:)
        i2 = index(name, '[')
        name = name(1:i2-1)
      elseif  (name == '#' .or. name == '#index') then
        call set_this_show_lat_header (line2, line3, 'Index', 'I', called_from_pipe_cmd, ix2-5)
        ix1 = ix2
        cycle    
      elseif  (name == '#branch') then
        call set_this_show_lat_header (line2, line3, 'Branch', 'I', called_from_pipe_cmd, ix2-5)
        ix1 = ix2
        cycle    
      elseif  (name == '#branch>>index') then
        call set_this_show_lat_header (line2, line3, 'Brnch>>Indx', 'I', called_from_pipe_cmd, ix2-5)
        ix1 = ix2
        cycle    
      elseif (name == 'x') then
        ix1 = ix2
        cycle    
      else
        name = ''
      endif

      ix = index(name, '.')
      if (ix == 0) ix = index(name, '_')
      n = len_trim(name)

      if (called_from_pipe_cmd) then
        call set_this_show_lat_header (line2, line3, name, col(i)%format, called_from_pipe_cmd)

      elseif (index(col(i)%format, 'A') /= 0) then
        line2(ix1:) = name

      elseif (ix == 0) then
        line2(ix2-n:) = name

      else
        if (ix2 - ix + 1 > 0) then
          line2(ix2-ix+1:) = name(1:ix-1)
        else
          line2 = name(1:ix-1)
        endif

        if (ix2 - n + ix > 0) then
          line3(ix2-n+ix:) = name(ix+1:)
        else
          line3 = name(ix+1:)
        endif
      endif

    else
      if (index(col(i)%name, '[') /= 0 .and. index(col(i)%label, '|') == 0) then
        i1 = index(col(i)%name, '[')
        i2 = index(col(i)%name, ']')
        name = upcase(col(i)%name(i1+1:i2-1))
        if (attribute_units(name) == 'rad') then
          col(i)%label = trim(col(i)%label) // '|[' // trim(phase_units_str) // ']'
          col(i)%scale_factor = phase_units
        endif
      endif

      name = col(i)%label

      ix = index(name, '|')

      if (called_from_pipe_cmd) then
        if (ix /= 0) name(ix:ix) = '_'
        call set_this_show_lat_header (line2, line3, name, col(i)%format, called_from_pipe_cmd)

      elseif (ix == 0) then
        if (index(col(i)%format, 'A') /= 0) then
          line2(ix1:) = name
        else
          j = len_trim(name)
          line2(ix2-j:) = name(1:j)
        endif

      else
        if (index(col(i)%format, 'A') /= 0) then
          line2(ix1:) = name(1:ix-1)
          line3(ix1:) = trim(name(ix+1:))
        else
          j = ix-1
          line2(ix2-j:) = name(1:j)
          name = name(ix+1:)
          j = len_trim(name)
          line3(ix2-j:) = name(1:j)
        endif
      endif
    endif

    ix1 = ix2
  enddo

  ! Collect lines

  if (print_header_lines) then
    if (called_from_pipe_cmd) then
      nl=nl+1; lines(nl) = line2
      nl=nl+1; lines(nl) = line3
    else
      nl=nl+1; lines(nl) = line1
      nl=nl+1; lines(nl) = line2
      if (line3 /= '#') then
        nl=nl+1; lines(nl) = line3
      endif
    endif
  endif

  !--------------------------------------------------------------------------------------------
  ! Loop over all rows

  s_last = real_garbage$
  ie0 = branch%n_ele_max

  row_loop: do ie = 0, branch%n_ele_max
    if (.not. picked_ele(ie)) cycle

    if (size(lines) < nl+100) call re_allocate (lines, 2*nl, .false.)

    ! Add separator line to distinguish lord vs slave elements

    if (ie0 <= branch%n_ele_track .and. ie > branch%n_ele_track .and. print_header_lines) then
      nl=nl+1; lines(nl) = 'Lord Elements:'
    endif
    ie0 = ie

    !---------------------------------------------------------------------------------------------
    ! Loop over all columns

    line = ''
    nc = 1
    ele => branch%ele(ie)
    n_zeros_found = 0
    n_remove = 0
    select case (where)
    case ('beginning');   s_now = ele%s_start
    case ('middle');      s_now = (ele%s_start + ele%s) / 2
    case default;         s_now = ele%s
    end select

    do i = 1, size(column)

      j = max(i-1, 1)
      if (i > 1 .and. col(j)%name /= '') then
        if (called_from_pipe_cmd) then
          if (col(j)%name /= 'x') then
            do i0 = nc, nc+col(j)%width
              if (line(i0:i0) /= ' ') exit
            enddo
            if (i0 > nc) line(nc:) = line(i0:)
            nc = len_trim(line) + 2
            line(nc-1:nc-1) = ';'
          endif
        else
          nc  = nc + col(i-1)%width
        endif
      endif

      name = col(i)%name
      if (name == '') cycle

      if (col(i)%remove_line_if_zero) n_remove = n_remove + 1

      if ((name(1:7) == 'ele::#[' .and. index(name, ']') /= 0) .or. &
              col_info(i)%attrib_type == is_string$ .or. col_info(i)%attrib_type == is_switch$) then
        sub_name = col_info(i)%attrib_name
        a_type = col_info(i)%attrib_type

        ! Note: a_type = real$ is handled later...
        select case (a_type)

        ! If recognized as a Bmad name.
        case (is_logical$, is_integer$, is_switch$)
          call pointer_to_attribute(ele, sub_name, .true., a_ptr, err, .false.)
          if (err) then
            write (line(nc:), col(i)%format, iostat = ios) replacement_for_blank
            cycle
          endif

          select case (a_type)
          case (is_logical$)
            if (associated(a_ptr%l)) then
              write (line(nc:), col(i)%format, iostat = ios) a_ptr%l
            else  ! Must be stored as a real
              write (line(nc:), col(i)%format, iostat = ios) is_true(a_ptr%r)
            endif
          case (is_real$)
            line(nc:) = '?????????' ! Should not be here. write (line(nc:), col(i)%format, iostat = ios) a_ptr%r
          case (is_integer$)
            if (associated(a_ptr%i)) then
              write (line(nc:), col(i)%format, iostat = ios) a_ptr%i
            else  ! Must be stored as a real
              write (line(nc:), col(i)%format, iostat = ios) nint(a_ptr%r)
            endif
          case (is_switch$)
            if (associated(a_ptr%r)) then
              r = a_ptr%r
            else
              r = a_ptr%i
            endif
            write (line(nc:), col(i)%format, iostat = ios) switch_attrib_value_name(sub_name, r, ele)
          end select

          if (line(nc:) == '') write (line(nc:), col(i)%format, iostat = ios) replacement_for_blank
          cycle

        case (is_string$)
          call string_attrib (sub_name, ele, line(nc:))
          if (line(nc:) == '') write (line(nc:), col(i)%format, iostat = ios) replacement_for_blank
          cycle
        end select
      endif

      ios = 0

      if (name == '#' .or. name == '#index') then
        if (ele%ix_branch /= ix_branch) then
          myname = ele_loc_name(ele, .true.)
          line(nc:) = adjustr(myname(1:col(i)%width))
        else
          write (line(nc:), col(i)%format, iostat = ios) ele%ix_ele
        endif

      elseif (name == '#branch') then
        write (line(nc:), col(i)%format, iostat = ios) ele%ix_branch

      elseif (name == 'ele::#[type]') then
        if (ele%type == '') then
        else
          write (line(nc:), col(i)%format, iostat = ios) ele%type
        endif

      elseif (name /= 'x') then
        write (myname, '(i0, a, i0)') ix_branch, '>>', ie
        call str_substitute (name, '#', trim(myname), ignore_escaped = .true.)
        call str_substitute (name, '\#', '#')
        ix = index(name, 'ele::')

        select case (where)
        case ('middle');    eval_pt = anchor_center$
        case ('beginning'); eval_pt = anchor_beginning$
        case default;       eval_pt = anchor_end$
        end select

        call tao_evaluate_expression (name, 1, .false., value, err, .false., &
                      dflt_component = tao_lat_type_name(lat_type), dflt_source = 'ele', dflt_ele = ele, &
                      dflt_uni = u%ix_uni, dflt_eval_point = eval_pt)

        if (err .or. .not. allocated(value) .or. size(value) /= 1) then
          if (col(i)%remove_line_if_zero) n_zeros_found = n_zeros_found + 1
          if (undef_uses_column_format .and. index(col(i)%format, 'A') == 0) then
            if (index(col(i)%format, 'I') /= 0) then
              write (line(nc:), col(i)%format, iostat = ios) 0
            else
              write (line(nc:), col(i)%format, iostat = ios) 0.0_rp
            endif
          else
            n = len(undef_str)
            k = min(n, col(i)%width - 1)
            j = nc + col(i)%width - k
            line(j:) = undef_str(n-k+1:n)
          endif

        elseif (col(i)%name == 'ele::#[state]') then
          write (line(nc:), col(i)%format, iostat = ios) coord_state_name(nint(value(1)))

        elseif (index(col(i)%format, 'L') /= 0) then
          if (value(1) == 0) then
            write (line(nc:), col(i)%format, iostat = ios) .false.
          else
            write (line(nc:), col(i)%format, iostat = ios) .true.
          endif

        elseif (index(col(i)%format, 'I') /= 0) then
          write (line(nc:), col(i)%format, iostat = ios) nint(value(1))
          if (col(i)%remove_line_if_zero .and. nint(value(1)) == 0) n_zeros_found = n_zeros_found + 1

        else
          r = value(1) * col(i)%scale_factor
          call write_real (line(nc:), col(i)%format, r)
          if (col(i)%remove_line_if_zero .and. value(1) == 0) then
            n_zeros_found = n_zeros_found + 1
          else
            col_info(i)%indent = nc
            col_info(i)%n_sum = col_info(i)%n_sum + 1
            col_info(i)%sum   = col_info(i)%sum + r
            col_info(i)%sum2  = col_info(i)%sum2 + r*r
            if (s_last /= real_garbage$ .and. ie <= branch%n_ele_track) then
              ds  = s_now - s_last
              z = (r + col_info(i)%val_last) / 2
              col_info(i)%int_s    = col_info(i)%int_s + ds
              col_info(i)%int_sum  = col_info(i)%int_sum + ds * z
              col_info(i)%int_sum2 = col_info(i)%int_sum2 + ds * z*z
              col_info(i)%val_last = r
            endif
          endif
        endif
      endif

      if (ios /= 0) then
        lines(1) = 'WIDTH TOO SMALL FOR NUMBER OR BAD FORMAT: ' // col(i)%format
        lines(2) = 'FOR DISPLAYING: ' // col(i)%name
        nl = 2
        return
      endif

    enddo  ! column loop

    if (n_remove > 0 .and. n_zeros_found == n_remove) cycle
    if (called_from_pipe_cmd) line(nc-1:nc-1) = ' '  ! Remove final ';'
    s_last = s_now

    nl=nl+1; lines(nl) = line

  enddo row_loop

  !-----------------------------------------------------------------------------------------------

  if (print_tail_lines) then
    nl=nl+1; lines(nl) = line2
    if (line3 /= '#') then
      nl=nl+1; lines(nl) = line3
    endif
    nl=nl+1; lines(nl) = line1
  endif

  if (limited .and. print_tail_lines) then
    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), '(a, i0)') &
          'NOTE: Since no range given, the number of elements shown is first 200 of ', branch%n_ele_track
  endif

  ! RMS

  if (print_rms) then
    lines(nl+1) = 'Average:'
    lines(nl+2) = 'RMS:'
    lines(nl+3) = 'N_count:'
    lines(nl+4) = 'Integrated Average:'
    lines(nl+5) = 'Integrated RMS:'

    do i = 1, size(column)
      n = col_info(i)%n_sum
      if (n == 0) cycle
      nc = col_info(i)%indent

      mean = col_info(i)%sum / n
      rms = sqrt(max(0.0_rp, col_info(i)%sum2 / n - mean**2))
      write (lines(nl+1)(nc:), col(i)%format) mean 
      write (lines(nl+2)(nc:), col(i)%format) rms

      myname = int_str(n)
      nc = len_trim(lines(nl+2)) - len_trim(myname) + 1
      lines(nl+3)(nc:) = myname    ! Right justified

      ds = col_info(i)%int_s
      if (ds /= 0) then
        mean = col_info(i)%int_sum / ds
        rms = sqrt(max(0.0_rp, col_info(i)%int_sum2 / ds - mean**2))
        write (lines(nl+4)(nc:), col(i)%format) mean 
        write (lines(nl+5)(nc:), col(i)%format) rms
      endif
    enddo

    nl=nl+5
  endif

  deallocate(picked_ele)

!----------------------------------------------------------------------
! merit (old style: top10)

case ('merit', 'top10')

  if (what2 == '') then
    call tao_show_constraints (0, 'TOP10')
    call tao_top10_merit_categories_print (0)
  elseif (index('-derivative', trim(what2)) == 1) then 
    call tao_top10_derivative_print ()
  elseif (index('-merit_only', trim(what2)) == 1) then
    nl=nl+1; write (lines(nl), '(a, es14.6)') 'Merit = ', tao_merit()
  else
    nl=1; lines(1) = 'UNKNOWN SWITCH: ' // what2
    return
  endif

!----------------------------------------------------------------------
! normal_form

case ('normal_form')

  nl=nl+1; lines(nl) = 'The "show normal_form" command is currently being reworked.'
  nl=nl+1; lines(nl) = 'Use "show chromaticity" for the chromaticity, momentum compaction and phase slip Taylor maps.'

!  select case(attrib0(1:5))
!    case ('dhdj ')
!      call type_taylors (bmad_nf%dhdj, max_order = n_order, lines = lines, n_lines = nl, clean = .true.)
!    case ('A    ')
!      call type_taylors (bmad_nf%A, max_order = n_order, lines = lines, n_lines = nl, clean = .true.)
!    case ('A_inv')
!      call type_taylors (bmad_nf%A_inv, max_order = n_order, lines = lines, n_lines = nl, clean = .true.)
!    case ('M    ')
!      call type_taylors (bmad_nf%M, max_order = n_order, lines = lines, n_lines = nl, clean = .true.)
!    case ('F    ')
!      call type_complex_taylors (bmad_nf%F, max_order = n_order, lines = lines, n_lines = nl, clean = .true.)
!    case ('L    ')
!      call type_complex_taylors (bmad_nf%L, max_order = n_order, lines = lines, n_lines = nl, clean = .true.)      
!    case ('h    ')
!      do i=1, size(bmad_nf%h)
!        write(lines(i),'(A,"   (",f20.10,", ",f20.10,")   ")') bmad_nf%h(i)%c, bmad_nf%h(i)%c_val
!      enddo
!      nl = size(bmad_nf%h)
!    case default 
!      nl=nl+1; lines(nl) = 'bad normal_form map: '//trim(attrib0)
!      nl=nl+1; lines(nl) = 'Must be one of: M A A_inv dhdj F L'
!  end select

!----------------------------------------------------------------------
! optimizer

case ('optimizer')

  do i = lbound(s%u, 1), ubound(s%u, 1)
    u => s%u(i)
    nl=nl+1; lines(nl) = 'Data Used:'
    if (size(s%u) > 1) then
      nl=nl+1; write(lines(nl), '(2(a, i0))') 'Universe: ', i, '  of: ', ubound(s%u, 1)
    endif
    do j = 1, u%n_d2_data_used
      if (u%d2_data(j)%name == ' ') cycle
      call tao_data_show_use (u%d2_data(j), lines, nl)
    enddo
  enddo

  nl=nl+1; lines(nl) = 'Variables Used:'
  do j = 1, s%n_v1_var_used
    if (s%v1_var(j)%name == ' ') cycle
    call tao_var_show_use (s%v1_var(j), lines, nl)
  enddo

  nl=nl+1; lines(nl) = ' '
  nl=nl+1; write(lines(nl), amt) 'optimizer:        ', quote(s%global%optimizer)
  call show_opt
  call out_io (s_blank$, r_name, lines(1:nl))
  nl = 0

!----------------------------------------------------------------------
! particle

case ('orbit')

  call tao_locate_elements (word1, u%ix_uni, eles, err)
  if (err) return
  ele => eles(1)%ele
  orb = u%model%tao_branch(ele%ix_branch)%orbit(ele%ix_ele)
  do i = 1, 6
    nl=nl+1; write(lines(nl), rmt) coord_name(i), orb%vec(i)
  enddo

!----------------------------------------------------------------------
! particle

case ('particle')

  nb = s%global%bunch_to_plot
  ix_branch = s%global%default_branch
  show_all = .false.
  show_lost = .false.
  ele_name = ''
  ix_ele = -1
  ix_p = 1

  do
    call tao_next_switch (what2, [character(16):: '-element', '-particle', '-bunch', '-lost', '-all'], .true., switch, err)
    if (err) return

    select case (switch)
    case ('') 
      exit

    case ('-lost') 
      show_lost = .true.

    case ('-all')
      show_all = .true.

    case ('-element')
      call tao_next_word(what2, ele_name)

      if (ele_name /= 'init') then
        call tao_pick_universe (ele_name, ele_name, picked_uni, err, ix_u)
        if (err) return
        call tao_locate_elements (ele_name, ix_u, eles, err)
        if (err) return
        ix_ele = eles(1)%ele%ix_ele
        ix_branch = eles(1)%ele%ix_branch
      endif

    case ('-particle')
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios) ix_p
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ PARTICLE INDEX!'
        return
      endif

    case ('-bunch')
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios) nb
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ BUNCH INDEX!'
        return
      endif

    case default
      call out_io (s_error$, r_name, 'EXTRA STUFF ON LINE: ' // attrib0)
      return
    end select

  enddo

  model_branch => u%model_branch(ix_branch)
  branch => u%model%lat%branch(ix_branch)

  ! show lost

  if (show_lost) then
    bunch => u%model_branch(ix_branch)%ele(lat%n_ele_track)%beam%bunch(nb)
    nl=nl+1; write(lines(nl), *) 'Bunch:', nb
    nl=nl+1; lines(nl) = 'Particles lost at:'
    nl=nl+1; lines(nl) = '    Ix Ix_Ele  Ele_Name '
    do i = 1, size(bunch%particle)
      if (bunch%particle(i)%state == alive$) cycle
      if (nl == size(lines)) call re_allocate (lines, nl+100, .false.)
      ie = bunch%particle(i)%ix_ele
      nl=nl+1; write(lines(nl), '(i6, i7, 2x, a)') i, ie, lat%ele(ie)%name
    enddo
    result_id = 'particle:lost'
    return
  endif

  ! check

  if (.not. allocated(model_branch%beam%beam_at_start%bunch)) then
    call out_io (s_error$, r_name, 'NO BEAM TRACKING HAS BEEN DONE.')
    return
  endif

  if (nb < 1 .or. nb > size(model_branch%beam%beam_at_start%bunch)) then
    call out_io (s_error$, r_name, 'BUNCH INDEX OUT OF RANGE: \i0\ ', i_array = [ nb ])
    return
  endif

  !

  if (ix_ele == -1) then
    bunch => model_branch%beam%beam_at_start%bunch(nb)

  else
    if (.not. allocated(model_branch%ele(ix_ele)%beam%bunch)) then
      call out_io (s_error$, r_name, 'BUNCH NOT ASSOCIATED WITH THIS ELEMENT.')
      return
    endif
    bunch => model_branch%ele(ix_ele)%beam%bunch(nb)
  endif

  ! show all

  if (show_all) then
    if (ix_ele == -1) then
      nl=nl+1; write(lines(nl), *) 'Initial Distribution'
    else
      nl=nl+1; write(lines(nl), *) 'Element:', ix_ele, '  ', branch%ele(ix_ele)%name
    endif
    nl=nl+1; write(lines(nl), '(a, 6(11x, a, 2x), (8x, a), (13x, a), 7x, a)') '  Ix', ' x', 'px', ' y', 'py', ' z', 'pz', 'dTime', 'State'
    do i = 1, size(bunch%particle)
      if (nl == size(lines)) call re_allocate (lines, nl+100, .false.)
      vec_in = bunch%particle(i)%vec
      if (bunch%particle(i)%beta == 0) then
        dt = 0
      else
        dt = -vec_in(5) / (c_light * bunch%particle(i)%beta)
      endif
      nl=nl+1; write(lines(nl), '(i6, 7es15.7, 2x, a)') i, (vec_in(j), j = 1, 6), dt, adjustr(coord_state_name(bunch%particle(i)%state))
    enddo
    result_id = 'particle:lost'
    return
  endif

  !

  if (ix_p < 1 .or. ix_p > size(bunch%particle)) then
    call out_io (s_error$, r_name, 'PARTICLE INDEX OUT OF RANGE: \i0\ ', i_array = [ ix_p ])
    return
  endif

  if (ix_ele == -1) then
    nl=nl+1; write(lines(nl), imt) 'Initial Distribtion'
  else
    nl=nl+1; write(lines(nl), imt) 'At lattice element:', ix_ele
  endif
  nl=nl+1; write(lines(nl), imt) 'Bunch:       ', nb
  nl=nl+1; write(lines(nl), imt) 'Particle:    ', ix_p
  nl=nl+1; write(lines(nl), lmt) 'Is Alive?    ', bunch%particle(ix_p)%state == alive$
  if (u%model%lat%branch(ix_branch)%param%particle == photon$) then
    nl=nl+1; write(lines(nl), rmt) 'Intensity_x: ', bunch%particle(ix_p)%field(1)**2
    nl=nl+1; write(lines(nl), rmt) 'Intensity_y: ', bunch%particle(ix_p)%field(2)**2
  else
    nl=nl+1; write(lines(nl), rmt) 'Charge:      ', bunch%particle(ix_p)%charge
  endif
  nl=nl+1; write(lines(nl), lmt) 'Coords: '
  nl=nl+1; write(lines(nl), '(a, 6es13.5)') '  ', bunch%particle(ix_p)%vec

!----------------------------------------------------------------------
! plot

case ('plot')

  ! Look for switches

  what = ''
  attrib0 = ''

  do
    call tao_next_switch (what2, [character(16) :: '-floor_plan', '-lat_layout', '-templates', &
                                     '-global', '-regions', '-plot_page', '-page'], .true., switch, err)
    if (err) return
    select case (switch)
    case ('') 
      exit

    case default
      if (switch(1:1) == '-') Then
        what = switch
      else 
        if (attrib0 /= '') then
          call out_io (s_error$, r_name, 'EXTRA STUFF ON LINE: ' // attrib0)
          return
        endif
        attrib0 = switch
      endif
    end select
  enddo

  ! Find particular plot

  if (attrib0 /= '') then

    call tao_find_plots (err, attrib0, 'BOTH', plot, print_flag = .false., only_visible = .false.)
    if (err) then
      nl = 1; lines(nl) = 'CANNOT FIND PLOT WITH NAME: ' // attrib0
      return
    endif

    found = .false.
    do j = 1, size(plot)
      p => plot(j)%p
      if (p%name == '') cycle

      if (found) then
        nl=nl+1; lines(nl) = ''
        nl=nl+1; lines(nl) = '======================================='
        nl=nl+1; lines(nl) = ''
      endif
      found = .true.

      nl=nl+1; lines(nl) = 'Plot:  ' // p%name
      if (associated(p%r)) then
        nl=nl+1; lines(nl) = 'Region:  ' // trim(p%r%name)
        nl=nl+1; write (lines(nl), lmt)  'Visible                = ', p%r%visible
        nl=nl+1; write (lines(nl), f3mt) 'Location [x1,x2,y1,y2] = ', p%r%location
      endif
      nl=nl+1; write(lines(nl), amt) 'x_axis_type          = ', quote(p%x_axis_type)
      nl=nl+1; write(lines(nl), lmt) 'autoscale_x          = ', p%autoscale_x
      nl=nl+1; write(lines(nl), lmt) 'autoscale_y          = ', p%autoscale_y
      nl=nl+1; write(lines(nl), lmt) 'autoscale_gang_x     = ', p%autoscale_gang_x
      nl=nl+1; write(lines(nl), lmt) 'autoscale_gang_y     = ', p%autoscale_gang_y
      nl=nl+1; write(lines(nl), imt) 'n_curve_pts          = ', p%n_curve_pts

      nl=nl+1; lines(nl) = 'Graphs:'
      do i = 1, size(p%graph)
        nl=nl+1; write(lines(nl), amt) '   ', quote(p%graph(i)%name)
      enddo
    enddo

    if (.not. found) then
      nl=1; lines(1) = 'This is not a name of a plot'
      result_id = 'ERROR'
    endif

    return
  endif

  ! Floor plan info

  select case (what)
  case ('-floor_plan', '-lat_layout')

    if (what == '-floor_plan') then
      shapes => s%plot_page%floor_plan%ele_shape
    else
      shapes => s%plot_page%lat_layout%ele_shape
    endif

    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Element Shapes:'
    nl=nl+1; lines(nl) = '                                                                                                                  Line_'
    nl=nl+1; lines(nl) = '                  Ele_ID                              Shape           Color           Size  Label    Draw  Multi  Width  Offset'
    nl=nl+1; lines(nl) = '                  ------------------------------      ----------      -------         ----  ------  -----  -----  -----  ------'

    do i = 1, size(shapes)
      shape => shapes(i)
      if (shape%ele_id == '') cycle
      nl=nl+1; write(lines(nl), '(a, i0, a, t19, a, t55, a, t71, a, t83, f8.1, 2x, a, t102, l5, l6, i7, 2x, f8.1)') &
                'ele_shape(', i, ') = ', quote(shape%ele_id), quote(shape%shape), quote(shape%color), &
                shape%size, quote(shape%label), shape%draw, shape%multi, shape%line_width, shape%offset
    enddo

    do i = 1, size(s%plot_page%pattern)
      pattern => s%plot_page%pattern(i)
      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = 'Shape Pattern Name: ' // trim(pattern%name)
      nl=nl+1; write (lines(nl), '(a, i0)') 'Line Width = ', pattern%line%width
      nl=nl+1; lines(nl) = '            s         x'
      do j = 1, size(pattern%pt)
        nl=nl+1; write (lines(nl), '(5x, 2f10.5)') pattern%pt(j)%s, pattern%pt(j)%y
      enddo
    enddo


    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Shape label text size scales as:'
    if (what == '-floor_plan') then
      nl=nl+1; lines(nl) = '    plot_page%text_height * plot_page%legend_text_scale * plot_page%floor_plan_text_scale'
    else
      nl=nl+1; lines(nl) = '    plot_page%text_height * plot_page%legend_text_scale * plot_page%lat_layout_text_scale'
    endif

    result_id = 'plot:floor_plan'
    return

  ! Plot_page parameters.

  case ('-page', '-plot_page', '-global')

    nl=nl+1; lines(nl) = 'plot_page parameters:'
    nl=nl+1; write(lines(nl), amt)  '  %title%string                  = ', quote(s%plot_page%title%string)
    nl=nl+1; write(lines(nl), '(a, 2f10.3, 2x, a, 2x, a)') &
                                    '  %title%x, y, units, justify    = ', s%plot_page%title%x, s%plot_page%title%y, &
                                                          quote(s%plot_page%title%units), quote(s%plot_page%title%justify)
    nl=nl+1; write(lines(nl), amt)  '  %subtitle%string               = ', quote(s%plot_page%subtitle%string)
    nl=nl+1; write(lines(nl), '(a, 2f10.3, 2x, a, 2x, a)') &
                                    '  %subtitle%x, y, units, justify = ', s%plot_page%title%x, s%plot_page%title%y, &
                                                          quote(s%plot_page%title%units), quote(s%plot_page%title%justify)
    nl=nl+1; write(lines(nl), imt)  '  %size                          = ', nint(s%plot_page%size)
    nl=nl+1; write(lines(nl), imt)  '  %n_curve_pts                   = ', s%plot_page%n_curve_pts
    nl=nl+1; write(lines(nl), '(a, 4(f0.3, 2x), a)') &
                                    '  %border                        = ', s%plot_page%border%x1, s%plot_page%border%x2, &
                                                     s%plot_page%border%y1, s%plot_page%border%y2, s%plot_page%border%units
    nl=nl+1; write(lines(nl), f3mt) '  %text_height                   = ', s%plot_page%text_height 
    nl=nl+1; write(lines(nl), f3mt) '  %main_title_text_scale         = ', s%plot_page%main_title_text_scale 
    nl=nl+1; write(lines(nl), f3mt) '  %graph_title_text_scale        = ', s%plot_page%graph_title_text_scale 
    nl=nl+1; write(lines(nl), f3mt) '  %axis_number_text_scale        = ', s%plot_page%axis_number_text_scale 
    nl=nl+1; write(lines(nl), f3mt) '  %axis_label_text_scale         = ', s%plot_page%axis_label_text_scale 
    nl=nl+1; write(lines(nl), f3mt) '  %key_table_text_scale          = ', s%plot_page%key_table_text_scale 
    nl=nl+1; write(lines(nl), '(a, f0.3, 3x, a)') &
                                    '  %legend_text_scale             = ', s%plot_page%legend_text_scale, &
                                                                            '! For legends, plot_page, and lat_layout' 
    nl=nl+1; write(lines(nl), f3mt) '  %floor_plan_shape_scale        = ', s%plot_page%floor_plan_shape_scale 
    nl=nl+1; write(lines(nl), f3mt) '  %floor_plan_text_scale         = ', s%plot_page%floor_plan_text_scale 
    nl=nl+1; write(lines(nl), f3mt) '  %lat_layout_shape_scale        = ', s%plot_page%lat_layout_shape_scale 
    nl=nl+1; write(lines(nl), f3mt) '  %lat_layout_text_scale         = ', s%plot_page%lat_layout_text_scale 
    nl=nl+1; write(lines(nl), lmt)  '  %delete_overlapping_plots      = ', s%plot_page%delete_overlapping_plots
    nl=nl+1; write(lines(nl), lmt)  '  %draw_graph_title_suffix       = ', s%plot_page%draw_graph_title_suffix


    result_id = 'plot:global'
    return 

  ! Template plot parameters

  case ('-templates')

    nl=nl+1; lines(nl) = ''
    nl=nl+1; lines(nl) = 'Templates:'
    nl=nl+1; lines(nl) = '   Plot Name                               Description        '
    nl=nl+1; lines(nl) = '   ----------------------------            -------------------'
    do i = 1, size(s%plot_page%template)
      p => s%plot_page%template(i)
      if (p%name == '') cycle
      if (p%name == 'scratch') cycle
      if (.not. p%list_with_show_plot_command) cycle
      nl=nl+1; write(lines(nl), '(3x, 2a)') p%name, trim(p%description)
    enddo

    result_id = 'plot:template'
    return 

  ! Plot regions

  case ('-regions', '')

    nl=nl+1; lines(nl) = ''
    if (s%global%plot_on) then
      nl=nl+1; lines(nl) = '                                               Location on Page'
      nl=nl+1; lines(nl) = 'Plot Region         <-->  Plot                 x1    x2    y1    y2     Visible'
      nl=nl+1; lines(nl) = '-----------               -----------------------------------------------------'
    else
      nl=nl+1; lines(nl) = 'Plot Region         <-->  Plot                 Visible'
      nl=nl+1; lines(nl) = '-----------               ----------------------------'
    endif
    found = .false.
    do i = 1, size(s%plot_page%region)
      region => s%plot_page%region(i)
      if (region%name == '') cycle
      if (.not. region%list_with_show_plot_command .and. region%plot%name == '' .and. what == '') then
        found = .true.
        cycle
      endif
      if (.not. s%global%plot_on) then
        nl=nl+1; write(lines(nl), '(a20, a, a21, l1)') region%name, '<-->  ', region%plot%name, region%visible
      elseif (region%plot%name == '') then
        nl=nl+1; write(lines(nl), '(a20, a, 18x, 4f6.2)') region%name, '<-->  ', region%location
      else
        nl=nl+1; write(lines(nl), '(a20, a, a18, 4f6.2, l10)') region%name, '<-->  ', region%plot%name, region%location, region%visible
      endif
    enddo

    if (found .and. what == '') then
      nl=nl+1; write(lines(nl), '(a)') '[Etc... In the interest of brevity, other plotless regions are not listed. Use "show plot -regions" for the entire list.]'
    endif

    result_id = 'plot:'
    return

  end select

!----------------------------------------------------------------------
! ptc

case ('ptc')

  what_to_show = ''

  do
    call tao_next_switch (what2, [character(24):: '-emittance'], .true., switch, err)
    if (err) return

    select case (switch)
    case ('')
      exit
    case ('-emittance')
      what_to_show = switch
    end select
  enddo

  if (.not. associated(branch%ptc%m_t_layout)) then
    call out_io (s_info$, r_name, 'Note: Creating PTC layout (equivalent to "ptc init").')
    call lat_to_ptc_layout (lat)
  endif


  call ptc_emit_calc (branch%ele(0), mode_ptc, sig_mat, orb)
  nl=nl+1; lines(nl) = '  Mode      Emit        Tune'
  nl=nl+1; write(lines(nl), '(1x, a, 2x, es15.7, f12.6)') 'A', mode_ptc%a%emittance, mode_ptc%a%tune
  nl=nl+1; write(lines(nl), '(1x, a, 2x, es15.7, f12.6)') 'B', mode_ptc%b%emittance, mode_ptc%b%tune
  nl=nl+1; write(lines(nl), '(1x, a, 2x, es15.7, f12.6)') 'C', mode_ptc%z%emittance, mode_ptc%z%tune
  nl=nl+1; lines(nl) = ''
  nl=nl+1; write(lines(nl), '(1x, a, 6es15.7)') 'Starting orbit:', orb%vec

!----------------------------------------------------------------------
! radiation_integrals

case ('radiation_integrals')

  ix_u = s%global%default_universe
  b_name = ''

  do
    call tao_next_switch (what2, [character(20):: '-branch'], .true., switch, err)
    if (err) return
    select case (switch)
    case ('')
      exit

    case ('-branch')
      call tao_next_word(what2, b_name)

    case default
      read (switch, *, iostat = ios) ix_u
      if (ios /= 0) then
        nl=1; lines(1) = 'BAD UNIVERSE NUMBER'
        return
      endif
      if (ix_u < lbound(s%u, 1) .or. ix_u > ubound(s%u, 1)) then
        nl=1; lines(1) = 'UNIVERSE NUMBER OUT OF RANGE'
        return
      endif
    end select
  enddo

  !

  u => s%u(ix_u)
  lat => u%model%lat

  if (b_name /= '') then
    branch => pointer_to_branch(b_name, u%model%lat)
    if (.not. associated(branch)) then
      nl=1; write(lines(1), *) 'Bad branch name or index: ', b_name
      return
    endif
    ix_branch = branch%ix_branch
  endif

  branch => lat%branch(ix_branch)
  model_branch => u%model_branch(ix_branch)
  tao_branch => u%model%tao_branch(ix_branch)

  design_lat => u%design%lat
  design_branch => design_lat%branch(ix_branch)
  design_tao_branch => u%design%tao_branch(ix_branch)

  nl = 0
  nl=nl+1; write(lines(nl), '(2(a, i0))') 'Universe: ', ix_u,      '  Of: ', ubound(s%u, 1)
  nl=nl+1; write(lines(nl), '(3(a, i0))') 'Branch:   ', ix_branch, '  In range: [0, ', ubound(lat%branch,1), ']'

  call radiation_integrals (lat, tao_branch%orbit, tao_branch%modes_ri, tao_branch%ix_rad_int_cache, ix_branch)

  mode_d => design_tao_branch%modes_ri
  mode_m => tao_branch%modes_ri
  l_lat = branch%param%total_length
  mode_d%momentum_compaction = mode_d%synch_int(1)/l_lat
  mode_m%momentum_compaction = mode_m%synch_int(1)/l_lat

  nt = branch%n_ele_track
  time1 = branch%ele(nt)%ref_time
  gamma = branch%ele(0)%value(e_tot$) / mass_of(branch%ele(0)%ref_species)

  fmt  = '(1x, a16, 2es13.5, 2x, 2es13.5, 2x, a)'
  fmt2 = '(1x, a16,  2f13.6, 2x,  2f13.6, 2x, a)'
  fmt3 = '(1x, a16,         28x, 2es13.5, 2x, a)'

  nl=nl+1; lines(nl) = ''
  nl=nl+1; write(lines(nl), '(23x, a)') '         X            |              Y'
  nl=nl+1; write(lines(nl), '(23x, a)') '  Model       Design  |       Model       Design'

  if (branch%param%geometry == closed$) then
    nl=nl+1; write(lines(nl), fmt2) 'J_damp', mode_m%a%j_damp, mode_d%a%j_damp, mode_m%b%j_damp, mode_d%b%j_damp, '! Damping Partition #'
    nl=nl+1; write(lines(nl), fmt) 'Emittance', mode_m%a%emittance, mode_d%a%emittance, mode_m%b%emittance, mode_d%b%emittance, '! Unnormalized'
    nl=nl+1; write(lines(nl), '(a43, 2x, 2es13.5)') 'Emit (photon vert opening angle ignored)', &
                                                                      mode_m%b%emittance_no_vert, mode_d%b%emittance_no_vert
    if (mode_m%b%alpha_damp /= 0) then
      nl=nl+1; write(lines(nl), fmt) 'Alpha_damp', mode_m%a%alpha_damp, mode_d%a%alpha_damp, mode_m%b%alpha_damp, mode_d%b%alpha_damp, '! Damping per turn'
      nl=nl+1; write(lines(nl), fmt) 'Damping_time', time1/mode_m%a%alpha_damp, time1/mode_d%a%alpha_damp, time1/mode_m%b%alpha_damp, time1/mode_d%b%alpha_damp, '! Sec'
    endif

    nl=nl+1; write(lines(nl), fmt) 'I4', mode_m%a%synch_int(4), mode_d%a%synch_int(4), mode_m%b%synch_int(4), mode_d%b%synch_int(4), '! Radiation Integral'
    nl=nl+1; write(lines(nl), fmt) 'I5', mode_m%a%synch_int(5), mode_d%a%synch_int(5), mode_m%b%synch_int(5), mode_d%b%synch_int(5), '! Radiation Integral'
    nl=nl+1; write(lines(nl), fmt3) 'I6b/gamma^2', mode_m%b%synch_int(6) / gamma**2, mode_d%b%synch_int(6) / gamma**2, '! Radiation Integral'
  else
    nl=nl+1; write(lines(nl), fmt) 'Final Emittance', mode_m%lin%a_emittance_end, mode_d%lin%a_emittance_end, mode_m%lin%b_emittance_end, mode_d%lin%b_emittance_end, '! Meters'
    nl=nl+1; write(lines(nl), fmt) 'I5*gamma^6', mode_m%lin%i5a_e6, mode_d%lin%i5a_e6, mode_m%lin%i5b_e6, mode_d%lin%i5b_e6, '! Linac Radiation Integral'
  endif

  fmt  = '(1x, a16, 2es13.5, 3x, a)'
  fmt2 = '(1x, a16, 2f13.7, 3x, a)'

  nl=nl+1; lines(nl) = ''
  nl=nl+1; write(lines(nl), '(23x, a)') '  Model       Design'

  if (branch%param%geometry == closed$) then
    if (mode_m%z%alpha_damp /= 0) then
      nl=nl+1; write(lines(nl), fmt) 'Sig_E/E:', mode_m%sigE_E, mode_d%sigE_E
      nl=nl+1; write(lines(nl), fmt) 'Sig_z:  ', mode_m%sig_z, mode_d%sig_z, '! Only calculated when RF is on'
      nl=nl+1; write(lines(nl), fmt) 'Energy Loss:', mode_m%e_loss, mode_d%e_loss, '! Energy_Loss (eV / Turn)'
      nl=nl+1; write(lines(nl), fmt) 'J_damp:', mode_m%z%j_damp, mode_d%z%j_damp, '! Longitudinal Damping Partition #'
      nl=nl+1; write(lines(nl), fmt) 'Alpha_damp:', mode_m%z%alpha_damp, mode_d%z%alpha_damp, '! Longitudinal Damping per turn'
      nl=nl+1; write(lines(nl), fmt) 'damp_time:', time1/mode_m%z%alpha_damp, time1/mode_d%z%alpha_damp, '! Longitudinal Damping time (sec)'
    endif

    nl=nl+1; write(lines(nl), fmt) 'Alpha_p:', mode_m%momentum_compaction, mode_d%momentum_compaction, '! Momentum Compaction'
    nl=nl+1; write(lines(nl), fmt) 'Eta_p:', mode_m%momentum_compaction - 1.0_rp/gamma**2, mode_d%momentum_compaction - 1.0_rp/gamma**2, '! Slip factor'
    if (mode_m%momentum_compaction > 0) then
      nl=nl+1; write(lines(nl), fmt) 'gamma_trans:', sqrt(1.0_rp/mode_m%momentum_compaction), sqrt(1.0_rp/mode_d%momentum_compaction), '! Gamma at transition'
    else
      nl=nl+1; write(lines(nl), amt) 'gamma_trans:  No transition since momentum compaction is negative.'
    endif

    nl=nl+1; write(lines(nl), fmt) 'I0:', mode_m%synch_int(0), mode_d%synch_int(0), '! Radiation Integral'
    nl=nl+1; write(lines(nl), fmt) 'I1:', mode_m%synch_int(1), mode_d%synch_int(1), '! Radiation Integral'
    nl=nl+1; write(lines(nl), fmt) 'I2:', mode_m%synch_int(2), mode_d%synch_int(2), '! Radiation Integral'
    nl=nl+1; write(lines(nl), fmt) 'I3:', mode_m%synch_int(3), mode_d%synch_int(3), '! Radiation Integral'
  else
    nl=nl+1; write(lines(nl), fmt) 'I1:', mode_m%synch_int(1), mode_d%synch_int(1), '! Radiation Integral'
    nl=nl+1; write(lines(nl), fmt) 'I2*gamma^4', mode_m%lin%i2_e4, mode_d%lin%i2_e4, '! Linac Radiation Integral'
    nl=nl+1; write(lines(nl), fmt) 'I3*gamma^7', mode_m%lin%i3_e7, mode_d%lin%i3_e7, '! Linac Radiation Integral'
  endif


!----------------------------------------------------------------------
! rampers

case ('rampers')

  show_energy = .false.

  do 
    call tao_next_switch (what2, [character(16):: '-universe', '-energy_show'], .true., switch, err)
    if (err) return
    if (switch == '') exit

    select case (switch)

    case ('-universe')
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios) ix
      u => tao_pointer_to_universe(ix)
      if (aname == '' .or. ios /= 0 .or. .not. associated(u)) then
        nl=1; lines(1) = 'CANNOT READ OR OUT-OF RANGE "-universe" ARGUMENT'
        return
      endif

    case ('-energy_show')
      show_energy = .true.
    end select
  enddo

  lat => u%model%lat

  do ie = lat%n_ele_track+1, lat%n_ele_max
    lord => lat%ele(ie)
    if (lord%key /= ramper$) cycle
    cntr => lord%control
    if (nl + size(cntr%ramp) + 20 > size(lines)) call re_allocate (lines, 2*(nl + size(cntr%ramp)), .false.)
    nl=nl+1; write(lines(nl), '(2a)') 'Ramper: ', ele_full_name(lord)
    do ix = 1, size(cntr%ramp)
      rmp => cntr%ramp(ix)
      if (allocated(rmp%stack)) then
        call split_expression_string (expression_stack_to_string(rmp%stack), 100, 5, alloc_lines)
      else  ! Spline
        call split_expression_string (knots_to_string(lord%control%x_knot, rmp%y_knot), 100, 5, alloc_lines, '),')
      endif
      value_here = ramper_value(lord, rmp, err_flag)
      nl=nl+1; write (lines(nl), '(4x, 4a, t34, es20.12, 4x, a)') trim(rmp%slave_name), '[', trim(rmp%attribute), ']', value_here, trim(alloc_lines(1))
      if (size(alloc_lines) > 1) then
        lines(nl) = trim(lines(nl)) // ' ...'
      endif
    enddo
    nl=nl+1; lines(nl) = ''
  enddo

  n_count = 0

  do ib = 0, ubound(lat%branch, 1)
    do ie = 0, lat%branch(ib)%n_ele_max
      ele => lat%branch(ib)%ele(ie)
      if (ele%n_lord_ramper == 0) cycle
      if (ele%n_lord_ramper == 1 .and. .not. show_energy) then
        rml => ele%control%ramper_lord(1)
        rmp => lat%ele(rml%ix_ele)%control%ramp(rml%ix_con)
        if (rmp%attribute == 'P0C' .or. rmp%attribute == 'E_TOT') then
          n_count = n_count + 1
          cycle
        endif
      endif

      if (nl + ele%n_lord_ramper + 20 > size(lines)) call re_allocate (lines, 2*nl, .false.)
      nl=nl+1; write(lines(nl), '(2a, 4x, a)') 'Slave: ', ele_full_name(ele), key_name(ele%key)
      cntr => ele%control
      do ix = 1, ele%n_lord_ramper
        rml => cntr%ramper_lord(ix)
        lord => lat%ele(rml%ix_ele)
        rmp => lord%control%ramp(rml%ix_con)
        if (rmp%attribute == 'P0C' .or. rmp%attribute == 'E_TOT') then
          n_count = n_count + 1
          if (.not. show_energy) cycle
        endif
        nl=nl+1; write(lines(nl), '(4x, a, t40, a)') ele_full_name(lord), rmp%attribute
      enddo
      nl=nl+1; lines(nl) = ''
    enddo
  enddo

  if (.not. show_energy) then
    nl=nl+1; lines(nl) = 'Number of elements with energy ramping: ' // int_str(n_count)
  endif

!----------------------------------------------------------------------
! spin

case ('spin')

  what_to_show = 'standard'
  show_mat = .false.
  sm => datum%spin_map
  sm%axis_input = spin_axis_struct()
  ele_ref_name = ''
  ele_ref => null()
  flip = .false.
  excite_zero = ''
  veto = ''

  do
    call tao_next_switch (what2, [character(24):: '-element', '-n_axis', '-l_axis', &
                            '-g_map', '-flip_n_axis', '-x_zero', '-y_zero', &
                            '-z_zero', '-ignore_kinetic', '-isf', '-spin_tune'], .true., switch, err)
    if (err) return

    select case (switch)
    case ('')
      exit
    case ('-element')
      what_to_show = 'element'
      call tao_next_word(what2, ele_ref_name)
      call upcase_string(ele_ref_name)

      ele_name = ele_ref_name
      if (what2(1:1) /= '-' .and. what2(1:1) /= ' ') then
        call tao_next_word(what2, ele_name)
        call upcase_string(ele_name)
      endif
    case ('-flip_n_axis')
      flip = .true.
    case ('-isf')
      what_to_show = 'isf'
    case ('-spin_tune')
      what_to_show = 'spin_tune'
    case ('-n_axis')
      read (what2, *, iostat = ios) sm%axis_input%n0
      if (ios /= 0) then
        nl=nl+1; lines(nl) = 'CANNOT PARSE N-AXIS: ' // what2
        return
      endif
      show_mat = .true.
      call word_read(what2, ' ,', word1, ix, delim, delim_found, what2) ! Strip Axis from what2
      call word_read(what2, ' ,', word1, ix, delim, delim_found, what2)
      call word_read(what2, ' ,', word1, ix, delim, delim_found, what2)
    case ('-l_axis')
      read (what2, *, iostat = ios) sm%axis_input%l
      if (ios /= 0) then
        nl=nl+1; lines(nl) = 'CANNOT PARSE L-AXIS: ' // what2
        return
      endif
      show_mat = .true.
      call word_read(what2, ' ,', word1, ix, delim, delim_found, what2) ! Strip Axis from what2
      call word_read(what2, ' ,', word1, ix, delim, delim_found, what2)
      call word_read(what2, ' ,', word1, ix, delim, delim_found, what2)
    case ('-q_map')
      nl=nl+1; lines(nl) = 'Note: "-q_map" now no longer needed or used.'
    case ('-g_map')
      show_mat = .true.
    case ('-ignore_kinetic')
      call word_read(what2, ' ', veto, ix, delim, delim_found, what2)
    case ('-x_zero')
      call word_read(what2, ' ', excite_zero(1), ix, delim, delim_found, what2)
    case ('-y_zero')
      call word_read(what2, ' ', excite_zero(2), ix, delim, delim_found, what2)
    case ('-z_zero')
      call word_read(what2, ' ', excite_zero(3), ix, delim, delim_found, what2)
    end select
  enddo

  show_q = (.not. show_mat)

  !

  if (.not. bmad_com%spin_tracking_on) call tao_spin_tracking_turn_on

  !

  if (what_to_show == 'isf') then
    if (branch%param%geometry == open$) then
      nl=nl+1; lines(nl) = 'No ISF for an open lattice!'
      return
    endif

    tao_branch%spin_map_valid = .false.
    if (.not. u%calc%one_turn_map) call tao_ptc_normal_form (.true., u%model, branch%ix_branch)

    ptc_nf  => tao_branch%ptc_normal_form
    do i = 1, 3
      ctaylor(i) = ptc_nf%isf%x(i)
    enddo

    call type_complex_taylors(ctaylor, out_type = 'SPIN')
    return
  endif
  
  if (what_to_show == 'spin_tune') then
    if (branch%param%geometry == open$) then
      nl=nl+1; lines(nl) = 'No spin tune for an open lattice!'
      return
    endif

    tao_branch%spin_map_valid = .false.
    if (.not. u%calc%one_turn_map) call tao_ptc_normal_form (.true., u%model, branch%ix_branch)

    ptc_nf  => tao_branch%ptc_normal_form
    ctaylor1 = ptc_nf%spin_tune

    call type_complex_taylors([ctaylor1], out_type = 'NONE')
    return   
  endif

  ! what_to_show = standard

  r = anomalous_moment_of(branch%param%particle) * branch%ele(1)%value(e_tot$) / mass_of(branch%param%particle)
  nl=nl+1; lines(nl) = 'a_anomalous_moment * gamma = ' // real_str(r, 6)
  nl=nl+1; lines(nl) = 'E_tot = ' // real_str(branch%ele(1)%value(e_tot$), 6)
  qs = branch%param%spin_tune/twopi
  nl=nl+1; write (lines(nl), '(a, f12.6, a)') 'Spin Tune', qs, 'rad/2pi'

  if (what_to_show == 'standard') then
    ele => branch%ele(0)

    nl=nl+1; write(lines(nl), lmt) 'bmad_com%spin_tracking_on                = ', bmad_com%spin_tracking_on

    if (branch%param%geometry == open$) then
      orb = tao_branch%orbit(0)
      nl=nl+1; lines(nl) = ''
      nl=nl+1; write(lines(nl), '(a, 3f12.8)') 'Beginning spin:', orb%spin
    else
      tao_branch%spin_map_valid = .false.
      call tao_spin_polarization_calc (branch, tao_branch, excite_zero, veto)

      !

      nl=nl+1; lines(nl) = ''
      nl=nl+1; write (lines(nl), '(a, es18.7)') 'spin_tune: ', tao_branch%spin%tune / twopi
      if (tao_branch%spin%valid) then
        if (tao_branch%spin%pol_rate_bks == 0) then
          nl=nl+1; lines(nl) = 'No bends or other radiation producing lattice elements detected!'
        else
          r = c_light * tao_branch%orbit(0)%beta / branch%param%total_length
          nl=nl+1; write(lines(nl), '(a, f12.8)')  'Polarization Limit ST:                   ', tao_branch%spin%pol_limit_st
          nl=nl+1; write(lines(nl), '(a, f12.8)')  'Polarization Limit DK:                   ', tao_branch%spin%pol_limit_dk
          nl=nl+1; write(lines(nl), '(a, 3f12.8)') 'Polarization Limits DK (a,b,c-modes):    ', tao_branch%spin%pol_limit_dk_partial
          nl=nl+1; write(lines(nl), '(a, 3f12.8)') 'Polarization Limits DK (bc,ac,ab-modes): ', tao_branch%spin%pol_limit_dk_partial2

          if (tao_branch%spin%pol_rate_bks == 0) then
            nl=nl+1; write(lines(nl), '(a)')    'Polarization Time BKS (minutes, turns): plarization rate is zero!'
          else
            x = 1.0_rp / tao_branch%spin%pol_rate_bks
            nl=nl+1; write(lines(nl), '(a, a12, es12.4)')    'Polarization Time BKS (minutes, turns): ', real_str(x/60.0_rp, 3), r*x
          endif

          if (1.0_rp / tao_branch%spin%depol_rate == 0) then
            nl=nl+1; write(lines(nl), '(a)')    'Depolarization Time (minutes, turns):   Depolarization rate is zero!'
          else
            x = 1.0_rp / tao_branch%spin%depol_rate
            nl=nl+1; write(lines(nl), '(a, a12, es12.4)')    'Depolarization Time (minutes, turns):   ', real_str(x/60.0_rp, 3), r*x
          endif

          if (tao_branch%spin%depol_rate_partial(1) == 0) then
            nl=nl+1; write(lines(nl), '(a)')   'Depolarization Time (a-mode) (minutes, turns): Depolarization rate is zero!'
          else
            x = 1 / tao_branch%spin%depol_rate_partial(1)
            nl=nl+1; write(lines(nl), '(a, a12, 3es12.4)')   'Depolarization Time (a-mode) (minutes, turns):', real_str(x/60.0_rp, 3), r*x
          endif

          if (tao_branch%spin%depol_rate_partial(2) == 0) then
            nl=nl+1; write(lines(nl), '(a)')   'Depolarization Time (b-mode) (minutes, turns): Depolarization rate is zero!'
          else
            x = 1 / tao_branch%spin%depol_rate_partial(2)
            nl=nl+1; write(lines(nl), '(a, a12, 3es12.4)')   'Depolarization Time (b-mode) (minutes, turns):', real_str(x/60.0_rp, 3), r*x
          endif

          if (tao_branch%spin%depol_rate_partial(3) == 0) then
            nl=nl+1; write(lines(nl), '(a)')   'Depolarization Time (c-mode) (minutes, turns): Depolarization rate is zero!'
          else
            x = 1 / tao_branch%spin%depol_rate_partial(3)
            nl=nl+1; write(lines(nl), '(a, a12, 3es12.4)')   'Depolarization Time (c-mode) (minutes, turns):', real_str(x/60.0_rp, 3), r*x
          endif

          if (tao_branch%spin%depol_rate_partial2(1) == 0) then
            nl=nl+1; write(lines(nl), '(a)')   'Depolarization Time (b&c modes) (minutes, turns): Depolarization rate is zero!'
          else
            x = 1 / tao_branch%spin%depol_rate_partial2(1)
            nl=nl+1; write(lines(nl), '(a, a12, 3es12.4)')   'Depolarization Time (b&c modes) (minutes, turns):', real_str(x/60.0_rp, 3), r*x
          endif

          if (tao_branch%spin%depol_rate_partial2(2) == 0) then
            nl=nl+1; write(lines(nl), '(a)')   'Depolarization Time (a&c modes) (minutes, turns): Depolarization rate is zero!'
          else
            x = 1 / tao_branch%spin%depol_rate_partial2(2)
            nl=nl+1; write(lines(nl), '(a, a12, 3es12.4)')   'Depolarization Time (a&c modes) (minutes, turns):', real_str(x/60.0_rp, 3), r*x
          endif

          if (tao_branch%spin%depol_rate_partial2(3) == 0) then
            nl=nl+1; write(lines(nl), '(a)')   'Depolarization Time (a&b modes) (minutes, turns): Depolarization rate is zero!'
          else
            x = 1 / tao_branch%spin%depol_rate_partial2(3)
            nl=nl+1; write(lines(nl), '(a, a12, 3es12.4)')   'Depolarization Time (a&b modes) (minutes, turns):', real_str(x/60.0_rp, 3), r*x
          endif

          nl=nl+1; write(lines(nl), '(a, a14)')            'Integral g^3 * b_hat * n_0:         ', real_str(tao_branch%spin%integral_bn, 5)
          nl=nl+1; write(lines(nl), '(a, a14)')            'Integral g^3 * b_hat * dn/ddelta:   ', real_str(tao_branch%spin%integral_bdn, 5)
          nl=nl+1; write(lines(nl), '(a, a14)')            'Integral g^3 (1 - 2(n * s_hat)/9):  ', real_str(tao_branch%spin%integral_1ns, 5)
          nl=nl+1; write(lines(nl), '(a, a14)')            'Integral g^3 * 11 (dn/ddelta)^2 / 9:', real_str(tao_branch%spin%integral_dn2, 5)
        endif

      else
        nl=nl+1; lines(nl) = 'Polarization calc not valid.'
      endif
    endif

    j = 0
    do iu = 1, size(s%u)
      u => s%u(iu)
      do i = 1, u%n_data_used
        d_ptr => u%data(i)
        if (.not. d_ptr%spin_map%valid) cycle
        if (nl+100 > size(lines)) call re_allocate(lines, nl+200, .false.)
        if (j == 0) then
          nl=nl+1; lines(nl) = ''
          nl=nl+1; lines(nl) = 'Spin G-matrices used in data:'
        else
          nl=nl+1; lines(nl) = ''
        endif
        j = j + 1
        ix = d_ptr%spin_map%ix_ref
        if (ix < 0) ix = d_ptr%spin_map%ix_ele
        nl=nl+1; write(lines(nl), '(2x, a, i0, a, i0)')     'Universe: ', d_ptr%spin_map%ix_uni, '  of: ', ubound(s%u, 1)
        nl=nl+1; write(lines(nl), '(2x, a, 2i6)')    'Ix_Ref, Ix_Ele:', ix, d_ptr%spin_map%ix_ele
        nl=nl+1; write(lines(nl), '(26x, a, 26x, a)') 'Initial', 'Final'
        nl=nl+1; write(lines(nl), '(2x, a, 3f12.8, 5x, 3f12.8)') 'L-axis: ', d_ptr%spin_map%axis0%l, d_ptr%spin_map%axis1%l
        nl=nl+1; write(lines(nl), '(2x, a, 3f12.8, 5x, 3f12.8)') 'N0-axis:', d_ptr%spin_map%axis0%n0, d_ptr%spin_map%axis1%n0
        nl=nl+1; write(lines(nl), '(2x, a, 3f12.8, 5x, 3f12.8)') 'M-axis: ', d_ptr%spin_map%axis0%m, d_ptr%spin_map%axis1%m
        nl=nl+1; write(lines(nl), '(2x, a)') '8x8 matrix:'
        do j = 1, 8
          nl=nl+1; write(lines(nl), '(5x, a)') reals_to_table_row(d_ptr%spin_map%mat8(j,:), 13, 7)
        enddo

        if (show_q) then
          nl=nl+1; lines(nl) = ''
          nl=nl+1; lines(nl) = 'Spin quaternion map - 1st order.'
          nl=nl+1; write (lines(nl), '(14x, a, 7x, 6(2x, a, 7x))') '0th order', 'dx  ', 'dpx', 'dy ', 'dpy', 'dz ', 'dpz'
          do ix = 0, 3
            nl=nl+1; write(lines(nl), '(i8, f12.6, 4x, a)') ix, d_ptr%spin_map%map1%spin_q(ix,0), reals_to_string(d_ptr%spin_map%map1%spin_q(ix,1:), 11, 1, 6, 6)
          enddo
          if (d_ptr%spin_map%ix_ref == d_ptr%spin_map%ix_ele) then
            dn_dpz = spin_dn_dpz_from_qmap(real(d_ptr%spin_map%map1%orb_mat, rp), real(d_ptr%spin_map%map1%spin_q, rp), dn_partial, dn_partial2, err)
            nl=nl+1; lines(nl) = ''
            nl=nl+1; write(lines(nl), '(3x, a, 3es14.6)') 'dn/dpz:', dn_dpz
          endif
        endif
      enddo
    enddo

  ! what_to_show = element
  else
    call tao_pick_universe (ele_name, ele_name, picked_uni, err, ix_u)
    if (err) return
    u => s%u(ix_u)
    call tao_locate_elements (ele_name, ix_u, eles, err)
    if (err) then
      nl = 0
      return
    endif
    ele => eles(1)%ele

    if (ele_ref_name /= '') then
      call tao_locate_elements (ele_ref_name, ix_u, eles, err)
      if (err) return
      ele_ref => eles(1)%ele
    else
      ele_ref => pointer_to_next_ele(ele, -1)
    endif

    ! If one-turn calc: Use PTC n0 to be consistant
    if (all(sm%axis_input%n0 == 0) .and. ele_ref%ix_ele /= ele%ix_ele) then
      if (all(u%model%tao_branch(ele%ix_branch)%orbit(ele_ref%ix_ele)%spin == 0)) call tao_lattice_calc(ok)
      sm%axis_input%n0 = u%model%tao_branch(ele%ix_branch)%orbit(ele_ref%ix_ele)%spin
      if (flip) sm%axis_input%n0 = -sm%axis_input%n0 
    endif

    if (show_mat) then
      datum%ix_branch = ix_branch
      call tao_spin_matrix_calc (datum, u, ele_ref, ele, excite_zero)
      if (.not. sm%valid) then
        nl=nl+1; lines(nl) = 'INVALID: ' // datum%why_invalid
        return
      endif

      if (all(sm%axis_input%n0 == 0) .and. ele_ref%ix_ele /= ele%ix_ele) then
        nl=nl+1; lines(nl) = 'NO SPIN AXIS COMPUTED.' 
        nl=nl+1; lines(nl) = 'TO TURN SPIN TRACKING ON FROM THE COMMAND LINE: "set bmad spin_tracking_on = T"'
        nl=nl+1; lines(nl) = 'TO TURN SPIN TRACKING ON IN THE LATTICE FILE: "bmad_com[spin_tracking_on] = T"'
        return
      endif
      nl=nl+1; write (lines(nl), '(23x, a, 34x, a)') 'Initial', 'Final'
      nl=nl+1; write (lines(nl), '(a, 3f12.8, 5x, 3f12.8)') 'L-axis:', sm%axis0%l, sm%axis1%l
      nl=nl+1; write (lines(nl), '(a, 3f12.8, 5x, 3f12.8)') 'N-axis:', sm%axis0%n0, sm%axis1%n0
      nl=nl+1; write (lines(nl), '(a, 3f12.8, 5x, 3f12.8)') 'M-axis:', sm%axis0%m, sm%axis1%m
      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = '8x8 Matrix:'
      do i = 1, 8
        nl=nl+1; write (lines(nl), '(5x, a)') reals_to_table_row(sm%mat8(i,:), 13, 7)
      enddo
    endif

    if (show_q) then
      call spin_concat_linear_maps(err, sm%map1, branch, ele_ref%ix_ele, ele%ix_ele, &
                                             orbit = tao_branch%orbit, excite_zero = excite_zero)
      if (err) return

      nl=nl+1; lines(nl) = ''
      nl=nl+1; lines(nl) = 'Spin quaternion map - 1st order' 
      nl=nl+1; write (lines(nl), '(15x, a, 7x, 6(4x, a, 7x))') '0th order', 'dx  ', 'dpx', 'dy ', 'dpy', 'dz ', 'dpz'
      do i = 0, 3
        nl=nl+1; write(lines(nl), '(i4, 2x, a, 2x, f14.8, 4x, a)') i, q_name(i), sm%map1%spin_q(i,0), reals_to_string(sm%map1%spin_q(i,1:), 13, 1, 8, 8)
      enddo
      if (ele_ref%ix_ele == ele%ix_ele) then
        dn_dpz = spin_dn_dpz_from_qmap(sm%map1%orb_mat, sm%map1%spin_q, dn_partial, dn_partial2, err)
        nl=nl+1; lines(nl) = ''
        nl=nl+1; write(lines(nl), '(3x, a, 3es14.6)') 'dn/dpz:', dn_dpz
      endif
    endif

    if (ele_ref%ix_ele == ele%ix_ele) then
      call spin_mat_to_eigen (sm%map1%orb_mat, sm%map1%spin_q, eval, evec, n0, n_eigen, err)
      if (err) return
      if (dot_product(n0, sm%axis0%n0) < 0) n_eigen = -n_eigen

      if (ele%ix_ele == ele_ref%ix_ele) then
        nl=nl+1; lines(nl) = ''
        nl=nl+1; lines(nl) = 'Eigen:'
        nl=nl+1; lines(nl) = '     |Eval|     E_val            x           px          y           py          z           pz              Sx          Sy          Sz'

        do i = 1, 6
          nl=nl+1; write (lines(nl), '(a, 2f10.6,     4x, 6f12.6, 4x, 3es12.4)', iostat = ios) 're', abs(eval(i)), real(eval(i),rp), real(evec(i,:),rp), real(n_eigen(i,:),rp)
          nl=nl+1; write (lines(nl), '(a, 10x, f10.6, 4x, 6f12.6, 4x, 3es12.4)', iostat = ios) 'im', aimag(eval(i)), aimag(evec(i,:)), aimag(n_eigen(i,:))
          nl=nl+1; lines(nl) = ''
        enddo

        nl=nl+1; lines(nl) = ''
        nl=nl+1; lines(nl) = 'Resonance strengths:'
        nl=nl+1; lines(nl) = '     Orb_Tune         frac(Qs+Q)    Xi_sum          frac(Qs-Q)   Xi_diff'

        do i = 1, 3
          j = 2 * i - 1
          q = atan2(aimag(eval(j)), real(eval(j),rp)) / twopi
          call spin_quat_resonance_strengths(evec(j,:), sm%map1%spin_q, xi_sum, xi_diff)
          nl=nl+1; write (lines(nl), '(a, f13.7, 2(f17.7, es13.5))') abc_name(i), &
                      q, modulo2(qs+q, 0.5_rp), xi_sum, modulo2(qs-q, 0.5_rp), xi_diff
        enddo
        nl=nl+1; lines(nl) = 'Note: "help show spin" will display information on this table.'
      endif
    endif
  endif

  result_id = 'spin:' // what_to_show
  return

!----------------------------------------------------------------------
! string

case ('string')

  nc = 0
  map_calc = u%calc%one_turn_map

  if (index(what2, 'chrom') /= 0) then
    s%com%force_chrom_calc = .true.
    s%u%calc%lattice = .true.
    if (index(what2, 'ptc') /= 0) u%calc%one_turn_map = .true.
    call tao_lattice_calc(ok)
  endif

  do
    ix = index(what2, '`')
    if (ix == 0) then
      line = line(1:nc) // what2
      exit
    endif

    line = line(1:nc) // what2(:ix-1)
    nc = nc + ix - 1
    what2 = what2(ix+1:)
    ix = index(what2, '`')
    if (ix == 0) then
      nl=nl+1; lines(nl) = 'UNMATCHED BACKTICK.'
      u%calc%one_turn_map = map_calc
      return
    endif

    str = what2(1:ix-1)
    what2 = what2(ix+1:)

    n_order = 14

    n = index(str, '@@')
    if (n /= 0) then
      if (.not. is_integer(str(n+2:), n_order)) then
        nl=nl+1; lines(nl) = 'Not an integer after "@@": ' // str(n+2:)
        u%calc%one_turn_map = map_calc
        return
      endif
      str = str(:n-1)
    endif

    ! If str evaluates to an element parameter. EG: "ele::34[name]".

    call tao_parse_element_param_str(err, str, uni_str, ele_name, param_name, loc, component)
    param_name = upcase(param_name)
    if (.not. err .and. loc == anchor_end$ .and. ele_name /= '' .and. component == ''.and. &
                                                            attribute_type(param_name) == is_string$) then
      u => tao_pointer_to_universe(uni_str)
      if (associated(u)) then
        call tao_locate_elements(ele_name, u%ix_uni, eles, err, multiple_eles_is_err = .true.)
        if (.not. err) then
          call string_attrib(param_name, eles(1)%ele, str)
          line = line(1:nc) // str
          nc = len_trim(line)
          cycle
        endif
      endif
    else
      call tao_evaluate_expression (str, 0, .false., value, err)
    endif

    if (err) then
      line = line(1:nc) // '??????'
      nc = len_trim(line)
      cycle
    endif

    !

    if (size(value) == 1) then
      line = line(1:nc) // real_str(value(1), n_order)
    else
      line = line(1:nc) // '[' // real_str(value(1), n_order)
      nc = len_trim(line)
      do i = 2, size(value)
        line = line(1:nc) // ', ' // real_str(value(i), n_order)
        nc = len_trim(line)
      enddo
      line = line(1:nc) // ']'
    endif

    nc = len_trim(line)
  enddo

  u%calc%one_turn_map = map_calc

  do
    ix = index(line, '\n')
    if (ix == 0) exit
    nl=nl+1; lines(nl) = line(:ix-1)
    line = line(ix+2:)
  enddo

  nl=nl+1; lines(nl) = line

!----------------------------------------------------------------------
! symbolic_numbers

case ('symbolic_numbers')

  what_to_show = 'tao'

  do
    call tao_next_switch (what2, [character(24):: '-physical_constants', '-lattice_constants'], .true., switch, err)
    if (err) return
    select case (switch)
    case ('');           exit
    case ('-physical_constants');    what_to_show = 'physical'
    case ('-lattice_constants');     what_to_show = 'lattice'
    end select
  enddo

  select case (what_to_show)
  case ('tao')
    if (allocated(s%com%symbolic_num)) then
      n = nl+10+size(s%com%symbolic_num); if (n > size(lines)) call re_allocate (lines, n, .false.)
      do i = 1, size(s%com%symbolic_num)
        nl=nl+1; write (lines(nl), '(2x, 2a, es22.15)') s%com%symbolic_num(i)%name, '=', s%com%symbolic_num(i)%value
      enddo
    else
      nl=nl+1; lines(nl) = 'No symbolic numbers yet defined. Use "show sym -phys" to show physical constants'
    endif

  case ('physical')
    n = nl+10+size(physical_const_list); if (n > size(lines)) call re_allocate (lines, n, .false.)
    do i = 1, size(physical_const_list)
      if (physical_const_list(i)%name == 'emass' .or. physical_const_list(i)%name == 'pmass') then
        nl=nl+1; write (lines(nl), '(2x, 2a, es22.15, a)') physical_const_list(i)%name, '=', physical_const_list(i)%value, &
                            '  ! For compatibility with MAD. Please avoid using this constant.'
      else
        nl=nl+1; write (lines(nl), '(2x, 2a, es22.15)') physical_const_list(i)%name, '=', physical_const_list(i)%value
      endif
    enddo

  case ('lattice')
    if (allocated(lat%constant)) then
      n = nl+10+size(lat%constant); if (n > size(lines)) call re_allocate (lines, n, .false.)
      do i = 1, size(lat%constant)
        nl=nl+1; write (lines(nl), '(2x, 2a, es22.15)') lat%constant(i)%name, '=', lat%constant(i)%value
      enddo
    else
      nl=nl+1; lines(nl) = 'No constants were defined in the lattice.'
    endif
  end select

!----------------------------------------------------------------------
! taylor_map

case ('taylor_map', 'matrix')

  ix_branch = s%global%default_branch
  where = 'ELE_LOCS' ! "S_LOCS", "ELE_LIST", "ELE_LOCS" 
  print_eigen = .false.
  angle_units = .false.
  force_use_ptc = .false.
  do_inverse = .false.
  clean = .true.
  output_type = 'TAYLOR_STANDARD'  ! "BMAD_LATTICE_FORMAT", "SCIBMAD", "RADIATION", "MATRIX", "TAYLOR_STANDARD"
  fmt = ''
  ele1_name = ''
  ele2_name = ''
  style = ''

  if (show_what == 'matrix') then
    n_order = 1
  else
    n_order = -1
  endif

  do
    call tao_next_switch (what2, [character(20):: '-order', '-s', '-ptc', '-eigen_modes', '-elements', &
              '-lattice_format', '-universe', '-angle_coordinates', '-number_format', '-inverse', &
              '-radiation', '-scibmad', '-noclean'], .true., switch, err)
    if (err) return
    if (switch == '') exit

    select case (switch)
    case ('-angle_coordinates')
      angle_units = .true.

    case ('-noclean')
      clean = .false.

    case ('-eigen_modes')
      print_eigen = .true.

    case ('-elements')
      output_type = 'ELEMENTS'
      call tao_next_word(what2, ele_name)
      if (n_order == -1) n_order = 1

    case ('-inverse')
      do_inverse = .true.

    case ('-lattice_format')
      output_type = 'BMAD_LATTICE_FORMAT'

    case ('-number_format')
      call tao_next_word(what2, fmt)

    case ('-order')
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios) n_order
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ ORDER NUMBER!'
        return
      endif
      if (n_order < 2 .and. output_type == 'TAYLOR_STANDARD') output_type = 'MATRIX'

    case ('-ptc')
      force_use_ptc = .true.

    case ('-radiation')
      output_type = 'RADIATION'

    case ('-s')
      where = "S_LOC"

    case ('-scibmad')
      output_type = 'SCIBMAD'
      style = 'SCIBMAD'
      call tao_next_word(what2, var_name)
      
    case ('-universe')
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios) ix
      u => tao_pointer_to_universe(ix)
      if (ix == 0 .or. ios /= 0 .or. .not. associated(u)) then
        nl=1; lines(1) = 'CANNOT READ OR OUT-OF RANGE "-universe" ARGUMENT'
        return
      endif
      lat => u%model%lat

    case default
      if (ele1_name == '') then
        ele1_name = switch
      elseif (ele2_name == '') then
        ele2_name = switch
      else
        nl=1; lines(1) = 'EXTRA STUFF ON LINE!'
        return
      endif
    end select
  enddo

  !

  if (n_order == -1) n_order = ptc_private%taylor_order_ptc

  if (where == "ELE_LIST" .and. output_type == "RADIATION") then
    nl=1; lines(1) = 'ERROR: "-radiation" AND "-elements" SWITCHES CANNOT BOTH BE PRESENT.'
    return
  endif

  if (where == "S_LOC" .and. output_type == "RADIATION") then
    nl=1; lines(1) = 'ERROR: "-radiation" AND "-s" SWITCHES CANNOT BOTH BE PRESENT.'
    return
  endif

  !---------------------------------------
  ! By s

  if (where == "S_LOC") then
    if (ele1_name == '') then
      s2 = lat%ele(lat%n_ele_track)%s
      s1 = 0
    else
      call tao_evaluate_expression(ele1_name, 1, .false., value, err)
      if (err) then
        nl=1; lines(1) = 'BAD S1 VALUE:' // ele1_name
        return
      endif
      s1 = value(1)
    endif

    if (ele2_name == '') then
      s2 = s1
      if (lat%param%geometry == open$) s1 = 0
    else 
      call tao_evaluate_expression(ele2_name, 1, .false., value, err)
      if (err) then
        nl=1; lines(1) = 'BAD S2 VALUE:' // ele2_name
        return
      endif
      s2 = value(1)
    endif

    call twiss_and_track_at_s (lat, s1, ele0, u%model%tao_branch(ix_branch)%orbit, orb, ix_branch)

    if (n_order > 1 .or. force_use_ptc) then
      call transfer_map_from_s_to_s (lat, taylor, s1, s2, orb, orb2, ix_branch = ix_branch, &
                                                        one_turn = .true., concat_if_possible = s%global%concatenate_maps)
      call taylor_to_mat6(taylor, orb%vec, vec0, mat6)
      ref_vec = taylor%ref
    else
      call mat6_from_s_to_s (lat, mat6, vec0, s1, s2, orb, orb2, ix_branch, one_turn = .true.)
      ref_vec = orb%vec
    endif

    if (do_inverse) then
      if (n_order > 1 .or. force_use_ptc) call taylor_inverse (taylor, taylor)
      call mat_inverse(mat6, mat6, vec0, vec0)
      ref_vec = orb2%vec
    endif

  !

  elseif (where == 'ELE_LOCS') then

    if (ele1_name == '') then
      ix2 = lat%n_ele_track
      ix1 = 0
    else
      call tao_locate_elements (ele1_name, u%ix_uni, eles, err)
      if (size(eles) > 1) then
        nl=1; lines(1) = 'MULTIPLE ELEMENTS BY THIS NAME: ' // ele1_name
        return
      endif
      if (err .or. size(eles) == 0) return
      ele => eles(1)%ele

      select case (ele%lord_status)
      case (super_lord$) 
        ele => pointer_to_slave (ele, ele%n_slave)
      case (overlay_lord$, multipass_lord$, girder_lord$, group_lord$)
        nl=1; lines(1) = 'LORD ELEMENT OF THIS TYPE (' // trim(control_name(ele%lord_status)) // &
                         ') DOES NOT HAVE A DEFINITE POSITION: ' // ele1_name
        return
      end select

      ix1 = ele%ix_ele
      ix_branch = ele%ix_branch
    endif

    branch => lat%branch(ix_branch)

    if (ele2_name == '' .and. ele1_name /= '') then
      ix2 = ix1
      if (lat%param%geometry == open$) then
        ix1 = 0
      else
        ix1 = ix2
      endif

    elseif (ele2_name /= '') then
      call tao_locate_elements (ele2_name, u%ix_uni, eles, err)
      if (size(eles) > 1) then
        nl=1; lines(1) = 'MULTIPLE ELEMENTS BY THIS NAME: ' // ele2_name
        return
      endif
      if (err .or. size(eles) == 0) return
      ele => eles(1)%ele

      select case (ele%lord_status)
      case (super_lord$) 
        ele => pointer_to_slave (ele, ele%n_slave)
      case (overlay_lord$, multipass_lord$, girder_lord$, group_lord$)
        nl=1; lines(1) = 'LORD ELEMENT OF THIS TYPE (' // trim(control_name(ele%lord_status)) // &
                         ') DOES NOT HAVE A DEFINITE POSITION: ' // ele1_name
        return
      end select

      ix2 = ele%ix_ele
      if (ele%ix_branch /= ix_branch) then
        nl=1; lines(1) = 'ELEMENTS ARE IN DIFFERENT LATTICE BRANCHES'
        return
      endif
    endif
  endif

    !---------------------------------------
    ! Radiation

  if (output_type == 'RADIATION') then
    if (.not. s%global%rf_on) then
      nl=nl+1; lines(nl) = 'Note: RF IS TURNED OFF!'
      nl=nl+1; lines(nl) = '   TO TURN ON USE: "set global rf_on = T"'
    endif

    if (.not. associated (branch%ptc%m_t_layout)) call lat_to_ptc_layout (lat)

    ! Use PTC or Bmad closed orbit for the radiation calculation?
    if (force_use_ptc) then
      call ptc_setup_map_with_radiation (rad_map, branch%ele(ix1), branch%ele(ix2), 1, .true.)
    else
      call ptc_setup_map_with_radiation (rad_map, branch%ele(ix1), branch%ele(ix2), &
                                         1, .true., orbit1 = u%model%tao_branch(ix_branch)%orbit(ix1))
    endif
    rmap => rad_map%sub_map
    nl=nl+1; write (lines(nl), '(a, 6es16.8)') 'Ref_orb_start: ', rmap(1)%fix0
    nl=nl+1; write (lines(nl), '(a, 6es16.8)') 'Ref_orb_end:   ', rmap(1)%fix
    nl=nl+1; lines(nl) = ''
    nl=nl+1; call mat_type (rad_map%nodamp_mat, 0, 'T-Matrix without Damping. Symplectic Error: ' // real_str(mat_symp_error(rad_map%nodamp_mat), 6), '(4x, 6es16.8)', lines(nl:), n)
    nl=nl+7; lines(nl) = ''
    nl=nl+1; call mat_type (rad_map%damp_mat, 0, 'D-Damping Matrix. Damp Factor: ' // real_str((1-determinant(rad_map%damp_mat))/10, 6), '(4x, 6es16.8)', lines(nl:), n)
    nl=nl+7; lines(nl) = ''
    nl=nl+1; call mat_type (rad_map%stoc_mat, 0, 'S-Radiation Matrix:', '(4x, 6es16.8)', lines(nl:), n)
    nl=nl+6
    return
  endif

  !---------------------------------------

  if (where == 'ELE_LOCS') then
    if (ele2_name == '') then
      nl=nl+1; lines(nl) = 'From: ' // trim(branch%ele(ix1)%name)
      nl=nl+1; lines(nl) = 'To:   ' // trim(branch%ele(ix2)%name)
    endif

    if (output_type /= 'MATRIX' .or. force_use_ptc) then
      call transfer_map_calc (lat, taylor, err, ix1, ix2, u%model%tao_branch(ix_branch)%orbit(ix1), ix_branch, &
                                                      one_turn = .true., concat_if_possible = s%global%concatenate_maps)
      if (err) then
        nl = 1; lines(1) = 'TAYLOR MAP TERM OVERFLOW.'
        return
      endif

      call taylor_to_mat6(taylor, u%model%tao_branch(ix_branch)%orbit(ix1)%vec, vec0, mat6)
      ref_vec = u%model%tao_branch(ix_branch)%orbit(ix1)%vec

    else
      call transfer_matrix_calc (lat, mat6, vec0, ix1, ix2, ix_branch, one_turn = .true.)
      ref_vec = u%model%tao_branch(ix_branch)%orbit(ix1)%vec
    endif
  endif

  if (do_inverse) then
    if (output_type /= 'MATRIX' .or. force_use_ptc) call taylor_inverse (taylor, taylor)
    call mat_inverse(mat6, mat6, vec0, vec0)
    ref_vec = u%model%tao_branch(ix_branch)%orbit(ix2)%vec
  endif

  ! Print results

  if (output_type == 'ELEMENTS') THEN
    call tao_locate_elements (ele_name, u%ix_uni, eles, err)
    if (err .or. size(eles) == 0) return

    do i = 1, size(eles)
      ele => eles(i)%ele
      select case (ele%key)
      case (group$, overlay$, girder$, ramper$); cycle
      end select
      if (i > 1) then
        nl=nl+1; lines(nl) = ''
      endif
      nl=nl+1; write (lines(nl), '(4a)') 'Element: ', trim(ele%name), ' ', ele_loc_name(ele, .false., '()')
      if (size(lines) < nl+400) call re_allocate(lines, 2*size(lines))

      if (output_type == 'MATRIX') then
        mat6 = ele%mat6
        vec0 = ele%vec0

        if (angle_units) then
          call mat6_to_taylor (vec0, mat6, taylor, ref_vec)
          call map_to_angle_coords (taylor, taylor)
          call taylor_to_mat6 (taylor, ref_vec, vec0, mat6)
        endif

        if (nl+10 > size(lines)) call re_allocate (lines, 2*nl, .false.)

        if (fmt /= '') then
          fmt2 = '(6' // trim(fmt) // ', a, ' // trim(fmt) // ')'
        elseif (any(abs(mat6) >= 1000)) then
          fmt2 = '(6es15.7, a, es15.7)'
        else
          fmt2 = '(6f15.8, a, es15.7)'
        endif

        do j = 1, 6
          nl=nl+1; write(lines(nl), fmt2, iostat = ios) mat6(j,:), '   : ', vec0(j)
        enddo

      else  ! output_type /= 'MATRIX'
        i0 = ele%ix_ele-1
        call transfer_map_calc (lat, taylor, err, i0, ele%ix_ele, u%model%tao_branch(ix_branch)%orbit(i0), ele%ix_branch)
        if (do_inverse) call taylor_inverse(taylor, taylor)
        call truncate_taylor_to_order (taylor, n_order, taylor)
        call type_taylors (taylor, n_order, alloc_lines, n, out_style = style, clean = clean, out_var_suffix = var_name)
        do j = 1, n
          nl=nl+1; lines(nl) = alloc_lines(j)
        enddo
      endif
    enddo

  else
    if (output_type /= 'MATRIX') then
      if (angle_units) call map_to_angle_coords (taylor, taylor)
      if (n_order > 1) call truncate_taylor_to_order (taylor, n_order, taylor)
      call type_taylors (taylor, n_order, lines, n_lines = nl, out_style = style, clean = clean, out_var_suffix = var_name)
      if (print_eigen) call taylor_to_mat6 (taylor, taylor%ref, vec0, mat6)

    elseif (output_type == 'BMAD_LATTICE_FORMAT') then
      call mat6_to_taylor (vec0, mat6, taylor, ref_vec)
      call type_taylors (taylor, n_order, lines, nl, out_style = style, clean = clean, out_var_suffix = var_name)

    else
      if (angle_units) then
        call mat6_to_taylor (vec0, mat6, taylor, ref_vec)
        call map_to_angle_coords (taylor, taylor)
        call taylor_to_mat6 (taylor, ref_vec, vec0, mat6)
      endif

      vec_in = 0
      if (n_order == 0) then 
        nl = nl+1; write(lines(nl), '(6f11.6)') vec0
      else
        if (fmt /= '') then
          fmt2 = '(6' // trim(fmt) // ', a, ' // trim(fmt) // ')'
        elseif (any(abs(mat6) >= 1000)) then
          fmt2 = '(6es15.7, a, es15.8)'
        else
          fmt2 = '(6f15.8, a, es15.8)'
        endif

        nl=nl+1; write (lines(nl), '(a, es10.2)') 'Symplectic Error:', mat_symp_error(mat6)
        nl=nl+1; lines(nl) = ''
        do i = 1, 6
          nl=nl+1; write(lines(nl), fmt2, iostat = ios) mat6(i,:), '   : ', vec0(i)
        enddo
      endif
    endif

    if (print_eigen) then
      call mat_eigen (mat6, eval, evec, err)
      nl=nl+1; lines(nl) = ''
      nl=nl+1; write(lines(nl), '(75x, a)') 'eVector'
      nl=nl+1; write(lines(nl), '(t11, a, t29, a, 3(15x, a, 14x, a))') '|eValue|', 'eValue', 'x', 'px', 'y', 'py', 'z', 'pz'

      if (fmt == '') then
        fmt2 = '(a, 8es16.8)'
        fmt3 = '(a, 16x, 8es16.8)'
      else
        fmt2 = '(a, 8' // trim(fmt) // ')'
        fmt3 = '(a, 16x, 8' // trim(fmt) // ')'
      endif

      do i = 1, 6
        nl=nl+1; write (lines(nl), fmt2, iostat = ios) 're', abs(eval(i)), real(eval(i), rp), real(evec(i,:), rp)
        nl=nl+1; write (lines(nl), fmt3, iostat = ios) 'im',              aimag(eval(i)), aimag(evec(i,:))
        nl=nl+1; lines(nl) = ''
      enddo
      nl = nl-1
    endif  
  endif

  !

  if (n_order > ptc_private%taylor_order_ptc) then
    nl=nl+1; write(lines(nl), '(a, i0)') 'NOTE: -order value is above order used in calculations which is ', ptc_private%taylor_order_ptc
    nl=nl+1; write(lines(nl), '(a, i0)') '      Thus the -order value is ignored.'
  endif

!----------------------------------------------------------------------
! track

case ('track')

  what_to_show = ''
  ele_name = ''
  velocity_fmt = ''
  e_field_fmt = ''
  b_field_fmt = ''
  position_fmt = '3pf14.6'
  s_fmt = 'f13.6'
  momentum_fmt = 'f14.8'
  t_fmt = ''
  spin_fmt = ''
  energy_fmt = ''
  twiss_fmt = ''
  disp_fmt = ''
  print_header_lines = .true.
  s1 = branch%ele(0)%s
  s2 = branch%ele(branch%n_ele_track)%s
  n_print = s%plot_page%n_curve_pts
  tao_lat => u%model
  branch => tao_lat%lat%branch(s%global%default_branch)
  lat_type = model$

  do 
    call tao_next_switch (what2, [character(16):: '-e_field', '-b_field', '-velocity', '-momentum', &
                '-energy', '-position', '-no_label_lines', '-s', '-spin', '-points', '-time', &
                '-range', '-twiss', '-dispersion', '-branch', '-universe', '-design', '-base', '-element'], &
                .false., switch, err)

    if (err) return
    if (switch == '') exit

    select case (switch)
    case ('-e_field')
      e_field_fmt = get_this_track_fmt(what2, 'es15.6', err); if (err) return
    case ('-b_field')
      b_field_fmt = get_this_track_fmt(what2, 'es15.6', err); if (err) return
    case ('-element')
      call tao_next_word(what2, ele_name)
    case ('-energy')
      energy_fmt = get_this_track_fmt(what2, 'es15.6', err); if (err) return
    case ('-velocity')
      velocity_fmt = get_this_track_fmt(what2, 'f13.8', err); if (err) return
    case ('-momentum')
      momentum_fmt = get_this_track_fmt(what2, momentum_fmt, err); if (err) return
    case ('-twiss')
      twiss_fmt = get_this_track_fmt(what2, 'f14.6', err); if (err) return
    case ('-dispersion')
      disp_fmt = get_this_track_fmt(what2, 'f14.6', err); if (err) return
    case ('-position')
      position_fmt = get_this_track_fmt(what2, position_fmt, err); if (err) return
    case ('-no_label_lines')
      print_header_lines = .false.
    case ('-points')
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios) n_print
      if (aname == '' .or. ios /= 0 .or. n_print < 1) then
        nl=1; lines(1) = 'CANNOT READ OR OUT-OF RANGE "-points" ARGUMENT'
        return
      endif
    case ('-range')
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios1) s1
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios2) s2
      if (aname == '' .or. ios1 /= 0 .or. ios2 /= 0) then
        nl=1; lines(1) = 'CANNOT READ OR OUT-OF RANGE "-points" ARGUMENT'
        return
      endif
    case ('-s')
      s_fmt = get_this_track_fmt(what2, s_fmt, err); if (err) return
    case ('-spin')
      spin_fmt = get_this_track_fmt(what2, 'f13.8', err); if (err) return
    case ('-time')
      t_fmt = get_this_track_fmt(what2, 'es15.6', err); if (err) return
    case ('-base')
      lat_type = base$
    case ('-branch')
      call tao_next_word(what2, aname)
      branch => pointer_to_branch(aname, lat)
      ix_branch = branch%ix_branch
      if (.not. associated(branch)) then
        nl=1; write(lines(1), *) 'Bad branch index:', ix_branch
        return
      endif
    case ('-design')
      lat_type = design$
    case ('-universe')
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios) ix
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ OR OUT-OF RANGE "-universe" argument'
        return
      endif
      call point_to_this_uni(ix, ix_branch, u, lat, branch, model_branch, tao_branch, design_tao_branch)
      if (.not. associated(u)) return
    end select
  enddo

  !

  if (ele_name /= '') then
    call tao_locate_elements (ele_name, u%ix_uni, eles, err, lat_type, multiple_eles_is_err = .true.)
    if (err) return
    ele => eles(1)%ele
    s1 = ele%s_start
    s2 = ele%s
    if (ele%key == beambeam$) what_to_show = 'beambeam'
  endif

  ! Beambeam element.

  if (what_to_show == 'beambeam') then
    call track1 (tao_branch%orbit(ele%ix_ele-1), ele, ele%branch%param, orb2, track, err)
    call re_allocate(lines, nl+track%n_pt+10)

    nl=nl+1; lines(nl) = '                      |                 Tracked particle (Laboratory Coordinates)                     |         Strong Beam (Lab Coords)                            | Particle - Beam distance'
    nl=nl+1; lines(nl) = '     s_lab     s_body |       x            px           y            py           z            pz     | slice     x_center     y_center      x_sigma      y_sigma   |         dx           dy'

    do i = 0, track%n_pt
      tp => track%pt(i)
      sb => tp%strong_beam
      if (sb%ix_slice == 0) then
        nl=nl+1; write (lines(nl), '(2f11.6, 1x, 6es13.5)') tp%s_lab, tp%s_body, tp%orb%vec
      else
        nl=nl+1; write (lines(nl), '(2f11.6, 1x, 6es13.5, i8, 2x, 4es13.5, 2x, 2es13.5)') tp%s_lab, tp%s_body, tp%orb%vec, &
                                  sb%ix_slice, sb%x_center, sb%y_center, sb%x_sigma, sb%y_sigma, sb%dx, sb%dy
      endif
    enddo
    return
  endif

  !

  if (print_header_lines) then
    line1 = '#   Ix'
    i1 = 7
    call write_track_header (line1, i1, s_fmt,        ['S'], err); if (err) return
    call write_track_header (line1, i1, t_fmt,        ['Time'], err); if (err) return
    call write_track_header (line1, i1, position_fmt, ['X', 'Y', 'Z'], err); if (err) return
    call write_track_header (line1, i1, velocity_fmt, ['Vx/c', 'Vy/c', 'Vs/c'], err); if (err) return
    call write_track_header (line1, i1, momentum_fmt, ['px', 'py', 'pz'], err); if (err) return
    call write_track_header (line1, i1, energy_fmt,   ['E_tot'], err); if (err) return
    call write_track_header (line1, i1, twiss_fmt,    ['Beta_a ', 'Alpha_a', 'Beta_b ', 'Alpha_b'], err); if (err) return
    call write_track_header (line1, i1, disp_fmt,     ['Eta_x ', 'Etap_x', 'Eta_y ', 'Etap_y'], err); if (err) return
    call write_track_header (line1, i1, spin_fmt,     ['Spin_x', 'Spin_y', 'Spin_z'], err); if (err) return
    call write_track_header (line1, i1, b_field_fmt,  ['Bx', 'By', 'Bz'], err); if (err) return
    call write_track_header (line1, i1, e_field_fmt,  ['Ex', 'Ey', 'Ez'], err); if (err) return
    nl=nl+1; lines(nl) = line1
  endif

  !    

  if (s1 < 0 .and. .not. tao_branch%plot_cache_valid) then
    nl = 1; lines(nl) = 'PROBLEM: "-s" BOUNDS NOT SPECIFIED AND NO PLOT DATA AVAILABLE. NO TABLE CAN BE GENERATED.'
    return
  endif

  call re_allocate(lines, nl+n_print+10)

  do_field = ((e_field_fmt /= '' .and. e_field_fmt /= 'no') .or. (b_field_fmt /= '' .and. b_field_fmt /= 'no'))

  do i = 1, n_print
    s_pos = s1 + (i - 1) * (s2 - s1) / max(1, (n_print - 1))
    call twiss_and_track_at_s (lat, s_pos, ele0, tao_branch%orbit, orbit, branch%ix_branch, err, (i /= 1), .false.)

    write (line1, '(i7)') i
    i1 = 7

    call write_track_info (line1, i1, s_fmt, [s_pos], err);  if (err) return
    call write_track_info (line1, i1, t_fmt, [real(rp):: orbit%t], err);  if (err) return
    call write_track_info (line1, i1, position_fmt, orbit%vec(1:5:2), err);  if (err) return

    if (orbit%beta == 0) then
      vec3 = 0
    else
      vec3(1:2) = [orbit%vec(2), orbit%vec(4)] / (1 + orbit%vec(6))
      vec3 = orbit%beta * [vec3(1), vec3(2), sqrt(max(0.0_rp, 1 - vec3(1)**2 - vec3(2)**2))]
    endif
    call write_track_info (line1, i1, velocity_fmt, vec3, err);  if (err) return

    call write_track_info (line1, i1, momentum_fmt, orbit%vec(2:6:2), err);  if (err) return
    call convert_pc_to((1 + orbit%vec(6)) * orbit%p0c,  orbit%species, e_tot = e_tot)
    call write_track_info (line1, i1, energy_fmt, [e_tot], err);  if (err) return
    call write_track_info (line1, i1, twiss_fmt, [ele0%a%beta, ele0%a%alpha, ele0%b%beta, ele0%b%alpha], err);  if (err) return
    call write_track_info (line1, i1, disp_fmt, [ele0%x%eta, ele0%x%etap, ele0%y%eta, ele0%y%etap], err);  if (err) return
    call write_track_info (line1, i1, spin_fmt, [orbit%spin], err);  if (err) return

    if (do_field) then
      call em_field_calc (ele0, branch%param, orbit%s-ele0%s_start, orbit, .false., field)
      call write_track_info (line1, i1, b_field_fmt, field%B, err);  if (err) return
      call write_track_info (line1, i1, e_field_fmt, field%E, err);  if (err) return
    endif

    nl=nl+1; lines(nl) = line1
  enddo

!----------------------------------------------------------------------
! tune

case ('tune')

  nl=nl+1; lines(nl) = 'Use "show universe" instead.'

!----------------------------------------------------------------------
! twiss
    
case ('twiss_and_orbit')

  tao_lat => u%model
  branch => tao_lat%lat%branch(s%global%default_branch)
  lat_type = model$
  attrib0 = ''

  do 

    call tao_next_switch (what2, [character(16):: '-branch', '-universe', '-design', '-base'], .true., switch, err)
    if (err) return
    if (switch == '') exit

    select case (switch)
    case ('-base')
      lat_type = base$

    case ('-branch')
      call tao_next_word(what2, aname)
      branch => pointer_to_branch(aname, lat)
      if (.not. associated(branch)) then
        nl=1; write(lines(1), *) 'Bad branch index:', ix_branch
        return
      endif

    case ('-design')
      lat_type = design$

    case ('-universe')
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios) ix
      u => tao_pointer_to_universe(ix)
      if (aname == '' .or. ios /= 0 .or. .not. associated(u)) then
        nl=1; lines(1) = 'CANNOT READ OR OUT-OF RANGE "-universe" argument'
        return
      endif

    case default
      if (what2 /= '') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON LINE: ' // what2)
        return
      endif
      attrib0 = switch
    end select

  enddo

  !

  tao_lat => tao_pointer_to_tao_lat (u, lat_type)
  lat => tao_lat%lat
  branch => lat%branch(branch%ix_branch)
  ix_branch = branch%ix_branch

  call string_trim(attrib0, attrib0, ix)
  if (ix == 0) then
    s_pos = 0
  else
    call tao_evaluate_expression (attrib0, 1, .false., value, err, .false., &
                                           dflt_component = tao_lat_type_name(lat_type), dflt_uni = u%ix_uni)

    if (err) then
      nl=1; lines(1) = 'NOT A REAL NUMBER: ' // attrib0
      return
    endif
    s_pos = value(1)
  endif

  call twiss_and_track_at_s (lat, s_pos, ele0, tao_lat%tao_branch(ix_branch)%orbit, orb, ix_branch, err)
  if (err) return 
  ele => ele0
  if (associated(ele%lord)) ele => ele0%lord    ! If ele0 is a slice slave

  nl=nl+1; write(lines(nl), '(a, f10.5)')     'At S =', s_pos
  nl=nl+1; write(lines(nl), '(3a, 2(i0, a))') 'In Element: ', trim(ele0%name), '  (', ele%ix_branch, '>>', ele%ix_ele, ')'

  call type_twiss (ele0, s%global%phase_units, lines = lines(nl+1:), n_lines = n)
  nl = nl + n

  if (branch%param%particle /= photon$) then
    ! When evaluating g, stay away from possible element edge fringe fields
    del = 1d-5
    s0 = max(ele%s_start+bmad_com%significant_length, s_pos - del)
    call twiss_and_track_at_s(lat, s0, ele3, tao_lat%tao_branch(ix_branch)%orbit, orb0, ix_branch, err)
    s2 = min(ele%s-bmad_com%significant_length, s_pos + del)
    call twiss_and_track_at_s(lat, s2, ele3, tao_lat%tao_branch(ix_branch)%orbit, orb2, ix_branch, err)
    if (s2 <= s0) then
      nl=nl+1; lines(nl) = 'Cannot evaluate synchrotron radiation parameters in lattice element of negligible width.'
    else
      dr(1:2) = orb2%vec(1:3:2) - orb0%vec(1:3:2)
      dr(3) = s2 - s0
      v0 = [orb0%vec(2), orb0%vec(4), sqrt((1 + orb0%vec(6))**2 - orb0%vec(2)**2 - orb0%vec(4)**2)] / (1 + orb0%vec(6))
      v2 = [orb2%vec(2), orb2%vec(4), sqrt((1 + orb2%vec(6))**2 - orb2%vec(2)**2 - orb2%vec(4)**2)] / (1 + orb2%vec(6))
      g_vec = v2 - v0
      g_vec = g_vec - dr * (dot_product(g_vec, dr) / dot_product(dr, dr))  ! Component of g_vec perpendicular to dr
      g_vec = -g_vec / norm2(dr)
      if (ele%key == sbend$) then
        g_vec(1:2) = g_vec(1:2) + ele%value(g$) * [cos(ele%value(ref_tilt$)), sin(ele%value(ref_tilt$))]
      endif
      g_bend = norm2(g_vec)

      e_tot = ele0%value(e_tot$)
      mc2 = mass_of(branch%param%particle)
      gamma = e_tot / mc2
      E_crit = 3 * h_bar_planck * c_light * g_bend * gamma**3 / 2
      E_ave = 8 * e_crit / (15 * sqrt(3.0_rp))
      c_gamma = 4 * pi * classical_radius_factor / (3 * mc2**4)
      P_gam = c_light * c_gamma * e_tot**4 * g_bend**2 / twopi 
      N_gam = 5 * sqrt(3.0_rp) * c_gamma * e_tot**4 * g_bend / (8 * pi * h_bar_planck * gamma**3)
      N_E2 = 55.0_rp * classical_radius_factor * h_bar_planck * c_light**2 * gamma**7 * g_bend**3 / (24 * sqrt(3.0_rp))
      H_a = ele0%a%gamma * ele0%a%eta**2 + 2 * ele0%a%alpha * ele0%a%eta * ele0%a%etap + ele0%a%beta * ele0%a%etap**2
      H_b = ele0%b%gamma * ele0%b%eta**2 + 2 * ele0%b%alpha * ele0%b%eta * ele0%b%etap + ele0%b%beta * ele0%b%etap**2

      fmt = '(2x, a, es16.8, t40, a)'
      nl=nl+1; lines(nl) = ''
      nl=nl+1; write (lines(nl), '(2x, a, 3f13.8, a)') '[g_x, g_y, g_z]:   ', g_vec,  '  ! g bending strength vector'
      nl=nl+1; write (lines(nl), '(2x, a, f13.8, t40, a)')  'g:                 ', g_bend, '  ! bending strength = 1/rho'
      nl=nl+1; write (lines(nl), fmt) 'E_crit:            ', E_crit, '  ! Critical photon energy (eV)'
      nl=nl+1; write (lines(nl), fmt) '<E_gam>:           ', E_ave,  '  ! Average photon energy (eV)'
      nl=nl+1; write (lines(nl), fmt) 'N_gam:             ', N_gam,  '  ! Photon emission rate per particle (1/sec)'
      nl=nl+1; write (lines(nl), fmt) 'P_gam:             ', P_gam,  '  ! Photon power emitted per particle (eV/sec)'
      nl=nl+1; write (lines(nl), fmt) 'N_gam * <E_gam^2>: ', N_E2,   '  ! #photons * Mean energy^2 (eV^2/sec)'
      nl=nl+1; write (lines(nl), fmt) 'H_a:               ', H_a,    '  ! a-mode curly H'
      nl=nl+1; write (lines(nl), fmt) 'H_b:               ', H_b,    '  ! b-mode curly H'
      nl=nl+1; write (lines(nl), fmt) 'g^3 * H_a          ', H_a * g_bend**3, '  ! Integrand of I_5a radiation integral'
      nl=nl+1; write (lines(nl), fmt) 'g^3 * H_b:         ', H_b * g_bend**3, '  ! Integrand of I_5b radiation integral'
      nl=nl+1; write (lines(nl), fmt) 'H_a * N_gam * <E^2>', H_a * N_E2
    endif
  endif

  fmt = '(2x, a, 3p2f11.4)'
  nl=nl+1; write(lines(nl), *) ''
  nl=nl+1; write(lines(nl), *)   'Orbit: [mm, mrad]'
  nl=nl+1; write(lines(nl), fmt) "X  X':", orb%vec(1), orb%vec(2)
  nl=nl+1; write(lines(nl), fmt) "Y  Y':", orb%vec(3), orb%vec(4)
  nl=nl+1; write(lines(nl), fmt) "Z  Z':", orb%vec(5), orb%vec(6)

!----------------------------------------------------------------------
! universe
! Note: Currently code for -element switch not implemented.
  
case ('universe')

  ix_u = s%global%default_universe
  b_name = ''

  do
    call tao_next_switch (what2, [character(20):: '-branch'], .true., switch, err)
    if (err) return
    select case (switch)
    case ('')
      exit

    case ('-branch')
      b_name = what2(1:ix)

    case default
      read (switch, *, iostat = ios) ix_u
      if (ios /= 0) then
        nl=1; lines(1) = 'BAD UNIVERSE NUMBER'
        return
      endif
      if (ix_u < lbound(s%u, 1) .or. ix_u > ubound(s%u, 1)) then
        nl=1; lines(1) = 'UNIVERSE NUMBER OUT OF RANGE'
        return
      endif
    end select
  enddo

  !

  u => s%u(ix_u)
  lat => u%model%lat

  if (b_name /= '') then
    branch => pointer_to_branch(b_name, u%model%lat)
    if (.not. associated(branch)) then
      nl=1; write(lines(1), *) 'Bad branch name or index: ', b_name
      return
    endif
    ix_branch = branch%ix_branch
  endif

  branch => lat%branch(ix_branch)
  model_branch => u%model_branch(ix_branch)
  tao_branch => u%model%tao_branch(ix_branch)
  tao_lat => tao_pointer_to_tao_lat (u, model$)

  design_lat => u%design%lat
  design_branch => design_lat%branch(ix_branch)
  design_tao_branch => u%design%tao_branch(ix_branch)
  nt = branch%n_ele_track
  species = branch%param%default_tracking_species

  nl = 0
  nl=nl+1; write(lines(nl), '(2(a, i0))') 'Universe: ', ix_u,      '  Of: ', ubound(s%u, 1)
  nl=nl+1; write(lines(nl), '(3(a, i0))') 'Branch:   ', ix_branch, '  In range: [0, ', ubound(lat%branch,1), ']'
  nl=nl+1; write(lines(nl), lmt) 'is_on:                    ', u%is_on
  nl=nl+1; write(lines(nl), amt) 'Lattice name:             ', quote(lat%lattice)
  nl=nl+1; write(lines(nl), amt) 'Machine name:             ', quote(lat%machine)
  nl=nl+1; write(lines(nl), amt) 'Used line(s) in lat file: ', quote(lat%use_name)
  nl=nl+1; write(lines(nl), amt) 'Lattice file name:        ', quote(lat%input_file_name)
  nl=nl+1; write(lines(nl), amt) 'Reference species:        ', species_name(branch%param%particle)

  if (species == ref_particle$ .or. species == anti_ref_particle$) then
    nl=nl+1; write(lines(nl), amt) 'Default tracking species: ', trim(species_name(species)), ' (', trim(species_name(default_tracking_species(branch%param))), ')'
  else
    nl=nl+1; write(lines(nl), amt) 'Default tracking species: ', species_name(species)
  endif

  if (branch%param%particle == photon$) then
    nl=nl+1; write(lines(nl), amt) 'photon_type:                 ', photon_type_name(lat%photon_type)
  else
    species = branch%ele(0)%ref_species
    gamma = branch%ele(0)%value(e_tot$) / mass_of(species)
    nl=nl+1; write(lines(nl), rmt) 'a_anomalous_moment * gamma   ', anomalous_moment_of(species) * gamma
  endif

  if (branch%param%geometry == closed$ .or. branch%param%particle == photon$) then
    nl=nl+1; write(lines(nl), rmt) 'Reference energy:            ', branch%ele(0)%value(e_tot$)
    nl=nl+1; write(lines(nl), rmt) 'Reference momentum:          ', branch%ele(0)%value(p0c$)
  else
    nl=nl+1; write(lines(nl), rmt) 'Starting reference energy:   ', branch%ele(0)%value(e_tot$)
    nl=nl+1; write(lines(nl), rmt) 'Starting reference momentum: ', branch%ele(0)%value(p0c$)
    nl=nl+1; write(lines(nl), rmt) 'Ending reference energy:     ', branch%ele(nt)%value(e_tot$)
    nl=nl+1; write(lines(nl), rmt) 'Ending reference momentum:   ', branch%ele(nt)%value(p0c$)
  endif

  nl=nl+1; write(lines(nl), lmt) 'Absolute_Time_Tracking:      ', bmad_com%absolute_time_tracking
  nl=nl+1; write(lines(nl), amt) 'Geometry:                    ', geometry_name(branch%param%geometry)
  nl=nl+1; write(lines(nl), lmt) 'global%rf_on:                ', s%global%rf_on
  nl=nl+1; write(lines(nl), imt) 'Elements used in tracking: From 1 through ', nt
  if (branch%n_ele_max > nt) then
    nl=nl+1; write(lines(nl), '(2(a, i0))') 'Lord elements:   ', nt+1, '  through ', branch%n_ele_max
  else
    nl=nl+1; write(lines(nl), '(a)') 'There are no Lord elements'
  endif

  nl=nl+1; write(lines(nl), '(a, f0.3)')       'Lattice branch length:      ', branch%param%total_length
  nl=nl+1; write(lines(nl), '(a, es13.6)')     'Lattice branch transit time:', branch%ele(nt)%ref_time - branch%ele(0)%ref_time
  nl=nl+1; write(lines(nl), '(a, 2(f0.6, a))') 'Lattice branch S-range:     [', &
                                                branch%ele(0)%s, ', ', branch%ele(nt)%s, ']'

  if (branch%param%geometry == open$ .and. tao_branch%track_state /= moving_forward$) then
    if (s%global%track_type == 'beam') then
      nl=nl+1; write(lines(nl), '(a, i0)') 'Tracking: Lost beam at:     ', tao_branch%track_state
    else
      nl=nl+1; write(lines(nl), '(a, i0)') 'Tracking: Lost particle at: ', tao_branch%track_state
    endif
  endif

  if (branch%param%geometry == closed$) then
    if (.not. branch%param%stable) then
      nl=nl+1; write(lines(nl), '(a, l)') 'Model lattice stability: ', branch%param%stable
      nl=nl+1; write(lines(nl), '(a, l)') 'Design lattice stability:', design_branch%param%stable
      result_id = 'universe:unstable'
      return
    endif
  endif

  if (allocated(lat%custom)) then
    nl=nl+1; lines(nl) = 'Custom lattice parameters defined in lattice file:'
    do i = 1, size(lat%custom)
      aname = attribute_name(def_parameter$, i+custom_attribute0$)
      if (aname(1:1) == '!') cycle
      nl= nl+1; write (lines(nl), rmt) '  parameter[' // trim(aname) // ']: ', lat%custom(i)
    enddo
  endif
 
  fmt  = '(1x, a16, 2es13.5, 2x, 2es13.5, 2x, a)'
  fmt2 = '(1x, a16,  2f13.6, 2x,  2f13.6, 2x, a)'
  fmt3 = '(1x, a16,         28x, 2es13.5, 2x, a)'
  phase_units = 1 / twopi
  time1 = branch%ele(nt)%ref_time

  !

  if (lat%param%geometry == closed$) then
    u%model%high_e_lat = u%model%lat
    ele2 => u%model%high_e_lat%branch(ix_branch)%ele(0)
    call emit_6d (ele2, .false., tao_branch%modes_6d, sig_mat, tao_branch%orbit, tao_lat%rad_int_by_ele_6d)
    call emit_6d (ele2, .true., tao_branch%modes_6d, sig_mat, tao_branch%orbit, tao_lat%rad_int_by_ele_6d)
    if (tao_branch%modes_6d%a%j_damp < 0 .or. tao_branch%modes_6d%b%j_damp < 0 .or. &
                                           (tao_branch%modes_6d%z%j_damp < 0 .and. rf_is_on(branch))) then
      call out_io (s_info$, r_name, &
        'Negative damping partition number detected therefore the lattice is unstable with radiation excitations.', &
        'Note1: This may not be a problem if the amount of radiation generated is low (like for protons).', &
        'Note2: Instability with respect to radiation excitations does not affect such things as the closed orbit calculation.')
    endif
    call chrom_calc (lat, s%global%delta_e_chrom, tao_branch%a%chrom, tao_branch%b%chrom, &
                                   pz = tao_branch%orbit(0)%vec(6), ix_branch = ix_branch, orb0 = tao_branch%orbit(0))

    mode_d => design_tao_branch%modes_6d
    mode_m => tao_branch%modes_6d
    mode_m%momentum_compaction = momentum_compaction(branch)

    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), '(23x, a)') '       A-Mode         |            B-Mode'
    nl=nl+1; write(lines(nl), '(23x, a)') '  Model       Design  |       Model       Design'

    nl=nl+1; write(lines(nl), fmt2) 'Q', phase_units*branch%ele(nt)%a%phi, &
          phase_units*design_branch%ele(nt)%a%phi, phase_units*branch%ele(nt)%b%phi, phase_units*design_branch%ele(nt)%b%phi,  '! Tune'
    nl=nl+1; write(lines(nl), fmt2) 'Chrom', tao_branch%a%chrom, design_tao_branch%a%chrom, tao_branch%b%chrom, design_tao_branch%b%chrom, '! dQ/(dE/E)'
    nl=nl+1; write(lines(nl), fmt2) 'J_damp', mode_m%a%j_damp, mode_d%a%j_damp, mode_m%b%j_damp, mode_d%b%j_damp, '! Damping Partition #'
    nl=nl+1; write(lines(nl), fmt) 'Emittance', mode_m%a%emittance, mode_d%a%emittance, mode_m%b%emittance, mode_d%b%emittance, '! Unnormalized'
    nl=nl+1; write(lines(nl), '(a43, 2x, 2es13.5)') 'Emit (photon vert opening angle ignored)', &
                                                                    mode_m%b%emittance_no_vert, mode_d%b%emittance_no_vert
    if (mode_m%b%alpha_damp /= 0) then
      nl=nl+1; write(lines(nl), fmt) 'Alpha_damp', mode_m%a%alpha_damp, mode_d%a%alpha_damp, mode_m%b%alpha_damp, mode_d%b%alpha_damp, '! Damping per turn'
      nl=nl+1; write(lines(nl), fmt) 'Damping_time', time1/mode_m%a%alpha_damp, time1/mode_d%a%alpha_damp, time1/mode_m%b%alpha_damp, time1/mode_d%b%alpha_damp, '! Sec'
    endif
  endif

  !

  fmt  = '(1x, a16, 2es13.5, 3x, a)'
  fmt2 = '(1x, a16, 2f13.7, 3x, a)'

  if (branch%param%geometry == closed$) then
    nl=nl+1; lines(nl) = ''
    nl=nl+1; write(lines(nl), '(23x, a)') '  Model       Design'

    call calc_z_tune(branch)
    if (abs(design_lat%z%tune/twopi)  > 1d-3 .and. abs(branch%z%tune/twopi) > 1d-3) then
      nl=nl+1; write(lines(nl), fmt) 'Z_tune:', -branch%z%tune/twopi, -design_lat%z%tune/twopi
    else
      if (.not. branch%z%stable) then
        str1 = '  Unstable'
      else
        write (str1, '(f13.7)') -branch%z%tune/twopi
      endif
      if (.not. design_lat%z%stable) then
        str2 = '  Unstable'
      else
        write (str2, '(f13.7)') -design_lat%z%tune/twopi
      endif
      nl=nl+1; write(lines(nl), '(1x, a16, 2a13, 3x, a)') 'Z_tune:', str1, str2, '! The design value is calculated with RF on'
    endif

    if (mode_m%z%alpha_damp /= 0) then
      nl=nl+1; write(lines(nl), fmt) 'Sig_E/E:', mode_m%sigE_E, mode_d%sigE_E
      nl=nl+1; write(lines(nl), fmt) 'Sig_z:', mode_m%sig_z, mode_d%sig_z, '! Only calculated when RF is on'
      nl=nl+1; write(lines(nl), fmt) 'Emittance_z:', mode_m%z%emittance, mode_d%z%emittance, '! Only calculated when RF is on'
      nl=nl+1; write(lines(nl), fmt) 'Energy Loss:', mode_m%e_loss, mode_d%e_loss, '! Energy_Loss (eV / Turn)'
      nl=nl+1; write(lines(nl), fmt) 'J_damp:', mode_m%z%j_damp, mode_d%z%j_damp, '! Longitudinal Damping Partition #'
      nl=nl+1; write(lines(nl), fmt) 'Alpha_damp:', mode_m%z%alpha_damp, mode_d%z%alpha_damp, '! Longitudinal Damping per turn'
      nl=nl+1; write(lines(nl), fmt) 'damp_time:', time1/mode_m%z%alpha_damp, time1/mode_d%z%alpha_damp, '! Longitudinal Damping time (sec)'
    endif

    nl=nl+1; write(lines(nl), fmt) 'Alpha_p:', mode_m%momentum_compaction, mode_d%momentum_compaction, '! Momentum Compaction'
    nl=nl+1; write(lines(nl), fmt) 'Eta_p:', mode_m%momentum_compaction - 1.0_rp/gamma**2, mode_d%momentum_compaction - 1.0_rp/gamma**2, '! Slip factor'
    if (mode_m%momentum_compaction > 0) then
      nl=nl+1; write(lines(nl), fmt) 'gamma_trans:', sqrt(1.0_rp/mode_m%momentum_compaction), sqrt(1.0_rp/mode_d%momentum_compaction), '! Gamma at transition'
    endif

    if (bmad_com%spin_tracking_on) then
      nl=nl+1; write(lines(nl), fmt) 'Spin Tune:', branch%param%spin_tune/twopi, design_branch%param%spin_tune/twopi, '! Spin Tune on Closed Orbit (Units of 2pi)'
    endif

    pz1 = 0;  pz2 = 0
    do i = 1, branch%n_ele_track
      pz1 = pz1 + branch%ele(i)%value(l$) * (tao_branch%orbit(i-1)%vec(6) + tao_branch%orbit(i)%vec(6)) / (2 * branch%param%total_length)
      pz2 = pz2 + branch%ele(i)%value(l$) * (design_tao_branch%orbit(i-1)%vec(6) + design_tao_branch%orbit(i)%vec(6)) / (2 * branch%param%total_length)
    enddo
    nl=nl+1; write(lines(nl), fmt) '<pz>:', pz1, pz2, '! Average closed orbit pz (momentum deviation)'
  endif

  nl=nl+1; lines(nl) = ''
  nl=nl+1; lines(nl) = '[Note: The emittance, energy loss, etc. calcs have changed. Use "help show universe" for more info.]'

!----------------------------------------------------------------------
! variable
    
case ('use')  

  nl=nl+1; lines(nl) = 'veto data *@*'

  do i = lbound(s%u, 1), ubound(s%u, 1)
    do j = 1, s%u(i)%n_d2_data_used
      d2_ptr => s%u(i)%d2_data(j)
      call re_allocate (lines, nl+size(d2_ptr%d1)+10, .false.)
      do k = lbound(d2_ptr%d1, 1), ubound(d2_ptr%d1, 1)
        d1_ptr => d2_ptr%d1(k)
        call location_encode(line, d1_ptr%d%useit_opt, d1_ptr%d%exists, lbound(d1_ptr%d, 1), err_flag = err)
        if (err .or. len_trim(line) + 50 > len(line)) then
          nl=nl+1; lines(nl) = 'veto data ' // trim(d1_ptr%name)
          do n = lbound(d1_ptr%d, 1), ubound(d1_ptr%d, 1)
            if (.not. d1_ptr%d(n)%useit_opt) cycle
            if (nl + 100 > size(lines)) call re_allocate(lines, nl+100, .false.)
            nl=nl+1; write(lines(nl), '(3a, i0, a)') 'restore data ', trim(tao_d2_d1_name(d1_ptr)), '[', n, ']'
          enddo
        else
          if (line == '') cycle
          nl=nl+1; write(lines(nl), '(5a)') 'use data ', trim(tao_d2_d1_name(d1_ptr)), '[', trim(line), ']'
        endif
      enddo
    enddo
  enddo
  nl=nl+1; lines(nl) = ''

  call re_allocate (lines, nl+s%n_v1_var_used+10, .false.)
  do i = 1, s%n_v1_var_used
    v1_ptr => s%v1_var(i)
    if (v1_ptr%name == ' ') cycle
    call re_allocate (lines, nl+200, .false.)
    call location_encode (line, v1_ptr%v%useit_opt, v1_ptr%v%exists, lbound(v1_ptr%v, 1), err_flag = err)
    if (err .or. len_trim(line) + 50 > len(line)) then
      nl=nl+1; lines(nl) = 'veto var ' // trim(v1_ptr%name)
      do j = lbound(v1_ptr%v, 1), ubound(v1_ptr%v, 1)
        if (.not. v1_ptr%v(j)%useit_opt) cycle
        if (nl + 100 > size(lines)) call re_allocate(lines, nl+100, .false.)
        nl=nl+1; write(lines(nl), '(3a, i0, a)') 'restore var ', trim(v1_ptr%name), '[', j, ']'
      enddo
    else
      nl=nl+1; write(lines(nl), '(5a)') 'use var ', trim(v1_ptr%name), '[', trim(line), ']'
    endif
  enddo

!----------------------------------------------------------------------
! value

case ('value')

  s_fmt = 'es24.16'
  line = ''

  do
    call tao_next_switch (what2, [character(20):: '#format'], .true., line2, err, print_err = .false.)
    if (err) return
    select case (line2)
    case ('#format')
      call tao_next_word(what2, s_fmt)
    case default
      line = trim(line) // trim(line2)
      if (what2 == '') exit
    end select
  enddo


  map_calc = u%calc%one_turn_map

  if (index(line, 'chrom') /= 0 .or. index(line, 'rad') /= 0) then
    if (index(line, 'chrom') /= 0) s%com%force_chrom_calc = .true.
    if (index(line, 'rad') /= 0) s%com%force_rad_int_calc = .true.
    s%u%calc%lattice = .true.
    if (index(line, 'ptc') /= 0) u%calc%one_turn_map = .true.
    call tao_lattice_calc(ok)
  endif

  call tao_evaluate_expression (line, 0, .false., value, err, .true., info)
  u%calc%one_turn_map = map_calc
  if (err) then
    if (index(what, '-f') /= 0) call out_io (s_warn$, r_name, &
                '"Note: -format" has been changed to "#format" to avoid confusion with any negative signs.')
    return
  endif

  if (size(value) == 1) then
    s_fmt = '(3x, ' // trim(s_fmt) // ')'
    nl=nl+1; write(lines(nl), s_fmt, iostat = ios) value(1)
    ! For some funny reason ios can be zero on a bad format so check for a star in the string.
    if (ios /= 0 .or. index(lines(nl), '*') /= 0) then
      write(lines(nl), '(3x, es24.16, a)') value(1), '  ! Note: Value/format mismatch detected'
    endif

  else
    if (is_integer(s_fmt(1:1)) .or. index(s_fmt, ',') /= 0) then   ! print value on a single line
      s_fmt = '(5x, ' // trim(s_fmt) // ')'    
      nl=nl+1; write(lines(nl), s_fmt, iostat = ios) value

    else
      s_fmt = '(i5, a, ' // trim(s_fmt) // ', 5x, a)'    
      call re_allocate (lines, size(value)+100, .false.)
      do i = 1, size(value)
        if (associated(info(i)%ele)) then
          nl=nl+1; write(lines(nl), s_fmt, iostat = ios) i, ':  ', value(i), ele_full_name(info(i)%ele)
        else
          nl=nl+1; write(lines(nl), s_fmt, iostat = ios) i, ':  ', value(i)
        endif
        ! For some funny reason ios can be zero on a bad format so check for a star in the string.
        if (ios /= 0 .or. index(lines(nl), '*') /= 0) then
          write(lines(nl), '(i5, a, es24.16, a)', iostat = ios) i, ':  ', value(i), '  ! Note: Value/format mismatch detected'
        endif
      enddo
    endif
  endif

!----------------------------------------------------------------------
! variable

case ('variables')

  good_opt_only = .false.
  bmad_format = .false.
  print_header_lines = .true.
  print_by_uni = .false.
  attrib0 = ''

  do
    call tao_next_switch (what2, [character(16):: '-bmad_format', '-good_opt_only', & 
                                                   '-no_label_lines', '-universe'], .true., switch, err)
    if (err) return

    select case (switch)  
    case ('') 
      exit
    case ('-bmad_format') 
      bmad_format = .true.
    case ('-good_opt_only') 
      good_opt_only = .true.
    case ('-no_label_lines')
      print_header_lines = .false.
    case ('-universe')
      call tao_next_word(what2, aname)
      call tao_pick_universe (trim(aname) // '@',  str, picked_uni, err)
      if (err) return
      print_by_uni = .true.
    case default
      if (attrib0 /= '') then
        call out_io (s_error$, r_name, 'EXTRA STUFF ON LINE: ' // attrib0)
        return
      endif
      attrib0 = switch
    end select

  enddo

  if (s%n_v1_var_used == 0) then
    nl=1; lines(1) = 'NO VARIABLES HAVE BEEN DEFINED IN THE INPUT FILES!'
    return 
  endif

  ! If 'n@' is present then write out stuff for universe n

  if (print_by_uni) then
    do ix_u = 1, size(s%u)
      if (.not. picked_uni(ix_u)) cycle

      if (print_header_lines) then
        nl=nl+1; lines(nl) = ''
        nl=nl+1; write(lines(nl), '(a, i4)') 'Variables controlling universe:', ix_u
        nl=nl+1; write(lines(nl), '(5x, a)') 'Name:'
      endif

      do i = 1, s%n_var_used
        if (.not. s%var(i)%exists) cycle
        found = .false.
        do j = 1, size(s%var(i)%slave)
          if (s%var(i)%slave(j)%ix_uni == ix_u) found = .true.
        enddo
        if (.not. found) cycle
        if (nl+10 > size(lines)) call re_allocate(lines, nl+100, .false.)
        nl=nl+1; write(lines(nl), '(5x, a25, a40)') tao_var1_name(s%var(i)), tao_var_attrib_name(s%var(i))
      enddo
    enddo

    result_id = 'variable:@'
    return
  endif

  ! If just "show var" then show all names

  if (attrib0 == '') then
    ! Bmad format
    if (bmad_format) then
      call tao_print_vars (0, 0, good_opt_only)
      result_id = 'variable:bmad'
      return
    endif

    if (print_header_lines) then
      nl=nl+1; write(lines(nl), '(7x, a, t50, a)') 'Name', 'Using for Optimization'
    endif
    do i = 1, s%n_v1_var_used
      v1_ptr => s%v1_var(i)
      if (v1_ptr%name == ' ') cycle
      call re_allocate (lines, nl+200, .false.)
      call location_encode (line, v1_ptr%v%useit_opt, v1_ptr%v%exists, lbound(v1_ptr%v, 1))
      nl=nl+1; write(lines(nl), '(4x, 2a, i0, a, i0, a, t50, a)') &
                      trim(v1_ptr%name), '[', lbound(v1_ptr%v, 1), ':', &
                      ubound(v1_ptr%v, 1), ']', trim(line)
      
    enddo

    result_id = 'variable:'
    return
  endif

  ! are we looking at a range of locations?

  call tao_find_var(err, attrib0, v1_array, v_array) 
  if (err) return
  n_size = 0
  if (allocated(v_array)) n_size = size(v_array)

  ! Bmad format
  if (bmad_format) then
    call tao_print_vars (0, 0, good_opt_only, .false., v_array)
    result_id = 'variable:bmad'
    return
  endif

  ! v_ptr is valid then show the variable info.

  if (n_size == 1) then
    v_ptr => v_array(1)%v
    if (.not. v_ptr%exists) then
      nl=nl+1; write(lines(nl), *) 'Variable does not exist: ' // tao_var1_name(v_ptr)
      return
    endif

    nl=nl+1; write(lines(nl), amt)  '%ele_name         = ', quote(v_ptr%ele_name)
    nl=nl+1; write(lines(nl), amt)  '%attrib_name      = ', quote(v_ptr%attrib_name)
    if (v_ptr%id /= '') then
      nl=nl+1; write(lines(nl), amt)  '%id                = ', quote(v_ptr%id)
    endif
    nl=nl+1; write(lines(nl), imt)  '%ix_attrib        = ', v_ptr%ix_attrib 
    nl=nl+1; write(lines(nl), imt)  '%ix_var           = ', v_ptr%ix_var
    nl=nl+1; write(lines(nl), imt)  '%ix_dvar          = ', v_ptr%ix_dvar           
    nl=nl+1; write(lines(nl), imt)  '%ix_v1            = ', v_ptr%ix_v1
    nl=nl+1; write(lines(nl), rmt)  '%model            = ', v_ptr%model_value
    nl=nl+1; write(lines(nl), rmt)  '%base             = ', v_ptr%base_value

    if (.not. allocated (v_ptr%slave)) then
      nl=nl+1; write(lines(nl), imt)  'this(:) -- Not associated!'
    else
      n = nl + 3*size(v_ptr%slave) + 100
      if (size(lines) < n) call re_allocate(lines, n)
      do i = 1, size(v_ptr%slave)
        if (v_ptr%slave(i)%ix_ele < 0) then  ! Happens for particle_start
          nl=nl+1; write(lines(nl), '(2(a, i0), 2a)')  '%slave(', i, '): Slave Element: ', &
                                                          v_ptr%slave(i)%ix_uni, '@', 'PARTICLE_START'
        else
          ele => s%u(v_ptr%slave(i)%ix_uni)%model%lat%branch(v_ptr%slave(i)%ix_branch)%ele(v_ptr%slave(i)%ix_ele)
          nl=nl+1; write(lines(nl), '(2(a, i0), 2a)')  '%slave(', i, '): Slave Element: ', &
                                                          v_ptr%slave(i)%ix_uni, '@', ele_full_name(ele)
        endif

        if (associated (v_ptr%slave(i)%model_value)) then
          nl=nl+1; write(lines(nl), irmt)  '%slave(', i, ')%Model_value: ', &
                                                            v_ptr%slave(i)%model_value
        else
          nl=nl+1; write(lines(nl), irmt)  '%slave(', i, ')%Model_value: <not associated>'
        endif
        if (associated (v_ptr%slave(i)%base_value)) then
          nl=nl+1; write(lines(nl), irmt)  '%slave(', i, ')%Base_value:  ', &
                                                            v_ptr%slave(i)%base_value
        else
          nl=nl+1; write(lines(nl), irmt)  '%slave(', i, ')%Base_value:  <not associated>'
        endif
      enddo
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
    nl=nl+1; write(lines(nl), amt)  '%merit_type       = ', quote(v_ptr%merit_type)
    nl=nl+1; write(lines(nl), rmt)  '%merit            = ', v_ptr%merit
    nl=nl+1; write(lines(nl), rmt)  '%dmerit_dvar      = ', v_ptr%dMerit_dVar
    nl=nl+1; write(lines(nl), rmt)  '%s                = ', v_ptr%s
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
    nl=nl+1; write(lines(nl), lmt)  '%key_bound        = ', v_ptr%key_bound

    result_id = 'variable:1:' // attrib0

  ! check if there is a variable number
  ! if no variable number requested, show a range

  elseif (n_size > 1) then

    nc = 40
    do i = 1, size(v_array)
      v_ptr => v_array(i)%v
      if (.not. v_ptr%exists) cycle
      nc = max(nc, len_trim(tao_var_attrib_name(v_ptr)))
    enddo

    if (print_header_lines) then
      write (line1, '(a, t30, a)') '  Variable', 'Slave Parameters'
      line1(nc+17:) = 'Meas         Model        Design  Useit_opt'
      nl=nl+1; lines(nl) = line1
    endif

    ! if a range is specified, show the variable range   
    do i = 1, size(v_array)
      v_ptr => v_array(i)%v
      if (.not. v_ptr%exists) cycle
      call re_allocate (lines, nl+200, .false.)
      nl=nl+1
      write(lines(nl), '(2x, 2a, i0, a, t30, a)') trim(v_ptr%v1%name), '[', v_ptr%ix_v1, ']', tao_var_attrib_name(v_ptr)
      write(lines(nl)(nc+9:), '(3es14.4, 7x, l)') v_ptr%meas_value, v_ptr%model_value, v_ptr%design_value, v_ptr%useit_opt
    enddo

    if (print_header_lines) then
      nl=nl+1; lines(nl) = line1
    endif

  else
    nl=1; lines(1) = 'Cannot find variables matching: ' // attrib0
    result_id = 'variable:?'
  endif

!----------------------------------------------------------------------
! version

case ('version')

  nl=nl+1; lines(nl) = 'Date: ' // tao_version_date

!----------------------------------------------------------------------
! wake_elements

case ('wake_elements')

  nl=nl+1; write(lines(nl), lmt) '  bmad_com%sr_wakes_on = ', bmad_com%sr_wakes_on
  nl=nl+1; write(lines(nl), lmt) '  bmad_com%lr_wakes_on = ', bmad_com%lr_wakes_on
  nl=nl+1; lines(nl) = ''
  nl=nl+1; lines(nl) = '                                           Short-Range    Short-Range   Long-'
  nl=nl+1; lines(nl) = 'Index  Element               Key           Longitudinal   Transverse    Range'
  do ib = 0, ubound(lat%branch, 1)
    branch => lat%branch(ib)

    do ie = 1, branch%n_ele_max
      ele => branch%ele(ie)

      if (ele%slave_status == super_slave$) cycle
      if (ele%slave_status == multipass_slave$) cycle
      if (.not. associated(ele%wake)) cycle
      wake => ele%wake

      nl=nl+1; write(lines(nl), '(a5, 2x, a20, 2x, a15, 3x, l1, 13x, l1, 13x, l1)') &
        ele_loc_name(ele, .true.), ele%name, key_name(ele%key), &
        allocated(wake%sr%long), allocated(wake%sr%trans), allocated(wake%lr%mode)
    enddo
  enddo

!----------------------------------------------------------------------
! wall

case ('wall')

  by_s = .false.
  ix_sec = -1
  angle = 0
  sub_name = ''
  attrib0 = ''


  do
    call tao_next_switch (what2, [character(16):: '-section', '-element', &
                                   '-angle', '-s', '-branch'], .true., switch, err)
    if (err) return
    if (switch == '') exit
    select case (switch)

    case ('-angle')
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios) angle
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ ANGLE.'
        return
      endif

    case ('-branch')
      call tao_next_word(what2, aname)
      branch => pointer_to_branch(aname, u%model%lat)
      if (.not. associated(branch)) then
        nl=1; write(lines(1), *) 'Bad branch index:', ix_branch
        return
      endif
      ix_branch = branch%ix_branch

    case ('-section')
      call tao_next_word(what2, aname)
      read (aname, *, iostat = ios) ix_sec
      if (ios /= 0) then
        nl=1; lines(1) = 'CANNOT READ SECTION INDEX.'
        return
      endif

    case ('-s')
      by_s = .true.

    case default
      attrib0 = trim(attrib0) // ' ' // trim(switch)
    end select
  enddo

  !-------

  if (.not. associated(branch%wall3d)) then
    nl=1; lines(nl) = 'No associated vacuum chamber wall.'
    result_id = 'wall:none'
    return
  endif

  wall => branch%wall3d(1)

  if (ix_sec > 0) then 
    if (ix_sec > size(wall%section)) then
      nl=1; write(lines(nl), '(a, i0)') 'Section index larger than number of sections: ', size(wall%section)
      result_id = 'wall:sec:large'
      return
    endif

    wall_sec => wall%section(ix_sec)
    ele => pointer_to_ele(lat, wall_sec%ix_ele, wall_sec%ix_branch)
    nl=nl+1; write(lines(nl), '(5a)')            'ele:    ', trim(ele%name), '   (', trim(ele_loc_name(ele)), ')'
    nl=nl+1; write(lines(nl), '(2a)')            'type:   ', trim(wall3d_section_type_name(wall_sec%type))
    nl=nl+1; write(lines(nl), '(a, f14.6)')      'S:      ', wall_sec%s
    nl=nl+1; write(lines(nl), '(3(a, f10.6))')  ' r0:     (', wall_sec%r0(1), ',', wall_sec%r0(2), ')'
    if (wall_sec%dr_ds == real_garbage$) then
      nl=nl+1; write(lines(nl), '(3(a, f10.6))')  'dr_ds:       Not-Set'
    else
      nl=nl+1; write(lines(nl), '(3(a, f10.6))')  'dr_ds: ', wall_sec%dr_ds
    endif

    do j = 1, size(wall_sec%v)
      v => wall_sec%v(j)
      nl=nl+1; write(lines(nl), '(a, i0, a, 5f11.6)') &
                            'v(', j, ') =', v%x, v%y, v%radius_x, v%radius_y, v%tilt
    enddo

    return

  endif

  !-------

  if (by_s) then
    ix2 = index(attrib0, ':')
    if (ix2 == 0) then
      nl=1; lines(nl) = 'NO ":" FOUND FOR RANGE SELECTION'
      return
    endif
    read (attrib0(1:ix2-1), *, iostat = ios1) s1
    read (attrib0(ix2+1:), *, iostat = ios2) s2
    if (ios1 /= 0 .or. ios2 /= 0) then
      nl=1; lines(1) = 'ERROR READING RANGE SELECTION: ' // attrib0
      return
    endif

    ix1 = bracket_index (s1 - 1d-10, wall%section%s, 1)
    ix2 = bracket_index (s2 + 1d-10, wall%section%s, 1) 
    ix1 = ix1 + 1

  elseif (attrib0 /= '') then
    ix2 = index(attrib0, ':')
    if (ix2 == 0) then
      nl=1; lines(nl) = 'NO ":" FOUND FOR RANGE SELECTION'
      return
    endif
    read (attrib0(1:ix2-1), *, iostat = ios1) ix1
    read (attrib0(ix2+1:), *, iostat = ios2) ix2
    if (ios1 /= 0 .or. ios2 /= 0) then
      nl=1; lines(1) = 'ERROR READING RANGE SELECTION: ' // attrib0
      return
    endif

  else
    ix1 = 1; ix2 = min(200, size(wall%section))
  endif

  nl=nl+1; lines(nl) = '    Ix             S    ix_ele                 Ele          Type   Radius (mm)'

  do i = ix1, ix2
    wall_sec => wall%section(i)
    ele => pointer_to_ele (lat, wall_sec%ix_ele, wall_sec%ix_branch)
    
    call calc_wall_radius (wall%section(i)%v, cos(angle), sin(angle), r, z)
    nl=nl+1; write(lines(nl), '(i6, f14.6, a10, a20, a14, f14.3)') i, wall_sec%s, &
                trim(ele_loc_name(ele)), trim(ele%name), trim(wall3d_section_type_name(wall_sec%type)), 1000*r
  enddo

  nl=nl+1; lines(nl) = '    Ix             S    ix_ele                 Ele          Type   Radius (mm)'

!----------------------------------------------------------------------
! wave

case ('wave')

  nl=nl+1; write(lines(nl), '(a, 2i5)') 'ix_a:', s%wave%ix_a1, s%wave%ix_a2
  nl=nl+1; write(lines(nl), '(a, 2i5)') 'ix_b:', s%wave%ix_b1, s%wave%ix_b2

  select case (s%wave%data_type)
  case ('orbit.x', 'orbit.y', 'eta.x', 'eta.y', 'beta.a', 'beta.b', 'ping_a.amp_x', 'ping_b.amp_y')
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
    nl=nl+1; lines(nl) = 'After Dat#    Norm_Kick     s-pos   ele@kick                                 phi'
    do i = 1, min(s%wave%n_kick, 20)
      wk => s%wave%kick(i)
      nl=nl+1; write(lines(nl), '(i9, f14.2, f10.2, 3x, a6, a30, f10.3)') wk%ix_dat_before_kick, 1e6*wk%amp, &
                                                          wk%s, ele_loc_name(wk%ele), wk%ele%name, wk%phi
    enddo

  case ('phase.a', 'phase.b', 'ping_a.phase_x', 'ping_b.phase_y')
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'A Region Sigma_Fit/Amp_Fit:  ', s%wave%rms_rel_a
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'B Region Sigma_Fit/Amp_Fit:  ', s%wave%rms_rel_b
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_Kick/Kick: ', s%wave%rms_rel_k
    nl=nl+1; write(lines(nl), '(a, f8.3)') 'Sigma_phi:       ', s%wave%rms_phi
    nl=nl+1; write(lines(nl), '(a, f8.3, a)') &
                                    'Chi_C:           ', s%wave%chi_c, ' [Figure of Merit]'
    nl=nl+1; lines(nl) = ' '
    nl=nl+1; lines(nl) = 'Normalized Kick = k * l * beta [dimensionless]'
    nl=nl+1; lines(nl) = '   where k = quadrupole gradient [rad/m^2].'
    nl=nl+1; lines(nl) = 'After Dat#    Norm_Kick     s-pos   ele@kick                                 phi'
    do i = 1, min(s%wave%n_kick, 20)
      wk => s%wave%kick(i)
      nl=nl+1; write(lines(nl), '(i9, f14.2, f10.2, 3x, a6, a30, f10.3)') wk%ix_dat_before_kick, wk%amp, &
                                                          wk%s, ele_loc_name(wk%ele), wk%ele%name, wk%phi
    enddo

  case ('ping_a.amp_sin_rel_y', 'ping_a.amp_cos_rel_y', 'ping_b.amp_sin_rel_x', 'ping_b.amp_cos_rel_x', &
        'ping_a.amp_sin_y', 'ping_a.amp_cos_y', 'ping_b.amp_sin_x', 'ping_b.amp_cos_x', 'cbar.11', 'cbar.12', 'cbar.22')
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

    nl=nl+1; lines(nl) = 'After Dat#   Norm_Kick     s-pos   ele@kick                         phi+    phi-   phi_a   phi_b'
    do i = 1, min(s%wave%n_kick, 20)
      wk => s%wave%kick(i)
      nl=nl+1; write(lines(nl), '(i9, f10.4, f10.2, 3x, a6, a30, 4f8.3)') wk%ix_dat_before_kick, &
            wk%amp, wk%s, ele_loc_name(wk%ele), wk%ele%name, wk%phi_s, wk%phi_r, (wk%phi_s+wk%phi_r)/2, (wk%phi_s-wk%phi_r)/2
    enddo

  end select

  if (s%wave%n_kick > 20) then
    nl=nl+1; lines(nl) = ' etc...'
  endif

!----------------------------------------------------------------------

case default

  nl=1; lines(1) = "INTERNAL ERROR, SHOULDN'T BE HERE!"

end select

!----------------------------------------------------------------------
!----------------------------------------------------------------------
contains

subroutine point_to_this_uni(ix_uni, ix_branch, u, lat, branch, model_branch, tao_branch, design_tao_branch)

type (tao_universe_struct), pointer :: u
type (tao_lattice_branch_struct), pointer :: tao_branch, design_tao_branch
type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (tao_model_branch_struct), pointer :: model_branch

integer ix_uni, ix_branch

!

u => tao_pointer_to_universe(ix_uni)
if (.not. associated(u)) return

lat => u%model%lat
branch => lat%branch(ix_branch)
model_branch => u%model_branch(ix_branch)
tao_branch => u%model%tao_branch(ix_branch)
design_tao_branch => u%design%tao_branch(ix_branch)

end subroutine point_to_this_uni

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! contains

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

nl=nl+1; write(lines(nl), '(a)') "  "
write (l1, '(a, 20x, a)') "Data Name", &
          "Data Type             |  Model Value  |  Design Value |  Base Value"
nl=nl+1; lines(nl) = l1

found_one = .false.
do i = 1, size(u%data)
  if (u%data(i)%ix_ele == ele%ix_ele .and. u%data(i)%ix_branch == ele%ix_branch) then
    found_one = .true.
    datum => u%data(i)
    nl=nl+1; write(lines(nl), "(a, t30, a20, 3(1x, es15.5))") &
                trim(tao_datum_name(datum)),  datum%data_type, datum%model_value, &
                datum%design_value, datum%base_value 
  endif
enddo

if (found_one) then
  nl=nl+1; lines(nl) = l1
else
  write(lines(nl), '(a)') "No data associated with this element."
endif

end subroutine show_ele_data

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! contains

subroutine show_opt ()

use geodesic_lm

implicit none

character(200), allocatable :: alloc_lines(:)

!

nl=nl+1; lines(nl) = 'Global optimization parameters (use "set global" to change):'
nl=nl+1; write(lines(nl), imt) '  %datum_err_messages_max        = ', s%global%datum_err_messages_max
nl=nl+1; write(lines(nl), rmt) '  %de_lm_step_ratio              = ', s%global%de_lm_step_ratio
nl=nl+1; write(lines(nl), rmt) '  %de_var_to_population_factor   = ', s%global%de_var_to_population_factor
nl=nl+1; write(lines(nl), rmt) '  %lm_opt_deriv_reinit           = ', s%global%lm_opt_deriv_reinit
nl=nl+1; write(lines(nl), rmt) '  %lmdif_eps                     = ', s%global%lmdif_eps
nl=nl+1; write(lines(nl), rmt) '  %lmdif_negligible_merit        = ', s%global%lmdif_negligible_merit
nl=nl+1; write(lines(nl), rmt) '  %merit_stop_value              = ', s%global%merit_stop_value
nl=nl+1; write(lines(nl), rmt) '  %dmerit_stop_value             = ', s%global%dmerit_stop_value
nl=nl+1; write(lines(nl), rmt) '  %svd_cutoff                    = ', s%global%svd_cutoff
nl=nl+1; write(lines(nl), rmt) '  %unstable_penalty              = ', s%global%unstable_penalty

nl=nl+1; write(lines(nl), imt) '  %n_top10_merit                 = ', s%global%n_top10_merit
nl=nl+1; write(lines(nl), imt) '  %n_opti_loops                  = ', s%global%n_opti_loops
nl=nl+1; write(lines(nl), imt) '  %n_opti_cycles                 = ', s%global%n_opti_cycles

nl=nl+1; write(lines(nl), amt) '  %optimizer                     = ', quote(s%global%optimizer)
nl=nl+1; write(lines(nl), amt) '  %var_out_file                  = ', quote(s%global%var_out_file)
nl=nl+1; write(lines(nl), lmt) '  %opt_match_auto_recalc         = ', s%global%opt_match_auto_recalc
nl=nl+1; write(lines(nl), lmt) '  %opti_write_var_file           = ', s%global%opti_write_var_file
nl=nl+1; write(lines(nl), lmt) '  %derivative_recalc             = ', s%global%derivative_recalc
nl=nl+1; write(lines(nl), lmt) '  %derivative_uses_design        = ', s%global%derivative_uses_design
nl=nl+1; write(lines(nl), lmt) '  %opt_with_ref                  = ', s%global%opt_with_ref 
nl=nl+1; write(lines(nl), lmt) '  %opt_with_base                 = ', s%global%opt_with_base
nl=nl+1; write(lines(nl), lmt) '  %optimizer_allow_user_abort    = ', s%global%optimizer_allow_user_abort
nl=nl+1; write(lines(nl), lmt) '  %optimizer_var_limit_warn      = ', s%global%optimizer_var_limit_warn
nl=nl+1; write(lines(nl), lmt) '  %svd_retreat_on_merit_increase = ', s%global%svd_retreat_on_merit_increase 
nl=nl+1; write(lines(nl), lmt) '  %var_limits_on                 = ', s%global%var_limits_on
nl=nl+1; write(lines(nl), lmt) '  %only_limit_opt_vars           = ', s%global%only_limit_opt_vars
nl=nl+1; lines(nl) = ''
nl=nl+1; lines(nl) = 'opti_de_param Parameters:'
nl=nl+1; write(lines(nl), rmt) '  %CR                   = ', opti_de_param%CR
nl=nl+1; write(lines(nl), rmt) '  %F                    = ', opti_de_param%F
nl=nl+1; write(lines(nl), rmt) '  %l_best               = ', opti_de_param%l_best
nl=nl+1; write(lines(nl), lmt) '  %binomial_cross       = ', opti_de_param%binomial_cross
nl=nl+1; write(lines(nl), lmt) '  %use_2nd_diff         = ', opti_de_param%use_2nd_diff
nl=nl+1; write(lines(nl), lmt) '  %randomize_F          = ', opti_de_param%randomize_F
nl=nl+1; write(lines(nl), lmt) '  %minimize_merit       = ', opti_de_param%minimize_merit
nl=nl+1; lines(nl) = ''
nl=nl+1; lines(nl) = 'geodesic_lm Parameters:'
call type_geodesic_lm (lines = alloc_lines, n_lines = n, prefix = '  %')
lines(nl+1:nl+n) = alloc_lines(1:n)
nl = nl + n

end subroutine show_opt

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! contains

subroutine show_ping(ix_uni)

type (tao_universe_struct), pointer :: u
integer ix_uni

u => tao_pointer_to_universe(ix_uni)
if (.not. associated(u)) return
nl=nl+1; lines(nl) = ''
nl=nl+1; write(lines(nl), rmt) 'ping_scale%a_mode_meas = ', u%ping_scale%a_mode_meas
nl=nl+1; write(lines(nl), rmt) 'ping_scale%a_mode_ref  = ', u%ping_scale%a_mode_ref
nl=nl+1; write(lines(nl), rmt) 'ping_scale%b_mode_meas = ', u%ping_scale%b_mode_meas
nl=nl+1; write(lines(nl), rmt) 'ping_scale%b_mode_ref  = ', u%ping_scale%b_mode_ref

end subroutine show_ping

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
! contains

function get_this_track_fmt (string, dflt_fmt, err) result (fmt)

character(*) string, dflt_fmt
character(16) fmt
integer ix
logical err

!

call string_trim(string, string, ix)
if (string(1:1) == '-' .or. ix == 0) then
  fmt = dflt_fmt
else
  fmt = string(1:ix)
  call string_trim(string(ix+1:), string, ix)
endif

err = .false.

end function get_this_track_fmt

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
! contains

subroutine write_track_header (line, ix_line, fmt, label, err)

implicit none

character(*) line, fmt, label(:)
character(16) code, str
integer ix_line, i, n, multiplyer, power, width, digits, ix
logical err

!

err = .false.
if (fmt == '' .or. fmt == 'no') return

call parse_fortran_format (fmt, multiplyer, power, code, width, digits)
if (code == '') then
  nl = 1; lines(1) = 'BAD FORMAT: ' // fmt
  err = .true.
  return
endif

do i = 1, size(label)
  if (power == 0 .or. code /= 'F') then
    str = label(i)
  else
    write (str, '(2a, i0, a)') trim(label(i)), ' (*1E', power, ')'
  endif
  n = len_trim(str)
  ix = ix_line + max(1, width-n-2)
  line(ix:) = str
  ix_line = ix_line + width
enddo

ix_line = ix_line + 2

end subroutine write_track_header

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
! contains

subroutine write_track_info (line, ix_line, fmt, value, err)

implicit none

character(*) line, fmt
character(8) code
real(rp) value(:)
integer ix_line, i, n, multiplyer, power, width, digits, ios
logical err

!

err = .false.
if (fmt == '' .or. fmt == 'no') return

call parse_fortran_format (fmt, multiplyer, power, code, width, digits)

do i = 1, size(value)
  write (line(ix_line+1:), '(' // fmt // ')', iostat = ios) value(i)
  if (ios /= 0) then
    nl = 1; lines(1) = 'BAD FORMAT: ' // fmt
    err = .true.
    return
  endif
  ix_line = ix_line + width
enddo

ix_line = ix_line + 2

end subroutine write_track_info

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
! contains

subroutine write_this_arg (nl, lines, str, value)
integer nl
character(*), allocatable :: lines(:)
character(*) str, value
!
select case (value)
case ('');                              return
case (present_str, 'BLUE', 'all');      nl=nl+1; write(lines(nl), '(2x, a)') str
case (negated_str, 'off', 'DEFAULT');   nl=nl+1; write(lines(nl), '(2x, 2a)') '-' // str
case default;                           nl=nl+1; write(lines(nl), '(2x, a, 2x, a)') str, trim(value)
end select
end subroutine write_this_arg

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
! contains

subroutine set_this_show_lat_header (line2, line3, who, fmt, called_from_pipe_cmd, ix)

character(*) line2, line3, who, fmt
integer, optional :: ix
logical called_from_pipe_cmd

!

if (called_from_pipe_cmd) then
  if (line2 /= '') then
    line2 = trim(line2) // ';'
    line3 = trim(line3) // ';'
  endif

  line2 = trim(line2) // who

  if (index(fmt, 'A') /= 0) then
    line3 = trim(line3) // 'STR'
  elseif (index(fmt, 'I') /= 0) then
    line3 = trim(line3) // 'INT'
  elseif (index(fmt, 'L') /= 0) then
    line3 = trim(line3) // 'LOGIC'
  else
    line3 = trim(line3) // 'REAL'
  endif

else
  line2(ix:) = who
endif

end subroutine set_this_show_lat_header 

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
! contains

subroutine write_real (line, fmt, value)

implicit none

real(rp) value

integer id, ip, ix, wid, pl, wid_want, pl_want, ios

character(16) fmt2, num_str
character(*) line, fmt

!

write (line, fmt, iostat = ios) value
if (ios /= 0) then
  do ix = 1, len(line)
    line(ix:ix) = '*'
  enddo
  return
endif

if (value == 0) return
if (line(1:1) == ' ') return
if (fmt(2:2) /= 'f' .and. fmt(2:2) /= 'F') return

! Always have a blank as first character.

id = index(fmt, '.')
ip = index(fmt, ')')

read (fmt(3:id-1), *) wid
read (fmt(id+1:ip-1), *) pl

wid_want = log10((1 + 1/10.0**(wid-2)) * abs(value)) + 1
if (value < 0) wid_want = wid_want + 1

pl_want = wid - wid_want - 2

if (pl_want == 0) then
  write (fmt2, '(a, i0, a, i0, a)') '(1x, f', wid, '.', pl_want, ')'
  write (line, fmt2) value
  line(wid+1:wid+1) = ' ' ! Get rid of '.'  
  return
endif

if (pl_want > 0) then  ! Can use F format.
  write (fmt2, '(a, i0, a, i0, a)') '(1x, f', wid-1, '.', pl_want, ')'
  write (line, fmt2) value
  return
endif

! Number is too large so switch to ES format.

wid_want = 4
if (value < 0) wid_want = wid_want + 1
if (abs(value) > 0.99d10) wid_want = wid_want + 1
if (abs(value) > 0.99d100) wid_want = wid_want + 1

if (wid < wid_want) return  ! Cannot do anything

pl = wid - wid_want - 1
if (pl < 0) pl = 0
write (fmt2, '(a, i0, a, i0, a)') '(es16.', pl, ')'
write (num_str, fmt2) value

id = index(num_str, '+')  ! Get rid of "+" sign in "E+nnn"
if (id /= 0) num_str = num_str(1:id-1) // num_str(id+1:)

id = index(num_str, 'E0')  ! Get rid of "0" in "E005"
if (id /= 0) num_str = num_str(1:id) // num_str(id+2:)

id = index(num_str, 'E0')  ! Get rid of "0" in "E005"
if (id /= 0) num_str = num_str(1:id) // num_str(id+2:)

call string_trim (num_str, num_str, ix)
line = ''
line(wid-ix+1:) = num_str

end subroutine write_real

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
! contains

function setup_lat_column(name, format, label, remove_line_if_zero, scale_factor) result (col)

type (show_lat_column_struct) col
character(*) :: name, format, label
character(4) descrip
character(20) fmt
logical :: remove_line_if_zero
real(rp) :: scale_factor
integer :: width, multi, power, digits

!

fmt = upcase(format)
col = show_lat_column_struct(name, fmt, label, remove_line_if_zero, scale_factor, 0)

call parse_fortran_format (fmt, multi, power, descrip, col%width, digits)
call upcase_string(descrip)
if (descrip == 'X') col%width = max(1, multi, col%width) 
if (fmt == 'A') col%width = 0

end function setup_lat_column

end subroutine tao_show_this
