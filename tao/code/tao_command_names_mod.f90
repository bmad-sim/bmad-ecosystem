!+
! Module tao_command_names_mod
!
! Central definitions of the command name lists used for command parsing.
! Hoisted here (rather than being literals at the match_word call sites) so the
! tab completion engine (tao_completion_mod) can offer the same names the
! parsers accept. Keep each list in sync with the corresponding select case.
!-

module tao_command_names_mod

use tao_struct
use quick_plot
use attribute_mod, only: switch_attrib_value_name

implicit none

! Top level Tao commands. Matched case sensitively in tao_command.

character(16), parameter :: tao_command_names(49) = [character(16):: &
                      'alias', 'call', 'change', 'clear', 'clip', 'continue', 'create', 'cut_ring', 'derivative', &
                      'end_file', 'exit', 'fixer', 'flatten', 'help', 'json', 'ls', 'misalign', 'pause', 'pipe', 'place', &
                      'plot', 'ptc', 'python', 'quit', 're_execute', 'read', 'regression', 'reinitialize', 'reset', &
                      'restore', 'run_optimizer', 'scale', 'set', 'show', 'single_mode', 'spawn', 'taper', &
                      'timer', 'use', 'veto', 'view', 'wave', 'write', 'x_axis', 'x_scale', 'xy_scale', &
                      'debug', 'verbose', 'tree']

! Commands accepted by the parser but not documented in command-list.tex.
! Tab completion omits these; the parser list above must still include them.

character(16), parameter :: tao_hidden_command_names(3) = [character(16):: 'debug', 'verbose', 'tree']

! "show <what>" names. See tao_show_this.

character(20), parameter :: tao_show_what_names(48) = [character(20):: 'alias', 'beam', 'branch', 'building_wall', &
        'chromaticity', 'constraints', 'control', 'curve', 'data', 'debug', &
        'derivative', 'dynamic_aperture', 'element', 'emittance', 'field', 'global', 'graph', &
        'history', 'hom', 'internal', 'key_bindings', 'lattice', 'matrix', 'merit', 'normal_form', &
        'optimizer', 'orbit', 'particle', 'plot', 'ptc', 'radiation_integrals', 'rampers', 'spin', 'string', &
        'symbolic_numbers', 'taylor_map',  'top10', &
        'track', 'tune', 'twiss_and_orbit', 'universe', 'use', 'value', 'variables', 'version', &
        'wake_elements', 'wall', 'wave']

! "set <target>" names. See the set case in tao_command.

character(20), parameter :: tao_set_target_names(33) = [character(20) :: 'branch', 'data', 'variable', 'lattice', &
      'universe', 'curve', 'graph', 'beam_init', 'wave', 'plot', 'bmad_com', 'element', 'opti_de_param', &
      'csr_param', 'floor_plan', 'lat_layout', 'geodesic_lm', 'default', 'key', 'particle_start', &
      'plot_page', 'ran_state', 'symbolic_number', 'beam', 'beam_start', 'dynamic_aperture', &
      'global', 'region', 'calculate', 'space_charge_com', 'ptc_com', 'tune', 'z_tune']

! "pipe <subcommand>" names. See tao_pipe_cmd.

character(40), parameter :: tao_pipe_cmd_names(114) = [character(40) :: &
          'beam', 'beam_init', 'branch1', 'bunch_comb', 'bunch_params', 'bunch1', 'bmad_com',&
          'building_wall_list', 'building_wall_graph', 'building_wall_point', 'building_wall_section', &
          'complete', 'constraints', 'da_params', 'da_aperture', &
          'data', 'data_d2_create', 'data_d2_destroy', 'data_d_array', 'data_d1_array', &
          'data_d2', 'data_d2_array', 'data_set_design_value', 'data_parameter', &
          'datum_create', 'datum_has_ele', 'derivative', &
          'ele:ac_kicker', 'ele:cartesian_map', 'ele:chamber_wall', 'ele:control_var', &
          'ele:cylindrical_map', 'ele:elec_multipoles', 'ele:floor', 'ele:gen_attribs', 'ele:gen_gradients', &
          'ele:grid_field', 'ele:head', 'ele:lord_slave', 'ele:mat6', 'ele:methods', &
          'ele:multipoles', 'ele:orbit', 'ele:param', 'ele:photon', 'ele:shape', 'ele:spin_taylor', 'ele:taylor', &
          'ele:twiss', 'ele:wake', 'ele:wall3d', &
          'em_field', 'enum', 'evaluate', 'floor_plan', 'floor_orbit', &
          'global', 'global:opti_de', 'global:optimization', 'global:ran_state', 'help', 'inum', &
          'lat_branch_list', 'lat_calc_done', 'lat_ele_list', 'lat_header', 'lat_list', 'lat_param_units', 'lord_control', &
          'matrix', 'merit', 'orbit_at_s', 'place_buffer', &
          'plot_curve', 'plot_curve_manage', 'plot_graph', 'plot_graph_manage', 'plot_histogram', &
          'plot_lat_layout', 'plot_line', 'plot_list', &
          'plot_symbol', 'plot_template_manage', 'plot_transfer', 'plot1', &
          'ptc_com', 'ring_general', &
          'shape_list', 'shape_manage', 'shape_pattern_list', 'shape_pattern_manage', 'shape_pattern_point_manage', 'shape_set', &
          'show', 'slave_control', 'space_charge_com', 'species_to_int', 'species_to_str', &
          'spin_invariant', 'spin_polarization', 'spin_resonance', 'super_universe', &
          'taylor_map', 'twiss_at_s', 'universe', &
          'var_v1_create', 'var_v1_destroy', 'var_create', 'var_general', 'var_v1_array', 'var_v_array', 'var', &
          'wall3d_radius', 'wave']

! Switch (flag) lists per completion context, used by tao_completion_mod for tab
! completion of "-switch" tokens and, where identical, read directly by the
! command parsers via tao_switches_for(). A context is a command name, or a
! command name plus a resolved subcommand name.
!
! Shared with the parser (single source of truth via tao_switches_for):
!   'show'         -- tao_show_cmd global switches
!   'show <what>'  -- every tao_next_switch case in tao_show_this
!   'pipe'         -- tao_pipe_cmd leading switches
! Completion-only (parser list intentionally differs, so kept separate):
!   'set'          -- parser also accepts the deprecated -lord_no_set
!   'change'/'place' -- parsed by index() matching, not tao_next_switch
!   'show merit'/'show top10' -- parsed by index() matching in tao_show_this

type tao_switch_set_struct
  character(tao_switch_name_len) :: context = ''
  character(tao_switch_name_len), allocatable :: switches(:)
end type

! Populated once by tao_switch_sets_init (below). Consumers should call
! tao_switches_for(context), which returns the pre-split switch array directly.

type (tao_switch_set_struct), allocatable, protected :: tao_switch_sets(:)

private :: sw_row

! Marker in tao_enum_value_names' ix_names(:) for enum values that have no index.

integer, parameter :: no_enum_index$ = -999999

contains

!------------------------------------------------------------------------------
!+
! Subroutine tao_enum_value_names (who, names, ix_names, ele)
!
! Allowed values of an enumerated Tao or Bmad parameter. This is the single
! source for both "pipe enum" and tab completion of "set ... = <value>".
!
! Input:
!   who         -- character(*): Enum name as accepted by "pipe enum": a Tao name
!                    like "track_type" or "symbol^type", anything containing "color",
!                    or a Bmad switch attribute name like "tracking_method".
!   ele         -- ele_struct, optional: For switch attributes, the element being
!                    set. Restricts the values to those valid for that element.
!   switch_attribs -- logical, optional: If False, do not consult Bmad's switch
!                    attribute table for names not in the Tao list. Default True.
!                    (The Bmad lookup reports unknown names, which is unwanted
!                    when the name is known to be a Tao struct component.)
!
! Output:
!   names(:)    -- character(*), allocatable: Value names. Not allocated if who is
!                    not a known enum.
!   ix_names(:) -- integer, allocatable: Index of each value, or no_enum_index$ for
!                    enums whose values have no index.
!-

subroutine tao_enum_value_names (who, names, ix_names, ele, switch_attribs)

type (ele_struct), optional, target :: ele
logical, optional :: switch_attribs
type (ele_struct), target :: dummy_ele

character(*) who
character(*), allocatable :: names(:)
integer, allocatable :: ix_names(:)

character(40), allocatable :: name_list(:)
character(40) tmp(300), nm
integer itmp(300), n, i

!

if (allocated(names)) deallocate (names)
if (allocated(ix_names)) deallocate (ix_names)
n = 0

if (index(who, 'color') /= 0) then
  do i = lbound(qp_color_name, 1), ubound(qp_color_name, 1)
    call add (qp_color_name(i), i)
  enddo
  call done
  return
endif

select case (who)
case ('axis^type')
  call add ('LINEAR', 1); call add ('LOG', 2)
case ('bounds')
  call add ('GENERAL', 1); call add ('ZERO_AT_END', 2); call add ('ZERO_SYMMETRIC', 3)
case ('building^constraint')
  call add ('none', 1); call add ('left_side', 2); call add ('right_side', 3)
case ('data^merit_type')
  call add_array (tao_data_merit_type_name)
case ('data_source')
  call add_array (tao_data_source_name)
case ('distribution_type')
  call add_array (beam_distribution_type_name)
case ('floor_plan_view_name')
  call add_array (tao_floor_plan_view_name)
case ('graph^type')
  call add_array (tao_graph_type_name)
case ('line^pattern', 'orbit_pattern')
  call add_array (qp_line_pattern_name)
case ('lord_status')
  call add ('Group_Lord', 4);     call add ('Super_Lord', 5);      call add ('Overlay_Lord', 6)
  call add ('Girder_Lord', 7);    call add ('Multipass_Lord', 8);  call add ('Not_a_Lord', 10)
  call add ('Control_Lord', 12);  call add ('Ramper_Lord', 13)
case ('optimizer')
  call add_array (tao_optimizer_name)
case ('orbit_lattice')
  call add ('model', 1); call add ('design', 2); call add ('base', 3)
case ('photon_type')
  call add_array (photon_type_name)
case ('plot^type')
  call add ('normal', 1); call add ('wave', 2)
case ('random_engine')
  call add ('pseudo', 1); call add ('quasi', 2)
case ('random_gauss_converter')
  call add ('exact', 1); call add ('quick', 2)
case ('shape^label')
  call add_array (tao_shape_label_name)
case ('shape^shape')
  call add_array (tao_shape_shape_name)
case ('slave_status')
  call add ('Minor_Slave', 1);  call add ('Super_Slave', 2);  call add ('Free', 3)
  call add ('Multipass_Slave', 9);  call add ('Slice_Slave', 11)
case ('fill_pattern')
  call add_array (qp_symbol_fill_pattern_name)
case ('symbol^type')
  call add_array (qp_symbol_type_name)
case ('track_type')
  call add ('single', no_enum_index$); call add ('beam', no_enum_index$)
case ('var^merit_type')
  call add_array (tao_var_merit_type_name)
case ('view')
  call add ('zx', no_enum_index$); call add ('xz', no_enum_index$); call add ('xy', no_enum_index$)
  call add ('yx', no_enum_index$); call add ('zy', no_enum_index$); call add ('yz', no_enum_index$)
case ('wave_data_type')
  call add_array (tao_wave_data_name)
case ('x_axis_type')
  call add_array (tao_x_axis_type_name)
case ('data_type_z')
  call add_array (tao_data_type_z_name)

case default
  ! A Bmad switch attribute.
  if (present(switch_attribs)) then
    if (.not. switch_attribs) return
  endif
  nm = upcase(who)
  if (nm == 'EVAL_POINT') nm = 'ELE_ORIGIN'  ! data%eval_point is not recognized by switch_attrib_value_name
  if (present(ele)) then
    nm = switch_attrib_value_name(nm, 1.0_rp, ele, name_list = name_list)
  else
    nm = switch_attrib_value_name(nm, 1.0_rp, dummy_ele, name_list = name_list)
  endif
  if (.not. allocated(name_list)) return   ! Unknown enum: names left unallocated.
  do i = lbound(name_list, 1), ubound(name_list, 1)
    if (name_list(i) == '' .or. index(name_list(i), '!') /= 0) cycle
    call add (name_list(i), i)
  enddo
end select

call done

!------------------------------------------
contains

subroutine add (name, ix)
character(*) name
integer ix
if (n >= size(tmp)) return
n = n + 1
tmp(n) = name
itmp(n) = ix
end subroutine add

subroutine add_array (arr)
character(*) arr(:)
integer ia
do ia = lbound(arr, 1), ubound(arr, 1)
  call add (arr(ia), ia)
enddo
end subroutine add_array

subroutine done
allocate (names(n), ix_names(n))
names = tmp(1:n)
ix_names = itmp(1:n)
end subroutine done

end subroutine tao_enum_value_names

!------------------------------------------------------------------------------
!+
! Function tao_switches_for (context) result (switches)
!
! Return the switch (flag) name array for a completion/parsing context, or a
! zero-length array if the context has no switches. The context is a command
! name, or "show <subcommand>".
!-

function tao_switches_for (context) result (switches)

character(*), intent(in) :: context
character(tao_switch_name_len), allocatable :: switches(:)
integer i

call tao_switch_sets_init()

do i = 1, size(tao_switch_sets)
  if (tao_switch_sets(i)%context == context) then
    switches = tao_switch_sets(i)%switches
    return
  endif
enddo

allocate (switches(0))

end function tao_switches_for

!------------------------------------------------------------------------------
!+
! Subroutine tao_switch_sets_init ()
!
! One-time build of the tao_switch_sets(:) table. Idempotent.
!-

subroutine tao_switch_sets_init ()

if (allocated(tao_switch_sets)) return

tao_switch_sets = [ &
  sw_row('change', [character(tao_switch_name_len):: '-silent', '-update', '-listing', '-branch', '-mask']), &
  sw_row('pipe',   [character(tao_switch_name_len):: '-append', '-write', '-noprint']), &
  sw_row('place',  [character(tao_switch_name_len):: '-no_buffer']), &
  sw_row('set',    [character(tao_switch_name_len):: '-update', '-mask', '-branch', '-listing', '-silent']), &
  sw_row('show',   [character(tao_switch_name_len):: '-append', '-write', '-noprint', '-no_err_out']), &
  sw_row('show beam', [character(tao_switch_name_len):: '-universe', '-lattice', '-comb', '-z']), &
  sw_row('show branch', [character(tao_switch_name_len):: '-universe']), &
  sw_row('show chromaticity', [character(tao_switch_name_len):: '-universe', '-taylor']), &
  sw_row('show curve', [character(tao_switch_name_len):: '-symbol', '-line', '-no_header']), &
  sw_row('show derivative', [character(tao_switch_name_len):: '-derivative_recalc']), &
  sw_row('show element', [character(tao_switch_name_len):: '-taylor', '-em_field', '-all', '-data', &
      '-design', '-no_slaves', '-wall', '-base', '-field', '-floor_coords', '-xfer_mat', '-ptc', &
      '-everything', '-attributes', '-no_super_slaves', '-radiation_kick', '-internal']), &
  sw_row('show emittance', [character(tao_switch_name_len):: '-universe', '-element', '-xmatrix', '-sigma_matrix']), &
  sw_row('show field', [character(tao_switch_name_len):: '-derivatives', '-grid_pt', '-percent_len', '-absolute_s']), &
  sw_row('show global', [character(tao_switch_name_len):: '-optimization', '-bmad_com', '-environment', &
      '-csr_param', '-space_charge_com', '-ran_state', '-ptc_com', '-internal']), &
  sw_row('show graph', [character(tao_switch_name_len):: '-debug', '-rms']), &
  sw_row('show history', [character(tao_switch_name_len):: '-no_num', '-all', '-filed']), &
  sw_row('show internal', [character(tao_switch_name_len):: '-pipe', '-control']), &
  sw_row('show lattice', [character(tao_switch_name_len):: '-branch', '-blank_replacement', '-lords', &
      '-center', '-middle', '-tracking_elements', '-0undef', '-beginning', '-pipe', '-no_label_lines', &
      '-no_tail_lines', '-custom', '-s', '-radiation_integrals', '-remove_line_if_zero', '-base', &
      '-design', '-floor_coords', '-orbit', '-attribute', '-all', '-no_slaves', '-energy', '-spin', &
      '-undef0', '-no_super_slaves', '-sum_radiation_integrals', '-python', '-universe', '-rms', &
      '-6d_radiation_integrals', '-ri_radiation_integrals']), &
  sw_row('show matrix', [character(tao_switch_name_len):: '-order', '-s', '-ptc', '-eigen_modes', &
      '-elements', '-lattice_format', '-universe', '-angle_coordinates', '-number_format', &
      '-inverse', '-radiation', '-scibmad', '-noclean']), &
  sw_row('show merit', [character(tao_switch_name_len):: '-derivative', '-merit_only']), &
  sw_row('show particle', [character(tao_switch_name_len):: '-element', '-particle', '-bunch', '-lost', '-all']), &
  sw_row('show plot', [character(tao_switch_name_len):: '-floor_plan', '-lat_layout', '-templates', &
      '-global', '-regions', '-plot_page', '-page']), &
  sw_row('show ptc', [character(tao_switch_name_len):: '-emittance']), &
  sw_row('show radiation_integrals', [character(tao_switch_name_len):: '-branch']), &
  sw_row('show rampers', [character(tao_switch_name_len):: '-universe', '-energy_show']), &
  sw_row('show spin', [character(tao_switch_name_len):: '-element', '-n_axis', '-l_axis', '-g_map', &
      '-flip_n_axis', '-x_zero', '-y_zero', '-z_zero', '-ignore_kinetic', '-isf', '-spin_tune']), &
  sw_row('show symbolic_numbers', [character(tao_switch_name_len):: '-physical_constants', '-lattice_constants']), &
  sw_row('show taylor_map', [character(tao_switch_name_len):: '-order', '-s', '-ptc', '-eigen_modes', &
      '-elements', '-lattice_format', '-universe', '-angle_coordinates', '-number_format', &
      '-inverse', '-radiation', '-scibmad', '-noclean']), &
  sw_row('show top10', [character(tao_switch_name_len):: '-derivative', '-merit_only']), &
  sw_row('show track', [character(tao_switch_name_len):: '-e_field', '-b_field', '-velocity', '-momentum', &
      '-energy', '-position', '-no_label_lines', '-s', '-spin', '-points', '-time', '-range', &
      '-twiss', '-dispersion', '-branch', '-universe', '-design', '-base', '-element']), &
  sw_row('show twiss_and_orbit', [character(tao_switch_name_len):: '-branch', '-universe', '-design', '-base']), &
  sw_row('show universe', [character(tao_switch_name_len):: '-branch']), &
  sw_row('show variables', [character(tao_switch_name_len):: '-bmad_format', '-good_opt_only', &
      '-no_label_lines', '-universe']), &
  sw_row('show wall', [character(tao_switch_name_len):: '-section', '-element', '-angle', '-s', '-branch'])]

end subroutine tao_switch_sets_init

!------------------------------------------------------------------------------
! Row constructor helper for tao_switch_sets_init, keeping each row on one line.

function sw_row (context, switches) result (set)
character(*), intent(in) :: context, switches(:)
type (tao_switch_set_struct) :: set
set%context = context
set%switches = switches
end function sw_row

end module
