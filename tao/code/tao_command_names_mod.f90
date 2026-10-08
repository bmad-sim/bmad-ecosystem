!+
! Module tao_command_names_mod
!
! Command name and switch lists shared by the command parsers and tab completion
! (tao_completion_mod). Keep each list in sync with the corresponding select case.
!-

module tao_command_names_mod

use tao_struct
use quick_plot
use attribute_mod, only: switch_attrib_value_name

implicit none

! Top level Tao commands. Matched case sensitively in tao_command. The hidden
! commands are accepted by the parser but not documented, so completion omits them.

character(16), parameter :: tao_visible_command_names(46) = [character(16):: &
                      'alias', 'call', 'change', 'clear', 'clip', 'continue', 'create', 'cut_ring', 'derivative', &
                      'end_file', 'exit', 'fixer', 'flatten', 'help', 'json', 'ls', 'misalign', 'pause', 'pipe', 'place', &
                      'plot', 'ptc', 'python', 'quit', 're_execute', 'read', 'regression', 'reinitialize', 'reset', &
                      'restore', 'run_optimizer', 'scale', 'set', 'show', 'single_mode', 'spawn', 'taper', &
                      'timer', 'use', 'veto', 'view', 'wave', 'write', 'x_axis', 'x_scale', 'xy_scale']

character(16), parameter :: tao_hidden_command_names(3) = [character(16):: 'debug', 'verbose', 'tree']

character(16), parameter :: tao_command_names(49) = [tao_visible_command_names, tao_hidden_command_names]

! "write <action>" names. See tao_write_cmd.

character(20), parameter :: tao_write_action_names(34) = [character(20):: &
              '3d_model', 'beam', 'bmad', 'blender', 'bunch_comb', 'covariance_matrix', 'curve', &
              'derivative_matrix', 'digested', 'elegant', 'field', &
              'gif', 'gif-l', 'hard', 'hard-l', 'mad', 'mad8', 'madx', 'matrix', &
              'namelist', 'opal', 'pals', 'pdf', 'pdf-l', 'plot_commands', 'ps', 'ps-l', 'ptc', &
              'sad', 'scibmad', 'spin_mat8', 'tao', 'variable', 'xsif']

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

! "set ptc_com <component>" names. tao_set_ptc_com_cmd validates against the full list
! before its select case; completion offers T/F for the logical ones.

character(24), parameter :: tao_set_ptc_com_logical_names(6) = [character(24):: 'exact_model', 'exact_misalign', &
          'use_orientation_patches', 'print_info_messages', 'pancake_symplectic', 'pancake_canonical']

character(24), parameter :: tao_set_ptc_com_names(10) = [character(24):: 'vertical_kick', 'cut_factor', &
          'max_fringe_order', 'old_integrator', tao_set_ptc_com_logical_names]

! "set beam <parameter>" names. See tao_set_beam_cmd, which accepts the deprecated
! names too (tao_set_beam_all_names); completion offers only the current ones.

character(24), parameter :: tao_set_beam_names(11) = [character(24):: 'beginning', 'comb_ds_save', 'always_reinit', &
          'track_start', 'track_end', 'beam_init_position_file', 'dump_file', 'dump_at', 'saved_at', &
          'add_saved_at', 'subtract_saved_at']

character(24), parameter :: tao_set_beam_deprecated_names(6) = [character(24):: 'beam_track_start', 'beam_track_end', &
          'beam_init_file_name', 'beam_saved_at', 'beam_dump_at', 'beam_dump_file']

character(24), parameter :: tao_set_beam_all_names(17) = [tao_set_beam_names, tao_set_beam_deprecated_names]

! "change <what>" names. The change case in tao_command matches these by abbreviation
! ("particle_start" may carry an "n@" prefix).

character(16), parameter :: tao_change_what_names(5) = [character(16):: 'element', 'variable', 'tune', 'z_tune', &
                                                                        'particle_start']

! "read <what>" names. See tao_read_cmd.

character(8), parameter :: tao_read_what_names(2) = [character(8):: 'lattice', 'ptc']

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

! Switch lists are in tao_switches_for (below). Switch names are character(28).

integer, parameter :: tao_switch_name_len = 28

! Marker in tao_enum_value_names' ix_names(:) for enum values that have no index.

integer, parameter :: no_enum_index$ = -999999

contains

!------------------------------------------------------------------------------
!+
! Subroutine tao_enum_value_names (who, names, ix_names, ele, switch_attribs)
!
! Allowed values of an enumerated Tao or Bmad parameter, for "pipe enum" and for
! tab completion of "set ... = <value>".
!
! Input:
!   who         -- character(*): Enum name as accepted by "pipe enum": a Tao name
!                    like "track_type" or "symbol^type", "prompt_color" (terminal colors),
!                    anything else containing "color" (plot colors), or a Bmad switch
!                    attribute name like "tracking_method".
!   ele         -- ele_struct, optional: For switch attributes, the element being
!                    set. Restricts the values to those valid for that element.
!   switch_attribs -- logical, optional: If False, do not consult Bmad's switch
!                    attribute table for names not in the Tao list (that lookup
!                    prints an error for unknown names). Default True.
!
! Output:
!   names(:)    -- character(*), allocatable: Value names. Not allocated if who is
!                    not a known enum.
!   ix_names(:) -- integer, allocatable: Index of each value, or no_enum_index$ for
!                    enums whose values have no index.
!-

subroutine tao_enum_value_names (who, names, ix_names, ele, switch_attribs)

type (ele_struct), optional, target :: ele
type (ele_struct), target :: dummy_ele
type (ele_struct), pointer :: ele_ptr
logical, optional :: switch_attribs

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

! The prompt takes terminal colors, not plot colors, so this must precede the generic test.

if (who == 'prompt_color') then
  do i = 1, size(terminal_color_name)
    call add (terminal_color_name(i), no_enum_index$)
  enddo
  call done
  return
endif

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
  ele_ptr => dummy_ele
  if (present(ele)) ele_ptr => ele
  nm = switch_attrib_value_name(nm, 1.0_rp, ele_ptr, name_list = name_list)
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
! Switch names for a command, or for "show <subcommand>". Zero-length if none.
! Shared by the parsers (tao_command, tao_show_cmd, tao_show_this, tao_pipe_cmd) and
! by tab completion. The 'place', 'show merit' and 'show top10' lists are
! completion-only since those parsers accept a different set. The set parser also
! accepts the deprecated '-lord_no_set', which is deliberately left out of 'set'.
!-

function tao_switches_for (context) result (switches)

character(*), intent(in) :: context
character(28), allocatable :: switches(:)  ! TODO match tao_switch_name_len (sorry, this is for cppbmad)
integer, parameter :: c = tao_switch_name_len

select case (context)
case ('change');                    switches = [character(c):: '-silent', '-update', '-listing', '-branch', '-mask']
case ('pipe');                      switches = [character(c):: '-append', '-write', '-noprint']
case ('place');                     switches = [character(c):: '-no_buffer']
case ('read');                      switches = [character(c):: '-universe', '-silent']
case ('set');                       switches = [character(c):: '-update', '-mask', '-branch', '-listing', '-silent']
case ('show');                      switches = [character(c):: '-append', '-write', '-noprint', '-no_err_out']
case ('show beam');                 switches = [character(c):: '-universe', '-lattice', '-comb', '-z']
case ('show branch');               switches = [character(c):: '-universe']
case ('show chromaticity');         switches = [character(c):: '-universe', '-taylor']
case ('show curve');                switches = [character(c):: '-symbol', '-line', '-no_header']
case ('show derivative');           switches = [character(c):: '-derivative_recalc']
case ('show element')
  switches = [character(c):: '-taylor', '-em_field', '-all', '-data', '-design', '-no_slaves', '-wall', '-base', '-field', &
                  '-floor_coords', '-xfer_mat', '-ptc', '-everything', '-attributes', '-no_super_slaves', &
                  '-radiation_kick', '-internal']
case ('show emittance');            switches = [character(c):: '-universe', '-element', '-xmatrix', '-sigma_matrix']
case ('show field');                switches = [character(c):: '-derivatives', '-grid_pt', '-percent_len', '-absolute_s']
case ('show global')
  switches = [character(c):: '-optimization', '-bmad_com', '-environment', '-csr_param', '-space_charge_com', '-ran_state', &
                  '-ptc_com', '-internal']
case ('show graph');                switches = [character(c):: '-debug', '-rms']
case ('show history');              switches = [character(c):: '-no_num', '-all', '-filed']
case ('show internal');             switches = [character(c):: '-pipe', '-control']
case ('show lattice')
  switches = [character(c):: '-branch', '-blank_replacement', '-lords', '-center', '-middle', '-tracking_elements', '-0undef', &
                  '-beginning', '-pipe', '-no_label_lines', '-no_tail_lines', '-custom', '-s', '-radiation_integrals', &
                  '-remove_line_if_zero', '-base', '-design', '-floor_coords', '-orbit', '-attribute', '-all', &
                  '-no_slaves', '-energy', '-spin', '-undef0', '-no_super_slaves', '-sum_radiation_integrals', &
                  '-python', '-universe', '-rms', '-6d_radiation_integrals', '-ri_radiation_integrals']
case ('show matrix', 'show taylor_map')
  switches = [character(c):: '-order', '-s', '-ptc', '-eigen_modes', '-elements', '-lattice_format', '-universe', &
                  '-angle_coordinates', '-number_format', '-inverse', '-radiation', '-scibmad', '-noclean']
case ('show merit', 'show top10');  switches = [character(c):: '-derivative', '-merit_only']
case ('show particle');             switches = [character(c):: '-element', '-particle', '-bunch', '-lost', '-all']
case ('show plot')
  switches = [character(c):: '-floor_plan', '-lat_layout', '-templates', '-global', '-regions', '-plot_page', '-page']
case ('show ptc');                  switches = [character(c):: '-emittance']
case ('show radiation_integrals');  switches = [character(c):: '-branch']
case ('show rampers');              switches = [character(c):: '-universe', '-energy_show']
case ('show spin')
  switches = [character(c):: '-element', '-n_axis', '-l_axis', '-g_map', '-flip_n_axis', '-x_zero', '-y_zero', '-z_zero', &
                  '-ignore_kinetic', '-isf', '-spin_tune']
case ('show symbolic_numbers');     switches = [character(c):: '-physical_constants', '-lattice_constants']
case ('show track')
  switches = [character(c):: '-e_field', '-b_field', '-velocity', '-momentum', '-energy', '-position', '-no_label_lines', '-s', &
                  '-spin', '-points', '-time', '-range', '-twiss', '-dispersion', '-branch', '-universe', '-design', &
                  '-base', '-element']
case ('show twiss_and_orbit');      switches = [character(c):: '-branch', '-universe', '-design', '-base']
case ('show universe');             switches = [character(c):: '-branch']
case ('show variables');            switches = [character(c):: '-bmad_format', '-good_opt_only', '-no_label_lines', '-universe']
case ('show wall');                 switches = [character(c):: '-section', '-element', '-angle', '-s', '-branch']
case default;                       allocate (switches(0))
end select

end function tao_switches_for

end module
