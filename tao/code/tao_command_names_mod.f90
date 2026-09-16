!+
! Module tao_command_names_mod
!
! Central definitions of the command name lists used for command parsing.
! Hoisted here (rather than being literals at the match_word call sites) so the
! tab completion engine (tao_completion_mod) can offer the same names the
! parsers accept. Keep each list in sync with the corresponding select case.
!-

module tao_command_names_mod

implicit none

! Top level Tao commands. Matched case sensitively in tao_command.

character(16), parameter :: tao_command_names(49) = [character(16):: &
                      'alias', 'call', 'change', 'clear', 'clip', 'continue', 'create', 'cut_ring', 'derivative', &
                      'end_file', 'exit', 'fixer', 'flatten', 'help', 'json', 'ls', 'misalign', 'pause', 'pipe', 'place', &
                      'plot', 'ptc', 'python', 'quit', 're_execute', 'read', 'regression', 'reinitialize', 'reset', &
                      'restore', 'run_optimizer', 'scale', 'set', 'show', 'single_mode', 'spawn', 'taper', &
                      'timer', 'use', 'veto', 'view', 'wave', 'write', 'x_axis', 'x_scale', 'xy_scale', &
                      'debug', 'verbose', 'tree']

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

character(40), parameter :: tao_pipe_cmd_names(113) = [character(40) :: &
          'beam', 'beam_init', 'branch1', 'bunch_comb', 'bunch_params', 'bunch1', 'bmad_com',&
          'building_wall_list', 'building_wall_graph', 'building_wall_point', 'building_wall_section', &
          'complete', 'constraints', 'da_params', 'da_aperture', &
          'data', 'data_d2_create', 'data_d2_destroy', 'data_d_array', 'data_d1_array', &
          'data_d2', 'data_d2_array', 'data_set_design_value', 'data_parameter', &
          'datum_create', 'datum_has_ele', 'derivative', &
          'ele:ac_kicker', 'ele:cartesian_map', 'ele:chamber_wall', 'ele:control_var', &
          'ele:cylindrical_map', 'ele:elec_multipoles', 'ele:floor', 'ele:gen_attribs', 'ele:gen_gradients', &
          'ele:grid_field', 'ele:head', 'ele:lord_slave', 'ele:mat6', 'ele:methods', &
          'ele:multipoles', 'ele:orbit', 'ele:param', 'ele:photon', 'ele:spin_taylor', 'ele:taylor', &
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

end module
