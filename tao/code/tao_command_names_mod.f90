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

! Switch (flag) lists per completion context, used by tao_completion_mod for tab
! completion of "-switch" tokens. A context is a command name, or a command name
! plus a resolved subcommand name. The two arrays are parallel: tao_switch_list(i)
! holds the blank-separated switches for tao_switch_context(i).
! Keep in sync with the switch parsing at:
!   'set', 'change', 'place', 'pipe'  -- tao_command / tao_pipe_cmd
!   'show'                            -- tao_show_cmd (pre-dispatch switches)
!   'show <what>'                     -- the corresponding case in tao_show_this

integer, parameter :: n_tao_switch_contexts = 34

character(28), parameter :: tao_switch_context(n_tao_switch_contexts) = [character(28):: &
    'change', 'pipe', 'place', 'set', 'show', &
    'show beam', 'show branch', 'show chromaticity', 'show curve', 'show derivative', &
    'show element', 'show emittance', 'show field', 'show global', 'show graph', &
    'show history', 'show internal', 'show lattice', 'show matrix', 'show merit', &
    'show particle', 'show plot', 'show ptc', 'show radiation_integrals', 'show rampers', &
    'show spin', 'show symbolic_numbers', 'show taylor_map', 'show top10', 'show track', &
    'show twiss_and_orbit', 'show universe', 'show variables', 'show wall']

character(400), parameter :: tao_switch_list(n_tao_switch_contexts) = [character(400):: &
    '-silent -update -listing -branch -mask', &
    '-append -write -noprint', &
    '-no_buffer', &
    '-update -mask -branch -listing -silent', &
    '-append -write -noprint -no_err_out', &
    '-universe -lattice -comb -z', &
    '-universe', &
    '-universe -taylor', &
    '-symbol -line -no_header', &
    '-derivative_recalc', &
    '-taylor -em_field -all -data -design -no_slaves -wall -base -field -floor_coords &
    &-xfer_mat -ptc -everything -attributes -no_super_slaves -radiation_kick -internal', &
    '-universe -element -xmatrix -sigma_matrix', &
    '-derivatives -grid_pt -percent_len -absolute_s', &
    '-optimization -bmad_com -environment -csr_param -space_charge_com -ran_state -ptc_com -internal', &
    '-debug -rms', &
    '-no_num -all -filed', &
    '-pipe -control', &
    '-branch -blank_replacement -lords -center -middle -tracking_elements -0undef -beginning &
    &-pipe -no_label_lines -no_tail_lines -custom -s -radiation_integrals -remove_line_if_zero &
    &-base -design -floor_coords -orbit -attribute -all -no_slaves -energy -spin -undef0 &
    &-no_super_slaves -sum_radiation_integrals -python -universe -rms -6d_radiation_integrals &
    &-ri_radiation_integrals', &
    '-order -s -ptc -eigen_modes -elements -lattice_format -universe -angle_coordinates &
    &-number_format -inverse -radiation -scibmad -noclean', &
    '-derivative -merit_only', &
    '-element -particle -bunch -lost -all', &
    '-floor_plan -lat_layout -templates -global -regions -plot_page -page', &
    '-emittance', &
    '-branch', &
    '-universe -energy_show', &
    '-element -n_axis -l_axis -g_map -flip_n_axis -x_zero -y_zero -z_zero -ignore_kinetic -isf -spin_tune', &
    '-physical_constants -lattice_constants', &
    '-order -s -ptc -eigen_modes -elements -lattice_format -universe -angle_coordinates &
    &-number_format -inverse -radiation -scibmad -noclean', &
    '-derivative -merit_only', &
    '-e_field -b_field -velocity -momentum -energy -position -no_label_lines -s -spin -points &
    &-time -range -twiss -dispersion -branch -universe -design -base -element', &
    '-branch -universe -design -base', &
    '-branch', &
    '-bmad_format -good_opt_only -no_label_lines -universe', &
    '-section -element -angle -s -branch']

end module
