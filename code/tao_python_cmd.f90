!+
! Subroutine tao_python_cmd (input_str)
!
! Print information in a form easily parsed by a scripting program like python.
!
! Note: The syntax for "variable list form" is:
!   <component_name>;<type>;<variable>;<component_value>
! <type> is one of:
!   STR
!   INT
!   REAL
!   LOGIC
! <variable> indicates if the component can be varied. It is one of:
!   T
!   F
!
! Input:
!   input_str  -- Character(*): What to show.
!-


subroutine tao_python_cmd (input_str)

use tao_mod
use tao_command_mod
use location_encode_mod

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
type (tao_plot_region_struct), pointer :: pr
type (tao_plot_struct), pointer :: p
type (tao_graph_struct), pointer :: g
type (tao_curve_struct), pointer :: cur
type (tao_lattice_struct), pointer :: tao_lat
type (tao_plot_region_struct), pointer :: region
type (tao_d2_data_array_struct), allocatable, save :: d2_array(:)
type (tao_d1_data_array_struct), allocatable, save :: d1_array(:)
type (tao_data_array_struct), allocatable, save :: d_array(:)
type (beam_struct), pointer :: beam
type (beam_init_struct), pointer :: beam_init
type (lat_struct), pointer :: lat
type (bunch_struct), pointer :: bunch
type (wake_lr_mode_struct), pointer :: lr_mode
type (ele_struct), pointer :: ele
type (coord_struct), target :: orb
type (ele_struct), target :: this_ele
type (bunch_params_struct), pointer :: bunch_params
type (bunch_params_struct), pointer :: bunch_p
type (ele_pointer_struct), allocatable, save :: eles(:)
type (branch_struct), pointer :: branch
type (tao_universe_branch_struct), pointer :: uni_branch
type (random_state_struct) ran_state
type (tao_scratch_space_struct), pointer :: ss
type (ele_attribute_struct) attrib

character(*) input_str
character(n_char_show), pointer :: li(:) 
character(24) imt, rmt, lmt, amt, iamt, vamt, vrmt
character(40) max_loc, loc_ele, name1(40), name2(40), a_name, name
character(200) line, file_name
character(20) cmd, command, who
character(20) :: r_name = 'tao_python_cmd'
character(20) :: cmd_names(24)= &
          ['help           ', 'global         ', 'plot_list      ', 'plot1          ', 'plot_graph     ', &
           'plot_curve     ', 'plot_line      ', 'plot_symbol    ', 'universe       ', 'var_general    ', &
           'var_v1         ', 'var1           ', 'data_general   ', 'data_d2        ', 'data_d1        ', &
           'data1          ', 'beam_init      ', 'bunch1         ', 'lat_general    ', 'branch1        ', &
           'lat_ele_list   ', 'lat_ele1       ', 'twiss_at_s     ', 'orbit_at_s     ']


real(rp) s_pos

integer :: i, j, ie, iu, md, nl, ct, n1, nl2, n, ix, ix2, iu_write
integer :: ix_ele, ix_ele1, ix_ele2, ix_branch, ix_universe
integer :: ios, n_loc
logical :: err, print_flag, opened, doprint, free

character(20) switch

!

line = input_str
doprint = .true.
opened = .false.

do
  call tao_next_switch (line, ['-append ', '-write  ', '-noprint'], .false., switch, err, ix)
  if (err) return
  if (switch == '') exit

  select case (switch)
  case ('-noprint')
    doprint = .false.

  case ('-append', '-write')
    call string_trim(line, line, ix)
    file_name = line(:ix)
    call string_trim(line(ix+1:), line, ix)

    iu_write = lunget()
    if (switch == '-append') then
      open (iu_write, file = file_name, position = 'APPEND', status = 'UNKNOWN', recl = 200)
    else
      open (iu_write, file = file_name, status = 'REPLACE', recl = 200)
    endif

    opened = .true.
  end select
enddo

call string_trim(line, line, ix)
cmd = line(1:ix)
call string_trim(line(ix+1:), line, ix)

call match_word (cmd, cmd_names, ix, matched_name = command)
if (ix == 0) then
  call out_io (s_error$, r_name, 'python what? "What" not recognized: ' // command)
  return
endif

if (ix < 0) then
  call out_io (s_error$, r_name, 'python what? Ambiguous command: ' // command)
  return
endif

amt = '(3a)'
imt = '(a,i0,a)'
rmt = '(a,es21.13,a)'
lmt = '(a,l1,a)'
vamt = '(a, i0, 3a)'
vrmt = '(a, i0, a, es21.13)'

nl = 0
ss => scratch
li => ss%lines

call re_allocate_lines (200)

select case (command)

!----------------------------------------------------------------------
! help
! returns list of "help xxx" topics

case ('help')

  call tao_help ('help-list', '', ss%lines, n)

  nl2 = 0
  do i = 1, n
    if (li(i) == '') cycle
    call string_trim(li(i), line, ix)
    nl=nl+1; name1(nl) = line(1:ix)
    call string_trim(line(ix+1:), line, ix)
    if (ix == 0) cycle
    nl2=nl2+1; name2(nl2) = line
  enddo

  li(1:nl) = name1(1:nl)
  li(nl+1:nl+nl2) = name2(1:nl2)
  nl = nl + nl2

!----------------------------------------------------------------------
! Global parameters
! Input syntax: 
!   python global
! Output syntax is variable list form. See documentation at beginning of this file.

case ('global')

  nl=nl+1; write (li(nl), rmt) 'y_axis_plot_dmin;REAL;T;',                s%global%y_axis_plot_dmin
  nl=nl+1; write (li(nl), rmt) 'lm_opt_deriv_reinit;REAL;T;',             s%global%lm_opt_deriv_reinit
  nl=nl+1; write (li(nl), rmt) 'de_lm_step_ratio;REAL;T;',                s%global%de_lm_step_ratio
  nl=nl+1; write (li(nl), rmt) 'de_var_to_population_factor;REAL;T;',     s%global%de_var_to_population_factor
  nl=nl+1; write (li(nl), rmt) 'lmdif_eps;REAL;T;',                       s%global%lmdif_eps
  nl=nl+1; write (li(nl), rmt) 'svd_cutoff;REAL;T;',                      s%global%svd_cutoff
  nl=nl+1; write (li(nl), rmt) 'unstable_penalty;REAL;T;',                s%global%unstable_penalty
  nl=nl+1; write (li(nl), rmt) 'merit_stop_value;REAL;T;',                s%global%merit_stop_value
  nl=nl+1; write (li(nl), rmt) 'random_sigma_cutoff;REAL;T;',             s%global%random_sigma_cutoff
  nl=nl+1; write (li(nl), rmt) 'delta_e_chrom;REAL;T;',                   s%global%delta_e_chrom
  nl=nl+1; write (li(nl), imt) 'n_opti_cycles;INT;T;',                    s%global%n_opti_cycles
  nl=nl+1; write (li(nl), imt) 'n_opti_loops;INT;T;',                     s%global%n_opti_loops
  nl=nl+1; write (li(nl), imt) 'phase_units;INT;T;',                      s%global%phase_units
  nl=nl+1; write (li(nl), imt) 'bunch_to_plot;INT;T;',                    s%global%bunch_to_plot
  nl=nl+1; write (li(nl), imt) 'random_seed;INT;T;',                      s%global%random_seed
  nl=nl+1; write (li(nl), imt) 'n_top10;INT;T;',                          s%global%n_top10
  nl=nl+1; write (li(nl), amt) 'random_engine;STR;T;',                    s%global%random_engine
  nl=nl+1; write (li(nl), amt) 'random_gauss_converter;STR;T;',           s%global%random_gauss_converter
  nl=nl+1; write (li(nl), amt) 'track_type;STR;T;',                       s%global%track_type
  nl=nl+1; write (li(nl), amt) 'prompt_string;STR;T;',                    s%global%prompt_string
  nl=nl+1; write (li(nl), amt) 'prompt_color;STR;T;',                     s%global%prompt_color
  nl=nl+1; write (li(nl), amt) 'optimizer;STR;T;',                        s%global%optimizer
  nl=nl+1; write (li(nl), amt) 'print_command;STR;T;',                    s%global%print_command
  nl=nl+1; write (li(nl), amt) 'var_out_file;STR;T;',                     s%global%var_out_file
  nl=nl+1; write (li(nl), lmt) 'initialized;LOGIC;T;',                    s%global%initialized
  nl=nl+1; write (li(nl), lmt) 'opt_with_ref;LOGIC;T;',                   s%global%opt_with_ref
  nl=nl+1; write (li(nl), lmt) 'opt_with_base;LOGIC;T;',                  s%global%opt_with_base
  nl=nl+1; write (li(nl), lmt) 'label_lattice_elements;LOGIC;T;',         s%global%label_lattice_elements
  nl=nl+1; write (li(nl), lmt) 'label_keys;LOGIC;T;',                     s%global%label_keys
  nl=nl+1; write (li(nl), lmt) 'derivative_recalc;LOGIC;T;',              s%global%derivative_recalc
  nl=nl+1; write (li(nl), lmt) 'derivative_uses_design;LOGIC;T;',         s%global%derivative_uses_design
  nl=nl+1; write (li(nl), lmt) 'init_plot_needed;LOGIC;T;',               s%global%init_plot_needed
  nl=nl+1; write (li(nl), lmt) 'orm_analysis;LOGIC;T;',                   s%global%orm_analysis
  nl=nl+1; write (li(nl), lmt) 'plot_on;LOGIC;T;',                        s%global%plot_on
  nl=nl+1; write (li(nl), lmt) 'lattice_calc_on;LOGIC;T;',                s%global%lattice_calc_on
  nl=nl+1; write (li(nl), lmt) 'svd_retreat_on_merit_increase;LOGIC;T;',  s%global%svd_retreat_on_merit_increase
  nl=nl+1; write (li(nl), lmt) 'stop_on_error;LOGIC;T;',                  s%global%stop_on_error
  nl=nl+1; write (li(nl), lmt) 'command_file_print_on;LOGIC;T;',          s%global%command_file_print_on
  nl=nl+1; write (li(nl), lmt) 'box_plots;LOGIC;T;',                      s%global%box_plots
  nl=nl+1; write (li(nl), lmt) 'beam_timer_on;LOGIC;T;',                  s%global%beam_timer_on
  nl=nl+1; write (li(nl), lmt) 'var_limits_on;LOGIC;T;',                  s%global%var_limits_on
  nl=nl+1; write (li(nl), lmt) 'only_limit_opt_vars;LOGIC;T;',            s%global%only_limit_opt_vars
  nl=nl+1; write (li(nl), lmt) 'optimizer_var_limit_warn;LOGIC;T;',       s%global%optimizer_var_limit_warn
  nl=nl+1; write (li(nl), lmt) 'rf_on;LOGIC;T;',                          s%global%rf_on
  nl=nl+1; write (li(nl), lmt) 'draw_curve_off_scale_warn;LOGIC;T;',      s%global%draw_curve_off_scale_warn
  nl=nl+1; write (li(nl), lmt) 'wait_for_cr_in_single_mode;LOGIC;T;',     s%global%wait_for_CR_in_single_mode
  nl=nl+1; write (li(nl), lmt) 'disable_smooth_line_calc;LOGIC;T;',       s%global%disable_smooth_line_calc
  nl=nl+1; write (li(nl), lmt) 'debug_on;LOGIC;T;',                       s%global%debug_on


!----------------------------------------------------------------------
! List of plot templates or plot regions.
! Input syntax:  
!   python plot_list <r/g>
! where "<r/g>" is:
!   "r"      ! list regions
!   "t"      ! list template plots 


case ('plot_list')
  if (line == 't') then
    do i = 1, size(s%plot_page%template)
      p => s%plot_page%template(i)
      if (p%phantom) cycle
      if (p%name == '') cycle
      if (p%name == 'scratch') cycle
      nl=nl+1; write (li(nl), '(i0, 2a)') i, ';', trim(p%name)
    enddo

  elseif (line == 'r') then
    do i = 1, size(s%plot_page%region)
      pr => s%plot_page%region(i)
      if (pr%name == '') cycle
      nl=nl+1; write (li(nl), '(i0, 5a, l1)') i, ';', trim(pr%name), ';', trim(pr%plot%name), ';', pr%visible
    enddo

  else
    nl=nl+1; li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Expect "r" or "t"')
  endif

!----------------------------------------------------------------------
! Info on a given plot.
! Input syntax:
!   python plot1 <name>
! <name> should be the region name if the plot is associated with a region.
! Output syntax is variable list form. See documentation at beginning of this file.

case ('plot1')

  call tao_find_plots (err, line, 'COMPLETE', plot, print_flag = .false.)
  if (err) then
    nl=nl+1; li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Expect "r" or "t" at end.')
    call end_stuff()
    return
  endif

  p => plot(1)%p

  n = 0
  if (allocated(p%graph)) n = size(p%graph)

  nl=nl+1; write (li(nl), imt) 'num_graphs;INT;T;',                       n
  do i = 1, n
    nl=nl+1; write (li(nl), vamt) 'graph[', i, '];STR;T;',              p%graph(i)%name
  enddo

  nl=nl+1; write (li(nl), amt) 'name;STR;T;',                             p%name
  nl=nl+1; write (li(nl), amt) 'description;STR;T;',                      p%description
  nl=nl+1; write (li(nl), amt) 'x_axis_type;STR;T;',                      p%x_axis_type
  nl=nl+1; write (li(nl), lmt) 'autoscale_x;LOGIC;T;',                    p%autoscale_x
  nl=nl+1; write (li(nl), lmt) 'autoscale_y;LOGIC;T;',                    p%autoscale_y
  nl=nl+1; write (li(nl), lmt) 'autoscale_gang_x;LOGIC;T;',               p%autoscale_gang_x
  nl=nl+1; write (li(nl), lmt) 'autoscale_gang_y;LOGIC;T;',               p%autoscale_gang_y
  nl=nl+1; write (li(nl), lmt) 'list_with_show_plot_command;LOGIC;T;',    p%list_with_show_plot_command
  nl=nl+1; write (li(nl), lmt) 'phantom;LOGIC;T;',                        p%phantom


!----------------------------------------------------------------------
! Graph
! Syntax:
!   python plot_graph <graph_name>
! <graph_name> is in the form:
!   <p_name>.<g_name>
! where 
!   <p_name> is the plot region name if from a region or the plot name if a template plot.
!   This name is obtained from the python plot_list command. 
!   <g_name> is the graph name obtained from the python plot1 command.

case ('plot_graph')

  call tao_find_plots (err, line, 'COMPLETE', graph = graph)

  if (err .or. .not. allocated(graph)) then
    nl=nl+1; li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Bad graph name')
    call end_stuff()
    return
  endif

  g => graph(1)%g

  n = 0
  if (allocated(g%curve)) n = size(g%curve)

  nl=nl+1; write (li(nl), imt) 'num_curves;INT;T;',                       n
  do i = 1, n
    nl=nl+1; write (li(nl), vamt) 'curve[', i, '];STR;T;',                g%curve(i)%name
  enddo

  nl=nl+1; write (li(nl), amt) 'name;STR;T;',                             g%name
  nl=nl+1; write (li(nl), amt) 'type;STR;T;',                             g%type
  nl=nl+1; write (li(nl), amt) 'title;STR;T;',                            g%title
  nl=nl+1; write (li(nl), amt) 'title_suffix;STR;F;',                     g%title_suffix
  nl=nl+1; write (li(nl), amt) 'component;STR;T;',                        g%component
  nl=nl+1; write (li(nl), amt) 'why_invalid;STR;F;',                      g%why_invalid
  nl=nl+1; write (li(nl), amt) 'floor_plan_view;STR;T;',                  g%floor_plan_view
  nl=nl+1; write (li(nl), amt) 'floor_plan_orbit_color;STR;T;',           g%floor_plan_orbit_color
  nl=nl+1; write (li(nl), rmt) 'x_axis_scale_factor;REAL;T;',             g%x_axis_scale_factor
  nl=nl+1; write (li(nl), rmt) 'symbol_size_scale;REAL;T;',               g%symbol_size_scale
  nl=nl+1; write (li(nl), rmt) 'floor_plan_rotation;REAL;T;',             g%floor_plan_rotation
  nl=nl+1; write (li(nl), rmt) 'floor_plan_orbit_scale;REAL;T;',          g%floor_plan_orbit_scale
  nl=nl+1; write (li(nl), imt) 'ix_branch;INT;T;',                        g%ix_branch
  nl=nl+1; write (li(nl), imt) 'ix_universe;INT;T;',                      g%ix_universe
  nl=nl+1; write (li(nl), lmt) 'clip;LOGIC;T;',                           g%clip
  nl=nl+1; write (li(nl), lmt) 'valid;LOGIC;F;',                          g%valid
  nl=nl+1; write (li(nl), lmt) 'y2_mirrors_y;LOGIC;T;',                   g%y2_mirrors_y
  nl=nl+1; write (li(nl), lmt) 'limited;LOGIC;F;',                        g%limited
  nl=nl+1; write (li(nl), lmt) 'draw_axes;LOGIC;T;',                      g%draw_axes
  nl=nl+1; write (li(nl), lmt) 'correct_xy_distortion;LOGIC;T;',          g%correct_xy_distortion
  nl=nl+1; write (li(nl), lmt) 'floor_plan_size_is_absolute;LOGIC;T;',    g%floor_plan_size_is_absolute
  nl=nl+1; write (li(nl), lmt) 'floor_plan_draw_only_first_pass;LOGIC;T;',  g%floor_plan_draw_only_first_pass
  nl=nl+1; write (li(nl), lmt) 'draw_curve_legend;LOGIC;T;',              g%draw_curve_legend
  nl=nl+1; write (li(nl), lmt) 'draw_grid;LOGIC;T;',                      g%draw_grid
  nl=nl+1; write (li(nl), lmt) 'visible;LOGIC;T;',                        g%visible
  nl=nl+1; write (li(nl), lmt) 'draw_only_good_user_data_or_vars;LOGIC;T;', g%draw_only_good_user_data_or_vars

  nl=nl+1; write (li(nl), amt) 'y.label;STR;T;',                         g%y%label
  nl=nl+1; write (li(nl), rmt) 'y.max;REAL;T;',                          g%y%max
  nl=nl+1; write (li(nl), rmt) 'y.min;REAL;T;',                          g%y%min
  nl=nl+1; write (li(nl), imt) 'y.major_div;INT;T;',                     g%y%major_div
  nl=nl+1; write (li(nl), imt) 'y.major_div_nominal;INT;T;',             g%y%major_div_nominal
  nl=nl+1; write (li(nl), imt) 'y.places;INT;T;',                        g%y%places
  nl=nl+1; write (li(nl), lmt) 'y.draw_label;LOGIC;T;',                  g%y%draw_label
  nl=nl+1; write (li(nl), lmt) 'y.draw_numbers;LOGIC;T;',                g%y%draw_numbers

  nl=nl+1; write (li(nl), amt) 'y2.label;STR;T;',                        g%y2%label
  nl=nl+1; write (li(nl), rmt) 'y2.max;REAL;T;',                         g%y2%max
  nl=nl+1; write (li(nl), rmt) 'y2.min;REAL;T;',                         g%y2%min
  nl=nl+1; write (li(nl), imt) 'y2.major_div;INT;T;',                    g%y2%major_div
  nl=nl+1; write (li(nl), imt) 'y2.major_div_nominal;INT;T;',            g%y2%major_div_nominal
  nl=nl+1; write (li(nl), imt) 'y2.places;INT;T;',                       g%y2%places
  nl=nl+1; write (li(nl), lmt) 'y2.draw_label;LOGIC;T;',                 g%y2%draw_label
  nl=nl+1; write (li(nl), lmt) 'y2.draw_numbers;LOGIC;T;',               g%y2%draw_numbers

!----------------------------------------------------------------------
! Curve information for a plot
! Input syntax:
!   pyton curve <curve_name>

case ('plot_curve')

  call tao_find_plots (err, line, 'COMPLETE', curve = curve)

  if (err .or. .not. allocated(curve)) then
    nl=nl+1; li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid curve')
    call end_stuff()
    return
  endif

  cur => curve(1)%c

  nl=nl+1; write (li(nl), amt) 'name;STR;T;',                             cur%name
  nl=nl+1; write (li(nl), amt) 'data_source;STR;T;',                      cur%data_source
  nl=nl+1; write (li(nl), amt) 'data_type_x;STR;T;',                      cur%data_type_x
  nl=nl+1; write (li(nl), amt) 'data_type_z;STR;T;',                      cur%data_type_z
  nl=nl+1; write (li(nl), amt) 'data_type;STR;T;',                        cur%data_type
  nl=nl+1; write (li(nl), amt) 'ele_ref_name;STR;T;',                     cur%ele_ref_name
  nl=nl+1; write (li(nl), amt) 'legend_text;STR;T;',                      cur%legend_text
  nl=nl+1; write (li(nl), amt) 'message_text;STR;T;',                     cur%message_text
  nl=nl+1; write (li(nl), rmt) 'y_axis_scale_factor;REAL;T;',             cur%y_axis_scale_factor
  nl=nl+1; write (li(nl), rmt) 's;REAL;F;',                               cur%s
  nl=nl+1; write (li(nl), rmt) 'z_color0;REAL;T;',                        cur%z_color0
  nl=nl+1; write (li(nl), rmt) 'z_color1;REAL;T;',                        cur%z_color0
  nl=nl+1; write (li(nl), imt) 'ix_universe;INT;T;',                      cur%ix_universe
  nl=nl+1; write (li(nl), imt) 'symbol_every;INT;T;',                     cur%symbol_every
  nl=nl+1; write (li(nl), imt) 'ix_branch;INT;T;',                        cur%ix_branch
  nl=nl+1; write (li(nl), imt) 'ix_ele_ref;INT;T;',                       cur%ix_ele_ref
  nl=nl+1; write (li(nl), imt) 'ix_ele_ref_track;INT;T;',                 cur%ix_ele_ref_track
  nl=nl+1; write (li(nl), imt) 'ix_bunch;INT;T;',                         cur%ix_bunch
  nl=nl+1; write (li(nl), lmt) 'use_y2;LOGIC;T;',                         cur%use_y2
  nl=nl+1; write (li(nl), lmt) 'draw_line;LOGIC;T;',                      cur%draw_line
  nl=nl+1; write (li(nl), lmt) 'draw_symbols;LOGIC;T;',                   cur%draw_symbols
  nl=nl+1; write (li(nl), lmt) 'draw_symbol_index;LOGIC;T;',              cur%draw_symbol_index
  nl=nl+1; write (li(nl), lmt) 'smooth_line_calc;LOGIC;T;',               cur%smooth_line_calc
  nl=nl+1; write (li(nl), lmt) 'use_z_color;LOGIC;T;',                    cur%use_z_color

  nl=nl+1; write (li(nl), imt)  'line.width;INT;T;',                      cur%line%width
  nl=nl+1; write (li(nl), amt)  'line.color;STR;T;',                      qp_color_name(cur%line%color)
  nl=nl+1; write (li(nl), amt)  'line.pattern;STR;T;',                    qp_line_pattern_name(cur%line%pattern)

  nl=nl+1; write (li(nl), amt)  'symbol.type;STR;T;',                     qp_symbol_type_name(cur%symbol%type)
  nl=nl+1; write (li(nl), rmt)  'symbol.height;REAL;T;',                  cur%symbol%height
  nl=nl+1; write (li(nl), amt)  'symbol.fill_pattern;STR;T;',             qp_fill_name(cur%symbol%fill_pattern)
  nl=nl+1; write (li(nl), imt)  'symbol.line_width;INT;T;',               cur%symbol%line_width

!----------------------------------------------------------------------
! Points used to construct a smooth line for a plot curve.
! Input syntax:
!   python plot_line <region_name>.<graph_name>.<curve_name>
! Note: The plot must come from a region, and not a template, since on plots associated with a resion of line data.

case ('plot_line')

  call tao_find_plots (err, line, 'COMPLETE', curve = curve)

  if (.not. allocated(curve) .or. size(curve) /= 1) then
    nl=nl+1; li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid curve')
    call end_stuff()
    return
  endif

  cur => curve(1)%c

  if (.not. allocated(cur%x_line)) then
    nl=nl+1; li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": No line associated with curve')
    call end_stuff()
    return
  endif
    
  call re_allocate_lines (nl+size(cur%x_line)+100)
  do i = 1, size(cur%x_line)
    nl=nl+1; write (li(nl), '(i0, 2(a, es21.13))') i, ';', cur%x_line(i), ';', cur%y_line(i)
  enddo

!----------------------------------------------------------------------
! Locations to draw symbols for a plot curve.
! Input syntax:
!   python plot_symbol <region_name>.<graph_name>.<curve_name>
! Note: The plot must come from a region, and not a template, since on plots associated with a resion of line data.

case ('plot_symbol')

  call tao_find_plots (err, line, 'COMPLETE', curve = curve)

  if (.not. allocated(curve) .or. size(curve) /= 1) then
    nl=nl+1; li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid curve')
    call end_stuff()
    return
  endif

  cur => curve(1)%c

  if (.not. allocated(cur%x_symb)) then
    nl=nl+1; li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": No line associated with curve')
    call end_stuff()
    return
  endif
    
  call re_allocate_lines (size(cur%x_symb)+100)
  do i = 1, size(cur%x_symb)
    nl=nl+1; write (li(nl), '(2(i0, a), 2(es21.13, a))') i, ';', cur%ix_symb(i), ';', cur%x_symb(i), ';', cur%y_symb(i)
  enddo

!----------------------------------------------------------------------
! Universe info
! Input syntax:
!   python universe <ix_universe>
! Use "python global" to get the number of universes.

case ('universe')

  u => point_to_uni(.false., err); if (err) return
  
  nl=nl+1; write (li(nl), imt) 'ix_uni;INT;F;',                           u%ix_uni
  nl=nl+1; write (li(nl), imt) 'n_d2_data_used;INT;F;',                   u%n_d2_data_used
  nl=nl+1; write (li(nl), imt) 'n_data_used;INT;F;',                      u%n_data_used
  nl=nl+1; write (li(nl), lmt) 'reverse_tracking;LOGIC;T;',               u%reverse_tracking
  nl=nl+1; write (li(nl), lmt) 'is_on;LOGIC;T;',                          u%is_on

!----------------------------------------------------------------------
! List of all variable v1 arrays
! Input syntax: 
!   python var_general
! Output syntax:
!   <v1_var name>;<v1_var%v lower bound>;<v1_var%v upper bound>

case ('var_general')

  do i = 1, s%n_v1_var_used
    v1_ptr => s%v1_var(i)
    if (v1_ptr%name == '') cycle
    nl=nl+1; write (li(nl), '(2a, 2(i0, a))') trim(v1_ptr%name), ';', lbound(v1_ptr%v, 1), ';', ubound(v1_ptr%v, 1)
  enddo

!----------------------------------------------------------------------
! List of variables in a given variable v1 array
! Input syntax: 
!   python var_v1 <v1_var>

case ('var_v1')

  call tao_find_var (err, line, v1_array = v1_array)

  if (err .or. .not. allocated(v1_array)) then
    nl=nl+1; li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid v1 name')
    call end_stuff()
    return
  endif

  v1_ptr => v1_array(1)%v1

  do i = lbound(v1_ptr%v, 1), ubound(v1_ptr%v, 1)
    v_ptr => v1_ptr%v(i)
    if (.not. v_ptr%exists) cycle
    nl=nl+1; write (li(nl), '(2a, i0, 5a, 3(es21.13, a), 2 (l1, a))') trim(v1_ptr%name), '[', &
                     v_ptr%ix_v1, '];', trim(v_ptr%ele_name), ';', trim(v_ptr%attrib_name), ';', &
                     v_ptr%meas_value, ';', v_ptr%model_value, ';', &
                     v_ptr%design_value, ';', v_ptr%good_user, ';', v_ptr%useit_opt
  enddo

!----------------------------------------------------------------------
! Info on an individual variable
! Input syntax: 
!   python var1 <var>
! Output syntax is variable list form. See documentation at beginning of this file.

case ('var1')

  call tao_find_var (err, line, v_array = v_array)

  if (.not. allocated(v_array) .or. size(v_array) /= 1) then
    nl=nl+1; li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid variable name')
    call end_stuff()
    return
  endif

  v_ptr => v_array(1)%v

  nl=nl+1; write (li(nl), rmt)  'model_value;REAL;T;',          v_ptr%model_value
  nl=nl+1; write (li(nl), rmt)  'base_value;REAL;T;',           v_ptr%base_value

  nl=nl+1; write (li(nl), amt) 'ele_name;STR;T;',                         v_ptr%ele_name
  nl=nl+1; write (li(nl), amt) 'attrib_name;STR;T;',                      v_ptr%attrib_name
  nl=nl+1; write (li(nl), imt) 'ix_v1;INT;F;',                            v_ptr%ix_v1
  nl=nl+1; write (li(nl), imt) 'ix_var;INT;F;',                           v_ptr%ix_var
  nl=nl+1; write (li(nl), imt) 'ix_dvar;INT;F;',                          v_ptr%ix_dvar
  nl=nl+1; write (li(nl), imt) 'ix_attrib;INT;F;',                        v_ptr%ix_attrib
  nl=nl+1; write (li(nl), imt) 'ix_key_table;INT;T;',                     v_ptr%ix_key_table
  nl=nl+1; write (li(nl), rmt) 'design_value;REAL;F;',                    v_ptr%design_value
  nl=nl+1; write (li(nl), rmt) 'scratch_value;REAL;F;',                   v_ptr%scratch_value
  nl=nl+1; write (li(nl), rmt) 'old_value;REAL;F;',                       v_ptr%old_value
  nl=nl+1; write (li(nl), rmt) 'meas_value;REAL;T;',                      v_ptr%meas_value
  nl=nl+1; write (li(nl), rmt) 'ref_value;REAL;T;',                       v_ptr%ref_value
  nl=nl+1; write (li(nl), rmt) 'correction_value;REAL;F;',                v_ptr%correction_value
  nl=nl+1; write (li(nl), rmt) 'high_lim;REAL;T;',                        v_ptr%high_lim
  nl=nl+1; write (li(nl), rmt) 'low_lim;REAL;T;',                         v_ptr%low_lim
  nl=nl+1; write (li(nl), rmt) 'step;REAL;T;',                            v_ptr%step
  nl=nl+1; write (li(nl), rmt) 'weight;REAL;T;',                          v_ptr%weight
  nl=nl+1; write (li(nl), rmt) 'delta_merit;REAL;F;',                     v_ptr%delta_merit
  nl=nl+1; write (li(nl), rmt) 'merit;REAL;F;',                           v_ptr%merit
  nl=nl+1; write (li(nl), rmt) 'dmerit_dvar;REAL;F;',                     v_ptr%dMerit_dVar
  nl=nl+1; write (li(nl), rmt) 'key_val0;REAL;F;',                        v_ptr%key_val0
  nl=nl+1; write (li(nl), rmt) 'key_delta;REAL;T;',                       v_ptr%key_delta
  nl=nl+1; write (li(nl), rmt) 's;REAL;F;',                               v_ptr%s
  nl=nl+1; write (li(nl), amt) 'merit_type;STR;T;',                       v_ptr%merit_type
  nl=nl+1; write (li(nl), lmt) 'exists;LOGIC;F;',                         v_ptr%exists
  nl=nl+1; write (li(nl), lmt) 'good_var;LOGIC;F;',                       v_ptr%good_var
  nl=nl+1; write (li(nl), lmt) 'good_user;LOGIC;T;',                      v_ptr%good_user
  nl=nl+1; write (li(nl), lmt) 'good_opt;LOGIC;T;',                       v_ptr%good_opt
  nl=nl+1; write (li(nl), lmt) 'good_plot;LOGIC;T;',                      v_ptr%good_plot
  nl=nl+1; write (li(nl), lmt) 'useit_opt;LOGIC;F;',                      v_ptr%useit_opt
  nl=nl+1; write (li(nl), lmt) 'useit_plot;LOGIC;F;',                     v_ptr%useit_plot
  nl=nl+1; write (li(nl), lmt) 'key_bound;LOGIC;T;',                      v_ptr%key_bound

!----------------------------------------------------------------------
! Data d2 info for a given universe.
! Input syntax:
!   python data_general <ix_universe>

case ('data_general')

  u => point_to_uni(.false., err); if (err) return

  do i = 1, u%n_d2_data_used
    d2_ptr => u%d2_data(i)
    if (d2_ptr%name == '') cycle
    nl=nl+1; write (li(nl), '(a)') d2_ptr%name
  enddo

!----------------------------------------------------------------------
! List of d1 arrays in a given data d2.
! Input syntax:
!   python data_d2 <d2_datum>
! <d2_datum> should be of the form 
!   <ix_uni>@<d2_datum_name>

case ('data_d2')


  call tao_find_data (err, line, d2_array = d2_array)

  if (err .or. .not. allocated(d2_array)) then
    nl=nl+1; li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid d2 data name')
    call end_stuff()
    return
  endif

  d2_ptr => d2_array(1)%d2


  do i = lbound(d2_ptr%d1, 1), ubound(d2_ptr%d1, 1)
    nl=nl+1; write (li(nl), '(a, i0, 2a)') 'd1[', i, '];STR;T;', d2_ptr%d1(i)%name
  enddo

  nl=nl+1; write (li(nl), amt) 'name;STR;T;',                             d2_ptr%name
  nl=nl+1; write (li(nl), amt) 'data_file_name;STR;F;',                   d2_ptr%data_file_name
  nl=nl+1; write (li(nl), amt) 'ref_file_name;STR;F;',                    d2_ptr%ref_file_name
  nl=nl+1; write (li(nl), amt) 'data_date;STR;T;',                        d2_ptr%data_date
  nl=nl+1; write (li(nl), amt) 'ref_date;STR;T;',                         d2_ptr%ref_date
  nl=nl+1; write (li(nl), imt) 'ix_uni;INT;F;',                           d2_ptr%ix_uni
  nl=nl+1; write (li(nl), imt) 'ix_data;INT;F;',                          d2_ptr%ix_data
  nl=nl+1; write (li(nl), imt) 'ix_ref;INT;F;',                           d2_ptr%ix_ref
  nl=nl+1; write (li(nl), lmt) 'data_read_in;LOGIC;F;',                   d2_ptr%data_read_in
  nl=nl+1; write (li(nl), lmt) 'ref_read_in;LOGIC;F;',                    d2_ptr%ref_read_in

!----------------------------------------------------------------------
! List of datums in a given data d1 array.
! Input syntax:
!   python data_d1 <ix_universe>@<d2_name>.<d1_datum>
! Use the "python data_d2 <name>" command to get a list of d1 arrays. 
! Use the "python data1" command to get detailed information on a particular datum.
! Example:
!   python data_d1 1@orbit.x

case ('data_d1')

  call tao_find_data (err, line, d1_array = d1_array)

  if (err .or. .not. allocated(d1_array)) then
    nl=nl+1; li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid d1 data name.')
    call end_stuff()
    return
  endif

  d1_ptr => d1_array(1)%d1
  nl=nl+1; write (li(nl), '(2a, 2(i0, a))') trim(d1_ptr%name), ';', lbound(d1_ptr%d, 1), ';', ubound(d1_ptr%d, 1)

!----------------------------------------------------------------------
! Individual datum info.
! Input syntax:
!   python data1 <ix_universe>@<d2_name>.<d1_datum>[<dat_index>]
! Use the "python data-d1" command to get detailed info on a specific d1 array.
! Output syntax is variable list form. See documentation at beginning of this file.
! Example:
!   python data_d1 1@orbit.x[10]

case ('data1')

  call tao_find_data (err, line, d_array = d_array)

  if (.not. allocated(d_array) .or. size(d_array) /= 1) then
    nl=nl+1; li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid datum name.')
    call end_stuff()
    return
  endif

  d_ptr => d_array(1)%d

  nl=nl+1; write (li(nl), amt) 'ele_name;STR;T;',                         d_ptr%ele_name
  nl=nl+1; write (li(nl), amt) 'ele_start_name;STR;T;',                   d_ptr%ele_start_name
  nl=nl+1; write (li(nl), amt) 'ele_ref_name;STR;T;',                     d_ptr%ele_ref_name  
  nl=nl+1; write (li(nl), amt) 'data_type;STR;T;',                        d_ptr%data_type     
  nl=nl+1; write (li(nl), amt) 'merit_type;STR;T;',                       d_ptr%merit_type    
  nl=nl+1; write (li(nl), amt) 'data_source;STR;T;',                      d_ptr%data_source   
  nl=nl+1; write (li(nl), imt) 'ix_bunch;INT;T;',                         d_ptr%ix_bunch      
  nl=nl+1; write (li(nl), imt) 'ix_branch;INT;T;',                        d_ptr%ix_branch     
  nl=nl+1; write (li(nl), imt) 'ix_ele;INT;T;',                           d_ptr%ix_ele        
  nl=nl+1; write (li(nl), imt) 'ix_ele_start;INT;T;',                     d_ptr%ix_ele_start  
  nl=nl+1; write (li(nl), imt) 'ix_ele_ref;INT;T;',                       d_ptr%ix_ele_ref    
  nl=nl+1; write (li(nl), imt) 'ix_ele_merit;INT;F;',                     d_ptr%ix_ele_merit  
  nl=nl+1; write (li(nl), imt) 'ix_d1;INT;F;',                            d_ptr%ix_d1         
  nl=nl+1; write (li(nl), imt) 'ix_data;INT;F;',                          d_ptr%ix_data       
  nl=nl+1; write (li(nl), imt) 'ix_dmodel;INT;F;',                        d_ptr%ix_dModel     
  nl=nl+1; write (li(nl), rmt) 'meas_value;REAL;T;',                      d_ptr%meas_value    
  nl=nl+1; write (li(nl), rmt) 'ref_value;REAL;T;',                       d_ptr%ref_value     
  nl=nl+1; write (li(nl), rmt) 'model_value;REAL;F;',                     d_ptr%model_value   
  nl=nl+1; write (li(nl), rmt) 'design_value;REAL;F;',                    d_ptr%design_value  
  nl=nl+1; write (li(nl), rmt) 'old_value;REAL;F;',                       d_ptr%old_value     
  nl=nl+1; write (li(nl), rmt) 'base_value;REAL;F;',                      d_ptr%base_value    
  nl=nl+1; write (li(nl), rmt) 'delta_merit;REAL;F;',                     d_ptr%delta_merit   
  nl=nl+1; write (li(nl), rmt) 'weight;REAL;T;',                          d_ptr%weight        
  nl=nl+1; write (li(nl), rmt) 'invalid_value;REAL;F;',                   d_ptr%invalid_value 
  nl=nl+1; write (li(nl), rmt) 'merit;REAL;F;',                           d_ptr%merit         
  nl=nl+1; write (li(nl), rmt) 's;REAL;F;',                               d_ptr%s             
  nl=nl+1; write (li(nl), lmt) 'exists;LOGIC;F;',                         d_ptr%exists        
  nl=nl+1; write (li(nl), lmt) 'good_model;LOGIC;F;',                     d_ptr%good_model
  nl=nl+1; write (li(nl), lmt) 'good_base;LOGIC;F;',                      d_ptr%good_base
  nl=nl+1; write (li(nl), lmt) 'good_design;LOGIC;F;',                    d_ptr%good_design
  nl=nl+1; write (li(nl), lmt) 'good_meas;LOGIC;T;',                      d_ptr%good_meas
  nl=nl+1; write (li(nl), lmt) 'good_ref;LOGIC;T;',                       d_ptr%good_ref
  nl=nl+1; write (li(nl), lmt) 'good_user;LOGIC;T;',                      d_ptr%good_user
  nl=nl+1; write (li(nl), lmt) 'good_opt;LOGIC;T;',                       d_ptr%good_opt
  nl=nl+1; write (li(nl), lmt) 'good_plot;LOGIC;T;',                      d_ptr%good_plot
  nl=nl+1; write (li(nl), lmt) 'useit_plot;LOGIC;F;',                     d_ptr%useit_plot
  nl=nl+1; write (li(nl), lmt) 'useit_opt;LOGIC;F;',                      d_ptr%useit_opt

!----------------------------------------------------------------------
! Lattice element list.
! Input syntax:
!   python lat_general <ix_universe>

case ('lat_general')

  u => point_to_uni(.false., err); if (err) return
  
  lat => u%design%lat
  do i = 0, ubound(lat%branch, 1)
    branch => lat%branch(i)
    nl=nl+1; write (li(nl), '(i0, 3a, 2(i0, a))') i, ';', trim(branch%name), ';', branch%n_ele_track, ';', branch%n_ele_max
  enddo

!----------------------------------------------------------------------
! Lattice element list.
! Input syntax:
!   python branch1 <ix_universe>@<ix_branch>

case ('branch1')

  u => point_to_uni(.false., err); if (err) return
  ix_branch = parse_branch(.false., err); if (err) return
  branch => u%design%lat%branch(ix_branch)

  nl=nl+1; write (li(nl), amt) 'name;STR;F;',                              branch%name
  nl=nl+1; write (li(nl), imt) 'ix_branch;INT;F;',                         branch%ix_branch
  nl=nl+1; write (li(nl), imt) 'ix_from_branch;INT;F;',                    branch%ix_from_branch
  nl=nl+1; write (li(nl), imt) 'ix_from_ele;INT;F;',                       branch%ix_from_ele

  nl=nl+1; write (li(nl), rmt) 'param.n_part;REAL;F;',                           branch%param%n_part
  nl=nl+1; write (li(nl), rmt) 'param.total_length;REAL;F;',                     branch%param%total_length
  nl=nl+1; write (li(nl), rmt) 'param.unstable_factor;REAL;F;',                  branch%param%unstable_factor
  nl=nl+1; write (li(nl), amt) 'param.particle;STR;F;',                          species_name(branch%param%particle)
  nl=nl+1; write (li(nl), imt) 'param.default_tracking_species;INT;F;',          branch%param%default_tracking_species
  nl=nl+1; write (li(nl), imt) 'param.geometry;INT;F;',                          branch%param%geometry
  nl=nl+1; write (li(nl), imt) 'param.ixx;INT;F;',                               branch%param%ixx
  nl=nl+1; write (li(nl), lmt) 'param.stable;LOGIC;F;',                          branch%param%stable
  nl=nl+1; write (li(nl), lmt) 'param.backwards_time_tracking;LOGIC;F;',         branch%param%backwards_time_tracking


!----------------------------------------------------------------------
! Lattice element list.
! Input syntax:
!   python lat_ele <branch_name>
! <branch_name> should have the form:
!   <ix_uni>@<ix_branch>

case ('lat_ele_list')

  u => point_to_uni(.true., err); if (err) return
  ix_branch = parse_branch(.false., err); if (err) return
  branch => u%design%lat%branch(ix_branch)

  do i = 0, branch%n_ele_max
    nl=nl+1; write (li(nl), '(i0, 2a)') i, ';', branch%ele(i)%name
  enddo

!----------------------------------------------------------------------
! parameters associated with given lattice element. 
! Input syntax: 
!   python lat_ele1 ix_universe@ix_branch>>ix_ele|which who
! where "which" is one of:
!   model
!   base
!   design
! and "who" is one of:
!   general         ! ele%xxx compnents where xxx is "simple" component (not a structure nor an array, nor allocatable, nor pointer).
!   parameters      ! parameters in ele%value array
!   multipole       ! nonzero multipole components.
!   floor           ! floor coordinates.
!   twiss           ! twiss parameters at exit end.
!   orbit           ! orbit at exit end.
! Example:
!   python lat_ele1 1@0>>547|design twiss

case ('lat_ele1')

  ix = index(line, ' ')
  call string_trim(line(ix:), who, ix2)
  who = line(1:ix)

  u => point_to_uni(.true., err); if (err) return
  tao_lat => point_to_tao_lat(err); if (err) return
  ele => point_to_ele(err); if (err) return

  select case (who)
  case ('general')
    nl=nl+1; write (li(nl), amt) 'name;STR;T;',                             ele%name
    nl=nl+1; write (li(nl), amt) 'type;STR;T;',                             ele%type
    nl=nl+1; write (li(nl), amt) 'alias;STR;T;',                            ele%alias
    nl=nl+1; write (li(nl), amt) 'component_name;STR;F;',                   ele%component_name
    nl=nl+1; write (li(nl), rmt) 'gamma_c;REAL;F;',                         ele%gamma_c
    nl=nl+1; write (li(nl), rmt) 's;REAL;F;',                               ele%s
    nl=nl+1; write (li(nl), rmt) 'ref_time;REAL;F;',                        ele%ref_time
    nl=nl+1; write (li(nl), amt) 'key;STR;F;',                              key_name(ele%key)
    nl=nl+1; write (li(nl), amt) 'sub_key;STR;F;',                          sub_key_name(ele%sub_key)
    nl=nl+1; write (li(nl), imt) 'ix_ele;INT;F;',                           ele%ix_ele
    nl=nl+1; write (li(nl), imt) 'ix_branch;INT;F;',                        ele%ix_branch
    nl=nl+1; write (li(nl), amt) 'slave_status;STR;F;',                     control_name(ele%slave_status)
    nl=nl+1; write (li(nl), imt) 'n_slave;INT;F;',                          ele%n_slave
    nl=nl+1; write (li(nl), imt) 'n_slave_field;INT;F;',                    ele%n_slave_field
    nl=nl+1; write (li(nl), amt) 'lord_status;STR;F;',                      control_name(ele%lord_status)
    nl=nl+1; write (li(nl), imt) 'n_lord;INT;F;',                           ele%n_lord
    nl=nl+1; write (li(nl), imt) 'n_lord_field;INT;F;',                     ele%n_lord_field
    nl=nl+1; write (li(nl), amt) 'mat6_calc_method;STR;T;',                 mat6_calc_method_name(ele%mat6_calc_method)
    nl=nl+1; write (li(nl), amt) 'tracking_method;STR;T;',                  tracking_method_name(ele%tracking_method)
    nl=nl+1; write (li(nl), amt) 'spin_tracking_method;STR;T;',             spin_tracking_method_name(ele%spin_tracking_method)
    nl=nl+1; write (li(nl), amt) 'ptc_integration_type;STR;T;',             ptc_integration_type_name(ele%ptc_integration_type)
    nl=nl+1; write (li(nl), amt) 'field_calc;STR;T;',                       field_calc_name(ele%field_calc)
    nl=nl+1; write (li(nl), amt) 'aperture_at;STR;T;',                      aperture_at_name(ele%aperture_at)
    nl=nl+1; write (li(nl), amt) 'aperture_type;STR;T;',                    aperture_type_name(ele%aperture_type)
    nl=nl+1; write (li(nl), imt) 'orientation;INT;T;',                      ele%orientation
    nl=nl+1; write (li(nl), lmt) 'symplectify;LOGIC;T;',                    ele%symplectify
    nl=nl+1; write (li(nl), lmt) 'mode_flip;LOGIC;F;',                      ele%mode_flip
    nl=nl+1; write (li(nl), lmt) 'multipoles_on;LOGIC;T;',                  ele%multipoles_on
    nl=nl+1; write (li(nl), lmt) 'scale_multipoles;LOGIC;T;',               ele%scale_multipoles
    nl=nl+1; write (li(nl), lmt) 'taylor_map_includes_offsets;LOGIC;T;',    ele%taylor_map_includes_offsets
    nl=nl+1; write (li(nl), lmt) 'field_master;LOGIC;T;',                   ele%field_master
    nl=nl+1; write (li(nl), lmt) 'is_on;LOGIC;T;',                          ele%is_on
    nl=nl+1; write (li(nl), lmt) 'csr_calc_on;LOGIC;T;',                    ele%csr_calc_on
    nl=nl+1; write (li(nl), lmt) 'offset_moves_aperture;LOGIC;T;',          ele%offset_moves_aperture

  case ('parameters')
    do i = 1, num_ele_attrib$
      attrib = attribute_info(ele, i)
      a_name = attrib%name
      if (a_name == null_name$) cycle
      if (attrib%type == private$) cycle
      free = attribute_free (ele, a_name, .false.)

      select case (attribute_type(a_name))
      case (is_logical$)
        nl=nl+1; write (li(nl), '(2a, l1, 2a)') trim(a_name), ';LOGIC;', free, ';', is_true(ele%value(i))
      case (is_integer$)
        nl=nl+1; write (li(nl), '(2a, l1, a, i0)') trim(a_name), ';INT;', free, ';', nint(ele%value(i))
      case (is_real$)
        nl=nl+1; write (li(nl), '(2a, l1, a, es21.13)') trim(a_name), ';REAL;', free, ';', ele%value(i)
      case (is_switch$)
        name = switch_attrib_value_name (a_name, ele%value(i), ele)
        nl=nl+1; write (li(nl), '(2a, l1, 2a)')  trim(a_name), ';STR;', free, ';', trim(name)
      end select
    enddo

  case ('multipole')
    if (associated(ele%a_pole)) then
      do i = 0, ubound(ele%a_pole, 1)
        if (ele%a_pole(i) /= 0) then
          nl=nl+1; write (li(nl), vrmt) 'a_pole[', i, '];REAL;T;', ele%a_pole(i) 
        endif
        if (ele%b_pole(i) /= 0) then
          nl=nl+1; write (li(nl), vrmt) 'b_pole[', i, '];REAL;T;', ele%b_pole(i) 
        endif
      enddo
    endif

    if (associated(ele%a_pole_elec)) then
      do i = 0, ubound(ele%a_pole_elec, 1)
        if (ele%a_pole_elec(i) /= 0) then
          nl=nl+1; write (li(nl), vrmt) 'a_pole_elec[', i, '];REAL;T;', ele%a_pole_elec(i) 
        endif
        if (ele%b_pole_elec(i) /= 0) then
          nl=nl+1; write (li(nl), vrmt) 'b_pole_elec[', i, '];REAL;T;', ele%b_pole_elec(i) 
        endif
      enddo
    endif

  case ('floor')
    nl=nl+1; write (li(nl), '(3(es21.13, a))') ele%floor%r(1), ';',ele%floor%r(2), ';', ele%floor%r(3) 
    nl=nl+1; write (li(nl), '(3(es21.13, a))') ele%floor%theta, ';',ele%floor%phi, ';', ele%floor%psi

  case ('twiss')
    call twiss_out (ele%a, 'a')
    call twiss_out (ele%b, 'b')

  case ('orbit')
    call orbit_out (tao_lat%lat_branch(ele%ix_branch)%orbit(ele%ix_ele))

  case default
    nl=nl+1; li(nl) = 'INVALID'
    call out_io (s_error$, r_name, 'python lat_ele1 <ele>|<which> <who>: Bad <who>: ' // who)
    return
  end select  

!----------------------------------------------------------------------
! Beam initialization parameters.
! Input syntax:
!   python beam_init ix_universe

case ('beam_init')

  u => point_to_uni(.true., err); if (err) return
  beam_init => u%beam%beam_init

  nl=nl+1; write (li(nl), amt) 'file_name;STR;F;',                         beam_init%file_name
  nl=nl+1; write (li(nl), rmt) 'sig_z_jitter;REAL;T;',                     beam_init%sig_z_jitter
  nl=nl+1; write (li(nl), rmt) 'sig_e_jitter;REAL;T;',                     beam_init%sig_e_jitter
  nl=nl+1; write (li(nl), imt) 'n_particle;INT;T;',                        beam_init%n_particle
  nl=nl+1; write (li(nl), lmt) 'renorm_center;LOGIC;T;',                   beam_init%renorm_center
  nl=nl+1; write (li(nl), lmt) 'renorm_sigma;LOGIC;T;',                    beam_init%renorm_sigma
  nl=nl+1; write (li(nl), amt) 'random_engine;STR;T;',                     beam_init%random_engine
  nl=nl+1; write (li(nl), amt) 'random_gauss_converter;STR;T;',            beam_init%random_gauss_converter
  nl=nl+1; write (li(nl), rmt) 'random_sigma_cutoff;REAL;T;',              beam_init%random_sigma_cutoff
  nl=nl+1; write (li(nl), rmt) 'a_norm_emit;REAL;T;',                      beam_init%a_norm_emit
  nl=nl+1; write (li(nl), rmt) 'b_norm_emit;REAL;T;',                      beam_init%b_norm_emit
  nl=nl+1; write (li(nl), rmt) 'a_emit;REAL;T;',                           beam_init%a_emit
  nl=nl+1; write (li(nl), rmt) 'b_emit;REAL;T;',                           beam_init%b_emit
  nl=nl+1; write (li(nl), rmt) 'dpz_dz;REAL;T;',                           beam_init%dPz_dz
  nl=nl+1; write (li(nl), rmt) 'dt_bunch;REAL;T;',                         beam_init%dt_bunch
  nl=nl+1; write (li(nl), rmt) 'sig_z;REAL;T;',                            beam_init%sig_z
  nl=nl+1; write (li(nl), rmt) 'sig_e;REAL;T;',                            beam_init%sig_e
  nl=nl+1; write (li(nl), rmt) 'bunch_charge;REAL;T;',                     beam_init%bunch_charge
  nl=nl+1; write (li(nl), imt) 'n_bunch;INT;T;',                           beam_init%n_bunch
  nl=nl+1; write (li(nl), amt) 'species;STR;T;',                           species_name(beam_init%species)
  nl=nl+1; write (li(nl), lmt) 'init_spin;LOGIC;T;',                       beam_init%init_spin
  nl=nl+1; write (li(nl), lmt) 'full_6d_coupling_calc;LOGIC;T;',           beam_init%full_6D_coupling_calc
  nl=nl+1; write (li(nl), lmt) 'use_lattice_center;LOGIC;T;',              beam_init%use_lattice_center
  nl=nl+1; write (li(nl), lmt) 'use_t_coords;LOGIC;T;',                    beam_init%use_t_coords
  nl=nl+1; write (li(nl), lmt) 'use_z_as_t;LOGIC;T;',                      beam_init%use_z_as_t

!----------------------------------------------------------------------
! Bunch parameters at the exit end of a given lattice element.
! Input syntax:
!   python bunch1 ix_universe@ix_branch>>ix_ele|which
! where "which" is one of:
!   model
!   base
!   design

case ('bunch1')  

  u => point_to_uni(.true., err); if (err) return
  tao_lat => point_to_tao_lat(err); if (err) return
  ele => point_to_ele(err); if (err) return

  bunch_params => tao_lat%lat_branch(ele%ix_branch)%bunch_params(ele%ix_ele)

  call twiss_out(bunch_params%x, 'x')
  call twiss_out(bunch_params%y, 'y')
  call twiss_out(bunch_params%z, 'z')
  call twiss_out(bunch_params%a, 'a')
  call twiss_out(bunch_params%b, 'b')
  call twiss_out(bunch_params%c, 'c')

  nl=nl+1; write (li(nl), rmt) 's;REAL;F;',                                bunch_params%s
  nl=nl+1; write (li(nl), rmt) 'charge_live;REAL;F;',                      bunch_params%charge_live
  nl=nl+1; write (li(nl), imt) 'n_particle_tot;INT;F;',                    bunch_params%n_particle_tot
  nl=nl+1; write (li(nl), imt) 'n_particle_live;INT;F;',                   bunch_params%n_particle_live
  nl=nl+1; write (li(nl), imt) 'n_particle_lost_in_ele;INT;F;',            bunch_params%n_particle_lost_in_ele

!----------------------------------------------------------------------
! Twiss at given s position
! Input syntax:
!   python twiss_at_s ix_uni@ix_branch>>s|which
! where "which" is one of:
!   model
!   base
!   design

case ('twiss_at_s')

  u => point_to_uni(.true., err); if (err) return
  tao_lat => point_to_tao_lat(err); if (err) return
  ix_branch = parse_branch(.true., err); if (err) return
  s_pos = parse_real(err); if (err) return

  call twiss_and_track_at_s (tao_lat%lat, s_pos, this_ele, tao_lat%lat_branch(ix_branch)%orbit, ix_branch = ix_branch)
  call twiss_out (this_ele%a, 'a')
  call twiss_out (this_ele%b, 'b')

!----------------------------------------------------------------------
! Twiss at given s position
! Input syntax:
!   python orbit_at_s ix_uni@ix_branch>>s|which
! where "which" is one of:
!   model
!   base
!   design

case ('orbit_at_s')

  u => point_to_uni(.true., err); if (err) return
  tao_lat => point_to_tao_lat(err); if (err) return
  ix_branch = parse_branch(.true., err); if (err) return
  s_pos = parse_real(err); if (err) return

  call twiss_and_track_at_s (tao_lat%lat, s_pos, orb = tao_lat%lat_branch(ix_branch)%orbit, orb_at_s = orb, ix_branch = ix_branch)
  call orbit_out (orb)

!----------------------------------------------------------------------

case default

  call out_io (s_error$, r_name, "python command internal error, shouldn't be here!")

end select

call end_stuff()

!----------------------------------------------------------------------
! return through scratch

contains

subroutine end_stuff()

scratch%n_lines = nl

if (doprint) call out_io (s_blank$, r_name, li(1:nl))

if (opened) then
  do i = 1, nl
    write (iu_write, '(a)') trim(li(i))
  enddo
  close (iu_write)
endif

end subroutine

!----------------------------------------------------------------------
! contains

function point_to_uni (has_ampersand, err) result (u)

type (tao_universe_struct), pointer :: u
logical has_ampersand, err
character(40) str

err = .false.

if (has_ampersand) then
  ix = index(line, '@')
  if (ix == 0) then
    nl=nl+1; li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Missing "@"')
    call end_stuff()
    err = .true.
    return
  endif
  str = line(1:ix-1)
  line = line(ix+1:)
else
  str = line
endif

read (str, *,  iostat = ios)  ix_universe
if (ios /= 0) ix_universe = -999

u => tao_pointer_to_universe(ix_universe)

if (.not. associated(u)) then
  nl=nl+1; li(nl) = 'INVALID'
  call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": bad universe index')
  call end_stuff()
  err = .true.
endif

end function point_to_uni

!----------------------------------------------------------------------
! contains

subroutine re_allocate_lines (n_lines)

integer n_lines

if (.not. allocated(ss%lines)) allocate (ss%lines(n_lines))
if (size(ss%lines) < n_lines) call re_allocate (ss%lines, n_lines)

li => ss%lines

end subroutine re_allocate_lines

!----------------------------------------------------------------------
! contains

function point_to_tao_lat (err) result (tao_lat)

type (tao_lattice_struct), pointer :: tao_lat
logical err

err = .false.

call string_trim(line, line, ix)
call string_trim(line(ix+1:), who, i)
line = line(1:ix)

ix = index(line, '|')
if (ix == 0) then
  nl=nl+1; li(nl) = 'INVALID'
  call out_io (s_error$, r_name, '"python ' // trim(command) // '" Expecting "|" character')
  err = .true.
  return
endif

select case (line(ix+1:))
case ('model')
  tao_lat => u%model
case ('base')
  tao_lat => u%base
case ('design')
  tao_lat => u%design
case default
  nl=nl+1; li(nl) = 'INVALID'
  call out_io (s_error$, r_name, 'python ' // trim(input_str) //  ': Expecting "|<which>" where <which> must be one of "model", "base", or "design"')
  err = .true.
end select

line = line(1:ix-1)

end function point_to_tao_lat

!----------------------------------------------------------------------
! contains

function point_to_ele (err) result (ele)

type (ele_struct), pointer :: ele
logical err

!

call lat_ele_locator (line, tao_lat%lat, eles, n_loc)

if (n_loc /= 1) then
  nl=nl+1; li(nl) = 'INVALID'
  call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Cannot locate element.')
  return
endif

ele => eles(1)%ele

end function point_to_ele

!----------------------------------------------------------------------
! contains

function parse_branch (has_separator, err) result (ix_branch)

integer ix_branch
logical has_separator, err
character(40) str

!

err = .false.

if (has_separator) then
  ix = index(line, '>>')

  if (ix == 0) then
    nl=nl+1; li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Missing ">>"')
    call end_stuff()
    err = .true.
    return
  endif

  str = line(1:ix-1)
  line = line(ix+2:)
else
  str = line
endif

read (str, *, iostat = ios) ix_branch
if (ios /= 0) ix_branch = -999

if (ix_branch < 0 .or. ix_branch > ubound(u%design%lat%branch, 1) .or. len_trim(str) == 0) then
  nl=nl+1; li(nl) = 'INVALID'
  call out_io (s_error$, r_name, '"python ' // trim(input_str) // '" missing or out of range branch index')
  call end_stuff()
  err = .true.
  return
endif

end function parse_branch

!----------------------------------------------------------------------
! contains

function parse_real (err) result (a_real)

real(rp) a_real
logical err

call string_to_real (line, real_garbage$, a_real, err)
if (err .or. a_real == real_garbage$) then
  nl=nl+1; li(nl) = 'INVALID'
  call out_io (s_error$, r_name, '"python ' // trim(input_str) // '" Bad real number')
  call end_stuff()
  return
endif

end function parse_real

!----------------------------------------------------------------------
! contains

subroutine orbit_out (orbit)

type (coord_struct) orbit

nl=nl+1; write (li(nl), rmt) 'x;REAL;F;',                                orbit%vec(1)
nl=nl+1; write (li(nl), rmt) 'px;REAL;F;',                               orbit%vec(2)
nl=nl+1; write (li(nl), rmt) 'y;REAL;F;',                                orbit%vec(3)
nl=nl+1; write (li(nl), rmt) 'py;REAL;F;',                               orbit%vec(4)
nl=nl+1; write (li(nl), rmt) 'z;REAL;F;',                                orbit%vec(5)
nl=nl+1; write (li(nl), rmt) 'pz;REAL;F;',                               orbit%vec(6)

nl=nl+1; write (li(nl), rmt) 'spin_x;REAL;F;',                           orbit%spin(1)
nl=nl+1; write (li(nl), rmt) 'spin_y;REAL;F;',                           orbit%spin(2)
nl=nl+1; write (li(nl), rmt) 'spin_z;REAL;F;',                           orbit%spin(3)

nl=nl+1; write (li(nl), rmt) 'field_x;REAL;F;',                          orbit%field(1)
nl=nl+1; write (li(nl), rmt) 'field_y;REAL;F;',                          orbit%field(2)

nl=nl+1; write (li(nl), rmt) 'phase_x;REAL;F;',                          orbit%phase(1)
nl=nl+1; write (li(nl), rmt) 'phase_y;REAL;F;',                          orbit%phase(2)

nl=nl+1; write (li(nl), rmt) 's;REAL;F;',                                orbit%s
nl=nl+1; write (li(nl), rmt) 't;REAL;F;',                                orbit%t
nl=nl+1; write (li(nl), rmt) 'charge;REAL;F;',                           orbit%charge
nl=nl+1; write (li(nl), rmt) 'path_len;REAL;F;',                         orbit%path_len
nl=nl+1; write (li(nl), rmt) 'p0c;REAL;F;',                              orbit%p0c
nl=nl+1; write (li(nl), rmt) 'beta;REAL;F;',                             orbit%beta
nl=nl+1; write (li(nl), imt) 'ix_ele;INT;F;',                            orbit%ix_ele
nl=nl+1; write (li(nl), amt) 'state;STR;F;',                             coord_state_name(orbit%state)
nl=nl+1; write (li(nl), imt) 'direction;INT;F;',                         orbit%direction
nl=nl+1; write (li(nl), amt) 'species;STR;F;',                           species_name(orbit%species)
nl=nl+1; write (li(nl), amt) 'location;STR;F;',                          location_name(orbit%location)

end subroutine orbit_out

!----------------------------------------------------------------------
! contains

subroutine twiss_out (twiss, suffix)

type (twiss_struct) twiss
character(*) suffix
character(20) fmt

fmt = '(3a, es21.13)'

nl=nl+1; write (li(nl), fmt) 'beta_', suffix, ';REAL;F;',                          twiss%beta
nl=nl+1; write (li(nl), fmt) 'alpha_', suffix, ';REAL;F;',                         twiss%alpha
nl=nl+1; write (li(nl), fmt) 'gamma_', suffix, ';REAL;F;',                         twiss%gamma
nl=nl+1; write (li(nl), fmt) 'phi_', suffix, ';REAL;F;',                           twiss%phi
nl=nl+1; write (li(nl), fmt) 'eta_', suffix, ';REAL;F;',                           twiss%eta
nl=nl+1; write (li(nl), fmt) 'etap_', suffix, ';REAL;F;',                          twiss%etap
nl=nl+1; write (li(nl), fmt) 'sigma_', suffix, ';REAL;F;',                         twiss%sigma
nl=nl+1; write (li(nl), fmt) 'sigma_p_', suffix, ';REAL;F;',                       twiss%sigma_p
nl=nl+1; write (li(nl), fmt) 'emit_', suffix, ';REAL;F;',                          twiss%emit
nl=nl+1; write (li(nl), fmt) 'norm_emit_', suffix, ';REAL;F;',                     twiss%norm_emit

end subroutine twiss_out

end subroutine tao_python_cmd
