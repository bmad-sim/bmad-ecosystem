!+
! Subroutine tao_python_cmd (input_str)
!
! Print information in a form easily parsed by a scripting program like python.
!
! Note: The syntax for "variable list form" is:
!   {component_name};{type};{variable};{component_value}
!
! {type} is one of:
!   STR
!   INT
!   REAL
!   LOGIC
!   ENUM
!
! {variable} indicates if the component can be varied. It is one of:
!   T
!   F
!
! Input:
!   input_str  -- Character(*): What to show.
!-

subroutine tao_python_cmd (input_str)

use tao_interface, dummy => tao_python_cmd
use tao_command_mod, only: tao_next_switch, tao_cmd_split
use tao_init_data_mod, only: tao_point_d1_to_data
use tao_init_variables_mod, only: tao_point_v1_to_var
use iso_c_binding
use tao_c_interface_mod, only: tao_c_interface_com
use location_encode_mod
use attribute_mod
use twiss_and_track_mod, only: twiss_and_track_at_s

implicit none

type (tao_universe_struct), pointer :: u
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_ptr
type (tao_data_struct), pointer :: d_ptr
type (tao_d2_data_struct), allocatable :: d2_temp(:)
type (tao_d1_data_struct), allocatable :: d1_temp(:)
type (tao_data_struct), allocatable :: d_temp(:)
type (tao_v1_var_array_struct), allocatable, save, target :: v1_array(:)
type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_var_struct), pointer :: v_ptr
type (tao_var_array_struct), allocatable, save, target :: v_array(:)
type (tao_v1_var_struct), allocatable :: v1_temp(:)
type (tao_var_struct), allocatable :: v_temp(:)
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
type (ele_struct), target :: this_ele
type (coord_struct), pointer :: orbit
type (coord_struct), target :: orb
type (bunch_params_struct), pointer :: bunch_params
type (bunch_params_struct), pointer :: bunch_p
type (ele_pointer_struct), allocatable, save :: eles(:)
type (branch_struct), pointer :: branch
type (tao_universe_branch_struct), pointer :: uni_branch
type (random_state_struct) ran_state
type (ele_attribute_struct) attrib

character(*) input_str
character(n_char_show), allocatable :: li(:) 
character(24) imt, rmt, lmt, amt, iamt, vamt, vrmt, vrmt2
character(40) max_loc, ele_name, name1(40), name2(40), a_name, name
character(200) line, file_name, all_who
character(20), allocatable :: name_list(:)
character(20) cmd, command, who, which, v_str
character(20) :: r_name = 'tao_python_cmd'
character(20) :: cmd_names(34)= [character(20) :: &
  'beam_init', 'branch1', 'bunch1', 'bmad_com', &
  'data_create', 'data_destroy', 'data_general', 'data_d2', 'data_d1', 'data1', &
  'enum', 'global', 'help', &
  'lat_ele_list', 'lat_ele1', 'lat_general', 'lat_list', 'lat_param_units', &
  'orbit_at_s', &
  'plot_list', 'plot1', 'plot_graph', 'plot_curve', 'plot_line', 'plot_symbol', &
  'species_to_int', 'species_to_str', 'twiss_at_s', 'universe', &
  'var_create', 'var_destroy', 'var_general', 'var_v1', 'var1']

real(rp) s_pos, value
real(rp), allocatable :: re_array(:)

integer :: i, j, k, ie, iu, nn, md, nl, ct, nl2, n, ix, ix2, iu_write, n1, n2, i1, i2
integer :: ix_ele, ix_ele1, ix_ele2, ix_branch, ix_d2, n_who
integer :: ios, n_loc, ix_line, n_d1, ix_min(20), ix_max(20), n_delta

logical :: err, print_flag, opened, doprint, free, matched, track_only, use_real_array_buffer

character(20) switch

!

line = input_str
doprint = .true.
opened = .false.
tao_c_interface_com%n_real = 0
tao_c_interface_com%n_int = 0

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
call string_trim(line(ix+1:), line, ix_line)

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
rmt = '(a,es24.16,a)'
lmt = '(a,l1,a)'
vamt = '(a, i0, 3a)'
vrmt = '(a, i0, a, es24.16)'
vrmt2 = '(a, i0, i0, a, es24.16)'

nl = 0
call re_allocate_lines (200)

select case (command)

!----------------------------------------------------------------------
! Beam initialization parameters.
! Command syntax:
!   python beam_init {ix_universe}
! where
!   {ix_universe} is a universe index.

case ('beam_init')

  u => point_to_uni(line, .false., err); if (err) return
  beam_init => u%beam%beam_init

  nl=incr(nl); write (li(nl), amt) 'file_name;STR;F;',                         beam_init%file_name
  nl=incr(nl); write (li(nl), rmt) 'sig_z_jitter;REAL;T;',                     beam_init%sig_z_jitter
  nl=incr(nl); write (li(nl), rmt) 'sig_e_jitter;REAL;T;',                     beam_init%sig_e_jitter
  nl=incr(nl); write (li(nl), imt) 'n_particle;INT;T;',                        beam_init%n_particle
  nl=incr(nl); write (li(nl), lmt) 'renorm_center;LOGIC;T;',                   beam_init%renorm_center
  nl=incr(nl); write (li(nl), lmt) 'renorm_sigma;LOGIC;T;',                    beam_init%renorm_sigma
  nl=incr(nl); write (li(nl), amt) 'random_engine;STR;T;',                     beam_init%random_engine
  nl=incr(nl); write (li(nl), amt) 'random_gauss_converter;STR;T;',            beam_init%random_gauss_converter
  nl=incr(nl); write (li(nl), rmt) 'random_sigma_cutoff;REAL;T;',              beam_init%random_sigma_cutoff
  nl=incr(nl); write (li(nl), rmt) 'a_norm_emit;REAL;T;',                      beam_init%a_norm_emit
  nl=incr(nl); write (li(nl), rmt) 'b_norm_emit;REAL;T;',                      beam_init%b_norm_emit
  nl=incr(nl); write (li(nl), rmt) 'a_emit;REAL;T;',                           beam_init%a_emit
  nl=incr(nl); write (li(nl), rmt) 'b_emit;REAL;T;',                           beam_init%b_emit
  nl=incr(nl); write (li(nl), rmt) 'dpz_dz;REAL;T;',                           beam_init%dPz_dz
  nl=incr(nl); write (li(nl), rmt) 'dt_bunch;REAL;T;',                         beam_init%dt_bunch
  nl=incr(nl); write (li(nl), rmt) 'sig_z;REAL;T;',                            beam_init%sig_z
  nl=incr(nl); write (li(nl), rmt) 'sig_e;REAL;T;',                            beam_init%sig_e
  nl=incr(nl); write (li(nl), rmt) 'bunch_charge;REAL;T;',                     beam_init%bunch_charge
  nl=incr(nl); write (li(nl), imt) 'n_bunch;INT;T;',                           beam_init%n_bunch
  nl=incr(nl); write (li(nl), amt) 'species;STR;T;',                           beam_init%species
  nl=incr(nl); write (li(nl), lmt) 'init_spin;LOGIC;T;',                       beam_init%init_spin
  nl=incr(nl); write (li(nl), lmt) 'full_6d_coupling_calc;LOGIC;T;',           beam_init%full_6D_coupling_calc
  nl=incr(nl); write (li(nl), lmt) 'use_particle_start_for_center;LOGIC;T;',   beam_init%use_particle_start_for_center
  nl=incr(nl); write (li(nl), lmt) 'use_t_coords;LOGIC;T;',                    beam_init%use_t_coords
  nl=incr(nl); write (li(nl), lmt) 'use_z_as_t;LOGIC;T;',                      beam_init%use_z_as_t

!----------------------------------------------------------------------
! Lattice element list.
! Command syntax:
!   python branch1 {ix_universe}@{ix_branch}
! where
!   {ix_universe} is a universe index
!   {ix_branch} is a lattice branch index

case ('branch1')

  u => point_to_uni(line, .true., err); if (err) return
  ix_branch = parse_branch(line, .false., err); if (err) return
  branch => u%design%lat%branch(ix_branch)

  nl=incr(nl); write (li(nl), amt) 'name;STR;F;',                              branch%name
  nl=incr(nl); write (li(nl), imt) 'ix_branch;INT;F;',                         branch%ix_branch
  nl=incr(nl); write (li(nl), imt) 'ix_from_branch;INT;F;',                    branch%ix_from_branch
  nl=incr(nl); write (li(nl), imt) 'ix_from_ele;INT;F;',                       branch%ix_from_ele

  nl=incr(nl); write (li(nl), rmt) 'param.n_part;REAL;F;',                           branch%param%n_part
  nl=incr(nl); write (li(nl), rmt) 'param.total_length;REAL;F;',                     branch%param%total_length
  nl=incr(nl); write (li(nl), rmt) 'param.unstable_factor;REAL;F;',                  branch%param%unstable_factor
  nl=incr(nl); write (li(nl), amt) 'param.particle;STR;F;',                          species_name(branch%param%particle)
  nl=incr(nl); write (li(nl), amt) 'param.default_tracking_species;INT;F;',          species_name(branch%param%default_tracking_species)
  nl=incr(nl); write (li(nl), imt) 'param.geometry;INT;F;',                          branch%param%geometry
  nl=incr(nl); write (li(nl), imt) 'param.ixx;INT;F;',                               branch%param%ixx
  nl=incr(nl); write (li(nl), lmt) 'param.stable;LOGIC;F;',                          branch%param%stable

!----------------------------------------------------------------------
! Bunch parameters at the exit end of a given lattice element.
! Command syntax:
!   python bunch1 {ix_universe}@{ix_branch}>>{ix_ele}|{which} coordinate
! where {which} is one of:
!   model
!   base
!   design
!
! Optional coordinate is one of:
! x, px, y, py, z, pz, 's', 't', 'charge', 'p0c'
! and will return an array. 

case ('bunch1')  

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, err, which, who); if (err) return
  ele => point_to_ele(line, err); if (err) return
  beam => u%uni_branch(ele%ix_branch)%ele(ele%ix_ele)%beam
  ! 
  
  !save_beam flag: u%uni_branch(<branch-index>)%ele(<ele-index>)%save_beam
  
  select case (who)
  case ('x', 'px', 'y', 'py', 'z', 'pz', 's', 't', 'charge', 'p0c') 
    call coord_out(beam, who)
    return
  case ('')
    bunch_params => tao_lat%tao_branch(ele%ix_branch)%bunch_params(ele%ix_ele)
    ! Continues below
  case default
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": coordinate not "x", "px", etc. ')
    call end_stuff()
    return
  end select

  call twiss_out(bunch_params%x, 'x', .true.)
  call twiss_out(bunch_params%y, 'y', .true.)
  call twiss_out(bunch_params%z, 'z', .true.)
  call twiss_out(bunch_params%a, 'a', .true.)
  call twiss_out(bunch_params%b, 'b', .true.)
  call twiss_out(bunch_params%c, 'c', .true.)

  ! Sigma matrix
  do i = 1, 6
    do j = 1,6
      nl=incr(nl); write (li(nl), vrmt2) 'sigma_', i, j, ';REAL;F;',           bunch_params%sigma(i,j)  
    enddo
  enddo

  ! Relative min, max, centroid
  do i = 1, 6
    nl=incr(nl); write (li(nl), vrmt) 'rel_min_', i, ';REAL;F;',               bunch_params%rel_min(i)
    nl=incr(nl); write (li(nl), vrmt) 'rel_max_', i, ';REAL;F;',               bunch_params%rel_max(i) 
    nl=incr(nl); write (li(nl), vrmt) 'centroid_vec_', i, ';REAL;F;',          bunch_params%centroid%vec(i) 
  enddo

  nl=incr(nl); write (li(nl), rmt) 'centroid_t;REAL;F;',                       bunch_params%centroid%t
  nl=incr(nl); write (li(nl), rmt) 'centroid_p0c;REAL;F;',                     bunch_params%centroid%p0c
  nl=incr(nl); write (li(nl), rmt) 'centroid_beta;REAL;F;',                    bunch_params%centroid%beta
  nl=incr(nl); write (li(nl), imt) 'ix_ele;INT;F;',                            bunch_params%centroid%ix_ele
  nl=incr(nl); write (li(nl), imt) 'direction;INT;F;',                         bunch_params%centroid%direction
  nl=incr(nl); write (li(nl), amt) 'species;STR;F;',                           species_name(bunch_params%centroid%species)
  nl=incr(nl); write (li(nl), amt) 'location;STR;F;',                          location_name(bunch_params%centroid%location)
  nl=incr(nl); write (li(nl), rmt) 's;REAL;F;',                                bunch_params%s
  nl=incr(nl); write (li(nl), rmt) 'charge_live;REAL;F;',                      bunch_params%charge_live
  nl=incr(nl); write (li(nl), imt) 'n_particle_tot;INT;F;',                    bunch_params%n_particle_tot
  nl=incr(nl); write (li(nl), imt) 'n_particle_live;INT;F;',                   bunch_params%n_particle_live
  nl=incr(nl); write (li(nl), imt) 'n_particle_lost_in_ele;INT;F;',            bunch_params%n_particle_lost_in_ele
  nl=incr(nl); write (li(nl), lmt) 'beam_saved;LOGIC;T;',                      allocated(beam%bunch)

!----------------------------------------------------------------------
! Bmad_com structure components
! Command syntax:
!   python bmad_com

case ('bmad_com')

  nl=incr(nl); write (li(nl), rmt) 'max_aperture_limit;REAL;T;',                 bmad_com%max_aperture_limit
  nl=incr(nl); write (li(nl), rmt) 'd_orb;REAL;T;',                              bmad_com%d_orb
  nl=incr(nl); write (li(nl), rmt) 'default_ds_step;REAL;T;',                    bmad_com%default_ds_step
  nl=incr(nl); write (li(nl), rmt) 'significant_length;REAL;T;',                 bmad_com%significant_length
  nl=incr(nl); write (li(nl), rmt) 'rel_tol_tracking;REAL;T;',                   bmad_com%rel_tol_tracking
  nl=incr(nl); write (li(nl), rmt) 'abs_tol_tracking;REAL;T;',                   bmad_com%abs_tol_tracking
  nl=incr(nl); write (li(nl), rmt) 'rel_tol_adaptive_tracking;REAL;T;',          bmad_com%rel_tol_adaptive_tracking
  nl=incr(nl); write (li(nl), rmt) 'abs_tol_adaptive_tracking;REAL;T;',          bmad_com%abs_tol_adaptive_tracking
  nl=incr(nl); write (li(nl), rmt) 'init_ds_adaptive_tracking;REAL;T;',          bmad_com%init_ds_adaptive_tracking
  nl=incr(nl); write (li(nl), rmt) 'min_ds_adaptive_tracking;REAL;T;',           bmad_com%min_ds_adaptive_tracking
  nl=incr(nl); write (li(nl), rmt) 'fatal_ds_adaptive_tracking;REAL;T;',         bmad_com%fatal_ds_adaptive_tracking
  nl=incr(nl); write (li(nl), rmt) 'autoscale_amp_abs_tol;REAL;T;',              bmad_com%autoscale_amp_abs_tol
  nl=incr(nl); write (li(nl), rmt) 'autoscale_amp_rel_tol;REAL;T;',              bmad_com%autoscale_amp_rel_tol
  nl=incr(nl); write (li(nl), rmt) 'autoscale_phase_tol;REAL;T;',                bmad_com%autoscale_phase_tol
  nl=incr(nl); write (li(nl), rmt) 'electric_dipole_moment;REAL;T;',             bmad_com%electric_dipole_moment
  nl=incr(nl); write (li(nl), rmt) 'ptc_cut_factor;REAL;T;',                     bmad_com%ptc_cut_factor
  nl=incr(nl); write (li(nl), rmt) 'sad_eps_scale;REAL;T;',                      bmad_com%sad_eps_scale
  nl=incr(nl); write (li(nl), rmt) 'sad_amp_max;REAL;T;',                        bmad_com%sad_amp_max
  nl=incr(nl); write (li(nl), imt) 'sad_n_div_max;INT;T;',                       bmad_com%sad_n_div_max
  nl=incr(nl); write (li(nl), imt) 'taylor_order;INT;T;',                        bmad_com%taylor_order
  nl=incr(nl); write (li(nl), imt) 'runge_kutta_order;INT;T;',                   bmad_com%runge_kutta_order
  nl=incr(nl); write (li(nl), imt) 'default_integ_order;INT;T;',                 bmad_com%default_integ_order
  nl=incr(nl); write (li(nl), imt) 'ptc_max_fringe_order;INT;T;',                bmad_com%ptc_max_fringe_order
  nl=incr(nl); write (li(nl), imt) 'max_num_runge_kutta_step;INT;T;',            bmad_com%max_num_runge_kutta_step
  nl=incr(nl); write (li(nl), lmt) 'rf_phase_below_transition_ref;LOGIC;T;',     bmad_com%rf_phase_below_transition_ref
  nl=incr(nl); write (li(nl), lmt) 'use_hard_edge_drifts;LOGIC;T;',              bmad_com%use_hard_edge_drifts
  nl=incr(nl); write (li(nl), lmt) 'sr_wakes_on;LOGIC;T;',                       bmad_com%sr_wakes_on
  nl=incr(nl); write (li(nl), lmt) 'lr_wakes_on;LOGIC;T;',                       bmad_com%lr_wakes_on
  nl=incr(nl); write (li(nl), lmt) 'mat6_track_symmetric;LOGIC;T;',              bmad_com%mat6_track_symmetric
  nl=incr(nl); write (li(nl), lmt) 'auto_bookkeeper;LOGIC;T;',                   bmad_com%auto_bookkeeper
  nl=incr(nl); write (li(nl), lmt) 'space_charge_on;LOGIC;T;',                   bmad_com%space_charge_on
  nl=incr(nl); write (li(nl), lmt) 'coherent_synch_rad_on;LOGIC;T;',             bmad_com%coherent_synch_rad_on
  nl=incr(nl); write (li(nl), lmt) 'spin_tracking_on;LOGIC;T;',                  bmad_com%spin_tracking_on
  nl=incr(nl); write (li(nl), lmt) 'backwards_time_tracking_on;LOGIC;T;',        bmad_com%backwards_time_tracking_on
  nl=incr(nl); write (li(nl), lmt) 'spin_sokolov_ternov_flipping_on;LOGIC;T;',   bmad_com%spin_sokolov_ternov_flipping_on
  nl=incr(nl); write (li(nl), lmt) 'radiation_damping_on;LOGIC;T;',              bmad_com%radiation_damping_on
  nl=incr(nl); write (li(nl), lmt) 'radiation_fluctuations_on;LOGIC;T;',         bmad_com%radiation_fluctuations_on
  nl=incr(nl); write (li(nl), lmt) 'conserve_taylor_maps;LOGIC;T;',              bmad_com%conserve_taylor_maps
  nl=incr(nl); write (li(nl), lmt) 'absolute_time_tracking_default;LOGIC;T;',    bmad_com%absolute_time_tracking_default
  nl=incr(nl); write (li(nl), lmt) 'convert_to_kinetic_momentum;LOGIC;T;',       bmad_com%convert_to_kinetic_momentum
  nl=incr(nl); write (li(nl), lmt) 'aperture_limit_on;LOGIC;T;',                 bmad_com%aperture_limit_on
  nl=incr(nl); write (li(nl), lmt) 'ptc_print_info_messages;LOGIC;T;',           bmad_com%ptc_print_info_messages
  nl=incr(nl); write (li(nl), lmt) 'debug;LOGIC;T;',                             bmad_com%debug

!----------------------------------------------------------------------
! List of datums in a given data d1 array.
! Command syntax:
!   python data_d1 {ix_universe}@{d2_name}.{d1_datum}
! Use the "python data_d2 {name}" command to get a list of d1 arrays. 
! Use the "python data1" command to get detailed information on a particular datum.
! Example:
!   python data_d1 1@orbit.x

case ('data_d1')

  call tao_find_data (err, line, d1_array = d1_array)

  if (err .or. .not. allocated(d1_array)) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid d1 data name.')
    call end_stuff()
    return
  endif

  d1_ptr => d1_array(1)%d1
  nl=incr(nl); write (li(nl), '(2a, 2(i0, a))') trim(d1_ptr%name), ';', lbound(d1_ptr%d, 1), ';', ubound(d1_ptr%d, 1)

!----------------------------------------------------------------------
! Create a d2 data structure along with associated d1 and data arrays.
!
! Command syntax:
!   python data_create {d2_name} {n_d1_data} {d_data_arrays_min_max}
! {d2_name} should be of the form {ix_uni}@{d2_datum_name}
! {n_d1_data} is the number of associated d1 data structures.
! {d_data_arrays_min_max} is an array of pairs of integers. The number of pairs is {n_d1_data}. 
!   The first number in the n^th pair gives the lower bound of the n^th d1 structure and the 
!   second number in the n^th pair gives the upper bound of the n^th d1 structure.
!
! The d1 structures created will be assigned initial names "1", "2", "3", etc.
!
! Example:
!   python data_create 2@orbit 2 0 45 1 47
! This example creates a d2 data structure called "orbit" with two d1 structures.
! The first d1 structure, assigned the name "1", has an associated data array with indexes in the range [0, 45].
! The second d1 structure, assigned the name "2", has an associated data arrray with indexes in the range [1, 47].
!
! Use the "set data" command to set a created datum parameters.
! Note: When setting multiple data parameters, temporarily toggle s%global%lattice_calc_on to False
!   ("set global lattice_calc_on = F") to prevent Tao trying to evaluate the partially created datum and
!   generating unwanted error messages.

case ('data_create')

  if (ix_line == 0) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": No d2 name given')
    call end_stuff()
    return
  endif

  name = line(1:ix_line)

  call string_trim (line(ix_line+1:), line, ix_line)
  if (.not. is_integer(line)) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Number of d1 arrays missing or invalid')
    call end_stuff()
    return
  endif
  read (line, *) n_d1
  call string_trim (line(ix_line+1:), line, ix_line)

  ix_min = 0; ix_max = 0 
  do i = 1, n_d1
    if (.not. is_integer(line)) exit
    read (line, *) ix_min(i)
    call string_trim (line(ix_line+1:), line, ix_line)
    if (.not. is_integer(line)) exit
    read (line, *) ix_max(i)
    call string_trim (line(ix_line+1:), line, ix_line)
  enddo

  if (ix_line /= 0 .or. i /= n_d1+1) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Malformed array of datum min/max for each d1.')
    call end_stuff()
    return
  endif

  ! Now create the d2 structure

  a_name = name
  ix = index(name, '@')
  if (ix == 0) then
    u => s%u(s%com%default_universe)
  else
    read (name(1:ix-1), *) iu
    u => s%u(iu)
    name = name(ix+1:)
  endif

  call tao_find_data(err, name, d2_array, print_err = .false.)
  if (size(d2_array) /= 0) then
    call destroy_this_data (a_name)
    call out_io (s_warn$, r_name, '"python ' // trim(input_str) // '": Data with this name already exists.', &
                                   'This old data has been destroyed to make room for the new data.')
!!    nl=incr(nl); li(nl) = 'INVALID'
!!    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Data with this name already exists.', &
!!                                   'Destroy with "python data_destroy" first.')
!!    call end_stuff()
!!    return
  endif


  if (allocated(u%d2_data)) then
    n2 = size(u%d2_data)
    if (u%n_d2_data_used + 1 > n2) then
      call move_alloc(u%d2_data, d2_temp)
      allocate (u%d2_data(n2+1))
      u%d2_data(1:n2) = d2_temp
    endif
  else
    allocate (u%d2_data(1))
  endif

  n_delta = sum(ix_max(1:n_d1)) - sum(ix_min(1:n_d1)) + n_d1

  if (allocated(u%d2_data)) then
    n = size(u%data)
    if (u%n_data_used + n_delta > n) then
      call move_alloc(u%data, d_temp)
      allocate (u%data(u%n_data_used + n_delta))
      u%data(1:u%n_data_used) = d_temp(1:u%n_data_used)
      do i = 1, size(u%data)
        u%data(i)%ix_data = i
      enddo
    endif
  else
    allocate (u%data(n_delta))
  endif

  i2 = 0   ! In case no d2 structures have yet been defined.

  do i = 1, u%n_d2_data_used
    n1 = size(u%d2_data(i)%d1)
    do j = 1, n1
      d1_ptr => u%d2_data(i)%d1(j)
      d1_ptr%d2 => u%d2_data(i)
      i1 = lbound(d1_ptr%d, 1)
      i1 = d1_ptr%d(i1)%ix_data
      i2 = ubound(d1_ptr%d, 1)
      i2 = d1_ptr%d(i2)%ix_data
      call tao_point_d1_to_data (d1_ptr, u%data(i1:i2), u%data(i1)%ix_d1)
    enddo
  enddo

  u%n_data_used = u%n_data_used + n_delta
  nn = u%n_d2_data_used + 1
  u%n_d2_data_used = nn
  u%d2_data(nn)%ix_d2_data = nn
  u%d2_data(nn)%name = name
  u%d2_data(nn)%ix_uni = iu
  if (allocated(u%d2_data(nn)%d1)) deallocate(u%d2_data(nn)%d1) ! Can happen if data has been destroyed.
  allocate (u%d2_data(nn)%d1(n_d1))

  do j = 1, n_d1
    d1_ptr => u%d2_data(nn)%d1(j)
    d1_ptr%d2 => u%d2_data(nn)
    write (d1_ptr%name, '(i0)') j
    i1 = i2 + 1
    i2 = i2 + 1 + ix_max(j) - ix_min(j)
    call tao_point_d1_to_data (d1_ptr, u%data(i1:i2), ix_min(j))
  enddo

!----------------------------------------------------------------------
! Destroy a d2 data structure along with associated d1 and data arrays.
! Command syntax:
!   python data_destroy {d2_datum}
! {d2_datum} should be of the form 
!   {ix_uni}@{d2_datum_name}

case ('data_destroy')

call destroy_this_data(line)

!----------------------------------------------------------------------
! List of d1 arrays in a given data d2.
! Command syntax:
!   python data_d2 {d2_datum}
! {d2_datum} should be of the form 
!   {ix_uni}@{d2_datum_name}

case ('data_d2')

  call tao_find_data (err, line, d2_array = d2_array)

  if (err .or. .not. allocated(d2_array)) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid d2 data name')
    call end_stuff()
    return
  endif

  d2_ptr => d2_array(1)%d2

  do i = lbound(d2_ptr%d1, 1), ubound(d2_ptr%d1, 1)
    nl=incr(nl); write (li(nl), '(a, i0, 2a)') 'd1[', i, '];STR;T;', d2_ptr%d1(i)%name
  enddo

  nl=incr(nl); write (li(nl), imt) 'ix_d2_data;INT;F;',                       d2_ptr%ix_d2_data
  nl=incr(nl); write (li(nl), amt) 'name;STR;T;',                             d2_ptr%name
  nl=incr(nl); write (li(nl), amt) 'data_file_name;STR;F;',                   d2_ptr%data_file_name
  nl=incr(nl); write (li(nl), amt) 'ref_file_name;STR;F;',                    d2_ptr%ref_file_name
  nl=incr(nl); write (li(nl), amt) 'data_date;STR;T;',                        d2_ptr%data_date
  nl=incr(nl); write (li(nl), amt) 'ref_date;STR;T;',                         d2_ptr%ref_date
  nl=incr(nl); write (li(nl), imt) 'ix_uni;INT;F;',                           d2_ptr%ix_uni
  nl=incr(nl); write (li(nl), imt) 'ix_data;INT;F;',                          d2_ptr%ix_data
  nl=incr(nl); write (li(nl), imt) 'ix_ref;INT;F;',                           d2_ptr%ix_ref
  nl=incr(nl); write (li(nl), lmt) 'data_read_in;LOGIC;F;',                   d2_ptr%data_read_in
  nl=incr(nl); write (li(nl), lmt) 'ref_read_in;LOGIC;F;',                    d2_ptr%ref_read_in

!----------------------------------------------------------------------
! Data d2 info for a given universe.
! Command syntax:
!   python data_general {ix_universe}

case ('data_general')

  u => point_to_uni(line, .false., err); if (err) return

  do i = 1, u%n_d2_data_used
    d2_ptr => u%d2_data(i)
    if (d2_ptr%name == '') cycle
    nl=incr(nl); write (li(nl), '(a)') d2_ptr%name
  enddo

!----------------------------------------------------------------------
! Individual datum info.
! Command syntax:
!   python data1 {ix_universe}@{d2_name}.{d1_datum}[{dat_index}]
! Use the "python data-d1" command to get detailed info on a specific d1 array.
! Output syntax is variable list form. See documentation at beginning of this file.
! Example:
!   python data_d1 1@orbit.x[10]

case ('data1')

  call tao_find_data (err, line, d_array = d_array)

  if (.not. allocated(d_array) .or. size(d_array) /= 1) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid datum name.')
    call end_stuff()
    return
  endif

  d_ptr => d_array(1)%d

  nl=incr(nl); write (li(nl), amt) 'ele_name;STR;T;',                         d_ptr%ele_name
  nl=incr(nl); write (li(nl), amt) 'ele_start_name;STR;T;',                   d_ptr%ele_start_name
  nl=incr(nl); write (li(nl), amt) 'ele_ref_name;STR;T;',                     d_ptr%ele_ref_name
  nl=incr(nl); write (li(nl), amt) 'data_type;STR;T;',                        d_ptr%data_type
  nl=incr(nl); write (li(nl), amt) 'merit_type;STR;T;',                       d_ptr%merit_type
  nl=incr(nl); write (li(nl), amt) 'data_source;STR;T;',                      d_ptr%data_source
  nl=incr(nl); write (li(nl), amt) 'eval_point;STR;T;',                       d_ptr%eval_point
  nl=incr(nl); write (li(nl), imt) 'ix_bunch;INT;T;',                         d_ptr%ix_bunch
  nl=incr(nl); write (li(nl), imt) 'ix_branch;INT;T;',                        d_ptr%ix_branch
  nl=incr(nl); write (li(nl), imt) 'ix_ele;INT;T;',                           d_ptr%ix_ele
  nl=incr(nl); write (li(nl), imt) 'ix_ele_start;INT;T;',                     d_ptr%ix_ele_start
  nl=incr(nl); write (li(nl), imt) 'ix_ele_ref;INT;T;',                       d_ptr%ix_ele_ref
  nl=incr(nl); write (li(nl), imt) 'ix_ele_merit;INT;F;',                     d_ptr%ix_ele_merit
  nl=incr(nl); write (li(nl), imt) 'ix_d1;INT;F;',                            d_ptr%ix_d1
  nl=incr(nl); write (li(nl), imt) 'ix_data;INT;F;',                          d_ptr%ix_data
  nl=incr(nl); write (li(nl), imt) 'ix_dmodel;INT;F;',                        d_ptr%ix_dModel
  nl=incr(nl); write (li(nl), rmt) 'meas_value;REAL;T;',                      d_ptr%meas_value
  nl=incr(nl); write (li(nl), rmt) 'ref_value;REAL;T;',                       d_ptr%ref_value
  nl=incr(nl); write (li(nl), rmt) 'model_value;REAL;F;',                     d_ptr%model_value
  nl=incr(nl); write (li(nl), rmt) 'design_value;REAL;F;',                    d_ptr%design_value
  nl=incr(nl); write (li(nl), rmt) 'old_value;REAL;F;',                       d_ptr%old_value
  nl=incr(nl); write (li(nl), rmt) 'base_value;REAL;F;',                      d_ptr%base_value
  nl=incr(nl); write (li(nl), rmt) 'delta_merit;REAL;F;',                     d_ptr%delta_merit
  nl=incr(nl); write (li(nl), rmt) 'weight;REAL;T;',                          d_ptr%weight
  nl=incr(nl); write (li(nl), rmt) 'invalid_value;REAL;F;',                   d_ptr%invalid_value
  nl=incr(nl); write (li(nl), rmt) 'merit;REAL;F;',                           d_ptr%merit
  nl=incr(nl); write (li(nl), rmt) 's;REAL;F;',                               d_ptr%s
  nl=incr(nl); write (li(nl), rmt) 's_offset;REAL;F;',                        d_ptr%s_offset
  nl=incr(nl); write (li(nl), lmt) 'exists;LOGIC;F;',                         d_ptr%exists
  nl=incr(nl); write (li(nl), lmt) 'good_model;LOGIC;F;',                     d_ptr%good_model
  nl=incr(nl); write (li(nl), lmt) 'good_base;LOGIC;F;',                      d_ptr%good_base
  nl=incr(nl); write (li(nl), lmt) 'good_design;LOGIC;F;',                    d_ptr%good_design
  nl=incr(nl); write (li(nl), lmt) 'good_meas;LOGIC;T;',                      d_ptr%good_meas
  nl=incr(nl); write (li(nl), lmt) 'good_ref;LOGIC;T;',                       d_ptr%good_ref
  nl=incr(nl); write (li(nl), lmt) 'good_user;LOGIC;T;',                      d_ptr%good_user
  nl=incr(nl); write (li(nl), lmt) 'good_opt;LOGIC;T;',                       d_ptr%good_opt
  nl=incr(nl); write (li(nl), lmt) 'good_plot;LOGIC;T;',                      d_ptr%good_plot
  nl=incr(nl); write (li(nl), lmt) 'useit_plot;LOGIC;F;',                     d_ptr%useit_plot
  nl=incr(nl); write (li(nl), lmt) 'useit_opt;LOGIC;F;',                      d_ptr%useit_opt

!----------------------------------------------------------------------
! List of possible values for enumerated numbers.
! Command syntax:
!   python enum {enum_name}
! Example:
!   python enum tracking_method

case ('enum')

  name = upcase(line)
  a_name = switch_attrib_value_name(name, 1.0_rp, this_ele, name_list = name_list)
  if (.not. allocated(name_list)) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid switch name.')
    call end_stuff()
    return
  endif

  n = 0
  do i = lbound(name_list, 1), ubound(name_list, 1)
    if (index(name_list(i), '!') == 0 .and. name_list(i) /= '') n = n + 1
  enddo

  do i = lbound(name_list, 1), ubound(name_list, 1)
    if (index(name_list(i), '!') /= 0 .or. name_list(i) == '') cycle
    nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(name_list(i))
  enddo

!----------------------------------------------------------------------
! Global parameters
! Command syntax: 
!   python global
! Output syntax is variable list form. See documentation at beginning of this file.

case ('global')

  nl=incr(nl); write (li(nl), rmt) 'lm_opt_deriv_reinit;REAL;T;',             s%global%lm_opt_deriv_reinit
  nl=incr(nl); write (li(nl), rmt) 'de_lm_step_ratio;REAL;T;',                s%global%de_lm_step_ratio
  nl=incr(nl); write (li(nl), rmt) 'de_var_to_population_factor;REAL;T;',     s%global%de_var_to_population_factor
  nl=incr(nl); write (li(nl), rmt) 'lmdif_eps;REAL;T;',                       s%global%lmdif_eps
  nl=incr(nl); write (li(nl), rmt) 'svd_cutoff;REAL;T;',                      s%global%svd_cutoff
  nl=incr(nl); write (li(nl), rmt) 'unstable_penalty;REAL;T;',                s%global%unstable_penalty
  nl=incr(nl); write (li(nl), rmt) 'merit_stop_value;REAL;T;',                s%global%merit_stop_value
  nl=incr(nl); write (li(nl), rmt) 'dmerit_stop_value;REAL;T;',               s%global%dmerit_stop_value
  nl=incr(nl); write (li(nl), rmt) 'random_sigma_cutoff;REAL;T;',             s%global%random_sigma_cutoff
  nl=incr(nl); write (li(nl), rmt) 'delta_e_chrom;REAL;T;',                   s%global%delta_e_chrom
  nl=incr(nl); write (li(nl), imt) 'n_opti_cycles;INT;T;',                    s%global%n_opti_cycles
  nl=incr(nl); write (li(nl), imt) 'n_opti_loops;INT;T;',                     s%global%n_opti_loops
  nl=incr(nl); write (li(nl), imt) 'phase_units;INT;T;',                      s%global%phase_units
  nl=incr(nl); write (li(nl), imt) 'bunch_to_plot;INT;T;',                    s%global%bunch_to_plot
  nl=incr(nl); write (li(nl), imt) 'random_seed;INT;T;',                      s%global%random_seed
  nl=incr(nl); write (li(nl), imt) 'n_top10_merit;INT;T;',                    s%global%n_top10_merit
  nl=incr(nl); write (li(nl), amt) 'random_engine;STR;T;',                    s%global%random_engine
  nl=incr(nl); write (li(nl), amt) 'random_gauss_converter;STR;T;',           s%global%random_gauss_converter
  nl=incr(nl); write (li(nl), amt) 'track_type;STR;T;',                       s%global%track_type
  nl=incr(nl); write (li(nl), amt) 'prompt_string;STR;T;',                    s%global%prompt_string
  nl=incr(nl); write (li(nl), amt) 'prompt_color;STR;T;',                     s%global%prompt_color
  nl=incr(nl); write (li(nl), amt) 'optimizer;STR;T;',                        s%global%optimizer
  nl=incr(nl); write (li(nl), amt) 'print_command;STR;T;',                    s%global%print_command
  nl=incr(nl); write (li(nl), amt) 'var_out_file;STR;T;',                     s%global%var_out_file
  nl=incr(nl); write (li(nl), lmt) 'opt_with_ref;LOGIC;T;',                   s%global%opt_with_ref
  nl=incr(nl); write (li(nl), lmt) 'opt_with_base;LOGIC;T;',                  s%global%opt_with_base
  nl=incr(nl); write (li(nl), lmt) 'label_lattice_elements;LOGIC;T;',         s%global%label_lattice_elements
  nl=incr(nl); write (li(nl), lmt) 'label_keys;LOGIC;T;',                     s%global%label_keys
  nl=incr(nl); write (li(nl), lmt) 'concatenate_maps;LOGIC;T;',               s%global%concatenate_maps
  nl=incr(nl); write (li(nl), lmt) 'derivative_recalc;LOGIC;T;',              s%global%derivative_recalc
  nl=incr(nl); write (li(nl), lmt) 'derivative_uses_design;LOGIC;T;',         s%global%derivative_uses_design
  nl=incr(nl); write (li(nl), lmt) 'orm_analysis;LOGIC;T;',                   s%global%orm_analysis
  nl=incr(nl); write (li(nl), lmt) 'plot_on;LOGIC;T;',                        s%global%plot_on
  nl=incr(nl); write (li(nl), lmt) 'lattice_calc_on;LOGIC;T;',                s%global%lattice_calc_on
  nl=incr(nl); write (li(nl), lmt) 'svd_retreat_on_merit_increase;LOGIC;T;',  s%global%svd_retreat_on_merit_increase
  nl=incr(nl); write (li(nl), lmt) 'stop_on_error;LOGIC;T;',                  s%global%stop_on_error
  nl=incr(nl); write (li(nl), lmt) 'command_file_print_on;LOGIC;T;',          s%global%command_file_print_on
  nl=incr(nl); write (li(nl), lmt) 'box_plots;LOGIC;T;',                      s%global%box_plots
  nl=incr(nl); write (li(nl), lmt) 'beam_timer_on;LOGIC;T;',                  s%global%beam_timer_on
  nl=incr(nl); write (li(nl), lmt) 'var_limits_on;LOGIC;T;',                  s%global%var_limits_on
  nl=incr(nl); write (li(nl), lmt) 'only_limit_opt_vars;LOGIC;T;',            s%global%only_limit_opt_vars
  nl=incr(nl); write (li(nl), lmt) 'optimizer_var_limit_warn;LOGIC;T;',       s%global%optimizer_var_limit_warn
  nl=incr(nl); write (li(nl), lmt) 'rf_on;LOGIC;T;',                          s%global%rf_on
  nl=incr(nl); write (li(nl), lmt) 'draw_curve_off_scale_warn;LOGIC;T;',      s%global%draw_curve_off_scale_warn
  nl=incr(nl); write (li(nl), lmt) 'wait_for_cr_in_single_mode;LOGIC;T;',     s%global%wait_for_CR_in_single_mode
  nl=incr(nl); write (li(nl), lmt) 'disable_smooth_line_calc;LOGIC;T;',       s%global%disable_smooth_line_calc
  nl=incr(nl); write (li(nl), lmt) 'debug_on;LOGIC;T;',                       s%global%debug_on


!----------------------------------------------------------------------
! returns list of "help xxx" topics
! Command syntax:
!   python help

case ('help')

  call tao_help ('help-list', '', li, n)

  nl2 = 0
  do i = 1, n
    if (li(i) == '') cycle
    call string_trim(li(i), line, ix)
    nl=incr(nl); name1(nl) = line(1:ix)
    call string_trim(line(ix+1:), line, ix)
    if (ix == 0) cycle
    nl2=nl2+1; name2(nl2) = line
  enddo

  li(1:nl) = name1(1:nl)
  li(nl+1:nl+nl2) = name2(1:nl2)
  nl = nl + nl2

!----------------------------------------------------------------------
! Lattice element list.
! Command syntax:
!   python lat_ele {branch_name}
! {branch_name} should have the form:
!   {ix_uni}@{ix_branch}

case ('lat_ele_list')

  u => point_to_uni(line, .true., err); if (err) return
  ix_branch = parse_branch(line, .false., err); if (err) return
  branch => u%design%lat%branch(ix_branch)

  call re_allocate_lines (branch%n_ele_max+100)

  do i = 0, branch%n_ele_max
    nl=incr(nl); write (li(nl), '(i0, 2a)') i, ';', branch%ele(i)%name
  enddo

!----------------------------------------------------------------------
! parameters associated with given lattice element. 
! Command syntax: 
!   python lat_ele1 {ix_universe}@{ix_branch}>>{ix_ele}|{which} {who}
! where {which} is one of:
!   model
!   base
!   design
! and {who} is one of:
!   general         ! ele%xxx compnents where xxx is "simple" component (not a structure nor an array, nor allocatable, nor pointer).
!   parameters      ! parameters in ele%value array
!   multipole       ! nonzero multipole components.
!   floor           ! floor coordinates.
!   twiss           ! twiss parameters at exit end.
!   orbit           ! orbit at exit end.
! Example:
!   python lat_ele1 1@0>>547|design twiss

case ('lat_ele1')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, err, which, who); if (err) return
  ele => point_to_ele(line, err); if (err) return

  select case (who)
  case ('general')
    nl=incr(nl); write (li(nl), amt) 'name;STR;T;',                             ele%name
    nl=incr(nl); write (li(nl), amt) 'type;STR;T;',                             ele%type
    nl=incr(nl); write (li(nl), amt) 'alias;STR;T;',                            ele%alias
    nl=incr(nl); write (li(nl), amt) 'component_name;STR;F;',                   ele%component_name
    nl=incr(nl); write (li(nl), rmt) 'gamma_c;REAL;F;',                         ele%gamma_c
    nl=incr(nl); write (li(nl), rmt) 's;REAL;F;',                               ele%s
    nl=incr(nl); write (li(nl), rmt) 'ref_time;REAL;F;',                        ele%ref_time
    nl=incr(nl); write (li(nl), amt) 'key;STR;F;',                              key_name(ele%key)
    nl=incr(nl); write (li(nl), amt) 'sub_key;STR;F;',                          sub_key_name(ele%sub_key)
    nl=incr(nl); write (li(nl), imt) 'ix_ele;INT;F;',                           ele%ix_ele
    nl=incr(nl); write (li(nl), imt) 'ix_branch;INT;F;',                        ele%ix_branch
    nl=incr(nl); write (li(nl), amt) 'slave_status;STR;F;',                     control_name(ele%slave_status)
    nl=incr(nl); write (li(nl), imt) 'n_slave;INT;F;',                          ele%n_slave
    nl=incr(nl); write (li(nl), imt) 'n_slave_field;INT;F;',                    ele%n_slave_field
    nl=incr(nl); write (li(nl), amt) 'lord_status;STR;F;',                      control_name(ele%lord_status)
    nl=incr(nl); write (li(nl), imt) 'n_lord;INT;F;',                           ele%n_lord
    nl=incr(nl); write (li(nl), imt) 'n_lord_field;INT;F;',                     ele%n_lord_field
    nl=incr(nl); write (li(nl), amt) 'mat6_calc_method;ENUM;T;',                mat6_calc_method_name(ele%mat6_calc_method)
    nl=incr(nl); write (li(nl), amt) 'tracking_method;ENUM;T;',                 tracking_method_name(ele%tracking_method)
    nl=incr(nl); write (li(nl), amt) 'spin_tracking_method;ENUM;T;',            spin_tracking_method_name(ele%spin_tracking_method)
    nl=incr(nl); write (li(nl), amt) 'ptc_integration_type;ENUM;T;',            ptc_integration_type_name(ele%ptc_integration_type)
    nl=incr(nl); write (li(nl), amt) 'field_calc;ENUM;T;',                      field_calc_name(ele%field_calc)
    nl=incr(nl); write (li(nl), amt) 'aperture_at;ENUM;T;',                     aperture_at_name(ele%aperture_at)
    nl=incr(nl); write (li(nl), amt) 'aperture_type;ENUM;T;',                   aperture_type_name(ele%aperture_type)
    nl=incr(nl); write (li(nl), imt) 'orientation;INT;T;',                      ele%orientation
    nl=incr(nl); write (li(nl), lmt) 'symplectify;LOGIC;T;',                    ele%symplectify
    nl=incr(nl); write (li(nl), lmt) 'mode_flip;LOGIC;F;',                      ele%mode_flip
    nl=incr(nl); write (li(nl), lmt) 'multipoles_on;LOGIC;T;',                  ele%multipoles_on
    nl=incr(nl); write (li(nl), lmt) 'scale_multipoles;LOGIC;T;',               ele%scale_multipoles
    nl=incr(nl); write (li(nl), lmt) 'taylor_map_includes_offsets;LOGIC;T;',    ele%taylor_map_includes_offsets
    nl=incr(nl); write (li(nl), lmt) 'field_master;LOGIC;T;',                   ele%field_master
    nl=incr(nl); write (li(nl), lmt) 'is_on;LOGIC;T;',                          ele%is_on
    nl=incr(nl); write (li(nl), lmt) 'csr_calc_on;LOGIC;T;',                    ele%csr_calc_on
    nl=incr(nl); write (li(nl), lmt) 'offset_moves_aperture;LOGIC;T;',          ele%offset_moves_aperture

  case ('parameters')
    do i = 1, num_ele_attrib$
      attrib = attribute_info(ele, i)
      a_name = attrib%name
      if (a_name == null_name$) cycle
      if (attrib%type == private$) cycle
      free = attribute_free (ele, a_name, .false.)
      if (which /= 'model') free = .false.

      select case (attribute_type(a_name))
      case (is_logical$)
        nl=incr(nl); write (li(nl), '(2a, l1, a, l1)') trim(a_name), ';LOGIC;', free, ';', is_true(ele%value(i))
      case (is_integer$)
        nl=incr(nl); write (li(nl), '(2a, l1, a, i0)') trim(a_name), ';INT;', free, ';', nint(ele%value(i))
      case (is_real$)
        nl=incr(nl); write (li(nl), '(2a, l1, a, es24.16)') trim(a_name), ';REAL;', free, ';', ele%value(i)
      case (is_switch$)
        name = switch_attrib_value_name (a_name, ele%value(i), ele)
        nl=incr(nl); write (li(nl), '(2a, l1, 2a)')  trim(a_name), ';ENUM;', free, ';', trim(name)
      end select
    enddo

  case ('multipole')
    if (which == 'model') then
      v_str = '];REAL;T;'
    else
      v_str = '];REAL;F;'
    endif

    if (associated(ele%a_pole)) then
      do i = 0, ubound(ele%a_pole, 1)
        if (ele%a_pole(i) /= 0) then
          nl=incr(nl); write (li(nl), vrmt) 'a_pole[', i, v_str, ele%a_pole(i) 
        endif
        if (ele%b_pole(i) /= 0) then
          nl=incr(nl); write (li(nl), vrmt) 'b_pole[', i, v_str, ele%b_pole(i) 
        endif
      enddo
    endif

    if (associated(ele%a_pole_elec)) then
      do i = 0, ubound(ele%a_pole_elec, 1)
        if (ele%a_pole_elec(i) /= 0) then
          nl=incr(nl); write (li(nl), vrmt) 'a_pole_elec[', i, v_str, ele%a_pole_elec(i) 
        endif
        if (ele%b_pole_elec(i) /= 0) then
          nl=incr(nl); write (li(nl), vrmt) 'b_pole_elec[', i, v_str, ele%b_pole_elec(i) 
        endif
      enddo
    endif

  case ('floor')
    nl=incr(nl); write (li(nl), '(3(es24.16, a))') ele%floor%r(1), ';',ele%floor%r(2), ';', ele%floor%r(3) 
    nl=incr(nl); write (li(nl), '(3(es24.16, a))') ele%floor%theta, ';',ele%floor%phi, ';', ele%floor%psi

  case ('twiss')
    free = attribute_free(ele, 'BETA_A', .false.) .and. (which == 'model')
    call twiss_out (ele%a, 'a', can_vary = free)
    call twiss_out (ele%b, 'b', can_vary = free)

  case ('orbit')
    call orbit_out (tao_lat%tao_branch(ele%ix_branch)%orbit(ele%ix_ele))

  case default
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, 'python lat_ele1 {ele}|{which} {who}: Bad {who}: ' // who)
    return
  end select  

!----------------------------------------------------------------------
! Lattice general
! Command syntax:
!   python lat_general {ix_universe}
!
! Output syntax:
!   branch_index;branch_name;n_ele_track;n_ele_max
case ('lat_general')

  u => point_to_uni(line, .false., err); if (err) return
  
  lat => u%design%lat
  do i = 0, ubound(lat%branch, 1)
    branch => lat%branch(i)
    nl=incr(nl); write (li(nl), '(i0, 3a, 2(i0, a))') i, ';', trim(branch%name), ';', branch%n_ele_track, ';', branch%n_ele_max
  enddo

!----------------------------------------------------------------------
! List of parameters at ends of lattice elements
! Command syntax:
!   python lat_list {ix_uni}@{ix_branch}>>{elements}|{which} {who}
! where 
!   {which} is one of:
!     model
!     base
!     design
!   {who} is a comma deliminated list of:
!     orbit.spin.1, orbit.spin.2, orbit.spin.3,
!     orbit.vec.1, orbit.vec.2, orbit.vec.3, orbit.vec.4, orbit.vec.5, orbit.vec.6,
!     orbit.t, orbit.beta,
!     orbit.state,     ! Note: state is an integer. alive$ = 1, anything else is lost.
!     orbit.energy, orbit.pc,
!     ele.a.beta, ele.a.alpha, ele.a.eta, ele.a.etap, ele.a.gamma, ele.a.phi,
!     ele.b.beta, ele.b.alpha, ele.b.eta, ele.b.etap, ele.b.gamma, ele.b.phi,
!     ele.x.eta, ele.x.etap,
!     ele.y.eta, ele.y.etap,
!     ele.s, ele.l
!     ele.e_tot, ele.p0c
!   {elements} is a string to match element names to. 
!     Use "*" to match to all elements.
!     Use the prefix "track:" to exclude lord elements.
! Note: To output through the real array buffer, add the prefix "real:" to {who}. In this
! case, {who} must only contain a single item
!
! Examples:
!   python lat_list 3@0>>track:Q*|base ele.s,orbit.vec.2
!   python lat_list 3@0>>Q*|base real:ele.s    ! Only a single item permitted with real buffer out.

case ('lat_list')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, err, who = all_who); if (err) return
  use_real_array_buffer = (all_who(1:5) == 'real:')
  if (use_real_array_buffer) then
    all_who = all_who(6:)
    call re_allocate(re_array, 1000)
  endif
  ix_branch = parse_branch(line, .true., err); if (err) return
  branch => u%design%lat%branch(ix_branch)
  ele_name = upcase(line)
  track_only = (ele_name(1:6) == 'TRACK:') 
  if (track_only) ele_name = ele_name(7:)

  n_who = 0
  do
    ix = index(all_who, ',')
    n_who = n_who + 1
    if (ix == 0) then
      name1(n_who) = all_who
      exit
    else
      name1(n_who) = all_who(1:ix-1)
      all_who = adjustl(all_who(ix+1:))
    endif
  enddo

  if (use_real_array_buffer .and. n_who /= 1) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, 'python lat_list: Number of "who" must be 1 for real buffered output.')
    return
  endif

  n = 0
  do ie = 1, branch%n_ele_max
    if (track_only .and. ie > branch%n_ele_track) cycle
    ele => branch%ele(ie)
    orbit => tao_lat%tao_branch(ix_branch)%orbit(ele%ix_ele)

    matched = match_ele_name(ele_name, ele, err); if (err) return
    if (.not. matched) cycle

    do i = 1, n_who
      select case (name1(i))
      case ('orbit.spin.1')
        value = orbit%spin(1)
      case ('orbit.spin.2')
        value = orbit%spin(2)
      case ('orbit.spin.3')
        value = orbit%spin(3)
      case ('orbit.vec.1')
        value = orbit%vec(1)
      case ('orbit.vec.2')
        value = orbit%vec(2)
      case ('orbit.vec.3')
        value = orbit%vec(3)
      case ('orbit.vec.4')
        value = orbit%vec(4)
      case ('orbit.vec.5')
        value = orbit%vec(5)
      case ('orbit.vec.6')
        value = orbit%vec(6)
      case ('orbit.t')
        value = orbit%t
      case ('orbit.beta')
        value = orbit%beta
      case ('orbit.state')
        value = orbit%state
      case ('orbit.energy')
        value = (1 + orbit%vec(6)) * orbit%p0c
      case ('orbit.pc')
        call convert_pc_to ((1 + orbit%vec(6)) * orbit%p0c, orbit%species, E_tot = value)
      case ('ele.a.beta')
        value = ele%a%beta
      case ('ele.a.alpha')
        value = ele%a%alpha
      case ('ele.a.eta')
        value = ele%a%eta
      case ('ele.a.etap')
        value = ele%a%etap
      case ('ele.a.gamma')
        value = ele%a%gamma
      case ('ele.a.phi')
        value = ele%a%phi
      case ('ele.b.beta')
        value = ele%b%beta
      case ('ele.b.alpha')
        value = ele%b%alpha
      case ('ele.b.eta')
        value = ele%b%eta
      case ('ele.b.etap')
        value = ele%b%etap
      case ('ele.b.gamma')
        value = ele%b%gamma
      case ('ele.b.phi')
        value = ele%b%phi
      case ('ele.e_tot')
        value = ele%value(e_tot$)
      case ('ele.p0c')
        value = ele%value(p0c$)
      case ('ele.x.eta')
        value = ele%x%eta
      case ('ele.x.etap')
        value = ele%x%etap
      case ('ele.y.eta')
        value = ele%y%eta
      case ('ele.y.etap')
        value = ele%y%etap
      case ('ele.s')
        value = ele%s
      case ('ele.l')
        value = ele%value(l$)
      case default
        nl=incr(nl); li(nl) = 'INVALID'
        call out_io (s_error$, r_name, 'python lat_ele1 {ele}|{which} {who}: Bad {who}: ' // who)
        return
      end select

      if (use_real_array_buffer) then
        n = n + 1
        if (n > size(re_array)) call re_allocate(re_array, 2*n)
        re_array(n) = value

      else
        if (who == 'orbit.state') then
          if (i == 1) then
            nl=incr(nl); write (li(nl), '(i0)') nint(value)
          else
            write (li(nl), '(2a, i0)') trim(li(nl)), ';', nint(value)
          endif
        else
          if (i == 1) then
            nl=incr(nl); write (li(nl), '(es24.16)') value
          else
            write (li(nl), '(2a, es24.16)') trim(li(nl)), ';', value
          endif
        endif
      endif
    enddo

  enddo

  if (use_real_array_buffer) then
    if (who == 'orbit.state') then
      if (.not. allocated(tao_c_interface_com%c_integer)) allocate (tao_c_interface_com%c_integer(n))
      if (size(tao_c_interface_com%c_integer) < n) then
        deallocate (tao_c_interface_com%c_integer)
        allocate (tao_c_interface_com%c_integer(n)) 
      endif

      tao_c_interface_com%n_int = n
      tao_c_interface_com%c_integer(1:n) = nint(re_array(1:n))

    else
      if (.not. allocated(tao_c_interface_com%c_real)) allocate (tao_c_interface_com%c_real(n))
      if (size(tao_c_interface_com%c_real) < n) then
        deallocate (tao_c_interface_com%c_real)
        allocate (tao_c_interface_com%c_real(n)) 
      endif

      tao_c_interface_com%n_real = n
      tao_c_interface_com%c_real(1:n) = re_array(1:n)
    endif
  endif

!----------------------------------------------------------------------
! Units of a parameter associated with a lattice or lattice element.
! Command syntax:
!   python lat_param_units {param_name}

case ('lat_param_units')

  name = upcase(line)
  a_name = attribute_units(name)
  nl=incr(nl); write(li(nl), '(a)') a_name

!----------------------------------------------------------------------
! Twiss at given s position.
! Command syntax:
!   python orbit_at_s {ix_uni}@{ix_branch}>>{s}|{which}
! where:
!   {which} is one of:
!     model
!     base
!     design
!   {s} is the longitudinal s-position.
case ('orbit_at_s')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, err); if (err) return
  ix_branch = parse_branch(line, .true., err); if (err) return
  s_pos = parse_real(line, err); if (err) return

  call twiss_and_track_at_s (tao_lat%lat, s_pos, orb = tao_lat%tao_branch(ix_branch)%orbit, orb_at_s = orb, ix_branch = ix_branch)
  call orbit_out (orb)

!----------------------------------------------------------------------
! List of plot templates or plot regions.
! Command syntax:  
!   python plot_list {r/g}
! where "{r/g}" is:
!   "r"      ! list regions
!   "t"      ! list template plots 


case ('plot_list')
  if (line == 't') then
    do i = 1, size(s%plot_page%template)
      p => s%plot_page%template(i)
      if (p%phantom) cycle
      if (p%name == '') cycle
      if (p%name == 'scratch') cycle
      nl=incr(nl); write (li(nl), '(i0, 2a)') i, ';', trim(p%name)
    enddo

  elseif (line == 'r') then
    do i = 1, size(s%plot_page%region)
      pr => s%plot_page%region(i)
      if (pr%name == '') cycle
      nl=incr(nl); write (li(nl), '(i0, 5a, l1)') i, ';', trim(pr%name), ';', trim(pr%plot%name), ';', pr%visible
    enddo

  else
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Expect "r" or "t"')
  endif

!----------------------------------------------------------------------
! Graph
! Syntax:
!   python plot_graph {graph_name}
! {graph_name} is in the form:
!   {p_name}.{g_name}
! where 
!   {p_name} is the plot region name if from a region or the plot name if a template plot.
!   This name is obtained from the python plot_list command. 
!   {g_name} is the graph name obtained from the python plot1 command.

case ('plot_graph')

  call tao_find_plots (err, line, 'COMPLETE', graph = graph)

  if (err .or. .not. allocated(graph)) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Bad graph name')
    call end_stuff()
    return
  endif

  g => graph(1)%g

  n = 0
  if (allocated(g%curve)) n = size(g%curve)

  nl=incr(nl); write (li(nl), imt) 'num_curves;INT;T;',                       n
  do i = 1, n
    nl=incr(nl); write (li(nl), vamt) 'curve[', i, '];STR;T;',                g%curve(i)%name
  enddo

  nl=incr(nl); write (li(nl), amt) 'name;STR;T;',                             g%name
  nl=incr(nl); write (li(nl), amt) 'type;STR;T;',                             g%type
  nl=incr(nl); write (li(nl), amt) 'title;STR;T;',                            g%title
  nl=incr(nl); write (li(nl), amt) 'title_suffix;STR;F;',                     g%title_suffix
  nl=incr(nl); write (li(nl), amt) 'component;STR;T;',                        g%component
  nl=incr(nl); write (li(nl), amt) 'why_invalid;STR;F;',                      g%why_invalid
  nl=incr(nl); write (li(nl), amt) 'floor_plan_view;STR;T;',                  g%floor_plan_view
  nl=incr(nl); write (li(nl), amt) 'floor_plan_orbit_color;STR;T;',           g%floor_plan_orbit_color
  nl=incr(nl); write (li(nl), rmt) 'x_axis_scale_factor;REAL;T;',             g%x_axis_scale_factor
  nl=incr(nl); write (li(nl), rmt) 'symbol_size_scale;REAL;T;',               g%symbol_size_scale
  nl=incr(nl); write (li(nl), rmt) 'floor_plan_rotation;REAL;T;',             g%floor_plan_rotation
  nl=incr(nl); write (li(nl), lmt) 'floor_plan_flip_axis_side;LOGIC;T;',      g%floor_plan_flip_label_side
  nl=incr(nl); write (li(nl), rmt) 'floor_plan_orbit_scale;REAL;T;',          g%floor_plan_orbit_scale
  nl=incr(nl); write (li(nl), imt) 'ix_branch;INT;T;',                        g%ix_branch
  nl=incr(nl); write (li(nl), imt) 'ix_universe;INT;T;',                      g%ix_universe
  nl=incr(nl); write (li(nl), lmt) 'clip;LOGIC;T;',                           g%clip
  nl=incr(nl); write (li(nl), lmt) 'valid;LOGIC;F;',                          g%valid
  nl=incr(nl); write (li(nl), lmt) 'y2_mirrors_y;LOGIC;T;',                   g%y2_mirrors_y
  nl=incr(nl); write (li(nl), lmt) 'limited;LOGIC;F;',                        g%limited
  nl=incr(nl); write (li(nl), lmt) 'draw_axes;LOGIC;T;',                      g%draw_axes
  nl=incr(nl); write (li(nl), lmt) 'correct_xy_distortion;LOGIC;T;',          g%correct_xy_distortion
  nl=incr(nl); write (li(nl), lmt) 'floor_plan_size_is_absolute;LOGIC;T;',    g%floor_plan_size_is_absolute
  nl=incr(nl); write (li(nl), lmt) 'floor_plan_draw_only_first_pass;LOGIC;T;',  g%floor_plan_draw_only_first_pass
  nl=incr(nl); write (li(nl), lmt) 'draw_curve_legend;LOGIC;T;',              g%draw_curve_legend
  nl=incr(nl); write (li(nl), lmt) 'draw_grid;LOGIC;T;',                      g%draw_grid
  nl=incr(nl); write (li(nl), lmt) 'draw_only_good_user_data_or_vars;LOGIC;T;', g%draw_only_good_user_data_or_vars

  nl=incr(nl); write (li(nl), amt) 'x.label;STR;T;',                         g%x%label
  nl=incr(nl); write (li(nl), rmt) 'x.max;REAL;T;',                          g%x%max
  nl=incr(nl); write (li(nl), rmt) 'x.min;REAL;T;',                          g%x%min
  nl=incr(nl); write (li(nl), imt) 'x.major_div;INT;T;',                     g%x%major_div
  nl=incr(nl); write (li(nl), imt) 'x.major_div_nominal;INT;T;',             g%x%major_div_nominal
  nl=incr(nl); write (li(nl), imt) 'x.places;INT;T;',                        g%x%places
  nl=incr(nl); write (li(nl), lmt) 'x.draw_label;LOGIC;T;',                  g%x%draw_label
  nl=incr(nl); write (li(nl), lmt) 'x.draw_numbers;LOGIC;T;',                g%x%draw_numbers

  nl=incr(nl); write (li(nl), amt) 'y.label;STR;T;',                         g%y%label
  nl=incr(nl); write (li(nl), rmt) 'y.max;REAL;T;',                          g%y%max
  nl=incr(nl); write (li(nl), rmt) 'y.min;REAL;T;',                          g%y%min
  nl=incr(nl); write (li(nl), imt) 'y.major_div;INT;T;',                     g%y%major_div
  nl=incr(nl); write (li(nl), imt) 'y.major_div_nominal;INT;T;',             g%y%major_div_nominal
  nl=incr(nl); write (li(nl), imt) 'y.places;INT;T;',                        g%y%places
  nl=incr(nl); write (li(nl), lmt) 'y.draw_label;LOGIC;T;',                  g%y%draw_label
  nl=incr(nl); write (li(nl), lmt) 'y.draw_numbers;LOGIC;T;',                g%y%draw_numbers

  nl=incr(nl); write (li(nl), amt) 'y2.label;STR;T;',                        g%y2%label
  nl=incr(nl); write (li(nl), rmt) 'y2.max;REAL;T;',                         g%y2%max
  nl=incr(nl); write (li(nl), rmt) 'y2.min;REAL;T;',                         g%y2%min
  nl=incr(nl); write (li(nl), imt) 'y2.major_div;INT;T;',                    g%y2%major_div
  nl=incr(nl); write (li(nl), imt) 'y2.major_div_nominal;INT;T;',            g%y2%major_div_nominal
  nl=incr(nl); write (li(nl), imt) 'y2.places;INT;T;',                       g%y2%places
  nl=incr(nl); write (li(nl), lmt) 'y2.draw_label;LOGIC;T;',                 g%y2%draw_label
  nl=incr(nl); write (li(nl), lmt) 'y2.draw_numbers;LOGIC;T;',               g%y2%draw_numbers

!----------------------------------------------------------------------
! Curve information for a plot
! Command syntax:
!   pyton plot_curve {curve_name}

case ('plot_curve')

  call tao_find_plots (err, line, 'COMPLETE', curve = curve)

  if (err .or. .not. allocated(curve)) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid curve')
    call end_stuff()
    return
  endif

  cur => curve(1)%c

  nl=incr(nl); write (li(nl), amt) 'name;STR;T;',                             cur%name
  nl=incr(nl); write (li(nl), amt) 'data_source;STR;T;',                      cur%data_source
  nl=incr(nl); write (li(nl), amt) 'data_type_x;STR;T;',                      cur%data_type_x
  nl=incr(nl); write (li(nl), amt) 'data_type_z;STR;T;',                      cur%data_type_z
  nl=incr(nl); write (li(nl), amt) 'data_type;STR;T;',                        cur%data_type
  nl=incr(nl); write (li(nl), amt) 'component;STR;T;',                        cur%component
  nl=incr(nl); write (li(nl), amt) 'ele_ref_name;STR;T;',                     cur%ele_ref_name
  nl=incr(nl); write (li(nl), amt) 'legend_text;STR;T;',                      cur%legend_text
  nl=incr(nl); write (li(nl), amt) 'message_text;STR;T;',                     cur%message_text
  nl=incr(nl); write (li(nl), amt) 'units;STR;T;',                            cur%units
  nl=incr(nl); write (li(nl), rmt) 'y_axis_scale_factor;REAL;T;',             cur%y_axis_scale_factor
  nl=incr(nl); write (li(nl), rmt) 's;REAL;F;',                               cur%s
  nl=incr(nl); write (li(nl), rmt) 'z_color0;REAL;T;',                        cur%z_color0
  nl=incr(nl); write (li(nl), rmt) 'z_color1;REAL;T;',                        cur%z_color1
  nl=incr(nl); write (li(nl), imt) 'ix_universe;INT;T;',                      cur%ix_universe
  nl=incr(nl); write (li(nl), imt) 'symbol_every;INT;T;',                     cur%symbol_every
  nl=incr(nl); write (li(nl), imt) 'ix_branch;INT;T;',                        cur%ix_branch
  nl=incr(nl); write (li(nl), imt) 'ix_ele_ref;INT;T;',                       cur%ix_ele_ref
  nl=incr(nl); write (li(nl), imt) 'ix_ele_ref_track;INT;T;',                 cur%ix_ele_ref_track
  nl=incr(nl); write (li(nl), imt) 'ix_bunch;INT;T;',                         cur%ix_bunch
  nl=incr(nl); write (li(nl), lmt) 'use_y2;LOGIC;T;',                         cur%use_y2
  nl=incr(nl); write (li(nl), lmt) 'draw_line;LOGIC;T;',                      cur%draw_line
  nl=incr(nl); write (li(nl), lmt) 'draw_symbols;LOGIC;T;',                   cur%draw_symbols
  nl=incr(nl); write (li(nl), lmt) 'draw_symbol_index;LOGIC;T;',              cur%draw_symbol_index
  nl=incr(nl); write (li(nl), lmt) 'smooth_line_calc;LOGIC;T;',               cur%smooth_line_calc
  nl=incr(nl); write (li(nl), lmt) 'use_z_color;LOGIC;T;',                    cur%use_z_color
  nl=incr(nl); write (li(nl), lmt) 'autoscale_z_color;LOGIC;T;',              cur%autoscale_z_color

  nl=incr(nl); write (li(nl), imt)  'line.width;INT;T;',                      cur%line%width
  nl=incr(nl); write (li(nl), amt)  'line.color;STR;T;',                      qp_color_name(cur%line%color)
  nl=incr(nl); write (li(nl), amt)  'line.pattern;STR;T;',                    qp_line_pattern_name(cur%line%pattern)

  nl=incr(nl); write (li(nl), amt)  'symbol.type;STR;T;',                     qp_symbol_type_name(cur%symbol%type)
  nl=incr(nl); write (li(nl), amt)  'symbol.color;STR;T;',                    qp_color_name(cur%symbol%color)
  nl=incr(nl); write (li(nl), rmt)  'symbol.height;REAL;T;',                  cur%symbol%height
  nl=incr(nl); write (li(nl), amt)  'symbol.fill_pattern;STR;T;',             qp_fill_name(cur%symbol%fill_pattern)
  nl=incr(nl); write (li(nl), imt)  'symbol.line_width;INT;T;',               cur%symbol%line_width

!----------------------------------------------------------------------
! Points used to construct a smooth line for a plot curve.
! Command syntax:
!   python plot_line {region_name}.{graph_name}.{curve_name} {x-or-y}
! Optional {x-or-y} may be set to "x" or "y" to get the smooth line points x or y component put into the real array buffer.
! Note: The plot must come from a region, and not a template, since no template plots have associated line data.
! Examples:
!   python plot_line r13.g.a       ! String array output.
!   python plot_line r13.g.a x     ! x-component of line points loaded into the real array buffer.
!   python plot_line r13.g.a y     ! y-component of line points loaded into the real array buffer.

case ('plot_line')

  call string_trim(line(ix_line+1:), who, ix2)
  line = line(1:ix_line)
  call tao_find_plots (err, line, 'COMPLETE', curve = curve)

  if (.not. allocated(curve) .or. size(curve) /= 1) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid curve')
    call end_stuff()
    return
  endif

  cur => curve(1)%c
  if (.not. allocated(cur%x_line)) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": No line associated with curve')
    call end_stuff()
    return
  endif
      
  n = size(cur%x_line)

  select case (who)
  case ('x', 'y')
    if (.not. allocated(tao_c_interface_com%c_real)) allocate (tao_c_interface_com%c_real(n))
    if (size(tao_c_interface_com%c_real) < n) then
      deallocate (tao_c_interface_com%c_real)
      allocate (tao_c_interface_com%c_real(n)) 
    endif

    tao_c_interface_com%n_real = n

    if (who == 'x') then
      tao_c_interface_com%c_real(1:n) = cur%x_line
    else
      tao_c_interface_com%c_real(1:n) = cur%y_line
    endif

  case ('')
    call re_allocate_lines (nl+n+100)
    do i = 1, n
      nl=incr(nl); write (li(nl), '(i0, 2(a, es24.16))') i, ';', cur%x_line(i), ';', cur%y_line(i)
    enddo

  case default
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": word after curve name not "x" nor "y"')
    call end_stuff()
  end select


!----------------------------------------------------------------------
! Locations to draw symbols for a plot curve.
! Command syntax:
!   python plot_symbol {region_name}.{graph_name}.{curve_name} {x-or-y}
! Optional {x-or-y} may be set to "x" or "y" to get the symbol x or y positions put into the real array buffer.
! Note: The plot must come from a region, and not a template, since no template plots have associated symbol data.
! Examples:
!   python plot_symbol r13.g.a       ! String array output.
!   python plot_symbol r13.g.a x     ! x-component of the symbol positions loaded into the real array buffer.
!   python plot_symbol r13.g.a y     ! y-component of the symbol positions loaded into the real array buffer.

case ('plot_symbol')

  call string_trim(line(ix_line+1:), who, ix2)
  line = line(1:ix_line)
  call tao_find_plots (err, line, 'COMPLETE', curve = curve)

  if (.not. allocated(curve) .or. size(curve) /= 1) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid curve')
    call end_stuff()
    return
  endif

  cur => curve(1)%c
  if (.not. allocated(cur%x_symb)) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": No line associated with curve')
    call end_stuff()
    return
  endif

  n = size(cur%x_symb)

  select case (who)
  case ('x', 'y')
    if (.not. allocated(tao_c_interface_com%c_real)) allocate (tao_c_interface_com%c_real(n))
    if (size(tao_c_interface_com%c_real) < n) then
      deallocate (tao_c_interface_com%c_real)
      allocate (tao_c_interface_com%c_real(n)) 
    endif

    tao_c_interface_com%n_real = n

    if (who == 'x') then
      tao_c_interface_com%c_real(1:n) = cur%x_symb
    else
      tao_c_interface_com%c_real(1:n) = cur%y_symb
    endif

  case ('')
    call re_allocate_lines (size(cur%x_symb)+100)
    do i = 1, size(cur%x_symb)
      nl=incr(nl); write (li(nl), '(2(i0, a), 2(es24.16, a))') i, ';', cur%ix_symb(i), ';', cur%x_symb(i), ';', cur%y_symb(i)
    enddo

  case default
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": word after curve name not "x" nor "y"')
    call end_stuff()
  end select

!----------------------------------------------------------------------
! Info on a given plot.
! Command syntax:
!   python plot1 {name}
! {name} should be the region name if the plot is associated with a region.
! Output syntax is variable list form. See documentation at beginning of this file.

case ('plot1')

  call tao_find_plots (err, line, 'COMPLETE', plot, print_flag = .false.)
  if (err) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Expect "r" or "t" at end.')
    call end_stuff()
    return
  endif

  p => plot(1)%p

  n = 0
  if (allocated(p%graph)) n = size(p%graph)

  nl=incr(nl); write (li(nl), imt) 'num_graphs;INT;T;',                       n
  do i = 1, n
    nl=incr(nl); write (li(nl), vamt) 'graph[', i, '];STR;T;',              p%graph(i)%name
  enddo

  nl=incr(nl); write (li(nl), amt) 'name;STR;T;',                             p%name
  nl=incr(nl); write (li(nl), amt) 'description;STR;T;',                      p%description
  nl=incr(nl); write (li(nl), amt) 'x_axis_type;STR;T;',                      p%x_axis_type
  nl=incr(nl); write (li(nl), lmt) 'autoscale_x;LOGIC;T;',                    p%autoscale_x
  nl=incr(nl); write (li(nl), lmt) 'autoscale_y;LOGIC;T;',                    p%autoscale_y
  nl=incr(nl); write (li(nl), lmt) 'autoscale_gang_x;LOGIC;T;',               p%autoscale_gang_x
  nl=incr(nl); write (li(nl), lmt) 'autoscale_gang_y;LOGIC;T;',               p%autoscale_gang_y
  nl=incr(nl); write (li(nl), lmt) 'list_with_show_plot_command;LOGIC;T;',    p%list_with_show_plot_command
  nl=incr(nl); write (li(nl), lmt) 'phantom;LOGIC;T;',                        p%phantom
  nl=incr(nl); write (li(nl), imt) 'n_curve_pts;INT;T;',                      p%n_curve_pts

!----------------------------------------------------------------------
! Convert species name to corresponding integer
! Command syntax:
!   python species_to_int {species_str}
! Example:
!   python species_to_int CO2++

case ('species_to_int')

  n = species_id(line)
  if (n == invalid$ .or. line == '') then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid species name.')
    call end_stuff()
    return
  endif

  nl=incr(nl); write (li(nl), '(i0)') n

!----------------------------------------------------------------------
! Convert species integer id to corresponding 
! Command syntax:
!   python species_to_str {species_int}
! Example:
!   python species_to_str -1     ! Returns 'Electron'

case ('species_to_str')

  n = string_to_int (line, 0, err)
  name = species_name(n)

  if (err .or. line == '' .or. name == invalid_name) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid species integer id number.')
    call end_stuff()
    return
  endif

  nl=incr(nl); write (li(nl), '(a)') trim(name)

!----------------------------------------------------------------------
! Twiss at given s position
! Command syntax:
!   python twiss_at_s {ix_uni}@{ix_branch}>>{s}|{which}
! where {which} is one of:
!   model
!   base
!   design

case ('twiss_at_s')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, err); if (err) return
  ix_branch = parse_branch(line, .true., err); if (err) return
  s_pos = parse_real(line, err); if (err) return

  call twiss_and_track_at_s (tao_lat%lat, s_pos, this_ele, tao_lat%tao_branch(ix_branch)%orbit, ix_branch = ix_branch)
  call twiss_out (this_ele%a, 'a')
  call twiss_out (this_ele%b, 'b')

!----------------------------------------------------------------------
! Universe info
! Command syntax:
!   python universe {ix_universe}
! Use "python global" to get the number of universes.

case ('universe')

  u => point_to_uni(line, .false., err); if (err) return
  
  nl=incr(nl); write (li(nl), imt) 'ix_uni;INT;F;',                           u%ix_uni
  nl=incr(nl); write (li(nl), imt) 'n_d2_data_used;INT;F;',                   u%n_d2_data_used
  nl=incr(nl); write (li(nl), imt) 'n_data_used;INT;F;',                      u%n_data_used
  nl=incr(nl); write (li(nl), lmt) 'reverse_tracking;LOGIC;T;',               u%reverse_tracking
  nl=incr(nl); write (li(nl), lmt) 'is_on;LOGIC;T;',                          u%is_on

!----------------------------------------------------------------------
! Create a v1 variable structure along with associated var array.
! Command syntax:
!   python var_create {v1_name} {n_var_min} {n_var_max}
! {n_var_min} and {n_var_max} are the lower and upper bounds of the var
! Example:
!   python var_create quad_k1 0 45
! This example creates a v1 var structure called "quad_k1" with an associated
! variable array that has the range [0, 45].
!
! Use the "set variable" command to set a created variable parameters.
! In particular, to slave a lattice parameter to a variable use the command:
!   set {v1_name}|ele_name = {lat_param}
! where {lat_param} is of the form {ix_uni}@{ele_name_or_location}{param_name}]
! Examples:
!   set quad_k1[2]|ele_name = 2@q01w[k1]
!   set quad_k1[2]|ele_name = 2@0>>10[k1]
! Note: When setting multiple variable parameters, temporarily toggle s%global%lattice_calc_on to False
!   ("set global lattice_calc_on = F") to prevent Tao trying to evaluate the partially created variable
!   and generating unwanted error messages.

case ('var_create')

  call tao_cmd_split (line, 3, name1, .true., err)

  if (err .or. .not. is_integer(name1(2)) .or. .not. is_integer(name1(3))) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Is Malformed')
    call end_stuff()
    return
  endif

  read (name1(2), *) ix_min(1)
  read (name1(3), *) ix_max(1)

  n1 = size(s%v1_var)
  if (s%n_v1_var_used + 1 > n1) then
    call move_alloc (s%v1_var, v1_temp)
    allocate (s%v1_var(s%n_v1_var_used + 1))
    s%v1_var(1:n1) = v1_temp
  endif

  n = size(s%var)
  n_delta = ix_max(1) + 1 - ix_min(1)
  if (s%n_var_used + n_delta  > n) then
    call move_alloc (s%var, v_temp)
    allocate (s%var(s%n_var_used+n_delta))
    s%var(1:n) = v_temp
    do k = s%n_var_used+1, size(s%var)
      s%var(k)%ix_var = k
    enddo
  endif

  i2 = 0   ! In case there are no v1 structures yet defined.

  do i = 1, s%n_v1_var_used
    v1_ptr => s%v1_var(i)
    i1 = lbound(v1_ptr%v, 1)
    i1 = v1_ptr%v(i1)%ix_var
    i2 = ubound(v1_ptr%v, 1)
    i2 = v1_ptr%v(i2)%ix_var
    call tao_point_v1_to_var (v1_ptr, s%var(i1:i2), s%var(i1)%ix_v1)
  enddo

  nn = s%n_v1_var_used + 1
  s%n_v1_var_used = nn
  s%v1_var(nn)%ix_v1_var = nn
  s%v1_var(nn)%name = name1(1)
  s%n_var_used = s%n_var_used + n_delta
  i1 = i2 + 1
  i2 = i2 + n_delta
  call tao_point_v1_to_var (s%v1_var(nn), s%var(i1:i2), ix_min(1))

!----------------------------------------------------------------------
! Destroy a v1 var structure along with associated var sub-array.
! Command syntax:
!   python var_destroy {v1_datum}

case ('var_destroy')

  call tao_find_var (err, line, v1_array = v1_array)
  if (err .or. .not. allocated(v1_array)) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid v1 var name')
    call end_stuff()
    return
  endif

  v1_ptr => v1_array(1)%v1
  i1 = lbound(v1_ptr%v, 1)
  i1 = v1_ptr%v(i1)%ix_var
  i2 = ubound(v1_ptr%v, 1)
  i2 = v1_ptr%v(i2)%ix_var

  n_delta = i2 + 1 - i1

  n = s%n_var_used

  do j = v1_ptr%ix_v1_var, s%n_v1_var_used - 1
    s%v1_var(j) = s%v1_var(j+1)
    v1_ptr => s%v1_var(j)
    i1 = v1_ptr%v(lbound(v1_ptr%v,1))%ix_var
    i2 = v1_ptr%v(ubound(v1_ptr%v,1))%ix_var
    s%var(i1-n_delta:i2-n_delta) = s%var(i1:i2)
    call tao_point_v1_to_var(v1_ptr, s%var(i1-n_delta:i2-n_delta), s%var(i1-n_delta)%ix_v1)
    do k = i1, i2
      s%var(i1-n_delta)%ix_var = i1 - n_delta
    enddo
  enddo

  s%n_v1_var_used = s%n_v1_var_used - 1
  s%n_var_used = s%n_var_used - n_delta

!----------------------------------------------------------------------
! List of all variable v1 arrays
! Command syntax: 
!   python var_general
! Output syntax:
!   {v1_var name};{v1_var%v lower bound};{v1_var%v upper bound}

case ('var_general')

  do i = 1, s%n_v1_var_used
    v1_ptr => s%v1_var(i)
    if (v1_ptr%name == '') cycle
    nl=incr(nl); write (li(nl), '(2a, 2(i0, a))') trim(v1_ptr%name), ';', lbound(v1_ptr%v, 1), ';', ubound(v1_ptr%v, 1)
  enddo

!----------------------------------------------------------------------
! List of variables in a given variable v1 array
! Command syntax: 
!   python var_v1 {v1_var}

case ('var_v1')

  call tao_find_var (err, line, v1_array = v1_array)

  if (err .or. .not. allocated(v1_array)) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid v1 name')
    call end_stuff()
    return
  endif

  v1_ptr => v1_array(1)%v1

  do i = lbound(v1_ptr%v, 1), ubound(v1_ptr%v, 1)
    v_ptr => v1_ptr%v(i)
    if (.not. v_ptr%exists) cycle
    nl=incr(nl); write (li(nl), '(2a, i0, 5a, 3(es24.16, a), 2 (l1, a))') trim(v1_ptr%name), '[', &
                     v_ptr%ix_v1, '];', trim(v_ptr%ele_name), ';', trim(v_ptr%attrib_name), ';', &
                     v_ptr%meas_value, ';', v_ptr%model_value, ';', &
                     v_ptr%design_value, ';', v_ptr%good_user, ';', v_ptr%useit_opt
  enddo

  nl=incr(nl); write (li(nl), imt) 'ix_v1_var;INT;F;',                       v1_ptr%ix_v1_var

!----------------------------------------------------------------------
! Info on an individual variable
! Command syntax: 
!   python var1 {var}
! Output syntax is variable list form. See documentation at beginning of this file.

case ('var1')

  call tao_find_var (err, line, v_array = v_array)

  if (.not. allocated(v_array) .or. size(v_array) /= 1) then
    nl=incr(nl); li(nl) = 'INVALID'
    call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid variable name')
    call end_stuff()
    return
  endif

  v_ptr => v_array(1)%v

  nl=incr(nl); write (li(nl), rmt)  'model_value;REAL;T;',          v_ptr%model_value
  nl=incr(nl); write (li(nl), rmt)  'base_value;REAL;T;',           v_ptr%base_value

  nl=incr(nl); write (li(nl), amt) 'ele_name;STR;T;',                         v_ptr%ele_name
  nl=incr(nl); write (li(nl), amt) 'attrib_name;STR;T;',                      v_ptr%attrib_name
  nl=incr(nl); write (li(nl), imt) 'ix_v1;INT;F;',                            v_ptr%ix_v1
  nl=incr(nl); write (li(nl), imt) 'ix_var;INT;F;',                           v_ptr%ix_var
  nl=incr(nl); write (li(nl), imt) 'ix_dvar;INT;F;',                          v_ptr%ix_dvar
  nl=incr(nl); write (li(nl), imt) 'ix_attrib;INT;F;',                        v_ptr%ix_attrib
  nl=incr(nl); write (li(nl), imt) 'ix_key_table;INT;T;',                     v_ptr%ix_key_table
  nl=incr(nl); write (li(nl), rmt) 'design_value;REAL;F;',                    v_ptr%design_value
  nl=incr(nl); write (li(nl), rmt) 'scratch_value;REAL;F;',                   v_ptr%scratch_value
  nl=incr(nl); write (li(nl), rmt) 'old_value;REAL;F;',                       v_ptr%old_value
  nl=incr(nl); write (li(nl), rmt) 'meas_value;REAL;T;',                      v_ptr%meas_value
  nl=incr(nl); write (li(nl), rmt) 'ref_value;REAL;T;',                       v_ptr%ref_value
  nl=incr(nl); write (li(nl), rmt) 'correction_value;REAL;F;',                v_ptr%correction_value
  nl=incr(nl); write (li(nl), rmt) 'high_lim;REAL;T;',                        v_ptr%high_lim
  nl=incr(nl); write (li(nl), rmt) 'low_lim;REAL;T;',                         v_ptr%low_lim
  nl=incr(nl); write (li(nl), rmt) 'step;REAL;T;',                            v_ptr%step
  nl=incr(nl); write (li(nl), rmt) 'weight;REAL;T;',                          v_ptr%weight
  nl=incr(nl); write (li(nl), rmt) 'delta_merit;REAL;F;',                     v_ptr%delta_merit
  nl=incr(nl); write (li(nl), rmt) 'merit;REAL;F;',                           v_ptr%merit
  nl=incr(nl); write (li(nl), rmt) 'dmerit_dvar;REAL;F;',                     v_ptr%dMerit_dVar
  nl=incr(nl); write (li(nl), rmt) 'key_val0;REAL;F;',                        v_ptr%key_val0
  nl=incr(nl); write (li(nl), rmt) 'key_delta;REAL;T;',                       v_ptr%key_delta
  nl=incr(nl); write (li(nl), rmt) 's;REAL;F;',                               v_ptr%s
  nl=incr(nl); write (li(nl), amt) 'merit_type;STR;T;',                       v_ptr%merit_type
  nl=incr(nl); write (li(nl), lmt) 'exists;LOGIC;F;',                         v_ptr%exists
  nl=incr(nl); write (li(nl), lmt) 'good_var;LOGIC;F;',                       v_ptr%good_var
  nl=incr(nl); write (li(nl), lmt) 'good_user;LOGIC;T;',                      v_ptr%good_user
  nl=incr(nl); write (li(nl), lmt) 'good_opt;LOGIC;T;',                       v_ptr%good_opt
  nl=incr(nl); write (li(nl), lmt) 'good_plot;LOGIC;T;',                      v_ptr%good_plot
  nl=incr(nl); write (li(nl), lmt) 'useit_opt;LOGIC;F;',                      v_ptr%useit_opt
  nl=incr(nl); write (li(nl), lmt) 'useit_plot;LOGIC;F;',                     v_ptr%useit_plot
  nl=incr(nl); write (li(nl), lmt) 'key_bound;LOGIC;T;',                      v_ptr%key_bound

!----------------------------------------------------------------------

case default

  call out_io (s_error$, r_name, "python command internal error, shouldn't be here!")

end select

call end_stuff()

!----------------------------------------------------------------------
! return through scratch

contains

subroutine end_stuff()

integer i

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

function point_to_uni (line, compound_word, err) result (u)

type (tao_universe_struct), pointer :: u
integer ix, ix_universe
logical compound_word, err
character(*) line
character(40) str

! A compound_word is something like "2@q10w" or "q10w". A non-compound word is something like "2" which 
! just represents a universe index.

nullify(u)
err = .false.

if (compound_word) then
  ix = index(line, '@')
  if (ix == 0) then
    str = '1' ! Use universe 1 as a default.
  else
    str = line(1:ix-1)
    line = line(ix+1:)
  endif
else
  ! In this case line is just a universe number
  str = line
endif

read (str, *,  iostat = ios)  ix_universe
if (ios /= 0) ix_universe = -999

u => tao_pointer_to_universe(ix_universe)

if (.not. associated(u)) then
  nl=incr(nl); li(nl) = 'INVALID'
  call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": bad universe index')
  call end_stuff()
  err = .true.
endif

end function point_to_uni

!----------------------------------------------------------------------
! contains

function incr(n) result (n1)

integer n, n1

n1 = n + 1
if (n1 > size(li)) call re_allocate_lines(int(1.5 * n1))

end function

!----------------------------------------------------------------------
! contains

subroutine re_allocate_lines (n_lines)

integer n_lines

if (.not. allocated(li)) allocate (li(n_lines))
if (size(li) < n_lines) call re_allocate (li, n_lines)

end subroutine re_allocate_lines

!----------------------------------------------------------------------
! contains

function point_to_tao_lat (line, err, which, who) result (tao_lat)

type (tao_lattice_struct), pointer :: tao_lat
logical err
character(*) line
character(*), optional :: which, who


err = .false.
nullify(tao_lat)

call string_trim(line, line, ix)
if (present(who)) call string_trim(line(ix+1:), who, i)
line = line(1:ix)

ix = index(line, '|')
if (ix == 0) then
  nl=incr(nl); li(nl) = 'INVALID'
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
  nl=incr(nl); li(nl) = 'INVALID'
  call out_io (s_error$, r_name, 'python ' // trim(input_str) //  ': Expecting "|<{hich}" where {which} must be one of "model", "base", or "design"')
  err = .true.
end select

if (present(which)) which = line(ix+1:)
line = line(1:ix-1)

end function point_to_tao_lat

!----------------------------------------------------------------------
! contains

function point_to_ele (line, err) result (ele)

type (ele_struct), pointer :: ele
character(*) line
logical err

!

err = .true.
nullify(ele)
call lat_ele_locator (line, tao_lat%lat, eles, n_loc)

if (n_loc /= 1) then
  nl=incr(nl); li(nl) = 'INVALID'
  call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Cannot locate element.')
  return
endif

ele => eles(1)%ele
err = .false.

end function point_to_ele

!----------------------------------------------------------------------
! contains

function parse_branch (line, has_separator, err) result (ix_branch)

integer ix_branch
logical has_separator, err
character(*) line
character(40) str

!

err = .false.

if (has_separator) then
  ix = index(line, '>>')

  if (ix == 0) then
    nl=incr(nl); li(nl) = 'INVALID'
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
  nl=incr(nl); li(nl) = 'INVALID'
  call out_io (s_error$, r_name, '"python ' // trim(input_str) // '" missing or out of range branch index')
  call end_stuff()
  err = .true.
  return
endif

end function parse_branch

!----------------------------------------------------------------------
! contains

function parse_real (line, err) result (a_real)

real(rp) a_real
logical err
character(*) line

a_real = string_to_real (line, real_garbage$, err)
if (err .or. a_real == real_garbage$) then
  nl=incr(nl); li(nl) = 'INVALID'
  call out_io (s_error$, r_name, '"python ' // trim(input_str) // '" Bad real number')
  call end_stuff()
  return
endif

end function parse_real

!----------------------------------------------------------------------
! contains

subroutine orbit_out (orbit)

type (coord_struct) orbit

nl=incr(nl); write (li(nl), rmt) 'x;REAL;F;',                                orbit%vec(1)
nl=incr(nl); write (li(nl), rmt) 'px;REAL;F;',                               orbit%vec(2)
nl=incr(nl); write (li(nl), rmt) 'y;REAL;F;',                                orbit%vec(3)
nl=incr(nl); write (li(nl), rmt) 'py;REAL;F;',                               orbit%vec(4)
nl=incr(nl); write (li(nl), rmt) 'z;REAL;F;',                                orbit%vec(5)
nl=incr(nl); write (li(nl), rmt) 'pz;REAL;F;',                               orbit%vec(6)

nl=incr(nl); write (li(nl), rmt) 'spin_x;REAL;F;',                           orbit%spin(1)
nl=incr(nl); write (li(nl), rmt) 'spin_y;REAL;F;',                           orbit%spin(2)
nl=incr(nl); write (li(nl), rmt) 'spin_z;REAL;F;',                           orbit%spin(3)

nl=incr(nl); write (li(nl), rmt) 'field_x;REAL;F;',                          orbit%field(1)
nl=incr(nl); write (li(nl), rmt) 'field_y;REAL;F;',                          orbit%field(2)

nl=incr(nl); write (li(nl), rmt) 'phase_x;REAL;F;',                          orbit%phase(1)
nl=incr(nl); write (li(nl), rmt) 'phase_y;REAL;F;',                          orbit%phase(2)

nl=incr(nl); write (li(nl), rmt) 's;REAL;F;',                                orbit%s
nl=incr(nl); write (li(nl), rmt) 't;REAL;F;',                                orbit%t
nl=incr(nl); write (li(nl), rmt) 'charge;REAL;F;',                           orbit%charge
nl=incr(nl); write (li(nl), rmt) 'path_len;REAL;F;',                         orbit%path_len
nl=incr(nl); write (li(nl), rmt) 'p0c;REAL;F;',                              orbit%p0c
nl=incr(nl); write (li(nl), rmt) 'beta;REAL;F;',                             orbit%beta
nl=incr(nl); write (li(nl), imt) 'ix_ele;INT;F;',                            orbit%ix_ele
nl=incr(nl); write (li(nl), amt) 'state;STR;F;',                             coord_state_name(orbit%state)
nl=incr(nl); write (li(nl), imt) 'direction;INT;F;',                         orbit%direction
nl=incr(nl); write (li(nl), amt) 'species;STR;F;',                           species_name(orbit%species)
nl=incr(nl); write (li(nl), amt) 'location;STR;F;',                          location_name(orbit%location)

end subroutine orbit_out

!----------------------------------------------------------------------
! contains

subroutine coord_out(beam, coordinate)
type (beam_struct), pointer :: beam
type (bunch_struct), pointer :: bunch
character(20) coordinate
integer :: i_vec

if (.not. allocated(beam%bunch)) then
    print *, 'no beam'
    return
endif

! TODO: generalize to n bunches
bunch => beam%bunch(1)

! Allocate scratch 
n = size(bunch%particle)
if (.not. allocated(tao_c_interface_com%c_real)) allocate (tao_c_interface_com%c_real(n))
if (size(tao_c_interface_com%c_real) < n) then
  deallocate (tao_c_interface_com%c_real)
  allocate (tao_c_interface_com%c_real(n)) 
endif
! Set scratch size
tao_c_interface_com%n_real = n

! Add data
select case (coordinate)
case ('x')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%vec(1)
case ('px')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%vec(2)
case ('y')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%vec(3)
case ('py')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%vec(4)
case ('z')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%vec(5)
case ('pz')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%vec(6)      
case ('s')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%s
case ('t')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%t
case ('charge')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%charge
case ('p0c')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%p0c

case default
  nl=incr(nl); li(nl) = 'INVALID'
  call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": coordinate not "x", "px", etc. ')
  call end_stuff()
  return
end select


end subroutine



!----------------------------------------------------------------------
! contains

subroutine twiss_out (twiss, suffix, emit_out, can_vary)

type (twiss_struct) twiss
character(*) suffix
character(20) fmt
character(8) v_str
logical, optional :: emit_out, can_vary

if (logic_option(.false., can_vary)) then
  v_str = ';REAL;T;'
else
  v_str = ';REAL;F;'
endif

fmt = '(3a, es24.16)'

nl=incr(nl); write (li(nl), fmt) 'beta_', suffix, v_str,                          twiss%beta
nl=incr(nl); write (li(nl), fmt) 'alpha_', suffix, v_str,                         twiss%alpha
nl=incr(nl); write (li(nl), fmt) 'gamma_', suffix, ';REAL;F;',                         twiss%gamma
nl=incr(nl); write (li(nl), fmt) 'phi_', suffix, v_str,                           twiss%phi
nl=incr(nl); write (li(nl), fmt) 'eta_', suffix, v_str,                           twiss%eta
nl=incr(nl); write (li(nl), fmt) 'etap_', suffix, v_str,                          twiss%etap

if (logic_option(.false., emit_out)) then
  nl=incr(nl); write (li(nl), fmt) 'sigma_', suffix, ';REAL;F;',                         twiss%sigma
  nl=incr(nl); write (li(nl), fmt) 'sigma_p_', suffix, ';REAL;F;',                       twiss%sigma_p
  nl=incr(nl); write (li(nl), fmt) 'emit_', suffix, ';REAL;F;',                          twiss%emit
  nl=incr(nl); write (li(nl), fmt) 'norm_emit_', suffix, ';REAL;F;',                     twiss%norm_emit
endif

end subroutine twiss_out

!----------------------------------------------------------------------
! contains

subroutine destroy_this_data(d_name)

character(*) d_name

call tao_find_data (err, d_name, d2_array = d2_array)
if (err .or. .not. allocated(d2_array)) then
  nl=incr(nl); li(nl) = 'INVALID'
  call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": Not a valid d2 data name')
  call end_stuff()
  return
endif

d2_ptr => d2_array(1)%d2
u => s%u(d2_ptr%ix_uni)
ix_d2 = d2_ptr%ix_d2_data

d1_ptr => d2_ptr%d1(1)
i1 = lbound(d1_ptr%d, 1)
i1 = d1_ptr%d(i1)%ix_data

n1 = size(d2_ptr%d1)
d1_ptr => d2_ptr%d1(n1)
i2 = ubound(d1_ptr%d, 1)
i2 = d1_ptr%d(i2)%ix_data

n_delta = i2 + 1 - i1

! Squeeze u%d2_data and u%data arrays

do i = ix_d2, u%n_d2_data_used - 1
  u%d2_data(i) = u%d2_data(i+1)
  do j = 1, size(u%d2_data(i)%d1)
    d1_ptr => u%d2_data(i)%d1(j)
    d1_ptr%d2 => u%d2_data(i)
    i1 = d1_ptr%d(lbound(d1_ptr%d,1))%ix_data
    i2 = d1_ptr%d(ubound(d1_ptr%d,1))%ix_data
    u%data(i1-n_delta:i2-n_delta) = u%data(i1:i2)
    call tao_point_d1_to_data(d1_ptr, u%data(i1-n_delta:i2-n_delta), u%data(i1-n_delta)%ix_d1)
    do k = i1, i2
      u%data(i1-n_delta)%ix_data = i1 - n_delta
    enddo
  enddo
enddo

u%n_d2_data_used = u%n_d2_data_used - 1
u%n_data_used = u%n_data_used - n_delta

end subroutine destroy_this_data

!----------------------------------------------------------------------
! contains

function match_ele_name (match_str, ele, err) result (is_a_match)

type (ele_struct) ele

integer ix, key

character(*) match_str
character(60) string

logical err, is_a_match

!

is_a_match = .false.
err = .false.

if (match_str == '*') then
  is_a_match = .true.
  return
endif

! key::name construct

string = match_str
ix = index(string, '::')
if (ix /= 0) then
  if (string(:ix-1) /= "*") then
    key = key_name_to_key_index (string(:ix-1), .true.)
    if (key < 1) then
      nl=incr(nl); li(nl) = 'INVALID'
      call out_io (s_error$, r_name, '"python ' // trim(input_str) // '": BAD ELEMENT KEY: ' // string(:ix-1))
      call end_stuff()
      err = .true.
      return
    endif

    if (ele%key /= key) then
      is_a_match = .false.
      return
    endif
  endif

  string = string(ix+2:)
endif

! 

if (index(string, "*") /= 0 .or. index(string, "%") /= 0) then
  is_a_match = match_wild(ele%name, string)
else
  is_a_match = (ele%name == string)
endif

end function match_ele_name

end subroutine tao_python_cmd
