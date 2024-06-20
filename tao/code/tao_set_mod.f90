module tao_set_mod

use tao_interface
use tao_data_and_eval_mod

implicit none

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_set_tune_cmd (branch_str, mask_str, print_list, qa_str, qb_str, delta_input)
!
! Routine to set the transverse tunes.
!
! Input:
!   branch_str    -- character(*): List of branches to apply tune set to.
!   mask_str      -- character(*): List of quadrupoles to veto.
!   print_list    -- logical: If True, print a list of elements varied and coefficients.
!   qa_str        -- character(*): Expression for Qa tune.
!   qb_str        -- character(*): Expression for Qb tune.
!   delta_input   -- logical: If true then qa_str and qb_str are deltas from present tune.
!-

subroutine tao_set_tune_cmd (branch_str, mask_str, print_list, qa_str, qb_str, delta_input)

implicit none

type (branch_pointer_struct), allocatable, target :: branches(:)
type (branch_struct), pointer :: branch
type (tao_universe_pointer_struct), allocatable, target :: unis(:)
type (tao_universe_struct), pointer :: u
type (ele_pointer_struct), allocatable :: eles(:)
type (lat_ele_order_struct) order
type (nametable_struct) nametab

real(rp) qa_val, qb_val
real(rp), allocatable :: dk1(:)

integer n, i, j, k, n_match
integer, allocatable :: indx(:)

logical print_list, delta_input, err

character(40) ele_name
character(*) branch_str, mask_str, qa_str, qb_str
character(*), parameter :: r_name = 'tao_set_tune_cmd'

! Evaluate expressions

call tao_pointer_to_branches(branch_str, branches, unis, err); if (err) return

do i = 1, size(branches)
  u => unis(i)%u
  branch => branches(i)%branch

  n = branch%n_ele_track
  qa_val = tao_evaluate_tune (qa_str, branch%ele(n)%a%phi/twopi, delta_input); if (qa_val == 0) return
  qb_val = tao_evaluate_tune (qb_str, branch%ele(n)%b%phi/twopi, delta_input); if (qb_val == 0) return

  call choose_quads_for_set_tune(branch, dk1, eles, mask_str, err)
  if (err) then
    call out_io (s_error$, r_name, &
      'CANNOT FIND A QUAD WITH BETA_A < BETA_B AND A QUAD WITH BETA_A > BETA_B (BOTH WITH NO TILT).')
    return
  endif

  if (print_list .and. i == 1) then
    call ele_order_calc (branch%lat, order)
    allocate(indx(size(dk1)))
    call nametable_init(nametab)
    do j = 1, size(dk1)
      ele_name = eles(j)%ele%name
      n = nametable_bracket_indexx (nametab, ele_name, n_match)
      if (n_match == 0) then
        call nametable_add (nametab, ele_name, nametab%n_max+1)
        indx(nametab%n_max) = j
        call out_io (s_blank$, r_name, ele_name // real_str(dk1(j), 4))

      else
        k = nametab%index(n)
        if (dk1(j) == dk1(indx(k))) cycle
        call out_io (s_blank$, r_name, ele_unique_name(eles(j)%ele, order) // real_str(dk1(j), 4))
      endif
    enddo
  endif

  err = .not. set_tune(twopi*qa_val, twopi*qb_val, dk1, eles, branch, u%model%tao_branch(branch%ix_branch)%orbit)
  u%calc%lattice = .true.
enddo

end subroutine tao_set_tune_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_set_z_tune_cmd (branch_str, q_str, delta_input)
!
! Routine to set the z-tune.
!
! Input:
!   branch_str    -- character(*): List of branches to apply tune set to.
!   q_str         -- character(*): Expression for Qc tune.
!   delta_input   -- logical: If true then qa_str and qb_str are deltas from present tune.
!-

subroutine tao_set_z_tune_cmd (branch_str, q_str, delta_input)

implicit none

type (branch_pointer_struct), allocatable, target :: branches(:)
type (branch_struct), pointer :: branch
type (tao_universe_pointer_struct), allocatable, target :: unis(:)
type (tao_universe_struct), pointer :: u

real(rp) q_val 
real(rp), allocatable :: dk1(:)

integer n, i
logical delta_input, err, ok

character(*) branch_str, q_str
character(*), parameter :: r_name = 'tao_set_z_tune_cmd'

! Evaluate expressions

call tao_pointer_to_branches(branch_str, branches, unis, err); if (err) return

do i = 1, size(branches)
  u => unis(i)%u
  branch => branches(i)%branch

  call calc_z_tune(branch)
  q_val = tao_evaluate_tune (q_str, branch%z%tune/twopi, delta_input); if (q_val == 0) return

  call set_z_tune(branch, twopi*q_val, ok)
  u%calc%lattice = .true.
enddo

end subroutine tao_set_z_tune_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_set_calculate_cmd (switch)
!
! Toggles off lattice calc and plotting.
!-

subroutine tao_set_calculate_cmd (switch)

character(*), optional :: switch
character(8) what
character(*), parameter :: r_name = 'tao_set_calculate_cmd'

!

what = 'toggle'

if (present(switch)) then
  if ((index('off', trim(switch)) == 1 .or. index('on', trim(switch)) == 1) .and. len_trim(switch) > 1) then
    what = switch
  elseif (switch /= '') then
    call out_io (s_error$, r_name, 'BAD SWITCH: ' // switch)
    return
  endif
endif

select case (what)
case ('toggle')
  s%global%lattice_calc_on = (.not. s%global%lattice_calc_on)
  s%global%plot_on = s%global%lattice_calc_on 

case ('off')
  s%global%lattice_calc_on = .false.
  s%global%plot_on = .false.

case ('on')
  s%global%lattice_calc_on = .true.
  s%global%plot_on = .true.
end select

end subroutine tao_set_calculate_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_set_key_cmd (key_str, cmd_str)
!
! Associates a command with a key press for single mode.
!
! Input:
!   key_str   -- character(*): keyboard key.
!   cmd_str   -- character(*): Command associated with key.
!-

subroutine tao_set_key_cmd (key_str, cmd_str)

integer i, n
character(*) key_str, cmd_str
character(*), parameter :: r_name = 'tao_set_key_cmd'

!

do i = 1, size(s%com%key)
  if (s%com%key(i)%name /= '' .and. s%com%key(i)%name /= key_str) cycle

  if (cmd_str == 'default') then
    if (s%com%key(i)%name /= key_str) then
      call out_io (s_error$, r_name, 'Key has not been set to begin with. Nothing to do.')
      return
    endif
    n = size(s%com%key)
    s%com%key(i:n) = [s%com%key(i+1:n), tao_alias_struct()]
  else
    s%com%key(i) = tao_alias_struct(key_str, cmd_str)
  endif

  return
enddo

call out_io (s_error$, r_name, 'KEY TABLE ARRAY OVERFLOW! PLEASE GET HELP!')

end subroutine tao_set_key_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_set_ran_state_cmd (state_string)
!
! Sets the random number generator state.
!
! Input:
!   state_string -- Character(*): Encoded random number state.
!-

subroutine tao_set_ran_state_cmd (state_string)

implicit none

type (random_state_struct) ran_state
character(*) state_string
character(100) state_str
character(*), parameter :: r_name = 'tao_set_ran_state_cmd'
integer ix, ios

!

call ran_default_state (get_state = ran_state)

call string_trim(state_string, state_str, ix)
read (state_str, *, iostat = ios) ran_state%ix 
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT READ FIRST RAN_STATE COMPONENT')
  return
endif

call string_trim(state_str(ix+1:), state_str, ix)
read (state_str, *, iostat = ios) ran_state%iy 
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT READ SECOND RAN_STATE COMPONENT')
  return
endif

call string_trim(state_str(ix+1:), state_str, ix)
read (state_str, *, iostat = ios) ran_state%number_stored 
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT READ THIRD RAN_STATE COMPONENT')
  return
endif

call string_trim(state_str(ix+1:), state_str, ix)
read (state_str, *, iostat = ios) ran_state%h_saved
if (ios /= 0 .or. ix == 0) then
  call out_io (s_error$, r_name, 'CANNOT READ FOURTH RAN_STATE COMPONENT')
  return
endif

call ran_default_state (set_state = ran_state)

end subroutine tao_set_ran_state_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_set_lattice_cmd (dest_lat, source_lat)
!
! Sets a lattice equal to another. This will also update the data structs
!
! Input:
!   dest_lat -- Character(*): Maybe: 'model', 'design', or 'base' with 
!                     optional '@n' at beginning to indicate the universe
!   source_lat  -- Character(*): Maybe: 'model', 'design', or 'base' 
!
!  Output:
!    s%u(n) -- lat_struct: changes specified lattice in specified universe 
!-

subroutine tao_set_lattice_cmd (dest_lat, source_lat)

implicit none

character(*) dest_lat, source_lat
character(16) dest1_name
character(*), parameter :: r_name = 'tao_set_lattice_cmd'

real(rp) source_val

integer i, j

logical, allocatable :: this_u(:)
logical err

! Lattice transfer

call tao_pick_universe (dest_lat, dest1_name, this_u, err)
if (err) return

do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. this_u(i)) cycle
  call set_lat (s%u(i))
  if (err) return
enddo

call tao_var_repoint()

! Variable transfer for those variables which vary parameters of the affected universe(s).
! This only needs to be done when dest_lat is a model lattice.

if (dest1_name == 'model') then
  do i = 1, s%n_var_used

    do j = 1, size(s%var(i)%slave)
      if (.not. this_u(s%var(i)%slave(j)%ix_uni)) cycle

      select case (source_lat)
      case ('model')
        source_val = s%var(i)%slave(j)%model_value
      case ('base')
        source_val = s%var(i)%slave(j)%base_value
      case ('design')
        source_val = s%var(i)%design_value
      end select

      call tao_set_var_model_value (s%var(i), source_val)
      exit
    enddo

  enddo
endif

!-------------------------------------------
contains

subroutine set_lat (u)

implicit none

type (tao_universe_struct), target :: u
type (tao_lattice_struct), pointer :: dest1_lat
type (tao_lattice_struct), pointer :: source1_lat
real(rp), pointer :: dest_data(:), source_data(:)
logical, pointer :: dest_good(:), source_good(:)
logical calc_ok

integer j, ib

!

err = .false.

select case (dest1_name)
case ('model')
  u%calc%lattice = .true.
  dest1_lat => u%model
  dest_data => u%data%model_value
  dest_good => u%data%good_model
case ('base')
  dest1_lat => u%base
  dest_data => u%data%base_value
  dest_good => u%data%good_base
case default
  call out_io (s_error$, r_name, 'BAD NAME: ' // dest_lat)
  err = .true.
  return
end select

select case (source_lat)
case ('model')
  ! make sure model data is up to date
  call tao_lattice_calc (calc_ok)
  source1_lat => u%model
  source_data => u%data%model_value
  source_good => u%data%good_model
case ('base')
  source1_lat => u%base
  source_data => u%data%base_value
  source_good => u%data%good_base
case ('design')
  source1_lat => u%design
  source_data => u%data%design_value
  source_good => u%data%good_design
case default
  call out_io (s_error$, r_name, 'BAD NAME: ' // source_lat)
  err = .true.
  return
end select

! dest_lat = source_lat will not mess up the pointers in s%var since both lattices have the same
! number of elements and therefore no reallocation needs to be done.

dest1_lat%lat = source1_lat%lat

do ib = 0, ubound(dest1_lat%tao_branch, 1)
  dest1_lat%tao_branch(ib)   = source1_lat%tao_branch(ib)
  do j = lbound(dest1_lat%tao_branch(ib)%bunch_params, 1), ubound(dest1_lat%tao_branch(ib)%bunch_params, 1)
    dest1_lat%tao_branch(ib)%bunch_params(j) = source1_lat%tao_branch(ib)%bunch_params(j)
  enddo
enddo

! Transfer the data

dest_data = source_data
dest_good = source_good

end subroutine set_lat

end subroutine tao_set_lattice_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_global_cmd (who, value_str)
!
! Routine to set global variables
! 
! Input:
!   who       -- Character(*): which global variable to set
!   value_str -- Character(*): Value to set to.
!
! Output:
!    s%global  -- Global variables structure.
!-

subroutine tao_set_global_cmd (who, value_str)

use bookkeeper_mod, only: set_on_off

implicit none

type (tao_global_struct) global, old_global
type (tao_universe_struct), pointer :: u

character(*) who, value_str
character(*), parameter :: r_name = 'tao_set_global_cmd'
character(len(value_str)+24) val

real(rp), allocatable :: set_val(:)
integer iu, ios, iuni, i, ix
logical err

namelist / params / global

! Special cases

if (who == 'phase_units') then
  call match_word (value_str, angle_units_name, ix)
  if (ix == 0) then
    call out_io (s_error$, r_name, 'BAD "phase_units" VALUE: ' // value_str)
  else
    s%global%phase_units = ix
  endif
  return
endif

if (who == 'silent_run') then  ! Old style syntax
  call out_io (s_error$, r_name, 'The "set global silent_run" command has been replaced by "set global quiet =  <logic>"')

  if (s%com%cmd_file_level == 0) then
    call out_io (s_error$, r_name, 'The "set global silent_run" command can only be used in command files.')
    return
  endif

  call match_word (value_str, [character(8):: 'true', 'false'], ix)
  if (ix == 1) then
    call tao_quiet_set('all')
  elseif (ix == 2) then
    call tao_quiet_set('off')
  else
    call out_io (s_error$, r_name, 'BAD "silent_run" VALUE: ' // value_str)
  endif
  return
endif

if (who == 'quiet') then
  if (s%com%cmd_file_level == 0) then
    call out_io (s_error$, r_name, 'The "set global quiet" command can only be used in command files.')
    return
  endif

  call tao_quiet_set (value_str)
  return
endif

! Surprisingly enough, a namelist read will ignore a blank value field so catch this problem here.

if (value_str == '') then
  call out_io (s_error$, r_name, 'SET VALUE IS BLANK!')
  return
endif

! open a scratch file for a namelist read

iu = tao_open_scratch_file (err);  if (err) return
ios = 0

select case (who)
case ('random_engine', 'random_gauss_converter', 'track_type', 'quiet', 'prompt_color'&
      'prompt_string', 'optimizer', 'print_command', 'var_out_file', 'history_file')
  val = quote(value_str)

case ('n_opti_cycles', 'n_opti_loops', 'phase_units', 'bunch_to_plot', &
      'random_seed', 'n_top10_merit', 'srdt_gen_n_slices', 'srdt_sxt_n_slices')
  call tao_evaluate_expression (value_str, 1, .false., set_val, err); if (err) return
  write (val, '(i0)', iostat = ios) nint(set_val(1))

case ('lm_opt_deriv_reinit', 'de_lm_step_ratio', 'de_var_to_population_factor', 'lmdif_eps', &
      'lmdif_negligible_merit', 'svd_cutoff', 'unstable_penalty', 'merit_stop_value', &
      'dmerit_stop_value', 'random_sigma_cutoff', 'delta_e_chrom')
  call tao_evaluate_expression (value_str, 1, .false., set_val, err); if (err) return
  write (val, '(es24.16)', iostat = ios) set_val(1)

case default
  val = value_str
end select

if (ios /= 0) then
  call out_io (s_error$, r_name, 'BAD NUMBER: ' // value_str)
  return
endif

write (iu, '(a)') '&params'
write (iu, '(a)') ' global%' // trim(who) // ' = ' // trim(val)
write (iu, '(a)') '/'
write (iu, *)
rewind (iu)
global = s%global  ! set defaults
read (iu, nml = params, iostat = ios)
close (iu, status = 'delete')

if (ios /= 0) then
  call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
  return
endif

call tao_data_check (err)
if (err) return

!

select case (who)
case ('optimizer')
  if (all(global%optimizer /= tao_optimizer_name)) then
    call out_io (s_error$, r_name, 'BAD OPTIMIZER NAME: ' // global%optimizer)
    return
  endif
case ('prompt_color')
  call upcase_string(global%prompt_color)
case ('random_seed')
  call ran_seed_put (global%random_seed)
case ('random_engine')
  call ran_engine (global%random_engine)
case ('random_gauss_converter', 'random_sigma_cutoff')
 call ran_gauss_converter (global%random_gauss_converter, global%random_sigma_cutoff)
case ('rf_on')
  do iuni = lbound(s%u, 1), ubound(s%u, 1)
    u => s%u(iuni)
    if (global%rf_on) then
      call set_on_off (rfcavity$, u%model%lat, on$)
    else
      call set_on_off (rfcavity$, u%model%lat, off$)
    endif
  enddo
  s%u%calc%lattice = .true.
case ('track_type')
  if (unquote(value_str) /= 'single' .and. unquote(value_str) /= 'beam') then
    call out_io (s_error$, r_name, 'BAD VALUE. MUST BE "single" OR "beam".')
    return
  endif
  s%u%calc%lattice = .true.
case ('srdt_gen_n_slices', 'srdt_sxt_n_slices', 'srdt_use_cache', 'init_lat_sigma_from_beam')
  s%u%calc%lattice = .true.
case ('symbol_import')
  if (global%symbol_import) then
    do iuni = lbound(s%u, 1), ubound(s%u, 1)
      call tao_symbol_import_from_lat (s%u(iuni)%model%lat)
    enddo
  endif
end select

!

s%global = global

!

select case (who)
case ('lattice_calc_on')
  if (s%global%lattice_calc_on) then
    do iuni = lbound(s%u, 1), ubound(s%u, 1)
      u => s%u(iuni)
      call lattice_bookkeeper(u%model%lat)
    enddo
  endif
end select

end subroutine tao_set_global_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_space_charge_com_cmd (who, value_str)
!
! Routine to set space_charge_com variables
! 
! Input:
!   who       -- Character(*): which space_charge_com variable to set
!   value_str -- Character(*): Value to set to.
!
! Output:
!    space_charge_com  -- space_charge_com variables structure.
!-

subroutine tao_set_space_charge_com_cmd (who, value_str)

implicit none

type (space_charge_common_struct) local_space_charge_com

character(*) who, value_str
character(len(value_str)+24) val
character(*), parameter :: r_name = 'tao_set_space_charge_com_cmd'

real(rp), allocatable :: set_val(:)
integer iu, ios
logical err

namelist / params / local_space_charge_com

! open a scratch file for a namelist read

select case (who)
case ('diagnostic_output_file')
  space_charge_com%diagnostic_output_file = unquote(value_str)
  return

case ('ds_track_step', 'beam_chamber_height', 'sigma_cutoff')
  call tao_evaluate_expression (value_str, 1, .false., set_val, err); if (err) return
  write (val, '(es24.16)', iostat = ios) set_val(1)

case ('n_bin', 'particle_bin_span', 'n_shield_images', 'sc_min_in_bin')
  call tao_evaluate_expression (value_str, 1, .false., set_val, err); if (err) return
  write (val, '(i0)', iostat = ios) nint(set_val(1))

case default  ! Is logical
  val = value_str
end select

iu = tao_open_scratch_file (err);  if (err) return

write (iu, '(a)') '&params'
if (who == 'sigma_cut') then   ! Old style
  write (iu, '(a)') ' local_space_charge_com%lsc_sigma_cut = ' // trim(val)
else
  write (iu, '(a)') ' local_space_charge_com%' // trim(who) // ' = ' // trim(val)
endif
write (iu, '(a)') '/'
rewind (iu)
local_space_charge_com = space_charge_com  ! set defaults
read (iu, nml = params, iostat = ios)
close (iu, status = 'delete')

if (ios /= 0) then
  call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
  return
endif

space_charge_com = local_space_charge_com
s%u%calc%lattice = .true.

end subroutine tao_set_space_charge_com_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_bmad_com_cmd (who, value_str)
!
! Routine to set bmad_com variables
! 
! Input:
!   who       -- Character(*): which bmad_com variable to set
!   value_str -- Character(*): Value to set to.
!-

subroutine tao_set_bmad_com_cmd (who, value_str)

implicit none

type (bmad_common_struct) this_bmad_com

character(*) who, value_str
character(*), parameter :: r_name = 'tao_set_bmad_com_cmd'

integer iu, ios
logical err

namelist / params / this_bmad_com

! open a scratch file for a namelist read

iu = tao_open_scratch_file (err);  if (err) return

write (iu, '(a)') '&params'
write (iu, '(a)') ' this_bmad_com%' // trim(who) // ' = ' // trim(value_str)
write (iu, '(a)') '/'
rewind (iu)
this_bmad_com = bmad_com  ! set defaults
read (iu, nml = params, iostat = ios)
close (iu, status = 'delete')

call tao_data_check (err)
if (err) return

if (ios == 0) then
  bmad_com = this_bmad_com
  s%u%calc%lattice = .true.
else
  call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
endif

end subroutine tao_set_bmad_com_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_ptc_com_cmd (who, value_str)
!
! Routine to set ptc_com variables
! 
! Input:
!   who       -- Character(*): which ptc_com variable to set
!   value_str -- Character(*): Value to set to.
!-

subroutine tao_set_ptc_com_cmd (who, value_str)

implicit none

character(*) who, value_str
character(*), parameter :: r_name = 'tao_set_ptc_com_cmd'

logical err

!

select case (who)
case ('vertical_kick');           call tao_set_real_value (ptc_com%vertical_kick, who, value_str, err)
case ('cut_factor');              call tao_set_real_value (ptc_com%cut_factor, who, value_str, err)
case ('max_fringe_order');        call tao_set_integer_value (ptc_com%max_fringe_order, who, value_str, err)
case ('old_integrator');          call tao_set_integer_value (ptc_com%old_integrator, who, value_str, err)
case ('exact_model');             call tao_set_logical_value (ptc_com%exact_model, who, value_str, err)
case ('exact_misalign');          call tao_set_logical_value (ptc_com%exact_misalign, who, value_str, err)
case ('use_orientation_patches'); call tao_set_logical_value (ptc_com%use_orientation_patches, who, value_str, err)
case ('print_info_messages');     call tao_set_logical_value (ptc_com%print_info_messages, who, value_str, err)
case default;                     call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
end select

end subroutine tao_set_ptc_com_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_geodesic_lm_cmd (who, value_str)
!
! Routine to set geodesic_lm variables
! 
! Input:
!   who       -- Character(*): which geodesic_lm variable to set
!   value_str -- Character(*): Value to set to.
!-

subroutine tao_set_geodesic_lm_cmd (who, value_str)

use geodesic_lm

implicit none

type (geodesic_lm_param_struct) this_geodesic_lm

character(*) who, value_str
character(*), parameter :: r_name = 'tao_set_geodesic_lm_cmd'

integer iu, ios
logical err

namelist / params / this_geodesic_lm

! open a scratch file for a namelist read

iu = tao_open_scratch_file (err);  if (err) return

write (iu, '(a)') '&params'
write (iu, '(a)') ' this_geodesic_lm%' // trim(who) // ' = ' // trim(value_str)
write (iu, '(a)') '/'
rewind (iu)
this_geodesic_lm = geodesic_lm_param  ! set defaults
read (iu, nml = params, iostat = ios)
close (iu, status = 'delete')

call tao_data_check (err)
if (err) return

if (ios == 0) then
  geodesic_lm_param = this_geodesic_lm
else
  call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
endif

end subroutine tao_set_geodesic_lm_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_opti_de_param_cmd (who, value_str)
!
! Routine to set opti_de_param variables
! 
! Input:
!   who       -- Character(*): which opti_de_param variable to set
!   value_str -- Character(*): Value to set to.
!-

subroutine tao_set_opti_de_param_cmd (who, value_str)

use opti_de_mod, only: opti_de_param

implicit none

character(*) who, value_str
character(*), parameter :: r_name = 'tao_set_opti_de_param_cmd'

integer iu, ios
logical err

namelist / params / opti_de_param

! open a scratch file for a namelist read

iu = tao_open_scratch_file (err);  if (err) return

write (iu, '(a)') '&params'
write (iu, '(a)') ' opti_de_param%' // trim(who) // ' = ' // trim(value_str)
write (iu, '(a)') '/'
rewind (iu)
read (iu, nml = params, iostat = ios)
close (iu, status = 'delete')

if (ios /= 0) then
  call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
endif

end subroutine tao_set_opti_de_param_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_wave_cmd (who, value_str, err)
!
! Routine to set wave variables
! 
! Input:
!   who       -- Character(*): which wave variable to set
!   value_str -- Character(*): Value to set to.
!
! Output:
!    err     -- logical: Set True if there is an error. False otherwise.
!    s%wave  -- Wave variables structure.
!-

subroutine tao_set_wave_cmd (who, value_str, err)

implicit none

type (tao_wave_struct) wave

character(*) who, value_str
character(*), parameter :: r_name = 'tao_set_wave_cmd'

real(rp) ix_a(2), ix_b(2)

integer iu, ios
logical err

namelist / params / ix_a, ix_b

! open a scratch file for a namelist read

iu = tao_open_scratch_file (err);  if (err) return

err = .true.

ix_a = [s%wave%ix_a1, s%wave%ix_a2]
ix_b = [s%wave%ix_b1, s%wave%ix_b2]

write (iu, '(a)') '&params'
write (iu, '(a)') trim(who) // ' = ' // trim(value_str)
write (iu, '(a)') '/'
rewind (iu)
wave = s%wave  ! set defaults
read (iu, nml = params, iostat = ios)
close (iu, status = 'delete')

if (ios /= 0) then
  call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
  return
endif

s%wave%ix_a1 = ix_a(1)
s%wave%ix_a2 = ix_a(2)
s%wave%ix_b1 = ix_b(1)
s%wave%ix_b2 = ix_b(2)

err = .false.

end subroutine tao_set_wave_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_beam_cmd (who, value_str, branch_str)
!
! Routine to set various beam parameters.
! 
! Input:
!   who         -- character(*): which parameter to set.
!   value_str   -- character(*): Value to set to.
!   branch_str  -- character(*): Branch to use. '' => branch 0.
!-

subroutine tao_set_beam_cmd (who, value_str, branch_str)

type (tao_universe_struct), pointer :: u
type (ele_pointer_struct), allocatable, target :: eles(:)
type (ele_struct), pointer :: ele
type (beam_struct), pointer :: beam
type (tao_beam_branch_struct), pointer :: bb

real(rp) value
integer ix, iu, n_loc, ie, ix_branch

logical, allocatable :: this_u(:)
logical err, logic, always_reinit

character(*) who, value_str, branch_str
character(20) switch, who2
character(*), parameter :: r_name = 'tao_set_beam_cmd'

!

call tao_pick_universe (unquote(who), who2, this_u, err); if (err) return

call match_word (who2, [character(32):: 'track_start', 'track_end', 'saved_at', 'comb_ds_save', 'comb_max_ds_save', &
                    'beam_track_start', 'beam_track_end', 'beam_init_file_name', 'beam_saved_at', &
                    'beginning', 'add_saved_at', 'subtract_saved_at', 'beam_init_position_file', &
                    'beam_dump_at', 'beam_dump_file', 'dump_at', 'dump_file', &
                    'always_reinit'], ix, matched_name=switch)

do iu = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. this_u(iu)) cycle
  u => s%u(iu)

  bb => u%model_branch(0)%beam

  if (switch /= 'beginning') then
    bb%init_starting_distribution = .true.
    u%calc%lattice = .true.
  endif

  select case (switch)
  case ('beginning')
    ele => tao_beam_track_endpoint (value_str, u%model%lat, '', 'BEGGINING', u)
    if (.not. associated(ele)) return
    beam => u%model_branch(ele%ix_branch)%ele(ele%ix_ele)%beam
    if (.not. allocated(beam%bunch)) then
      call out_io (s_error$, r_name, 'BEAM NOT SAVED AT: ' // who, 'NOTHING DONE.')
      err = .true.
      return
    endif
    bb%beam_at_start = beam
    u%calc%lattice = .true.

  case ('comb_ds_save')
    call tao_set_real_value (value, switch, value_str, err)
    u%model%tao_branch(:)%comb_ds_save = value

  case ('comb_max_ds_save')
    call tao_set_real_value (value, switch, value_str, err)
    u%model%tao_branch(:)%comb_max_ds_save = value

  case ('always_reinit')
    call tao_set_logical_value (u%beam%always_reinit, switch, value_str, err)

  case ('track_start', 'beam_track_start', 'track_end', 'beam_track_end')
    ele => tao_beam_track_endpoint (value_str, u%model%lat, branch_str, upcase(switch), u)
    if (.not. associated(ele)) return

    bb => u%model_branch(ele%ix_branch)%beam

    if (switch == 'track_start' .or. switch == 'beam_track_start') then
      bb%track_start = value_str
      bb%ix_track_start = ele%ix_ele
    else
      bb%track_end = value_str
      bb%ix_track_end = ele%ix_ele
    endif

  case ('beam_init_position_file', 'beam_init_file_name')
    if (switch == 'beam_init_file_name') call out_io (s_warn$, r_name, 'Note: "beam_init_file_name" has been renamed to "beam_init_position_file".')
    bb%beam_init%position_file = value_str
    bb%init_starting_distribution = .true.
  case ('dump_file', 'beam_dump_file')
    u%beam%dump_file = value_str

  case ('dump_at', 'beam_dump_at')
    call tao_locate_elements (value_str, u%ix_uni, eles, err, ignore_blank = .true.)
    if (err) then
      call out_io (s_error$, r_name, 'BAD DUMP_AT STRING: ' // value_str)
      return
    endif
    u%beam%dump_at = value_str

    do ix = 0, ubound(u%model_branch, 1)
      u%model_branch(ix)%ele(:)%save_beam_to_file = .false.
    enddo

    ! Note: Beam will automatically be dump at fork elements and at the ends of the beam tracking.
    do ix = 1, size(eles)
      ele => eles(ix)%ele
      u%model_branch(ele%ix_branch)%ele(ele%ix_ele)%save_beam_to_file = .true.
    enddo

  case ('saved_at', 'beam_saved_at')
    call tao_locate_elements (value_str, u%ix_uni, eles, err, ignore_blank = .true.)
    if (err) then
      call out_io (s_error$, r_name, 'BAD SAVED_AT STRING: ' // value_str)
      return
    endif
    u%beam%saved_at = value_str

    do ix = 0, ubound(u%model_branch, 1)
      u%model_branch(ix)%ele(:)%save_beam_internally = .false.
    enddo

    ! Note: Beam will automatically be saved at fork elements and at the ends of the beam tracking.
    do ix = 1, size(eles)
      ele => eles(ix)%ele
      u%model_branch(ele%ix_branch)%ele(ele%ix_ele)%save_beam_internally = .true.
    enddo

  case ('add_saved_at', 'subtract_saved_at')
    call tao_locate_elements (value_str, u%ix_uni, eles, err)
    if (err) then
      call out_io (s_error$, r_name, 'BAD SAVED_AT STRING: ' // value_str)
      return
    endif

    logic = (switch == 'add_saved_at')
    do ix = 1, size(eles)
      ele => eles(ix)%ele
      u%model_branch(ele%ix_branch)%ele(ele%ix_ele)%save_beam_internally = logic
    enddo

  case default
    call out_io (s_error$, r_name, 'PARAMETER NOT RECOGNIZED: ' // who2, &
                                   '[DID YOU WANT "SET BEAM_INIT" INSTEAD OF "SET BEAM"?]')
    return
  end select
enddo

end subroutine tao_set_beam_cmd 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_beam_init_cmd (who, value_str, branch_str)
!
! Routine to set beam_init variables
! 
! Input:
!   who         -- character(*): which beam_init variable to set
!   value_str   -- character(*): Value to set to.
!   branch_str  -- character(*): Branch to use. '' => branch 0
!
! Output:
!    s%beam_init  -- Beam_init variables structure.
!-

subroutine tao_set_beam_init_cmd (who, value_str, branch_str)

implicit none

type (beam_init_struct) beam_init
type (tao_universe_struct), pointer :: u
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), pointer :: ele
type (tao_beam_branch_struct), pointer :: bb

character(*) who, value_str, branch_str
character(40) who2
character(*), parameter :: r_name = 'tao_set_beam_init_cmd'

real(rp), allocatable :: set_val(:)
real(rp) r_val
integer i, ix, iu, ios, ib, n_loc
logical err, eval_err
logical, allocatable :: picked_uni(:)

character(40) name, switch

namelist / params / beam_init

! get universe

call tao_pick_universe (unquote(who), who2, picked_uni, err)
call downcase_string(who2)

! Special cases not associated with the beam_init structure

call match_word (who2, [character(32):: 'track_start', 'track_end', 'saved_at', 'beam_saved_at', &
                    'beam_track_start', 'beam_track_end', 'beam_init_file_name', &
                    'beginning', 'add_saved_at', 'subtract_saved_at', 'beam_init_position_file', &
                    'beam_dump_at', 'beam_dump_file', 'dump_at', 'dump_file', &
                    'always_reinit'], ix, matched_name=switch)

if (ix > 0) then
  call tao_set_beam_cmd(who, value_str, branch_str)
  return
endif

! Beam_init. open a scratch file for a namelist read

eval_err = .false.
if (who2 == 'sig_e') who2 = 'sig_pz'
iu = tao_open_scratch_file (err);  if (err) return

write (iu, '(a)') '&params'

if (is_real(value_str, real_num = r_val) .or. is_logical(value_str)) then
  if (who2 == 'n_particle') then  ! Sometimes people use "1E3" for the value
    write (iu, '(a, i0)') ' beam_init%' // trim(who2) // ' = ', nint(r_val)
  else
    write (iu, '(2a)') ' beam_init%' // trim(who2) // ' = ', trim(value_str)
  endif

elseif (who(1:17) == 'distribution_type') then  ! Value is a vector so quote() function is not good here.
  write (iu, '(2a)') ' beam_init%' // trim(who2) // ' = ', value_str

else
  select case (who2)
  case ('random_engine')
    call downcase_string(value_str)
    call match_word(unquote(value_str), random_engine_name, ix, .true., .false.)
    if (ix == 0) then
      call out_io (s_error$, r_name, 'INVALID RANDOM_ENGINE VALUE: ' // value_str)
      return
    endif
    write (iu, '(2a)') ' beam_init%' // trim(who2) // ' = ', quote(value_str)

  case ('random_gauss_converter')
    call downcase_string(value_str)
    call match_word(unquote(value_str), random_gauss_converter_name, ix, .true., .false.)
    if (ix == 0) then
      call out_io (s_error$, r_name, 'INVALID RANDOM_GAUSS_CONVERTER VALUE: ' // value_str)
      return
    endif
    write (iu, '(2a)') ' beam_init%' // trim(who2) // ' = ', quote(value_str)

  case ('species')
    if (species_id(unquote(value_str)) == invalid$) then
      call out_io (s_error$, r_name, 'INVALID SPECIES VALUE: ' // value_str)
      return
    endif
    write (iu, '(2a)') ' beam_init%' // trim(who2) // ' = ', quote(value_str)

  case ('position_file')
    write (iu, '(2a)') ' beam_init%' // trim(who2) // ' = ', quote(value_str)

  case ('center_jitter', 'center', 'spin')
    write (iu, '(2a)') ' beam_init%' // trim(who2) // ' = ', value_str

  case default
    ! If tao_evaluate_expression fails then the root cause may be that the User is trying
    ! something like "set beam_init beam_saved_at = END" so the basic problem is that beam_saved_at
    ! is not a valid beam_init component. So delay error messages until we know for sure.
    call tao_evaluate_expression (value_str, 1, .false., set_val, eval_err, print_err = .false.)
    if (eval_err) then
      write (iu, '(a)') ' beam_init%' // trim(who2) // ' = 0'  ! For a test
    else
      write (iu, '(a, es23.15)') ' beam_init%' // trim(who2) // ' = ', set_val(1)
    endif
  end select
endif

write (iu, '(a)') '/'

!

do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. picked_uni(i)) cycle

  rewind (iu)
  u => s%u(i)
  bb => u%model_branch(0)%beam
  beam_init = bb%beam_init  ! set defaults
  read (iu, nml = params, iostat = ios)
  if (ios /= 0 .or. eval_err) then
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'BAD BEAM_INIT COMPONENT: ' // who2)
    else
      call tao_evaluate_expression (value_str, 1, .false., set_val, eval_err, print_err = .true.) ! Print error message
    endif
    exit
  endif

  select case (who2)
  case ('a_emit');        beam_init%a_norm_emit = 0
  case ('b_emit');        beam_init%b_norm_emit = 0
  case ('a_norm_emit');   beam_init%a_emit = 0
  case ('b_norm_emit');   beam_init%b_emit = 0
  end select

  bb%beam_init = beam_init
  bb%init_starting_distribution = .true.  ! Force reinit
  if (bb%beam_init%use_particle_start) then
    bb%beam_init%center = u%model%lat%particle_start%vec
    bb%beam_init%spin = u%model%lat%particle_start%spin
  endif
  u%calc%lattice = .true.
  u%beam%track_beam_in_universe = .true.
enddo

close (iu, status = 'delete') 
deallocate (picked_uni)

end subroutine tao_set_beam_init_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_particle_start_cmd (who, value_str)
!
! Routine to set particle_start variables.
! 
! Input:
!   who       -- Character(*): which particle_start variable to set
!   value_str -- Character(*): Value to set to.
!
! Output:
!    s%particle_start  -- Beam_start variables structure.
!-

subroutine tao_set_particle_start_cmd (who, value_str)

type (tao_universe_struct), pointer :: u
type (all_pointer_struct), allocatable :: a_ptr(:)
type (tao_d2_data_array_struct), allocatable :: d2_array(:)

real(rp), allocatable :: set_val(:)

integer ix, iu

character(*) who, value_str
character(40) who2, name

character(*), parameter :: r_name = 'tao_set_particle_start_cmd'

logical, allocatable :: this_u(:)
logical err

!

call tao_pick_universe (who, who2, this_u, err); if (err) return
call string_trim (upcase(who2), who2, ix)

if (who2 == '*') then
  call tao_evaluate_expression (value_str, 6, .false., set_val, err); if (err) return
else
  call tao_evaluate_expression (value_str, 1, .false., set_val, err); if (err) return
endif

!

do iu = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. this_u(iu)) cycle
  u => s%u(iu)

  if (who2 == '*') then
    u%model%lat%particle_start%vec = set_val

  else
    call pointers_to_attribute (u%model%lat, 'PARTICLE_START', who2, .true., a_ptr, err, .true.)
    if (err) return

    if (u%model%lat%param%geometry == closed$) then
      ! All parameters free to vary if multi-turn orbit data is being collected.
      if (.not. s%com%multi_turn_orbit_is_plotted) then
        if (who2 == 'PZ') then
          if (s%global%rf_on) then
            call out_io (s_warn$, r_name, 'Setting particle_start[pz] will not affect lattice calculations since the rf is on and the lattice is closed.', &
                              'Note: "set global rf_on = F" or "set branch 0 geometry = open" will change this situation.')
          endif
        else
          call out_io (s_error$, r_name, 'Setting particle_start non-pz parameter will not affect lattice calculations since the lattice is closed.', &
                              'Note: "set branch 0 geometry = open" will change this situation.')
        endif
      endif
    endif

    ! Set value

    a_ptr(1)%r = set_val(1)
  endif

  call tao_set_flags_for_changed_attribute (u, 'PARTICLE_START')
enddo

end subroutine tao_set_particle_start_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Subroutine tao_set_plot_page_cmd (component, value_str, value_str2)
!
!  Set various aspects of the plotting window
!
! Input:
!   component     -- Character(*): Which component to set.
!   value_str     -- Character(*): What value to set to.
!   value_str2    -- Character(*): 2nd value if component is an array.
!
!  Output:
!    s%plot       -- tao_plotting_struct:
!-

subroutine tao_set_plot_page_cmd (component, value_str, value_str2)

use tao_input_struct, only: tao_plot_page_input, tao_set_plotting
use tao_plot_window_mod, only: tao_destroy_plot_window, tao_create_plot_window

implicit none

type (tao_plot_page_input) plot_page

character(*) component, value_str
character(*), optional :: value_str2
character(24) :: r_name = 'tao_set_plot_page_cmd'

real(rp) x, y
integer iu, ios
logical err


namelist / params / plot_page

! Special cases

select case (component)

case ('title')
  s%plot_page%title%string = trim(value_str)
  return

case ('subtitle')
  s%plot_page%subtitle%string = trim(value_str)
  return

end select

! For everything else...
! open a scratch file for a namelist read

iu = tao_open_scratch_file (err);  if (err) return

write (iu, '(a)') '&params'

select case (component)
case ('subtitle%string', 'subtitle%units', 'subtitle%justify', 'plot_display_type', &
      'title%string', 'title%units', 'title%justify')
  write (iu, '(a)') ' plot_page%' // trim(component) // ' = ' // quote(value_str)
case ('size')
  write (iu, '(a)') ' plot_page%' // trim(component) // ' = ' // trim(value_str) // ',' // trim(value_str2)
case default
  write (iu, '(a)') ' plot_page%' // trim(component) // ' = ' // trim(value_str)
end select

write (iu, '(a)') '/'
rewind (iu)

call tao_set_plotting (plot_page, s%plot_page, .false., .true.)

read (iu, nml = params, iostat = ios)
close (iu, status = 'delete')

if (ios /= 0) then
  call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
  return
endif

call tao_set_plotting (plot_page, s%plot_page, .false.)

if (component == 'size') then
  call tao_destroy_plot_window()
  call tao_create_plot_window()
endif

end subroutine tao_set_plot_page_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_curve_cmd (curve_name, component, value_str)
!
! Routine to set var values.
!
! Input:
!   curve_name -- Character(*): Which curve to set.
!   component  -- Character(*): Which component to set.
!   value_str  -- Character(*): What value to set it to.
!-

subroutine tao_set_curve_cmd (curve_name, component, value_str)

implicit none

type (tao_curve_array_struct), allocatable :: curve(:)
type (tao_graph_array_struct), allocatable :: graph(:)
type (lat_struct), pointer :: lat

integer i, j, ios, i_uni
integer, allocatable :: ix_ele(:)

character(*) curve_name, component, value_str
character(*), parameter :: r_name = 'tao_set_curve_cmd'

logical err

!

call tao_find_plots (err, curve_name, 'BOTH', curve = curve, blank_means_all = .true., only_visible = .false.)
if (err) return

if (size(curve) == 0) then
  call out_io (s_error$, r_name, 'CURVE OR GRAPH NOT SPECIFIED')
  return
else
  do i = 1, size(curve)
    call set_this_curve (curve(i)%c)
  enddo
endif

!---------------------------------------------
contains

subroutine set_this_curve (this_curve)

type (tao_curve_struct) this_curve
type (tao_graph_struct), pointer :: this_graph
type (tao_universe_struct), pointer :: u
type (tao_model_branch_struct), pointer :: model_branch
type (ele_pointer_struct), allocatable :: eles(:)

integer ix, i_branch
logical err, is_int
character(40) name, comp

!

i_branch = tao_branch_index(this_curve%ix_branch)
i_uni = tao_universe_index(tao_curve_ix_uni(this_curve))

this_graph => this_curve%g
this_graph%p%default_plot = .true.   ! Plot has been modified

! if the universe is changed then need to check ele_ref

comp = component
ix = index(comp, '.')
if (ix /= 0) comp(ix:ix) = '%'
select case (comp)

case ('ele_ref_name', 'ix_ele_ref')
  is_int = is_integer(value_str, ix)
  if (value_str == '' .or. (is_int .and. ix < 0)) then
    this_curve%ele_ref_name = ''
    this_curve%ix_ele_ref = -1

  else
    call tao_locate_elements (value_str, i_uni, eles, err, ignore_blank = .true.)
    if (size(eles) == 0) return
    this_curve%ele_ref_name = upcase(value_str)
    this_curve%ix_ele_ref = eles(1)%ele%ix_ele
    this_curve%ix_branch  = eles(1)%ele%ix_branch
    call tao_ele_to_ele_track (i_uni, i_branch, this_curve%ix_ele_ref, this_curve%ix_ele_ref_track)
  endif

case ('name')
  this_curve%name = value_str
  
case ('ix_universe')
  call tao_set_integer_value (this_curve%ix_universe, component, value_str, err, -2, ubound(s%u, 1))
  if (err) return
  call tao_locate_elements (this_curve%ele_ref_name, tao_curve_ix_uni(this_curve), eles, err, ignore_blank = .true.)
  if (size(eles) == 0) return
  this_curve%ix_ele_ref = eles(1)%ele%ix_ele
  this_curve%ix_branch  = eles(1)%ele%ix_branch
  call tao_ele_to_ele_track (tao_curve_ix_uni(this_curve), this_curve%ix_branch, &
                                     this_curve%ix_ele_ref, this_curve%ix_ele_ref_track)

case ('ix_branch') 
  call tao_set_integer_value (this_curve%ix_branch, component, value_str, err, -1, ubound(s%u(i_uni)%model%lat%branch, 1))

case ('ix_bunch')
  u => tao_pointer_to_universe (tao_curve_ix_uni(this_curve))
  if (.not. associated(u)) return
  call tao_set_integer_value (this_curve%ix_bunch, component, value_str, err, 0, u%model_branch(0)%beam%beam_init%n_bunch)

case ('symbol_every')
  call tao_set_integer_value (this_curve%symbol_every, component, value_str, err, 0, 1000000)

case ('symbol_size')
  call tao_set_real_value (this_curve%symbol%height, component, value_str, err)

case ('symbol_color', 'symbol%color')
  call tao_set_switch_value (ix, component, value_str, qp_color_name, lbound(qp_color_name,1), err, this_curve%symbol%color)

case ('symbol_type', 'symbol%type')
  call tao_set_switch_value (ix, component, value_str, qp_symbol_type_name, lbound(qp_symbol_type_name,1), err, this_curve%symbol%type)

case ('symbol_fill_pattern', 'symbol%fill_pattern')
  call tao_set_switch_value (ix, component, value_str, qp_symbol_fill_pattern_name, &
                                                              lbound(qp_symbol_fill_pattern_name,1), err, this_curve%symbol%fill_pattern)

case ('symbol_height', 'symbol%height')
  call tao_set_real_value (this_curve%symbol%height, component, value_str, err)

case ('symbol_line_width', 'symbol%line_width')
  call tao_set_integer_value (this_curve%symbol%line_width, component, value_str, err)

case ('smooth_line_calc')
  call tao_set_logical_value (this_curve%smooth_line_calc, component, value_str, err)

case ('line_color', 'line%color')
  call tao_set_switch_value (ix, component, value_str, qp_color_name, lbound(qp_color_name,1), err, this_curve%line%color)

case ('line_width', 'line%width')
  call tao_set_integer_value (this_curve%line%width, component, value_str, err)

case ('line_pattern', 'line%pattern')
  call tao_set_switch_value (ix, component, value_str, qp_line_pattern_name, lbound(qp_line_pattern_name,1), err, this_curve%line%pattern)

case ('component')
  this_curve%component = unquote(value_str)

case ('draw_error_bars')
  call tao_set_logical_value (this_curve%draw_error_bars, component, value_str, err)

case ('draw_line')
  call tao_set_logical_value (this_curve%draw_line, component, value_str, err)

case ('draw_symbols')
  call tao_set_logical_value (this_curve%draw_symbols, component, value_str, err)

case ('draw_symbol_index')
  call tao_set_logical_value (this_curve%draw_symbol_index, component, value_str, err)

case ('use_y2')
  call tao_set_logical_value (this_curve%use_y2, component, value_str, err)

case ('data_source')
  this_curve%data_source = unquote(value_str)

case ('data_index')
  this_curve%data_index = unquote(value_str)

case ('data_type')
  this_curve%data_type = unquote(value_str)

case ('data_type_x')
  this_curve%data_type_x = unquote(value_str)

case ('legend_text')
  this_curve%legend_text = unquote(value_str)

case ('z_color%data_type', 'data_type_z')
  this_curve%z_color%data_type = unquote(value_str)

case ('z_color%is_on', 'use_z_color')
  call tao_set_logical_value (this_curve%z_color%is_on, component, value_str, err)
  
case ('z_color%autoscale', 'autoscale_z_color')
  call tao_set_logical_value (this_curve%z_color%autoscale, component, value_str, err)  

case ('z_color%min', 'z_color0')
  call tao_set_real_value (this_curve%z_color%min, component, value_str, err, dflt_uni = i_uni)

case ('z_color%max', 'z_color1')
  call tao_set_real_value (this_curve%z_color%max, component, value_str, err, dflt_uni = i_uni) 

case ('hist%number')
  this_curve%hist%width = 0
  call tao_set_integer_value (this_curve%hist%number, component, value_str, err, min_val = 0)

case ('hist%density_normalized')
  call tao_set_logical_value (this_curve%hist%density_normalized, component, value_str, err)
  
case ('hist%weight_by_charge')
  call tao_set_logical_value (this_curve%hist%weight_by_charge, component, value_str, err)
  
case ('hist%center')  
  call tao_set_real_value (this_curve%hist%center, component, value_str, err, dflt_uni = i_uni)
  
case ('hist%width')  
  this_curve%hist%number = 0
  call tao_set_real_value (this_curve%hist%width, component, value_str, err, dflt_uni = i_uni)  
  
case ('orbit')
  call string_trim (value_str, value_str, ix)
  call tao_set_real_value (this_curve%orbit%x, component, value_str(1:ix), err, dflt_uni = i_uni)
  if (err) return
  call string_trim (value_str(ix+1:), value_str, ix)
  if (ix == 0) return
  call tao_set_real_value (this_curve%orbit%y, component, value_str(1:ix), err, dflt_uni = i_uni)
  if (err) return
  call string_trim (value_str(ix+1:), value_str, ix)
  if (ix == 0) return
  call tao_set_real_value (this_curve%orbit%t, component, value_str(1:ix), err, dflt_uni = i_uni)

case ('orbit%x')
  call tao_set_real_value (this_curve%orbit%x, component, value_str, err, dflt_uni = i_uni)
case ('orbit%y')
  call tao_set_real_value (this_curve%orbit%y, component, value_str, err, dflt_uni = i_uni)
case ('orbit%t')
  call tao_set_real_value (this_curve%orbit%t, component, value_str, err, dflt_uni = i_uni)

case ('y_axis_scale_factor')
  call tao_set_real_value (this_curve%y_axis_scale_factor, component, value_str, err, dflt_uni = i_uni)

case default
  call out_io (s_error$, r_name, "BAD CURVE COMPONENT")
  return

end select

! Enable

if (this_graph%type == 'phase_space') then
  model_branch => s%u(i_uni)%model_branch(i_branch)
  if (.not. model_branch%ele(this_curve%ix_ele_ref)%save_beam_internally) then
    s%u(i_uni)%calc%lattice = .true.
    model_branch%ele(this_curve%ix_ele_ref)%save_beam_internally = .true.
  endif
endif

end subroutine set_this_curve

end subroutine tao_set_curve_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_plot_cmd (plot_name, component, value_str)
!
! Routine to set plot parameters.
!
! Input:
!   plot_name  -- character(*): Which plot to set.
!   component  -- character(*): Which component to set.
!   value_str  -- character(*): What value to set it to.
!
!  Output:
!-

subroutine tao_set_plot_cmd (plot_name, component, value_str)

implicit none

type (tao_plot_array_struct), allocatable :: plot(:)
type (tao_universe_struct), pointer :: u
type (tao_plot_struct), pointer :: p

character(*) plot_name, component, value_str
character(40) comp, sub_comp
character(*), parameter :: r_name = 'tao_set_plot_cmd'

integer iset, iw, iu
integer i, j, ix, ios
logical err_flag, found

!

call tao_find_plots (err_flag, plot_name, 'BOTH', plot = plot, only_visible = .false.)
if (err_flag) return

if (size(plot) == 0) then
  call out_io (s_error$, r_name, 'PLOT OR PLOT NOT SPECIFIED')
  return
endif

! And set

comp = component
ix = index(component, '%')
if (component(1:17) == 'floor_plan_orbit_') ix = 17   ! Old style.

if (ix /= 0) then  ! Split on '%'
  comp = component(:ix-1)
  sub_comp = component(ix+1:)
endif

found = .false.

do i = 1, size(plot)
  p => plot(i)%p
  p%default_plot = .false.  ! Plot has been modified

  select case (comp)

  case ('autoscale_x')
    call tao_set_logical_value (p%autoscale_x, component, value_str, err_flag)

  case ('autoscale_y')
    call tao_set_logical_value (p%autoscale_y, component, value_str, err_flag)

  case ('autoscale_gang_x')
    call tao_set_logical_value (p%autoscale_gang_x, component, value_str, err_flag)

  case ('autoscale_gang_y')
    call tao_set_logical_value (p%autoscale_gang_y, component, value_str, err_flag)

  case ('description')
    p%description = value_str

  case ('component')
    do j = 1, size(p%graph)
      if (.not. allocated(p%graph(j)%curve)) cycle
      p%graph(j)%curve%component = unquote(value_str)
    enddo

  case ('n_curve_pts')
    call tao_set_integer_value (p%n_curve_pts, component, value_str, err_flag)

  case ('name')
    p%name = value_str

  case ('visible')
    if (.not. associated(p%r)) cycle
    call tao_set_logical_value (p%r%visible, component, value_str, err_flag)
    call tao_turn_on_special_calcs_if_needed_for_plotting()
    found = .true.

  case ('x')
    if (allocated(p%graph)) then
      do j = 1, size(p%graph)
        call tao_set_qp_axis_struct('x', sub_comp, p%graph(j)%x, value_str, err_flag)
        if (err_flag) exit
      enddo
    endif

  case ('x_axis_type')
    call tao_set_switch_value (ix, component, value_str, tao_x_axis_type_name, lbound(tao_x_axis_type_name,1), err_flag)
    if (.not. err_flag) p%x_axis_type = tao_x_axis_type_name(ix)

  case default
    call out_io (s_error$, r_name, "BAD PLOT COMPONENT: " // component)
    return
      
  end select

enddo

!

if (comp == 'visible' .and. .not. found) then
  call out_io (s_error$, r_name, 'NO PLOT ASSOCIATED WITH: ' // plot_name)
endif

end subroutine tao_set_plot_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_region_cmd (region_name, component, value_str)
!
! Routine to set region parameters.
!
! Input:
!   region_name -- character(*): Which region to set.
!   component   -- character(*): Which component to set.
!   value_str   -- character(*): What value to set it to.
!
!  Output:
!-

subroutine tao_set_region_cmd (region_name, component, value_str)

implicit none

type (tao_universe_struct), pointer :: u
type (tao_plot_region_struct), pointer :: r

character(*) region_name, component, value_str
character(*), parameter :: r_name = 'tao_set_region_cmd'

integer iset, iw, iu
integer i, j, ix, ios
logical err_flag, found

!

call tao_find_plot_region (err_flag, region_name, r)
if (err_flag) return

if (.not. associated(r)) then
  call out_io (s_error$, r_name, 'BAD REGION NAME.')
  return
endif

! And set

select case (component)
case ('x1')
  call tao_set_real_value (r%location(1), component, value_str, err_flag)
case ('x2')
  call tao_set_real_value (r%location(2), component, value_str, err_flag)
case ('y1')
  call tao_set_real_value (r%location(3), component, value_str, err_flag)
case ('y2')
  call tao_set_real_value (r%location(4), component, value_str, err_flag)
case ('visible')
  call tao_set_logical_value (r%visible, component, value_str, err_flag)
end select

end subroutine tao_set_region_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_graph_cmd (graph_name, component, value_str)
!
! Routine to set var values.
!
! Input:
!   graph_name -- Character(*): Which graph to set.
!   component  -- Character(*): Which component to set.
!   value_str  -- Character(*): What value to set it to.
!
!  Output:
!-

subroutine tao_set_graph_cmd (graph_name, component, value_str)

implicit none

type (tao_plot_array_struct), allocatable :: plot(:)
type (tao_graph_array_struct), allocatable :: graph(:)

character(*) graph_name, component, value_str
character(*), parameter :: r_name = 'tao_set_graph_cmd'

integer i, j, ios
logical err

! 'BOTH' was 'REGION'. Not sure why.

call tao_find_plots (err, graph_name, 'BOTH', plot = plot, graph = graph, only_visible = .false.)
if (err) return

if (size(graph) > 0) then
  do i = 1, size(graph)
    call set_this_graph (graph(i)%g)
  enddo
elseif (size(plot) > 0) then
  do i = 1, size(plot)
    do j = 1, size(plot(i)%p%graph)
      call set_this_graph (plot(i)%p%graph(j))
    enddo
  enddo
else
  call out_io (s_error$, r_name, 'GRAPH OR PLOT NOT SPECIFIED')
  return
endif

!---------------------------------------------
contains

subroutine set_this_graph (this_graph)

type (tao_graph_struct) this_graph
type (tao_universe_struct), pointer :: u
character(40) comp, sub_comp
character(200) value
integer iset, iw, ix, ix2
logical logic, error

!

value = unquote(value_str)

comp = component
sub_comp = ''

select case (comp)
case ('floor_plan_size_is_absolute');      comp = 'floor_plan%size_is_absolute'
case ('floor_plan_draw_only_first_pass');  comp = 'floor_plan%draw_only_first_pass'
case ('floor_plan_flip_label_side');       comp = 'floor_plan%flip_label_side'
case ('floor_plan_rotation');              comp = 'floor_plan%rotation'
case ('floor_plan_scale');                 comp = 'floor_plan%scale'
case ('floor_plan_color');                 comp = 'floor_plan%color'
case ('floor_plan_pattern');               comp = 'floor_plan%pattern'
case ('floor_plan_view');                  comp = 'floor_plan%view'
case ('floor_plan_width');                 comp = 'floor_plan%width'
case ('correct_xy_distortion');            comp = 'floor_plan%correct_distortion'
end select

ix = max(index(comp, '%'), index(comp, '.'), index(comp, '('))
if (ix /= 0) then
  sub_comp = comp(ix+1:)
  comp = comp(:ix-1)
endif

u => tao_pointer_to_universe(this_graph%ix_universe, .true.)
this_graph%p%default_plot = .false. ! Plot has been modified

select case (comp)
case ('text_legend')
  ix = len_trim(sub_comp)
  call tao_set_integer_value (iw, component, sub_comp(:ix-1), error, 1, size(this_graph%text_legend));  if (error) return
  this_graph%text_legend(iw) = value

case ('allow_wrap_around')
  call tao_set_logical_value (this_graph%allow_wrap_around, component, value, error)

case ('box')
  ix = index(component, '(')
  if (ix == 0) then
    read (value, *, iostat = ios) this_graph%box
  else
    ix2 = index(component, ')')
    if (ix2 == 0) then
      call out_io (s_error$, r_name, 'MISSING ")" CHARACTER IN "BOX(N)" CONSTRUCT.')
      return
    endif
    error = (.not. is_integer(component(ix+1:ix2-1), iw))
    if (error .or. iw < 0 .or. iw > 4) then
      call out_io (s_error$, r_name, 'BAD BOX INDEX.')
      return
    endif
    read (value, *, iostat = ios) this_graph%box(iw)
  endif
  if (ios /= 0) then
    call out_io (s_error$, r_name, 'BAD BOX INTEGER VALUE.')
  endif

case ('component')
  if (allocated(this_graph%curve)) this_graph%curve%component = value_str
case ('clip')
  call tao_set_logical_value (this_graph%clip, component, value, error)
case ('curve_legend_origin')
  call tao_set_qp_point_struct (comp, sub_comp, this_graph%curve_legend_origin, value, error, u%ix_uni)
case ('draw_axes')
  call tao_set_logical_value (this_graph%draw_axes, component, value, error)
case ('draw_title')
  call tao_set_logical_value (this_graph%draw_title, component, value, error)
case ('draw_curve_legend')
  call tao_set_logical_value (this_graph%draw_curve_legend, component, value, error)
case ('draw_grid')
  call tao_set_logical_value (this_graph%draw_grid, component, value, error)
case ('draw_only_good_user_data_or_vars')
  call tao_set_logical_value (this_graph%draw_only_good_user_data_or_vars, component, value, error)
case ('floor_plan')
  select case (sub_comp)
  case ('correct_distortion')
    call tao_set_logical_value(this_graph%floor_plan%correct_distortion, component, value, error)
  case ('size_is_absolute')
    call tao_set_logical_value(this_graph%floor_plan%size_is_absolute, component, value, error)
  case ('draw_only_first_pass')
    call tao_set_logical_value(this_graph%floor_plan%draw_only_first_pass, component, value, error)
  case ('flip_label_side')
    call tao_set_logical_value(this_graph%floor_plan%flip_label_side, component, value, error)
  case ('rotation')
    call tao_set_real_value(this_graph%floor_plan%rotation, component, value, error, dflt_uni = u%ix_uni)
  case ('orbit_scale', 'scale')
    call tao_set_real_value(this_graph%floor_plan%orbit_scale, component, value, error, dflt_uni = u%ix_uni)
  case ('orbit_color', 'color')
    this_graph%floor_plan%orbit_color = value
  case ('orbit_lattice', 'lattice')
    this_graph%floor_plan%orbit_lattice = value
  case ('orbit_pattern', 'pattern')
    this_graph%floor_plan%orbit_pattern = value
  case ('view')
    if (.not. any(value == tao_floor_plan_view_name)) then
      call out_io(s_info$, r_name, "Valid floor_plan_view settings are: 'xy', 'zx', etc.")
      return
    endif
    this_graph%floor_plan%view = value
  case ('orbit_width', 'width')
    call tao_set_integer_value(this_graph%floor_plan%orbit_width, component, value, error, 1, 999)
  case default
    call out_io (s_error$, r_name, "BAD GRAPH floor_plan SUB-COMPONENT: " // sub_comp)
    return
  end select
case ('ix_universe')
  if (this_graph%type == 'floor_plan' .or. this_graph%type == 'lat_layout') then
    call tao_set_integer_value (this_graph%ix_universe, component, value, error, -2, ubound(s%u, 1))
    u => tao_pointer_to_universe(this_graph%ix_universe, .true.)
  else
    call tao_set_integer_value (this_graph%ix_universe, component, value, error, -1, ubound(s%u, 1))
    u => tao_pointer_to_universe(this_graph%ix_universe)
  endif
case ('ix_branch')
  call tao_set_integer_value (this_graph%ix_branch, component, value, error, -1, ubound(u%model%lat%branch, 1))
case ('margin')
  call tao_set_qp_rect_struct (comp, sub_comp, this_graph%margin, value, error, u%ix_uni)
case ('name')
  this_graph%name = value_str
case ('scale_margin')
  call tao_set_qp_rect_struct (comp, sub_comp, this_graph%scale_margin, value, error, u%ix_uni)
case ('symbol_size_scale')
  call tao_set_real_value(this_graph%symbol_size_scale, component, value, error, dflt_uni = u%ix_uni)
case ('text_legend_origin')
  call tao_set_qp_point_struct (comp, sub_comp, this_graph%text_legend_origin, value, error, u%ix_uni)
case ('title')
  this_graph%title = value
case ('type')
  this_graph%type = value
case ('x')
  call tao_set_qp_axis_struct (comp, sub_comp, this_graph%x, value, error, u%ix_uni)
case ('y')
  call tao_set_qp_axis_struct (comp, sub_comp, this_graph%y, value, error, u%ix_uni)
case ('x2')
  call tao_set_qp_axis_struct (comp, sub_comp, this_graph%x2, value, error, u%ix_uni)
case ('y2')
  call tao_set_qp_axis_struct (comp, sub_comp, this_graph%y2, value, error, u%ix_uni)
case ('y2_mirrors_y')
  call tao_set_logical_value(this_graph%y2_mirrors_y, component, value, error)
case ('x_axis_scale_factor')
  call tao_set_real_value(this_graph%x_axis_scale_factor, component, value, error, dflt_uni = u%ix_uni)
case default
  call out_io (s_error$, r_name, "BAD GRAPH COMPONENT: " // component)
  return
end select

if (associated(u)) u%calc%lattice = .true.

end subroutine set_this_graph

end subroutine tao_set_graph_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_var_cmd (var_str, value_str)
!
! Routine to set var values.
!
! Input:
!   var_str  -- Character(*): Which var name to set.
!   value_str  -- Character(*): What value to set it to.
!
!  Output:
!-

subroutine tao_set_var_cmd (var_str, value_str)

implicit none

type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_real_pointer_struct), allocatable    :: r_var(:), r_set(:)
type (tao_logical_array_struct), allocatable :: l_var(:), l_set(:)
type (tao_var_array_struct), allocatable, target :: v_var(:)
type (tao_string_array_struct), allocatable :: s_var(:), s_set(:)
type (tao_universe_struct), pointer :: u
type (ele_pointer_struct), allocatable :: eles(:)
type (all_pointer_struct) a_ptr
type (tao_var_struct), pointer :: v_ptr

real(rp), allocatable :: r_value(:)
real(rp) value

integer i, j, ix, np, n_loc
integer, allocatable :: u_pick(:)

character(*) var_str, value_str
character(*), parameter :: r_name = 'tao_set_var_cmd'
character(20) set_is, component
character(40) ele_name, attrib_name
character(len(value_str)) val

logical err, l_value, err_flag

! Decode variable component to set.

call tao_find_var (err, var_str, v_array = v_var, re_array=r_var, log_array=l_var, &
                                                 str_array = s_var, component = component)
if (err) return

select case (component)
case ('base')
  call out_io (s_error$, r_name, &
        'VARIABLES IN THE BASE LATTICE ARE NOT ALLOWED TO BE SET DIRECTLY SINCE DEPENDENT', &
        'PARAMETERS (LIKE THE TWISS PARAMETERS) IN THE BASE LATTICE ARE NEVER COMPUTED.', &
        'USE THE "SET LATTICE BASE = ..." COMMAND INSTEAD.')
  return

case ('ele_name', 'attrib_name', 'design', 'old', 'base_value', &
      'design_value', 'old_value', 'merit', 'delta_merit', 'exists', 'good_var', 'useit_opt', &
      'useit_plot')
  call out_io (s_error$, r_name, 'VARIABLE ATTRIBUTE NOT SETTABLE: ' // component)
  return
end select

! A logical value_str is either a logical or an array of datum values.

if (size(l_var) > 0) then
  if (is_logical(value_str)) then
    read (value_str, *) l_value
    do i = 1, size(l_var)
      l_var(i)%l = l_value
    enddo

  else
    call tao_find_var (err, value_str, log_array=l_set)
    if (size(l_set) /= size(l_var)) then
      call out_io (s_error$, r_name, 'ARRAY SIZES ARE NOT THE SAME')
      return
    endif
    do i = 1, size(l_var)
      l_var(i)%l = l_set(i)%l
    enddo
  endif

! A string set
! If value_str has "|" then it must be a datum array

elseif (size(s_var) /= 0) then
  if (index(var_str, '|merit_type') /= 0) then
    if (index(value_str, '|') == 0) then
      if (all (value_str /= tao_var_merit_type_name)) then
        call out_io (s_error$, r_name, 'BAD VARIABLE MERIT_TYPE NAME:' // value_str)
        return
      endif
      do i = 1, size(s_var)
        s_var(i)%s = value_str
      enddo

    else
      call tao_find_var (err, value_str, str_array=s_set)
      if (size(l_set) /= size(l_var)) then
        call out_io (s_error$, r_name, 'ARRAY SIZES ARE NOT THE SAME')
        return
      endif
      do i = 1, size(s_var)
        s_var(i)%s = s_set(i)%s
      enddo
    endif
  endif

! Only possibility left is real. The value_str might be a number or it might 
! be a mathematical expression involving datum values or array of values.

elseif (size(r_var) /= 0) then
  call tao_evaluate_expression (value_str, size(r_var),  .false., r_value, err, dflt_source = 'var')
  if (err) then
    call out_io (s_error$, r_name, 'BAD SET VALUE ' // value_str)
    return
  endif

  do i = 1, size(r_var)
    if (component == 'model') then
      call tao_set_var_model_value (v_var(i)%v, r_value(i))
    else
      r_var(i)%r = r_value(i)
    endif
  enddo

! Else must be an error

else
  call out_io (s_error$, r_name, 'NOTHING TO SET!')
  return

endif

call tao_set_var_useit_opt()
call tao_setup_key_table ()

end subroutine tao_set_var_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_branch_cmd (branch_str, component_str, value_str)
!
! Routine to set lattice branch values.
!
! Input:
!   branch_str      -- character(*): Which branch to set.
!   component_str   -- character(*): Which branch parameter to set.
!   value_str       -- character(*): What value to set it to.
!-

subroutine tao_set_branch_cmd (branch_str, component_str, value_str)

implicit none

integer i
logical, allocatable :: this_u(:)
logical err

character(*) branch_str, component_str, value_str
character(*), parameter :: r_name = 'tao_set_branch_cmd'
character(40) b_str

!

call tao_pick_universe (branch_str, b_str, this_u, err)
if (err) return

do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. this_u(i)) cycle
  call set_this_branch(s%u(i), err)
  s%u(i)%calc%lattice = .true.
  if (err) return
enddo

!--------------------------------------------
contains

subroutine set_this_branch(u, err)

type (tao_universe_struct), target :: u
type (branch_struct), pointer :: branch
integer ix
logical err
character(40) c_str

!

err = .true.

branch => pointer_to_branch(b_str, u%model%lat)
if (.not. associated(branch)) then
  call out_io (s_error$, r_name, 'BAD BRANCH NAME OR INDEX: ' // b_str)
  return
endif

!

call match_word (component_str, [character(28):: 'particle', 'default_tracking_species', 'geometry', 'live_branch'], &
                                                                                                    ix, matched_name = c_str)
if (ix < 1) then
  call out_io (s_error$, r_name, 'BAD BRANCH COMPONENT NAME: ' // component_str)
  return
endif

select case (c_str)
case ('particle')
  ix = species_id(value_str)
  if (ix == invalid$ .or. ix == ref_particle$ .or. ix == anti_ref_particle$) then
    call out_io (s_error$, r_name, 'INVALID REFERENCE PARTICLE SPECIES: ' // value_str)
    return
  endif
  branch%param%particle = ix

case ('default_tracking_species')
  ix = species_id(value_str)
  if (ix == invalid$) then
    call out_io (s_error$, r_name, 'INVALID DEFAULT TRACKING SPECIES: ' // value_str)
    return
  endif
  branch%param%default_tracking_species = ix

case ('geometry')
  call tao_set_switch_value (ix, c_str, value_str, geometry_name(1:), 1, err)
  if (err) return
  branch%param%geometry = ix
  if (ix == open$) u%model%lat%particle_start = u%model%tao_branch(branch%ix_branch)%orbit(0)
  if (ix == closed$) s%com%force_chrom_calc = .true.

case ('live_branch')
  call tao_set_logical_value (branch%param%live_branch, c_str, value_str, err)
  if (err) return

end select

err = .false.

end subroutine set_this_branch

end subroutine tao_set_branch_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_data_cmd (who_str, value_str, silent)
!
! Routine to set data values.
!
! Input:
!   who_str   -- Character(*): Which data component(s) to set.
!   value_str -- Character(*): What value to set it to.
!-

subroutine tao_set_data_cmd (who_str, value_str, silent)

implicit none

type (tao_data_array_struct), allocatable    :: d_dat(:)
type (tao_real_pointer_struct), allocatable  :: r_dat(:)
type (tao_integer_array_struct), allocatable :: int_dat(:), int_value(:)
type (tao_logical_array_struct), allocatable :: l_dat(:), l_value(:)
type (tao_string_array_struct), allocatable :: s_dat(:), s_value(:)
type (tao_universe_struct), pointer :: u
type (tao_data_struct), pointer :: d
type (branch_struct), pointer :: branch
type (ele_pointer_struct), allocatable :: eles(:)

real(rp), allocatable :: r_value(:)

integer i, ix, i1, n_loc, ib, ie
integer, allocatable :: int_save(:)

character(*) who_str, value_str
character(20) component
character(*), parameter :: r_name = 'tao_set_data_cmd'
character(100) :: why_invalid, tmpstr
character, allocatable :: s_save(:)

logical, optional :: silent
logical err, l1


! Decode data component to set.

call tao_find_data (err, who_str, d_array = d_dat, re_array=r_dat, &
          log_array=l_dat, str_array = s_dat, int_array = int_dat, component = component)
if (err) return

select case (component)
  case ('model', 'base', 'design', 'old', 'model_value', 'base_value', 'design_value', &
             'old_value', 'invalid', 'delta_merit', 'merit', 'exists', 'good_base ', &
             'useit_opt ', 'useit_plot', 'ix_d1', 'good_model', 'good_design')
  call out_io (s_error$, r_name, 'DATUM ATTRIBUTE NOT SETTABLE: ' // component)
  return
end select

!------------------------
! A logical value_str is either a logical or an array of datum values.

if (size(l_dat) /= 0) then
  if (is_logical(value_str)) then
    read (value_str, *) l1
    do i = 1, size(l_dat)
      l_dat(i)%l = l1
    enddo

  else
    call tao_find_data (err, value_str, log_array=l_value)
    if (size(l_value) /= size(l_dat) .and. size(l_value) /= 1) then
      call out_io (s_error$, r_name, 'ARRAY SIZES ARE NOT THE SAME')
      return
    endif
    do i = 1, size(l_dat)
      if (size(l_value) == 1) then
        l_dat(i)%l = l_value(1)%l
      else
        l_dat(i)%l = l_value(i)%l
      endif
    enddo
  endif

!------------------------
! An integer value_str is either an integer or an array of datum values.

elseif (size(int_dat) /= 0) then

  allocate (int_save(size(int_dat)))  ! Used to save old values in case of error

  if (is_integer(value_str)) then
    read (value_str, *) i1
    do i = 1, size(int_dat)
      int_save(i) = int_dat(i)%i
      int_dat(i)%i = i1
    enddo

  elseif (component == 'eval_point' .and. index(value_str, '|') == 0) then
    call match_word (value_str, anchor_pt_name(1:), i1, can_abbreviate = .false.)
    if (i1 == 0) then
      call out_io (s_error$, r_name, 'eval_point setting is "beginning", "center", or "end".')
      return
    endif
    do i = 1, size(int_dat)
      int_save(i) = int_dat(i)%i
      int_dat(i)%i = i1
    enddo

  else
    call tao_find_data (err, value_str, int_array=int_value)
    if (size(int_value) /= size(int_dat) .and. size(int_dat) /= 1) then
      call out_io (s_error$, r_name, 'ARRAY SIZES ARE NOT THE SAME')
      return
    endif
    do i = 1, size(int_dat)
      if (size(int_dat) == 1) then
        int_save(i) = int_dat(1)%i
        int_dat(i)%i = int_value(1)%i
      else
        int_save(i) = int_dat(i)%i
        int_dat(i)%i = int_value(i)%i
      endif
    enddo
  endif

  if (component == 'ix_ele' .or. component == 'ix_ele_start' .or. component == 'ix_ele_ref') then
    do i = 1, size(int_dat)
      u => s%u(d_dat(i)%d%d1%d2%ix_universe)
      branch => u%design%lat%branch(d_dat(i)%d%ix_branch)
      ie = int_dat(i)%i
      if (ie < 0) ie = -1
      tmpstr = ''

      if (ie > branch%n_ele_max) then
        int_dat(i)%i = int_save(i)
        call out_io (s_error$, r_name, 'ELEMENT INDEX OUT OF RANGE.')
        return
      endif

      if (component == 'ix_ele') then
        if (ie > -1) tmpstr = branch%ele(ie)%name
        d_dat(i)%d%ele_name = upcase(tmpstr)   ! Use temp due to bug on Windows
        d_dat(i)%d%ix_ele = ie
      elseif (component == 'ix_ele_start') then
        if (ie > -1) tmpstr = branch%ele(ie)%name
        d_dat(i)%d%ele_start_name = upcase(tmpstr)   ! Use temp due to bug on Windows
        d_dat(i)%d%ix_ele_start = ie
      else
        if (ie > -1) tmpstr = branch%ele(ie)%name
        d_dat(i)%d%ele_ref_name = upcase(tmpstr)   ! Use temp due to bug on Windows
        d_dat(i)%d%ix_ele_ref = ie
      endif
    enddo

  elseif (component == 'ix_branch') then
    do i = 1, size(int_dat)
      u => s%u(d_dat(i)%d%d1%d2%ix_universe)
      ib = int_dat(i)%i
      if (ib < 0 .or. ib > ubound(u%design%lat%branch, 1)) then
        int_dat(i)%i = int_save(i)
        call out_io (s_error$, r_name, 'ELEMENT INDEX OUT OF RANGE.')
        return
      endif
      d_dat(i)%d%ele_name = ''
      d_dat(i)%d%ix_ele = -1
      d_dat(i)%d%ele_ref_name = ''
      d_dat(i)%d%ix_ele_ref = -1
      d_dat(i)%d%ele_start_name = ''
      d_dat(i)%d%ix_ele_start = -1
    enddo
  endif

!------------------------
! A string:

elseif (size(s_dat) /= 0) then

  allocate (s_save(size(s_dat)))  ! Used to save old values in case of error

  ! If value_string has "|" then it must be a datum array

  if (index(value_str, '|') /= 0) then
    call tao_find_data (err, value_str, str_array=s_value)
    if (size(s_value) /= size(s_dat) .and. size(s_value) /= 1) then
      call out_io (s_error$, r_name, 'ARRAY SIZES ARE NOT THE SAME')
      return
    endif

    do i = 1, size(s_dat)
      tmpstr = s_value(i)%s
      s_save(i) = tmpstr
      s_dat(i)%s = tmpstr   ! Use temp due to bug on Windows
    enddo

  else
    if (component == 'merit_type' .and. all(value_str /= tao_data_merit_type_name)) then
      call out_io (s_error$, r_name, 'BAD DATA MERIT_TYPE NAME:' // value_str)
      return
    endif

    do i = 1, size(s_dat)
      s_save(i) = value_str
      ! Setting s_data(i)%s for %data_type does not work since %data_type is a var length string.
      select case (component)
      case ('data_type');       d_dat(i)%d%data_type = value_str
      case default;             s_dat(i)%s = value_str
      end select
    enddo
  endif

  !

  if (component == 'ele_name' .or. component == 'ele_start_name' .or. component == 'ele_ref_name') then
    do i = 1, size(d_dat)
      u => s%u(d_dat(i)%d%d1%d2%ix_universe)
      call upcase_string (s_dat(i)%s)
      if (s_dat(i)%s /= '') then
        call lat_ele_locator (s_dat(i)%s, u%design%lat, eles, n_loc)
        if (n_loc == 0) then
          call out_io (s_error$, r_name, 'ELEMENT NOT LOCATED: ' // s_dat(i)%s)
          s_dat(i)%s = s_save(i)
          return
        endif
        if (n_loc > 1) then
          call out_io (s_error$, r_name, 'MULTIPLE ELEMENT OF THE SAME NAME EXIST: ' // s_dat(i)%s)
          s_dat(i)%s = s_save(i)
          return
        endif
      endif

      if (component == 'ele_name') then
        if (s_dat(i)%s == '') then
          d_dat(i)%d%ix_ele = -1
          d_dat(i)%d%ele_name = ''
        else
          d_dat(i)%d%ix_ele    = eles(1)%ele%ix_ele
          d_dat(i)%d%ele_name  = eles(1)%ele%name
          if (d_dat(i)%d%ix_branch /= eles(1)%ele%ix_branch) then
            d_dat(i)%d%ele_ref_name = ''
            d_dat(i)%d%ix_ele_ref = -1
            d_dat(i)%d%ele_start_name = ''
            d_dat(i)%d%ix_ele_start = -1
          endif
          d_dat(i)%d%ix_branch = eles(1)%ele%ix_branch
        endif

      elseif (component == 'ele_start_name') then
        if (s_dat(i)%s == '') then
          d_dat(i)%d%ix_ele_start = -1
          d_dat(i)%d%ele_start_name = ''
        else
          if (d_dat(i)%d%ix_branch /= eles(1)%ele%ix_branch) then
            s_dat(i)%s = s_save(i)
            call out_io (s_error$, r_name, 'START_ELEMENT IS IN DIFFERENT BRANCH FROM ELEMENT.')
            return
          endif
          d_dat(i)%d%ix_ele_start   = eles(1)%ele%ix_ele
          d_dat(i)%d%ele_start_name = eles(1)%ele%name
        endif

      else
        if (s_dat(i)%s == '') then
          d_dat(i)%d%ix_ele_ref = -1
          d_dat(i)%d%ele_ref_name = ''
        else
          if (d_dat(i)%d%ix_branch /= eles(1)%ele%ix_branch) then
            s_dat(i)%s = s_save(i)
            call out_io (s_error$, r_name, 'REF_ELEMENT IS IN DIFFERENT BRANCH FROM ELEMENT.')
            return
          endif
          d_dat(i)%d%ix_ele_ref   = eles(1)%ele%ix_ele
          d_dat(i)%d%ele_ref_name = eles(1)%ele%name
        endif
      endif
    enddo
  endif

!------------------------
! Only possibility left is real. The value_str might be a number or it might 
! be a mathematical expression involving datum values or array of values.

elseif (size(r_dat) /= 0) then
  call tao_evaluate_expression (value_str, size(r_dat), .false., r_value, err, dflt_source = 'data')
  if (err) then
    call out_io (s_error$, r_name, 'BAD SET VALUE ' // value_str)
    return
  endif

  do i = 1, size(r_dat)
    r_dat(i)%r = r_value(i)
    if (component == 'meas') d_dat(i)%d%good_meas = .true.
    if (component == 'ref')  d_dat(i)%d%good_ref = .true.
    if (component == 'base') d_dat(i)%d%good_base = .true.
  enddo

else
  call out_io (s_error$, r_name, 'LEFT HAND SIDE MUST POINT TO A SCALAR OR ARRAY OF DATA COMPONENTS.')
  return
endif

!----------------------
! If the "exists" component has been set (used by gui interface) check if the datum is truely valid.

select case (component)
case ('ele_name', 'ele_start_name', 'ele_ref_name', 'data_type', 'data_source', 'ix_uni', &
      'ix_branch', 'ix_ele', 'ix_ele_start', 'ix_ele_ref', 'eval_point', 'exists')
  do i = 1, size(d_dat)
    d => d_dat(i)%d
    if (component == 'exists' .and. .not. d%exists) cycle

    d%exists = tao_data_sanity_check(d, .not. logic_option(.false., silent), '')
    if (.not. d%exists) cycle

    u => s%u(d%d1%d2%ix_universe) 
    call tao_evaluate_a_datum (d, u, u%model, d%model_value, d%good_model, why_invalid)
    if (.not. d%good_model) then
      call out_io (s_error$, r_name, 'Datum is not valid since: ' // why_invalid)
    endif
    if (d%good_model) call tao_evaluate_a_datum (d, u, u%design, d%design_value, d%good_design, why_invalid)
  enddo
end select

! End stuff

call tao_set_data_useit_opt()

end subroutine tao_set_data_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_default_cmd (who_str, value_str)
!
! Routine to set default values.
!
! Input:
!   who_str   -- Character(*): Which default component(s) to set.
!   value_str -- Character(*): What value to set it to.
!
!  Output:
!-

subroutine tao_set_default_cmd (who_str, value_str)

implicit none

type (tao_universe_struct), pointer :: u
integer ix, iu
logical err

character(*) who_str, value_str
character(16) switch
character(*), parameter :: r_name = 'tao_set_default_cmd'
!

call match_word (who_str, ['universe', 'branch  '], ix, matched_name=switch)
if (ix < 1) then
  call out_io (s_error$, r_name, 'BAD DEFAULT NAME: ' // who_str)
  return
endif

select case (switch)
case ('universe')
  call tao_set_integer_value (s%global%default_universe, 'UNIVERSE', value_str, err, lbound(s%u, 1), ubound(s%u, 1))
  if (err) return
  call tao_turn_on_special_calcs_if_needed_for_plotting()
  u => tao_pointer_to_universe(-1)
  if (s%global%default_branch > ubound(u%model%lat%branch, 1)) then
    call out_io (s_error$, r_name, 'DEFAULT_BRANCH VALUE NOW OUT OF RANGE: ' // int_str(s%global%default_branch), &
                                   'SETTING TO ZERO')
    s%global%default_branch = 0
  endif

case ('branch')
  u => tao_pointer_to_universe(-1)
  call tao_set_integer_value (s%global%default_branch, 'BRANCH', value_str, err, 0, ubound(u%model%lat%branch, 1))
  if (err) return

end select

end subroutine tao_set_default_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_dynamic_aperture_cmd (who, value_str)
!
! Sets dynamic aperture parameters.
!
! Input:
!   who       -- Character(*): which parameter to set.
!   value_str -- Character(*): Value to set to.
!-

subroutine tao_set_dynamic_aperture_cmd (who, value_str)

type (tao_universe_struct), pointer :: u
type (aperture_param_struct) ap_param

real(rp) pz(20), a_emit, b_emit, ellipse_scale
integer i, iu, ios, iof, n, n2

character(*) who, value_str
character(20) who2
character(*), parameter :: r_name = 'tao_set_dynamic_aperture_cmd'
logical, allocatable :: this_u(:)
logical err

namelist / params / ap_param, pz, a_emit, b_emit, ellipse_scale

!

call tao_pick_universe (unquote(who), who2, this_u, err); if (err) return

iof = tao_open_scratch_file (err);  if (err) return

write (iof, '(a)') '&params'
if (who2(1:2) == 'pz') then
  write (iof, '(a)') trim(who2) // ' = ' // trim(value_str)
else
  write (iof, '(a)') ' ap_param%' // trim(who2) // ' = ' // trim(value_str)
endif
write (iof, '(a)') '/'

do iu = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. this_u(iu)) cycle
  u => s%u(iu)

  n = size(u%dynamic_aperture%pz)
  pz = -2
  pz(1:n) = u%dynamic_aperture%pz
  a_emit = u%dynamic_aperture%a_emit
  b_emit = u%dynamic_aperture%b_emit
  ellipse_scale = u%dynamic_aperture%ellipse_scale
  ap_param = u%dynamic_aperture%param
  rewind(iof)
  read (iof, nml = params, iostat = ios)

  if (ios /= 0) then
    call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
    close (iof, status = 'delete') 
    return
  endif

  u%dynamic_aperture%param = ap_param
  u%dynamic_aperture%a_emit = a_emit
  u%dynamic_aperture%b_emit = b_emit
  u%dynamic_aperture%ellipse_scale = ellipse_scale

  n = 0
  do i = 1, size(pz)
    if (pz(i) <= -1) cycle
    n = n + 1
    pz(n) = pz(i)
  enddo

  call re_allocate(u%dynamic_aperture%pz, n)
  u%dynamic_aperture%pz = pz(1:n)

  if (allocated(u%dynamic_aperture%scan)) deallocate (u%dynamic_aperture%scan)
  allocate(u%dynamic_aperture%scan(n))
enddo

close (iof, status = 'delete') 

end subroutine tao_set_dynamic_aperture_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_universe_cmd (uni, who, what)
!
! Sets a universe on or off, or sets the recalculate or twiss_calc logicals, etc.
!
! Input:
!   uni     -- Character(*): which universe; 0 => current viewed universe
!   who     -- Character(*): "on", "off", "recalculate", "dynamic_aperture_calc", "one_turn_map_calc", or "twiss_calc"
!   what    -- Character(*): "on" or "off" for who = "dynamic_aperture_calc", "one_turn_map_calc" or "twiss_calc".
!-

subroutine tao_set_universe_cmd (uni, who, what)

implicit none

type (tao_universe_struct), pointer :: u

integer i, iu

character(*) uni, who, what
character(*), parameter :: r_name = "tao_set_universe_cmd"
character(40) str

logical err, mat6_toggle, calc_ok
logical, allocatable :: picked(:)


! Pick universe

call downcase_string(what)

call tao_pick_universe(uni, str, picked, err, dflt_uni = -1, pure_uni = .true.)
if (err) return
if (str /= '') then
  call out_io (s_error$, r_name, 'MALFORMED "SET UNIVERSE" COMMAND.')
endif

s%u%picked_uni = .false.  ! Used by recalculate command.

do iu = 1, ubound(s%u, 1)
  if (.not. picked(iu)) cycle
  u => s%u(iu)

  ! Twiss calc.
  ! "mat6_recalc" is old style

  if (index('twiss_calc', trim(who)) == 1 .or. index('mat6_recalc', trim(who)) == 1) then
    if (what == 'on') then
      u%calc%twiss = .true.
      u%calc%lattice = .true.
    elseif (what == 'off') then
      u%calc%twiss = .false.
    else
      call out_io (s_error$, r_name, 'Syntax is: "set universe <uni_num> twiss_calc on/off"')
      return
    endif
    cycle
  endif

  ! Track calc
  ! "track_recalc" is old style.

  if (index('track_calc', trim(who)) == 1 .or. index('track_recalc', trim(who)) == 1) then
    if (what == 'on') then
      u%calc%track = .true.
      u%calc%lattice = .true.
    elseif (what == 'off') then
      u%calc%track = .false.
    else
      call out_io (s_error$, r_name, 'Syntax is: "set universe <uni_num> track_calc on/off"')
      return
    endif

    cycle
  endif

  ! Dynamic aperture calc.

  if (index('dynamic_aperture_calc', trim(who)) == 1) then
    if (what == 'on') then
      u%calc%dynamic_aperture = .true.
      u%calc%lattice = .true.
    elseif (what == 'off') then
      u%calc%dynamic_aperture = .false.
    else
      call out_io (s_error$, r_name, 'Syntax is: "set universe <uni_num> dynamic_aperture_calc on/off"')
      return
    endif

    cycle
  endif  

  ! One turn map calc.

  if ('one_turn_map_calc' == trim(who)) then
    if (what == 'on' .or. index('true', trim(what)) == 1) then
      u%calc%one_turn_map = .true.
      u%calc%lattice = .true.
    elseif (what == 'off' .or. index('false', trim(what)) == 1) then
      u%calc%one_turn_map = .false.
    else
      call out_io (s_error$, r_name, 'Syntax is: "set universe <uni_num> one_turn_map_calc on/off"')
      return
    endif

    cycle
  endif  
    
  ! Recalc.
  ! If universe is off then turn on and mark it to be turned off after a recalc.

  if (what /= '') then
    call out_io (s_error$, r_name, 'Extra stuff on line. Nothing done.')
    return
  endif

  if (index('recalculate', trim(who)) == 1) then
    u%calc%lattice = .true.
    if (.not. u%is_on) then
      u%is_on = .true.
      u%picked_uni = .true.
    endif
    cycle
  endif

  !

  if (who == 'on') then
    u%is_on = .true.
  elseif (who == 'off') then
    u%is_on = .false.
  else
    call out_io (s_error$, r_name, 'Choices are: "on", "off", "recalculate", "track_recalc", "twiss_calc", etc.')
    return
  endif

enddo

!

call tao_set_data_useit_opt()
call tao_lattice_calc (calc_ok)
  
do iu = 1, ubound(s%u, 1)
  u => s%u(iu)
  if (u%picked_uni) u%is_on = .false.
enddo

end subroutine tao_set_universe_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_elements_cmd (ele_list, attribute, value, update, lord_set)
!
! Sets element parameters.
!
! Input:
!   ele_list   -- Character(*): which elements.
!   attribute  -- Character(*): Attribute to set.
!   value      -- Character(*): Value to set.
!-

subroutine tao_set_elements_cmd (ele_list, attribute, value, update, lord_set)

use attribute_mod, only: attribute_type

implicit none

type (ele_pointer_struct), allocatable :: eles(:), v_eles(:)
type (tao_universe_struct), pointer :: u
type (all_pointer_struct) a_ptr
type (tao_lattice_struct), pointer :: tao_lat

real(rp), allocatable :: set_val(:)
integer i, ix, ix2, j, n_uni, n_set, n_eles, lat_type

character(*) ele_list, attribute, value
character(*), parameter :: r_name = "tao_set_elements_cmd"
character(100) val_str

logical update, lord_set
logical is_on, err, mat6_toggle

! Find elements

call tao_locate_all_elements (ele_list, eles, err)
if (err) return
if (size(eles) == 0) then
  call out_io (s_error$, r_name, 'CANNOT FIND ANY ELEMENTS CORRESPONDING TO: ' // ele_list)
  return
endif

!-----------
! The first complication is that what is being set can be a logical, switch (like an element's tracking_method), etc.
! So use set_ele_attribute to do the set.
! But set_ele_attribute does not know about Tao syntax so it may have problems evaluateing the value string.
! And set_ele_attribute cannot handle the situation where there is an array of set values.
! How to handle this depends upon what type of attribute it is.

! If a real attribute then use tao_evaluate_expression to evaluate.
! If attribute_type returns invalid_name$ then assume attribute is a controller variable which are always real.

if (attribute_type(upcase(attribute)) == is_real$ .or. attribute_type(upcase(attribute)) == invalid_name$) then
  ! Important to use "size(eles)" as 2nd arg instead of "0" since if value is something like "ran()" then
  ! want a an array of set_val values with each value different.
  call tao_evaluate_expression (value, size(eles), .false., set_val, err)
  if (err) return

  if (size(eles) /= size(set_val)) then
    call out_io (s_error$, r_name, 'SIZE OF VALUE ARRAY NOT EQUAL TO THE SIZE OF THE ELEMENTS TO BE SET.', &
                                   'NOTHING DONE.')
    return
  endif

  n_set = 0
  do i = 1, size(eles)
    call pointer_to_attribute(eles(i)%ele, attribute, .true., a_ptr, err, err_print_flag = .false.)
    if (err) cycle
    if (.not. associated(a_ptr%r)) then
      call out_io (s_error$, r_name, 'STRANGE ERROR: PLEASE CONTACT HELP.')
      return
    endif

    call set_ele_real_attribute (eles(i)%ele, attribute, set_val(i), err, .false.)
    if (.not. err) n_set = n_set + 1
    call tao_set_flags_for_changed_attribute (s%u(eles(i)%id), eles(i)%ele%name, eles(i)%ele, a_ptr%r)
  enddo

  if (n_set == 0) then
    i = size(eles)
    call set_ele_real_attribute (eles(i)%ele, attribute, set_val(i), err, .true.)
    call out_io (s_error$, r_name, 'NOTHING SET.')
  endif

  do i = lbound(s%u, 1), ubound(s%u, 1)
    u => s%u(i)
    if (.not. u%calc%lattice .or. .not. s%global%lattice_calc_on) cycle
    call lattice_bookkeeper (u%model%lat)
  enddo

  call tao_var_check(eles, attribute, update)

  return

! If there is a "ele::" construct in the value string...

elseif (index(value, 'ele::') /= 0) then

  val_str = value
  u => tao_pointer_to_universe(val_str)

  if (val_str(1:5) /= 'ele::') then
    call out_io (s_error$, r_name, 'CANNOT PARSE SET VALUE: ' // value, &
                                   'PLEASE CONTACT BMAD MAINTAINER DAVID SAGAN.')
    return
  endif
  val_str = val_str(6:)

  lat_type = model$
  if (index(val_str, '|model') /= 0) then
    lat_type = model$
  elseif (index(val_str, '|design') /= 0) then
    lat_type = design$
  elseif (index(val_str, '|base') /= 0) then
    lat_type = base$
  elseif (index(val_str, '|') /= 0) then
    call out_io (s_error$, r_name, 'BAD SET VALUE: ' // value)
    return
  endif

  ix = index(val_str, '|')
  if (ix /= 0) val_str = val_str(:ix-1)
  tao_lat => tao_pointer_to_tao_lat(u, lat_type)

  ix = index(val_str, '[')
  if (ix == 0) then
    call out_io (s_error$, r_name, 'BAD SET VALUE: ' // value)
    return
  endif

  call lat_ele_locator (val_str(:ix-1), tao_lat%lat, v_eles, n_eles, err); if (err) return
  ix2 = len_trim(val_str)
  call pointer_to_attribute (v_eles(1)%ele, val_str(ix+1:ix2-1), .false., a_ptr, err); if (err) return
  val_str = all_pointer_to_string (a_ptr, err = err)
  if (err) then
    call out_io (s_error$, r_name, 'STRANGE SET VALUE: ' // value)
    return
  endif

! If the value string does not have "ele::" then just 
! assume that set_ele_attribute will be able to evaluate the value string.

else
  val_str = value
endif

! When a wild card is used so there are multiple elements involved, an error
! generated by some, but not all elements is not considered a true error.
! For example: "set ele * csr_calc_on = t" is not valid for markers.

n_set = 0
do i = 1, size(eles)
  u => s%u(eles(i)%id)
  call set_ele_attribute (eles(i)%ele, trim(attribute) // '=' // trim(val_str), err, .false., lord_set)
  call tao_set_flags_for_changed_attribute (u, eles(i)%ele%name, eles(i)%ele)
  if (.not. err) n_set = n_set + 1
enddo

if (n_set == 0) then
  if (lord_set) then
    call out_io (s_error$, r_name, &
          'NOTHING SET. REMEMBER: "-lord_no_set" USAGE WILL PREVENT LORD ELEMENT ATTRIBUTES FROM BEING SET.')
  else
    call out_io (s_error$, r_name, 'NOTHING SET. EG:')
    u => s%u(eles(1)%id)
    call set_ele_attribute (eles(1)%ele, trim(attribute) // '=' // trim(val_str),  err, .true., lord_set)
  endif
  return
endif

! End stuff

if (n_set /= size(eles)) then
  call out_io (s_info$, r_name, 'Note: \i0\ elements (out of \i0\) set.', i_array = [n_set, size(eles)])
endif

do i = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(i)
  if (.not. u%calc%lattice .or. .not. s%global%lattice_calc_on) cycle
  call lattice_bookkeeper (u%model%lat)
enddo

call tao_var_check(eles, attribute, update)

end subroutine tao_set_elements_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_logical_value (var, var_str, value_str, error)
!
! Subroutine to read and set the value of an logical varialbe.
!
! If the value is out of the range [min_val, max_val] then an error message will
! be generated and the variable will not be set.
!
! Input:
!   var_str   -- Character(*): Used for error messages.
!   value_str -- Character(*): String with encoded value.
!
! Output:
!   var   -- Logical: Variable to set.
!   error -- Logical: Set True on an error. False otherwise.
!-

subroutine tao_set_logical_value (var, var_str, value_str, error)

implicit none

logical var, ix
integer ios

character(*) var_str, value_str
character(*), parameter :: r_name = 'tao_set_logical_value'
logical error

!

error = .true.
read (value_str, '(l)', iostat = ios) ix

if (ios /= 0 .or. len_trim(value_str) == 0) then
  call out_io (s_error$, r_name, 'BAD ' // trim(var_str) // ' VALUE.')
  return
endif

var = ix      
error = .false.

end subroutine tao_set_logical_value 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_integer_value (var, var_str, value_str, error, min_val, max_val, print_err)
!
! Subroutine to read and set the value of an integer varialbe.
!
! If the value is out of the range [min_val, max_val] then an error message will
! be generated and the variable will not be set.
!
! Input:
!   var_str   -- Character(*): Used for error messages.
!   value_str -- Character(*): String with encoded value.
!   min_val   -- Integer, optional: Minimum value. 
!   max_val   -- Integer, optional: Maximum value.
!   print_err -- logical, optional: If True, print error message. Default is true
!
! Output:
!   var   -- Integer: Variable to set.
!   error -- Logical: Set True on an error. False otherwise.
!-

subroutine tao_set_integer_value (var, var_str, value_str, error, min_val, max_val, print_err)

implicit none

integer var
integer, optional :: min_val, max_val
integer ios, ix

character(*) var_str, value_str
character(*), parameter :: r_name = 'tao_set_integer_value'
logical error
logical, optional :: print_err

!

error = .true.
read (value_str, *, iostat = ios) ix

if (ios /= 0 .or. len_trim(value_str) == 0) then
  if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'BAD ' // trim(var_str) // ' VALUE.')
  return
endif

if (present(min_val)) then
  if (ix < min_val) then
  if (logic_option(.true., print_err)) call out_io (s_error$, r_name, trim(var_str) // ' VALUE TOO SMALL.')
    return 
  endif
endif

if (present(max_val)) then
  if (ix > max_val) then
  if (logic_option(.true., print_err)) call out_io (s_error$, r_name, trim(var_str) // ' VALUE TOO LARGE.')
    return 
  endif
endif

var = ix
error = .false.

end subroutine tao_set_integer_value

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_switch_value (switch_val, err_str, value_str, name_list, l_bound, error, switch_str)
!
! Routine to set the value of an integer switch.
!
! If the value is out of the range [min_val, max_val] then an error message will
! be generated and the switch will not be set.
!
! Input:
!   err_str       -- character(*): Used for error messages.
!   value_str     -- character(*): String with encoded value.
!   name_list(:)  -- character(*): Names to match to.
!   l_bound       -- integer: Lower bound to name_list(:) array.
!
! Output:
!   switch_val  -- integer: Parameter to set. Not set if there is an error.
!   error       -- logical: Set True on an error. False otherwise.
!   switch_str  -- character(*), optional: Set to the string representation of switch_val. 
!                   Not set if there is an error.
!-

subroutine tao_set_switch_value (switch_val, err_str, value_str, name_list, l_bound, error, switch_str)

implicit none

integer switch_val, l_bound
integer ios, ix

character(*) err_str, value_str
character(*) name_list(l_bound:)
character(*), optional :: switch_str
character(*), parameter :: r_name = 'tao_set_switch_value'
logical error

!

error = .true.

call match_word(unquote(value_str), name_list, ix, .false., .true.)

if (ix == 0) then
  call out_io (s_error$, r_name, trim(err_str) // ' VALUE IS UNKNOWN.')
  return 
endif

if (ix < 0) then
  call out_io (s_error$, r_name, trim(err_str) // ' ABBREVIATION MATCHES MULTIPLE NAMES.')
  return 
endif

switch_val = ix + (l_bound - 1)
if (present(switch_str)) switch_str = downcase(name_list(switch_val))

error = .false.

end subroutine tao_set_switch_value

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_real_value (var, var_str, value_str, error, min_val, max_val, dflt_uni)
!
! Subroutine to read and set the value of a real variable.
!
! If the value is out of the range [min_val, max_val] then an error message will
! be generated and the variable will not be set.
!
! Input:
!   var_str   -- Character(*): Used for error messages.
!   value_str -- Character(*): String with encoded value.
!   min_val   -- real(rp), optional: Minimum value. 
!   max_val   -- real(rp), optional: Maximum value.
!   dflt_uni  -- integer, optional: Default universe used to evaluate parameters.
!
! Output:
!   var   -- real(rp): Variable to set.
!   error -- Logical: Set True on an error. False otherwise.
!-

subroutine tao_set_real_value (var, var_str, value_str, error, min_val, max_val, dflt_uni)

implicit none

real(rp) var, var_value
real(rp), allocatable :: var_array(:)
real(rp), optional :: min_val, max_val
integer, optional :: dflt_uni
integer ios

character(*) var_str, value_str
character(*), parameter :: r_name = 'tao_set_real_value'
logical error

!

call tao_evaluate_expression (value_str, 1, .false., var_array, error, .true., dflt_uni = dflt_uni)
if (error) return

var_value = var_array(1)
error = .true.

if (present(min_val)) then
  if (var_value < min_val) then
    call out_io (s_error$, r_name, trim(var_str) // ' VALUE OUT OF RANGE.')
    return
  endif
endif

if (present(max_val)) then
  if (var_value > max_val) then
    call out_io (s_error$, r_name, trim(var_str) // ' VALUE OUT OF RANGE.')
    return
  endif
endif

var = var_value
error = .false.

end subroutine tao_set_real_value

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_drawing_cmd (drawing, component, value_str)
!
! Routine to set floor_plan and lat_layout parameters.
! 
! Input:
!   drawing   -- tao_drawing_struct: s%plot_page%floor_plan or s%plot_page%lat_layout.
!   component -- Character(*): Which shape component to set.
!   value_str -- Character(*): Value to set to.
!
! Output:
!    s%shape  -- Shape variables structure.
!-

subroutine tao_set_drawing_cmd (drawing, component, value_str)

use tao_input_struct

implicit none

type (tao_ele_shape_input), target :: ele_shape(50)
type (tao_drawing_struct), target :: drawing
type (tao_ele_shape_input), pointer :: es_in
type (tao_ele_shape_struct), pointer :: es

character(*) component, value_str
character(*), parameter :: r_name = 'tao_set_drawing_cmd'

integer i, ix, n, iu, ios

logical err, needs_quotes

namelist / params / ele_shape

! Init

n = size(drawing%ele_shape)
ele_shape = tao_ele_shape_input()
ele_shape(1:n) = tao_ele_shape_struct_to_input(drawing%ele_shape)

! Setup

needs_quotes = .false.
ix = index(component, '%')

if (ix /= 0) then
  select case (component(ix+1:))
  case ('ele_name', 'name')
    needs_quotes = .true.
    component = 'ele_id%' // component(ix+1:)
  case ('shape', 'color', 'label', 'ele_id')
    needs_quotes = .true.
  end select

  ! Something like 30*"XXX" does not need quotes
  n = len_trim(value_str)
  if (value_str(n:n) == "'" .or. value_str(n:n) == '"') needs_quotes = .false.
endif

! open a scratch file for a namelist read

iu = tao_open_scratch_file (err);  if (err) return

write (iu, '(a)') '&params'
if (needs_quotes) then
  write (iu, '(a)') trim(component) // ' = "' // trim(value_str) // '"'
else
  write (iu, '(a)') trim(component) // ' = ' // trim(value_str)
endif
write (iu, '(a)') '/'
write (iu, *)
rewind (iu)
read (iu, nml = params, iostat = ios)
close (iu, status = 'delete')

if (ios /= 0) then
  call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
  return
endif

! Transfer to drawing

do n = size(ele_shape), 1, -1
  if (ele_shape(n)%ele_id /= '') exit
enddo

if (n > size(drawing%ele_shape)) then
  deallocate(drawing%ele_shape)
  allocate (drawing%ele_shape(n))
endif

do i = 1, n
  drawing%ele_shape(i) = tao_ele_shape_input_to_struct(ele_shape(i))
  call tao_shape_init (drawing%ele_shape(i), err, .true.);  if (err) return
enddo

end subroutine tao_set_drawing_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_symbolic_number_cmd (sym_str, num_str, val)
!
! Associates a given symbol with a given number.
! Note: Either num_str or val argument must be present. 
!
! Input:
!   sym_str     -- character(*): Symbol.
!   num_str     -- character(*), optional: Symbol value expression.
!   val         -- real(rp), optional: Value of symbol
!-

subroutine tao_set_symbolic_number_cmd (sym_str, num_str, val)

type (named_number_struct), allocatable :: sym_temp(:)

integer i, n
real(rp), optional :: val
real(rp), allocatable :: value(:)
logical err, ok

character(*) sym_str
character(*), optional :: num_str
character(40) s_str
character(*), parameter :: r_name = 'tao_set_symbolic_number_cmd'

!

s_str = adjustl(sym_str)

do i = 1, size(physical_const_list)
  if (s_str == physical_const_list(i)%name) then
    call out_io (s_error$, r_name, 'NAME MATCHES NAME OF A PHYSICAL CONSTANT. SET IGNORED.')
    return
  endif
enddo

if (present(val)) then
  allocate(value(1))
  value(1) = val

else
  if (index(num_str, 'chrom') /= 0) then
    s%com%force_chrom_calc = .true.
    s%u%calc%lattice = .true.
    call tao_lattice_calc(ok)
  endif

  call tao_evaluate_expression (num_str, 1, .false., value, err); if (err) return
endif

!

if (allocated(s%com%symbolic_num)) then
  n = size(s%com%symbolic_num)
  do i = 1, n
    if (s_str == s%com%symbolic_num(i)%name) exit
  enddo

  if (i == n + 1) then
    call move_alloc (s%com%symbolic_num, sym_temp)
    allocate (s%com%symbolic_num(n+1))
    s%com%symbolic_num(1:n) = sym_temp
  endif

  s%com%symbolic_num(i)%name = s_str
  s%com%symbolic_num(i)%value = value(1)

else
  allocate (s%com%symbolic_num(1)) 
  s%com%symbolic_num(1)%name = s_str
  s%com%symbolic_num(1)%value = value(1)
endif

end subroutine tao_set_symbolic_number_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_qp_rect_struct (qp_rect_name, component, qp_rect, value, error, ix_uni)
!
! Routine to set qp_rect_names of a qp_rect_struct.
!
! Input:
!   qp_rect_name    -- character(*): qp_rect name. Used for error messages.
!   component       -- character(*): qp_rect component name.
!   qp_rect         -- qp_rect_struct: qp_rect_struct with component to modify
!   value           -- character(*): Component value.
!
! Output:
!   qp_rect         -- qp_rect_struct: qp_rect_struct with changed component value.
!   error           -- logical: Set true if there is an error. False otherwise.
!   ix_uni          -- integer, optional: Tao universe number in case the value depends upon
!                       a parameter of a particular universe.
!-

subroutine tao_set_qp_rect_struct (qp_rect_name, component, qp_rect, value, error, ix_uni)

type (qp_rect_struct) qp_rect

real(rp) :: values(4)
integer, optional :: ix_uni
integer ios
character(*) qp_rect_name, component, value
character(*), parameter :: r_name = 'tao_set_qp_rect_struct '
logical error

!

select case (component)
case ('x1');  call tao_set_real_value(qp_rect%x1, component, value, error, dflt_uni = ix_uni)
case ('x2');  call tao_set_real_value(qp_rect%x2, component, value, error, dflt_uni = ix_uni)
case ('y1');  call tao_set_real_value(qp_rect%y1, component, value, error, dflt_uni = ix_uni)
case ('y2');  call tao_set_real_value(qp_rect%y2, component, value, error, dflt_uni = ix_uni)

case ('')
  read (value, *, iostat = ios) values
  if (ios /= 0) then
    call out_io (s_error$, r_name, "Need four values: " // value)
    error = .true.
    return
  endif
  qp_rect%x1 = values(1)
  qp_rect%x2 = values(2)
  qp_rect%y1 = values(3)
  qp_rect%y2 = values(4)
  error = .false.

case default
  call out_io (s_error$, r_name, "BAD QP_RECT COMPONENT: " // component)
  error = .true.
  return
end select

end subroutine tao_set_qp_rect_struct

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_qp_axis_struct (qp_axis_name, component, qp_axis, value, error, ix_uni)
!
! Routine to set qp_axis_names of a qp_axis_struct.
!
! Input:
!   qp_axis_name    -- character(*): qp_axis name. Used for error messages.
!   component       -- character(*): qp_axis component name.
!   qp_axis         -- qp_axis_struct: qp_axis_struct with component to modify
!   value           -- character(*): Component value.
!
! Output:
!   qp_axis         -- qp_axis_struct: qp_axis_struct with changed component value.
!   error           -- logical: Set true if there is an error. False otherwise.
!   ix_uni          -- integer, optional: Tao universe number in case the value depends upon
!                       a parameter of a particular universe.
!-

subroutine tao_set_qp_axis_struct (qp_axis_name, component, qp_axis, value, error, ix_uni)

use quick_plot, only: qp_string_to_enum

type (qp_axis_struct) qp_axis
character(*) component, value, qp_axis_name
character(40) val
character(*), parameter :: r_name = 'tao_set_qp_axis_struct '
integer, optional :: ix_uni
integer indx
logical error

! Note: Setting %tick_min, %tick_max, and %dtick are not implemented since they are considered dependent parameters.

select case (component)
case ('bounds')
  val = unquote(upcase(value))
  select case (val)
  case ('ZERO_AT_END', 'ZERO_SYMMETRIC', 'GENERAL', 'EXACT')
    qp_axis%bounds = val
    error = .false.
  case default
    call out_io (s_error$, r_name, 'BAD AXIS BOUNDS NAME:' // value, 'POSSIBILITIES ARE: "ZERO_AT_END", "ZERO_SYMMETRIC", OR "GENERAL".')
    error = .true.
  end select

case ('min')
  call tao_set_real_value (qp_axis%min, qp_axis_name, value, error, dflt_uni = ix_uni)
case ('max')
  call tao_set_real_value (qp_axis%max, qp_axis_name, value, error, dflt_uni = ix_uni)
case ('number_offset')
  call tao_set_real_value (qp_axis%number_offset, qp_axis_name, value, error, dflt_uni = ix_uni)
case ('label_offset')
  call tao_set_real_value (qp_axis%label_offset, qp_axis_name, value, error, dflt_uni = ix_uni)
case ('major_tick_len')
  call tao_set_real_value (qp_axis%major_tick_len, qp_axis_name, value, error, dflt_uni = ix_uni)
case ('minor_tick_len')
  call tao_set_real_value (qp_axis%minor_tick_len, qp_axis_name, value, error, dflt_uni = ix_uni)

case ('label_color')
  indx = qp_string_to_enum(value, 'color', -1)
  if (indx < 1) then
    call out_io (s_error$, r_name, 'BAD COLOR NAME: ' // value)
    error = .true.
  else
    qp_axis%label_color = component
    error = .false.
  endif

case ('major_div')
  call tao_set_integer_value (qp_axis%major_div, qp_axis_name, value, error, 1)
  ! If %major_div_nominal is positive, setting %major_div does not make sense.
  ! So if %major_div_nominal is positive, set %major_div_nominal to the set value.
  if (.not. error .and. qp_axis%major_div_nominal > 0) qp_axis%major_div_nominal = qp_axis%major_div

case ('major_div_nominal')
  call tao_set_integer_value (qp_axis%major_div_nominal, qp_axis_name, value, error)
case ('minor_div')
  call tao_set_integer_value (qp_axis%minor_div, qp_axis_name, value, error, 0)
case ('minor_div_max')
  call tao_set_integer_value (qp_axis%minor_div_max, qp_axis_name, value, error, 1)
case ('places')
  call tao_set_integer_value (qp_axis%places, qp_axis_name, value, error)
case ('tick_side')
  call tao_set_integer_value (qp_axis%tick_side, qp_axis_name, value, error, -1, 1)
case ('number_side')
  call tao_set_integer_value (qp_axis%number_side, qp_axis_name, value, error, -1, 1)
case ('label')
  qp_axis%label = value
  error = .false.

case ('type')
  val = unquote(upcase(value))
  select case (val)
  case ('LOG', 'LINEAR')
    qp_axis%type = val
    error = .false.
  case default
    call out_io (s_error$, r_name, 'BAD AXIS TYPE NAME:' // value, 'POSSIBILITIES ARE: "LOG", OR "LINEAR".')
    error = .true.
  end select

case ('draw_label')
  call tao_set_logical_value (qp_axis%draw_label, qp_axis_name, value, error)

case ('draw_numbers')
  call tao_set_logical_value (qp_axis%draw_numbers, qp_axis_name, value, error)

case default
  call out_io (s_error$, r_name, "BAD QP_AXIS COMPONENT " // component)
  error = .true.
  return
end select

end subroutine tao_set_qp_axis_struct

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_qp_point_struct (qp_point_name, component, qp_point, value, error, ix_uni)
!
! Routine to set qp_point_names of a qp_point_struct.
!
! Input:
!   qp_point_name   -- character(*): qp_point name. Used for error messages.
!   component       -- character(*): qp_point component name.
!   qp_point        -- qp_point_struct: qp_point_struct with component to modify
!   value           -- character(*): Component value.
!
! Output:
!   qp_point        -- qp_point_struct: qp_point_struct with changed component value.
!   error           -- logical: Set true if there is an error. False otherwise.
!   ix_uni          -- integer, optional: Tao universe number in case the value depends upon
!                       a parameter of a particular universe.
!-

subroutine tao_set_qp_point_struct (qp_point_name, component, qp_point, value, error, ix_uni)

type (qp_point_struct) qp_point
character(*) component, value, qp_point_name
character(*), parameter :: r_name = 'tao_set_qp_point_struct '
integer, optional :: ix_uni
integer iu, ios
logical error

namelist / params / qp_point

!

select case (component)
case ('x')
  call tao_set_real_value(qp_point%x, qp_point_name, value, error, dflt_uni = ix_uni)
case ('y')
  call tao_set_real_value(qp_point%y, qp_point_name, value, error, dflt_uni = ix_uni)
case ('units')
  qp_point%units = value
  error = .false.
case ('')
  iu = tao_open_scratch_file (error);  if (error) return
  write (iu, '(a)') '&params'
  write (iu, '(a)') ' qp_point = ' // trim(value)
  write (iu, '(a)') '/'
  write (iu, *)

  rewind (iu)
  read (iu, nml = params, iostat = ios)
  close (iu, status = 'delete')

  if (ios /= 0) then
    call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
    return
  endif

case default
  call out_io (s_error$, r_name, "BAD GRAPH QP_POINT COMPONENT " // component)
  error = .true.
  return
end select

end subroutine tao_set_qp_point_struct

end module tao_set_mod
