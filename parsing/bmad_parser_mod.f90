!+
! Module bmad_parser_mod
!
! This module is a collection of helper routines used by bmad_parser and bmad_parser2.
! The routines in this module are specifically taylored for bmad_parser and
! bmad_parser2 and cannot, in general, be used otherwise.
!-

module bmad_parser_mod

use bmad_parser_struct
use superimpose_mod
use binary_parser_mod
use wake_mod
use bookkeeper_mod 
use wall3d_mod
use random_mod
use taylor_mod
implicit none

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_set_attribute (how, ele, delim, delim_found, err_flag, pele, check_free, heterogeneous_ele_list, set_field_master)
!
! Subroutine used by bmad_parser and bmad_parser2 to get the value of
! an attribute from the input file and set the appropriate value in an element.
!
! This subroutine is not intended for general use.
!
! Input:
!   how              -- integer: Either def$ if the element is being construct from scratch or redef$ if the element has
!                         already been formed and this is part of a "ele_name[attrib_name] = value" construct.
!   ele              -- ele_struct: Element whose attribute this is.
!   check_free       -- logical, optional: If present and True then an error will be generated
!                          if the attribute is not free to vary. Used by bmad_parser2.
!   heterogeneous_ele_list 
!                    -- logical, optional: If True (default = False), we are parsing something like something like 
!                           "*[tracking_method] = runge_kutta". 
!                         In this case, runge_kutta may not be valid for this ele but this is not an error.
!   set_field_master -- logical, optional: If True (the default) set ele%field_master = T if the
!                         attribute to be set is something like DB_FIELD. If this routine is being
!                         called post lattice parsing, setting ele%field_master is *not* wanted.
!
! Output
!   delim          -- character(1): Delimiter found where the parsing of the input line stops.
!   delim_found    -- logical: Delimiter found? False if end of input command.
!   err_flag       -- logical: Set True if there is a problem parsing the input.
!   pele           -- parser_ele_struct, optional: Structure to hold additional 
!                       information that cannot be stored in the ele argument.
!-

subroutine parser_set_attribute (how, ele, delim, delim_found, err_flag, pele, check_free, heterogeneous_ele_list, set_field_master)

use photon_reflection_mod, only: finalize_reflectivity_tables

implicit none

type (lat_struct), pointer :: lat
type (parser_ele_struct), optional :: pele
type (ele_struct), target ::  ele
type (ele_struct), pointer ::  ele2
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), target, save ::  ele0
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: bele
type (all_pointer_struct), allocatable :: a_ptrs(:)
type (all_pointer_struct) a_ptr
type (wall3d_struct), pointer :: wall3d_arr(:), wall3d
type (wall3d_section_struct), pointer :: section
type (wall3d_vertex_struct), pointer :: v_ptr
type (cylindrical_map_struct), pointer :: cl_map
type (cartesian_map_term1_struct), pointer :: ct_term
type (cartesian_map_term1_struct), allocatable :: ct_terms(:)
type (gen_grad_map_struct), pointer :: gg_map
type (grid_field_struct), pointer :: g_field
type (cartesian_map_struct), pointer :: ct_map
type (ac_kicker_struct), pointer :: ac
type (photon_element_struct), pointer :: ph
type (photon_reflect_table_struct), allocatable :: rt_save(:)
type (photon_reflect_table_struct), pointer :: rt

real(rp) kx, ky, kz, tol, value, coef, r_vec(10), r0(2), vec(100)
real(rp), allocatable :: table(:,:)
real(rp), pointer :: r_ptr

integer i, i2, j, k, n, na, ne, nn, nt, ix_word, how, ix_word1, ix_word2, ios, ix, iy, i_out, ix_coef, switch
integer expn(6), ix_attrib, i_section, ix_v, ix_sec, i_ptr, i_term, ib, ie, im
integer ix_bounds(2), iy_bounds(2), i_vec(2), n_sec, key

character(40) :: word, str_ix, attrib_word, word2, name
character(40), allocatable :: name_list(:)
character(1) delim, delim1, delim2
character(80) str, err_str
character(200) line

logical, target :: delim_found, err_flag, logic, set_done, end_of_file, do_evaluate, hetero_list
logical is_attrib, err_flag2, old_style_input, ok, err
logical, optional :: check_free, heterogeneous_ele_list, set_field_master

! Get next WORD.
! If an overlay or group element then word is just an attribute to control
! [except for a "GROUP[COMMAND] = 0.343" redef construct]

err_flag = .true.  ! assume the worst
call get_next_word (word, ix_word, ':, =(){', delim, delim_found, call_check = .true., err_flag = err); if (err) return
lat => ele%branch%lat

if (ele%key == def_particle_start$ .and. word == 'SIG_E') word = 'SIG_PZ'

! Taylor

hetero_list = logic_option(.false., heterogeneous_ele_list)

if ((ele%key == taylor$ .or. ele%key == hybrid$) .and. delim == '{' .and. word == '') then

  call get_next_word (word, ix_word, ':, =(){', delim, delim_found, call_check = .true.)

  call match_word (word, ['XX', 'XY', 'XZ', 'YX', 'YY', 'YZ', 'ZX', 'ZY', 'ZZ'], i_out, .true., .false.)
  if (i_out > 0) then
    call parser_error ('OLD STYLE, ROTATION MATRIX BASED SPIN TAYLOR MAP FOR ' // ele%name, &
                       'MUST BE CONVERTED TO NEW STYLE QUATERNON BASED MAP.')
    return
  endif

  call match_word (word, ['S1', 'SX', 'SY', 'SZ'], i_out, .true., .false.)
  if (i_out > 0) then
    i_out = i_out + 99 ! Make i_out not in range [1:6]
  else
    read (word, *, iostat = ios) i_out
    if (delim /= ':' .or. ix_word == 0 .or. ios /= 0) then
      call parser_error ('BAD "OUT" COMPONENT: ' // word, 'IN TERM FOR TAYLOR ELEMENT: ' // ele%name)
      return
    endif
  endif

  call parse_evaluate_value (ele%name, coef, lat, delim, delim_found, err_flag, ',|', ele);  if (err_flag) return
  delim2 = delim   ! Save
  expn = 0

  ! Need to check for "{1: xxx |}" case where there are no exponents.
  if (delim2 == '|') then
    call get_next_word (word, ix_word, '}', delim, delim_found)
    ! If there are exponents then rewind the get_next_word call.
    if (ix_word /= 0 .or. delim /= '}') then
      bp_com%parse_line = trim(word) // delim // bp_com%parse_line
      delim = delim2
    endif
  endif

  ! Parse exponents

  do j = 1, 100 
    if (delim == '}') exit
    call parser_get_integer (n, word, ix_word, delim, delim_found, err_flag, 'BAD EXPONENT');   if (err_flag) return
    if (.not. expect_one_of ('} ', .true., ele%name, delim, delim_found)) return
    if (delim2 == ',') then
      select case (j)
      case (6);      if( .not. expect_one_of ('}', .true., ele%name, delim, delim_found)) return
      case default;  if (.not. expect_one_of (' ', .true., ele%name, delim, delim_found)) return
      end select
      expn(j) = n
    else
      ! Where, for example, n = 34, must separate into 3 and 4.
      do
        nn = modulo(n, 10)
        if (nn < 1 .or. nn > 6) then
          call parser_error ('BAD EXPONENT VALUE FOR TAYLOR ELEMENT: ' // ele%name, 'CANNOT PARSE: ' // str)
          return
        endif
        expn(nn) = expn(nn) + 1
        n = (n - nn) / 10
        if (n == 0) exit
      enddo
    endif
  enddo

  call add_this_taylor_term (ele, i_out, coef, expn)

  call get_next_word (word, ix_word, '},', delim, delim_found)
  if (ix_word /= 0 .or. (delim_found .and. delim /= ',')) then
    call parser_error ('BAD TERM ENDING FOR TAYLOR ELEMENT: ' // ele%name, 'CANNOT PARSE: ' // str)
    return
  endif

  return
endif  ! Taylor term

! Overlay, ramper, or group.
! Redef of "slave(i)%y_knot(j)" must be handled here since the y_knot info is in pele (unlike the x_knot info).

if ((ele%key == overlay$ .or. ele%key == group$ .or. ele%key == ramper$) .and. (word /= 'X_KNOT' .or. how /= redef$)) then
  i = attribute_index(ele, word)       ! general attribute search

  if (how == redef$ .and. word == 'SLAVE') then
    if (.not. expect_this ('(', .true., .false., 'NO "(" AFTER: ' // word, ele, delim, delim_found)) return
    call parser_get_integer (n, word, ix_word, delim, delim_found, err_flag, 'BAD INDEX');  if (err_flag) return
    if (.not. expect_this (')', .true., .false., 'NO "(...)" AFTER: ' // word, ele, delim, delim_found)) return
    if (n < 1 .or. n > size(pele%control)) then
      call parser_error ('SLAVE INDEX OUT OF RANGE: ' // int_str(n))
      return
    endif
    if (.not. expect_this ('%', .false., .false., 'NO "%" AFTER: ' // trim(word) // '()', ele, delim, delim_found)) return
    call get_next_word (word2, ix_word, '(]=', delim, delim_found)
    if (word2 /= 'Y_KNOT') THEN
      call parser_error ('MALFORMED SLAVE PARAMETER REDEF')
      return
    endif
    if (.not. expect_this ('(', .true., .false., 'MALFORMED SLAVE PARAMETER REDEF', ele, delim, delim_found)) return
    call parser_get_integer (ix, word, ix_word, delim, delim_found, err_flag, 'BAD INDEX');  if (err_flag) return
    if (.not. expect_this (')=', .true., .false., 'NO "(...)" AFTER: ' // word, ele, delim, delim_found)) return
    if (ix < 1 .or. ix > size(pele%control(n)%y_knot)) then
      call parser_error ('Y_KNOT INDEX OUT OF RANGE: ' // int_str(ix))
      return
    endif
    call parse_evaluate_value (trim(ele%name) // ' ' // word, pele%control(n)%y_knot(ix), lat, delim, delim_found, err_flag, ele = ele)
    return
  endif

  select case (i)
  case (type$, alias$, descrip$, gang$, is_on$, interpolation$)
    ! Handled below

  case (var$)
    if (how == redef$ .or. allocated(ele%control%var)) then
      call parser_error ('RESETTING VAR = {...} IS NOT PERMITTED', 'FOR: ' // ele%name)
      return
    endif
    call get_overlay_group_names(ele, lat, pele, delim, delim_found, .true., err_flag)
    if (size(ele%control%var) > 0) pele%default_attrib = ele%control%var(1)%name
    return

  case (x_knot$)
    if (.not. parse_real_list2 (lat, 'ERROR PARSING X_KNOT POINTS FOR: ' // ele%name, ele%control%x_knot, n, delim, delim_found, 10, '{', ',', '}')) return
    call re_allocate(ele%control%x_knot, n)
    if (.not. expect_one_of (', ', .false., ele%name, delim, delim_found)) return
    err_flag = .false.
    return

  case default
    ! Parse old style control var syntax: "i > num_ele_attrib$" handles accordion_edge for example.

    is_attrib = (attribute_index(0, word) > 0 .or. (ele%key == group$ .and. word == 'COMMAND'))
    if (how == def$ .and. .not. allocated(ele%control%var) .and. (i < 1 .or. i > num_ele_attrib$) .and. is_attrib) then 
      pele%default_attrib = word
      allocate (ele%control)
      allocate (ele%control%var(1))
      if (ele%key == group$) word = 'COMMAND'
      ele%control%var(1)%name = word
      i = 1 + var_offset$
    endif

    !

    if (i < 1) then
      if (hetero_list) then
        err_flag = .false.
        return
      endif
      call parser_error ('BAD OVERLAY/GROUP ATTRIBUTE: ' // word, 'FOR: ' // ele%name)
      return
    endif

    value = 0
    if (delim == '=') then  ! value
      call parse_evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag, ele = ele)
      if (err_flag) return
    endif

    call pointer_to_attribute (ele, word, .true., a_ptr, err_flag, .true.)
    if (err_flag) then
      call parser_error ('')
      return
    endif

    a_ptr%r = value

    if (attrib_free_problem(word)) return
    
    err_flag = .false.
    return
  end select
endif   ! Overlay, Ramper, or Group

! L_pole, N_pole for wiggler/undulator are deprecated in favor of L_period, N_period.

if (ele%key == wiggler$ .or. ele%key == undulator$) then
  if (word == 'L_POLE' .or. word == 'N_POLE') then
    if (.not. expect_one_of ('=', .true., ele%name, delim, delim_found)) return
    call parse_evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag, ele = ele)
    if (err_flag) return
    if (word == 'L_POLE') then
      ele%value(l_period$) = 2.0_rp * value
    else
      ele%value(n_period$) = 0.5_rp * value
    endif
    return
  endif
endif

! Other old-style conversions

if (ele%key == beambeam$) then
  select case (word)
  case ('BETA_A');    word = 'BETA_A_STRONG'
  case ('BETA_B');    word = 'BETA_B_STRONG'
  case ('ALPHA_A');   word = 'ALPHA_A_STRONG'
  case ('ALPHA_B');   word = 'ALPHA_B_STRONG'
  end select
endif

if (word == 'SPINOR_POLARIZATION' .or. word == 'SPINOR_PHI' .or. word == 'SPINOR_THETA' .or. word == 'SPINOR_XI') then
  call parser_error ('DUE TO BOOKKEEPING COMPLICATIONS, THE OLD SPINOR ATTRIBUTES NO LONGER EXIST: ' // word, &
                     'PLEASE CONVERT TO SPIN_X, SPIN_Y, SPIN_Z COMPONENTS.', 'FOR ELEMENT: ' // ele%name)
  return
endif

! Setting n_ref_pass and multipass_ref_energy is no longer valid and
! will be ignored for backwards compatibility.
if (word == 'N_REF_PASS' .or. word == 'MULTIPASS_REF_ENERGY') then
  call parser_error(quote(word) // ' IS NOT SETTABLE. PLEASE REMOVE FROM LATTICE FILE. PARSING WILL PROCEED AS NORMAL.', &
                    'FOR ELEMENT: ' // ele%name, level = s_warn$)
  call parse_evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag, ele = ele) 
  return
endif

!

key = ele%key

select case (word)
case ('HIGHER_ORDER_FRINGE_TYPE')
  call parser_error ('Note: HIGHER_ORDER_FRINGE_TYPE is now no longer used and will be ignored.', &
                     'Remove reference to this parameter to avoid this warning.', level = s_warn$)
  call get_switch (attrib_word, fringe_type_name(1:), ix, err_flag, ele, delim, delim_found)
  return

case ('SPACE_CHARGE_ON')
  call parser_error ('Note: "bmad_com[SPACE_CHARGE_ON]" has been renamed "bmad_com[HIGH_ENERGY_SPACE_CHARGE_ON]"', &
                     'Will run normally...', level = s_warn$)
  word = 'HIGH_ENERGY_SPACE_CHARGE_ON'
  key = def_bmad_com$

case ('COHERENT_SYNCH_RAD_ON')
  call parser_error ('Note: "bmad_com[COHERENT_SYNCH_RAD_ON]" has been renamed "bmad_com[CSR_AND_SPACE_CHARGE_ON]"', &
                     'Will run normally...', level = s_warn$)
  word = 'CSR_AND_SPACE_CHARGE_ON'

case ('X_PITCH_MULT')
  call parser_error ('X_PITCH_MULT no longer supported (necessitated extremely complicated bookkeeping). Will use X_PITCH instead.', &
                     'Will run normally...', level = s_warn$)
  word = 'X_PITCH'

case ('Y_PITCH_MULT')
  call parser_error ('Y_PITCH_MULT no longer supported (necessitated extremely complicated bookkeeping). Will use Y_PITCH instead.', &
                     'Will run normally...', level = s_warn$)
  word = 'Y_PITCH'
end select                    

! Particle_start and bmad_com elements can have attributes that are not part of the element so
! Need to use pointers_to_attribute.

! For historical reasons, a few paramter[...] parameters are actually in bmad_com.

if (word == 'REF')    word = 'REFERENCE' ! allowed abbrev
if (key == rfcavity$ .and. word == 'LAG') word = 'PHI0'   ! For MAD compatibility
if (key == def_parameter$) then
  select case (word)
  case ('ABSOLUTE_TIME_TRACKING', 'APERTURE_LIMIT_ON', 'ELECTRIC_DIPOLE_MOMENT') 
    key = def_bmad_com$  ! "Parameter[absolute_time_tracking]", etc. is deprecated
  case ('PTC_EXACT_MODEL', 'EXACT_MODEL')
    key = def_ptc_com$
    word = 'EXACT_MODEL'
  case ('PTC_EXACT_MISALIGN', 'EXACT_MISALIGN')
    key = def_ptc_com$
    word = 'EXACT_MISALIGN'
  end select
endif

if (key == sbend$ .or. key == rbend$) then
  if (word == 'G_ERR')       word = 'DG'
  if (word == 'B_FIELD_ERR') word = 'DB_FIELD'
endif

if (key == def_particle_start$ .or. key == def_bmad_com$ .or. key == def_space_charge_com$ .or. key == def_ptc_com$) then
  name = ele%name

  if (word(1:4) == 'PTC_') then    ! For backwards compatibility
    name = 'PTC_COM'

  elseif (ele%name == 'PARAMETER') then
    if (word == 'EXACT_MODEL' .or. word == 'EXACT_MISALIGN') then
      name = 'PTC_COM'
    else
      name = 'BMAD_COM'
    endif

  elseif (word == 'SIGMA_CUTOFF') then
    word = 'LSC_SIGMA_CUTOFF'
  endif

  if (delim == '(') then
    ix = index(bp_com%parse_line, '=')
    if (ix == 0) then
      call parser_error ('MALFORMED BMAD_COM OR PARAMETER SET')
      return
    endif
    word = trim(word) // '(' // bp_com%parse_line(:ix-1)
    delim = '='
    call string_trim(bp_com%parse_line(ix+1:), bp_com%parse_line, ix)
  endif

  ! USE_HARD_EDGE_DRIFTS does not exist anymore. Will ignore to preserve backwards compatibility.
  if (word == 'USE_HARD_EDGE_DRIFTS') then
    call parser_get_logical (word, logic, ele%name, delim, delim_found, err_flag) ! Parse rest of line & ignore.
    return
  endif

  if (ele%key == def_parameter$ .and. word == 'APERTURE_LIMIT_ON') then
    call parser_error ('SYNTAX HAS CHANGED: PARAMETER[APERTURE_LIMIT_ON] = ... NEEDS TO BE REPLACED BY BMAD_COM[APERTURE_LIMIT_ON] = ...', &
                       'THIS IS A WARNING ONLY. THE PROGRAM WILL RUN NORMALLY.', level = s_warn$)
  endif

  if (word == 'D_ORB') then
    if (.not. parse_real_list (lat, trim(ele%name) // ' ' // word, bmad_com%d_orb, .true., delim, delim_found)) return
    bp_com%extra%d_orb_set = .true.
    return
  endif

  if (word == 'SPACE_CHARGE_MESH_SIZE') then
    if (.not. parse_integer_list (trim(ele%name) // ' ' // word, lat, space_charge_com%space_charge_mesh_size, .true., delim, delim_found)) return
    bp_com%extra%space_charge_mesh_size_set = .true.
    return
  endif

  if (word == 'CSR3D_MESH_SIZE') then
    if (.not. parse_integer_list (trim(ele%name) // ' ' // word, lat, space_charge_com%csr3d_mesh_size, .true., delim, delim_found)) return
    bp_com%extra%csr3d_mesh_size_set = .true.
    return
  endif

  if (word == 'DIAGNOSTIC_OUTPUT_FILE') then
    call get_next_word (space_charge_com%diagnostic_output_file, ix_word, ',', delim, delim_found)
    bp_com%extra%diagnostic_output_file_set = .true.
    return
  endif

  !

  call pointers_to_attribute (lat, name, word, .false., a_ptrs, err_flag, .false.)
  if (err_flag .or. size(a_ptrs) == 0) then
    call parser_error ('BAD ATTRIBUTE: ' // word, 'FOR ELEMENT: ' // ele%name)
    return
  endif

  if (associated(a_ptrs(1)%r)) then
    call parse_evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag, ele = ele) 
    if (err_flag) return
    a_ptrs(1)%r = value
    if (associated(a_ptrs(1)%r, bmad_com%max_aperture_limit))              bp_com%extra%max_aperture_limit_set          = .true.
    if (associated(a_ptrs(1)%r, bmad_com%default_ds_step))                 bp_com%extra%default_ds_step_set             = .true.
    if (associated(a_ptrs(1)%r, bmad_com%significant_length))              bp_com%extra%significant_length_set          = .true.
    if (associated(a_ptrs(1)%r, bmad_com%rel_tol_tracking))                bp_com%extra%rel_tol_tracking_set            = .true.
    if (associated(a_ptrs(1)%r, bmad_com%abs_tol_tracking))                bp_com%extra%abs_tol_tracking_set            = .true.
    if (associated(a_ptrs(1)%r, bmad_com%rel_tol_adaptive_tracking))       bp_com%extra%rel_tol_adaptive_tracking_set   = .true.
    if (associated(a_ptrs(1)%r, bmad_com%abs_tol_adaptive_tracking))       bp_com%extra%abs_tol_adaptive_tracking_set   = .true.
    if (associated(a_ptrs(1)%r, bmad_com%init_ds_adaptive_tracking))       bp_com%extra%init_ds_adaptive_tracking_set   = .true.
    if (associated(a_ptrs(1)%r, bmad_com%min_ds_adaptive_tracking))        bp_com%extra%min_ds_adaptive_tracking_set    = .true.
    if (associated(a_ptrs(1)%r, bmad_com%fatal_ds_adaptive_tracking))      bp_com%extra%fatal_ds_adaptive_tracking_set  = .true.
    if (associated(a_ptrs(1)%r, bmad_com%autoscale_amp_abs_tol))           bp_com%extra%autoscale_amp_abs_tol_set       = .true.
    if (associated(a_ptrs(1)%r, bmad_com%autoscale_amp_rel_tol))           bp_com%extra%autoscale_amp_rel_tol_set       = .true.
    if (associated(a_ptrs(1)%r, bmad_com%autoscale_phase_tol))             bp_com%extra%autoscale_phase_tol_set         = .true.
    if (associated(a_ptrs(1)%r, bmad_com%electric_dipole_moment))          bp_com%extra%electric_dipole_moment_set      = .true.
    if (associated(a_ptrs(1)%r, bmad_com%synch_rad_scale))                 bp_com%extra%synch_rad_scale_set      = .true.
    if (associated(a_ptrs(1)%r, bmad_com%sad_eps_scale))                   bp_com%extra%sad_eps_scale_set               = .true.
    if (associated(a_ptrs(1)%r, bmad_com%sad_amp_max))                     bp_com%extra%sad_amp_max_set                 = .true.

    if (associated(a_ptrs(1)%r, space_charge_com%ds_track_step))           bp_com%extra%ds_track_step_set               = .true.
    if (associated(a_ptrs(1)%r, space_charge_com%dt_track_step))           bp_com%extra%dt_track_step_set               = .true.
    if (associated(a_ptrs(1)%r, space_charge_com%cathode_strength_cutoff)) bp_com%extra%cathode_strength_cutoff_set     = .true.
    if (associated(a_ptrs(1)%r, space_charge_com%rel_tol_tracking))        bp_com%extra%sc_rel_tol_tracking_set         = .true.
    if (associated(a_ptrs(1)%r, space_charge_com%abs_tol_tracking))        bp_com%extra%sc_abs_tol_tracking_set         = .true.
    if (associated(a_ptrs(1)%r, space_charge_com%beam_chamber_height))     bp_com%extra%beam_chamber_height_set         = .true.
    if (associated(a_ptrs(1)%r, space_charge_com%lsc_sigma_cutoff))        bp_com%extra%lsc_sigma_cutoff_set            = .true.
    if (associated(a_ptrs(1)%r, space_charge_com%particle_sigma_cutoff))   bp_com%extra%particle_sigma_cutoff_set       = .true.

    if (associated(a_ptrs(1)%r, ptc_com%vertical_kick))                    bp_com%extra%vertical_kick_set               = .true.
    if (associated(a_ptrs(1)%r, ptc_com%cut_factor))                       bp_com%extra%cut_factor_set                  = .true.

  elseif (associated(a_ptrs(1)%i)) then
    call parse_evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag, ele = ele) 
    if (err_flag) return
    if (associated(a_ptrs(1)%i, lat%particle_start%direction) .and. nint(value) /= -1 .and. nint(value) /= 1) then
      call parser_error ('VALUE OF PARTICLE_START[DIRECTION] MUST BE -1 OR 1.')
      return
    endif
    a_ptrs(1)%i = nint(value)
    if (associated(a_ptrs(1)%i, bmad_com%taylor_order))                   bp_com%extra%taylor_order_set                    = .true.
    if (associated(a_ptrs(1)%i, bmad_com%default_integ_order))            bp_com%extra%default_integ_order_set             = .true.
    if (associated(a_ptrs(1)%i, bmad_com%runge_kutta_order))              bp_com%extra%runge_kutta_order_set               = .true.
    if (associated(a_ptrs(1)%i, bmad_com%sad_n_div_max))                  bp_com%extra%sad_n_div_max_set                   = .true.
    if (associated(a_ptrs(1)%i, bmad_com%max_num_runge_kutta_step))       bp_com%extra%max_num_runge_kutta_step_set        = .true.

    if (associated(a_ptrs(1)%i, space_charge_com%n_bin))                   bp_com%extra%n_bin_set                       = .true.
    if (associated(a_ptrs(1)%i, space_charge_com%particle_bin_span))       bp_com%extra%particle_bin_span_set           = .true.
    if (associated(a_ptrs(1)%i, space_charge_com%n_shield_images))         bp_com%extra%n_shield_images_set             = .true.
    if (associated(a_ptrs(1)%i, space_charge_com%sc_min_in_bin))           bp_com%extra%sc_min_in_bin_set               = .true.

    if (associated(a_ptrs(1)%i, ptc_com%max_fringe_order))                 bp_com%extra%max_fringe_order_set            = .true.
    if (associated(a_ptrs(1)%i, ptc_com%old_integrator))                   bp_com%extra%old_integrator_set              = .true.

  elseif (associated(a_ptrs(1)%l)) then
    if (associated(a_ptrs(1)%l, bmad_com%auto_bookkeeper)) a_ptrs(1)%l => logic  ! Auto_bookkeeper must not be set.
    call parser_get_logical (word, a_ptrs(1)%l, ele%name, delim, delim_found, err_flag)
    if (err_flag) return
    if (associated(a_ptrs(1)%l, bmad_com%absolute_time_ref_shift))        bp_com%extra%absolute_time_ref_shift_set         = .true.
    if (associated(a_ptrs(1)%l, bmad_com%rf_phase_below_transition_ref))  bp_com%extra%rf_phase_below_transition_ref_set   = .true.
    if (associated(a_ptrs(1)%l, bmad_com%sr_wakes_on))                    bp_com%extra%sr_wakes_on_set                     = .true.
    if (associated(a_ptrs(1)%l, bmad_com%lr_wakes_on))                    bp_com%extra%lr_wakes_on_set                     = .true.
    if (associated(a_ptrs(1)%l, bmad_com%high_energy_space_charge_on))    bp_com%extra%high_energy_space_charge_on_set     = .true.
    if (associated(a_ptrs(1)%l, bmad_com%csr_and_space_charge_on))        bp_com%extra%csr_and_space_charge_on_set         = .true.
    if (associated(a_ptrs(1)%l, bmad_com%spin_tracking_on))               bp_com%extra%spin_tracking_on_set                = .true.
    if (associated(a_ptrs(1)%l, bmad_com%radiation_damping_on))           bp_com%extra%radiation_damping_on_set            = .true.
    if (associated(a_ptrs(1)%l, bmad_com%radiation_zero_average))         bp_com%extra%radiation_zero_average_set          = .true.
    if (associated(a_ptrs(1)%l, bmad_com%radiation_fluctuations_on))      bp_com%extra%radiation_fluctuations_on_set       = .true.
    if (associated(a_ptrs(1)%l, bmad_com%conserve_taylor_maps))           bp_com%extra%conserve_taylor_maps_set            = .true.
    if (associated(a_ptrs(1)%l, bmad_com%absolute_time_tracking))         bp_com%extra%absolute_time_tracking_set          = .true.
    if (associated(a_ptrs(1)%l, bmad_com%convert_to_kinetic_momentum))    bp_com%extra%convert_to_kinetic_momentum_set     = .true.
    if (associated(a_ptrs(1)%l, bmad_com%aperture_limit_on))              bp_com%extra%aperture_limit_on_set               = .true.
    if (associated(a_ptrs(1)%l, bmad_com%debug))                          bp_com%extra%debug_set                           = .true.

    if (associated(a_ptrs(1)%l, space_charge_com%lsc_kick_transverse_dependence)) bp_com%extra%lsc_kick_transverse_dependence_set = .true.
    if (associated(a_ptrs(1)%l, space_charge_com%debug))                   bp_com%extra%sc_debug_set                    = .true.

    if (associated(a_ptrs(1)%l, ptc_com%use_orientation_patches))          bp_com%extra%use_orientation_patches_set     = .true.
    if (associated(a_ptrs(1)%l, ptc_com%print_info_messages))              bp_com%extra%print_info_messages_set         = .true.
    if (associated(a_ptrs(1)%l, ptc_com%exact_model))                      bp_com%extra%exact_model_set                 = .true.
    if (associated(a_ptrs(1)%l, ptc_com%exact_misalign))                   bp_com%extra%exact_misalign_set              = .true.
    if (associated(a_ptrs(1)%l, ptc_com%translate_patch_drift_time))       bp_com%extra%translate_patch_drift_time_set  = .true.

  else
    call parser_error ('BOOKKEEPING ERROR. PLEASE CONTACT A BMAD MAINTAINER!')
  endif

  return
endif

! Redefinition of Long-range "wake()", "r_custom()", etc.
! Old style wiggler "term()" handled below.

if (delim == '(' .and. .not. (word == 'TERM' .and. how == def$)) then
  word2 = trim(word) // '('
  call get_next_word (word, ix_word, '=', delim, delim_found)
  word = trim(word2) // word

  if (.not. delim_found) then
    call parser_error ('NO "=" SIGN FOUND', 'FOR ELEMENT: ' // ele%name)
    return
  endif

  call pointer_to_attribute (ele, word, how == def$, a_ptr, err_flag, .false.)

  if (err_flag .or. (.not. associated(a_ptr%r) .and. .not. associated(a_ptr%i) .and. .not. associated(a_ptr%l))) then
    call parser_error ('BAD ATTRIBUTE: ' // word, 'FOR ELEMENT: ' // ele%name)
    return
  endif

  if (associated(a_ptr%r)) then
    call parse_evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag, ele = ele)
    a_ptr%r = value
  elseif (associated(a_ptr%i)) then
    if (index(word, '%MASTER_PARAMETER') /= 0) then
      call get_next_word (word2, ix_word, ',', delim, delim_found)
      ix = attribute_index(ele, word2)
      if (ix < 1) then
        call parser_error ('BAD MASTER_PARAMETER NAME FOR ELEMENT: ' // ele%name)
        return
      endif
      a_ptr%i = ix
    elseif (index(word, '%ELE_ANCHOR_PT') /= 0) then
      call get_switch ('WALL ELE_ANCHOR_PT', anchor_pt_name(1:), a_ptr%i, err_flag, ele, delim, delim_found)
    else
      call parse_evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag, ele = ele)
      a_ptr%i = nint(value)
    endif
  else
    call parser_get_logical (word, a_ptr%l, ele%name, delim, delim_found, err_flag)
  endif
  return
endif

! "WALL%" redef

if (word(1:5) == 'WALL%') then

  select case (word(6:))

  ! Section redef

  case ('SECTION')

    if (delim /= '(') then
      call parser_error ('MALFORMED WALL COMPONENT REDEF IN ELEMENT: ' // ele%name)
      return
    endif

    ix_sec = evaluate_array_index (err_flag, ')', word2, '(=', delim)
    n_sec = -1
    if (associated(ele%wall3d)) n_sec = size(ele%wall3d(1)%section)
    if (err_flag .or. ix_sec < 0 .or. ix_sec > n_sec) then
      call parser_error('BAD ' // trim(word) // ' INDEX', 'FOR ELEMENT: ' // ele%name)
      return
    endif
    section => ele%wall3d(1)%section(ix_sec)

    if (word2 == '%S' .and. delim == '=') then
      r_ptr => section%s

    elseif (word2 == '%DR_DS' .and. delim == '=') then
      r_ptr => section%dr_ds

    elseif (word2 == '%V' .and. delim == '(') then
      ix_v = evaluate_array_index (err_flag, ')', word, '=', delim)
      if (err_flag .or. ix_v < 0 .or. ix_v > size(section%v)) then
        call parser_error('BAD VERTEX INDEX',  'FOR ELEMENT: ' // ele%name)
        return
      endif
      v_ptr => section%v(ix_v)

      select case (word)
      case ('%X');        r_ptr => v_ptr%x
      case ('%Y');        r_ptr => v_ptr%y
      case ('%RADIUS_X'); r_ptr => v_ptr%radius_x
      case ('%RADIUS_Y'); r_ptr => v_ptr%radius_y
      case ('%TILT');     r_ptr => v_ptr%tilt
      case default
        call parser_error('BAD WALL SECTION VERTEX COMPONENT: ' // word, 'FOR ELEMENT: ' // ele%name)
        return
      end select
    else
      call parser_error('BAD WALL SECTION COMPONENT: ' // word2, 'FOR ELEMENT: ' // ele%name)
      return
    endif

    call parse_evaluate_value (ele%name, r_ptr, lat, delim, delim_found, err_flag, ele = ele)

  ! Not recognized

  case default
    call parser_error ('BAD WALL COMPONENT REDEF: ' // word, 'IN ELEMENT: ' // ele%name)
  end select

  return
endif

! if not an overlay/group then see if it is an ordinary attribute.
! if not an ordinary attribute then might be a superimpose switch

if (ix_word == 0) then  ! no word
  call parser_error  ('"," NOT FOLLOWED BY ATTRIBUTE NAME FOR: ' // ele%name)
  return
endif


select case (word)
case ('ELE_BEGINNING', 'ELE_CENTER', 'END_END', 'REF_BEGINNING', 'REF_CENTER', 'REF_END')
  call parser_error ('OLD SUPERPOSITION SYNTAX: ' // word, &
              'PLEASE CONVERT (SEE THE BMAD MANUAL)', 'WARNING ONLY, PROGRAM WILL RUN NORMALLY...', level = s_warn$)
end select

select case (word)
case ('TILT')
  if (ele%key == sbend$ .or. ele%key == rbend$ .or. ele%key == rf_bend$) then
    call parser_error ('BENDS HAVE A "REF_TILT" ATTRIBUTE BUT NOT A "TILT" ATTRIBUTE.')
  endif

case ('DPHI0')
  call parser_error ('THE ATTRIBUTE NAME "DPHI0" HAS BEEN CHANGED TO "PHI0_MULTIPASS"', &
                     'PLEASE MAKE THE CHANGE IN THE LATTICE FILE.', &
                     '[THIS IS A WARNING ONLY. THIS PROGRAM WILL RUN NORMALLY]', level = s_warn$)
  word = 'PHI0_MULTIPASS'

case ('REL_TRACKING_CHARGE') 
  call parser_error ('THE ATTRIBUTE NAME "REL_TRACKING_CHARGE" HAS BEEN CHANGED TO "DEFAULT_TRACKING_SPECIES"', &
                     'PLEASE MAKE THE CHANGE IN THE LATTICE FILE.')

case ('RADIUS')
  call parser_error ('THE ATTRIBUTE "RADIUS" HAS BEEN CHANGED TO "R0_MAG"', &
                     'PLEASE MAKE THE CHANGE IN THE LATTICE FILE.', &
                     '[THIS IS A WARNING ONLY. THIS PROGRAM WILL RUN NORMALLY]', level = s_warn$)
  word = 'R0_MAG'

case ('FIELD')
  call parser_error ('THE "FIELD = {...}" SYNTAX HAS BEEN CHANGED TO "GRID_FIELD = {...} AND/OR CYLINDRICAL_MAP = {...}.', &
                     'NOTE: THIS INVOLVES MORE THAN CHANGING "FIELD" TO "GRID_FIELD" OR "CYLINDRICAL_MAP".', &
                     'PLEASE SEE THE BMAD MANUAL FOR MORE DETAILS')

case ('REF_BEGINNING')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  pele%ref_pt = anchor_beginning$
  err_flag = .false.
  return

case ('REF_CENTER')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  pele%ref_pt = anchor_center$
  err_flag = .false.
  return

case ('REF_END')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  pele%ref_pt = anchor_end$
  err_flag = .false.
  return

case ('ELE_BEGINNING')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  pele%ele_pt = anchor_beginning$
  err_flag = .false.
  return

case ('ELE_CENTER')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  pele%ele_pt = anchor_center$
  err_flag = .false.
  return

case ('ELE_END')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  pele%ele_pt = anchor_end$
  err_flag = .false.
  return
end select

if (word(1:16) == 'CUSTOM_ATTRIBUTE') then
  read(word(17:), *, iostat = ios) k
  if (ios /= 0) then
    call parser_error ('BAD NUMBER FOR: PARAMETER[' // trim(word) // ']')
    err_flag = .true.
  endif

  call get_next_word (str, ix_word, ',= ', delim, delim_found, .false.) 
  str = unquote(str)
  call set_custom_attribute_name(str, err_flag, k)
  if (err_flag) call parser_error ('CANNOT SET PARAMETER[' // trim(word) // ']')
  return
endif

!-------------------------------------------------------------------
! "SURFACE" is old style

ix_attrib = attribute_index(ele, word, attrib_word)
if (attrib_free_problem(word)) return

if (word == "SURFACE") then
  ix_attrib = 999  
  attrib_word = word
endif

if (ix_attrib < 1) then
  call pointer_to_attribute(ele, word, .true., a_ptr, err_flag, .false.)
  if (associated(a_ptr%r) .or. associated(a_ptr%i) .or. associated(a_ptr%l)) then
    attrib_word = word
  else
    if (ele%key == drift$ .and. (word == 'HKICK' .or. word == 'VKICK' .or. &
          word == 'BL_HKICK' .or. word == 'BL_VKICK')) then
      call parser_error ('BAD ATTRIBUTE: ' // word, 'FOR ELEMENT: ' // ele%name, &
                        'ONE SOLUTION IS TO MAKE THIS DRIFT A "PIPE" ELEMENT.')
    else
      call parser_error ('BAD ATTRIBUTE NAME: ' // word, 'FOR ELEMENT: ' // ele%name)
    endif
    return
  endif
endif

! ac_kicker amp_vs_time

if (attrib_word == 'AMP_VS_TIME') then
  ac => ele%ac_kick
  if (.not. parse_real_lists (lat, ele, trim(ele%name) // ' AMP_VS_TIME', table, 2, delim, delim_found)) return
  if (.not. expect_one_of (', ', .false., ele%name, delim, delim_found)) return
  n = size(table, 1)
  allocate (ac%amp_vs_time(n))
  do i = 1, n
    ac%amp_vs_time(i)%time = table(i,1)
    ac%amp_vs_time(i)%amp  = table(i,2)
    ac%amp_vs_time(i)%spline%x0 = table(i,1)
    ac%amp_vs_time(i)%spline%y0 = table(i,2)
  enddo
  call spline_akima (ac%amp_vs_time%spline, ok)
  if (.not. ok) call parser_error ('ERROR CREATING SPLINE FOR AC_KICKER AMP_VS_TIME CURVE.', 'FOR ELEMENT: ' // ele%name)
  err_flag = .false.
  return
endif

if (attrib_word == 'FREQUENCIES') then
  ac => ele%ac_kick
  if (.not. parse_real_lists (lat, ele, trim(ele%name) // ' FREQUENCIES', table, 3, delim, delim_found)) return
  if (.not. expect_one_of (', ', .false., ele%name, delim, delim_found)) return
  n = size(table, 1)
  allocate (ac%frequency(n))
  do i = 1, n
    ac%frequency(i)%f    = table(i,1)
    ac%frequency(i)%amp  = table(i,2)
    ac%frequency(i)%phi  = table(i,3)
  enddo
  err_flag = .false.
  return
endif

! ac_kicker frequencies

if (attrib_word == 'FREQUENCIES') then

  return
endif

! wall cross-section definition

if (attrib_word == 'WALL') then

  i_section = 0
  if (.not. expect_this ('=', .true., .true., 'AFTER "WALL"', ele, delim, delim_found)) return
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)

  ! "ele1[wall] = ele2[wall]" construct

  if (delim == '[') then
    ele2 => parser_find_ele_for_attrib_transfer ('WALL', word)
    if (.not. associated(ele2%wall3d)) then
      call parser_error ('NO WALL ASSOCIATED WITH LATTICE ELEMENT: ' // word)
      return
    endif
    call transfer_wall3d(ele2%wall3d, ele%wall3d)
    return
  endif

  !

  if (.not. expect_this ('{', .true., .true., 'AFTER "WALL"', ele, delim, delim_found)) return

  ! Loop wall3d_struct components.

  if (associated (ele%wall3d)) then
    n = size(ele%wall3d)
    wall3d_arr => ele%wall3d
    allocate (ele%wall3d(n+1))
    do i = 1, n
      wall3d => ele%wall3d(i)
      n_sec = size(wall3d_arr(i)%section)
      allocate(wall3d%section(n_sec))
      do ix_sec = 1, n_sec
        nn = size(wall3d%section(ix_sec)%v)
        allocate(wall3d%section(ix_sec)%v(nn))
      enddo
      wall3d = wall3d_arr(i)
      wall3d%n_link = 1
    enddo
    call unlink_wall3d (wall3d_arr)
    wall3d => ele%wall3d(n+1)
  else
    allocate (ele%wall3d(1))
    wall3d => ele%wall3d(1)
  endif

  ! Can imagine in the future that an element could have different types of walls.
  ! Right now this is not true.

  if (ele%key == mask$ .or. ele%key == diffraction_plate$) then
    wall3d%type = mask_plate$
  else
    wall3d%type = chamber_wall$
  endif

  wall3d_loop: do

    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

    ! Possible "}" is end of wall 
    if (delim /= '}' .and. word == '') exit
    if (.not. expect_this ('=', .true., .false., 'AFTER ' // trim(word) // ' IN WALL CONSTRUCT', ele, delim, delim_found)) return

    select case (word)

    case ('NAME')
      call bmad_parser_string_attribute_set (ele, word, delim, delim_found, str_out = wall3d%name)

    case ('OPAQUE_MATERIAL') 
      call bmad_parser_string_attribute_set (ele, word, delim, delim_found, str_out = wall3d%opaque_material)

    case ('CLEAR_MATERIAL') 
      call bmad_parser_string_attribute_set (ele, word, delim, delim_found, str_out = wall3d%clear_material)

    case ('THICKNESS') 
      call parse_evaluate_value (ele%name, wall3d%thickness, lat, delim, delim_found, err_flag, ',}', ele)
      if (err_flag) return

    case ('ELE_ANCHOR_PT')
      call get_switch ('WALL ELE_ANCHOR_PT', anchor_pt_name(1:), wall3d%ele_anchor_pt, err_flag2, ele, delim, delim_found)
      if (err_flag2) return

    case ('SUPERIMPOSE')
      call parser_get_logical ('WALL SUPERIMPOSE', wall3d%superimpose, ele%name, delim, delim_found, err_flag2); if (err_flag2) return

    ! Must be "section = {"

    case ('SECTION')

      ! Read in section

      if (.not. expect_this ('{', .false., .true., 'AFTER "SECTION =" IN WALL CONSTRUCT', ele, delim, delim_found)) return

      i_section = i_section + 1
      ix_v = 0
      call re_allocate (wall3d%section, i_section)
      section => wall3d%section(i_section)

      wall3d_section_loop: do

        call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

        ! Possible "}" is end of wall 
        if (delim /= '}' .and. word == '') exit
        if (word == 'V') then
          if (.not. expect_this ('(', .true., .false., 'AFTER ' // trim(word) // ' IN WALL CONSTRUCT', ele, delim, delim_found)) return
        else
          if (.not. expect_this ('=', .true., .false., 'AFTER ' // trim(word) // ' IN WALL CONSTRUCT', ele, delim, delim_found)) return
        endif

        select case (word)

        case ('TYPE') 
          call get_switch ('WALL TYPE', wall3d_section_type_name(1:), section%type, err_flag2, ele, delim, delim_found)
          if (err_flag2) return

        case ('MATERIAL') 
          call bmad_parser_string_attribute_set (ele, word, delim, delim_found, str_out = section%material)

        case ('THICKNESS')
          call parse_evaluate_value (trim(ele%name) // ' ' // word, section%thickness, lat, delim, delim_found, err_flag, ',}', ele)
          if (err_flag) return
          if (ele%key == capillary$) ele%value(l$) = section%s

        case ('S')
          call parse_evaluate_value (trim(ele%name) // ' ' // word, section%s, lat, delim, delim_found, err_flag, ',}', ele)
          if (err_flag) return
          if (ele%key == capillary$) ele%value(l$) = section%s

        case ('DR_DS') 
          call parse_evaluate_value (trim(ele%name) // ' ' // word, section%dr_ds, lat, delim, delim_found, err_flag, ',}', ele)
          if (err_flag) return
                  
        case ('ABSOLUTE_VERTICES') 
          call parser_get_logical (word, logic, ele%name, delim, delim_found, err_flag)
          if (err_flag) return
          if (logic) then
            section%vertices_state = absolute$
          else
            section%vertices_state = relative$
          endif

        case ('X0') 
          call parse_evaluate_value (trim(ele%name) // ' ' // word, section%r0(1), lat, delim, delim_found, err_flag, ',}', ele)
          if (err_flag) return

        case ('Y0') 
          call parse_evaluate_value (trim(ele%name) // ' ' // word, section%r0(2), lat, delim, delim_found, err_flag, ',}', ele)
          if (err_flag) return

        case ('R0')
          if (.not. parse_real_list (lat, trim(ele%name) // ' SECTION R0', section%r0, .true., delim, delim_found)) return
          if (.not. expect_one_of (',}', .false., ele%name, delim, delim_found)) return

        ! Parse "V() = ..." constructs.

        case ('V')

          ix_v = ix_v + 1
          section%n_vertex_input = ix_v
          call re_allocate (section%v, ix_v)

          call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
          read (word, *, iostat = ios) j 
          if (ios /= 0 .or. ix_v /= j) then
            call parser_error ('BAD OR OUT OF ORDER WALL SECTION VERTEX INDEX NUMBER FOR: ' // ele%name)
            return
          endif

          if (.not. expect_this (')={', .true., .false., 'AFTER "V(n)" IN WALL CONSTRUCT', ele, delim, delim_found)) return

          call parse_evaluate_value (trim(ele%name), section%v(ix_v)%x, lat, delim, delim_found, err_flag, ',', ele)
          if (err_flag) return

          call parse_evaluate_value (trim(ele%name), section%v(ix_v)%y, lat, delim, delim_found, err_flag, ',}', ele)
          if (err_flag) return

          if (delim == ',') then
            call parse_evaluate_value (trim(ele%name), section%v(ix_v)%radius_x, lat, delim, delim_found, err_flag, ',}', ele)
            if (err_flag) return
          endif

          if (delim == ',') then
            call parse_evaluate_value (trim(ele%name), section%v(ix_v)%radius_y, lat, delim, delim_found, err_flag, ',}', ele)
            if (err_flag) return
          endif

          if (delim == ',') then
            call parse_evaluate_value (trim(ele%name), section%v(ix_v)%tilt, lat, delim, delim_found, err_flag, '}', ele)
            if (err_flag) return
          endif

          call get_next_word (word, ix_word, '{},()=', delim, delim_found)
          if (word /= '' .or. (delim /= '}' .and. delim /= ',')) then
            call parser_error ('BAD SYNTAX IN WALL SECTION DEFINITION FOR ELEMENT: ' // ele%name)
            return
          endif

        case default
          call parser_error ('WALL SECTION COMPONENT NOT RECOGNIZED: ' // word, 'FOR ELEMENT: ' // ele%name)
          return
        end select   ! section components

        if (.not. expect_one_of (',}', .true., ele%name, delim, delim_found)) return
        if (delim == '}') then
          if (.not. expect_one_of(',}', .false., ele%name, delim, delim_found)) return
          exit
        endif
      enddo wall3d_section_loop

    case default
      call parser_error ('WALL COMPONENT NOT RECOGNIZED: ' // word, 'FOR ELEMENT: ' // ele%name)
      return
    end select   ! wall components

    if (.not. expect_one_of (',}', .true., ele%name, delim, delim_found)) return
    if (delim == '}') exit

  enddo wall3d_loop

  ! Next thing on line should be either a "," or end-of-line

  logic = expect_one_of(', ', .false., ele%name, delim, delim_found)
  return

endif

!-------------------------------
! Reflecting Surface

select case (attrib_word)
case ('ENERGY_PROBABILITY_CURVE')
  ph => ele%photon
  nt = 0
  if (.not. allocated(ph%init_energy_prob)) allocate(ph%init_energy_prob(100))
  if (.not. expect_this ('={', .true., .true., 'AFTER ' // quote(attrib_word), ele, delim, delim_found)) return
  call parser_call_check(word, ix_word, delim, delim_found)
  do
    nt = nt + 1
    if (nt > size(ph%init_energy_prob)) call reallocate_spline(ph%init_energy_prob, 2*nt)
    if (.not. parser_fast_real_read(vec(:2), ele, ' ,}', delim, '', .true.)) return
    ph%init_energy_prob(nt)%x0 = vec(1); ph%init_energy_prob(nt)%y0 = vec(2)
    if (delim == '}') exit
  enddo

  call reallocate_spline(ph%init_energy_prob, nt)
  call spline_akima(ph%init_energy_prob, ok)
  call re_allocate(ph%integrated_init_energy_prob, nt)
  ph%integrated_init_energy_prob(1) = 0
  do i = 2, nt
    ph%integrated_init_energy_prob(i) = ph%integrated_init_energy_prob(i-1) + &
                      spline1(ph%init_energy_prob(i-1), ph%init_energy_prob(i-1)%x1, -1)
  enddo

  if (.not. expect_one_of(', ', .false., ele%name, delim, delim_found)) return
  err_flag = .false.
  return

case ('REFLECTIVITY_TABLE')
  ph => ele%photon
  nt = 0
  if (allocated(ph%reflectivity_table)) deallocate(ph%reflectivity_table)
  if (.not. expect_this ('={', .true., .true., 'AFTER ' // quote(attrib_word), ele, delim, delim_found)) return

  do
    nt = nt + 1
    if (nt == 1) then
      allocate(ph%reflectivity_table(nt))
    else
      call move_alloc(ph%reflectivity_table, rt_save)
      allocate (ph%reflectivity_table(nt))
      ph%reflectivity_table(1:nt-1) = rt_save
    endif
    rt => ph%reflectivity_table(nt)
    if (.not. expect_this ('{', .false., .true., 'AFTER ' // quote(attrib_word), ele, delim, delim_found)) return
    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
    if (word /= 'ANGLES') then
      call parser_error ('EXPECTING "ANGLES" ATTRIBUTE IN REFLECTIVITY_TABLE CONSTRUCT FOR ELEMENT: ' // ele%name)
      return
    endif
    if (.not. expect_this ('=', .true., .false., 'AFTER ' // quote(attrib_word), ele, delim, delim_found)) return
    if (.not. parse_real_list2(lat, trim(ele%name) // ' ANGLES', rt%angle, na, delim, delim_found)) return
    if (.not. expect_this (',', .false., .false., 'AFTER ' // quote(attrib_word), ele, delim, delim_found)) return
    ne = 0
    do
      ne = ne + 1
      call re_allocate(rt%energy, ne)
      call re_allocate2d(rt%p_reflect, na, ne)
      call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
      if (word /= 'P_REFLECT') then
        call parser_error ('EXPECTING "P_REFLECT" ATTRIBUTE IN REFLECTIVITY_TABLE CONSTRUCT FOR ELEMENT: ' // ele%name)
        return
      endif
      if (.not. expect_this ('=', .true., .false., 'AFTER ' // quote(attrib_word), ele, delim, delim_found)) return
      call parse_evaluate_value (ele%name, rt%energy(ne), lat, delim, delim_found, err_flag, ',', ele);  if (err_flag) return
      if (.not. parse_real_list (lat, trim(ele%name) // ' P_REFLECT', rt%p_reflect(:,ne), .true., delim, delim_found)) return
      rt%max_energy = vec(1)
      if (.not. expect_one_of(',}', .false., ele%name, delim, delim_found)) return
      if (delim == '}') exit
    enddo
    allocate(rt%int1(ne))
    if (.not. expect_one_of(',}', .false., ele%name, delim, delim_found)) return
    if (delim == '}') exit
  enddo

  call finalize_reflectivity_tables (ph%reflectivity_table, .false.)

  if (.not. expect_one_of(', ', .false., ele%name, delim, delim_found)) return
  err_flag = .false.
  return

!

case ('SURFACE', 'PIXEL', 'DISPLACEMENT', 'H_MISALIGN', 'SEGMENTED')
  ph => ele%photon

  name = attrib_word
  if (attrib_word == 'SURFACE') then
    if (ele%key == detector$) name = 'PIXEL'
  elseif (attrib_word /= 'PIXEL') then
    call match_word (attrib_word, surface_grid_type_name(1:), ph%grid%type)
  endif

  if (.not. expect_this ('=', .true., .true., 'AFTER ' // quote(attrib_word), ele, delim, delim_found)) return
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)

  ! "ele1[surface] = ele2[surface]" construct

  if (delim == '[') then
    ele2 => parser_find_ele_for_attrib_transfer (attrib_word, word)
    if (.not. associated(ele2%photon)) then
      call parser_error ('NO ' // trim(attrib_word) // ' ASSOCIATED WITH LATTICE ELEMENT: ' // word)
      return
    endif

    if (attrib_word == 'PIXEL') then
      ph%pixel = ele2%photon%pixel
    else
      ph%grid = ele2%photon%grid
    endif

    return
  endif

  !

  if (.not. expect_this ('{', .true., .true., 'AFTER ' // quote(attrib_word), ele, delim, delim_found)) return

  if (attrib_word == 'SURFACE') then   ! Old style
!!!    call parser_error ('Old style "SURFACE" construct. Please change this (see Bmad manual).', level = s_warn$)
    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
    ! Expect "GRID = {" 
    if (word /= 'GRID') then
      call parser_error ('EXPECT "GRID" AFTER "SURFACE" BUT GOT: ' // word, 'FOR: ' // ele%name)
      return
    endif
    if (.not. expect_this ('={', .true., .true., 'AFTER "GRID"', ele, delim, delim_found)) return
  endif

  ix_bounds = int_garbage$; iy_bounds = int_garbage$

  do
    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
    if (word /= 'PT') then
      if (.not. expect_this ('=', .true., .false., 'AFTER ' // trim(word) // ' IN ' // trim(attrib_word) // ' CONSTRUCT', ele, delim, delim_found)) return
    endif

    select case (word)
    case ('TYPE')   ! This is old style.
      call get_switch ('SURFACE GRID TYPE', surface_grid_type_name(1:), ph%grid%type, err_flag2, ele, delim, delim_found)
      if (err_flag2) return
      bp_com%parse_line = delim // bp_com%parse_line

    case ('ACTIVE')
      call parser_get_logical (word, ph%grid%active, ele%name, delim, delim_found, err_flag2); if (err_flag2) return

    case ('DR')
      if (name == 'PIXEL') then
        if (.not. parse_real_list (lat, trim(ele%name) // ' GRID DR', ph%pixel%dr, .true., delim, delim_found)) return
      else
        if (.not. parse_real_list (lat, trim(ele%name) // ' GRID DR', ph%grid%dr, .true., delim, delim_found)) return
      endif

    case ('R0')
      if (name == 'PIXEL') then
        if (.not. parse_real_list (lat, trim(ele%name) // ' GRID R0', ph%pixel%r0, .true., delim, delim_found)) return
      else
        if (.not. parse_real_list (lat, trim(ele%name) // ' GRID R0', ph%grid%r0, .true., delim, delim_found)) return
      endif

    case ('IX_BOUNDS', 'IY_BOUNDS')
      if (.not. parse_integer_list (trim(ele%name) // ' GRID ' // trim(word), lat, i_vec, .true., delim, delim_found)) return
      if (word == 'IX_BOUNDS') ix_bounds = i_vec
      if (word == 'IY_BOUNDS') iy_bounds = i_vec

      if (any(ix_bounds /= int_garbage$) .and. any(iy_bounds /= int_garbage$)) then
        if (any(ix_bounds == int_garbage$) .or. any(iy_bounds == int_garbage$) .or. &
            ix_bounds(1) > ix_bounds(2) .or. iy_bounds(1) > iy_bounds(2)) then
          call parser_error ('SURFACE GRID X/IY_BOUNDS NOT PROPERLY SET', trim(ele%name))
          return
        endif

        if (name == 'PIXEL') then
          if (allocated (ph%pixel%pt)) deallocate (ph%pixel%pt)
          allocate (ph%pixel%pt(ix_bounds(1):ix_bounds(2), iy_bounds(1):iy_bounds(2)))
        else
          if (allocated (ph%grid%pt)) deallocate (ph%grid%pt)
          allocate (ph%grid%pt(ix_bounds(1):ix_bounds(2), iy_bounds(1):iy_bounds(2)))
        endif
      endif

    case ('PT')
      bp_com%parse_line = delim // bp_com%parse_line
      if (.not. parse_integer_list (trim(ele%name) // ' GRID PT', lat, i_vec, .true., delim, delim_found)) return

      if (.not. allocated(ph%grid%pt)) then
        call parser_error ('IX_BOUNDS OR IY_BOUNDS MISSING WHEN CONSTRUCTING: ' // attrib_word, 'FOR: ' // ele%name)
        return
      endif

      if (any(i_vec < lbound(ph%grid%pt)) .or. any(i_vec > ubound(ph%grid%pt))) then
        call parser_error ('PT(I,J) INDEX OUT OF BOUNDS WHEN CONSTRUCTING: ' // attrib_word, 'FOR: ' // ele%name)
        return
      endif

      if (.not. expect_this ('=', .false., .false., 'IN GRID PT', ele, delim, delim_found)) return

      if (ph%grid%type == h_misalign$) then
        if (.not. parse_real_list (lat, trim(ele%name) // 'IN GRID PT', r_vec(1:4), .true., delim, delim_found)) return
        ph%grid%pt(i_vec(1), i_vec(2))%orientation = surface_orientation_struct(r_vec(1), r_vec(2), r_vec(3), r_vec(4))

      elseif (ph%grid%type == displacement$) then
        r_vec(1:4) = real_garbage$
        if (.not. parse_real_list (lat, trim(ele%name) // 'IN GRID PT', r_vec(1:4), .false., delim, delim_found, num_found = n)) return
        if (n /= 1 .and. n /= 3 .and. n /= 4) then
          call parser_error ('NUMBER OF PT(I,J) VALUES NOT 1, 3, NOR 4 FOR SURFACE GRID OF: ' // ele%name)
          return
        endif
        ph%grid%pt(i_vec(1), i_vec(2))%z0 = r_vec(1)
        ph%grid%pt(i_vec(1), i_vec(2))%dz_dx = r_vec(2)
        ph%grid%pt(i_vec(1), i_vec(2))%dz_dy = r_vec(3)
        ph%grid%pt(i_vec(1), i_vec(2))%d2z_dxdy = r_vec(4)
      elseif (ph%grid%type == not_set$) then
        call parser_error ('THE SURFACE GRID TYPE MUST BE SET BEFORE SETTING A TABLE OF SURFACE GRID "PT" POINTS.', &
                           'FOR: ' // ele%name)
        return
      else
        call parser_error ('A TABLE OF SURFACE GRID "PT" POINTS IS NOT ALLOWED IF THE GRID TYPE IS', &
                           'SOMETHING OTHER THAN "OFFSET" OR "H_MISALIGN" FOR: ' // ele%name)
        return
      endif

    case default
      call parser_error ('GRID COMPONENT NOT RECOGNIZED: ' // word, 'FOR ELEMENT: ' // ele%name)
      return
    end select

    if (.not. expect_one_of (',}', .false., ele%name, delim, delim_found)) return

    call string_trim(bp_com%parse_line, bp_com%parse_line, ix)
    if (word == 'PT' .and. delim == ',' .and. bp_com%parse_line(1:1) == '}') then
      delim = '}'
      bp_com%parse_line = bp_com%parse_line(2:)
    endif

    if (delim == '}') then
      if (attrib_word == 'SURFACE') then
        if (.not. expect_one_of (',}', .false., ele%name, delim, delim_found)) return
      endif
      exit
    endif
  enddo

  if (.not. expect_one_of(', ', .false., ele%name, delim, delim_found)) return
  err_flag = .false.
  return

!------------------------
! Curvature

case ('CURVATURE')
  ph => ele%photon

  if (.not. expect_this ('=', .true., .true., 'AFTER ' // quote(attrib_word), ele, delim, delim_found)) return
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)

  if (delim == '[') then
    ele2 => parser_find_ele_for_attrib_transfer (attrib_word, word)
    if (.not. associated(ele2%photon)) then
      call parser_error ('NO ' // trim(attrib_word) // ' ASSOCIATED WITH LATTICE ELEMENT: ' // word)
      return
    endif
    ph%curvature = ele2%photon%curvature
  endif

  if (.not. expect_this ('{', .true., .true., 'AFTER ' // quote(attrib_word), ele, delim, delim_found)) return

  do
    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
    call pointer_to_attribute (ele, 'CURVATURE%' // word, .false., a_ptr, err_flag)
    if (err_flag) then
      call parser_error ('BAD CURVATURE PARAMETER: ' // word, 'FOR: ' // ele%name)
      return
    endif
    call parse_evaluate_value (trim(ele%name) // ' ' // attrib_word, a_ptr%r, &
                                                             lat, delim, delim_found, err_flag, ele = ele)

    if (.not. expect_one_of (',}', .true., ele%name, delim, delim_found)) return
    if (delim == '}') exit
  enddo

  if (.not. expect_one_of(', ', .false., ele%name, delim, delim_found)) return
  err_flag = .false.
  return
end select

!-------------------------------

if (attrib_word == 'SR_WAKE') then
  if (.not. expect_this ('=', .true., .true., 'AFTER "SR_WAKE"', ele, delim, delim_found)) return
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)
  ! ele1[sr_wake] = ele2[sr_wake] construct.
  if (delim == '[') then
    ele2 => parser_find_ele_for_attrib_transfer ('SR_WAKE', word); if (err_flag) return
    if (.not. associated(ele%wake)) allocate (ele%wake)
    if (.not. associated(ele2%wake)) then
      call parser_error ('SR_WAKE NOT DEFINED FOR: ' // ele2%name)
      return
    endif
    ele%wake%sr = ele2%wake%sr
  ! "ele1[sr_wake] = call::..." or "ele1: ..., sr_wake = {...}, ..." construct.
  else
    if (word /= 'CALL::') then
      if (.not. expect_this ('{', .true., .true., 'AFTER "SR_WAKE"', ele, delim, delim_found)) return
    endif
    call parser_read_sr_wake (ele, delim, delim_found, err_flag)
  endif

  return
endif

!-------------------------------

if (attrib_word == 'LR_WAKE') then
  if (.not. expect_this ('=', .true., .true., 'AFTER "LR_WAKE"', ele, delim, delim_found)) return
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)
  ! ele1[lr_wake] = ele2[lr_wake] construct.
  if (delim == '[') then
    ele2 => parser_find_ele_for_attrib_transfer ('LR_WAKE', word); if (err_flag) return
    if (.not. associated(ele%wake)) allocate (ele%wake)
    if (.not. associated(ele2%wake)) then
      call parser_error ('LR_WAKE NOT DEFINED FOR: ' // ele2%name)
      return
    endif
    ele%wake%lr = ele2%wake%lr
  ! "ele1[lr_wake] = call::..." or "ele1: ..., lr_wake = {...}, ..." construct.
  else
    if (word /= 'CALL::') then
      if (.not. expect_this ('{', .true., .true., 'AFTER "LR_WAKE"', ele, delim, delim_found)) return
    endif
    call parser_read_lr_wake (ele, delim, delim_found, err_flag)
  endif

  return
endif

!-------------------------------
! Converter distribution

if (attrib_word == 'DISTRIBUTION') then
  if (.not. expect_this ('=', .true., .true., 'AFTER "CARTESIAN_MAP"', ele, delim, delim_found)) return
  call converter_distribution_parser (ele, delim, delim_found, err_flag)
  return
endif

!-------------------------------
! Cartesian_map field

if (attrib_word == 'CARTESIAN_MAP') then

  if (.not. expect_this ('=', .true., .true., 'AFTER "CARTESIAN_MAP"', ele, delim, delim_found)) return
  call get_next_word (word, ix_word, ':[],(){}', delim, delim_found, call_check = .true.)

  ! "ele1[cartesian_map] = ele2[cartesian_map]" construct

  if (delim == '[') then
    ele2 => parser_find_ele_for_attrib_transfer ('CARTESIAN_MAP', word)
    if (err_flag) return
    if (.not. associated(ele2%cartesian_map)) then
      call parser_error ('NO CARTESIAN_MAP ASSOCIATED WITH LATTICE ELEMENT: ' // word)
      return
    endif
    call transfer_fieldmap(ele2, ele, cartesian_map$)
    return
  endif

  !

  if (associated(ele%cartesian_map)) then
    i_ptr = size(ele%cartesian_map) + 1
    ele0%cartesian_map => ele%cartesian_map
    allocate(ele%cartesian_map(i_ptr))
    do i = 1, i_ptr-1
     ele%cartesian_map(i) = ele0%cartesian_map(i)
    enddo
  else
    allocate(ele%cartesian_map(1))
    i_ptr = 1
  endif

  ! "ele1[cartesian_map] = call::..." or "ele1: ..., cartesian_map = {...}, ..." construct.

  if (.not. expect_this ('{', .true., .true., 'AFTER "CARTESIAN_MAP"', ele, delim, delim_found)) return
  allocate (ele%cartesian_map(i_ptr)%ptr)
  call parse_cartesian_map(ele%cartesian_map(i_ptr), ele, lat, delim, delim_found, err_flag)

  if (ele%key == wiggler$ .or. ele%key == undulator$) ele%field_calc = fieldmap$
  return
endif

!-------------------------------
! Cylindrical_map field

if (attrib_word == 'CYLINDRICAL_MAP') then

  if (.not. expect_this ('=', .true., .true., 'AFTER "CYLINDRICAL_MAP"', ele, delim, delim_found)) return
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)

  ! "ele1[cylindrical_map] = ele2[cylindrical_map]" construct

  if (delim == '[') then
    ele2 => parser_find_ele_for_attrib_transfer ('CYLINDRICAL_MAP', word)
    if (err_flag) return
    if (.not. associated(ele2%cylindrical_map)) then
      call parser_error ('NO CYLINDRICAL_MAP ASSOCIATED WITH LATTICE ELEMENT: ' // word)
      return
    endif
    call transfer_fieldmap(ele2, ele, cylindrical_map$)
    return
  endif

  if (associated(ele%cylindrical_map)) then
    i_ptr = size(ele%cylindrical_map) + 1
    ele0%cylindrical_map => ele%cylindrical_map
    allocate(ele%cylindrical_map(i_ptr))
    do i = 1, i_ptr-1
     ele%cylindrical_map(i) = ele0%cylindrical_map(i)
    enddo
  else
    allocate(ele%cylindrical_map(1))
    i_ptr = 1
  endif

  if (.not. expect_this ('{', .true., .true., 'AFTER "CYLINDRICAL_MAP"', ele, delim, delim_found)) return
  allocate (ele%cylindrical_map(i_ptr)%ptr)
  cl_map => ele%cylindrical_map(i_ptr)
  if (ele%key == lcavity$ .or. ele%key == rfcavity$) cl_map%harmonic = 1 ! Default
  call parse_cylindrical_map(cl_map, ele, lat, delim, delim_found, err_flag)

  if (ele%key == wiggler$ .or. ele%key == undulator$) ele%field_calc = fieldmap$
  return
endif

!-------------------------------
! grid_field field

if (attrib_word == 'GRID_FIELD') then

  ! Note: get_next_word will change "call::" to "hdf5" or "binary" if appropriate.
  if (.not. expect_this ('=', .true., .true., 'AFTER "GRID_FIELD"', ele, delim, delim_found)) return
  call get_next_word (word, ix_word, ':[],(){}', delim, delim_found, call_check = .true.)

  ! "ele1[grid_field] = ele2[grid_field]" construct

  if (delim == '[') then
    ele2 => parser_find_ele_for_attrib_transfer ('GRID_FIELD', word)
    if (err_flag) return
    if (.not. associated(ele2%grid_field)) then
      call parser_error ('NO GRID_FIELD ASSOCIATED WITH LATTICE ELEMENT: ' // word)
      return
    endif
    call transfer_fieldmap(ele2, ele, grid_field$)
    return
  endif

  if (word /= 'hdf5') then
    if (associated(ele%grid_field)) then
      i_ptr = size(ele%grid_field) + 1
      ele0%grid_field => ele%grid_field
      allocate(ele%grid_field(i_ptr))
      do i = 1, i_ptr-1
       ele%grid_field(i) = ele0%grid_field(i)
      enddo
      deallocate (ele0%grid_field)
    else
      allocate(ele%grid_field(1))
      i_ptr = 1
    endif
  endif

  if (word == 'binary') then
    call get_next_word (line, ix, ', ', delim, delim_found, .false.)
    call parser_file_stack('push', line, err = err_flag, open_file = .false.); if (err_flag) return
    call read_binary_grid_field(bp_com%current_file%full_name, ele, ele%grid_field(i_ptr), err_flag)
    call parser_file_stack('pop')
    if (err_flag) then
      call parser_error ('ERROR READING BINARY GRID_FIELD FILE.')
      return
    endif
  elseif (word == 'hdf5') then
    call get_next_word (line, ix, ', ', delim, delim_found, .false.)
    call parser_file_stack('push', line, err = err_flag, open_file = .false.); if (err_flag) return
    call hdf5_read_grid_field(bp_com%current_file%full_name, ele, ele%grid_field, err_flag, combine = .true.)
    call parser_file_stack('pop')
    if (err_flag) then
      call parser_error ('ERROR READING HDF5 GRID_FIELD FILE.')
      return
    endif
  else
    if (.not. expect_this ('{', .true., .true., 'AFTER "GRID_FIELD"', ele, delim, delim_found)) return
    allocate (ele%grid_field(i_ptr)%ptr)
    g_field => ele%grid_field(i_ptr)
    if (ele%key == lcavity$ .or. ele%key == rfcavity$) g_field%harmonic = 1 ! Default

    call parse_grid_field(g_field, ele, lat, delim, delim_found, err_flag)
  endif

  if (ele%key == wiggler$ .or. ele%key == undulator$) ele%field_calc = fieldmap$
  return
endif

!-------------------------------
! Gen_Grad_field field

if (attrib_word == 'GEN_GRAD_MAP') then

  if (.not. expect_this ('=', .true., .true., 'AFTER "GEN_GRAD_MAP"', ele, delim, delim_found)) return
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)

  ! "ele1[gen_grad_map] = ele2[gen_grad_map]" construct

  if (delim == '[') then
    ele2 => parser_find_ele_for_attrib_transfer ('GEN_GRAD_MAP', word)
    if (err_flag) return
    if (.not. associated(ele2%gen_grad_map)) then
      call parser_error ('NO GEN_GRAD_MAP ASSOCIATED WITH LATTICE ELEMENT: ' // word)
      return
    endif
    call transfer_fieldmap(ele2, ele, gen_grad_map$)
    return
  endif

  if (associated(ele%gen_grad_map)) then
    i_ptr = size(ele%gen_grad_map) + 1
    ele0%gen_grad_map => ele%gen_grad_map
    allocate(ele%gen_grad_map(i_ptr))
    allocate(ele%gen_grad_map(i_ptr)%gg(0))
    do i = 1, i_ptr-1
      ele%gen_grad_map(i) = ele0%gen_grad_map(i)
    enddo
    deallocate(ele0%gen_grad_map)
  else
    allocate(ele%gen_grad_map(1))
    allocate(ele%gen_grad_map(1)%gg(0))
    i_ptr = 1
  endif

  if (.not. expect_this ('{', .true., .true., 'AFTER "GEN_GRAD_MAP"', ele, delim, delim_found)) return
  gg_map => ele%gen_grad_map(i_ptr)

  call parse_gen_grad_map(gg_map, ele, lat, delim, delim_found, err_flag)

  if (ele%key == wiggler$ .or. ele%key == undulator$) ele%field_calc = fieldmap$
  return
endif

!------------------------------
! wiggler term attribute

if (ix_attrib == term$ .and. (ele%key == wiggler$ .or. ele%key == undulator$)) then

  err_flag = .true. ! assume the worst

  if (delim /= '(') then   ! ) then
    call parser_error ('"TERM" FOR A WIGGLER NOT FOLLOWED BY A "(" FOR: ' // ele%name)  ! )
    return
  endif

  call parser_get_integer (ix, word, ix_word, delim, delim_found, err_flag, 'BAD WIGGLER "TERM(IX)" CONSTRUCT'); if (err_flag) return

  if (delim /= ')') then
    call parser_error ('CANNOT FIND CLOSING ")" for a "TERM(i)" FOR A WIGGLER"', 'FOR: ' // ele%name)
    return
  endif

  write (str_ix, '(a, i3, a)') 'TERM(', ix, ')'

  if (.not. associated(ele%cartesian_map)) then
    allocate(ele%cartesian_map(1))
    ct_map => ele%cartesian_map(1)
    allocate(ct_map%ptr)
    allocate(ct_map%ptr%term(ix))
    ct_map%ptr%file = bp_com%line2_file_name
    ct_map%master_parameter = polarity$
  else
    ct_map => ele%cartesian_map(1)
    if (ix > size(ct_map%ptr%term)) then
      call move_alloc (ct_map%ptr%term, ct_terms)
      allocate (ct_map%ptr%term(ix))
      ct_map%ptr%term(1:size(ct_terms)) = ct_terms
      deallocate (ct_terms)
    endif
  endif

  ! 1) chop "=", 2) chop to "{", 3) chop to "}", 4) chop to "," if it exists

  call get_next_word (word, ix_word1, ':,={}', delim1, delim_found, .true.) 
  call get_next_word (word, ix_word2, ':,={}', delim2, delim_found, .true., call_check = .true.)  

  if (delim1 /= '=' .or. delim2 /= '{' .or. ix_word1 /= 0 .or. ix_word2 /= 0) then
    call parser_error ('CONFUSED SYNTAX FOR TERM IN WIGGLER: ' // ele%name, str_ix)
    return
  endif

  err_str = trim(ele%name) // ' ' // str_ix
  ct_term => ct_map%ptr%term(ix)

  call parse_evaluate_value (err_str, ct_term%coef, lat, delim, delim_found, err_flag, ',', ele);   if (err_flag) return
  call parse_evaluate_value (err_str, ct_term%kx, lat, delim, delim_found, err_flag, ',', ele);     if (err_flag) return
  call parse_evaluate_value (err_str, ct_term%ky, lat, delim, delim_found, err_flag, ',', ele);     if (err_flag) return
  call parse_evaluate_value (err_str, ct_term%kz, lat, delim, delim_found, err_flag, ',', ele);     if (err_flag) return
  call parse_evaluate_value (err_str, ct_term%phi_z, lat, delim, delim_found, err_flag, ',}', ele); if (err_flag) return

  old_style_input = .true.
  ct_term%family = family_y$

  if (delim == ',') then
    ct_term%x0 = ct_term%phi_z
    call parse_evaluate_value (err_str, ct_term%y0, lat, delim, delim_found, err_flag, ',', ele); if (err_flag) return
    call parse_evaluate_value (err_str, ct_term%phi_z, lat, delim, delim_found, err_flag, ',', ele); if (err_flag) return
    call get_switch ('FAMILY', ['Y ', 'X ', 'QU', 'SQ'], ct_term%family, err_flag, ele, delim, delim_found); if (err_flag) return
    if (.not. expect_this ('}', .true., .false., 'AFTER "FAMILY" SWITCH', ele, delim, delim_found)) return
    old_style_input = .false.
    call parser_error ('"HYBRID" STYLE WIGGLER TERMS DEPRECATED. PLEASE CONVERT TO CARTESIAN_MAP FORM.', level = s_warn$)
  endif

  kx = ct_term%kx
  ky = ct_term%ky
  kz = ct_term%kz
  tol = 1d-5 * (kx**2 + ky**2 + kz**2)

  if (abs(ky**2 - kx**2 - kz**2) < tol) then
    ct_term%form = hyper_y$
    ky = sign_of(ky, .false.) * sqrt(kx**2 + kz**2)

    if (old_style_input) then
      if (ct_term%kx == 0) ct_term%kx = 1d-30  ! Something small to prevent divide by zero problems.
    endif

  elseif (abs(ky**2 + kx**2 - kz**2) < tol) then
    ct_term%form = hyper_xy$
    kz = sign_of(kz, .false.) * sqrt(kx**2 + ky**2)

    if (old_style_input) then
      ct_term%coef = ct_term%coef * ct_term%kz / ct_term%ky
      if (ct_term%kx == 0) ct_term%kx = 1d-30  ! Something small to prevent divide by zero problems.
      if (ct_term%ky == 0) ct_term%ky = 1d-30  ! Something small to prevent divide by zero problems.
    endif

  elseif (abs(ky**2 - kx**2 + kz**2) < tol) then
    ct_term%form = hyper_x$
    kx = sign_of(kx, .false.) * sqrt(ky**2 + kz**2)

    if (old_style_input) then
      ct_term%coef = ct_term%coef * ct_term%kx / ct_term%ky
      if (ct_term%ky == 0) ct_term%ky = 1d-30  ! Something small to prevent divide by zero problems.
    endif

  else
    call parser_error ('WIGGLER TERM DOES NOT HAVE CONSISTANT Kx, Ky, and Kz', &
                  'FOR WIGGLER: ' // ele%name // '  ' // str_ix)
    err_flag = .true.
    return
  endif

  call get_next_word (word, ix_word,  ':,=()', delim,  delim_found, .true.)  
  if (ix_word /= 0) then
    call parser_error ('BAD SYNTAX FOR WIGGLER: ' // ele%name, str_ix)
    err_flag = .true.
    return
  endif

  ele%field_calc = fieldmap$
  return

endif

! Check that next delim is a "=". 
! If not, it might be a flag attribute or an attribute that has a default value.

if (delim /= '=')  then
  err_flag = .false.

  if (ele%key == multipole$ .and. ix_attrib >= t0$ .and. attrib_word(1:1) == 'T') then
    ele%b_pole(ix_attrib-t0$) = pi / (2*(ix_attrib-t0$) + 2)
    return
  endif

  if (attrib_word == 'TILT') then
    select case (ele%key)
    case (quadrupole$, sol_quad$) 
      ele%value(tilt$) = pi / 4.0_rp
      return
    case (sextupole$) 
      ele%value(tilt$) = pi / 6.0_rp
      return
    case (octupole$) 
      ele%value(tilt$) = pi / 8.0_rp
      return
    case default
      call parser_error ('SORRY I''M NOT PROGRAMMED TO USE A "TILT" DEFAULT' // &
                         'FOR A: ' // key_name(ele%key), 'FOR: ' // ele%name)
      err_flag = .true.
      return
    end select
  endif

  if (ele%key == sbend$ .or. ele%key == rbend$) then
    select case (ix_attrib)
    case (fint$)
      ele%value(fint$) = 0.5_rp
      return
    case (fintx$)
      ele%value(fintx$) = 0.5_rp
      return
    end select
  endif

  select case (attrib_word)

  case ('SUPERIMPOSE')
    ele%lord_status = super_lord$
    pele%superposition_has_been_set = .true.

  case default
    call parser_error ('EXPECTING "=" AFTER ATTRIBUTE: ' // word,  'FOR ELEMENT: ' // ele%name)
    err_flag = .true.
  end select

  return
endif

! get the value of the attribute.
! Stuff like TYPE, ALIAS, and DESCRIP attributes are special because their "values"
! are character strings

select case (attrib_word)

case ('REFERENCE')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  call get_next_word(pele%ref_name, ix_word,  '=,', delim, delim_found, .true.)

case ('OFFSET')
  call parse_evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag, ele = ele)
  if (err_flag) return
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  pele%offset = value

case ('FIELD_OVERLAPS')

  ! If pele is not present then bmad_parser2 is the parser and this is an element in the lattice.
  ! In this case, simple call create_field_overlap directly.

  call get_list_of_names (ele, 'FIELD_OVERLAPS', name_list, delim, delim_found, err_flag); if (err_flag) return
  nn = size(name_list)

  if (present(pele)) then
    n = 0
    if (allocated(pele%field_overlaps)) n = size(pele%field_overlaps)
    call re_allocate (pele%field_overlaps, n+nn)
    pele%field_overlaps(n+1:n+nn) = name_list

  else
    do i = 1, n
      call create_field_overlap (ele%branch%lat, ele%name, name_list(i), err_flag)
      if (err_flag) then
        call parser_error ('OVERLAP ELEMENT: ' // name_list(i), 'NOT FOUND FOR OVERLAPPING ELEMENT: ' // ele%name)
      endif
    enddo
  endif

case('TYPE', 'ALIAS', 'DESCRIP', 'SR_WAKE_FILE', 'LR_WAKE_FILE', 'LATTICE', 'TO', 'MACHINE', &
     'TO_LINE', 'TO_ELEMENT', 'CRYSTAL_TYPE', 'MATERIAL_TYPE', 'ORIGIN_ELE', 'PHYSICAL_SOURCE')
  call bmad_parser_string_attribute_set (ele, attrib_word, delim, delim_found, pele = pele)

case ('REF_ORBIT')
  if (.not. parse_real_list (lat, ele%name // ' REF_ORBIT', ele%taylor%ref, .true., delim, delim_found)) return
  if (.not. expect_one_of (', ', .false., ele%name, delim, delim_found)) return

case ('TAYLOR_ORDER')
  call parser_get_integer (ix, word, ix_word, delim, delim_found, err_flag); if (err_flag) return
  if (ix <= 0) then
    call parser_error ('TAYLOR_ORDER IS LESS THAN 1')
    return
  endif
  ptc_private%taylor_order_saved = ix
  lat%input_taylor_order = ix

case ('RUNGE_KUTTA_ORDER')
  call parser_get_integer (ix, word, ix_word, delim, delim_found, err_flag); if (err_flag) return
  if (ix /= 2 .and. ix /= 4) then
    call parser_error ('RUNGE_KUTTA_ORDER NOT EQUAL TO 2 OR 4')
    return
  endif
  bmad_com%runge_kutta_order = ix
  bp_com%extra%runge_kutta_order_set = .true.

case ('SYMPLECTIFY') 
  if (how == def$ .and. (delim == ',' .or. .not. delim_found)) then
    ele%symplectify = .true.
  else
    call parser_get_logical (attrib_word, ele%symplectify, ele%name, delim, delim_found, err_flag); if (err_flag) return
  endif
  
case ('IS_ON')
  call parser_get_logical (attrib_word, ele%is_on, ele%name, delim, delim_found, err_flag)

case ('SUPERIMPOSE')
  call parser_get_logical (attrib_word, logic, ele%name, delim, delim_found, err_flag); if (err_flag) return
  if (logic) then
    ele%lord_status = super_lord$
  else
    ele%lord_status = not_a_lord$
  endif
  pele%superposition_has_been_set = .true.

case ('APERTURE_AT')
  call get_switch (attrib_word, aperture_at_name(1:), ele%aperture_at, err_flag, ele, delim, delim_found); if (err_flag) return

case ('APERTURE_TYPE')
  call get_switch (attrib_word, aperture_type_name(1:), ele%aperture_type, err_flag, ele, delim, delim_found); if (err_flag) return

case ('CAVITY_TYPE')
  call get_switch (attrib_word, cavity_type_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(cavity_type$) = ix

case ('COUPLER_AT')
  call get_switch (attrib_word, end_at_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(coupler_at$) = ix

case ('CREATE_JUMBO_SLAVE')
  if (.not. present(pele)) call parser_error ('INTERNAL ERROR...')
  call parser_get_logical (attrib_word, pele%create_jumbo_slave, ele%name, delim, delim_found, err_flag); if (err_flag) return

case ('CSR_METHOD')
  call get_switch (attrib_word, csr_method_name(1:), switch, err_flag, ele, delim, delim_found)
  if (err_flag) return
  ele%csr_method = switch

case ('DEFAULT_TRACKING_SPECIES')
  call get_next_word (word, ix_word, ':,=(){}', delim, delim_found, .false.)
  ix = species_id(word)
  if (ix == invalid$) then
    call parser_error ('INVALID PARTICLE SPECIES: ' // word)
    return
  endif

  ele%value(default_tracking_species$) = ix
  j = nint(ele%value(ix_branch$)) 
  if (j >= 0) lat%branch(j)%param%default_tracking_species = ix 

case ('ELE_ORIGIN')
  call get_switch (attrib_word, anchor_pt_name(1:), pele%ele_pt, err_flag, ele, delim, delim_found); if (err_flag) return

case ('ENERGY_DISTRIBUTION')
  call get_switch (attrib_word, distribution_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(energy_distribution$) = ix

case ('EXACT_MULTIPOLES')
  call get_switch (attrib_word, exact_multipoles_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(exact_multipoles$) = ix

case ('FIELD_CALC')
  call get_switch (attrib_word, field_calc_name(1:), ele%field_calc, err_flag, ele, delim, delim_found); if (err_flag) return

case ('FIELD_MASTER')
  call parser_get_logical (attrib_word, ele%field_master, ele%name, delim, delim_found, err_flag); if (err_flag) return

case ('FRINGE_AT')
  call get_switch (attrib_word, end_at_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(fringe_at$) = ix

case ('FRINGE_TYPE')
  call get_switch (attrib_word, fringe_type_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  if (.not. valid_fringe_type(ele, ix)) then
    call parser_error ('NOT A VALID FRINGE_TYPE: ' // word, &
                       'FOR: ' // trim(ele%name), 'WHICH IS A: ' // key_name(ele%key))
    return
  endif
  ele%value(fringe_type$) = ix

case ('GEOMETRY')
  call get_switch (attrib_word, geometry_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(geometry$) = ix
  j = nint(ele%value(ix_branch$)) 
  if (j >= 0) lat%branch(j)%param%geometry = ix

case ('INTERPOLATION')
  if (attrib_word == 'spline') then
    call parser_error ('Setting "interpolation = spline" replaced by "interpolation = cubic".', &
                       'Please revise the lattice file.', level = s_warn$)
  endif
  call get_switch (attrib_word, interpolation_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(interpolation$) = ix

case ('LATTICE_TYPE')   ! Old style
  call parser_error ('PARAMETER[LATTICE_TYPE] IS OLD SYNTAX.', &
                     'PLEASE REPLACE WITH PARAMETER[GEOMETRY] = OPEN/CLOSED')
  call get_switch (attrib_word, lattice_type_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(geometry$) = ix

case ('LIVE_BRANCH')
  call get_logical_real (attrib_word, ele%value(live_branch$), err_flag); if (err_flag) return
  j = nint(ele%value(ix_branch$)) 
  if (j >= 0) lat%branch(j)%param%live_branch = is_true(ele%value(live_branch$))

case ('MAT6_CALC_METHOD')
  call get_switch (attrib_word, mat6_calc_method_name(1:), switch, err_flag, ele, delim, delim_found); if (err_flag) return
  if (.not. valid_mat6_calc_method (ele, not_set$, switch)) then
    if (hetero_list) then
      err_flag = .false.
    else
      err_flag = .true.
      call parser_error ('NOT A VALID MAT6_CALC_METHOD: ' // mat6_calc_method_name(switch), &
                         'FOR: ' // trim(ele%name), 'WHICH IS A: ' // key_name(ele%key))
    endif
    return
  endif
  ele%mat6_calc_method = switch

case ('MODE')
  call get_switch (attrib_word, mode_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(mode$) = ix

case ('OFFSET_MOVES_APERTURE')
  call parser_get_logical (attrib_word, ele%offset_moves_aperture, ele%name, delim, delim_found, err_flag); if (err_flag) return

case ('ORIGIN_ELE_REF_PT')
  call get_switch (attrib_word, ref_pt_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(origin_ele_ref_pt$) = ix

case ('PARTICLE')
  call get_next_word (word, ix_word, ':,=(){}', delim, delim_found, .false.)
  ix = species_id(word)
  if (ix == invalid$ .or. ix == ref_particle$ .or. ix == anti_ref_particle$) then
    call parser_error ('INVALID REFERENCE PARTICLE SPECIES: ' // word)
    return
  endif

  ele%ref_species = ix
  if (ele%key == def_parameter$) lat%param%particle = ix 

case ('PHOTON_TYPE')
  call get_switch (attrib_word, photon_type_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  lat%photon_type = ix   ! photon_type has been set.

case ('PTC_FRINGE_GEOMETRY')
  call get_switch (attrib_word, ptc_fringe_geometry_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(ptc_fringe_geometry$) = ix

case ('PTC_INTEGRATION_TYPE')
  call get_switch (attrib_word, ptc_integration_type_name(1:), ele%ptc_integration_type, err_flag, ele, delim, delim_found); if (err_flag) return

case ('PTC_FIELD_GEOMETRY')
  call get_switch (attrib_word, ptc_field_geometry_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(ptc_field_geometry$) = ix

case ('REF_ORIGIN')
  call get_switch (attrib_word, anchor_pt_name(1:), pele%ref_pt, err_flag, ele, delim, delim_found); if (err_flag) return

case ('REF_COORDS')
  call get_switch (attrib_word, ref_coords_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  if (ix == no_end$) then
    call parser_error ('"REF_COORDS = NO_END" NOW SHOULD BE "USER_SETS_LENGTH = T". PLEASE CHANGE.', level = s_warn$)
    ele%value(user_sets_length$) = 1
  else
    ele%value(ref_coords$) = ix
  endif

case ('REF_ORBIT_FOLLOWS')
  call get_switch (attrib_word, ref_orbit_follows_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(ref_orbit_follows$) = ix

case ('SCALE_MULTIPOLES')
  call parser_get_logical (attrib_word, ele%scale_multipoles, ele%name, delim, delim_found, err_flag); if (err_flag) return

case ('SPATIAL_DISTRIBUTION')
  call get_switch (attrib_word, distribution_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(spatial_distribution$) = ix

case ('SPECIES_OUT')
  call get_next_word (word, ix_word, ':,=(){}', delim, delim_found, .false.)
  ix = species_id(word)
  if (ix == invalid$ .or. ix == ref_particle$ .or. ix == anti_ref_particle$) then
    call parser_error ('INVALID SPECIES_OUT: ' // word)
    return
  endif
  ele%converter%species_out = ix

case ('SPECIES_STRONG')
  call get_next_word (word, ix_word, ':,=(){}', delim, delim_found, .false.)
  ix = species_id(word)
  if (ix == invalid$ .or. ix == ref_particle$ .or. ix == anti_ref_particle$) then
    call parser_error ('INVALID SPECIES_STRONG: ' // word)
    return
  endif
  ele%value(species_strong$) = ix

case ('SPIN_TRACKING_METHOD')
  if (attrib_word == 'BMAD_STANDARD') then
    call parser_error ('SPIN_TRACKING_METHOD = BMAD_STANDARD NOW NO LONGER VALID.', &
                     'PLEASE REPLACE WITH SPIN_TRACKING_METHOD = TRACKING.', &
                     'THIS PROGRAM WILL RUN NORMALLY...', level = s_warn$)
    attrib_word = 'TRACKING'
  endif
  call get_switch (attrib_word, spin_tracking_method_name(1:), switch, err_flag, ele, delim, delim_found)
  if (err_flag) return
  if (.not. valid_spin_tracking_method (ele, switch)) then
    if (hetero_list) then
      err_flag = .false.
    else
      call parser_error ('NOT A VALID SPIN_TRACKING_METHOD: ' // word, &
                         'FOR: ' // trim(ele%name), 'WHICH IS A: ' // key_name(ele%key))
    endif
    return
  endif
  ele%spin_tracking_method = switch

case ('TAYLOR_MAP_INCLUDES_OFFSETS')
  call parser_get_logical (attrib_word, ele%taylor_map_includes_offsets, ele%name, delim, delim_found, err_flag); if (err_flag) return

case ('TRACKING_METHOD')
  call get_switch (attrib_word, tracking_method_name(1:), switch, err_flag, ele, delim, delim_found)
  if (err_flag) return
  if (.not. valid_tracking_method (ele, not_set$, switch)) then
    if (hetero_list) then
      err_flag = .false.
    else
      call parser_error ('NOT A VALID TRACKING_METHOD: ' // bp_com%last_word, &
                         'FOR: ' // trim(ele%name), 'WHICH IS A: ' // key_name(ele%key))
    endif
    return
  endif
  ele%tracking_method = switch

case ('SPACE_CHARGE_METHOD')
  call get_switch (attrib_word, space_charge_method_name(1:), switch, err_flag, ele, delim, delim_found)
  if (err_flag) return
  ele%space_charge_method = switch

case ('VELOCITY_DISTRIBUTION')
  call get_switch (attrib_word, distribution_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(velocity_distribution$) = ix

case ('WRAP_SUPERIMPOSE')
  call parser_get_logical (attrib_word, pele%wrap_superimpose, ele%name, delim, delim_found, err_flag); if (err_flag) return


!------------------------------------------------
case default   ! normal attribute

  if (ele%key == def_line$) then
    select case (attrib_word)
    case ('CBAT_11', 'CMAT_12', 'CMAT_21', 'CMAT_22', 'P0C', 'E_TOT', 'ETA_X', 'ETA_Y', &
          'ETAP_X', 'ETAP_Y', 'ALPHA_A', 'ALPHA_B', 'BETA_A', 'BETA_B', 'PHI_A', 'PHI_B')
      ele%value(inherit_from_fork$) = false$
    end select
  endif

  ! attrib_word = "x_limit" for example will generate an error here but this is not a true error.
  call pointer_to_attribute (ele, attrib_word, .true., a_ptr, err_flag, .false.)
  
  select case (attribute_type(attrib_word))
  case (is_logical$)
    if (associated (a_ptr%l)) then
      call parser_get_logical (trim(ele%name) // ' ' // attrib_word, a_ptr%l, ele%name, delim, delim_found, err_flag)
    else
      call get_logical_real (attrib_word, ele%value(ix_attrib), err_flag)
    endif
    if (err_flag) return

  case (is_integer$)
    if (associated (a_ptr%i)) then
      call parser_get_integer (a_ptr%i, word, ix_word, delim, delim_found, err_flag, trim(ele%name) // ' ' // attrib_word)
      call set_flags_for_changed_attribute (ele, a_ptr%i, set_dependent = (bp_com%parser_name == 'bmad_parser2'))
    else
      call parse_evaluate_value (trim(ele%name) // ' ' // word, ele%value(ix_attrib), lat, delim, delim_found, err_flag, ele = ele)
      call set_flags_for_changed_attribute (ele, ele%value(ix_attrib), set_dependent = (bp_com%parser_name == 'bmad_parser2'))
    endif
    if (err_flag) return
    

  case default
    call parse_evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag, ele = ele)
    if (err_flag) return

    ! multipole attribute?
    if (ele%key == hybrid$ .and. is_attribute(ix_attrib, multipole$)) then
      ele%vec0(ix_attrib-a0$) = value
    elseif (ele%key == hybrid$ .and. is_attribute(ix_attrib, elec_multipole$)) then
      i = 1 + (ix_attrib - a0_elec$ - 1) / 6
      j = ix_attrib - a0_elec$ - 6 * (i - 1)
      ele%mat6(i,j) = value
    elseif (is_attribute(ix_attrib, multipole$) .and. attrib_word(1:4) /= 'CURV') then  
        if (.not. associated(ele%a_pole)) call multipole_init (ele, magnetic$)
        if (ix_attrib >= b0$) then
          ele%b_pole(ix_attrib-b0$) = value
        else
          ele%a_pole(ix_attrib-a0$) = value
        endif
    ! Electric multipole attribute
    elseif (is_attribute(ix_attrib, elec_multipole$)) then
        if (.not. associated(ele%a_pole_elec)) call multipole_init (ele, electric$)
        if (ix_attrib >= b0_elec$) then
          ele%b_pole_elec(ix_attrib-b0_elec$) = value
        else
          ele%a_pole_elec(ix_attrib-a0_elec$) = value
        endif
    !
    elseif (attrib_word == 'RAN_SEED') then
      bp_com%extra%ran_seed = nint(value)
      call ran_seed_put (bp_com%extra%ran_seed)  ! init random number generator
    elseif (attrib_word == 'APERTURE') then
      ele%value(x1_limit$) = value
      ele%value(x2_limit$) = value
      ele%value(y1_limit$) = value
      ele%value(y2_limit$) = value
    elseif (attrib_word == 'X_LIMIT') then
      ele%value(x1_limit$) = value
      ele%value(x2_limit$) = value
    elseif (attrib_word == 'Y_LIMIT') then
      ele%value(y1_limit$) = value
      ele%value(y2_limit$) = value
    else
      if (err_flag .or. .not. associated(a_ptr%r)) then
        call parser_error ('BAD ATTRIBUTE: ' // attrib_word, 'FOR ELEMENT: ' // ele%name)
        return
      endif
      a_ptr%r = value
      call set_flags_for_changed_attribute (ele, a_ptr, set_dependent = (bp_com%parser_name == 'bmad_parser2'))

      if (logic_option(.true., set_field_master)) then
        ix = len_trim(attrib_word)
        if (ix > 9 .and. index(attrib_word, '_GRADIENT') == ix-8) ele%field_master = .true.
        if (ix > 6 .and. index(attrib_word, '_FIELD') == ix-5) ele%field_master = .true.
        if (ix > 10 .and. index(attrib_word, '_FIELD_ERR') == ix-9) ele%field_master = .true.
        if (attrib_word(1:3) == 'BL_') ele%field_master = .true.
        if (ele%key == elseparator$ .and. attrib_word == 'VOLTAGE') ele%field_master = .true.
        if (ele%key == elseparator$ .and. attrib_word == 'E_FIELD') ele%field_master = .true.
      endif

      !

      select case (attrib_word)
      case ('CMAT_11', 'CMAT_12', 'CMAT_21', 'CMAT_22')
        coef = 1 - determinant(ele%c_mat)
        if (coef >= 0) ele%gamma_c = sqrt(coef)

      case ('E_TOT')
        if (ele%key == def_parameter$) then
          lat%ele(0)%value(e_tot$) = value
          lat%ele(0)%value(p0c$) = -1
        else
          ele%value(p0c$) = -1
        endif

        branch => pointer_to_branch(ele%name, lat, parameter_is_branch0 = .true.)
        if (associated(branch)) then
          branch%ele(0)%value(e_tot$) = value
          call set_flags_for_changed_attribute (branch%ele(0), branch%ele(0)%value(e_tot$), &
                                                        set_dependent = (bp_com%parser_name == 'bmad_parser2'))
        endif

      case ('ENERGY')    ! Only in def_mad_beam
        lat%ele(0)%value(e_tot$) = 1d9 * value
        lat%ele(0)%value(p0c$) = -1

      case ('PARTICLE')
        if (ele%key == def_mad_beam$) then
          ele2 => lat%ele(ele%ix_ele+1)   ! Points to def_parameter element
          ele2%ref_species = ele%ref_species
        endif

      case ('P0C')
        if (ele%key == def_parameter$) then
          lat%ele(0)%value(p0c$) = value
          lat%ele(0)%value(e_tot$) = -1
        else
          ele%value(e_tot$) = -1
        endif

        branch => pointer_to_branch(ele%name, lat, parameter_is_branch0 = .true.)
        if (associated(branch)) then
          branch%ele(0)%value(p0c$) = value
          call set_flags_for_changed_attribute (branch%ele(0), branch%ele(0)%value(p0c$), &
                                                    set_dependent = (bp_com%parser_name == 'bmad_parser2'))
        endif

      case ('PC')    ! Only in def_mad_beam
        lat%ele(0)%value(p0c$) = 1d9 * value
        ele%value(e_tot$) = -1

      case ('LR_FREQ_SPREAD')
        call randomize_lr_wake_frequencies (ele, set_done)
        if (set_done) call bp_set_ran_status

      case ('N_PART')
        branch => pointer_to_branch(ele%name, lat, parameter_is_branch0 = .true.)
        if (associated(branch)) branch%param%n_part = value

      case ('RF_FREQUENCY')
        if (ele%key == rfcavity$) ele%value(harmon$) = 0

      case ('HARMON')
        ele%value(rf_frequency$) = 0

      end select       ! attrib_word

    endif

  end select  ! attribute_type(attrib_word)

end select

err_flag = .false.

!--------------------------------------------------------
contains

function parser_find_ele_for_attrib_transfer (attribute, word) result (target_ele)

type (ele_struct), pointer :: target_ele
integer n
character(*) attribute, word
character(40) word2

!

nullify(target_ele)

call get_next_word (word2, ix_word, '[],(){}', delim2, delim_found, call_check = .true.)
if (delim2 /= ']' .or. word2 /= attribute) then
  call parser_error ('BAD ' // attribute // ' CONSTRUCT')
  return
endif

if (.not. expect_this (' ', .false., .false., '', ele, delim, delim_found)) return
call lat_ele_locator (word, lat, eles, n, err_flag)

if (err_flag .or. n /= 1) then
  call parser_error ('LATTICE ELEMENT NOT FOUND: ' // word)
  return
endif

target_ele => eles(1)%ele

end function parser_find_ele_for_attrib_transfer

!--------------------------------------------------------
! contains

function attrib_free_problem (attrib_name) result (is_problem)

type (ele_attribute_struct) attrib_info
type (all_pointer_struct) a_ptr

character(*) attrib_name
logical is_problem, is_free

! Attributes may be definitely free, definitely dependent, or may be free or
! dependent depending upon the state of other element parameters.

! If not check_free then at least check if it is a dependent attribute.

is_problem = .false.

attrib_info = attribute_info(ele, attribute_index(ele, attrib_name))
if (attrib_info%state == dependent$) then
  if (.not. hetero_list) then
    call parser_error ('DEPENDENT ATTRIBUTE NOT FREE TO BE SET: ' // attrib_name, 'FOR: ' // ele%name)
  endif
  is_problem = .true.
  return
endif

if (logic_option(.false., check_free)) then
  is_free = attribute_free (ele, attrib_name, .false.)
  if (.not. is_free) then
    call pointer_to_attribute(ele, attrib_name, .true., a_ptr, err_flag, .false.)
    call set_flags_for_changed_attribute (ele, a_ptr%r, .true.)
  endif
endif

end function attrib_free_problem

!--------------------------------------------------------
! contains

subroutine get_logical_real (name, logic_real, err)

character(*) name
real(rp) logic_real
logical this_logical, err

!

call parser_get_logical (name, this_logical, ele%name, delim, delim_found, err)
if (err) return

if (this_logical) then
  logic_real = 1
else
  logic_real = 0
endif

err = .false.

end subroutine get_logical_real

end subroutine parser_set_attribute 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine get_called_file (delim, call_file, err)

implicit none

character(1) delim
character(*) call_file

integer ix_word, ix, n
logical delim_found, finished, err

!

err = .true.

if (delim /= ',')  call parser_error ('"CALL" NOT FOLLOWED BY COMMA', stop_here = .true.)
call get_next_word(call_file, ix_word, ':=,', delim, delim_found, .true.)

if (ix_word == 0) then
  call parser_error ('NOTHING AFTER "CALL"', stop_here = .true.)
  return
elseif (index('FILENAME', call_file(:ix_word)) /= 1) then
  call parser_error ('INVALID "CALL" COMMAND', stop_here = .true.)
  return
elseif (delim /= '=') then
  call parser_error ('NO "=" AFTER "FILENAME"', stop_here = .true.)
  return
endif

call get_next_word(call_file, ix_word, ',', delim, delim_found, .false.)
if (ix_word == 0) then
  call parser_error ('NO FILE NAME SPECIFIED', stop_here = .true.)
  return
endif

if (call_file(1:1) == '"') then
  call_file = call_file(2:)
  ix = index(call_file, '"')
  if (ix == 0 .or. ix /= len_trim(call_file)) then
    call parser_error ('MISSING DOUBLE QUOTE MARK (") FOR CALL STATEMENT', stop_here = .true.)
    return
  endif
  call_file(ix:ix) = ' '
endif

if (call_file(1:1) == "'") then
  call_file = call_file(2:)
  ix = index(call_file, "'")
  if (ix == 0 .or. ix /= len_trim(call_file)) then
    call parser_error ("MISSING SINGLE QUOTE MARK (') FOR CALL STATEMENT", stop_here = .true.)
    return
  endif
  call_file(ix:ix) = ' '
endif

call parser_file_stack ('push', call_file, finished, err) ! err gets set here

end subroutine get_called_file

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine add_this_taylor_term (ele, i_out, coef, expn) 
!
! Subroutine used by bmad_parser and bmad_parser2 to parse the input file.
! This subroutine is not intended for general use.
!-

subroutine add_this_taylor_term (ele, i_out, coef, expn)

implicit none

type (ele_struct), target :: ele
type (taylor_struct), pointer :: taylor

real(rp) coef
integer i, j, i_out, expn(6)

!

if (i_out >= 100) then
  i_out = i_out - 100
  taylor => ele%spin_taylor(i_out) 
else
  if (i_out < 1 .or. i_out > 6) then
    call parser_error ('"OUT" VALUE IN TAYLOR TERM NOT IN RANGE (1 - 6)', &
                  'FOR TAYLOR ELEMENT: ' // ele%name)
    return
  endif
  taylor => ele%taylor(i_out)
endif

call add_taylor_term (taylor, coef, expn, .true.)

end subroutine add_this_taylor_term

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_call_check(word, ix_word, delim, delim_found, err_flag))
!
! Routine to check if there is a "call::XXX" construct in the input stream.
!-

subroutine parser_call_check(word, ix_word, delim, delim_found, err_flag)

implicit none

integer ix, ix_word

logical delim_found
logical, optional :: err_flag

character(*) word, delim
character(6) str
character(20) suffix
character(n_parse_line) line

!


call string_trim(bp_com%parse_line, bp_com%parse_line, ix)
call str_upcase (str, bp_com%parse_line(1:6))
if (str /= 'CALL::') return

bp_com%parse_line = bp_com%parse_line(7:)
call word_read (bp_com%parse_line, ',} ',  line, ix_word, delim, delim_found, bp_com%parse_line)    
ix = str_last_in_set(line, '.')
suffix = ''
if (ix /= 0 .and. ix > len_trim(line)-10) suffix = line(ix:)

if (suffix == '.h5' .or. suffix == '.hdf5') then
  word = 'hdf5'
  bp_com%parse_line = trim(line) // delim // bp_com%parse_line  ! Put line back on parse line.
  return
elseif (suffix == '.bin') then
  word = 'binary'
  bp_com%parse_line = trim(line) // delim // bp_com%parse_line  ! Put line back on parse line.
  return
else
  bp_com%parse_line = delim // bp_com%parse_line  ! Put delim back on parse line.
  call parser_file_stack ('push_inline', line)
  if (bp_com%fatal_error_flag) then
    if (present(err_flag)) err_flag = .true.
    word = ''
    return
  endif
endif

end subroutine parser_call_check

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_next_word (word, ix_word, delim_list, delim, delim_found, upper_case_word, call_check, err_flag)
!
! Subroutine to get the next word from the input stream.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
! Input:
!   word            -- Character(*): Word returned
!   delim_list      -- Character(*): List of valid delimiters
!   upper_case_word -- Logical, optional: if True then convert word to 
!                       upper case. Default is True.
!   call_check      -- Logical, optional: If present and True then check for 'call::<filename>' construct.
!                         Default is False.
!
! Output
!   ix_word     -- Integer: length of word argument
!   delim       -- Character(1): Actual delimiter found
!   delim_found -- Logical: Set true if a delimiter found. A delimiter
!                    may not be found if the end of the line is reached first.
!   err_flag    -- logical, optional: Set True if there is an error. False otherwise.
!-

subroutine get_next_word (word, ix_word, delim_list, delim, delim_found, upper_case_word, call_check, err_flag)

implicit none

integer ix_a, ix_word

character(*) word, delim_list, delim

integer n, ix

logical delim_found, end_of_file
logical, optional :: upper_case_word, call_check, err_flag

character(n_parse_line) line

! Possible inline call...

if (present(err_flag)) err_flag = .false.

if (logic_option(.false., call_check)) then
  call parser_call_check(word, ix_word, delim, delim_found, err_flag)
  if (logic_option(.false.,err_flag) .or. word == 'hdf5' .or. word == 'binary') return
endif

! Check for continuation character and, if found, then load more characters
! into the parse line from the lattice file. 
! If the input is not from a file then skip this.

if (bp_com%input_from_file) then 
  do
    n = len_trim(bp_com%parse_line)
    if (n == 0 .or. n > 90) exit

    if (bp_com%parse_line(n:n) == '&') then
      call load_parse_line('continue', n, end_of_file, err_flag = err_flag); if (logic_option(.false., err_flag)) return

    else if (index(',({[=', bp_com%parse_line(n:n)) /= 0) then
      call load_parse_line('continue', n+2, end_of_file, err_flag = err_flag); if (logic_option(.false., err_flag)) return

    elseif (index(',)}]=', bp_com%next_line_from_file(1:1)) /= 0) then
      call load_parse_line('continue', n+2, end_of_file, err_flag = err_flag); if (logic_option(.false., err_flag)) return

    else
      if (.not. bp_com%inline_call_active) exit
      ! If in an inline called file then make sure the rest of the file is blank and
      ! return to the calling file
      call load_parse_line('continue', n+2, end_of_file, err_flag = err_flag); if (logic_option(.false., err_flag)) return
      if (bp_com%parse_line(n+1:) /= '') then
        call string_trim (bp_com%parse_line(n+1:), line, ix)
        call str_upcase (line(1:10), line(1:10))
        if (line /= 'END_FILE') THEN
          call parser_error ('EXTRA STUFF IN FILE')
          if (present(err_flag)) err_flag = .true.
        endif
      endif
      bp_com%parse_line(n+1:) = ''
    endif

    if (end_of_file) call parser_file_stack ('pop')
  enddo
endif

! Get the first word in bp_com%parse_line

call word_read (bp_com%parse_line, delim_list, word, ix_word, delim, delim_found, bp_com%parse_line)

if (len(word) < ix_word) then
  call parser_error ('BAD WORD: ' // bp_com%parse_line)
  if (present(err_flag)) err_flag = .true.
  ix_word = len(word)
endif

if (present(upper_case_word)) then
  if (upper_case_word) call str_upcase (word, word)
else
  call str_upcase (word, word)
endif

! Note: "var := num" is old-style variable definition syntax.
! If delim is ":" and next char is "=" then use "=" as the delim

if (delim == ':' .and. index(delim_list, '=') /= 0 .and. bp_com%parse_line(1:1) == '=') then
  delim = '='
  bp_com%parse_line = bp_com%parse_line(2:)
endif

bp_com%last_word = word

end subroutine get_next_word

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_file_stack (how, file_name_in, finished, err, open_file)
!
! Subroutine to keep track of the files that are opened for reading.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_file_stack (how, file_name_in, finished, err, open_file)

implicit none

integer i, ix, ios, n, n_file
integer, pointer :: i_level
character(*) how
character(*), optional :: file_name_in
character(200) file_name, basename, file_name2

logical, optional :: finished, err, open_file
logical found_it, is_relative, valid, err_flag

! "Init" means init

i_level => bp_com%i_file_level
if (present(err)) err = .true.

if (how == 'init') then
  i_level = 0
  call fullfilename ('./', bp_com%file(0)%dir)
  if (present(err)) err = .false.
  if (.not. allocated(bp_com%lat_file_names)) allocate(bp_com%lat_file_names(100))
  bp_com%next_line_from_file = str_garbage$
  bp_com%inline_call_active = .false.
  return
endif

! "push" means open a file and put its name on the stack.

if (present(finished)) finished = .false.

select case (how)
case ('push', 'push_inline')

  i_level = i_level + 1    ! number of files currently open
  if (i_level > f_maxx) then
    call parser_error ('CALL NESTING GREATER THAN 20 LEVELS')
    if (global_com%exit_on_error) call err_exit
  endif

  bp_com%file(i_level)%input_line1_saved  = bp_com%input_line1
  bp_com%file(i_level)%input_line2_saved  = bp_com%input_line2
  bp_com%file(i_level)%rest_of_line_saved = bp_com%rest_of_line
  bp_com%input_line1  = ''
  bp_com%input_line2  = ''
  bp_com%rest_of_line = ''

  bp_com%file(i_level)%next_line_from_file_saved = bp_com%next_line_from_file
  bp_com%file(i_level)%ios_next_line_from_file_saved = bp_com%ios_next_line_from_file
  bp_com%next_line_from_file = str_garbage$

  if (how == 'push_inline') then
    bp_com%file(i_level)%parse_line_saved = bp_com%parse_line
    bp_com%file(i_level)%inline_call_active = .true.    
    bp_com%parse_line = '&'
    bp_com%inline_call_active = .true.
  endif

  bp_com%current_file => bp_com%file(i_level)

  if (i_level == 1) then   ! if we are just starting out then init some vars.
    bp_com%num_lat_files = 0           ! total number of files opened
    bp_com%error_flag = .false.  ! set to true on an error
    bp_com%current_file%full_name = ' '
    bp_com%input_line_meaningful = .false.
  endif

  call fullfilename (file_name_in, file_name2, valid)
  if (.not. valid) then
    call parser_error ('MALFORMED FILE NAME: ' // file_name_in, stop_here = .true.)
    if (global_com%exit_on_error) call err_exit
    do i = 1, i_level-1
      close (bp_com%file(i_level)%f_unit)
    enddo
    return
  endif

  ix = splitfilename (file_name2, bp_com%file(i_level)%dir, basename, is_relative)

  if (bp_com%use_local_lat_file) then
    inquire (file = basename, exist = found_it, name = file_name2)
    if (found_it) bp_com%file(i_level)%dir = bp_com%file(0)%dir
  else
    found_it = .false.
  endif

  if (is_relative .and. .not. found_it) then
    call append_subdirectory (trim(bp_com%file(i_level-1)%dir), bp_com%file(i_level)%dir, bp_com%file(i_level)%dir, err_flag)
    if (err_flag) call parser_error ('BAD DIRECTORY SYNTAX FOR: ' // file_name, stop_here = .true.)
    call append_subdirectory (bp_com%file(i_level-1)%dir, file_name2, file_name2, err_flag)
  endif

  inquire (file = file_name2, exist = found_it, name = file_name)

  bp_com%file(i_level)%full_name = file_name
  bp_com%file(i_level)%f_unit = lunget()

  ! Note: open_file will be False when the file is a binary file.

  if (logic_option(.true., open_file)) then
    open (bp_com%file(i_level)%f_unit, file = file_name, status = 'OLD', action = 'READ', iostat = ios)
    if (ios /= 0 .or. .not. found_it) then
      bp_com%current_file => bp_com%file(i_level-1)  ! For warning
      if (i_level == 1)  then !
        call parser_error ('UNABLE TO OPEN FILE: ' // file_name_in, &
                           '(FULL NAME: ' // trim(file_name) // ')', stop_here = .true.)
      else
        call parser_error ('UNABLE TO OPEN FILE: ' // file_name, &
                           'THIS FROM THE LOGICAL FILE NAME: ' // file_name_in, stop_here = .true.)
      endif
      do i = 1, i_level-1
        close (bp_com%file(i_level)%f_unit)
      enddo
      return
    endif
  endif

  bp_com%current_file%i_line = 0

  n = size(bp_com%lat_file_names)
  n_file = bp_com%num_lat_files + 1
  if (n < n_file) call re_allocate (bp_com%lat_file_names, n + 100)
  bp_com%num_lat_files = n_file
  inquire (file = file_name, name = bp_com%lat_file_names(n_file))

  ! Note: The same file may be validly called multiple times if it is an inline file.
  ! EG: A wall file called inline.
  ! Therefore the warning is disabled.

! "pop" means close the current file and pop its name off the stack

case ('pop')
  close (unit = bp_com%current_file%f_unit)
  i_level = i_level - 1
  if (i_level < 0) then
    call parser_error ('BAD "RETURN"')
    return
  elseif (i_level > 0) then
    bp_com%current_file => bp_com%file(i_level)
  else    ! i_level == 0
    if (present(finished)) finished = .true.
  endif

  bp_com%input_line1  = bp_com%file(i_level+1)%input_line1_saved
  bp_com%input_line2  = bp_com%file(i_level+1)%input_line2_saved
  bp_com%rest_of_line = bp_com%file(i_level+1)%rest_of_line_saved
  bp_com%next_line_from_file = bp_com%file(i_level+1)%next_line_from_file_saved
  bp_com%ios_next_line_from_file = bp_com%file(i_level+1)%ios_next_line_from_file_saved

  if (bp_com%inline_call_active) then
    bp_com%parse_line = trim(bp_com%parse_line) // ' ' // bp_com%file(i_level+1)%parse_line_saved
    bp_com%inline_call_active = bp_com%file(i_level+1)%inline_call_active
  endif

  bp_com%inline_call_active = .false.

! Programming error

case default
  call parser_error ('INTERNAL ERROR IN PARSER_FILE_STACK SUBROUTINE!')
  if (global_com%exit_on_error) call err_exit
end select

if (present(err)) err = .false.

end subroutine parser_file_stack

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine load_parse_line (action, ix_start, end_of_file, err_flag) 
!
! Subroutine to load characters from the input file.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
! Input:
!   action        -- character(*): 'continue', 'new_command', or 'init'
!   ix_start      -- integer: Index in bp_com%parse_line string where to append stuff.
!
! Output:
!   end_of_file       -- logical: End of file reached?
!   err_flag          -- logical, optional: Set True if there is an error. False otherwise
!   bp_com%parse_line -- string to append to.
!-

recursive subroutine load_parse_line (action, ix_start, end_of_file, err_flag)

implicit none

integer ix_start, ix, n, ios

character(*) action
character(n_parse_line+20) :: line
character(1) quote_mark
character(1), parameter :: tab = achar(9)

logical :: end_of_file, flush_this
logical, optional :: err_flag

! action = 'init'

if (action == 'init') then
  bp_com%parse_line = ''
  bp_com%rest_of_line = ''
  bp_com%next_line_from_file = str_garbage$
  bp_com%ios_next_line_from_file = 0  
  return
endif

!

end_of_file = .false.
flush_this = .false.
if (present(err_flag)) err_flag = .false.

! If 'new_command' then will need to flush lines that are part of the rest of the current command. 
! This will happen when there has been an error and the entire command was not parsed.

if (action == 'new_command' .and. bp_com%rest_of_line == '') then
  n = len_trim(bp_com%parse_line)
  if (n /= 0) then
    if (index (',+-*/({[=&', bp_com%parse_line(n:n)) /= 0) flush_this = .true.
  endif
endif

!

do
  ! Read a line or use bp_com%rest_of_line if it exists
  if (bp_com%rest_of_line /= '') then
    line = bp_com%rest_of_line
    bp_com%rest_of_line = ''

  else
    if (bp_com%next_line_from_file == str_garbage$) then
      read (bp_com%current_file%f_unit, '(a)', iostat = bp_com%ios_next_line_from_file) bp_com%next_line_from_file
      if (bp_com%ios_next_line_from_file == 0) then
        call string_trim(bp_com%next_line_from_file, bp_com%next_line_from_file, ix)
      else
        bp_com%next_line_from_file = ''
      endif
    endif

    line = bp_com%next_line_from_file
    ios = bp_com%ios_next_line_from_file

    if (ios < 0) then
      end_of_file = .true.
      if (bp_com%parse_line /= '' .and. .not. bp_com%inline_call_active) then
        call parser_error ('FILE ENDED BEFORE PARSING FINISHED', stop_here = .true.)
        if (present(err_flag)) err_flag = .true.
      endif
      bp_com%next_line_from_file = str_garbage$
      bp_com%ios_next_line_from_file = 0  
      return
    elseif (ios > 0) then
      call parser_error ('ERROR READING INPUT LINE.', '[PERHAPS THE LINE HAS TOO MANY CHARACTERS.]', line)
      if (present(err_flag)) err_flag = .true.
      return
    endif
    read (bp_com%current_file%f_unit, '(a)', iostat = bp_com%ios_next_line_from_file) bp_com%next_line_from_file
    if (bp_com%ios_next_line_from_file == 0) then
      call string_trim(bp_com%next_line_from_file, bp_com%next_line_from_file, ix)
    else
      bp_com%next_line_from_file = ''
    endif
    bp_com%current_file%i_line = bp_com%current_file%i_line + 1
  endif

  ! %input_line1 and %input_line2 are for error messages if needed.
  ! Only the input string being parsed is saved in these lines.
  ! 'new_command' action means we are loading a new input string so start from scratch.
  ! 'continue' action means keep the existing input string.

  if (action == 'continue') then
    bp_com%input_line1 = bp_com%input_line2
    bp_com%input_line2 = line
    bp_com%line1_file_name = bp_com%line2_file_name
    write (bp_com%line2_file_name, '(2a, i0)') trim(bp_com%current_file%full_name), ':', bp_com%current_file%i_line
  elseif (action == 'new_command') then
    bp_com%input_line1 = ''
    bp_com%input_line2 = line
    bp_com%line1_file_name = ''
    write (bp_com%line2_file_name, '(2a, i0)') trim(bp_com%current_file%full_name), ':', bp_com%current_file%i_line
  else
    call parser_error ('INTERNAL ERROR #4: CALL HELP')
    if (present(err_flag)) err_flag = .true.
    if (global_com%exit_on_error) call err_exit
  endif

  ! Look for '!' or ';' except in quoted text.

  quote_mark = ''  
  do ix = 1, len_trim(line)
    select case (line(ix:ix))
    case ('"')
      if (ix == 1) then
        quote_mark = '"'
      elseif (line(ix-1:ix-1) /= '\' .and. quote_mark == '"') then !" 
        quote_mark = ''
      else
        quote_mark = '"'
      endif

    case ("'")
      if (ix == 1) then
        quote_mark = "'"
      elseif (line(ix-1:ix-1) /= '\' .and. quote_mark == "'") then !" 
        quote_mark = ''
      else
        quote_mark = "'"
      endif

    ! strip off comments
    case ('!')
      if (quote_mark /= '') cycle
      line = line(:ix-1)
      exit

    ! semi-colon delimiter means that we need to split the line
    ! and save the 2nd piece for the next time around.
    case (';')
      if (quote_mark /= '') cycle
      if (ix == 1) then
        call string_trim(line(ix+1:), bp_com%rest_of_line, n)
        line = ' '
      elseif (ix > 1) then
        call string_trim(line(ix+1:), bp_com%rest_of_line, n)
        line = line(:ix-1)
      else
        bp_com%rest_of_line = ''
      endif
      exit
    end select
  enddo

  ! if the command line is blank then go back for more input
  call string_trim (line, line, ix)
  if (ix /= 0 .or. bp_com%rest_of_line /= '') exit
enddo

! now simply append the line to %parse_line starting at ix_start

call str_substitute (line)
line = adjustl(line)
if (len_trim(line) + ix_start > n_parse_line) then
  call parser_error ('INPUT LINE HAS TOO MANY CHARACTERS:', line)
  if (present(err_flag)) err_flag = .true.
  return
endif

bp_com%parse_line(ix_start:) = line

! Flush this line if needed

if (flush_this) call load_parse_line (action, ix_start, end_of_file)

end subroutine load_parse_line

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function evaluate_array_index (err_flag, delim_list1, word2, delim_list2, delim2) result (this_index)
!
! Function of evaluate the index of an array. Typically the text being parsed looks like:
!      "5) = ..."         or
!      "6).COMP = ..."
!
! Input:
!   delim_list1 -- Character(1): Delimitor after the integer. Normally ')'.
!   delim_list2 -- Character(*): Delimitor list to mark the end of word2. Normally '='.
!
! Output:
!   err_flag    -- Logical: Set True if there is an error. False otherwise.
!   word2       -- Character(*): Word found after delim1. Normally this should be blank.
!   delim2      -- Character(1): Actual delimitor found after word2.
!   this_index  -- Integer: Integer value
!-

function evaluate_array_index (err_flag, delim_list1, word2, delim_list2, delim2) result (this_index)

implicit none

character(*) delim_list1, word2, delim_list2, delim2
integer this_index, ix_word
character(1) delim
character(20) word
logical err_flag, delim_found


! Init

err_flag = .true.
this_index = -1

! Get integer

call get_next_word (word, ix_word, delim_list1, delim, delim_found)
if (.not. delim_found) return

if (.not. is_integer(word)) return
read (word, *) this_index

! Get word after integer

call get_next_word (word2, ix_word, delim_list2, delim2, delim_found)
if (delim_found) err_flag = .false.

end function evaluate_array_index

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function evaluate_logical (word, iostat) result (this_logic)
!
! Function of convert a string into a logical value.
! Accepted possibilities are:
!   .TRUE.  .FALSE. 
!    TRUE    FALSE
!    T       F
!
! Input:
!   word   -- Character(*): Input string.
!
! Output:
!   this_logic -- Logical: Result.
!   iostat     -- Integer: Status: Returns 0 if conversion successful. 
!-

function evaluate_logical (word, iostat) result (this_logic)

implicit none

character(*), intent(in) :: word
character(len(word)+8) wd
logical this_logic
integer, intent(out) :: iostat
integer i

!

iostat = -1
this_logic = .false.  ! To avoid uninit compiler warnings.
if (word == '') return

call str_upcase(wd, word)

do i = 1, len(word)
  if (wd(i:i) == '') cycle

  if (wd(i:i+6) == '.TRUE. ' .or. wd(i:i+4) == 'TRUE ' .or. wd(i:i+1) == 'T ') then
    this_logic = .true.
    iostat = 0
  elseif (wd(i:i+7) == '.FALSE. ' .or. wd(i:i+5) == 'FALSE ' .or. wd(i:i+1) == 'F ') then
    this_logic = .false.
    iostat = 0
  endif
  return
enddo

end function evaluate_logical

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------      
!+
! Subroutine parse_evaluate_value (err_str, value, lat, delim, delim_found, err_flag, end_delims, ele, string_out, string_in)
!
! This routine evaluates as a real number the characters at the beginning of bp_com%parse_line.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
! Input:
!   err_str     -- character(*): String to print as part of error message if there is an error.
!   lat         -- lat_struct: 
!   end_delims  -- character(*), optional: List of delimiters that should be present after section of line used for evaluation.
!   ele         -- ele_struct, optional: Used to evaluat "%[...]" constructs.
!   string_in   -- character(*), optional: If present then use this string as input instead of reading from the lattice file.
!
! Output:
!   value       -- real(rp):
!   delim       -- character(1): Actual delimiter found. Set to blank is no delim found
!   delim_found -- logical: Set False if end-of-line found instead of a delimiter.
!   err_flag    -- logical: Set True if there is an error. False otherwise.
!   string_out  -- character(:), allocatable, optional: If present then parsed expression is returned but not evaluated.
!                   Useful for group and overlay control expressions.
!-

subroutine parse_evaluate_value (err_str, value, lat, delim, delim_found, err_flag, end_delims, ele, string_out, string_in)

use expression_mod

implicit none

type (lat_struct)  lat
type (expression_atom_struct), allocatable :: stk(:)
type (ele_struct), optional, target :: ele

real(rp) value

integer i, ix_word, ix_str, n_parens, n_stk

character(*) err_str
character(*), optional :: end_delims
character(*), optional :: string_in
character(:), allocatable, optional :: string_out
character(1) delim
character(:), allocatable :: str
character(100) err_str2
character(200) word, str_in

logical delim_found, ran_function_pending
logical err_flag, err, call_check

! Get string

if (present(string_in)) str_in = string_in
call_check = .true.
allocate(character(100):: str)
str = ''
ix_str = 0
n_parens = 0
err_flag = .true.

do
  ! Include "+-" as delims to avoid error with sub-expression exceeding 90 characters and with ending "&" continuation char.
  if (present(string_in)) then
    call word_read (str_in, '(),:}+-|', word, ix_word, delim, delim_found, str_in)
  else
    call get_next_word (word, ix_word, '(),:}+-|', delim, delim_found, upper_case_word = .false., call_check = call_check)
  endif
  call_check = .false.
  str = str(1:ix_str) // word
  ix_str = ix_str + ix_word
  if (.not. delim_found) exit

  select case (delim)
  case (',', ')')
    if (n_parens == 0) exit
  case (':', '|', '}')
    exit
  case default
    ! Nothing to do
  end select

  ix_str = ix_str + 1
  str(ix_str:ix_str) = delim
  if (delim == '(') n_parens = n_parens + 1
  if (delim == ')') n_parens = n_parens - 1
  if (n_parens < 0) then
    call parser_error ('CLOSING PARENS ")" FOUND WITHOUT MATCHING OPENING "(" PARENS.', 'FOR: ' // err_str)
    return
  endif

enddo

! Check Parens
if (n_parens /= 0) then
  call parser_error ('OPENING PARENS "(" FOUND WITHOUT MATCHING CLOSING ")" PARENS.', 'FOR: ' // err_str)
  return
endif

! Check that final delim matches.

if (present(end_delims)) then
  if (.not. delim_found) then
    call parser_error ('NO DELIMITOR AFTER STRING: ' // quote(str), 'FOR: ' // err_str)
    return
  elseif (index(end_delims, delim) == 0) then
    call parser_error ('BAD DELIMITOR: ' // quote(delim) // ' AFTER STRING: ' // quote(str), 'FOR: ' // err_str)
    return
  endif
endif

! If the string_out argument is present then just return the string for later evaluation

if (present(string_out)) then
  if (.not. allocated(string_out)) allocate(character(len(str)):: string_out)
  string_out = str
  err_flag = .false.
  return
endif

! If expression is just a number then evaluate and return

if (is_real(str)) then
  read (str, *) value
  err_flag = .false.
  return
endif

! Make a stack

call expression_string_to_stack(str, stk, n_stk,  err, err_str2)
if (err) then
  call parser_error (err_str2, 'FOR: ' // err_str)
  if (err_str2 == 'MALFORMED EXPRESSION') bp_com%parse_line = ''
  return
endif

do i = 1, n_stk
  select case (stk(i)%type)
  case (ran$, ran_gauss$)
    call bp_set_ran_status
  case (variable$)
    call word_to_value (stk(i)%name, lat, stk(i)%value, err, ele)
    if (err) return
  case (species_const$)
    stk(i)%value = species_id(stk(i)%name)
    if (stk(i)%value == invalid$) then
      call parser_error ('INVALID PARTICLE SPECIES: ' // stk(i)%name)
      return
    endif
  end select
enddo

! Evaluate

value = expression_stack_value (stk, err, err_str2)
if (err) then
  call parser_error (err_str2, 'IN EXPRESSION: ' // expression_stack_to_string(stk), 'FOR: ' // err_str)
endif

err_flag = .false.

end subroutine parse_evaluate_value

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine word_to_value (word, lat, value, err_flag, ele)

implicit none

type (lat_struct), target ::  lat
type (ele_struct), optional, target :: ele
type (all_pointer_struct), allocatable :: ptr(:)
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_attribute_struct) attrib_info

integer i, ix1, ix2, ix_word, ios, ix, n_loc, ix_attrib
integer i_const, ixm, ixm2
real(rp) value
real(rp), pointer :: v(:)
character(*) word
character(40) attrib_name, ele_name
character(80) var_name
logical err_flag, err

! Word may be something like "0>>3[k1]".

err_flag = .true.

! See if this is numeric

if (is_real(word, real_num = value)) then
  err_flag = .false.
  return
endif

! If not numeric...

ix_word = len_trim(word)
var_name = upcase(word)

! If word does not have a "[...]" then it must be a variable

ix1 = index(var_name, '[')
if (ix1 == 0) then
  call find_index (var_name, bp_com2%const%name, bp_com2%const%index, bp_com%i_const_tot, i)
  if (i == 0) then
    do i = 1, size(old_style_physical_const_list)
      if (var_name == upcase(old_style_physical_const_list(i)%name)) then
        call parser_error ('DEPRECATED OLD STYLE PHYSICAL CONSTANT: ' // old_style_physical_const_list(i)%name, &
                           'PLEASE REPLACE WITH: '// physical_const_list(i+24)%name, level = s_warn$)
        value = old_style_physical_const_list(i)%value
        err_flag = .false.
        return
      endif
    enddo

    call parser_error ('VARIABLE USED BUT NOT YET DEFINED: ' // word, 'WILL TREAT AS ZERO.', level = s_warn$)
    value = 0
    err_flag = .false.
    ! To prevent multiple error messages define this variable.
    bp_com%i_const_tot = bp_com%i_const_tot + 1
    if (bp_com%i_const_tot > size(bp_com2%const%name)) call reallocate_bp_com_const()
    i_const = bp_com%i_const_tot
    bp_com2%const(i_const)%name = var_name
    bp_com2%const(i_const)%value = 0
    call find_index (var_name, bp_com2%const%name, bp_com2%const%index, i_const-1, ixm, ixm2)
    bp_com2%const(ixm2+1:i_const)%index = bp_com2%const(ixm2:i_const-1)%index
    bp_com2%const(ixm2)%index = i_const

  else
    value = bp_com2%const(i)%value
    err_flag = .false.
  endif

  return
endif

! Here if word does have a "[...]" then is a element attribute

ele_name = var_name(:ix1-1)    ! name of attribute

ix2 = index(var_name, ']')
attrib_name = var_name(ix1+1:ix2-1)

if (attrib_name == 'S' .and. bp_com%parser_name == 'bmad_parser') then
  call parser_error ('"S" ATTRIBUTE CAN ONLY BE USED WITH BMAD_PARSER2')
  return
endif

if (ele_name == '%') then
  if (present(ele)) then
    call re_allocate_eles(eles, 1)
    eles(1)%ele => ele
    err_flag = .false.
  else
    err_flag = .true.
    call parser_error ('"%[...]" CONSTRUCT USED IN EXPRESSION WHERE THERE IS NO ASSOCIATED LATTICE ELEMENT.')
    return
  endif
endif

! Apertures, etc.

select case (attrib_name)
case ('APERTURE', 'X_LIMIT', 'Y_LIMIT')
  if (ele_name /= '%') call lat_ele_locator (ele_name, lat, eles, n_loc, err_flag)

  if (.not. err_flag .and. size(eles) > 0) then
    v => eles(1)%ele%value
    if (attrib_name == 'APERTURE') then
      if (v(x1_limit$) /= v(x2_limit$) .or. v(x1_limit$) /= v(y1_limit$) .or. &
          v(x1_limit$) /= v(y2_limit$)) then
        err_flag = .true.
      else
        value = v(x1_limit$)
      endif
    elseif (attrib_name == 'X_LIMIT') then
      if (v(x1_limit$) /= v(x2_limit$)) then
        err_flag = .true.
      else
        value = v(x1_limit$)
      endif
    elseif (attrib_name == 'Y_LIMIT') then
      if (v(y1_limit$) /= v(y2_limit$)) then
        err_flag = .true.
      else
        value = v(y1_limit$)
      endif
    endif
  endif

  ! If size(eles) > 1 then there must be more than one element of the same name.

  if (.not. err_flag .and. size(eles) > 1) then
    do i = 2, size(eles)
      if (eles(i)%ele%value(x1_limit$) == eles(1)%ele%value(x1_limit$) .and. &
          eles(i)%ele%value(y1_limit$) == eles(1)%ele%value(y1_limit$)) cycle
      call parser_error (&
            'MULTIPLE ELEMENTS OF THE SAME NAME BUT WITH DIFFERENT ATTRIBUTE VALUES REFERENCED IN: ' // word)
      err_flag = .true.
      exit
    enddo
  endif

  return

! Everything else

case default
  if (ele_name == '%') then
    call re_allocate (ptr, 1)
    call pointer_to_attribute (ele, attrib_name, .false., ptr(1), err, .false., ix_attrib)
  else
    call pointers_to_attribute (lat, ele_name, attrib_name, .false., ptr, err, .false., eles, ix_attrib)
  endif

  if (err .or. size(ptr) == 0) then
    call parser_error('BAD ATTRIBUTE: ' // word)
    return
  elseif (associated(ptr(1)%r)) then
    value = ptr(1)%r
  elseif (associated(ptr(1)%i)) then
    value = ptr(1)%i
  else  ! Must
    call parser_error('ATTRIBUTE IS NOT REAL OR INTEGER: ' // word)
    return
  endif

  ! If this is bmad_parser, and not bmad_parser2, then dependent attributes have not been set and cannot
  ! be used.

  if (ix_attrib > 0 .and. associated(eles(1)%ele) .and. bp_com%parser_name == 'bmad_parser') then
    attrib_info = attribute_info (eles(1)%ele, ix_attrib)
    if (attrib_info%state == dependent$) then
      call parser_error ('DEPENDENT ATTRIBUTE IS NOT CALCULATED BEFORE LATTICE EXPANSION AND', &
                         'THEREFORE CANNOT BE USED BEFORE ANY EXPAND_LATTICE COMMAND: ' // word)
      return
    endif
  endif

  ! If size(ptr) > 1 then there must be more than one element of the same name.

  if (size(ptr) > 1) then
    do i = 2, size(ptr)
      if (ptr(i)%r /= ptr(1)%r) then
        call parser_error (&
              'MULTIPLE ELEMENTS OF THE SAME NAME BUT WITH DIFFERENT ATTRIBUTE VALUES REFERENCED IN: ' // word)
        return
      endif
    enddo
  endif

end select

err_flag = .false.

end subroutine word_to_value

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_add_constant (word, lat, redef_is_error)

implicit none

type (lat_struct) lat
type (expression_atom_struct), allocatable :: temp_const(:)
character(*) word
character(1) delim
real(rp) old_val
integer i, i_const, ixm, ixm2, n
logical delim_found, err_flag, redef_is_error

!

call find_index (word, bp_com2%const%name, bp_com2%const%index, bp_com%i_const_tot, i)
if (i /= 0) then
  old_val = bp_com2%const(i)%value
  call parse_evaluate_value (word, bp_com2%const(i)%value, lat, delim, delim_found, err_flag)

  if (redef_is_error) then
    if (bp_com2%const(i)%value == old_val) then
      call parser_error ('CONSTANTS ARE NOT ALLOWED TO BE REDEFINED: ' // word, &
                         'BUT SINCE OLD_VALUE = NEW_VALUE THIS IS ONLY A WARNING...', &
                         'USE "REDEF:" CONSTRUCT TO GET AROUND THIS (BUT NOT RECOMMENDED).', level = s_warn$)
    else
      call parser_error ('CONSTANTS ARE NOT ALLOWED TO BE REDEFINED: ' // word, &
                         'USE "REDEF:" CONSTRUCT TO GET AROUND THIS (BUT NOT RECOMMENDED).')
    endif
  endif

else
  bp_com%i_const_tot = bp_com%i_const_tot + 1
  if (bp_com%i_const_tot > size(bp_com2%const%name)) call reallocate_bp_com_const()
  i_const = bp_com%i_const_tot
  bp_com2%const(i_const)%name = word
  ! Reindex.
  call find_index (word, bp_com2%const%name, bp_com2%const%index, i_const-1, ixm, ixm2)
  bp_com2%const(ixm2+1:i_const)%index = bp_com2%const(ixm2:i_const-1)%index
  bp_com2%const(ixm2)%index = i_const
  ! Evaluate
  call parse_evaluate_value (bp_com2%const(i_const)%name, bp_com2%const(i_const)%value, lat, delim, delim_found, err_flag)
  ! Put in lat%constant(:) array
  n = 0
  if (allocated(lat%constant)) n = size(lat%constant)
  call move_alloc(lat%constant, temp_const)
  allocate (lat%constant(n+1))
  if (allocated(temp_const)) lat%constant(1:n) = temp_const
  lat%constant(n+1)%name = bp_com2%const(i_const)%name
  lat%constant(n+1)%value = bp_com2%const(i_const)%value
endif

!

if (.not. verify_valid_name(word, len_trim(word), .true.)) then
  call parser_error ('MALFORMED CONSTANT NAME: ' // word)
  return
endif

if (delim_found .and. .not. err_flag) call parser_error  &
                  ('EXTRA CHARACTERS ON RHS: ' // bp_com%parse_line,  &
                   'FOR CONSTANT: ' // bp_com2%const(i_const)%name)

end subroutine parser_add_constant

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine bmad_parser_string_attribute_set (ele, attrib_name, delim, delim_found, pele, str_out)

implicit none

type (ele_struct)  ele
type (parser_ele_struct), optional :: pele

integer ix, ix_word, n

character(*) attrib_name
character(*), optional :: str_out
character(40)  word
character(1)   delim, str_end
character(200) type_name

logical delim_found, err

!

call string_trim(bp_com%parse_line, bp_com%parse_line, ix)

str_end = bp_com%parse_line(1:1)
if (str_end == '"' .or. str_end == "'") then
  bp_com%parse_line = bp_com%parse_line(2:)
  ix = index(bp_com%parse_line, str_end)
  if (ix == 0) then
    call parser_error ('MISSING ENDING QUOTE MARK FOR: ' // attrib_name,  &
                        'FOR ELEMENT: ' // ele%name)
    type_name = ' '
  else
    type_name = bp_com%parse_line(1:ix-1)
    bp_com%parse_line = bp_com%parse_line(ix+1:)
    call get_next_word (word, ix_word, ',=', delim, delim_found, .true.)
    if (ix_word /= 0) call parser_error ( &
              'EXTRA CHARACTERS FOUND AFTER VALUE OF: ' // attrib_name, &
              'I DO NOT UNDERSTAND: ' // word,  &
              'FOR ELEMENT: ' // ele%name)
  endif
else
  call get_next_word (type_name, ix_word, ',= ', delim, delim_found, .false.)
endif

!--------------

select case (attrib_name)
case ('ALIAS')
  ele%alias = type_name
case ('CRYSTAL_TYPE', 'MATERIAL_TYPE', 'PHYSICAL_SOURCE')
  ele%component_name = type_name
case ('DESCRIP')
  if (.not. associated(ele%descrip)) allocate (ele%descrip) 
  ele%descrip = type_name
case ('LATTICE')
  ele%branch%lat%lattice = type_name
case ('MACHINE')
  ele%branch%lat%machine = type_name
case ('LR_WAKE_FILE') 
  call parser_read_old_format_lr_wake (ele, type_name)
case ('MATERIAL', 'CLEAR_MATERIAL', 'OPAQUE_MATERIAL')
  str_out = type_name
case ('TO', 'TO_LINE', 'ORIGIN_ELE')
  ele%component_name = type_name
  call upcase_string (ele%component_name)
case ('TO_ELEMENT')
  pele%ele_name = type_name
  call upcase_string (pele%ele_name)
case ('SR_WAKE_FILE') 
  call parser_read_old_format_sr_wake (ele, type_name)
case ('TYPE')
  ele%type = type_name
case default
  call parser_error ('INTERNAL ERROR IN BMAD_PARSER_STRING_ATTRIBUTE_SET: I NEED HELP!')
  if (global_com%exit_on_error) call err_exit
end select

end subroutine bmad_parser_string_attribute_set

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_read_sr_wake (ele, delim, delim_found, err_flag)
!
! Subroutine to read in a short-range wake field from an external file.
! This subroutine is used by bmad_parser and bmad_parser2.
!
! Input:
!   ele          -- Ele_struct: Element containing wake structure.
!
! Output:
!   ele -- Ele_struct: Element with wake information.
!     %wake%sr      -- Short-range wake potential.
!-
      
subroutine parser_read_sr_wake (ele, delim, delim_found, err_flag)

implicit none

type (ele_struct) ele
type (lat_struct), pointer :: lat
type (wake_sr_mode_struct), target :: trans(100), long(100)
type (wake_sr_mode_struct), pointer :: srm
type (wake_sr_struct), pointer :: wake_sr

integer itrans, ilong, ix_word, which

logical delim_found, err_flag

character(40) err_str, word, attrib_name
character(1) delim

! Init

if (.not. associated(ele%wake)) allocate (ele%wake)
if (.not. allocated(ele%wake%lr%mode)) allocate (ele%wake%lr%mode(0))
if (allocated(ele%wake%sr%long))  deallocate (ele%wake%sr%long)
if (allocated(ele%wake%sr%trans)) deallocate (ele%wake%sr%trans)

lat => ele%branch%lat
wake_sr => ele%wake%sr
trans = wake_sr_mode_struct()
long = wake_sr_mode_struct()
err_flag = .true.

! get data

itrans = 0
ilong = 0

do
  call get_next_word (attrib_name, ix_word, '{}=,()', delim, delim_found, call_check = .true.)
  if (.not. expect_this ('=', .true., .false., 'IN SR_WAKE DEFINITION', ele, delim, delim_found)) return

  select case (attrib_name)
  case ('Z_MAX')
    call parse_evaluate_value (err_str, wake_sr%z_max, lat, delim, delim_found, err_flag, ',}', ele);  if (err_flag) return
    if (delim == '}') exit
    cycle
  case ('Z_SCALE')
    call parse_evaluate_value (err_str, wake_sr%z_scale, lat, delim, delim_found, err_flag, ',}', ele);  if (err_flag) return
    if (delim == '}') exit
    cycle
  case ('AMP_SCALE')
    call parse_evaluate_value (err_str, wake_sr%amp_scale, lat, delim, delim_found, err_flag, ',}', ele);  if (err_flag) return
    if (delim == '}') exit
    cycle
  case ('SCALE_WITH_LENGTH')
    call parser_get_logical (attrib_name, wake_sr%scale_with_length, ele%name, delim, delim_found, err_flag);  if (err_flag) return
    if (.not. expect_one_of (',}', .true., ele%name, delim, delim_found)) return
    if (delim == '}') exit
    cycle
  case ('LONGITUDINAL')
    ilong = ilong + 1
    srm => long(ilong)
    which = longitudinal_mode$
  case ('TRANSVERSE')
    itrans = itrans + 1
    srm => trans(itrans)
    which = -1  ! Not longitudinal. That is, transverse.
  case default
    call parser_error ('UNKNOWN SR_WAKE COMPONENT: ' // attrib_name, 'FOR ELEMENT: ' // ele%name)
    return
  end select

  if (.not. expect_this ('{', .false., .false., 'AFTER "' // trim(attrib_name) // ' =" IN SR_WAKE DEFINITION', ele, delim, delim_found)) return

  err_str = trim(ele%name) // ' SR_WAKE ' // attrib_name
  call parse_evaluate_value (err_str, srm%amp, lat, delim, delim_found, err_flag, ',', ele);  if (err_flag) return
  call parse_evaluate_value (err_str, srm%damp, lat, delim, delim_found, err_flag, ',', ele);  if (err_flag) return
  call parse_evaluate_value (err_str, srm%k, lat, delim, delim_found, err_flag, ',', ele);  if (err_flag) return
  call parse_evaluate_value (err_str, srm%phi, lat, delim, delim_found, err_flag, ',', ele);  if (err_flag) return
  if (which == longitudinal_mode$) then
    call get_switch ('POSITION_DEPENDENCE', sr_longitudinal_position_dep_name, srm%position_dependence, err_flag, ele, delim, delim_found)
  else
    call get_switch ('POLARIZATION', sr_transverse_polarization_name, srm%polarization, err_flag, ele, delim, delim_found)
    if (.not. expect_one_of (',', .true., ele%name, delim, delim_found)) return
    call get_switch ('POSITION_DEPENDENCE', sr_transverse_position_dep_name, srm%position_dependence, err_flag, ele, delim, delim_found)
  endif

  if (.not. expect_one_of ('}', .true., ele%name, delim, delim_found)) return
  if (.not. expect_one_of (',}', .false., ele%name, delim, delim_found)) return
  if (delim == '}') exit
enddo

!

if (.not. expect_one_of (', ', .false., ele%name, delim, delim_found)) return

allocate (ele%wake%sr%long(ilong))
ele%wake%sr%long = long(1:ilong)

allocate (ele%wake%sr%trans(itrans))
ele%wake%sr%trans = trans(1:itrans)

err_flag = .false.

end subroutine parser_read_sr_wake

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_read_lr_wake (ele, delim, delim_found, err_flag)
!
! Subroutine to read in a long-range wake field from an external file.
! This subroutine is used by bmad_parser and bmad_parser2.
!
! Input:
!   ele          -- Ele_struct: Element containing wake structure.
!
! Output:
!   ele          -- Ele_struct: Element with wake information.
!     %wake%lr     -- Long-range wake.
!-
      
subroutine parser_read_lr_wake (ele, delim, delim_found, err_flag)

implicit none

type (ele_struct) ele
type (lat_struct), pointer :: lat
type (wake_lr_mode_struct), target :: lr_mode(100)
type (wake_lr_mode_struct), pointer :: lrm
type (wake_lr_struct), pointer :: wake_lr

integer iterm, ix_word

logical delim_found, err_flag, set_done

character(40) err_str, word, attrib_name
character(1) delim

! Init

if (.not. associated(ele%wake)) allocate (ele%wake)
if (.not. allocated(ele%wake%sr%long))  allocate (ele%wake%sr%long(0))
if (.not. allocated(ele%wake%sr%trans)) allocate (ele%wake%sr%trans(0))
if (allocated(ele%wake%lr%mode)) deallocate (ele%wake%lr%mode)

lat => ele%branch%lat
wake_lr => ele%wake%lr
lr_mode = wake_lr_mode_struct()
err_flag = .true.

! get data

iterm = 0

do
  call get_next_word (attrib_name, ix_word, '{}=,()', delim, delim_found, call_check = .true.)
  if (.not. expect_this ('=', .true., .false., 'IN LR_WAKE DEFINITION', ele, delim, delim_found)) return

  select case (attrib_name)
  case ('AMP_SCALE')
    call parse_evaluate_value (err_str, wake_lr%amp_scale, lat, delim, delim_found, err_flag, ',}', ele);  if (err_flag) return
    if (delim == '}') exit
    cycle
  case ('TIME_SCALE')
    call parse_evaluate_value (err_str, wake_lr%time_scale, lat, delim, delim_found, err_flag, ',}', ele);  if (err_flag) return
    if (delim == '}') exit
    cycle
  case ('FREQ_SPREAD')
    call parse_evaluate_value (err_str, wake_lr%freq_spread, lat, delim, delim_found, err_flag, ',}', ele);  if (err_flag) return
    if (delim == '}') exit
    cycle
  case ('T_REF')
    call parse_evaluate_value (err_str, wake_lr%t_ref, lat, delim, delim_found, err_flag, ',}', ele);  if (err_flag) return
    if (delim == '}') exit
    cycle
  case ('SELF_WAKE_ON')
    call parser_get_logical (attrib_name, wake_lr%self_wake_on, ele%name, delim, delim_found, err_flag);  if (err_flag) return
    if (.not. expect_one_of (',}', .true., ele%name, delim, delim_found)) return
    if (delim == '}') exit
    cycle
  case ('MODE')
  case default
    call parser_error ('UNKNOWN LR_WAKE COMPONENT: ' // attrib_name, 'FOR ELEMENT: ' // ele%name)
    return
  end select

  if (.not. expect_this ('{', .false., .false., 'AFTER "MODE =" IN LR_WAKE DEFINITION', ele, delim, delim_found)) return
  iterm = iterm + 1
  lrm => lr_mode(iterm)

  err_str = trim(ele%name) // ' LR_WAKE MODE'
  call parse_evaluate_value (err_str, lrm%freq_in, lat, delim, delim_found, err_flag, ',', ele);  if (err_flag) return
  lrm%freq = lrm%freq_in
  call parse_evaluate_value (err_str, lrm%r_over_q, lat, delim, delim_found, err_flag, ',', ele);  if (err_flag) return
  call parse_evaluate_value (err_str, lrm%damp, lat, delim, delim_found, err_flag, ',', ele);  if (err_flag) return
  call parse_evaluate_value (err_str, lrm%phi, lat, delim, delim_found, err_flag, ',', ele);  if (err_flag) return

  call parser_get_integer (lrm%m, word, ix_word, delim, delim_found, err_flag, 'BAD LR_WAKE M MODE VALUE')

  if (.not. expect_this (',', .true., .false., 'AFTER M MODE VALUE', ele, delim, delim_found)) return
  call get_next_word (attrib_name, ix_word, ',{}=', delim, delim_found)
  if (index('UNPOLARIZED', trim(upcase(attrib_name))) == 1 .and. attrib_name /= '') then
    lrm%polarized = .false.
  else
    lrm%polarized = .true.
    bp_com%parse_line = trim(attrib_name) // delim // bp_com%parse_line
    call parse_evaluate_value (err_str, lrm%angle, lat, delim, delim_found, err_flag, ',}', ele);  if (err_flag) return
  endif

  if (delim == '}') then
    if (.not. expect_one_of (',}', .false., ele%name, delim, delim_found)) return
    if (delim == '}') exit
    cycle
  endif

  call parse_evaluate_value (err_str, lrm%b_sin, lat, delim, delim_found, err_flag, ',', ele);  if (err_flag) return
  call parse_evaluate_value (err_str, lrm%b_cos, lat, delim, delim_found, err_flag, ',', ele);  if (err_flag) return
  call parse_evaluate_value (err_str, lrm%a_sin, lat, delim, delim_found, err_flag, ',', ele);  if (err_flag) return
  call parse_evaluate_value (err_str, lrm%a_cos, lat, delim, delim_found, err_flag, '}', ele);  if (err_flag) return

  if (.not. expect_one_of (',}', .false., ele%name, delim, delim_found)) return
  if (delim == '}') exit
enddo

!

if (.not. expect_one_of (', ', .false., ele%name, delim, delim_found)) return

allocate (wake_lr%mode(iterm))
wake_lr%mode = lr_mode(1:iterm)

if (wake_lr%freq_spread /= 0) then
  call randomize_lr_wake_frequencies (ele, set_done)
  if (set_done) call bp_set_ran_status
endif

err_flag = .false.

end subroutine parser_read_lr_wake

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_read_old_format_lr_wake (ele, lr_file_name)
!
! Subroutine to read in a long-range wake field from an external file.
! This subroutine is used by bmad_parser and bmad_parser2.
!
! Input:
!   ele          -- Ele_struct: Element containing wake structure.
!   lr_file_name -- Character(*):  Name of long-range wake field file.
!
! Output:
!   ele -- Ele_struct: Element with wake information.
!     %wake%lr%mode(:)       -- Long-range wake potential.
!-
      
subroutine parser_read_old_format_lr_wake (ele, lr_file_name)

implicit none

type lr_wake_input_struct
  real(rp) :: freq  = real_garbage$ ! Actual Frequency in Hz.
  real(rp) :: R_over_Q = 0          ! Strength in V/C/m^(2*m_mode).
  real(rp) :: Q = 0                 ! Quality factor.
  integer :: m = 0                  ! Mode order (1 = dipole, 2 = quad, etc.)
  character(16) :: angle = ''       ! polarization angle (radians/2pi).
  real(rp) :: b_sin = 0, b_cos = 0, a_sin = 0, a_cos = 0, t_ref = 0
end type

type (ele_struct) ele
type (lr_wake_input_struct) lr(500)
integer n_row, iu, i, j, ios
character(*) lr_file_name
character(200) full_file_name
logical set_done, finished, err

namelist / long_range_modes / lr

! Init

if (.not. associated(ele%wake)) allocate (ele%wake)
if (.not. allocated(ele%wake%sr%long))  allocate (ele%wake%sr%long(0))
if (.not. allocated(ele%wake%sr%trans)) allocate (ele%wake%sr%trans(0))
if (allocated(ele%wake%lr%mode)) deallocate (ele%wake%lr%mode)

! get data

call parser_file_stack ('push', lr_file_name, finished, err)
if (err) return
iu = bp_com%current_file%f_unit

ele%wake%lr%file = lr_file_name

lr = lr_wake_input_struct()
read (iu, nml = long_range_modes, iostat = ios)
call parser_file_stack ('pop')

if (ios > 0 .or. lr(1)%freq == -1) then
  call parser_error ('CANNOT READ LONG_RANGE_MODES NAMELIST FOR ELEMENT: ' // ele%name, & 
                'FROM FILE: '// full_file_name)
  return
endif

! Transfer info to ele structure.

n_row = count(lr%freq /= real_garbage$)
allocate (ele%wake%lr%mode(n_row))
j = 0
do i = 1, size(lr)
  if (lr(i)%freq == real_garbage$) cycle
  if (lr(i)%freq == 0) lr(i)%freq = -1

  j = j + 1
  ele%wake%lr%mode(j)%freq_in   = lr(i)%freq
  ele%wake%lr%mode(j)%freq      = lr(i)%freq
  ele%wake%lr%mode(j)%r_over_q  = lr(i)%r_over_q
  ele%wake%lr%mode(j)%phi       = 0
  ele%wake%lr%mode(j)%Q         = lr(i)%Q
  ele%wake%lr%mode(j)%m         = lr(i)%m
  ele%wake%lr%mode(j)%b_sin     = lr(i)%b_sin
  ele%wake%lr%mode(j)%b_cos     = lr(i)%b_cos
  ele%wake%lr%mode(j)%a_sin     = lr(i)%a_sin
  ele%wake%lr%mode(j)%a_cos     = lr(i)%a_cos

  call downcase_string(lr(i)%angle)
  if (lr(i)%angle == '') then
    call parser_error ('LONG_RANGE_MODE ANGLE IS MISSING. MUST BE NUMBER OR "UNPOLARIZED"', & 
                  'FOR ELEMENT: ' // ele%name, &
                  'IN FILE: ' // full_file_name)
    cycle
  endif

  if (index('unpolarized', trim(lr(j)%angle)) == 1) then
    ele%wake%lr%mode(j)%polarized = .false.
    ele%wake%lr%mode(j)%angle     = 0
  else
    ele%wake%lr%mode(j)%polarized = .true.
    read (lr(j)%angle, *, iostat = ios) ele%wake%lr%mode(j)%angle
    if (ios /= 0) then
      call parser_error ('BAD LONG_RANGE_MODE ANGLE.', &
                    'FOR ELEMENT: ' // ele%name, &
                    'IN FILE: ' // full_file_name)
      cycle
    endif
  endif
enddo

call randomize_lr_wake_frequencies (ele, set_done)
if (set_done) call bp_set_ran_status

end subroutine parser_read_old_format_lr_wake

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_read_old_format_sr_wake (ele, sr_file_name)
!
! Subroutine to read in a short-range wake field from an external file.
! This subroutine is used by bmad_parser and bmad_parser2.
!
! Input:
!   ele          -- Ele_struct: Element containing wake structure.
!   sr_file_name -- Character(*):  Name of short-range wake field file.
!
! Output:
!   ele -- Ele_struct: Element with wake information.
!     %wake%sr%table(:)       -- Short-range wake potential.
!     %wake%sr%long(:)  -- Short-range wake potential.
!     %wake%sr%trans(:) -- Short-range wake potential.
!-

subroutine parser_read_old_format_sr_wake (ele, sr_file_name)

implicit none

type (ele_struct) ele
type (wake_sr_mode_struct), target :: longitudinal(100), transverse(100)
type (wake_sr_mode_struct), pointer :: sr(:), sr1

real(rp) z_max
integer n, j, iu, ios, ix, i, ixx

character(*) sr_file_name
character(140) line, line_in
character(200) full_file_name

logical in_namelist, finished, err

character(16), parameter :: old_sr_position_dependence_name(3) = [character(16):: 'none', 'linear_leading', 'linear_trailing']

! init

if (.not. associated(ele%wake))   allocate (ele%wake)
if (.not. allocated(ele%wake%lr%mode)) allocate (ele%wake%lr%mode(0))
if (allocated(ele%wake%sr%long))  deallocate (ele%wake%sr%long)
if (allocated(ele%wake%sr%trans)) deallocate (ele%wake%sr%trans)

! Open file

ele%wake%sr%file = sr_file_name
call parser_file_stack ('push', sr_file_name, finished, err)
if (err) return
iu = bp_com%current_file%f_unit

! Get data

longitudinal = wake_sr_mode_struct()
longitudinal%phi = real_garbage$

transverse = wake_sr_mode_struct()
transverse%phi = real_garbage$

z_max = real_garbage$
in_namelist = .false.

do
  read (iu, '(a)', iostat = ios) line_in
  if (ios /= 0) then
    call parser_error ('END OF FILE REACHED BEFORE SHORT_RANGE_MODES NAMELIST PARSED.', &
                       'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
    return
  endif

  line = line_in
  ix = index(line, '!')
  if (ix /= 0) line = line(1:ix-1)
  call string_trim(line, line, ixx)
  if (ixx == 0) cycle
  call downcase_string (line) 

  if (in_namelist) then
    if (line(1:1) == '/') exit
  else
    if (line == '&short_range_modes') then
      in_namelist = .true.
      cycle
    endif
  endif

  if (line(1:12) == 'longitudinal') then
    call string_trim(line(13:), line, ixx)
    sr => longitudinal
    if (.not. get_this_sr1 (sr, sr1)) return
    sr1%position_dependence = none$    ! Default
    sr1%polarization = none$             ! Default

    if (.not. expect_equal_sign()) return
    if (.not. get_this_param (sr1%amp)) return
    if (.not. get_this_param (sr1%damp)) return
    if (.not. get_this_param (sr1%k)) return
    if (.not. get_this_param (sr1%phi)) return
    if (.not. get_this_switch (sr1%polarization, sr_transverse_polarization_name)) return
    if (.not. get_this_switch (sr1%position_dependence, old_sr_position_dependence_name)) return
    if (.not. expect_nothing ()) return

  elseif (line(1:10) == 'transverse') then
    call string_trim(line(11:), line, ixx)
    sr => transverse
    if (.not. get_this_sr1 (sr, sr1)) return
    sr1%position_dependence = leading$     ! Default

    if (.not. expect_equal_sign()) return
    if (.not. get_this_param (sr1%amp)) return
    if (.not. get_this_param (sr1%damp)) return
    if (.not. get_this_param (sr1%k)) return
    if (.not. get_this_param (sr1%phi)) return
    if (.not. get_this_switch (sr1%polarization, sr_transverse_polarization_name)) return
    if (.not. get_this_switch (sr1%position_dependence, old_sr_position_dependence_name)) return
    if (.not. expect_nothing ()) return

  elseif (line(1:ixx) == 'z_max') then
    call string_trim(line(ixx+1:), line, ixx)
    if (.not. expect_equal_sign()) return
    if (.not. get_this_param (z_max)) return
    if (.not. expect_nothing ()) return

  else
    call parser_error ('BAD PARAMETER NAME IN SHORT_RANGE_MODES NAMELIST: ' // line_in, &
                       '(MUST BE "LONGITUDINAL", "TRANSVERSE", OR "Z_MAX")', &
                       'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
    return
  endif
enddo

call parser_file_stack ('pop')

! Put data in element

n = count(longitudinal%phi /= real_garbage$)
allocate (ele%wake%sr%long(n))
ele%wake%sr%long = longitudinal(1:n)
if (any(longitudinal(1:n)%phi == real_garbage$)) call parser_error ( &
    'JUMBLED INDEX FOR LONGITUDINAL SHORT_RANGE_MODES FROM FILE: ' // full_file_name, &
    'FOR ELEMENT: ' // ele%name)

do i = 1, n
  sr1 => ele%wake%sr%long(i)
  if (sr1%polarization == none$ .and. sr1%position_dependence == none$) then
    sr1%position_dependence = none$
  elseif (sr1%polarization == x_polarization$ .and. sr1%position_dependence == leading$) then
    sr1%position_dependence = x_leading$
  elseif (sr1%polarization == x_polarization$ .and. sr1%position_dependence == trailing$) then
    sr1%position_dependence = x_trailing$
  elseif (sr1%polarization == y_polarization$ .and. sr1%position_dependence == leading$) then
    sr1%position_dependence = y_leading$
  elseif (sr1%polarization == y_polarization$ .and. sr1%position_dependence == trailing$) then
    sr1%position_dependence = y_trailing$
  else
    call parser_error ('MISMATCHED POLARIZATION AND POSITION_DEPENDENCE FOR SHORT-RANG WAKE LONGITUDINAL MODE.', &
                       'FOR MODE WITH INDEX \i0\ OF ELEMENT: ' // ele%name, &
                       'FROM FILE: ' // full_file_name, i_array = [i])
  endif
enddo

n = count(transverse%phi /= real_garbage$)
allocate (ele%wake%sr%trans(n))
ele%wake%sr%trans = transverse(1:n)
if (any(transverse(1:n)%phi == real_garbage$)) call parser_error ( &
    'JUMBLED INDEX FOR TRANSVERSE SHORT_RANGE_MODES FROM FILE: ' // full_file_name, &
    'FOR ELEMENT: ' // ele%name)

ele%wake%sr%long%phi  = ele%wake%sr%long%phi / twopi
ele%wake%sr%trans%phi = ele%wake%sr%trans%phi / twopi

ele%wake%sr%z_max = z_max
if (z_max == real_garbage$) call parser_error ( &
    'Z_MAX NOT SET FOR SHORT_RANGE_MODES FROM FILE: ' // full_file_name, &
    'FOR ELEMENT: ' // ele%name)

!-------------------------------------------------------------------------
contains

function expect_equal_sign () result (is_ok)

logical is_ok

is_ok = .false.

if (line(1:1) /= '=') then
  call parser_error ('MISING "=" SIGN IN SHORT_RANGE_MODES NAMELIST: ' // line_in, &
                     'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

call string_trim(line(2:), line, ixx)
is_ok = .true.

end function expect_equal_sign

!-------------------------------------------------------------------------
! contains

function get_this_sr1 (sr, sr1) result (is_ok)

type (wake_sr_mode_struct), pointer :: sr(:), sr1
integer ixp, ios, n
logical is_ok


!
is_ok = .false.

if (line(1:1) /= '(') then
  call parser_error ('MISING "(" IN SHORT_RANGE_MODES NAMELIST: ' // line_in, &
                     'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

ixp = index(line, ')') 
if (ixp == 0) then
  call parser_error ('MISING "(" IN SHORT_RANGE_MODES NAMELIST: ' // line_in, &
                     'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

read(line(2:ixp-1), *, iostat = ios) n
if (ios /= 0 .or. n < 1 .or. n > size(sr)) then
  call parser_error ('BAD WAKE MODE INDEX IN SHORT_RANGE_MODES NAMELIST: ' // line_in, &
                     'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

call string_trim (line(ixp+1:), line, ixx)
sr1 => sr(n)

is_ok = .true.

end function get_this_sr1

!-------------------------------------------------------------------------
! contains

function get_this_param (param) result (is_ok)

real(rp) param
integer ios
logical is_ok

!

is_ok = .false.

if (ixx == 0) then
  call parser_error ('MISING NUMBER IN SHORT_RANGE_MODES NAMELIST: ' // line_in, &
                     'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

read (line(1:ixx), *, iostat = ios) param

if (ios /= 0) then
  call parser_error ('MALFORMED NUMBER IN SHORT_RANGE_MODES NAMELIST: ' // line_in, &
                     'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

call string_trim(line(ixx+1:), line, ixx)

is_ok = .true.

end function get_this_param

!-------------------------------------------------------------------------
! contains

function get_this_switch (switch, names) result (is_ok)

integer switch, ix
character(*) :: names(:)
logical is_ok

is_ok = .false.

if (ixx == 0) then
  is_ok = .true.
  return
endif

call match_word (line(1:ixx), names, switch)

if (switch < 1) then
  call parser_error ('BAD SWITCH NAME IN SHORT_RANGE_MODES NAMELIST: ' // line_in, &
                     'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif

call string_trim(line(ixx+1:), line, ixx)

is_ok = .true.

end function get_this_switch

!-------------------------------------------------------------------------
! contains

function expect_nothing () result (is_ok)

logical is_ok

is_ok = .false.

if (line /= '') then    
  call parser_error ('EXTRA STUFF ON LINE: ' // line_in, &
                     'FROM FILE: ' // full_file_name, 'FOR ELEMENT: ' // ele%name)
  return
endif


is_ok = .true.

end function expect_nothing

end subroutine parser_read_old_format_sr_wake

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_list_of_names (ele, err_str, delim, delim_found)
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-
      
subroutine get_list_of_names (ele, err_str, name_list, delim, delim_found, err_flag)

implicit none

type (ele_struct)  ele

character(*) err_str
character(*), allocatable :: name_list(:)
character(1) delim
character(40) word

integer ix_word, n_name 

logical delim_found, curly_parens_found, err_flag

! Opening "}" is optional if there is only one word

err_flag = .true.

call get_next_word (word, ix_word, '{,}()', delim, delim_found, .true.)
curly_parens_found = (delim == '{')
if (curly_parens_found .and. ix_word /= 0) then
  call parser_error ('ERROR PARSING: ' // err_str, 'FOUND STUFF BEFORE OPENING "{"')
  return
endif

if (curly_parens_found) call get_next_word (word, ix_word, '{,}()', delim, delim_found, .true.)

n_name = 0
call re_allocate (name_list, 10)

do
  n_name = n_name + 1
  if (n_name > size(name_list)) call re_allocate(name_list, n_name+10)
  name_list(n_name) = word

  if ((delim == '}' .and. .not. curly_parens_found) .or. (.not. delim_found .and. curly_parens_found)) then
    call parser_error ('ERROR PARSING: ' // err_str, 'MISMATCHED {} BRACES')
    return
  endif

  if (delim_found .and. delim /= '}' .and. delim /= ',') then
    call parser_error ('ERROR PARSING: ' // err_str, 'MALFORMED STATEMENT')
    return
  endif

  if (delim == '}' .or. .not. curly_parens_found) then 
    if (delim == '}') call get_next_word (word, ix_word, ',=:', delim, delim_found, .true.)
    exit
  endif

  call get_next_word (word, ix_word, '{,}()', delim, delim_found, .true.)
enddo

call re_allocate(name_list, n_name)
err_flag = .false.

end subroutine get_list_of_names

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_overlay_group_names (ele, lat, pele, delim, delim_found, is_control_var_list)
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
! Input:
!   is_control_var_list   -- logical: If True then parsing "var = {...}" list.
!                                     If False then parsing "group/overlay/girder = {...}" list.
!-
      
subroutine get_overlay_group_names (ele, lat, pele, delim, delim_found, is_control_var_list, err_flag)

implicit none

type my_knot_struct
  real(rp), allocatable :: y(:)
end type

type (ele_struct)  ele
type (parser_ele_struct), target :: pele
type (lat_struct)  lat
type (parser_controller_struct), pointer :: pc
type (var_length_string_struct), allocatable :: expression(:), etemp(:)
type (my_knot_struct), allocatable :: knot(:), ktemp(:)

real(rp) value

integer ix_word, ix, n_slave, i, j, k, n, i_last
                           
character(1) delim
character(40) word_in, word
character(40), allocatable :: name(:), attrib_name(:)
character(200) err_str

logical :: is_control_var_list, err_flag
logical delim_found, err, end_of_file, ele_names_only

!

err_flag = .true.
allocate (name(40), attrib_name(40), expression(40), knot(40))

call get_next_word (word_in, ix_word, '{,}', delim, delim_found, .true.)
if (delim /= '{' .or. ix_word /= 0) then
  call parser_error  ('BAD SYNTAX FOR ' // trim(upcase(control_name(ele%lord_status))), &
                      'EXPECTING A "{" AFTER "=")', 'FOR ELEMENT: ' // ele%name)
  return
endif

ele_names_only = (is_control_var_list .or. ele%key == girder$)

! loop over all names in "{...}" list

n_slave = 0
do 
  call get_next_word (word_in, ix_word, '{,}/:=', delim, delim_found, .true.)
  if (bp_com%parse_line(1:1) == ':') then   ! Element is something like "rfcavity::*" in which case the break is in the wrong place
    bp_com%parse_line = bp_com%parse_line(2:)
    call get_next_word (word, ix, '{,}/:=', delim, delim_found, .true.)
    word_in = trim(word_in) // '::' // word
    ix_word = ix_word + 2 + ix
  endif

  if (delim == ':' .and. ele%key == girder$) pele%is_range = .true.

  ! A variable list may be empty for rampers using expressions with ran() or ran_gauss().
  ! In this case, create a variable named NULL.

  if (word_in == '' .and. delim == '}' .and. n_slave == 0 .and. is_control_var_list) word_in = 'NULL'

  n_slave = n_slave + 1
  n = size(name)
  if (n_slave > n) then
    call re_allocate(name, 2*n_slave)
    call re_allocate(attrib_name, 2*n_slave)
    call move_alloc(expression, etemp)
    allocate (expression(2*n_slave))
    expression(1:n) = etemp
    call move_alloc(knot, ktemp)
    allocate(knot(2*n_slave))
    knot(1:n) = ktemp
  endif

  word = word_in
  allocate(character(1):: expression(n_slave)%str)
  expression(n_slave)%str = ''

  j = index(word, '[')
  if (j > 0) then
    k = index(word, ']')
    if (k <= j+1) then
      call parser_error ('BAD ATTRIBUTE SYNTAX: ' // word_in, 'FOR: ' // ele%name)
      return
    elseif (word(k+1:) /= '') then
      call parser_error ('MALFORMED SLAVE NAME: ' // word_in, 'FOR: ' // ele%name)
      return
    endif
    attrib_name(n_slave) = word(j+1:k-1)
    word = word(:j-1)
  else
    attrib_name(n_slave) = blank_name$
  endif

  name(n_slave) = word

  if (word == '') then
    if (is_control_var_list) then
      call parser_error ('VARIABLE NAME MISSING WHEN PARSING LORD: ' // ele%name)
      return
    else
      !! call parser_error ('No slave elements defined for lord: ' // ele%name, level = s_warn$)
      n_slave = n_slave - 1
    endif
  endif

  ! If ele_names_only = True then evaluating "var = {...}" construct or is a girder.
  ! In this case, there are no expressions

  bp_com%parse_line = adjustl(bp_com%parse_line)

  if (delim == '/' .or. (delim == ':' .and. ele%key /= girder$)) then

    ! Not expecting this delim if ele_names_only
    if (ele_names_only) then
      call parser_error ('BAD VAR = {...} CONSTRUCT.', 'FOR: ' // ele%name)
      return


    ! If a list of knot y-values
    elseif (bp_com%parse_line(1:1) == '{') then   ! y_knot array
      if (.not. parse_real_list2 (lat, 'ERROR PARSING Y KNOT POINTS FOR: ' // ele%name, knot(n_slave)%y, n, delim, delim_found, 10, '{', ',', '}')) return
      if (.not. expect_one_of ('},', .false., ele%name, delim, delim_found)) return
      call re_allocate(knot(n_slave)%y, n)

    ! Parse expression
    else
      ! Just parse and not evaluate.
      call parse_evaluate_value (trim(ele%name), value, lat, delim, delim_found, err, ',}', ele, expression(n_slave)%str)
      if (err) then
        call parser_error ('BAD EXPRESSION: ' // word_in,  'FOR ELEMENT: ' // ele%name)
        return
      endif
    endif
  endif

  if (delim == '}') then
    call get_next_word (word, ix_word, ',=:', delim, delim_found, .true.)
    exit
  elseif (delim /= ',' .and. .not. (delim == ':' .and. ele%key == girder$)) then
    call parser_error ('MALFORMED DEFINITION OF SLAVE: ' // word_in, 'FOR CONTROLLER: ' // ele%name)
    return
  endif
                        
enddo

!

if (is_control_var_list) then
  ! Note: ele%control is already allocated
  allocate(ele%control%var(n_slave))
  ele%control%var%name = name(1:n_slave)
  err_flag = .false.
  return
endif

!

if (ele%lord_status /= girder_lord$) allocate(ele%control)
allocate (pele%control(n_slave))
pele%control%name = name(1:n_slave)
pele%control%attrib_name = attrib_name(1:n_slave)
i_last = -1

if (ele%lord_status /= girder_lord$) then
  do i = 1, n_slave
    pc => pele%control(i)

    if (expression(i)%str == '' .and. .not. allocated(knot(i)%y)) then  ! Use default which is derived from last slave.
      if (i_last == -1) then
        expression(i)%str = '1'
      else
        knot(i)%y = knot(i_last)%y
      endif
    endif

    if (allocated(knot(i)%y)) then
      i_last = i
      pc%y_knot = knot(i)%y
    else        ! Expression string to stack
      i_last = -1
      call expression_string_to_stack (expression(i)%str, pc%stack, pc%n_stk, err, err_str)
      if (err) then
        call parser_error (err_str, 'FOR ELEMENT: ' // ele%name, 'EXPRESSION: ' // trim(expression(i)%str))
        return
      endif
    endif
  enddo
endif

err_flag = .false.

end subroutine get_overlay_group_names

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function verify_valid_name (name, ix_name, pure_name) result (is_valid)
!
! Routine to check if a name is well formed. Examples:
!   "0>>Q0"                           -- Invalid
!   "Q1##1"                           -- Invalid (double hash not accepted).
!   "Q2A_C.\7#"                       -- Pure name (no "[", "]", "(", ")", "%" characters present).
!   "Q3[GRID_FIELD(1)%FIELD_SCALE]"   -- Valid but not a pure name.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
!
! Input:
!   name        -- character(*): Name(1:ix_name) is the string to check.
!   ix_name     -- integer: Number of characters in the name.
!   pure_name   -- logical, optional: If tu
!
! Output:
!   is_valid    -- logical: True if name is well formed. False otherwise.   
!-

function verify_valid_name (name, ix_name, pure_name) result (is_valid)

implicit none

integer i, ix_name, ix1, ix2

character(*) name
character(*), parameter :: letters = '\ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
character(*), parameter :: valid_chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ\0123456789_[]().#%'
character(*), parameter :: pure_chars  = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ\0123456789_.#'
character(1), parameter :: tab = achar(9)

logical, optional :: pure_name
logical OK, is_valid

! Check for blank spaces

is_valid = .false.

do i = 1, min(ix_name, len(name))
  if (name(i:i) == ' ' .or. name(i:i) == tab) then
    call parser_error  ('NO DELIMITER BETWEEN NAMES: ' // name)
    return
  endif
enddo

! Check for name too long

if (ix_name > len(name)) then
   call parser_error ('NAME TOO LONG: ' // name)
   ix_name = len(name)      ! chop name
  return
endif

! Check for name too short

if (ix_name == 0) then
  call parser_error ('BLANK NAME')
  return
endif

! Check for invalid characters in name.
! Example valid: "Q1[GRID_FIELD(1)%FIELD_SCALE]"

OK = .true.
if (index(letters, name(1:1)) == 0) OK = .false.

if (logic_option(.false., pure_name)) then
  do i = 1, ix_name
    if (index(pure_chars, name(i:i)) == 0) OK = .false.
  enddo
else
  do i = 1, ix_name
    if (index(valid_chars, name(i:i)) == 0) OK = .false.
  enddo
endif

if (.not. OK) then
  call parser_error ('INVALID NAME: UNRECOGNIZED CHARACTERS IN: ' // name)
  return
endif

! Check for non-matched "(" ")" pairs

ix1 = index(name, '(')
ix2 = index(name, ')')
if (ix1 /= 0 .or. ix2 /= 0) then
  if (ix1 == 0) then
    call parser_error ('UNMATCHED PARENTHESIS: ' // name)
    return
  endif
  if (ix2 <= ix1+1) then
    call parser_error  ('INVALID: REVERSED PARENTHESES: ' // name)
    return
  endif
  if (index(name(ix1+1:), '(') /= 0 .or. index(name(ix2+1:), ')') /=  0) then
    call parser_error ('INVALID: BAD PARENTHESES: ' // name)
    return
  endif
endif

! Check for non matched "[" "]" pairs

ix1 = index(name, '[')
ix2 = index(name, ']')

if (ix1 /= 0 .or. ix2 /= 0) then
  if (ix1 == 0) then 
    call parser_error ('UNMATCHED BRACKET: ' // name)
    return
  endif

  if (ix2 <= ix1+1) then
    call parser_error  ('INVALID: REVERSED BRACKETS: ' // name)
    return
  endif

  if (index(name(ix1+1:), '[') /= 0 .or. index(name(ix2+1:), ']') /=  0) then
    call parser_error ('INVALID: BAD BRACKETS: ' // name)
    return
  endif

  if (ix2 /= len(name)) then
    if (name(ix2+1:ix2+1) /= ' ') then
      call parser_error  ('INVALID: SOMETHING AFTER CLOSING "]" BRACKET: ' // name)
      return
    endif
  endif

endif

! check for more than 40 characters

if ((ix1 == 0 .and. ix_name > 40) .or. (ix1 > 41 .or. ix2 - ix1 > 41)) then
  call parser_error ('NAME HAS > 40 CHARACTERS: ' // name)
  return
endif

is_valid = .true.

end function verify_valid_name

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_error (what1, what2, what3, what4, seq, pele, stop_here, level, r_array, i_array)
!
! Routine to print an error message generated when parsing a lattice.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
! Input:
!   what1       -- character(*): First line in error message.
!   what2       -- character(*), optional: Second line in error message.
!   what3       -- character(*), optional: Third line in error message.
!   what4       -- character(*), optional: Fourth line in error message.
!   seq         -- seq_struct, optional: Used when error is generated during reading of a lattice file.
!                     Contains information such as file name, and line number where error was detected.
!   pele        -- parser_ele_struct, optional: Used when error is associated with a lattice element.
!                     Contains information on the lattice element.
!   stop_here   -- logical, optional: If present and True then immediately stop. Used with severe errors.
!   level       -- integer, optional: Possibilities are:
!                     s_error$   -- Parser error (default). At end parser error argument will be set True.
!                     s_warn$    -- Warning level.
!                     s_info$    -- Informational message to be printed.
!   r_array(:)  -- real(rp), optional: Real numbers to be encoded in error message. See out_io doc.
!   i_array(:)  -- integer, optional: Integer numbers to be encoded in error message. See out_io doc.
!-

subroutine parser_error (what1, what2, what3, what4, seq, pele, stop_here, level, r_array, i_array)

implicit none

type (seq_struct), optional :: seq
type (parser_ele_struct), optional :: pele

character(*) what1
character(*), optional :: what2, what3, what4
character(160) lines(12)
character(16), parameter :: r_name = 'parser_error'
real(rp), optional :: r_array(:)
integer, optional :: level
integer nl, err_level
integer, optional :: i_array(:)
logical, optional :: stop_here

! bp_com%error_flag is a common logical used so program will stop at end of parsing

err_level = integer_option(s_error$, level)

if (bp_com%print_err) then

  nl = 0

  if (bp_com%parser_name == '') then
    select case (err_level)
    case (s_info$)
      nl=nl+1; lines(nl) = 'Note: ' // trim(what1)
    case (s_warn$)
      nl=nl+1; lines(nl) = 'WARNING: ' // trim(what1)
    case (s_error$)
      nl=nl+1; lines(nl) = 'ERROR: ' // trim(what1)
    end select

  else
    select case (err_level)
    case (s_info$)
      nl=nl+1; lines(nl) = 'Note from: ' // trim(bp_com%parser_name) // ': ' // trim(what1)
    case (s_warn$)
      nl=nl+1; lines(nl) = 'WARNING IN ' // trim(bp_com%parser_name) // ': ' // trim(what1)
    case (s_error$)
      nl=nl+1; lines(nl) = 'ERROR IN ' // trim(bp_com%parser_name) // ': ' // trim(what1)
    end select
  endif

  if (present(what2)) then
    nl=nl+1; lines(nl) = '     ' // trim(what2)
    if (lines(nl) == '') nl = nl - 1 
  endif

  if (present(what3)) then
    nl=nl+1; lines(nl) = '     ' // trim(what3)
    if (lines(nl) == '') nl = nl - 1 
  endif

  if (present(what4)) then
    nl=nl+1; lines(nl) = '     ' // trim(what4)
    if (lines(nl) == '') nl = nl - 1 
  endif

  ! If bp_com%num_lat_files = 0 then no parser init has been done

  if (bp_com%num_lat_files /= 0) then
    if (present(seq)) then
      nl=nl+1; lines(nl) = '      IN FILE: ' // trim(seq%file_name)
      nl=nl+1; write (lines(nl), '(a, i0)') '      AT LINE: ', seq%ix_file_line
    elseif (bp_com%current_file%full_name /= '') then
      if (bp_com%input_line_meaningful) then
        nl=nl+1; lines(nl) = '      IN FILE: ' // trim(bp_com%current_file%full_name)
        nl=nl+1; write (lines(nl), '(a, i0)') '      AT OR BEFORE LINE: ', bp_com%current_file%i_line
      else
        if (bp_com%i_file_level > 1) then
          nl=nl+1; lines(nl) = '     IN FILE: ' // trim(bp_com%file(bp_com%i_file_level)%full_name)
        endif
        nl=nl+1; lines(nl) = '     ROOT FILE: ' // trim(bp_com%file(1)%full_name)
      endif
    endif

    if (bp_com%input_line_meaningful) then
      if (len_trim(bp_com%input_line1) /= 0) then
        nl=nl+1; lines(nl) = '     ' // trim(bp_com%input_line1)
      endif
      if (len_trim(bp_com%input_line2) /= 0) then
        nl=nl+1; lines(nl) = '     ' // trim(bp_com%input_line2)
      endif
    endif

    if (present(pele)) then
      nl=nl+1; lines(nl) = '      ELEMENT DEFINED IN FILE: ' // trim(pele%lat_file)
      nl=nl+1; write (lines(nl), '(a, i0)') '      AT LINE: ', pele%ix_line_in_file
    endif
  endif

  nl=nl+1; lines(nl) = ''

  call out_io (err_level, r_name, lines(1:nl), r_array = r_array, i_array = i_array, insert_tag_line = .false.)

endif

! Warnings do not result in bp_com%error_flag being set. Just no digested file is generated.

if (err_level == s_warn$) then
  bp_com%write_digested = .false.
elseif (err_level == s_error$) then
  bp_com%error_flag = .true.
  if (logic_option(.false., stop_here)) then
    if (global_com%exit_on_error) stop
    bp_com%fatal_error_flag = .true.
  endif
endif

end subroutine parser_error

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine add_this_multipass (lat, m_slaves, lord_in)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: slave, lord, slave2, lord2, ele
type (ele_struct), optional :: lord_in
type (branch_struct), pointer :: branch
type (lat_ele_loc_struct) m_slaves(:)

integer i, j, k, n, i1, ix, ixc, ixic, ix_lord, ixb, ixb2, ix_n
integer n_multipass, ic, ix_l1, ix_l0, ix_pass, n_links, lmax
integer, allocatable :: indx(:)

character(40) base_name
character(40), allocatable :: names(:)
character(100) slave2_name

! Make sure that all slaves are of the same lord/slave status

slave => pointer_to_ele (lat, m_slaves(1))
do i = 2, size(m_slaves)
  slave2 => pointer_to_ele (lat, m_slaves(i))

  if (slave2%lord_status /= slave%lord_status) then
    ele => slave
    if (slave2%lord_status == super_lord$) ele => slave2
    if (ele%lord_status == super_lord$) then
      call parser_error ('A SUPERPOSITION IN A MULTIPASS REGION HAS CAUSED A CONFLICT NEAR ELEMENT: ' // ele%name)
    else
      call parser_error ('MULTIPASS SETUP ERROR. PLEASE REPORT.')
    endif
    return
  endif

  if (slave2%slave_status /= slave%slave_status) then
    ele => slave
    if (slave2%slave_status == super_slave$) ele => slave2
    if (ele%slave_status == super_slave$) then
      ele => pointer_to_lord(ele, 1)
      call parser_error ('A SUPERPOSITION IN A MULTIPASS REGION HAS CAUSED A CONFLICT NEAR ELEMENT: ' // ele%name)
    else
      call parser_error ('MULTIPASS SETUP ERROR. PLEASE REPORT.')
    endif
    return
  endif
enddo

! Count slaves.
! If n_multi > lat%n_ele_track we are looking at cloning a super_lord which should not happen.
! If n_multi = 1 then, to symplify the lattice, do not create a lord

n_multipass = size(m_slaves)
if (n_multipass == 1) return

! setup multipass_lord

call new_control (lat, ix_lord)
lord => lat%ele(ix_lord)

if (present(lord_in)) then
  lord = lord_in   ! Use lord_in as template
else
  lord = pointer_to_ele(lat, m_slaves(1))  ! Set attributes equal to first slave.
endif

call set_ele_name(lord, lord%name)
lord%logic = .false.  ! So parser_add_superimpose will not try to use as ref ele.
lord%lord_status = multipass_lord$
lord%n_slave = 0
lord%ix1_slave = 0
call add_lattice_control_structs (lord, n_add_slave = n_multipass)

! Multipass_lord does not have reference energy or s_position or map bookkeeping. 

lord%bookkeeping_state%ref_energy = ok$   
lord%bookkeeping_state%s_position = ok$   
lord%bookkeeping_state%mat6       = ok$   

! A multipass lord defaults to multipass_ref_energy = first_pass if neither multipass_ref_energy, p0c and e_tot are not set.

if (lord%value(p0c$) == 0 .and. lord%value(e_tot$) == 0) then
  lord%value(multipass_ref_energy$) = first_pass$
else
  lord%value(multipass_ref_energy$) = user_set$ 
endif

! Setup bookkeeping between lord and slaves

do i = 1, n_multipass
  slave => pointer_to_ele (lat, m_slaves(i))
  ixc = i + lord%ix1_slave - 1
  lat%control(ixc)%lord%ix_ele = ix_lord
  lat%control(ixc)%slave = lat_ele_loc_struct(slave%ix_ele, slave%ix_branch)
  if (slave%n_lord /= 0) then
    call parser_error ('INTERNAL ERROR: CONFUSED MULTIPASS SETUP.', 'PLEASE GET EXPERT HELP!')
    if (global_com%exit_on_error) call err_exit
  endif
  call set_ele_name (slave, trim(slave%name) // '\' // int_str(i))   ! '
  call add_lattice_control_structs (slave, n_add_lord = 1)
  slave%slave_status = multipass_slave$
  ixic = slave%ic1_lord
  lat%ic(ixic) = ixc
  ! If slave is a super_lord then create the appropriate super_slave names

  allocate (names(slave%n_slave), indx(slave%n_slave))
  do j = 1, slave%n_slave
    if (slave%orientation == 1) then
      slave2 => pointer_to_slave(slave, j)
    else
      slave2 => pointer_to_slave(slave, slave%n_slave+1-j)
    endif

    slave2_name = ''
    lmax = len(slave2%name) - 2
    do k = 1, slave2%n_lord
      lord2 => pointer_to_lord(slave2, k)
      if (lord2%n_lord > 0) lord2 => pointer_to_lord(lord2, 1)
      slave2_name = trim(slave2_name) // trim(lord2%name) // '\'     ! '
      if (len_trim(slave2_name) > lmax) exit
    enddo
    if (len_trim(slave2_name) > lmax) slave2_name = slave2_name(1:lmax) // '\ '
    call set_ele_name (slave2, trim(slave2_name) // int_str(i))
    

    names(j) = slave2_name
  enddo

  ! If a name is not unique then add "#NNN" suffix
  call indexer(names, indx)
  j = 0
  outer_loop: do 
    j = j + 1
    if (j >= slave%n_slave) exit
    if (names(indx(j)) /= names(indx(j+1))) cycle  ! Is unique name
    do i1 = 0, slave%n_slave
      if (j + i1 > slave%n_slave) exit outer_loop
      if (names(indx(j)) /= names(indx(j+i1))) then
        j = j + i1
        exit
      endif

      if (slave%orientation == 1) then
        slave2 => pointer_to_slave(slave, indx(j+i1))
      else
        slave2 => pointer_to_slave(slave, slave%n_slave+1-indx(j+i1))
      endif
      call set_ele_name (slave2, trim(slave2%name) // '#' // int_str(i1+1))
    enddo
  enddo outer_loop

  deallocate(names, indx)

enddo

!

call control_bookkeeper (lat, lord)

end subroutine add_this_multipass

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
! If this is a drift multipass whose multipass_slave elements are the result
! of splitting a drift with superposition then make sure that all split drift elements 
! of the lattice with the same base name have a name of the form "<base_name>#<n>" where
! <n> is an index from 1 for the first split drift.
!-

subroutine drift_multipass_name_correction (lat)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, slave, lord, lord2
type (branch_struct), pointer :: branch

integer i, j, k, ie, ixb, ixb2, ix_pass, n_links, ix_n

character(40) base_name

!

do ie = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(ie)
  ixb = index(lord%name, '#') - 1
  if (lord%key /= drift$ .or. ixb < 1) cycle

  ix_n = 0
  base_name = lord%name(1:ixb) 

  do i = 0, ubound(lat%branch, 1)
    branch => lat%branch(i)
    do j = 1, branch%n_ele_track
      ele => branch%ele(j)
      if (ele%key /= drift$) cycle
      ixb2 = index(ele%name, '#') - 1
      if (ixb2 /= ixb) cycle
      if (base_name(1:ixb) /= ele%name(1:ixb)) cycle
      ! super_slave drifts are temporary constructs that need to be ignored.
      ! This routine will be called later to correct the name of such elements.
      if (ele%slave_status == super_slave$) cycle 
      if (ele%slave_status == multipass_slave$) then
        call multipass_chain (ele, ix_pass, n_links)
        if (ix_pass /= 1) cycle  ! Only do renaming once
        lord2 => pointer_to_lord(ele, 1)
        ix_n = ix_n + 1
        call set_ele_name (lord2, base_name(1:ixb) // '#' // int_str(ix_n))
        do k = 1, lord2%n_slave
          slave => pointer_to_slave(lord2, k)
          call set_ele_name (slave, base_name(1:ixb) // '#' // int_str(ix_n) // '\' // int_str(k))      !'
        enddo
      else
        ix_n = ix_n + 1
        call set_ele_name (ele, base_name(1:ixb) // '#' // int_str(ix_n))
      endif
    enddo
  enddo

enddo

end subroutine drift_multipass_name_correction

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine reallocate_bp_com_const()

implicit none

type (bp_const_struct), allocatable :: var_temp(:)
integer n

!

call move_alloc(bp_com2%const, var_temp)
n = bp_com%i_const_tot+100
allocate (bp_com2%const(n))
n = size(var_temp)
bp_com2%const(1:n) = var_temp

end subroutine reallocate_bp_com_const

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser (but not bmad_parser2)
! This subroutine is not intended for general use.
!-

subroutine parser_add_superimpose (branch, super_ele_in, pele, in_lat, plat)

implicit none

type (ele_struct) super_ele_in
type (ele_struct), save :: super_ele_saved, super_ele
type (ele_struct), pointer :: ref_ele, ele, slave, lord, super_ele_out, ele_at_s
type (ele_pointer_struct), allocatable :: eles(:)
type (parser_lat_struct), target :: plat
type (parser_ele_struct) :: pele
type (parser_ele_struct), pointer :: ref_pele
type (multipass_all_info_struct) m_info
type (lat_struct) :: in_lat
type (branch_struct), target :: branch
type (branch_struct), pointer :: ref_branch
type (lat_struct), pointer :: lat

integer ix, i, j, k, it, nn, i_ele, ib, il
integer n_con, ix_branch, n_loc, n_loc0, ix_insert, ix_pass

character(40) name, ref_name
character(80) line

logical have_inserted, found, err_flag, err

! init

if (.not. bp_com%do_superimpose) return

call init_ele(super_ele_saved)
call init_ele(super_ele)

super_ele = super_ele_in
super_ele%logic = .false.
call settable_dep_var_bookkeeping (super_ele)

super_ele_saved = super_ele     ! in case super_ele_in changes
lat => branch%lat
pele%ix_super_ref_multipass = 0

! Find all matches
! If no refrence element then this implies superposition on branch 0.

if (pele%ref_name == blank_name$) pele%ref_name = '0>>BEGINNING'

call lat_ele_locator (pele%ref_name, lat, eles, n_loc, err)
if (err) then
  call parser_error ('MALFORMED SUPERIMPOSE REFERENCE ELEMENT NAME: ' // pele%ref_name, &
                     'FOR SUPERPOSITION OF: ' // super_ele_saved%name, pele = pele)
  return
endif

if (n_loc > 0) then
  if (eles(1)%ele%key == drift$ .and. index(pele%ref_name, '##') /= 0) then
    call parser_error ('SUPERPOSITION WITH THE REFERENCE ELEMENT BEING THE N^TH DRIFT OF A CERTAIN NAME: ' // trim(pele%ref_name), &
                       'NOT ALLOWED SINCE THIS CAN LEAD TO UNEXPECTED BEHAVIOR.')
    return
  endif
endif

! Throw out super_slave elements.
! Throw out elements that are the same physical element.

n_loc0 = n_loc
i = 1
do
  if (i > n_loc) exit

  if (eles(i)%ele%slave_status == super_slave$) then
    eles(i:n_loc-1) = eles(i+1:n_loc)  ! Remove
    n_loc = n_loc - 1
    cycle
  endif

  if (pele%ix_super_ref_multipass /= 0 .and. eles(i)%ele%iyy /= pele%ix_super_ref_multipass) then
    eles(i:n_loc-1) = eles(i+1:n_loc)  ! Remove
    n_loc = n_loc - 1
    cycle
  endif

  i = i + 1
enddo

if (n_loc == 0 .and. n_loc0 > 0) then
  call parser_error ('SUPERPOSITION WITH THE REFERENCE ELEMENT BEING A SUPER SLAVE NOT ALLOWED: ' // trim(pele%ref_name))
  return
endif


! Group and overlay elements may have not yet been transfered from in_lat to lat.
! So search in_lat for a match if there has not been a match using lat. 
! Note: Using a group or overlay as the reference element is only valid if there is only one slave element.

if (n_loc == 0) then
  call lat_ele_locator (pele%ref_name, in_lat, eles, n_loc, err) ! Notice: Searching in_lat
  if (err) then
    call parser_error ('MALFORMED SUPERIMPOSE REFERENCE NAME: ' // pele%ref_name)
    return
  endif
  ! Since superposition in bmad_parser is done branch-by-branch, n_loc = 0 might not be an error.
  ! This will get checked later.
  if (n_loc == 0) return
  if (n_loc > 1) then
    call parser_error ('CONFUSED SUPERPOSITION. PLEASE CONTACT A BMAD MAINTAINER.')
    return
  endif

  ! n_loc = 1 case
  if (eles(1)%ele%key /= group$ .and. eles(1)%ele%key /= overlay$) return ! Should be handled in another branch.

  ref_ele => eles(1)%ele
  ref_pele => plat%ele(ref_ele%ixx)
  ref_name = ref_pele%control(1)%name

  do i = 2, size(ref_pele%control)
    if (ref_pele%control(i)%name /= ref_name) then
      call parser_error ('SUPERPOSITION USING A GROUP OR OVERLAY AS A REFERENCE ELEMENT IS ONLY ALLOWED', &
                         'WHEN THE GROUP OR OVERLAY CONTROLS A SINGLE ELEMENT.', &
                         'FOR SUPERPOSITION OF: ' // super_ele_saved%name, pele = pele)
      return
    endif
  enddo 

  call lat_ele_locator (ref_name, lat, eles, n_loc, err)
  if (err) then
    call parser_error ('MALFORMED OVERLAY OR GROUP SLAVE NAME: ' // ref_name)
    return
  endif
endif

! If a ref element is outside, and the superposition offset puts the super_ele into a multipass
! region, must add further superpositions to keep the multipass regions on separate passes looking the same.
! shift the ref element to the multipass region.

if (n_loc == 1) then
  ref_ele => eles(1)%ele
  ref_branch => pointer_to_branch(ref_ele)
  if (ref_ele%iyy == 0 .and. ref_branch%ix_branch == branch%ix_branch) then
    call compute_super_lord_s (eles(1)%ele, super_ele, pele, ix_insert)
    ele_at_s => pointer_to_element_at_s (branch, super_ele%s_start, .true., err_flag)
    if (err_flag) then
      call parser_error ('PROBLEM SUPERIMPOSING: ' // super_ele%name)
    else
      if (ele_at_s%slave_status == super_slave$) ele_at_s => pointer_to_lord(ele_at_s, 1)
      if (ele_at_s%iyy /= 0) then  ! If in multipass region...
        pele%ref_name = ele_at_s%name
        pele%ref_pt = anchor_end$
        pele%ele_pt = anchor_end$
        pele%offset = super_ele%s - ele_at_s%s
        pele%ix_super_ref_multipass = ele_at_s%iyy
        ! multipass bookkeeping has not been done yet and ele_at_s%name is
        ! the name of the multipass_lord (something like "Q1". 
        call lat_ele_locator (ele_at_s%name, lat, eles, n_loc, err)

        i = 1 ! throw out elements that are the same physical element
        do
          if (i > n_loc) exit
          if (eles(i)%ele%iyy /= pele%ix_super_ref_multipass .or. eles(i)%ele%slave_status == super_slave$) then
            eles(i:n_loc-1) = eles(i+1:n_loc)  ! Remove
            n_loc = n_loc - 1
          else
            i = i + 1
          endif
        enddo
      endif
    endif
  endif
endif

! Tag reference elements using %logic flag which is not otherwise used during parsing.

branch%ele(:)%logic = .false.
lat%ele(:)%logic = .false.  ! For lord elements
 
do i = 1, n_loc
  eles(i)%ele%logic = .true. ! Tag reference element.
enddo

! Note: branch%n_ele_max will vary when an element is superimposed.

do 
  have_inserted = .false.
  i_ele = -1

  do 
    i_ele = i_ele + 1
    if (i_ele > branch%n_ele_track) exit
    call do_this_superimpose(lat, branch%ele(i_ele), have_inserted)
  enddo

  i_ele = lat%n_ele_track
  do
    i_ele = i_ele + 1
    if (i_ele > lat%n_ele_max) exit
    call do_this_superimpose(lat, lat%ele(i_ele), have_inserted)
  enddo

  if (.not. have_inserted) exit
enddo

!-------------------------------------------------------------
contains

subroutine do_this_superimpose (lat, start_ele, have_inserted)

type (lat_struct), target :: lat
type (ele_struct), target :: start_ele
type (ele_struct), pointer :: ref_ele
logical have_inserted

!

ref_ele => start_ele
if (.not. ref_ele%logic) return
ref_ele%logic = .false.  ! So only use this reference once

if (ref_ele%slave_status == super_slave$) then
  call parser_error ('SUPERIMPOSE REFERENCE ELEMENT MUST NOT BE A SUPER_SLAVE: ' // ref_ele%name)
  return
endif

if (ref_ele%key == group$ .or. ref_ele%key == overlay$) return
if (ref_ele%key == girder$) return
if (ref_ele%slave_status == super_slave$) return

ref_branch => pointer_to_branch(ref_ele)
if (ref_branch%ix_branch /= branch%ix_branch) return

call compute_super_lord_s (ref_ele, super_ele, pele, ix_insert)
super_ele%iyy = ref_ele%iyy   ! Multipass info
call check_for_superimpose_problem (branch, super_ele, err_flag, ref_ele, pele%wrap_superimpose); if (err_flag) return
call string_trim(super_ele_saved%name, super_ele_saved%name, ix)
super_ele%name = super_ele_saved%name(:ix)            
call add_superimpose (lat, super_ele, branch%ix_branch, err_flag, super_ele_out, save_null_drift = .true., &
      create_jumbo_slave = pele%create_jumbo_slave, ix_insert = ix_insert, mangle_slave_names = .false., wrap = pele%wrap_superimpose)
if (err_flag) then
  bp_com%error_flag = .true.
  return
endif
call control_bookkeeper (lat, super_ele_out)

call s_calc (lat)
have_inserted = .true.   

end subroutine do_this_superimpose

end subroutine parser_add_superimpose

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser2 (but not bmad_parser).
! This subroutine is not intended for general use.
!-

subroutine parser2_add_superimpose (lat, super_ele_in, pele, in_lat)

implicit none

type (ele_struct) super_ele_in
type (ele_struct), save :: super_ele_saved, super_ele
type (ele_struct), pointer :: ref_ele, ele
type (ele_pointer_struct), allocatable :: eles(:)
type (parser_ele_struct) :: pele
type (parser_ele_struct), pointer :: ref_pele
type (lat_struct), optional :: in_lat
type (branch_struct), pointer :: branch
type (branch_struct), pointer :: ref_branch
type (lat_struct), target :: lat

integer ix, i, j, k, it, nn, i_ele, n_super
integer n_con, ix_branch, n_loc, ix_insert, ix_pass

character(40) name, ref_name
character(80) line

logical have_inserted, found, err_flag, err

! init

if (.not. bp_com%do_superimpose) return

call settable_dep_var_bookkeeping (super_ele_in)

call init_ele(super_ele_saved)
call init_ele(super_ele)

super_ele = super_ele_in
super_ele%logic = .false.
super_ele_saved = super_ele     ! in case super_ele_in changes

! Find all matches
! If no refrence element then this implies superposition on branch 0.

if (pele%ref_name == blank_name$) pele%ref_name = '0>>BEGINNING'

call lat_ele_locator (pele%ref_name, lat, eles, n_loc, err)
if (err) then
  call parser_error ('MALFORMED SUPERIMPOSE REFERENCE ELEMENT NAME: ' // pele%ref_name, &
                     'FOR SUPERPOSITION OF: ' // super_ele_saved%name, pele = pele)
  return
endif

! Throw out super_slave elements.

i = 1
do
  if (i > n_loc) exit

  if (eles(i)%ele%slave_status == super_slave$) then
    eles(i:n_loc-1) = eles(i+1:n_loc)  ! Remove
    n_loc = n_loc - 1
    cycle
  endif

  i = i + 1
enddo

! The superimpose may reference elements that were defined but not used in the construction of the lattice
! by bmad_parser. If so, this is not an error and just do not do the superposition.

if (n_loc == 0 .and. present(in_lat)) then
  call lat_ele_locator (pele%ref_name, in_lat, eles, n_loc, err) ! Notice: Searching in_lat
  if (err) then
    call parser_error ('MALFORMED SUPERIMPOSE REFERENCE NAME: ' // pele%ref_name)
    return
  endif
  if (n_loc > 0) return
endif

! Check

if (n_loc == 0) then
  call parser_error ('NO MATCH FOR REFERENCE ELEMENT: ' //  pele%ref_name, &      
                     'FOR SUPERPOSITION OF: ' // super_ele%name, pele = pele)
  return
endif

! If a ref element is outside is outside, and the superposition offset puts the super_ele into a multipass
! region, must add further superpositions to keep the multipass regions on separate passes looking the same.
! shift the ref element to the multipass region.

! Tag reference elements using %logic flag which is not otherwise used during parsing.

do i = 0, ubound(lat%branch, 1)
  lat%branch(i)%ele(:)%logic = .false.
  lat%branch(i)%ele(:)%iyy = 0
enddo

do i = 1, n_loc
  eles(i)%ele%logic = .true. ! Tag reference element.
enddo

n_super = 0
do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  i_ele = -1
  do 
    i_ele = i_ele + 1
    if (i_ele > branch%n_ele_max) exit
    n_super = n_super + 1
    call do2_this_superimpose(lat, branch%ele(i_ele), have_inserted, pele, n_super)
  enddo
enddo

!-------------------------------------------------------------
contains

subroutine do2_this_superimpose (lat, ref0_ele, have_inserted, pele, n_super)

type (lat_struct), target :: lat
type (ele_struct), target :: ref0_ele
type (ele_struct), pointer :: ref_ele, ele, slave, slave2, ele_at_s, super_ele_out
type (branch_struct), pointer :: branch
type (parser_ele_struct), target :: pele
type (parser_ele_struct) :: pele2
type (multipass_all_info_struct) m_info
type (lat_ele_loc_struct), allocatable :: m_slaves(:)

logical have_inserted
integer i, j, k, ib, n_slave, n_super

!

ref_ele => ref0_ele
if (.not. ref_ele%logic) return
ref_ele%logic = .false.  ! So only use this reference once

if (ref_ele%slave_status == super_slave$) then
  call parser_error ('SUPERIMPOSE REFERENCE ELEMENT MUST NOT BE A SUPER_SLAVE: ' // ref_ele%name)
  return
endif

! Note: Using a group or overlay as the reference element is only valid if there is only one slave element.

if (ref_ele%key == group$ .or. ref_ele%key == overlay$) then
  slave => pointer_to_slave(ref_ele, 1)

  do i = 2, ref_ele%n_slave
    slave2 => pointer_to_slave(ref_ele, i)
    if (slave2%ix_ele /= slave%ix_ele .or. slave2%ix_branch /= slave%ix_branch) then
      call parser_error ('SUPERPOSITION USING A GROUP OR OVERLAY AS A REFERENCE ELEMENT IS ONLY ALLOWED', &
                         'WHEN THE GROUP OR OVERLAY CONTROLS A SINGLE ELEMENT.', &
                         'FOR SUPERPOSITION OF: ' // super_ele_saved%name, pele = pele)
      return
    endif
  enddo

  ref_ele => slave
endif

!

pele2 = pele

! If a ref element is outside is outside, and the superposition offset puts the super_ele into a multipass
! region, must use the appropriate multipass_lord as the reference element so that a superposition is done
! one for each pass.

if (.not. associated (pointer_to_multipass_lord(ref_ele))) then
  branch => pointer_to_branch(ref_ele)
  call compute_super_lord_s (ref_ele, super_ele, pele2, ix_insert)
  ele_at_s => pointer_to_element_at_s (branch, super_ele%s_start, .true., err_flag)
  if (ele_at_s%slave_status == super_slave$) ele_at_s => pointer_to_lord(ele_at_s, 1)
  if (ele_at_s%slave_status == multipass_slave$) then
    ref_ele => ele_at_s
    pele2%ref_pt = anchor_end$
    pele2%ele_pt = anchor_end$
    pele2%offset = super_ele%s - ele_at_s%s
  endif
endif

if (ref_ele%slave_status == multipass_slave$) then
  ref_ele => pointer_to_lord(ref_ele, 1)
endif

if (ref_ele%key == girder$) return
if (ref_ele%slave_status == super_slave$) return

! If superimposing on a multipass_lord
! then the superposition must be done at all multipass locations.

if (ref_ele%lord_status == multipass_lord$) then
  ref_ele%iyy = -1  ! So we can find it again after it has moved during superposition

  n_slave = ref_ele%n_slave
  allocate (m_slaves(n_slave))
  do j = 1, n_slave
    ! Find where ref_ele is
    do k = lat%n_ele_track+1, lat%n_ele_max
      ref_ele => lat%ele(k)
      if (ref_ele%iyy == -1) exit
    enddo
    if (j == n_slave) ref_ele%iyy = 0  ! No longer needed. Important to reset to zero.
    slave => pointer_to_slave(ref_ele, j)
    branch => pointer_to_branch(slave)
    call compute_super_lord_s (slave, super_ele, pele2, ix_insert)
    call add_superimpose (lat, super_ele, branch%ix_branch, err_flag, super_ele_out, save_null_drift = .false., &
                 create_jumbo_slave = pele2%create_jumbo_slave, ix_insert = ix_insert, mangle_slave_names = .false., wrap = pele2%wrap_superimpose) 
    if (err_flag) bp_com%error_flag = .true.
    super_ele_out%iyy = n_super  ! Is unique
  enddo

  ! If the super_lords have a single super_slave and the super_slave
  ! has only a single super_lord, the super_lords
  ! can be eliminated and the created multipass_lord can control the
  ! super_slaves directly.

  j = 0

  do ib = 0, ubound(lat%branch, 1)
    branch => lat%branch(ib)
    do i = 1, branch%n_ele_max
      if (branch%ele(i)%iyy == n_super) then
        j = j + 1
        m_slaves(j) = ele_loc(branch%ele(i))
      endif
    enddo
  enddo

  ele => pointer_to_ele (lat, m_slaves(1))
  if (ele%lord_status == super_lord$ .and. ele%n_slave == 1) then
    slave => pointer_to_slave(ele, 1)
    if (slave%n_lord == 1) then
      do i = 1, n_slave
        ele => pointer_to_ele (lat, m_slaves(i))
        ele%ix_ele = -1 ! Mark for deletion
        ele => pointer_to_slave(ele, 1)
        ele%name = super_ele_saved%name
        m_slaves(i) = ele_loc(ele)
      enddo
    endif
  endif

  ! Remove eles marked for deletion. But first shift m_slaves list

  do i = 1, size(m_slaves)
    ele => pointer_to_ele (lat, m_slaves(i))
    m_slaves(i)%ix_ele = m_slaves(i)%ix_ele - count(lat%branch(ele%ix_branch)%ele(1:ele%ix_ele)%key == -1)
  enddo

  call remove_eles_from_lat (lat, .false.)

  ! Add a multipass_lord to control the created super_lords.

  call add_this_multipass (lat, m_slaves, super_ele_saved) 

  deallocate (m_slaves)

!-----------------------
! Else not superimposing on a multipass_lord ...
! [Note: Only will superimpose on a multipass_lord in bmad_parser2.]

else
  ref_branch => pointer_to_branch(ref_ele)
  if (ref_branch%ix_branch /= branch%ix_branch) return

  call compute_super_lord_s (ref_ele, super_ele, pele2, ix_insert)
  super_ele%iyy = ref_ele%iyy   ! Multipass info
  call check_for_superimpose_problem (branch, super_ele, err_flag, ref_ele, pele%wrap_superimpose); if (err_flag) return
  call string_trim(super_ele_saved%name, super_ele_saved%name, ix)
  super_ele%name = super_ele_saved%name(:ix)            
  call add_superimpose (lat, super_ele, branch%ix_branch, err_flag, super_ele_out, save_null_drift = .true., &
              create_jumbo_slave = pele2%create_jumbo_slave, ix_insert = ix_insert, mangle_slave_names = .false., wrap = pele2%wrap_superimpose)
  if (err_flag) bp_com%error_flag = .true.
  call control_bookkeeper (lat, super_ele_out)
endif

!---------------------

call s_calc (lat)
have_inserted = .true.   

end subroutine do2_this_superimpose

end subroutine parser2_add_superimpose

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine compute_super_lord_s (ref_ele, super_ele, pele, ix_insert)

implicit none

type (ele_struct), target :: ref_ele, super_ele
type (ele_struct), pointer :: slave
type (parser_ele_struct) pele
type (branch_struct), pointer :: branch

integer i, ix, ix_insert, ele_pt, ref_pt

real(rp) s_ref_begin, s_ref_end, s0
logical reflected_or_reversed

! Find the reference point on the element being superimposed.
! ref_ele%select has been set such that if True then reference element was in
! a reflected or reversed line.

ix_insert = -1
reflected_or_reversed = ref_ele%select

super_ele%orientation = ref_ele%orientation
if (reflected_or_reversed) then
  super_ele%s = -pele%offset
else
  super_ele%s = pele%offset
endif

ele_pt = pele%ele_pt
if (reflected_or_reversed) then
  if (ele_pt == anchor_beginning$) then
    ele_pt = anchor_end$
  elseif (ele_pt == anchor_end$) then
    ele_pt = anchor_beginning$
  endif
endif

if (ele_pt == anchor_beginning$) then
  super_ele%s = super_ele%s + super_ele%value(l$)
elseif (ele_pt == anchor_center$ .or. ele_pt == not_set$) then
  super_ele%s = super_ele%s + 0.5_rp * super_ele%value(l$)
elseif (ele_pt /= anchor_end$) then
  call parser_error ('ERROR IN COMPUTE_SUPER_LORD_S: CONTROL #1 INTERNAL ERROR!')
  if (global_com%exit_on_error) call err_exit
endif

! Find the refernce point in the lattice.

select case (ref_ele%key)
case (overlay$, group$, girder$)
  s_ref_begin = 1d10
  s_ref_end = 0
  do i = 1, ref_ele%n_slave
    slave => pointer_to_slave(ref_ele, i)
    s_ref_begin = min(s_ref_begin, slave%s_start)
    s_ref_end = max(s_ref_end, slave%s)
  enddo
case (ramper$)
  call parser_error ('SUPERPOSING: ' // super_ele%name, 'UPON RAMPER' // pele%ref_name)
  return
case default
  s_ref_begin = ref_ele%s_start
  s_ref_end = ref_ele%s
end select

! Now compute the s position at the end of the element and put it in ele%s.

ref_pt = pele%ref_pt
if (reflected_or_reversed) then
  if (ref_pt == anchor_beginning$) then
    ref_pt = anchor_end$
  elseif (ref_pt == anchor_end$) then
    ref_pt = anchor_beginning$
  endif
endif

if (ref_pt == anchor_beginning$) then
  super_ele%s = super_ele%s + s_ref_begin
elseif (ref_pt == anchor_center$ .or. ref_pt == not_set$) then
  super_ele%s = super_ele%s + 0.5_rp * (s_ref_begin + s_ref_end)
elseif (ref_pt == anchor_end$) then
  super_ele%s = super_ele%s + s_ref_end
else
  call parser_error ('ERROR IN COMPUTE_SUPER_LORD_S: CONTROL #2 INTERNAL ERROR!')
  if (global_com%exit_on_error) call err_exit
endif

! A superimpose can wrap around the beginning or the end of the lattice. 
! This is done independent of the geometry. The reason why this is geometry 
! independent is that it is sometimes convenient to treat a closed lattice as open.

branch => pointer_to_branch(ref_ele)
s0 = branch%ele(0)%s

if (pele%wrap_superimpose) then
  if (super_ele%s > branch%ele(branch%n_ele_track)%s) then
    super_ele%s = super_ele%s - branch%param%total_length
  elseif (super_ele%s < s0) then
    super_ele%s = super_ele%s + branch%param%total_length
  endif

  super_ele%s_start = super_ele%s - super_ele%value(l$)
  if (super_ele%s_start < s0) then
    super_ele%s_start = super_ele%s_start + branch%param%total_length
  endif
endif

! The "nominal" insert point is at the downstream end of element with index ix_insert.
! ix_insert is used for positioning zero length super_lords in case
! They are next to an element which is also zero length.

! First special case: A null_ele that has been moved to the lord list

if (ref_ele%key == null_ele$ .and. ref_ele%ix_branch == 0 .and. ref_ele%ix_ele > branch%lat%n_ele_track) then
  ix_insert = -1  ! Will be ignored 

elseif (ref_ele%value(l$) == 0 .and. super_ele%value(l$) == 0 .and. pele%offset == 0) then
  if (ref_ele%ix_ele == 0) then  ! Automatically must be at downstream end.
    ix_insert = 0
  elseif (pele%ref_pt == anchor_beginning$) then
    ix_insert = ref_ele%ix_ele - 1
  elseif (pele%ref_pt == anchor_end$) then
    ix_insert = ref_ele%ix_ele
  endif

! For elements with finite length just return the ref element index.
! This will be useful in superpositions near elements that have negative length (EG a patch).
else
  slave => ref_ele
  if (slave%n_slave /= 0) slave => pointer_to_slave(slave, 1) ! Does not matter which slave chosen
  if (slave%n_slave /= 0) slave => pointer_to_slave(slave, 1) ! Does not matter which slave chosen
  ix_insert = slave%ix_ele
endif

end subroutine compute_super_lord_s

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine check_for_superimpose_problem (branch, super_ele, err_flag, ref_ele, wrap)
!
! Subroutine to check if there is a problem superimposing an element when there is multipass.
! In particular will check that:
!   1) If the ref_ele is part of a multipass region then super_ele must be superimposed
!      within the region.
! Or:
!   2) If the ref_ele is not part of a multipass region then super_ele must also not
!      be part of a multipass region.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine check_for_superimpose_problem (branch, super_ele, err_flag, ref_ele, wrap)

implicit none

type (ele_struct) super_ele
type (ele_struct), optional :: ref_ele
type (ele_struct), pointer :: ele, ele1, ele2, ele_stop
type (branch_struct) :: branch
real(rp) eps
logical err_flag, wrap
integer ix1, ix2


! Check for out-of-bounds.
! If wrap = False then out-of-bounds is not an error.

err_flag = .true.
eps = bmad_com%significant_length

ele1 => pointer_to_element_at_s (branch, super_ele%s_start + eps, .true., err_flag)
if (err_flag) then
  if (wrap) then
    call parser_error ('BAD SUPERIMPOSE OF: ' // super_ele%name, 'UPSTREAM ELEMENT EDGE OUT OF BOUNDS.')
    return
  else
    ele1 => branch%ele(0)
  endif
endif

ele2 => pointer_to_element_at_s (branch, super_ele%s - eps, .false., err_flag, print_err = wrap)
if (err_flag) then
  if (wrap) then
    call parser_error ('BAD SUPERIMPOSE OF: ' // super_ele%name, 'DOWNSTREAM ELEMENT EDGE OUT OF BOUNDS.')
    return
  else
    ele2 => branch%ele(branch%n_ele_track)
  endif
endif

! Ref ele check.

if (ele1%slave_status == super_slave$) ele1 => pointer_to_lord(ele1, 1)
if (ele2%slave_status == super_slave$) ele2 => pointer_to_lord(ele2, 1)

if (present(ref_ele)) then
  if (ref_ele%iyy /= 0) then     ! Ref element in multipass region
    if (abs(super_ele%value(l$)) < eps .and. (ele1%iyy /= 0 .or. ele2%iyy /= 0)) return ! At multipass edge is OK
    if (ele1%iyy == 0 .or. ele2%iyy == 0) then
      call parser_error ('SUPERIMPOSE OF: ' // super_ele%name, &
                         'WITH REFERENCE ELEMENT: ' // trim(ref_ele%name) // ' IN BRANCH \i0\ ', &
                         'USES MULTIPASS REFERENCE ELEMENT BUT OFFSET PLACES IT OUT OF THE MULTIPASS REGION!', &
                          i_array = [branch%ix_branch])
      return
    endif

  else
    if (abs(super_ele%value(l$)) < eps .and. (ele1%iyy == 0 .or. ele2%iyy == 0)) return  ! At multipass edge is OK
    if (ele1%iyy /= 0 .or. ele2%iyy /= 0) then
      call parser_error ('SUPERIMPOSE OF: ' // super_ele%name, &
                         'WITH REFERENCE ELEMENT: ' // trim(ref_ele%name) // ' IN BRANCH \i0\ ', &
                         'USES NON-MULTIPASS REFERENCE ELEMENT BUT OFFSET PLACES IT IN A MULTIPASS REGION!', &
                          i_array = [branch%ix_branch])
      return
    endif
  endif
endif

err_flag = .false.

end subroutine check_for_superimpose_problem 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine get_sequence_args (seq_name, arg_list, delim, err_flag)
!
! Subroutine to get the argument list for a replacement_line or a list.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine get_sequence_args (seq_name, arg_list, delim, err_flag)

implicit none

integer n_arg, ix_word

character(*), allocatable :: arg_list(:)
character(*) seq_name
character(1) delim
character(40) name(20), word

logical delim_found, err_flag

!

n_arg = 0
err_flag = .true.

do
  call get_next_word (word, ix_word, '(): =,', delim, delim_found, .true.)
  if (ix_word == 0 .or. delim == '( :=') then
    call parser_error ('BAD ARGUMENT LIST FOR: ', seq_name)
    return
  endif
  n_arg = n_arg + 1
  name(n_arg) = word
  if (delim == ')') exit
enddo

err_flag = .false.
if (allocated(arg_list)) deallocate(arg_list)
allocate (arg_list(n_arg))
arg_list = name(1:n_arg)

end subroutine get_sequence_args

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parse_line_or_list (sequence, iseq_tot, lat, top_level)
!
! Subroutine to parse a sequence.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

recursive subroutine parse_line_or_list (sequence, iseq_tot, lat, top_level)

implicit none

type (seq_struct), allocatable, target :: sequence(:)
type (seq_struct), pointer :: seq, sub_seq
type (seq_ele_struct), pointer :: seq_ele_arr(:)
type (seq_ele_struct), pointer :: seq_ele_arr2(:)
type (seq_ele_struct), pointer :: s_ele
type (lat_struct) lat

integer ix_ele, iseq_tot, ix_word, ix_word2, ix, i, n, ios, rcount
integer, save :: ix_internal = 0

character(1) delim, c_delim
character(40) str, word, word2
character(n_parse_line) parse_line_saved

logical delim_found, replacement_line_here, c_delim_found
logical err_flag, top_level

! init

allocate (seq_ele_arr(ubound(lat%ele, 1)))

! save info on what file we are parsing for error messages.

seq => sequence(iseq_tot)
seq%file_name = bp_com%current_file%full_name 
seq%ix_file_line = bp_com%current_file%i_line

! first thing should be a "("

call get_next_word(word, ix_word, ':=(),', delim, delim_found, .true.)

if (delim /= '(') call parser_error ('EXPECTING "(", GOT: ' // delim, 'FOR LINE: ' // seq%name)
if (ix_word /= 0)  call parser_error ('EXTRANEOUS STUFF BEFORE "(", GOT: ' // word, 'FOR LINE: ' // seq%name)

! now parse list proper

ix_ele = 1
s_ele => seq_ele_arr(ix_ele)

do 

  call get_next_word (word, ix_word, ':=(,)[]*@', delim, delim_found, .true.)

  s_ele%rep_count = 1

  if (word(1:1) == '-') then
    s_ele%ele_order_reflect = .true.
    word = word(2:)
    ix_word = ix_word - 1
  endif

  if (word(1:1) == '-') then
    s_ele%ele_orientation = -1
    word = word(2:)
    ix_word = ix_word - 1
  endif

  if (delim == '*') then    ! E.g. '-3*(A,B)'
    ! Evaluate the rep count.
    read (word, *, iostat = ios) rcount
    if (ix_word == 0 .or. ios /= 0) then
      call parser_error ('MALFORMED REPETION COUNT FOUND IN SEQUENCE: ' // seq%name)
      return
    endif
    s_ele%rep_count = rcount
    call get_next_word (word, ix_word, ':=(,)[]*@', delim, delim_found, .true.)
  endif

  s_ele%name = word
  if (word /= ' ') then
    if (.not. verify_valid_name (word, ix_word)) return
  endif

  ! Check for line slice or tag

  if (delim == '@') then   ! "tag@line_name" syntax
    s_ele%tag = s_ele%name
    call get_next_word (s_ele%name, ix_word, '[]:=(,)', delim, delim_found, .true.)
    if (.not. verify_valid_name (s_ele%name, ix_word)) return
  endif

  if (delim == '[') then
    call get_next_word (s_ele%slice_start, ix_word, '[]:=(,)', delim, delim_found, .true.)
    select case (delim)

    case (']')  ! Old style "line_name[tag]" tag syntax
      s_ele%tag = s_ele%slice_start
      s_ele%slice_start = ''
      call parser_error ('OLD STYLE "LINE_NAME[TAG]" LINE TAG SYNTAX DETECTED.', &
                         'PLEASE CONVERT TO NEW STYLE "TAG@LINE_NAME" SYNTAX.', &
                         'WILL RUN AS NORMAL FOR NOW...', level = s_warn$)
    case (':')
      call get_next_word (s_ele%slice_end, ix_word, '[]:=(,)', delim, delim_found, .true.)
      if (delim /= ']') then
        call parser_error ('NO MATCHING "]" FOUND FOR OPENING "[" IN SEQUENCE: ' // seq%name)
        return
      endif
    case default
      call parser_error ('NO MATCHING "]" FOUND FOR OPENING "[" IN SEQUENCE: ' // seq%name)
      return
    end select

    call get_next_word (word, ix_word, '[]:=(,)', delim, delim_found, .true.)
    if (ix_word > 0) then
      call parser_error ('ILLEGAL CHARACTERS AFTER CLOSING "]" FOUND IN SEQUENCE: ' // seq%name)
      return
    endif   
  endif

  ! Check for a subline or replacement line.
  ! If there is one then save as an internal sequence.

  replacement_line_here = .false.

  if (delim == '(') then ! subline or replacement line

    ! if a subline...
    if (s_ele%name == '') then  
      ix_internal = ix_internal + 1
      write (str, '(a, i3.3)') '#Internal', ix_internal   ! unique name 
      s_ele%name = str
      iseq_tot = iseq_tot + 1
      if (iseq_tot > size(sequence)) call reallocate_sequence(sequence, 2*iseq_tot)
      sub_seq => sequence(iseq_tot) 
      sub_seq%name = str
      sub_seq%type = seq%type
      sub_seq%multipass = seq%multipass
      if (sub_seq%type == replacement_line$) then
        ix = size (seq%dummy_arg)
        allocate (sub_seq%dummy_arg(ix), sub_seq%corresponding_actual_arg(ix), s_ele%actual_arg(ix))
        sub_seq%dummy_arg = seq%dummy_arg
        s_ele%actual_arg = seq%dummy_arg
      endif
      bp_com%parse_line = '(' // bp_com%parse_line
      call parse_line_or_list (sequence, iseq_tot, lat, .false.)

    ! else this is a replacement line
    else    
      replacement_line_here = .true.
      call get_sequence_args (s_ele%name, s_ele%actual_arg, delim, err_flag)
      if (err_flag) return
    endif

    call get_next_word(word, ix_word, ':=(),', delim, delim_found, .true.)
    if (word /= ' ') call parser_error &
              ('NO COMMA AFTER SUBLINE OR REPLACEMENT LINE. FOUND: ' // &
               word, 'IN THE SEQUENCE: ' // seq%name)
  endif

  if (s_ele%name == ' ') call parser_error ('SUB-ELEMENT NAME IS BLANK FOR LINE/LIST: ' // seq%name)

  ! if a replacement line then look for element in argument list

  s_ele%ix_arg = 0
  if (seq%type == replacement_line$) then
    do i = 1, size(seq%dummy_arg)
      if (seq%dummy_arg(i) == s_ele%name) then
        s_ele%ix_arg = i
        exit
      endif
    enddo
  endif

  ! 

  n = size(seq_ele_arr)
  ix_ele = ix_ele + 1

  if (ix_ele > n) then
    allocate (seq_ele_arr2(n))      
    seq_ele_arr2 = seq_ele_arr(1:n)
    deallocate (seq_ele_arr)
    allocate (seq_ele_arr(n+1000))
    seq_ele_arr(1:n) = seq_ele_arr2
    deallocate(seq_ele_arr2)
  endif

  s_ele => seq_ele_arr(ix_ele)

  if (delim == ')') exit

  if (delim /= ',') then
    call parser_error ('EXPECTING "," GOT: ' // delim, 'FOR LINE: ' // seq%name)
    exit
  endif
         
enddo

! make sure there is nothing else if at top level

if (top_level) then
  call get_next_word(word, ix_word, ':=() ', delim, delim_found, .true.)
  if (delim_found .or. ix_word /= 0) call parser_error  &
        ('EXTRA CHARACTERS AFTER CLOSING ")"',  'FOR LINE: ' // seq%name)
endif

! transfer

ix_ele = ix_ele - 1
allocate (seq%ele(ix_ele))

do i = 1, ix_ele
  seq%ele(i) = seq_ele_arr(i)
enddo

deallocate (seq_ele_arr)

end subroutine parse_line_or_list

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine allocate_plat (plat, n_ele_max) 
!
! Subroutine to allocate allocatable array sizes.
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

Subroutine allocate_plat (plat, n_ele_max)

implicit none

type (parser_lat_struct) plat, temp_plat

integer i, n_now, n_ele_max

! assume all the arrays have the same size

if (allocated(plat%ele)) then
  n_now = ubound(plat%ele, 1)
  call move_alloc (plat%ele, temp_plat%ele)
  allocate (plat%ele(0:n_ele_max))
  plat%ele(0:n_now) = temp_plat%ele
  deallocate (temp_plat%ele)

else
  allocate (plat%ele(0:n_ele_max))
  n_now = -1
endif

! %ixx is used as a pointer from the in_lat%ele array to the plat%ele array

do i = n_now+1, ubound(plat%ele, 1)
  plat%ele(i)%ele_name = ''
  plat%ele(i)%ref_name = blank_name$
  plat%ele(i)%ref_pt  = not_set$
  plat%ele(i)%ele_pt  = not_set$
  plat%ele(i)%offset  = 0
  plat%ele(i)%create_jumbo_slave = .false.
enddo

end subroutine allocate_plat

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_add_lord (lord_lat, n_ele_max, plat, lat, check_lat)
!
! Subroutine to add overlay, group, and girder lords to the lattice.
! For overlays and groups: If multiple elements have the same name then 
! use all of them.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
! Input:
!   lord_lat    -- lat_struct: List of lord elements to add to lat.
!   n_ele_max   -- integer: lord elements in lord_lat are in range [1:n_ele_max].
!   plat        -- parser_lat_struct: Extra info needed to place the lord elements
!
! Output:
!   lat         -- lat_struct: Lattice to add lord elements to.
!   check_lat   -- lat_struct, optional: If slave elements of a lord are not in lat but
!                   are in check_lat, do not issue error message about slave elements
!                   not found. 
!-

subroutine parser_add_lord (lord_lat, n_ele_max, plat, lat, check_lat)

implicit none

type multi_ele_pointer_struct
  type (ele_pointer_struct), allocatable :: eles(:)
  integer n_loc
end type  

type (lat_struct), target :: lord_lat, lat
type (lat_struct), optional :: check_lat
type (ele_struct), pointer :: lord, slave, ele, g_lord, g_slave0, g_slave1
type (parser_lat_struct), target :: plat
type (parser_ele_struct), pointer :: pele
type (control_struct), allocatable, target :: cs(:)
type (branch_struct), pointer :: branch
type (multi_ele_pointer_struct), allocatable :: m_eles(:)
type (ele_pointer_struct), allocatable :: in_eles(:)
type (parser_controller_struct), pointer :: pc
type (control_struct), pointer :: con

integer ip, n_in, ic, ig, k, ix, ib, ie_start, n_list, ns, ixs, ii, ix_end, n_ele_max
integer n_slave, nn, n_loc, n_names, n_in_loc, ix1_slave
integer ix_lord, k_slave, ix_ele_now, ix_super_lord_end

character(60) err_str
character(40) slave_name, attrib_name

logical err, created_girder_lord, err_flag, matched_to_drift, have_ignored_a_drift
logical have_wrapped

! loop over lord elements

main_loop: do n_in = 1, n_ele_max

  lord => lord_lat%ele(n_in)  ! next lord to add
  pele => plat%ele(lord%ixx)
  
  !-----------------------------------------------------
  ! overlays, groups, and rampers

  select case (lord%key)

  case (ramper$)
    if (allocated(cs)) deallocate(cs)
    allocate (cs(size(pele%control)))

    nn = size(pele%control)
    do ip = 1, nn
      pc => pele%control(ip)

      if (allocated(pc%y_knot)) then
        cs(ip)%y_knot = pc%y_knot
      else
        call reallocate_expression_stack (cs(ip)%stack, pc%n_stk)
        cs(ip)%stack = pc%stack(1:pc%n_stk)
      endif

      ! Note: It is not possible to check attrib_name for a misspelling since the ramper may control
      ! an overlay or group variable that has a non-standard variable name.
      attrib_name = pc%attrib_name
      if (attrib_name == blank_name$) attrib_name = pele%default_attrib
      cs(ip)%ix_attrib = attribute_index(0, attrib_name)
      cs(ip)%attribute = attrib_name
      cs(ip)%slave_name = pc%name
    enddo

    ! Create the ramper

    call new_control (lat, ix_lord, lord%name)  ! get index in lat where lord goes
    lat%ele(ix_lord) = lord

    call create_ramper (lat%ele(ix_lord), cs(1:nn), err)

  case (overlay$, group$)
 
    ! If a slave name does not match any name in lord_lat then this is an error (to catch typos).
    ! If a slave is defined in lord_lat but not present in lat then the slave is ignored.
    ! If all slave elements are defined in lord_lat, but are not present in lat, then
    ! this lord can be ignored.
    ! Variation used by bmad_parser2: Check for missing slaves in check_lat instead of lord_lat.

    n_slave = 0

    if (allocated(m_eles)) deallocate (m_eles)
    allocate (m_eles(size(pele%control)))

    do ip = 1, size(pele%control)
      pc => pele%control(ip)
      call lat_ele_locator (pc%name, lat, m_eles(ip)%eles, m_eles(ip)%n_loc, err)
      n_loc = m_eles(ip)%n_loc
      n_slave = n_slave + n_loc

      if (n_loc == 0) then
        if (present(check_lat)) then
          call lat_ele_locator (pc%name, check_lat, in_eles, n_in_loc, err)
        else
          call lat_ele_locator (pc%name, lord_lat, in_eles, n_in_loc, err)
        endif

        if (n_loc == 0 .and. n_in_loc == 0) then
          call parser_error ('CANNOT FIND SLAVE ELEMENT FOR ' // trim(upcase(control_name(lord%lord_status))) // &
                                                   ' ELEMENT: ' // lord%name, 'CANNOT FIND: '// pc%name, pele = pele)
          cycle main_loop
        endif
      endif
    enddo

    if (n_slave == 0) cycle main_loop

    ! Create the lord(s)

    if (is_true(lord%value(gang$))) then
      call make_this_overlay_group_lord(0, lord, lat, n_slave, cs, err_flag, pele, m_eles)
    else
      do ip = 2, size(pele%control)
        if (m_eles(1)%n_loc /= m_eles(ip)%n_loc) then
          call parser_error ('IN OVERLAY OR GROUP ELEMENT: ' // lord%name, &
                    'WITH GANG = FALSE NEED ALL SLAVES WITH A GIVEN NAME TO HAVE THE SAME NUMBER OF', &
                    'ELEMENTS IN THE LATTICE. BUT ' // trim(pele%control(1)%name) // ' HAS \i0\ ELEMENTS', &
                    'WHILE ' // trim(pele%control(ip)%name) // ' HAS \i0\ ELEMENTS', &
                    i_array = [m_eles(1)%n_loc, m_eles(ip)%n_loc])
          cycle main_loop
        endif
      enddo

      do nn = 1, m_eles(1)%n_loc
        call make_this_overlay_group_lord(nn, lord, lat, n_slave, cs, err_flag, pele, m_eles)
        if (err_flag) exit
      enddo
    endif

  !-----------------------------------------------------
  ! girder
  ! Create an girder element for each lattice section that matches the slave list names.
  ! If no lattice element names match any names on the girder slave list then assume 
  ! this girder is for a different lattice and ignore the girder. If some do match and 
  ! some don't then flag this as an error.

  case (girder$)

    if (pele%is_range .and. size(pele%control) /= 2) then
      call parser_error ('GIRDER HAS BAD "ELE_START:ELE_END" RANGE CONSTRUCT. ' // ele%name)
      cycle
    endif

    ! Loop over all elements in the lattice.

    if (allocated(cs)) deallocate(cs)
    allocate (cs(size(pele%control)))

    created_girder_lord = .false.

    branch_loop: do ib = 0, ubound(lat%branch, 1)
      branch => lat%branch(ib)
      ie_start = 1
      ix1_slave = -1

      ! Loop over all possible first elements

      ele_loop: do ie_start = 1, branch%n_ele_track


        ! Loop over girder slave list and see if this section matches.

        ix_ele_now = ie_start
        ix_super_lord_end = -1   ! Not in a super_lord
        ixs = 1                  ! Index of girder slave element we are looking for.
        n_slave = 0              ! Number of actual slaves found.
        matched_to_drift = .false.
        have_ignored_a_drift = .false.
        have_wrapped = .false.

        if (size(pele%control) == 0) then
          call parser_error ('GIRDER DOES NOT HAVE ANY ELEMENTS TO SUPPORT: ' // lord%name)
          cycle main_loop
        endif

        ! loop over all girder slaves and see if lattice eles match slaves.

        slave_loop: do
          if (n_slave > 0 .and. cs(1)%slave%ix_ele == ix1_slave) cycle ele_loop  ! Can happen with superposition
          if (ixs > size(pele%control)) exit

          ! Wrap around the origin if needed.
          if (ix_ele_now > branch%n_ele_track) then
            ix_ele_now = ix_ele_now - branch%n_ele_track
            have_wrapped = .true.
          endif

          if (have_wrapped .and. ix_ele_now == ie_start) then
            call parser_error ('GIRDER SLAVE ELEMENT NOT FOUND: ' // slave_name)
            cycle main_loop
          endif

          slave_name = pele%control(ixs)%name

          ele => pointer_to_ele (lat, ix_ele_now, ib)

          if (girder_match_slave_element(ele, slave_name, ele, n_slave, cs, ix_super_lord_end, ixs, ix_ele_now, pele)) cycle

          ! Here if no match. 
          ! If not previously in a super lord then there is no overall match to the slave list.

          if (ele%slave_status == super_slave$ .and. ix_super_lord_end > -1) then
            if (ix_ele_now == ix_super_lord_end) ix_super_lord_end = -1
            ix_ele_now = ix_ele_now + 1
            cycle
          endif

          ! No match to the slave list. If a marker or drift then ignore except 
          ! if this is the first slave in which case there is no match.

          if ((ele%key == marker$ .or. ele%key == drift$) .and. ixs > 1) then
            if (ix_ele_now == ix_super_lord_end) ix_super_lord_end = -1
            ix_ele_now = ix_ele_now + 1
            if (ele%key == drift$) have_ignored_a_drift = .true.
            cycle
          endif

          ! Match failed. Start again 

          cycle ele_loop

        enddo slave_loop

        if (matched_to_drift .and. have_ignored_a_drift .and. .not. pele%is_range) cycle  ! matching rules violated.

        ! create the girder element

        if (n_slave == 0) then
          call parser_error ('LIST OF GIRDER SLAVES IN LATTICE FILE DOES NOT INCLUDE A NON-DRIFT ELEMENT: ' // &
                              lord%name, level = s_warn$)
          cycle main_loop
        endif

        call new_control (lat, ix_lord, lord%name)
        call create_girder (lat, ix_lord, cs(1:n_slave), lord, err)
        if (err) then
          call parser_error ('ERROR CONSTRUCTING GIRDER')
          lat%n_ele_max = lat%n_ele_max - 1
          return
        endif

        created_girder_lord = .true.
        ix1_slave = cs(1)%slave%ix_ele

      enddo ele_loop

    enddo branch_loop

    if (.not. created_girder_lord) then
      call parser_error ('FAILED TO FIND REGION IN LATTICE FOR CREATING GIRDER: ' // &
                          lord%name, level = s_warn$)
    endif

  end select

enddo main_loop

call control_bookkeeper (lat)

!-------------------------------------------------------------------------
contains

recursive function girder_match_slave_element (ele, slave_name, slave, n_slave, cs, &
                                        ix_super_lord_end, ixs, ix_ele_now, pele) result (is_matched)

type (ele_struct), target :: ele, slave
type (ele_struct), pointer :: lord, slave1, slave2, ele2
type (control_struct), allocatable, target :: cs(:)
type (control_struct), allocatable :: cs_temp(:)
type (parser_ele_struct) :: pele

integer n_slave, ixs, ix_ele_now, ix_super_lord_end
integer ii, ls, n, ie
logical is_matched, add_slave

character(*) slave_name

! If element is already slaved to a girder must match to the girder (an element cannot have multiple girder lords).
! Need to avoid the situation where "gird1:gird2" matches to elements inbetween gird1 and gird2 even if
! the element is slaved to a third girder.

lord => pointer_to_girder(ele)
if (associated(lord)) then
  ele2 => lord
else
  ele2 => ele
endif

! Try to match

is_matched = match_wild(ele2%name, slave_name)

if (is_matched .or. (pele%is_range .and. ixs == 2 .and. (ele2%slave_status == free$ .or. ele2%slave_status == minor_slave$))) then

  if (is_matched) ixs = ixs + 1  ! Next element in list

  if (ele2%key == drift$) then
    if (is_matched) matched_to_drift = .true.

  else ! If not a drift then ele will be a girder slave
    ! First check for duplicates. For example, an element superimposed over a marker
    ! marker will be duplicated unless a check is made.
    add_slave = .true.
    do ie = 1, n_slave
      if (ele2%ix_ele == cs(ie)%slave%ix_ele .and. ele2%ix_branch == cs(ie)%slave%ix_branch) add_slave = .false.
    enddo

    if (add_slave) then
      n_slave = n_slave + 1
      if (size(cs) < n_slave) then  ! Can happen if there is a range.
        n = size(cs)
        call move_alloc(cs, cs_temp)
        allocate (cs(2*n))
        cs(1:n) = cs_temp
      endif
      cs(n_slave)%slave = lat_ele_loc_struct(ele2%ix_ele, ele2%ix_branch)
    endif
  endif

  is_matched = .true.

  ! If a super_lord the logic here is complicated by the fact that 
  ! elements can, for example, be completely contained within another element.
  ! Note: slave can be a super_lord if ele is a multipass_lord

  if (ele2%lord_status == super_lord$) then
    slave2 => pointer_to_slave(ele2, ele2%n_slave)
    ix_super_lord_end = ix_far_index(ix_ele_now-1, ix_super_lord_end, slave2%ix_ele)
  endif

  if (slave%lord_status == super_lord$) then
    slave2 => pointer_to_slave(slave, slave%n_slave)
    ix_super_lord_end = ix_far_index(ix_ele_now-1, ix_super_lord_end, slave2%ix_ele)
  endif

  ! If in a super_slave region then need to recheck the current slave against the next girder slave name.

  if (ix_super_lord_end == -1 .or. pele%is_range) ix_ele_now = ix_ele_now + 1

  ! If match to a girder then move pointers to element after last girder slave

  if (ele2%key == girder$) then
    call find_element_ends (ele2, slave1, slave2)
    ix_ele_now = slave2%ix_ele + 1
    ix_super_lord_end = -1
  endif

  return
endif

! Since ele does not match, look for match at a super or multipass lord of this element.
! Girder lords have been handled above.

do ii = 1, ele%n_lord
  lord => pointer_to_lord (ele, ii)
  ls = lord%lord_status
  if (ls /= super_lord$ .and. ls /= multipass_lord$) cycle
  is_matched = girder_match_slave_element(lord, slave_name, ele, n_slave, cs, ix_super_lord_end, ixs, ix_ele_now, pele)
  if (is_matched) return
enddo

end function girder_match_slave_element

!-------------------------------------------------------------------------
! contains

! Function to return the index that is farthest (reached last) when moving
! from ix_now in a positive direction. Tricky part is if there is wrap around.
! An index of -1 means that the corresponding point does not exist.

function ix_far_index (ix_now, ix1, ix2) result (ix_far)

implicit none

integer ix_now, ix1, ix2, ix_far

!

if (ix1 < 0) then      ! Point 1 does not exist so point 2 wins by default
  ix_far = ix2
elseif (ix2 < 0) then  ! Point 2 does not exist so point 1 wins by default
  ix_far = ix1
elseif (ix1 < ix_now .and. ix2 < ix_now) then  ! both wrapped case
  ix_far = max(ix1, ix2)
elseif (ix1 > ix_now .and. ix2 > ix_now) then ! No wrap ("normal") case
  ix_far = max(ix1, ix2)
else                      ! One is wrapped but not the other case
  ix_far = min(ix1, ix2)
endif

end function ix_far_index

!-------------------------------------------------------------------------
! contains

subroutine make_this_overlay_group_lord (ix_pick, lord, lat, n_slave, cs, err_flag, pele, m_eles)

type (ele_struct) :: lord
type (lat_struct) lat
type (control_struct), allocatable, target :: cs(:)
type (parser_ele_struct), target :: pele
type (multi_ele_pointer_struct), allocatable :: m_eles(:)
type (parser_controller_struct), pointer :: pc
type (all_pointer_struct) a_ptr

integer n_slave
integer k, n, ix, ip, iv, ix_pick, ix_lord
logical err_flag, err
character(40) attrib_name

!

err_flag = .true.

n = n_slave * size(lord%control%var)
if (allocated(cs)) deallocate(cs)
allocate (cs(n))

! Slave setup

n_slave = 0 ! number of slaves found
do ip = 1, size(pele%control)

  pc => pele%control(ip)

  do iv = 1, size(lord%control%var)
    ! Only when attrib_name == '*' is the loop executed multiple times.
    if (pc%attrib_name /= '*' .and. iv > 1) exit

    select case (pc%attrib_name)
    case (blank_name$);   attrib_name = pele%default_attrib
    case ('*');           attrib_name = lord%control%var(iv)%name
    case default;         attrib_name = pc%attrib_name
    end select

    ! There might be more than 1 element with same name. 
    ! Loop over all elements whose name matches name.
    ! Put the info into the cs structure.

    do k = 1, m_eles(ip)%n_loc
      if (ix_pick /= 0 .and. k /= ix_pick) cycle
      slave => pointer_to_ele (lat, m_eles(ip)%eles(k)%loc)
      n_slave = n_slave + 1
      if (allocated(pc%y_knot)) then
        cs(n_slave)%y_knot = pc%y_knot
      else
        call reallocate_expression_stack (cs(n_slave)%stack, pc%n_stk)
        cs(n_slave)%stack = pc%stack(1:pc%n_stk)
      endif
      cs(n_slave)%slave = lat_ele_loc_struct(slave%ix_ele, slave%ix_branch)
      cs(n_slave)%lord%ix_ele = -1             ! dummy value
      call pointer_to_attribute (slave, attrib_name, .true., a_ptr, err, .false., ix_attrib = ix)
      
      ! If attribute not found it may be a special attribute like accordion_edge$.
      ! A special attribute will have ix > num_ele_attrib$
      if (.not. associated(a_ptr%r) .and. lord%key == group$) then
        ix = attribute_index(lord, attrib_name)
        if (ix <= num_ele_attrib$) ix = 0  ! Mark as not valid
      endif
      cs(n_slave)%ix_attrib = ix
      cs(n_slave)%attribute = attrib_name
      if (ix < 1 .and. .not. associated(a_ptr%r)) then
        call parser_error ('IN OVERLAY OR GROUP ELEMENT: ' // lord%name, 'ATTRIBUTE: ' // attrib_name, &
                      'IS NOT A VALID ATTRIBUTE OF: ' // slave%name, pele = pele)
        return
      endif
      if (iv > 1) call parser_transfer_control_struct(cs(n_slave), cs(n_slave), lord, iv)
    enddo
  enddo

enddo

! create the lord

if (n_slave == 0) return   ! If the lord has no slaves then do not create a lord.

call new_control (lat, ix_lord, lord%name)  ! get index in lat where lord goes
lat%ele(ix_lord) = lord

select case (lord%key)
case (overlay$)
  call create_overlay (lat%ele(ix_lord), cs(1:n_slave), err)
case (group$)
  call create_group (lat%ele(ix_lord), cs(1:n_slave), err)
end select
if (err) return

lat%ele(ix_lord)%value(gang$) = lord%value(gang$)

! Finish.

err_flag = .false.

end subroutine make_this_overlay_group_lord

end subroutine parser_add_lord

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine drift_and_pipe_track_methods_adjustment(lat)
!
! Drift and pipe elements can be used in both photon and non-photon lines.
! A problem occures if, for example, a lattice file with both photon and
! non-photon lines contains a line like:
!   drift::*[tracking_method] = taylor
! So this routine resets drift and pipe tracking_method and mat6_calc_method
! parameters in photon lines to bmad_standard if needed.
!
! Input:
!   lat   -- lat_struct: Lattice
!
! Output:
!   lat   -- lat_struct: Lattice with tracking methods adjusted if needed.
!-

subroutine drift_and_pipe_track_methods_adjustment(lat)

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, lord, lord2

integer ib, ie, j

! Only adjust for drift and pipe elements. 
! If there is a problem with other types of elements we don't want to cover up a problem.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  if (branch%param%particle /= photon$) cycle

  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    if (ele%key /= drift$ .and. ele%key /= pipe$) cycle

    if (.not. valid_tracking_method(ele, photon$, ele%tracking_method)) then
      ele%tracking_method = bmad_standard$
      do j = 1, ele%n_lord
        lord => pointer_to_lord(ele, j)
        if (lord%lord_status /= super_lord$ .and. lord%lord_status /= multipass_lord$) cycle
        lord%tracking_method = bmad_standard$
        if (lord%slave_status == multipass_slave$) then
          lord2 => pointer_to_lord(lord, 1)
          lord2%tracking_method = bmad_standard$
        endif
      enddo
    endif

    if (.not. valid_mat6_calc_method(ele, photon$, ele%mat6_calc_method)) then
      ele%mat6_calc_method = bmad_standard$
      do j = 1, ele%n_lord
        lord => pointer_to_lord(ele, j)
        if (lord%lord_status /= super_lord$ .and. lord%lord_status /= multipass_lord$) cycle
        lord%mat6_calc_method = bmad_standard$
        if (lord%slave_status == multipass_slave$) then
          lord2 => pointer_to_lord(lord, 1)
          lord2%mat6_calc_method = bmad_standard$
        endif
      enddo
    endif
  enddo
enddo

end subroutine drift_and_pipe_track_methods_adjustment

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine settable_dep_var_bookkeeping (ele)
!
! Subroutine to initialize dependent variables in an element.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine settable_dep_var_bookkeeping (ele)

use random_mod

implicit none

type (ele_struct),  target :: ele
type (branch_struct), pointer :: branch
type (surface_grid_struct), pointer :: grid
type (surface_grid_pt_struct), pointer :: pt

real(rp) a, rr, v_inv_mat(4,4), eta_vec(4)

integer n, i, j, i0, i1, j0, j1
logical kick_set, length_set, set_done, err_flag
logical b_field_set, g_set

! Wall3d init.

if (associated(ele%wall3d)) then
  do n = 1, size(ele%wall3d)
    call wall3d_initializer (ele%wall3d(n), err_flag)
    if (err_flag) then
      call parser_error ('WALL INIT ERROR FOR ELEMENT: ' // ele%name)
      return
    endif
  enddo
endif

! Surface init

if (associated(ele%photon)) then
  grid => ele%photon%grid
  select case (grid%type)
  case (segmented$)
    do i = lbound(grid%pt, 1), ubound(grid%pt, 1)
    do j = lbound(grid%pt, 2), ubound(grid%pt, 2)
      call init_surface_segment (ele%photon, i, j)
    enddo
    enddo

  case (displacement$)
    do i = lbound(grid%pt, 1), ubound(grid%pt, 1)
    do j = lbound(grid%pt, 2), ubound(grid%pt, 2)
      pt => grid%pt(i,j)

      pt%x0 = i * grid%dr(1) + grid%r0(1)
      pt%y0 = j * grid%dr(2) + grid%r0(2)

      i0 = i - 1; i1 = i + 1
      j0 = j - 1; j1 = j + 1
      if (i == lbound(grid%pt, 1)) i0 = i
      if (i == ubound(grid%pt, 1)) i1 = i
      if (j == lbound(grid%pt, 2)) j0 = j
      if (j == ubound(grid%pt, 2)) j1 = j
      
      if (pt%dz_dx == real_garbage$) then
        pt%dz_dx = (grid%pt(i1,j)%z0 - grid%pt(i0,j)%z0) / ((i1-i0)*grid%dr(1))
        pt%dz_dy = (grid%pt(i,j1)%z0 - grid%pt(i,j0)%z0) / ((j1-j0)*grid%dr(2))
      endif

      if (pt%d2z_dxdy == real_garbage$) then
        pt%d2z_dxdy = (grid%pt(i1,j1)%z0 - grid%pt(i1,j0)%z0 - grid%pt(i0,j1)%z0 + grid%pt(i0,j0)%z0) / &
                                                            ((i1-i0)*grid%dr(1)*(j1-j0)*grid%dr(2))
      endif
    enddo
    enddo

  end select
endif

! Aperture init

if (ele%aperture_type == auto_aperture$) then
  call aperture_bookkeeper (ele)
endif

! Note: If an attribute has a value of real_garbage$ then that attribute has not been set by the user.

kick_set = (ele%value(hkick$) /= 0) .or. (ele%value(vkick$) /= 0)

select case (ele%key)

! Taylor map gets unit spin quaternion if nothing has been set for spin.

case (taylor$)

  if (size(ele%spin_taylor(0)%term) == 0 .and. size(ele%spin_taylor(1)%term) == 0 .and. &
                size(ele%spin_taylor(2)%term) == 0 .and. size(ele%spin_taylor(3)%term) == 0) then
    call add_taylor_term(ele%spin_taylor(0), 1.0_rp, [0, 0, 0, 0, 0, 0])
  endif

!------------------
case (beginning_ele$)

  if (ele%a%beta /= 0) ele%a%gamma = (1 + ele%a%alpha**2) / ele%a%beta
  if (ele%b%beta /= 0) ele%b%gamma = (1 + ele%b%alpha**2) / ele%b%beta

  ele%gamma_c = sqrt(1 - ele%c_mat(1,1)*ele%c_mat(2,2) + ele%c_mat(1,2)*ele%c_mat(2,1))

  call make_v_mats (ele, v_inv_mat = v_inv_mat)
  eta_vec = matmul (v_inv_mat, [ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap])

  ele%a%eta  = eta_vec(1)
  ele%a%etap = eta_vec(2)
  ele%b%eta  = eta_vec(3)
  ele%b%etap = eta_vec(4)

!------------------
! Convert rbends to sbends and evaluate G if needed.
! Needed is the length and either: angle, G, or rho.
! Note: l -> l_chord for rbends has already been done.

case (sbend$, rbend$, rf_bend$) 

  b_field_set = (ele%value(b_field$) /= 0 .or. ele%value(db_field$) /= 0)
  g_set = (ele%value(g$) /= 0 .or. ele%value(dg$) /= 0)
  if ((ele%value(angle$) /= 0 .or. ele%value(rho$) /= 0 .or. g_set) .and. &
                                      .not. b_field_set) ele%value(b_field$) = real_garbage$

  if (ele%key /= rf_bend$) ele%sub_key = ele%key  ! Save sbend/rbend input type.

  ! Only one of b_field, g, or rho may be set.
  ! B_field may not be set for an rbend since, in this case, L is not computable (we don't know the ref energy).

  if (b_field_set .and. ele%key == rbend$) call parser_error &
          ("B_FIELD NOT SETTABLE FOR AN RBEND (USE AN SBEND INSTEAD): " // ele%name)

  if (b_field_set .and. g_set) call parser_error &
          ('BOTH G (OR DG) AND B_FIELD (OR DB_FIELD) SET FOR A BEND: ' // ele%name)

  if (b_field_set .and. ele%value(rho$) /= 0) call parser_error &
          ('BOTH RHO AND B_FIELD (OR DB_FIELD) SET FOR A BEND: ' // ele%name)

  if (ele%value(g$) /= 0 .and. ele%value(rho$) /= 0) &
            call parser_error ('BOTH G AND RHO SPECIFIED FOR BEND: ' // ele%name)

  ! if rho is set then this gives g

  if (ele%value(l$) /= 0 .and. ele%value(angle$) /= 0 .and. ele%value(g$) /= 0) &
                      call parser_error ('ANGLE, G, AND L ARE ALL SPECIFIED FOR BEND: ' // ele%name)
  if (ele%value(l$) /= 0 .and. ele%value(angle$) /= 0 .and. ele%value(rho$) /= 0) &
                      call parser_error ('ANGLE, RHO, AND L ARE ALL SPECIFIED FOR BEND: ' // ele%name)

  if (ele%value(rho$) /= 0) ele%value(g$) = 1 / ele%value(rho$)

  ! If g and angle are set then this determines l

  if (ele%value(g$) /= 0 .and. ele%value(angle$) /= 0) ele%value(l$) = ele%value(angle$) / ele%value(g$)

  if (ele%value(angle$) /= 0 .and. ele%value(l$) == 0 .and. ele%value(l_chord$) == 0) then
    call parser_error ('THE BENDING ANGLE IS NONZERO IN A ZERO LENGTH BEND! ' // ele%name)
  endif

  ! Convert an rbend to an sbend

  if (ele%key == rbend$) then
    ! Note: L must be zero if g and angle have both been specified and are non-zero.
    if (ele%value(l$) == 0 .and. ele%value(l_chord$) /= 0) then
      if (ele%value(angle$) /= 0) then
        ele%value(l$) = ele%value(l_chord$) * ele%value(angle$) / (2.0_rp * sin(0.5_rp*ele%value(angle$)))
      elseif (ele%value(g$) /= 0) then
        a = 0.5_rp * ele%value(l_chord$) * ele%value(g$)
        if (abs(a) >= 1) then
          call parser_error ('G * L FOR RBEND IS TOO LARGE TO BE PHYSICAL! ' // ele%name)
          return
        endif
        a = 2 * asin(a)
        ele%value(l$) = ele%value(l_chord$) * a / (2.0_rp * sin(0.5_rp*a))
      else  ! g and angle are zero.
        ele%value(l$) = ele%value(l_chord$)
      endif
    endif

    if (ele%value(l$) /= 0 .and. ele%value(angle$) /= 0) then
      ele%value(g$) = ele%value(angle$) / ele%value(l$) 
    elseif (ele%value(g$) /= 0) then
      ele%value(angle$) = ele%value(g$) * ele%value(l$) 
    endif

    ele%value(e1$) = ele%value(e1$) + 0.5_rp * ele%value(angle$)
    ele%value(e2$) = ele%value(e2$) + 0.5_rp * ele%value(angle$)
    ele%key = sbend$
  endif

  ! 

  if (ele%value(angle$) /= 0) ele%value(g$) = ele%value(angle$) / ele%value(l$) 

  ! If fintx or hgapx are real_garbage then they have not been set.
  ! If so, set their valuse to fint and hgap.

  if (ele%key /= rf_bend$) then
    if (ele%value(hgapx$) == real_garbage$) ele%value(hgapx$) = ele%value(hgap$)
    if (ele%value(fintx$) == real_garbage$) ele%value(fintx$) = ele%value(fint$)
  endif

!------------------
! Accept use of Voltage for lcavities and vary the mode frequencies.

case (lcavity$) 

  if (ele%value(voltage$) /= 0) then
    if (ele%value(gradient$) /= 0) call parser_error &
                ('BOTH VOLTAGE AND GRADIENT NON-ZERO FOR A LCAVITY:', ele%name)
    if (ele%value(l$) == 0) then
      ele%value(gradient$) = 0
    else
      ele%value(gradient$) = ele%value(voltage$) / ele%value(l$)
    endif
  endif

  if (ele%value(voltage_err$) /= 0) then
    if (ele%value(gradient_err$) /= 0) call parser_error &
                ('BOTH VOLTAGE_ERR AND GRADIENT_ERR NON-ZERO FOR A LCAVITY:', ele%name)
    if (ele%value(l$) == 0) then
      ele%value(gradient_err$) = 0
    else
      ele%value(gradient_err$) = ele%value(voltage_err$) / ele%value(l$)
    endif
  endif

!------------------

case (crab_cavity$)

  if (ele%value(voltage$) /= 0) then
    if (ele%value(gradient$) /= 0) call parser_error &
                ('BOTH VOLTAGE AND GRADIENT NON-ZERO FOR A CRAB_CAVITY:', ele%name)
    if (ele%value(l$) == 0) then
      ele%value(gradient$) = 0
    else
      ele%value(gradient$) = ele%value(voltage$) / ele%value(l$)
    endif
  endif

!------------------

case (rfcavity$) 

  if (ele%value(rf_frequency$) /= 0 .and. ele%value(harmon$) /= 0) call parser_error &
              ('BOTH RF_FREQUENCY AND HARMON SET FOR RFCAVITY: ' // ele%name, &
               'SETTING OF HARMON WILL BE IGNORED!', level = s_warn$)

  if (ele%value(voltage$) /= 0) then
    if (ele%value(gradient$) /= 0) call parser_error &
                ('BOTH VOLTAGE AND GRADIENT NON-ZERO FOR A RFCAVITY:', ele%name)
    if (ele%value(l$) == 0) then
      ele%value(gradient$) = 0
    else
      ele%value(gradient$) = ele%value(voltage$) / ele%value(l$)
    endif
  endif


!------------------
! for a planar and helical wigglers n_pole is a dependent attribute

case (wiggler$, undulator$)
  if (ele%field_calc == int_garbage$) ele%field_calc = planar_model$

  ! Default tracking_method for map type elements is Taylor. Bmad_standard is not possible.
  if (ele%field_calc /= planar_model$ .and. ele%field_calc /= helical_model$) then
    if (ele%tracking_method == bmad_standard$) ele%tracking_method = taylor$
    if (ele%mat6_calc_method == bmad_standard$) ele%mat6_calc_method = taylor$
  endif

  if (ele%value(l_period$) == 0 .and. ele%value(n_period$) /= 0) then
    ele%value(l_period$) = ele%value(l$) / ele%value(n_period$) 
  endif

!------------------
case (match$)
  ele%value(phase_trombone_input$)  = ele%value(phase_trombone$)
  ele%value(match_end_input$)       = ele%value(match_end$)
  ele%value(match_end_orbit_input$) = ele%value(match_end_orbit$)

!------------------

case (quadrupole$)
  if (ele%field_master .and. (ele%value(k1$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM FOR ELEMENT: ' // ele%name, &
       'BOTH STRENGTH (K1, HKICK, ETC.) AND FIELD SET FOR A QUAD.')

!------------------

case (solenoid$)
  if (ele%field_master .and. (ele%value(ks$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM FOR ELEMENT: ' // ele%name, &
       'BOTH STRENGTH (KS, HKICK, ETC.) AND FIELD SET FOR A SOLENOID.')

!------------------

case (sol_quad$)
  if (ele%field_master .and. (ele%value(ks$) /= 0 .or. &
                            ele%value(k1$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM FOR ELEMENT: ' // ele%name, &
       'BOTH STRENGTH (K1, HKICK, ETC.) AND FIELD SET FOR A SOL_QUAD.')

!------------------

case (sextupole$)
  if (ele%field_master .and. (ele%value(k2$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM FOR ELEMENT: ' // ele%name, &
       'BOTH STRENGTH (K2, HKICK, ETC.) AND FIELD SET FOR A SEXTUPOLE.')

!------------------

case (octupole$)
  if (ele%field_master .and. (ele%value(k3$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM FOR ELEMENT: ' // ele%name, &
       'BOTH STRENGTH (K3, HKICK, ETC.) AND FIELD SET FOR A OCTUPOLE.')

!------------------

case (hkicker$, vkicker$)
  if (ele%field_master .and. (ele%value(kick$) /= 0 .or. kick_set)) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM FOR ELEMENT: ' // ele%name, &
       'BOTH STRENGTH AND BL_KICK SET FOR A H/VKICKER.')

!------------------

case (elseparator$)
  if (ele%field_master .and. kick_set) call parser_error &
      ('INDEPENDENT VARIABLE PROBLEM FOR ELEMENT: ' // ele%name, &
       'BOTH KICK (HKICK OR VKICK) AND E_FIELD OR VOLTAGE SET FOR A ELSEPARATOR.')

  if (ele%field_master) then
    if (ele%value(voltage$) /= 0 .and. ele%value(e_field$) /= 0) call parser_error &
              ('INDEPENDENT VARIABLE PROBLEM FOR ELSEPARATOR: ' // ele%name, &
               'BOTH VOLTAGE AND E_FIELD SET FOR THIS ELEMENT.')

    if (ele%value(voltage$) /= 0) then
      if (ele%value(gap$) == 0) then
        call parser_error ('FOR ELSEPARATOR: ' // ele%name, 'VOLTAGE IS SET BUT GAP IS NOT!')
      else
        ele%value(e_field$) = ele%value(voltage$) / ele%value(gap$)
      endif
    endif
  endif

!------------------

case (e_gun$)
  if (ele%value(gradient$) == 0 .and. ele%value(l$) /= 0) ele%value(gradient$) = ele%value(voltage$) / ele%value(l$)

end select

end subroutine settable_dep_var_bookkeeping 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine form_digested_bmad_file_name (lat_file, digested_file, full_lat_file, use_line)
!
! Subroutine to form the standard name of the Bmad digested file. 
! The standard digested file name has the suffix added to the file name:
!     suffix = '.digested' + bmad_inc_version$ 
! Exception: If the use_line argument is present and not blank, the suffix will be:
!     suffix = '.' + use_line + '.digested' + bmad_inc_version$ 
!   
!
! Input:
!   lat_file -- Character(*): Input lattice file name.
!   use_line -- Character(*), optional: Line used for lattice expansion. If not present
!                 or blank, the line used is the one that was specified in the lattice file.
!
! Output:
!   digested_file -- Character(200): Name of the digested file.
!   full_lat_file -- Character(200), optional: Input lattice file name with full directory.
!                       Can be used for error messages.
!-

subroutine form_digested_bmad_file_name (lat_file, digested_file, full_lat_file, use_line)

character(*) lat_file, digested_file
character(*), optional :: full_lat_file, use_line
character(200) name, full_name

integer ix

! Singular case

if (lat_file == '') then
  digested_file = ''
  return
endif

! Get the full_lat_file name

call fullfilename (lat_file, name)
inquire (file = name, name = full_name)  ! full input file_name
if (present (full_lat_file)) full_lat_file = full_name

! Construct the digested_file name

if (present(use_line)) then
  if (use_line /= '') then
    write (digested_file, '(4a, i0)') trim(full_name), '.', trim(use_line), '.digested', bmad_inc_version$ 
    return
  endif
endif

write (digested_file, '(2a, i0)') trim(full_name), '.digested', bmad_inc_version$ 

end subroutine form_digested_bmad_file_name

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_add_branch (fork_ele, lat, sequence, in_name, in_indexx, &
!                                  seq_name, seq_indexx, no_end_marker, in_lat, plat, created_new_branch)
!
! Subroutine to do line expansion.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_add_branch (fork_ele, lat, sequence, seq_name, &
                                    seq_indexx, no_end_marker, in_lat, plat, created_new_branch, new_branch_name)

implicit none

type (lat_struct), target :: lat, in_lat
type (parser_lat_struct) plat
type (ele_struct) fork_ele
type (ele_struct), pointer :: target_ele
type (seq_struct), allocatable, target :: sequence(:)
type (branch_struct), pointer :: branch

integer, allocatable :: seq_indexx(:), in_indexx(:)
integer i, j, nb, n_ele_use, n, ix, key

character(*), optional :: new_branch_name
character(*), allocatable ::  seq_name(:)
character(40) branch_name

logical created_new_branch, no_end_marker

!

created_new_branch = .true.
nb = ubound(lat%branch, 1)

if (is_false(fork_ele%value(new_branch$))) then ! Branch back if
  do i = 0, nb - 1
    branch => lat%branch(i)
    if (branch%name /= fork_ele%component_name) cycle
    fork_ele%value(ix_to_branch$) = i
    created_new_branch = .false.
    if (present(new_branch_name)) new_branch_name = ''
  enddo
endif

branch_name = fork_ele%component_name
if (present(new_branch_name)) new_branch_name = branch_name
fork_ele%component_name = plat%ele(fork_ele%ixx)%ele_name  ! Substitute element name for line name.

if (created_new_branch) then
  call parser_expand_line (1, branch_name, sequence, seq_name, seq_indexx, no_end_marker, n_ele_use, lat, in_lat)
  if (bp_com%error_flag) return
  nb = nb + 1
  fork_ele%value(ix_to_branch$) = nb
  branch => lat%branch(nb)

  branch%ix_from_branch     = fork_ele%ix_branch
  branch%ix_from_ele        = fork_ele%ix_ele
  n = branch%ix_from_branch

  ! Need to know if the reference energy needs to be specified for the new branch which depends upon if the 
  ! "to" ele is the Beginning element. But the "to" element may not yet exist (may be superimposed later).
  ! So just do something crude for now.
  if (fork_ele%component_name == '' .or. fork_ele%component_name == 'BEGINNING') then
    branch%ix_to_ele = 0
  else
    branch%ix_to_ele = -1
  endif
endif


end subroutine parser_add_branch

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_identify_fork_to_element (lat)
!
! Routine to identify the elements the forks in a lattice are branching to.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_identify_fork_to_element (lat)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: target_ele
type (ele_struct), pointer :: fork_ele
type (branch_struct), pointer :: branch
type (ele_pointer_struct), allocatable :: eles(:)

integer ib, ie, j, n_loc
logical err

character(40) name

!

do ib = 0, ubound(lat%branch, 1)

  do ie = 1, lat%branch(ib)%n_ele_max

    fork_ele => lat%branch(ib)%ele(ie)
    if (fork_ele%key /= fork$ .and. fork_ele%key /= photon_fork$) cycle

    branch => lat%branch(nint(fork_ele%value(ix_to_branch$)))
    nullify(target_ele)
    name = fork_ele%component_name

    if (name == '') then
      if (nint(fork_ele%value(direction$)) == 1) then
        target_ele => branch%ele(0)
      else
        target_ele => branch%ele(branch%n_ele_track)
      endif

    else
      call lat_ele_locator (name, lat, eles, n_loc, err, ix_dflt_branch = branch%ix_branch)
      if (n_loc > 1) then
        call parser_error('DUPLICATE TO_ELEMENT: ' // name, 'FOR FORK ELEMENT: ' // fork_ele%name)
        return
      endif

      if (n_loc == 0) then
        call parser_error('TO_ELEMENT NOT FOUND: ' // name, 'FOR FORK ELEMENT: ' // fork_ele%name)
        return
      endif

      target_ele => eles(1)%ele
    endif


    fork_ele%value(ix_to_element$) = target_ele%ix_ele

    select case (target_ele%key)
    case (marker$, fork$, photon_fork$, fiducial$, beginning_ele$)
    case default
      call parser_error('TO_ELEMENT: ' // name, 'FOR FORK ELEMENT: ' // fork_ele%name, &
                        'IS NOT A ZERO-LENGTH MARKER-LIKE ELEMENT')
    end select

  enddo

  branch => lat%branch(ib)
  if (branch%ix_from_branch > -1) then
    fork_ele => lat%branch(branch%ix_from_branch)%ele(branch%ix_from_ele)
    branch%ix_to_ele = nint(fork_ele%value(ix_to_element$))
  endif

enddo

end subroutine parser_identify_fork_to_element

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_expand_line (lat, line_name, sequence, &
!               seq_name, seq_indexx, no_end_marker, n_ele_expand, lat, in_lat, expanded_line)
!
! Subroutine to do line expansion.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!
! Note: Either lat and in_lat must be present or expanded_line must be present.
!
! Input:
!   i_lev         -- integer: Subsequence level. 1 => Root level.
!   line_name     -- character(*): Root line to expand.
!   sequence(:)   -- seq_struct: Array of sequencies.
!   seq_name(:)   -- character(*): Array of sequence names.
!   seq_indexx(:) -- integer: Index array for the sequence names.
!   no_end_marker -- logical: Put a marker named "end" at the end of the branch?
!   lat           -- lat_struct, optional: Lattice to put the expanded line
!   in_lat        -- lat_struct, optional: Lattice with array of defined elements.
!
! Output:
!   n_ele_expand     -- integer: Number of elements in the finished line.
!   lat              -- lat_struct, optional: Lattice with new line. Except if expanded_line is present.
!   expanded_line(:) -- base_line_ele_struct, optional: If present, lat argument will be
!                         ignored and the expanded line will be put into expanded_line.
!-

recursive subroutine parser_expand_line (i_lev, line_name, sequence, &
               seq_name, seq_indexx, no_end_marker, n_ele_expand, lat, in_lat, expanded_line)

implicit none

type (lat_struct), optional, target :: lat, in_lat
type (ele_struct), pointer :: ele_line(:), ele, ele2
type (seq_struct), allocatable, target :: sequence(:)
type (seq_ele_struct), pointer :: s_ele, this_seq_ele
type (seq_ele_struct), target :: dummy_seq_ele
type (seq_struct), pointer :: seq, seq2
type (base_line_ele_struct), allocatable, target ::  base_line(:), sub_line(:)
type (base_line_ele_struct), optional, allocatable :: expanded_line(:)
type (base_line_ele_struct), pointer :: b_ele
type (branch_struct), pointer :: branch

integer, allocatable :: seq_indexx(:)
integer iseq_tot, i_lev, ix_seq, n_ele_expand, rep_count, ix_ele, n_ele2
integer i, j, k, n, ix, ix_multipass, ix_branch, flip, ixs_start, ixs_end

character(*), allocatable ::  seq_name(:)
character(*) line_name
character(40) name
character(40), allocatable :: my_line(:)
character(100) err_line2

logical no_end_marker

! find line corresponding to the "use" statement.

iseq_tot = size(seq_indexx)
allocate (base_line(100))
err_line2 = ''
if (bp_com%detected_expand_lattice_cmd) err_line2 = &
                      'NOTE: ERROR HAPPENS AFTER AN EXPAND_LATTICE COMMAND. THIS MAY BE OF SIGNIFICANCE.'

call find_index (line_name, seq_name, seq_indexx, iseq_tot, ix_seq)
if (ix_seq == 0) then
  if (i_lev == 1) then
    call parser_error ('CANNOT FIND DEFINITION OF LINE IN "USE" STATEMENT: ' // line_name, err_line2, stop_here = .true.)
  else
    call parser_error ('CANNOT FIND DEFINITION OF SUBLINE: ' // line_name, err_line2, stop_here = .true.)
  endif
  return
endif

if (i_lev == 1 .and. sequence(ix_seq)%type /= line$) then
  call parser_error ('NAME IN "USE" STATEMENT IS NOT A LINE!', stop_here = .true.)
  return
endif

seq => sequence(ix_seq)
rep_count = seq%ele(1)%rep_count
ix_ele =  1              ! we start at the beginning
n_ele_expand = 0
if (seq%active) then
  call parser_error('INFINITE RECURSION: A SUBLINE OF LINE: ' // trim(line_name) // ' IS THIS SAME LINE!', stop_here = .true.)
  return
endif
seq%active = .true.
         
sequence(:)%ix_list = 1  ! Init. Used for replacement list index
sequence(:)%list_upcount = 0 ! Count up

!-------------------------------------------------------------------------
! Expand base line...

line_expansion: do

  ! If rep_count is zero then we are finished with the current element.

  if (rep_count == 0) then      ! goto next element in the sequence
    ! Goto the next element 
    ix_ele = ix_ele + 1

    ! Check if off the end of the current line...
    if (ix_ele <= size(seq%ele)) then  ! Nope. Still have more element to process
      rep_count = seq%ele(ix_ele)%rep_count           ! set the rep_count for the next ele.
      if (rep_count == 0) cycle                       ! For "0*sub_line" construct.

    ! If we have got to the end of the current line then pop the stack back to
    ! the next lower level.
    else
      exit
    endif
  endif

  rep_count = rep_count - 1

  ! if s_ele is a dummy arg then get corresponding actual arg.

  s_ele => seq%ele(ix_ele)  ! next element, line, or list

  ix = s_ele%ix_arg
  if (ix /= 0) then  ! it is a dummy argument.
    name = seq%corresponding_actual_arg(ix)
    s_ele => dummy_seq_ele
    s_ele%name = name
    s_ele%ele_orientation = seq%ele(ix_ele)%ele_orientation

    call find_index (name, seq_name, seq_indexx, iseq_tot, j)
    if (j > 0) then  ! if a sequence
      s_ele%ix_ele = j
      s_ele%type = sequence(j)%type
    else  ! Must be an element
      call find_index (name, in_lat%nametable, j)
      if (j == 0) then  ! if not an element then I don't know what it is
        call parser_error ('CANNOT FIND DEFINITION FOR: ' // name, &
                          'IN LINE: ' // seq%name, err_line2, seq = seq)
        if (global_com%exit_on_error) call err_exit
        return
      endif
      s_ele%ix_ele = j 
      s_ele%type = element$
    endif
    
  endif

  ! Select type

  select case (s_ele%type)

  ! If an element

  case (element$, list$) 
    if (s_ele%type == list$) then
      seq2 => sequence(s_ele%ix_ele)
      j = seq2%ix_list
      this_seq_ele => seq2%ele(j)
      seq2%list_upcount = seq2%list_upcount + 1
      if (seq2%list_upcount == this_seq_ele%rep_count) then
        seq2%ix_list = seq2%ix_list + 1
        if (seq2%ix_list > size(seq2%ele(:))) seq2%ix_list = 1
        seq2%list_upcount = 0
      endif
    else
      if (s_ele%tag /= '') then
        call parser_error ('ELEMENTS IN A LINE OR LIST ARE NOT ALLOWED TO HAVE A TAG.', &
                      'FOUND ILLEGAL TAG FOR ELEMENT: ' // s_ele%name, &
                      'IN THE LINE/LIST: ' // seq%name, seq = seq)
      endif
      this_seq_ele => s_ele
    endif

    if (this_seq_ele%ix_ele < 1) call parser_error('NOT A DEFINED ELEMENT: ' // &
                          s_ele%name, 'IN THE LINE/LIST: ' // seq%name, err_line2, seq = seq)

    
    if (n_ele_expand+10 > size(base_line)) call reallocate_base_line(base_line, 2*n_ele_expand)

    n_ele_expand = n_ele_expand + 1
    base_line(n_ele_expand)%ix_ele_in_in_lat = this_seq_ele%ix_ele

    base_line(n_ele_expand)%name = this_seq_ele%name
    base_line(n_ele_expand)%orientation = this_seq_ele%ele_orientation
    base_line(n_ele_expand)%ele_order_reflect = (this_seq_ele%ele_orientation == -1)

  ! if a line:
  !     a) move pointer on current level past line element
  !     b) go to the next higher level
  !     c) initialize pointers for the higher level to use the line

  case (line$, replacement_line$)
    if (i_lev > 100) then
      call parser_error ('NESTED LINES EXCEED STACK DEPTH! SUSPECT INFINITE LOOP!')
      if (global_com%exit_on_error) call err_exit
      return
    endif

    seq2 => sequence(s_ele%ix_ele)
    if (s_ele%type == replacement_line$) then
      if (size(seq2%dummy_arg) /= size(s_ele%actual_arg)) then
        call parser_error ('WRONG NUMBER OF ARGUMENTS FOR REPLACEMENT LINE: ' // &
            s_ele%name, 'WHEN USED IN LINE: ' // seq%name, seq = seq)
        return
      endif
      arg_loop: do i = 1, size(seq2%dummy_arg)
        seq2%corresponding_actual_arg(i) = s_ele%actual_arg(i)
        if (allocated(seq%dummy_arg)) then
          do j = 1, size(seq%dummy_arg)
            if (seq2%corresponding_actual_arg(i) == seq%dummy_arg(j)) then
              seq2%corresponding_actual_arg(i) = seq%corresponding_actual_arg(j)
              cycle arg_loop
            endif
          enddo
        endif
        name = seq2%corresponding_actual_arg(i)
      enddo arg_loop
    endif

    call parser_expand_line (i_lev+1, s_ele%name, sequence, seq_name, &
                                            seq_indexx, no_end_marker, n_ele2, lat, in_lat, sub_line)
    if (bp_com%fatal_error_flag) return
    ixs_start = find_slice_edge(s_ele%slice_start, 1, sub_line(1:n_ele2), s_ele)
    ixs_end = find_slice_edge(s_ele%slice_end, n_ele2, sub_line(1:n_ele2), s_ele)
    if (ixs_start > ixs_end) then
      call parser_error ('FOR SLICE OF LINE: ' // s_ele%name, &
                         'STARTING SLICE POSITION AT: ' // s_ele%slice_start, &
                         'IS PAST ENDING SLICE POSITION AT: ' // s_ele%slice_end)
      ixs_start = ixs_end-1  ! Just to be able to limp along.
    endif
    n_ele2 = ixs_end + 1 - ixs_start
    sub_line(1:n_ele2) = sub_line(ixs_start:ixs_end)

    call reallocate_base_line (base_line, 2*(n_ele_expand+n_ele2))

    do i = 1, n_ele2
      b_ele => base_line(i+n_ele_expand)
      if (s_ele%ele_order_reflect) then
        b_ele = sub_line(n_ele2-i+1)
      else
        b_ele = sub_line(i)
      endif

      b_ele%orientation = b_ele%orientation * s_ele%ele_orientation
      b_ele%ele_order_reflect = (b_ele%ele_order_reflect .neqv. s_ele%ele_order_reflect)

      if (seq2%multipass .and. b_ele%ix_multi == 0) b_ele%ix_multi = i + 1000000 * seq2%index

      if (b_ele%tag /= '' .and. s_ele%tag /= '') then
        b_ele%tag =  trim(s_ele%tag) // '.' // trim(b_ele%tag)
      elseif (s_ele%tag /= '') then
        b_ele%tag = s_ele%tag
      endif
    enddo

    n_ele_expand = n_ele_expand + n_ele2

  case default
    call parser_error ('INTERNAL SEQUENCE ERROR!')

  end select

enddo line_expansion

! Transfer the ele information from the in_lat to lat and
! do the bookkeeping for settable dependent variables.

if (present(expanded_line)) then  ! Used by girders and sublines.
  if (allocated(expanded_line)) deallocate(expanded_line)
  allocate(expanded_line(n_ele_expand))
  expanded_line = base_line(1:n_ele_expand)
  seq%active = .false.
  return
endif

! Stop here if there has been an error

if (bp_com%error_flag) return

!

if (lat%n_ele_max == 0) then
  ix_branch = 0
else
  ix_branch = ubound(lat%branch, 1) + 1
endif

if (ix_branch == 0) then  ! Main branch
  call allocate_lat_ele_array(lat, n_ele_expand+1)
  lat%n_ele_track = n_ele_expand
  lat%n_ele_max   = n_ele_expand
  ele_line => lat%ele
else                    ! branch line
  call allocate_branch_array (lat, ix_branch)
  call allocate_lat_ele_array(lat, n_ele_expand+1, ix_branch)
  ele_line => lat%branch(ix_branch)%ele
endif

branch => lat%branch(ix_branch)

do i = 1, n_ele_expand
  ele_line(i) = in_lat%ele(base_line(i)%ix_ele_in_in_lat) 
  ele_line(i)%name              = base_line(i)%name
  ele_line(i)%iyy               = base_line(i)%ix_multi
  ele_line(i)%orientation       = base_line(i)%orientation
  ele_line(i)%select            = base_line(i)%ele_order_reflect
  ele_line(i)%lord_status       = not_a_lord$  ! In case element is also being superimposed.
  if (base_line(i)%tag /= '') ele_line(i)%name = trim(base_line(i)%tag) // '.' // ele_line(i)%name
  call settable_dep_var_bookkeeping (ele_line(i))
enddo

call init_ele(ele_line(0), beginning_ele$, ix_branch, 0, lat%branch(ix_branch))
call set_ele_defaults (ele_line(0))   ! Defaults for beginning_ele element
ele_line(0)%name = 'BEGINNING'
ele_line(0)%orientation = ele_line(1)%orientation

deallocate(base_line)

! Add End marker and make sure it's orientation is consistant

if (.not. no_end_marker) then
  n_ele_expand = n_ele_expand + 1
  ele => ele_line(n_ele_expand)
  ele%name = 'END'
  ele%key = marker$
  call set_ele_defaults (ele)
  flip = 1
  do j = n_ele_expand-1, 0, -1
    ele2 => ele_line(j)
    if (ele2%key == patch$ .or. ele2%key == floor_shift$) then
      if (patch_flips_propagation_direction (ele2%value(x_pitch$), ele2%value(y_pitch$))) flip = -flip
      cycle
    endif
    exit
  enddo
  ele%orientation = ele2%orientation * flip
endif

! Branch info

branch%n_ele_track = n_ele_expand
branch%n_ele_max   = n_ele_expand
branch%ix_branch   = ix_branch
branch%name        = line_name

seq%active = .false.

!-------------------------------------------------
contains

subroutine reallocate_base_line(base_line, n_ele)
type (base_line_ele_struct), allocatable :: base_line(:), base_temp(:)
integer n_ele, n

!

n = size(base_line)
if (n_ele <= n) return

call move_alloc (base_line, base_temp)
allocate (base_line(1:n_ele))
base_line(1:n) = base_temp(1:n)
deallocate (base_temp)

end subroutine reallocate_base_line

!-------------------------------------------------
! contains

function find_slice_edge(slice_edge, ix_default, this_line, s_ele) result (ix_slice)

type (base_line_ele_struct) ::  this_line(:)
type (seq_ele_struct) s_ele
integer i, ix, ix_default, ix_slice, n_count
character(*) slice_edge
character(40) name

!

ix_slice = ix_default
if (slice_edge == '') return

ix = index(slice_edge, '##')
if (ix == 0) then
  n_count = 1
  name = slice_edge
else
  if (.not. is_integer(slice_edge(ix+2:), n_count)) then
    call parser_error ('EXPECTING INTEGER AFTER SLICE EDGE NAME: ' // slice_edge, &
                       'FOR SLICE: ' // trim(s_ele%name) // '[' // trim(s_ele%slice_start) // ':' // trim(s_ele%slice_end) // ']')
    return
  endif
  name = slice_edge(:ix-1)
endif

do i = 1, size(this_line)
  if (this_line(i)%name /= name) cycle
  n_count = n_count - 1
  if (n_count > 0) cycle
  ix_slice = i
  return
enddo

call parser_error ('CANNOT FIND SLICE EDGE NAME IN LIST OF LINE ELEMENT: ' // slice_edge, &
                   'FOR SLICE: ' // trim(s_ele%name) // '[' // trim(s_ele%slice_start) // ':' // trim(s_ele%slice_end) // ']')

end function find_slice_edge

end subroutine parser_expand_line

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine bp_set_ran_status
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine bp_set_ran_status

if (bp_com%extra%ran_seed == 0) then
  bp_com%extra%undeterministic_ran_function_called = .true.
endif

end subroutine bp_set_ran_status
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parser_debug_print_info (lat, debug_line, sequence)
!
! Subroutine to remove all null_ele elements.
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is not intended for general use.
!-

subroutine parser_debug_print_info (lat, debug_line, sequence)

type (lat_struct) lat
type (seq_struct), target, optional :: sequence(:)
type (seq_struct), pointer :: seq
integer i, j, ix
logical found
character(*) debug_line

!

found = .false.
call str_upcase (debug_line, debug_line)

if (index(debug_line, 'TIME') /= 0) then
  print *
  print *, '----------------------------------------'
  print '(a, f12.2)', 'Parse time (min):               ', (bp_com%time1 - bp_com%time0) / 60
  print '(a, f12.2)', 'Lattice Bookkeeping time (min): ', (bp_com%time2 - bp_com%time1) / 60
  print '(a, f12.2)', 'Make_mat6 calc time (min):      ', (bp_com%time3 - bp_com%time2) / 60
  found = .true.
endif

if (index(debug_line, 'SEQ') /= 0 .and. present(sequence)) then
  print *
  print *, '----------------------------------------'
  print *, 'Number of sequences:', size(sequence)
  print *, 'Use Line: ', trim(lat%use_name)
  do i = 1, size(sequence)
    seq => sequence(i)
    print *
    print *, '----------------------------------'
    print '(a, i3, 2x, a, 2x, a)', 'Sequence: ', i, trim(seq%name), trim(this_type(seq%type))
    print '(1x, a, 2x, a)', 'In file:', trim(seq%file_name)
    print *, 'Multipass:', seq%multipass
    do j = 1, size(seq%ele)
      print '(i5, 2x, a16, 2x, i1, 2x, a)', j, &
              this_type(seq%ele(j)%type), seq%ele(j)%ele_orientation, trim(seq%ele(j)%name)
    enddo
  enddo
  found = .true.
endif

if (index(debug_line, 'CONST') /= 0) then
  print *
  print *, '----------------------------------------'
  print *, 'Number of Defined Constants:', bp_com%i_const_tot - bp_com%i_const_init
  do i = bp_com%i_const_init+1, bp_com%i_const_tot
    print '(i6, 2x, a, es18.10)', i, bp_com2%const(i)%name, bp_com2%const(i)%value
  enddo
  found = .true.
endif

if (index(debug_line, 'SLAVE') /= 0) then
  print *
  print *, '----------------------------------------'
  print *, 'Number of Elements in Tracking Lattice:', lat%n_ele_track
  do i = 1, lat%n_ele_track
    print *, '-------------'
    print *, 'Ele #', i
    call type_ele (lat%ele(i), .false., 0, .false., 0, .true., .true., .false., all$, .true.)
  enddo
  found = .true.
endif

if (index(debug_line, 'LORD') /= 0) then
  print *
  print *, '----------------------------------------'
  print *, 'LORD elements: ', lat%n_ele_max - lat%n_ele_track
  do i = lat%n_ele_track+1, lat%n_ele_max
    print *, '-------------'
    print *, 'Ele #', i
    call type_ele (lat%ele(i), .false., 0, .false., 0, .true., .true., .false., all$, .true.)
  enddo
  found = .true.
endif

if (index(debug_line, 'LATTICE') /= 0) then  
  print *
  print *, '----------------------------------------'
  print *, 'Lattice Used: ', lat%use_name
  print *, 'Number of lattice elements:', lat%n_ele_track
  print *, 'List:                                 Key                 Length         S'
  do i = 1, lat%n_ele_track
    print '(i4, 2a, 3x, a, 2f10.2)', i, ') ', lat%ele(i)%name(1:30),  &
      key_name(lat%ele(i)%key), lat%ele(i)%value(l$), lat%ele(i)%s
  enddo
  print *, '---- Lord Elements ----'
  do i = lat%n_ele_track+1, lat%n_ele_max
    print '(2x, i4, 2a, 3x, a, 2f10.2)', i, ') ', lat%ele(i)%name(1:30),  &
           key_name(lat%ele(i)%key), lat%ele(i)%value(l$), lat%ele(i)%s
  enddo
  found = .true.
endif

ix = index(debug_line, 'ELE')
if (ix /= 0) then
  print *
  print *, '----------------------------------------'
  call string_trim (debug_line(ix+3:), debug_line, ix)
  do
    if (ix == 0) exit
    read (debug_line, *) i
    print *
    print *, '----------------------------------------'
    print *, 'Element #', i
    call type_ele (lat%ele(i), .false., 0, .true., 0, .true., .true., .true., all$, .true.)
    call string_trim (debug_line(ix+1:), debug_line, ix)
  enddo
  found = .true.
endif

if (index(debug_line, 'PARTICLE_START') /= 0) then
  print *
  print *, '----------------------------------------'
  print *, 'particle_start:'
  print '(3x, 6es13.4)', lat%particle_start%vec      
  found = .true.
endif

!--------------------------------

if (.not. found) then
  print *, 'BAD PARSER_DEBUG LINE: ' // trim(debug_line)
endif

!--------------------------------------------------------------
contains

function this_type (ix) result (type_str)
integer ix, k
character(20) type_str

select case (ix)
case (line$)
  type_str = 'Line'
case (list$)
  type_str = 'List'
case (replacement_line$)
  type_str = 'Replacement_Line'
case (element$)
  type_str = 'Element'
case (0)
  type_str = 'Zero!!'
case default
  type_str = '???'
end select

end function

end subroutine parser_debug_print_info

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parse_cartesian_map (ct_map, ele, lat, delim, delim_found, err_flag)
!
! Subroutine to parse a "cartesian_map = {}" construct
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is private to bmad_parser_mod.
! This must read in:
! {type = ,
!    dr = , 
!    r0 = , 
!    pt(i,j,k) = ( (ex_re, ex_im), .... (bz_re, bz_im) ) 
!    .
!    .
!    . ) },
!-

subroutine parse_cartesian_map (ct_map, ele, lat, delim, delim_found, err_flag)

implicit none

type (cartesian_map_struct) :: ct_map
type (ele_struct), target :: ele
type (ele_struct), pointer :: match_ele
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (cartesian_map_term1_struct), allocatable :: term(:)
type (cartesian_map_term1_struct), pointer :: tm

real(rp) kx, ky, kz, tol

complex(rp), pointer :: c_ptr(:)

integer n, ix_word, i_term, ib, ie, im, ix

character(80) err_str
character(40) word, word2, name, attrib_name
character(1) delim, delim2

logical err_flag, delim_found

!

name = 'xxx'
err_flag = .true.
ct_map%ptr%file = bp_com%line2_file_name   ! In case there are no terms

if (ele%key == wiggler$) ct_map%master_parameter = polarity$

!

do

  ! Read attriubute
  call get_next_word (attrib_name, ix_word, '{}=,()', delim, delim_found, call_check = .true.)
  if (.not. expect_this ('=', .true., .false., 'IN CARTESIAN_MAP DEFINITION', ele, delim, delim_found)) return

  select case (attrib_name)

  case ('FIELD_SCALE')
    call parse_evaluate_value (ele%name, ct_map%field_scale, lat, delim, delim_found, err_flag, ',}', ele)

  case ('R0')
    if (.not. equal_sign_here(ele, delim)) return
    if (.not. parse_real_list (lat, trim(ele%name) // 'GRID_FIELD', ct_map%r0, .true., delim, delim_found)) return
    if (.not. expect_one_of (',}', .false., ele%name, delim, delim_found)) return


  case ('ELE_ANCHOR_PT', 'FIELD_TYPE', 'MASTER_PARAMETER')
    ! Expect "<component> = "
    if (delim /= '=') then
      call parser_error ('NO "=" SIGN FOUND AFTER ' // attrib_name,  &
                         'IN ELEMENT: ' // ele%name)
      return
    endif
    call get_next_word (word2, ix_word, ',}', delim, delim_found)

    !

    select case (attrib_name)
 
    case ('MASTER_PARAMETER')
      if (word2 == 'NONE') then
        ix = 0
      else
        ix = attribute_index(ele, word2)
        if (ix < 1) then
          call parser_error ('BAD NAME FOR "MASTER_PARAMETER = <NAME>" CONSTRUCT', &
                             'FOUND IN ELEMENT: ' // ele%name)
          return
        endif
      endif
      ct_map%master_parameter = ix

    case ('ELE_ANCHOR_PT')
      call match_word(word2, anchor_pt_name(1:), ct_map%ele_anchor_pt, can_abbreviate = .false., matched_name = name)
  
    case ('FIELD_TYPE')
      call match_word(word2, em_field_type_name(1:2), ct_map%field_type, can_abbreviate = .false., matched_name = name)
  
    end select

    !

    if (name == '') then
      call parser_error ('UNKNKOWN ' // trim(attrib_name) // ' VALUE:' // word2, &
                         'IN ELEMENT: ' // ele%name)
      return        
    endif      
      
  case ('TERM')

    if (.not. allocated(ct_map%ptr%term)) then
      allocate (ct_map%ptr%term(1))
      tm => ct_map%ptr%term(1)
      ! Set %file to be the last called file:<line_number>. 
      ct_map%ptr%file = bp_com%line2_file_name
    else
      n = size(ct_map%ptr%term) + 1
      call move_alloc(ct_map%ptr%term, term)
      allocate (ct_map%ptr%term(n))
      ct_map%ptr%term(1:n-1) = term
      deallocate (term)
      tm => ct_map%ptr%term(n)
    endif

    err_str = trim(ele%name) // ' CARTESIAN_MAP TERM'

    if (.not. expect_this ('{', .false., .false., 'AFTER "TERM =" IN CARTESIAN_MAP DEFINITION', ele, delim, delim_found)) return
    call parse_evaluate_value (err_str, tm%coef, lat, delim, delim_found, err_flag, ',', ele);  if (err_flag) return
    call parse_evaluate_value (err_str, tm%kx, lat, delim, delim_found, err_flag, ',', ele);    if (err_flag) return
    call parse_evaluate_value (err_str, tm%ky, lat, delim, delim_found, err_flag, ',', ele);    if (err_flag) return
    call parse_evaluate_value (err_str, tm%kz, lat, delim, delim_found, err_flag, ',', ele);    if (err_flag) return
    call parse_evaluate_value (err_str, tm%x0, lat, delim, delim_found, err_flag, ',', ele);    if (err_flag) return
    call parse_evaluate_value (err_str, tm%y0, lat, delim, delim_found, err_flag, ',', ele);    if (err_flag) return
    call parse_evaluate_value (err_str, tm%phi_z, lat, delim, delim_found, err_flag, ',', ele); if (err_flag) return
    call get_switch ('FAMILY', ['Y ', 'X ', 'QU', 'SQ'], tm%family, err_flag, ele, delim, delim_found); if (err_flag) return
    if (.not. expect_this ('}', .true., .false., 'AFTER "FAMILY" SWITCH', ele, delim, delim_found)) return
    if (.not. expect_one_of(',}', .false., ele%name, delim, delim_found)) return

    kx = tm%kx
    ky = tm%ky
    kz = tm%kz
    tol = 1d-5 * (kx**2 + ky**2 + kz**2)

    if (abs(ky**2 - kx**2 - kz**2) < tol) then
      tm%form = hyper_y$
      ky = sign_of(ky, .false.) * sqrt(kx**2 + kz**2)

    elseif (abs(ky**2 + kx**2 - kz**2) < tol) then
      tm%form = hyper_xy$
      kz = sign_of(kz, .false.) * sqrt(kx**2 + ky**2)

    elseif (abs(ky**2 - kx**2 + kz**2) < tol) then
      tm%form = hyper_x$
      kx = sign_of(kx, .false.) * sqrt(ky**2 + kz**2)

    else
      call parser_error ('CARTESIAN_MAP TERM DOES NOT HAVE CONSISTANT Kx, Ky, and Kz', &
                    'FOR ELEMENT: ' // ele%name)
      err_flag = .true.
      return
    endif

  case default
    if (attrib_name == '') then
      call parser_error ('MANGLED CARTESIAN_MAP DEFINITION FOR ELEMENT: ' // ele%name)
    else
      call parser_error ('UNKNOWN CARTESIAN_MAP COMPONENT: ' // attrib_name, 'FOR ELEMENT: ' // ele%name)
    endif
    return

  end select

  ! Possible "}" is end of mode
  if (delim == '}') exit

enddo

!

if (.not. expect_one_of (', ', .false., ele%name, delim, delim_found)) return
err_flag = .false.

! Check if data has already been read in for another element.
! If so, save space by pointing to the data.

call find_matching_fieldmap(ct_map%ptr%file, ele, cartesian_map$, match_ele, im)
if (im > 0) then
  deallocate(ct_map%ptr)
  ct_map%ptr => match_ele%cartesian_map(im)%ptr
  ct_map%ptr%n_link = ct_map%ptr%n_link + 1        
endif

end subroutine parse_cartesian_map

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! parse_cylindrical_map (cl_map, ele, lat, delim, delim_found, err_flag)
!
! Subroutine to parse a "cylindrical_map = {}" construct
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is private to bmad_parser_mod.
! This must read in:
! {type = ,
!    dr = , 
!    r0 = , 
!    pt(i,j,k) = ( (ex_re, ex_im), .... (bz_re, bz_im) ) 
!    .
!    .
!    . ) },
!-

subroutine parse_cylindrical_map (cl_map, ele, lat, delim, delim_found, err_flag)

implicit none

type (cylindrical_map_struct), pointer :: cl_map
type (ele_struct), target :: ele
type (ele_struct), pointer :: match_ele
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch

real(rp), allocatable :: array(:)

complex(rp), pointer :: c_ptr(:)

integer ix_word, i_term, ib, ie, im, ix

character(1) delim, delim2
character(40) word, word2, name, attrib_name

logical err_flag, delim_found, file_name_set

! Init

err_flag = .true.
allocate (array(1024))
file_name_set = .false.
cl_map%ptr%file = bp_com%line2_file_name  ! In case there are no terms

!

do

  ! Read attriubute
  call get_next_word (attrib_name, ix_word, '{}=,()', delim, delim_found, call_check = .true.)

  select case (attrib_name)

  case ('DPHI0_REF')
    call parser_error ('THE ATTRIBUTE NAME "DPHI0_REF" HAS BEEN CHANGED TO "PHI0_FIELDMAP"', &
                       'PLEASE MAKE THE CHANGE IN THE LATTICE FILE.')

  case ('PHI0_FIELDMAP')
    call parse_evaluate_value (ele%name, cl_map%phi0_fieldmap, lat, delim, delim_found, err_flag, ',}', ele)

  case ('THETA0_AZIMUTH')
    call parse_evaluate_value (ele%name, cl_map%theta0_azimuth, lat, delim, delim_found, err_flag, ',}', ele)

  case ('FIELD_SCALE')
    call parse_evaluate_value (ele%name, cl_map%field_scale, lat, delim, delim_found, err_flag, ',}', ele)

  case ('DZ')            
    call parse_evaluate_value (ele%name, cl_map%dz, lat, delim, delim_found, err_flag, ',}', ele)

  case ('R0')
    if (.not. equal_sign_here(ele, delim)) return
    if (.not. parse_real_list (lat, trim(ele%name) // ' GRID_FIELD', cl_map%r0, .true., delim, delim_found)) return
    if (.not. expect_one_of (',}', .false., ele%name, delim, delim_found)) return

  case ('M')
   call parser_get_integer (cl_map%m, word, ix_word, delim, delim_found, err_flag, 'BAD CYLINDRICAL_MAP M CONSTRUCT', 'IN ELEMENT: ' // ele%name)

  case ('HARMONIC')
    call parser_get_integer (cl_map%harmonic, word, ix_word, delim, delim_found, err_flag, 'BAD CYLINDRICAL_MAP HARMONIC CONSTRUCT', 'IN ELEMENT: ' // ele%name)

  case ('MASTER_PARAMETER')
    call get_next_word (word, ix_word, ',}', delim, delim_found)
    if (word == 'NONE') then
      ix = 0
    else
      ix = attribute_index(ele, word)
      if (ix < 1) then
        call parser_error ('BAD NAME FOR "MASTER_PARAMETER = <NAME>" CONSTRUCT', &
                             'FOUND IN ELEMENT: ' // ele%name)
        return
      endif
    endif
    cl_map%master_parameter = ix


  case ('ELE_ANCHOR_PT')
    ! Expect "<component> = "
    if (delim /= '=') then
      call parser_error ('NO "=" SIGN FOUND AFTER ELE_ANCHOR_PT ' // attrib_name,  &
                         'IN ELEMENT: ' // ele%name)
      return
    endif
    call get_next_word (word2, ix_word, ',}', delim, delim_found)

    ! Evaluate string into integer.

    call match_word(word2, anchor_pt_name(1:), cl_map%ele_anchor_pt, can_abbreviate = .false., matched_name = name)
  
    if (name == '') then
      call parser_error ('UNKNKOWN ELE_ANCHOR_PT ' // trim(word) // ': ' // word2, &
                         'FOUND IN ELEMENT: ' // ele%name)
      return        
    endif      
      
  case ('E_COEF_RE', 'E_COEF_IM', 'B_COEF_RE', 'B_COEF_IM')

    if (.not. file_name_set) then
      ! Set %file to be the last called file:<line_number>. 
      cl_map%ptr%file = bp_com%line2_file_name
      file_name_set = .true.
    endif

    ! Expect "("
    call get_next_word (word, ix_word, ',({', delim, delim_found)
    if (word /= '' .or. delim /= '(') then
      call parser_error ('NO "(" FOUND AFTER "' // trim(attrib_name) // ' =" ', &
                           'IN ELEMENT: ' // ele%name)
      return
    endif

    ! Read list of values.
    call re_allocate(array, 1024, .false.)
    do i_term = 1, 100000
      call get_next_word (word, ix_word, '{},()', delim, delim_found)
      if ((delim /= ',' .and. delim /= ')') .or. .not. is_real(word)) then
        call parser_error ('ERROR PARSING CYLINDRICAL_MAP ARRAY: ' // word2, &
                             'IN ELEMENT: ' // ele%name)
        return
      endif
      if (i_term > size(array)) call re_allocate(array, 2*size(array))
      read (word, *) array(i_term)
      if (delim == ')') exit
    enddo

    if (allocated(cl_map%ptr%term)) then
      if (size(cl_map%ptr%term) /= i_term) then
        call parser_error ('ARRAY SIZE MISMATCH FOR: ' // word2, &
                           'IN CYLINDRICAL_MAP DEFINITION IN ELEMENT: ' // ele%name)
        return
      endif
    else
      allocate(cl_map%ptr%term(i_term))
    endif

    select case (attrib_name)
    case ('E_COEF_RE', 'E_COEF_IM'); c_ptr => cl_map%ptr%term%e_coef 
    case ('B_COEF_RE', 'B_COEF_IM'); c_ptr => cl_map%ptr%term%b_coef
    end select

    if (attrib_name(8:9) == 'RE') then
      if (any(real(c_ptr) /= 0)) then
        call parser_error ('DUPLICATE ARRAY FOR: ' // attrib_name, &
                           'IN CYLINDRICAL_MAP IN ELEMENT: ' // ele%name)
        return
      endif
      c_ptr = c_ptr + array(1:i_term)

    else
      if (any(aimag(c_ptr) /= 0)) then
        call parser_error ('DUPLICATE ARRAY FOR: ' // attrib_name, &
                           'IN CYLINDRICAL_MAP IN ELEMENT: ' // ele%name)
        return
      endif
      c_ptr = c_ptr + i_imag * array(1:i_term)
    endif

    ! Expect "," or "}"
    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
    if (word /= '' .or. (delim /= ',' .and. delim /= '}')) then
      call parser_error ('BAD ' // trim(attrib_name) // ' = (...) CONSTRUCT', &
                           'FOUND IN CYLINDRICAL_MAP DEFINITION IN ELEMENT: ' // ele%name)
      return
    endif

  case default
    if (attrib_name == '') then
      call parser_error ('MANGLED CYLINDRICAL_MAP DEFINITION FOR ELEMENT: ' // ele%name)
    else
      call parser_error ('UNKNOWN CYLINDRICAL_MAP COMPONENT: ' // attrib_name, &
                         'FOR ELEMENT: ' // ele%name)
    endif
    return
    
  end select

  ! Possible "}" is end of mode
  if (delim == '}') exit

enddo

! Get final separator after grid construct.
 
if (.not. expect_one_of (', ', .false., ele%name, delim, delim_found)) return

deallocate(array)
err_flag = .false.

! Check if data has already been read in for another element.
! If so, save space by pointing to the data.

call find_matching_fieldmap(cl_map%ptr%file, ele, cylindrical_map$, match_ele, im)
if (im > 0) then
  deallocate(cl_map%ptr)
  cl_map%ptr => match_ele%cylindrical_map(im)%ptr
  cl_map%ptr%n_link = cl_map%ptr%n_link + 1        
endif

end subroutine parse_cylindrical_map

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! parse_grid_field (grid, ele, lat, delim, delim_found, err_flag)
!
! Subroutine to parse a "grid_field = {}" construct
!
! This subroutine is used by bmad_parser and bmad_parser2.
! This subroutine is private to bmad_parser_mod.
! This must read in:
! {type = ,
!    dr = , 
!    r0 = , 
!    pt(i,j,k) = ( (ex_re, ex_im), .... (bz_re, bz_im) ) 
!    .
!    .
!    . ) },
!-

subroutine parse_grid_field(g_field, ele, lat, delim, delim_found, err_flag)

type grid_pt_struct
  integer :: ix(3) = [1, 1, 1]
  complex(rp) :: field(6) = 0
end type

type (grid_field_struct), pointer :: g_field
type (ele_struct) :: ele
type (ele_struct), pointer :: match_ele
type (lat_struct),  target :: lat
type (branch_struct), pointer :: branch
type (grid_pt_struct), allocatable :: array(:), array2(:)

character(1) delim, delim2
character(40) :: word, word2, name
character(200) line

integer ix_word, ix_word2, ix
integer pt_counter, n, i, ib, ie, im, ix0, ix1, iy0, iy1, iz0, iz1
integer grid_dim,  num_dr, num_r0, ios

logical delim_found, delim_found2, err_flag, err_flag2

! Init. Last thing read in was initial "{"

allocate(array(1024))
pt_counter = 0
err_flag = .true.
g_field%ptr%file = bp_com%line2_file_name    ! In case there are no terms

do
  ! Read attriubute
  call get_next_word (word, ix_word, '{}=,()', delim, delim_found, call_check = .true.)
  if (word == '' .and. delim == '{') word = '{'

  select case (word)

  case ('PHI0_FIELDMAP')
    call parse_evaluate_value (ele%name, g_field%phi0_fieldmap, lat, delim, delim_found, err_flag, ',}', ele)

  case ('FIELD_SCALE')
    call parse_evaluate_value (ele%name, g_field%field_scale, lat, delim, delim_found, err_flag, ',}', ele)

  case ('HARMONIC')
    call parser_get_integer (g_field%harmonic, word, ix_word, delim, delim_found, err_flag, 'CANNOT READ GRID_FIELD HARMONIC NUMBER', 'IN ELEMENT: ' // ele%name)

  case ('INTERPOLATION_ORDER')
    call parser_get_integer (g_field%interpolation_order, word, ix_word, delim, delim_found, err_flag, 'CANNOT READ INTERPOLATION_ORDER NUMBER', 'IN ELEMENT: ' // ele%name)

  case ('MASTER_PARAMETER')
    call get_next_word (word, ix_word, ',}', delim, delim_found)
    if (word == 'NONE') then
      ix = 0
    else
      ix = attribute_index(ele, word)
      if (ix < 1) then
        call parser_error ('UNKNOWN ELEMENT PARAMETER NAME FOR GRID_FIELD MASTER_PARAMETER: ' // word, &
                           'THIS IS NOT A PARAMETER DEFINED FOR THIS TYPE OF ELEMENT: ' // key_name(ele%key), &
                           'FOUND IN ELEMENT: ' // ele%name)
      endif
    endif
    g_field%master_parameter = ix

  case ('FIELD_TYPE', 'ELE_ANCHOR_PT', 'GEOMETRY')
    if (.not. equal_sign_here(ele, delim)) return

    call get_next_word (word2, ix_word, ',}', delim, delim_found)
    ! Check to see if this is a valid type by checking against grid_field_geometry_name(:)

    if (word == 'FIELD_TYPE') then
      call match_word(word2, em_field_type_name, g_field%field_type, can_abbreviate = .false., matched_name = name)
    elseif (word == 'GEOMETRY') then
      call match_word(word2, grid_field_geometry_name(1:), g_field%geometry, can_abbreviate = .false., matched_name = name)
    else
      call match_word(word2, anchor_pt_name(1:), g_field%ele_anchor_pt, can_abbreviate = .false., matched_name = name)
    endif
  
    if (name == '') then
      call parser_error ('UNKNKOWN GRID_FIELD ' // trim(word) // ': ' // word2, &
                         'FOUND IN GRID_FIELD DEFINITION FOR ELEMENT: ' // ele%name)
      return        
    endif      
      
  case ('CURVED_COORDS', 'CURVED_REF_FRAME')  ! 'curved_coords' is old style.
    if (.not. equal_sign_here(ele, delim)) return
    call get_next_word (word2, ix_word, ':,=()', delim, delim_found, .true.)
    g_field%curved_ref_frame = evaluate_logical (word2, ios)
    if (ios /= 0 .or. ix_word == 0) then
      call parser_error ('BAD GRID_FIELD CURVED_REF_FRAME SETTING ' // word2, 'FOR: ' // ele%name)
    endif

  case ('R0')
    if (.not. equal_sign_here(ele, delim)) return
    if (.not. parse_real_list (lat, trim(ele%name) // ' GRID_FIELD', g_field%r0, .false., delim, delim_found)) return
    if (.not. expect_one_of (',}', .false., ele%name, delim, delim_found)) return

    case ('DR')
    if (.not. equal_sign_here(ele, delim)) return
    ! expect ( 1.) or (1. , 2.) or (1., 2., 3.)
    if (.not. parse_real_list (lat, trim(ele%name) // ' GRID', g_field%dr, .false., delim, delim_found)) return
    call get_next_word (word, ix_word, ',}', delim, delim_found)     
    if (word /= '') then
      call parser_error ('BAD INPUT AFTER DR DEFINITION: ' // word , &
                                 'FOUND IN GRID_FIELD DEFINITION FOR ELEMENT: ' // ele%name)
      return
    end if

  case ('PT')

    if (pt_counter == 0) then
      ! Set %file to be the last called file:<line_number>.
      g_field%ptr%file = bp_com%line2_file_name
    endif

    ! Increment 
    pt_counter = pt_counter + 1
    ! Reallocate temporary structure if needed
    n = size(array)
    if (pt_counter > n) then
      call move_alloc(array, array2)
      allocate(array(2*n))
      array(1:n) = array2
      deallocate(array2)
    end if

    ! Get indices
    bp_com%parse_line = delim // bp_com%parse_line
    if (.not. parse_integer_list (trim(ele%name) // ' GRID_FIELD PT', lat, array(pt_counter)%ix, .false., delim, delim_found)) return
      
    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
    call get_next_word (word2, ix_word2, '{}=,()', delim2, delim_found2)
    if ((word /= '') .or. (word2 /= '') .or. (delim /= '=') .or. (delim2 /= '(')) then
      call parser_error ('BAD GRID_FIELD PT CONSTRUCT, NO  = "(" ', &
                 'FOUND IN GRID_FIELD DEFINITION FOR ELEMENT: ' // ele%name)
      return
    end if
    ! Get as many field components as listed
    do i = 1, 6
      call parse_complex_component(array(pt_counter)%field(i), array(pt_counter), delim, err_flag2)
      if (err_flag2) return
      if (delim == ')') exit
      if (delim /= ',') then
        call parser_error ('BAD GRID_FIELD PT CONSTRUCT, NO "," BETWEEN FIELD COMPONENTS', &
            'FOUND IN GRID_FIELD DEFINITION FOR ELEMENT: ' // ele%name)
        return
      end if
    end do

    select case (i)
    case (3)
      if (g_field%field_type == mixed$) then
        call parser_error ('FIELD_GRID WITH FIELD_TYPE = MIXED IN ELEMENT: ' // ele%name, 'MUST SPECIFY BOTH E AND B FIELDS.')
        return
      endif        
    case (6)
      if (g_field%field_type /= mixed$) then
        call parser_error ('FIELD_GRID WITH FIELD_TYPE = ELECTRIC OR MAGNETIC IN ELEMENT: ' // ele%name, 'CANNOT SPECIFY BOTH E AND B FIELDS.')
        return
      endif
    case default
      call parser_error ('FIELD_GRID IN ELEMENT: ' // ele%name, 'DOES NOT HAVE THE CORRECT NUMBER OF FIELD COMPOENTS IN PT = (...) CONSTRUCT')
      return
    end select

    ! Expect , or }
    if (.not. expect_one_of(',}', .false., ele%name, delim, delim_found)) return

  case ('{')
    ! Set %file to be the last called file:<line_number>.
    g_field%ptr%file = bp_com%line2_file_name

    do pt_counter = 1, 10000000
      ! Reallocate temporary structure if needed
      n = size(array)
      if (pt_counter > n) then
        call move_alloc(array, array2)
        allocate(array(2*n))
        array(1:n) = array2
        deallocate(array2)
      end if

      ! Get indices

      if (g_field%geometry == rotationally_symmetric_rz$) then
        if (.not. parser_fast_integer_read(array(pt_counter)%ix(1:2), ele, ':', 'FIELD_GRID POINT TABLE')) return
      else
        if (.not. parser_fast_integer_read(array(pt_counter)%ix, ele, ':', 'FIELD_GRID POINT TABLE')) return
      endif

      ! Get as many field components as listed

      if (g_field%field_type == mixed$) then
        if (.not. parser_fast_complex_read(array(pt_counter)%field, ele, delim, 'FIELD_GRID POINT TABLE')) return
      else
        if (.not. parser_fast_complex_read(array(pt_counter)%field(1:3), ele, delim, 'FIELD_GRID POINT TABLE')) return
      endif

      ! Allow extra comma "...(0 0), }" at end of field point list.

      call string_trim (bp_com%parse_line, bp_com%parse_line, ix)
      if (delim == ',' .and. bp_com%parse_line(1:1) == '}') then
        delim = '}'
        bp_com%parse_line = bp_com%parse_line(2:)
      endif

      if (delim == '}') exit
    enddo

    ! Expect , or }
    if (.not. expect_one_of(',}', .false., ele%name, delim, delim_found)) return

  case default
    if (word == '') then
      call parser_error ('MANGLED GRID_FIELD DEFINITION FOR ELEMENT: ' // ele%name)
    else
      call parser_error ('UNKNOWN GRID_FIELD COMPONENT: ' // word, &
                         'FOR ELEMENT: ' // ele%name)
    endif
    return
    
  end select 


  ! Allow extra comma
  call string_trim (bp_com%parse_line, bp_com%parse_line, ix)
  if (delim == ',' .and. bp_com%parse_line(1:1) == '}') then
    delim = '}'
    bp_com%parse_line = bp_com%parse_line(2:)
  endif

  if (delim == '}') exit   

enddo

! Get final separator after grid_field construct.
 
if (.not. expect_one_of (', ', .false., ele%name, delim, delim_found)) return

! Clear pts

if (allocated(g_field%ptr%pt)) deallocate(g_field%ptr%pt)

! Allocate grid_field for different dimensions

grid_dim = grid_field_dimension(g_field%geometry)

if (grid_dim < 1 .or. grid_dim > 3) then
  call parser_error ('BAD GRID_FIELD DIMENSION', &
             'FOUND IN GRID_FIELD DEFINITION FOR ELEMENT: ' // ele%name)
  return
endif

ix0 = minval(array(1:pt_counter)%ix(1))
ix1 = maxval(array(1:pt_counter)%ix(1))
iy0 = minval(array(1:pt_counter)%ix(2))
iy1 = maxval(array(1:pt_counter)%ix(2))
iz0 = minval(array(1:pt_counter)%ix(3))
iz1 = maxval(array(1:pt_counter)%ix(3))

allocate(g_field%ptr%pt(ix0:ix1, iy0:iy1, iz0:iz1))

n = (ix1+1-ix0) * (iy1+1-iy0) * (iz1+1-iz0)
if (n /= pt_counter) then
  call out_io (s_warn$, bp_com%parser_name, &
                 'Note: Number of grid_field points (\i0\) in the file not equal to grid_field array size (\i0\:\i0\, \i0\:\i0\, \i0\:\i0\).', &
                 'for element: ' // ele%name, &
                 i_array = [pt_counter, ix0, ix1, iy0, iy1, iz0, iz1])
endif

! Assign grid_field values
do i = 1, pt_counter
  ix1 = array(i)%ix(1)
  iy1 = array(i)%ix(2)
  iz1 = array(i)%ix(3)
  if (g_field%field_type == magnetic$) then
    g_field%ptr%pt(ix1, iy1, iz1)%B = array(i)%field(1:3)
  else
    g_field%ptr%pt(ix1, iy1, iz1)%E = array(i)%field(1:3)
    g_field%ptr%pt(ix1, iy1, iz1)%B = array(i)%field(4:6)
  endif
end do

! Clear temporary array

deallocate(array)

! Check if grid_field data has already been read in for another element.
! If so, save space by pointing to the existing grid.

call find_matching_fieldmap(g_field%ptr%file, ele, grid_field$, match_ele, im)
if (im > 0) then
  deallocate(g_field%ptr)
  g_field%ptr => match_ele%grid_field(im)%ptr
  g_field%ptr%n_link = g_field%ptr%n_link + 1        
endif

err_flag = .false.

!-----------------------------------------------------------
contains 

! subroutine parse_complex_component(complex_component, delim, err_flag)
! looks for (x, y) or x followed by , or ) 
! returns complex field_component and next delim, which should be , or )

subroutine parse_complex_component(complex_component, array_pt, delim, err_flag)

type (grid_pt_struct) array_pt

real(rp) x, y
complex(rp) complex_component
integer ix_word
logical delim_found, err_flag

character(1) delim
character(40) word

!

err_flag = .true.

! Expect "(" for complex, "," for real in the middle of the list, and ")" at the end of the list

call string_trim(bp_com%parse_line, bp_com%parse_line, ix_word)
if (bp_com%parse_line(1:1) == '(') then
  bp_com%parse_line = bp_com%parse_line(2:)
  call get_this_value(x, array_pt, ',', delim, delim_found, err_flag); if (err_flag) return
  call get_this_value(y, array_pt, ')', delim, delim_found, err_flag); if (err_flag) return
  complex_component = cmplx(x, y, rp)
  call get_next_word (word, ix_word, ',)', delim, delim_found)

else
  call get_this_value(x, array_pt, ',)', delim, delim_found, err_flag); if (err_flag) return
  complex_component = cmplx(x, 0.0_rp, rp)
endif

err_flag = .false.

end subroutine parse_complex_component

!-----------------------------------------------------------
! contains 

subroutine get_this_value (val, array_pt, delim_list, delim, delim_found, err_flag)

type (grid_pt_struct) array_pt

real(rp) val
integer ix_word
logical delim_found, err_flag

character(*) delim_list, delim
character(40) :: word

!

err_flag = .false.
call get_next_word (word, ix_word, delim_list, delim, delim_found)
if (is_real(word, real_num = val)) return

call parser_error ('BAD FIELD VALUE IN GRID_FIELD PT: ' // word, &
                   'AT GRID POINT INDEX: ' // int_str(array_pt%ix(1)) // ', ' // int_str(array_pt%ix(2)) // &
                                                                         ', ' // int_str(array_pt%ix(3)), &
                   '[NOTE: EXPRESSIONS FOR GRID_FIELD PT FIELD VALUES NOT IMPLEMENTED.]')
err_flag = .true.

end subroutine get_this_value

end subroutine parse_grid_field

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine parse_gen_grad_map (gg_map, ele, lat, delim, delim_found, err_flag)
!
! Subroutine to parse a "gen_grad_map = {}" construct
!-

subroutine parse_gen_grad_map (gg_map, ele, lat, delim, delim_found, err_flag)

implicit none

type (gen_grad_map_struct), pointer :: gg_map
type (gen_grad1_struct), allocatable :: gg(:)
type (gen_grad1_struct), pointer :: gg1
type (ele_struct), target :: ele
type (ele_struct), pointer :: match_ele
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (em_taylor_term_struct), allocatable :: term(:)
type (em_taylor_term_struct), pointer :: tm

real(rp) coef, deriv(0:50), z(1)

integer i, j, nn, n, i_ib, ie, im, ix, ios, nder, n_gg
integer lb, i_out, ix_word, iz, iz_here

character(80) err_str
character(40) word, word2, name, attrib_name
character(1) delim, delim2

logical err_flag, delim_found, valid

! Init

name = 'xxx'
err_flag = .true.
gg_map%file = bp_com%line2_file_name

!

do

  ! Read attriubute
  call get_next_word (attrib_name, ix_word, '{}=,()', delim, delim_found, call_check = .true.)
  if (attrib_name == '' .and. delim == '{') attrib_name = '{'

  select case (attrib_name)

  case ('FIELD_SCALE')
    call parse_evaluate_value (ele%name, gg_map%field_scale, lat, delim, delim_found, err_flag, ',}', ele)

  case ('R0')
    if (.not. equal_sign_here(ele, delim)) return
    if (.not. parse_real_list (lat, trim(ele%name) // ' GEN_GRAD_MAP R0', gg_map%r0, .true., delim, delim_found)) return
    if (.not. expect_one_of (',}', .false., ele%name, delim, delim_found)) return

  case ('DZ')
    call parse_evaluate_value (ele%name, gg_map%dz, lat, delim, delim_found, err_flag, ',}', ele)

  case ('ELE_ANCHOR_PT', 'FIELD_TYPE', 'MASTER_PARAMETER')
    ! Expect "<component> = "
    if (delim /= '=') then
      call parser_error ('NO "=" SIGN FOUND AFTER ' // attrib_name,  &
                         'IN ELEMENT: ' // ele%name)
      return
    endif
    call get_next_word (word2, ix_word, ',}', delim, delim_found)

    !

    select case (attrib_name)
 
    case ('MASTER_PARAMETER')
      if (word2 == 'NONE') then
        ix = 0
      else
        ix = attribute_index(ele, word2)
        if (ix < 1) then
          call parser_error ('BAD NAME FOR "MASTER_PARAMETER = <NAME>" CONSTRUCT', &
                               'FOUND IN ELEMENT: ' // ele%name)
          return
        endif
      endif
      gg_map%master_parameter = ix

    case ('ELE_ANCHOR_PT')
      call match_word(word2, anchor_pt_name(1:), gg_map%ele_anchor_pt, can_abbreviate = .false., matched_name = name)
  
    case ('FIELD_TYPE')
      call match_word(word2, em_field_type_name(1:2), gg_map%field_type, can_abbreviate = .false., matched_name = name)
  
    end select

    !

    if (name == '') then
      call parser_error ('UNKNKOWN ' // trim(attrib_name) // ' VALUE:' // word2, &
                         'IN ELEMENT: ' // ele%name)
      return        
    endif      
      
  case ('CURVED_REF_FRAME')    ! 'curved_coords' is old style.
    if (.not. equal_sign_here(ele, delim)) return
    call get_next_word (word2, ix_word, ':,=()', delim, delim_found, .true.)
    gg_map%curved_ref_frame = evaluate_logical (word2, ios)
    if (ios /= 0 .or. ix_word == 0) then
      call parser_error ('BAD GEN_GRAD_MAP CURVED_REF_FRAME SETTING ' // word2, 'FOR: ' // ele%name)
    endif

  case ('CURVE')
    if (.not. expect_this('={', .true., .false., 'NO "={" AFTER "CURVE" IN GEN_GRAD_MAP', ele, delim, delim_found)) return   
    n_gg = size(gg_map%gg) + 1
    call move_alloc(gg_map%gg, gg)
    allocate (gg_map%gg(n_gg))
    gg_map%gg(1:n_gg-1) = gg
    deallocate(gg)
    gg1 => gg_map%gg(n_gg)

    do
      call get_next_word (attrib_name, ix_word, '{}=,()', delim, delim_found, call_check = .true.)
      select case (attrib_name)

      case ('M')
        call parser_get_integer (gg1%m, word, ix_word, delim, delim_found, err_flag, 'BAD GEN_GRAD%GG%M CONSTRUCT', 'IN ELEMENT: ' // ele%name)

      case ('KIND')
        call get_next_word (word, ix_word, ',}', delim, delim_found)
        call match_word(word, ['COS', 'SIN'], gg1%sincos, can_abbreviate = .false., matched_name = word)
        select case (word)
        case ('SIN');   gg1%sincos = sin$
        case ('COS');   gg1%sincos = cos$
        case default
          call parser_error ('BAD GEN_GRAD_MAP TYPE = <SIN-OR-COS>" CONSTRUCT', 'FOUND IN ELEMENT: ' // ele%name)
          return
        end select

      case ('DERIVS')
        iz = int_garbage$
        nder = -1

        if (.not. expect_this ('={', .true., .false., 'NO "={" AFTER "DERIVES" IN GEN_GRAD_MAP', ele, delim, delim_found)) return

        do
          if (.not. parser_fast_real_read(z, ele, ':', delim, 'GEN_GRAD_MAP DERIVS Z-POSITION')) return
          iz_here = nint(z(1)/gg_map%dz)

          if (iz == int_garbage$) then
            iz = iz_here

            if (n_gg == 1) then
              gg_map%iz0 = iz_here

            else
              if (gg_map%iz0 /= iz_here) then
                call parser_error ('LOWER BOUND INDEX IN GEN_GRAD_MAP DERIVS TABLE IS DIFFERENT FROM LOWER BOUND INDEX IN PRIOR DERIVS TABLE', &
                                   'FOR ELEMENT: ' // ele%name)
                return
              endif
            endif

          else
            iz = iz + 1
            if (iz /= iz_here) then
              call parser_error ('GEN_GRAD_MAP DERIVS TABLE INDEXES NOT IN CORRECT ORDER. EXPECTED: ' // int_str(iz) // ' BUT GOT: ' // int_str(iz_here), &
                                 'FOR ELEMENT: ' // ele%name)
              return
            endif

            if (gg_map%iz1 /= int_garbage$ .and. iz > gg_map%iz1) then
              call parser_error ('UPPER BOUND INDEX IN GEN_GRAD_MAP DERIVS TABLE IS GREATER THAN UPPER BOUND INDEX IN PRIOR DERIVS TABLE', &
                                 'FOR ELEMENT: ' // ele%name)
              return
            endif
          endif          

          if (nder == -1) then
            if (.not. parser_fast_real_read (deriv, ele, ',}', delim, 'GEN_GRAD_MAP DERIVS TABLE', .false., nder)) return
            nder = nder - 1   ! Since derivs are indexed from 0.
          else
            if (.not. parser_fast_real_read (deriv(0:nder), ele, ',}', delim, 'GEN_GRAD_MAP DERIVS TABLE')) return
          endif

          if (.not. allocated(gg1%deriv)) allocate (gg1%deriv(gg_map%iz0:gg_map%iz0+1000, 0:2*nder+1))
          if (iz > ubound(gg1%deriv,1)) call re_allocate2d(gg1%deriv, ubound(gg1%deriv,1)+1000, 2*nder+1, lb1 = gg_map%iz0, lb2 = 0)

          gg1%deriv(iz,0:nder) = deriv(0:nder)
          gg1%n_deriv_max = nder

          if (delim == '}') exit
        enddo

        if (gg_map%iz1 == int_garbage$) gg_map%iz1 = iz
        call re_allocate2d(gg1%deriv, gg_map%iz1, 2*nder+1, lb1 = gg_map%iz0, lb2 = 0)
        gg1%n_deriv_max = nder

        if (iz /= gg_map%iz1) then
          call parser_error ('ENDING IZ-INDEX IN GEN_GRAD_MAP DERIVS TABLE NOT IS DIFFERENT FROM PRIOR DERIVS TABLE', 'FOR ELEMENT: ' // ele%name)
          return
        endif

        if (.not. expect_one_of ('},', .false., ele%name, delim, delim_found)) return
      end select

      if (delim == '}') exit
      if (.not. expect_one_of (',}', .true., ele%name, delim, delim_found)) return
    enddo

    if (.not. expect_one_of (',} ', .false., ele%name, delim, delim_found)) return

  case default
    if (attrib_name == '') then
      call parser_error ('MANGLED GRID_FIELD DEFINITION FOR ELEMENT: ' // ele%name)
    else
      call parser_error ('UNKNOWN GRID_FIELD COMPONENT: ' // attrib_name, &
                         'FOR ELEMENT: ' // ele%name)
    endif
    return

  end select

  ! Possible "}" is end of mode
  if (delim == '}') exit

enddo

! Extend derivatives to form interpolating spline polynomial

do i = 1, size(gg_map%gg)
  gg1 => gg_map%gg(i)
  n = gg1%n_deriv_max

  do iz = gg_map%iz0, gg_map%iz1-1
    call n_spline_create(gg1%deriv(iz,0:n), gg1%deriv(iz+1,0:n), gg_map%dz, gg1%deriv(iz,:))
  enddo
enddo

! Get final separator after gen_grad_map construct.
 
if (.not. expect_one_of (', ', .false., ele%name, delim, delim_found)) return
err_flag = .false.

end subroutine parse_gen_grad_map

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function parse_integer_list (err_str, lat, int_array, exact_size, delim, delim_found, open_delim, 
!                                       separator, close_delim, default_value) result (is_ok)
!
! Routine to parse a list of integers of the form:
!    open_delim integer_1 separator integer_2 . . . close_delim
! Example:   "(1.2, 2.3, 4.4, 8.5)"
! 
! Similar to parse_integer_list2 except does not use allocatable array.
! See parse_integer_list2 for more details
!-

function parse_integer_list (err_str, lat, int_array, exact_size, delim, delim_found, open_delim, &
                                      separator, close_delim, default_value) result (is_ok)

implicit none

type (lat_struct) lat

integer int_array(:)
integer, optional :: default_value
integer, allocatable :: vec(:)

integer num_found

character(*) err_str, delim
character(*), optional :: open_delim, separator, close_delim

logical is_ok, exact_size, delim_found

!

is_ok = .false.
if (.not. parse_integer_list2 (err_str, lat, vec, num_found, delim, delim_found, size(int_array), &
                               open_delim, separator, close_delim, default_value)) return

if (num_found > size(int_array) .or. (exact_size .and. num_found < size(int_array))) then
  call parser_error (err_str)
  return
endif

int_array(1:num_found) = vec(1:num_found)

is_ok = .true.

end function parse_integer_list

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function parse_integer_list2 (err_str, lat, int_array, num_found, delim, delim_found, num_expected, 
!                                        open_delim, separator, close_delim, default_value) result (is_ok)
!
! Routine to parse a list of integers of the form
!    open_delim integer_1 separator integer_2 . . . close_delim
! Example:   (1, 2, 4, 8) 
!
! Input:
!   err_str       -- character(*): Error string to print if there is an error. 
!   lat           -- lat_struct: lattice
!   int_array(:)  -- Integer, allocatable: the array to be read in 
!
! Optional: 
!   num_expected = 1     -- integer: number of expected arguments. Used to initialize int_array.
!   open_delim   = '('   -- character(1): opening delimeter.
!   separator    = ','   -- character(1): separating character.
!   close_delim  = ')'   -- character(1): closing delimeter.
!   default_value = 0    -- real(rp): inital assignment of int_array elements.
!
! Output:
!   is_ok                   -- logical: Set True if everything is ok.
!   int_array(1:num_found)  -- integer(rp): Array of values.
!   num_found               -- integer: number of elements.
!   delim                   -- character(1): Delimiter found where the parsing of the input line stops.
!   delim_found             -- logical: Delimiter found? False if end of input command.
!-

function parse_integer_list2 (err_str, lat, int_array, num_found, delim, delim_found, num_expected, &
                                       open_delim, separator, close_delim, default_value) result (is_ok)


type (lat_struct) lat

integer, allocatable :: int_array(:)
integer :: num_found
integer, optional :: num_expected, default_value
character(*) err_str
character(*), optional :: open_delim, close_delim, separator
logical is_ok

! Local
integer num_expect
character(1) delim, op_delim, cl_delim, sep
character(40) :: word
real(rp) rval
integer  ix_word
logical delim_found, err_flag

! Optional arguments

is_ok = .false.
num_expect = integer_option (1, num_expected)
op_delim = '('
cl_delim = ')'
sep      = ','
if (present(open_delim)) op_delim = open_delim
if (present(close_delim)) cl_delim = close_delim
if (present(separator)) sep = separator

! Expect op_delim
if (op_delim /= '') then
  call get_next_word (word, ix_word, op_delim, delim, delim_found)
  if ((word /= '') .or. (delim /= op_delim)) then
    call parser_error (err_str)
    return
  end if
endif

! Initial allocation
call re_allocate(int_array, num_expected, .false.)
int_array = integer_option(0, default_value)

! counter
num_found = 0

! Get integers
do 

  call parse_evaluate_value ('BAD NUMBER IN: ' // err_str, rval, lat, delim, delim_found, err_flag, sep // cl_delim)
  if (err_flag) return
  if (abs(rval - nint(rval)) > 1d-10) then
    call parser_error ('BAD INTEGER NUMBER IN: ' // err_str)
    return
   end if    

  num_found = num_found + 1
  if (size(int_array) < num_found) then
    call re_allocate (int_array, 2*num_found, .false.)
    int_array(num_found:2*num_found) = integer_option(0, default_value)
  endif
  
  int_array(num_found) = nint(rval)
  
  ! Exit if cl_delim is found
  if (delim == cl_delim) exit
  
  ! Check separator
  if (delim /= sep) then
    call parser_error ('BAD SEPARATOR IN: ' // err_str)
    return  
  end if
  
end do

is_ok = .true.

end function parse_integer_list2

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function parse_real_list (lat, err_str, real_array, exact_size, delim, delim_found, open_delim, 
!                                separator, close_delim, default_value, num_found) result (is_ok)
!
! Routine to parse a list of reals of the form:
!    open_delim real_1 separator real_2 . . . close_delim
! Example:   "(1.2, 2.3, 4.4, 8.5)"
! 
! Similar to parse_real_list2 except does not use allocatable array.
! Also see: parse_real_lists.
!
! Input:
!   lat           -- lat_struct: Lattice
!   err_str       -- character(*): Error string to print if there is an error. 
!   exact_size    --
!   open_delim    --
!   separator     --
!   close_delim   --
!   default_value --
!
! Output:
!   real_array    --
!   delim         --
!   delim_found   --
!   num_found     --
!-

function parse_real_list (lat, err_str, real_array, exact_size, delim, delim_found, open_delim, &
                               separator, close_delim, default_value, num_found) result (is_ok)

implicit none

type (lat_struct) lat

real(rp) real_array(:)
real(rp), optional :: default_value
real(rp), allocatable :: vec(:)

integer, optional :: num_found
integer num_here

character(*) err_str, delim
character(*), optional :: open_delim, separator, close_delim

logical is_ok, exact_size, delim_found

!

is_ok = .false.
if (.not. parse_real_list2 (lat, err_str, vec, num_here, delim, delim_found, size(real_array), &
                          open_delim, separator, close_delim, default_value)) return

if (num_here > size(real_array) .or. (exact_size .and. num_here < size(real_array))) then
  call parser_error (err_str)
  return
endif

real_array(1:num_here) = vec(1:num_here)
if (present(num_found)) num_found = num_here

is_ok = .true.

end function parse_real_list

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function parse_real_lists (lat, err_str, table, size2, delim, delim_found) result (is_ok)
!
! Routine to parse a list of reals of the form:
!    {(re_11, re_12, ..., re_1<size2>), (re_21, re_22, ...), ...} 
! And re_IJ is put in table(I,J).
! size2 is the size of the inner array.
! The size of the outer array can be anything.
! 
! Similar to parse_real_list2 except does not use allocatable array.
! Also see: parse_real_lists.
!-

function parse_real_lists (lat, ele, err_str, table, size2, delim, delim_found) result (is_ok)

implicit none

type (lat_struct) lat
type (ele_struct) ele

real(rp), allocatable :: vec(:)
real(rp), allocatable :: table(:,:)

integer size2
integer nn, num_found

character(*) err_str
character(1) delim

logical is_ok, delim_found

!

is_ok = .false.
if (.not. allocated(table)) allocate (table(100,size2))

if (.not. expect_one_of ('{', .false., ele%name, delim, delim_found)) return
nn = 0
do
  if (.not. parse_real_list2 (lat, err_str, vec, num_found, delim, delim_found, size2, '(', ',', ')')) return
  if (num_found /= size2) then
    call parser_error (err_str)
    return
  endif
  nn = nn + 1
  if (nn > size(table, 1)) call re_allocate2d(table, 2*nn, size2)
  table(nn,:) = vec
  if (.not. expect_one_of (',}', .false., ele%name, delim, delim_found)) return
  if (delim == '}') exit
enddo

call re_allocate2d(table, nn, size2)

is_ok = .true.

end function parse_real_lists

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function parse_real_list2 (lat, err_str, real_array, num_found, delim, delim_found, num_expected, 
!                            open_delim, separator, close_delim, default_value) result (is_ok)
!
! Routine to parse a list of reals of the form:
!    open_delim real_1 separator real_2 . . . close_delim
! Example:   "(1.2, 2.3, 4.4, 8.5)"
!
! Also see:
!   pase_real_list
!   parse_real_lists.
!
! Input:
!  lat            -- lat_struct: lattice
!  err_str        -- character(*): Error string to print if there is an error. 
!  real_array     -- real(rp), allocatable: the array to be read in 
!
!   Optional: 
!   num_expected = 10    -- integer, optional: number of expected arguments
!                             Used to initialize real_array
!   open_delim   = '('   -- character(1), optional: opening delimeter
!   separator    = ','   -- character(1), optional: separating character
!   close_delim  = ')'   -- character(1), optional: closing delimeter
!   default_value = 0.0_rp -- real(rp) : inital assignment of real_array elements
!
! Output:
!   is_ok                   -- logical: Set True if everything is ok
!   real_array(1:num_found) -- real(rp) : Array of values
!   num_found               -- integer : number of elements
!   delim                   -- character(1): Delimiter found where the parsing of the input line stops.
!   delim_found             -- logical: Delimiter found? False if end of input command.
!-

function parse_real_list2 (lat, err_str, real_array, num_found, delim, delim_found, num_expected, &
          open_delim, separator, close_delim, default_value) result (is_ok)

! Arguments

type (lat_struct) lat

real(rp), allocatable :: real_array(:)
real(rp), optional :: default_value

integer :: num_found
integer, optional :: num_expected

logical is_ok

character(*) err_str
character(1), optional :: open_delim, close_delim, separator

! Local

real(rp) :: default_val, value

integer num_expect
integer  ix_word

character(1) delim, op_delim, cl_delim, sep
character(40) :: word

logical delim_found, err_flag

! Optional arguments

is_ok = .false.
num_expect = integer_option(10, num_expected)
default_val = real_option(0.0_rp, default_value)

op_delim = '('
cl_delim = ')'
sep      = ','
if (present(open_delim)) op_delim = open_delim
if (present(close_delim)) cl_delim = close_delim
if (present(separator)) sep = separator

! Expect op_delim
call get_next_word (word, ix_word, op_delim, delim, delim_found)
if ((word /= '') .or. (delim /= op_delim)) then
  call parser_error ('BAD OPENING DELIMITER IN: ' // err_str)
  return
end if

! Initial allocation
call re_allocate(real_array, num_expect, .false.)
real_array = default_val

! Get reals

num_found = 0

do 
  call parse_evaluate_value ('BAD REAL NUMBER IN: ' // err_str, value, lat, delim, delim_found, err_flag, sep // cl_delim)
  if (err_flag) return
  ! real is found
  num_found = num_found + 1
  ! reallocate if needed  
  if (size(real_array) < num_found) then
    call re_allocate (real_array, 2*num_found, .false.)
    real_array(num_found:2*num_found) = default_val
  endif

  ! Set value
   real_array(num_found) = value
  
  ! Exit if cl_delim is found
  if (delim == cl_delim) exit
  
  ! Check separator
  if (delim /= sep) then
    call parser_error ('BAD SEPARATOR IN: ' // err_str)
    return  
  end if
end do

is_ok = .true.
call re_allocate(real_array, num_found)

end function parse_real_list2

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

subroutine parser_get_integer (int_val, word, ix_word, delim, delim_found, err, str1, str2)

integer int_val, ix_word
logical err
character(*) word, delim
character(*), optional :: str1, str2
logical delim_found

!

call get_next_word (word, ix_word, ':,=(){} ', delim, delim_found, .true.)
if (.not. is_integer(word) ) then
  if (present(str1)) then
    call parser_error (str1, str2)
  else
    call parser_error ('INTEGER EXPECTED, I DO NOT UNDERSTAND: ' // word)
  endif
  err = .true.

else
  read (word, *) int_val
  err = .false.
endif

end subroutine parser_get_integer

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

subroutine parser_get_logical (attrib_name, this_logic, ele_name, delim, delim_found, err)

type (ele_struct) ele
integer ix_word, ios
logical delim_found, err, this_logic
character(*) attrib_name, ele_name, delim
character(40) word

!

call get_next_word (word, ix_word, ':,=()', delim, delim_found, .true.)
this_logic = evaluate_logical (word, ios)
if (ios /= 0 .or. ix_word == 0) then
  call parser_error ('BAD "' // trim(attrib_name) // '" SWITCH FOR: ' // ele_name, 'I DO NOT UNDERSTAND: ' // word)
  err = .true.
else
  err = .false.
endif

end subroutine parser_get_logical

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!+
! Function expect_this (expecting, check_delim, call_check, err_str, ele, delim, delim_found) result (is_ok)
!
! Checks that the next character or characters in the parse stream corresponds to the 
! characters in the expecting argument. For example, if expecting is ')={' these three characters 
! should be the next non-blank characters in the parse stream.
!
! Also see: expect_one_of
!
! Input:
!   expecting(*)  -- character: list of characters that are expected to be next in the parse stream.
!   check_delim   -- logical: If True then use delim argument as first token to check.
!                      A blank character indicates end of command is expected.
!   call_check    -- Logical: If True then check for 'call::<filename>' construct.
!   err_str       -- character(*): String used for error messages.
!   ele           -- ele_struct: Element parameters being parsed.
!
! Output:
!   delim         -- character(*): Final delim
!   delim_found   -- logical: Is there a final delim (as opposed to end of command).
!-

function expect_this (expecting, check_delim, call_check, err_str, ele, delim, delim_found) result (is_ok)

implicit none

type (ele_struct) ele
character(*) expecting, err_str
character(1) delim
character(40) word
logical is_ok, delim_found, check_delim, call_check
integer ix, ix_word

!

is_ok = .false.

do ix = 1, len(expecting)
  if (ix == 1 .and. check_delim) then
    word = ''
  else
    call get_next_word (word, ix_word, expecting(ix:ix), delim, delim_found, call_check = call_check)
  endif

  if (expecting(ix:ix) == ' ') then
    if (delim_found .or. word /= '') then
      call parser_error ('EXTRA STUFF ON LINE.', err_str)
      return
    endif
  elseif (delim /= expecting(ix:ix) .or. word /= '') then
    call parser_error ('NO "' // expecting // '" FOUND ' // err_str, 'FOR ELEMENT: ' // ele%name)
    return
  endif
enddo

is_ok = .true.

end function expect_this

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

subroutine get_switch (name, name_list, switch, err, ele, delim, delim_found)

type (ele_struct) :: ele
type (all_pointer_struct), allocatable :: ptr(:)

character(*) name, name_list(:)
character(1) delim
character(60) word
character(40) ele_name, attrib_name 
character(:), allocatable :: line
integer i, this_switch, switch, ix_word, ixp, ixp2
logical err, delim_found

!

err = .true.
call get_next_word (word, ix_word, ':,=(){}', delim, delim_found, .true.)

!

if (name == 'FIELD_CALC') then
  if (word == 'GRID' .or. word == 'MAP') then
    
    call parser_error ('FIELD_CALC = ' // word, 'HAS BEEN CHANGED TO FIELD_CALC = FIELDMAP', &
                       'Program will execute as normal...', &
                       '[But eventually this warning will be converted to an error. You have been warned!]', level = s_warn$)
    word = 'FIELDMAP'
  endif
endif

! If word is something like "q1[tracking_method]" then need retrieve this value.

ixp = index(word, '[')
if (ixp == 0) then
  call match_word (word, name_list, this_switch, can_abbreviate = .false.)
  if (this_switch < 1) then
    line = trim(name_list(1))
    do  i = 2, size(name_list)
      if (name_list(i) == null_name$) cycle
      if (upcase(name_list(i)) == 'GARBAGE!') cycle
      line = line // ', ' // trim(name_list(i))
    enddo
    call parser_error ('BAD "' // trim(name) // '" SWITCH FOR: ' // ele%name, 'I DO NOT UNDERSTAND: ' // word, &
                       'POSSIBILITIES ARE: ' // line)
    return
  else
    switch = this_switch
  endif

else
  ixp2 = len_trim(word)
  if (word(ixp2:ixp2) /= ']' .or. word(ixp+1:ixp2-1) /= name) then
    call parser_error ('BAD "' // trim(name) // '" SWITCH FOR: ' // ele%name, 'I DO NOT UNDERSTAND: ' // word)
    return
  endif

  ele_name = word(:ixp-1)
  attrib_name = word(ixp+1:ixp2-1)
  call pointers_to_attribute (ele%branch%lat, ele_name, attrib_name, .true., ptr, err, .false.)
  if (size(ptr) == 0) then
    call parser_error ('NO ELEMENT FOUND TO EVALUATE: ' // word, 'EVALUATING SWITCH IN ELEMENT: ' // ele%name)
    return
  endif
  if (size(ptr) > 1) then
    call parser_error ('MULTIPLE ELEMENTS FOUND FOR EVALUATING: ' // word, 'EVALUATING SWITCH IN ELEMENT: ' // ele%name)
    return
  endif
  if (associated(ptr(1)%i)) switch = ptr(1)%i
  if (associated(ptr(1)%r)) switch = nint(ptr(1)%r)
endif

err = .false.

end subroutine get_switch

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!+
! Function expect_one_of (delim_list, check_input_delim, ele_name, delim, delim_found) result (is_ok)
!
! Routine to check either that the current delimitor or the next character in the parse stream is the 
! expected delimitor.
! This routine is used for Bmad lattice file parsing and is not meant for general use.
!
! Also see: expect_this
!
! Input:
!   delim_list  -- character(*): List of expected (valid) delimitors. If list contains a space character
!                    then no delimitor (indicating the end of the command) is a valid possibility.
!   check_input_delim 
!               -- logical: If True, then check if delim argument is in the delim_list. 
!                    If False, check that the next character in the parse stream is an expected delimitor.
!   ele_name    -- character(*): Lattice element under construction. Used for error messages.
!   delim       -- character(1): Current delimitor that will be checked if check_input_delim = .true.
!
! Output:
!   delim       -- character(1): Next delim if check_input_delim = False.
!-

function expect_one_of (delim_list, check_input_delim, ele_name, delim, delim_found) result (is_ok)

type (ele_struct) ele
integer ix_word
character(*) delim_list, ele_name
character(1) delim
character(40) word
logical check_input_delim, delim_found, is_ok, must_have_delim

!

is_ok = .false.
must_have_delim = (index(delim_list, ' ') == 0)

if (check_input_delim) then
  if ((must_have_delim .and. .not. delim_found) .or. &
                        (delim /= '' .and. index(delim_list, delim) == 0)) then
    if (ele_name(1:1) == '!') then  ! Indicates is not an element
      call parser_error ('BAD DELIMITOR', 'FOR: ' // ele_name)
    else
      call parser_error ('BAD DELIMITOR', 'FOR ELEMENT: ' // ele_name)
    endif
    return
  endif

else
  call get_next_word (word, ix_word, '{}=,()[]', delim, delim_found)
  if (word /= '' .or. (must_have_delim .and. .not. delim_found) .or. &
                      (delim /= '' .and. index(delim_list, delim) == 0)) then
    if (ele_name(1:1) == '!') then  ! Indicates is not an element
      call parser_error ('BAD DELIMITOR', 'FOR: ' // ele_name)
    else
      call parser_error ('BAD DELIMITOR', 'FOR ELEMENT: ' // ele_name)
    endif
    return
  endif
endif

is_ok = .true.

end function expect_one_of

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

function equal_sign_here(ele, delim) result (is_here)

type (ele_struct) ele
character(40) word
character(1) delim
logical is_here

! Expect "<component> = "

is_here = .false.

if (delim /= '=') then
call parser_error ('NO "=" SIGN FOUND AFTER: ' // word,  &
                   'IN GRID_FIELD STRUCTURE IN ELEMENT: ' // ele%name)
  return
endif

is_here = .true.

end function equal_sign_here

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!+
! Subroutine parser_print_line(end_of_file)
!
! This routine is called when a print statement is found in the lattice file.
!-

subroutine parser_print_line(lat, end_of_file)

implicit none

type (lat_struct) lat
integer ix, ix2, n
real(rp) value
logical end_of_file, err_flag, delim_found
character(1) delim
character(280) print_line
character(*), parameter :: r_name = 'parser_print_line'

!

call string_trim (bp_com%input_line2, print_line, ix) ! To strip off initial "print"
print_line = print_line(ix+2:)

do
  ix = index(print_line, '`')
  if (ix == 0) exit
  ix2 = index(print_line(ix+1:), '`')
  if (ix2 == 0) exit
  ix2 = ix + ix2
  call parse_evaluate_value ('PRINT STATEMENT EVALUATION', value, lat, delim, delim_found, err_flag, string_in = print_line(ix+1:ix2-1))
  if (err_flag) then
    print_line = print_line(:ix-1) // '???' // print_line(ix2+1:)
  else
    print_line = print_line(:ix-1) // trim(real_to_string(value, 20, 14)) // print_line(ix2+1:)
  endif
enddo

if (.not. allocated(lat%print_str)) then
  allocate (lat%print_str(1))
else
  call re_allocate(lat%print_str, size(lat%print_str)+1)
endif
n = size(lat%print_str)
lat%print_str(n) = trim(print_line) ! To save in digested file
call out_io (s_important$, r_name, 'Message in Lattice File: ' // print_line, insert_tag_line = .false.)

end subroutine parser_print_line

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!+
! Subroutine parser_init_custom_elements (lat)
!
!-

subroutine parser_init_custom_elements (lat)

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
integer i, n
logical err

! Init custom stuff.

do n = 0, ubound(lat%branch, 1)
  branch => lat%branch(n)
  do i = 1, branch%n_ele_max
    ele => branch%ele(i)
    if (ele%key == custom$ .or. ele%tracking_method == custom$ .or. &
        ele%mat6_calc_method == custom$ .or. ele%field_calc == custom$ .or. &
        ele%aperture_type == custom_aperture$) then
      call init_custom (ele, err)
      if (err) bp_com%error_flag = .true.
    endif
  enddo
enddo

end subroutine parser_init_custom_elements

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!+
! Subroutine reallocate_sequence (sequence, n_seq)
!-

subroutine reallocate_sequence (sequence, n_seq)

implicit none

type (seq_struct), target, allocatable :: sequence(:), temp_seq(:)
integer n_seq, n

!

n = size(sequence)
call move_alloc (sequence, temp_seq)
allocate(sequence(n_seq))
sequence(1:n) = temp_seq

end subroutine reallocate_sequence

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!+
! Subroutine parse_superimpose_command(lat, delim)
!-

subroutine parse_superimpose_command(lat, ele, pele, delim)

type (lat_struct) lat
type (ele_struct) ele
type (parser_ele_struct) pele

integer ix_word, n
logical delim_found, err_flag

character(*) delim
character(40) var, value


!

err_flag = .false.

ele%lord_status = super_lord$
pele%superposition_command_here = .true.

if (delim /= ',') then
  call parser_error ('MISSING COMMA IN "SUPERIMPOSE" STATEMENT.')
  return
endif

do
  call get_next_word(var, ix_word, '=,', delim, delim_found, .true.)
  if (.not. expect_one_of('=', .true., '!SUPERIMPOSE STATEMENT', delim, delim_found)) return

  select case (var)
  case ('CREATE_JUMBO_SLAVE')
    call parser_get_logical (value, pele%create_jumbo_slave, ele%name, delim, delim_found, err_flag)
  case ('REF')
    call get_next_word(pele%ref_name, ix_word,  '=,', delim, delim_found, .true.)
  case ('ELEMENT')
    call get_next_word(pele%ele_name, ix_word,  '=,', delim, delim_found, .true.)
  case ('ELE_ORIGIN')
    call get_switch (value, anchor_pt_name(1:), pele%ele_pt, err_flag, ele, delim, delim_found)
  case ('REF_ORIGIN')
    call get_switch (value, anchor_pt_name(1:), pele%ref_pt, err_flag, ele, delim, delim_found)
  case ('OFFSET')
    call parse_evaluate_value ('SUPERIMPOSE STATEMENT OFFSET', pele%offset, lat, delim, delim_found, err_flag)
  case default
    call parser_error ('UNKNOWN PARAMETER OF SUPERIMPOSE COMMAND: ' // var)
    return
  end select

  if (err_flag) return
  if (.not. delim_found) exit
enddo

end subroutine parse_superimpose_command

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine init_surface_segment (phot, ix, iy)
!
! Routine to init the componentes in ele%photon%grid%pt(ix,iy) for use with segmented surface calculations.
!
! Input:
!   phot    -- photon_element_struct: Surface structure.
!   ix, iy  -- integer: index of grid point to init.
!-

subroutine init_surface_segment (phot, ix, iy)

type (photon_element_struct), target :: phot
type (surface_grid_pt_struct), pointer :: pt

real(rp) zt, x0, y0, dx, dy, coef_xx, coef_xy, coef_yy, coef_diag, g(3), gs
integer ix, iy

!

pt => phot%grid%pt(ix, iy)

x0 = ix * phot%grid%dr(1) + phot%grid%r0(1)
y0 = iy * phot%grid%dr(2) + phot%grid%r0(2)

pt%x0 = x0
pt%y0 = y0
pt%z0 = 0

pt%dz_dx = 0
pt%dz_dy = 0
coef_xx = 0; coef_xy = 0; coef_yy = 0

do ix = 0, ubound(phot%curvature%xy, 1)
do iy = 0, ubound(phot%curvature%xy, 2) - ix
  if (phot%curvature%xy(ix, iy) == 0) cycle
  pt%z0 = pt%z0 - phot%curvature%xy(ix, iy) * x0**ix * y0**iy
  if (ix > 0) pt%dz_dx = pt%dz_dx - ix * phot%curvature%xy(ix, iy) * x0**(ix-1) * y0**iy
  if (iy > 0) pt%dz_dy = pt%dz_dy - iy * phot%curvature%xy(ix, iy) * x0**ix * y0**(iy-1)
  if (ix > 1) coef_xx = coef_xx - ix * (ix-1) * phot%curvature%xy(ix, iy) * x0**(ix-2) * y0**iy / 2
  if (iy > 1) coef_yy = coef_yy - iy * (iy-1) * phot%curvature%xy(ix, iy) * x0**ix * y0**(iy-2) / 2
  if (ix > 0 .and. iy > 0) coef_xy = coef_xy - ix * iy * phot%curvature%xy(ix, iy) * x0**(ix-1) * y0**(iy-1)
enddo
enddo

g = phot%curvature%elliptical
if (g(3) /= 0) then
  zt = sqrt(1 - (x0 * g(1))**2 - (y0 * g(2))**2)
  pt%z0 = pt%z0 + sqrt_one(-(g(1) * x0)**2 - (g(2) * y0)**2) / g(3)
  pt%dz_dx = pt%dz_dx - x0 * g(1)**2 / (g(3) * zt)
  pt%dz_dy = pt%dz_dy - y0 * g(2)**2 / (g(3) * zt)
  coef_xx = coef_xx - (g(1)**2 / zt - (x0 * g(1)**2)**2 / zt**3) / (2 * g(3))
  coef_yy = coef_yy - (g(2)**2 / zt - (y0 * g(2)**2)**2 / zt**3) / (2 * g(3))
  coef_xy = coef_xy - (x0 * y0 * (g(1) * g(2))**2 / zt**3) / (g(3))
endif

gs = phot%curvature%spherical
if (gs /= 0) then
  zt = sqrt(1 - (x0 * gs)**2 - (y0 * gs)**2)
  pt%z0 = pt%z0 + sqrt_one(-(gs * x0)**2 - (gs * y0)**2) / gs
  pt%dz_dx = pt%dz_dx - x0 * gs**2 / (gs * zt)
  pt%dz_dy = pt%dz_dy - y0 * gs**2 / (gs * zt)
  coef_xx = coef_xx - (gs**2 / zt - (x0 * gs**2)**2 / zt**3) / (2 * gs)
  coef_yy = coef_yy - (gs**2 / zt - (y0 * gs**2)**2 / zt**3) / (2 * gs)
  coef_xy = coef_xy - (x0 * y0 * (gs * gs)**2 / zt**3) / (gs)
endif

! Correct for fact that segment is supported at the corners of the segment and the segment is flat.
! This correction only affects z0 and not the slopes

dx = phot%grid%dr(1) / 2
dy = phot%grid%dr(2) / 2
coef_xx = coef_xx * dx**2
coef_xy = coef_xy * dx * dy
coef_yy = coef_yy * dy**2
coef_diag = coef_xx + coef_yy - abs(coef_xy)

if (abs(coef_diag) > abs(coef_xx) .and. abs(coef_diag) > abs(coef_yy)) then
  pt%z0 = pt%z0 + coef_diag
else if (abs(coef_xx) > abs(coef_yy)) then
  pt%z0 = pt%z0 + coef_xx
else
  pt%z0 = pt%z0 + coef_yy
endif

end subroutine init_surface_segment 

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine parser_transfer_control_struct (con_in, con_out, lord, ix_var)
!
! Routine to transfer the information from an input control_struct (which stores
! the user input parameters) to a control_struct that will be stored in the lat%control
! or lord%control%ramp for a ramper.
!
! Input:
!   con_in    -- control_struct: Input control structure.
!   lord      -- ele_struct: Lord element associated with the control_struct.
!   ix_var    -- integer:  If an expression stack evaluates to a constant, this routine will 
!                 modify the expression stack to evaluate to the value of: 
!                     lord%control%var(ix_var) * constant
!
! Output:
!   con_out   -- control_struct: Output control structure.
!-

subroutine parser_transfer_control_struct (con_in, con_out, lord, ix_var)

type (control_struct) con_in, con_out, con0
type (ele_struct) lord
integer ix_var, is, n, iv, atype
logical err, var_found

!

con0 = con_in  ! In case con_in and con_out actual arguments are the same.

if (allocated(con0%stack)) then
  do is = 1, size(con0%stack)
    if (con0%stack(is)%type == end_stack$) exit
  enddo
  call reallocate_expression_stack(con_out%stack, is-1)

  con_out%stack = con0%stack(1:is-1)

  ! Convert variable$ type to ramper variable index if name matches an ramper variable name.

  do is = 1, size(con_out%stack)
    if (con_out%stack(is)%type == end_stack$) exit
    if (con_out%stack(is)%type /= variable$) cycle
    do iv = 1, size(lord%control%var)
      if (upcase(con_out%stack(is)%name) /= lord%control%var(iv)%name) cycle
      con_out%stack(is)%type = iv + var_offset$
      exit
    enddo
  enddo

  ! Convert a stack of an expression "expres" to "express * control_var(1)" if "expres" does not use any control vars.
  ! Exception: Ramper slave with ran() or ran_gauss() in "expres".

  var_found = .false.
  do is = 1, size(con_out%stack)
    atype = con_out%stack(is)%type
    if (atype == ran$ .or. atype == ran_gauss$) var_found = .true. 
    if (.not. is_attribute (atype, all_control_var$)) cycle
    if (atype == end_stack$) exit
    var_found = .true.
    exit
  enddo

  if (.not. var_found) then
    if (size(con_out%stack) == 1 .and. con_out%stack(1)%name == '1' .or. con_out%stack(1)%name == '1.0') then
      con_out%stack(1) = expression_atom_struct(lord%control%var(ix_var)%name, ix_var+var_offset$, 0.0_rp)
    else
      n = size(con_out%stack)
      call reallocate_expression_stack(con_out%stack, n+2)
      con_out%stack(n+1) = expression_atom_struct(lord%control%var(ix_var)%name, ix_var+var_offset$, 0.0_rp)
      con_out%stack(n+2) = expression_atom_struct('', times$, 0.0_rp)
    endif
  endif

  ! Evaluate any variable values.

  do is = 1, size(con_out%stack)
    select case (con_out%stack(is)%type)
    case (ran$, ran_gauss$)
      if (lord%key == ramper$) cycle
      call parser_error ('RANDOM NUMBER FUNCITON MAY NOT BE USED WITH A GROUP, OR OVERLAY', &
                         'FOR ELEMENT: ' // lord%name)
      if (global_com%exit_on_error) call err_exit
      return
    case (variable$)
      call word_to_value (con_out%stack(is)%name, lord%branch%lat, con_out%stack(is)%value, err)
      if (err) then
        call parser_error ('ERROR CONVERTING WORD TO VALUE: ' // con_out%stack(is)%name, &
                           'FOR ELEMENT: ' // lord%name)
        return
      endif
      ! Variables in the arithmetic expression are immediately evaluated and never reevaluated.
      ! If the variable is an element attribute (looks like: "ele_name[attrib_name]") then this may
      ! be confusing if the attribute value changes later. To avoid some (but not all) confusion, 
      ! turn the variable into a numeric$ so the output from the type_ele routine looks "sane".
      if (index(con_out%stack(is)%name, '[') /= 0) then
        con_out%stack(is)%type = numeric$
        con_out%stack(is)%name = ''
      endif
    end select
  enddo

else
  con_out%y_knot = con0%y_knot
  if (size(con_out%y_knot) /= size(lord%control%x_knot)) then
    call parser_error ('NUMBER OF Y_SPLINE POINTS CONTROLLING PARAMETER: ' // con_out%attribute, &
                       'IS NOT THE SAME AS THE NUMBER OF X_SPLINE POINTS FOR ELEMENT: ' // lord%name)
    return
  endif
endif

end subroutine parser_transfer_control_struct 

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Function parser_fast_integer_read (int_vec, ele, delim_wanted, err_str)  result (is_ok)
!-

function parser_fast_integer_read (int_vec, ele, delim_wanted, err_str) result (is_ok)

type (ele_struct) ele
integer int_vec(:)
integer i, n, ix_word
logical is_ok, delim_found, err
character(*) err_str
character(40) word
character(1) delim_wanted, delim

!

is_ok = .false.

n = size(int_vec)
do i = 1, n
  call get_next_word (word, ix_word, delim_wanted // ' ', delim, delim_found, err_flag = err)
  if (err) then
    call parser_error ('ERROR IN ' // err_str, 'IN ELEMENT: ' // ele%name)
    return
  endif

  if (.not. is_integer(word, int_vec(i))) then
    call parser_error ('ERROR READING INTEGER IN ' // err_str, 'IN ELEMENT: ' // ele%name)
    return
  endif

  if (i == n .and. delim == delim_wanted) then
    is_ok = .true.
    return
  endif

  if (delim == delim_wanted) then
    call parser_error ('NOT ENOUGH INTEGERS FOR ' // err_str, 'IN ELEMENT: ' // ele%name)
    return
  endif
enddo

call parser_error ('EXPECTING DELIMITOR: ' // delim_wanted, &
                   'AFTER READING ' // int_str(n) // ' INTEGERS IN ' // err_str, &
                   'IN ELEMENT: ' // ele%name)

end function parser_fast_integer_read 

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Function parser_fast_complex_read (cmplx_vec, ele, delim, err_str)  result (is_ok)
!
! Routine to read an array of complex numbers. 
!
! This routine assumes that the array values are pure numbers in the form "<re>" or "(<re> <im>)" 
! where <re> and <im> are real numbers (not expressions) and there are no commas except possibly 
! at the end of the array. 
!
! Input:
!   ele             -- ele_struct: Lattice element associated with the array. Used for error messages.
!   err_str         -- character(*): String used when printing error messages identifying where in
!                         the lattice file the error is occuring.
!
! Output:
!   cmplx_vec(:)    -- complex(rp): Complex vector.
!   delim           -- character(1): Delimitor at end of array. Must be "," or "}"
!   is_ok           -- logical: True if everything OK. False otherwise.
!-

function parser_fast_complex_read (cmplx_vec, ele, delim, err_str) result (is_ok)

type (ele_struct) ele
complex(rp) cmplx_vec(:)
real(rp) rval, ival
integer i, n, ix, ios1, ios2, ix_word
logical is_ok, delim_found, err, left_parens_here
character(*) err_str
character(40) word
character(1) delim

!

is_ok = .false.
call string_trim(bp_com%parse_line, bp_com%parse_line, ix)
if (bp_com%parse_line(1:1) == '(') then
  delim = '('
  bp_com%parse_line = bp_com%parse_line(2:)
endif

n = size(cmplx_vec)
do i = 1, n
  left_parens_here = (delim == '(')
  call get_next_word (word, ix_word, ' ,()}', delim, delim_found, err_flag = err)
  if (err .or. (left_parens_here .and. delim /= ' ') .or. (.not. left_parens_here .and. delim == ')')) then
    call parser_error ('ERROR IN ' // err_str, 'IN ELEMENT: ' // ele%name)
    return
  endif

  read (word, *, iostat = ios1) rval

  if (left_parens_here) then
    call get_next_word (word, ix_word, ',()', delim, delim_found, err_flag = err)
    read (word, *, iostat = ios2) ival
    if (ios1 /= 0 .or. ios2 /= 0 .or. err .or. delim /= ')') then
      call parser_error ('ERROR READING COMPLEX VALUE IN ' // err_str, 'IN ELEMENT: ' // ele%name)
      return
    endif

    cmplx_vec(i) = cmplx(rval, ival)

    call string_trim(bp_com%parse_line, bp_com%parse_line, ix)
    if (bp_com%parse_line(1:1) == ',' .or. bp_com%parse_line(1:1) == '}' .or. bp_com%parse_line(1:1) == '(') then
      delim = bp_com%parse_line(1:1)
      bp_com%parse_line = bp_com%parse_line(2:)
    else
      delim = ' '
    endif

  else
    if (ios1 /= 0 .or. err .or. delim == ')') then
      call parser_error ('ERROR READING VALUE IN ' // err_str, 'IN ELEMENT: ' // ele%name)
      return
    endif
    cmplx_vec(i) = rval
  endif

  if ((delim == ',' .or. delim == '}') .and. i == n) then
    is_ok = .true.
    return
  endif

  if (i == n .or. (delim == ',' .or. delim == '}')) then
    call parser_error ('NOT ENOUGH VALUES FOR ' // err_str, 'IN ELEMENT: ' // ele%name)
    return
  endif
enddo

call parser_error ('EXPECTING DELIMITOR: "," or "}"', &
                   'AFTER READING ' // int_str(n) // ' COMPLEXS IN ' // err_str, &
                   'IN ELEMENT: ' // ele%name)

end function parser_fast_complex_read 

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Function parser_fast_real_read (real_vec, ele, end_delims, delim, err_str, exact_size, n_real)  result (is_ok)
!
! Routine to read an array of real numbers. 
!
! This routine assumes that the array values are pure numbers in the form "<re1> <re2> ...," 
! where <re1>, <re2>, etc. are real numbers (not expressions) and there are no commas except possibly, 
! at the end of the array.
!
! Input:
!   ele             -- ele_struct: Lattice element associated with the array. Used for error messages.
!   end_delims      -- character(*): List of possible ending delimitors.
!   err_str         -- character(*): String used when printing error messages identifying where in
!                         the lattice file the error is occuring.
!   exact_size      -- logical, optional: If True, 
!
! Output:
!   real_vec(:)     -- complex(rp): Complex vector.
!   delim           -- character(1): Delimitor at end of array.
!   is_ok           -- logical: True if everything OK. False otherwise.
!   n_real          -- integer, optional: Number of elements found
!-

function parser_fast_real_read (real_vec, ele, end_delims, delim, err_str, exact_size, n_real) result (is_ok)

type (ele_struct) ele
real(rp) real_vec(:)
real(rp) val

integer, optional :: n_real
integer i, n, ix, ios, ix_word

logical is_ok, delim_found, err, exact
logical, optional :: exact_size

character(*) err_str, end_delims
character(40) word
character(1) delim

!

is_ok = .false.
call string_trim(bp_com%parse_line, bp_com%parse_line, ix)
exact = logic_option(.true., exact_size)

n = size(real_vec)
do i = 1, n
  call get_next_word (word, ix_word, ' ' // end_delims, delim, delim_found, err_flag = err)
  if (err) then
    call parser_error ('ERROR IN ' // err_str, 'IN ELEMENT: ' // ele%name)
    return
  endif
  if (.not. delim_found) then
    call parser_error ('ERROR IN ' // err_str, 'MISSING END DELIMITOR', 'IN ELEMENT: ' // ele%name)
    return
  endif

  read (word, *, iostat = ios) val

  if (ios /= 0 .or. err) then
    call parser_error ('ERROR READING VALUE IN ' // err_str, 'IN ELEMENT: ' // ele%name)
    return
  endif
  real_vec(i) = val

  if (delim /= ' ' .and. (.not. exact .or. i == n)) then
    is_ok = .true.
    if (present(n_real)) n_real = i
    return
  endif

  if ((exact .and. i == n) .or. delim /= ' ') then
    call parser_error ('NOT ENOUGH VALUES FOR ' // err_str, 'IN ELEMENT: ' // ele%name)
    return
  endif
enddo

call parser_error ('EXPECTING DELIMITOR TO BE ONE OF: ' // end_delims, &
                   'AFTER READING ' // int_str(n) // ' REALS IN ' // err_str, &
                   'IN ELEMENT: ' // ele%name)

end function parser_fast_real_read 

end module
