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
!                         attribute to be set is something like B_FIELD_ERR. If this routine is being
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

use bmad_parser_mod, dummy => parser_set_attribute
use wall3d_mod, dummy2 => parser_set_attribute
use random_mod
       
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
type (photon_surface_struct), pointer :: surf
type (cylindrical_map_struct), pointer :: cl_map
type (cartesian_map_term1_struct), pointer :: ct_term
type (cartesian_map_term1_struct), allocatable :: ct_terms(:)
type (grid_field_struct), pointer :: g_field
type (taylor_field_struct), pointer :: t_field
type (cartesian_map_struct), pointer :: ct_map
type (ac_kicker_struct), pointer :: ac

real(rp) kx, ky, kz, tol, value, coef, r_vec(10), r0(2)
real(rp), allocatable :: table(:,:)
real(rp), pointer :: r_ptr

integer i, i2, j, k, n, nn, ix_word, how, ix_word1, ix_word2, ios, ix, i_out, ix_coef, switch
integer expn(6), ix_attrib, i_section, ix_v, ix_sec, i_ptr, i_term, ib, ie, im
integer ix_bounds(2), iy_bounds(2), i_vec(2), n_sec, key

character(40) :: word, str_ix, attrib_word, word2, name
character(40), allocatable :: name_list(:)
character(1) delim, delim1, delim2
character(80) str, err_str, line

logical delim_found, err_flag, logic, set_done, end_of_file, do_evaluate, hetero_list
logical is_attrib, err_flag2, old_style_input, ok
logical, optional :: check_free, heterogeneous_ele_list, set_field_master

! Get next WORD.
! If an overlay or group element then word is just an attribute to control
! [except for a "GROUP[COMMAND] = 0.343" redef construct]

err_flag = .true.  ! assume the worst
call get_next_word (word, ix_word, ':, =(){', delim, delim_found, call_check = .true.)
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

! overlay or group

if (ele%key == overlay$ .or. ele%key == group$) then
  i = attribute_index(ele, word)       ! general attribute search

  select case (i)
  case (type$, alias$, descrip$)
    call bmad_parser_string_attribute_set (ele, word, delim, delim_found)
    err_flag = .false.
    return

  case (var$)
    if (how == redef$ .or. allocated(ele%control%var)) then
      call parser_error ('RESETTING VAR = {...} IS NOT PERMITTED', 'FOR: ' // ele%name)
      return
    endif
    call get_overlay_group_names(ele, lat, pele, delim, delim_found, .true., err_flag)
    pele%default_attrib = ele%control%var(1)%name
    return

  case (gang$)
    call get_logical_real ('GANG', ele%value(gang$), err_flag)
    return

  case (x_knot$)
    if (.not. parse_real_list2 (lat, 'ERROR PARSING X_KNOT POINTS FOR: ' // ele%name, ele%control%x_knot, n, delim, delim_found, 10, '{', ',', '}')) return
    call re_allocate(ele%control%x_knot, n)
    if (.not. expect_one_of (', ', .false., ele%name, delim, delim_found)) return
    err_flag = .false.
    return

  end select

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
endif   ! Overlay or Group

! L_pole, N_pole for wiggler/undulator are deprecated in favor of L_period, N_period.

if (ele%key == wiggler$ .or. ele%key == undulator$) then
  if (word == 'L_POLE' .or. word == 'N_POLE') then
    if (.not. expect_one_of ('=', .true., ele%name, delim, delim_found)) return
    call parse_evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag, ele = ele)
    if (err_flag) return
    if (word == 'L_POLE') then
      ele%value(l_period$) = 2 * value
    else
      ele%value(n_period$) = value / 2
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

! For historical reasons, a few paramter[...] parameters are actually in bmad_com.

key = ele%key
if (ele%key == def_parameter$ .and. word == 'APERTURE_LIMIT_ON') key = def_bmad_com$
if (ele%key == def_parameter$ .and. word == 'ELECTRIC_DIPOLE_MOMENT') key = def_bmad_com$
if (ele%key == def_parameter$ .and. word == 'PTC_CUT_FACTOR') key = def_bmad_com$
if (ele%key == def_parameter$ .and. word == 'USE_HARD_EDGE_DRIFTS') key = def_bmad_com$

if (word == 'SPINOR_POLARIZATION' .or. word == 'SPINOR_PHI' .or. word == 'SPINOR_THETA' .or. word == 'SPINOR_XI') then
  call parser_error ('DUE TO BOOKKEEPING COMPLICATIONS, THE OLD SPINOR ATTRIBUTES NO LONGER EXIST: ' // word, &
                     'PLEASE CONVERT TO SPIN_X, SPIN_Y, SPIN_Z COMPONENTS.', 'FOR ELEMENT: ' // ele%name)
  return
endif

if (word == 'N_REF_PASS') word = 'MULTIPASS_REF_ENERGY'  ! This will be ignored...

if (word == 'SPACE_CHARGE_ON') then
  call parser_error ('Note: "bmad_com[SPACE_CHARGE_ON]" has been renamed "bmad_com[HIGH_ENERGY_SPACE_CHARGE_ON]"', &
                     'Will run normally...', level = s_warn$)
  word = 'HIGH_ENERGY_SPACE_CHARGE_ON'
endif

if (word == 'COHERENT_SYNCH_RAD_ON') then
  call parser_error ('Note: "bmad_com[COHERENT_SYNCH_RAD_ON]" has been renamed "bmad_com[CSR_AND_SPACE_CHARGE_ON]"', &
                     'Will run normally...', level = s_warn$)
  word = 'CSR_AND_SPACE_CHARGE_ON'
endif

! particle_start and bmad_com element can have attributes that are not part of the element so
! Need to use pointers_to_attribute.

if (key == def_particle_start$ .or. key == def_bmad_com$) then
  name = ele%name
  if (ele%name == 'PARAMETER') name = 'BMAD_COM'

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

  call pointers_to_attribute (lat, name, word, .false., a_ptrs, err_flag, .false.)
  if (err_flag .or. size(a_ptrs) == 0) then
    call parser_error ('BAD ATTRIBUTE: ' // word, 'FOR ELEMENT: ' // ele%name)
    return
  endif

  if (ele%key == def_parameter$ .and. word == 'APERTURE_LIMIT_ON') then
    call parser_error ('SYNTAX HAS CHANGED: PARAMETER[APERTURE_LIMIT_ON] = ... NEEDS TO BE REPLACED BY BMAD_COM[APERTURE_LIMIT_ON] = ...', &
                       'THIS IS A WARNING ONLY. THE PROGRAM WILL RUN NORMALLY.', level = s_warn$)
  endif

  if (name == 'D_ORB') then
    if (.not. parse_real_list (lat, trim(ele%name) // ' ' // word, bmad_com%d_orb, .true., delim, delim_found)) return
    bp_com%extra%d_orb_set = .true.
    return
  elseif (name == 'SPACE_CHARGE_MESH_SIZE') then
    if (.not. parse_integer_list (trim(ele%name) // ' ' // word, lat, bmad_com%space_charge_mesh_size, .true., delim, delim_found)) return
    bp_com%extra%space_charge_mesh_size_set = .true.
    return
  endif

  if (associated(a_ptrs(1)%r)) then
    call parse_evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag, ele = ele) 
    if (err_flag) return
    a_ptrs(1)%r = value
    if (associated(a_ptrs(1)%r, bmad_com%max_aperture_limit))             bp_com%extra%max_aperture_limit_set              = .true.
    if (associated(a_ptrs(1)%r, bmad_com%default_ds_step))                bp_com%extra%default_ds_step_set                 = .true.
    if (associated(a_ptrs(1)%r, bmad_com%significant_length))             bp_com%extra%significant_length_set              = .true.
    if (associated(a_ptrs(1)%r, bmad_com%rel_tol_tracking))               bp_com%extra%rel_tol_tracking_set                = .true.
    if (associated(a_ptrs(1)%r, bmad_com%abs_tol_tracking))               bp_com%extra%abs_tol_tracking_set                = .true.
    if (associated(a_ptrs(1)%r, bmad_com%rel_tol_adaptive_tracking))      bp_com%extra%rel_tol_adaptive_tracking_set       = .true.
    if (associated(a_ptrs(1)%r, bmad_com%abs_tol_adaptive_tracking))      bp_com%extra%abs_tol_adaptive_tracking_set       = .true.
    if (associated(a_ptrs(1)%r, bmad_com%init_ds_adaptive_tracking))      bp_com%extra%init_ds_adaptive_tracking_set       = .true.
    if (associated(a_ptrs(1)%r, bmad_com%min_ds_adaptive_tracking))       bp_com%extra%min_ds_adaptive_tracking_set        = .true.
    if (associated(a_ptrs(1)%r, bmad_com%fatal_ds_adaptive_tracking))     bp_com%extra%fatal_ds_adaptive_tracking_set      = .true.
    if (associated(a_ptrs(1)%r, bmad_com%autoscale_amp_abs_tol))          bp_com%extra%autoscale_amp_abs_tol_set           = .true.
    if (associated(a_ptrs(1)%r, bmad_com%autoscale_amp_rel_tol))          bp_com%extra%autoscale_amp_rel_tol_set           = .true.
    if (associated(a_ptrs(1)%r, bmad_com%autoscale_phase_tol))            bp_com%extra%autoscale_phase_tol_set             = .true.
    if (associated(a_ptrs(1)%r, bmad_com%electric_dipole_moment))         bp_com%extra%electric_dipole_moment_set          = .true.
    if (associated(a_ptrs(1)%r, bmad_com%ptc_cut_factor))                 bp_com%extra%ptc_cut_factor_set                  = .true.
    if (associated(a_ptrs(1)%r, bmad_com%sad_eps_scale))                  bp_com%extra%sad_eps_scale_set                   = .true.
    if (associated(a_ptrs(1)%r, bmad_com%sad_amp_max))                    bp_com%extra%sad_amp_max_set                     = .true.

  elseif (associated(a_ptrs(1)%i)) then
    call parse_evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag, ele = ele) 
    if (err_flag) return
    if (associated(a_ptrs(1)%i, lat%particle_start%direction) .and. nint(value) /= -1 .and. nint(value) /= 1) then
      call parser_error ('VALUE OF BEAM_SART[DIRECTION] MUST BE -1 OR 1.')
      return
    endif
    a_ptrs(1)%i = nint(value)
    if (associated(a_ptrs(1)%i, bmad_com%taylor_order))                   bp_com%extra%taylor_order_set                    = .true.
    if (associated(a_ptrs(1)%i, bmad_com%default_integ_order))            bp_com%extra%default_integ_order_set             = .true.
    if (associated(a_ptrs(1)%i, bmad_com%ptc_max_fringe_order))           bp_com%extra%ptc_max_fringe_order_set            = .true.
    if (associated(a_ptrs(1)%i, bmad_com%runge_kutta_order))              bp_com%extra%runge_kutta_order_set               = .true.
    if (associated(a_ptrs(1)%i, bmad_com%sad_n_div_max))                  bp_com%extra%sad_n_div_max_set                   = .true.
    if (associated(a_ptrs(1)%i, bmad_com%max_num_runge_kutta_step))       bp_com%extra%max_num_runge_kutta_step_set        = .true.

  elseif (associated(a_ptrs(1)%l)) then
    call parser_get_logical (word, a_ptrs(1)%l, ele%name, delim, delim_found, err_flag)
    if (err_flag) return
    if (associated(a_ptrs(1)%l, bmad_com%rf_phase_below_transition_ref))  bp_com%extra%rf_phase_below_transition_ref_set   = .true.
    if (associated(a_ptrs(1)%l, bmad_com%use_hard_edge_drifts))           bp_com%extra%use_hard_edge_drifts_set            = .true.
    if (associated(a_ptrs(1)%l, bmad_com%sr_wakes_on))                    bp_com%extra%sr_wakes_on_set                     = .true.
    if (associated(a_ptrs(1)%l, bmad_com%lr_wakes_on))                    bp_com%extra%lr_wakes_on_set                     = .true.
    if (associated(a_ptrs(1)%l, bmad_com%mat6_track_symmetric))           bp_com%extra%mat6_track_symmetric_set            = .true.
    if (associated(a_ptrs(1)%l, bmad_com%auto_bookkeeper))                bp_com%extra%auto_bookkeeper_set                 = .true.
    if (associated(a_ptrs(1)%l, bmad_com%csr_and_space_charge_on))        bp_com%extra%csr_and_space_charge_on_set         = .true.
    if (associated(a_ptrs(1)%l, bmad_com%spin_tracking_on))               bp_com%extra%spin_tracking_on_set                = .true.
    if (associated(a_ptrs(1)%l, bmad_com%backwards_time_tracking_on))     bp_com%extra%backwards_time_tracking_on_set      = .true.
    if (associated(a_ptrs(1)%l, bmad_com%radiation_damping_on))           bp_com%extra%radiation_damping_on_set            = .true.
    if (associated(a_ptrs(1)%l, bmad_com%radiation_fluctuations_on))      bp_com%extra%radiation_fluctuations_on_set       = .true.
    if (associated(a_ptrs(1)%l, bmad_com%conserve_taylor_maps))           bp_com%extra%conserve_taylor_maps_set            = .true.
    if (associated(a_ptrs(1)%l, bmad_com%absolute_time_tracking_default)) bp_com%extra%absolute_time_tracking_default_set  = .true.
    if (associated(a_ptrs(1)%l, bmad_com%convert_to_kinetic_momentum))    bp_com%extra%convert_to_kinetic_momentum_set     = .true.
    if (associated(a_ptrs(1)%l, bmad_com%aperture_limit_on))              bp_com%extra%aperture_limit_on_set               = .true.
    if (associated(a_ptrs(1)%l, bmad_com%ptc_print_info_messages))        bp_com%extra%ptc_print_info_messages_set         = .true.
    if (associated(a_ptrs(1)%l, bmad_com%debug))                          bp_com%extra%debug_set                           = .true.

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

  if (err_flag .or. .not. associated(a_ptr%r)) then
    call parser_error ('BAD ATTRIBUTE: ' // word, 'FOR ELEMENT: ' // ele%name)
    return
  endif

  call parse_evaluate_value (trim(ele%name) // ' ' // word, value, lat, delim, delim_found, err_flag, ele = ele)
  a_ptr%r = value
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
case ('X_RAY_LINE_LEN')
  call parser_error ('X_RAY_LINE_LEN NO LONGER SUPPORTED. WILL BE IGNORED. PLEASE USE A PHOTON_FORK ELEMENT INSTEAD', level = s_warn$)
  err_flag = .false.
  return

case ('TILT')
  if (ele%key == sbend$ .or. ele%key == rbend$) then
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

word = parser_translate_attribute_name (ele%key, word)

ix_attrib = attribute_index(ele, word, attrib_word)
if (attrib_free_problem(word)) return

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
  allocate (ac%frequencies(n))
  do i = 1, n
    ac%frequencies(i)%f    = table(i,1)
    ac%frequencies(i)%amp  = table(i,2)
    ac%frequencies(i)%phi  = table(i,3)
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
          call parser_get_logical (word, section%absolute_vertices_input, ele%name, delim, delim_found, err_flag)
          if (err_flag) return

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

if (attrib_word == 'SURFACE') then
  surf => ele%photon%surface

  if (.not. expect_this ('=', .true., .true., 'AFTER "SURFACE"', ele, delim, delim_found)) return
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)

  ! "ele1[surface] = ele2[surface]" construct

  if (delim == '[') then
    ele2 => parser_find_ele_for_attrib_transfer ('SURFACE', word)
    if (.not. associated(ele2%photon)) then
      call parser_error ('NO SURFACE ASSOCIATED WITH LATTICE ELEMENT: ' // word)
      return
    endif
    ele%photon%surface = ele2%photon%surface
    return
  endif

  !

  if (.not. expect_this ('{', .true., .true., 'AFTER "SURFACE"', ele, delim, delim_found)) return

  surface_loop: do

    ! Expect "GRID ={" or "TYPE ="

    call get_next_word (word, ix_word, '{}=,()', delim, delim_found)

    select case (word)

    case ('GRID')

      if (.not. expect_this ('={', .true., .true., 'AFTER "GRID"', ele, delim, delim_found)) return
      ix_bounds = int_garbage$; iy_bounds = int_garbage$

      do
        call get_next_word (word, ix_word, '{}=,()', delim, delim_found)
        if (word /= 'PT') then
          if (.not. expect_this ('=', .true., .false., 'AFTER ' // trim(word) // ' IN SURFACE CONSTRUCT', ele, delim, delim_found)) return
        endif

        select case (word)
        case ('DR')
          if (.not. parse_real_list (lat, trim(ele%name) // ' GRID DR', surf%grid%dr, .true., delim, delim_found)) return

        case ('R0')
          if (.not. parse_real_list (lat, trim(ele%name) // ' GRID R0', surf%grid%r0, .true., delim, delim_found)) return

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
            if (allocated (surf%grid%pt)) deallocate (surf%grid%pt)
            allocate (surf%grid%pt(ix_bounds(1):ix_bounds(2), iy_bounds(1):iy_bounds(2)))
          endif

        case ('PT')
          bp_com%parse_line = delim // bp_com%parse_line
          if (.not. parse_integer_list (trim(ele%name) // ' GRID PT', lat, i_vec, .true., delim, delim_found)) return
          if (.not. allocated(surf%grid%pt)) then
            call parser_error ('SURFACE PT_MAX MISSING', 'FOR: ' // ele%name)
            return
          endif
          if (any(i_vec < 0) .or. any(i_vec > ubound(surf%grid%pt))) then
            call parser_error ('SURFACE PT(I,J) INDEX OUT OF BOUNDS', 'FOR: ' // ele%name)
            return
          endif
          if (.not. expect_this ('=', .false., .false., 'IN GRID PT', ele, delim, delim_found)) return
          if (.not. parse_real_list (lat, trim(ele%name) // 'IN GRID PT', r_vec(1:4), .true., delim, delim_found)) return
          surf%grid%pt(i_vec(1), i_vec(2))%orientation = surface_orientation_struct(r_vec(1), r_vec(2), r_vec(3), r_vec(4))

        case ('TYPE')
          call get_switch ('SURFACE GRID TYPE', surface_grid_type_name(1:), surf%grid%type, err_flag2, ele, delim, delim_found)
          if (err_flag2) return
          bp_com%parse_line = delim // bp_com%parse_line

        case default
          call parser_error ('GRID COMPONENT NOT RECOGNIZED: ' // word, 'FOR ELEMENT: ' // ele%name)
          return
        end select

        if (.not. expect_one_of (',}', .false., ele%name, delim, delim_found)) return
        if (delim == '}') then
          if (.not. expect_one_of (',}', .false., ele%name, delim, delim_found)) return
          exit
        endif

      enddo

    case default
      call parser_error ('UNKNOWN SURFACE COMPONENT: ' // word2, 'FOR: ' // ele%name)
      return
    end select

    if (delim == '}') exit

  enddo surface_loop

  if (.not. expect_one_of(', ', .false., ele%name, delim, delim_found)) return
  err_flag = .false.
  return
endif

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
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)

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

  if (word(1:8) == 'binary::') then
    call parser_file_stack('push', word(9:), err = err_flag, open_file = .false.); if (err_flag) return
    call read_binary_cartesian_map(word(9:), ele, ele%cartesian_map(i_ptr), err_flag)
    call parser_file_stack('pop')
    if (err_flag) then
      call parser_error ('ERROR READING BINARY CARTESIAN_MAP FILE.')
      return
    endif
  else
    if (.not. expect_this ('{', .true., .true., 'AFTER "CARTESIAN_MAP"', ele, delim, delim_found)) return
    allocate (ele%cartesian_map(i_ptr)%ptr)
    call parse_cartesian_map(ele%cartesian_map(i_ptr), ele, lat, delim, delim_found, err_flag)
  endif

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

  if (word(1:8) == 'binary::') then
    call parser_file_stack('push', word(9:), err = err_flag, open_file = .false.); if (err_flag) return
    call read_binary_cylindrical_map(word(9:), ele, ele%cylindrical_map(i_ptr), err_flag)
    call parser_file_stack('pop')
    if (err_flag) then
      call parser_error ('ERROR READING BINARY CYLINDIRCAL_MAP FILE.')
      return
    endif
  else
    if (.not. expect_this ('{', .true., .true., 'AFTER "CYLINDRICAL_MAP"', ele, delim, delim_found)) return
    allocate (ele%cylindrical_map(i_ptr)%ptr)
    cl_map => ele%cylindrical_map(i_ptr)
    if (ele%key == lcavity$ .or. ele%key == rfcavity$) cl_map%harmonic = 1 ! Default
    call parse_cylindrical_map(cl_map, ele, lat, delim, delim_found, err_flag)
  endif

  if (ele%key == wiggler$ .or. ele%key == undulator$) ele%field_calc = fieldmap$
  return
endif

!-------------------------------
! grid_field field

if (attrib_word == 'GRID_FIELD') then

  if (.not. expect_this ('=', .true., .true., 'AFTER "GRID_FIELD"', ele, delim, delim_found)) return
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)

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

  if (word(1:6) /= 'hdf5::') then
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

  if (word(1:8) == 'binary::') then
    call parser_file_stack('push', word(9:), err = err_flag, open_file = .false.); if (err_flag) return
    call read_binary_grid_field(word(9:), ele, ele%grid_field(i_ptr), err_flag)
    call parser_file_stack('pop')
    if (err_flag) then
      call parser_error ('ERROR READING BINARY GRID_FIELD FILE.')
      return
    endif
  elseif (word(1:6) == 'hdf5::') then
    call parser_file_stack('push', word(7:), err = err_flag, open_file = .false.); if (err_flag) return
    call hdf5_read_grid_field(word(7:), ele, ele%grid_field, err_flag, combine = .true.)
    call parser_file_stack('pop')
    if (err_flag) then
      call parser_error ('ERROR READING BINARY GRID_FIELD FILE.')
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
! Taylor_field field

if (attrib_word == 'TAYLOR_FIELD') then

  if (.not. expect_this ('=', .true., .true., 'AFTER "TAYLOR_FIELD"', ele, delim, delim_found)) return
  call get_next_word (word, ix_word, '[],(){}', delim, delim_found, call_check = .true.)

  ! "ele1[taylor_field] = ele2[taylor_field]" construct

  if (delim == '[') then
    ele2 => parser_find_ele_for_attrib_transfer ('TAYLOR_FIELD', word)
    if (err_flag) return
    if (.not. associated(ele2%taylor_field)) then
      call parser_error ('NO TAYLOR_FIELD ASSOCIATED WITH LATTICE ELEMENT: ' // word)
      return
    endif
    call transfer_fieldmap(ele2, ele, taylor_field$)
    return
  endif

  if (associated(ele%taylor_field)) then
    i_ptr = size(ele%taylor_field) + 1
    ele0%taylor_field => ele%taylor_field
    allocate(ele%taylor_field(i_ptr))
    do i = 1, i_ptr-1
     ele%taylor_field(i) = ele0%taylor_field(i)
    enddo
  else
    allocate(ele%taylor_field(1))
    i_ptr = 1
  endif

  if (word(1:8) == 'binary::') then
    call parser_file_stack('push', word(9:), err = err_flag, open_file = .false.); if (err_flag) return
    call read_binary_taylor_field(word(9:), ele, ele%taylor_field(i_ptr), err_flag)
    call parser_file_stack('pop')
    if (err_flag) then
      call parser_error ('ERROR READING BINARY TAYLOR_FIELD FILE.')
      return
    endif
  else
    if (.not. expect_this ('{', .true., .true., 'AFTER "TAYLOR_FIELD"', ele, delim, delim_found)) return
    allocate (ele%taylor_field(i_ptr)%ptr)
    t_field => ele%taylor_field(i_ptr)

    call parse_taylor_field(t_field, ele, lat, delim, delim_found, err_flag)
  endif

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
    ky = sign_of(ky) * sqrt(kx**2 + kz**2)

    if (old_style_input) then
      if (ct_term%kx == 0) ct_term%kx = 1d-30  ! Something small to prevent divide by zero problems.
    endif

  elseif (abs(ky**2 + kx**2 - kz**2) < tol) then
    ct_term%form = hyper_xy$
    kz = sign_of(kz) * sqrt(kx**2 + ky**2)

    if (old_style_input) then
      ct_term%coef = ct_term%coef * ct_term%kz / ct_term%ky
      if (ct_term%kx == 0) ct_term%kx = 1d-30  ! Something small to prevent divide by zero problems.
      if (ct_term%ky == 0) ct_term%ky = 1d-30  ! Something small to prevent divide by zero problems.
    endif

  elseif (abs(ky**2 - kx**2 + kz**2) < tol) then
    ct_term%form = hyper_x$
    kx = sign_of(kx) * sqrt(ky**2 + kz**2)

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
      ele%value(tilt$) = pi / 4
      return
    case (sextupole$) 
      ele%value(tilt$) = pi / 6
      return
    case (octupole$) 
      ele%value(tilt$) = pi / 8
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
      ele%value(fint$) = 0.5
      return
    case (fintx$)
      ele%value(fintx$) = 0.5
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

case ('PTC_MAX_FRINGE_ORDER')
  call parser_error ('PLEASE CONVERT "PARAMETER[PTC_MAX_FRINGE_ORDER]" TO "BMAD_COM[PTC_MAX_FRINGE_ORDER]"', level = s_warn$)
  call parser_get_integer (bmad_com%ptc_max_fringe_order, word, ix_word, delim, delim_found, err_flag); if (err_flag) return
  bp_com%extra%ptc_max_fringe_order_set = .true.

case ('TAYLOR_ORDER')
  call parser_get_integer (ix, word, ix_word, delim, delim_found, err_flag); if (err_flag) return
  if (ix <= 0) then
    call parser_error ('TAYLOR_ORDER IS LESS THAN 1')
    return
  endif
  ptc_com%taylor_order_saved = ix
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

case ('USE_HARD_EDGE_DRIFTS')
  call parser_error ('PLEASE CONVERT "PARAMETER[USE_HARD_EDGE_DRIFTS]" TO "BMAD_COM[USE_HARD_EDGE_DRIFTS]"', level = s_warn$)
  call parser_get_logical (attrib_word, bmad_com%use_hard_edge_drifts, ele%name, delim, delim_found, err_flag)
  bp_com%extra%use_hard_edge_drifts_set = .true.

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

case ('ABSOLUTE_TIME_TRACKING')
  call parser_get_logical (attrib_word, lat%absolute_time_tracking, ele%name, delim, delim_found, err_flag); if (err_flag) return

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

case ('HIGH_ENERGY_SPACE_CHARGE_ON')
  call get_logical_real (attrib_word, ele%value(high_energy_space_charge_on$), err_flag); if (err_flag) return
  j = nint(ele%value(ix_branch$)) 
  if (j >= 0) lat%branch(j)%param%high_energy_space_charge_on = is_true(ele%value(high_energy_space_charge_on$))

case ('HIGHER_ORDER_FRINGE_TYPE')
  call get_switch (attrib_word, higher_order_fringe_type_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(higher_order_fringe_type$) = ix

case ('INTERPOLATION')
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

case ('PTC_EXACT_MODEL')
  call parser_get_logical (attrib_word, logic, ele%name, delim, delim_found, err_flag); if (err_flag) return
  call set_ptc (exact_modeling = logic)

case ('PTC_EXACT_MISALIGN')
  call parser_get_logical (attrib_word, logic, ele%name, delim, delim_found, err_flag)
  if (err_flag) return
  call set_ptc (exact_misalign = logic)

case ('PTC_FRINGE_GEOMETRY')
  call get_switch (attrib_word, ptc_fringe_geometry_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(ptc_fringe_geometry$) = ix

case ('PTC_INTEGRATION_TYPE')
  call get_switch (attrib_word, ptc_integration_type_name(1:), ele%ptc_integration_type, err_flag, ele, delim, delim_found); if (err_flag) return

case ('PTC_FIELD_GEOMETRY')
  call get_switch (attrib_word, ptc_field_geometry_name(1:), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(ptc_field_geometry$) = ix

  if (ele%key == sbend$ .and. ix == true_rbend$) then
    call parser_error ('TRUE_RBEND IS NOT A VALID PTC_FIELD_GEOMETRY VALUE FOR AN SBEND')
    return
  endif

case ('REF_ORIGIN')
  call get_switch (attrib_word, anchor_pt_name(1:), pele%ref_pt, err_flag, ele, delim, delim_found); if (err_flag) return

case ('REF_COORDINATES')
  call get_switch (attrib_word, end_at_name(1:2), ix, err_flag, ele, delim, delim_found); if (err_flag) return
  ele%value(ref_coordinates$) = ix

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
      call ran_seed_put (nint(value))  ! init random number generator
      if (nint(value) == 0) then  ! Using system clock -> Not determinisitc.
        bp_com%extra%deterministic = 0
      else
        bp_com%extra%deterministic = 2
      endif
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
      call set_flags_for_changed_attribute (ele, a_ptr, .false.)

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

        branch => pointer_to_branch(ele%name, lat, .true.)
        if (associated(branch)) then
          branch%ele(0)%value(e_tot$) = value
          call set_flags_for_changed_attribute (branch%ele(0), branch%ele(0)%value(e_tot$), .false.)
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

        branch => pointer_to_branch(ele%name, lat, .true.)
        if (associated(branch)) then
          branch%ele(0)%value(p0c$) = value
          call set_flags_for_changed_attribute (branch%ele(0), branch%ele(0)%value(p0c$), .false.)
        endif

      case ('PC')    ! Only in def_mad_beam
        lat%ele(0)%value(p0c$) = 1d9 * value
        ele%value(e_tot$) = -1

      case ('LR_FREQ_SPREAD')
        call randomize_lr_wake_frequencies (ele, set_done)
        if (set_done) call bp_set_ran_status

      case ('N_PART')
        branch => pointer_to_branch(ele%name, lat, .true.)
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
character(*) attrib_name
logical is_problem, is_free

! Attributes may be definitely free, definitely dependent, or may be free or
! dependent depending upon the state of other element parameters.

! If not check_free then at least check if it is a dependent attribute.

is_problem = .false.

if (logic_option(.false., check_free)) then
  is_free = attribute_free (ele, attrib_name, bp_com%print_err .and. .not. heterogeneous_ele_list)
  if (.not. is_free) then
    if (.not. heterogeneous_ele_list) call parser_error ('ATTRIBUTE NOT FREE TO BE SET: ' // attrib_name, 'FOR: ' // ele%name)
    err_flag = .true.
    is_problem = .true.
  endif
else
  attrib_info = attribute_info(ele, attribute_index(ele, attrib_name))
  if (attrib_info%state == dependent$) then
    if (.not. hetero_list) then
      call parser_error ('DEPENDENT ATTRIBUTE NOT FREE TO BE SET: ' // attrib_name, 'FOR: ' // ele%name)
    endif
    is_problem = .true.
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
