!+
! Subroutine tao_evaluate_a_datum (datum, u, tao_lat, datum_value, valid_value, why_invalid, called_from_lat_calc)
!
! Subroutine to put the proper data in the specified datum
!
! Input:
!   datum          -- Tao_data_struct: What type of datum
!   u              -- Tao_universe_struct: Which universe to use.
!   tao_lat        -- Tao_lattice_struct: Lattice to use.
!   called_from_lat_calc -- logical, optional: Default is False. If true, prevents infinite loop of this
!                             routine calling tao_lattice_calc
!
! Output:
!   datum          -- Tao_data_struct: 
!     %ix_ele_merit   -- For max/min type constraints: Place where value is max/min. 
!   datum_value   -- Real(rp): Value of the datum.
!   valid_value   -- Logical: Set false when there is a problem. Set true otherwise.
!   why_invalid   -- Character(*), optional: Tells why datum value is invalid.
!-

recursive subroutine tao_evaluate_a_datum (datum, u, tao_lat, datum_value, valid_value, why_invalid, called_from_lat_calc)

use tao_data_and_eval_mod, dummy => tao_evaluate_a_datum
use pointer_lattice, only: operator(.sub.)
use ptc_interface_mod, only: taylor_inverse
use ptc_layout_mod, only: normal_form_rd_terms
use measurement_mod, only: to_orbit_reading, to_eta_reading, ele_is_monitor
use expression_mod, only: numeric$

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct) datum
type (tao_lattice_branch_struct), pointer :: tao_branch
type (tao_building_wall_section_struct), pointer :: section
type (tao_building_wall_point_struct) :: pt
type (tao_data_struct), pointer :: dptr
type (tao_data_array_struct), allocatable :: d_array(:)
type (tao_lattice_struct), target :: tao_lat
type (tao_expression_info_struct), allocatable :: info(:)
type (lat_struct), pointer :: lat
type (tao_dynamic_aperture_struct), pointer :: da
type (aperture_scan_struct), pointer :: scan
type (normal_modes_struct) mode
type (ele_struct), pointer :: ele, ele_start, ele_ref, ele2
type (ele_pointer_struct), allocatable :: eles(:)
type (coord_struct), pointer :: orb0, orbit(:), orb
type (coord_struct) :: orb_at_s, orb1
type (bpm_phase_coupling_struct) bpm_data
type (taylor_struct), save :: taylor_save(6), taylor(6) ! Saved taylor map
type (floor_position_struct) floor
type (branch_struct), pointer :: branch
type (bunch_params_struct), pointer :: bunch_params(:)
type (bmad_normal_form_struct), pointer :: bmad_nf
type (ptc_normal_form_struct), pointer :: ptc_nf
type (taylor_struct), pointer :: taylor_ptr
type (complex_taylor_struct), pointer :: complex_taylor_ptr
type (all_pointer_struct) a_ptr
type (rad_int_branch_struct), pointer :: branch_ri, branch_6d
type (c_taylor), pointer :: phase_map
type (twiss_struct), pointer :: z0, z1, z2

real(rp) datum_value, mat6(6,6), vec0(6), angle, px, py, vec2(2)
real(rp) eta_vec(4), v_mat(4,4), v_inv_mat(4,4), a_vec(4), mc2, charge
real(rp) beta_gamma, one_pz, xi_sum, xi_diff, w0_mat(3,3), w_mat(3,3), vec3(3), value, s_len, n0(3)
real(rp) dz, dx, cos_theta, sin_theta, zz_pt, xx_pt, zz0_pt, xx0_pt, dpz, s_offset
real(rp) zz_center, xx_center, xx_wall, phase, amp, dalpha, dbeta, aa, bb, g2
real(rp) xx_a, xx_b, dxx1, dzz1, drad, ang_a, ang_b, ang_c, dphi, amp_a, amp_b
real(rp), allocatable :: value_vec(:)
real(rp), allocatable :: expression_value_vec(:)
real(rp) theta, phi, psi

complex(rp) eval(6), evec(6,6), n_eigen(6,3)
complex(rp) temp_cplx

integer i, j, jj, k, m, n, k_old, ix, ie, is, iz, ix_ele, ix_start, ix_ref, ie0, ie1
integer n_size, ix0, which, expnt(6), n_track, n_max, ix_branch, expo(6), n_da

character(*), optional :: why_invalid
character(6) expn_str
character(16) constraint
character(20) :: r_name = 'tao_evaluate_a_datum'
character(40) head_data_type, sub_data_type, data_source, name, dflt_dat_index
character(100) data_type, str
character(:), allocatable :: e_str

logical, optional :: called_from_lat_calc
logical found, valid_value, err, err1, err2, taylor_is_complex, use_real_part, term_found, ok
logical particle_lost, exterminate, printit, twiss_at_ele
logical, allocatable :: good(:)

! If does not exist

valid_value = .false.
datum%why_invalid = ''
twiss_at_ele = .true.

if (.not. datum%exists) then
  datum_value = real_garbage$
  call tao_set_invalid(datum, 'Datum does not exist.', why_invalid)
  return
endif

! To save time, don't evaluate if unnecessary when the running an optimizer.
! Exception: When there are datums that use expressions, things are 
!   complicated so don't try to save time in this case.

if (s%com%optimizer_running .and. .not. datum%useit_opt .and. .not. s%com%have_datums_using_expressions) then
  datum_value = 0
  return
endif

! See if there is a hook for this datum

if (associated(tao_hook_evaluate_a_datum_ptr)) then
  call tao_hook_evaluate_a_datum_ptr (found, datum, u, tao_lat, datum_value, valid_value, why_invalid)
  if (found) return
endif

! Set ix_ele, ix_ref, and ix_start. 
! Note: To check that a start element was set, need to look at datum%ele_start_name, not ix_start.

data_source = datum%data_source
data_type = datum%data_type  ! Needed since %data_type is a var length str so evaluating something like %data_type(1:10) is problematical
head_data_type = datum%data_type
lat => tao_lat%lat

if (head_data_type == 'null') then
  datum_value = 0
  call tao_set_invalid (datum, 'Datum data_type is set to "null".', why_invalid)
  valid_value = .false.
  return
endif

ele => tao_pointer_to_datum_ele (lat, datum%ele_name, datum%ix_ele, datum, valid_value, why_invalid)
if (.not. valid_value) return
ix_ele = -1
ix_branch = datum%ix_branch
if (associated(ele)) then
  ix_ele = tao_tracking_ele_index(ele, datum, ix_branch)
  if (ele%a%beta == 0) twiss_at_ele = .false.
endif

ele_ref => tao_pointer_to_datum_ele (lat, datum%ele_ref_name, datum%ix_ele_ref, datum, valid_value, why_invalid)
if (.not. valid_value) return
ix_ref = -1
if (associated(ele_ref)) ix_ref = tao_tracking_ele_index(ele_ref, datum)

ele_start => tao_pointer_to_datum_ele (lat, datum%ele_start_name, datum%ix_ele_start, datum, valid_value, why_invalid)
if (.not. valid_value) return
ix_start = ix_ele
if (associated(ele_start)) ix_start = tao_tracking_ele_index(ele_start, datum)

! Some inits

valid_value = .false.

datum_value = 0           ! default
datum%ix_ele_merit = -1   ! default

branch => lat%branch(ix_branch)
tao_branch => tao_lat%tao_branch(ix_branch)
orbit => tao_branch%orbit
bunch_params => tao_branch%bunch_params

n_track = branch%n_ele_track
n_max   = branch%n_ele_max

call re_allocate2 (value_vec, 0, n_track, .false.) ! Make sure is large enough if used.
call re_allocate2 (good,      0, n_track, .false.) ! Make sure is large enough if used.

ix = index(head_data_type, '.')
if (head_data_type(1:11) == 'expression:') then
  head_data_type = 'expression:'
elseif (ix /= 0) then
  sub_data_type  = head_data_type(ix+1:)
  head_data_type = head_data_type(1:ix) 
endif

select case (data_source)
case ('lat', 'beam')
  ! Valid data source
case default
  if ( head_data_type /= 'expression:') then
    call tao_set_invalid (datum, 'UNKNOWN DATA_SOURCE: ' // data_source, why_invalid, .true.)
    return
  endif
end select

if (index(head_data_type, 'stable') == 0 .and. head_data_type /= 'expression:') then
  if (associated(ele)) then
    if (orbit(ix_ele)%state /= alive$) then
      call tao_set_invalid (datum, 'UNSTABLE ORBIT AT EVALUATION POINT', why_invalid)
      return
    endif
  endif

  if (associated(ele_ref)) then
    if (orbit(ix_ref)%state /= alive$) then
      call tao_set_invalid (datum, 'UNSTABLE ORBIT AT REFERENCE EVALUATION POINT', why_invalid)
      return
    endif
  endif

  if (associated(ele_start)) then
    if (orbit(ix_start)%state /= alive$) then
      call tao_set_invalid (datum, 'UNSTABLE ORBIT AT START EVALUATION POINT', why_invalid)
      return
    endif
  endif
endif

!

if (.not. twiss_at_ele .or. (.not. tao_branch%twiss_valid .and. data_source == 'lat')) then
  select case (head_data_type)
  case ('alpha.', 'apparent_emit.', 'norm_apparent_emit.', 'beta.', 'bpm_eta.', 'bpm_phase.', 'cbar.', 'chrom.', &
        'chrom_ptc.', 'curly_h.', 'damp.', 'deta_ds.', 'emit.', 'norm_emit.', 'eta.', 'etap.', 'gamma.', &
        'phase.', 'phase_frac.', 'phase_frac_diff', 'ping_a.', 'ping_b.', 'rad_int.', 'srdt.', 'tune.')
    call tao_set_invalid (datum, 'UNSTABLE 1-TURN MATRIX', why_invalid)
    return
  end select
endif

! ele_ref must not be specified for some data types. Check this.

select case (head_data_type)
case ('wall.')
  if (datum%ele_ref_name /= '') then
    call tao_set_invalid (datum, 'SPECIFYING ELE_REF NOT VALID', why_invalid, .true.)
    return
 endif
end select


! ele_start must not be specified for some data types. Check this.

if (data_type(1:11) == 'periodic.tt' .or. data_type == 'sigma.pz') then
  if (datum%ele_start_name /= '') then
    call tao_set_invalid (datum, 'SPECIFYING ELE_START NOT VALID', why_invalid, .true.)
    return
  endif
endif

if (data_type(1:11) == 'periodic.tt') then
  if (datum%ele_ref_name /= '') then
    call tao_set_invalid (datum, 'SPECIFYING ELE_REF NOT VALID', why_invalid, .true.)
    return
  endif
endif

if (tao_branch%track_state /= moving_forward$ .and. ix_ele >= tao_branch%track_state) then
  if ((data_source == 'beam' .and. head_data_type /= 'n_particle_loss') .or. &
                         head_data_type(1:4) == 'bpm_' .or. head_data_type == 'orbit.') then
    call tao_set_invalid (datum, 'CANNOT EVALUATE DUE TO PARTICLE LOSS.', why_invalid)
    return
  endif
endif

if (data_source == 'beam' .and. .not. s%com%have_tracked_beam) then
  call tao_set_invalid (datum, 'DATA_SOURCE FOR DATUM SET TO "beam". BUT NO BEAM TRACKING HAS BEEN DONE!', why_invalid, err_level = s_warn$)
  return
endif

!-------------------------------------------------------------
! Case where evaluation point not at the end of the element.

if (head_data_type /= 'expression:' .and. (datum%s_offset /= 0 .or. datum%eval_point /= anchor_end$)) then
  if (data_source /= 'lat') then
    call tao_set_invalid (datum, 'CANNOT USE A BEAM DATA_SOURCE WITH A FINITE S_OFFSET OR EVAL_POINT = CENTER.', why_invalid, .true.)
    return
  endif

  if (datum%ele_start_name /= '') then
    call out_io (s_warn$, r_name, 'If there is an evaluation range (that is, ele_start is set), s_offset and', & 
                                  ' eval_point are ignored. For datum:' // tao_datum_name(datum))
  else
    datum_value = tao_evaluate_datum_at_s (datum, tao_lat, ele, ele_ref, valid_value, str, exterminate)
    if (.not. valid_value) call tao_set_invalid (datum, str, why_invalid, exterminate)
    return
  endif
endif

!---------------------------------------------------

select case (head_data_type)

!-----------

case ('alpha.')

  select case (data_type)

  case ('alpha.a')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%a%alpha, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%alpha, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = valid_value .and. (bunch_params(ix_ele)%a%norm_emit /= 0)
      if (.not. valid_value) call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
    endif
    
  case ('alpha.b')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%b%alpha, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%alpha, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = valid_value .and. (bunch_params(ix_ele)%b%norm_emit /= 0)
      if (.not. valid_value) call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
    endif

  case ('alpha.z')
    if (data_source == 'lat') goto 9001  ! Error message and return
    call tao_load_this_datum (bunch_params(:)%z%alpha, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID.', why_invalid, .true.)
    return

  end select

!-----------

case ('apparent_emit.', 'norm_apparent_emit.')

  select case (data_type)

  case ('apparent_emit.x', 'norm_apparent_emit.x')
    do i = ix_start, ix_ele
      if (data_source == 'lat') then
        value_vec(i) = tao_lat_emit_calc (x_plane$, apparent_emit$, branch%ele(i), tao_branch%modes_6d)
      else
        value_vec(i) = tao_beam_emit_calc (x_plane$, apparent_emit$, branch%ele(i), bunch_params(i))
      endif
    enddo

    if (ix_ref > -1) then
      if (data_source == 'lat') then
        value_vec(ix_ref) = tao_lat_emit_calc (x_plane$, apparent_emit$, branch%ele(ix_ref), tao_branch%modes_6d)
      else
        value_vec(ix_ref) = tao_beam_emit_calc (x_plane$, apparent_emit$, branch%ele(i), bunch_params(ix_ref))
      endif
    endif

    call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

    if (data_type == 'norm_apparent_emit.x') then
      beta_gamma = ele%value(p0c$) /  mass_of(branch%param%particle)
      datum_value = datum_value * beta_gamma
    endif

  case ('apparent_emit.y', 'norm_apparent_emit.y')
    do i = ix_start, ix_ele
      if (data_source == 'lat') then
        value_vec(i) = tao_lat_emit_calc (y_plane$, apparent_emit$, branch%ele(i), tao_branch%modes_6d)
      else
        value_vec(i) = tao_beam_emit_calc (y_plane$, apparent_emit$, branch%ele(i), bunch_params(i))
      endif
    enddo

    if (ix_ref > -1) then
      if (data_source == 'lat') then
        value_vec(ix_ref) = tao_lat_emit_calc (y_plane$, apparent_emit$, branch%ele(ix_ref), tao_branch%modes_6d)
        call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      else
        value_vec(ix_ref) = tao_beam_emit_calc (y_plane$, apparent_emit$, branch%ele(i), bunch_params(ix_ref))
        call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      endif
    endif


    if (data_type == 'norm_apparent_emit.y') then
      beta_gamma = ele%value(p0c$) /  mass_of(branch%param%particle)
      datum_value = datum_value * beta_gamma
    endif


  case default
    call tao_set_invalid (datum, 'UNKNOWN DATUM TYPE: ' // data_type, why_invalid, .true.)
    return

  end select

!-----------

case ('beta.')

  select case (data_type)

  case ('beta.x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%x%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      if (bunch_params(ix_ele)%x%norm_emit == 0) then
        valid_value = .false.
        call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
      endif
    else
      call tao_set_invalid (datum, 'DATA_SOURCE: ' // trim(data_source) // ' INVALID WITH: beta.x DATA_TYPE', why_invalid, .true.)
    endif
    
  case ('beta.y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%y%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      if (bunch_params(ix_ele)%y%norm_emit == 0) then
        valid_value = .false.
        call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
      endif
    else
      call tao_set_invalid (datum, 'DATA_SOURCE: ' // trim(data_source) // ' INVALID WITH: beta.y DATA_TYPE', why_invalid, .true.)
    endif

  case ('beta.z')
    if (data_source == 'lat') goto 9001  ! Error message and return
    call tao_load_this_datum (bunch_params(:)%z%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    if (bunch_params(ix_ele)%z%norm_emit == 0) then
      valid_value = .false.
      call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
    endif

  case ('beta.a')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%a%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      if (bunch_params(ix_ele)%a%norm_emit == 0) then
        valid_value = .false.
        call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
      endif
    endif
    
  case ('beta.b')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%b%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      if (bunch_params(ix_ele)%b%norm_emit == 0) then
        valid_value = .false.
        call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
      endif
    endif

  case ('beta.c')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%z%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%c%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      if (bunch_params(ix_ele)%a%norm_emit == 0) then
        valid_value = .false.
        call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
      endif
    endif
    
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('bpm_orbit.')

  select case (data_type)
  case ('bpm_orbit.x')
    which = x_plane$
  case ('bpm_orbit.y')
    which = y_plane$
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  end select

  if (data_source == 'beam') orbit => bunch_params%centroid

  valid_value = .true.
  particle_lost = .false.

  do i = ix_start, ix_ele
    if (i /= ix_ele .and. .not. ele_is_monitor(branch%ele(i), .false.)) then
      value_vec(i) = 0
      cycle
    endif
    call to_orbit_reading (orbit(i), branch%ele(i), which, s%com%add_measurement_noise, value_vec(i), err)
    particle_lost = particle_lost .or. (tao_branch%track_state /= moving_forward$ .and. i > tao_branch%track_state)
    valid_value = valid_value .and. .not. err
  enddo

  if (ix_ref > -1) then
    call to_orbit_reading (orbit(ix_ref), branch%ele(ix_ref), which, s%com%add_measurement_noise, value_vec(ix_ref), err)
    particle_lost = particle_lost .or. (tao_branch%track_state /= moving_forward$ .and. ix_ref > tao_branch%track_state)
    valid_value = valid_value .and. .not. err
  endif

  if (particle_lost) then
    call tao_set_invalid (datum, 'CANNOT EVALUATE DUE TO PARTICLE LOSS', why_invalid)
    return
  elseif (.not. valid_value) then
    call tao_set_invalid (datum, 'NO VALID MONITOR ELEMENT', why_invalid, .true.)
    return
  endif

  call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)


!-----------

case ('bpm_eta.')

  select case (data_type)
  case ('bpm_eta.x')
    which = x_plane$
  case ('bpm_eta.y')
    which = y_plane$
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  end select

  if (data_source == 'beam') goto 9000  ! Set error message and return
  vec2 = [ele%x%eta, ele%y%eta]
  call to_eta_reading (vec2, ele, which, s%com%add_measurement_noise, datum_value, err)
  valid_value = .not. err

  if (ix_ref > -1) then
    vec2 = [ele_ref%x%eta, ele_ref%y%eta]
    call to_eta_reading (vec2, ele_ref, which, s%com%add_measurement_noise, value, err)
    datum_value = datum_value - value
    valid_value = valid_value .and. .not. err
  endif


!-----------

case ('bpm_phase.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value, why_invalid, datum); if (.not. valid_value) return

  select case (data_type)
  case ('bpm_phase.a')
    datum_value = bpm_data%phi_a
  case ('bpm_phase.b')
    datum_value = bpm_data%phi_b
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    valid_value = .false.
    return
  end select

  if (ix_ref > -1 .and. valid_value) then
    call tao_to_phase_and_coupling_reading (ele_ref, bpm_data, valid_value, why_invalid, datum); if (.not. valid_value) return
    select case (data_type)
    case ('bpm_phase.a')
      datum_value = datum_value - bpm_data%phi_a
    case ('bpm_phase.b')
      datum_value = datum_value - bpm_data%phi_b
    end select
  endif

!-----------

case ('bpm_k.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value, why_invalid, datum); if (.not. valid_value) return

  select case (data_type)
  case ('bpm_k.22a')
    datum_value = bpm_data%k_22a
  case ('bpm_k.12a')
    datum_value = bpm_data%k_12a
  case ('bpm_k.11b')
    datum_value = bpm_data%k_11b
  case ('bpm_k.12b')
    datum_value = bpm_data%k_12b
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    valid_value = .false.
    return
  end select

  if (ix_ref > -1 .and. valid_value) then
    call tao_to_phase_and_coupling_reading (ele_ref, bpm_data, valid_value, why_invalid, datum); if (.not. valid_value) return
    select case (data_type)
    case ('bpm_k.22a')
      datum_value = datum_value - bpm_data%k_22a
    case ('bpm_k.12a')
      datum_value = datum_value - bpm_data%k_12a
    case ('bpm_k.11b')
      datum_value = datum_value - bpm_data%k_11b
    case ('bpm_k.12b')
      datum_value = datum_value - bpm_data%k_12b
    end select
  endif

!-----------

case ('bpm_cbar.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value, why_invalid, datum); if (.not. valid_value) return

  select case (data_type)
  case ('bpm_cbar.22a')
    datum_value = bpm_data%cbar22_a
  case ('bpm_cbar.12a')
    datum_value = bpm_data%cbar12_a
  case ('bpm_cbar.11b')
    datum_value = bpm_data%cbar11_b
  case ('bpm_cbar.12b')
    datum_value = bpm_data%cbar12_b
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    valid_value = .false.
    return
  end select

  if (ix_ref > -1 .and. valid_value) then
    call tao_to_phase_and_coupling_reading (ele_ref, bpm_data, valid_value, why_invalid, datum); if (.not. valid_value) return
    select case (data_type)
    case ('bpm_cbar.22a')
      datum_value = datum_value - bpm_data%cbar22_a
    case ('bpm_cbar.12a')
      datum_value = datum_value - bpm_data%cbar12_a
    case ('bpm_cbar.11b')
      datum_value = datum_value - bpm_data%cbar11_b
    case ('bpm_cbar.12b')
      datum_value = datum_value - bpm_data%cbar12_b
    end select
  endif

!-----------

case ('bunch_charge.')

  call tao_load_this_datum (bunch_params(:)%charge_live, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  if (data_type == 'bunch_charge.live_relative') then
    charge = bunch_params(ele%ix_ele)%charge_tot
    if (charge == 0) then
      call tao_set_invalid (datum, 'BUNCH HAS NO CHARGE FOR EVALUATING A DATUM OF TYPE "bunch_charge_live.percent', why_invalid)
      valid_value = .false.
      return
    endif
    datum_value = datum_value / charge

  elseif (data_type /= 'bunch_charge.live') then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    valid_value = .false.
    return
  endif

!-----------

case ('bunch_max.', 'bunch_min.')
  if (data_source /= 'beam') goto 9000  ! Set error message and return
  select case (data_type(11:))
  case ('x');  i = 1
  case ('px'); i = 2
  case ('y');  i = 3
  case ('py'); i = 4
  case ('z');  i = 5
  case ('pz'); i = 6
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  end select
  
  select case (data_type(1:10))
  case ('bunch_max.')
    call tao_load_this_datum (bunch_params(:)%rel_max(i), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
  case ('bunch_min.')
    call tao_load_this_datum (bunch_params(:)%rel_min(i), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
  end select

!-----------

case ('c_mat.', 'cmat.')

  if (data_type(1:5) == 'c_mat') then
    data_type = 'cmat' // data_type(6:)
    call out_io (s_warn$, r_name, 'Note: "c_mat" data type is now called "cmat"')
  endif

  select case (data_type)

  case ('cmat.11')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_load_this_datum (branch%ele(:)%c_mat(1,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('cmat.12')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_load_this_datum (branch%ele(:)%c_mat(1,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('cmat.21')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_load_this_datum (branch%ele(:)%c_mat(2,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('cmat.22')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_load_this_datum (branch%ele(:)%c_mat(2,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('cbar.')

  select case (data_type)

  case ('cbar.11')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%cbar(1,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('cbar.12')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%cbar(1,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('cbar.21')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%cbar(2,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('cbar.22')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%cbar(2,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('chrom.')
  
  if (data_source == 'beam') goto 9000  ! Set error message and return

  ! Can happen with a command like "show lat -attribute chrom.a" that the chromaticity has not been computed.
  ! So try to compute it if needed.
  if (.not. tao_lat%chrom_calc_ok) then
    ok = .false.
    if (.not. logic_option(.false., called_from_lat_calc)) then ! Try calling tao_lattice_calc.
      s%com%force_chrom_calc = .true.
      s%u%calc%lattice = .true.
      call tao_lattice_calc(ok)
    endif
  endif

  if (.not. allocated(tao_lat%low_E_lat%branch) .or. .not. tao_lat%chrom_calc_ok) then
    if (branch%param%unstable_factor == 0) then
      call tao_set_invalid (datum, 'Chrom calculation problem.', why_invalid)
    else
      call tao_set_invalid (datum, 'Unstable lattice.', why_invalid)
    endif
    return
  endif

  !----

  select case (data_type)

  case ('chrom.dtune.a', 'chrom.a')
    if (data_type == 'chrom.dtune.a') call out_io (s_warn$, r_name, '"chrom.dtune.a" IS DEPRECATED. PLEASE CHANGE TO "chrom.a".')
    datum_value = tao_branch%a%chrom
    valid_value = .true.

  case ('chrom.dtune.b', 'chrom.b')
    if (data_type == 'chrom.dtune.b') call out_io (s_warn$, r_name, '"chrom.dtune.b" IS DEPRECATED. PLEASE CHANGE TO "chrom.b".')
    datum_value = tao_branch%b%chrom
    valid_value = .true.

  case ('chrom.dbeta.a')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = tao_lat%lat%ele(i)%a%dbeta_dpz / tao_lat%lat%ele(i)%a%beta
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.dbeta.b')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = tao_lat%lat%ele(i)%b%dbeta_dpz / tao_lat%lat%ele(i)%b%beta
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif
  
  case ('chrom.dphi.a')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        dpz = tao_branch%high_E_orb(i)%vec(6) - tao_branch%low_E_orb(i)%vec(6)
        value_vec(i) = (tao_lat%high_E_lat%branch(ix_branch)%ele(i)%a%phi - tao_lat%low_E_lat%branch(ix_branch)%ele(i)%a%phi)/ dpz
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.dphi.b')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        dpz = tao_branch%high_E_orb(i)%vec(6) - tao_branch%low_E_orb(i)%vec(6)
        value_vec(i) = (tao_lat%high_E_lat%branch(ix_branch)%ele(i)%b%phi - tao_lat%low_E_lat%branch(ix_branch)%ele(i)%b%phi)/ dpz
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.deta.x')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        dpz = tao_branch%high_E_orb(i)%vec(6) - tao_branch%low_E_orb(i)%vec(6)
        value_vec(i) = tao_lat%lat%ele(i)%x%deta_dpz
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.deta.y')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        dpz = tao_branch%high_E_orb(i)%vec(6) - tao_branch%low_E_orb(i)%vec(6)
        value_vec(i) = tao_lat%lat%ele(i)%y%deta_dpz
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.detap.x')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        dpz = tao_branch%high_E_orb(i)%vec(6) - tao_branch%low_E_orb(i)%vec(6)
        value_vec(i) = tao_lat%lat%ele(i)%x%detap_dpz
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.detap.y')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        dpz = tao_branch%high_E_orb(i)%vec(6) - tao_branch%low_E_orb(i)%vec(6)
        value_vec(i) = tao_lat%lat%ele(i)%y%detap_dpz
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.w.a', 'chrom.w.b')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        if (data_type == 'chrom.w.a') then
          z0 => branch%ele(i)%a
        else
          z0 => branch%ele(i)%b
        endif
        bb = z0%dbeta_dpz / z0%beta
        aa = z0%dalpha_dpz - z0%alpha * bb
        value_vec(i) = sqrt(aa**2 + bb**2)
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif      

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('chrom_ptc.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  ptc_nf => tao_branch%ptc_normal_form

  if (.not. ptc_nf%valid_map) then
    if (.not. u%calc%one_turn_map) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE ONE_TURN_MAP_CALC IS NOT SET TO TRUE.', why_invalid)
    elseif (branch%param%geometry /= closed$) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE LATTICE GEOMETRY IS NOT CLOSED.', why_invalid)
    else
      call tao_set_invalid (datum, '?????', why_invalid)
    endif
    return
  endif

  select case (data_type(1:12))
  case ('chrom_ptc.a.')
    phase_map => ptc_nf%phase(1)
  case ('chrom_ptc.b.')
    phase_map => ptc_nf%phase(2)
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  end select

  if (.not. is_integer(data_type(13:), n)) then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  endif

  expo = 0
  expo(6) = n 

  datum_value = real(phase_map .sub. expo)
  valid_value = .true.

!-----------

case ('curly_h.')

  if (data_source == 'beam')  goto 9000  ! Set error message and return

  select case (data_type)
  case ('curly_h.a')
    do i = ix_start, ix_ele
      ele => branch%ele(i)
      value_vec(i) = ele%a%gamma * ele%a%eta**2 + 2 * ele%a%alpha * ele%a%eta * ele%a%etap + ele%a%beta * ele%a%etap**2
    enddo
    if (ix_ref > -1) then
      ele => branch%ele(ix_ref)
      value_vec(ix_ref) = ele%a%gamma * ele%a%eta**2 + 2 * ele%a%alpha * ele%a%eta * ele%a%etap + ele%a%beta * ele%a%etap**2
    endif
    call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('curly_h.b')
    do i = ix_start, ix_ele
      ele => branch%ele(i)
      value_vec(i) = ele%b%gamma * ele%b%eta**2 + 2 * ele%b%alpha * ele%b%eta * ele%b%etap + ele%b%beta * ele%b%etap**2
    enddo
    if (ix_ref > -1) then
      ele => branch%ele(ix_ref)
      value_vec(ix_ref) = ele%b%gamma * ele%b%eta**2 + 2 * ele%b%alpha * ele%b%eta * ele%b%etap + ele%b%beta * ele%b%etap**2
    endif
    call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  end select

!-----------

case ('damp.')

  if (data_source == 'beam') goto 9000  ! Set error message and return

  select case (data_type)

  case ('damp.j_a')
    datum_value = tao_branch%modes_6d%a%j_damp
    valid_value = .true.

  case ('damp.j_b')
    datum_value = tao_branch%modes_6d%b%j_damp
    valid_value = .true.

  case ('damp.j_z')
    datum_value = tao_branch%modes_6d%z%j_damp
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('deta_ds.')

  select case (data_type)

  case ('deta_ds.a')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%deta_ds, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%a%deta_ds, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('deta_ds.b')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%deta_ds, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%b%deta_ds, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('deta_ds.x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%x%deta_ds, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%x%deta_ds, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('deta_ds.y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%y%deta_ds, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%y%deta_ds, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('dpx_dx') 
  if (data_source == 'lat') then
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) = tao_branch%lat_sigma(ix_ref)%mat(1,2) / tao_branch%lat_sigma(ix_ref)%mat(1,1)
      value_vec(ix_ele) = tao_branch%lat_sigma(ix_ele)%mat(1,2) / tao_branch%lat_sigma(ix_ele)%mat(1,1)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (tao_branch%lat_sigma%mat(1,2) / tao_branch%lat_sigma%mat(1,1), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif
  else
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) = bunch_params(ix_ref)%sigma(1,2) / bunch_params(ix_ref)%sigma(1,1)
      value_vec(ix_ele) = bunch_params(ix_ele)%sigma(1,2) / bunch_params(ix_ele)%sigma(1,1)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    else
      call tao_load_this_datum (bunch_params%sigma(1,2) / bunch_params%sigma(1,1), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif
  endif

case ('dpy_dy') 
  if (data_source == 'lat') then
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) =  tao_branch%lat_sigma(ix_ref)%mat(3,4) / tao_branch%lat_sigma(ix_ref)%mat(3,3)
      value_vec(ix_ele) = tao_branch%lat_sigma(ix_ele)%mat(3,4) / tao_branch%lat_sigma(ix_ele)%mat(3,3)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (tao_branch%lat_sigma%mat(3,4) / tao_branch%lat_sigma%mat(3,3), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif
  else
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) =  bunch_params(ix_ref)%sigma(3,4) / bunch_params(ix_ref)%sigma(3,3)
      value_vec(ix_ele) = bunch_params(ix_ele)%sigma(3,4) / bunch_params(ix_ele)%sigma(3,3)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    else
      call tao_load_this_datum (bunch_params%sigma(3,4) / bunch_params%sigma(3,3), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif
  endif

case ('dpz_dz') 
  if (data_source == 'lat') then
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) =  tao_branch%lat_sigma(ix_ref)%mat(5,6) / tao_branch%lat_sigma(ix_ref)%mat(5,5)
      value_vec(ix_ele) = tao_branch%lat_sigma(ix_ele)%mat(5,6) / tao_branch%lat_sigma(ix_ele)%mat(5,5)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (tao_branch%lat_sigma%mat(5,6) / tao_branch%lat_sigma%mat(5,5), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif
  else
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) =  bunch_params(ix_ref)%sigma(5,6) / bunch_params(ix_ref)%sigma(5,5)
      value_vec(ix_ele) = bunch_params(ix_ele)%sigma(5,6) / bunch_params(ix_ele)%sigma(5,5)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    else
      call tao_load_this_datum (bunch_params%sigma(5,6) / bunch_params%sigma(5,5), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif
  endif

!-----------

case ('dynamic_aperture.')
  da => u%dynamic_aperture
  if (allocated(da%scan)) then
    if (.not. allocated(da%scan(1)%point)) then
      call tao_set_invalid (datum, 'DYNAMIC APERTURE NOT CALCULATED', why_invalid)
      return
    endif
  else
    call tao_set_invalid (datum, 'DYNAMIC APERTURE NOT CALCULATED', why_invalid)
    return
  endif

  if (da%a_emit <= 0 .or. da%b_emit <= 0) then
    call tao_set_invalid (datum, 'A_EMIT OR B_EMIT NOT SET IN TAO_DYNAMIC_APERTURE STRUCTURE.', why_invalid)
    return
  endif

  n_da = size(da%scan)

  if (.not. is_integer(sub_data_type, n)) then
    call tao_set_invalid (datum, 'MALFORMED DATA_TYPE: ' // quote(data_type) // '. ' // quote(sub_data_type) // ' IS NOT AN INTEGER.', why_invalid, .true.)
    return
  endif

  if (n < 1 .or. n > n_da) then
    call tao_set_invalid (datum, 'SCAN INDEX OUT OF RANGE FOR DATA_TYPE: ' // quote(data_type), why_invalid, .true.)
    return
  endif

  scan => da%scan(n)
  if (da%param%start_ele == '') then
    ele => lat%ele(0)
  else
    call lat_ele_locator (da%param%start_ele, lat, eles, n)
    ele => eles(1)%ele
  endif

  datum_value = 1d100   ! Something large
  do j = 1, size(scan%point)
    orb1 = scan%ref_orb
    orb1%vec(1:4) = [scan%point(j)%x, 0.0_rp, scan%point(j)%y, 0.0_rp]
    call orbit_amplitude_calc (ele, orb1, amp_a, amp_b)
    amp = sqrt(2 * (amp_a / da%a_emit + amp_b / da%b_emit))
    datum_value = min(datum_value, amp)
  enddo

  valid_value = .true.

!-----------

case ('e_tot_ref')
  if (data_source == 'beam') goto 9000  ! Set error message and return
  call tao_load_this_datum (branch%ele(:)%value(e_tot$), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('element_attrib.')

  name = upcase(data_type(16:))
  value_vec = 0
  good = .false.

  do i = ix_start, ix_ele
    call pointer_to_attribute (branch%ele(i), name, .false., a_ptr, err, .false.)
    if (.not. associated (a_ptr%r)) cycle
    value_vec(i) = a_ptr%r
    good(i) = .true.
  enddo

  if (all(.not. good(ix_start:ix_ele))) then
    call tao_set_invalid (datum, 'CANNOT EVALUATE DATUM WITH DATA_TYPE = "' // trim(data_type) // '" AT ASSOCIATED ELEMENT', why_invalid, .true.)
    return
  endif

  if (ix_ref > -1) then
    call pointer_to_attribute (ele_ref, name, .false., a_ptr, err, .false.)
    if (associated (a_ptr%r)) then
      value_vec(ix_ref) = a_ptr%r
      good(ix_ref) = .true.
    else
      if (.not. good(ix_ref)) then
        call tao_set_invalid (datum, 'CANNOT EVALUATE DATUM WITH DATA_TYPE = "' // trim(data_type) // '" AT REFERENCE ELEMENT', why_invalid, .true.)
        return
      endif
    endif
  endif

  call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, good)

!-----------

case ('emit.', 'norm_emit.')

  if (associated(ele)) then
    beta_gamma = ele%value(p0c$) / mass_of(branch%param%particle)
  else
    beta_gamma = 0
  endif

  select case (data_type)

  case ('emit.x', 'norm_emit.x')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = tao_lat_emit_calc (x_plane$, projected_emit$, branch%ele(i), tao_branch%modes_6d)
      enddo
      if (ix_ref > -1) then
        value_vec(ix_ref) = tao_lat_emit_calc (x_plane$, projected_emit$, branch%ele(ix_ref), tao_branch%modes_6d)
      endif
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%x%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif

    if (data_type == 'norm_emit.x') datum_value = datum_value * beta_gamma

  case ('emit.y', 'norm_emit.y')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = tao_lat_emit_calc (y_plane$, projected_emit$, branch%ele(i), tao_branch%modes_6d)
      enddo
      if (ix_ref > -1) then
        value_vec(ix_ref) = tao_lat_emit_calc (y_plane$, projected_emit$, branch%ele(ix_ref), tao_branch%modes_6d)
      endif
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%y%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif

    if (data_type == 'norm_emit.y') datum_value = datum_value * beta_gamma

  case ('emit.z', 'norm_emit.z')
    if (data_source == 'lat') then
      goto 9001   ! Error message and return
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%z%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif

    if (data_type == 'norm_emit.y') datum_value = datum_value * beta_gamma

  case ('emit.a', 'norm_emit.a')
    if (data_source == 'lat') then
      if (lat%param%geometry == open$ .and. ix_ele > -1) then
        s%com%rad_int_ri_calc_on = .true.
        if (.not. allocated(tao_lat%rad_int_by_ele_ri%branch)) then
          call out_io (s_error$, r_name, 'tao_lat%rad_int_by_ele_ri%branch not allocated')
          return
        endif
        
        call tao_load_this_datum (tao_lat%rad_int_by_ele_ri%branch(ix_branch)%ele%lin_norm_emit_a, &
                                ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
        datum_value = datum_value / beta_gamma
      else
        datum_value = tao_branch%modes_6d%a%emittance
        valid_value = .true.
      endif
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif

    if (data_type == 'norm_emit.a') datum_value = datum_value * beta_gamma
    
  case ('emit.b', 'norm_emit.b')
    if (data_source == 'lat') then
      if (lat%param%geometry == open$ .and. ix_ele > -1) then
        s%com%rad_int_ri_calc_on = .true.
        if (.not. allocated(tao_lat%rad_int_by_ele_ri%branch)) then
          call out_io (s_error$, r_name, 'tao_lat%rad_int_by_ele_ri%branch not allocated')
          return
        endif
        call tao_load_this_datum (tao_lat%rad_int_by_ele_ri%branch(ix_branch)%ele%lin_norm_emit_b, &
                                ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
        datum_value = datum_value / beta_gamma
      else
        datum_value = tao_branch%modes_6d%b%emittance
        valid_value = .true.
      endif
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif

    if (data_type == 'norm_emit.b') datum_value = datum_value * beta_gamma

  case ('emit.c', 'norm_emit.c')
    if (data_source == 'lat') then
      goto 9001     ! Error message and return
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%c%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif

    if (data_type == 'norm_emit.y') datum_value = datum_value * beta_gamma

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('eta.')

  select case (data_type)

  case ('eta.a')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%a%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('eta.b')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%b%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('eta.x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%x%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%x%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('eta.y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%y%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%y%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('eta.z')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%z%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%z%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('etap.')

  select case (data_type)

  case ('etap.a')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%a%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('etap.b')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%b%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('etap.x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%x%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%x%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('etap.y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%y%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%y%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('expression:', 'expression.')

  ! Something like "expression:max(data::orbit.x|model)" needs dflt_dat_indx = "*"
  !! write (dflt_dat_index, '(i0)') datum%ix_d1    ???? why was this used previously?
  dflt_dat_index = '*'

  e_str = datum%data_type(12:)
  do
    ix = index(e_str, 'ele::#[')
    if (ix == 0) exit
    if (ix_ele == -1) then
      call tao_set_invalid (datum, 'NO ASSOCIATED ELEMENT' // e_str, why_invalid, .true.)
      return
    endif
    e_str = e_str(1:ix+4) // trim(ele_loc_name(ele)) // e_str(ix+6:)
  enddo

  printit = (s%com%n_err_messages_printed < s%global%datum_err_messages_max) 
  call tao_evaluate_expression (e_str, 0, .false., expression_value_vec, err, printit, info, &
                  datum%stack, tao_lat%name, datum%data_source, ele_ref, ele_start, ele, &
                  dflt_dat_index, u%ix_uni, datum%eval_point, datum%s_offset, datum = datum)
  if (err) then
    call tao_set_invalid (datum, 'CANNOT EVALUATE EXPRESSION: ' // e_str, why_invalid)
    return
  endif

  select case (datum%merit_type)
  case ('min');      datum_value = minval(expression_value_vec)
  case ('max');      datum_value = maxval(expression_value_vec)
  case ('abs_min');  datum_value = minval(abs(expression_value_vec))
  case ('abs_max');  datum_value = maxval(abs(expression_value_vec))
  case ('max-min');  datum_value = maxval(expression_value_vec) - minval(expression_value_vec)

  case ('integral', 'average', 'rms')
    s_offset = 0
    if (.not. associated(ele_start)) then
      call tao_set_invalid (datum, 'ELE_START NOT SET. THIS IS NEEDED WHEN MERIT_TYPE IS SET TO "integral", "average", OR "rms".', why_invalid)
      return
    endif

    do i = 1, size(info)
      j = i + ele_start%ix_ele - 1
      if (j > branch%n_ele_track) then
        j = j - branch%n_ele_track - 1
        s_offset = branch%param%total_length
      endif
      info(i)%s = tao_datum_s_position(datum, branch%ele(j)) + s_offset
    enddo

    if (j /= ele%ix_ele) then
      call out_io (s_error$, r_name, 'BOOKKEEPING ERROR IN EVALUATING INTEGRAL/AVERAGE OF EXPRESSION.', &
                                     'PLEASE REPORT!')
      call tao_set_invalid (datum, 'CANNOT EVALUATE EXPRESSION: ' // datum%data_type, why_invalid)
    endif

    datum_value = tao_datum_integrate(datum, branch, info(:)%s, expression_value_vec, valid_value, why_invalid)
    return

  case ('target')
    if (size(expression_value_vec) /= 1) then
      call tao_set_invalid (datum, 'MERIT_TYPE IS SET TO "TARGET" BUT DATUM DOES NOT EVALUATE TO A SINGLE NUMBER!', why_invalid, .true.)
      return
    endif
    datum_value = expression_value_vec(1)
  case default
    call out_io (s_error$, r_name, &
                'SINCE THIS DATUM: ' // tao_datum_name(datum), &
                'SPECIFIES A RANGE OF ELEMENTS, THEN THIS MERIT_TYPE: ' // datum%merit_type, &
                'IS NOT VALID. VALID MERIT_TYPES ARE MIN, MAX, ABS_MIN, AND ABS_MAX.')
    call tao_set_invalid (datum, 'MERIT_TYPE: ' // quote(datum%merit_type) // ' IS NOT VALID WHEN THERE IS AN EVALUATION RANGE', why_invalid, .true.)
    return
  end select

  ! Make sure that any datums used in the expression have already been evaluated.
  do i = 1, size(datum%stack)
    if (datum%stack(i)%type /= numeric$) cycle
    call tao_find_data (err, datum%stack(i)%name, d_array = d_array, print_err = .false.)
    if (err .or. size(d_array) == 0) cycle  ! Err -> This is not associated then not a datum.
    dptr => d_array(1)%d
    if (dptr%d1%d2%ix_universe < u%ix_uni) cycle ! OK
    if (dptr%d1%d2%ix_universe == u%ix_uni .and. dptr%ix_data < datum%ix_data) cycle
    call out_io (s_error$, r_name, 'DATUM: ' // tao_datum_name(datum), &
                    'WHICH IS OF TYPE EXPRESSION:' // datum%data_type, &
                    'THE EXPRESSION HAS A COMPONENT: ' // datum%stack(i)%name, &
                    'AND THIS COMPONENT IS EVALUATED AFTER THE EXPRESSION!', &
                    'TO FIX: MOVE THE EXPRESSION DATUM TO BE AFTER THE COMPONENT DATUM IN THE FILE THAT DEFINES THE DATA.')
    return
  enddo
  valid_value = .true.

!-----------

case ('floor.')

  select case (data_type)

  case ('floor.x')
    call tao_load_this_datum (branch%ele(:)%floor%r(1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('floor.y')
    call tao_load_this_datum (branch%ele(:)%floor%r(2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('floor.z')
    call tao_load_this_datum (branch%ele(:)%floor%r(3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('floor.theta')
    call tao_load_this_datum (branch%ele(:)%floor%theta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('floor.phi')
    call tao_load_this_datum (branch%ele(:)%floor%phi, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('floor.psi')
    call tao_load_this_datum (branch%ele(:)%floor%psi, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('floor_actual.')

  value_vec(ix_ele) = tao_ele_geometry_with_misalignments(datum, ele, valid_value, why_invalid)
  if (.not. valid_value) return
  if (associated(ele_ref)) value_vec(ix_ref) = tao_ele_geometry_with_misalignments(datum, ele_ref, valid_value, why_invalid)
  if (.not. valid_value) return

  if (associated(ele_start)) then
    do ie = ix_start, ix_ele - 1
      value_vec(ie) = tao_ele_geometry_with_misalignments(datum, branch%ele(ie), valid_value, why_invalid)
      if (.not. valid_value) return
    enddo
  endif

  call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('floor_orbit.')

  value_vec(ix_ele) = tao_eval_floor_orbit (datum, ele, orbit(ix_ele), bunch_params(ix_ele), valid_value, why_invalid)
  if (.not. valid_value) return
  if (associated(ele_ref)) value_vec(ix_ref) = tao_eval_floor_orbit (datum, ele_ref, orbit(ix_ref), bunch_params(ix_ref), valid_value, why_invalid)
  if (.not. valid_value) return

  if (associated(ele_start)) then
    do ie = ix_start, ix_ele - 1
      value_vec(ie) = tao_eval_floor_orbit (datum, branch%ele(ie), orbit(ie), bunch_params(ie), valid_value, why_invalid)
      if (.not. valid_value) return
    enddo
  endif

  call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('gamma.')

  select case (data_type)

  case ('gamma.a')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%a%gamma, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%gamma, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = valid_value .and. (bunch_params(ix_ele)%a%norm_emit /= 0)
      call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
    endif
    
  case ('gamma.b')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%b%gamma, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%gamma, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = valid_value .and. (bunch_params(ix_ele)%b%norm_emit /= 0)
      call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
    endif

  case ('gamma.z')
    if (data_source == 'lat') goto 9001  ! Error message and return
    call tao_load_this_datum (bunch_params(:)%z%gamma, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('k.')

  select case (data_type)

  case ('k.11b')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%k_11a, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  case ('k.12a')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%k_12a, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  case ('k.12b')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%k_12b, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  case ('k.22a')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%k_22b, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('momentum')
  if (data_source == 'beam') goto 9000  ! Set error message and return
  call tao_load_this_datum (branch%ele(0:n_track)%value(p0c$) * (1+orbit(0:n_track)%vec(6)), &
                            ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('momentum_compaction')
  if (data_source == 'beam') goto 9000  ! Set error message and return

  if (ix_ref < 0) then
    ix_ref = 0
    ele_ref => lat%branch(branch%ix_branch)%ele(ix_ref)
  endif

  g2 = (mass_of(ele_ref%ref_species) / ele_ref%value(E_tot$))**2   ! 1/gamma^2

  orb0 => orbit(ix_ref)
  call make_v_mats (ele_ref, v_mat, v_inv_mat)
  eta_vec = [ele_ref%a%eta, ele_ref%a%etap, ele_ref%b%eta, ele_ref%b%etap]
  eta_vec = matmul (v_mat, eta_vec)
  one_pz = 1 + orb0%vec(6)
  eta_vec(2) = eta_vec(2) * one_pz + orb0%vec(2) / one_pz
  eta_vec(4) = eta_vec(4) * one_pz + orb0%vec(4) / one_pz

  call transfer_matrix_calc (lat, mat6, vec0, ix_ref, ix_start, branch%ix_branch)

  do i = ix_start, ix_ele
    s_len = branch%ele(i)%s - branch%ele(ix_ref)%s
    if (s_len == 0) then
      value_vec(i) = 0
    else
      value_vec(i) = g2 - (sum(mat6(5,1:4) * eta_vec) + mat6(5,6)) / s_len
    endif
    if (i /= ix_ele) mat6 = matmul(branch%ele(i+1)%mat6, mat6)
  enddo
  call tao_load_this_datum (value_vec, null(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('momentum_compaction_ptc.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  ptc_nf => tao_branch%ptc_normal_form

  if (.not. ptc_nf%valid_map) then
    if (.not. u%calc%one_turn_map) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE ONE_TURN_MAP_CALC IS NOT SET TO TRUE.', why_invalid)
    elseif (branch%param%geometry /= closed$) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE LATTICE GEOMETRY IS NOT CLOSED.', why_invalid)
    else
      call tao_set_invalid (datum, '?????', why_invalid)
    endif
    return
  endif

  if (.not. is_integer(data_type(25:), n)) then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  endif

  expo = 0
  expo(6) = n 

  datum_value = real(ptc_nf%path_length .sub. expo) / branch%param%total_length
  valid_value = .true.

!-----------

case ('n_particle_loss')
  if (data_source /= 'beam') goto 9001  ! Set error message and return
  if (ix_ele < 0) ix_ele = branch%n_ele_track
  datum_value = sum(bunch_params(ix_ref+1:ix_ele)%n_particle_lost_in_ele)
  valid_value = .true.

!--------------

case ('normal.')

  ! Fetches normal_form components.
  if (data_source == 'beam') goto 9000  ! Set error message and return

  ptc_nf => tao_branch%ptc_normal_form
  bmad_nf => tao_branch%bmad_normal_form

  ! Do nothing it the map wasn't made
  if (.not. ptc_nf%valid_map) then
    call tao_set_invalid (datum, 'PTC one-turn map not calculated.', why_invalid)
    return
  endif


  if (.not. associated(bmad_nf%ele_origin)) then
    ! Get resonant driving terms
    call normal_form_rd_terms(ptc_nf%one_turn_map, bmad_nf, rf_is_on(branch))
    bmad_nf%ele_origin => ptc_nf%ele_origin
  endif

  ! Expect: taylor.#.######
  ! Example: normal.dhdj.2.000001 is the b-mode chromaticity
  !          head   sub
  ! Get position of first number. 
  iz = index(sub_data_type, '.') + 1
  
  if(sub_data_type(1:2) == 'h.') then
    valid_value = .true.
    term_found = .false.
    do i=1, size(bmad_nf%h)
      if(sub_data_type(3:8) == bmad_nf%h(i)%id) then
        temp_cplx = bmad_nf%h(i)%c_val
        term_found = .true.
      endif
    enddo
    if(term_found) then
      select case (sub_data_type(10:10))
      case('r')
        datum_value = real(temp_cplx)
      case('i')
        datum_value = aimag(temp_cplx)
      case('a')
        datum_value = abs(temp_cplx)
      case default
        call tao_set_invalid (datum, 'Data_type not ending in .r, .i, or .a.', why_invalid, .true.)
        valid_value = .false.
        return
      end select
    else
      call tao_set_invalid (datum, 'Data_type not found in normal_form_struct', why_invalid, .true.)
      valid_value = .false.
      return
    endif
  else
    i = tao_read_phase_space_index (sub_data_type, iz, .false.)
    if (i == 0) then
      call tao_set_invalid (datum, 'Bad phase space index.', why_invalid, .true.)
      return
    endif
    ! Point to taylor
    taylor_is_complex = .false.
    if (sub_data_type(1:5) == 'dhdj.') then
      taylor_ptr => bmad_nf%dhdj(i)
    else if (sub_data_type(1:2) == 'A.') then
      taylor_ptr => bmad_nf%A(i)
    else if (sub_data_type(1:6) == 'A_inv.') then
      taylor_ptr => bmad_nf%A_inv(i)
    else if (sub_data_type(1:2) == 'M.') then
      taylor_ptr => bmad_nf%M(i)
    else if (sub_data_type(1:4) == 'ReF.') then
      taylor_is_complex = .true.
      use_real_part = .true.
      complex_taylor_ptr => bmad_nf%f(i)
    else if (sub_data_type(1:4) == 'ImF.') then
      taylor_is_complex = .true.
      use_real_part = .false.
      complex_taylor_ptr => bmad_nf%F(i)
    else if (sub_data_type(1:4) == 'ReL.') then
      taylor_is_complex = .true.
      use_real_part = .true.
      complex_taylor_ptr => bmad_nf%L(i)
    else if (sub_data_type(1:4) == 'ImL.') then
      taylor_is_complex = .true.
      use_real_part = .false.
      complex_taylor_ptr => bmad_nf%L(i)
    endif
   
    ! Check for second dot
    if (sub_data_type(iz+1:iz+1) /= '.') then
     call tao_set_invalid (datum, 'Missing dot "." in data_type', why_invalid, .true.)
     call out_io (s_error$, r_name, 'data_type: '//trim(data_type) )
     call out_io (s_error$, r_name, 'expect dot: ', sub_data_type(1:iz)//'.######' )
    endif
   
    ! Get exponent
    expn_str = sub_data_type(iz+2:iz+7)
    expnt = 0
    do j = 1, 6
      if (expn_str(j:j) == ' ') exit
      expnt(j) = index('0123456789', expn_str(j:j)) - 1
    enddo
    
    ! Coefficient
    if (taylor_is_complex) then
      if (use_real_part) then
        datum_value = real(complex_taylor_coef(complex_taylor_ptr, expnt))
      else
        datum_value = aimag(complex_taylor_coef(complex_taylor_ptr, expnt))
      endif
    else
      datum_value = taylor_coef(taylor_ptr, expnt)
    endif
    valid_value = .true.  
  endif

!-----------

case ('orbit.')

  if (tao_branch%track_state /= moving_forward$ .and. ix_ele > tao_branch%track_state) then
    valid_value = .false.
    call tao_set_invalid (datum, 'Particle lost.', why_invalid)
    return
  endif

  select case (data_type)

  case ('orbit.energy', 'orbit.e_tot', 'orbit.kinetic')  ! orbit.e_tot is old style
    if (ix_ref > -1) then
      if (data_source == 'beam') then
        orb => bunch_params(ix_ref)%centroid
      else
        orb => orbit(ix_ref)
      endif
      if (orb%state == not_set$) goto 7000  ! Set error message and return
      value_vec(ix_ref) = (1 + orb%vec(6)) * orb%p0c / orb%beta
    endif

    do i = ix_start, ix_ele
      if (data_source == 'beam') then
        orb => bunch_params(i)%centroid
      else
        orb => orbit(i)
      endif
      if (orb%state == not_set$) goto 7000  ! Set error message and return
      call convert_pc_to ((1 + orb%vec(6))*orb%p0c, orb%species, e_tot = value_vec(i))
    enddo

    call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

    if (data_type == 'orbit.kinetic') datum_value = datum_value - mass_of(orb%species)

  case ('orbit.x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%vec(1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%vec(3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.z')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(5), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%vec(5), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.px')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%vec(2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.py')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(4), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%vec(4), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.pz')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(6), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%vec(6), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.amp_a')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%amp_a, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)

  case ('orbit.amp_b')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%amp_b, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)

  case ('orbit.norm_amp_a')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%amp_na, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)

  case ('orbit.norm_amp_b')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%amp_nb, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('pc')
  if (ix_ref > -1) then
    if (data_source == 'beam') then
      value_vec(ix_ref) = (1 + bunch_params(ix_ref)%centroid%vec(6)) * bunch_params(ix_ref)%centroid%p0c
    else
      value_vec(ix_ref) = (1 + orbit(ix_ref)%vec(6)) * orbit(ix_ref)%p0c
    endif
  endif

  if (data_source == 'beam') then
    do i = ix_start, ix_ele
      value_vec(i) = (1 + bunch_params(i)%centroid%vec(6)) * bunch_params(ix_ref)%centroid%p0c
    enddo
  else
    do i = ix_start, ix_ele
      value_vec(i) = (1 + orbit(i)%vec(6)) * orbit(ix_ref)%p0c
    enddo
  endif

  call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('periodic.')

  ix = index(data_type(10:), '.') + 9
  select case (data_type(1:ix))

  case ('periodic.tt.')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    if (lat%param%geometry /= closed$ .and. .not. associated(ele_ref)) then
      call tao_set_invalid (datum, 'LATTICE MUST BE CIRCULAR FOR A DATUM LIKE: ' // data_type, why_invalid)
      call err_exit
    endif

    call transfer_map_calc (lat, taylor, err, ix_ele, ix_ele, orbit(ix_ele), branch%ix_branch, &
                                                       one_turn = .true., concat_if_possible = s%global%concatenate_maps)
    if (err) then
      call tao_set_invalid (datum, 'CANNOT INVERT MAP', why_invalid)
      return
    endif

    do i = 1, 4
      call add_taylor_term (taylor(i), -1.0_rp, taylor_expn([i]))
    enddo
    call taylor_inverse (taylor, taylor, err)
    if (err) then
      call tao_set_invalid (datum, 'CANNOT INVERT MAP', why_invalid)
      return
    endif

    expnt = 0
    i = tao_read_phase_space_index (data_type, 13, .false.)
    if (i == 0) then
      call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
      return
    endif

    do j = 14, 24
      if (data_type(j:j) == ' ') exit
      k = tao_read_phase_space_index (data_type, j, .false.)
      if (k == 0) then
        call tao_set_invalid (datum, 'BAD DATA_TYPE = "' // trim(data_type), why_invalid, .true.)
        return
      endif
      expnt(k) = expnt(k) + 1
    enddo

    datum_value = taylor_coef (taylor(i), expnt)
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('phase.', 'phase_frac.')

  select case (data_type)

  case ('phase.a', 'phase_frac.a')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    if (ix_ref < 0) then
      datum_value = ele%a%phi
    else
      datum_value = ele%a%phi - ele_ref%a%phi
      if (ix_ref > ix_ele .and. branch%param%geometry == closed$) then
        dphi = branch%ele(n_track)%a%phi - branch%ele(0)%a%phi
        if (2*datum_value < -dphi) datum_value = datum_value + dphi
      endif
    if (data_type == 'phase_frac.a') datum_value = modulo2(datum_value, pi)
    endif
    valid_value = .true.

  case ('phase.b', 'phase_frac.b')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    if (ix_ref < 0) then
      datum_value = ele%b%phi
    else
      datum_value = ele%b%phi - ele_ref%b%phi
      if (ix_ref > ix_ele .and. branch%param%geometry == closed$) then
        dphi = branch%ele(n_track)%b%phi - branch%ele(0)%b%phi
        if (2*datum_value < -dphi) datum_value = datum_value + dphi
      endif
    endif
    if (data_type == 'phase_frac.b') datum_value = modulo2(datum_value, pi)
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('phase_frac_diff')
  if (data_source == 'beam') goto 9000  ! Set error message and return

  if (ix_ref < 0) then
    px = ele%a%phi 
    py = ele%b%phi 
  else
    px = ele%a%phi - ele_ref%a%phi
    py = ele%b%phi - ele_ref%b%phi
    if (ix_ref > ix_ele .and. branch%param%geometry == closed$) then
      dphi = branch%ele(n_track)%a%phi - branch%ele(0)%a%phi
      if (2*px < -dphi) px = px + dphi
      dphi = branch%ele(n_track)%b%phi - branch%ele(0)%b%phi
      if (2*px < -dphi) py = py + dphi
    endif
  endif

  datum_value = modulo2 (px - py, pi)
  valid_value = .true.

!-----------

case ('photon.')

  select case (data_type)

  case ('photon.intensity_x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%centroid%field(1)**2, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%field(1)**2, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('photon.intensity_y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%centroid%field(2)**2, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%field(2)**2, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('photon.intensity')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%centroid%field(1)**2+bunch_params(:)%centroid%field(2)**2, &
                                                                      ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%field(1)**2 + orbit(:)%field(2)**2, ele_ref, ele_start, ele, &
                                                                           datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('photon.phase_x')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_load_this_datum (orbit(:)%phase(1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('photon.phase_y')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_load_this_datum (orbit(:)%phase(2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
  end select

!-----------

case ('ping_a.')
  select case (data_type)
  case ('ping_a.amp_x')
    datum_value = ele%gamma_c * sqrt(ele%a%beta)
    if (associated(ele_ref)) datum_value = datum_value - ele%gamma_c * sqrt(ele_ref%a%beta)
    valid_value = .true.

  case ('ping_a.phase_x')
    datum_value = ele%a%phi
    if (associated(ele_ref)) datum_value = datum_value - ele_ref%a%phi
    valid_value = .true.

  case ('ping_a.amp_y')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    datum_value = sqrt(ele%b%beta * (scratch%cc(ix_ele)%cbar(1,2)**2 + scratch%cc(ix_ele)%cbar(2,2)**2))
    if (associated(ele_ref)) then
      datum_value = datum_value - sqrt(ele_ref%b%beta * (scratch%cc(ix_ref)%cbar(1,2)**2 + scratch%cc(ix_ref)%cbar(2,2)**2))
    endif
    valid_value = .true.

  case ('ping_a.phase_y')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    datum_value = ele%a%phi + atan2(scratch%cc(ix_ele)%cbar(1,2), - scratch%cc(ix_ele)%cbar(2,2))
    if (associated(ele_ref)) then
      datum_value = datum_value - ele_ref%a%phi - atan2(scratch%cc(ix_ref)%cbar(1,2), - scratch%cc(ix_ref)%cbar(2,2))
    endif
    valid_value = .true.

  case ('ping_a.amp_sin_y')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    amp = sqrt(ele%b%beta * (scratch%cc(ix_ele)%cbar(1,2)**2 + scratch%cc(ix_ele)%cbar(2,2)**2))
    phase = ele%a%phi + atan2(scratch%cc(ix_ele)%cbar(1,2), -scratch%cc(ix_ele)%cbar(2,2))
    datum_value = amp * sin(phase)
    if (associated(ele_ref)) then
      amp = sqrt(ele_ref%b%beta * (scratch%cc(ix_ref)%cbar(1,2)**2 + scratch%cc(ix_ref)%cbar(2,2)**2))
      phase = ele_ref%a%phi + atan2(scratch%cc(ix_ref)%cbar(1,2), -scratch%cc(ix_ref)%cbar(2,2))
      datum_value = datum_value - amp * sin(phase)
    endif
    valid_value = .true.

  case ('ping_a.amp_cos_y')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    amp = sqrt(ele%b%beta * (scratch%cc(ix_ele)%cbar(1,2)**2 + scratch%cc(ix_ele)%cbar(2,2)**2))
    phase = ele%a%phi + atan2(scratch%cc(ix_ele)%cbar(1,2), -scratch%cc(ix_ele)%cbar(2,2))
    datum_value = amp * cos(phase)
    if (associated(ele_ref)) then
      amp = sqrt(ele_ref%b%beta * (scratch%cc(ix_ref)%cbar(1,2)**2 + scratch%cc(ix_ref)%cbar(2,2)**2))
      phase = ele_ref%a%phi + atan2(scratch%cc(ix_ref)%cbar(1,2), -scratch%cc(ix_ref)%cbar(2,2))
      datum_value = datum_value - amp * cos(phase)
    endif
    valid_value = .true.

  case ('ping_a.amp_sin_rel_y')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    datum_value = -sqrt(ele%b%beta) * scratch%cc(ix_ele)%cbar(1,2)

    if (associated(ele_ref)) then
      datum_value = datum_value + sqrt(ele_ref%b%beta) * scratch%cc(ix_ref)%cbar(1,2)
    endif
    valid_value = .true.

  case ('ping_a.amp_cos_rel_y')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    datum_value = -sqrt(ele%b%beta) * scratch%cc(ix_ele)%cbar(2,2)

    if (associated(ele_ref)) then
      datum_value = datum_value + sqrt(ele_ref%b%beta) * scratch%cc(ix_ref)%cbar(2,2)
    endif
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
  end select

!-----------

case ('ping_b.')
  select case (data_type)
  case ('ping_b.amp_y')
    datum_value = ele%gamma_c * sqrt(ele%b%beta)
    if (associated(ele_ref)) datum_value = datum_value - ele%gamma_c * sqrt(ele_ref%b%beta)
    valid_value = .true.

  case ('ping_b.phase_y')
    datum_value = ele%b%phi
    if (associated(ele_ref)) datum_value = datum_value - ele_ref%b%phi
    valid_value = .true.

  case ('ping_b.amp_x')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    datum_value = sqrt(ele%a%beta * (scratch%cc(ix_ele)%cbar(1,2)**2 + scratch%cc(ix_ele)%cbar(1,1)**2))
    if (associated(ele_ref)) then
      datum_value = datum_value - sqrt(ele_ref%a%beta * (scratch%cc(ix_ref)%cbar(1,2)**2 + scratch%cc(ix_ref)%cbar(1,1)**2))
    endif
    valid_value = .true.

  case ('ping_b.phase_x')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    datum_value = ele%b%phi + atan2(scratch%cc(ix_ele)%cbar(1,2), scratch%cc(ix_ele)%cbar(1,1))
    if (associated(ele_ref)) then
      datum_value = datum_value - ele_ref%b%phi - atan2(scratch%cc(ix_ref)%cbar(1,2), scratch%cc(ix_ref)%cbar(1,1))
    endif
    valid_value = .true.

  case ('ping_b.amp_sin_x')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    amp = sqrt(ele%a%beta * (scratch%cc(ix_ele)%cbar(1,2)**2 + scratch%cc(ix_ele)%cbar(1,1)**2))
    phase = ele%b%phi + atan2(scratch%cc(ix_ele)%cbar(1,2), scratch%cc(ix_ele)%cbar(1,1))
    datum_value = amp * sin(phase)
    if (associated(ele_ref)) then
      amp = sqrt(ele_ref%a%beta * (scratch%cc(ix_ref)%cbar(1,2)**2 + scratch%cc(ix_ref)%cbar(1,1)**2))
      phase = ele_ref%b%phi + atan2(scratch%cc(ix_ref)%cbar(1,2), scratch%cc(ix_ref)%cbar(1,1))
      datum_value = datum_value - amp * sin(phase)
    endif
    valid_value = .true.

  case ('ping_b.amp_cos_x')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    amp = sqrt(ele%a%beta * (scratch%cc(ix_ele)%cbar(1,2)**2 + scratch%cc(ix_ele)%cbar(1,1)**2))
    phase = ele%b%phi + atan2(scratch%cc(ix_ele)%cbar(1,2), scratch%cc(ix_ele)%cbar(1,1))
    datum_value = amp * cos(phase)
    if (associated(ele_ref)) then
      amp = sqrt(ele_ref%a%beta * (scratch%cc(ix_ref)%cbar(1,2)**2 + scratch%cc(ix_ref)%cbar(1,1)**2))
      phase = ele_ref%b%phi + atan2(scratch%cc(ix_ref)%cbar(1,2), scratch%cc(ix_ref)%cbar(1,1))
      datum_value = datum_value - amp * cos(phase)
    endif
    valid_value = .true.

  case ('ping_b.amp_sin_rel_x')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    datum_value = -sqrt(ele%a%beta) * scratch%cc(ix_ele)%cbar(1,2)
    if (associated(ele_ref)) then
      datum_value = datum_value + sqrt(ele_ref%a%beta) * scratch%cc(ix_ref)%cbar(1,2)
    endif
    valid_value = .true.

  case ('ping_b.amp_cos_rel_x')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    datum_value = sqrt(ele%a%beta) * scratch%cc(ix_ele)%cbar(1,1)
    if (associated(ele_ref)) then
      datum_value = datum_value - sqrt(ele_ref%a%beta) * scratch%cc(ix_ref)%cbar(1,1)
    endif
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
  end select

!-----------

case ('r.')
  if (data_source == 'beam') goto 9000  ! Set error message and return

  i = tao_read_phase_space_index (data_type, 3, .false.)
  j = tao_read_phase_space_index (data_type, 4, .false.)
  if (i == 0 .or. j == 0 .or. len_trim(data_type) /= 4) then
    call tao_set_invalid (datum, 'BAD DATA_TYPE = "' // trim(data_type), why_invalid, .true.)
    return
  endif

  if (ix_ref < 0) ix_ref = 0
  call transfer_matrix_calc (lat, mat6, vec0, ix_ref, ix_start, branch%ix_branch)

  k = ix_start
  do 
    value_vec(k) = mat6(i, j)
    if (k == ix_ele) exit
    k = k + 1
    if (k > n_track) k = 0
    mat6 = matmul(branch%ele(k)%mat6, mat6)
  enddo

  call tao_load_this_datum (value_vec, NULL(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('r56_compaction')
  if (data_source == 'beam') goto 9000  ! Set error message and return

  if (ix_ref < 0) then
    ix_ref = 0
    ele_ref => lat%branch(branch%ix_branch)%ele(ix_ref)
  endif

  orb0 => orbit(ix_ref)
  call make_v_mats (ele_ref, v_mat, v_inv_mat)
  eta_vec = [ele_ref%a%eta, ele_ref%a%etap, ele_ref%b%eta, ele_ref%b%etap]
  eta_vec = matmul (v_mat, eta_vec)
  one_pz = 1 + orb0%vec(6)
  eta_vec(2) = eta_vec(2) * one_pz + orb0%vec(2) / one_pz
  eta_vec(4) = eta_vec(4) * one_pz + orb0%vec(4) / one_pz

  call transfer_matrix_calc (lat, mat6, vec0, ix_ref, ix_start, branch%ix_branch)

  do i = ix_start, ix_ele
    value_vec(i) = sum(mat6(5,1:4) * eta_vec) + mat6(5,6)
    if (i /= ix_ele) mat6 = matmul(branch%ele(i+1)%mat6, mat6)
  enddo
  call tao_load_this_datum (value_vec, null(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('rad_int.')

  if (ix_ref > -1 .or. ix_ele > -1) then
    if (ix_ele < 0) ix_ele = branch%n_ele_track
    if (ix_ref < 0) ix_ref = 0
  endif

  if (data_source == 'beam') goto 9000  ! Set error message and return

  if (.not. tao_lat%rad_int_calc_ok .or. .not. tao_lat%emit_6d_calc_ok) then
    if (.not. logic_option(.false., called_from_lat_calc)) then ! Try calling tao_lattice_calc.
      s%com%force_rad_int_calc = .true.
      u%calc%lattice = .true.
      call tao_lattice_calc(ok)
    endif

    if (.not. tao_lat%rad_int_calc_ok .or. .not. tao_lat%emit_6d_calc_ok) then
      call tao_set_invalid (datum, 'Radiation integral calc failed.', why_invalid)
      return
    endif
  endif

  branch_ri => tao_lat%rad_int_by_ele_ri%branch(ix_branch)
  branch_6d => tao_lat%rad_int_by_ele_6d%branch(ix_branch)

  select case (data_type)
  case ('rad_int.i0')
    s%com%rad_int_6d_calc_on = .true.
    if (ix_ele > -1) then
      datum_value = sum(branch_6d%ele(ix_ref:ix_ele)%i0)
    else
      datum_value = tao_branch%modes_6d%synch_int(0)
    endif

  case ('rad_int.i1')
    s%com%rad_int_6d_calc_on = .true.
    if (ix_ele > -1) then
      datum_value = sum(branch_6d%ele(ix_ref:ix_ele)%i1)
    else
      datum_value = tao_branch%modes_6d%synch_int(1)
    endif

  case ('rad_int.i2')
    s%com%rad_int_6d_calc_on = .true.
    if (ix_ele > -1) then
      datum_value = sum(branch_6d%ele(ix_ref:ix_ele)%i2)
    else
      datum_value = tao_branch%modes_6d%synch_int(2)
    endif

  case ('rad_int.i3')
    s%com%rad_int_6d_calc_on = .true.
    if (ix_ele > -1) then
      datum_value = sum(branch_6d%ele(ix_ref:ix_ele)%i3)
    else
      datum_value = tao_branch%modes_6d%synch_int(3)
    endif

  case ('rad_int.i2_e4')
    s%com%rad_int_ri_calc_on = .true.
    if (ix_ele > -1) then
      datum_value = sum(branch_ri%ele(ix_ref:ix_ele)%lin_i2_e4)
    else
      datum_value = tao_branch%modes_ri%lin%i2_e4
    endif

  case ('rad_int.i3_e7')
    s%com%rad_int_ri_calc_on = .true.
    if (ix_ele > -1) then
      datum_value = sum(branch_ri%ele(ix_ref:ix_ele)%lin_i3_e7)
    else
      datum_value = tao_branch%modes_ri%lin%i3_e7
    endif

  case ('rad_int.i4a')
    s%com%rad_int_ri_calc_on = .true.
    if (ix_ele > -1) then
      datum_value = sum(branch_ri%ele(ix_ref:ix_ele)%i4a)
    else
      datum_value = tao_branch%modes_ri%a%synch_int(4)
    endif

  case ('rad_int.i4b')
    s%com%rad_int_ri_calc_on = .true.
    if (ix_ele > -1) then
      datum_value = sum(branch_ri%ele(ix_ref:ix_ele)%i4b)
    else
      datum_value = tao_branch%modes_ri%b%synch_int(4)
    endif

  case ('rad_int.i4z')
    s%com%rad_int_ri_calc_on = .true.
    if (ix_ele > -1) then
      datum_value = sum(branch_ri%ele(ix_ref:ix_ele)%i4z)
    else
      datum_value = tao_branch%modes_ri%z%synch_int(4)
    endif

  case ('rad_int.i5a')
    s%com%rad_int_ri_calc_on = .true.
    if (ix_ele > -1) then
      datum_value = sum(branch_ri%ele(ix_ref:ix_ele)%i5a)
    else
      datum_value = tao_branch%modes_ri%a%synch_int(5)
    endif

  case ('rad_int.i5a_e6')
    s%com%rad_int_ri_calc_on = .true.
    if (ix_ele > -1) then
      datum_value = sum(branch_ri%ele(ix_ref:ix_ele)%lin_i5a_e6)
    else
      datum_value = tao_branch%modes_ri%lin%i5a_e6
    endif

  case ('rad_int.i5b')
    s%com%rad_int_ri_calc_on = .true.
    if (ix_ele > -1) then
      datum_value = sum(branch_ri%ele(ix_ref:ix_ele)%i5b)
    else
      datum_value = tao_branch%modes_ri%b%synch_int(5)
    endif

  case ('rad_int.i5b_e6')
    s%com%rad_int_ri_calc_on = .true.
    if (ix_ele > -1) then
      datum_value = sum(branch_ri%ele(ix_ref:ix_ele)%lin_i5b_e6)
    else
      datum_value = tao_branch%modes_ri%lin%i5b_e6
    endif

  case ('rad_int.i6b')
    s%com%rad_int_ri_calc_on = .true.
    if (ix_ele > -1) then
      datum_value = sum(branch_ri%ele(ix_ref:ix_ele)%i6b)
    else
      datum_value = tao_branch%modes_ri%b%synch_int(6)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

  valid_value = .true.

!-----------
! Basically, 'rad_int1_ri.' and 'rad_int1_6d.' are used for internal use.

case ('rad_int1.', 'rad_int1_ri.', 'rad_int1_6d.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  if (ix_ele < 0) return

  select case (head_data_type)
  case ('rad_int1.')
    branch_ri => tao_lat%rad_int_by_ele_ri%branch(ix_branch)
    branch_6d => tao_lat%rad_int_by_ele_6d%branch(ix_branch)
  case ('rad_int1_ri.')
    branch_ri => tao_lat%rad_int_by_ele_ri%branch(ix_branch)
    branch_6d => tao_lat%rad_int_by_ele_ri%branch(ix_branch)
    data_type = 'rad_int1.' // data_type(13:)
  case ('rad_int1_6d.')
    branch_ri => tao_lat%rad_int_by_ele_6d%branch(ix_branch)
    branch_6d => tao_lat%rad_int_by_ele_6d%branch(ix_branch)
    data_type = 'rad_int1.' // data_type(13:)
  end select

  select case (data_type)
  case ('rad_int1.i0')
    datum_value = branch_6d%ele(ix_ele)%i0
    if (ix_ref > -1) datum_value = datum_value - branch_6d%ele(ix_ref)%i0

  case ('rad_int1.i1')
    datum_value = branch_6d%ele(ix_ele)%i1
    if (ix_ref > -1) datum_value = datum_value - branch_6d%ele(ix_ref)%i1

  case ('rad_int1.i2')
    datum_value = branch_6d%ele(ix_ele)%i2
    if (ix_ref > -1) datum_value = datum_value - branch_6d%ele(ix_ref)%i2

  case ('rad_int1.i3')
    datum_value = branch_6d%ele(ix_ele)%i3
    if (ix_ref > -1) datum_value = datum_value - branch_6d%ele(ix_ref)%i3

  case ('rad_int1.i2_e4')
    datum_value = branch_ri%ele(ix_ele)%lin_i2_e4
    if (ix_ref > -1) datum_value = datum_value - branch_ri%ele(ix_ref)%lin_i2_e4

  case ('rad_int1.i3_e7')
    datum_value = branch_ri%ele(ix_ele)%lin_i3_e7
    if (ix_ref > -1) datum_value = datum_value - branch_ri%ele(ix_ref)%lin_i3_e7

  case ('rad_int1.i4a')
    datum_value = branch_ri%ele(ix_ele)%i4a
    if (ix_ref > -1) datum_value = datum_value - branch_ri%ele(ix_ref)%i4a

  case ('rad_int1.i5a')
    datum_value = branch_ri%ele(ix_ele)%i5a
    if (ix_ref > -1) datum_value = datum_value - branch_ri%ele(ix_ref)%i5a

  case ('rad_int1.i5a_e6')
    datum_value = branch_ri%ele(ix_ele)%lin_i5a_e6
    if (ix_ref > -1) datum_value = datum_value - branch_ri%ele(ix_ref)%lin_i5a_e6

  case ('rad_int1.i4b')
    datum_value = branch_ri%ele(ix_ele)%i4b
    if (ix_ref > -1) datum_value = datum_value - branch_ri%ele(ix_ref)%i4b

  case ('rad_int1.i5b')
    datum_value = branch_ri%ele(ix_ele)%i5b
    if (ix_ref > -1) datum_value = datum_value - branch_ri%ele(ix_ref)%i5b

  case ('rad_int1.i5b_e6')
    datum_value = branch_ri%ele(ix_ele)%lin_i5b_e6
    if (ix_ref > -1) datum_value = datum_value - branch_ri%ele(ix_ref)%lin_i5b_e6

  case ('rad_int1.i6b')
    datum_value = branch_ri%ele(ix_ele)%i6b
    if (ix_ref > -1) datum_value = datum_value - branch_ri%ele(ix_ref)%i6b

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

  valid_value = .true.

!-----------

case ('ref_time')
    call tao_load_this_datum (branch%ele(:)%ref_time, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    
!-----------

case ('rel_floor.')

  select case (data_type)

  case ('rel_floor.x', 'rel_floor.y', 'rel_floor.z')

    if (ix_ref < 0) ele_ref => lat%branch(branch%ix_branch)%ele(0)

    call floor_angles_to_w_mat (-ele_ref%floor%theta, -ele_ref%floor%phi, -ele_ref%floor%psi, w0_mat)

    i = ix_start
    do 
      ele2 => branch%ele(i)
      vec3 = ele2%floor%r - ele_ref%floor%r
      vec3 = matmul (w0_mat, vec3)
      select case (data_type)
      case ('rel_floor.x')
        value_vec(i) = vec3(1)
      case ('rel_floor.y')
        value_vec(i) = vec3(2)
      case ('rel_floor.z')
        value_vec(i) = vec3(3)
      end select
      if (i == ix_ele) exit
      i = i + 1
      if (i > n_track) i = 0
    enddo

    call tao_load_this_datum (value_vec, null(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('rel_floor.theta', 'rel_floor.phi', 'rel_floor.psi')

    if (ix_ref < 0) ele_ref => lat%branch(branch%ix_branch)%ele(0)

    call floor_angles_to_w_mat (-ele_ref%floor%theta, -ele_ref%floor%phi, -ele_ref%floor%psi, w0_mat)

    i = ix_start
    do 
      ele2 => branch%ele(i)
      call floor_angles_to_w_mat (ele2%floor%theta, ele2%floor%phi, ele2%floor%psi, w_mat)
      w_mat = matmul (w0_mat, w_mat)
      call floor_w_mat_to_angles (w_mat, theta, phi, psi)

      select case (data_type)
      case ('rel_floor.theta')
        value_vec(i) = theta
      case ('rel_floor.phi')
        value_vec(i) = phi
      case ('rel_floor.psi')
        value_vec(i) = psi
      end select
      if (i == ix_ele) exit
      i = i + 1
      if (i > n_track) i = 0
    enddo

    call tao_load_this_datum (value_vec, null(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('s_position') 
  if (data_source == 'beam') goto 9000  ! Set error message and return
  if (ix_ref >= 0) then
    datum_value = ele%s - ele_ref%s
  else
    datum_value = ele%s 
  endif
  valid_value = .true.

!-----------

case ('sigma.')

  ! Looks for numbers: e.g. sigma.13
  i = index('123456', data_type(7:7))
  j = index('123456', data_type(8:8))
  if (i > 0 .and. j > 0) then
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%lat_sigma%mat(i,j), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(i,j), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif
    return
  endif

  select case (data_type)

  case ('sigma.x')
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%lat_sigma%mat(1,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      datum_value = sqrt(datum_value)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(1,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      datum_value = sqrt(datum_value)
    endif

  case ('sigma.px')  
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%lat_sigma%mat(2,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      datum_value = sqrt(datum_value)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(2,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      datum_value = sqrt(datum_value)
    endif

  case ('sigma.y')  
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%lat_sigma%mat(3,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      datum_value = sqrt(datum_value)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(3,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      datum_value = sqrt(datum_value)
    endif
    
  case ('sigma.py')  
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%lat_sigma%mat(4,4), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      datum_value = sqrt(datum_value)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(4,4), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      datum_value = sqrt(datum_value)
    endif
    
  case ('sigma.z')
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%lat_sigma%mat(5,5), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      datum_value = sqrt(datum_value)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(5,5), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      datum_value = sqrt(datum_value)
    endif

  case ('sigma.pz')  
    if (data_source == 'lat') then
      if (lat%param%geometry == closed$) then
        call tao_load_this_datum (tao_branch%lat_sigma%mat(6,6), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
        datum_value = sqrt(datum_value)
      else
        if (ix_ele == -1) ix_ele = branch%n_ele_track
        branch_ri => tao_lat%rad_int_by_ele_ri%branch(ix_branch)
        datum_value = branch_ri%ele(ix_ele)%lin_sig_E / ele%value(E_tot$)
        if (ix_ref > 0) datum_value = datum_value - branch_ri%ele(ix_ref)%lin_sig_E / ele_ref%value(E_tot$)
        valid_value = .true.
      endif
    else
      call tao_load_this_datum (bunch_params(:)%sigma(6,6), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      datum_value = sqrt(datum_value)
    endif
    
  case ('sigma.xy')  
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%lat_sigma%mat(1,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(1,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif

  case ('sigma.Lxy')  
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%lat_sigma%mat(1,4) - tao_branch%lat_sigma%mat(2,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(1,4) - bunch_params(:)%sigma(2,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('slip_factor_ptc.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  ptc_nf => tao_branch%ptc_normal_form

  if (.not. ptc_nf%valid_map) then
    if (.not. u%calc%one_turn_map) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE ONE_TURN_MAP_CALC IS NOT SET TO TRUE.', why_invalid)
    elseif (branch%param%geometry /= closed$) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE LATTICE GEOMETRY IS NOT CLOSED.', why_invalid)
    else
      call tao_set_invalid (datum, '?????', why_invalid)
    endif
    return
  endif

  if (.not. is_integer(data_type(17:), n)) then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  endif

  expo = 0
  expo(6) = n 

  datum_value = -real(ptc_nf%phase(3) .sub. expo) / branch%param%total_length
  valid_value = .true.

!-----------

case ('spin.')

  if (.not. bmad_com%spin_tracking_on) call tao_spin_tracking_turn_on()

  select case (data_type)

  case ('spin.x', 'spin.y', 'spin.z', 'spin.amp')
    do i = ix_start, ix_ele
      if (data_source == 'beam') then
        vec3 = bunch_params(i)%centroid%spin
      else
        vec3 = orbit(i)%spin
      endif

      select case (data_type)
      case ('spin.x');    value_vec(i) = vec3(1)
      case ('spin.y');    value_vec(i) = vec3(2)
      case ('spin.z');    value_vec(i) = vec3(3)
      case ('spin.amp');  value_vec(i) = norm2(vec3)
      end select
    enddo

    if (ix_ref > -1) then
      if (data_source == 'beam') then
        vec3 = bunch_params(ix_ref)%centroid%spin
      else
        vec3 = orbit(ix_ref)%spin
      endif

      select case (data_type)
      case ('spin.x');    value_vec(ix_ref) = vec3(1)
      case ('spin.y');    value_vec(ix_ref) = vec3(2)
      case ('spin.z');    value_vec(ix_ref) = vec3(3)
      case ('spin.amp');  value_vec(ix_ref) = norm2(vec3)
      end select
    endif

    call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    
  case ('spin.depolarization_rate', 'spin.polarization_rate', 'spin.polarization_limit')
    if (.not. tao_branch%spin%valid) call tao_spin_polarization_calc(branch, tao_branch, err_flag = err)
    valid_value = tao_branch%spin%valid

    if (.not. valid_value) then
      call tao_set_invalid (datum, 'ERROR IN SPIN POLARIZAITON CALC.', why_invalid, .false.)
      return
    endif
    select case (data_type)
    case ('spin.depolarization_rate')
      datum_value = tao_branch%spin%depol_rate
    case ('spin.polarization_rate')
      datum_value = tao_branch%spin%pol_rate_bks
    case ('spin.polarization_limit')
      datum_value = tao_branch%spin%pol_limit_dk
    end select

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('spin_dn_dpz.')

  if (data_source == 'beam') goto 9000  ! Set error message and return

  if (.not. bmad_com%spin_tracking_on) call tao_spin_tracking_turn_on()

  call tao_spin_polarization_calc(branch, tao_branch)

  valid_value = (tao_branch%spin_ele(ix_start)%valid .and. tao_branch%spin_ele(ix_ele)%valid)
  if (ix_ref > -1) valid_value = (valid_value .and. tao_branch%spin_ele(ix_ref)%valid)

  if (.not. valid_value) then
     call tao_set_invalid (datum, 'ERROR IN SPIN POLARIZAITON CALC.', why_invalid, .false.)
     return
  endif

  select case (data_type)
  case ('spin_dn_dpz.x')
    do i = ix_start, ix_ele
      value_vec(i) = tao_branch%spin_ele(i)%dn_dpz%vec(1)
    enddo
    value_vec(ix_ele) = tao_branch%spin_ele(ix_ele)%dn_dpz%vec(1)
  case ('spin_dn_dpz.y')
    do i = ix_start, ix_ele
      value_vec(i) = tao_branch%spin_ele(i)%dn_dpz%vec(2)
    enddo
    value_vec(ix_ele) = tao_branch%spin_ele(ix_ele)%dn_dpz%vec(2)
  case ('spin_dn_dpz.z')
    do i = ix_start, ix_ele
      value_vec(i) = tao_branch%spin_ele(i)%dn_dpz%vec(3)
    enddo
    value_vec(ix_ele) = tao_branch%spin_ele(ix_ele)%dn_dpz%vec(3)
  case ('spin_dn_dpz.amp')
    do i = ix_start, ix_ele
      value_vec(i) = norm2(tao_branch%spin_ele(i)%dn_dpz%vec)
    enddo
    value_vec(ix_ele) = norm2(tao_branch%spin_ele(ix_ele)%dn_dpz%vec)
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  end select

  call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)


!-----------

case ('spin_g_matrix.')

  call tao_spin_matrix_calc (datum, u, ele_ref, ele)
  valid_value = datum%spin_map%valid

  if (.not. valid_value) then
    call tao_set_invalid (datum, datum%why_invalid)
    return
  endif

  select case (data_type)
  case ('spin_g_matrix.11');  datum_value = datum%spin_map%mat8(7,1)
  case ('spin_g_matrix.12');  datum_value = datum%spin_map%mat8(7,2)
  case ('spin_g_matrix.13');  datum_value = datum%spin_map%mat8(7,3)
  case ('spin_g_matrix.14');  datum_value = datum%spin_map%mat8(7,4)
  case ('spin_g_matrix.15');  datum_value = datum%spin_map%mat8(7,5)
  case ('spin_g_matrix.16');  datum_value = datum%spin_map%mat8(7,6)
  case ('spin_g_matrix.21');  datum_value = datum%spin_map%mat8(8,1)
  case ('spin_g_matrix.22');  datum_value = datum%spin_map%mat8(8,2)
  case ('spin_g_matrix.23');  datum_value = datum%spin_map%mat8(8,3)
  case ('spin_g_matrix.24');  datum_value = datum%spin_map%mat8(8,4)
  case ('spin_g_matrix.25');  datum_value = datum%spin_map%mat8(8,5)
  case ('spin_g_matrix.26');  datum_value = datum%spin_map%mat8(8,6)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    valid_value = .false.
    return
  end select

!-----------

case ('spin_res.')

  if (data_source == 'beam') goto 9000  ! Set error message and return

  if (.not. bmad_com%spin_tracking_on) then
    call out_io (s_info$, r_name, 'Note: Turning on spin tracking (setting: bmad_com%spin_tracking_on = T)')
    bmad_com%spin_tracking_on = .true.
  endif

  call tao_spin_matrix_calc (datum, u, ele, ele)

  call spin_mat_to_eigen (datum%spin_map%map1%orb_mat, datum%spin_map%map1%spin_q, eval, evec, n0, n_eigen, err)
  if (err) then
    call tao_set_invalid (datum, 'ERROR CONVERTING SPIN/ORBIT 1-TURN MATRIX TO EIGEN VALUES.', why_invalid)
    return
  endif

  j = index('abc', data_type(10:10))
  if (j == 0) then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  endif

  call spin_quat_resonance_strengths(evec(2*j-1,:), datum%spin_map%map1%spin_q, xi_sum, xi_diff)

  select case (data_type(11:))
  case ('.sum');   datum_value = xi_sum
  case ('.diff');  datum_value = xi_diff
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  end select

  valid_value = .true.

!-----------

case ('spin_map_ptc.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  ptc_nf => tao_branch%ptc_normal_form

  if (.not. ptc_nf%valid_map) then
    if (.not. u%calc%one_turn_map) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE ONE_TURN_MAP_CALC IS NOT SET TO TRUE.', why_invalid)
    elseif (branch%param%geometry /= closed$) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE LATTICE GEOMETRY IS NOT CLOSED.', why_invalid)
    else
      call tao_set_invalid (datum, '?????', why_invalid)
    endif
    return
  endif

  if (.not. is_integer(data_type(14:), n)) then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  elseif (n > 999999 .or. n < 0) then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  endif

  do i = 1, 6
    expo(7-i) = mod(n, 10)
    n = n / 10
  enddo
  
  datum_value = real(ptc_nf%spin_tune .sub. expo)
  valid_value = .true.

!-----------

case ('spin_tune')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  if (.not. bmad_com%spin_tracking_on) call tao_spin_tracking_turn_on()

  select case (data_type)
  case ('spin_tune')
    datum_value = branch%param%spin_tune
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  end select

!-----------

case ('spin_tune_ptc.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  ptc_nf => tao_branch%ptc_normal_form

  if (.not. ptc_nf%valid_map) then
    if (.not. u%calc%one_turn_map) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE ONE_TURN_MAP_CALC IS NOT SET TO TRUE.', why_invalid)
    elseif (branch%param%geometry /= closed$) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE LATTICE GEOMETRY IS NOT CLOSED.', why_invalid)
    else
      call tao_set_invalid (datum, '?????', why_invalid)
    endif
    return
  endif

  if (.not. is_integer(data_type(15:), n)) then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  endif

  expo = 0
  expo(6) = n 

  datum_value = real(ptc_nf%spin_tune .sub. expo)
  valid_value = .true.

!-----------

case ('srdt.')
  select case(sub_data_type(1:6))
  case('h00111');  temp_cplx = tao_branch%srdt%h00111
  case('h00201');  temp_cplx = tao_branch%srdt%h00201
  case('h00220');  temp_cplx = tao_branch%srdt%h00220
  case('h00310');  temp_cplx = tao_branch%srdt%h00310
  case('h00400');  temp_cplx = tao_branch%srdt%h00400
  case('h10002');  temp_cplx = tao_branch%srdt%h10002
  case('h10020');  temp_cplx = tao_branch%srdt%h10020
  case('h10110');  temp_cplx = tao_branch%srdt%h10110
  case('h10200');  temp_cplx = tao_branch%srdt%h10200
  case('h11001');  temp_cplx = tao_branch%srdt%h11001
  case('h11110');  temp_cplx = tao_branch%srdt%h11110
  case('h11200');  temp_cplx = tao_branch%srdt%h11200
  case('h20001');  temp_cplx = tao_branch%srdt%h20001
  case('h20020');  temp_cplx = tao_branch%srdt%h20020
  case('h20110');  temp_cplx = tao_branch%srdt%h20110
  case('h20200');  temp_cplx = tao_branch%srdt%h20200
  case('h21000');  temp_cplx = tao_branch%srdt%h21000
  case('h22000');  temp_cplx = tao_branch%srdt%h22000
  case('h30000');  temp_cplx = tao_branch%srdt%h30000
  case('h31000');  temp_cplx = tao_branch%srdt%h31000
  case('h40000');  temp_cplx = tao_branch%srdt%h40000
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    term_found = .false.
    valid_value = .false.
    return
  end select

  valid_value = .true.
  select case (sub_data_type(8:8))
  case('r')
    datum_value = real(temp_cplx)
  case('i')
    datum_value = aimag(temp_cplx)
  case('a')
    datum_value = abs(temp_cplx)
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID (DATA_TYPE NOT ENDING IN .r, .i, or .a).', why_invalid, .true.)
    valid_value = .false.
  end select

!-----------

case ('time')
  if (data_source == 'beam') then
    call tao_load_this_datum (real(bunch_params%centroid%t, rp), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
  else
    call tao_load_this_datum (real(orbit(:)%t, rp), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  endif

!-----------

case ('tune.')

  if (data_source == 'beam') goto 9000  ! Set error message and return

  select case (data_type)
  case ('tune.a')
    datum_value = branch%ele(branch%n_ele_track)%a%phi
    valid_value = .true.

  case ('tune.b')
    datum_value = branch%ele(branch%n_ele_track)%b%phi
    valid_value = .true.

  case ('tune.z')
    call calc_z_tune (branch)
    datum_value = -branch%z%tune
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  end select

!-----------

case ('t.', 'tt.')
  if (branch%ix_branch /= 0) then
    call out_io (s_fatal$, r_name, 'TRANSFER MAP CALC NOT YET MODIFIED FOR BRANCHES.')
    return
  endif
  if (data_source == 'beam') goto 9000  ! Set error message and return

  expnt = 0
  i = tao_read_phase_space_index (sub_data_type, 1, .false.)
  do j = 2, 20
    if (sub_data_type(j:j) == ' ') exit
    k = tao_read_phase_space_index (sub_data_type, j, .false.); if (k == 0) exit
    expnt(k) = expnt(k) + 1
  enddo

  if (i == 0 .or. k == 0) then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  endif

  if (ix_ref < 0) ix_ref = 0

  ! Computation if there is no range

  if (ix_start == ix_ele) then
    if (s%com%ix_ref_taylor /= ix_ref .or. s%com%ix_ele_taylor /= ix_ele) then
      ix0 = s%com%ix_ele_taylor
      if (s%com%ix_ref_taylor == ix_ref .and. ix_ele > ix0) then
        call transfer_map_calc (lat, taylor_save, err, ix0, ix_ele, orbit(ix0), &
                                                  unit_start = .false., concat_if_possible = s%global%concatenate_maps)
      else
        call transfer_map_calc (lat, taylor_save, err, ix_ref, ix_ele, orbit(ix_ref), concat_if_possible = s%global%concatenate_maps)
      endif

      if (err) then
        call tao_set_invalid (datum, 'MAP TERM OVERFLOW', why_invalid)
        return
      endif

      s%com%ix_ref_taylor = ix_ref
      s%com%ix_ele_taylor = ix_ele
    endif
    datum_value = taylor_coef (taylor_save(i), expnt)
    valid_value = .true.

  ! Here if there is a range.
  else
    k = ix_start
    call transfer_map_calc (lat, taylor, err, ix_ref, k, orbit(ix_ref), concat_if_possible = s%global%concatenate_maps)
    if (err) then
      call tao_set_invalid (datum, 'MAP TERM OVERFLOW', why_invalid)
      return
    endif

    do
      value_vec(k) = taylor_coef (taylor(i), expnt)
      if (k == ix_ele) exit
      k_old = k
      k = k + 1
      if (k > branch%n_ele_track) k = 0
      call transfer_map_calc (lat, taylor, err, k_old, k, unit_start = .false., concat_if_possible = s%global%concatenate_maps)
      if (err) then
        call tao_set_invalid (datum, 'MAP TERM OVERFLOW', why_invalid)
        return
      endif
    enddo
    call tao_load_this_datum (value_vec, NULL(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  endif

!-----------

case ('unstable.')

  select case (data_type)

  case ('unstable.eigen', 'unstable.eigen.a', 'unstable.eigen.b', 'unstable.eigen.c')
    call transfer_matrix_calc (lat, mat6, vec0, 0, branch%n_ele_track, branch%ix_branch, one_turn = .true.)
    call mat_eigen (mat6, eval, evec, err)
    if (err) then
      call tao_set_invalid (datum, 'CANNOT COMPUTE EIGENVALUES FOR TRANSFER MATRIX', why_invalid)
      return
    endif

    select case (data_type)
    case ('unstable.eigen');    datum_value = maxval(abs(eval))
    case ('unstable.eigen.a');  datum_value = max(abs(eval(1)), abs(eval(2)))
    case ('unstable.eigen.b');  datum_value = max(abs(eval(3)), abs(eval(4)))
    case ('unstable.eigen.c');  datum_value = max(abs(eval(5)), abs(eval(6)))
    end select
    valid_value = .true.
    

  case ('unstable.orbit')

    if (data_source == 'beam') then
      ie0 = u%model_branch(ix_branch)%beam%ix_track_start

      if (datum%ele_name == '') then
        ie1 = u%model_branch(ix_branch)%beam%ix_track_end
      else
        ie1 = ix_ele
      endif

      if (ie0 == not_set$) then
        call tao_set_invalid (datum, 'NO TRACKING DONE IN BRANCH', why_invalid)
        return
      endif

      if (ie1 > ie0) then
        n = ie1 - ie0
      else
        n = branch%n_ele_track - ie0 + ie1
      endif

      datum_value = 0
      do j = 1, branch%n_ele_track
        jj = j + ie0
        if (jj > branch%n_ele_track) jj = jj - branch%n_ele_track
        datum_value = datum_value + (n - j + 1) * bunch_params(jj)%n_particle_lost_in_ele
        if (jj == ie1) exit
      enddo
      datum_value = datum_value / bunch_params(ie0)%n_particle_tot
      datum%ix_ele_merit = -1

    elseif (lat%param%geometry == open$) then
      if (datum%ele_name == '') ix_ele = branch%n_ele_track
      iz = tao_branch%track_state
      if (iz /= moving_forward$ .and. iz <= ix_ele) then
        datum_value = 1 + ix_ele - iz
        if (orbit(iz)%s < branch%ele(iz)%s) then
          orb => orbit(iz-1)
          datum%ix_ele_merit = iz - 1
          if (branch%ele(iz)%value(L$) /= 0) then
            ! Add s_rel/L 
            datum_value = datum_value + (branch%ele(iz)%s - orbit(iz)%s)/branch%ele(iz)%value(L$)
          endif
        else
          datum_value = datum_value - 0.5
          orb => orbit(iz)
          datum%ix_ele_merit = iz
        endif

        datum_value = datum_value + 0.5 * tanh(lat%param%unstable_factor)
      endif

    else   ! closed geometry
      if (tao_branch%track_state == moving_forward$) then
        datum_value = 0
      else
        datum_value = 1
      endif
      datum%ix_ele_merit = 0
    endif

    valid_value = .true.

  case ('unstable.ring', 'unstable.lattice')
    if (data_source == 'beam') goto 9000  ! Set error message and return

    if (data_type == 'unstable.ring') then
      call out_io (s_error$, r_name, '"unstable.ring" has been replaced by "unstable.lattice". Please change this in your input file.')
    endif

    if (lat%param%geometry == closed$ .and. tao_branch%track_state /= moving_forward$) then
      datum_value = 1
    else
      datum_value = lat%param%unstable_factor
      ! unstable_penalty is needed since at the metastable borderline the growth rate is zero.
      if (.not. lat%param%stable) datum_value = datum_value + s%global%unstable_penalty
    endif

    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('velocity', 'velocity.')

  if (tao_branch%track_state /= moving_forward$ .and. ix_ele > tao_branch%track_state) then
    valid_value = .false.
    call tao_set_invalid (datum, 'Particle lost.', why_invalid)
    return
  endif

  select case (data_type)

  case ('velocity')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('velocity.x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(2)*(1+bunch_params%centroid%vec(6))*bunch_params%centroid%beta, &
                        ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%vec(2)*(1+orbit(:)%vec(6))*orbit(:)%beta, &
                        ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('velocity.y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(4)*(1+bunch_params%centroid%vec(6))*bunch_params%centroid%beta, &
                        ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%vec(4)*(1+orbit(:)%vec(6))*orbit(:)%beta, &
                        ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('velocity.z')
    if (data_source == 'beam') then
      call tao_load_this_datum (sqrt(1 - (bunch_params%centroid%vec(2)*(1+bunch_params%centroid%vec(6)))**2 - (bunch_params%centroid%vec(4)*(1+bunch_params%centroid%vec(6)))**2)*(1+bunch_params%centroid%vec(6))*bunch_params%centroid%beta, &
                        ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (sqrt(1 - (orbit(:)%vec(2)*(1+orbit(:)%vec(6)))**2 - (orbit(:)%vec(4)*(1+orbit(:)%vec(6)))**2)*(1+orbit(:)%vec(6))*orbit(:)%beta, &
                        ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
  end select

!-----------

case ('wall.')
  if (data_source == 'beam') goto 9000  ! Set error message and return

  constraint = data_type(6:)
  zz0_pt = 0
  datum_value = 1e10   ! Something large
  found = .false.

  if (.not. allocated(s%building_wall%section)) then
    valid_value = .false.
    call tao_set_invalid (datum, 'No building wall sections defined.', why_invalid)
    return
  endif

  do i = 1, size(s%building_wall%section)
    section => s%building_wall%section(i)
    if (section%constraint /= constraint) cycle
    do ie = ix_start, ix_ele
      ele => branch%ele(ie)

      ! (zz, xx) is the local reference coordinate system at the element

      do is = 1, size(section%point)
        pt = tao_oreint_building_wall_pt(section%point(is))
        dz = pt%z - ele%floor%r(3); dx = pt%x - ele%floor%r(1)
        cos_theta = cos(ele%floor%theta); sin_theta = sin(ele%floor%theta)
        zz_pt =  dz * cos_theta + dx * sin_theta
        xx_pt = -dz * sin_theta + dx * cos_theta

        if (is == 1) then
          zz0_pt = zz_pt
          xx0_pt = xx_pt
          cycle
        endif


        if (pt%radius == 0 .or. zz_pt == 0 .or. zz0_pt == 0) then
          ! The perpendicular to the machine line intersects the segment if
          ! zz_pt and zz0_pt have a different signs.
          if (zz_pt * zz0_pt > 0) cycle
          xx_wall = (xx0_pt * zz_pt - xx_pt * zz0_pt) / (zz_pt - zz0_pt)

        else  ! Circular arc
          dz = pt%z_center - ele%floor%r(3); dx = pt%x_center - ele%floor%r(1)
          zz_center =  dz * cos_theta + dx * sin_theta
          xx_center = -dz * sin_theta + dx * cos_theta
          drad = pt%radius**2 -zz_center**2
          if (drad <= 0) cycle
          drad = sqrt(drad)
          ! There are two possible points at (0, xx_a) and (0, xx_b) where the perpendicular 
          ! to the machine line intersects the arc.
          xx_a = xx_center - drad
          xx_b = xx_center + drad
          dxx1 = xx_pt - xx0_pt
          dzz1 = zz_pt - zz0_pt
          ang_a = dzz1 * (xx_a - xx0_pt) + dxx1 * zz0_pt
          ang_b = dzz1 * (xx_b - xx0_pt) + dxx1 * zz0_pt
          ang_c = dzz1 * (xx_center - xx0_pt) - dxx1 * (zz_center - zz0_pt)
          ! A point is within the circular arc if it and the center point are on opposite
          ! sides of the chord from (zz0_pt, xx0_pt) to (zz_pt, xx_pt). 
          ! This assumes the arc is less than 180^deg.
          ! It should not be that both intersection points are within the arc.
          if (ang_a * ang_c < 0) then
            xx_wall = xx_a
          elseif (ang_b * ang_c < 0) then
            xx_wall = xx_b
          else
            cycle
          endif
        endif

        if (data_type =='wall.right_side')  xx_wall = -xx_wall
        datum_value = min(datum_value, xx_wall)
        valid_value = .true.

        zz0_pt = zz_pt
        xx0_pt = xx_pt
      enddo
    enddo

  enddo

  if (.not. valid_value) call tao_set_invalid (datum, 'No wall section found in the transverse plane of the evaluation point.', why_invalid)

!-----------

case ('wire.')  
  if (data_source == 'lat') goto 9001  ! Error message and return
  read (data_type(6:), '(a)') angle
  datum_value = tao_do_wire_scan (ele, angle, u%model_branch(branch%ix_branch)%ele(ix_ele)%beam)
  valid_value = .true.
  
case default
  call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
  return
end select

!-----------------------------------------------------------------------
! End stuff

if (datum%ix_ele_merit > -1) then
  datum%s = tao_datum_s_position(datum, branch%ele(datum%ix_ele_merit))
elseif (associated(ele)) then
  datum%s = tao_datum_s_position(datum, ele)
else
  datum%s = real_garbage$
endif

if (valid_value) datum%err_message_printed = .false.  ! Reset

return

!----------------------------------------------------------------------

7000 continue
call tao_set_invalid (datum, 'PARTICLE SPECIES TYPE NOT SET ??!! PLEASE SEEK HELP!', why_invalid)
return

9000 continue
call tao_set_invalid (datum, 'DATA_SOURCE = "beam" NOT VALID FOR THIS DATA_TYPE: ' // datum%data_type, why_invalid, .true.)
return

9001 continue
call tao_set_invalid (datum, 'DATA_SOURCE = "lat" NOT VALID FOR THIS DATA_TYPE: ' // datum%data_type, why_invalid, .true.)
return

9100 continue
call tao_set_invalid (datum, 'DATA_TYPE: ' // quote(datum%data_type) // ' NOT APPLICABLE TO A LATTICE BRANCH WITH AN OPEN GEOMETRY.', why_invalid, .true.)
return

9101 continue
call tao_set_invalid (datum, 'DATA_TYPE: ' // quote(datum%data_type) // ' NOT APPLICABLE TO A LATTICE BRANCH WITH A CLOSED GEOMETRY.', why_invalid, .true.)
return

end subroutine tao_evaluate_a_datum

