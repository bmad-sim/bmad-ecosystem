!+
! Subroutine pointer_to_attribute (ele, attrib_name, do_allocation, a_ptr, err_flag, err_print_flag, ix_attrib)
!
! Returns a pointer to an attribute of an element ele with attribute name attrib_name.
! Note: Use attribute_free to see if the attribute may be varied independently.
! Note: This routine will not work on bmad_com components. Rather use pointers_to_attribute.
! Note: Alternatively consider the routines:
!     pointers_to_attribute
!     set_ele_attribute
!     value_of_attribute
!
! Input:
!   ele             -- Ele_struct: After this routine finishes Ptr_attrib 
!                        will point to a variable within this element.
!   attrib_name     -- Character(40): Name of attribute. Must be uppercase.
!                       For example: "HKICK".
!   do_allocation   -- Logical: If True then do an allocation if needed.
!                       EG: The multipole An and Bn arrays need to be allocated
!                       before their use.
!   err_print_flag  -- Logical, optional: If present and False then suppress
!                       printing of an error message on error.
!
! Output:
!   a_ptr      -- all_pointer_struct: Pointer to the attribute. 
!     %r           -- pointer to real attribute. Nullified if error or attribute is not real.               
!     %i           -- pointer to integer attribute. Nullified if error or attribute is not integer.
!     %l           -- pointer to logical attribute. Nullified if error or attribute is not logical.               
!   err_flag   -- Logical: Set True if attribtute not found. False otherwise.
!   ix_attrib  -- Integer, optional: If applicable, this is the index to the 
!                     attribute in the ele%value(:), ele%control%var(:), ele%a_pole(:) or ele%b_pole(:) arrays.
!                     Set to 0 if not in any of these arrays.
!-

subroutine pointer_to_attribute (ele, attrib_name, do_allocation, a_ptr, err_flag, err_print_flag, ix_attrib)

use bmad_interface, except_dummy => pointer_to_attribute

implicit none

type (ele_struct), target :: ele
type (wake_lr_mode_struct), allocatable :: lr_mode(:)
type (cartesian_map_struct), pointer :: ct_map
type (cylindrical_map_struct), pointer :: cl_map
type (grid_field_struct), pointer :: g_field
type (taylor_field_struct), pointer :: t_field
type (all_pointer_struct) a_ptr
type (branch_struct), pointer :: branch
type (lat_struct), pointer :: lat

real(rp), pointer :: ptr_attrib, r(:,:,:)

integer, optional :: ix_attrib
integer ix_d, n, ios, n_lr_mode, ix_a, ix1, ix2, n_cc, n_coef, n_v, ix, iy, i, j, ivec(3)
integer expn(6)
integer lb0(3), ub0(3), lb(3), ub(3)
character(*) attrib_name
character(40) a_name
character(40) str
character(24) :: r_name = 'pointer_to_attribute'

logical err_flag, do_allocation, do_print, err, out_of_bounds
logical, optional :: err_print_flag

! init check

err_flag = .true.
out_of_bounds = .false.
branch => pointer_to_branch(ele)

nullify (a_ptr%r, a_ptr%i, a_ptr%l)

do_print = logic_option (.true., err_print_flag)
call str_upcase (a_name, attrib_name)
if (present(ix_attrib)) ix_attrib = 0

!--------------------
! If a controller with a defined list of variables

if (associated (ele%control)) then

  if (len(a_name) > 4) then
    if (a_name(1:4) == 'OLD_') then
      do i = 1, size(ele%control%var)
        if (ele%control%var(i)%name /= a_name(5:)) cycle
        a_ptr%r => ele%control%var(i)%old_value
        if (present(ix_attrib)) ix_attrib = old_control_var_offset$ + i
        err_flag = .false.
        return
      enddo
      goto 9000 ! Error message and return
    endif
  endif

  do i = 1, size(ele%control%var)
    if (ele%control%var(i)%name /= a_name) cycle
    a_ptr%r => ele%control%var(i)%value
    if (present(ix_attrib)) ix_attrib = var_offset$ + i
    err_flag = .false.
    return
  enddo

endif

! r_custom(...)

if (a_name(1:9) == 'R_CUSTOM(') THEN
  ix_d = index(a_name, ')')
  if (ix_d == 0) goto 9000 ! Error message and return
  str = a_name(10:ix_d-1) // ', -9999, 0, 0'  
  read (str, *, iostat = ios) ivec
  if (ios /= 0 .or. ivec(1) == -9999) goto 9000 ! ivec(1) must be present
  lb0 = 0; ub0 = 0
  if (associated(ele%r)) lb0 = lbound(ele%r)
  if (associated(ele%r)) ub0 = ubound(ele%r)
  if (ivec(2) == -9999) ivec(2) = 0
  if (ivec(3) == -9999) ivec(3) = 0

  lb = min(lb0, ivec)
  ub = max(ub0, ivec)
  if (associated(ele%r)) then
    if (.not. all(lb == lb0) .or. .not. all (ub == ub0)) then
      if (.not. do_allocation) goto 9110
      r => ele%r
      allocate(ele%r(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))
      ele%r = 0
      ele%r(lb0(1):ub0(1), lb0(2):ub0(2), lb0(3):ub0(3)) = r
      deallocate(r)
    endif
  else
    if (.not. do_allocation) goto 9110
    allocate(ele%r(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))
    ele%r = 0
  endif

  a_ptr%r => ele%r(ivec(1),ivec(2),ivec(3))
endif

!--------------------
! Check to see if the attribute is a long-range wake

if (a_name(1:3) == 'LR(' .or. a_name(1:13) == 'LR_WAKE%MODE(') then
  if (.not. associated (ele%wake)) then
    if (.not. do_allocation) goto 9100
    call init_wake (ele%wake, 0, 0, n)
  endif

  if (a_name(1:3) == 'LR(') then
    n = get_cross_index(a_name, 3, err, 1, 1000);  if (err) goto 9140
  else
    n = get_cross_index(a_name, 13, err, 1, 1000);  if (err) goto 9140
  endif

  n_lr_mode = size(ele%wake%lr%mode)
  if (n_lr_mode < n) then
    if (.not. do_allocation) goto 9100
    allocate (lr_mode(n_lr_mode))
    lr_mode = ele%wake%lr%mode
    deallocate (ele%wake%lr%mode)
    allocate (ele%wake%lr%mode(n))
    ele%wake%lr%mode = wake_lr_mode_struct ()
    ele%wake%lr%mode(1:n_lr_mode) = lr_mode
    deallocate (lr_mode)
  endif

  select case (a_name(ix_d+2:))
  case ('FREQ_IN');         a_ptr%r => ele%wake%lr%mode(n)%freq_in
  case ('FREQ');            a_ptr%r => ele%wake%lr%mode(n)%freq
  case ('R_OVER_Q');        a_ptr%r => ele%wake%lr%mode(n)%r_over_q
  case ('DAMP');            a_ptr%r => ele%wake%lr%mode(n)%damp
  case ('PHI');             a_ptr%r => ele%wake%lr%mode(n)%phi
  case ('POLAR_ANGLE');     a_ptr%r => ele%wake%lr%mode(n)%angle
  case ('POLARIZED');       a_ptr%l => ele%wake%lr%mode(n)%polarized
  case default; goto 9000
  end select    

  err_flag = .false.
  return
endif

!--------------------
! cartesian_map

if (a_name(1:14) == 'CARTESIAN_MAP(') then
  if (.not. associated(ele%cartesian_map)) goto 9130
  n_cc = get_cross_index(a_name, 14, err, 1, size(ele%cartesian_map))
  if (err) goto 9140
  ct_map => ele%cartesian_map(n_cc)

  select case (a_name)
  case ('%FIELD_SCALE');    a_ptr%r => ct_map%field_scale
  case ('%R0(1)');          a_ptr%r => ct_map%r0(1)
  case ('%R0(2)');          a_ptr%r => ct_map%r0(2)
  case ('%R0(3)');          a_ptr%r => ct_map%r0(3)
  case default;           goto 9000
  end select

  err_flag = .false.
  return

endif

!--------------------
! cylindrical_map

if (a_name(1:16) == 'CYLINDRICAL_MAP(') then
  if (.not. associated(ele%cylindrical_map)) goto 9130
  n_cc = get_cross_index(a_name, 16, err, 1, size(ele%cylindrical_map))
  if (err) goto 9140
  cl_map => ele%cylindrical_map(n_cc)

  select case (a_name)
  case ('%PHI0_FIELDMAP');  a_ptr%r => cl_map%phi0_fieldmap
  case ('%THETA0_AZIMUTH'); a_ptr%r => cl_map%theta0_azimuth
  case ('%FIELD_SCALE');    a_ptr%r => cl_map%field_scale
  case ('%DZ');             a_ptr%r => cl_map%dz
  case ('%R0(1)');          a_ptr%r => cl_map%r0(1)
  case ('%R0(2)');          a_ptr%r => cl_map%r0(2)
  case ('%R0(3)');          a_ptr%r => cl_map%r0(3)
  case default;           goto 9000
  end select

  err_flag = .false.
  return

endif

!--------------------
! grid_field

if (a_name(1:11) == 'GRID_FIELD(') then
  if (.not. associated(ele%grid_field)) goto 9130
  n_cc = get_cross_index(a_name, 11, err, 1, size(ele%grid_field))
  if (err) goto 9140
  g_field => ele%grid_field(n_cc)

  select case (a_name)
  case ('%PHI0_FIELDMAP');  a_ptr%r => g_field%phi0_fieldmap
  case ('%FIELD_SCALE');    a_ptr%r => g_field%field_scale
  case ('%R0(1)');          a_ptr%r => g_field%r0(1)
  case ('%R0(2)');          a_ptr%r => g_field%r0(2)
  case ('%R0(3)');          a_ptr%r => g_field%r0(3)
  case ('%DR(1)');          a_ptr%r => g_field%dr(1)
  case ('%DR(2)');          a_ptr%r => g_field%dr(2)
  case ('%DR(3)');          a_ptr%r => g_field%dr(3)
  case default;           goto 9000
  end select

  err_flag = .false.
  return

endif

!--------------------
! grid_field

if (a_name(1:13) == 'TAYLOR_FIELD(') then
  if (.not. associated(ele%taylor_field)) goto 9130
  n_cc = get_cross_index(a_name, 13, err, 1, size(ele%taylor_field))
  if (err) goto 9140
  t_field => ele%taylor_field(n_cc)

  select case (a_name)
  case ('%FIELD_SCALE');    a_ptr%r => t_field%field_scale
  case ('%DZ');             a_ptr%r => t_field%dz
  case ('%R0(1)');          a_ptr%r => t_field%r0(1)
  case ('%R0(2)');          a_ptr%r => t_field%r0(2)
  case ('%R0(3)');          a_ptr%r => t_field%r0(3)
  case default;           goto 9000
  end select

  err_flag = .false.
  return

endif

!--------------------
! wall3d section

if (a_name(1:12) == 'WALL%SECTION') then
  if (.not. associated(ele%wall3d)) goto 9210
  n_cc = get_cross_index(a_name, 13, err, 1, size(ele%wall3d(1)%section))
  if (err) goto 9200

  if (a_name == 'S') then
    if (n_cc == 1) goto 9210  ! must have s = 0
    a_ptr%r => ele%wall3d(1)%section(n_cc)%s
    err_flag = .false.
    return
  endif

  if (a_name(1:11) == 'WALL%DR_DS') then
    a_ptr%r => ele%wall3d(1)%section(n_cc)%dr_ds
    err_flag = .false.
    return
  endif

  if (a_name(1:1) == 'V') then
    n_v = get_cross_index(a_name, 2, err, 1, size(ele%wall3d(1)%section(n_cc)%v))
    if (err) goto 9200

    select case (a_name)
    case ('%X');        a_ptr%r => ele%wall3d(1)%section(n_cc)%v(n_v)%x
    case ('%Y');        a_ptr%r => ele%wall3d(1)%section(n_cc)%v(n_v)%y
    case ('%RADIUS_X'); a_ptr%r => ele%wall3d(1)%section(n_cc)%v(n_v)%radius_x
    case ('%RADIUS_Y'); a_ptr%r => ele%wall3d(1)%section(n_cc)%v(n_v)%radius_y
    case ('%TILT');     a_ptr%r => ele%wall3d(1)%section(n_cc)%v(n_v)%tilt
    case default;       goto 9200
    end select

    err_flag = .false.
    return
  endif

  goto 9200
endif

! Special cases

select case (a_name)
case ('BETA_A');          a_ptr%r => ele%a%beta
case ('BETA_B');          a_ptr%r => ele%b%beta
case ('ALPHA_A');         a_ptr%r => ele%a%alpha
case ('ALPHA_B');         a_ptr%r => ele%b%alpha
case ('GAMMA_A');         a_ptr%r => ele%a%gamma
case ('GAMMA_B');         a_ptr%r => ele%b%gamma
case ('PHI_A');           a_ptr%r => ele%a%phi
case ('PHI_B');           a_ptr%r => ele%b%phi
case ('ETA_A');           a_ptr%r => ele%a%eta
case ('ETA_B');           a_ptr%r => ele%b%eta
case ('ETA_X');           a_ptr%r => ele%x%eta
case ('ETA_Y');           a_ptr%r => ele%y%eta
case ('ETA_Z');           a_ptr%r => ele%z%eta
case ('ETAP_A');          a_ptr%r => ele%a%etap
case ('ETAP_B');          a_ptr%r => ele%b%etap
case ('ETAP_X');          a_ptr%r => ele%x%etap
case ('ETAP_Y');          a_ptr%r => ele%y%etap
case ('ETAP_Z');          a_ptr%r => ele%z%etap
case ('CMAT_11');         a_ptr%r => ele%c_mat(1,1)
case ('CMAT_12');         a_ptr%r => ele%c_mat(1,2)
case ('CMAT_21');         a_ptr%r => ele%c_mat(2,1)
case ('CMAT_22');         a_ptr%r => ele%c_mat(2,2)
case ('X_POSITION');      a_ptr%r => ele%floor%r(1)
case ('Y_POSITION');      a_ptr%r => ele%floor%r(2)
case ('Z_POSITION');      a_ptr%r => ele%floor%r(3)
case ('THETA_POSITION');  a_ptr%r => ele%floor%theta
case ('PHI_POSITION');    a_ptr%r => ele%floor%phi
case ('PSI_POSITION');    a_ptr%r => ele%floor%psi
case ('S');               a_ptr%r => ele%s
case ('LORD_STATUS');     a_ptr%i => ele%lord_status
case ('SLAVE_STATUS');    a_ptr%i => ele%slave_status
case ('REF_TIME');        a_ptr%r => ele%ref_time
case ('KEY');             a_ptr%i => ele%key
case ('N_SLAVE');         a_ptr%i => ele%n_slave
case ('N_LORD');          a_ptr%i => ele%n_lord
case ('LR_FREQ_SPREAD', 'LR_SELF_WAKE_ON', 'LR_WAKE%AMP_SCALE', 'LR_WAKE%TIME_SCALE', &
      'LR_WAKE%FREQ_SPREAD', 'LR_WAKE%SELF_WAKE_ON', &
      'SR_WAKE%SCALE_WITH_LENGTH', 'SR_WAKE%AMP_SCALE', 'SR_WAKE%Z_SCALE')
  if (.not. associated(ele%wake)) then
    if (.not. do_allocation) goto 9100
    call init_wake (ele%wake, 0, 0, 0, .true.)
  endif
  select case (a_name)
  case ('LR_SELF_WAKE_ON', 'LR_WAKE%SELF_WAKE_ON')
    a_ptr%l => ele%wake%lr%self_wake_on
  case ('LR_WAKE%AMP_SCALE')
    a_ptr%r => ele%wake%lr%amp_scale
  case ('LR_WAKE%TIME_SCALE')
    a_ptr%r => ele%wake%lr%time_scale
  case ('LR_FREQ_SPREAD', 'LR_WAKE%FREQ_SPREAD')
    a_ptr%r => ele%wake%lr%freq_spread
  case ('SR_WAKE%AMP_SCALE')
    a_ptr%r => ele%wake%sr%amp_scale
  case ('SR_WAKE%Z_SCALE')
    a_ptr%r => ele%wake%sr%z_scale
  case ('SR_WAKE%SCALE_WITH_LENGTH')
    a_ptr%l => ele%wake%sr%scale_with_length
  end select

case ('SPHERICAL_CURVATURE')
  a_ptr%r => ele%photon%surface%spherical_curvature
case ('ELLIPTICAL_CURVATURE_X')
  a_ptr%r => ele%photon%surface%elliptical_curvature(1)
case ('ELLIPTICAL_CURVATURE_Y')
  a_ptr%r => ele%photon%surface%elliptical_curvature(2)
case ('ELLIPTICAL_CURVATURE_Z')
  a_ptr%r => ele%photon%surface%elliptical_curvature(3)
end select

if (a_name(1:11) == 'CURVATURE_X' .and. a_name(13:14) == '_Y' .and. a_name(16:) == '') then
  ix = index('0123456789', a_name(12:12)) - 1
  iy = index('0123456789', a_name(15:15)) - 1
  if (ix == -1 .or. iy == -1) goto 9000 ! Error message and return
  if (ix > ubound(ele%photon%surface%curvature_xy, 1)) goto 9000 ! Error message and return
  if (iy > ubound(ele%photon%surface%curvature_xy, 2)) goto 9000 ! Error message and return
  a_ptr%r => ele%photon%surface%curvature_xy(ix,iy)
  err_flag = .false.
  return
endif

if (a_name(1:5) == "XMAT_") then
  if (len(a_name) >= 7) then
    ix1 = index('123456', a_name(6:6))
    ix2 = index('123456', a_name(7:7))
    if (ix1 > 0 .and. ix2 > 0) then
      a_ptr%r => ele%mat6(ix1,ix2)
      err_flag = .false.
      return
    endif
  endif
  goto 9000 ! Error message and return
endif

if (a_name(1:5) == 'VEC0_') then
  if (len(a_name) >= 6) then
    ix1 = index('123456', a_name(6:6))
    if (ix1 > 0) then
      a_ptr%r => ele%vec0(ix1)
      err_flag = .false.
      return
    endif
  endif
  goto 9000 ! Error message and return
endif

if (associated(a_ptr%r) .or. associated(a_ptr%l) .or. associated(a_ptr%i)) then
  err_flag = .false.
  return
endif

! Must be an indexed attribute

ix_a = attribute_index (ele, a_name)
if (present(ix_attrib)) ix_attrib = ix_a
if (ix_a < 1) goto 9000 ! Error message and return

! "Normal" indexed attribute

if (ix_a > 0 .and. ix_a <= num_ele_attrib$) then
  a_ptr%r => ele%value(ix_a)
  err_flag = .false.
  return
endif

! Custom attribute

if (ix_a > custom_attribute0$ .and. ix_a <= custom_attribute0$+custom_attribute_num$) then
  if (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) return

  n = ix_a - custom_attribute0$
  if (ele%key == def_parameter$) then
    lat => branch%lat
    if (.not. allocated(lat%custom)) then
      if (.not. do_allocation) return
      call re_allocate(lat%custom, custom_attribute_ubound_index(ele%key), .true., 0.0_rp)
    else
      if (size(lat%custom) < n .and. .not. do_allocation) return
      if (size(lat%custom) < n) call re_allocate(lat%custom, custom_attribute_ubound_index(ele%key), .true., 0.0_rp)
    endif
    a_ptr%r => lat%custom(n)
  else
    if (.not. associated(ele%custom)) then
      if (.not. do_allocation) return
      call re_associate(ele%custom, custom_attribute_ubound_index(ele%key), .true., 0.0_rp)
    else
      if (size(ele%custom) < n .and. .not. do_allocation) return
      if (size(ele%custom) < n) call re_associate(ele%custom, custom_attribute_ubound_index(ele%key), .true., 0.0_rp)
    endif
    a_ptr%r => ele%custom(n)
  endif
  err_flag = .false.
  return
endif

! Taylor term?

if (a_name(1:2) == 'TT') then
  n = index('123456', a_name(3:3))
  if (.not. associated(ele%taylor(1)%term)) then
    if (.not. do_allocation) return
    do i = 1, 6
      call init_taylor_series(ele%taylor(i), 0)
    enddo
  endif

  expn = 0
  do i = 4, len_trim(a_name)
    j = index('123456', a_name(i:i))
    expn(j) = expn(j) + 1
  enddo

  i = taylor_term_index(ele%taylor(n), expn, do_allocation)
  if (i /= 0) then
    a_ptr%r => ele%taylor(n)%term(i)%coef
    err_flag = .false.
  endif
  return
endif

! Magnetic Multipole

if (ix_a >= a0$ .and. ix_a <= b21$) then
  if (.not. associated(ele%a_pole)) then
    if (do_allocation) then
      call multipole_init (ele, magnetic$)
    else
      if (do_print) call out_io (s_error$, r_name, 'MULTIPOLE NOT ALLOCATED FOR ELEMENT: ' // ele%name)
      return
    endif
  endif

  if (ix_a >= b0$) then
    a_ptr%r => ele%b_pole(ix_a-b0$)
  else
    a_ptr%r => ele%a_pole(ix_a-a0$)
  endif

  err_flag = .false.
  return
endif

! Electric Multipole 

if (ix_a >= a0_elec$ .and. ix_a <= b21_elec$) then   
  if (.not. associated(ele%a_pole_elec)) then
    if (do_allocation) then
      call multipole_init (ele, electric$)
    else
      if (do_print) call out_io (s_error$, r_name, 'MULTIPOLE NOT ALLOCATED FOR ELEMENT: ' // ele%name)
      return
    endif
  endif

  if (ix_a >= b0_elec$) then
    a_ptr%r => ele%b_pole_elec(ix_a-b0_elec$)
  else
    a_ptr%r => ele%a_pole_elec(ix_a-a0_elec$)
  endif

  err_flag = .false.
  return
endif

! Special cases.

select case (a_name)
! attrib_type = is_real$
! attrib_type = is_logical$
case ('MATCH_END');                      a_ptr%r => ele%value(match_end$)
case ('MATCH_END_ORBIT');                a_ptr%r => ele%value(match_end_orbit$)
case ('FLEXIBLE');                       a_ptr%r => ele%value(flexible$)
case ('X_REF');                          a_ptr%r => ele%taylor(1)%ref
case ('PX_REF');                         a_ptr%r => ele%taylor(2)%ref
case ('Y_REF');                          a_ptr%r => ele%taylor(3)%ref
case ('PY_REF');                         a_ptr%r => ele%taylor(4)%ref
case ('Z_REF');                          a_ptr%r => ele%taylor(5)%ref
case ('PZ_REF');                         a_ptr%r => ele%taylor(6)%ref
case ('SYMPLECTIFY');                    a_ptr%l => ele%symplectify
case ('ABSOLUTE_TIME_TRACKING');         a_ptr%l => branch%lat%absolute_time_tracking
case ('TAYLOR_MAP_INCLUDES_OFFSETS');    a_ptr%l => ele%taylor_map_includes_offsets
case ('OFFSET_MOVES_APERTURE');          a_ptr%l => ele%offset_moves_aperture
case ('FIELD_MASTER');                   a_ptr%l => ele%field_master
case ('SCALE_MULTIPOLES');               a_ptr%l => ele%scale_multipoles
case ('MULTIPOLES_ON');                  a_ptr%l => ele%multipoles_on
case ('IS_ON');                          a_ptr%l => ele%is_on
case ('HIGH_ENERGY_SPACE_CHARGE_ON');    a_ptr%l => branch%param%high_energy_space_charge_on
!  attrib_type = is_integer$
case ('N_SLICE');                        a_ptr%r => ele%value(n_slice$)
case ('MULTIPASS_REF_ENERGY');           a_ptr%r => ele%value(multipass_ref_energy$)
case ('DIRECTION');                      a_ptr%r => ele%value(direction$)
case ('N_CELL');                         a_ptr%r => ele%value(n_cell$)
case ('IX_TO_BRANCH');                   a_ptr%r => ele%value(ix_to_branch$)
case ('IX_TO_ELEMENT');                  a_ptr%r => ele%value(ix_to_element$)
case ('NUM_STEPS');                      a_ptr%r => ele%value(num_steps$)
case ('INTEGRATOR_ORDER');               a_ptr%r => ele%value(integrator_order$)
!  attrib_type = is_switch$
case ('APERTURE_AT');                    a_ptr%i => ele%aperture_at
case ('APERTURE_TYPE');                  a_ptr%i => ele%aperture_type
case ('COUPLER_AT');                     a_ptr%r => ele%value(coupler_at$)
case ('CSR_METHOD');                     a_ptr%i => ele%csr_method
case ('DEFAULT_TRACKING_SPECIES');       a_ptr%i => branch%param%default_tracking_species
case ('FIELD_CALC');                     a_ptr%i => ele%field_calc
case ('FRINGE_TYPE');                    a_ptr%r => ele%value(fringe_type$)
case ('GEOMETRY');                       a_ptr%r => ele%value(geometry$)
case ('LIVE_BRANCH');                    a_ptr%r => ele%value(live_branch$)
case ('FRINGE_AT');                      a_ptr%r => ele%value(fringe_at$)
case ('HIGHER_ORDER_FRINGE_TYPE');       a_ptr%r => ele%value(higher_order_fringe_type$)
case ('MAT6_CALC_METHOD');               a_ptr%i => ele%mat6_calc_method
case ('MODE');                           a_ptr%r => ele%value(mode$)
case ('ORIGIN_ELE_REF_PT');              a_ptr%r => ele%value(origin_ele_ref_pt$)
case ('PTC_FIELD_GEOMETRY');             a_ptr%r => ele%value(ptc_field_geometry$)
case ('PTC_INTEGRATION_TYPE');           a_ptr%i => ele%ptc_integration_type
case ('PTC_FRINGE_GEOMETRY');            a_ptr%r => ele%value(ptc_fringe_geometry$)
case ('REF_ORBIT_FOLLOWS');              a_ptr%r => ele%value(ref_orbit_follows$)
case ('REF_COORDINATES');                a_ptr%r => ele%value(ref_coordinates$)
case ('SPACE_CHARGE_METHOD');            a_ptr%i => ele%space_charge_method
case ('SPIN_TRACKING_METHOD');           a_ptr%i => ele%spin_tracking_method
case ('TRACKING_METHOD');                a_ptr%i => ele%tracking_method
case ('REF_SPECIES');                    a_ptr%i => ele%ref_species
case ('PARTICLE')
  if (ele%key == def_line$) then
    a_ptr%i => ele%ref_species
  else
    a_ptr%i => branch%param%particle
  endif

! No corresponding attribute in element.
case ('TAYLOR_ORDER')
case ('PTC_EXACT_MODEL')
case ('PTC_EXACT_MISALIGN')
case ('PTC_MAX_FRINGE_ORDER')
case ('UPSTREAM_ELE_DIR')
case ('DOWNSTREAM_ELE_DIR')
end select

if (associated(a_ptr%r) .or. associated(a_ptr%i) .or. associated(a_ptr%l)) then
  err_flag = .false.
else
  goto 9000
endif

return

!----------------------------------------
! Error message and return

9000 continue
if (do_print) call out_io (s_error$, r_name, &
          'INVALID ATTRIBUTE: ' // a_name, 'FOR THIS ELEMENT: ' // ele%name)
return

!----------------------------------------
9100 continue
if (do_print) call out_io (s_error$, r_name, &
                 'WAKE ATTRIBUTE NOT ALLOCATED: ' // a_name, &
                 'FOR THIS ELEMENT: ' // ele%name)
return

!----------------------------------------
9110 continue
if (do_print) call out_io (s_error$, r_name, &
                 'R_CUSTOM ATTRIBUTE NOT ALLOCATED: ' // a_name, &
                 'FOR THIS ELEMENT: ' // ele%name)
return

!----------------------------------------
9130 continue
if (do_print) then
  if (out_of_bounds) then
    call out_io (s_error$, r_name, &
        'INDEX OUT OF BOUNDS IN ATTRIBUTE: ' // a_name, &
        'FOR THIS ELEMENT: ' // ele%name)
  else
    call out_io (s_error$, r_name, &
        'MALFORMED ATTRIBUTE: ' // a_name, &
        'FOR THIS ELEMENT: ' // ele%name)
  endif
endif
return

!----------------------------------------
9140 continue
if (do_print) call out_io (s_error$, r_name, &
                 '(EM) FIELD ATTRIBUTE NOT ALLOCATED: ' // a_name, &
                 'FOR THIS ELEMENT: ' // ele%name)
return

!----------------------------------------
9200 continue
if (do_print) then
  if (out_of_bounds) then
    call out_io (s_error$, r_name, &
        'INDEX OUT OF BOUNDS IN ATTRIBUTE: ' // a_name, &
        'FOR THIS ELEMENT: ' // ele%name)
  else
    call out_io (s_error$, r_name, &
        'MALFORMED ATTRIBUTE: ' // a_name, &
        'FOR THIS ELEMENT: ' // ele%name)
  endif
endif
return

!----------------------------------------
9210 continue
if (do_print) call out_io (s_error$, r_name, &
        'CROSS-SECTION NOT DEFINED SO CANNOT SET ATTRIBUTE: ' // a_name, &
        'FOR THIS ELEMENT: ' // ele%name)
return

!---------------------------------------------------------------
contains

!+
! Function reads number of the form "...(num)" and checks to
! see if num is between n_min and n_max.
! Function also chops "...(num)" from name.
!-

function get_cross_index(name, ix_name, err, n_min, n_max) result (ixc)

character(*) name

integer ix_name, n_min, n_max, ixc, ios

logical err

!

err = .true.

if (name(ix_name:ix_name) /= '(') return
name = name(ix_name+1:)

ix = index(name, ')')
if (ix < 2) return

read (name(1:ix-1), *, iostat = ios) ixc
if (ios /= 0 .or. name(1:ix-1) == '') return
name = name(ix+1:)

if (ixc < n_min .or. ixc > n_max) then
  out_of_bounds = .true.
  return
endif

err = .false.

end function

end subroutine
