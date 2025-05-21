!+
! Subroutine pointer_to_attribute (ele, attrib_name, do_allocation, a_ptr, err_flag, err_print_flag, ix_attrib, do_unlink)
!
! Returns a pointer to an attribute of an element ele with attribute name attrib_name.
! Note: Use attribute_free to see if the attribute may be varied independently.
! Note: This routine will not work on bmad_com components. Rather use pointers_to_attribute.
!
! Note: To save memory, ele%cartesian_map (and other field maps), can point to the same memory location as the 
! Cartesian maps of other elements. This linkage is not desired if the attribute to be pointed to is varied. 
! In this case, the do_unlink argumnet should be set to True.
!
! Note: Alternatively consider the routines:
!     pointers_to_attribute
!     set_ele_attribute
!     value_of_attribute
!
! Input:
!   ele             -- ele_struct: After this routine finishes Ptr_attrib 
!                        will point to a variable within this element.
!   attrib_name     -- character(40): Name of attribute. Must be uppercase.
!                       For example: "HKICK".
!   do_allocation   -- logical: If True then do an allocation if needed.
!                       EG: The multipole An and Bn arrays need to be allocated
!                       before their use.
!   err_print_flag  -- logical, optional: If present and False then suppress
!                       printing of an error message on error.
!   do_unlink       -- logical, optional: Default False. If True and applicable, unlink the structure containing the attribute.
!                       See above for details.
!
! Output:
!   a_ptr      -- all_pointer_struct: Pointer to the attribute. 
!     %r           -- pointer to real attribute. Nullified if error or attribute is not real.               
!     %i           -- pointer to integer attribute. Nullified if error or attribute is not integer.
!     %l           -- pointer to logical attribute. Nullified if error or attribute is not logical.               
!   err_flag   -- logical: Set True if attribtute not found. False otherwise.
!   ix_attrib  -- integer, optional: If applicable, this is the index to the 
!                     attribute in the ele%value(:), ele%control%var(:), ele%a_pole(:) or ele%b_pole(:) arrays.
!                     Set to 0 if not in any of these arrays.
!-

subroutine pointer_to_attribute (ele, attrib_name, do_allocation, a_ptr, err_flag, err_print_flag, ix_attrib, do_unlink)

use bmad_interface, except_dummy => pointer_to_attribute

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: slave
type (wake_lr_mode_struct), allocatable :: lr_mode(:)
type (cartesian_map_struct), pointer :: ct_map
type (cartesian_map_term_struct), pointer :: ct_ptr
type (cartesian_map_term1_struct), pointer :: ct_term
type (cylindrical_map_struct), pointer :: cl_map
type (grid_field_struct), pointer :: g_field
type (gen_grad_map_struct), pointer :: gg_map
type (all_pointer_struct) a_ptr
type (branch_struct), pointer :: branch
type (lat_struct), pointer :: lat
type (control_struct), pointer :: ctl
type (control_ramp1_struct), pointer :: rmp
type (taylor_struct), pointer :: tlr

real(rp), pointer :: ptr_attrib, r(:,:,:)

integer, optional :: ix_attrib
integer ix_d, n, ios, n_lr_mode, ix_a, ix1, ix2, n_cc, n_coef, n_v, ix, iy, i, j, ivec(3), ixs, i0, nt
integer expn(6)
integer lb0(3), ub0(3), lb(3), ub(3)
character(*) attrib_name
character(40) a_name
character(40) str
character(24) :: r_name = 'pointer_to_attribute'

logical err_flag, do_allocation, do_print, err, out_of_bounds
logical, optional :: err_print_flag, do_unlink

! init check

err_flag = .true.
out_of_bounds = .false.
branch => pointer_to_branch(ele)

nullify (a_ptr%r, a_ptr%i, a_ptr%l, a_ptr%r1, a_ptr%i1)

do_print = logic_option (.true., err_print_flag)
call str_upcase (a_name, attrib_name)
if (present(ix_attrib)) ix_attrib = 0

!

if (ele%key == def_ptc_com$) then
  select case (a_name)
  case ('VERTICAL_KICK');                             a_ptr%r => ptc_com%vertical_kick
  case ('CUT_FACTOR');                                a_ptr%r => ptc_com%cut_factor
  case ('TRANSLATE_PATCH_DRIFT_TIME');                a_ptr%l => ptc_com%translate_patch_drift_time
  case ('PRINT_INFO_MESSAGES');                       a_ptr%l => ptc_com%print_info_messages
  case ('USE_ORIENTATION_PATCHES');                   a_ptr%l => ptc_com%use_orientation_patches
  case ('OLD_INTEGRATOR');                            a_ptr%i => ptc_com%old_integrator
  case ('EXACT_MODEL', 'PTC_EXACT_MODEL');            a_ptr%l => ptc_com%exact_model
  case ('EXACT_MISALIGN', 'PTC_EXACT_MISALIGN');      a_ptr%l => ptc_com%exact_misalign
  case ('MAX_FRINGE_ORDER', 'PTC_MAX_FRINGE_ORDER');  a_ptr%i => ptc_com%max_fringe_order
  case default; goto 9000
  end select
  err_flag = .false.
  return
endif

!--------------------
! If a controller with a defined list of variables
! Note: ele%control or ele%control%var may not be allocated during parsing.

if ((ele%key == ramper$ .or. ele%key == overlay$ .or. ele%key == group$) .and. associated(ele%control)) then

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

  if (a_name(1:6) == 'SLAVE(') then
    n = get_this_index(a_name, 6, err, 1, ele%n_slave); if (err) goto 9130
    if (a_name(1:1) /= '%') goto 9000

    if (ele%key == ramper$) then
      rmp => ele%control%ramp(n)

      if (allocated(rmp%y_knot)) then
        if (a_name(1:8) /= '%Y_KNOT(') goto 9300
        n = get_this_index(a_name, 8, err, 1, size(rmp%y_knot)); if (err) goto 9130
        a_ptr%r => rmp%y_knot(n)
        err_flag = .false.
        return
      endif

      a_name = a_name(2:)  ! Remove '%'
      do i = 1, size(rmp%stack)
        if (upcase(rmp%stack(i)%name) /= a_name) cycle
        a_ptr%r => rmp%stack(i)%value
        err_flag = .false.
        return
      enddo

    else
      slave => pointer_to_slave(ele, n, ctl)

      if (allocated(ctl%y_knot)) then
        if (a_name(1:8) /= '%Y_KNOT(') goto 9300
        n = get_this_index(a_name, 8, err, 1, size(ctl%y_knot)); if (err) goto 9130
        a_ptr%r => ctl%y_knot(n)
        err_flag = .false.
        return
      endif

      a_name = a_name(2:)  ! Remove '%'
      do i = 1, size(ctl%stack)
        if (upcase(ctl%stack(i)%name) /= a_name) cycle
        a_ptr%r => ctl%stack(i)%value
        err_flag = .false.
        return
      enddo
    endif

    goto 9310
  endif

  if (a_name(1:7) == 'X_KNOT(') then
    if (.not. allocated(ele%control%x_knot)) goto 9320
    n = get_this_index(a_name, 7, err, 1, size(ele%control%x_knot)); if (err) goto 9130
    a_ptr%r => ele%control%x_knot(n)
    err_flag = .false.
    return
  endif
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
    call init_wake (ele%wake, 0, 0, 0, n)
  endif

  if (a_name(1:3) == 'LR(') then
    n = get_this_index(a_name, 3, err, 1, 1000);  if (err) goto 9140
  else
    n = get_this_index(a_name, 13, err, 1, 1000);  if (err) goto 9140
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

  select case (a_name)
  case ('%FREQ_IN');         a_ptr%r => ele%wake%lr%mode(n)%freq_in
  case ('%FREQ');            a_ptr%r => ele%wake%lr%mode(n)%freq
  case ('%R_OVER_Q');        a_ptr%r => ele%wake%lr%mode(n)%r_over_q
  case ('%DAMP');            a_ptr%r => ele%wake%lr%mode(n)%damp
  case ('%PHI');             a_ptr%r => ele%wake%lr%mode(n)%phi
  case ('%POLAR_ANGLE');     a_ptr%r => ele%wake%lr%mode(n)%angle
  case ('%POLARIZED');       a_ptr%l => ele%wake%lr%mode(n)%polarized
  case default; goto 9000
  end select    

  err_flag = .false.
  return
endif

!--------------------
! cartesian_map

if (a_name(1:14) == 'CARTESIAN_MAP(') then
  if (.not. associated(ele%cartesian_map)) goto 9130
  n_cc = get_this_index(a_name, 14, err, 1, size(ele%cartesian_map))
  if (err) goto 9140
  ct_map => ele%cartesian_map(n_cc)
  if (.not. associated(ct_map%ptr)) return

  ct_ptr => ct_map%ptr
  if (logic_option(.false., do_unlink) .and. ct_map%ptr%n_link > 1) then
    ct_map%ptr%n_link = ct_map%ptr%n_link - 1
    allocate(ct_map%ptr)
    ct_map%ptr = ct_ptr
    ct_map%ptr%n_link = 1
  endif

  if (a_name(1:3) == '%T(' .or. a_name(1:6) == '%TERM(') then
    nt = get_this_index(a_name, index(a_name, '('), err, 1, size(ct_map%ptr%term))
    if (err) goto 9140
    ct_term => ct_map%ptr%term(nt)
    select case (a_name)
    case ('%A');            a_ptr%r => ct_term%coef
    case ('%KX', '%K_X');   a_ptr%r => ct_term%kx
    case ('%KY', '%K_Y');   a_ptr%r => ct_term%ky
    case ('%KZ', '%K_Z');   a_ptr%r => ct_term%kz
    case ('%X0', '%X_0');   a_ptr%r => ct_term%x0
    case ('%Y0', '%Y_0');   a_ptr%r => ct_term%y0
    case ('%PHI_Z');        a_ptr%r => ct_term%phi_z
    case default;           goto 9000
    end select

  else
    select case (a_name)
    case ('%FIELD_SCALE');      a_ptr%r => ct_map%field_scale
    case ('%R0');               a_ptr%r1 => ct_map%r0
    case ('%R0(1)');            a_ptr%r => ct_map%r0(1)
    case ('%R0(2)');            a_ptr%r => ct_map%r0(2)
    case ('%R0(3)');            a_ptr%r => ct_map%r0(3)
    case ('%MASTER_PARAMETER'); a_ptr%i => ct_map%master_parameter
    case default;           goto 9000
    end select
  endif

  err_flag = .false.
  return

endif

!--------------------
! cylindrical_map

if (a_name(1:16) == 'CYLINDRICAL_MAP(') then
  if (.not. associated(ele%cylindrical_map)) goto 9130
  n_cc = get_this_index(a_name, 16, err, 1, size(ele%cylindrical_map))
  if (err) goto 9140
  cl_map => ele%cylindrical_map(n_cc)

  select case (a_name)
  case ('%PHI0_FIELDMAP');    a_ptr%r  => cl_map%phi0_fieldmap
  case ('%THETA0_AZIMUTH');   a_ptr%r  => cl_map%theta0_azimuth
  case ('%FIELD_SCALE');      a_ptr%r  => cl_map%field_scale
  case ('%DZ');               a_ptr%r  => cl_map%dz
  case ('%R0');               a_ptr%r1 => cl_map%r0
  case ('%R0(1)');            a_ptr%r  => cl_map%r0(1)
  case ('%R0(2)');            a_ptr%r  => cl_map%r0(2)
  case ('%R0(3)');            a_ptr%r  => cl_map%r0(3)
  case ('%MASTER_PARAMETER'); a_ptr%i  => cl_map%master_parameter
  case default;           goto 9000
  end select

  err_flag = .false.
  return

endif

!--------------------
! gen_grad_map

if (a_name(1:13) == 'GEN_GRAD_MAP(') then
  if (.not. associated(ele%gen_grad_map)) goto 9130
  n_cc = get_this_index(a_name, 13, err, 1, size(ele%gen_grad_map))
  if (err) goto 9140
  gg_map => ele%gen_grad_map(n_cc)

  select case (a_name)
  case ('%FIELD_SCALE');      a_ptr%r => gg_map%field_scale
  case ('%DZ');               a_ptr%r => gg_map%dz
  case ('%R0(1)');            a_ptr%r => gg_map%r0(1)
  case ('%R0(2)');            a_ptr%r => gg_map%r0(2)
  case ('%R0(3)');            a_ptr%r => gg_map%r0(3)
  case ('%MASTER_PARAMETER'); a_ptr%i => gg_map%master_parameter
  case default;           goto 9000
  end select

  err_flag = .false.
  return

endif

!--------------------
! grid_field

if (a_name(1:11) == 'GRID_FIELD(') then
  if (.not. associated(ele%grid_field)) goto 9130
  n_cc = get_this_index(a_name, 11, err, 1, size(ele%grid_field))
  if (err) goto 9140
  g_field => ele%grid_field(n_cc)

  select case (a_name)
  case ('%INTERPOLATION_ORDER');  a_ptr%i => g_field%interpolation_order
  case ('%HARMONIC');             a_ptr%i => g_field%harmonic
  case ('%GEOMETRY');             a_ptr%i => g_field%geometry
  case ('%ELE_ANCHOR_PT');        a_ptr%i => g_field%ele_anchor_pt
  case ('%PHI0_FIELDMAP');        a_ptr%r => g_field%phi0_fieldmap
  case ('%FIELD_SCALE');          a_ptr%r => g_field%field_scale
  case ('%R0(1)');                a_ptr%r => g_field%r0(1)
  case ('%R0(2)');                a_ptr%r => g_field%r0(2)
  case ('%R0(3)');                a_ptr%r => g_field%r0(3)
  case ('%DR(1)');                a_ptr%r => g_field%dr(1)
  case ('%DR(2)');                a_ptr%r => g_field%dr(2)
  case ('%DR(3)');                a_ptr%r => g_field%dr(3)
  case ('%MASTER_PARAMETER');     a_ptr%i => g_field%master_parameter
  case default;                   goto 9000
  end select

  err_flag = .false.
  return

endif

!--------------------
! wall3d section

if (a_name(1:12) == 'WALL%SECTION') then
  if (.not. associated(ele%wall3d)) goto 9210
  n_cc = get_this_index(a_name, 13, err, 1, size(ele%wall3d(1)%section))
  if (err) goto 9130

  if (a_name == '%S') then
    if (n_cc == 1) goto 9210  ! must have s = 0
    a_ptr%r => ele%wall3d(1)%section(n_cc)%s
    err_flag = .false.
    return
  endif

  if (a_name == '%DR_DS') then
    a_ptr%r => ele%wall3d(1)%section(n_cc)%dr_ds
    err_flag = .false.
    return
  endif

  if (a_name(1:2) == '%V') then
    n_v = get_this_index(a_name, 2, err, 1, size(ele%wall3d(1)%section(n_cc)%v))
    if (err) goto 9130

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

  goto 9130
endif

!---------------
! AC_kicker

if (a_name(1:12) == 'AMP_VS_TIME(') then
  if (.not. associated(ele%ac_kick)) goto 9400
  if (.not. allocated(ele%ac_kick%amp_vs_time)) goto 9410
  n = get_this_index(a_name, 12, err, 1, size(ele%ac_kick%amp_vs_time)); if (err) goto 9420

  select case (a_name)
  case ('%TIME');   a_ptr%r => ele%ac_kick%amp_vs_time(n)%time
  case ('%AMP');    a_ptr%r => ele%ac_kick%amp_vs_time(n)%amp
  case default;     goto 9430
  end select

  err_flag = .false.
  return
endif

if (a_name(1:12) == 'FREQUENCIES(') then
  if (.not. associated(ele%ac_kick)) goto 9400
  if (.not. allocated(ele%ac_kick%frequency)) goto 9450
  n = get_this_index(a_name, 12, err, 1, size(ele%ac_kick%frequency)); if (err) goto 9460

  select case (a_name)
  case ('%FREQ'); a_ptr%r => ele%ac_kick%frequency(n)%f
  case ('%AMP');  a_ptr%r => ele%ac_kick%frequency(n)%amp
  case ('%PHI');  a_ptr%r => ele%ac_kick%frequency(n)%phi
  case default;   goto 9470
  end select

  err_flag = .false.
  return
endif

!---------------
! Special cases

if (ele%key == rbend$) then   ! Note: Rbend elements only exist during lattice parsing.
  select case (a_name)
  case ('L');       a_ptr%r => ele%value(l_chord$)
  case ('L_CHORD'); a_ptr%r => ele%value(l_chord$)
  case ('L_ARC');   a_ptr%r => ele%value(l$)
  end select

  if (associated(a_ptr%r)) then
    err_flag = .false.
    return
  endif
endif

select case (a_name)
case ('G_ERR');           a_ptr%r => ele%value(dg$)        ! Old name
case ('B_FIELD_ERR');     a_ptr%r => ele%value(db_field$)  ! Old name
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

case ('DETA_A_DS');       a_ptr%r => ele%a%deta_ds
case ('DETA_B_DS');       a_ptr%r => ele%b%deta_ds
case ('DETA_X_DS');       a_ptr%r => ele%x%deta_ds
case ('DETA_Y_DS');       a_ptr%r => ele%y%deta_ds
case ('DETA_Z_DS');       a_ptr%r => ele%z%deta_ds

case ('DETA_DPZ_A');      a_ptr%r => ele%a%deta_dpz
case ('DETA_DPZ_B');      a_ptr%r => ele%b%deta_dpz
case ('DETA_DPZ_X');      a_ptr%r => ele%x%deta_dpz
case ('DETA_DPZ_Y');      a_ptr%r => ele%y%deta_dpz
case ('DETA_DPZ_Z');      a_ptr%r => ele%z%deta_dpz

case ('DETAP_DPZ_A');      a_ptr%r => ele%a%detap_dpz
case ('DETAP_DPZ_B');      a_ptr%r => ele%b%detap_dpz
case ('DETAP_DPZ_X');      a_ptr%r => ele%x%detap_dpz
case ('DETAP_DPZ_Y');      a_ptr%r => ele%y%detap_dpz
case ('DETAP_DPZ_Z');      a_ptr%r => ele%z%detap_dpz

case ('DBETA_DPZ_A');     a_ptr%r => ele%a%dbeta_dpz
case ('DBETA_DPZ_B');     a_ptr%r => ele%b%dbeta_dpz
case ('DALPHA_DPZ_A');    a_ptr%r => ele%a%dalpha_dpz
case ('DALPHA_DPZ_B');    a_ptr%r => ele%b%dalpha_dpz
case ('CMAT_11');         a_ptr%r => ele%c_mat(1,1)
case ('CMAT_12');         a_ptr%r => ele%c_mat(1,2)
case ('CMAT_21');         a_ptr%r => ele%c_mat(2,1)
case ('CMAT_22');         a_ptr%r => ele%c_mat(2,2)
case ('MODE_FLIP');       a_ptr%l => ele%mode_flip
case ('X_POSITION');      a_ptr%r => ele%floor%r(1)
case ('Y_POSITION');      a_ptr%r => ele%floor%r(2)
case ('Z_POSITION');      a_ptr%r => ele%floor%r(3)
case ('THETA_POSITION');  a_ptr%r => ele%floor%theta
case ('PHI_POSITION');    a_ptr%r => ele%floor%phi
case ('PSI_POSITION');    a_ptr%r => ele%floor%psi
case ('S');               a_ptr%r => ele%s
case ('LORD_STATUS');     a_ptr%i => ele%lord_status
case ('SLAVE_STATUS');    a_ptr%i => ele%slave_status
case ('ORIENTATION');     a_ptr%i => ele%orientation
case ('REF_TIME');        a_ptr%r => ele%ref_time
case ('KEY');             a_ptr%i => ele%key
case ('N_SLAVE');         a_ptr%i => ele%n_slave
case ('N_LORD');          a_ptr%i => ele%n_lord
case ('LR_FREQ_SPREAD', 'LR_SELF_WAKE_ON', 'LR_WAKE%AMP_SCALE', 'LR_WAKE%TIME_SCALE', &
      'LR_WAKE%FREQ_SPREAD', 'LR_WAKE%SELF_WAKE_ON', &
      'SR_WAKE%SCALE_WITH_LENGTH', 'SR_WAKE%AMP_SCALE', 'SR_WAKE%Z_SCALE', &
      'SR_WAKE%Z_LONG%SMOOTHING_SIGMA')
  if (.not. associated(ele%wake)) then
    if (.not. do_allocation) goto 9100
    call init_wake (ele%wake, 0, 0, 0, 0, .true.)
  endif
  select case (a_name)
  case ('SR_WAKE%Z_LONG%SMOOTHING_SIGMA');             a_ptr%r => ele%wake%sr%z_long%smoothing_sigma
  case ('SR_WAKE%AMP_SCALE');                          a_ptr%r => ele%wake%sr%amp_scale
  case ('SR_WAKE%Z_SCALE');                            a_ptr%r => ele%wake%sr%z_scale
  case ('SR_WAKE%SCALE_WITH_LENGTH');                  a_ptr%l => ele%wake%sr%scale_with_length
  case ('LR_SELF_WAKE_ON', 'LR_WAKE%SELF_WAKE_ON');    a_ptr%l => ele%wake%lr%self_wake_on
  case ('LR_WAKE%AMP_SCALE');                          a_ptr%r => ele%wake%lr%amp_scale
  case ('LR_WAKE%TIME_SCALE');                         a_ptr%r => ele%wake%lr%time_scale
  case ('LR_FREQ_SPREAD', 'LR_WAKE%FREQ_SPREAD');      a_ptr%r => ele%wake%lr%freq_spread
  end select

case ('H_MISALIGN%ACTIVE');     a_ptr%l => ele%photon%h_misalign%active
case ('DISPLACEMENT%ACTIVE');   a_ptr%l => ele%photon%displacement%active
case ('SEGMENTED%ACTIVE');      a_ptr%l => ele%photon%segmented%active

case ('H_MISALIGN%DR');         a_ptr%r1 => ele%photon%h_misalign%dr
case ('DISPLACEMENT%DR');       a_ptr%r1 => ele%photon%displacement%dr
case ('SEGMENTED%DR');          a_ptr%r1 => ele%photon%segmented%dr
case ('PIXEL%DR');              a_ptr%r1 => ele%photon%pixel%dr

case ('H_MISALIGN%R0');         a_ptr%r1 => ele%photon%h_misalign%r0
case ('DISPLACEMENT%R0');       a_ptr%r1 => ele%photon%displacement%r0
case ('SEGMENTED%R0');          a_ptr%r1 => ele%photon%segmented%r0
case ('PIXEL%R0');              a_ptr%r1 => ele%photon%pixel%r0
end select

if (a_name(1:11) == 'CURVATURE_X' .and. a_name(13:14) == '_Y' .and. a_name(16:) == '') then  ! Deprecated syntax
  ix = index('0123456789', a_name(12:12)) - 1
  iy = index('0123456789', a_name(15:15)) - 1
  if (ix == -1 .or. iy == -1) goto 9000 ! Error message and return
  if (ix > ubound(ele%photon%curvature%xy, 1)) goto 9000 ! Error message and return
  if (iy > ubound(ele%photon%curvature%xy, 2)) goto 9000 ! Error message and return
  a_ptr%r => ele%photon%curvature%xy(ix,iy)
  err_flag = .false.
  return
endif

if (a_name(1:10) == 'CURVATURE%') then
  select case (a_name(11:))
  case ('SPHERICAL')
    a_ptr%r => ele%photon%curvature%spherical
  case ('ELLIPTICAL_X')
    a_ptr%r => ele%photon%curvature%elliptical(1)
  case ('ELLIPTICAL_Y')
    a_ptr%r => ele%photon%curvature%elliptical(2)
  case ('ELLIPTICAL_Z')
    a_ptr%r => ele%photon%curvature%elliptical(3)
  case default
    if (a_name(11:11) /= 'X' .or. a_name(13:13) /= 'Y' .or. a_name(15:15) /= ' ') goto 9000
    ix = index('0123456789', a_name(12:12)) - 1
    iy = index('0123456789', a_name(14:14)) - 1
    if (ix == -1 .or. iy == -1) goto 9000 ! Error message and return
    if (ix > ubound(ele%photon%curvature%xy, 1)) goto 9000 ! Error message and return
    if (iy > ubound(ele%photon%curvature%xy, 2)) goto 9000 ! Error message and return
    a_ptr%r => ele%photon%curvature%xy(ix,iy)
  end select
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

if (associated(a_ptr%r) .or. associated(a_ptr%l) .or. associated(a_ptr%i) .or. &
    associated(a_ptr%r1) .or. associated(a_ptr%i1)) then
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

! Custom attribute.
! If a super_slave has multiple super_lords, it is not clear what do do. 
! In this case, return with err_flag = True.

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
  if (a_name(3:3) == 'S') then
    ixs = index('1XYZ', a_name(4:4))
    if (ixs == 0) return
    if (.not. associated(ele%spin_taylor(0)%term)) then
      if (.not. do_allocation) return
      do i = 0, 3
        call init_taylor_series(ele%spin_taylor(i), 0)
      enddo
    endif
    tlr => ele%spin_taylor(ixs-1)
    i0 = 5

  else
    n = index('123456', a_name(3:3))
    if (.not. associated(ele%taylor(1)%term)) then
      if (.not. do_allocation) return
      do i = 1, 6
        call init_taylor_series(ele%taylor(i), 0)
      enddo
    endif
    tlr => ele%taylor(n)
    i0 = 4
  endif

  expn = 0
  do i = i0, len_trim(a_name)
    j = index('123456', a_name(i:i))
    expn(j) = expn(j) + 1
  enddo

  i = taylor_term_index(tlr, expn, do_allocation)
  if (i /= 0) then
    a_ptr%r => tlr%term(i)%coef
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
case ('MATRIX');                         a_ptr%r => ele%value(matrix$)
case ('KICK0');                          a_ptr%r => ele%value(kick0$)
case ('FLEXIBLE');                       a_ptr%r => ele%value(flexible$)
case ('MODE_FLIP0');                     a_ptr%r => ele%value(mode_flip0$)
case ('MODE_FLIP1');                     a_ptr%r => ele%value(mode_flip1$)
case ('X_REF');                          a_ptr%r => ele%taylor(1)%ref
case ('PX_REF');                         a_ptr%r => ele%taylor(2)%ref
case ('Y_REF');                          a_ptr%r => ele%taylor(3)%ref
case ('PY_REF');                         a_ptr%r => ele%taylor(4)%ref
case ('Z_REF');                          a_ptr%r => ele%taylor(5)%ref
case ('PZ_REF');                         a_ptr%r => ele%taylor(6)%ref
case ('SYMPLECTIFY');                    a_ptr%l => ele%symplectify
case ('ABSOLUTE_TIME_TRACKING');         a_ptr%l => bmad_com%absolute_time_tracking
case ('ABSOLUTE_TIME_REF_SHIFT');        a_ptr%l => bmad_com%absolute_time_ref_shift
case ('TAYLOR_MAP_INCLUDES_OFFSETS');    a_ptr%l => ele%taylor_map_includes_offsets
case ('OFFSET_MOVES_APERTURE');          a_ptr%l => ele%offset_moves_aperture
case ('FIELD_MASTER');                   a_ptr%l => ele%field_master
case ('SCALE_MULTIPOLES');               a_ptr%l => ele%scale_multipoles
case ('MULTIPOLES_ON');                  a_ptr%l => ele%multipoles_on
case ('IS_ON');                          a_ptr%l => ele%is_on
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
case ('MAT6_CALC_METHOD');               a_ptr%i => ele%mat6_calc_method
case ('MODE');                           a_ptr%r => ele%value(mode$)
case ('ORIGIN_ELE_REF_PT');              a_ptr%r => ele%value(origin_ele_ref_pt$)
case ('PTC_INTEGRATION_TYPE');           a_ptr%i => ele%ptc_integration_type
case ('PTC_FRINGE_GEOMETRY');            a_ptr%r => ele%value(ptc_fringe_geometry$)
case ('PTC_FIELD_GEOMETRY');             a_ptr%r => ele%value(ptc_field_geometry$)
case ('REF_ORBIT_FOLLOWS');              a_ptr%r => ele%value(ref_orbit_follows$)
case ('REF_COORDS');                     a_ptr%r => ele%value(ref_coords$)
case ('SCATTER_METHOD');                 a_ptr%r => ele%value(scatter_method$)
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
case ('UPSTREAM_ELE_DIR')
case ('DOWNSTREAM_ELE_DIR')
end select

if (associated(a_ptr%r) .or. associated(a_ptr%i) .or. associated(a_ptr%l) .or. &
    associated(a_ptr%r1) .or. associated(a_ptr%i1)) then
  err_flag = .false.
else
  goto 9000
endif

return

!----------------------------------------
! Error message and return

9000 continue
if (do_print) call out_io (s_error$, r_name, &
          'INVALID ATTRIBUTE: ' // attrib_name, 'FOR THIS ELEMENT: ' // ele%name)
return

!----------------------------------------
9100 continue
if (do_print) call out_io (s_error$, r_name, &
                 'WAKE ATTRIBUTE NOT ALLOCATED: ' // attrib_name, &
                 'FOR THIS ELEMENT: ' // ele%name)
return

!----------------------------------------
9110 continue
if (do_print) call out_io (s_error$, r_name, &
                 'R_CUSTOM ATTRIBUTE NOT ALLOCATED: ' // attrib_name, &
                 'FOR THIS ELEMENT: ' // ele%name)
return

!----------------------------------------
9130 continue
if (do_print) then
  if (out_of_bounds) then
    call out_io (s_error$, r_name, &
        'INDEX OUT OF BOUNDS IN ATTRIBUTE: ' // attrib_name, &
        'FOR THIS ELEMENT: ' // ele%name)
  else
    call out_io (s_error$, r_name, &
        'MALFORMED ATTRIBUTE: ' // attrib_name, &
        'FOR THIS ELEMENT: ' // ele%name)
  endif
endif
return

!----------------------------------------
9140 continue
if (do_print) call out_io (s_error$, r_name, &
                 '(EM) FIELD ATTRIBUTE NOT ALLOCATED: ' // attrib_name, &
                 'FOR THIS ELEMENT: ' // ele%name)
return

!----------------------------------------
9200 continue
if (do_print) call out_io (s_error$, r_name, &
        'BAD VERTEX COMPONENT IN ATTRIBUTE: ' // attrib_name, &
        'FOR THIS ELEMENT: ' // ele%name)
return

!----------------------------------------
9210 continue
if (do_print) call out_io (s_error$, r_name, &
        'CROSS-SECTION NOT DEFINED SO CANNOT SET ATTRIBUTE: ' // attrib_name, &
        'FOR THIS ELEMENT: ' // ele%name)
return

!----------------------------------------
9300 continue
if (do_print) call out_io (s_error$, r_name, &
        'ATTRIBUTE: ' // trim(attrib_name) // ' IS NOT A KNOT POINT', &
        'FOR THIS ELEMENT: ' // ele%name)
return

!----------------------------------------
9310 continue
if (do_print) call out_io (s_error$, r_name, &
        'ATTRIBUTE: ' // trim(attrib_name) // ' IS NOT A NAMED PARAMETER OF CONTROL EXPRESSION.', &
        'FOR THIS ELEMENT: ' // ele%name)
return

!----------------------------------------
9320 continue
if (do_print) call out_io (s_error$, r_name, &
        'ATTRIBUTE: ' // trim(attrib_name) // ' IS NOT VALID SINCE CONTROLLER DOES NOT USE KNOTS.', &
        'FOR THIS ELEMENT: ' // ele%name)
return

!----------------------------------------
9400 continue
if (do_print) call out_io (s_error$, r_name, &
        'ATTRIBUTE: ' // trim(attrib_name) // ' IS NOT VALID FOR THIS ELEMENT: ' // ele%name)
return

!----------------------------------------
9410 continue
if (do_print) call out_io (s_error$, r_name, &
        'ELEMENT: ' // trim(ele%name) // ' USES "FREQUENCIES" TO SPECIFY WAVEFORM. NOT "AMP_VS_TIME"')
return

!----------------------------------------
9420 continue
if (do_print) call out_io (s_error$, r_name, &
        'ATTRIBUTE: ' // trim(attrib_name) // ' HAS INDEX OUT OF RANGE. VALID RANGE IS FROM 1 TO ', &
                                                                      int_str(size(ele%ac_kick%amp_vs_time)), &
        'FOR ELEMENT: ' // ele%name)
return

!----------------------------------------
9430 continue
if (do_print) call out_io (s_error$, r_name, &
        'ATTRIBUTE: ' // trim(attrib_name) // ' IS NOT A VALID AMP_VS_TIME COMPONENT.', &
        'FOR ELEMENT: ' // ele%name)
return

!----------------------------------------
9450 continue
if (do_print) call out_io (s_error$, r_name, &
        'ELEMENT: ' // trim(ele%name) // ' USES "AMP_VS_TIME" TO SPECIFY WAVEFORM. NOT "FREQUENCIES"')
return

!----------------------------------------
9460 continue
if (do_print) call out_io (s_error$, r_name, &
        'ATTRIBUTE: ' // trim(attrib_name) // ' HAS INDEX OUT OF RANGE. VALID RANGE IS FROM 1 TO ', &
                                                                      int_str(size(ele%ac_kick%frequency)), &
        'FOR ELEMENT: ' // ele%name)
return

!----------------------------------------
9470 continue
if (do_print) call out_io (s_error$, r_name, &
        'ATTRIBUTE: ' // trim(attrib_name) // ' IS NOT A VALID FREQUENCIES COMPONENT.', &
        'FOR ELEMENT: ' // ele%name)
return

!---------------------------------------------------------------
contains

!+
! Function reads number of the form "...(num)" and checks to
! see if num is between n_min and n_max.
! Function also chops "...(num)" from name.
!-

function get_this_index(name, ix_name, err, n_min, n_max) result (ixc)

character(*) name

integer ix_name, n_min, n_max, ixc, ios

logical err

!

err = .true.
ixc = int_garbage$

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
