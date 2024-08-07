!+
! Subroutine ele_geometry (floor_start, ele, floor_end, len_scale, ignore_patch_err)
!
! Routine to calculate the non-misaligned global (floor) coordinates of an element at the downstream end
! (if len_scale = 1) given the non-misaligned global coordinates of the preceeding element. 
!
! The coordinates are computed without misalignments. That is, the coordinates are the "laboratory" 
! coordinates and not the "body" coordinates. To compute coordinates with misalignments, use
! the routine ele_geometry_with_misalignments.
!
! Note: For crystal and mirror elements, len_scale ~ 0.5 means compute the floor coordinates of the
! element surface.
!
! Note: For a floor_position element, floor_end is independent of floor_start.
!
! Input:
!   floor_start      -- Starting floor coordinates at upstream end.
!                         Not used for fiducial and girder elements.
!   ele              -- Ele_struct: Element to propagate the geometry through.
!   len_scale        -- Real(rp), optional: factor to scale the length of the element.
!                          1.0_rp => Output is geometry at end of element (default).
!                          0.5_rp => Output is geometry at center of element. 
!                         -1.0_rp => Used to propagate geometry in reverse.
!   ignore_patch_err -- logical, optional: If present and True, ignore flexible patch errors.
!                         This is used by ele_compute_ref_energy_and_time to suppress unnecessary messages.
!
! Output:
!   floor_end        -- floor_position_struct, optional: Output floor position. If not present then 
!                         ele%floor will be used and ele%bookkeeping_state%floor_position will be set to ok$.
!     %r(3)              -- X, Y, Z Floor position at end of element
!     %w(3,3)            -- W matrix corresponding to orientation angles
!     %theta, phi, %psi  -- Orientation angles 
!-

recursive subroutine ele_geometry (floor_start, ele, floor_end, len_scale, ignore_patch_err)

use bmad_interface, dummy => ele_geometry

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0, ele00, slave0, slave1, ele2, this_ele, slave
type (floor_position_struct), optional, target :: floor_end
type (floor_position_struct) :: floor_start, this_floor, old_floor, floor0
type (floor_position_struct), pointer :: floor
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_pointer_struct), allocatable, target :: chain_ele(:)
type (lat_param_struct) param

real(rp), optional :: len_scale
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx), dtheta
real(rp) r0(3), w0_mat(3,3), rot_angle, graze_angle_in, graze_angle_out
real(rp) chord_len, angle, ang, leng, rho, len_factor
real(rp) theta, phi, psi, tlt, dz(3), z0(3), z_cross(3), eps, signif(6)
real(rp) :: w_mat(3,3), w_mat_inv(3,3), s_mat(3,3), r_vec(3), t_mat(3,3)

integer i, j, k, n, ie, key, n_loc, ix_pass, n_links, ix_pole_max, ib_to, ix, iv(6)

logical err, doit, finished, has_multipole_rot_tilt, ele_floor_geometry_calc
logical, optional :: ignore_patch_err

character(*), parameter :: r_name = 'ele_geometry'

! Only set ele%bookkeeping_state%floor_position = ok$ if computing ele%floor.

if (present(floor_end)) then
  floor => floor_end
  ele_floor_geometry_calc = .false.
else
  floor => ele%floor
  ele_floor_geometry_calc = .true.
  ele%bookkeeping_state%floor_position = ok$
endif

! Custom geometry

if (associated(ele_geometry_hook_ptr)) then
  call ele_geometry_hook_ptr (floor_start, ele, floor_end, finished, len_scale)
  if (finished) return
endif

! Init.

len_factor = ele%orientation * real_option(1.0_rp, len_scale)

if (ele%key == crystal$) then
  if (ele%value(graze_angle_in$) /= 0) then
    graze_angle_in  = ele%value(graze_angle_in$) 
    graze_angle_out = ele%value(graze_angle_out$) 
  else
    graze_angle_in  = ele%value(bragg_angle_in$) 
    graze_angle_out = ele%value(bragg_angle_out$) 
  endif
endif

old_floor = floor
floor0 = floor_start  ! Use floor0 in case actual args floor = floor_start

theta   = floor0%theta
phi     = floor0%phi
psi     = floor0%psi
w_mat   = floor0%w

leng = ele%value(l$) * len_factor

key = ele%key
if ((key == sbend$ .or. ele%key == rf_bend$) .and. (leng == 0 .or. ele%value(g$) == 0)) key = drift$

! Fiducial, floor_shift and girder elements.
! Note that fiducial, and girder elements are independent of floor0

if (key == fiducial$ .or. key == girder$ .or. key == floor_shift$) then
  ele0 => null()

  ! Fiducial with no origin ele: Use global origin.
  if (key == fiducial$ .and. ele%component_name == '') then
    call mat_make_unit(w_mat)
    r0 = 0

  ! Girder with global origin
  elseif (key == girder$ .and. ele%component_name == 'GLOBAL_COORDINATES') then
    call mat_make_unit(w_mat)
    r0 = 0

  ! Girder uses center of itself by default.
  else  if (key == girder$ .and. ele%component_name == '') then
    call find_element_ends (ele, slave0, slave1)
    r0 = (slave0%floor%r + slave1%floor%r) / 2
    ! 
    w0_mat = slave0%floor%w
    dz = r0 - slave0%floor%r
    z0 = w0_mat(:,3)
    z_cross = cross_product(z0, dz)
    if (all(dz == 0) .or. sum(abs(z_cross)) <= 1d-14 * sum(abs(dz))) then
      w_mat = w0_mat
    else
      angle = atan2(sqrt(dot_product(z_cross, z_cross)), dot_product(z0, dz))
      call axis_angle_to_w_mat (z_cross, angle, w_mat)
      w_mat = matmul (w_mat, w0_mat)
    endif


  else  ! Must have a reference element
    if (ele%component_name == '') then  ! Must be a floor_shift element
      ele0 => pointer_to_next_ele(ele, -1)

    else
      call lat_ele_locator (ele%component_name, ele%branch%lat, eles, n_loc, err)

      ! If multiple matches but one match is a slave of ele (which must be a girder), this is OK.
      if (n_loc == 0) then
        call out_io (s_fatal$, r_name, 'ORIGIN_ELE NAME: ' // ele%component_name,  &
                                       'FOR ELEMENT: ' // ele%name, &
                                       'DOES NOT MATCH ANY ELEMENT!')
        if (global_com%exit_on_error) call err_exit
        return
      endif

      if (n_loc > 0) then
        n = 0
        do i = 1, n_loc
          do j = 1, ele%n_slave
            slave => pointer_to_slave(ele, j)
            if (.not. (ele_loc(eles(i)%ele) == ele_loc(slave))) cycle
            select case (n)
            case (0);     n = j
            case default; n = -1  ! Mark that there are multiple slave element matches
            end select
          enddo
        enddo
        if (n > 0) then
          n_loc = 1
          eles(1)%ele => pointer_to_slave(ele, n)
        endif
        if (n_loc > 1) then
          call out_io (s_fatal$, r_name, 'ORIGIN_ELE: ' // ele%component_name,  &
                                         'FOR ELEMENT: ' // ele%name, &
                                         'IS NOT UNIQUE!')
          if (global_com%exit_on_error) call err_exit
          return
        endif
      endif

      if (ele%lord_status == multipass_lord$ .or. ele%key == ramper$) then
        call out_io (s_fatal$, r_name, 'ORIGIN_ELE: ' // ele%component_name,  &
                                       'FOR ELEMENT: ' // ele%name, &
                                       'IS A MULTIPASS_LORD OR RAMPER WHICH DOES NOT HAVE A UNIQUE POSITION')
        if (global_com%exit_on_error) call err_exit
        return
      endif

      if (ele%n_slave > 1 .and. (ele%key == overlay$ .or. ele%key == group$)) then
        call out_io (s_fatal$, r_name, 'ORIGIN_ELE: ' // ele%component_name,  &
                                       'FOR ELEMENT: ' // ele%name, &
                                       'IS AN OVERLAY OR GROUP ELEMENT WHICH HAS MORE THAN ONE SLAVE SO THERE IS NO UNIQUE POSITION')
        if (global_com%exit_on_error) call err_exit
        return        
      endif

      ele0 => eles(1)%ele
    endif

    select case (stream_ele_end(nint(ele%value(origin_ele_ref_pt$)), ele%orientation))
    case (upstream_end$)
      call ele_geometry (ele0%floor, ele0, this_floor, -1.0_rp)
      w_mat = this_floor%w
      r0 = this_floor%r

    case (center_pt$)
      select case (ele0%key)
      case (crystal$, mirror$, multilayer_mirror$)
        ele00 => pointer_to_next_ele(ele0, -1)
        w_mat = ele00%floor%w

        select case (ele0%key)
        case (crystal$)
          if (ele0%value(tilt_corr$) /= 0) then
            t_mat = w_mat_for_tilt(ele0%value(tilt_corr$))
            w_mat = matmul(w_mat, t_mat)
          endif
          rot_angle = ele0%value(bragg_angle_in$)
          if (ele0%value(graze_angle_in$) /= 0) rot_angle = ele0%value(graze_angle_in$)
          if (ele0%value(b_param$) < 0) rot_angle = rot_angle - pi/2  ! Bragg
        case (mirror$, multilayer_mirror$)
          rot_angle = ele0%value(graze_angle$) - pi/2
        end select

        s_mat = w_mat_for_bend_angle (rot_angle, ele0%value(ref_tilt_tot$))
        w_mat = matmul (w_mat, s_mat)
        r0 = ele00%floor%r

      case default
        call ele_geometry (ele0%floor, ele0, this_floor, -0.5_rp)
        w_mat = this_floor%w
        r0 = this_floor%r
      end select

    case (downstream_end$)
      w_mat = ele0%floor%w
      r0 = ele0%floor%r

    case default
      call out_io (s_fatal$, r_name, 'ORIGIN_ELE_REF_PT NOT VALID FOR ELEMENT: ' // ele%name)
      if (global_com%exit_on_error) call err_exit
    end select
  endif

  ! Now offset from origin pt. 
  ! Fiducial and floor_shift elements are not allowed to be be turned off.

  if (ele%key == girder$ .and. .not. ele%is_on) then
    r_vec = 0
    theta = 0
  elseif (ele%key == floor_shift$) then
    r_vec = [ele%value(x_offset$), ele%value(y_offset$), ele%value(z_offset$)]
    theta = ele%value(x_pitch$);  phi = ele%value(y_pitch$);  psi = ele%value(tilt$)
  else
    r_vec = [ele%value(dx_origin$), ele%value(dy_origin$), ele%value(dz_origin$)]
    theta = ele%value(dtheta_origin$);  phi = ele%value(dphi_origin$); psi = ele%value(dpsi_origin$)
  endif

  floor%r = r0 + matmul(w_mat, r_vec)
  call floor_angles_to_w_mat (theta, phi, psi, s_mat)
  w_mat = matmul(w_mat, s_mat)
  
  ! Update floor angles
  floor%w = w_mat

  if (ele%key == fiducial$) then
    ! floor0 for a fiducial element is meaningless. Instead use approximate floor coords to try to
    ! have the angles come out without annoying factors of 2pi offsets
    this_floor%theta = theta; this_floor%phi = phi; this_floor%psi = psi
    if (associated(ele0)) then
      this_floor%theta = this_floor%theta + ele0%floor%theta
      this_floor%phi   = this_floor%phi   + ele0%floor%phi
      this_floor%psi   = this_floor%psi   + ele0%floor%psi
    endif
    call update_floor_angles(floor, this_floor)
  else
    call update_floor_angles(floor, floor0)
  endif

  call end_bookkeeping(ele, old_floor, floor)
  return
endif   ! Fiducial, girder, floor_shift

!---------------------------
! General case where layout is not in the horizontal plane
! Note: 

has_multipole_rot_tilt = .false.
if (key == multipole$) then
  call multipole_ele_to_kt (ele, .true., ix_pole_max, knl, tilt)
  if (knl(0) /= 0 .and. tilt(0) /= 0) has_multipole_rot_tilt = .true.
endif

if (((key == mirror$  .or. key == sbend$ .or. key == rf_bend$ .or. key == multilayer_mirror$) .and. &
         ele%value(ref_tilt_tot$) /= 0) .or. phi /= 0 .or. psi /= 0 .or. key == patch$ .or. &
         key == crystal$ .or. has_multipole_rot_tilt) then

  select case (key)

  ! sbend and multipole

  case (sbend$, rf_bend$, multipole$)
    if (key == sbend$ .or. key == rf_bend$) then
      angle = leng * dble(ele%value(g$))
      tlt = ele%value(ref_tilt_tot$)
      rho = 1.0_dp / ele%value(g$)
      r_vec = [rho * cos_one(angle), 0.0_dp, rho * sin(angle)]
    else
      angle = knl(0) * len_factor
      tlt = tilt(0)
      r_vec = 0
    endif

    s_mat = w_mat_for_bend_angle (angle, tlt, r_vec)
    call floor_angles_to_w_mat (theta, phi, psi, w_mat)
    floor%r = floor0%r + matmul(w_mat, r_vec)
    w_mat = matmul (w_mat, s_mat)

  ! mirror, multilayer_mirror, crystal
  ! Note: The reference frame is discontinuous.

  case (mirror$, multilayer_mirror$, crystal$)
    
    ! Len_factor of 1/2 means compute geometry of surface
    if (len_factor > 0.75) then
      len_factor = 1.0_rp
    elseif (len_factor > 0.25) then
      len_factor = 0.5_rp
    elseif (len_factor > -0.50) then
      len_factor = 0
    else
      len_factor = -1.0_rp
    endif

    if (ele%key == crystal$) then
      select case (nint(ele%value(ref_orbit_follows$)))
      case (bragg_diffracted$)
        if (len_factor == 0.5_rp) then
          angle = graze_angle_in
        else
          angle = len_factor * (graze_angle_in + graze_angle_out)
        endif
      case (forward_diffracted$, undiffracted$)
        angle = 0
      end select

      if (ele%value(b_param$) > 0 .and. len_factor /= 0.5_rp) then  ! Laue
        ! %l_ref is with respect to the body coords
        r_vec = len_factor * ele%photon%material%l_ref
        if (len_factor > 0) then  ! Forward propagation -> Express r_vec in entrance coords.
          ang = len_factor * graze_angle_in
          r_vec = [cos(ang) * r_vec(1) - sin(ang) * r_vec(3), r_vec(2), sin(ang) * r_vec(1) + cos(ang) * r_vec(3)]
        else                      ! Express r_vec in exit coords
          ang = angle - len_factor * graze_angle_in
          r_vec = [cos(ang) * r_vec(1) - sin(ang) * r_vec(3), r_vec(2), sin(ang) * r_vec(1) + cos(ang) * r_vec(3)]
        endif
      else  ! Bragg
        r_vec = 0
      endif

    else   ! Non-crystal
      angle = 2 * len_factor * ele%value(graze_angle$)
      r_vec = 0
    endif

    !

    call floor_angles_to_w_mat (theta, phi, psi, w_mat)
    floor%r = floor0%r + matmul(w_mat, r_vec)

    if (len_factor == 0.5_rp) then
      if (ele%key == crystal$) then
        s_mat = w_mat_for_bend_angle(angle-pi/2, ele%value(ref_tilt_tot$) + ele%value(tilt_corr$), r_vec)
      else
        s_mat = w_mat_for_bend_angle(angle-pi/2, ele%value(ref_tilt_tot$), r_vec)
      endif
    else
      s_mat = w_mat_for_bend_angle(angle, ele%value(ref_tilt_tot$), r_vec)
    endif
    w_mat = matmul (w_mat, s_mat)

  ! patch

  case (patch$)

    ! Flexible bookkeeping to compute offsets and pitches. 
    ! Only needs to be done if element is part of a lattice (as opposed to, for example, a slice_slave) and,
    ! if part of a multipass region, only needs to be done if this is a first pass slave.

    if (is_true(ele%value(flexible$)) .and. ele%ix_ele > 0) then
      doit = .true.
      if (ele%lord_status == multipass_lord$) doit = .false.
      call multipass_chain (ele, ix_pass, n_links)
      if (ix_pass > 1) doit = .false.

      if (doit) then
        ! In the case that the next element has zero length (and so does not affect the floor position),
        ! and does not yet have a well defined position, look at the element after.
        ele2 => ele
        do
          ! If ele2 is a super_slave of a super_lord that is also a multipass_slave there is a potential problem 
          ! in that, currently, multipass_chain will bomb if not all of the associated multipass_slave elements are
          ! also super_lords (this can happen if, say, only one of the multipass_slave elements is superimposed upon).
          ! Using use_super_lord = T here avoids this problem.
          ele2 => pointer_to_next_ele(ele2, 1)
          call multipass_chain (ele2, ix_pass, n_links, chain_ele, use_super_lord = .true.)
          if (ix_pass > 1) then
            ele2 => chain_ele(1)%ele
            if (ele2%lord_status == super_lord$) ele2 => pointer_to_slave(ele2, 1)
          endif
          if (ele2%bookkeeping_state%floor_position /= stale$) exit
          if (ele2%value(l$) /= 0 .or. ele2%key == patch$) exit
        enddo

        if (ele2%bookkeeping_state%floor_position == stale$ .and. .not. logic_option(.false., ignore_patch_err)) then
          call out_io (s_fatal$, r_name, 'FOR FLEXIBLE PATCH: ' // trim(ele%name) // '  ' // trim(ele_loc_name(ele, parens = "()")), &
                                         '"DOWNSTREAM" FIDUCIAL ELEMENT: ' // trim(ele2%name) // '  ' // trim(ele_loc_name(ele2, parens = "()")), &
                                         ' DOES NOT HAVE A WELL DEFINED POSITION')
          if (global_com%exit_on_error) call err_exit
        endif

        call ele_geometry (ele2%floor, ele2, this_floor, -1.0_rp) 
        
        ele0 => pointer_to_next_ele(ele, -1)
        w_mat_inv = transpose(ele0%floor%w)
        w_mat = this_floor%w
        w_mat = matmul(w_mat_inv, w_mat)
        call floor_w_mat_to_angles (w_mat, ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$))
        r_vec = matmul(w_mat_inv, this_floor%r - ele0%floor%r)
        ele%value(x_offset$) = r_vec(1)
        ele%value(y_offset$) = r_vec(2)
        ele%value(z_offset$) = r_vec(3)
        ele%value(l$) = patch_length(ele)

        iv = [x_offset$, y_offset$, z_offset$, x_pitch$, y_pitch$, tilt$]
        signif(1:3) = bmad_com%significant_length;  signif(4:6) = 1d-12
        if (ele_floor_geometry_calc .and. ele_value_has_changed(ele, iv, signif, .true.)) then
          call set_ele_status_stale (ele, s_position_group$)
        endif

        ! Transfer offsets and pitches if patch is a multipass_slave.
        if (ele%slave_status == multipass_slave$) then
          call multipass_chain (ele, ix_pass, n_links, chain_ele)
          do i = 1, n_links
            if (i == 1) then
              ele2 => pointer_to_lord(ele, 1)
            else
              ele2 => chain_ele(i)%ele
            endif
            ele2%value(x_offset$) = ele%value(x_offset$)
            ele2%value(y_offset$) = ele%value(y_offset$)
            ele2%value(z_offset$) = ele%value(z_offset$)
            ele2%value(x_pitch$)  = ele%value(x_pitch$)
            ele2%value(y_pitch$)  = ele%value(y_pitch$)
            ele2%value(tilt$)     = ele%value(tilt$)
            ele2%value(l$)        = ele%value(l$)
          enddo
        endif

      endif
    endif

    ! Now calculate the geometry.

    r_vec = [ele%value(x_offset$), ele%value(y_offset$), ele%value(z_offset$)]

    call floor_angles_to_w_mat (theta, phi, psi, w_mat)
    if (len_factor < 0) then
      call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$), w_mat_inv = s_mat)
      w_mat = matmul(w_mat, s_mat)
      floor%r = floor0%r - matmul(w_mat, r_vec)
    else
      floor%r = floor0%r + matmul(w_mat, r_vec)
      call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$), s_mat)
      w_mat = matmul(w_mat, s_mat)
    endif

  ! Everything else. Just a translation

  case default
    floor%r = floor0%r + w_mat(:,3) * leng

  end select 

! Simple case where the local reference frame stays in the horizontal plane.

else

  select case (key)
  case (sbend$, rf_bend$)
    angle = leng * ele%value(g$)
    chord_len = 2 * ele%value(rho$) * sin(angle/2)
  case (multipole$)
    angle = knl(0)
    chord_len = 0
  case (mirror$, multilayer_mirror$)
    angle = 2 * ele%value(graze_angle$) * len_factor
    chord_len = 0
  case default
    angle = 0
    chord_len = leng
  end select

  theta = theta - angle / 2
  floor%r(1) = floor0%r(1) + chord_len * sin(theta)
  floor%r(2) = floor0%r(2)
  floor%r(3) = floor0%r(3) + chord_len * cos(theta)
  theta = theta - angle / 2

  call rotate_mat(w_mat, y_axis$, -angle)

endif

! Update floor angles

floor%w = w_mat 
call update_floor_angles(floor, floor0)
call end_bookkeeping(ele, old_floor, floor)

!-------------------------------------------------------------------------------------------
contains

subroutine end_bookkeeping(ele, old_floor, floor)

type (ele_struct), target :: ele
type (floor_position_struct) old_floor, floor
type (ele_struct), pointer :: lord, slave, ele2
type (lat_struct), pointer :: lat
type (ele_pointer_struct), allocatable, target :: chain_ele(:)

integer k, ib_to, ix, ix_pass, n_links, ie

! End bookkeeping. Only set ele%bookkeeping_state if computing 
! ele%floor (ele_floor_geometry_calc = T) and element is associated with a lattice...

eps = bmad_com%significant_length
if (ele_floor_geometry_calc .and. associated(ele%branch) .and. (any(abs(floor%r - old_floor%r) > eps) .or. &
              abs(floor%theta - old_floor%theta) > eps .or. abs(floor%phi - old_floor%phi) > eps .or. &
              abs(floor%psi - old_floor%psi) > eps)) then

  lat => ele%branch%lat

  ! If there is a girder element then *_tot attributes need to be recomputed.
  do k = 1, ele%n_lord
    lord => pointer_to_lord(ele, k)
    if (lord%lord_status == girder_lord$) ele%bookkeeping_state%control = stale$
  enddo

  if (ele%key == girder$) then
    do k = 1, ele%n_slave
      slave => pointer_to_slave(ele, k)
      slave%bookkeeping_state%control = stale$
      lat%branch(slave%ix_branch)%param%bookkeeping_state%control = stale$
    enddo
  endif

  ! Fork target branch only needs to be recomputed if target branch index is greater than present branch.
  if (ele%key == fork$ .or. ele%key == photon_fork$) then
    ib_to = nint(ele%value(ix_to_branch$))
    if (ib_to > ele%ix_branch) then
      ix = nint(ele%value(ix_to_element$))
      lat%branch(ib_to)%ele(ix)%bookkeeping_state%floor_position = stale$
      lat%branch(ib_to)%param%bookkeeping_state%floor_position = stale$
    endif
  endif

  call multipass_chain(ele, ix_pass, n_links, chain_ele, use_super_lord = .true.)
  if (ix_pass > 0) then
    do k = ix_pass+1, n_links
      this_ele => chain_ele(k)%ele
      this_ele%bookkeeping_state%floor_position = stale$
      lat%branch(this_ele%ix_branch)%param%bookkeeping_state%floor_position = stale$
      if (this_ele%lord_status == super_lord$) then
        do ie = 1, this_ele%n_slave
          ele2 => pointer_to_slave(this_ele, ie)
          ele2%bookkeeping_state%floor_position = stale$
        enddo
      endif
    enddo
  endif

  if (ele%slave_status == super_slave$) then
    do k = 1, ele%n_lord
      lord => pointer_to_lord(ele, k)
      if (lord%lord_status /= super_lord$) exit
      lord%bookkeeping_state%floor_position = stale$
      if (lord%slave_status == multipass_slave$) then
        lord => pointer_to_lord(lord, 1)  ! multipass lord
        lord%bookkeeping_state%floor_position = stale$
      endif
    enddo

  elseif (ele%slave_status == multipass_slave$) then
    lord => pointer_to_lord(ele, 1)
    lord%bookkeeping_state%floor_position = stale$
  endif
endif

end subroutine end_bookkeeping

end subroutine ele_geometry
