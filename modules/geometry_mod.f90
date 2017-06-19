module geometry_mod

use multipass_mod
use lat_ele_loc_mod
use rotation_3d_mod

implicit none

contains

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine lat_geometry (lat)
!
! Routine to calculate the physical placement of all the elements in a lattice.
! That is, the layout on the floor. This is the same as the MAD convention.
!
! Note: This routine does NOT update %ele(i)%s. To do this call s_calc.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat -- lat_struct: The lattice.
!     %ele(0)%floor  -- Floor_position_struct: The starting point for the calculations.
!
! Output:
!   lat -- lat_struct: The lattice.
!     %ele(i)%floor --  floor_position_struct: Floor position.
!-

subroutine lat_geometry (lat)

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, lord, slave, b_ele, ele0, ele2
type (branch_struct), pointer :: branch
type (floor_position_struct) dummy

integer i, i2, n, ix2, ie, ib, ie0
logical stale

character(16), parameter :: r_name = 'lat_geometry'

!

do n = 0, ubound(lat%branch, 1)
  branch => lat%branch(n)

  if (bmad_com%auto_bookkeeper) then
    stale = .true.
  else
    if (branch%param%bookkeeping_state%floor_position /= stale$) cycle
    stale = .false.
  endif

  branch%param%bookkeeping_state%floor_position = ok$
  branch%ele(0)%value(floor_set$) = false$

  ! If there are fiducial elements then survey the fiducial regions

  do i = 1, branch%n_ele_track
    ele => branch%ele(i)
    if (ele%key /= fiducial$) cycle
    call ele_geometry (dummy, ele, ele%floor, set_ok = .true.)

    do i2 = i+1, branch%n_ele_track
      ele2 => branch%ele(i2)
      if (ele2%key == patch$ .and. is_true(ele2%value(flexible$))) exit
      if (ele2%key == fiducial$) then
        call out_io (s_fatal$, r_name, 'FIDUCIAL ELEMENTS IN A BRANCH MUST BE SEPARATED BY A FLEXIBLE PATCH')
        if (global_com%exit_on_error) call err_exit
        exit
      endif
      call ele_geometry (branch%ele(i2-1)%floor, ele2, ele2%floor, set_ok = .true.)
    enddo

    branch%ele(i-1)%floor = ele%floor  ! Save time

    do i2 = i-1, 1, -1
      ele2 => branch%ele(i2)
      if (ele2%key == patch$ .and. is_true(ele2%value(flexible$))) exit
      if (ele2%key == fiducial$) then
        call out_io (s_fatal$, r_name, 'FIDUCIAL ELEMENTS IN A BRANCH MUST BE SEPARATED BY A FLEXIBLE PATCH')
        if (global_com%exit_on_error) call err_exit
        exit
      endif
      call ele_geometry (ele2%floor, ele2, branch%ele(i2-1)%floor, -1.0_rp)
      if (i2 - 1 == 0) branch%ele(i2-1)%value(floor_set$) = true$
      branch%ele(i2-1)%bookkeeping_state%floor_position = ok$
    enddo
  enddo

  ! Transfer info from the from_branch element if that element exists.

  if (branch%ix_from_branch > -1 .and. (stale .or. branch%ele(0)%bookkeeping_state%floor_position == stale$)) then
    b_ele => pointer_to_ele (lat, branch%ix_from_ele, branch%ix_from_branch)
    ie0 = nint(b_ele%value(ix_to_element$))
    call ele_geometry (b_ele%floor, b_ele, branch%ele(ie0)%floor)
    branch%ele(ie0)%bookkeeping_state%floor_position = ok$
    stale = .true.
  else
    ie0 = 0
  endif

  if (branch%ele(ie0)%bookkeeping_state%floor_position == stale$) branch%ele(ie0)%bookkeeping_state%floor_position = ok$

  do i = ie0+1, branch%n_ele_track
    call propagate_geometry(i, 1, stale)
  enddo

  do i = ie0-1, 0, -1
    call propagate_geometry(i, -1, stale)
  enddo

enddo

! put info in super_lords and multipass_lords

lat%lord_state%floor_position = ok$
lat%param%bookkeeping_state%floor_position = ok$

if (bmad_com%auto_bookkeeper) lat%ele(lat%n_ele_track+1:lat%n_ele_max)%bookkeeping_state%floor_position = stale$

do i = lat%n_ele_track+1, lat%n_ele_max  
  lord => lat%ele(i)

  if (lord%bookkeeping_state%floor_position /= stale$) cycle
  lord%bookkeeping_state%floor_position = ok$

  if (lord%n_slave == 0) cycle

  select case (lord%lord_status)
  case (super_lord$)
    slave => pointer_to_slave(lord, lord%n_slave) ! Last slave is at exit end.
    lord%floor = slave%floor
  case (multipass_lord$)
    slave => pointer_to_slave(lord, 1)
    lord%floor = slave%floor
  case (girder_lord$)
    call girder_lord_geometry (lord)
  end select

enddo

!---------------------------------------------------
contains

subroutine propagate_geometry (ie, dir, stale)

type (floor_position_struct) floor0
type (ele_pointer_struct), allocatable :: chain_ele(:)

integer ie, dir, ix, ix_pass, n_links, k
logical stale

!

ele => branch%ele(ie)

if (ele%key == patch$ .and. is_true(ele%value(flexible$))) then
  if (dir == -1 .or. ie == branch%n_ele_track) then
    call out_io (s_fatal$, r_name, 'CONFUSION! PLEASE CONTACT DAVID SAGAN!')
    if (global_com%exit_on_error) call err_exit
    return
  endif
  ! If the position of the element just after a flexible patch is stale, 
  ! the patch geometry should be recomputed just to be on the safe side.
  if (branch%ele(ie+1)%bookkeeping_state%floor_position /= ok$) stale = .true.  
endif

!

if (.not. stale .and. ele%bookkeeping_state%floor_position /= stale$) return
if (ele%ix_ele == 0) ele%value(floor_set$) = true$

floor0 = ele%floor

if (dir == 1) then
  call ele_geometry (branch%ele(ie-1)%floor, ele, ele%floor, 1.0_rp, set_ok = .true.)
else
  call ele_geometry (branch%ele(ie+1)%floor, branch%ele(ie+1), ele%floor, -1.0_rp)
endif

stale = (.not. (ele%floor == floor0))

! target branch only needs to be recomputed if target branch index is greater than present branch.

if (ele%key == fork$ .or. ele%key == photon_fork$) then
  ib = nint(ele%value(ix_to_branch$))
  if (ib > n) then
    ix = nint(ele%value(ix_to_element$))
    lat%branch(ib)%ele(ix)%bookkeeping_state%floor_position = stale$
    lat%branch(ib)%param%bookkeeping_state%floor_position = stale$
  endif
endif

!

call multipass_chain(ele, ix_pass, n_links, chain_ele)

if (ix_pass > 0) then
  do k = ix_pass+1, n_links
    chain_ele(k)%ele%bookkeeping_state%floor_position = stale$
  enddo

  if (ele%slave_status == super_slave$) then
    do k = 1, ele%n_lord
      lord => pointer_to_lord(ele, k)
      if (lord%lord_status /= super_lord$) exit
      lord%bookkeeping_state%floor_position = stale$
      lord => pointer_to_lord(lord, 1)  ! multipass lord
      lord%bookkeeping_state%floor_position = stale$
    enddo
  else
    lord => pointer_to_lord(ele, 1)
    lord%bookkeeping_state%floor_position = stale$
  endif
endif

!

ele%bookkeeping_state%floor_position = ok$

end subroutine propagate_geometry

!---------------------------------------------------
! contains

! When girders themselves have girder lords, must do computation in order: Lord before slave.

recursive subroutine girder_lord_geometry (ele)
type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
integer i

!

do i = 1, ele%n_lord
  lord => pointer_to_lord(ele, i)
  if (lord%key /= girder$) cycle
  if (lord%bookkeeping_state%floor_position /= stale$) cycle
  call girder_lord_geometry (lord)
enddo

call ele_geometry (dummy, ele, ele%floor, set_ok = .true.)

end subroutine girder_lord_geometry

end subroutine lat_geometry

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine ele_geometry (floor0, ele, floor, len_scale, set_ok)
!
! Routine to calculate the global (floor) coordinates of an element given the
! global coordinates of the preceeding element. This is the same as the MAD convention.
!
! floor0 will correspond to the coordinates at the upstream end of the element
! and floor will correspond to the coordinates at the downstream end.
!
! Note: For floor_position element, floor is independent of floor0.
!
! Modules Needed:
!   use bmad
!
! Input:
!   floor0        -- Starting floor coordinates at upstream end.
!                      Not used for fiducial and girder elements.
!   ele           -- Ele_struct: Element to propagate the geometry through.
!   len_scale     -- Real(rp), optional: factor to scale the length of the element.
!                       1.0_rp => Output is geometry at end of element (default).
!                       0.5_rp => Output is geometry at center of element. [Cannot be used for crystals.]
!                      -1.0_rp => Used to propagate geometry in reverse.
!   set_ok        -- logical, optional: If present and True, set ele%bookkeeping_state%floor_position = T.
!
! Output:
!   floor       -- floor_position_struct: Floor position at downstream end.
!     %r(3)              -- X, Y, Z Floor position at end of element
!     %w(3,3)            -- W matrix corresponding to orientation angles
!     %theta, phi, %psi  -- Orientation angles 
!-

recursive subroutine ele_geometry (floor0, ele, floor, len_scale, set_ok)

use multipole_mod
use multipass_mod

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0, ele00, slave0, slave1, ele2
type (floor_position_struct) floor0, floor, floor_ref, floor_saved
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_pointer_struct), allocatable, target :: chain_ele(:)
type (lat_param_struct) param

real(rp), optional :: len_scale
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx), dtheta
real(rp) r0(3), w0_mat(3,3), rot_angle
real(rp) chord_len, angle, ang, leng, rho, len_factor
real(rp) theta, phi, psi, tlt, dz(3), z0(3), z_cross(3)
real(rp) :: s_ang, c_ang, w_mat(3,3), w_mat_inv(3,3), s_mat(3,3), r_vec(3), t_mat(3,3)

integer i, key, n_loc, ix_pass, n_links, ix_pole_max

logical err, calc_done, doit, finished
logical, optional :: set_ok

character(*), parameter :: r_name = 'ele_geometry'

! Custom geometry

if (logic_option(.false., set_ok)) ele%bookkeeping_state%floor_position = ok$

call ele_geometry_hook (floor0, ele, floor, finished, len_scale)
if (finished) return

! Init.

len_factor = ele%orientation * real_option(1.0_rp, len_scale)

if (ele%key /= girder$) then
  floor   = floor0
endif

theta   = floor0%theta
phi     = floor0%phi
psi     = floor0%psi
w_mat   = floor0%w


knl  = 0   ! initialize
tilt = 0  

leng = ele%value(l$) * len_factor

key = ele%key
if (key == sbend$ .and. (leng == 0 .or. ele%value(g$) == 0)) key = drift$

if (key == multipole$) then
  call multipole_ele_to_kt (ele, .true., ix_pole_max, knl, tilt)
endif

! Fiducial, floor_shift and girder elements.
! Note that fiducial, and girder elements are independent of floor0

if (key == fiducial$ .or. key == girder$ .or. key == floor_shift$) then

  ele0 => null()
  if (ele%component_name /= '') then
    call lat_ele_locator (ele%component_name, ele%branch%lat, eles, n_loc, err)
    if (n_loc /= 1) then
      call out_io (s_fatal$, r_name, 'ORIGIN_ELE: ' // ele%component_name,  &
                                     'FOR ELEMENT: ' // ele%name, &
                                     'IS NOT UNIQUE!')
      if (global_com%exit_on_error) call err_exit
      return
    endif

    ele0 => eles(1)%ele
  elseif (key == floor_shift$) then
    ele0 => pointer_to_next_ele(ele, -1)
  endif

  if (associated(ele0)) then
    calc_done = .false.
    select case (stream_ele_end(nint(ele%value(origin_ele_ref_pt$)), ele%orientation))
    case (upstream_end$)
      call ele_geometry (ele0%floor, ele0, floor_ref, -1.0_rp)

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
          if (ele0%value(b_param$) < 0) rot_angle = rot_angle - pi/2  ! Bragg
        case (mirror$, multilayer_mirror$)
          rot_angle = ele0%value(graze_angle$) - pi/2
        end select

        s_mat = w_mat_for_x_pitch (-rot_angle)
        w_mat = matmul (w_mat, s_mat)

        if (ele0%value(ref_tilt_tot$) /= 0) then
          t_mat = w_mat_for_tilt(ele0%value(ref_tilt_tot$))
          w_mat = matmul (t_mat, w_mat)
          t_mat(1,2) = -t_mat(1,2); t_mat(2,1) = -t_mat(2,1) ! form inverse
          w_mat = matmul (w_mat, t_mat)
        endif
        r0 = ele00%floor%r
        calc_done = .true.

      case default
        call ele_geometry (ele0%floor, ele0, floor_ref, -0.5_rp)
      end select

    case (downstream_end$)
      floor_ref = ele0%floor
    case default
      call out_io (s_fatal$, r_name, 'ORIGIN_ELE_REF_PT NOT VALID FOR ELEMENT: ' // ele%name)
      if (global_com%exit_on_error) call err_exit
    end select

    if (.not. calc_done) then
      w_mat = floor_ref%w
      r0 = floor_ref%r
    endif

  ! Fiducial with no origin ele: Use global origin.
  elseif (key == fiducial$) then
    call mat_make_unit(w_mat)
    r0 = 0

  ! Girder uses center of itself by default.
  else  ! must be a girder
    floor_ref = floor  ! Save
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
    ! floor0 for a fiducial element is meaningless. Instead use approximate floor_ref to try to
    ! have the angles come out without annoying factors of 2pi offsets
    floor_ref%theta = theta; floor_ref%phi = phi; floor_ref%psi = psi
    if (associated(ele0)) then
      floor_ref%theta = floor_ref%theta + ele0%floor%theta
      floor_ref%phi   = floor_ref%phi   + ele0%floor%phi
      floor_ref%psi   = floor_ref%psi   + ele0%floor%psi
    endif
    call update_floor_angles(floor, floor_ref)
  else
    call update_floor_angles(floor, floor0)
  endif

  if (any(floor%r /= floor_ref%r) .or. floor%theta /= floor_ref%theta .or. &
      floor%phi /= floor_ref%phi .or. floor%psi /= floor_ref%psi) then
    call set_ele_status_stale (ele, control_group$, .false.)
  endif
  return
endif

! General case where layout is not in the horizontal plane
! Note: 

if (((key == mirror$  .or. key == sbend$ .or. key == multilayer_mirror$) .and. &
         ele%value(ref_tilt_tot$) /= 0) .or. phi /= 0 .or. psi /= 0 .or. key == patch$ .or. &
         key == crystal$ .or. (key == multipole$ .and. knl(0) /= 0 .and. tilt(0) /= 0)) then

  !

  select case (key)

  ! sbend and multipole

  case (sbend$, multipole$)
    if (key == sbend$) then
      angle = leng * dble(ele%value(g$))
      tlt = ele%value(ref_tilt_tot$)
      rho = 1.0_dp / ele%value(g$)
      s_ang = sin(angle); c_ang = cos(angle)
      r_vec = [rho * (c_ang - 1), 0.0_dp, rho * s_ang]
    else
      angle = knl(0) * len_factor
      tlt = tilt(0)
      r_vec = 0
    endif
    ! By definition, positive angle is equivalent to negative x_pitch
    s_mat = w_mat_for_x_pitch(-angle)

    if (tlt /= 0) then
      t_mat = w_mat_for_tilt (tlt)

      r_vec = matmul (t_mat, r_vec)

      s_mat = matmul (t_mat, s_mat)
      t_mat(1,2) = -t_mat(1,2); t_mat(2,1) = -t_mat(2,1) ! form inverse
      s_mat = matmul (s_mat, t_mat)
    endif

    call floor_angles_to_w_mat (theta, phi, psi, w_mat)
    floor%r = floor%r + matmul(w_mat, r_vec)
    w_mat = matmul (w_mat, s_mat)

  ! mirror, multilayer_mirror, crystal

  case (mirror$, multilayer_mirror$, crystal$)
    
    if (ele%key == crystal$) then   ! Laue
      select case (nint(ele%value(ref_orbit_follows$)))
      case (bragg_diffracted$)
        angle = len_factor * (ele%value(bragg_angle_in$) + ele%value(bragg_angle_out$))
      case (forward_diffracted$, undiffracted$)
        angle = 0
      end select
      if (ele%value(b_param$) > 0) then
        ! %l_ref is with respect to the body coords
        r_vec = len_factor * ele%photon%material%l_ref
        if (len_factor > 0) then  ! Forward propagation -> Express r_vec in entrance coords.
          ang = len_factor * ele%value(bragg_angle_in$)
          r_vec = [cos(ang) * r_vec(1) - sin(ang) * r_vec(3), r_vec(2), sin(ang) * r_vec(1) + cos(ang) * r_vec(3)]
        else                      ! Express r_vec in exit coords
          ang = angle - len_factor * ele%value(bragg_angle_in$)
          r_vec = [cos(ang) * r_vec(1) - sin(ang) * r_vec(3), r_vec(2), sin(ang) * r_vec(1) + cos(ang) * r_vec(3)]
        endif
      else  ! Bragg
        r_vec = 0
      endif

    else   ! Non-crystal
      angle = 2 * len_factor * ele%value(graze_angle$)
      r_vec = 0
    endif

    ! By definition, positive angle is equivalent to negative x_pitch
    s_mat = w_mat_for_x_pitch (-angle)

    tlt = ele%value(ref_tilt_tot$)
    if (tlt /= 0) then
      t_mat = w_mat_for_tilt(tlt)
      r_vec = matmul(t_mat, r_vec)
      s_mat = matmul (t_mat, s_mat)
      t_mat(1,2) = -t_mat(1,2); t_mat(2,1) = -t_mat(2,1) ! form inverse
      s_mat = matmul (s_mat, t_mat)
    endif

    call floor_angles_to_w_mat (theta, phi, psi, w_mat)
    floor%r = floor%r + matmul(w_mat, r_vec)
    w_mat = matmul (w_mat, s_mat)

  ! patch

  case (patch$)

    ! Flexible bookkeeping to compute offsets and pitches. 
    ! Only needs to be done if element is part of a lattice (as opposed to, for example, a slice_slave) and,
    ! if part of a multipass region, only needs to be done if this is a first pass slave.

    if (is_true(ele%value(flexible$)) .and. ele%ix_ele > 0) then
      doit = .true.
      if (ele%lord_status == multipass_lord$) doit = .false.
      call multipass_chain (ele, ix_pass, n_links, chain_ele)
      if (ix_pass > 1) doit = .false.

      if (doit) then
        ele2 => pointer_to_next_ele(ele, 1)
        call multipass_chain (ele2, ix_pass, n_links, chain_ele)
        if (ix_pass > 0) then
          ele2 => chain_ele(1)%ele
        endif

        if (ele2%bookkeeping_state%floor_position == stale$) then
          call out_io (s_fatal$, r_name, 'ELEMENT AFTER FLEXIBLE PATCH: ' // trim(ele%name) // &
                                                          '  (' // trim(ele_loc_to_string(ele)) // ')', &
                                         'DOES NOT HAVE A WELL DEFINED POSITION')
          if (global_com%exit_on_error) call err_exit
        endif

        call ele_geometry (ele2%floor, ele2, floor_ref, -1.0_rp) 
        
        ele0 => pointer_to_next_ele(ele, -1)
        w_mat_inv = transpose(ele0%floor%w)
        w_mat = floor_ref%w
        w_mat = matmul(w_mat_inv, w_mat)
        call floor_w_mat_to_angles (w_mat, ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$))
        r_vec = matmul(w_mat_inv, floor_ref%r - ele0%floor%r)
        ele%value(x_offset$) = r_vec(1)
        ele%value(y_offset$) = r_vec(2)
        ele%value(z_offset$) = r_vec(3)
        w_mat_inv = transpose(w_mat)
        ele%value(l$) = w_mat_inv(3,1) * ele%value(x_offset$) + w_mat_inv(3,2) * ele%value(y_offset$) + &
                        w_mat_inv(3,3) * ele%value(z_offset$)
        if (ele_value_has_changed(ele, [l$], [bmad_com%significant_length], .true.)) then
          call set_ele_status_stale (ele, s_position_group$)
        endif

        ! Transfer offsets and pitches if patch is part of a multipass retion
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
      floor%r = floor%r - matmul(w_mat, r_vec)
    else
      floor%r = floor%r + matmul(w_mat, r_vec)
      call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$), s_mat)
      w_mat = matmul(w_mat, s_mat)
    endif

  ! everything else. Just a translation

  case default
    !call floor_angles_to_w_mat (theta, phi, psi, w_mat)
    floor%r = floor%r + w_mat(:,3) * leng

  end select 

! Simple case where the local reference frame stays in the horizontal plane.

else

  select case (key)
  case (sbend$)
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
  floor%r(1) = floor%r(1) + chord_len * sin(theta)
  floor%r(3) = floor%r(3) + chord_len * cos(theta)
  theta = theta - angle / 2

  call rotate_mat(w_mat, y_axis$, -angle)

endif

! Update floor angles

floor%w = w_mat 
call update_floor_angles(floor, floor0)

end subroutine ele_geometry

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine floor_angles_to_w_mat (theta, phi, psi, w_mat, w_mat_inv)
!
! Routine to construct the W matrix that specifies the orientation of an element
! in the global "floor" coordinates. See the Bmad manual for more details.
!
! Modules needed:
!   use geometry_mod
!
! Input:
!   theta -- Real(rp): Azimuth angle.
!   phi   -- Real(rp): Pitch angle.
!   psi   -- Real(rp): Roll angle.
!
! Output:
!   w_mat(3,3)     -- Real(rp), optional: Orientation matrix.
!   w_mat_inv(3,3) -- Real(rp), optional: Inverse Orientation matrix.
!-

subroutine floor_angles_to_w_mat (theta, phi, psi, w_mat, w_mat_inv)

real(rp), optional :: w_mat(3,3), w_mat_inv(3,3)
real(rp) theta, phi, psi
real(rp) s_the, c_the, s_phi, c_phi, s_psi, c_psi

!

s_the = sin(theta); c_the = cos(theta)
s_phi = sin(phi);   c_phi = cos(phi)
s_psi = sin(psi);   c_psi = cos(psi)

if (present(w_mat)) then
  w_mat(1,1) =  c_the * c_psi - s_the * s_phi * s_psi
  w_mat(1,2) = -c_the * s_psi - s_the * s_phi * c_psi
  w_mat(1,3) =  s_the * c_phi
  w_mat(2,1) =  c_phi * s_psi
  w_mat(2,2) =  c_phi * c_psi
  w_mat(2,3) =  s_phi 
  w_mat(3,1) = -s_the * c_psi - c_the * s_phi * s_psi
  w_mat(3,2) =  s_the * s_psi - c_the * s_phi * c_psi 
  w_mat(3,3) =  c_the * c_phi
endif

if (present(w_mat_inv)) then
  w_mat_inv(1,1) =  c_the * c_psi - s_the * s_phi * s_psi
  w_mat_inv(1,2) =  c_phi * s_psi 
  w_mat_inv(1,3) = -s_the * c_psi - c_the * s_phi * s_psi 
  w_mat_inv(2,1) = -c_the * s_psi - s_the * s_phi * c_psi 
  w_mat_inv(2,2) =  c_phi * c_psi
  w_mat_inv(2,3) =  s_the * s_psi - c_the * s_phi * c_psi
  w_mat_inv(3,1) =  s_the * c_phi
  w_mat_inv(3,2) =  s_phi
  w_mat_inv(3,3) =  c_the * c_phi
endif

end subroutine floor_angles_to_w_mat 

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine floor_w_mat_to_angles (w_mat, theta, phi, psi, floor0)
!
! Routine to construct the angles that define the orientation of an element
! in the global "floor" coordinates from the W matrix. See the Bmad manual for more details.
!
! Modules needed:
!   use geometry_mod
!
! Input:
!   w_mat(3,3) -- Real(rp): Orientation matrix.
!   floor0     -- floor_position_struct, optional: There are two solutions related by:
!                   [theta, phi, psi] & [pi+theta, pi-phi, pi+psi]
!                 If floor0 is present, choose the solution "nearest" the angles in floor0.
!
! Output:
!   theta -- Real(rp): Azimuth angle.
!   phi   -- Real(rp): Pitch angle.
!   psi   -- Real(rp): Roll angle.
!-

subroutine floor_w_mat_to_angles (w_mat, theta, phi, psi, floor0)

type (floor_position_struct), optional :: floor0
type (floor_position_struct) f0
real(rp) theta, phi, psi, w_mat(3,3)
real(rp) diff1(3), diff2(3)

! special degenerate case

if (abs(w_mat(1,3)) + abs(w_mat(3,3)) < 1d-12) then 
  ! Note: Only theta +/- psi is well defined here so this is rather arbitrary.
  if (present(floor0)) then
    theta = floor0%theta
  else
    theta = 0
  endif

  if (w_mat(2,3) > 0) then
    phi = pi/2
    psi = atan2(-w_mat(3,1), w_mat(1,1)) - theta
  else
    phi = -pi/2
    psi = atan2(w_mat(3,1), w_mat(1,1)) + theta
  endif

! normal case

else 
  if (present(floor0)) f0 = floor0      ! In case actual theta, phi, psi args are floor%theta, etc.
  theta = atan2 (w_mat(1,3), w_mat(3,3))
  phi = atan2 (w_mat(2,3), sqrt(w_mat(1,3)**2 + w_mat(3,3)**2))
  psi = atan2 (w_mat(2,1), w_mat(2,2))

  if (present(floor0)) then
    diff1 = [modulo2(theta-f0%theta, pi), modulo2(phi-f0%phi, pi), modulo2(psi-f0%psi, pi)]
    diff2 = [modulo2(pi+theta-f0%theta, pi), modulo2(pi-phi-f0%phi, pi), modulo2(pi+psi-f0%psi, pi)]
    if (sum(abs(diff2)) < sum(abs(diff1))) diff1 = diff2
    theta = diff1(1) + f0%theta
    phi   = diff1(2) + f0%phi
    psi   = diff1(3) + f0%psi
  else
    theta = theta - twopi * nint((theta ) / twopi)
  endif

endif

end subroutine floor_w_mat_to_angles 

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Function patch_flips_propagation_direction (x_pitch, y_pitch) result (is_flip)
!
! Routine to determine if the propagation direction is flipped in a patch.
! This is true if the tranformation matrix element S(3,3) = cos(x_pitch) * cos(y_pitch) 
! is negative.
!
! Module needed:
!   use geometry_mod
!
! Input:
!   x_pitch   -- Real(rp): Rotaion around y-axis
!   y_pitch   -- Real(rp): Rotation around x-axis.
!
! Output:
!   is_flip -- Logical: True if patch does a flip
!-

function patch_flips_propagation_direction (x_pitch, y_pitch) result (is_flip)

real(rp) x_pitch, y_pitch
logical is_flip

!

is_flip = (cos(x_pitch) * cos(y_pitch) < 0)

end function patch_flips_propagation_direction 

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Function coords_relative_to_floor (floor0, dr, theta, phi, psi) result (floor1)
!
! Starting from a given reference frame specified by its orientation and
! position in the global (floor) coordinates, and given a shift in position
! and angular orientation with respect to this reference frame, return the 
! resulting reference frame orientation and position.
!
! Also see: coords_floor_to_relative
!
! Module needed:
!   use geometry_mod
!
! Input:
!   floor0   -- floor_position_struct: Initial reference frame.
!   dr(3)    -- real(rp): (x, y, z) positional shift of the reference frame.
!   theta, phi, psi
!            -- real(rp), optional: Angular shift of the reference frame. See the 
!                 Bmad manual on the Global Coordinate system for more details.
!                 All angles must either be absent or present.
!
! Output:
!   floor1   -- floor_position_struct: Shifted reference frame.
!-

function coords_relative_to_floor (floor0, dr, theta, phi, psi) result (floor1)

type (floor_position_struct) floor0, floor1
real(rp) dr(3)
real(rp), optional :: theta, phi, psi
real(rp) w_mat(3,3), w0_mat(3,3)

!
floor1%r = matmul(floor0%w, dr) + floor0%r

if (present(theta)) then
  call floor_angles_to_w_mat (theta, phi, psi, w_mat)
  floor1%W = matmul(floor0%W, w_mat)
  call update_floor_angles (floor1, floor0)
endif

end function coords_relative_to_floor

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Function coords_floor_to_relative (floor0, global_position, calculate_angles, is_delta_position) result (local_position)
!
! Returns local floor position relative to floor0 given a global floor position.
! This is an essentially an inverse of routine coords_relative_to_floor.
!
! Input:
!   floor0            -- floor_position_struct: reference position
!   global_position   -- floor_position_struct: global position 
!   calculate_angles  -- logical, optional: calculate angles for local_position 
!                          Default: True.
!                          False returns local_position angles (%theta, %phi, %psi) = 0.
!   is_delta_position -- logical, optional: If True then treat global_position%r as a difference
!                           position in global space and only rotate the position but not shift it.
!                           Default: False.
!
! Output:
!  local_position -- floor_position_struct: position relative to floor0
!-

function coords_floor_to_relative (floor0, global_position, calculate_angles, is_delta_position) result (local_position)

type (floor_position_struct) floor0, global_position, local_position
real(rp) :: w0_mat_T(3,3), w_mat(3,3)
logical, optional :: calculate_angles, is_delta_position

! transpose
w0_mat_T = transpose(floor0%W)

! Solve for r_local = [x, y, z]_local
   
if (logic_option(.false., is_delta_position)) then
  local_position%r = matmul(w0_mat_T, global_position%r)
else
  local_position%r = matmul(w0_mat_T, global_position%r - floor0%r)
endif

local_position%w =  matmul(w0_mat_T, global_position%w)

! If angles are not needed, just return zeros; 
if (logic_option(.true., calculate_angles)) then
  call update_floor_angles(local_position)
else
  local_position%theta = 0
  local_position%phi = 0
  local_position%psi = 0
endif 


end function coords_floor_to_relative

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function coords_floor_to_local_curvilinear  (global_position, ele, status, w_mat) result(local_position)
!
! Given a position in global coordinates, return local curvilinear coordinates defined by ele.
!
! If the calculated longitudinal s places global_position outside of ele, the status argument
! will be set appropriately and local_position is not will be defined.
!
! Angular orientation is ignored.
!
! Module needed:
!   use geometry_mod
!
! Input:
!   global_position -- floor_position_struct: %r = [X, Y, Z] position in global coordinates
!   ele             -- ele_struct: element to find local coordinates of
!
! Result:
!   local_position  -- floor_position_struct: %r = [x, y, s] position in local curvilinear coordinates
!                        with s relative to start of element.
!   status          -- logical: inside$: s is inside ele
!                               upstream_end$: s is before element's entrance
!                               downstream_end$: s is beyond element's end
!   w_mat(3,3)      -- real(rp) (optional): W matrix at s, to transform vectors. 
!                                  v_global = w_mat.v_local
!                                  v_local = transpose(w_mat).v_global
!-  

function coords_floor_to_local_curvilinear (global_position, ele, status, w_mat) result(local_position)

use nr, only: zbrent

type (floor_position_struct) :: global_position, local_position
type (ele_struct)   :: ele
type (floor_position_struct) :: floor0
real(rp), optional :: w_mat(3,3)
real(rp) s_local, z_min, z_max, theta, dtheta, r1
integer :: status

! sbend case

if (ele%key == sbend$ .and. ele%value(g$) /= 0) then
  floor0 = ele%branch%ele(ele%ix_ele-1)%floor        ! Get floor0 from previous element
  local_position = coords_floor_to_relative (floor0, global_position)
  theta = atan(local_position%r(3) / ele%value(rho$))
  r1 = sqrt(local_position%r(3)**2 + ele%value(rho$)**2)
  dtheta = abs(atan(sqrt(local_position%r(1)**2 + local_position%r(2)**2) / r1) + 1d-6)
  z_min = (theta - dtheta) * ele%value(rho$)
  z_max = (theta + dtheta) * ele%value(rho$)
  local_position%r(3) = zbrent (delta_s_func, z_min, z_max, 1d-9)

! Straight line case

else
  if (ele%key == patch$ .and. ele%orientation == -1) then ! rotations at exit end
    floor0 = ele%branch%ele(ele%ix_ele-1)%floor        ! Get floor0 from previous element
    local_position = coords_floor_to_relative (floor0, global_position)
  else
    floor0 = ele%floor
    local_position = coords_floor_to_relative (floor0, global_position)
    local_position%r(3) = local_position%r(3) + ele%value(l$)
  endif
endif

! Inside or outside?

if (local_position%r(3) < 0) then
  status = upstream_end$
elseif (local_position%r(3) > ele%value(l$)) then
  status = downstream_end$
else
  status = inside$
endif

! Optionally return w_mat

if (present(w_mat)) then
  w_mat = local_position%W
endif

!-------------------------------------------------------------------------
contains

function delta_s_func (this_len) result (delta_s)

type (floor_position_struct) floor_at_s
real(rp), intent(in)  :: this_len
real(rp) :: delta_s

!

call ele_geometry(floor0, ele, floor_at_s, this_len /ele%value(l$))  
local_position = coords_floor_to_relative (floor_at_s, global_position, calculate_angles = .false.)
delta_s = local_position%r(3)

end function delta_s_func

end function coords_floor_to_local_curvilinear

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function coords_floor_to_curvilinear (floor_coords, ele0, ele1, status, w_mat) result(local_coords)
!
! Given a position in global "floor" coordinates, return local curvilinear (ie element) coordinates 
! for an appropriate element, ele1, near ele0. That is, the s-position of local_coords will be within
! the longitudinal extent of ele1.
!
! There may not be a corresponding local_coords. This may happen if the lattice is open or
! may happen near patch elemnts that have a non-zero x_pitch or y_pitch. 
! In these cases the status argument will be set accorrdingly.
!
! If there is no corresponding local_coords an "approximate" local_coords will
! be returned that may be used, for example, to do plotting but should not be used in calculations.
!
! Angular orientation of floor_coords is ignored.
!
! Module needed:
!   use geometry_mod
!
! Input:
!   floor_coords  -- floor_position_struct: %r = [X, Y, Z] position in global coordinates
!   ele0          -- ele_struct: Element to start the search at.
!
! Result:
!   local_coords  -- floor_position_struct: %r = [x, y, s] position in curvilinear coordinates
!                      with respect to ele1 with s relative to start the lattice branch.
!   ele1          -- ele_struct, pointer: Element that local_coords is with respect to.
!   status        -- logical: ok$             -> Local_coords found.
!                             patch_problem$  -> No solution due to a patch element.
!                             outside$        -> Outside of lattice ends (for open lattices).
!   w_mat(3,3)    -- real(rp) (optional): W matrix at s, to transform vectors from floor to local. 
!                      w_mat will only be well defined if status = ok$
!-  

function coords_floor_to_curvilinear (floor_coords, ele0, ele1, status, w_mat) result(local_coords)

type (floor_position_struct) floor_coords, local_coords
type (ele_struct), target :: ele0
type (ele_struct), pointer :: ele1
type (branch_struct), pointer :: branch

real(rp), optional :: w_mat(3,3)
real(rp) w_mat0(3,3), w_mat1(3,3)
integer status, last_direction, ix_ele, n_try
character(*), parameter :: r_name = 'coords_floor_to_curvilinear'

! Loop over neighboring elements until an encompassing one is found

branch => ele0%branch
ix_ele = ele0%ix_ele
last_direction = 0
n_try = 0

do
  ele1 => branch%ele(ix_ele)
  local_coords = coords_floor_to_local_curvilinear (floor_coords, ele1, status) 
  local_coords%r(3) = local_coords%r(3) + ele1%s_start
  if (status == inside$) exit
  n_try = n_try + 1

  ! No good so look for a problem and keep on searching if warranted

  if (status == upstream_end$) then
    if (n_try >= branch%n_ele_track .or. (ix_ele == 1 .and. branch%param%geometry == open$)) then
      status = outside$
      return
    endif

    if (last_direction == 1) then ! endless loop detected
      call no_convergence_calc (branch, ix_ele-1, ix_ele, status)
      return
    endif

    ! Try previous element
    ix_ele = ix_ele - 1
    if (ix_ele == 0) ix_ele = branch%n_ele_track
    last_direction = -1
    cycle

  else if (status == downstream_end$) then
    if (n_try >= branch%n_ele_track .or. (ix_ele == branch%n_ele_track .and. branch%param%geometry == open$)) then
      status = outside$
      return
    endif

    if (last_direction == -1) then ! endless loop detected
      call no_convergence_calc (branch, ix_ele, ix_ele+1, status)
      return
    endif

    ! Try next element
    ix_ele = ix_ele + 1
    if (ix_ele == branch%n_ele_track) ix_ele = 1
    last_direction = 1
    cycle

  else 
    ! This element contains position
    exit
  end if  
enddo

! Optionally return rotation matrix
if (present(w_mat) ) then
  w_mat = matmul(transpose(local_coords%W), floor_coords%W)
endif

status = ok$

!---------------------------------------------------------------------------
contains

subroutine no_convergence_calc (branch, ie1, ie2, status)

type (branch_struct), target :: branch
type (ele_struct), pointer :: elem1, elem2
type (floor_position_struct) pos1, pos2

integer ie1, ie2, status, stat

!

status = patch_problem$

elem1 => branch%ele(ie1)
elem2 => branch%ele(ie2)

if (elem1%key /= patch$ .and. elem2%key /= patch$) then
  call out_io (s_fatal$, r_name, 'CANNOT FIND CORRESPONDING LOCAL POSITION')
  if (global_com%exit_on_error) call err_exit
  return
endif

if (elem1%key == patch$ .and. ele1%orientation == 1) then
  pos1 = coords_floor_to_local_curvilinear (floor_coords, branch%ele(ie1-1), stat)
  pos1%r(3) = floor_coords%r(3) - elem1%value(z_offset$)
else
  pos1 = coords_floor_to_local_curvilinear (floor_coords, elem1, stat)
  pos1%r(3) = pos1%r(3) - elem1%value(l$)
endif

pos2 = coords_floor_to_local_curvilinear (floor_coords, elem2, stat)

if (abs(pos1%r(3)) < abs(pos2%r(3))) then
  local_coords = pos1
  ele1 => elem1
else
  local_coords = pos2
  ele1 => elem2
endif

local_coords%r(3) = elem1%s

end subroutine no_convergence_calc

end function coords_floor_to_curvilinear

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function coords_local_curvilinear_to_floor (local_position, ele, 
!                              in_ele_frame, w_mat, calculate_angles) result (global_position)
!
! Given a position local to ele, return global floor coordinates.
!
! Input:
!   local_position  -- floor_position_struct: Floor position in local curvilinear coordinates,
!                                             where %r = [x, y, s_local], 
!   ele             -- ele_struct: element that local_position coordinates are relative to.
!   in_ele_frame    -- logical, optional :: local_position is in ele body frame and includes misalignments.
!                               Default: .false.
!
! Result:
!   global_position -- floor_position_struct: Position in global coordinates.
!                       %r and %w
!   w_mat(3,3)      -- real(rp), optional: W matrix at s, to transform vectors. 
!                                  v_global = w_mat . v_local
!                                  v_local = transpose(w_mat) . v_global
!   
!   calculate_angles  -- logical, optional: calculate angles for global_position 
!                          Default: True.
!                          False returns local_position angles (%theta, %phi, %psi) = 0.
!-  

function coords_local_curvilinear_to_floor (local_position, ele, &
                              in_ele_frame, w_mat, calculate_angles) result (global_position)

type (floor_position_struct) :: local_position, global_position, p
type (ele_struct)  ele
real(rp) :: L_save
real(rp) :: w_mat_local(3,3), L_vec(3), S_mat(3,3), s
real(rp), optional :: w_mat(3,3)
logical, optional :: in_ele_frame
logical, optional :: calculate_angles

! Set x and y for floor offset 

p = local_position
 
if (logic_option(.false., in_ele_frame)) then
   ! General geometry with possible misalignments
   p = coords_element_frame_to_local(p, ele, w_mat = S_mat)
   
elseif (ele%key == sbend$) then
  ! Curved geometry, no misalignments. Get relative to ele's end
  s = p%r(3)
  p%r(3) = 0
  p = bend_shift(p, ele%value(g$), ele%value(L$) - s, w_mat = S_mat, tilt=ele%value(ref_tilt_tot$) )
  
else
   ! Element has Cartesian geometry, no misalignments. 
   ! position is relative to ele's exit: 
  p%r(3) = p%r(3) - ele%value(L$)
  call mat_make_unit(S_mat)
endif 

! Get global floor coordinates
global_position%r = matmul(ele%floor%w, p%r) + ele%floor%r
global_position%w = matmul(ele%floor%w, p%w)

! If angles are not needed, just return zeros; 
if (logic_option(.true., calculate_angles)) then
  call update_floor_angles(global_position)
else
  global_position%theta = 0
  global_position%phi = 0
  global_position%psi = 0
endif 

! Optionally return w_mat used in these transformations
if (present(w_mat) ) then
  w_mat = global_position%w
endif 

end function coords_local_curvilinear_to_floor

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function coords_element_frame_to_local(body_position, ele, w_mat, calculate_angles) result (local_position)
!
! Returns Cartesian coordinates relative to ele%floor (end of element).
!
! Input:
!   body_position   -- floor_position_struct: element body frame coordinates.
!     %r                [x, y, s] position with s = Position from beginning of element .
!   ele             -- ele_struct: element that local_position coordinates are relative to.
!   calculate_angles  -- logical, optional: calculate angles for local_position 
!                          Default: True.
!                          False returns local_position angles (%theta, %phi, %psi) = 0.
!  
! Output         
!   local_position  -- floor_position_struct: Cartesian coordinates relative to end of the element (ele%floor).
!
!   w_mat(3,3)      -- real(rp), optional: W matrix at to transform vectors. 
!                                  v_local     = w_mat . v_ele_frame
!                                  v_ele_frame = transpose(w_mat) . v_local
!-

function coords_element_frame_to_local(body_position, ele, w_mat, calculate_angles) result(local_position)

type (floor_position_struct) :: body_position, local_position
type (ele_struct) :: ele
real(rp) :: s, g, theta, Sb(3,3), Lb(3)
real(rp) ::  L_mis(3), S_mis(3,3) , S_mat0(3,3)
real(rp), optional :: w_mat(3,3)
logical, optional :: calculate_angles
!

local_position = body_position
s = body_position%r(3)

if (ele%key == sbend$) then
  ! Get coords relative to center
  theta = ele%value(g$)*s

  ! In ele frame at s relative to beginning. 
  ! Center coordinates at s, and move to center frame by ds = L/2 -s
  local_position%r(3) = 0
  local_position = bend_shift(local_position, ele%value(g$), ele%value(L$)/2 -s)
  
  ! Put into tilted frame
  call rotate_vec(local_position%r, z_axis$, ele%value(ref_tilt_tot$))
  call rotate_mat(local_position%W, z_axis$, ele%value(ref_tilt_tot$))

  ! Misalign
  call ele_misalignment_L_S_calc(ele, L_mis, S_mis)
  local_position%r = matmul(s_mis, local_position%r) + L_mis
  local_position%W = matmul(s_mis, local_position%W)

  ! Transform from center frame to element end frame
  local_position = bend_shift(local_position, ele%value(g$), ele%value(L$)/2, w_mat = Sb, tilt=ele%value(ref_tilt_tot$))
 
  if (present(w_mat)) then
    ! Initial rotation to ele's center frame and the tilt
    S_mat0 = w_mat_for_x_pitch (-(theta-ele%value(angle$)/2))
    call rotate_mat(S_mat0, z_axis$, ele%value(ref_tilt_tot$))
    w_mat = matmul(S_mis, S_mat0)
    w_mat = matmul(Sb, w_mat)
  endif

else
  if (ele%value(x_pitch_tot$) /= 0 .or. ele%value(y_pitch_tot$) /= 0 .or. ele%value(tilt_tot$) /= 0) then
    ! need to rotate element frame to align with local frame. 
    call floor_angles_to_w_mat (ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$), S_mis)
    local_position%r(3) = local_position%r(3) - ele%value(L$)/2 ! Angle offsets are relative to the center of the element
    local_position%r = matmul(S_mis, local_position%r)
    local_position%w = matmul(S_mis, local_position%w)
    local_position%r(3)  = local_position%r(3) - ele%value(L$)/2 ! Shift relative to end of element for output
  else
    ! Just shift relative to end
    call mat_make_unit(S_mis)
    local_position%r(3) = local_position%r(3)  - ele%value(L$)
  endif
  
  ! Add offsets
  local_position%r = local_position%r + [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)]
  
  if (present(w_mat)) w_mat = S_mis
endif

! If angles are not needed, just return zeros; 
if (logic_option(.true., calculate_angles)) then
  call update_floor_angles(local_position)
else
  local_position%theta = 0
  local_position%phi = 0
  local_position%psi = 0
endif 

end function coords_element_frame_to_local

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function coords_curvilinear_to_floor (xys, branch, err_flag) result (global)
!
! Routine to find the global position of a local (x, y, s) position.
!
! Input:
!   xys(3)      -- real(rp): (x, y, s) position vector.
!   branch      -- branch_struct: Lattice branch that defines the local reference coordinates.
!
! Output:
!   global      -- floor_position_struct: Global floor position corresponding to (x, y, s)
!   err_flag    -- logical: Set True if global floor position cannot be computed.
!-

function coords_curvilinear_to_floor (xys, branch, err_flag) result (global)

type (branch_struct), target :: branch
type (ele_struct), pointer :: ele
type (floor_position_struct) global, local

real(rp) xys(3)
integer ix_ele
logical err_flag

!

ix_ele = element_at_s (branch%lat, xys(3), .true., branch%ix_branch, err_flag)
if (err_flag) return
ele => branch%ele(ix_ele)

local%r = [xys(1), xys(2), xys(3) - ele%s_start]
global = coords_local_curvilinear_to_floor (local, ele)

end function coords_curvilinear_to_floor

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function w_mat_for_x_pitch (x_pitch, return_inverse) result (w_mat)
! 
! Routine to return the transformation matrix for an x_pitch.
!
! Module needed:
!   use geometry_mod
!
! Input:
!   x_pitch        -- real(rp): pitch angle
!   return_inverse -- logical, optional: If True, return the inverse matrix. Default is False.
! Output:
!   w_mat(3,3)     -- real(rp): Transformation matrix.
!-   

Function w_mat_for_x_pitch (x_pitch, return_inverse) result (w_mat)

real(rp) x_pitch, c_ang, s_ang
real(rp) :: w_mat(3,3)
logical, optional :: return_inverse

! An x_pitch corresponds to a rotation around the y axis.

c_ang = cos(x_pitch); s_ang = sin(x_pitch)

if (logic_option(.false., return_inverse)) then
  w_mat(1,:) = [ c_ang, 0.0_rp,  -s_ang]
  w_mat(2,:) = [0.0_rp, 1.0_rp,  0.0_rp]
  w_mat(3,:) = [ s_ang, 0.0_rp,   c_ang]
else
  w_mat(1,:) = [ c_ang, 0.0_rp,   s_ang]
  w_mat(2,:) = [0.0_rp, 1.0_rp,  0.0_rp]
  w_mat(3,:) = [-s_ang, 0.0_rp,   c_ang]
endif

end function w_mat_for_x_pitch

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function w_mat_for_y_pitch (y_pitch, return_inverse) result (w_mat)
! 
! Routine to return the transformation matrix for an y_pitch.
!
! Module needed:
!   use geometry_mod
!
! Input:
!   y_pitch        -- real(rp): pitch angle
!   return_inverse -- logical, optional: If True, return the inverse matrix. Default is False.
!
! Output:
!   w_mat(3,3)     -- real(rp): Transformation matrix.
!-   

Function w_mat_for_y_pitch (y_pitch, return_inverse) result (w_mat)

real(rp) y_pitch, c_ang, s_ang
real(rp) :: w_mat(3,3)
logical, optional :: return_inverse

! An y_pitch corresponds to a rotation around the y axis.

c_ang = cos(y_pitch); s_ang = sin(y_pitch)

if (logic_option(.false., return_inverse)) then
  w_mat(1,:) = [1.0_rp,  0.0_rp, 0.0_rp]
  w_mat(2,:) = [0.0_rp,  c_ang,  -s_ang]
  w_mat(3,:) = [0.0_rp,  s_ang,   c_ang]
else
  w_mat(1,:) = [1.0_rp,  0.0_rp, 0.0_rp]
  w_mat(2,:) = [0.0_rp,  c_ang,   s_ang]
  w_mat(3,:) = [0.0_rp, -s_ang,   c_ang]
endif

end function w_mat_for_y_pitch

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function w_mat_for_tilt (tilt, return_inverse) result (w_mat)
! 
! Routine to return the transformation matrix for an tilt.
!
! Module needed:
!   use geometry_mod
!
! Input:
!   tilt           -- real(rp): pitch angle
!   return_inverse -- logical, optional: If True, return the inverse matrix. Default is False.
!
! Output:
!   w_mat(3,3)     -- real(rp): Transformation matrix.
!-   

Function w_mat_for_tilt (tilt, return_inverse) result (w_mat)

real(rp) tilt, c_ang, s_ang
real(rp) :: w_mat(3,3)
logical, optional :: return_inverse

! An tilt corresponds to a rotation around the y axis.

c_ang = cos(tilt); s_ang = sin(tilt)

if (logic_option(.false., return_inverse)) then
  w_mat(1,:) = [ c_ang,  s_ang,  0.0_dp ]
  w_mat(2,:) = [-s_ang,  c_ang,  0.0_dp ]
  w_mat(3,:) = [0.0_dp,  0.0_dp, 1.0_dp ]
else
  w_mat(1,:) = [c_ang,  -s_ang,  0.0_dp ]
  w_mat(2,:) = [s_ang,   c_ang,  0.0_dp ]
  w_mat(3,:) = [0.0_dp,  0.0_dp, 1.0_dp ]
endif

end function w_mat_for_tilt

!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!+
! Subroutine ele_misalignment_L_S_calc (ele, L_mis, S_mis)
! 
! Calculates transformation vector L_mis and matrix S_mis due to misalignments for an ele
! Used to transform coordinates and vectors relative to the center of the element
!
! Module needed:
!   use geometry_mod
!
! Input:
!   ele       -- real(rp): Element
!
! Output:
!   L_mis(3)  -- real(rp): Misalignment vector relative to center of element
!   S_mis(3)  -- real(rp): Misalignment matrix relative to center of element
!
!-

subroutine ele_misalignment_L_S_calc (ele, L_mis, S_mis)

type(ele_struct) :: ele 
real(rp) :: Lc(3), Sb(3,3), s0
real(rp) :: L_mis(3), S_mis(3,3)

!

L_mis = [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)]

select case(ele%key)
case(sbend$)
  ! L_mis at ele center:
  ! L_mis = L_offsets + [Rz(roll) - 1] . Rz(tilt) . Ry(bend_angle/2) . rho . [cos(bend_angle/2) -1, 0, sin(bend_angle/2)]
  if (ele%value(roll_tot$) /= 0) then
    Lc = ele%value(rho$) * [cos(ele%value(angle$)/2) - 1, 0.0_rp, sin(ele%value(angle$)/2)]
    call rotate_vec(Lc, y_axis$, ele%value(angle$)/2)  ! rotate to entrance about y axis by half angle 
    call rotate_vec(Lc, z_axis$, ele%value(ref_tilt_tot$)) ! rotate about z axis by tilt
    L_mis = L_mis - Lc
    call rotate_vec(Lc, z_axis$, ele%value(roll_tot$))        ! rotate about z axis for roll
    L_mis = L_mis + Lc
  endif

  ! S_mis at ele center
  call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(roll_tot$), s_mis)
  
case default
    call floor_angles_to_w_mat (ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$), S_mis)
end select

end subroutine ele_misalignment_L_S_calc

!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!+
! Function bend_shift(position1, g, delta_s, w_mat, tilt) result(position2)
!
! Function to shift frame of reference within a bend with curvature g and tilt.
! Note: position2%theta, %phi, and %psi are not calculated. 
! 
! Module needed:
!   use geometry_mod
!
! Input:
!   position1    -- floor_position_struct: Position of particle in frame 1 coordinates (Caretesian).
!   g            -- real(rp): curvature (1/rho)
!   delta_s      -- real(rp): relative s-position of frame 2 relative to frame 1
!   tilt         -- real(rp), optional: tilt. Default: 0
!
! Result:
!   position2    -- floor_position_struct: Coordinates relative to frame 2
!   w_mat(3,3)   -- real(rp), optional: W matrix used in the transformation   
!-

function bend_shift (position1, g, delta_s, w_mat, tilt) result(position2)

type (floor_position_struct) :: position1, position2
real(rp) :: g, delta_s, S_mat(3,3), L_vec(3), tlt, angle
real(rp), optional :: w_mat(3,3), tilt

!

angle = delta_s * g

if (angle == 0) then
  position2 = position1
  position2%r(3) = position2%r(3) - delta_s
  if (present(w_mat)) call mat_make_unit(w_mat)
  return
endif

!

tlt = real_option(0.0_rp, tilt)
call mat_make_unit(S_mat)

if (tlt /= 0) then
  call rotate_mat(S_mat, z_axis$, -tilt)
  call rotate_mat(S_mat, y_axis$,  angle)
  call rotate_mat(S_mat, z_axis$,  tilt)
else
  call rotate_mat(S_mat, y_axis$, angle)
endif

L_vec = [cos(angle)-1, 0.0_rp, -sin(angle)]/g
if (present(tilt)) call rotate_vec(L_vec, z_axis$, tilt)

position2%r = matmul(S_mat, position1%r) + L_vec
position2%w = matmul(S_mat, position1%w)

if (present(w_mat)) w_mat = s_mat

end function bend_shift

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Helper routine to calculate floor angles from its W matrix
!-

subroutine update_floor_angles (floor, floor0)

type(floor_position_struct) :: floor
type(floor_position_struct), optional :: floor0

!

if (present(floor0)) then
  call floor_w_mat_to_angles (floor%W, floor%theta, floor%phi, floor%psi, floor0)
else
  call floor_w_mat_to_angles (floor%W, floor%theta, floor%phi, floor%psi)
endif

end subroutine update_floor_angles

end module
