module lat_geometry_mod

use bmad_struct
use bmad_interface
use lat_ele_loc_mod
use rotation_3d_mod

contains

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine lat_geometry (lat)
!
! Subroutine to calculate the physical placement of all the elements in a lattice.
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
!       %r(3)              -- X, Y, Z Floor position at end of element
!       %theta, phi, %psi  -- Orientation angles 
!-

subroutine lat_geometry (lat)

use multipass_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, lord, slave, b_ele, ele2, ele0
type (branch_struct), pointer :: branch
type (ele_pointer_struct), allocatable :: chain_ele(:)
type (floor_position_struct) dummy

real(rp) w_mat(3,3), w_mat_inv(3,3), r_vec(3)

integer i, i2, n, ix2, ie, ib, ix_pass
logical stale, stale_lord

character(16), parameter :: r_name = 'lat_geometry'

!

stale_lord = .false.

do n = 0, ubound(lat%branch, 1)
  branch => lat%branch(n)

  if (bmad_com%auto_bookkeeper) then
    stale = .true.
  else
    if (branch%param%bookkeeping_state%floor_position /= stale$) cycle
    stale = .false.
  endif

  branch%param%bookkeeping_state%floor_position = ok$

  ! If there are fiducial elements then survey the fiducial regions

  do i = 1, branch%n_ele_track
    ele => branch%ele(i)
    if (ele%key /= fiducial$) cycle
    call ele_geometry (dummy, ele, ele%floor)

    do i2 = i+1, branch%n_ele_track
      ele2 => branch%ele(i2)
      if (ele2%key == patch$ .and. ele2%value(flexible$) /= 0) exit
      if (ele2%key == fiducial$) then
        call out_io (s_fatal$, r_name, 'FIDUCIAL ELEMENTS IN A BRANCH MUST BE SEPARATED BY A FLEXIBLE PATCH')
        if (global_com%exit_on_error) call err_exit
        exit
      endif
      call ele_geometry (branch%ele(i2-1)%floor, ele2, ele2%floor)
    enddo

    branch%ele(i-1)%floor = ele%floor  ! Save time

    do i2 = i-1, 1, -1
      ele2 => branch%ele(i2)
      if (ele2%key == patch$ .and. ele2%value(flexible$) /= 0) exit
      if (ele2%key == fiducial$) then
        call out_io (s_fatal$, r_name, 'FIDUCIAL ELEMENTS IN A BRANCH MUST BE SEPARATED BY A FLEXIBLE PATCH')
        if (global_com%exit_on_error) call err_exit
        exit
      endif
      call ele_geometry (ele2%floor, ele2, branch%ele(i2-1)%floor, -1.0_rp)
    enddo
  enddo

  ! Transfer info from the from_branch element if that element exists.

  if (branch%ix_from_branch > -1 .and. (stale .or. branch%ele(0)%bookkeeping_state%floor_position == stale$)) then
    b_ele => pointer_to_ele (lat, branch%ix_from_ele, branch%ix_from_branch)
    call ele_geometry (b_ele%floor, b_ele, branch%ele(0)%floor, treat_as_patch = .true.)
    stale = .true.
  endif

  if (branch%ele(0)%bookkeeping_state%floor_position == stale$) branch%ele(0)%bookkeeping_state%floor_position = ok$

  do i = 1, branch%n_ele_track
    ele => branch%ele(i)
    if (.not. stale .and. ele%bookkeeping_state%floor_position /= stale$) cycle

    if (ele%key == patch$ .and. ele%value(flexible$) /= 0) then
      ele2 => branch%ele(i+1)
      call multipass_chain (ele2, ix_pass, chain_ele = chain_ele)
      if (ix_pass > 0) then
        ele2 => chain_ele(1)%ele
        if (ele2%ix_ele /= 0) ele2 => pointer_to_next_ele(ele2, -1)
        ele%floor = ele2%floor
      endif

      ele0 => branch%ele(i-1)
      call floor_angles_to_w_mat (ele0%floor%theta, ele0%floor%phi, ele0%floor%psi, w_mat_inv = w_mat_inv)
      call floor_angles_to_w_mat (ele%floor%theta, ele%floor%phi, ele%floor%psi, w_mat)
      w_mat = matmul(w_mat_inv, w_mat)
      call floor_w_mat_to_angles (w_mat, 0.0_rp, ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$))
      r_vec = matmul(w_mat_inv, ele%floor%r - ele0%floor%r)
      ele%value(x_offset$) = r_vec(1)
      ele%value(y_offset$) = r_vec(2)
      ele%value(z_offset$) = r_vec(3)
      stale = .false.
    else
      stale = .true.
    endif

    call ele_geometry (branch%ele(i-1)%floor, ele, ele%floor)

    ! target branch only needs to be recomputed if target branch index is greater than present branch.

    if (ele%key == fork$ .or. ele%key == photon_fork$) then
      ib = nint(ele%value(ix_to_branch$))
      if (ib > n) lat%branch(ib)%ele(0)%bookkeeping_state%floor_position = stale$
    endif

    if (ele%n_lord > 0) then
      call set_lords_status_stale (ele, floor_position_group$)
      stale_lord = .true.
    endif
  enddo

enddo

! put info in super_lords and multipass_lords

lat%lord_state%floor_position = ok$
lat%param%bookkeeping_state%floor_position = ok$

if (.not. stale_lord) return
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
    lord%bookkeeping_state%control = stale$
  end select

enddo

!---------------------------------------------------
contains

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

call ele_geometry (dummy, ele, ele%floor)

end subroutine girder_lord_geometry

end subroutine lat_geometry

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine ele_geometry (floor0, ele, floor, len_scale, treat_as_patch)
!
! Subroutine to calculate the global (floor) coordinates of an element given the
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
!   floor0          -- Starting floor coordinates at upstream end.
!                        Not used for fiducial and girder elements.
!   ele             -- Ele_struct: Element to propagate the geometry through.
!   len_scale       -- Real(rp), optional: factor to scale the length of the element.
!                         1.0_rp => Output is geometry at end of element (default).
!                         0.5_rp => Output is geometry at center of element. [Cannot be used for crystals.]
!                        -1.0_rp => Used to propagate geometry in reverse.
!   treat_as_patch  -- Logical, option: If present and True then treat the element
!                        like a patch element. This is used by branch and photon_branch
!                        elements for constructing the coordinates of the "to" lattice branch.
!
! Output:
!   floor       -- floor_position_struct: Floor position at downstream end.
!     %r(3)              -- X, Y, Z Floor position at end of element
!     %theta, phi, %psi  -- Orientation angles 
!-

recursive subroutine ele_geometry (floor0, ele, floor, len_scale, treat_as_patch)

use multipole_mod

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0, ele00, slave0, slave1
type (floor_position_struct) floor0, floor, floor_ref
type (ele_pointer_struct), allocatable :: eles(:)
type (lat_param_struct) param

real(rp), optional :: len_scale
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx), dtheta
real(rp) r0(3), w0_mat(3,3), rot_angle
real(rp) chord_len, angle, ang, leng, rho, len_factor
real(rp) theta, phi, psi, tlt, dz(3), z0(3), z_cross(3)
real(rp) :: s_ang, c_ang, w_mat(3,3), s_mat(3,3), r_vec(3), t_mat(3,3)

integer i, key, n_loc

logical has_nonzero_pole, err, calc_done
logical, optional :: treat_as_patch

character(16), parameter :: r_name = 'ele_geometry'

! Init

ele%bookkeeping_state%floor_position = ok$
len_factor = ele%orientation * real_option(1.0_rp, len_scale)

floor   = floor0
theta   = floor0%theta
phi     = floor0%phi
psi     = floor0%psi

knl  = 0   ! initialize
tilt = 0  

leng = ele%value(l$) * len_factor

key = ele%key
if (key == sbend$ .and. (leng == 0 .or. ele%value(g$) == 0)) key = drift$
if (logic_option(.false., treat_as_patch)) key = patch$

if (key == multipole$) then
  call multipole_ele_to_kt (ele, param, .true., has_nonzero_pole, knl, tilt)
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
        call floor_angles_to_w_mat (ele00%floor%theta, ele00%floor%phi, ele00%floor%psi, w_mat)

        select case (ele0%key)
        case (crystal$)
          if (ele0%value(tilt_corr$) /= 0) then
            call w_mat_for_tilt(t_mat, ele0%value(tilt_corr$))
            w_mat = matmul(w_mat, t_mat)
          endif
          rot_angle = ele0%value(bragg_angle_in$) 
          if (ele0%value(b_param$) < 0) rot_angle = rot_angle - pi/2  ! Bragg
        case (mirror$, multilayer_mirror$)
          rot_angle = ele0%value(graze_angle$) - pi/2
        end select

        call w_mat_for_x_pitch (s_mat, -rot_angle)
        w_mat = matmul (w_mat, s_mat)

        if (ele0%value(ref_tilt_tot$) /= 0) then
          call w_mat_for_tilt(t_mat, ele0%value(ref_tilt_tot$))
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
      call floor_angles_to_w_mat (floor_ref%theta, floor_ref%phi, floor_ref%psi, w_mat)
      r0 = floor_ref%r
    endif

  ! Fiducial with no origin ele: Use global origin.
  elseif (key == fiducial$) then
    call mat_make_unit(w_mat)
    r0 = 0

  ! Girder uses center of itself by default.
  else  ! must be a girder
    call find_element_ends (ele, slave0, slave1)
    r0 = (slave0%floor%r + slave1%floor%r) / 2
    ! 
    call floor_angles_to_w_mat (slave0%floor%theta, slave0%floor%phi, slave0%floor%psi, w0_mat)
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
  if (ele%key == floor_shift$) then
    r_vec = [ele%value(x_offset$), ele%value(y_offset$), ele%value(z_offset$)]
    theta = ele%value(x_pitch$);  phi = ele%value(y_pitch$);  psi = ele%value(tilt$)
  else
    r_vec = [ele%value(dx_origin$), ele%value(dy_origin$), ele%value(dz_origin$)]
    theta = ele%value(dtheta_origin$);  phi = ele%value(dphi_origin$); psi = ele%value(dpsi_origin$)
  endif

  floor%r = r0 + matmul(w_mat, r_vec)
  call floor_angles_to_w_mat (theta, phi, psi, s_mat)
  w_mat = matmul(w_mat, s_mat)
  call floor_w_mat_to_angles (w_mat, 0.0_rp, floor%theta, floor%phi, floor%psi, floor0)

  return
endif

! General case where layout is not in the horizontal plane
! Note: 

if (((key == mirror$  .or. key == sbend$ .or. key == multilayer_mirror$) .and. &
         ele%value(ref_tilt$) /= 0) .or. phi /= 0 .or. psi /= 0 .or. key == patch$ .or. &
         key == crystal$ .or. (key == multipole$ .and. knl(0) /= 0 .and. tilt(0) /= 0)) then

  call floor_angles_to_w_mat (theta, phi, psi, w_mat)

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
    call w_mat_for_x_pitch(s_mat, -angle)

    if (tlt /= 0) then
      call w_mat_for_tilt (t_mat, tlt)

      r_vec = matmul (t_mat, r_vec)

      s_mat = matmul (t_mat, s_mat)
      t_mat(1,2) = -t_mat(1,2); t_mat(2,1) = -t_mat(2,1) ! form inverse
      s_mat = matmul (s_mat, t_mat)
    endif

    floor%r = floor%r + matmul(w_mat, r_vec)
    w_mat = matmul (w_mat, s_mat)

    call floor_w_mat_to_angles (w_mat, 0.0_rp, theta, phi, psi, floor0)

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
    call w_mat_for_x_pitch (s_mat, -angle)

    tlt = ele%value(ref_tilt_tot$)
    if (tlt /= 0) then
      call w_mat_for_tilt(t_mat, tlt)
      r_vec = matmul(t_mat, r_vec)
      s_mat = matmul (t_mat, s_mat)
      t_mat(1,2) = -t_mat(1,2); t_mat(2,1) = -t_mat(2,1) ! form inverse
      s_mat = matmul (s_mat, t_mat)
    endif

    floor%r = floor%r + matmul(w_mat, r_vec)
    w_mat = matmul (w_mat, s_mat)

    call floor_w_mat_to_angles (w_mat, 0.0_rp, theta, phi, psi, floor0)

  ! patch

  case (patch$)

    r_vec = [ele%value(x_offset$), ele%value(y_offset$), ele%value(z_offset$)]

    if (len_factor < 0) then
      call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$), w_mat_inv = s_mat)
      w_mat = matmul(w_mat, s_mat)
      floor%r = floor%r - matmul(w_mat, r_vec)
    else
      floor%r = floor%r + matmul(w_mat, r_vec)
      call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$), s_mat)
      w_mat = matmul(w_mat, s_mat)
    endif
     
    call floor_w_mat_to_angles (w_mat, 0.0_rp, theta, phi, psi, floor0)

  ! everything else. Just a translation

  case default
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

endif

!

floor%theta = theta
floor%phi   = phi
floor%psi   = psi

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
!   use lat_geometry_mod
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

implicit none

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
! Subroutine floor_w_mat_to_angles (w_mat, theta0, theta, phi, psi, floor0)
!
! Routine to construct the angles that define the orientation of an element
! in the global "floor" coordinates from the W matrix. See the Bmad manual for more details.
!
! Modules needed:
!   use lat_geometry_mod
!
! Input:
!   w_mat(3,3) -- Real(rp): Orientation matrix.
!   theta0     -- Real(rp): Reference azimuth angle. The output theta will be in the range:
!                   [theta0 - pi, theta0 + pi]. Theta0 is used to keep track of the total
!                   winding angle. If you care, set theta0 to the old value of theta.
!                   If you don't care, set theta0 to 0.0_rp.
!   floor0     -- floor_position_struct, optional: There are two solutions related by:
!                   [theta, phi, psi] & [pi+theta, pi-phi, pi+psi]
!                 If floor0 is present, choose the solution "nearest" the angles in floor0.
!                 If floor0 is present then theta0 is ignored.
!
! Output:
!   theta -- Real(rp): Azimuth angle.
!   phi   -- Real(rp): Pitch angle.
!   psi   -- Real(rp): Roll angle.
!-

subroutine floor_w_mat_to_angles (w_mat, theta0, theta, phi, psi, floor0)

implicit none

type (floor_position_struct), optional :: floor0
type (floor_position_struct) f0
real(rp) theta0, theta, phi, psi, w_mat(3,3)
real(rp) diff1(3), diff2(3)

! special degenerate case

if (abs(w_mat(1,3)) + abs(w_mat(3,3)) < 1e-12) then 
  ! Note: Only theta +/- psi is well defined here so this is rather arbitrary.
  if (present(floor0)) then
    theta = floor0%theta
  else
    theta = theta0  
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
    theta = theta - twopi * nint((theta - theta0) / twopi)
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
!   use lat_geometry_mod
!
! Input:
!   x_pitch   -- Real(rp): Rotaion around y-axis
!   y_pitch   -- Real(rp): Rotation around x-axis.
!
! Output:
!   is_flip -- Logical: True if patch does a flip
!-

function patch_flips_propagation_direction (x_pitch, y_pitch) result (is_flip)

implicit none

real(rp) x_pitch, y_pitch
logical is_flip

!

is_flip = (cos(x_pitch) * cos(y_pitch) < 0)

end function patch_flips_propagation_direction 

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Function local_to_floor (floor0, dr, theta, phi, psi) result (floor1)
!
! Starting from a given reference frame specified by its orientation and
! position in the global (floor) coordinates, and given a shift in position
! and angular orientation with respect to this reference frame, return the 
! resulting reference frame orientation and position.
!
! Also see: floor_to_local
!
! Module needed:
!   use lat_geometry_mod
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

function local_to_floor (floor0, dr, theta, phi, psi) result (floor1)

implicit none

type (floor_position_struct) floor0, floor1
real(rp) dr(3)
real(rp), optional :: theta, phi, psi
real(rp) w_mat(3,3), w0_mat(3,3)

!

call floor_angles_to_w_mat (floor0%theta, floor0%phi, floor0%psi, w0_mat)

floor1%r = matmul(w0_mat, dr) + floor0%r

if (present(theta)) then
  call floor_angles_to_w_mat (theta, phi, psi, w_mat)
  w_mat = matmul(w0_mat, w_mat)
  call floor_w_mat_to_angles (w_mat, 0.0_rp, floor1%theta, floor1%phi, floor1%psi, floor0)
endif

end function local_to_floor

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Function floor_to_local (floor0, global_position, calculate_angles, is_delta_position) result (local_position)
!
! Returns local floor position relative to floor0 given a global floor position.
! This is an essentially an inverse of routine local_to_floor.
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

function floor_to_local (floor0, global_position, calculate_angles, is_delta_position) result (local_position)

implicit none

type (floor_position_struct) floor0, global_position, local_position
real(rp) :: w0_mat(3,3), w_mat(3,3)
logical, optional :: calculate_angles, is_delta_position

! Get w0_mat and invert

call floor_angles_to_w_mat (floor0%theta,floor0%phi, floor0%psi, w0_mat)
w0_mat = transpose(w0_mat)

! Solve for r_local = [x, y, z]_local
   
if (logic_option(.false., is_delta_position)) then
  local_position%r = matmul(w0_mat, global_position%r)
else
  local_position%r = matmul(w0_mat, global_position%r - floor0%r)
endif

! If angles are not needed, just return zeros; 
if (.not. logic_option(.true., calculate_angles) ) then
  local_position%theta = 0
  local_position%phi = 0
  local_position%psi = 0
  return
endif 

call floor_angles_to_w_mat (global_position%theta, global_position%phi, global_position%psi, w_mat)
w_mat = matmul(w0_mat, w_mat) ! Remember that w0_mat is actually w0_mat^T
call floor_w_mat_to_angles (w_mat, 0.0_rp, local_position%theta, local_position%phi, local_position%psi)

end function floor_to_local

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function position_in_local_frame  (global_position, ele, status, w_mat) result(local_position)
!
! Given a position in global coordinates, return local curvilinear coordinates in ele
!   relative to floor0
!
! Module needed:
!   use lat_geometry_mod
!
! Input:
!   global_position -- floor_position_struct: [X, Y, Z] position in global coordinates
!   ele             -- ele_struct: element to find local coordinates of
!
! Result:
!   local_position  -- floor_position_struct: [x, y, s] position in local curvilinear coordinates
!   status          -- logical: inside$: s is inside ele
!                               upstream_end$: s is before element's entrance
!                               downstream_end$: s is beyond element's end
!   w_mat(3,3)      -- real(rp) (optional): W matrix at s, to transform vectors. 
!                                  v_global = w_mat.v_local
!                                  v_local = transpose(w_mat).v_global
!       
!-  

function position_in_local_frame (global_position, ele, status, w_mat) result(local_position)

use nr, only: zbrent

implicit none

type (floor_position_struct) :: global_position, local_position
type (ele_struct)   :: ele
type (floor_position_struct) :: floor0, floor_at_s
real(rp) :: L_save, s_local, r_global(3)
real(rp), optional :: w_mat(3,3)
integer :: status
logical  :: err

!

status = inside$

! Save ele's L. We will vary this. 
L_save = ele%value(L$)

if (associated (ele%branch) ) then
  ! Get floor0 from previous element
  floor0 = ele%branch%ele(ele%ix_ele-1)%floor
else
  ! ele is without a lat. Propagate backwards to get floor0
  ele%value(L$) = -L_save
  call ele_geometry(ele%floor, ele, floor0)
  ele%value(L$) = L_save
endif


! Check to see if position is within 0 < s < ele%value(L$)
local_position = floor_to_local (floor0, global_position)
if (local_position%r(3) < 0) then
  status = upstream_end$
  return
endif
local_position = floor_to_local (ele%floor, global_position)
if (local_position%r(3) > 0) then
  status = downstream_end$
  return
endif


! Find s_local between 0 and L_save
s_local = zbrent(delta_s_in_ele_for_zbrent, 0.0_rp, L_save, 1d-9) 

! Restore ele's length
ele%value(L$) = L_save

! r_local was calculated in the zbrent function. Add in s_local
local_position%r(3) = s_local

! Optionally return w_mat
if (present(w_mat) ) then
  call floor_angles_to_w_mat (local_position%theta, local_position%phi, local_position%psi, w_mat)
endif

!-------------------------------------------------------------------------
contains

!+
! function for zbrent to calculate s
!
! r_global = w_mat(s_local).[x_local, y_local, 0] + floor%r (s_local)
! Invert:
! => transpose(w_mat) . ( r_global - floor%r (s_local) )
!     == [x_local, y_local, 0]  when w_mat and floor% are calculated from correct s. 
!-

function delta_s_in_ele_for_zbrent (this_s)

real(rp), intent(in)  :: this_s
real(rp) :: delta_s_in_ele_for_zbrent

!Vary L  
ele%value(L$) = this_s  

!Get floor_at_s
call ele_geometry(floor0, ele, floor_at_s)

!Get local coordinates   
local_position = floor_to_local (floor_at_s, global_position, calculate_angles = .false.)

delta_s_in_ele_for_zbrent = local_position%r(3)

end function  delta_s_in_ele_for_zbrent
end function position_in_local_frame


!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function position_in_global_frame (local_position, ele, w_mat) result (global_position)
!
! Given a position local to ele, return global floor coordinates
!
! Input:
!   local_position  -- floor_position_struct: Floor position in local curvilinear coordinates.
!     %r(3)             -- Position from beginning of element.
!   ele             -- ele_struct: element that local_position coordinates are relative to.
!
! Result:
!   global_position -- floor_position_struct: Position in global coordinates.
!   w_mat(3,3)      -- real(rp), optional: W matrix at s, to transform vectors. 
!                                  v_global = w_mat . v_local
!                                  v_local = transpose(w_mat) . v_global
!       
!-  

function position_in_global_frame (local_position, ele, w_mat) result (global_position)

implicit none

type (floor_position_struct) :: local_position, global_position, floor, floor0
type (ele_struct)  ele
real(rp) :: L_save
real(rp) :: dr(3)
real(rp), optional :: w_mat(3,3)

! Set x and y for floor offset 

dr = local_position%r
 
if (ele%key == sbend$ .or. ele%key == rbend$) then
  ! Element has a curved geometry. Shorten ele
  L_save = ele%value(L$)
  ele%value(L$) = local_position%r(3)
  
  ! calculate floor from previous element
  if (associated (ele%branch) ) then
    ! Get floor0 from previous element
    floor0 = ele%branch%ele(ele%ix_ele-1)%floor
  else
    ! ele is without a lat. Propagate backwards to get floor0
    ele%value(L$) = -L_save
    call ele_geometry(ele%floor, ele, floor0)
    ele%value(L$) = L_save
  endif
    
  call ele_geometry(floor0, ele, floor)
  ! position is exactly at ele's exit now. 
  dr(3) = 0
  
  !Restore ele's length
  ele%value(L$) = L_save

else
   ! Element has Cartesian geometry. 
   floor = ele%floor
      
   ! position is relative to ele's exit: 
   dr(3) =  local_position%r(3) - ele%value(L$) 
endif 

! Get global floor coordinates
global_position = local_to_floor (floor, dr, local_position%theta, local_position%phi, local_position%psi)

! Optionally return w_mat
if (present(w_mat) ) then
  call floor_angles_to_w_mat (global_position%theta, global_position%phi,global_position%psi, w_mat)
endif 

end function position_in_global_frame

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine switch_local_positions (position0, ele0, ele_try, position1, ele1, ww_mat)
! 
! Subroutine to take a local position0 in ele0 and find a local position near ele_try.
! If this position is beyond the bounds of ele_try, neighboring elements will be 
! stepped to until a containing element is found.
! Optionally returns the ww_mat = W1^T.W0 matrix needed to rotate vectors:  
!      W0.v0 = W1.v1 => W1^T.W0.v0 = v1
!
! Input:
!   position0   -- floor_position_struct: local position in ele0
!   ele0        -- ele_struct: Element that position0 is local to
!   ele_try     -- ele_struct: Element to try to find a local position
!
! Output: 
!   position1   -- floor_position_struct: local position in ele1
!   ele1        --  ele_struct, pointer :: element that contains position1
!   ww_mat(3,3) -- real(rp), optional: W1^T.W0 matrix
!
!-
subroutine switch_local_positions (position0, ele0, ele_try, position1, ele1, ww_mat)

implicit none

type (floor_position_struct) :: position0, position1, global_position
type (ele_struct) :: ele0
type (ele_struct), target :: ele_try
type (ele_struct), pointer :: ele1
real(rp), optional :: ww_mat(3,3)
real(rp) :: w_mat0(3,3),  w_mat1(3,3)
integer :: ix_ele, status

character(30), parameter :: r_name = 'switch_local_positions'

!

! Make sure ele_try has a branch
if (.not. associated (ele_try%branch) ) then
      call out_io (s_fatal$, r_name, 'ELE_TRY HAS NO ASSOCIATED BRANCH')
      if (global_com%exit_on_error) call err_exit
endif

!
ix_ele = ele_try%ix_ele

! Get global position
global_position = position_in_global_frame (position0, ele0)

! Loop over neighboring elements until an encompassing one is found
do
  position1 = position_in_local_frame  (global_position, ele_try%branch%ele(ix_ele), status) 
  if (status == upstream_end$) then
    ! Try previous element
    ix_ele = ix_ele -1
    cycle
  else if (status == downstream_end$) then
    ! Try next element
    ix_ele = ix_ele + 1
    cycle
  else if (status == inside$) then
    ! This element contains position
    ele1 => ele_try%branch%ele(ix_ele)
    exit
  end if  
enddo

! Optionally return rotation matrix
if (present(ww_mat) ) then
  call floor_angles_to_w_mat (global_position%theta, global_position%phi,global_position%psi, w_mat0)
  call floor_angles_to_w_mat (position1%theta, position1%phi, position1%psi, w_mat1)
  ww_mat = matmul( transpose(w_mat1), w_mat0)
endif
  
end subroutine switch_local_positions

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine w_mat_for_x_pitch (w_mat, x_pitch)
! 
! Routine to return the transformation matrix for an x_pitch.
!
! Module needed:
!   use lat_geometry_mod
!
! Input:
!   x_pitch     -- real(rp): pitch angle
!
! Output:
!   w_mat(3,3)  -- real(rp): Transformation matrix.
!-   

Subroutine w_mat_for_x_pitch (w_mat, x_pitch)

implicit none

real(rp) w_mat(3,3), x_pitch, c_ang, s_ang

! An x_pitch corresponds to a rotation around the y axis.

c_ang = cos(x_pitch); s_ang = sin(x_pitch)

w_mat(1,:) = [ c_ang, 0.0_rp,   s_ang]
w_mat(2,:) = [0.0_rp, 1.0_rp,  0.0_rp]
w_mat(3,:) = [-s_ang, 0.0_rp,   c_ang]

end subroutine w_mat_for_x_pitch

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine w_mat_for_y_pitch (w_mat, y_pitch)
! 
! Routine to return the transformation matrix for an y_pitch.
!
! Module needed:
!   use lat_geometry_mod
!
! Input:
!   y_pitch     -- real(rp): pitch angle
!
! Output:
!   w_mat(3,3)  -- real(rp): Transformation matrix.
!-   

Subroutine w_mat_for_y_pitch (w_mat, y_pitch)

implicit none

real(rp) w_mat(3,3), y_pitch, c_ang, s_ang

! An y_pitch corresponds to a rotation around the y axis.

c_ang = cos(y_pitch); s_ang = sin(y_pitch)

w_mat(2,:) = [1.0_rp,  0.0_rp, 0.0_rp]
w_mat(1,:) = [0.0_rp,  c_ang,   s_ang]
w_mat(3,:) = [0.0_rp, -s_ang,   c_ang]

end subroutine w_mat_for_y_pitch

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine w_mat_for_tilt (w_mat, tilt)
! 
! Routine to return the transformation matrix for an tilt.
!
! Module needed:
!   use lat_geometry_mod
!
! Input:
!   tilt     -- real(rp): pitch angle
!
! Output:
!   w_mat(3,3)  -- real(rp): Transformation matrix.
!-   

Subroutine w_mat_for_tilt (w_mat, tilt)

implicit none

real(rp) w_mat(3,3), tilt, c_ang, s_ang

! An tilt corresponds to a rotation around the y axis.

c_ang = cos(tilt); s_ang = sin(tilt)

w_mat(1,:) = [c_ang,  -s_ang,  0.0_dp ]
w_mat(2,:) = [s_ang,   c_ang,  0.0_dp ]
w_mat(3,:) = [0.0_dp,  0.0_dp, 1.0_dp ]

end subroutine w_mat_for_tilt

end module
