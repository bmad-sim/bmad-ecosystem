module lat_geometry_mod

use bmad_struct
use bmad_interface
use lat_ele_loc_mod
use nr, only: zbrent

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
!       %a                -- X position at end of element
!       %b                -- Y position at end of element
!       %z                -- Z position at end of element
!       %theta            -- Orientation angle at end of element in X-Z plane
!       %phi              -- Elevation angle.
!       %psi              -- Roll angle.
!-

subroutine lat_geometry (lat)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, lord, slave, b_ele
type (branch_struct), pointer :: branch

integer i, n, ix2, ie, ib
logical stale, stale_lord

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

  ! Transfer info from branch element if that element exists.

  if (branch%ix_from_branch > -1 .and. (stale .or. branch%ele(0)%bookkeeping_state%floor_position == stale$)) then
    b_ele => pointer_to_ele (lat, branch%ix_from_ele, branch%ix_from_branch)
    branch%ele(0)%floor = b_ele%floor

    if (nint(b_ele%value(direction$)) == -1) then
      branch%ele(0)%floor%theta = modulo2(branch%ele(0)%floor%theta + pi, pi)
      branch%ele(0)%floor%phi   = -branch%ele(0)%floor%phi
      branch%ele(0)%floor%psi   = -branch%ele(0)%floor%psi
    endif
    stale = .true.
  endif

  if (branch%ele(0)%bookkeeping_state%floor_position == stale$) branch%ele(0)%bookkeeping_state%floor_position = ok$

  do i = 1, branch%n_ele_track
    ele => branch%ele(i)
    if (.not. stale .and. ele%bookkeeping_state%floor_position /= stale$) cycle
    call ele_geometry (branch%ele(i-1)%floor, ele, ele%floor)
    stale = .true.
    if (ele%key == branch$ .or. ele%key == photon_branch$) then
      ib = nint(ele%value(ix_branch_to$))
      lat%branch(ib)%ele(0)%bookkeeping_state%floor_position = stale$
    endif
    if (ele%n_lord > 0) then
      call set_lords_status_stale (ele, lat, floor_position_group$)
      stale_lord = .true.
    endif
    ele%bookkeeping_state%floor_position = ok$
  enddo

enddo

! put info in super_lords and multipass_lords

lat%lord_state%floor_position = ok$
lat%param%bookkeeping_state%floor_position = ok$

if (.not. stale_lord) return

do i = lat%n_ele_track+1, lat%n_ele_max  
  lord => lat%ele(i)
  if (.not. bmad_com%auto_bookkeeper .and. lord%bookkeeping_state%floor_position /= stale$) cycle
  if (lord%n_slave == 0) cycle

  select case (lord%lord_status)
  case (super_lord$)
    slave => pointer_to_slave(lord, lord%n_slave) ! Last slave is at exit end.
    lord%floor = slave%floor
  case (multipass_lord$)
    slave => pointer_to_slave(lord, 1)
    lord%floor = slave%floor
  end select

  lord%bookkeeping_state%floor_position = ok$
enddo

end subroutine lat_geometry

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine ele_geometry (floor0, ele, floor)
!
! Subroutine to calculate the global (floor) coordinates of an element given the
! global coordinates of the preceeding element. This is the same as the MAD convention.
!
! Modules Needed:
!   use bmad
!
! Input:
!   floor0 -- Starting floor coordinates.
!   ele    -- Ele_struct: Element to propagate the geometry through.
!
! Output:
!   floor -- floor_position_struct: Floor position at the exit end of ele.
!        %a                -- X position at end of element
!        %b                -- Y position at end of element
!        %z                -- Z position at end of element
!        %theta            -- Orientation angle at end of element in X-Z plane
!        %phi              -- Elevation angle.
!        %psi              -- Roll angle.
!-

subroutine ele_geometry (floor0, ele, floor)

use multipole_mod

implicit none

type (ele_struct) ele
type (floor_position_struct) floor0, floor
integer i, key

real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
real(dp) chord_len, angle, leng, rho
real(dp) pos(3), theta, phi, psi, tlt
real(dp), save :: old_theta = 100  ! garbage number
real(dp), save :: old_phi, old_psi
real(dp), save :: s_ang, c_ang
real(dp), save :: w_mat(3,3), s_mat(3,3), r_mat(3), t_mat(3,3)
real(dp), parameter :: twopi_dp = 2 * 3.14159265358979

logical has_nonzero_pole

! floor_position element: Floor position is already set in the element.

if (ele%key == floor_position$) return

! init
! old_theta is used to tell if we have to reconstruct the w_mat

pos   = [floor0%x, floor0%y, floor0%z ]
theta = floor0%theta
phi   = floor0%phi
psi   = floor0%psi

knl  = 0   ! initialize
tilt = 0  

leng = ele%value(l$)
key = ele%key
if (key == sbend$ .and. (leng == 0 .or. ele%value(g$) == 0)) key = drift$

if (key == multipole$) then
  call multipole_ele_to_kt (ele, positron$, .true., has_nonzero_pole, knl, tilt)
endif

! General case where layout is not in the horizontal plane

if (((key == mirror$  .or. key == crystal$ .or. key == sbend$) .and. ele%value(tilt_tot$) /= 0) .or. &
         phi /= 0 .or. psi /= 0 .or. key == patch$ .or. &
         (key == multipole$ .and. knl(0) /= 0 .and. tilt(0) /= 0)) then

  if (old_theta /= theta .or. old_phi /= phi .or. old_psi /= psi) then
    call floor_angles_to_w_mat (theta, phi, psi, w_mat)
  endif

  !

  select case (key)

  ! sbend and multipole

  case (sbend$, multipole$)
    if (key == sbend$) then
      angle = leng * dble(ele%value(g$))
      tlt = ele%value(tilt_tot$)
      rho = 1.0_dp / ele%value(g$)
      s_ang = sin(angle); c_ang = cos(angle)
      r_mat = [rho * (c_ang - 1), 0.0_dp, rho * s_ang ]
    else
      angle = knl(0)
      tlt = tilt(0)
      s_ang = sin(angle); c_ang = cos(angle)
      r_mat = 0
    endif

    s_mat(1,:) = [c_ang,  0.0_dp, -s_ang ]
    s_mat(2,:) = [0.0_dp, 1.0_dp,  0.0_dp ]
    s_mat(3,:) = [s_ang,  0.0_dp,  c_ang ]

    if (tlt /= 0) then
      s_ang = sin(tlt); c_ang = cos(tlt)
      t_mat(1,:) = [c_ang,  -s_ang,  0.0_dp ]
      t_mat(2,:) = [s_ang,   c_ang,  0.0_dp ]
      t_mat(3,:) = [0.0_dp,  0.0_dp, 1.0_dp ]

      r_mat = matmul (t_mat, r_mat)

      s_mat = matmul (t_mat, s_mat)
      t_mat(1,2) = -t_mat(1,2); t_mat(2,1) = -t_mat(2,1) ! form inverse
      s_mat = matmul (s_mat, t_mat)
    endif

    pos = pos + matmul(w_mat, r_mat)
    w_mat = matmul (w_mat, s_mat)

  ! mirror

  case (mirror$, crystal$)
    
    if (ele%key == mirror$) then
      angle = 2 * ele%value(graze_angle$)
    else
      angle = ele%value(graze_angle_in$) + ele%value(graze_angle_out$)
    endif

    tlt = ele%value(tilt_tot$)
    s_ang = sin(angle); c_ang = cos(angle)

    s_mat(1,:) = [c_ang,  0.0_dp, -s_ang ]
    s_mat(2,:) = [0.0_dp, 1.0_dp,  0.0_dp ]
    s_mat(3,:) = [s_ang,  0.0_dp,  c_ang ]

    if (tlt /= 0) then
      s_ang = sin(tlt); c_ang = cos(tlt)
      t_mat(1,:) = [c_ang,  -s_ang,  0.0_dp ]
      t_mat(2,:) = [s_ang,   c_ang,  0.0_dp ]
      t_mat(3,:) = [0.0_dp,  0.0_dp, 1.0_dp ]

      s_mat = matmul (t_mat, s_mat)
      t_mat(1,2) = -t_mat(1,2); t_mat(2,1) = -t_mat(2,1) ! form inverse
      s_mat = matmul (s_mat, t_mat)
    endif

    w_mat = matmul (w_mat, s_mat)

  ! patch

  case (patch$)

    if (ele%value(translate_after$) == 0) then
      r_mat = [ele%value(x_offset$), ele%value(y_offset$), ele%value(z_offset$) ]
      pos = pos + matmul(w_mat, r_mat)
    endif

    angle = ele%value(tilt$)
    if (angle /= 0) then
      s_ang = sin(angle); c_ang = cos(angle)
      s_mat(1,:) = [c_ang,  -s_ang,  0.0_dp ]
      s_mat(2,:) = [s_ang,   c_ang,  0.0_dp ]
      s_mat(3,:) = [0.0_dp,  0.0_dp, 1.0_dp ]
      w_mat = matmul(w_mat, s_mat)
    endif

    angle = ele%value(y_pitch_tot$)           ! 
    if (angle /= 0) then
      s_ang = sin(angle); c_ang = cos(angle)
      s_mat(1,:) = [1.0_dp,  0.0_dp, 0.0_dp ]
      s_mat(2,:) = [0.0_dp,  c_ang,  s_ang ]
      s_mat(3,:) = [0.0_dp, -s_ang,  c_ang ]
      w_mat = matmul(w_mat, s_mat)
    endif

    angle = ele%value(x_pitch_tot$)            ! x_pitch is negative MAD yrot
    if (angle /= 0) then
      s_ang = sin(angle); c_ang = cos(angle)
      s_mat(1,:) = [ c_ang,  0.0_dp, s_ang ]
      s_mat(2,:) = [ 0.0_dp, 1.0_dp, 0.0_dp ]
      s_mat(3,:) = [-s_ang,  0.0_dp, c_ang ]
      w_mat = matmul(w_mat, s_mat)
    endif
     
    if (ele%value(translate_after$) /= 0) then
      r_mat = [ele%value(x_offset$), ele%value(y_offset$), ele%value(z_offset$) ]
      pos = pos + matmul(w_mat, r_mat)
    endif

  ! everything else. Just a translation

  case default
    pos = pos + w_mat(:,3) * leng

  end select

  ! if there has been a rotation calculate new theta, phi, and psi

  if (key == sbend$ .or. key == patch$ .or. key == multipole$ .or. key == mirror$) then
    call floor_w_mat_to_angles (w_mat, floor0%theta, theta, phi, psi)
  endif

  old_theta = theta
  old_phi   = phi
  old_psi   = psi

! Simple case where the local reference frame stays in the horizontal plane.

else

  select case (key)
  case (sbend$)
    angle = leng * ele%value(g$)
    chord_len = 2 * leng * sin(angle/2) / angle
  case (multipole$)
    angle = knl(0)
    chord_len = 0
  case (mirror$)
    angle = 2 * ele%value(graze_angle$)
    chord_len = 0
  case (crystal$)
    angle = ele%value(graze_angle_in$) + ele%value(graze_angle_out$)
    chord_len = 0
  case default
    angle = 0
    chord_len = leng
  end select

  theta = theta - angle / 2
  pos(1) = pos(1) + chord_len * sin(theta)
  pos(3) = pos(3) + chord_len * cos(theta)
  theta = theta - angle / 2

endif

!

floor%x = pos(1)
floor%y = pos(2)
floor%z = pos(3)
floor%theta = theta
floor%phi   = phi
floor%psi   = psi

end subroutine ele_geometry

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine floor_angles_to_w_mat (theta, phi, psi, w_mat)
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
!   w_mat(3,3) -- Real(rp): Orientation matrix.
!-

subroutine floor_angles_to_w_mat (theta, phi, psi, w_mat)

implicit none

real(rp) theta, phi, psi, w_mat(3,3)
real(rp) s_the, c_the, s_phi, c_phi, s_psi, c_psi

!

s_the = sin(theta); c_the = cos(theta)
s_phi = sin(phi);   c_phi = cos(phi)
s_psi = sin(psi);   c_psi = cos(psi)
w_mat(1,1) =  c_the * c_psi - s_the * s_phi * s_psi
w_mat(1,2) = -c_the * s_psi - s_the * s_phi * c_psi
w_mat(1,3) =  s_the * c_phi
w_mat(2,1) =  c_phi * s_psi
w_mat(2,2) =  c_phi * c_psi
w_mat(2,3) =  s_phi 
w_mat(3,1) = -s_the * c_psi - c_the * s_phi * s_psi
w_mat(3,2) =  s_the * s_psi - c_the * s_phi * c_psi 
w_mat(3,3) =  c_the * c_phi

end subroutine floor_angles_to_w_mat 

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine floor_w_mat_to_angles (w_mat, theta0, theta, phi, psi)
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
!
! Output:
!   theta -- Real(rp): Azimuth angle.
!   phi   -- Real(rp): Pitch angle.
!   psi   -- Real(rp): Roll angle.
!-

subroutine floor_w_mat_to_angles (w_mat, theta0, theta, phi, psi)

implicit none

real(rp) theta0, theta, phi, psi, w_mat(3,3)

if (abs(w_mat(1,3)) + abs(w_mat(3,3)) < 1e-12) then ! special degenerate case
  ! Note: Only theta +/- psi is well defined here so this is rather arbitrary.
  theta = theta0  
  if (w_mat(2,3) > 0) then
    phi = pi/2
    psi = atan2(-w_mat(3,1), w_mat(1,1)) - theta
  else
    phi = -pi/2
    psi = atan2(w_mat(3,1), w_mat(1,1)) + theta
  endif
else  ! normal case
  theta = atan2 (w_mat(1,3), w_mat(3,3))
  theta = theta - twopi * nint((theta - theta0) / twopi)
  phi = atan2 (w_mat(2,3), sqrt(w_mat(1,3)**2 + w_mat(3,3)**2))
  psi = atan2 (w_mat(2,1), w_mat(2,2))
endif

end subroutine floor_w_mat_to_angles 

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine shift_reference_frame (floor0, dr, theta, phi, psi, floor1)
!
! Starting from a given reference frame specified by its orientation and
! position in the global (floor) coordinates, and given a shift in position
! and angular orientation with respect to this reference frame, return the 
! resulting reference frame orientation and position.
!
! Module needed:
!   use lat_geometry_mod
!
! Input:
!   floor0   -- floor_position_struct: Initial reference frame.
!   dr(3)    -- real(rp): (x, y, z) positional shift of the reference frame.
!   theta, phi, psi
!            -- real(rp): Angular shift of the reference frame. See the 
!                 Bmad manual on the Global Coordinate system for more details.
!
! Output:
!   floor1   -- floor_position_struct: Shifted reference frame.
!-

subroutine shift_reference_frame (floor0, dr, theta, phi, psi, floor1)

implicit none

type (floor_position_struct) floor0, floor1
real(rp) dr(3), theta, phi, psi
real(rp) r(3), w_mat(3,3), w0_mat(3,3)

!

call floor_angles_to_w_mat (floor0%theta, floor0%phi, floor0%psi, w0_mat)

r = matmul(w0_mat, dr) + [floor0%x, floor0%y, floor0%z]
floor1%x = r(1)
floor1%y = r(2)
floor1%z = r(3)

call floor_angles_to_w_mat (theta, phi, psi, w_mat)
w_mat = matmul(w0_mat, w_mat)
call floor_w_mat_to_angles (w_mat, 0.0_rp, floor1%theta, floor1%phi, floor1%psi)


end subroutine shift_reference_frame



!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Function floor_to_local (floor0, r_global) result (r_local)
!
! Returns local position relative to floor0 given a global position
!
! Input:
!   floor0        -- floor_position_struct: Coordinates to convert to.
!   r_global(3)   -- real(rp): [X, Y, Z] global position 
!
! Output:
!   r_local(3)    -- real(rp): [x, y, z] local position
!
!-

function floor_to_local (floor0, r_global) result (r_local)

implicit none

type (floor_position_struct) floor0
real(rp) :: r_local(3), r_global(3), w_mat(3,3)

!  Get w_mat
call floor_angles_to_w_mat (floor0%theta,floor0%phi, floor0%psi, w_mat)

!Solve for r_local = [x, y, z]_local
r_local = matmul(transpose(w_mat), r_global - [floor0%x, floor0%y, floor0%z])

end function floor_to_local


!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function position_in_local_frame  (r_global, floor0, ele) result(r_local)
!
! Given a point in global coordinates, return local curvilinear coordinates in ele
!   to floor0
!
! Module needed:
!   nr, only: zbrent
!
! Input:
!   r_global    -- real(rp) (3): [X, Y, Z] position in global coordinates
!   floor0      -- floor_position_struct: floor coordinates at the beginning of ele
!   ele         -- ele_struct: element to find local coordinates of
!
! Result:
!   r_local     -- real(rp) (3): [x, y, s] position in local curvilinear coordinates
!   err         -- logical: position is outside ele's bounds.  

function position_in_local_frame (r_global, floor0, ele, err) result(r_local)

implicit none

type (ele_struct)   :: ele
type (floor_position_struct) :: floor0, floor_at_s
real(rp) :: L_save, s_local, r_global(3), r_local(3) 
real(rp) :: w_mat(3,3)
logical  :: err

!

err = .false. 

! Check to see if position is within 0 < s < ele%value(L$)
r_local = floor_to_local (floor0, r_global)
if (r_local(3) < 0) then
  err = .true.
  return
endif
r_local = floor_to_local (ele%floor, r_global)
if (r_local(3) > 0) then
  err = .true.
  return
endif

! Save ele's L. We will vary this. 
L_save = ele%value(L$)

! Find s_local between 0 and L_save
s_local = zbrent(delta_s_in_ele_for_zbrent, 0.0_rp, L_save, 1d-9) 

! Restore ele's length
ele%value(L$) = L_save

! r_local was calculated in the zbrent function. Add in s_local
r_local = [r_local(1), r_local(2), s_local]

contains

  !
  ! function for zbrent to calculate s
  !
  !
  ! r_global = w_mat(s_local).[x_local, y_local, 0] 
  !             + [floor%x, floor%y, floor%z](s_local)
  ! Invert:
  ! => transpose(w_mat) . ( r_global - [floor%x, floor%y, floor%z](s_local) )
  !     == [x_local, y_local, 0]  when w_mat and floor% are calculated from correct s. 

  function delta_s_in_ele_for_zbrent (this_s)
  
  real(rp), intent(in)  :: this_s
  real(rp) :: delta_s_in_ele_for_zbrent

  !Vary L  
  ele%value(L$) = this_s  
  
  !Get floor_at_s
  call ele_geometry(floor0, ele, floor_at_s)
  
  !Get local coordinates   
  r_local = floor_to_local (floor_at_s, r_global)

  delta_s_in_ele_for_zbrent = r_local(3)
  
  end function  delta_s_in_ele_for_zbrent
end function position_in_local_frame


end module
