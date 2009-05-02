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

#include "CESR_platform.inc"

subroutine lat_geometry (lat)

use bmad_struct
use bmad_interface, except_dummy => lat_geometry

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: line

integer i, n, ix2, ie

!

do i = 1, lat%n_ele_track
  call ele_geometry (lat%ele(i-1)%floor, lat%ele(i), lat%ele(i)%floor)
enddo

do n = 1, ubound(lat%branch, 1)
  line => lat%branch(n)

  call pointer_to_ele (lat, line%ix_from_branch, line%ix_from_ele, ele)
  line%ele(0)%floor = ele%floor

  if (ele%key == photon_branch$ .or. ele%key == photon_branch$) then
    if (nint(ele%value(direction$)) == -1) then
      line%ele(0)%floor%theta = -line%ele(0)%floor%theta
      line%ele(0)%floor%phi   = -line%ele(0)%floor%phi
    endif
  endif

  do i = 1, line%n_ele_track
    call ele_geometry (line%ele(i-1)%floor, line%ele(i), lat%ele(i)%floor)
  enddo

enddo

! put info in superimpose lords

do i = lat%n_ele_track+1, lat%n_ele_max  
  ele => lat%ele(i)
  if (ele%key == super_lord$) then
    ix2 = ele%ix2_slave
    ie = lat%control(ix2)%ix_slave
    ele%floor = lat%ele(ie)%floor
  endif
enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine ele_geometry (floor0, ele, floor)
!
! Subroutine to calculate the physical (floor) placement of an element given the
! placement of the preceeding element. This is the same as the MAD convention.
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

use bmad_struct
use bmad_interface, except => ele_geometry
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
real(dp), save :: s_ang, c_ang, s_the, c_the, s_phi, c_phi, s_psi, c_psi
real(dp), save :: w_mat(3,3), s_mat(3,3), r_mat(3), t_mat(3,3)
real(dp), parameter :: twopi_dp = 2 * 3.14159265358979

! init
! old_theta is used to tell if we have to reconstruct the w_mat

pos   = (/ floor0%x, floor0%y, floor0%z /)
theta = floor0%theta
phi   = floor0%phi
psi   = floor0%psi

knl  = 0   ! initialize
tilt = 0  

leng = ele%value(l$)
key = ele%key
if (key == sbend$ .and. (leng == 0 .or. ele%value(g$) == 0)) key = drift$

if (key == multipole$) then
  call multipole_ele_to_kt (ele, positron$, knl, tilt, .true.)
endif

! General case where layout is not in the horizontal plane

if (phi /= 0 .or. psi /= 0 .or. key == patch$ .or. &
             (key == multipole$ .and. knl(0) /= 0 .and. tilt(0) /= 0) .or. &
             (key == sbend$ .and. ele%value(tilt_tot$) /= 0)) then

  if (old_theta /= theta .or. old_phi /= phi .or. old_psi /= psi) then
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
      r_mat = (/ rho * (c_ang - 1), 0.0_dp, rho * s_ang /)
    else
      angle = knl(0)
      tlt = tilt(0)
      s_ang = sin(angle); c_ang = cos(angle)
      r_mat = 0
    endif

    s_mat(1,:) = (/ c_ang,  0.0_dp, -s_ang /)
    s_mat(2,:) = (/ 0.0_dp, 1.0_dp,  0.0_dp /)
    s_mat(3,:) = (/ s_ang,  0.0_dp,  c_ang /) 

    if (tlt /= 0) then
      s_ang = sin(tlt); c_ang = cos(tlt)
      t_mat(1,:) = (/ c_ang,  -s_ang,  0.0_dp /)
      t_mat(2,:) = (/ s_ang,   c_ang,  0.0_dp /)
      t_mat(3,:) = (/ 0.0_dp,  0.0_dp, 1.0_dp /)

      r_mat = matmul (t_mat, r_mat)

      s_mat = matmul (t_mat, s_mat)
      t_mat(1,2) = -t_mat(1,2); t_mat(2,1) = -t_mat(2,1) ! form inverse
      s_mat = matmul (s_mat, t_mat)
    endif

    pos = pos + matmul(w_mat, r_mat)
    w_mat = matmul (w_mat, s_mat)

  ! mirror

  case (mirror$)
    
    angle = 2 * ele%value(graze_angle$)
    tlt = ele%value(tilt_tot$)
    s_ang = sin(angle); c_ang = cos(angle)

    s_mat(1,:) = (/ c_ang,  0.0_dp, -s_ang /)
    s_mat(2,:) = (/ 0.0_dp, 1.0_dp,  0.0_dp /)
    s_mat(3,:) = (/ s_ang,  0.0_dp,  c_ang /) 

    if (tlt /= 0) then
      s_ang = sin(tlt); c_ang = cos(tlt)
      t_mat(1,:) = (/ c_ang,  -s_ang,  0.0_dp /)
      t_mat(2,:) = (/ s_ang,   c_ang,  0.0_dp /)
      t_mat(3,:) = (/ 0.0_dp,  0.0_dp, 1.0_dp /)

      s_mat = matmul (t_mat, s_mat)
      t_mat(1,2) = -t_mat(1,2); t_mat(2,1) = -t_mat(2,1) ! form inverse
      s_mat = matmul (s_mat, t_mat)
    endif

    w_mat = matmul (w_mat, s_mat)

  ! patch

  case (patch$)

    angle = ele%value(x_pitch_tot$)            ! x_pitch is negative MAD yrot
    if (angle /= 0) then
      s_ang = sin(angle); c_ang = cos(angle)
      s_mat(1,:) = (/  c_ang,  0.0_dp, s_ang /)
      s_mat(2,:) = (/  0.0_dp, 1.0_dp, 0.0_dp /)
      s_mat(3,:) = (/ -s_ang,  0.0_dp, c_ang /) 
      w_mat = matmul(w_mat, s_mat)
    endif
     
    angle = ele%value(y_pitch_tot$)           ! 
    if (angle /= 0) then
      s_ang = sin(angle); c_ang = cos(angle)
      s_mat(1,:) = (/ 1.0_dp,  0.0_dp, 0.0_dp /)
      s_mat(2,:) = (/ 0.0_dp,  c_ang,  s_ang /)
      s_mat(3,:) = (/ 0.0_dp, -s_ang,  c_ang /) 
      w_mat = matmul(w_mat, s_mat)
    endif

    angle = ele%value(tilt$)
    if (angle /= 0) then
      s_ang = sin(angle); c_ang = cos(angle)
      s_mat(1,:) = (/ c_ang,  -s_ang,  0.0_dp /)
      s_mat(2,:) = (/ s_ang,   c_ang,  0.0_dp /)
      s_mat(3,:) = (/ 0.0_dp,  0.0_dp, 1.0_dp /) 
      w_mat = matmul(w_mat, s_mat)
    endif

    r_mat = (/ ele%value(x_offset$), ele%value(y_offset$), &
                                                    ele%value(z_offset$) /)
    pos = pos + matmul(w_mat, r_mat)

  ! everything else. Just a translation

  case default
    pos = pos + w_mat(:,3) * leng

  end select

  ! if there has been a rotation calculate new theta, phi, and psi

  if (key == sbend$ .or. key == patch$ .or. key == multipole$) then
    if (abs(w_mat(1,3)) + abs(w_mat(3,3)) < 1e-12) then ! special degenerate case
      ! Note: Only theta +/- psi is well defined here so this is rather arbitrary.
      theta = floor0%theta  
      if (w_mat(2,3) > 0) then
        phi = pi/2
        psi = atan2(-w_mat(3,1), w_mat(1,1)) - theta
      else
        phi = -pi/2
        psi = atan2(w_mat(3,1), w_mat(1,1)) + theta
      endif
    else  ! normal case
      theta = atan2 (w_mat(1,3), w_mat(3,3))
      theta = theta - twopi_dp * &
                       nint((theta - floor0%theta) / twopi_dp)
      phi = atan2 (w_mat(2,3), sqrt(w_mat(1,3)**2 + w_mat(3,3)**2))
      psi = atan2 (w_mat(2,1), w_mat(2,2))
    endif
  endif

  old_theta = theta
  old_phi   = phi
  old_psi   = psi

! Simple case where the local reference frame stays in the horizontal plane.

else

  select case (key)
  case (sbend$)
    angle = leng * dble(ele%value(g$))
    chord_len = 2 * leng * sin(angle/2) / angle
  case (multipole$)
    angle = knl(0)
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

end subroutine
