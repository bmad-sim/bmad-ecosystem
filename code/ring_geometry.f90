!+
! Subroutine ring_geometry (ring)
!
! Subroutine to calculate the physical placement of all the elements in a ring.
! That is, the layout on the floor. 
! This is the same as the MAD convention. See the MAD manual for more details.
!
! Note: At present this routine assumes no  vertical bends. That is, 
! y_position is always 0 and the ring in in the X-Z plane.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring -- ring_struct: The ring
!
! Output:
!   ring
!     %ele_(i)
!       %x_position        -- X position at end of element
!       %y_position        -- Y position at end of element
!       %z_position        -- Z position at end of element
!       %theta_position    -- Orientation angle at end of element in X-Z plane
!
! Note: The starting point is taken to be the position in RING%ELE_(0)
!-

#include "CESR_platform.inc"

subroutine ring_geometry (ring)

  use bmad_struct
  use bmad_interface
  use multipole_mod

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele

  integer i, key

  real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
  real(8) chord_len, angle, leng, old_theta, rho
  real(8) s_ang, c_ang, s_the, c_the, s_phi, c_phi, s_psi, c_psi
  real(8) pos(3), theta, phi, psi, tlt
  real(8) w_mat(3,3), s_mat(3,3), r_mat(3), t_mat(3,3)
  real(8) :: twopi_8 = 2 * 3.14159265358979

! init
! old_theta is used to tell if we have to reconstruct the w_mat

  pos(1) = ring%ele_(0)%position%x
  pos(2) = ring%ele_(0)%position%y
  pos(3) = ring%ele_(0)%position%z
  theta = ring%ele_(0)%position%theta
  phi   = ring%ele_(0)%position%phi
  psi   = ring%ele_(0)%position%psi

  old_theta = 100  ! garbage number

!

  do i = 1, ring%n_ele_ring

    ele => ring%ele_(i)
    leng = ele%value(l$)
    key = ele%key
    if (key == sbend$ .and. (leng == 0 .or. ele%value(g$) == 0)) key = drift$

    if (key == ab_multipole$ .or. key == multipole$) then
      call multipole_ele_to_kt (ele, ring%param%particle, knl, tilt, .true.)
      key = multipole$
    endif

! General case where layout is not in the horizontal plane

    if (phi /= 0 .or. psi /= 0 .or. key == patch$ .or. &
               (key == multipole$ .and. knl(0) /= 0 .and. tilt(0) /= 0) .or. &
               (key == sbend$ .and. ele%value(tilt_tot$) /= 0)) then

      if (old_theta /= theta) then
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
          rho = 1.0_8 / ele%value(g$)
          s_ang = sin(angle); c_ang = cos(angle)
          r_mat = (/ rho * (c_ang - 1), 0.0_8, rho * s_ang /)
        else
          angle = knl(0)
          tlt = tilt(0)
          s_ang = sin(angle); c_ang = cos(angle)
          r_mat = 0
        endif
        s_mat(1,:) = (/ c_ang, 0.0_8, -s_ang /)
        s_mat(2,:) = (/ 0.0_8, 1.0_8,  0.0_8 /)
        s_mat(3,:) = (/ s_ang, 0.0_8,  c_ang /) 
        pos = pos + matmul(w_mat, r_mat)

        if (tlt /= 0) then
          s_ang = sin(tlt); c_ang = cos(tlt)
          t_mat(1,:) = (/ c_ang, -s_ang, 0.0_8 /)
          t_mat(2,:) = (/ s_ang,  c_ang, 0.0_8 /)
          t_mat(3,:) = (/ 0.0_8,  0.0_8, 1.0_8 /)
          s_mat = matmul (t_mat, s_mat)
          t_mat(1,2) = -t_mat(1,2); t_mat(2,1) = -t_mat(2,1) ! form inverse
          s_mat = matmul (s_mat, t_mat)
        endif

        w_mat = matmul (w_mat, s_mat)

! patch

      case (patch$)

        angle = ele%value(x_pitch$)            ! x_pitch is negative MAD yrot
        if (angle /= 0) then
          s_ang = sin(angle); c_ang = cos(angle)
          s_mat(1,:) = (/  c_ang, 0.0_8, s_ang /)
          s_mat(2,:) = (/  0.0_8, 1.0_8, 0.0_8 /)
          s_mat(3,:) = (/ -s_ang, 0.0_8, c_ang /) 
          w_mat = matmul(w_mat, s_mat)
        endif
     
        angle = ele%value(y_pitch$)            ! 
        if (angle /= 0) then
          s_ang = sin(angle); c_ang = cos(angle)
          s_mat(1,:) = (/ 1.0_8,  0.0_8, 0.0_8 /)
          s_mat(2,:) = (/ 0.0_8,  c_ang, s_ang /)
          s_mat(3,:) = (/ 0.0_8, -s_ang, c_ang /) 
          w_mat = matmul(w_mat, s_mat)
        endif

        angle = ele%value(tilt$)
        if (angle /= 0) then
          s_ang = sin(angle); c_ang = cos(angle)
          s_mat(1,:) = (/ c_ang, -s_ang, 0.0_8 /)
          s_mat(2,:) = (/ s_ang,  c_ang, 0.0_8 /)
          s_mat(3,:) = (/ 0.0_8,  0.0_8, 1.0_8 /) 
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
        theta = atan2 (w_mat(1,3), w_mat(3,3))
        theta = theta - twopi_8 * &
                       nint((theta - ring%ele_(i-1)%position%theta) / twopi_8)
        phi = atan2 (w_mat(2,3), sqrt(w_mat(1,3)**2 + w_mat(3,3)**2))
        psi = atan2 (w_mat(2,1), w_mat(2,2))
      endif

      old_theta = theta

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

    ring%ele_(i)%position%x = pos(1)
    ring%ele_(i)%position%y = pos(2)
    ring%ele_(i)%position%z = pos(3)
    ring%ele_(i)%position%theta = theta
    ring%ele_(i)%position%phi   = phi
    ring%ele_(i)%position%psi   = psi

  enddo

end subroutine
