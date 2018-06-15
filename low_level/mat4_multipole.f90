!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mat4_multipole (knl, tilt, n, orbit, kick_mat)
!
! Subroutine to find the kick matrix (Jacobian) due to a multipole.
! This routine is not meant for general use.
!
! Input:
!   orbit   -- Coord_struct: coordinates of particle
!   knl     -- Real(rp): Strength of multipole
!   tilt    -- Real(rp): Tilt of multipole
!
! Output:
!   kick_mat(4,4) -- Real(rp): Kick matrix (Jacobian) at orbit.
!-


subroutine mat4_multipole (knl, tilt, n, orbit, kick_mat)

use equal_mod, dummy => mat4_multipole

implicit none

type (coord_struct) orbit

real(rp) x_pos, y_pos, x, y, knl, tilt
real(rp) sin_ang, cos_ang, mat(2,2), rot(2,2)
real(rp) kick_mat(4,4)

integer m, n

! init

kick_mat = 0
forall (m = 1:4) kick_mat(m,m) = 1

x_pos = orbit%vec(1)
y_pos = orbit%vec(3)
         
! simple case

if (knl == 0 .or. (x_pos == 0 .and. y_pos == 0 .and. n /= 1)) then
  kick_mat(2:4:2, 1:3:2) = 0
  return
endif

! get position of particle in frame of multipole

if (tilt == 0) then
  x = x_pos
  y = y_pos
else
  sin_ang = sin(tilt)
  cos_ang = cos(tilt)
  x =  x_pos * cos_ang + y_pos * sin_ang
  y = -x_pos * sin_ang + y_pos * cos_ang
endif

! compute kick matrix

mat = 0

do m = 0, n, 2
  mat(1,1) = mat(1,1) +  knl * (n-m) * c_multi(n, m) * mexp(x, n-m-1) * mexp(y, m)
  mat(1,2) = mat(1,2) +  knl * m * c_multi(n, m) * mexp(x, n-m) * mexp (y, m-1)
enddo

do m = 1, n, 2
  mat(2,1) = mat(2,1) +  knl * (n-m) * c_multi(n, m) * mexp(x, n-m-1) * mexp(y, m)
  mat(2,2) = mat(2,2) +  knl * m * c_multi(n, m) * mexp(x, n-m) * mexp(y, m-1)
enddo

! transform back to lab frame

if (tilt /= 0) then
  rot(1,1) =  cos_ang
  rot(1,2) = -sin_ang
  rot(2,1) =  sin_ang
  rot(2,2) =  cos_ang
  mat = matmul(rot, mat)
  rot(1,2) =  sin_ang
  rot(2,1) = -sin_ang
  mat = matmul (mat, rot)
endif

kick_mat(2,1) = mat(1,1)
kick_mat(2,3) = mat(1,2)
kick_mat(4,1) = mat(2,1)
kick_mat(4,3) = mat(2,2)

end subroutine mat4_multipole

