!+
! Subroutine multipole_kick (knl, tilt, n, coord)
!
! Subroutine to put in the kick due to a multipole.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!                          
! Input:
!   knl   -- Real: Multipole strength (mad units).
!   tilt  -- Real: Multipole tilt.
!   n     -- Real: Multipole order.
!   coord -- Coord_struct:
!     %x%pos -- X position.
!     %y%pos -- Y position.
!
! Output:
!   coord -- Coord_struct: 
!     %x%vel -- X kick.
!     %y%vel -- Y kick.
!-

subroutine multipole_kick (knl, tilt, n, coord)

  use bmad_struct
  implicit none

  type (coord_struct)  coord

  real knl, tilt, c(0:n_pole_maxx, 0:n_pole_maxx), x, y, sin_ang, cos_ang
  real mexp, x_vel, y_vel

  integer n, m

  logical init_needed / .true. /

! init

  if (init_needed) then
    call multipole_c_init (c, n_pole_maxx)
    init_needed = .false.
  endif

! simple case

  if (knl == 0) return

! normal case

  if (tilt == 0) then
    x = coord%x%pos
    y = coord%y%pos
  else
    sin_ang = sin(tilt)
    cos_ang = cos(tilt)
    x =  coord%x%pos * cos_ang + coord%y%pos * sin_ang
    y = -coord%x%pos * sin_ang + coord%y%pos * cos_ang
  endif

  x_vel = 0
  y_vel = 0

  do m = 0, n, 2
    x_vel = x_vel + knl * c(n, m) * mexp(x, n-m) * mexp(y, m)
  enddo

  do m = 1, n, 2
    y_vel = y_vel + knl * c(n, m) * mexp(x, n-m) * mexp(y, m)
  enddo

  if (tilt == 0) then
    coord%x%vel = coord%x%vel + x_vel
    coord%y%vel = coord%y%vel + y_vel
  else
    coord%x%vel = coord%x%vel + x_vel * cos_ang - y_vel * sin_ang
    coord%y%vel = coord%y%vel + x_vel * sin_ang + y_vel * cos_ang
  endif

end subroutine
