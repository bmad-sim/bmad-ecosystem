!+
! Subroutine aperture_info (lattice, orbit, ix_lost, plane_lost)
!
! Subroutine to show where in an orbit a particle hit an aperture and was lost.
!
! The element in the lattice where the particle is lost is:
!   lattice%param%ix_lost
! If the particle is lost at the exit end of the element then
!   ix_lost = lattice%param%ix_lost
! If the particle is lost at the entrance end of the element then
!   ix_lost = lattice%param%ix_lost - 1
!
! plane_lost tells in what plane the aperture was hit. If the
! aperture is elliptical then what plane was the "major contributor" 
! for the particle getting lost is used. Note: In this case it is possible
! that:
!   orbit(ix_lost)%vec(1) < x_apperture and simultaneously
!   orbit(ix_lost)%vec(3) < y_apperture 
!
! Modules needed:
!   use bmad
!
! Input:
!   lattice   -- ring_struct: lattice tracked through.
!   orbit(0:) -- Coord_struct: coordinates.
!
! Output:
!   ix_lost    -- Integer: Index in orbit array of where particle is lost.
!   plane_lost -- Integer: Plane where particle is lost:
!                     x_plane$ or y_plane$
!-

subroutine aperture_info (lattice, orbit, ix_lost, plane_lost)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) lattice
  type (coord_struct) orbit(0:)

  real(rp) x_lim, y_lim
  integer ix_lost, plane_lost, ix
  character(16) :: r_name = 'aperture_info'

!

  if (.not. lattice%param%lost) then
    call out_io (s_fatal$, r_name, 'PARTICLE NOT LOST!')
    return
  endif

  ix = lattice%param%ix_lost

  select case (lattice%param%end_lost_at)
  case (entrance_end$)
    ix_lost = lattice%param%ix_lost - 1
  case (exit_end$)
    ix_lost = lattice%param%ix_lost 
  case default
    call out_io (s_abort$, r_name, 'INTERNAL ERROR')
  end select

  call check_aperture_limit (orbit(ix_lost), lattice%ele_(ix), &
                                                   lattice%param, plane_lost)

end subroutine
