!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! subroutine track_ray_to_wall (ray, lat, inside, outside, 
!                                hit_flag, track_max)
!
! subroutine to propagate a synch radiation ray until it hits
!    a wall
!
! Modules needed:
!   use sr_mod
!
! Input:
!   ray    -- ray_struct: synch radiation ray with starting
!                         parameters set
!   lat   -- lat_struct: with twiss propagated and mat6s made
!   inside  -- wall_struct: inside wall
!   outside -- wall_struct: outside wall
!   track_max -- real(rp), optional: Maximum length in m to track
!                                    the ray
!
! Output:
!   ray    -- ray_struct: synch radiation ray propagated to wall
!   hit_flag -- logical, optional: true if wall was hit,
!                                  false if track_max was reached first
!-

subroutine track_ray_to_wall (ray, lat, inside, outside, hit_flag, track_max)

  use sr_struct
  use sr_interface

  implicit none

  type (lat_struct), target :: lat
  type (ray_struct), target :: ray
  type (wall_struct) inside, outside

  logical, optional :: hit_flag
  real(rp), optional :: track_max

  integer ix_in, ix_out

  real(rp) s_next

  logical is_hit

! init

  if (present(hit_flag)) hit_flag = .true.  ! assume that we will hit

! ix_in and ix_out are the next inside and outside wall points that
! are at or just "downstream" of the ray.

  call get_initial_pt (ray, inside, ix_in, lat)
  call get_initial_pt (ray, outside, ix_out, lat)

! propagation loop:
! Propagate the ray. Figure out how far to advance in s.
! Do not advance past the next wall point (either inside or outside).

  do

    if (ray%direction == 1) then
      s_next = min(inside%pt(ix_in)%s, outside%pt(ix_out)%s, &
                                                 ray%now%vec(5) + 1.0)
    else
      s_next = max(inside%pt(ix_in)%s, outside%pt(ix_out)%s, &
                                                 ray%now%vec(5) - 1.0)
    endif

    call propagate_ray (ray, s_next, lat)

! See if we are outside the beam pipe.
! If so we calculate the exact hit spot where the ray crossed the
! wall boundry and return

    call hit_spot_calc (ray, inside, ix_in, is_hit, lat)
    if (is_hit) return

    call hit_spot_calc (ray, outside, ix_out, is_hit, lat)
    if (is_hit) return

    if (present(track_max)) then
      if (ray%track_len .ge. track_max) then
        hit_flag = .false.
        return
      endif
    endif

    if (ray%now%vec(5) == inside%pt(ix_in)%s) then
      call next_pt (ray, inside, ix_in)
    endif

    if (ray%now%vec(5) == outside%pt(ix_out)%s) then
      call next_pt (ray, outside, ix_out)
    endif

  enddo

end subroutine
