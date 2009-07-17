!+
! Gets next point in increasing or decreasing s depending upon ray%direction.
! -

subroutine next_pt (ray, wall, ix_wall, passed_end)

  use synrad_struct
  use synrad_interface, except => next_pt

  implicit none

  type (ray_struct) ray
  type (wall_struct) wall

  integer ix_wall
  logical passed_end
  
  integer direct

! check

  if (wall%pt(wall%ix_pt)%ix_pt /= ix_wall) then
    print *, 'ERROR IN NEXT_PT: INTERNAL ERROR'
    call err_exit
  endif

! wrap around cases

  direct = ray%direction

	passed_end = .false.

  if (ix_wall == 0 .and. direct == -1) then
    ix_wall = wall%n_pt_tot
    wall%ix_pt = wall%n_pt_tot
    passed_end = .true.
    return
  endif

  if (ix_wall == wall%n_pt_tot .and. direct == 1) then
    ix_wall = 0
    wall%ix_pt = 0
	passed_end = .true.
    return
  endif

! normal

  wall%ix_pt = wall%ix_pt + direct
  ix_wall = wall%pt(wall%ix_pt)%ix_pt


end subroutine
