!+
! subroutine get_initial_pt (ray, wall, ix_wall, lat)
!
! Routine to get the initial point on the wall to track the ray to.
!
! Modules needed:
!   use synrad_mod
!
! Input:
!   ray    -- ray_struct:
!   wall   -- wall_struct: An outside or inside wall.
!   lat    -- lat_struct: Lattice.
!
! Output:
!   ix_wall -- Integer: wall%pt(:) index.
!-

subroutine get_initial_pt (ray, wall, ix_wall, lat)

  use synrad_struct
  use synrad_interface, except => get_initial_pt

  implicit none

  type (ray_struct) ray
  type (wall_struct) wall
  type (lat_struct) lat

  integer ix_wall, ix, ix0, ix1, ix2


  if (wall%n_pt_tot == 0) then
    print *, 'There are no points in the wall!'
    print *, 'You should check the wall first with check_wall!'
    call err_exit
  endif


! point ix_wall is at or just downstream of ray%now%vec(5).

! edge cases

  if (ray%now%vec(5) == lat%param%total_length) then
    ix_wall = wall%n_pt_tot
    wall%ix_pt = wall%n_pt_tot
    return
  endif

  if (ray%now%vec(5) == 0) then
    ix_wall = 0
    wall%ix_pt = 0
    return
  endif

! normal case. divide and conquer.

  ix0 = 0
  ix2 = wall%n_pt_tot

  do
    ix1 = (ix2 + ix0) / 2
    ix = wall%pt(ix1)%ix_pt
    if (wall%pt(ix)%s < ray%now%vec(5)) then
      ix0 = ix1
    elseif (wall%pt(ix)%s > ray%now%vec(5)) then
      ix2 = ix1
    elseif (ray%direction == 1) then   ! here wall%pt(ix)%s == ray%now%vec(5)
      ix2 = ix1
    else
      ix0 = ix1
    endif
    if (ix2 - ix0 == 1) then
      if (ray%direction == 1) then
        ix_wall = wall%pt(ix2)%ix_pt
        wall%ix_pt = ix2
      else
        ix_wall = wall%pt(ix0)%ix_pt
        wall%ix_pt = ix0
      endif
      return
    endif
  enddo

end subroutine
