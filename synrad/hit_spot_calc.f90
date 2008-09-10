subroutine hit_spot_calc (ray, wall, ix_wall, has_hit, lat, circular)

  use sr_struct
  use sr_interface

  implicit none

  type (lat_struct) lat
  type (ray_struct) :: ray, ray0, ray1, ray2
  type (wall_struct), target :: wall
  type (wall_pt_struct), pointer :: pt(:)
  type (wall_pt_struct), pointer :: pt0

  integer ix_wall, ix0, ix1, ix2, i

  real(rp) dx_wall, ds_wall, del_s, s1, denom
  real(rp) del0, del1, del2, x, r_wall

  logical has_hit, circular

! init

  pt => wall%pt
  pt0 => pt(ix_wall)

  has_hit = .false.  ! assume no hit

! if the ray is on the opposite side of the beam centerline from the wall there
! is no hit

  if (sign(1.0, ray%now%vec(1)) /= sign(1.0, pt0%x)) return

! figure out if there is a hit.
! if there is no alley here then we have not hit if x_pos < x_wall

  if (pt0%type == no_alley$) then
    ray%alley_status = no_local_alley$
    ix1 = ix_wall - ray%direction
    if (ray%now%vec(5) == pt0%s) then
      x = pt0%x
    else
      x = pt0%x + (pt(ix1)%x - pt0%x) * &       ! interpolate to get x
                    (ray%now%vec(5) - pt0%s) / (pt(ix1)%s - pt0%s)
    endif
    if (abs(ray%now%vec(1)) < abs(x)) return  ! return if no hit.

! Here for an alley...
! Note. To simplify the program we do not allow a calculation of whether
! we have hit if we are in an ally but between wall points (that is, in
! an alley the wall points cannot be further apart than 1 meter).

  elseif (ray%now%vec(5) == pt0%s) then   ! in alley at a wall pt

    if (pt0%type == open_end$) then
      if (ray%direction == pt0%closed_end_direct) then
        if (abs(ray%now%vec(1)) .gt. abs(pt0%x)) then
          ray%alley_status = in_alley$
        else
          ray%alley_status = out_of_alley$
        endif
        ix1 = ix_wall - 1
        if (pt(ix1)%type == open_end$) then
          if (abs(ray%now%vec(1)) > max(abs(pt(ix1)%x), abs(pt0%x))) return
          if (abs(ray%now%vec(1)) < min(abs(pt(ix1)%x), abs(pt0%x))) return
        else
          return
        endif
      else
        if (ray%alley_status == in_alley$) then
          ix1 = ix_wall + ray%direction
          if (pt(ix1)%type == open_end$) return
          if (abs(ray%now%vec(1)) > abs(pt0%x)) return
        else
          ix1 = ix_wall - ray%direction
          if (pt(ix1)%type == open_end$) return
          if (abs(ray%now%vec(1)) < abs(pt0%x)) return
        endif
      endif

    elseif (pt0%type == closed_end$) then
      if (ray%direction /= pt0%closed_end_direct) return
      if (ray%alley_status /= in_alley$) return
      if (abs(ray%now%vec(1)) < abs(pt0%x)) then
        ix1 = ix_wall + ray%direction
        if (pt(ix1)%type == closed_end$ .and. &
             abs(ray%now%vec(1)) < abs(pt(ix1)%x)) return
      else
        ix1 = ix_wall - ray%direction
        if (pt(ix1)%type == closed_end$ .and. &
             abs(ray%now%vec(1)) > abs(pt(ix1)%x)) return
      endif

    elseif (pt0%type == inner_wall$) then
      if (ray%alley_status == in_alley$) return
      ix1 = ix_wall - ray%direction
      if (abs(ray%now%vec(1)) < abs(pt0%x)) return

    elseif (pt0%type == middle_wall$) then
      if (ray%alley_status /= in_alley$) return
      ix1 = ix_wall + ray%direction
      if (abs(ray%now%vec(1)) > abs(pt0%x)) return

    elseif (pt0%type == outer_wall$) then
      if (ray%alley_status /= in_alley$) return
      ix1 = ix_wall - ray%direction
      if (abs(ray%now%vec(1)) < abs(pt0%x)) return

    else
      type *, 'ERROR IN HIT_SPOT_CALC: ALLEY ERROR'
      call err_exit
    endif

  else
    if (pt0%type /= no_alley$) then
      type *, 'ERROR IN HIT_SPOT_CALC: CALCULATION IN ALLEY BUT BETWEEN POINTS.'
      print *, "ray%now%vec(5): <", ray%now%vec(5), &
           ">  pt0%s: <", pt0%s, ">"
      print *, "ray%now%vec(1): <", ray%now%vec(1), &
           ">  pt0%x: <", pt0%x, ">"
      call err_exit
    endif
  endif

  ix0 = min(ix_wall, ix1)
  ix2 = max(ix_wall, ix1)
  has_hit = .true.
  ray%ix_wall_pt = ix2

! Here if we have a hit.
! We need to find where exactly the ray hit the wall.
! ray%r_wall is the percentage distance along the wall piece from pt(ix0)
! where the ray hits.
! ray%r_wall = 0.0 => the hit is at pt(ix0)
! ray%r_wall = 1.0 => the hit is at pt(ix2)
! we need to iterate in a bend since the wall is actually curved.

  ray0 = ray
  ray1 = ray
  ray2 = ray

  if (ray%now%vec(5) < ray%old%vec(5)) then
    ray2%now = ray%old
  elseif (ray%now%vec(5) > ray%old%vec(5)) then
    ray0%now = ray%old
  endif

  dx_wall = pt(ix2)%x - pt(ix0)%x
  ds_wall = pt(ix2)%s - pt(ix0)%s
  denom = sqrt (dx_wall**2 + ds_wall**2)

  del0 = (dx_wall*(ray0%now%vec(5) - pt(ix0)%s) - &
       ds_wall*(ray0%now%vec(1) - pt(ix0)%x)) / denom
  del2 = (dx_wall*(ray2%now%vec(5) - pt(ix0)%s) - &
       ds_wall*(ray2%now%vec(1) - pt(ix0)%x)) / denom

  do i = 1, 20

    if (abs(del0) < 1.0e-4) then
      ray1 = ray0
      exit
    elseif (abs(del2) < 1.0e-4) then
      ray1 = ray2
      exit
    endif

    if (i == 20) then
      print *, 'ERROR IN HIT_SPOT_CALC: CALCULATION IS NOT CONVERGING'
!      print *, del0,del1
      call err_exit
    endif

    s1 = (del2 * ray0%now%vec(5) - del0 * ray2%now%vec(5)) / (del2 - del0)
    if (s1 < min(pt(ix0)%s, pt(ix2)%s) .or. s1 > max(pt(ix0)%s, pt(ix2)%s)) then
      type *, 'ERROR IN HIT_SPOT_CALC: INTERPOLATION ERROR'
      call err_exit
    endif

    if (s1 < ray1%now%vec(5)) then
      ray1%direction = -1
    else
      ray1%direction = +1
    endif
    call propagate_ray (ray1, s1, lat, .false.,circular)

    del1 = (dx_wall*(ray1%now%vec(5) - pt(ix0)%s) - &
                        ds_wall*(ray1%now%vec(1) - pt(ix0)%x)) / denom

    if (s1 < ray0%now%vec(5)) then
      ray2 = ray0; del2 = del0
      ray0 = ray1; del0 = del1
    elseif (s1 > ray2%now%vec(5)) then
      ray0 = ray2; del0 = del2
      ray2 = ray1; del2 = del1
    elseif (sign(1.0, del0) == sign(1.0, del1)) then
      ray0 = ray1; del0 = del1
    else
      ray2 = ray1; del2 = del1
    endif

  enddo

! cleanup

  r_wall = (dx_wall*(ray1%now%vec(1) - pt(ix0)%x) + &
                ds_wall*(ray1%now%vec(5) - pt(ix0)%s)) / denom**2
  if (r_wall > 1 .and. r_wall <  1.0001) r_wall = 1
  if (r_wall < 0 .and. r_wall > -0.0001) r_wall = 0

  if (r_wall > 1 .or. r_wall < 0) then
    type *, 'ERROR IN HIT_SPOT_CALC: R_WALL OUT OF BOUNDS.', r_wall
    call err_exit
  endif

  ray%now = ray1%now
  ray%r_wall = r_wall

  ! We assume that the travel length cannot be greater 
  ! then half the circumference.

  del_s = ray%now%vec(5) - ray%start%vec(5)
  if (del_s*ray%direction < 0) then
    ray%track_len = lat%param%total_length - abs(del_s)
    ray%crossed_end = .true.
!    print *, "ray from ", lat%ele(ray%ix_source)%name, "crossed end."
  else
    ray%track_len = abs(del_s)
    ray%crossed_end = .false.
  endif

  ray%wall => wall

end subroutine
