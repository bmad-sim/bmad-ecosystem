subroutine propagate_ray (ray, s_end, lat, stop_at_extremum)

  use synrad_struct
  use synrad_interface, except => propagate_ray

  implicit none

  type (lat_struct), target :: lat
  type (ray_struct), target :: ray
  type (ele_struct), pointer :: ele

  real(rp) s_end, s_next, del_s, s_target
  real(rp) g, new_x, theta0, theta1, c_t0, c_t1

  logical stop_at_extremum

  ! find the target

  s_target = s_end

  if ((s_target - ray%now%s) * ray%direction < 0) s_target = &
          s_target + ray%direction * lat%param%total_length

  if (abs(s_target - ray%now%s) > 200) then
    print *, ' ERROR IN PROPAGATE_RAY: TRYING TO PROPAGATE TOO FAR.'
    print *, '      ', ray%now%s, s_end, s_target
    if (global_com%exit_on_error) call err_exit
  endif

  ! update old (but only if we have moved or gone through the IP)

  if (s_target == ray%now%s .and. &
            abs(s_end - ray%now%s) /= lat%param%total_length) return

  ray%old = ray%now

  ! Propagate the ray until we get to s_end

  propagation_loop: do

    ! First some bookkeeping...
    ! If we are crossing over to a new element then update ray%ix_ele.
    ! Additionally, if we cross the lat end we need to 
    ! reset ray%now%s and ray%ix_ele.
    ! Note: Since we can be going "backwards" to find a shadow then 
    ! ray%crossed_end can toggle from true to false.

    if (ray%direction == 1) then
      do
        if (ray%now%s .ge. lat%ele(ray%ix_ele)%s) then
          ray%ix_ele = ray%ix_ele + 1
          if (ray%ix_ele > lat%n_ele_track) then
            ray%ix_ele = 1
            ray%now%s = ray%now%s - lat%param%total_length
            ray%now%s = ray%now%s
            s_target = s_target - lat%param%total_length
            ray%crossed_end = .not. ray%crossed_end
          endif
        elseif (ray%now%s .lt. lat%ele(ray%ix_ele-1)%s) then
          ray%ix_ele = ray%ix_ele - 1
          if (ray%ix_ele == 0) then
            print *, 'ERROR IN PROPAGATE_RAY: INTERNAL + ERROR'
            if (global_com%exit_on_error) call err_exit
          endif
        else
          exit
        endif
      enddo

    else   ! direction = -1
      do
        if (ray%now%s .le. lat%ele(ray%ix_ele-1)%s) then
          ray%ix_ele = ray%ix_ele - 1
          if (ray%ix_ele .le. 0) then
            ray%ix_ele = lat%n_ele_track
            ray%now%s = ray%now%s + lat%param%total_length
            ray%now%s = ray%now%s
            s_target = s_target + lat%param%total_length
            ray%crossed_end = .not. ray%crossed_end
          endif
        elseif (ray%now%s .gt. lat%ele(ray%ix_ele)%s) then
          ray%ix_ele = ray%ix_ele + 1
          if (ray%ix_ele == lat%n_ele_track+1) then
            print *, 'ERROR IN PROPAGATE_RAY: INTERNAL - ERROR'
            if (global_com%exit_on_error) call err_exit
          endif
        else
          exit
        endif
      enddo
    endif

    if (ray%direction == 1) then
      s_next = min (s_target, lat%ele(ray%ix_ele)%s)
    else
      s_next = max (s_target, lat%ele(ray%ix_ele-1)%s)
    endif

    del_s = s_next - ray%now%s
    ele => lat%ele(ray%ix_ele)

    ! Now that the step length has been calculated, propagate the ray.

    if (ele%key == patch$) then



    ! In a bend: Exact formula is:
    !   new_x = (rho * (cos(theta0) - cos(theta1)) + ray%now%vec(1) * cos(theta0)) / cos(theta1)
    ! where theta = v_x / v_z
    ! We stop then theta = 0 since that is an extremum in x.

    elseif (ele%key == sbend$ .and. ele%value(g$) /= 0) then
      ! Rotate to ele coords if needed
      if (ele%value(ref_tilt_tot$) /= 0) call tilt_coords (ele%value(ref_tilt_tot$), ray%now%vec)
      g = ele%value(g$)
      theta0 = atan2(ray%now%vec(2), ray%now%vec(6))
      theta1 = theta0 + del_s * g
      ! if theta has changed sign then there is an extremum 
      if (stop_at_extremum .and. theta0 * theta1 < 0) then 
        del_s = -theta0 / g            ! step to extremum
        theta1 = 0
        s_next = ray%now%s + del_s
        s_target = s_next
      endif
      c_t0 = -(theta0**2)/2 + theta0**4/24
      c_t1 = -(theta1**2)/2 + theta1**4/24
      new_x = ((c_t0 - c_t1) / g + ray%now%vec(1) * cos(theta0)) / cos(theta1)
      ray%now%vec(1) = new_x
      ray%now%vec(2) = sin(theta1)
      ray%now%vec(3) = ray%now%vec(3) + del_s * ray%now%vec(4)
      ray%now%vec(6) = cos(theta1)
      if (ele%value(ref_tilt_tot$) /= 0) call tilt_coords (-ele%value(ref_tilt_tot$), ray%now%vec)

    ! Else not a patch or bend...

    else
      ray%now%vec(1) = ray%now%vec(1) + del_s * ray%now%vec(2)
      ray%now%vec(3) = ray%now%vec(3) + del_s * ray%now%vec(4)
    endif

    ray%track_len = ray%track_len + abs(del_s)
    ray%now%s = s_next
    ray%now%s      = s_next

    if (ray%crossed_end .and. lat%param%geometry == open$) return

    if (s_next == s_target) return

  enddo propagation_loop

end subroutine
