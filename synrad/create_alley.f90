!+
! this routine marks the wall points according to whether they are part of an
! alley.
!
! To make the calculations tractible, there is a required restriction: 
! A "vertical" line (line drawn at constant s) must go through the wall
! exactly once or three times or be tangent to some wall segment. There
! is no restriction on x values except that a wall cannot intersect
! itself.
!-

subroutine create_alley (wall)

  use synrad_struct
  use synrad_interface, except => create_alley
  use nr
  use sim_utils

  implicit none

  type (wall_struct), target :: wall
  type (wall_pt_struct), pointer :: pt(:)

  integer i, j, i1, i2, i3, i4, i5, i6
  integer closed_end_direct

  real(rp) s_min_local, s_max_local, x_wall

! loop over all points until we come to a switchback which indicates
! we are in an alley.
! i1 is the beginning of the alley with i2 being the first point such that
!     s_wall(i1) >= s_min_local == s_wall(i4).
! i2 is the first point after i1 where s_wall is a local maximum == s_max_local.
!    [Note: There might be multiple points with s = s_max_local.]
! i3 is the point after i2 at which the wall starts going backward (s_wall
!    starts decreasing).
! i4 is the first point after i3 where s_wall is a local minimum.
!    [Note: There might be multiple points with s = s_min_local.]
! i5 is the point after i3 where the wall starts to go forward again.
! i6 is the end of the alley. It is the first point after i5 with:
!    s_wall(i6) >= s_max_local == s_wall(i2).

  pt => wall%pt

  i6 = 0
  do

    i3 = i6
    do
      i3 = i3 + 1
      if (pt(i3)%s < pt(i3-1)%s) exit

      ! At end cleanup: possible_alley -> no_alley$

      if (i3 == wall%n_pt_tot) then
        do i = 1, wall%n_pt_tot
          if (pt(i)%type == possible_alley$) pt(i)%type = no_alley$
        enddo
        return
      endif

    enddo

    s_max_local = pt(i3-1)%s

! go forward to 2nd switchback of the alley

    i5 = i3
    do
      i5 = i5 + 1
      if (pt(i5)%s > pt(i5-1)%s) exit
    enddo

    s_min_local = pt(i5-1)%s

! find the beginning of the second switchback

    i4 = i3
      do
        if (pt(i4)%s == s_min_local) exit
        i4 = i4 + 1
      enddo

! find the alley orientation

    i = i3
    do
      i = i - 1
      if (pt(i)%s < pt(i3)%s) then
        x_wall = pt(i)%x + (pt(i+1)%x - pt(i)%x) * &
              (pt(i3)%s - pt(i)%s) / (pt(i+1)%s - pt(i)%s)
        if (x_wall < pt(i3)%x) then
          closed_end_direct = -1
        else
          closed_end_direct = +1
        endif
        exit
      endif
    enddo

! find the beginning of the alley

    i1 = i3
    do
      if (pt(i1)%s == s_max_local) i2 = i1
      if (pt(i1-1)%s < s_min_local) exit
      i1 = i1 - 1
    enddo

! find the end of the alley

    i6 = i5
    do
      if (pt(i6)%s < pt(i6-1)%s) then
        print *, 'ERROR IN CREATE_ALLEY:'
        print *, '      I CANNOT HANDLE A DOUBLE BACKTRACK ALLEY'
        print *, '      ON THE WALL: ', wall_name(wall%side)
        print *, '      BETWEEN POINTS:', i6-1, i6
        do j = i1, i6-1
          print '(i6, 2x, a, f11.4)', j, wall%pt(j)%name, wall%pt(j)%s
        enddo
        call err_exit
      endif
      if (pt(i6)%s > s_max_local) exit
      i6 = i6 + 1
    enddo

! mark sections of the alley

    do i = i1, i6-1
      if (pt(i)%type /= possible_alley$) then
        print *, 'ERROR IN CREATE_ALLEY:'
        print *, '      FOUND A WALL ALLEY WHERE THERE SHOULD BE NONE.'
        print *, '      ON THE WALL: ', wall_name(wall%side)
        print *, '      BETWEEN POINTS:', i1, i6-1
        do j = max(1,i1-1), i6
          print '(i6, 2x, a, f11.4)', j, wall%pt(j)%name, wall%pt(j)%s
        enddo
        call err_exit
      endif
    enddo

    if (closed_end_direct == +1) then
      pt(i1:i2-1)%type = outer_wall$
      pt(i2:i3-1)%type = closed_end$
      pt(i3:i4-1)%type = middle_wall$
      pt(i4:i5-1)%type = open_end$
      pt(i5:i6-1)%type = inner_wall$
    else
      pt(i1:i2-1)%type = inner_wall$
      pt(i2:i3-1)%type = open_end$
      pt(i3:i4-1)%type = middle_wall$
      pt(i4:i5-1)%type = closed_end$
      pt(i5:i6-1)%type = outer_wall$
    endif

    pt(i1:i6-1)%closed_end_direct = closed_end_direct

! order the points in assending s

    call indexx(pt(i1:i6)%s, pt(i1:i6)%ix_pt)
    pt(i1:i6)%ix_pt = pt(i1:i6)%ix_pt + (i1 - 1)

  enddo

end subroutine
