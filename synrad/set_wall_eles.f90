!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine set_wall_eles (walls, lat)

  use sr_struct
  use sr_interface

  implicit none

  type (walls_struct), target :: walls
  type (wall_struct), pointer :: wall
  type (lat_struct) :: lat
  type (ele_struct) :: ele

  integer i, is, ix_ele

  !
  do i = 1,2
    if (i == 1) wall => walls%positive_x_wall
    if (i == 2) wall => walls%negative_x_wall
    do is = 1, wall%n_seg_tot

      ! find the element at the s midpoint of the wall segment
      call ele_at_s (lat, wall%seg(is)%s_mid, ix_ele)
      wall%seg(is)%ix_ele = ix_ele

      ! Get the betas at the s midpoint
      call twiss_and_track_at_s (lat, wall%seg(is)%s_mid, ele)
      wall%seg(is)%a = ele%a
      wall%seg(is)%b = ele%b
      
    end do

  end do


end subroutine set_wall_eles
