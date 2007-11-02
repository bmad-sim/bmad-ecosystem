!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine init_wall_ends (walls)

  use sr_struct
  use sr_interface

  implicit none

  type (walls_struct) walls
	real(rp) end_n_x, end_p_x
	real(rp) end_n_s, end_p_s
	
! assuming both walls start at s=0
  walls%initial_end%ix_pt = -1
  walls%initial_end%s = 0
  walls%initial_end%x = walls%positive_x_wall%pt(0)%x
  walls%initial_end%s_mid = 0
  walls%initial_end%len = abs(walls%positive_x_wall%pt(0)%x - &
  	walls%negative_x_wall%pt(0)%x)
  walls%initial_end%x_mid = walls%negative_x_wall%pt(0)%x + &
  	walls%initial_end%len/2.
  walls%initial_end%sr_power%power = 0
  walls%initial_end%sr_power%power_per_len = 0
  walls%initial_end%sr_power%power_per_area = 0
  walls%initial_end%sr_power%ix_ele_source = 0
  walls%initial_end%sr_power%s_source = 0
  walls%initial_end%sr_power%n_source = 0
  nullify (walls%initial_end%sr_power%sources)
  
  end_n_s = walls%negative_x_wall%seg(walls%negative_x_wall%n_seg_tot)%s
  end_p_s = walls%positive_x_wall%seg(walls%positive_x_wall%n_seg_tot)%s
  end_n_x = walls%negative_x_wall%seg(walls%negative_x_wall%n_seg_tot)%x
  end_p_x = walls%positive_x_wall%seg(walls%positive_x_wall%n_seg_tot)%x

  walls%final_end%ix_pt = -1
  walls%final_end%s = end_n_s
  walls%final_end%x = end_n_x
  walls%final_end%s_mid = (end_n_s + end_p_s)/2.
  walls%final_end%len = sqrt( (end_p_s - end_n_s)**2 + (end_p_x - end_n_x)**2) 

  walls%final_end%x_mid = (end_n_x + end_p_x)/2.
  walls%final_end%sr_power%power = 0
  walls%final_end%sr_power%power_per_len = 0
  walls%final_end%sr_power%power_per_area = 0
  walls%final_end%sr_power%ix_ele_source = 0
  walls%final_end%sr_power%s_source = 0
  walls%final_end%sr_power%n_source = 0
  nullify (walls%final_end%sr_power%sources)


end subroutine init_wall_ends
