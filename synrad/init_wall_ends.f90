subroutine init_wall_ends (walls)

  use synrad_struct
  use synrad_interface, except => init_wall_ends

  implicit none

  type (walls_struct) walls
	real(rp) end_n_x, end_p_x
	real(rp) end_n_s, end_p_s
	
! assuming both walls start at s=0
  walls%start_end%ix_pt = -1
  walls%start_end%s = 0
  walls%start_end%x = walls%positive_x_wall%pt(0)%x
  walls%start_end%s_mid = 0
  walls%start_end%len = abs(walls%positive_x_wall%pt(0)%x - &
  	walls%negative_x_wall%pt(0)%x)
  walls%start_end%x_mid = walls%negative_x_wall%pt(0)%x + &
  	walls%start_end%len/2.
  walls%start_end%power%power_tot = 0
  walls%start_end%power%power_per_len = 0
  walls%start_end%power%power_per_area = 0
  walls%start_end%power%ix_ele_source = 0
  walls%start_end%power%s_source = 0
  walls%start_end%power%n_source = 0
  nullify (walls%start_end%power%sources)
  
  end_n_s = walls%negative_x_wall%seg(walls%negative_x_wall%n_seg_tot)%s
  end_p_s = walls%positive_x_wall%seg(walls%positive_x_wall%n_seg_tot)%s
  end_n_x = walls%negative_x_wall%seg(walls%negative_x_wall%n_seg_tot)%x
  end_p_x = walls%positive_x_wall%seg(walls%positive_x_wall%n_seg_tot)%x

  walls%exit_end%ix_pt = -1
  walls%exit_end%s = end_n_s
  walls%exit_end%x = end_n_x
  walls%exit_end%s_mid = (end_n_s + end_p_s)/2.
  walls%exit_end%len = sqrt( (end_p_s - end_n_s)**2 + (end_p_x - end_n_x)**2) 

  walls%exit_end%x_mid = (end_n_x + end_p_x)/2.
  walls%exit_end%power%power_tot = 0
  walls%exit_end%power%power_per_len = 0
  walls%exit_end%power%power_per_area = 0
  walls%exit_end%power%ix_ele_source = 0
  walls%exit_end%power%s_source = 0
  walls%exit_end%power%n_source = 0
  nullify (walls%exit_end%power%sources)


end subroutine init_wall_ends
