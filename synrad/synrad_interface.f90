module synrad_interface

use synrad_struct
use bmad_interface

interface
  function theta_floor (s, lat, theta_base) result (theta_fl)
    use synrad_struct
    implicit none
    type (lat_struct) lat
    real(rp) s, theta_fl
    real(rp), optional :: theta_base
  end function
end interface

interface
  subroutine break_wall_into_segments (wall, seg_len_max)
    use synrad_struct
    implicit none
    type (wall_struct) wall
    real(rp) seg_len_max
  end subroutine
end interface

interface
  subroutine calculate_synrad_power (lat, orb, direction, power, &
       walls, gen, use_ele_ix)
    use synrad_struct
    implicit none
    type (lat_struct), target :: lat
    type (coord_struct) orb(0:)
    type (walls_struct), target :: walls
    type (synrad_param_struct) gen
    type (ele_power_struct) power(:)
    integer direction
    integer use_ele_ix
  end subroutine
end interface

interface
  subroutine ele_synrad_power (lat, ie, orb, direction, power, &
       walls, gen)
    use synrad_struct
    implicit none
    type (lat_struct), target :: lat
    type (coord_struct) orb(0:)
    type (walls_struct), target :: walls
    type (synrad_param_struct) gen
    type (ele_power_struct) power(:)
    integer direction, ie
  end subroutine
end interface

interface
  subroutine hit_spot_calc (ray, wall, ix_wall, has_hit, lat)
    use synrad_struct
    implicit none
    type (lat_struct) lat
    type (ray_struct) :: ray
    type (wall_struct), target :: wall
    integer ix_wall
    logical has_hit
  end subroutine
end interface

interface
  subroutine check_end (point, ix, string_in)
    use synrad_struct
    implicit none
    type (outline_pt_struct) point(:)
    integer ix
    character(*) string_in
  end subroutine
end interface

interface
  subroutine check_wall (wall, s_lat)
    use synrad_struct, only: wall_struct, rp
    implicit none
    type (wall_struct) wall
    real(rp) s_lat
  end subroutine
end interface

interface
  subroutine track_ray_to_wall (ray, lat, walls, &
       hit_flag, track_max)
    use synrad_struct
    implicit none
    type (lat_struct), target :: lat
    type (ray_struct), target :: ray
    type (walls_struct), target :: walls
    logical, optional :: hit_flag
    real(rp), optional :: track_max
  end subroutine
end interface

interface
  subroutine create_alley (wall)
    use synrad_struct
    implicit none
    type (wall_struct), target :: wall
  end subroutine
end interface

interface
  subroutine delete_overlapping_wall_points (wall)
    use synrad_struct
    implicit none
    type (wall_struct), target :: wall
  end subroutine
end interface

interface
  subroutine init_ray (ray, lat, ix_ele, l_offset, orb, direction)
    use synrad_struct
    implicit none
    type (lat_struct), target :: lat
    type (coord_struct) orb(0:*)
    type (ray_struct) ray
    real(rp) l_offset
    integer direction
    integer ix_ele
  end subroutine
end interface

interface
  subroutine outline_concat (outline1, outline2, outline3)
    use synrad_struct
    implicit none
    type (outline_struct) outline1
    type (outline_struct) outline2
    type (outline_struct) outline3
  end subroutine
end interface

interface
  subroutine outline_reverse (outline1, outline2)
    use synrad_struct
    implicit none
    type (outline_struct) outline1
    type (outline_struct) outline2
  end subroutine
end interface

interface
  subroutine propagate_ray (ray, s_end, lat, stop_at_extremum)
    use synrad_struct
    implicit none
    type (lat_struct), target :: lat
    type (ray_struct), target :: ray
    real(rp) s_end
    logical stop_at_extremum
  end subroutine
end interface

interface
  subroutine synrad_read_vac_wall_geometry (wall_file, component_file, dflt_dir, s_lat, walls)
    use synrad_struct, only: walls_struct, rp
    implicit none
    character(*) wall_file, component_file, dflt_dir
    type (walls_struct) walls
    real(rp) s_lat
    logical type_warning
  end subroutine
end interface

interface
  subroutine seg_power_calc (rays, i_ray, walls, lat, gen, power)
    use synrad_struct
    implicit none
    type (ray_struct) :: rays(:)
    type (walls_struct), target :: walls
    type (synrad_param_struct) gen
    type (lat_struct)          lat
    type (ele_power_struct)     power
    integer i_ray
  end subroutine
end interface

interface
  subroutine get_initial_pt (ray, wall, ix_wall, lat)
    use synrad_struct
    implicit none
    type (ray_struct) ray
    type (wall_struct) wall
    type (lat_struct) lat
    integer ix_wall
  end subroutine
end interface

interface
  subroutine init_wall (wall)
    use synrad_struct
    implicit none
    type (wall_struct) wall
  end subroutine
end interface

interface
  subroutine init_wall_ends (walls)
    use synrad_struct
    implicit none
    type (walls_struct) walls
  end subroutine
end interface

interface
  subroutine next_pt (ray, wall, ix_wall, passed_end)
    use synrad_struct
    implicit none
    type (ray_struct) ray
    type (wall_struct) wall
    integer ix_wall
	  logical passed_end
  end subroutine
end interface

! This is to suppress the ranlib "has no symbols" message
integer, private :: private_dummy

end module
