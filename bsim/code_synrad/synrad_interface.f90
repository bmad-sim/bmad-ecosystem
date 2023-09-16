module synrad_interface

use synrad_struct

interface
  subroutine set_wall_eles (wall, branch)
    import
    implicit none
    type (wall_struct) wall
    type (branch_struct) branch
  end subroutine
end interface

interface
  function ray_is_outside_wall (ray, walls) result (is_outside)
    import
    implicit none
    type (ray_struct) ray
    type (walls_struct), target :: walls
    logical is_outside
  end function
end interface

interface
  function theta_floor (s, branch, theta_base) result (theta_fl)
    import
    implicit none
    type (branch_struct) branch
    real(rp) s, theta_fl
    real(rp), optional :: theta_base
  end function
end interface

interface
  subroutine break_wall_into_segments (wall, seg_len_max, branch, seg_len_phantom_max)
    import
    implicit none
    type (wall_struct), target :: wall
    type (branch_struct) branch
    real(rp) seg_len_max
    real(rp), optional :: seg_len_phantom_max
  end subroutine
end interface

interface
  subroutine calculate_synrad_power (branch, orb, direction, power, walls, gen, ix_ele1, ix_ele2)
    import
    implicit none
    type (branch_struct), target :: branch
    type (coord_struct) orb(0:)
    type (walls_struct), target :: walls
    type (synrad_param_struct) gen
    type (ele_power_struct) power(:)
    integer direction
    integer ix_ele1, ix_ele2
  end subroutine
end interface

interface
  subroutine ele_synrad_power (branch, ie, orb, direction, power, walls, gen)
    import
    implicit none
    type (branch_struct), target :: branch
    type (coord_struct) orb(0:)
    type (walls_struct), target :: walls
    type (synrad_param_struct) gen
    type (ele_power_struct) power(:)
    integer direction, ie
  end subroutine
end interface

interface
  subroutine check_end (point, ix, string_in)
    import
    implicit none
    type (outline_pt_struct) point(:)
    integer ix
    character(*) string_in
  end subroutine
end interface

interface
  subroutine check_wall (wall, s_lat, lat_type)
    import
    implicit none
    type (wall_struct), target :: wall
    real(rp) s_lat
    integer lat_type
  end subroutine
end interface

interface
  subroutine track_ray_to_wall (ray, walls)
    import
    implicit none
    type (ray_struct), target :: ray
    type (walls_struct), target :: walls
  end subroutine
end interface

interface
  subroutine synrad_adjust_wall_points (wall, branch)
    import
    implicit none
    type (wall_struct), target :: wall
    type (branch_struct), target :: branch
  end subroutine
end interface

interface
  subroutine init_ray (ray, branch, ix_ele, l_offset, orb, direction)
    import
    implicit none
    type (branch_struct), target :: branch
    type (coord_struct) orb(0:*)
    type (ray_struct) ray
    real(rp) l_offset
    integer direction
    integer ix_ele
  end subroutine
end interface

interface
  subroutine outline_concat (outline1, outline2, outline3)
    import
    implicit none
    type (outline_struct) outline1
    type (outline_struct) outline2
    type (outline_struct) outline3
  end subroutine
end interface

interface
  subroutine outline_reverse (outline1, outline2)
    import
    implicit none
    type (outline_struct) outline1
    type (outline_struct) outline2
  end subroutine
end interface

interface
  subroutine synrad_read_vac_wall_geometry (wall_file, seg_len_max, branch, walls, err_flag, seg_len_phantom_max)
    import
    implicit none
    character(*) wall_file
    type (branch_struct) branch
    type (walls_struct), target :: walls
    real(rp) seg_len_max
    real(rp), optional :: seg_len_phantom_max
    logical, optional :: err_flag
  end subroutine
end interface

interface
  subroutine seg_power_calc (fan, i_ray, walls, wall_side, branch, gen, power)
    import
    implicit none
    type (ray_struct), target :: fan(:)
    type (walls_struct), target :: walls
    type (synrad_param_struct) gen
    type (branch_struct) :: branch
    type (ele_power_struct) power
    integer i_ray, wall_side
  end subroutine
end interface

interface
  subroutine synrad_custom_seg_calc (wall, ray, seg, frac_illum)
    import
    implicit none
    type (wall_struct) wall
    type (ray_struct) ray
    type (wall_seg_struct) seg
    real(rp) frac_illum
  end subroutine
end interface

interface
  subroutine synrad_setup_walls (walls, branch, seg_len_max, seg_len_phantom_max)
    import
    implicit none
    type (walls_struct), target :: walls
    type (branch_struct) branch
    real(rp), optional :: seg_len_phantom_max
    real(rp) seg_len_max
  end subroutine
end interface

interface
  subroutine get_initial_wall_pt (ray, wall, ix_pt)
    import
    implicit none
    type (ray_struct) ray
    type (wall_struct) wall
    integer ix_pt
  end subroutine
end interface

interface
  subroutine init_wall_ends (walls)
    import
    implicit none
    type (walls_struct) walls
  end subroutine
end interface

! This is to suppress the ranlib "has no symbols" message
integer, private :: private_dummy

end module
