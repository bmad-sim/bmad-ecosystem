module sr_interface

  use sr_struct

  interface
    subroutine write_power_results (wall, ring, gen_params)
      use sr_struct
      implicit none
      type (wall_struct), target :: wall
      type (general_param_struct) gen_params
      type (ring_struct) ring
    end subroutine
  end interface

  interface
    subroutine break_wall_into_segments (wall, seg_len_max)
      use sr_struct
      implicit none
      type (wall_struct) wall
      real(rp) seg_len_max
    end subroutine
  end interface

  interface
    subroutine calculate_sr_power (ring, orb, direction, power, &
                                                     inside, outside, gen)
      use sr_struct
      implicit none
      type (ring_struct), target :: ring
      type (coord_struct) orb(0:*)
      type (wall_struct) inside
      type (wall_struct) outside
      type (general_param_struct) gen
      type (ele_power_struct) power(*)
      integer direction
    end subroutine
  end interface

  interface
    subroutine ele_sr_power (ring, ie, orb, direction, power, &
                                                     inside, outside, gen)
      use sr_struct
      implicit none
      type (ring_struct), target :: ring
      type (coord_struct) orb(0:*)
      type (wall_struct) inside
      type (wall_struct) outside
      type (general_param_struct) gen
      type (ele_power_struct) power(*)
      integer direction, ie
    end subroutine
  end interface

  interface
    subroutine hit_spot_calc (ray, wall, ix_wall, has_hit, ring)
      use sr_struct
      implicit none
      type (ring_struct) ring
      type (ray_struct) :: ray
      type (wall_struct), target :: wall
      integer ix_wall
      logical has_hit
    end subroutine
  end interface

  interface
    subroutine check_end (point_, ix, string_in)
      use sr_struct
      implicit none
      type (outline_pt_struct) point_(*)
      integer ix
      character*(*) string_in
    end subroutine
  end interface

  interface
    subroutine check_wall (wall, ring)
      use sr_struct
      implicit none
      type (wall_struct) wall
      type (ring_struct) ring
    end subroutine
  end interface

  interface
    subroutine track_ray_to_wall (ray, ring, inside, outside, &
                                               hit_flag, track_max)
      use sr_struct
      implicit none
      type (ring_struct), target :: ring
      type (ray_struct), target :: ray
      type (wall_struct) inside
      type (wall_struct) outside
      logical, optional :: hit_flag
      real(rp), optional :: track_max
    end subroutine
  end interface

  interface
    subroutine convert_blanks_to_underscore (string_in, string_out)
      use sr_struct
      implicit none
      character*(*) string_in
      character*(*) string_out
    end subroutine
  end interface

  interface
    subroutine create_alley (wall)
      use sr_struct
      implicit none
      type (wall_struct), target :: wall
    end subroutine
  end interface

  interface
    subroutine delete_overlaping_wall_points (wall)
      use sr_struct
      implicit none
      type (wall_struct), target :: wall
    end subroutine
  end interface

  interface
    subroutine init_ray (ray, ring, ix_ele, l_offset, orb, direction)
      use sr_struct
      implicit none
      type (ring_struct), target :: ring
      type (coord_struct) orb(0:*)
      type (ray_struct) ray
      real(rp) l_offset
      integer direction
      integer ix_ele
    end subroutine
  end interface

  interface
    subroutine outline_concat (outline1, outline2, outline3)
      use sr_struct
      implicit none
      type (outline_struct) outline1
      type (outline_struct) outline2
      type (outline_struct) outline3
    end subroutine
  end interface

  interface
    subroutine outline_reverse (outline1, outline2)
      use sr_struct
      implicit none
      type (outline_struct) outline1
      type (outline_struct) outline2
    end subroutine
  end interface

  interface
    subroutine propagate_ray (ray, s_end, ring)
      use sr_struct
      implicit none
      type (ring_struct), target :: ring
      type (ray_struct), target :: ray
      real(rp) s_end
    end subroutine
  end interface

  interface
    subroutine seg_power_calc (rays, i_ray, inside, outside, ring, gen, power)
      use sr_struct
      implicit none
      type (ray_struct) :: rays(*)
      type (wall_struct)          inside
      type (wall_struct)          outside
      type (general_param_struct) gen
      type (ring_struct)          ring
      type (ele_power_struct)     power
      integer i_ray
    end subroutine
  end interface

  interface
    subroutine get_initial_pt (ray, wall, ix_wall, ring)
      use sr_struct
      implicit none
      type (ray_struct) ray
      type (wall_struct) wall
      type (ring_struct) ring
      integer ix_wall
    end subroutine
  end interface

  interface
    subroutine init_wall (wall)
      use sr_struct
      implicit none
      type (wall_struct) wall
    end subroutine
  end interface

  interface
    subroutine next_pt (ray, wall, ix_wall)
      use sr_struct
      implicit none
      type (ray_struct) ray
      type (wall_struct) wall
      integer ix_wall
    end subroutine
  end interface

  interface
    subroutine write_power_header (iu, file, gen_params)
      use sr_struct
      implicit none
      type (general_param_struct) gen_params
      character*(*) file
      integer iu
    end subroutine
  end interface

end module
