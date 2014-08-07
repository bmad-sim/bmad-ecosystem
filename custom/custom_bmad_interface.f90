module custom_bmad_interface

interface 
  subroutine check_aperture_limit_custom (orb, ele, particle_at, param, err_flag)
    use bmad_struct, only: coord_struct, ele_struct, lat_param_struct
    implicit none
    type (coord_struct) :: orb
    type (ele_struct) :: ele
    type (lat_param_struct) :: param
    integer particle_at
    logical err_flag
  end subroutine

  subroutine wall_hit_handler_custom (orb, ele, s, t)
    use bmad_struct, only: coord_struct, ele_struct, rp
    implicit none
    type (coord_struct) :: orb
    type (ele_struct) :: ele
    real(rp) s, t
  end subroutine

  subroutine em_field_custom (ele, param, s_rel, t_rel, orb, local_ref_frame, field, calc_dfield, err_flag)
    use bmad_struct
    implicit none
    type (ele_struct) :: ele
    type (lat_param_struct) param
    type (coord_struct), intent(in) :: orb
    real(rp), intent(in) :: s_rel, t_rel
    logical local_ref_frame
    type (em_field_struct) :: field
    logical, optional :: err_flag
    logical, optional :: calc_dfield
  end subroutine

  subroutine radiation_integrals_custom (lat, ir, orb, err_flag)
    use bmad_struct, only: lat_struct, coord_struct
    implicit none
    type (lat_struct) lat
    type (coord_struct) orb(0:)
    integer ir
    logical err_flag
  end subroutine

  subroutine init_custom (ele, err_flag)
    use bmad_struct, only: ele_struct
    implicit none
    type (ele_struct), target :: ele
    logical err_flag
  end subroutine

  subroutine make_mat6_custom (ele, param, start_orb, end_orb, err_flag)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (ele_struct), target :: ele
    type (coord_struct) :: start_orb, end_orb
    type (lat_param_struct) param
    logical err_flag
  end subroutine

  subroutine make_mat6_custom2 (ele, param, start_orb, end_orb, err_flag)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (ele_struct), target :: ele
    type (coord_struct) :: start_orb, end_orb
    type (lat_param_struct) param
    logical err_flag
  end subroutine

  subroutine track1_custom (start_orb, ele, param, end_orb, track, err_flag)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct, track_struct
    implicit none
    type (coord_struct) :: start_orb
    type (coord_struct) :: end_orb
    type (ele_struct) :: ele
    type (lat_param_struct) :: param
    logical err_flag
    type (track_struct), optional :: track
  end subroutine

  subroutine track1_custom2 (start_orb, ele, param, end_orb, track, err_flag)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct, track_struct
    implicit none
    type (coord_struct) :: start_orb
    type (coord_struct) :: end_orb
    type (ele_struct) :: ele
    type (lat_param_struct) :: param
    logical err_flag
    type (track_struct), optional :: track
  end subroutine

  subroutine track1_bunch_custom (bunch_start, lat, ele, bunch_end, err_flag)
    use bmad_struct, only: lat_struct, ele_struct
    use beam_def_struct, only: bunch_struct
    implicit none
    type (bunch_struct) bunch_start, bunch_end
    type (lat_struct), target :: lat
    type (ele_struct) :: ele
    logical err_flag
  end subroutine

  subroutine track1_postprocess (start_orb, ele, param, end_orb)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct, track_struct
    implicit none
    type (coord_struct) :: start_orb
    type (coord_struct) :: end_orb
    type (ele_struct) :: ele
    type (lat_param_struct) :: param
  end subroutine

  subroutine track1_spin_custom (start_orb, ele, param, end_orb, err_flag, track)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct, track_struct
    implicit none
    type (coord_struct) :: start_orb
    type (coord_struct) :: end_orb
    type (ele_struct) :: ele
    type (lat_param_struct) :: param
    type (track_struct), optional :: track
    logical err_flag
  end subroutine
end interface

! This is to suppress the ranlib "has no symbols" message
integer, private :: private_dummy

end module
