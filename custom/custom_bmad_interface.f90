module custom_bmad_interface

interface 
  subroutine check_aperture_limit_custom (orb, ele, at, param)
    use bmad_struct, only: coord_struct, ele_struct, lat_param_struct
    implicit none
    type (coord_struct) :: orb
    type (ele_struct) :: ele
    type (lat_param_struct) :: param
    integer at
  end subroutine
end interface

interface 
  subroutine em_field_custom (ele, param, s_rel, t_rel, orb, local_ref_frame, field, calc_dfield)
    use bmad_struct
    implicit none
    type (ele_struct) :: ele
    type (lat_param_struct) param
    type (coord_struct), intent(in) :: orb
    real(rp), intent(in) :: s_rel, t_rel
    logical local_ref_frame
    type (em_field_struct), intent(out) :: field
    logical, optional :: calc_dfield
  end subroutine
end interface

interface
  subroutine radiation_integrals_custom (lat, ir, orb)
    use bmad_struct, only: lat_struct, coord_struct
    implicit none
    type (lat_struct) lat
    type (coord_struct) orb(0:)
    integer ir
  end subroutine
end interface

interface
  subroutine init_custom (ele)
    use bmad_struct, only: ele_struct
    implicit none
    type (ele_struct), target :: ele
  end subroutine
end interface

interface
  subroutine make_mat6_custom (ele, param, start, end)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (ele_struct), target :: ele
    type (coord_struct) :: start, end
    type (lat_param_struct) param
  end subroutine
end interface

interface
  subroutine track1_custom (start, ele, param, end, track)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct, track_struct
    implicit none
    type (coord_struct) :: start
    type (coord_struct) :: end
    type (ele_struct) :: ele
    type (lat_param_struct) :: param
    type (track_struct), optional :: track
  end subroutine
end interface

interface
  subroutine track1_bunch_custom (bunch_start, lat, ele, bunch_end)
    use bmad_struct, only: lat_struct, ele_struct
    use beam_def_struct, only: bunch_struct
    implicit none
    type (bunch_struct) bunch_start, bunch_end
    type (lat_struct), target :: lat
    type (ele_struct) :: ele
  end subroutine
end interface

! This is to suppress the ranlib "has no symbols" message
integer, private :: private_dummy

end module
