#include "CESR_platform.inc"

module custom_bmad_interface

interface
  subroutine emit_calc_custom (lat, ir, i2, i3, i5a, i5b)
    use bmad_struct, only: lat_struct, rp
    type (lat_struct) lat
    integer ir
    real(rp) i2, i3, i5a, i5b
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
  subroutine make_mat6_custom (ele, param, start, end)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (ele_struct), target :: ele
    type (coord_struct) :: start, end
    type (lat_param_struct) param
  end subroutine
end interface

interface
  subroutine track1_custom (start, ele, param, end)
    use bmad_struct, only: ele_struct, coord_struct, lat_param_struct
    implicit none
    type (coord_struct) :: start
    type (coord_struct) :: end
    type (ele_struct) :: ele
    type (lat_param_struct) :: param
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
