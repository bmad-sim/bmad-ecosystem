module da_program_mod

use bmad
use dynamic_aperture_mod

implicit none


type da_common_struct
  ! User settable
  logical :: ramping_on = .false.
  real(rp) :: ramping_start_time = 0
  ! Internal params
  integer :: n_ramper_loc = 0
  type (lat_ele_loc_struct), allocatable :: ramper(:)    ! Ramper element locations.
end type


type (da_common_struct), save :: da_com

!

end module
