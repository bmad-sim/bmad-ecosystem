!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_lat_bookkeeper (u, tao_lat)
!
! This will make sure all bookkeeping is up to date.
!
! Input:
!  u            -- tao_universe_struct
!  lat_name     -- Integer: Which lattice
!
! Output:
!  tao_lat      -- lat_struct
!-

subroutine tao_lat_bookkeeper (u, tao_lat)

use bmad
use tao_struct

implicit none

type (tao_universe_struct), target :: u
type (tao_lattice_struct) :: tao_lat

integer i, j

character(20) :: r_name = "tao_lat_bookkeeper"

! Setup from common if it exists

if (associated(u%common)) then

  ! First put in the common values

  do i = 1, s%n_var_used
    do j = 1, size(s%var(i)%slave)
      s%var(i)%slave(j)%model_value = s%var(i)%common_slave%model_value
      s%var(i)%slave(j)%base_value  = s%var(i)%common_slave%base_value
    enddo
  enddo

  ! Then put in the values for this universe

  do i = 1, s%n_var_used
    do j = 1, size(s%var(i)%slave)
      if (s%var(i)%slave(j)%ix_uni /= u%ix_uni) cycle
      s%var(i)%slave(j)%model_value = s%var(i)%model_value
      s%var(i)%slave(j)%base_value = s%var(i)%base_value
    enddo
  enddo

endif

! Do bookkeeping

call lattice_bookkeeper (tao_lat%lat)

end subroutine tao_lat_bookkeeper

