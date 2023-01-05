!+
! Subroutine tao_lat_bookkeeper (u, tao_lat, err_flag)
!
! This will make sure all bookkeeping is up to date.
!
! Input:
!  u            -- tao_universe_struct:
!  tao_lat      -- tao_lattice_struct: 
!
! Output:
!  err_flag     -- logical: Set True if there is a problem. False otherwise.
!-

subroutine tao_lat_bookkeeper (u, tao_lat, err_flag)

use bmad
use tao_struct

implicit none

type (tao_universe_struct), target :: u
type (tao_lattice_struct) :: tao_lat

integer i, j
logical err_flag

character(*), parameter :: r_name = "tao_lat_bookkeeper"

! Do bookkeeping

call lattice_bookkeeper (tao_lat%lat, err_flag)

if (err_flag .and. s%com%optimizer_running) then
  do i = 1, size(u%data)
    if (u%data(i)%data_type == 'unstable.lattice') return
  enddo

  call out_io (s_error$, r_name, 'CANNOT DO LATTICE BOOKKEEPING.', &
                                 'ADD AN "unstable.lattice" DATUM TO DRIVE THE LATTICE TOWARDS STABILITY.')
endif

end subroutine tao_lat_bookkeeper

