!+
! Subroutine tao_hook_merit_data (i_uni, j_data, data)
! 
! Dummy routine that needs to be over written in order to implement a
! custom merit calculation for data .
!
! Input:
!   i_uni   -- Integer: Universe index for the datiable in the s%u(:) array.
!   j_data  -- Integer: Index of the data in the u%data(:) array.
!   data    -- Tao_dat_struct: Data whose contribution to the merit function
!               is to be calculated.
!
! Output:
!   data     -- Tao_dat_struct: Data whose contribution to the merit function
!                is to be calculated.
!     %merit -- Real(rp): Contribution to the merit function.
!-

subroutine tao_hook_merit_data (i_uni, j_data, data)

use tao_mod

implicit none

type (tao_data_struct) data

integer, intent(in) :: i_uni, j_data
character(20) :: r_name = 'tao_hook_merit_data'

!

call out_io (s_error$, r_name, 'THIS ROUTINE SHOULD NOT HAVE BEEN CALLED')
call out_io (s_error$, r_name, 'MERIT_TYPE NOT RECOGNIZED FOR DATA: ' &
        // data%name, 'MERIT_TYPE: ' // data%merit_type)

stop

end subroutine