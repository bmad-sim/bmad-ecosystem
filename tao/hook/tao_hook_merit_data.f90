!+
! Subroutine tao_hook_merit_data (i_uni, j_data, data, valid_value_set)
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
!   data        -- Tao_dat_struct: Data whose contribution to the merit function
!                   is to be calculated.
!     %merit    -- Real(rp): Contribution to the merit function.
!   valid_value_set  -- logical: Set True if this routine properly handles this datum. False otherwise.
!-

subroutine tao_hook_merit_data (i_uni, j_data, data, valid_value_set)

use tao_interface, dummy => tao_hook_merit_data

implicit none

type (tao_data_struct) data

integer, intent(in) :: i_uni, j_data
character(*), parameter :: r_name = 'tao_hook_merit_data'
logical valid_value_set

!

valid_value_set = .false.

end subroutine
