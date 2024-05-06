!+
! Subroutine transfer_wake (wake_in, wake_out)
!
! Subroutine to transfer the wake info from one struct to another.
!
! Input:
!   wake_in -- Wake_struct, pointer: Input wake.
!
! Output:
!   wake_out -- Wake_struct, pointer: Output wake.
!-

subroutine transfer_wake (wake_in, wake_out)

use bmad_routine_interface, except_dummy => transfer_wake

implicit none

type (wake_struct), pointer :: wake_in, wake_out
integer n_sr_long, n_sr_trans, n_sr_time, n_lr_mode, i

!

if (associated (wake_in)) then
  n_sr_long   = size(wake_in%sr%long)
  n_sr_trans  = size(wake_in%sr%trans)
  n_sr_time   = size(wake_in%sr_time)
  n_lr_mode   = size(wake_in%lr%mode)

  call init_wake (wake_out, n_sr_long, n_sr_trans, n_sr_time, n_lr_mode, .true.)
  wake_out = wake_in

else
  if (associated(wake_out)) call init_wake (wake_out, 0, 0, 0, 0)
endif

end subroutine transfer_wake

