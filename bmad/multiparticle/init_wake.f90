!+
! Subroutine init_wake (wake, n_sr_long, n_sr_trans, n_sr_z, n_lr_mode, always_allocate)
!
! Subroutine to initialize a wake struct.
! If the wake is allocated, All components are always allocated even when the size is zero.
!
! Input:
!   n_sr_long           -- Integer: Number of terms: wake%sr%long.
!   n_sr_trans          -- Integer: Number of terms: wake%sr%trans.
!   n_sr_z              -- Integer: Number of terms: wake%sr%z.
!   n_lr_mode           -- Integer: Number of terms: wake%lr%mode.
!   always_allocate     -- logical, optional: If present and True then allways allocate wake
!                           even if n_lr_mode, etc. are all 0. Default is False.
!
! Output:
!   wake -- Wake_struct, pointer: Initialized structure. 
!-

subroutine init_wake (wake, n_sr_long, n_sr_trans, n_sr_z, n_lr_mode, always_allocate)

use bmad_struct

implicit none

type (wake_struct), pointer :: wake
integer n_sr_long, n_sr_trans, n_sr_z, n_lr_mode
integer i
logical, optional :: always_allocate

! Deallocate wake if all inputs are zero.

if (n_sr_long == 0 .and. n_sr_trans == 0 .and. n_sr_z == 0 .and. n_lr_mode == 0 .and. .not. logic_option(.false., always_allocate)) then
  if (associated(wake)) deallocate (wake)
  return
endif

!

if (associated (wake)) then
  if (size(wake%sr%long) /= n_sr_long) then
    deallocate (wake%sr%long)
    allocate (wake%sr%long(n_sr_long))
  endif

  if (size(wake%sr%trans) /= n_sr_trans) then
    deallocate (wake%sr%trans)
    allocate (wake%sr%trans(n_sr_trans))
  endif

  if (size(wake%sr%z_long%w) /= n_sr_z) then
    deallocate (wake%sr%z_long%w, wake%sr%z_long%fw, wake%sr%z_long%fbunch, wake%sr%z_long%w_out)
    allocate (wake%sr%z_long%w(n_sr_z), wake%sr%z_long%fw(n_sr_z), wake%sr%z_long%fbunch(n_sr_z), wake%sr%z_long%w_out(n_sr_z))
  endif

  if (size(wake%lr%mode) /= n_lr_mode) then
    deallocate (wake%lr%mode)
    allocate (wake%lr%mode(n_lr_mode))
  endif

else
  allocate (wake)
  allocate (wake%sr%long(n_sr_long))
  allocate (wake%sr%trans(n_sr_trans))
  allocate (wake%sr%z_long%w(n_sr_z), wake%sr%z_long%fw(n_sr_z), wake%sr%z_long%fbunch(n_sr_z), wake%sr%z_long%w_out(n_sr_z))
  allocate (wake%lr%mode(n_lr_mode))
endif

end subroutine init_wake

