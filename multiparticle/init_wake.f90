!+
! Subroutine init_wake (wake, n_sr_long, n_sr_trans, n_lr_mode, n_lr_position_array, always_allocate)
!
! Subroutine to initialize a wake struct.
! If the wake is allocated, All components are always allocated even when the size is zero.
!
! Input:
!   n_sr_long           -- Integer: Number of terms: wake%nr(n_sr_long).
!   n_sr_trans          -- Integer: Number of terms: wake%nr(n_sr_trans).
!   n_lr_mode           -- Integer: Number of terms: wake%nr(n_lr_mode)
!   n_lr_position_array -- Integer: Number of terms: wake%nr_formula(n_lr_position_array)
!   always_allocate     -- logical, optional: If present and True then allways allocate wake
!                           even if n_lr_mode, etc. are all 0. Default is False.
!
! Output:
!   wake -- Wake_struct, pointer: Initialized structure. 
!-

subroutine init_wake (wake, n_sr_long, n_sr_trans, n_lr_mode, n_lr_position_array, always_allocate)

use bmad_struct

implicit none

type (wake_struct), pointer :: wake
integer n_sr_long, n_sr_trans, n_lr_mode, n_lr_position_array
integer i
logical, optional :: always_allocate

! Deallocate wake if all inputs are zero.

if (n_sr_long == 0 .and. n_sr_trans == 0 .and. n_lr_mode == 0 .and. n_lr_position_array == 0) then
  if (logic_option(.false., always_allocate)) then
    if (.not. associated(wake)) allocate (wake)
    allocate (wake%sr_long%mode(n_sr_long))
    allocate (wake%sr_trans%mode(n_sr_trans))
    allocate (wake%lr_mode(n_lr_mode))
    allocate (wake%lr_position_array(n_lr_position_array))
  else
    if (associated(wake)) deallocate (wake)
  endif
  return
endif

!

if (associated (wake)) then
  if (size(wake%sr_long%mode) /= n_sr_long) then
    deallocate (wake%sr_long%mode)
    allocate (wake%sr_long%mode(n_sr_long))
  endif

  if (size(wake%sr_trans%mode) /= n_sr_trans) then
    deallocate (wake%sr_trans%mode)
    allocate (wake%sr_trans%mode(n_sr_trans))
  endif

  if (size(wake%lr_mode) /= n_lr_mode) then
    deallocate (wake%lr_mode)
    allocate (wake%lr_mode(n_lr_mode))
  endif

  if (size(wake%lr_position_array) /= n_lr_position_array) then
    deallocate (wake%lr_position_array)
    allocate (wake%lr_position_array(n_lr_position_array))
  endif

else
  allocate (wake)
  allocate (wake%sr_long%mode(n_sr_long))
  allocate (wake%sr_trans%mode(n_sr_trans))
  allocate (wake%lr_mode(n_lr_mode))
  allocate (wake%lr_position_array(n_lr_position_array))
endif

end subroutine init_wake

