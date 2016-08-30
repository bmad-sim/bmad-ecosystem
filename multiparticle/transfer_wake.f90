!+
! Subroutine transfer_wake (wake_in, wake_out)
!
! Subroutine to transfer the wake info from one struct to another.
!
! Modules needed:
!   use bmad
!
! Input:
!   wake_in -- Wake_struct, pointer: Input wake.
!
! Output:
!   wake_out -- Wake_struct, pointer: Output wake.
!-

subroutine transfer_wake (wake_in, wake_out)

use basic_bmad_interface, except_dummy => transfer_wake

implicit none

type (wake_struct), pointer :: wake_in, wake_out
type (wake_lr_spline_struct), pointer :: lrp_in, lrp_out
integer n_sr_long, n_sr_trans, n_lr_mode, n_lr_pa, i

!

if (associated (wake_in)) then
  n_sr_long   = size(wake_in%sr_long%mode)
  n_sr_trans  = size(wake_in%sr_trans%mode)
  n_lr_mode   = size(wake_in%lr_mode)
  n_lr_pa     = size(wake_in%lr_spline)

  call init_wake (wake_out, n_sr_long, n_sr_trans, n_lr_mode, n_lr_pa, .true.)

  do i = 1, n_lr_pa
    lrp_in  => wake_in%lr_spline(i)
    lrp_out => wake_out%lr_spline(i)

    if (allocated(lrp_in%bunch) .and. .not. allocated(lrp_out%bunch)) then
      allocate (lrp_out%bunch(size(lrp_in%bunch)))
    elseif (.not. allocated(lrp_in%bunch) .and. allocated(lrp_out%bunch)) then
      deallocate(lrp_out%bunch)
    elseif (allocated(lrp_in%bunch)) then
      if (size(lrp_out%bunch) /= size(lrp_in%bunch)) then
        deallocate (lrp_out%bunch)
        allocate (lrp_out%bunch(size(lrp_in%bunch)))
      endif
    endif

    if (allocated(lrp_in%spline) .and. .not. allocated(lrp_out%spline)) then
      allocate (lrp_out%spline(size(lrp_in%spline)))
    elseif (.not. allocated(lrp_in%spline) .and. allocated(lrp_out%spline)) then
      deallocate(lrp_out%spline)
    elseif (allocated(lrp_in%spline)) then
      if (size(lrp_out%spline) /= size(lrp_in%spline)) then
        deallocate (lrp_out%spline)
        allocate (lrp_out%spline(size(lrp_in%spline)))
      endif
    endif
  enddo

  wake_out = wake_in

else
  if (associated(wake_out)) call init_wake (wake_out, 0, 0, 0, 0)
endif

end subroutine transfer_wake

