!+
! subroutine tao_set_data_useit_opt (data)
!
! Routine to set the %data%useit_opt logicals which determine which data can be used by an optimizer.
!
! Input:
!   data(:)   -- Tao_data_struct, optional: Data to work on.
!                 Default is all data in all universes.
!-

subroutine tao_set_data_useit_opt (data)

use tao_interface, dummy => tao_set_data_useit_opt

implicit none

type (tao_data_struct), optional :: data(:)

integer i, j

! 

if (present(data)) then
  call set_this_data (data)
else
  do i = lbound(s%u, 1), ubound(s%u, 1)
    call set_this_data (s%u(i)%data)
  enddo
endif

!-----------------------------------------------------------------------
contains

subroutine set_this_data (d)

type (tao_data_struct) :: d(:)
integer ix_uni

!

if (size(d) == 0) return
ix_uni = d(1)%d1%d2%ix_universe

if (s%u(ix_uni)%is_on) then
  d%useit_opt = d%good_opt .and. d%exists .and. d%good_user .and. d%good_meas
  if (s%global%opt_with_ref) then
    d%useit_opt = d%useit_opt .and. d%good_ref
  endif
else   ! data in off universes do not get used in optimizations.
  d%useit_opt = .false.
endif

end subroutine set_this_data

end subroutine tao_set_data_useit_opt
