!+
! subroutine tao_set_data_useit_opt ()
!
! Routine to set the %data%useit_opt logicals which
! determine which data will be used by an optimizer.
!
! Input/Output:
!   s 		-- type_super_universe_struct
!-

subroutine tao_set_data_useit_opt ()

use tao_mod

implicit none

type (tao_data_struct), pointer :: d(:)

integer i, j

!

do i = 1, size(s%u)
  d => s%u(i)%data
  d(:)%useit_opt = d(:)%good_opt .and. d(:)%exists .and. &
                              d(:)%good_user .and. d(:)%good_meas
  if (s%global%opt_with_ref) d(:)%useit_opt = &
                          d(:)%useit_opt .and. d(:)%good_ref
enddo

end subroutine tao_set_data_useit_opt
