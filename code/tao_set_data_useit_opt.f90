!+
! subroutine tao_set_data_useit_opt (s)
!
! Routine to set the %data%useit_opt logicals which
! determine which data will be used by an optimizer.
!
! Input/Output:
!   s 		-- type_super_universe_struct
!-

subroutine tao_set_data_useit_opt (s)

use tao_mod

implicit none

type (tao_super_universe_struct), target :: s
type (tao_data_struct), pointer :: d(:)

integer i, j, k
logical good_opt

!

do i = 1, size(s%u)

!return if no data
  if (.not. associated (s%u(i)%d2_data)) return 

  do j = 1, size(s%u(i)%d2_data)
    good_opt = s%u(i)%d2_data(j)%good_opt
    do k = 1, size(s%u(i)%d2_data(j)%d1)
      d => s%u(i)%d2_data(j)%d1(k)%d
      d(:)%useit_opt = good_opt .and. d(:)%exists .and. &
                    d(:)%good_user .and. d(:)%good_data

      if (s%global%opt_with_ref) d(:)%useit_opt = &
                            d(:)%useit_opt .and. d(:)%good_ref
    enddo
  enddo
enddo

end subroutine tao_set_data_useit_opt
