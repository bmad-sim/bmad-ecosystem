!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_data_check (err)
!
! Routine to do some checking of data.
!-

subroutine tao_data_check (err)

use tao_interface, dummy => tao_data_check

implicit none

type (tao_data_struct), pointer :: datum
integer iu, id
logical err
character(16) :: r_name = 'tao_data_check'

!

err = .false.

do iu = lbound(s%u, 1), ubound(s%u, 1)
  do id = 1, size(s%u(iu)%data)
    datum => s%u(iu)%data(id)
    if (datum%merit_type(1:4) == 'int_' .and. &
            (s%global%opt_with_ref .or. s%global%opt_with_base)) then
      call out_io (s_error$, r_name, &
                        'BAD DATUM INTEGRATION FOR: ' // tao_datum_name(datum))
      err = .true.
    endif
  enddo
enddo

end subroutine tao_data_check
