!+
! Subroutine tao_hook_does_data_exist (datum)
!
! This routine is for use with custom data types. 
! This routine sets datum%exists to true for datums with no
! associated datum%ele_name. This is necessary since otherwise
! tao_init_global_and_universes will flag the datum as
! not existing.
!
! Input:
!   datum        -- tao_data_struct: the current datum to examine
!-

subroutine tao_hook_does_data_exist (datum)

use tao_struct
use tao_interface

implicit none

type (tao_data_struct) datum

select case(datum%data_type)

!   case('data_type')
!   datum%exists = .true.

   case default
   !do nothing

end select

end subroutine
