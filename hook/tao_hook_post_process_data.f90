!+
! Subroutine tao_hook_post_process_data ()
!
! Here place anything that needs to be done to the data, or anything else in the
! superuniverse after the data arrays have been loaded.
!
! For example, BPM resolution can be handled here.
!
!-


subroutine tao_hook_post_process_data ()

use bmad_struct
use bmad_interface
use tao_mod
use tao_data_mod
implicit none

type(ring_struct), pointer :: lat
type(ele_struct), pointer :: ele
type(modes_struct), pointer :: mode
!type(tao_real_array_struct), pointer :: r(:)
type(tao_data_array_struct), allocatable :: d(:)
integer i,j
logical err
character data_name*6
real(rp) datum_value
integer foo1

lat=>s%u(1)%model%lat
mode=>s%u(1)%model%modes

data_name = 'ibsT.x'
call tao_find_data(err,data_name, d_array=d)
if(err .eq. .true.) THEN
  WRITE(*,*) "tao_find_data error!"
  STOP
endif

do i = 1, size(d)-1
  ele=>lat%ele_( d(i)%d%ix_ele )

  d(i)%d%model_value=(mode%a%emittance*ele%x%beta) / &
    (mode%a%emittance*ele%x%beta + &
    (ele%x%eta**2)*(mode%sigE_E**2))

enddo

end subroutine tao_hook_post_process_data
