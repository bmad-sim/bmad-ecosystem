!+
! Function num_lords (slave, lord_type) result (num)
!
! Routine to return the number of lords of a given type for a given lattice element.
!
! Input:
!   slave     -- ele_struct: Slave element.
!   lord_type -- integer: Type of lord. super_lord$, multipass_lord$, girder_lord$, 
!                   group_lord$, overlay_lord$, and governor$ (= group + overlay + control + girder)
!
! Output:
!   num       -- integer: Number of lords of the given type.
!-

function num_lords (slave, lord_type) result (num)

use bmad_interface, except_dummy => num_lords

implicit none

type (ele_struct), target :: slave
type (ele_struct), pointer :: lord
integer lord_type, num
integer i

!

num = 0
do i = 1, slave%n_lord
  lord => pointer_to_lord(slave, i)
  if (lord_type == governor$) then
    select case (lord%lord_status)
    case (group_lord$, overlay_lord$, control_lord$, girder_lord$); num = num + 1
    end select
  else
    if (lord%lord_status == lord_type) num = num + 1
  endif
enddo

end function num_lords
