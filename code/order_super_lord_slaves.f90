!+
! Subroutine order_super_lord_slaves (ring, ix_lord)
!
! Subroutine to make the slave elements of a super_lord in order.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ring    -- Ring_struct: Ring.
!   ix_lord -- Integer: Index of lord element.
!
! Output
!   ring -- Ring_struct: Ring with fixed controls.
!-

!$Id$
!$Log$
!Revision 1.3  2002/01/08 21:44:42  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:56  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine order_super_lord_slaves (ring, ix_lord)

  use bmad_struct
  use bmad_interface
  use nr

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele
  type (control_struct), allocatable :: cs_(:)

  integer i, ix, ix_lord, ix1, ix2, ns
  integer, allocatable :: ixx(:), iyy(:)

  real, allocatable :: s_rel(:)

!

  ele => ring%ele_(ix_lord)
  ix1 = ele%ix1_slave; ix2 = ele%ix2_slave

  if (ele%control_type /= super_lord$) then
    type *, 'ERROR IN ORDER_SUPER_LORD_SLAVES: ELEMENT NOT A SUPER_LORD'
    call err_exit
  endif

!

  ns = ele%n_slave
  allocate (s_rel(ns), ixx(ns), iyy(ns), cs_(ns))

  do i = ix1, ix2
    ix = ring%control_(i)%ix_slave
    s_rel(i+1-ix1) = ring%ele_(ix)%s - ele%s
  enddo

  where (s_rel > ring%param%total_length / 2) &
                    s_rel = s_rel - ring%param%total_length

  where (s_rel < ring%param%total_length / 2) &
                    s_rel = s_rel + ring%param%total_length

  call indexx (s_rel, ixx)
  cs_ = ring%control_(ix1:ix2) 

  do i = 1, ns
    ring%control_(i+ix1-1) = cs_(ixx(i))
    iyy(ixx(i)) = i
  enddo
  
  do i = 1, ring%n_ic_array
    if (ring%ic_(i) >= ix1 .and. ring%ic_(i) <= ix2) then
      ring%ic_(i) = iyy(ring%ic_(i)+1-ix1) + ix1 - 1
    endif
  enddo

  deallocate (s_rel)

end subroutine
