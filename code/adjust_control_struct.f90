!+
! Subroutine adjust_control_struct (ring, ix_ele)
! 
! Subroutine to adjust the control structure of a ring so that extra control
! elements can be added.
!
! Modules to use:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ring   -- Ring_struct: ring whose control structure needs fixing
!     %ele_(ix_ele)%n_slave -- Increase this to reserve more room in the
!                              ring%control_(:) array.
!     %ele_(ix_ele)%n_lord  -- Increase this to reserve more room in the 
!                              ring%ic_(:) array.
!   ix_ele -- Integer: Index of element that needs extra control elements.
!
! Output:
!   ring -- Ring_struct: Ring with control structure fixed.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:47  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine adjust_control_struct (ring, ix_ele)

  use bmad_struct
                                                 
  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele

  integer ix_ele, n_add, iz, i2, ic


! fix slave problems

  ele => ring%ele_(ix_ele)
  n_add = ele%n_slave - (ele%ix2_slave - ele%ix1_slave + 1) 

  if (n_add < 0) then
    type *, 'ERROR IN ADJUST_CONTROL_STRUCT: N_SLAVE < CURRENT ALLOCATION'
    type *, '      FOR: ', ele%name, ele%n_slave, ele%ix2_slave, ele%ix1_slave
    call err_exit
  endif

  if (n_add > 0) then
    iz = ring%n_control_array
    i2 = ele%ix2_slave

    if (i2 < 0) then
      ele%ix1_slave = iz + 1
      ele%ix2_slave = iz + n_add
    else
      ring%control_(i2+1+n_add:iz+n_add) = ring%control_(i2+1:iz)
      ring%control_(i2+1:i2+n_add)%ix_lord = ix_ele
      ring%control_(i2+1:i2+n_add)%ix_slave = 0
      ring%control_(i2+1:i2+n_add)%ix_attrib = 0
      ring%control_(i2+1:i2+n_add)%coef = 0
      where (ring%ele_%ix1_slave > i2) ring%ele_%ix1_slave = &
                                            ring%ele_%ix1_slave + n_add
      where (ring%ele_%ix2_slave >= i2) ring%ele_%ix2_slave = &
                                            ring%ele_%ix2_slave + n_add
      where (ring%ic_ > i2) ring%ic_ = ring%ic_ + n_add
    endif

    ring%n_control_array = iz + n_add

  endif
                                        
! fix lord problems

  n_add = ele%n_lord - (ele%ic2_lord - ele%ic1_lord + 1) 

  if (n_add < 0) then
    type *, 'ERROR IN ADJUST_CONTROL_STRUCT: N_LORD < CURRENT ALLOCATION'
    type *, '      FOR: ', ele%name, ele%n_lord, ele%ic2_lord, ele%ic1_lord
    call err_exit
  endif

  if (n_add > 0) then
    ic = ring%n_ic_array
    i2 = ele%ic2_lord

    if (i2 < 0) then
      ele%ic1_lord = ic + 1
      ele%ic2_lord = ic + n_add
    else
      ring%ic_(i2+1+n_add:ic+n_add) = ring%ic_(i2+1:ic)
      ring%ic_(i2+1:i2+n_add) = 0
      where (ring%ele_%ic1_lord > i2) ring%ele_%ic1_lord = &
                                            ring%ele_%ic1_lord + n_add
      where (ring%ele_%ic2_lord >= i2) ring%ele_%ic2_lord = &
                                            ring%ele_%ic2_lord + n_add
    endif

    ring%n_ic_array = ic + n_add

  endif

end subroutine
