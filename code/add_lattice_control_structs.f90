!+
! Subroutine add_lattice_control_structs (ring, ix_ele)
! 
! Subroutine to adjust the control structure of a ring so that extra control
! elements can be added.
!
! Modules to use:
!   use bmad
!
! Input:
!   ring   -- Ring_struct: ring whose control structure needs fixing
!     %ele_(ix_ele)  -- This could be a new element or an existing element
!                       that needs more control info.
!     %ele_(ix_ele)%n_slave -- Increase this to reserve more room in the
!                              ring%control_(:) array.
!     %ele_(ix_ele)%n_lord  -- Increase this to reserve more room in the 
!                              ring%ic_(:) array.
!   ix_ele -- Integer: Index of element that needs extra control elements.
!
! Output:
!   ring -- Ring_struct: Ring with control structure fixed.
!-

#include "CESR_platform.inc"

subroutine add_lattice_control_structs (ring, ix_ele)

  use bmad_struct
  use bmad_interface, except => add_lattice_control_structs

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele

  integer ix_ele, n_add, n_con, i2, n_con2, n_ic, n_ic2


! fix slave problems

  ele => ring%ele_(ix_ele)
  n_add = ele%n_slave - (ele%ix2_slave - ele%ix1_slave + 1) 

  if (n_add < 0) then
    print *, 'ERROR IN ADD_LATTICE_CONTROL_STRUCTS: N_SLAVE < CURRENT ALLOCATION'
    print *, '      FOR: ', ele%name, ele%n_slave, ele%ix2_slave, ele%ix1_slave
    call err_exit
  endif

  if (n_add > 0) then

    n_con = ring%n_control_max
    i2 = ele%ix2_slave
    n_con2 = ring%n_control_max + n_add
    if (n_con2 > size(ring%control_)) &
                        ring%control_ => reallocate(ring%control_, n_con2+500)
    ring%control_(ring%n_control_max+1:) = control_struct(0.0_rp, 0, 0, 0)

    if (i2 < 0) then
      ele%ix1_slave = n_con + 1
      ele%ix2_slave = n_con + n_add
    else
      ring%control_(i2+1+n_add:n_con+n_add) = ring%control_(i2+1:n_con)
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

    ring%n_control_max = n_con2

  endif
                                        
! fix lord problems

  n_add = ele%n_lord - (ele%ic2_lord - ele%ic1_lord + 1) 

  if (n_add < 0) then
    print *, 'ERROR IN ADD_LATTICE_CONTROL_STRUCTS: N_LORD < CURRENT ALLOCATION'
    print *, '      FOR: ', ele%name, ele%n_lord, ele%ic2_lord, ele%ic1_lord
    call err_exit
  endif

  if (n_add > 0) then

    n_ic = ring%n_ic_max
    n_ic2 = n_ic + n_add
    if (n_ic2 > size(ring%ic_)) call re_associate(ring%ic_, n_ic2+500)

    i2 = ele%ic2_lord

    if (i2 < 0) then
      ele%ic1_lord = n_ic + 1
      ele%ic2_lord = n_ic + n_add
    else
      ring%ic_(i2+1+n_add:n_ic+n_add) = ring%ic_(i2+1:n_ic)
      ring%ic_(i2+1:i2+n_add) = 0
      where (ring%ele_%ic1_lord > i2) ring%ele_%ic1_lord = &
                                            ring%ele_%ic1_lord + n_add
      where (ring%ele_%ic2_lord >= i2) ring%ele_%ic2_lord = &
                                            ring%ele_%ic2_lord + n_add
    endif

    ring%n_ic_max = n_ic + n_add

  endif

end subroutine
