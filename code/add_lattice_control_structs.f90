!+
! Subroutine add_lattice_control_structs (lat, ix_ele)
! 
! Subroutine to adjust the control structure of a lat so that extra control
! elements can be added.
!
! Modules to use:
!   use bmad
!
! Input:
!   lat   -- lat_struct: lat whose control structure needs fixing
!     %ele(ix_ele)  -- This could be a new element or an existing element
!                       that needs more control info.
!     %ele(ix_ele)%n_slave -- Increase this to reserve more room in the
!                              lat%control(:) array.
!     %ele(ix_ele)%n_lord  -- Increase this to reserve more room in the 
!                              lat%ic(:) array.
!   ix_ele -- Integer: Index of element that needs extra control elements.
!
! Output:
!   lat -- lat_struct: Lat with control structure fixed.
!-

#include "CESR_platform.inc"

subroutine add_lattice_control_structs (lat, ix_ele)

  use bmad_struct
  use bmad_interface, except => add_lattice_control_structs

  implicit none

  type (lat_struct), target :: lat
  type (ele_struct), pointer :: ele

  integer ix_ele, n_add, n_con, i2, n_con2, n_ic, n_ic2


! fix slave problems

  ele => lat%ele(ix_ele)
  n_add = ele%n_slave - (ele%ix2_slave - ele%ix1_slave + 1) 

  if (n_add < 0) then
    print *, 'ERROR IN ADD_LATTICE_CONTROL_STRUCTS: N_SLAVE < CURRENT ALLOCATION'
    print *, '      FOR: ', ele%name, ele%n_slave, ele%ix2_slave, ele%ix1_slave
    call err_exit
  endif

  if (n_add > 0) then

    n_con = lat%n_control_max
    i2 = ele%ix2_slave
    n_con2 = lat%n_control_max + n_add
    if (n_con2 > size(lat%control)) &
                        lat%control => reallocate(lat%control, n_con2+500)
    lat%control(lat%n_control_max+1:) = control_struct(0.0_rp, 0, 0, 0)

    if (i2 < 0) then
      ele%ix1_slave = n_con + 1
      ele%ix2_slave = n_con + n_add
    else
      lat%control(i2+1+n_add:n_con+n_add) = lat%control(i2+1:n_con)
      lat%control(i2+1:i2+n_add)%ix_lord = ix_ele
      lat%control(i2+1:i2+n_add)%ix_slave = 0
      lat%control(i2+1:i2+n_add)%ix_attrib = 0
      lat%control(i2+1:i2+n_add)%coef = 0
      where (lat%ele%ix1_slave > i2) lat%ele%ix1_slave = &
                                            lat%ele%ix1_slave + n_add
      where (lat%ele%ix2_slave >= i2) lat%ele%ix2_slave = &
                                            lat%ele%ix2_slave + n_add
      where (lat%ic > i2) lat%ic = lat%ic + n_add
    endif

    lat%n_control_max = n_con2

  endif
                                        
! fix lord problems

  n_add = ele%n_lord - (ele%ic2_lord - ele%ic1_lord + 1) 

  if (n_add < 0) then
    print *, 'ERROR IN ADD_LATTICE_CONTROL_STRUCTS: N_LORD < CURRENT ALLOCATION'
    print *, '      FOR: ', ele%name, ele%n_lord, ele%ic2_lord, ele%ic1_lord
    call err_exit
  endif

  if (n_add > 0) then

    n_ic = lat%n_ic_max
    n_ic2 = n_ic + n_add
    if (n_ic2 > size(lat%ic)) call re_associate(lat%ic, n_ic2+500)

    i2 = ele%ic2_lord

    if (i2 < 0) then
      ele%ic1_lord = n_ic + 1
      ele%ic2_lord = n_ic + n_add
    else
      lat%ic(i2+1+n_add:n_ic+n_add) = lat%ic(i2+1:n_ic)
      lat%ic(i2+1:i2+n_add) = 0
      where (lat%ele%ic1_lord > i2) lat%ele%ic1_lord = &
                                            lat%ele%ic1_lord + n_add
      where (lat%ele%ic2_lord >= i2) lat%ele%ic2_lord = &
                                            lat%ele%ic2_lord + n_add
    endif

    lat%n_ic_max = n_ic + n_add

  endif

end subroutine
