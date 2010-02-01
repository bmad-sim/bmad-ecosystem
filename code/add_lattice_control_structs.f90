!+
! Subroutine add_lattice_control_structs (lat, ele)
! 
! Subroutine to adjust the control structure of a lat so that extra control
! elements can be added.
!
! Modules to use:
!   use bmad
!
! Input:
!   lat   -- lat_struct: lat whose control structure needs fixing
!   ele -- Ele_struct: Element that needs extra control elements.
!          This could be a new element or an existing element that needs more control info.
!     %n_slave -- Increase this to reserve more room in the
!                              lat%control(:) array.
!     %n_lord  -- Increase this to reserve more room in the 
!                              lat%ic(:) array.
!
! Output:
!   lat -- lat_struct: Lat with control structure fixed.
!-

subroutine add_lattice_control_structs (lat, ele)

  use bmad_struct
  use bmad_interface, except_dummy => add_lattice_control_structs

  implicit none

  type (lat_struct), target :: lat
  type (ele_struct) ele
  type (branch_struct), pointer :: branch

  integer n_add, n_con, i2, n_con2, n_ic, n_ic2, ib

  ! fix slave problems

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
    if (n_con2 > size(lat%control)) call reallocate_control(lat, nint(1.2*n_con2)+10)
    lat%control(lat%n_control_max+1:) = control_struct(0.0_rp, 0, 0, 0, 0)

    ! If no existing slaves of this lord then just put the new lat%control elements
    ! at the end of the array.

    if (i2 < 0) then  
      ele%ix1_slave = n_con + 1
      ele%ix2_slave = n_con + n_add
      lat%control(n_con+1:n_con+n_add)%ix_lord = ele%ix_ele

    ! Else we need to make room in lat%control for the new slaves by moving
    ! a slice of lat%control.

    else
      lat%control(i2+1+n_add:n_con+n_add) = lat%control(i2+1:n_con)
      lat%control(i2+1:i2+n_add)%ix_lord = ele%ix_ele
      lat%control(i2+1:i2+n_add)%ix_slave = -1
      lat%control(i2+1:i2+n_add)%ix_branch = 0
      lat%control(i2+1:i2+n_add)%ix_attrib = 0
      lat%control(i2+1:i2+n_add)%coef = 0
      do ib = 0, ubound(lat%branch, 1)
        branch => lat%branch(ib)
        where (branch%ele%ix1_slave > i2) branch%ele%ix1_slave = branch%ele%ix1_slave + n_add
        where (branch%ele%ix2_slave >= i2) branch%ele%ix2_slave = branch%ele%ix2_slave + n_add
      enddo
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

    i2 = ele%ic2_lord

    if (i2 < 0) then
      ele%ic1_lord = n_ic + 1
      ele%ic2_lord = n_ic + n_add
    else
      if (n_ic+n_add > size(lat%ic)) call re_allocate (lat%ic, nint(1.2*(n_ic+n_add)) + 10)
      lat%ic(i2+1+n_add:n_ic+n_add) = lat%ic(i2+1:n_ic)
      lat%ic(i2+1:i2+n_add) = 0
      do ib = 0, ubound(lat%branch, 1)
        branch => lat%branch(ib)
        where (branch%ele%ic1_lord > i2) branch%ele%ic1_lord = branch%ele%ic1_lord + n_add
        where (branch%ele%ic2_lord >= i2) branch%ele%ic2_lord = branch%ele%ic2_lord + n_add
      enddo
    endif

    lat%n_ic_max = n_ic + n_add

  endif

end subroutine
