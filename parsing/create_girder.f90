!+
! Subroutine create_girder (lat, ix_girder, contrl, girder_info, err_flag)
!
! Subroutine to add the controller information to slave elements of
! an girder_lord.
!
! Input:
!   lat          -- lat_struct: Lat to modify.
!   ix_girder    -- Integer: Index of girder element.
!   contrl(:)    -- control_struct: Array of elements that are supported by the girder.
!     slave%ix_ele       -- Integer: Index in lat%branch()%ele() of element controlled.
!     %ix_branch      -- Integer: Branch index.  
!   girder_info  -- ele_struct: Element containing attributes to be transfered to the Girder element:
!                       girder_info%name        
!                       girder_info%alias
!                       girder_info%descrip
!                       girder_info%value(:)
!
! Output:
!   lat    -- lat_struct: Modified lattice.
!
! Example: Create the Girder supporting elements 
! lat%ele(10) and lat%ele(12)
!
!   call new_control (lat, ix_ele, ele_name)        ! get IX_ELE index of the girder.
!   control(1:2)%ix_branch = 0
!   control(1)%ix_ele = 10; control(2)%ix_ele = 12
!   call create_girder (lat, ix_ele, control(1:2))  ! create the girder
!-

subroutine create_girder (lat, ix_girder, contrl, girder_info, err_flag)

use bmad_interface, except_dummy => create_girder

implicit none

type (lat_struct), target :: lat
type (ele_struct) :: girder_info
type (control_struct)  contrl(:)
type (ele_struct), pointer ::  slave, slave0, girder_ele
type (ele_pointer_struct), allocatable :: eles(:)

integer, intent(in) :: ix_girder
integer i, j, ix, ix2, ixc, n_con2, n_loc
integer ixs, idel, n_slave

logical err_flag
character(*), parameter :: r_name = 'create_girder'

! girder_info stuff

girder_ele => lat%ele(ix_girder)

girder_ele%name           = girder_info%name
girder_ele%alias          = girder_info%alias
girder_ele%value          = girder_info%value
girder_ele%component_name = girder_info%component_name

if (associated(girder_info%descrip)) then
  allocate (girder_ele%descrip)
  girder_ele%descrip = girder_info%descrip
endif

! Check that origin element exists and is unique.

err_flag = .true.

if (girder_ele%component_name /= '' .and. girder_ele%component_name /= 'GLOBAL_COORDINATES') then
  call lat_ele_locator (girder_ele%component_name, lat, eles, n_loc)

  if (n_loc == 0) then
    call out_io (s_fatal$, r_name, &
                 'GIRDER ORIGIN ELEMENT DOES NOT EXIST: ' // girder_ele%component_name, &
                 'FOR GIRDER: ' // girder_ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  if (n_loc > 1) then
    call out_io (s_fatal$, r_name, &
                 'GIRDER ORIGIN ELEMENT IS NOT UNIQUE: ' // girder_ele%component_name, &
                 'FOR GIRDER: ' // girder_ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif
endif

! Mark element as an girder lord

n_slave = size (contrl)
ix = lat%n_control_max
n_con2 = ix + n_slave

if (n_con2 > size(lat%control)) call reallocate_control (lat, n_con2+500)

do j = 1, n_slave
  lat%control(ix+j)           = contrl(j)
  lat%control(ix+j)%lord      = lat_ele_loc_struct(ix_girder, 0)
  lat%control(ix+j)%ix_attrib = 0
enddo

girder_ele%n_slave = n_slave
girder_ele%ix1_slave = ix + 1
girder_ele%lord_status = girder_lord$
girder_ele%key = girder$
lat%n_control_max = n_con2

! Loop over all slaves
! Free elements convert to overlay slaves.

do i = 1, girder_ele%n_slave

  slave => pointer_to_slave(girder_ele, i)
  if (slave%slave_status == free$) slave%slave_status = minor_slave$

  ! You cannot control super_slaves, group or overlay elements

  if (slave%slave_status == super_slave$ .or. slave%key == group$ .or. &
                                    slave%key == overlay$ .or. slave%key == ramper$) then
    call out_io (s_fatal$, r_name, &
                 'GIRDER IS NOT ALLOWED TO CONTROL SUPER_SLAVE, GROUP OR OVERLAY ELEMENTS.', &
                 'BUT GIRDER: ' // girder_ele%name, 'IS TRYING TO CONTROL: ' // slave%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  ! update controller info for the slave ele

  call add_lattice_control_structs (slave, n_add_lord = 1)
  lat%ic(slave%ic1_lord+slave%n_lord-1) = girder_ele%ix1_slave + i - 1
enddo

!

call find_element_ends (girder_ele, slave0, slave)
girder_ele%value(l$) = slave%s - slave0%s
if (girder_ele%value(l$) < 0) girder_ele%value(l$) = girder_ele%value(l$) + slave0%branch%param%total_length

err_flag = .false.

end subroutine
