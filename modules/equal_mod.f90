module equal_mod

use bmad_core_struct_mod

interface assignment (=)
  module procedure ele_equal_ele
  module procedure lat_equal_lat 
  module procedure lat_vec_equal_lat_vec 
  module procedure branch_equal_branch
end interface

interface operator (+)
  module procedure em_field_plus_em_field
end interface

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Function em_field_plus_em_field (field1, field2) result (field_tot)
!
! Subroutine to add fields.
!
! Note: This subroutine is called by the overloaded plus sign:
!		field_tot = field1 + field2 
!
! Input:
!   field1 -- em_field_struct: Input field
!   field2 -- em_field_struct: Input field
!
! Output:
!   field_tot -- em_field_struct: Combined field.
!-

function em_field_plus_em_field (field1, field2) result (field_tot)

type (em_field_struct), intent(in) :: field1, field2
type (em_field_struct) field_tot

!

field_tot%e = field1%e + field2%e
field_tot%b = field1%b + field2%b

field_tot%de = field1%de + field2%de
field_tot%db = field1%db + field2%db

end function em_field_plus_em_field 

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine ele_equal_ele (ele_out, ele_in)
!
! Subroutine that is used to set one element equal to another. 
! This routine takes care of the pointers in ele_out. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		ele_out = ele_in
!
! Input:
!   ele_in -- Ele_struct: Input element.
!
! Output:
!   ele_out -- Ele_struct: Output element.
!-

subroutine ele_equal_ele (ele_out, ele_in)

implicit none
	
type (ele_struct), intent(inout) :: ele_out
type (ele_struct), intent(in) :: ele_in
type (ele_struct) ele_save

integer i, j, n, n2, ub(2), ub1

! 1) Save ele_out pointers in ele_save
! 2) Set ele_out = ele_in.

call transfer_ele (ele_out, ele_save)
call transfer_ele (ele_in, ele_out)

! ele_out%ix_ele and ele_out%ix_branch should not change.
! ele_out%branch should not change if ele_out is a component of a lat_struct.
!   Otherwise ele_out%lat should point to ele_in%lat (For cases where ele_out 
!   represents a sliced piece of ele_in)

ele_out%ix_ele    = ele_save%ix_ele    ! This should not change.
ele_out%ix_branch = ele_save%ix_branch ! This should not change.
if (ele_out%ix_ele > -1) then          ! If part of a lattice...
  ele_out%branch => ele_save%branch     !   then ele_out%branch should not change.
endif

! Transfer pointer info.
! When finished ele_out's pointers will be pointing to a different memory
! location from ele_in's so that the elements are separate.
! Exceptions: %em_field%mode%cylindrical_map, %em_field%mode%grid.

! %cartesian_map exception: The problem with having ele_out%cartesian_map and ele_in%cartesian_map point
! to the same memory location is when we have a periodic_wiggler and ele_out is not a slave of ele_in. 
! In this case, the wiggler field depends upon the setting of
! ele%value(b_max$) and ele%value(l_pole$) so sharing the same memeory location would
! lead to trouble if these attributes are modified in one element but not the other.

! If the memory allocated for the wiggler field for ele_out and ele_in are different
! then must adjust the number of links and deallocate if necessary.

! Note: A periodic_type wiggler always has one cartesian_map createded in attribute_bookkeeper.

if ((ele_out%key == wiggler$ .or. ele_out%key == undulator$) .and. ele_out%sub_key == periodic_type$ .and. &
    ele_save%slave_status /= super_slave$ .and. ele_save%slave_status /= multipass_slave$ .and. &
    ele_save%slave_status /= slice_slave$) then

  if (associated(ele_save%cartesian_map)) then
    if (size(ele_save%cartesian_map) /= 1 .or. .not. associated(ele_in%cartesian_map)) then
      call unlink_fieldmap (cartesian_map = ele_save%cartesian_map)
    endif
  endif

  nullify(ele_out%cartesian_map)

  if (associated(ele_in%cartesian_map)) then
    n2 = size(ele_in%cartesian_map(1)%ptr%term)  ! Should be 1

    if (associated(ele_save%cartesian_map)) then
      ele_out%cartesian_map => ele_save%cartesian_map
      if (associated(ele_out%cartesian_map(1)%ptr, ele_in%cartesian_map(1)%ptr)) &
                            ele_out%cartesian_map(1)%ptr%n_link = ele_out%cartesian_map(1)%ptr%n_link - 1
    else
      allocate(ele_out%cartesian_map(1))
    endif

    ele_out%cartesian_map = ele_in%cartesian_map
    allocate(ele_out%cartesian_map(1)%ptr)
    allocate(ele_out%cartesian_map(1)%ptr%term(n2))
    ele_out%cartesian_map(1)%ptr%term = ele_in%cartesian_map(1)%ptr%term
    write (ele_out%cartesian_map(1)%ptr%file, '(a, i0, a, i0)') 'ele_equal_ele:', &
                                                       ele_out%ix_branch, '>>', ele_out%ix_ele ! Unique name
  endif

else
  ele_out%cartesian_map => ele_save%cartesian_map ! Reinstate for transfer call 
  call transfer_fieldmap (ele_in, ele_out, cartesian_map$)
endif

! %cylindrical_map, etc.

ele_out%cylindrical_map => ele_save%cylindrical_map ! Reinstate for transfer call 
call transfer_fieldmap (ele_in, ele_out, cylindrical_map$)

ele_out%grid_field => ele_save%grid_field ! Reinstate for transfer call 
call transfer_fieldmap (ele_in, ele_out, grid_field$)

ele_out%taylor_field => ele_save%taylor_field ! Reinstate for transfer call
call transfer_fieldmap (ele_in, ele_out, taylor_field$)

! %rad_int_cache

if (associated(ele_in%rad_int_cache)) then
  if (associated (ele_save%rad_int_cache)) then
      ele_out%rad_int_cache => ele_save%rad_int_cache
  else
    allocate (ele_out%rad_int_cache)
  endif
  ele_out%rad_int_cache = ele_in%rad_int_cache
else
  if (associated (ele_save%rad_int_cache)) deallocate (ele_save%rad_int_cache)
endif

! %r

if (associated(ele_in%r)) then
  if (associated (ele_save%r)) then
    if (all(lbound(ele_save%r) == lbound(ele_in%r)) .and. &
        all(ubound(ele_save%r) == ubound(ele_in%r)) ) then
      ele_out%r => ele_save%r
    else
      deallocate (ele_save%r)
      allocate (ele_out%r(lbound(ele_in%r,1):ubound(ele_in%r,1), &
                       lbound(ele_in%r,2):ubound(ele_in%r,2), &
                       lbound(ele_in%r,3):ubound(ele_in%r,3)))
    endif
  else
    allocate (ele_out%r(lbound(ele_in%r,1):ubound(ele_in%r,1), &
                     lbound(ele_in%r,2):ubound(ele_in%r,2), &
                     lbound(ele_in%r,3):ubound(ele_in%r,3)))
  endif
  ele_out%r = ele_in%r
else
  if (associated (ele_save%r)) deallocate (ele_save%r)
endif

! %photon

if (associated(ele_in%photon)) then
  ele_out%photon => ele_save%photon  ! reinstate
  if (.not. associated(ele_out%photon)) allocate(ele_out%photon)

  if (allocated (ele_in%photon%surface%grid%pt)) then
    ub = ubound(ele_in%photon%surface%grid%pt)
    if (allocated (ele_out%photon%surface%grid%pt)) then
      if (any(ub /= ubound(ele_out%photon%surface%grid%pt))) deallocate (ele_out%photon%surface%grid%pt)
    endif
    if (.not. allocated (ele_out%photon%surface%grid%pt)) allocate (ele_out%photon%surface%grid%pt(0:ub(1), 0:ub(2)))
  else
    if (allocated(ele_out%photon%surface%grid%pt)) deallocate (ele_out%photon%surface%grid%pt)
  endif

  ele_out%photon = ele_in%photon
else
  if (associated (ele_save%photon)) deallocate (ele_save%photon)
endif

! %control_var

if (associated(ele_in%control_var)) then
  n = size(ele_in%control_var)
  ele_out%control_var => ele_save%control_var   ! reinstate
  if (associated(ele_out%control_var)) then
    if (size(ele_out%control_var) /= n) deallocate(ele_out%control_var)
  endif
  if (.not. associated(ele_out%control_var)) allocate(ele_out%control_var(n))
  ele_out%control_var = ele_in%control_var

else
  if (associated (ele_save%control_var)) deallocate(ele_save%control_var)
endif

! %taylor

do i = 1, 6
  ele_out%taylor(i)%term => ele_save%taylor(i)%term ! reinstate
  ele_out%taylor(i) = ele_in%taylor(i)      ! use overloaded taylor_equal_taylor
enddo

! %spin_taylor

do i = 1, 3
do j = 1, 3
  ele_out%spin_taylor(i,j)%term => ele_save%spin_taylor(i,j)%term ! reinstate
  ele_out%spin_taylor(i,j) = ele_in%spin_taylor(i,j)      ! use overloaded taylor_equal_taylor
enddo
enddo

! %wall3d

ele_out%wall3d => ele_save%wall3d        ! reinstate
call transfer_wall3d (ele_in%wall3d, ele_out%wall3d)

! %a_pole, and %b_pole

if (associated(ele_in%a_pole)) then
  ele_out%a_pole => ele_save%a_pole   ! reinstate
  ele_out%b_pole => ele_save%b_pole   ! reinstate
  call multipole_init (ele_out, magnetic$)
  ele_out%a_pole = ele_in%a_pole
  ele_out%b_pole = ele_in%b_pole
else
  if (associated (ele_save%a_pole)) deallocate (ele_save%a_pole, ele_save%b_pole)
endif

! %a_pole_elec, and %b_pole_elec

if (associated(ele_in%a_pole_elec)) then
  ele_out%a_pole_elec => ele_save%a_pole_elec   ! reinstate
  ele_out%b_pole_elec => ele_save%b_pole_elec   ! reinstate
  call multipole_init (ele_out, electric$)
  ele_out%a_pole_elec = ele_in%a_pole_elec
  ele_out%b_pole_elec = ele_in%b_pole_elec
else
  if (associated (ele_save%a_pole_elec)) deallocate (ele_save%a_pole_elec, ele_save%b_pole_elec)
endif

! %descrip

if (associated(ele_in%descrip)) then
  if (associated (ele_save%descrip)) then
    ele_out%descrip => ele_save%descrip
  else
    allocate (ele_out%descrip)
  endif
  ele_out%descrip = ele_in%descrip
else
  if (associated (ele_save%descrip)) deallocate (ele_save%descrip)
endif

! %mode3

if (associated(ele_in%mode3)) then
  if (associated (ele_save%mode3)) then
    ele_out%mode3 => ele_save%mode3
  else
    allocate (ele_out%mode3)
  endif
  ele_out%mode3 = ele_in%mode3
else
  if (associated (ele_save%mode3)) deallocate (ele_save%mode3)
endif

! %space_charge

if (associated(ele_in%space_charge)) then
  if (associated (ele_save%space_charge)) then
    ele_out%space_charge => ele_save%space_charge
  else
    allocate (ele_out%space_charge)
  endif
  ele_out%space_charge = ele_in%space_charge
else
  if (associated (ele_save%space_charge)) deallocate (ele_save%space_charge)
endif

! %wake

ele_out%wake => ele_save%wake  ! reinstate
call transfer_wake (ele_in%wake, ele_out%wake)

! %ptc_genfield%fields are hard because it involves pointers in PTC.
! just kill the ptc_genfield in ele_out for now.

if (associated(ele_save%ptc_genfield%field)) call kill_ptc_genfield (ele_save%ptc_genfield%field)
if (associated(ele_out%ptc_genfield%field)) nullify (ele_out%ptc_genfield%field)

end subroutine ele_equal_ele

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine lat_equal_lat (lat_out, lat_in)
!
! Subroutine that is used to set one lat equal to another. 
! This routine takes care of the pointers in lat_in. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		lat_out = lat_in
!
! Input:
!   lat_in -- lat_struct: Input lat.
!
! Output:
!   lat_out -- lat_struct: Output lat.
!-

subroutine lat_equal_lat (lat_out, lat_in)

implicit none

type (lat_struct), intent(inout), target :: lat_out
type (lat_struct), intent(in), target :: lat_in
type (branch_struct), pointer :: branch_out
type (control_struct), pointer :: c_in, c_out
integer i, n, ie, n_out, n_in

! Kill allociated PTC layouts if they exist

call kill_ptc_layouts(lat_out)

! If the element arrays have not been initialized in lat_in then deallocate lat_out.

if (.not. allocated (lat_in%branch) .or. .not. associated(lat_in%ele)) then
  call deallocate_lat_pointers (lat_out)
  return
endif

! Care must be taken here since lat%ele points to the same memory as lat%branch(0).
! First take care of the branch lines.

n = ubound(lat_in%branch, 1)
call allocate_branch_array (lat_out, n)

do i = 0, n
  call allocate_lat_ele_array (lat_out, ubound(lat_in%branch(i)%ele, 1), i)
  branch_out => lat_out%branch(i)
  branch_out =  lat_in%branch(i)
  branch_out%lat => lat_out
  do ie = 0, ubound(branch_out%ele, 1)
    branch_out%ele(ie)%ix_ele = ie
    branch_out%ele(ie)%ix_branch = i
    branch_out%ele(ie)%branch => branch_out
  enddo
enddo

lat_out%ele_init = lat_in%ele_init
nullify(lat_out%ele_init%branch)

! Handle lat%control array

if (allocated (lat_in%control)) then
  n = size(lat_in%control)
  if (.not. allocated(lat_out%control)) allocate(lat_out%control(n))
  if (size(lat_in%control) /= size(lat_out%control)) then
    deallocate (lat_out%control)
    allocate (lat_out%control(n))
  endif

  do i = 1, size(lat_in%control)
    c_in => lat_in%control(i); c_out => lat_out%control(i)
    if (allocated(c_in%stack)) then
      n = size(c_in%stack)
      if (.not. allocated(c_out%stack)) allocate(c_out%stack(n))
      if (size(c_out%stack) /= size(c_in%stack)) then
        deallocate (c_out%stack)
        allocate (c_out%stack(n))
      endif

    else
      if (allocated(c_out%stack)) deallocate(c_out%stack)
    endif
  enddo

  lat_out%control = lat_in%control
else
  if (allocated(lat_out%control)) deallocate (lat_out%control)
endif

! handle lat%ic array

if (allocated(lat_in%ic)) then
  call re_allocate(lat_out%ic, size(lat_in%ic))
  lat_out%ic = lat_in%ic
else
  if (allocated(lat_out%ic)) deallocate (lat_out%ic)
endif

! lat%attribute_alias

if (allocated(lat_in%attribute_alias)) then
  n = size(lat_in%attribute_alias)
  if (.not. allocated(lat_out%attribute_alias)) allocate (lat_out%attribute_alias(n))
  if (size(lat_out%attribute_alias) /= n) then
    deallocate(lat_out%attribute_alias)
    allocate (lat_out%attribute_alias(n))
  endif
  lat_out%attribute_alias = lat_in%attribute_alias
else
  if (allocated(lat_out%attribute_alias)) deallocate(lat_out%attribute_alias)
endif

! non-pointer transfer

call transfer_lat_parameters (lat_in, lat_out)

end subroutine lat_equal_lat

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine lat_vec_equal_lat_vec (lat1, lat2)
!
! Subroutine that is used to set one lat vector equal to another. 
! This routine takes care of the pointers in lat1. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		lat1(:) = lat2(:)
!
! Input:
!   lat2(:) -- lat_struct: Input lat vector.
!
! Output:
!   lat1(:) -- lat_struct: Output lat vector.
!-

subroutine lat_vec_equal_lat_vec (lat1, lat2)

implicit none
	
type (lat_struct), intent(inout) :: lat1(:)
type (lat_struct), intent(in) :: lat2(:)

integer i

! error check

if (size(lat1) /= size(lat2)) then
  print *, 'ERROR IN lat_vec_equal_lat_vec: ARRAY SIZES ARE NOT THE SAME!'
  if (global_com%exit_on_error) call err_exit
endif

! transfer

do i = 1, size(lat1)
  call lat_equal_lat (lat1(i), lat2(i))
enddo

end subroutine lat_vec_equal_lat_vec 

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine branch_equal_branch (branch1, branch2)
!
! Subroutine that is used to set one branch equal to another. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		branch1 = branch2
!
! Input:
!   branch2 -- branch_struct: Input branch.
!
! Output:
!   branch1 -- branch_struct: Output branch.
!-

subroutine branch_equal_branch (branch1, branch2)

implicit none
	
type (branch_struct), intent(inout) :: branch1
type (branch_struct), intent(in) :: branch2
integer i

!

branch1%name           = branch2%name
branch1%ix_branch      = branch2%ix_branch
branch1%ix_from_branch = branch2%ix_from_branch
branch1%ix_from_ele    = branch2%ix_from_ele
branch1%n_ele_track    = branch2%n_ele_track
branch1%n_ele_max      = branch2%n_ele_max

call allocate_element_array (branch1%ele, ubound(branch2%ele, 1))
do i = 0, ubound(branch2%ele, 1)
  branch1%ele(i)  = branch2%ele(i)
enddo

branch1%param          = branch2%param
branch1%a              = branch2%a
branch1%b              = branch2%b
branch1%z              = branch2%z
branch1%ele%ix_branch  = branch2%ix_branch
branch1%wall3d         => branch2%wall3d   

end subroutine branch_equal_branch

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine coord_equal_coord (coord1, coord2)
!
! Subroutine that is used to set one coord equal to another. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		coord1 = coord2
!
! Input:
!   coord2 -- coord_struct: Input coord.
!
! Output:
!   coord1 -- coord_struct: Output coord.
!-

elemental subroutine coord_equal_coord (coord1, coord2)

implicit none
	
type (coord_struct), intent(inout) :: coord1
type (coord_struct), intent(in) :: coord2

!

coord1%vec = coord2%vec
coord1%spin = coord2%spin
 
end subroutine coord_equal_coord

end module

