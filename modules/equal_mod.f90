module equal_mod

use bmad_utils_mod

interface assignment (=)
  module procedure ele_equal_ele
  module procedure ele_vec_equal_ele_vec
  module procedure lat_equal_lat 
  module procedure lat_vec_equal_lat_vec 
  module procedure branch_equal_branch
  module procedure wall3d_equal_wall3d
end interface

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine ele_equal_ele (ele1, ele2)
!
! Subroutine that is used to set one element equal to another. 
! This routine takes care of the pointers in ele1. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		ele1 = ele2
!
! Input:
!   ele2 -- Ele_struct: Input element.
!
! Output:
!   ele1 -- Ele_struct: Output element.
!-

subroutine ele_equal_ele (ele1, ele2)

use tpsalie_analysis 
use multipole_mod

implicit none
	
type (ele_struct), intent(inout) :: ele1
type (ele_struct), intent(in) :: ele2
type (ele_struct) ele_save

integer i, n2

! 1) Save ele1 pointers in ele_save
! 2) Set ele1 = ele2.

call transfer_ele (ele1, ele_save)
call transfer_ele (ele2, ele1)

! ele1%ix_ele and ele1%ix_branch should not change.
! ele1%lat should not change if ele1 is a component of a lat_struct.
!   Otherwise ele1%lat should point to ele2%lat (For cases where ele1 
!   represents a sliced piece of ele2)

ele1%ix_ele    = ele_save%ix_ele    ! This should not change.
ele1%ix_branch = ele_save%ix_branch ! This should not change.
if (ele1%ix_ele > -1) then          ! If part of a lattice...
  ele1%lat      => ele_save%lat     !   then ele1%lat should not change.
endif
  
! Transfer pointer info.
! When finished ele1's pointers will be pointing to a different memory
! location from ele2's so that the elements are separate.
! Exceptions: %em_field%mode%map, %em_field%mode%grid and %wig.

! %wig%term. 
! The exception for having ele1%wig and ele2%wig point to the same memory location 
! is if we have a periodic_wiggler and ele1 is not a slave or ele2. 
! In this case, the wiggler field depends upon the setting of
! ele%value(b_max$) and ele%value(l_pole$) so sharing the same memeory location would
! lead to trouble if these attributes are modified in one element but not the other.

if (ele1%key == wiggler$ .and. ele1%sub_key == periodic_type$ .and. &
    ele_save%slave_status /= super_slave$ .and. ele_save%slave_status /= multipass_slave$ .and. &
    ele_save%slave_status /= slice_slave$) then

  nullify(ele1%wig)

  if (associated(ele_save%wig) .and. .not. associated(ele2%wig)) then
    ele_save%wig%n_link = ele_save%wig%n_link - 1
    if (ele_save%wig%n_link == 0) deallocate (ele_save%wig)

  elseif (.not. associated(ele_save%wig) .and. associated(ele2%wig)) then
    allocate(ele1%wig)
    n2 = size(ele2%wig%term)
    allocate(ele1%wig%term(n2))
    ele1%wig = ele2%wig

  elseif (associated(ele_save%wig) .and. associated(ele2%wig)) then
    n2 = size(ele2%wig%term)
    if (associated(ele_save%wig, ele2%wig)) then
      ele_save%wig%n_link = ele_save%wig%n_link - 1
      allocate(ele1%wig)
      allocate(ele1%wig%term(n2))
    elseif (size(ele_save%wig%term) /= n2) then
      ele_save%wig%n_link = ele_save%wig%n_link - 1
      if (ele_save%wig%n_link == 0) deallocate (ele_save%wig)
      allocate(ele1%wig)
      allocate(ele1%wig%term(n2))
    else
      ele1%wig => ele_save%wig
    endif
    ele1%wig = ele2%wig
  endif

else
  ele1%wig => ele_save%wig ! Reinstate for trnasfer call 
  call transfer_wig (ele2%wig, ele1%wig)
endif

! If the memory allocated for the wiggler field for ele1 and ele2 are different
! then must adjust the number of links and deallocate if necessary.

! %em_field

ele1%em_field => ele_save%em_field ! Reinstate for transfer call 
call transfer_em_field (ele2%em_field, ele1%em_field)

! %const

if (associated(ele2%const)) then
  if (associated (ele_save%const)) then
    if (size(ele_save%const) == size(ele2%const)) then
      ele1%const => ele_save%const
    else
      deallocate (ele_save%const)
      allocate (ele1%const(size(ele2%const)))
    endif
  else
    allocate (ele1%const(size(ele2%const)))
  endif
  ele1%const = ele2%const
else
  if (associated (ele_save%const)) deallocate (ele_save%const)
endif

! %r

if (associated(ele2%r)) then
  if (associated (ele_save%r)) then
    if (all(lbound(ele_save%r) == lbound(ele2%r)) .and. &
        all(ubound(ele_save%r) == ubound(ele2%r)) ) then
      ele1%r => ele_save%r
    else
      deallocate (ele_save%r)
      allocate (ele1%r(lbound(ele2%r,1):ubound(ele2%r,1), &
                       lbound(ele2%r,2):ubound(ele2%r,2), &
                       lbound(ele2%r,3):ubound(ele2%r,3)))
    endif
  else
    allocate (ele1%r(lbound(ele2%r,1):ubound(ele2%r,1), &
                     lbound(ele2%r,2):ubound(ele2%r,2), &
                     lbound(ele2%r,3):ubound(ele2%r,3)))
  endif
  ele1%r = ele2%r
else
  if (associated (ele_save%r)) deallocate (ele_save%r)
endif

! %taylor

do i = 1, 6
  ele1%taylor(i)%term => ele_save%taylor(i)%term ! reinstate
  ele1%taylor(i) = ele2%taylor(i)      ! use overloaded taylor_equal_taylor
enddo

! %wall3d

ele1%wall3d%section => ele_save%wall3d%section  ! reinstate
ele1%wall3d = ele2%wall3d                       ! use overloaded wall3d_equal_wall3d

! %a_pole, and %b_pole

if (associated(ele2%a_pole)) then
  ele1%a_pole => ele_save%a_pole   ! reinstate
  ele1%b_pole => ele_save%b_pole   ! reinstate
  call multipole_init (ele1)
  ele1%a_pole = ele2%a_pole
  ele1%b_pole = ele2%b_pole
else
  if (associated (ele_save%a_pole)) deallocate (ele_save%a_pole, ele_save%b_pole)
endif

! %descrip

if (associated(ele2%descrip)) then
  if (associated (ele_save%descrip)) then
    ele1%descrip => ele_save%descrip
  else
    allocate (ele1%descrip)
  endif
  ele1%descrip = ele2%descrip
else
  if (associated (ele_save%descrip)) deallocate (ele_save%descrip)
endif

! %mode3

if (associated(ele2%mode3)) then
  if (associated (ele_save%mode3)) then
    ele1%mode3 => ele_save%mode3
  else
    allocate (ele1%mode3)
  endif
  ele1%mode3 = ele2%mode3
else
  if (associated (ele_save%mode3)) deallocate (ele_save%mode3)
endif

! %space_charge

if (associated(ele2%space_charge)) then
  if (associated (ele_save%space_charge)) then
    ele1%space_charge => ele_save%space_charge
  else
    allocate (ele1%space_charge)
  endif
  ele1%space_charge = ele2%space_charge
else
  if (associated (ele_save%space_charge)) deallocate (ele_save%space_charge)
endif

! %rf_wake

ele1%rf_wake => ele_save%rf_wake  ! reinstate
call transfer_rf_wake (ele2%rf_wake, ele1%rf_wake)

! %ptc_genfields are hard because it involves pointers in PTC.
! just kill the ptc_genfield in ele1 for now.

if (associated(ele_save%ptc_genfield)) call kill_ptc_genfield (ele_save%ptc_genfield)
if (associated(ele1%ptc_genfield)) nullify (ele1%ptc_genfield)

end subroutine ele_equal_ele

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine ele_vec_equal_ele_vec (ele1, ele2)
!
! Subroutine that is used to set one element vector equal to another.
! This routine takes care of the pointers in ele1.
!
! Note: This subroutine is called by the overloaded equal sign:
!               ele1(:) = ele2(:)
!
! Input:
!   ele2(:) -- Ele_struct: Input element vector.
!
! Output:
!   ele1(:) -- Ele_struct: Output element vector.
!-

subroutine ele_vec_equal_ele_vec (ele1, ele2)

implicit none

type (ele_struct), intent(inout) :: ele1(:)
type (ele_struct), intent(in) :: ele2(:)

integer i

! error check

if (size(ele1) /= size(ele2)) then
  print *, 'ERROR IN ELE_VEC_EQUAL_ELE_VEC: ARRAY SIZES ARE NOT THE SAME!'
  call err_exit
endif

! transfer

do i = 1, size(ele1)
  call ele_equal_ele (ele1(i), ele2(i))
enddo

end subroutine ele_vec_equal_ele_vec

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
type (lat_struct), intent(in) :: lat_in

integer i, n, ie, n_out, n_in

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
  lat_out%branch(i) = lat_in%branch(i)
  do ie = 0, ubound(lat_out%branch(i)%ele, 1)
    lat_out%branch(i)%ele(ie)%ix_ele = ie
    lat_out%branch(i)%ele(ie)%ix_branch = i
    lat_out%branch(i)%ele(ie)%lat => lat_out
  enddo
enddo

lat_out%ele_init = lat_in%ele_init
lat_out%ele_init%lat => lat_out

! Handle lat%control array

if (allocated (lat_in%control)) then
  n = size(lat_in%control)
  if (.not. allocated(lat_out%control)) allocate(lat_out%control(n))
  if (size(lat_in%control) /= size(lat_out%control)) then
    deallocate (lat_out%control)
    allocate (lat_out%control(n))
  endif
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

lat_out%wall3d = lat_in%wall3d

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
  call err_exit
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

!

branch1%name           = branch2%name
branch1%ix_branch      = branch2%ix_branch
branch1%ix_from_branch = branch2%ix_from_branch
branch1%ix_from_ele    = branch2%ix_from_ele
branch1%n_ele_track    = branch2%n_ele_track
branch1%n_ele_max      = branch2%n_ele_max
call allocate_element_array (branch1%ele, ubound(branch2%ele, 1))
branch1%ele            = branch2%ele
branch1%param          = branch2%param
branch1%a              = branch2%a
branch1%b              = branch2%b
branch1%z              = branch2%z
branch1%ele%ix_branch  = branch2%ix_branch
branch1%wall3d         = branch2%wall3d   

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

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine wall3d_equal_wall3d (wall3d_out, wall3d_in)
!
! Subroutine that is used to set one wall3d equal to another. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		wall3d_out = wall3d_in
!
! Input:
!   wall3d_in -- wall3d_struct: Input wall3d.
!
! Output:
!   wall3d_out -- wall3d_struct: Output wall3d.
!-

elemental subroutine wall3d_equal_wall3d (wall3d_out, wall3d_in)

implicit none
	
type (wall3d_struct), intent(inout) :: wall3d_out
type (wall3d_struct), intent(in) :: wall3d_in

integer i, n_sec, nv

!

if (associated(wall3d_in%section)) then
  n_sec = size(wall3d_in%section)
  if (associated(wall3d_out%section)) then
    if (size(wall3d_out%section) /= n_sec) deallocate (wall3d_out%section)
  endif
  if (.not. associated(wall3d_out%section)) allocate(wall3d_out%section(n_sec))

  do i = 1, n_sec
    if (allocated(wall3d_in%section(i)%v)) then
      nv = size(wall3d_in%section(i)%v)
      if (allocated(wall3d_out%section(i)%v)) then
        if (size(wall3d_out%section(i)%v) /= nv) deallocate(wall3d_out%section(i)%v)
      endif
      if (.not. allocated(wall3d_out%section(i)%v)) allocate(wall3d_out%section(i)%v(nv))
    else
      if (allocated(wall3d_out%section(i)%v)) deallocate(wall3d_out%section(i)%v)
    endif
    wall3d_out%section(i) = wall3d_in%section(i)
  enddo 

else
  if (associated(wall3d_out%section)) deallocate(wall3d_out%section)
endif
 
end subroutine wall3d_equal_wall3d

end module

