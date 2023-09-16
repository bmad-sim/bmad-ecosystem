!+
! Function pointer_to_next_ele (this_ele, offset, skip_beginning, follow_fork) result (next_ele)
!
! Function to return a pointer to the next (if offset = 1) tracking element relative to this_ele.
! 
! If the this_ele is a super_lord element, the appropriate element in the tracking 
! part of the lattice is returned.
!
! If this_ele is a lord element but not a super_lord then this is an error since
! it is not clear what to do in this case.
!
! This routine will always wrap around between branch end and branch beginning
! independent of the branch geometry. That is, it is the responsibility of the
! calling routine to check this if important.
!
! Notice that the first element in a lattice is the beginning element with index 0.
!
! Also see:
!   iterator_to_ele_in_branch
!
! Input:
!   this_ele       -- ele_struct: Starting element.
!   offset         -- integer, optional: +1 -> return next element, +2 -> element 
!                       after that, etc. Can be negative. Default = +1.
!   skip_beginning -- logical, optional: If True then skip beginning element #0
!                       when wrapping around. Default is False.
!   follow_fork    -- logical, optional: If True then fork at any fork element.
!                       Default is False.
!
!   next_ele -- ele_struct, pointer: Element after this_ele (if offset = 1).
!                Nullified if there is an error. EG bad this_ele.
!-

function pointer_to_next_ele (this_ele, offset, skip_beginning, follow_fork) result (next_ele)

use equal_mod, dummy => pointer_to_next_ele

implicit none

type (ele_struct), target :: this_ele
type (ele_struct), pointer :: next_ele
type (ele_struct), pointer :: an_ele
type (branch_struct), pointer :: branch

integer, optional :: offset
integer i, ix_ele, n_off
logical, optional :: skip_beginning, follow_fork

! If element is not associated with a lat_struct then do nothing.

next_ele => null()

if (.not. associated(this_ele%branch)) return   ! Not associated with a lat_struct
if (this_ele%ix_ele < 0 .or. this_ele%ix_ele > this_ele%branch%n_ele_max) return  ! Slice_slave?

! Compute offset

n_off = integer_option(+1, offset)
if (n_off == 0) then
  next_ele => this_ele
  return
endif

! If a lord then point to a slave. 
! Initially point to the first or last slave element so that this routine
! will return a pointer to an element that is not a slave of this_ele.

if (this_ele%ix_ele > this_ele%branch%n_ele_track) then  ! Is a lord
  if (this_ele%lord_status /= super_lord$) return    ! Error
  if (n_off > 0) then
    an_ele => pointer_to_slave(this_ele, this_ele%n_slave)
  else
    an_ele => pointer_to_slave(this_ele, 1)
  endif
else
  an_ele => this_ele
endif

! Apply offset.
! If follow_fork = True then must check all elements in between for a possible fork element.

ix_ele = an_ele%ix_ele
branch => an_ele%branch

if (n_off > 0) then

  do i = 1, abs(n_off)

    if (logic_option(.false., follow_fork) .and. (an_ele%key == fork$ .or. an_ele%key == photon_fork$)) then
      ix_ele = nint(an_ele%value(ix_to_element$))
      branch => an_ele%branch%lat%branch(nint(an_ele%value(ix_to_branch$)))
    else
      ix_ele = ix_ele + 1
    endif

    if (ix_ele == branch%n_ele_track + 1) then
      if (logic_option(.false., skip_beginning)) then
        ix_ele = 1
      else
        ix_ele = 0
      endif
    endif

    an_ele => branch%ele(ix_ele)

  enddo

! Going backwards

else
  do i = 1, abs(n_off)

    if (logic_option(.false., follow_fork) .and. nint(an_ele%value(direction$)) == -1 .and. & 
                            (an_ele%key == fork$ .or. an_ele%key == photon_fork$)) then
      ix_ele = nint(an_ele%value(ix_to_element$))
      branch => an_ele%branch%lat%branch(nint(an_ele%value(ix_to_branch$)))
    else
      ix_ele = ix_ele - 1
    endif

    if (logic_option(.false., skip_beginning)) then
      if (ix_ele == 0) ix_ele = branch%n_ele_track
    else
      if (ix_ele == -1) ix_ele = branch%n_ele_track
    endif

    an_ele => branch%ele(ix_ele)

  enddo
endif

next_ele => an_ele

end function pointer_to_next_ele
