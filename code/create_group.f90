!+
! Subroutine create_group (ring, ix_ele, contrl)
!
! Subroutine to create a group control element.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring          -- Ring_struct: Lattice
!     %ele_(ix_ele)%value(command$) -- Real(rp): Initial command value.
!   ix_ele        -- Integer: Index of group element (see below)
!   contrl(:)     -- Control_struct: What to control.
!     %ix_slave       -- Integer: Index in RING%ELE_() of element controlled.
!     %ix_attrib      -- Integer: Index in %VALUE() array of
!                                 attribute controlled.
!     %coef           -- Real(rp): Coefficient.
!
! Output:
!     ring -- Ring_struct: Appropriate values are set in the RING structure.
!
! Note: Use NEW_CONTROL to get an index for the group element
!
! Example:
!   call new_control (ring, ix_ele)   ! get IX_ELE index
!   ring%ele_(ix_ele)%name = 'GROUP1' ! group name
!   ring%ele_(ix_ele)%value(command$) = 0 ! start at zero
!   n_control = 2               ! control 2 elements
!
!   contrl(1)%ix_slave = 10   ! RING%ELE_(10) is Q01W say.
!   contrl(1)%ix_attrib = k1$ ! The group controls the quadrupole strength.
!   contrl(1)%coef = 0.1      ! A change in the group value of 1 produces
!                             !    a change of 0.1 in k1 of element 10.
!
!   contrl(2)%ix_slave = 790  ! RING%ELE_(790) is Q01E say.
!   contrl(2)%ix_attrib = k1$ ! The group controls the quadrupole strength.
!   contrl(2)%coef = -0.1     ! make changes antisymmetric.
!
!   call create_group (ring, ix_ele, 2, contrl)  ! create the group
!
! Notes:
!   A) The value of the group is stored in RING%ELE_(IX_ELE)%VALUE(COMMAND$)
!   B) Only changes from the previous value are significant. The
!      old value is stored in RING%ELE_(IX_ELE)%VALUE(OLD_COMMAND$).
!   C) Use CONTROL_BOOKKEEPER to update the attributes of the elements
!      controlled by the group element.
!   D) Use RING_MAKE_MAT6 to update the attributes AND remake MAT6 for the
!     elements controlled by the group element.
!
! Proceeding with the previous example:
!   ix_ele1 = ix_ele                        ! save the group index
!   ring%ele_(ix_ele1)%value(command$) = 10 ! put in a value for the group
!   call ring_make_mat6 (ring, ix_ele1)     ! update the k1's for the 2 quads
!                                           !   AND remake the MAT6's
!   call control_bookkeeper (ring, ix_ele1) ! use this instead
!                                           !   to only update the k1's
!   ring%ele_(ix_ele1)%value(command$) = 13 ! put in a new value for the group
!   call ring_make_mat6 (ring, ix_ele1)     ! the change now in k1's is
!                                           !   based upon the delta: 13 - 10
!
! Note: You can control an element's position by setting:
!       CONTRL(i)%IX_ATTRIB = START_EDGE$      or
!                           = END_EDGE$        or
!                           = ACCORDION_EDGE$  or
!                           = SYMMETRIC_EDGE$
!
! %IX_ATTRIB = START_EDGE$ and %IX_ATTRIB = END_EDGE$ controls the
! placement of the edges of an element keeping the ring total length invariant.
! this is done by lengthening and shortening the elements to either side
! keeping the total ring length invariant.
!
! %IX_ATTRIB = ACCORDION_EDGE$ and %IX_ATTRIB = SYMMETRIC_EDGE$ moves both
! the start and end edges simultaneously.
! ACCORDION_EDGE$ moves the edges antisymmetrically (and thus there is a length
! change of the element).
! SYMMETRIC_EDGE$ moves the edges symmetrically creating a z-offset with no
! length change.
!-

#include "CESR_platform.inc"

subroutine create_group (ring, ix_ele, contrl)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct)  ring
  type (control_struct)  contrl(:)

  integer i, ix_ele, ixa, n_control, n_con
  integer ix1, ix2, ix_min, ix_max, ixe

! init

  n_control = size(contrl)
  ring%ele_(ix_ele)%control_type = group_lord$
  ring%ele_(ix_ele)%key = group$
  n_con = ring%n_control_max
  ring%ele_(ix_ele)%ix1_slave = n_con + 1

! loop over all controlled elements

  do i = 1, n_control

! For position control: We need to figure out the elements that
! need to be controlled.
! Find beginning and ending positions of element
! if a super_lord then we must go to the slave elements to find the ends
! else not a super lord so finding the ends is simple

    ixe = contrl(i)%ix_slave
    ixa = contrl(i)%ix_attrib

    if (ixa == start_edge$ .or. ixa == end_edge$ .or. &
                                              ixa == accordion_edge$) then

      if (ring%ele_(ixe)%control_type == super_lord$) then
        ix_min = ring%control_(ring%ele_(ixe)%ix1_slave)%ix_slave
        ix_max = ring%control_(ring%ele_(ixe)%ix2_slave)%ix_slave
      elseif (ixe < ring%n_ele_use) then
        ix_min = ixe
        ix_max = ixe
      else
        print *, 'ERROR IN CREATE_GROUP: A GROUP IS NOT ALLOWED TO CONTROL'
        print *, '      A ', control_name(ring%ele_(ixe)%control_type)
        print *, '      YOU TRIED TO CONTROL: ', ring%ele_(ixe)%name
        call err_exit
      endif

! now that we have the ends we find the elements to either side whose length
! the group can adjust

      ix1 = ix_min - 1
      do 
        if (ring%ele_(ix1)%value(l$) == 0) then
          ix1 = ix1 - 1
        else
          exit
        endif
      enddo

      ix2 = ix_max + 1 
      do
        if (ring%ele_(ix2)%value(l$) == 0) then
          ix2 = ix2 + 1
        else
          exit
        endif
      enddo

      if (ixa == start_edge$ .or. ixa == accordion_edge$ .or. &
                                       ixa == symmetric_edge$) then
        if (ix1 < 1) then
          print *, 'ERROR IN CREATE_GROUP: START_EDGE OF CONTROLED'
          print *, '      ELEMENT IS AT BEGINNING OF RING AND CANNOT BE'
          print *, '      VARIED FOR GROUP: ', ring%ele_(ix_ele)%name
          call err_exit
        endif
      endif

      if (ixa == end_edge$ .or. ixa == accordion_edge$ .or. &
                                        ixa == symmetric_edge$) then
        if (ix2 > ring%n_ele_use) then
          print *, 'ERROR IN CREATE_GROUP: END_EDGE OF CONTROLED'
          print *, '      ELEMENT IS AT END OF RING AND CANNOT BE'
          print *, '      VARIED FOR GROUP: ', ring%ele_(ix_ele)%name
          call err_exit
        endif
      endif

! put in coefficients

      select case (ixa)

      case (start_edge$)
        call bookit (ix1, 1)
        call bookit (ix_min, -1)

      case (end_edge$)
        call bookit (ix_max, 1)
        call bookit (ix2, -1)

      case (accordion_edge$)
        call bookit (ix1, -1)
        if (ix_min == ix_max) then
          call bookit (ix_min, 2)
        else
          call bookit (ix_min, 1)
          call bookit (ix_max, 1)
        endif
        call bookit (ix2, -1)

      case (symmetric_edge$)
        call bookit (ix1, 1)
        call bookit (ix2, -1)

      end select

! for all else without position control the group setup is simple.

    else

      n_con = n_con + 1
      if (n_con > size(ring%control_)) &
                      ring%control_ => reallocate (ring%control_, n_con+500)
      ring%control_(n_con) = contrl(i)
      ring%control_(n_con)%ix_lord = ix_ele

    endif

  enddo

! final bookkeping

  ring%ele_(ix_ele)%ix2_slave = n_con
  ring%ele_(ix_ele)%n_slave = n_con - ring%ele_(ix_ele)%ix1_slave + 1
  ring%n_control_max = n_con

!---------------------------------------------------------------------------

contains

subroutine bookit (i_ele, scale)

  integer scale, i_ele

  n_con = n_con + 1
  if (n_con > size(ring%control_)) &
                      ring%control_ => reallocate (ring%control_, n_con+500)
  ring%control_(n_con)%ix_lord = ix_ele
  ring%control_(n_con)%ix_slave = i_ele
  ring%control_(n_con)%ix_attrib = l$
  ring%control_(n_con)%coef = scale * contrl(i)%coef

end subroutine

end subroutine
