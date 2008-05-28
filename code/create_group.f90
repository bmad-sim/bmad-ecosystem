!+
! Subroutine create_group (lat, ix_lord, contrl, err, err_print_flag)
!
! Subroutine to create a group control element.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat          -- lat_struct: Lattice
!     %ele(ix_lord)%value(command$) -- Real(rp): Initial command value.
!   ix_lord       -- Integer: Index of group lord element (see below).
!   contrl(:)     -- Control_struct: What to control.
!     %ix_slave       -- Integer: Index in LAT%ele() of element controlled.
!     %ix_attrib      -- Integer: Index in %VALUE() array of
!                                 attribute controlled.
!     %coef           -- Real(rp): Coefficient.
!   err            -- Logical: Set True if an attribute is not free to be controlled.
!   err_print_flag -- Logical, optional: If present and False then supress                                
!                       printing of an error message if attribute is not free.  
!
! Output:
!     lat -- lat_struct: Appropriate values are set in the LAT structure.
!
! Note: Use NEW_CONTROL to get an index for the group element
!
! Example:
!   call new_control (lat, ix_lord)   ! get IX_LORD index
!   lat%ele(ix_lord)%name = 'GROUP1' ! group name
!   lat%ele(ix_lord)%value(command$) = 0 ! start at zero
!   n_control = 2               ! control 2 elements
!
!   contrl(1)%ix_slave = 10   ! LAT%ele(10) is Q01W say.
!   contrl(1)%ix_attrib = k1$ ! The group controls the quadrupole strength.
!   contrl(1)%coef = 0.1      ! A change in the group value of 1 produces
!                             !    a change of 0.1 in k1 of element 10.
!
!   contrl(2)%ix_slave = 790  ! LAT%ele(790) is Q01E say.
!   contrl(2)%ix_attrib = k1$ ! The group controls the quadrupole strength.
!   contrl(2)%coef = -0.1     ! make changes antisymmetric.
!
!   call create_group (lat, ix_lord, 2, contrl)  ! create the group
!
! Notes:
!   A) The value of the group is stored in LAT%ele(IX_LORD)%VALUE(COMMAND$)
!   B) Only changes from the previous value are significant. The
!      old value is stored in LAT%ele(IX_LORD)%VALUE(OLD_COMMAND$).
!   C) Use CONTROL_BOOKKEEPER to update the attributes of the elements
!      controlled by the group element.
!   D) Use lat_make_mat6 to update the attributes AND remake MAT6 for the
!     elements controlled by the group element.
!
! Proceeding with the previous example:
!   ix_lord1 = ix_lord                       ! save the group index
!   lat%ele(ix_lord1)%value(command$) = 10   ! put in a value for the group
!   call lat_make_mat6 (lat, ix_lord1)       ! update the k1's for the 2 quads
!                                            !   AND remake the MAT6's
!   call control_bookkeeper (lat, ix_lord1)  ! use this instead
!                                            !   to only update the k1's
!   lat%ele(ix_lord1)%value(command$) = 13   ! put in a new value for the group
!   call lat_make_mat6 (lat, ix_lord1)       ! the change now in k1's is
!                                            !   based upon the delta: 13 - 10
!
! Note: You can control an element's position by setting:
!       CONTRL(i)%IX_ATTRIB = START_EDGE$      or
!                           = END_EDGE$        or
!                           = ACCORDION_EDGE$  or
!                           = SYMMETRIC_EDGE$
!
! %IX_ATTRIB = START_EDGE$ and %IX_ATTRIB = END_EDGE$ controls the
! placement of the edges of an element keeping the lat total length invariant.
! this is done by lengthening and shortening the elements to either side
! keeping the total lat length invariant.
!
! %IX_ATTRIB = ACCORDION_EDGE$ and %IX_ATTRIB = SYMMETRIC_EDGE$ moves both
! the start and end edges simultaneously.
! ACCORDION_EDGE$ moves the edges antisymmetrically (and thus there is a length
! change of the element).
! SYMMETRIC_EDGE$ moves the edges symmetrically creating a z-offset with no
! length change.
!-

#include "CESR_platform.inc"

subroutine create_group (lat, ix_lord, contrl, err, err_print_flag)

  use bmad_struct
  use bmad_interface, except_dummy => create_group

  implicit none

  type (lat_struct)  lat
  type (control_struct)  contrl(:)

  integer i, ix_lord, ixa, n_control, n_con
  integer ix1, ix2, ix_min, ix_max, ixe

  logical err, free
  logical, optional :: err_print_flag

! init

  call check_controller_controls (contrl, lat%ele(ix_lord)%name, err)
  if (err) return

  n_control = size(contrl)
  lat%ele(ix_lord)%control_type = group_lord$
  lat%ele(ix_lord)%key = group$
  n_con = lat%n_control_max
  lat%ele(ix_lord)%ix1_slave = n_con + 1

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

      if (lat%ele(ixe)%control_type == super_lord$) then
        ix_min = lat%control(lat%ele(ixe)%ix1_slave)%ix_slave
        ix_max = lat%control(lat%ele(ixe)%ix2_slave)%ix_slave
      elseif (ixe < lat%n_ele_track) then
        ix_min = ixe
        ix_max = ixe
      else
        print *, 'ERROR IN CREATE_GROUP: A GROUP IS NOT ALLOWED TO CONTROL'
        print *, '      A ', control_name(lat%ele(ixe)%control_type)
        print *, '      YOU TRIED TO CONTROL: ', lat%ele(ixe)%name
        call err_exit
      endif

      ! now that we have the ends we find the elements to either side whose length
      ! the group can adjust

      ix1 = ix_min - 1
      do 
        if (lat%ele(ix1)%value(l$) == 0) then
          ix1 = ix1 - 1
        else
          exit
        endif
      enddo

      ix2 = ix_max + 1 
      do
        if (lat%ele(ix2)%value(l$) == 0) then
          ix2 = ix2 + 1
        else
          exit
        endif
      enddo

      if (ixa == start_edge$ .or. ixa == accordion_edge$ .or. &
                                       ixa == symmetric_edge$) then
        if (ix1 < 1) then
          print *, 'ERROR IN CREATE_GROUP: START_EDGE OF CONTROLED'
          print *, '      ELEMENT IS AT BEGINNING OF LAT AND CANNOT BE'
          print *, '      VARIED FOR GROUP: ', lat%ele(ix_lord)%name
          call err_exit
        endif
      endif

      if (ixa == end_edge$ .or. ixa == accordion_edge$ .or. &
                                        ixa == symmetric_edge$) then
        if (ix2 > lat%n_ele_track) then
          print *, 'ERROR IN CREATE_GROUP: END_EDGE OF CONTROLED'
          print *, '      ELEMENT IS AT END OF LAT AND CANNOT BE'
          print *, '      VARIED FOR GROUP: ', lat%ele(ix_lord)%name
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

    ! x_limit and y_limit

    elseif (ixa == x_limit$ .or. ixa == y_limit$) then

      if (n_con+2 > size(lat%control)) call reallocate_control (lat, n_con+100)
      lat%control(n_con+1) = contrl(i)
      lat%control(n_con+2) = contrl(i)
      lat%control(n_con+1)%ix_lord = ix_lord
      lat%control(n_con+2)%ix_lord = ix_lord
      if (ixa == x_limit$) then
        lat%control(n_con+1)%ix_attrib = x1_limit$
        lat%control(n_con+2)%ix_attrib = x2_limit$
      else
        lat%control(n_con+1)%ix_attrib = y1_limit$
        lat%control(n_con+2)%ix_attrib = y2_limit$
      endif
      n_con = n_con + 2

    ! x_limit and y_limit

    elseif (ixa == aperture$) then

      if (n_con+4 > size(lat%control)) call reallocate_control (lat, n_con+100)
      lat%control(n_con+1:n_con+4) = contrl(i)
      lat%control(n_con+1:n_con+4)%ix_lord = ix_lord
      lat%control(n_con+1)%ix_attrib = x1_limit$
      lat%control(n_con+2)%ix_attrib = x2_limit$
      lat%control(n_con+3)%ix_attrib = y1_limit$
      lat%control(n_con+4)%ix_attrib = y2_limit$
      n_con = n_con + 4

    ! For all else without position control the group setup is simple.

    else

      free = attribute_free (contrl(i)%ix_slave, ixa, lat, err_print_flag)
      err = err .or. .not. free
      n_con = n_con + 1
      if (n_con > size(lat%control)) call reallocate_control (lat, n_con+100)
      lat%control(n_con) = contrl(i)
      lat%control(n_con)%ix_lord = ix_lord

    endif

  enddo

! final bookkeping

  lat%ele(ix_lord)%ix2_slave = n_con
  lat%ele(ix_lord)%n_slave = n_con - lat%ele(ix_lord)%ix1_slave + 1
  lat%n_control_max = n_con

!---------------------------------------------------------------------------

contains

subroutine bookit (i_ele, scale)

  integer scale, i_ele

  n_con = n_con + 1
  if (n_con > size(lat%control)) call reallocate_control (lat, n_con+100)
  lat%control(n_con)%ix_lord = ix_lord
  lat%control(n_con)%ix_slave = i_ele
  lat%control(n_con)%ix_attrib = l$
  lat%control(n_con)%coef = scale * contrl(i)%coef

end subroutine

end subroutine
