!+
! Subroutine SPLIT_RING (RING, S_SPLIT, IX_SPLIT, SPLIT_DONE)
!
! Subroutine to split a ring at a point. Subroutine will not split the ring
! if the split would create a "runt" element with length less than 10um
!
! Note: split_ring does NOT call make_mat6. The Twiss parameters are also
!       not recomputed.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     RING    -- Ring_struct: Original ring structure.
!     S_SPLIT -- Real: longitudinal distance at which ring is to be split.
!
! Output:
!     RING       -- Ring_struct: Modified ring structure.
!     IX_SPLIT   -- Integer: Index of element just before the split.
!     SPLIT_DONE -- Logical: .true. if ring was split
!-


subroutine split_ring (ring, s_split, ix_split, split_done)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct)  ring
  type (ele_struct)  ele

  real s_split, len_orig, len1, len2, coef1, coef2, angle0, coef_old

  integer i, j, k, ix, ix1, ix_del, ix1_del, ix2_del, ixx1
  integer ix_split, ix_lord, ixc, ix_attrib, ix_super_lord
  integer icon, ix2, inc

  logical split_done

! check for s_split out of bounds

  if (s_split < 0 .or. s_split > ring%param%total_length) then
    type *, 'ERROR IN SPLIT_RING: S_SPLIT NOT VALID: ', s_split
    call err_exit
  endif

! find where to split

  do ix_split = 1, ring%n_ele_ring
    if (abs(ring%ele_(ix_split)%s - s_split) < 1.0e-5) then
      split_done = .false.
      return
    endif
    if (ring%ele_(ix_split)%s > s_split) exit
  enddo

  split_done = .true.
  ele = ring%ele_(ix_split)
  len_orig = ele%value(l$)
  len2 = ring%ele_(ix_split)%s - s_split
  len1 = len_orig - len2

! insert a new element
! note: ring.control_().ix_ele pointers to ix_split will now point to ix_split+1

  ele%value(l$) = 0       ! so no s recalc with insert_element
  call insert_element (ring, ele, ix_split)

  ix = len_trim(ele%name)
  ring%ele_(ix_split)%name   = ele%name(:ix) // '\1'
  ring%ele_(ix_split+1)%name = ele%name(:ix) // '\2'

! if an sbend:
!     1) renormalize the angles
!     2) zero the face angles next to the split

  if (ele%key == sbend$) then
    angle0 = ele%value(angle$)
    ring%ele_(ix_split)%value(angle$) = angle0 * len1 / len_orig
    ring%ele_(ix_split)%value(e2$) = 0.0
    ring%ele_(ix_split+1)%value(angle$) = angle0 * len2 / len_orig
    ring%ele_(ix_split+1)%value(e1$) = 0.0
  endif                       

! hkicks and vkicks get distributed in proportion to the length

  ring%ele_(ix_split)%value(hkick$) = ele%value(hkick$) * len1 / len_orig
  ring%ele_(ix_split)%value(vkick$) = ele%value(vkick$) * len1 / len_orig

  ring%ele_(ix_split+1)%value(hkick$) = ele%value(hkick$) * len2 / len_orig
  ring%ele_(ix_split+1)%value(vkick$) = ele%value(vkick$) * len2 / len_orig

! put in correct lengths and s positions

  ring%ele_(ix_split)%value(l$) = len1
  ring%ele_(ix_split)%s = s_split
  ring%ele_(ix_split+1)%value(l$) = len2

!-------------------------------------------------------------
! Now to correct the slave/lord bookkeeping...

  ix_super_lord = 0   ! no super lord made yet.

! a free drift or free bend needs nothing more.

  if (ele%key == drift$ .and. ele%control_type == free$) goto 8000
  if (ele%key == sbend$ .and. ele%control_type == free$) goto 8000

! If we have split a super_slave we need to make a 2nd control list for one
! of the split elements (can't have both split elements using the same list).
! Also: Redo the control list for the lord elements.

! A split bend element gets handeled the same as a split super_slave.
! We do not try to create a super_lord for a split bend because of the
! problem of handling the face angles.

! Note: If the super_slave is a bend then we are in trouble since the face
! angles cannot be handled correctly.

  if (ele%control_type == super_slave$ .and. ele%key == sbend$) then
    type *, 'ERROR IN SPLIT_RING: BEND ELEMENT TO SPLIT IS A SUPER_SLAVE!'
    type *, '      I DO NOT KNOW HOW TO HANDLE THIS!'
    type *, '      PLEASE SEAK EXPERT (HUMAN) HELP!'
    call err_exit
  endif

  if (ele%control_type == super_slave$ .or. ele%key == sbend$) then

    if (ele%n_lord == 0) goto 8000  ! nothing to do for free sbend$

    ixc = ring%n_ic_array
    ring%ele_(ix_split)%ic1_lord = ixc + 1
    ring%ele_(ix_split)%ic2_lord = ixc + ele%n_lord
    ring%n_ic_array = ixc + ele%n_lord

    do j = 1, ele%n_lord

      ix = ring%ele_(ix_split+1)%ic1_lord - 1
      icon = ring%ic_(ix + j)

      coef_old = ring%control_(icon)%coef
      ix_attrib = ring%control_(icon)%ix_attrib
      ix_lord = ring%control_(icon)%ix_lord

      if (ele%control_type == super_slave$ .or.  &
            ix_attrib == hkick$ .or. ix_attrib == vkick$) then
        coef1 = coef_old * len1 / len_orig
        coef2 = coef_old * len2 / len_orig
      else
        coef1 = coef_old
        coef2 = coef_old
      endif

      ring%control_(icon)%coef = coef2

      ring%ele_(ix_lord)%n_slave = ring%ele_(ix_lord)%n_slave + 1
      call adjust_control_struct (ring, ix_lord)

      ix2 = ring%ele_(ix_lord)%ix2_slave
      ring%control_(ix2)%ix_slave = ix_split
      ring%control_(ix2)%ix_attrib = ix_attrib
      ring%control_(ix2)%coef = coef1
      ring%ic_(ixc+j) = ix2

      if (ring%ele_(ix_lord)%control_type == super_lord$) &
                    call order_super_lord_slaves (ring, ix_lord)

    enddo

    goto 8000   ! and return

  endif   ! split element is a super_slave

! Here if a free or overlay element, not a bend.
! Need to make a super lord to control the split elements.

  ix_super_lord = ring%n_ele_max + 1
  if (ix_super_lord > n_ele_maxx) then
    type *, 'ERROR IN SPLIT_RING: NOT ENOUGH RING ELEMENTS!!!'
    type *, '      YOU NEED TO INCREASE N_ELE_MAXX IN BMAD_STRUCT!!!'
    call err_exit
  endif
                  
  ring%n_ele_max = ix_super_lord
  ring%ele_(ix_super_lord) = ele
  ring%ele_(ix_super_lord)%control_type = super_lord$
  ring%ele_(ix_super_lord)%value(l$) = len_orig
  ixc = ring%n_control_array + 1
  ring%ele_(ix_super_lord)%ix1_slave = ixc
  ring%ele_(ix_super_lord)%ix2_slave = ixc + 1
  ring%ele_(ix_super_lord)%n_slave = 2
  ring%n_control_array = ixc + 1
  ring%control_(ixc)%ix_lord = ix_super_lord
  ring%control_(ixc)%ix_slave = ix_split
  ring%control_(ixc)%coef = len1 / len_orig
  ring%control_(ixc+1)%ix_lord = ix_super_lord
  ring%control_(ixc+1)%ix_slave = ix_split + 1
  ring%control_(ixc+1)%coef = len2 / len_orig

! overlay lord elements of the split element must now point towards the
! super lord

  do i = ele%ic1_lord, ele%ic2_lord
    j = ring%ic_(i)
    ix_lord = ring%control_(j)%ix_lord
    do k = ring%ele_(ix_lord)%ix1_slave, ring%ele_(ix_lord)%ix2_slave
      if (ring%control_(k)%ix_slave == ix_split+1) then
        ring%control_(k)%ix_slave = ix_super_lord
      endif
    enddo
  enddo

! split elements must now be pointing towards their lord

  ring%ele_(ix_split)%control_type = super_slave$
  inc = ring%n_ic_array + 1
  ring%ele_(ix_split)%ic1_lord = inc
  ring%ele_(ix_split)%ic2_lord = inc
  ring%ele_(ix_split)%n_lord = 1
  ring%n_ic_array = inc
  ring%ic_(inc) = ixc

  ring%ele_(ix_split+1)%control_type = super_slave$
  inc = ring%n_ic_array + 1
  ring%ele_(ix_split+1)%ic1_lord = inc
  ring%ele_(ix_split+1)%ic2_lord = inc
  ring%ele_(ix_split+1)%n_lord = 1
  ring%n_ic_array = inc
  ring%ic_(inc) = ixc + 1

! last details:
!     1) Groups that point to the split element must be redirected to the lord
!     2) Call control_bookkeeper to remake the split elements
!     3) And return

8000  continue

  do i = ring%n_ele_ring+1, ring%n_ele_max
    if (ring%ele_(i)%control_type == group_lord$) then
      do j = ring%ele_(i)%ix1_slave, ring%ele_(i)%ix2_slave
        if (ring%control_(j)%ix_slave == ix_split+1) then
          if (ring%control_(j)%ix_attrib == l$) then
            type *, 'WARNING IN SPLIT_RING: GROUP: ', ring%ele_(i)%name
            type *, '        CONTROLS L$ OF SPLIT ELEMENT: ', ele%name
          elseif (ix_super_lord /= 0) then
            ring%control_(j)%ix_slave = ix_super_lord
          else
            type *, 'ERROR IN SPLIT_RING: GROUP: ', ring%ele_(i)%name
            type *, '      CONTROLS SPLIT ELEMENT: ', ele%name
            type *, '      BUT NO LORD WAS MADE!'
            call err_exit
          endif
        endif
      enddo
    endif
  enddo


  if (ring%n_control_array > n_control_maxx) then
    type *, 'ERROR IN SPLIT_RING: NOT ENOUGH CONTROL ELEMENTS !!!'
    type *, '      YOU NEED TO INCREASE N_CONTROL_MAXX IN BMAD_STRUCT !!!'
    call err_exit
  endif

  if (ele%control_type /= free$) then
    call control_bookkeeper (ring, ix_split)
    call control_bookkeeper (ring, ix_split+1)
  endif

end subroutine
