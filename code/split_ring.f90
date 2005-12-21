!+
! Subroutine split_ring (ring, s_split, ix_split, split_done)
!
! Subroutine to split a ring at a point. Subroutine will not split the ring
! if the split would create a "runt" element with length less than 10um
!
! Note: split_ring does NOT call make_mat6. The Twiss parameters are also
!       not recomputed.
!
! Modules Needed:
!   use bmad
!
! Input:
!     ring    -- Ring_struct: Original ring structure.
!     s_split -- Real(rp): Position at which ring is to be split.
!
! Output:
!     ring       -- Ring_struct: Modified ring structure.
!     ix_split   -- Integer: Index of element just before the split.
!     split_done -- Logical: True if ring was split.
!-

#include "CESR_platform.inc"

subroutine split_ring (ring, s_split, ix_split, split_done)

  use bmad_struct
  use bmad_interface, except => split_ring
  use bookkeeper_mod, only: control_bookkeeper
  use nrutil, only: reallocate

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), save :: ele
  type (ele_struct), pointer :: ele1, ele2

  real(rp) s_split, len_orig, len1, len2, coef1, coef2, coef_old

  integer i, j, k, ix
  integer ix_split, ix_lord, ixc, ix_attrib, ix_super_lord
  integer icon, ix2, inc, nr, n_ic2, ct

  logical split_done
  character(16) :: r_name = "split_ring"

! Check for s_split out of bounds.

  nr = ring%n_ele_use
  if (s_split < ring%ele_(0)%s .or. s_split > ring%ele_(nr)%s) then
    print *, 'ERROR IN SPLIT_RING: POSITION OF SPLIT NOT WITHIN RING: ', s_split
    call err_exit
  endif

! Find where to split.

  do ix_split = 0, ring%n_ele_use   
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

! there is a problem with custom elements in that we don't know which
! attributes (if any) scale with length.

  if (ele%key == custom$) then
    print *, "ERROR IN SPLIT_RING: I DON'T KNOW HOW TO SPLIT A CUSTOM ELEMENT!"
    call err_exit
  endif

! Insert a new element.
! Note: Any ring%control_()%ix_ele pointing to ix_split will now 
!  point to ix_split+1.

  ele%value(l$) = 0       ! so no s recalc with insert_element
  call insert_element (ring, ele, ix_split)
  ele1 => ring%ele_(ix_split)
  ele2 => ring%ele_(ix_split+1)

  ix = len_trim(ele%name)
  ele1%name = ele%name(:ix) // '\1'
  ele2%name = ele%name(:ix) // '\2'

! kill any talyor series

  if (associated (ele1%taylor(1)%term)) call kill_taylor (ele1%taylor)
  if (associated (ele2%taylor(1)%term)) call kill_taylor (ele2%taylor)

! put in correct lengths and s positions

  ele1%value(l$) = len1
  ele1%s = s_split
  ele2%value(l$) = len2

! correct aperture limits

  select case (ele%aperture_at)
  case (entrance_end$) 
    ele2%value(x_limit$) = 0
    ele2%value(y_limit$) = 0
    ele2%value(aperture$) = 0

  case (exit_end$)
    ele1%value(x_limit$) = 0
    ele1%value(y_limit$) = 0
    ele1%value(aperture$) = 0

  case (both_ends$)
    ele1%aperture_at = entrance_end$
    ele2%aperture_at = exit_end$

  case (no_end$)

  case default
    call out_io (s_abort$, r_name, 'BAD ELE%APERTURE_AT VALUE: ' // &
                                     element_end_name(ele%aperture_at))
    call err_exit
  end select

! correct coupler location

  if (ele%key == lcavity$) then
    select case (ele%coupler_at)
    case (entrance_end$) 
      ele2%value(coupler_strength$) = 0

    case (exit_end$)
      ele1%value(coupler_strength$) = 0

    case (both_ends$)
      ele1%coupler_at = entrance_end$
      ele2%coupler_at = exit_end$

    case (no_end$)

    case default
      call out_io (s_abort$, r_name, 'BAD ELE%COUPLER_AT VALUE: ' // &
                                       element_end_name(ele%coupler_at))
      call err_exit
    end select
  endif

!-------------------------------------------------------------
! Now to correct the slave/lord bookkeeping...

  ix_super_lord = 0   ! no super lord made yet.

! a free drift needs nothing more.

  if (ele%key == drift$ .and. ele%control_type == free$) goto 8000

! If we have split a super_slave we need to make a 2nd control list for one
! of the split elements (can't have both split elements using the same list).
! Also: Redo the control list for the lord elements.

  if (ele%control_type == super_slave$) then

    if (ele%n_lord == 0) goto 8000  ! nothing to do for free element

    ixc = ring%n_ic_max
    n_ic2 = ixc + ele%n_lord
    ele1%ic1_lord = ixc + 1
    ele1%ic2_lord = n_ic2
    ring%n_ic_max = n_ic2
    if (n_ic2 > size(ring%ic_)) ring%ic_ => reallocate (ring%ic_, n_ic2+500)

    do j = 1, ele%n_lord

      ix = ele2%ic1_lord - 1
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
      call add_lattice_control_structs (ring, ix_lord)

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

! Here if a free or overlay element
! Need to make a super lord to control the split elements.

  ix_super_lord = ring%n_ele_max + 1
  if (ix_super_lord > ubound(ring%ele_, 1)) then
    call allocate_ring_ele_(ring)
    ele1 => ring%ele_(ix_split)
    ele2 => ring%ele_(ix_split+1)
  endif            
  ring%n_ele_max = ix_super_lord
  ring%ele_(ix_super_lord) = ele
  ring%ele_(ix_super_lord)%control_type = super_lord$
  ring%ele_(ix_super_lord)%value(l$) = len_orig
  ixc = ring%n_control_max + 1
  if (ixc+1 > size(ring%control_)) &
                      ring%control_ => reallocate (ring%control_, ixc+500)
  ring%ele_(ix_super_lord)%ix1_slave = ixc
  ring%ele_(ix_super_lord)%ix2_slave = ixc + 1
  ring%ele_(ix_super_lord)%n_slave = 2
  ring%n_control_max = ixc + 1
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

  if (ring%n_ic_max+2 > size(ring%ic_)) &
                    ring%ic_ => reallocate (ring%ic_, ring%n_ic_max+500)

  ele1%control_type = super_slave$
  inc = ring%n_ic_max + 1
  ele1%ic1_lord = inc
  ele1%ic2_lord = inc
  ele1%n_lord = 1
  ring%n_ic_max = inc
  ring%ic_(inc) = ixc

  ele2%control_type = super_slave$
  inc = ring%n_ic_max + 1
  ele2%ic1_lord = inc
  ele2%ic2_lord = inc
  ele2%n_lord = 1
  ring%n_ic_max = inc
  ring%ic_(inc) = ixc + 1

! last details:
!     1) Groups that point to the split element must be redirected to the lord
!     2) Call control_bookkeeper to remake the split elements
!     3) And return

8000  continue

  do i = ring%n_ele_use+1, ring%n_ele_max
    ct = ring%ele_(i)%control_type
    if (ct == group_lord$ .or. ct == i_beam_lord$) then
      do j = ring%ele_(i)%ix1_slave, ring%ele_(i)%ix2_slave
        if (ring%control_(j)%ix_slave == ix_split+1) then
          if (ring%control_(j)%ix_attrib == l$) then
            print *, 'WARNING IN SPLIT_RING: GROUP: ', ring%ele_(i)%name
            print *, '        CONTROLS L$ OF SPLIT ELEMENT: ', ele%name
          elseif (ix_super_lord /= 0) then
            ring%control_(j)%ix_slave = ix_super_lord
          else
            print *, 'ERROR IN SPLIT_RING: GROUP: ', ring%ele_(i)%name
            print *, '      CONTROLS SPLIT ELEMENT: ', ele%name
            print *, '      BUT NO LORD WAS MADE!'
            call err_exit
          endif
        endif
      enddo
    endif
  enddo

  call control_bookkeeper (ring, ix_split)
  call control_bookkeeper (ring, ix_split+1)

  call check_ring_controls (ring, .true.)

end subroutine
