!+
! Subroutine add_superimpose (ring, super_ele, ix_super)
!
! Subroutine to make a superimposed element. If the element can be inserted
! into the ring without making a super_lord element then this will be done.
!
! Modules Needed:
!   use bmad
!
! Input:
!     ring      -- Ring_struct: Ring to modify
!     super_ele -- Ele_struct: Element to superimpose
!         %s       -- Position of end of element.
!                      Negative distances mean distance from the end.
!
!
! Output:
!     ring     -- Ring_struct: Modified ring.
!     ix_super -- Integer: Index where element is put
!-

#include "CESR_platform.inc"

subroutine add_superimpose (ring, super_ele, ix_super)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct)  ring
  type (ele_struct)  super_ele, sup_ele, slave_ele
  type (control_struct)  sup_con(100)

  real(rdef) s1, s2, length, s0

  integer j, jj, k, idel, ix, n, i2, ic
  integer ix1_split, ix2_split, ix_super, ix_super_con
  integer ix_slave, ixn, ixc, superimpose_key, ix_slave_name

  logical setup_lord, split1_done, split2_done

!-------------------------------------------------------------------------
! We need a copy of super_ele since the actual argument may be in the ring
! and split_ring can then overwrite it.

  sup_ele = super_ele
  ix_slave_name = 0

! s1 is the left edge of the superimpose.

  s0 = ring%ele_(0)%s   ! normally this is zero.
  s1 = sup_ele%s - sup_ele%value(l$)
  if (s1 < s0) s1 = s1 + ring%param%total_length

!-------------------------------------------------------------------------
! if element has zero length then just insert it in the regular ring list

  if (sup_ele%value(l$) == 0) then
    call split_ring (ring, s1, ix1_split, split1_done)
    call insert_element (ring, sup_ele, ix1_split+1)
    ix_super = ix1_split + 1
    ring%ele_(ix_super)%control_type = free$
    return
  endif

!-------------------------------------------------------------------------
! Split ring at begining and end of the superimpose.
! the correct order of splitting is important since we are creating elements
! so that the numbering of the elments after the split changes.
! SPLIT_RING adds "\1" and "\2" to the names of the split elements.
! For aesthetic reasons remove the last digit of the names

  s2 = sup_ele%s                 
  if (s2 < s1) then   ! if superimpose wraps around 0 ...
    call split_ring (ring, s2, ix2_split, split2_done)
    if (split2_done) call delete_last_chars (ix2_split)
    call split_ring (ring, s1, ix1_split, split1_done)
    if (split1_done) call delete_last_chars (ix1_split)
  else                ! no wrap case
    call split_ring (ring, s1, ix1_split, split1_done)
    if (split1_done) call delete_last_chars (ix1_split)
    call split_ring (ring, s2, ix2_split, split2_done)
    if (split2_done) call delete_last_chars (ix2_split)
  endif

  call delete_double_slash (ix1_split)
  call delete_double_slash (ix1_split+1)
  call delete_double_slash (ix2_split)
  call delete_double_slash (ix2_split+1)

! if element overlays a drift then just insert it in the regular ring list

  if (ix2_split == ix1_split + 1 .and.  &
                              ring%ele_(ix2_split)%key == drift$) then
    ix_super = ix2_split
    ring%ele_(ix_super) = sup_ele
    ring%ele_(ix_super)%control_type = free$
    return
  endif

! Only possibility left means we have to set up a super_lord element for the
! superposition

  ix_super = ring%n_ele_max + 1
  ring%n_ele_max = ix_super
  if (ring%n_ele_max > n_ele_maxx) then
    print *, 'ERROR IN ADD_SUPERIMPOSE: NOT ENOUGH RING ELEMENTS!!!'
    print *, '      YOU NEED TO INCREASE N_ELE_MAXX IN BMAD_STRUCT !!!'
  endif
  ring%ele_(ix_super) = sup_ele
  ring%ele_(ix_super)%control_type = super_lord$

  ix_super_con = 0
  length = sup_ele%value(l$)

  ix_slave = ix1_split

! Go through the list of elements being superimposed upon.
! Zero length elements (markers and multipoles) do not get involved here.

  do 

    ix_slave = ix_slave + 1
    if (ix_slave == ring%n_ele_ring + 1) ix_slave = 1
    slave_ele = ring%ele_(ix_slave)
    if (slave_ele%value(l$) == 0) goto 8000     ! skip rest of loop

! Do we need to set up a super lord to control this slave element?

    if (slave_ele%control_type == overlay_slave$) then
      setup_lord = .true.
    elseif (slave_ele%control_type == super_slave$) then
      setup_lord = .false.
    elseif (slave_ele%key == drift$) then
      setup_lord = .false.
    else
      setup_lord = .true.
    endif

! if yes then create the super lord element

    if (setup_lord) then
      ixn = ring%n_ele_max + 1
      if (ring%n_ele_max > n_ele_maxx) then
        print *, 'ERROR IN ADD_SUPERIMPOSE: NOT ENOUGH RING ELEMENTS!!!'
        print *, '      YOU NEED TO INCREASE N_ELE_MAXX IN BMAD_STRUCT!!!'
      endif
      ring%ele_(ixn) = slave_ele
      ring%ele_(ixn)%control_type = super_lord$
      ring%n_ele_max = ixn
      ixc = ring%n_control_array + 1
      ring%ele_(ixn)%ix1_slave = ixc
      ring%ele_(ixn)%ix2_slave = ixc
      ring%ele_(ixn)%n_slave = 1
      ring%control_(ixc)%ix_lord = ixn
      ring%control_(ixc)%ix_slave = ix_slave
      ring%control_(ixc)%coef = 1.0
      ring%n_control_array = ixc

      do j = ring%ele_(ixn)%ic1_lord, ring%ele_(ixn)%ic2_lord
        jj = ring%ic_(j)
        ring%control_(jj)%ix_slave = ixn
      enddo

      ic = ring%n_ic_array + 1
      ring%ele_(ix_slave)%ic1_lord = ic
      ring%ele_(ix_slave)%ic2_lord = ic + 1
      ring%ele_(ix_slave)%n_lord = 2
      ring%n_ic_array = ic + 1
      ring%ic_(ic) = ixc 

    else
      ring%ele_(ix_slave)%n_lord = slave_ele%n_lord + 1
      call adjust_control_struct (ring, ix_slave)
    endif

    if (slave_ele%key == drift$) then
      ix_slave_name = ix_slave_name + 1
      n = log10(1.001*ix_slave_name) + 1
      write (ring%ele_(ix_slave)%name, '(2a, i<n>)') &
                                   trim(sup_ele%name), '\', ix_slave_name
    else
      ring%ele_(ix_slave)%name = trim(ring%ele_(ix_slave)%name) //  &
                                                         '\' // sup_ele%name
    endif

    call delete_double_slash (ix_slave)

    ring%ele_(ix_slave)%control_type = super_slave$

! add control info for main super lord to list

    ix_super_con = ix_super_con + 1
    sup_con(ix_super_con)%ix_slave = ix_slave
    sup_con(ix_super_con)%ix_lord = ix_super
    sup_con(ix_super_con)%coef = slave_ele%value(l$) / length
    sup_con(ix_super_con)%ix_attrib = 0

! change the element key

    ring%ele_(ix_slave)%key = superimpose_key(slave_ele%key, sup_ele%key)
    if (ring%ele_(ix_slave)%key <= 0) then
      print *, 'ERROR IN ADD_SUPERIMPOSE: BAD SUPERIMPOSE FOR ',  &
                                        sup_ele%name
      print *, '      SUPERIMPOSED UPON: ', slave_ele%name
      call err_exit                    
    endif

8000    continue
    if (ix_slave == ix2_split) exit

  enddo

! Special case where elements on either side of the superimpose have the same
! name

  if (split1_done .and. split2_done .and. &
                ring%ele_(ix1_split)%name == ring%ele_(ix2_split+1)%name) then
    ring%ele_(ix1_split)%name = trim(ring%ele_(ix1_split)%name) // '1'
    ring%ele_(ix2_split+1)%name = trim(ring%ele_(ix2_split+1)%name) // '2'
  endif

! transfer control info from sup_con array

  ixc = ring%n_control_array + 1
  ring%ele_(ix_super)%ix1_slave = ixc
  ring%ele_(ix_super)%ix2_slave = ixc + ix_super_con - 1
  ring%ele_(ix_super)%n_slave = ix_super_con

  do k = 1, ix_super_con
    ring%control_(k+ixc-1) = sup_con(k)
    ix_slave = ring%control_(k+ixc-1)%ix_slave
    i2 = ring%ele_(ix_slave)%ic2_lord
    ring%ic_(i2) = k+ixc-1
  enddo

  ring%n_control_array = ring%ele_(ix_super)%ix2_slave

  if (ring%n_control_array > n_control_maxx) then
    print *, 'ERROR IN ADD_SUPERIMPOSE: NOT ENOUGH CONTROL ELEMENTS !!!'
    print *, '      YOU NEED TO INCREASE N_CONTROL_MAXX IN BMAD_STRUCT !!!'
    call err_exit
  endif

! order slave elements in the super_lord list to be in the correct order

  call order_super_lord_slaves (ring, ix_super)

!------------------------------------------------------------------------------

contains

subroutine delete_last_chars (ix_split)

  integer ix_split

  ix = len_trim(ring%ele_(ix_split)%name) - 1
  ring%ele_(ix_split)%name = ring%ele_(ix_split)%name(1:ix)
  ix = len_trim(ring%ele_(ix_split+1)%name) - 1
  ring%ele_(ix_split+1)%name = ring%ele_(ix_split+1)%name(1:ix)

end subroutine

!------------------------------------------------------------------------------

subroutine delete_double_slash(ix_ele)

  integer ix_ele

  ix = index(ring%ele_(ix_ele)%name, '\\')
  if (ix /= 0) ring%ele_(ix_ele)%name = ring%ele_(ix_ele)%name(:ix) // &
                                              ring%ele_(ix_ele)%name(ix+2:)

end subroutine

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

integer function superimpose_key (key1, key2)

  use bmad_struct
  use bmad_interface

  implicit none

  integer key1, key2

!

  if (key1 == drift$) then
    superimpose_key = key2
  elseif (key1 == key2) then
    superimpose_key = key1
  elseif (key2 == drift$) then
    superimpose_key = key1
  elseif ((key1 == quadrupole$  .or. key1 == solenoid$ .or.  &
          key1 == sol_quad$) .and. (key2 == quadrupole$  .or.  &
          key2 == solenoid$ .or. key2 == sol_quad$)) then
    superimpose_key = sol_quad$
  else
    superimpose_key = -1
  endif

end function
