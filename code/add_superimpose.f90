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
  use bmad_interface, except => add_superimpose
  use nrutil, only: reallocate

  implicit none

  type (ring_struct)  ring
  type (ele_struct)  super_ele, sup_ele, slave_ele, drift
  type (control_struct)  sup_con(100)

  real(rp) s1, s2, length, s1_lat, s2_lat

  integer j, jj, k, ix, n, i2, ic, n_con
  integer ix1_split, ix2_split, ix_super, ix_super_con
  integer ix_slave, ixn, ixc, superimpose_key, ix_slave_name

  logical setup_lord, split1_done, split2_done

  character(20) fmt
  character(20) :: r_name = "add_superimpose"

!-------------------------------------------------------------------------
! We need a copy of super_ele since the actual argument may be in the ring
! and split_ring can then overwrite it.

  sup_ele = super_ele
  ix_slave_name = 0

  call init_ele (drift)
  drift%key = drift$

! s1 is the left edge of the superimpose.
! s2 is the right edge of the superimpose.
! For a ring a superimpose can wrap around the ends of the lattice.

  s1_lat = ring%ele_(0)%s               ! normally this is zero.
  s2_lat = ring%ele_(ring%n_ele_use)%s

  s1 = sup_ele%s - sup_ele%value(l$)
  s2 = sup_ele%s                 

  if (s1 < s1_lat .and. ring%param%lattice_type == circular_lattice$) &
                                      s1 = s1 + ring%param%total_length

  if (s2 < s1_lat .or. s1 > s2_lat) then
    call out_io (s_abort$, r_name, 'SUPERIMPOSE POSITION BEYOUND END OF LATTICE')
    call out_io (s_blank$, r_name, 'LEFT EDGE: \F10.1\ ', s1)
    call out_io (s_blank$, r_name, 'RIGHT EDGE:\F10.1\ ', s2)
    call err_exit
  endif
 
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

  if (s2 < s1) then   ! if superimpose wraps around 0 ...
    call split_ring (ring, s2, ix2_split, split2_done)
    if (split2_done) call delete_last_chars (ix2_split)
    call split_ring (ring, s1, ix1_split, split1_done)
    if (split1_done) call delete_last_chars (ix1_split)

  else                ! no wrap case
    if (s1 < s1_lat) then  ! superimpose off end case
      if (ring%ele_(1)%key /= drift$) then
        length = s1_lat - s1
        drift%value(l$) = length
        call insert_element (ring, drift, 1)
        s1 = s1_lat
        s2 = s2 + length
        s2_lat = s2_lat + length
      endif
      ix1_split = 0
    else
      call split_ring (ring, s1, ix1_split, split1_done)
      if (split1_done) call delete_last_chars (ix1_split)
    endif

    if (s2 > s2_lat) then  ! superimpose off end case
      if (ring%ele_(ring%n_ele_use)%key /= drift$) then
        drift%value(l$) = s2 - s2_lat
        call insert_element (ring, drift, ring%n_ele_use + 1)
        s2_lat = s2
      endif
      ix2_split = ring%n_ele_use
    else
      call split_ring (ring, s2, ix2_split, split2_done)
      if (split2_done) call delete_last_chars (ix2_split)
    endif

    if (s1 < s1_lat) ring%ele_(1)%value(l$) = ring%ele_(1)%s - s1
    if (s2 > s2_lat) then
      n = ring%n_ele_use
      ring%ele_(n)%value(l$) = s2 - ring%ele_(n-1)%s
    endif

  endif

! SPLIT_RING adds "\1" and "\2" to the names of the split elements.
! For aesthetic reasons remove the last digit of the names.

  call delete_double_slash (ix1_split)
  call delete_double_slash (ix1_split+1)
  call delete_double_slash (ix2_split)
  call delete_double_slash (ix2_split+1)

! if element overlays a drift then just insert it in the regular ring list

  if (ix2_split == ix1_split + 1 .and. ring%ele_(ix2_split)%key == drift$) then
    ix_super = ix2_split
    ring%ele_(ix_super) = sup_ele
    ring%ele_(ix_super)%control_type = free$
    return
  endif

! Only possibility left means we have to set up a super_lord element for the
! superposition

  ix_super = ring%n_ele_max + 1
  ring%n_ele_max = ix_super
  if (ring%n_ele_max > ubound(ring%ele_, 1)) call allocate_ring_ele_(ring)
  ring%ele_(ix_super) = sup_ele
  ring%ele_(ix_super)%control_type = super_lord$

  ix_super_con = 0
  length = sup_ele%value(l$)

  ix_slave = ix1_split

! Go through the list of elements being superimposed upon.
! Zero length elements (markers and multipoles) do not get involved here.

  do 

    ix_slave = ix_slave + 1
    if (ix_slave == ix2_split + 1) exit

    if (ix_slave == ring%n_ele_use + 1) ix_slave = 1
    slave_ele = ring%ele_(ix_slave)
    if (slave_ele%value(l$) == 0) cycle

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
      if (ixn > ubound(ring%ele_, 1)) call allocate_ring_ele_(ring)
      ring%ele_(ixn) = slave_ele
      ring%ele_(ixn)%control_type = super_lord$
      ring%n_ele_max = ixn
      ixc = ring%n_control_max + 1
      if (ixc > size(ring%control_)) ring%control_ => reallocate(ring%control_, ixc+500)
      ring%ele_(ixn)%ix1_slave = ixc
      ring%ele_(ixn)%ix2_slave = ixc
      ring%ele_(ixn)%n_slave = 1
      ring%control_(ixc)%ix_lord = ixn
      ring%control_(ixc)%ix_slave = ix_slave
      ring%control_(ixc)%coef = 1.0
      ring%n_control_max = ixc

      do j = ring%ele_(ixn)%ic1_lord, ring%ele_(ixn)%ic2_lord
        jj = ring%ic_(j)
        ring%control_(jj)%ix_slave = ixn
      enddo

      ic = ring%n_ic_max + 1
      ring%ele_(ix_slave)%ic1_lord = ic
      ring%ele_(ix_slave)%ic2_lord = ic + 1
      ring%ele_(ix_slave)%n_lord = 2
      ring%n_ic_max = ic + 1
      if (ic+1 > size(ring%ic_)) ring%ic_ => reallocate (ring%ic_, ic+500)
      ring%ic_(ic) = ixc 

    else
      ring%ele_(ix_slave)%n_lord = slave_ele%n_lord + 1
      call add_lattice_control_structs (ring, ix_slave)
    endif

    if (slave_ele%key == drift$) then
      ix_slave_name = ix_slave_name + 1
      write (ring%ele_(ix_slave)%name, '(2a, i0)') &
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

  enddo

! Special case where elements on either side of the superimpose have the same
! name

  if (split1_done .and. split2_done .and. &
                ring%ele_(ix1_split)%name == ring%ele_(ix2_split+1)%name) then
    ring%ele_(ix1_split)%name = trim(ring%ele_(ix1_split)%name) // '1'
    ring%ele_(ix2_split+1)%name = trim(ring%ele_(ix2_split+1)%name) // '2'
  endif

! transfer control info from sup_con array

  ixc = ring%n_control_max + 1
  n_con = ixc + ix_super_con - 1
  if (n_con > size(ring%control_)) ring%control_ => reallocate(ring%control_, n_con+500) 
  ring%ele_(ix_super)%ix1_slave = ixc
  ring%ele_(ix_super)%ix2_slave = n_con
  ring%ele_(ix_super)%n_slave = ix_super_con

  do k = 1, ix_super_con
    ring%control_(k+ixc-1) = sup_con(k)
    ix_slave = ring%control_(k+ixc-1)%ix_slave
    i2 = ring%ele_(ix_slave)%ic2_lord
    ring%ic_(i2) = k+ixc-1
  enddo

  ring%n_control_max = n_con

! order slave elements in the super_lord list to be in the correct order

  call s_calc (ring)  ! just in case superimpose extended before beginning of lattice.
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
!+
! Function superimpose_key (key1, key2) result (key12)
!
! Function to decide what the element key (key12) should be when
! an element with key1 is superimpsed upon with an element
! with key2. 
!
! This function is not meant for general use.
!-

function superimpose_key (key1, key2) result (key12)

  use bmad_struct
  use bmad_interface

  implicit none

  integer key1, key2, key12

!

  key12 = -1  ! Default if no superimpse possible

  if (key1 == key2) then
    key12 = key1
    return
  endif

  select case (key1)
  case (drift$)
    key12 = key2
    return

  case (kicker$, rcollimator$, monitor$, instrument$)
    key12 = key2

  case (quadrupole$,  solenoid$, sol_quad$) 
    select case (key2)
    case (quadrupole$);    key12 = sol_quad$
    case (solenoid$);      key12 = sol_quad$
    case (sol_quad$);      key12 = sol_quad$
    case (bend_sol_quad$); key12 = bend_sol_quad$
    case (sbend$);         key12 = bend_sol_quad$
    end select

  case (bend_sol_quad$)
    select case (key2)
    case (quadrupole$);    key12 = bend_sol_quad$
    case (solenoid$);      key12 = bend_sol_quad$
    case (sol_quad$);      key12 = bend_sol_quad$
    case (sbend$);         key12 = bend_sol_quad$
    end select
  end select

!

  select case (key2)
  case (drift$, kicker$, rcollimator$, monitor$, instrument$)
    key12 = key1
  end select

end function
