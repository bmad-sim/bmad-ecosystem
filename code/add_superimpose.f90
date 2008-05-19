!+
! Subroutine add_superimpose (lat, super_ele, ix_super)
!
! Subroutine to make a superimposed element. If the element can be inserted
! into the lat without making a super_lord element then this will be done.
!
! Modules Needed:
!   use bmad
!
! Input:
!     lat      -- lat_struct: Lat to modify
!     super_ele -- Ele_struct: Element to superimpose
!         %s       -- Position of end of element.
!                      Negative distances mean distance from the end.
!
!
! Output:
!     lat     -- lat_struct: Modified lat.
!     ix_super -- Integer: Index where element is put
!-

#include "CESR_platform.inc"

subroutine add_superimpose (lat, super_ele, ix_super)

  use bmad_struct
  use bmad_interface, except_dummy => add_superimpose

  implicit none

  type (lat_struct)  lat
  type (ele_struct)  super_ele
  type (ele_struct), save :: sup_ele, slave_ele, drift
  type (control_struct)  sup_con(100)

  real(rp) s1, s2, length, s1_lat, s2_lat

  integer i, j, jj, k, ix, n, i2, ic, n_con
  integer ix1_split, ix2_split, ix_super, ix_super_con
  integer ix_slave, ixn, ixc, ix_slave_name

  logical setup_lord, split1_done, split2_done, all_drift

  character(20) fmt
  character(20) :: r_name = "add_superimpose"

!-------------------------------------------------------------------------
! We need a copy of super_ele since the actual argument may be in the lat
! and split_lat can then overwrite it.

  call init_ele (sup_ele)
  call init_ele (slave_ele)
  call init_ele (drift)
  drift%key = drift$

  sup_ele = super_ele
  ix_slave_name = 0

! s1 is the left edge of the superimpose.
! s2 is the right edge of the superimpose.
! For a lat a superimpose can wrap around the ends of the lattice.

  s1_lat = lat%ele(0)%s               ! normally this is zero.
  s2_lat = lat%ele(lat%n_ele_track)%s

  s1 = sup_ele%s - sup_ele%value(l$)
  s2 = sup_ele%s                 

  if (s1 < s1_lat) then
    if (lat%param%lattice_type == linear_lattice$) call out_io (s_warn$, &
           r_name, 'superimpose is being wrapped around for: ' // sup_ele%name)
    s1 = s1 + lat%param%total_length
  endif

  if (s2 < s1_lat .or. s1 > s2_lat) then
    call out_io (s_abort$, r_name, 'SUPERIMPOSE POSITION BEYOUND END OF LATTICE')
    call out_io (s_blank$, r_name, 'LEFT EDGE: \F10.1\ ', s1)
    call out_io (s_blank$, r_name, 'RIGHT EDGE:\F10.1\ ', s2)
    call err_exit
  endif
 
!-------------------------------------------------------------------------
! if element has zero length then just insert it in the tracking part of the lattice list

  if (sup_ele%value(l$) == 0) then
    call split_lat (lat, s1, ix1_split, split1_done)
    call insert_element (lat, sup_ele, ix1_split+1)
    ix_super = ix1_split + 1
    lat%ele(ix_super)%control_type = free$
    return
  endif

!-------------------------------------------------------------------------
! Split lat at begining and end of the superimpose.
! the correct order of splitting is important since we are creating elements
! so that the numbering of the elments after the split changes.

  if (s2 < s1) then   ! if superimpose wraps around 0 ...
    call split_lat (lat, s2, ix2_split, split2_done)
    if (split2_done) call delete_last_chars (ix2_split)
    call split_lat (lat, s1, ix1_split, split1_done)
    if (split1_done) call delete_last_chars (ix1_split)

  else                ! no wrap case
    if (s1 < s1_lat) then  ! superimpose off end case
      if (lat%ele(1)%key /= drift$) then
        length = s1_lat - s1
        drift%value(l$) = length
        call insert_element (lat, drift, 1)
        s1 = s1_lat
        s2 = s2 + length
        s2_lat = s2_lat + length
      endif
      ix1_split = 0
    else
      call split_lat (lat, s1, ix1_split, split1_done)
      if (split1_done) call delete_last_chars (ix1_split)
    endif

    if (s2 > s2_lat) then  ! superimpose off end case
      if (lat%ele(lat%n_ele_track)%key /= drift$) then
        drift%value(l$) = s2 - s2_lat
        call insert_element (lat, drift, lat%n_ele_track + 1)
        s2_lat = s2
      endif
      ix2_split = lat%n_ele_track
    else
      call split_lat (lat, s2, ix2_split, split2_done)
      if (split2_done) call delete_last_chars (ix2_split)
    endif

    if (s1 < s1_lat) lat%ele(1)%value(l$) = lat%ele(1)%s - s1
    if (s2 > s2_lat) then
      n = lat%n_ele_track
      lat%ele(n)%value(l$) = s2 - lat%ele(n-1)%s
    endif

  endif

! split_lat adds "\1" and "\2" to the names of the split elements.
! For aesthetic reasons remove the last digit of the names.

  call delete_double_slash (ix1_split)
  call delete_double_slash (ix1_split+1)
  call delete_double_slash (ix2_split)
  call delete_double_slash (ix2_split+1)

! If element overlays only drifts then just 
! insert it in the tracking part of the lat list.

  all_drift = (ix2_split > ix1_split)
  do i = ix1_split+1, ix2_split
    if (lat%ele(i)%key /= drift$) all_drift = .false.
    if (.not. all_drift) exit
  enddo

  if (all_drift) then  
    do i = ix2_split, ix1_split+2, -1 ! remove all drifts but one
      call remove_ele_from_lat(lat, i)
    enddo
    ix_super = ix1_split + 1
    lat%ele(ix_super) = sup_ele
    lat%ele(ix_super)%control_type = free$
    return
  endif

! Only possibility left means we have to set up a super_lord element for the
! superposition

  ix_super = lat%n_ele_max + 1
  lat%n_ele_max = ix_super
  if (lat%n_ele_max > ubound(lat%ele, 1)) call allocate_ele_array(lat%ele)
  lat%ele(ix_super) = sup_ele
  lat%ele(ix_super)%control_type = super_lord$

  ix_super_con = 0
  length = sup_ele%value(l$)

  ix_slave = ix1_split

! Go through the list of elements being superimposed upon.
! Zero length elements (markers and multipoles) do not get involved here.

  do 

    ix_slave = ix_slave + 1
    if (ix_slave == ix2_split + 1) exit

    if (ix_slave == lat%n_ele_track + 1) ix_slave = 1
    slave_ele = lat%ele(ix_slave)
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
      ixn = lat%n_ele_max + 1
      if (ixn > ubound(lat%ele, 1)) call allocate_ele_array(lat%ele)
      lat%ele(ixn) = slave_ele
      lat%ele(ixn)%control_type = super_lord$
      lat%n_ele_max = ixn
      ixc = lat%n_control_max + 1
      if (ixc > size(lat%control)) call reallocate_control(lat, ixc+100)
      lat%ele(ixn)%ix1_slave = ixc
      lat%ele(ixn)%ix2_slave = ixc
      lat%ele(ixn)%n_slave = 1
      lat%control(ixc)%ix_lord = ixn
      lat%control(ixc)%ix_slave = ix_slave
      lat%control(ixc)%coef = 1.0
      lat%n_control_max = ixc

      do j = lat%ele(ixn)%ic1_lord, lat%ele(ixn)%ic2_lord
        jj = lat%ic(j)
        lat%control(jj)%ix_slave = ixn
      enddo

      ic = lat%n_ic_max + 1
      lat%ele(ix_slave)%ic1_lord = ic
      lat%ele(ix_slave)%ic2_lord = ic + 1
      lat%ele(ix_slave)%n_lord = 2
      lat%n_ic_max = ic + 1
      lat%ic(ic) = ixc 

    else
      lat%ele(ix_slave)%n_lord = slave_ele%n_lord + 1
      call add_lattice_control_structs (lat, ix_slave)
    endif

    if (slave_ele%key == drift$) then
      ix_slave_name = ix_slave_name + 1
      write (lat%ele(ix_slave)%name, '(2a, i0)') &
                                   trim(sup_ele%name), '\', ix_slave_name
    else
      lat%ele(ix_slave)%name = trim(lat%ele(ix_slave)%name) //  &
                                                         '\' // sup_ele%name
    endif

    call delete_double_slash (ix_slave)

    lat%ele(ix_slave)%control_type = super_slave$

! add control info for main super lord to list

    ix_super_con = ix_super_con + 1
    sup_con(ix_super_con)%ix_slave = ix_slave
    sup_con(ix_super_con)%ix_lord = ix_super
    sup_con(ix_super_con)%coef = slave_ele%value(l$) / length
    sup_con(ix_super_con)%ix_attrib = 0

! change the element key

    call calc_superimpose_key(slave_ele, sup_ele, lat%ele(ix_slave))
    if (lat%ele(ix_slave)%key <= 0) then
      print *, 'ERROR IN ADD_SUPERIMPOSE: BAD SUPERIMPOSE FOR ',  &
                                        sup_ele%name
      print *, '      SUPERIMPOSED UPON: ', slave_ele%name
      call err_exit                    
    endif

  enddo

! Special case where elements on either side of the superimpose have the same
! name

  if (split1_done .and. split2_done .and. &
                lat%ele(ix1_split)%name == lat%ele(ix2_split+1)%name) then
    lat%ele(ix1_split)%name = trim(lat%ele(ix1_split)%name) // '1'
    lat%ele(ix2_split+1)%name = trim(lat%ele(ix2_split+1)%name) // '2'
  endif

! transfer control info from sup_con array

  ixc = lat%n_control_max + 1
  n_con = ixc + ix_super_con - 1
  if (n_con > size(lat%control)) call reallocate_control(lat, n_con+500) 
  lat%ele(ix_super)%ix1_slave = ixc
  lat%ele(ix_super)%ix2_slave = n_con
  lat%ele(ix_super)%n_slave = ix_super_con

  do k = 1, ix_super_con
    lat%control(k+ixc-1) = sup_con(k)
    ix_slave = lat%control(k+ixc-1)%ix_slave
    i2 = lat%ele(ix_slave)%ic2_lord
    lat%ic(i2) = k+ixc-1
  enddo

  lat%n_control_max = n_con

! order slave elements in the super_lord list to be in the correct order

  call s_calc (lat)  ! just in case superimpose extended before beginning of lattice.
  call order_super_lord_slaves (lat, ix_super)

!------------------------------------------------------------------------------
contains

subroutine delete_last_chars (ix_split)

  integer ix_split

  ix = len_trim(lat%ele(ix_split)%name) - 1
  lat%ele(ix_split)%name = lat%ele(ix_split)%name(1:ix)
  ix = len_trim(lat%ele(ix_split+1)%name) - 1
  lat%ele(ix_split+1)%name = lat%ele(ix_split+1)%name(1:ix)

end subroutine

!------------------------------------------------------------------------------
! contains

subroutine delete_double_slash(ix_ele)

  integer ix_ele

  ix = index(lat%ele(ix_ele)%name, '\\')
  if (ix /= 0) lat%ele(ix_ele)%name = lat%ele(ix_ele)%name(:ix) // &
                                              lat%ele(ix_ele)%name(ix+2:)

end subroutine

end subroutine
