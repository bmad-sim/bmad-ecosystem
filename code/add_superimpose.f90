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
!   lat       -- lat_struct: Lat to modify.
!   super_ele -- Ele_struct: Element to superimpose.
!         %s       -- Position of end of element.
!                      Negative distances mean distance from the end.
!
!
! Output:
!   lat      -- lat_struct: Modified lat.
!   ix_super -- Integer: Index where element is put.
!-

#include "CESR_platform.inc"

subroutine add_superimpose (lat, super_ele, ix_super)

  use bmad_struct
  use bmad_interface, except_dummy => add_superimpose

  implicit none

  type (lat_struct), target :: lat
  type (ele_struct)  super_ele
  type (ele_struct), save :: super_saved, slave_saved, drift, null_ele
  type (ele_struct), pointer :: slave, lord
  type (control_struct)  sup_con(100)

  real(rp) s1, s2, length, s1_lat, s2_lat

  integer i, j, jj, k, ix, n, i2, ic, n_con
  integer ix1_split, ix2_split, ix_super, ix_super_con
  integer ix_slave, ixn, ixc, ix_1lord, n_ele_max_old

  logical setup_lord, split1_done, split2_done, all_drift

  character(100) name
  character(20) fmt
  character(20) :: r_name = "add_superimpose"

  !-------------------------------------------------------------------------
  ! Check for negative length

  if (super_ele%value(l$) < 0) then
    call out_io (s_abort$, r_name, &
                    'Superposition of element with negative length not allowed!', &
                    'Element: ' // super_ele%name, &
                    'Length: \es10.2\ ', r_array = (/ super_ele%value(l$) /) )
    call err_exit
  endif

  ! We need a copy of super_ele since the actual argument may be in the lat
  ! and split_lat can then overwrite it.

  call init_ele (super_saved)
  call init_ele (slave_saved)
  call init_ele (drift)
  drift%key = drift$

  super_saved = super_ele

  ! s1 is the left edge of the superimpose.
  ! s2 is the right edge of the superimpose.
  ! For a lat a superimpose can wrap around the ends of the lattice.

  n_ele_max_old = lat%n_ele_max

  s1_lat = lat%ele(0)%s                 ! normally this is zero.
  s2_lat = lat%ele(lat%n_ele_track)%s

  s1 = super_saved%s - super_saved%value(l$)
  s2 = super_saved%s                 

  if (s1 < s1_lat) then
    if (lat%param%lattice_type == linear_lattice$) call out_io (s_warn$, &
           r_name, 'superimpose is being wrapped around for: ' // super_saved%name)
    s1 = s1 + lat%param%total_length
  endif

  if (s2 < s1_lat .or. s1 > s2_lat) then
    call out_io (s_abort$, r_name, &
      'SUPERIMPOSE POSITION BEYOUND END OF LATTICE', &
      'LEFT EDGE: \F10.1\ ', &
      'RIGHT EDGE:\F10.1\ ', r_array = (/ s1, s2 /))
    call err_exit
  endif
 
  !-------------------------------------------------------------------------
  ! If element has zero length then just insert it in the tracking part 
  ! of the lattice list.

  if (super_saved%value(l$) == 0) then
    call split_lat (lat, s1, ix1_split, split1_done, check_controls = .false.)
    call insert_element (lat, super_saved, ix1_split+1)
    ix_super = ix1_split + 1
    lat%ele(ix_super)%lord_status  = free$
    lat%ele(ix_super)%slave_status = free$
    call adjust_slave_names ()
    return
  endif

  !-------------------------------------------------------------------------
  ! Split lat at begining and end of the superimpose.
  ! the correct order of splitting is important since we are creating elements
  ! so that the numbering of the elments after the split changes.

  if (s2 < s1) then     ! if superimpose wraps around 0 ...
    call split_lat (lat, s2, ix2_split, split2_done, .false., .false.)
    call split_lat (lat, s1, ix1_split, split1_done, .false., .false.)

  else                  ! no wrap case
    if (s1 < s1_lat) then    ! superimpose off end case
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
      call split_lat (lat, s1, ix1_split, split1_done, .false., .false.)
    endif

    if (s2 > s2_lat) then    ! superimpose off end case
      if (lat%ele(lat%n_ele_track)%key /= drift$) then
        drift%value(l$) = s2 - s2_lat
        call insert_element (lat, drift, lat%n_ele_track + 1)
        s2_lat = s2
      endif
      ix2_split = lat%n_ele_track
    else
      call split_lat (lat, s2, ix2_split, split2_done, .false., .false.)
    endif

    if (s1 < s1_lat) lat%ele(1)%value(l$) = lat%ele(1)%s - s1
    if (s2 > s2_lat) then
      n = lat%n_ele_track
      lat%ele(n)%value(l$) = s2 - lat%ele(n-1)%s
    endif

  endif

  ! If there are null_ele elements in the superimpose region then just move them
  ! out of the way to the lord section of the lattice. This prevents unnecessary
  ! splitting.

  i = ix1_split
  do
    i = i + 1
    if (i > lat%n_ele_track) i = 0
    if (lat%ele(i)%key == null_ele$) then
      call new_control (lat, ix)
      lat%ele(ix) = lat%ele(i)  ! copy null_ele
      do ic = lat%ele(i)%ic1_lord, lat%ele(i)%ic2_lord
        j = lat%ic(ic)
        lat%control(j)%ix_slave = ix ! point to new null_ele.
      enddo
      lat%ele(i)%key = -1  ! Mark old null_ele for deletion
      call remove_eles_from_lat (lat)
      i = i - 1
      if (ix2_split > i) ix2_split = ix2_split - 1
    endif
    if (i == ix2_split) exit
  enddo

  ! If element overlays only drifts then just 
  ! insert it in the tracking part of the lat list.

  all_drift = (ix2_split > ix1_split)
  do i = ix1_split+1, ix2_split
    if (lat%ele(i)%key /= drift$) all_drift = .false.
    if (lat%ele(i)%slave_status /= free$) all_drift = .false.
    if (.not. all_drift) exit
  enddo

  if (all_drift) then  
    do i = ix1_split+2, ix2_split    ! remove all drifts but one
      lat%ele(i)%key = -1    ! mark for deletion
    enddo
    call remove_eles_from_lat(lat)    ! And delete
    ix_super = ix1_split + 1
    lat%ele(ix_super) = super_saved
    lat%ele(ix_super)%lord_status  = free$
    lat%ele(ix_super)%slave_status = free$
    ! If a single drift was split give the runt drifts on either end 
    ! Unique names by adding "#1" and "#2" suffixes.
    if (split1_done .and. split2_done) then
      if (lat%ele(ix_super-1)%name == lat%ele(ix_super+1)%name .and. &
                                    lat%ele(ix_super-1)%key == drift$) then
        lat%ele(ix_super-1)%name = trim(lat%ele(ix_super-1)%name) // '#1'
        lat%ele(ix_super+1)%name = trim(lat%ele(ix_super+1)%name) // '#2'
      endif
    endif
    return
  endif

  ! Only possibility left means we have to set up a super_lord element for the
  ! superposition

  ix_super = lat%n_ele_max + 1
  lat%n_ele_max = ix_super
  if (lat%n_ele_max > ubound(lat%ele, 1)) call allocate_lat_ele_array(lat)
  lat%ele(ix_super) = super_saved
  lat%ele(ix_super)%lord_status = super_lord$

  ix_super_con = 0
  length = super_saved%value(l$)

  ix_slave = ix1_split

  ! Go through the list of elements being superimposed upon.
  ! Zero length elements (markers and multipoles) do not get involved here.

  do 

    ix_slave = ix_slave + 1
    if (ix_slave == ix2_split + 1) exit
    if (ix_slave == lat%n_ele_track + 1) ix_slave = 1

    slave => lat%ele(ix_slave)
    slave_saved = slave
    if (slave_saved%value(l$) == 0) cycle

    ! Do we need to set up a super lord to control this slave element?

    if (slave%slave_status == overlay_slave$) then
      setup_lord = .true.
    elseif (slave%slave_status == super_slave$) then
      setup_lord = .false.
    elseif (slave%key == drift$ .and. slave%slave_status /= multipass_slave$) then
      setup_lord = .false.
    else
      setup_lord = .true.
    endif

    ! if yes then create the super lord element

    if (setup_lord) then
      ixn = lat%n_ele_max + 1
      if (ixn > ubound(lat%ele, 1)) call allocate_lat_ele_array(lat)
      lat%ele(ixn) = slave_saved
      lat%ele(ixn)%lord_status = super_lord$
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
      slave%ic1_lord = ic
      slave%ic2_lord = ic + 1
      slave%n_lord = 2
      lat%n_ic_max = ic + 1
      lat%ic(ic) = ixc 

    else
      slave%n_lord = slave_saved%n_lord + 1
      call add_lattice_control_structs (lat, ix_slave)
    endif

    slave%slave_status = super_slave$

    ! add control info for main super lord to list

    ix_super_con = ix_super_con + 1
    sup_con(ix_super_con)%ix_slave = ix_slave
    sup_con(ix_super_con)%ix_lord = ix_super
    sup_con(ix_super_con)%coef = slave_saved%value(l$) / length
    sup_con(ix_super_con)%ix_attrib = 0

    ! change the element key

    call calc_superimpose_key(slave_saved, super_saved, slave)
    if (slave%key <= 0) then
      call out_io (s_abort$, r_name, (/ &
              'ELEMENT: ' // trim(super_saved%name), &
              'OF TYPE: ' // key_name(super_saved%key), &
              'IS TO BE SUPERIMPOSED UPON: ' // trim(slave_saved%name), &
              'OF TYPE: ' // key_name(slave_saved%key), &
              'I DO NOT KNOW HOW TO DO THIS!' /) )
      call err_exit                    
    endif

  enddo

  ! Special case where elements on either side of the superimpose have the same
  ! name

  if (split1_done .and. split2_done .and. &
                lat%ele(ix1_split)%name == lat%ele(ix2_split+1)%name) then
    lat%ele(ix1_split)%name = trim(lat%ele(ix1_split)%name) // '#1'
    lat%ele(ix2_split+1)%name = trim(lat%ele(ix2_split+1)%name) // '#2'
    call delete_underscore (lat%ele(ix1_split))
    call delete_underscore (lat%ele(ix2_split+1))
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
  call adjust_slave_names

!------------------------------------------------------------------------------
contains

! Modify: "#\" -> "\"
!         "##" -> "#"

subroutine delete_underscore(ele)

  implicit none

  type (ele_struct) ele
  integer ix

  !

  ix = index(ele%name, '#\')  ! '
  if (ix /= 0) ele%name = ele%name(1:ix-1) // ele%name(ix+1:)

  ix = index(ele%name, '##')
  if (ix /= 0) ele%name = ele%name(1:ix-1) // ele%name(ix+1:)

end subroutine

!------------------------------------------------------------------------------
! contains

! Adjust the names of the slaves

subroutine adjust_slave_names ()

do i = n_ele_max_old+1, lat%n_ele_max
  lord => lat%ele(i)
  if (lord%lord_status /= super_lord$) cycle
  ix_1lord = 0
  do j = lord%ix1_slave, lord%ix2_slave
    slave => lat%ele(lat%control(j)%ix_slave)
    if (slave%n_lord == 1) then
      ix_1lord = ix_1lord + 1
      write (slave%name, '(2a, i0)') trim(lord%name), '#', ix_1lord
    else
      name = ''
      do k = slave%ic1_lord, slave%ic2_lord
        ix = lat%control(lat%ic(k))%ix_lord
        name = trim(name) //  '\' // lat%ele(ix)%name !'
      enddo
      slave%name = name(2:len(slave%name))
    endif
  enddo
enddo

end subroutine

end subroutine
