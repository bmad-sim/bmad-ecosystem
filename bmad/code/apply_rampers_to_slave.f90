!+
! Subroutine apply_rampers_to_slave (slave, ramper, err_flag)
!
! Routine to apply the ramper elements of a lattice to a particular element.
! Also see: apply_all_rampers.
!
! Input:
!   slave       -- ele_struct: Element to apply ramper elements to.
!   ramper(:)   -- ele_pointer_struct: Pointers to ramper elements in the lattice to use.
!
! Output:
!   err_flag    -- logical: Set true if there is an error. False otherwise.
!-

subroutine apply_rampers_to_slave (slave, ramper, err_flag)

use bmad, except_dummy => apply_rampers_to_slave

implicit none

type (ele_struct), target :: slave
type (ele_struct), pointer :: rmp, slave2
type (ele_pointer_struct), target :: ramper(:)
type (lat_struct), pointer :: lat
type (control_ramp1_struct), pointer :: r1
type (ele_pointer_struct), allocatable :: slave_list(:)

integer iv, key, ix, ir, n, is
logical err_flag, ok, found

character(100) err_str
character(40) name
character(*), parameter :: r_name = 'apply_rampers_to_slave'

! Init

err_flag = .false.
if (size(ramper) == 0) return
lat => ramper(1)%ele%branch%lat

! Bookkeeping for ramper controlling ramper.

do ix = 1, size(ramper)
  ramper(ix)%ele%select = .false.   ! Bookkeeping has not yet been done.
enddo

do ix = 1, size(ramper)
  call this_ramper_bookkeeper(ramper(ix)%ele, ramper, lat)
enddo

! Calculate slave values

do ir = 1, size(ramper)
  rmp => ramper(ir)%ele
  if (.not. rmp%is_on) cycle

  do iv = 1, size(rmp%control%ramp)
    r1 => rmp%control%ramp(iv)

    ! slave%key = int_garbage$ is used by the controller_function_plot program to bypass 
    ! some of the bookkeeping of this routine.

    if (slave%key /= int_garbage$) then
      ix = index(r1%slave_name, '::')
      if (ix == 0) then
        key = 0
        name = r1%slave_name
      else
        key = key_name_to_key_index(r1%slave_name(1:ix-1), .true.)
        name = r1%slave_name(ix+2:)
      endif

      if (r1%is_controller) then
        slave2 => pointer_to_ele(lat, r1%slave)
        ! In case elements have shifted in the lattice, check that slave2 is the correct element
        if (.not. associated(slave2)) slave2 => lat%ele(0)   ! Just to point to something
        if (slave2%name /= r1%slave_name) then
          slave2 => pointer_to_ele(lat, r1%slave_name)
          r1%slave = ele_loc(slave2)
        endif
        !
        call get_slave_list(slave2, slave_list, n)
        found = .false.
        do is = 1, n
          if (slave_list(is)%ele%name /= slave%name) cycle
          found = .true.
          exit
        enddo
        if (.not. found) then
          r1%value = real_garbage$  ! This ramper does not control this slave.
          cycle
        endif

      elseif ((key /= 0 .and. key /= slave%key) .or. .not. match_wild(slave%name, name)) then
        r1%value = real_garbage$  ! This ramper does not control this slave.
        cycle
      endif
    endif

    call this_slave_bookkeeper(rmp, slave, r1)
  enddo
enddo

if (slave%key == int_garbage$) return
call attribute_bookkeeper(slave, .true.)

!-------------------------------------------------------------------
contains

recursive subroutine this_ramper_bookkeeper (this_ramp, ramper, lat)

type (ele_struct) this_ramp
type (ele_pointer_struct) ramper(:)
type (ele_struct), pointer :: lord, slave
type (lat_struct) lat
type (control_ramp1_struct), pointer :: r1

integer ir, is

! Nothing to do if bookkeeping has been done.

if (this_ramp%select) return 

! Ramper lord bookkeeping

do ir = 1, this_ramp%n_lord
  lord => pointer_to_lord(this_ramp, ir)
  if (lord%key /= ramper$) cycle
  call this_ramper_bookkeeper(this_ramp, ramper, lat)
enddo

! Ramper slave bookkeeping

if (this_ramp%is_on) then
  do is = 1, size(this_ramp%control%ramp)
    r1 => this_ramp%control%ramp(is)
    ! No bookkeeping needed if the slave name has wild card characters.
    if (index(r1%slave_name, '*') /= 0 .or. index(r1%slave_name,'%') /= 0) cycle
    do ir = 1, size(ramper)
      if (ramper(ir)%ele%name /= r1%slave_name) cycle
      slave => pointer_to_ele(lat, r1%slave)
      call this_slave_bookkeeper(ramper(ir)%ele, slave, r1)
    enddo
  enddo
endif

! And bookkeeping has been done

this_ramp%select = .true.   

end subroutine this_ramper_bookkeeper

!-------------------------------------------------------------------
! contains

subroutine this_slave_bookkeeper (this_ramp, slave, r1)

type (ele_struct), target :: this_ramp, slave
type (ele_struct), pointer :: slave2
type (control_ramp1_struct) r1
type (all_pointer_struct) a_ptr

logical err_flag

!

if (allocated(r1%stack)) then
  r1%value = expression_stack_value(r1%stack, err_flag, err_str, this_ramp%control%var, .false.)
  if (err_flag) then
    call out_io (s_error$, r_name, err_str, ' OF RAMPER: ' // this_ramp%name)
    err_flag = .true.
    return
  endif

elseif (allocated(r1%y_knot)) then
  r1%value = knot_interpolate(this_ramp%control%x_knot, r1%y_knot, &
          this_ramp%control%var(1)%value, nint(this_ramp%value(interpolation$)), err_flag)
  if (err_flag) then
    call out_io (s_error$, r_name, 'VARIABLE VALUE (\es12.4\) OF RAMPER ELEMENT: ' // this_ramp%name, &
                                   'IS OUTSIDE OF SPLINE KNOT RANGE OF SLAVE: ' // r1%slave_name)
    return
  endif
endif

! slave%key = int_garbage$ is used by the controller_function_plot program to bypass 
! some of the bookkeeping of this routine.

if (slave%key == int_garbage$) return

if (r1%is_controller) then
  slave2 => pointer_to_ele(lat, r1%slave)
else
  slave2 => slave  
endif

call pointer_to_attribute (slave2, r1%attribute, .true., a_ptr, err_flag, .false.)
if (err_flag .or. .not. associated(a_ptr%r)) then
  r1%value = real_garbage$
  return
endif

a_ptr%r = r1%value
call set_flags_for_changed_attribute (slave2, a_ptr%r, .true.)
if (r1%is_controller) call control_bookkeeper(lat, slave2)

end subroutine this_slave_bookkeeper

end subroutine
