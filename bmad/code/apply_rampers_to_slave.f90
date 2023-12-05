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
type (ele_struct), pointer :: rmp
type (ele_pointer_struct), target :: ramper(:)
type (lat_struct), pointer :: lat
type (control_ramp1_struct), pointer :: r1

integer iv, key, ix, ir
logical err_flag, ok

character(100) err_str
character(40) name
character(*), parameter :: r_name = 'apply_rampers_to_slave'

! Init

err_flag = .false.
if (size(ramper) == 0) return
lat => ramper(1)%ele%branch%lat

do ix = 1, size(ramper)
  ramper(ix)%ele%select = .false.
enddo

! Bookkeeping for ramper controlling ramper.

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

      if ((key /= 0 .and. key /= slave%key) .or. .not. match_wild(slave%name, name)) then
        r1%value = real_garbage$  ! This ramper does not control this slave.
        cycle
      endif
    endif

    call this_slave_bookkeeper(rmp, slave, r1)
  enddo
enddo

if (slave%key /= int_garbage$) call attribute_bookkeeper(slave, .true.)

!-------------------------------------------------------------------
contains

recursive subroutine this_ramper_bookkeeper (this_ramp, ramper, lat)

type (ele_struct) this_ramp
type (ele_pointer_struct) ramper(:)
type (ele_struct), pointer :: lord, slave
type (lat_struct) lat
type (control_ramp1_struct), pointer :: r1

integer ir, is

!

if (this_ramp%select) return

do ir = 1, this_ramp%n_lord
  lord => pointer_to_lord(this_ramp, ir)
  if (lord%key /= ramper$) cycle
  call this_ramper_bookkeeper(this_ramp, ramper, lat)
enddo

if (this_ramp%is_on) then
  do is = 1, size(this_ramp%control%ramp)
    r1 => this_ramp%control%ramp(is)
    if (index(r1%slave_name, '*') /= 0 .or. index(r1%slave_name,'%') /= 0) cycle
    do ir = 1, size(ramper)
      if (ramper(ir)%ele%name /= r1%slave_name) cycle
      slave => pointer_to_ele(lat, r1%slave)
      call this_slave_bookkeeper(ramper(ir)%ele, slave, r1)
    enddo
  enddo
endif

this_ramp%select = .true.

end subroutine this_ramper_bookkeeper

!-------------------------------------------------------------------
! contains

subroutine this_slave_bookkeeper (this_ramp, slave, r1)

type (ele_struct) this_ramp, slave
type (control_ramp1_struct) r1
type (all_pointer_struct) a_ptr

logical err_flag

! slave%key = int_garbage$ is used by the controller_function_plot program to bypass 
! some of the bookkeeping of this routine.

if (slave%key /= int_garbage$) then
  call pointer_to_attribute (slave, r1%attribute, .true., a_ptr, err_flag, .false.)
  if (err_flag .or. .not. associated(a_ptr%r)) then
    r1%value = real_garbage$
    return
  endif
endif

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

if (slave%key /= int_garbage$) then
  a_ptr%r = r1%value
  call set_flags_for_changed_attribute (slave, a_ptr%r, .true.)
endif

end subroutine this_slave_bookkeeper

end subroutine
