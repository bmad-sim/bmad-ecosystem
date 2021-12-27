!+
! Subroutine choose_quads_for_set_tune (lat, dk1, mask)
!
! Routine to choose assign weights for quadrupole changes in a lattice when varying the tune.
! The output of this routine, dk1, can be used in the set_tune routine.
!
! Input:
!   lat         -- lat_struct: lattice.
!   mask        -- character(*), optional: If present, assign weight of zero for all quads that
!                   do not match. That is, no variation for non-matching quads.
!
! Output:
!   dk1(:)      -- real(rp): Weights for the quadrupoles. All values will be +1 or -1.
!-

subroutine choose_quads_for_set_tune (lat, dk1, mask)

use bmad_interface, dummy => choose_quads_for_set_tune

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), pointer :: slave
type (control_struct), pointer :: ctl

real(rp), allocatable, intent(inout) :: dk1(:)

integer i, j, is, n_loc

character(*), optional :: mask

logical found

! find which quads to change

if (.not. allocated(dk1)) allocate(dk1(lat%n_ele_max))
dk1 = 0

lat%ele%select = .true.
if (present(mask)) then
  call lat_ele_locator (mask, lat, eles, n_loc)
  do i = 1, n_loc
    eles(i)%ele%select = .false.
  enddo  
endif

do i = 1, lat%n_ele_max
  ele => lat%ele(i)
  if (.not. ele%select) cycle

  if (ele%key == quadrupole$ .and. attribute_free(ele, 'K1',err_print_flag = .false.) .and. &
                                                                      abs(ele%value(tilt$)) < 0.01) then
    if (ele%a%beta > ele%b%beta) then
      dk1(i) = +1
    else
      dk1(i) = -1
    endif
  endif

  if (ele%lord_status == overlay_lord$) then
    found = .false.
    do is = 1, ele%n_slave    
      slave => pointer_to_slave(ele, is, ctl)
      if (ctl%ix_attrib == k1$ .and. slave%key == quadrupole$ .and. slave%value(tilt$) == 0) then
        found = .true.
        exit
      endif
    enddo
    if (.not. found) cycle
    if (slave%a%beta > slave%b%beta) then
      dk1(i) = +1
    else
      dk1(i) = -1
    endif
  endif
enddo

end subroutine choose_quads_for_set_tune
