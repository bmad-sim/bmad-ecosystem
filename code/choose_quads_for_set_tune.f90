!+
! Subroutine choose_quads_for_set_tune (lat, dk1, regex_mask)
!
! Routine to choose assign weights for quadrupole changes in a lattice when varying the tune.
! The output of this routine, dk1, can be used in the set_tune routine.
!
! Input:
!   lat         -- lat_struct: lattice.
!   regex_mask  -- character(*), optional: If present, assign weight of zero for all quads that
!                   do not match this regular expression. That is, no variation for non-matching quads.
!
! Output:
!   dk1(:)      -- real(rp): Weights for the quadrupoles.
!-

subroutine choose_quads_for_set_tune (lat, dk1, regex_mask)

use bmad_interface, dummy => choose_quads_for_set_tune

implicit none

type (lat_struct) lat
type (ele_struct), pointer :: slave
type (control_struct), pointer :: ctl

real(rp), allocatable, intent(inout) :: dk1(:)

integer i, j, is

character(*), optional :: regex_mask
character(40) :: r_mask 

logical found

! find which quads to change

if (.not. allocated(dk1)) allocate(dk1(lat%n_ele_max))

if (present(regex_mask)) then
  call upcase_string(regex_mask)
  r_mask = regex_mask
else
  r_mask = ''
endif

do i = 1, lat%n_ele_max
  if (lat%ele(i)%key == quadrupole$ .and. &
    attribute_free(lat%ele(i), 'K1',err_print_flag = .false.) .and. &
    abs(lat%ele(i)%value(tilt$)) < 0.01) then
    if (.not. match_reg(lat%ele(i)%name, r_mask)) cycle ! if no mask provided, mask set to '', thereby always matching
    if (lat%ele(i)%a%beta > lat%ele(i)%b%beta) then
      dk1(i) = +1
    else
      dk1(i) = -1
    endif
  else
    dk1(i) = 0
  endif

  if (lat%ele(i)%key == match$)then
    dk1(i) = 1 !If there is a match element we will use it to qtune 
    cycle
  endif

  if (lat%ele(i)%lord_status == overlay_lord$) then
    found = .false.
    do is = 1, lat%ele(i)%n_slave    
      slave => pointer_to_slave(lat%ele(i), is, ctl)
      if (ctl%ix_attrib == k1$ .and. slave%key == quadrupole$ .and. slave%value(tilt$) == 0) then
        if (.not. match_reg(slave%name, r_mask)) cycle
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
