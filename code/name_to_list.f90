!+
! Subroutine name_to_list (lat, ele_names)
!
! Routine to mark the elements in lat whose name matches the names in ELE_NAMES.
! This routine is typiclly used with make_hybrid_lat.
!
! Input:
!   lat         -- lat_struct: Input lat.
!   ele_names(:) -- Character(*): list of element names. Wild card
!                     characters may be used.
!
! Output:
!   lat%branch(:)%ele(:)%select  -- Set True if there is a match. False otherwise
!
! Example: The following makes a list of element whose name begin with "Q" or "B"
!
!   ele_names(1) = 'Q*'  
!   ele_names(2) = 'B*'
!-

subroutine name_to_list (lat, ele_names)

use bmad_interface, except_dummy => name_to_list

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele

integer ib, n, m

character(*) ele_names(:)

logical searching

! match

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do n = 1, branch%n_ele_max
    ele => branch%ele(n)
    ele%select = .false.

    do m = 1, size(ele_names)
      if (.not. match_wild(ele%name, ele_names(m))) cycle
      ele%select = .true.
      exit
    enddo

  enddo
enddo

end subroutine
