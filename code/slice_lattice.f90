!+
! Subroutine slice_lattice (lat, orbit, ele_list, error)
!
! Routine to discard from the lattice all elements not in ele_list.
! Note controllers that control elements that remain will not be cut.
!
! For each branch where there are elements to be deleted and where the reference energy has been computed:
!   1) The Twiss and reference energy parameters from the first non-deleted element are 
!       transferred to the beginning element.
!   2) The beginning betatron phase is set to zero.
!   3) The branch geometry is set to open.
!   
!
! Note: Not modified is:
!   1) Beginning s (longitudinal position) value.
!   2) Beginning floor position.
!
! Input:
!   lat           -- lat_struct: Lattice to slice.
!   orbit(:)      -- coord_struct: Orbit to transfer to lat%particle_start.
!   ele_list      -- character(*): List of elements to retain. See the documentation for
!                     the lat_ele_locator routine for the syntax of the list.
!
! Output:
!   lat           -- lat_struct: Lattice with unwanted elements sliced out.
!   error         -- logical: Set True if there is an error Set False if not.
!-

subroutine slice_lattice (lat, orbit, ele_list, error)

use bmad, dummy => slice_lattice

implicit none

type (lat_struct), target :: lat
type (coord_struct) :: orbit(0:)
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele0, ele1
type (ele_pointer_struct), allocatable :: eles(:)

integer ie, ib, n_loc
logical error, err

character(*) ele_list
character(*), parameter :: r_name = 'slice_lattice'

!

error = .true.

call lat_ele_locator (ele_list, lat, eles, n_loc, err, above_ubound_is_err = .false.)
if (err) return
if (n_loc == 0) then
  call out_io (s_error$, r_name, 'NO LATTICE ELEMENTS FOUND: ' // ele_list, 'LATTICE NOT SLICED.')
  return
endif

! Use ele%ixx = -1 to tag elements to be deleted which is everything not in eles list.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  branch%ele(1:)%ixx = -1
enddo

do ie = 1, n_loc
  eles(ie)%ele%ixx = 0  ! Do not delete
enddo

! Now go through and save all controllers that can be saved

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do ie = 1, branch%n_ele_max
    if (branch%ele(ie)%ixx == -1) cycle
    call add_back_controllers(branch%ele(ie))
  enddo
enddo

! Transfer Twiss from first non-deleted element back to beginning element.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do ie = 1, branch%n_ele_track
    if (branch%ele(ie)%ixx == -1) cycle
    ele0 => branch%ele(0)
    ele1 => branch%ele(ie-1)
    if (ele1%value(e_tot$) <= 0) exit  ! Energy has not been computed
    call transfer_twiss (ele1, ele0)
    ele0%value(e_tot$)       = ele1%value(e_tot$)
    ele0%value(e_tot_start$) = ele0%value(e_tot$)
    ele0%value(p0c$)         = ele1%value(p0c$)
    ele0%value(p0c_start$)   = ele0%value(p0c$)
    ele0%a%phi = 0
    ele0%b%phi = 0
    ele0%z%phi = 0
    call set_flags_for_changed_attribute(ele0, ele1%value(p0c$))
    branch%param%geometry = open$
    if (ib == 0) lat%particle_start = orbit(ie-1)
    exit
  enddo
enddo

! And remove

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do ie = 1, branch%n_ele_max
    if (branch%ele(ie)%ixx == -1) branch%ele(ie)%key = -1
  enddo
enddo

call remove_eles_from_lat (lat)
call lattice_bookkeeper (lat)

error = .false.

!-------------------------------------------------------------------------------------
contains

recursive subroutine add_back_controllers (ele)

type (ele_struct) ele
type (ele_struct), pointer :: ele2
integer i

do i = 1, ele%n_lord
  ele2 => pointer_to_lord (ele, i)
  ele2%ixx = 0
  call add_back_controllers(ele2)
enddo

end subroutine add_back_controllers

end subroutine slice_lattice
