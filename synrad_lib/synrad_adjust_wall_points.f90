!+
! Subroutine synrad_adjust_wall_points (wall, branch)
!
! Routine to delete overlapping wall points and add points as needed
! for the triangle construction
!
! Input:
!   wall    -- wall_struct: Wall with points
!   branch  -- branch_struct: lattice branch.
! 
! Output:
!   wall    -- wall_struct: Wall with points deleted and added as needed.
!-

subroutine synrad_adjust_wall_points (wall, branch)

use synrad_struct
use synrad_interface, except => synrad_adjust_wall_points

implicit none

type (wall_struct), target :: wall
type (branch_struct), target :: branch
type (ele_struct), pointer :: ele

real(rp) theta0, dtheta
integer ix_pt, ix1, ix2, k_max, j, k, n_bend
character(*), parameter :: r_name = 'synrad_adjust_wall_points'

! Delete ovlapping wall points

ix_pt = 1
do while (ix_pt <= wall%n_pt_max)
  if (wall%pt(ix_pt)%x == wall%pt(ix_pt-1)%x .and. wall%pt(ix_pt)%s == wall%pt(ix_pt-1)%s) then
    call out_io (s_info$, r_name, &
            'Note: Deleting overlapping point on: ' // wall_name(wall%side), &
            '\i8\ \2f12.6\ ' // trim(wall%pt(ix_pt)%name), &
            r_array = [wall%pt(ix_pt)%s, wall%pt(ix_pt)%x], i_array = [ix_pt])
    wall%pt(ix_pt:wall%n_pt_max-1) = wall%pt(ix_pt+1:wall%n_pt_max)
    wall%n_pt_max = wall%n_pt_max - 1
  else
    ix_pt = ix_pt + 1
  endif
enddo

! Add wall points at bend boundaries to prevent:
!   1) A region between points spanning multiple bends
!   2) A region between points having a bend angle of more than 1 radian.
!      1 radian is rather arbitrary. In theory can handle up to pi radian bends.

ix_pt = 1
pt_loop: do while (ix_pt <= wall%n_pt_max)
  ix1 = element_at_s (branch%lat, wall%pt(ix_pt-1)%s, .true., branch%ix_branch)
  ix2 = element_at_s (branch%lat, wall%pt(ix_pt)%s, .false., branch%ix_branch)
  n_bend = 0
  theta0 = branch%ele(ix1-1)%floor%theta
  do j = ix1, ix2
    ele => branch%ele(j)
    if (ele%key /= sbend$) cycle
    n_bend = n_bend + 1

    if (n_bend == 2) then
      call add_this_pt (branch%ele(j-1)%s, ix_pt)
      cycle pt_loop
    endif

    dtheta = abs(ele%floor%theta - theta0)
    if (dtheta > 1) then
      k_max = int(dtheta)
      do k = 1, k_max
        call add_this_pt (branch%ele(j-1)%s + k * ele%value(l$) / (k_max + 1), ix_pt)
      enddo
      cycle pt_loop
    endif

  enddo

  ix_pt = ix_pt + 1

enddo pt_loop

!------------------------------------------------------------------------------------
contains

subroutine add_this_pt (s, ix_pt)

type (wall_pt_struct), allocatable :: pt_temp(:)
real(rp) s, s0, s1
integer ix_pt, np

!

np = wall%n_pt_max
wall%n_pt_max = np + 1

if (np+1 > ubound(wall%pt, 1)) then
  call move_alloc (wall%pt, pt_temp)
  allocate (wall%pt(0:np + 10))
  wall%pt(0:np) = pt_temp(0:np)
  deallocate (pt_temp)
endif

wall%pt(ix_pt+1:np+1) = wall%pt(ix_pt:np)
wall%pt(ix_pt) = wall%pt(ix_pt+1)
wall%pt(ix_pt)%name = trim(wall%pt(ix_pt)%name) // '#Added'
s0 = wall%pt(ix_pt-1)%s
s1 = wall%pt(ix_pt+1)%s
wall%pt(ix_pt)%x = (wall%pt(ix_pt+1)%x * (s - s0) + wall%pt(ix_pt-1)%x * (s1 - s)) / (s1 - s0)
wall%pt(ix_pt)%s = s
ix_pt = ix_pt + 1

end subroutine

end subroutine
