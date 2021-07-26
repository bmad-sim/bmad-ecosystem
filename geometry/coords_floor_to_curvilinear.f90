!+
! Function coords_floor_to_curvilinear (floor_coords, ele0, ele1, status, w_mat) result (local_coords)
!
! Given a position in global "floor" coordinates, return local curvilinear (laboratory) coordinates 
! for an appropriate element, ele1, near ele0. That is, the s-position of local_coords will be within
! the longitudinal extent of ele1.
!
! There may not be a corresponding local_coords. This may happen if the lattice is open.
! The status argument will be set accorrdingly.
!
! If there is no corresponding local_coords an "approximate" local_coords will
! be returned that may be used, for example, to do plotting but should not be used in calculations.
!
! Note: Angular orientation of floor_coords is ignored.
!
! Input:
!   floor_coords  -- floor_position_struct: %r = [X, Y, Z] position in global coordinates
!   ele0          -- ele_struct: Element to start the search at.
!
! Output:
!   local_coords  -- floor_position_struct: %r = [x, y, s] position in curvilinear coordinates
!                      with respect to ele1 with s relative to start the lattice branch.
!   ele1          -- ele_struct, pointer: Element that local_coords is with respect to.
!   status        -- logical: ok$             -> Local_coords found.
!                             patch_problem$  -> No solution due to a patch element.
!                             outside$        -> Outside of lattice ends (for open lattices).
!   w_mat(3,3)    -- real(rp) (optional): W matrix at s, to transform vectors from floor to local. 
!                      w_mat will only be well defined if status = ok$
!-  

function coords_floor_to_curvilinear (floor_coords, ele0, ele1, status, w_mat) result (local_coords)

use bmad_interface, dummy => coords_floor_to_curvilinear

implicit none

type (floor_position_struct) floor_coords, local_coords
type (ele_struct), target :: ele0
type (ele_struct), pointer :: ele1
type (branch_struct), pointer :: branch

real(rp), optional :: w_mat(3,3)
real(rp) w_mat0(3,3), w_mat1(3,3), ds_now, ds_old

integer status, this_stat, last_direction, ix_ele, n_try
character(*), parameter :: r_name = 'coords_floor_to_curvilinear'

! Loop over neighboring elements until an encompassing one is found

branch => ele0%branch
ix_ele = ele0%ix_ele
last_direction = 0
n_try = 0
ds_now = 0
status = ok$

do
  n_try = n_try + 1
  ele1 => branch%ele(ix_ele)
  local_coords = coords_floor_to_local_curvilinear (floor_coords, ele1, this_stat) 

  if (ele1%key == patch$) then
    local_coords%r(3) = ele1%value(downstream_coord_dir$) * local_coords%r(3) + ele1%s
  elseif (ele1%orientation == 1) then
    local_coords%r(3) = local_coords%r(3) + ele1%s_start
  else
    local_coords%r(3) = ele1%value(l$) - local_coords%r(3) + ele1%s_start
  endif

  ! No good so look for a problem and keep on searching if warranted

  ds_old = ds_now

  if (this_stat == upstream_end$) then
    ds_now = local_coords%r(3) - ele1%s_start
   
    if (n_try > branch%n_ele_track .or. (ix_ele == 1 .and. branch%param%geometry == open$)) then
      status = outside$
      return
    endif

    if (last_direction == 1) then ! endless loop detected
      call no_convergence_calc (branch, ix_ele-1, ix_ele, ds_old, ds_now, status)
      return
    endif

    ! Try previous element
    ix_ele = ix_ele - 1
    if (ix_ele < 1) ix_ele = branch%n_ele_track
    last_direction = -1
    cycle

  !

  else if (this_stat == downstream_end$) then
    ds_now = local_coords%r(3) - ele1%s

    if (n_try > branch%n_ele_track .or. (ix_ele == branch%n_ele_track .and. branch%param%geometry == open$)) then
      status = outside$
      return
    endif

    if (last_direction == -1) then ! endless loop detected
      call no_convergence_calc (branch, ix_ele, ix_ele+1, ds_now, ds_old, status)
      return
    endif

    ! Try next element
    ix_ele = ix_ele + 1
    if (ix_ele > branch%n_ele_track) ix_ele = 1
    last_direction = 1
    cycle

  ! This element contains position
  else 
    exit
  end if  
enddo

! Optionally return rotation matrix

if (present(w_mat) ) then
  w_mat = matmul(transpose(local_coords%W), floor_coords%W)
endif

!---------------------------------------------------------------------------
contains

subroutine no_convergence_calc (branch, ie1, ie2, ds1, ds2, status)

type (branch_struct), target :: branch
type (ele_struct), pointer :: elem1, elem2
type (floor_position_struct) pos1, pos2

real(rp) ds1, ds2
integer ie1, ie2, status, this_stat

! Check if roundoff errors are throwing off the calculation

elem1 => branch%ele(ie1)
elem2 => branch%ele(ie2)

if (10*abs(ds1) < bmad_com%significant_length .and. 10*abs(ds2) < bmad_com%significant_length) then
  if (abs(ds1) < abs(ds2)) then
    local_coords = coords_floor_to_local_curvilinear (floor_coords, elem1, this_stat)
  else
    local_coords = coords_floor_to_local_curvilinear (floor_coords, elem2, this_stat)
  endif
  status = ok$
  return
endif

!

status = patch_problem$

if (elem1%key /= patch$ .and. elem2%key /= patch$) then
  call out_io (s_fatal$, r_name, 'CANNOT FIND CORRESPONDING LOCAL POSITION')
  if (global_com%exit_on_error) call err_exit
  return
endif

if (elem1%key == patch$ .and. ele1%orientation == 1) then
  pos1 = coords_floor_to_local_curvilinear (floor_coords, branch%ele(ie1-1), this_stat)
  pos1%r(3) = floor_coords%r(3) - elem1%value(z_offset$)
else
  pos1 = coords_floor_to_local_curvilinear (floor_coords, elem1, this_stat)
  pos1%r(3) = pos1%r(3) - elem1%value(l$)
endif

pos2 = coords_floor_to_local_curvilinear (floor_coords, elem2, this_stat)

if (abs(pos1%r(3)) < abs(pos2%r(3))) then
  local_coords = pos1
  ele1 => elem1
else
  local_coords = pos2
  ele1 => elem2
endif

local_coords%r(3) = elem1%s

end subroutine no_convergence_calc

end function coords_floor_to_curvilinear
