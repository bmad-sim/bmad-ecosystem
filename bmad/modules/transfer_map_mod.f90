module transfer_map_mod

use element_at_s_mod
use taylor_mod

private transfer_this_map, transfer_this_mat

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+         
! Subroutine transfer_map_from_s_to_s (lat, t_map, s1, s2, ref_orb_in, ref_orb_out, ix_branch, 
!                                          one_turn, unit_start, err_flag, concat_if_possible, spin_map)
!
! Subroutine to calculate the transfer map between longitudinal positions s1 to s2.
!
! If s2 < s1 and lat%param%geometry is closed$ then the
! calculation will 'wrap around' the lattice end.
! For example, if s1 = 900 and s2 = 10 then the t_map is the map from
! element 900 to the lattice end plus from 0 through 10.
!
! If s2 < s1 and lat%param%geometry is open$ then the inverse of the forward map of s2 -> s1 is computed.
!
! If s2 = s1 then you get the unit map except if one_turn = True and the lattice is circular.
!
! Input:
!   lat         -- lat_struct: Lattice used in the calculation.
!   t_map(6)    -- taylor_struct: Initial map (used when unit_start = False)
!   s1          -- real(rp), optional: Element start position for the calculation.
!                    Default is 0.
!   s2          -- real(rp), optional: Element end position for the calculation.
!                    Default is lat%param%total_length.
!   ref_orb_in  -- coord_struct, optional: Reference orbit/particle at s1 around which the map is made.
!                    This arg is needed if: unit_start = True or particle is not the same as the reference 
!                    particle of the lattice.
!   ix_branch   -- integer, optional: Lattice branch index. Default is 0 (main branch).
!   one_turn    -- logical, optional: If present and True, and s1 = s2, and the lattice
!                    is circular: Construct the one-turn map from s1 back to s1.
!                    Otherwise t_map is unchanged or the unit map if unit_start = T.
!                    Default = False.
!   unit_start  -- logical, optional: If present and False then t_map will be
!                    used as the starting map instead of the unit map.
!                    Default = True
!   concat_if_possible
!               -- logical, optional: If present and True then use map concatenation rather than tracking 
!                    if a map is present for a given lattice element. See above. Default is False.
!   spin_map(4) -- taylor_struct, optional: Initial spin map.
!
! Output:
!   t_map(6)    -- taylor_struct: Transfer map.
!   ref_orb_out -- coord_struct, optional: Ending coordinates of the reference orbit.
!                    This is also the actual orbit of particle 
!   err_flag    -- logical, optional: Set true if there is an error. False otherwise.
!   spin_map(4) -- taylor_struct, optional: Final spin map. Only computed if bmad_com%spin_tracking_on = T.
!-

subroutine transfer_map_from_s_to_s (lat, t_map, s1, s2, ref_orb_in, ref_orb_out, ix_branch, &
                                            one_turn, unit_start, err_flag, concat_if_possible, spin_map)

use ptc_interface_mod, only: concat_taylor, taylor_inverse

implicit none

type (lat_struct), target :: lat
type (coord_struct), optional :: ref_orb_in, ref_orb_out
type (coord_struct) orb
type (taylor_struct) :: t_map(:)
type (taylor_struct) a_map(6), s_map(4)
type (taylor_struct), optional :: spin_map(:)
type (branch_struct), pointer :: branch

real(rp), intent(in), optional :: s1, s2
real(rp) ss1, ss2, v6(6)

integer, optional :: ix_branch
integer ix_br, ix_ele

logical, optional :: one_turn, unit_start, err_flag, concat_if_possible
logical unit_start_this, error_flag, do_spin

character(*), parameter :: r_name = 'transfer_map_from_s_to_s'

!

ix_br = integer_option(0, ix_branch)
branch => lat%branch(ix_br)
do_spin = present(spin_map) .and. bmad_com%spin_tracking_on

unit_start_this = logic_option(.true., unit_start)
if (present(err_flag)) err_flag = .true.

call check_if_s_in_bounds (branch, real_option(0.0_rp, s1), error_flag, ss1)
if (error_flag) return

call check_if_s_in_bounds (branch, real_option(branch%param%total_length, s2), error_flag, ss2)
if (error_flag) return

ix_ele = element_at_s (branch, ss1, .true.)

if (unit_start_this) then
  call taylor_make_unit (t_map)
  if (present(ref_orb_in)) t_map%ref = ref_orb_in%vec
  if (do_spin) call taylor_make_unit(spin_map)
endif

if (present(ref_orb_in)) then
  call init_coord(orb, ref_orb_in, branch%ele(ix_ele), inside$)
else
  v6 = 0
  call init_coord(orb, v6, branch%ele(ix_ele), inside$)
endif

! One turn or not calc?

if (ss1 == ss2 .and. (.not. logic_option (.false., one_turn) .or. branch%param%geometry == open$)) then
  if (present(err_flag)) err_flag = .false.
  return
endif

! Normal case

if (ss1 < ss2) then 
  call transfer_this_map (t_map, branch, ss1, ss2, error_flag, orb, do_spin, concat_if_possible, spin_map)
  if (error_flag) return

! For a circular lattice push through the origin.

elseif (branch%param%geometry == closed$) then
  call transfer_this_map (t_map, branch, ss1, branch%param%total_length, error_flag, orb, do_spin, concat_if_possible, spin_map)
  if (error_flag) return
  call transfer_this_map (t_map, branch, 0.0_rp, ss2, error_flag, orb, do_spin, concat_if_possible, spin_map)
  if (error_flag) return

! For an open lattice compute the backwards map

else
  if (unit_start_this) then
    call transfer_this_map (t_map, branch, ss2, ss1, error_flag, orb, do_spin, concat_if_possible, spin_map)
    if (error_flag) return
    call taylor_inverse (t_map, t_map)
    if (do_spin) call taylor_inverse(spin_map, spin_map)
  else  
    call taylor_make_unit (a_map)
    if (do_spin) call taylor_make_unit(s_map)
    call transfer_this_map (a_map, branch, ss2, ss1, error_flag, orb, do_spin, concat_if_possible, s_map)
    if (error_flag) return
    call taylor_inverse (a_map, a_map)
    call concat_taylor (t_map, a_map, t_map)
    call kill_taylor (a_map)
    if (do_spin) then
      call taylor_inverse (s_map, s_map)
      call concat_taylor (spin_map, s_map, spin_map)
      call kill_taylor (s_map)
    endif
  endif

endif

if (present(ref_orb_out)) ref_orb_out = orb
if (present(err_flag)) err_flag = .false.

end subroutine transfer_map_from_s_to_s

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine transfer_this_map (t_map, branch, s_1, s_2, error_flag, orb, do_spin, concat_if_possible, spin_map)
!
! Private subroutine used by transfer_map_from_s_to_s
!-

subroutine transfer_this_map (t_map, branch, s_1, s_2, error_flag, orb, do_spin, concat_if_possible, spin_map)

use ptc_interface_mod, only: taylor_propagate1, concat_ele_taylor

implicit none

type (branch_struct), target :: branch
type (taylor_struct) :: t_map(:)
type (taylor_struct), optional :: spin_map(:)
type (ele_struct), pointer :: ele
type (ele_struct), pointer :: runt => null()
type (ele_struct), target, save :: runt_save
type (ele_struct), target :: runt_nosave
type (coord_struct) :: orb
type (coord_struct) orb2

real(rp) s_1, s_2, s_now, s_end, ds

integer i, ix_ele

logical do_spin, create_it, track_upstream_end, track_downstream_end, track_entire_ele
logical runt_points_to_new, error_flag, include_next_ele
logical, optional :: concat_if_possible
logical, save :: old_track_end = .false.

! Init
! Want to get the whole lattice if [s_1, s_2] spans the entire lattice

error_flag = .false.

if (s_1 == branch%ele(0)%s) then
  ix_ele = 1
else
  ix_ele = element_at_s (branch, s_1, .true.)
endif

s_now = s_1

! Loop over all the element to track through.

do
  ele => branch%ele(ix_ele)
  s_end = min(s_2, ele%s)

  track_upstream_end   = (s_now == branch%ele(ix_ele)%s_start) 
  track_downstream_end = (s_end == ele%s)
  track_entire_ele     = (track_upstream_end .and. track_downstream_end)

  ds = s_end - s_now

  ! runt points at ele if we are tracking through the entire element and
  ! at runt_save if only a partial track. We do this since we will mangle
  ! the element with a partial track.

  runt_points_to_new = .false.

  if (track_entire_ele) then
    runt => ele
  elseif (.not. global_com%mp_threading_is_safe) then
    runt => runt_nosave
    runt = ele
  else if (.not. associated(runt, ele) .or. .not. associated(runt, runt_save)) then ! partial track
    call transfer_ele (ele, runt_save, .true.)
    runt => runt_save
    runt_points_to_new = .true.
  endif

  ! We only need to do the "split" bookkeeping if we are only partially tracking
  ! through the element and we have not done the bookkeeping before.

  if (.not. track_entire_ele) then
    create_it = .false.

    if (.not. global_com%mp_threading_is_safe .or. ds /= runt%value(l$) .or. runt_points_to_new) then
      create_it = .true.
    elseif (ele%key == sbend$ .or. ele%key == rf_bend$) then
      if (track_upstream_end .or. track_downstream_end .or. old_track_end) create_it = .true.
    elseif (.not. ele_has_constant_ds_dt_ref(ele)) then
      create_it = .true.
    endif

    if (create_it) then
      call create_element_slice (runt, ele, ds, s_now-branch%ele(ix_ele)%s_start, &
                                     branch%param, track_upstream_end, track_downstream_end, error_flag)
      if (error_flag) exit
    endif
  endif

  ! Now for the integration step

  if (track_entire_ele .and. logic_option(.false., concat_if_possible)  .and. associated(ele%taylor(1)%term)) then
    call concat_ele_taylor (t_map, ele, error_flag, spin_map); if (error_flag) return
    call init_coord(orb, t_map%ref, ele, downstream_end$)
  else
    call taylor_propagate1 (t_map, runt, branch%param, error_flag, orb, spin_map); if (error_flag) return
    call init_coord(orb, t_map%ref, runt, downstream_end$)
  endif

  ! Save the present integration step parameters so that if this routine
  ! is called in the future we can tell if the saved taylor t_map is still valid.

  old_track_end = track_upstream_end .or. track_downstream_end

  ! Are we done?
  ! Include any zero length elements at end of region.

  include_next_ele = .false.
  if (ix_ele + 1 <= branch%n_ele_track) include_next_ele = (branch%ele(ix_ele+1)%s <= s_end)
  if (.not. include_next_ele .and. abs(s_end - s_2) < bmad_com%significant_length) exit

  ! We are not done so move to the next element.

  s_now = s_end
  ix_ele = ix_ele + 1

enddo

! Cleanup

if (.not. global_com%mp_threading_is_safe) then
  call deallocate_ele_pointers (runt_nosave)
endif

end subroutine transfer_this_map

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+         
! Subroutine mat6_from_s_to_s (lat, mat6, vec0, s1, s2, ref_orb_in, ref_orb_out, 
!                                                      ix_branch, one_turn, unit_start, err_flag, ele_save)
!
! Subroutine to calculate the transfer map between longitudinal positions
! s1 to s2.
!
! If s2 < s1 and lat%param%geometry is closed$ then the
! calculation will "wrap around" the lattice end.
! For example, if s1 = 900 and s2 = 10 then the xfer_mat is the matrix from
! element 900 to the lattice end plus from 0 through 10.
!
! If s2 < s1 and lat%param%geometry is open$ then the inverse matrix to the s2 -> s1 matrix is computed.
!
! If s2 = s1 then you get the unit matrix except if one_turn = True and the lattice is circular.
!
! Important Note: To save time, it is assumed that the transfer matrices stored in the
! elements (ele%mat6) are consistant with the reference orbit given by the "orbit" argument.
! If this is not true, lat_make_mat6 must be called before hand.
!
! Input:
!   lat         -- Lat_struct: Lattice used in the calculation.
!   mat6(6,6)   -- Real(rp): Initial matrix (used when unit_start = False)
!   vec0(6)     -- Real(rp): Initial 0th order map (used when unit_start = False)
!   s1          -- Real(rp), optional: Element start index for the calculation.
!                    Default is 0.
!   s2          -- Real(rp), optional: Element end index for the calculation.
!                    Default is lat%param%total_length.
!   ref_orb_in  -- coord_struct, optional: Starting coordinates of the reference orbit about 
!                    which mat6 is made. Default is to make it around the zero-orbit.
!   ix_branch   -- Integer, optional: Lattice branch index. Default is 0 (main branch).
!   one_turn    -- Logical, optional: If present and True, and s1 = s2, and the lattice
!                    is circular: Construct the one-turn matrix from s1 back to s1.
!                    Otherwise mat6 is unchanged or the unit map if unit_start = T.
!                    Default = False.
!   unit_start  -- Logical, optional: If present and False then mat6 will be
!                    used as the starting matrix instead of the unit matrix.
!                    Default = True
!   ele_save    -- ele_struct, optional: Used as scratch space to save values between calls.
!                    Only useful if the next call has s1 equal to the present s2. 
!
! Output:
!   mat6(6,6)   -- Real(rp): Transfer matrix.
!   vec0(6)     -- Real(rp): 0th order part of the map.
!   ref_orb_out -- coord_struct, optional: Ending coordinates of the reference orbit.
!                    This is also the actual orbit of particle 
!   err_flag    -- Logical, optional: Set True if there is an error. False otherwise.
!   ele_save    -- ele_struct, optional: Used as scratch space to save values between calls.
!                    Only useful if the next call has s1 equal to the present s2. 
!-

subroutine mat6_from_s_to_s (lat, mat6, vec0, s1, s2, ref_orb_in, ref_orb_out, &
                                                     ix_branch, one_turn, unit_start, err_flag, ele_save)

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (coord_struct), optional :: ref_orb_in, ref_orb_out
type (coord_struct) orb
type (ele_struct), optional, target :: ele_save

real(rp) mat6(:,:), vec0(:)
real(rp), intent(in), optional :: s1, s2
real(rp) ss1, ss2, v6(6)

integer, optional :: ix_branch
integer ix_ele

logical, optional :: one_turn, unit_start, err_flag
logical error_flag

!

if (present(err_flag)) err_flag = .true.
branch => lat%branch(integer_option(0, ix_branch))

call check_if_s_in_bounds (branch, real_option(0.0_rp, s1), error_flag, ss1)
if (error_flag) return

call check_if_s_in_bounds (branch, real_option(branch%param%total_length, s2), error_flag, ss2)
if (error_flag) return

if (logic_option(.true., unit_start)) then
  call mat_make_unit (mat6)
  vec0 = 0
endif

if (present(ref_orb_in)) then
  orb = ref_orb_in
else
  ix_ele = element_at_s (branch, s1, .true.)
  v6 = 0
  call init_coord(orb, v6, branch%ele(ix_ele), inside$)
endif

! One turn or not calc?

if (ss1 == ss2 .and. (.not. logic_option (.false., one_turn) .or. branch%param%geometry == open$)) then
  if (present(err_flag)) err_flag = .false.
  return
endif

! Normal case

if (ss1 < ss2) then
  call transfer_this_mat (mat6, vec0, branch, ss1,  ss2, error_flag, orb, ele_save)
  if (error_flag) return

! For a circular lattice push through the origin.

elseif (branch%param%geometry == closed$) then
  call transfer_this_mat (mat6, vec0, branch, ss1,  branch%param%total_length, error_flag, orb)
  if (error_flag) return
  call transfer_this_mat (mat6, vec0, branch, 0.0_rp, ss2, error_flag, orb, ele_save)
  if (error_flag) return

! For an open lattice compute the backwards matrix

else
  call transfer_this_mat (mat6, vec0, branch, ss2, ss1, error_flag, orb, ele_save)
  if (error_flag) return
  call mat_inverse (mat6, mat6)
  vec0 = -matmul(mat6, vec0)
endif

if (present(err_flag)) err_flag = .false.
if (present(ref_orb_out)) ref_orb_out = orb

end subroutine mat6_from_s_to_s

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine transfer_this_mat (mat6, vec0, branch, s_1, s_2, error_flag, orb, ele_save)
!
! Private subroutine used by mat6_from_s_to_s
!-

subroutine transfer_this_mat (mat6, vec0, branch, s_1, s_2, error_flag, orb, ele_save)

implicit none

type (branch_struct), target :: branch
type (ele_struct), pointer :: ele
type (ele_struct), pointer :: runt
type (ele_struct), target :: runt_nosave
type (ele_struct), optional, target :: ele_save
type (coord_struct) :: orb, orb2

real(rp) mat6(:,:), vec0(:)
real(rp) s_1, s_2, s_end, s_now, ds

integer ix_ele

logical track_upstream_end, track_downstream_end, track_entire_ele
logical use_saved, error_flag, include_next_ele

! Init

if (s_1 == branch%ele(0)%s) then
  ix_ele = 1
else
  ix_ele = element_at_s (branch, s_1, .true.)
endif

s_now = s_1

! Loop over all the element to track through.

do
  ele => branch%ele(ix_ele)
  s_end = min(s_2, ele%s)

  track_upstream_end   = (s_now == branch%ele(ix_ele)%s_start) 
  track_downstream_end       = (s_end == ele%s)
  track_entire_ele = (track_upstream_end .and. track_downstream_end)

  ds = s_end - s_now

  ! runt points at ele if we are tracking through the entire element and
  ! at runt_save if only a partial track. We do this since we will mangle
  ! the element with a partial track.

  use_saved = .false.
  if (present(ele_save)) then
    if (ele_save%s == s_now .and. ele_save%ix_ele == ele%ix_ele) use_saved = .true.
  endif

  if (track_entire_ele) then
    runt => ele
    call track1 (orb, runt, branch%param, orb2)

  elseif (use_saved) then
    runt => ele_save
    call transfer_ele (ele, runt, .true.)
    call create_element_slice (runt, ele, ds, s_now-branch%ele(ix_ele)%s_start, &
                           branch%param, track_upstream_end, track_downstream_end, error_flag, runt)
    if (error_flag) exit
    call make_mat6 (runt, branch%param, orb, orb2)

  else
    runt => runt_nosave
    call transfer_ele (ele, runt, .true.)
    call create_element_slice (runt, ele, ds, s_now-branch%ele(ix_ele)%s_start, &
                                      branch%param, track_upstream_end, track_downstream_end, error_flag)
    if (error_flag) exit
    call make_mat6 (runt, branch%param, orb, orb2)
  endif

  ! Now for the integration step

  mat6 = matmul (runt%mat6, mat6)
  vec0 = matmul (runt%mat6, vec0) + runt%vec0

  ! Are we done?
  ! Include any zero length elements at end of region.

  orb = orb2
  include_next_ele = .false.
  if (ix_ele + 1 <= branch%n_ele_track) include_next_ele = (branch%ele(ix_ele+1)%s <= s_end)
  if (.not. include_next_ele .and. abs(s_end - s_2) < bmad_com%significant_length) exit

  ! We are not done so move to the next element.

  s_now = s_end
  ix_ele = ix_ele + 1
enddo

! Cleanup

if (.not. global_com%mp_threading_is_safe) then
  call deallocate_ele_pointers (runt_nosave)
endif

end subroutine transfer_this_mat

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+         
! Subroutine concat_transfer_mat (mat_1, vec_1, mat_0, vec_0, mat_out, vec_out)
!
! Routine to concatinate two linear maps:
!   mat_out = matmul(mat_1, mat_0)
!   vec_out = matmul(mat_1, vec_0) + vec_1
!
! Input:
!   mat_1(6,6), vec_1(6) -- Real(rp): Map from s1 to s2
!   mat_0(6,6), vec_0(6) -- Real(rp): Map from s0 to s1
!
! Output:
!   mat_out(6,6), vec_out(6) -- Real(rp): Map from s0 to s2
!-

subroutine concat_transfer_mat (mat_1, vec_1, mat_0, vec_0, mat_out, vec_out)

implicit none

real(rp) mat_1(6,6), vec_1(6), mat_0(6,6), vec_0(6), mat_out(6,6), vec_out(6)
real(rp) mat6(6,6), vec0(6)

!

mat6 = matmul(mat_1, mat_0)
mat_out = mat6

vec0 = matmul(mat_1, vec_0) + vec_1
vec_out = vec0

end subroutine concat_transfer_mat 

end module
