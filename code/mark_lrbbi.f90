!+
! subroutine MARK_LRBBI(master_lat, master_lat_oppos, lat, crossings)
!
! Subroutine to insert markers at parasitic crossing point.
!
! Modules Needed:
!   use bmad
!
! Input:
!   master_lat        -- Lat struct: Lat with markers at LRBBI locations.
!   master_lat_oppos  -- Lat struct: Lat for oppositely circulating
!                         particles, with markers at all LRBBI locations.
!   lat(:)            -- Lat struct: Each bunch has its own lat(i)
!                         with markers at parasitic crossings that bunch sees.
!   crossings(:,:)     -- Real(rp): First column is the position (as a
!                         fraction/percent of the total lat length) of a beam-
!                         beam interaction seen by a bunch. Second column is
!                         the index of the bunch that sees the crossings in the
!                         first column. (Third column should be an index
!                         indicating which crossing this is: 1rst, 2nd, etc,
!                         but it is not used or changed here.) The fourth and
!                         fifth columns should be empty and will hold,
!                         respectively, the index of the crossings in the lat,
!                         and the index in the master_lat.
!
! Output:
!   lat(:)            -- Lat struct: Lat for each bunch with
!                         markers placed where parasitic crossings are seen
!                         by that bunch.
!   master_lat        -- Lat struct: Master_lat with markers placed at every
!                         parasitic crossing (seen by any bunch).
!   master_lat_oppos  -- Lat struct: Master_lat_oppos with markers placed
!                         at every parasitic crossings (seen by any bunch).
!
!
! Note: The elements placed at the parasitic crossing sites are simply markers
!       with unit 6x6 matrices. The fourth and fifth columns of crossings hold
!       the indices of the location of all inserted markers (named
!       'LRBBI_MARKER'). lat_make_mat6 is called for all lattices.
!-

#include "CESR_platform.inc"

subroutine mark_lrbbi (master_lat, master_lat_oppos, lat, crossings)

  use bmad_struct
  use bmad_interface, except => mark_lrbbi

  implicit none

  type (lat_struct), dimension(:) :: lat
  type (lat_struct) :: master_lat, master_lat_oppos
  type (coord_struct), allocatable, save :: orbit(:), orbit_oppos(:)
  type (ele_struct) :: insert_ele

  real(rp), dimension(:,:) :: crossings
  real(rp) :: smallest, s_split

  integer :: i, j, k, loc_smallest, ix_split
  integer :: ix_ele, master_ix_split, master_ix_ele, lat_index, total
  integer :: index1, index2
  integer, dimension(1) :: minloc_array
	
  logical :: split_done, split_done_1, split_done_2, split_done_3

!

  total = size(crossings, 1)
  
  call init_ele(insert_ele)
  insert_ele%key = marker$
  insert_ele%name = "LRBBI_MARKER"
  call mat_make_unit(insert_ele%mat6)
        
!---------------------------------------------------------------------------

  do k = 1, total
    smallest = minval(crossings(k:total,1))
    minloc_array = minloc(crossings(k:total,1))
    loc_smallest = k - 1 + minloc_array(1)
    index1 = crossings(loc_smallest, 2)
    index2 = crossings(loc_smallest, 3)
    if (loc_smallest .ne. k) then
      crossings(loc_smallest,1) = crossings(k,1)
      crossings(loc_smallest, 2) = crossings(k,2)
      crossings(loc_smallest, 3) = crossings(k,3)
      crossings(k,1) = smallest
      crossings(k,2) = index1
      crossings(k,3) = index2
    endif
  enddo
  
  j_loop: do j = 1, total
     
    s_split = crossings(j,1) * master_lat%param%total_length
    lat_index = crossings(j,2)
     
    if (s_split == 0) cycle

    if (j .gt. 1 .and. crossings(j,1) == crossings(j-1,1)) then
       call split_lat(lat(lat_index), s_split, ix_split, split_done)
       if (.not. split_done) cycle j_loop
       insert_ele%s = s_split
       ix_ele = ix_split + 1
       crossings(j,4) = ix_ele
       crossings(j,5) = master_ix_ele                 
       call insert_element(lat(lat_index), insert_ele, ix_ele)
       cycle j_loop		
     endif

     ix_split = 0
     
     call split_lat(master_lat_oppos, s_split, master_ix_split,split_done_1)
     call split_lat(master_lat, s_split, master_ix_split, split_done_2)
     call split_lat(lat(lat_index), s_split, ix_split, split_done_3)
     
     insert_ele%s = s_split
     
     ix_ele = ix_split + 1
     crossings(j, 4) = ix_ele
     
     master_ix_ele = master_ix_split + 1
     crossings(j,5) = master_ix_ele

     if (split_done_2) then
        call insert_element(master_lat, insert_ele, master_ix_ele)
     endif

     if (split_done_1) then
        call insert_element(master_lat_oppos, insert_ele, master_ix_ele)
     endif

     if (split_done_3) then
        call insert_element(lat(lat_index), insert_ele, ix_ele)
     endif

  enddo j_loop

!

  call reallocate_coord (orbit,master_lat%n_ele_max)
  call reallocate_coord (orbit_oppos,master_lat%n_ele_max)

  call closed_orbit_calc (master_lat, orbit, 4)
  call closed_orbit_calc(master_lat_oppos, orbit_oppos, 4)
  
  call lat_make_mat6 (master_lat, -1, orbit)
  call lat_make_mat6 (master_lat_oppos, -1, orbit_oppos)
  
  do i = 1, size(lat)
     call closed_orbit_calc(lat(i), orbit, 4)
     call lat_make_mat6(lat(i), -1, orbit)
  enddo

!--------------------------------------------------------------------------


end subroutine
