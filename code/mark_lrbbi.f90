!+
! subroutine MARK_LRBBI(master_ring, master_ring_oppos, ring, crossings)
!
!	Input:
!		master_ring						-- Ring struct: Ring with markers at all LRBBI crossings
!		master_ring_oppos  -- Ring struct: Ring for oppositely circulating
!																particles, with markers at all LRBBI crossings.
!   ring									-- Array of ring structs: Each bunch has its own ring(i)
!                         with markers at all LRBBI crossings that bunch sees.
!		crossings							-- Array(i,j): First column is the position (as a
!																fraction/percent of the total ring length) of a beam-
!																beam interaction seen by a bunch. Second column is
!																the index of the bunch that sees the crossings in the
!																first column. (Third column should be an index
!																indicating which crossing this is: 1rst, 2nd, etc,
!																but it is not used or changed here.) The fourth and
!																fifth columns should be empty and will hold,
!																respectively, the index of the crossings in the ring,
!																and the index in the master_ring.
!
!	Output:			
!		ring									-- Array of ring structs: Ring for each bunch with
!																markers placed where parasitic crossings are seen
!																by that bunch.
!   master_ring						-- Ring struct: Master_ring with markers placed at every 
!                         parasitic crossing (seen by any bunch).
!   master_ring_oppos		-- Ring struct: Master_ring_oppos with markers placed
!                         at every parasitic crossings (seen by any bunch).
!   ix_LRBBI							-- Array(i,j): First column (i) is the index of the ring
!																(i.e., bunch), second column (j) is the index of a
!																beam-beam element's positions in ring(i).
!		master_ix_LRBBI    -- Array(i,j): First column (i) is the index of the
!																ring, second column (j) is the index of a beam-beam
!																element (seen by the ith bunch) in the master_ring 
!																and master_ring_oppos. This index is used to 
!																calculate sigmas and offsets.
!
!
! Note: The elements placed at the parasitic crossing sites are simply markers 
!       with unit 6x6 matrices. ix_LRBBI and master_ix_LRBBI hold the
!				indices of the location of all inserted markers (named 'LRBBI_MARKER').
!       Ring_make_mat6 is called for all rings.
!-

subroutine MARK_LRBBI(master_ring, master_ring_oppos, ring, crossings)

  use bmad_struct                      
  use bmad_interface

  implicit none

	type (ring_struct), dimension(:) :: ring
	type (ring_struct) :: master_ring, master_ring_oppos
  type (coord_struct) :: orbit_(0:n_ele_maxx), orbit_oppos_(0:n_ele_maxx)
  type (ele_struct), dimension(:), allocatable ::  insert_ele

	real, dimension(:,:) :: crossings
	real :: smallest, s_split

  integer :: i, j, k, m, ierr, end, loc_smallest, ix_split
	integer :: ix_ele, master_ix_split, master_ix_ele, ring_index, total
	integer :: index1, index2, index3, index4
	integer, dimension(1) :: minloc_array
	
	logical :: ok, split_done, split_done_1, split_done_2, split_done_3

!

	total = size(crossings, 1)

	allocate(insert_ele(1:total), stat=ierr)
	if(ierr .ne. 0) then
		print*, "INSERT_ELE: ALLOCATION REQUEST DENIED."
		call err_exit
	endif

!---------------------------------------------------------------------------

	end = total

	do k = 1, end
		smallest = minval(crossings(k:end,1))
		minloc_array = minloc(crossings(k:end,1))
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

	jay: do j = 1, end

		s_split = crossings(j,1) * master_ring%param%total_length
		ring_index = crossings(j,2)

		if (s_split == 0) cycle

		if (j .gt. 1 .and. crossings(j,1) == crossings(j-1,1)) then
			call split_ring(ring(ring_index), s_split, ix_split, split_done)
			if (split_done == .false.) cycle jay

			call init_ele(insert_ele(j))
			insert_ele(j)%key = marker$
			insert_ele(j)%name = "LRBBI_MARKER"
			call mat_unit(insert_ele(j)%mat6, 6, 6)

			insert_ele(j)%s = s_split
     	ix_ele = ix_split + 1
			crossings(j,4) = ix_ele
			crossings(j,5) = master_ix_ele

 		  call insert_element(ring(ring_index), insert_ele(j), ix_ele)

			cycle jay		
		endif

    ix_split = 0

    call split_ring(master_ring_oppos, s_split, master_ix_split,split_done_1)
 		call split_ring(master_ring, s_split, master_ix_split, split_done_2)
 		call split_ring(ring(ring_index), s_split, ix_split, split_done_3)
 
		call init_ele(insert_ele(j))
		insert_ele(j)%key = marker$
		insert_ele(j)%name = "LRBBI_MARKER"
		call mat_unit(insert_ele(j)%mat6, 6, 6)

		insert_ele(j)%s = s_split
     
		ix_ele = ix_split + 1
		crossings(j, 4) = ix_ele

		master_ix_ele = master_ix_split + 1
		crossings(j,5) = master_ix_ele

		if (split_done_2 == .true.) then
			call insert_element(master_ring, insert_ele(j), master_ix_ele)
		endif

		if (split_done_1 == .true.) then
			call insert_element(master_ring_oppos, insert_ele(j), master_ix_ele)
    endif

		if (split_done_3 == .true.) then
 	  	call insert_element(ring(ring_index), insert_ele(j), ix_ele)
    endif

  enddo jay

	call closed_orbit_at_start(master_ring, orbit_(0), 4, .true.)
	call track_all(master_ring, orbit_)
	call closed_orbit_at_start(master_ring_oppos, orbit_oppos_(0), 4, .true.)
	call track_all(master_ring_oppos, orbit_oppos_)

 	call ring_make_mat6(master_ring, -1, orbit_)
 	call ring_make_mat6(master_ring_oppos, -1, orbit_oppos_)

	do i = 1, size(ring)
		call closed_orbit_at_start(ring(i), orbit_(0), 4, .true.)
		call track_all(ring(i), orbit_)
 		call ring_make_mat6(ring(i), -1, orbit_)
  enddo

!--------------------------------------------------------------------------

	deallocate(insert_ele, stat=ierr)
  if (ierr .ne. 0) then
    print*, "INSERT_ELE: DEALLOCATION REQUEST DENIED."
		call err_exit
  endif

end subroutine
