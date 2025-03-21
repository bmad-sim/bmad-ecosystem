!+
! Subroutine lat_make_mat6 (lat, ix_ele, ref_orb, ix_branch, err_flag)
!
! Subroutine to make the first order transfer map for an element or elements:
!   r_out = M * r_in + vec0
! M is the 6x6 linear transfer matrix (Jacobian) about the 
! reference orbit ref_orb.
!
! If the element lat%ele(ix_ele) is a lord element then the martices of 
! all the slave elements will be recomputed.
!
! Input:
!   lat         -- lat_struct: Lat containing the elements.
!   ix_ele      -- Integer, optional: Index of the element. If not present
!                    or negative, the matrices for all elements will be calculated.
!   ref_orb(0:) -- Coord_struct, optional: Coordinates of the reference orbit
!                   around which the matrix is calculated. If not present 
!                   then the referemce is taken to be the origin.
!   ix_branch   -- Integer, optional: Branch index. Default is 0 (main lattice).
!                   -1 => All branches/all elements (ref_orb & ix_ele will be ignored).
!
! Output:
!   lat        -- lat_struct:
!     ele(:)%mat6  -- Real(rp): 1st order (Jacobian) 6x6 transfer matrix.
!     ele(:)%vec0  -- Real(rp): 0th order transfer vector.
!   err_flag    -- Logical, optional: True if there is an error. False otherwise.
!-

recursive subroutine lat_make_mat6 (lat, ix_ele, ref_orb, ix_branch, err_flag)

use bmad_interface, except_dummy => lat_make_mat6

implicit none
                                       
type (lat_struct), target :: lat
type (coord_struct), optional :: ref_orb(0:)
type (coord_struct) orb_start
type (ele_struct), pointer :: ele, slave, lord, slave0, slave1
type (branch_struct), pointer :: branch

real(rp), pointer :: mat6(:,:), vec0(:)

integer, optional :: ix_ele, ix_branch
integer i, j, i0, i1, ie, ild, n_taylor, i_ele, i_branch, ix_slave, species

logical, optional :: err_flag
logical transferred, zero_orbit, err, finished

character(*), parameter :: r_name = 'lat_make_mat6'

! 

if (associated(lat_make_mat6_hook_ptr)) then
  finished = .false.
  call lat_make_mat6_hook_ptr(finished, lat, ix_ele, ref_orb, ix_branch, err_flag)
  if (finished) return
endif

!

if (present(err_flag)) err_flag = .true.

i_ele = integer_option (-1, ix_ele)
i_branch = integer_option (0, ix_branch)

if (i_branch == -1) then
  do i = 0, ubound(lat%branch, 1)
    call lat_make_mat6 (lat, ix_branch = i, err_flag = err)
    if (err) return
  enddo
  if (present(err_flag)) err_flag = .false.
  return
endif

! Error check

branch => lat%branch(i_branch)
if (present(ref_orb)) then
  species = ref_orb(0)%species
else
  species = default_tracking_species(branch%param)
endif

if (i_ele == 0 .or. i_ele > branch%n_ele_max) then
  call out_io (s_fatal$, r_name, 'ELEMENT INDEX OUT OF BOUNDS: \i0\ ', i_ele)
  if (global_com%exit_on_error) call err_exit
  return
endif

if (present(ref_orb)) then
  if (ubound(ref_orb, 1) < branch%n_ele_track) then
    call out_io (s_fatal$, r_name, 'REF_ORB(:) ARRAY SIZE IS TOO SMALL!')
    if (global_com%exit_on_error) call err_exit
    return
  endif
endif

if (bmad_com%auto_bookkeeper) then
  call lat_compute_ref_energy_and_time (lat, err)
  if (err) return
endif

!--------------------------------------------------------------
! Make entire lat if i_ele < 0.
! First do the inter-element bookkeeping.

if (i_ele < 0) then         

  ! Is the reference orbit zero?

  zero_orbit = .true.
  if (present(ref_orb)) then
    do i = 0, branch%n_ele_track
      if (any(ref_orb(i)%vec /= 0)) then
        zero_orbit = .false.
        exit
      endif
    enddo
  endif

  if (bmad_com%auto_bookkeeper) call control_bookkeeper (lat)

  ! Now make the transfer matrices.
  ! For speed if a element needs a taylor series then check if we can use
  ! one from a previous element.

  n_taylor = 0  ! number of taylor map found

  do i = 1, branch%n_ele_track
    ele => branch%ele(i)

    ! Check if transfer matrix needs to be recomputed

    if (.not. bmad_com%auto_bookkeeper .and. ele%bookkeeping_state%mat6 /= stale$) then
      if (present(ref_orb)) then
        if (all(ref_orb(i-1)%vec == ele%map_ref_orb_in%vec)) cycle
      else
        if (all(ele%map_ref_orb_in%vec == 0)) cycle
      endif
    endif

    ! Find orbit at start of element.

    if (zero_orbit) then 
      orb_start%vec = 0  ! Default if no previous slave found.
      if (ele%slave_status == super_slave$) then
        do ild = 1, ele%n_lord
          lord => pointer_to_lord(ele, ild, ix_slave_back = ix_slave)
          if (lord%lord_status /= super_lord$) cycle
          if (ix_slave == 1) cycle  ! If first one then no preceeding slave
          slave0 => pointer_to_slave(lord, ix_slave-1) ! slave before this element.
          orb_start = slave0%map_ref_orb_out
          exit
        enddo
      endif
      call init_coord (orb_start, orb_start%vec, ele, upstream_end$, species)

    else  ! else ref_orb must be present
      orb_start = ref_orb(i-1)
    endif

    ! If a Taylor map is needed then check for an appropriate map from a previous element.

    transferred = .false.

    if (.not. associated(ele%taylor(1)%term) .and. ((ele%mat6_calc_method == taylor$) .or. &
                                                               (ele%tracking_method == taylor$))) then
      do j = 1, n_taylor
        ie = branch%ele(j)%iyy   ! Keep track of where Taylor maps are.
        if (.not. equivalent_taylor_attributes (ele, branch%ele(ie))) cycle
        if (any(orb_start%vec /= branch%ele(ie)%taylor%ref)) cycle
        call transfer_ele_taylor (branch%ele(ie), ele)
        transferred = .true.
        exit
      enddo
    endif

    ! call make_mat6 for this element
    ! For consistancy, if no orbit is given, the starting coords in a super_slave
    ! will be taken as the ending coords of the previous super_slave.

    if (zero_orbit) then 
      call make_mat6(ele, branch%param, orb_start, err_flag = err)
    else  ! else ref_orb must be present
      if (ref_orb(i)%state /= alive$ .or. ref_orb(i-1)%state /= alive$) then
        ele%mat6 = 0
        ele%vec0 = 0
        ele%map_ref_orb_in%vec = real_garbage$
        cycle
      endif
      call make_mat6(ele, branch%param, ref_orb(i-1), err_flag = err)
    endif
    if (err) return

    ! save this taylor in the list if it is a new one. 

    if (associated(ele%taylor(1)%term) .and. .not. transferred) then
      n_taylor = n_taylor + 1
      branch%ele(n_taylor)%iyy = i  ! Keep track of where Taylor maps are.
    endif

    call set_lords_status_stale (ele, mat6_group$)
    ele%bookkeeping_state%mat6 = ok$
  enddo

  if (branch%param%bookkeeping_state%mat6 == stale$) branch%param%bookkeeping_state%mat6 = ok$

  ! calc super_lord matrices

  do i = branch%n_ele_track+1, branch%n_ele_max
    lord => branch%ele(i)
    if (lord%lord_status /= super_lord$) cycle
    slave0 => pointer_to_slave(lord, 1)
    slave1 => pointer_to_slave(lord, lord%n_slave)
    i0 = slave0%ix_ele; i1 = slave1%ix_ele

    if (.not. bmad_com%auto_bookkeeper .and. lord%bookkeeping_state%mat6 /= stale$) then
      if (present(ref_orb)) then
        if (all(ref_orb(i0-1)%vec == lord%map_ref_orb_in%vec)) cycle
      else
        if (all(lord%map_ref_orb_in%vec == 0)) cycle
      endif
    endif


    if (zero_orbit .or. i_branch /= slave0%ix_branch) then
      call make_mat6(lord, lat%branch(slave0%ix_branch)%param, err_flag = err)
    else
      call make_mat6(lord, lat%branch(slave0%ix_branch)%param, ref_orb(i0-1), err_flag = err)
    endif
    if (err) return
  enddo 

  lat%lord_state%mat6 = ok$

  if (present(err_flag)) err_flag = .false.
  return

endif

!-----------------------------------------------------------
! otherwise make a single element

ele => branch%ele(i_ele)
if (bmad_com%auto_bookkeeper) call control_bookkeeper (lat, ele)

! Check if transfer matrix needs to be recomputed

if (.not. bmad_com%auto_bookkeeper .and. ele%bookkeeping_state%mat6 /= stale$ .and. &
                                         i_ele <= branch%n_ele_track) then
  if (present(ref_orb)) then
    if (all(ref_orb(i_ele-1)%vec == ele%map_ref_orb_in%vec)) then
      return
      if (present(err_flag)) err_flag = .false.
    endif
  else
    if (all(ele%map_ref_orb_in%vec == 0)) then
      return
      if (present(err_flag)) err_flag = .false.
    endif
  endif
endif

! For an element in the tracking part of the lattice

if (i_ele <= branch%n_ele_track) then
  if (present(ref_orb)) then
     call make_mat6(ele, branch%param, ref_orb(i_ele-1), err_flag = err_flag)
  else
     call make_mat6(ele, branch%param, err_flag = err_flag)
  endif
  return
endif                        

! for a control element

if (ele%lord_status == super_lord$) then
  mat6 => ele%mat6
  call mat_make_unit(mat6)
  vec0 => ele%vec0
  vec0 = 0
endif

do i = 1, ele%n_slave
  slave => pointer_to_slave(ele, i)

  if (present(ref_orb)) then
    call lat_make_mat6 (lat, slave%ix_ele, ref_orb, slave%ix_branch, err_flag = err_flag)
  else
    call lat_make_mat6 (lat, slave%ix_ele, ix_branch = slave%ix_branch, err_flag = err_flag)
  endif

  if (ele%lord_status == super_lord$) then
    mat6 = matmul(slave%mat6, mat6)
    vec0 = matmul(slave%mat6, vec0) + slave%vec0
  endif
enddo

end subroutine
