!+
! Subroutine lat_make_mat6 (lat, ix_ele, coord)
!
! Subroutine to make the first order transfer map for an element:
!   r_out = M * r_in + vec0
! M is the 6x6 linear transfer matrix (Jacobian) about the 
! reference orbit coord.
!
! If the element lat%ele(ix_ele) is a lord element then the martices of 
! all the slave elements will be recomputed.
!
! The routine will also call control_bookkeeper to make sure that all
! lord/slave dependencies are correct.
!
! Moudules Needed:
!   use bmad
!
! Input:
!   lat       -- lat_struct: Lat containing the elements.
!   ix_ele    -- Integer, optional: Index of the element. if not present
!                   or negative then the entire lattice will be made.
!   coord(0:) -- Coord_struct, optional: Coordinates of the reference orbit
!                   around which the matrix is calculated. If not present 
!                   then the referemce is taken to be the origin.
!
! Output:
!   lat        -- lat_struct:
!     ele(:)%mat6  -- Real(rp): 1st order (Jacobian) 6x6 transfer matrix.
!     ele(:)%vec0  -- Real(rp): 0th order transfer vector.
!-

#include "CESR_platform.inc"

recursive subroutine lat_make_mat6 (lat, ix_ele, coord)

  use bmad_struct
  use bmad_utils_mod
  use bmad_interface, except => lat_make_mat6
  use bookkeeper_mod, only: control_bookkeeper

  implicit none
                                         
  type (lat_struct), target :: lat
  type (coord_struct), optional, volatile :: coord(0:)
  type (coord_struct) orb_start, orb_end
  type (ele_struct), pointer :: ele

  integer, optional :: ix_ele
  integer i, j, ie, i1, n_taylor, i_ele
  integer, save, allocatable :: ix_taylor(:)

  logical transferred, want_taylor

! Error check

  if (.not. allocated(ix_taylor)) allocate(ix_taylor(200))

  i_ele = integer_option (-1, ix_ele)

  if (i_ele == 0 .or. i_ele > lat%n_ele_max) then
    print *, 'ERROR IN lat_make_mat6: ELEMENT INDEX OUT OF BOUNDS:', i_ele
    if (bmad_status%exit_on_error) call err_exit
    return
  endif

  if (present(coord)) then
    if (ubound(coord, 1) < lat%n_ele_track) then
      print *, 'ERROR IN lat_make_mat6: coord(:) ARGUMENT SIZE IS TOO SMALL!'
      call err_exit
    endif
  endif

  if (bmad_com%auto_bookkeeper) call compute_reference_energy (lat)

!--------------------------------------------------------------
! Make entire lat if i_ele < 0.
! First do the inter-element bookkeeping.

  if (i_ele < 0) then         

    if (bmad_com%auto_bookkeeper) call control_bookkeeper (lat)

! Now make the transfer matrices.
! For speed if a element needs a taylor series then check if we can use
! one from a previous element.

! For consistancy, if no orbit is given, the starting coords in a super_slave
! will be taken as the ending coords of the previous super_slave.

    n_taylor = 0  ! number of taylor series found
    orb_end%vec = 0

    do i = 1, lat%n_ele_track

      ele => lat%ele(i)
      want_taylor = (ele%mat6_calc_method == taylor$) .or. &
                    (ele%mat6_calc_method == symp_map$) .or. &
                    (ele%tracking_method == taylor$) .or. &
                    (ele%tracking_method == symp_map$)

      if (want_taylor) then
        transferred = .false.
        if (.not. associated(ele%taylor(1)%term)) then
          do j = 1, n_taylor
            ie = ix_taylor(j)
            if (.not. equivalent_eles (ele, lat%ele(ie))) cycle
            if (present(coord)) then
              if (any(coord(i-1)%vec /= coord(ie-1)%vec)) cycle
            endif
            call transfer_ele_taylor (lat%ele(ie), ele, &
                                                 lat%ele(ie)%taylor_order)
            transferred = .true.
            exit
          enddo
        endif
      endif

      if (present(coord)) then
        call make_mat6(ele, lat%param, coord(i-1), coord(i), .true.)
      elseif (ele%control_type == super_slave$) then
        orb_start = orb_end
        call make_mat6(ele, lat%param, orb_start, orb_end)
      else
        call make_mat6(ele, lat%param)
        ! Reset orb_end if not in a superposition block.
        if (ele%value(l$) /= 0) orb_end%vec = 0  
      endif

      ! save this taylor in the list if it is a new one. 

      if (want_taylor .and. .not. transferred) then
        n_taylor = n_taylor + 1
        if (n_taylor > size(ix_taylor)) &
                         call re_allocate (ix_taylor, 2*size(ix_taylor))
        ix_taylor(n_taylor) = i
      endif

    enddo

    return

  endif

!-----------------------------------------------------------
! otherwise make a single element

  call control_bookkeeper (lat, i_ele)

! For an element in the tracking part of the lattice

  if (i_ele <= lat%n_ele_track) then
     if (present(coord)) then
        call make_mat6(lat%ele(i_ele), lat%param, &
                                  coord(i_ele-1), coord(i_ele), .true.)
     else
        call make_mat6(lat%ele(i_ele), lat%param)
     endif

    return
  endif                        

! for a control element

  do i1 = lat%ele(i_ele)%ix1_slave, lat%ele(i_ele)%ix2_slave
    i = lat%control(i1)%ix_slave
     if (present(coord)) then
        call lat_make_mat6 (lat, i, coord)
     else
        call lat_make_mat6 (lat, i)
     endif
  enddo

end subroutine
