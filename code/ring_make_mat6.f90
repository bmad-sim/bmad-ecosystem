!+
! Subroutine ring_make_mat6 (ring, ix_ele, coord_)
!
! Subroutine to make the first order transfer map:
!   r_out = M * r_in + vec0
! M is the 6x6 linear transfer matrix (Jacobian) about the 
! reference orbit coord_.
!
! The routine will also call control_bookkeeper to make sure that all
! lord/slave dependencies are correct.
!
! Moudules Needed:
!   use bmad
!
! Input:
!   ring       -- Ring_struct: Ring containing the elements.
!   ix_ele     -- Integer: Index of the element. if < 0 then entire
!                    ring will be made. In this case group elements will
!                    be made up last.
!   coord_(0:) -- Coord_struct, optional: Coordinates of the reference orbit
!                   around which the matrix is calculated. If not present 
!                   then the referemce is taken to be the origin.
!
! Output:
!   ring        -- ring_struct:
!     ele_(:)%mat6  -- Real(rp): 1st order (Jacobian) 6x6 transfer matrix.
!     ele_(:)%vec0  -- Real(rp): 0th order transfer vector.
!-

#include "CESR_platform.inc"

recursive subroutine ring_make_mat6 (ring, ix_ele, coord_)

  use bmad_struct
  use bmad_utils_mod
  use bmad_interface, except => ring_make_mat6
  use bookkeeper_mod, only: control_bookkeeper

  implicit none
                                         
  type (ring_struct), target :: ring
  type (coord_struct), optional, volatile :: coord_(0:)
  type (ele_struct), pointer :: ele

  integer i, j, ie, ix_ele, i1, ix_taylor(100), n_taylor

! Error check

  if (ix_ele == 0 .or. ix_ele > ring%n_ele_max) then
    print *, 'ERROR IN RING_MAKE_MAT6: ELEMENT INDEX OUT OF BOUNDS:', ix_ele
    if (bmad_status%exit_on_error) call err_exit
    return
  endif

  if (present(coord_)) then
    if (ubound(coord_, 1) < ring%n_ele_use) then
      print *, 'ERROR IN RING_MAKE_MAT6: COORD_(:) ARGUMENT SIZE IS TOO SMALL!'
      call err_exit
    endif
  endif

  call compute_element_energy (ring)

!--------------------------------------------------------------
! make entire ring if ix_ele < 0
! first do the inter element bookkeeping

  if (ix_ele < 0) then         

    call control_bookkeeper (ring)

! now make the transfer matrices.
! for speed if an element needs a taylor series then check if we can use
! one from a previous element.

    n_taylor = 0  ! number of taylor series found

    do i = 1, ring%n_ele_use

      ele => ring%ele_(i)

      if (ele%mat6_calc_method == taylor$ .and. ele%key == wiggler$) then
        if (.not. associated(ele%taylor(1)%term)) then
          do j = 1, n_taylor
            ie = ix_taylor(j)
            if (.not. equivalent_eles (ele, ring%ele_(ie))) cycle
            if (present(coord_)) then
              if (any(coord_(i-1)%vec /= coord_(ie-1)%vec)) cycle
            endif
            call transfer_ele_taylor (ring%ele_(ie), ele, &
                                                 ring%ele_(ie)%taylor_order)
            exit
          enddo
        endif
        n_taylor = n_taylor + 1
        ix_taylor(n_taylor) = i
      endif

      if (present(coord_)) then
        call make_mat6(ele, ring%param, coord_(i-1), coord_(i), .true.)
      else
        call make_mat6(ele, ring%param)
      endif
    enddo

    return

  endif

!-----------------------------------------------------------
! otherwise make a single element

  call control_bookkeeper (ring, ix_ele)

! for a regular element

  if (ix_ele <= ring%n_ele_use) then
     if (present(coord_)) then
        call make_mat6(ring%ele_(ix_ele), ring%param, &
                                  coord_(ix_ele-1), coord_(ix_ele), .true.)
     else
        call make_mat6(ring%ele_(ix_ele), ring%param)
     endif

    return
  endif                        

! for a control element

  do i1 = ring%ele_(ix_ele)%ix1_slave, ring%ele_(ix_ele)%ix2_slave
    i = ring%control_(i1)%ix_slave
     if (present(coord_)) then
        call ring_make_mat6 (ring, i, coord_)
     else
        call ring_make_mat6 (ring, i)
     endif
  enddo

end subroutine
