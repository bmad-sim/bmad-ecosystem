!+
! Subroutine ring_make_mat6 (ring, ix_ele, coord_)
!
! Subroutine to make the 6x6 linear transfer matrix for an element or
! elements in a ring.
!
! Moudules Needed:
!   use bmad
!
! Input:
!   ring     -- Ring_struct: Ring containing the elements.
!   ix_ele   -- Integer: Index of the element. if < 0 then entire
!                    ring will be made. In this case group elements will
!                    be made up last.
!   coord_(0:n_ele_maxx) 
!            -- Coord_struct, optional: Coordinates of the 
!                   nominal orbit around which the matrix is calculated. 
!                   If not present then the orbit is taken to be the orign.
!
! Output:
!   ring    -- ring_struct:
!     %ele_(i)%mat6 -- 6x6 transfer matrices.
!-

!$Id$
!$Log$
!Revision 1.8  2003/01/27 14:40:42  dcs
!bmad_version = 56
!
!Revision 1.7  2002/11/17 01:01:43  dcs
!compiler bug fix
!
!Revision 1.6  2002/11/06 06:48:32  dcs
!Changed arg array
!
!Revision 1.5  2002/06/13 14:54:28  dcs
!Interfaced with FPP/PTC
!
!Revision 1.4  2002/02/23 20:32:23  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/16 21:04:18  helms
!Fixed problem with passing optional arguments.
!
!Revision 1.2  2001/09/27 18:31:57  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


recursive subroutine ring_make_mat6 (ring, ix_ele, coord_)

  use bmad_struct
  use bmad_interface, only: control_bookkeeper, make_mat6

  implicit none
                                         
  type (ring_struct), target :: ring
  type (coord_struct), optional, volatile :: coord_(0:)
  type (ele_struct), pointer :: ele
  type (coord_struct) c1

  integer i, j, k, ie, ix_ele, i1, i2, i3, ix1, ix2, ix3
  integer ix_taylor(100), n_taylor, istat


!

  if (ring%param%lattice_type == linac_lattice$)  call compute_element_energy (ring)

! Check Energy

  if (abs(ring%param%beam_energy-1d9*ring%param%energy) > &
                                             1e-5*ring%param%beam_energy) then
    print *, 'ERROR IN RING_MAKE_MAT6:'
    print *, '      RING%PARAM%ENERGY AND RING%PARAM%BEAM_ENERGY DO NOT MATCH'
    print *, '      ', ring%param%energy, ring%param%beam_energy 
    call err_exit
  endif

! Error check

  if (ix_ele == 0 .or. ix_ele > ring%n_ele_max) then
    type *, 'ERROR IN RING_MAKE_MAT6: ELEMENT INDEX OUT OF BOUNDS:', ix_ele
    if (bmad_status%exit_on_error) call err_exit
    return
  endif

!--------------------------------------------------------------
! make entire ring if ix_ele < 0
! first do the inter element bookkeeping

  if (ix_ele < 0) then         

    do i = ring%n_ele_ring+1, ring%n_ele_max
      if (ring%ele_(i)%control_type /= group_lord$)  &
                                 call control_bookkeeper (ring, i)
    enddo

    do i = ring%n_ele_ring+1, ring%n_ele_max
      if (ring%ele_(i)%control_type == group_lord$)  &
                                 call control_bookkeeper (ring, i)
    enddo

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
            do k = 1, 6
              deallocate (ele%taylor(k)%term, stat = istat)
              allocate (ele%taylor(k)%term(size(ring%ele_(ie)%taylor(k)%term)))
              ele%taylor(k)%term = ring%ele_(ie)%taylor(k)%term
            enddo
            exit
          enddo
        endif
        n_taylor = n_taylor + 1
        ix_taylor(n_taylor) = i
      endif

      if (present(coord_)) then
        call make_mat6(ele, ring%param, coord_(i-1), c1)
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

  if (ix_ele <= ring%n_ele_ring) then
     if (present(coord_)) then
        call make_mat6(ring%ele_(ix_ele), ring%param, coord_(ix_ele-1), c1)
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

!----------------------------------------------------------------
contains

function equivalent_eles(ele1, ele2) result (equiv)

  type (ele_struct) ele1, ele2
  logical equiv
  integer it

!

  equiv = .false.

  if (any(ele1%value /= ele2%value)) return
  if (associated(ele1%wig_term) .xor. associated(ele2%wig_term)) return
  if (ele1%num_steps /= ele2%num_steps) return
  if (ele1%integration_order /= ele2%integration_order) return

  if (associated(ele1%wig_term)) then
    if (size(ele1%wig_term) /= size(ele2%wig_term)) return
    do it = 1, size(ele1%wig_term)
      if (ele1%wig_term(it)%coef /= ele2%wig_term(it)%coef) return
      if (ele1%wig_term(it)%kx /= ele2%wig_term(it)%kx) return
      if (ele1%wig_term(it)%ky /= ele2%wig_term(it)%ky) return
      if (ele1%wig_term(it)%kz /= ele2%wig_term(it)%kz) return
      if (ele1%wig_term(it)%phi_z /= ele2%wig_term(it)%phi_z) return
    enddo
  endif

  equiv = .true.

end function


end subroutine
