!+
! subroutine MAKE_LRBBI(master_ring_oppos, ring, ix_LRBBI, master_ix_LRBBI_oppos)
!
! This subroutine takes rings with markers for LRBBI's (like those created
!   in MARK_LRBBI) and turns the marked locations into LRBBI elements. The
!   master_ring and master_ring_oppos supply the data for calculating sigmas
!   and offsets.
!
! Modules Needed:
!   use bmad
!
! Input:
!   master_ring_oppos    -- Ring struct: Ring for oppositely circulating
!                           particles, with markers at parasitic crossings.
!   ring(:)              -- Ring struct: Each bunch has its own ring(i)
!                           with markers at all crossings that bunch sees.
!   ix_LRBBI(:,:)        -- Real(rp): First index (i) is the index of the ring
!                           (i.e., bunch), second index (j) is the index of a
!                           beam-beam element's position in ring(i).
!   master_ix_LRBBI_oppos(:,:) -- Real(rp): First index (i) is the index of the
!                           ring, second index (j) is the index of a beam-beam
!                           element (seen by the ith bunch) in the master_ring_oppos.
!                           This index is used to calculate sigmas and offsets.
!
! Output:
!   ring(:)              -- Ring struct: Ring for each bunch with
!                           markers updated to beam-beam elements.
!
! NOTE: The following attributes of the LRBBI elements are changed:
!         %name = "LRBBI"
!         %key = beambeam
!         %value(charge) = -1
!         %value(n_slice) = 1
!         %value(sig_x), (sig_y), (sig_z), (x_offset), (y_offset)
!           (These are assigned values from master_ring and master_ring_oppos.)
!    All other values (i.e., %s_pos) must be set elsewhere!! Also, this does
!     NOT call ring_make_mat6.
!-

#include "CESR_platform.inc"

subroutine MAKE_LRBBI(master_ring_oppos, ring, ix_LRBBI, master_ix_LRBBI_oppos)

  use bmad_struct
  use bmad_interface                     

  implicit none

  type (ring_struct), dimension(:) :: ring
  type (ring_struct) :: master_ring_oppos
  type (coord_struct), allocatable, save :: orbit_oppos_(:)
  type (modes_struct) :: mode_oppos

  real(rp) :: beta_x, beta_y, eta_x, eta_y
  real(rp) :: sigma_z
  real(rp) :: e_spread, emit_x, emit_y
  real(rp) :: sigma_x(master_ring_oppos%n_ele_ring), sigma_y(master_ring_oppos%n_ele_ring)

  integer :: i,j
  integer, dimension(:,:) :: ix_LRBBI
  integer, dimension(:,:) :: master_ix_LRBBI_oppos
  
! init

  call reallocate_coord (orbit_oppos_, master_ring_oppos%n_ele_max)

! calc master_ring_opps parameters

  call twiss_at_start(master_ring_oppos)
  call twiss_propagate_all(master_ring_oppos)
  call closed_orbit_at_start(master_ring_oppos, orbit_oppos_(0), 4, .true.)
  if (.not. bmad_status%ok) return
  call track_all(master_ring_oppos, orbit_oppos_)
  call twiss_at_start(master_ring_oppos)
  if (.not. bmad_status%ok) return
  call twiss_propagate_all(master_ring_oppos)
  call radiation_integrals(master_ring_oppos, orbit_oppos_, mode_oppos)

!

  do i = 1, size(ring)
    do j = 1, size(ix_LRBBI,2)
      if (ix_LRBBI(i,j) == 0) cycle
      ring(i)%ele_(ix_LRBBI(i,j))%name = "LRBBI"
      ring(i)%ele_(ix_LRBBI(i,j))%key = beambeam$
      ring(i)%ele_(ix_LRBBI(i,j))%value(charge$) = -1
      ring(i)%ele_(ix_LRBBI(i,j))%value(n_slice$) = 1
    enddo
  enddo

  e_spread = mode_oppos%sigE_E
  sigma_z = mode_oppos%sig_z

  emit_x = mode_oppos%a%emittance
  emit_y = emit_x * 0.01
  do j = 1, master_ring_oppos%n_ele_ring
    beta_x = master_ring_oppos%ele_(j)%x%beta
    beta_y = master_ring_oppos%ele_(j)%y%beta
    eta_x = master_ring_oppos%ele_(j)%x%eta
    eta_y = master_ring_oppos%ele_(j)%y%eta
    sigma_x(j) = sqrt(emit_x * beta_x + (eta_x**2) * ((e_spread)**2))
    sigma_y(j) = sqrt(emit_y * beta_y + (eta_y**2) *((e_spread)**2))
  enddo

  do i = 1, size(ring)

    do j = 1, size(ix_LRBBI, 2)
      if (ix_lrbbi(i,j) == 0) cycle
      ring(i)%ele_(ix_lrbbi(i,j))%value(sig_x$) = sigma_x(master_ix_LRBBI_oppos(i,j))
      ring(i)%ele_(ix_lrbbi(i,j))%value(sig_y$) = sigma_y(master_ix_LRBBI_oppos(i,j))
      ring(i)%ele_(ix_lrbbi(i,j))%value(sig_z$) = sigma_z
    enddo
     
    do j = 1, size(ix_LRBBI, 2)
      if (ix_lrbbi(i,j) == 0) cycle
      ring(i)%ele_(ix_lrbbi(i,j))%value(x_offset$) = orbit_oppos_(master_ix_lrbbi_oppos(i,j))%vec(1)
      ring(i)%ele_(ix_lrbbi(i,j))%value(y_offset$) = orbit_oppos_(master_ix_lrbbi_oppos(i,j))%vec(3)
    enddo
        
  enddo

end subroutine
