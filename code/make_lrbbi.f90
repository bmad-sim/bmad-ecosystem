!+
! subroutine MAKE_LRBBI(master_ring, master_ring_oppos, ring, &
!										                   			ix_LRBBI, master_ix_LRBBI)
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
!   master_ring          -- Ring struct: Ring with markers at LRBBI locations.
!   master_ring_oppos    -- Ring struct: Ring for oppositely circulating
!                           particles, with markers at parasitic crossings.
!   ring(:)              -- Ring struct: Each bunch has its own ring(i)
!                           with markers at all crossings that bunch sees.
!   ix_LRBBI(:,:)        -- Real(rdef): First index (i) is the index of the ring
!                           (i.e., bunch), second index (j) is the index of a
!                           beam-beam element's position in ring(i).
!   master_ix_LRBBI(:,:) -- Real(rdef): First index (i) is the index of the
!                           ring, second index (j) is the index of a beam-beam
!                           element (seen by the ith bunch) in the master_ring
!                           and master_ring_oppos. This index is used to
!                           calculate sigmas and offsets.
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

!$Id$
!$Log$
!Revision 1.5  2003/01/27 14:40:37  dcs
!bmad_version = 56
!
!Revision 1.4  2002/02/23 20:32:18  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/08 21:44:39  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:53  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine MAKE_LRBBI(master_ring, master_ring_oppos, ring, &
													ix_LRBBI, master_ix_LRBBI)

  use bmad_struct
  use bmad_interface                     

  implicit none

	type (ring_struct), dimension(:) :: ring
	type (ring_struct) :: master_ring, master_ring_oppos
  type (coord_struct) :: orbit_(0:n_ele_maxx), orbit_oppos_(0:n_ele_maxx)
	type (modes_struct) :: mode

	real, dimension(0:n_ele_maxx) :: beta_x, beta_y, eta_x, eta_y
	real, dimension(0:n_ele_maxx) :: sigma_x, sigma_y, sigma_z 
	real :: e_spread, emit_x, emit_y

	integer :: i,j
  integer, dimension(:,:) :: ix_LRBBI
  integer, dimension(:,:) :: master_ix_LRBBI

!

	call closed_orbit_at_start(master_ring, orbit_(0), 4, .true.)
	call track_all(master_ring, orbit_)
	call twiss_at_start(master_ring)
	call twiss_propagate_all(master_ring)

	do i = 1, size(ring)
		do j = 1, size(ix_LRBBI,2)
			if (ix_LRBBI(i,j) == 0) cycle
			ring(i)%ele_(ix_LRBBI(i,j))%name = "LRBBI"
			ring(i)%ele_(ix_LRBBI(i,j))%key = beambeam$
			ring(i)%ele_(ix_LRBBI(i,j))%value(charge$) = -1
			ring(i)%ele_(ix_LRBBI(i,j))%value(n_slice$) = 1
		enddo
	enddo

  do j = 1, master_ring%n_ele_ring
    beta_x(j) = master_ring%ele_(j)%x%beta
    beta_y(j) = master_ring%ele_(j)%y%beta

    eta_x(j) = master_ring%ele_(j)%x%eta
    eta_y(j) = master_ring%ele_(j)%y%eta
  enddo

  call radiation_integrals(master_ring, orbit_, mode)
  e_spread = mode%sig_E

  do j = 1, master_ring%n_ele_ring
    sigma_z(j) = mode%sig_z
  enddo

  emit_x = mode%a%emittance
  emit_y = emit_x * 0.01

  do j = 1, master_ring%n_ele_ring
    sigma_x(j) = sqrt(emit_x * beta_x(j) + (eta_x(j)**2) * ((e_spread)**2))
    sigma_y(j) = sqrt(emit_y * beta_y(j) + (eta_y(j)**2) *((e_spread)**2))
  enddo

	do i = 1, size(ring)
		do j = 1, size(ix_LRBBI, 2)
			if (ix_lrbbi(i,j) == 0) cycle
  		ring(i)%ele_(ix_lrbbi(i,j))%value(sig_x$) = sigma_x(master_ix_LRBBI(i,j))
  		ring(i)%ele_(ix_lrbbi(i,j))%value(sig_y$) = sigma_y(master_ix_LRBBI(i,j))
  		ring(i)%ele_(ix_lrbbi(i,j))%value(sig_z$) = sigma_z(master_ix_LRBBI(i,j))
    enddo

		call twiss_at_start(master_ring_oppos)
		call twiss_propagate_all(master_ring_oppos)
		call closed_orbit_at_start(master_ring_oppos, orbit_oppos_(0), 4, .true.)
		call track_all(master_ring_oppos, orbit_oppos_)

		do j = 1, size(ix_LRBBI, 2)
			if (ix_lrbbi(i,j) == 0) cycle
  		ring(i)%ele_(ix_lrbbi(i,j))%value(x_offset$) = orbit_oppos_(master_ix_lrbbi(i,j))%vec(1)
  		ring(i)%ele_(ix_lrbbi(i,j))%value(y_offset$) = orbit_oppos_(master_ix_lrbbi(i,j))%vec(3)
		enddo
        
  enddo


end subroutine
