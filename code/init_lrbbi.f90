!+
!
! Subroutine INIT_LRBBI(ring, oppos_ring, LRBBI_ele, ix_LRBBI, ix_oppos)
!
! Subroutine to calculate the basic parameters of a beambeam element. 
! Initializes the element and establishes the following values:
!     key, value(n_slice), value(charge), value(current -- provided 
!     ring%param%n_part has been specified), x_offset, y_offset,  
!     sigma_x, sigma_y, sigma_z. 
! Names the element 'LRBBI'.
! DOES NOT establish any other values, including %s (s_positions)!!
!                                                                     
!
! Modules needed:
!   use bmad
!
! Input:
!   ring            -- Ring_struct: ring in which LRBBI's are inserted
!   oppos_ring      -- Ring_sturct: ring for the oppositely circulating beam
!   LRBBI_ele       -- Ele_struct: ele_struct for the lrbbi element
!   ix_LRBBI        -- Integer: index of LRBBI element in ring (probably the 
!                        same as ix_split from split_ring since the
!                        LRBBI element has zero length)
!   ix_oppos
!
! Output:
!   LRBBI_ele       -- Ele_struct: updated lrbbi structure 
!
!-  

!$Id$
!$Log$
!Revision 1.3  2002/02/23 20:32:16  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:52  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine init_LRBBI(ring, oppos_ring, LRBBI_ele, ix_LRBBI, ix_oppos)

  use bmad
  
  implicit none

  type (ring_struct)  ring
  type (ring_struct) oppos_ring
  type (ele_struct)     LRBBI_ele
  integer, intent(in) :: ix_LRBBI, ix_oppos 

  type (modes_struct)  mode
  type (coord_struct)  orbit_(0:n_ele_maxx)
  type (coord_struct)  oppos_orbit_(0:n_ele_maxx)
  
  real(rdef), dimension(1:n_ele_maxx) :: beta_x, beta_y, eta_x, eta_y
  real(rdef), dimension(1:n_ele_maxx) :: sigma_x, sigma_y, sigma_z

  real(rdef) :: e_spread, emit_x, emit_y, current
  integer :: i

  character*16 :: call_it

!

  call init_ele(LRBBI_ele)

  call_it = 'lrbbi'
  LRBBI_ele%name = call_it
  LRBBI_ele%key = beambeam$
  LRBBI_ele%value(charge$) = -1
  LRBBI_ele%value(n_slice$) = 1

  current = ring%param%n_part*e_charge*c_light/(ring%param%total_length)

  call twiss_at_start(ring)
  call twiss_propagate_all(ring)
  call closed_orbit_at_start(ring, orbit_(0), 4, .true.)
  call track_all(ring, orbit_)

  do i=1, ring%n_ele_ring
    beta_x(i) = ring%ele_(i)%x%beta
    beta_y(i) = ring%ele_(i)%y%beta

    eta_x(i) = ring%ele_(i)%x%eta
    eta_y(i) = ring%ele_(i)%y%eta
  enddo

! Get additional needed parameter values.

  call radiation_integrals(ring, orbit_, mode)
  e_spread = mode%sig_E

  do i=1, ring%n_ele_ring
    sigma_z(i) = mode%sig_z
  enddo

  emit_x = mode%a%emittance !a=horizontal-like. Use %b% for vertical-like.      

! Calculate sigmas for the beambeam element.

  emit_y = emit_x * 0.01

  do i=1, ring%n_ele_ring
    sigma_x(i) = sqrt(emit_x * beta_x(i) + (eta_x(i)**2) * ((e_spread)**2))
    sigma_y(i) = sqrt(emit_y * beta_y(i) + (eta_y(i)**2) *((e_spread)**2))
  enddo

  LRBBI_ele%value(sig_x$) = sigma_x(ix_LRBBI)
  LRBBI_ele%value(sig_y$) = sigma_y(ix_LRBBI)
  LRBBI_ele%value(sig_z$) = sigma_z(ix_LRBBI)

  call twiss_at_start(oppos_ring)
  call twiss_propagate_all(oppos_ring)
  call closed_orbit_at_start(oppos_ring, oppos_orbit_(0), 4, .true.)
  call track_all(oppos_ring, oppos_orbit_)

  LRBBI_ele%value(x_offset$) = oppos_orbit_(ix_oppos)%vec(1)
  LRBBI_ele%value(y_offset$) = oppos_orbit_(ix_oppos)%vec(3)

end subroutine
