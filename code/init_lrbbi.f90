!+
! Subroutine init_lrbbi(ring, oppos_ring, lrbbi_ele, ix_lrbbi, ix_oppos)
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
!   LRBBI_ele       -- Ele_struct: updated lrbbi structure.
!-  

#include "CESR_platform.inc"

subroutine init_LRBBI(ring, oppos_ring, LRBBI_ele, ix_LRBBI, ix_oppos)

  use bmad_struct
  use bmad_interface
  
  implicit none

  type (ring_struct)  ring
  type (ring_struct) oppos_ring
  type (ele_struct)     LRBBI_ele
  integer, intent(in) :: ix_LRBBI, ix_oppos 

  type (modes_struct)  mode
  type (coord_struct), allocatable, save ::  orbit_(:)
  type (coord_struct), allocatable, save ::  oppos_orbit_(:)
  
  real(rp) :: beta_x, beta_y, eta_x, eta_y, sigma_z
  real(rp) :: sigma_x(ring%n_ele_max), sigma_y(ring%n_ele_max)

  real(rp) :: e_spread, emit_x, emit_y, current
  integer :: i

  character*16 :: call_it

! init

  call init_ele(LRBBI_ele)

  call_it = 'lrbbi'
  LRBBI_ele%name = call_it
  LRBBI_ele%key = beambeam$
  LRBBI_ele%value(charge$) = -1
  LRBBI_ele%value(n_slice$) = 1

  current = ring%param%n_part*e_charge*c_light/(ring%param%total_length)

  call reallocate_coord (orbit_, ring%n_ele_max)
  call reallocate_coord (oppos_orbit_, ring%n_ele_max)

!

  call twiss_at_start(ring)
  call twiss_propagate_all(ring)
  call closed_orbit_at_start(ring, orbit_(0), 4, .true.)
  call track_all(ring, orbit_)

! Get additional needed parameter values.

  call radiation_integrals(ring, orbit_, mode)

  e_spread = mode%sigE_E
  sigma_z = mode%sig_z
  emit_x = mode%a%emittance !a=horizontal-like. Use %b% for vertical-like.
  emit_y = emit_x * 0.01

! Calculate sigmas for the beambeam element.

  do i=1, ring%n_ele_ring
    beta_x = ring%ele_(i)%x%beta
    beta_y = ring%ele_(i)%y%beta
    eta_x = ring%ele_(i)%x%eta
    eta_y = ring%ele_(i)%y%eta
    sigma_x(i) = sqrt(emit_x * beta_x + (eta_x**2) * ((e_spread)**2))
    sigma_y(i) = sqrt(emit_y * beta_y + (eta_y**2) *((e_spread)**2))
  enddo

  LRBBI_ele%value(sig_x$) = sigma_x(ix_LRBBI)
  LRBBI_ele%value(sig_y$) = sigma_y(ix_LRBBI)
  LRBBI_ele%value(sig_z$) = sigma_z

  call twiss_at_start(oppos_ring)
  call twiss_propagate_all(oppos_ring)
  call closed_orbit_at_start(oppos_ring, oppos_orbit_(0), 4, .true.)
  call track_all (oppos_ring, oppos_orbit_)

  LRBBI_ele%value(x_offset$) = oppos_orbit_(ix_oppos)%vec(1)
  LRBBI_ele%value(y_offset$) = oppos_orbit_(ix_oppos)%vec(3)

end subroutine
