!+
! Subroutine ring_to_quad_calib (ring, cesr, k_theory, k_base, len_quad,
!                               cu_per_k_gev, quad_rot, dk_gev_dcu, cu_theory)
!
! This subroutine returns the calibration constants for the CESR quads.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ring  -- Ring_struct: Ring with lattice loaded.
!   cesr  -- Cesr_struct: Locations of the quads. 
!              Need previous call to bmad_to_cesr.
!
! Output:
! k_theory(0:120)     -- Real: Theory K of quad_i,
! k_base(0:120)       -- Real: Extrapolated K for zero cu command
! cu_per_k_gev(0:120) -- Real: CU to K*GEV calibration
! len_quad(0:120)     -- Real: Length of quad_i
! quad_rot(0:120)     -- Real: Quad rotation angle in degrees
! dk_gev_dcu(0:120)   -- Real: Derivative of K*GEV vs CU curve.
! cu_theory(0:120)    -- Integer: Scaler needed to get K_THEORY.
!
! Notes:
!
! 0) Without corrections
!         CU_THEORY = (K_THEORY - K_BASE) * GEV * CU_PER_K_GEV
!    In actuality there are corrections due to the calibration of CU with
!    current so the above equation will be off slightly.
!
! 1)      Index     Quad
!             0     REQ W
!            99     REQ E
!           101    QADD 50 W
!           102    QADD 47 W
!           103    QADD 47 E
!           104    QADD 50 E
!
! 2) Nonzero K_BASE is due to using current from the dipoles in a quadrupole.
!
! 3) DK_GEV_DCU = 1 / CU_PER_K_GEV except if there are iron saturation effects:
!
!    DK_GEV_DCU is the derivative when the actual k = K_THEORY.
!    I.e. DK_GEV_DCU is the tangent of the curve.
!
!    CU_PER_K_GEV is the average slope from K_BASE to K_THEORY.
!    I.e. CU_PER_K_GEV is the chord between 2 points on the curve.
!-

subroutine ring_to_quad_calib (ring, cesr, k_theory, k_base,  &
                 len_quad, cu_per_k_gev, quad_rot, dk_gev_dcu, cu_theory)

  use bmad_struct
  implicit none

  record /cesr_struct/ cesr
  record /ring_struct/ ring
  real energy, k_theory(0:*), k_base(0:*), len_quad(0:*)
  real cu_per_k_gev(0:*), dk_gev_dcu(0:*), quad_rot(0:*)
  integer cindex, rindex, cu_theory(0:*)

! init  &
  do cindex = 0, 120
    quad_rot(cindex) = 0.0
  enddo

! read lattice file
  energy = ring%param%energy

  do cindex = 0, 120
    rindex = cesr%quad_(cindex)%ix_ring
    if(rindex/=0) then
      k_theory(cindex) = ring%ele_(rindex)%value(k1$)
      len_quad(cindex) = ring%ele_(rindex)%value(l$)
      quad_rot(cindex) = ring%ele_(rindex)%value(tilt$)*(180./pi)
    endif
  enddo
      
! convert k_theory to scalar computer units given the specified design energy

  call k_to_quad_calib (k_theory, energy, cu_theory, k_base, dk_gev_dcu,  &
                                                                cu_per_k_gev)

end subroutine
