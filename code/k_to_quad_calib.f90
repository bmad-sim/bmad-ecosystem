!+
! Subroutine k_to_quad_calib (k_theory, energy, cu_theory, k_base,
!                                                  dk_gev_dcu, cu_per_k_gev)
!
! This subroutine returns the calibration constants for the CESR quads.
! See also: QUAD_CALIB.
!
! Input:
!   energy          -- Real: Energy in GeV
!   k_theory(0:120) -- Real: Theory K of quad_i,
!
! Output:
!   k_base(0:120)       -- Real: Extrapolated K for zero CU command
!   cu_per_k_gev(0:120) -- Real: CU to K*GEV calibration
!   dk_gev_dcu(0:120)   -- Real: Derivative of K*GEV vs CU curve.
!   cu_theory(0:120)    -- Integer: Scaler needed to get K_THEORY.
!
! Notes:
!                                           
! 1)      Index     Quad                    
!             0     REQ W
!            99     REQ E
!           101    QADD 49 W
!           102    QADD ...
!
! 2) Nonzero K_BASE is due to using current from the dipoles in a quadrupole.
!
! 3) DK_GEV_DCU = 1 / CU_PER_K_GEV except if there are iron saturation effects:
!
!    DK_GEV_DCU is the derivative when the k = K_THEORY.
!    I.e. DK_GEV_DCU is the tangent of the curve:
!         dK = dCU * (dK_GeV_dCU / GeV)
!
!    CU_PER_K_GEV is the average slope from K_BASE to K_THEORY.
!    I.e. CU_PER_K_GEV is the chord between 2 points on the curve:
!         CU_THEORY = (K_THEORY - K_BASE) * GEV * CU_PER_K_GEV
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:52  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine k_to_quad_calib (k_theory, energy, cu_theory, k_base,  &
                                                   dk_gev_dcu, cu_per_k_gev)

  implicit none

  real energy, theory
  real k_theory(0:*), dk_gev_dcu(0:*)
  real cu_per_k_gev(0:120), k_base(0:120)
  real satfac(0:120), gain_(0:120), offset_(0:120)

  integer lun, lunget, i, cu_theory(0:*)

  namelist / quadrupole_cal / cu_per_k_gev, k_base

  namelist / quad_scaler_corrections / gain_, offset_

  character*60 file_scaler, file_cal
  	
! init

  cu_per_k_gev = 0
  k_base       = 0.0
  gain_        = 1.0

! Current to CU corrections

  lun = lunget()
  file_scaler = 'CESR_CONST:QUADRUPOLE_SCALER.CAL'
  file_scaler = FullFileName(file_scaler)
  open (lun, file = file_scaler,  &
                             status = 'old', readonly, shared, err = 9200)
  read (lun, nml=quad_scaler_corrections, err = 9210)
  close (lun)

! setup calibrations with constants read from QUAD.CALIBRATION
! A non-zero k_base is due to dipole current through the magnet.
! Thus k_base scales with energy.

  lun = lunget()
  file_cal = 'CESR_CONST:QUADRUPOLE.CAL'
  file_cal = FullFileName(file_cal)
  open (lun, file = file_cal, status = 'old', readonly, shared, err = 9000)
  read (lun, nml=quadrupole_cal, err = 9010)
  close (lun)

  do i = 0, 120
    if(cu_per_k_gev(i) /= 0) dk_gev_dcu(i) = 1 / cu_per_k_gev(i)
  enddo

! flip signs for vertically focusing quads

  do i = 1, 120
    if (k_theory(i) < 0) then
      k_base(i)       = -k_base(i)
      cu_per_k_gev(i) = -cu_per_k_gev(i)
      dk_gev_dcu(i)   = -dk_gev_dcu(i)
    endif
  enddo

! Give CU values with current calibration corrections

  do i = 0, 120
    if (gain_(i) /= 0)  then
      theory = (k_theory(i)-k_base(i)) * energy * cu_per_k_gev(i)
      cu_theory(i) = nint((theory - offset_(i)) / gain_(i))
      dk_gev_dcu(i) = dk_gev_dcu(i) * gain_(i)
      cu_per_k_gev(i) = cu_per_k_gev(i) / gain_(i)
      if (cu_per_k_gev(i) /= 0) k_base(i) = &
                    k_theory(i) - cu_theory(i) / (energy * cu_per_k_gev(i))
    endif
  enddo

  return

! For errors...

9000  type *
  type *, 'ERROR IN K_TO_QUAD_CALIB: UNABLE TO OPEN PARAMETER FILE:'
  type *, '      ', trim(file_cal)
  call err_exit

9010  type *
  type *, 'ERROR IN K_TO_QUAD_CALIB: UNABLE TO READ PARAMETERS IN FILE:'
  type *, '      ', trim(file_cal)
  rewind (lun)
  read (lun, nml=quadrupole_cal)
  call err_exit

9200  type *
  type *, 'ERROR IN K_TO_QUAD_CALIB: UNABLE TO OPEN CORRECTION FILE:'
  type *, '     ', trim(file_scaler)
  call err_exit

9210  type *
  type *, 'ERROR IN K_TO_QUAD_CALIB: UNABLE TO READ PARAMETERS IN FILE:'
  type *, '     ', trim(file_scaler)
  call err_exit

end subroutine
              
