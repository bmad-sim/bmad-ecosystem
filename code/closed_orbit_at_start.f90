!+   
! ***********************************************************************
! ***********************************************************************
! ***************                                         ***************
! *************** OBSOLETE: USE CLOSED_ORBIT_CALC INSTEAD ***************
! ***************                                         ***************
! ***********************************************************************
! ***********************************************************************
!                        
! Subroutine closed_orbit_at_start (ring, co, i_dim, iterate)
!
! Subroutine to calculate the closed orbit at the beginning of the ring.
! The Subroutine works by first guessing the closed orbit based upon 
! the 1-turn matrix and tracking for 1-turn.
!
! If ITERATE is True the subroutine iteratively uses the current guess
! for the closed orbit to converge to the correct answer. If False
! then the first guess is returned as the closed orbit. ITERATE should
! always be set True unless you *know* you have a linear lattice and you
! need to save time.
!
! Note: This routine uses the 1-turn matrix ring%param%t1_no_RF or 
! ring%param%t1_with_RF in the computations. If you have changed conditions 
! significantly enough you might want to force a remake of the 1-turn matrices
! by calling clear_ring_1turn_mats.
!
! Note: See also closed_orbit_from_tracking as an alternative method
! of finding the closed orbit.
!
! Note: Use TRACK_ALL to propagate the closed orbit through the rest
! of the ring.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring   -- Ring_struct: Ring to track through.
!   i_dim  -- Integer: Dimensions to use
!     = 4  Transverse closed orbit at constant energy (dE/E = CO.Z.VEL)
!     = 6  Full closed orbit for 6x6 matrix.
!   iterate -- Logical: See above for details. This should always be set
!              True unless you know what you are doing.  
!   bmad_status -- Bmad status common block
!         %exit_on_error -- If True then subroutine will terminate program
!                           if the orbit does not converge.
!         %type_out      -- If True then the subroutine will type out
!                           a warning message if the orbit does not converge.
!
! Output:
!   co          -- Coord_struct: Closed orbit at the starting point.
!   bmad_status -- Bmad status common block
!       %ok         -- Set False if orbit does not converge, True otherwise.
!-

#include "CESR_platform.inc"

subroutine closed_orbit_at_start (ring, co, i_dim, iterate)

  use bmad_struct
  use bmad_interface, except => closed_orbit_at_start
  use bookkeeper_mod, only: set_on_off, save_state$, restore_state$, off$

  implicit none

  type (ring_struct)  ring
  type (coord_struct)  orbit_end, del_co, del_orb
  type (coord_struct)  co
  type (coord_struct), allocatable, save :: orbit_(:)

  real(rp) s_mat(6,6), mat2(6,6), mat(6,6)
  real(rp) :: amp_co, amp_del, factor = 1.0, t1(6,6)

  integer i, n, n_ele, i_dim

  logical iterate
  logical fluct_saved, aperture_saved

! Init.
! allocate storage
! Random fluctuations must be turned off to find the closed orbit.

  call reallocate_coord (orbit_, ring%n_ele_max)

  fluct_saved = bmad_com%radiation_fluctuations_on
  bmad_com%radiation_fluctuations_on = .false.  

  aperture_saved = ring%param%aperture_limit_on
  ring%param%aperture_limit_on = .false.

! Error check

  bmad_status%ok = .true.

  n = i_dim

  if (i_dim == 4) then
    if (all(ring%param%t1_no_RF == 0)) &
                call transfer_matrix_calc (ring, .false., ring%param%t1_no_RF)
    t1 = ring%param%t1_no_RF

  elseif (i_dim == 6) then
    if (all(ring%param%t1_with_RF == 0)) &
                call transfer_matrix_calc (ring, .true., ring%param%t1_with_RF)
    t1 = ring%param%t1_with_RF

    if (t1(6,5) == 0) then
      print *, 'ERROR IN CLOSED_ORBIT_AT_START: CANNOT DO FULL 6-DIMENSIONAL'
      print *, '      CALCULATION WITH NO RF VOLTAGE!'
      bmad_status%ok = .false. 
      call err_exit
    endif

  else 
    print *, 'ERROR IN CLOSED_ORBIT_AT_START: BAD I_DIM ARGUMENT:', i_dim
    bmad_status%ok = .false. 
    call err_exit
  endif
          
! init

  call mat_make_unit(s_mat(1:n,1:n))
  s_mat(2,2) = -1.0
  s_mat(4,4) = -1.0

  orbit_(0)%vec(:) = 0.0
  orbit_(0)%vec(6) = co%vec(6)

  n_ele = ring%n_ele_use

! Turn off RF voltage if i_dim == 4 (for constant delta_E)

  if (n == 4) then
    call set_on_off (rfcavity$, ring, save_state$)
    call set_on_off (rfcavity$, ring, off$)
  endif

!---------------------------------------------------------------------
!----------------------------------------------------------------------
! X = (Mat-1)^-1 * orbit_end

  call track_all (ring, orbit_)
  orbit_end = orbit_(n_ele)
  call mat_make_unit (mat(1:n,1:n))
  mat(1:n,1:n) = mat(1:n,1:n) - t1(1:n,1:n)

!----------------------------------------------------------------------
! Now find closed orbit at start

  call mat_inverse(mat(1:n,1:n), mat2(1:n,1:n))
  co%vec(1:n) = matmul(mat2(1:n,1:n), orbit_end%vec(1:n))

  if (n == 4) co%vec(5) = 0

!--------------------------------------------------------------------------
! Because of nonlinearities we may need to iterate to find the solution
         
  if (iterate) then
                                                   
    do i = 1, 100

      orbit_(0) = co
      call track_all (ring, orbit_)

      del_orb%vec(1:n) = orbit_(n_ele)%vec(1:n) - co%vec(1:n)
      del_co%vec(1:n) = matmul(mat2(1:n,1:n), del_orb%vec(1:n)) 
      co%vec(1:n) = co%vec(1:n) + factor * del_co%vec(1:n)

      amp_co = sum(abs(co%vec(1:n)))
      amp_del = sum(abs(del_co%vec(1:n)))

      if (amp_del < amp_co * bmad_com%rel_tolerance + bmad_com%abs_tolerance) exit

    enddo
                              
    if (i == 101) then
      if (bmad_status%type_out) print *,  &
          'ERROR IN CLOSED_ORBIT_AT_START: NONLINEAR ORBIT NOT CONVERGING!'
      if (bmad_status%exit_on_error) call err_exit
      bmad_status%ok = .false.
    endif

  endif

! return rf cavities to original state

  if (n == 4) call set_on_off (rfcavity$, ring, restore_state$)

  bmad_com%radiation_fluctuations_on = fluct_saved  ! restore state
  ring%param%aperture_limit_on = aperture_saved

end subroutine




