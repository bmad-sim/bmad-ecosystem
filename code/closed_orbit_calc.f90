!+                           
! Subroutine closed_orbit_calc (ring, closed_orb, i_dim)
!
! Subroutine to calculate the closed orbit at the beginning of the ring.
! Closed_orbit_calc uses the 1-turn transfer matrix to converge upon a  
! solution. closed_orb(0) is used as an initial guess to the closed orbit.
!
! Note: This routine uses the 1-turn matrix ring%param%t1_mat4 or 
! ring%param%t1_mat6 in the computations. If you have changed conditions 
! significantly enough you might want to force a remake of the 1-turn matrices
! by calling clear_ring_1turn_mats.
!
! Note: See also closed_orbit_from_tracking as an alternative method
! of finding the closed orbit.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring          -- Ring_struct: Ring to track through.
!   closed_orb(0) -- Coord_struct, allocatable: Initial Guess.
!     %vec(6)          This is used as the energy around which the closed orbit
!                      is calculated if i_dim = 4.   
!   i_dim         -- Integer: Dimensions to use
!                   = 4  Transverse closed orbit at constant energy 
!                        (dE/E = closed_orb(0)%vec(6))
!                   = 6  Full closed orbit for 6x6 matrix.
!   bmad_status -- Bmad status common block
!     %exit_on_error -- If True then subroutine will terminate program
!                         if the orbit does not converge.
!     %type_out      -- If True then the subroutine will type out
!                         a warning message if the orbit does not converge.
!
! Output:
!   closed_orb(:) -- Coord_struct, allocatable: Closed orbit.
!   bmad_status   -- Bmad status common block
!     %ok          -- Set False if orbit does not converge, True otherwise.
!-

#include "CESR_platform.inc"

subroutine closed_orbit_calc (ring, closed_orb, i_dim)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct)  ring
  type (coord_struct)  del_co, del_orb
  type (coord_struct), allocatable ::  closed_orb(:)

  real(rp) s_mat(6,6), mat1(6,6), mat2(6,6), mat(6,6)
  real(rp) :: amp_co, amp_del, t1(6,6)

  integer i, j, n, n_ele, i_dim, i_max

!----------------------------------------------------------------------
! allocate storage

  call reallocate_coord (closed_orb, ring%n_ele_max)

! Error check

  bmad_status%ok = .true.

  n = i_dim

  if (i_dim == 4) then
    if (all(ring%param%t1_mat4 == 0)) &
                        call one_turn_matrix (ring, ring%param%t1_mat4)
    t1(1:4,1:4) = ring%param%t1_mat4
    closed_orb(0)%vec(5) = 0

  elseif (i_dim == 6) then
    if (all(ring%param%t1_mat6 == 0)) &
                        call one_turn_matrix (ring, ring%param%t1_mat6)
    t1 = ring%param%t1_mat6

    if (t1(6,5) == 0) then
      print *, 'ERROR IN CLOSED_ORBIT_CALC: CANNOT DO FULL 6-DIMENSIONAL'
      print *, '      CALCULATION WITH NO RF VOLTAGE!'
      bmad_status%ok = .false. 
      call err_exit
    endif

  else 
    print *, 'ERROR IN CLOSED_ORBIT_CALC: BAD I_DIM ARGUMENT:', i_dim
    bmad_status%ok = .false. 
    call err_exit
  endif
          
! init

  call mat_make_unit(s_mat(1:n,1:n))
  s_mat(2,2) = -1.0
  s_mat(4,4) = -1.0

  n_ele = ring%n_ele_use

! Turn off RF voltage if i_dim == 4 (for constant delta_E)

  if (n == 4) then
    ring%ele_(:)%internal_logic = ring%ele_(:)%is_on
    call set_on (rfcavity$, ring, .false.)
  endif

!----------------------------------------------------------------------
! d_orb = (T-1)^-1 * orbit_end

  call mat_make_unit (mat(1:n,1:n))
  mat(1:n,1:n) = mat(1:n,1:n) - t1(1:n,1:n)
  call mat_inverse(mat(1:n,1:n), mat2(1:n,1:n))

!--------------------------------------------------------------------------
! Because of nonlinearities we may need to iterate to find the solution

  i_max = 100  
  do i = 1, i_max

    call track_all (ring, closed_orb)

    del_orb%vec(1:n) = closed_orb(n_ele)%vec(1:n) - closed_orb(0)%vec(1:n)
    del_co%vec(1:n) = matmul(mat2(1:n,1:n), del_orb%vec(1:n)) 

    amp_co = sum(abs(closed_orb(0)%vec(1:n)))
    amp_del = sum(abs(del_co%vec(1:n)))
                                            
    if (amp_del < amp_co*1.0e-5 + 1.0e-8) exit

    closed_orb(0)%vec(1:n) = closed_orb(0)%vec(1:n) + del_co%vec(1:n)
            
    if (i == i_max) then
      if (bmad_status%type_out) print *,  &
          'ERROR IN CLOSED_ORBIT_CALC: NONLINEAR ORBIT NOT CONVERGING!'
      if (bmad_status%exit_on_error) call err_exit
      bmad_status%ok = .false.
    endif

  enddo

! return rf cavities to original state

  if (n == 4) ring%ele_(:)%is_on = ring%ele_(:)%internal_logic

end subroutine
