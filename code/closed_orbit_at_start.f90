!+                           
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
! Note: See also closed_orbit_from_tracking as an alternative method
! of finding the closed orbit.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ring   -- Ring_struct: Ring
!     %ele_(0).mat6      -- This must contain the 6x6 1-turn matrix
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
!
! Note:
!     1) You can use TWISS_AT_START if I_DIM = 4 to calculate the 1-turn
!         matrix in RING%ELE_(0)%MAT6. NOTE: This does not work with I_DIM = 6.
!     2) Use CLOSED_ORBIT_AT_START to find the beginning closed orbit.
!     3) Use TRACK_ALL to propagate the closed orbit through the rest
!         of the ring.
!-

!$Id$
!$Log$
!Revision 1.4  2002/01/16 21:04:17  helms
!Fixed problem with passing optional arguments.
!
!Revision 1.3  2002/01/08 21:44:38  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:49  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine closed_orbit_at_start (ring, co, i_dim, iterate)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct)  ring
  type (coord_struct)  orbit_end, del_co, del_orb
  type (coord_struct)  orbit_end_e_(0:n_ele_maxx), co
  type (coord_struct)  orbit_(0:n_ele_maxx)

  real s_mat(6,6), mat1(6,6), mat2(6,6), mat(6,6)
  real amp_co, amp_del, factor / 1.0 /

  integer i, j, n, n_ele, i_dim

  logical iterate, is_on(0:n_ele_maxx)

! Error check

  bmad_status%ok = .true.

  if (ring%ele_(0)%mat6(6,5) == 0 .and. i_dim == 6) then
    type *, 'ERROR IN CLOSED_ORBIT_AT_START: CANNOT DO FULL 6-DIMENSIONAL'
    type *, '      CALCULATION WITH NO RF VOLTAGE!'
    bmad_status%ok = .false. 
    call err_exit
  endif

  if (i_dim /= 4 .and. i_dim /= 6) then
    type *, 'ERROR IN CLOSED_ORBIT_AT_START: BAD I_DIM ARGUMENT:', i_dim
    bmad_status%ok = .false. 
    call err_exit
  endif
          
! init

  n = i_dim
  call mat_unit(s_mat, n, 6)
  s_mat(2,2) = -1.0
  s_mat(4,4) = -1.0

  orbit_(0)%vec(:) = 0.0
  orbit_(0)%vec(6) = co%vec(6)

  orbit_end_e_(0)%vec(:) = 0.0
  orbit_end_e_(0)%vec(6) = co%vec(6)

  n_ele = ring%n_ele_use

! Turn off RF voltage if i_dim == 4 (for constant delta_E)

  if (n == 4) then
    is_on = .false.
    do i = 1, n_ele
      if (ring%ele_(i)%key == rfcavity$) then
        is_on(i) = ring%ele_(i)%is_on
        ring%ele_(i)%is_on = .false.
      endif
    enddo
  endif

!---------------------------------------------------------------------
! with e/w symmetry

  if (ring%param%symmetry == ew_antisymmetry$) then

    call track_all (ring, orbit_)

     do i = 1, ring%n_ele_use
       if (ring%ele_(i)%key == elseparator$) then ! set kicks as in east
         if (ring%ele_(i)%type(1:4) == 'ANTI') then
           ring%ele_(i)%value(hkick$) = -ring%ele_(i)%value(hkick$)
           ring%ele_(i)%value(vkick$) = -ring%ele_(i)%value(vkick$)
         endif
       endif !separators
     end do

     call track_all (ring, orbit_end_e_)

     do i = 1, ring%n_ele_use  ! set separator kicks back to values in west
       if (ring%ele_(i)%key == elseparator$) then
         if (ring%ele_(i)%type(1:4) == 'ANTI') then
           ring%ele_(i)%value(hkick$) = -ring%ele_(i)%value(hkick$)
           ring%ele_(i)%value(vkick$) = -ring%ele_(i)%value(vkick$)
         endif
       endif !separators
     end do

! symmetric ring then
!       X = (Mat*S_Mat-S_Mat*Mat)^-1 * (S_mat * orbit_w - orbit_e)

    orbit_(n_ele)%x%vel = -orbit_(n_ele)%x%vel
    orbit_(n_ele)%y%vel = -orbit_(n_ele)%y%vel
    orbit_end%vec = orbit_(n_ele)%vec - orbit_end_e_(n_ele)%vec

    mat(1:n,1:n) = matmul(ring%ele_(0)%mat6(1:n,1:n), s_mat(1:n,1:n)) - &
                  matmul (s_mat(1:n,1:n), ring%ele_(0)%mat6(1:n,1:n))
                   
!----------------------------------------------------------------------
! No symmetry: X = (Mat-1)^-1 * orbit_end

  else
    call track_all (ring, orbit_)
    orbit_end = orbit_(n_ele)
    call mat_unit (mat, n, 6)
    mat(1:n,1:n) = mat(1:n,1:n) - ring%ele_(0)%mat6(1:n,1:n)
  endif

!----------------------------------------------------------------------
! Now find closed orbit at start

  call mat_inverse(mat(1:n,1:n), mat2(1:n,1:n))
  co%vec(1:n) = matmul(mat2(1:n,1:n), orbit_end%vec(1:n))

  if (n == 4) co%z%pos = 0

!--------------------------------------------------------------------------
! Because of nonlinearities we may need to iterate to find the solution
!  But iteration makes no sense in the case of ew_antisymmetry
         
  if (iterate .and. ring%param%symmetry /= ew_antisymmetry$) then
                                                   
    do i = 1, 100

      orbit_(0) = co
      call track_all (ring, orbit_)

      del_orb%vec(1:n) = orbit_(n_ele)%vec(1:n) - co%vec(1:n)
      del_co%vec(1:n) = matmul(mat2(1:n,1:n), del_orb%vec(1:n)) 
      co%vec(1:n) = co%vec(1:n) + factor * del_co%vec(1:n)

      amp_co = sum(abs(co%vec(1:n)))
      amp_del = sum(abs(del_co%vec(1:n)))
                                            
      if (amp_del < amp_co*1.0e-5 + 1.0e-8) exit

    enddo
                              
    if (i == 101) then
      if (bmad_status%type_out) type *,  &
          'ERROR IN CLOSED_ORBIT_AT_START: NONLINEAR ORBIT NOT CONVERGING!'
      if (bmad_status%exit_on_error) call err_exit
      bmad_status%ok = .false.
    endif

  endif

! return rf cavities to original state

  if (n == 4) then
    where (is_on(1:n_ele)) ring%ele_(1:n_ele)%is_on = .true.
  endif

end subroutine
