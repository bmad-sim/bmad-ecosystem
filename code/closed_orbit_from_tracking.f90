!+
! Subroutine closed_orbit_from_tracking (ring, closed_orb_, i_dim, 
!                                                eps_rel, eps_abs, init_guess)
!
! Subroutine to find the closed orbit via tracking. 
!
! See closed_orbit_at_start for an alternative routine.
! This routine has the advantage of not using the 1-turn transfer matrix 
! (which is calculated from the transfer matrices of the elements).
! This routine is probably slower than closed_orbit_at_start (especially for 
! linear lattices) but (hopefully) will be better behaved when there
! are strong nonlinearities. Possible convergence porblems with this 
! routine occurs when the initial guess is far from the true closed orbit.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring     -- Ring_struct: Ring to track through.
!   i_dim    -- Integer: 
!       = 2,4  Transverse closed orbit at constant energy (dE/E = CO.Z.VEL)
!       = 6    Full closed orbit using the entire transfer 6x6 matrix.
!   eps_rel(:) -- Real(rdef): relative allowed error.
!   eps_abs(:) -- Real(rdef): absolute allowed error.
!   init_guess -- [Optional] Coord_struct: Starting guess for the closed 
!                orbit at the start of the ring. If not present then
!                the origin will be used. 
!   bmad_status -- BMAD status common block:
!         %exit_on_error -- If True then subroutine will terminate program
!                           if the orbit does not converge.
!         %type_out      -- If True then the subroutine will type out
!                           a warning message if the orbit does not converge.
!
!
! Output:
!   closed_orb_(0:n_ele_maxx) -- Coord_struct: closed orbit.
!   bmad_status -- BMAD status common block:
!       %ok         -- Set False if orbit does not converge.
!-

!$Id$
!$Log$
!Revision 1.11  2003/01/27 14:40:32  dcs
!bmad_version = 56
!
!Revision 1.10  2002/12/13 16:23:32  dcs
!*** empty log message ***
!
!Revision 1.8  2002/11/14 16:07:39  dcs
!Fixed bug to overwrite closed_orb_(0)%vec(5) when i_dim = 4
!
!Revision 1.7  2002/11/06 06:48:31  dcs
!Changed arg array
!
!Revision 1.5  2002/06/13 14:54:23  dcs
!Interfaced with FPP/PTC
!
!Revision 1.4  2002/02/23 20:32:13  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/16 21:04:17  helms
!Fixed problem with passing optional arguments.
!
!Revision 1.2  2001/09/27 18:31:49  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc" 

subroutine closed_orbit_from_tracking (ring, closed_orb_, i_dim, &
                                                 eps_rel, eps_abs, init_guess)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) ring
  type (coord_struct), intent(inout) :: closed_orb_(0:)
  type (coord_struct) :: orb_(0:n_ele_maxx), start(100), end(100)
  type (coord_struct), optional :: init_guess

  real(rdef), intent(in) :: eps_rel(:), eps_abs(:)

  real(rdef) amp_co, amp_del, orb_diff(6), amp(6)
  real(rdef) mat6(6,6), mat6_inv(6,6), mat6_unit(6,6)
  real(rdef) start_mat(6,6), end_mat(6,6)
  real(rdef) :: error

  integer i_dim, i, i1, i2, j, k, jmax, n_ele, j0, jj, nd, nnd, msk(6)

  logical is_on(0:n_ele_maxx)
  logical :: debug = .false.

! make a unit matrix

  call mat_unit (mat6_unit, 6, 6)

  n_ele = ring%n_ele_ring
  nd = i_dim

! Turn off RF voltage if i_dim == 4 (for constant delta_E)

  if (nd ==2 .or. nd == 4) then
    is_on = .false.
    do i = 1, n_ele
      if (ring%ele_(i)%key == rfcavity$) then
        is_on(i) = ring%ele_(i)%is_on
        ring%ele_(i)%is_on = .false.
      endif
    enddo
  elseif (nd /= 6) then
    print *, 'ERROR IN CLOSED_ORBIT_FROM_TRACKING: BAD "I_DIM":', nd
    call err_exit    
  endif

! determine initial guess

  if (present(init_guess)) then
    closed_orb_(0) = init_guess
  else
    closed_orb_(0)%vec = 0
  endif

!-------------------
! iterate until we find the closed orbit

  jmax = 100
  do j = 1, jmax

! track the current guess

    call track_all (ring, closed_orb_)
    if (ring%param%lost) exit

! save start and end coords.

    start(j)%vec = closed_orb_(0)%vec
    end(j)%vec   = closed_orb_(n_ele)%vec

! is this good enough? if so return.

    orb_diff = closed_orb_(n_ele)%vec - closed_orb_(0)%vec
    amp = max(abs(closed_orb_(0)%vec), abs(closed_orb_(n_ele)%vec))

    if (all( (abs(orb_diff(1:nd)) < eps_abs(1:nd)) .or. &
             (abs(orb_diff(1:nd)) < eps_rel(1:nd) * amp(1:nd)) )) then
      if (nd == 2 .or. nd == 4) where (is_on(1:n_ele)) ring%ele_(1:n_ele)%is_on = .true.
      return
    endif

! if we have enough data then form the transfer matrix and invert to get 
! the new closed orbit guess. 
! do not use the matrix if it is very non-symplectic.
!
! msk is used to mask out any plans where closed_orb_ has converged to
! the closed orbit. If we do not do this then the computation can blow up.

  if (j > nd) then

    nnd = 0  ! number of unmasked dimensions
    do i = 2, nd, 2
      i1 = i-1
      i2 = i
      if (.not. all( (abs(orb_diff(i1:i2)) < eps_abs(i1:i2)) .or. &
                   (abs(orb_diff(i1:i2)) < eps_rel(i1:i2) * amp(i1:i2)) )) then
        msk(nnd+1:nnd+2) = (/ i1, i2 /)
        nnd = nnd + 2
      endif
     enddo

    do i = 1, nd
      j0 = j - nd
      jj = j0 + i
      end_mat(1:nnd,i) = end(jj)%vec(msk(1:nnd)) - end(j0)%vec(msk(1:nnd))
      start_mat(1:nnd, i) = &
                    start(jj)%vec(msk(1:nnd)) - start(j0)%vec(msk(1:nnd))
    enddo

    call mat_det(start_mat(1:nnd,1:nnd), error, nnd, nnd)
    if (debug) print *, 'det:', error
    call mat_inverse (start_mat(1:nnd,1:nnd), start_mat(1:nnd, 1:nnd))
    mat6 = matmul(end_mat(1:nnd,1:nnd), start_mat(1:nnd,1:nnd))
    call mat_symp_check (mat6(1:nnd,1:nnd), error)
    if (debug) print *, 'error:', error
    if (debug) print '(i4, a, 3p6f11.5)', j, ':', (closed_orb_(0)%vec(i), i = 1, nd)

    if (error < 0.1) then  ! if error is low use matrix
      mat6(1:nnd,1:nnd) = mat6_unit(1:nnd,1:nnd) - mat6(1:nnd,1:nnd)
      call mat_inverse (mat6(1:nnd,1:nnd), mat6_inv(1:nnd,1:nnd))
      closed_orb_(0)%vec(msk(1:nnd)) = closed_orb_(0)%vec(msk(1:nnd)) + &
                        matmul(mat6_inv(1:nnd,1:nnd), orb_diff(msk(1:nnd)))
      cycle
    endif

  endif


! if we are here then we did not make a guess using the matrix.
! The new guess is the average of the start and end orbits.

    closed_orb_(0)%vec(1:nd) = &
              (closed_orb_(0)%vec(1:nd) + closed_orb_(n_ele)%vec(1:nd)) / 2

  enddo

! error

  if (bmad_status%type_out) then
    if (j == jmax+1) then
      print *, 'ERROR IN CLOSED_ORBIT_FROM_TRACKING: ORBIT NOT CONVERGING!'
    else
      print *, 'ERROR IN CLOSED_ORBIT_FROM_TRACKING: PROBLEM TRACKING ORBIT.'
    endif
  endif

  if (bmad_status%exit_on_error) call err_exit
  bmad_status%ok = .false.

end subroutine
