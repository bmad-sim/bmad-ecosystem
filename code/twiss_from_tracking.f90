!+
! Subroutine twiss_from_tracking (ring, closed_orb_, d_orb, error)
!
! Subroutine to compute from tracking, for every element in the ring, 
! the Twiss parameters and the transfer matrices. 
! Because of nonlinearities, and the fact that d_orb is finite,
! the computed 6x6 matrices will not be exactly symplectic. To handle this
! these matrices are symplecitified.
!     
! See twiss_at_start and twiss_propagate_all for alternative routines.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring    -- Ring_struct: Ring to track through.
!   closed_orb_(0:n_ele_maxx) -- Coord_struct: closed orbit.
!   d_orb(6) -- Coord_struct: Vector of offsets to use. All 6 components 
!                must be non-zero.
!   bmad_status -- BMAD status common block
!         %exit_on_error -- If True then subroutine will terminate program
!                           if the orbit does not converge.
!         %type_out      -- If True then the subroutine will type out
!                           a warning message if the orbit does not converge.
!
! Output:
!   error   -- Real(rdef): A measure of how symplectic the constructed 
!                                                           matrices were.
!              before symplecitification. See mat_symp_check for more details.
!   ring    -- Ring_struct:
!     ele_(:)%mat6  -- 6x6 transfer matrices. The for
!     %x%beta, etc  -- Twiss parameters.
!   bmad_status  -- BMAD status common block
!       %ok         -- Set False if orbit does not converge.
!-

!$Id$
!$Log$
!Revision 1.7  2002/11/06 06:48:32  dcs
!Changed arg array
!
!Revision 1.6  2002/08/20 20:34:55  dcs
!symp_lie_bmad / symp_lie_ptc added
!
!Revision 1.5  2002/07/16 21:33:58  dcs
!*** empty log message ***
!
!Revision 1.4  2002/02/23 20:32:29  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2001/10/02 18:49:13  rwh24
!More compatibility updates; also added many explicit variable declarations.
!
!Revision 1.2  2001/09/27 18:32:00  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine twiss_from_tracking (ring, closed_orb_, d_orb, error)

  use bmad

  implicit none

  type (ring_struct), intent(inout) :: ring
  type (coord_struct), intent(in) :: closed_orb_(0:), d_orb
  real(rdef), intent(out) :: error

  type multi_orb_struct
    type (coord_struct) orb(0:n_ele_maxx)
  end type
  type (multi_orb_struct) mo(6)

  real(rdef) mat(6,6), mat1(6,6), mat0(6,6), mat_inv(6,6), mat6_unit(6,6)
  integer i, j, n_ele

  logical is_on(0:n_ele_maxx)

! make a unit matrix

  call mat_unit (mat6_unit, 6, 6)

  n_ele = ring%n_ele_max

! Turn off RF voltage 

  is_on = .false.
  do i = 1, n_ele
    if (ring%ele_(i)%key == rfcavity$) then
      is_on(i) = ring%ele_(i)%is_on
      ring%ele_(i)%is_on = .false.
    endif
  enddo

! track offset particles

  do i = 1, 6
    mo(i)%orb(0) = closed_orb_(0) 
    mo(i)%orb(0)%vec(i) = mo(i)%orb(0)%vec(i) + d_orb%vec(i)
    call track_all (ring, mo(i)%orb)
    if (.not. bmad_status%ok) return
  enddo

! find 1-turn matrix and twiss at starting point

  do i = 1, 6
    mat(1:6, i) = (mo(i)%orb(n_ele)%vec - closed_orb_(n_ele)%vec) / d_orb%vec(i)
  enddo

  call mat_symp_check (mat, error)
  call mat_symplectify (mat, ring%ele_(0)%mat6)
  call twiss_from_mat6 (ring%ele_(0)%mat6, ring%ele_(0), &
                           ring%param%stable, ring%param%growth_rate)
  ring%x%tune = ring%ele_(0)%x%phi
  ring%y%tune = ring%ele_(0)%y%phi

  if (.not. ring%param%stable) return

! compute the element transfer matrices.
! mat0 is the transfer matrix from the start to the beginning of element #j
! mat  is the transfer matrix from the start to the end of element #j
! mat1 is the transfer matrix through element #j

  call mat_unit (mat0, 6, 6)  

  do j = 1, n_ele

    do i = 1, 6
      mat(1:6, i) = (mo(i)%orb(j)%vec - closed_orb_(j)%vec) / d_orb%vec(i)
    enddo

    call mat_symp_conj (mat0, mat_inv, 6, 6)   ! symp_conj is the inverse
    mat1 = matmul (mat, mat_inv)
    call mat_symplectify (mat1, ring%ele_(j)%mat6)
    mat0 = matmul (ring%ele_(j)%mat6, mat0)

  enddo

! now compute the twiss parameters

  call twiss_propagate_all (ring)

! turn RF back on

  where (is_on(1:n_ele)) ring%ele_(1:n_ele)%is_on = .true.

end subroutine
