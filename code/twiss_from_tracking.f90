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
!   ring            -- Ring_struct: Ring to track through.
!   closed_orb_(0:) -- Coord_struct, allocatable: closed orbit.
!   d_orb(6)        -- Coord_struct: Vector of offsets to use. 
!                       All 6 components must be non-zero.
!   bmad_status      -- BMAD status common block
!         %exit_on_error -- If True then subroutine will terminate program
!                           if the orbit does not converge.
!         %type_out      -- If True then the subroutine will type out
!                           a warning message if the orbit does not converge.
!
! Output:
!   ring        -- Ring_struct: Structure holding the Twiss parameters.
!   error       -- Real(rp): A measure of how symplectic the constructed 
!                   matrices were before symplecitification. 
!                   See mat_symp_check for more details.
!   bmad_status -- BMAD status common block
!       %ok         -- Set False if orbit does not converge.
!-

#include "CESR_platform.inc"

subroutine twiss_from_tracking (ring, closed_orb_, d_orb, error)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct), intent(inout) :: ring
  type (coord_struct), intent(in), allocatable :: closed_orb_(:)
  type (coord_struct), intent(in) :: d_orb
  real(rp), intent(out) :: error

  type multi_orb_struct
    type (coord_struct), pointer :: orb(:) => null()
  end type
  type (multi_orb_struct), save :: mo(6)

  real(rp) mat(6,6), mat1(6,6), mat0(6,6), mat_inv(6,6), mat6_unit(6,6)
  integer i, j, n_ele

!

  do i=1,6
    call reallocate_coord_pointer (mo(i)%orb, ring%n_ele_max)
  enddo

! make a unit matrix

  call mat_make_unit (mat6_unit)

  n_ele = ring%n_ele_max

! Turn off RF voltage 

  ring%ele_(:)%internal_logic = ring%ele_(:)%is_on
  call set_on (rfcavity$, ring, .false.)

! track offset particles

  do i = 1, 6
    mo(i)%orb(0) = closed_orb_(0) 
    mo(i)%orb(0)%vec(i) = mo(i)%orb(0)%vec(i) + d_orb%vec(i)
    call track_all (ring, mo(i)%orb)
    if (.not. bmad_status%ok) return
  enddo

! find 1-turn matrix and twiss at starting point

  do i = 1, 6
    mat(1:6, i) = (mo(i)%orb(n_ele)%vec - closed_orb_(n_ele)%vec) / &
                                                                d_orb%vec(i)
  enddo

  call mat_symp_check (mat, error)
  call mat_symplectify (mat, ring%param%t1_mat6)
  call twiss_from_mat6 (ring%param%t1_mat6, ring%ele_(0), &
                           ring%param%stable, ring%param%growth_rate)
  ring%x%tune = ring%ele_(0)%x%phi
  ring%y%tune = ring%ele_(0)%y%phi

  if (.not. ring%param%stable) return

! compute the element transfer matrices.
! mat0 is the transfer matrix from the start to the beginning of element #j
! mat  is the transfer matrix from the start to the end of element #j
! mat1 is the transfer matrix through element #j

  call mat_make_unit (mat0)  

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

  ring%ele_(:)%is_on = ring%ele_(:)%internal_logic

end subroutine
