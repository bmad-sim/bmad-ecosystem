!+
! Subroutine twiss_from_tracking (ring, ref_orb0, error, d_orb)
!
! Subroutine to compute, from tracking, and for every element in the ring, 
! the Twiss parameters and the transfer matrices about some reference orbit.
!
! The is done by tracking 12 particle offset from the reference orbit 
! by +/- d_orb.
!
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
!   ref_orb0        -- Coord_struct: Reference orbit at ring%ele_(0).
!   d_orb(6)        -- real(rp), optional: Vector of offsets to use. 
!                       If not present or zero bmad_com%d_orb(:) will be used.
!   bmad_status      -- BMAD status common block
!     %type_out      -- If True then the subroutine will type out
!                         a warning message if a particle is lost in tracking.
!
! Output:
!   ring        -- Ring_struct: Structure holding the Twiss parameters.
!     %mat6       -- transfer matrices.
!
!   error       -- Real(rp): A measure of how symplectic the constructed 
!                   matrices were before symplecitification. 
!                   See mat_symp_check for more details.
!   bmad_status -- BMAD status common block
!       %ok         -- Set False if particle lost while tracking.
!-

#include "CESR_platform.inc"

subroutine twiss_from_tracking (ring, ref_orb0, error, d_orb)

  use bmad_struct
  use bmad_interface
  use bookkeeper_mod

  implicit none

  type (ring_struct), intent(inout) :: ring
  type (coord_struct), intent(in) :: ref_orb0

  real(rp), optional, intent(in) :: d_orb(:)
  real(rp), intent(out) :: error
  real(rp) delta(0:6)

  type multi_orb_struct
    type (coord_struct), allocatable :: orb(:) 
  end type
  type (multi_orb_struct), save :: mo(-6:6)

  real(rp) mat(6,6), mat1(6,6), mat0(6,6), mat_inv(6,6), mat6_unit(6,6)
  integer i, j, n_ele

! Init

  do i = -6, 6
    call reallocate_coord (mo(i)%orb, ring%n_ele_max)
  enddo

  call mat_make_unit (mat6_unit)
  n_ele = ring%n_ele_ring

  call set_on_off (rfcavity$, ring, save_state$)
  call set_on_off (rfcavity$, ring, off$)

  delta(0) = 0
  delta(1:6) = bmad_com%d_orb(1:6)
  if (present(d_orb)) where (d_orb(1:6) /= 0) delta(1:6) = d_orb(1:6)

! track offset particles

  do i = 1, 6
    mo(i)%orb(0) = ref_orb0 
    mo(i)%orb(0)%vec(i) = ref_orb0%vec(i) + delta(i)
    call track_all (ring, mo(i)%orb)
    if (.not. bmad_status%ok) return
    mo(-i)%orb(0) = ref_orb0 
    mo(-i)%orb(0)%vec(i) = ref_orb0%vec(i) - delta(i)
    call track_all (ring, mo(-i)%orb)
    if (.not. bmad_status%ok) return
    mat(:,i) = (mo(i)%orb(n_ele)%vec - mo(-i)%orb(n_ele)%vec) / (2*delta(i))
  enddo

! find 1-turn matrix and twiss at starting point

  call mat_symp_check (mat, error)
  call mat_symplectify (mat, mat)
  call twiss_from_mat6 (mat, ring%ele_(0), &
                           ring%param%stable, ring%param%growth_rate)
  ring%x%tune = ring%ele_(0)%x%phi
  ring%y%tune = ring%ele_(0)%y%phi
  ring%param%t1_no_RF = mat

  if (.not. ring%param%stable) return

! compute the element transfer matrices.
! mat0 is the transfer matrix from the start to the beginning of element #j
! mat  is the transfer matrix from the start to the end of element #j
! mat1 is the transfer matrix through element #j

  call mat_make_unit (mat0)  

  do j = 1, n_ele

    do i = 1, 6
      mat(:,i) = (mo(i)%orb(j)%vec - mo(-i)%orb(j)%vec) / (2*delta(i))
    enddo

    call mat_symp_conj (mat0, mat_inv)   ! symp_conj is the inverse
    mat1 = matmul (mat, mat_inv)
    call mat_symplectify (mat1, ring%ele_(j)%mat6)
    mat0 = matmul (ring%ele_(j)%mat6, mat0)

  enddo

! Now compute the twiss parameters.
! And turn the RF back on.

  call twiss_propagate_all (ring)
  call set_on_off (rfcavity$, ring, restore_state$)

end subroutine
