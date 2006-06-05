!+
! Subroutine twiss_at_start (ring)
!
! Subroutine to calculate, for a circular machine, the closed 1-turn 
! solution for the Twiss parameters at the start of the ring.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring   -- Ring_struct: Ring
!   bmad_status -- BMAD Common block status structure
!     %type_out  -- Logical: If .true. then will type a message for
!                       non ok$ STATUS
!
! Output:
!   ring
!     %param%t1_no_RF --  Note: Only the linear part is computed.
!     %ele_(0)%x      -- "a" mode Twiss parameters at the start of the ring.
!     %ele_(0)%y      -- "b" mode Twiss parameters at the start of the ring.
!     %ele_(0)%c_mat  -- Coupling matrix.
!     %x%tune         -- Fractional part of the tune in radians
!     %y%tune         -- Fractional part of the tune in radians
!     %param%stable   -- Set true or false.
!     %param%growth_rate -- unstable growth rate (= 0 if stable)
! 
!   bmad_status  -- BMAD Common block status structure
!     %ok            -- Logical: .True. if everything is OK, False otherwise.
!     %status        -- Integer: Calculation status.
!                        See MAT_SYMP_DECOUPLE for for more info
!-

#include "CESR_platform.inc"

subroutine twiss_at_start (ring)

  use bmad_struct
  use bmad_interface, except => twiss_at_start

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele

  real(rp) eta_vec(4), t0_4(4,4), mat6(6,6), error, map0(4)

  integer i, j, n, iu, n_lines

  logical :: debug = .false. 

  character(100), pointer :: lines(:)

! init one turn. T0 is the transverse part of the matrix

  bmad_status%ok = .false.             ! assume the worst

  call mat_make_unit (t0_4)       ! form unit matrix
  eta_vec = 0
  map0 = 0

! Propagate the transfer map around ring. 
! Since the RF is taken to be off we use a trick so we only have to multiply
! 4x4 matrices.

  if (debug) then
    iu = lunget()
    open (iu, file = 'twiss_at_start.dat')
  endif

  do n = 1, ring%n_ele_use
    ele => ring%ele_(n)
    eta_vec = matmul (ele%mat6(1:4,1:4), eta_vec)
    eta_vec = eta_vec + ele%mat6(1:4,6)
    map0 = matmul (ele%mat6(1:4,1:4), map0) + ele%vec0(1:4)
    t0_4 = matmul (ele%mat6(1:4,1:4), t0_4)
    if (debug) then
      write (iu, *) '!------------------------------------', n
      call type2_ele (ele, lines, n_lines, .false., 0, .false., 0)
      do i = 1, n_lines
        write (iu, *) lines(i)
      enddo
      deallocate (lines)
      call mat_symp_check (t0_4, error)
      write (iu, *) 'Symplectic Check:', error
      do i = 1, 4
        write (iu, '(4f15.10, 5x, 2f15.10)') (t0_4(i, j), j = 1, 4), &
                                                      eta_vec(i), map0(i)
      enddo
    endif
  enddo

  if (debug) close (iu)

! Put 1-turn matrix into ring%param%t1_no_RF

  call mat_make_unit (mat6)
  mat6(1:4,1:4) = t0_4

  call mat6_dispersion (eta_vec, mat6) ! dispersion to %mat6
  ring%param%t1_no_RF = mat6

! compute twiss parameters

  call twiss_from_mat6 (mat6, map0, ring%ele_(0), &
                                  ring%param%stable, ring%param%growth_rate)
  ring%x%tune = ring%ele_(0)%x%phi
  ring%y%tune = ring%ele_(0)%y%phi
  ring%ele_(0)%x%phi = 0
  ring%ele_(0)%y%phi = 0

end subroutine
