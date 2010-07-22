!+
! Subroutine twiss_at_start (lat, ix_branch)
!
! Subroutine to calculate, for a circular machine, the closed 1-turn 
! solution for the Twiss parameters at the start of the lat.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat         -- lat_struct: Lat
!   ix_branch   -- Integer, option: Branch to use. Default is 0 (main branch).
!   bmad_status -- BMAD Common block status structure
!     %type_out  -- Logical: If .true. then will type a message for
!                       non ok$ STATUS
!
! Output:
!   lat
!     %param%t1_no_RF --  Note: Only the linear part is computed.
!     %ele(0)%a      -- "a" mode Twiss parameters at the start of the lat.
!     %ele(0)%b      -- "b" mode Twiss parameters at the start of the lat.
!     %ele(0)%c_mat  -- Coupling matrix.
!     %a%tune         -- Fractional part of the tune in radians
!     %b%tune         -- Fractional part of the tune in radians
!     %param%stable   -- Set true or false.
!     %param%unstable_factor -- unstable growth rate (= 0 if stable)
! 
!   bmad_status  -- BMAD Common block status structure
!     %ok            -- Logical: .True. if everything is OK, False otherwise.
!     %status        -- Integer: Calculation status.
!                        See MAT_SYMP_DECOUPLE for for more info
!-

subroutine twiss_at_start (lat, ix_branch)

use bmad_struct
use bookkeeper_mod, except_dummy => twiss_at_start

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch

real(rp) eta_vec(4), t0_4(4,4), mat6(6,6), map0(4)

integer, optional :: ix_branch
integer i, j, n, iu, n_lines

logical :: debug = .false. 
logical saved_state

character(200), pointer :: lines(:)

! init one turn. T0 is the transverse part of the matrix

bmad_status%ok = .false.             ! assume the worst

call mat_make_unit (t0_4)       ! form unit matrix
eta_vec = 0
map0 = 0

branch => lat%branch(integer_option(0, ix_branch))

! Propagate the transfer map around branch. 
! Since the RF is taken to be off we use a trick so we only have to multiply
! 4x4 matrices.

if (debug) then
  iu = lunget()
  open (iu, file = 'twiss_at_start.dat')
endif

call set_on_off (rfcavity$, lat, save_state$)
call set_on_off (rfcavity$, lat, off$, use_ref_orb = .true.)

do n = 1, branch%n_ele_track
  ele => branch%ele(n)
  eta_vec = matmul (ele%mat6(1:4,1:4), eta_vec)
  eta_vec = eta_vec + ele%mat6(1:4,6)
  map0 = matmul (ele%mat6(1:4,1:4), map0) + ele%vec0(1:4)
  t0_4 = matmul (ele%mat6(1:4,1:4), t0_4)
  if (debug) then
    write (iu, *) '!------------------------------------', n
    call type2_ele (ele, lines, n_lines, .false., 0, .false., 0)
    do i = 1, n_lines
      write (iu, '(a)') lines(i)
    enddo
    deallocate (lines)
    write (iu, *) 'Symplectic Check:', mat_symp_error(t0_4)

    do i = 1, 4
      write (iu, '(4f18.13)') (t0_4(i, j), j = 1, 4)
    enddo
    do i = 1, 4
      write (iu, '(es20.12, 2x, es20.12)') eta_vec(i), map0(i)
    enddo
  endif
enddo

call set_on_off (rfcavity$, lat, restore_state$, use_ref_orb = .true.)

if (debug) close (iu)

! Put 1-turn matrix into branch%param%t1_no_RF

call mat_make_unit (mat6)
mat6(1:4,1:4) = t0_4

call mat6_dispersion (eta_vec, mat6) ! dispersion to %mat6
branch%param%t1_no_RF = mat6

! compute twiss parameters

call twiss_from_mat6 (mat6, map0, branch%ele(0), branch%param%stable, branch%param%unstable_factor)

if (integer_option(0, ix_branch) == 0) then
  lat%a%tune = lat%ele(0)%a%phi
  lat%b%tune = lat%ele(0)%b%phi
endif

branch%ele(0)%a%phi = 0
branch%ele(0)%b%phi = 0

end subroutine
