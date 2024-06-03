!+
! Subroutine twiss_at_start (lat, status, ix_branch, type_out)
!
! Subroutine to calculate, for a circular machine, the closed 1-turn 
! solution for the Twiss parameters at the start of the lat.
!
! Input:
!   lat         -- lat_struct: Lat
!   ix_branch   -- Integer, optional: Branch to use. Default is 0 (main branch).
!   type_out    -- logical, optional: If True (the default), print an error message
!                    If the 1-turn matrix is unstable.
!
! Output:
!   lat         -- Lat_struct: Lattice with twiss parameters computed.
!     %param%t1_no_RF --  Note: Only the linear part is computed.
!     %ele(0)%a      -- "a" mode Twiss parameters at the start of the lat.
!     %ele(0)%b      -- "b" mode Twiss parameters at the start of the lat.
!     %ele(0)%c_mat  -- Coupling matrix.
!     %a%tune         -- Fractional part of the tune in radians
!     %b%tune         -- Fractional part of the tune in radians
!     %a%stable       -- Set True or False.
!     %b%stable       -- Set True or False.
!     %param%stable   -- Set True or False.
!     %param%unstable_factor -- unstable growth rate (= 0 if stable)
!   status      -- Integer, optional: Calculation status:
!                       ok$, in_stop_band$, unstable$, or non_symplectic$
!-

subroutine twiss_at_start (lat, status, ix_branch, type_out)

use bookkeeper_mod, except_dummy => twiss_at_start

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (ele_struct), target :: drift_ele
type (branch_struct), pointer :: branch

real(rp) eta_vec(5), t0_5(5,5), mat6(6,6), map0(5), m56

integer, optional, intent(in) :: ix_branch
integer, optional, intent(out) ::status
integer i, j, n, iu, n_lines, stat

logical, optional :: type_out
logical :: debug = .false. 
logical saved_state

character(200), allocatable :: lines(:)

! Init one turn. T0 is the transverse part of the matrix

call mat_make_unit (t0_5)       ! form unit matrix
eta_vec = 0
map0 = 0
m56 = 0

branch => lat%branch(integer_option(0, ix_branch))


! Propagate the transfer map around branch. 
! Since the RF is taken to be off we use a trick so we only have to multiply
! 5x5 matrices.

if (debug) then
  iu = lunget()
  open (iu, file = 'twiss_at_start.dat')
endif

!

do n = 1, branch%n_ele_track
  ele => branch%ele(n)
  if (ele%key == rfcavity$) then
    drift_ele%map_ref_orb_out = ele%map_ref_orb_in
    call mat_make_unit(drift_ele%mat6)
    call track_a_drift(drift_ele%map_ref_orb_out, ele%value(l$), drift_ele%mat6, .true.)
    drift_ele%vec0 = 0
    ele => drift_ele
  endif

  eta_vec = matmul (ele%mat6(1:5,1:5), eta_vec) + ele%mat6(1:5,6)
  map0 = matmul (ele%mat6(1:5,1:5), map0) + ele%vec0(1:5)
  t0_5 = matmul (ele%mat6(1:5,1:5), t0_5)
  if (debug) then
    write (iu, *) '!------------------------------------', n
    call type_ele (ele, .false., 0, .false., 0, lines = lines, n_lines = n_lines)
    do i = 1, n_lines
      write (iu, '(a)') trim(lines(i))
    enddo
    write (iu, *) 'Symplectic Check:', mat_symp_error(t0_5(1:4,1:4))
    do i = 1, 5
      write (iu, '(5x, 5f14.8, 6x, 5f14.8)') t0_5(i,:), ele%mat6(i,1:5)
    enddo
    write (iu, '(a, 6es18.10)') 'Map_ref_orb:', ele%map_ref_orb_in%vec
    write (iu, '(a, 5es20.12)') 'Eta: ', eta_vec
    write (iu, '(a, 5es20.12)') 'Map0:', map0
  endif
enddo

if (debug) close (iu)

! Put 1-turn matrix into branch%param%t1_no_RF

call mat_make_unit (mat6)
mat6(1:5,1:5) = t0_5
mat6(1:5, 6) = eta_vec
branch%param%t1_no_RF = mat6

! Compute twiss parameters

call twiss_from_mat6 (mat6, branch%ele(1)%map_ref_orb_in%vec, branch%ele(0), branch%param%stable, &
                                   branch%param%unstable_factor, stat, logic_option(.true., type_out))
if (present(status)) status = stat

lat%a%tune = branch%ele(0)%a%phi
lat%b%tune = branch%ele(0)%b%phi
lat%a%stable = (stat == ok$)
lat%b%stable = (stat == ok$)

branch%ele(0)%a%phi = 0
branch%ele(0)%b%phi = 0

end subroutine
