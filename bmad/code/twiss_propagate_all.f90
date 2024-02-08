!+
! Subroutine twiss_propagate_all (lat, ix_branch, err_flag, ie_start, ie_end, zero_uncalculated)
!
! Subroutine to propagate the twiss, coupling, and dispersion parameters from 
! the start to the end of a lattice branch.
!
!
! Input:
!   lat                -- lat_struct: lattice.
!     %branch(ix_branch)%ele(0) -- Branch beginning element with the starting parameters.
!   ix_branch          -- integer, optional: Branch index. Default is 0 (main lattice).
!   ie_start           -- integer, optional: Starting element index. Default is 0.
!                           Note: The first element at which the Twiss parameters are calculated is ie_start+1.
!   ie_end             -- integer, optional: Ending element index, Default is branch%n_ele_track.
!   zero_uncalculated  -- logical, optional: Set to zero Twiss parameters not calculated in 
!                           range [ie_start, ie_end]? Default is True.
!
! Output:
!   lat          -- lat_struct: Lattice with parameters computed for the branch.
!   err_flag     -- logical, optional: Set True if there is an error. False otherwise.
!-

subroutine twiss_propagate_all (lat, ix_branch, err_flag, ie_start, ie_end, zero_uncalculated)

use bmad_interface, except_dummy => twiss_propagate_all

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: lord, slave, ele, ele0

real(rp) v_inv_mat(4,4), eta_vec(4)

integer n, n_track, i_start, i_end
integer, optional :: ix_branch, ie_start, ie_end

logical, optional :: err_flag, zero_uncalculated
logical err

character(*), parameter :: r_name = 'twiss_propagate_all'

! Twiss parameters for photons do not make much sense so do not bother to calculate anything.

branch => lat%branch(integer_option(0, ix_branch))
if (branch%param%particle == photon$) return

! Make sure gamma for ele(0) is correct.

i_start = integer_option(0, ie_start)
n_track = branch%n_ele_track
i_end = integer_option(n_track, ie_end)

ele => branch%ele(i_start)
if (ele%a%beta /= 0) ele%a%gamma = (1 + ele%a%alpha**2) / ele%a%beta
if (ele%b%beta /= 0) ele%b%gamma = (1 + ele%b%alpha**2) / ele%b%beta

call make_v_mats (ele, v_inv_mat = v_inv_mat)
eta_vec = [ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap]
eta_vec = matmul (v_inv_mat, eta_vec)
ele%a%eta  = eta_vec(1)
ele%a%etap = eta_vec(2)
ele%b%eta  = eta_vec(3)
ele%b%etap = eta_vec(4)

ele%x%deta_ds = ele%x%etap
ele%y%deta_ds = ele%y%etap
ele%a%deta_ds = eta_vec(2)
ele%b%deta_ds = eta_vec(4)

! Propagate twiss

if (present(err_flag)) err_flag = .true.

do n = i_start+1, i_end
  ele => branch%ele(n)
  call twiss_propagate1 (branch%ele(n-1), ele, err)
  if (err) return
enddo

if (logic_option(.true., zero_uncalculated)) then
  do n = 0, i_start-1
    branch%ele(n)%a = twiss_struct()
    branch%ele(n)%b = twiss_struct()
    branch%ele(n)%x = xy_disp_struct()
    branch%ele(n)%y = xy_disp_struct()
  enddo

  do n = i_end+1, n_track
    branch%ele(n)%a = twiss_struct()
    branch%ele(n)%b = twiss_struct()
    branch%ele(n)%x = xy_disp_struct()
    branch%ele(n)%y = xy_disp_struct()
  enddo
endif

! Make sure final mode is same as initial mode

if (branch%param%geometry == closed$ .and. i_start == 0 .and. i_end == n_track) then
  if (branch%ele(0)%mode_flip .neqv. branch%ele(n_track)%mode_flip) then
    call do_mode_flip (branch%ele(n_track), err)
    if (err) call out_io (s_error$, r_name,  'CANNOT MAKE FINAL FLIP STATE EQUAL TO THE INITIAL')
  endif
endif

! Super_lord elements get the twiss parameters at the exit end

do n = lat%n_ele_track + 1, lat%n_ele_max
  lord => lat%ele(n)
  if (lord%n_slave == 0) cycle
  select case (lord%lord_status) 
  case (super_lord$, overlay_lord$, group_lord$)
    slave => pointer_to_slave(lord, lord%n_slave)
    if (slave%ix_branch /= branch%ix_branch) cycle
    call transfer_twiss (slave, lord)
  end select
enddo

if (present(err_flag)) err_flag = .false.

end subroutine
