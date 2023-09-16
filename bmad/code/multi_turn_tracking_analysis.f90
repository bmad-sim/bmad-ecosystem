!+
! Subroutine multi_turn_tracking_analysis (track, i_dim, track0, ele,
!                                                stable, growth_rate, chi, err_flag)
!
! Subroutine to analyze multi-turn tracking data. 
! Calculated is the symplectified 1-turn matrix and the Twiss parameters. 
! If you just want the 1-turn matrix see the routine:
!        multi_turn_tracking_to_mat.
!
! Input:
!   track(:) -- Cooord_struct: multi-turn tracking data to analyze.
!                track(i) is the particle position at a given point
!                in the lat on the i^th turn.
!   i_dim    -- Integer: number of dimensions used in the tracking: 2, or 4.
!
! Output:
!   track0      -- Coord_struct: Closed orbit.
!   ele         -- ele_struct: structure holding the 1-turn matrix and Twiss parameters.
!     %mat6        -- Symplectified 1-turn matrix. If you want the
!                      true non-symplectified 1-turn matrix use the routine
!                      multi_turn_tracking_to_mat.
!     %a%beta, etc -- a-mode beta,  etc.  
!     %a%phi       -- a-mode fractional tune in radians.
!     %a%sigma     -- a-mode amplitude = sqrt(ele%a%beta * ele%a%sigma)
!     %c_mat       -- c coupling matrix (only with i_dim = 4)
!   stable      -- Logical: Is motion stable?
!   growth_rate -- Real(rp): Unstable growth rate (= 0 if stable).
!   chi         -- Real(rp): How symplectic the computed 1-turn matrix is.
!                    See mat_symp_check for more details.
!   err_flag    -- Logical: Set true if there is an error. False otherwise.
!-

subroutine multi_turn_tracking_analysis (track, i_dim, track0, ele, &
                                                  stable, growth_rate, chi, err_flag)

use bmad_interface, except_dummy => multi_turn_tracking_analysis

implicit none

type (coord_struct), intent(in) :: track(:)
type (coord_struct), intent(out) :: track0
type (ele_struct) :: ele

real(rp), intent(out) :: growth_rate, chi
real(rp) a_vec(4), v_mat(4,4), v_inv_mat(4,4), map0(6)

integer, intent(in) :: i_dim
logical, intent(out) :: stable, err_flag

integer stat, i

! get 1-turn matrix and symplectify

err_flag = .true.

call mat_make_unit (ele%mat6)
call multi_turn_tracking_to_mat (track, i_dim, ele%mat6(1:i_dim,1:i_dim), map0, track0, chi)
ele%vec0 = track0%vec
call mat_symplectify (ele%mat6(1:i_dim,1:i_dim), ele%mat6(1:i_dim,1:i_dim))

! get twiss parameters

if (i_dim == 2) then

  call twiss_from_mat2 (ele%mat6(1:2,1:2), ele%a, stat, .false.)

  if (stat /= ok$) then
    stable = .false.
    growth_rate = sqrt(abs(ele%mat6(1,1) + ele%mat6(2,2)) - 2)
  else
    stable = .true.
    growth_rate = 0
    ele%a%sigma = (ele%a%gamma * sum(track(:)%vec(1)**2) + &
                   ele%a%alpha * sum(track(:)%vec(1) * track(:)%vec(2)) + &
                   ele%a%beta *  sum(track(:)%vec(2)**2)) / size(track)
  endif

elseif (i_dim == 4) then

  call twiss_from_mat6 (ele%mat6, track0%vec, ele, stable, growth_rate, stat, .true.)
  if (stat /= ok$) return

  call make_v_mats (ele, v_mat, v_inv_mat)

  ele%a%sigma = 0
  ele%b%sigma = 0
  do i = 1, size(track)
    a_vec = matmul(v_inv_mat, track(i)%vec(1:4))
    ele%a%sigma = ele%a%sigma + &
                   ele%a%gamma * a_vec(1)**2 + &
                   ele%a%alpha * a_vec(1) * a_vec(2) + &
                   ele%a%beta *  a_vec(2)**2 
    ele%b%sigma = ele%b%sigma + &
                   ele%b%gamma * a_vec(3)**2 + &
                   ele%b%alpha * a_vec(3) * a_vec(4) + &
                   ele%b%beta *  a_vec(4)**2 
  enddo
  ele%a%sigma = ele%a%sigma / size(track)
  ele%b%sigma = ele%b%sigma / size(track)

endif

err_flag = .false.

end subroutine
