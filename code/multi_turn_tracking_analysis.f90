!+
! Subroutine multi_turn_tracking_analysis (track, i_dim, track0, ele,
!                                                stable, growth_rate, chi)
!
! Subroutine to analyze multi-turn tracking data. 
! Calculated is the symplectified 1-turn matrix and the Twiss parameters. 
! If you just want the 1-turn matrix see the routine:
!        multi_turn_tracking_to_mat.
!
! Modules needed:
!   use bmad
!
! Input:
!   track(:) -- Cooord_struct: multi-turn tracking data to analyze.
!                track(i) is the particle position at a given point
!                in the ring on the i^th turn.
!   i_dim    -- Integer: number of dimensions used in the tracking: 2, or 4.
!
! Output:
!   track0 -- Coord_struct: Closed orbit.
!   ele    -- ele_struct: structure holding the 1-turn matrix and 
!                                                      Twiss parameters.
!     %mat6        -- Symplectified 1-turn matrix. If you want the
!                      true non-symplectified 1-turn matrix use the routine
!                      multi_turn_tracking_to_mat.
!     %x%beta, etc -- a-mode beta,  etc.  
!     %x%phi       -- a-mode fractional tune in radians.
!     %x%sigma     -- a-mode amplitude = sqrt(ele%x%beta * ele%x%sigma)
!     %c_mat       -- c coupling matrix (only with i_dim = 4)
!   stable -- Logical: Is motion stable?
!   growth_rate -- Real(rdef): Unstable growth rate (= 0 if stable).
!   chi    -- Real(rdef): How symplectic the computed 1-turn matrix is.
!              See mat_symp_check for more details.
!-

#include "CESR_platform.inc"
                       
subroutine multi_turn_tracking_analysis (track, i_dim, track0, ele, &
                                                 stable, growth_rate, chi)

  use bmad_struct
  use bmad_interface

  implicit none

  type (coord_struct), intent(in) :: track(:)
  type (coord_struct), intent(out) :: track0
  type (ele_struct), intent(out) :: ele

  real(rdef), intent(out) :: growth_rate, chi
  real(rdef) a_vec(4), v_mat(4,4), v_inv_mat(4,4)

  integer, intent(in) :: i_dim
  logical, intent(out) :: stable

  real(rdef) det
  integer stat, i

! get 1-turn matrix and symplectify

  call mat_make_unit (ele%mat6)
  call multi_turn_tracking_to_mat (track, i_dim, ele%mat6(1:i_dim,1:i_dim), &
                                                                 track0, chi)
  call mat_symplectify (ele%mat6(1:i_dim,1:i_dim), ele%mat6(1:i_dim,1:i_dim))

! get twiss parameters

  if (i_dim == 2) then

    call twiss_from_mat2 (ele%mat6(1:2,1:2), det, ele%x, stat, 1e-4_rdef, .false.)

    if (stat == unstable$) then
      stable = .false.
      growth_rate = sqrt(abs(ele%mat6(1,1) + ele%mat6(2,2)) - 2)
    else
      stable = .true.
      growth_rate = 0
      ele%x%sigma = (ele%x%gamma * sum(track(:)%vec(1)**2) + &
                     ele%x%alpha * sum(track(:)%vec(1) * track(:)%vec(2)) + &
                     ele%x%beta *  sum(track(:)%vec(2)**2)) / size(track)
    endif

  elseif (i_dim == 4) then

    call twiss_from_mat6 (ele%mat6, ele, stable, growth_rate)
    call make_v_mats (ele, v_mat, v_inv_mat)

    ele%x%sigma = 0
    ele%y%sigma = 0
    do i = 1, size(track)
      a_vec = matmul(v_inv_mat, track(i)%vec(1:4))
      ele%x%sigma = ele%x%sigma + &
                     ele%x%gamma * a_vec(1)**2 + &
                     ele%x%alpha * a_vec(1) * a_vec(2) + &
                     ele%x%beta *  a_vec(2)**2 
      ele%y%sigma = ele%y%sigma + &
                     ele%y%gamma * a_vec(3)**2 + &
                     ele%y%alpha * a_vec(3) * a_vec(4) + &
                     ele%y%beta *  a_vec(4)**2 
    enddo
    ele%x%sigma = ele%x%sigma / size(track)
    ele%y%sigma = ele%y%sigma / size(track)

  endif

end subroutine
