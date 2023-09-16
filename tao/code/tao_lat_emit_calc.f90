!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_lat_emit_calc (plane, emit_type, ele, modes) result (emit)
!
! Routine to calculate the emittance.
! The "projected" emittance is:
!   emit = sqrt(sigma_xx * sigma_pp - sigma_xp^2)
! Where the sigmas are calculated including coupling effects but assuming that sigma_pz = 0.
! The "apparent" emittance is:
!   emit = sqrt(sigma_xx - sigma_xp^2 / sigma_pp) / beta
!
! Input:
!   plane     -- Integer: x_plane$ or y_plane$.
!   emit_type -- Integer: Either projected_emit$ or apparent_emit$
!   ele       -- ele_struct: Element holding the Twiss and coupling parameters.
!   modes     -- normal_modes_struct: Structure holding the emittances
!
! Output:
!   emit -- Real(rp): emittance.
!-

function tao_lat_emit_calc (plane, emit_type, ele, modes) result (emit)

use tao_struct

implicit none

type (ele_struct) ele
type (normal_modes_struct) modes

real(rp) s_mat(4,4), v_mat(4,4), emit
real(rp), save :: a_mat(4,4) = 0
integer plane, emit_type

!

a_mat(1,1) =  modes%a%emittance * ele%a%beta
a_mat(1,2) = -modes%a%emittance * ele%a%alpha
a_mat(2,2) =  modes%a%emittance * ele%a%gamma
a_mat(2,1) = a_mat(1,2)

a_mat(3,3) =  modes%b%emittance * ele%b%beta
a_mat(3,4) = -modes%b%emittance * ele%b%alpha
a_mat(4,4) =  modes%b%emittance * ele%b%gamma
a_mat(4,3) = a_mat(3,4)

call make_v_mats (ele, v_mat)
s_mat = matmul(matmul(v_mat, a_mat), transpose(v_mat))

if (plane == x_plane$) then
  if (emit_type == projected_emit$) then
    emit = sqrt(s_mat(1,1) * s_mat(2,2) - s_mat(1,2)**2)
  elseif (emit_type == apparent_emit$) then
    emit = s_mat(1,1) / ele%a%beta
  else
    call err_exit    
  endif

elseif (plane == y_plane$) then
  if (emit_type == projected_emit$) then
    emit = sqrt(s_mat(3,3) * s_mat(4,4) - s_mat(3,4)**2)
  elseif (emit_type == apparent_emit$) then
    emit = s_mat(3,3) / ele%b%beta
  else
    call err_exit    
  endif

else
  call err_exit
endif

end function tao_lat_emit_calc
