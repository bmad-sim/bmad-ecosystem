!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_beam_emit_calc (plane, emit_type, ele, bunch_params) result (emit)
!
! Routine to calculate the emittance from beam parameters.
! The "projected" emittance is:
!   emit = sqrt(sigma_xx * sigma_pp - sigma_xp^2)
! Where the sigmas are calculated including coupling effects but assuming that sigma_pz = 0.
! The "apparent" emittance is:
!   emit = sqrt(sigma_xx - sigma_xp^2 / sigma_pp) / beta
!
! Input:
!   plane        -- Integer: x_plane$ or y_plane$.
!   emit_type    -- Integer: Either projected_emit$ or apparent_emit$
!   ele          -- ele_struct: Element.
!   bunch_params -- bunch_params_struct: Bunch sigma matrix
!
! Output:
!   emit -- Real(rp): emittance.
!-

function tao_beam_emit_calc (plane, emit_type, ele, bunch_params) result (emit)

use tao_struct

implicit none

type (ele_struct) ele
type (bunch_params_struct) bunch_params
real(rp) emit
integer plane, emit_type

!

if (plane == x_plane$) then
  if (emit_type == projected_emit$) then
    emit = bunch_params%x%emit
  elseif (emit_type == apparent_emit$) then
    if (bunch_params%sigma(6,6) == 0) then
      emit = bunch_params%sigma(1,1) / ele%a%beta
    else
      emit = (bunch_params%sigma(1,1) - bunch_params%sigma(1,6)**2 / bunch_params%sigma(6,6)) / ele%a%beta
    endif
  else
    call err_exit    
  endif

elseif (plane == y_plane$) then
  if (emit_type == projected_emit$) then
    emit = bunch_params%y%emit
  elseif (emit_type == apparent_emit$) then
    if (bunch_params%sigma(6,6) == 0) then
      emit = bunch_params%sigma(3,3) / ele%b%beta
    else
      emit = (bunch_params%sigma(3,3) - bunch_params%sigma(3,6)**2 / bunch_params%sigma(6,6)) / ele%b%beta
    endif
  else
    call err_exit    
  endif

else
  call err_exit
endif

end function tao_beam_emit_calc
