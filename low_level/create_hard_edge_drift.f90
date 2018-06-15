!+
! Subroutine create_hard_edge_drift (ele_in, which_end, drift_ele)
!
! Routine to create the drift element for the end drifts of an element.
! The end drifts are present, for example, when doing runge_kutta tracking
! through an rf_cavity with field_calc = bmad_standard. In this case, the
! field model is a pi-wave hard-edge resonator whose length may not match
! the length of the element and so particles must be drifted from the edge
! of the element to the edge of the field model.
!
! Input:
!   ele_in     -- ele_struct: Input element.
!   which_end  -- Integer: Which end is being created. upstream_end$ or downstream_end$.
!                   For an Lcavity one can have differences in reference energy.
!
! Output:
!   drift_ele  -- ele_struct: drift elment.
!-

subroutine create_hard_edge_drift (ele_in, which_end, drift_ele)

use equal_mod, dummy => create_hard_edge_drift
implicit none

type (ele_struct) ele_in, drift_ele
real(rp) E_tot, p0c
integer which_end

!

select case (which_end)
case (upstream_end$)
  e_tot = ele_in%value(e_tot_start$)
  p0c   = ele_in%value(p0c_start$)
  drift_ele%name                   = 'drift1_' // ele_in%name(1:33)
case (downstream_end$)
  e_tot = ele_in%value(e_tot$)
  p0c   = ele_in%value(p0c$)
  drift_ele%name                   = 'drift2_' // ele_in%name(1:33)
case default
  if (global_com%exit_on_error) call err_exit
end select

drift_ele%key                    = drift$

drift_ele%value                  = 0
drift_ele%value(p0c$)            = p0c
drift_ele%value(e_tot$)          = e_tot
drift_ele%value(p0c_start$)      = p0c
drift_ele%value(e_tot_start$)    = e_tot
drift_ele%value(l$)              = (ele_in%value(l$) - hard_edge_model_length(ele_in)) / 2 
drift_ele%value(ds_step$)        = drift_ele%value(l$)
drift_ele%value(num_steps$)      = 1
drift_ele%value(delta_ref_time$) = drift_ele%value(l$) * e_tot / (c_light * p0c)
drift_ele%orientation            = ele_in%orientation

end subroutine create_hard_edge_drift 
