!+
! Subroutine track_a_patch (ele, orbit, drift_to_exit, s_ent, ds_ref,  track_spin, mat6, make_matrix)
! 
! Bmad_standard routine to track through a patch element.
!
! The steps for tracking are:
!   1) Transform from entrance to exit coordinates.
!   2) Drift particle from the entrance to the exit coordinants.
!
! Input:
!   ele           -- ele_struct: patch element.
!   orbit         -- coord_struct: Starting phase space coords
!   drift_to_exit -- logical, optional: If False then do not drift the particle from
!                      beginning to end face. Also do not correct for a reference energy shift.
!                      Default is True. 
!   track_spin    -- logical, optional: If True rotate the spin vector appropriately. 
!                       If ele%spin_tracking_method = symp_lie_ptc -> default = True. Else default = False.
!   make_matrix   -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: Coords after applying a patch transformation.
!   s_ent      -- real(rp), optional: Longitudinal coordinate of the initial particle 
!                   position in the frame of reference of the face where the particle exits.
!                   For a patch with positive z_offset and all other attributes zero, s_ent = -z_offset.
!   ds_ref     -- real(rp), optional: Distance reference particle travels from entrance to exit.
!   mat6(6,6)  -- Real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_patch (ele, orbit, drift_to_exit, s_ent, ds_ref, track_spin, mat6, make_matrix)

use bmad_interface, dummy => track_a_patch

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit, orb_in

real(rp), pointer :: v(:)
real(rp) p_vec(3), r_vec(3), rel_p, ww(3,3), beta0, ds0, s_vec(3), r_s
real(rp) pz_out, mc2, pz_in, beta_ref, dps_dpx, dps_dpy, dps_dpz, dp_ratio
real(rp), optional :: mat6(6,6), s_ent, ds_ref

integer dir

logical, optional :: drift_to_exit, track_spin, make_matrix

character(*), parameter :: r_name = 'track_a_patch'

! Transform to exit face coords.

v => ele%value

if (logic_option(.true., make_matrix)) orb_in = orbit
rel_p = 1 + orbit%vec(6)
p_vec = [orbit%vec(2), orbit%vec(4), sqrt(rel_p**2 - orbit%vec(2)**2 - orbit%vec(4)**2)]

! For other types of elements, the ele%orientation is the same as the upstream and downstream 
! elements orientation. For a patch this is not necessarily true which is why a patch element
! needs to store the upstream and downstream orientations.

! orbit%direction * orbit%time_dir * ele%orientation = -1 means we are going from the exit end to the entrance end.
! In this case, must reverse the order of application of the offsets and pitches.

if (orbit%direction * orbit%time_dir * ele%orientation == 1) then
  if (ele%orientation == 1) then ! Entering from upstream
    p_vec(3) = p_vec(3) * ele%value(upstream_coord_dir$) * orbit%direction  
  else    ! Entering from downstream
    p_vec(3) = p_vec(3) * ele%value(downstream_coord_dir$) * orbit%direction  
  endif
  r_vec = [orbit%vec(1) - v(x_offset$), orbit%vec(3) - v(y_offset$), -v(z_offset$)]
  if (v(x_pitch$) /= 0 .or. v(y_pitch$) /= 0 .or. v(tilt$) /= 0) then
    call floor_angles_to_w_mat (v(x_pitch$), v(y_pitch$), v(tilt$), w_mat_inv = ww)
    p_vec = matmul(ww, p_vec)
    r_vec = matmul(ww, r_vec)
    orbit%vec(2) = p_vec(1)
    orbit%vec(4) = p_vec(2)
    if (logic_option((ele%spin_tracking_method /= symp_lie_ptc$), track_spin)) then
      orbit%spin = matmul(ww, orbit%spin)
    endif
  else
    call mat_make_unit (ww)
  endif
  ds0 = ww(3,1) * v(x_offset$) + ww(3,2) * v(y_offset$) + ww(3,3) * v(z_offset$)
  orbit%vec(5) = orbit%vec(5) + orbit%beta * c_light * v(t_offset$)

else
  if (ele%orientation == 1) then ! Entering from downstream
    p_vec(3) = p_vec(3) * ele%value(downstream_coord_dir$) * orbit%direction  
  else    ! Entering from upstream
    p_vec(3) = p_vec(3) * ele%value(upstream_coord_dir$) * orbit%direction  
  endif
  r_vec = [orbit%vec(1), orbit%vec(3), 0.0_rp]
  if (v(x_pitch$) /= 0 .or. v(y_pitch$) /= 0 .or. v(tilt$) /= 0) then
    call floor_angles_to_w_mat (v(x_pitch$), v(y_pitch$), v(tilt$), w_mat = ww)
    p_vec = matmul(ww, p_vec)
    r_vec = matmul(ww, r_vec)
    orbit%vec(2) = p_vec(1)
    orbit%vec(4) = p_vec(2)
    if (logic_option((ele%spin_tracking_method /= symp_lie_ptc$), track_spin)) then
      orbit%spin = matmul(ww, orbit%spin)
    endif
  else
    call mat_make_unit (ww)
  endif
  r_vec = r_vec + [v(x_offset$), v(y_offset$), v(z_offset$)]
  ds0 = ww(1,3) * v(x_offset$) + ww(2,3) * v(y_offset$) + ww(3,3) * v(z_offset$)
  orbit%vec(5) = orbit%vec(5) - orbit%beta * c_light * v(t_offset$)
endif

!

orbit%vec(1) = r_vec(1)
orbit%vec(3) = r_vec(2)

if (present(s_ent))  s_ent = r_vec(3)
if (present(ds_ref)) ds_ref = ds0

! Drift to exit face.
! Notice that the drift distance is -r_vec(3). 

if (logic_option(.true., drift_to_exit)) then
  ! Set track edge so that energy correction does not ignore an energy shift.
  if (orbit%direction * orbit%time_dir == 1) then
    call reference_energy_correction (ele, orbit, first_track_edge$)
  else
    call reference_energy_correction (ele, orbit, second_track_edge$)
  endif
  beta0 = v(p0c$) / v(e_tot$)
  orbit%vec(1) = orbit%vec(1) - r_vec(3) * p_vec(1) / p_vec(3)
  orbit%vec(3) = orbit%vec(3) - r_vec(3) * p_vec(2) / p_vec(3)
  orbit%vec(5) = orbit%vec(5) + r_vec(3) * rel_p / p_vec(3) + orbit%time_dir * v(l$) * orbit%beta / beta0
  orbit%t = orbit%t - r_vec(3) * rel_p / (p_vec(3) * orbit%beta * c_light)
  orbit%s = orbit%s + v(l$)
endif

! Matrix

if (logic_option(.false., make_matrix)) then

  if (ele%field_calc == custom$) then
    call out_io (s_fatal$, r_name, 'MAT6_CALC_METHOD=BMAD_STANDARD CANNOT HANDLE FIELD_CALC=CUSTOM', &
                                   'FOR PATCH ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  pz_out = p_vec(3)
  r_s = r_vec(3)
  mc2 = mass_of(orb_in%species)

  ! Coordinate transform before drift

  rel_p = 1 + orb_in%vec(6)
  if (ele%orientation*orb_in%direction == 1) then
    dir = ele%value(upstream_coord_dir$) * orb_in%direction
  else
    dir = ele%value(downstream_coord_dir$) * orb_in%direction
  endif
  pz_in = sqrt(rel_p**2 - orb_in%vec(2)**2 - orb_in%vec(4)**2) * dir
  beta_ref = v(p0c$) / v(e_tot$)

  dps_dpx = ww(3,1) - ww(3,3) * orb_in%vec(2) / pz_in
  dps_dpy = ww(3,2) - ww(3,3) * orb_in%vec(4) / pz_in
  dps_dpz = ww(3,3) * rel_p / pz_in

  mat6(1,1) = ww(1,1) - ww(3,1) * orbit%vec(2) / pz_out
  mat6(1,2) = -r_s * ((ww(1,1) - ww(1,3) * orb_in%vec(2) / pz_in) / pz_out - orbit%vec(2) * dps_dpx / pz_out**2) 
  mat6(1,3) = ww(1,2) - ww(3,2) * orbit%vec(2) / pz_out
  mat6(1,4) = -r_s * ((ww(1,2) - ww(1,3) * orb_in%vec(4) / pz_in) / pz_out - orbit%vec(2) * dps_dpy / pz_out**2)
  mat6(1,6) = -r_s * (ww(1,3) * rel_p / (pz_out * pz_in) - orbit%vec(2) * dps_dpz / pz_out**2)

  mat6(2,2) = ww(1,1) - ww(1,3) * orb_in%vec(2) / pz_in
  mat6(2,4) = ww(1,2) - ww(1,3) * orb_in%vec(4) / pz_in
  mat6(2,6) = ww(1,3) * rel_p / pz_in

  mat6(3,1) = ww(2,1) - ww(3,1) * orbit%vec(4) / pz_out
  mat6(3,2) = -r_s * ((ww(2,1) - ww(2,3) * orb_in%vec(2) / pz_in) / pz_out - orbit%vec(4) * dps_dpx / pz_out**2)
  mat6(3,3) = ww(2,2) - ww(3,2) * orbit%vec(4) / pz_out
  mat6(3,4) = -r_s * ((ww(2,2) - ww(2,3) * orb_in%vec(4) / pz_in) / pz_out - orbit%vec(4) * dps_dpy / pz_out**2) 
  mat6(3,6) = -r_s * (ww(2,3) * rel_p / (pz_out * pz_in) - orbit%vec(4) * dps_dpz / pz_out**2)

  mat6(4,2) = ww(2,1) - ww(2,3) * orb_in%vec(2) / pz_in
  mat6(4,4) = ww(2,2) - ww(2,3) * orb_in%vec(4) / pz_in
  mat6(4,6) = ww(2,3) * rel_p / pz_in

  mat6(5,1) = ww(3,1) * rel_p / pz_out
  mat6(5,2) = -r_s * rel_p * dps_dpx / pz_out**2
  mat6(5,3) = ww(3,2) * rel_p / pz_out
  mat6(5,4) = -r_s * rel_p * dps_dpy / pz_out**2
  mat6(5,6) = v(t_offset$) * c_light * mc2**2 * orb_in%beta**3 / (v(p0c_start$)**2 * rel_p**3) + &
              r_s / pz_out - r_s * rel_p * dps_dpz / pz_out**2 + &
              v(l$) * mc2**2 * orbit%beta**3 / (rel_p**3 * v(p0c_start$)**2 * beta_ref)

  ! Energy offset

  if (v(p0c_start$) /= v(p0c$)) then
    dp_ratio = v(p0c_start$) / v(p0c$)
    mat6(2,:) = mat6(2,:) * dp_ratio
    mat6(4,:) = mat6(4,:) * dp_ratio
    mat6(6,:) = mat6(6,:) * dp_ratio
  endif

endif

end subroutine track_a_patch
