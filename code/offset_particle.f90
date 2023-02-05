!+
! Subroutine offset_particle (ele, set, orbit, set_tilt, set_hvkicks, drift_to_edge, &
!                                                           s_pos, s_out, set_spin, mat6, make_matrix, spin_qrot)
!
! Routine to transform a particles's coordinates between laboratory and element body coordinates.
! Additionally, this routine can:
!   a) Apply the half kicks due to multipole and kick attributes.
!   b) drift to the actual misaligned edge of the element.
!
! set = set$:
!    Transforms from lab to element coords. 
!    If s_pos is not present:
!      If coord%direction*coord%time_dir = +1 -> Assume the particle is at the upstream (-S) end.
!      If coord%direction*coord%time_dir = -1 -> Assume the particle is at the downstream (+S) end.
!
! set = unset$:
!    Transforms from element body to lab coords.
!    If s_pos is not present:
!      If coord%direction*coord%time_dir = +1 -> Assume the particle is at the downstream (+S) end.
!      If coord%direction*coord%time_dir = -1 -> Assume the particle is at the upstream (-S) end.
!
! Note: If ele%orientation = -1 then the upstream end is the exit end of the element and 
!   the downstream end is the entrance end of the element.
!
! Note: There are no element coordinates associated with a patch element so this routine will do nothing in this case.
! Note: x/y_offset_mult are included in sad_mult offsets.
!
! Input:
!   ele               -- Ele_struct: Element
!   set               -- Logical: 
!                       T (= set$)   -> Translate from lab coords to the local element coords.
!                       F (= unset$) -> Translate back from element to lab coords.
!   orbit             -- Coord_struct: Coordinates of the particle.
!   set_tilt          -- Logical, optional: Default is True.
!                       T -> Rotate using ele%value(tilt$) and ele%value(roll$) for sbends.
!                       F -> Do not rotate
!   set_hvkicks       -- Logical, optional: Default is True.
!                       T -> Apply 1/2 any hkick or vkick.
!   drift_to_edge     -- Logical, optional: Default is True if s_pos is absent and False if s_pos is present.
!                       T -> Particle will be propagated from where the particle is to:
!                            With set = set$: The upsteam body edge.
!                            With set = unset$: The nominal (lab coords) downstream edge of the element.
!                       F -> Do not propagate. 
!   s_pos             -- Real(rp), optional: Longitudinal particle position:
!                         If set = set$: Relative to upstream end (in lab coords).
!                         If set = unset$: Relative to entrance end (in body coords).
!   set_spin          -- Logical, optional: Default if False.
!                         Rotate spin coordinates? Also bmad_com%spin_tracking_on must be T to rotate.
!   mat6(6,6)         -- Real(rp), optional: Transfer matrix before off setting.
!   make_matrix       -- logical, optional: Propagate the transfer matrix? Default is false.
!                                               
! Output:
!     orbit           -- coord_struct: Coordinates of particle.
!                         If set = set$: In body coords.
!                         If set = unset$: In lab coords.
!     s_out           -- real(rp), optional: Longitudinal particle position. 
!                         If set = set$: Relative to entrance end (in body coords).
!                         If set = unset$: Relative to upstream end (in lab coords).
!     mat6(6,6)       -- real(rp), optional: Transfer matrix transfer matrix after offsets applied.
!     spin_qrot(0:3)  -- real(rp), optional: Spin rotation quaternion
!-

subroutine offset_particle (ele, set, orbit, set_tilt, set_hvkicks, drift_to_edge, &
                                                            s_pos, s_out, set_spin, mat6, make_matrix, spin_qrot)

use bmad_interface, except_dummy => offset_particle

implicit none

type (ele_struct) :: ele
type (coord_struct), intent(inout) :: orbit
type (em_field_struct) field
type (floor_position_struct) position

real(rp), optional :: s_pos, s_out, mat6(6,6), spin_qrot(0:3)
real(rp) rel_p, knl(0:n_pole_maxx), tilt(0:n_pole_maxx), dx, f, B_factor, ds_center
real(rp) angle, xp, yp, x_off, y_off, z_off, off(3), m_trans(3,3), pz, s_target, s_lab, s_body
real(rp) beta_ref, charge_dir, dz, rel_tracking_charge, rtc, Ex, Ey, kx, ky, length, s_pos0
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), ws(3,3), L_mis(3), p_vec0(3), p_vec(3), ref_tilt
real(rp) L_half, ds

integer particle, sign_z_vel
integer n, ix_pole_max

logical, intent(in) :: set
logical, optional, intent(in) :: set_tilt, set_spin
logical, optional, intent(in) :: set_hvkicks, drift_to_edge
logical, optional :: make_matrix
logical set_hv, set_t, set_hv1, set_hv2, do_drift, set_spn, is_misaligned

!

length = ele%value(l$)
L_half = 0.5_rp * length
sign_z_vel = ele%orientation * orbit%direction

! set:   s_pos is in lab coords
! unset: s_pos is in body coords.

if (present(s_pos)) then
  s_pos0 = s_pos
elseif ((set .and. orbit%direction * orbit%time_dir == 1) .or. (.not. set .and. sign_z_vel* orbit%time_dir == -1)) then
  s_pos0 = 0
else
  s_pos0 = length
endif

! ds_center is distance to the center of the element from the particle position

if (set) then
  ds_center = L_half - s_pos0                       ! Lab coords: Nominal center - S_pos
  s_target = (1 - sign_z_vel * orbit%time_dir) * L_half              ! Position to drift to (body coords).
else
  ds_center = L_half - s_pos0                       ! Body coords: Body center - s_pos
  s_target = (1 + orbit%direction * orbit%time_dir) * L_half         ! Position to drift to (lab coords).
endif

!

if (ele%key == patch$) then
  if (present(s_out)) s_out = s_pos0
  return
endif

if (present(spin_qrot)) spin_qrot = [1, 0, 0, 0]

!---------------------------------------------------------------         

rel_p = (1 + orbit%vec(6))

set_hv     = logic_option (.true., set_hvkicks) .and. ele%is_on .and. &
                   (has_kick_attributes(ele%key) .or. has_hkick_attributes(ele%key))
set_t      = logic_option (.true., set_tilt) .and. has_orientation_attributes(ele)
do_drift   = logic_option (.not. present(s_pos), drift_to_edge) .and. has_orientation_attributes(ele)
set_spn    = (logic_option (.false., set_spin) .and. bmad_com%spin_tracking_on) .or. present(spin_qrot)

rel_tracking_charge = rel_tracking_charge_to_mass (orbit, ele%ref_species)
charge_dir = rel_tracking_charge * sign_z_vel * orbit%time_dir 

if (set_hv) then
  select case (ele%key)
  case (elseparator$, kicker$, hkicker$, vkicker$)
    set_hv1 = .false.
    set_hv2 = .true.
  case default
    set_hv1 = .true.
    set_hv2 = .false.
  end select
else
  set_hv1 = .false.
  set_hv2 = .false.
endif

B_factor = ele%value(p0c$) / (charge_of(ele%ref_species) * c_light)
is_misaligned = .false.

!----------------------------------------------------------------
! Set...

if (set) then

  ! Set: Offset and pitch

  if (has_orientation_attributes(ele)) then
    x_off = ele%value(x_offset_tot$)
    y_off = ele%value(y_offset_tot$)
    z_off = ele%value(z_offset_tot$)
    xp    = ele%value(x_pitch_tot$)
    yp    = ele%value(y_pitch_tot$)
    ref_tilt = ele%value(ref_tilt_tot$)

    if (ele%key == sad_mult$) then
      x_off = x_off + ele%value(x_offset_mult$)
      y_off = y_off + ele%value(y_offset_mult$)
    endif

    is_misaligned = (x_off /= 0 .or. y_off /= 0 .or. z_off /= 0 .or. xp /= 0 .or. yp /= 0 .or. &
                                (ele%key == sbend$ .and. (ref_tilt /= 0 .or. ele%value(roll$) /= 0)))
  endif

  if (is_misaligned) then
    position%r = [orbit%vec(1), orbit%vec(3), 0.0_rp]
    call mat_make_unit (position%w)

    if (ele%key == sbend$ .and. (ele%value(g$) /= 0 .or. ref_tilt /= 0 .or. ele%value(roll$) /= 0)) then
      position = bend_shift(position, ele%value(g$), ele%orientation*ds_center, ref_tilt = ref_tilt)

      call ele_misalignment_L_S_calc(ele, L_mis, ws)
      ws = transpose(ws)
      position%r = matmul(ws, position%r - L_mis)
      position%w = matmul(ws, position%w)

      if (ref_tilt /= 0) then
        position = bend_shift(position, ele%value(g$), -L_half, ref_tilt = ref_tilt)
        ws = w_mat_for_tilt(-ref_tilt)
        position%r = matmul(ws, position%r)
        position%w = matmul(ws, position%w)
        position = bend_shift(position, ele%value(g$), L_half)
      endif        

      if (do_drift) then
        position = bend_shift(position, ele%value(g$), s_target-L_half)
        position%r(3) = position%r(3) + s_target
      else
        if (ele%value(g$) == 0) then
          ds = -ele%orientation*ds_center
        else
          ds = ele%value(rho$) * atan2(position%r(3)*ele%value(g$), 1.0_rp + position%r(1)*ele%value(g$))
        endif
        position = bend_shift(position, ele%value(g$), ds)
        position%r(3) = position%r(3) + L_half + ds
      endif

    ! Else not a bend or a bend with zero bending angle

    else
      position%r  = position%r - [x_off, y_off, z_off+ele%orientation*ds_center]
      call floor_angles_to_w_mat (xp, yp, 0.0_rp, w_mat_inv = ws)
      position%r = matmul(ws, position%r)
      position%w = matmul(ws, position%w)
      position%r(3) = position%r(3) + L_half
    endif

    s_body = position%r(3)

    pz = rel_p**2 - orbit%vec(2)**2 - orbit%vec(4)**2
    if (pz <= 0) then
      orbit%state = lost_pz_aperture$
    else
      p_vec0 = [orbit%vec(2), orbit%vec(4), sign_z_vel * sqrt(pz)]
      p_vec = matmul(position%w, p_vec0)
      orbit%vec(2:4:2) = p_vec(1:2)
      orbit%vec(1:3:2) = position%r(1:2)
      if (logic_option(.false., make_matrix)) call apply_offsets_to_matrix (p_vec0, p_vec, position%w, mat6)
    endif


    if (set_spn) orbit%spin = matmul(position%w, orbit%spin)
    if (present(spin_qrot)) spin_qrot = quat_mul(w_mat_to_quat(position%w), spin_qrot)
  endif    ! has nonzero offset or pitch

  if (.not. is_misaligned) then
    select case (ele%orientation)
    case (1);  s_body = s_pos0
    case (-1); s_body = length - s_pos0
    end select
  endif

  ! When drifting the reference particle does not move! That is, the reference time is the time the 
  ! referece particle reaches the *nominal* edge of the element and this is independent of any misalignments.

  if (do_drift .and. orbit%state == alive$) then
    call track_a_drift (orbit, s_target - s_body, mat6, make_matrix, ele%orientation, include_ref_motion = .false.)
  endif

  !

  if (present(s_out)) then
    if (do_drift) then
      s_out = s_target
    else
      s_out = s_body
    endif
  endif

  ! Set: HV kicks for quads, etc. but not hkicker, vkicker, elsep and kicker elements.
  ! HV kicks must come after z_offset but before any tilts are applied.
  ! Note: Change in %vel is NOT dependent upon energy since we are using
  ! canonical momentum.
  ! Note: Since this is applied before tilt_coords, kicks are independent of any tilt.

  if (set_hv1) then
    orbit%vec(2) = orbit%vec(2) + 0.5_rp * charge_dir * ele%value(hkick$)
    orbit%vec(4) = orbit%vec(4) + 0.5_rp * charge_dir * ele%value(vkick$)
    if (set_spn) call rotate_spin_given_field (orbit, sign_z_vel, 0.5_rp*B_factor*[ele%value(vkick$), -ele%value(hkick$), 0.0_rp], qrot = spin_qrot)
  endif

  ! Set: Tilt

  if (set_t .and. ele%key /= sbend$ .and. ele%value(tilt_tot$) /= 0) then
    call tilt_coords (ele%value(tilt_tot$), orbit%vec, mat6, make_matrix)
    if (set_spn) call rotate_spin ([0.0_rp, 0.0_rp, -ele%value(tilt_tot$)], orbit%spin, qrot = spin_qrot)
  endif

  ! Set: HV kicks for kickers and separators only.
  ! Note: Since this is applied after tilt_coords, kicks are dependent on any tilt.

  if (set_hv2) then
    if (ele%key == elseparator$) then
      rtc = 0.5_rp * abs(rel_tracking_charge) * sign(1, charge_of(orbit%species))
      orbit%vec(2) = orbit%vec(2) + rtc * ele%value(hkick$)
      orbit%vec(4) = orbit%vec(4) + rtc * ele%value(vkick$)
      if (set_spn .and. ele%value(e_field$) /= 0) call rotate_spin_given_field (orbit, sign_z_vel, &
                                      EL = 0.5_rp * [ele%value(hkick$), ele%value(vkick$), 0.0_rp] * ele%value(p0c$), qrot = spin_qrot)
    elseif (ele%key == hkicker$) then
      orbit%vec(2) = orbit%vec(2) + 0.5_rp * charge_dir * ele%value(kick$)
      if (set_spn) call rotate_spin_given_field (orbit, sign_z_vel, 0.5_rp*B_factor*[0.0_rp, -ele%value(kick$), 0.0_rp], qrot = spin_qrot)
    elseif (ele%key == vkicker$) then
      orbit%vec(4) = orbit%vec(4) + 0.5_rp * charge_dir * ele%value(kick$)
      if (set_spn) call rotate_spin_given_field (orbit, sign_z_vel, 0.5_rp*B_factor*[ele%value(kick$), 0.0_rp, 0.0_rp], qrot = spin_qrot)
    else
      orbit%vec(2) = orbit%vec(2) + 0.5_rp * charge_dir * ele%value(hkick$)
      orbit%vec(4) = orbit%vec(4) + 0.5_rp * charge_dir * ele%value(vkick$)
      if (set_spn) call rotate_spin_given_field (orbit, sign_z_vel, 0.5_rp*B_factor*[ele%value(vkick$), -ele%value(hkick$), 0.0_rp], qrot = spin_qrot)
    endif
  endif

!----------------------------------------------------------------
! Unset... 

else

  ! Unset: HV kicks for kickers and separators only.

  if (set_hv2) then
    if (ele%key == elseparator$) then
      rtc = 0.5_rp * abs(rel_tracking_charge) * sign(1, charge_of(orbit%species))
      orbit%vec(2) = orbit%vec(2) + rtc * ele%value(hkick$)
      orbit%vec(4) = orbit%vec(4) + rtc * ele%value(vkick$)
      if (set_spn .and. ele%value(e_field$) /= 0) call rotate_spin_given_field (orbit, sign_z_vel, &
                                         EL = 0.5_rp * [ele%value(hkick$), ele%value(vkick$), 0.0_rp] * ele%value(p0c$), qrot = spin_qrot)
    elseif (ele%key == hkicker$) then
      orbit%vec(2) = orbit%vec(2) + 0.5_rp * charge_dir * ele%value(kick$)
      if (set_spn) call rotate_spin_given_field (orbit, sign_z_vel, 0.5_rp*B_factor*[0.0_rp, -ele%value(kick$), 0.0_rp], qrot = spin_qrot)
    elseif (ele%key == vkicker$) then
      orbit%vec(4) = orbit%vec(4) + 0.5_rp * charge_dir * ele%value(kick$)
      if (set_spn) call rotate_spin_given_field (orbit, sign_z_vel, 0.5_rp*B_factor*[ele%value(kick$), 0.0_rp, 0.0_rp], qrot = spin_qrot)
    else
      orbit%vec(2) = orbit%vec(2) + 0.5_rp * charge_dir * ele%value(hkick$)
      orbit%vec(4) = orbit%vec(4) + 0.5_rp * charge_dir * ele%value(vkick$)
      if (set_spn) call rotate_spin_given_field (orbit, sign_z_vel, 0.5_rp*B_factor*[ele%value(vkick$), -ele%value(hkick$), 0.0_rp], qrot = spin_qrot)
    endif
  endif

  ! Unset: Tilt & Roll

  if (set_t .and. ele%key /= sbend$) then
    call tilt_coords (-ele%value(tilt_tot$), orbit%vec, mat6, make_matrix)
    if (set_spn) call rotate_spin ([0.0_rp, 0.0_rp, ele%value(tilt_tot$)], orbit%spin, qrot = spin_qrot)
  endif

  ! UnSet: HV kicks for quads, etc. but not hkicker, vkicker, elsep and kicker elements.
  ! HV kicks must come after z_offset but before any tilts are applied.
  ! Note: Change in %vel is NOT dependent upon energy since we are using
  ! canonical momentum.

  if (set_hv1) then
    orbit%vec(2) = orbit%vec(2) + 0.5_rp * charge_dir * ele%value(hkick$)
    orbit%vec(4) = orbit%vec(4) + 0.5_rp * charge_dir * ele%value(vkick$)
    if (set_spn) call rotate_spin_given_field (orbit, sign_z_vel, 0.5_rp*B_factor*[ele%value(vkick$), -ele%value(hkick$), 0.0_rp], qrot = spin_qrot)
  endif

  ! Unset: Offset and pitch

  if (has_orientation_attributes(ele)) then
    x_off = ele%value(x_offset_tot$)
    y_off = ele%value(y_offset_tot$)
    z_off = ele%value(z_offset_tot$)
    xp    = ele%value(x_pitch_tot$)
    yp    = ele%value(y_pitch_tot$)
    ref_tilt = ele%value(ref_tilt_tot$)

    if (ele%key == sad_mult$) then
      x_off = x_off + ele%value(x_offset_mult$)
      y_off = y_off + ele%value(y_offset_mult$)
    endif

    is_misaligned = (x_off /= 0 .or. y_off /= 0 .or. z_off /= 0 .or. xp /= 0 .or. yp /= 0 .or. &
                                (ele%key == sbend$ .and. (ref_tilt /= 0 .or. ele%value(roll$) /= 0)))
  endif

  if (is_misaligned) then
    position%r = [orbit%vec(1), orbit%vec(3), 0.0_rp]
    call mat_make_unit (position%w)

    if (ele%key == sbend$ .and. (ele%value(g$) /= 0 .or. ref_tilt /= 0 .or. ele%value(roll$) /= 0)) then

      position = bend_shift(position, ele%value(g$), ds_center)

      if (ref_tilt /= 0) then
        position = bend_shift(position, ele%value(g$), -L_half)
        ws = w_mat_for_tilt(ref_tilt)
        position%r = matmul(ws, position%r)
        position%w = matmul(ws, position%w)
        position = bend_shift(position, ele%value(g$), L_half, ref_tilt = ref_tilt)
      endif

      call ele_misalignment_L_S_calc(ele, L_mis, ws)
      position%r = matmul(ws, position%r) + L_mis
      position%w = matmul(ws, position%w)

      ! Coordinates when drifting to the bend ends are the coordinates at the end. 
      ! Coordinates when not drifting are the coordinates at the longitudinal position of the point.

      if (do_drift) then
        position = bend_shift(position, ele%value(g$), sign_z_vel*L_half, ref_tilt = ref_tilt)
        position%r(3) = position%r(3) + (sign_z_vel+1.0_rp)*L_half
      else
        if (ele%value(g$) == 0) then
          ds = -ds_center
        else
          ds = ele%value(rho$) * atan2(position%r(3)*ele%value(g$), 1.0_rp + position%r(1)*ele%value(g$))
        endif
        position = bend_shift(position, ele%value(g$), ds, ref_tilt = ref_tilt)
        position%r(3) = position%r(3) + L_half + ds
      endif

    ! Else not a bend or a bend with zero bending angle

    else
      position%r(3) = position%r(3) - ds_center
      call floor_angles_to_w_mat (xp, yp, 0.0_rp, w_mat = ws)
      position%r = matmul(ws, position%r)
      position%w = matmul(ws, position%w)
      position%r  = position%r + [x_off, y_off, z_off+L_half]
    endif

    if (ele%orientation == -1) position%r(3) = length - position%r(3)
    s_lab = position%r(3)

    pz = rel_p**2 - orbit%vec(2)**2 - orbit%vec(4)**2
    if (pz <= 0) then
      orbit%state = lost_pz_aperture$
    else
      p_vec0 = [orbit%vec(2), orbit%vec(4), sign_z_vel * sqrt(pz)]
      p_vec = matmul(position%w, p_vec0)
      orbit%vec(2:4:2) = p_vec(1:2)
      orbit%vec(1:3:2) = position%r(1:2)
      if (logic_option(.false., make_matrix)) call apply_offsets_to_matrix (p_vec0, p_vec, position%w, mat6)
    endif

    if (set_spn) orbit%spin = matmul(position%w, orbit%spin)
    if (present(spin_qrot)) spin_qrot = quat_mul(w_mat_to_quat(position%w), spin_qrot)
  endif    ! has nonzero offset or pitch

  if (.not. is_misaligned) then
    select case (ele%orientation)
    case (1);  s_lab = s_pos0
    case (-1); s_lab = length - s_pos0
    end select
  endif

  ! When drifting the reference particle does not move! That is, the reference time is the time the 
  ! referece particle reaches the *nominal* edge of the element and this is independent of any misalignments.

  if (do_drift .and. orbit%state == alive$) then
    call track_a_drift (orbit, s_target-s_lab, mat6, make_matrix, +1, include_ref_motion = .false.)
  endif

  !

  if (present(s_out)) then
    if (do_drift) then
      s_out = s_target
    else
      s_out = s_lab
    endif
  endif

endif

!--------------------------------------------------------------------------
contains

subroutine apply_offsets_to_matrix (p_vec0, p_vec, ww, mat6)

real(rp) p_vec0(3), p_vec(3), ww(3,3), mat6(6,6)
real(rp) dmat(6,6), pz_in, pz_out

!

pz_out = p_vec(3)
pz_in = p_vec0(3)

call mat_make_unit(dmat)

dmat(1,1) = ww(1,1) - ww(3,1) * p_vec(1) / pz_out
dmat(1,3) = ww(1,2) - ww(3,2) * p_vec(1) / pz_out

dmat(2,2) = ww(1,1) - ww(1,3) * p_vec0(1) / pz_in
dmat(2,4) = ww(1,2) - ww(1,3) * p_vec0(2) / pz_in
dmat(2,6) = ww(1,3) * rel_p / pz_in

dmat(3,1) = ww(2,1) - ww(3,1) * p_vec(2) / pz_out
dmat(3,3) = ww(2,2) - ww(3,2) * p_vec(2) / pz_out

dmat(4,2) = ww(2,1) - ww(2,3) * p_vec0(1) / pz_in
dmat(4,4) = ww(2,2) - ww(2,3) * p_vec0(2) / pz_in
dmat(4,6) = ww(2,3) * rel_p / pz_in

dmat(5,1) = ww(3,1) * rel_p / pz_out
dmat(5,3) = ww(3,2) * rel_p / pz_out

mat6 = matmul (dmat, mat6)

end subroutine apply_offsets_to_matrix

end subroutine offset_particle
