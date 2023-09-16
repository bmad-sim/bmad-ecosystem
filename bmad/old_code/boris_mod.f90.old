!+
! Module boris_mod
!
! Module to do Boris integration tracking.
!
! Reference: 
!  "Efficiency of a Boris-like Integration Scheme with Spatial Stepping", 
!   P. H. Stolz et al., Physical Review Special Topics.
!   5, 094001 (2002).
!
! When comparing the paper to this module remember that the paper uses a
! coordinate system:
!   (x, P_x, y, P_y, ct, U)
! While the Bmad coordinate system is:
!   (x, P_x/P_0, y, P_y/P_0, beta*c*dt, dP/P_0)
!-

module boris_mod

use em_field_mod
use fringe_mod

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine track1_boris (start_orb, ele, param, end_orb, err_flag, track, s_start, s_end)
!
! Subroutine to do Boris tracking. For more information on Boris tracking 
! see the boris_mod module documentation.
! 
! Modules needed:
!   use bmad
!
! Input: 
!   start_orb  -- Coord_struct: start_orbing coords.
!   ele        -- Ele_struct: Element to track through.
!     %tracking_method -- Determines which subroutine to use to calculate the 
!                         field. Note: BMAD does no supply em_field_custom.
!                           == custom$ then use em_field_custom
!                           /= custom$ then use em_field
!     %value(ds_step$) -- Step size.
!   param      -- lat_param_struct: Beam parameters.
!     %particle    -- Particle type [positron$, or electron$]
!   s_start    -- Real, optional: Starting point.
!   s_end      -- Real, optional: Ending point.
!
! Output:
!   end        -- Coord_struct: Ending coords.
!   err_flag   -- Logical: Set True if there is an error. False otherwise.
!   track      -- Track_struct, optional: Structure holding the track information.
!                  When tracking through multiple elements, the trajectory in an element
!                  is appended to the existing trajectory. To reset: Set track%n_pt = -1.
!-

subroutine track1_boris (start_orb, ele, param, end_orb, err_flag, track, s_start, s_end)

implicit none

type (coord_struct), intent(in) :: start_orb
type (coord_struct), intent(out) :: end_orb
type (ele_struct) ele
type (lat_param_struct) param
type (track_struct), optional :: track
type (fringe_edge_info_struct) fringe_info

real(rp), optional, intent(in) :: s_start, s_end
real(rp) s1, s2, ds, s, beta, s_edge_track, s_target
real(rp) dref_time, s0, ds_ref, beta0

integer i, n_step

character(16), parameter :: r_name = 'track1_boris'

logical err_flag, track_spin

! Boris is not able to handle a zero length element with a non-zero multipole.

if (ele%key /= patch$ .and. ele%value(l$) == 0) then
  call track_a_zero_length_element (start_orb, ele, param, end_orb, err_flag, track)
  return
endif

! init

err_flag = .true.

if (present(s_start)) then
  s1 = s_start
else
  s1 = 0
endif

if (present(s_end)) then
  s2 = s_end
else 
  s2 = ele%value(l$)
endif

! go to local coords

end_orb = start_orb
end_orb%s = s1 + ele%s_start

if (ele%key == patch$) then
  call track_a_patch (ele, end_orb, .false., s0, ds_ref)
  beta0 = ele%value(p0c$) / ele%value(e_tot$)
  end_orb%vec(5) = end_orb%vec(5) + ds_ref * end_orb%beta / beta0 + s0
else
  if (start_orb%direction == -1) then
    s0 = s1
    s1 = s2
    s2 = s0
  endif
endif

! If the element is using a hard edge model then need to stop at the hard edges
! to apply the appropriate hard edge kick.

call calc_next_fringe_edge (ele, s_edge_track, fringe_info, end_orb, .true.)

call compute_even_steps (ele%value(ds_step$), s2-s1, bmad_com%default_ds_step, ds, n_step)

call reference_energy_correction (ele, end_orb, first_track_edge$)

call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false.)

! if we are saving the trajectory then allocate enough space in the arrays

if (present(track)) then
  call save_a_step (track, ele, param, .true., end_orb, s1, .true.)
endif

! track through the body

s = s1

do 

  do
    if (abs(s - s_edge_track) > bmad_com%significant_length .or. .not. associated(fringe_info%hard_ele)) exit
    track_spin = (ele%spin_tracking_method == tracking$ .and. ele%field_calc == bmad_standard$)
    call apply_element_edge_kick (end_orb, fringe_info, ele, param, track_spin)
    call calc_next_fringe_edge (ele, s_edge_track, fringe_info, end_orb)
  enddo

  if (abs(s - s2) < bmad_com%significant_length) exit

  s_target = min(s2, s_edge_track)
  call compute_even_steps (ele%value(ds_step$), s_target-s, bmad_com%default_ds_step, ds, n_step)

  call track1_boris_partial (end_orb, ele, param, s, ds, end_orb)
  s = s + ds

  if (present(track)) call save_a_step (track, ele, param, .true., end_orb, s, .true.)
  
enddo

! back to lab coords

call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false.)

call reference_energy_correction (ele, end_orb, second_track_edge$)

! The z value computed in Boris tracking is off for elements where the particle changes energy is not 
! constant (see Boris tracking for more details). In this case make the needed correction.
! dref_time is reference time for transversing the element under the assumption, used by Boris tracking, that 
! the reference velocity is constant and equal to the velocity at the final enegy.

beta0 = ele%value(p0c$) / ele%value(e_tot$)
dref_time = ele%value(l$) / (beta0 * c_light)
end_orb%vec(5) = end_orb%vec(5) + (ele%value(delta_ref_time$) - dref_time) * end_orb%beta * c_light

err_flag = .false.

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine track1_boris_partial (start, ele, param, s, ds, end)
!
! Subroutine to track 1 step using boris tracking.
! This subroutine is used by track1_boris.
!
! Note: Coordinates are with respect to the element coordinate frame.
!
! For more information on Boris tracking see the boris_mod module documentation.
!
! Modules needed:
!   use bmad
!   use em_field_mod
!
! Input:
!   start -- Coord_struct: Starting coordinates.
!   ele   -- Ele_struct: Element that we are tracking through.
!   param -- lat_param_struct: 
!     %spin_tracking_on -- If True then also track the spin
!   s     -- Real(rp): Starting point relative to element beginning.
!   t     -- Real(rp): Particle time
!   ds    -- Real(rp): step size
!
! Output:
!   end   -- Coord_struct: Ending coordinates.
!-

subroutine track1_boris_partial (start, ele, param, s, ds, end)

implicit none

type (ele_struct) :: ele
type (lat_param_struct) param
type (coord_struct) :: start, end  
type (em_field_struct) :: field

real(rp), intent(in) :: s, ds
real(rp) :: f, p_z, d2, alpha, dxv, dyv, dt2_f, charge, U_tot, p_tot, ds2
real(rp) :: r(3,3), w(3), ex, ey, ex2, ey2, exy, bz, bz2, mass, old_beta, beta
real(rp) :: p2, dt, beta_ref, p2_z

character(*), parameter :: r_name = 'track1_boris_partial'

!


charge = charge_of(start%species)
mass = mass_of(start%species) / ele%value(p0c$)

end = start
ds2 = end%direction * ds / 2
beta_ref = ele%value(p0c$) / ele%value(e_tot$)

! 1) Push the position 1/2 step

p_tot = 1 + end%vec(6)
p2_z = p_tot**2 - end%vec(2)**2 - end%vec(4)**2
if (p2_z < 0 .or. p_tot < 0) then
  end%state = lost_pz_aperture$
  return
endif
p_z = sqrt(p2_z) * end%direction * ele%orientation
dt2_f = abs(ds2 / p_z)
U_tot = sqrt (p_tot**2 + mass**2)
old_beta = p_tot / U_tot  ! particle velocity: v/c
dt = dt2_f * p_tot / (old_beta * c_light)

end%vec(1) = end%vec(1) + dt2_f * end%vec(2) 
end%vec(3) = end%vec(3) + dt2_f * end%vec(4)
end%vec(5) = end%vec(5) + (ds2 * old_beta/beta_ref - dt2_f * p_tot) 

end%s = end%s + ds2
end%t = end%t + dt

! 2) Evaluate the fields.

call em_field_calc (ele, param, s+ds2, end, .true., field)

! 2.5) Push the spin 1/2 step
! This uses the momentum at the beginning and the fields at (ds2)

if (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$) then
  call rotate_spin_a_step (end, field, ele, ds2)
endif

! 3) Push the momenta a 1/2 step using only the "b" term.

f = ds2 * charge * c_light / ele%value(p0c$)

end%vec(2) = end%vec(2) - field%B(2) * f
end%vec(4) = end%vec(4) + field%B(1) * f
U_tot = U_tot + field%e(3) * f / c_light
p_tot = sqrt (U_tot**2 - mass**2)
p2_z = p_tot**2 - end%vec(2)**2 - end%vec(4)**2
if (p2_z < 0 .or. p_tot < 0) then
  end%state = lost_pz_aperture$
  return
endif
p_z = sqrt(p2_z) * end%direction * ele%orientation

! 4) Push the momenta a full step using the "R" matrix.

d2 = ds2 * charge * c_light / (p_z * ele%value(p0c$)) 

if (field%e(1) == 0 .and. field%e(2) == 0) then
  if (field%B(3) /= 0) then
    d2 = d2 * field%B(3)
    alpha = 2 * d2 / (1 + d2**2)
    dxv = -d2 * end%vec(2) + end%vec(4)
    dyv = -end%vec(2) - d2 * end%vec(4)
    end%vec(2) = end%vec(2) + alpha * dxv
    end%vec(4) = end%vec(4) + alpha * dyv
  endif
else
  ex = field%e(1) / c_light;     ex2 = ex**2
  ey = field%e(2) / c_light;     ey2 = ey**2
  bz = field%B(3);               bz2 = bz**2
  exy = ex * ey
  alpha = 2 * d2 / (1 + d2**2 * (bz2 - ex2 - ey2))
  r(1,1:3) = [d2 * (ex2 - bz2), bz + d2*exy,      ex + d2*bz*ey    ]
  r(2,1:3) = [-bz + d2*exy,     d2 * (ey2 - bz2), ey - d2*bz*ex    ]
  r(3,1:3) = [ex - d2*bz*ey,    ey + d2*bz*ex,    d2 * (ex2 + ey2) ]

  w = [end%vec(2), end%vec(4), U_tot ]
  w = w + alpha * matmul(r, w)
  end%vec(2:4:2) = w(1:2); U_tot = w(3)
endif

! 5) Push the momenta a 1/2 step using only the "b" term.

end%vec(2) = end%vec(2) - field%B(2) * f
end%vec(4) = end%vec(4) + field%B(1) * f
U_tot = U_tot + field%e(3) * f / c_light
p_tot = sqrt (U_tot**2 - mass**2)

! Since beta changes in steps 2-5 we need to update vec(5)

beta = p_tot / U_tot
end%vec(5) = end%vec(5) * beta / old_beta

! 6) Push the position 1/2 step.

p_z = sqrt(p_tot**2 - end%vec(2)**2 - end%vec(4)**2)
dt2_f = ds2 / p_z

end%vec(1) = end%vec(1) + dt2_f * end%vec(2) 
end%vec(3) = end%vec(3) + dt2_f * end%vec(4)
end%vec(5) = end%vec(5) + ds2 * (beta/beta_ref - p_tot/p_z) 
end%vec(6) = p_tot - 1

end%s = end%s + ds2
dt = ds2 * p_tot / (p_z * beta * c_light)
end%t = end%t + dt
end%beta = beta

! 6.5) Push the spin 1/2 step

if (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$) then
  call rotate_spin_a_step (end, field, ele, ds2)
endif

end subroutine

end module
