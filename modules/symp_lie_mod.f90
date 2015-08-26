module symp_lie_mod

use bmad_struct
use bmad_interface
use make_mat6_mod
use em_field_mod
use random_mod
use track1_mod

type save_coef_struct
  real(rp) :: coef = 0, dx_coef = 0, dy_coef = 0
end type

type wiggler_computations_struct
  type (save_coef_struct) Ax, dint_Ax_dy, Ay, dint_Ay_dx, dAz_dx, dAz_dy
  real(rp) :: c_x = 0, s_x = 0, c_y = 0, s_y = 0, c_z = 0, s_z = 0, s_over_kx = 0, one_minus_c_over_ky
end type

private save_coef_struct, wiggler_computations_struct

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine symp_lie_bmad (ele, param, start_orb, end_orb, calc_mat6, track, offset_ele)
!
! Subroutine to track through an element (which gives the 0th order map) 
! and optionally make the 6x6 transfer matrix (1st order map) as well.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele        -- Ele_struct: Element with transfer matrix
!   param      -- lat_param_struct: Parameters are needed for some elements.
!   start_orb  -- Coord_struct: Coordinates at the beginning of element. 
!   calc_mat6  -- Logical: If True then make the 6x6 transfer matrix.
!   offset_ele -- Logical, optional: Offset the element using ele%value(x_offset$), etc.
!                   Default is True.
!
! Output:
!   ele        -- Ele_struct: Element with transfer matrix.
!     %mat6(6,6)  -- 6x6 transfer matrix.
!     %vec0(6)    -- 0th order part of the transfer matrix.
!   end_orb    -- Coord_struct: Coordinates at the end of element.
!   track      -- Track_struct, optional: Structure holding the track information.
!                   When tracking through multiple elements, the trajectory in an element
!                   is appended to the existing trajectory. To reset: Set track%n_pt = -1.
!-

subroutine symp_lie_bmad (ele, param, start_orb, end_orb, calc_mat6, track, offset_ele)

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: field_ele
type (ele_pointer_struct), allocatable :: field_eles(:)
type (coord_struct) :: start_orb, end_orb, start2_orb
type (lat_param_struct)  param
type (track_struct), optional :: track
type (wig_term_struct), pointer :: wig_term(:)

type (wiggler_computations_struct), allocatable, target :: tm(:)
type (wig_term_struct), pointer :: wt

real(rp) rel_E, rel_E2, rel_E3, ds, ds2, s, m6(6,6), kmat6(6,6), Ax_saved, Ay_saved
real(rp) g_x, g_y, k1, k1_norm, k1_skew, x_q, y_q, ks_tot_2, ks, dks_ds, z_patch
real(rp), pointer :: mat6(:,:)
real(rp), parameter :: z0 = 0, z1 = 1
real(rp) gamma_0, fact_d, fact_f, this_ran, g2, g3
real(rp) dE_p, dpx, dpy, mc2, z_offset, orientation, rel_tracking_charge, charge_dir
real(rp), parameter :: rad_fluct_const = 55 * classical_radius_factor * h_bar_planck * c_light / (24 * sqrt_3)
real(rp), allocatable :: dz_offset(:)

integer i, n_step, n_field

integer num_wig_terms  ! number of wiggler terms

logical calc_mat6, calculate_mat6, err, do_offset
logical, optional :: offset_ele

character(16) :: r_name = 'symp_lie_bmad'

! init

calculate_mat6 = (calc_mat6 .or. synch_rad_com%i_calc_on)

do_offset = logic_option (.true., offset_ele)
rel_E = (1 + start_orb%vec(6))
rel_E2 = rel_E**2
rel_E3 = rel_E**3

start2_orb = start_orb
end_orb = start_orb
end_orb%s = ele%s - ele%value(l$)
mat6 => ele%mat6

err = .false.

orientation = ele%orientation * start_orb%direction
rel_tracking_charge = relative_tracking_charge(start_orb, param)
charge_dir = rel_tracking_charge * orientation

! element offset 

if (calculate_mat6) then
  call drift_mat6_calc (mat6, ele%value(z_offset_tot$), ele, param, end_orb)
endif

if (do_offset) call offset_particle (ele, param, set$, end_orb)

! init

call compute_even_steps (ele%value(ds_step$), ele%value(l$), bmad_com%default_ds_step, ds, n_step)
ds2 = ds / 2

s = 0   ! longitudianl position

! radiation damping and fluctuations...
! The same kick is applied over the entire wiggler to save time.

if ((bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on)) then

  mc2 = mass_of(param%particle)
  gamma_0 = ele%value(e_tot$) / mass_of(param%particle)

  fact_d = 0
  if (bmad_com%radiation_damping_on) fact_d = 2 * classical_radius_factor * gamma_0**3 * ds / (3 * mc2)

  fact_f = 0
  if (bmad_com%radiation_fluctuations_on) then
    fact_f = sqrt(rad_fluct_const * ds * gamma_0**5) / mc2 
  endif

endif

!------------------------------------------------------------------
! select the element

select case (ele%key)

!------------------------------------------------------------------
! Wiggler

Case (wiggler$, undulator$)

  call get_field_ele_list (ele, field_eles, dz_offset, n_field)
  do i = 1, n_field
    field_ele => field_eles(i)%ele
    if (field_ele%key == ele%key) exit
  enddo

  wig_term => field_ele%wig%term
  z_offset = dz_offset(i)

  num_wig_terms = size(wig_term)

  allocate (tm(num_wig_terms))

  call calc_wig_coefs (calculate_mat6)
  call update_wig_y_terms (err); if (err) return
  call update_wig_x_terms (err); if (err) return
  call update_wig_s_terms

  if (present(track)) then
    call save_this_track_pt (s)
  endif

  ! Symplectic integration gives a kick to the physical momentum (but not the canonical momentum) at the ends of an 
  ! element if there is a finite vector potential. This is necessary for symplecticity but can complicate comparisons
  ! with runge_kutta. That is, only set bmad_com%convert_to_kinetic_momentum = True for testing purposes.

  if (bmad_com%convert_to_kinetic_momentum) then
    end_orb%vec(2) = end_orb%vec(2) + Ax()
    end_orb%vec(4) = end_orb%vec(4) + Ay()
  endif

  ! loop over all steps

  do i = 1, n_step

    ! s half step

    s = s + ds2
    call update_wig_s_terms

    ! Drift_1 = (P_x - Ax)^2 / (2 * (1 + dE))

    call apply_wig_exp_int_ax (-1, calculate_mat6)

    call apply_p_x (calculate_mat6)
    call update_wig_x_terms (err); if (err) return

    call apply_wig_exp_int_ax (+1, calculate_mat6)

    ! Drift_2 = (P_y - Ay)^2 / (2 * (1 + dE))

    call apply_wig_exp_int_ay (-1, calculate_mat6)

    call apply_p_y (calculate_mat6)
    call update_wig_y_terms (err); if (err) return

    call apply_wig_exp_int_ay (+1, calculate_mat6)

    ! Kick = Az

    dpx = dAz_dx()
    dpy = dAz_dy()
    end_orb%vec(2) = end_orb%vec(2) + ds * dpx
    end_orb%vec(4) = end_orb%vec(4) + ds * dpy

    call radiation_kick()

    if (calculate_mat6) then
      mat6(2,1:6) = mat6(2,1:6) + ds * (dAz_dx__dx() * mat6(1,1:6) + dAz_dx__dy() * mat6(3,1:6))
      mat6(4,1:6) = mat6(4,1:6) + ds * (dAz_dy__dx() * mat6(1,1:6) + dAz_dy__dy() * mat6(3,1:6))
    endif 

    ! Drift_2

    call apply_wig_exp_int_ay (-1, calculate_mat6)

    call apply_p_y (calculate_mat6)
    call update_wig_y_terms (err); if (err) return

    call apply_wig_exp_int_ay (+1, calculate_mat6)

    ! Drift_1

    call apply_wig_exp_int_ax (-1, calculate_mat6)

    call apply_p_x (calculate_mat6)
    call update_wig_x_terms (err); if (err) return

    call apply_wig_exp_int_ax (+1, calculate_mat6)

    ! s half step

    s = s + ds2
    call update_wig_s_terms

    if (present(track)) call save_this_track_pt (s)

  enddo

  ! Correction for finite vector potential at exit end.

  if (bmad_com%convert_to_kinetic_momentum) then
    end_orb%vec(2) = end_orb%vec(2) - Ax()
    end_orb%vec(4) = end_orb%vec(4) - Ay()
  endif

!----------------------------------------------------------------------------
! rf cavity

case (lcavity$, rfcavity$)

  ! loop over all steps

  do i = 1, n_step

    ! s half step

    s = s + ds2
!    call rf_drift1 (calculate_mat6)
!    call rf_drift2 (calculate_mat6)
!    call rf_kick (calculate_mat6)
    call radiation_kick()
!    call rf_drift2 (calculate_mat6)
!    call rf_drift1 (calculate_mat6)

    s = s + ds2

    if (present(track)) call save_this_track_pt (s)

  enddo

!----------------------------------------------------------------------------
! solenoid, quadrupole, sol_quad, or bend_sol_quad
! Notice that we don't have to worry about a solenoid hard edge since we are
! using canonical coords so the hard edge is automatically included.

case (bend_sol_quad$, solenoid$, quadrupole$, sol_quad$)

  g_x = 0
  g_y = 0
  x_q = 0
  y_q = 0
  dks_ds = 0
  k1_norm = 0
  k1_skew = 0
  ks = 0

  select case (ele%key)
  case (bend_sol_quad$)
    g_x = ele%value(g$) * cos (ele%value(bend_tilt$)) * charge_dir
    g_y = ele%value(g$) * sin (ele%value(bend_tilt$)) * charge_dir
    k1_norm = ele%value(k1$) * cos (2 * ele%value(quad_tilt$)) * charge_dir
    k1_skew = ele%value(k1$) * sin (2 * ele%value(quad_tilt$)) * charge_dir
    x_q = ele%value(x_quad$)
    y_q = ele%value(y_quad$)
    ks = ele%value(ks$) * rel_tracking_charge
    dks_ds = ele%value(dks_ds$) * charge_dir
  case (solenoid$)
    ks = ele%value(ks$) * rel_tracking_charge
  case (quadrupole$)
    k1_norm = ele%value(k1$) * charge_dir 
  case (sol_quad$)
    k1_norm = ele%value(k1$) * charge_dir
    ks = ele%value(ks$) * rel_tracking_charge
  end select

  call hard_multipole_edge_kick (ele, param, first_track_edge$, end_orb, mat6, calculate_mat6)
  call soft_quadrupole_edge_kick (ele, param, first_track_edge$, end_orb, mat6, calculate_mat6)

  ! loop over all steps

  do i = 1, n_step

    s = s + ds2
    ks_tot_2 = (ks + dks_ds * s) / 2

    call bsq_drift1 (calculate_mat6)
    call bsq_drift2 (calculate_mat6)
    call bsq_kick (calculate_mat6)
    call radiation_kick()
    call bsq_drift2 (calculate_mat6)
    call bsq_drift1 (calculate_mat6)

    s = s + ds2
    ks_tot_2 = (ks + dks_ds * s) / 2

    if (present(track)) call save_this_track_pt (s)

  enddo

  call soft_quadrupole_edge_kick (ele, param, second_track_edge$, end_orb, mat6, calculate_mat6)
  call hard_multipole_edge_kick (ele, param, second_track_edge$, end_orb, mat6, calculate_mat6)

!----------------------------------------------------------------------------
! unknown element

case default

  call out_io (s_fatal$, r_name, 'TRACKING NOT YET IMPLEMENTED FOR: ' // key_name(ele%key), &
                                 'FOR ELEMENT: ', ele%name)

end select

!----------------------------------------------------------------------------
! element offset

if (calculate_mat6) then
  call drift_mat6_calc (m6, -ele%value(z_offset_tot$), ele, param, end_orb)
  mat6(1,1:6) = mat6(1,1:6) + m6(1,2) * mat6(2,1:6) + m6(1,6) * mat6(6,1:6)
  mat6(3,1:6) = mat6(3,1:6) + m6(3,4) * mat6(4,1:6) + m6(3,6) * mat6(6,1:6)
  mat6(5,1:6) = mat6(5,1:6) + m6(5,2) * mat6(2,1:6) + m6(5,4) * mat6(4,1:6) + m6(5,6) * mat6(6,1:6)

  if (ele%value(tilt_tot$) /= 0) call tilt_mat6 (mat6, ele%value(tilt_tot$))
  call mat6_add_pitch (ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%orientation, mat6)
endif

if (do_offset) call offset_particle (ele, param, unset$, end_orb)

! Correct z-position for wigglers, etc. 

z_patch = ele%value(delta_ref_time$) * c_light * end_orb%beta - ele%value(l$)
end_orb%vec(5) = end_orb%vec(5) + z_patch
end_orb%t = start2_orb%t + ele%value(delta_ref_time$) + (start2_orb%vec(5) - end_orb%vec(5)) / &
                                                                            (end_orb%beta * c_light)

if (calculate_mat6) then
  mat6(5,6) = mat6(5,6) + ele%value(delta_ref_time$) * c_light * &
                (mass_of(param%particle)/ele%value(p0c$))**2 * (end_orb%beta / (1 + end_orb%vec(6)))**3
endif

! calc vec0

if (calculate_mat6) then
  ele%vec0(1:5) = end_orb%vec(1:5) - matmul (mat6(1:5,1:6), start_orb%vec)
  ele%vec0(6) = 0
endif

!

end_orb%s = ele%s
end_orb%p0c = ele%value(p0c$)

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
contains

subroutine err_set (err, plane)

logical err
integer plane

!

call out_io (s_warn$, r_name, &
                'FLOATING OVERFLOW IN TRACKING:' // ele%name, &
                'PARTICLE WILL BE TAGGED AS LOST.')

if (plane == x_plane$) then
  end_orb%vec(1) = sign(2 * bmad_com%max_aperture_limit, end_orb%vec(1))
  if (end_orb%vec(1) > 0) then ; end_orb%state = lost_pos_x_aperture$
  else;                          end_orb%state = lost_neg_x_aperture$
  endif
else
  end_orb%vec(3) = sign(2 * bmad_com%max_aperture_limit, end_orb%vec(3))
  if (end_orb%vec(3) > 0) then ; end_orb%state = lost_pos_y_aperture$
  else;                          end_orb%state = lost_neg_y_aperture$
  endif
endif

err = .true.

end subroutine err_set

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine save_this_track_pt (s)

real(rp) s, s_sav
integer ix

!

call save_a_step (track, ele, param, .true., s, end_orb, s_sav)

if (calculate_mat6) then
  ix = track%n_pt
  track%map(ix)%mat6 = mat6
  if (ele%value(tilt_tot$) /= 0) call tilt_mat6 (track%map(ix)%mat6, ele%value(tilt_tot$))
  call mat6_add_pitch (ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%orientation, track%map(ix)%mat6)
  track%map(ix)%vec0 = track%orb(ix)%vec - matmul (track%map(ix)%mat6, start_orb%vec)
endif

if (bmad_com%convert_to_kinetic_momentum) then
  ix = track%n_pt
  track%orb(ix)%vec(2) = track%orb(ix)%vec(2) - Ax()
  track%orb(ix)%vec(4) = track%orb(ix)%vec(4) - Ay()
endif

end subroutine save_this_track_pt

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine apply_p_x (do_mat6)

logical do_mat6

end_orb%vec(1) = end_orb%vec(1) + ds2 * end_orb%vec(2) / rel_E
end_orb%vec(5) = end_orb%vec(5) - ds2 * end_orb%vec(2)**2 / (2*rel_E2)
end_orb%s = end_orb%s + ds2

if (do_mat6) then
  mat6(1,1:6) = mat6(1,1:6) + (ds2 / rel_E)               * mat6(2,1:6) - (ds2*end_orb%vec(2)/rel_E2)    * mat6(6,1:6) 
  mat6(5,1:6) = mat6(5,1:6) - (ds2*end_orb%vec(2)/rel_E2) * mat6(2,1:6) + (ds2*end_orb%vec(2)**2/rel_E3) * mat6(6,1:6)
endif

end subroutine apply_p_x 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine apply_p_y (do_mat6)

logical do_mat6

end_orb%vec(3) = end_orb%vec(3) + ds2 * end_orb%vec(4) / rel_E
end_orb%vec(5) = end_orb%vec(5) - ds2 * end_orb%vec(4)**2 / (2*rel_E2)

if (do_mat6) then
  mat6(3,1:6) = mat6(3,1:6) + (ds2 / rel_E)               * mat6(4,1:6) - (ds2*end_orb%vec(4)/rel_E2)    * mat6(6,1:6) 
  mat6(5,1:6) = mat6(5,1:6) - (ds2*end_orb%vec(4)/rel_E2) * mat6(4,1:6) + (ds2*end_orb%vec(4)**2/rel_E3) * mat6(6,1:6)
endif      

end subroutine apply_p_y

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine rf_drift1 (do_mat6)

logical do_mat6

! Drift_1 = (P_x - Ax)**2 / (2 * (1 + dE))

end_orb%vec(2) = end_orb%vec(2) 
end_orb%vec(4) = end_orb%vec(4) 

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) 
  mat6(4,1:6) = mat6(4,1:6) 
endif      

!

call apply_p_x (do_mat6)

!

end_orb%vec(2) = end_orb%vec(2) 
end_orb%vec(4) = end_orb%vec(4) 

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6)
  mat6(4,1:6) = mat6(4,1:6)
endif  

end subroutine rf_drift1

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine bsq_drift1 (do_mat6)

logical do_mat6

! Drift_1 = (P_x - Ax)**2 / (2 * (1 + dE))

end_orb%vec(2) = end_orb%vec(2) + end_orb%vec(3) * ks_tot_2   !  vec(2) - Ax
end_orb%vec(4) = end_orb%vec(4) + end_orb%vec(1) * ks_tot_2   !  vec(4) - dint_Ax_dy

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) + ks_tot_2 * mat6(3,1:6)
  mat6(4,1:6) = mat6(4,1:6) + ks_tot_2 * mat6(1,1:6)
endif      

!

call apply_p_x (do_mat6)

!

end_orb%vec(2) = end_orb%vec(2) - end_orb%vec(3) * ks_tot_2   !  vec(2) + Ax
end_orb%vec(4) = end_orb%vec(4) - end_orb%vec(1) * ks_tot_2   !  vec(4) + dint_Ax_dy

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) - ks_tot_2 * mat6(3,1:6)
  mat6(4,1:6) = mat6(4,1:6) - ks_tot_2 * mat6(1,1:6)
endif  

end subroutine bsq_drift1

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine bsq_drift2 (do_mat6)

logical do_mat6

! Drift_2 = (P_y - Ay)**2 / (2 * (1 + dE))

end_orb%vec(2) = end_orb%vec(2) - end_orb%vec(3) * ks_tot_2   !  vec(2) - dint_Ay_dx
end_orb%vec(4) = end_orb%vec(4) - end_orb%vec(1) * ks_tot_2   !  vec(4) - Ay

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) - ks_tot_2 * mat6(3,1:6)
  mat6(4,1:6) = mat6(4,1:6) - ks_tot_2 * mat6(1,1:6)
endif      

!

call apply_p_y (do_mat6)

!

end_orb%vec(2) = end_orb%vec(2) + end_orb%vec(3) * ks_tot_2   !  vec(2) + dint_Ay_dx
end_orb%vec(4) = end_orb%vec(4) + end_orb%vec(1) * ks_tot_2   !  vec(4) + Ay

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) + ks_tot_2 * mat6(3,1:6)
  mat6(4,1:6) = mat6(4,1:6) + ks_tot_2 * mat6(1,1:6)
endif  

end subroutine bsq_drift2

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine bsq_kick (do_mat6)

logical do_mat6

dpx = k1_norm * (x_q - end_orb%vec(1)) - k1_skew * end_orb%vec(3) - g_x
end_orb%vec(2) = end_orb%vec(2) + ds * dpx  ! dAz_dx
              
dpy = k1_norm * (end_orb%vec(3) - y_q) - k1_skew * end_orb%vec(1) - g_y
end_orb%vec(4) = end_orb%vec(4) + ds * dpy  ! dAz_dy

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) - ds * k1_norm * mat6(1,1:6) - ds * k1_skew * mat6(3,1:6)
  mat6(4,1:6) = mat6(4,1:6) - ds * k1_skew * mat6(1,1:6) + ds * k1_norm * mat6(3,1:6)
endif 

end subroutine bsq_kick

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine apply_wig_exp_int_ax (sgn, do_mat6)

integer sgn
logical do_mat6

Ax_saved = Ax()
end_orb%vec(2) = end_orb%vec(2) + sgn * Ax_saved
end_orb%vec(4) = end_orb%vec(4) + sgn * dint_Ax_dy()

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) + sgn * (Ax__dx()         * mat6(1,1:6) + Ax__dy()         * mat6(3,1:6))
  mat6(4,1:6) = mat6(4,1:6) + sgn * (dint_Ax_dy__dx() * mat6(1,1:6) + dint_Ax_dy__dy() * mat6(3,1:6))
endif      

end subroutine apply_wig_exp_int_ax

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine apply_wig_exp_int_ay (sgn, do_mat6)

integer sgn
logical do_mat6

Ay_saved = Ay()
end_orb%vec(2) = end_orb%vec(2) + sgn * dint_Ay_dx()
end_orb%vec(4) = end_orb%vec(4) + sgn * Ay_saved

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) + sgn * (dint_Ay_dx__dx() * mat6(1,1:6) + dint_Ay_dx__dy() * mat6(3,1:6))
  mat6(4,1:6) = mat6(4,1:6) + sgn * (Ay__dx()         * mat6(1,1:6) + Ay__dy()         * mat6(3,1:6))
endif      

end subroutine apply_wig_exp_int_ay

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine calc_wig_coefs (do_mat6)

type (wiggler_computations_struct), pointer :: tmj
real(rp) factor, coef
integer j
logical do_mat6

factor = rel_tracking_charge * c_light * field_ele%value(polarity$) / ele%value(p0c$)

do j = 1, num_wig_terms
  wt => wig_term(j)
  tmj => tm(j)
  coef = factor * wt%coef 

  select case (wt%type)
  case (hyper_y$)
    tmj%Ax%coef         =  coef * wt%kz / wt%ky**2
    tmj%Ay%coef         =  0
    tmj%dint_Ax_dy%coef =  tmj%Ax%coef * wt%ky               ! / wt%kx
    tmj%dint_Ay_dx%coef =  0
    tmj%dAz_dx%coef     = -coef * orientation * (wt%kx / wt%ky)**2
    tmj%dAz_dy%coef     = -coef * orientation * wt%kx / wt%ky
    tmj%Ax%dx_coef         = -tmj%Ax%coef * wt%kx
    tmj%Ax%dy_coef         =  tmj%Ax%coef * wt%ky
    tmj%Ay%dx_coef         =  0
    tmj%Ay%dy_coef         =  0
    tmj%dint_Ax_dy%dx_coef =  tmj%Ax%coef * wt%ky
    tmj%dint_Ax_dy%dy_coef =  tmj%dint_Ax_dy%coef * wt%ky   ! / wt%kx
    tmj%dint_Ay_dx%dx_coef =  0
    tmj%dint_Ay_dx%dy_coef =  0
    tmj%dAz_dx%dx_coef     = -tmj%dAz_dx%coef * wt%kx
    tmj%dAz_dx%dy_coef     =  tmj%dAz_dx%coef * wt%ky
    tmj%dAz_dy%dx_coef     =  tmj%dAz_dy%coef * wt%kx
    tmj%dAz_dy%dy_coef     =  tmj%dAz_dy%coef * wt%ky

  case (hyper_xy$)
    tmj%Ax%coef         =  coef * wt%ky / wt%kz**2
    tmj%Ay%coef         = -coef * wt%kx / wt%kz**2
    tmj%dint_Ax_dy%coef =  tmj%Ax%coef * wt%ky              ! / wt%kx
    tmj%dint_Ay_dx%coef =  tmj%Ay%coef * wt%kx              ! / wt%ky
    tmj%dAz_dx%coef     =  0
    tmj%dAz_dy%coef     =  0
    tmj%Ax%dx_coef         =  tmj%Ax%coef * wt%kx
    tmj%Ax%dy_coef         =  tmj%Ax%coef * wt%ky
    tmj%Ay%dx_coef         =  tmj%Ay%coef * wt%kx
    tmj%Ay%dy_coef         =  tmj%Ay%coef * wt%ky
    tmj%dint_Ax_dy%dx_coef =  tmj%Ax%coef * wt%ky
    tmj%dint_Ax_dy%dy_coef =  tmj%dint_Ax_dy%coef * wt%ky   ! / wt%kx
    tmj%dint_Ay_dx%dx_coef =  tmj%dint_Ay_dx%coef * wt%kx   ! / wt%ky
    tmj%dint_Ay_dx%dy_coef =  tmj%Ay%coef * wt%kx
    tmj%dAz_dx%dx_coef     =  0
    tmj%dAz_dx%dy_coef     =  0
    tmj%dAz_dy%dx_coef     =  0
    tmj%dAz_dy%dy_coef     =  0

  case (hyper_x$)
    tmj%Ax%coef         =  0
    tmj%Ay%coef         = -coef * wt%kz / wt%kx**2 
    tmj%dint_Ax_dy%coef =  0
    tmj%dint_Ay_dx%coef =  tmj%Ay%coef * wt%kx             ! / wt%ky
    tmj%dAz_dx%coef     = -coef * orientation * wt%ky / wt%kx
    tmj%dAz_dy%coef     =  coef * orientation * (wt%ky / wt%kx)**2
    tmj%Ax%dx_coef         =  0
    tmj%Ax%dy_coef         =  0
    tmj%Ay%dx_coef         =  tmj%Ay%coef * wt%kx
    tmj%Ay%dy_coef         =  tmj%Ay%coef * wt%ky
    tmj%dint_Ax_dy%dx_coef =  0
    tmj%dint_Ax_dy%dy_coef =  0
    tmj%dint_Ay_dx%dx_coef =  tmj%dint_Ay_dx%coef * wt%kx  ! / wt%ky
    tmj%dint_Ay_dx%dy_coef =  tmj%Ay%coef * wt%kx
    tmj%dAz_dx%dx_coef     =  tmj%dAz_dx%coef * wt%kx
    tmj%dAz_dx%dy_coef     = -tmj%dAz_dx%coef * wt%ky
    tmj%dAz_dy%dx_coef     =  tmj%dAz_dy%coef * wt%kx
    tmj%dAz_dy%dy_coef     =  tmj%dAz_dy%coef * wt%ky

  case (hyper_y_old$, hyper_xy_old$, hyper_x_old$)
    if (abs(wt%kx) < 1d-30) wt%kx = 1d-30   ! To prevent 0/0
    tmj%Ax%coef         =  0
    tmj%Ay%coef         = -coef * wt%kz / (wt%kx * wt%ky) 
    tmj%dint_Ax_dy%coef =  0
    tmj%dint_Ay_dx%coef =  tmj%Ay%coef * wt%kx             ! / wt%ky
    tmj%dAz_dx%coef     = -coef * orientation
    tmj%dAz_dy%coef     = -coef * orientation * wt%ky / wt%kx
    tmj%Ax%dx_coef         =  0
    tmj%Ax%dy_coef         =  0
    tmj%Ay%dx_coef         =  tmj%Ay%coef * wt%kx
    tmj%Ay%dy_coef         =  tmj%Ay%coef * wt%ky
    tmj%dint_Ax_dy%dx_coef =  0
    tmj%dint_Ax_dy%dy_coef =  0
    tmj%dint_Ay_dx%dx_coef =  tmj%dint_Ay_dx%coef * wt%kx  ! / wt%ky
    tmj%dint_Ay_dx%dy_coef =  tmj%Ay%coef * wt%kx
    tmj%dAz_dx%dx_coef     =  tmj%dAz_dx%coef * wt%kx
    tmj%dAz_dx%dy_coef     =  tmj%dAz_dx%coef * wt%ky
    tmj%dAz_dy%dx_coef     =  tmj%dAz_dy%coef * wt%kx
    tmj%dAz_dy%dy_coef     =  tmj%dAz_dy%coef * wt%ky

    select case (wt%type)
    case (hyper_y_old$)
      tmj%dint_Ay_dx%dx_coef =  -tmj%dint_Ay_dx%dx_coef
      tmj%dAz_dx%dx_coef     =  -tmj%dAz_dx%dx_coef
    case (hyper_x_old$)
      tmj%dAz_dy%coef        =  -tmj%dAz_dy%coef
      tmj%dAz_dx%dy_coef     =  -tmj%dAz_dx%dy_coef
      tmj%dAz_dy%dy_coef     =  -tmj%dAz_dy%dy_coef
      tmj%dAz_dy%dx_coef     =  -tmj%dAz_dy%dx_coef
    end select

  end select
enddo

end subroutine calc_wig_coefs

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine update_wig_x_terms (err)

type (wiggler_computations_struct), pointer :: tmj

real(rp) kxx, kxx0
integer j
logical err

!

do j = 1, num_wig_terms
  wt => wig_term(j)
  tmj => tm(j)

  kxx0 = wt%kx * end_orb%vec(1)
  kxx = kxx0 + wt%phi_x

  select case (wt%type)
  case (hyper_y$, hyper_y_old$)
    tmj%c_x = cos(kxx)
    tmj%s_x = sin(kxx)

    if (abs(kxx0) < 1e-4) then
      tmj%s_over_kx = end_orb%vec(1) * ((1 - kxx0**2/6) * cos(wt%phi_x) + (-kxx0/2 + kxx0**3/24) * sin(wt%phi_x))
    else
      tmj%s_over_kx = (tmj%s_x - sin(wt%phi_x)) / wt%kx
    endif

  !
  case (hyper_xy$, hyper_x$, hyper_xy_old$, hyper_x_old$)
    if (abs(kxx) > 30 .or. abs(wt%phi_x) > 30) then
      call err_set (err, x_plane$)
      return
    endif

    tmj%c_x = cosh(kxx)
    tmj%s_x = sinh(kxx)

    if (abs(kxx0) < 1e-4) then
      tmj%s_over_kx = end_orb%vec(1) * ((1 + kxx0**2/6) * cosh(wt%phi_x) + &
                                        kxx0 * (0.5_rp + kxx0**2/24) * sinh(wt%phi_x))
    else
      tmj%s_over_kx = (tmj%s_x - sinh(wt%phi_x)) / wt%kx
    endif
  end select

enddo

end subroutine update_wig_x_terms

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine update_wig_y_terms (err)

type (wiggler_computations_struct), pointer :: tmj

real(rp) kyy, kyy0
integer j
logical err

!

do j = 1, num_wig_terms
  wt => wig_term(j)
  tmj => tm(j)

  kyy0 = wt%ky * end_orb%vec(3)
  kyy = kyy0 + wt%phi_y

  select case (wt%type)
  case (hyper_y$, hyper_xy$, hyper_y_old$, hyper_xy_old$)
    if (abs(kyy) > 30 .or. abs(wt%phi_y) > 30) then
      call err_set (err, y_plane$)
      return
    endif

    tmj%c_y = cosh(kyy)
    tmj%s_y = sinh(kyy)

    if (abs(kyy0) < 1e-4) then
      tmj%one_minus_c_over_ky = end_orb%vec(3) * ((kyy0/2 + kyy0**3/24) * cosh(wt%phi_y) + (1 + kyy0**2/6) * sinh(wt%phi_y))
    else
      tmj%one_minus_c_over_ky = (tmj%c_y - cosh(wt%phi_y)) / wt%ky
    endif

  !
  case (hyper_x$, hyper_x_old$)
    tmj%c_y = cos(kyy)
    tmj%s_y = sin(kyy)

    if (abs(kyy0) < 1e-4) then
      tmj%one_minus_c_over_ky = end_orb%vec(3) * ((kyy0/2 - kyy0**3/24) * cos(wt%phi_y) + (1 - kyy0**2/6) * sin(wt%phi_y))
    else
      tmj%one_minus_c_over_ky = (cos(wt%phi_y) - tmj%c_y) / wt%ky
    endif
  end select
enddo

end subroutine update_wig_y_terms

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine update_wig_s_terms

real(rp) spz_offset
real(rp) kzz(1:num_wig_terms)

!

spz_offset = s + z_offset
if (orientation == -1) spz_offset = ele%value(l$) - spz_offset
kzz(1:num_wig_terms) = wig_term(1:num_wig_terms)%kz * spz_offset + wig_term(1:num_wig_terms)%phi_z

tm(1:num_wig_terms)%c_z = cos(kzz(1:num_wig_terms))
tm(1:num_wig_terms)%s_z = sin(kzz(1:num_wig_terms))

end subroutine update_wig_s_terms

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function Ax() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%Ax%coef * tm(j)%c_x * tm(j)%c_y * tm(j)%s_z
enddo

end function Ax

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function Ay() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%Ay%coef * tm(j)%s_x * tm(j)%s_y * tm(j)%s_z
enddo

end function Ay

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dAx_dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%Ax%coef * tm(j)%s_x * tm(j)%c_y * tm(j)%s_z
enddo

end function dAx_dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dAy_dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%Ay%coef * tm(j)%c_x * tm(j)%s_y * tm(j)%s_z
enddo

end function dAy_dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dAx_dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%Ax%coef * tm(j)%c_x * tm(j)%s_y * tm(j)%s_z
enddo

end function dAx_dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dAy_dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%Ay%coef * tm(j)%s_x * tm(j)%c_y * tm(j)%s_z
enddo

end function dAy_dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dint_Ax_dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%dint_Ax_dy%coef * tm(j)%s_over_kx * tm(j)%s_y * tm(j)%s_z
enddo

end function dint_Ax_dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dint_Ay_dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%dint_Ay_dx%coef * tm(j)%c_x * tm(j)%one_minus_c_over_ky * tm(j)%s_z
enddo

end function dint_Ay_dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dAz_dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%dAz_dx%coef * tm(j)%c_x * tm(j)%c_y * tm(j)%c_z
enddo

end function dAz_dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dAz_dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%dAz_dy%coef * tm(j)%s_x * tm(j)%s_y * tm(j)%c_z
enddo

end function dAz_dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dint_Ax_dy__dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%dint_Ax_dy%dx_coef * tm(j)%c_x * tm(j)%s_y * tm(j)%s_z
enddo

end function dint_Ax_dy__dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dint_Ay_dx__dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%dint_Ay_dx%dx_coef * tm(j)%s_x * tm(j)%one_minus_c_over_ky * tm(j)%s_z
enddo

end function dint_Ay_dx__dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dint_Ax_dy__dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%dint_Ax_dy%dy_coef * tm(j)%s_over_kx * tm(j)%c_y * tm(j)%s_z
enddo

end function dint_Ax_dy__dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dint_Ay_dx__dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%dint_Ay_dx%dy_coef * tm(j)%c_x * tm(j)%s_y * tm(j)%s_z
enddo

end function dint_Ay_dx__dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function Ax__dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%Ax%dx_coef * tm(j)%s_x * tm(j)%c_y * tm(j)%s_z
enddo

end function Ax__dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function Ay__dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%Ay%dx_coef * tm(j)%c_x * tm(j)%s_y * tm(j)%s_z
enddo

end function Ay__dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function Ax__dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%Ax%dy_coef * tm(j)%c_x * tm(j)%s_y * tm(j)%s_z
enddo

end function Ax__dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function Ay__dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%Ay%dy_coef * tm(j)%s_x * tm(j)%c_y * tm(j)%s_z
enddo

end function Ay__dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dAz_dx__dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%dAz_dx%dx_coef * tm(j)%s_x * tm(j)%c_y * tm(j)%c_z
enddo

end function dAz_dx__dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dAz_dx__dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%dAz_dx%dy_coef * tm(j)%c_x * tm(j)%s_y * tm(j)%c_z
enddo

end function dAz_dx__dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dAz_dy__dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%dAz_dy%dx_coef * tm(j)%c_x * tm(j)%s_y * tm(j)%c_z
enddo

end function dAz_dy__dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dAz_dy__dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%dAz_dy%dy_coef * tm(j)%s_x * tm(j)%c_y * tm(j)%c_z
enddo

end function dAz_dy__dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine radiation_kick()

type (ele_struct) :: temp_ele
type (ele_struct), pointer :: ele0

! Test if kick should be applied

if (.not. bmad_com%radiation_damping_on .and. .not. bmad_com%radiation_fluctuations_on) return

! g2 and g3 radiation integrals can be computed from the change in momentum.

g2 = dpx**2 + dpy**2
g3 = g2 * sqrt(g2)

! synch_rad_com%scale is normally 1 but can be set by a program for testing purposes.

call ran_gauss (this_ran)
dE_p = (1 + end_orb%vec(6)) * (fact_d * g2 + fact_f * sqrt(g3) * this_ran) * synch_rad_com%scale 

! And kick the particle.

end_orb%vec(2) = end_orb%vec(2) - (end_orb%vec(2) - Ax_saved) * dE_p
end_orb%vec(4) = end_orb%vec(4) - (end_orb%vec(4) - Ay_saved) * dE_p
end_orb%vec(6) = end_orb%vec(6) - dE_p * (1 + end_orb%vec(6))

! synch_ran_com%i_calc_on is, by default, False but a program can set this to True for testing purposes.

if (synch_rad_com%i_calc_on) then
  synch_rad_com%i2 = synch_rad_com%i2 + g2 * ds
  synch_rad_com%i3 = synch_rad_com%i3 + g3 * ds
  if (associated(ele%branch)) then
    temp_ele%mat6 = mat6
    temp_ele%vec0(1:5) = end_orb%vec(1:5) - matmul (mat6(1:5,1:6), start_orb%vec)
    temp_ele%vec0(6) = 0
    temp_ele%map_ref_orb_in = start_orb
    temp_ele%map_ref_orb_out = end_orb
    ele0 => ele%branch%ele(ele%ix_ele-1)
    call twiss_propagate1 (ele0, temp_ele)
    synch_rad_com%i5a = synch_rad_com%i5a + g3 * ds * (temp_ele%a%gamma * temp_ele%a%eta**2 + &
          2 * temp_ele%a%alpha * temp_ele%a%eta * temp_ele%a%etap + temp_ele%a%beta * temp_ele%a%etap**2)
    synch_rad_com%i5b = synch_rad_com%i5b + g3 * ds * (temp_ele%b%gamma * temp_ele%b%eta**2 + &
          2 * temp_ele%b%alpha * temp_ele%b%eta * temp_ele%b%etap + temp_ele%b%beta * temp_ele%b%etap**2)
  endif
endif

end subroutine radiation_kick

end subroutine

end module
