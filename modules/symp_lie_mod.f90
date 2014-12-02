module symp_lie_mod

use bmad_struct
use bmad_interface
use make_mat6_mod
use em_field_mod   
use random_mod
use track1_mod

type save_coef_struct
  real(rp) coef, dx_coef, dy_coef
end type

type wiggler_computations_struct
  type (save_coef_struct) a_y, dint_a_y_dx, da_z_dx, da_z_dy
  real(rp) c_x, s_x, c_y, s_y, c_z, s_z, s_x_kx, s_y_ky, c1_ky2
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
! Convention: Start and end p_y and p_x coordinates are the field free momentum.
! That is, at the start the coordinates are transformed by:
!   (p_x, p_y) -> (p_x + A_x, p_y + A_y)
! and at the end there is a transformation:
!   (p_x, p_y) -> (p_x - A_x, p_y - A_y)
! Where (A_x, A_y) components of the magnetic vector potential.
! If the start and end coordinates are in field free regions then (A_x, A_y) will be zero
! and the transformations will not affect the result. 
! The reason for this convention is to be able to compute the local bending radius via 
! tracking. Also this convention gives more "intuative" results when, say, using
! a single wiggler term as a "toy" model for a wiggler.
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

type (wiggler_computations_struct), allocatable :: tm(:)
type (wig_term_struct), pointer :: wt

real(rp) rel_E, rel_E2, rel_E3, ds, ds2, s, m6(6,6), kmat6(6,6)
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

if (present(track)) then
  call init_saved_orbit (track, n_step)
  call save_this_track_pt (0.0_rp)
endif

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

  call update_wig_coefs (calculate_mat6)
  call update_wig_y_terms (err); if (err) return
  call update_wig_x_s_terms (err); if (err) return

  ! Correction for finite vector potential at the entrance end

  call apply_wig_exp_int_ay(1, calculate_mat6)

  ! loop over all steps

  do i = 1, n_step

    ! s half step

    s = s + ds2

    ! Drift_1 = P_x^2 / (2 * (1 + dE))
    ! Note: We are using the gauge where A_x = 0.

    call apply_p_x (calculate_mat6)
    call update_wig_x_s_terms (err); if (err) return

    ! Drift_2 = (P_y - a_y)^2 / (2 * (1 + dE))

    call apply_wig_exp_int_ay (-1, calculate_mat6)

    call apply_p_y (calculate_mat6)
    call update_wig_y_terms (err); if (err) return

    call apply_wig_exp_int_ay (+1, calculate_mat6)

    ! Kick = a_z

    dpx = da_z_dx()
    dpy = da_z_dy()
    end_orb%vec(2) = end_orb%vec(2) + ds * dpx
    end_orb%vec(4) = end_orb%vec(4) + ds * dpy

    call radiation_kick()

    if (calculate_mat6) then
      mat6(2,1:6) = mat6(2,1:6) + ds * da_z_dx__dx() * mat6(1,1:6) + ds * da_z_dx__dy() * mat6(3,1:6)
      mat6(4,1:6) = mat6(4,1:6) + ds * da_z_dy__dx() * mat6(1,1:6) + ds * da_z_dy__dy() * mat6(3,1:6)
    endif 

    ! Drift_2

    call apply_wig_exp_int_ay (-1, calculate_mat6)

    call apply_p_y (calculate_mat6)
    call update_wig_y_terms (err); if (err) return

    call apply_wig_exp_int_ay (+1, calculate_mat6)

    ! Drift_1

    call apply_p_x (calculate_mat6)

    ! s half step

    s = s + ds2

    if (present(track)) call save_this_track_pt (s)

  enddo

  ! Correction for finite vector potential at exit end

  call update_wig_coefs (calculate_mat6)
  call update_wig_y_terms (err); if (err) return
  call update_wig_x_s_terms (err); if (err) return

  call apply_wig_exp_int_ay(-1, calculate_mat6)

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
  mat6(1,1:6) = mat6(1,1:6) + (ds2 / rel_E)           * mat6(2,1:6) - (ds2*end_orb%vec(2)/rel_E2)    * mat6(6,1:6) 
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
  mat6(3,1:6) = mat6(3,1:6) + (ds2 / rel_E)           * mat6(4,1:6) - (ds2*end_orb%vec(4)/rel_E2)    * mat6(6,1:6) 
  mat6(5,1:6) = mat6(5,1:6) - (ds2*end_orb%vec(4)/rel_E2) * mat6(4,1:6) + (ds2*end_orb%vec(4)**2/rel_E3) * mat6(6,1:6)
endif      

end subroutine apply_p_y

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine rf_drift1 (do_mat6)

logical do_mat6

! Drift_1 = (P_x - a_x)**2 / (2 * (1 + dE))

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

! Drift_1 = (P_x - a_x)**2 / (2 * (1 + dE))

end_orb%vec(2) = end_orb%vec(2) + end_orb%vec(3) * ks_tot_2   !  vec(2) - a_x
end_orb%vec(4) = end_orb%vec(4) + end_orb%vec(1) * ks_tot_2   !  vec(4) - dint_a_x_dy

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) + ks_tot_2 * mat6(3,1:6)
  mat6(4,1:6) = mat6(4,1:6) + ks_tot_2 * mat6(1,1:6)
endif      

!

call apply_p_x (do_mat6)

!

end_orb%vec(2) = end_orb%vec(2) - end_orb%vec(3) * ks_tot_2   !  vec(2) + a_x
end_orb%vec(4) = end_orb%vec(4) - end_orb%vec(1) * ks_tot_2   !  vec(4) + dint_a_x_dy

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

! Drift_2 = (P_y - a_y)**2 / (2 * (1 + dE))

end_orb%vec(2) = end_orb%vec(2) - end_orb%vec(3) * ks_tot_2   !  vec(2) - dint_a_y_dx
end_orb%vec(4) = end_orb%vec(4) - end_orb%vec(1) * ks_tot_2   !  vec(4) - a_y

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) - ks_tot_2 * mat6(3,1:6)
  mat6(4,1:6) = mat6(4,1:6) - ks_tot_2 * mat6(1,1:6)
endif      

!

call apply_p_y (do_mat6)

!

end_orb%vec(2) = end_orb%vec(2) + end_orb%vec(3) * ks_tot_2   !  vec(2) + dint_a_y_dx
end_orb%vec(4) = end_orb%vec(4) + end_orb%vec(1) * ks_tot_2   !  vec(4) + a_y

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
end_orb%vec(2) = end_orb%vec(2) + ds * dpx  ! da_z_dx
              
dpy = k1_norm * (end_orb%vec(3) - y_q) - k1_skew * end_orb%vec(1) - g_y
end_orb%vec(4) = end_orb%vec(4) + ds * dpy  ! da_z_dy

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) - ds * k1_norm * mat6(1,1:6) - ds * k1_skew * mat6(3,1:6)
  mat6(4,1:6) = mat6(4,1:6) - ds * k1_skew * mat6(1,1:6) + ds * k1_norm * mat6(3,1:6)
endif 

end subroutine bsq_kick

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine apply_wig_exp_int_ay (sgn, do_mat6)

integer sgn
logical do_mat6

end_orb%vec(2) = end_orb%vec(2) + sgn * dint_a_y_dx()
end_orb%vec(4) = end_orb%vec(4) + sgn * a_y()

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) + sgn * (dint_a_y_dx__dx() * mat6(1,1:6) + dint_a_y_dx__dy() * mat6(3,1:6))
  mat6(4,1:6) = mat6(4,1:6) + sgn * (a_y__dx()         * mat6(1,1:6) + a_y__dy()         * mat6(3,1:6))
endif      

end subroutine apply_wig_exp_int_ay

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine update_wig_coefs (do_mat6)

real(rp) factor, coef
integer j
logical do_mat6

factor = charge_dir * c_light / ele%value(p0c$)

do j = 1, num_wig_terms
  wt => wig_term(j)
  coef = factor * wt%coef * field_ele%value(polarity$)
  tm(j)%a_y%coef         = -coef * wt%kz      ! / (wt%kx * wt%ky)
  tm(j)%dint_a_y_dx%coef = -coef * wt%kz      ! / wt%ky**2
  tm(j)%da_z_dx%coef     = -coef 
  tm(j)%da_z_dy%coef     = -coef * wt%ky      ! / wt%kx
  if (wt%type == hyper_x$) then
    tm(j)%da_z_dy%coef     = -tm(j)%da_z_dy%coef
    tm(j)%dint_a_y_dx%coef = -tm(j)%dint_a_y_dx%coef 
  endif
enddo

if (.not. do_mat6) return

do j = 1, num_wig_terms
  wt => wig_term(j)
  tm(j)%a_y%dx_coef = tm(j)%a_y%coef
  tm(j)%a_y%dy_coef = tm(j)%a_y%coef
  tm(j)%dint_a_y_dx%dx_coef = tm(j)%dint_a_y_dx%coef * wt%kx
  tm(j)%dint_a_y_dx%dy_coef = tm(j)%dint_a_y_dx%coef
  tm(j)%da_z_dx%dx_coef = tm(j)%da_z_dx%coef * wt%kx
  tm(j)%da_z_dx%dy_coef = tm(j)%da_z_dx%coef * wt%ky
  tm(j)%da_z_dy%dx_coef = tm(j)%da_z_dy%coef
  tm(j)%da_z_dy%dy_coef = tm(j)%da_z_dy%coef * wt%ky
  
  if (wt%type == hyper_y$) then
    tm(j)%dint_a_y_dx%dx_coef = -tm(j)%dint_a_y_dx%dx_coef 
    tm(j)%da_z_dx%dx_coef     = -tm(j)%da_z_dx%dx_coef 
  elseif (wt%type == hyper_x$) then
    tm(j)%dint_a_y_dx%dy_coef = -tm(j)%dint_a_y_dx%dy_coef
    tm(j)%da_z_dx%dy_coef     = -tm(j)%da_z_dx%dy_coef      
  endif
enddo


end subroutine update_wig_coefs

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine update_wig_y_terms (err)

real(rp) kyy(1:num_wig_terms)
integer j
logical err

real(rp) exp_term
real(rp) exp_term_inv

real(rp) hypr_kyy(1:num_wig_terms)
real(rp) trig_kyy(1:num_wig_terms)

real(rp) exp_kyy(1:num_wig_terms)
real(rp) exp_kyy_inv(1:num_wig_terms)
real(rp) sin_kyy(1:num_wig_terms)
real(rp) cos_kyy(1:num_wig_terms)

integer n_trig, n_hypr

kyy(1:num_wig_terms) = wig_term(1:num_wig_terms)%ky * end_orb%vec(3)
n_hypr = 0
n_trig = 0
do j = 1, num_wig_terms
  wt => wig_term(j)
  if (abs(kyy(j)) < 1e-20) then
    !do nothing
  elseif (wt%type == hyper_y$ .or. wt%type == hyper_xy$) then
    if (abs(kyy(j)) > 30) then
      call err_set (err, y_plane$)
      return
    endif
    n_hypr = n_hypr + 1
    hypr_kyy(n_hypr) = kyy(j)    
  else
    n_trig = n_trig + 1
    trig_kyy(n_trig) = kyy(j)
  endif
enddo

exp_kyy(1:n_hypr) = EXP(hypr_kyy(1:n_hypr))
exp_kyy_inv(1:n_hypr) = 1.0/exp_kyy(1:n_hypr)

! !DIR$ vector always
! do j=1,n_trig
!   sin_kyy(j) = SIN(trig_kyy(j))
!   cos_kyy(j) = COS(trig_kyy(j))
! enddo
sin_kyy(1:n_trig) = SIN(trig_kyy(1:n_trig))
cos_kyy(1:n_trig) = COS(trig_kyy(1:n_trig))

n_hypr = 0
n_trig = 0
do j = 1, num_wig_terms
  wt => wig_term(j)
  if (abs(kyy(j)) < 1e-20) then
    tm(j)%c_y = 1
    tm(j)%s_y = kyy(j)
    tm(j)%s_y_ky = end_orb%vec(3)
    tm(j)%c1_ky2 = end_orb%vec(3)**2 / 2
    if (wt%type == hyper_x$) then
      tm(j)%c1_ky2 = -tm(j)%c1_ky2 
    endif
  elseif (wt%type == hyper_y$ .or. wt%type == hyper_xy$) then
    n_hypr = n_hypr + 1
    tm(j)%c_y = (exp_kyy(n_hypr)+exp_kyy_inv(n_hypr))*0.5d0 !  tm(j)%c_y = cosh(kyy)
    tm(j)%s_y = (exp_kyy(n_hypr)-exp_kyy_inv(n_hypr))*0.5d0 !  tm(j)%s_y = sinh(kyy)
    tm(j)%c1_ky2 = (tm(j)%c_y-1) / wt%ky**2   !  tm(j)%c1_ky2 = 2 * sinh(kyy(j)/2)**2 / wt%ky**2
    tm(j)%s_y_ky = tm(j)%s_y / wt%ky
  else
    n_trig = n_trig + 1
    tm(j)%c_y = cos_kyy(n_trig)
    tm(j)%s_y = sin_kyy(n_trig)
    tm(j)%s_y_ky = tm(j)%s_y / wt%ky
    tm(j)%c1_ky2 = -1*(1.0d0 - cos_kyy(n_trig)) / wt%ky**2
  endif
enddo

end subroutine update_wig_y_terms

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine update_wig_x_s_terms (err)

integer j
logical err

real(rp) kxx(1:num_wig_terms)
real(rp) kzz(1:num_wig_terms)

real(rp) hypr_kxx(1:num_wig_terms)
real(rp) trig_kxx(1:num_wig_terms)

real(rp) exp_kxx(1:num_wig_terms)
real(rp) exp_kxx_inv(1:num_wig_terms)
real(rp) cos_kxx(1:num_wig_terms)
real(rp) sin_kxx(1:num_wig_terms)

real(rp) exp_term
real(rp) exp_term_inv

real(rp) spz_offset

integer n_hypr, n_trig

kxx(1:num_wig_terms) = wig_term(1:num_wig_terms)%kx * end_orb%vec(1)
n_hypr = 0
n_trig = 0
do j = 1, num_wig_terms
  wt => wig_term(j)
  if (abs(kxx(j)) < 1e-20) then
    !do nothing
  elseif (wt%type == hyper_x$ .or. wt%type == hyper_xy$) then
    if (abs(kxx(j)) > 30) then
      call err_set (err, x_plane$)
      return
    endif
    n_hypr = n_hypr + 1
    hypr_kxx(n_hypr) = kxx(j)
  else
    n_trig = n_trig + 1
    trig_kxx(n_trig) = kxx(j)
  endif
enddo

exp_kxx(1:n_hypr) = EXP(hypr_kxx(1:n_hypr))
exp_kxx_inv(1:n_hypr) = 1.0/exp_kxx(1:n_hypr)

sin_kxx(1:n_trig) = sin(trig_kxx(1:n_trig))
cos_kxx(1:n_trig) = cos(trig_kxx(1:n_trig))

n_hypr = 0
n_trig = 0
do j = 1, num_wig_terms
  wt => wig_term(j)

  if (abs(kxx(j)) < 1e-20) then
    tm(j)%c_x = 1
    tm(j)%s_x = kxx(j)
    tm(j)%s_x_kx = end_orb%vec(1)
  elseif (wt%type == hyper_x$ .or. wt%type == hyper_xy$) then
    n_hypr = n_hypr + 1
    tm(j)%c_x = (exp_kxx(n_hypr)+exp_kxx_inv(n_hypr))*0.5d0 ! tm(j)%c_x = cosh(kxx(j))
    tm(j)%s_x = (exp_kxx(n_hypr)-exp_kxx_inv(n_hypr))*0.5d0 ! tm(j)%s_x = sinh(kxx(j))
    tm(j)%s_x_kx = tm(j)%s_x / wt%kx
  else
    n_trig = n_trig + 1
    tm(j)%c_x = cos_kxx(n_trig)
    tm(j)%s_x = sin_kxx(n_trig)
    tm(j)%s_x_kx = tm(j)%s_x / wt%kx
  endif
enddo

! update s-terms
spz_offset = s + z_offset
if (orientation == -1) spz_offset = ele%value(l$) - spz_offset
kzz(1:num_wig_terms) = wig_term(1:num_wig_terms)%kz * spz_offset + wig_term(1:num_wig_terms)%phi_z

! !DIR$ vector always
! do j=1,num_wig_terms
!   tm(j)%c_z = cos(kzz(j))
!   tm(j)%s_z = sin(kzz(j))
! enddo
tm(1:num_wig_terms)%c_z = cos(kzz(1:num_wig_terms))
tm(1:num_wig_terms)%s_z = sin(kzz(1:num_wig_terms))

end subroutine update_wig_x_s_terms

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function a_y() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%a_y%coef * tm(j)%s_x_kx * tm(j)%s_y_ky * tm(j)%s_z
enddo

end function a_y

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function da_y_dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%a_y%coef * tm(j)%c_x * tm(j)%s_y_ky * tm(j)%s_z
enddo

end function da_y_dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function da_y_dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%a_y%coef * tm(j)%s_x_kx * tm(j)%c_y * tm(j)%s_z
enddo

end function da_y_dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dint_a_y_dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%dint_a_y_dx%coef * tm(j)%c_x * tm(j)%c1_ky2 * tm(j)%s_z
enddo

end function dint_a_y_dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function da_z_dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%da_z_dx%coef * tm(j)%c_x * tm(j)%c_y * tm(j)%c_z
enddo

end function da_z_dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function da_z_dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%da_z_dy%coef * tm(j)%s_x_kx * tm(j)%s_y * tm(j)%c_z
enddo

end function da_z_dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dint_a_y_dx__dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%dint_a_y_dx%dx_coef * tm(j)%s_x * tm(j)%c1_ky2 * tm(j)%s_z
enddo

end function dint_a_y_dx__dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dint_a_y_dx__dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%dint_a_y_dx%dy_coef * tm(j)%c_x * tm(j)%s_y_ky * tm(j)%s_z
enddo

end function dint_a_y_dx__dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function a_y__dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%a_y%dx_coef * tm(j)%c_x * tm(j)%s_y_ky * tm(j)%s_z
enddo

end function a_y__dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function a_y__dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%a_y%dy_coef * tm(j)%s_x_kx * tm(j)%c_y * tm(j)%s_z
enddo

end function a_y__dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function da_z_dx__dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%da_z_dx%dx_coef * tm(j)%s_x * tm(j)%c_y * tm(j)%c_z
enddo

end function da_z_dx__dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function da_z_dx__dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%da_z_dx%dy_coef * tm(j)%c_x * tm(j)%s_y * tm(j)%c_z
enddo

end function da_z_dx__dy

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function da_z_dy__dx() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%da_z_dy%dx_coef * tm(j)%c_x * tm(j)%s_y * tm(j)%c_z
enddo

end function da_z_dy__dx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function da_z_dy__dy() result (value)

real(rp) value
integer j

!

value = 0
do j = 1, size(wig_term)
  value = value + tm(j)%da_z_dy%dy_coef * tm(j)%s_x_kx * tm(j)%c_y * tm(j)%c_z
enddo

end function da_z_dy__dy

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

end_orb%vec(2) = end_orb%vec(2) * (1 - dE_p)
end_orb%vec(4) = end_orb%vec(4) * (1 - dE_p)
end_orb%vec(6) = end_orb%vec(6)  - dE_p * (1 + end_orb%vec(6))

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
