#include "CESR_platform.inc"

module symp_lie_mod

use bmad_struct
use bmad_interface
use make_mat6_mod
use em_field_mod   
use random_mod

type save_coef_struct
  real(rp) coef, dx_coef, dy_coef
end type

type save_computations_struct
  type (save_coef_struct) a_y, dint_a_y_dx, da_z_dx, da_z_dy
  real(rp) c_x, s_x, c_y, s_y, c_z, s_z, s_x_kx, s_y_ky, c1_ky2
end type

private save_coef_struct, save_computations_struct

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine symp_lie_bmad (ele, param, start, end, calc_mat6, track, offset_ele)
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
!   start      -- Coord_struct: Coordinates at the beginning of element. 
!   calc_mat6  -- Logical: If True then make the 6x6 transfer matrix.
!   track      -- Track_struct, optional: Structure holding the track information.
!     %save_track -- Logical: Set True if track is to be saved.
!   offset_ele -- Logical, optional: Offset the element using ele%value(x_offset$), etc.
!                   Default is True.
!
! Output:
!   ele    -- Ele_struct: Element with transfer matrix.
!     %mat6(6,6) -- 6x6 transfer matrix.
!     %vec0(6)   -- 0th order part of the transfer matrix.
!   end    -- Coord_struct: Coordinates at the end of element.
!   track      -- Track_struct, optional: Structure holding the track information.
!-

subroutine symp_lie_bmad (ele, param, start, end, calc_mat6, track, offset_ele)

implicit none

type (ele_struct), target :: ele
type (coord_struct) :: start, end
type (lat_param_struct)  param
type (track_struct), optional :: track

type (save_computations_struct), allocatable, save :: tm(:)
type (wig_term_struct), pointer :: wt

real(rp) rel_E, rel_E2, rel_E3, ds, ds2, s, m6(6,6), x_pitch, y_pitch
real(rp) g_x, g_y, k1_norm, k1_skew, x_q, y_q, ks_tot_2, ks, dks_ds
real(rp), pointer :: mat6(:,:)
real(rp), parameter :: z0 = 0, z1 = 1
real(rp) gamma_0, fact_d, fact_f, this_ran, g2, g3
real(rp) dE_p, dpx, dpy
integer i, j, n_step

logical calc_mat6, calculate_mat6, err, save_track, do_offset
logical, optional :: offset_ele

character(16) :: r_name = 'symp_lie_bmad'

! init

calculate_mat6 = (calc_mat6 .or. synch_rad_com%i_calc_on)

do_offset = logic_option (.true., offset_ele)
rel_E = (1 + start%vec(6))
rel_E2 = rel_E**2
rel_E3 = rel_E**3

end = start
err = .false.

x_pitch = ele%value(x_pitch_tot$)
y_pitch = ele%value(y_pitch_tot$)

save_track = present(track)
if (save_track) save_track = track%save_track 

! element offset 

if (calculate_mat6) then
  mat6 => ele%mat6
  call drift_mat6_calc (mat6, ele%value(s_offset_tot$), end%vec)
endif

if (do_offset) call offset_particle (ele, param, end, set$, set_canonical = .false.)

! init

call compute_even_steps (ele%value(ds_step$), ele%value(l$), &
                                    bmad_com%default_ds_step, ds, n_step)
ds2 = ds / 2

s = 0   ! longitudianl position

if (save_track) then
  call allocate_saved_orbit (track, n_step)
  track%n_pt = n_step
  call save_this_track_pt (0, 0.0_rp)
endif

! radiation damping and fluctuations...
! The same kick is applied over the entire wiggler to save time.

if ((bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on)) then

  gamma_0 = ele%value(E_TOT$) / m_electron

  fact_d = 0
  if (bmad_com%radiation_damping_on) fact_d = 2 * r_e * gamma_0**3 * ds / 3

  fact_f = 0
  if (bmad_com%radiation_fluctuations_on) then
!    call ran_gauss (this_ran)
!    fact_f = sqrt(rad_fluct_const * ele%value(l$) * gamma_0**5) * this_ran / n_step
    fact_f = sqrt(rad_fluct_const * ds * gamma_0**5) 
  endif

endif

!------------------------------------------------------------------
! select the element

select case (ele%key)

!------------------------------------------------------------------
! Wiggler

Case (wiggler$)

  if (.not. allocated(tm)) then
    allocate (tm(size(ele%wig_term)))
  elseif (size(tm) < size(ele%wig_term)) then
    deallocate(tm)
    allocate (tm(size(ele%wig_term)))
  endif

  call update_wig_coefs (calculate_mat6)
  call update_wig_y_terms (err); if (err) return

  ! loop over all steps

  do i = 1, n_step

    ! s half step

    s = s + ds2

    ! Drift_1 = P_x^2 / (2 * (1 + dE))

    call apply_p_x (calculate_mat6)

    ! Drift_2 = (P_y - a_y)**2 / (2 * (1 + dE))

    call update_wig_x_s_terms (err); if (err) return
    call apply_wig_exp_int_ay (-1, calculate_mat6)
    call apply_p_y (calculate_mat6)
    call update_wig_y_terms (err); if (err) return
    call apply_wig_exp_int_ay (+1, calculate_mat6)

    ! Kick = a_z

    dpx = da_z_dx()
    dpy = da_z_dy()
    end%vec(2) = end%vec(2) + ds * dpx
    end%vec(4) = end%vec(4) + ds * dpy

    call radiation_kick()

    if (calculate_mat6) then
      mat6(2,1:6) = mat6(2,1:6) + ds * da_z_dx__dx() * mat6(1,1:6) + &
                                                ds * da_z_dx__dy() * mat6(3,1:6)
      mat6(4,1:6) = mat6(4,1:6) + ds * da_z_dy__dx() * mat6(1,1:6) + &
                                                ds * da_z_dy__dy() * mat6(3,1:6)
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

    if (save_track) call save_this_track_pt (i, s)

  enddo

  ! z_patch: This should have been computed if doing tracking with an offset.

  if (ele%value(z_patch$) == 0 .and. do_offset) then
    call out_io (s_fatal$, r_name, 'WIGGLER Z_PATCH VALUE HAS NOT BEEN COMPUTED!')
    call err_exit 
  endif

  end%vec(5) = end%vec(5) - ele%value(z_patch$)

!----------------------------------------------------------------------------
! bend_sol_quad

case (bend_sol_quad$, solenoid$, quadrupole$)

  if (ele%key == bend_sol_quad$) then
    g_x = ele%value(g$) * cos (ele%value(bend_tilt$))
    g_y = ele%value(g$) * sin (ele%value(bend_tilt$))
    k1_norm = ele%value(k1$) * cos (2 * ele%value(quad_tilt$))
    k1_skew = ele%value(k1$) * sin (2 * ele%value(quad_tilt$))
    x_q = ele%value(x_quad$)
    y_q = ele%value(y_quad$)
    ks = ele%value(ks$)
    dks_ds = ele%value(dks_ds$)
  elseif (ele%key == solenoid$) then
    g_x = 0
    g_y = 0
    k1_norm = 0
    k1_skew = 0
    x_q = 0
    y_q = 0
    ks = ele%value(ks$)
    dks_ds = 0
  elseif (ele%key == quadrupole$) then
    g_x = 0
    g_y = 0
    k1_norm = ele%value(k1$) 
    k1_skew = 0
    x_q = 0
    y_q = 0
    ks = 0
    dks_ds = 0
  endif

! loop over all steps

  do i = 1, n_step

    s = s + ds2
    ks_tot_2 = (ks + dks_ds * s) / 2

    call bsq_drift1 (calculate_mat6)
    call bsq_drift2 (calculate_mat6)
    call bsq_kick (calculate_mat6)
    call bsq_drift2 (calculate_mat6)
    call bsq_drift1 (calculate_mat6)

    s = s + ds2
    ks_tot_2 = (ks + dks_ds * s) / 2

    if (save_track) call save_this_track_pt (i, s)

  enddo

!----------------------------------------------------------------------------
! unknown element

case default

  print *, 'ERROR IN SYMP_LIE_BMAD: NOT YET IMPLEMENTED:', ele%key
  print *, '      FOR ELEMENT: ', ele%name
  call err_exit

end select

! element offset

if (calculate_mat6) then
  call drift_mat6_calc (m6, -ele%value(s_offset_tot$), end%vec)
  mat6(1,1:6) = mat6(1,1:6) + m6(1,2) * mat6(2,1:6) + m6(1,6) * mat6(6,1:6)
  mat6(3,1:6) = mat6(3,1:6) + m6(3,4) * mat6(4,1:6) + m6(3,6) * mat6(6,1:6)
  mat6(5,1:6) = mat6(5,1:6) + m6(5,2) * mat6(2,1:6) + m6(5,4) * mat6(4,1:6) + &
                                                                 m6(5,6) * mat6(6,1:6)

  if (ele%value(tilt_tot$) /= 0) call tilt_mat6 (mat6, ele%value(tilt_tot$))
  if (x_pitch /= 0 .or. y_pitch /= 0) call mat6_add_pitch (ele, mat6)
endif

if (do_offset) call offset_particle (ele, param, end, unset$, set_canonical = .false.)

! Correct for finite pitches & calc vec0

if (calculate_mat6) then
  ele%vec0(1:5) = end%vec(1:5) - matmul (mat6(1:5,1:6), start%vec)
  ele%vec0(6) = 0
endif

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
contains

subroutine err_set (err, plane)

logical err
integer plane

!

print *, 'ERROR IN SYMP_LIE_BMAD: FLOATING OVERFLOW IN WIGGLER TRACKING.'
print *, '      PARTICLE WILL BE TAGGED AS LOST.'
param%plane_lost_at = plane
param%lost = .true.
end%vec(1) = 2 * bmad_com%max_aperture_limit
end%vec(3) = 2 * bmad_com%max_aperture_limit
err = .true.

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine save_this_track_pt (ix, s)

real(rp) s
integer ix

!

track%pt(ix)%s = s
track%pt(ix)%orb = end
call offset_particle (ele, param, track%pt(ix)%orb, unset$, set_canonical = .false.)
  
if (calculate_mat6) track%pt(ix)%mat6 = mat6

if (ele%value(tilt_tot$) /= 0) call tilt_mat6 (track%pt(ix)%mat6, ele%value(tilt_tot$))
if (x_pitch /= 0 .or. y_pitch /= 0) call mat6_add_pitch (ele, track%pt(ix)%mat6)

if (calculate_mat6) then
  track%pt(ix)%vec0(1:5) = track%pt(ix)%orb%vec(1:5) - matmul (mat6(1:5,1:6), start%vec)
  track%pt(ix)%vec0(6) = 0
endif
 
end subroutine  

!----------------------------------------------------------------------------
! contains

subroutine apply_p_x (do_mat6)

logical do_mat6

end%vec(1) = end%vec(1) + ds2 * end%vec(2) / rel_E
end%vec(5) = end%vec(5) - ds2 * end%vec(2)**2 / (2*rel_E2)

if (do_mat6) then
  mat6(1,1:6) = mat6(1,1:6) + (ds2 / rel_E)           * mat6(2,1:6) - (ds2*end%vec(2)/rel_E2)    * mat6(6,1:6) 
  mat6(5,1:6) = mat6(5,1:6) - (ds2*end%vec(2)/rel_E2) * mat6(2,1:6) + (ds2*end%vec(2)**2/rel_E3) * mat6(6,1:6)
endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine apply_p_y (do_mat6)

logical do_mat6

end%vec(3) = end%vec(3) + ds2 * end%vec(4) / rel_E
end%vec(5) = end%vec(5) - ds2 * end%vec(4)**2 / (2*rel_E2)

if (do_mat6) then
  mat6(3,1:6) = mat6(3,1:6) + (ds2 / rel_E)           * mat6(4,1:6) - (ds2*end%vec(4)/rel_E2)    * mat6(6,1:6) 
  mat6(5,1:6) = mat6(5,1:6) - (ds2*end%vec(4)/rel_E2) * mat6(4,1:6) + (ds2*end%vec(4)**2/rel_E3) * mat6(6,1:6)
endif      

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine bsq_drift1 (do_mat6)

logical do_mat6

! Drift_1 = (P_x - a_x)**2 / (2 * (1 + dE))

end%vec(2) = end%vec(2) + end%vec(3) * ks_tot_2   !  vec(2) - a_x
end%vec(4) = end%vec(4) + end%vec(1) * ks_tot_2   !  vec(4) - dint_a_x_dy

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) + ks_tot_2 * mat6(3,1:6)
  mat6(4,1:6) = mat6(4,1:6) + ks_tot_2 * mat6(1,1:6)
endif      

!

call apply_p_x (do_mat6)

!

end%vec(2) = end%vec(2) - end%vec(3) * ks_tot_2   !  vec(2) + a_x
end%vec(4) = end%vec(4) - end%vec(1) * ks_tot_2   !  vec(4) + dint_a_x_dy

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) - ks_tot_2 * mat6(3,1:6)
  mat6(4,1:6) = mat6(4,1:6) - ks_tot_2 * mat6(1,1:6)
endif  

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine bsq_drift2 (do_mat6)

logical do_mat6

! Drift_2 = (P_y - a_y)**2 / (2 * (1 + dE))

end%vec(2) = end%vec(2) - end%vec(3) * ks_tot_2   !  vec(2) - dint_a_y_dx
end%vec(4) = end%vec(4) - end%vec(1) * ks_tot_2   !  vec(4) - a_y

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) - ks_tot_2 * mat6(3,1:6)
  mat6(4,1:6) = mat6(4,1:6) - ks_tot_2 * mat6(1,1:6)
endif      

!

call apply_p_y (do_mat6)

!

end%vec(2) = end%vec(2) + end%vec(3) * ks_tot_2   !  vec(2) + dint_a_y_dx
end%vec(4) = end%vec(4) + end%vec(1) * ks_tot_2   !  vec(4) + a_y

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) + ks_tot_2 * mat6(3,1:6)
  mat6(4,1:6) = mat6(4,1:6) + ks_tot_2 * mat6(1,1:6)
endif  

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine bsq_kick (do_mat6)

logical do_mat6

end%vec(2) = end%vec(2) + ds * &  ! da_z_dx
              (k1_norm * (x_q - end%vec(1)) - k1_skew * end%vec(3) - g_x)    
end%vec(4) = end%vec(4) + ds * &  ! da_z_dy
              (k1_norm * (end%vec(3) - y_q) - k1_skew * end%vec(1) - g_y)    

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) - ds * k1_norm * mat6(1,1:6) - ds * k1_skew * mat6(3,1:6)
  mat6(4,1:6) = mat6(4,1:6) - ds * k1_skew * mat6(1,1:6) + ds * k1_norm * mat6(3,1:6)
endif 

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine apply_wig_exp_int_ay (sgn, do_mat6)

integer sgn
logical do_mat6

end%vec(2) = end%vec(2) + sgn * dint_a_y_dx()
end%vec(4) = end%vec(4) + sgn * a_y()

if (do_mat6) then
  mat6(2,1:6) = mat6(2,1:6) + sgn * &
          (dint_a_y_dx__dx() * mat6(1,1:6) + dint_a_y_dx__dy() * mat6(3,1:6))
  mat6(4,1:6) = mat6(4,1:6) + sgn * &
          (a_y__dx()         * mat6(1,1:6) + a_y__dy()         * mat6(3,1:6))
endif      

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine update_wig_coefs (do_mat6)

real(rp) factor, coef
logical do_mat6

factor = c_light / ele%value(p0c$)

do j = 1, size(ele%wig_term)
  wt => ele%wig_term(j)
  coef = factor * wt%coef * ele%value(polarity$)
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

do j = 1, size(ele%wig_term)
  wt => ele%wig_term(j)
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


end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine update_wig_y_terms (err)

real(rp) kyy
logical err

do j = 1, size(ele%wig_term)
  wt => ele%wig_term(j)
  kyy = wt%ky * end%vec(3)
  if (abs(kyy) < 1e-20) then
    tm(j)%c_y = 1
    tm(j)%s_y = kyy
    tm(j)%s_y_ky = end%vec(3)
    tm(j)%c1_ky2 = end%vec(3)**2 / 2
    if (wt%type == hyper_x$) tm(j)%c1_ky2 = -tm(j)%c1_ky2 
  elseif (wt%type == hyper_y$ .or. wt%type == hyper_xy$) then
    if (abs(kyy) > 30) then
      call err_set (err, y_plane$)
      return
    endif
    tm(j)%c_y = cosh(kyy)
    tm(j)%s_y = sinh(kyy)
    tm(j)%s_y_ky = tm(j)%s_y / wt%ky
    tm(j)%c1_ky2 = 2 * sinh(kyy/2)**2 / wt%ky**2
  else
    tm(j)%c_y = cos(kyy)
    tm(j)%s_y = sin(kyy)
    tm(j)%s_y_ky = tm(j)%s_y / wt%ky
    tm(j)%c1_ky2 = -2 * sin(kyy/2)**2 / wt%ky**2
  endif
enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine update_wig_x_s_terms (err)

real(rp) kxx, kzz
logical err

!

do j = 1, size(ele%wig_term)
  wt => ele%wig_term(j)

  kxx = wt%kx * end%vec(1)
  if (abs(kxx) < 1e-20) then
    tm(j)%c_x = 1
    tm(j)%s_x = kxx
    tm(j)%s_x_kx = end%vec(1)
  elseif (wt%type == hyper_x$ .or. wt%type == hyper_xy$) then
    if (abs(kxx) > 30) then
      call err_set (err, x_plane$)
      return
    endif
    tm(j)%c_x = cosh(kxx)
    tm(j)%s_x = sinh(kxx)
    tm(j)%s_x_kx = tm(j)%s_x / wt%kx
  else
    tm(j)%c_x = cos(kxx)
    tm(j)%s_x = sin(kxx)
    tm(j)%s_x_kx = tm(j)%s_x / wt%kx
  endif

  kzz = wt%kz * s + wt%phi_z
  tm(j)%c_z = cos(kzz)
  tm(j)%s_z = sin(kzz)

enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function a_y() result (value)

real(rp) value

value = 0
do j = 1, size(ele%wig_term)
  value = value + tm(j)%a_y%coef * tm(j)%s_x_kx * tm(j)%s_y_ky * tm(j)%s_z
enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dint_a_y_dx() result (value)

real(rp) value

value = 0
do j = 1, size(ele%wig_term)
  value = value + tm(j)%dint_a_y_dx%coef * tm(j)%c_x * tm(j)%c1_ky2 * tm(j)%s_z
enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function da_z_dx() result (value)

real(rp) value

value = 0
do j = 1, size(ele%wig_term)
  value = value + tm(j)%da_z_dx%coef * tm(j)%c_x * tm(j)%c_y * tm(j)%c_z
enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function da_z_dy() result (value)

real(rp) value

value = 0
do j = 1, size(ele%wig_term)
  value = value + tm(j)%da_z_dy%coef * tm(j)%s_x_kx * tm(j)%s_y * tm(j)%c_z
enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dint_a_y_dx__dx() result (value)

real(rp) value

value = 0
do j = 1, size(ele%wig_term)
  value = value + tm(j)%dint_a_y_dx%dx_coef * tm(j)%s_x * tm(j)%c1_ky2 * tm(j)%s_z
enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function dint_a_y_dx__dy() result (value)

real(rp) value

value = 0
do j = 1, size(ele%wig_term)
  value = value + tm(j)%dint_a_y_dx%dy_coef * tm(j)%c_x * tm(j)%s_y_ky * tm(j)%s_z
enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function a_y__dx() result (value)

real(rp) value

value = 0
do j = 1, size(ele%wig_term)
  value = value + tm(j)%a_y%dx_coef * tm(j)%c_x * tm(j)%s_y_ky * tm(j)%s_z
enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function a_y__dy() result (value)

real(rp) value

value = 0
do j = 1, size(ele%wig_term)
  value = value + tm(j)%a_y%dy_coef * tm(j)%s_x_kx * tm(j)%c_y * tm(j)%s_z
enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function da_z_dx__dx() result (value)

real(rp) value

value = 0
do j = 1, size(ele%wig_term)
  value = value + tm(j)%da_z_dx%dx_coef * tm(j)%s_x * tm(j)%c_y * tm(j)%c_z
enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function da_z_dx__dy() result (value)

real(rp) value

value = 0
do j = 1, size(ele%wig_term)
  value = value + tm(j)%da_z_dx%dy_coef * tm(j)%c_x * tm(j)%s_y * tm(j)%c_z
enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function da_z_dy__dx() result (value)

real(rp) value

value = 0
do j = 1, size(ele%wig_term)
  value = value + tm(j)%da_z_dy%dx_coef * tm(j)%c_x * tm(j)%s_y * tm(j)%c_z
enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

function da_z_dy__dy() result (value)

real(rp) value

value = 0
do j = 1, size(ele%wig_term)
  value = value + tm(j)%da_z_dy%dy_coef * tm(j)%s_x_kx * tm(j)%c_y * tm(j)%c_z
enddo

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine radiation_kick()

type (ele_struct), save :: temp_ele

!

if ((bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on)) then

  g2 = dpx**2 + dpy**2
  g3 = g2 * sqrt(g2)

  !  dE_p = (1 + end%vec(6)) * (fact_d * g2 + fact_f * sqrt(g3)) * synch_rad_com%scale
  call ran_gauss (this_ran)
  dE_p = (1 + end%vec(6)) * (fact_d * g2 + fact_f * sqrt(g3) * this_ran) * &
                                                                      synch_rad_com%scale 

  end%vec(2) = end%vec(2) * (1 - dE_p)
  end%vec(4) = end%vec(4) * (1 - dE_p)
  end%vec(6) = end%vec(6)  - dE_p * (1 + end%vec(6))

  if (synch_rad_com%i_calc_on) then
    synch_rad_com%i2 = synch_rad_com%i2 + g2 * ds
    synch_rad_com%i3 = synch_rad_com%i3 + g3 * ds
    temp_ele%mat6 = mat6
    temp_ele%vec0(1:5) = end%vec(1:5) - matmul (mat6(1:5,1:6), start%vec)
    temp_ele%vec0(6) = 0
    temp_ele%map_ref_orb_in = start
    temp_ele%map_ref_orb_out = end
    call twiss_propagate1 (synch_rad_com%ele0, temp_ele)
    synch_rad_com%i5a = synch_rad_com%i5a + g3 * ds * &
                  (temp_ele%a%gamma * temp_ele%a%eta**2 + &
                  2 * temp_ele%a%alpha * temp_ele%a%eta * temp_ele%a%etap + &
                  temp_ele%a%beta * temp_ele%a%etap**2)
    synch_rad_com%i5b = synch_rad_com%i5b + g3 * ds * &
                  (temp_ele%b%gamma * temp_ele%b%eta**2 + &
                  2 * temp_ele%b%alpha * temp_ele%b%eta * temp_ele%b%etap + &
                  temp_ele%b%beta * temp_ele%b%etap**2)
  endif

endif


end subroutine

end subroutine

end module
