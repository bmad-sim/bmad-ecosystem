#include "CESR_platform.inc"

module symp_lie_mod

  use bmad_struct
  use bmad_interface
  use make_mat6_mod
  use em_field_mod   

contains

!+
! Subroutine symp_lie_bmad (ele, param, start, end, calc_mat6, track)
!
! Subroutine to track through an element (which gives the 0th order 
! taylor series) and optionally make the 6x6 transfer matrix (1st order 
! taylor series) as well.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele       -- Ele_struct: Element with transfer matrix
!   param     -- Param_struct: Parameters are needed for some elements.
!   start     -- Coord_struct: Coordinates at the beginning of element. 
!   calc_mat6 -- Logical: If True then make the 6x6 transfer matrix.
!   track      -- Track_struct: Structure holding the track information.
!     %save_track -- Logical: Set True if track is to be saved.
!
! Output:
!   ele    -- Ele_struct: Element with transfer matrix.
!     %mat6  -- 6x6 transfer matrix.
!   end    -- Coord_struct: Coordinates at the end of element.
!   track      -- Track_struct: Structure holding the track information.
!-

subroutine symp_lie_bmad (ele, param, start, end, calc_mat6, track)

  implicit none

  type save_coef_struct
    real(rp) coef, dx_coef, dy_coef
  end type

  type save_computations_struct
    type (save_coef_struct) a_y, dint_a_y_dx, da_z_dx, da_z_dy
    real(rp) c_x, s_x, c_y, s_y, c_z, s_z, s_x_kx, s_y_ky, c1_ky2
  end type

  type (save_computations_struct), allocatable, save :: tm(:)

  type (ele_struct), target :: ele
  type (coord_struct) :: start, end, start0
  type (param_struct)  param
  type (wig_term_struct), pointer :: wt
  type (track_struct) track

  real(rp) rel_E, rel_E2, rel_E3, ds, ds2, s, m6(6,6)
  real(rp) g_x, g_y, k1_norm, k1_skew, x_q, y_q, ks_tot_2
  real(rp), pointer :: mat6(:,:)
  real(rp), parameter :: z0 = 0, z1 = 1

  integer i, j, n_step

  logical calc_mat6

! If z_patch has not been calculated we need to track the reference orbit.

  if (ele%key == wiggler$ .and. any(start%vec /= 0) .and. &
                                      ele%value(z_patch$) == 0) then
    start0%vec = 0
    call track_it (start0, .false., .false.)
  endif

  call track_it (start, .true., calc_mat6)

!----------------------------------------------------------------------------
contains

subroutine track_it (start, real_track, local_calc_mat6)

  type (coord_struct) start
  logical real_track, local_calc_mat6, err

! init

  if (local_calc_mat6) mat6 => ele%mat6

  rel_E = (1 + start%vec(6))
  rel_E2 = rel_E**2
  rel_E3 = rel_E**3

  end = start
  err = .false.

! element offset 

  if (local_calc_mat6) then
    call drift_mat6_calc (mat6, ele%value(s_offset_tot$), end%vec)
  endif

  if (real_track) call offset_particle (ele, param, end, set$, set_canonical = .false.)

! init

  call compute_even_steps (ele%value(ds_step$), ele%value(l$), &
                                      bmad_com%default_ds_step, ds, n_step)
  ds2 = ds / 2

  s = 0   ! longitudianl position

  if (track%save_track) then
    call allocate_saved_orbit (track, n_step)
    track%pt(0)%s = 0
    track%pt(0)%orb = start
    track%n_pt = n_step
    if (local_calc_mat6) track%pt(0)%mat6 = mat6
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

    call update_wig_coefs (local_calc_mat6)
    call update_wig_y_terms (err); if (err) return

! loop over all steps

    do i = 1, n_step

! s half step

      s = s + ds2

! Drift_1 = P_x^2 / (2 * (1 + dE))

      call apply_p_x (local_calc_mat6)

! Drift_2 = (P_y - a_y)**2 / (2 * (1 + dE))

      call update_wig_x_s_terms (err); if (err) return
      call apply_wig_exp_int_ay (-1, local_calc_mat6)
      call apply_p_y (local_calc_mat6)
      call update_wig_y_terms (err); if (err) return
      call apply_wig_exp_int_ay (+1, local_calc_mat6)

! Kick = a_z

      end%vec(2) = end%vec(2) + ds * da_z_dx()
      end%vec(4) = end%vec(4) + ds * da_z_dy()

      if (local_calc_mat6) then
        mat6(2,1:6) = mat6(2,1:6) + ds * da_z_dx__dx() * mat6(1,1:6) + ds * da_z_dx__dy() * mat6(3,1:6)
        mat6(4,1:6) = mat6(4,1:6) + ds * da_z_dy__dx() * mat6(1,1:6) + ds * da_z_dy__dy() * mat6(3,1:6)
      endif 

! Drift_2

      call apply_wig_exp_int_ay (-1, local_calc_mat6)
      call apply_p_y (local_calc_mat6)
      call update_wig_y_terms (err); if (err) return
      call apply_wig_exp_int_ay (+1, local_calc_mat6)

! Drift_1

      call apply_p_x (local_calc_mat6)

! s half step

      s = s + ds2

      if (track%save_track) then
        track%pt(i)%s = s
        track%pt(i)%orb = end
        if (local_calc_mat6) track%pt(i)%mat6 = mat6
      endif

    enddo

! z_patch

  if (all(start%vec == 0)) then
    ele%value(z_patch$) = end%vec(5)
    end%vec(5) = 0
  else
    end%vec(5) = end%vec(5) - ele%value(z_patch$)
  endif

!----------------------------------------------------------------------------
! bend_sol_quad

  case (bend_sol_quad$)

    g_x = ele%value(g$) * cos (ele%value(bend_tilt$))
    g_y = ele%value(g$) * sin (ele%value(bend_tilt$))
    k1_norm = ele%value(k1$) * cos (2 * ele%value(quad_tilt$))
    k1_skew = ele%value(k1$) * sin (2 * ele%value(quad_tilt$))
    x_q = ele%value(x_quad$)
    y_q = ele%value(y_quad$)

! loop over all steps

    do i = 1, n_step

      s = s + ds2
      ks_tot_2 = (ele%value(ks$) + ele%value(dks_ds$) * s) / 2

      call bsq_drift1 (local_calc_mat6)
      call bsq_drift2 (local_calc_mat6)
      call bsq_kick (local_calc_mat6)
      call bsq_drift2 (local_calc_mat6)
      call bsq_drift1 (local_calc_mat6)

      s = s + ds2
      ks_tot_2 = (ele%value(ks$) + ele%value(dks_ds$) * s) / 2

      if (track%save_track) then
        track%pt(i)%s = s
        track%pt(i)%orb = end
        if (local_calc_mat6) track%pt(i)%mat6 = mat6
      endif

    enddo

!----------------------------------------------------------------------------
! unknown element

  case default

    print *, 'ERROR IN SYMP_LIE_BMAD: NOT YET IMPLEMENTED:', ele%key
    print *, '      FOR ELEMENT: ', ele%name
    call err_exit

  end select

! element offset

  if (local_calc_mat6) then
    call drift_mat6_calc (m6, -ele%value(s_offset_tot$), end%vec)
    mat6(1,1:6) = mat6(1,1:6) + m6(1,2) * mat6(2,1:6) + m6(1,6) * mat6(6,1:6)
    mat6(3,1:6) = mat6(3,1:6) + m6(3,4) * mat6(4,1:6) + m6(3,6) * mat6(6,1:6)
    mat6(5,1:6) = mat6(5,1:6) + m6(5,2) * mat6(2,1:6) + m6(5,4) * mat6(4,1:6) + m6(5,6) * mat6(6,1:6)

    if (ele%value(tilt_tot$) /= 0) call tilt_mat6 (mat6, ele%value(tilt_tot$))
  endif

  if (real_track) call offset_particle (ele, param, end, unset$, set_canonical = .false.)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine err_set (err)

  logical err

!

  print *, 'ERROR IN SYMP_LIE_BMAD: FLOATING OVERFLOW IN WIGGLER TRACKING.'
  print *, '      PARTICLE WILL BE TAGGED AS LOST.'
  param%lost = .true.
  end%vec(1) = 2 * bmad_com%max_aperture_limit
  end%vec(3) = 2 * bmad_com%max_aperture_limit
  err = .true.

end subroutine

!----------------------------------------------------------------------------
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
        call err_set (err)
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
        call err_set (err)
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

end subroutine



end module
