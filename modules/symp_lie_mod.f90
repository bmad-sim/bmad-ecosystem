#include "CESR_platform.inc"

module symp_lie_mod

  use bmad_struct
  use bmad_interface
  use make_mat6_mod
  use em_field_mod   ! For the track_com_struct

contains

!+
! Subroutine symp_lie_bmad (ele, param, start, end, calc_mat6)
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
!
! Output:
!   ele    -- Ele_struct: Element with transfer matrix.
!     %mat6  -- 6x6 transfer matrix.
!   end    -- Coord_struct: Coordinates at the end of element.
!
! To save the orbit through the ele there is a global variable:
!   track_com  -- Symp_lie_com_struct: Global variable
!     %save_track -- Logical: Set True to save the orbit
!     %s(:)       -- real(rp): S-positions
!     %orb(:)     -- Coord_struct: orbit.
!-

subroutine symp_lie_bmad (ele, param, start, end, calc_mat6)

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
  type (coord_struct) :: start, end, orb, start0
  type (param_struct)  param
  type (wig_term_struct), pointer :: wt

  real(rp) rel_E, rel_E2, rel_E3, ds, ds2, s, m6(6,6)
  real(rp), pointer :: mat6(:,:)
  real(rp), parameter :: z0 = 0, z1 = 1

  integer i, j, n

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

subroutine track_it (start, real_track, calc_mat6)

  type (coord_struct) start
  logical real_track, calc_mat6

! init

  if (calc_mat6) mat6 => ele%mat6

  rel_E = (1 + start%vec(6))
  rel_E2 = rel_E**2
  rel_E3 = rel_E**3

  end = start

! element offset 

  if (calc_mat6) then
    call drift_mat6_calc (mat6, ele%value(s_offset$), end%vec)
  endif

  if (real_track) call offset_particle (ele, param, end, set$, set_canonical = .false.)

!------------------------------------------------------------------
! select the element

  select case (ele%key)

! Wiggler

  Case (wiggler$)

    if (.not. allocated(tm) .or. size(tm) < size(ele%wig_term)) then
      if (allocated(tm)) deallocate(tm)
      allocate (tm(size(ele%wig_term)))
    endif

    ds = ele%value(l$) / ele%num_steps
    ds2 = ds / 2

    s = 0   ! longitudianl position

    call update_coefs
    call update_y_terms

    if (track_com%save_track) then
      n = ele%num_steps 
      if (associated(track_com%s)) then
        if (size(track_com%s) /= n + 1) then
          deallocate (track_com%s, track_com%orb)
          allocate (track_com%s(0:n), track_com%orb(0:n))
        endif
      else
        allocate (track_com%s(0:n), track_com%orb(0:n))
      endif
      track_com%s(0) = 0
      track_com%orb(0) = start
    endif

! loop over all steps

    do i = 1, ele%num_steps

! check for overflow

      if (maxval(abs(ele%wig_term(:)%kx * end%vec(1))) > 30 .or. &
                  maxval(abs(ele%wig_term(:)%ky * end%vec(3))) > 30) then
        print *, 'ERROR IN SYMP_LIE_BMAD: ', &
                                     'FLOATING OVERFLOW IN WIGGLER TRACKING.'
        print *, '      PARTICLE WILL BE TAGGED AS LOST.'
        param%lost = .true.
        return
      endif

! s half step

      s = s + ds2

! Drift_1 = P_x^2 / (2 * (1 + dE))

      end%vec(1) = end%vec(1) + ds2 * end%vec(2) / rel_E
      end%vec(5) = end%vec(5) - ds2 * end%vec(2)**2 / (2*rel_E2)

      if (calc_mat6) then
        mat6(1,1:6) = mat6(1,1:6) + (ds2 / rel_E)          * mat6(2,1:6) - (ds2*end%vec(2)/rel_E2)    * mat6(6,1:6) 
        mat6(5,1:6) = mat6(5,1:6) - (ds2*end%vec(2)/rel_E2) * mat6(2,1:6) + (ds2*end%vec(2)**2/rel_E3) * mat6(6,1:6)
      endif

! Drift_2 = (P_y - a_y)**2 / (2 * (1 + dE))

      call update_x_s_terms

      end%vec(2) = end%vec(2) - dint_a_y_dx()
      end%vec(4) = end%vec(4) - a_y()

      if (calc_mat6) then
        mat6(2,1:6) = mat6(2,1:6) - dint_a_y_dx__dx() * mat6(1,1:6) - dint_a_y_dx__dy() * mat6(3,1:6)
        mat6(4,1:6) = mat6(4,1:6) - a_y__dx()         * mat6(1,1:6) - a_y__dy()         * mat6(3,1:6)
      endif      

!

      end%vec(3) = end%vec(3) + ds2 * end%vec(4) / rel_E
      end%vec(5) = end%vec(5) - ds2 * end%vec(4)**2 / (2*rel_E2)

      if (calc_mat6) then
        mat6(3,1:6) = mat6(3,1:6) + (ds2 / rel_E)          * mat6(4,1:6) - (ds2*end%vec(4)/rel_E2)    * mat6(6,1:6) 
        mat6(5,1:6) = mat6(5,1:6) - (ds2*end%vec(4)/rel_E2) * mat6(4,1:6) + (ds2*end%vec(4)**2/rel_E3) * mat6(6,1:6)
      endif      

!

      call update_y_terms

      end%vec(2) = end%vec(2) + dint_a_y_dx()
      end%vec(4) = end%vec(4) + a_y()

      if (calc_mat6) then
        mat6(2,1:6) = mat6(2,1:6) + dint_a_y_dx__dx() * mat6(1,1:6) + dint_a_y_dx__dy() * mat6(3,1:6)
        mat6(4,1:6) = mat6(4,1:6) + a_y__dx()         * mat6(1,1:6) + a_y__dy()         * mat6(3,1:6)
      endif      

! Kick = a_z

      end%vec(2) = end%vec(2) + ds * da_z_dx()
      end%vec(4) = end%vec(4) + ds * da_z_dy()

      if (calc_mat6) then
        mat6(2,1:6) = mat6(2,1:6) + ds * da_z_dx__dx() * mat6(1,1:6) + ds * da_z_dx__dy() * mat6(3,1:6)
        mat6(4,1:6) = mat6(4,1:6) + ds * da_z_dy__dx() * mat6(1,1:6) + ds * da_z_dy__dy() * mat6(3,1:6)
      endif 

! Drift_2

      end%vec(2) = end%vec(2) - dint_a_y_dx()
      end%vec(4) = end%vec(4) - a_y()

      if (calc_mat6) then
        mat6(2,1:6) = mat6(2,1:6) - dint_a_y_dx__dx() * mat6(1,1:6) - dint_a_y_dx__dy() * mat6(3,1:6)
        mat6(4,1:6) = mat6(4,1:6) - a_y__dx()         * mat6(1,1:6) - a_y__dy()         * mat6(3,1:6)
      endif      

!

      end%vec(3) = end%vec(3) + ds2 * end%vec(4) / rel_E
      end%vec(5) = end%vec(5) - ds2 * end%vec(4)**2 / (2*rel_E2)

      if (calc_mat6) then
        mat6(3,1:6) = mat6(3,1:6) + (ds2 / rel_E)          * mat6(4,1:6) - (ds2*end%vec(4)/rel_E2)    * mat6(6,1:6) 
        mat6(5,1:6) = mat6(5,1:6) - (ds2*end%vec(4)/rel_E2) * mat6(4,1:6) + (ds2*end%vec(4)**2/rel_E3) * mat6(6,1:6)
      endif      

!

      call update_y_terms

      end%vec(2) = end%vec(2) + dint_a_y_dx()
      end%vec(4) = end%vec(4) + a_y()

      if (calc_mat6) then
        mat6(2,1:6) = mat6(2,1:6) + dint_a_y_dx__dx() * mat6(1,1:6) + dint_a_y_dx__dy() * mat6(3,1:6)
        mat6(4,1:6) = mat6(4,1:6) + a_y__dx()         * mat6(1,1:6) + a_y__dy()         * mat6(3,1:6)
      endif      

! Drift_1

      end%vec(1) = end%vec(1) + ds2 * end%vec(2) / rel_E
      end%vec(5) = end%vec(5) - ds2 * end%vec(2)**2 / (2*rel_E2)

      if (calc_mat6) then
        mat6(1,1:6) = mat6(1,1:6) + (ds2 / rel_E)          * mat6(2,1:6) - (ds2*end%vec(2)/rel_E2)    * mat6(6,1:6) 
        mat6(5,1:6) = mat6(5,1:6) - (ds2*end%vec(2)/rel_E2) * mat6(2,1:6) + (ds2*end%vec(2)**2/rel_E3) * mat6(6,1:6)
      endif      

! s half step

      s = s + ds2

      if (track_com%save_track) then
        track_com%s(i) = s
        track_com%orb(i) = end
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
! unknown element

  case default

    print *, 'ERROR IN CALC_MAT6_SYMP_LIE_BMAD: NOT YET IMPLEMENTED:', ele%key
    print *, '      FOR ELEMENT: ', ele%name
    call err_exit

  end select

! element offset

  if (calc_mat6) then
    call drift_mat6_calc (m6, -ele%value(s_offset$), end%vec)
    mat6(1,1:6) = mat6(1,1:6) + m6(1,2) * mat6(2,1:6) + m6(1,6) * mat6(6,1:6)
    mat6(3,1:6) = mat6(3,1:6) + m6(3,4) * mat6(4,1:6) + m6(3,6) * mat6(6,1:6)
    mat6(5,1:6) = mat6(5,1:6) + m6(5,2) * mat6(2,1:6) + m6(5,4) * mat6(4,1:6) + m6(5,6) * mat6(6,1:6)

    if (ele%value(tilt$) /= 0) call tilt_mat6 (mat6, ele%value(tilt$))
  endif

  if (real_track) call offset_particle (ele, param, end, unset$, set_canonical = .false.)

end subroutine

!----------------------------------------------------------------------------
! contains

subroutine update_coefs

  real(rp) factor, coef

  factor = c_light / param%beam_energy

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

  if (.not. calc_mat6) return

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
! contains

subroutine update_y_terms

  real(rp) kyy

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
! contains

subroutine update_x_s_terms

  real(rp) kxx, kzz

  do j = 1, size(ele%wig_term)
    wt => ele%wig_term(j)

    kxx = wt%kx * end%vec(1)
    if (abs(kxx) < 1e-20) then
      tm(j)%c_x = 1
      tm(j)%s_x = kxx
      tm(j)%s_x_kx = end%vec(1)
    elseif (wt%type == hyper_x$ .or. wt%type == hyper_xy$) then
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
! contains

function a_y() result (value)

  real(rp) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%a_y%coef * tm(j)%s_x_kx * tm(j)%s_y_ky * tm(j)%s_z
  enddo

end function

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
! contains

function da_z_dx() result (value)

  real(rp) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%da_z_dx%coef * tm(j)%c_x * tm(j)%c_y * tm(j)%c_z
  enddo

end function

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
! contains

function dint_a_y_dx__dx() result (value)

  real(rp) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%dint_a_y_dx%dx_coef * tm(j)%s_x * tm(j)%c1_ky2 * tm(j)%s_z
  enddo

end function

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
! contains

function a_y__dx() result (value)

  real(rp) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%a_y%dx_coef * tm(j)%c_x * tm(j)%s_y_ky * tm(j)%s_z
  enddo

end function

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
! contains

function da_z_dx__dx() result (value)

  real(rp) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%da_z_dx%dx_coef * tm(j)%s_x * tm(j)%c_y * tm(j)%c_z
  enddo

end function

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
! contains

function da_z_dy__dx() result (value)

  real(rp) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%da_z_dy%dx_coef * tm(j)%c_x * tm(j)%s_y * tm(j)%c_z
  enddo

end function

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
