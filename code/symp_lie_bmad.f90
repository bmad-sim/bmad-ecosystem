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
!-

#include "CESR_platform.inc"

subroutine symp_lie_bmad (ele, param, start, end, calc_mat6)

  use bmad_struct
  use bmad_interface
  use make_mat6_mod

  implicit none

  type save_coef_struct
    real(rdef) coef, dx_coef, dy_coef
  end type

  type save_computations_struct
    type (save_coef_struct) a_y, dint_a_y_dx, da_z_dx, da_z_dy
    real(rdef) c_x, s_x, c_y, s_y, c_z, s_z, s_x_kx, s_y_ky, c1_ky2
  end type

  type (save_computations_struct), allocatable, save :: tm(:)

  type (ele_struct), target :: ele
  type (coord_struct) :: start, end, orb
  type (param_struct)  param
  type (wig_term_struct), pointer :: wt

  real(rdef) rel_E, rel_E2, rel_E3, ds, ds2, s, m6(6,6)
  real(rdef), pointer :: mat6(:,:)
  real(rdef), parameter :: z0 = 0, z1 = 1

  integer i, j

  logical calc_mat6

! init

  if (calc_mat6) mat6 => ele%mat6

  rel_E = (1 + start%z%vel)
  rel_E2 = rel_E**2
  rel_E3 = rel_E**3

  end = start

! element offset 

  if (calc_mat6) then
    call drift_mat6_calc (mat6, ele%value(s_offset$), end%vec)
    call mat6_add_tilt_at_end (mat6, ele%value(tilt$))
  endif

  call offset_particle (ele, param, end, set$, set_canonical = .false.)

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

    do i = 1, ele%num_steps

! check for overflow

      if (maxval(abs(ele%wig_term(:)%kx * end%x%pos)) > 30 .or. &
                  maxval(abs(ele%wig_term(:)%ky * end%y%pos)) > 30) then
        print *, 'ERROR IN SYMP_LIE_BMAD: ', &
                                     'FLOATING OVERFLOW IN WIGGLER TRACKING.'
        print *, '      PARTICLE WILL BE TAGGED AS LOST.'
        param%lost = .true.
        return
      endif

! s half step

      s = s + ds2

! Drift_1 = P_x^2 / (2 * (1 + dE))

      end%x%pos = end%x%pos + ds2 * end%x%vel / rel_E
      end%z%pos = end%z%pos - ds2 * end%x%vel**2 / (2*rel_E2)

      if (calc_mat6) then
        mat6(1,1:6) = mat6(1,1:6) + (ds2 / rel_E)          * mat6(2,1:6) - (ds2*end%x%vel/rel_E2)    * mat6(6,1:6) 
        mat6(5,1:6) = mat6(5,1:6) - (ds2*end%x%vel/rel_E2) * mat6(2,1:6) + (ds2*end%x%vel**2/rel_E3) * mat6(6,1:6)
      endif

! Drift_2 = (P_y - a_y)**2 / (2 * (1 + dE))

      call update_x_s_terms

      end%x%vel = end%x%vel - dint_a_y_dx()
      end%y%vel = end%y%vel - a_y()

      if (calc_mat6) then
        mat6(2,1:6) = mat6(2,1:6) - dint_a_y_dx__dx() * mat6(1,1:6) - dint_a_y_dx__dy() * mat6(3,1:6)
        mat6(4,1:6) = mat6(4,1:6) - a_y__dx()         * mat6(1,1:6) - a_y__dy()         * mat6(3,1:6)
      endif      

!

      end%y%pos = end%y%pos + ds2 * end%y%vel / rel_E
      end%z%pos = end%z%pos - ds2 * end%y%vel**2 / (2*rel_E2)

      if (calc_mat6) then
        mat6(3,1:6) = mat6(3,1:6) + (ds2 / rel_E)          * mat6(4,1:6) - (ds2*end%y%vel/rel_E2)    * mat6(6,1:6) 
        mat6(5,1:6) = mat6(5,1:6) - (ds2*end%y%vel/rel_E2) * mat6(4,1:6) + (ds2*end%y%vel**2/rel_E3) * mat6(6,1:6)
      endif      

!

      call update_y_terms

      end%x%vel = end%x%vel + dint_a_y_dx()
      end%y%vel = end%y%vel + a_y()

      if (calc_mat6) then
        mat6(2,1:6) = mat6(2,1:6) + dint_a_y_dx__dx() * mat6(1,1:6) + dint_a_y_dx__dy() * mat6(3,1:6)
        mat6(4,1:6) = mat6(4,1:6) + a_y__dx()         * mat6(1,1:6) + a_y__dy()         * mat6(3,1:6)
      endif      

! Kick = a_z

      end%x%vel = end%x%vel + ds * da_z_dx()
      end%y%vel = end%y%vel + ds * da_z_dy()

      if (calc_mat6) then
        mat6(2,1:6) = mat6(2,1:6) + ds * da_z_dx__dx() * mat6(1,1:6) + ds * da_z_dx__dy() * mat6(3,1:6)
        mat6(4,1:6) = mat6(4,1:6) + ds * da_z_dy__dx() * mat6(1,1:6) + ds * da_z_dy__dy() * mat6(3,1:6)
      endif 

! Drift_2

      end%x%vel = end%x%vel - dint_a_y_dx()
      end%y%vel = end%y%vel - a_y()

      if (calc_mat6) then
        mat6(2,1:6) = mat6(2,1:6) - dint_a_y_dx__dx() * mat6(1,1:6) - dint_a_y_dx__dy() * mat6(3,1:6)
        mat6(4,1:6) = mat6(4,1:6) - a_y__dx()         * mat6(1,1:6) - a_y__dy()         * mat6(3,1:6)
      endif      

!

      end%y%pos = end%y%pos + ds2 * end%y%vel / rel_E
      end%z%pos = end%z%pos - ds2 * end%y%vel**2 / (2*rel_E2)

      if (calc_mat6) then
        mat6(3,1:6) = mat6(3,1:6) + (ds2 / rel_E)          * mat6(4,1:6) - (ds2*end%y%vel/rel_E2)    * mat6(6,1:6) 
        mat6(5,1:6) = mat6(5,1:6) - (ds2*end%y%vel/rel_E2) * mat6(4,1:6) + (ds2*end%y%vel**2/rel_E3) * mat6(6,1:6)
      endif      

!

      call update_y_terms

      end%x%vel = end%x%vel + dint_a_y_dx()
      end%y%vel = end%y%vel + a_y()

      if (calc_mat6) then
        mat6(2,1:6) = mat6(2,1:6) + dint_a_y_dx__dx() * mat6(1,1:6) + dint_a_y_dx__dy() * mat6(3,1:6)
        mat6(4,1:6) = mat6(4,1:6) + a_y__dx()         * mat6(1,1:6) + a_y__dy()         * mat6(3,1:6)
      endif      

! Drift_1

      end%x%pos = end%x%pos + ds2 * end%x%vel / rel_E
      end%z%pos = end%z%pos - ds2 * end%x%vel**2 / (2*rel_E2)

      if (calc_mat6) then
        mat6(1,1:6) = mat6(1,1:6) + (ds2 / rel_E)          * mat6(2,1:6) - (ds2*end%x%vel/rel_E2)    * mat6(6,1:6) 
        mat6(5,1:6) = mat6(5,1:6) - (ds2*end%x%vel/rel_E2) * mat6(2,1:6) + (ds2*end%x%vel**2/rel_E3) * mat6(6,1:6)
      endif      

! s half step

      s = s + ds2

    enddo

!----------------------------------------------------------------------------
! unknown element

  case default

    type *, 'ERROR IN CALC_MAT6_SYMP_LIE_BMAD: NOT YET IMPLEMENTED:', ele%key
    type *, '      FOR ELEMENT: ', ele%name
    call err_exit

  end select

! element offset

  if (calc_mat6) then
    call drift_mat6_calc (m6, -ele%value(s_offset$), end%vec)
    mat6(1,1:6) = mat6(1,1:6) + m6(1,2) * mat6(2,1:6) + m6(1,6) * mat6(6,1:6)
    mat6(3,1:6) = mat6(3,1:6) + m6(3,4) * mat6(4,1:6) + m6(3,6) * mat6(6,1:6)
    mat6(5,1:6) = mat6(5,1:6) + m6(5,2) * mat6(2,1:6) + m6(5,4) * mat6(4,1:6) + m6(5,6) * mat6(6,1:6)

    call mat6_add_tilt_at_end (mat6, -ele%value(tilt$))
  endif

  call offset_particle (ele, param, end, unset$, set_canonical = .false.)

!----------------------------------------------------------------------------
contains

subroutine update_coefs

  real(rdef) factor

  factor = c_light / param%beam_energy

  do j = 1, size(ele%wig_term)
    wt => ele%wig_term(j)
    tm(j)%a_y%coef     = -factor * wt%coef * wt%kz          ! / (wt%kx * wt%ky)
    tm(j)%dint_a_y_dx%coef = -factor * wt%coef * wt%kz      ! / wt%ky**2
    tm(j)%da_z_dx%coef = -factor * wt%coef 
    tm(j)%da_z_dy%coef = -factor * wt%coef * wt%ky          ! / wt%kx
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

  real(rdef) kyy

  do j = 1, size(ele%wig_term)
    wt => ele%wig_term(j)
    kyy = wt%ky * end%y%pos
    if (abs(kyy) < 1e-20) then
      tm(j)%c_y = 1
      tm(j)%s_y = kyy
      tm(j)%s_y_ky = end%y%pos
      tm(j)%c1_ky2 = end%y%pos**2 / 2
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

  real(rdef) kxx, kzz

  do j = 1, size(ele%wig_term)
    wt => ele%wig_term(j)

    kxx = wt%kx * end%x%pos
    if (abs(kxx) < 1e-20) then
      tm(j)%c_x = 1
      tm(j)%s_x = kxx
      tm(j)%s_x_kx = end%x%pos
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

  real(rdef) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%a_y%coef * tm(j)%s_x_kx * tm(j)%s_y_ky * tm(j)%s_z
  enddo

end function

!----------------------------------------------------------------------------
! contains

function dint_a_y_dx() result (value)

  real(rdef) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%dint_a_y_dx%coef * tm(j)%c_x * tm(j)%c1_ky2 * tm(j)%s_z
  enddo

end function

!----------------------------------------------------------------------------
! contains

function da_z_dx() result (value)

  real(rdef) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%da_z_dx%coef * tm(j)%c_x * tm(j)%c_y * tm(j)%c_z
  enddo

end function

!----------------------------------------------------------------------------
! contains

function da_z_dy() result (value)

  real(rdef) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%da_z_dy%coef * tm(j)%s_x_kx * tm(j)%s_y * tm(j)%c_z
  enddo

end function

!----------------------------------------------------------------------------
! contains

function dint_a_y_dx__dx() result (value)

  real(rdef) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%dint_a_y_dx%dx_coef * tm(j)%s_x * tm(j)%c1_ky2 * tm(j)%s_z
  enddo

end function

!----------------------------------------------------------------------------
! contains

function dint_a_y_dx__dy() result (value)

  real(rdef) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%dint_a_y_dx%dy_coef * tm(j)%c_x * tm(j)%s_y_ky * tm(j)%s_z
  enddo

end function

!----------------------------------------------------------------------------
! contains

function a_y__dx() result (value)

  real(rdef) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%a_y%dx_coef * tm(j)%c_x * tm(j)%s_y_ky * tm(j)%s_z
  enddo

end function

!----------------------------------------------------------------------------
! contains

function a_y__dy() result (value)

  real(rdef) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%a_y%dy_coef * tm(j)%s_x_kx * tm(j)%c_y * tm(j)%s_z
  enddo

end function

!----------------------------------------------------------------------------
! contains

function da_z_dx__dx() result (value)

  real(rdef) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%da_z_dx%dx_coef * tm(j)%s_x * tm(j)%c_y * tm(j)%c_z
  enddo

end function

!----------------------------------------------------------------------------
! contains

function da_z_dx__dy() result (value)

  real(rdef) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%da_z_dx%dy_coef * tm(j)%c_x * tm(j)%s_y * tm(j)%c_z
  enddo

end function

!----------------------------------------------------------------------------
! contains

function da_z_dy__dx() result (value)

  real(rdef) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%da_z_dy%dx_coef * tm(j)%c_x * tm(j)%s_y * tm(j)%c_z
  enddo

end function

!----------------------------------------------------------------------------
! contains

function da_z_dy__dy() result (value)

  real(rdef) value

  value = 0
  do j = 1, size(ele%wig_term)
    value = value + tm(j)%da_z_dy%dy_coef * tm(j)%s_x_kx * tm(j)%c_y * tm(j)%c_z
  enddo

end function

end subroutine

