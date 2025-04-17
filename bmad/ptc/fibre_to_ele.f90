!+
! Subroutine fibre_to_ele (ptc_fibre, branch, ix_ele, err_flag, from_mad)
!
! Routine to transfer parameters from a PTC fibre to one, or more if needed, bmad lattice elements.
! For example, a fibre with a patch needs multiple elements.
!
! Input:
!   ptc_fibre     -- fibre: PTC fibre.
!   ix_ele        -- integer: Index in ele(:) array of element last used.
!   branch        -- branch_struct: branch containing elements.
!   from_mad      -- logical, optional: If True, ignore PTC specific parameters like integrator_order.
!                      Default is False. True is used when the fibre has been created via MAD. In this
!                      case, the PTC specific parameters may not have good values.
!
! Output:
!   branch        -- branch_struct: branch containing elements.
!   ix_ele        -- integer: Index to element created (upper index if more than one created).
!   err_flag      -- logical: Set true if there is an error. False otherwise.
!-

! To do: lcavity energy change !?
! open or closed geometry?
! Energy patch


subroutine fibre_to_ele (ptc_fibre, branch, ix_ele, err_flag, from_mad)

use madx_ptc_module, ptc_pi => pi, ptc_twopi => twopi
use bmad_routine_interface, except_dummy => fibre_to_ele

implicit none

type (fibre), target :: ptc_fibre
type (ele_struct), pointer :: ele
type (branch_struct) branch
type (fibre), pointer :: fib
type (element), pointer :: mag
type (work) wk

real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
real(rp) knl(0:n_pole_maxx), tn(0:n_pole_maxx), ele_tilt, this_kick

integer ix_ele, nmul, i, ix

logical err_flag
logical, optional :: from_mad

character(*), parameter :: r_name = 'fibre_to_ele'
character(40) name

! Init

fib => ptc_fibre
mag => fib%mag

err_flag = .false.

! Pure patch or marker

if (fib%mag%kind == kind0) then
  select case (fib%patch%patch)
  case (0)  ! marker
    call ele_out (fib%mag%name, err_flag)
  case (1)  ! Entrance patch
    call patch_out (fib%mag%name, fib%patch%a_d, fib%patch%a_ang)
  case (2)  ! Exit patch
    call patch_out (fib%mag%name, fib%patch%b_d, fib%patch%b_ang)
  case default
    call out_io (s_error$, r_name, 'I DO NOT KNOW HOW TO HANDLE PATCH TYPE: \i0\ ', int(fib%patch%patch))
    err_flag = .true.
  end select
  return
endif
  
! Not a pure patch nor a marker

if (fib%patch%patch == 1 .or. fib%patch%patch == 3) then
  call patch_out ('ENT_PATCH_' // fib%mag%name, fib%patch%a_d, fib%patch%a_ang)
endif

call ele_out (fib%mag%name, err_flag)

if (fib%patch%patch == 2 .or. fib%patch%patch == 3) then
  call patch_out ('EXIT_PATCH_' // fib%mag%name, fib%patch%b_d, fib%patch%b_ang)
endif

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
contains

subroutine patch_out (name, offset, angles)

real(rp) offset(3), angles(3)
real(rp) w_mat(3,3), off(3)

character(*) name

!

ix_ele = ix_ele + 1
ele => branch%ele(ix_ele)
call init_ele(ele, patch$, -1, ix_ele, branch)

ele%name = name

w_mat = matmul(matmul(w_mat_for_y_pitch(-angles(2)), w_mat_for_x_pitch(-angles(1))), w_mat_for_tilt(angles(3)))
call floor_w_mat_to_angles (w_mat, ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$))

off = matmul(transpose(w_mat), offset)
ele%value(x_offset$) = off(1)
ele%value(y_offset$) = off(2)
ele%value(z_offset$) = off(3)

call set_energy (ele)

end subroutine patch_out

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
! contains

subroutine ele_out (name, err_flag)

type (all_pointer_struct) a_ptr

real(rp) ab_max, ab
real(rp), parameter :: cm1 = 0.01  ! 1 cm radius
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)

integer i, np, ix_max

logical err_flag, err

character(*) name

!

ix_ele = ix_ele + 1
ele => branch%ele(ix_ele)

np = 0
a_pole = 0
b_pole = 0

if (associated(mag%an)) then
  np = size(mag%an)
  a_pole(0:np-1) = mag%an
  b_pole(0:np-1) = mag%bn
endif

! Sbend

if (mag%p%b0 /= 0) then
  call init_ele(ele, sbend$, -1, ix_ele, branch)
  if (associated(mag%k16)) then
    ele%sub_key = rbend$
  else
    ele%sub_key = sbend$
  endif

  ele%value(g$) = fib%mag%p%b0
  ele%value(angle$) = ele%value(g$) * ele%value(l$)
  ele%value(hgap$) = fib%mag%hgap(1)
  ele%value(fint$) = fib%mag%fint(1)
  ele%value(hgapx$) = fib%mag%hgap(2)
  ele%value(fintx$) = fib%mag%fint(2)
  ele%value(ref_tilt$) = fib%mag%p%tiltd

  if (nint(ele%value(ptc_field_geometry$)) == straight$) then
    ele%value(e1$) = fib%mag%p%edge(1) + ele%value(angle$)/2
    ele%value(e2$) = fib%mag%p%edge(2) + ele%value(angle$)/2
  else
    ele%value(e1$) = fib%mag%p%edge(1)
    ele%value(e2$) = fib%mag%p%edge(2)
  endif

  if (np >= 2) ele%value(k1$) = mag%bn(2)
  if (np >= 3) ele%value(k2$) = 2 * mag%bn(3)

! Marker

elseif (fib%mag%kind == kind0) then
  call init_ele(ele, marker$, -1, ix_ele, branch)

! Elseparator

elseif (fib%mag%kind == kind15) then
  call init_ele(ele, elseparator$, -1, ix_ele, branch)

  this_kick = fib%mag%volt * 1d6 / ele%value(e_tot$)
  if (fib%charge < 0) this_kick = -this_kick
  ele_tilt = fib%mag%p%tiltd - ele%value(tilt_tot$)
  ele%value(hkick$) = -this_kick * sin(ele_tilt)
  ele%value(vkick$) = -this_kick * cos(ele_tilt)

! Multipole

elseif (associated(mag%k3)) then
  call init_ele(ele, ab_multipole$, -1, ix_ele, branch)
  ele%a_pole = a_pole
  ele%b_pole = b_pole

! Solenoid

elseif (associated(mag%s5)) then
  call init_ele(ele, solenoid$, -1, ix_ele, branch)
  ele%value(ks$) = fib%mag%b_sol

! RFcavity

elseif (associated(mag%c4)) then
  call init_ele(ele, rfcavity$, -1, ix_ele, branch)

  ele%value(rf_frequency$) = fib%mag%freq
  ele%value(voltage$) = 1d6 * fib%mag%p%ld * fib%mag%volt

  if (ele%key == lcavity$) then
    ele%value(phi0$) = pi/2 - fib%mag%lag/twopi
  else
    ele%value(phi0$) = fib%mag%lag/twopi
  endif

! Collimator

elseif (associated(mag%d0)) then
  if (associated (mag%p%aperture)) then
    if (mag%p%aperture%kind == 1) then
      call init_ele(ele, ecollimator$, -1, ix_ele, branch)
    elseif (mag%p%aperture%kind == 2) then
      call init_ele(ele, rcollimator$, -1, ix_ele, branch)
    else
      call init_ele(ele, drift$, -1, ix_ele, branch)
      call out_io (s_error$, r_name, 'UNKNOWN PTC APERTURE TYPE: \i0\ ', &
                                     'FOR ELEMENT: ' // name, i_array = [fib%mag%p%aperture%kind])
      err_flag = .true.
    endif
  else
    call init_ele(ele, drift$, -1, ix_ele, branch)
  endif

! k2:  Non-exact, mad_model = 1
! k16: Exact, mad_model = 1 (Drift-Kick)
! t7:  Exact or Non-exact, mad_model = 2 (Matrix-Kick)
! t6:  Exact or Non-exact, mad_model = 3 (Delta_Dependent_Matrix-Kick]

elseif (associated(mag%k2) .or. associated(mag%k16) .or. associated(mag%t7) .or. associated(mag%t6)) then
  ab_max = 0
  ix_max = 1
  do i = 2, np
    ab = cm1**(i-1) * abs(mag%bn(i))
    if (ab <= ab_max) cycle
    ab_max = ab
    ix_max = i
  enddo

  select case (ix_max)
  case(2)
    call init_ele(ele, quadrupole$, -1, ix_ele, branch)
    ele%value(k1$) = mag%bn(2)
    b_pole(1) = 0

  case(3)
    call init_ele(ele, sextupole$, -1, ix_ele, branch)
    ele%value(k2$) = 2 * mag%bn(3)
    b_pole(2) = 0

  case(4)
    call init_ele(ele, octupole$, -1, ix_ele, branch)
    ele%value(k3$) = 6 * mag%bn(4)
    b_pole(3) = 0

  case default
    call init_ele(ele, kicker$, -1, ix_ele, branch)
    ele%value(hkick$) = -fib%mag%p%ld * mag%bn(1)
    ele%value(vkick$) =  fib%mag%p%ld * mag%an(1)
    a_pole(0) = 0
    b_pole(0) = 0
  end select

  if (any(ele%a_pole /= 0) .or. any(ele%b_pole /= 0)) then
    allocate (ele%a_pole(0:n_pole_maxx), ele%b_pole(0:n_pole_maxx))
    ele%a_pole = a_pole
    ele%b_pole = b_pole
  endif

endif

call set_energy (ele)
ele%name = name
ele%value(l$) = fib%mag%p%ld


if (.not. logic_option(.false., from_mad)) then
  call pointer_to_attribute(ele, 'NUM_STEPS', .false., a_ptr, err, .false.)
  if (associated(a_ptr%r)) then
    a_ptr%r = real(fib%mag%p%nst, rp)
    if (a_ptr%r == 0) then
      ele%value(ds_step$) = 0
    else
      ele%value(ds_step$) = ele%value(l$) / a_ptr%r
    endif
  endif
endif

! If integrator_order is defined for this element then update

if (.not. logic_option(.false., from_mad)) then
  call pointer_to_attribute(ele, 'INTEGRATOR_ORDER', .false., a_ptr, err, .false.)
  if (associated(a_ptr%r)) a_ptr%r = fib%mag%p%method
endif

!

if (ele%key /= sbend$) ele%value(tilt$) = fib%mag%p%tiltd

! kicks

! Fringes

! Apertures

if (associated(mag%p%aperture)) then
  select case (mag%p%aperture%kind)
  case (1)
    ele%aperture_type = elliptical$
    ele%value(x1_limit$) = mag%p%aperture%r(1)
    ele%value(x2_limit$) = mag%p%aperture%r(1)
    ele%value(y1_limit$) = mag%p%aperture%r(2)
    ele%value(y2_limit$) = mag%p%aperture%r(2)
  case (2)
    ele%aperture_type = rectangular$
    ele%value(x1_limit$) = mag%p%aperture%x + mag%p%aperture%dx
    ele%value(x2_limit$) = mag%p%aperture%x - mag%p%aperture%dx
    ele%value(y1_limit$) = mag%p%aperture%y + mag%p%aperture%dy
    ele%value(y2_limit$) = mag%p%aperture%y - mag%p%aperture%dy
  case default
    call out_io (s_error$, r_name, 'UNKNOWN PTC APERTURE TYPE: \i0\ ', &
                                   'FOR ELEMENT: ' // name, i_array = [fib%mag%p%aperture%kind])
    err_flag = .true.
  end select

endif

end subroutine ele_out

!------------------------------------------------------------------------
! contains

subroutine set_energy (ele)
type (ele_struct) ele

wk = fib
ele%value(E_tot_start$) = 1d9 * wk%energy
ele%value(E_tot$)       = 1d9 * wk%energy
ele%value(p0c_start$)   = 1d9 * wk%p0c
ele%value(p0c$)         = 1d9 * wk%p0c

end subroutine set_energy

end subroutine fibre_to_ele

