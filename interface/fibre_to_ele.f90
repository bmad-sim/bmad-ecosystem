!+
! Subroutine fibre_to_ele (ptc_fibre, ele_array, ix_ele, err_flag)
!
! Routine to transfer parameters from a PTC fibre to one, or more if needed, bmad lattice elements.
! For example, a fibre with a patch needs multiple elements.
!
! Module Needed:
!   use ptc_layout_mod
!
! Input:
!   ptc_fibre     -- fibre: PTC fiber.
!   ix_ele        -- integer: Index in ele(:) array of element to use
!   ele_array(0:) -- ele_struct: Array of Bmad elements. Indexed from 0.
!
! Output:
!   ele_array(0:) -- ele_struct: Array of Bmad elements. Indexed from 0.
!   ix_ele        -- integer: Index to next unused element.
!   err_flag      -- logical: Set true if there is an error. False otherwise.
!-

! To do: lcavity energy change !?
! open or closed geometry?
! Energy patch


subroutine fibre_to_ele (ptc_fibre, ele_array, ix_ele, err_flag)

use madx_ptc_module, ptc_pi => pi, ptc_twopi => twopi
use bmad

type (fibre), target :: ptc_fibre
type (ele_struct), target :: ele_array(0:)
type (ele_struct), pointer :: ele
type (fibre), pointer :: fib
type (element), pointer :: mag
type (work) wk

real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
real(rp) knl(0:n_pole_maxx), tn(0:n_pole_maxx), ele_tilt, this_kick

integer ix_ele, nmul, i, ix

logical err_flag

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
    call ele_out (fib%mag%name)
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

call ele_out (fib%mag%name)

if (fib%patch%patch == 2 .or. fib%patch%patch == 3) then
  call patch_out ('EXIT_PATCH_' // fib%mag%name, fib%patch%b_d, fib%patch%b_ang)
endif

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
contains

subroutine patch_out (name, offset, angles)

character(*) name
real(rp) offset(3), angles(3)

!

ele => ele_array(ix_ele)
ele%name = name
ix_ele = ix_ele + 1

call set_energy (ele)

end subroutine patch_out

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
! contains

subroutine ele_out (name)

character(*) name

!

ele => ele_array(ix_ele)
ele%name = name
ix_ele = ix_ele + 1

call set_energy (ele)

!

if (associated(mag%d0)) then
  ele%key = drift$
else
endif

!  call is_associated (associated(mag%k2),     '%k2     [Integrator]:')
!  call is_associated (associated(mag%k3),     '%k3     [Thin Kick]:')
!  call is_associated (associated(mag%c4),     '%c4     [Cavity]:')
!  call is_associated (associated(mag%s5),     '%s5     [Solenoid]:')
!  call is_associated (associated(mag%t6),     '%t6     [Integrator thick,slow]:')
!  call is_associated (associated(mag%t7),     '%t7     [Integrator thick,fast]:')
!  call is_associated (associated(mag%s8),     '%s8     [Normal SMI]:')
!  call is_associated (associated(mag%s9),     '%s9     [Skew SMI]:')
!  call is_associated (associated(mag%tp10),   '%tp10   [Sector Teapot]:')
!  call is_associated (associated(mag%mon14),  '%mon14  [Monitor/Instrument]:')
!  call is_associated (associated(mag%sep15),  '%sep15  [Monitor/Instrument]:')
!  call is_associated (associated(mag%k16),    '%k16    [Exact Straight Integrator]:')
!  call is_associated (associated(mag%enge17), '%enge17 [Solenoid Sixtrack style]:')
!  call is_associated (associated(mag%rcol18), '%rcol18 [Rcollimator]:')
!  call is_associated (associated(mag%ecol19), '%ecol19 [Ecollimator]:')
!  call is_associated (associated(mag%cav21),  '%cav21  [Cavity. Traveling Wave]:')
!  call is_associated (associated(mag%wi),     '%wi     [Wiggler]:')
!  call is_associated (associated(mag%pa),     '%pa     [General B]:')
!  call is_associated (associated(mag%he22),   '%he22   [Helical Dipole]:')



call update_this_real (ele%value(l$), fib%mag%p%ld)

if (attribute_name(ele, num_steps$) == 'NUM_STEPS') then
  call update_this_real (ele%value(num_steps$), real(fib%mag%p%nst, rp))
  if (ele%value(num_steps$) == 0) then
    ele%value(ds_step$) = 0
  else
    ele%value(ds_step$) = ele%value(l$) / ele%value(num_steps$)
  endif
endif

! If integrator_order is defined for this element then update

name = attribute_name(ele, integrator_order$)
if (name(1:1) /= '!') call update_this_real (ele%value(integrator_order$), real(fib%mag%p%method, rp))

! Multipole

a_pole = 0
b_pole = 0
nmul = min(fib%mag%p%nmul, n_pole_maxx+1)
a_pole(0:nmul-1) = fib%mag%an(1:nmul)
b_pole(0:nmul-1) = fib%mag%bn(1:nmul)

call multipole_ab_to_kt (a_pole, b_pole, knl, tn)

! Electric Multipole

if (associated(fib%mag%tp10)) then
  if (associated(fib%mag%tp10%ae)) then
    if (.not. associated(ele%a_pole_elec)) call elec_multipole_init(ele)
    nmul = min(size(fib%mag%tp10%ae), n_pole_maxx+1)
    ele%a_pole_elec(0:nmul-1) = 1d9 * VOLT_C * fib%mag%tp10%ae(1:nmul)
    ele%b_pole_elec(0:nmul-1) = 1d9 * VOLT_C * fib%mag%tp10%be(1:nmul)
  endif
endif

!

if (ele%key == sbend$) then
  call update_this_real (ele%value(ref_tilt_tot$), fib%mag%p%tiltd)
else
  call update_this_real (ele%value(tilt_tot$), fib%mag%p%tiltd)
endif

!

select case (ele%key)
case (ab_multipole$)
  ele%a_pole = a_pole
  ele%b_pole = b_pole

case (drift$)

! Use dsin & dcos due to bug in ifort 13.1 compiler. 
! Can rename to sin & cos when bug is fixed.

case (elseparator$)
  this_kick = fib%mag%volt * 1d6 / ele%value(e_tot$)
  if (fib%charge < 0) this_kick = -this_kick
  ele_tilt = fib%mag%p%tiltd - ele%value(tilt_tot$)
  call update_this_real (ele%value(hkick$), -this_kick * dsin(ele_tilt))
  call update_this_real (ele%value(vkick$), -this_kick * dcos(ele_tilt))

case (hkicker$)
  call update_this_real (ele%value(kick$), knl(1))

case (lcavity$, rfcavity$)
  call update_this_real (ele%value(rf_frequency$), fib%mag%freq)

  select case (nint(ele%value(cavity_type$)))
  case (traveling_wave$);     call update_this_real (ele%value(voltage$), fib%mag%volt*1d6)
  case (standing_wave$);      call update_this_real (ele%value(voltage$), fib%mag%volt*0.5d6)
  case (ptc_standard$);       call update_this_real (ele%value(voltage$), fib%mag%volt*1d6)
  end select

  if (ele%key == lcavity$) then
    call update_this_real (ele%value(phi0$), pi/2 - fib%mag%lag/twopi)
  else
    call update_this_real (ele%value(phi0$), fib%mag%lag/twopi)
  endif

case (multipole$)
  ele%a_pole = knl
  ele%b_pole = tn

case (octupole$)
  call update_this_real (ele%value(k3$), knl(3))
  call update_this_real (ele%value(tilt$), tn(3))
  knl(3) = 0
  tn(3) = 0

case (quadrupole$)
  call update_this_real (ele%value(k1$), knl(1))
  call update_this_real (ele%value(tilt$), tn(1))
  knl(1) = 0
  tn(1) = 0

case (sbend$)
  call update_this_real (ele%value(g$), fib%mag%p%b0)
  call update_this_real (ele%value(angle$), ele%value(g$) * ele%value(l$))
  ix = nint(ele%value(ptc_field_geometry$))
  if (ix == straight$ .or. ix == true_rbend$) then
    call update_this_real (ele%value(e1$), fib%mag%p%edge(1) + ele%value(angle$)/2)
    call update_this_real (ele%value(e2$), fib%mag%p%edge(2) + ele%value(angle$)/2)
  else
    call update_this_real (ele%value(e1$), fib%mag%p%edge(1))
    call update_this_real (ele%value(e2$), fib%mag%p%edge(2))
  endif
  call update_this_real (ele%value(hgap$), fib%mag%hgap)
  call update_this_real (ele%value(fint$), fib%mag%fint)

case (sextupole$)
  call update_this_real (ele%value(k2$), knl(2))
  call update_this_real (ele%value(tilt$), tn(2))
  knl(2) = 0
  tn(2) = 0

case (solenoid$)
  call update_this_real (ele%value(ks$), fib%mag%b_sol)

case (sol_quad$)
  call update_this_real (ele%value(ks$), fib%mag%b_sol)
  call update_this_real (ele%value(k1$), knl(1))
  call update_this_real (ele%value(tilt$), tn(1))
  knl(1) = 0
  tn(1) = 0

case (wiggler$, undulator$)


case default
end select

! multipoles

if (any(knl /= 0)) then
endif


! kicks

! Fringes

end subroutine ele_out

!------------------------------------------------------------------------
! contains

subroutine set_energy (ele)
type (ele_struct) ele

wk = fib
ele%value(p0c_start$) = wk%p0c
ele%value(E_tot_start$) = wk%energy
ele%value(p0c$) = wk%p0c
ele%value(E_tot$) = wk%energy

end subroutine set_energy

!------------------------------------------------------------------------
! contains

subroutine update_this_real (var, value)

real(rp) var, value

if (var == value) return

var = value
call set_flags_for_changed_attribute (ele, var)

end subroutine update_this_real

!------------------------------------------------------------------------
! contains

subroutine update_this_integer (var, value)

integer var, value

if (var == value) return

var = value
call set_flags_for_changed_attribute (ele, var)

end subroutine update_this_integer

end subroutine fibre_to_ele

