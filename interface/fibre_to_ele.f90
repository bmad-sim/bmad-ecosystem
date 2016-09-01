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
use bmad, except_dummy => fibre_to_ele

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

ele => ele_array(ix_ele)

ele%name = name
ix_ele = ix_ele + 1

w_mat = matmul(matmul(w_mat_for_y_pitch(-angles(1)), w_mat_for_x_pitch(-angles(2))), w_mat_for_tilt(angles(3)))
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

character(*) name

real(rp) ab_max, ab
real(rp), parameter :: cm1 = 0.01  ! 1 cm radius
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
integer i, n, ix_max
logical err_flag

!

ele => ele_array(ix_ele)
ele%name = name
ix_ele = ix_ele + 1

call set_energy (ele)

!

if (mag%p%b0 /= 0) then
  ele%key = sbend$

elseif (associated(mag%k3)) then
  ele%key = multipole$

elseif (associated(mag%s5)) then
  ele%key = solenoid$

elseif (associated(mag%c4)) then
  ele%key = rfcavity$


! elseparator ???
! taylor ???

elseif (associated(mag%d0)) then
  if (associated (mag%p%aperture)) then
    if (mag%p%aperture%kind == 1) then
      ele%key = ecollimator$
    elseif (mag%p%aperture%kind == 2) then
      ele%key = rcollimator$
    else
      ele%key = drift$
      call out_io (s_error$, r_name, 'UNKNOWN PTC APERTURE TYPE: \i0\ ', &
                                     'FOR ELEMENT: ' // name, i_array = [fib%mag%p%aperture%kind])
      err_flag = .true.
    endif
  else
    ele%key = drift$
  endif

elseif (associated(mag%k2)) then
  n = size(mag%an)
  ab_max = 0
  ix_max = 1
  do i = 2, n
    ab = cm1**(i-1) * abs(mag%bn(i))
    if (ab <= ab_max) cycle
    ab_max = ab
    ix_max = i
  enddo


  a_pole(0:n-1) = mag%an
  b_pole(0:n-1) = mag%bn

  select case (i)
  case(2)
    ele%key = quadrupole$
    ele%value(k1$) = mag%bn(2)
    b_pole(1) = 0

  case(3)
    ele%key = sextupole$
    ele%value(k1$) = mag%bn(3)
    b_pole(2) = 0

  case(4)
    ele%key = sextupole$
    ele%value(k1$) = mag%bn(4)
    b_pole(3) = 0

  case default
    ele%key = kicker$
    ele%value(hkick$) = mag%bn(1)
    ele%value(vkick$) = mag%an(1)
    a_pole(0) = 0
    b_pole(0) = 0
  end select

  if (any(ele%a_pole /= 0) .or. any(ele%b_pole /= 0)) then
    allocate (ele%a_pole(0:n_pole_maxx), ele%b_pole(0:n_pole_maxx))
    ele%a_pole = a_pole
    ele%b_pole = b_pole
  endif

endif

ele%value(l$) = fib%mag%p%ld

if (attribute_name(ele, num_steps$) == 'NUM_STEPS') then
  ele%value(num_steps$) = real(fib%mag%p%nst, rp)
  if (ele%value(num_steps$) == 0) then
    ele%value(ds_step$) = 0
  else
    ele%value(ds_step$) = ele%value(l$) / ele%value(num_steps$)
  endif
endif

! If integrator_order is defined for this element then update

name = attribute_name(ele, integrator_order$)
if (name(1:1) /= '!') ele%value(integrator_order$) = fib%mag%p%method

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
  ele%value(ref_tilt_tot$) = fib%mag%p%tiltd
else
  ele%value(tilt_tot$) = fib%mag%p%tiltd
endif

!

select case (ele%key)
case (ab_multipole$)
  ele%a_pole = a_pole
  ele%b_pole = b_pole

case (drift$)

! 

case (elseparator$)
  this_kick = fib%mag%volt * 1d6 / ele%value(e_tot$)
  if (fib%charge < 0) this_kick = -this_kick
  ele_tilt = fib%mag%p%tiltd - ele%value(tilt_tot$)
  ele%value(hkick$) = -this_kick * sin(ele_tilt)
  ele%value(vkick$) = -this_kick * cos(ele_tilt)

case (hkicker$)
  ele%value(kick$) = knl(1)

case (lcavity$, rfcavity$)
  ele%value(rf_frequency$) = fib%mag%freq

  select case (nint(ele%value(cavity_type$)))
  case (traveling_wave$);     ele%value(voltage$) = fib%mag%volt*1d6
  case (standing_wave$);      ele%value(voltage$) = fib%mag%volt*0.5d6
  case (ptc_standard$);       ele%value(voltage$) = fib%mag%volt*1d6
  end select

  if (ele%key == lcavity$) then
    ele%value(phi0$) = pi/2 - fib%mag%lag/twopi
  else
    ele%value(phi0$) = fib%mag%lag/twopi
  endif

case (multipole$)
  ele%a_pole = knl
  ele%b_pole = tn

case (octupole$)
  ele%value(k3$) = knl(3)
  ele%value(tilt$) = tn(3)
  knl(3) = 0
  tn(3) = 0

case (quadrupole$)
  ele%value(k1$) = knl(1)
  ele%value(tilt$) = tn(1)
  knl(1) = 0
  tn(1) = 0

case (sbend$)
  ele%value(g$) = fib%mag%p%b0
  ele%value(angle$) = ele%value(g$) * ele%value(l$)
  ix = nint(ele%value(ptc_field_geometry$))
  if (ix == straight$ .or. ix == true_rbend$) then
    ele%value(e1$) = fib%mag%p%edge(1) + ele%value(angle$)/2
    ele%value(e2$) = fib%mag%p%edge(2) + ele%value(angle$)/2
  else
    ele%value(e1$) = fib%mag%p%edge(1)
    ele%value(e2$) = fib%mag%p%edge(2)
  endif
  ele%value(hgap$) = fib%mag%hgap
  ele%value(fint$) = fib%mag%fint

case (sextupole$)
  ele%value(k2$) = knl(2)
  ele%value(tilt$) = tn(2)
  knl(2) = 0
  tn(2) = 0

case (solenoid$)
  ele%value(ks$) = fib%mag%b_sol

case (sol_quad$)
  ele%value(ks$) = fib%mag%b_sol
  ele%value(k1$) = knl(1)
  ele%value(tilt$) = tn(1)
  knl(1) = 0
  tn(1) = 0

case (wiggler$, undulator$)


case default
end select

! Multipoles

if (any(knl /= 0)) then
endif


! kicks

! Fringes

! Apertures

if (associated(mag%p%aperture)) then
  select case (mag%p%aperture%kind)
  case (1)
    ele%aperture_type = elliptical$
  case (2)
    ele%aperture_type = elliptical$
  case default
    call out_io (s_error$, r_name, 'UNKNOWN PTC APERTURE TYPE: \i0\ ', &
                                   'FOR ELEMENT: ' // name, i_array = [fib%mag%p%aperture%kind])
    err_flag = .true.
  end select

  ele%value(x1_limit$) = mag%p%aperture%x + mag%p%aperture%dx
  ele%value(x2_limit$) = mag%p%aperture%x - mag%p%aperture%dx
  ele%value(y1_limit$) = mag%p%aperture%y + mag%p%aperture%dy
  ele%value(y2_limit$) = mag%p%aperture%y - mag%p%aperture%dy
endif

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

end subroutine fibre_to_ele

