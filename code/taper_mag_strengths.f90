!+
! Subroutine taper_mag_strengths (lat, ref_lat, vary_solenoid)
!
! Routine to "correct" magnet strengths in a ring to counteract local radiation damping energy offsets.
!
! Routine will scale magnet strengths around a ring according to the local closed orbit momentum 
! so that the normalized strength at the closed orbit momentum is equal to the unshifted normalized 
! strength at the reference momentum. Varied will be "multipolar" like strengths. In particular, varied will be:
!   DG, K1, K2, K3, HKICK, VKICK, KICK, a_pole(:), b_pole(:), KS (if vary_solenoid = T).
! Notice that something like wiggler strengths are not varied.
!
! Radiation damping needs to be turned on before calling this routine.
!
! Input:
!   lat             -- lat_struct: Lattice to vary.
!   ref_lat         -- lat_struct, optional: Reference lattice. If not present, lat will be used as the ref.
!   vary_solenoid   -- logical, optional: Vary solenoid strengths? Default is False.
!
! Output:
!   lat             -- lat_struct: Lattice with magnet strengths varied.
!-

subroutine taper_mag_strengths (lat, ref_lat, vary_solenoid)

use bmad, dummy => taper_mag_strengths

implicit none

type (lat_struct), target :: lat, lat0
type (lat_struct), optional, target :: ref_lat
type (branch_struct), pointer :: branch, branch0
type (ele_struct), pointer :: ele
type (coord_struct), allocatable :: closed_orb(:)

real(rp) tol, max_change, pz_ave, weight
integer k, ib, ie
logical, optional :: vary_solenoid
logical err_flag, vary_ks

character(*), parameter :: r_name = 'taper_mag_strengths'

!

tol = 1e-6
vary_ks = logic_option(.false., vary_solenoid)

if (present(ref_lat)) then
  lat0 = ref_lat
else
  lat0 = lat
endif

branch_loop: do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  branch0 => lat0%branch(ib)
  if (branch%param%geometry == open$) cycle

  do k = 1, 10
    call closed_orbit_calc (lat, closed_orb, 6, +1, ib, err_flag)
    if (err_flag) then
      call out_io (s_error$, r_name, 'CLOSED ORBIT CALC FAILED ON LOOP: ' // int_str(k))
      return
    endif

    pz_ave = 0
    weight = 0 
    do ie = 1, branch%n_ele_track
      ele => branch%ele(ie)
      if (ele%key /= sbend$) cycle
      pz_ave = pz_ave + 0.5_rp * ele%value(angle$) * (closed_orb(ie-1)%vec(6) + closed_orb(ie)%vec(6))
      weight = weight + ele%value(angle$)
    enddo
    pz_ave = pz_ave / weight

    max_change = 0
    do ie = 1, branch%n_ele_track
      max_change = max(max_change, taper_this_ele (branch%ele(ie), lat0, closed_orb(ie), pz_ave, vary_ks))
    enddo

    call lattice_bookkeeper(lat)
    if (max_change < tol .and. k > 1) cycle branch_loop
  enddo

  call out_io (s_warn$, r_name, 'CALCULATION NOT CONVERGING. MAX CHANGE ON LAST LOOP IS: ' // real_str(max_change, 3))
enddo branch_loop

!-----------------------------------------------
contains

recursive function taper_this_ele (ele_in, lat0, orb, pz_ave, vary_ks) result (change)

type (ele_struct), target :: ele_in
type (lat_struct) lat0
type (ele_struct), pointer :: ele, ele0, lord
type (coord_struct) :: orb

real(rp) change, pz_ave
integer i, ip
logical vary_ks

!

ele => ele_in
change = 0

if (ele%slave_status == super_slave$) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%lord_status /= super_lord$) return
    change = max(change, taper_this_ele(lord, lat0, orb, pz_ave, vary_ks))
  enddo
  return
endif

if (ele%slave_status == multipass_slave$) ele => pointer_to_lord(ele, 1)
ele0 => pointer_to_ele(lat0, ele)

if (has_attribute(ele, 'DG'))    change = max(change, taper_this_attrib(ele, ele%value(dg$), ele0%value(dg$), orb, pz_ave, ele%value(g$)))
if (has_attribute(ele, 'K1'))    change = max(change, taper_this_attrib(ele, ele%value(k1$), ele0%value(k1$), orb, pz_ave))
if (has_attribute(ele, 'K2'))    change = max(change, taper_this_attrib(ele, ele%value(k2$), ele0%value(k2$), orb, pz_ave))
if (has_attribute(ele, 'K3'))    change = max(change, taper_this_attrib(ele, ele%value(k3$), ele0%value(k3$), orb, pz_ave))
if (has_attribute(ele, 'HKICK')) change = max(change, taper_this_attrib(ele, ele%value(hkick$), ele0%value(hkick$), orb, pz_ave))
if (has_attribute(ele, 'VKICK')) change = max(change, taper_this_attrib(ele, ele%value(vkick$), ele0%value(vkick$), orb, pz_ave))
if (has_attribute(ele, 'KICK'))  change = max(change, taper_this_attrib(ele, ele%value(kick$), ele0%value(kick$), orb, pz_ave))

if (associated(ele%a_pole)) then
  do ip = 0, ubound(ele%a_pole, 1)
    change = max(change, taper_this_attrib(ele, ele%a_pole(ip), ele0%a_pole(ip), orb, pz_ave))
    change = max(change, taper_this_attrib(ele, ele%b_pole(ip), ele0%b_pole(ip), orb, pz_ave))
  enddo
endif

if (vary_ks .and. has_attribute(ele, 'KS')) change = max(change, taper_this_attrib(ele, ele%value(ks$), ele0%value(ks$), orb, pz_ave))

end function taper_this_ele

!-----------------------------------------------
! contains

function taper_this_attrib(ele, attrib, attrib0, orb, pz_ave, g0) result (change)

type (ele_struct) :: ele
type (coord_struct) orb
real(rp), optional :: g0
real(rp) attrib, attrib0, pz_ave, change, base
real(rp) f, a_new

!

change = 0
base = real_option(0.0_rp, g0)

if (attrib0 + base == 0) return
a_new = (attrib0 + base) * (1 + orb%vec(6)) / (1 + pz_ave) - base
change = (a_new - attrib) / (attrib0 + base)
change = abs(change)
attrib = a_new

call set_flags_for_changed_attribute (ele, attrib)

end function taper_this_attrib

end subroutine
