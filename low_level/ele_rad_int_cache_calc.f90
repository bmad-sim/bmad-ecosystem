!+
! Subroutine ele_rad_int_cache_calc (ele)
! 
! Routine to calculate the radiation integrals for a lattice element needed when tracking with radiation.
!
! Input:
!   ele          -- ele_struct: Lattice element 
!
! Output:
!   ele          -- ele_struct: Lattice element 
!     %rad_int_cache
!-

subroutine ele_rad_int_cache_calc (ele)

use pointer_to_ele_mod, dummy => ele_rad_int_cache_calc

implicit none

type (ele_struct) ele
type (ele_struct), pointer :: field_ele
type (coord_struct) start0_orb, start_orb, end_orb
type (track_struct) track
type (branch_struct), pointer :: branch

real(rp) g2, g3
real(rp), parameter :: del_orb = 1d-4

integer j, track_method_saved
logical err_flag

!

if (ele%key /= wiggler$ .and. ele%key /= undulator$ .and. ele%key /= em_field$) return
if (ele%value(l$) == 0) return

field_ele => pointer_to_field_ele(ele, 1)
if (field_ele%field_calc == planar_model$ .or. field_ele%field_calc == helical_model$) return

if (.not. associated(ele%rad_int_cache)) allocate (ele%rad_int_cache)
if (.not. ele%rad_int_cache%stale .and. all(ele%map_ref_orb_in%vec == ele%rad_int_cache%rm0%ref_orb)) return

! If %map_ref_orb_in has not been set (can happen before the closed orbit is 
! calculated or during optimization), use %time_ref_orb_in.

if (ele%map_ref_orb_in%state == not_set$ .or. ele%map_ref_orb_in%p0c /= ele%value(p0c$)) then
  start0_orb = ele%time_ref_orb_in
else
  start0_orb = ele%map_ref_orb_in
endif
ele%rad_int_cache%rm0%ref_orb = start0_orb%vec

branch => pointer_to_branch(ele)

track%n_pt = -1
track_method_saved = ele%tracking_method
if (ele%tracking_method == taylor$) ele%tracking_method = runge_kutta$
call track1 (start0_orb, ele, branch%param, end_orb, track, err_flag, .true.)
call calc_this_g (ele, branch, track, ele%rad_int_cache%g2_0, ele%rad_int_cache%g3_0)

do j = 1, 4
  start_orb = start0_orb
  start_orb%vec(j) = start_orb%vec(j) + del_orb
  track%n_pt = -1
  call track1 (start_orb, ele, branch%param, end_orb, track, err_flag, .true.)
  call calc_this_g (ele, branch, track, g2, g3)
  ele%rad_int_cache%dg2_dorb(j) = (g2 - ele%rad_int_cache%g2_0) / del_orb
  ele%rad_int_cache%dg3_dorb(j) = (g3 - ele%rad_int_cache%g3_0) / del_orb
enddo

ele%rad_int_cache%stale = .false.
ele%tracking_method = track_method_saved

!-------------------------------------------------------
contains

subroutine calc_this_g (ele, branch, track, g2, g3)

type (ele_struct) ele
type (branch_struct) branch
type (track_struct) track
real(rp) g2, g3, g2_here, g3_here, g(3), f, s0
integer j, n1

! g2 is the average g^2 over the element for an on-energy particle.

track%pt(:)%orb%vec(6) = 0  ! on-energy

g2 = 0; g3 = 0

n1 = track%n_pt
s0 = ele%s_start

do j = 0, n1

  call g_bending_strength_from_em_field (ele, branch%param, track%pt(j)%orb%s - s0, track%pt(j)%orb, .false., g)

  g2_here = g(1)**2 + g(2)**2 ! = g_x^2 + g_y^2
  g3_here = sqrt(g2_here)**3

  if (j == 0) then
    g2 = g2 + g2_here * track%pt(1)%s_body / 2.0_rp
    g3 = g3 + g3_here * track%pt(1)%s_body / 2.0_rp
  elseif (j == n1) then
    g2 = g2 + g2_here * (track%pt(n1)%s_body - track%pt(n1-1)%s_body) / 2.0_rp
    g3 = g3 + g3_here * (track%pt(n1)%s_body - track%pt(n1-1)%s_body) / 2.0_rp
  else
    g2 = g2 + g2_here * (track%pt(j-1)%s_body - track%pt(j+1)%s_body) / 2.0_rp
    g3 = g3 + g3_here * (track%pt(j-1)%s_body - track%pt(j+1)%s_body) / 2.0_rp
  endif


enddo

g2 = g2 / track%pt(n1)%s_body
g3 = g3 / track%pt(n1)%s_body

end subroutine calc_this_g

end subroutine
