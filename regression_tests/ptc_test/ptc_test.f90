!+
! Program ptc_test
!
! This program is part of the Bmad regression testing suite.
!-

program ptc_test

use bmad, dummy => dp
use s_def_kind, dummy2 => dp
use ptc_layout_mod
use ptc_map_with_radiation_mod
use duan_zhe_map, only: grnf_zhe, get_seed, set_seed

implicit none

type (lat_struct), target :: lat, lat2, lat3
type (ele_struct), pointer :: ele, ele2
type (ele_struct) ele0
type (coord_struct), allocatable :: orbit_track(:)
type (coord_struct) start_orb, end_orb, end_orb1, end_orb2, end_orb1p, end_orb2p
type (coord_struct) end_orb1t, end_orb2t, closed_orb, orbit
type (normal_modes_struct) mode
type (taylor_struct) bmad_taylor(6), orb_taylor(6), spin_taylor(0:3), tlr
type (real_8) y8(6)
type (probe_8) prb
type (info), pointer :: info0
type (internal_state) ptc_state
type (branch_struct), pointer :: branch, branch2
type (track_struct) orb_track
type (em_field_struct) field
type (ptc_rad_map_struct) rad_map

real(rp) diff_mat(6,6), diff_vec(6), vec(6)
real(rp) vec_bmad(6), vec_ptc(6), vec_bmad2(6), beta0, beta1 
real(rp) m6_to_ptc(6,6), m6_to_bmad(6,6), m6(6,6), sigma_mat(6,6)
real(dp) b_pot, e_pot, b_field_ptc(3), e_field_ptc(3), a(3), da(3,2), x(6), z

integer i, j, k, seed0, ixd(3)
character(80) str
logical err_flag, survey_needed

namelist / params / start_orb

!

open (1, file = 'output.now')

!----------------------------------------------------------
! Check information passing between bmad element and associated ptc fibre

bmad_com%auto_bookkeeper = .false.

call bmad_parser ('diff_test.bmad', lat)
call lattice_bookkeeper (lat)
call lat_to_ptc_layout(lat)

branch => lat%branch(0)
do i = 1, branch%n_ele_track
  ele => branch%ele(i)
  ele0 = ele
  call vary_ele_attributes(ele)

  call update_fibre_from_ele (ele, survey_needed)
  call update_ele_from_fibre (ele0)
  str = 'NO-DIFF'
  call check_if_ele_different(ele0, ele, str)

  write (1, '(6a)') '"IN-OUT:', trim(ele%name), '" STR  "', trim(str), '"'
enddo

!----------------------------

call bmad_parser ('tiny_ring.bmad', lat)
call lat_to_ptc_layout (lat)

call ptc_setup_map_with_radiation (rad_map, lat%ele(0), lat%ele(0), map_order = 1)

call init_coord (start_orb, [1, 2, 3, 4, 5, 6] * 1e-3_rp, lat%ele(0), upstream_end$)
orbit = start_orb

call get_seed(seed0)
call init_coord (start_orb, rad_map%ref0, lat%ele(0), upstream_end$)
orbit = start_orb
call ptc_track_map_with_radiation (orbit, rad_map, .true., .true.)

call set_seed(seed0)
do i = 1, 6
  vec(i) = grnf_zhe()
enddo
vec = matmul(rad_map%stoc_mat, vec)
vec = vec + rad_map%ref1

write (1, '(a, 6es16.8)') '"Stoc-Track"      ABS 1E-13', orbit%vec
write (1, '(a, 6es16.8)') '"Diff-Stoc-Track" ABS 1E-20', orbit%vec - vec

call ptc_track_map_with_radiation (orbit, rad_map, .true., .false.)
end_orb%vec = start_orb%vec - rad_map%ref0
end_orb%vec = matmul(rad_map%nodamp_mat, end_orb%vec)
end_orb%vec = matmul(rad_map%damp_mat, end_orb%vec)
end_orb%vec = end_orb%vec + rad_map%ref1

write (1, '(a, 6f16.10)') '"Damp-Track"      ABS 2E-10', orbit%vec
write (1, '(a, 6f16.10)') '"Diff-Damp-Track" ABS 2E-10', orbit%vec - end_orb%vec


!----------------------------

!! call ptc_invariant_spin_field (lat%ele(0), 1, closed_orb)


!----------------------------

call bmad_parser ('figure_8.bmad', lat)
call lat_to_ptc_layout (lat)

bmad_com%spin_tracking_on = .true.
call ptc_one_turn_map_at_ele (lat%ele(2), orbit%vec, prb, ptc_state)

orb_taylor = prb%x
spin_taylor = prb%q%x
call taylor_clean(orb_taylor)
call taylor_clean(spin_taylor)

ixd = [1, 4, 5]  ! To cut down on output
do i = 1, 3
  j = ixd(i)
  call sort_taylor_terms(orb_taylor(ixd(i)), tlr)
  do k = 1, size(tlr%term)
    write (1, '(a, i0, a, i0, a, f20.12, 6i3)') '"OrbMap-', ixd(i), '-', k, '" ABS 1E-7', &
                                                                           tlr%term(k)%coef, tlr%term(k)%expn
  enddo
enddo

ixd = [1, 2, 0]  ! To cut down on output
do i = 1, 2
  j = ixd(i)
  call sort_taylor_terms(spin_taylor(ixd(i)), tlr)
  do k = 1, size(tlr%term)
    write (1, '(a, i0, a, i0, a, 2x, a, 6i3)') '"SpinMap-', ixd(i), '-', k, '" ABS 1E-7', &
                                                real_str(tlr%term(k)%coef,14), tlr%term(k)%expn
  enddo
enddo

call ptc_one_turn_mat_and_closed_orbit_calc (lat%branch(0))
info0 => lat%ele(0)%ptc_fibre%i
write (1, '(a, 6f16.10)') '"1Turn-Fix" ABS 1E-10', info0%fix
do i = 1, 6
  write (1, '(a, i0, a, 6es18.10)') '"1Turn-Mat', i, '" ABS 1E10', info0%m(i,:)
enddo

!----------------------------

call bmad_parser ('tiny_ring.bmad', lat)
call lat_to_ptc_layout (lat)

call ptc_emit_calc (lat%ele(0), mode, sigma_mat, closed_orb)

write (1, '(a, 3es16.8)') '"layout-tune" REL 1E-8', mode%a%tune, mode%b%tune, mode%z%tune
write (1, '(a, 3es16.8)') '"layout-emit" REL 1E-3', mode%a%emittance, mode%b%emittance, mode%z%emittance
write (1, '(a, 3es16.8)') '"layout-damp" REL 3E-5', mode%a%alpha_damp, mode%b%alpha_damp, mode%z%alpha_damp
write (1, '(a, 6es16.8)') '"layout-orb"  ABS 1E-14', closed_orb%vec
do i = 1, 6
  write (1, '(a, i0, a, 6es16.8)') '"layout-sigma', i, '" ABS 2E-12', sigma_mat(i,:)
enddo

bmad_com%radiation_damping_on = .false.

!-------------------------------------------------------

call bmad_parser ('single_quad.bmad', lat)

ele => lat%ele(1)
call ele_to_fibre(ele, ele%ptc_fibre, .true., err_flag)
x(1) = lat%particle_start%vec(1)
x(3) = lat%particle_start%vec(3)
z = lat%particle_start%vec(5)
call b_e_field(ele%ptc_fibre%mag%ab, x, z, e_pot, e_field_ptc, b_pot, b_field_ptc, a, da)
call em_field_calc (ele, lat%param, z, lat%particle_start, .true., field)

write (1, '(a, 3f14.6)') '"B PTC:" ABS 1E-20  ', b_field_ptc * ele%value(p0c$) / c_light
write (1, '(a, 3f14.6)') '"B Bmad:" ABS 1E-20 ', field%b
write (1, '(a, 3f14.6)') '"B Diff:" ABS 1E-20 ', b_field_ptc * ele%value(p0c$) / c_light - field%b

write (1, '(a, 3f14.6)') '"E PTC:" ABS 1E-20  ', e_field_ptc * ele%value(p0c$)
write (1, '(a, 3f14.6)') '"E Bmad:" ABS 1E-20 ', field%e
write (1, '(a, 3f14.6)') '"E Diff:" ABS 1E-20 ', e_field_ptc * ele%value(p0c$) - field%e

!----------------------------------------------------------
! Check information passing via flat file.

call bmad_parser ('ptc_test.bmad', lat)
!!
!!call kill_ptc_layouts(lat)
!!call branch_to_ptc_m_u(branch, .false.)
!!call write_ptc_flat_file_lattice ('ptc_test.flat', branch)
!!call ptc_read_flat_file (['ptc_test.flat'], err_flag, lat2)
!!do i = 1, branch%n_ele_track
!!  call check_if_diffrent
!!

!------------------------------------------------------------------------
! Tracking tests
! ele(1) is a wiggler, ele(2:4) is the same wiggler split by a marker.

open (2, file = 'ptc_test.bmad')
read (2, nml = params)
close (2)

lat%ele(3)%ix_ele = -1
call remove_eles_from_lat (lat)

call init_coord(start_orb, lat%particle_start, lat%ele(0), downstream_end$)
call track1 (start_orb, lat%ele(1), lat%param, end_orb1)
call track1 (start_orb, lat%ele(2), lat%param, end_orb2)
call track1 (end_orb2, lat%ele(3), lat%param, end_orb2)
call lat_make_mat6(lat, -1)
lat%ele(3)%vec0 = lat%ele(3)%vec0 + matmul(lat%ele(3)%mat6, lat%ele(2)%vec0)
lat%ele(3)%mat6 = matmul(lat%ele(3)%mat6, lat%ele(2)%mat6)

diff_mat = lat%ele(3)%mat6 - lat%ele(1)%mat6
diff_vec = lat%ele(3)%vec0 - lat%ele(1)%vec0

write (1, *)
write (1, '(a, es20.10)') '"Bmad:vec(1)" ABS  2E-16', end_orb1%vec(1)
write (1, '(a, es20.10)') '"Bmad:vec(2)" ABS  2E-16', end_orb1%vec(2)
write (1, '(a, es20.10)') '"Bmad:vec(3)" ABS  2E-16', end_orb1%vec(3)
write (1, '(a, es20.10)') '"Bmad:vec(4)" ABS  2E-16', end_orb1%vec(4)
write (1, '(a, es20.10)') '"Bmad:vec(5)" ABS  4E-16', end_orb1%vec(5)
write (1, '(a, es20.10)') '"Bmad:vec(6)" ABS  2E-16', end_orb1%vec(6)
write (1, '(a, es20.10)') '"Bmad:orb%t " ABS  2E-16', end_orb1%t

write (1, *)
write (1, '(a, 6es10.2)') '"Bmad2:dvec"   ABS  1E-11', end_orb1%vec - end_orb2%vec
write (1, '(a, es10.2)')  '"Bmad2:dmat"   ABS  1E-11', maxval(abs(diff_mat))
write (1, '(a, es10.2)')  '"Bmad2:dvec0"  ABS  1E-11', maxval(abs(diff_vec))
write (1, '(a, es10.2)')  '"Bmad2:dorb%t" ABS  1E-22', &
          (end_orb1%t - lat%ele(1)%ref_time) - (end_orb2%t - lat%ele(3)%ref_time)

!

lat2 = lat
lat2%ele%tracking_method = symp_lie_ptc$
lat2%ele%mat6_calc_method = symp_lie_ptc$

call track1 (start_orb, lat2%ele(1), lat2%param, end_orb1p)
call track1 (start_orb, lat2%ele(2), lat2%param, end_orb2p)
call track1 (end_orb2p, lat2%ele(3), lat2%param, end_orb2p)
call lat_make_mat6(lat2, -1)
lat2%ele(3)%vec0 = lat2%ele(3)%vec0 + matmul(lat2%ele(3)%mat6, lat2%ele(2)%vec0)
lat2%ele(3)%mat6 = matmul(lat2%ele(3)%mat6, lat2%ele(2)%mat6)

diff_mat = lat2%ele(3)%mat6 - lat%ele(1)%mat6
diff_vec = lat2%ele(3)%vec0 - lat%ele(1)%vec0

write (1, *)
write (1, '(a, 6es10.2)') '"PTC1:dvec"   ABS  5E-09', end_orb1%vec - end_orb1p%vec
write (1, '(a, es10.2)')  '"PTC1:dorb%t" ABS  1E-20', end_orb1%t - end_orb1p%t
write (1, '(a, 6es10.2)') '"PTC2:dvec"   ABS  5E-09', end_orb1%vec - end_orb2p%vec
write (1, '(a, es10.2)')  '"PTC2:dmat"   ABS  2E-06', maxval(abs(diff_mat))
write (1, '(a, es10.2)')  '"PTC2:dvec0"  ABS  2E-11', maxval(abs(diff_vec))
write (1, '(a, es10.2)')  '"PTC2:dorb%t" ABS  1E-22', &
          (end_orb1%t - lat%ele(1)%ref_time) - (end_orb2%t - lat%ele(3)%ref_time)

!

lat3 = lat2
lat3%ele%tracking_method = taylor$
lat3%ele%mat6_calc_method = taylor$

call track1 (start_orb, lat3%ele(1), lat3%param, end_orb1t)
call track1 (start_orb, lat3%ele(2), lat3%param, end_orb2t)
call track1 (end_orb2t, lat3%ele(3), lat3%param, end_orb2t)
call lat_make_mat6(lat3, -1)
lat3%ele(3)%vec0 = lat3%ele(3)%vec0 + matmul(lat3%ele(3)%mat6, lat3%ele(2)%vec0)
lat3%ele(3)%mat6 = matmul(lat3%ele(3)%mat6, lat3%ele(2)%mat6)

diff_mat = lat3%ele(3)%mat6 - lat2%ele(1)%mat6
diff_vec = lat3%ele(3)%vec0 - lat2%ele(1)%vec0

write (1, *)
write (1, '(a, 6es10.2)') '"TAYLOR1:dvec"   ABS  1E-9', end_orb1p%vec - end_orb1t%vec
write (1, '(a, es10.2)')  '"TAYLOR1:dorb%t" ABS  1E-17', end_orb1p%t - end_orb1t%t
write (1, '(a, 6es10.2)') '"TAYLOR2:dvec"   ABS  1E-09', end_orb1p%vec - end_orb2t%vec
write (1, '(a, es10.2)')  '"TAYLOR2:dmat"   ABS  1E-11', maxval(abs(diff_mat))
write (1, '(a, es10.2)')  '"TAYLOR2:dvec0"  ABS  1E-11', maxval(abs(diff_vec))
write (1, '(a, es10.2)')  '"TAYLOR2:dorb%t" ABS  1E-22', &
          (end_orb1p%t - lat2%ele(1)%ref_time) - (end_orb2p%t - lat2%ele(3)%ref_time)

! Map translation

call alloc(y8)
y8 = lat3%ele(1)%taylor
bmad_taylor%ref = lat3%ele(1)%taylor%ref
bmad_taylor = y8
bmad_taylor = bmad_taylor - lat3%ele(1)%taylor

do i = 1, 6
  diff_vec = maxval(abs(bmad_taylor(i)%term(:)%coef))
enddo

!write (1, '(a, 6es10.2)') '"map_convert" ABS 1E-15', diff_vec

!-----------------------------------------------------------------
contains

subroutine vary_ele_attributes(ele)

type (ele_struct) ele
type (ele_attribute_struct) attrib
character(40) name
character(40), allocatable :: name_list(:)

!

do j = 1, num_ele_attrib$
  attrib = attribute_info(ele, j)
  if (attrib%state == does_not_exist$ .or. attrib%state == private$) cycle
  select case (attrib%kind)
  case (is_logical$)
    if (ele%value(j) == 0) then
      ele%value(j) = 1
    else
      ele%value(j) = 0
    endif
  case (is_integer$)
    ele%value(j) = 2 * ele%value(j) + 1
  case (is_real$)
    ele%value(j) = 1.4 * ele%value(j) + 1
  case (is_switch$)
    if (nint(ele%value(j)) < 2) then
      ele%value(j) = nint(ele%value(j)) + 1
    else
      ele%value(j) = nint(ele%value(j)) - 1
    endif
  end select

  call set_flags_for_changed_attribute(ele, ele%value(j))
enddo

if (has_attribute(ele, 'A0') .or. ele%key == multipole$) then
  call multipole_init(ele, all$)
  ele%b_pole(4) = ele%b_pole(3) + 1
endif

call attribute_bookkeeper (ele, .true.)

end subroutine vary_ele_attributes

!-----------------------------------------------------------------
! contains

subroutine check_if_ele_different(ele, ele2, str)

type (ele_struct) ele, ele2
type (ele_attribute_struct) attrib
real(rp) a_pole, b_pole, a_pole2, b_pole2
integer j
character(*) str

!

do j = 1, num_ele_attrib$
  attrib = attribute_info(ele, j)
  if (attrib%name == null_name$) cycle
  if (attrib%state == private$) cycle
  call check_if_value_different (str, ele, j, ele%value(j), ele2%value(j))
enddo

do j = 0, n_pole_maxx
  a_pole = 0; b_pole = 0
  a_pole2 = 0; b_pole2 = 0

  if (associated(ele%a_pole)) then 
    a_pole = ele%a_pole(j); b_pole = ele%a_pole(j)
  endif
  if (associated(ele2%a_pole)) then
    a_pole2 = ele2%a_pole(j); b_pole2 = ele2%a_pole(j)
  endif

  call check_if_value_different (str, ele, j+a0$, a_pole, a_pole2)
  call check_if_value_different (str, ele, j+b0$, b_pole, b_pole2)
enddo

end subroutine check_if_ele_different

!-----------------------------------------------------------------
! contains

subroutine check_if_value_different (str, ele, ix_attrib, val, val2)

type (ele_struct) ele
real(rp) val, val2
integer ix_attrib
character(*) str

! drift like elements in ptc only use one integration step independent
! of Bmad so ignore any of these differences.

select case (ele%key)
case (drift$, rcollimator$, ecollimator$, monitor$, instrument$, pipe$) 
  if (ix_attrib == ds_step$ .or. ix_attrib == num_steps$) return
end select

if (val == 0 .and. val2 == 0) return

if (val * val2 <= 0) then  ! opposite sign
  if (abs(val) + abs(val2) < 2e-10) return
else
  if (abs(val-val2) / abs(val+val2) < 1e-10) return
endif

write (str, '(a, 2es20.10)') trim(attribute_name(ele, ix_attrib)), val, val2

end subroutine check_if_value_different

end program
