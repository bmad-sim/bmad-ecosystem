!+
! Program ptc_test
!
! This program is part of the Bmad regression testing suite.
!-

program ptc_test

use bmad
use polymorphic_taylor
use ptc_layout_mod

implicit none

type (lat_struct), target :: lat, lat2, lat3
type (ele_struct), pointer :: ele, ele2
type (coord_struct) start_orb, end_orb1, end_orb2, end_orb1p, end_orb2p
type (coord_struct) end_orb1t, end_orb2t, closed_orb
type (normal_modes_struct) mode
type (taylor_struct) bmad_taylor(6)
type (real_8) y8(6)
type (branch_struct), pointer :: branch, branch2
type (track_struct) track

real(rp) diff_mat(6,6), diff_vec(6)
real(rp) vec_bmad(6), vec_ptc(6), vec_bmad2(6), beta0, beta1 
real(rp) m6_to_ptc(6,6), m6_to_bmad(6,6), m6(6,6), sigma_mat(6,6)
real(rp) a_pole, b_pole, a_pole2, b_pole2

integer i, j
character(80) str
logical err_flag

namelist / params / start_orb

!

open (1, file = 'output.now')

!-------------------------------------------------------

call bmad_parser ('track.bmad', lat)

call track1 (lat%beam_start, lat%ele(1), lat%param, end_orb1, track, err_flag)

print *

do i = 0, track%n_pt
  print '(i4, f10.6, 4x, 6f10.6)', i, track%orb(i)%s, track%orb(i)%vec
enddo

stop

!----------------------------------------------------------
! Check information passing between bmad element and associated ptc fibre
! Procedure: 
!   0) branch%ele(i) and branch2%ele(i) are the same type of element but
!      branch2 elements have zero parameter values.
!   1) Setup ptc layout using branch2.
!   2) Point branch%ele to same fibre as branch%ele
!   3) Transfer attribute values from branch%ele -> ptc fibre -> branch2%ele
!   4) Check that branch%ele and branch2%ele have same attribute values.

bmad_com%use_hard_edge_drifts = .false.

call bmad_parser ('ptc_test.bmad', lat)

branch => lat%branch(1)
branch2 => lat%branch(2)

call branch_to_ptc_m_u (branch2, .false.)

do i = 1, branch%n_ele_track
  ele => branch%ele(i)
  ele2 => branch2%ele(i)
  ele%ptc_fibre => ele2%ptc_fibre
  call update_ptc_fibre_from_bmad (ele)
  call update_bmad_ele_from_ptc (ele2)
enddo

call lattice_bookkeeper (lat)

do i = 1, branch%n_ele_track
  ele => branch%ele(i)
  ele2 => branch2%ele(i)
  str = 'NO-DIFF'
  call check_if_ele_different(ele, ele2)

  write (1, '(6a)') '"IN-OUT:', trim(ele%name), '" STR  "', trim(str), '"'
enddo

!----------------------------------------------------------
! Check information passing via flat file.

call kill_ptc_layouts(lat)
call branch_to_ptc_m_u(branch, .false.)
call write_ptc_flat_file_lattice ('ptc_test.flat', branch)
call ptc_read_flat_file (['ptc_test.flat'], err_flag, lat2)
!!do i = 1, branch%n_ele_track
!!  call check_if_diffrent
!!

!------------------------------------------------------------------------
! Tracking tests
! ele(1) is a wiggler, ele(2:4) is the same wiggler split by a marker.

open (2, file = 'ptc_test.bmad')
read (2, nml = params)
close (2)

lat%ele(3)%key = -1
call remove_eles_from_lat (lat)

call track1 (start_orb, lat%ele(1), lat%param, end_orb1)
call track1 (start_orb, lat%ele(2), lat%param, end_orb2)
call track1 (end_orb2, lat%ele(3), lat%param, end_orb2)
call lat_make_mat6(lat, -1)
lat%ele(3)%vec0 = lat%ele(3)%vec0 + matmul(lat%ele(3)%mat6, lat%ele(2)%vec0)
lat%ele(3)%mat6 = matmul(lat%ele(3)%mat6, lat%ele(2)%mat6)

diff_mat = lat%ele(3)%mat6 - lat%ele(1)%mat6
diff_vec = lat%ele(3)%vec0 - lat%ele(1)%vec0

write (1, *)
write (1, '(a, es20.10)') '"Bmad:vec(1)" REL  1E-10', end_orb1%vec(1)
write (1, '(a, es20.10)') '"Bmad:vec(2)" REL  1E-10', end_orb1%vec(2)
write (1, '(a, es20.10)') '"Bmad:vec(3)" REL  1E-10', end_orb1%vec(3)
write (1, '(a, es20.10)') '"Bmad:vec(4)" REL  1E-10', end_orb1%vec(4)
write (1, '(a, es20.10)') '"Bmad:vec(5)" REL  1E-10', end_orb1%vec(5)
write (1, '(a, es20.10)') '"Bmad:vec(6)" REL  1E-10', end_orb1%vec(6)
write (1, '(a, es20.10)') '"Bmad:orb%t " REL  1E-10', end_orb1%t

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

! Vector translation

vec_bmad = [0.01, 0.02, 0.03, 0.04, 1.0, 2.0]
beta0 = 0.6
beta1 = 0.7

call vec_bmad_to_ptc (vec_bmad, beta0, vec_ptc, m6_to_ptc)
call vec_ptc_to_bmad (vec_ptc, beta0, vec_bmad2, m6_to_bmad)

m6 = matmul(m6_to_ptc, m6_to_bmad)
forall (j = 1:6) m6(j,j) = m6(j,j) - 1

write (1, *) 
write (1, '(a, 6es10.2)') '"vec_convert" ABS 1E-15', vec_bmad2 - vec_bmad
write (1, '(a, es10.2)')  '"mat_convert" ABS 1E-15', maxval(abs(m6))

! Map translation

call real_8_init(y8)
call taylor_to_real_8 (lat3%ele(1)%taylor, beta0, beta1, y8)
bmad_taylor%ref = lat3%ele(1)%taylor%ref
call real_8_to_taylor (y8, beta0, beta1, bmad_taylor)
bmad_taylor = bmad_taylor - lat3%ele(1)%taylor

do i = 1, 6
  diff_vec = maxval(abs(bmad_taylor(i)%term(:)%coef))
enddo

write (1, '(a, 6es10.2)') '"map_convert" ABS 1E-15', diff_vec

!----------------------------
! Layout test

call bmad_parser ('figure_8.bmad', lat)
call lat_to_ptc_layout (lat)
call ptc_emit_calc (lat%ele(0), mode, sigma_mat, closed_orb)

write (1, '(a, 3es16.8)') '"layout-tune" REL 1E-8', mode%a%tune, mode%b%tune, mode%z%tune
write (1, '(a, 3es16.8)') '"layout-emit" REL 1E-8', mode%a%emittance, mode%b%emittance, mode%z%emittance
write (1, '(a, 3es16.8)') '"layout-damp" REL 1E-8', mode%a%alpha_damp, mode%b%alpha_damp, mode%z%alpha_damp
write (1, '(a, 6es16.8)') '"layout-orb"  REL 1E-8', closed_orb%vec
do i = 1, 6
  write (1, '(a, i0, a, 6es16.8)') '"layout-sigma', i, '" REL 1E-8', sigma_mat(i,:)
enddo


!-----------------------------------------------------------------
contains

subroutine check_if_ele_different(ele, ele2)

type (ele_struct) ele, ele2

!!!  do j = 1, num_ele_attrib$
!!!    if (j == ds_step$) cycle
!!!    if (j == num_steps$) cycle
!!!    if (j == integrator_order$) cycle
!!!    call check_if_different (str, ele, j, ele%value(j), ele2%value(j))
!!!  enddo
!!!  do j = 0, n_pole_maxx
!!!    if (associated(ele%a_pole)) then
!!!      a_pole = ele%a_pole(j); b_pole = ele%a_pole(j)
!!!    else
!!!      a_pole = 0; b_pole = 0
!!!    endif
!!!
!!!    if (associated(ele2%a_pole)) then
!!!      a_pole2 = ele2%a_pole(j); b_pole2 = ele2%a_pole(j)
!!!    else
!!!      a_pole2 = 0; b_pole2 = 0
!!!    endif
!!!
!!!    call check_if_different (str, ele, j+a0$, a_pole, a_pole2)
!!!    call check_if_different (str, ele, j+b0$, b_pole, b_pole2)
!!!  enddo

end subroutine check_if_ele_different

!-----------------------------------------------------------------
! contains

subroutine check_if_different (str, ele, ix_attrib, val, val2)

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

end subroutine

end program
