!+
! Program test_bmad_to_c
!
! Program to test the conversion routines between bmad and C/C++.
! The general idea is for each structure type:
!
! 1) From the F main program, pass a structure to a C++ routine which 
! translates the structure to C++ and checks that all the components have the
! correct value. This checks that a C++ routine can do the F -> C++ conversion.
!
! 2) The C++ routine then passes the C++ structure to a F routine. The F
! routine translates the structure to F and checks that all the components
! have the correct value. This checks that a F routine can do the C++ -> F conversion.
! 
! 3) The F routine converts a structure from F to C++ and sends it back to the calling
! C++ routine that checks that the conversion is good. This checks that a F routine
! can do the F -> C++ conversion.
!
! 4) The C++ routine converts the C++ structure to F and sends it back to the
! main program which checks it. This checks that a C++ routine can do the 
! C++ -> F conversion.
!-

module test_mod

use bmad
use spin_mod

type (coord_struct), save ::           coord_in, coord_out
type (twiss_struct), save ::           twiss_in, twiss_out
type (xy_disp_struct), save ::         xy_disp_in, xy_disp_out
type (floor_position_struct), save ::  floor_position_in, floor_position_out
type (wig_term_struct), save ::        wig_term_in, wig_term_out
type (taylor_term_struct), save ::     taylor_term_in, taylor_term_out
type (taylor_struct), save ::          taylor_in, taylor_out
type (sr_table_wake_struct), save ::   sr_table_wake_in, sr_table_wake_out
type (sr_mode_wake_struct), save ::    sr_mode_wake_in, sr_mode_wake_out
type (lr_wake_struct), save ::         lr_wake_in, lr_wake_out
type (wake_struct), save, target ::    wake_in, wake_out
type (control_struct), save ::         control_in, control_out
type (lat_param_struct), save ::       param_in, param_out
type (anormal_mode_struct), save ::    amode_in, amode_out
type (linac_normal_mode_struct), save :: linac_mode_in, linac_mode_out
type (normal_modes_struct), save ::    modes_in, modes_out
type (bmad_common_struct), save ::     bmad_com_in, bmad_com_out
type (em_field_struct), save ::        em_field_in, em_field_out
type (ele_struct), save ::             ele_in, ele_out
type (mode_info_struct), save ::       mode_info_in, mode_info_out
type (lat_struct), save ::             lat_in, lat_out

logical, save :: all_ok = .true.

!-----------------------

contains

subroutine init_all_structs

implicit none

real(rp) mat6_a(6,6), mat6_b(6,6), mat3_a(3,3), mat3_b(3,3), mat3_c(3,3)
real(rp) vec6_a(6), vec6_b(6), vec3_a(3), vec3_b(3), vec3_c(3)
real(rp) mat2_a(2,2)
complex(rp) spinor2_a(2), spinor2_b(2)
integer i, j
logical :: T = .true., F = .false.

!

forall(i = 1:6, j = 1:6) 
  mat6_a(i, j) = 10*(i-1) + j 
  mat6_b(i, j) = 10*(i-1) + j  + 100
end forall

forall(i = 1:3, j = 1:3) 
  mat3_a(i, j) = 10*(i-1) + j 
  mat3_b(i, j) = 10*(i-1) + j  + 100
  mat3_c(i, j) = 10*(i-1) + j  + 200
end forall

forall(i = 1:2, j = 1:2) 
  mat2_a(i, j) = 10*(i-1) + j 
end forall

forall (i = 1:6)
  vec6_a(i) = i + 100
  vec6_b(i) = i + 200
end forall

forall (i = 1:2)
  spinor2_a(i) = i + 100
  spinor2_b(i) = i + 200
end forall

forall (i = 1:3)
  vec3_a(i) = i + 100
  vec3_b(i) = i + 200
  vec3_c(i) = i + 300
end forall

!

coord_in  = coord_struct(vec6_a, spinor2_a, 1.0_rp)
coord_out = coord_struct(vec6_b, spinor2_b, 2.0_rp)

twiss_in  = twiss_struct(1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp, 5.0_rp, 6.0_rp, &
                                      7.0_rp, 8.0_rp, 9.0_rp, 10.0_rp)
twiss_out = twiss_struct(9.0_rp, 8.0_rp, 7.0_rp, 6.0_rp, 5.0_rp, &
                                      4.0_rp, 3.0_rp, 2.0_rp, 1.0_rp, 0.0_rp)

xy_disp_in  = xy_disp_struct(1.0_rp, 2.0_rp)
xy_disp_out = xy_disp_struct(2.0_rp, 1.0_rp)

floor_position_in  = floor_position_struct(1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp, 5.0_rp, 6.0_rp)
floor_position_out = floor_position_struct(6.0_rp, 5.0_rp, 4.0_rp, 3.0_rp, 2.0_rp, 1.0_rp)

wig_term_in  = wig_term_struct(1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp, 5.0_rp, 6)
wig_term_out = wig_term_struct(6.0_rp, 5.0_rp, 4.0_rp, 3.0_rp, 2.0_rp, 1)

taylor_term_in  = taylor_term_struct(1.0_rp, (/ 2, 3, 4, 5, 6, 7 /))
taylor_term_out = taylor_term_struct(7.0_rp, (/ 6, 5, 4, 3, 2, 1/))

allocate(taylor_in%term(2), taylor_out%term(2))

taylor_in%ref = -1
taylor_in%term(1) = taylor_term_struct(1.0_rp, (/ 2, 3, 4, 5, 6, 7 /))
taylor_in%term(2) = taylor_term_struct(8.0_rp, (/ 9, 10, 11, 12, 13, 14 /)) 

taylor_out%ref = 1.0_rp
taylor_out%term(1) = taylor_term_struct(-1.0_rp, (/ -2, -3, -4, -5, -6, -7 /))
taylor_out%term(2) = taylor_term_struct(-8.0_rp, (/ -9, -10, -11, -12, -13, -14 /)) 

sr_table_wake_in  = sr_table_wake_struct (1.0_rp, 2.0_rp, 3.0_rp)
sr_table_wake_out = sr_table_wake_struct (3.0_rp, 2.0_rp, 1.0_rp)

sr_mode_wake_in  = sr_mode_wake_struct (21.0_rp, 22.0_rp, 23.0_rp, 24.0_rp, 25.0_rp, &
                                                     26.0_rp, 27.0_rp, 28.0_rp)
sr_mode_wake_out = sr_mode_wake_struct (31.0_rp, 32.0_rp, 33.0_rp, 34.0_rp, 35.0_rp, &
                                                     36.0_rp, 37.0_rp, 38.0_rp)

lr_wake_in  = lr_wake_struct(1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp, 5.0_rp, &
                             6.0_rp, 7.0_rp, 8.0_rp, 9.0_rp, 10.0_rp, 11, .true.)
lr_wake_out = lr_wake_struct(11.0_rp, 10.0_rp, 9.0_rp, 8.0_rp, 7.0_rp, 6, 5.0_rp, &
                              4.0_rp, 3.0_rp, 2.0_rp, 1, .false.)

allocate (wake_in%sr_table(0:1), wake_in%sr_mode_long(2), wake_in%sr_mode_trans(0), wake_in%lr(1))
allocate (wake_out%sr_table(0), wake_out%sr_mode_long(0), wake_out%sr_mode_trans(2), wake_out%lr(2))

wake_in%sr_file = "ABCD"
wake_in%lr_file = "XYZZY"
wake_in%z_sr_mode_max = 100
wake_in%sr_table(0) = sr_table_wake_in
wake_in%sr_table(1) = sr_table_wake_out
wake_in%sr_mode_long(1) = sr_mode_wake_in
wake_in%sr_mode_long(2) = sr_mode_wake_out
wake_in%lr(1) = lr_wake_in

wake_out%sr_file = "abcd"
wake_out%lr_file = "xyzzy"
wake_out%z_sr_mode_max = 101
wake_out%sr_mode_trans(1) = sr_mode_wake_in
wake_out%sr_mode_trans(2) = sr_mode_wake_out
wake_out%lr(1) = lr_wake_struct (-1.0_rp, -2.0_rp, -3.0_rp, -4.0_rp, -5.0_rp, &
                          -6.0_rp, -7.0_rp, -8.0_rp, -9.0_rp, -10.0_rp, 11, .true.)
wake_out%lr(2) = lr_wake_struct (-11.0_rp, -12.0_rp, -13.0_rp, -14.0_rp, -15.0_rp, &
                                -16.0_rp, -17.0_rp, -18.0_rp, -19.0_rp, -20.0_rp, 21, .false.)

control_in  = control_struct(1.0_rp, 2, 3, 4, 5)
control_out = control_struct(5.0_rp, 4, 3, 2, 1)

param_in  = lat_param_struct (1.0_rp, 2.0_rp, &
      3.0_rp, mat6_a, mat6_b, 11, 12, 13, 14, 15, 16, T, F, T)
param_out = lat_param_struct (11.0_rp, 12.0_rp, &
      13.0_rp, mat6_b, mat6_a, 111, 112, 113, 114, 115, 116, F, F, F)

amode_in  = anormal_mode_struct (1.0_rp, (/2.0_rp, 3.0_rp/), 4.0_rp, 5.0_rp, 6.0_rp, 7.0_rp)
amode_out  = anormal_mode_struct (11.0_rp, (/12.0_rp, 13.0_rp/), 14.0_rp, 15.0_rp, 16.0_rp, 17.0_rp)

linac_mode_in  = linac_normal_mode_struct (1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp, 5.0_rp, 6.0_rp, 7.0_rp)
linac_mode_out = linac_normal_mode_struct (11.0_rp, 12.0_rp, 13.0_rp, 14.0_rp, 15.0_rp, 16.0_rp, 17.0_rp)

modes_in  = normal_modes_struct ((/-1.0_rp, 1.0_rp, 2.0_rp, 3.0_rp/), 4.0_rp, 5.0_rp, 6.0_rp, 7.0_rp, &
                  amode_in, amode_out, amode_in, linac_mode_in)
modes_out = normal_modes_struct ((/20.0_rp, 11.0_rp, 12.0_rp, 13.0_rp/), 14.0_rp, 15.0_rp, 16.0_rp, 17.0_rp, &
                  amode_out, amode_out, amode_in, linac_mode_out)

bmad_com_in  = bmad_common_struct (2.0_rp, vec6_a, 3_rp, 4.0_rp, 5.0_rp, 6.0_rp, 7.0_rp, 8.0_rp, 9.0_rp, &
                                10, 11, T, F, T, F, T, F, T, F, T, T, F, T)
bmad_com_out = bmad_common_struct(12.0_rp, vec6_b, 13_rp, 14.0_rp, 15.0_rp, 16.0_rp, 17.0_rp, &
                                18.0_rp, 19.0_rp, 20, 21, T, T, F, F, T, T, F, F, T, F, T, F)

em_field_in  = em_field_struct (vec3_a, vec3_b, vec3_c, mat3_a, mat3_b, mat3_c, 77)
em_field_out = em_field_struct (vec3_c, vec3_b, vec3_a, mat3_c, mat3_b, mat3_a, -77)

mode_info_in   = mode_info_struct (1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp, 5.0_rp)
mode_info_out  = mode_info_struct (-1.0_rp, -2.0_rp, -3.0_rp, -4.0_rp, -5.0_rp)

call init_ele (ele_in)
call init_ele (ele_out)

ele_in%name = "1234"
ele_in%type = "abcd"
ele_in%alias = "5678"
ele_in%component_name = "efgh" 
ele_in%value(1) = -1
ele_in%value(n_attrib_maxx) = -n_attrib_maxx
ele_in%gen0 = vec6_a
ele_in%vec0 = vec6_b
ele_in%mat6 = mat6_a
ele_in%c_mat = mat2_a
ele_in%gamma_c = 51
ele_in%s = 52
ele_in%ref_time = 53
ele_in%key = 54
ele_in%sub_key = 55
ele_in%lord_status = 561
ele_in%slave_status = 562
ele_in%ix_value = 57
ele_in%n_slave = 58
ele_in%ix1_slave = 59
ele_in%ix2_slave = 60
ele_in%n_lord = 61
ele_in%ic1_lord = 62
ele_in%ic2_lord = 63
ele_in%ix_pointer = 64
ele_in%ixx = 65
ele_in%ix_ele = 66 
ele_in%mat6_calc_method = 67
ele_in%tracking_method = 68
ele_in%field_calc = 69
ele_in%ref_orbit = 72
ele_in%taylor_order = 73
ele_in%aperture_at = 74
ele_in%aperture_type = 75
ele_in%symplectify = T
ele_in%mode_flip = F
ele_in%multipoles_on = T
ele_in%map_with_offsets = F
ele_in%field_master = T
ele_in%is_on = F
ele_in%old_is_on = T
ele_in%logic = F
ele_in%on_a_girder = T
ele_in%csr_calc_on = F
ele_in%offset_moves_aperture = T

ele_in%floor = floor_position_in
ele_in%a = twiss_in
ele_in%b = twiss_out
ele_in%z = twiss_in
ele_in%x = xy_disp_in
ele_in%y = xy_disp_out

ele_out = ele_in
ele_out%ix_ele = ele_in%ix_ele
deallocate(ele_out%r)

call multipole_init (ele_in)
allocate (ele_in%r(2,3), ele_in%const(6), ele_in%descrip, ele_in%wig_term(2), ele_in%gen_field)

ele_in%descrip = "descrip"
ele_in%const = -vec6_b
ele_in%wake => wake_in
ele_in%wig_term(1) = wig_term_in
ele_in%wig_term(2) = wig_term_out
ele_in%taylor(1) = taylor_in
ele_in%taylor(2) = taylor_out
ele_in%r = mat3_c(1:2,1:3)
forall (i = 0:n_pole_maxx)
  ele_in%a_pole(i) = -i - 100
  ele_in%b_pole(i) = -i - 200
end forall

lat_in%name = "abc"
lat_in%lattice = "def"
lat_in%input_file_name = "123"
lat_in%title = "title"
lat_in%version = 21
lat_in%n_ele_track = 22
lat_in%n_ele_max = 1
lat_in%n_control_max = 2
lat_in%n_ic_max = 6
lat_in%input_taylor_order = 34

allocate (lat_in%ele(0:1))
lat_in%ele(0) = ele_in
lat_in%ele(1) = ele_out
allocate (lat_in%control(lat_in%n_control_max))
lat_in%control(1) = control_in
lat_in%control(2) = control_out
allocate(lat_in%ic(6))
lat_in%ic = nint(vec6_a)
lat_in%ele_init = ele_in
lat_in%ele_init%name = 'ele_init'
lat_in%ele(0)%ix_ele = 0
lat_in%ele(1)%ix_ele = 1

lat_out = lat_in

end subroutine

end module

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

program test_fortran_and_cpp

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (coord_struct)          coord1, coord2
type (twiss_struct)          twiss1, twiss2
type (xy_disp_struct)        xy_disp1, xy_disp2
type (floor_position_struct) floor_position1, floor_position2
type (wig_term_struct)       wig_term1, wig_term2
type (taylor_term_struct)    taylor_term1, taylor_term2
type (taylor_struct)         taylor1, taylor2
type (sr_table_wake_struct)  sr_table_wake1, sr_table_wake2
type (sr_mode_wake_struct)   sr_mode_wake1, sr_mode_wake2
type (lr_wake_struct)        lr_wake1, lr_wake2
type (wake_struct)           wake1, wake2
type (control_struct)        control1, control2
type (lat_param_struct)      param1, param2
type (anormal_mode_struct)   amode1, amode2
type (linac_normal_mode_struct)     linac_mode1, linac_mode2
type (normal_modes_struct)   modes1, modes2
type (em_field_struct)       em_field1, em_field2
type (ele_struct)            ele1, ele2
type (mode_info_struct)      mode_info1, mode_info2
type (lat_struct)           lat1, lat2

integer :: c_ok = 1

!------------------------------------------------
! init

call init_all_structs

!------------------------------------------------
! Coord Check

coord1 = coord_in

call test_c_coord (coord1, coord2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (coord2 == coord_out) then
  print *, 'C_side_convert: Coord C to F: OK'
else
  print *, 'C_SIDE_CONVERT: COORD C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! Twiss Check

twiss1 = twiss_in

print *
call test_c_twiss (twiss1, twiss2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (twiss2 == twiss_out) then
  print *, 'C_side_convert: Twiss C to F: OK'
else
  print *, 'C_SIDE_CONVERT: TWISS C TO F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! Xy_disp Check

xy_disp1 = xy_disp_in

print *
call test_c_xy_disp (xy_disp1, xy_disp2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (xy_disp2 == xy_disp_out) then
  print *, 'C_side_convert: Xy_disp C to F: OK'
else
  print *, 'C_SIDE_CONVERT: XY_DISP C TO F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! floor_position Check

floor_position1 = floor_position_in

print *
call test_c_floor_position (floor_position1, floor_position2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (floor_position2 == floor_position_out) then
  print *, 'C_side_convert: floor_position C to F: OK'
else
  print *, 'C_SIDE_CONVERT: floor_position C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! wig_term Check

wig_term1 = wig_term_in

print *
call test_c_wig_term (wig_term1, wig_term2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (wig_term2 == wig_term_out) then
  print *, 'C_side_convert: wig_term C to F: OK'
else
  print *, 'C_SIDE_CONVERT: wig_term C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! taylor_term Check

taylor_term1 = taylor_term_in

print *
call test_c_taylor_term (taylor_term1, taylor_term2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (taylor_term2 == taylor_term_out) then
  print *, 'C_side_convert: taylor_term C to F: OK'
else
  print *, 'C_SIDE_CONVERT: taylor_term C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! taylor Check

taylor1 = taylor_in

print *
call test_c_taylor (taylor1, taylor2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (taylor2 == taylor_out) then
  print *, 'C_side_convert: taylor C to F: OK'
else
  print *, 'C_SIDE_CONVERT: taylor C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! sr_table_wake Check

sr_table_wake1 = sr_table_wake_in

print *
call test_c_sr_table_wake (sr_table_wake1, sr_table_wake2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (sr_table_wake2 == sr_table_wake_out) then
  print *, 'C_side_convert: sr_table_wake C to F: OK'
else
  print *, 'C_SIDE_CONVERT: sr_table_wake C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! sr_mode_wake Check

sr_mode_wake1 = sr_mode_wake_in

print *
call test_c_sr_mode_wake (sr_mode_wake1, sr_mode_wake2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (sr_mode_wake2 == sr_mode_wake_out) then
  print *, 'C_side_convert: sr_mode_wake C to F: OK'
else
  print *, 'C_SIDE_CONVERT: sr_mode_wake C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! lr_wake Check

lr_wake1 = lr_wake_in

print *
call test_c_lr_wake (lr_wake1, lr_wake2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (lr_wake2 == lr_wake_out) then
  print *, 'C_side_convert: lr_wake C to F: OK'
else
  print *, 'C_SIDE_CONVERT: lr_wake C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! wake Check

wake1 = wake_in

print *
call test_c_wake (wake1, wake2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (wake2 == wake_out) then
  print *, 'C_side_convert: wake C to F: OK'
else
  print *, 'C_SIDE_CONVERT: wake C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! control Check

control1 = control_in

print *
call test_c_control (control1, control2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (control2 == control_out) then
  print *, 'C_side_convert: control C to F: OK'
else
  print *, 'C_SIDE_CONVERT: control C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! param Check

param1 = param_in

print *
call test_c_param (param1, param2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (param2 == param_out) then
  print *, 'C_side_convert: param C to F: OK'
else
  print *, 'C_SIDE_CONVERT: param C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! amode Check

amode1 = amode_in

print *
call test_c_amode (amode1, amode2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (amode2 == amode_out) then
  print *, 'C_side_convert: amode C to F: OK'
else
  print *, 'C_SIDE_CONVERT: amode C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! linac_mode Check

linac_mode1 = linac_mode_in

print *
call test_c_linac_mode (linac_mode1, linac_mode2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (linac_mode2 == linac_mode_out) then
  print *, 'C_side_convert: linac_mode C to F: OK'
else
  print *, 'C_SIDE_CONVERT: linac_mode C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! modes Check

modes1 = modes_in

print *
call test_c_modes (modes1, modes2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (modes2 == modes_out) then
  print *, 'C_side_convert: modes C to F: OK'
else
  print *, 'C_SIDE_CONVERT: modes C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! bmad_com Check

bmad_com = bmad_com_in

print *
call test_c_bmad_com (c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (bmad_com == bmad_com_out) then
  print *, 'C_side_convert: bmad_com C to F: OK'
else
  print *, 'C_SIDE_CONVERT: bmad_com C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! em_field Check

em_field1 = em_field_in

print *
call test_c_em_field (em_field1, em_field2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (em_field2 == em_field_out) then
  print *, 'C_side_convert: em_field C to F: OK'
else
  print *, 'C_SIDE_CONVERT: em_field C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! ele Check

ele1 = ele_in
ele1%ix_ele = ele_in%ix_ele
ele1%gen_field => ele_in%gen_field

print *
call test_c_ele (ele1, ele2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (ele2 == ele_out) then
  print *, 'C_side_convert: ele C to F: OK'
else
  print *, 'C_SIDE_CONVERT: ele C to F: FAILED!!'
  call print_eq_ele (ele2, ele_out)
  all_ok = .false.
endif

!------------------------------------------------
! mode_info Check

mode_info1 = mode_info_in

print *
call test_c_mode_info (mode_info1, mode_info2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (mode_info2 == mode_info_out) then
  print *, 'C_side_convert: mode_info C to F: OK'
else
  print *, 'C_SIDE_CONVERT: mode_info C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! lat Check

lat1 = lat_in

print *
call test_C_lat (lat1, lat2, c_ok)
all_ok = all_ok .and. f_logic(c_ok)

if (lat2 == lat_out) then
  print *, 'C_side_convert: lat C to F: OK'
else
  print *, 'C_SIDE_CONVERT: lat C to F: FAILED!!'
  all_ok = .false.
endif

!------------------------------------------------
! All Check

print *
if (all_ok) then
  print *, 'Convert check: All OK'
else
  print *, 'CONVERT CHECK: ALL NOT OK!!'
endif

end program

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_coord (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (coord_struct) f1, f2

!

call coord_to_f (c1, f1)

if (f1 == coord_in) then
  print *, 'F_side_convert: Coord C to F: OK'
else
  print *, 'F_SIDE_CONVERT: COORD C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = coord_out
call coord_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_twiss (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (twiss_struct) f1, f2

!

call twiss_to_f (c1, f1)

if (f1 == twiss_in) then
  print *, 'F_side_convert: Twiss C to F: OK'
else
  print *, 'F_SIDE_CONVERT: TWISS C TO F: FAILED!!'
  all_ok = .false.
endif


f2 = twiss_out
call twiss_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_xy_disp (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (xy_disp_struct) f1, f2

!

call xy_disp_to_f (c1, f1)

if (f1 == xy_disp_in) then
  print *, 'F_side_convert: Xy_disp C to F: OK'
else
  print *, 'F_SIDE_CONVERT: XY_DISP C TO F: FAILED!!'
  all_ok = .false.
endif


f2 = xy_disp_out
call xy_disp_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_floor_position (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (floor_position_struct) f1, f2

!

call floor_position_to_f (c1, f1)

if (f1 == floor_position_in) then
  print *, 'F_side_convert: floor_position C to F: OK'
else
  print *, 'F_SIDE_CONVERT: floor_position C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = floor_position_out
call floor_position_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_wig_term (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (wig_term_struct) f1, f2

!

call wig_term_to_f (c1, f1)

if (f1 == wig_term_in) then
  print *, 'F_side_convert: wig_term C to F: OK'
else
  print *, 'F_SIDE_CONVERT: wig_term C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = wig_term_out
call wig_term_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_taylor_term (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (taylor_term_struct) f1, f2

!

call taylor_term_to_f (c1, f1)

if (f1 == taylor_term_in) then
  print *, 'F_side_convert: taylor_term C to F: OK'
else
  print *, 'F_SIDE_CONVERT: taylor_term C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = taylor_term_out
call taylor_term_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_taylor (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (taylor_struct) f1, f2

!

call taylor_to_f (c1, f1)

if (f1 == taylor_in) then
  print *, 'F_side_convert: taylor C to F: OK'
else
  print *, 'F_SIDE_CONVERT: taylor C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = taylor_out
call taylor_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_sr_table_wake (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (sr_table_wake_struct) f1, f2

!

call sr_table_wake_to_f (c1, f1)

if (f1 == sr_table_wake_in) then
  print *, 'F_side_convert: sr_table_wake C to F: OK'
else
  print *, 'F_SIDE_CONVERT: sr_table_wake C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = sr_table_wake_out
call sr_table_wake_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_sr_mode_wake (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (sr_mode_wake_struct) f1, f2

!

call sr_mode_wake_to_f (c1, f1)

if (f1 == sr_mode_wake_in) then
  print *, 'F_side_convert: sr_mode_wake C to F: OK'
else
  print *, 'F_SIDE_CONVERT: sr_mode_wake C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = sr_mode_wake_out
call sr_mode_wake_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_lr_wake (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (lr_wake_struct) f1, f2

!

call lr_wake_to_f (c1, f1)

if (f1 == lr_wake_in) then
  print *, 'F_side_convert: lr_wake C to F: OK'
else
  print *, 'F_SIDE_CONVERT: lr_wake C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = lr_wake_out
call lr_wake_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_wake (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (wake_struct) f1, f2

!

call wake_to_f (c1, f1)

if (f1 == wake_in) then
  print *, 'F_side_convert: wake C to F: OK'
else
  print *, 'F_SIDE_CONVERT: wake C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = wake_out
call wake_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_control (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (control_struct) f1, f2

!

call control_to_f (c1, f1)

if (f1 == control_in) then
  print *, 'F_side_convert: control C to F: OK'
else
  print *, 'F_SIDE_CONVERT: control C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = control_out
call control_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_param (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (lat_param_struct) f1, f2

!

call param_to_f (c1, f1)

if (f1 == param_in) then
  print *, 'F_side_convert: param C to F: OK'
else
  print *, 'F_SIDE_CONVERT: param C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = param_out
call param_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_amode (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (anormal_mode_struct) f1, f2

!

call amode_to_f (c1, f1)

if (f1 == amode_in) then
  print *, 'F_side_convert: amode C to F: OK'
else
  print *, 'F_SIDE_CONVERT: amode C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = amode_out
call amode_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_linac_mode (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (linac_normal_mode_struct) f1, f2

!

call linac_mode_to_f (c1, f1)

if (f1 == linac_mode_in) then
  print *, 'F_side_convert: linac_mode C to F: OK'
else
  print *, 'F_SIDE_CONVERT: linac_mode C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = linac_mode_out
call linac_mode_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_modes (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (normal_modes_struct) f1, f2

!

call modes_to_f (c1, f1)

if (f1 == modes_in) then
  print *, 'F_side_convert: modes C to F: OK'
else
  print *, 'F_SIDE_CONVERT: modes C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = modes_out
call modes_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_bmad_com (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2

!

call bmad_com_to_f (c1)

if (bmad_com == bmad_com_in) then
  print *, 'F_side_convert: bmad_com C to F: OK'
else
  print *, 'F_SIDE_CONVERT: bmad_com C TO F: FAILED!!'
  all_ok = .false.
endif

bmad_com = bmad_com_out
call bmad_com_to_c (c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_em_field (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (em_field_struct) f1, f2

!

call em_field_to_f (c1, f1)

if (f1 == em_field_in) then
  print *, 'F_side_convert: em_field C to F: OK'
else
  print *, 'F_SIDE_CONVERT: em_field C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = em_field_out
call em_field_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_ele (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (ele_struct) f1, f2

!

call ele_to_f (c1, f1)

if (f1 == ele_in) then
  print *, 'F_side_convert: ele C to F: OK'
else
  print *, 'F_SIDE_CONVERT: ele C TO F: FAILED!!'
  call print_eq_ele (f1, ele_in)
  all_ok = .false.
endif

f2 = ele_out
f2%ix_ele = ele_out%ix_ele
f2%gen_field => ele_out%gen_field
call ele_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_mode_info (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (mode_info_struct) f1, f2

!

call mode_info_to_f (c1, f1)

if (f1 == mode_info_in) then
  print *, 'F_side_convert: mode_info C to F: OK'
else
  print *, 'F_SIDE_CONVERT: mode_info C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = mode_info_out
call mode_info_to_c (f2, c2)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine test_f_lat (c1, c2)

use fortran_and_cpp
use equality_mod
use test_mod

implicit none

type (c_dummy_struct) c1, c2
type (lat_struct) f1, f2

!

call lat_to_f (c1, f1)

if (f1 == lat_in) then
  print *, 'F_side_convert: lat C to F: OK'
else
  print *, 'F_SIDE_CONVERT: lat C TO F: FAILED!!'
  all_ok = .false.
endif

f2 = lat_out
call lat_to_c (f2, c2)

end subroutine

