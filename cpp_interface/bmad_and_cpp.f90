subroutine init_coord_struct (f_coord)

use bmad_struct

implicit none

type (coord_struct), pointer :: f_coord
allocate (f_coord)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine coord_to_c (f_coord, c_coord)
!
! Subroutine to convert a Bmad coord_struct to a C++ C_coord.
!
! Input:
!   f_coord -- Coord_struct: Input Bmad coord_struct.
!
! Output:
!   c_coord -- c_dummy_struct: Output C_coord.
!-

subroutine coord_to_c (f_coord, c_coord)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (coord_struct) f_coord
type (c_dummy_struct) c_coord

call coord_to_c2 (c_coord, f_coord%vec, f_coord%s, f_coord%t, &
          real(f_coord%spin(1)), aimag(f_coord%spin(1)), real(f_coord%spin(2)), aimag(f_coord%spin(2)), &
          f_coord%e_field_x, f_coord%e_field_y, f_coord%phase_x, f_coord%phase_y)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine coord_to_f2 (f_coord, vec, s, t, sp1_re, sp1_im, sp2_re, sp2_im, e_field_x, e_field_y, phase_x, phase_y)
!
! Subroutine used by coord_to_f to convert a C++ C_coord into
! a Bmad coord_struct. This routine is not for general use.
!-

subroutine coord_to_f2 (f_coord, vec, s, t, sp1_re, sp1_im, sp2_re, sp2_im, e_field_x, e_field_y, phase_x, phase_y)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (coord_struct) f_coord
real(rp) vec(6), e_field_x, e_field_y, phase_x, phase_y, s, t, sp1_re, sp1_im, sp2_re, sp2_im

f_coord = coord_struct(vec, s, t, [cmplx(sp1_re, sp1_im), cmplx(sp2_re, sp2_im)], &
                                    e_field_x, e_field_y, phase_x, phase_y)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine twiss_to_c (f_twiss, c_twiss)
!
! Subroutine to convert a Bmad twiss_struct to a C++ C_twiss.
!
! Input:
!   f_twiss -- Twiss_struct: Input Bmad twiss_struct.
!
! Output:
!   c_twiss -- c_dummy_struct: Output C_twiss.
!-

subroutine twiss_to_c (f_twiss, c_twiss)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (twiss_struct), target :: f_twiss
type (twiss_struct), pointer :: f
type (c_dummy_struct) c_twiss

f => f_twiss
call twiss_to_c2 (c_twiss, f%beta, f%alpha, f%gamma, &
                    f%phi, f%eta, f%etap, f%sigma, f%sigma_p, f%emit, f%norm_emit)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine twiss_to_f2 (f_twiss, beta, alpha, gamma, phi, eta, etap, 
!                                           sigma, sigma_p, emit, norm_emit)
!
! Subroutine used by twiss_to_f to convert a C++ C_twiss into
! a Bmad twiss_struct. This routine is not for general use.
!-

subroutine twiss_to_f2 (f_twiss, beta, alpha, gamma, phi, eta, etap, &
                                             sigma, sigma_p, emit, norm_emit)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (twiss_struct) f_twiss
real(rp) beta, alpha, gamma, phi, eta, etap, sigma, emit, sigma_p, norm_emit

f_twiss = twiss_struct(beta, alpha, gamma, phi, eta, etap, &
                                             sigma, sigma_p, emit, norm_emit)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine xy_disp_to_c (f_xy_disp, c_xy_disp)
!
! Subroutine to convert a Bmad xy_disp_struct to a C++ C_xy_disp.
!
! Input:
!   f_xy_disp -- Xy_disp_struct: Input Bmad xy_disp_struct.
!
! Output:
!   c_xy_disp -- c_dummy_struct: Output C_xy_disp.
!-

subroutine xy_disp_to_c (f_xy_disp, c_xy_disp)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (xy_disp_struct), target :: f_xy_disp
type (xy_disp_struct), pointer :: f
type (c_dummy_struct) c_xy_disp

f => f_xy_disp
call xy_disp_to_c2 (c_xy_disp, f%eta, f%etap)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine xy_disp_to_f2 (f_xy_disp, eta, etap)
!
! Subroutine used by xy_disp_to_f to convert a C++ C_xy_disp into
! a Bmad xy_disp_struct. This routine is not for general use.
!-

subroutine xy_disp_to_f2 (f_xy_disp, eta, etap)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (xy_disp_struct) f_xy_disp
real(rp) eta, etap

f_xy_disp = xy_disp_struct(eta, etap)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine floor_position_to_c (f_floor_position, c_floor_position)
!
! Subroutine to convert a Bmad floor_position_struct to a C++ C_floor_position.
!
! Input:
!   f_floor_position -- Floor_position_struct: Input. 
!
! Output:
!   c_floor_position -- c_dummy_struct: Output C_floor_position.
!-

subroutine floor_position_to_c (f_floor_position, c_floor_position)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (floor_position_struct), target :: f_floor_position
type (floor_position_struct), pointer :: f
type (c_dummy_struct) c_floor_position

f => f_floor_position
call floor_position_to_c2 (c_floor_position, f%x, f%y, f%z, &
                                                       f%theta, f%phi, f%psi)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine floor_position_to_f2 (f_floor_position, x, y, z, theta, phi, psi)
!
! Subroutine used by floor_position_to_f to convert a C++ C_floor_position into
! a Bmad floor_position_struct. This routine is not for general use.
!-

subroutine floor_position_to_f2 (f_floor_position, x, y, z, theta, phi, psi)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (floor_position_struct) f_floor_position
real(rp) x, y, z, theta, phi, psi

f_floor_position = floor_position_struct(x, y, z, theta, phi, psi)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine wig_term_to_c (f_wig_term, c_wig_term)
!
! Subroutine to convert a Bmad wig_term_struct to a C++ C_wig_term.
!
! Input:
!   f_wig_term -- Wig_term_struct: Input Bmad wig_term_struct.
!
! Output:
!   c_wig_term -- c_dummy_struct: Output C_wig_term.
!-

subroutine wig_term_to_c (f_wig_term, c_wig_term)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (wig_term_struct), target :: f_wig_term
type (wig_term_struct), pointer :: f
type (c_dummy_struct) c_wig_term

f => f_wig_term
call wig_term_to_c2 (c_wig_term, f%coef, f%kx, f%ky, f%kz, f%phi_z, f%type)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine wig_term_to_f2 (f_wig_term, coef, kx, ky, kz, phi_z, tp)
!
! Subroutine used by wig_term_to_f to convert a C++ C_wig_term into
! a Bmad wig_term_struct. This routine is not for general use.
!-

subroutine wig_term_to_f2 (f_wig_term, coef, kx, ky, kz, phi_z, tp)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (wig_term_struct) f_wig_term
real(rp) coef, kx, ky, kz, phi_z
integer tp

f_wig_term = wig_term_struct(coef, kx, ky, kz, phi_z, tp)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine taylor_term_to_c (f_taylor_term, c_taylor_term)
!
! Subroutine to convert a Bmad taylor_term_struct to a C++ C_taylor_term.
!
! Input:
!   f_taylor_term -- Taylor_term_struct: Input Bmad taylor_term_struct.
!
! Output:
!   c_taylor_term -- c_dummy_struct: Output C_taylor_term.
!-

subroutine taylor_term_to_c (f_taylor_term, c_taylor_term)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (taylor_term_struct), target :: f_taylor_term
type (taylor_term_struct), pointer :: f
type (c_dummy_struct) c_taylor_term

f => f_taylor_term
call taylor_term_to_c2 (c_taylor_term, f%coef, f%expn)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine taylor_term_to_f2 (f_taylor_term, coef, exp)
!
! Subroutine used by taylor_term_to_f to convert a C++ C_taylor_term into
! a Bmad taylor_term_struct. This routine is not for general use.
!-

subroutine taylor_term_to_f2 (f_taylor_term, coef, exp)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (taylor_term_struct) f_taylor_term
real(rp) coef
integer exp(6)

f_taylor_term = taylor_term_struct(coef, exp)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine taylor_to_c (f_taylor, c_taylor)
!
! Subroutine to convert a Bmad taylor_struct to a C++ C_taylor.
!
! Input:
!   f_taylor -- Taylor_struct: Input Bmad taylor_struct.
!
! Output:
!   c_taylor -- c_dummy_struct: Output C_taylor.
!-

subroutine taylor_to_c (f_taylor, c_taylor)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (taylor_struct), target :: f_taylor
type (taylor_struct), pointer :: f
type (c_dummy_struct) c_taylor
integer i, n_term


f => f_taylor

n_term = 0
if (associated(f%term)) n_term = size(f%term)

call taylor_to_c2 (c_taylor, n_term, f%ref)

do i = 1, n_term
  call taylor_term_in_taylor_to_c2 (c_taylor, i, f%term(i)%coef, f%term(i)%expn)
end do

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine taylor_to_f2 (f_taylor, n_term, ref)
!
! Subroutine used by taylor_to_f to convert a C++ C_taylor into
! a Bmad taylor_struct. This routine is not for general use.
!-

subroutine taylor_to_f2 (f_taylor, n_term, ref)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (taylor_struct) f_taylor
real(rp) ref
integer n_term

call init_taylor_series (f_taylor, n_term)
f_taylor%ref = ref

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine taylor_term_in_taylor_to_f2 (f_taylor, it, coef, exp)
!
! Subroutine used by taylor_to_f to convert a C++ C_taylor into
! a Bmad taylor_struct. This routine is not for general use.
!-

subroutine taylor_term_in_taylor_to_f2 (f_taylor, it, coef, exp)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (taylor_struct) f_taylor
real(rp) coef
integer it, exp(6)

f_taylor%term(it) = taylor_term_struct(coef, exp)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine rf_wake_sr_table_to_c (f_rf_wake_sr_table, c_rf_wake_sr_table)
!
! Subroutine to convert a Bmad rf_wake_sr_table_struct to a C++ C_rf_wake_sr_table.
!
! Input:
!   f_rf_wake_sr_table -- rf_wake_sr_table_struct: Input Bmad rf_wake_sr_table_struct.
!
! Output:
!   c_rf_wake_sr_table -- c_dummy_struct: Output C_rf_wake_sr_table.
!-

subroutine rf_wake_sr_table_to_c (f_rf_wake_sr_table, c_rf_wake_sr_table)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (rf_wake_sr_table_struct), target :: f_rf_wake_sr_table
type (rf_wake_sr_table_struct), pointer :: f
type (c_dummy_struct) c_rf_wake_sr_table

f => f_rf_wake_sr_table
call rf_wake_sr_table_to_c2 (c_rf_wake_sr_table, f%z, f%long, f%trans)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine rf_wake_sr_table_to_f2 (f_rf_wake_sr_table, z, long, trans)
!
! Subroutine used by rf_wake_sr_table_to_f to convert a C++ C_rf_wake_sr_table into
! a Bmad rf_wake_sr_table_struct. This routine is not for general use.
!-

subroutine rf_wake_sr_table_to_f2 (f_rf_wake_sr_table, z, long, trans)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (rf_wake_sr_table_struct) f_rf_wake_sr_table
real(rp) z, long, trans

f_rf_wake_sr_table = rf_wake_sr_table_struct(z, long, trans)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine rf_wake_sr_mode_to_c (f_rf_wake_sr_mode, c_rf_wake_sr_mode)
!
! Subroutine to convert a Bmad rf_wake_sr_mode_struct to a C++ C_rf_wake_sr_mode.
!
! Input:
!   f_rf_wake_sr_mode -- rf_wake_sr_mode_struct: Input Bmad rf_wake_sr_mode_struct.
!
! Output:
!   c_rf_wake_sr_mode -- c_dummy_struct: Output C_rf_wake_sr_mode.
!-

subroutine rf_wake_sr_mode_to_c (f_rf_wake_sr_mode, c_rf_wake_sr_mode)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (rf_wake_sr_mode_struct), target :: f_rf_wake_sr_mode
type (rf_wake_sr_mode_struct), pointer :: f
type (c_dummy_struct) c_rf_wake_sr_mode

f => f_rf_wake_sr_mode
call rf_wake_sr_mode_to_c2 (c_rf_wake_sr_mode, f%amp, f%damp, f%k, f%phi, &
                          f%b_sin, f%b_cos, f%a_sin, f%a_cos)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine rf_wake_sr_mode_to_f2 (f_rf_wake_sr_mode, amp, damp, freq, phi, 
!                                     b_sin, b_cos, a_sin, a_cos)
!
! Subroutine used by rf_wake_sr_mode_to_f to convert a C++ C_rf_wake_sr_mode into
! a Bmad rf_wake_sr_mode_struct. This routine is not for general use.
!-

subroutine rf_wake_sr_mode_to_f2 (f_rf_wake_sr_mode, amp, damp, freq, phi, &
                                        b_sin, b_cos, a_sin, a_cos)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (rf_wake_sr_mode_struct) f_rf_wake_sr_mode
real(rp) amp, damp, freq, phi, b_sin, b_cos, a_sin, a_cos

f_rf_wake_sr_mode = rf_wake_sr_mode_struct(amp, damp, freq, phi, &
                                      b_sin, b_cos, a_sin, a_cos)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine rf_wake_lr_to_c (f_rf_wake_lr, c_rf_wake_lr)
!
! Subroutine to convert a Bmad rf_wake_lr_struct to a C++ C_rf_wake_lr.
!
! Input:
!   f_rf_wake_lr -- Rf_wake_lr_struct: Input Bmad rf_wake_lr_struct.
!
! Output:
!   c_rf_wake_lr -- c_dummy_struct: Output C_rf_wake_lr.
!-

subroutine rf_wake_lr_to_c (f_rf_wake_lr, c_rf_wake_lr)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (rf_wake_lr_struct), target :: f_rf_wake_lr
type (rf_wake_lr_struct), pointer :: f
type (c_dummy_struct) c_rf_wake_lr

f => f_rf_wake_lr
call rf_wake_lr_to_c2 (c_rf_wake_lr, f%freq, f%freq_in, f%R_over_Q, f%q, f%angle, &
         f%b_sin, f%b_cos, f%a_sin, f%a_cos, f%t_ref, f%m, f%polarized)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine rf_wake_lr_to_f2 (f_rf_wake_lr, freq, freq_in, r_over_q, q, angle, &
!                                   n_sin, n_cos, s_cos, s_sin, t_ref, m, polarized)
!
! Subroutine used by rf_wake_lr_to_f to convert a C++ C_rf_wake_lr into
! a Bmad rf_wake_lr_struct. This routine is not for general use.
!-

subroutine rf_wake_lr_to_f2 (f_rf_wake_lr, freq, freq_in, r_over_q, q, angle, &
                                   n_sin, n_cos, s_cos, s_sin, t_ref, m, polarized)


use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (rf_wake_lr_struct) f_rf_wake_lr
real(rp) freq, freq_in, r_over_q, q, n_sin, n_cos, s_cos, s_sin, s_ref
real(rp) angle, t_ref
integer m
logical polarized

f_rf_wake_lr = rf_wake_lr_struct(freq, freq_in, r_over_q, q, angle, &
                                   n_sin, n_cos, s_cos, s_sin, t_ref, m, polarized)

end subroutine


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine rf_wake_to_c (f_rf_wake, c_rf_wake)
!
! Subroutine to convert a Bmad rf_wake_struct to a C++ C_rf_wake.
!
! Input:
!   f_rf_wake -- Rf_wake_struct: Input Bmad rf_wake_struct.
!
! Output:
!   c_rf_wake -- c_dummy_struct: Output C_rf_wake.
!-

subroutine rf_wake_to_c (f_rf_wake, c_rf_wake)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (rf_wake_struct), target :: f_rf_wake
type (rf_wake_struct), pointer :: f
type (rf_wake_sr_mode_struct), pointer :: sr_mode
type (c_dummy_struct) c_rf_wake
integer i, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr

!

f => f_rf_wake
n_sr_table       = size(f%sr_table)
n_sr_mode_long  = size(f%sr_mode_long)
n_sr_mode_trans = size(f%sr_mode_trans)
n_lr        = size(f%lr)

call rf_wake_to_c2 (c_rf_wake, c_str(f%sr_file), c_str(f%lr_file), f%z_sr_mode_max, &
                                             n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr)

do i = 0, n_sr_table-1
  call rf_wake_sr_table_in_rf_wake_to_c2 (c_rf_wake, i, f%sr_table(i)%z, f%sr_table(i)%long, f%sr_table(i)%trans)
enddo

do i = 1, n_sr_mode_long
  sr_mode => f%sr_mode_long(i)
  call rf_wake_sr_mode_long_in_rf_wake_to_c2 (c_rf_wake, i, sr_mode%amp, sr_mode%damp, sr_mode%k, &
                  sr_mode%phi, sr_mode%b_sin, sr_mode%b_cos, sr_mode%a_sin, sr_mode%a_cos)
enddo

do i = 1, n_sr_mode_trans
  sr_mode => f%sr_mode_trans(i)
  call rf_wake_sr_mode_trans_in_rf_wake_to_c2 (c_rf_wake, i, sr_mode%amp, sr_mode%damp, sr_mode%k, &
                  sr_mode%phi, sr_mode%b_sin, sr_mode%b_cos, sr_mode%a_sin, sr_mode%a_cos)
enddo

do i = 1, n_lr
  call rf_wake_lr_in_rf_wake_to_c2 (c_rf_wake, i, f%lr(i)%freq, f%lr(i)%freq_in, &
         f%lr(i)%r_over_q, f%lr(i)%Q, f%lr(i)%angle, f%lr(i)%b_sin, &
         f%lr(i)%b_cos, f%lr(i)%a_sin, f%lr(i)%a_cos, f%lr(i)%t_ref, f%lr(i)%m, &
         f%lr(i)%polarized)
enddo 

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine rf_wake_to_f2 (f_rf_wake, sr_file, n_srf, lr_file, n_lrf, z_sr_mode_max,
!                                            n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr)
!
! Subroutine used by rf_wake_to_f to convert a C++ C_rf_wake into
! a Bmad rf_wake_struct. This routine is not for general use.
!-

subroutine rf_wake_to_f2 (f_rf_wake, sr_file, n_srf, lr_file, n_lrf, z_sr_mode_max, &
                                         n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (rf_wake_struct) f_rf_wake
integer n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr, n_srf, n_lrf
real(rp) z_sr_mode_max
character(n_srf) :: sr_file
character(n_lrf) :: lr_file

!

f_rf_wake%sr_file = sr_file
f_rf_wake%lr_file = lr_file
f_rf_wake%z_sr_mode_max = z_sr_mode_max
allocate (f_rf_wake%sr_table(0:n_sr_table-1))
allocate (f_rf_wake%sr_mode_long(n_sr_mode_long))
allocate (f_rf_wake%sr_mode_trans(n_sr_mode_trans))
allocate (f_rf_wake%lr(n_lr))

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine rf_wake_sr_table_in_rf_wake_to_f2 (f_rf_wake, it, z, long, trans)
!
! Subroutine used by rf_wake_to_f to convert a C++ C_rf_wake into
! a Bmad rf_wake_struct. This routine is not for general use.
!-

subroutine rf_wake_sr_table_in_rf_wake_to_f2 (f_rf_wake, it, z, long, trans)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (rf_wake_struct) f_rf_wake
real(rp) z, long, trans
integer it

f_rf_wake%sr_table(it) = rf_wake_sr_table_struct(z, long, trans)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine rf_wake_sr_mode_long_in_rf_wake_to_f2 (f_rf_wake, it, amp, damp, freq, phi, &
!                                        b_sin, b_cos, a_sin, a_cos)
!
! Subroutine used by rf_wake_to_f to convert a C++ C_rf_wake into
! a Bmad rf_wake_struct. This routine is not for general use.
!-

subroutine rf_wake_sr_mode_long_in_rf_wake_to_f2 (f_rf_wake, it, amp, damp, freq, phi, &
                                        b_sin, b_cos, a_sin, a_cos)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (rf_wake_struct) f_rf_wake
real(rp) amp, damp, freq, phi, b_sin, b_cos, a_sin, a_cos
integer it

f_rf_wake%sr_mode_long(it) = rf_wake_sr_mode_struct(amp, damp, freq, phi, &
                                        b_sin, b_cos, a_sin, a_cos)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine rf_wake_sr_mode_trans_in_rf_wake_to_f2 (f_rf_wake, it, amp, damp, freq, phi, &
!                                        b_sin, b_cos, a_sin, a_cos)
!
! Subroutine used by rf_wake_to_f to convert a C++ C_rf_wake into
! a Bmad rf_wake_struct. This routine is not for general use.
!-

subroutine rf_wake_sr_mode_trans_in_rf_wake_to_f2 (f_rf_wake, it, amp, damp, freq, phi, &
                                        b_sin, b_cos, a_sin, a_cos)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (rf_wake_struct) f_rf_wake
real(rp) amp, damp, freq, phi, b_sin, b_cos, a_sin, a_cos
integer it

f_rf_wake%sr_mode_trans(it) = rf_wake_sr_mode_struct(amp, damp, freq, phi, &
                                        b_sin, b_cos, a_sin, a_cos)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine rf_wake_lr_in_rf_wake_to_f2 (f_rf_wake, it, freq, freq_in, r_over_q, q, angle, &
!                                     n_sin, n_cos, s_cos, s_sin, t_ref, m, polarized)
!
! Subroutine used by rf_wake_to_f to convert a C++ C_rf_wake into
! a Bmad rf_wake_struct. This routine is not for general use.
!-

subroutine rf_wake_lr_in_rf_wake_to_f2 (f_rf_wake, it, freq, freq_in, r_over_q, q, angle, &
                                    n_sin, n_cos, s_cos, s_sin, t_ref, m, polarized)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (rf_wake_struct) f_rf_wake
real(rp) freq, freq_in, r_over_q, q, n_sin, n_cos, s_cos, s_sin, angle, t_ref
integer it, m
logical polarized

f_rf_wake%lr(it) = rf_wake_lr_struct(freq, freq_in, r_over_q, q, angle, &
                              n_sin, n_cos, s_cos, s_sin, t_ref, m, polarized)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine control_to_c (f_control, c_control)
!
! Subroutine to convert a Bmad control_struct to a C++ C_control.
!
! Input:
!   f_control -- Control_struct: Input Bmad control_struct.
!
! Output:
!   c_control -- c_dummy_struct: Output C_control.
!-

subroutine control_to_c (f_control, c_control)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (control_struct), target :: f_control
type (control_struct), pointer :: f
type (c_dummy_struct) c_control

f => f_control
call control_to_c2 (c_control, f%coef, f%ix_lord, f%ix_slave, f%ix_branch, f%ix_attrib)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine control_to_f2 (f_control, coef, ix_lord, ix_slave, ix_attrib)
!
! Subroutine used by control_to_f to convert a C++ C_control into
! a Bmad control_struct. This routine is not for general use.
!-

subroutine control_to_f2 (f_control, coef, ix_lord, ix_slave, ix_branch, ix_attrib)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (control_struct) f_control
real(rp) coef
integer ix_lord, ix_slave, ix_branch, ix_attrib

f_control = control_struct(coef, ix_lord, ix_slave, ix_branch, ix_attrib)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine lat_param_to_c (f_lat_param, c_lat_param)
!
! Subroutine to convert a Bmad lat_param_struct to a C++ C_lat_param.
!
! Input:
!   f_lat_param -- lat_param_struct: Input Bmad lat_param_struct.
!
! Output:
!   c_lat_param -- c_dummy_struct: Output C_lat_param.
!-

subroutine lat_param_to_c (f_lat_param, c_lat_param)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (lat_param_struct), target :: f_lat_param
type (lat_param_struct), pointer :: f
type (c_dummy_struct) c_lat_param

f => f_lat_param

call lat_param_to_c2 (c_lat_param, f%n_part, f%total_length, f%unstable_factor, &
      mat2arr(f%t1_with_RF), mat2arr(f%t1_no_RF), &
      f%particle, f%ix_lost, f%end_lost_at, f%plane_lost_at, f%lattice_type, &
      f%ixx, c_logic(f%stable), c_logic(f%aperture_limit_on), c_logic(f%lost))

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine lat_param_to_f2 (f_lat_param, n_part, total_length, &
!      growth_rate, m1, m2, particle, ix_lost, end_lost_at, plane_lost_at, &
!      lat_type, ixx, stable, ap_limit_on, lost)
!
! Subroutine used by lat_param_to_f to convert a C++ C_lat_param into
! a Bmad lat_param_struct. This routine is not for general use.
!-

subroutine lat_param_to_f2 (f_lat_param, n_part, total_length, &
      growth_rate, m1, m2, particle, ix_lost, end_lost_at, plane_lost_at, &
      lat_type, ixx, stable, ap_limit_on, lost) 

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (lat_param_struct) f_lat_param
real(rp) n_part, total_length, growth_rate
real(rp) m1(36), m2(36)
integer particle, ix_lost, end_lost_at, lat_type, ixx, stable, &
        ap_limit_on, lost, plane_lost_at

! Added status component

!f_lat_param = lat_param_struct(n_part, total_length, growth_rate, &
!      arr2mat(m1, 6, 6), arr2mat(m2, 6, 6), particle, ix_lost, end_lost_at, &
!      plane_lost_at, lat_type, ixx, f_logic(stable), f_logic(ap_limit_on), &
!      f_logic(lost))

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine anormal_mode_to_c (f_anormal_mode, c_anormal_mode)
!
! Subroutine to convert a Bmad anormal_mode_struct to a C++ C_anormal_mode.
!
! Input:
!   f_anormal_mode -- Anormal_mode_struct: Input Bmad anormal_mode_struct.
!
! Output:
!   c_anormal_mode -- c_dummy_struct: Output C_anormal_mode.
!-

subroutine anormal_mode_to_c (f_anormal_mode, c_anormal_mode)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (anormal_mode_struct), target :: f_anormal_mode
type (anormal_mode_struct), pointer :: f
type (c_dummy_struct) c_anormal_mode

f => f_anormal_mode
call anormal_mode_to_c2 (c_anormal_mode, f%emittance, f%synch_int(4), f%synch_int(5), f%synch_int(6), &
                            f%j_damp, f%alpha_damp, f%chrom, f%tune)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine anormal_mode_to_f2 (f_anormal_mode, emit, s_int4, s_int5, s_int6, j_damp, a_damp, chrom, tune)
!
! Subroutine used by anormal_mode_to_f to convert a C++ C_anormal_mode into
! a Bmad anormal_mode_struct. This routine is not for general use.
!-

subroutine anormal_mode_to_f2 (f_anormal_mode, emit, s_int4, s_int5, s_int6, j_damp, a_damp, chrom, tune)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (anormal_mode_struct) f_anormal_mode
real(rp) emit, s_int4, s_int5, s_int6, j_damp, a_damp, chrom, tune

f_anormal_mode = anormal_mode_struct(emit, [s_int4, s_int5, s_int6], j_damp, a_damp, chrom, tune)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine linac_normal_mode_to_c (f_linac_normal_mode, c_linac_normal_mode)
!
! Subroutine to convert a Bmad linac_normal_mode_struct to a C++ C_linac_normal_mode.
!
! Input:
!   f_linac_normal_mode -- Linac_normal_mode_struct: Input Bmad linac_normal_mode_struct.
!
! Output:
!   c_linac_normal_mode -- c_dummy_struct: Output C_linac_normal_mode.
!-

subroutine linac_normal_mode_to_c (f_linac_normal_mode, c_linac_normal_mode)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (linac_normal_mode_struct), target :: f_linac_normal_mode
type (linac_normal_mode_struct), pointer :: f
type (c_dummy_struct) c_linac_normal_mode

f => f_linac_normal_mode
call linac_normal_mode_to_c2 (c_linac_normal_mode, f%i2_E4, f%i3_E7, f%i5a_E6, f%i5b_E6, &
                                f%sig_E1, f%a_emittance_end, f%b_emittance_end)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine linac_normal_mode_to_f2 (f_linac_normal_mode, i2, i3, i5a, i5b, sig_e, ea, eb)
!
! Subroutine used by linac_normal_mode_to_f to convert a C++ C_linac_normal_mode into
! a Bmad linac_normal_mode_struct. This routine is not for general use.
!-

subroutine linac_normal_mode_to_f2 (f_linac_normal_mode, i2, i3, i5a, i5b, sig_e, ea, eb)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (linac_normal_mode_struct) f_linac_normal_mode
real(rp) i2, i3, i5a, i5b, sig_e, ea, eb

f_linac_normal_mode = linac_normal_mode_struct(i2, i3, i5a, i5b, sig_e, ea, eb)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine normal_modes_to_c (f_normal_modes, c_normal_modes)
!
! Subroutine to convert a Bmad normal_modes_struct to a C++ C_normal_modes.
!
! Input:
!   f_normal_modes -- normal_modes_struct: Input Bmad normal_modes_struct.
!
! Output:
!   c_normal_modes -- c_dummy_struct: Output C_normal_modes.
!-

subroutine normal_modes_to_c (f_normal_modes, c_normal_modes)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (normal_modes_struct), target :: f_normal_modes
type (normal_modes_struct), pointer :: f
type (c_dummy_struct) c_normal_modes

f => f_normal_modes
call normal_modes_to_c2 (c_normal_modes, f%synch_int, f%sige_e, f%sig_z, f%e_loss, f%rf_voltage, &
                           f%pz_aperture, f%a, f%b, f%z, f%lin)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine normal_modes_to_f2 (f_normal_modes, synch_int, sige, sig_z, e_loss, pz_aperture, a, b, z, lin)
!
! Subroutine used by normal_modes_to_f to convert a C++ C_normal_modes into
! a Bmad normal_modes_struct. This routine is not for general use.
!-

subroutine normal_modes_to_f2 (f_normal_modes, synch_int, sige, sig_z, e_loss, rf_volt, pz, a, b, z, lin)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (normal_modes_struct) f_normal_modes
type (c_dummy_struct) a, b, z, lin
real(rp) synch_int(0:3), sige, sig_z, e_loss, rf_volt, pz

!

call anormal_mode_to_f (a, f_normal_modes%a)
call anormal_mode_to_f (b, f_normal_modes%b)
call anormal_mode_to_f (z, f_normal_modes%z)
call linac_normal_mode_to_f (lin, f_normal_modes%lin)

f_normal_modes = normal_modes_struct(synch_int, sige, sig_z, e_loss, rf_volt, pz, &
                              f_normal_modes%a, f_normal_modes%b, f_normal_modes%z, f_normal_modes%lin)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine bmad_com_to_c (c_bmad_com)
!
! Subroutine to convert the Bmad bmad_common_struct common block to 
! a C++ C_bmad_com.
!
! Output:
!   c_bmad_com -- c_dummy_struct: Output C_bmad_com.
!-

subroutine bmad_com_to_c (c_bmad_com)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (bmad_common_struct) :: f
type (c_dummy_struct) c_bmad_com

f = bmad_com
call bmad_com_to_c2 (c_bmad_com, f%max_aperture_limit, f%d_orb, &
      f%default_ds_step, f%significant_length, &
      f%rel_tolerance, f%abs_tolerance, &
      f%rel_tol_adaptive_tracking, f%abs_tol_adaptive_tracking, &
      f%taylor_order, f%default_integ_order, &
      c_logic(f%canonical_coords), &
      c_logic(f%sr_wakes_on), c_logic(f%lr_wakes_on), &
      c_logic(f%mat6_track_symmetric), c_logic(f%auto_bookkeeper), & 
      c_logic(f%space_charge_on), c_logic(f%coherent_synch_rad_on), &
      c_logic(f%spin_tracking_on), &
      c_logic(f%radiation_damping_on), c_logic(f%radiation_fluctuations_on), &
      c_logic(f%conserve_taylor_maps))

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine bmad_com_to_f2 (max_ap, orb, ds_step, signif, rel, abs, rel_track, 
!              abs_track, taylor_ord, dflt_integ, cc, sr, lr, sym,
!              a_book, tsc_on, csr_on, st_on, rad_d, rad_f, ref_e, conserve_t)
!
! Subroutine used by bmad_com_to_f to transfer the data from a C++ 
! C_bmad_com variable into the Bmad bmad_com_stuct common block.
! This routine is not for general use.
!-

subroutine bmad_com_to_f2 (max_ap, orb, ds_step, signif, rel, abs, rel_track, &
        abs_track, taylor_ord, dflt_integ, cc, sr, lr, sym, &
        a_book, tsc_on, csr_on, st_on, rad_d, rad_f, conserve_t)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

real(rp) orb(6), max_ap, rel, abs, rel_track, abs_track, ds_step, signif
integer taylor_ord, dflt_integ, cc, sr, lr, sym
integer st_on, rad_d, rad_f, a_book, tsc_on, csr_on
integer conserve_t

bmad_com = bmad_common_struct(max_ap, orb, ds_step, signif, &
    rel, abs, rel_track, abs_track, taylor_ord, dflt_integ, &
    f_logic(cc), f_logic(sr), f_logic(lr), f_logic(sym), &
    f_logic(a_book), f_logic(tsc_on), f_logic(csr_on), f_logic(st_on), &
    f_logic(rad_d), f_logic(rad_f), f_logic(conserve_t))

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine em_field_to_c (f_em_field, c_em_field)
!
! Subroutine to convert a Bmad em_field_struct to a C++ C_em_field.
!
! Input:
!   f_em_field -- Em_field_struct: Input Bmad em_field_struct.
!
! Output:
!   c_em_field -- c_dummy_struct: Output C_em_field.
!-

subroutine em_field_to_c (f_em_field, c_em_field)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (em_field_struct), target :: f_em_field
type (em_field_struct), pointer :: f
type (c_dummy_struct) c_em_field

f => f_em_field

call em_field_to_c2 (c_em_field, f%E, f%b, mat2arr(f%dE), mat2arr(f%dB))

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine em_field_to_f2 (f_em_field, e, b, de, db)
!
! Subroutine used by em_field_to_f to convert a C++ C_em_field into
! a Bmad em_field_struct. This routine is not for general use.
!-

subroutine em_field_to_f2 (f_em_field, e, b, de, db)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (em_field_struct) f_em_field
real(rp) e(3), b(3), de(9), db(9)

f_em_field = em_field_struct(e, b, arr2mat(de, 3, 3), arr2mat(db, 3, 3))

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ele_to_c (f_ele, c_ele)
!
! Subroutine to convert a Bmad ele_struct to a C++ C_ele.
!
! Input:
!   f_ele -- Ele_struct: Input Bmad ele_struct.
!
! Output:
!   c_ele -- c_dummy_struct: Output C_ele.
!-

subroutine ele_to_c (f_ele, c_ele)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (ele_struct), target :: f_ele
type (ele_struct), pointer :: f
type (wig_term_struct), pointer :: w
type (c_dummy_struct) c_ele

real(rp) m6(36), value(0:n_attrib_maxx)
real(rp), allocatable :: r_arr(:)
integer i, nr1, nr2, nd, n_wig
character(200) descrip

!

f => f_ele

nr1 = 0; nr2 = 0
!if (associated(f%r)) then
!  nr1 = size(f%r, 1)
!  nr2 = size(f%r, 2)
!  allocate (r_arr(nr1*nr2))
!  r_arr = mat2arr (f%r)
!else
  allocate (r_arr(0))
!endif

descrip =  ' '
if (associated(f%descrip)) descrip = f%descrip

n_wig = 0
if (associated(f%wig_term)) n_wig = size(f%wig_term)

value(0) = 0
value(1:) = f%value

call ele_to_c2 (c_ele, c_str(f%name), c_str(f%type), c_str(f%alias), &
      c_str(f%component_name), c_str(descrip), f%a, f%b, f%z, f%x, f%y, &
      f%floor, f%map_ref_orb_in, f%map_ref_orb_out, f%gen_field, &
      f%taylor(1), f%taylor(2), f%taylor(3), f%taylor(4), f%taylor(5), f%taylor(6), &
      n_wig, value, f%gen0, f%vec0, mat2arr(f%mat6), mat2arr(f%c_mat), f%gamma_c, f%s, f%ref_time, &
      r_arr, nr1, nr2, f%a_pole, f%b_pole, r_size(f%a_pole), f%const, r_size(f%const), &
      f%key, f%sub_key, f%ix_ele, f%ix_branch, f%ix_value, &
      f%slave_status, f%n_slave, f%ix1_slave, f%ix2_slave, &
      f%lord_status, f%n_lord, f%ic1_lord, f%ic2_lord, &
      f%ix_pointer, f%ixx, &
      f%mat6_calc_method, f%tracking_method, f%field_calc, f%ref_orbit, &
      f%aperture_at, f%aperture_type, &
      f%symplectify, f%mode_flip, f%multipoles_on, f%scale_multipoles, f%map_with_offsets, &
      f%field_master, f%reversed, f%is_on, f%old_is_on, f%logic, f%bmad_logic, f%on_a_girder, &
      f%csr_calc_on, f%offset_moves_aperture)

if (associated(f%r)) deallocate(r_arr)

do i = 1, n_wig
  w => f_ele%wig_term(i)
  call wig_term_in_ele_to_c2 (c_ele, i, w%coef, w%kx, w%ky, w%kz, w%phi_z, w%type)
end do

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ele_to_f2 (f, nam, n_nam, typ, n_typ, ali, n_ali, component_nam, ...)
!
! Subroutine used by ele_to_f to convert a C++ C_ele into
! a Bmad ele_struct. This routine is not for general use.
!-

subroutine ele_to_f2 (f, nam, n_nam, typ, n_typ, ali, n_ali, component_nam, n_component_nam, des, n_des, &
    a, b, z, x, y, &
    floor, ref_orb_in, ref_orb_out, gen_f, &
    tlr1, tlr2, tlr3, tlr4, tlr5, tlr6, &
    n_wig, value, gen0, vec0, mat6, c_mat, &
    gamma_c, s, ref_t, r_arr, nr1, nr2, &
    a_pole, b_pole, n_ab, const, n_const, &
    key, sub_key, ix_ele, ix_branch, ix_value, &
    slave_status, n_slave, ix1_slave, ix2_slave, &
    lord_status, n_lord, ic1_lord, ic2_lord, &
    ix_point, ixx, mat6_meth, tracking_meth, field_calc, ref_orb, &
    aperture_at, aperture_type, symp, mode_flip, multi_on, scale_multi, &
    map_with_off, field_master, reversed, is_on, old_is_on, logic, bmad_logic, &
    girder, csr_calc, offset_moves_ap)   

use fortran_and_cpp
use multipole_mod
use bmad_struct
use bmad_interface

implicit none

type (ele_struct) f
type (c_dummy_struct) a, b, x, y, z, floor, ref_orb_in, ref_orb_out
type (c_dummy_struct) tlr1, tlr2, tlr3, tlr4, tlr5, tlr6
type (genfield), target :: gen_f

integer n_nam, nr1, nr2, n_ab, n_const, key, sub_key, lord_status, slave_status
integer ix2_slave, n_lord, ic1_lord, ic2_lord, ix_point, ixx, ix_ele, mat6_meth, tracking_meth, field_calc
integer ref_orb, aperture_at, symp, mode_flip, multi_on, map_with_off, field_master
integer reversed, is_on, old_is_on, logic, bmad_logic, girder, csr_calc, n_typ, n_ali, n_component_nam, n_des, ix_branch
integer n_wig, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr, aperture_type, offset_moves_ap
integer ix_value, n_slave, ix1_slave, scale_multi

real(rp) value(n_attrib_maxx), gen0(6), vec0(6), mat6(36), c_mat(4), gamma_c, s, ref_t
real(rp) a_pole(n_ab), b_pole(n_ab), r_arr(nr1*nr2), const(n_const)

character(n_nam)  nam
character(n_typ)  typ
character(n_ali)  ali
character(n_component_nam) component_nam
character(n_des)  des

!

call str_upcase(f%name, nam)
f%type  = typ
f%alias = ali
call str_upcase (f%component_name, component_nam)

if (n_des == 0) then
  if (associated (f%descrip)) deallocate (f%descrip)
else
  if (.not. associated(f%descrip)) allocate (f%descrip)
  f%descrip = des
endif

call twiss_to_f (a, f%a)
call twiss_to_f (b, f%b)
call twiss_to_f (z, f%z)
call xy_disp_to_f (x, f%x)
call xy_disp_to_f (y, f%y)
call floor_position_to_f (floor, f%floor)
call coord_to_f (ref_orb_in, f%map_ref_orb_in)
call coord_to_f (ref_orb_out, f%map_ref_orb_out)

f%gen_field => gen_f

call taylor_to_f (tlr1, f%taylor(1))
call taylor_to_f (tlr2, f%taylor(2))
call taylor_to_f (tlr3, f%taylor(3))
call taylor_to_f (tlr4, f%taylor(4))
call taylor_to_f (tlr5, f%taylor(5))
call taylor_to_f (tlr6, f%taylor(6))

if (n_wig == 0) then
  if (associated(f%wig_term)) deallocate (f%wig_term)
else
  if (associated(f%wig_term)) then
    if (size(f%wig_term) /= n_wig) then
      deallocate (f%wig_term)
      allocate (f%wig_term(n_wig))
    endif
  else
    allocate (f%wig_term(n_wig))
  endif
endif

f%value                 = value
f%gen0                  = gen0
f%vec0                  = vec0
f%mat6                  = arr2mat(mat6, 6, 6)
f%c_mat                 = arr2mat(c_mat, 2, 2) 
f%gamma_c               = gamma_c
f%s                     = s
f%ref_time              = ref_t

!if (nr1 == 0 .or. nr2 == 0) then
  if (associated(f%r)) deallocate (f%r)
!else
!  if (.not. associated(f%r)) then
!    allocate (f%r(nr1, nr2))
!  elseif ((size(f%r, 1) /= nr1) .or. (size(f%r, 2) /= nr2)) then
!    deallocate (f%r)
!    allocate (f%r(nr1, nr2))
!  endif
!  f%r =  arr2mat(r_arr, nr1, nr2)
!endif

if (n_ab == 0) then
  if (associated (f%a_pole)) deallocate (f%a_pole, f%b_pole)
else
  call multipole_init(f)
  f%a_pole = a_pole
  f%b_pole = b_pole
endif

if (n_const == 0) then
  if (associated (f%const)) deallocate (f%const)
else
  if (associated(f%const)) then
    if (size(f%const) /= n_const) then
      deallocate (f%const)
      allocate (f%const(n_const))
    endif
  else
    allocate (f%const(n_const))
  endif
  f%const = const
endif

f%key                   = key
f%sub_key               = sub_key
f%ix_ele                = ix_ele
f%ix_branch             = ix_branch
f%ix_value              = ix_value

f%slave_status          = slave_status
f%n_slave               = n_slave
f%ix1_slave             = ix1_slave
f%ix2_slave             = ix2_slave

f%lord_status           = lord_status
f%n_lord                = n_lord
f%ic1_lord              = ic1_lord
f%ic2_lord              = ic2_lord

f%ix_pointer            = ix_point
f%ixx                   = ixx

! f%status !

f%mat6_calc_method      = mat6_meth
f%tracking_method       = tracking_meth
f%field_calc            = field_calc
f%ref_orbit             = ref_orb
f%aperture_at           = aperture_at
f%aperture_type         = aperture_type
f%symplectify           = f_logic(symp)
f%mode_flip             = f_logic(mode_flip)
f%multipoles_on         = f_logic(multi_on)
f%scale_multipoles      = f_logic(scale_multi)
f%map_with_offsets      = f_logic(map_with_off)
f%field_master          = f_logic(field_master)
f%reversed              = f_logic(reversed)
f%is_on                 = f_logic(is_on)
f%old_is_on             = f_logic(old_is_on)
f%logic                 = f_logic(logic)
f%bmad_logic            = f_logic(bmad_logic)
f%on_a_girder           = f_logic(girder)
f%csr_calc_on           = f_logic(csr_calc)
f%offset_moves_aperture = f_logic(offset_moves_ap)

! pointer stuff

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine wig_term_in_ele_to_f2 (f_ele, it, coef, kx, ky, kz, phi_z, tp)
!
! Subroutine used by ele_to_f to convert a C++ C_ele into
! a Bmad ele_struct. This routine is not for general use.
!-


subroutine wig_term_in_ele_to_f2 (f_ele, it, coef, kx, ky, kz, phi_z, tp)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (ele_struct) f_ele
real(rp) coef, kx, ky, kz, phi_z
integer it, tp

f_ele%wig_term(it) = wig_term_struct(coef, kx, ky, kz, phi_z, tp)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine mode_info_to_c (f_mode_info, c_mode_info)
!
! Subroutine to convert a Bmad mode_info_struct to a C++ C_mode_info.
!
! Input:
!   f_mode_info -- Mode_info_struct: Input Bmad mode_info_struct.
!
! Output:
!   c_mode_info -- c_dummy_struct: Output C_mode_info.
!-

subroutine mode_info_to_c (f_mode_info, c_mode_info)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (mode_info_struct), target :: f_mode_info
type (mode_info_struct), pointer :: f
type (c_dummy_struct) c_mode_info

f => f_mode_info
call mode_info_to_c2 (c_mode_info, f%tune, f%emit, f%chrom, f%sigma, f%sigmap)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine mode_info_to_f2 (f_mode_info, tune, emit, chrom, sigma, sigamp)
!
! Subroutine used by mode_info_to_f to convert a C++ C_mode_info into
! a Bmad mode_info_struct. This routine is not for general use.
!-

subroutine mode_info_to_f2 (f_mode_info, tune, emit, chrom, sigma, sigmap)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (mode_info_struct) f_mode_info
real(rp) tune, emit, chrom, sigma, sigmap

f_mode_info = mode_info_struct(tune, emit, chrom, sigma, sigmap)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine lat_to_c (f_lat, C_lat)
!
! Subroutine to convert a Bmad lat_struct to a C++ C_lat.
!
! Input:
!   f_lat -- lat_struct: Input Bmad lat_struct.
!
! Output:
!   C_lat -- c_dummy_struct: Output C_lat.
!-

subroutine lat_to_c (f_lat, C_lat)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (lat_struct), target :: f_lat
type (lat_struct), pointer :: f
type (c_dummy_struct) C_lat, c_ele
integer i, n_con, n_ele, n_ic

!

print *, 'LAT_STRUCT CONVERSION BETWEEN C++/FORTRAN NOT YET IMPLEMENTED!'
!!call err_exit

!

f => f_lat

n_con = 0
if (allocated(f%control)) n_con = size(f%control)

n_ele = 0
if (associated(f%ele)) n_ele = size(f%ele)

n_ic = 0
if (allocated(f%ic)) n_ic = size(f%ic)

call lat_to_c2 (C_lat, c_str(f%use_name), c_str(f%lattice), &
      c_str(f%input_file_name), c_str(f%title), f%a, f%b, f%z, f%param, &
      f%version, f%n_ele_track, f%n_ele_max, n_ele, &
      f%n_control_max, f%n_ic_max, f%input_taylor_order, &
      f%ele_init, n_con, f%ic, n_ic)

do i = 1, n_con
  call control_from_lat_to_c2 (C_lat, i, f%control(i))
enddo

do i = 0, f%n_ele_max
  call ele_from_lat_to_c2 (C_lat, i, f%ele(i))
enddo


end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine lat_to_f2 (f, use_name, n_name, lat, n_lat, file, n_file, title, &
!    n_title, x, y, z, param, ver, n_use, n_max, nc_max, n_ic_max, &
!    tlr_ord, ele_init, n_con, ic, n_ic)
!
! Subroutine used by lat_to_f to convert a C++ C_lat into
! a Bmad lat_struct. This routine is not for general use.
!-

subroutine lat_to_f2 (f, use_name, n_name, lat, n_lat, file, n_file, title, &
    n_title, x, y, z, param, ver, n_use, n_max, n_maxx, nc_max, n_ic_max, &
    tlr_ord, ele_init, n_con, ic, n_ic)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (lat_struct) f
type (c_dummy_struct) x, y, z, param, ele_init
integer n_ic_max, n_ic, n_name, n_lat, n_file, n_title, ic(n_ic)
integer ver, n_use, n_max, nc_max, tlr_ord, n_con, n_maxx

character(n_name) use_name
character(n_lat) lat
character(n_file) file
character(n_title) title

!

print *, 'LAT_STRUCT CONVERSION BETWEEN C++/FORTRAN NOT YET IMPLEMENTED!'
!!call err_exit

!

f%use_name         = use_name
f%lattice          = lat
f%input_file_name  = file
f%title            = title
f%version          = ver
f%n_ele_track        = n_use
f%n_ele_track       = n_use
f%n_ele_max        = n_max
f%n_control_max    = nc_max
f%n_ic_max         = n_ic_max
f%input_taylor_order = tlr_ord
call mode_info_to_f (x, f%a)
call mode_info_to_f (y, f%b)
call mode_info_to_f (z, f%z)
call lat_param_to_f (param, f%param)
call ele_to_f (ele_init, f%ele_init)
call allocate_lat_ele_array(f, n_maxx)
allocate (f%control(n_con))
allocate (f%ic(n_ic))
f%ic = ic

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ele_from_lat_to_f2 (f_lat, it, c_ele)
!
! Subroutine used by lat_to_f to convert a C++ C_lat into
! a Bmad lat_struct. This routine is not for general use.
!-

subroutine ele_from_lat_to_f2 (f_lat, it, c_ele)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (lat_struct) f_lat
type (c_dummy_struct) c_ele
integer it

call ele_to_f (c_ele, f_lat%ele(it))

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine control_from_lat_to_f2 (f_lat, it, c_control)
!
! Subroutine used by lat_to_f to convert a C++ C_lat into
! a Bmad lat_struct. This routine is not for general use.
!-

subroutine control_from_lat_to_f2 (f_lat, it, c_control)

use fortran_and_cpp
use bmad_struct
use bmad_interface

implicit none

type (lat_struct) f_lat
type (c_dummy_struct) c_control
integer it

call control_to_f (c_control, f_lat%control(it))

end subroutine
