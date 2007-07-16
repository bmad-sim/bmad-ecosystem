module bmad_and_cpp

use bmad_struct
use bmad_interface

type c_dummy_struct
  real(rp) dummy
end type

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function c_logic (logic) result (c_log)
!
! Function to convert from a fortran logical to a C logical.
!
! Modules needed:
!   use bmad_and_cpp
!
! Input:
!   logic -- Logical: Fortran logical.
!
! Output:
!   c_log -- Integer: C logical.
!-

function c_logic (logic) result (c_log)

implicit none

logical logic
integer c_log

!

if (logic) then
  c_log = 1
else
  c_log = 0
endif

end function

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function f_logic (logic) result (f_log)
!
! Function to convert from a fortran logical to a C logical.
!
! Modules needed:
!   use bmad_and_cpp
!
! Input:
!   logic -- Integer: C logical.
!
! Output:
!   f_log -- Logical: Fortran logical.
!-

function f_logic (logic) result (f_log)

implicit none

integer logic
logical f_log

!

if (logic == 0) then
  f_log = .false.
else
  f_log = .true.
endif

end function

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function r_size (ptr) result (this_size)
!
! Function to return the size of a real pointer.
! If the pointer is not associated then 0 is returned.
!
! Modules needed:
!  use bmad_and_cpp
!
! Input:
!   ptr(:) -- Real(rp), pointer: Pointer to an array.
!
! Output:
!   this_size -- Integer: Size of array. 0 if not associated.
!-

function r_size (ptr) result (this_size)

real(rp), pointer :: ptr(:)
integer this_size

this_size = 0
if (associated(ptr)) this_size = size(ptr)

end function

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function i_size (ptr) result (this_size)
!
! Function to return the size of an integer pointer.
! If the pointer is not associated then 0 is returned.
!
! Modules needed:
!  use bmad_and_cpp
!
! Input:
!   ptr(:) -- Integer, pointer: Pointer to an array.
!
! Output:
!   this_size -- Integer: Size of array. 0 if not associated.
!-


function i_size (ptr) result (this_size)

integer, pointer :: ptr(:)
integer this_size

this_size = 0
if (associated(ptr)) this_size = size(ptr)

end function

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function c_str (str) result (c_string)
!
! Function to append a null (0) character at the end of a string (trimmed
! of trailing blanks) so it will look like a C character array. 
!
! Modules needed:
!  use bmad_and_cpp
!
! Input:
!   str   -- Character(*): Input character string
!
! Output:
!   c_str -- Character(*): String with a null put just after the last
!             non-blank character.
!-

function c_str (str) result (c_string)

character(*) str
character(len_trim(str)+1) c_string

c_string = trim(str) // char(0)

end function

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function mat2arr (mat) result (arr)
!
! Function to take a matrix and turn it into an array:
!   arr(n2*(i-1) + j) = mat(i,j)
! where n2 = size(mat,2).
! This is used for passing matrices to C++ routines.
!
! Modules needed:
!  use bmad_and_cpp
!
! Input:
!   mat(:,:)  -- Real(rp): Input matrix
!
! Output:
!   arr(:)   -- Real(rp): Output array 
!-

function mat2arr (mat) result (arr)

real(rp) mat(:,:)
real(rp) arr(size(mat))
integer i, j, n1, n2

n1 = size(mat, 1); n2 = size(mat, 2)
forall (i = 1:n1, j = 1:n2) arr(n2*(i-1) + j) = mat(i,j)
 
end function

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function arr2mat (arr, n1, n2) result (mat)
!
! Function to take a an array and turn it into a matrix:
!   mat(i,j) = arr(n2*(i-1) + j) 
! This is used for getting matrices from C++ routines.
!
! Modules needed:
!  use bmad_and_cpp
!
! Input:
!   arr(:)   -- Real(rp): Input array.
!   n1       -- Integer: Size of first mat index.
!   n2       -- Integer: Size of second mat index.
!
! Output:
!   mat(n1,n2)  -- Real(rp): Output matrix
!-

function arr2mat (arr, n1, n2) result (mat)

integer i, j, n1, n2
real(rp) arr(:)
real(rp) mat(n1,n2)

forall (i = 1:n1, j = 1:n2) mat(i,j) = arr(n2*(i-1) + j) 
 
end function

end module

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

use bmad_and_cpp

implicit none

type (coord_struct) f_coord
type (c_dummy_struct) c_coord

call coord_to_c2 (c_coord, f_coord%vec)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine coord_to_f2 (f_coord, vec)
!
! Subroutine used by coord_to_f to convert a C++ C_coord into
! a Bmad coord_struct. This routine is not for general use.
!-

subroutine coord_to_f2 (f_coord, vec)

use bmad_and_cpp

implicit none

type (coord_struct) f_coord
real(rp) vec(6)

f_coord%vec = vec

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine orbit_to_c (f_orbit, c_orbit)
!
! Subroutine to convert an a Bmad orbit_struct to a C++ C_orbit.
!
! Input:
!   f_orbit -- Orbit_struct: Input Bmad orbit_struct.
!
! Output:
!   c_orbit -- c_dummy_struct: Output C_orbit.
!-

subroutine orbit_to_c (f_orbit, c_orbit)

use bmad_and_cpp

implicit none

type (orbit_struct) f_orbit
type (c_dummy_struct) c_orbit
integer i

call orbit_to_c2 (c_orbit, size(f_orbit%at))
do i = 0, ubound(f_orbit%at, 1)
  call coord_in_orbit_to_c2 (c_orbit, i, f_orbit%at(i)%vec)
enddo

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine orbit_to_f2 (f_orbit, n_orb)
!
! Subroutine used by orbit_to_f to convert a C++ C_orbit into
! a Bmad orbit_struct. This routine is not for general use.
!-

subroutine orbit_to_f2 (f_orbit, n_orb)

use bmad_and_cpp

implicit none

type (orbit_struct) f_orbit
type (coord_struct) coord0
integer n_orb

if (.not. allocated(f_orbit%at)) then
  allocate (f_orbit%at(0:n_orb-1))
elseif (size(f_orbit%at) < n_orb) then
  deallocate (f_orbit%at)
  allocate (f_orbit%at(0:n_orb-1))
endif

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine coord_in_orbit_to_f2 (f_orbit, ix, vec)
!
! Subroutine used by orbit_to_f to convert a C++ C_orbit into
! a Bmad orbit_struct. This routine is not for general use.
!-

subroutine coord_in_orbit_to_f2 (f_orbit, ix, vec)

use bmad_and_cpp

implicit none

type (orbit_struct) f_orbit
real(rp) vec(6)
integer ix

f_orbit%at(ix)%vec = vec

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

use bmad_and_cpp

implicit none

type (twiss_struct), target :: f_twiss
type (twiss_struct), pointer :: f
type (c_dummy_struct) c_twiss

f => f_twiss
call twiss_to_c2 (c_twiss, f%beta, f%alpha, f%gamma, &
                                      f%phi, f%eta, f%etap, f%sigma, f%emit)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine twiss_to_f2 (f_twiss, beta, alpha, gamma, phi, eta, etap, sigma, emit)
!
! Subroutine used by twiss_to_f to convert a C++ C_twiss into
! a Bmad twiss_struct. This routine is not for general use.
!-

subroutine twiss_to_f2 (f_twiss, beta, alpha, gamma, phi, eta, etap, sigma, emit)

use bmad_and_cpp

implicit none

type (twiss_struct) f_twiss
real(rp) beta, alpha, gamma, phi, eta, etap, sigma, emit

f_twiss = twiss_struct(beta, alpha, gamma, phi, eta, etap, sigma, emit)

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

use bmad_and_cpp

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

use bmad_and_cpp

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

use bmad_and_cpp

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

use bmad_and_cpp

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

use bmad_and_cpp

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

use bmad_and_cpp

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

use bmad_and_cpp

implicit none

type (taylor_term_struct), target :: f_taylor_term
type (taylor_term_struct), pointer :: f
type (c_dummy_struct) c_taylor_term

f => f_taylor_term
call taylor_term_to_c2 (c_taylor_term, f%coef, f%exp)

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

use bmad_and_cpp

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

use bmad_and_cpp

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
  call taylor_term_in_taylor_to_c2 (c_taylor, i, f%term(i)%coef, f%term(i)%exp)
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

use bmad_and_cpp

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

use bmad_and_cpp

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
! Subroutine sr_table_wake_to_c (f_sr_table_wake, c_sr_table_wake)
!
! Subroutine to convert a Bmad sr_table_wake_struct to a C++ C_sr_table_wake.
!
! Input:
!   f_sr_table_wake -- sr_table_wake_struct: Input Bmad sr_table_wake_struct.
!
! Output:
!   c_sr_table_wake -- c_dummy_struct: Output C_sr_table_wake.
!-

subroutine sr_table_wake_to_c (f_sr_table_wake, c_sr_table_wake)

use bmad_and_cpp

implicit none

type (sr_table_wake_struct), target :: f_sr_table_wake
type (sr_table_wake_struct), pointer :: f
type (c_dummy_struct) c_sr_table_wake

f => f_sr_table_wake
call sr_table_wake_to_c2 (c_sr_table_wake, f%z, f%long, f%trans)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine sr_table_wake_to_f2 (f_sr_table_wake, z, long, trans)
!
! Subroutine used by sr_table_wake_to_f to convert a C++ C_sr_table_wake into
! a Bmad sr_table_wake_struct. This routine is not for general use.
!-

subroutine sr_table_wake_to_f2 (f_sr_table_wake, z, long, trans)

use bmad_and_cpp

implicit none

type (sr_table_wake_struct) f_sr_table_wake
real(rp) z, long, trans

f_sr_table_wake = sr_table_wake_struct(z, long, trans)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine sr_mode_wake_to_c (f_sr_mode_wake, c_sr_mode_wake)
!
! Subroutine to convert a Bmad sr_mode_wake_struct to a C++ C_sr_mode_wake.
!
! Input:
!   f_sr_mode_wake -- sr_mode_wake_struct: Input Bmad sr_mode_wake_struct.
!
! Output:
!   c_sr_mode_wake -- c_dummy_struct: Output C_sr_mode_wake.
!-

subroutine sr_mode_wake_to_c (f_sr_mode_wake, c_sr_mode_wake)

use bmad_and_cpp

implicit none

type (sr_mode_wake_struct), target :: f_sr_mode_wake
type (sr_mode_wake_struct), pointer :: f
type (c_dummy_struct) c_sr_mode_wake

f => f_sr_mode_wake
call sr_mode_wake_to_c2 (c_sr_mode_wake, f%amp, f%damp, f%k, f%phi, &
                          f%norm_sin, f%norm_cos, f%skew_sin, f%skew_cos)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine sr_mode_wake_to_f2 (f_sr_mode_wake, amp, damp, freq, phi, 
!                                     norm_sin, norm_cos, skew_sin, skew_cos)
!
! Subroutine used by sr_mode_wake_to_f to convert a C++ C_sr_mode_wake into
! a Bmad sr_mode_wake_struct. This routine is not for general use.
!-

subroutine sr_mode_wake_to_f2 (f_sr_mode_wake, amp, damp, freq, phi, &
                                        norm_sin, norm_cos, skew_sin, skew_cos)

use bmad_and_cpp

implicit none

type (sr_mode_wake_struct) f_sr_mode_wake
real(rp) amp, damp, freq, phi, norm_sin, norm_cos, skew_sin, skew_cos

f_sr_mode_wake = sr_mode_wake_struct(amp, damp, freq, phi, &
                                      norm_sin, norm_cos, skew_sin, skew_cos)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine lr_wake_to_c (f_lr_wake, c_lr_wake)
!
! Subroutine to convert a Bmad lr_wake_struct to a C++ C_lr_wake.
!
! Input:
!   f_lr_wake -- Lr_wake_struct: Input Bmad lr_wake_struct.
!
! Output:
!   c_lr_wake -- c_dummy_struct: Output C_lr_wake.
!-

subroutine lr_wake_to_c (f_lr_wake, c_lr_wake)

use bmad_and_cpp

implicit none

type (lr_wake_struct), target :: f_lr_wake
type (lr_wake_struct), pointer :: f
type (c_dummy_struct) c_lr_wake

f => f_lr_wake
call lr_wake_to_c2 (c_lr_wake, f%freq, f%freq_in, f%R_over_Q, f%q, f%angle, &
         f%norm_sin, f%norm_cos, f%skew_sin, f%skew_cos, f%m, f%polarized)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine lr_wake_to_f2 (f_lr_wake, freq, freq_in, r_over_q, q, angle, &
!                                   n_sin, n_cos, s_cos, s_sin, m, polarized)
!
! Subroutine used by lr_wake_to_f to convert a C++ C_lr_wake into
! a Bmad lr_wake_struct. This routine is not for general use.
!-

subroutine lr_wake_to_f2 (f_lr_wake, freq, freq_in, r_over_q, q, angle, &
                                   n_sin, n_cos, s_cos, s_sin, m, polarized)


use bmad_and_cpp

implicit none

type (lr_wake_struct) f_lr_wake
real(rp) freq, freq_in, r_over_q, q, n_sin, n_cos, s_cos, s_sin, s_ref
real(rp) angle
integer m
logical polarized

f_lr_wake = lr_wake_struct(freq, freq_in, r_over_q, q, angle, &
                                   n_sin, n_cos, s_cos, s_sin, m, polarized)

end subroutine


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine wake_to_c (f_wake, c_wake)
!
! Subroutine to convert a Bmad wake_struct to a C++ C_wake.
!
! Input:
!   f_wake -- Wake_struct: Input Bmad wake_struct.
!
! Output:
!   c_wake -- c_dummy_struct: Output C_wake.
!-

subroutine wake_to_c (f_wake, c_wake)

use bmad_and_cpp

implicit none

type (wake_struct), target :: f_wake
type (wake_struct), pointer :: f
type (sr_mode_wake_struct), pointer :: sr_mode
type (c_dummy_struct) c_wake
integer i, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr

!

f => f_wake
n_sr_table       = size(f%sr_table)
n_sr_mode_long  = size(f%sr_mode_long)
n_sr_mode_trans = size(f%sr_mode_trans)
n_lr        = size(f%lr)

call wake_to_c2 (c_wake, c_str(f%sr_file), c_str(f%lr_file), f%z_sr_mode_max, &
                                             n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr)

do i = 0, n_sr_table-1
  call sr_table_wake_in_wake_to_c2 (c_wake, i, f%sr_table(i)%z, f%sr_table(i)%long, f%sr_table(i)%trans)
enddo

do i = 1, n_sr_mode_long
  sr_mode => f%sr_mode_long(i)
  call sr_mode_long_wake_in_wake_to_c2 (c_wake, i, sr_mode%amp, sr_mode%damp, sr_mode%k, &
                  sr_mode%phi, sr_mode%norm_sin, sr_mode%norm_cos, sr_mode%skew_sin, sr_mode%skew_cos)
enddo

do i = 1, n_sr_mode_trans
  sr_mode => f%sr_mode_trans(i)
  call sr_mode_trans_wake_in_wake_to_c2 (c_wake, i, sr_mode%amp, sr_mode%damp, sr_mode%k, &
                  sr_mode%phi, sr_mode%norm_sin, sr_mode%norm_cos, sr_mode%skew_sin, sr_mode%skew_cos)
enddo

do i = 1, n_lr
  call lr_wake_in_wake_to_c2 (c_wake, i, f%lr(i)%freq, f%lr(i)%freq_in, &
         f%lr(i)%r_over_q, f%lr(i)%Q, f%lr(i)%angle, f%lr(i)%norm_sin, &
         f%lr(i)%norm_cos, f%lr(i)%skew_sin, f%lr(i)%skew_cos, f%lr(i)%m, &
         f%lr(i)%polarized)
enddo 

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine wake_to_f2 (f_wake, sr_file, n_srf, lr_file, n_lrf, z_sr_mode_max,
!                                            n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr)
!
! Subroutine used by wake_to_f to convert a C++ C_wake into
! a Bmad wake_struct. This routine is not for general use.
!-

subroutine wake_to_f2 (f_wake, sr_file, n_srf, lr_file, n_lrf, z_sr_mode_max, &
                                         n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr)

use bmad_and_cpp

implicit none

type (wake_struct) f_wake
integer n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr, n_srf, n_lrf
real(rp) z_sr_mode_max
character(n_srf) :: sr_file
character(n_lrf) :: lr_file

!

f_wake%sr_file = sr_file
f_wake%lr_file = lr_file
f_wake%z_sr_mode_max = z_sr_mode_max
allocate (f_wake%sr_table(0:n_sr_table-1))
allocate (f_wake%sr_mode_long(n_sr_mode_long))
allocate (f_wake%sr_mode_trans(n_sr_mode_trans))
allocate (f_wake%lr(n_lr))

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine sr_table_wake_in_wake_to_f2 (f_wake, it, z, long, trans)
!
! Subroutine used by wake_to_f to convert a C++ C_wake into
! a Bmad wake_struct. This routine is not for general use.
!-

subroutine sr_table_wake_in_wake_to_f2 (f_wake, it, z, long, trans)

use bmad_and_cpp

implicit none

type (wake_struct) f_wake
real(rp) z, long, trans
integer it

f_wake%sr_table(it) = sr_table_wake_struct(z, long, trans)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine sr_mode_long_wake_in_wake_to_f2 (f_wake, it, amp, damp, freq, phi, &
!                                        norm_sin, norm_cos, skew_sin, skew_cos)
!
! Subroutine used by wake_to_f to convert a C++ C_wake into
! a Bmad wake_struct. This routine is not for general use.
!-

subroutine sr_mode_long_wake_in_wake_to_f2 (f_wake, it, amp, damp, freq, phi, &
                                        norm_sin, norm_cos, skew_sin, skew_cos)

use bmad_and_cpp

implicit none

type (wake_struct) f_wake
real(rp) amp, damp, freq, phi, norm_sin, norm_cos, skew_sin, skew_cos
integer it

f_wake%sr_mode_long(it) = sr_mode_wake_struct(amp, damp, freq, phi, &
                                        norm_sin, norm_cos, skew_sin, skew_cos)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine sr_mode_trans_wake_in_wake_to_f2 (f_wake, it, amp, damp, freq, phi, &
!                                        norm_sin, norm_cos, skew_sin, skew_cos)
!
! Subroutine used by wake_to_f to convert a C++ C_wake into
! a Bmad wake_struct. This routine is not for general use.
!-

subroutine sr_mode_trans_wake_in_wake_to_f2 (f_wake, it, amp, damp, freq, phi, &
                                        norm_sin, norm_cos, skew_sin, skew_cos)

use bmad_and_cpp

implicit none

type (wake_struct) f_wake
real(rp) amp, damp, freq, phi, norm_sin, norm_cos, skew_sin, skew_cos
integer it

f_wake%sr_mode_trans(it) = sr_mode_wake_struct(amp, damp, freq, phi, &
                                        norm_sin, norm_cos, skew_sin, skew_cos)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine lr_wake_in_wake_to_f2 (f_wake, it, freq, freq_in, r_over_q, q, angle, &
!                                          n_sin, n_cos, s_cos, s_sin, m, polarized)
!
! Subroutine used by wake_to_f to convert a C++ C_wake into
! a Bmad wake_struct. This routine is not for general use.
!-

subroutine lr_wake_in_wake_to_f2 (f_wake, it, freq, freq_in, r_over_q, q, angle, &
                                          n_sin, n_cos, s_cos, s_sin, m, polarized)

use bmad_and_cpp

implicit none

type (wake_struct) f_wake
real(rp) freq, freq_in, r_over_q, q, n_sin, n_cos, s_cos, s_sin, angle
integer it, m
logical polarized

f_wake%lr(it) = lr_wake_struct(freq, freq_in, r_over_q, q, angle, &
                              n_sin, n_cos, s_cos, s_sin, m, polarized)

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

use bmad_and_cpp

implicit none

type (control_struct), target :: f_control
type (control_struct), pointer :: f
type (c_dummy_struct) c_control

f => f_control
call control_to_c2 (c_control, f%coef, f%ix_lord, f%ix_slave, f%ix_attrib)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine control_to_f2 (f_control, coef, ix_lord, ix_slave, ix_attrib)
!
! Subroutine used by control_to_f to convert a C++ C_control into
! a Bmad control_struct. This routine is not for general use.
!-

subroutine control_to_f2 (f_control, coef, ix_lord, ix_slave, ix_attrib)

use bmad_and_cpp

implicit none

type (control_struct) f_control
real(rp) coef
integer ix_lord, ix_slave, ix_attrib

f_control = control_struct(coef, ix_lord, ix_slave, ix_attrib)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine param_to_c (f_param, c_param)
!
! Subroutine to convert a Bmad lat_param_struct to a C++ C_param.
!
! Input:
!   f_param -- lat_param_struct: Input Bmad lat_param_struct.
!
! Output:
!   c_param -- c_dummy_struct: Output C_param.
!-

subroutine param_to_c (f_param, c_param)

use bmad_and_cpp

implicit none

type (lat_param_struct), target :: f_param
type (lat_param_struct), pointer :: f
type (c_dummy_struct) c_param

f => f_param

call param_to_c2 (c_param, f%n_part, f%total_length, f%growth_rate, &
      mat2arr(f%t1_with_RF), mat2arr(f%t1_no_RF), &
      f%particle, f%ix_lost, f%end_lost_at, f%lattice_type, &
      f%ixx, c_logic(f%stable), c_logic(f%aperture_limit_on), c_logic(f%lost))

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine param_to_f2 (f_param, n_part, total_length, &
!      growth_rate, m1, m2, particle, ix_lost, end_lost_at, &
!      lat_type, ixx, stable, ap_limit_on, lost)
!
! Subroutine used by param_to_f to convert a C++ C_param into
! a Bmad lat_param_struct. This routine is not for general use.
!-

subroutine param_to_f2 (f_param, n_part, total_length, &
      growth_rate, m1, m2, particle, ix_lost, end_lost_at, &
      lat_type, ixx, stable, ap_limit_on, lost) 

use bmad_and_cpp

implicit none

type (lat_param_struct) f_param
real(rp) n_part, total_length, growth_rate
real(rp) m1(36), m2(36)
integer particle, ix_lost, end_lost_at, lat_type, ixx, stable, &
        ap_limit_on, lost

f_param = lat_param_struct(n_part, total_length, growth_rate, &
      arr2mat(m1, 6, 6), arr2mat(m2, 6, 6), particle, ix_lost, end_lost_at, &
      lat_type, ixx, f_logic(stable), f_logic(ap_limit_on), &
      f_logic(lost))

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine amode_to_c (f_amode, c_amode)
!
! Subroutine to convert a Bmad anormal_mode_struct to a C++ C_amode.
!
! Input:
!   f_amode -- Anormal_mode_struct: Input Bmad anormal_mode_struct.
!
! Output:
!   c_amode -- c_dummy_struct: Output C_amode.
!-

subroutine amode_to_c (f_amode, c_amode)

use bmad_and_cpp

implicit none

type (anormal_mode_struct), target :: f_amode
type (anormal_mode_struct), pointer :: f
type (c_dummy_struct) c_amode

f => f_amode
call amode_to_c2 (c_amode, f%emittance, f%synch_int(4), f%synch_int(5), &
                                    f%j_damp, f%alpha_damp, f%chrom, f%tune)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine amode_to_f2 (f_amode, emit, i4, i5, j_damp, a_damp, chrom, tune)
!
! Subroutine used by amode_to_f to convert a C++ C_amode into
! a Bmad anormal_mode_struct. This routine is not for general use.
!-

subroutine amode_to_f2 (f_amode, emit, i4, i5, j_damp, a_damp, chrom, tune)

use bmad_and_cpp

implicit none

type (anormal_mode_struct) f_amode
real(rp) emit, i4, i5, j_damp, a_damp, chrom, tune

f_amode = anormal_mode_struct(emit, (/i4, i5/), j_damp, a_damp, chrom, tune)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine linac_mode_to_c (f_linac_mode, c_linac_mode)
!
! Subroutine to convert a Bmad linac_normal_mode_struct to a C++ C_linac_mode.
!
! Input:
!   f_linac_mode -- Linac_normal_mode_struct: Input Bmad linac_normal_mode_struct.
!
! Output:
!   c_linac_mode -- c_dummy_struct: Output C_linac_mode.
!-

subroutine linac_mode_to_c (f_linac_mode, c_linac_mode)

use bmad_and_cpp

implicit none

type (linac_normal_mode_struct), target :: f_linac_mode
type (linac_normal_mode_struct), pointer :: f
type (c_dummy_struct) c_linac_mode

f => f_linac_mode
call linac_mode_to_c2 (c_linac_mode, f%i2_E4, f%i3_E7, f%i5a_E6, f%i5b_E6, &
                                f%sig_E1, f%a_emittance_end, f%b_emittance_end)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine linac_mode_to_f2 (f_linac_mode, i2, i3, i5a, i5b, sig_e, ea, eb)
!
! Subroutine used by linac_mode_to_f to convert a C++ C_linac_mode into
! a Bmad linac_normal_mode_struct. This routine is not for general use.
!-

subroutine linac_mode_to_f2 (f_linac_mode, i2, i3, i5a, i5b, sig_e, ea, eb)

use bmad_and_cpp

implicit none

type (linac_normal_mode_struct) f_linac_mode
real(rp) i2, i3, i5a, i5b, sig_e, ea, eb

f_linac_mode = linac_normal_mode_struct(i2, i3, i5a, i5b, sig_e, ea, eb)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine modes_to_c (f_modes, c_modes)
!
! Subroutine to convert a Bmad normal_modes_struct to a C++ C_modes.
!
! Input:
!   f_modes -- normal_modes_struct: Input Bmad normal_modes_struct.
!
! Output:
!   c_modes -- c_dummy_struct: Output C_modes.
!-

subroutine modes_to_c (f_modes, c_modes)

use bmad_and_cpp

implicit none

type (normal_modes_struct), target :: f_modes
type (normal_modes_struct), pointer :: f
type (c_dummy_struct) c_modes

f => f_modes
call modes_to_c2 (c_modes, f%synch_int(1), f%synch_int(2), f%synch_int(3), &
            f%sige_e, f%sig_z, f%e_loss, f%pz_aperture, f%a, f%b, f%z, f%lin)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine modes_to_f2 (f_modes, i1, i2, i3, sige, sig_z, e_loss, 
!                                                  pz_aperture, a, b, z, lin)
!
! Subroutine used by modes_to_f to convert a C++ C_modes into
! a Bmad normal_modes_struct. This routine is not for general use.
!-

subroutine modes_to_f2 (f_modes, i1, i2, i3, sige, sig_z, e_loss, &
                                                          pz, a, b, z, lin)

use bmad_and_cpp

implicit none

type (normal_modes_struct) f_modes
type (c_dummy_struct) a, b, z, lin
real(rp) i1, i2, i3, sige, sig_z, e_loss, pz

!

call amode_to_f (a, f_modes%a)
call amode_to_f (b, f_modes%b)
call amode_to_f (z, f_modes%z)
call linac_mode_to_f (lin, f_modes%lin)

f_modes = normal_modes_struct((/i1, i2, i3/), sige, sig_z, e_loss, pz, &
                          f_modes%a, f_modes%b, f_modes%z, f_modes%lin)

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

use bmad_and_cpp

implicit none

type (bmad_common_struct) :: f
type (c_dummy_struct) c_bmad_com

f = bmad_com
call bmad_com_to_c2 (c_bmad_com, &
      f%max_aperture_limit, f%d_orb, f%grad_loss_sr_wake, &
      f%default_ds_step, f%rel_tolerance, f%abs_tolerance, &
      f%taylor_order, f%default_integ_order, &
      c_logic(f%canonical_coords), c_logic(f%use_liar_lcavity), &
      c_logic(f%sr_wakes_on), c_logic(f%lr_wakes_on), &
      c_logic(f%mat6_track_symmetric), c_logic(f%auto_bookkeeper), & 
      c_logic(f%trans_space_charge_on), c_logic(f%coherent_synch_rad_on), &
      c_logic(f%spin_tracking_on), &
      c_logic(f%radiation_damping_on), c_logic(f%radiation_fluctuations_on), &
      c_logic(f%compute_ref_energy))

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine bmad_com_to_f2 (max_ap, orb, grad_loss, ds_step, rel, 
!                abs, taylor_ord, dflt_integ, cc, liar, &
!                sr, lr, sym, a_book, tsc_on, csr_on, st_on, rad_d, rad_f, ref_e)
!
! Subroutine used by bmad_com_to_f to transfer the data from a C++ 
! C_bmad_com variable into the Bmad bmad_com_stuct common block.
! This routine is not for general use.
!-

subroutine bmad_com_to_f2 (max_ap, orb, grad_loss, ds_step, rel, &
                abs, taylor_ord, dflt_integ, cc, liar, &
                sr, lr, sym, a_book, tsc_on, csr_on, st_on, rad_d, rad_f, ref_e)

use bmad_and_cpp

implicit none

real(rp) orb(6), max_ap, grad_loss, rel, abs, ds_step
integer taylor_ord, dflt_integ, cc, liar, sr, lr, sym
integer st_on, rad_d, rad_f, ref_e, a_book, tsc_on, csr_on

bmad_com = bmad_common_struct(max_ap, orb, grad_loss, ds_step, &
    rel, abs, taylor_ord, dflt_integ, &
    f_logic(cc), f_logic(liar), f_logic(sr), f_logic(lr), f_logic(sym), &
    f_logic(a_book), f_logic(tsc_on), f_logic(csr_on), f_logic(st_on), &
    f_logic(rad_d), f_logic(rad_f), f_logic(ref_e))

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

use bmad_and_cpp

implicit none

type (em_field_struct), target :: f_em_field
type (em_field_struct), pointer :: f
type (c_dummy_struct) c_em_field

f => f_em_field

call em_field_to_c2 (c_em_field, f%E, f%b, f%kick, &
            mat2arr(f%dE), mat2arr(f%dB), mat2arr(f%dkick), f%type)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine em_field_to_f2 (f_em_field, e, b, k, de, db, dk, tp)
!
! Subroutine used by em_field_to_f to convert a C++ C_em_field into
! a Bmad em_field_struct. This routine is not for general use.
!-

subroutine em_field_to_f2 (f_em_field, e, b, k, de, db, dk, tp)

use bmad_and_cpp

implicit none

type (em_field_struct) f_em_field
real(rp) e(3), b(3), k(3), de(9), db(9), dk(9)
integer tp

f_em_field = em_field_struct(e, b, k, &
               arr2mat(de, 3, 3), arr2mat(db, 3, 3), arr2mat(dk, 3, 3), tp)

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

use bmad_and_cpp

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
if (associated(f%r)) then
  nr1 = size(f%r, 1)
  nr2 = size(f%r, 2)
  allocate (r_arr(nr1*nr2))
  r_arr = mat2arr (f%r)
else
  allocate (r_arr(0))
endif

descrip =  ' '
if (associated(f%descrip)) descrip = f%descrip

n_wig = 0
if (associated(f%wig_term)) n_wig = size(f%wig_term)

value(0) = 0
value(1:) = f%value


call ele_to_c2 (c_ele, c_str(f%name), c_str(f%type), c_str(f%alias), &
      c_str(f%attribute_name), f%x, f%y, f%a, f%b, f%z, f%floor, value, f%gen0, &
      f%vec0, mat2arr(f%mat6), mat2arr(f%c_mat), f%gamma_c, f%s, &
      r_arr, nr1, nr2, f%a_pole, f%b_pole, r_size(f%a_pole), f%const, r_size(f%const), &
      c_str(descrip), f%gen_field, f%taylor(1), f%taylor(2), f%taylor(3), &
      f%taylor(4), f%taylor(5), f%taylor(6), f%wake, &
      c_logic(associated(f%wake)), n_wig, f%key, &
      f%sub_key, f%control_type, f%ix_value, f%n_slave, f%ix1_slave, &
      f%ix2_slave, f%n_lord, f%ic1_lord, f%ic2_lord, f%ix_pointer, f%ixx, &
      f%ix_ele, f%mat6_calc_method, f%tracking_method, f%field_calc, &
      f%num_steps, f%integrator_order, f%ptc_kind, f%taylor_order, &
      f%aperture_at, f%coupler_at, f%symplectify, f%mode_flip, &
      f%multipoles_on, f%map_with_offsets, &
      f%field_master, f%is_on, f%old_is_on, f%logic, f%on_an_i_beam, &
      f%csr_calc_on)

if (associated(f%r)) deallocate(r_arr)

do i = 1, n_wig
  w => f_ele%wig_term(i)
  call wig_term_in_ele_to_c2 (c_ele, i, &
                      w%coef, w%kx, w%ky, w%kz, w%phi_z, w%type)
end do

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine ele_to_f2 (f, nam, n_nam, typ, n_typ, ali, n_ali, attrib, &
!    n_attrib, x, y, a, b, z, floor, val, g0, v0, m6, c2, gam, s, r_arr, nr1, nr2, &
!    a_pole, b_pole, n_ab, const, n_const, des, n_des, gen, tlr1, tlr2, tlr3, tlr4, &
!    tlr5, tlr6, wake, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr, n_wig, key, sub, &
!    con_tp, ixv, nsl, ix1s, ix2s, nlrd, ic1_l, ic2_l, &
!    ixp, ixx, ixe, m6_meth, tk_meth, f_calc, steps, int_ord, &
!    ptc, tlr_ord, aperture_at, coupler_at, symp, mode, mult, ex_rad,  &
!    f_master, on, intern, logic, i_beam, csr_calc)
!
! Subroutine used by ele_to_f to convert a C++ C_ele into
! a Bmad ele_struct. This routine is not for general use.
!-

subroutine ele_to_f2 (f, nam, n_nam, typ, n_typ, ali, n_ali, attrib, &
    n_attrib, x, y, a, b, z, floor, val, g0, v0, m6, c2, gam, s, r_arr, nr1, nr2, &
    a_pole, b_pole, n_ab, const, n_const, des, n_des, gen, tlr1, tlr2, tlr3, tlr4, &
    tlr5, tlr6, wake, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr, n_wig, key, sub, &
    con_tp, ixv, nsl, ix1s, ix2s, &
    nlrd, ic1_l, ic2_l, ixp, ixx, ixe, m6_meth, tk_meth, f_calc, steps, &
    int_ord, ptc, tlr_ord, aperture_at, coupler_at, symp, mode, mult, ex_rad, &
    f_master, on, intern, logic, i_beam, csr_calc)   

use bmad_and_cpp
use multipole_mod

implicit none

type (ele_struct) f
type (c_dummy_struct) a, b, x, y, z, floor, wake, wig
type (c_dummy_struct) tlr1, tlr2, tlr3, tlr4, tlr5, tlr6
type (genfield), target :: gen

integer n_nam, nr1, nr2, n_ab, n_const, key, sub, con_tp, ixv, nsl, ix1s, &
    ix2s, nlrd, ic1_l, ic2_l, ixp, ixx, ixe, m6_meth, tk_meth, f_calc, steps, &
    int_ord, ptc, tlr_ord, aperture_at, symp, mode, mult, ex_rad, f_master, &
    on, intern, logic, i_beam, csr_calc, n_typ, n_ali, n_attrib, n_des, &
    n_wig, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr, coupler_at

real(rp) val(n_attrib_maxx), g0(6), v0(6), m6(36), c2(4), gam, s
real(rp) a_pole(n_ab), b_pole(n_ab), r_arr(nr1*nr2), const(n_const)

character(n_nam)  nam
character(n_typ)  typ
character(n_ali)  ali
character(n_attrib) attrib
character(n_des)  des

!

f%name  = nam
f%type  = typ
f%alias = ali
f%attribute_name = attrib

call twiss_to_f (a, f%a)
call twiss_to_f (b, f%b)
call twiss_to_f (z, f%z)
call xy_disp_to_f (x, f%x)
call xy_disp_to_f (y, f%y)
call floor_position_to_f (floor, f%floor)

f%value               = val
f%gen0                = g0
f%vec0                = v0
f%mat6                = arr2mat(m6, 6, 6)
f%c_mat               = arr2mat(c2, 2, 2) 
f%gamma_c             = gam
f%s                   = s
f%key                 = key
f%sub_key             = sub
f%control_type        = con_tp
f%ix_value            = ixv
f%n_slave             = nsl
f%ix1_slave           = ix1s
f%ix2_slave           = ix2s
f%n_lord              = nlrd
f%ic1_lord            = ic1_l
f%ic2_lord            = ic2_l
f%ix_pointer          = ixp
f%ixx                 = ixx
f%ix_ele              = ixe
f%mat6_calc_method    = m6_meth
f%tracking_method     = tk_meth
f%field_calc          = f_calc
f%num_steps           = steps
f%integrator_order     = int_ord
f%ptc_kind            = ptc
f%taylor_order        = tlr_ord
f%aperture_at         = aperture_at
f%coupler_at          = coupler_at
f%symplectify         = f_logic(symp)
f%mode_flip           = f_logic(mode)
f%multipoles_on       = f_logic(mult)
f%map_with_offsets  = f_logic(ex_rad)
f%field_master        = f_logic(f_master)
f%is_on               = f_logic(on)
f%old_is_on      = f_logic(intern)
f%logic               = f_logic(logic)
f%on_an_i_beam        = f_logic(i_beam)
f%csr_calc_on         = f_logic(csr_calc)

! pointer stuff

if (nr1 == 0 .or. nr2 == 0) then
  if (associated(f%r)) deallocate (f%r)
else
  if (.not. associated(f%r)) then
    allocate (f%r(nr1, nr2))
  elseif ((size(f%r, 1) /= nr1) .or. (size(f%r, 2) /= nr2)) then
    deallocate (f%r)
    allocate (f%r(nr1, nr2))
  endif
  f%r =  arr2mat(r_arr, nr1, nr2)
endif

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

if (n_des == 0) then
  if (associated (f%descrip)) deallocate (f%descrip)
else
  if (.not. associated(f%descrip)) allocate (f%descrip)
  f%descrip = des
endif

f%gen_field => gen

call taylor_to_f (tlr1, f%taylor(1))
call taylor_to_f (tlr2, f%taylor(2))
call taylor_to_f (tlr3, f%taylor(3))
call taylor_to_f (tlr4, f%taylor(4))
call taylor_to_f (tlr5, f%taylor(5))
call taylor_to_f (tlr6, f%taylor(6))

call init_wake (f%wake, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr)
if (associated(f%wake)) call wake_to_f (wake, f%wake)

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

use bmad_and_cpp

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

use bmad_and_cpp

implicit none

type (mode_info_struct), target :: f_mode_info
type (mode_info_struct), pointer :: f
type (c_dummy_struct) c_mode_info

f => f_mode_info
call mode_info_to_c2 (c_mode_info, f%tune, f%emit, f%chrom)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine mode_info_to_f2 (f_mode_info, tune, emit, chrom)
!
! Subroutine used by mode_info_to_f to convert a C++ C_mode_info into
! a Bmad mode_info_struct. This routine is not for general use.
!-

subroutine mode_info_to_f2 (f_mode_info, tune, emit, chrom)

use bmad_and_cpp

implicit none

type (mode_info_struct) f_mode_info
real(rp) tune, emit, chrom

f_mode_info = mode_info_struct(tune, emit, chrom)

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

use bmad_and_cpp

implicit none

type (lat_struct), target :: f_lat
type (lat_struct), pointer :: f
type (c_dummy_struct) C_lat, c_ele
integer i, n_con, n_ele

!

f => f_lat

n_con = 0
if (associated(f%control)) n_con = size(f%control)

n_ele = 0
if (associated(f%ele)) n_ele = size(f%ele)

call lat_to_c2 (C_lat, c_str(f%name), c_str(f%lattice), &
      c_str(f%input_file_name), c_str(f%title), f%a, f%b, f%z, f%param, &
      f%version, f%n_ele_track, f%n_ele_max, n_ele, &
      f%n_control_max, f%n_ic_max, f%input_taylor_order, &
      f%ele_init, n_con, f%ic, i_size(f%ic))

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
! Subroutine lat_to_f2 (f, name, n_name, lat, n_lat, file, n_file, title, &
!    n_title, x, y, z, param, ver, n_use, n_max, nc_max, n_ic_max, &
!    tlr_ord, ele_init, n_con, ic, n_ic)
!
! Subroutine used by lat_to_f to convert a C++ C_lat into
! a Bmad lat_struct. This routine is not for general use.
!-

subroutine lat_to_f2 (f, name, n_name, lat, n_lat, file, n_file, title, &
    n_title, x, y, z, param, ver, n_use, n_max, n_maxx, nc_max, n_ic_max, &
    tlr_ord, ele_init, n_con, ic, n_ic)

use bmad_and_cpp

implicit none

type (lat_struct) f
type (c_dummy_struct) x, y, z, param, ele_init
integer n_ic_max, n_ic, n_name, n_lat, n_file, n_title, ic(n_ic)
integer ver, n_use, n_max, nc_max, tlr_ord, n_con, n_maxx

character(n_name) name
character(n_lat) lat
character(n_file) file
character(n_title) title

!

f%name             = name
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
call param_to_f (param, f%param)
call ele_to_f (ele_init, f%ele_init)
call allocate_lat_ele(f, n_maxx)
allocate (f%control(n_con))
allocate (f%ic(n_ic))
f%ic = ic
f%E_TOT => f%ele(0)%value(E_TOT$)

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

use bmad_and_cpp

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

use bmad_and_cpp

implicit none

type (lat_struct) f_lat
type (c_dummy_struct) c_control
integer it

call control_to_f (c_control, f_lat%control(it))

end subroutine



