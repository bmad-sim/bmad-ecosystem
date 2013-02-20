
module bmad_cpp_test_mod

use bmad_cpp_convert_mod
use equality_mod

contains

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_coord (ok)

implicit none

type(coord_struct), target :: f_coord, f2_coord
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_coord (c_coord, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_coord
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_coord_test_pattern (f2_coord, 1)

call test_c_coord(c_loc(f2_coord), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_coord_test_pattern (f_coord, 4)
if (f_coord == f2_coord) then
  print *, 'coord: C side convert C->F: Good'
else
  print *, 'coord: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_coord

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_coord (c_coord, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_coord
type(coord_struct), target :: f_coord, f2_coord
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call coord_to_f (c_coord, c_loc(f_coord))

call set_coord_test_pattern (f2_coord, 2)
if (f_coord == f2_coord) then
  print *, 'coord: F side convert C->F: Good'
else
  print *, 'coord: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_coord_test_pattern (f2_coord, 3)
call coord_to_c (c_loc(f2_coord), c_coord)

end subroutine test2_f_coord

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_coord_test_pattern (F, ix_patt)

implicit none

type(coord_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%vec,1); lb1 = lbound(F%vec,1) - 1
  rhs = 100 + jd1 + 1 + offset
  F%vec(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%s = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%t = rhs
!! f_side.test_pat[complex, 1, NOT]
do jd1 = 1, size(F%spin,1); lb1 = lbound(F%spin,1) - 1
  rhs = 100 + jd1 + 4 + offset
  F%spin(jd1+lb1) = cmplx(rhs, 100+rhs)
enddo
!! f_side.test_pat[real, 0, NOT]
rhs = 5 + offset; F%e_field_x = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 6 + offset; F%e_field_y = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 7 + offset; F%phase_x = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 8 + offset; F%phase_y = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 9 + offset; F%charge = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 10 + offset; F%p0c = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 11 + offset; F%beta = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 12 + offset; F%ix_ele = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 13 + offset; F%state = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 14 + offset; F%location = rhs

end subroutine set_coord_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_coord_array (ok)

implicit none

type(coord_array_struct), target :: f_coord_array, f2_coord_array
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_coord_array (c_coord_array, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_coord_array
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_coord_array_test_pattern (f2_coord_array, 1)

call test_c_coord_array(c_loc(f2_coord_array), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_coord_array_test_pattern (f_coord_array, 4)
if (f_coord_array == f2_coord_array) then
  print *, 'coord_array: C side convert C->F: Good'
else
  print *, 'coord_array: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_coord_array

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_coord_array (c_coord_array, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_coord_array
type(coord_array_struct), target :: f_coord_array, f2_coord_array
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call coord_array_to_f (c_coord_array, c_loc(f_coord_array))

call set_coord_array_test_pattern (f2_coord_array, 2)
if (f_coord_array == f2_coord_array) then
  print *, 'coord_array: F side convert C->F: Good'
else
  print *, 'coord_array: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_coord_array_test_pattern (f2_coord_array, 3)
call coord_array_to_c (c_loc(f2_coord_array), c_coord_array)

end subroutine test2_f_coord_array

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_coord_array_test_pattern (F, ix_patt)

implicit none

type(coord_array_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[type, 1, ALLOC]

if (ix_patt < 3) then
  if (allocated(F%orb)) deallocate (F%orb)
else
  if (.not. allocated(F%orb)) allocate (F%orb(-1:1))
  do jd1 = 1, size(F%orb,1); lb1 = lbound(F%orb,1) - 1
    call set_coord_test_pattern (F%orb(jd1+lb1), ix_patt+jd1)
  enddo
endif

end subroutine set_coord_array_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_bpm_phase_coupling (ok)

implicit none

type(bpm_phase_coupling_struct), target :: f_bpm_phase_coupling, f2_bpm_phase_coupling
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_bpm_phase_coupling (c_bpm_phase_coupling, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_bpm_phase_coupling
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_bpm_phase_coupling_test_pattern (f2_bpm_phase_coupling, 1)

call test_c_bpm_phase_coupling(c_loc(f2_bpm_phase_coupling), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_bpm_phase_coupling_test_pattern (f_bpm_phase_coupling, 4)
if (f_bpm_phase_coupling == f2_bpm_phase_coupling) then
  print *, 'bpm_phase_coupling: C side convert C->F: Good'
else
  print *, 'bpm_phase_coupling: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_bpm_phase_coupling

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_bpm_phase_coupling (c_bpm_phase_coupling, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_bpm_phase_coupling
type(bpm_phase_coupling_struct), target :: f_bpm_phase_coupling, f2_bpm_phase_coupling
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call bpm_phase_coupling_to_f (c_bpm_phase_coupling, c_loc(f_bpm_phase_coupling))

call set_bpm_phase_coupling_test_pattern (f2_bpm_phase_coupling, 2)
if (f_bpm_phase_coupling == f2_bpm_phase_coupling) then
  print *, 'bpm_phase_coupling: F side convert C->F: Good'
else
  print *, 'bpm_phase_coupling: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_bpm_phase_coupling_test_pattern (f2_bpm_phase_coupling, 3)
call bpm_phase_coupling_to_c (c_loc(f2_bpm_phase_coupling), c_bpm_phase_coupling)

end subroutine test2_f_bpm_phase_coupling

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_bpm_phase_coupling_test_pattern (F, ix_patt)

implicit none

type(bpm_phase_coupling_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%k_22a = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%k_12a = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%k_11b = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%k_12b = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 5 + offset; F%cbar22_a = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 6 + offset; F%cbar12_a = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 7 + offset; F%cbar11_b = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 8 + offset; F%cbar12_b = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 9 + offset; F%phi_a = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 10 + offset; F%phi_b = rhs

end subroutine set_bpm_phase_coupling_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_wig_term (ok)

implicit none

type(wig_term_struct), target :: f_wig_term, f2_wig_term
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_wig_term (c_wig_term, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_wig_term
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_wig_term_test_pattern (f2_wig_term, 1)

call test_c_wig_term(c_loc(f2_wig_term), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_wig_term_test_pattern (f_wig_term, 4)
if (f_wig_term == f2_wig_term) then
  print *, 'wig_term: C side convert C->F: Good'
else
  print *, 'wig_term: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_wig_term

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_wig_term (c_wig_term, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_wig_term
type(wig_term_struct), target :: f_wig_term, f2_wig_term
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call wig_term_to_f (c_wig_term, c_loc(f_wig_term))

call set_wig_term_test_pattern (f2_wig_term, 2)
if (f_wig_term == f2_wig_term) then
  print *, 'wig_term: F side convert C->F: Good'
else
  print *, 'wig_term: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_wig_term_test_pattern (f2_wig_term, 3)
call wig_term_to_c (c_loc(f2_wig_term), c_wig_term)

end subroutine test2_f_wig_term

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_wig_term_test_pattern (F, ix_patt)

implicit none

type(wig_term_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%coef = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%kx = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%ky = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%kz = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 5 + offset; F%phi_z = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 6 + offset; F%type = rhs

end subroutine set_wig_term_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_wig (ok)

implicit none

type(wig_struct), target :: f_wig, f2_wig
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_wig (c_wig, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_wig
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_wig_test_pattern (f2_wig, 1)

call test_c_wig(c_loc(f2_wig), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_wig_test_pattern (f_wig, 4)
if (f_wig == f2_wig) then
  print *, 'wig: C side convert C->F: Good'
else
  print *, 'wig: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_wig

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_wig (c_wig, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_wig
type(wig_struct), target :: f_wig, f2_wig
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call wig_to_f (c_wig, c_loc(f_wig))

call set_wig_test_pattern (f2_wig, 2)
if (f_wig == f2_wig) then
  print *, 'wig: F side convert C->F: Good'
else
  print *, 'wig: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_wig_test_pattern (f2_wig, 3)
call wig_to_c (c_loc(f2_wig), c_wig)

end subroutine test2_f_wig

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_wig_test_pattern (F, ix_patt)

implicit none

type(wig_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[integer, 0, NOT]
rhs = 1 + offset; F%n_link = rhs
!! f_side.test_pat[type, 1, ALLOC]

if (ix_patt < 3) then
  if (allocated(F%term)) deallocate (F%term)
else
  if (.not. allocated(F%term)) allocate (F%term(-1:1))
  do jd1 = 1, size(F%term,1); lb1 = lbound(F%term,1) - 1
    call set_wig_term_test_pattern (F%term(jd1+lb1), ix_patt+jd1)
  enddo
endif

end subroutine set_wig_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_rf_wake_sr_table (ok)

implicit none

type(rf_wake_sr_table_struct), target :: f_rf_wake_sr_table, f2_rf_wake_sr_table
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_rf_wake_sr_table (c_rf_wake_sr_table, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_rf_wake_sr_table
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_rf_wake_sr_table_test_pattern (f2_rf_wake_sr_table, 1)

call test_c_rf_wake_sr_table(c_loc(f2_rf_wake_sr_table), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_rf_wake_sr_table_test_pattern (f_rf_wake_sr_table, 4)
if (f_rf_wake_sr_table == f2_rf_wake_sr_table) then
  print *, 'rf_wake_sr_table: C side convert C->F: Good'
else
  print *, 'rf_wake_sr_table: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_rf_wake_sr_table

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_rf_wake_sr_table (c_rf_wake_sr_table, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_rf_wake_sr_table
type(rf_wake_sr_table_struct), target :: f_rf_wake_sr_table, f2_rf_wake_sr_table
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call rf_wake_sr_table_to_f (c_rf_wake_sr_table, c_loc(f_rf_wake_sr_table))

call set_rf_wake_sr_table_test_pattern (f2_rf_wake_sr_table, 2)
if (f_rf_wake_sr_table == f2_rf_wake_sr_table) then
  print *, 'rf_wake_sr_table: F side convert C->F: Good'
else
  print *, 'rf_wake_sr_table: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_rf_wake_sr_table_test_pattern (f2_rf_wake_sr_table, 3)
call rf_wake_sr_table_to_c (c_loc(f2_rf_wake_sr_table), c_rf_wake_sr_table)

end subroutine test2_f_rf_wake_sr_table

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_rf_wake_sr_table_test_pattern (F, ix_patt)

implicit none

type(rf_wake_sr_table_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%z = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%long = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%trans = rhs

end subroutine set_rf_wake_sr_table_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_rf_wake_sr_mode (ok)

implicit none

type(rf_wake_sr_mode_struct), target :: f_rf_wake_sr_mode, f2_rf_wake_sr_mode
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_rf_wake_sr_mode (c_rf_wake_sr_mode, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_rf_wake_sr_mode
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_rf_wake_sr_mode_test_pattern (f2_rf_wake_sr_mode, 1)

call test_c_rf_wake_sr_mode(c_loc(f2_rf_wake_sr_mode), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_rf_wake_sr_mode_test_pattern (f_rf_wake_sr_mode, 4)
if (f_rf_wake_sr_mode == f2_rf_wake_sr_mode) then
  print *, 'rf_wake_sr_mode: C side convert C->F: Good'
else
  print *, 'rf_wake_sr_mode: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_rf_wake_sr_mode

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_rf_wake_sr_mode (c_rf_wake_sr_mode, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_rf_wake_sr_mode
type(rf_wake_sr_mode_struct), target :: f_rf_wake_sr_mode, f2_rf_wake_sr_mode
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call rf_wake_sr_mode_to_f (c_rf_wake_sr_mode, c_loc(f_rf_wake_sr_mode))

call set_rf_wake_sr_mode_test_pattern (f2_rf_wake_sr_mode, 2)
if (f_rf_wake_sr_mode == f2_rf_wake_sr_mode) then
  print *, 'rf_wake_sr_mode: F side convert C->F: Good'
else
  print *, 'rf_wake_sr_mode: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_rf_wake_sr_mode_test_pattern (f2_rf_wake_sr_mode, 3)
call rf_wake_sr_mode_to_c (c_loc(f2_rf_wake_sr_mode), c_rf_wake_sr_mode)

end subroutine test2_f_rf_wake_sr_mode

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_rf_wake_sr_mode_test_pattern (F, ix_patt)

implicit none

type(rf_wake_sr_mode_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%amp = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%damp = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%k = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%phi = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 5 + offset; F%b_sin = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 6 + offset; F%b_cos = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 7 + offset; F%a_sin = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 8 + offset; F%a_cos = rhs

end subroutine set_rf_wake_sr_mode_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_rf_wake_lr (ok)

implicit none

type(rf_wake_lr_struct), target :: f_rf_wake_lr, f2_rf_wake_lr
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_rf_wake_lr (c_rf_wake_lr, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_rf_wake_lr
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_rf_wake_lr_test_pattern (f2_rf_wake_lr, 1)

call test_c_rf_wake_lr(c_loc(f2_rf_wake_lr), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_rf_wake_lr_test_pattern (f_rf_wake_lr, 4)
if (f_rf_wake_lr == f2_rf_wake_lr) then
  print *, 'rf_wake_lr: C side convert C->F: Good'
else
  print *, 'rf_wake_lr: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_rf_wake_lr

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_rf_wake_lr (c_rf_wake_lr, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_rf_wake_lr
type(rf_wake_lr_struct), target :: f_rf_wake_lr, f2_rf_wake_lr
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call rf_wake_lr_to_f (c_rf_wake_lr, c_loc(f_rf_wake_lr))

call set_rf_wake_lr_test_pattern (f2_rf_wake_lr, 2)
if (f_rf_wake_lr == f2_rf_wake_lr) then
  print *, 'rf_wake_lr: F side convert C->F: Good'
else
  print *, 'rf_wake_lr: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_rf_wake_lr_test_pattern (f2_rf_wake_lr, 3)
call rf_wake_lr_to_c (c_loc(f2_rf_wake_lr), c_rf_wake_lr)

end subroutine test2_f_rf_wake_lr

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_rf_wake_lr_test_pattern (F, ix_patt)

implicit none

type(rf_wake_lr_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%freq = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%freq_in = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%r_over_q = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%q = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 5 + offset; F%angle = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 6 + offset; F%b_sin = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 7 + offset; F%b_cos = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 8 + offset; F%a_sin = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 9 + offset; F%a_cos = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 10 + offset; F%t_ref = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 11 + offset; F%m = rhs
!! f_side.test_pat[logical, 0, NOT]
rhs = 12 + offset; F%polarized = (modulo(rhs, 2) == 0)

end subroutine set_rf_wake_lr_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_rf_wake (ok)

implicit none

type(rf_wake_struct), target :: f_rf_wake, f2_rf_wake
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_rf_wake (c_rf_wake, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_rf_wake
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_rf_wake_test_pattern (f2_rf_wake, 1)

call test_c_rf_wake(c_loc(f2_rf_wake), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_rf_wake_test_pattern (f_rf_wake, 4)
if (f_rf_wake == f2_rf_wake) then
  print *, 'rf_wake: C side convert C->F: Good'
else
  print *, 'rf_wake: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_rf_wake

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_rf_wake (c_rf_wake, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_rf_wake
type(rf_wake_struct), target :: f_rf_wake, f2_rf_wake
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call rf_wake_to_f (c_rf_wake, c_loc(f_rf_wake))

call set_rf_wake_test_pattern (f2_rf_wake, 2)
if (f_rf_wake == f2_rf_wake) then
  print *, 'rf_wake: F side convert C->F: Good'
else
  print *, 'rf_wake: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_rf_wake_test_pattern (f2_rf_wake, 3)
call rf_wake_to_c (c_loc(f2_rf_wake), c_rf_wake)

end subroutine test2_f_rf_wake

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_rf_wake_test_pattern (F, ix_patt)

implicit none

type(rf_wake_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[character, 0, NOT]
do jd1 = 1, len(F%sr_file)
  F%sr_file(jd1:jd1) = char(ichar("a") + modulo(100+1+offset+jd1, 26))
enddo
!! f_side.test_pat[character, 0, NOT]
do jd1 = 1, len(F%lr_file)
  F%lr_file(jd1:jd1) = char(ichar("a") + modulo(100+2+offset+jd1, 26))
enddo
!! f_side.test_pat[type, 1, PTR]

if (ix_patt < 3) then
  if (associated(F%sr_table)) deallocate (F%sr_table)
else
  if (.not. associated(F%sr_table)) allocate (F%sr_table(-1:1))
  do jd1 = 1, size(F%sr_table,1); lb1 = lbound(F%sr_table,1) - 1
    call set_rf_wake_sr_table_test_pattern (F%sr_table(jd1+lb1), ix_patt+jd1)
  enddo
endif
!! f_side.test_pat[type, 1, PTR]

if (ix_patt < 3) then
  if (associated(F%sr_mode_long)) deallocate (F%sr_mode_long)
else
  if (.not. associated(F%sr_mode_long)) allocate (F%sr_mode_long(-1:1))
  do jd1 = 1, size(F%sr_mode_long,1); lb1 = lbound(F%sr_mode_long,1) - 1
    call set_rf_wake_sr_mode_test_pattern (F%sr_mode_long(jd1+lb1), ix_patt+jd1)
  enddo
endif
!! f_side.test_pat[type, 1, PTR]

if (ix_patt < 3) then
  if (associated(F%sr_mode_trans)) deallocate (F%sr_mode_trans)
else
  if (.not. associated(F%sr_mode_trans)) allocate (F%sr_mode_trans(-1:1))
  do jd1 = 1, size(F%sr_mode_trans,1); lb1 = lbound(F%sr_mode_trans,1) - 1
    call set_rf_wake_sr_mode_test_pattern (F%sr_mode_trans(jd1+lb1), ix_patt+jd1)
  enddo
endif
!! f_side.test_pat[type, 1, PTR]

if (ix_patt < 3) then
  if (associated(F%lr)) deallocate (F%lr)
else
  if (.not. associated(F%lr)) allocate (F%lr(-1:1))
  do jd1 = 1, size(F%lr,1); lb1 = lbound(F%lr,1) - 1
    call set_rf_wake_lr_test_pattern (F%lr(jd1+lb1), ix_patt+jd1)
  enddo
endif
!! f_side.test_pat[real, 0, NOT]
rhs = 11 + offset; F%z_sr_mode_max = rhs

end subroutine set_rf_wake_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_em_field_map_term (ok)

implicit none

type(em_field_map_term_struct), target :: f_em_field_map_term, f2_em_field_map_term
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_em_field_map_term (c_em_field_map_term, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_em_field_map_term
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_em_field_map_term_test_pattern (f2_em_field_map_term, 1)

call test_c_em_field_map_term(c_loc(f2_em_field_map_term), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_em_field_map_term_test_pattern (f_em_field_map_term, 4)
if (f_em_field_map_term == f2_em_field_map_term) then
  print *, 'em_field_map_term: C side convert C->F: Good'
else
  print *, 'em_field_map_term: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_em_field_map_term

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_em_field_map_term (c_em_field_map_term, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_em_field_map_term
type(em_field_map_term_struct), target :: f_em_field_map_term, f2_em_field_map_term
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call em_field_map_term_to_f (c_em_field_map_term, c_loc(f_em_field_map_term))

call set_em_field_map_term_test_pattern (f2_em_field_map_term, 2)
if (f_em_field_map_term == f2_em_field_map_term) then
  print *, 'em_field_map_term: F side convert C->F: Good'
else
  print *, 'em_field_map_term: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_em_field_map_term_test_pattern (f2_em_field_map_term, 3)
call em_field_map_term_to_c (c_loc(f2_em_field_map_term), c_em_field_map_term)

end subroutine test2_f_em_field_map_term

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_em_field_map_term_test_pattern (F, ix_patt)

implicit none

type(em_field_map_term_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[complex, 0, NOT]
rhs = 1 + offset; F%e_coef = cmplx(rhs, 100+rhs)
!! f_side.test_pat[complex, 0, NOT]
rhs = 2 + offset; F%b_coef = cmplx(rhs, 100+rhs)

end subroutine set_em_field_map_term_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_em_field_map (ok)

implicit none

type(em_field_map_struct), target :: f_em_field_map, f2_em_field_map
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_em_field_map (c_em_field_map, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_em_field_map
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_em_field_map_test_pattern (f2_em_field_map, 1)

call test_c_em_field_map(c_loc(f2_em_field_map), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_em_field_map_test_pattern (f_em_field_map, 4)
if (f_em_field_map == f2_em_field_map) then
  print *, 'em_field_map: C side convert C->F: Good'
else
  print *, 'em_field_map: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_em_field_map

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_em_field_map (c_em_field_map, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_em_field_map
type(em_field_map_struct), target :: f_em_field_map, f2_em_field_map
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call em_field_map_to_f (c_em_field_map, c_loc(f_em_field_map))

call set_em_field_map_test_pattern (f2_em_field_map, 2)
if (f_em_field_map == f2_em_field_map) then
  print *, 'em_field_map: F side convert C->F: Good'
else
  print *, 'em_field_map: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_em_field_map_test_pattern (f2_em_field_map, 3)
call em_field_map_to_c (c_loc(f2_em_field_map), c_em_field_map)

end subroutine test2_f_em_field_map

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_em_field_map_test_pattern (F, ix_patt)

implicit none

type(em_field_map_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[character, 0, NOT]
do jd1 = 1, len(F%file)
  F%file(jd1:jd1) = char(ichar("a") + modulo(100+1+offset+jd1, 26))
enddo
!! f_side.test_pat[integer, 0, NOT]
rhs = 2 + offset; F%n_link = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 3 + offset; F%ele_anchor_pt = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%dz = rhs
!! f_side.test_pat[type, 1, ALLOC]

if (ix_patt < 3) then
  if (allocated(F%term)) deallocate (F%term)
else
  if (.not. allocated(F%term)) allocate (F%term(-1:1))
  do jd1 = 1, size(F%term,1); lb1 = lbound(F%term,1) - 1
    call set_em_field_map_term_test_pattern (F%term(jd1+lb1), ix_patt+jd1)
  enddo
endif

end subroutine set_em_field_map_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_em_field_grid_pt (ok)

implicit none

type(em_field_grid_pt_struct), target :: f_em_field_grid_pt, f2_em_field_grid_pt
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_em_field_grid_pt (c_em_field_grid_pt, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_em_field_grid_pt
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_em_field_grid_pt_test_pattern (f2_em_field_grid_pt, 1)

call test_c_em_field_grid_pt(c_loc(f2_em_field_grid_pt), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_em_field_grid_pt_test_pattern (f_em_field_grid_pt, 4)
if (f_em_field_grid_pt == f2_em_field_grid_pt) then
  print *, 'em_field_grid_pt: C side convert C->F: Good'
else
  print *, 'em_field_grid_pt: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_em_field_grid_pt

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_em_field_grid_pt (c_em_field_grid_pt, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_em_field_grid_pt
type(em_field_grid_pt_struct), target :: f_em_field_grid_pt, f2_em_field_grid_pt
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call em_field_grid_pt_to_f (c_em_field_grid_pt, c_loc(f_em_field_grid_pt))

call set_em_field_grid_pt_test_pattern (f2_em_field_grid_pt, 2)
if (f_em_field_grid_pt == f2_em_field_grid_pt) then
  print *, 'em_field_grid_pt: F side convert C->F: Good'
else
  print *, 'em_field_grid_pt: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_em_field_grid_pt_test_pattern (f2_em_field_grid_pt, 3)
call em_field_grid_pt_to_c (c_loc(f2_em_field_grid_pt), c_em_field_grid_pt)

end subroutine test2_f_em_field_grid_pt

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_em_field_grid_pt_test_pattern (F, ix_patt)

implicit none

type(em_field_grid_pt_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[complex, 1, NOT]
do jd1 = 1, size(F%e,1); lb1 = lbound(F%e,1) - 1
  rhs = 100 + jd1 + 1 + offset
  F%e(jd1+lb1) = cmplx(rhs, 100+rhs)
enddo
!! f_side.test_pat[complex, 1, NOT]
do jd1 = 1, size(F%b,1); lb1 = lbound(F%b,1) - 1
  rhs = 100 + jd1 + 2 + offset
  F%b(jd1+lb1) = cmplx(rhs, 100+rhs)
enddo

end subroutine set_em_field_grid_pt_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_em_field_grid (ok)

implicit none

type(em_field_grid_struct), target :: f_em_field_grid, f2_em_field_grid
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_em_field_grid (c_em_field_grid, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_em_field_grid
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_em_field_grid_test_pattern (f2_em_field_grid, 1)

call test_c_em_field_grid(c_loc(f2_em_field_grid), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_em_field_grid_test_pattern (f_em_field_grid, 4)
if (f_em_field_grid == f2_em_field_grid) then
  print *, 'em_field_grid: C side convert C->F: Good'
else
  print *, 'em_field_grid: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_em_field_grid

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_em_field_grid (c_em_field_grid, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_em_field_grid
type(em_field_grid_struct), target :: f_em_field_grid, f2_em_field_grid
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call em_field_grid_to_f (c_em_field_grid, c_loc(f_em_field_grid))

call set_em_field_grid_test_pattern (f2_em_field_grid, 2)
if (f_em_field_grid == f2_em_field_grid) then
  print *, 'em_field_grid: F side convert C->F: Good'
else
  print *, 'em_field_grid: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_em_field_grid_test_pattern (f2_em_field_grid, 3)
call em_field_grid_to_c (c_loc(f2_em_field_grid), c_em_field_grid)

end subroutine test2_f_em_field_grid

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_em_field_grid_test_pattern (F, ix_patt)

implicit none

type(em_field_grid_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[character, 0, NOT]
do jd1 = 1, len(F%file)
  F%file(jd1:jd1) = char(ichar("a") + modulo(100+1+offset+jd1, 26))
enddo
!! f_side.test_pat[integer, 0, NOT]
rhs = 2 + offset; F%type = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 3 + offset; F%ele_anchor_pt = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 4 + offset; F%n_link = rhs
!! f_side.test_pat[type, 3, ALLOC]
if (ix_patt < 3) then
  if (allocated(F%pt)) deallocate (F%pt)
else
  if (.not. allocated(F%pt)) allocate (F%pt(-1:1, 2, 1))
  do jd1 = 1, size(F%pt,1); lb1 = lbound(F%pt,1) - 1
  do jd2 = 1, size(F%pt,2); lb2 = lbound(F%pt,2) - 1
  do jd3 = 1, size(F%pt,3); lb3 = lbound(F%pt,3) - 1
    call set_em_field_grid_pt_test_pattern (F%pt(jd1+lb1,jd2+lb2,jd3+lb3), ix_patt+jd1+2*jd2+3*jd3)
  enddo
  enddo
  enddo
endif
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%dr,1); lb1 = lbound(F%dr,1) - 1
  rhs = 100 + jd1 + 9 + offset
  F%dr(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%r0,1); lb1 = lbound(F%r0,1) - 1
  rhs = 100 + jd1 + 10 + offset
  F%r0(jd1+lb1) = rhs
enddo

end subroutine set_em_field_grid_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_em_field_mode (ok)

implicit none

type(em_field_mode_struct), target :: f_em_field_mode, f2_em_field_mode
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_em_field_mode (c_em_field_mode, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_em_field_mode
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_em_field_mode_test_pattern (f2_em_field_mode, 1)

call test_c_em_field_mode(c_loc(f2_em_field_mode), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_em_field_mode_test_pattern (f_em_field_mode, 4)
if (f_em_field_mode == f2_em_field_mode) then
  print *, 'em_field_mode: C side convert C->F: Good'
else
  print *, 'em_field_mode: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_em_field_mode

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_em_field_mode (c_em_field_mode, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_em_field_mode
type(em_field_mode_struct), target :: f_em_field_mode, f2_em_field_mode
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call em_field_mode_to_f (c_em_field_mode, c_loc(f_em_field_mode))

call set_em_field_mode_test_pattern (f2_em_field_mode, 2)
if (f_em_field_mode == f2_em_field_mode) then
  print *, 'em_field_mode: F side convert C->F: Good'
else
  print *, 'em_field_mode: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_em_field_mode_test_pattern (f2_em_field_mode, 3)
call em_field_mode_to_c (c_loc(f2_em_field_mode), c_em_field_mode)

end subroutine test2_f_em_field_mode

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_em_field_mode_test_pattern (F, ix_patt)

implicit none

type(em_field_mode_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[integer, 0, NOT]
rhs = 1 + offset; F%m = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 2 + offset; F%harmonic = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%f_damp = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%dphi0_ref = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 5 + offset; F%stored_energy = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 6 + offset; F%phi0_azimuth = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 7 + offset; F%field_scale = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 8 + offset; F%master_scale = rhs
!! f_side.test_pat[type, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%map)) deallocate (F%map)
else
  if (.not. associated(F%map)) allocate (F%map)
  rhs = 9 + offset
  call set_em_field_map_test_pattern (F%map, ix_patt)
endif
!! f_side.test_pat[type, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%grid)) deallocate (F%grid)
else
  if (.not. associated(F%grid)) allocate (F%grid)
  rhs = 11 + offset
  call set_em_field_grid_test_pattern (F%grid, ix_patt)
endif

end subroutine set_em_field_mode_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_em_fields (ok)

implicit none

type(em_fields_struct), target :: f_em_fields, f2_em_fields
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_em_fields (c_em_fields, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_em_fields
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_em_fields_test_pattern (f2_em_fields, 1)

call test_c_em_fields(c_loc(f2_em_fields), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_em_fields_test_pattern (f_em_fields, 4)
if (f_em_fields == f2_em_fields) then
  print *, 'em_fields: C side convert C->F: Good'
else
  print *, 'em_fields: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_em_fields

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_em_fields (c_em_fields, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_em_fields
type(em_fields_struct), target :: f_em_fields, f2_em_fields
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call em_fields_to_f (c_em_fields, c_loc(f_em_fields))

call set_em_fields_test_pattern (f2_em_fields, 2)
if (f_em_fields == f2_em_fields) then
  print *, 'em_fields: F side convert C->F: Good'
else
  print *, 'em_fields: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_em_fields_test_pattern (f2_em_fields, 3)
call em_fields_to_c (c_loc(f2_em_fields), c_em_fields)

end subroutine test2_f_em_fields

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_em_fields_test_pattern (F, ix_patt)

implicit none

type(em_fields_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[type, 1, ALLOC]

if (ix_patt < 3) then
  if (allocated(F%mode)) deallocate (F%mode)
else
  if (.not. allocated(F%mode)) allocate (F%mode(-1:1))
  do jd1 = 1, size(F%mode,1); lb1 = lbound(F%mode,1) - 1
    call set_em_field_mode_test_pattern (F%mode(jd1+lb1), ix_patt+jd1)
  enddo
endif

end subroutine set_em_fields_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_floor_position (ok)

implicit none

type(floor_position_struct), target :: f_floor_position, f2_floor_position
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_floor_position (c_floor_position, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_floor_position
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_floor_position_test_pattern (f2_floor_position, 1)

call test_c_floor_position(c_loc(f2_floor_position), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_floor_position_test_pattern (f_floor_position, 4)
if (f_floor_position == f2_floor_position) then
  print *, 'floor_position: C side convert C->F: Good'
else
  print *, 'floor_position: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_floor_position

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_floor_position (c_floor_position, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_floor_position
type(floor_position_struct), target :: f_floor_position, f2_floor_position
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call floor_position_to_f (c_floor_position, c_loc(f_floor_position))

call set_floor_position_test_pattern (f2_floor_position, 2)
if (f_floor_position == f2_floor_position) then
  print *, 'floor_position: F side convert C->F: Good'
else
  print *, 'floor_position: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_floor_position_test_pattern (f2_floor_position, 3)
call floor_position_to_c (c_loc(f2_floor_position), c_floor_position)

end subroutine test2_f_floor_position

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_floor_position_test_pattern (F, ix_patt)

implicit none

type(floor_position_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%r,1); lb1 = lbound(F%r,1) - 1
  rhs = 100 + jd1 + 1 + offset
  F%r(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%theta = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%phi = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%psi = rhs

end subroutine set_floor_position_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_space_charge (ok)

implicit none

type(space_charge_struct), target :: f_space_charge, f2_space_charge
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_space_charge (c_space_charge, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_space_charge
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_space_charge_test_pattern (f2_space_charge, 1)

call test_c_space_charge(c_loc(f2_space_charge), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_space_charge_test_pattern (f_space_charge, 4)
if (f_space_charge == f2_space_charge) then
  print *, 'space_charge: C side convert C->F: Good'
else
  print *, 'space_charge: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_space_charge

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_space_charge (c_space_charge, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_space_charge
type(space_charge_struct), target :: f_space_charge, f2_space_charge
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call space_charge_to_f (c_space_charge, c_loc(f_space_charge))

call set_space_charge_test_pattern (f2_space_charge, 2)
if (f_space_charge == f2_space_charge) then
  print *, 'space_charge: F side convert C->F: Good'
else
  print *, 'space_charge: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_space_charge_test_pattern (f2_space_charge, 3)
call space_charge_to_c (c_loc(f2_space_charge), c_space_charge)

end subroutine test2_f_space_charge

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_space_charge_test_pattern (F, ix_patt)

implicit none

type(space_charge_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[type, 0, NOT]
call set_coord_test_pattern (F%closed_orb, ix_patt)
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%kick_const = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%sig_x = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%sig_y = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 5 + offset; F%phi = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 6 + offset; F%sin_phi = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 7 + offset; F%cos_phi = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 8 + offset; F%sig_z = rhs

end subroutine set_space_charge_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_xy_disp (ok)

implicit none

type(xy_disp_struct), target :: f_xy_disp, f2_xy_disp
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_xy_disp (c_xy_disp, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_xy_disp
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_xy_disp_test_pattern (f2_xy_disp, 1)

call test_c_xy_disp(c_loc(f2_xy_disp), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_xy_disp_test_pattern (f_xy_disp, 4)
if (f_xy_disp == f2_xy_disp) then
  print *, 'xy_disp: C side convert C->F: Good'
else
  print *, 'xy_disp: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_xy_disp

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_xy_disp (c_xy_disp, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_xy_disp
type(xy_disp_struct), target :: f_xy_disp, f2_xy_disp
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call xy_disp_to_f (c_xy_disp, c_loc(f_xy_disp))

call set_xy_disp_test_pattern (f2_xy_disp, 2)
if (f_xy_disp == f2_xy_disp) then
  print *, 'xy_disp: F side convert C->F: Good'
else
  print *, 'xy_disp: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_xy_disp_test_pattern (f2_xy_disp, 3)
call xy_disp_to_c (c_loc(f2_xy_disp), c_xy_disp)

end subroutine test2_f_xy_disp

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_xy_disp_test_pattern (F, ix_patt)

implicit none

type(xy_disp_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%eta = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%etap = rhs

end subroutine set_xy_disp_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_twiss (ok)

implicit none

type(twiss_struct), target :: f_twiss, f2_twiss
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_twiss (c_twiss, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_twiss
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_twiss_test_pattern (f2_twiss, 1)

call test_c_twiss(c_loc(f2_twiss), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_twiss_test_pattern (f_twiss, 4)
if (f_twiss == f2_twiss) then
  print *, 'twiss: C side convert C->F: Good'
else
  print *, 'twiss: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_twiss

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_twiss (c_twiss, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_twiss
type(twiss_struct), target :: f_twiss, f2_twiss
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call twiss_to_f (c_twiss, c_loc(f_twiss))

call set_twiss_test_pattern (f2_twiss, 2)
if (f_twiss == f2_twiss) then
  print *, 'twiss: F side convert C->F: Good'
else
  print *, 'twiss: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_twiss_test_pattern (f2_twiss, 3)
call twiss_to_c (c_loc(f2_twiss), c_twiss)

end subroutine test2_f_twiss

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_twiss_test_pattern (F, ix_patt)

implicit none

type(twiss_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%beta = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%alpha = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%gamma = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%phi = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 5 + offset; F%eta = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 6 + offset; F%etap = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 7 + offset; F%sigma = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 8 + offset; F%sigma_p = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 9 + offset; F%emit = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 10 + offset; F%norm_emit = rhs

end subroutine set_twiss_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_mode3 (ok)

implicit none

type(mode3_struct), target :: f_mode3, f2_mode3
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_mode3 (c_mode3, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_mode3
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_mode3_test_pattern (f2_mode3, 1)

call test_c_mode3(c_loc(f2_mode3), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_mode3_test_pattern (f_mode3, 4)
if (f_mode3 == f2_mode3) then
  print *, 'mode3: C side convert C->F: Good'
else
  print *, 'mode3: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_mode3

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_mode3 (c_mode3, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_mode3
type(mode3_struct), target :: f_mode3, f2_mode3
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call mode3_to_f (c_mode3, c_loc(f_mode3))

call set_mode3_test_pattern (f2_mode3, 2)
if (f_mode3 == f2_mode3) then
  print *, 'mode3: F side convert C->F: Good'
else
  print *, 'mode3: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_mode3_test_pattern (f2_mode3, 3)
call mode3_to_c (c_loc(f2_mode3), c_mode3)

end subroutine test2_f_mode3

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_mode3_test_pattern (F, ix_patt)

implicit none

type(mode3_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 2, NOT]
do jd1 = 1, size(F%v,1); lb1 = lbound(F%v,1) - 1
do jd2 = 1, size(F%v,2); lb2 = lbound(F%v,2) - 1
  rhs = 100 + jd1 + 10*jd2 + 1 + offset
  F%v(jd1+lb1,jd2+lb2) = rhs
enddo; enddo
!! f_side.test_pat[type, 0, NOT]
call set_twiss_test_pattern (F%a, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_twiss_test_pattern (F%b, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_twiss_test_pattern (F%c, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_twiss_test_pattern (F%x, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_twiss_test_pattern (F%y, ix_patt)

end subroutine set_mode3_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_bookkeeping_state (ok)

implicit none

type(bookkeeping_state_struct), target :: f_bookkeeping_state, f2_bookkeeping_state
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_bookkeeping_state (c_bookkeeping_state, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_bookkeeping_state
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_bookkeeping_state_test_pattern (f2_bookkeeping_state, 1)

call test_c_bookkeeping_state(c_loc(f2_bookkeeping_state), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_bookkeeping_state_test_pattern (f_bookkeeping_state, 4)
if (f_bookkeeping_state == f2_bookkeeping_state) then
  print *, 'bookkeeping_state: C side convert C->F: Good'
else
  print *, 'bookkeeping_state: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_bookkeeping_state

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_bookkeeping_state (c_bookkeeping_state, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_bookkeeping_state
type(bookkeeping_state_struct), target :: f_bookkeeping_state, f2_bookkeeping_state
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call bookkeeping_state_to_f (c_bookkeeping_state, c_loc(f_bookkeeping_state))

call set_bookkeeping_state_test_pattern (f2_bookkeeping_state, 2)
if (f_bookkeeping_state == f2_bookkeeping_state) then
  print *, 'bookkeeping_state: F side convert C->F: Good'
else
  print *, 'bookkeeping_state: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_bookkeeping_state_test_pattern (f2_bookkeeping_state, 3)
call bookkeeping_state_to_c (c_loc(f2_bookkeeping_state), c_bookkeeping_state)

end subroutine test2_f_bookkeeping_state

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_bookkeeping_state_test_pattern (F, ix_patt)

implicit none

type(bookkeeping_state_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[integer, 0, NOT]
rhs = 1 + offset; F%attributes = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 2 + offset; F%control = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 3 + offset; F%floor_position = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 4 + offset; F%s_position = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 5 + offset; F%ref_energy = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 6 + offset; F%mat6 = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 7 + offset; F%rad_int = rhs

end subroutine set_bookkeeping_state_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_rad_int_ele_cache (ok)

implicit none

type(rad_int_ele_cache_struct), target :: f_rad_int_ele_cache, f2_rad_int_ele_cache
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_rad_int_ele_cache (c_rad_int_ele_cache, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_rad_int_ele_cache
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_rad_int_ele_cache_test_pattern (f2_rad_int_ele_cache, 1)

call test_c_rad_int_ele_cache(c_loc(f2_rad_int_ele_cache), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_rad_int_ele_cache_test_pattern (f_rad_int_ele_cache, 4)
if (f_rad_int_ele_cache == f2_rad_int_ele_cache) then
  print *, 'rad_int_ele_cache: C side convert C->F: Good'
else
  print *, 'rad_int_ele_cache: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_rad_int_ele_cache

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_rad_int_ele_cache (c_rad_int_ele_cache, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_rad_int_ele_cache
type(rad_int_ele_cache_struct), target :: f_rad_int_ele_cache, f2_rad_int_ele_cache
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call rad_int_ele_cache_to_f (c_rad_int_ele_cache, c_loc(f_rad_int_ele_cache))

call set_rad_int_ele_cache_test_pattern (f2_rad_int_ele_cache, 2)
if (f_rad_int_ele_cache == f2_rad_int_ele_cache) then
  print *, 'rad_int_ele_cache: F side convert C->F: Good'
else
  print *, 'rad_int_ele_cache: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_rad_int_ele_cache_test_pattern (f2_rad_int_ele_cache, 3)
call rad_int_ele_cache_to_c (c_loc(f2_rad_int_ele_cache), c_rad_int_ele_cache)

end subroutine test2_f_rad_int_ele_cache

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_rad_int_ele_cache_test_pattern (F, ix_patt)

implicit none

type(rad_int_ele_cache_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%orb0,1); lb1 = lbound(F%orb0,1) - 1
  rhs = 100 + jd1 + 1 + offset
  F%orb0(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%g2_0 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%g3_0 = rhs
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%dg2_dorb,1); lb1 = lbound(F%dg2_dorb,1) - 1
  rhs = 100 + jd1 + 4 + offset
  F%dg2_dorb(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%dg3_dorb,1); lb1 = lbound(F%dg3_dorb,1) - 1
  rhs = 100 + jd1 + 5 + offset
  F%dg3_dorb(jd1+lb1) = rhs
enddo
!! f_side.test_pat[logical, 0, NOT]
rhs = 6 + offset; F%stale = (modulo(rhs, 2) == 0)

end subroutine set_rad_int_ele_cache_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_wall3d_vertex (ok)

implicit none

type(wall3d_vertex_struct), target :: f_wall3d_vertex, f2_wall3d_vertex
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_wall3d_vertex (c_wall3d_vertex, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_wall3d_vertex
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_wall3d_vertex_test_pattern (f2_wall3d_vertex, 1)

call test_c_wall3d_vertex(c_loc(f2_wall3d_vertex), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_wall3d_vertex_test_pattern (f_wall3d_vertex, 4)
if (f_wall3d_vertex == f2_wall3d_vertex) then
  print *, 'wall3d_vertex: C side convert C->F: Good'
else
  print *, 'wall3d_vertex: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_wall3d_vertex

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_wall3d_vertex (c_wall3d_vertex, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_wall3d_vertex
type(wall3d_vertex_struct), target :: f_wall3d_vertex, f2_wall3d_vertex
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call wall3d_vertex_to_f (c_wall3d_vertex, c_loc(f_wall3d_vertex))

call set_wall3d_vertex_test_pattern (f2_wall3d_vertex, 2)
if (f_wall3d_vertex == f2_wall3d_vertex) then
  print *, 'wall3d_vertex: F side convert C->F: Good'
else
  print *, 'wall3d_vertex: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_wall3d_vertex_test_pattern (f2_wall3d_vertex, 3)
call wall3d_vertex_to_c (c_loc(f2_wall3d_vertex), c_wall3d_vertex)

end subroutine test2_f_wall3d_vertex

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_wall3d_vertex_test_pattern (F, ix_patt)

implicit none

type(wall3d_vertex_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%x = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%y = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%radius_x = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%radius_y = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 5 + offset; F%tilt = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 6 + offset; F%angle = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 7 + offset; F%x0 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 8 + offset; F%y0 = rhs

end subroutine set_wall3d_vertex_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_wall3d_section (ok)

implicit none

type(wall3d_section_struct), target :: f_wall3d_section, f2_wall3d_section
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_wall3d_section (c_wall3d_section, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_wall3d_section
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_wall3d_section_test_pattern (f2_wall3d_section, 1)

call test_c_wall3d_section(c_loc(f2_wall3d_section), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_wall3d_section_test_pattern (f_wall3d_section, 4)
if (f_wall3d_section == f2_wall3d_section) then
  print *, 'wall3d_section: C side convert C->F: Good'
else
  print *, 'wall3d_section: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_wall3d_section

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_wall3d_section (c_wall3d_section, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_wall3d_section
type(wall3d_section_struct), target :: f_wall3d_section, f2_wall3d_section
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call wall3d_section_to_f (c_wall3d_section, c_loc(f_wall3d_section))

call set_wall3d_section_test_pattern (f2_wall3d_section, 2)
if (f_wall3d_section == f2_wall3d_section) then
  print *, 'wall3d_section: F side convert C->F: Good'
else
  print *, 'wall3d_section: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_wall3d_section_test_pattern (f2_wall3d_section, 3)
call wall3d_section_to_c (c_loc(f2_wall3d_section), c_wall3d_section)

end subroutine test2_f_wall3d_section

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_wall3d_section_test_pattern (F, ix_patt)

implicit none

type(wall3d_section_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[integer, 0, NOT]
rhs = 1 + offset; F%type = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%s = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 3 + offset; F%n_vertex_input = rhs
!! f_side.test_pat[type, 1, ALLOC]

if (ix_patt < 3) then
  if (allocated(F%v)) deallocate (F%v)
else
  if (.not. allocated(F%v)) allocate (F%v(-1:1))
  do jd1 = 1, size(F%v,1); lb1 = lbound(F%v,1) - 1
    call set_wall3d_vertex_test_pattern (F%v(jd1+lb1), ix_patt+jd1)
  enddo
endif
!! f_side.test_pat[real, 0, NOT]
rhs = 6 + offset; F%x0 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 7 + offset; F%y0 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 8 + offset; F%dx0_ds = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 9 + offset; F%dy0_ds = rhs
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%x0_coef,1); lb1 = lbound(F%x0_coef,1) - 1
  rhs = 100 + jd1 + 10 + offset
  F%x0_coef(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%y0_coef,1); lb1 = lbound(F%y0_coef,1) - 1
  rhs = 100 + jd1 + 11 + offset
  F%y0_coef(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 0, NOT]
rhs = 12 + offset; F%dr_ds = rhs
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%p1_coef,1); lb1 = lbound(F%p1_coef,1) - 1
  rhs = 100 + jd1 + 13 + offset
  F%p1_coef(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%p2_coef,1); lb1 = lbound(F%p2_coef,1) - 1
  rhs = 100 + jd1 + 14 + offset
  F%p2_coef(jd1+lb1) = rhs
enddo

end subroutine set_wall3d_section_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_wall3d_crotch (ok)

implicit none

type(wall3d_crotch_struct), target :: f_wall3d_crotch, f2_wall3d_crotch
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_wall3d_crotch (c_wall3d_crotch, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_wall3d_crotch
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_wall3d_crotch_test_pattern (f2_wall3d_crotch, 1)

call test_c_wall3d_crotch(c_loc(f2_wall3d_crotch), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_wall3d_crotch_test_pattern (f_wall3d_crotch, 4)
if (f_wall3d_crotch == f2_wall3d_crotch) then
  print *, 'wall3d_crotch: C side convert C->F: Good'
else
  print *, 'wall3d_crotch: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_wall3d_crotch

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_wall3d_crotch (c_wall3d_crotch, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_wall3d_crotch
type(wall3d_crotch_struct), target :: f_wall3d_crotch, f2_wall3d_crotch
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call wall3d_crotch_to_f (c_wall3d_crotch, c_loc(f_wall3d_crotch))

call set_wall3d_crotch_test_pattern (f2_wall3d_crotch, 2)
if (f_wall3d_crotch == f2_wall3d_crotch) then
  print *, 'wall3d_crotch: F side convert C->F: Good'
else
  print *, 'wall3d_crotch: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_wall3d_crotch_test_pattern (f2_wall3d_crotch, 3)
call wall3d_crotch_to_c (c_loc(f2_wall3d_crotch), c_wall3d_crotch)

end subroutine test2_f_wall3d_crotch

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_wall3d_crotch_test_pattern (F, ix_patt)

implicit none

type(wall3d_crotch_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[integer, 0, NOT]
rhs = 1 + offset; F%location = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 2 + offset; F%ix_section = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 3 + offset; F%ix_v1_cut = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 4 + offset; F%ix_v2_cut = rhs
!! f_side.test_pat[type, 0, NOT]
call set_wall3d_section_test_pattern (F%section, ix_patt)

end subroutine set_wall3d_crotch_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_wall3d (ok)

implicit none

type(wall3d_struct), target :: f_wall3d, f2_wall3d
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_wall3d (c_wall3d, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_wall3d
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_wall3d_test_pattern (f2_wall3d, 1)

call test_c_wall3d(c_loc(f2_wall3d), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_wall3d_test_pattern (f_wall3d, 4)
if (f_wall3d == f2_wall3d) then
  print *, 'wall3d: C side convert C->F: Good'
else
  print *, 'wall3d: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_wall3d

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_wall3d (c_wall3d, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_wall3d
type(wall3d_struct), target :: f_wall3d, f2_wall3d
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call wall3d_to_f (c_wall3d, c_loc(f_wall3d))

call set_wall3d_test_pattern (f2_wall3d, 2)
if (f_wall3d == f2_wall3d) then
  print *, 'wall3d: F side convert C->F: Good'
else
  print *, 'wall3d: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_wall3d_test_pattern (f2_wall3d, 3)
call wall3d_to_c (c_loc(f2_wall3d), c_wall3d)

end subroutine test2_f_wall3d

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_wall3d_test_pattern (F, ix_patt)

implicit none

type(wall3d_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[integer, 0, NOT]
rhs = 1 + offset; F%n_link = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 2 + offset; F%priority = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 3 + offset; F%ele_anchor_pt = rhs
!! f_side.test_pat[type, 0, NOT]
call set_wall3d_crotch_test_pattern (F%crotch, ix_patt)
!! f_side.test_pat[type, 1, ALLOC]

if (ix_patt < 3) then
  if (allocated(F%section)) deallocate (F%section)
else
  if (.not. allocated(F%section)) allocate (F%section(-1:1))
  do jd1 = 1, size(F%section,1); lb1 = lbound(F%section,1) - 1
    call set_wall3d_section_test_pattern (F%section(jd1+lb1), ix_patt+jd1)
  enddo
endif

end subroutine set_wall3d_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_taylor_term (ok)

implicit none

type(taylor_term_struct), target :: f_taylor_term, f2_taylor_term
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_taylor_term (c_taylor_term, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_taylor_term
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_taylor_term_test_pattern (f2_taylor_term, 1)

call test_c_taylor_term(c_loc(f2_taylor_term), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_taylor_term_test_pattern (f_taylor_term, 4)
if (f_taylor_term == f2_taylor_term) then
  print *, 'taylor_term: C side convert C->F: Good'
else
  print *, 'taylor_term: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_taylor_term

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_taylor_term (c_taylor_term, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_taylor_term
type(taylor_term_struct), target :: f_taylor_term, f2_taylor_term
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call taylor_term_to_f (c_taylor_term, c_loc(f_taylor_term))

call set_taylor_term_test_pattern (f2_taylor_term, 2)
if (f_taylor_term == f2_taylor_term) then
  print *, 'taylor_term: F side convert C->F: Good'
else
  print *, 'taylor_term: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_taylor_term_test_pattern (f2_taylor_term, 3)
call taylor_term_to_c (c_loc(f2_taylor_term), c_taylor_term)

end subroutine test2_f_taylor_term

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_taylor_term_test_pattern (F, ix_patt)

implicit none

type(taylor_term_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%coef = rhs
!! f_side.test_pat[integer, 1, NOT]
do jd1 = 1, size(F%expn,1); lb1 = lbound(F%expn,1) - 1
  rhs = 100 + jd1 + 2 + offset
  F%expn(jd1+lb1) = rhs
enddo

end subroutine set_taylor_term_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_taylor (ok)

implicit none

type(taylor_struct), target :: f_taylor, f2_taylor
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_taylor (c_taylor, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_taylor
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_taylor_test_pattern (f2_taylor, 1)

call test_c_taylor(c_loc(f2_taylor), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_taylor_test_pattern (f_taylor, 4)
if (f_taylor == f2_taylor) then
  print *, 'taylor: C side convert C->F: Good'
else
  print *, 'taylor: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_taylor

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_taylor (c_taylor, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_taylor
type(taylor_struct), target :: f_taylor, f2_taylor
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call taylor_to_f (c_taylor, c_loc(f_taylor))

call set_taylor_test_pattern (f2_taylor, 2)
if (f_taylor == f2_taylor) then
  print *, 'taylor: F side convert C->F: Good'
else
  print *, 'taylor: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_taylor_test_pattern (f2_taylor, 3)
call taylor_to_c (c_loc(f2_taylor), c_taylor)

end subroutine test2_f_taylor

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_taylor_test_pattern (F, ix_patt)

implicit none

type(taylor_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%ref = rhs
!! f_side.test_pat[type, 1, PTR]

if (ix_patt < 3) then
  if (associated(F%term)) deallocate (F%term)
else
  if (.not. associated(F%term)) allocate (F%term(-1:1))
  do jd1 = 1, size(F%term,1); lb1 = lbound(F%term,1) - 1
    call set_taylor_term_test_pattern (F%term(jd1+lb1), ix_patt+jd1)
  enddo
endif

end subroutine set_taylor_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_control (ok)

implicit none

type(control_struct), target :: f_control, f2_control
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_control (c_control, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_control
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_control_test_pattern (f2_control, 1)

call test_c_control(c_loc(f2_control), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_control_test_pattern (f_control, 4)
if (f_control == f2_control) then
  print *, 'control: C side convert C->F: Good'
else
  print *, 'control: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_control

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_control (c_control, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_control
type(control_struct), target :: f_control, f2_control
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call control_to_f (c_control, c_loc(f_control))

call set_control_test_pattern (f2_control, 2)
if (f_control == f2_control) then
  print *, 'control: F side convert C->F: Good'
else
  print *, 'control: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_control_test_pattern (f2_control, 3)
call control_to_c (c_loc(f2_control), c_control)

end subroutine test2_f_control

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_control_test_pattern (F, ix_patt)

implicit none

type(control_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%coef = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 2 + offset; F%ix_lord = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 3 + offset; F%ix_slave = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 4 + offset; F%ix_branch = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 5 + offset; F%ix_attrib = rhs

end subroutine set_control_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_lat_param (ok)

implicit none

type(lat_param_struct), target :: f_lat_param, f2_lat_param
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_lat_param (c_lat_param, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_lat_param
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_lat_param_test_pattern (f2_lat_param, 1)

call test_c_lat_param(c_loc(f2_lat_param), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_lat_param_test_pattern (f_lat_param, 4)
if (f_lat_param == f2_lat_param) then
  print *, 'lat_param: C side convert C->F: Good'
else
  print *, 'lat_param: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_lat_param

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_lat_param (c_lat_param, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_lat_param
type(lat_param_struct), target :: f_lat_param, f2_lat_param
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call lat_param_to_f (c_lat_param, c_loc(f_lat_param))

call set_lat_param_test_pattern (f2_lat_param, 2)
if (f_lat_param == f2_lat_param) then
  print *, 'lat_param: F side convert C->F: Good'
else
  print *, 'lat_param: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_lat_param_test_pattern (f2_lat_param, 3)
call lat_param_to_c (c_loc(f2_lat_param), c_lat_param)

end subroutine test2_f_lat_param

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_lat_param_test_pattern (F, ix_patt)

implicit none

type(lat_param_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%n_part = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%total_length = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%unstable_factor = rhs
!! f_side.test_pat[real, 2, NOT]
do jd1 = 1, size(F%t1_with_rf,1); lb1 = lbound(F%t1_with_rf,1) - 1
do jd2 = 1, size(F%t1_with_rf,2); lb2 = lbound(F%t1_with_rf,2) - 1
  rhs = 100 + jd1 + 10*jd2 + 4 + offset
  F%t1_with_rf(jd1+lb1,jd2+lb2) = rhs
enddo; enddo
!! f_side.test_pat[real, 2, NOT]
do jd1 = 1, size(F%t1_no_rf,1); lb1 = lbound(F%t1_no_rf,1) - 1
do jd2 = 1, size(F%t1_no_rf,2); lb2 = lbound(F%t1_no_rf,2) - 1
  rhs = 100 + jd1 + 10*jd2 + 5 + offset
  F%t1_no_rf(jd1+lb1,jd2+lb2) = rhs
enddo; enddo
!! f_side.test_pat[real, 0, NOT]
rhs = 6 + offset; F%rel_tracking_charge = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 7 + offset; F%particle = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 8 + offset; F%geometry = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 9 + offset; F%ixx = rhs
!! f_side.test_pat[logical, 0, NOT]
rhs = 10 + offset; F%stable = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 11 + offset; F%aperture_limit_on = (modulo(rhs, 2) == 0)
!! f_side.test_pat[type, 0, NOT]
call set_bookkeeping_state_test_pattern (F%bookkeeping_state, ix_patt)

end subroutine set_lat_param_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_mode_info (ok)

implicit none

type(mode_info_struct), target :: f_mode_info, f2_mode_info
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_mode_info (c_mode_info, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_mode_info
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_mode_info_test_pattern (f2_mode_info, 1)

call test_c_mode_info(c_loc(f2_mode_info), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_mode_info_test_pattern (f_mode_info, 4)
if (f_mode_info == f2_mode_info) then
  print *, 'mode_info: C side convert C->F: Good'
else
  print *, 'mode_info: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_mode_info

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_mode_info (c_mode_info, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_mode_info
type(mode_info_struct), target :: f_mode_info, f2_mode_info
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call mode_info_to_f (c_mode_info, c_loc(f_mode_info))

call set_mode_info_test_pattern (f2_mode_info, 2)
if (f_mode_info == f2_mode_info) then
  print *, 'mode_info: F side convert C->F: Good'
else
  print *, 'mode_info: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_mode_info_test_pattern (f2_mode_info, 3)
call mode_info_to_c (c_loc(f2_mode_info), c_mode_info)

end subroutine test2_f_mode_info

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_mode_info_test_pattern (F, ix_patt)

implicit none

type(mode_info_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%tune = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%emit = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%chrom = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%sigma = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 5 + offset; F%sigmap = rhs

end subroutine set_mode_info_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_pre_tracker (ok)

implicit none

type(pre_tracker_struct), target :: f_pre_tracker, f2_pre_tracker
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_pre_tracker (c_pre_tracker, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_pre_tracker
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_pre_tracker_test_pattern (f2_pre_tracker, 1)

call test_c_pre_tracker(c_loc(f2_pre_tracker), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_pre_tracker_test_pattern (f_pre_tracker, 4)
if (f_pre_tracker == f2_pre_tracker) then
  print *, 'pre_tracker: C side convert C->F: Good'
else
  print *, 'pre_tracker: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_pre_tracker

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_pre_tracker (c_pre_tracker, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_pre_tracker
type(pre_tracker_struct), target :: f_pre_tracker, f2_pre_tracker
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call pre_tracker_to_f (c_pre_tracker, c_loc(f_pre_tracker))

call set_pre_tracker_test_pattern (f2_pre_tracker, 2)
if (f_pre_tracker == f2_pre_tracker) then
  print *, 'pre_tracker: F side convert C->F: Good'
else
  print *, 'pre_tracker: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_pre_tracker_test_pattern (f2_pre_tracker, 3)
call pre_tracker_to_c (c_loc(f2_pre_tracker), c_pre_tracker)

end subroutine test2_f_pre_tracker

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_pre_tracker_test_pattern (F, ix_patt)

implicit none

type(pre_tracker_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[integer, 0, NOT]
rhs = 1 + offset; F%who = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 2 + offset; F%ix_ele_start = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 3 + offset; F%ix_ele_end = rhs
!! f_side.test_pat[character, 0, NOT]
do jd1 = 1, len(F%input_file)
  F%input_file(jd1:jd1) = char(ichar("a") + modulo(100+4+offset+jd1, 26))
enddo

end subroutine set_pre_tracker_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_anormal_mode (ok)

implicit none

type(anormal_mode_struct), target :: f_anormal_mode, f2_anormal_mode
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_anormal_mode (c_anormal_mode, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_anormal_mode
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_anormal_mode_test_pattern (f2_anormal_mode, 1)

call test_c_anormal_mode(c_loc(f2_anormal_mode), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_anormal_mode_test_pattern (f_anormal_mode, 4)
if (f_anormal_mode == f2_anormal_mode) then
  print *, 'anormal_mode: C side convert C->F: Good'
else
  print *, 'anormal_mode: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_anormal_mode

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_anormal_mode (c_anormal_mode, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_anormal_mode
type(anormal_mode_struct), target :: f_anormal_mode, f2_anormal_mode
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call anormal_mode_to_f (c_anormal_mode, c_loc(f_anormal_mode))

call set_anormal_mode_test_pattern (f2_anormal_mode, 2)
if (f_anormal_mode == f2_anormal_mode) then
  print *, 'anormal_mode: F side convert C->F: Good'
else
  print *, 'anormal_mode: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_anormal_mode_test_pattern (f2_anormal_mode, 3)
call anormal_mode_to_c (c_loc(f2_anormal_mode), c_anormal_mode)

end subroutine test2_f_anormal_mode

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_anormal_mode_test_pattern (F, ix_patt)

implicit none

type(anormal_mode_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%emittance = rhs
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%synch_int,1); lb1 = lbound(F%synch_int,1) - 1
  rhs = 100 + jd1 + 2 + offset
  F%synch_int(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%j_damp = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%alpha_damp = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 5 + offset; F%chrom = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 6 + offset; F%tune = rhs

end subroutine set_anormal_mode_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_linac_normal_mode (ok)

implicit none

type(linac_normal_mode_struct), target :: f_linac_normal_mode, f2_linac_normal_mode
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_linac_normal_mode (c_linac_normal_mode, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_linac_normal_mode
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_linac_normal_mode_test_pattern (f2_linac_normal_mode, 1)

call test_c_linac_normal_mode(c_loc(f2_linac_normal_mode), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_linac_normal_mode_test_pattern (f_linac_normal_mode, 4)
if (f_linac_normal_mode == f2_linac_normal_mode) then
  print *, 'linac_normal_mode: C side convert C->F: Good'
else
  print *, 'linac_normal_mode: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_linac_normal_mode

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_linac_normal_mode (c_linac_normal_mode, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_linac_normal_mode
type(linac_normal_mode_struct), target :: f_linac_normal_mode, f2_linac_normal_mode
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call linac_normal_mode_to_f (c_linac_normal_mode, c_loc(f_linac_normal_mode))

call set_linac_normal_mode_test_pattern (f2_linac_normal_mode, 2)
if (f_linac_normal_mode == f2_linac_normal_mode) then
  print *, 'linac_normal_mode: F side convert C->F: Good'
else
  print *, 'linac_normal_mode: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_linac_normal_mode_test_pattern (f2_linac_normal_mode, 3)
call linac_normal_mode_to_c (c_loc(f2_linac_normal_mode), c_linac_normal_mode)

end subroutine test2_f_linac_normal_mode

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_linac_normal_mode_test_pattern (F, ix_patt)

implicit none

type(linac_normal_mode_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%i2_e4 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%i3_e7 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%i5a_e6 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%i5b_e6 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 5 + offset; F%sig_e1 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 6 + offset; F%a_emittance_end = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 7 + offset; F%b_emittance_end = rhs

end subroutine set_linac_normal_mode_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_normal_modes (ok)

implicit none

type(normal_modes_struct), target :: f_normal_modes, f2_normal_modes
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_normal_modes (c_normal_modes, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_normal_modes
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_normal_modes_test_pattern (f2_normal_modes, 1)

call test_c_normal_modes(c_loc(f2_normal_modes), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_normal_modes_test_pattern (f_normal_modes, 4)
if (f_normal_modes == f2_normal_modes) then
  print *, 'normal_modes: C side convert C->F: Good'
else
  print *, 'normal_modes: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_normal_modes

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_normal_modes (c_normal_modes, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_normal_modes
type(normal_modes_struct), target :: f_normal_modes, f2_normal_modes
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call normal_modes_to_f (c_normal_modes, c_loc(f_normal_modes))

call set_normal_modes_test_pattern (f2_normal_modes, 2)
if (f_normal_modes == f2_normal_modes) then
  print *, 'normal_modes: F side convert C->F: Good'
else
  print *, 'normal_modes: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_normal_modes_test_pattern (f2_normal_modes, 3)
call normal_modes_to_c (c_loc(f2_normal_modes), c_normal_modes)

end subroutine test2_f_normal_modes

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_normal_modes_test_pattern (F, ix_patt)

implicit none

type(normal_modes_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%synch_int,1); lb1 = lbound(F%synch_int,1) - 1
  rhs = 100 + jd1 + 1 + offset
  F%synch_int(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%sige_e = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%sig_z = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%e_loss = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 5 + offset; F%rf_voltage = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 6 + offset; F%pz_aperture = rhs
!! f_side.test_pat[type, 0, NOT]
call set_anormal_mode_test_pattern (F%a, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_anormal_mode_test_pattern (F%b, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_anormal_mode_test_pattern (F%z, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_linac_normal_mode_test_pattern (F%lin, ix_patt)

end subroutine set_normal_modes_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_em_field (ok)

implicit none

type(em_field_struct), target :: f_em_field, f2_em_field
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_em_field (c_em_field, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_em_field
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_em_field_test_pattern (f2_em_field, 1)

call test_c_em_field(c_loc(f2_em_field), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_em_field_test_pattern (f_em_field, 4)
if (f_em_field == f2_em_field) then
  print *, 'em_field: C side convert C->F: Good'
else
  print *, 'em_field: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_em_field

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_em_field (c_em_field, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_em_field
type(em_field_struct), target :: f_em_field, f2_em_field
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call em_field_to_f (c_em_field, c_loc(f_em_field))

call set_em_field_test_pattern (f2_em_field, 2)
if (f_em_field == f2_em_field) then
  print *, 'em_field: F side convert C->F: Good'
else
  print *, 'em_field: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_em_field_test_pattern (f2_em_field, 3)
call em_field_to_c (c_loc(f2_em_field), c_em_field)

end subroutine test2_f_em_field

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_em_field_test_pattern (F, ix_patt)

implicit none

type(em_field_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%e,1); lb1 = lbound(F%e,1) - 1
  rhs = 100 + jd1 + 1 + offset
  F%e(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%b,1); lb1 = lbound(F%b,1) - 1
  rhs = 100 + jd1 + 2 + offset
  F%b(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 2, NOT]
do jd1 = 1, size(F%de,1); lb1 = lbound(F%de,1) - 1
do jd2 = 1, size(F%de,2); lb2 = lbound(F%de,2) - 1
  rhs = 100 + jd1 + 10*jd2 + 3 + offset
  F%de(jd1+lb1,jd2+lb2) = rhs
enddo; enddo
!! f_side.test_pat[real, 2, NOT]
do jd1 = 1, size(F%db,1); lb1 = lbound(F%db,1) - 1
do jd2 = 1, size(F%db,2); lb2 = lbound(F%db,2) - 1
  rhs = 100 + jd1 + 10*jd2 + 4 + offset
  F%db(jd1+lb1,jd2+lb2) = rhs
enddo; enddo

end subroutine set_em_field_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_track_map (ok)

implicit none

type(track_map_struct), target :: f_track_map, f2_track_map
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_track_map (c_track_map, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_track_map
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_track_map_test_pattern (f2_track_map, 1)

call test_c_track_map(c_loc(f2_track_map), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_track_map_test_pattern (f_track_map, 4)
if (f_track_map == f2_track_map) then
  print *, 'track_map: C side convert C->F: Good'
else
  print *, 'track_map: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_track_map

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_track_map (c_track_map, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_track_map
type(track_map_struct), target :: f_track_map, f2_track_map
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call track_map_to_f (c_track_map, c_loc(f_track_map))

call set_track_map_test_pattern (f2_track_map, 2)
if (f_track_map == f2_track_map) then
  print *, 'track_map: F side convert C->F: Good'
else
  print *, 'track_map: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_track_map_test_pattern (f2_track_map, 3)
call track_map_to_c (c_loc(f2_track_map), c_track_map)

end subroutine test2_f_track_map

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_track_map_test_pattern (F, ix_patt)

implicit none

type(track_map_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%vec0,1); lb1 = lbound(F%vec0,1) - 1
  rhs = 100 + jd1 + 1 + offset
  F%vec0(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 2, NOT]
do jd1 = 1, size(F%mat6,1); lb1 = lbound(F%mat6,1) - 1
do jd2 = 1, size(F%mat6,2); lb2 = lbound(F%mat6,2) - 1
  rhs = 100 + jd1 + 10*jd2 + 2 + offset
  F%mat6(jd1+lb1,jd2+lb2) = rhs
enddo; enddo

end subroutine set_track_map_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_track (ok)

implicit none

type(track_struct), target :: f_track, f2_track
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_track (c_track, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_track
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_track_test_pattern (f2_track, 1)

call test_c_track(c_loc(f2_track), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_track_test_pattern (f_track, 4)
if (f_track == f2_track) then
  print *, 'track: C side convert C->F: Good'
else
  print *, 'track: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_track

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_track (c_track, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_track
type(track_struct), target :: f_track, f2_track
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call track_to_f (c_track, c_loc(f_track))

call set_track_test_pattern (f2_track, 2)
if (f_track == f2_track) then
  print *, 'track: F side convert C->F: Good'
else
  print *, 'track: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_track_test_pattern (f2_track, 3)
call track_to_c (c_loc(f2_track), c_track)

end subroutine test2_f_track

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_track_test_pattern (F, ix_patt)

implicit none

type(track_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[type, 1, ALLOC]

if (ix_patt < 3) then
  if (allocated(F%orb)) deallocate (F%orb)
else
  if (.not. allocated(F%orb)) allocate (F%orb(-1:1))
  do jd1 = 1, size(F%orb,1); lb1 = lbound(F%orb,1) - 1
    call set_coord_test_pattern (F%orb(jd1+lb1), ix_patt+jd1)
  enddo
endif
!! f_side.test_pat[type, 1, ALLOC]

if (ix_patt < 3) then
  if (allocated(F%field)) deallocate (F%field)
else
  if (.not. allocated(F%field)) allocate (F%field(-1:1))
  do jd1 = 1, size(F%field,1); lb1 = lbound(F%field,1) - 1
    call set_em_field_test_pattern (F%field(jd1+lb1), ix_patt+jd1)
  enddo
endif
!! f_side.test_pat[type, 1, ALLOC]

if (ix_patt < 3) then
  if (allocated(F%map)) deallocate (F%map)
else
  if (.not. allocated(F%map)) allocate (F%map(-1:1))
  do jd1 = 1, size(F%map,1); lb1 = lbound(F%map,1) - 1
    call set_track_map_test_pattern (F%map(jd1+lb1), ix_patt+jd1)
  enddo
endif
!! f_side.test_pat[real, 0, NOT]
rhs = 7 + offset; F%ds_save = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 8 + offset; F%n_pt = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 9 + offset; F%n_bad = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 10 + offset; F%n_ok = rhs

end subroutine set_track_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_synch_rad_common (ok)

implicit none

type(synch_rad_common_struct), target :: f_synch_rad_common, f2_synch_rad_common
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_synch_rad_common (c_synch_rad_common, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_synch_rad_common
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_synch_rad_common_test_pattern (f2_synch_rad_common, 1)

call test_c_synch_rad_common(c_loc(f2_synch_rad_common), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_synch_rad_common_test_pattern (f_synch_rad_common, 4)
if (f_synch_rad_common == f2_synch_rad_common) then
  print *, 'synch_rad_common: C side convert C->F: Good'
else
  print *, 'synch_rad_common: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_synch_rad_common

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_synch_rad_common (c_synch_rad_common, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_synch_rad_common
type(synch_rad_common_struct), target :: f_synch_rad_common, f2_synch_rad_common
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call synch_rad_common_to_f (c_synch_rad_common, c_loc(f_synch_rad_common))

call set_synch_rad_common_test_pattern (f2_synch_rad_common, 2)
if (f_synch_rad_common == f2_synch_rad_common) then
  print *, 'synch_rad_common: F side convert C->F: Good'
else
  print *, 'synch_rad_common: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_synch_rad_common_test_pattern (f2_synch_rad_common, 3)
call synch_rad_common_to_c (c_loc(f2_synch_rad_common), c_synch_rad_common)

end subroutine test2_f_synch_rad_common

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_synch_rad_common_test_pattern (F, ix_patt)

implicit none

type(synch_rad_common_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%scale = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%i2 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%i3 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%i5a = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 5 + offset; F%i5b = rhs
!! f_side.test_pat[logical, 0, NOT]
rhs = 6 + offset; F%i_calc_on = (modulo(rhs, 2) == 0)

end subroutine set_synch_rad_common_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_bmad_common (ok)

implicit none

type(bmad_common_struct), target :: f_bmad_common, f2_bmad_common
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_bmad_common (c_bmad_common, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_bmad_common
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_bmad_common_test_pattern (f2_bmad_common, 1)

call test_c_bmad_common(c_loc(f2_bmad_common), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_bmad_common_test_pattern (f_bmad_common, 4)
if (f_bmad_common == f2_bmad_common) then
  print *, 'bmad_common: C side convert C->F: Good'
else
  print *, 'bmad_common: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_bmad_common

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_bmad_common (c_bmad_common, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_bmad_common
type(bmad_common_struct), target :: f_bmad_common, f2_bmad_common
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call bmad_common_to_f (c_bmad_common, c_loc(f_bmad_common))

call set_bmad_common_test_pattern (f2_bmad_common, 2)
if (f_bmad_common == f2_bmad_common) then
  print *, 'bmad_common: F side convert C->F: Good'
else
  print *, 'bmad_common: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_bmad_common_test_pattern (f2_bmad_common, 3)
call bmad_common_to_c (c_loc(f2_bmad_common), c_bmad_common)

end subroutine test2_f_bmad_common

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_bmad_common_test_pattern (F, ix_patt)

implicit none

type(bmad_common_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%max_aperture_limit = rhs
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%d_orb,1); lb1 = lbound(F%d_orb,1) - 1
  rhs = 100 + jd1 + 2 + offset
  F%d_orb(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%default_ds_step = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%significant_length = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 5 + offset; F%rel_tol_tracking = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 6 + offset; F%abs_tol_tracking = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 7 + offset; F%rel_tol_adaptive_tracking = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 8 + offset; F%abs_tol_adaptive_tracking = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 9 + offset; F%init_ds_adaptive_tracking = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 10 + offset; F%min_ds_adaptive_tracking = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 11 + offset; F%taylor_order = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 12 + offset; F%default_integ_order = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 13 + offset; F%ptc_max_fringe_order = rhs
!! f_side.test_pat[logical, 0, NOT]
rhs = 14 + offset; F%use_hard_edge_drifts = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 15 + offset; F%sr_wakes_on = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 16 + offset; F%lr_wakes_on = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 17 + offset; F%mat6_track_symmetric = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 18 + offset; F%auto_bookkeeper = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 19 + offset; F%space_charge_on = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 20 + offset; F%coherent_synch_rad_on = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 21 + offset; F%spin_tracking_on = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 22 + offset; F%radiation_damping_on = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 23 + offset; F%radiation_fluctuations_on = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 24 + offset; F%conserve_taylor_maps = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 25 + offset; F%absolute_time_tracking_default = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 26 + offset; F%rf_auto_scale_phase_default = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 27 + offset; F%rf_auto_scale_amp_default = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 28 + offset; F%use_ptc_layout_default = (modulo(rhs, 2) == 0)

end subroutine set_bmad_common_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_rad_int1 (ok)

implicit none

type(rad_int1_struct), target :: f_rad_int1, f2_rad_int1
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_rad_int1 (c_rad_int1, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_rad_int1
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_rad_int1_test_pattern (f2_rad_int1, 1)

call test_c_rad_int1(c_loc(f2_rad_int1), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_rad_int1_test_pattern (f_rad_int1, 4)
if (f_rad_int1 == f2_rad_int1) then
  print *, 'rad_int1: C side convert C->F: Good'
else
  print *, 'rad_int1: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_rad_int1

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_rad_int1 (c_rad_int1, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_rad_int1
type(rad_int1_struct), target :: f_rad_int1, f2_rad_int1
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call rad_int1_to_f (c_rad_int1, c_loc(f_rad_int1))

call set_rad_int1_test_pattern (f2_rad_int1, 2)
if (f_rad_int1 == f2_rad_int1) then
  print *, 'rad_int1: F side convert C->F: Good'
else
  print *, 'rad_int1: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_rad_int1_test_pattern (f2_rad_int1, 3)
call rad_int1_to_c (c_loc(f2_rad_int1), c_rad_int1)

end subroutine test2_f_rad_int1

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_rad_int1_test_pattern (F, ix_patt)

implicit none

type(rad_int1_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[real, 0, NOT]
rhs = 1 + offset; F%i0 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 2 + offset; F%i1 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 3 + offset; F%i2 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 4 + offset; F%i3 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 5 + offset; F%i4a = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 6 + offset; F%i4b = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 7 + offset; F%i4z = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 8 + offset; F%i5a = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 9 + offset; F%i5b = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 10 + offset; F%i6b = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 11 + offset; F%lin_i2_e4 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 12 + offset; F%lin_i3_e7 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 13 + offset; F%lin_i5a_e6 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 14 + offset; F%lin_i5b_e6 = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 15 + offset; F%lin_norm_emit_a = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 16 + offset; F%lin_norm_emit_b = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 17 + offset; F%n_steps = rhs

end subroutine set_rad_int1_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_rad_int_all_ele (ok)

implicit none

type(rad_int_all_ele_struct), target :: f_rad_int_all_ele, f2_rad_int_all_ele
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_rad_int_all_ele (c_rad_int_all_ele, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_rad_int_all_ele
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_rad_int_all_ele_test_pattern (f2_rad_int_all_ele, 1)

call test_c_rad_int_all_ele(c_loc(f2_rad_int_all_ele), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_rad_int_all_ele_test_pattern (f_rad_int_all_ele, 4)
if (f_rad_int_all_ele == f2_rad_int_all_ele) then
  print *, 'rad_int_all_ele: C side convert C->F: Good'
else
  print *, 'rad_int_all_ele: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_rad_int_all_ele

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_rad_int_all_ele (c_rad_int_all_ele, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_rad_int_all_ele
type(rad_int_all_ele_struct), target :: f_rad_int_all_ele, f2_rad_int_all_ele
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call rad_int_all_ele_to_f (c_rad_int_all_ele, c_loc(f_rad_int_all_ele))

call set_rad_int_all_ele_test_pattern (f2_rad_int_all_ele, 2)
if (f_rad_int_all_ele == f2_rad_int_all_ele) then
  print *, 'rad_int_all_ele: F side convert C->F: Good'
else
  print *, 'rad_int_all_ele: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_rad_int_all_ele_test_pattern (f2_rad_int_all_ele, 3)
call rad_int_all_ele_to_c (c_loc(f2_rad_int_all_ele), c_rad_int_all_ele)

end subroutine test2_f_rad_int_all_ele

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_rad_int_all_ele_test_pattern (F, ix_patt)

implicit none

type(rad_int_all_ele_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[type, 1, ALLOC]

if (ix_patt < 3) then
  if (allocated(F%ele)) deallocate (F%ele)
else
  if (.not. allocated(F%ele)) allocate (F%ele(-1:1))
  do jd1 = 1, size(F%ele,1); lb1 = lbound(F%ele,1) - 1
    call set_rad_int1_test_pattern (F%ele(jd1+lb1), ix_patt+jd1)
  enddo
endif

end subroutine set_rad_int_all_ele_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_ele (ok)

implicit none

type(ele_struct), target :: f_ele, f2_ele
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_ele (c_ele, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_ele
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_ele_test_pattern (f2_ele, 1)

call test_c_ele(c_loc(f2_ele), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_ele_test_pattern (f_ele, 4)
if (f_ele == f2_ele) then
  print *, 'ele: C side convert C->F: Good'
else
  print *, 'ele: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_ele

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_ele (c_ele, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_ele
type(ele_struct), target :: f_ele, f2_ele
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call ele_to_f (c_ele, c_loc(f_ele))

call set_ele_test_pattern (f2_ele, 2)
if (f_ele == f2_ele) then
  print *, 'ele: F side convert C->F: Good'
else
  print *, 'ele: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_ele_test_pattern (f2_ele, 3)
call ele_to_c (c_loc(f2_ele), c_ele)

end subroutine test2_f_ele

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_ele_test_pattern (F, ix_patt)

implicit none

type(ele_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[character, 0, NOT]
do jd1 = 1, len(F%name)
  F%name(jd1:jd1) = char(ichar("a") + modulo(100+1+offset+jd1, 26))
enddo
!! f_side.test_pat[character, 0, NOT]
do jd1 = 1, len(F%type)
  F%type(jd1:jd1) = char(ichar("a") + modulo(100+2+offset+jd1, 26))
enddo
!! f_side.test_pat[character, 0, NOT]
do jd1 = 1, len(F%alias)
  F%alias(jd1:jd1) = char(ichar("a") + modulo(100+3+offset+jd1, 26))
enddo
!! f_side.test_pat[character, 0, NOT]
do jd1 = 1, len(F%component_name)
  F%component_name(jd1:jd1) = char(ichar("a") + modulo(100+4+offset+jd1, 26))
enddo
!! f_side.test_pat[character, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%descrip)) deallocate (F%descrip)
else
  if (.not. associated(F%descrip)) allocate (F%descrip)
  do jd1 = 1, len(F%descrip)
    F%descrip(jd1:jd1) = char(ichar("a") + modulo(100+5+offset+jd1, 26))
  enddo
endif
!! f_side.test_pat[type, 0, NOT]
call set_twiss_test_pattern (F%a, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_twiss_test_pattern (F%b, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_twiss_test_pattern (F%z, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_xy_disp_test_pattern (F%x, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_xy_disp_test_pattern (F%y, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_bookkeeping_state_test_pattern (F%bookkeeping_state, ix_patt)
!! f_side.test_pat[type, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%em_field)) deallocate (F%em_field)
else
  if (.not. associated(F%em_field)) allocate (F%em_field)
  rhs = 13 + offset
  call set_em_fields_test_pattern (F%em_field, ix_patt)
endif
!! f_side.test_pat[type, 0, NOT]
call set_floor_position_test_pattern (F%floor, ix_patt)
!! f_side.test_pat[type, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%mode3)) deallocate (F%mode3)
else
  if (.not. associated(F%mode3)) allocate (F%mode3)
  rhs = 16 + offset
  call set_mode3_test_pattern (F%mode3, ix_patt)
endif
!! f_side.test_pat[type, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%rad_int_cache)) deallocate (F%rad_int_cache)
else
  if (.not. associated(F%rad_int_cache)) allocate (F%rad_int_cache)
  rhs = 18 + offset
  call set_rad_int_ele_cache_test_pattern (F%rad_int_cache, ix_patt)
endif
!! f_side.test_pat[type, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%rf_wake)) deallocate (F%rf_wake)
else
  if (.not. associated(F%rf_wake)) allocate (F%rf_wake)
  rhs = 20 + offset
  call set_rf_wake_test_pattern (F%rf_wake, ix_patt)
endif
!! f_side.test_pat[type, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%space_charge)) deallocate (F%space_charge)
else
  if (.not. associated(F%space_charge)) allocate (F%space_charge)
  rhs = 22 + offset
  call set_space_charge_test_pattern (F%space_charge, ix_patt)
endif
!! f_side.test_pat[type, 1, NOT]
do jd1 = 1, size(F%taylor,1); lb1 = lbound(F%taylor,1) - 1
  rhs = 100 + jd1 + 24 + offset
  call set_taylor_test_pattern (F%taylor(jd1+lb1), ix_patt+jd1)
enddo
!! f_side.test_pat[type, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%wall3d)) deallocate (F%wall3d)
else
  if (.not. associated(F%wall3d)) allocate (F%wall3d)
  rhs = 25 + offset
  call set_wall3d_test_pattern (F%wall3d, ix_patt)
endif
!! f_side.test_pat[type, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%wig)) deallocate (F%wig)
else
  if (.not. associated(F%wig)) allocate (F%wig)
  rhs = 27 + offset
  call set_wig_test_pattern (F%wig, ix_patt)
endif
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%value,1); lb1 = lbound(F%value,1) - 1
  rhs = 100 + jd1 + 29 + offset
  F%value(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%old_value,1); lb1 = lbound(F%old_value,1) - 1
  rhs = 100 + jd1 + 30 + offset
  F%old_value(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%gen0,1); lb1 = lbound(F%gen0,1) - 1
  rhs = 100 + jd1 + 31 + offset
  F%gen0(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%vec0,1); lb1 = lbound(F%vec0,1) - 1
  rhs = 100 + jd1 + 32 + offset
  F%vec0(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 2, NOT]
do jd1 = 1, size(F%mat6,1); lb1 = lbound(F%mat6,1) - 1
do jd2 = 1, size(F%mat6,2); lb2 = lbound(F%mat6,2) - 1
  rhs = 100 + jd1 + 10*jd2 + 33 + offset
  F%mat6(jd1+lb1,jd2+lb2) = rhs
enddo; enddo
!! f_side.test_pat[real, 2, NOT]
do jd1 = 1, size(F%c_mat,1); lb1 = lbound(F%c_mat,1) - 1
do jd2 = 1, size(F%c_mat,2); lb2 = lbound(F%c_mat,2) - 1
  rhs = 100 + jd1 + 10*jd2 + 34 + offset
  F%c_mat(jd1+lb1,jd2+lb2) = rhs
enddo; enddo
!! f_side.test_pat[real, 0, NOT]
rhs = 35 + offset; F%gamma_c = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 36 + offset; F%s = rhs
!! f_side.test_pat[real, 0, NOT]
rhs = 37 + offset; F%ref_time = rhs
!! f_side.test_pat[real, 3, PTR]
if (ix_patt < 3) then
  if (associated(F%r)) deallocate (F%r)
else
  if (.not. associated(F%r)) allocate (F%r(-1:1, 2, 1))
  do jd1 = 1, size(F%r,1); lb1 = lbound(F%r,1) - 1
  do jd2 = 1, size(F%r,2); lb2 = lbound(F%r,2) - 1
  do jd3 = 1, size(F%r,3); lb3 = lbound(F%r,3) - 1
    rhs = 100 + jd1 + 10*jd2 + 100*jd3 + 38 + offset
    F%r(jd1+lb1,jd2+lb2,jd3+lb3) = rhs
  enddo; enddo; enddo
endif
!! f_side.test_pat[real, 1, PTR]

if (ix_patt < 3) then
  if (associated(F%a_pole)) deallocate (F%a_pole)
else
  if (.not. associated(F%a_pole)) allocate (F%a_pole(-1:1))
  do jd1 = 1, size(F%a_pole,1); lb1 = lbound(F%a_pole,1) - 1
    rhs = 100 + jd1 + 42 + offset
    F%a_pole(jd1+lb1) = rhs
  enddo
endif
!! f_side.test_pat[real, 1, PTR]

if (ix_patt < 3) then
  if (associated(F%b_pole)) deallocate (F%b_pole)
else
  if (.not. associated(F%b_pole)) allocate (F%b_pole(-1:1))
  do jd1 = 1, size(F%b_pole,1); lb1 = lbound(F%b_pole,1) - 1
    rhs = 100 + jd1 + 44 + offset
    F%b_pole(jd1+lb1) = rhs
  enddo
endif
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%map_ref_orb_in,1); lb1 = lbound(F%map_ref_orb_in,1) - 1
  rhs = 100 + jd1 + 46 + offset
  F%map_ref_orb_in(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%map_ref_orb_out,1); lb1 = lbound(F%map_ref_orb_out,1) - 1
  rhs = 100 + jd1 + 47 + offset
  F%map_ref_orb_out(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%time_ref_orb_in,1); lb1 = lbound(F%time_ref_orb_in,1) - 1
  rhs = 100 + jd1 + 48 + offset
  F%time_ref_orb_in(jd1+lb1) = rhs
enddo
!! f_side.test_pat[real, 1, NOT]
do jd1 = 1, size(F%time_ref_orb_out,1); lb1 = lbound(F%time_ref_orb_out,1) - 1
  rhs = 100 + jd1 + 49 + offset
  F%time_ref_orb_out(jd1+lb1) = rhs
enddo
!! f_side.test_pat[integer, 0, NOT]
rhs = 50 + offset; F%key = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 51 + offset; F%sub_key = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 52 + offset; F%ix_ele = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 53 + offset; F%ix_branch = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 54 + offset; F%ix_value = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 55 + offset; F%slave_status = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 56 + offset; F%n_slave = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 57 + offset; F%ix1_slave = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 58 + offset; F%ix2_slave = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 59 + offset; F%lord_status = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 60 + offset; F%n_lord = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 61 + offset; F%ic1_lord = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 62 + offset; F%ic2_lord = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 63 + offset; F%ix_pointer = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 64 + offset; F%ixx = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 65 + offset; F%iyy = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 66 + offset; F%mat6_calc_method = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 67 + offset; F%tracking_method = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 68 + offset; F%spin_tracking_method = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 69 + offset; F%ptc_integration_type = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 70 + offset; F%field_calc = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 71 + offset; F%aperture_at = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 72 + offset; F%aperture_type = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 73 + offset; F%orientation = rhs
!! f_side.test_pat[logical, 0, NOT]
rhs = 74 + offset; F%symplectify = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 75 + offset; F%mode_flip = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 76 + offset; F%multipoles_on = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 77 + offset; F%scale_multipoles = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 78 + offset; F%map_with_offsets = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 79 + offset; F%field_master = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 80 + offset; F%is_on = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 81 + offset; F%old_is_on = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 82 + offset; F%logic = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 83 + offset; F%bmad_logic = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 84 + offset; F%on_a_girder = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 85 + offset; F%csr_calc_on = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 86 + offset; F%offset_moves_aperture = (modulo(rhs, 2) == 0)

end subroutine set_ele_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_branch (ok)

implicit none

type(branch_struct), target :: f_branch, f2_branch
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_branch (c_branch, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_branch
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_branch_test_pattern (f2_branch, 1)

call test_c_branch(c_loc(f2_branch), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_branch_test_pattern (f_branch, 4)
if (f_branch == f2_branch) then
  print *, 'branch: C side convert C->F: Good'
else
  print *, 'branch: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_branch

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_branch (c_branch, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_branch
type(branch_struct), target :: f_branch, f2_branch
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call branch_to_f (c_branch, c_loc(f_branch))

call set_branch_test_pattern (f2_branch, 2)
if (f_branch == f2_branch) then
  print *, 'branch: F side convert C->F: Good'
else
  print *, 'branch: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_branch_test_pattern (f2_branch, 3)
call branch_to_c (c_loc(f2_branch), c_branch)

end subroutine test2_f_branch

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_branch_test_pattern (F, ix_patt)

implicit none

type(branch_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[character, 0, NOT]
do jd1 = 1, len(F%name)
  F%name(jd1:jd1) = char(ichar("a") + modulo(100+1+offset+jd1, 26))
enddo
!! f_side.test_pat[integer, 0, NOT]
rhs = 2 + offset; F%ix_branch = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 3 + offset; F%ix_root_branch = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 4 + offset; F%ix_from_branch = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 5 + offset; F%ix_from_ele = rhs
!! f_side.test_pat[integer, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%n_ele_track)) deallocate (F%n_ele_track)
else
  if (.not. associated(F%n_ele_track)) allocate (F%n_ele_track)
  rhs = 6 + offset
  F%n_ele_track = rhs
endif
!! f_side.test_pat[integer, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%n_ele_max)) deallocate (F%n_ele_max)
else
  if (.not. associated(F%n_ele_max)) allocate (F%n_ele_max)
  rhs = 8 + offset
  F%n_ele_max = rhs
endif
!! f_side.test_pat[type, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%a)) deallocate (F%a)
else
  if (.not. associated(F%a)) allocate (F%a)
  rhs = 10 + offset
  call set_mode_info_test_pattern (F%a, ix_patt)
endif
!! f_side.test_pat[type, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%b)) deallocate (F%b)
else
  if (.not. associated(F%b)) allocate (F%b)
  rhs = 12 + offset
  call set_mode_info_test_pattern (F%b, ix_patt)
endif
!! f_side.test_pat[type, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%z)) deallocate (F%z)
else
  if (.not. associated(F%z)) allocate (F%z)
  rhs = 14 + offset
  call set_mode_info_test_pattern (F%z, ix_patt)
endif
!! f_side.test_pat[type, 1, PTR]

if (ix_patt < 3) then
  if (associated(F%ele)) deallocate (F%ele)
else
  if (.not. associated(F%ele)) allocate (F%ele(-1:1))
  do jd1 = 1, size(F%ele,1); lb1 = lbound(F%ele,1) - 1
    call set_ele_test_pattern (F%ele(jd1+lb1), ix_patt+jd1)
  enddo
endif
!! f_side.test_pat[type, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%param)) deallocate (F%param)
else
  if (.not. associated(F%param)) allocate (F%param)
  rhs = 18 + offset
  call set_lat_param_test_pattern (F%param, ix_patt)
endif
!! f_side.test_pat[type, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%wall3d)) deallocate (F%wall3d)
else
  if (.not. associated(F%wall3d)) allocate (F%wall3d)
  rhs = 20 + offset
  call set_wall3d_test_pattern (F%wall3d, ix_patt)
endif

end subroutine set_branch_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_lat (ok)

implicit none

type(lat_struct), target :: f_lat, f2_lat
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_lat (c_lat, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_lat
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_lat_test_pattern (f2_lat, 1)

call test_c_lat(c_loc(f2_lat), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_lat_test_pattern (f_lat, 4)
if (f_lat == f2_lat) then
  print *, 'lat: C side convert C->F: Good'
else
  print *, 'lat: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_lat

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_lat (c_lat, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_lat
type(lat_struct), target :: f_lat, f2_lat
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call lat_to_f (c_lat, c_loc(f_lat))

call set_lat_test_pattern (f2_lat, 2)
if (f_lat == f2_lat) then
  print *, 'lat: F side convert C->F: Good'
else
  print *, 'lat: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_lat_test_pattern (f2_lat, 3)
call lat_to_c (c_loc(f2_lat), c_lat)

end subroutine test2_f_lat

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_lat_test_pattern (F, ix_patt)

implicit none

type(lat_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[character, 0, NOT]
do jd1 = 1, len(F%use_name)
  F%use_name(jd1:jd1) = char(ichar("a") + modulo(100+1+offset+jd1, 26))
enddo
!! f_side.test_pat[character, 0, NOT]
do jd1 = 1, len(F%lattice)
  F%lattice(jd1:jd1) = char(ichar("a") + modulo(100+2+offset+jd1, 26))
enddo
!! f_side.test_pat[character, 0, NOT]
do jd1 = 1, len(F%input_file_name)
  F%input_file_name(jd1:jd1) = char(ichar("a") + modulo(100+3+offset+jd1, 26))
enddo
!! f_side.test_pat[character, 0, NOT]
do jd1 = 1, len(F%title)
  F%title(jd1:jd1) = char(ichar("a") + modulo(100+4+offset+jd1, 26))
enddo
!! f_side.test_pat[character, 1, ALLOC]
if (ix_patt < 3) then
  if (allocated(F%attribute_alias)) deallocate (F%attribute_alias)
else
  if (.not. allocated(F%attribute_alias)) allocate (F%attribute_alias(3))
  do jd1 = 1, 3
  do jd = 1, len(F%attribute_alias)
    F%attribute_alias(jd1)(jd:jd) = char(ichar("a") + modulo(100+5+offset+10*jd+jd1, 26))
  enddo; enddo
endif
!! f_side.test_pat[type, 0, NOT]
call set_mode_info_test_pattern (F%a, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_mode_info_test_pattern (F%b, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_mode_info_test_pattern (F%z, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_lat_param_test_pattern (F%param, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_bookkeeping_state_test_pattern (F%lord_state, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_ele_test_pattern (F%ele_init, ix_patt)
!! f_side.test_pat[type, 1, PTR]

if (ix_patt < 3) then
  if (associated(F%ele)) deallocate (F%ele)
else
  if (.not. associated(F%ele)) allocate (F%ele(-1:1))
  do jd1 = 1, size(F%ele,1); lb1 = lbound(F%ele,1) - 1
    call set_ele_test_pattern (F%ele(jd1+lb1), ix_patt+jd1)
  enddo
endif
!! f_side.test_pat[type, 1, ALLOC]

if (ix_patt < 3) then
  if (allocated(F%branch)) deallocate (F%branch)
else
  if (.not. allocated(F%branch)) allocate (F%branch(-1:1))
  do jd1 = 1, size(F%branch,1); lb1 = lbound(F%branch,1) - 1
    call set_branch_test_pattern (F%branch(jd1+lb1), ix_patt+jd1)
  enddo
endif
!! f_side.test_pat[type, 1, ALLOC]

if (ix_patt < 3) then
  if (allocated(F%control)) deallocate (F%control)
else
  if (.not. allocated(F%control)) allocate (F%control(-1:1))
  do jd1 = 1, size(F%control,1); lb1 = lbound(F%control,1) - 1
    call set_control_test_pattern (F%control(jd1+lb1), ix_patt+jd1)
  enddo
endif
!! f_side.test_pat[type, 0, NOT]
call set_coord_test_pattern (F%beam_start, ix_patt)
!! f_side.test_pat[type, 0, NOT]
call set_pre_tracker_test_pattern (F%pre_tracker, ix_patt)
!! f_side.test_pat[integer, 0, NOT]
rhs = 21 + offset; F%version = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 22 + offset; F%n_ele_track = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 23 + offset; F%n_ele_max = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 24 + offset; F%n_control_max = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 25 + offset; F%n_ic_max = rhs
!! f_side.test_pat[integer, 0, NOT]
rhs = 26 + offset; F%input_taylor_order = rhs
!! f_side.test_pat[integer, 1, ALLOC]

if (ix_patt < 3) then
  if (allocated(F%ic)) deallocate (F%ic)
else
  if (.not. allocated(F%ic)) allocate (F%ic(-1:1))
  do jd1 = 1, size(F%ic,1); lb1 = lbound(F%ic,1) - 1
    rhs = 100 + jd1 + 27 + offset
    F%ic(jd1+lb1) = rhs
  enddo
endif
!! f_side.test_pat[logical, 0, NOT]
rhs = 29 + offset; F%absolute_time_tracking = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 30 + offset; F%rf_auto_scale_phase = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 31 + offset; F%rf_auto_scale_amp = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, NOT]
rhs = 32 + offset; F%use_ptc_layout = (modulo(rhs, 2) == 0)

end subroutine set_lat_test_pattern

end module
