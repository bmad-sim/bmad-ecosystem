
module bmad_cpp_test_mod

use bmad_cpp_convert_mod
use equality_mod

contains

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_my (ok)

implicit none

type(my_struct), target :: f_my, f2_my
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_my (c_my, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_my
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_my_test_pattern (f2_my, 1)

call test_c_my(c_loc(f2_my), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_my_test_pattern (f_my, 4)
if (f_my == f2_my) then
  print *, 'my: C side convert C->F: Good'
else
  print *, 'my: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_my

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_my (c_my, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_my
type(my_struct), target :: f_my, f2_my
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call my_to_f (c_my, c_loc(f_my))

call set_my_test_pattern (f2_my, 2)
if (f_my == f2_my) then
  print *, 'my: F side convert C->F: Good'
else
  print *, 'my: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_my_test_pattern (f2_my, 3)
call my_to_c (c_loc(f2_my), c_my)

end subroutine test2_f_my

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_my_test_pattern (F, ix_patt)

implicit none

type(my_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[integer, 0, NOT]
rhs = 1 + offset; F%a = rhs

end subroutine set_my_test_pattern

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test1_f_ttt (ok)

implicit none

type(ttt_struct), target :: f_ttt, f2_ttt
logical(c_bool) c_ok
logical ok

interface
  subroutine test_c_ttt (c_ttt, c_ok) bind(c)
    import c_ptr, c_bool
    type(c_ptr), value :: c_ttt
    logical(c_bool) c_ok
  end subroutine
end interface

!

ok = .true.
call set_ttt_test_pattern (f2_ttt, 1)

call test_c_ttt(c_loc(f2_ttt), c_ok)
if (.not. f_logic(c_ok)) ok = .false.

call set_ttt_test_pattern (f_ttt, 4)
if (f_ttt == f2_ttt) then
  print *, 'ttt: C side convert C->F: Good'
else
  print *, 'ttt: C SIDE CONVERT C->F: FAILED!'
  ok = .false.
endif

end subroutine test1_f_ttt

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine test2_f_ttt (c_ttt, c_ok) bind(c)

implicit  none

type(c_ptr), value ::  c_ttt
type(ttt_struct), target :: f_ttt, f2_ttt
logical(c_bool) c_ok

!

c_ok = c_logic(.true.)
call ttt_to_f (c_ttt, c_loc(f_ttt))

call set_ttt_test_pattern (f2_ttt, 2)
if (f_ttt == f2_ttt) then
  print *, 'ttt: F side convert C->F: Good'
else
  print *, 'ttt: F SIDE CONVERT C->F: FAILED!'
  c_ok = c_logic(.false.)
endif

call set_ttt_test_pattern (f2_ttt, 3)
call ttt_to_c (c_loc(f2_ttt), c_ttt)

end subroutine test2_f_ttt

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine set_ttt_test_pattern (F, ix_patt)

implicit none

type(ttt_struct) F
integer ix_patt, offset, jd, jd1, jd2, jd3, lb1, lb2, lb3, rhs

!

offset = 100 * ix_patt

!! f_side.test_pat[logical, 0, NOT]
rhs = 1 + offset; F%i0 = (modulo(rhs, 2) == 0)
!! f_side.test_pat[logical, 0, PTR]
if (ix_patt < 3) then
  if (associated(F%ip0)) deallocate (F%ip0)
else
  if (.not. associated(F%ip0)) allocate (F%ip0)
  rhs = 2 + offset
  F%ip0 = (modulo(rhs, 2) == 0)
endif
!! f_side.test_pat[logical, 0, ALLOC]
if (ix_patt < 3) then
  if (allocated(F%ia0)) deallocate (F%ia0)
else
  if (.not. allocated(F%ia0)) allocate (F%ia0)
  rhs = 4 + offset
  F%ia0 = (modulo(rhs, 2) == 0)
endif
!! f_side.test_pat[logical, 1, NOT]
do jd1 = 1, size(F%i1,1); lb1 = lbound(F%i1,1) - 1
  rhs = 100 + jd1 + 6 + offset
  F%i1(jd1+lb1) = (modulo(rhs, 2) == 0)
enddo
!! f_side.test_pat[logical, 1, PTR]

if (ix_patt < 3) then
  if (associated(F%ip1)) deallocate (F%ip1)
else
  if (.not. associated(F%ip1)) allocate (F%ip1(-1:1))
  do jd1 = 1, size(F%ip1,1); lb1 = lbound(F%ip1,1) - 1
    rhs = 100 + jd1 + 7 + offset
    F%ip1(jd1+lb1) = (modulo(rhs, 2) == 0)
  enddo
endif
!! f_side.test_pat[logical, 1, ALLOC]

if (ix_patt < 3) then
  if (allocated(F%ia1)) deallocate (F%ia1)
else
  if (.not. allocated(F%ia1)) allocate (F%ia1(-1:1))
  do jd1 = 1, size(F%ia1,1); lb1 = lbound(F%ia1,1) - 1
    rhs = 100 + jd1 + 9 + offset
    F%ia1(jd1+lb1) = (modulo(rhs, 2) == 0)
  enddo
endif
!! f_side.test_pat[logical, 2, NOT]
do jd1 = 1, size(F%i2,1); lb1 = lbound(F%i2,1) - 1
do jd2 = 1, size(F%i2,2); lb2 = lbound(F%i2,2) - 1
  rhs = 100 + jd1 + 10*jd2 + 11 + offset
  F%i2(jd1+lb1,jd2+lb2) = (modulo(rhs, 2) == 0)
enddo; enddo
!! f_side.test_pat[logical, 2, PTR]

if (ix_patt < 3) then
  if (associated(F%ip2)) deallocate (F%ip2)
else
  if (.not. associated(F%ip2)) allocate (F%ip2(-1:1, 2))
  do jd1 = 1, size(F%ip2,1); lb1 = lbound(F%ip2,1) - 1
  do jd2 = 1, size(F%ip2,2); lb2 = lbound(F%ip2,2) - 1
    rhs = 100 + jd1 + 10*jd2 + 12 + offset
    F%ip2(jd1+lb1,jd2+lb2) = (modulo(rhs, 2) == 0)
  enddo; enddo
endif
!! f_side.test_pat[logical, 2, ALLOC]

if (ix_patt < 3) then
  if (allocated(F%ia2)) deallocate (F%ia2)
else
  if (.not. allocated(F%ia2)) allocate (F%ia2(-1:1, 2))
  do jd1 = 1, size(F%ia2,1); lb1 = lbound(F%ia2,1) - 1
  do jd2 = 1, size(F%ia2,2); lb2 = lbound(F%ia2,2) - 1
    rhs = 100 + jd1 + 10*jd2 + 15 + offset
    F%ia2(jd1+lb1,jd2+lb2) = (modulo(rhs, 2) == 0)
  enddo; enddo
endif
!! f_side.test_pat[logical, 3, NOT]
do jd1 = 1, size(F%i3,1); lb1 = lbound(F%i3,1) - 1
do jd2 = 1, size(F%i3,2); lb2 = lbound(F%i3,2) - 1
do jd3 = 1, size(F%i3,3); lb3 = lbound(F%i3,3) - 1
  rhs = 100 + jd1 + 10*jd2 + 100*jd3 + 18 + offset
  F%i3(jd1+lb1,jd2+lb2,jd3+lb3) = (modulo(rhs, 2) == 0)
enddo; enddo; enddo
!! f_side.test_pat[logical, 3, PTR]
if (ix_patt < 3) then
  if (associated(F%ip3)) deallocate (F%ip3)
else
  if (.not. associated(F%ip3)) allocate (F%ip3(-1:1, 2, 1))
  do jd1 = 1, size(F%ip3,1); lb1 = lbound(F%ip3,1) - 1
  do jd2 = 1, size(F%ip3,2); lb2 = lbound(F%ip3,2) - 1
  do jd3 = 1, size(F%ip3,3); lb3 = lbound(F%ip3,3) - 1
    rhs = 100 + jd1 + 10*jd2 + 100*jd3 + 19 + offset
    F%ip3(jd1+lb1,jd2+lb2,jd3+lb3) = (modulo(rhs, 2) == 0)
  enddo; enddo; enddo
endif
!! f_side.test_pat[logical, 3, ALLOC]
if (ix_patt < 3) then
  if (allocated(F%ia3)) deallocate (F%ia3)
else
  if (.not. allocated(F%ia3)) allocate (F%ia3(-1:1, 2, 1))
  do jd1 = 1, size(F%ia3,1); lb1 = lbound(F%ia3,1) - 1
  do jd2 = 1, size(F%ia3,2); lb2 = lbound(F%ia3,2) - 1
  do jd3 = 1, size(F%ia3,3); lb3 = lbound(F%ia3,3) - 1
    rhs = 100 + jd1 + 10*jd2 + 100*jd3 + 23 + offset
    F%ia3(jd1+lb1,jd2+lb2,jd3+lb3) = (modulo(rhs, 2) == 0)
  enddo; enddo; enddo
endif

end subroutine set_ttt_test_pattern

end module
