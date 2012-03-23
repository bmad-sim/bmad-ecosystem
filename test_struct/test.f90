module cpp_convert_mod

use, intrinsic :: iso_c_binding

type :: zzz_struct 
  integer i, j
end type

interface 
  subroutine zzz_to_f (c_zzz, f_zzz) bind(c)
    use, intrinsic :: iso_c_binding
    type (c_ptr), value :: c_zzz, f_zzz
  end subroutine
end interface

contains

!------------------------------------------------------

subroutine zzz_to_c (f_zzz, c_zzz) bind(c)

implicit none

type (c_ptr), value :: f_zzz
type (c_ptr), value :: c_zzz
type (zzz_struct), pointer :: fp_zzz

integer(c_int), target :: c_int_vec(2)

interface 
  subroutine zzz_to_c2 (c_zzz, c_int_vec) bind(c)
    use, intrinsic :: iso_c_binding
    type (c_ptr), value :: c_zzz
    type (c_ptr), value :: c_int_vec
  end subroutine
end interface

call c_f_pointer (f_zzz, fp_zzz)

c_int_vec = [fp_zzz%i, fp_zzz%j]
call zzz_to_c2 (c_zzz, c_loc(c_int_vec))

end subroutine zzz_to_c

!------------------------------------------------------

subroutine zzz_to_f2 (f_zzz, c_int_vec) bind(c)

implicit none

type (c_ptr), value :: f_zzz
integer(c_int) :: c_int_vec(2)
type (zzz_struct), pointer :: fp_zzz

call c_f_pointer (f_zzz, fp_zzz)

fp_zzz%i = c_int_vec(1)
fp_zzz%j = c_int_vec(2)

end subroutine zzz_to_f2

end module

!------------------------------------------------------
!------------------------------------------------------
!------------------------------------------------------

module test_mod

use cpp_convert_mod

contains

!------------------------------------------------------

subroutine test1_f_zzz (ok)

implicit none

type (zzz_struct), target :: f1_zzz
integer c_ok
logical ok

interface 
  subroutine test_c_zzz (c_zzz, c_ok) bind(c)
    use, intrinsic :: iso_c_binding
    type (c_ptr), value :: c_zzz
    integer(c_int) c_ok
  end subroutine
end interface

!

f1_zzz%i = 8
f1_zzz%j = 999

call test_c_zzz(c_loc(f1_zzz), c_ok)

print *, 'C_side_convert C->F: ', f1_zzz%i, f1_zzz%j
ok = .true.

end subroutine test1_f_zzz

!------------------------------------------------------

subroutine test2_f_zzz (c_zzz) bind(c)

implicit  none

type (c_ptr), value ::  c_zzz
type (zzz_struct), target :: f2_zzz

!

call zzz_to_f (c_zzz, c_loc(f2_zzz))

print *, 'F_side_convert C->F:', f2_zzz%i, f2_zzz%j

f2_zzz%i = 4
f2_zzz%j = 5

call zzz_to_c (c_loc(f2_zzz), c_zzz)

end subroutine test2_f_zzz

end module test_mod

!------------------------------------------------------
!------------------------------------------------------
!------------------------------------------------------

program test_main

use test_mod

implicit none

logical ok, all_ok

!

all_ok = .true.
call test1_f_zzz(ok); if (.not. ok) all_ok = .false.

end program

