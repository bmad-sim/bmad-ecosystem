module test_mod

use, intrinsic :: iso_c_binding

type :: zzz_struct 
  integer i, j
end type

contains

subroutine test_sub (cptr) bind(c)

interface 
  subroutine my_c(cp) bind(c)
    import :: c_ptr
    type(c_ptr), value :: cp
  end subroutine
end interface

type(c_ptr)  cptr
type (zzz_struct), target :: ts

!

call my_c(cptr)

end subroutine

!------------------------------------------------------

subroutine test_f1_zzz (ok)

use test_mod

implicit none

type (zzz_struct), target :: zzz1, zzz2
logical ok

!

zzz1%i = 8
zzz1%j = 999

call test_c_zzz(c_loc(zzz1), c_loc(zzz2))

print *, 'test_f1: ', zzz2%i, zzz2%j
ok = .true.

!------------------------------------------------------

subroutine test_f2_zzz (c_zzz1, c_zzz2) bind(c)

implicit  none

type (c_dummy_struct) c_zzz1, c_zzz2

!

end module

!------------------------------------------------------
!------------------------------------------------------
!------------------------------------------------------

program test_main

use test_mod

implicit none

logical ok, all_ok

!

all_ok = .true.
call test_f1_zzz(ok); if (.not. ok) all_ok = .false.

end program

