module test_mod

use, intrinsic :: iso_c_binding

type :: test_struct 
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
type (test_struct), target :: ts

!

call my_c(cptr)

end subroutine

end module

!------------------------------------------------------

program test

use test_mod

implicit none

type (test_struct), target :: ts
type (test_struct), pointer :: ts2
type(c_ptr)  cptr

!

ts%i = 8
ts%j = 999

cptr = c_loc(ts)

call test_sub(cptr)

call c_f_pointer (cptr, ts2)

print *, '!: ', ts2%i, ts2%j

end program
