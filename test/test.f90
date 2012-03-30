program test_main

implicit none

interface
  subroutine c_test (i, r, i1, r1, str) bind(c)
    use iso_c_binding
    implicit none
    integer(c_int) i, i1(*)
    real(c_double) r, r1(*)
    character(c_char), intent(in) :: str(*)
  end subroutine
end interface

integer i, i1(2)
real(8), pointer :: rp, rp1(:)
real(8), target :: r, r1(3)
character(8) str

!

i = 7
i1 = [11, 22]
r = 45
r1 = [100, 101, 102]
rp => r
rp1 => r1

call c_test (i, rp, i1, rp1, str)

end program

