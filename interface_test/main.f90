
program interface_test

use bmad_cpp_test_mod

logical ok, all_ok

!

all_ok = .true.
call test1_f_my(ok); if (.not. ok) all_ok = .false.
call test1_f_ttt(ok); if (.not. ok) all_ok = .false.

if (all_ok) then
  print *, 'Bottom Line: Everything OK!'
else
  print *, 'BOTTOM LINE: PROBLEMS FOUND!'
endif

end program
