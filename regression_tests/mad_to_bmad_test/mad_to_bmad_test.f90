! Program to test mad_to_bmad translation

program mad_to_bmad_test

use bmad

implicit none

logical err


!

call system_command ('rm ', err)
call system_command ('python3', err)

end program
