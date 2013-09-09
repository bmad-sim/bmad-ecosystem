integer function ttychr()

use sim_utils 
use input_mod

implicit none

character*1 ic,ignore_this(2)

! wrapper to get 1 char, non-blocking  

ignore_this(1)=' '
call get_a_char(ic,.false.,ignore_this)  !dont wait, ignore blank
ttychr=ichar(ic)

end function
