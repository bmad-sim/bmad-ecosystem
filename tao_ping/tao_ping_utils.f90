module tao_ping_utils

use tao_ping_struct

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

subroutine ping_read_parameters ()

implicit none

type (ping_param_struct) param

namelist / ping_params / param

!

open (1, file = 'tao_ping.init', status = 'old')
read (1, nml = ping_params)
close (1)

ping_s%param = param


end subroutine

end module
