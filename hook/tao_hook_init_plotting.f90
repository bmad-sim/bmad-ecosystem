!+
! Subroutine tao_hook_init_plotting (place) 
!
! Hook routine to initialize plotting.
!
! Input:
!   place(:)    -- tao_place_input: Placement of plots
!
! Output:
!   place(:)    -- tao_place_input: Possibly modified placement of plots.
!-

subroutine tao_hook_init_plotting (place)

use tao_input_struct
use tao_mod, dummy => tao_hook_init_plotting

implicit none

type (tao_place_input) place(:)

!

end subroutine
