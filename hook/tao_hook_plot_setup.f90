!+
! Subroutine tao_hook_plot_setup (place)
!
! Subroutine for custimizing setting up the plot data.
!
! Input:
!   place(:)    -- tao_place_input: Placement of plots
!
! Output:
!   place(:)    -- tao_place_input: Possibly modified placement of plots.
!-

subroutine tao_hook_plot_setup (place)

use tao_input_struct
use tao_mod, dummy => tao_hook_plot_setup

implicit none

type (tao_place_input) place(:)

! This is just a dummy subroutine used as the default.

end subroutine
