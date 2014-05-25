subroutine init_wall (wall)

use synrad_struct
use synrad_interface, except => init_wall

implicit none

type (wall_struct) wall

integer i

!

wall%seg(:)%power%power_tot = 0
wall%seg(:)%power%power_per_len = 0
wall%seg(:)%power%power_per_area = 0
wall%seg(:)%power%n_source = 0
wall%seg(:)%power%main_source%ix_ele = 0
wall%seg(:)%power%main_source%power_per_len = 0
wall%seg(:)%power%main_source%s = 0

end subroutine
