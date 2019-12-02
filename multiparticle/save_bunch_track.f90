!+
! Subroutine save_bunch_track (bunch, ele, s_travel)
!
! Routine is called by track1_bunch_csr to allow recoding of bunch statistics.
! This file contains a dummy routine that does nothing. That is, to be useful,
! a program must be linked with a custom save_bunch_track routine.
! 
! Input:
!   bunch    -- bunch_struct: A bunch of particles.
!   ele      -- ele_struct: The lattice element the bunch is going through.
!   s_travel -- Real(rp): Longitudinal distance from beginning of element.
!-

subroutine save_bunch_track (bunch, ele, s_travel)

use bmad_struct, only: bunch_struct, ele_struct, rp

implicit none

type (bunch_struct) bunch
type (ele_struct) ele
real(rp) s_travel

!

end subroutine
