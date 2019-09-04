!+
! Subroutine ptc_ran_seed_put (iseed)
!
! Routine to set the PTC random number seed used for the radiation excitation part of the map.
!
! Input:
!   iseed -- Integer: 0 -> Use system clock.
!-

subroutine ptc_ran_seed_put (iseed)

use sim_utils
use duan_zhe_map, only: zhe_iseed
use gauss_dis, only: gaussian_seed

implicit none

integer is, iseed, v(10)

! Note: There are actually two independent seeds used by PTC

if (iseed == 0) then
  call date_and_time (values = v)
  is = v(1) + v(2) + 11*v(3) + 111*v(5) + 1111*v(6) + 11111*v(7) + 111111*v(8)
else
  is = iseed
endif

zhe_iseed = is
call gaussian_seed(is)  ! Actually used also for PTC non-gaussian random number generator

end subroutine ptc_ran_seed_put

