!+
! Function tao_srdt_calc_needed (data_type, data_source) result (do_rad_int)
! 
! Routine decide if a datum or plot curve needs the summation RDTs
! to be evaluated. 0 = not needed, 1 = first order, 2 = second order
!-

function tao_srdt_calc_needed (data_type, data_source) result (do_srdt)

use srdt_mod

implicit none

character(*) data_type, data_source
integer do_srdt

!

do_srdt = 0
if (data_source /= 'lat') return

if (data_type(1:5)  == 'srdt.') then
  if(any(data_type(6:11) == srdt_first)) then
    do_srdt = 1
  elseif(any(data_type(6:11) == srdt_second)) then
    do_srdt = 2
  endif
endif

end function tao_srdt_calc_needed
