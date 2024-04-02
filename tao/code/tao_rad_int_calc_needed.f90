!+
! Function tao_rad_int_calc_needed (data_type, data_source) result (do_rad_int)
! 
! Routine decide if a datum or plot curve needs the radiation integrals 
! to be evaluated.
!-

function tao_rad_int_calc_needed (data_type, data_source) result (do_rad_int)

use tao_interface, dummy => tao_rad_int_calc_needed

implicit none

character(*) data_type, data_source
integer n
logical do_rad_int

!

do_rad_int = .false.

if (data_source /= 'lat') return
n = len(data_type)

if (data_type == 'sigma.pz') do_rad_int = .true. 
if (data_type(1:min(5,n))  == 'emit.') do_rad_int = .true. 
if (data_type(1:min(10,n)) == 'norm_emit.') do_rad_int = .true. 
if (data_type(1:min(7,n))  == 'rad_int') do_rad_int = .true.
if (data_type(1:min(16,n)) == 'apparent_rad_int') do_rad_int = .true.
if (data_type(1:min(11,n)) == 'expression:') then
  if (index(data_type, 'sigma.pz') /= 0) do_rad_int = .true.
  if (index(data_type, 'emit.') /= 0) do_rad_int = .true.
  if (index(data_type, 'rad_int') /= 0) do_rad_int = .true.
endif

if (tao_lat_sigma_calc_needed(data_type, data_source)) do_rad_int = .true.

end function tao_rad_int_calc_needed
