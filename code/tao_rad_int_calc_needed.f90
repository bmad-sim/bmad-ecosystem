!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
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
logical do_rad_int

!

do_rad_int = .false.

if (data_source /= 'lat') return

if (data_type  == 'sigma.pz') do_rad_int = .true. 
if (data_type(1:5)  == 'emit.') do_rad_int = .true. 
if (data_type(1:10) == 'norm_emit.') do_rad_int = .true. 
if (data_type(1:7)  == 'rad_int') do_rad_int = .true.
if (data_type(1:16)  == 'apparent_rad_int') do_rad_int = .true.
if (data_type(1:11) == 'expression:') then
  if (index(data_type, 'sigma.pz') /= 0) do_rad_int = .true.
  if (index(data_type, 'emit.') /= 0) do_rad_int = .true.
  if (index(data_type, 'rad_int') /= 0) do_rad_int = .true.
endif

if (tao_beam_sigma_calc_needed(data_type, data_source)) do_rad_int = .true.

end function tao_rad_int_calc_needed
