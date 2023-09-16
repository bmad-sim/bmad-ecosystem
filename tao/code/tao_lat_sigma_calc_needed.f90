!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_lat_sigma_calc_needed (data_type, data_source) result (do_lat_sigma)
! 
! Routine decide if a datum or plot curve needs a lattice based beam sigma matrix calculation.
!-

function tao_lat_sigma_calc_needed (data_type, data_source) result (do_lat_sigma)


implicit none

character(*) data_type, data_source
logical do_lat_sigma

!

do_lat_sigma = .false.

if (data_source /= 'lat') return
if (len(data_type) < 6) return
if (data_type(1:6)  == 'sigma.') do_lat_sigma = .true. 

end function tao_lat_sigma_calc_needed
