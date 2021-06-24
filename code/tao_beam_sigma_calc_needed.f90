!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_beam_sigma_calc_needed (data_type, data_source) result (do_beam_sigma)
! 
! Routine decide if a datum or plot curve needs a beam sigma calculation.
!-

function tao_beam_sigma_calc_needed (data_type, data_source) result (do_beam_sigma)


implicit none

character(*) data_type, data_source
logical do_beam_sigma

!

do_beam_sigma = .false.

if (data_source /= 'lat') return
if (len(data_type) < 6) return
if (data_type(1:6)  == 'sigma.') do_beam_sigma = .true. 

end function tao_beam_sigma_calc_needed
