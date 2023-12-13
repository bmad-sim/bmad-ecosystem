!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_chrom_calc_needed (data_type, data_source) result (do_chrom)
! 
! Routine decide if a datum or plot curve needs the chromaticity calculation.
!-

function tao_chrom_calc_needed (data_type, data_source) result (do_chrom)

implicit none

character(*) data_type, data_source
logical do_chrom

!

do_chrom = .false.

if (data_source /= 'lat') return
if (len(data_type) < 6) return
if (data_type(1:6)  == 'chrom.') do_chrom = .true. 

end function tao_chrom_calc_needed
