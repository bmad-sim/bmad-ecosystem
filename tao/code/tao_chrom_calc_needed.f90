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

! Note: The data_type may be an expression.

do_chrom = .false.

if (data_source /= 'lat') return
if (len(data_type) < 6) return
if (index(data_type, 'chrom.') /= 0) do_chrom = .true. 
if (index(data_type, '_dpz') /= 0 .and. index(data_type, 'spin') == 0) do_chrom = .true. 

end function tao_chrom_calc_needed
