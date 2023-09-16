!+
! Function tao_spin_matrices_calc_needed (data_type, data_source) result (do_calc)
! 
! Routine decide if a datum needs spin related matrices calculated.
!-

function tao_spin_matrices_calc_needed (data_type, data_source) result (do_calc)

implicit none

character(*) data_type, data_source
logical do_calc

!

do_calc = .false.
if (len(data_type) < 14) return
if (data_type(1:14)  == 'spin_g_matrix.') do_calc = .true. 

end function tao_spin_matrices_calc_needed
