!+
! Function tao_one_turn_map_calc_needed (data_type, data_source) result (do_one_turn_map)
! 
! Routine decide if a datum or plot curve needs the PTC one turn map normal form calculation.
!-

function tao_one_turn_map_calc_needed (data_type, data_source) result (do_one_turn_map)

implicit none

character(*) data_type, data_source
logical do_one_turn_map

!

do_one_turn_map = .false.
if (data_source /= 'lat') return

call setit (data_type, 'chrom_ptc.', do_one_turn_map)
call setit (data_type, 'momentum_compaction_ptc.', do_one_turn_map) 
call setit (data_type, 'slip_factor_ptc.', do_one_turn_map) 
call setit (data_type, 'spin_tune_ptc.', do_one_turn_map) 
call setit (data_type, 'normal.', do_one_turn_map) 

!------------------------------------------------
contains

subroutine setit(data_type, who, do_one_turn_map)

character(*) data_type, who
integer n
logical do_one_turn_map

!

n = len(who)
if (len(data_type) < n) return
if (data_type(1:n) == who) do_one_turn_map = .true.

end subroutine setit

end function tao_one_turn_map_calc_needed
