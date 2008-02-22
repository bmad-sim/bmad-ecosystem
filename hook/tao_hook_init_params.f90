!+
! Subroutine tao_hook_init_params (global, is_set)
!
! This routine is part of a collection of hook routines 
! used to bypass the use of an initialization file.
!-

subroutine tao_hook_init_params (global, is_set)

use tao_mod
use tao_input_struct, only: tao_global_struct

implicit none

type (tao_global_struct) global
logical is_set

!

is_set = .false.

end subroutine



!+
! Subroutine tao_hook_init_connected_uni (ix_universe, connect, is_set, ended)
!-

subroutine tao_hook_init_connected_uni (ix_universe, connect, is_set, ended)

use tao_mod
use tao_input_struct, only: tao_connected_uni_input

implicit none

integer ix_universe
type (tao_connected_uni_input) connect
logical is_set, ended

!

is_set = .false.

end subroutine




!+
! Subroutine tao_hook_count_d2_data (n_d2_data, is_set) 
!-

subroutine tao_hook_count_d2_data (n_d2_data, is_set) 

use tao_mod

implicit none

integer n_d2_data(lbound(s%u, 1):)
logical is_set

!

is_set = .false.

end subroutine




!+
! Subroutine tao_hook_count_var(n_var, is_set) 
!-

subroutine tao_hook_count_var(n_var, is_set) 

use tao_mod

implicit none

integer n_var
logical is_set

!

is_set = .false.

end subroutine

!+
!
!-

subroutine tao_hook_read_d1_data (d1_data, data, ix_d1_data, ix_min_data, ix_max_data, &
               default_weight, default_data_type, default_data_source, &
               use_same_lat_eles_as, search_for_lat_eles, is_set)

use tao_mod
use tao_input_struct, only: tao_d1_data_input, tao_data_input, n_data_minn, &
                            n_data_maxx

implicit none

type (tao_d1_data_input) d1_data
type (tao_data_input) data(n_data_minn:n_data_maxx)
integer ix_d1_data, ix_min_data, ix_max_data
real(rp) default_weight        ! default merit function weight
character(*) default_data_type, default_data_source
character(*) use_same_lat_eles_as, search_for_lat_eles
logical is_set

!

is_set = .false.

end subroutine


!+
!
!-

subroutine tao_hook_read_d2_data (d2_data, n_d1_data, &
                             default_merit_type, universe, is_set, ended)

use tao_mod
use tao_input_struct, only: tao_d2_data_input

implicit none

type (tao_d2_data_input) d2_data
integer n_d1_data
character(*) default_merit_type, universe
logical is_set, ended

!

is_set = .false.

end subroutine


!+
!
!-

subroutine tao_hook_read_var (v1_var, var, default_weight, default_step, default_key_delta, &
                ix_min_var, ix_max_var, default_universe, default_attribute, &
                default_low_lim, default_high_lim, default_merit_type, &
                use_same_lat_eles_as, search_for_lat_eles, default_key_bound, &
                is_set, ended)

use tao_mod
use tao_input_struct, only: n_var_minn, n_var_maxx, tao_v1_var_input, tao_var_input

implicit none

type (tao_v1_var_input) v1_var
type (tao_var_input) var(n_var_minn:n_var_maxx)
real(rp) default_weight, default_step, default_key_delta, &
                default_low_lim, default_high_lim
integer ix_min_var, ix_max_var
character(*) default_universe, default_attribute, default_merit_type
character(*) use_same_lat_eles_as, search_for_lat_eles
character(*) default_key_bound
logical is_set, ended

!

is_set = .false.

end subroutine
 

!+
!
!-

subroutine tao_hook_init_beam (ix_universe, beam0_file, ix_track_start, ix_track_end, &
            beam_all_file, beam_init, save_beam_at, is_set, ended)

use tao_mod
use tao_input_struct, only: beam_init_struct

implicit none

type (beam_init_struct) beam_init
integer ix_universe, ix_track_start, ix_track_end
character(*) beam0_file, beam_all_file
character(*) save_beam_at(:)
logical is_set, ended

!

is_set = .false.

end subroutine

