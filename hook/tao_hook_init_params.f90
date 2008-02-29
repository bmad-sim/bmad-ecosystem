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
! Subroutine tao_hook_init_data (is_set) 
!-

subroutine tao_hook_init_data (is_set) 

use tao_mod

implicit none

logical is_set

!

is_set = .false.

end subroutine




!+
! Subroutine tao_hook_init_var(is_set) 
!-

subroutine tao_hook_init_var(is_set) 

use tao_mod

implicit none

logical is_set

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

