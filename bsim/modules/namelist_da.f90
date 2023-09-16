module namelist_da

! provides:
! tracking_method, n_adts, n_turn, n_angle, track_dims, dE, init_len,
! adts_x_min, adts_x_max, adts_y_min, adts_y_max

use precision_def

implicit none

integer, parameter :: max_dE = 11

integer tracking_method
integer n_adts
integer n_turn
integer n_angle
integer track_dims
real(rp) dE(max_dE)
real(rp) init_len
real(rp) adts_x_min, adts_x_max
real(rp) adts_y_min, adts_y_max

namelist / da /     tracking_method, &
                    n_adts, &
                    n_turn, &          !Number of turns particle must survive
                    n_angle, &
                    track_dims, &            !either 4 or 6
                    dE, &
                    adts_x_min, &
                    adts_x_max, &
                    adts_y_min, &
                    adts_y_max, &
                    init_len

end module
