program lux

use lux_module

implicit none

type (lux_common_struct), target :: lux_com
type (lux_param_struct) lux_param
type (lux_output_data_struct) lux_data
type (lat_struct), pointer :: lat

integer nx, ny

!------------------------------------------

call lux_init (lux_param, lux_com)
call lux_init_data (lux_param, lux_com, lux_data)
call lux_track_photons (lux_param, lux_com, lux_data)
call lux_write_data (lux_param, lux_com, lux_data)

end program
