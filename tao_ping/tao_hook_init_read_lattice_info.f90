subroutine tao_hook_init_read_lattice_info (init_file)

use tao_ping_utils

implicit none

character(*) init_file
character(40) :: r_name = 'tao_hook_init_read_lattice_info'

! 

tao_com%init_lat_file = ping_s%param%lattice_file
tao_com%n_universes = 1
tao_com%combine_consecutive_elements_of_like_name = .false.
tao_com%unique_name_suffix = ''
tao_com%aperture_limit_on = ''

end subroutine
