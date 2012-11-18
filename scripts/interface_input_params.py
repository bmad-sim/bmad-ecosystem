# List of files containing definitions of the structures to setup interfaces for.

struct_def_files = ['../bmad/modules/bmad_struct.f90']

# List of use statements needed in various Fortran modules

use_statements = ['use bmad_struct']

# List of structures to setup interfaces for.
# List must be in ordered such that if struct A is a component of struct B,
# then A must be before B in the list.

struct_list = [
    'coord_struct',
    'coord_array_struct',
    'bpm_phase_coupling_struct',
    'wig_term_struct',
    'wig_struct',
    'rf_wake_sr_table_struct'
]

to_do_struct_list = [
    'rf_wake_sr_mode_struct',
    'rf_wake_lr_struct',
    'rf_wake_struct',
    'em_field_map_term_struct',
    'em_field_map_struct',
    'em_field_grid_pt_struct',
    'em_field_grid_struct',
    'em_field_mode_struct',
    'em_fields_struct',
    'floor_position_struct',
    'space_charge_struct',
    'xy_disp_struct',
    'mode3_struct',
    'bookkeeper_status_struct',
    'rad_int_ele_cache_struct',
    'ele_struct',
    'control_struct',
    'lat_param_struct',
    'mode_info_struct',
    'branch_struct',
    'pre_tracker_struct',
    'lat_struct',
    'anormal_mode_struct',
    'linac_normal_mode_struct',
    'normal_modes_struct',
    'em_field_struct',
    'track_map_struct',
    'track_struct',
    'synch_rad_common_struct',
    'bmad_common_struct',
    'rad_int1_struct',
    'rad_int_all_ele_struct'
]

# List of sub-structures to ignore.
# That is, do not translate these sub-structure components.

component_ignore_list = set(['fibre', 'genfield'])

# Directory where the output is put

output_dir = 'interface_test'
test_dir = 'interface_test'

# Function to customize the interface code.

def customize(struct_definitions):
  pass
