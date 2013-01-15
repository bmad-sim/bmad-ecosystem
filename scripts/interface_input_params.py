# List of files containing definitions of the structures to setup interfaces for.

struct_def_files = ['../bmad/modules/bmad_struct.f90', 
                    '../bmad/modules/twiss_mod.f90', 
                    '../bmad/modules/bmad_taylor_mod.f90'
]

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
    'rf_wake_sr_table_struct',
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
    'twiss_struct',
    'mode3_struct',
    'bookkeeping_state_struct',
    'rad_int_ele_cache_struct',
    'wall3d_vertex_struct',
    'wall3d_section_struct',
    'wall3d_crotch_struct',
    'wall3d_struct',
    'taylor_term_struct',
    'taylor_struct',
    'control_struct',
    'lat_param_struct',
    'mode_info_struct',
    'pre_tracker_struct',
    'anormal_mode_struct',
    'linac_normal_mode_struct',
    'normal_modes_struct',
    'em_field_struct',
    'track_map_struct',
    'track_struct',
    'synch_rad_common_struct',
    'bmad_common_struct',
    'rad_int1_struct',
    'rad_int_all_ele_struct',
    'ele_struct',
    'branch_struct',
    'lat_struct'
]

# List of structure components to not translate.
# Can specify these using the syntax:
#   <component_struct_name>        or
#   <struct>%<component_name>

component_no_translate_list = set([
    'fibre', 
    'genfield',
    'ptc_branch1_info_struct',
    'layout',
    'ele_struct%lord'
])

# List of structure components links:
# Structure components that are just links to other structures are handled differently.
#   1) Not call delete in the C++ destructor.
#   2) Ignore in Fortran and C++ equality tests (could go around in circles).
#   3) Do not create a test pattern in interface test code.

interface_ignore_list = set([
  'ele_struct%branch',
  'branch_struct%lat'
])


# Translations on C++ side to avoid clash with reserved words

c_side_name_translation = {
    'rf_wake_sr_table_struct%long' : 'long_wake',
    'rf_wake_sr_table_struct%trans' : 'trans_wake'
}

# Directory where the output is put

equality_mod_dir  = '../bmad/modules'
test_dir          = 'interface_test'
code_dir          = 'code'

# Lower bounds for allocatable and pointer arrays on the fortran side

def f_side_lbound (id_name):
  if id_name == 'branch%ele':
    return '0'
  else:
    return '1'

# Function to customize the interface code.

def customize(struct_definitions):
  pass
