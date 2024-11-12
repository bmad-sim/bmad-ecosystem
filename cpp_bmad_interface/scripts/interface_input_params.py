# List of files containing definitions of the structures to setup interfaces for.

struct_def_files = [
    '../sim_utils/math/spline_mod.f90',
    '../bmad/modules/bmad_struct.f90', 
    '../bmad/modules/taylor_mod.f90',
    '../bmad/space_charge/csr_and_space_charge_mod.f90',
    '../bmad/modules/complex_taylor_mod.f90',
]

# List of use statements needed in various Fortran modules.

conversion_use_statements = ['use bmad_struct']
equality_use_statements = ['use bmad_struct']
test_use_statements = []

# List of structures to setup interfaces for.
# List must be in ordered such that if struct A is a component of struct B,
# then A must be before B in the list.

struct_list = [
    'spline_struct',
    'spin_polar_struct',
    'ac_kicker_time_struct',
    'ac_kicker_freq_struct',
    'ac_kicker_struct',
    'interval1_coef_struct',
    'photon_reflect_table_struct',
    'photon_reflect_surface_struct',
    'coord_struct',
    'coord_array_struct',
    'bpm_phase_coupling_struct',
    'expression_atom_struct',
    'wake_sr_z_long_struct',
    'wake_sr_mode_struct',
    'wake_sr_struct',
    'wake_lr_mode_struct',
    'wake_lr_struct',
    'lat_ele_loc_struct',
    'wake_struct',
    'taylor_term_struct',
    'taylor_struct',
    'em_taylor_term_struct',
    'em_taylor_struct',
    'cartesian_map_term1_struct',
    'cartesian_map_term_struct',
    'cartesian_map_struct',
    'cylindrical_map_term1_struct',
    'cylindrical_map_term_struct',
    'cylindrical_map_struct',
    'grid_field_pt1_struct',
    'grid_field_pt_struct',
    'grid_field_struct',
    'floor_position_struct',
    'high_energy_space_charge_struct',
    'xy_disp_struct',
    'twiss_struct',
    'mode3_struct',
    'bookkeeping_state_struct',
    'rad_map_struct',
    'rad_map_ele_struct',
    'gen_grad1_struct',
    'gen_grad_map_struct',
    'surface_segmented_pt_struct',
    'surface_segmented_struct',
    'surface_h_misalign_pt_struct',
    'surface_h_misalign_struct',
    'surface_displacement_pt_struct',
    'surface_displacement_struct',
    'target_point_struct',
    'surface_curvature_struct',
    'photon_target_struct',
    'photon_material_struct',
    'pixel_pt_struct',
    'pixel_detec_struct',
    'photon_element_struct',
    'wall3d_vertex_struct',
    'wall3d_section_struct',
    'wall3d_struct',
    'ramper_lord_struct',
    'control_struct',
    'control_var1_struct',
    'control_ramp1_struct',
    'controller_struct',
    'ellipse_beam_init_struct',
    'kv_beam_init_struct',
    'grid_beam_init_struct',
    'beam_init_struct',
    'lat_param_struct',
    'mode_info_struct',
    'pre_tracker_struct',
    'anormal_mode_struct',
    'linac_normal_mode_struct',
    'normal_modes_struct',
    'em_field_struct',
    'strong_beam_struct',
    'track_point_struct',
    'track_struct',
    'space_charge_common_struct',
    'bmad_common_struct',
    'rad_int1_struct',
    'rad_int_branch_struct',
    'rad_int_all_ele_struct',
    'ele_struct',
    'complex_taylor_term_struct',
    'complex_taylor_struct',
    'branch_struct',
    'lat_struct',
    'bunch_struct',
    'bunch_params_struct',
    'beam_struct',
    'aperture_point_struct',
    'aperture_param_struct',
    'aperture_scan_struct',
]

# List of structure components to not translate.
# Can specify these using the syntax:
#   <component_struct_name>        or
#   <struct>%<component_name>

component_no_translate_list = set([
  'fibre', 
  'ptc_branch1_info_struct',
  'layout',
  'exact_bend_multipole_struct',
  'branch_struct%ptc',
  'ele_struct%lord',
  'ele_struct%branch',
  'ele_struct%converter',
  'ele_struct%multipole_cache',
  'ele_struct%foil',
  'lat_struct%nametable',
  'branch_struct%lat',
  'normal_form_struct',
  'grid_field_struct%bi_coef',
  'grid_field_struct%tri_coef',
])

# List of structure components links:
# Structure components that are just links to other structures are handled differently.
#   1) No call to delete in the C++ destructor.
#   2) Ignore in Fortran and C++ equality tests (could go around in circles).
#   3) Do not create a test pattern in interface test code.

interface_ignore_list = set([
  'ele_struct%branch',
  'branch_struct%lat',
  'pixel_grid_struct',
])

# List of structure components that are structures and are defined externally.
# There are no such structures for the cpp_bmad_interface library but there
# are for the cpp_tao_interface library.

structs_defined_externally = set([])

# Translations on C++ side to avoid clash with reserved words

c_side_name_translation = {
    'wake_sr_struct%long' : 'long_wake',
    'wake_sr_struct%trans' : 'trans_wake',
}

# Include header files for main header file

include_header_files = [
  '#include "bmad_enums.h"',
  '#include "bmad_std_typedef.h"',
]

# Directory where the output is put

equality_mod_dir  = '../bmad/modules'
equality_mod_file = 'equality_mod'
test_dir          = 'interface_test'
code_dir          = 'code'

# Lower bounds for allocatable and pointer arrays on the fortran side

def f_side_lbound (id_name):
  if id_name == 'branch%ele':
    return '0'
  else:
    return '1'

# custom C++ side init

c_custom_constructors = {
    'ele%key' : 'key(key_)',
    'ele%value' : 'value(double(0), Bmad::NUM_ELE_ATTRIB+1)',
    'ele%old_value' : 'old_value(double(0), Bmad::NUM_ELE_ATTRIB+1)', 
    'ele%ix_ele' : 'ix_ele(-1)',
    'ele%slave_status' : 'slave_status(Bmad::FREE)', 
    'ele%ix2_slave' : 'ix2_slave(-1)',
    'ele%lord_status' : 'lord_status(Bmad::NOT_A_LORD)', 
    'ele%ic2_lord' : 'ic2_lord(-1)',
    'ele%mat6_calc_method' : 'mat6_calc_method(Bmad::BMAD_STANDARD)', 
    'ele%tracking_method' : 'tracking_method(Bmad::BMAD_STANDARD)',
    'ele%spin_tracking_method' : 'spin_tracking_method(Bmad::BMAD_STANDARD)',
    'ele%field_calc' : 'field_calc(Bmad::BMAD_STANDARD)',
    'ele%ptc_integration_type' : 'ptc_integration_type(Bmad::MATRIX_KICK)',
    'ele%aperture_at' : 'aperture_at(Bmad::DOWNSTREAM_END)', 
    'ele%aperture_type' : 'aperture_type(Bmad::RECTANGULAR)',
    'ele%multipoles_on' : 'multipoles_on(true)', 
    'ele%scale_multipoles' : 'scale_multipoles(true)', 
    'ele%map_with_offsets' : 'map_with_offsets(true)',
    'ele%is_on' : 'is_on(true)', 
    'ele%csr_calc_on' : 'csr_calc_on(true)',
    'ele%orientation' : 'orientation(1)',
    'floor_position%w' : 'w(Real_ARRAY(0.0, 3), 3)',
    'aperture_param%max_angle' : 'max_angle(Bmad::PI)',
    'bmad_common%space_charge_mesh_size' : 'space_charge_mesh_size(32, 3)',
    'rad_map%xfer_damp_mat' : 'xfer_damp_mat(Real_ARRAY(0.0, 6), 6)',
}

#-----------------------------------------------
# Function to customize the interface code.

def customize(struct_definitions):

  for struct in struct_definitions:

    if struct.f_name == 'ele_struct': 
      struct.c_constructor_arg_list = 'const int key_ = 0'
      struct.c_extra_methods = '''
  void class_init (const int key_) {
    key = key_;

    if (key == Bmad::LCAVITY) {
      value[Bmad::COUPLER_AT] = Bmad::DOWNSTREAM_END;
      value[Bmad::FIELD_AUTOSCALE] = 1;
      value[Bmad::N_CELL] = 1;
    }

    if (key == Bmad::RFCAVITY) {
      value[Bmad::COUPLER_AT] = Bmad::DOWNSTREAM_END;
      value[Bmad::FIELD_AUTOSCALE] = 1;
      value[Bmad::N_CELL] = 1;
    }

    if (key == Bmad::RBEND || key == Bmad::SBEND) {
      value[Bmad::FRINGE_AT] = Bmad::BOTH_ENDS;
      value[Bmad::FRINGE_TYPE] = Bmad::BASIC_BEND;
      value[Bmad::PTC_FIELD_GEOMETRY] = Bmad::SECTOR;
    }
  }
'''

      struct.c_constructor_body = '''\
    {
      class_init(key);
    }
'''
    for arg in struct.arg:

      id_name = struct.short_name + '%' + arg.f_name 

      if id_name in c_custom_constructors: 
        arg.c_side.constructor = c_custom_constructors[id_name]

      # ele%value and ele%old_value must be handled specially since the array sizes are different
      # between Fortran and C++.

      if id_name == 'ele%value' or id_name == 'ele%old_value':
        arg.c_side.to_c2_set = '''\
  C.NAME[0] = 0;
  for (unsigned int i = 1; i < Bmad::NUM_ELE_ATTRIB+1; i++) C.NAME[i] = z_NAME[i-1];
'''.replace('NAME', arg.f_name)

        arg.c_side.test_pat = '''\
  C.NAME[0] = 0;
  for (unsigned int i = 1; i < Bmad::NUM_ELE_ATTRIB+1; i++)
    {int rhs = 100 + i + XXX + offset; C.NAME[i] = rhs;}
'''.replace('NAME', arg.f_name)

        arg.f_side.to_f2_trans = 'F%NAME = z_NAME(2:num_ele_attrib$+1)'.replace('NAME', arg.f_name)
