debug = True

# List of files containing definitions of the structures to setup interfaces for.

struct_def_files = ['scripts/test_struct_defs.f90']

# List of use statements needed in various Fortran modules

conversion_use_statements = ['use bmad_struct']
equality_use_statements = ['use bmad_struct']
test_use_statements = []

# List of structures to setup interfaces for.
# List must be in ordered such that if struct A is a component of struct B,
# then A must be before B in the list.

struct_list = ['wake_lr_struct', 'wake_struct']

# List of sub-structures to ignore.
# That is, do not translate these sub-structure components.

component_no_translate_list = set(['fibre', 'genfield'])
interface_ignore_list = set(['fibre', 'genfield'])
structs_defined_externally = set([])

# Directory where the output is put

equality_mod_dir  = '.'
equality_mod_file = 'test_equality_mod'
test_dir          = 'interface_test'
code_dir          = 'code'

# Lower bounds for allocatable and pointer arrays on the fortran side

def f_side_lbound (id_name):
  if id_name == 'branch%ele':
    return '0'
  else:
    return '1'

# Function to customize the interface code.

c_side_name_translation = {
    'rf_wake_sr_table_struct%long' : 'long_wake',
    'rf_wake_sr_table_struct%trans' : 'trans_wake'
}

# Include header files for main header file

include_header_files = [
  '#include "bmad_enums.h"',
  '#include "bmad_std_typedef.h"',
]

#
def customize(struct_definitions):
  pass
