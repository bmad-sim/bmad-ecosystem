# List of files containing definitions of the structures to setup interfaces for.

struct_def_files = ['interface_test/test_struct_defs.f90']

# List of use statements needed in various Fortran modules

use_statements = ['use test_mod']

# List of structures to setup interfaces for.
# List must be in ordered such that if struct A is a component of struct B,
# then A must be before B in the list.

struct_list = ['z_struct', 'ttt_struct']

# List of sub-structures to ignore.
# That is, do not translate these sub-structure components.

component_ignore_list = set(['fibre', 'genfield'])

# Directory where the output is put

output_dir = 'interface_test'
test_dir = 'interface_test'

# Function to customize the interface code.

def customize(struct_definitions):
  pass
