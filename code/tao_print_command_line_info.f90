!+
! Subroutine tao_print_command_line_info
!
! Routine to print a list of the command line options.
!-

subroutine tao_print_command_line_info

use sim_utils

implicit none

character(40), parameter :: r_name = 'tao_print_command_line_info'

!

call out_io (s_blank$, r_name, [ &
        'Syntax:                                                                                      ', &
        '  <path-to-tao-exe-directory>/tao {OPTIONS}                                                  ', &
        'Options are:                                                                                 ', &
        '  -beam <beam_file>                 # Beam init params (beam size, starting point, etc.)     ', &
        '  -beam_init_file_name <file_name>  # File with beam initial particle positions              ', &
        '  -beam_all <all_beam_file>         # Beam info from previous tracking                       ', &
        '  -building_wall <wall_file>        # Define the building tunnel wall                        ', &
        '  -color_prompt                     # Set color of prompt string to blue                     ', &
        '  -data <data_file>                 # Define data for plotting and optimization              ', &
        '  -debug                            # Debug mode for Wizards                                 ', &
        '  -disable_smooth_line_calc         # Disable the smooth line calc used in plotting          ', &
        '  -geometry <width>x<height>        # Plot window geometry (pixels)                          ', &
        '  -help                             # Display this list of command line options              ', &
        '  -hook_init_file <init_file>       # Init file for hook routines (Default = tao_hook.init)  ', &
        '  -init <tao_init_file>             # Tao init file                                          ', &
        '  -lat <bmad_lattice_file>          # Bmad lattice file                                      ', &
        '  -lat xsif::<xsif_lattice_file>    # XSIF lattice file                                      ', &
        '  -log_startup                      # Write startup debugging info                           ', &
        '  -no_stopping                      # For debugging: Prevents Tao from exiting on errors     ', &
        '  -noinit                           # Do not use Tao init file                               ', &
        '  -noplot                           # Do not open a plotting window                          ', &
        '  -plot <plot_file>                 # Plotting initialization file                           ', &
        '  -rf_on                            # Keep RF on (Default is to turn off)                    ', &
        '  -silent_run                       # Suppress terminal output when running a command file?  ', &
        '  -slice_lattice <ele_list>         # Discards elements from lattice that are not in the list', &
        '  -startup <starup_command_file>    # Commands to run after parsing Tao init file            ', &
        '  -var <var_file>                   # Define variables for plotting and optimization         '])

end subroutine tao_print_command_line_info

