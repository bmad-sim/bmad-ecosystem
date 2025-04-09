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
        'Syntax:                                                                                         ', &
        '  tao {OPTIONS}                                      ! If the tao exe is in your PATH           ', &
        '  <path-to-tao-exe-directory>/tao {OPTIONS}                                                     ', &
        'Options are:                                                                                    ', &
        '  -beam_file <file_name>               # File containing the tao_beam_init namelist.            ', &
        '  -beam_init_position_file <file_name> # File containing initial particle positions.            ', &
        '  -building_wall_file <file_name>      # Define the building tunnel wall                        ', &
        '  -command <command_string>            # Commands to run after startup file commands            ', &
        '  -data_file <file_name>               # Define data for plotting and optimization              ', &
        '  -debug                               # Debug mode for Wizards                                 ', &
        '  -disable_smooth_line_calc            # Disable the smooth line calc used in plotting          ', &
        '  -external_plotting                   # Tells Tao that plotting is done externally to Tao.     ', &
        '  -geometry <width>x<height>           # Plot window geometry (pixels)                          ', &
        '  -help                                # Display this list of command line options              ', &
        '  -hook_init_file <file_name>          # Init file for hook routines (Default = tao_hook.init)  ', &
        '  -init_file <file_name>               # Tao init file                                          ', &
        '  -lattice_file <file_name>            # Bmad lattice file                                      ', &
        '  -log_startup                         # Write startup debugging info                           ', &
        '  -no_stopping                         # For debugging: Prevents Tao from exiting on errors     ', &
        '  -noinit                              # Do not use Tao init file.                              ', &
        '  -noplot                              # Do not open a plotting window                          ', &
        '  -nostartup                           # Do not open a startup command file                     ', &
        '  -no_rad_int                          # Do not do any radiation integrals calculations.        ', &
        '  -plot_file <file_name>               # Plotting initialization file                           ', &
        '  -prompt_color <color>                # Set color of prompt string. Default is blue.           ', &
        '  -reverse                             # Reverse lattice element order?                         ', &
        '  -rf_on                               # Use "--rf_on" to turn off RF (default is now RF on)    ', &
        '  -quiet <level>                       # Suppress terminal output when running a command file?  ', &
        '                                       #  Levels: "all" (default), "warnings".                  ', &
        '  -slice_lattice <ele_list>            # Discards elements from lattice that are not in the list', &
        '  -start_branch_at <ele_name>          # Start lattice branch at element.                       ', &
        '  -startup_file <file_name>            # Commands to run after parsing Tao init file            ', &
        '  -symbol_import                       # Import symbols defined in lattice files(s)?            ', &
        '  -var_file <file_name>                # Define variables for plotting and optimization         ', &
        'Note: Use "--" instead of "-" prefix to negate an option. See the Tao manual for details.       '])

end subroutine tao_print_command_line_info

