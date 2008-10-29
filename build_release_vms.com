
$ TRIGGER_FILE = "cesr_mount:[000000.bctl]cesr_VMS_builds.flag" !this is the folder where to look for the file
$ LOG_DIRECTORY= "cesr2f$dkb200:[cesr_libs.log]"
$ LIBS_DIR = "cesr2f$dkb200:[cesr_libs" !where the zip/unzip files are. Note, "]" is missing
$ PACKAGES = "packages_2008_0601_d"
$ TEMP_OUTPUT_FILE = "cesr2f$dkb200:[cesr_libs]temp_output.dat" !temporary file, do not alter
$ SPACE_NEEDED = 1600000 !Require 1.6 Gb for the VMS build (need both zip and build, so double)

$ start:
$     file = f$search(TRIGGER_FILE)
$     if file .nes. "" !means file is present (make sure to delete file afterwards)
$     then
$         write sys$output "File found "+f$time()
$     !-------------------------------------------------------------------------------------
$     ! Step 1: Get the information about available space, determine if enough space
$         pipe show device cesr2f$dkb200 > 'TEMP_OUTPUT_FILE'
$         open/error=open_error temp_output 'TEMP_OUTPUT_FILE'
$         read_loop_1:
$             read/end_of_file=end_read_1 temp_output line
$             if f$locate("CESR2F",line) .eq. f$length(line)
$             then
$                 goto read_loop_1
$             else !means this is the line
$                 availablespace = f$extract(61,10,line)
$                 availablespace = f$edit(availablespace,"TRIM")
$                 freespace = f$integer(availablespace)
$                 write sys$output "Available space: "'freespace'
$                 write sys$output "Space required for build: "'SPACE_NEEDED' 
$                 goto end_read_1
$             endif
$         end_read_1:
$             close temp_output
$             delete 'TEMP_OUTPUT_FILE';*
$             if (freespace .lt. SPACE_NEEDED)
$             then
$                 write sys$output "Not enough space. Exiting"
$                 exit
$             endif
$     !------------------------------------------------------------------------------------
$     !Step 2: Extract the release name from the trigger file.
$         open/error=open_error temp_output 'TRIGGER_FILE'
$         read_loop_2:
$             read/end_of_file=end_read_2 temp_output line
$             release_name = f$edit(line,"TRIM") !assumes first line has release name
$             goto read_loop_2 !This is in case more parameters needed, then 
$         end_read_2:
$             close temp_output
$             delete 'TRIGGER_FILE';*
$     !-------------------------------------------------------------------------------------
$     !Step 3: Run the unzip command, save output to a file and look for errors
$         pipe unzip 'LIBS_DIR'.zip]'release_name'.zip /dir='LIBS_DIR'.'release_name'] > 'TEMP_OUTPUT_FILE'
$         open/error=open_error temp_output 'TEMP_OUTPUT_FILE'
$         read_loop_3:
$             read/end_of_file=end_read_3 temp_output line
$             if f$locate("ERROR",line) .eq. f$length(line)
$             then
$                 goto read_loop_3
$             else
$                 write sys$output "ERROR detected during UNZIP command. Exiting"
$                 exit
$             endif
$         end_read_3:
$             close temp_output
$             delete 'TEMP_OUTPUT_FILE';*  
$     !-------------------------------------------------------------------------------------
$     !Step 4: Run the build command, save output to a file and look for errors
$         pipe @'LIBS_DIR'.'release_name'.util]do_release_vms [cesr_libs.'release_name'] 'release_name' 'PACKAGES' > 'LOG_DIRECTORY''release_name'.log
$         open/error=open_error temp_output 'LOG_DIRECTORY''release_name'.log
$         read_loop_4:
$             read/end_of_file=end_read_4 temp_output line
$             if f$locate("ERROR",line) .eq. f$length(line)
$             then
$                 goto read_loop_4
$             else
$                 write sys$output "A build ERROR detected. Exiting"
$!                 exit !uncomment this line if error detection OK
$             endif
$         end_read_4:
$             close temp_output

$     write sys$output "Build completed successfully at "+f$time()
$     endif
$     goto start
$ open_error:
$     write sys&output "Unable to open the file"
$     exit

