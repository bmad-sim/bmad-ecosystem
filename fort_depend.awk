######################################################################
#                                                                    #
# A fortran 90 dependency parser                                     #
#                                                                    #
# Author:  M. Palmer   7/27/00                                       #
#                                                                    #
# Modifications:  MAP  02/22/01  - Fix depend_match function         #
#                                - Fix parsing of multiple use       #
#                                  statements on a single line.      #
#                 MAP  02/23/01  - Set up paths to modules so that   #
#                                  the $(localmod) is included.      #
#                 MAP  09/23/01  - Set up paths to modules in the    #
#                                  $(MOD_OUT_DIR) area.              #
#                                                                    #
# WARNINGS:  Never loop over a variable in a main check that calls a #
#            function that loops over the same variable.  Variables  #
#            have global extent in awk (or is this a bug?!?!?)!      #
#                                                                    #
# FIX NEEDED:  include statements need to have their paths searched  #
#              properly.  makedepend would do this properly if all   #
#              include files were handled by #include                #
######################################################################

BEGIN{
# Initializations
	START = 0;
	IN_INTBLK = 0;
  target = "";
	inc_depend = "";
  mod_depend = "";

# Ignore upper versus lower case since we are parsing Fortran
	IGNORECASE = 1;

# DEBUG flag
DEBUG = 0;
}

############################
# Look for "use" statments #
############################
/use/ {

	if (DEBUG == 1) printf("In use block:  %s\n",$0);

# Strip comments from the directive line
	stripcomments($0);
	
# Parse the directives from the non-comment portion of the current line
	ndir = split(dir_line,direct,";");
	
	if (DEBUG == 1) printf("Number of directives:  %d\n",ndir);
	
# Loop over directives and determine module dependencies
	for(j=1; j <= ndir; j++) {

#    printf("Index is %d; Directive is %d\n", j, ndir);

		nfields = split(direct[j],fields," ");
    
	  if (DEBUG == 1) printf("Number of fields for directive %d:  %d\n",j,nfields);

		if ((nfields >= 2) && (tolower(fields[1]) == "use")) {
			fieldlen = length(fields[2]);
			
#     Test for comma at end of module name
			if (substr(fields[2],fieldlen,1) == ",") fieldlen -= 1;
			
#     Generate the module name and add to the dependency list
			module_name=substr(fields[2],1,fieldlen);
			module_lcname=tolower(module_name);
			module_file=sprintf("%s.mod",module_lcname);

    	if (DEBUG == 1) printf("Module Name:  %s\n", module_name);
			if (depend_match(mod_depend,module_file) == 0) {
				if (DEBUG == 1) printf("Add module:  %s\n", module_file);
				mod_depend = sprintf("%s %s", mod_depend, module_file);
				if (DEBUG == 1) printf("DEPEND:  %s\n",mod_depend);
			}
		}

	}
}

# Look for include files (Fortran-style only - let CPP handle any 
# #include statements) 
/include/ {

#	printf("In include block:  %s\n",$0);

# Strip comments from the directive line
	stripcomments($0);

# Parse the fields from the non-comment portion of the current line
	nfields = split(dir_line,fields," ");
	if ((nfields == 2) && (tolower(fields[1]) == "include")) {

#   Get the include file name and add to the dependency list
#   Make sure to remove quotes from around the include file name
		inclen = length(fields[2]);
    inc_name = substr(fields[2],2,inclen-2);
		if (depend_match(inc_depend,inc_name) == 0) {		
			inc_depend = sprintf("%s %s", inc_depend, inc_name);
		}
	}
}

END{
#	print "In end block";
#	printdepend(START,depend);
	printf("%s %s\n",mod_depend,inc_depend);
}

########################
# FUNCTION DEFINITIONS #
########################
# Function to print a dependency line
function printdepend(doprint,depend) {
	if (doprint == 1) {
		printf("%s\n",depend);
	}
	START = 1;
}


# Function to strip comments from a directive line
function stripcomments(line_in) {

# parse the line looking for a comment
	icomment = match(line_in,"!");

# Load dir_line with the comment-stripped line
# NOTE:  dir_line is a global variable
	if (icomment != 0) {
		dir_line = substr($0,1,icomment-1);
	} else {
		dir_line = line_in;
	}
}


# Function to search for a match to existing 
# dependencies in list
function depend_match(list,module) {

	if (DEBUG == 1) printf("In depend_match function\n");

# Loop through the individual items and compare
	ndepend = split(list,depends," ");
	for (i=1; i<= ndepend; i++) {
		
		if (DEBUG == 1) {
			printf("COMPARE:  %s  %s\n",module,depends[i]);
		}
		if (depends[i] == module) return 1;
	}
	
	return 0;
}
			
