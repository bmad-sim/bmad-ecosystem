BEGIN{
	IGNORECASE = 1;
  in_comment = 0;
}


# Look for Fortran programs
/program/ {
	if (tolower($1) == "program") {
    status = get_basename();
  }
}



# For C programs need to be careful about comments
/\/\*/ {
  in_comment = 1;
}

/\*\// {
  in_comment = 0;
}


# Look for C/C++ main functions
/main/ {

# Check for leading C++ style comment
# Check that not in a C comment
# Check for leading quote for a string
  if ((substr($1,1,2) != "//") && 
      (in_comment != 1)        &&
      (substr($1,1,1) != "\"")    ) {
    if (($1 == "main")  ||
        (substr($1,1,5) == "main(") ||
        ($2 == "main")  ||
        (substr($2,1,5) == "main(") ||
        ($3 == "main")  ||
        (substr($3,1,5) == "main(")   ) {

      file_type = get_filetype();
      if ((file_type == "c")   ||
          (file_type == "cc")  ||
          (file_type == "cpp") ||
          (file_type == "cxx")  ) {
        status = get_basename();
#      Diagnostic printout
#      print $0;
      }
    }
  }
}



# Function definitions

# Process FILENAME for basename
# Recall that fields returned from split are numbered 1...N
function get_basename() {
  nfields   = split(FILENAME,file_fields, /\//);
  if (nfields < 1) {
    return 1;
  }
  file_name = file_fields[nfields];
  ndummy    = split(file_name,base_fields, /\./);
  printf("%s ", base_fields[1]);	
  return 0;
}

# Process FILENAME for file type
# Recall that fields returned from split are numbered 1...N
function get_filetype() {
  nfields   = split(FILENAME,file_fields, /\//);
  if (nfields < 1) {
    return 1;
  }
  file_name = file_fields[nfields];
  nfields   = split(file_name,base_fields, /\./);
  if (nfields != 2) {
    return 1;
  }	
  return base_fields[2];
}
