$!################################################################
$!
$! Check for the existence of a directory and create as necessary
$!
$! Mark Palmer  Jan 3, 2002
$!
$!################################################################
$ set noverify
$!
$! Check for valid input
$!
$ if ("''p1'" .eqs. "") 
$ then 
$   write sys$output "No directory specified..."
$   exit
$ endif
$!
$!  Generate the string to parse for directory name
$!
$ inlen = f$length(p1)
$ start = f$locate("[",p1)
$ dirname = f$extract(start,inlen-1,p1)
$!
$! Get the base directory name that we want to create
$!
$dir_name:
$  pos = f$locate(".",dirname)
$  len = f$length(dirname)
$!  write sys$output pos," ",len
$  if pos .LT. len 
$  then 
$    dirname = f$extract(pos+1,len,dirname)
$!    write sys$output "DEBUG dirname:  ",dirname
$    goto dir_name
$  else
$    end     = inlen - len
$    path    = f$extract(0,inlen-len-2,p1) + "]"
$    dirname = dirname + ".dir"
$    search_name = path + dirname
$  endif
$!
$! write sys$output "Checking for directory file:  ", search_name
$!
$! Check for directory and create if not present.  Note that this 
$! com file presumes that higher directories in the hierarchy are
$! already present!
$!
$ if f$search(search_name) .eqs. ""
$ then
$   write sys$output "Creating:  ", search_name
$   create/dir 'p1'
$ endif
$!
$! $Id$
$!
$! $Log$
$! Revision 1.1  2002/01/09 16:29:57  cesrulib
$! First pass at VMS library build infrastructure
$!
$!
