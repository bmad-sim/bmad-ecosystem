$!#####################################################################
$!
$! Builds packages-release on VMS.
$!
$! Usage:  @do_packages_vms <dir> <rel> <MAKE_SRC>
$!
$! Example:
$!  @do_packages_vms [CESR_LIBS.packages_2004_0621_D] PACKAGES_2004_0621_D [CESR_LIBS.current.GMAKE]
$!
$!#####################################################################
$ tp1=0
$ tp2=0
$ tp3=0
$ if ("''p1'".eqs."") then tp1 = 1
$ if ("''p2'".eqs."") then tp2 = 1
$ if ("''p3'".eqs."") then tp3 = 1
$ if (tp1.or.tp2.or.tp3)
$ then  
$   write sys$output "Incorrect inputs specified..."
$   exit
$ endif
$!--------------------------------------------------
$!set up some parameters
$!--------------------------------------------------
$ par1="RELEASE=''p2'"
$ if ("''p3'".eqs."")
$ then
$   par2=""
$ else
$   par2="MAKE_SRC=''p3'"
$ endif
$ par4="[CESR_LIBS.current.util]makedir"
$!
$!
$!--------------------------------------------------
$!  Move to the specified directory
$!--------------------------------------------------
$!
$ crsdir=F$DIRECTORY()
$! write sys$output "bubub ''crsdir'"
$ set def 'p1'
$!
$!--------------------------------------------------
$! building PGPLOT
$!--------------------------------------------------
$!
$ @'par4' [.PGBIN]
$ set def [.PGBIN]
$ backup/replace [-.PGPLOT]drivers.list []
$ @[-.PGPLOT.SYS_VMS]install [-.PGPLOT]
$ @[-.PGPLOT.SYS_VMS]install [-.PGPLOT] cpg
$ set def 'p1'
$!
$!--------------------------------------------------
$! building forest
$!--------------------------------------------------
$!
$ gmake --directory=[.forest] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building RECIPES_F_90
$!--------------------------------------------------
$!  
$ gmake --directory=[.RECIPES_F_90] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building  RECIPES_C_ANSI
$!--------------------------------------------------
$!  
$ gmake --directory=[.RECIPES_C_ANSI] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building xsif
$!--------------------------------------------------
$!  
$ gmake --directory=[.xsif] "''par1'" "''par2'"
$!
$ set def 'crsdir'
$ EXIT






