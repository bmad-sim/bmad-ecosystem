$!#####################################################################
$!
$! Builds cesr-release on VMS.
$!
$! Usage:  @do_release_vms <dir> <rel> <pkg>
$!
$! Example:
$!  @do_release_vms [CESR_LIBS.CESR_2004_0621_D] CESR_2004_0621_D packages_2004_0608_D
$!
$!#####################################################################
$ tp1=0
$ tp2=0
$ tp3=0
$ if ("''p1'".eqs."") then tp1 = 1
$ if ("''p2'".eqs."") then tp2 = 1
$ if ("''p3'".eqs."") then tp3 = 1
$ if (tp1.or.tp2)
$ then  
$   write sys$output "Incorrect inputs specified..."
$   exit
$ endif
$ par2="PBASE=[CESR_LIBS.''p3'"
$ if (tp3) then par2="PBASE=[CESR_LIBS.PACKAGES"
$!
$!--------------------------------------------------
$!  Move to the specified directory
$!--------------------------------------------------
$!
$ crsdir=F$DIRECTORY()
$ set def 'p1'
$!
$!--------------------------------------------------
$! building c_utils
$!--------------------------------------------------
$!
$ par1="RELEASE=''p2'"
$ gmake --directory=[.cvssrc.c_utils] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building cesr_utils
$!--------------------------------------------------
$!  
$ gmake --directory=[.cvssrc.cesr_utils] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building RECIPES_F_90_LEPP
$!--------------------------------------------------
$!  
$ gmake --directory=[.cvssrc.RECIPES_F_90_LEPP] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building dcslib
$!--------------------------------------------------
$!  
$ gmake --directory=[.cvssrc.dcslib] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building bmad
$!--------------------------------------------------
$!  
$ gmake --directory=[.cvssrc.bmad] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building CesrBPM
$!--------------------------------------------------
$!  
$ gmake --directory=[.cvssrc.CesrBPM] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building bmadz
$!--------------------------------------------------
$!  
$ gmake --directory=[.cvssrc.bmadz] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building odysseus
$!--------------------------------------------------
$!  
$ gmake --directory=[.cvssrc.odysseus] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building toad
$!--------------------------------------------------
$!  
$! gmake --directory=[.cvssrc.toad] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building nonlin_bpm
$!--------------------------------------------------
$!  
$ gmake --directory=[.cvssrc.nonlin_bpm] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building BeamInstSupport
$!--------------------------------------------------
$!  
$ gmake --directory=[.cvssrc.BEAMINSTSUPPORT] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building MPMnet
$!--------------------------------------------------
$!  
$ gmake --directory=[.cvssrc.MPMnet] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$ set def 'crsdir'
$ EXIT















