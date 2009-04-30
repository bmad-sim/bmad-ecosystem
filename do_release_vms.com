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
$ write sys$output "========== Building c_utils =========="
$ par1="RELEASE=''p2'"
$ gmake --directory=[.src.c_utils] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building RECIPES_F_90_LEPP
$!--------------------------------------------------
$!  
$ write sys$output "========== Building RECIPES_F_90_LEPP =========="
$ gmake --directory=[.src.RECIPES_F_90_LEPP] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building sim_utils
$!--------------------------------------------------
$!
$ write sys$output "========== Building sim_utils ==========="
$ gmake --directory=[.src.sim_utils] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building MPMnet
$!--------------------------------------------------
$!  
$ write sys$output "========== Building MPMnet =========="
$ gmake --directory=[.src.mpmnet] "''par1'" "''par2'"
$ gmake --directory=[.src.mpmnet] --file=MAKEFILE_MPMEQ.VMS  "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building bmad
$!--------------------------------------------------
$!  
$ write sys$output "========== Building bmad =========="
$ gmake --directory=[.src.bmad] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building cesr_utils
$!--------------------------------------------------
$!
$ write sys$output "========== Building cesr_utils =========="  
$ gmake --directory=[.src.cesr_utils] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building mpm_utils
$!--------------------------------------------------
$!  
$ write sys$output "========== Building mpm_utils =========="
$ gmake --directory=[.src.mpm_utils] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building nonlin_bpm
$!--------------------------------------------------
$!  
$ write sys$output "========== Building nonlin_bpm =========="
$ gmake --directory=[.src.nonlin_bpm] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building toad
$!--------------------------------------------------
$!  
$! write sys$output "========== Building toad =========="
$! gmake --directory=[.src.toad] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building odysseus
$!--------------------------------------------------
$!  
$! write sys$output "========== Building odysseus =========="
$! gmake --directory=[.src.odysseus] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building CesrBPM
$!--------------------------------------------------
$!  
$ write sys$output "========== Building CesrBPM =========="
$ gmake --directory=[.src.CesrBPM] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building dcslib
$!--------------------------------------------------
$!  
$! write sys$output "========== Building dcslib =========="
$! gmake --directory=[.src.dcslib] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building pecklib
$!--------------------------------------------------
$! write sys$output "========== Building pecklib =========="
$! gmake --directory=[.src.pecklib] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building bmadz
$!--------------------------------------------------
$!  
$! write sys$output "========== Building bmadz =========="
$! gmake --directory=[.src.bmadz] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$! building BeamInstSupport
$!--------------------------------------------------
$!  
$! write sys$output "========== Building BeamInstSupport =========="
$! gmake --directory=[.src.BEAMINSTSUPPORT] "''par1'" "''par2'"
$!
$!--------------------------------------------------
$ set def 'crsdir'
$ EXIT















