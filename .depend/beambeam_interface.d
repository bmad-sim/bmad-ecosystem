# DO NOT DELETE

$(LIBRARY)(beambeam_interface$(PREC)$(DBG).o) beambeam_interface.d: /nfs/acc/libs/Linux_i686_intel/devel/src/include/CESR_platform.h
$(LIBRARY)(beambeam_interface$(PREC)$(DBG).o) beambeam_interface.d $(MOD_OUT_DIR)/beambeam_interface.mod:  /nfs/acc/libs/Linux_i686_intel/devel/modules/bmad_struct.mod /nfs/acc/libs/Linux_i686_intel/devel/modules/bmad_interface.mod ../modules/bmadz_mod.mod ../modules/bmadz_interface.mod $(MOD_OUT_DIR)/bsim_interface.mod $(MOD_OUT_DIR)/scan_parameters.mod
