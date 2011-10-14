#--------------------------------------------------------------
# acc_vars_csh
#
# Wrapper script for C-shell compatibility with the ACCelerator
# libraries management and software build system environment
# setup script.
# 
# See the file acc_vars_sh for usage notes.
#--------------------------------------------------------------

# Temporary for cleanliness.
source ~cesrulib/bin/clean_env.csh


setenv SETUP_SCRIPTS_DIR /nfs/acc/libs/util

#---------------------------------------------------------------
# Filename constants for obtaining environment differences.
#---------------------------------------------------------------
setenv ENV_USE_SNAPSHOTS "Y"
setenv ENV_SETUP_FILE "${HOME}/.ACC_env_var_set.tmp"
setenv ENV_SNAPSHOT_1 "${HOME}/.ACC_env_snapshot_1.tmp"
setenv ENV_SNAPSHOT_2 "${HOME}/.ACC_env_snapshot_2.tmp"


#---------------------------------------------------------------
# The main environment setup script, in Bourne shell syntax.
# All side effects are generated in a subshell and logged to
# disk.  The changes effected are then sourced into the
# environment that started _this_ script.
#---------------------------------------------------------------
`which bash` ${SETUP_SCRIPTS_DIR}/acc_vars.sh
unsetenv SET_SCRIPTS_DIR


#---------------------------------------------------------------
# Incorporate the difference in the representative environment
# provided by the above script into the active C-shell variant.
# For each line that has changed or that is new in the most
# recent environment snapshot, add quotes to the value of
# each variable's contents to safeguard against whitespace.
#
#    For environment snapshots generated in Bourne compatible
#    shells being converted to syntax appropriate for 
#    C-type shells:
#
#  1st 'sed' - replace  (=)   with  ( ")
#  2nd 'sed' - replace  (> )  with  (setenv )
#  3rd 'sed' - replace   LF   with  (")LF
#---------------------------------------------------------------
diff ${ENV_SNAPSHOT_1} ${ENV_SNAPSHOT_2} | grep \> | sed 's/=/ \"/g' | sed 's/> /setenv /g' | sed 's/$/\"/g' > ${ENV_SETUP_FILE}
source ${ENV_SETUP_FILE}


unsetenv ENV_SETUP_FILE
unsetenv ENV_SNAPSHOT_1
unsetenv ENV_SNAPSHOT_2
unsetenv ENV_USE_SNAPSHOTS
unsetenv ACC_RELEASE_REQUEST
