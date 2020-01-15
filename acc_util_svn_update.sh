#!/bin/bash

# Script to update the ACC util directory from the ACC Repository 

#set -x

# Function for kerberos ticket expiration issue
func_kerberos_setup () {
    kdestroy -A
    kinit -k -t /home/$USER/etc/$USER-keytab $USER
    klist
}

# Function to set the login enironment
func_acc_login_setup () {
case $(uname -n) in  
    cesr*) . /nfs/cesr/online/acc_control/program_info/shell/cesr_online.bashrc;;  
    *) . /nfs/acc/libs/cesr/cesr_online.bashrc;;  
esac 
}

# Function to update the directories
func_svn_up_dirs () {
    [ $(pwd) == "${ACC_UTIL}" ] \
	|| cd ${ACC_UTIL} \
	&& echo -e '\nUpdating '${ACC_UTIL}
    svn up

    [ $(pwd) == "${ACC_BUILD_SYSTEM}" ] \
	|| cd ${ACC_BUILD_SYSTEM} \
	&& echo -e '\nUpdating '${ACC_BUILD_SYSTEM}
    svn up
}

# Confirm that cesrulib is running the script
if [[ $(id) == *cesrulib* ]] ; then
    func_kerberos_setup
else    
    echo -e '\nThis script can only be run by cesrulib\n'
    exit 1
fi

# Check and enable the ACC environment
[ "${ACC_UTIL}" ] \
    || func_acc_login_setup

# Check the system and update for the Repository 
case $(uname -n) in 
    lnx6166*)
	func_svn_up_dirs
	[ $? == 0 ] && \
	    echo -e '\nThis script must now be run on a cesr system - e.g. "ssh cesrulib@cesrshell"\n'
	;;
    cesr*)
	func_svn_up_dirs
        [ $? == 0 ] && \
	    echo -e '\nThis script must now be run on lnx6166 - e.g. "ssh cesrulib@lnx6166"\nif not done already...\n'
	;;
    *)
	echo -e '\nThis script can only be run on lnx6166 or a cesr system.\n'
	exit 1
	;;
esac


