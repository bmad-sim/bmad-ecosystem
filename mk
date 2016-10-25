#!/bin/sh
#
# Frontend script to do a production build.
# This script removes the need for symlinks, for Windows builds.

export THIS_SCRIPT='mk'
${ACC_ROOT_DIR}/util/mk-mkd $1

#
