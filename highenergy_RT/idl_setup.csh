# C shell commands to define IDL environment variables and aliases.
#
# This script can be used by IDL users who use csh as their interactive shell
# to define the environment variables and aliases required by IDL
# related commands (idl, idlde, idlhelp) if the symbolic links to
# the default directory (/usr/local/exelis/idl) are not being used.
#
# csh users should source idl_setup from their .cshrc files,
# using the following command:
#
#    source /Applications/exelis/idl83/bin/idl_setup
#

#setenv EXELIS_DIR /Applications/exelis
#setenv IDL_DIR /Applications/exelis/idl83
setenv EXELIS_DIR /usr/local/exelis
setenv IDL_DIR /usr/local/exelis/idl85
alias exelislicense $IDL_DIR/bin/exelislicense
alias idl $IDL_DIR/bin/idl
if ( -x $IDL_DIR/bin/idlde ) alias idlde $IDL_DIR/bin/idlde
if ( -x $IDL_DIR/bin/idlhelp ) alias idlhelp $IDL_DIR/bin/idlhelp
if ( -x $IDL_DIR/bin/idlrpc ) alias idlrpc $IDL_DIR/bin/idlrpc

