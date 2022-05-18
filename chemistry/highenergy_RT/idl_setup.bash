# Bash shell commands to define IDL environment variables and aliases.
#
# This script can be used by IDL users who use Bash as their interactive shell
# to define the environment variables and aliases required by IDL
# related commands (idl, idlde, idlhelp) if the symbolic links to
# the default directory (/usr/local/harris/idl) are not being used.
#
# Bash users should run idl_setup from their .profile file 
# using the following command:
#
#    . /home/kschwarz/idl87/bin/idl_setup.bash
#
EXELIS_DIR=/home/kschwarz
IDL_DIR=/home/kschwarz/idl87
export IDL_DIR EXELIS_DIR
alias harrislicense=$IDL_DIR/bin/harrislicense
alias idl=$IDL_DIR/bin/idl
if [ -x $IDL_DIR/bin/idlde ]; then
  alias idlde=$IDL_DIR/bin/idlde
fi
if [ -x $IDL_DIR/bin/idlhelp ]; then
  alias idlhelp=$IDL_DIR/bin/idlhelp
fi
if [ -x $IDL_DIR/bin/idlrpc ]; then
  alias idlrpc=$IDL_DIR/bin/idlrpc
fi
if [ -x $IDL_DIR/bin/idltaskengine ]; then 
  alias idltaskengine=$IDL_DIR/bin/idltaskengine
fi

