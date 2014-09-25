#!/bin/bash
# Setup script to set the appropriate environment variaibles for nustar_stray_light
#
# . setup_stray_light.sh
#
# Brian Grefenstette, August 2014, Caltech
#

#############################################

# Change the following lines for your own build!!! ###

# Set the IDL startup script
STRAY_DIR=$HOME/science/local/git/nustar_stray_light


# Set the path tot eh nustar-idl libraries:
export NUSTAR_IDL_ROOT=$HOME/science/local/git/nustar-idl

#############################################


# Do not change below.

# Check to make sure this the startup file and libraries exists:

export IDL_STARTUP=$STRAY_DIR/startup.pro
if [ ! -f $IDL_STARTUP ]
then
    echo "Change the IDL_STARTUP environment variable in $0 first!"
    echo "$IDL_STARTUP does not exist!"
    exit 1
fi

if [ ! -d $NUSTAR_IDL_ROOT ]
then
    echo "Set the NUSTAR_IDL_ROOT environment variable in $0 first!"
    echo "$NSUTAR_IDL_ROOT does not exist!"
    exit 1
fi

