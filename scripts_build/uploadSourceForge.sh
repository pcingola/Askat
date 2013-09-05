#!/bin/sh -e

#-------------------------------------------------------------------------------
#
# The code is in GIT, but since we cannot upload binary distribution (GIT 
# only allows code), we upload the binaries to SourceForge.
#
# We'll probably end up moving the whole project to SourceForge...
#
#																Pablo Cingolani
#-------------------------------------------------------------------------------

# Compile
./scripts_build/make.sh

# Go to home dir
cd $HOME/askat

# Create ZIP file
./scripts_build/zip.sh 

# Upload (use SnpEff site)
scp $HOME/askat.zip pcingola,snpeff@frs.sourceforge.net:/home/frs/project/s/sn/snpeff/

