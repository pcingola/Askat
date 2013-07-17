#!/bin/sh

# Build
rm -vf Askat.jar
cd $HOME/workspace/Askat/
ant

# Run script
cd -
./scripts/geno_cov.sh
#./scripts/karim4k.sh
#./scripts/karim4k_intervals.sh
#./scripts/t2d2.sh
