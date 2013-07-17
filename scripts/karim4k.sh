#!/bin/sh

# Delete old files
rm -rvf karim4k.block.*

# Copy Karim's kinship matrix
cp karim4k/kinshipMatrix.karim.RData karim4k.block.1_100.kinship.RData

#	-d \
#	-p 2 \
# Run using block size 10
java -Xmx4G -jar Askat.jar \
	-v \
	-noDep \
	-kin all \
	-sb 10 \
	-maxMaf 1.0 \
	karim4k \
	2>&1 | tee karim4k.out

