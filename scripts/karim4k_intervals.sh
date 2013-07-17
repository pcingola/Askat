#!/bin/sh

# Copy Karim's kinship matrix
cp karim4k/kinshipMatrix.karim.RData karim4k.block.1_100.kinship.RData

#	-d \
# Run using block size 10
java -Xmx4G -jar Askat.jar \
	-v \
	-noDep \
	-kin all \
	-sb 10 \
	-maxMaf 1.0 \
	-i karim4k.bed \
	-p 2 \
	karim4k \
	2>&1 | tee karim4k.out

