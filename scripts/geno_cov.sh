#!/bin/sh

# Delete old files
rm -rvf geno_cov.block.*

# Run askat
java -Xmx4G -jar Askat.jar \
	-v \
	-p 1 \
	-noDep \
	-kin all \
	-sb 20 \
	-maxMaf 1.0 \
	geno_cov \
	2>&1 | tee geno_cov.out

