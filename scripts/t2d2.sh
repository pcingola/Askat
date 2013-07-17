#!/bin/sh

# Delete old files
rm -rvf replicate.head.block.*

# Run wrapper
java -Xmx4G -jar Askat.jar \
	-v \
	-noDep \
	-kin all \
	-sb 100 \
	-maxMaf 1.0 \
	-p 5 \
	replicate.head \
	2>&1 | tee replicate.head.out

