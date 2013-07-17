#!/bin/sh

cd ..
tar -cvhzf askat.tgz \
	askat/Askat.jar \
	askat/fastlmmc \
	askat/r/ \
	askat/scripts/ \
	askat/data/ \
	askat/karim4k.* 
cd -

