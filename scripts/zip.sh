#!/bin/sh

cd ..

rm -vf askat.zip

zip -r askat.zip \
	askat/Askat.jar \
	askat/r/ \
	askat/scripts/ \
	askat/data/ \
	askat/bed/ \
	askat/karim4k.* 

cd -

