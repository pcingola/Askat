#!/bin/sh

for out in *.out
do
	base=`basename $out var0.01.out`
	cat $out | grep ASKAT_RESUTS: | cut -f 8,12 | sed "s/-/:/" | tr ":" "\t" | cut -f 2,5 | sort -n > $base.txt
done
