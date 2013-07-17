#!/bin/sh

SNPEFF_SCRIPTS=$HOME/snpEff/scripts/

#---
# Show a report on each file
#---
echo "<table> "
echo "<tr> <th> Analysis </th> <th> P-Values </th> <th> P-Values [Db] </th> <th> QQ-plot </th> <th> Hitogram p-values [Db] </th> </tr>"
for f in *.out
do
	echo "<tr> "

	echo "<td> <b> $f </b> </td>"

	# Get a file showing "chr \t pos \t pvalue", sorted by chr:pos
	cat $f | grep ASKAT_RESUTS: | cut -f 8,12 | tr ":" "\t" | sed "s/-/|/" | tr "|" "\t" | cut -f 1,2,5 | sort -k1n -k2n > pvalues.txt
	cat pvalues.txt | cut -f 3 | $SNPEFF_SCRIPTS/db.pl | grep -v NA > pvalues_db.txt

	#---
	# Show P-value 
	#---
	cat pvalues.txt | cut -f 3 | $SNPEFF_SCRIPTS/plot.pl "p-value_$f" > /dev/null
	echo "<td> <a href=p-value_$f.png> <img width=300 src=p-value_$f.png> </a> </td>"

	cat pvalues_db.txt | $SNPEFF_SCRIPTS/plot.pl "p-value_db_$f" > /dev/null
	echo "<td> <a href=p-value_db_$f.png> <img width=300 src=p-value_db_$f.png> </a> </td>"

	#---
	# QQ plot
	#---
	cat pvalues.txt | cut -f 3 | $SNPEFF_SCRIPTS/qqplot.pl "qq_$f" > /dev/null
	echo "<td> <a href=qq_$f.png> <img width=300 src=qq_$f.png> </a> </td>"

	#---
	# Histogram
	#---
	cat pvalues_db.txt | $SNPEFF_SCRIPTS/hist.pl "hist_$f" > /dev/null
	echo "<td> <a href=hist_$f.png> <img width=300 src=hist_$f.png> </a> </td>"

	echo "</tr> "
done
echo "</table> "

#rm pvalues.txt pvalues_db.txt

