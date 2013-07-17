#!/usr/bin/perl

#-------------------------------------------------------------------------------
#
# Find the minimum p-value and position from ASKAT results
#
# Udage: cat resutls.txt | ./minPval.pl
#
#																Pablo Cingolani
#-------------------------------------------------------------------------------

$pvalMin = 1;

# Read STDIN
while( $l = <STDIN> ) {
	(@t) = split /\t/, $l;

	# Parse p-value and position
	($pval, $pos) = ($t[2], $t[4]);

	# Find minimum p-value
	if( $pval < $pvalMin ) { ($pvalMin, $posMin) = ($pval, $pos); }
}
print "pValue\t$pvalMin\tpos\t$posMin\n";
