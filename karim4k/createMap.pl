#!/usr/bin/perl


for( $i=1 ; $i <= 4000 ; $i++ ) { 
	$pos = $i * 100;
	print "1 snp_$i 0 $pos\n";
}
