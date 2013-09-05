#!/usr/bin/perl

#-------------------------------------------------------------------------------
#Make sure TPED file has unique IDs
#
#
#                                                               pcingola 2013
#-------------------------------------------------------------------------------

for( $lineNum=1 ; $l = <STDIN> ; $lineNum++ ) {
	chomp $l;
	@t = ();
	@t = split /\s+/, $l;
	$t[1] = "id_$lineNum";			# Change ID to 'unique' one.
	print join("\t", @t) . "\n";
}

