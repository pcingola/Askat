#!/usr/bin/perl


while( $l = <STDIN> ) {
	chomp $l;
	@t = ();
	@t = split / /, $l;

	
	print "$t[0] ";

	for( $i=1 ; $i <=3  ; $i++ )	{ 
		if( $t[$i] > 0 )	{ print "$t[0]_$t[$i] "; }
		else				{ print "0 "; }
	}

	print "$t[4] ";
	print "$t[5] ";

	for( $i=6 ; $i <= $#t ; $i++ )	{ 

		# Show mutation
		if( $t[$i] == 0 )		{ print "A A"; }
		elsif( $t[$i] == 1 )	{ 
			if( rand() < 0.5 )	{ print "A T"; }
			else				{ print "T A"; }
		}
		elsif( $t[$i] == 2 )	{ print "T T"; }


		# Show space or newline
		if ($i < $#t)	{ print " "; }
		else 			{ print "\n"; }
	}
}
