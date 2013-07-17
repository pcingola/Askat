
files <- c('ADCY5_D.txt')
files <- c('ADCY5_D.txt',  'ADCY5_M.txt',  'ADCY5_P.txt',  'TCF7L2_D.txt',  'TCF7L2_M.txt',  'TCF7L2_P.txt')

IDX_TH <- 100;	# All values below this index should be true positives

for( file in files ) {
	pngFile <- paste( file , '.png', sep="");
	png(pngFile, width = 1024/2, height = 1024);

	par( mfrow=c(2,1) );

	d <- read.table(file, header = FALSE, col.names = c('pos' , 'pvalue'));
	cat('\n', file, '\tNumber of pvalues:', length(d$pvalue), '\n');

	#---
	# Convert zero p-values to non-zero values
	#---
	pval <- d$pvalue
	idx <- 1:length(pval)
	nonNull <- pval > 0
	minp <- min( pval[ nonNull ] )
	pval[ ! nonNull ] <- minp;

	# Plot pvalues
	title <- paste( file, 'p-values [db]');
	q <- -10 * log(pval) / log(10);
	plot( idx, q, main=title, xlab='Index', ylab='p-values [db]' )
	abline( v=IDX_TH , col='blue' )

	th <- mean(q) + 3 * sd(q)
	abline( h=th , col='red' )

	# Stats
	tp <- (q >= th) & (idx <= IDX_TH)
	fn <- (q < th) & (idx <= IDX_TH)
	fp <- (q >= th) & (idx > IDX_TH)
	tn <- (q < th) & (idx > IDX_TH)
	cat('True positives  :', sum(tp), '\n')
	cat('False positives :', sum(fp), '\n')
	cat('False negatives :', sum(fn), '\n')
	cat('True negatives  :', sum(tn), '\n')
	points( idx[tp], q[tp], col='green')
	points( idx[fn], q[fn], col='blue')
	points( idx[fp], q[fp], col='red')

	# QQ-plot
	title <- paste( file, ' QQ-plot');
	y <- sort(pval)
	n <- length(y);
	x <- (1:n) / (n+1)
    lx <- -log10( x )
    ly <- -log10( y )
	range <- c(0 , max(lx, ly) );
	plot( lx, ly, xlim=range, ylim=range, main=title, xlab="-Log[ rank / (N+1) ]", ylab="-Log[ p ]" );
	abline( 0 , 1 , col='red');

	dev.off();
}
