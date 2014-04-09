#-------------------------------------------------------------------------------
# ASKAT Package
#
# By: Karim Oualkacha
# Modifications by: Pablo Cingolani
#
#
# 2013-01:
# --------
#
#	- FaST-LMM code moved to kinship.r
#
# Code review:
# ------------
#
#	- Joined all functions into one file (askar.r)
#
#	- Call to fastLmm should be in one function. It is now split in two or 
#     more, thus creating several "global" names and variables (not modular).
#
#	- In VC.FaST.LMM you write to fastlmmOutFileName and immediately read 
#     the same file. Why!?!? (you already have the data, why do you need 
#     to write/read the same file?).
#
#		# Write results to 'fastlmmOutFileName'
#		output = scan(paste(genoOutFileName, sep=""), what="character", sep="\n")
#		write.table(output[2], file = fastlmmOutFileName, quote = FALSE, sep = " ", eol = "\r\n", na = "NA", row.names = FALSE, col.names = FALSE)
#		# Read results from fastlmmOutFileName. Why do you write and then read the same file?
#		VC.FaST.LMM = read.table(fastlmmOutFileName, header=FALSE)
#
#		Changed that code by this one liner:
#
#		VC.FaST.LMM = read.table(genoOutFileName, header=TRUE)
#
#	- In VC.FaST.LMM: Side effect deletes a file created outside the 
#     function. This is bad practice, we should change it.
#
#	- Mixed assignment operators: '=' and '<-' are both used. R programmers 
#     are quite 'orthodox' and usually prefer '<-'.
#
#	- Changed 'system("rm ...") by 'unlink(...)' which is system independent.
#
#	- Changed 'paste("...", sep="")' by just "...". It has the same effect, so
#     why the 'paste' function anyway?
#
#	- In VC.FaST.LMM, 'ans' is created just to return a list (removed).
#
#	- When writing pheno.txt and TFAM files, the separator was "    " 
#     (i.e.four spaces). Changed it to single space. Apparently it only made 
#     files larger and marginally slower to parse (when running this 
#     millons of times, everything counts).
#
#	- Added return statement to clarify the code.
#
#	- Show ASKAT results in one line (using ASKAT_RESULTS) to make it easier
#	  for other programs to parse this output.
#
#
#																		2012
#-------------------------------------------------------------------------------

library(GenABEL)
library(nFactors)
library(CompQuadForm)

#-------------------------------------------------------------------------------
# ASKAT Main function
#
# ped    : Is the pedigree data file. It has subject IDs as first column (IDs 
#          should be differents for all subjects), phenotype as second column 
#          and region-based SNPs that will be analized together  
#
# fastlmm: Results form fastlmm program (see kinship.r)
#
#-------------------------------------------------------------------------------
ASKAT <- function(ped, fastlmm) {

	#####  STEP 1: Check for Missing Data and construction of pedigree without Missing Data   #####

	Y.trait = ped[,2]
	X = as.matrix(ped[,3:dim(ped)[2]])

	##### STEP 2: Under the NULL we call FaST-LMM to estimate the VCs ####
	estim.sigma.RG = as.numeric(as.character(fastlmm$nullGeneticVar));
	estim.sigma.e = as.numeric(as.character(fastlmm$nullResidualVar));
		
	S = fastlmm$S
	U = fastlmm$U
    #print(U[U!=0])
	##### STEP 3: Calculation of weights matrix W and the matrix K =GWG #####
	freq.MAF = apply(X, 2, mean)/2
	#Geno.no.sparse = which(!freq.MAF==0)
    #X = X[,Geno.no.sparse]
    #freq.MAF = apply(X, 2, mean)/2

	#Attempt to reduce time to O(N^2), using PCA trick
  
  	if( length(freq.MAF) == 1){
    	w <- dbeta(freq.MAF, 1, 25)
    	K.sqrt <- w * t(X)
  	} else {
    	w <- vector(length = length(freq.MAF))
    	for (i in 1:length(freq.MAF)){
    	  w[i] <- dbeta(freq.MAF[i], 1, 25)
    	}
    	w <- diag(w)
    	K.sqrt <- w %*% t(X)
  	}
	
##### STEP 2: ASKAT score test statistic calculations #####


	Gamma = estim.sigma.RG[1] / estim.sigma.e[1] 
	UT<-t(U)
	D.0 <- (Gamma * S) + diag(1, dim(X)[1], dim(X)[1])
	inv.sqrt.D.0 <- diag(1/sqrt(diag(D.0)))
	inv.D.0 = diag(1/diag(D.0))
	un.n <- c(rep(1,dim(U)[1]))
	Z <- 1/(( t(un.n) %*% (U %*% (inv.D.0 %*% (UT %*% un.n))))[1,1]) 
	X.tilde <- inv.sqrt.D.0 %*% (UT %*% un.n)
	Y.tilde <- inv.sqrt.D.0 %*% (UT %*% Y.trait)
	s2 = estim.sigma.e
	P.0.tilde = (diag(1, dim(U)[1], dim(U)[2])) - (X.tilde %*% Z) %*% ((t(X.tilde)))
	#K.tilde2 <- inv.sqrt.D.0 %*% (UT %*% (K %*% (U %*% inv.sqrt.D.0))) 
	#W1 <- P.0.tilde %*% K.tilde2 %*% P.0.tilde
	#PDU <- P.0.tilde %*% inv.sqrt.D.0 %*% UT
	#PDUT <- t(PDU)
	#W1 <- P.0.tilde %*% inv.sqrt.D.0 %*% UT %*% K %*% U %*% inv.sqrt.D.0 %*% P.0.tilde #55 s
	RM <- ((K.sqrt %*% U) %*% inv.sqrt.D.0) %*%  P.0.tilde# 44 s
	W <- RM %*% t(RM)

	#P.0.tilde is symmetric
	#res <- P.0.tilde %*% Y.tilde
	Q <- (t(Y.tilde) %*% t(RM) %*% ((RM %*% Y.tilde)))/(2 * s2)
	#print("time for Q computation, after associativity update")
	#print(proc.time()-ptm)

	#print("Q, Q2, Q-Q2"); print(c(Q, Q2, Q-Q2))


	#W1 = P.0.tilde %*% K.tilde
	#W1 = W1 %*% P.0.tilde/2
	#K.tilde = inv.sqrt.D.0 %*% UT

	#P.0.tilde = (diag(1, dim(U)[1], dim(U)[2])) - (X.tilde %*% Z) %*% ((t(X.tilde)))
	#W1 = P.0.tilde %*% K.tilde2 %*% P.0.tilde/2#UPD: on average 0.5 s gain
   
  	out = Get_PValue.Modif(W/2, Q)
  	pvalue.davies = out$p.value
  	lambda = out$lambda

	return( list(pvalue.ASKAT = pvalue.davies, Q.ASKAT = Q, Polygenic.VC = estim.sigma.RG, Env.VC = estim.sigma.e, lambda = lambda) );
}

#-------------------------------------------------------------------------------
# Get lambda
#-------------------------------------------------------------------------------
Get_Lambda <- function (K) {
	out.s <- eigen(K, symmetric = TRUE, only.values = TRUE)
	lambda1 <- out.s$values
	IDX1 <- which(lambda1 >= 0)
	IDX2 <- which(lambda1 > mean(lambda1[IDX1])/1e+05)
	if (length(IDX2) == 0) { fatalError("No Eigenvalue is bigger than 0!!") }
	lambda <- lambda1[IDX2]
	return( lambda );
}

#-------------------------------------------------------------------------------
# Get p-value
# Note: This function was taken from SKAT package
#-------------------------------------------------------------------------------
Get_PValue.Modif <- function(K, Q){
	lambda <- Get_Lambda(K)
	n1 <- length(Q)
	p.val <- rep(0, n1)
	p.val.liu <- rep(0, n1)
	is_converge <- rep(0, n1)
	for (i in 1:n1) {
		out <- davies(Q[i], lambda, acc = 10^(-6))
		p.val[i] <- out$Qq
		p.val.liu[i] <- liu(Q[i], lambda)
		is_converge[i] <- 1

		if (length(lambda) == 1) {
			p.val[i] <- p.val.liu[i]
		} else if (out$ifault != 0) {
			is_converge[i] <- 0
		}

		if (p.val[i] > 1 || p.val[i] < 0) {
			is_converge[i] <- 0
			p.val[i] <- p.val.liu[i]
		}
	}

	return(list(p.value = p.val, p.val.liu = p.val.liu, is_converge = is_converge, lambda = lambda))
}

#-------------------------------------------------------------------------------
# Fatal error
#-------------------------------------------------------------------------------
fatalError <- function(errStr)	{ 
	errStr <- paste('FATAL ERROR:',  errStr, "\n");
	stop( errStr ); 
}

#-------------------------------------------------------------------------------
# Main program
#-------------------------------------------------------------------------------

# Example of command line arguments
# You can just uncomment one of these lines for debugging.
#
# cmdLineArgs <- c( 'geno_cov.block.1_1.1.askat', 'geno_cov.tfam', 'geno_cov.block.1_1.kinship.RData', 20, './tmp_zzz', './fastlmmc', 'TRUE' )
# cmdLineArgs <- c( 'karim4k.block.1_100.1.askat', 'karim4k.tfam', 'karim4k.block.1_100.kinship.RData', 10, 'karim4k.block.1_100.1', 'fastlmmc', 'TRUE' )
#

#---
# Default values
#---
path.FastLmm   <- 'fastlmmc'
debug          <- FALSE
exitAfterError <- !debug		# Exit after any error, (unless we are in debug mode)
tmpDir         <- '.'

#---
# Parse  command line arguments
#---
if( !exists('cmdLineArgs') )    { cmdLineArgs <- commandArgs(trailingOnly = TRUE); }

# Stop if there are no command line argument
if( length(cmdLineArgs) < 1 ) { fatalError('No command line arguments!\n'); }

dataFileStr     <- cmdLineArgs[1];
tfamFile        <- cmdLineArgs[2];
kinshipFile     <- cmdLineArgs[3];
subBlockSize    <- as.integer( cmdLineArgs[4] );
onlyOnce        <- (cmdLineArgs[5] == 'TRUE') || (cmdLineArgs[5] == 'T')
debug           <- debug || onlyOnce;		# Set debug mode

cat("ASKAT arguments:\n");
cat("\tData file/s         : ", dataFileStr , "\n" );
cat("\tTFAM file           : ", tfamFile , "\n" );
cat("\tKinship matrix file : ", kinshipFile , "\n" );
cat("\tSub-block size      : ", subBlockSize , "\n" );
cat("\tTemporal dir        : ", tmpDir , "\n" );

#---
# Should we exit immediatly after any error or warning?
#---
if( exitAfterError ) {
	options(
		#warn=2,	# This will change all warnings into errors, so warnings will also be handled like errors
        error= quote({ traceback(); q("no",status=1,FALSE) })	# Make sure we get an error code after a 'stop'
	);
}

#---
# Read files
#---
if( debug )	{ cat('Loading tfam file: ', tfamFile , '\n' ); }
tfam <- read.csv(tfamFile, sep="", header=FALSE, col.names=c('familyId','individualId', 'paternalId', 'maternalId', 'sex', 'phenotype') );

if( debug )	{ cat('Loading kinship & FaST-LMM file: ', kinshipFile , '\n' ); }
load(kinshipFile);


# More than one file (comma separated list of files)
dataFiles <- unlist( strsplit(dataFileStr , ",") )

for( dataFile in dataFiles ) {
	cat("Data file : ", dataFile , "\n" );

	# Read data file
	if( debug )	{ cat('Load data file: ', dataFile, '\n' ); }
	dat  <- read.csv(dataFile, sep="", header=FALSE );


	#---
	# Iterate for each sub-block
	#---
	snpIdx <- 5:dim(dat)[2];							# Columns having SNP data
	ped12 <- tfam[,c(2,6)]; 							# First two columns of data structure (inividualID and phenotype)
	sbIdx <- seq( 1, dim(dat)[1] , by = subBlockSize );	# SubBlock indeces

	# Iterate on every sub-block
	for( i in sbIdx )  {
		# Create pedigree matrix for ASKAT function
		maxBlock <- min( dim(dat)[1] , i+subBlockSize-1 );
		snpsBlock <- i:maxBlock;
		ped <- cbind( ped12, t(dat[snpsBlock,snpIdx]) );

		if( debug )	{ cat('Iterating on sub-block: ', paste( dat[i,1], ':', dat[i,4], ' - ' , dat[maxBlock,1], ':', dat[maxBlock,4], sep="") , '\n' ); }

		# Call ASKAT
		results <- ASKAT(ped, fastlmm)

		# Show results
		cat("\nASKAT_RESUTS:"
			, "p-value:" , results$pvalue.ASKAT 
			, "chr:pos:", paste( dat[i,1], ':', dat[i,4], ' - ' , dat[maxBlock,1], ':', dat[maxBlock,4], sep="")
			, "Block:", dataFile
			, "Sub-Block:", paste( i, ' - ', maxBlock, sep="" )
			, "Id:", paste( dat[i,2], ' - ', dat[maxBlock,2], sep="" )
			, "Q:", results$Q.ASKAT 
			, "Polygenic.VC:", results$Polygenic.VC
			, "Env.VC:", results$Env.VC
			, "lambda:", results$lambda
			, "\n"
			, sep="\t" 
			);

		if( onlyOnce )	{ 
			# Excecute only one sub-block? => Stop now
			fatalError('Execute onlyOnce is set. Stopping after first iteration.\n'); 
		}
	}
}

