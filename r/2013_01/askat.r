#-------------------------------------------------------------------------------
# ASKAT Package
#
# By: Karim Oualkacha
# Modifications by: Pablo Cingolani
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
# Ped    : Is the Pedigree data file. It has subject IDs as first column (IDs 
#          should be differents for all subjects), phenotype as second column 
#          and region-based SNPs that will be analized together  
#
# tfam   : This is the TFAM matrix (PLINK's TFAM format file)
#
# kin1   : The kinship matrix. If the kinship matrix noted is calculated from 
#          GenABEL we should convert it to the kinship matrix using this 
#          command: 
#              kin1 = diagReplace(kin1, upper=TRUE)
#
# Missing: A logical parameter, if TRUE, it means that there is missing values 
#          in the Pedigree data set: ASKAT function moves out subjects with 
#          missing data
#
#-------------------------------------------------------------------------------
ASKAT <- function(Ped, tfam, kin1, simFileName, phenoFileName, tfamFile ) {

	#####  STEP 1: Check for Missing Data and construction of Pedigree without Missing Data   #####

	# Directory for eigenvalue output (Fast-LMM)
	eigenDir <- tmpFile("ASKAT_FaSTLMM");

	#
	Y.trait = Ped[,2]
	X = as.matrix(Ped[,3:dim(Ped)[2]])

	##### STEP 2: Under the NULL we call FaST-LMM to estimate the VCs ####
	res.VC.FaST = VC.FaST.LMM(Ped, tfam, eigenDir, simFileName, phenoFileName, tfamFile)
	Estim.Sigma.RG = res.VC.FaST$Polygenic
	Estim.Sigma.e = res.VC.FaST$Env
	S = res.VC.FaST$S
	U = res.VC.FaST$U

	##### STEP 3: Calculation of weights matrix W and the matrix K =GWG #####
	freq.MAF = apply(X, 2, mean)/2

	if( length(freq.MAF) == 1){
		w = (dbeta(freq.MAF, 1, 25))^2
		K = w * X %*% t(X)
	} else {
		w = vector(length = length(freq.MAF))
		for (i in 1:length(freq.MAF)){
			w[i] = (dbeta(freq.MAF[i], 1, 25))^2
		}
		w = diag(w)
		K = X %*% w %*% t(X)
	}

	##### STEP 4: ASKAT score test statistic calculations #####
	Gamma = Estim.Sigma.RG / Estim.Sigma.e
	D.0 = (Gamma * S)  + diag(1, dim(X)[1], dim(X)[1])
	inv.sqrt.D.0 = diag(1/sqrt(diag(D.0)))

	K.tilde = inv.sqrt.D.0 %*% t(U)

	un.n = c(rep(1,dim(U)[1]))
	X.tilde = K.tilde %*% un.n
	Y.tilde = K.tilde %*% Y.trait

	K.tilde =  K.tilde %*% K %*% t(K.tilde)

	P.0.tilde = diag(1, dim(U)[1], dim(U)[2]) - ( X.tilde %*% solve( t(X.tilde) %*% X.tilde ) %*% t(X.tilde) )
	res = P.0.tilde %*% Y.tilde
	s2 = Estim.Sigma.e

	Q = t(res) %*% K.tilde
	Q = Q %*% res/(2 * s2)
	W1 = P.0.tilde %*% K.tilde
	W1 = W1 %*% P.0.tilde/2
	out = Get_PValue.Modif(W1, Q)
	pvalue.davies = out$p.value
	lambda = out$lambda

	return( list(pvalue.ASKAT = pvalue.davies, Q.ASKAT = Q, Polygenic.VC = Estim.Sigma.RG, Env.VC = Estim.Sigma.e, lambda = lambda, tmp.kinshipFileName <- simFileName) );
}

#-------------------------------------------------------------------------------
# Write a TPED file for FastLmm
# NOTE: Since FaST-LMM doesn't actually use this file, we can use any random data
#-------------------------------------------------------------------------------
Geno.FaST.LMM <- function(geno, tpedFile){

#	# Create a new vector
#	geno.FaST = vector(length = 2 * length(geno))
#
#	# Code SNPs as hom/het using 'T' and 'A' (4 and 1)
#	# So we code {'T/T', 'A/T', 'A/A'} as {'4/4', '1/4', '1/1'}
#	for (i in 1:length(geno)){
#
#		cat("DEBUG CODE!!!!\n");
#		geno[i] = floor( 3 * runif(200) )
#
#		if (geno[i] == 0){
#			geno.FaST[((2*(i-1))+1)] = 4
#			geno.FaST[(2*i)] = 4
#		}
#
#		if (geno[i] == 1){
#			geno.FaST[((2*(i-1))+1)] = 1
#			geno.FaST[(2*i)] = 4
#		}
#
#		if (geno[i] == 2){
#			geno.FaST[((2*(i-1))+1)] = 1
#			geno.FaST[(2*i)] = 1
#		}
#	}

	# File already exists: Don't even bother.
	if( file.exists(tpedFile) )	{ 
		cat('File ', tpedFile ,' already exists. Nothing done.\n');
		return; 
	}

	geno.FaST = floor( 2 * runif(2 * length(geno)) ) + 1;	# A random vector of '1' and '2' (representing 'A' and 'C'). Note: Any pther combination should also do the trick.
	geno.test = c(1, "snp0", 0, 0, 1, geno.FaST)
	geno.test = rbind(geno.test)

	if( debug )	{ cat('Writing tped matrix to file: ', tpedFile , '\tsize: ', dim(geno.test), '\n' ); }
	cat('GENO.TEST:', geno.test , '\n');
	write.matrix(geno.test, file = tpedFile, sep="")
}

#-------------------------------------------------------------------------------
# Get lambda
#-------------------------------------------------------------------------------
Get_Lambda <- function (K) {
	out.s <- eigen(K, symmetric = TRUE)
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
# Invoke FastLmm and get VC
#
# Note: In order to invoke Fast-LMM we must create some temporal files
#       After invoking the program, we have to read and parse the result files
#-------------------------------------------------------------------------------
VC.FaST.LMM <- function(Ped, tfam, eigenDir, kinshipFileName, phenoFileName, tfamFile ){
	# TMP File names to be used
	genoName <- tmpFile("geno_test", sep="_");				# FaSTLMM version 2.03 doesn't allow any '/' in this name (bug)
	genoTfamFileName <- paste(genoName,".tfam", sep="");
	genoTpedFileName <- paste(genoName,".tped", sep="");
	genoOutFileName <- paste(genoName,".out.txt", sep="");
	fastlmmOutFileName <- tmpFile("OUTFaST-LMM.txt");

	#---
	# Create TMP files used by the program
	#---

	# Create TPED file
	Y.trait = Ped[,2]
	FID = tfam[,1]
	IID = tfam[,2]
	Geno = Ped[,3]
	Geno.FaST.LMM(Geno, genoTpedFileName)

	# Copy TFAM file
	if( ! file.exists(genoTfamFileName) ) {
		file.copy(tfamFile, genoTfamFileName);
	}

	#---
	# Invoke Fast-LMM
	#---

	# Create Fast-LMM command line and execute it
	fastlmmcCmd <- paste( path.FastLmm		# Full path to fastlmm binary
		, "-tfile", genoName				# basename for PLINK's transposed .tfam and .tped files
		, "-sim" , kinshipFileName			# file containing the genetic similarity matrix
		, "-eigenOut", eigenDir				# save the spectral decomposition object to the directoryname
		, "-pheno", phenoFileName			# name of phenotype file
		, "-out", genoOutFileName			# name of output file
		, "-maxThreads", 1					# Number of threads to use (by the math library). We perform paralelization, so we need to set this to one
		, "-mpheno 1"						# index for phenotype in -pheno file to process, starting at 1 for the first phenotype column
		);
	if( debug )	{ cat('Execute system command: ', fastlmmcCmd , '\n' ); }
	retCode <- system(fastlmmcCmd)
	if( retCode != 0 )	fatalError( paste("Cannot execute command\n\t", fastlmmcCmd) );

	#---
	# Read Fast-LMM results (from TMP files)
	#---

	# Read results from 'genoOutFileName'
	if( debug )	{ cat('Reading VC.FaST.LMM table from: ', genoOutFileName , '\n' ); }
	VC.FaST.LMM = read.table(genoOutFileName, header=TRUE)
	Estim.Sigma.RG = VC.FaST.LMM$NullGeneticVar
	Estim.Sigma.e = VC.FaST.LMM$NullResidualVar
	cat('\tEstim.Sigma.RG  : ', Estim.Sigma.RG, '\n');
	cat('\tEstim.Sigma.e   : ', Estim.Sigma.e, '\n');

	# Read 'S' matrix
	p = dim(Ped)[1]
	sFileName <- paste(eigenDir,"S.bin", sep="/")
	if( debug )	{ cat('Reading read.SVD.bin S-matrix from: ', sFileName , '\n' ); }
	read.SVD.bin = file(sFileName, "rb")
	S = readBin(read.SVD.bin, "numeric", n=p, endian="little")
	close(read.SVD.bin)
	S = diag(sort(S, decreasing = T))

	# Read 'U' matrix
	uFileName <- paste(eigenDir,"U.bin", sep="/")
	if( debug )	{ cat('Reading upper read.SVD.bin U-matrix from: ', uFileName , '\n' ); }
	read.SVD.bin = file(uFileName, "rb")
	U = readBin(read.SVD.bin, "numeric", n=p*p, endian="little")
	close(read.SVD.bin)

	U = matrix(U,p,p, byrow=F)
	U = U[ ,ncol(U):1]

	# Remove TMP files and dirs
	# WARNING: A function should NOT have side effects (e.g. deleting a file created somewhere else)
	if( !debug) unlink( c(fastlmmOutFileName, genoOutFileName ) );
	if( !debug) unlink( eigenDir, recursive = TRUE)
	# unlink( kinshipFileName );	# May be we can re-use the kinship matrix file in the next iteration

	# Create results list
	return( list(Polygenic = Estim.Sigma.RG, Env = Estim.Sigma.e, S = S, U = U) );
}

#-------------------------------------------------------------------------------
# Create a name for a temporal file
#-------------------------------------------------------------------------------
tmpFile <- function(name, sep="/")	{ return( paste(tmpDir, name, sep=sep) ); }

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
debug          <- TRUE
exitAfterError <- !debug		# Exit after any error, (unless we are in debug mode)
tmpDir         <- '.'

#---
# Parse  command line arguments
#---
if( !exists('cmdLineArgs') )    { cmdLineArgs <- commandArgs(trailingOnly = TRUE); }

# Stop if there are no command line argument
if( length(cmdLineArgs) < 1 ) { fatalError('No command line arguments!\n'); }

dataFile        <- cmdLineArgs[1];
tfamFile        <- cmdLineArgs[2];
kinshipFile     <- cmdLineArgs[3];
subBlockSize    <- as.integer( cmdLineArgs[4] );
tmpDir          <- cmdLineArgs[5];
path.FastLmm    <- cmdLineArgs[6];
onlyOnce        <- (cmdLineArgs[7] == 'TRUE') || (cmdLineArgs[7] == 'T')
simFileName     <- cmdLineArgs[8];
phenoFileName   <- cmdLineArgs[9];
debug           <- debug || onlyOnce;		# Set debug mode

cat("ASKAT arguments:\n");
cat("\tData file           : ", dataFile , "\n" );
cat("\tTFAM file           : ", tfamFile , "\n" );
cat("\tKinship matrix file : ", kinshipFile , "\n" );
cat("\tSub-block size      : ", subBlockSize , "\n" );
cat("\tTemporal dir        : ", tmpDir , "\n" );
cat("\tFast-LMM path       : ", path.FastLmm , "\n" );
cat("\tOnly once           : ", onlyOnce , "\n" );
cat("\tSimilarity file     : ", simFileName , "\n" );
cat("\tPhenotype file      : ", phenoFileName , "\n" );

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
if( debug )	{ cat('Load data file: ', dataFile, '\n' ); }
dat  <- read.csv(dataFile, sep="", header=FALSE );

if( debug )	{ cat('Loading tfam file: ', tfamFile , '\n' ); }
tfam <- read.csv(tfamFile, sep="", header=FALSE, col.names=c('familyId','individualId', 'paternalId', 'maternalId', 'sex', 'phenotype') );

if( debug )	{ cat('Loading kinship file: ', kinshipFile , '\n' ); }
load(kinshipFile);

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
	results <- ASKAT(ped, tfam, kinshipMatrix, simFileName, phenoFileName, tfamFile)

	# Show results
	cat("\nASKAT_RESUTS:"
		, "p-value:" , results$pvalue.ASKAT 
		, "Block:", dataFile
		, "Sub-Block:", paste( i, ' - ', maxBlock, sep="" )
		, "chr:pos:", paste( dat[i,1], ':', dat[i,4], ' - ' , dat[maxBlock,1], ':', dat[maxBlock,4], sep="")
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

# Remove tmp kinship matrix file (used for fastlmm)
if( !debug) unlink( results$tmp.kinshipFileName );

