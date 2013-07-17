#-------------------------------------------------------------------------------
# KINSHIP matrix calculation using GenABEL package
#
# Here we convert PLINK tped and tfam files to GenABEL data and we calculate 
# kinship matrix using "ibs" function. We need tped and tfam files
#
# Note: tped should have subject alleles coded by letters, EX: A,C,G,T (see 
#       page 233 of GenABEL tutorial)
#
# Author     : Karim Oualkacha
# Adapted by : Pablo Cinoglani
#
#																			2012
#-------------------------------------------------------------------------------

library(GenABEL)
library(nFactors)
library(CompQuadForm)

#-------------------------------------------------------------------------------
# Create kinship file
#-------------------------------------------------------------------------------
createKinship <- function( tfam , tpedFile ) {
	# Convert TFAM to GenABEL phenotype
	# Save as genable phenotype file
	# It contains: subject IDs, SEX and the quantitaive phenotype. 
	# This file should have headers for each column: "id   sex   pheno"
	cat('Writing GenABEL phenotype file\n');
	sexOri <- (2 - tfam$sex);			# Note: Sex has to be coded as 0=female and 1=male, instead of 1=male, 2=female)
	sex <- pmin(1, pmax(0,sexOri) );	# Make sure it is always {0,1}
	if( sum(sex!=sexOri) )	{ cat("\n\nWARNING: Sex had to be adjusted to {0,1}. This will probably create some problems\n\n"); }
	genPhen <- data.frame( id=tfam$individualId, sex=sex, pheno=tfam$phenotype );
	write.table( genPhen, file=genabelPhenFile, col.names=TRUE, row.names=FALSE, quote=FALSE);

	#---
	# Converts PLINK tped file to GenABEL genotype
	#---
	cat('Converting TPED to GenABEL file:', genabelGenFile,'\n');
	convert.snp.tped(tped = tpedFile, tfam = tfamFile, out = genabelGenFile, strand = "+")

	#---
	# Calculate and save kinship matrix
	#---

	# Load GenABEL data
	cat('Loafind GenABEL file: ', genabelGenFile,'\n');
	data.GenABEL <- load.gwaa.data(phe = genabelPhenFile, gen = genabelGenFile, force = T)

	# Next command calculates the kinship matrix using all autosomal markers that we have in the tped file
	cat('Calculating Kinship matrix (IBS)\n');
	kinshipMatrix = ibs(data.GenABEL[, autosomal(data.GenABEL)], weight = "freq")
	cat('Calculating diagReplace on Kinship matrix\n');
	kinshipMatrix = diagReplace(kinshipMatrix, upper=TRUE)

	kinshipMatrix;
}

#-------------------------------------------------------------------------------
# Create 'pheno.txt' file (used by FaST-LMM)
#-------------------------------------------------------------------------------
createPheno <- function(phenoFile, tfam){
	# Create 'pheno' file
	pheno.test = data.frame(tfam$familyId, tfam$individualId, tfam$phenotype)
	colnames(pheno.test) = c( 'familyId', 'individualId', 'phenotype');
	cat('Writing phenotype to file: ', phenoFile , '\n' );
	write.table(pheno.test, file = phenoFile, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

#-------------------------------------------------------------------------------
# Create 'sim' file for FaST-LMM
#
# From FastLMM's handbook:
#		Instead of SNP data from which genetic similarities are computed, the user may provide
#		the genetic similarities directly using the –sim <filename> option. The file containing
#		the genetic similarities should be tab delimited and have both row and column labels for
#		the IDs (family ID and individual ID separated by a space). The value in the top-left
#		corner of the file should be var.
#
# Arguments:
#
# tfam            : This is the TFAM matrix (PLINK's TFAM format file)
#
# kinshipMatrix   : The kinship matrix. If the kinship matrix noted is calculated from 
#          GenABEL we should convert it to the kinship matrix using this 
#          command: 
#              kin1 = diagReplace(kin1, upper=TRUE)
#
#-------------------------------------------------------------------------------
createSim <- function(simFile, tfam, kinshipMatrix ) {
	cat('Creating simmilarity\n' );

	# Create 'simmilarity' file for FaST-LMM (simFile).
	kin.FaST = 2 * kinshipMatrix
	n.col = paste(tfam[,1] , tfam[,2])		# Create labels 'family_ID individual_ID'
	colnames(kin.FaST) <- n.col
	rownames(kin.FaST) <- n.col

	# Write simmilarity file
	cat('Writing simmilarity matrix to file: ', simFile , '\n' );
	write.table(c('var\t'), file = simFile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, eol = "");	# This is just to create the required 'var'
	write.table(kin.FaST, file = simFile, quote = FALSE, sep = "\t", row.names = TRUE , col.names = TRUE , append=TRUE); # Now we dump the data
}

#-------------------------------------------------------------------------------
# Create a TPED file for FastLmm
# NOTE: Since FaST-LMM doesn't actually use this file, we can use any random data
#-------------------------------------------------------------------------------
createTpedFastlmm <- function(numIndividuals, tpedFile){
	cat('Creating TPED file (for FaST-LMM)\n' );

	# File already exists: Don't even bother.
	if( file.exists(tpedFile) )	{ 
		if( debug) { cat('File ', tpedFile ,' already exists. Nothing done.\n'); }
		return; 
	}

	geno.FaST = floor( 2 * runif(2 * numIndividuals) ) + 1;	# A random vector of '1' and '2' (representing 'A' and 'C'). Note: Any pther combination should also do the trick.
	geno.test = c(1, "snp0", 0, 0, 1, geno.FaST)
	geno.test = rbind(geno.test)

	if( debug )	{ cat('Writing tped matrix to file: ', tpedFile , '\tsize: ', dim(geno.test), '\n' ); }
	write.matrix(geno.test, file = tpedFile, sep="")
}

#-------------------------------------------------------------------------------
# Fatal error
#-------------------------------------------------------------------------------
fatalError <- function(errStr)	{ 
	errStr <- paste('FATAL ERROR:',  errStr, "\n");
	stop( errStr ); 
}

#-------------------------------------------------------------------------------
# Invoke FastLmm and get VC
#
# Under the NULL we call FaST-LMM to estimate the VCs 
#
# Note: In order to invoke Fast-LMM we must create some temporal files
#       After invoking the program, we have to read and parse the result files
#-------------------------------------------------------------------------------
invokeFastlmm <- function(tfam, simFile, phenoFile, tfamFile) {
	cat('Invoking FaST-LMM\n' );

	# Directory for eigenvalue output (Fast-LMM)
	eigenDir <- tmpFile("ASKAT_FaSTLMM");

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
	numIndividuals = dim(tfam)[1]
	createTpedFastlmm(numIndividuals, genoTpedFileName)

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
		, "-sim" , simFile					# file containing the genetic similarity matrix
		, "-eigenOut", eigenDir				# save the spectral decomposition object to the directoryname
		, "-pheno", phenoFile				# name of phenotype file
		, "-out", genoOutFileName			# name of output file
		, "-mpheno 1"						# index for phenotype in -pheno file to process, starting at 1 for the first phenotype column
		);
	if( debug )	{ cat('Execute system command: ', fastlmmcCmd , '\n' ); }
	retCode <- system(fastlmmcCmd)
	if( retCode != 0 )	fatalError( paste("Cannot execute command\n\t", fastlmmcCmd) );

	#---
	# Read Fast-LMM results (from TMP files)
	#---

	# Read results from 'genoOutFileName'
	if( debug )	{ cat('Reading FaST-LMM table from: ', genoOutFileName , '\n' ); }
	res.fastlmm = read.table(genoOutFileName, header=TRUE)
	nullGeneticVar = res.fastlmm$NullGeneticVar
	nullResidualVar = res.fastlmm$NullResidualVar

	# Read 'S' matrix
	p = dim(tfam)[1]
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
	if( !debug) unlink( c(fastlmmOutFileName, genoTfamFileName, genoTpedFileName, genoOutFileName, fastlmmOutFileName ) );
	if( !debug) unlink( eigenDir, recursive = TRUE)

	# Create results list
	return( list(nullGeneticVar = nullGeneticVar, nullResidualVar = nullResidualVar, S = S, U = U) );
}

#-------------------------------------------------------------------------------
# Create a name for a temporal file
#-------------------------------------------------------------------------------
tmpFile <- function(name, sep="/")	{ return( paste(tmpDir, name, sep=sep) ); }

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

debug          <- FALSE
exitAfterError <- !debug        # Exit after any error, (unless we are in debug mode)

if( exitAfterError ) {
	# Make sure we get an error code after a 'stop' (and we quit the process)
    options( error= quote({ traceback(); q("no",status=1,FALSE) })	);
} else {
	# Make sure we get a stack trace
    options( error= traceback );
}

#---
# Parse command line arguments
# Debug: 
#	cmdLineArgs <- c( 'geno_cov.tped', 'geno_cov.tfam', 'geno_cov.genable.gen', 'geno_cov.genable.phen', 'geno_cov.kinship.RData', 'geno_cov.sim', 'geno_cov.pheno.txt' );
#---
if( !exists('cmdLineArgs') )	{ cmdLineArgs <- commandArgs(trailingOnly = TRUE); }

tpedFile		<- cmdLineArgs[1];
tfamFile		<- cmdLineArgs[2];
genabelGenFile	<- cmdLineArgs[3];
genabelPhenFile	<- cmdLineArgs[4];
kinshipFile		<- cmdLineArgs[5];
simFile			<- cmdLineArgs[6];
phenoFile	    <- cmdLineArgs[7];
path.FastLmm    <- cmdLineArgs[8];

cat("Kinship arguments:\n");
cat("\tInput TPED file                 :", tpedFile, "\n");
cat("\tInput TFAM file                 :", tfamFile, "\n");
cat("\tOutput GenABEL 'genotype' file  :", genabelGenFile, "\n");
cat("\tOutput GenABEL 'phenotype' file :", genabelPhenFile, "\n");
cat("\tOutput Kinship matrix file      :", kinshipFile, "\n");
cat("\tOutput SIM matrix file          :", simFile, "\n");
cat("\tOutput pheno file               :", phenoFile, "\n");
cat("\tFast-LMM path                   :", path.FastLmm , "\n" );

#---
# TMP dir (form tpedFile)
#---
len <- nchar(tpedFile);
extLen <- nchar( ".tped" );
if( len <= 5 )	{ fatalError(paste("TPED file", tpedFile, "does not have '.tped' extension\n")); }
ext <- tolower( substr( tpedFile, len - extLen + 1 , len ) );
if( ext != ".tped" )	{ fatalError(paste("TPED file", tpedFile, "does not have '.tped' extension\n")); }
tmpDir <- substr( tpedFile, 0 , len - extLen );
dir.create(tmpDir, showWarnings = FALSE);

#---
# Read TFAM
#---
cat('Reading TFAM file\n');
tfam <- read.csv(tfamFile, sep="", header=FALSE, col.names=c('familyId','individualId', 'paternalId', 'maternalId', 'sex', 'phenotype') );

#---
# Create (or load) kinship file
#---
if( file.exists(kinshipFile) ) {
	cat("Kinship file '", kinshipFile ,"' exists. Loading.\n");
	load( kinshipFile )
} else {
	kinshipMatrix <- createKinship( tfam, tpedFile );
}

#---
# Create SIM file (for FaST-LMM)
#---
createSim(simFile, tfam, kinshipMatrix );

#---
# Create pheno.txt file (for FaST-LMM)
#---
createPheno(phenoFile, tfam);

#---
# Invoke FaST-LMM
#---

fastlmm = invokeFastlmm(tfam, simFile, phenoFile, tfamFile );

#---
# Save kinship as RData file
#---
cat('Saving results to file', kinshipFile, '\n');
save( kinshipMatrix, fastlmm, file=kinshipFile );

if( !debug) unlink( tmpDir, recursive = TRUE)
