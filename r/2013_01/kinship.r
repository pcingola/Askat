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
createKinship <- function(kinshipFile, tfam ) {
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
	convert.snp.tped(tped = tpedFile, tfam = tfamFile, out = genabelGenFile, strand = "+")

	#---
	# Calculate and save kinship matrix
	#---

	# Load GenABEL data
	data.GenABEL <- load.gwaa.data(phe = genabelPhenFile, gen = genabelGenFile, force = T)

	# Next command calculates the kinship matrix using all autosomal markers that we have in the tped file
	cat('Calculating Kinship matrix (IBS)\n');
	kinshipMatrix = ibs(data.GenABEL[, autosomal(data.GenABEL)], weight = "freq")
	kinshipMatrix = diagReplace(kinshipMatrix, upper=TRUE)

	# Save kinship as RData file
	cat('Saving file', kinshipFile, '\n');
	save( kinshipMatrix, file=kinshipFile );
	kinshipMatrix;
}

#-------------------------------------------------------------------------------
# Create 'pheno.txt' file (used by FaST-LMM)
#-------------------------------------------------------------------------------
createPheno <- function(phenoFileName, tfam){
	# Create 'pheno' file
	pheno.test = data.frame(tfam$familyId, tfam$individualId, tfam$phenotype)
	colnames(pheno.test) = c( 'familyId', 'individualId', 'phenotype');
	cat('Writing phenotype to file: ', phenoFileName , '\n' );
	write.table(pheno.test, file = phenoFileName, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
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
createSim <- function(simFileName, tfam, kinshipMatrix ) {

	# Create 'simmilarity' file for FaST-LMM (simFileName).
	kin.FaST = 2 * kinshipMatrix
	n.col = paste(tfam[,1] , tfam[,2])		# Create labels 'family_ID individual_ID'
	colnames(kin.FaST) <- n.col
	rownames(kin.FaST) <- n.col

	# Write simmilarity file
	cat('Writing simmilarity matrix to file: ', simFileName , '\n' );
	write.table(c('var\t'), file = simFileName, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, eol = "");	# This is just to create the required 'var'
	write.table(kin.FaST, file = simFileName, quote = FALSE, sep = "\t", row.names = TRUE , col.names = TRUE , append=TRUE); # Now we dump the data
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

debug          <- FALSE
exitAfterError <- !debug        # Exit after any error, (unless we are in debug mode)

if( exitAfterError ) {
    options(
        #warn=2,	# This will change all warnings into errors, so warnings will also be handled like errors
        error= quote({ traceback(); q("no",status=1,FALSE) })	# Make sure we get an error code after a 'stop'
    );
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
phenoFileName	<- cmdLineArgs[7];

cat("Kinship arguments:\n");
cat("\tInput TPED file                 :", tpedFile, "\n");
cat("\tInput TFAM file                 :", tfamFile, "\n");
cat("\tOutput GenABEL 'genotype' file  :", genabelGenFile, "\n");
cat("\tOutput GenABEL 'phenotype' file :", genabelPhenFile, "\n");
cat("\tOutput Kinship matrix file      :", kinshipFile, "\n");
cat("\tOutput SIM matrix file          :", simFile, "\n");
cat("\tOutput pheno file               :", phenoFileName, "\n");

#---
# Read TFAM
#---
cat('Reading TFAM file\n');
tfam <- read.csv(tfamFile, sep="", header=FALSE, col.names=c('familyId','individualId', 'paternalId', 'maternalId', 'sex', 'phenotype') );

# Create (or load) kinship file
if( file.exists(kinshipFile) ) {
	cat("Kinship file '", kinshipFile ,"' exists. Loading.\n");
	load( kinshipFile )
} else {
	kinshipMatrix <- createKinship(kinshipFile, tfam );
}

# Create SIM file (for FaST-LMM)
createSim(simFile, tfam, kinshipMatrix );

# Create pheno.txt file (for FaST-LMM)
createPheno(phenoFileName, tfam);
