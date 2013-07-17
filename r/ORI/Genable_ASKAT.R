#######################################################################################
#                  KINSHIP matrix calculation using GenABEL package                   #
#######################################################################################
library(GenABEL)
library(nFactors)
library(CompQuadForm)

# Here we convert PLINK tped and tfam files to GenABEL data and we calculate the the kinship matrix using "ibs" function 
# We need tped and tfam files
# tped should have subject alleles coded by letters, EX: A,C,G,T (see page 233 of GenABEL tutorial)

# path = pathway for ".tped" and ".tfam" files
path = "C:/Users/karim.oualkacha/Desktop/STT66943_Genetique/Rare_Variants/ASKAT_Rpackage/askat_test/askat_test/"

# Next command converts PLINK tped and tfam files to GenABEL data
convert.snp.tped(tped = paste(path,"geno_cov.tped",sep=""), tfam = paste(path, "geno_cov.tfam",sep=""), out = "gen0tped.raw", strand = "+")

# To load the data in the next step, we need the phenotype file: In this example it called "GenableTrait.dat", and it contains: subject IDs, SEX and the quantitaive phenotype. This file should have headers for each column: id   sex   pheno 
data.GenABEL <- load.gwaa.data(phe = "GenableTrait.dat", gen = "gen0tped.raw", force = T)

phenodata = phdata(data.GenABEL)  
matrix.Genotypes = as.numeric(gtdata(data.GenABEL))

# Next command calculates the kinship matrix using all autosomal markers that we have in the tped file
kin1 = ibs(data.GenABEL[, autosomal(data.GenABEL)], weight = "freq")
kin1 = diagReplace(kin1, upper=TRUE)

# Next command gives the Pedigree file as is used in ASKAT function (see "ASKAT.R" file)
Ped = cbind(phenodata[,c(1,3)],matrix.Genotypes)

# Next step: Now you can run ASKAT

#path = "your pathway to the R files below"
source(file = paste(path, "/ASKAT.R", sep=""))
source(file = paste(path, "/Geno.FaST.LMM.R", sep=""))
source(file = paste(path, "/Get_Lambda.R", sep=""))
source(file = paste(path, "/Get_PValue.Modif.R", sep=""))
source(file = paste(path, "/pheno.geno.kin.Without.NA.R", sep=""))
source(file = paste(path, "/VC.FaST.LMM.R", sep=""))

Missing = FALSE
ASKAT(Ped, kin1, Missing)


