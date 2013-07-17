library(CompQuadForm)
library(nFactors)

###############################
path1 = '.';

source(file = paste(path1, "/ASKAT.R", sep=""))
source(file = paste(path1, "/Geno.FaST.LMM.R", sep=""))
source(file = paste(path1, "/Get_Lambda.R", sep=""))
source(file = paste(path1, "/Get_PValue.Modif.R", sep=""))
source(file = paste(path1, "/pheno.geno.kin.Without.NA.R", sep=""))
source(file = paste(path1, "/VC.FaST.LMM.R", sep=""))

load(paste(path1, "/kinsip4000SNPs.RData", sep=""))

Ped.pheno.Geno = read.table(paste(path1, "/Ped4ASKAT.dat", sep=""), header=FALSE)
Y.trait = Ped.pheno.Geno[,2]
Missing = FALSE
GenoAll = Ped.pheno.Geno[,3:dim(Ped.pheno.Geno)[2]]

IDs = c(1:length(Y.trait))
wp1 = seq(0, dim(GenoAll)[2], by = 10)

p.value.ASKAT = vector(length = (length(wp1)-1))

#for (i in 1:(length(wp1)-1)){
for (i in c(1, 101, 201, 301, 401) ) {
	Ped = cbind(IDs, Y.trait, GenoAll[,(wp1[i]+1):wp1[i+1]])
	resultats = ASKAT(Ped, kin1, Missing)
	p.value.ASKAT[i] = resultats$pvalue.ASKAT
	cat('Index: ', i, '\t[', wp1[i]+1, ',' , wp1[i+1] , ']\t' , 'pvalue: ', p.value.ASKAT[i] ,'\n')
}

p.value.ASKAT
write.table(p.value.ASKAT, file = "pvalueASKAT1.dat", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

log.pval.Obs1 = sort(- log10(p.value.ASKAT), decreasing = T)
p= length(log.pval.Obs1)
log.pval.Exp = - log10(c(1:p)/(p+1))

pdf("QQplotSimDATA_ASKAT1.pdf")
plot(log.pval.Exp, log.pval.Obs1, xlim = c(0,13), ylim = c(0,13), ylab = expression("-log10(Obs p-value)"), xlab = "-log10(Exp p-value)", main= paste("ASKAT, ", p, " Windows, 10 SNPs in each W", sep=""))
abline(0,1)
dev.off()

