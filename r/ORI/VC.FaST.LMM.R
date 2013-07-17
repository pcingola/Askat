VC.FaST.LMM <-
function(Ped){
Y.trait = Ped[,2]
IID = Ped[,1]
Geno = Ped[,3]
FID = c(rep(1,length(IID)))
Geno.FaST.LMM(Geno)

pheno.test = cbind(FID, IID, Y.trait)
write.table(pheno.test, file = paste("pheno.txt", sep=""), sep = "     ", quote = FALSE, row.names = FALSE, col.names = TRUE)
geno.tfam = cbind(pheno.test[,1], pheno.test[,2], c(rep(0,length(Y.trait))), c(rep(0,length(Y.trait))), c(rep(0,length(Y.trait))), Y.trait)
write.table(geno.tfam, file = paste("geno_test.tfam", sep=""), sep = "     ", quote = FALSE, row.names = FALSE, col.names = FALSE)
system(paste("fastlmmc -tfile geno_test -sim kin.FaST.data.miss.txt -eigenOut OURSKAT_FaST -pheno pheno.txt -mpheno 1", sep=""))
output = scan(paste("geno_test.out.txt", sep=""), what="character", sep="\n")
write.table(output[2], file = paste("OUTFaST-LMM.txt", sep=""), quote = FALSE, sep = " ", eol = "\r\n", na = "NA", row.names = FALSE, col.names = FALSE)
VC.FaST.LMM = read.table(paste("OUTFaST-LMM.txt", sep=""), header=FALSE)
Estim.Sigma.RG = VC.FaST.LMM[1,19]
Estim.Sigma.e = VC.FaST.LMM[1,20]
pvalue.FaST.LMM = VC.FaST.LMM[1,9]
system("rm OUTFaST-LMM.txt")
system("rm geno_test.out.txt")
system("rm pheno.txt")
system("rm geno_test.tfam")
system("rm geno_test.tped")
system("rm kin.FaST.data.miss.txt")

ans = list(Polygenic = Estim.Sigma.RG, Env = Estim.Sigma.e, pvalue.FaST.LMM = pvalue.FaST.LMM)
ans
}

