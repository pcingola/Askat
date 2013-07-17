pheno.geno.kin.Without.NA <-
function(Ped.All, kin2){
pheno.geno0 = cbind(Ped.All, c(1:dim(Ped.All)[1]))
cont1 = 0
for (i in 2:dim(Ped.All)[2]){
                            pheno.geno1 = pheno.geno0
                            pheno.geno0 = pheno.geno0[which(!is.na(pheno.geno0[,(i-cont1)])),(1:dim(pheno.geno0)[2])]
                            if (dim(pheno.geno0)[1] < 300){
                                                            pheno.geno0 = pheno.geno1[-(i-cont1)]
                                                            cont1 = cont1 + 1
                                                            }
                            }
pheno.geno = pheno.geno0[,1:(dim(Ped.All)[2] - cont1)]

kin.FaST = 2 * kin2

ord.kin2.No.Miss = as.integer(pheno.geno0[,dim(pheno.geno0)[2]])
kin2.data.miss0 = kin2[ord.kin2.No.Miss,ord.kin2.No.Miss]

n.col = Ped.All[ord.kin2.No.Miss,1]
kin.FaST.data.miss = kin.FaST[ord.kin2.No.Miss,ord.kin2.No.Miss]
kin.FaST.data.miss1 = cbind(n.col,kin.FaST.data.miss)
kin.FaST.data.miss1 = rbind(c("var",n.col),kin.FaST.data.miss1)
file.name.data = paste("kin.FaST.data.miss.txt",sep="")
write(t(kin.FaST.data.miss1[ ,1:(dim(kin.FaST.data.miss1)[1])]), file = file.name.data, ncolumns = dim(kin.FaST.data.miss1)[1], sep= "\t")

res = list(kin2.data.miss = kin2.data.miss0, pheno.Geno.No.NA = pheno.geno, ord.kin2.No.Miss = ord.kin2.No.Miss, W = kin.FaST)
}

