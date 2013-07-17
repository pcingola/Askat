Geno.FaST.LMM <-
function(Geno){
Geno.FaST = vector(length = 2*length(Geno))
for (i in 1:length(Geno)){
if (Geno[i] == 0){
                 Geno.FaST[((2*(i-1))+1)] = 4
                 Geno.FaST[(2*i)] = 4
                 }
if (Geno[i] == 1){
                 Geno.FaST[((2*(i-1))+1)] = 1
                 Geno.FaST[(2*i)] = 4
                 }
if (Geno[i] == 2){
                 Geno.FaST[((2*(i-1))+1)] = 1
                 Geno.FaST[(2*i)] = 1
                 }
}
geno.test = c(1, "snp0", 0, 0, 1, Geno.FaST)
geno.test = rbind(geno.test)
write.matrix(geno.test, file = paste("geno_test.tped", sep=""), sep="")
}

