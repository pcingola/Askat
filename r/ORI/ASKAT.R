ASKAT <-
function(Ped, kin1, Missing){
# Missing is a logical parameter, if TRUE, it means that there is missing values in the Pedigree data set: ASKAT function moves out subjects with missing data
# kin1 is the kinship matrix. If the kinship matrix noted is calculated from GenABEL we should convert it to the kinship matrix using this command: kin1 = diagReplace(kin1, upper=TRUE)
# Ped is the Pedigree data file: it has subject IDs as first column (IDs should be differents for all subjects), phenotype as second column and region-based SNPs that will be analized together  

#####  STEP 1: Check for Missing Data and construction of Pedigree without Missing Data   #####
##### Also creation of the file "kin.FaST.data.miss.txt" which will be needed by FaST-LMM #####
if (Missing == "TRUE"){
                      data.Without.NA = pheno.geno.kin.Without.NA(Ped, kin1)
                      Ped = data.Without.NA$pheno.Geno.No.NA
                      }
if (Missing == "FALSE"){
                       kin.FaST = 2 * kin1
                       n.col = Ped[,1]
                       kin.FaST = cbind(n.col,kin.FaST)
                       kin.FaST = rbind(c("var",n.col),kin.FaST)
                       file.name.data = paste("kin.FaST.data.miss.txt",sep="")
                       write(t(kin.FaST[ ,1:(dim(kin.FaST)[1])]), file = file.name.data, ncolumns = dim(kin.FaST)[1], sep= "\t")
                       }
#####
p = dim(Ped)[1]
Y.trait = Ped[,2]
X = as.matrix(Ped[,3:dim(Ped)[2]])

##### STEP 2: Under the NULL we call FaST-LMM to estimate the VCs ####
res.VC.FaST = VC.FaST.LMM(Ped)
Estim.Sigma.RG = res.VC.FaST$Polygenic
Estim.Sigma.e = res.VC.FaST$Env
pvalue.FaST.LMM = res.VC.FaST$pvalue.FaST.LMM

read.SVD.bin = file(paste("OURSKAT_FaST/S.bin",sep=""), "rb")
S = readBin(read.SVD.bin, "numeric", n=p, endian="little")
close(read.SVD.bin)
S = diag(sort(S, decreasing = T))
read.SVD.bin = file(paste("OURSKAT_FaST/U.bin",sep=""), "rb")
U = readBin(read.SVD.bin, "numeric", n=p*p, endian="little")
close(read.SVD.bin)
U = matrix(U,p,p, byrow=F)
U = U[ ,ncol(U):1]
system("rm -r OURSKAT_FaST")

##### STEP 3: Calculation of weights matrix W and the matrix K =GWG #####
freq.MAF = apply(X, 2, mean)/2

if( length(freq.MAF) == 1){
                          w = (dbeta(freq.MAF, 1, 25))^2
                          K = w * X %*% t(X)
                          } else
                                {
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
out = Get_PValue.Modif(W1, Q)    ### This function was taken from SKAT package
pvalue.davies = out$p.value
lambda = out$lambda

resultats = list(pvalue.ASKAT = pvalue.davies, Q.ASKAT = Q, Polygenic.VC = Estim.Sigma.RG, Env.VC = Estim.Sigma.e, lambda = lambda)
resultats
}

