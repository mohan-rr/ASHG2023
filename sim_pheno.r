#simpheno
#install.packages("PhenotypeSimulator")
#Biocmanager::install("snpStats")
library(PhenotypeSimulator)

# load("../geno_matrix.Rdata")
# load("../genedrop_array.Rdata")
# dim geno_array: 227, 33, 400
load("../ped_kin_mx.Rdata")
# "kin_mx"

correlation<-matrix(NA, 7,7)
correlation[1,]<-c(1,0.47,0.31,0.41,0.4,0.42,0.3)
correlation[2,]<-c(0.47,1,0.27,0.42,0.48,0.49,0.3)
correlation[3,]<-c(0.31,0.27,1,0.51,0.29,0.24,0.3)
correlation[4,]<-c(0.41,0.42,0.51,1,0.46,0.42,0.3)
correlation[5,]<-c(0.40,0.48,0.29,0.46,1,0.88,0.3)
correlation[6,]<-c(0.42,0.49,0.24,0.42,0.88,1,0.3)
correlation[7,]<-c(0.3,0.3,0.3,0.3,0.3,0.3,1)

sim_pheno <- function(kinvar){
  N = 227; P = 7
  corrBG   <- correlatedBgEffects(N = N , P = P , corr_mat = correlation)
  genBG    <- geneticBgEffects(P, N, kin_mx, shared=TRUE)
  
  genBG  <- rescaleVariance(genBG$shared, kinvar)$component
  corrBG <- rescaleVariance(corrBG$correlatedBg, 1-kinvar)$component
  pheno  <- scale(genBG + corrBG)
  #last column is the binary eye disease, use threshold to turn binary
  pheno[,P] <- ifelse(pheno[,P] > 0, 1, 0)
  return(pheno)
}


# #--------------------------
# #----------------simulations
# #---------------------------
kinvar <- seq(0.01,0.10,0.01)
n_ind  <- 227
Nreps  <- 400
P      <- 7

pheno_array <- array(NA, dim=c(n_ind, P, Nreps, length(kinvar)))
                     
for (j in 1:Nreps){
    print(j)
    for (i in 1:length(kinvar)){
        print(i)
        pheno_matrix = sim_pheno(kinvar[i])
        pheno_array[,,j,i] = pheno_matrix
  }
}
save(pheno_array, file="pheno_array.Rdata")
