# remove.packages(GENLIB)
# install.packages()

library(GENLIB)
qqq <- read.csv("EyeDis227.csv")
qqq <- gen.genealogy(qqq)
founders <- gen.founder(qqq)
chr = 1

for (i in 1:chr){
  # load ("chr_params.Rdata")
  
  BP_len   = 290000000
  nmarker  = 100
  pos      = 1:nmarker * (BP_len %% nmarker)
  # load("founders.Rdata")
  # founder  = 
  nfounder = length(founders)
  founders = rep(founders, each=2)
  founders = founders * c(1, -1)
  nfounder = length(founders)
  
  MAF  <- runif(nmarker, 0.15, 0.45)
  MAF_ <- rep(MAF, each=2)
  MAF_ <- floor(MAF_ * nfounder)
  
  geno_matrix <- matrix(-9, nrow=nfounder, ncol= nmarker*2 + 6)
  geno_matrix[, 7:ncol(geno_matrix)] <- sapply(MAF_, function(x) sample(c(rep(1, x), rep(0, nfounder - x)), nfounder)) 
  geno_matrix[,1] <- founders
  geno_matrix[,2] <- founders
  geno_matrix[,3] <- 0
  geno_matrix[,4] <- 0
  geno_matrix[,5] <- 1
  geno_matrix[,6] <- 1
  write.table(geno_matrix, file="./founder_geno/chr1.ped", quote=F, row.names=F, col.names=F, sep='\t')
  
  mapfile <- matrix(nrow=nmarker, ncol=4)
  mapfile[,1] <- 1
  mapfile[,2] <- 1:nmarker
  mapfile[,3] <- 0
  mapfile[,4] <- pos
  write.table(mapfile, file="./founder_geno/chr1.map", quote=F, row.names=F, col.names=F, sep='\t')
}

gen.drop(qqq, model_params=c(3.28, 1.96), cM_len=c(328, 196), BP_len=BP_len, mapfile_path = "./founder_geno/chr1.map", 
         pedfile_path = "./founder_geno/chr1.ped", out = "./proband_geno/chr1")

for (i in 1:chr){
  gen.drop = function (gen, pro=NULL, ancestors=NULL, model =1, model_params, cM_len, 
                       BP_len, nsimul = 100, physical_map_Mo = NULL, physical_map_Fa = NULL, 
                       mapfile_path, pedfile_path, out = NULL, seed=0)
}
#combine into one array and save