load ("chr_params.Rdata")

BP_len   = 290000000
nmarker  = 100
pos      = 1:nmarker * (BP_len %% nmarker)
load("founders.Rdata")
# founder  = 
nfounder = length(founder)
founder  = rep(founder, each=2)
founder  = founder * c(1, -1)
nfounder = length(founders)

MAF  <- runif(nmarker, 0.15, 0.45)
MAF_ <- rep(MAF, each=2)
MAF_ <- floor(MAF_ * nfounder * 2)

geno_matrix <- matrix(-9, nrow=nfounder, ncol= nmarker + 6)
geno_matrix[, 7:ncol(geno_matrix)] <- sapply(MAF_, function(x) sample(c(rep(1, x), rep(0, nfounder*2 - x)), nfounder*2)) 
geno_matrix[,1] <- founders
geno_matrix[,2] <- founders
geno_matrix[,3] <- 0
geno_matrix[,4] <- 0
geno_matrix[,5] <- 1
geno_matrix[,6] <- 1
write.table(mapfile, "./founder_geno/chr1.map", quote=F, row.names=F, col.names=F, sep='\t')

mapfile <- matrix(nrow=nmarker, ncol=4)
mapfile[,1] <- 1
mapfile[,2] <- 1:nmarker
mapfile[,3] <- 0
mapfile[,4] <- pos
write.table(mapfile, "./founder_geno/chr1.map", quote=F, row.names=F, col.names=F, sep='\t')