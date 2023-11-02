#simpheno
library(Matrix)
#install.packages("PhenotypeSimulator")
#Biocmanager::install("snpStats")
library(PhenotypeSimulator)

#read in and clean real data

#convert compound genotype into 0,1,2
#missing data encoded as 0
convert_genotype <- function(x){
  missing <- x==0
  #get counts of alleles from vector of compound genotypes
  alleles <- table(unlist(strsplit(x, split="")))
  alleles <- alleles[order(alleles)]
  alleles <- alleles[names(alleles) != 0]
  A1 <- names(alleles)[1]
  A2 <- names(alleles)[2]
  minor_homo <- x==paste0(A1, A1)
  major_homo <- x==paste0(A2, A2)
  
  x <- rep(1, length(x))
  x[minor_homo] <- 2
  x[major_homo] <- 0
  x[missing]    <- 0
  x
}

#read in real data.
#will use the real genotypes with simulated phenotypes
real_data <- read.csv("../mohan_real_data/Merged_Information_Final.csv")
real_data <- real_data[real_data$centre == 1, ]

#Remove QC SNPs
myvars<-names(real_data) %in% c("rs63750847","rs3213787","rs9621532","rs4151667","rs9332739")
real_data<-real_data[!myvars]

#create geno matrix
geno_matrix <- sapply(real_data[,c(17:50)], convert_genotype)
colinear = sapply(1:ncol(geno_matrix), function(x) rankMatrix(geno_matrix[,-x])[[1]])
colinear = which(colinear==rankMatrix(geno_matrix)[[1]])[1]
geno_matrix <- geno_matrix[, -colinear]
rownames(geno_matrix) <- real_data$pid
save(geno_matrix, file="geno_matrix.Rdata")
# correlation<-matrix(NA, 7,7)
# correlation[1,]<-c(1,0.47,0.31,0.41,0.4,0.42,0.3)
# correlation[2,]<-c(0.47,1,0.27,0.42,0.48,0.49,0.3)
# correlation[3,]<-c(0.31,0.27,1,0.51,0.29,0.24,0.3)
# correlation[4,]<-c(0.41,0.42,0.51,1,0.46,0.42,0.3)
# correlation[5,]<-c(0.40,0.48,0.29,0.46,1,0.88,0.3)
# correlation[6,]<-c(0.42,0.49,0.24,0.42,0.88,1,0.3)
# correlation[7,]<-c(0.3,0.3,0.3,0.3,0.3,0.3,1)
# 
# 
# sim_pheno <- function(SNPvar, causalSNPs, p_trait, p_ind_gen, p_trait_ind_gen){
#   shared   = 1 - p_ind_gen; N = nrow(geno_matrix); P = 7
#   X_causal <- geno_matrix[, causalSNPs]
#   corrBG   <- correlatedBgEffects(N = N , P = P , corr_mat = correlation)
#   genFixed <- geneticFixedEffects(N = N , P = P, X_causal = X_causal,
#                                   pTraitsAffected = p_trait,
#                                   pIndependentGenetic = p_ind_gen,
#                                   pTraitIndependentGenetic = p_trait_ind_gen)
#   
#   genFixed_s <- rescaleVariance(genFixed$shared, shared * SNPvar)$component
#   corrBG_    <- rescaleVariance(corrBG$correlatedBg, 1-SNPvar)$component
#   
#   if (p_ind_gen == 0){
#     pheno <- scale(genFixed_s + corrBG_)
#   } else {
#     genFixed_i <- rescaleVariance(genFixed$independent, (1-shared)*SNPvar)$component
#     pheno <- scale(genFixed_s + genFixed_i + corrBG_)
#   }
#   
#   #last column is the binary eye disease, use threshold to turn binary
#   pheno[,P] <- ifelse(pheno[,P] > 0, 1, 0)
#   return(pheno)
# }
# 
# 
# #--------------------------
# #----------------simulations
# #---------------------------
# SNPvar <- seq(0.01,0.05,0.01)
# Nsnps  <- ncol(geno_matrix)
# n_ind  <- nrow(geno_matrix)
# Nreps  <- 400
# P      <- 7
# 
# ####################
# #    scenario 1 ####
# ####################
# 
# #all 7 phenotypes associated w genetic effect
# # 3 causal SNPs randomly selected
# # no independent SNP effects (all shared)
# # SNP var 0.01-0.05
# 
# n_causal    <- 3
# n_notcausal <- Nsnps - n_causal
# 
# pheno_array <- array(NA, dim=c(n_ind, P, Nreps, 5))
# causal_array<- array(NA, dim=c(n_causal, Nreps, 5))
#                      
# for (i in 1:5){
#   print(i)
#   for (j in 1:Nreps){
#     causalSNPs   = sample(1:ncol(geno_matrix), n_causal)
#     causalSNPs_  = colnames(geno_matrix)[causalSNPs]
#     SNPvar_      = SNPvar[i]
#     pheno_matrix = sim_pheno(SNPvar = SNPvar_, causalSNPs = causalSNPs_, p_trait = 1,
#                              p_ind_gen = 0, p_trait_ind_gen = 0)
#     pheno_array[,,j,i] = pheno_matrix
#     causal_array[,j,i] = causalSNPs
#   }
# }
# save(pheno_array, causal_array, file="scenario1_phenos.Rdata")
# ####################
# #    scenario 2 ####
# ####################
# 
# #all 7 phenotypes associated w genetic effect
# # 15 causal SNPs randomly selected
# # no independent SNP effects (all shared)
# # SNP var 0.01-0.05
# 
# n_causal    <- 15
# pheno_array <- array(NA, dim=c(n_ind, P, Nreps, 5))
# causal_array<- array(NA, dim=c(n_causal, Nreps, 5))
#                      
# for (i in 1:5){
#   for (j in 1:Nreps){
#      causalSNPs   = sample(1:ncol(geno_matrix), n_causal)
#      causalSNPs_  = colnames(geno_matrix)[causalSNPs]
#      SNPvar_      = SNPvar[i]
#      pheno_matrix = sim_pheno(SNPvar = SNPvar_, causalSNPs = causalSNPs_, p_trait = 1,
#                               p_ind_gen = 0, p_trait_ind_gen = 0)
#      pheno_array[,,j,i] = pheno_matrix
#      causal_array[,j,i] = causalSNPs
#   }
# }                   
# save(pheno_array, causal_array, file="scenario2_phenos.Rdata")
# 
# ####################
# #    scenario 3 ####
# ####################
# 
# #4/7 phenotypes associated w genetic effect
# # 3 causal SNPs randomly selected
# # no independent SNP effects (all shared)
# # SNP var 0.01-0.05
# 
# n_causal    <- 3
# pheno_array <- array(NA, dim=c(n_ind, P, Nreps, 5))
# causal_array<- array(NA, dim=c(n_causal, Nreps, 5))
#                      
# for (i in 1:5){
#  for (j in 1:Nreps){
#    causalSNPs   = sample(1:ncol(geno_matrix), n_causal)
#    causalSNPs_  = colnames(geno_matrix)[causalSNPs]
#    SNPvar_      = SNPvar[i]
#    pheno_matrix = sim_pheno(SNPvar = SNPvar_, causalSNPs = causalSNPs_, p_trait = 0.5,
#                             p_ind_gen = 0, p_trait_ind_gen = 0)
#    pheno_array[,,j,i] = pheno_matrix
#    causal_array[,j,i] = causalSNPs
#  }
# }                   
# save(pheno_array, causal_array, file="scenario3_phenos.Rdata")
# 
# ####################
# #    scenario 4 ####
# ####################
# 
# # all 7 phenotypes associated w genetic effect
# # 5 causal SNPs randomly selected
# #2 of the phenos are only affected by independent effects
# #genetic effect partitioned into 0.4 independent, 0.6 shared
# # SNP var 0.01-0.05
# 
# n_causal    <- 5
# pheno_array <- array(NA, dim=c(n_ind, P, Nreps, 5))
# causal_array<- array(NA, dim=c(n_causal, Nreps, 5))
#                      
# for (i in 1:5){
#  for (j in 1:Nreps){
#    causalSNPs   = sample(1:ncol(geno_matrix), n_causal)
#    causalSNPs_  = colnames(geno_matrix)[causalSNPs]
#    SNPvar_      = SNPvar[i]
#    pheno_matrix = sim_pheno(SNPvar = SNPvar_, causalSNPs = causalSNPs_, p_trait = 1,
#                             p_ind_gen = 0.4, p_trait_ind_gen = 0.2)
#    pheno_array[,,j,i] = pheno_matrix
#    causal_array[,j,i] = causalSNPs
#  }
# }                   
# save(pheno_array, causal_array, file="scenario4_phenos.Rdata")
# 
# ####################
# #    scenario 5 ####
# ####################
# 
# # all 7 phenotypes associated w genetic effect
# # 15 causal SNPs randomly selected
# #2 of the phenos are only affected by independent effects
# #genetic effect partitioned into 0.4 independent, 0.6 shared
# # SNP var 0.01-0.05
# 
# n_causal    <- 15
# pheno_array <- array(NA, dim=c(n_ind, P, Nreps, 5))
# causal_array<- array(NA, dim=c(n_causal, Nreps, 5))
#                      
#  for (i in 1:5){
#    for (j in 1:Nreps){
#      causalSNPs   = sample(1:ncol(geno_matrix), n_causal)
#      causalSNPs_  = colnames(geno_matrix)[causalSNPs]
#      SNPvar_      = SNPvar[i]
#      pheno_matrix = sim_pheno(SNPvar = SNPvar_, causalSNPs = causalSNPs_, p_trait = 1,
#                               p_ind_gen = 0.4, p_trait_ind_gen = 0.2)
#      pheno_array[,,j,i] = pheno_matrix
#      causal_array[,j,i] = causalSNPs
#    }
#  }                   
#  save(pheno_array, causal_array, file="scenario5_phenos.Rdata")
#                      
# ####################
# #    scenario 6 ####
# ####################
# 
# #all 7 phenotypes associated w genetic effect
# # 25 causal SNPs randomly selected
# # no independent SNP effects (all shared)
# # SNP var 0.01-0.05
# 
# n_causal    <- 25
# pheno_array <- array(NA, dim=c(n_ind, P, Nreps, 5))
# causal_array<- array(NA, dim=c(n_causal, Nreps, 5))
#                     
# for (i in 1:5){
#   for (j in 1:Nreps){
#     causalSNPs   = sample(1:ncol(geno_matrix), n_causal)
#     causalSNPs_  = colnames(geno_matrix)[causalSNPs]
#     SNPvar_      = SNPvar[i]
#     pheno_matrix = sim_pheno(SNPvar = SNPvar_, causalSNPs = causalSNPs_, p_trait = 1,
#                              p_ind_gen = 0, p_trait_ind_gen = 0)
#     pheno_array[,,j,i] = pheno_matrix
#     causal_array[,j,i] = causalSNPs
#   }
# }                   
# save(pheno_array, causal_array, file="scenario6_phenos.Rdata")                     
