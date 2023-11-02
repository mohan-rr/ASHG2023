#run association on simulated phenotype matrices
library(MDMR)
library(CompQuadForm)
source("./GAMuT-functions.R")
source("RA_functions.R")

print("check1")
# load data
#############
# args = commandArgs(trailingOnly = TRUE)
job_num = sample(1:120, 1)
print("job num: ")
print(job_num)
# job_num  = as.integer(args[1])
#job num will be 1:120
job_num   = job_num - 1
scenario  = (job_num %/% 20) + 1 # {1,2,3,4,5,6}
SNPvar    = (job_num %% 20 %/% 4) + 1 #{1,2,3,4,5}
reps      = job_num %% 20 %% 4 #{0,1,2,3}
print(reps)
rep_range = (reps*100 + 1):((reps + 1)*100)
Nreps = 100

load("geno_matrix.Rdata")
n_snps  = ncol(geno_matrix)

pheno_file = paste0("scenario",scenario,"_phenos.Rdata")
load(pheno_file)
#has 2 objects:
# pheno_array  : dim (n_ind, P, 400, 5)
# causal_array : dim (n_causal, 400, 5)
print(rep_range)
print(SNPvar)
pheno_array = pheno_array[,,rep_range,SNPvar] # now has dim (n_ind, P, Nreps)
causal_array= causal_array[, rep_range, SNPvar] # dim (n_causal, Nreps)
n_causal    = dim(causal_array)[1]
n_notcausal = n_snps - n_causal
print("check2")

#-------------
#functions to run assoc tests
MDMR_pval <- function(pheno_matrix, geno_matrix){
  D <- dist(pheno_matrix, method="euclidean")
  p <- mdmr(X=geno_matrix, D=D)
  p <- p$pv$analytic.pvals
  p
}

#GAMUT
gamut_pval <- function(pheno_matrix, geno_matrix) {
  pheno_ <- as.matrix(apply(pheno_matrix, 2, scale, scale=T, center=T))
  linear_pheno <- linear_GAMuT_pheno(pheno_)
  Yc = linear_pheno$Kc           # linear kernel similarity matrix
  lambda_Y = linear_pheno$ev_Kc  # eigenvalues of Yc
  # sample Major allele frequency  of each causal variant in the sample
  MAF    = colMeans(geno_matrix)/2                       
  betas  = rep(1/length(MAF),length(MAF))
  w_geno = geno_matrix %*% diag(betas)      # Weighted common causal variants
  c_geno = scale(w_geno, center=T, scale=F) # centered genotype matrix
  #construct the weighted linear kernel matrix for genotypes
  linear_geno <- linear_GAMuT_geno(c_geno)
  X_geno      <- linear_geno$Lc    # linear kernel similarity matrix
  lambda_X    <- linear_geno$ev_Lc  
  
  p <- TestGAMuT(Yc,lambda_Y,X_geno ,lambda_X)
  return(p)
}

#MANOVA
manova_pval <- function(pheno_matrix, geno_matrix) {
  p <- vector("numeric", ncol(geno_matrix))
  for (j in 1:ncol(geno_matrix)){
    manova_ <- manova(pheno_matrix~factor(geno_matrix[,j]))
    p[j] <- as.numeric(summary(manova_)$stats[,6][1])
  } 
  return(p)
}

reverse_regression <- function(geno, pheno)  # one SNP for each test
{
  nn = nrow(pheno)
  P  = rowsum(rep(1,nn), geno)/nn
  DF = ncol(pheno)
  yy = rowsum(pheno, geno)
  w  = as.vector(c(1-P[1], P[3]-P[1], -(1-P[3])) %*% yy)
  S = as.vector((w %*% solve(var(pheno), w)))/nn/prod(1-P)
  p.value = 1-pchisq(S, DF)
  return(p.value)
}
#-------------

#-------------run assoc tests and store pvals, power, type1
power = list(gamut  = array(NA, dim=c(Nreps)),
             MDMR_o = array(NA, dim=c(Nreps)),
             MDMR   = array(NA, dim=c(n_causal, Nreps)),
             revreg = array(NA, dim=c(n_causal, Nreps)),
             manova = array(NA, dim=c(n_causal, Nreps)),
             revAB  = array(NA, dim=c(n_causal, Nreps)))

type1 = list(MDMR   = array(NA, dim=c(n_notcausal, Nreps)),
             revreg = array(NA, dim=c(n_notcausal, Nreps)),
             manova = array(NA, dim=c(n_notcausal, Nreps)),
             revAB  = array(NA, dim=c(n_notcausal, Nreps)))
print("check3")
for (j in 1:Nreps){
  pheno_matrix = as.matrix(pheno_array[,,j])
  causalSNPs   = as.vector(causal_array[,j])
  p_mdmr    = MDMR_pval(pheno_matrix, geno_matrix)
  p_mdmr_o  = p_mdmr[1]
  p_mdmr    = p_mdmr[2:length(p_mdmr)]
  p_rev     = apply(X=geno_matrix, MARGIN=2, FUN=reverse_regression, pheno=pheno_matrix)
  p_manova  = manova_pval(pheno_matrix, geno_matrix)
  p_RAB     = apply(X=geno_matrix, MARGIN=2, FUN=RA_assoc_indep, y=pheno_matrix, z=NULL, HWE=T)
  
  if (j==1) print("check4")
  power[["gamut"]][j]  = gamut_pval(pheno_matrix, geno_matrix)
  power[["MDMR_o"]][j] = p_mdmr_o
  power[["MDMR"]][,j]  = p_mdmr[causalSNPs]
  power[["revreg"]][,j]= p_rev[causalSNPs]
  power[["manova"]][,j]= p_manova[causalSNPs]
  power[["revAB"]][,j] = p_RAB[causalSNPs]
  
  type1[["MDMR"]][,j]   = p_mdmr[-causalSNPs]
  type1[["revreg"]][,j] = p_rev[-causalSNPs]
  type1[["manova"]][,j] = p_manova[-causalSNPs]
  type1[["revAB"]] [,j] = p_RAB[-causalSNPs]
  if (j==1) print("check5")
}

print("check6")
alpha = 0.05/n_snps

power_vec = c(mean(power[["gamut"]] < alpha),
              mean(power[["MDMR_o"]]< alpha),
              mean(power[["MDMR"]]  < alpha),
              mean(power[["revreg"]]< alpha),
              mean(power[["manova"]]< alpha),
              mean(power[["revAB"]] < alpha))
names(power_vec) = c("gamut", "MDMR_o", "MDMR", "revreg", "manova", "revAB")

type1_vec = c(mean(type1[["MDMR"]] < alpha),
              mean(type1[["revreg"]]< alpha),
              mean(type1[["manova"]]< alpha),
              mean(type1[["revAB"]] < alpha))
names(type1_vec)= c("MDMR", "revreg", "manova", "revAB")

filename1 = paste0("power_scenario", scenario,"_snpvar", SNPvar, "_",reps+1,".Rdata") 
filename2 = paste0("type1_scenario", scenario,"_snpvar", SNPvar, "_",reps+1,".Rdata")

save(power, power_vec, file=filename1)
save(type1, type1_vec, file=filename2)
print("\ndone")





