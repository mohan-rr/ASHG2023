library(Rcpp)

sourceCpp("founder_genotypes.cpp")
#founder_genos(sample pool size, founder vector, num markers, MAF for all markers, lambda, out)
founder_genos(1000, as.integer(1:20), 10000, 0.30, 20, "test")