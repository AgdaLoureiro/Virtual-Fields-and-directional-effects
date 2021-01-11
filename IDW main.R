rm(list = ls())
source("idws.R")

anis_forte =  c(30, 0.3)
anis_fraco =  c(30, 0.9)
anis = matrix(data = NA, nrow = 3, ncol = 2)
anis[1,] = anis_forte
anis[2,] = anis_fraco

beta_fraco = c(2,0.001,0.001)
beta_forte = c(2,0.002,0.002)

beta = matrix(data = NA, nrow = 2, ncol = 3)
beta[1,] = beta_forte
beta[2,] = beta_fraco

for (i in 1:nrow(beta)) {
  for (j in 1:nrow(anis)) {
    idws(trend = beta[i,], anispar = anis[j,])
  }
}


