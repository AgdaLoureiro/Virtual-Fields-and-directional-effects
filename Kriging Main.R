rm(list = ls())
source("interpo_2.R")


nug_a1 = 0.2
nug_b1 = 0.16
nug_c1 = 0.42
nug_d1 = 0.17
psil_a1 = 2.2 
psil_b1 = 0.1
psil_c1 = 1
psil_d1 = 0.08  
range_a1 = 900
range_b1 = 400
range_c1 = 900
range_d1 = 400

chute_inicial_dense = c(nug_a1, nug_b1, nug_c1, nug_d1,
                  psil_a1, psil_b1, psil_c1, psil_d1,
                  range_a1, range_b1, range_c1, range_d1)

nug_a2 = 0.2
nug_b2 = 0.16
nug_c2 = 0.6
nug_d2 = 0.17
psil_a2 = 4 
psil_b2 = 2
psil_c2 = 5
psil_d2 = 6  
range_a2 = 600 
range_b2 = 200
range_c2 = 1000
range_d2 = 400

chute_inicial_sparse = c(nug_a2, nug_b2, nug_c2, nug_d2,
                         psil_a2, psil_b2, psil_c2, psil_d2,
                         range_a2, range_b2, range_c2, range_d2)
anis_forte =  c(30, 0.3)
anis_fraco =  c(30, 0.9)
anis = matrix(data = NA, nrow = 2, ncol = 2)
anis[1,] = anis_forte
anis[2,] = anis_fraco

beta_fraco = c(2,0.001,0.001)
beta_forte = c(2,0.002,0.002)

beta = matrix(data = NA, nrow = 2, ncol = 3)
beta[1,] = beta_forte
beta[2,] = beta_fraco

for (i in 1:nrow(beta)) {
  for (j in 1:nrow(anis)) {
    interpo(trend = beta[i,], anispar = anis[j,], chute_ini_1 = chute_inicial_dense,
            chute_ini_2 = chute_inicial_sparse)
  }
}


