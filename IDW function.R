  idws = function(trend,anispar){ 
  #We will use the following packages for this first part
  library(pacman)
  pacman::p_load(gstat, sp, dplyr,raster, lattice, ggplot2, RStoolbox, dplyr)
  
  
  beta = trend
  anis = anispar
  
  #We will first generate x and y 
  
  # virtual data
  x.range <- c(5.5, 1400.5)
  y.range <- c(5.5, 1450.5)
  area <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 10), 
                      y = seq(from = y.range[1], to = y.range[2], by = 10))
  
  nsim = 100
  #We firt define a model of semivariogram and its parameters
  set.seed(93)
  model <- gstat(formula=z ~ x+y, 
                 locations = ~ x+y,
                 dummy=T,#Inconditional simulation term 
                 beta= beta, 
                 model=vgm(nugget =0, psill= 0.28,
                           range= 400, model='Sph', 
                           anis = anis),
                 nmin = 4, #min number of nearest observations
                 nmax=100,
                 weights = T,#max nummber of nearest observations
                 data = area)
  
  #Now we perform the Simulation
  t1 = Sys.time()
  set.seed(1993)
  #number of simulations that will be done is seted in this part of the code
  #simulation uses the model made in the previous step
  vf <- predict(model, newdata=area, nsim= nsim)
  t2 = Sys.time()
  print (t2-t1)
  
  
  #In a computer with 8 cores and 8GB of RAM, it takes around 21 minutes
  
  #####################Sampling###################################
  
  
  #Now we sample from the VF
  
  #For this task, I used the following packages
  pacman::p_load(Rcpp, hydroGOF, raster, rgdal, sp)
  
  #Now we set x and y locations
  #I used two sample densities, so:
  
  #first sampling (dense with 100 units of spacing)
  Coords_am1 <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 100), 
                            y = seq(from = y.range[1], to = y.range[2], by = 100))
  
  #second sampling (sparse with 200 units of spacing)
  Coords_am2 <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 200), 
                            y = seq(from = y.range[1], to = y.range[2], by = 200))
  
  #third, we eliminate coords_am1 and coords_am2 from the VF
  val = vf
  val1 <- dplyr::anti_join(vf, Coords_am1)
  val_tot <- dplyr::anti_join(val1, Coords_am2) 
  
  #And then we extract the values of each VF
  #For dense sampling
  dense_sampling = dplyr::inner_join(vf, Coords_am1)
  
  #For sparse sampling
  sparse_sampling = dplyr::inner_join(vf, Coords_am2)
  
  #For this purpose I draw a boundary of the VF in QGIS:
  boundary <- shapefile("./data/boundary/cotorno.shp")
  r = raster(boundary, res = 10) 
  rp = rasterize(boundary, r, 0) 
  grid <- as(rp, "SpatialPixelsDataFrame")
  
  mapas_idw_dense = list()
  mapas_idw_sparse = list()
  
  idw1 = list()
  idw2 = list()
  for(i in 3:length(vf)) {
    
    #dense
    data1 <- dense_sampling[, c(1,2,i)]
    names(data1) <- c("x", "y", "z")
    
    sp::coordinates(data1) = ~x+y
    proj4string(data1) <- CRS(proj4string(grid))
    
    #Cross Validation
    
  
    eval.idw1 = data.frame(expand.grid (k = c(seq (0.5, 4, by = 0.5)),
                                       me_train = NA,
                                       mae_train = NA,
                                       rmse_train = NA,
                                       ave_train = NA,
                                       r2_train= NA))
    set.seed(42)
    for (j in 1:nrow(eval.idw1)){
      model1 = gstat::krige.cv(z  ~ 1, locations = data1, nfold = nrow(data1),set = list(idp = j))
      
      eval.idw1$me_train[j] = hydroGOF::me.data.frame(model1$var1.pred, model1$observed)
      eval.idw1$mae_train[j] = hydroGOF::mae(model1$var1.pred, model1$observed)
      eval.idw1$rmse_train[j] = hydroGOF::rmse(model1$var1.pred, model1$observed)
      eval.idw1$r2_train[j] = hydroGOF::br2(model1$var1.pred, model1$observed)
      eval.idw1$ave_train[j] = hydroGOF::NSE(model1$var1.pred, model1$observed)
  
      }  
    
    library(dplyr)
    
    eval.idw1 = eval.idw1 %>% dplyr::select(k, me_train, mae_train, rmse_train, ave_train) %>% group_by(k) %>%
      summarise(me_train = mean(me_train),
                mae_train = mean(mae_train),
                rmse_train = mean(rmse_train),
                ave_train = mean(ave_train)) %>%
      arrange(desc(ave_train))
    
    idw_exp1 = eval.idw1 %>%  filter(ave_train == max(ave_train)) %>%  dplyr::select(k) 
    
    model1 = gstat::krige.cv(z  ~ 1, locations = data1, nfold = nrow(data1), set = list(idp = idw_exp1$k))
    
    ##Accessing Metrics
    idw1[[i]] = data.frame(mae = hydroGOF::mae(model1$var1.pred, model1$observed),
                          me = hydroGOF::me.data.frame(model1$var1.pred, model1$observed),
                          rmse = hydroGOF::rmse(model1$var1.pred, model1$observed),
                          r2 = hydroGOF::br2(model1$var1.pred, model1$observed),
                          ave = hydroGOF::NSE(model1$var1.pred, model1$observed),
                          willmott = hydroGOF::md(model1$var1.pred, model1$observed)) 
  
    map1 <- idw(z~1, data1, grid, idp = idw_exp1$k)
  
    
    mapa_idw1 <- raster(map1)
    proj4string(mapa_idw1) <- CRS("+init=epsg:32722")
    
    mapas_idw_dense[[i]] = mapa_idw1
    
    
    
    ####sparse
    data2 <- sparse_sampling[, c(1,2,i)]
    names(data2) <- c("x", "y", "z")
    
    sp::coordinates(data2) = ~x+y
    proj4string(data2) <- CRS(proj4string(grid))
    
    #Cross Validation
    
    
    eval.idw2 = data.frame(expand.grid (k = c(seq (0.5, 4, by = 0.5)),
                                        me_train = NA,
                                        mae_train = NA,
                                        rmse_train = NA,
                                        ave_train = NA,
                                        r2_train= NA))
    set.seed(42)
    for (j in 1:nrow(eval.idw2)){
      model2 = gstat::krige.cv(z  ~ 1, locations = data2, nfold = nrow(data2),set = list(idp = j))
      
      eval.idw2$me_train[j] = hydroGOF::me.data.frame(model2$var1.pred, model2$observed)
      eval.idw2$mae_train[j] = hydroGOF::mae(model2$var1.pred, model2$observed)
      eval.idw2$rmse_train[j] = hydroGOF::rmse(model2$var1.pred, model2$observed)
      eval.idw2$r2_train[j] = hydroGOF::br2(model2$var1.pred, model2$observed)
      eval.idw2$ave_train[j] = hydroGOF::NSE(model2$var1.pred, model2$observed)
   
      
    }  
    
    library(dplyr)
    
    eval.idw2 = eval.idw2 %>% dplyr::select(k, me_train, mae_train, rmse_train, ave_train) %>% group_by(k) %>%
      summarise(me_train = mean(me_train),
                mae_train = mean(mae_train),
                rmse_train = mean(rmse_train),
                ave_train = mean(ave_train)) %>%
      arrange(desc(ave_train))
    
    idw_exp2 = eval.idw2 %>%  filter(ave_train == max(ave_train)) %>%  dplyr::select(k) 
    
    model2 = gstat::krige.cv(z  ~ 1, locations = data2, nfold = nrow(data2), set = list(idp = idw_exp2$k))
    
    ##Accessing Metrics
    idw2[[i]] = data.frame(mae = hydroGOF::mae(model2$var1.pred, model2$observed),
                           me = hydroGOF::me.data.frame(model2$var1.pred, model2$observed),
                           rmse = hydroGOF::rmse(model2$var1.pred, model2$observed),
                           r2 = hydroGOF::br2(model2$var1.pred, model2$observed),
                           ave = hydroGOF::NSE(model2$var1.pred, model2$observed),
                           willmott = hydroGOF::md(model2$var1.pred, model2$observed)) 
    
    map2 <- idw(z~1, data2, grid, idp = idw_exp2$k)
    
    
    mapa_idw2 <- raster(map2)
    proj4string(mapa_idw2) <- CRS("+init=epsg:32722")
    
    mapas_idw_sparse[[i]] = mapa_idw2
    
    
    
   }
  
  name = paste("output/RData/idw_beta", beta[1], beta[2], beta[3], "anis", anispar[1], anispar[2], ".RData",sep = "_")
  
  save(list = ls(all.names = TRUE), file = name, envir = 
         environment())  


}



