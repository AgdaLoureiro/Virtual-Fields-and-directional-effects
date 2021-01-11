interpo = function(trend,anispar, chute_ini_1, chute_ini_2){
  
  #This code is made to generate a virtual field (VF) using Gaussian Simulation, 
  #sample the VF, interpolate it by Kriging (MoM and REML) and IDW,
  #and to validate interpolation performance
  
  
  #####################Generation of VF###################################
  #We will start by VF generetion
  
  beta = trend
  anis = anispar
  
  #Initial values for semivariogram parameters
  nug_a1 = chute_ini_1[1]
  nug_b1 = chute_ini_1[2]
  nug_c1 = chute_ini_1[3]
  nug_d1 = chute_ini_1[4]
  psil_a1 = chute_ini_1[5] 
  psil_b1 = chute_ini_1[6]
  psil_c1 = chute_ini_1[7]
  psil_d1 = chute_ini_1[8]  
  range_a1 = chute_ini_1[9]
  range_b1 = chute_ini_1[10]
  range_c1 = chute_ini_1[11]
  range_d1 = chute_ini_1[12]
  
  
  nug_a2 = chute_ini_2[1]
  nug_b2 = chute_ini_2[2]
  nug_c2 = chute_ini_2[3]
  nug_d2 = chute_ini_2[4]
  psil_a2 = chute_ini_2[5] 
  psil_b2 = chute_ini_2[6]
  psil_c2 = chute_ini_2[7]
  psil_d2 = chute_ini_2[8]  
  range_a2 = chute_ini_2[9]
  range_b2 = chute_ini_2[10]
  range_c2 = chute_ini_2[11]
  range_d2 = chute_ini_2[12]
  
  
  
  #We will use the following packages for this first part
  library(pacman)
  pacman::p_load(gstat, sp, dplyr,raster, lattice, ggplot2, RStoolbox)
  
  #We will first generate x and y 
  
  # virtual data
  x.range <- c(5.5, 1400.5)
  y.range <- c(5.5, 1450.5)
  area <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 10), 
                      y = seq(from = y.range[1], to = y.range[2], by = 10))
  
  nsim = 100
  #We first define a model of semivariogram and its parameters
  set.seed(93)
  model <- gstat(formula=z ~ x+y, 
                 locations = ~ x+y,
                 dummy= T,#Inconditional simulation term 
                 beta= beta, 
                 model=vgm(nugget = 0, psill= 0.28,
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
  
  #####################Interpolation###################################
  #The necessary packages:
  pacman::p_load(Metrics, geoR, gstat, hydroGOF)
  
  #We have two sample densities. We will perform interpolation for each separately
  #In this part I will divide Geostatistical interpolation into:
  #Interpolation with MoM
  #Interpolation with REML
  #I will also perform interpolation recognizing the following conditions:
  #Condition A - I will not verify if data have anisotropy or trend during semivariogram modeling and kriging
  #Condition B - I will verify only the presence of trend, and use the residual semivariogram
  #Condition C - I will verify only the presence of anisotropy, and model the directional semivariogram
  #Condition D - I will verify trend and anisotropy
  
  #For all interpolations we will use a boundary and get as output a 
  #archive in raster format. 
  
  #For this purpose I draw a boundary of the VF in QGIS:
  boundary <- shapefile("./data/boundary/cotorno.shp")
  r = raster(boundary, res = 10) 
  rp = rasterize(boundary, r, 0) 
  grid <- as(rp, "SpatialPixelsDataFrame")
  
  grid_ge <- as.data.frame(x = rp, xy = T, na.rm = T)
  
  grid.geodata <- as.geodata(grid_ge, data.col = 3,# ID
                             coords.col = c(1:2))
  
  ######----Interpolation----
  #outputs----
  #We fist create a object to store outputs
  #metrics
  mom_a1 = list()
  mom_b1 = list()
  mom_c1 = list()
  mom_d1 = list()
  reml_a1 = list()
  reml_b1 = list()
  reml_c1 = list()
  reml_d1 = list()
  
  mom_a2 = list()
  mom_b2 = list()
  mom_c2 = list()
  mom_d2 = list()
  reml_a2 = list()
  reml_b2 = list()
  reml_c2 = list()
  reml_d2 = list()
  
  #comparing maps
  
  mom_a1_r = data.frame(MAE = NA,
                      RMSE = NA,
                      AVE = NA,
                      ME = NA)
  mom_b1_r = data.frame(MAE = NA,
                      RMSE = NA,
                      AVE = NA,
                      ME = NA)
  mom_c1_r = data.frame(MAE = NA,
                      RMSE = NA,
                      AVE = NA,
                      ME = NA)
  mom_d1_r = data.frame(MAE = NA,
                      RMSE = NA,
                      AVE = NA,
                      ME = NA)
  reml_a1_r = data.frame(MAE = NA,
                       RMSE = NA,
                       AVE = NA,
                       ME = NA)
  reml_b1_r = data.frame(MAE = NA,
                       RMSE = NA,
                       AVE = NA,
                       ME = NA)
  reml_c1_r = data.frame(MAE = NA,
                       RMSE = NA,
                       AVE = NA,
                       ME = NA)
  reml_d1_r = data.frame(MAE = NA,
                       RMSE = NA,
                       AVE = NA,
                       ME = NA)
  
  mom_a2_r = data.frame(MAE = NA,
                      RMSE = NA,
                      AVE = NA,
                      ME = NA)
  mom_b2_r = data.frame(MAE = NA,
                      RMSE = NA,
                      AVE = NA,
                      ME = NA)
  mom_c2_r = data.frame(MAE = NA,
                      RMSE = NA,
                      AVE = NA,
                      ME = NA)
  mom_d2_r = data.frame(MAE = NA,
                      RMSE = NA,
                      AVE = NA,
                      ME = NA)
  reml_a2_r = data.frame(MAE = NA,
                       RMSE = NA,
                       AVE = NA,
                       ME = NA)
  reml_b2_r = data.frame(MAE = NA,
                       RMSE = NA,
                       AVE = NA,
                       ME = NA)
  reml_c2_r = data.frame(MAE = NA,
                       RMSE = NA,
                       AVE = NA,
                       ME = NA)
  reml_d2_r = data.frame(MAE = NA,
                       RMSE = NA,
                       AVE = NA,
                       ME = NA)
  
  
  
  #Semivariogram parameters
  
  fits_mom_a1 = list()
  fits_mom_b1 = list()
  fits_mom_c1 = list()
  fits_mom_d1 = list()
  fits_reml_a1 = list()
  fits_reml_b1 = list()
  fits_reml_c1 = list()
  fits_reml_d1 = list()
  
  
  fits_mom_a2 = list()
  fits_mom_b2 = list()
  fits_mom_c2 = list()
  fits_mom_d2 = list()
  fits_reml_a2 = list()
  fits_reml_b2 = list()
  fits_reml_c2 = list()
  fits_reml_d2 = list()
  
  
  #Interpolated surfaces
  mapa_mom_a1 = list()
  mapa_mom_b1 = list()
  mapa_mom_c1 = list()
  mapa_mom_d1 = list()
  mapa_reml_a1 = list()
  mapa_reml_b1 = list()
  mapa_reml_c1 = list()
  mapa_reml_d1 = list()
  
  mapa_mom_a2 = list()
  mapa_mom_b2 = list()
  mapa_mom_c2 = list()
  mapa_mom_d2 = list()
  mapa_reml_a2 = list()
  mapa_reml_b2 = list()
  mapa_reml_c2 = list()
  mapa_reml_d2 = list()
  
  #Interpo----
  #Now, to ease go through all realizations (Simulations) of the VF,
  #we will use loops

  for(j in 1:4){
    for(i in 3:102) {
      vf_i = vf[,c(1,2,i)]
      vf_i = rasterFromXYZ(vf_i)
      names(vf_i) <- "sim1"
      data1 <- dense_sampling[, c(1,2,i)]
      data2 <-sparse_sampling[, c(1,2,i)]
      names(data1) <- c("x", "y", "z")
      names(data2) <- c("x", "y", "z")
      
      data1_ge = as.geodata(data1, coords.col = c(1,2),
                        data.col = "z")
      data2_ge = as.geodata(data2, coords.col = c(1,2),
                            data.col = "z")
      
      sp::coordinates(data1) = ~x+y
      sp::coordinates(data2) = ~x+y
      proj4string(data1) <- CRS(proj4string(grid))
      proj4string(data2) <- CRS(proj4string(grid))
      print(i)
      
      if (j == 1) {
        #Condition A
        #dense sampling
        semivar1 = gstat::variogram(z ~ 1, data = data1, cutoff= 1200,
                                    width= 100)  
        
        fit.mom1 = gstat::fit.variogram(semivar1, 
                                        vgm(nugget = nug_a1,
                                            psill = psil_a1,
                                            range = range_a1,
                                            model = "Sph"), 
                                        debug.level = 0)
        
        
        val.mom1 = gstat::krige.cv(z ~ 1, locations = data1, model = fit.mom1)
        
        fit.reml1 <- likfit(data1_ge, ini = c(psill = fit.mom1$psill[2],fit.mom1$range[2]), 
                            fix.nugget = T, method = "REML",cov.model = "spherical")
        
        val.reml1 = xvalid(data1_ge, model=fit.reml1)
        
        
        #sparse sampling
        semivar2 = gstat::variogram(z ~ 1, data = data2, cutoff= 2200,
                                    width= 200)  
        
        fit.mom2 = gstat::fit.variogram(semivar2, 
                                        vgm(nugget = nug_a2,
                                            psill = psil_a2,
                                            range = range_a2,
                                            model = "Sph"), 
                                        debug.level = 0)
        
        
        val.mom2 = gstat::krige.cv(z ~ 1, locations = data2, model = fit.mom2)
        
        fit.reml2 <- likfit(data2_ge, ini = c(psill = fit.mom2$psill[2],fit.mom2$range[2]), 
                            fix.nugget = T, method = "REML",cov.model = "spherical")
        
        fits_mom_a1[[i]] = data.frame(fit.mom1)
        fits_reml_a1[[i]] = (fit.reml1)
        
        fits_mom_a2[[i]] = data.frame(fit.mom2)
        fits_reml_a2[[i]] = (fit.reml2)
        
        val.reml2 = xvalid(data2_ge, model=fit.reml2)
        
        
        #prediction
        #dense sampling
        map_mom1 <- raster(krige(z~1, data1, grid, model = fit.mom1))
        proj4string(map_mom1) <- CRS("+init=epsg:32722")
        mapa_mom_a1[[i]] = map_mom1
        
        
        k.df1 <-  krige.conv(geodata = data1_ge, locations = grid.geodata$coords,
                             krige = krige.control(type.krige = "OK",
                                                   obj.model = fit.reml1))
        map_reml1 = data.frame(coords = grid_ge[,1:2], predict = k.df1$predict) 
        map_reml1 <- rasterFromXYZ(map_reml1)
        proj4string(map_reml1) <- CRS("+init=epsg:32722")
        mapa_reml_a1[[i]] = map_reml1
        
        
        #sparse sampling
        map_mom2 <- raster(krige(z~1, data2, grid, model = fit.mom2))
        proj4string(map_mom2) <- CRS("+init=epsg:32722")
        mapa_mom_a2[[i]] = map_mom2
        
        
        k.df2 <-  krige.conv(geodata = data2_ge, locations = grid.geodata$coords,
                             krige = krige.control(type.krige = "OK",
                                                   obj.model = fit.reml2))
        
        map_reml2 = data.frame(coords = grid_ge[,1:2], predict = k.df2$predict) 
        map_reml2 <- rasterFromXYZ(map_reml2)
        proj4string(map_reml2) <- CRS("+init=epsg:32722")
        mapa_reml_a2[[i]] = map_reml2
        
        
        
        
        ##Accessing Metrics
        mom_a1[[i]] = data.frame(mae = hydroGOF::mae(val.mom1$var1.pred, val.mom1$observed),
                                 me = hydroGOF::me.data.frame(val.mom1$var1.pred, val.mom1$observed),
                                 rmse = hydroGOF::rmse(val.mom1$var1.pred, val.mom1$observed),
                                 r2 = hydroGOF::br2(val.mom1$var1.pred, val.mom1$observed),
                                 ave = hydroGOF::NSE(val.mom1$var1.pred, val.mom1$observed),
                                 willmott = hydroGOF::md(val.mom1$var1.pred, val.mom1$observed))
       
        reml_a1[[i]] = data.frame(mae = hydroGOF::mae(val.reml1$predicted, val.reml1$data),
                                  me = hydroGOF::me.data.frame(val.reml1$predicted, val.reml1$data),
                                  rmse = hydroGOF::rmse(val.reml1$predicted, val.reml1$data),
                                  r2 = hydroGOF::br2(val.reml1$predicted, val.reml1$data),
                                  ave = hydroGOF::NSE(val.reml1$predicted, val.reml1$data),
                                  willmott = hydroGOF::md(val.reml1$predicted, val.reml1$data))
        
        
       
        mom_a2[[i]] = data.frame(mae = hydroGOF::mae(val.mom2$var1.pred, val.mom2$observed),
                                 me = hydroGOF::me.data.frame(val.mom2$var1.pred, val.mom2$observed),
                                 rmse = hydroGOF::rmse(val.mom2$var1.pred, val.mom2$observed),
                                 r2 = hydroGOF::br2(val.mom2$var1.pred, val.mom2$observed),
                                 ave = hydroGOF::NSE(val.mom2$var1.pred, val.mom2$observed),
                                 willmott = hydroGOF::md(val.mom2$var1.pred, val.mom2$observed))
        
        reml_a2[[i]] = data.frame(mae = hydroGOF::mae(val.reml2$predicted, val.reml2$data),
                                  me = hydroGOF::me.data.frame(val.reml2$predicted, val.reml2$data),
                                  rmse = hydroGOF::rmse(val.reml2$predicted, val.reml2$data),
                                  r2 = hydroGOF::br2(val.reml2$predicted, val.reml2$data),
                                  ave = hydroGOF::NSE(val.reml2$predicted, val.reml2$data),
                                  willmott = hydroGOF::md(val.reml2$predicted, val.reml2$data))
        
        mom_a1_r[i-2,1] = hydroGOF::mae((as.data.frame(map_mom1$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_a1_r[i-2,4] = hydroGOF::me.data.frame((as.data.frame(map_mom1$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_a1_r[i-2,2] = hydroGOF::rmse((as.data.frame(map_mom1$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_a1_r[i-2,3] = hydroGOF::NSE((as.data.frame(map_mom1$var1.pred)), (as.data.frame(vf_i$sim1)))
        
        reml_a1_r[i-2,1] = hydroGOF::mae((as.data.frame(map_reml1$predict)), (as.data.frame(vf_i$sim1)))
        reml_a1_r[i-2,4] = hydroGOF::me.data.frame((as.data.frame(map_reml1$predict)), (as.data.frame(vf_i$sim1)))
        reml_a1_r[i-2,2] = hydroGOF::rmse((as.data.frame(map_reml1$predict)), (as.data.frame(vf_i$sim1)))
        reml_a1_r[i-2,3] = hydroGOF::NSE((as.data.frame(map_reml1$predict)), (as.data.frame(vf_i$sim1)))
        
        mom_a2_r[i-2,1] = hydroGOF::mae((as.data.frame(map_mom2$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_a2_r[i-2,4] = hydroGOF::me.data.frame((as.data.frame(map_mom2$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_a2_r[i-2,2] = hydroGOF::rmse((as.data.frame(map_mom2$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_a2_r[i-2,3] = hydroGOF::NSE((as.data.frame(map_mom2$var1.pred)), (as.data.frame(vf_i$sim1)))
        
        reml_a2_r[i-2,1] = hydroGOF::mae((as.data.frame(map_reml2$predict)), (as.data.frame(vf_i$sim1)))
        reml_a2_r[i-2,4] = hydroGOF::me.data.frame((as.data.frame(map_reml2$predict)), (as.data.frame(vf_i$sim1)))
        reml_a2_r[i-2,2] = hydroGOF::rmse((as.data.frame(map_reml2$predict)), (as.data.frame(vf_i$sim1)))
        reml_a2_r[i-2,3] = hydroGOF::NSE((as.data.frame(map_reml2$predict)), (as.data.frame(vf_i$sim1)))
        
        
      }
      if (j == 2){
        #Condition B
        #dense sampling
        semivar1 = gstat::variogram(z ~ x+y, data = data1, cutoff= 1200,
                                    width= 100)   
        
        
        fit.mom1 = fit.variogram(semivar1, 
                                 vgm(nugget = nug_b1,
                                     psill = psil_b1,
                                     range = range_b1,
                                     model = "Sph")) 
        
        
        val.mom1 = krige.cv(z ~ x+y, locations = data1, model = fit.mom1)
        
        fit.reml1 <- likfit(data1_ge, ini = c(psill = fit.mom1$psill[2],fit.mom1$range[2]), 
                            fix.nugget = T, method = "REML",cov.model = "spherical", trend = "1st")
        
        val.reml1 = xvalid(data1_ge, model=fit.reml1)
        
        
        fits_mom_b1[[i]] = data.frame(fit.mom1)
        fits_reml_b1[[i]] = (fit.reml1)
        
        #sparse sampling
        semivar2 = gstat::variogram(z ~ x+y, data = data2, cutoff= 1200,
                                    width= 200)   
        
        
        fit.mom2 = fit.variogram(semivar2, 
                                 vgm(nugget = nug_b2,
                                     psill = psil_b2,
                                     range = range_b2,
                                     model = "Sph")) 
        
        
        val.mom2 = krige.cv(z ~ x+y, locations = data2, model = fit.mom2)
        
        fit.reml2 <- likfit(data2_ge, ini = c(psill = fit.mom2$psill[2],fit.mom2$range[2]), 
                            fix.nugget = T, method = "REML",cov.model = "spherical", trend = "1st")
        
        val.reml2 = xvalid(data2_ge, model=fit.reml2)
        
       
        
        fits_mom_b2[[i]] = data.frame(fit.mom2)
        fits_reml_b2[[i]] = (fit.reml2)
        
        
        
        #dense sampling
        map_mom1 <- raster(krige(z~x+y, data1, grid, model = fit.mom1))
        proj4string(map_mom1) <- CRS("+init=epsg:32722")
        mapa_mom_b1[[i]] = map_mom1
        
        
        
        
        k.df1 <-  krige.conv(data1_ge, locations = grid.geodata$coords,
                             krige = krige.control(obj.model = fit.reml1, trend.d = "1st",
                                                   trend.l = "1st"))
        
        map_reml1 = data.frame(coords = grid_ge[,1:2], predict = k.df1$predict) 
        map_reml1 <- rasterFromXYZ(map_reml1)
        proj4string(map_reml1) <- CRS("+init=epsg:32722")
        mapa_reml_b1[[i]] = map_reml1
        
        
        
        #sparse sampling
        map_mom2 <- raster(krige(z~x+y, data2, grid, model = fit.mom2))
        proj4string(map_mom2) <- CRS("+init=epsg:32722")
        mapa_mom_b2[[i]] = map_mom2
        
        
        
        k.df2 <-  krige.conv(data2_ge, locations = grid.geodata$coords,
                             krige = krige.control(obj.model = fit.reml2, trend.d = "1st",
                                                   trend.l = "1st"))
        
        map_reml2 = data.frame(coords = grid_ge[,1:2], predict = k.df2$predict) 
        map_reml2 <- rasterFromXYZ(map_reml2)
        proj4string(map_reml2) <- CRS("+init=epsg:32722")
        mapa_reml_b2[[i]] = map_reml2
        
        
        
        ##Accessing Metrics
        mom_b1[[i]] = data.frame(mae = hydroGOF::mae(val.mom1$var1.pred, val.mom1$observed),
                                 me = hydroGOF::me.data.frame(val.mom1$var1.pred, val.mom1$observed),
                                 rmse = hydroGOF::rmse(val.mom1$var1.pred, val.mom1$observed),
                                 r2 = hydroGOF::br2(val.mom1$var1.pred, val.mom1$observed),
                                 ave = hydroGOF::NSE(val.mom1$var1.pred, val.mom1$observed),
                                 willmott = hydroGOF::md(val.mom1$var1.pred, val.mom1$observed))
        
        reml_b1[[i]] = data.frame(mae = hydroGOF::mae(val.reml1$predicted, val.reml1$data),
                                  me = hydroGOF::me.data.frame(val.reml1$predicted, val.reml1$data),
                                  rmse = hydroGOF::rmse(val.reml1$predicted, val.reml1$data),
                                  r2 = hydroGOF::br2(val.reml1$predicted, val.reml1$data),
                                  ave = hydroGOF::NSE(val.reml1$predicted, val.reml1$data),
                                  willmott = hydroGOF::md(val.reml1$predicted, val.reml1$data))
        
        
        
        ##Accessing Metrics
        
        mom_b2[[i]] = data.frame(mae = hydroGOF::mae(val.mom2$var1.pred, val.mom2$observed),
                                 me = hydroGOF::me.data.frame(val.mom2$var1.pred, val.mom2$observed),
                                 rmse = hydroGOF::rmse(val.mom2$var1.pred, val.mom2$observed),
                                 r2 = hydroGOF::br2(val.mom2$var1.pred, val.mom2$observed),
                                 ave = hydroGOF::NSE(val.mom2$var1.pred, val.mom2$observed),
                                 willmott = hydroGOF::md(val.mom2$var1.pred, val.mom2$observed))
        
        reml_b2[[i]] = data.frame(mae = hydroGOF::mae(val.reml2$predicted, val.reml2$data),
                                  me = hydroGOF::me.data.frame(val.reml2$predicted, val.reml2$data),
                                  rmse = hydroGOF::rmse(val.reml2$predicted, val.reml2$data),
                                  r2 = hydroGOF::br2(val.reml2$predicted, val.reml2$data),
                                  ave = hydroGOF::NSE(val.reml2$predicted, val.reml2$data),
                                  willmott = hydroGOF::md(val.reml2$predicted, val.reml2$data))
        
        
        mom_b1_r[i-2,1] = hydroGOF::mae((as.data.frame(map_mom1$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_b1_r[i-2,4] = hydroGOF::me.data.frame((as.data.frame(map_mom1$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_b1_r[i-2,2] = hydroGOF::rmse((as.data.frame(map_mom1$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_b1_r[i-2,3] = hydroGOF::NSE((as.data.frame(map_mom1$var1.pred)), (as.data.frame(vf_i$sim1)))
        
        reml_b1_r[i-2,1] = hydroGOF::mae((as.data.frame(map_reml1$predict)), (as.data.frame(vf_i$sim1)))
        reml_b1_r[i-2,4] = hydroGOF::me.data.frame((as.data.frame(map_reml1$predict)), (as.data.frame(vf_i$sim1)))
        reml_b1_r[i-2,2] = hydroGOF::rmse((as.data.frame(map_reml1$predict)), (as.data.frame(vf_i$sim1)))
        reml_b1_r[i-2,3] = hydroGOF::NSE((as.data.frame(map_reml1$predict)), (as.data.frame(vf_i$sim1)))
        
        mom_b2_r[i-2,1] = hydroGOF::mae((as.data.frame(map_mom2$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_b2_r[i-2,4] = hydroGOF::me.data.frame((as.data.frame(map_mom2$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_b2_r[i-2,2] = hydroGOF::rmse((as.data.frame(map_mom2$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_b2_r[i-2,3] = hydroGOF::NSE((as.data.frame(map_mom2$var1.pred)), (as.data.frame(vf_i$sim1)))
        
        reml_b2_r[i-2,1] = hydroGOF::mae((as.data.frame(map_reml2$predict)), (as.data.frame(vf_i$sim1)))
        reml_b2_r[i-2,4] = hydroGOF::me.data.frame((as.data.frame(map_reml2$predict)), (as.data.frame(vf_i$sim1)))
        reml_b2_r[i-2,2] = hydroGOF::rmse((as.data.frame(map_reml2$predict)), (as.data.frame(vf_i$sim1)))
        reml_b2_r[i-2,3] = hydroGOF::NSE((as.data.frame(map_reml2$predict)), (as.data.frame(vf_i$sim1)))
        
      }
      if (j == 3){
        #Condition C
        #dense sampling
        semivar1 <-  gstat::variogram(z ~ 1, data = data1, cutoff= 1200,
                                      width= 100, alpha =30) 
        
        
        fit.mom1 = fit.variogram(semivar1, 
                                 vgm(nugget = nug_c1,
                                     psill = psil_c1,
                                     range = range_c1,
                                     model = "Sph",
                                     anis = anis)) 
        
        
        
        val.mom1 = krige.cv(z ~ 1, locations = data1, model = fit.mom1)
        
        fit.reml1 <- likfit(data1_ge, ini = c(psill = fit.mom1$psill[2],fit.mom1$range[2]), 
                            fix.nugget = T, method = "REML",cov.model = "spherical",fix.psiA = FALSE, fix.psiR = FALSE)
        
        val.reml1 = xvalid(data1_ge, model=fit.reml1)
        
        
        
        
       
        fits_mom_c1[[i]] = data.frame(fit.mom1)
        fits_reml_c1[[i]] = (fit.reml1)
        
        #sparse sampling
        semivar2 <-  gstat::variogram(z ~ 1, data = data2, cutoff= 1200,
                                      width= 200, alpha =30) 
        
        
        fit.mom2 = fit.variogram(semivar2, 
                                 vgm(nugget = nug_c2,
                                     psill = psil_c2,
                                     range = range_c2,
                                     model = "Sph",
                                     anis = anis)) 
        
        
        
        val.mom2 = krige.cv(z ~ 1, locations = data2, model = fit.mom2)
        
        fit.reml2 <- likfit(data2_ge, ini = c(psill = fit.mom2$psill[2],fit.mom2$range[2]), 
                            fix.nugget = T, method = "REML",cov.model = "spherical",fix.psiA = FALSE, fix.psiR = FALSE)
        
        val.reml2 = xvalid(data2_ge, model=fit.reml2)
         
        
       
        fits_mom_c2[[i]] = data.frame(fit.mom2)
        fits_reml_c2[[i]] = (fit.reml2)
        
        
        #prediction
        #dense sampling
        
        map_mom1 <- raster(krige(z~1, data1, grid, model = fit.mom1))
        proj4string(map_mom1) <- CRS("+init=epsg:32722")
        mapa_mom_c1[[i]] = map_mom1
        
        
        
        k.df1 <-  krige.conv(geodata = data1_ge, locations = grid.geodata$coords,
                             krige = krige.control(type.krige = "OK",
                                                   obj.model = fit.reml1,
                                                   aniso.pars = c(30,2)))
        
        map_reml1 = data.frame(coords = grid_ge[,1:2], predict = k.df1$predict) 
        map_reml1 <- rasterFromXYZ(map_reml1)
        proj4string(map_reml1) <- CRS("+init=epsg:32722")
        mapa_reml_c1[[i]] = map_reml1
        
        
        
        #sparse sampling
        
        map_mom2 <- raster(krige(z~1, data2, grid, model = fit.mom2))
        proj4string(map_mom2) <- CRS("+init=epsg:32722")
        mapa_mom_c2[[i]] = map_mom2
        
        
        k.df2 <-  krige.conv(geodata = data2_ge, locations = grid.geodata$coords,
                             krige = krige.control(type.krige = "OK",
                                                   obj.model = fit.reml2,
                                                   aniso.pars = c(30,2)))
        
        map_reml2 = data.frame(coords = grid_ge[,1:2], predict = k.df2$predict) 
        map_reml2 <- rasterFromXYZ(map_reml2)
        proj4string(map_reml2) <- CRS("+init=epsg:32722")
        mapa_reml_c2[[i]] = map_reml2
        
        
        ##Accessing Metrics
        mom_c1[[i]] = data.frame(mae = hydroGOF::mae(val.mom1$var1.pred, val.mom1$observed),
                                 me = hydroGOF::me.data.frame(val.mom1$var1.pred, val.mom1$observed),
                                 rmse = hydroGOF::rmse(val.mom1$var1.pred, val.mom1$observed),
                                 r2 = hydroGOF::br2(val.mom1$var1.pred, val.mom1$observed),
                                 ave = hydroGOF::NSE(val.mom1$var1.pred, val.mom1$observed),
                                 willmott = hydroGOF::md(val.mom1$var1.pred, val.mom1$observed))
        
        reml_c1[[i]] = data.frame(mae = hydroGOF::mae(val.reml1$predicted, val.reml1$data),
                                  me = hydroGOF::me.data.frame(val.reml1$predicted, val.reml1$data),
                                  rmse = hydroGOF::rmse(val.reml1$predicted, val.reml1$data),
                                  r2 = hydroGOF::br2(val.reml1$predicted, val.reml1$data),
                                  ave = hydroGOF::NSE(val.reml1$predicted, val.reml1$data),
                                  willmott = hydroGOF::md(val.reml1$predicted, val.reml1$data))
        
        
        
        ##Accessing Metrics
        
        mom_c2[[i]] = data.frame(mae = hydroGOF::mae(val.mom2$var1.pred, val.mom2$observed),
                                 me = hydroGOF::me.data.frame(val.mom2$var1.pred, val.mom2$observed),
                                 rmse = hydroGOF::rmse(val.mom2$var1.pred, val.mom2$observed),
                                 r2 = hydroGOF::br2(val.mom2$var1.pred, val.mom2$observed),
                                 ave = hydroGOF::NSE(val.mom2$var1.pred, val.mom2$observed),
                                 willmott = hydroGOF::md(val.mom2$var1.pred, val.mom2$observed))
        
        reml_c2[[i]] = data.frame(mae = hydroGOF::mae(val.reml2$predicted, val.reml2$data),
                                  me = hydroGOF::me.data.frame(val.reml2$predicted, val.reml2$data),
                                  rmse = hydroGOF::rmse(val.reml2$predicted, val.reml2$data),
                                  r2 = hydroGOF::br2(val.reml2$predicted, val.reml2$data),
                                  ave = hydroGOF::NSE(val.reml2$predicted, val.reml2$data),
                                  willmott = hydroGOF::md(val.reml2$predicted, val.reml2$data))
        
        mom_c1_r[i-2,1] = hydroGOF::mae((as.data.frame(map_mom1$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_c1_r[i-2,4] = hydroGOF::me.data.frame((as.data.frame(map_mom1$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_c1_r[i-2,2] = hydroGOF::rmse((as.data.frame(map_mom1$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_c1_r[i-2,3] = hydroGOF::NSE((as.data.frame(map_mom1$var1.pred)), (as.data.frame(vf_i$sim1)))
        
        reml_c1_r[i-2,1] = hydroGOF::mae((as.data.frame(map_reml1$predict)), (as.data.frame(vf_i$sim1)))
        reml_c1_r[i-2,4] = hydroGOF::me.data.frame((as.data.frame(map_reml1$predict)), (as.data.frame(vf_i$sim1)))
        reml_c1_r[i-2,2] = hydroGOF::rmse((as.data.frame(map_reml1$predict)), (as.data.frame(vf_i$sim1)))
        reml_c1_r[i-2,3] = hydroGOF::NSE((as.data.frame(map_reml1$predict)), (as.data.frame(vf_i$sim1)))
        
        mom_c2_r[i-2,1] = hydroGOF::mae((as.data.frame(map_mom2$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_c2_r[i-2,4] = hydroGOF::me.data.frame((as.data.frame(map_mom2$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_c2_r[i-2,2] = hydroGOF::rmse((as.data.frame(map_mom2$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_c2_r[i-2,3] = hydroGOF::NSE((as.data.frame(map_mom2$var1.pred)), (as.data.frame(vf_i$sim1)))
        
        reml_c2_r[i-2,1] = hydroGOF::mae((as.data.frame(map_reml2$predict)), (as.data.frame(vf_i$sim1)))
        reml_c2_r[i-2,4] = hydroGOF::me.data.frame((as.data.frame(map_reml2$predict)), (as.data.frame(vf_i$sim1)))
        reml_c2_r[i-2,2] = hydroGOF::rmse((as.data.frame(map_reml2$predict)), (as.data.frame(vf_i$sim1)))
        reml_c2_r[i-2,3] = hydroGOF::NSE((as.data.frame(map_reml2$predict)), (as.data.frame(vf_i$sim1)))
        
      }
      if (j == 4){
        #Condition D
        #dense sampling
        semivar1 <- gstat::variogram(z ~ x+y, data = data1, cutoff= 1200,
                                     width= 100, alpha =30) 
        
        fit.mom1 =  fit.variogram(semivar1,
                                  vgm(nugget = nug_d1,
                                      psill = psil_d1,
                                      range = range_d1,
                                      model = "Sph",
                                      anis = anis)) 
        
        
        val.mom1 = krige.cv(z ~ x+y, locations = data1, model = fit.mom1)
        
        fit.reml1 <- likfit(data1_ge, ini = c(psill = fit.mom1$psill[2],fit.mom1$range[2]), 
                            fix.nugget = T, method = "REML",cov.model = "spherical", trend = "1st",
                            fix.psiA = FALSE, fix.psiR = FALSE)
        
        val.reml1 = xvalid(data1_ge, model=fit.reml1)
        
        
        fits_mom_d1[[i]] = data.frame(fit.mom1)
        fits_reml_d1[[i]] = (fit.reml1)
        
        #sparse sampling
        semivar2 <- gstat::variogram(z ~ x+y, data = data2, cutoff= 1200,
                                     width= 200, alpha =30) 
        
        fit.mom2 =  fit.variogram(semivar2,
                                  vgm(nugget = nug_d2,
                                      psill = psil_d2,
                                      range = range_d2,
                                      model = "Sph",
                                      anis = anis)) 
        
        
        val.mom2 = krige.cv(z ~ x+y, locations = data2, model = fit.mom2)
        
        fit.reml2 <- likfit(data2_ge, ini = c(psill = fit.mom2$psill[2],fit.mom2$range[2]), 
                            fix.nugget = T, method = "REML",cov.model = "spherical", trend = "1st",
                            fix.psiA = FALSE, fix.psiR = FALSE)
        
        val.reml2 = xvalid(data2_ge, model=fit.reml2)
        
       
        fits_mom_d2[[i]] = data.frame(fit.mom2)
        fits_reml_d2[[i]] = (fit.reml2)
        
        
        #prediction
        #dense sampling
        map_mom1 <- raster(krige(z~x+y, data1, grid, model = fit.mom1))
        proj4string(map_mom1) <- CRS("+init=epsg:32722")
        mapa_mom_d1[[i]] = map_mom1
        
        
        k.df1 <-  krige.conv(geodata = data1_ge, locations = grid.geodata$coords,
                             krige = krige.control(obj.model = fit.reml1, trend.d = "1st",
                                                   trend.l = "1st",aniso.pars = c(30,2)))
        
        map_reml1 = data.frame(coords = grid_ge[,1:2], predict = k.df1$predict) 
        map_reml1 <- rasterFromXYZ(map_reml1)
        proj4string(map_reml1) <- CRS("+init=epsg:32722")
        mapa_reml_d1[[i]] = map_reml1
        
        
        
        #sparse sampling
        
        map_mom2 <- raster(krige(z~x+y, data2, grid, model = fit.mom2))
        proj4string(map_mom2) <- CRS("+init=epsg:32722")
        mapa_mom_d2[[i]] = map_mom2
        
        
        
        k.df2 <-  krige.conv(geodata = data2_ge, locations = grid.geodata$coords,
                             krige = krige.control(obj.model = fit.reml2, trend.d = "1st",
                                                   trend.l = "1st",aniso.pars = c(30,2)))
        
        map_reml2 = data.frame(coords = grid_ge[,1:2], predict = k.df2$predict) 
        map_reml2 <- rasterFromXYZ(map_reml2)
        proj4string(map_reml2) <- CRS("+init=epsg:32722")
        mapa_reml_d2[[i]] = map_reml2
        
        
        
        ##Accessing Metrics
        mom_d1[[i]] = data.frame(mae = hydroGOF::mae(val.mom1$var1.pred, val.mom1$observed),
                                 me = hydroGOF::me.data.frame(val.mom1$var1.pred, val.mom1$observed),
                                 rmse = hydroGOF::rmse(val.mom1$var1.pred, val.mom1$observed),
                                 r2 = hydroGOF::br2(val.mom1$var1.pred, val.mom1$observed),
                                 ave = hydroGOF::NSE(val.mom1$var1.pred, val.mom1$observed),
                                 willmott = hydroGOF::md(val.mom1$var1.pred, val.mom1$observed))
        
        reml_d1[[i]] = data.frame(mae = hydroGOF::mae(val.reml1$predicted, val.reml1$data),
                                  me = hydroGOF::me.data.frame(val.reml1$predicted, val.reml1$data),
                                  rmse = hydroGOF::rmse(val.reml1$predicted, val.reml1$data),
                                  r2 = hydroGOF::br2(val.reml1$predicted, val.reml1$data),
                                  ave = hydroGOF::NSE(val.reml1$predicted, val.reml1$data),
                                  willmott = hydroGOF::md(val.reml1$predicted, val.reml1$data))
        
        
        
        ##Accessing Metrics
        
        mom_d2[[i]] = data.frame(mae = hydroGOF::mae(val.mom2$var1.pred, val.mom2$observed),
                                 me = hydroGOF::me.data.frame(val.mom2$var1.pred, val.mom2$observed),
                                 rmse = hydroGOF::rmse(val.mom2$var1.pred, val.mom2$observed),
                                 r2 = hydroGOF::br2(val.mom2$var1.pred, val.mom2$observed),
                                 ave = hydroGOF::NSE(val.mom2$var1.pred, val.mom2$observed),
                                 willmott = hydroGOF::md(val.mom2$var1.pred, val.mom2$observed))
        
        reml_d2[[i]] = data.frame(mae = hydroGOF::mae(val.reml2$predicted, val.reml2$data),
                                  me = hydroGOF::me.data.frame(val.reml2$predicted, val.reml2$data),
                                  rmse = hydroGOF::rmse(val.reml2$predicted, val.reml2$data),
                                  r2 = hydroGOF::br2(val.reml2$predicted, val.reml2$data),
                                  ave = hydroGOF::NSE(val.reml2$predicted, val.reml2$data),
                                  willmott = hydroGOF::md(val.reml2$predicted, val.reml2$data))
        
        mom_d1_r[i-2,1] = hydroGOF::mae((as.data.frame(map_mom1$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_d1_r[i-2,4] = hydroGOF::me.data.frame((as.data.frame(map_mom1$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_d1_r[i-2,2] = hydroGOF::rmse((as.data.frame(map_mom1$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_d1_r[i-2,3] = hydroGOF::NSE((as.data.frame(map_mom1$var1.pred)), (as.data.frame(vf_i$sim1)))
        
        reml_d1_r[i-2,1] = hydroGOF::mae((as.data.frame(map_reml1$predict)), (as.data.frame(vf_i$sim1)))
        reml_d1_r[i-2,4] = hydroGOF::me.data.frame((as.data.frame(map_reml1$predict)), (as.data.frame(vf_i$sim1)))
        reml_d1_r[i-2,2] = hydroGOF::rmse((as.data.frame(map_reml1$predict)), (as.data.frame(vf_i$sim1)))
        reml_d1_r[i-2,3] = hydroGOF::NSE((as.data.frame(map_reml1$predict)), (as.data.frame(vf_i$sim1)))
        
        mom_d2_r[i-2,1] = hydroGOF::mae((as.data.frame(map_mom2$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_d2_r[i-2,4] = hydroGOF::me.data.frame((as.data.frame(map_mom2$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_d2_r[i-2,2] = hydroGOF::rmse((as.data.frame(map_mom2$var1.pred)), (as.data.frame(vf_i$sim1)))
        mom_d2_r[i-2,3] = hydroGOF::NSE((as.data.frame(map_mom2$var1.pred)), (as.data.frame(vf_i$sim1)))
        
        reml_d2_r[i-2,1] = hydroGOF::mae((as.data.frame(map_reml2$predict)), (as.data.frame(vf_i$sim1)))
        reml_d2_r[i-2,4] = hydroGOF::me.data.frame((as.data.frame(map_reml2$predict)), (as.data.frame(vf_i$sim1)))
        reml_d2_r[i-2,2] = hydroGOF::rmse((as.data.frame(map_reml2$predict)), (as.data.frame(vf_i$sim1)))
        reml_d2_r[i-2,3] = hydroGOF::NSE((as.data.frame(map_reml2$predict)), (as.data.frame(vf_i$sim1)))
        
      }
      
    }
  }
  
  name = paste("output/RData/beta", beta[1], beta[2], beta[3], "anis", anispar[1], anispar[2], ".RData",sep = "_")
  
  save(list = ls(all.names = TRUE), file = name, envir = 
         environment())
  
}


