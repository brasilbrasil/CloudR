
sre_projection=function(NewData, ExtremCond){
  if (is.data.frame(NewData) | is.matrix(NewData)) {
    out <- rep(1, nrow(NewData))
    for (j in 1:ncol(NewData)) {
      out <- out * as.numeric(NewData[, j] >= ExtremCond[j, 
                                                         1] & NewData[, j] <= ExtremCond[j, 2])
    }
  }
  if (inherits(NewData, "Raster")) {
    out <- reclassify(raster:::subset(NewData, 1, drop = TRUE), 
                      c(-Inf, Inf, 1))
    for (j in 1:nlayers(NewData)) {
      out <- out * (raster:::subset(NewData, j, drop = TRUE) >= 
                      ExtremCond[j, 1]) * (raster:::subset(NewData, 
                                                           j, drop = TRUE) <= ExtremCond[j, 2])
    }
    out <- raster:::subset(out, 1, drop = TRUE)
  }
  return(out)
}

sre_dataframe=function(Response,explanatory, NewData, Quant){  
  require(biomod2)
  
  lout <- list()
  if (is.data.frame(Response) | is.matrix(Response)) {
    resp.names <- colnames(Response)[3]
    occ.pts <- which(Response[, 3] == 1)
    extrem.cond <- t(apply(as.data.frame(explanatory[occ.pts,]), 2, quantile, probs = c(0 + Quant, 1 - Quant), na.rm = TRUE))
    lout[[1]] <- sre_projection(NewData, extrem.cond)
  }  
  lout <- stack(lout)
  if (nlayers(lout) == 1) {
    lout <- raster:::subset(lout, 1, drop = TRUE)
  }
  names(lout) <- resp.names
  return(lout)
}

BS_sre=function(Response,explanatory, NewData, pct_data, reps){
  #pct_data
  if (reps==1){
    Quant=0.005 #0.0001 #for 
  }else{
    Quant=0.0001 #0.0001 #for     
  }
  #subsample points
  nr<-dim(Response)[1]
  if (pct_data<=1){
    nr_s<-dim(Response)[1]*pct_data    
  }else{
    nr_s<-pct_data    
  }
  Response_sample=Response[sample.int(nr,nr_s),]
  explanatory_sample=explanatory[sample.int(nr,nr_s),]
  
  #apply sre
  CE=sre_dataframe(Response_sample,explanatory_sample, NewData, Quant)
  return(CE)
}

###parallel computation

Parallel_fx=function(x){
  library(rgdal)
  cat('doing rep ', x,'\n')
  out=BS_sre(Response_var, explanatory, biovars2000, pct_data, reps)  
  return(out)
}

parallel_fx_run=function(reps, fx){  
  require(snowfall)
  # Init Snowfall with settings from sfCluster
  cpucores=as.integer(Sys.getenv('NUMBER_OF_PROCESSORS'))
  sfInit( parallel=TRUE, cpus=cpucores) # 
  sfExportAll( except=c( "working_dir" ) )
  jnk=c(1:reps)
  all_reps=sfLapply(x=jnk,fun=fx)
  sfRemoveAll( except=c( "rep_FX" ) )
  sfStop()
  cat('summarizing all BS CEs for ',sp_str,'\n')
  #Sum all rasters
  all_reps=Reduce("+",all_reps) #much faster than sum of stack or do.call
  all_reps=all_reps/reps #divide by number of reps
  return(all_reps)
}

Process_BS_SRE_data=function(all_reps,plot_var){
  jpeg_name=paste("jpg_outputs/", sp_str, "_",plot_var,".jpg", sep = "")
  jpeg(jpeg_name,
       width = 10, height = 10, units = "in",
       pointsize = 12, quality = 90, bg = "white", res = 300)
  plot(all_reps)
  dev.off()  
  tifname=paste0("tifs/", sp_str, "_",plot_var,".tif") 
  writeRaster(all_reps, tifname, format="GTiff", overwrite=TRUE)
}



my_sre_BS=function(Response_var, sp_str, reps, biovars2000, biovars2100, pct_data, overwrite=F){
  #response_var is dataframe with XY of occurrence
  #sp_str is name of species  
  sp_str0=sp_str #debug
  if (dim(Response_var)[1]<min_n_points){
    cat("not enough points for modeling ", sp_str, "\n")
    break #FIX THIS
  }
  
  out_raster_name=paste('tifs/', sp_str,"_response_zones.tif", sep = "")
  BS_out_raster_name=paste0("tifs/", sp_str, "_BS_CE_change.tif") #
  #BS_out_raster_name=paste0("jpg_outputs/", sp_str, "_BS_CE_change.jpg") #
  
  if ((reps==1 & file.exists(out_raster_name)==F)|(reps>1 & file.exists(BS_out_raster_name)==F)|overwrite==T){    
    explanatory=extract(biovars2000, Response_var[,1:2])
    n_points=dim(Response_var)[1]
    #present
    if (reps==1){
      pred <- sre_dataframe(Response_var, explanatory, biovars2000, Quant = 0.001)
      SRE_raster_present=subset(pred, 3)
      plot_var="SRE_present"    
      
    }else{
      Parallel_fx=function(x){
        library(rgdal)
        cat('doing rep ', x,'\n')
        out=BS_sre(Response_var, explanatory, biovars2000, pct_data = pct_data)  
        return(out)
      }
      SRE_raster_present=parallel_fx_run(reps, fx=Parallel_fx)
      plot_var="BS_SRE_present"    
    }
    Process_BS_SRE_data(SRE_raster_present, plot_var)
    
    #future
    if (reps==1){
      pred <- sre_dataframe(Response_var, explanatory, biovars2100, Quant = 0.001)
      SRE_raster_future=subset(pred, 3)
      plot_var="SRE_future"    
      
    }else{
      Parallel_fx=function(x){
        library(rgdal)
        cat('doing rep ', x,'\n')
        out=BS_sre(Response_var, explanatory, biovars2100, pct_data = pct_data)  
        return(out)
      }
      SRE_raster_future=parallel_fx_run(reps, fx=Parallel_fx)
      plot_var="BS_SRE_future"    
    }
    Process_BS_SRE_data(SRE_raster_future, plot_var)
    
    sp_str=sp_str0 #debug    
    if (reps==1){
      jnk=SRE_raster_future*10
      BIN_dif=SRE_raster_present+jnk
      m  =  c(9.9,  10.1,  3, 10.9, 11.1, 2)
      rclmat  =  matrix(m,  ncol=3,  byrow=TRUE)
      resp_zone  =  reclassify(BIN_dif,  rclmat)
      
      mypalette_numbers=c(0, 1, 2, 3)
      mypalette=c("Grey", "Red", "Green", "Yellow")
      resp_zone_names0=c("Micro refugia", "Tolerate", "Migrate")
      
      jnk=unique(resp_zone)
      zones_present=jnk[jnk>0]
      zones_present=zones_present[zones_present<=3]
      resp_zone_colors=mypalette[zones_present+1]
      resp_zone_names=resp_zone_names0[zones_present]
      
      jpeg_name=paste('jpg_outputs/', sp_str,"_response_zones.jpg", sep = "")
      jpeg(jpeg_name,
           width = 10, height = 8, units = "in",
           pointsize = 12, quality = 90, bg = "white", res = 300)
      plot(resp_zone,  col=mypalette, legend=F)
      #legend("bottomleft",legend = c("Micro refugia", "Tolerate", "Migrate"), col = mypalette[2:4],pch = 16)
      legend("bottomleft",legend = resp_zone_names, col = resp_zone_colors,pch = 16)
      title(paste(sp_name, " response zones (n=", n_points, ")", sep=""))
      dev.off()
      writeRaster(resp_zone, out_raster_name, format="GTiff", overwrite=TRUE) 
      
    }else{
      jnk=SRE_raster_future-SRE_raster_present
      plot_var="BS_CE_change"
      Process_BS_SRE_data(jnk, plot_var)
    }
  }else{
    cat(sp_str," already done", "\n")
  }
  
}
