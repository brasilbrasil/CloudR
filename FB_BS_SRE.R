rm(list = ls()) #remove all past worksheet variables
source(paste0("C:/Users/lfortini/","directory_registry.r"))

library(raster)
library(stringr)
library(biomod2)

###USER CONFIGURATION
overwrite=1
min_n_points=2
working_dir=paste0(analysisDir,"/global_clim/") #this is where r will create a mirror folder structure with the 
setwd(working_dir)

clim_data_dir0=paste0(bioclimData2013Dir,"all_baseline/500m/")  
clim_data_dir1=paste0(bioclimData2013Dir,"all_future/500m/")

clim_dirs=c(clim_data_dir0, clim_data_dir1)
clim_names=c('biovars2000', 'biovars2100')
env_var_files=c("bio1", "bio7", "bio12", "bio15") 
i=1
for (clim_data_dir in clim_dirs){
  predictors = raster( paste(clim_data_dir, env_var_files[1], ".tif", sep=""))
  for (jj in 2:length(env_var_files)){
    temp=raster(paste(clim_data_dir, env_var_files[jj], ".tif", sep=""))
    predictors = addLayer(predictors, temp)
  }
  names(predictors)<- env_var_files
  assign(clim_names[i],predictors)
  rm("predictors" ,"temp") 
  i=i+1
}

source("SRE_functions.r")

dir.create(paste0(working_dir, "tifs/"),showWarnings=F)
dir.create(paste0(working_dir, 'jpg_outputs/'),showWarnings=F)


#############
#forest birds
#############
all_spp=c('Akekee','Akiapolauu','Apapane','Akohekohe', 'Hawaii_Akepa', 'Iiwi', 'Hawaii_Creeper', 'Hawaii_Elepaio', 'Kauai_Elepaio', 'Amakihi', 'Maui_Alauahio', 'Hawaii_Amakihi', 'Maui_Parrotbill', 'Oahu_Amakihi',  'Omao', 'Oahu_Elepaio', 'Palila', 'Elepaio', 'Puaiohi', 'Kauai_Amakihi', 'Anianiau', 'Akikiki')

source(paste0("C:/Users/lfortini/","directory_registry.r"))
spdata_dir=paste0(resultsDir,'/necessary_run_data/')
csv_dir=paste(spdata_dir,"single_sp_CSVs/", sep="")
pct_data=2 #if this is <1 value will be used as proportion of total points to be sampled at each boot strap; if >1, value will be used directly as n to sample

for (sp in all_spp){
  ptm0 <- proc.time()
  cat('\n',sp,'modeling...')
  spp_data=read.csv(paste(csv_dir,sp,'_pres_abs.csv', sep = "")) #FB_data_points4_PAandA
  sp_data=data.frame(cbind(x=spp_data$X, y=spp_data$Y, pres=spp_data$pa))
  Response_var=sp_data[sp_data[,"pres"]==1,]
  reps=100 #200 for BS
  sp_str=sp
  sp_name=sp
  my_sre_BS(Response_var, sp_str=sp, reps, biovars2000, biovars2100,pct_data, overwrite=F)
  ptm1=proc.time() - ptm0
  jnk=as.numeric(ptm1[3])
  jnk=jnk/60
  cat('\n','It took ', jnk, "minutes to model", sp_name)
}


