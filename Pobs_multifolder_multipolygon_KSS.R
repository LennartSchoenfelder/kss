library(dplyr)
library(ncdf4)
library(tictoc)
library(sf)
library(xts)
library(stringr)


# #this script clips gridded Precipitation from netCDF files to subbasins in multi-polygons
# output is a "Pobs.txt" for direct use in a HYPE model setup
##  hist_CNRM_CCLM_RR_daily_mm_1990.nc is an example filename

#set to where delineation is located
dir.gis = "D:/OneDrive - NTNU/HYPE_T/GIS/"
#dir.nc = "C:/Users/schoe/Desktop/KSS_1785_97095/KSSData/" #all folders with scenarios and variables in this directory
#set work directory to path where climatedata is located
dir.nc = "H:/RR/"

varname_nc = "RR"

#create vector with scenario names 
v.scenarios = c("hist", "rcp45", "rcp85")

setwd(dir.nc)
# shell.exec(getwd())
# import shapefile (polygons) that delineate the subcatchments (exportable from WHIST)
#see WHIST from SMHI for HYPE
#Polygons can be in any projection, projection has to be changed to nc file projection
#subbasins <- readOGR(dsn = "C:/Users/lennarts/OneDrive - SINTEF/SINTEF/PHD/HYPE_T/GIS", layer= "Gaula_delineation.shp", encoding = "UTF8")
subbasins = st_read(dsn = paste0(dir.gis,"Gaula_delineation_2.gpkg"))
st_crs(subbasins)
subbasins <- st_transform(subbasins, crs = 32633)
subbasins$ROWNR = seq(1:nrow(subbasins))
st_crs(subbasins)
getwd()
############################################################################
############################################################################

# # #load first NC PREC data to create spatial structure of dataset
setwd(dir.nc)
v.folder = list.dirs(path = ".", full.names = F, recursive = F)
init.filelist <- list.files(path = paste0("./",v.folder[1], "/", v.scenarios[1],"/"), pattern = "*.nc")

nc <- nc_open(paste0("./",v.folder[1], "/", v.scenarios[1],"/",init.filelist[1]), verbose = F, write = FALSE)


nc1 = nc_open("./CNRM_CCLM/hist/hist_CNRM_CCLM_RR_daily_1971_v4.nc", verbose = F, write = FALSE)
nc2 = nc_open("./CNRM_CCLM/rcp45/rcp45_CNRM_CCLM_RR_daily_2006_v4.nc", verbose = F, write = FALSE)
nc3 = nc_open("./CNRM_CCLM/rcp85/rcp85_CNRM_CCLM_RR_daily_2006_v4.nc", verbose = F, write = FALSE)
names(nc1$var)[5]
names(nc2$var)[5]
names(nc3$var)[5]
#check if projection is correct and find variable name
show(nc)
nc$dim$time
############################################################################
############################################################################
#select variable
nc$var$precipitation__map_hist_daily
varname = "precipitation__map_hist_daily"

#lonlat by values for wgs84 utm 33 projection from VALS
lat <- ncvar_get(nc,"Yc")
lon <- ncvar_get(nc,"Xc")
lonlat <- as.data.frame(merge(lon,lat))
rm(lat, lon)
nc_close(nc)
# #check the cornerpoints of your domain (optional)
# #left top
# lonlat[1,]
# #right top
# lonlat[nc$dim$X$len,]
# #left bottom
# lonlat[(nc$dim$Y$len-1)*nc$dim$X$len+1,]
# #right bottom
# lonlat[nc$dim$X$len*nc$dim$Y$len,]
#create matrix of indices, this will be needed to access the NC file rapidly
lonind <- as.double(rep(1:nc$dim$Xc$len, times = nc$dim$Yc$len))
latind <- as.double(rep(1:nc$dim$Yc$len, each = nc$dim$Xc$len))
lonlatind <- cbind (lonind, latind)
lonlatind <- as.data.frame(lonlatind)
coords_ind_mat <- cbind(lonlat,lonlatind)
rm(lonind, latind)
# ## check coordinates
# # #left top
# coords_ind_mat[1,]
# #right top
# coords_ind_mat[nc$dim$X$len,]
# #left bottom
# coords_ind_mat[(nc$dim$Y$len-1)*nc$dim$X$len+1,]
# #right bottom
# coords_ind_mat[nc$dim$X$len*nc$dim$Y$len,]
lonlat_mat = as.matrix(lonlat)

# # # create & write dataframe
sf_Prec_points = st_as_sf(coords_ind_mat, coords = c("x","y"), crs = 32633)
st_crs(sf_Prec_points)
plot(head(sf_Prec_points))
class(sf_Prec_points)

#write_sf(sf_Prec_points, dsn="Pre_nodes_KSS.gpkg")
shell.exec(getwd())
#create function to retrieve dataframe of indices of temp/prec points within a subbasin from the subbasin ID
indexlist = list()

# i = 195 is sybbasins ID 119
for (i in 1:nrow(subbasins)) {
  spdf_idx <- sf_Prec_points[subbasins[i,], ]
  if (nrow(spdf_idx)==0) {
    subid_center = st_centroid(subbasins[i,])
    nearest = st_nearest_feature(subid_center, sf_Prec_points)
    nearest = sf_Prec_points[nearest,]
    idx_list = as.matrix(st_drop_geometry(nearest))
  } else  {
    idx_list <- as.matrix(st_drop_geometry(spdf_idx)) 
  }
  colnames(idx_list) <- c("lon", "lat")
  indexlist[[i]] = idx_list
  print(i)
}


#check if all subbasins have at least one prec/temp point within them, also assign closest point
#if ever used when a raster has a very coarse resolution compared to the catchment size, an interpolation scheme should be implemented here
# #get list entries with no lat/lon entries or no gridnodes within subbasin polygon respectively
# cond_small <- sapply(indexlist, function(x) length(x) < 1)
# small_subID <- indexlist[cond_small]
# small_subID <- names(small_subID)
# rm(cond_small)

# #manual fixing of subcatchments that do not contain any precipitaton grid points
# #check in GIS system first if small subcatchments are meaningful to have (maybe only sliver polygons due to cutting)
# #add entries such as "indexlist$`119`" according to subbasinIDs in vector "small_subID"
# this is poorly hardcoded for the example of two points replacing the empty indexlist of subID 119
# # #1 (example)
# indexlist$`119` <- indexlist$`1`
# indexlist$`119`[1,] <- c(70,28)
# indexlist$`119`[2,] <- c(70,29)
# indexlist$`119` <- indexlist$`119`[1:2,]
# rm(small_subID)
#function to get mean of each subcatchment for given timestep by position in indexlist

fn_get_Prec_avg <- function(subbasinindex){
  avg_Prec_idx <- round(mean(Prec_ts[indexlist[[subbasinindex]]]), digits = 2)
  return(avg_Prec_idx )
}

#index vector from 1 to number of subbasins
index <- c(1:length(indexlist))

save.image("2023-06-11.RData")
#load("2023-05-26.RData")


getwd()
setwd(dir.nc)
shell.exec(getwd())
# #get total #iterations n
# tic()
# for (n_folder in 1:length(v.folder)) {
#   #loop through folders
#   setwd(paste0(dir.nc, v.folder[n_folder]))
#         for (j in 1:length(v.scenarios)) {
#           setwd(paste0(dir.nc,v.folder[n_folder],"/", v.scenarios[j]))
#           filelist =  list.files(pattern = "*.nc")
#           #filelist = filelist[str_detect(filelist,paste(scenario, collapse = "|"))]
#           for (filenumber in 1:length(filelist)) {
#             nc <- nc_open(filelist[filenumber], verbose = F, write = FALSE)
#             n = n + length(nc$dim$time$vals)
#             nc_close(nc)
#           }
#         }
# }
# #initialize progress bar
# pb <- txtProgressBar(min = 0, max = n, style = 3)
# toc()

#filelist = filelist[str_detect(filelist, paste(scenario, collapse = "|"))]

getwd()
setwd(dir.nc)
m=0
tic()
n_folder = 1
n_scenario = 1
for (n_folder in 1:length(v.folder)) {  #loop through model folders
  setwd(paste0(dir.nc, v.folder[n_folder]))
  print(getwd())
  for (n_scenario in 1:length(v.scenarios)) {  #loop through scenarios
    setwd(paste0(dir.nc, v.folder[n_folder], "/", v.scenarios[n_scenario]))
    print(getwd())
    filelist =  list.files(pattern = "*.nc", recursive = T)
    for (filenumber in 1:length(filelist)) {
      nc <- nc_open(filelist[filenumber], verbose = F, write = FALSE)
      if (filenumber == 1) {
        timeseries <- as.Date("1970-01-01") + nc$dim$time$vals
      }
      else {
        timeseries_temp <- as.Date("1970-01-01") + nc$dim$time$vals
        timeseries <- c(timeseries, timeseries_temp)
      }
      nc_close(nc)
    }
    for (filenumber in 1:length(filelist)) {
      # for (filenumber in 15:1){
      nc <- nc_open(filelist[filenumber], verbose = F, write = FALSE)
      ################AAAAA
      names(nc$var)[5]
      varname = names(nc$var)[5]
      print(paste("File opened: ",filelist[filenumber] ))
      #m = m + nc$dim$time$len
      for (ts in 1:nc$dim$time$len) {
        #loop through timesteps
        if (ts == 1 &&
            filenumber == 1) {
          #initialize df and timeseries
          Prec_ts <-
            ncvar_get(
              nc,
              varid = varname,
              start = c(1, 1, ts),
              count = c(-1,-1, 1),
              verbose = FALSE,
              signedbyte = TRUE,
              collapse_degen = TRUE
            )
          Prec_ts = Prec_ts * 0.1
          df_Prec <- sapply(index, fn_get_Prec_avg)
          df_Prec <- as.data.frame(t(df_Prec))
          colnames(df_Prec) <- names(indexlist)
          m = nc$dim$time$len
          
        }
        else {
          #extract data frome timestep and add to dataframe
          Prec_ts <-
            ncvar_get(
              nc,
              varid = varname,
              start = c(1, 1, ts),
              count = c(-1,-1, 1),
              verbose = FALSE,
              signedbyte = TRUE,
              collapse_degen = TRUE
            )
          Prec_ts = Prec_ts * 0.1
          new_row <- sapply(index, fn_get_Prec_avg)
          df_Prec <- rbind(df_Prec, new_row)
          
          if (ts == nc$dim$time$len &&
              filenumber == length(filelist)) {
            #end condition, write dataframe to file
            df_Prec <- cbind(as.character(timeseries), df_Prec)
            colnames(df_Prec)[1] <- "SUBID"
            setwd(dir.gis)
            filename = paste0("POBS_",
                              v.scenarios[n_scenario],
                              str_remove(v.folder[n_folder], "./"),
                              "_Gaula_247.txt")
            write.table(
              df_Prec,
              file = filename,
              quote = F,
              row.names = F,
              sep = "\t"
            )
            print(paste("File written: ", filename))
            
          }
        }
      }
      nc_close(nc)
    }
  }
}

toc()
