########################################################################################################
########################################################################################################
#
#
#
#                                        ALS Clipping
#                                     By: Hannah Redford
#
#
#
#
# The purpose of this script is to clip the plots according to precise x/y locations. 
# The way these are found is at first through sub meter x/y coordinates from a 
# JAVAD GNSS receiver. These are not always totally accurate - so the TLS point
# clouds are registered to the ALS point cloud and the center coordinate is taken 
# to clip again. 
# 
#
#
# Inputs:
#   - Base data frame of the plot name, x/y coordinates, and dependent variables of oak counts, conifer counts, and shrub percent cover
#   - LAZ catalog data set of the 2019 ALS Yosemite Acquisition
#
# Outputs: 
#   - LAZ files of the plots
#
#
########################################################################################################
########################################################################################################
  



{ library(raster)
  library(dplyr)
  library(lidR)
  library(sp)
  library(ggplot2)
  library(rgdal)
  library(mapview)
  library(gridExtra)
  library(rgl)
  library(ggplot2)
  library(png)
  library(grid)
  library(sf)
}

  # Path to Lidar data
  YoseNP.laz <- "G:\\GDrive_Backup_NOYOSE_08042025\\YOSE_2019_ALS" # Load / Read LAZ files
  sarr.laz = "F:\\ssarr_raw"

  outputfilepath = "G:\\REGEN\\Plot_Level_Datasets\\PlotPatches_NormClipped"
  patches = st_read("G:\\REGEN\\Park_Level_Datasets\\PlotPatches\\PlotPatches.shp")
  patches <- st_zm(patches, drop = TRUE, what = "ZM")
  
  
  ################################################################################
  laz <- readLAScatalog(c(YoseNP.laz,sarr.laz))
  plot(laz, mapview = TRUE)
  
  for (i in nrow(patches)) {

    # Point Cloud Image
    lazclip <- clip_roi(laz, patches[i,])
    
    #plot(lazclip)
    #plot(lazclip, color = "Classification")
    
    #lazclip <- normalize_height(lazclip, knnidw())
    
    #plot(lazclip)

    writeLAS(lazclip, sprintf("%s\\%s.laz", outputfilepath, i)) # Write the laz file
    
    rm(lazclip)
    
  }
 

