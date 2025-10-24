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
}

  
  plots <- read.csv("Plot_Level_Datasets\\BaseDF_10-18_USEABLE.csv") #Read the X/Y dataset
  head(plots) # find the x/y columns
  p.size <- 10 # set plot size to clip 
  
  # Path to Lidar data
  YoseNP.laz <- "G:\\GDrive_Backup_NOYOSE_08042025\\YOSE_2019_ALS" # Load / Read LAZ files
  sarr.laz = "F:\\ssarr_raw"

  outputfilepath = "Plot_Level_Datasets\\NormalizedClipped"
  
  
  ################################################################################
  laz <- readLAScatalog(YoseNP.laz)
  sarr.laz <- readLAScatalog(sarr.laz)
  plot(laz, mapview = TRUE)
  plot(sarr.laz, mapview = TRUE)

  plot.list <- plots$plot
  plot.list = plot.list[c(125:131)]
  
  for (i in plot.list) {
    
    print(i)
    plot <- filter(plots, plot == i) # Select the plot you want to clip
    x <- as.numeric(plot$x) # extract x coordinate
    y <- as.numeric(plot$y) # extract y coordinate
    
    # Point Cloud Image
    lazclip <- clip_circle(sarr.laz, xcenter = x, ycenter = y, radius = p.size)
    
    #plot(lazclip)
    plot(lazclip, color = "Classification")
    
    lazclip <- normalize_height(lazclip, knnidw())
    lazclip <- clip_circle(lazclip, xcenter = x, ycenter = y, radius = p.size)
    
    plot(lazclip)

    writeLAS(lazclip, sprintf("%s\\%s_Reclip_10m.las", outputfilepath, i)) # Write the laz file
    
    rm(lazclip)
    
  }
 

