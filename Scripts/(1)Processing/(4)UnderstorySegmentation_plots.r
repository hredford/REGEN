########################################################################################################
########################################################################################################
#
#
#
#                               Understory Vegetation Segmentation
#                                     By: Hannah Redford
#
#
#
#
# The purpose of this script is to clip out the understory under a selected height value. 
#
# Inputs:
#   - ALS LAZ files
#
# Outputs: 
#   - Understory point cloud of all the plots 
#
########################################################################################################
########################################################################################################




{
  library(lidRplugins)
  library(terra)
  library(dplyr)
  library(lidR)
  library(sp)
  library(rgdal)
  library(mapview)
  library(rgl)
  library(colorspace)
  library(TreeLS)
  library(nabor)
  library(sf)
  library(lubridate)
  library(spatialEco)
  library(tidyverse)
  library(data.table)
}

setwd("G:\\REGEN\\")

tiles = list.files("Plot_Level_Datasets\\NormalizedClipped", 
                   patter = ".las",
                   full.names = TRUE)
names = list.files("Plot_Level_Datasets\\NormalizedClipped", 
                   pattern = ".las", 
                   full.names = FALSE)

base.folder = "G:\\REGEN\\"

uFolderName = "Plot_Level_Products\\UnderstoryPointClouds"

dir.create(paste0(base.folder, uFolderName)) # Create a folder that contains the classified tiles

#--------------------------- USER INPUTS -----------------------------#

# set maximum tree height - this will remove all TAOs above the set height AND their understory 
maxZ = 16L 

#---------------------------- 99th percentile Function -------------------------------#

z99 = function(z){
  q99 = quantile(z, probs = 0.99)
  return(q99)}

z95 = function(z){
  q95 = quantile(z, probs = 0.95)
  return(q95)}

###############################################################################################
#
#                        Loop through plots and get metrics
#
###############################################################################################

# unique_to_folder1 = paste0(unique_to_folder1, ".laz")
# names = names[basename(names) %in% unique_to_folder1]
# 
# tiles <- tiles[basename(tiles) %in% unique_to_folder1]


for (i in c(1:length(tiles))) {
  
  
  ###########################################################################
  #                         Pre-Processing
  ###########################################################################  
  
  #-------------------- initialize the necessary information ---------------------------# 
  name = tools::file_path_sans_ext(names[i])
  
  lazclip <- tryCatch(
    {
      readLAS(tiles[i])  # Attempt to read the .laz file
    },
    error = function(e) {
      # Print a warning message and return NULL
      message(paste("Error reading file:", name))
      message("Skipping this file due to error: ", e$message)
      return(NULL)
    }
  )
  
  if (length(lazclip$X) < 1) {
    next
  }
  
  # Proceed with processing if the file is successfully read
  message(paste("Processing file:", name))
  
  #plot(lazclip)
  
  # normalize the point cloud and re-clip to correct size
  lazclip <- normalize_height(lazclip, knnidw())
  
  #plot(lazclip)
  ###########################################################################
  #                    Tall Tree and Snag Segmentation
  ########################################################################### 
  
  #------------------------- Canopy Vegetation Metrics ----------------------------#
  
  # Canopy height model
  chm <- rasterize_canopy(lazclip, 0.5, pitfree(subcircle = 0.2)) 
  plot(chm)
  
  # locate trees and filter to trees only over the tallest tree value (16m)
  ttops <- locate_trees(lazclip, lmf(ws = 3))
  ttops = filter(ttops, Z > maxZ)
  
  #plot(chm)
  #plot(sf::st_geometry(ttops), add = TRUE, pch = 3)
  
  if (length(ttops$treeID) > 0) { # the ttops variable must have more than 1 tree for this next aprt
    
    #-------------- Segment Tree Columns --------------#
    print("Tree Found")
    
    # segmentation inputs
    ker <- matrix(1,3,3) # set moving window 
    filter1 <- filter_poi(lazclip, Z <= max(ttops$Z)) # filter the points to the height of the segmented trees
    vegrast<- grid_canopy(filter1, res = .3, p2r(0.5)) # new higher resolution CHM
    vegrast <- raster::focal(vegrast, w = ker, fun = mean, na.rm = TRUE) # smooth the canopy height model
    #plot(vegrast)
    
    # Segment the trees using the high-resolution CHM
    # this likely oversegments the tall trees, we don't care, we just want them gone
    veght2 <- tryCatch(segment_trees(filter1, silva2016(vegrast, ttops)))
    #plot(veght2, color = "treeID")
    
    
    # Assign classification to the segmented trees
    veght2@data[treeID > 0, Classification := 21]
    
    # filter the point cloud for only the segmented trees
    all.stems = filter_poi(veght2, Classification == 21L)
    #plot(all.stems)
    
    # filter the point cloud for the non-tree vegetation
    nostems = filter_poi(veght2, Classification < 20)
    nostems = filter_poi(nostems, Z < maxZ)
    #plot(nostems, bg = "white")
    
    #plot(nostems)
    #rgl::axes3d(col='white')
    
    
    #----------------- 99th percentile of understory at 20m  -------------------#
    
    p99 = grid_metrics(nostems, ~z99(Z), res = 20) # rasterize the height of the understory canopy layer 
    #plot(p99)
    
    veght2 = merge_spatial(veght2, p99, "p99") # add p99 as an attribute to the point cloud
    #plot(veght2, color = "treeID")
    
    
    
    #---------- Identify Tall Trees and remove them from the point cloud ------------# 
    
    # Extract only the points with a treeID
    trees <- veght2@data[!is.na(treeID)]  
    
    
    # Calculate 95th percentile for Z for each treeID (using only points 
    # above p99 without trees which we calculated in the previous step)
    q_vals <- trees[Z < p99, list(q = quantile(Z, 0.95, na.rm = TRUE)), by = treeID]
    
    # merge the qvals with the trees data table
    trees <- merge(trees, q_vals, by = "treeID", all.x = TRUE)
    
    # Classify points as tree or non-tree points (only for tree points)
    trees[, Classification := ifelse(Z > q, 21, 31)]
    
    # 7. Merge the updated tree classifications back into the full point cloud
    veght2@data <- merge(veght2@data, 
                         trees[, .(X, Y, Z, Classification, treeID)], 
                         by = c("X", "Y", "Z", "treeID"), 
                         all.x = TRUE)
    
    # SInce the merge step split the classificaitons, merge those by using the 
    # new classification attribute (Classification.y) as the primary layer
    veght2@data[, Classification := coalesce(Classification.y, Classification.x)]
    
    # filter the points to create a NOTREE layer
    notree = filter_poi(veght2, Classification != 21)
    #plot(notree)
    
    
    notree@data[treeID > 0, treeID := 0] # remove tree IDs
    
  } else{
    notree = lazclip
  }
  
  
  #----------------- Filter any outlying points that got missed ----------------#
  
  notree = nnFilter(notree, d = 2, n = 5)
  #plot(notree)
  
  
  ########################### Snag Removal ####################################
  print("Snag Removal")
  
  #-------------- Segment TAOs from the understory point cloud -----------------#
  # Locate the trees
  ttops2 <- locate_trees(notree, lmf(ws = 2))
  
  # make sure there are no points over 16m
  filter1 <- filter_poi(notree, Z <= maxZ)
  
  # create a CHM
  vegrast<- grid_canopy(filter1, res = .3, p2r(0.5))
  #plot(vegrast)
  
  # Segment trees
  veght2 <- tryCatch(segment_trees(filter1, silva2016(vegrast, ttops2)))
  #plot(veght2, color = "treeID")
  veght2@data[treeID > 0, Classification := 23] # classify the trees
  segmentedTAOs = filter_poi(veght2, Classification == 23) # create a point cloud of only the TAOS
  #plot(segmentedTAOs)
  
  # List tree IDs
  treeIds = unique(veght2$treeID)
  treeIds = na.omit(treeIds)
  
  #--------------- Run the Height to Width Ratio for Each TAO ------------------#
  
  
  # Extract treeIDs
  valid_trees <- veght2@data[Z > 0.5, .SD[.N > 5], by = treeID]
  
  # Calculate height, width, and ratio for each treeID
  ratios <- valid_trees[, {
    maxh <- max(Z, na.rm = TRUE)
    widthx <- max(X) - min(X)
    widthy <- max(Y) - min(Y)
    width <- mean(c(widthx, widthy))
    ratio <- width / maxh
    .(ratio = ratio, maxh = maxh)
  }, by = treeID]
  
  # Identify trees with ratio < 0.25
  treeIDs_low_ratio <- ratios[ratio < 0.25, treeID]
  
  # Update Classification for points belonging to identified trees
  veght2@data[treeID %in% treeIDs_low_ratio & Z > 1, Classification := as.integer(30)]
  veght2@data$Classification <- as.integer(veght2@data$Classification)
  
  nosnags = filter_poi(veght2, Classification != 30) 
  
  writeLAS(nosnags, paste0(base.folder, uFolderName,"\\", name, "_UnderstoyOnly.laz"))
  
  #plot(segmentedTAOs)
  rm(lazclip)
  rm(nosnags)
  rm(veght2)
  rm(notree)
  rm(nostems)
  rm(segmentedTAOs)
  rm(all.stems)
  rm(vegrast)
  rm(valid_trees)
  gc()
  
}

