########################################################################################################
########################################################################################################
#
#
#
#                              Building Understory Only Point Clouds
#                                     By: Hannah Redford
#
#
#
#
# The purpose of this script is to remove all tall trees while preserving the understory. 
#
# 1. It starts by normalizing and filtering the points clouds, then identifies tall trees using
#    a Local Maximum Filter (LMF) algorithm and Filter (LMF) algorithm and segements them
#    using silva2016 (a seed + voronoi tessellation method). 
#
# 2. Then, the 99th percentile height is calculated using grid_metrics at a 20m resolution.
#    The output for this step is a raster. 
#
# 3. Next, the script clips out each segmented TAO and the understory vegetation. 
#    The point cloud is clipped at the height of the surrounding vegetation from the 
#    previous step. Then, it is clipped once more at the 95th percentile of the 
#    remaining height. 
#
# 4. Finally, the remaining TAOs are identified. These will be smaller, and can 
#    represent anything from regeneration and shrubs to small snags and parts of
#    tall trees that didn't get clipped from the previous step. The width and height 
#    of these TAos is calculated, and the ratio (0.25) from these measurements determines
#    whether or not to keep them in the point cloud.  
#
# The resulting point clouds will be input through other scripts to calculate metrics. 
#
#
# Inputs:
#   - Understory ONLY Laz files
#   - 
#
# Outputs: 
#   - Metric raster stacks for each plot
#
#
########################################################################################################
########################################################################################################

library(lidRplugins)

{
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
  library(pbapply)
}

################################################################################
#                         INPUTS
################################################################################


input_folder = "G:\\REGEN\\Plot_Level_Datasets\\"
outputfilepath = "G:\\REGEN\\Plot_Level_Products\\"

plot.csv = read.csv(paste0(input_folder, "BaseDF_10-18_USEABLE.csv"))
als.path = paste0(outputfilepath, "UnderstoryPointClouds")

output_folder_name = "UnderstoryMetrics_101725"
output_folder = paste0(outputfilepath, output_folder_name)
dir.create(output_folder)


heightstrata <- c("55to80",
                  "80to90",
                  "90to99")

plotsize <- 10

################################################################################
# Functions
################################################################################

################################################################################
# For Loop
################################################################################

# SAFER PARALLEL APPROACH FOR WINDOWS
# Uses parallel package with explicit cluster management

library(parallel)
library(data.table)

# Detect available cores (leave one free)
n_cores <- max(1, detectCores() - 1)
n_cores <- min(n_cores, 5)  # Limit to 4 to avoid memory issues

cat(sprintf("Using %d cores for parallel processing\n", n_cores))

# Define the processing function for a single plot
process_plot <- function(i, plot.csv, als.path, heightstrata, output_folder) {
  
  # Load required packages within each worker
  suppressPackageStartupMessages({
    library(lidR)
    library(terra)
    library(raster)
    library(data.table)
  })
  
  row <- plot.csv[i,]
  plot.name <- row$plot
  message(paste("Processing file:", plot.name))
  
  tryCatch({
    lazclip <- list.files(als.path, pattern = plot.name, full.names = TRUE)
    lazclip <- readLAS(lazclip)
    
    # Filter snags and normalize if needed
    lazclip <- TreeLS::nnFilter(lazclip, d = 1, n = 10)
    
    # Create output raster list
    metric_rasters <- list()
    
    ref_raster <- rast(grid_metrics(lazclip, ~quantile(Z, probs = 0.25), 20))
    
    print("Height Metrics")
    
    for (h in seq_along(heightstrata)) {
      high <- as.numeric(sub(".*o", "", heightstrata[h]))
      low <- as.numeric(sub("t.*", "", heightstrata[h]))
      
      filter_noise <- function(las, low_val, high_val) {
        # Store thresholds as attributes so formulas can access them
        las <- add_attribute(las, low_val, "low_thresh")
        las <- add_attribute(las, high_val, "high_thresh")
        
        lowPC <- grid_metrics(las, ~quantile(Z, probs = low_thresh[1] / 100), 20)
        highPC <- grid_metrics(las, ~quantile(Z, probs = high_thresh[1] / 100), 20)
        las <- merge_spatial(las, lowPC, "low")
        las <- merge_spatial(las, highPC, "high")
        las <- filter_poi(las, Z > low)
        las <- filter_poi(las, Z < high)
        return(las)
      }
      
      hclip <- filter_noise(lazclip, low, high)
      
      if (length(hclip@data$X) == 0) {
        next  # Skip to next iteration
      } else {
        # Create a temporary column for the threshold comparison
        # This avoids scope issues with formulas in parallel processing
        hclip <- add_attribute(hclip, low, "low_threshold")
        p25 <- rast(grid_metrics(hclip, ~quantile(Z, probs = 0.25), 20))
        p25 <- extend(p25, ref_raster)
        names(p25) <- paste0("p25_", heightstrata[h])
        
        p75 <- rast(grid_metrics(hclip, ~quantile(Z, probs = 0.75), 20))
        p75 <- extend(p75, ref_raster)
        names(p75) <- paste0("p75_", heightstrata[h])
        
        p95 <- rast(grid_metrics(hclip, ~quantile(Z, probs = 0.95), 20))
        p95 <- extend(p95, ref_raster)
        names(p95) <- paste0("p95_", heightstrata[h])
        
        iqr <- p75 - p25
        names(iqr) <- paste0("iqr_", heightstrata[h])
        
        iqr95 <- p95 - p25
        names(iqr95) <- paste0("iqr95_", heightstrata[h])
        
        propor <- rast(grid_metrics(hclip, ~mean(Z > low_threshold[1]), 20))
        propor <- extend(propor, ref_raster)
        names(propor) <- paste0("propor_", heightstrata[h])
        
        hrumple <- rast(grid_metrics(hclip, ~rumple_index(X,Y,Z), 20))
        hrumple <- extend(hrumple, ref_raster)
        names(hrumple) <- paste0("rumple_", heightstrata[h])
        
        hmeanz <- rast(grid_metrics(hclip, ~mean(Z), 20))
        hmeanz <- extend(hmeanz, ref_raster)
        names(hmeanz) <- paste0("zmean_", heightstrata[h])
        
        hsdz <- rast(grid_metrics(hclip, ~sd(Z), 20))
        hsdz <- extend(hsdz, ref_raster)
        names(hsdz) <- paste0("zsd_", heightstrata[h])
        
        httops <- locate_trees(hclip, lmf(ws = 2))
        htaos_grid <- rast(ext(ref_raster), resolution = 20, crs = crs(ref_raster))
        httops_sv <- vect(httops)
        htaos <- rasterize(httops_sv, htaos_grid, field = 1, fun = "sum", background = 0)
        htaos <- resample(htaos, ref_raster, method = "near")
        
        if (length(hclip@data$X) > 21) {
          
          eigenMetrics <- point_eigenvalues(hclip, k = 20, metrics = TRUE)
          
          # Add the metrics to the point cloud
          hclip <- add_attribute(hclip, eigenMetrics$eigen_largest, "emax")
          hclip <- add_attribute(hclip, eigenMetrics$eigen_medium, "emid")
          hclip <- add_attribute(hclip, eigenMetrics$eigen_smallest, "emin")
          hclip <- add_attribute(hclip, eigenMetrics$curvature, "curv")
          hclip <- add_attribute(hclip, eigenMetrics$linearity, "lin")
          hclip <- add_attribute(hclip, eigenMetrics$planarity, "plan")
          hclip <- add_attribute(hclip, eigenMetrics$sphericity, "sphere")
          hclip <- add_attribute(hclip, eigenMetrics$anisotropy, "anisotr")
          
          # extract metrics to rasters
          emax <- rast(grid_metrics(hclip, ~mean(emax, na.rm = TRUE), 20))
          emid <- rast(grid_metrics(hclip, ~mean(emid, na.rm = TRUE), 20))
          emin <- rast(grid_metrics(hclip, ~mean(emin, na.rm = TRUE), 20))
          curvature <- rast(grid_metrics(hclip, ~mean(curv, na.rm = TRUE), 20))
          linearity <- rast(grid_metrics(hclip, ~mean(lin, na.rm = TRUE), 20))
          planarity <- rast(grid_metrics(hclip, ~mean(plan, na.rm = TRUE), 20))
          sphericity <- rast(grid_metrics(hclip, ~mean(sphere, na.rm = TRUE), 20))
          anisotropy <- rast(grid_metrics(hclip, ~mean(anisotr, na.rm = TRUE), 20))
          
          # Change the extent
          emax <- extend(emax, ref_raster)
          emid <- extend(emid, ref_raster)
          emin <- extend(emin, ref_raster)
          curvature <- extend(curvature, ref_raster)
          linearity <- extend(linearity, ref_raster)
          planarity <- extend(planarity, ref_raster)
          sphericity <- extend(sphericity, ref_raster)
          anisotropy <- extend(anisotropy, ref_raster)
          
          # Change the names
          names(emax) <- paste0("emax_", heightstrata[h])
          names(emid) <- paste0("emid_", heightstrata[h])
          names(emin) <- paste0("emin_", heightstrata[h])
          names(curvature) <- paste0("curv_", heightstrata[h])
          names(linearity) <- paste0("linear_", heightstrata[h])
          names(planarity) <- paste0("planar_", heightstrata[h])
          names(sphericity) <- paste0("sphere_", heightstrata[h])
          names(anisotropy) <- paste0("anisotr_", heightstrata[h])
          names(htaos) <- paste0("taos_", heightstrata[h])
        }
        
        # Add these rasters to the master list immediately
        metric_rasters <- c(metric_rasters, list(
          p25, p75, p95, iqr, iqr95, propor, hrumple, hmeanz, hsdz, htaos,
          emax, emid, emin, curvature, linearity, planarity, sphericity, anisotropy
        ))
      }
    }
    
    ############################ TAO Metrics #####################################
    
    print("TAO Metric Calculation")
    chm <- grid_canopy(lazclip, res = 2, p2r())
    ker <- matrix(1,3,3)
    vegrast <- raster::focal(chm, w = ker, fun = "mean", na.rm = TRUE)
    
    vegrast <- grid_canopy(lazclip, res = .3, p2r(0.5))
    ttops <- locate_trees(lazclip, lmf(ws = 2))
    vegrast <- raster::focal(vegrast, w = ker, fun = mean, na.rm = TRUE)
    veght2 <- tryCatch(segment_trees(lazclip, silva2016(vegrast, ttops)))
    
    veght2@data[treeID > 0, Classification := 24]
    
    if (all(is.na(veght2@data$treeID))) {
      
      ref_raster <- metric_rasters[[1]]
      r0 <- rast(ref_raster)
      r0 <- setValues(r0, 0)
      
      emax_TAO <- r0
      emid_TAO <- r0
      emin_TAO <- r0
      curv_TAO <- r0
      lin_TAO <- r0
      plan_TAO <- r0
      spher_TAO <- r0
      aniso_TAO <- r0
      
      names(emax_TAO) <- "emax_TAO"
      names(emid_TAO) <- "emid_TAO"
      names(emin_TAO) <- "emin_TAO"
      names(curv_TAO) <- "curv_TAO"
      names(lin_TAO) <- "lin_TAO"
      names(plan_TAO) <- "plan_TAO"
      names(spher_TAO) <- "spher_TAO"
      names(aniso_TAO) <- "aniso_TAO"
      
      r0 <- rast(ref_raster)
      r0 <- setValues(r0, 0)
      totalNtaos <- r0
      names(totalNtaos) <- "totalNtaos"
      values(totalNtaos) <- length(ttops$treeID)
      
      metric_rasters <- c(metric_rasters, list(
        emax_TAO, emid_TAO, emin_TAO, curv_TAO,
        lin_TAO, plan_TAO, spher_TAO, aniso_TAO, totalNtaos
      ))
      
    } else {
      
      segmentedTAOs <- filter_poi(veght2, Classification == 24)
      k_eig <- 21L
      
      grid20 <- rast(ref_raster)
      
      npt <- npoints(segmentedTAOs)
      pid_vals <- as.integer(seq_len(npt))
      segmentedTAOs <- add_attribute(segmentedTAOs, pid_vals, "pid")
      
      tree_ids <- sort(unique(segmentedTAOs@data$treeID))
      chunks <- lapply(tree_ids, function(id) filter_poi(segmentedTAOs, treeID == id))
      names(chunks) <- as.character(tree_ids)
      
      k_eig <- as.integer(k_eig)
      eig_fun <- function(las_i) {
        np <- npoints(las_i)
        if (np < 2L) return(NULL)
        k_use <- k_eig
        if (np - 1L < k_use) k_use <- max(1L, np - 1L)
        e <- point_eigenvalues(las_i, k = k_use, metrics = TRUE)
        data.frame(
          pid    = las_i@data$pid,
          treeID = las_i@data$treeID[1],
          emax   = e$eigen_largest,
          emid   = e$eigen_medium,
          emin   = e$eigen_smallest,
          curv   = e$curvature,
          lin    = e$linearity,
          plan   = e$planarity,
          spher  = e$sphericity,
          aniso  = e$anisotropy
        )
      }
      res_list <- lapply(chunks, eig_fun)
      res_list <- res_list[!vapply(res_list, is.null, TRUE)]
      res_all  <- rbindlist(res_list, use.names = TRUE)
      
      xy <- cbind(segmentedTAOs@data$X, segmentedTAOs@data$Y)
      cell_idx <- cellFromXY(grid20, xy)
      
      dt <- data.table(
        pid  = segmentedTAOs@data$pid,
        cell = cell_idx
      )
      setkey(dt, pid)
      
      setDT(res_all)
      setkey(res_all, pid)
      dt <- res_all[dt, nomatch = 0L]
      dt <- dt[!is.na(cell)]
      
      per_cell_tree <- dt[, .(
        emax  = mean(emax,  na.rm = TRUE),
        emid  = mean(emid,  na.rm = TRUE),
        emin  = mean(emin,  na.rm = TRUE),
        curv  = mean(curv,  na.rm = TRUE),
        lin   = mean(lin,   na.rm = TRUE),
        plan  = mean(plan,  na.rm = TRUE),
        spher = mean(spher, na.rm = TRUE),
        aniso = mean(aniso, na.rm = TRUE)
      ), by = .(cell, treeID)]
      
      per_cell_sum <- per_cell_tree[, .(
        emax  = sum(emax,  na.rm = TRUE),
        emid  = sum(emid,  na.rm = TRUE),
        emin  = sum(emin,  na.rm = TRUE),
        curv  = sum(curv,  na.rm = TRUE),
        lin   = sum(lin,   na.rm = TRUE),
        plan  = sum(plan,  na.rm = TRUE),
        spher = sum(spher, na.rm = TRUE),
        aniso = sum(aniso, na.rm = TRUE)
      ), by = cell]
      
      mk_layer <- function(cells, vals, name) {
        r <- rast(ref_raster)
        r <- setValues(r, rep(0, ncell(r)))
        if (length(cells)) r[cells] <- vals
        names(r) <- name
        r
      }
      
      emax_TAO  <- mk_layer(per_cell_sum$cell, per_cell_sum$emax,  "emax_TAO")
      emid_TAO  <- mk_layer(per_cell_sum$cell, per_cell_sum$emid,  "emid_TAO")
      emin_TAO  <- mk_layer(per_cell_sum$cell, per_cell_sum$emin,  "emin_TAO")
      curv_TAO  <- mk_layer(per_cell_sum$cell, per_cell_sum$curv,  "curv_TAO")
      lin_TAO   <- mk_layer(per_cell_sum$cell, per_cell_sum$lin,   "lin_TAO")
      plan_TAO  <- mk_layer(per_cell_sum$cell, per_cell_sum$plan,  "plan_TAO")
      spher_TAO <- mk_layer(per_cell_sum$cell, per_cell_sum$spher, "spher_TAO")
      aniso_TAO <- mk_layer(per_cell_sum$cell, per_cell_sum$aniso, "aniso_TAO")
      
      r0 <- rast(ref_raster)
      r0 <- setValues(r0, 0)
      totalNtaos <- r0
      names(totalNtaos) <- "totalNtaos"
      values(totalNtaos) <- length(ttops$treeID)
      
      metric_rasters <- c(metric_rasters, list(
        emax_TAO, emid_TAO, emin_TAO, curv_TAO,
        lin_TAO, plan_TAO, spher_TAO, aniso_TAO, totalNtaos
      ))
      
      align_to_ref <- function(r, ref) {
        r <- extend(r, ref) 
        r
      }
      
      metric_rasters <- lapply(metric_rasters, align_to_ref, ref = ref_raster)
      
      ok <- vapply(metric_rasters, function(x) compareGeom(x, ref_raster, stopOnError = FALSE), TRUE)
      stopifnot(all(ok))
    }
    
    # Stack and write out
    if (length(metric_rasters) > 0) {
      rstack <- rast(metric_rasters)
      raster_out_dir <- file.path(output_folder)
      out_path <- file.path(raster_out_dir, paste0(plot.name, "_stack.tif"))
      writeRaster(rstack, out_path, overwrite = TRUE)
    }
    
    # Force garbage collection
    rm(lazclip, metric_rasters, rstack)
    gc()
    
    return(list(success = TRUE, plot = plot.name))
    
  }, error = function(e) {
    return(list(success = FALSE, plot = plot.name, error = as.character(e)))
  })
}

# Create cluster explicitly
cl <- makeCluster(n_cores, type = "PSOCK")

# Export necessary variables to cluster
clusterExport(cl, c("plot.csv", "als.path", "heightstrata", "output_folder"))

# Set RNG for reproducibility
clusterSetRNGStream(cl, 123)

# Run parallel processing with error handling
cat("Starting parallel processing...\n")
results <- parLapply(cl, 1:nrow(plot.csv), process_plot,
                     plot.csv = plot.csv,
                     als.path = als.path,
                     heightstrata = heightstrata,
                     output_folder = output_folder)

# Stop cluster
stopCluster(cl)

# Check results
successful <- sapply(results, function(x) x$success)
cat(sprintf("\nProcessing complete: %d/%d plots successful\n", 
            sum(successful), length(successful)))

if (any(!successful)) {
  cat("\nFailed plots:\n")
  failed <- results[!successful]
  for (f in failed) {
    cat(sprintf("  - %s: %s\n", f$plot, f$error))
  }
}





