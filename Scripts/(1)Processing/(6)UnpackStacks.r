# Unpacking the Understory Rasters 

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
  library(stringr)
  library(tools)
}

input_folder = "G:\\REGEN\\"

plot.csv = read.csv(paste0(input_folder, "Plot_Level_Datasets\\BaseDF_10-18_USEABLE.csv"))
stack.path = list.files(paste0(input_folder, "Plot_Level_Products\\UnderstoryMetrics_101725"), full.names = TRUE, pattern = ".tif$")
shp = vect(paste0(input_folder, "Plot_Level_Datasets\\Shapefile\\BaseDF_shp.shp"))

seasonDiff = list.files("E:\\ARCHIVE_03-17-25\\YOSE_Regen_Sampling\\Datasets\\Sentinel_Image", pattern = "seasonDiff-", full.names = TRUE)
seasonDiffs = list.files("E:\\ARCHIVE_03-17-25\\YOSE_Regen_Sampling\\Datasets\\Sentinel_Image", pattern = "Sum", full.names = TRUE)
seasonDiffw = list.files("E:\\ARCHIVE_03-17-25\\YOSE_Regen_Sampling\\Datasets\\Sentinel_Image", pattern = "Win", full.names = TRUE)

yose.naip = rast("E:\\ARCHIVE_03-17-25\\YOSE_Regen_Sampling\\Datasets\\NAIP\\yose_2020_utm11.tif")

outputfilepath = "G:\\REGEN\\Plot_Level_Datasets\\"
outputfilename = "responsematrix_10-17-25_2.csv"

################################################################################

# seasonDiff1 = rast(seasonDiff[[1]])
# seasonDiff1 = terra::project(seasonDiff1, shp)
# 
# seasonDiffs = rast(seasonDiffs[[1]])
# seasonDiffs = terra::project(seasonDiffs, shp)
# 
# seasonDiffw = rast(seasonDiffw[[1]])
# seasonDiffw = terra::project(seasonDiffw, shp)

################################################################################


stack_paths <- stack.path

# Optional: plot IDs from filenames (or replace with your own 'plot_names' vector of length 121)
plot_ids <- basename(file_path_sans_ext(stack_paths))

# 1) Pass 1 — compute per-plot, per-layer means (as named numeric vectors)
means_list <- vector("list", length(stack_paths))
names(means_list) <- plot_ids

for (i in seq_along(stack_paths)) {
  pth <- stack_paths[i]
  plot_id <- plot_ids[i]
  
  # read one plot stack
  x <- try(rast(pth), silent = TRUE)
  if (inherits(x, "try-error")) {
    warning("Failed to read: ", pth)
    means_list[[i]] <- numeric(0)
    next
  }
  
  # get means per layer (global works chunked; memory-safe)
  gm <- try(global(x, "mean", na.rm = TRUE), silent = TRUE)
  if (inherits(gm, "try-error") || is.null(gm) || !nrow(gm)) {
    warning("Failed to compute means for: ", pth)
    means_list[[i]] <- numeric(0)
    next
  }
  
  # turn into a named numeric vector: names = layer names, values = means
  rn <- rownames(gm)
  v  <- gm$mean
  names(v) <- rn
  means_list[[i]] <- v
}

# 2) Build the union of all metric (layer) names across plots
all_metrics <- sort(unique(unlist(lapply(means_list, names))))

# 3) Allocate result matrix and fill (rows = plots, cols = metrics)
res_mat <- matrix(NA_real_, nrow = length(stack_paths), ncol = length(all_metrics))
rownames(res_mat) <- plot_ids
colnames(res_mat) <- all_metrics

for (i in seq_along(means_list)) {
  v <- means_list[[i]]
  if (!length(v)) next
  res_mat[i, names(v)] <- v
}

# 4) Final data.frame (plot id column + metric columns)
all_data <- data.frame(plot = rownames(res_mat), res_mat, row.names = NULL, check.names = FALSE)
all_data$plot <- sub("_stack$", "", all_data$plot)

################################################################################
# Add NAIP and Sentinel
################################################################################
# 
plotsh = all_data
plotsh = merge(plot.csv, plotsh, by = "plot")
# 
# plotsh$red = 0
# plotsh$green = 0
# plotsh$blue = 0
# plotsh$nir = 0
# 
# for (i in c(1:length(plotsh$y))) {
#   
#   row = plotsh[i,]
#   plotID = row$plot
#   x = row$x
#   y = row$y
#   
#   shp_select = shp[shp$plot == plotID,]
#   
#   buff = buffer(shp_select, 10)
#   buff = terra::project(buff, yose.naip)
#   #yose.naip = extend(yose.naip, buff)
#   
#   clip.naip = crop(yose.naip, buff)
#   #plot(clip.naip)
#   
#   plotsh[i,]$red = mean(values(clip.naip$yose_2020_utm11_1))
#   plotsh[i,]$green = mean(values(clip.naip$yose_2020_utm11_2))
#   plotsh[i,]$blue = mean(values(clip.naip$yose_2020_utm11_3))
#   plotsh[i,]$nir = mean(values(clip.naip$yose_2020_utm11_4))
#   
# }

shp <- vect(plotsh, geom = c("x", "y"),
            crs = crs(shp))

yosejunendvi = rast("G:\\REGEN\\Park_Level_Datasets\\Sentinel\\june_ndir.tif")
yosejunendir = rast("G:\\REGEN\\Park_Level_Datasets\\Sentinel\\june_ndvi.tif")

yoseoctndvi = rast("G:\\REGEN\\Park_Level_Datasets\\Sentinel\\oct_ndvi.tif")
yoseoctndir = rast("G:\\REGEN\\Park_Level_Datasets\\Sentinel\\oct_ndir.tif")

tkjunendvi = rast("G:\\REGEN\\Park_Level_Datasets\\Sentinel\\TK_juneNDVI.tif")
tkjunendir = rast("G:\\REGEN\\Park_Level_Datasets\\Sentinel\\TK_juneNDIR.tif")

tkoctndvi = rast("G:\\REGEN\\Park_Level_Datasets\\Sentinel\\TK_octNDVI.tif")
tkoctndir = rast("G:\\REGEN\\Park_Level_Datasets\\Sentinel\\TK_octNDIR.tif")


ext_buf5 <- function(r, pts) {
  terra::extract(r, terra::project(pts, r), buffer = 10, fun = mean, na.rm = TRUE, ID = FALSE)[[1]]
}


for (nm in c("junendvi","junendir","octndvi","octndir")) {
  if (!nm %in% names(shp)) shp[[nm]] <- NA_real_
}

n <- nrow(shp)
idx_tk   <- intersect(111:122, seq_len(n))
idx_yose <- setdiff(seq_len(n), idx_tk)

# ---- YOSE rows ----
if (length(idx_yose)) {
  shp$junendvi[idx_yose] <- ext_buf5(yosejunendvi, shp[idx_yose, ])
  shp$junendir[idx_yose] <- ext_buf5(yosejunendir, shp[idx_yose, ])
  shp$octndvi[idx_yose]  <- ext_buf5(yoseoctndvi, shp[idx_yose, ])
  shp$octndir[idx_yose]  <- ext_buf5(yoseoctndir, shp[idx_yose, ])
}

# ---- TK rows (118:129) ----
if (length(idx_tk)) {
  shp$junendvi[idx_tk] <- ext_buf5(tkjunendvi, shp[idx_tk, ])
  shp$junendir[idx_tk] <- ext_buf5(tkjunendir, shp[idx_tk, ])
  shp$octndvi[idx_tk]  <- ext_buf5(tkoctndvi, shp[idx_tk, ])
  shp$octndir[idx_tk]  <- ext_buf5(tkoctndir, shp[idx_tk, ])
}


df = as.data.frame(shp)
df$seasDiff_ndvi = df$junendvi-df$octndvi
df$seasDiff_ndir = df$junendir-df$octndir

write.csv(df, sprintf("%s\\%s", outputfilepath, outputfilename)) 












































