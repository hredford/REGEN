library(raster)
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
library(TreeLS)

aoi <- raster("D:\\SEFS502\\AOI_Raster.tif") # Load AOI File (UTM 11)
path = "E:\\YOSE_Regen_Sampling\\Datasets\\TLS_LAZ\\RawTLS\\2023_TLS"

tls = list.files("E:\\YOSE_Regen_Sampling\\Datasets\\TLS_LAZ\\RawTLS\\2023_TLS")
outputfilepath = "E:\\YOSE_Regen_Sampling\\AnalysisMaterials\\Sep4_run_CURRENT\\TLS_Normalized"


################################################################################

tls = tls[70:86]

for (i in tls) {
  
  lazclip = readTLS(paste0(path,'\\', i))
  lazclip <- tlsNormalize(lazclip, keep_ground = TRUE)
  lazclip = filter_poi(lazclip, Z > 0)
  
  #plot(lazclip)
  
  writeLAS(lazclip, sprintf("%s\\%s_TLSnorm_0910.las", outputfilepath, i)) # Write the laz file
  
  rm(lazclip)
  
}
