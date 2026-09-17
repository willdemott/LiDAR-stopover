# Modeling dependency of stopover intensity on LiDAR-derived Forest Characteristics
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠓⠒⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⣀⠀⠀⠀⠀⠀⢠⢤⣤⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠔⠒⠒⠲⠎⠀⠀⢹⡃⢀⣀⠀⠑⠃⠀⠈⢀⠔⠒⢢⠀⠀⠀⡖⠉⠉⠉⠒⢤⡀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠔⠚⠙⠒⠒⠒⠤⡎⠀⠀⠀⠀⢀⣠⣴⣦⠀⠈⠘⣦⠑⠢⡀⠀⢰⠁⠀⠀⠀⠑⠰⠋⠁⠀⠀⠀⠀⠀⠈⢦⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⣸⠁⠀⠀⠀⠀⠀⠀⢰⠃⠀⣀⣀⡠⣞⣉⡀⡜⡟⣷⢟⠟⡀⣀⡸⠀⡎⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⣻⠀⠀⠀⠀
#⢰⠂⠀⠀⠀⠀⠀⠀⠀⣗⠀⠀⢀⣀⣀⣀⣀⣀⣓⡞⢽⡚⣑⣛⡇⢸⣷⠓⢻⣟⡿⠻⣝⢢⠀⢇⣀⡀⠀⠀⠀⢈⠗⠒⢶⣶⣶⡾⠋⠉⠀⠀⠀⠀⠀
#⠈⠉⠀⠀⠀⠀⠀⢀⠀⠈⠒⠊⠻⣷⣿⣚⡽⠃⠉⠀⠀⠙⠿⣌⠳⣼⡇⠀⣸⣟⡑⢄⠘⢸⢀⣾⠾⠥⣀⠤⠖⠁⠀⠀⠀⢸⡇⠀⠀⠀⠀⠀⢀⠀⠀
#⠀⠀⠀⢰⢆⠀⢀⠏⡇⠀⡀⠀⠀⠀⣿⠉⠀⠀⠀⠀⠀⠀⠀⠈⢧⣸⡇⢐⡟⠀⠙⢎⢣⣿⣾⡷⠊⠉⠙⠢⠀⠀⠀⠀⠀⢸⡇⢀⠀⠀⠀⠀⠈⠣⡀
#⠀⠀⠀⠘⡌⢣⣸⠀⣧⢺⢃⡤⢶⠆⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⣟⠋⢀⠔⣒⣚⡋⠉⣡⠔⠋⠉⢰⡤⣇⠀⠀⠀⠀⢸⡇⡇⠀⠀⠀⠀⠀⠀⠸
#⠀⠀⠀⠀⠑⢄⢹⡆⠁⠛⣁⠔⠁⠀⣿⠀⠀⠀⠀⢸⡇⠀⠀⠀⠀⠀⣿⢠⡷⠋⠁⠀⠈⣿⡇⠀⠀⠀⠈⡇⠉⠀⠀⠀⠀⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠑⣦⡔⠋⠁⠀⠀⠀⣿⠀⠀⢠⡀⢰⣼⡇⠀⡀⠀⠀⣿⠀⠁⠀⠀⠀⠀⣿⣷⠀⠀⠀⠀⡇⠀⠀⢴⣤⠀⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⢰⣿⡇⠀⠀⠀⠀⠀⣿⡀⠀⢨⣧⡿⠋⠀⠘⠛⠀⠀⣿⠀⠀⢀⠀⠀⠀⣿⣿⠀⠀⠀⠀⢲⠀⠀⠀⠀⠀⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⢸⣿⡇⠀⠀⠀⠀⠀⢸⡧⡄⠀⠹⣇⡆⠀⠀⠀⠀⠀⣿⠀⢰⣏⠀⣿⣸⣿⣿⠀⠀⠀⠀⣼⠀⠀⠰⠗⠀⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⢸⣿⡇⠀⠀⠀⠀⠀⢸⡇⣷⣛⣦⣿⢀⠈⠑⠀⢠⡆⣿⠐⢠⣟⠁⢸⠸⣿⣿⢱⣤⢀⠀⣼⠀⠀⢀⠀⠀⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⢸⣿⡇⠀⢀⠀⠀⠀⢸⡇⠘⠫⣟⡇⠊⣣⠘⠛⣾⡆⢿⠀⠙⣿⢀⣘⡃⣿⣿⡏⠉⠒⠂⡿⠀⠰⣾⡄⠀⢸⡟⣽⣀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠸⣿⡇⠀⠘⣾⠀⠀⢸⡇⢸⣇⡙⠣⠀⣹⣇⠀⠈⠧⢀⣀⣀⡏⣸⣿⣇⢹⣿⡇⢴⣴⣄⣀⡀⢰⣿⡇⠀⢸⣇⢿⡿⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠓⠁⠈⠻⢷⠾⠦⠤⠬⣅⣹⣿⣖⣶⣲⣈⡥⠤⠶⡖⠛⠒⠛⠁⠉⠛⠮⠐⢛⡓⠒⢛⠚⠒⠒⠒⠛⣚⣫⡼⠿⠿⣯⠛⠤⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠉⠉⠉⠉⡉⠉⠁⠀⠀⠘⠓⠀⠀⠀⠀⠀⣀⣞⡿⡉⠉⠉⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠹⣶⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀

###################################################################
# Load Packages
{
library(terra)
library(mgcv)
library(mapview)
library(tidyverse)
library(viridisLite)
library(marginaleffects)
library(patchwork)
}

##################################################################
# LiDAR DATA
{
p95_file <- rast("C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/Northeast_epsg6350_p95_10m_rounded.tif")

nlcd <- rast("C:/Users/willd/Desktop/Everything/GIS 5453/GIS_5453/NLCD_TIFS/Annual_NLCD_LndCov_2024_CU_C1V1.tif")

# crop NLCD in its own CRS then align to the 10 m LiDAR grid
foot <- project(as.polygons(ext(p95_file), crs = crs(p95_file)), crs(nlcd), filename = "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/foot.tif", overwrite = T)
nlcd <- crop(nlcd, foot)

# Forest-only
forest <- project(nlcd %in% c(41, 42, 43, 90), p95_file, method = "near", filename = "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/forest.tif", overwrite = T)

# p95 masked at 10 m and aggregated
p95_for <- mask(p95_file, forest, maskvalues = c(0, NA))

# Aggregate
p95      <- aggregate(p95_for, 30, fun = "mean", na.rm = T)
sdp95    <- aggregate(p95_for, 30, fun = "sd",   na.rm = T)
}
##################################################################
# RADAR DATA
{
box_dir <- "C:/Users/willd/Box/"
setwd(box_dir)

# Every station in the list
stations <- unique(read.csv("./DFF_resources/DFF_radar_BCR_list.csv")$radar)

spring_files <- sprintf("./DFF_ppi_years/%s_2013-24_%s_bioRad_ppi_years_AGL_inverse_Refit.tif", stations, "spring")
Test<-rast(spring_files[1])
autumn_files <- sprintf("./DFF_ppi_years/%s_2013-24_%s_bioRad_ppi_years_AGL_inverse_Refit.tif", stations, "autumn")

# VIR mosaic onto that grid; mean(rast(f)...) averages across all years, while fun = "mean" averages overlapping values in the mosaic
spring_mean <- mosaic(sprc(lapply(spring_files, \(f) mean(rast(f), na.rm = T))), fun = "mean")
autumn_mean <- mosaic(sprc(lapply(autumn_files, \(f) mean(rast(f), na.rm = T))), fun = "mean")

# 2024 only
spring_mean <- mosaic(sprc(lapply(spring_files, \(f) rast(f)[["Mean_2024"]])), fun = "mean")
autumn_mean <- mosaic(sprc(lapply(autumn_files, \(f) rast(f)[["Mean_2024"]])), fun = "mean")

# Align CRS, crop to LiDAR area
spring_vir <- project(spring_mean, p95, method = "bilinear")
autumn_vir <- project(autumn_mean, p95, method = "bilinear")
}
##################################################################
# Stack and turn into df for modeling
{
spring_stk <- c(spring_vir, p95, sdp95)
autumn_stk <- c(autumn_vir, p95, sdp95)

names(spring_stk) <- names(autumn_stk) <- c("vir", "p95", "sdp95")

# keep only cells where every layer is valid
spring_stk <- mask(spring_stk, noNA(spring_stk, falseNA = T))
autumn_stk <- mask(autumn_stk, noNA(autumn_stk, falseNA = T))

grn <- colorRampPalette(c("#c7e9c0", "#74c476", "#31a354", "#00441b"))

# Quick map to visualize our data spatially
vir_cap <- clamp(spring_stk[["vir"]], upper = 7000, values = FALSE)
mapview(vir_cap, layer.name = "Spring VIR",
        col.regions = viridis, na.color = "transparent") +
  mapview(spring_stk[["p95"]], layer.name = "Canopy height p95 (m)",
          col.regions = grn, na.color = "transparent") +
  mapview(spring_stk[["sdp95"]], layer.name = "Canopy height SD (m)",
          col.regions = viridis, na.color = "transparent")

spring_d <- as.data.frame(spring_stk, xy = T, na.rm = T)
autumn_d <- as.data.frame(autumn_stk, xy = T, na.rm = T)

nrow(spring_d) 
nrow(autumn_d)
}
##################################################################
# Create alternative stopover responses
{
thr <- global(spring_mean, \(v, ...) quantile(v, 0.90, na.rm = TRUE))[[1]]
  
spring_d$hot <- as.integer(spring_d$vir >= thr)
  
mean(spring_d$hot)
}
##################################################################
# Model
{
# Remove outliers
spring_d <- spring_d[spring_d$vir <= 7000, ]

spring_m_bam <- bam(vir ~ s(p95, k = 20) + s(sdp95, k = 20) + s(x, y, k = 800), family = Gamma(link = "log"), data = spring_d, method = "fREML", discrete = T)
summary(spring_m_bam)

spring_m_hot <- bam(hot ~ s(p95, k = 20) + s(sdp95, k = 20) + s(x, y, k = 800), family = binomial, data = spring_d, method = "fREML", discrete = T)
summary(spring_m_hot)
}
#################################################################
# Inspect Models

g1 <- plot_predictions(spring_m_bam, condition = "p95", type = "response") +
  labs(x = "p95 (m)", y = "Predicted VIR")

g2 <- plot_predictions(spring_m_bam, condition = "sdp95", type = "response") +
  labs(x = "sd p95 (m)", y = "Predicted VIR")

h1 <- plot_predictions(spring_m_hot, condition = "p95", type = "response") +
  labs(x = "p95 (m)", y = "P(hotspot)")

h2 <- plot_predictions(spring_m_hot, condition = "sdp95", type = "response") +
  labs(x = "sd p95 (m)", y = "P(hotspot)")

(g1 | g2) / (h1 | h2)

