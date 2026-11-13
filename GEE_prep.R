# Stopover intensity vs LiDAR-derived forest structure: GEE PREP (run order 1 of 5)
# Makes the shared 300 m grid (lidar_300m.tif), the radar station table (stations.csv), radar coverage (radar_300m.tif), and the GEE export region (gee_region.shp in NE_gee/, made only when not found there)
# Next: upload gee_region.shp (with its .shx, .dbf, and .prj) as the GEE table asset gee_region, run chili_dw_export.js, and put the exported tiles in NE_gee/

###################################################################
# Load Packages
{
  library(terra)
  library(dplyr)
}



##################################################################
# SETUP: folders, the 10 m LiDAR grid, and the shared 300 m grid
{
  # Folders: radar rasters, inputs and outputs, and the GEE exports (with the GEE region)
  radar_dir <- "C:/Users/willd/Box/"
  out_dir <- "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/"
  gee_dir <- paste0(out_dir, "NE_gee/")
  dir.create(gee_dir, showWarnings = F)
  
  # LiDAR datasets (folder + file-name prefix); where they overlap, the last one listed is used
  lidar_src <- c("C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/JT_LiDAR/Northeast_epsg6350_",
                 "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/Burch_LiDAR/Burch/mosaic_epsg6350_")
  
  # Processing: small blocks, large GDAL cache (the JT files are stored in 512 x 512 tiles)
  terraOptions(memmax = 2)
  gdalCache(12000)
  
  # p95 of both LiDAR datasets as one 10 m virtual mosaic (the reference grid)
  p95 <- vrt(paste0(lidar_src, "p95_10m.tif"))
  
  # Shared 300 m grid (the 10 m grid aggregated by 30), saved once as the maximum p95 of each 300 m cell
  f_lid <- paste0(out_dir, "lidar_300m.tif")
  if (!file.exists(f_lid)) aggregate(p95, 30, fun = "max", na.rm = T, filename = f_lid)
  grid300 <- rast(f_lid)
  
  # LiDAR footprint (bounding box) for projecting
  foot <- as.polygons(ext(p95), crs = crs(p95))
}



##################################################################
# RADAR STATIONS (saved as stations.csv in out_dir): stations with spring VIR inside the LiDAR footprint, numbered 1, 2, ... (the ids in station.tif)
{
  # Yearly VIR raster (one layer per year, 2013-2024) of station s and season
  fp <- function(s, season) sprintf("%sDFF_ppi_years/%s_2013-24_%s_bioRad_ppi_years_AGL_inverse_Refit.tif", radar_dir, s, season)
  
  f_sta <- paste0(out_dir, "stations.csv")
  if (!file.exists(f_sta)) {
    # Spring VIR cells of station s inside the LiDAR footprint (crop() errors when there is no overlap, so tryCatch gives 0)
    n_vir <- function(s) {
      r <- rast(fp(s, "spring"))
      tryCatch(sum(!is.na(values(crop(r, project(foot, crs(r)))))), error = function(e) 0)
    }
    
    # NOAA NCEI station list (ICAO code, latitude, longitude, site elevation in feet) -> DFF stations with VIR inside the footprint -> id and elevation (m)
    read.fwf("https://www.ncei.noaa.gov/access/homr/file/nexrad-stations.txt", widths = c(-9, 4, -93, 9, -1, 10, -1, 6), skip = 2,
             col.names = c("icao", "lat", "lon", "elev_ft"), strip.white = T) %>%
      filter(icao %in% read.csv(paste0(radar_dir, "DFF_resources/DFF_radar_BCR_list.csv"))$radar) %>%
      filter(sapply(icao, n_vir) > 0) %>%
      mutate(id = row_number(), elev_m = elev_ft * 0.3048) %>%
      write.csv(f_sta, row.names = F)
  }
  
  # Station table
  sta <- read.csv(f_sta)
}



##################################################################
# RADAR COVERAGE (saved as radar_300m.tif in out_dir): 1 = 300 m cell with VIR from any station, any season, any year; NA = not covered
{
  # Check: CRS, resolution, and grid origin of every station raster; one row per distinct combination, with its number of files (n)
  # mosaic() here and in stopover_lidar_covariates.R combines the station rasters as one grid, so 1 row is expected
  lapply(c(fp(sta$icao, "spring"), fp(sta$icao, "autumn")), function(f) {
    r <- rast(f)
    data.frame(crs = crs(r, proj = T), res = paste(round(res(r), 6), collapse = " x "), origin = paste(round(origin(r), 3), collapse = ", "))
  }) %>%
    bind_rows() %>%
    count(crs, res, origin) %>%
    print()
  
  f_rad <- paste0(out_dir, "radar_300m.tif")
  if (!file.exists(f_rad)) {
    c(fp(sta$icao, "spring"), fp(sta$icao, "autumn")) %>%
      lapply(rast) %>% sprc() %>% mosaic(fun = "max") %>% # all stations and seasons, one layer per year
      max(na.rm = T) %>%  # any year
      project(grid300, method = "bilinear") %>% # 300 m grid (bilinear, as the VIR responses)
      noNA(falseNA = T) %>% # 1 = covered, NA = not covered
      writeRaster(f_rad)
  }
  rad <- rast(f_rad)
}



##################################################################
# GEE REGION (saved as gee_region.shp in NE_gee/ when not found there): 300 m cells with LiDAR inside radar coverage, as polygons grown by one 300 m cell and simplified to 300 m
# Growing and simplifying cuts the vertices GEE has to handle (complex regions made the exports fail) while every LiDAR cell stays inside the region
{
  f_reg <- paste0(gee_dir, "gee_region.shp")
  if (!file.exists(f_reg)) {
    mask(rad, grid300) %>%   # covered cells with LiDAR
      as.polygons() %>%      # one polygon (all cells = 1)
      buffer(300) %>%
      simplifyGeom(300) %>%
      writeVector(f_reg)
  }
  
  # Check: vertices, area (km2), and map of the LiDAR cells inside radar coverage (maximum p95) with the region outline in red
  region <- vect(f_reg)
  cat("gee_region:", nrow(geom(region)), "vertices,", round(expanse(region, unit = "km")), "km2\n")
  plot(mask(grid300, rad), main = "LiDAR inside radar coverage (max p95 per 300 m cell); GEE region in red")
  lines(region, col = "red")
}