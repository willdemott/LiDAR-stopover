# Stopover intensity vs LiDAR-derived forest structure: 300 m COVARIATES
# Each covariate is its own 300 m raster on one shared grid (EPSG:6350, 300 m cells on the LiDAR 10 m grid origin), saved to cov_dir as <name>.tif
# LiDAR covariates use the kept forest pixels of each 300 m cell (woody, flown leaf-off, p95 >= h_forest); GEE and radar covariates describe the whole cell
# A chunk skips its work when its file exists: delete a file to rebuild that covariate (after a parameter change, also delete the masks in out_dir it depends on)
# Pipes (%>%) read left to right: take a raster or table, do one step, pass the result to the next step

###################################################################
# Load Packages
{
  library(terra)
  library(dplyr)
  library(bioRad)
}



##################################################################
# SETUP: folders, parameters, the 10 m LiDAR grid, and the shared 300 m grid
{
  # Folders: radar rasters, inputs and saved masks (LiDAR, GEE exports, first-leaf layer), and the 300 m covariates
  radar_dir <- "C:/Users/willd/Box/"
  out_dir <- "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/"
  cov_dir <- paste0(out_dir, "cov300/")
  dir.create(cov_dir, showWarnings = F)
  
  # File of a 300 m covariate, e.g. cf("p95")
  cf <- function(name) paste0(cov_dir, name, ".tif")
  
  # LiDAR datasets (folder + file-name prefix); where they overlap, the last one listed is used
  lidar_src <- c("C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/JT_LiDAR/Northeast_epsg6350_",
                 "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/Burch_LiDAR/Burch/mosaic_epsg6350_")
  
  # Parameters
  woody_min <- 75    # % of a 300 m cell's 10 m pixels classed woody (WorldCover tree cover or shrubland)
  h_scrub   <- 1     # p95 (m) from which a pixel counts as vegetation (scrub)
  h_forest  <- 5     # p95 (m) from which a pixel counts as forest
  leaf_lag  <- 14    # days after the local first-leaf day from which a flight counts as leaf-on
  leafoff   <- 1001  # MMDD from which a flight counts as leaf-off again (Oct 1)
  brk <- c(1, 2, 3, 4, 5, 6, 18, 30, 42, Inf)   # height layer breaks (m), shared by both datasets' V bins
  
  # Processing: small blocks, large GDAL cache (the JT files are stored in 512 x 512 tiles)
  terraOptions(memmax = 2)
  gdalCache(12000)
  
  # LiDAR metric v of all datasets as one 10 m virtual mosaic, e.g. lr("p95")
  lr <- function(v) vrt(paste0(lidar_src, v, "_10m.tif"))
  p95 <- lr("p95")
  
  # Any file(s) set to the 10 m LiDAR grid (virtual): tiles are stitched and every 10 m layer lines up with p95
  va <- function(f) vrt(f, options = c("-te", xmin(p95), ymin(p95), xmax(p95), ymax(p95), "-tr", res(p95)))
  
  # Shared 300 m grid (the 10 m grid aggregated by 30), saved once as the maximum p95 of each 300 m cell
  f_lid <- paste0(out_dir, "lidar_300m.tif")
  if (!file.exists(f_lid)) aggregate(p95, 30, fun = "max", na.rm = T, filename = f_lid)
  grid300 <- rast(f_lid)
  
  # LiDAR footprint (bounding box)
  foot <- as.polygons(ext(p95), crs = crs(p95))
}



##################################################################
# RADAR STATIONS (saved as cov300/stations.csv): stations with spring VIR inside the LiDAR footprint, numbered 1, 2, ... (the ids in station.tif)
{
  # Yearly VIR raster (one layer per year, 2013-2024) of station s and season
  fp <- function(s, season) sprintf("%sDFF_ppi_years/%s_2013-24_%s_bioRad_ppi_years_AGL_inverse_Refit.tif", radar_dir, s, season)
  
  f_sta <- paste0(cov_dir, "stations.csv")
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
  
  # Station table, and station points in the LiDAR CRS
  sta <- read.csv(f_sta)
  sta_v <- vect(sta, geom = c("lon", "lat"), crs = "EPSG:4326") %>% project(crs(p95))
}



##################################################################
# RADAR COVERAGE (saved as radar_300m.tif in out_dir): 1 = 300 m cell with VIR from any station, either season, any year; NA = not covered
{
  f_rad <- paste0(out_dir, "radar_300m.tif")
  if (!file.exists(f_rad)) {
    c(fp(sta$icao, "spring"), fp(sta$icao, "autumn")) %>%
      lapply(rast) %>% sprc() %>% mosaic(fun = "max") %>%   # all stations and seasons, one layer per year
      max(na.rm = T) %>%                                     # any year
      project(grid300, method = "bilinear") %>%              # 300 m grid (bilinear, as the VIR responses)
      noNA(falseNA = T) %>%                                  # 1 = covered, NA = not covered
      writeRaster(f_rad)
  }
  rad <- rast(f_rad)
}



##################################################################
# WOODY CELLS (saved as woody_cells_300m.tif in out_dir): 1 = 300 m cell inside radar coverage with at least woody_min % woody 10 m pixels
{
  # 10 m woody pixels (100 = WorldCover tree cover or shrubland, 0 = other): GEE tiles stitched on the 10 m grid
  woody10 <- va(list.files(out_dir, pattern = "^gee_woody_10m.*\\.tif$", full.names = T))
  
  f_cell <- paste0(out_dir, "woody_cells_300m.tif")
  if (!file.exists(f_cell)) {
    # % woody of each 300 m cell (mean of the 0/100 pixels; 75% = at least 675 of the 900 pixels)
    w300 <- aggregate(woody10, 30, fun = "mean", na.rm = T)
    ifel(w300 >= woody_min & rad == 1, 1, NA, filename = f_cell, datatype = "INT1U")
  }
}



##################################################################
# KEPT PIXELS (saved as veg_10m.tif in out_dir): the 10 m pixels behind the LiDAR covariates; 2 = forest (p95 >= h_forest), 1 = scrub (h_scrub <= p95 < h_forest), NA = not used
# Kept = woody pixel, inside a woody 300 m cell, p95 >= h_scrub, and both flights (MinTime, MaxTime) leaf-off; leaf-on = on/after the local first-leaf day + leaf_lag and before leafoff
{
  f_veg <- paste0(out_dir, "veg_10m.tif")
  if (!file.exists(f_veg)) {
    # Local first-leaf day of year (USA-NPN 30-year average, ~4 km): cropped to the LiDAR footprint, projected to the LiDAR CRS (small file), set to the 10 m grid
    fl <- rast(paste0(out_dir, "si-x__30yr_avg_4k_leaf.tif")) %>%
      crop(ext(project(foot, crs(.))) + 0.1) %>%
      project(crs(p95), filename = tempfile(fileext = ".tif")) %>%
      sources() %>%
      va()
    
    # Day of year of a flight date (YYMMDD)
    doy <- function(t) c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)[t %% 10000 %/% 100] + t %% 100
    
    # One pass over the 10 m grid: w = woody, k = woody cell, h = p95, t0 / t1 = first / last flight date, f = first-leaf day (May 1 = day 121 where missing)
    lapp(c(woody10, va(f_cell), p95, lr("MinTime"), lr("MaxTime"), fl),
         function(w, k, h, t0, t1, f) {
           f <- ifelse(is.na(f), 121, f) + leaf_lag
           on0 <- doy(t0) >= f & t0 %% 10000 < leafoff
           on1 <- doy(t1) >= f & t1 %% 10000 < leafoff
           ifelse(w == 100 & k == 1 & h >= h_scrub & !on0 & !on1, ifelse(h >= h_forest, 2, 1), NA)
         },
         filename = f_veg, wopt = list(datatype = "INT1U"))
  }
  veg <- rast(f_veg)
}



##################################################################
# CANOPY HEIGHT (p95.tif, sdp95.tif): mean and SD of p95 (m) over each 300 m cell's forest pixels
{
  if (!all(file.exists(cf(c("p95", "sdp95"))))) {
    tf <- tempfile(fileext = ".tif")
    h <- mask(p95, veg, maskvalues = c(NA, 1), filename = tf)   # forest pixels only (10 m, temporary file)
    h %>% aggregate(30, fun = "mean", na.rm = T) %>% setNames("p95") %>% writeRaster(cf("p95"), overwrite = T)
    h %>% aggregate(30, fun = "sd", na.rm = T) %>% setNames("sdp95") %>% writeRaster(cf("sdp95"), overwrite = T)
    unlink(tf)
  }
}



##################################################################
# RUMPLE INDEX (rumple.tif): mean over each 300 m cell's forest pixels of 1 / cos(slope of the 10 m p95 surface) = canopy surface area / flat area (1 = flat canopy)
{
  if (!file.exists(cf("rumple"))) {
    tf <- tempfile(fileext = c("_slope.tif", "_rumple.tif"))
    terrain(p95, "slope", unit = "radians", filename = tf[1]) %>%                # slope of the canopy surface (3 x 3 neighborhood)
      c(veg) %>%
      lapp(function(s, v) ifelse(v == 2, 1 / cos(s), NA), filename = tf[2]) %>%  # forest pixels only
      aggregate(30, fun = "mean", na.rm = T) %>%
      setNames("rumple") %>% writeRaster(cf("rumple"))
    unlink(tf)
  }
}



##################################################################
# POINT DENSITY (pd.tif): mean returns per m2 (NumPoints / 100 m2) over each 300 m cell's forest pixels
{
  if (!file.exists(cf("pd"))) {
    tf <- tempfile(fileext = ".tif")
    c(lr("NumPoints"), veg) %>%
      lapp(function(n, v) ifelse(v == 2, n / 100, NA), filename = tf) %>%
      aggregate(30, fun = "mean", na.rm = T) %>%
      setNames("pd") %>% writeRaster(cf("pd"))
    unlink(tf)
  }
}



##################################################################
# FLIGHT MONTH (yymm.tif): most common month-year (YYMM) of the first flight (MinTime) over each 300 m cell's forest pixels
{
  if (!file.exists(cf("yymm"))) {
    tf <- tempfile(fileext = ".tif")
    c(lr("MinTime"), veg) %>%
      lapp(function(t, v) ifelse(v == 2, t %/% 100, NA), filename = tf) %>%
      aggregate(30, fun = "modal", na.rm = T) %>%
      setNames("yymm") %>% writeRaster(cf("yymm"))
    unlink(tf)
  }
}



##################################################################
# FLIGHT YEAR (year.tif): year (YY) of each 300 m cell's most common flight month
{
  if (!file.exists(cf("year"))) {
    (rast(cf("yymm")) %/% 100) %>% setNames("year") %>% writeRaster(cf("year"))
  }
}



##################################################################
# HEIGHT LAYERS (used by the VOXEL, NRD, and FHD chunks): % of all returns in each height layer, on the 10 m grid
# Both datasets' V bins (JT: 1 m to 6 m, then 3 m; Burch: 1 m to 6 m, then 4 m) nest in the layers 1-2, ..., 5-6, 6-18, 18-30, 30-42, 42+ m
{
  # V bins: file, dataset (position in lidar_src), bin top (m), and height layer (1 = 0-1 m, 2 = 1-2 m, ..., 10 = 42+ m)
  vb <- lapply(seq_along(lidar_src), function(d) {
    data.frame(file = list.files(dirname(lidar_src[d]), pattern = paste0("^", basename(lidar_src[d]), "V[0-9]+_10m\\.tif$"), full.names = T), ds = d)
  }) %>%
    bind_rows() %>%
    mutate(top = as.numeric(sub(".*_V([0-9]+)_10m\\.tif$", "\\1", file)),
           layer = cut(top, c(0, brk), labels = F))
  
  # Check: V bins per dataset (rows) and height layer (columns); a 0 means that dataset has no bins in that layer
  print(table(dataset = vb$ds, layer = vb$layer))
  
  # Height layer names (above 1 m) and widths (m; the open top layer runs to the highest V bin)
  vl <- paste0("v", head(brk, -1), "_", brk[-1])
  vw <- diff(pmin(brk, max(vb$top)))
  
  # % of all returns in height layer L as one 10 m raster: each dataset's bins in L summed (saved once in out_dir when there are several), datasets joined (the later dataset wins where they overlap)
  layer10 <- function(L) {
    lapply(seq_along(lidar_src), function(d) {
      b <- vb$file[vb$ds == d & vb$layer == L]
      s <- paste0(out_dir, "vsum_ds", d, "_layer", L, "_10m.tif")
      if (length(b) > 1 & !file.exists(s)) app(rast(b), "sum", filename = s)
      if (length(b) > 1) s else b
    }) %>%
      unlist() %>%
      va()
  }
}



##################################################################
# VOXEL (voxel.tif, one layer per height layer above 1 m): mean % of all returns in each height layer over each 300 m cell's kept pixels (forest and scrub)
{
  if (!file.exists(cf("voxel"))) {
    tf <- tempfile(fileext = ".tif")
    lapply(2:length(brk), function(L) {
      layer10(L) %>%
        mask(veg, filename = tf, overwrite = T) %>%   # kept pixels only (10 m, temporary file)
        aggregate(30, fun = "mean", na.rm = T)
    }) %>%
      rast() %>% setNames(vl) %>% writeRaster(cf("voxel"))
    unlink(tf)
  }
}



##################################################################
# SHRUB LAYER (nrd.tif): mean % of all returns at 1-2 m (height layer 2 = V bin V2 of both datasets) over each 300 m cell's forest pixels
{
  if (!file.exists(cf("nrd"))) {
    tf <- tempfile(fileext = ".tif")
    c(layer10(2), veg) %>%
      lapp(function(a, v) ifelse(v == 2, a, NA), filename = tf) %>%   # 1-2 m of forest pixels
      aggregate(30, fun = "mean", na.rm = T) %>%
      setNames("nrd") %>% writeRaster(cf("nrd"))
    unlink(tf)
  }
}



##################################################################
# VERTICAL DIVERSITY (fhd.tif, sdfhd.tif): mean and SD of FHD over each 300 m cell's forest pixels
# FHD of a pixel = Shannon index of its returns over 1 m slices from 1 m up (each layer's returns spread evenly over its width vw),
# divided by log(M), M = meters from 1 m to the top of its highest layer with returns (1 = returns spread evenly up to that top)
{
  if (!all(file.exists(cf(c("fhd", "sdfhd"))))) {
    tf <- tempfile(fileext = ".tif")
    h <- c(rast(lapply(2:length(brk), layer10)), veg) %>%
      lapp(function(...) {
        m <- cbind(...)                                            # one row per pixel: % of returns per height layer, then veg
        if (all(is.na(m[, ncol(m)]))) return(rep(NA, nrow(m)))     # skip blocks without kept pixels
        v <- m[, -ncol(m), drop = F]
        tot <- rowSums(v)
        M <- cumsum(vw)[max.col((v > 0) * 1, ties.method = "last")]
        fhd <- (log(tot) - rowSums(ifelse(v > 0, v * log(v), 0)) / tot + drop(v %*% log(vw)) / tot) / log(M)
        ifelse(m[, ncol(m)] == 2 & tot > 0 & M > 1, fhd, NA)        # forest pixels only
      }, filename = tf)
    h %>% aggregate(30, fun = "mean", na.rm = T) %>% setNames("fhd") %>% writeRaster(cf("fhd"), overwrite = T)
    h %>% aggregate(30, fun = "sd", na.rm = T) %>% setNames("sdfhd") %>% writeRaster(cf("sdfhd"), overwrite = T)
    unlink(tf)
  }
}



##################################################################
# HEAT LOAD (chili.tif): mean CHILI (0 = cool to 255 = warm) of each 300 m cell, from the 90 m GEE export
{
  if (!file.exists(cf("chili"))) {
    vrt(list.files(out_dir, pattern = "^gee_chili_90m.*\\.tif$", full.names = T)) %>%   # GEE tiles stitched
      project(grid300, method = "average") %>%                                         # mean over each 300 m cell
      setNames("chili") %>% writeRaster(cf("chili"))
  }
}



##################################################################
# DISTANCE TO WATER (waterdist.tif): mean distance (m) to the nearest water of each 300 m cell, from the 30 m GEE export
{
  if (!file.exists(cf("waterdist"))) {
    vrt(list.files(out_dir, pattern = "^gee_waterdist_30m.*\\.tif$", full.names = T)) %>%
      project(grid300, method = "average") %>%
      setNames("waterdist") %>% writeRaster(cf("waterdist"))
  }
}



##################################################################
# ELEVATION (elev.tif): mean elevation (m) of each 300 m cell, from the 30 m GEE export (3DEP)
{
  if (!file.exists(cf("elev"))) {
    vrt(list.files(out_dir, pattern = "^gee_elev_30m.*\\.tif$", full.names = T)) %>%
      project(grid300, method = "average") %>%
      setNames("elev") %>% writeRaster(cf("elev"))
  }
}



##################################################################
# FOREST TYPE (fortype.tif): most common TreeMap forest type code (FIA FORTYPCD) of each 300 m cell's TreeMap forest pixels, from the 30 m GEE export
# Values up to 100 are not forest type codes (the lowest is 101, jack pine; e.g. 0 or -1 for non-forest) and are set to NA before the mode; cells without forest pixels stay NA
{
  if (!file.exists(cf("fortype"))) {
    tf <- tempfile(fileext = ".tif")
    vrt(list.files(out_dir, pattern = "^gee_fortype_30m.*\\.tif$", full.names = T)) %>%
      classify(cbind(-Inf, 100, NA), filename = tf, datatype = "INT2S") %>%   # non-forest values -> NA (30 m, temporary file)
      project(grid300, method = "mode") %>%
      setNames("fortype") %>% writeRaster(cf("fortype"))
    unlink(tf)
  }
}



##################################################################
# FOREST TYPE CLASS (ftype.tif): codes with at least 500 300 m cells of forest LiDAR keep their own class (most common first); all other codes = Other
# Names and colors: treemap_fortypcd.csv in out_dir (code, TreeMap color, name)
{
  if (!file.exists(cf("ftype"))) {
    # 300 m cells with forest LiDAR (p95) per forest type code -> codes with at least 500 cells, most common first
    ft <- read.csv(paste0(out_dir, "treemap_fortypcd.csv")) %>%
      inner_join(freq(mask(rast(cf("fortype")), rast(cf("p95")))), by = c("code" = "value")) %>%
      filter(count >= 500) %>%
      arrange(desc(count))
    
    # Class of each 300 m cell (class i = row i of ft; last class = Other), labeled code + name
    ftype <- classify(rast(cf("fortype")), cbind(ft$code, 1:nrow(ft)), others = nrow(ft) + 1)
    levels(ftype) <- data.frame(value = 1:(nrow(ft) + 1), ftype = c(paste(ft$code, ft$name), "Other"))
    ftype %>% setNames("ftype") %>% writeRaster(cf("ftype"))
  }
}



##################################################################
# RADAR DISTANCE (raddist.tif): distance (km) from each 300 m cell to the nearest station
{
  if (!file.exists(cf("raddist"))) {
    (distance(rast(grid300), sta_v) / 1000) %>% setNames("raddist") %>% writeRaster(cf("raddist"))
  }
}



##################################################################
# NEAREST STATION (station.tif): id (stations.csv) of each 300 m cell's nearest station
{
  if (!file.exists(cf("station"))) {
    lapply(1:nrow(sta), function(i) distance(rast(grid300), sta_v[i])) %>%   # distance to each station, one layer per station
      rast() %>% which.min() %>%                                            # layer of the nearest station = its id
      setNames("station") %>% writeRaster(cf("station"))
  }
}



##################################################################
# RADAR BEAM HEIGHT (beamagl.tif): height (m) of the nearest station's lowest beam (0.5 degrees) above the ground
# = station site elevation + beam height above the radar at that distance - cell elevation
{
  if (!file.exists(cf("beamagl"))) {
    site <- classify(rast(cf("station")), cbind(sta$id, sta$elev_m))                                              # site elevation (m)
    beam <- lapp(rast(cf("raddist")) * 1000, function(d) beam_height(beam_range(d, 0.5), 0.5, lat = mean(sta$lat)))  # beam height above the radar (m)
    (site + beam - rast(cf("elev"))) %>% setNames("beamagl") %>% writeRaster(cf("beamagl"))
  }
}



##################################################################
# CHECK: all 300 m covariates as one stack (fails if any file is off the shared grid), with summaries and maps
{
  cov300 <- rast(list.files(cov_dir, pattern = "\\.tif$", full.names = T))
  print(cov300)
  print(summary(cov300))
  plot(cov300)
}