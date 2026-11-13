# Stopover intensity vs LiDAR-derived forest structure: DATA PREP, 300 m COVARIATES and RESPONSES (run order 3 of 5)
# Each covariate (and each season's responses) is its own 300 m raster on one shared grid (EPSG:6350, 300 m cells on the LiDAR 10 m grid origin), cached in ./cov_dir as <name>.tif
# All cached layers are written as ONE stacked raster, stack_300m.tif in out_dir (with its .aux.xml, which holds the forest type and station labels): the only input of stopover_lidar_responses.R and stopover_lidar_model.R
# Needs: lidar_300m.tif, stations.csv (stopover_lidar_gee_prep.R), veg_10m.tif (stopover_lidar_mask.R), and the GEE exports (NE_gee/)
# LiDAR covariates only contain non-leaf forest pixels of each 300 m cell (woody, flown leaf-off, p95 >= h_forest)

###################################################################
# Load Packages
{
  library(terra)
  library(dplyr)
  library(bioRad)
  library(mapview)
}



##################################################################
# SETUP: load in folders, parameters, the 10 m LiDAR grid, the shared 300 m grid, the radar stations, and the kept pixels
{
  # Folders: radar rasters, inputs and starter-rasters (LiDAR, masks), the GEE exports, the directory for the 300m covariates and responses, and the summed V bins
  radar_dir <- "C:/Users/willd/Box/"
  out_dir <- "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/"
  gee_dir <- paste0(out_dir, "NE_gee/")
  cov_dir <- paste0(out_dir, "cov300/")
  vsum_dir <- paste0(out_dir, "vsum/")
  dir.create(cov_dir, showWarnings = F)
  dir.create(vsum_dir, showWarnings = F)
  
  # Function to call a 300 m covariate once made, e.g. cf("p95")
  cf <- function(name) paste0(cov_dir, name, ".tif")
  
  # LiDAR datasets (folder + file-name prefix); where they overlap, the last one listed is used
  lidar_src <- c("C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/JT_LiDAR/Northeast_epsg6350_",
                 "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/Burch_LiDAR/Burch/mosaic_epsg6350_")
  
  # Parameters: height layer breaks (m), shared by both datasets' V bins; seasons and radar years (layer 1 = 2013) of the responses
  brk <- c(1, 2, 3, 4, 5, 6, 18, 30, 42, Inf)
  seasons <- c("spring", "autumn")
  yrs <- 2013:2024
  
  # Processing: small blocks, large GDAL cache (the JT files are stored in 512 x 512 tiles)
  terraOptions(memmax = 2)
  gdalCache(12000)
  
  # A function that pulls LiDAR metrics from both datasets into one 10m virtual mosaic, e.g. lr("p95")
  lr <- function(v) vrt(paste0(lidar_src, v, "_10m.tif"))
  # We'll use p95 as a CRS reference-raster
  p95 <- lr("p95")
  
  # 'virtual alignment' function; tiles are stitched and every 10 m layer lines up with p95
  va <- function(f) vrt(f, options = c("-te", xmin(p95), ymin(p95), xmax(p95), ymax(p95), "-tr", res(p95)))
  
  # Shared 300 m grid: maximum p95 of each 300 m cell (lidar_300m.tif, made by stopover_lidar_gee_prep.R)
  grid300 <- rast(paste0(out_dir, "lidar_300m.tif"))
  
  # Radar stations (stations.csv, made by stopover_lidar_gee_prep.R; id = value in station.tif), and station points in the LiDAR CRS
  sta <- read.csv(paste0(out_dir, "stations.csv"))
  sta_v <- vect(sta, geom = c("lon", "lat"), crs = "EPSG:4326") %>% project(crs(p95))
  
  # Yearly VIR raster (one layer per year, 2013-2024) of station s and season
  fp <- function(s, season) sprintf("%sDFF_ppi_years/%s_2013-24_%s_bioRad_ppi_years_AGL_inverse_Refit.tif", radar_dir, s, season)
  
  # Kept 10 m pixels, the final mask (veg_10m.tif, made by stopover_lidar_mask.R): 2 = forest, 1 = scrub, NA = not used
  veg <- rast(paste0(out_dir, "veg_10m.tif"))
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
    terrain(p95, "slope", unit = "radians", filename = tf[1]) %>% # slope of the canopy surface (3 x 3 neighborhood)
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
  
  # % of all returns in height layer L as one 10 m raster: each dataset's bins in L summed (saved once in vsum_dir when there are several), datasets joined (the later dataset wins where they overlap)
  layer10 <- function(L) {
    lapply(seq_along(lidar_src), function(d) {
      b <- vb$file[vb$ds == d & vb$layer == L]
      s <- paste0(vsum_dir, "vsum_ds", d, "_layer", L, "_10m.tif")
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
# FHD of a pixel = Shannon index of its returns over 1 m slices from 1 m up (each layer's returns spread evenly over its width vw), divided by log(M), M = meters from 1 m to the top of its highest layer with returns (1 = returns spread evenly up to that top)
{
  if (!all(file.exists(cf(c("fhd", "sdfhd"))))) {
    tf <- tempfile(fileext = ".tif")
    h <- c(rast(lapply(2:length(brk), layer10)), veg) %>%
      lapp(function(...) {
        m <- cbind(...) # one row per pixel: % of returns per height layer, then veg
        if (all(is.na(m[, ncol(m)]))) return(rep(NA, nrow(m))) # skip blocks without kept pixels
        v <- m[, -ncol(m), drop = F]
        tot <- rowSums(v)
        M <- cumsum(vw)[max.col((v > 0) * 1, ties.method = "last")]
        fhd <- (log(tot) - rowSums(ifelse(v > 0, v * log(v), 0)) / tot + drop(v %*% log(vw)) / tot) / log(M)
        ifelse(m[, ncol(m)] == 2 & tot > 0 & M > 1, fhd, NA) # forest pixels only
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
    vrt(list.files(gee_dir, pattern = "^gee_chili_90m.*\\.tif$", full.names = T)) %>%   # GEE tiles stitched
      project(grid300, method = "average") %>% # mean over each 300 m cell
      setNames("chili") %>% writeRaster(cf("chili"))
  }
}



##################################################################
# DISTANCE TO WATER (waterdist.tif): mean distance (m) to the nearest water of each 300 m cell, from the 30 m GEE export
{
  if (!file.exists(cf("waterdist"))) {
    vrt(list.files(gee_dir, pattern = "^gee_waterdist_30m.*\\.tif$", full.names = T)) %>%
      project(grid300, method = "average") %>%
      setNames("waterdist") %>% writeRaster(cf("waterdist"))
  }
}



##################################################################
# ELEVATION (elev.tif): mean elevation (m) of each 300 m cell, from the 30 m GEE export (3DEP)
{
  if (!file.exists(cf("elev"))) {
    vrt(list.files(gee_dir, pattern = "^gee_elev_30m.*\\.tif$", full.names = T)) %>%
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
    vrt(list.files(gee_dir, pattern = "^gee_fortype_30m.*\\.tif$", full.names = T)) %>%
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
      rast() %>% which.min() %>%  # layer of the nearest station = its id
      setNames("station") %>% writeRaster(cf("station"))
  }
}



##################################################################
# RADAR BEAM HEIGHT (beamagl.tif): height (m) of the nearest station's lowest beam (0.5 degrees) above the ground
# = station site elevation + beam height above the radar at that distance - cell elevation
{
  if (!file.exists(cf("beamagl"))) {
    site <- classify(rast(cf("station")), cbind(sta$id, sta$elev_m)) # site elevation (m)
    beam <- lapp(rast(cf("raddist")) * 1000, function(d) beam_height(beam_range(d, 0.5), 0.5, lat = mean(sta$lat)))  # beam height above the radar (m)
    (site + beam - rast(cf("elev"))) %>% setNames("beamagl") %>% writeRaster(cf("beamagl"))
  }
}



##################################################################
# HOT-YEAR CUTOFFS (saved as cutoffs.csv in cov_dir): 90th percentile of each station's VIR in each year and season, over its full range (the station's own raster: uncropped and unmasked)
{
  f_cut <- paste0(cov_dir, "cutoffs.csv")
  if (!file.exists(f_cut)) {
    lapply(seasons, function(s) {
      lapply(sta$id, function(i) {
        q90 <- apply(values(rast(fp(sta$icao[i], s))), 2, quantile, 0.9, na.rm = T)   # one value per year (layer)
        data.frame(season = s, id = i, icao = sta$icao[i], year = yrs, q90 = q90)
      }) %>% bind_rows()
    }) %>%
      bind_rows() %>%
      write.csv(f_cut, row.names = F)
  }
  
  # One row per season, station, and year
  cutoffs <- read.csv(f_cut)
  print(head(cutoffs))
}



##################################################################
# RESPONSES (resp_<season>.tif, layers vir_<season>, hot_<season>, nyr_<season>): vir = mean VIR of the radar years within the cell's LiDAR year +-1;
# hot = 1 if the cell's VIR reached its nearest station's cutoff in all 3 of those years; nyr = years tested (VIR and a cutoff), so hot = 1 needs nyr = 3
{
  # Radar years within the LiDAR year +-1 of each 300 m cell: one layer per radar year (1 = within, 0 = not)
  win <- rast(lapply(yrs, function(y) abs(y - 2000 - rast(cf("year"))) <= 1))
  
  for (s in seasons) {
    if (!file.exists(cf(paste0("resp_", s)))) {
      # Yearly VIR on the 300 m grid (one layer per year; mean where station ranges overlap; bilinear), kept in the LiDAR window only
      vir <- fp(sta$icao, s) %>%
        lapply(rast) %>% sprc() %>% mosaic(fun = "mean") %>%
        project(grid300, method = "bilinear") %>%
        mask(win, maskvalues = c(0, NA))
      
      # Hot-year cutoff of each cell's nearest station: one layer per year
      cut_r <- rast(lapply(yrs, function(y) {
        rcl <- cutoffs %>% filter(season == s, year == y) %>% select(id, q90) %>% as.matrix()   # station id -> cutoff
        classify(rast(cf("station")), rcl)
      }))
      
      # Hot years: 1 = VIR >= cutoff, 0 = below, NA = no VIR or no cutoff that year
      hot <- vir >= cut_r
      
      # Responses: mean VIR of the window years, hot in all 3 window years (1/0), years tested
      c(mean(vir, na.rm = T), (sum(hot, na.rm = T) == 3) * 1, sum(!is.na(hot))) %>%
        setNames(paste0(c("vir_", "hot_", "nyr_"), s)) %>%
        writeRaster(cf(paste0("resp_", s)))
    }
  }
}



##################################################################
# STACK (saved as stack_300m.tif in out_dir, rewritten on every run so it always matches the cached layers): every 300 m layer in cov_dir as ONE raster
# The station layer gets its ICAO codes as labels (stations.csv), so the model and figure scripts need only the stack
{
  stk <- rast(list.files(cov_dir, pattern = "\\.tif$", full.names = T))   # fails if any file is off the shared grid
  levels(stk[["station"]]) <- data.frame(id = sta$id, station = sta$icao)
  writeRaster(stk, paste0(out_dir, "stack_300m.tif"), overwrite = T)
  
  # Check: layers and their summaries
  print(stk)
  print(summary(stk))
}