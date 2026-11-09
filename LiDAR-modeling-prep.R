# Stopover intensity vs LiDAR-derived forest structure: DATA PREP
# 10 m LiDAR + GEE layers + radar VIR -> 300 m covariates and responses -> spring_d.csv and autumn_d.csv
# Study pixels = 10 m LiDAR pixels inside radar coverage, classed WorldCover woody (tree cover or shrubland), in 300 m cells with >= woody_min % woody pixels, and with p95 >= h_scrub
# Forest = study pixels with p95 >= h_forest; scrub = study pixels with p95 below h_forest; only pixels flown leaf-off are used
# Covariates from the 10 m pixels use forest pixels only (elevation and radar beam height describe the whole 300 m cell); scrub pixels count only in the height-layer (voxel) layers
# Slow steps are saved in out_dir and read back on later runs (delete a saved file to rebuild it)

###################################################################
# Load Packages
{
  library(terra)
  library(mapview)
  library(viridisLite)
  library(ggplot2)
  library(ggridges)
  library(bioRad)
}



##################################################################
# SETUP: directories, parameters, saved files, and the 10 m LiDAR grid
{
  # Where radar files are kept
  radar_dir <- "C:/Users/willd/Box/"
  
  # Where the GEE exports, first-leaf layer, saved files, and model dataframes are kept
  out_dir <- "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/"
  
  # LiDAR datasets: folder + file-name prefix of each; where they overlap, the last one listed is used
  # Des Moines: "C:/Users/willd/Box/LiDAR/IA_lidar_4_freda/mosaic_epsg6344_"
  lidar_src <- c("C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/JT_LiDAR/Northeast_epsg6350_", "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/Burch_LiDAR/Burch/mosaic_epsg6350_")
  
  # Saved files
  f_cov   <- paste0(out_dir, "covariates_300m.tif") # 300 m covariates
  
  # Leaf-on = flown on/after the local first-leaf day + leaf_lag (days) and before leafoff (MMDD)
  leaf_lag <- 14 # 2 weeks
  leafoff <- 1001 #oct 1st
  
  # Day of year at the start of each month (flight date MMDD -> day of year)
  m0 <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
  
  # Woody cutoff (% of a 300 m cell's 10 m pixels classed WorldCover tree cover or shrubland) and height cutoffs (p95, m)
  woody_min <- 75 #percent
  h_scrub <- 1 # meters
  h_forest <- 5 # meters
  
  # Height layer breaks (m) chosen by both datasets' V bins
  brk <- c(1, 2, 3, 4, 5, 6, 18, 30, 42, Inf)
  vl <- paste0("v", head(brk, -1), "_", brk[-1])
  
  # 10 m LiDAR grid: p95 of all dataset(s) as one virtual mosaic (vrt() = very really terrific).
  # This object will be used as the reference for aligning CRS between datasets
  p95_file <- vrt(paste0(lidar_src, "p95_10m.tif"))
  
  # While were at it, lets make a function that can easily turn any of our lidar data  ('v') into that virtual raster.
  lr <- function(v) vrt(paste0(lidar_src, v, "_10m.tif"))
  
  # We want to make sure these virtual rasters are aligned geographically, so lets make one last function that does just that. 
  # 'te' = target extent, 'tr' = target resolution
  va <- function(f) vrt(f, options = c("-te", xmin(p95_file), ymin(p95_file), xmax(p95_file), ymax(p95_file), "-tr", res(p95_file)))
  
  # Now we can make our voxel objects that we'll work with. 'vb' is a reference dataframe for V bins of both datasets (burch/JTs): file, dataset (position in lidar_src), bin top (m), and height layer (1 = 0-1 m, then one per pair of brk edges)
  vb <- do.call(rbind, lapply(seq_along(lidar_src), function(i) {
    f <- list.files(dirname(lidar_src[i]), pattern = paste0("^", basename(lidar_src[i]), "V[0-9]+_10m\\.tif$"), full.names = T)
    data.frame(file = f, ds = i, top = as.numeric(sub(".*_V([0-9]+)_10m\\.tif$", "\\1", f)))
  }))
  
  vb$layer <- cut(vb$top, c(0, brk), labels = F)
  
  # Bin -> layer matrix (1 = bin is in that layer): a dataset's bins %*% G = its % of returns per height layer
  G <- outer(vb$layer, 1:length(brk), "==") * 1
  
  # Width (m) of each height layer above 1 m, used to weight the layers in FHD (1-6 m layers = 1 m, 6-42 m layers = 12 m); the open 42+ m layer runs to the top of the highest V bin
  vw <- diff(pmin(brk, max(vb$top)))
  
}



##################################################################
# LiDAR FOOTPRINT, RADAR STATIONS, and GEE REGION
# GEE region = 300 m cells with LiDAR inside radar coverage, saved as gee_region.shp (lon/lat) for chili_dw_export.js. After running this chunk, use gee_region.shp to export covariates from code.googleearthengine.com, and drag them into out_dir before continuing with script.
{
  setwd(radar_dir)
  
  # Stations
  stations <- unique(read.csv("./DFF_resources/DFF_radar_BCR_list.csv")$radar)
  
  # Function to reference the file of a station-year-season's raster
  fp <- function(s, season) sprintf("./DFF_ppi_years/%s_2013-24_%s_bioRad_ppi_years_AGL_inverse_Refit.tif", s, season)
  
  # LiDAR footprint/bounding box
  foot_p <- as.polygons(ext(p95_file), crs = crs(p95_file))
  
  # Spring VIR values of each station inside the LiDAR footprint (crop() errors when there's no overlap, so tryCatch returns 0); stations with values become those included
  n_values <- sapply(stations, function(s) {
    r <- rast(fp(s, "spring"))
    tryCatch(sum(!is.na(values(crop(r, project(foot_p, crs(r)))))), error = function(e) 0)
  })
  mos <- stations[n_values > 0]
  
  # 300 m cells with LiDAR: maximum p95 of each 300 m cell
  f_lid   <- paste0(out_dir, "lidar_300m.tif")
  if (!file.exists(f_lid)) 
    aggregate(p95_file, 30, fun = "max", na.rm = T, filename = f_lid)
  lid <- rast(f_lid)
  
  # Radar coverage: any VIR value of those stations (either season, any year); 300 m cells inside it (bilinear, as the VIR in STOPOVER RESPONSES), saved
  rcov <- max(mosaic(sprc(lapply(c(fp(mos, "spring"), fp(mos, "autumn")), rast)), fun = "max"), na.rm = T)
  f_rad   <- paste0(out_dir, "radar_300m.tif")
  rad <- ifel(!is.na(project(rcov, lid, method = "bilinear")), 1, NA,
              filename = f_rad, overwrite = T)
  
  # GEE region: LiDAR inside radar coverage as polygons, grown by one 300 m cell for covariate edge-case calculations, and simplified to 300 m; area (km2) and map (red outline)
  region <- simplifyGeom(buffer(as.polygons(mask(rad, lid)), 300), 300)  
  print(sum(expanse(region, unit = "km")))
  plot(lid)
  plot(region, add = T, border = "red", main = "p95")
  
  # Saved in lon/lat (upload the .shp, .shx, .dbf, .prj as a GEE table asset)
  if(!file.exists(paste0(out_dir, "gee_region.shp"))) 
    writeVector(project(region, "EPSG:4326"), 
                paste0(out_dir, "gee_region.shp"))
  
  # 10 m grid origin for GEE covariate exporting (XMIN, YMAX)
  print(c(xmin(p95_file), ymax(p95_file)), digits = 15)
}



##################################################################
# Loading in GEE layers: woody (10 m) for the study mask, CHILI, distance to water, and forest type
# (elevation, 30 m, is averaged to the 300 m grid in STOPOVER RESPONSES)
# GEE splits large exports into several tiles; va() stitches all tiles of a variable into one virtual raster on the 10 m LiDAR grid
{
  # Woody: 100 = WorldCover tree cover or shrubland, 0 = any other class (10 m, one value per LiDAR pixel)
  woody10 <- va(list.files(out_dir, pattern = "^gee_woody_10m.*\\.tif$", full.names = T))
  
  gee <- c(
    va(list.files(out_dir, pattern = "^gee_chili_90m.*\\.tif$", 
                  full.names = T)), 
    va(list.files(out_dir, pattern = "^gee_waterdist_30m.*\\.tif$", 
                  full.names = T)), 
    va(list.files(out_dir, pattern = "^gee_fortype_30m.*\\.tif$", 
                  full.names = T)))
  
  names(gee) <- c("chili", "waterdist", "fortype")
  
  compareGeom(c(woody10, gee), p95_file)
}



##################################################################
# WOODY STUDY MASK (saved as study_10m.tif): 1 = 10 m pixel classed woody, inside radar coverage, in a 300 m cell with at least woody_min % woody 10 m pixels
{
  # Saved files: 300 m cells kept (small) and 10 m study pixels
  f_cell  <- paste0(out_dir, "woody_cells_300m.tif")
  f_study <- paste0(out_dir, "study_10m.tif")
  
  # 300 m cells kept: mean of the 0/100 woody pixels = % woody (75% = at least 675 of the 900 pixels), at least woody_min, and inside radar coverage (rad, same 300 m grid)
  if (!file.exists(f_cell))
    lapp(c(aggregate(woody10, 30, fun = "mean", na.rm = T), rad), 
         function(w, r) ifelse(w >= woody_min & r == 1, 1, NA), 
         filename = f_cell, wopt = list(datatype = "INT1U"))
  
  # 10 m study pixels: woody pixels (100) inside kept 300 m cells (cell mask set to the 10 m grid)
  if (!file.exists(f_study))
    lapp(c(woody10, va(f_cell)), 
         function(w, k) ifelse(w == 100 & k == 1, 1, NA), 
         filename = f_study, wopt = list(datatype = "INT1U"))
}



##################################################################
# FLIGHT DATES of the study pixels: histogram of flight dates and map of flight month
{
  # MinTime (YYMMDD), study mask, and p95 of a regular sample of 10 m pixels, as a small raster; MinTime kept for study pixels only
  s <- spatSample(c(lr("MinTime"), va(f_study), p95_file), 1e6, method = "regular", as.raster = T)
  s <- mask(s[[1]], s[[2]] == 1 & s[[3]] >= h_scrub, maskvalues = c(0, NA))
  
  # Histogram of flight dates by week; every date set in 2000 so month/day line up across years (220415 -> Apr 15)
  hist(as.POSIXct(sprintf("2000%04d", as.integer(values(s %% 10000))),
                  format = "%Y%m%d", tz = "UTC"),
       breaks = "weeks", freq = T, format = "%b %d", 
       main = "Flight dates of study pixels (MinTime)", xlab = "")
  
  # Map of flight month (MMDD classes)
  plot(s %% 10000, type = "interval", 
       breaks = c(1:12 * 100 + 0.5, 1231.5), 
       col = c("#08306B", "#2171B5", "#9ECAE1", "#A1D99B", "#41AB5D",
               "#006D2C", "#00441B", "#006D2C", "#41AB5D", "#DD8D3C",
               "#6BAED6", "#4292C6"),
       plg = list(legend = month.abb), 
       main = "Flight month of study pixels (MinTime)")
}



#################################################################
# LEAF-OUT LAG FIGURE: how many study points up to the local first-leaf date + custom leaf_lag are kept?
{
  # Regular sample of 10 m pixels in the woody study mask, as points
  pts <- spatSample(va(f_study), 1e6, method = "regular", as.points = T, na.rm = T)
  
  # At the points: coordinates (x, y), p95, flight dates (YYMMDD), and % of all returns in the 1 m layers V2-V6 (1-2 ... 5-6 m)
  x <- data.frame(crds(pts), sapply(c("p95", "MinTime", "MaxTime", paste0("V", 2:6)), function(v) terra::extract(lr(v), pts, ID = F)[, 1]))
  
  # Local first-leaf day of year at the points; May 1 (121) where missing, based on LEAF-OFF DATE FILTER
  fl <- rast(paste0(out_dir, "si-x__30yr_avg_4k_leaf.tif"))
  x$fl <- terra::extract(fl, project(pts, crs(fl)), ID = F)[, 1]
  x$fl[is.na(x$fl)] <- 121
  
  # Study pixels only: p95 >= h_scrub and both flight dates (the pixels the LEAF-OFF DATE FILTER tests)
  x <- x[which(x$p95 >= h_scrub & !is.na(x$MinTime) & !is.na(x$MaxTime)), ]
  
  # Day of year of MinTime and MaxTime, days from first leaf to the last flight, and % of all returns at 1-6 m
  x$doy0 <- m0[x$MinTime %% 10000 %/% 100] + x$MinTime %% 100
  x$doy1 <- m0[x$MaxTime %% 10000 %/% 100] + x$MaxTime %% 100
  x$d <- x$doy1 - x$fl
  x$r16 <- x$V2 + x$V3 + x$V4 + x$V5 + x$V6
  
  # Figure pixels: forest, last flight 4 weeks before to 6 weeks after first leaf, all flights within 7 days (same year); their week from first leaf
  i <- which(x$p95 >= h_forest & x$d >= -28 & x$d < 42 & x$MaxTime %/% 10000 == x$MinTime %/% 10000 & x$doy1 - x$doy0 <= 7 & !is.na(x$r16))
  wk <- cut(x$d[i], seq(-28, 42, 7), right = F)
  
  # Weekly mean % of returns at 1-6 m; per week: forest pixels, flight dates, and that mean
  wm <- tapply(x$r16[i], wk, mean)
  print(data.frame(pixels = c(table(wk)), flight_dates = c(tapply(x$MaxTime[i], wk, function(t) length(unique(t)))), mean_pct_1_6m = c(round(wm, 2))))
  
  # Per flight date (MaxTime): mean days from first leaf, mean % of returns at 1-6 m, and forest pixels
  fd <- aggregate(cbind(d, r16) ~ MaxTime, data = x[i, ], FUN = mean)
  fd$n <- aggregate(d ~ MaxTime, data = x[i, ], FUN = length)$d
  
  # flight dates (gray, size ~ forest pixels) and weekly means (black); dotted line = first leaf; dashed line and shading = dropped as leaf-on
  plot(fd$d, fd$r16, type = "n", xlim = c(-28, 42), 
       xlab = "Days from local first leaf to last flight", 
       ylab = "% of returns at 1-6 m (study forest pixels)", 
       main = "A. Lower-canopy returns")
  rect(leaf_lag, par("usr")[3], par("usr")[2], par("usr")[4], 
       col = adjustcolor("orange", 0.2), border = NA)
  abline(v = c(0, leaf_lag), lty = c(3, 2))
  points(fd$d, fd$r16, pch = 16, col = "grey60", 
         cex = 0.5 + 2 * sqrt(fd$n / max(fd$n)))
  lines(seq(-25, 38, 7), wm, type = "b", pch = 15, lwd = 2)
  legend("topleft", c("Flight date (size ~ forest pixels)", 
                      "Weekly mean", "Dropped as leaf-on"), 
         pch = c(16, 15, 15), pt.cex = c(1, 1, 2),
         col = c("grey60", "black", adjustcolor("orange", 0.4)), 
         lty = c(NA, 1, NA), bty = "n", cex = 0.8)
  
  # Map: sampled study pixels kept (grey), with those dropped as leaf-on at leaf_lag (same test as panel B) on top in red
  x$drop <- (x$doy0 >= x$fl + leaf_lag & x$MinTime %% 10000 < leafoff) | (x$doy1 >= x$fl + leaf_lag & x$MaxTime %% 10000 < leafoff)
  
  plot(x$x, x$y, pch = 16, cex = 0.2, col = "grey70", asp = 1, 
       xlab = "", ylab = "", 
       main = paste0("Study pixels dropped as leaf-on (red), lag ",
                     leaf_lag, " days"))
  points(x$x[x$drop], x$y[x$drop], pch = 16, cex = 0.2, col = "red")
}



##################################################################
# LEAF-OFF DATE FILTER (saved as keep_study_10m.tif) for the study pixels:
{
  f_keep  <- paste0(out_dir, "keep_study_10m.tif")  # 10 m study pixels: 1 = flown leaf-off, 0 = flown leaf-on
  
  if (!file.exists(f_keep)) {
    # USA-NPN 30-year average first-leaf day of year, downloaded into out_dir from https://geoserver.usanpn.org/geoserver/wcs?service=WCS&version=2.0.1&request=GetCoverage&coverageId=si-x__30yr_avg_4k_leaf&format=geotiff
    fl <- rast(paste0(out_dir, "si-x__30yr_avg_4k_leaf.tif"))
    
    # Cropped to the LiDAR footprint, re-projected to the LiDAR CRS (small temporary file), then set to the 10 m grid (virtual)
    fl <- va(sources(project(crop(fl, ext(project(foot_p, crs(fl))) + 0.1), crs(p95_file), filename = tempfile(fileext = ".tif"))))
    
    # Keep test for a block of pixels (t0 = MinTime, t1 = MaxTime, f = first-leaf day, s = study mask, h = p95); only study pixels with flight dates are tested
    keep_fun <- function(t0, t1, f, s, h) {
      out <- rep(NA_real_, length(t0))
      i <- which(s == 1 & h >= h_scrub & !is.na(t0) & !is.na(t1))
      f <- f[i]
      f[is.na(f)] <- 121 # May 1 where the first-leaf layer has no value
      f <- f + leaf_lag  # spring cutoff: local first leaf + lag
      on0 <- m0[t0[i] %% 10000 %/% 100] + t0[i] %% 100 >= f & t0[i] %% 10000 < leafoff # MinTime leaf-on
      on1 <- m0[t1[i] %% 10000 %/% 100] + t1[i] %% 100 >= f & t1[i] %% 10000 < leafoff # MaxTime leaf-on
      out[i] <- as.numeric(!(on0 | on1))
      out
    }
    
    # One pass (study mask set to the 10 m grid)
    lapp(c(lr("MinTime"), lr("MaxTime"), fl, va(f_study), p95_file), keep_fun, filename = f_keep, wopt = list(datatype = "INT1U"))
  }
  
  keep <- rast(f_keep)
}



##################################################################
# MAIN PASS + 300 m COVARIATES, band by band: each band (band_rows = rows of 300 m cells run per loop) is read only over the columns with kept study pixels, its 10 m layers go to a temporary file, and its 300 m covariates are saved to out_dir/cov300_bands/ (saved bands are skipped, so a stopped run resumes)
# rumple index: mean over the cell's forest pixels of 1 / cos(slope of the 10 m p95 surface) = canopy surface area / flat area (1 = flat canopy)
{
  
  # Rows of 300 m cells per band (100 = 30 km)
  band_rows <- 100
  
  # Block of 10 m pixels; m = matrix with one row per pixel; columns: V bins (rows of vb), p95, keep, CHILI, distance to water, forest type, NumPoints, MinTime
  # m is not reshaped on purpose: app() first tries the function on single pixels and, if that works, calls it once per pixel (very slow);
  # a single pixel (plain vector) makes it fail at nrow(m), so app() calls it once per block of rows
  main_fun <- function(m) {
    
    p <- nrow(vb)
    out <- matrix(NA_real_, nrow(m), 9 + length(vl))
    
    # Kept study pixels (rows of the block); all other rows stay NA
    i <- which(m[, p + 2] == 1)
    if (length(i) == 0) return(out)
    m <- m[i, , drop = F]
    
    # % of all returns per height layer (column 1 = 0-1 m, then the vl layers) of each dataset; later datasets win where they overlap
    x <- matrix(NA, nrow(m), ncol(G))
    for (d in seq_along(lidar_src)) {
      xd <- m[, which(vb$ds == d), drop = F] %*% G[which(vb$ds == d), , drop = F]
      x <- ifelse(is.na(xd), x, xd)
    }
    v <- x[, -1, drop = F] # height layers above 1 m
    
    # FHD weighted by layer width: Shannon index of v over 1 m slices, each layer's returns spread evenly over its vw meters (same as -sum(p * log(p / vw))),
    # standardized by log(M), M = meters from 1 m to the top of the highest layer with returns (1 = returns spread evenly over every meter up to that top)
    tot <- rowSums(v)
    M <- cumsum(vw)[max.col((v > 0) * 1, ties.method = "last")]
    fhd <- ifelse(tot > 0 & M > 1, (log(tot) - rowSums(ifelse(v > 0, v * log(v), 0)) / tot + drop(v %*% log(vw)) / tot) / log(M), NA)
    
    # 1 for forest pixels, NA for scrub pixels (multiplies the forest-only layers)
    f <- ifelse(m[, p + 1] >= h_forest, 1, NA)
    
    out[i, ] <- cbind(
      # p95 (m)
      f * m[, p + 1], 
      # FHD
      f * fhd, 
      # NRD: share of all returns at 1-3 m (1-2 m layer + 2-3 m layer)
      f * (x[, 2] + x[, 3]),
      # point density (returns per m2 = NumPoints / 100 m2)
      f * m[, p + 6] / 100,
      # CHILI
      f * m[, p + 3], 
      # distance to water (m)
      f * m[, p + 4], 
      # TreeMap forest type (FIA code)
      f * ifelse(m[, p + 5] > 0, m[, p + 5], NA), 
      # flight year (YY of MinTime)
      f * (m[, p + 7] %/% 10000),
      # flight month-year (YYMM of MinTime)
      f * (m[, p + 7] %/% 100), 
      # % of all returns per height layer
      v
    ) 
    
    out
  }
  
  # Small processing blocks; large GDAL block cache (the JT files are stored in 512 x 512 tiles)
  terraOptions(memmax = 2)
  gdalCache(12000)
  
  # Bands not yet saved, then their mosaic (skipped once covariates_300m.tif exists)
  if (!file.exists(f_cov)) {
    # Main pass inputs: V bins set to the 10 m grid (virtual), p95, keep, GEE layers, NumPoints, MinTime
    x <- c(rast(lapply(vb$file, va)), p95_file, keep, gee, lr("NumPoints"), lr("MinTime"))
    
    # 300 m cells with any kept study pixel (1); same 300 m grid as the covariates
    work <- aggregate(keep, 30, fun = "max", na.rm = T)
    
    # Folder for the band files
    dir.create(paste0(out_dir, "cov300_bands"), showWarnings = F)
    
    for (b in seq(1, nrow(work), by = band_rows)) {
      # Band file, and the band's rows of 300 m cells
      f_band <- paste0(out_dir, "cov300_bands/band_", b, ".tif")
      n <- min(band_rows, nrow(work) - b + 1)
      
      # Columns of 300 m cells with kept study pixels in the band; skip the band if it is saved already or has none
      cols <- which(colSums(matrix(values(work, row = b, nrows = n), nrow = n, byrow = T) == 1, na.rm = T) > 0)
      if (file.exists(f_band) | length(cols) == 0) next
      
      # Window: the band's rows over the column range with kept study pixels, inside the 10 m grid (the previous band's window is removed first)
      e <- intersect(ext(work, cellFromRowCol(work, c(b, b + n - 1), range(cols))), ext(p95_file))
      window(x) <- NULL
      window(x) <- e
      
      # 10 m layers of the window (temporary file)
      f10 <- tempfile(fileext = ".tif")
      m10 <- app(x, main_fun, filename = f10, wopt = list(names = c("p95", "fhd", "nrd", "pd", "chili", "waterdist", "fortype", "year", "yymm", vl)))
      
      # Rumple index of each 10 m forest pixel: canopy surface area / flat area = 1 / cos(slope of the p95 surface, from its 3 x 3 neighborhood); 1 = flat canopy
      # p95 cropped to the window plus one pixel on each side (so edge pixels have all 8 neighbors), slope (radians) computed, then read over the window only
      f_rum <- c(tempfile(fileext = ".tif"), tempfile(fileext = ".tif"), tempfile(fileext = ".tif"))
      slp <- terrain(crop(p95_file, e + 10, filename = f_rum[1]), "slope", unit = "radians", filename = f_rum[2])
      window(slp) <- e
      rum10 <- lapp(c(slp, m10[["p95"]]), function(s, h) ifelse(is.na(h), NA, 1 / cos(s)), filename = f_rum[3])
      
      # 300 m covariates of the band, saved; temporary 10 m files deleted
      band <- c(aggregate(m10[[c("p95", "fhd", "nrd", "pd", "chili", "waterdist", vl)]], 30, fun = "mean", na.rm = T),
                aggregate(m10[[c("p95", "fhd")]], 30, fun = "sd", na.rm = T),
                aggregate(m10[[c("fortype", "year", "yymm")]], 30, fun = "modal", na.rm = T),
                aggregate(rum10, 30, fun = "mean", na.rm = T))
      names(band) <- c("p95", "fhd", "nrd", "pd", "chili", "waterdist", vl, "sdp95", "sdfhd", "fortype", "year", "yymm", "rumple")
      writeRaster(band, f_band)
      unlink(c(f10, f_rum))
      
      # Progress
      cat("300 m rows", b, "-", b + n - 1, "of", nrow(work), "done:", format(Sys.time(), "%H:%M"), "\n")
    }
    
    # 300 m covariates: mosaic of the band files (keeps the band layer names), saved
    writeRaster(merge(sprc(list.files(paste0(out_dir, "cov300_bands"), full.names = T))), f_cov)
  }
  cov300 <- rast(f_cov)
  
  # Summaries
  print(summary(cov300[[c("p95", "sdp95", "rumple", "fhd", "sdfhd", "nrd", "pd", "chili", "waterdist")]]))
  hist(cov300[[c("p95", "sdp95", "rumple", "fhd", "sdfhd", "nrd", "pd", "chili", "waterdist")]])
}



##################################################################
# FOREST TYPE GROUPS: 300 m cells per TreeMap forest type code, and per FIA forest type group (bar chart)
{
  # Get ordered table of forest-types and counts
  forest_counts <- freq(cov300[["fortype"]])
  forest_counts[order(-forest_counts$count),]
  
  # FIA forest type groups: each group code is the first forest type code of its range
  grp <- c("100" = "White/red/jack pine", "120" = "Spruce/fir", "140" = "Longleaf/slash pine", "150" = "Tropical softwoods",
           "160" = "Loblolly/shortleaf pine", "170" = "Other eastern softwoods", "180" = "Pinyon/juniper", "200" = "Douglas-fir",
           "220" = "Ponderosa pine", "240" = "Western white pine", "260" = "Fir/spruce/mountain hemlock", "280" = "Lodgepole pine",
           "300" = "Hemlock/Sitka spruce", "320" = "Western larch", "340" = "Redwood", "360" = "Other western softwoods",
           "370" = "California mixed conifer", "380" = "Exotic softwoods", "390" = "Other softwoods", "400" = "Oak/pine",
           "500" = "Oak/hickory", "600" = "Oak/gum/cypress", "700" = "Elm/ash/cottonwood", "800" = "Maple/beech/birch",
           "900" = "Aspen/birch", "910" = "Alder/maple", "920" = "Western oak", "940" = "Tanoak/laurel",
           "950" = "Other western hardwoods (retired)", "960" = "Other hardwoods", "970" = "Woodland hardwoods",
           "980" = "Tropical hardwoods", "990" = "Exotic hardwoods", "999" = "Nonstocked")
  
  # Group of each forest type code: the largest group code at or below it (e.g. 515 -> 500, 962 -> 960, 995 -> 990)
  forest_counts$group <- grp[findInterval(forest_counts$value, as.numeric(names(grp)))]
  
  # 300 m cells per group, sorted
  par(mar = c(4, 12, 2, 1))
  barplot(sort(tapply(forest_counts$count, forest_counts$group, sum)), horiz = T, las = 1, cex.names = 0.8,
          xlab = "300 m cells", main = "Forest type group (TreeMap FORTYPCD) of the 300 m cells")
  par(mar = c(5, 4, 4, 2) + 0.1)
}



##################################################################
# STOPOVER RESPONSES and model dataframes (saved to out_dir as <season>_d.csv)
# Each season keeps its own radar coverage; hot = 1 when a cell's VIR reaches the top-10% cutoff of its nearest station's full range in all 3 years of its LiDAR year +-1
{
  # Forest type table: every TreeMap forest type code (FIA FORTYPCD) in the 300 m cells, most common first, with its TreeMap color and name
  # (from the TreeMap FORTYPCD class table, " / " written as "/"); read.csv keeps the "#" of the colors (its comment.char is "")
  ft <- read.csv(text = "
code,col,name
801,#a2c9eb,Sugar maple/beech/yellow birch
503,#f22ab1,White oak/red oak/hickory
515,#eeafa3,Chestnut oak/black oak/scarlet oak
502,#bccc4b,Chestnut oak
506,#718640,Yellow-poplar/white oak/northern red oak
505,#77989d,Northern red oak
805,#6ab7f2,Hard maple/basswood
401,#ec42f6,Eastern white pine/northern red oak/white ash
103,#efdbcc,Eastern white pine
504,#f6e095,White oak
520,#cce5b9,Mixed upland hardwoods
708,#e4c3c0,Red maple/lowland
516,#9bd763,Cherry/white ash/yellow-poplar
510,#5fa7cc,Scarlet oak
519,#e88fbb,Red maple/oak
802,#3ff451,Black cherry
809,#b3714c,Red maple/upland
962,#f5f169,Other hardwoods
508,#c376e4,Sweetgum/yellow-poplar
702,#b1bef2,River birch/sycamore
608,#54c7ef,Sweetbay/swamp tupelo/red maple
105,#d0ce83,Eastern hemlock
163,#63f3ac,Virginia pine
167,#bba847,Pitch pine
405,#37999a,Virginia pine/southern red oak
409,#82eb3e,Other pine/hardwood
701,#cecceb,Black ash/American elm/red maple
606,#60dad1,Atlantic white cedar
517,#b838ee,Elm/ash/black locust
161,#97f2ad,Loblolly pine
705,#c4ec84,Sycamore/pecan/American elm
706,#efadec,Sugarberry/hackberry/elm/green ash
511,#9ae6e8,Yellow-poplar
104,#a8a9f2,Eastern white pine/eastern hemlock
512,#def3b1,Black walnut
901,#d28f25,Aspen
402,#7e9f81,Eastern redcedar/hardwood
501,#b8db98,Post oak/blackjack oak
102,#c765ec,Red pine
127,#f2c531,Northern white cedar
602,#c8ed2d,Sweetgum/Nuttall oak/willow oak
703,#f077ef,Cottonwood
513,#b88bf2,Black locust
507,#9d4f8d,Sassafras/persimmon
903,#dd82c7,Gray birch
999,#d5cc36,Nonstocked
707,#da23cf,Silver maple/American elm
171,#95eacd,Eastern redcedar
122,#9d86a6,White spruce
406,#ed54dd,Loblolly pine/hardwood
381,#e9e2eb,Scotch pine
509,#7cb133,Bur oak
165,#e9c991,Table Mountain pine
384,#eee1b3,Norway spruce
905,#e3f2e7,Pin cherry
404,#5cd76e,Shortleaf pine/oak
121,#47d0b6,Balsam fir
162,#d45549,Shortleaf pine
123,#a5f77a,Red spruce
125,#64e1f7,Black spruce
126,#afa9b0,Tamarack
995,#e67774,Other exotic hardwoods
704,#969aca,Willow
902,#f59550,Paper birch
385,#e4db79,Introduced larch
166,#ddbef2,Pond pine
601,#ed8a9c,Swamp chestnut oak/cherrybark oak
124,#dcf4d9,Red spruce/balsam fir
709,#bf90e1,Cottonwood/willow
514,#a5f031,Southern scrub oak
904,#c5f2a0,Balsam poplar
101,#6e26ec,Jack pine
605,#f0bd53,Overcup oak/water hickory
609,#8e6a31,Baldcypress/pondcypress
")
  
  # Forest types kept as their own class: at least 500 300 m cells (rows of ft stay in count order); all other types become "Other"
  ft <- ft[ft$code %in% subset(freq(cov300[["fortype"]]), count >= 500)$value, ]
  
  # Forest type class of each 300 m cell (class i = row i of ft; last class = Other); label = code + name
  ftype <- classify(cov300[["fortype"]], cbind(ft$code, 1:nrow(ft)), others = nrow(ft) + 1)
  levels(ftype) <- data.frame(value = 1:(nrow(ft) + 1), ftype = c(paste(ft$code, ft$name), "Other"))
  
  # Yearly VIR of the stations per season (one layer per year, 2013-2024; mean where station ranges overlap), resampled onto the 300 m grid
  yr_spring <- project(mosaic(sprc(lapply(fp(mos, "spring"), rast)), fun = "mean"), cov300, method = "bilinear")
  yr_autumn <- project(mosaic(sprc(lapply(fp(mos, "autumn"), rast)), fun = "mean"), cov300, method = "bilinear")
  
  # Radar stations from the NOAA NCEI station list
  # ICAO code, latitude, longitude, and site elevation (feet); stations in mos only (same order)
  nx <- read.fwf("https://www.ncei.noaa.gov/access/homr/file/nexrad-stations.txt", widths = c(-9, 4, -93, 9, -1, 10, -1, 6), skip = 2,  col.names = c("icao", "lat", "lon", "elev_ft"), strip.white = T)
  nx <- nx[match(mos, nx$icao), ]
  sta <- project(vect(nx, geom = c("lon", "lat"), crs = "EPSG:4326"), crs(cov300))
  
  # Distance (km) from each 300 m cell to each station (one layer per station), distance to the nearest station, the nearest station, and the site elevation (m) of that nearest station
  d_sta <- rast(lapply(1:nrow(sta), function(i) distance(rast(cov300[["p95"]]), sta[i]) / 1000))
  raddist <- min(d_sta)
  near <- which.min(d_sta)
  site <- classify(near, cbind(1:nrow(sta), nx$elev_ft * 0.3048))
  
  # Hot-year cutoffs of each season: 90th percentile of each station's VIR in each year over its full range (the station's own raster: uncropped and unmasked);
  # one layer per year, each 300 m cell holding the cutoff of its nearest station (made before the station codes are added to near, so classify() reads station numbers)
  cutoffs <- lapply(c(spring = "spring", autumn = "autumn"), function(season) {
    # 90th percentiles: rows = years (2013-2024), columns = stations (mos order)
    q90 <- sapply(fp(mos, season), function(f) apply(values(rast(f)), 2, quantile, 0.9, na.rm = T))
    rast(lapply(1:nrow(q90), function(i) classify(near, cbind(1:nrow(sta), q90[i, ]))))
  })
  
  # Nearest station as a class (station code), for the station random effect in the models
  levels(near) <- data.frame(value = 1:nrow(sta), station = mos)
  
  # Elevation (m): mean of the 30 m 3DEP pixels in each 300 m cell (GEE splits large exports into several files; all files are joined)
  elev <- resample(vrt(list.files(out_dir, pattern = "^gee_elev_30m.*\\.tif$", full.names = T)), cov300, method = "average")
  
  # Radar beam height above ground (m): nearest station's site elevation + height of its lowest beam (0.5 degrees) above the radar at that distance - cell elevation
  beamagl <- site + lapp(raddist * 1000, function(d) beam_height(beam_range(d, 0.5), 0.5, lat = mean(nx$lat))) - elev
  
  # Check: nearest radar station of each 300 m cell (sets the cell's hot-year cutoffs) and stations (red)
  plot(near, main = "Nearest radar station of each 300 m cell (colors) and stations (red)")
  plot(sta, add = T, col = "red", pch = 16)
  
  build <- function(yr, cutoffs, season) {
    
    # Each year's VIR kept only in cells where that year (layer i = 2012 + i) is within the LiDAR year +-1
    yr <- rast(lapply(1:nlyr(yr), function(i) ifel(abs(12 + i - cov300[["year"]]) <= 1, yr[[i]], NA)))
    
    # VIR and model covariates (incl. nearest station); model cells (ok) have all of them (this season's radar coverage only)
    base <- c(mean(yr, na.rm = T), cov300[[c("p95", "sdp95", "rumple", "fhd", "sdfhd", "nrd", "pd", "chili", "waterdist")]], raddist, elev, beamagl, cov300[["fortype"]], ftype, near)
    ok <- noNA(base, falseNA = T)
    
    # Hot years (one layer per year): 1 = VIR >= the nearest station's cutoff (top 10% of that station's full range), 0 = below, NA = no VIR or no cutoff that year
    hot <- yr >= cutoffs
    
    # Responses, covariates, hot in all 3 years (1 = hot in each year of the LiDAR year +-1, 0 = not), and years tested (years with VIR and a cutoff) of the model cells
    out <- mask(c(base, (sum(hot, na.rm = T) == 3) * 1, sum(!is.na(hot))), ok)
    
    names(out) <- c("vir", "p95", "sdp95", "rumple", "fhd", "sdfhd", "nrd", "pd", "chili", "waterdist", "raddist", "elev", "beamagl", "fortype", "ftype", "station", "hot", "nyr")
    
    # Model dataframe (one row per model cell), saved to out_dir
    write.csv(as.data.frame(out, xy = T, na.rm = T), paste0(out_dir, season, "_d.csv"), row.names = F)
    
    # Model cells as a raster
    out
    
  }
  
  spring <- build(yr_spring, cutoffs$spring, "spring")
  autumn <- build(yr_autumn, cutoffs$autumn, "autumn")
  
}



##################################################################
# MAPS of the spring model cells: responses and covariates (deciles as color breaks); mean % of returns per height layer (voxel map)
{
  
  # Mapview basemap options, getting rid of watermarks
  options("mapviewBasemaps" = c("Esri.WorldGrayCanvas", "CartoDB.DarkMatter", "OpenStreetMap", "Esri.WorldImagery", "OpenTopoMap"))
  
  # Map layer names (no "/" or "\": mapview uses each name in an HTML dependency name, which cannot contain them)
  lyr <- c(hot = "Hot in all 3 years (top 10% of the nearest station's full range)", p95 = "Canopy height p95 (m)", 
           sdp95 = "Canopy height SD (m)", rumple = "Rumple index (canopy surface per flat area)", fhd = "Vertical diversity (FHD)",
           sdfhd = "FHD SD (horizontal heterogeneity)", nrd = "Shrub layer (share of returns at 1-3 m)", 
           chili = "Heat load (CHILI)", 
           waterdist = "Distance to water (m)", pd = "Point density (returns per m2)", elev = "Elevation (m)")
  
  # Main map: spring VIR shown; other layers hidden (turn on in the layer menu); forest type classes as polygons so the legend shows their names
  map <- mapview(spring[["vir"]], layer.name = "Spring VIR", col.regions = turbo, 
                 at = unique(quantile(values(spring[["vir"]]), seq(0, 1, 0.1), na.rm = T)), 
                 na.color = "transparent", maxpixels = 1.4e6)
  
  for (v in names(lyr)) {
    map <- map + mapview(
      spring[[v]], layer.name = lyr[[v]], col.regions = turbo, 
      at = unique(quantile(values(spring[[v]]), seq(0, 1, 0.1), na.rm = T)), 
      na.color = "transparent", maxpixels = 1.4e6, hide = T)
  }
  
  # Forest type colors: TreeMap color of each class (ft, from STOPOVER RESPONSES) and grey for Other, named by class label
  fcol <- setNames(c(ft$col, "#bdbdbd"), c(paste(ft$code, ft$name), "Other"))
  
  # Forest type classes as polygons (one per class), rows sorted alphabetically so the map (colors used row by row) and the legend (colors given to classes in alphabetical order) agree
  ft_poly <- as.polygons(spring[["ftype"]])
  map <- map + mapview(
    ft_poly[order(ft_poly$ftype), ], zcol = "ftype", 
    col.regions = unname(fcol[sort(ft_poly$ftype)]), 
    alpha.regions = 1, layer.name = "Forest type class", 
    lwd = 0, hide = T)
  
  # LiDAR flight month of the spring model cells (most common MinTime month-year of the cell's forest pixels) as polygons labeled YYYY-MM; rows sorted so the colors (row by row) and the legend (alphabetical = chronological) agree
  fl_poly <- as.polygons(mask(cov300[["yymm"]], spring[["vir"]]))
  fl_poly$date <- sprintf("20%02d-%02d", fl_poly$yymm %/% 100, fl_poly$yymm %% 100)
  map <- map + mapview(
    fl_poly[order(fl_poly$date), ], zcol = "date", 
    col.regions = turbo(nrow(fl_poly)), 
    alpha.regions = 1, layer.name = "LiDAR flight month (YYYY-MM)", 
    lwd = 0, hide = T)
  print(map)
  
  # Voxel map: mean % of all returns per height layer (study pixels) in the spring model cells; one map layer per height layer
  vox <- mask(cov300[[vl]], spring[["vir"]])
  
  vlab <- sub("-Inf", "+", paste0("% of returns ", head(brk, -1), "-", brk[-1], " m"))
  
  vmap <- mapview(
    vox[[1]], layer.name = vlab[1], col.regions = viridis, 
    at = unique(quantile(values(vox[[1]]), seq(0, 1, 0.1), na.rm = T)), 
    na.color = "transparent", maxpixels = 1.4e6)
  
  for (k in 2:nlyr(vox)) {
    vmap <- vmap + mapview(
      vox[[k]], layer.name = vlab[k], col.regions = viridis, 
      at = unique(quantile(values(vox[[k]]), seq(0, 1, 0.1), na.rm = T)), 
      na.color = "transparent", maxpixels = 1.4e6, hide = T)
  }
  
  print(vmap)
}



##################################################################
# VOXEL RIDGELINES BY FOREST TYPE: mean voxel profile (% of all returns per m of height) of the model cells in each forest type class
# x = height above ground, with each height layer drawn over its true height range (1-2 ... 5-6 m are 1 m wide; 6-18, 18-30, 30-42 m are 12 m wide; the open 42+ m layer is drawn as 42-54 m)
# Each ridge is a step line: its height over a layer = the class's mean % of returns in that layer divided by the layer's width (% of returns per m of height); the 0-1 m layer is not kept in cov300, so profiles start at 1 m
# Voxel layers = mean of the cell's kept study pixels (forest and scrub); spring model cells
{
  # Forest type class, canopy height (p95, forest pixels), and % of returns per height layer of the model cells; classes in FIA code order (lowest code at the top, Other at the bottom)
  vox_d <- as.data.frame(c(spring[["ftype"]], spring[["p95"]], cov300[[vl]]), na.rm = T)
  vox_d$ftype <- factor(as.character(vox_d$ftype), levels = rev(sort(unique(as.character(vox_d$ftype)))))
  
  # Bottom and top (m) of each height layer; the open top layer (42+ m) drawn to 54 m
  bot <- head(brk, -1)
  top <- pmin(brk[-1], 54)
  
  # Mean % of returns per height layer in each forest type class
  m_ft <- aggregate(vox_d[, vl], list(class = vox_d$ftype), mean)
  
  # Ridge points: two per class and layer (the layer's bottom and top), both at the layer's % of returns per m (layer value / layer width; 42+ m divided by its drawn width, 12 m)
  pr <- data.frame()
  for (i in 1:nrow(m_ft)) {
    for (j in seq_along(vl)) {
      pr <- rbind(pr, data.frame(class = m_ft$class[i], x = c(bot[j], top[j]), value = m_ft[i, vl[j]] / (top[j] - bot[j])))
    }
  }
  
  # Ridgelines (tallest step = 1.5 row spacings); colors = forest type colors (fcol, from the MAPS chunk)
  g <- ggplot(pr, aes(x = x, y = class, height = value, fill = class)) +
    geom_ridgeline(scale = 1.5 / max(pr$value), alpha = 0.8) +
    scale_fill_manual(values = fcol, guide = "none") +
    scale_x_continuous(breaks = c(1, 6, 18, 30, 42, 54)) +
    labs(x = "Height above ground (m)", y = NULL, title = "Mean voxel profile (% of all returns per m of height) by forest type")
  print(g)
}



##################################################################
# VOXEL RIDGELINES BY CANOPY HEIGHT: the same profiles by canopy height class (p95 of the cell's forest pixels, 5 m classes; lowest class at the bottom)
{
  # Canopy height class of each model cell
  vox_d$hclass <- cut(vox_d$p95, c(5, 10, 15, 20, 25, 30, Inf), right = F, labels = c("5-10 m", "10-15 m", "15-20 m", "20-25 m", "25-30 m", "30+ m"))
  
  # Mean % of returns per height layer in each canopy height class
  m_h <- aggregate(vox_d[, vl], list(class = vox_d$hclass), mean)
  
  # Ridge points: two per class and layer (the layer's bottom and top), both at the layer's % of returns per m (layer value / layer width)
  pr <- data.frame()
  for (i in 1:nrow(m_h)) {
    for (j in seq_along(vl)) {
      pr <- rbind(pr, data.frame(class = m_h$class[i], x = c(bot[j], top[j]), value = m_h[i, vl[j]] / (top[j] - bot[j])))
    }
  }
  
  # Ridgelines (tallest step = 1.5 row spacings); viridis colors from short (dark) to tall (yellow)
  g <- ggplot(pr, aes(x = x, y = class, height = value, fill = class)) +
    geom_ridgeline(scale = 1.5 / max(pr$value), alpha = 0.8) +
    scale_fill_manual(values = viridis(nlevels(vox_d$hclass)), guide = "none") +
    scale_x_continuous(breaks = c(1, 6, 18, 30, 42, 54)) +
    labs(x = "Height above ground (m)", y = "Canopy height (p95)", title = "Mean voxel profile (% of all returns per m of height) by canopy height")
  print(g)
}