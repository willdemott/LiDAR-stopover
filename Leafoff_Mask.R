# Stopover intensity vs LiDAR-derived forest structure: LEAF-OFF MASK (run order 2 of 5)
# Makes the final mask used in data prep: woody 300 m cells (woody_cells_300m.tif) and the kept 10 m pixels behind every LiDAR covariate (veg_10m.tif)
# The green-up figures (flight dates, leaf-out lag) justify the leaf-on window and run before the slow KEPT PIXELS pass, so leaf_lag can be settled first (a changed leaf_lag needs veg_10m.tif deleted)
# Needs: lidar_300m.tif and radar_300m.tif (stopover_lidar_gee_prep.R), the GEE woody tiles (NE_gee/), and the USA-NPN first-leaf layer (si-x__30yr_avg_4k_leaf.tif)

###################################################################
# Load Packages
{
  library(terra)
  library(dplyr)
}



##################################################################
# SETUP: folders, parameters, and the 10 m layers on the LiDAR grid
{
  # Folders: inputs and outputs (LiDAR, first-leaf layer, masks), and the GEE exports
  out_dir <- "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/"
  gee_dir <- paste0(out_dir, "NE_gee/")
  
  # LiDAR datasets (folder + file-name prefix); where they overlap, the last one listed is used
  lidar_src <- c("C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/JT_LiDAR/Northeast_epsg6350_",
                 "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/Burch_LiDAR/Burch/mosaic_epsg6350_")
  
  # Parameters
  woody_min <- 75    # % of a 300 m cell's 10 m pixels classed woody (WorldCover tree cover or shrubland)
  h_scrub   <- 1     # p95 (m) from which a pixel counts as vegetation (scrub)
  h_forest  <- 5     # p95 (m) from which a pixel counts as forest
  leaf_lag  <- 14    # days after the local first-leaf day from which a flight counts as leaf-on
  leafoff   <- 1001  # MMDD from which a flight counts as leaf-off again (Oct 1)
  
  # Processing: small blocks, large GDAL cache (the JT files are stored in 512 x 512 tiles)
  terraOptions(memmax = 2)
  gdalCache(12000)
  
  # LiDAR metric v of both datasets as one 10 m virtual mosaic, e.g. lr("p95"); p95 is the reference grid
  lr <- function(v) vrt(paste0(lidar_src, v, "_10m.tif"))
  p95 <- lr("p95")
  
  # Any file(s) set to the 10 m LiDAR grid (virtual): tiles are stitched and every 10 m layer lines up with p95
  va <- function(f) vrt(f, options = c("-te", xmin(p95), ymin(p95), xmax(p95), ymax(p95), "-tr", res(p95)))
  
  # LiDAR footprint (bounding box) for projecting
  foot <- as.polygons(ext(p95), crs = crs(p95))
  
  # 10 m woody pixels (100 = WorldCover tree cover or shrubland, 0 = other): GEE tiles stitched on the 10 m grid
  woody10 <- va(list.files(gee_dir, pattern = "^gee_woody_10m.*\\.tif$", full.names = T))
  
  # Local first-leaf day of year (USA-NPN 30-year average, ~4 km): cropped to the LiDAR footprint, projected to the LiDAR CRS (small file), set to the 10 m grid
  fl <- rast(paste0(out_dir, "si-x__30yr_avg_4k_leaf.tif")) %>%
    crop(ext(project(foot, crs(.))) + 0.1) %>%
    project(crs(p95), filename = tempfile(fileext = ".tif")) %>%
    sources() %>%
    va()
  
  # Day of year of a flight date (YYMMDD)
  doy <- function(t) c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)[t %% 10000 %/% 100] + t %% 100
}



##################################################################
# WOODY CELLS (saved as woody_cells_300m.tif in out_dir): 1 = 300 m cell inside radar coverage with at least woody_min % woody 10 m pixels
{
  f_cell <- paste0(out_dir, "woody_cells_300m.tif")
  if (!file.exists(f_cell)) {
    # % woody of each 300 m cell (mean of the 0/100 pixels; 75% = at least 675 of the 900 pixels), kept inside radar coverage (radar_300m.tif)
    w300 <- aggregate(woody10, 30, fun = "mean", na.rm = T)
    ifel(w300 >= woody_min & rast(paste0(out_dir, "radar_300m.tif")) == 1, 1, NA, filename = f_cell, datatype = "INT1U")
  }
}



##################################################################
# GREEN-UP SAMPLE (used by the two figures below): regular sample of 10 m pixels, as a small raster
# Kept where the leaf-on test of KEPT PIXELS applies: woody pixel, inside a woody 300 m cell, p95 >= h_scrub
{
  s <- c(p95, woody10, va(f_cell), fl, lr("MinTime"), lr("MaxTime"), lr("V2"), lr("V3"), lr("V4"), lr("V5"), lr("V6")) %>%
    setNames(c("p95", "woody", "cell", "fl", "MinTime", "MaxTime", "V2", "V3", "V4", "V5", "V6")) %>%
    spatSample(1e6, method = "regular", as.raster = T)
  s <- mask(s, s[["woody"]] == 100 & s[["cell"]] == 1 & s[["p95"]] >= h_scrub, maskvalues = c(0, NA))
}



##################################################################
# FLIGHT DATES of the sampled pixels: histogram of flight dates and map of flight month (first flight, MinTime)
{
  # Histogram of flight dates by week; every date set in 2000 so month/day line up across years (220415 -> Apr 15)
  hist(as.POSIXct(sprintf("2000%04d", as.integer(values(s[["MinTime"]] %% 10000, na.rm = T))), format = "%Y%m%d", tz = "UTC"),
       breaks = "weeks", freq = T, format = "%b %d", main = "Flight dates of the sampled pixels (MinTime)", xlab = "")
  
  # Map of flight month (MMDD classes)
  plot(s[["MinTime"]] %% 10000, type = "interval", breaks = c(1:12 * 100 + 0.5, 1231.5),
       col = c("#08306B", "#2171B5", "#9ECAE1", "#A1D99B", "#41AB5D", "#006D2C", "#00441B", "#006D2C", "#41AB5D", "#DD8D3C", "#6BAED6", "#4292C6"),
       plg = list(legend = month.abb), main = "Flight month of the sampled pixels (MinTime)")
}



##################################################################
# LEAF-OUT LAG FIGURE: why spring flights up to the local first-leaf day + leaf_lag are kept; forest pixels (p95 >= h_forest) of the green-up sample
# A: % of all returns at 1-6 m vs days from the local first-leaf day to the pixel's last flight (MaxTime: the date the spring leaf-on test turns on); pixels flown over more than 7 days left out
# B: % of forest pixels dropped as leaf-on for lags of 0-35 days (same test as KEPT PIXELS)
{
  # Forest pixels of the sample with both flight dates: first-leaf day (May 1 = day 121 where missing, as in KEPT PIXELS), day of year of the first and last flight,
  # days from first leaf to the last flight (negative = before first leaf), and % of all returns at 1-6 m
  x <- as.data.frame(s, na.rm = F) %>%
    filter(p95 >= h_forest, !is.na(MinTime), !is.na(MaxTime)) %>%
    mutate(fl = ifelse(is.na(fl), 121, fl),
           doy0 = doy(MinTime),
           doy1 = doy(MaxTime),
           d = doy1 - fl,
           r16 = V2 + V3 + V4 + V5 + V6)
  
  # Panel A pixels: last flight 4 weeks before to 6 weeks after first leaf, all flights within 7 days of the same year; their week from first leaf
  a <- x %>%
    filter(d >= -28, d < 42, MaxTime %/% 10000 == MinTime %/% 10000, doy1 - doy0 <= 7, !is.na(r16)) %>%
    mutate(wk = cut(d, seq(-28, 42, 7), right = F))
  
  # Per week: forest pixels, flight dates, and mean % of returns at 1-6 m (few flight dates in a week = a rise may be one collection, not leaf-out)
  wm <- a %>%
    group_by(wk, .drop = F) %>%
    summarise(pixels = n(), flight_dates = n_distinct(MaxTime), mean_pct_1_6m = round(mean(r16), 2))
  print(wm)
  
  # Per flight date (MaxTime): mean days from first leaf, mean % of returns at 1-6 m, and forest pixels
  fd <- a %>%
    group_by(MaxTime) %>%
    summarise(d = mean(d), r16 = mean(r16), n = n())
  
  # Panel B: % of forest pixels dropped as leaf-on for lags of 0-35 days: MinTime or MaxTime on/after first leaf + lag and before leafoff
  dropped <- sapply(setNames(0:35, 0:35), function(L) {
    100 * mean((x$doy0 >= x$fl + L & x$MinTime %% 10000 < leafoff) | (x$doy1 >= x$fl + L & x$MaxTime %% 10000 < leafoff))
  })
  print(round(dropped[c("0", "7", "14", "21", "28")], 1))
  
  # Figure: two panels side by side
  par(mfrow = c(1, 2))
  
  # A: flight dates (gray, size ~ forest pixels) and weekly means (black); dotted line = first leaf; dashed line and shading = dropped as leaf-on
  plot(fd$d, fd$r16, type = "n", xlim = c(-28, 42), xlab = "Days from local first leaf to last flight", ylab = "% of returns at 1-6 m (forest pixels)", main = "A. Lower-canopy returns")
  rect(leaf_lag, par("usr")[3], par("usr")[2], par("usr")[4], col = adjustcolor("orange", 0.2), border = NA)   # shading from leaf_lag to the right edge of the plot
  abline(v = c(0, leaf_lag), lty = c(3, 2))
  points(fd$d, fd$r16, pch = 16, col = "grey60", cex = 0.5 + 2 * sqrt(fd$n / max(fd$n)))
  lines(seq(-25, 38, 7), wm$mean_pct_1_6m, type = "b", pch = 15, lwd = 2)
  legend("topleft", c("Flight date (size ~ forest pixels)", "Weekly mean", "Dropped as leaf-on"), pch = c(16, 15, 15), pt.cex = c(1, 1, 2),
         col = c("grey60", "black", adjustcolor("orange", 0.4)), lty = c(NA, 1, NA), bty = "n", cex = 0.8)
  
  # B: % of forest pixels dropped vs lag; point and dashed line = leaf_lag
  plot(0:35, dropped, type = "l", lwd = 2, xlab = "Days added to local first leaf (lag)", ylab = "% of forest pixels dropped as leaf-on", main = "B. Data lost to the leaf-on filter")
  abline(v = leaf_lag, lty = 2)
  points(leaf_lag, dropped[as.character(leaf_lag)], pch = 16)
  text(leaf_lag, dropped[as.character(leaf_lag)], paste0(round(dropped[as.character(leaf_lag)], 1), "%"), pos = 4)
  
  par(mfrow = c(1, 1))
}



##################################################################
# KEPT PIXELS = the final mask (saved as veg_10m.tif in out_dir): the 10 m pixels behind the LiDAR covariates; 2 = forest (p95 >= h_forest), 1 = scrub (h_scrub <= p95 < h_forest), NA = not used
# Kept = woody pixel, inside a woody 300 m cell, p95 >= h_scrub, and both flights (MinTime, MaxTime) leaf-off; leaf-on = on/after the local first-leaf day + leaf_lag and before leafoff
{
  f_veg <- paste0(out_dir, "veg_10m.tif")
  if (!file.exists(f_veg)) {
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
}