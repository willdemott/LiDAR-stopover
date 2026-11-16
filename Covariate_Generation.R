# Stopover intensity vs LiDAR-derived forest structure: DATA PREP, 300 m COVARIATES and RESPONSES (run order 3 of 5)
# Each covariate (and each season's responses) is its own 300 m raster on one shared grid (EPSG:6350, 300 m cells on the LiDAR 10 m grid origin), cached in ./cov_dir as <name>.tif
# All cached layers are written as ONE stacked raster, stack_300m.tif in out_dir (with its .aux.xml, which holds the habitat class and station labels): the only input of stopover_lidar_responses.R and stopover_lidar_model.R
# Needs: lidar_300m.tif, stations.csv (stopover_lidar_gee_prep.R), veg_10m.tif (stopover_lidar_mask.R), and the GEE exports (NE_gee/; every layer but woody already as 300 m cell values)
# LiDAR covariates only contain non-leaf forest pixels of each 300 m cell (woody, flown leaf-off, p95 >= h_forest); the habitat class (ftype) also uses the scrub pixels (p95 1-5 m)
# Each cell's station = the radar whose VIR raster covers it (the nearer one where two overlap): its VIR, distance, beam height, and hot cutoff all come from that radar

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
  
  # Parameters: height layer breaks (m), shared by both datasets' V bins; seasons and radar years (layer 1 = 2013) of the responses;
  # largest yearly VIR kept (higher values set to NA before the hot cutoffs, the hot test, and the mean VIR); % of a cell's kept pixels that must be scrub for the Shrub habitat class (more than)
  brk <- c(1, 2, 3, 4, 5, 6, 18, 30, 42, Inf)
  seasons <- c("spring", "autumn")
  yrs <- 2013:2024
  vir_cap <- 7000
  shrub_min <- 50
  
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
  
  # GEE 300 m layer, e.g. g300("elev"): tiles of NE_gee/gee_<name>_300m*.tif stitched on the shared 300 m grid (the GEE cells are the grid's cells, so nothing is aggregated or resampled)
  g300 <- function(name) {
    vrt(list.files(gee_dir, pattern = paste0("^gee_", name, "_300m.*\\.tif$"), full.names = T),
        options = c("-te", xmin(grid300), ymin(grid300), xmax(grid300), ymax(grid300), "-tr", res(grid300)))
  }
  
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
# HEAT LOAD (chili.tif): mean CHILI (0 = cool to 255 = warm) of each 300 m cell's 10 m pixels (GEE, US 10 m CHILI)
{
  if (!file.exists(cf("chili"))) {
    g300("chili") %>% setNames("chili") %>% writeRaster(cf("chili"))
  }
}



##################################################################
# DISTANCE TO WATER (waterdist.tif): mean distance (m) of each 300 m cell's 10 m pixels to the nearest water (GEE: NHD and WorldCover water)
if (!file.exists(cf("waterdist"))) {
  vrt(list.files(gee_dir, pattern = "^gee_waterdist_30m.*\\.tif$", full.names = T),   # 30 m GEE tiles (means of 10 m distances) on the shared grid's extent
      options = c("-te", xmin(grid300), ymin(grid300), xmax(grid300), ymax(grid300), "-tr", 30, 30)) %>%
    aggregate(10, fun = "mean", na.rm = T) %>%   # 10 x 10 30 m pixels = one 300 m cell
    setNames("waterdist") %>% writeRaster(cf("waterdist"))
}



##################################################################
# DISTANCE TO FOREST EDGE (edgedist.tif): mean distance (m) of each 300 m cell's 10 m pixels to the nearest opening
# (GEE: non-woody WorldCover land other than water, in patches of at least 0.5 ha; woody = tree cover or shrubland)
{
  if (!file.exists(cf("edgedist"))) {
    vrt(list.files(gee_dir, pattern = "^gee_edgedist_30m.*\\.tif$", full.names = T),   # 30 m GEE tiles (means of 10 m distances) on the shared grid's extent
        options = c("-te", xmin(grid300), ymin(grid300), xmax(grid300), ymax(grid300), "-tr", 30, 30)) %>%
      aggregate(10, fun = "mean", na.rm = T) %>%   # 10 x 10 30 m pixels = one 300 m cell
      setNames("edgedist") %>% writeRaster(cf("edgedist"))
  }
}



##################################################################
# ELEVATION (elev.tif): mean elevation (m) of each 300 m cell's 10 m pixels (GEE, 3DEP 10 m)
{
  if (!file.exists(cf("elev"))) {
    g300("elev") %>% setNames("elev") %>% writeRaster(cf("elev"))
  }
}



##################################################################
# FOREST TYPE GROUP (fortype.tif): FIA forest type group code (e.g. 500 = oak/hickory) most common among each 300 m cell's TreeMap forest pixels (GEE, TreeMap 2020); NA = no forest pixels
{
  if (!file.exists(cf("fortype"))) {
    g300("fortype") %>% setNames("fortype") %>% writeRaster(cf("fortype"), datatype = "INT2S")
  }
}



##################################################################
# SCRUB SHARE (scrub.tif): % of each 300 m cell's kept pixels that are scrub (veg = 1: h_scrub <= p95 < h_forest); used for the Shrub habitat class
{
  if (!file.exists(cf("scrub"))) {
    tf <- tempfile(fileext = ".tif")
    lapp(veg, function(v) ifelse(v == 1, 100, 0), filename = tf) %>%   # 100 = scrub, 0 = forest, NA = not kept (10 m, temporary file)
      aggregate(30, fun = "mean", na.rm = T) %>%
      setNames("scrub") %>% writeRaster(cf("scrub"))
    unlink(tf)
  }
}



##################################################################
# HABITAT CLASS (ftype.tif): Shrub when more than shrub_min % of a cell's kept pixels are scrub; otherwise its FIA forest type group
# Groups with at least 500 non-shrub 300 m cells of forest LiDAR keep their own class (most common first); all other groups = Other; NA = no TreeMap forest pixels and not Shrub
{
  if (!file.exists(cf("ftype"))) {
    # FIA forest type groups: each group code is the first forest type code of its range
    grp <- c("100" = "White/red/jack pine", "120" = "Spruce/fir", "140" = "Longleaf/slash pine", "150" = "Tropical softwoods", "160" = "Loblolly/shortleaf pine", "170" = "Other eastern softwoods", "180" = "Pinyon/juniper", "200" = "Douglas-fir", "220" = "Ponderosa pine", "240" = "Western white pine", "260" = "Fir/spruce/mountain hemlock", "280" = "Lodgepole pine", "300" = "Hemlock/Sitka spruce", "320" = "Western larch", "340" = "Redwood", "360" = "Other western softwoods", "370" = "California mixed conifer", "380" = "Exotic softwoods", "390" = "Other softwoods", "400" = "Oak/pine", "500" = "Oak/hickory", "600" = "Oak/gum/cypress", "700" = "Elm/ash/cottonwood", "800" = "Maple/beech/birch", "900" = "Aspen/birch", "910" = "Alder/maple", "920" = "Western oak", "940" = "Tanoak/laurel", "950" = "Other western hardwoods (retired)", "960" = "Other hardwoods", "970" = "Woodland hardwoods", "980" = "Tropical hardwoods", "990" = "Exotic hardwoods", "999" = "Nonstocked")
    
    # Shrub cells (1 = Shrub, 0 = not, NA = no kept pixels), and the forest type group of every other cell
    shrub <- rast(cf("scrub")) > shrub_min
    fgrp <- mask(rast(cf("fortype")), shrub, maskvalues = 1)
    
    # 300 m cells with forest LiDAR (p95) per group, not counting Shrub cells -> groups with at least 500 cells, most common first
    ft <- freq(mask(fgrp, rast(cf("p95")))) %>%
      filter(count >= 500) %>%
      arrange(desc(count))
    
    # Class of each 300 m cell: class i = row i of ft, then Other, then Shrub; labeled group code + name
    # (mask: cells without TreeMap forest pixels stay NA; terra versions differ in whether classify's 'others' also fills NA cells)
    ftype <- classify(fgrp, cbind(ft$value, 1:nrow(ft)), others = nrow(ft) + 1) %>% mask(fgrp)
    ftype <- ifel(shrub == 1, nrow(ft) + 2, ftype)
    levels(ftype) <- data.frame(value = 1:(nrow(ft) + 2), ftype = c(paste(ft$value, grp[as.character(ft$value)]), "Other", "Shrub (LiDAR p95 1-5 m)"))
    ftype %>% setNames("ftype") %>% writeRaster(cf("ftype"))
    
    # Check: 300 m cells with forest LiDAR per habitat class
    print(freq(mask(ftype, rast(cf("p95")))))
  }
}



##################################################################
# RADAR STATION AND DISTANCE (station.tif, raddist.tif): id (stations.csv) of the station whose VIR raster covers each 300 m cell, and the distance (km) to it
# A station covers the cells where its raster has VIR in any year of either season (bilinear, as in RESPONSES); where two stations cover a cell, the nearer one
{
  if (!all(file.exists(cf(c("station", "raddist"))))) {
    d <- lapply(sta$id, function(i) {
      spr <- max(rast(fp(sta$icao[i], "spring")), na.rm = T) %>% project(grid300, method = "bilinear")   # spring VIR in any year
      aut <- max(rast(fp(sta$icao[i], "autumn")), na.rm = T) %>% project(grid300, method = "bilinear")   # autumn VIR in any year
      mask(distance(rast(grid300), sta_v[i]), max(spr, aut, na.rm = T))   # distance (m) to station i, inside its coverage only
    }) %>%
      rast()   # one layer per station (layer i = station id i)
    which.min(d) %>% setNames("station") %>% writeRaster(cf("station"), overwrite = T)   # layer of the nearest covering station = its id
    (min(d, na.rm = T) / 1000) %>% setNames("raddist") %>% writeRaster(cf("raddist"), overwrite = T)
  }
}



##################################################################
# RADAR BEAM HEIGHT (beamagl.tif): height (m) of the cell's station's lowest beam (0.5 degrees) above the ground
# = station site elevation + beam height above the radar at that distance - cell elevation
{
  if (!file.exists(cf("beamagl"))) {
    site <- classify(rast(cf("station")), cbind(sta$id, sta$elev_m)) # site elevation (m)
    beam <- lapp(rast(cf("raddist")) * 1000, function(d) beam_height(beam_range(d, 0.5), 0.5, lat = mean(sta$lat)))  # beam height above the radar (m)
    (site + beam - rast(cf("elev"))) %>% setNames("beamagl") %>% writeRaster(cf("beamagl"))
  }
}



##################################################################
# HOT-YEAR CUTOFFS (saved as cutoffs.csv in cov_dir): 90th percentile of each station's VIR in each year and season, over its full range (the station's own raster: uncropped and unmasked), yearly VIR above vir_cap removed first
{
  f_cut <- paste0(cov_dir, "cutoffs.csv")
  if (!file.exists(f_cut)) {
    lapply(seasons, function(s) {
      lapply(sta$id, function(i) {
        v <- values(rast(fp(sta$icao[i], s)))   # one column per year (layer)
        v[v > vir_cap] <- NA                     # yearly VIR above vir_cap removed
        q90 <- apply(v, 2, quantile, 0.9, na.rm = T)
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
# hot = 1 if the cell's VIR reached its station's cutoff in all 3 of those years; nyr = years tested (VIR and a cutoff), so hot = 1 needs nyr = 3
# Yearly VIR above vir_cap is removed first, so it counts neither for hot nor for the mean (a removed year lowers nyr)
{
  # Radar years within the LiDAR year +-1 of each 300 m cell: one layer per radar year (1 = within, 0 = not)
  win <- rast(lapply(yrs, function(y) abs(y - 2000 - rast(cf("year"))) <= 1))
  
  for (s in seasons) {
    if (!file.exists(cf(paste0("resp_", s)))) {
      # Yearly VIR on the 300 m grid (one layer per year): each cell from its own station's raster (station.tif), values above vir_cap removed, bilinear; kept in the LiDAR window only
      vir <- lapply(sta$id, function(i) {
        rast(fp(sta$icao[i], s)) %>%
          clamp(upper = vir_cap, values = F) %>%                     # yearly VIR above vir_cap -> NA
          project(grid300, method = "bilinear") %>%
          mask(rast(cf("station")) == i, maskvalues = c(0, NA))     # cells of station i only
      }) %>%
        sprc() %>% merge() %>%   # stations joined (their cells do not overlap)
        mask(win, maskvalues = c(0, NA))
      
      # Hot-year cutoff of each cell's station: one layer per year
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