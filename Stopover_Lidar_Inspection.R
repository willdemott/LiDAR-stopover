# Stopover intensity vs LiDAR-derived forest structure: RESPONSES and MODEL DATAFRAMES
# Reads the 300 m covariates made by stopover_lidar_covariates.R (cov300/), adds the radar responses on the same grid, and writes spring_d.csv and autumn_d.csv
# vir = mean VIR of the radar years within the cell's LiDAR year +-1; hot = 1 if the cell's VIR reached the top 10% of its nearest station's full range in all 3 of those years
# A chunk skips its work when its file exists: delete a file to rebuild it
# Pipes (%>%) read left to right: take a raster or table, do one step, pass the result to the next step

###################################################################
# Load Packages
{
  library(terra)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
  library(mapview)
  library(viridisLite)
}



##################################################################
# SETUP: folders, radar files, years, model covariates, and the 300 m covariates
{
  # Folders: radar rasters, inputs and model dataframes, the 300 m covariates, and the 300 m responses
  radar_dir <- "C:/Users/willd/Box/"
  out_dir <- "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/"
  cov_dir <- paste0(out_dir, "cov300/")
  resp_dir <- paste0(out_dir, "resp300/")
  dir.create(resp_dir, showWarnings = F)
  
  # File of a 300 m covariate, e.g. cf("p95"), and of a season's 300 m responses, e.g. rf("spring")
  cf <- function(name) paste0(cov_dir, name, ".tif")
  rf <- function(season) paste0(resp_dir, season, ".tif")
  
  # Yearly VIR raster (one layer per year) of station s and season
  fp <- function(s, season) sprintf("%sDFF_ppi_years/%s_2013-24_%s_bioRad_ppi_years_AGL_inverse_Refit.tif", radar_dir, s, season)
  
  # Seasons, radar years (layer 1 = 2013), and height layer breaks (m) of the voxel layers (as in the covariate script)
  seasons <- c("spring", "autumn")
  yrs <- 2013:2024
  brk <- c(1, 2, 3, 4, 5, 6, 18, 30, 42, Inf)
  
  # Covariates in the model dataframes (year and yymm = LiDAR flight year and month, YY and YYMM; not model covariates)
  vars <- c("p95", "sdp95", "rumple", "fhd", "sdfhd", "nrd", "pd", "chili", "waterdist", "raddist", "elev", "beamagl", "fortype", "ftype", "station", "year", "yymm")
  
  # 300 m covariates as one stack (layer names = file names; voxel.tif adds one layer per height layer), stations (id = value in station.tif), and voxel layer names
  cov300 <- rast(list.files(cov_dir, pattern = "\\.tif$", full.names = T))
  sta <- read.csv(paste0(cov_dir, "stations.csv"))
  vl <- names(rast(cf("voxel")))
}



##################################################################
# HOT-YEAR CUTOFFS (saved as resp300/cutoffs.csv): 90th percentile of each station's VIR in each year and season, over its full range (the station's own raster: uncropped and unmasked)
{
  f_cut <- paste0(resp_dir, "cutoffs.csv")
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
# RESPONSES (saved as resp300/<season>.tif): vir, hot, and nyr of each 300 m cell
# Only the radar years within the cell's LiDAR year +-1 count; nyr = years tested (VIR and a cutoff), so hot = 1 needs nyr = 3
{
  # Radar years within the LiDAR year +-1 of each 300 m cell: one layer per radar year (1 = within, 0 = not)
  win <- rast(lapply(yrs, function(y) abs(y - 2000 - cov300[["year"]]) <= 1))
  
  for (s in seasons) {
    if (!file.exists(rf(s))) {
      # Yearly VIR on the 300 m grid (one layer per year; mean where station ranges overlap; bilinear), kept in the LiDAR window only
      vir <- fp(sta$icao, s) %>%
        lapply(rast) %>% sprc() %>% mosaic(fun = "mean") %>%
        project(cov300[["p95"]], method = "bilinear") %>%
        mask(win, maskvalues = c(0, NA))
      
      # Hot-year cutoff of each cell's nearest station: one layer per year
      cut_r <- rast(lapply(yrs, function(y) {
        rcl <- cutoffs %>% filter(season == s, year == y) %>% select(id, q90) %>% as.matrix()   # station id -> cutoff
        classify(cov300[["station"]], rcl)
      }))
      
      # Hot years: 1 = VIR >= cutoff, 0 = below, NA = no VIR or no cutoff that year
      hot <- vir >= cut_r
      
      # Responses: mean VIR of the window years, hot in all 3 window years (1/0), years tested
      c(mean(vir, na.rm = T), (sum(hot, na.rm = T) == 3) * 1, sum(!is.na(hot))) %>%
        setNames(c("vir", "hot", "nyr")) %>%
        writeRaster(rf(s))
    }
  }
}



##################################################################
# MODEL DATAFRAMES (saved to out_dir as <season>_d.csv): one row per model cell = 300 m cell with every model covariate and response present
{
  for (s in seasons) {
    d <- c(cov300[[vars]], rast(rf(s))) %>%
      as.data.frame(xy = T, na.rm = T) %>%                         # model cells
      rename(id = station) %>%
      left_join(select(sta, id, station = icao), by = "id") %>%    # nearest station's ICAO code
      select(-id)
    write.csv(d, paste0(out_dir, s, "_d.csv"), row.names = F)
    
    # Check: model cells, cells with all 3 years tested, and cells hot in all 3 years
    cat(s, ":", nrow(d), "model cells;", sum(d$nyr == 3), "with 3 years tested;", sum(d$hot), "hot in all 3 years\n")
  }
}



##################################################################
# LIDAR TIMING BY STATION: how much the LiDAR flight dates of the model cells vary within each radar station's area
# Each cell's responses use the radar years within its own LiDAR year +-1, so cells of one station with different flight years are compared against different radar years
{
  # Model cells of both seasons: station, flight year, and flight month
  lt <- lapply(seasons, function(s) read.csv(paste0(out_dir, s, "_d.csv")) %>% mutate(season = s)) %>%
    bind_rows() %>%
    mutate(year = 2000 + year, month = yymm %% 100)
  
  # Per season and station: model cells, flight years, first and last flight year, share of cells in the most common year, and flight months (year-month combinations)
  lt %>%
    group_by(season, station) %>%
    summarise(cells = n(), years = n_distinct(year), first = min(year), last = max(year),
              top_year_share = round(max(table(year)) / n(), 2), months = n_distinct(yymm), .groups = "drop") %>%
    arrange(season, desc(years)) %>%
    print(n = Inf)
  
  # Spring model cells per flight year in each station's area, colored by flight month
  g <- lt %>%
    filter(season == "spring") %>%
    count(station, year, month) %>%
    ggplot(aes(x = year, y = n, fill = factor(month, levels = 1:12, labels = month.abb))) +
    geom_col() +
    facet_wrap(~ station, scales = "free_y") +
    scale_fill_viridis_d(drop = F) +
    labs(x = "LiDAR flight year", y = "Spring model cells", fill = "Flight month", title = "LiDAR flight dates of the model cells by nearest radar station")
  print(g)
}



##################################################################
# FOREST TYPE GROUPS: 300 m cells with forest LiDAR per FIA forest type group (bar chart)
{
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
  
  # Cells per forest type code -> group (the largest group code at or below the code, e.g. 515 -> 500) -> cells per group
  g <- freq(mask(cov300[["fortype"]], cov300[["p95"]])) %>%
    mutate(group = grp[findInterval(value, as.numeric(names(grp)))]) %>%
    group_by(group) %>%
    summarise(cells = sum(count)) %>%
    ggplot(aes(x = cells, y = reorder(group, cells))) +
    geom_col() +
    labs(x = "300 m cells", y = NULL, title = "Forest type group (TreeMap FORTYPCD) of the 300 m cells with forest LiDAR")
  print(g)
}



##################################################################
# SPRING MODEL CELLS: covariates and responses of the spring model cells as one stack (used by the maps and ridgelines below)
{
  sp <- c(cov300, rast(rf("spring")))
  sp <- mask(sp, noNA(sp[[c(vars, "vir", "hot", "nyr")]], falseNA = T))
  
  # Forest type colors: TreeMap color of every code (treemap_fortypcd.csv) and grey for Other, named by class label
  ft <- read.csv(paste0(out_dir, "treemap_fortypcd.csv"))
  fcol <- c(setNames(ft$col, paste(ft$code, ft$name)), Other = "#bdbdbd")
}



##################################################################
# MAPS of the spring model cells: responses and covariates (deciles as color breaks), forest type, LiDAR flight month, and the voxel layers
{
  # Mapview basemap options, getting rid of watermarks
  options("mapviewBasemaps" = c("Esri.WorldGrayCanvas", "CartoDB.DarkMatter", "OpenStreetMap", "Esri.WorldImagery", "OpenTopoMap"))
  
  # Map layer names (no "/" or "\": mapview uses each name in an HTML dependency name, which cannot contain them)
  lyr <- c(hot = "Hot in all 3 years (top 10% of the nearest station's full range)", p95 = "Canopy height p95 (m)",
           sdp95 = "Canopy height SD (m)", rumple = "Rumple index (canopy surface per flat area)", fhd = "Vertical diversity (FHD)",
           sdfhd = "FHD SD (horizontal heterogeneity)", nrd = "Shrub layer (share of returns at 1-2 m)",
           chili = "Heat load (CHILI)", waterdist = "Distance to water (m)", pd = "Point density (returns per m2)",
           elev = "Elevation (m)")
  
  # Main map: spring VIR shown; other layers hidden (turn on in the layer menu)
  map <- mapview(sp[["vir"]], layer.name = "Spring VIR", col.regions = turbo,
                 at = unique(quantile(values(sp[["vir"]]), seq(0, 1, 0.1), na.rm = T)),
                 na.color = "transparent", maxpixels = 1.4e6)
  
  for (v in names(lyr)) {
    map <- map + mapview(sp[[v]], layer.name = lyr[[v]], col.regions = turbo,
                         at = unique(quantile(values(sp[[v]]), seq(0, 1, 0.1), na.rm = T)),
                         na.color = "transparent", maxpixels = 1.4e6, hide = T)
  }
  
  # Forest type classes as polygons (one per class, so the legend shows names); rows sorted so the colors (row by row) and the legend (alphabetical) agree
  ft_poly <- as.polygons(sp[["ftype"]])
  map <- map + mapview(ft_poly[order(ft_poly$ftype), ], zcol = "ftype", col.regions = unname(fcol[sort(ft_poly$ftype)]),
                       alpha.regions = 1, layer.name = "Forest type class", lwd = 0, hide = T)
  
  # LiDAR flight month as polygons labeled YYYY-MM (alphabetical = chronological, so colors and legend agree)
  fl_poly <- as.polygons(sp[["yymm"]])
  fl_poly$date <- sprintf("20%02d-%02d", fl_poly$yymm %/% 100, fl_poly$yymm %% 100)
  map <- map + mapview(fl_poly[order(fl_poly$date), ], zcol = "date", col.regions = turbo(nrow(fl_poly)),
                       alpha.regions = 1, layer.name = "LiDAR flight month (YYYY-MM)", lwd = 0, hide = T)
  print(map)
  
  # Voxel map: mean % of all returns per height layer (kept pixels); one map layer per height layer
  vlab <- sub("-Inf", "+", paste0("% of returns ", head(brk, -1), "-", brk[-1], " m"))
  vmap <- mapview(sp[[vl[1]]], layer.name = vlab[1], col.regions = viridis,
                  at = unique(quantile(values(sp[[vl[1]]]), seq(0, 1, 0.1), na.rm = T)),
                  na.color = "transparent", maxpixels = 1.4e6)
  for (k in 2:length(vl)) {
    vmap <- vmap + mapview(sp[[vl[k]]], layer.name = vlab[k], col.regions = viridis,
                           at = unique(quantile(values(sp[[vl[k]]]), seq(0, 1, 0.1), na.rm = T)),
                           na.color = "transparent", maxpixels = 1.4e6, hide = T)
  }
  print(vmap)
}



##################################################################
# VOXEL PROFILES: data for the ridgelines below
# Each class's ridge is a step line: over each height layer, its height = the class's mean % of returns in that layer / the layer's width (% of returns per m of height)
# x = height above ground; each layer drawn over its true range (1 m layers to 6 m, then 12 m layers; the open 42+ m layer drawn as 42-54 m); the 0-1 m layer is not kept, so profiles start at 1 m
{
  # Forest type class, canopy height, and voxel layers of the spring model cells
  vox_d <- as.data.frame(sp[[c("ftype", "p95", vl)]], na.rm = T)
  
  # Ridge points of a table with a 'class' column: mean of each voxel layer per class, then two points per layer (its bottom and top, both at % of returns per m)
  ridge_points <- function(d) {
    d %>%
      group_by(class) %>%
      summarise(across(all_of(vl), mean)) %>%
      pivot_longer(all_of(vl), names_to = "layer", values_to = "pct") %>%
      mutate(bot = head(brk, -1)[match(layer, vl)], top = pmin(brk[-1], 54)[match(layer, vl)], value = pct / (top - bot)) %>%
      pivot_longer(c(bot, top), values_to = "x")
  }
}



##################################################################
# VOXEL RIDGELINES BY FOREST TYPE: mean voxel profile of each forest type class (lowest FIA code at the top, Other at the bottom)
{
  pr <- vox_d %>%
    mutate(class = factor(as.character(ftype), levels = rev(sort(unique(as.character(ftype)))))) %>%
    ridge_points()
  
  # Ridgelines (tallest step = 1.5 row spacings); colors = forest type colors
  g <- ggplot(pr, aes(x = x, y = class, height = value, fill = class)) +
    geom_ridgeline(scale = 1.5 / max(pr$value), alpha = 0.8) +
    scale_fill_manual(values = fcol, guide = "none") +
    scale_x_continuous(breaks = c(1, 6, 18, 30, 42, 54)) +
    labs(x = "Height above ground (m)", y = NULL, title = "Mean voxel profile (% of all returns per m of height) by forest type")
  print(g)
}



##################################################################
# VOXEL RIDGELINES BY CANOPY HEIGHT: mean voxel profile of each canopy height class (p95 of the cell's forest pixels, 5 m classes; lowest class at the bottom)
{
  pr <- vox_d %>%
    mutate(class = cut(p95, c(5, 10, 15, 20, 25, 30, Inf), right = F, labels = c("5-10 m", "10-15 m", "15-20 m", "20-25 m", "25-30 m", "30+ m"))) %>%
    ridge_points()
  
  # Ridgelines (tallest step = 1.5 row spacings); viridis colors from short (dark) to tall (yellow)
  g <- ggplot(pr, aes(x = x, y = class, height = value, fill = class)) +
    geom_ridgeline(scale = 1.5 / max(pr$value), alpha = 0.8) +
    scale_fill_manual(values = viridis(nlevels(pr$class)), guide = "none") +
    scale_x_continuous(breaks = c(1, 6, 18, 30, 42, 54)) +
    labs(x = "Height above ground (m)", y = "Canopy height (p95)", title = "Mean voxel profile (% of all returns per m of height) by canopy height")
  print(g)
}