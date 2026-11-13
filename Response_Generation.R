# Stopover intensity vs LiDAR-derived forest structure: FIGURES of the model cells (run order 4 of 5)
# Reads the one stacked raster made by stopover_lidar_covariates.R (stack_300m.tif: the 300 m covariates and the spring and autumn responses)
# vir = mean VIR of the radar years within the cell's LiDAR year +-1; hot = 1 if the cell's VIR reached the top 10% of its station's full range in all 3 of those years
# A cell's station = the radar whose VIR raster covers it; habitat class (ftype) = FIA forest type group, or Shrub when most of the cell's kept LiDAR pixels are scrub (p95 1-5 m)

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
# SETUP: folder, seasons, model covariates, and the stacked raster
{
  # Folder with the stacked raster
  out_dir <- "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/"
  
  # Seasons, and height layer breaks (m) of the voxel layers (as in the covariate script)
  seasons <- c("spring", "autumn")
  brk <- c(1, 2, 3, 4, 5, 6, 18, 30, 42, Inf)
  
  # Covariates of the model cells (year and yymm = LiDAR flight year and month, YY and YYMM; not model covariates)
  vars <- c("p95", "sdp95", "rumple", "fhd", "sdfhd", "nrd", "pd", "chili", "waterdist", "edgedist", "raddist", "elev", "beamagl", "ftype", "station", "year", "yymm")
  
  # Stacked raster (300 m covariates; responses vir_<season>, hot_<season>, nyr_<season>; habitat class and station as labeled classes), and voxel layer names (v1_2 ... v42_Inf)
  stk <- rast(paste0(out_dir, "stack_300m.tif"))
  vl <- grep("^v[0-9]", names(stk), value = T)
}



##################################################################
# SAMPLE POINTS PER RADAR: model rows (300 m cells) per station and season; each cell counts for its own station (the radar its VIR comes from)
{
  g <- lapply(seasons, function(s) {
    as.data.frame(stk[[c(vars, paste0(c("vir_", "hot_", "nyr_"), s))]], na.rm = T) %>% mutate(season = s)   # model cells of season s (every covariate and response present)
  }) %>%
    bind_rows() %>%
    count(season, station) %>%
    ggplot(aes(x = n, y = reorder(station, n), fill = season)) +
    geom_col(position = position_dodge()) +
    scale_fill_manual(values = c(spring = "#1b9e77", autumn = "#d95f02")) +
    labs(x = "Model cells (300 m)", y = NULL, fill = "Season", title = "Sample points per radar station")
  print(g)
}



##################################################################
# HABITAT CLASSES: 300 m cells with forest LiDAR per habitat class (FIA forest type group, Other, or Shrub; bar chart)
{
  g <- freq(mask(stk[["ftype"]], stk[["p95"]])) %>%
    ggplot(aes(x = count, y = reorder(value, count))) +
    geom_col() +
    labs(x = "300 m cells", y = NULL, title = "Habitat class of the 300 m cells with forest LiDAR")
  print(g)
}



##################################################################
# SEASONAL MODEL CELLS: covariates and responses of the spring (sp) and autumn (fa) model cells, one stack per season (used by the maps and ridgelines below)
# Autumn has the larger radar coverage, so the autumn stack (fa) is the source of the covariate layers in the maps and ridgelines
{
  # Covariates, voxel layers, and the season's responses (renamed vir, hot, nyr), masked to the season's model cells (every covariate and response present)
  sp <- c(stk[[c(vars, vl)]], setNames(stk[[c("vir_spring", "hot_spring", "nyr_spring")]], c("vir", "hot", "nyr")))
  fa <- c(stk[[c(vars, vl)]], setNames(stk[[c("vir_autumn", "hot_autumn", "nyr_autumn")]], c("vir", "hot", "nyr")))
  sp <- mask(sp, noNA(sp[[c(vars, "vir", "hot", "nyr")]], falseNA = T))
  fa <- mask(fa, noNA(fa[[c(vars, "vir", "hot", "nyr")]], falseNA = T))
  
  # Habitat class colors: one color per class, in class order (as in the forest type curves of stopover_lidar_model.R), grey for Other; named by class label
  lv <- levels(stk[["ftype"]])[[1]]$ftype
  fcol <- setNames(hcl.colors(length(lv), "Dark 3"), lv)
  fcol["Other"] <- "#bdbdbd"
  
  # Habitat classes of the autumn model cells as polygons, one per class (used by both maps below)
  ft_poly <- as.polygons(fa[["ftype"]])
}



##################################################################
# MAPS of the model cells: responses (spring and autumn) and covariates (deciles as color breaks), habitat class, LiDAR flight month, and the voxel layers
# Covariate, habitat class, flight month, and voxel layers come from the autumn model cells (fa), which have the larger radar coverage
{
  # Mapview basemap options, getting rid of watermarks
  options("mapviewBasemaps" = c("Esri.WorldGrayCanvas", "CartoDB.DarkMatter", "OpenStreetMap", "Esri.WorldImagery", "OpenTopoMap"))
  
  # Map layer names (no "/" or "\": mapview uses each name in an HTML dependency name, which cannot contain them)
  lyr <- c(p95 = "Canopy height p95 (m)", sdp95 = "Canopy height SD (m)", rumple = "Rumple index (canopy surface per flat area)",
           fhd = "Vertical diversity (FHD)", sdfhd = "FHD SD (horizontal heterogeneity)",
           nrd = "Shrub layer (share of returns at 1-2 m)", chili = "Heat load (CHILI)", waterdist = "Distance to water (m)",
           edgedist = "Distance to forest edge (m)", pd = "Point density (returns per m2)", elev = "Elevation (m)")
  
  # Main map: spring and autumn VIR and hotspot layers (shown); the covariate, habitat class, and flight month layers added below are hidden (turn on in the layer menu)
  map <- mapview(sp[["vir"]], layer.name = "Spring VIR", col.regions = turbo,
                 at = unique(quantile(values(sp[["vir"]]), seq(0, 1, 0.1), na.rm = T)),
                 na.color = "transparent", maxpixels = 1.4e6) +
    mapview(fa[["vir"]], layer.name = "Autumn VIR", col.regions = turbo,
            at = unique(quantile(values(fa[["vir"]]), seq(0, 1, 0.1), na.rm = T)),
            na.color = "transparent", maxpixels = 1.4e6) +
    mapview(sp[["hot"]], layer.name = "Spring Hotspot in all 3 years", col.regions = c("gray55", "firebrick"),
            at = c(-0.5, 0.5, 1.5), na.color = "transparent", maxpixels = 1.4e6) +
    mapview(fa[["hot"]], layer.name = "Autumn Hotspot in all 3 years", col.regions = c("gray55", "firebrick"),
            at = c(-0.5, 0.5, 1.5), na.color = "transparent", maxpixels = 1.4e6)
  
  # Covariate layers of the autumn model cells (hidden)
  for (v in names(lyr)) {
    map <- map + mapview(fa[[v]], layer.name = lyr[[v]], col.regions = turbo,
                         at = unique(quantile(values(fa[[v]]), seq(0, 1, 0.1), na.rm = T)),
                         na.color = "transparent", maxpixels = 1.4e6, hide = T)
  }
  
  # Habitat classes of the autumn model cells (polygons, one per class, so the legend shows names); rows sorted so the colors (row by row) and the legend (alphabetical) agree
  map <- map + mapview(ft_poly[order(ft_poly$ftype), ], zcol = "ftype", col.regions = unname(fcol[sort(ft_poly$ftype)]),
                       alpha.regions = 1, layer.name = "Habitat class", lwd = 0, hide = T)
  
  # LiDAR flight month of the autumn model cells as polygons labeled YYYY-MM (alphabetical = chronological, so colors and legend agree)
  fl_poly <- as.polygons(fa[["yymm"]])
  fl_poly$date <- sprintf("20%02d-%02d", fl_poly$yymm %/% 100, fl_poly$yymm %% 100)
  map <- map + mapview(fl_poly[order(fl_poly$date), ], zcol = "date", col.regions = turbo(nrow(fl_poly)),
                       alpha.regions = 1, layer.name = "LiDAR flight month (YYYY-MM)", lwd = 0, hide = T)
  print(map)
  
  # Voxel map of the autumn model cells: mean % of all returns per height layer (kept pixels); one map layer per height layer
  vlab <- sub("-Inf", "+", paste0("% of returns ", head(brk, -1), "-", brk[-1], " m"))
  vmap <- mapview(fa[[vl[1]]], layer.name = vlab[1], col.regions = viridis,
                  at = unique(quantile(values(fa[[vl[1]]]), seq(0, 1, 0.1), na.rm = T)),
                  na.color = "transparent", maxpixels = 1.4e6)
  for (k in 2:length(vl)) {
    vmap <- vmap + mapview(fa[[vl[k]]], layer.name = vlab[k], col.regions = viridis,
                           at = unique(quantile(values(fa[[vl[k]]]), seq(0, 1, 0.1), na.rm = T)),
                           na.color = "transparent", maxpixels = 1.4e6, hide = T)
  }
  print(vmap)
}



##################################################################
# HABITAT MAP of the model cells: spring and autumn VIR (shown), plus one layer per habitat class of the autumn model cells (hidden; turn each on in the layer menu)
{
  hmap <- mapview(sp[["vir"]], layer.name = "Spring VIR", col.regions = turbo,
                  at = unique(quantile(values(sp[["vir"]]), seq(0, 1, 0.1), na.rm = T)),
                  na.color = "transparent", maxpixels = 1.4e6) +
    mapview(fa[["vir"]], layer.name = "Autumn VIR", col.regions = turbo,
            at = unique(quantile(values(fa[["vir"]]), seq(0, 1, 0.1), na.rm = T)),
            na.color = "transparent", maxpixels = 1.4e6)
  
  # One layer per habitat class, in its class color (layer names without "/", which mapview cannot use in a layer name)
  for (k in sort(ft_poly$ftype)) {
    hmap <- hmap + mapview(ft_poly[ft_poly$ftype == k, ], col.regions = fcol[[k]], alpha.regions = 1, lwd = 0,
                           layer.name = gsub("/", "-", k), hide = T)
  }
  print(hmap)
}



##################################################################
# VOXEL PROFILES: data for the ridgelines below
# Each class's ridge is a step line: over each height layer, its height = the class's mean % of returns in that layer / the layer's width (% of returns per m of height)
# x = height above ground; each layer drawn over its true range (1 m layers to 6 m, then 12 m layers; the open 42+ m layer drawn as 42-54 m); the 0-1 m layer is not kept, so profiles start at 1 m
{
  # Habitat class, canopy height, and voxel layers of the autumn model cells (the larger radar coverage)
  vox_d <- as.data.frame(fa[[c("ftype", "p95", vl)]], na.rm = T)
  
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
# VOXEL RIDGELINES BY HABITAT CLASS: mean voxel profile of each habitat class (lowest FIA group code at the top; Other and Shrub at the bottom)
{
  pr <- vox_d %>%
    mutate(class = factor(as.character(ftype), levels = rev(sort(unique(as.character(ftype)))))) %>%
    ridge_points()
  
  # Ridgelines (tallest step = 1.5 row spacings); colors = habitat class colors
  g <- ggplot(pr, aes(x = x, y = class, height = value, fill = class)) +
    geom_ridgeline(scale = 1.5 / max(pr$value), alpha = 0.8) +
    scale_fill_manual(values = fcol, guide = "none") +
    scale_x_continuous(breaks = c(1, 6, 18, 30, 42, 54)) +
    labs(x = "Height above ground (m)", y = NULL, title = "Mean voxel profile (% of all returns per m of height) by habitat class")
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