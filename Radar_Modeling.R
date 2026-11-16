# Stopover intensity vs LiDAR-derived forest structure: STATION MODELS (run after stopover_lidar_covariates.R)
# One VIR GAM per season for each of the most-sampled radar stations, fitted to that station's cells only (the cells whose VIR comes from that radar)
# Product: two plotly figures (spring, autumn) of predicted VIR vs each covariate, each with a dropdown to choose the station
# Model cells, covariates, and formula as in stopover_lidar_model.R, without the station random effect (one station per model)

###################################################################
# Load Packages
{
  library(mgcv)
  library(tidyverse)
  library(plotly)
}



##################################################################
# SETUP: folder, stations, season colors, covariates, and the model formula
{
  # Folder with the stacked raster (stack_300m.tif)
  out_dir <- "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/"
  
  # Stations with their own models (the most model cells), and season colors (as in stopover_lidar_model.R)
  stations <- c("KCCX", "KBGM", "KRLX", "KFCX", "KENX", "KBOX")
  cols <- c(spring = "#1b9e77", autumn = "#d95f02")
  
  # Covariates: column name = panel label
  covs <- c(p95 = "Canopy height p95 (m)", sdp95 = "SD of canopy height (m)", rumple = "Rumple index",
            fhd = "Vertical diversity FHD (point density-adjusted)", sdfhd = "SD of FHD (point density-adjusted)",
            nrd = "Shrub layer: % of returns at 1-2 m", chili = "Heat load, CHILI (cool to warm)",
            waterdist = "Distance to water (m)", edgedist = "Distance to forest edge (m)",
            beamagl = "Radar beam height above ground (m)", elev = "Elevation (m)")
  
  # Formula: habitat class + one smooth per covariate (k = 20) + canopy height x vertical diversity (as in stopover_lidar_model.R, without s(station))
  fml <- as.formula(paste("vir ~ ftype +", paste0("s(", names(covs), ", k = 20)", collapse = " + "), "+ ti(p95, fhd)"))
}



##################################################################
# LOAD DATA: VIR model cells of each season (as in stopover_lidar_model.R), listed stations only
{
  # Point density adjustment (as in stopover_lidar_model.R): the part of y predicted by point density removed and the mean added back
  adj_pd <- function(y, pd) y - fitted(bam(y ~ s(pd))) + mean(y)
  
  # Stacked raster (stopover_lidar_covariates.R): 300 m covariates and responses; habitat class and station as labeled classes
  stk <- terra::rast(paste0(out_dir, "stack_300m.tif"))
  
  # One table per season: cells with every covariate and response (suffix _<season> dropped), FHD and SD of FHD adjusted over all cells (as in the main models), then VIR > 0 cells of the listed stations
  dat <- lapply(c(spring = "spring", autumn = "autumn"), function(s) {
    terra::as.data.frame(stk[[c(names(covs), "pd", "ftype", "station", paste0(c("vir_", "hot_", "nyr_"), s))]], na.rm = T) %>%
      rename_with(function(n) sub(paste0("_", s, "$"), "", n)) %>%
      mutate(fhd = adj_pd(fhd, pd), sdfhd = adj_pd(sdfhd, pd)) %>%
      filter(vir > 0, station %in% stations)
  })
  
  # Check: model cells per station and season
  bind_rows(dat, .id = "season") %>% count(station, season) %>% print()
}



##################################################################
# MODELS: one GAM per station and season (Gamma, log link; fREML, discrete), named e.g. "KBGM spring"
{
  runs <- expand.grid(season = names(dat), station = stations, stringsAsFactors = F) %>% mutate(model = paste(station, season))
  mods <- lapply(1:nrow(runs), function(i) {
    bam(fml, family = Gamma(link = "log"), data = filter(dat[[runs$season[i]]], station == runs$station[i]), method = "fREML", discrete = T)
  }) %>%
    setNames(runs$model)
}



##################################################################
# CURVES: predicted VIR vs each covariate, one curve per model and covariate
# Each curve: the covariate over its 5th-95th percentile, every other covariate at its mean, the most common habitat class (all from that model's cells); band = 95% CI of the predicted mean VIR
{
  grid <- expand.grid(covariate = names(covs), model = names(mods), stringsAsFactors = F)
  eff <- lapply(1:nrow(grid), function(i) {
    m <- mods[[grid$model[i]]]
    v <- grid$covariate[i]
    # Typical cell (every covariate at its mean, the most common habitat class) repeated 100 times, then v over its 5th-95th percentile
    nd <- m$model %>%
      summarise(across(all_of(names(covs)), mean), ftype = factor(names(which.max(table(ftype))), levels(ftype))) %>%
      slice(rep(1, 100))
    nd[[v]] <- seq(quantile(m$model[[v]], 0.05), quantile(m$model[[v]], 0.95), length.out = 100)
    # Prediction and its SE on the log scale, back-transformed to VIR (exp)
    p <- predict(m, nd, se.fit = T)
    data.frame(model = grid$model[i], covariate = v, x = nd[[v]], est = exp(p$fit), lo = exp(p$fit - 1.96 * p$se.fit), hi = exp(p$fit + 1.96 * p$se.fit))
  }) %>%
    bind_rows() %>%
    mutate(season = word(model, 2))
}



##################################################################
# FIGURES, one per season: predicted VIR vs each covariate (one panel per covariate) for the station chosen in the dropdown
{
  for (s in names(cols)) {
    # One panel per covariate: a curve and its 95% CI band per station (trace name = model, e.g. "KBGM spring"), in the season's color
    fig <- lapply(names(covs), function(v) {
      filter(eff, covariate == v, season == s) %>%
        plot_ly(x = ~x, name = ~model, color = I(cols[[s]]), showlegend = F) %>%
        add_ribbons(ymin = ~lo, ymax = ~hi, alpha = 0.2, line = list(width = 0)) %>%
        add_lines(y = ~est) %>%
        layout(xaxis = list(title = list(text = covs[[v]])))
    }) %>%
      subplot(nrows = 3, titleX = T, margin = c(0.03, 0.03, 0.03, 0.09))
    
    # Model of each trace (its name), used to show one station's traces at a time
    nm <- sapply(plotly_build(fig)$x$data, function(t) t$name)
    
    # First station shown at the start; dropdown = one button per station (shows that station's traces, hides the rest)
    fig %>%
      style(visible = F, traces = which(!startsWith(nm, stations[1]))) %>%
      layout(title = list(text = paste0(s, ": predicted VIR (typical cell of the station; band = 95% CI)")),
             updatemenus = list(list(x = 0, y = 1.08, xanchor = "left", buttons = lapply(stations, function(st) {
               list(label = st, method = "restyle", args = list("visible", startsWith(nm, st)))
             })))) %>%
      print()
  }
}