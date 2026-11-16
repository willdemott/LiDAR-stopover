# Stopover intensity vs LiDAR-derived forest structure: HABITAT CLASS MODELS (run after stopover_lidar_covariates.R)
# One VIR GAM per season for each habitat class (ftype: FIA forest type group, or Shrub), fitted to that class's cells only (cells of every radar station)
# Product: two plotly figures (spring, autumn) of predicted VIR vs each covariate, each with a dropdown to choose the habitat class
# Model cells, covariates, and formula as in stopover_lidar_model.R, without the habitat class term (one class per model); the station random effect stays (a class spans several radars)

###################################################################
# Load Packages
{
  library(mgcv)
  library(tidyverse)
  library(plotly)
}



##################################################################
# SETUP: folder, habitat classes, season colors, covariates, and the model formula
{
  # Folder with the stacked raster (stack_300m.tif)
  out_dir <- "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/"
  
  # Habitat classes with their own models (first word of the class label: FIA forest type group code, or Shrub; Other, which pools the rare groups, left out), and season colors (as in stopover_lidar_model.R)
  habitats <- c("500", "800", "100", "700", "400", "600", "160", "960", "120", "900", "Shrub")
  cols <- c(spring = "#1b9e77", autumn = "#d95f02")
  
  # Covariates: column name = panel label
  covs <- c(p95 = "Canopy height p95 (m)", sdp95 = "SD of canopy height (m)", rumple = "Rumple index",
            fhd = "Vertical diversity FHD (point density-adjusted)", sdfhd = "SD of FHD (point density-adjusted)",
            nrd = "Shrub layer: % of returns at 1-2 m", chili = "Heat load, CHILI (cool to warm)",
            waterdist = "Distance to water (m)", edgedist = "Distance to forest edge (m)",
            beamagl = "Radar beam height above ground (m)", elev = "Elevation (m)")
  
  # Formula: one smooth per covariate (k = 20) + canopy height x vertical diversity + station random effect (as in stopover_lidar_model.R, without ftype)
  fml <- as.formula(paste("vir ~", paste0("s(", names(covs), ", k = 20)", collapse = " + "), "+ ti(p95, fhd) + s(station, bs = 're')"))
}



##################################################################
# LOAD DATA: VIR model cells of each season (as in stopover_lidar_model.R), listed habitat classes only
{
  # Point density adjustment (as in stopover_lidar_model.R): the part of y predicted by point density removed and the mean added back
  adj_pd <- function(y, pd) y - fitted(bam(y ~ s(pd))) + mean(y)
  
  # Stacked raster (stopover_lidar_covariates.R): 300 m covariates and responses; habitat class and station as labeled classes
  stk <- terra::rast(paste0(out_dir, "stack_300m.tif"))
  
  # One table per season: cells with every covariate and response (suffix _<season> dropped), FHD and SD of FHD adjusted over all cells (as in the main models),
  # habitat = first word of the class label (e.g. "500", "Shrub"), then VIR > 0 cells of the listed classes
  dat <- lapply(c(spring = "spring", autumn = "autumn"), function(s) {
    terra::as.data.frame(stk[[c(names(covs), "pd", "ftype", "station", paste0(c("vir_", "hot_", "nyr_"), s))]], na.rm = T) %>%
      rename_with(function(n) sub(paste0("_", s, "$"), "", n)) %>%
      mutate(fhd = adj_pd(fhd, pd), sdfhd = adj_pd(sdfhd, pd), habitat = sub(" .*", "", ftype)) %>%
      filter(vir > 0, habitat %in% habitats)
  })
  
  # Check: model cells per habitat class and season
  bind_rows(dat, .id = "season") %>% count(habitat, season) %>% print()
}



##################################################################
# MODELS: one GAM per habitat class and season (Gamma, log link; fREML, discrete), named e.g. "500 spring"
# drop.unused.levels = F keeps every station as a level of s(station), so a class seen by only one radar still fits (stations without cells of the class stay at 0)
{
  runs <- expand.grid(season = names(dat), habitat = habitats, stringsAsFactors = F) %>% mutate(model = paste(habitat, season))
  mods <- lapply(1:nrow(runs), function(i) {
    bam(fml, family = Gamma(link = "log"), data = filter(dat[[runs$season[i]]], habitat == runs$habitat[i]), method = "fREML", discrete = T, drop.unused.levels = F)
  }) %>%
    setNames(runs$model)
}



##################################################################
# CURVES: predicted VIR vs each covariate, one curve per model and covariate
# Each curve: the covariate over its 5th-95th percentile, every other covariate at its mean, the most common station (all from that model's cells); band = 95% CI of the predicted mean VIR
# (the station only moves the whole curve up or down: on the log scale s(station) adds the same amount at every covariate value)
{
  grid <- expand.grid(covariate = names(covs), model = names(mods), stringsAsFactors = F)
  eff <- lapply(1:nrow(grid), function(i) {
    m <- mods[[grid$model[i]]]
    v <- grid$covariate[i]
    # Typical cell (every covariate at its mean, the most common station) repeated 100 times, then v over its 5th-95th percentile
    nd <- m$model %>%
      summarise(across(all_of(names(covs)), mean), station = factor(names(which.max(table(station))), levels(station))) %>%
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
# FIGURES, one per season: predicted VIR vs each covariate (one panel per covariate) for the habitat class chosen in the dropdown
{
  for (s in names(cols)) {
    # One panel per covariate: a curve and its 95% CI band per habitat class (trace name = model, e.g. "500 spring"), in the season's color
    fig <- lapply(names(covs), function(v) {
      filter(eff, covariate == v, season == s) %>%
        plot_ly(x = ~x, name = ~model, color = I(cols[[s]]), showlegend = F) %>%
        add_ribbons(ymin = ~lo, ymax = ~hi, alpha = 0.2, line = list(width = 0)) %>%
        add_lines(y = ~est) %>%
        layout(xaxis = list(title = list(text = covs[[v]])))
    }) %>%
      subplot(nrows = 3, titleX = T, margin = c(0.03, 0.03, 0.03, 0.09))
    
    # Model of each trace (its name), used to show one habitat class's traces at a time
    nm <- sapply(plotly_build(fig)$x$data, function(t) t$name)
    
    # First habitat class shown at the start; dropdown = one button per habitat class, labeled with its full class name (shows that class's traces, hides the rest)
    fig %>%
      style(visible = F, traces = which(nm != paste(habitats[1], s))) %>%
      layout(title = list(text = paste0(s, ": predicted VIR (typical cell of the habitat class, at its most common station; band = 95% CI)")),
             updatemenus = list(list(x = 0, y = 1.08, xanchor = "left", buttons = lapply(habitats, function(h) {
               list(label = grep(paste0("^", h, " "), levels(dat[[s]]$ftype), value = T), method = "restyle", args = list("visible", nm == paste(h, s)))
             })))) %>%
      print()
  }
}