# Stopover intensity vs LiDAR-derived forest structure: MODELING
# Fits GAMs to spring_d.csv and autumn_d.csv (from stopover_lidar_responses.R), then makes the effect, interaction, and forest type figures and the model stats
# vir = mean VIR (Gamma, log link); hot = 1 if hot in all 3 years of the LiDAR year +-1 (binomial); both with forest type, one smooth per covariate, ti(p95, fhd), and a station random effect
# Covariates are listed once (covs, in SETUP): the formulas, figures, and stats follow that list
# Pipes (%>%) read left to right: take a table, do one step, pass the result to the next step

###################################################################
# Load Packages
{
  library(mgcv)
  library(tidyverse)
}



##################################################################
# SETUP: folder, covariates, formulas, model runs, and figure settings
{
  # Folder with the model dataframes; the model stats are written here too
  out_dir <- "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/"
  
  # Covariates: column name, axis label, and whether it is a forest-structure covariate (structure covariates get per-forest-type curves)
  covs <- tribble(
    ~name,       ~label,                                              ~structure,
    "p95",       "Canopy height\np95 (m)",                            TRUE,
    "sdp95",     "SD of canopy\nheight (m)",                          TRUE,
    "rumple",    "Rumple index (canopy\nsurface / flat area)",        TRUE,
    "fhd",       "Vertical diversity\nFHD (point density-adjusted)",  TRUE,
    "sdfhd",     "SD of FHD\n(point density-adjusted)",               TRUE,
    "nrd",       "Shrub layer: share of\nreturns at 1-2 m (%)",       TRUE,
    "chili",     "Heat load, CHILI\n(cool to warm)",                  FALSE,
    "waterdist", "Distance to\nwater (m)",                            FALSE,
    "beamagl",   "Radar beam height\nabove ground (m)",               FALSE,
    "elev",      "Elevation (m)",                                     FALSE)
  
  # Formulas: forest type (fixed effect) + one smooth per covariate (k = 20) + canopy height x vertical diversity interaction + station random effect
  rhs <- paste("ftype +", paste0("s(", covs$name, ", k = 20)", collapse = " + "), "+ ti(p95, fhd) + s(station, bs = 're')")
  fml <- list(vir = as.formula(paste("vir ~", rhs)), hot = as.formula(paste("hot ~", rhs)))
  fam <- list(vir = Gamma(link = "log"), hot = binomial())
  print(fml)
  
  # Model runs: one row per season and response; model name = "season response"
  runs <- expand.grid(response = c("vir", "hot"), season = c("spring", "autumn"), stringsAsFactors = F) %>%
    mutate(model = paste(season, response))
  
  # Largest VIR used in the VIR models
  vir_cap <- 7000
  
  # Forest types (FIA codes): base type of the figures, and the types with their own curves in the forest type curves
  base_code <- "801"
  type_codes <- c("801", "503", "515")
  
  # Radar stations with their own curves in the station curves (ICAO codes, e.g. c("KENX", "KBGM")); NULL = every station; stations not listed share one group, Other
  sta_pick <- NULL
  
  # Response labels, and season colors (colorblind-friendly)
  ylab <- c(vir = "Predicted VIR", hot = "P(hot all 3 years)")
  cols <- c(spring = "#1b9e77", autumn = "#d95f02")
}



##################################################################
# LOAD DATA: model dataframes (ftype and station read as factors), FHD and SD of FHD adjusted for point density, and the forest type groups of the per-type curves
{
  # Point density adjustment: the part of y predicted by point density (a smooth fit across the model cells) removed and the mean added back,
  # so a value says how a cell compares with cells of the same point density (pd itself is not a model covariate)
  adj_pd <- function(y, pd) y - fitted(bam(y ~ s(pd))) + mean(y)
  
  # Model dataframes, one per season; ft3 = listed type or Other, ordered with Other first (s(v, by = ft3) then adds one difference curve per listed type to the shared curve)
  dat <- lapply(c(spring = "spring", autumn = "autumn"), function(s) {
    read.csv(paste0(out_dir, s, "_d.csv"), stringsAsFactors = T) %>%
      mutate(fhd = adj_pd(fhd, pd),
             sdfhd = adj_pd(sdfhd, pd),
             ft3 = ordered(if_else(substr(ftype, 1, 3) %in% type_codes, substr(ftype, 1, 3), "Other"), levels = c("Other", type_codes)),
             st2 = factor(if (is.null(sta_pick)) as.character(station) else if_else(station %in% sta_pick, as.character(station), "Other")))
  })
  
  # Full label of the base forest type, e.g. "801 Sugar maple/beech/yellow birch"
  base_type <- grep(paste0("^", base_code, " "), levels(dat$spring$ftype), value = T)
  
  # Model data of season s and response k: VIR models use cells with 0 < VIR <= vir_cap (Gamma needs positive values); hot models use cells with all 3 years tested
  model_data <- function(s, k) if (k == "vir") filter(dat[[s]], vir > 0, vir <= vir_cap) else filter(dat[[s]], nyr == 3)
  
  # Check: cells per season and model, hot cells, and correlations of the adjusted FHD and SD of FHD with point density (near 0)
  bind_rows(dat, .id = "season") %>%
    group_by(season) %>%
    summarise(cells = n(), vir_model = sum(vir > 0 & vir <= vir_cap), hot_model = sum(nyr == 3), hot_1 = sum(hot[nyr == 3]),
              cor_fhd_pd = cor(fhd, pd), cor_sdfhd_pd = cor(sdfhd, pd)) %>%
    print()
  plot(dat$spring$pd, dat$spring$fhd, pch = ".", xlab = "Point density (returns per m2)", ylab = "FHD (point density-adjusted)")
}



##################################################################
# POINT DENSITY CHECK: percentiles of the adjusted SD of FHD per point density class (should no longer shift with point density), and its map (spring)
{
  dat$spring %>%
    mutate(pd_class = cut(pd, c(0, 5, 10, 20, 30, 40, Inf))) %>%
    group_by(pd_class) %>%
    summarise(cells = n(), q05 = quantile(sdfhd, 0.05), q50 = quantile(sdfhd, 0.5), q95 = quantile(sdfhd, 0.95)) %>%
    print()
  
  # Map with decile colors, to compare with the raw map across LiDAR collection boundaries
  sd_adj <- terra::rast(dat$spring[, c("x", "y", "sdfhd")], type = "xyz", crs = "EPSG:6350")
  print(mapview::mapview(sd_adj, col.regions = viridisLite::turbo, at = unique(quantile(dat$spring$sdfhd, seq(0, 1, 0.1))),
                         layer.name = "SD of FHD (point density-adjusted)", na.color = "transparent"))
}



##################################################################
# MODELS: one GAM per run (fREML, discrete), named "spring vir", "spring hot", "autumn vir", "autumn hot"
{
  mods <- lapply(1:nrow(runs), function(i) {
    bam(fml[[runs$response[i]]], family = fam[[runs$response[i]]], data = model_data(runs$season[i], runs$response[i]), method = "fREML", discrete = T)
  }) %>%
    setNames(runs$model)
  
  # Summaries (also used by MODEL STATS) and concurvity (GAM version of collinearity: near 1 = a smooth is largely explained by the other terms)
  sums <- lapply(mods, summary)
  for (m in names(mods)) {
    cat("\n=====", m, "=====\n")
    print(sums[[m]])
    print(round(concurvity(mods[[m]], full = T)["worst", ], 2))
  }
}



##################################################################
# PREDICTION HELPERS (used by the figures)
{
  # n rows of a typical cell of model m: every covariate at its mean, forest type 'type' (the base type unless given)
  typical <- function(m, n, type = base_type) {
    nd <- m$model[rep(1, n), ]
    for (u in covs$name) nd[[u]] <- mean(m$model[[u]])
    nd$ftype[] <- type
    nd
  }
  
  # Prediction for the rows of nd, with a 95% CI of the curve's shape: link scale without s(station) (= an average station),
  # SE from the smooths that use the covariates in 'terms' only (as plot.gam with seWithMean = FALSE), back-transformed (exp for VIR, inverse logit for hot)
  shape_ci <- function(m, nd, terms) {
    p <- predict(m, nd, exclude = "s(station)")
    X <- predict(m, nd, type = "lpmatrix")
    X[, -unlist(lapply(m$smooth, function(sm) if (any(sm$term %in% terms)) sm$first.para:sm$last.para))] <- 0
    se <- sqrt(rowSums((X %*% vcov(m)) * X))
    data.frame(est = m$family$linkinv(p), lo = m$family$linkinv(p - 1.96 * se), hi = m$family$linkinv(p + 1.96 * se))
  }
  
  # Panel labels: response and covariate labels as factors (panel order = covs order)
  lab <- function(d) d %>% mutate(response = factor(ylab[response], levels = ylab), covariate = factor(covariate, levels = covs$name, labels = covs$label))
}



##################################################################
# EFFECT FIGURES, one per season: predicted VIR (rows 1-2) and P(hot all 3 years) (rows 3-4) vs each covariate
# Each line: the covariate over its 5th-95th percentile, every other covariate at its mean, the base forest type, an average station; band = 95% CI of the curve's shape; rug = observed values
{
  # One curve per model and covariate
  eff_grid <- expand.grid(model = names(mods), covariate = covs$name, stringsAsFactors = F) %>% left_join(runs, by = "model")
  eff <- lapply(1:nrow(eff_grid), function(i) {
    m <- mods[[eff_grid$model[i]]]
    v <- eff_grid$covariate[i]
    nd <- typical(m, 100)
    nd[[v]] <- seq(quantile(m$model[[v]], 0.05), quantile(m$model[[v]], 0.95), length.out = 100)
    shape_ci(m, nd, v) %>% mutate(season = eff_grid$season[i], response = eff_grid$response[i], covariate = v, x = nd[[v]])
  }) %>%
    bind_rows()
  
  # Rug: observed values of each covariate in each model, every 0.1% from the 5th to 95th percentile
  rug <- lapply(1:nrow(eff_grid), function(i) {
    data.frame(season = eff_grid$season[i], response = eff_grid$response[i], covariate = eff_grid$covariate[i],
               x = quantile(mods[[eff_grid$model[i]]]$model[[eff_grid$covariate[i]]], seq(0.05, 0.95, 0.001)))
  }) %>%
    bind_rows()
  
  for (s in c("spring", "autumn")) {
    g <- ggplot(lab(filter(eff, season == s)), aes(x = x, y = est)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = cols[[s]], alpha = 0.2) +
      geom_line(colour = cols[[s]], linewidth = 0.8) +
      geom_rug(data = lab(filter(rug, season == s)), aes(x = x), inherit.aes = F, sides = "b", colour = cols[[s]], alpha = 0.1) +
      facet_wrap(vars(response, covariate), ncol = 5, scales = "free") +
      scale_x_continuous(n.breaks = 4) +
      labs(x = NULL, y = NULL, title = paste0(s, ": predicted VIR (rows 1-2) and P(hot all 3 years) (rows 3-4)"))
    print(g)
  }
}



##################################################################
# INTERACTION FIGURES, one per season: predicted VIR (top) and P(hot all 3 years) (bottom) vs FHD at low, median, and high canopy height
# Canopy heights = 10th, 50th, and 90th percentile of p95 (both seasons pooled); FHD over its 5th-95th percentile; everything else as in the effect figures; band = shape CI of s(p95), s(fhd), and ti(p95, fhd)
{
  pq <- quantile(c(dat$spring$p95, dat$autumn$p95), c(0.1, 0.5, 0.9), names = F)
  
  # 100 FHD values at each of the 3 canopy heights, per model
  int <- lapply(1:nrow(runs), function(i) {
    m <- mods[[runs$model[i]]]
    nd <- typical(m, 300)
    nd$p95 <- rep(pq, each = 100)
    nd$fhd <- rep(seq(quantile(m$model$fhd, 0.05), quantile(m$model$fhd, 0.95), length.out = 100), 3)
    shape_ci(m, nd, c("p95", "fhd")) %>%
      mutate(season = runs$season[i], response = runs$response[i], x = nd$fhd,
             h = factor(paste0("p95 = ", round(nd$p95, 1), " m"), levels = paste0("p95 = ", round(pq, 1), " m")))
  }) %>%
    bind_rows()
  
  for (s in c("spring", "autumn")) {
    g <- int %>%
      filter(season == s) %>%
      mutate(response = factor(ylab[response], levels = ylab)) %>%
      ggplot(aes(x = x, y = est)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = cols[[s]], alpha = 0.2) +
      geom_line(colour = cols[[s]], linewidth = 0.8) +
      facet_grid(response ~ h, scales = "free_y") +
      labs(x = covs$label[covs$name == "fhd"], y = NULL, title = paste0(s, ": canopy height x vertical diversity"))
    print(g)
  }
}



##################################################################
# FOREST TYPE FIGURE: each forest type class vs the base type, by season
# Two otherwise identical cells (covariates at their means, same station) that differ only in forest type: ratio of predicted VIR (left) and odds ratio of hot in all 3 years (right); bars = 95% CIs; 1 = no difference
{
  ftc <- lapply(1:nrow(runs), function(i) {
    m <- mods[[runs$model[i]]]
    lv <- setdiff(levels(m$model$ftype), base_type)              # every other class
    nd <- typical(m, length(lv) + 1)
    nd$ftype[] <- c(base_type, lv)                                 # row 1 = base type, then one row per other class
    X <- predict(m, nd, type = "lpmatrix")
    X <- sweep(X[-1, , drop = F], 2, X[1, ])                       # each class minus the base type (link scale)
    est <- drop(X %*% coef(m))
    se <- sqrt(rowSums((X %*% vcov(m)) * X))
    data.frame(season = runs$season[i], response = runs$response[i], term = lv, est = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
  }) %>%
    bind_rows()
  
  # Class order: mean VIR ratio over both seasons (highest at the top)
  ord <- ftc %>% filter(response == "vir") %>% group_by(term) %>% summarise(r = mean(est)) %>% arrange(r) %>% pull(term)
  
  g <- ftc %>%
    mutate(term = factor(term, levels = ord),
           response = factor(response, levels = c("vir", "hot"), labels = paste0(c("VIR ratio", "Odds ratio, hot all 3 years"), " (vs ", base_code, ")"))) %>%
    ggplot(aes(x = est, y = term, colour = season)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5)) +
    facet_wrap(~ response, scales = "free_x") +
    scale_colour_manual(values = cols) +
    labs(x = NULL, y = NULL, colour = "Season")
  print(g)
}



##################################################################
# FOREST TYPE CURVES, one figure per model: predicted VIR or P(hot all 3 years) vs each structure covariate for the listed forest types (type_codes)
# Forest-type models, used only here: each main model plus s(v, by = ft3, k = 10) for every structure covariate; each type's curve = shared curve + its difference curve
# Lines cover the 5th-95th percentile of the covariate; band = 95% CI of the curve's shape (level differences between types: see the forest type figure)
{
  # Difference curves added to the main formulas
  by_terms <- as.formula(paste(". ~ . +", paste0("s(", covs$name[covs$structure], ", by = ft3, k = 10)", collapse = " + ")))
  
  # Forest-type models (same data, family, and settings as the main models; slow, like the main fits)
  ftm <- lapply(1:nrow(runs), function(i) {
    bam(update(fml[[runs$response[i]]], by_terms), family = fam[[runs$response[i]]], data = model_data(runs$season[i], runs$response[i]), method = "fREML", discrete = T)
  }) %>%
    setNames(runs$model)
  
  # Labels of the listed types, and their TreeMap colors (as in the maps)
  ft_lab <- sapply(type_codes, function(t) grep(paste0("^", t, " "), levels(dat$spring$ftype), value = T))
  ft_col <- read.csv(paste0(out_dir, "treemap_fortypcd.csv")) %>%
    filter(code %in% type_codes) %>%
    mutate(label = paste(code, name)) %>%
    pull(col, name = label)
  
  # One curve per model, structure covariate, and listed type
  cur_grid <- expand.grid(model = names(ftm), covariate = covs$name[covs$structure], type = type_codes, stringsAsFactors = F) %>% left_join(runs, by = "model")
  ftcur <- lapply(1:nrow(cur_grid), function(i) {
    m <- ftm[[cur_grid$model[i]]]
    v <- cur_grid$covariate[i]
    nd <- typical(m, 100, type = ft_lab[[cur_grid$type[i]]])
    nd$ft3[] <- cur_grid$type[i]
    nd[[v]] <- seq(quantile(m$model[[v]], 0.05), quantile(m$model[[v]], 0.95), length.out = 100)
    shape_ci(m, nd, v) %>% mutate(model = cur_grid$model[i], response = cur_grid$response[i], covariate = v, type = ft_lab[[cur_grid$type[i]]], x = nd[[v]])
  }) %>%
    bind_rows()
  
  for (mn in names(ftm)) {
    g <- ftcur %>%
      filter(model == mn) %>%
      mutate(covariate = factor(covariate, levels = covs$name, labels = gsub("\n", " ", covs$label))) %>%
      ggplot(aes(x = x, y = est, colour = type)) +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = type), alpha = 0.2, colour = NA, show.legend = F) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~ covariate, scales = "free") +
      scale_colour_manual(values = ft_col) +
      scale_fill_manual(values = ft_col, guide = "none") +
      labs(x = NULL, y = ylab[[runs$response[runs$model == mn]]], colour = NULL, title = paste0(mn, ": by forest type"))
    print(g)
  }
}



##################################################################
# STATION CURVES, one figure per model: predicted VIR or P(hot all 3 years) vs each structure covariate for each radar station (sta_pick)
# Station models, used only here: each main model plus s(v, by = st2, k = 10, m = 1) for every structure covariate = model "GI" of Pedersen et al. (2019):
#   one shared curve per covariate plus one deviation curve per station; m = 1 penalizes each deviation toward a flat line (no difference from the shared curve)
# Each station's line = shared curve + its deviation, drawn over that station's own 5th-95th percentile of the covariate (stations with < 100 cells skipped)
# s(station) is left out of the lines, so they compare curve shapes, not station levels; no CI bands (too many lines)
{
  # Deviation curves added to the main formulas
  st_terms <- as.formula(paste(". ~ . +", paste0("s(", covs$name[covs$structure], ", by = st2, k = 10, m = 1)", collapse = " + ")))
  
  # Station models (same data, family, and settings as the main models; slow, like the main fits)
  stm <- lapply(1:nrow(runs), function(i) {
    bam(update(fml[[runs$response[i]]], st_terms), family = fam[[runs$response[i]]], data = model_data(runs$season[i], runs$response[i]), method = "fREML", discrete = T)
  }) %>%
    setNames(runs$model)
  
  # Deviation terms of every station model: edf (1 = straight-line deviation, more = curved) and p (deviation different from 0), saved to out_dir
  st_trm <- lapply(names(stm), function(mn) {
    st <- summary(stm[[mn]])$s.table
    st <- st[grepl("st2", rownames(st)), , drop = F]
    data.frame(model = mn, term = rownames(st), edf = round(st[, 1], 2), p = signif(st[, 4], 3), row.names = NULL)
  }) %>%
    bind_rows()
  print(st_trm)
  write.csv(st_trm, paste0(out_dir, "station_term_stats.csv"), row.names = F)
  
  # One line per model, structure covariate, and station
  st_grid <- expand.grid(model = names(stm), covariate = covs$name[covs$structure], st = unique(c(levels(dat$spring$st2), levels(dat$autumn$st2))), stringsAsFactors = F) %>%
    left_join(runs, by = "model")
  stcur <- lapply(1:nrow(st_grid), function(i) {
    m <- stm[[st_grid$model[i]]]
    v <- st_grid$covariate[i]
    xs <- m$model[[v]][m$model$st2 == st_grid$st[i]]   # the station's values of v
    if (length(xs) < 100) return(NULL)
    nd <- typical(m, 100)
    nd$st2[] <- st_grid$st[i]
    nd[[v]] <- seq(quantile(xs, 0.05), quantile(xs, 0.95), length.out = 100)
    shape_ci(m, nd, v) %>% mutate(model = st_grid$model[i], covariate = v, station = st_grid$st[i], x = nd[[v]])
  }) %>%
    bind_rows()
  
  for (mn in names(stm)) {
    g <- stcur %>%
      filter(model == mn) %>%
      mutate(covariate = factor(covariate, levels = covs$name, labels = gsub("\n", " ", covs$label))) %>%
      ggplot(aes(x = x, y = est, colour = station)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~ covariate, scales = "free") +
      labs(x = NULL, y = ylab[[runs$response[runs$model == mn]]], colour = "Station", title = paste0(mn, ": by radar station"))
    print(g)
  }
}



##################################################################
# MODEL STATS (saved to out_dir)
# fit = one row per model (cells, % deviance explained, adjusted R2); trm = one row per smooth term (edf, test statistic, p; the statistic is F for the VIR models and Chi.sq for the hot models)
{
  fit <- lapply(names(sums), function(mn) {
    data.frame(model = mn, n = sums[[mn]]$n, dev_expl = round(100 * sums[[mn]]$dev.expl, 1), r2_adj = round(sums[[mn]]$r.sq, 3))
  }) %>%
    bind_rows()
  
  trm <- lapply(names(sums), function(mn) {
    st <- sums[[mn]]$s.table
    data.frame(model = mn, term = rownames(st), edf = round(st[, 1], 2), stat = round(st[, 3], 1), p = signif(st[, 4], 3), row.names = NULL)
  }) %>%
    bind_rows()
  
  print(fit)
  print(trm)
  write.csv(fit, paste0(out_dir, "model_fit_stats.csv"), row.names = F)
  write.csv(trm, paste0(out_dir, "model_term_stats.csv"), row.names = F)
}