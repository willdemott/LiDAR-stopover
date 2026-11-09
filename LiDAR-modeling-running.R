# Stopover intensity vs LiDAR-derived forest structure: MODELING
# Fits GAMs to the dataframes saved by stopover_lidar_prep.R; makes the effect figures and model stats

###################################################################
# Load Packages
{
  library(mgcv)
  library(tidyverse)
  library(patchwork)
}



##################################################################
# SETUP: directories and parameters
{
  # Where the model dataframes are kept, and model stats are written
  out_dir <- "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/"
  
  # Largest VIR used in the VIR models
  vir_cap <- 7000
  
  # Covariates to plot, with axis labels
  xlab <- c(p95   = "Canopy height\np95 (m)",
            sdp95 = "SD of canopy\nheight (m)",
            rumple = "Rumple index (canopy\nsurface / flat area)",
            fhd   = "Vertical diversity\nFHD (point density-adjusted)",
            sdfhd = "SD of FHD\n(point density-adjusted)",
            nrd   = "Shrub layer: share of\nreturns at 1-3 m (%)",
            chili = "Heat load, CHILI\n(cool to warm)",
            waterdist = "Distance to\nwater (m)",
            beamagl = "Radar beam height\nabove ground (m)",
            elev = "Elevation (m)")
  
  # Season colors (colorblind-friendly)
  cols <- c(spring = "#1b9e77", autumn = "#d95f02")
}



##################################################################
# LOAD DATA: model dataframes (ftype = forest type class and station = nearest radar station, read as factors); FHD and SD of FHD adjusted for point density
{
  spring_d <- read.csv(paste0(out_dir, "spring_d.csv"), stringsAsFactors = T)
  autumn_d <- read.csv(paste0(out_dir, "autumn_d.csv"), stringsAsFactors = T)
  
  # Forest type for the by-type curves: 801, 503, 515, or Other; ordered with Other first, so s(v, by = ft3) adds a difference curve for each of the three types to the shared curve s(v)
  spring_d$ft3 <- ordered(ifelse(substr(spring_d$ftype, 1, 3) %in% c("801", "503", "515"), substr(spring_d$ftype, 1, 3), "Other"), levels = c("Other", "801", "503", "515"))
  autumn_d$ft3 <- ordered(ifelse(substr(autumn_d$ftype, 1, 3) %in% c("801", "503", "515"), substr(autumn_d$ftype, 1, 3), "Other"), levels = c("Other", "801", "503", "515"))
  
  # FHD and SD of FHD adjusted for point density (pd, returns per m2): the part predicted by pd (smooth fit across the model cells) is removed and the overall mean added back,
  # so each value says how a cell's FHD compares with cells of the same point density; pd itself is not a model variable
  spring_d$fhd   <- spring_d$fhd   - fitted(bam(fhd ~ s(pd), data = spring_d))   + mean(spring_d$fhd)
  spring_d$sdfhd <- spring_d$sdfhd - fitted(bam(sdfhd ~ s(pd), data = spring_d)) + mean(spring_d$sdfhd)
  autumn_d$fhd   <- autumn_d$fhd   - fitted(bam(fhd ~ s(pd), data = autumn_d))   + mean(autumn_d$fhd)
  autumn_d$sdfhd <- autumn_d$sdfhd - fitted(bam(sdfhd ~ s(pd), data = autumn_d)) + mean(autumn_d$sdfhd)
  
  # Check: adjusted FHD and SD of FHD should no longer change with point density (correlations near 0)
  print(round(c(fhd = cor(spring_d$fhd, spring_d$pd), sdfhd = cor(spring_d$sdfhd, spring_d$pd)), 3))
  plot(spring_d$pd, spring_d$fhd, pch = ".", xlab = "Point density (returns per m2)", ylab = "FHD (point density-adjusted)")
  
  # Cells per season (each season keeps its own radar coverage, so the counts can differ)
  print(c(spring = nrow(spring_d), autumn = nrow(autumn_d)))
}

# SD of FHD percentiles per point density class (the same check as for FHD)
print(aggregate(spring_d$sdfhd, list(pd_class = cut(spring_d$pd, c(0, 5, 10, 20, 30, 40, Inf))), quantile, c(0.05, 0.5, 0.95)))

# Map of the point density-adjusted SD of FHD (same decile colors as the prep map), to compare with the raw map across the collection boundary
sd_adj <- terra::rast(spring_d[, c("x", "y", "sdfhd")], type = "xyz", crs = "EPSG:6350")
print(mapview(sd_adj, col.regions = viridisLite::turbo, at = unique(quantile(spring_d$sdfhd, seq(0, 1, 0.1))),
                       layer.name = "SD of FHD (point density-adjusted)", na.color = "transparent"))

##################################################################
# MODELS
# vir = mean VIR within the LiDAR year +-1; hot / nyr = proportion of tested years (VIR and a cutoff, within the LiDAR year +-1) in which the cell's VIR reached the top 10% of its nearest station's full range
# ftype = fixed effect of forest type class; ti(p95, fhd) = canopy height x vertical diversity interaction (main effects stay in s(p95) and s(fhd)); s(station, bs = "re") = random effect of the nearest radar station
# nrd = share of all returns at 1-3 m (shrub layer, forest pixels)
# fhd, sdfhd = FHD and SD of FHD adjusted for point density (see LOAD DATA); pd is not in the models
# rumple = rumple index: canopy surface area / flat area of the 10 m p95 surface (1 = flat canopy; higher = more jagged canopy)
# beamagl = height (m) of the nearest radar's lowest beam above the ground (from distance, radar site elevation, and cell elevation): how well the radar sees low-flying birds
# elev = mean elevation (m) of the cell; with beamagl in the models, the elevation curve compares cells with the same beam height above ground
{
  f_vir <- vir ~ ftype +
    s(p95, k = 20) + s(sdp95, k = 20) + s(rumple, k = 20) +
    s(fhd, k = 20) + s(sdfhd, k = 20) +
    ti(p95, fhd) +
    s(nrd, k = 20) +
    s(chili, k = 20) + s(waterdist, k = 20) + 
    s(beamagl, k = 20) + s(elev, k = 20) + s(station, bs = "re")
  
  f_hot <- hot / nyr ~ ftype +
    s(p95, k = 20) + s(sdp95, k = 20) + s(rumple, k = 20) +
    s(fhd, k = 20) + s(sdfhd, k = 20) +
    ti(p95, fhd) +
    s(nrd, k = 20) +
    s(chili, k = 20) + s(waterdist, k = 20) + 
    s(beamagl, k = 20) + s(elev, k = 20) + s(station, bs = "re")
  
  # VIR models: Gamma (log link); cells with VIR > 0 (Gamma needs positive values) and <= vir_cap
  # Hot-year models: binomial; cells with at least 1 year of data; weights = number of years (binomial trials)
  mods <- list(spring = list(vir = bam(f_vir, family = Gamma(link = "log"), 
                                       data = spring_d[spring_d$vir > 0 & spring_d$vir <= vir_cap, ], 
                                       method = "fREML", discrete = T),
                             
                             hot = bam(f_hot, family = binomial, data = spring_d[spring_d$nyr > 0, ], 
                                       weights = nyr, method = "fREML", discrete = T)),
               
               autumn = list(vir = bam(f_vir, family = Gamma(link = "log"), 
                                       data = autumn_d[autumn_d$vir > 0 & autumn_d$vir <= vir_cap, ], 
                                       method = "fREML", discrete = T),
                             
                             hot = bam(f_hot, family = binomial, data = autumn_d[autumn_d$nyr > 0, ], 
                                       weights = nyr, method = "fREML", discrete = T)))
  
  # Summary and concurvity of each model (concurvity = GAM version of collinearity: values near 1 mean a smooth is largely explained by the other terms)
  for (s in names(mods)) {
    for (k in names(mods[[s]])) {
      cat("\n=====", s, k, "=====\n")
      print(summary(mods[[s]][[k]]))
      print(round(concurvity(mods[[s]][[k]], full = T)["worst", ], 2))
    }
  }
}



#################################################################
# EFFECT FIGURES (raw scale), one figure per season: rows 1-2 = VIR model, rows 3-4 = hot-year model; one panel per covariate (5 per row)
# Each line = predicted VIR, or predicted probability that a year is a hot year (top 10% of its nearest station's full range), as one covariate changes and all others are fixed:
#   other covariates at their means, forest type 801 (sugar maple/beech/yellow birch), and the station random effect s(station) left out (= an average station)
# Each line covers the 5th-95th percentile of the covariate's values; rug = observed values
# band = 95% CI of the curve's shape: only the smooths that use the covariate count, as in plot.gam(seWithMean = FALSE); the uncertainty of the average station's level is left out
{
  for (s in c("spring", "autumn")) {
    panels <- list()
    for (k in c("vir", "hot")) {
      m <- mods[[s]][[k]]
      for (v in names(xlab)) {
        # Observed values, every 0.1% from the 5th to 95th percentile (rug)
        q <- quantile(m$model[[v]], seq(0.05, 0.95, by = 0.001))
        
        # Prediction rows: other covariates at their means, forest type 801; v = 100 values over its 5th-95th percentile
        nd <- m$model[rep(1, 100), ]
        for (u in names(xlab)) nd[[u]] <- mean(m$model[[u]])
        nd$ftype[] <- "801 Sugar maple/beech/yellow birch"
        nd[[v]] <- seq(min(q), max(q), length.out = 100)
        
        # Prediction on the link scale (log VIR or log-odds) without s(station)
        p <- predict(m, nd, exclude = "s(station)")
        
        # Standard error of the curve's shape: only the columns of the smooths that use v count (s(v), plus ti(p95, fhd) for p95 and fhd)
        X <- predict(m, nd, type = "lpmatrix")
        X[, -unlist(lapply(m$smooth, function(sm) if (any(sm$term %in% v)) sm$first.para:sm$last.para))] <- 0
        se <- sqrt(rowSums((X %*% vcov(m)) * X))
        
        # Back-transformed (linkinv: exp for VIR, inverse logit for the hot-year probability), with 95% CI
        pr <- data.frame(x = nd[[v]], est = m$family$linkinv(p), lo = m$family$linkinv(p - 1.96 * se), hi = m$family$linkinv(p + 1.96 * se))
        
        g <- ggplot(pr, aes(x = x, y = est)) +
          geom_ribbon(aes(ymin = lo, ymax = hi), fill = cols[[s]], alpha = 0.2) +
          geom_line(colour = cols[[s]], linewidth = 0.8) +
          geom_rug(data = data.frame(x = q), aes(x = x), inherit.aes = F, sides = "b", colour = cols[[s]], alpha = 0.1) +
          scale_x_continuous(n.breaks = 4) +
          labs(x = xlab[[v]], y = if (k == "vir") "Predicted VIR" else "P(hot year)")
        
        panels[[length(panels) + 1]] <- g
      }
    }
    
    print(wrap_plots(panels, ncol = 5) + plot_annotation(title = paste0(s, ": predicted VIR (rows 1-2) and P(hot year) (rows 3-4)")))
  }
}



#################################################################
# INTERACTION FIGURES (raw scale), one figure per season: canopy height (p95) x vertical diversity (FHD)
# Predicted VIR (top row) and P(hot year) (bottom row) vs FHD at low, median, and high canopy height (10th, 50th, 90th percentile of p95, both seasons pooled), one panel per p95 level
# Other covariates at their means, forest type 801, s(station) left out (= an average station); FHD covers the 5th-95th percentile
# band = 95% CI of the surface's shape (s(p95), s(fhd), and ti(p95, fhd) only, as in the effect figures)
{
  pq <- quantile(c(spring_d$p95, autumn_d$p95), c(0.1, 0.5, 0.9), names = F)
  for (s in c("spring", "autumn")) {
    panels_int <- list()
    for (k in c("vir", "hot")) {
      m <- mods[[s]][[k]]
      q <- quantile(m$model$fhd, c(0.05, 0.95))
      
      # Prediction rows: 100 FHD values at each of the 3 canopy heights; other covariates at their means, forest type 801
      nd <- m$model[rep(1, 300), ]
      for (u in names(xlab)) nd[[u]] <- mean(m$model[[u]])
      nd$ftype[] <- "801 Sugar maple/beech/yellow birch"
      nd$p95 <- rep(pq, each = 100)
      nd$fhd <- rep(seq(q[1], q[2], length.out = 100), 3)
      
      # Prediction on the link scale without s(station); standard error of the surface's shape (columns of s(p95), s(fhd), and ti(p95, fhd) only)
      p <- predict(m, nd, exclude = "s(station)")
      X <- predict(m, nd, type = "lpmatrix")
      X[, -unlist(lapply(m$smooth, function(sm) if (any(sm$term %in% c("p95", "fhd"))) sm$first.para:sm$last.para))] <- 0
      se <- sqrt(rowSums((X %*% vcov(m)) * X))
      
      # Back-transformed with 95% CI; panel label of each p95 level, in increasing order
      pr <- data.frame(x = nd$fhd, h = factor(paste0("p95 = ", round(nd$p95, 1), " m"), levels = paste0("p95 = ", round(pq, 1), " m")),
                       est = m$family$linkinv(p), lo = m$family$linkinv(p - 1.96 * se), hi = m$family$linkinv(p + 1.96 * se))
      
      g <- ggplot(pr, aes(x = x, y = est)) +
        geom_ribbon(aes(ymin = lo, ymax = hi), fill = cols[[s]], alpha = 0.2) +
        geom_line(colour = cols[[s]], linewidth = 0.8) +
        facet_wrap(~ h) +
        labs(x = xlab[["fhd"]], y = if (k == "vir") "Predicted VIR" else "P(hot year)")
      
      panels_int[[k]] <- g
    }
    
    print(wrap_plots(panels_int, ncol = 1) + plot_annotation(title = paste0(s, ": canopy height x vertical diversity")))
  }
}



#################################################################
# FOREST TYPE FIGURE: each forest type class vs 801 (sugar maple/beech/yellow birch, the most common class), by season
# Two otherwise identical cells (other covariates at their means, same station) that differ only in forest type
# Left: ratio of predicted VIR; right: odds ratio of a hot year; points = estimates, bars = 95% CIs; 1 = no difference (dashed line)
{
  pr <- data.frame() # estimates, both responses and seasons
  for (k in c("vir", "hot")) {
    for (s in c("spring", "autumn")) {
      m <- mods[[s]][[k]]
      
      # Forest type classes other than 801 (801 = the base class of this figure)
      lv <- setdiff(levels(m$model$ftype), "801 Sugar maple/beech/yellow birch")
      
      # Rows: the base cell (801), then one row per other forest type class; other covariates at their means
      nd <- m$model[rep(1, length(lv) + 1), ]
      for (u in names(xlab)) nd[[u]] <- mean(m$model[[u]])
      nd$ftype[] <- c("801 Sugar maple/beech/yellow birch", lv)
      
      # Difference of each row from the base cell on the link scale (log VIR or log-odds), and its standard error
      X <- predict(m, nd, type = "lpmatrix")
      X <- sweep(X[-1, , drop = F], 2, X[1, ])
      est <- drop(X %*% coef(m))
      se <- sqrt(rowSums((X %*% vcov(m)) * X))
      
      # Back-transformed: VIR ratio or odds ratio, with 95% CI
      pr <- rbind(pr, data.frame(term = lv, est = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se), season = s,
                                 resp = if (k == "vir") "VIR ratio (vs 801)" else "Odds ratio, hot year (vs 801)"))
    }
  }
  
  # Forest type rows sorted by their mean VIR ratio over both seasons (highest at the top); VIR panel on the left
  v <- pr[pr$resp == "VIR ratio (vs 801)", ]
  pr$term <- factor(pr$term, levels = names(sort(tapply(v$est, v$term, mean))))
  pr$resp <- factor(pr$resp, levels = unique(pr$resp))
  
  g <- ggplot(pr, aes(x = est, y = term, colour = season)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5)) +
    facet_wrap(~ resp, scales = "free_x") +
    scale_colour_manual(values = cols) +
    labs(x = NULL, y = NULL, colour = "Season")
  
  print(g)
}



#################################################################
# FOREST TYPE CURVES, one figure per season and response: predicted VIR or P(hot year) vs each forest-structure covariate for forest types 801, 503, and 515
# Separate forest-type models, used only here: each main model plus s(v, by = ft3, k = 10) for the six forest-structure covariates (ft3 ordered, Other first: no curve for Other, one difference curve for each of 801, 503, 515)
# Each type's curve = shared curve s(v) + its difference curve; ftype holds each type's level; other covariates at their means, s(station) left out (= an average station)
# The main models (and every other figure) keep one shared curve per covariate
# Lines cover the 5th-95th percentile of the covariate (all model cells); band = 95% CI of the curve's shape, as in the effect figures (level differences between types: see the forest type figure)
{
  # Forest type labels and colors (same labels and TreeMap colors as the prep maps)
  ft_lab <- c("801" = "801 Sugar maple/beech/yellow birch", "503" = "503 White oak/red oak/hickory", "515" = "515 Chestnut oak/black oak/scarlet oak")
  ft_col <- c("801 Sugar maple/beech/yellow birch" = "#a2c9eb", "503 White oak/red oak/hickory" = "#f22ab1", "515 Chestnut oak/black oak/scarlet oak" = "#eeafa3")
  
  for (s in c("spring", "autumn")) {
    for (k in c("vir", "hot")) {
      # Forest-type model: the main model refit with the difference curves added (same data, family, and settings; slow, like the main fits)
      m <- update(mods[[s]][[k]], . ~ . + s(p95, by = ft3, k = 10) + s(sdp95, by = ft3, k = 10) + s(rumple, by = ft3, k = 10) +
                    s(fhd, by = ft3, k = 10) + s(sdfhd, by = ft3, k = 10) + s(nrd, by = ft3, k = 10))
      pr <- data.frame() # curves of all covariates and types
      for (v in c("p95", "sdp95", "rumple", "fhd", "sdfhd", "nrd")) {
        for (f in names(ft_lab)) {
          # Prediction rows: other covariates at their means, forest type f; v = 100 values over its 5th-95th percentile
          nd <- m$model[rep(1, 100), ]
          for (u in names(xlab)) nd[[u]] <- mean(m$model[[u]])
          nd$ftype[] <- ft_lab[[f]]
          nd$ft3[] <- f
          nd[[v]] <- seq(quantile(m$model[[v]], 0.05), quantile(m$model[[v]], 0.95), length.out = 100)
          
          # Prediction on the link scale without s(station); standard error of the curve's shape (smooths that use v, including the type's difference curve)
          p <- predict(m, nd, exclude = "s(station)")
          X <- predict(m, nd, type = "lpmatrix")
          X[, -unlist(lapply(m$smooth, function(sm) if (any(sm$term %in% v)) sm$first.para:sm$last.para))] <- 0
          se <- sqrt(rowSums((X %*% vcov(m)) * X))
          
          # Back-transformed with 95% CI; covariate label on one line (panel title)
          pr <- rbind(pr, data.frame(covariate = gsub("\n", " ", xlab[[v]]), x = nd[[v]], type = ft_lab[[f]],
                                     est = m$family$linkinv(p), lo = m$family$linkinv(p - 1.96 * se), hi = m$family$linkinv(p + 1.96 * se)))
        }
      }
      
      # One panel per covariate, in the order above (own x and y ranges), one line and band per forest type (legend shows the lines only)
      pr$covariate <- factor(pr$covariate, levels = unique(pr$covariate))
      g <- ggplot(pr, aes(x = x, y = est, colour = type)) +
        geom_ribbon(aes(ymin = lo, ymax = hi, fill = type), alpha = 0.2, colour = NA, show.legend = F) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~ covariate, scales = "free") +
        scale_colour_manual(values = ft_col) +
        scale_fill_manual(values = ft_col, guide = "none") +
        labs(x = NULL, y = if (k == "vir") "Predicted VIR" else "P(hot year)", colour = NULL,
             title = paste0(s, ": ", if (k == "vir") "predicted VIR" else "P(hot year)", " by forest type"))
      
      print(g)
    }
  }
}



#################################################################
# MODEL STATS (saved to out_dir)
# fit = one row per model (sample size, % deviance explained, adjusted R2)
# trm = one row per smooth term (edf, test statistic, p); column 3 of s.table is F for the Gamma models and Chi.sq for the binomial models
{
  fit <- data.frame()
  trm <- data.frame()
  for (s in names(mods)) {
    for (k in names(mods[[s]])) {
      m <- summary(mods[[s]][[k]])
      fit <- rbind(fit, data.frame(model = paste(s, k), n = m$n, dev_expl = round(100 * m$dev.expl, 1), r2_adj = round(m$r.sq, 3)))
      trm <- rbind(trm, data.frame(model = paste(s, k), term = rownames(m$s.table), edf = round(m$s.table[, 1], 2), stat = round(m$s.table[, 3], 1), p = signif(m$s.table[, 4], 3), row.names = NULL))
    }
  }
  
  print(fit)
  print(trm)
  write.csv(fit, paste0(out_dir, "model_fit_stats.csv"), row.names = F)
  write.csv(trm, paste0(out_dir, "model_term_stats.csv"), row.names = F)
}
