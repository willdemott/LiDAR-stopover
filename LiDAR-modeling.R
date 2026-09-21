# Modeling dependency of stopover intensity on LiDAR-derived Forest Characteristics
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠓⠒⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⣀⠀⠀⠀⠀⠀⢠⢤⣤⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠔⠒⠒⠲⠎⠀⠀⢹⡃⢀⣀⠀⠑⠃⠀⠈⢀⠔⠒⢢⠀⠀⠀⡖⠉⠉⠉⠒⢤⡀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠔⠚⠙⠒⠒⠒⠤⡎⠀⠀⠀⠀⢀⣠⣴⣦⠀⠈⠘⣦⠑⠢⡀⠀⢰⠁⠀⠀⠀⠑⠰⠋⠁⠀⠀⠀⠀⠀⠈⢦⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⣸⠁⠀⠀⠀⠀⠀⠀⢰⠃⠀⣀⣀⡠⣞⣉⡀⡜⡟⣷⢟⠟⡀⣀⡸⠀⡎⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⣻⠀⠀⠀⠀
#⢰⠂⠀⠀⠀⠀⠀⠀⠀⣗⠀⠀⢀⣀⣀⣀⣀⣀⣓⡞⢽⡚⣑⣛⡇⢸⣷⠓⢻⣟⡿⠻⣝⢢⠀⢇⣀⡀⠀⠀⠀⢈⠗⠒⢶⣶⣶⡾⠋⠉⠀⠀⠀⠀⠀
#⠈⠉⠀⠀⠀⠀⠀⢀⠀⠈⠒⠊⠻⣷⣿⣚⡽⠃⠉⠀⠀⠙⠿⣌⠳⣼⡇⠀⣸⣟⡑⢄⠘⢸⢀⣾⠾⠥⣀⠤⠖⠁⠀⠀⠀⢸⡇⠀⠀⠀⠀⠀⢀⠀⠀
#⠀⠀⠀⢰⢆⠀⢀⠏⡇⠀⡀⠀⠀⠀⣿⠉⠀⠀⠀⠀⠀⠀⠀⠈⢧⣸⡇⢐⡟⠀⠙⢎⢣⣿⣾⡷⠊⠉⠙⠢⠀⠀⠀⠀⠀⢸⡇⢀⠀⠀⠀⠀⠈⠣⡀
#⠀⠀⠀⠘⡌⢣⣸⠀⣧⢺⢃⡤⢶⠆⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⣟⠋⢀⠔⣒⣚⡋⠉⣡⠔⠋⠉⢰⡤⣇⠀⠀⠀⠀⢸⡇⡇⠀⠀⠀⠀⠀⠀⠸
#⠀⠀⠀⠀⠑⢄⢹⡆⠁⠛⣁⠔⠁⠀⣿⠀⠀⠀⠀⢸⡇⠀⠀⠀⠀⠀⣿⢠⡷⠋⠁⠀⠈⣿⡇⠀⠀⠀⠈⡇⠉⠀⠀⠀⠀⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠑⣦⡔⠋⠁⠀⠀⠀⣿⠀⠀⢠⡀⢰⣼⡇⠀⡀⠀⠀⣿⠀⠁⠀⠀⠀⠀⣿⣷⠀⠀⠀⠀⡇⠀⠀⢴⣤⠀⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⢰⣿⡇⠀⠀⠀⠀⠀⣿⡀⠀⢨⣧⡿⠋⠀⠘⠛⠀⠀⣿⠀⠀⢀⠀⠀⠀⣿⣿⠀⠀⠀⠀⢲⠀⠀⠀⠀⠀⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⢸⣿⡇⠀⠀⠀⠀⠀⢸⡧⡄⠀⠹⣇⡆⠀⠀⠀⠀⠀⣿⠀⢰⣏⠀⣿⣸⣿⣿⠀⠀⠀⠀⣼⠀⠀⠰⠗⠀⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⢸⣿⡇⠀⠀⠀⠀⠀⢸⡇⣷⣛⣦⣿⢀⠈⠑⠀⢠⡆⣿⠐⢠⣟⠁⢸⠸⣿⣿⢱⣤⢀⠀⣼⠀⠀⢀⠀⠀⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⢸⣿⡇⠀⢀⠀⠀⠀⢸⡇⠘⠫⣟⡇⠊⣣⠘⠛⣾⡆⢿⠀⠙⣿⢀⣘⡃⣿⣿⡏⠉⠒⠂⡿⠀⠰⣾⡄⠀⢸⡟⣽⣀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠸⣿⡇⠀⠘⣾⠀⠀⢸⡇⢸⣇⡙⠣⠀⣹⣇⠀⠈⠧⢀⣀⣀⡏⣸⣿⣇⢹⣿⡇⢴⣴⣄⣀⡀⢰⣿⡇⠀⢸⣇⢿⡿⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠓⠁⠈⠻⢷⠾⠦⠤⠬⣅⣹⣿⣖⣶⣲⣈⡥⠤⠶⡖⠛⠒⠛⠁⠉⠛⠮⠐⢛⡓⠒⢛⠚⠒⠒⠒⠛⣚⣫⡼⠿⠿⣯⠛⠤⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠉⠉⠉⠉⡉⠉⠁⠀⠀⠘⠓⠀⠀⠀⠀⠀⣀⣞⡿⡉⠉⠉⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠹⣶⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀

###################################################################
# Load Packages
{
library(terra)
library(mgcv)
library(mapview)
library(tidyverse)
library(viridisLite)
library(marginaleffects)
library(patchwork)
}

##################################################################
# LiDAR DATA
{
#p95_file <- rast("C:/Users/willd/Box/LiDAR/IA_lidar_4_freda/mosaic_epsg6344_p95_10m.tif")

p95_file <- rast("C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/Northeast_epsg6350_p95_10m_rounded.tif")

nlcd <- rast("C:/Users/willd/Desktop/Everything/GIS 5453/GIS_5453/NLCD_TIFS/Annual_NLCD_LndCov_2024_CU_C1V1.tif")

# crop NLCD in its own CRS then align to the 10 m LiDAR grid
foot <- project(as.polygons(ext(p95_file), crs = crs(p95_file)), crs(nlcd))
nlcd <- crop(nlcd, foot)

# Forest-only
forest <- project(nlcd %in% c(41, 42, 43, 90), p95_file, method = "near", filename = "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/forest.tif", overwrite = T)

# p95 masked at 10 m and aggregated
p95_for <- mask(p95_file, forest, maskvalues = c(0, NA))

# Aggregate
p95      <- aggregate(p95_for, 30, fun = "mean", na.rm = T)
sdp95    <- aggregate(p95_for, 30, fun = "sd",   na.rm = T)
}
##################################################################
# RADAR STATIONS
{
box_dir <- "C:/Users/willd/Box/"
setwd(box_dir)

stations <- unique(read.csv("./DFF_resources/DFF_radar_BCR_list.csv")$radar)
fp <- \(s, season) sprintf("./DFF_ppi_years/%s_2013-24_%s_bioRad_ppi_years_AGL_inverse_Refit.tif", s, season)

sr <- lapply(stations, \(s) rast(fp(s, "spring")))
cen <- vect(data.frame(sta = stations, x = sapply(sr, \(r) mean(ext(r)[1:2])), y = sapply(sr, \(r) mean(ext(r)[3:4]))), geom = c("x", "y"), crs = crs(sr[[1]]))

# station ID layer on the LiDAR grid: each cell owned by its nearest station
foot_p <- as.polygons(ext(p95), crs = crs(p95))
vor <- project(voronoi(cen, bnd = as.polygons(ext(cen) + 1e6, crs = crs(cen))), crs(p95))
sta <- rasterize(vor, p95, field = "sta")

keep <- vor$sta[is.related(vor, foot_p, "intersects")]
mos <- stations[sapply(sr, \(r) !is.null(intersect(ext(r), ext(project(foot_p, crs(r))))))]
}
##################################################################
# Stack, stopover responses, and df for modeling
{
hotc <- \(r) c(sum(rast(Map(\(l, q) l >= q, as.list(r), unlist(global(r, \(v, ...) quantile(v, 0.9, na.rm = T))))), na.rm = T), sum(!is.na(r)))

build <- function(season) {
  yr <- project(mosaic(sprc(lapply(fp(mos, season), rast)), fun = "mean"), p95, method = "bilinear")
  base <- c(mean(yr, na.rm = T), p95, sdp95, sta)
  names(base) <- c("vir", "p95", "sdp95", "sta")
  ok <- noNA(base, falseNA = T)
  yr <- mask(yr, ok)
  
  # range-wide: top 10% of model cells within each year
  hot <- hotc(yr)
  
  # within-station: top 10% of that station's model cells within each year
  hots <- mosaic(sprc(lapply(keep, \(s) hotc(mask(project(rast(fp(s, season)), p95, method = "bilinear"), mask(ok, vor[vor$sta == s, ]))))), fun = "max")
  
  out <- mask(c(base, hot, hots), ok)
  names(out) <- c("vir", "p95", "sdp95", "sta", "hot", "nyr", "hots", "nyrs")
  
  d <- as.data.frame(out, xy = T, na.rm = T)
  d$sta <- factor(d$sta)
  list(r = out, d = d)
}

spring <- build("spring")
autumn <- build("autumn")

spring_d <- spring$d
autumn_d <- autumn$d

grn <- colorRampPalette(c("#c7e9c0", "#74c476", "#31a354", "#00441b"))

# Quick map to visualize our data spatially
vir_cap <- clamp(spring$r[["vir"]], upper = 7000, values = FALSE)
mapview(vir_cap, layer.name = "Spring VIR",
        col.regions = viridis, na.color = "transparent") +
  mapview(spring$r[["p95"]], layer.name = "Canopy height p95 (m)",
          col.regions = grn, na.color = "transparent") +
  mapview(spring$r[["sdp95"]], layer.name = "Canopy height SD (m)",
          col.regions = viridis, na.color = "transparent") +
  mapview(spring$r[["hot"]], layer.name = "Hot years (range-wide)",
          col.regions = magma, na.color = "transparent") +
  mapview(spring$r[["hots"]], layer.name = "Hot years (within-station)",
          col.regions = magma, na.color = "transparent")

nrow(spring_d)
nrow(autumn_d)
}
##################################################################
# Model
{
vir_cap <- 7000

sm <- \(f, fam, d) bam(f, family = fam, data = d, weights = wt, method = "fREML", discrete = T)

rhs <- "s(p95, k = 20) + s(sdp95, k = 20) + s(sta, bs = 're') + s(x, y, k = 800)"
ff <- \(lhs) as.formula(paste(lhs, "~", rhs))
f_vir <- ff("vir")
f_hot <- ff("hot / nyr")
f_hots <- ff("hots / nyrs")

mods <- lapply(list(spring = spring_d, autumn = autumn_d), \(d) {
  dv <- transform(d[d$vir > 0 & d$vir <= vir_cap, ], wt = 1)
  dh <- transform(d[d$nyr > 0, ], wt = nyr)
  ds <- transform(d[d$nyrs > 0, ], wt = nyrs)
  list(vir = sm(f_vir, Gamma(link = "log"), dv), hot = sm(f_hot, binomial, dh), hots = sm(f_hots, binomial, ds))
})

lapply(mods, \(m) lapply(m, summary))
}
#################################################################
# Inspect Models

ylab <- c(vir = "Predicted VIR", hot = "P(hot year)", hots = "P(hot year, within-station)")

pnl <- function(s, k, v) {
  p <- plot_predictions(mods[[s]][[k]], condition = v, type = "response") + labs(x = paste(v, "(m)"), y = ylab[[k]])
  if (k != "vir") p <- p + scale_y_continuous(labels = scales::percent)
  p
}

grd <- expand.grid(v = c("p95", "sdp95"), k = names(ylab), stringsAsFactors = F)

fig_spring <- wrap_plots(Map(\(v, k) pnl("spring", k, v), grd$v, grd$k), ncol = 2) + plot_annotation(title = "Spring")
fig_autumn <- wrap_plots(Map(\(v, k) pnl("autumn", k, v), grd$v, grd$k), ncol = 2) + plot_annotation(title = "Autumn")

fig_spring
fig_autumn

#################################################################
# Model stats

ms <- function(f) do.call(rbind, lapply(names(mods), \(s) do.call(rbind, lapply(names(mods[[s]]), \(k) f(s, k, summary(mods[[s]][[k]]))))))

fit <- ms(\(s, k, m) data.frame(season = s, model = k, n = m$n, dev_expl = round(100 * m$dev.expl, 1), r2_adj = round(m$r.sq, 3)))

# column 3 of s.table is F for the Gamma model and Chi.sq for the binomials
trm <- ms(\(s, k, m) data.frame(season = s, model = k, term = rownames(m$s.table), edf = round(m$s.table[, 1], 2), stat = round(m$s.table[, 3], 1), p = signif(m$s.table[, 4], 3), row.names = NULL))

fit
trm

out_dir <- "C:/Users/willd/Desktop/Everything/WLFW Project/LiDAR-stopover/"
write.csv(fit, paste0(out_dir, "model_fit_stats.csv"), row.names = F)
write.csv(trm, paste0(out_dir, "model_term_stats.csv"), row.names = F)