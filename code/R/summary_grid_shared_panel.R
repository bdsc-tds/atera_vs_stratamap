#!/usr/bin/env Rscript
# Figure 5: 5-panel base-R-graphics summary grid [N genes | N cells | nCount/cell | nCount/nGenes/cell
# | ARI], self (solid) vs. shared_panel_genes (hatched) dodge convention throughout, one shared
# y-axis (technology row) across all 5 panels via a single png() device + layout().
#
# Xenium 5k and Xenium MM (biomarkers) don't have a native shared_panel_genes computation (their own
# panels barely overlap the breast+100 gene list), so both reuse their own panel_genes stats as a
# proxy - except panels 3/4/5 (nCount, nCount/nGenes, ARI) deliberately omit Xenium 5k's hatched bar:
# restricting its cells to the shared-panel gene subset filters out enough of them that per-cell
# stats/clustering computed on the remainder would be misleading. Xenium MM's 96-gene panel doesn't
# have this issue, so it keeps all 5 panels. n_genes/n_cells (panels 1/2) stay honest either way and
# show both technologies' proxy.
#
# Also writes summary_grid_data.csv - the exact per-panel values plotted, including which
# (technology, panel) combinations were intentionally left out.
#
# Usage: Rscript summary_grid_shared_panel.R --ngenes-csv <ngenes_barplot.csv> \
#   --ncells-csv <ncells_barplot.csv> --overall-csv <sensitivity_overall_all_technologies.csv> \
#   --ari-csv <separation_ari_comparison.csv> --out-dir <dir>

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
})
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--ngenes-csv", type = "character"),
  make_option("--ncells-csv", type = "character"),
  make_option("--overall-csv", type = "character"),
  make_option("--ari-csv", type = "character"),
  make_option("--out-dir", type = "character")
)))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

techs <- TECH_ORDER
n <- length(techs)
y_pos <- rev(seq_len(n))
ylim <- c(0.3, n + 0.7)
TOP_MAR <- 6
BOT_MAR <- 9
label_levels <- vapply(techs, label_for_tech, character(1))
tech_color <- unname(TECH_FAMILY_COLORS[vapply(techs, tech_family, character(1))])

format_k <- function(x) vapply(x, function(v) {
  if (v >= 1000) paste0(formatC(round(v / 1000), format = "d", big.mark = ","), "K")
  else formatC(v, format = "d")
}, character(1))
nice_log_ticks <- function(max_val) {
  p <- 0:ceiling(log10(max_val))
  ticks <- 10^p
  ticks[ticks <= max_val * 1.05]
}

draw_bar_h <- function(yc, thickness, val, col, hatch) {
  if (is.na(val)) return(invisible())
  rect(0, yc - thickness / 2, val, yc + thickness / 2, col = col, border = "black", lwd = 0.8)
  if (hatch) {
    rect(0, yc - thickness / 2, val, yc + thickness / 2,
         density = 10, angle = 45, col = "white", border = NA, lwd = 1.8)
    rect(0, yc - thickness / 2, val, yc + thickness / 2, col = NA, border = "black", lwd = 0.8)
  }
}

circle_xy <- function(xc, yc, radius_in, n_pts = 120) {
  usr <- par("usr"); pin <- par("pin")
  rx <- radius_in * diff(usr[1:2]) / pin[1]
  ry <- radius_in * diff(usr[3:4]) / pin[2]
  theta <- seq(0, 2 * pi, length.out = n_pts)
  list(x = xc + rx * cos(theta), y = yc + ry * sin(theta))
}
draw_lollipop_h <- function(yc, val, col, hatch, radius_in = 0.10, lwd = 3) {
  if (is.na(val)) return(invisible())
  segments(0, yc, val, yc, col = col, lwd = lwd)
  circ <- circle_xy(val, yc, radius_in)
  polygon(circ$x, circ$y, col = col, border = NA)
  if (hatch) polygon(circ$x, circ$y, col = "white", border = NA, density = 14, angle = 45, lwd = 1.8)
  polygon(circ$x, circ$y, col = NA, border = "black", lwd = lwd * 0.35)
}

# ============================================================ Panels 1/2: N genes / N cells
ngenes_long <- read.csv(opt$`ngenes-csv`, stringsAsFactors = FALSE)
ncells_long <- read.csv(opt$`ncells-csv`, stringsAsFactors = FALSE)
get_long <- function(df, tech, type, value_col) {
  row <- df[df$technology == tech & df$gene_set_type == type, ]
  if (nrow(row) == 0) return(NA_real_)
  row[[value_col]][1]
}
ngenes_self <- vapply(techs, function(t) get_long(ngenes_long, t, "Self genes", "n_genes"), numeric(1))
ngenes_shared <- vapply(techs, function(t) get_long(ngenes_long, t, "Shared panel genes", "n_genes"), numeric(1))
ncells_self <- vapply(techs, function(t) get_long(ncells_long, t, "Self genes", "n_cells"), numeric(1))
ncells_shared <- vapply(techs, function(t) get_long(ncells_long, t, "Shared panel genes", "n_cells"), numeric(1))

draw_panel_ngenes <- function() {
  tf <- function(v) log10(pmax(v, 1))
  xlim <- c(0, max(tf(ngenes_self), na.rm = TRUE) * 1.08)
  par(mar = c(BOT_MAR, 17, TOP_MAR, 1))
  plot(NA, ylim = ylim, xlim = xlim, axes = FALSE, xlab = "", ylab = "", yaxs = "i")
  for (i in seq_len(n)) {
    col <- tech_color[i]; y <- y_pos[i]
    draw_lollipop_h(y + 0.20, tf(ngenes_self[i]), col, hatch = FALSE)
    draw_lollipop_h(y - 0.20, tf(ngenes_shared[i]), col, hatch = TRUE)
  }
  axis(2, at = y_pos, labels = label_levels, las = 1, cex.axis = 2.8, lwd = 0, lwd.ticks = 0)
  xt <- nice_log_ticks(max(ngenes_self, na.rm = TRUE))
  axis(1, at = tf(xt), labels = format_k(xt), cex.axis = 2.5, las = 1)
  mtext("N genes (log10 scale)", side = 3, line = 1.2, cex = 1.5, font = 2)
  box()
}

draw_panel_ncells <- function() {
  xlim <- c(0, max(ncells_self, na.rm = TRUE) * 1.08)
  par(mar = c(BOT_MAR, 1, TOP_MAR, 1))
  plot(NA, ylim = ylim, xlim = xlim, axes = FALSE, xlab = "", ylab = "", yaxs = "i", xaxs = "i")
  for (i in seq_len(n)) {
    col <- tech_color[i]; y <- y_pos[i]
    draw_bar_h(y + 0.20, 0.36, ncells_self[i], col, hatch = FALSE)
    draw_bar_h(y - 0.20, 0.36, ncells_shared[i], col, hatch = TRUE)
  }
  ct <- pretty(c(0, max(ncells_self, na.rm = TRUE)), n = 5)
  ct <- ct[ct >= 0]
  axis(1, at = ct, labels = format_k(ct), cex.axis = 2.5, las = 1)
  mtext("N cells (linear scale)", side = 3, line = 1.2, cex = 1.5, font = 2)
  box()
}

# ============================================================ Panels 3/4: N counts/cell, N counts/N genes/cell
# Xenium MM (biomarkers) reuses its own panel_genes stats as its shared-panel proxy (96-gene panel,
# no real filtering effect). Xenium 5k deliberately gets NO such proxy here - restricting its cells
# to the shared-panel subset drops a meaningful fraction of them (--min-count-in-subset filtering in
# STAGE 4b), so nCount/nCount-per-nGenes computed on that reduced, biased subpopulation would be
# misleading; n_genes/n_cells (panels 1/2 above) stay honest either way and keep the 5k proxy.
overall_all <- read.csv(opt$`overall-csv`, stringsAsFactors = FALSE)
shared_proxy_mm <- overall_all %>% filter(technology == "xenium_biomarkers", gene_set == "panel_genes") %>%
  mutate(gene_set = "shared_panel_genes")
overall_all <- bind_rows(overall_all, shared_proxy_mm)
overall_all <- overall_all %>%
  mutate(iqr = p75 - p25, whisk_lo = pmax(min, p25 - 1.5 * iqr), whisk_hi = pmin(max, p75 + 1.5 * iqr))

get_stats <- function(tech, metric, gene_set) {
  row <- overall_all %>% filter(technology == tech, metric == !!metric, gene_set == !!gene_set)
  if (nrow(row) == 0) return(NULL)
  row
}
draw_box_h <- function(yc, thickness, s, col, hatch, tf) {
  wlo <- tf(s$whisk_lo); lo <- tf(s$p25); med <- tf(s$median); hi <- tf(s$p75); whi <- tf(s$whisk_hi)
  segments(wlo, yc, lo, yc); segments(hi, yc, whi, yc)
  segments(wlo, yc - thickness * 0.3, wlo, yc + thickness * 0.3)
  segments(whi, yc - thickness * 0.3, whi, yc + thickness * 0.3)
  rect(lo, yc - thickness / 2, hi, yc + thickness / 2, col = col, border = "black", lwd = 0.8)
  if (hatch) {
    rect(lo, yc - thickness / 2, hi, yc + thickness / 2,
         density = 10, angle = 45, col = "white", border = NA, lwd = 1.8)
    rect(lo, yc - thickness / 2, hi, yc + thickness / 2, col = NA, border = "black", lwd = 0.8)
  }
  segments(med, yc - thickness / 2, med, yc + thickness / 2, lwd = 2.5)
}

draw_panel_box <- function(metric, log_x, title) {
  all_rows <- overall_all %>% filter(metric == !!metric)
  tf <- if (log_x) function(v) log10(pmax(v, 1)) else identity
  raw_lo <- if (log_x) pmax(all_rows$whisk_lo, 1) else all_rows$whisk_lo
  xlim <- range(c(tf(raw_lo), tf(all_rows$whisk_hi)), na.rm = TRUE)
  par(mar = c(BOT_MAR, 1, TOP_MAR, 1))
  plot(NA, ylim = ylim, xlim = xlim, axes = FALSE, xlab = "", ylab = "", yaxs = "i")
  for (i in seq_len(n)) {
    tech <- techs[i]; col <- tech_color[i]; y <- y_pos[i]
    s_self <- get_stats(tech, metric, "all"); s_shared <- get_stats(tech, metric, "shared_panel_genes")
    if (!is.null(s_self)) draw_box_h(y + 0.20, 0.36, s_self, col, hatch = FALSE, tf = tf)
    if (!is.null(s_shared)) draw_box_h(y - 0.20, 0.36, s_shared, col, hatch = TRUE, tf = tf)
  }
  if (log_x) {
    xt <- nice_log_ticks(max(all_rows$whisk_hi, na.rm = TRUE))
    axis(1, at = tf(xt), labels = formatC(xt, format = "d", big.mark = ","), cex.axis = 2.5, las = 1)
  } else {
    axis(1, cex.axis = 2.5, las = 1)
  }
  for (j in seq_along(title)) mtext(title[j], side = 3, line = 1.2 + (length(title) - j) * 1.4, cex = 1.5, font = 2)
  box()
}

# ============================================================ Panel 5: ARI
ari <- read.csv(opt$`ari-csv`, stringsAsFactors = FALSE)
get_ari <- function(tech, ari_type) {
  row <- ari %>% filter(technology == !!tech, ari_type == !!ari_type)
  if (nrow(row) == 0) return(NA_real_)
  row$ari[1]
}
draw_panel5_ari <- function() {
  xlim <- c(0, max(ari$ari, na.rm = TRUE) * 1.08)
  par(mar = c(BOT_MAR, 1, TOP_MAR, 2))
  plot(NA, ylim = ylim, xlim = xlim, axes = FALSE, xlab = "", ylab = "", yaxs = "i", xaxs = "i")
  for (i in seq_len(n)) {
    tech <- techs[i]; col <- tech_color[i]; y <- y_pos[i]
    draw_bar_h(y + 0.20, 0.36, get_ari(tech, "Self genes"), col, hatch = FALSE)
    # Xenium 5k excluded here - same reasoning as panels 3/4 above (shared-panel restriction drops
    # too many of its cells for the resulting clustering ARI to be meaningful).
    shared_ari <- if (tech == "xenium_5k") NA_real_ else get_ari(tech, "Shared panel genes")
    draw_bar_h(y - 0.20, 0.36, shared_ari, col, hatch = TRUE)
  }
  axis(1, cex.axis = 2.5, las = 1)
  mtext("Bioconservation", side = 3, line = 2.6, cex = 1.5, font = 2)
  mtext("(ARI cell-type vs clustering)", side = 3, line = 1.2, cex = 1.5, font = 2)
  box()
}

# ============================================================ CSV: exactly what's plotted above,
# one row per (technology, panel, gene_set) - so the "5k excluded from panels 3/4/5" choice is
# visible in the data too, not just implied by an absent bar.
csv_rows <- list(
  data.frame(technology = techs, label = label_levels, panel = "n_genes", gene_set = "self", value = ngenes_self),
  data.frame(technology = techs, label = label_levels, panel = "n_genes", gene_set = "shared_panel_genes", value = ngenes_shared),
  data.frame(technology = techs, label = label_levels, panel = "n_cells", gene_set = "self", value = ncells_self),
  data.frame(technology = techs, label = label_levels, panel = "n_cells", gene_set = "shared_panel_genes", value = ncells_shared)
)
for (metric in c("nCount", "nCount_per_nGenes")) {
  for (gs in c("all", "shared_panel_genes")) {
    for (i in seq_len(n)) {
      s <- get_stats(techs[i], metric, gs)
      if (is.null(s)) next
      csv_rows[[length(csv_rows) + 1]] <- data.frame(
        technology = techs[i], label = label_levels[i], panel = metric,
        gene_set = if (gs == "all") "self" else gs,
        value = s$median[1], p25 = s$p25[1], p75 = s$p75[1], whisk_lo = s$whisk_lo[1], whisk_hi = s$whisk_hi[1]
      )
    }
  }
}
for (i in seq_len(n)) {
  tech <- techs[i]
  csv_rows[[length(csv_rows) + 1]] <- data.frame(technology = tech, label = label_levels[i], panel = "ari",
                                                   gene_set = "self", value = get_ari(tech, "Self genes"))
  if (tech != "xenium_5k") {
    csv_rows[[length(csv_rows) + 1]] <- data.frame(technology = tech, label = label_levels[i], panel = "ari",
                                                     gene_set = "shared_panel_genes", value = get_ari(tech, "Shared panel genes"))
  }
}
csv_out <- bind_rows(csv_rows)
write.csv(csv_out, file.path(opt$`out-dir`, "summary_grid_data.csv"), row.names = FALSE)

# ============================================================ Assemble
out_path <- file.path(opt$`out-dir`, "summary_grid_shared_panel.png")
png(out_path, width = 25, height = 8, units = "in", res = 300)
layout(matrix(1:5, nrow = 1), widths = c(1.6, 1, 1, 1, 1))
draw_panel_ngenes()
draw_panel_ncells()
draw_panel_box("nCount", log_x = TRUE, title = c("N counts / cell", "(log10 scale)"))
draw_panel_box("nCount_per_nGenes", log_x = FALSE, title = "N counts / N genes / cell")
draw_panel5_ari()
dev.off()
cat("Saved:", out_path, "\n")
