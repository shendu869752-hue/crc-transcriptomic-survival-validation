suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) != 1L) stop("Cannot resolve script path")
cli <- commandArgs(trailingOnly = TRUE)
root <- if (length(cli) >= 1L && dir.exists(cli[[1L]])) {
  gsub("\\\\", "/", cli[[1L]])
} else {
  gsub("\\\\", "/", getwd())
}
if (!file.exists(file.path(root, "scripts", "build_submission_fig3.R"))) {
  script_file <- gsub("\\\\", "/", sub("^--file=", "", script_arg))
  root <- file.path(dirname(script_file), "..")
}
run_key <- "b2e747538453cf33b4321109514bdc2f06e38e688e14d1ee434c1677b4f6657c"
result_dir <- file.path(root, "runs", run_key, "results")
out_dir <- file.path(root, "submission_package", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scores <- fread(file.path(result_dir, "all_risk_scores.csv"))
stats <- fread(file.path(result_dir, "cohort_performance.csv"))
cohort_order <- data.frame(
  cohort = c("TCGA-COAD", "GSE14333", "GSE17536", "GSE17537"),
  endpoint = c("OS", "DFS", "OS", "OS"),
  stringsAsFactors = FALSE
)

km_data <- scores[scores$cohort %in% cohort_order$cohort, ]
if (!setequal(unique(km_data$cohort), cohort_order$cohort)) {
  stop("Submission Fig3 does not contain exactly the four external cohorts")
}
km_max_year <- max(km_data$time / 365.25)
km_break_by <- if (km_max_year <= 6) 1 else if (km_max_year <= 12) 2 else 5
km_x_limit <- ceiling(km_max_year / km_break_by) * km_break_by
km_x_breaks <- seq(0, km_x_limit, by = km_break_by)

curve_rows <- list()
risk_rows <- list()
for (nm in cohort_order$cohort) {
  d <- km_data[km_data$cohort == nm, ]
  d$group <- factor(d$group, levels = c("Low", "High"))
  fit <- survfit(Surv(time / 365.25, status) ~ group, data = d)
  sm <- summary(fit, censored = TRUE)
  current <- data.frame(
    time_years = sm$time,
    survival = sm$surv,
    n_censor = sm$n.censor,
    group = sub("group=", "", sm$strata),
    cohort = nm,
    endpoint = unique(d$endpoint),
    stringsAsFactors = FALSE
  )
  baseline <- data.frame(
    time_years = 0,
    survival = 1,
    n_censor = 0L,
    group = c("Low", "High"),
    cohort = nm,
    endpoint = unique(d$endpoint),
    stringsAsFactors = FALSE
  )
  curve_rows[[nm]] <- rbind(baseline, current)
  for (risk_group in c("Low", "High")) {
    group_time <- d$time[d$group == risk_group] / 365.25
    risk_rows[[paste(nm, risk_group, sep = "::")]] <- data.frame(
      cohort = nm,
      endpoint = unique(d$endpoint),
      time_years = km_x_breaks,
      group = risk_group,
      n_at_risk = vapply(km_x_breaks, function(tt) sum(group_time >= tt), integer(1)),
      stringsAsFactors = FALSE
    )
  }
}
km_curve <- rbindlist(curve_rows)
km_curve$group <- factor(km_curve$group, levels = c("Low", "High"))
km_risk <- rbindlist(risk_rows)
km_risk$group <- factor(km_risk$group, levels = c("Low", "High"))

ann <- stats[match(cohort_order$cohort, stats$cohort), ]
p_text <- ifelse(ann$logrank_p < 0.001, "P<0.001", sprintf("P=%.3f", ann$logrank_p))
ann$text <- sprintf(
  "n=%d; events=%d\nHR/SD %.2f (%.2f-%.2f)\nlog-rank %s",
  ann$n, ann$events, ann$HR_per_SD, ann$lower95, ann$upper95, p_text
)

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      axis.title = element_text(face = "plain"),
      legend.title = element_blank()
    )
}

curve_plots <- list()
risk_plots <- list()
for (i in seq_len(nrow(cohort_order))) {
  nm <- cohort_order$cohort[[i]]
  endpoint_i <- cohort_order$endpoint[[i]]
  curve_i <- km_curve[km_curve$cohort == nm, ]
  risk_i <- km_risk[km_risk$cohort == nm, ]
  ann_i <- ann[ann$cohort == nm, , drop = FALSE]
  curve_plots[[nm]] <- ggplot(curve_i, aes(time_years, survival, color = group, linetype = group)) +
    geom_step(linewidth = 0.72) +
    geom_point(
      data = curve_i[curve_i$n_censor > 0, ],
      aes(shape = group), size = 0.8, stroke = 0.4, show.legend = FALSE
    ) +
    geom_text(
      data = ann_i,
      aes(x = 0.02 * km_x_limit, y = 0.05, label = text),
      inherit.aes = FALSE, hjust = 0, vjust = 0,
      size = 2.05, family = "Arial", color = "#222222"
    ) +
    scale_color_manual(values = c(Low = "#2166AC", High = "#B2182B")) +
    scale_linetype_manual(values = c(Low = "solid", High = "22")) +
    scale_shape_manual(values = c(Low = 3, High = 4)) +
    scale_x_continuous(
      limits = c(0, km_x_limit), breaks = km_x_breaks,
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.25),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = paste0(nm, " - ", endpoint_i),
      x = NULL, y = "Survival probability", tag = LETTERS[[i]]
    ) +
    theme_pub(8) +
    theme(
      text = element_text(family = "Arial"), legend.position = "none",
      plot.title = element_text(face = "bold", size = 8),
      plot.tag = element_text(face = "bold", size = 9),
      plot.margin = margin(4, 5, 0, 5)
    )
  risk_plots[[nm]] <- ggplot(risk_i, aes(time_years, group, label = n_at_risk, color = group)) +
    geom_text(
      data = risk_i[risk_i$time_years > 0, ],
      size = 2.15, family = "Arial", show.legend = FALSE
    ) +
    geom_text(
      data = risk_i[risk_i$time_years == 0, ],
      aes(x = time_years + 0.12), hjust = 0,
      size = 2.15, family = "Arial", show.legend = FALSE
    ) +
    scale_color_manual(values = c(Low = "#2166AC", High = "#B2182B")) +
    scale_x_continuous(
      limits = c(0, km_x_limit), breaks = km_x_breaks,
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(x = "Time (years)", y = "At risk") +
    theme_pub(7) +
    theme(
      text = element_text(family = "Arial"), legend.position = "none",
      axis.line.y = element_blank(), axis.ticks.y = element_blank(),
      plot.margin = margin(0, 5, 4, 5)
    )
}

draw_figure <- function() {
  grid::grid.newpage()
  outer_layout <- grid::grid.layout(
    nrow = 3L, ncol = 2L,
    heights = grid::unit.c(
      grid::unit(0.22, "inches"),
      grid::unit(1, "null"), grid::unit(1, "null")
    )
  )
  grid::pushViewport(grid::viewport(layout = outer_layout))
  grid::pushViewport(grid::viewport(layout.pos.row = 1L, layout.pos.col = 1:2))
  grid::grid.text(
    "Risk groups", x = 0.26, y = 0.52,
    gp = grid::gpar(fontfamily = "Arial", fontsize = 7, fontface = "bold")
  )
  grid::grid.segments(
    x0 = 0.36, x1 = 0.43, y0 = 0.52, y1 = 0.52,
    gp = grid::gpar(col = "#2166AC", lwd = 1.4, lty = 1)
  )
  grid::grid.text("Low", x = 0.47, y = 0.52, gp = grid::gpar(fontfamily = "Arial", fontsize = 7))
  grid::grid.segments(
    x0 = 0.55, x1 = 0.62, y0 = 0.52, y1 = 0.52,
    gp = grid::gpar(col = "#B2182B", lwd = 1.4, lty = 2)
  )
  grid::grid.text("High", x = 0.67, y = 0.52, gp = grid::gpar(fontfamily = "Arial", fontsize = 7))
  grid::popViewport()
  for (i in seq_len(nrow(cohort_order))) {
    nm <- cohort_order$cohort[[i]]
    grid::pushViewport(grid::viewport(
      layout.pos.row = ((i - 1L) %/% 2L) + 2L,
      layout.pos.col = ((i - 1L) %% 2L) + 1L
    ))
    pair_layout <- grid::grid.layout(
      nrow = 2L, ncol = 1L,
      heights = grid::unit(c(3.2, 1.0), "null")
    )
    grid::pushViewport(grid::viewport(layout = pair_layout))
    print(curve_plots[[nm]], vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 1L))
    print(risk_plots[[nm]], vp = grid::viewport(layout.pos.row = 2L, layout.pos.col = 1L))
    grid::popViewport(2L)
  }
  grid::popViewport()
}

grDevices::cairo_pdf(file.path(out_dir, "Fig3_multicohort_KM.pdf"), width = 7.2, height = 7.4)
draw_figure()
dev.off()
png(
  file.path(out_dir, "Fig3_multicohort_KM.png"),
  width = 7.2, height = 7.4, units = "in", res = 300, type = "cairo-png"
)
draw_figure()
dev.off()

writeLines(
  c(
    "# Submission figure provenance",
    "",
    paste0("Certified source run: `", run_key, "`."),
    "",
    "Fig3 was regenerated from the certified `all_risk_scores.csv` and `cohort_performance.csv` files.",
    "The only intended presentation change from the certified run figure is standard P-value typography (`P<0.001` or `P=...`).",
    "No cohort, endpoint, survival time, event indicator, risk group, hazard ratio, confidence interval, censoring mark, or risk-table count was changed."
  ),
  file.path(out_dir, "README.md")
)

message(out_dir)
