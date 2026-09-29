# ---------------------------------------------------------------------------
# helpers.R - shared helpers for the plastisphere resistome case.
#
# The first chunk of plastpath_case.Rmd runs:
#     my_group <- 1             # your group number, 1 to 8
#     source("R/helpers.R")
#
# What you get:
#   env_levels, env_colors      fixed order and colour of the environments for all plots
#   read_metadata()             reads metadata/sample_metadata.csv with the right factor levels
#   select_group(meta, group)   keeps the samples of one group (metadata/group_design.csv)
#   keep_good_hits(x, ...)      identity/coverage filter for KMA hits
#   to_matrix(x, feature, value, samples)   long table -> sample x feature matrix (zeros kept)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

if (!file.exists("metadata/sample_metadata.csv")) {
  stop("Cannot find metadata/sample_metadata.csv. Open plastpath_case.Rproj so that the working directory is the case folder.")
}
if (!exists("my_group")) my_group <- 1

env_levels <- c("Lier river loc1", "Lier river loc2", "Raw WW", "Treated WW", "Blank control")
env_colors <- c("Lier river loc1" = "#2a9d8f", "Lier river loc2" = "#7fb069",
                "Raw WW" = "#d1495b", "Treated WW" = "#3d5a80", "Blank control" = "grey55")

read_metadata <- function() {
  read_csv("metadata/sample_metadata.csv", show_col_types = FALSE) |>
    mutate(environment = factor(environment, levels = env_levels),
           env_type    = factor(env_type, levels = c("River", "Raw WW", "Treated WW", "Control")),
           duration    = factor(paste(duration_weeks, "weeks"), levels = c("2 weeks", "4 weeks")),
           plastic     = factor(plastic, levels = c("PP", "PVC", "PEHD")),
           month       = factor(month, levels = c("June", "August", "September")))
}

# All raw WW and treated WW samples + the river subset of the group.
# A group can have several rows in group_design.csv (e.g. one location in both months).
select_group <- function(meta, group = my_group) {
  design <- read_csv("metadata/group_design.csv", show_col_types = FALSE)
  d <- design[design$group == group, ]
  if (nrow(d) == 0) stop("group must be one of: ", paste(unique(design$group), collapse = ", "))
  meta |>
    filter(env_type %in% c("Raw WW", "Treated WW") |
             paste(environment, month) %in% paste(d$river_site, d$river_month)) |>
    droplevels()
}

# Keep only hits where the CARD reference gene is well covered and well matched.
keep_good_hits <- function(x, min_identity = 80, min_coverage = 80) {
  x |> filter(template_identity >= min_identity, template_coverage >= min_coverage)
}

# Long table -> matrix with samples in rows and features (genes, taxa, ...) in columns.
# Samples in `samples` that have no rows in `x` get a row of zeros (they are NOT dropped).
to_matrix <- function(x, feature = "gene", value = "rpkm", samples = unique(x$sample)) {
  wide <- x |>
    filter(sample %in% samples) |>
    group_by(sample, .data[[feature]]) |>
    summarise(v = sum(.data[[value]]), .groups = "drop") |>
    pivot_wider(names_from = all_of(feature), values_from = v, values_fill = 0)
  m <- as.matrix(wide[, -1, drop = FALSE])
  rownames(m) <- wide$sample
  missing <- setdiff(samples, rownames(m))
  if (length(missing) > 0) {
    m <- rbind(m, matrix(0, nrow = length(missing), ncol = ncol(m), dimnames = list(missing, colnames(m))))
  }
  m[samples, , drop = FALSE]
}

theme_set(theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank()))
