# plotting_helpers.R
require(ggplot2)

# Rescale ggplot tick marks for variables --------------------------------------

rescaleTicks <- function(dataVar, dataVarScl, sep=3, ax='x', sd.scale=FALSE){

  data_mean <- mean(dataVar, na.rm=TRUE)
  scale_factor <- if (sd.scale) sd(dataVar, na.rm=TRUE) else 1
  
  dataVar_min <- min(dataVar, na.rm=TRUE)
  dataVar_max <- max(dataVar, na.rm=TRUE)
  
  ticks <- seq(floor(dataVar_min), ceiling(dataVar_max), by=sep)
  
  breaks <- (ticks - data_mean) / scale_factor
  
  if (ax == 'x'){
    outAx <- scale_x_continuous(breaks=breaks, labels=ticks)
  } else if (ax == 'y'){
    outAx <- scale_y_continuous(breaks=breaks, labels=ticks)
  } else {
    stop("Invalid axis argument")
  }
  
  return(outAx) 
}

# Helpers for working with the `params` table returned by standardise_list().
# Columns used: variable, new_variable, center, scale.
# Convention:  model-scale value = (original - center) / scale
#
# Pass a `params` table that has been filtered to ONE data frame (the one the
# model was fitted on); otherwise a variable matches several rows and
# get_scaling() stops with an error rather than guessing.


# --- Lookup ------------------------------------------------------------------

get_scaling <- function(params, var) {
  p <- params[params$new_variable == var, , drop = FALSE]
  if (nrow(p) != 1) {
    stop("Expected exactly one row for '", var, "' in `params`, found ", nrow(p),
         ". Filter `params` to a single data frame first.")
  }
  list(variable = p$variable, center = p$center, scale = p$scale)
}


# --- Converting vectors --------------------------------------------------------

# model scale -> original units
to_original <- function(x, var, params) {
  s <- get_scaling(params, var)
  x * s$scale + s$center
}

# original units -> model scale
to_model <- function(x, var, params) {
  s <- get_scaling(params, var)
  (x - s$center) / s$scale
}


# --- Converting data frames ------------------------------------------------------

# For every model-scale column of `data` that appears in `params`, add a column
# in original units, named by the ORIGINAL variable name. If that column already
# exists (e.g. raw data), it is checked against the back-transformed values
# instead of being overwritten, so a wrong `params` table fails loudly.
# Call this BEFORE converting any model-scale column to a factor.
add_original_units <- function(data, params,
                               vars = intersect(params$new_variable, names(data))) {
  for (v in vars) {
    s   <- get_scaling(params, v)
    new <- data[[v]] * s$scale + s$center
    if (s$variable %in% names(data)) {
      if (!isTRUE(all.equal(data[[s$variable]], new, check.attributes = FALSE))) {
        stop("'", s$variable, "' already exists in `data` but does not match '", v,
             "' back-transformed with `params`. Is `params` for the right data frame?")
      }
    } else {
      data[[s$variable]] <- new
    }
  }
  data
}


# --- Labels and factors ----------------------------------------------------------

# Consistent number formatting for labels and legends
fmt_num <- function(x, digits = 2) {
  format(round(x, digits), nsmall = digits, trim = TRUE)
}

# Numeric -> factor with levels in NUMERIC order, labelled in original units.
# (factor(round(x, 2)) would sort the labels as text, e.g. "-0.20" < "-0.40".)
num_factor <- function(x, digits = 2) {
  lv   <- sort(unique(x))
  labs <- fmt_num(lv, digits)
  if (anyDuplicated(labs)) {
    stop("Rounding to ", digits, " dp makes the factor labels non-unique; ",
         "increase `digits`.")
  }
  factor(x, levels = lv, labels = labs)
}


# --- Prediction grid pre-flight check ---------------------------------------------

# marginaleffects crosses every combination in `variables` with every row of
# `newdata`, so memory scales with combinations x rows. Print it before running.
check_grid <- function(variables, newdata, max_rows = Inf) {
  n_combo <- prod(lengths(variables))
  n_rows  <- n_combo * nrow(newdata)
  message(sprintf("Prediction grid: %s combinations x %s rows = %s rows",
                  format(n_combo, big.mark = ","),
                  format(nrow(newdata), big.mark = ","),
                  format(n_rows, big.mark = ",")))
  if (n_rows > max_rows) {
    stop("Grid has ", format(n_rows, big.mark = ","), " rows, above max_rows = ",
         format(max_rows, big.mark = ","), ".")
  }
  invisible(n_rows)
}


# --- Axis-side alternative ----------------------------------------------------------

# Direct replacement for rescaleTicks(): keeps the data on the model scale but
# labels the axis in original units, using the stored centre/scale rather than
# recomputing mean()/sd(). Prefer add_original_units() when you can.
scale_orig <- function(var, params, ax = c("x", "y"), digits = 2, n = 5, by = NULL, ...) {
  ax <- match.arg(ax)
  s  <- get_scaling(params, var)
  breaks_fn <- function(lims) {
    raw <- lims * s$scale + s$center
    tk  <- if (is.null(by)) pretty(raw, n = n) else
      seq(ceiling(min(raw) / by) * by, floor(max(raw) / by) * by, by = by)
    (tk - s$center) / s$scale
  }
  labels_fn <- function(b) fmt_num(b * s$scale + s$center, digits)
  f <- if (ax == "x") ggplot2::scale_x_continuous else ggplot2::scale_y_continuous
  f(breaks = breaks_fn, labels = labels_fn, ...)
}

# to_model_cols ----------------
to_model_cols <- function(data, params) {
  for (v in intersect(params$new_variable, names(data))) {
    s <- get_scaling(params, v)
    if (s$variable %in% names(data)) data[[v]] <- (data[[s$variable]] - s$center) / s$scale
  }
  data
}

# hi_lo ------------
# Median split of the raw data, labelled by the 10% / 90% quantile values
hi_lo <- function(x, q) {
  if (length(q) != 2 || any(q < 0 | q > 1) || q[1] >= q[2]) {
    stop("`q` must be a numeric vector of length 2 with 0 <= q[1] < q[2] <= 1")
  }
  
  q <- quantile(x, c(q[1], q[2]), na.rm = TRUE)
  num_factor(ifelse(x <= median(x, na.rm = TRUE), q[[1]], q[[2]]))
}


# plot_slice ---------------------------
# Plot a slice of the prediction data
plot_slice <- function(pred, data, x, y, colour_v, facet_v, facet_labeller,
                       q = c(0.3, 0.9),
                       option = "viridis",
                       colour_title = colour_v,
                       colour_breaks = ggplot2::waiver()) {
  pts <- data
  pts[[facet_v]] <- hi_lo(pts[[facet_v]], q)
  
  stopifnot(
    "facet levels in `data` not found in `pred`; check `q` and `num_factor` digits" =
      all(levels(pts[[facet_v]]) %in% levels(pred[[facet_v]]))
  )
  
  ggplot(pred) +
    
    geom_point(data = pts,
               aes(x = .data[[x]], y = .data[[y]], colour = .data[[colour_v]]),
               alpha = 0.55, size = 2) +
    scale_color_viridis(option = option, discrete = FALSE, begin = 0.2, end = 0.8,
                        breaks = colour_breaks, labels = fmt_num,
                        guide = guide_legend(title = colour_title, reverse = TRUE)) +
    ggnewscale::new_scale_colour() +
    
    geom_ribbon(aes(x = .data[[x]], ymin = .data[["conf.low"]], ymax = .data[["conf.high"]],
                    fill = .data[[colour_v]]), alpha = 0.15) +
    geom_line(aes(x = .data[[x]], y = .data[["estimate"]], colour = .data[[colour_v]]),
              linewidth = 1) +
    scale_color_viridis(option = option, discrete = TRUE, begin = 0.2, end = 0.8,
                        guide = guide_legend(title = colour_title, reverse = TRUE)) +
    scale_fill_viridis(option = option, discrete = TRUE, begin = 0.2, end = 0.8,
                       guide = guide_legend(title = colour_title, reverse = TRUE)) +

    facet_grid(rows = vars(.data[[facet_v]]), labeller = facet_labeller)
}


# save_fig ----------------
# Save to file (only if `saveplots` is TRUE). `plot` defaults to the last plot,
# ggplot2::get_last_plot(), exactly as in ggsave().
save_fig <- function(file_name, width, height, out_path, sub_ext = TRUE, extensions = c('svg', 'pdf'), plot = NULL) {

  if (is.null(plot)) plot <- ggplot2::get_last_plot()
  if (is.null(plot)) stop("No `plot` supplied and no last plot available.")
  for (ext in extensions) {
    if (sub_ext) {
      if (!dir.exists(file.path(out_path, ext))){dir.create(file.path(outFigPath, ext))}
      sub_folder <- ext
    }
    ggsave(filename = paste0(file_name, ".", ext), plot = plot, width = width, height = height,
           path = file.path(out_path, sub_folder))
  }
  invisible(plot)
}