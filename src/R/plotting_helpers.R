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
scale_orig <- function(var, params, ax = c("x", "y"), digits = 2, n = 5, ...) {
  ax <- match.arg(ax)
  s  <- get_scaling(params, var)
  breaks_fn <- function(lims) {
    (pretty(lims * s$scale + s$center, n = n) - s$center) / s$scale
  }
  labels_fn <- function(b) fmt_num(b * s$scale + s$center, digits)
  f <- if (ax == "x") ggplot2::scale_x_continuous else ggplot2::scale_y_continuous
  f(breaks = breaks_fn, labels = labels_fn, ...)
}