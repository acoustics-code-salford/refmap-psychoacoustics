# modelling_helpers.R
require(datawizard)

# Beta distribution transform from arbitrary scale to (0, 1) or back again -----

betaTransform <- function(x, scaleMin, scaleMax,
                          direction=c('forward', 'reverse')[1],
                          squeeze=c('none', 'basic', 'smithson')[1],
                          basic_val=1e-6){
  "Squeeze algorithm 'smithson' is from:
  Smithson, M & Verkuilen, J, 2006. A better lemon squeezer? Maximum-likelihood
  regression with beta-distributed dependent variables. Psychological Methods,
  11(1), 54-71."
  
  if (direction == 'forward') {
    # transform from original scale to (0,1)
    x_scale <- (x - scaleMin) / (scaleMax - scaleMin)
    
    if (squeeze == 'smithson') {
      
      N = length(x_scale)
      x_out <- (x_scale*(N - 1) + 0.5) / N
      
    } else if (squeeze == 'basic') {
      
      x_out <- x_scale*(1 - 2*basic_val) + basic_val
      
    } else {
      
      x_out <- x_scale
      
    }
    
  } else if (direction == 'reverse') {
    # transform from (0, 1) back to original scale
    x_scale <- x * (scaleMax - scaleMin) + scaleMin
    
    if (squeeze == 'smithson') {
      
      N = length(x_scale)
      x_out <- (x_scale*N - 0.5)/(N - 1)
      
    } else if (squeeze == 'basic') {
      
      x_out <- (x_scale - basic_val)/(1 - 2*basic_val)
      
    } else {
      
      x_out <- x_scale
      
    }
    
  } else {
    
    stop("Invalid direction argument: must be 'forward' or 'reverse'")
    
  }
  
  return(x_out)
  
}



# Scaling multiple dataframes --------------------------------------------------

#' Standardise selected columns across a list of data frames
#'
#' Applies \code{datawizard::standardise()} to user-specified columns in each
#' data frame of a list, adding the result as new columns with "Scl"
#' appended to the original column name. Standardisation can be switched on
#' or off per variable via the \code{scale} argument.
#'
#' @param df_list A list of data frames.
#' @param vars    Character vector of column names to (optionally) standardise.
#' @param scale   Controls which variables get full standardisation vs.
#'                centering only:
#'                \itemize{
#'                  \item A single logical (default \code{TRUE}) applied to all \code{vars}.
#'                  \item An unnamed logical vector the same length as \code{vars},
#'                        matched by position.
#'                  \item A named logical vector, e.g. \code{c(height = TRUE, weight = FALSE)}.
#'                        Any variable in \code{vars} not named here defaults to \code{TRUE}.
#'                }
#'                Variables with \code{scale = TRUE} are run through
#'                \code{datawizard::standardise()} (centered AND scaled to SD = 1).
#'                Variables with \code{scale = FALSE} are run through
#'                \code{datawizard::center()} instead, so they are still
#'                mean-centered but left on their original scale (SD unchanged).
#'                Either way the new "Scl" column is always created.
#' @param ...     Additional arguments passed on to \code{datawizard::standardise()}
#'                and \code{datawizard::center()} (e.g. \code{robust}).
#'                Note: only arguments valid for both functions can safely be
#'                passed this way (see Details).
#'
#' @return A list of data frames (same length/order as \code{df_list}), each
#'         with new "<var>Scl" columns added.
#'
#' @examples
#' \dontrun{
#' library(datawizard)
#'
#' df1 <- data.frame(height = rnorm(10, 170, 10), weight = rnorm(10, 70, 15))
#' df2 <- data.frame(height = rnorm(10, 165, 8),  weight = rnorm(10, 65, 12))
#'
#' dfs <- list(df1, df2)
#'
#' # Standardise both variables in every data frame
#' out <- standardise_list(dfs, vars = c("height", "weight"))
#'
#' # Fully standardise height, but only center weight (weightScl keeps its
#' # original SD, just shifted to a mean of 0)
#' out2 <- standardise_list(dfs, vars = c("height", "weight"),
#'                           scale = c(height = TRUE, weight = FALSE))
#' }
#'
#' @export
standardise_list <- function(df_list, vars, scale = TRUE, ...) {
  
  # --- Input checks --------------------------------------------------------
  if (!is.list(df_list) || !all(vapply(df_list, is.data.frame, logical(1)))) {
    stop("`df_list` must be a list of data frames.")
  }
  if (!is.character(vars) || length(vars) == 0) {
    stop("`vars` must be a non-empty character vector of column names.")
  }
  if (!is.logical(scale)) {
    stop("`scale` must be logical (TRUE/FALSE), optionally named.")
  }
  
  # --- Expand `scale` into a named logical vector aligned with `vars` -----
  if (length(scale) == 1) {
    scale <- stats::setNames(rep(scale, length(vars)), vars)
  } else {
    if (is.null(names(scale))) {
      if (length(scale) != length(vars)) {
        stop("`scale` must be length 1, named, or the same length as `vars`.")
      }
      names(scale) <- vars
    }
    # Any vars not explicitly named in `scale` default to TRUE
    missing_vars <- setdiff(vars, names(scale))
    if (length(missing_vars) > 0) {
      scale <- c(scale, stats::setNames(rep(TRUE, length(missing_vars)), missing_vars))
    }
    scale <- scale[vars]  # align order, drop any extras not in `vars`
  }
  
  # --- Apply to each data frame --------------------------------------------
  lapply(df_list, function(df) {
    
    present_vars <- intersect(vars, names(df))
    absent_vars  <- setdiff(vars, names(df))
    
    if (length(absent_vars) > 0) {
      warning("Columns not found and skipped: ", paste(absent_vars, collapse = ", "))
    }
    
    for (v in present_vars) {
      new_name <- paste0(v, "Scl")
      
      if (isTRUE(scale[[v]])) {
        df[[new_name]] <- datawizard::standardise(df[[v]], ...)
      } else {
        df[[new_name]] <- datawizard::center(df[[v]], ...)
      }
    }
    
    df
  })
}

# Prune terms ------------------------
prune_terms <- function(pop_level, drop, sanitize_fn = identity, respect_hierarchy = TRUE) {
  
  parse_term <- function(term) {
    sep <- if (grepl("*", term, fixed = TRUE)) "*"
    else if (grepl(":", term, fixed = TRUE)) ":"
    else NULL
    raw <- if (is.null(sep)) trimws(term) else trimws(strsplit(term, sep, fixed = TRUE)[[1]])
    list(vars_raw = raw, vars = sanitize_fn(raw), notation = if (is.null(sep)) "main" else sep)
  }
  term_implies <- function(term, candidate) {
    if (term$notation == "*") length(candidate) <= length(term$vars) && all(candidate %in% term$vars)
    else setequal(term$vars, candidate)
  }
  full_expansion <- function(term) {
    n <- length(term$vars)
    idx_sets <- unlist(lapply(seq_len(n), function(k) utils::combn(n, k, simplify = FALSE)), recursive = FALSE)
    list(pieces = vapply(idx_sets, function(idx)
      if (length(idx) == 1) term$vars_raw[idx] else paste(term$vars_raw[idx], collapse = ":"), character(1)),
      sets   = lapply(idx_sets, function(idx) term$vars[idx]))
  }
  
  pop_parsed  <- lapply(pop_level, parse_term)
  drop_parsed <- lapply(drop, parse_term)
  
  exact_match <- vapply(drop_parsed, function(dp) {
    hit <- which(vapply(pop_parsed, function(pp) setequal(pp$vars, dp$vars), logical(1)))
    if (length(hit) == 0) NA_integer_ else hit[1]
  }, integer(1))
  
  subset_match <- lapply(seq_along(drop), function(i) {
    if (!is.na(exact_match[i])) return(NULL)
    dp <- drop_parsed[[i]]
    hit <- which(vapply(pop_parsed, function(pp)
      pp$notation == "*" && length(dp$vars) < length(pp$vars) && all(dp$vars %in% pp$vars), logical(1)))
    if (length(hit) == 0) NULL else hit
  })
  
  kept <- pop_level
  additions <- character(0)
  unresolved <- character(0)
  
  # --- exact-match removals, retaining implied lower-order pieces not covered elsewhere
  remove_idx <- unique(stats::na.omit(exact_match))
  if (length(remove_idx) > 0) {
    kept <- pop_level[-remove_idx]
    kept_parsed <- pop_parsed[-remove_idx]
    for (i in seq_along(drop)) {
      if (is.na(exact_match[i])) next
      pop_term <- pop_parsed[[exact_match[i]]]
      if (pop_term$notation != "*") next
      surgical <- drop_parsed[[i]]$notation != "*"
      exp <- full_expansion(pop_term)
      for (j in seq_along(exp$pieces)) {
        if (length(exp$sets[[j]]) == length(pop_term$vars)) next
        if (any(vapply(kept_parsed, term_implies, logical(1), candidate = exp$sets[[j]]))) next
        if (!surgical) next
        if (!(exp$pieces[j] %in% kept) && !(exp$pieces[j] %in% additions)) additions <- c(additions, exp$pieces[j])
      }
    }
  }
  
  # --- subset-match (de-bundling) requests
  for (i in seq_along(drop)) {
    if (!is.na(exact_match[i]) || is.null(subset_match[[i]])) next
    if (respect_hierarchy) { unresolved <- c(unresolved, drop[i]); next }
    
    if (length(subset_match[[i]]) > 1) {
      warning("'", drop[i], "' matches multiple '*' terms; using the first.")
    }
    target_str <- pop_level[subset_match[[i]][1]]
    exp <- full_expansion(pop_parsed[[subset_match[[i]][1]]])
    exclude <- vapply(exp$sets, setequal, logical(1), y = drop_parsed[[i]]$vars)
    remaining_pieces <- exp$pieces[!exclude]
    remaining_sets   <- exp$sets[!exclude]
    
    other_terms <- lapply(kept[kept != target_str], parse_term)
    final_pieces <- remaining_pieces[!vapply(remaining_sets, function(s)
      any(vapply(other_terms, term_implies, logical(1), candidate = s)), logical(1))]
    
    kept <- c(kept[kept != target_str], final_pieces)
    message("respect_hierarchy = FALSE: de-bundled '", target_str, "' -> kept: ",
            paste(final_pieces, collapse = ", "), " (removed only '", drop[i], "')")
  }
  
  truly_missing <- drop[is.na(exact_match) & vapply(subset_match, is.null, logical(1))]
  if (length(truly_missing) > 0) warning("Not found in pop_level, NOT removed: ", paste(truly_missing, collapse = ", "))
  if (length(unresolved) > 0) {
    warning("Not removed (respect_hierarchy = TRUE, default): ", paste(unresolved, collapse = ", "),
            " would require partially de-bundling a '*' term into a non-hierarchical formula. ",
            "Pass respect_hierarchy = FALSE to allow this explicitly.")
  }
  
  kept <- c(kept, additions)   # <-- the fix: additions was never reaching the return value before this
  if (length(additions) > 0) message("Auto-retained implied term(s): ", paste(additions, collapse = ", "))
  kept
}