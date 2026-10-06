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

#' Standardise and/or centre selected columns across a list of data frames
#'
#' Applies \code{datawizard::standardise()} to user-specified columns in each
#' data frame of a list, adding the result as new columns (or overwriting the
#' originals) according to \code{append}. Each variable can be fully
#' standardised (\code{scale = TRUE}) or centered only (\code{scale = FALSE}),
#' and a variable can be given both treatments at once. The centering and
#' scaling values used are returned alongside the amended data frames.
#'
#' @param df_list A list of data frames. Names are used for the \code{df_name}
#'                column of \code{params}. If the list is unnamed but written
#'                inline (\code{list(df1, df2)}), the object names are used. A
#'                pre-built unnamed list gives \code{NA} names: name it first,
#'                e.g. \code{list(df1 = df1, df2 = df2)} or
#'                \code{mget(c("df1", "df2"))}.
#' @param vars    Column names to process. Either a character vector (use
#'                \code{scale} to control standardise vs. centre), or a named
#'                list that groups the variables by treatment, e.g.
#'                \code{list(standardise = c("a", "b"), centre = c("c", "d"))}.
#'                Accepted group names are \code{standardise}/\code{standardize}
#'                (centered and scaled) and \code{centre}/\code{center}
#'                (centered only). In the list form a variable may appear in
#'                both groups, giving it both versions; this requires a named
#'                \code{append} so the two outputs get different names. When
#'                \code{vars} is a list, \code{scale} must not be supplied.
#'                Non-numeric columns are skipped with a warning.
#' @param scale   (Character-vector form of \code{vars} only.) Controls full
#'                standardisation vs. centering only:
#'                \itemize{
#'                  \item A single logical (default \code{TRUE}) applied to all \code{vars}.
#'                  \item An unnamed logical vector the same length as \code{vars},
#'                        matched by position.
#'                  \item A named logical vector, e.g. \code{c(height = TRUE, weight = FALSE)}.
#'                        Any variable in \code{vars} not named here defaults to \code{TRUE}.
#'                  \item A character vector of variable names, e.g.
#'                        \code{c("age", "rt")}: only these are standardised and
#'                        all other \code{vars} are centered only. Every name must
#'                        be in \code{vars}, otherwise an error is raised.
#'                }
#'                Passed per variable to the \code{scale} argument of
#'                \code{datawizard::standardise()}.
#' @param two_sd  Logical. If \code{TRUE}, standardised variables are divided
#'                by two SDs (or MADs). Ignored for centered-only versions,
#'                because datawizard would otherwise halve them.
#' @param append  Controls the names of the output columns, following
#'                datawizard's convention:
#'                \itemize{
#'                  \item A single string (default \code{"Scl"}) is used as a
#'                        suffix for every new column, e.g. \code{heightScl}.
#'                  \item \code{TRUE} uses datawizard's default suffix \code{"_z"}.
#'                  \item \code{FALSE} (or \code{""}) overwrites the original columns.
#'                  \item A named character vector gives a separate suffix per
#'                        treatment, e.g.
#'                        \code{c(standardise = "_z", centre = "_c")}. Names are
#'                        \code{standardise}/\code{standardize} and
#'                        \code{centre}/\code{center}.
#'                }
#'                Two outputs may never share a name, so a variable requested
#'                in both groups needs the named form with different suffixes.
#' @param ...     Further arguments passed to \code{datawizard::standardise()}.
#'                Only \code{robust}, \code{weights}, \code{reference} and
#'                \code{verbose} are accepted; anything else is an error, so
#'                that a misspelled argument cannot be silently ignored. Note
#'                that \code{weights} and \code{reference} are applied
#'                identically to every data frame.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{data}{List of data frames (same length/order/names as \code{df_list})
#'                 with the new columns added (or the original columns
#'                 overwritten if \code{append = FALSE}).}
#'     \item{params}{Data frame with one row per variable per treatment per
#'                   data frame: \code{df_index} (position in the list),
#'                   \code{df_name} (list name, \code{NA} if unnamed),
#'                   \code{variable}, \code{new_variable}, \code{scaled},
#'                   \code{robust}, \code{two_sd}, \code{center}, \code{scale}.
#'                   \code{scale} is the effective divisor (already doubled when
#'                   \code{two_sd} was applied), so
#'                   \code{(x - center) / scale} reproduces the new column.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Mixed treatments in one call
#' out <- standardise_list(
#'   list(df1, df2),
#'   vars = list(standardise = c("age", "rt"), centre = "income")
#' )
#'
#' # Both a centred and a standardised version of the same variable
#' out <- standardise_list(
#'   list(df1, df2),
#'   vars   = list(standardise = "age", centre = "age"),
#'   append = c(standardise = "_z", centre = "_c")
#' )   # -> age_z and age_c
#' }
#'
#' @export
standardise_list <- function(df_list, vars, scale = TRUE, two_sd = FALSE,
                             append = "Scl", ...) {
  
  # --- df_list checks ------------------------------------------------------
  if (!is.list(df_list) || !all(vapply(df_list, is.data.frame, logical(1)))) {
    stop("`df_list` must be a list of data frames.")
  }
  if (length(df_list) == 0) {
    stop("`df_list` is empty.")
  }
  
  # --- Build `spec`: one row per (variable, treatment) ---------------------
  group_scale <- c(standardise = TRUE, standardize = TRUE,
                   centre = FALSE, center = FALSE)
  
  if (is.list(vars)) {
    # Grouped list form; the same variable may appear in both groups
    if (!missing(scale)) {
      stop("`scale` cannot be set when `vars` is a list; the group names ",
           "('standardise' / 'centre') determine it.")
    }
    groups <- names(vars)
    if (is.null(groups) || !all(groups %in% names(group_scale))) {
      stop("When `vars` is a list, its names must be from: ",
           "'standardise', 'centre'.")
    }
    # unlist()/data.frame() would silently coerce e.g. numbers, so check first
    if (!all(vapply(vars, function(g) is.null(g) || is.character(g), logical(1)))) {
      stop("Each group in `vars` must be a character vector of column names.")
    }
    spec <- do.call(rbind, lapply(seq_along(vars), function(k) {
      if (length(vars[[k]]) == 0) return(NULL)
      data.frame(variable = vars[[k]], scale = group_scale[[groups[k]]],
                 stringsAsFactors = FALSE)
    }))
    if (is.null(spec)) {
      stop("`vars` contains no variable names.")
    }
    if (anyNA(spec$variable) || any(spec$variable == "")) {
      stop("`vars` must not contain NA or empty names.")
    }
    if (anyDuplicated(spec[c("variable", "scale")])) {
      stop("A variable appears more than once within the same group of `vars`.")
    }
    
  } else {
    # Character-vector form; `scale` controls the treatment
    if (!is.character(vars) || length(vars) == 0) {
      stop("`vars` must be a non-empty character vector (or a named list of them).")
    }
    if (anyNA(vars) || any(vars == "")) {
      stop("`vars` must not contain NA or empty names.")
    }
    if (anyDuplicated(vars)) {
      stop("A variable appears more than once in `vars`: ",
           paste(unique(vars[duplicated(vars)]), collapse = ", "),
           ". To get both a centred and a standardised version of a variable, ",
           "pass `vars` as a list, e.g. list(standardise = \"x\", centre = \"x\").")
    }
    
    # `scale` as a character vector = the variables to standardise; every other
    # variable in `vars` is centered only. Names are validated against `vars`
    # so that a typo cannot silently turn standardisation into centering.
    if (is.character(scale)) {
      unknown <- setdiff(scale, vars)
      if (length(unknown) > 0) {
        stop("`scale` names variables that are not in `vars`: ",
             paste(unknown, collapse = ", "))
      }
      scale <- stats::setNames(vars %in% scale, vars)
    }
    if (!is.logical(scale)) {
      stop("`scale` must be logical (TRUE/FALSE, optionally named) or a ",
           "character vector of the variables to standardise.")
    }
    # An NA would be treated as FALSE downstream (centering only) without warning
    if (anyNA(scale)) {
      stop("`scale` must not contain NA.")
    }
    if (!is.null(names(scale)) && anyDuplicated(names(scale))) {
      stop("`scale` has duplicated names: ",
           paste(unique(names(scale)[duplicated(names(scale))]), collapse = ", "))
    }
    
    # Expand `scale` into a named logical vector aligned with `vars`
    if (length(scale) == 1 && is.null(names(scale))) {
      scale <- stats::setNames(rep(scale, length(vars)), vars)
    } else {
      if (is.null(names(scale))) {
        if (length(scale) != length(vars)) {
          stop("`scale` must be length 1, named, or the same length as `vars`.")
        }
        names(scale) <- vars
      }
      extra <- setdiff(names(scale), vars)
      if (length(extra) > 0) {
        stop("`scale` names variables that are not in `vars`: ",
             paste(extra, collapse = ", "))
      }
      missing_vars <- setdiff(vars, names(scale))
      if (length(missing_vars) > 0) {
        scale <- c(scale, stats::setNames(rep(TRUE, length(missing_vars)), missing_vars))
      }
      scale <- scale[vars]
    }
    spec <- data.frame(variable = vars, scale = unname(scale),
                       stringsAsFactors = FALSE)
  }
  
  # A non-logical two_sd (e.g. "yes") would be treated as FALSE without warning
  if (!is.logical(two_sd) || length(two_sd) != 1 || is.na(two_sd)) {
    stop("`two_sd` must be a single TRUE or FALSE.")
  }
  
  # standardize.numeric() swallows unknown arguments in its own `...`, so a
  # misspelled argument (e.g. `robst = TRUE`) would be silently ignored.
  allowed_dots <- c("robust", "weights", "reference", "verbose")
  dots <- list(...)
  if (length(dots) > 0) {
    dot_names <- names(dots)
    if (is.null(dot_names)) dot_names <- rep("", length(dots))
    bad <- !(dot_names %in% allowed_dots)
    if (any(bad)) {
      stop("Unsupported argument(s) in `...`: ",
           paste(ifelse(dot_names[bad] == "", "<unnamed>", dot_names[bad]),
                 collapse = ", "),
           ". Allowed: ", paste(allowed_dots, collapse = ", "), ".")
    }
  }
  
  # --- Resolve `append` into a suffix for each row of `spec` ---------------
  # Single value: TRUE -> "_z", FALSE or "" -> overwrite, string -> suffix.
  # Named character vector: a separate suffix per treatment.
  spec$treatment <- ifelse(spec$scale, "standardise", "centre")
  
  if (!is.null(names(append))) {
    canon <- c(standardise = "standardise", standardize = "standardise",
               centre = "centre", center = "centre")
    if (!is.character(append) || anyNA(append) || !all(names(append) %in% names(canon))) {
      stop("A named `append` must be a character vector with names from ",
           "'standardise' and 'centre', e.g. c(standardise = \"_z\", centre = \"_c\").")
    }
    names(append) <- canon[names(append)]
    if (anyDuplicated(names(append))) {
      stop("`append` gives more than one suffix for the same treatment.")
    }
    spec$suffix <- unname(append[spec$treatment])
    if (anyNA(spec$suffix)) {
      stop("`append` has no suffix for: ",
           paste(unique(spec$treatment[is.na(spec$suffix)]), collapse = ", "), ".")
    }
  } else if (isTRUE(append)) {
    spec$suffix <- "_z"
  } else if (isFALSE(append)) {
    spec$suffix <- ""
  } else if (is.character(append) && length(append) == 1 && !is.na(append)) {
    spec$suffix <- append
  } else {
    stop("`append` must be TRUE, FALSE, a single string, or a named character ",
         "vector of suffixes (one per treatment).")
  }
  
  spec$new_name <- paste0(spec$variable, spec$suffix)
  
  # Two outputs sharing a name would silently overwrite each other
  if (anyDuplicated(spec$new_name)) {
    stop("These output column names would be produced more than once: ",
         paste(unique(spec$new_name[duplicated(spec$new_name)]), collapse = ", "),
         ". A variable requested in both groups needs a named `append` with ",
         "different suffixes, e.g. append = c(standardise = \"_z\", centre = \"_c\").")
  }
  
  # --- Data frame names ----------------------------------------------------
  # Name of each data frame (NA if unnamed) plus a label used in warnings
  df_names <- names(df_list)
  if (is.null(df_names)) df_names <- rep(NA_character_, length(df_list))
  df_names[!is.na(df_names) & df_names == ""] <- NA_character_
  
  # If the list was written inline, e.g. standardise_list(list(df1, df2), ...),
  # fill any missing names from the symbols used in that call. This cannot
  # work when a pre-built list object is passed in (its element names are
  # not stored anywhere).
  list_expr <- substitute(df_list)
  if (anyNA(df_names) && is.call(list_expr) &&
      identical(list_expr[[1]], quote(list))) {
    inferred <- vapply(
      as.list(list_expr)[-1],
      function(a) if (is.symbol(a)) as.character(a) else NA_character_,
      character(1)
    )
    if (length(inferred) == length(df_names)) {
      df_names[is.na(df_names)] <- inferred[is.na(df_names)]
    }
  }
  df_labels <- ifelse(is.na(df_names), as.character(seq_along(df_list)), df_names)
  
  # --- Apply to each data frame --------------------------------------------
  results <- lapply(seq_along(df_list), function(i) {
    
    # Always read from the untouched original, so that overwriting a column
    # for one treatment can never feed into the other treatment.
    df_orig <- df_list[[i]]
    df      <- df_orig
    
    uvars        <- unique(spec$variable)
    present_vars <- intersect(uvars, names(df_orig))
    absent_vars  <- setdiff(uvars, names(df_orig))
    
    if (length(absent_vars) > 0) {
      warning("Data frame '", df_labels[i], "': columns not found and skipped: ",
              paste(absent_vars, collapse = ", "))
    }
    
    # standardise() silently returns factors/characters/dates unchanged
    is_num <- vapply(present_vars, function(v) is.numeric(df_orig[[v]]), logical(1))
    if (any(!is_num)) {
      warning("Data frame '", df_labels[i], "': non-numeric columns skipped: ",
              paste(present_vars[!is_num], collapse = ", "))
      present_vars <- present_vars[is_num]
    }
    
    rows <- which(spec$variable %in% present_vars)
    pars <- vector("list", length(rows))
    
    for (j in seq_along(rows)) {
      r        <- rows[j]
      v        <- spec$variable[r]
      new_name <- spec$new_name[r]
      do_scale <- spec$scale[r]
      use_2sd  <- do_scale && two_sd
      
      if (new_name != v && new_name %in% names(df_orig)) {
        warning("Data frame '", df_labels[i], "': column '", new_name,
                "' already exists and will be overwritten.")
      }
      
      out <- datawizard::standardise(
        df_orig[[v]],
        scale = do_scale,
        two_sd = use_2sd,
        add_transform_class = FALSE,
        ...
      )
      
      # standardize.numeric() stores the values it used as attributes
      ctr <- attr(out, "center")
      scl <- attr(out, "scale")
      rob <- attr(out, "robust")
      
      # An all-NA/Inf column is returned by datawizard without attributes
      # (.process_std_center() returns NULL), so fall back accordingly.
      # Otherwise center/scale are always numeric (scale = 1 when scale = FALSE).
      if (is.null(ctr))  ctr <- NA_real_
      if (is.null(scl))  scl <- if (do_scale) NA_real_ else 1
      if (is.null(rob))  rob <- FALSE
      if (use_2sd)       scl <- 2 * scl   # effective divisor
      
      df[[new_name]] <- as.vector(out)   # drop attributes from the new column
      
      pars[[j]] <- data.frame(
        df_index     = i,
        df_name      = df_names[i],
        variable     = v,
        new_variable = new_name,
        scaled       = do_scale,
        robust       = rob,
        two_sd       = use_2sd,
        center       = as.numeric(ctr),
        scale        = as.numeric(scl),
        stringsAsFactors = FALSE
      )
    }
    
    list(data = df, params = do.call(rbind, pars))
  })
  
  data <- lapply(results, `[[`, "data")
  if (any(!is.na(df_names))) names(data) <- ifelse(is.na(df_names), "", df_names)
  
  params <- do.call(rbind, lapply(results, `[[`, "params"))
  if (!is.null(params)) rownames(params) <- NULL
  
  list(data = data, params = params)
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
  n_drop <- length(drop)
  
  exact_match <- vapply(drop_parsed, function(dp) {
    hit <- which(vapply(pop_parsed, function(pp) setequal(pp$vars, dp$vars), logical(1)))
    if (length(hit) == 0) NA_integer_ else hit[1]
  }, integer(1))
  
  subset_hits <- function(dp) {
    which(vapply(pop_parsed, function(pp)
      pp$notation == "*" && length(dp$vars) < length(pp$vars) && all(dp$vars %in% pp$vars), logical(1)))
  }
  
  kept <- pop_level
  additions <- character(0)
  handled <- !is.na(exact_match)   # tracks which `drop` entries are fully accounted for
  
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
        piece_set <- exp$sets[[j]]
        if (length(piece_set) == length(pop_term$vars)) next   # the full term itself
        
        # NEW: does another drop request in THIS SAME CALL name this exact piece?
        other_hit <- which(vapply(seq_len(n_drop), function(k)
          k != i && is.na(exact_match[k]) && setequal(drop_parsed[[k]]$vars, piece_set), logical(1)))
        if (length(other_hit) > 0) { handled[other_hit] <- TRUE; next }   # consumed there - don't retain
        
        if (any(vapply(kept_parsed, term_implies, logical(1), candidate = piece_set))) next
        if (!surgical) next
        if (!(exp$pieces[j] %in% kept) && !(exp$pieces[j] %in% additions)) additions <- c(additions, exp$pieces[j])
      }
    }
  }
  
  unresolved <- character(0)
  for (i in seq_along(drop)) {
    if (handled[i]) next
    hit <- subset_hits(drop_parsed[[i]])
    if (length(hit) == 0) next   # genuinely not found - falls to truly_missing below
    
    if (respect_hierarchy) { unresolved <- c(unresolved, drop[i]); handled[i] <- TRUE; next }
    
    if (length(hit) > 1) warning("'", drop[i], "' matches multiple '*' terms; using the first.")
    target_str <- pop_level[hit[1]]
    exp <- full_expansion(pop_parsed[[hit[1]]])
    exclude <- vapply(exp$sets, setequal, logical(1), y = drop_parsed[[i]]$vars)
    remaining_pieces <- exp$pieces[!exclude]
    remaining_sets   <- exp$sets[!exclude]
    
    other_terms <- lapply(kept[kept != target_str], parse_term)
    final_pieces <- remaining_pieces[!vapply(remaining_sets, function(s)
      any(vapply(other_terms, term_implies, logical(1), candidate = s)), logical(1))]
    
    kept <- c(kept[kept != target_str], final_pieces)
    handled[i] <- TRUE
    message("respect_hierarchy = FALSE: de-bundled '", target_str, "' -> kept: ",
            paste(final_pieces, collapse = ", "), " (removed only '", drop[i], "')")
  }
  
  truly_missing <- drop[!handled]
  if (length(truly_missing) > 0) warning("Not found in pop_level, NOT removed: ", paste(truly_missing, collapse = ", "))
  if (length(unresolved) > 0) {
    warning("Not removed (respect_hierarchy = TRUE, default): ", paste(unresolved, collapse = ", "),
            " would require partially de-bundling a '*' term into a non-hierarchical formula. ",
            "Pass respect_hierarchy = FALSE to allow this explicitly.")
  }
  
  kept <- c(kept, additions)
  if (length(additions) > 0) message("Auto-retained implied term(s): ", paste(additions, collapse = ", "))
  kept
}