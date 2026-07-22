#' Average predictions with a hybrid CI: response-scale point estimate,
#' link-scale-derived (bounded) interval shape re-centred on that estimate.
#'
#' Wraps marginaleffects::avg_predictions(). Pass the same arguments you would
#' normally pass to avg_predictions() (model, newdata, variables, by, vcov, ...),
#' EXCEPT do not supply 'type' or 'transform' -- those are set internally.
#' 'by' must be an explicit character vector of grouping columns (not by=TRUE),
#' since the two internal calls need to be merged row-for-row.
#'
#' This is a pragmatic display fix, not a formally derived CI -- it does not
#' correspond to inverting a known test statistic. Treat it as such if reporting.
#'
#' @param ... arguments passed through to marginaleffects::avg_predictions()
#' @param recentre if TRUE (default), return the merged/recentred table.
#'   If FALSE, return list(response = ..., link = ...) with both raw calls,
#'   for inspection/debugging.
#' @param est_col,low_col,high_col column names in the avg_predictions() output
#'   holding the point estimate / lower / upper CI. Check names(your_result) for
#'   your installed marginaleffects version if the defaults don't match --
#'   these are typically "estimate", "conf.low", "conf.high" in the underlying
#'   data frame even though the printed table shows "Estimate", "2.5 %", "97.5 %".
#'
#' @return data.frame with the 'by' columns, response-scale 'estimate', and
#'   recentred 'conf.low'/'conf.high'; other columns (p-value, S, etc.) from
#'   the response-scale call are carried through unchanged.
avg_predictions_hybrid <- function(..., recentre = TRUE,
                                   est_col = "estimate",
                                   low_col = "conf.low",
                                   high_col = "conf.high") {
  
  args <- list(...)
  
  if (!is.null(args$type) || !is.null(args$transform)) {
    stop("Do not supply 'type' or 'transform' directly: avg_predictions_hybrid() ",
         "sets these internally (type='response' for the estimate; ",
         "type='link' for the CI shape).")
  }
  
  by_arg <- args$by
  if (is.null(by_arg) || isTRUE(by_arg) || isFALSE(by_arg) || !is.character(by_arg)) {
    stop("avg_predictions_hybrid() requires 'by' to be an explicit character vector ",
         "of column names, so the two internal calls can be merged unambiguously.")
  }
  by_cols <- by_arg
  
  resp <- do.call(marginaleffects::avg_predictions, c(args, list(type = "response")))
  link <- do.call(marginaleffects::avg_predictions, c(args, list(type = "link")))
  
  if (!recentre) return(list(response = resp, link = link))
  
  stopifnot(all(c(by_cols, est_col) %in% names(resp)),
            all(c(by_cols, est_col, low_col, high_col, "std.error") %in% names(link)))
  
  # rename immediately on entry to m, so nothing downstream is ambiguous
  resp_se <- resp[, c(by_cols, "std.error")]
  names(resp_se)[names(resp_se) == "std.error"] <- "std.error.response_delta"
  
  link_renamed <- link[, c(by_cols, est_col, low_col, high_col, "std.error")]
  names(link_renamed)[names(link_renamed) == "std.error"] <- "std.error.link_delta"
  
  m <- merge(resp[, c(by_cols, est_col)], link_renamed,
             by = by_cols, suffixes = c(".resp", ".link"))
  m <- merge(m, resp_se, by = by_cols)
  
  est_resp <- m[[paste0(est_col, ".resp")]]
  est_link <- m[[paste0(est_col, ".link")]]
  
  # critical value implied by marginaleffects' own link-scale interval --
  # inherits whatever conf_level/df was actually used, not hardcoded
  mult <- (m[[high_col]] - est_link) / m$std.error.link_delta
  
  target_link <- qlogis(est_resp)
  
  m$estimate  <- est_resp
  m$conf.low  <- plogis(target_link - mult * m$std.error.link_delta)
  m$conf.high <- plogis(target_link + mult * m$std.error.link_delta)
  
  # diagnostic only -- collapses the asymmetric interval to a single number,
  # using the SAME critical value as the interval itself (mult), so it
  # stays correct at any conf_level rather than assuming 95%
  m$std.error.hybrid_implied <- (m$conf.high - m$conf.low) / (2 * mult)
  m$lower_halfwidth <- m$estimate - m$conf.low
  m$upper_halfwidth <- m$conf.high - m$estimate
  m$asymmetry_ratio  <- m$upper_halfwidth / m$lower_halfwidth
  
  bad <- !is.finite(m$conf.low) | !is.finite(m$conf.high) |
    m$conf.low > m$estimate | m$conf.high < m$estimate
  if (any(bad)) {
    warning(sum(bad), " row(s) failed the sanity check (non-finite bound, or ",
            "estimate outside its own interval) -- inspect these rows.")
  }
  
  drop_cols <- c(by_cols, est_col, low_col, high_col, "std.error")
  extra <- setdiff(names(resp), drop_cols)
  
  out <- m[, c(by_cols, "estimate", "conf.low", "conf.high",
               "std.error.response_delta", "std.error.link_delta",
               "std.error.hybrid_implied",
               "lower_halfwidth", "upper_halfwidth", "asymmetry_ratio")]
  if (length(extra) > 0) {
    out <- merge(out, resp[, c(by_cols, extra)], by = by_cols)
  }
  
  out <- out[match(interaction(resp[by_cols]), interaction(out[by_cols])), ]
  rownames(out) <- NULL
  
  out
}