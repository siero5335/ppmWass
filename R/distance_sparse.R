#' Sparse exact capped-cost transport
#'
#' Since total mass is one, minimize cost by maximizing savings 1-C over
#' partial flows (row sums <= a, column sums <= b). Edges with C=1 contribute
#' no savings. Connected components can therefore be optimized independently.
#' @keywords internal
#' @noRd
sparse_exact_ot <- function(mz_p, mz_q, w_p, w_q, width, tolerance,
                            solver = NULL, diagnostics = FALSE) {
  if (length(width) != 1L || !is.finite(width) || width <= 0) {
    stop("Sparse OT saturation width must be positive and finite.", call. = FALSE)
  }
  ip <- order(mz_p[w_p > 0]); iq <- order(mz_q[w_q > 0])
  a <- cbind(mz_p[w_p > 0][ip], w_p[w_p > 0][ip])
  b <- cbind(mz_q[w_q > 0][iq], w_q[w_q > 0][iq])
  e <- sparse_edges(a[,1], b[,1], width)
  groups <- split(seq_along(e$i), e$component)
  saving <- 0; general <- 0L; attempts <- list()
  if (is.null(solver)) solver <- function(...) transport::transport(...)
  details <- list(edges = length(e$i), components = length(groups),
                  analytic_components = 0L, general_components = 0L)
  for (k in groups) {
    ai <- unique(e$i[k]); bj <- unique(e$j[k])
    if (length(ai) == 1L || length(bj) == 1L) {
      o <- k[order(e$cost[k])]
      cap <- if (length(ai) == 1L) a[ai,2] else b[bj,2]
      weights <- if (length(ai) == 1L) b[e$j[o],2] else a[e$i[o],2]
      flow <- pmin(weights, pmax(0, cap - c(0, utils::head(cumsum(weights), -1))))
      saving <- saving + sum(flow * (1-e$cost[o]))
      details$analytic_components <- details$analytic_components + 1L
    } else {
      general <- general + 1L
      details$general_components <- general
      na <- length(ai); nb <- length(bj)
      C <- matrix(1, na+1L, nb+1L); C[na+1L,] <- 0
      C[cbind(match(e$i[k], ai), match(e$j[k], bj))] <- e$cost[k]
      wa <- c(a[ai,2], sum(b[bj,2])); wb <- c(b[bj,2], sum(a[ai,2]))
      total <- sum(wa); wa <- wa/total; wb <- wb/total
      solver_error <- NA_character_
      plan <- tryCatch(solver(a=wa, b=wb, costm=C), error=function(err) {
        solver_error <<- conditionMessage(err)
        NULL
      })
      assessment <- assess_ot_plan(plan, C, wa, wb, tolerance=tolerance)
      if (!is.na(solver_error)) assessment$validation_error <- "solver_error"
      if (diagnostics) attempts[[length(attempts)+1L]] <- ot_attempt_row(
        paste0("sparse_component_", general), "transport_exact_sparse_component",
        NULL, NULL, assessment, solver_error
      )
      if (!isTRUE(assessment$accepted)) {
        return(list(total=NA_real_, attempts=attempts, details=details))
      }
      real <- plan$from <= na & plan$to <= nb
      saving <- saving + total * sum(plan$mass[real] *
        (1-C[cbind(plan$from[real], plan$to[real])]))
    }
  }
  list(total=1-saving, attempts=attempts, details=details)
}
