#' @keywords internal
#' @noRd
with_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (has_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  set.seed(seed)
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  force(expr)
}

#' @keywords internal
#' @noRd
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

#' @keywords internal
#' @noRd
resolve_parallel_cores <- function(n_cores) {
  n_cores <- n_cores %||% 1L
  cores_avail <- parallel::detectCores()
  if (!is.finite(cores_avail) || cores_avail < 1) cores_avail <- n_cores
  max(1L, min(as.integer(n_cores), as.integer(cores_avail)))
}

#' @keywords internal
#' @noRd
discrete_palette_presets <- function() {
  list(
    bright = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"),
    pastel = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5"),
    muted = c("#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#edc949", "#af7aa1", "#ff9da7"),
    vivid = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf"),
    gray = c("#111111", "#555555", "#999999", "#cccccc")
  )
}

#' @keywords internal
#' @noRd
resolve_discrete_palette <- function(n, palette, palette_name) {
  if (!is.null(palette)) {
    if (length(palette) < 1) stop("palette must have at least 1 color.")
    if (length(palette) >= n) return(palette[1:n])
    return(grDevices::colorRampPalette(palette)(n))
  }
  presets <- discrete_palette_presets()
  if (is.null(palette_name)) palette_name <- "bright"
  if (!palette_name %in% names(presets)) {
    stop("Unknown palette_name. Use one of: ", paste(names(presets), collapse = ", "))
  }
  pal <- presets[[palette_name]]
  if (length(pal) >= n) return(pal[1:n])
  grDevices::colorRampPalette(pal)(n)
}
