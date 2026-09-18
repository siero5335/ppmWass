#' @keywords internal
#' @noRd
is_near_zero <- function(x, tol = .Machine$double.eps) {
  is.finite(x) & abs(x) <= tol
}

#' @keywords internal
#' @noRd
# Reference implementation retained for regression tests against the
# optimized align_spectra() path.
align_spectra_reference <- function(specA, specB, ppm = 20) {
  if (nrow(specA) == 0 || nrow(specB) == 0) {
    return(list(p = numeric(0), q = numeric(0), mz = numeric(0)))
  }

  mzA <- specA[, 1]
  intA <- specA[, 2]
  mzB <- specB[, 1]
  intB <- specB[, 2]

  all_mz <- sort(unique(c(mzA, mzB)))

  n <- length(all_mz)
  merged_mz <- vector("list", n)
  k <- 0
  i <- 1
  while (i <= n) {
    current_mz <- all_mz[i]
    tol <- current_mz * ppm * 1e-6
    j <- i
    while (j <= n && abs(all_mz[j] - current_mz) <= tol) {
      j <- j + 1
    }
    k <- k + 1
    merged_mz[[k]] <- mean(all_mz[i:(j - 1)])
    i <- j
  }
  merged_mz <- unlist(merged_mz[1:k], use.names = FALSE)

  p <- numeric(length(merged_mz))
  q <- numeric(length(merged_mz))

  for (i in seq_along(merged_mz)) {
    mz <- merged_mz[i]
    tol <- mz * ppm * 1e-6

    idxA <- which(abs(mzA - mz) <= tol)
    if (length(idxA) > 0) p[i] <- sum(intA[idxA])

    idxB <- which(abs(mzB - mz) <= tol)
    if (length(idxB) > 0) q[i] <- sum(intB[idxB])
  }

  if (sum(p) > 0) p <- p / sum(p)
  if (sum(q) > 0) q <- q / sum(q)

  list(p = p, q = q, mz = merged_mz)
}

#' @keywords internal
#' @noRd
align_spectra <- function(specA, specB, ppm = 20) {
  if (nrow(specA) == 0 || nrow(specB) == 0) {
    return(list(p = numeric(0), q = numeric(0), mz = numeric(0)))
  }

  mzA <- specA[, 1]
  intA <- specA[, 2]
  mzB <- specB[, 1]
  intB <- specB[, 2]

  ordA <- order(mzA)
  ordB <- order(mzB)
  mzA <- mzA[ordA]
  intA <- intA[ordA]
  mzB <- mzB[ordB]
  intB <- intB[ordB]

  all_mz <- sort(unique(c(mzA, mzB)))

  n <- length(all_mz)
  merged_mz <- vector("list", n)
  k <- 0
  i <- 1
  while (i <= n) {
    current_mz <- all_mz[i]
    tol <- current_mz * ppm * 1e-6
    j <- i
    while (j <= n && abs(all_mz[j] - current_mz) <= tol) {
      j <- j + 1
    }
    k <- k + 1
    merged_mz[[k]] <- mean(all_mz[i:(j - 1)])
    i <- j
  }
  merged_mz <- unlist(merged_mz[1:k], use.names = FALSE)

  tol_vec <- merged_mz * ppm * 1e-6

  csA <- c(0, cumsum(intA))
  csB <- c(0, cumsum(intB))

  lower <- merged_mz - tol_vec
  upper <- merged_mz + tol_vec

  leftA <- findInterval(lower, mzA, left.open = TRUE) + 1
  rightA <- findInterval(upper, mzA)

  leftB <- findInterval(lower, mzB, left.open = TRUE) + 1
  rightB <- findInterval(upper, mzB)

  p <- numeric(length(merged_mz))
  q <- numeric(length(merged_mz))

  okA <- leftA <= rightA
  okB <- leftB <= rightB

  if (any(okA)) {
    p[okA] <- csA[rightA[okA] + 1] - csA[leftA[okA]]
  }
  if (any(okB)) {
    q[okB] <- csB[rightB[okB] + 1] - csB[leftB[okB]]
  }

  if (sum(p) > 0) p <- p / sum(p)
  if (sum(q) > 0) q <- q / sum(q)

  list(p = p, q = q, mz = merged_mz)
}
