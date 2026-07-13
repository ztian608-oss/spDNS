#' Transform SF expression to a bounded 0-1 scale
#' @keywords internal
#' @noRd
transform_sf_expression_01 <- function(sf_counts) {
  if (is.null(rownames(sf_counts)) || is.null(colnames(sf_counts)))
    stop("sf_counts must have row and column names.", call. = FALSE)
  lib_size <- Matrix::colSums(sf_counts)
  if (any(!is.finite(lib_size)) || any(lib_size <= 0))
    stop("Every cell must have a positive finite library size.", call. = FALSE)
  sf_data_log <- log1p(sweep(sf_counts, 2, lib_size, "/") * 10000)
  z_cdf01_row <- function(x) {
    x <- as.numeric(x); mu <- mean(x, na.rm = TRUE); sigma <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(mu) || !is.finite(sigma) || sigma == 0) return(rep(0.5, length(x)))
    y <- stats::pnorm((x - mu) / sigma); y[!is.finite(y)] <- 0.5; y
  }
  out <- t(apply(sf_data_log, 1, z_cdf01_row))
  dimnames(out) <- dimnames(sf_counts)
  out
}
