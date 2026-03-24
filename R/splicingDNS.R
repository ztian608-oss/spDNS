#' Expand an SF-gene network to SF-splicing event network
#'
#' Converts a directed network of splicing factor (SF) to gene edges into
#' a directed network of SF to splicing event edges via a gene-event mapping.
#'
#' @param sf_gene_network data.frame with at least columns `SF` and `target_gene`.
#' @param gene_event_map data.frame with at least columns `gene` and `splicing_event`.
#'
#' @return data.frame with columns `SF`, `target_gene`, `splicing_event`.
#' @export
expand_network_to_splicing <- function(sf_gene_network, gene_event_map) {
  if (!all(c("SF", "target_gene") %in% colnames(sf_gene_network))) {
    stop("sf_gene_network must contain columns: SF, target_gene")
  }
  if (!all(c("gene", "splicing_event") %in% colnames(gene_event_map))) {
    stop("gene_event_map must contain columns: gene, splicing_event")
  }

  sf_gene_network <- unique(sf_gene_network[, c("SF", "target_gene")])
  gene_event_map <- unique(gene_event_map[, c("gene", "splicing_event")])

  out <- merge(
    sf_gene_network,
    gene_event_map,
    by.x = "target_gene",
    by.y = "gene",
    all.x = FALSE,
    all.y = FALSE
  )

  out <- unique(out[, c("SF", "target_gene", "splicing_event")])
  rownames(out) <- NULL
  out
}

#' Compute single-cell SF activity scores on an SF-event directed network
#'
#' For SF \eqn{s} and cell \eqn{c}, score is:
#' \deqn{A_{s,c} = z(E_{s,c}) \times \frac{\sum_{e \in T(s)} w_{s,e}\tilde{\psi}_{e,c}}{\sum_{e \in T(s)} |w_{s,e}|}}
#' where \eqn{E_{s,c}} is SF expression, \eqn{\tilde{\psi}_{e,c}} is centered PSI,
#' and \eqn{w_{s,e}} is optional edge weight (defaults to 1).
#'
#' @param sf_expression matrix-like SF expression (rows = SF, columns = cells).
#' @param psi matrix-like PSI (rows = splicing events, columns = cells).
#' @param sf_event_network data.frame with columns `SF`, `splicing_event`
#'   and optional `weight`.
#' @param center_psi logical; center each event PSI by row median.
#' @param min_events minimum number of mapped events required for an SF.
#' @param na_to_zero logical; if TRUE, missing PSI is set to 0 after centering.
#'
#' @return matrix SF x cells of activity scores.
#' @export
compute_sf_activity_score <- function(
    sf_expression,
    psi,
    sf_event_network,
    center_psi = TRUE,
    min_events = 3,
    na_to_zero = TRUE
) {
  if (!all(c("SF", "splicing_event") %in% colnames(sf_event_network))) {
    stop("sf_event_network must contain columns: SF, splicing_event")
  }

  sf_expression <- as.matrix(sf_expression)
  psi <- as.matrix(psi)

  common_cells <- intersect(colnames(sf_expression), colnames(psi))
  if (length(common_cells) == 0) {
    stop("sf_expression and psi share no common cell names")
  }
  sf_expression <- sf_expression[, common_cells, drop = FALSE]
  psi <- psi[, common_cells, drop = FALSE]

  net <- unique(sf_event_network[, intersect(c("SF", "splicing_event", "weight"), colnames(sf_event_network)), drop = FALSE])
  if (!"weight" %in% colnames(net)) {
    net$weight <- 1
  }

  keep <- net$SF %in% rownames(sf_expression) & net$splicing_event %in% rownames(psi)
  net <- net[keep, , drop = FALSE]
  if (nrow(net) == 0) {
    stop("No SF-event edges overlap with sf_expression and psi row names")
  }

  sf_degree <- table(net$SF)
  valid_sf <- names(sf_degree)[sf_degree >= min_events]
  net <- net[net$SF %in% valid_sf, , drop = FALSE]

  psi_used <- psi[unique(net$splicing_event), , drop = FALSE]
  if (center_psi) {
    psi_used <- psi_used - matrixStats::rowMedians(psi_used, na.rm = TRUE)
  }
  if (na_to_zero) {
    psi_used[is.na(psi_used)] <- 0
  }

  sf_levels <- sort(unique(net$SF))
  event_levels <- rownames(psi_used)

  edge_sf <- match(net$SF, sf_levels)
  edge_ev <- match(net$splicing_event, event_levels)
  W <- Matrix::sparseMatrix(
    i = edge_sf,
    j = edge_ev,
    x = net$weight,
    dims = c(length(sf_levels), length(event_levels)),
    dimnames = list(sf_levels, event_levels)
  )

  weighted_psi <- W %*% psi_used
  denom <- Matrix::rowSums(abs(W))
  denom[denom == 0] <- NA_real_
  weighted_psi <- weighted_psi / denom

  sf_exp_used <- sf_expression[sf_levels, , drop = FALSE]
  sf_exp_z <- t(scale(t(sf_exp_used)))
  sf_exp_z[is.na(sf_exp_z)] <- 0

  score <- sf_exp_z * as.matrix(weighted_psi)
  score
}

#' Differential SF activity analysis between two conditions
#'
#' @param sf_activity matrix SF x cells from `compute_sf_activity_score`.
#' @param condition vector/factor of condition labels for each cell.
#' @param contrast length-2 character vector specifying comparison order.
#' @param method one of `wilcox` or `t.test`.
#'
#' @return data.frame with per-SF differential statistics.
#' @export
differential_sf_activity <- function(
    sf_activity,
    condition,
    contrast = NULL,
    method = c("wilcox", "t.test")
) {
  method <- match.arg(method)
  sf_activity <- as.matrix(sf_activity)

  if (length(condition) != ncol(sf_activity)) {
    stop("condition length must match number of columns in sf_activity")
  }
  condition <- as.character(condition)
  if (is.null(contrast)) {
    contrast <- unique(condition)
  }
  if (length(contrast) != 2) {
    stop("contrast must contain exactly two condition labels")
  }

  idx_a <- which(condition == contrast[1])
  idx_b <- which(condition == contrast[2])
  if (length(idx_a) < 2 || length(idx_b) < 2) {
    stop("Each condition in contrast must have at least 2 cells")
  }

  res <- lapply(seq_len(nrow(sf_activity)), function(i) {
    x <- sf_activity[i, idx_a]
    y <- sf_activity[i, idx_b]
    if (method == "wilcox") {
      p <- suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value)
    } else {
      p <- stats::t.test(x, y)$p.value
    }
    data.frame(
      SF = rownames(sf_activity)[i],
      mean_a = mean(x, na.rm = TRUE),
      mean_b = mean(y, na.rm = TRUE),
      delta = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
      p_value = p
    )
  })

  res <- do.call(rbind, res)
  res$FDR <- p.adjust(res$p_value, method = "BH")
  res <- res[order(res$FDR, -abs(res$delta)), ]
  rownames(res) <- NULL
  res
}

#' Detect condition-specific rewiring in SF-event regulation
#'
#' Rewiring is quantified as the difference in SF-expression vs event-PSI
#' correlation between two conditions for each directed edge.
#'
#' @param sf_expression matrix SF x cells.
#' @param psi matrix event x cells.
#' @param sf_event_network data.frame with `SF`, `splicing_event`.
#' @param condition condition label per cell.
#' @param contrast length-2 character vector.
#' @param cor_method correlation method passed to `stats::cor`.
#'
#' @return data.frame with edge-level rewiring score.
#' @export
detect_sf_rewiring <- function(
    sf_expression,
    psi,
    sf_event_network,
    condition,
    contrast = NULL,
    cor_method = "spearman"
) {
  sf_expression <- as.matrix(sf_expression)
  psi <- as.matrix(psi)

  common_cells <- Reduce(intersect, list(colnames(sf_expression), colnames(psi), names(condition)))
  if (length(common_cells) < 4) {
    stop("Not enough shared cells across sf_expression, psi, and condition")
  }

  sf_expression <- sf_expression[, common_cells, drop = FALSE]
  psi <- psi[, common_cells, drop = FALSE]
  condition <- as.character(condition[common_cells])

  if (is.null(contrast)) {
    contrast <- unique(condition)
  }
  if (length(contrast) != 2) {
    stop("contrast must contain exactly two conditions")
  }

  idx_a <- which(condition == contrast[1])
  idx_b <- which(condition == contrast[2])

  net <- sf_event_network[, c("SF", "splicing_event"), drop = FALSE]
  net <- unique(net)
  net <- net[net$SF %in% rownames(sf_expression) & net$splicing_event %in% rownames(psi), , drop = FALSE]

  res <- lapply(seq_len(nrow(net)), function(i) {
    sf <- net$SF[i]
    ev <- net$splicing_event[i]
    x_a <- sf_expression[sf, idx_a]
    y_a <- psi[ev, idx_a]
    x_b <- sf_expression[sf, idx_b]
    y_b <- psi[ev, idx_b]

    cor_a <- suppressWarnings(stats::cor(x_a, y_a, method = cor_method, use = "pairwise.complete.obs"))
    cor_b <- suppressWarnings(stats::cor(x_b, y_b, method = cor_method, use = "pairwise.complete.obs"))

    data.frame(
      SF = sf,
      splicing_event = ev,
      cor_a = cor_a,
      cor_b = cor_b,
      rewiring = cor_a - cor_b
    )
  })

  out <- do.call(rbind, res)
  out <- out[order(-abs(out$rewiring)), ]
  rownames(out) <- NULL
  out
}
