# Directed bipartite null-model helpers added by spDNS. The edge-swap logic
# preserves source out-degree and target in-degree.
is_directed_bipartite_network <- function(Network) {
  isTRUE(attr(Network, "directed_bipartite"))
}

rewire_directed_bipartite_once <- function(edge_df, n_swap = NULL) {
  left_nodes <- attr(edge_df, "left_nodes")
  right_nodes <- attr(edge_df, "right_nodes")
  edge_df <- unique(edge_df[, 1:2, drop = FALSE])
  names(edge_df) <- c("source", "target")
  if (is.null(left_nodes)) left_nodes <- unique(edge_df$source)
  if (is.null(right_nodes)) right_nodes <- unique(edge_df$target)
  if (nrow(edge_df) < 2L) return(.copy_bipartite_attributes(
    structure(edge_df, directed_bipartite = TRUE, left_nodes = left_nodes, right_nodes = right_nodes), edge_df))
  if (is.null(n_swap)) n_swap <- max(1000L, nrow(edge_df) * 10L)
  edge_key <- paste(edge_df$source, edge_df$target, sep = "\t")
  for (i in seq_len(n_swap)) {
    pick <- sample.int(nrow(edge_df), 2L)
    e1 <- edge_df[pick[1L], , drop = FALSE]; e2 <- edge_df[pick[2L], , drop = FALSE]
    if (e1$source == e2$source || e1$target == e2$target) next
    new1 <- c(e1$source, e2$target); new2 <- c(e2$source, e1$target)
    key1 <- paste(new1, collapse = "\t"); key2 <- paste(new2, collapse = "\t")
    if (key1 %in% edge_key || key2 %in% edge_key) next
    edge_df[pick[1L], ] <- new1; edge_df[pick[2L], ] <- new2
    edge_key[pick] <- c(key1, key2)
  }
  attr(edge_df, "directed_bipartite") <- TRUE
  attr(edge_df, "left_nodes") <- sort(unique(left_nodes))
  attr(edge_df, "right_nodes") <- sort(unique(right_nodes))
  edge_df
}

randomize_directed_bipartite_network <- function(Network, n.edge = NULL) {
  left_nodes <- attr(Network, "left_nodes")
  right_nodes <- attr(Network, "right_nodes")
  edge_df <- unique(Network[, 1:2, drop = FALSE])
  names(edge_df) <- c("source", "target")
  attr(edge_df, "left_nodes") <- left_nodes %||% sort(unique(edge_df$source))
  attr(edge_df, "right_nodes") <- right_nodes %||% sort(unique(edge_df$target))
  if (is.null(n.edge)) n.edge <- nrow(edge_df)
  if (!identical(as.integer(n.edge), as.integer(nrow(edge_df)))) {
    stop("A degree-preserving bipartite randomization must retain the observed edge count.", call. = FALSE)
  }
  rewire_directed_bipartite_once(edge_df)
}

.random_bipartite_edge_pool <- function(Network, minimum_edges) {
  n_rep <- max(1L, ceiling(minimum_edges / nrow(Network)))
  out <- do.call(rbind, lapply(seq_len(n_rep), function(i)
    randomize_directed_bipartite_network(Network)))
  rownames(out) <- NULL
  attr(out, "directed_bipartite") <- TRUE
  attr(out, "left_nodes") <- attr(Network, "left_nodes")
  attr(out, "right_nodes") <- attr(Network, "right_nodes")
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x
