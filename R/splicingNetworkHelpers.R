#' Expand SF-gene network to SF-event network
#'
#' @param sf_gene_network data.frame with columns SF and target_gene
#' @param gene_event_map data.frame with columns gene and splicing_event
#'
#' @return data.frame with columns source and target (SF -> splicing_event)
#' @export
expand_network_to_events <- function(sf_gene_network, gene_event_map) {
  required_sf <- c("SF", "target_gene")
  required_map <- c("gene", "splicing_event")
  if (!all(required_sf %in% colnames(sf_gene_network))) {
    stop("sf_gene_network must contain columns: SF, target_gene")
  }
  if (!all(required_map %in% colnames(gene_event_map))) {
    stop("gene_event_map must contain columns: gene, splicing_event")
  }

  sf_gene_network <- unique(sf_gene_network[, required_sf, drop = FALSE])
  gene_event_map <- unique(gene_event_map[, required_map, drop = FALSE])

  expanded <- merge(
    sf_gene_network,
    gene_event_map,
    by.x = "target_gene",
    by.y = "gene",
    all.x = FALSE,
    all.y = FALSE
  )

  expanded <- unique(expanded[, c("SF", "splicing_event"), drop = FALSE])
  colnames(expanded) <- c("source", "target")
  rownames(expanded) <- NULL

  attr(expanded, "directed_bipartite") <- TRUE
  attr(expanded, "left_nodes") <- sort(unique(expanded$source))
  attr(expanded, "right_nodes") <- sort(unique(expanded$target))
  expanded
}

is_directed_bipartite_network <- function(Network) {
  isTRUE(attr(Network, "directed_bipartite"))
}

rewire_directed_bipartite_once <- function(edge_df, n_swap = NULL) {
  edge_df <- unique(edge_df[, 1:2, drop = FALSE])
  colnames(edge_df)[1:2] <- c("source", "target")

  if (nrow(edge_df) < 2) {
    return(edge_df)
  }

  if (is.null(n_swap)) {
    n_swap <- max(1000, nrow(edge_df) * 10)
  }

  left_nodes <- attr(edge_df, "left_nodes")
  right_nodes <- attr(edge_df, "right_nodes")
  if (is.null(left_nodes)) left_nodes <- unique(edge_df$source)
  if (is.null(right_nodes)) right_nodes <- unique(edge_df$target)

  edge_key <- paste(edge_df$source, edge_df$target, sep = "\t")

  for (i in seq_len(n_swap)) {
    pick <- sample.int(nrow(edge_df), 2)
    e1 <- edge_df[pick[1], , drop = FALSE]
    e2 <- edge_df[pick[2], , drop = FALSE]

    if (e1$source == e2$source || e1$target == e2$target) {
      next
    }

    new1 <- c(e1$source, e2$target)
    new2 <- c(e2$source, e1$target)

    if (!(new1[1] %in% left_nodes && new1[2] %in% right_nodes &&
          new2[1] %in% left_nodes && new2[2] %in% right_nodes)) {
      next
    }

    key1 <- paste(new1[1], new1[2], sep = "\t")
    key2 <- paste(new2[1], new2[2], sep = "\t")
    if (key1 %in% edge_key || key2 %in% edge_key) {
      next
    }

    old1 <- paste(e1$source, e1$target, sep = "\t")
    old2 <- paste(e2$source, e2$target, sep = "\t")

    edge_df[pick[1], c("source", "target")] <- new1
    edge_df[pick[2], c("source", "target")] <- new2

    edge_key[match(old1, edge_key)] <- key1
    edge_key[match(old2, edge_key)] <- key2
  }

  attr(edge_df, "directed_bipartite") <- TRUE
  attr(edge_df, "left_nodes") <- sort(unique(left_nodes))
  attr(edge_df, "right_nodes") <- sort(unique(right_nodes))
  edge_df
}

#' Generate directed degree-preserving random network for SF->event bipartite graph
#'
#' @param Network two-column edge data.frame (source,target), marked by
#'   attr(Network, "directed_bipartite") = TRUE
#' @param n.edge number of random edges to return
#'
#' @return randomized edge data.frame preserving source out-degree and target in-degree
#' @export
randomize_directed_bipartite_network <- function(Network, n.edge = NULL) {
  edge_df <- unique(Network[, 1:2, drop = FALSE])
  colnames(edge_df)[1:2] <- c("source", "target")

  if (is.null(n.edge)) {
    n.edge <- nrow(edge_df)
  }

  if (nrow(edge_df) == 0) {
    return(edge_df)
  }

  n_copy <- ceiling(n.edge / nrow(edge_df))
  random_list <- vector("list", n_copy)
  for (i in seq_len(n_copy)) {
    random_list[[i]] <- rewire_directed_bipartite_once(edge_df)
  }

  random_edges <- do.call(rbind, random_list)
  random_edges <- random_edges[sample.int(nrow(random_edges), n.edge, replace = n.edge > nrow(random_edges)), , drop = FALSE]
  rownames(random_edges) <- NULL

  attr(random_edges, "directed_bipartite") <- TRUE
  attr(random_edges, "left_nodes") <- sort(unique(edge_df$source))
  attr(random_edges, "right_nodes") <- sort(unique(edge_df$target))
  random_edges
}
