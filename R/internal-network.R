.copy_bipartite_attributes <- function(from, to) {
  if (isTRUE(attr(from, "directed_bipartite"))) {
    attr(to, "directed_bipartite") <- TRUE
    attr(to, "left_nodes") <- attr(from, "left_nodes")
    attr(to, "right_nodes") <- attr(from, "right_nodes")
  }
  to
}

.prepare_sf_event_network <- function(Network, node_names) {
  if ((!is.data.frame(Network) && !is.matrix(Network)) || ncol(Network) < 2L)
    stop("Network must contain at least two columns.", call. = FALSE)
  net <- as.data.frame(Network[, 1:2, drop = FALSE], stringsAsFactors = FALSE)
  names(net) <- c("source", "target")
  net$source <- as.character(net$source)
  net$target <- as.character(net$target)
  valid <- !is.na(net$source) & !is.na(net$target) & nzchar(net$source) & nzchar(net$target)
  net <- unique(net[valid, , drop = FALSE])
  if (!nrow(net)) stop("No valid SF-to-event edges remain.", call. = FALSE)
  if (length(setdiff(unique(net$source), node_names)) ||
      length(setdiff(unique(net$target), node_names)))
    stop("All source and target nodes must occur in counts/data row names.", call. = FALSE)
  overlap <- intersect(unique(net$source), unique(net$target))
  if (length(overlap)) stop("SF source and splicing-event target node sets must not overlap: ",
                            paste(overlap, collapse = ", "), call. = FALSE)
  attr(net, "directed_bipartite") <- TRUE
  attr(net, "left_nodes") <- sort(unique(net$source))
  attr(net, "right_nodes") <- sort(unique(net$target))
  rownames(net) <- NULL
  net
}

.validate_scDNS_inputs <- function(counts, data, GroupLabel, network.type) {
  if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
      is.null(rownames(data)) || is.null(colnames(data)))
    stop("counts and data must have row and column names.", call. = FALSE)
  if (length(GroupLabel) != ncol(data))
    stop("GroupLabel length must equal the number of cells.", call. = FALSE)
  if (network.type == "sf_event") {
    if (!identical(rownames(counts), rownames(data)) || !identical(colnames(counts), colnames(data)))
      stop("counts and data must have identical row names, column names, and ordering.", call. = FALSE)
    vals <- as.numeric(data)
    if (any(!is.finite(vals))) stop("sf_event data must contain only finite values.", call. = FALSE)
    if (any(vals < 0 | vals > 1)) stop("sf_event data values must be within 0-1.", call. = FALSE)
  }
  invisible(TRUE)
}
