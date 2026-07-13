args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
old <- readRDS(args[1]); new <- readRDS(args[2])

canonical_network <- function(x) {
  x <- as.data.frame(x)
  x <- x[order(as.character(x[[1]]), as.character(x[[2]])), , drop = FALSE]
  rownames(x) <- NULL
  attributes(x)[c("directed_bipartite", "left_nodes", "right_nodes")] <- NULL
  x
}
canonical_zscore <- function(x) {
  x <- as.data.frame(x)
  if ("node_type" %in% names(x)) x$node_type <- NULL
  key <- if ("Gene" %in% names(x)) as.character(x$Gene) else rownames(x)
  x <- x[order(key), , drop = FALSE]
  rownames(x) <- NULL
  x
}
canonical_matrix <- function(x) {
  x <- as.matrix(x)
  x[order(rownames(x)), order(colnames(x)), drop = FALSE]
}
assert_equal <- function(a, b, label) {
  result <- all.equal(a, b, tolerance = 1e-8, check.attributes = TRUE)
  if (!isTRUE(result)) stop(label, " differs: ", paste(result, collapse = "; "))
  message(label, ": equal within tolerance 1e-8")
}

assert_equal(canonical_network(old$Network), canonical_network(new$Network), "Network")
assert_equal(canonical_zscore(old$Zscore), canonical_zscore(new$Zscore), "Zscore")
assert_equal(canonical_matrix(old$scZscore), canonical_matrix(new$scZscore), "scZscore")
