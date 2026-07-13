args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
set.seed(123)
n_gene <- 70L
n_cell <- 80L
genes <- paste0("G", seq_len(n_gene))
cells <- paste0("cell", seq_len(n_cell))
counts <- matrix(rpois(n_gene * n_cell, 8) + 1, nrow = n_gene,
                 dimnames = list(genes, cells))
data <- log1p(sweep(counts, 2, colSums(counts), "/") * 10000)
network <- do.call(rbind, lapply(seq_len(n_gene), function(i) data.frame(
  source = genes[i],
  target = genes[((i - 1L) * 3L + 1:6) %% n_gene + 1L],
  stringsAsFactors = FALSE
)))
group <- rep(c("A", "B"), each = n_cell / 2L)
saveRDS(list(counts = counts, data = data, network = network,
             group = group, genes = genes), args[1])
