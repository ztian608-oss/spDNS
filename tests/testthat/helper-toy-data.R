make_toy_spdns <- function(seed = 123L, n_cell = 60L) {
  set.seed(seed)
  sf <- paste0("SF", 1:4)
  event <- paste0("EV", 1:12)
  cells <- paste0("cell", seq_len(n_cell))
  group <- rep(c("A", "B"), each = n_cell / 2)
  sf_counts <- matrix(rpois(length(sf) * n_cell, 8), nrow = length(sf),
                      dimnames = list(sf, cells))
  event_support <- matrix(rpois(length(event) * n_cell, 5) + 1,
                          nrow = length(event), dimnames = list(event, cells))
  sf_data <- spDNS:::transform_sf_expression_01(sf_counts)
  event_psi <- matrix(stats::runif(length(event) * n_cell), nrow = length(event),
                      dimnames = list(event, cells))
  net <- data.frame(
    source = rep(sf, each = 4),
    target = unlist(lapply(seq_along(sf), function(i)
      event[((i - 1) * 3 + 0:3) %% length(event) + 1])),
    stringsAsFactors = FALSE
  )
  list(counts = rbind(sf_counts, event_support), data = rbind(sf_data, event_psi),
       network = net, group = group, sf = sf, event = event, cells = cells)
}
