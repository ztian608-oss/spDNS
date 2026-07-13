test_that("sf_event object validates matrices and topology", {
  toy <- make_toy_spdns()
  obj <- CreatScDNSobject(toy$counts, toy$data, toy$network, toy$group,
                          network.type = "sf_event", parallel.sz = 1)
  expect_identical(obj@Other$network_type, "sf_event")
  expect_identical(obj@Other$source_nodes, unique(obj@Network$source))
  expect_true(isTRUE(attr(obj@Network, "directed_bipartite")))
  expect_setequal(attr(obj@Network, "left_nodes"), toy$sf)
  expect_setequal(attr(obj@Network, "right_nodes"), toy$event)
  expect_error(CreatScDNSobject(toy$counts, toy$data[, rev(toy$cells)],
                                toy$network, toy$group, network.type = "sf_event"),
               "identical")
  bad <- toy$data; bad[1, 1] <- 2
  expect_error(CreatScDNSobject(toy$counts, bad, toy$network, toy$group,
                                network.type = "sf_event"), "0-1")
  overlap <- toy$network; overlap$target[1] <- overlap$source[1]
  expect_error(CreatScDNSobject(toy$counts, toy$data, overlap, toy$group,
                                network.type = "sf_event"), "must not overlap")
})

test_that("gene_gene remains the default mode", {
  set.seed(123)
  counts <- matrix(rpois(6 * 20, 5), 6, dimnames = list(paste0("G", 1:6), paste0("C", 1:20)))
  data <- log1p(counts)
  net <- data.frame(from = paste0("G", 1:5), to = paste0("G", 2:6))
  obj <- CreatScDNSobject(counts, data, net, rep(c("A", "B"), each = 10), parallel.sz = 1)
  expect_identical(obj@Other$network_type, "gene_gene")
  expect_false(isTRUE(attr(obj@Network, "directed_bipartite")))
})
