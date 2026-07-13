test_that("bipartite edge swaps preserve named degrees and direction", {
  toy <- make_toy_spdns()
  net <- spDNS::CreatScDNSobject(toy$counts, toy$data, toy$network, toy$group,
                                 network.type = "sf_event", parallel.sz = 1)@Network[, 1:2]
  attr(net, "directed_bipartite") <- TRUE
  attr(net, "left_nodes") <- toy$sf
  attr(net, "right_nodes") <- toy$event
  set.seed(123)
  random_net <- spDNS:::randomize_directed_bipartite_network(net)
  original_source_degree <- table(net$source)
  random_source_degree <- table(random_net$source)
  expect_equal(random_source_degree[names(original_source_degree)], original_source_degree)
  original_target_degree <- table(net$target)
  random_target_degree <- table(random_net$target)
  expect_equal(random_target_degree[names(original_target_degree)], original_target_degree)
  expect_true(all(random_net$source %in% unique(net$source)))
  expect_true(all(random_net$target %in% unique(net$target)))
  expect_false(any(random_net$source %in% unique(net$target)))
  expect_false(any(random_net$target %in% unique(net$source)))
})
