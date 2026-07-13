test_that("the five-step sf_event workflow returns typed SF-focused outputs", {
  skip_on_cran()
  toy <- make_integration_spdns()
  set.seed(123)
  obj <- CreatScDNSobject(toy$counts, toy$data, toy$network, toy$group,
                          network.type = "sf_event", k = 5, n.grid = 12,
                          n.coarse = 6, parallel.sz = 1, n.randNet = 5000)
  obj <- scDNS_1_CalDivs(obj)
  expect_true(isTRUE(attr(obj@Network, "directed_bipartite")))
  obj <- scDNS_2_creatNEAModel_v2(obj, n.randNet = 5000, repTime = 1)
  obj <- scDNS_3_GeneZscore_v2(obj)
  expect_gt(nrow(obj@Zscore), 0)
  expect_true("node_type" %in% names(obj@Zscore))
  expect_setequal(unique(obj@Zscore$node_type), c("SF", "splicing_event"))
  obj <- scDNS_4_scContribution(obj, q.th = 1)
  expect_identical(rownames(obj@scZscore), unique(obj@Network$source))
  expect_identical(colnames(obj@scZscore), colnames(toy$data))
})
