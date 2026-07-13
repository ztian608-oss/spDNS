test_that("SF transformation is bounded, stable and row-wise", {
  x <- matrix(c(1, 3, 7, 9, 4, 4, 4, 4, 9, 7, 3, 1), nrow = 3, byrow = TRUE,
              dimnames = list(c("varying", "constant", "balance"), paste0("c", 1:4)))
  y <- spDNS:::transform_sf_expression_01(x)
  expect_identical(dimnames(y), dimnames(x))
  expect_true(all(is.finite(y)))
  expect_true(all(y >= 0 & y <= 1))
  expect_equal(unname(y["constant", ]), rep(0.5, 4))
  expect_identical(order(y["varying", ]), order(x["varying", ]))
})
