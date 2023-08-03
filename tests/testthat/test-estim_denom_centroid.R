test_that("denom_centro", {
  n1 = sample(1:30,1)
  denoms = sapply(1:2, function(i){estim_denom_centroid(density_ratio = rep(1,30), cell_numbers = c(rep(1,n1), rep(2,30-n1)), cell = i)})
  expect_equal(denoms, c(n1,30-n1))
})
