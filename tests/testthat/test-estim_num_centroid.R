test_that("num_centro", {
  n1 = sample(1:30,1)
  nums = sapply(1:2, function(i){estim_num_centroid(data = matrix(1, ncol=30), density_ratio = rep(1,30), cell_numbers = c(rep(1,n1), rep(2,30-n1)), cell = i)})
  expect_equal(nums, c(n1,30-n1))
})
