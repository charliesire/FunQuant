test_that("distinct_prototypes_1", {
  expect_equal(distinct_prototypes(list(0,0,2)), FALSE)
})

test_that("distinct_prototypes_2", {
  expect_equal(distinct_prototypes(list(0,1,2)), TRUE)
})
