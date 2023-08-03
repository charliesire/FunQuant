test_that("sort_proto", {
  expect_equal(sort_prototypes(as.list(sample(1:10))), as.list(1:10))
})
