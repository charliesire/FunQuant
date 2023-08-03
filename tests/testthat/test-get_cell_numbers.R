test_that("length cell numbers", {
  data = array(1:10^3, dim = c(2,5,2,5,10))
  expect_equal(length(get_cell_numbers(data = data,prototypes = list(array(20,dim = c(2,5,2,5)), array(100,dim = c(2,5,2,5))))), 10)
})

test_that("cell numbers", {
  data = array(1:30)
  expect_equal(get_cell_numbers(data = data,prototypes = list(1.1, 20)), c(rep(1,10), rep(2,20)))
})
