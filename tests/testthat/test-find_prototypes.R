test_that("kmeans", {
  data = matrix(runif(30), ncol = 15)
  protos_1 = find_prototypes(data = data, starting_proto = list(matrix(c(0.3,0.3)), matrix(c(0.7,0.7))), budget = 20)$prototypes
  protos_2 = kmeans(t(data),centers = rbind(c(0.3,0.3),c(0.7,0.7)),iter.max = 20,algorithm = "Lloyd")$centers
  dimnames(protos_2) = NULL
  expect_equal(t(do.call("cbind", protos_1)), protos_2)
})


test_that("kmeans_is", {
  data = matrix(runif(60), ncol = 30)
  data_2 = cbind(data,data[,16:30])
  protos_1 = find_prototypes(data = data, starting_proto = list(matrix(c(0.3,0.3)), matrix(c(0.7,0.7))), density_ratio = c(rep(2/3, 15), rep(4/3,15)),budget = 20)$prototypes
  protos_2 = kmeans(t(data_2),centers = rbind(c(0.3,0.3),c(0.7,0.7)),iter.max = 20, algorithm = "Lloyd")$centers
  dimnames(protos_2) = NULL
  expect_equal(t(do.call("cbind", protos_1)), protos_2)
})

test_that("errors", {
  data = matrix(runif(60), ncol = 30)
  expect_error(find_prototypes(data = data, starting_proto = list(matrix(c(0.3,0.3)), matrix(c(7,7)))), "One of the initial Voronoi cells is empty")
})



