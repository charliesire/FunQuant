test_that("compute_centro", {
  data = matrix(runif(30), ncol = 15)
  prototypes = list(matrix(c(0.3,0.3)), matrix(c(0.7,0.7)))
  cell_numbers = get_cell_numbers(data,prototypes)
  centro = compute_centroids_and_proba(data,cell_numbers = cell_numbers, density_ratio = rep(1,15))$centroids
  centro_2 = lapply(1:2, function(i){apply(data[,cell_numbers==i], 1, mean)})
  expect_equal(unlist(centro), unlist(centro_2))
})
