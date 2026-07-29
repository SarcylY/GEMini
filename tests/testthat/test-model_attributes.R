test_that("converting model indices to string", {
  expect_equal(model_numtostring(c(1, 4, 7)), "AAA")
  expect_equal(model_numtostring(c(2, 4, 7)), "BAA")
  expect_equal(model_numtostring(c(2, 6, 7)), "BCA")
  expect_equal(model_numtostring(c(3, 6, 9)), "CCC")
  expect_equal(model_numtostring(c(3, 4, 9, 10, 15, 17, 21, 24, 27)), "CACACBCCC")
})

test_that("model index replacement", {
  expect_equal(replace_ind(c(1, 4, 7), 1), c(1, 4, 7))
  expect_equal(replace_ind(c(1, 4, 7), 2), c(2, 4, 7))
  expect_equal(replace_ind(c(1, 4, 7), 5), c(1, 5, 7))
  expect_equal(replace_ind(c(1, 4, 7), 8), c(1, 4, 8))
  expect_equal(replace_ind(c(2, 5, 9), 3), c(3, 5, 9))
  expect_equal(replace_ind(c(3, 4, 9, 10, 15, 17, 21, 24, 27), 16), c(3, 4, 9, 10, 15, 16, 21, 24, 27))
})
