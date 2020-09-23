context("check_scalar/vector/matrix")

test_that("check_scalar errors on missing input", {
  ## check_scalar <- svcommon:::check_scalar
  err <- "object 'a' not found"
  expect_error(check_scalar(a = a), regexp = err)
})

test_that("check_scalar errors on incorrect input", {
  ## check_scalar <- svcommon:::check_scalar
  err <- "'a' must be a numeric scalar."
  expect_error(check_scalar(a = 1:2), regexp = err)
  expect_error(check_scalar(a = t(3)), regexp = err)
  expect_error(check_scalar(a = "a"), regexp = err)
})

test_that("check_scalar works on correct input", {
  ## check_scalar <- svcommon:::check_scalar
  expect_equal(check_scalar(a = a, default = 5), 5)
  expect_equal(check_scalar(a = 4, default = 5), 4)
})

test_that("check_vector errors on missing input", {
  ## check_vector <- svcommon:::check_vector
  err <- "object 'a' not found"
  expect_error(check_vector(a = a), regexp = err)
  expect_error(check_vector(a = a, len = 2), regexp = err)
})

test_that("check_vector errors on incorrect input", {
  ## check_vector <- svcommon:::check_vector
  err <- "'a' must be a numeric vector."
  expect_error(check_vector(a = t(3)), regexp = err)
  expect_error(check_vector(a = "a"), regexp = err)
  err <- "'a' has incorrect length."
  expect_error(check_vector(a = a, len = 2, default = 1:3), regexp = err)
  expect_error(check_vector(a = 1:2, len = 1), regexp = err)
  expect_error(check_vector(a = 1, len = 2), regexp = err)
})

test_that("check_vector works on correct input", {
  ## check_vector <- svcommon:::check_vector
  expect_equal(check_vector(a = a, default = 5), 5)
  expect_equal(check_vector(a = a, len = 2, default = 6), c(6,6))
  expect_equal(check_vector(a = 1:3, len = 3), 1:3)
  expect_equal(check_vector(a = 1:3, default = 5), 1:3)
  expect_equal(check_vector(a = 1, len = 3, promote = TRUE), rep(1, 3))
})

test_that("check_matrix errors on missing input", {
  ## check_matrix <- svcommon:::check_matrix
  err <- "object 'a' not found"
  expect_error(check_matrix(a = a), regexp = err)
  expect_error(check_matrix(a = a, dim = 1:2), regexp = err)
})

test_that("check_matrix errors on incorrect input", {
  ## check_matrix <- svcommon:::check_matrix
  err <- "'a' must be a numeric matrix."
  expect_error(check_matrix(a = 1:2), regexp = err)
  expect_error(check_matrix(a = as.array(1:3)), regexp = err)
  expect_error(check_matrix(a = "a"), regexp = err)
  err <- "'a' has incorrect dimensions."
  expect_error(check_matrix(a = a, dim = 1:2, default = t(1)), regexp = err)
  expect_error(check_matrix(a = t(1:2), dim = 2:1), regexp = err)
})

test_that("check_matrix works on correct input", {
  ## check_matrix <- svcommon:::check_matrix
  expect_equal(check_matrix(a = a, default = t(5)), t(5))
  expect_equal(check_matrix(a = a, dim = 1:2, default = 6), t(c(6,6)))
  expect_equal(check_matrix(a = t(5+1:3), dim = c(1,3)), t(5+1:3))
  expect_equal(check_matrix(a = t(1:3), default = 5), t(1:3))
  expect_equal(check_matrix(a = 1:3, promote = TRUE), as.matrix(1:3))
  expect_equal(check_matrix(a = 2:5, dim = c(4,1), promote = TRUE),
               as.matrix(2:5))
  expect_equal(check_matrix(a = 2:5, dim = c(4,2), promote = TRUE),
               cbind(2:5, 2:5))
})
