context("sun.rise.set()")
library(LakeMetabolizer)

testthat::test_that("Return object class and format", {
  
  # Run sun.rise.set() with test data
  dates <- as.POSIXlt(
    c("2020-08-22", "2020-08-23", "2020-08-24"),
    tz = "America/Chicago")
  r <- sun.rise.set(dates, lat = 41.8781)
  
  # Returned object is a data frame
  expect_equal(class(r), "data.frame")
  
  # Index order is preserved from original implmentation
  expect_equal("sunrise", colnames(r)[1])
  expect_equal("sunset", colnames(r)[2])
  
  # Input time zone attribute is returned
  expect_true(attr(r$sunrise, "tzone") == "America/Chicago")
  
})
