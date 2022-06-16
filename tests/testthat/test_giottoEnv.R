### TESTS FOR GIOTTO ENVIRONMENT INSTALLATION
# --------------------------------------------------------------

# remove any preexisting giotto environment
if (checkGiottoEnvironment() == T) {
  removeGiottoEnvironment()
}

test_that("No Giotto environment exists", {
  expect_false(checkGiottoEnvironment())
})

# install environment

installGiottoEnvironment()

test_that("Environment now exists", {
  expect_true(checkGiottoEnvironment())
})


