python_path = NULL 
if(is.null(python_path)) {
  installGiottoEnvironment()
}

# best way to store and access visium data? can you use wget in
# github actions and would that slow down the workflow significantly?

# expr_path = system.file("extdata", "x", package = 'Giotto')
# loc_path = system.file("extdata", "x", package = 'Giotto')
# meta_path = system.file("extdata", "x", package = 'Giotto')

### TESTS FOR 10X VISIUM MOUSE BRAIN DATASET
# --------------------------------------------------------------