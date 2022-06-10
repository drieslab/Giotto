python_path = NULL 
if(is.null(python_path)) {
  installGiottoEnvironment()
}

data_path = system.file("extdata/", "visium_brain_data/", package = 'Giotto')

### TESTS FOR 10X VISIUM MOUSE BRAIN DATASET
# --------------------------------------------------------------

# CREATE VISIUM OBJECT
object = createGiottoVisiumObject(visium_dir = data_path,
                                  expr_data = 'raw',
                                  png_name = 'tissue_lowres_image.png',
                                  gene_column_index = 2)
