

require(testthat)


# Load subobjects
ex = GiottoData::loadSubObjectMini('exprObj')
sl = GiottoData::loadSubObjectMini('spatLocsObj')
cm = GiottoData::loadSubObjectMini('cellMetaObj')
fm = GiottoData::loadSubObjectMini('featMetaObj')
sn = GiottoData::loadSubObjectMini('spatialNetworkObj')
enr = GiottoData::loadSubObjectMini('spatEnrObj')
nn = GiottoData::loadSubObjectMini('nnNetObj')
gpoly = GiottoData::loadSubObjectMini('giottoPolygon')
gpoints = GiottoData::loadSubObjectMini('giottoPoints')





### TESTS FOR GOBJECT FUNCTIONALITY ####
## ------------------------------------------------------------------------ ##

test = giotto()

test_that('Gobject can be generated', {
  expect_s4_class(test, 'giotto')
  expect_true(validObject(test))
})

test_that('giottoInstructions is automatically generated', {
  expect_s3_class(instructions(test), 'giottoInstructions')
})


#### Default Spat/Feat Setting ####

test_that('Error is thrown when no data or actives set', {
  expect_warning(set_default_spat_unit(test), regexp = 'No default for spat_unit could be set')
  expect_warning(set_default_feat_type(test), regexp = 'No default for feat_type could be set')
})

test_that('No actives return NULL when checked', {
  expect_null(activeSpatUnit(test))
  expect_null(activeFeatType(test))
})

test_that('Actives can be set', {
  activeSpatUnit(test) = 'aggregate'
  activeFeatType(test) = 'rna'

  expect_no_warning(set_default_spat_unit(test)) # no warning despite lack of expr or spat_info
  expect_no_warning(set_default_feat_type(test)) # no warning despite lack of expr or spat_info

  expect_identical(set_default_spat_unit(test), 'aggregate')
  expect_identical(set_default_feat_type(test), 'rna')
})

test_that('Providing spat_unit returns unmodified', {
  expect_identical('test_value', set_default_spat_unit(spat_unit = 'test_value'))
})

test_that('Providing feat_type returns unmodified', {
  expect_identical('test_value', set_default_feat_type(feat_type = 'test_value'))
})



#### Aggregate Initialization ####

test_that('Expression initiates ID slots', {

  test_ex = setExpression(test, ex)

  expect_identical(spatIDs(test_ex), spatIDs(ex))
  expect_identical(featIDs(test_ex), featIDs(ex))

  expect_true(inherits(test_ex@cell_ID, 'list'))
  expect_true(inherits(test_ex@cell_ID$aggregate, 'character'))
  expect_true(inherits(test_ex@feat_ID, 'list'))
  expect_true(inherits(test_ex@feat_ID$rna, 'character'))
})



test_that('Expression initiates metadata slots', {

  test_ex = setExpression(test, ex)

  expect_identical(spatIDs(test_ex), pDataDT(test_ex)$cell_ID)
  expect_identical(featIDs(test_ex), fDataDT(test_ex)$feat_ID)

  expect_true(inherits(test_ex@cell_metadata$aggregate$rna, 'cellMetaObj'))
  expect_true(inherits(test_ex@feat_metadata$aggregate$rna, 'featMetaObj'))

  expect_true(inherits(test_ex@cell_metadata$aggregate$rna[], 'data.table'))
  expect_true(inherits(test_ex@feat_metadata$aggregate$rna[], 'data.table'))
})


test_that('Expression sets active spat_unit and feat_type', {
  test_ex = setExpression(test, ex)

  # check in instructions settings
  expect_identical(activeSpatUnit(test_ex), spatUnit(ex))
  expect_identical(activeFeatType(test_ex), featType(ex))

  # check output from default setting
  expect_identical(set_default_spat_unit(test_ex), spatUnit(ex))
  expect_identical(set_default_feat_type(test_ex), featType(ex))
})


test_that('expression_feats slot is set by expression', {

  # test single
  test_ex = setExpression(test, ex, feat_type = 'test_feat')
  expect_identical(test_ex@expression_feat, 'test_feat')
  # test multiple
  test_ex = setExpression(test_ex, ex) # default 'rna'
  expect_identical(test_ex@expression_feat, c('rna', 'test_feat'))
})




#### Subcellular Initialization ####


test_that('Spatial info initiates spat_ID slot', {

  test_si = setPolygonInfo(test, gpoly)

  expect_identical(spatIDs(test_si), spatIDs(gpoly))

  expect_true(inherits(test_si@cell_ID, 'list'))
  expect_true(inherits(test_si@cell_ID$aggregate, 'character'))

})

test_that('Spatial info sets active spat_unit', {

  test_si = setPolygonInfo(test, gpoly)

  expect_identical(activeSpatUnit(test_si), 'aggregate')
  expect_identical(set_default_spat_unit(test_si), 'aggregate')

})


test_that('Feature info initiates feat_ID slot', {

  featType(gpoints) = 'test_feat'
  test_fi = setFeatureInfo(test, gpoints)

  expect_identical(activeFeatType(test_fi), 'test_feat')
  expect_identical(set_default_feat_type(test_fi), 'test_feat')

})


test_that('Spat and Feat info initiates cell_metadata slot', {

  test_sf = setFeatureInfo(test, gpoints)

  expect_null(list_cell_metadata(test_sf))
  expect_null(list_feat_metadata(test_sf))

  test_sf = setPolygonInfo(test_sf, gpoly)

  expect_s3_class(list_cell_metadata(test_sf), 'data.table')
  expect_s3_class(list_feat_metadata(test_sf), 'data.table')

})


test_that('expression_feats slot is set by feature_info', {

  # test single
  test_fi = setFeatureInfo(test, gpoints, feat_type = 'test_feat')
  expect_identical(test_fi@expression_feat, 'test_feat')
  # test multiple
  test_fi = setFeatureInfo(test_fi, gpoints) # default 'rna'
  expect_identical(test_fi@expression_feat, c('rna', 'test_feat'))
})




#### ID interaction interactions ####

test_that('cell_ID from spatial_info is overwritten by expression', {

  expected_IDs = spatIDs(ex)

  test_int = setPolygonInfo(test, gpoly)
  expect_false(identical(spatIDs(test_int), expected_IDs))

  test_int = setExpression(test, ex)
  expect_identical(spatIDs(test_int), expected_IDs)

})


test_that('feat_ID from feat_info is overwritten by expression', {

  expected_IDs = featIDs(ex)

  test_int = setFeatureInfo(test, gpoints)
  test_IDs = expect_warning(featIDs(test_int)) # no spat_unit info
  expect_false(identical(test_IDs, expected_IDs))

  test_int = setExpression(test, ex)
  expect_identical(featIDs(test_int), expected_IDs)

})



### GOBJECT ASSEMBLY: VIZGEN ####

#### Step-wise ####

# test_that('Vizgen mini can be assembled with no errors', {
#   a = giotto()
# })
















