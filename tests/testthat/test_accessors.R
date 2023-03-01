

if(!require(remotes)){
  install.packages('R.utils', repos = 'http://cran.us.r-project.org')
}

if(!require(remotes)){
  install.packages('remotes', repos = 'http://cran.us.r-project.org')
}

if(!require(GiottoData)){
  library(remotes)
  install_github('drieslab/GiottoData')
}


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






### TESTS EXTERNAL SETTERS ####
## ------------------------------------------------------------------------ ##

#### Native setting ####






test_that('Native exprObj can be set directly', {

  test_ex = expect_no_condition(setExpression(test, ex))

  avail_ex = list_expression(test_ex)
  expect_s3_class(avail_ex, 'data.table')

  expect_identical(avail_ex$spat_unit, 'aggregate')
  expect_identical(avail_ex$feat_type, 'rna')
  expect_identical(avail_ex$name, 'normalized')

})

## ------------------------------------------------------------------------ ##


test_that('Native feature info is set directly', {

  test_fi = expect_no_condition(setFeatureInfo(test, gpoints))

  avail_fi = list_feature_info(test_fi)
  expect_s3_class(avail_fi, 'data.table')

  expect_identical(avail_fi$feat_info, 'rna')
})


test_that('Native feature info is set with lists', { # issues currently happen with unnamed lists

  # assign names by list names
  test_fi = setFeatureInfo(test, gpoints = list(rna = gpoints,
                                                protein = gpoints))
  avail_fi = list_feature_info(test_fi)
  expect_identical(avail_fi$feat_info, c('rna', 'protein'))
  expect_identical(test_fi@feat_info$rna@feat_type, 'rna')
  expect_identical(test_fi@feat_info$protein@feat_type, 'protein')

  # assign names by feat_type tag
  gp_1 = gp_2 = gpoints
  featType(gp_1) = 'new1'
  featType(gp_2) = 'new2'
  test_fi = setFeatureInfo(test, gpoints = list(gp_1, gp_2))

  avail_fi = list_feature_info(test_fi)
  expect_identical(avail_fi$feat_info, c('new1', 'new2'))
  expect_identical(test_fi@feat_info$new1@feat_type, 'new1')
  expect_identical(test_fi@feat_info$new2@feat_type, 'new2') # passes, but the message is wrong

})


## ------------------------------------------------------------------------ ##

test_that('Native spatial info is set directly', {

  test_si = expect_no_error(setPolygonInfo(test, gpoly))

  avail_si = list_spatial_info(test_si)
  expect_s3_class(avail_si, 'data.table')

  expect_identical(avail_si$spat_info, 'aggregate')
})

test_that('Spatlocs is also set if centroids are available', {

  # spat_unit (polygon_name) not explicitly set
  test_si = setPolygonInfo(test, gpoly)

  avail_spatlocs = list_spatial_locations(test_si)
  expect_s3_class(avail_spatlocs, 'data.table')

  expect_identical(avail_spatlocs$spat_unit, 'aggregate')
  expect_identical(avail_spatlocs$name, 'raw')

  # spat_unit (polygon_name) explicitly set
  test_si = setPolygonInfo(test, gpoly, polygon_name = 'test_unit')

  avail_spatlocs = list_spatial_locations(test_si)
  expect_s3_class(avail_spatlocs, 'data.table')

  expect_identical(avail_spatlocs$spat_unit, 'test_unit')
  expect_identical(avail_spatlocs$name, 'raw')
})

## ------------------------------------------------------------------------ ##


## ------------------------------------------------------------------------ ##


## ------------------------------------------------------------------------ ##




## ------------------------------------------------------------------------ ##


test_that('Spatlocs setting requires expression', {

  expect_error(setSpatialLocations(test, sl),
               regexp = 'Add expression')



})

# set expression first
test_sl = setExpression(test, ex)
test_that('Native spatlocs is set', {

  test_sl = suppressWarnings(setSpatialLocations(test_sl, sl))

  avail_sl = list_spatial_locations(test_sl)
  expect_s3_class(avail_sl, 'data.table') # should now exist

  expect_identical(avail_sl$spat_unit, 'aggregate')
  expect_identical(avail_sl$name, 'raw')
})


test_that('Native spatLocsObj is set with user specified nesting', {

  # add needed spat_unit in expression first
  test_sl = setExpression(test_sl, ex, spat_unit = 'new')

  # whether this will work with just spatial_info already tested if
  # test that centroids are set as spatlocs already passed

  test_sl = suppressWarnings(setSpatialLocations(test_sl, sl, spat_unit = 'new'))
  test_sl = suppressWarnings(setSpatialLocations(test_sl, sl, name = 'new'))

  avail_sl = list_spatial_locations(test_sl)
  expect_equal(nrow(avail_sl), 2L)

  set_sl_a = getSpatialLocations(test_sl, spat_unit = 'new')
  set_sl_b = getSpatialLocations(test_sl, name = 'new')

})


test_that('Spatlocs missing spat_unit in expr and spatial_info throws error', {
  # available spat unit in expression is only 'aggregate'
  test_sl = expect_error(setSpatialLocations(test_sl, sl, spat_unit = 'new'),
                         regexp = 'No expression')
})


test_that('Spatlocs spatID mismatch throws error', {
  test_sl = setPolygonInfo(test_sl, gpoly, polygon_name = 'new')
  # in spat_unit 'new', spatIDs have more entries (poly info) than the spatlocs
  # which are based on the later aggregated expression information
  expect_error(setSpatialLocations(test_sl, sl, spat_unit = 'new'),
               regexp = 'between spatial and')
})



# test_that('Native spatLocsObj is set with lists', {
#
# })


test_that('Native spatLocsObj can be removed', {

  test_sl = setSpatialLocations(test_sl, sl)

  test_sl = setSpatialLocations(test_sl,
                                spatlocs = NULL,
                                spat_unit = 'aggregate',
                                name = 'raw')

  expect_null(test_sl@spatial_locs)

})

rm(test_sl)
