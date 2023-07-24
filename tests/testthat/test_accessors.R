
require(testthat)

# Load subobjects
ex = GiottoData::loadSubObjectMini('exprObj')
sl = GiottoData::loadSubObjectMini('spatLocsObj')
cm = GiottoData::loadSubObjectMini('cellMetaObj')
fm = GiottoData::loadSubObjectMini('featMetaObj')
sn = GiottoData::loadSubObjectMini('spatialNetworkObj')
enr = GiottoData::loadSubObjectMini('spatEnrObj')
dr = GiottoData::loadSubObjectMini('dimObj')
nn = GiottoData::loadSubObjectMini('nnNetObj')
gpoly = GiottoData::loadSubObjectMini('giottoPolygon')
gpoints = GiottoData::loadSubObjectMini('giottoPoints')


# create test object
test = giotto()




ex1 = ex2 = ex
objName(ex1) = 'data1'
featType(ex1) = 'protein'
objName(ex2) = 'data2'
featType(ex2) = 'feat3'
spatUnit(ex2) = 'nucleus'



# GETTERS ####

## missing cases ####

test_that('Not found exprObj returns error', {
  expect_error(
    getExpression(test, spat_unit = 'none', feat_type = 'none', values = 'raw')
  )
})




# SETTERS ####
## ------------------------------------------------------------------------ ##

#### setting: expression ####

ex1 = ex2 = ex
objName(ex1) = 'data1'
featType(ex1) = 'protein'
objName(ex2) = 'data2'
featType(ex2) = 'feat3'
spatUnit(ex2) = 'nucleus'


test_that('Single: exprObj can be set', {

  test_ex = setExpression(test, ex)

  avail_ex = list_expression(test_ex)
  expect_s3_class(avail_ex, 'data.table')

  expect_identical(avail_ex$spat_unit, 'aggregate')
  expect_identical(avail_ex$feat_type, 'rna')
  expect_identical(avail_ex$name, 'normalized')

})

test_that('List: exprObj can be set', {

  test_ex = setExpression(test, list(ex, ex1, ex2))

  avail_ex = list_expression(test_ex)
  expect_s3_class(avail_ex, 'data.table')

  expect_identical(avail_ex$spat_unit, c('aggregate', 'aggregate', 'nucleus'))
  expect_identical(avail_ex$feat_type, c('rna', 'protein', 'feat3'))
  expect_identical(avail_ex$name, c('normalized', 'data1', 'data2'))

})

test_that('Non-native throws error', {
  test_ex = expect_error(setExpression(test, ex[]), regexp = 'Only exprObj')
})

spatUnit(ex) = 'aggregate'
spatUnit(ex2) = 'nucleus'
test_ex = setExpression(test, ex) # spat_unit aggregate
test_ex2 = setExpression(test_ex, ex2) # spat_unit nucleus


## setting: cell metadata ####




## setting: feat metadata ####




## setting: spatial locations ####

sl1 = sl2 = sl
objName(sl1) = 'data1'
objName(sl2) = 'data2'
spatUnit(sl2) = 'nucleus'


test_that('Single: spatLocsObj can be set', {

  test_ex = setSpatialLocations(test_ex, sl)

  avail_ex = list_spatial_locations(test_ex)
  expect_s3_class(avail_ex, 'data.table')

  expect_identical(avail_ex$spat_unit, spatUnit(sl))
  expect_identical(avail_ex$name, objName(sl))

})

test_that('List: spatLocsObj can be set', {
  # setup
  test_ex = setExpression(test_ex, ex2)

  # tests
  test_ex = setSpatialLocations(test_ex, list(sl, sl1, sl2))

  avail_ex = list_spatial_locations(test_ex)
  expect_s3_class(avail_ex, 'data.table')

  expect_identical(avail_ex$spat_unit, c('aggregate', 'aggregate', 'nucleus'))
  expect_identical(avail_ex$name, c('raw', 'data1', 'data2'))

})

test_that('Non-native throws error', {
  expect_error(setSpatialLocations(test_ex, sl[]), regexp = 'Only spatLocsObj')
})




## setting: spatial networks ####

sn1 = sn2 = sn
objName(sn1) = 'data1'
objName(sn2) = 'data2'
spatUnit(sl) = 'aggregate'
spatUnit(sl2) = 'nucleus'
spatUnit(sn2) = 'nucleus'

test_that('Spatial network requires matching spatial locations', {
  expect_error(setSpatialNetwork(test_ex, sn2), regexp = 'Add spatial location') # none
  test_ex = setSpatialLocations(test_ex, sl)
  expect_error(setSpatialNetwork(test_ex, sn2), regexp = 'Matching') # no match (nucleus vs aggregate)

  # test_ex2 contains info with spat_unit = 'nucleus'
  test_ex2 = setSpatialLocations(test_ex2, sl2)
  expect_no_error(setSpatialNetwork(test_ex2, sn2)) # should now work with correct sl spat_unit
})

test_that('Single: spatialNetworkObj can be set', {
  test_ex = setSpatialLocations(test_ex, sl)
  test_ex = setSpatialNetwork(test_ex, sn)

  avail_ex = list_spatial_networks(test_ex)
  expect_s3_class(avail_ex, 'data.table')

  expect_identical(avail_ex$spat_unit, spatUnit(sn))
  expect_identical(avail_ex$name, objName(sn))

})

test_that('List: spatialNetworkObj can be set', {
  # setup
  test_ex = setSpatialLocations(test_ex, sl)
  test_ex = setExpression(test_ex, ex2)
  test_ex = setSpatialLocations(test_ex, sl2)

  # tests
  test_ex = setSpatialNetwork(test_ex, list(sn, sn1, sn2))

  avail_ex = list_spatial_networks(test_ex)
  expect_s3_class(avail_ex, 'data.table')

  expect_identical(avail_ex$spat_unit, c('aggregate', 'aggregate', 'nucleus'))
  expect_identical(avail_ex$name, c(objName(sn), 'data1', 'data2'))

})

test_that('Non-native throws error', {
  test_ex = setSpatialLocations(test_ex, sl)
  expect_error(setSpatialNetwork(test_ex, sn[]), regexp = 'Only spatialNetworkObj')
})




## setting: spatial enrichment ####

enr1 = enr2 = enr
objName(enr1) = 'data1'
objName(enr2) = 'data2'
spatUnit(sl) = 'aggregate'
spatUnit(sl2) = 'nucleus'
spatUnit(enr2) = 'nucleus'

test_that('Spatial enrichment requires matching spatial locations', {
  expect_error(setSpatialEnrichment(test_ex, enr2), regexp = 'Add spatial location') # none
  test_ex = setSpatialLocations(test_ex, sl)
  expect_error(setSpatialEnrichment(test_ex, enr2), regexp = 'Matching') # no match (nucleus vs aggregate)

  # test_ex2 contains info with spat_unit = 'nucleus'
  test_ex2 = setSpatialLocations(test_ex2, sl2)
  expect_no_error(setSpatialEnrichment(test_ex2, enr2)) # should now work with correct sl spat_unit
})

test_that('Single: spatEnrObj can be set', {
  test_ex = setSpatialLocations(test_ex, sl)
  test_ex = setSpatialEnrichment(test_ex, enr)

  avail_se = list_spatial_enrichments(test_ex)
  expect_s3_class(avail_se, 'data.table')

  expect_identical(avail_se$spat_unit, spatUnit(enr))
  expect_identical(avail_se$name, objName(enr))

})

test_that('List: spatEnrObj can be set', {
  # setup
  test_ex = setSpatialLocations(test_ex, sl)
  test_ex = setExpression(test_ex, ex2)
  test_ex = setSpatialLocations(test_ex, sl2)

  # tests
  test_ex = setSpatialEnrichment(test_ex, list(enr, enr1, enr2))

  avail_se = list_spatial_enrichments(test_ex)
  expect_s3_class(avail_se, 'data.table')

  expect_identical(avail_se$spat_unit, c('aggregate', 'aggregate', 'nucleus'))
  expect_identical(avail_se$name, c(objName(enr), 'data1', 'data2'))

})

test_that('Non-native throws error', {
  test_ex = setSpatialLocations(test_ex, sl)
  expect_error(setSpatialEnrichment(test_ex, enr[]), regexp = 'Only spatEnrObj')
})




## setting: dim reduction ####

dr1 = dr2 = dr
objName(dr1) = 'data1'
featType(dr1) = 'test_feat'
objName(dr2) = 'data2'
spatUnit(dr2) = 'nucleus'
featType(dr2) = 'test_feat'

test_that('Dim red requires matching expression', {
  expect_error(setDimReduction(test, dr2), regexp = 'Add expression') # none
  expect_error(setDimReduction(test_ex, dr2), regexp = 'Matching') # no match (nucleus vs aggregate)

  # test_ex2 contains info with spat_unit = 'nucleus'
  expect_error(setDimReduction(test_ex2, dr2), regexp = 'Matching') # wrong feat
  featType(dr2) = 'feat3'
  expect_no_error(setDimReduction(test_ex2, dr2))
})

test_that('Single: dimObj can be set', {
  test_ex = setSpatialLocations(test_ex, sl)
  test_ex = setDimReduction(test_ex, dr)

  avail_dr = list_dim_reductions(test_ex)
  expect_s3_class(avail_dr, 'data.table')

  expect_identical(avail_dr$spat_unit, spatUnit(dr))
  expect_identical(avail_dr$name, objName(dr))

})

test_that('List: dimObj can be set', {
  # setup
  featType(ex2) = 'test_feat'
  test_ex = setExpression(test_ex, ex2)
  featType(ex2) = 'rna'
  test_ex = setExpression(test_ex, ex2)
  featType(ex) = 'test_feat'
  test_ex = setExpression(test_ex, ex)

  # tests
  test_ex = setDimReduction(test_ex, list(dr, dr1, dr2))

  avail_dr = list_dim_reductions(test_ex)
  expect_s3_class(avail_dr, 'data.table')

  expect_identical(avail_dr$spat_unit, c('aggregate', 'aggregate', 'nucleus'))
  expect_identical(avail_dr$feat_type, c('rna', 'test_feat', 'test_feat'))
  expect_identical(avail_dr$name, c(objName(dr), 'data1', 'data2'))

})

test_that('Non-native throws error', {
  test_ex = setSpatialLocations(test_ex, sl)
  expect_error(setDimReduction(test_ex, dr[]), regexp = 'Only dimObj')
})





## setting: nearest networks ####

nn1 = nn2 = nn
objName(nn1) = 'data1'
featType(nn1) = 'test_feat'
objName(nn2) = 'data2'
spatUnit(nn2) = 'nucleus'
featType(nn2) = 'test_feat'

test_that('Nearest neighbors requires matching dimreduction', {
  expect_error(setNearestNetwork(test, nn2), regexp = 'Add dimension reduction') # none
  test_ex = setDimReduction(test_ex, dr) # no match (nucleus vs aggregate)

  # test_ex2 contains info with spat_unit = 'nucleus'
  expect_error(setNearestNetwork(test_ex, nn2), regexp = 'Matching') # wrong feat
  featType(dr2) = 'test_feat'
  featType(ex2) = 'test_feat'

  test_ex = setExpression(test_ex, ex2)
  test_ex = setDimReduction(test_ex, dr2)
  expect_no_error(setNearestNetwork(test_ex, nn2))
})

# test_that('Single: nnNetObj can be set', {
#   test_ex = setSpatialLocations(test_ex, sl)
#   test_ex = setDimReduction(test_ex, nn)
#
#   avail_nn = list_dim_reductions(test_ex)
#   expect_s3_class(avail_nn, 'data.table')
#
#   expect_identical(avail_nn$spat_unit, spatUnit(nn))
#   expect_identical(avail_nn$name, objName(nn))
#
# })
#
# test_that('List: nnNetObj can be set', {
#   # setup
#   featType(ex2) = 'test_feat'
#   test_ex = setExpression(test_ex, ex2)
#   featType(ex2) = 'rna'
#   test_ex = setExpression(test_ex, ex2)
#   featType(ex) = 'test_feat'
#   test_ex = setExpression(test_ex, ex)
#
#   # tests
#   test_ex = setNearestNetwork(test_ex, list(nn, nn1, nn2))
#
#   avail_nn = list_nearest_networks(test_ex)
#   expect_s3_class(avail_nn, 'data.table')
#
#   expect_identical(avail_nn$spat_unit, c('aggregate', 'aggregate', 'nucleus'))
#   expect_identical(avail_nn$feat_type, c('rna', 'test_feat', 'test_feat'))
#   expect_identical(avail_nn$name, c(objName(nn), 'data1', 'data2'))
#
# })

test_that('Non-native throws error', {
  test_ex = setSpatialLocations(test_ex, sl)
  expect_error(setDimReduction(test_ex, nn[]), regexp = 'Only dimObj')
})











## ------------------------------------------------------------------------ ##



test_that('Native feature info is set directly', {

  test_fi = expect_no_error(setFeatureInfo(test, gpoints))

  avail_fi = list_feature_info(test_fi)
  expect_s3_class(avail_fi, 'data.table')

  expect_identical(avail_fi$feat_info, 'rna')
})


test_that('Native feature info is set with lists', { # issues currently happen with unnamed lists

  # assign names by list names - this now happens through read fxns only
  # test_fi = setFeatureInfo(test, x = list(rna = gpoints,
  #                                         protein = gpoints2))
  # avail_fi = list_feature_info(test_fi)
  # expect_identical(avail_fi$feat_info, c('rna', 'protein'))
  # expect_identical(test_fi@feat_info$rna@feat_type, 'rna')
  # expect_identical(test_fi@feat_info$protein@feat_type, 'protein')

  # assign names by feat_type tag
  gp_1 = gp_2 = gpoints
  featType(gp_1) = 'new1'
  featType(gp_2) = 'new2'
  test_fi = setFeatureInfo(test, x = list(gp_1, gp_2))

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
  test_si = setPolygonInfo(test, gpoly, centroids_to_spatlocs = TRUE)

  avail_spatlocs = list_spatial_locations(test_si)
  expect_s3_class(avail_spatlocs, 'data.table')

  expect_identical(avail_spatlocs$spat_unit, 'aggregate')
  expect_identical(avail_spatlocs$name, 'raw')

  # spat_unit (polygon_name) explicitly set
  test_si = setPolygonInfo(test, gpoly, name = 'test_unit', centroids_to_spatlocs = TRUE)

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
  test_sl = setPolygonInfo(test_sl, gpoly, name = 'new')
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
                                x = NULL,
                                spat_unit = 'aggregate',
                                name = 'raw')

  expect_null(test_sl@spatial_locs)

})

rm(test_sl)
