
require(testthat)

# Load subobjects
ex = GiottoData::loadSubObjectMini('exprObj')
sl = GiottoData::loadSubObjectMini('spatLocsObj')
cm = GiottoData::loadSubObjectMini('cellMetaObj')
fm = GiottoData::loadSubObjectMini('featMetaObj')
sn = GiottoData::loadSubObjectMini('spatialNetworkObj')
enr = GiottoData::loadSubObjectMini('spatEnrObj')
nn = GiottoData::loadSubObjectMini('nnNetObj')
dr = GiottoData::loadSubObjectMini('dimObj')
gpoly = GiottoData::loadSubObjectMini('giottoPolygon')
gpoints = GiottoData::loadSubObjectMini('giottoPoints')



#  TEST SUBOBJECT CREATION ####
## ------------------------------------------------------------------------ ##

# giottoPolygon ####
sv = gpoly[]
DT = spatVector_to_dt(sv)
DF = data.table::setDF(DT)
gp_IDs = spatIDs(gpoly)

test_that('giottoPolygon is created from SpatVector', {
  gp = createGiottoPolygonsFromDfr(sv)
  expect_no_error(validObject(gp))
  expect_s4_class(gp, 'giottoPolygon')
  expect_setequal(gp_IDs, spatIDs(gp))
})

test_that('giottoPolygon is created from data.frame', {
  gp = createGiottoPolygonsFromDfr(DT)
  expect_no_error(validObject(gp))
  expect_s4_class(gp, 'giottoPolygon')
  expect_setequal(gp_IDs, spatIDs(gp))
})

test_that('giottoPolygon is created from data.table', {
  gp = createGiottoPolygonsFromDfr(DF)
  expect_no_error(validObject(gp))
  expect_s4_class(gp, 'giottoPolygon')
  expect_setequal(gp_IDs, spatIDs(gp))
})

# TODO need the file uploaded to do this easily
# test_that('giottoPolygon is created from maskfile', {
#   gp = createGiottoPolygonsFromMask()
# })



# exprObj ####

a = as.array(ex[])
m = as.matrix(ex[])
dgC = ex[]
ex_IDs = spatIDs(ex)

test_that('exprObj is created from array', {
  exprObj = create_expr_obj(exprMat = a)
  expect_no_error(validObject(exprObj))
  expect_s4_class(exprObj, 'exprObj')
  expect_setequal(ex_IDs, spatIDs(exprObj))
})

test_that('exprObj is created from matrix', {
  exprObj = create_expr_obj(exprMat = m)
  expect_no_error(validObject(exprObj))
  expect_s4_class(exprObj, 'exprObj')
  expect_setequal(ex_IDs, spatIDs(exprObj))
})

test_that('exprObj is created from dgCMatrix', {
  exprObj = create_expr_obj(exprMat = dgC)
  expect_no_error(validObject(exprObj))
  expect_s4_class(exprObj, 'exprObj')
  expect_setequal(ex_IDs, spatIDs(exprObj))
})



















# Nearest Network ####





# SETUP
ig = nn[]
nnDT = nnDT_min = data.table::setDT(igraph::as_data_frame(nn[]))
nnDT_min[, c('weight', 'shared', 'rank') := NULL]


## nnNetObj ####

### input evaluation ####

test_that('Eval of nnNetObj returns nnNetObj', {
  expect_identical(evaluate_nearest_networks(nn), nn)
})






## data.table ####

### input evaluation ####

test_that('Eval with missing info throws error', {
  nnDT_min = data.table::copy(nnDT_min)
  nnDT_min[, distance := NULL]
  expect_error(evaluate_nearest_networks(nnDT_min), regexp = 'Unable to coerce')
})

test_that('Eval with minimum info works', {
  ig_new = expect_no_error(evaluate_nearest_networks(nnDT_min))
  expect_s3_class(ig_new, 'igraph')
  expect_true(all(c('distance', 'weight') %in% igraph::list.edge.attributes(ig_new)))
  expect_true('name' %in% igraph::list.vertex.attributes(ig_new))
})






## igraph ####

### input evaluation ####
test_that('Eval of igraph returns igraph', {
  expect_s3_class(evaluate_nearest_networks(ig), 'igraph')
})

test_that('Eval of igraph with no distance attr fails', {
  ig_nodist = igraph::delete_edge_attr(ig, 'distance')
  expect_error(evaluate_nearest_networks(ig_nodist), regexp = 'distance')
})

test_that('Eval of igraph with no name attr fails', {
  ig_noname = igraph::delete_vertex_attr(ig, 'name')
  expect_error(evaluate_nearest_networks(ig_noname), regexp = 'name')
})

test_that('Eval of minimal igraph adds weight attr', {
  ig_min = igraph::delete_edge_attr(ig, 'weight')
  ig_min = evaluate_nearest_networks(ig_min)
  expect_s3_class(ig_min, 'igraph')
  expect_true(all(c('distance', 'weight') %in% igraph::list.edge.attributes(ig_min)))
  expect_true('name' %in% igraph::list.vertex.attributes(ig_min))
  expect_true('weight' %in% igraph::list.edge.attributes(ig_min))
})






## list reading ####


test_that('Read returns nnNetObj list directly', {
  # this warning can't be tested properly
  read_list = suppressWarnings(readNearestNetData(list(nn, nn)))
  expect_true(is.list(read_list))
  expect_true(all(sapply(read_list, class) == 'nnNetObj'))
})


test_that('Depth 1 works', {
  read_list = readNearestNetData(list(ig, ig))
  expect_true(all(sapply(read_list, featType) == 'rna'))
  expect_true(all(sapply(read_list, spatUnit) == 'cell'))
  expect_identical(sapply(read_list, objName), c('nn_1', 'nn_2'))
  expect_identical(sapply(read_list, function(x) x@nn_type), c('nn_1', 'nn_2'))
})


test_that('Depth 2 works', {
  read_list = readNearestNetData(list(test_feat = list(ig,ig),
                                         list(test = ig)))
  expect_identical(sapply(read_list, featType), c('test_feat', 'test_feat', 'feat_2'))
  expect_identical(sapply(read_list, spatUnit), c('cell', 'cell', 'cell'))
  expect_identical(sapply(read_list, objName), c('nn_1', 'nn_2', 'test'))
  expect_identical(sapply(read_list, function(x) x@nn_type), c('nn_1', 'nn_2', 'test'))
})

test_that('Depth 3 works', {
  read_list = readNearestNetData(list(test_unit = list(test_feat = list(a = ig, ig),
                                                          list(ig)),
                                         list(list(b = ig))))
  expect_identical(sapply(read_list, spatUnit), c('test_unit', 'test_unit', 'test_unit', 'unit_2'))
  expect_identical(sapply(read_list, featType), c('test_feat', 'test_feat', 'feat_2', 'feat_1'))
  expect_identical(sapply(read_list, objName), c('a', 'nn_2', 'nn_1', 'b'))
  expect_identical(sapply(read_list, function(x) x@nn_type), c('a', 'nn_2', 'nn_1', 'b'))
})

test_that('Depth 4 works', {
  read_list = readNearestNetData(list(test_unit = list(test_feat = list(list(a = ig),
                                                                           test_meth2 = list(x = ig)),
                                                          list(test_meth = list(ig))),
                                         list(list(list(b = ig)))))
  expect_identical(sapply(read_list, spatUnit), c('test_unit', 'test_unit', 'test_unit', 'unit_2'))
  expect_identical(sapply(read_list, featType), c('test_feat', 'test_feat', 'feat_2', 'feat_1'))
  expect_identical(sapply(read_list, objName), c('a', 'x', 'nn_1', 'b'))
  expect_identical(sapply(read_list, function(x) x@nn_type), c('method_1', 'test_meth2', 'test_meth', 'method_1'))
})










# dimension reduction ####





### list reading ####
drm = dr[]

test_that('Read returns dimObj list directly', {
  read_list = suppressWarnings(readDimReducData(list(dr, dr)))
  expect_true(is.list(read_list))
  expect_true(all(sapply(read_list, class) == 'dimObj'))
})


test_that('Depth 1 works', {
  read_list = readDimReducData(list(drm, drm))
  expect_true(all(sapply(read_list, featType) == 'rna'))
  expect_true(all(sapply(read_list, spatUnit) == 'cell'))
  expect_identical(sapply(read_list, objName), c('dimRed_1', 'dimRed_2'))
  expect_identical(sapply(read_list, function(x) x@reduction_method), c('dimRed_1', 'dimRed_2'))
})


test_that('Depth 2 works', {
  read_list = readDimReducData(list(test_feat = list(drm,drm),
                                         list(test = drm)))
  expect_identical(sapply(read_list, featType), c('test_feat', 'test_feat', 'feat_2'))
  expect_identical(sapply(read_list, spatUnit), c('cell', 'cell', 'cell'))
  expect_identical(sapply(read_list, objName), c('dimRed_1', 'dimRed_2', 'test'))
  expect_identical(sapply(read_list, function(x) x@reduction_method), c('dimRed_1', 'dimRed_2', 'test'))
})

test_that('Depth 3 works', {
  read_list = readDimReducData(list(test_unit = list(test_feat = list(a = drm, drm),
                                                          list(drm)),
                                         list(list(b = drm))))
  expect_identical(sapply(read_list, spatUnit), c('test_unit', 'test_unit', 'test_unit', 'unit_2'))
  expect_identical(sapply(read_list, featType), c('test_feat', 'test_feat', 'feat_2', 'feat_1'))
  expect_identical(sapply(read_list, objName), c('a', 'dimRed_2', 'dimRed_1', 'b'))
  expect_identical(sapply(read_list, function(x) x@reduction_method), c('a', 'dimRed_2', 'dimRed_1', 'b'))
})

test_that('Depth 4 works', {
  read_list = readDimReducData(list(test_unit = list(test_feat = list(list(a = drm),
                                                                           test_meth2 = list(x = drm)),
                                                          list(test_meth = list(drm))),
                                         list(list(list(b = drm)))))
  expect_identical(sapply(read_list, spatUnit), c('test_unit', 'test_unit', 'test_unit', 'unit_2'))
  expect_identical(sapply(read_list, featType), c('test_feat', 'test_feat', 'feat_2', 'feat_1'))
  expect_identical(sapply(read_list, objName), c('a', 'x', 'dimRed_1', 'b'))
  expect_identical(sapply(read_list, function(x) x@reduction_method), c('method_1', 'test_meth2', 'test_meth', 'method_1'))
})




# spatial enrichment ####



## list reading ####

enrDT = enr[]

test_that('Read returns spatEnrObj list directly', {
  read_list = suppressWarnings(readSpatEnrichData(list(enr, enr)))
  expect_true(is.list(read_list))
  expect_true(all(sapply(read_list, class) == 'spatEnrObj'))
})


test_that('Depth 1 works', {
  read_list = readSpatEnrichData(list(enrDT, enrDT))
  expect_true(all(sapply(read_list, featType) == 'rna'))
  expect_true(all(sapply(read_list, spatUnit) == 'cell'))
  expect_identical(sapply(read_list, objName), c('enr_1', 'enr_2'))
  expect_identical(sapply(read_list, function(x) x@method), c('enr_1', 'enr_2'))
})


test_that('Depth 2 works', {
  read_list = readSpatEnrichData(list(test_feat = list(enrDT,enrDT),
                                            list(test = enrDT)))
  expect_identical(sapply(read_list, featType), c('test_feat', 'test_feat', 'feat_2'))
  expect_identical(sapply(read_list, spatUnit), c('cell', 'cell', 'cell'))
  expect_identical(sapply(read_list, objName), c('enr_1', 'enr_2', 'test'))
  expect_identical(sapply(read_list, function(x) x@method), c('enr_1', 'enr_2', 'test'))
})

test_that('Depth 3 works', {
  read_list = readSpatEnrichData(list(test_unit = list(test_feat = list(a = enrDT, enrDT),
                                                             list(enrDT)),
                                            list(list(b = enrDT))))
  expect_identical(sapply(read_list, spatUnit), c('test_unit', 'test_unit', 'test_unit', 'unit_2'))
  expect_identical(sapply(read_list, featType), c('test_feat', 'test_feat', 'feat_2', 'feat_1'))
  expect_identical(sapply(read_list, objName), c('a', 'enr_2', 'enr_1', 'b'))
  expect_identical(sapply(read_list, function(x) x@method), c('a', 'enr_2', 'enr_1', 'b'))
})


test_that('Depth 4 works', {
  read_list = readSpatEnrichData(list(test_unit = list(test_feat = list(list(a = enrDT),
                                                                              test_meth2 = list(x = enrDT)),
                                                             list(test_meth = list(enrDT))),
                                            list(list(list(b = enrDT)))))
  expect_identical(sapply(read_list, spatUnit), c('test_unit', 'test_unit', 'test_unit', 'unit_2'))
  expect_identical(sapply(read_list, featType), c('test_feat', 'test_feat', 'feat_2', 'feat_1'))
  expect_identical(sapply(read_list, objName), c('a', 'x', 'enr_1', 'b'))
  expect_identical(sapply(read_list, function(x) x@method), c('method_1', 'test_meth2', 'test_meth', 'method_1'))
})







# spatial network ####


## list reading ####

snDT = sn[]

test_that('Read returns spatialNetworkObj list directly', {
  read_list = suppressWarnings(readSpatNetData(list(sn, sn)))
  expect_true(is.list(read_list))
  expect_true(all(sapply(read_list, class) == 'spatialNetworkObj'))
})


test_that('Depth 1 works', {
  read_list = readSpatNetData(list(snDT, snDT))
  expect_true(all(sapply(read_list, spatUnit) == 'cell'))
  expect_identical(sapply(read_list, objName), c('sn_1', 'sn_2'))
  expect_identical(sapply(read_list, function(x) x@method), c('sn_1', 'sn_2'))
})


test_that('Depth 2 works', {
  read_list = readSpatNetData(list(test_unit = list(snDT,snDT),
                                         list(test = snDT)))
  expect_identical(sapply(read_list, spatUnit), c('test_unit', 'test_unit', 'unit_2'))
  expect_identical(sapply(read_list, objName), c('sn_1', 'sn_2', 'test'))
  expect_identical(sapply(read_list, function(x) x@method), c('sn_1', 'sn_2', 'test'))
})




# spatial locations ####


## list reading ####

slDT = sl[]

test_that('Read returns spatLocsObj list directly', {
  read_list = suppressWarnings(readSpatLocsData(list(sl, sl)))
  expect_true(is.list(read_list))
  expect_true(all(sapply(read_list, class) == 'spatLocsObj'))
})

test_that('Depth 1 works', {
  read_list = readSpatLocsData(list(slDT, slDT))
  expect_true(all(sapply(read_list, spatUnit) == 'cell'))
  expect_identical(sapply(read_list, objName), c('coord_1', 'coord_2'))
})


test_that('Depth 2 works', {
  read_list = readSpatLocsData(list(test_unit = list(slDT,slDT),
                                         list(test = slDT)))
  expect_identical(sapply(read_list, spatUnit), c('test_unit', 'test_unit', 'unit_2'))
  expect_identical(sapply(read_list, objName), c('coord_1', 'coord_2', 'test'))
})













# expression ####


## list reading ####

exMat = ex[]

test_that('Read returns dimObj list directly', {
  read_list = suppressWarnings(readExprData(list(ex, ex)))
  expect_true(is.list(read_list))
  expect_true(all(sapply(read_list, class) == 'exprObj'))
})


test_that('Depth 1 works', {
  read_list = readExprData(list(exMat, exMat))
  expect_true(all(sapply(read_list, featType) == 'rna'))
  expect_true(all(sapply(read_list, spatUnit) == 'cell'))
  expect_identical(sapply(read_list, objName), c('data_1', 'data_2'))
})


test_that('Depth 2 works', {
  read_list = readExprData(list(test_feat = list(exMat,exMat),
                                            list(test = exMat)))
  expect_identical(sapply(read_list, featType), c('test_feat', 'test_feat', 'feat_2'))
  expect_identical(sapply(read_list, spatUnit), c('cell', 'cell', 'cell'))
  expect_identical(sapply(read_list, objName), c('data_1', 'data_2', 'test'))
})

test_that('Depth 3 works', {
  read_list = readExprData(list(test_unit = list(test_feat = list(a = exMat, exMat),
                                                             list(exMat)),
                                            list(list(b = exMat))))
  expect_identical(sapply(read_list, spatUnit), c('test_unit', 'test_unit', 'test_unit', 'unit_2'))
  expect_identical(sapply(read_list, featType), c('test_feat', 'test_feat', 'feat_2', 'feat_1'))
  expect_identical(sapply(read_list, objName), c('a', 'data_2', 'data_1', 'b'))
})


# cell metadata ####

## list reading ####





# feat metadata ####

## list reading ####









