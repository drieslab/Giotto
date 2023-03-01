

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


# nnNetObj ####

## input evaluation ####

test_that('Eval of nnNetObj returns nnNetObj', {
  expect_identical(evaluate_nearest_networks(nn), nn)
})






# data.table ####

## input evaluation ####

test_that('Eval with missing info throws error', {
  nnDT_min = data.table::copy(nnDT_min)
  nnDT_min[, distance := NULL]
  expect_error(evaluate_nearest_networks(nnDT_min), regexp = 'Unable to coerce')
})

test_that('Eval with minimum info works', {
  ig_new = expect_no_error(evaluate_nearest_networks(nn_network = nnDT_min))
  expect_s3_class(ig_new, 'igraph')
  expect_true(all(c('distance', 'weight') %in% igraph::list.edge.attributes(ig_new)))
  expect_true('name' %in% igraph::list.vertex.attributes(ig_new))
})






# igraph ####

## input evaluation ####
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






