

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



### TEST SUBOBJECT CREATION ####
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











