#' Create epigenetic giotto object
#'
#' @param fragmentsPath A character vector containing the paths to the input files to use to generate the ArrowFiles. These files can be in one of the following formats: (i) scATAC tabix files, (ii) fragment files, or (iii) bam files.
#' @param genome A string indicating the default genome to be used for all ArchR functions. Currently supported values include "hg19","hg38","mm9", and "mm10". This value is stored as a global environment variable, not part of the ArchRProject. This can be overwritten on a per-function basis using the given function's geneAnnotationand genomeAnnotation parameter. For something other than one of the currently supported, see createGeneAnnnotation() and createGenomeAnnnotation()
#' @param expression expression information
#' @param expression_feat available features (e.g. atac, rna, ...)
#' @param spatial_locs data.table or data.frame with coordinates for cell centroids
#' @param verbose Default = TRUE
#' @param sampleNames A character vector containing the names to assign to the samples
#' @param minFrags minFrags 
#' @param maxFrags maxFrags 
#' @param minTSS minTSS
#' @param addTileMat addTileMat 
#' @param addGeneScoreMat addGeneScoreMat
#' @param offsetPlus offsetPlus
#' @param offsetMinus offsetMinus
#' @param TileMatParams TileMatParams
#' @param force TRUE
#' @param outputDirectory outputDirectory 
#' @param copyArrows copyArrows 
#' @param useMatrix useMatrix 
#' @param name name
#' @param ... Additional params passed to `createGiottoObject`
#'
#' @return A Giotto object with at least an atac or epigenetic modality
#' @export
#'
createEpiObject <- function(fragmentsPath,
                            genome = c('hg19', 'hg38', 'mm9', 'mm10'),
                            expression = NULL,
                            expression_feat = 'atac',
                            spatial_locs = NULL,
                            sampleNames = 'sample1',
                            minFrags = 0,
                            maxFrags = 1e+07,
                            minTSS = 0,
                            addTileMat = TRUE,
                            addGeneScoreMat = TRUE,
                            offsetPlus = 0,
                            offsetMinus = 0,
                            TileMatParams = list(tileSize = 5000),
                            force = TRUE,
                            outputDirectory = getwd(),
                            copyArrows = FALSE,
                            useMatrix = "TileMatrix", 
                            name = "IterativeLSI", 
                            verbose = TRUE,
                            ...) {
    
    if(!requireNamespace('ArchR')) {
        wrap_msg('ArchR is needed. Install the package using remotes::install_github("GreenleafLab/ArchR")')
    } else {require(ArchR)}
    
    ## Add reference genome
    wrap_msg('Loading reference genome')
    ArchR::addArchRGenome(genome)
    
    # Creating Arrow Files
    wrap_msg('Creating Arrow files')
    ArrowFiles <- ArchR::createArrowFiles(
        inputFiles = fragmentsPath,
        sampleNames = sampleNames,
        minFrags = minFrags,
        maxFrags = maxFrags,
        minTSS = minTSS,
        addTileMat = addTileMat,
        addGeneScoreMat = addGeneScoreMat,
        offsetPlus = offsetPlus,
        offsetMinus = offsetMinus,
        TileMatParams = TileMatParams,
        force = force
    )
    
    # Creating an ArchRProject 
    wrap_msg('Creating ArchRProject')
    proj <- ArchR::ArchRProject(
        ArrowFiles = ArrowFiles, 
        outputDirectory = getwd(),
        copyArrows = FALSE
    )
    
    # Data normalization and dimensionality reduction 
    wrap_msg('Running dimension reduction')
    proj <- addIterativeLSI(
        ArchRProj = proj,
        useMatrix = "TileMatrix", 
        name = "IterativeLSI"
    )
    
    # extract GeneScoreMatrix
    GeneScoreMatrix_summarizedExperiment <- ArchR::getMatrixFromProject(proj)
    GeneScoreMatrix <- slot(slot(GeneScoreMatrix_summarizedExperiment, 'assays'), 'data')[['GeneScoreMatrix']]
    
    ## get cell names
    cell_names <- colnames(GeneScoreMatrix)
    cell_names <- gsub(paste0(sampleNames,'#'),'',cell_names)
    cell_names <- gsub('-1','',cell_names)
    
    ## get gene names
    gene_names <- slot(GeneScoreMatrix_summarizedExperiment,'elementMetadata')[['name']]
    
    ## replace colnames with cell names
    colnames(GeneScoreMatrix) <- cell_names
    
    ## replace rownames with gene names
    rownames(GeneScoreMatrix) <- gene_names
    
    wrap_msg('Writting local GeneScoreMatrix.csv file')
    data.table::fwrite(as.data.frame(GeneScoreMatrix), 
                       "GeneScoreMatrix.csv",
                       row.names = TRUE)
    
    # filter spatial coordinates
    wrap_msg('Writting filtered coordinates spatial_coords_filtered.csv')
    x <- read.csv(spatial_locs, header = FALSE)
    rownames(x) <- x$V1
    x <- x[cell_names,c(1,3,4)]
    colnames(x) <- c('cell_ID', 'sdimx', 'sdimy')
    write.table(x, 'spatial_coords_filtered.csv',sep = ',', quote = FALSE, row.names = FALSE)
    
    # Creating GiottoObject 
    wrap_msg('Creating GiottoObject')
    
    if(!is.null(expression)) {
        gobject <- createGiottoObject(expression = list("GeneScoreMatrix.csv", expression),
                                      expression_feat = expression_feat,
                                      spatial_locs = 'spatial_coords_filtered.csv',
                                      ...)
    } else {
        gobject <- createGiottoObject(expression = "GeneScoreMatrix.csv",
                                      expression_feat = expression_feat,
                                      spatial_locs = 'spatial_coords_filtered.csv',
                                      ...)
    }
    
    return(gobject)
}
