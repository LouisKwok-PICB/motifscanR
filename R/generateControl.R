generate_control_regions_helper <- function(regions, chrom.size, n.random, genes=NULL,
                                            random.seed=NULL){
  set.seed(random.seed)
  control_regions_df <- data.frame(chromosome = rep(regions@seqnames, each = n.random),
                                   start_old = rep(regions@ranges@start, each = n.random),
                                   width = rep(regions@ranges@width, each = n.random),
                                   strand = rep(regions@strand, each = n.random))
  if (is.null(genes)) {
    control_regions_df$start <- apply(control_regions_df, 1, function(x){
      sample(1:(chrom.size[x['chromosome']] - as.integer(x['width'])), 1)})
  } else{
    control_regions_df$start <- apply(control_regions_df, 1, function(x){
      genes_chrom <- genes[genes@seqnames==x['chromosome']]
      
      if (length(genes_chrom) == 0) {
        return(NA)
      }
      
      distance_tss <- distance2nearestTss(as.integer(x['start_old']), genes_chrom)
      if (is.null(distance_tss)) {
        distance_tss <- sample(10000:100000, 1)
      }
      
      start_temp <- randomStartFromGene(genes_chrom, distance_tss, random.seed)
      while (start_temp < 1| start_temp > (chrom.size[x['chromosome']] - as.integer(x['width']))) {
        start_temp <- startFromGene(genes_chrom, distance_tss)
      }
      start_legal <- start_temp
      return(start_legal)
    })
  }
  control_regions_df <- na.omit(control_regions_df)
  return(GRanges(seqnames = control_regions_df$chromosome,
                 ranges = IRanges(start = control_regions_df$start,
                                  width = control_regions_df$width),
                 strand = control_regions_df$strand))
}

#' generateControl
#'
#' @description function to generate random control region for given regions.
#'
#' @param regions \code{\link[GenomicRanges]{GenomicRanges}} for generating
#' control regions
#' @param genome BSgenome object, or short string signifying genome build
#' recognized by \code{\link[BSgenome]{getBSgenome}}.
#' @param gene.txdb \code{\link[GenomicFeatures]{TxDb}} gene annotation for
#' control regions, See Details.(default=NULL)
#' @param n.random number for generating random control.For each region, generating
#' n.random control regions.(default=5)
#' @param random.seed seed for random programs
#'
#' @return \code{\link[GenomicRanges]{GenomicRanges}} random control regions
#' with n.random times the length of the input regions
#'
#' @export
#'
#' @details If the gene.txdb is specified,  control regions will be
#' generated from the same TSS distance to random gene of  given region.
#' Otherwise, the control region will be generated  by random location
#' with the same chromosome and region width of given regions.
#'
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'
#' example_motifs <- getJasparMotifs(species = "Homo sapiens",
#'                                   collection = "CORE")
#' # Make a set of input regions
#' Input <- GenomicRanges::GRanges(seqnames = c("chr1","chr2","chr2"),
#'                                 ranges = IRanges::IRanges(start = c(76585873,42772928,
#'                                                                     100183786),
#'                                                           width = 500))
#' # Make a set of control regions by generateControl
#' Control <- generateControl(Input, genome=BSgenome.Hsapiens.UCSC.hg19,
#'                            gene.txdb=TxDb.Hsapiens.UCSC.hg19.knownGene)
#'
#' # Scan motif for example motifs
#' motif_ix_input <- motifScan(example_motifs, Input, genome = "BSgenome.Hsapiens.UCSC.hg19")
#' motif_ix_control <- motifScan(example_motifs, Control, genome = "BSgenome.Hsapiens.UCSC.hg19")
#'
#' # Find Enrichment motif of input by control
#' Enrichment_result <- motifEnrichment(motif_ix_input, motif_ix_control)
#'
setGeneric("generateControl",
           function(regions, ...) standardGeneric("generateControl"))

#' @describeIn generateControl GenomicRanges
#' @export
setMethod("generateControl", signature(regions = "GenomicRanges"),
          function(regions,
                   genome,
                   gene.txdb = NULL,
                   n.random = 5,
                   random.seed = NULL) {
            genome <- validate_genome_input(genome)
            chrom.size <- seqlengths(genome)
            genes <- NULL
            if (class(gene.txdb) == "TxDb") {
              genes <- transcripts(gene.txdb)
            }
            if (class(genes) == "GRanges" | is.null(gene.txdb)) {
              generate_control_regions_helper(regions, chrom.size, n.random,
                                              genes, random.seed)
            }else{
              stop('gene.txdb is not supported.')
            }
          })

