motifScan_helper <- function(pwms, seqs, genome, bg, p.cutoff, out, ranges,
                             thread=1, random.seed, cutoff.matrix.loc) {
  if (p.cutoff %ni% c(0.01, 0.001, 1e-04, 1e-05, 1e-06)) {
    stop("p.cutoff should be one of [0.01, 0.001, 1e-04, 1e-05, 1e-06]!")
  }
  convert_pwm_message <- "pwm convertion section"
  message(gettextf("\n################ %s ################", convert_pwm_message))
  message("converting PWM...")
  motif_mats <- convert_pwms(pwms, bg, "log")
  for (pwm_name in names(pwms)) {
    pwms[[pwm_name]]@profileMatrix <- motif_mats[[pwm_name]]
  }
  
  load_cutoff_mat_message <- "motif scan section"
  message(gettextf("\n################ %s ################", load_cutoff_mat_message))
  message("Loading background cutoff matrix...")
  if (out == "matches") {
    tmp_out <- get_motif_ix(pwms, seqs, genome, p.cutoff, thread, random.seed, cutoff.matrix.loc)
    if (is.null(ranges)) {
      out <- SummarizedExperiment(assays =
                                    list(motifScans = as(tmp_out,
                                                         "lMatrix")),
                                  colData =
                                    DataFrame(name = name(pwms),
                                              row.names = names(pwms)))
    } else {
      out <- SummarizedExperiment(assays =
                                    list(motifScans = as(tmp_out,
                                                         "lMatrix")),
                                  rowRanges = ranges,
                                  colData =
                                    DataFrame(name = name(pwms),
                                              row.names = names(pwms)))
    }
  } else if (out == "scores") {
    tmp_out <- get_motif_ix_plus(pwms, seqs, genome, p.cutoff, thread, random.seed, cutoff.matrix.loc)
    tmp_out$motifScans <- as(tmp_out$motifScans, "lMatrix")
    if (is.null(ranges)) {
      out <- SummarizedExperiment(assays = tmp_out,
                                  colData =
                                    DataFrame(name = name(pwms),
                                              row.names = names(pwms)))
    } else {
      out <- SummarizedExperiment(assays = tmp_out,
                                  rowRanges = ranges,
                                  colData =
                                    DataFrame(name = name(pwms),
                                              row.names = names(pwms)))
    }
  } else {
    tmp_out <- get_motif_positions(pwms, seqs, genome, p.cutoff, thread, random.seed, cutoff.matrix.loc)
    cl <- makeCluster(thread)
    if (is.null(ranges)) {
      out <- parallel::parLapply(cl, 1:length(motif_mats), function(x) {
        IRangesList(lapply(seq_along(seqs),
                           function(y){
                             tmp <- IRanges(start = tmp_out$motifLocations[[y,x]] + 1,
                                            width = ncol(motif_mats[[x]]))
                             mcols(tmp) <- DataFrame(strand = tmp_out$motifStrands[[y,x]],
                                                     score = tmp_out$motifScores[[y,x]])
                             tmp
                           }))
      })
      names(out) <- names(pwms)
    } else {
      out <- parallel::parLapply(cl, 1:length(motif_mats), function(x) {
        mx_ix <- unlist(lapply(seq_along(seqs), function(y){rep(y, length(tmp_out$motifLocations[[y,x]]))}))
        GRanges(seqnames(ranges)[mx_ix],
                IRanges(start =
                          start(ranges[mx_ix] + 1) +
                          unlist(tmp_out$motifLocations[, x]),
                        width = ncol(motif_mats[[x]])),
                strand = unlist(tmp_out$motifStrands[,x]),
                score = unlist(tmp_out$motifScores[,x])
        )
      })
      names(out) <- names(pwms)
      out <- GRangesList(out)
    }
    stopCluster(cl)
  }
  message("ALL DONE! Congratulation!")
  return(out)
}



#' motifScan
#'
#' @description function to find motif matches
#'
#' @param pwms either \code{\link[TFBSTools]{PFMatrix}},
#' \code{\link[TFBSTools]{PFMatrixList}}, \code{\link[TFBSTools]{PWMatrix}},
#' \code{\link[TFBSTools]{PWMatrixList}}
#' @param subject either \code{\link[GenomicRanges]{GenomicRanges}},
#' \code{\link[Biostrings]{DNAStringSet}}, \code{\link[Biostrings]{DNAString}},
#'  or character vector
#' @param genome BSgenome object, or \code{\link[Rsamtools]{FaFile}}, or short
#' string signifying genome build recognized by \code{\link[BSgenome]{getBSgenome}}.
#' Only required if subject is \code{\link[GenomicRanges]{GenomicRanges}} or
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}} or if bg is set
#'  to "genome". if the bg is set to "genome" and genome is \code{\link[Rsamtools]{FaFile}},
#'  this function could only work on linux
#' @param bg background nucleotide frequencies. Default is to compute based on
#'  genome, i.e. the specific genome being evaluated. See Details.
#' @param out what to return? see value section
#' @param thread thread for running motifscan
#' @param random.seed seed number for random program
#' @param p.cutoff p-value cutoff for returning motifs, should be one of
#' 0.01, 0.001, 1e-04, 1e-05, 1e-06.(default=1e-04)
#' @param ranges if subject is not GenomicRanges or RangedSummarizedExperiment,
#'  these ranges can be used to specify what ranges the input sequences
#'  correspond to. These ranges will be incorporated into the
#'  SummarizedExperiment output if out is "matches" or "scores" or will be used
#'  to give absolute positions of motifs if out is "positions"
#' @param cutoff.matrix.loc the location of local motif score cutoff file, if the file is
#'  not present, motifscanR will generate by itself and save in the current
#'  working directory as './species_collect_cutoff_motifs_matrix.Rdata'(default),
#'  and the user could specify a specific file directory with 'save_path', 
#'  or specify the cutoff matrix file by user himself with parameter cutoff.matrix.name
#' @param cutoff.matrix.name the name of local motif score cutoff file, if the file is
#'  not present, motifscanR will generate by itself and save in the current
#'  working directory as './species_collect_cutoff_motifs_matrix.Rdata'(default),
#'  and the user could specify a specific file directory with 'save_path',
#'  then save as './[cutoff.matrix.name]_cutoff_motifs_matrix.Rdata, so user
#'  could use the next time if the pwms and the genome are same.
#' @param ... additional arguments depending on inputs
#'
#' @details Background nucleotide frequencies can be set to "genome" for
#'  using the genomice frequencies (in which case a genome must be specified),
#'  "subject" to use the subject sequences or ranges for computing
#'  the nucleotide frequencies, "even" for using 0.25 for each base,
#'  or a numeric vector with A, C, G, and T frequencies.
#'
#' @return Either returns a SummarizedExperiment with a sparse matrix with
#'  values set to TRUE for a match (if out == 'matches'), a
#'  SummarizedExperiment with a matches matrix as well as matrices with the
#'  maximum motif score and total motif counts (if out == 'scores'), or a
#'  \code{\link[GenomicRanges]{GenomicRangesList}} or a list of
#'  \code{\link[IRanges]{IRangesList}} with all the positions of matches
#'  (if out == 'positions')
#'
#' @export
#'
#' @examples
#' example_motifs <- getJasparMotifs(species = "Homo sapiens",
#'                                   collection = "CORE")
#'
#' # Make a set of peaks
#' peaks <- GenomicRanges::GRanges(seqnames = c("chr1","chr2","chr2"),
#'                 ranges = IRanges::IRanges(start = c(76585873,42772928,
#'                                           100183786),
#'                                           width = 500))
#'
#' # Scan motif for example motifs
#' motif_ix <- motifScan(example_motifs, peaks, genome = "BSgenome.Hsapiens.UCSC.hg19")
#'
setGeneric("motifScan",
           function(pwms, subject, ...) standardGeneric("motifScan"))

#' @describeIn motifScan PWMatrixList/DNAStringSet
#' @export
setMethod("motifScan", signature(pwms = "PWMatrixList",
                                 subject = "DNAStringSet"),
          function(pwms,
                   subject,
                   genome = NULL,
                   bg = c("genome","subject","even"),
                   out = c("matches", "scores", "positions"),
                   p.cutoff = 1e-04, ranges = NULL,
                   thread = 1, random.seed = NULL, 
                   cutoff.matrix.loc  = './', cutoff.matrix.name = NULL) {
            message('Checking input arugements...')
            out <- match.arg(out)
            
            if (is.numeric(bg)){
              bg <- check_bg(bg)
            } else{
              bg_method <- match.arg(bg)
              bg <- get_bg(bg_method, subject, genome)
            }
            
            seqs <- as.character(subject)
            if (is.null(cutoff.matrix.name)) {
              cutoff.matrix.loc <- path.connect(cutoff.matrix.loc, get.locfile(pwms, genome))
            }else{
              cutoff.matrix.loc <- path.connect(cutoff.matrix.loc, 
                                                paste0(cutoff.matrix.name, 
                                                       '_cutoff_motifs_matrix.Rdata'))
            }
            motifScan_helper(pwms, seqs, genome, bg, p.cutoff, out, ranges,
                             thread, random.seed, cutoff.matrix.loc)
          })

#' @describeIn motifScan PWMatrixList/character
#' @export
setMethod("motifScan", signature(pwms = "PWMatrixList",
                                 subject = "character"),
          function(pwms,
                   subject,
                   genome = NULL,
                   bg = c("genome","subject","even"),
                   out = c("matches", "scores",  "positions"),
                   p.cutoff = 1e-04, ranges = NULL,
                   thread = 1, random.seed = NULL, 
                   cutoff.matrix.loc  = './', cutoff.matrix.name = NULL) {
            message('Checking input arugements...')
            out <- match.arg(out)
            
            if (is.numeric(bg)){
              bg <- check_bg(bg)
            } else{
              bg_method <- match.arg(bg)
              bg <- get_bg(bg_method, subject, genome)
            }
            if (is.null(cutoff.matrix.name)) {
              cutoff.matrix.loc <- path.connect(cutoff.matrix.loc, get.locfile(pwms, genome))
            }else{
              cutoff.matrix.loc <- path.connect(cutoff.matrix.loc, 
                                                paste0(cutoff.matrix.name, 
                                                       '_cutoff_motifs_matrix.Rdata'))
            }
            motifScan_helper(pwms, subject, genome, bg, p.cutoff, out, ranges,
                             thread, random.seed, cutoff.matrix.loc)
          })

#' @describeIn motifScan PWMatrixList/DNAString
#' @export
setMethod("motifScan", signature(pwms = "PWMatrixList",
                                 subject = "DNAString"),
          function(pwms,
                   subject,
                   genome = NULL,
                   bg = c("genome","subject","even"),
                   out = c("matches", "scores", "positions"),
                   p.cutoff = 1e-04, ranges = NULL,
                   thread = 1, random.seed = NULL, 
                   cutoff.matrix.loc  = './', cutoff.matrix.name = NULL) {
            message('Checking input arugements...')
            out <- match.arg(out)
            
            if (is.numeric(bg)){
              bg <- check_bg(bg)
            } else{
              bg_method <- match.arg(bg)
              bg <- get_bg(bg_method, subject, genome)
            }
            
            seqs <- as.character(subject)
            if (is.null(cutoff.matrix.name)) {
              cutoff.matrix.loc <- path.connect(cutoff.matrix.loc, get.locfile(pwms, genome))
            }else{
              cutoff.matrix.loc <- path.connect(cutoff.matrix.loc, 
                                                paste0(cutoff.matrix.name, 
                                                       '_cutoff_motifs_matrix.Rdata'))
            }
            motifScan_helper(pwms, seqs, genome, bg, p.cutoff, out, ranges,
                             thread, random.seed, cutoff.matrix.loc)
          })

#' @describeIn motifScan PWMatrixList/GenomicRanges
#' @export
setMethod("motifScan", signature(pwms = "PWMatrixList",
                                 subject = "GenomicRanges"),
          function(pwms,
                   subject,
                   genome = GenomeInfoDb::genome(subject),
                   bg = c("genome","subject","even"),
                   out = c("matches", "scores", "positions"),
                   p.cutoff = 1e-04, thread = 1,
                   random.seed = NULL,
                   cutoff.matrix.loc  = './', cutoff.matrix.name = NULL) {
            message('Checking input arugements...')
            out <- match.arg(out)
            GenomicRanges::strand(subject) <- "+"
            subject_seq <- subject
            GenomicRanges::start(subject_seq) <- GenomicRanges::start(subject_seq) + 1
            genome <- validate_genome_input(genome)
            seqs <- BSgenome::getSeq(genome, subject_seq)
            
            if (is.numeric(bg)){
              bg <- check_bg(bg)
            } else{
              bg_method <- match.arg(bg)
              bg <- get_bg(bg_method, seqs, genome)
            }
            
            seqs <- as.character(seqs)
            
            if (is.null(cutoff.matrix.name)) {
              cutoff.matrix.loc <- path.connect(cutoff.matrix.loc, get.locfile(pwms, genome))
            }else{
              cutoff.matrix.loc <- path.connect(cutoff.matrix.loc, 
                                                paste0(cutoff.matrix.name, 
                                                       '_cutoff_motifs_matrix.Rdata'))
            }
            motifScan_helper(pwms, seqs, genome, bg, p.cutoff, out, subject,
                             thread, random.seed, cutoff.matrix.loc)
          })

#' @describeIn motifScan PWMatrixList/RangedSummarizedExperiment
#' @export
setMethod("motifScan", signature(pwms = "PWMatrixList",
                                 subject = "RangedSummarizedExperiment"),
          function(pwms, subject,
                   genome = GenomeInfoDb::genome(subject),
                   bg = c("genome","subject","even"),
                   out = c("matches", "scores", "positions"),
                   p.cutoff = 1e-04, thread = 1,
                   random.seed = NULL, 
                   cutoff.matrix.loc  = './', cutoff.matrix.name = NULL) {
            message('Checking input arugements...')
            out <- match.arg(out)
            motifScan(pwms, rowRanges(subject),
                      genome = genome,
                      bg = bg,
                      out = out,
                      p.cutoff = p.cutoff,
                      thread, random.seed,
                      cutoff.matrix.loc, 
                      cutoff.matrix.name)
          })

#' @describeIn motifScan PWMatrixList/BSGenomeViews
#' @export
setMethod("motifScan", signature(pwms = "PWMatrixList",
                                 subject = "BSgenomeViews"),
          function(pwms,
                   subject,
                   bg = c("genome","subject","even"),
                   out = c("matches", "scores", "positions"),
                   p.cutoff = 1e-04, thread = 1,
                   random.seed = NULL, 
                   cutoff.matrix.loc  = './', cutoff.matrix.name = NULL) {
            message('Checking input arugements...')
            out <- match.arg(out)
            seqs <- as.character(subject)
            ranges <- BSgenome::granges(subject)
            
            if (is.numeric(bg)){
              bg <- check_bg(bg)
            } else{
              bg_method <- match.arg(bg)
              bg <- get_bg(bg_method, subject, BSgenome::subject(subject))
            }
            
            seqs <- as.character(subject)
            
            if (is.null(cutoff.matrix.name)) {
              cutoff.matrix.loc <- path.connect(cutoff.matrix.loc, get.locfile(pwms, genome))
            }else{
              cutoff.matrix.loc <- path.connect(cutoff.matrix.loc, 
                                                paste0(cutoff.matrix.name, 
                                                       '_cutoff_motifs_matrix.Rdata'))
            }
            motifScan_helper(pwms, seqs, genome, bg,
                             p.cutoff, out, ranges,
                             thread, random.seed, cutoff.matrix.loc)
          })

### PFMatrixList ---------------------------------------------------------------


#' @describeIn motifScan PFMatrixList/ANY
#' @export
setMethod("motifScan", signature(pwms = "PFMatrixList", subject = "ANY"),
          function(pwms,
                   subject,
                   ...) {
            
            pwms_list <- do.call(PWMatrixList, lapply(pwms, toPWM))
            motifScan(pwms_list,
                      subject,
                      ...)
          })

# Single PWM input -------------------------------------------------------------

#' @describeIn motifScan PWMatrix/ANY
#' @export
setMethod("motifScan", signature(pwms = "PWMatrix", subject = "ANY"),
          function(pwms,
                   subject,
                   ...) {
            
            pwms_list <- PWMatrixList(pwms)
            motifScan(pwms_list,
                      subject,
                      ...)
          })


# Single PFM -------------------------------------------------------------------

#' @describeIn motifScan PFMatrix/ANY
#' @export
setMethod("motifScan",
          signature(pwms = "PFMatrix", subject = "ANY"),
          function(pwms,
                   subject,
                   ...) {
            
            pwms_list <- PWMatrixList(toPWM(pwms, pseudocounts=0.001))
            motifScan(pwms_list,
                      subject,
                      ...)
          })
