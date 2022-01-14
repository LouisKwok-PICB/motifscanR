# Not exported

# Not-in operator --------------------------------------------------------------

"%ni%" <- Negate("%in%")

## getNucFreqs function to find background sequence --------------------------


setGeneric("getNucFreqs",
           function(subject, ...) standardGeneric("getNucFreqs"))

setMethod("getNucFreqs", signature(subject = "BSgenome"),
          function(subject) {
            param <- new("BSParams", X = subject, FUN = letterFrequency)
            nucFreqs <- colSums(do.call(rbind,
                                        bsapply(param,
                                                letters =
                                                  c("A", "C", "G", "T"))))
            nucFreqs <- nucFreqs/sum(nucFreqs)
            return(nucFreqs)
          })

#' @importFrom Rsamtools scanFaIndex scanFa
setMethod("getNucFreqs", signature(subject = "FaFile"),
          function(subject) {
            nucFreqs <- c("A" = 0, "C" = 0, "G" = 0, "T" = 0)
            chroms <- scanFaIndex(subject)
            for (i in seq_along(chroms)){
              nucFreqs <- nucFreqs +
                letterFrequency(scanFa(subject, chroms[i])[[1]],
                                letters = c("A", "C", "G", "T"))
            }
            nucFreqs <- nucFreqs/sum(nucFreqs)
            return(nucFreqs)
          })

setMethod("getNucFreqs", signature(subject = "DNAStringSet"),
          function(subject) {
            nucFreqs <- colSums(letterFrequency(subject,
                                                c("A", "C", "G", "T")))
            nucFreqs <- nucFreqs/sum(nucFreqs)
            return(nucFreqs)
          })

setMethod("getNucFreqs", signature(subject = "DNAString"),
          function(subject) {
            nucFreqs <- letterFrequency(subject, c("A", "C", "G", "T"))
            nucFreqs <- nucFreqs/sum(nucFreqs)
            return(nucFreqs)
          })

setMethod("getNucFreqs", signature(subject = "BSgenomeViews"),
          function(subject) {
            nucFreqs <- letterFrequency(subject, c("A", "C", "G", "T"))
            nucFreqs <- nucFreqs/sum(nucFreqs)
            return(nucFreqs)
          })

setMethod("getNucFreqs", signature(subject = "character"),
          function(subject) {
            if (length(subject) == 1) {
              return(getNucFreqs(DNAString(subject)))
            } else {
              return(getNucFreqs(DNAStringSet(subject)))
            }
          })


## validate genome input

setGeneric("validate_genome_input",
           function(genome) standardGeneric("validate_genome_input"))

setMethod("validate_genome_input", signature(genome = "FaFile"),
          function(genome) {
            return(genome)
          })

setMethod("validate_genome_input", signature(genome = "BSgenome"),
          function(genome) {
            return(genome)
          })

setMethod("validate_genome_input", signature(genome = "character"),
          function(genome) {
            return(getBSgenome(genome))
          })

setMethod("validate_genome_input", signature(genome = "character"),
          function(genome) {
            if (any(is.na(genome)))
              stop("No genome provided")
            if (length(genome) > 1){
              stopifnot(all(genome == genome[[1]]))
              genome <- genome[[1]]
            }
            return(getBSgenome(genome))
          })

setMethod("validate_genome_input", signature(genome = "ANY"),
          function(genome) {
            stop("genome input must be a BSgenome, FaFile ",
                 "or a string recognized by getBSgenome")
          })


## toPWM by log wight matrix -----------------------------------------------
setGeneric("toPWM",
           function(x, ...) standardGeneric("toPWM"))

setMethod("toPWM", "PFMatrix",
          function(x, pseudocounts=0.001, bg=NULL){
            if(is.null(bg))
              bg <- bg(x)
            pwmMatrix <- toPWM(as.matrix(x),
                               pseudocounts=pseudocounts,
                               bg=bg)
            pwm <- PWMatrix(ID=ID(x), name=name(x),
                            matrixClass=matrixClass(x),
                            strand=strand(x), bg=bg,
                            tags=tags(x), profileMatrix=pwmMatrix,
                            pseudocounts=pseudocounts)
            pwm
          }
)

setMethod("toPWM", "PFMatrixList",
          function(x, pseudocounts=0.001,
                   bg=NULL){
            ans <- lapply(x, toPWM, pseudocounts=pseudocounts,
                          bg=bg)
            ans <- do.call(PWMatrixList, ans)
            return(ans)
          })

setMethod("toPWM", "matrix",
          function(x, pseudocounts=0.001,
                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25)){
            x <- apply(x, 2, function(x){x/sum(x)})
            x[x==0] <- pseudocounts
            x <- apply(x, 2, function(x){x/sum(x)})
            ans <- log(apply(x, 2, function(x){x/bg}))
            return(ans)
          }
)


## convert_pwm to adjust background in pwm model -------------------------------

convert_pwms <- function(pwms, bg_freqs, out_type) {
  stopifnot(inherits(pwms, "PWMatrixList"))
  lapply(pwms, convert_pwm, bg_freqs, out_type)
}

convert_pwm <- function(pwm, bg_freqs, out_type) {
  type <- pwmType(pwm)
  out <- as.matrix(pwm)
  if (type == "prob") {
    norm_mat <- matrix(bg_freqs, nrow = 4, ncol = length(pwm),
                       byrow = FALSE)
    if (out_type == "log") {
      out <- log(as.matrix(pwm)/norm_mat)
    }else{
      out <- log2(as.matrix(pwm)/norm_mat)
    }
  } else if (type == "log2") {
    norm_mat <- matrix(log2(bg_freqs) - log2(bg(pwm)), nrow = 4,
                       ncol = length(pwm),
                       byrow = FALSE)
    out <- as.matrix(pwm) - norm_mat
    if (out_type == "log") {
      out <- log(2 ** out)
    }
  } else if (type == "log") {
    norm_mat <- matrix(log(bg_freqs) - log(bg(pwm)), nrow = 4,
                       ncol = length(pwm),
                       byrow = FALSE)
    out <- as.matrix(pwm) - norm_mat
    if (out_type == "log2") {
      out <- log2(exp(out))
    }
  }
  return(out)
}

# make sure background is correct ----------------------------------------------
check_bg <- function(bg_freqs){
  if (length(bg_freqs) != 4)
    stop("Invalid background frequencies -- should be length 4")
  if (!all.equal(sum(bg_freqs),1) || min(bg_freqs) <= 0){
    stop("Invalid background frequencies. Should sum to 1")
  }

  if (!is.null(names(bg_freqs))){
    if (!all(names(c("A","C","G","T") %in% bg_freqs))){
      stop("Background nucleotide frequencies have names that ",
           "don't match nucleotides! (A,C,G,T)")
    } else{
      bg_freqs <- bg_freqs[c("A","C","G","T")]
    }
  }
  return(bg_freqs)
}

get.locfile <- function(pwms, genome){
  tags <- pwms[[1]]@tags
  if (inherits(genome, "BSgenome")) {
    return(paste0(c('BSgenome',
                    genome@metadata$genome,
                    gsub(' ', '_', tags$tax_group),
                    tags$collection,
                    'cutoff_motifs_matrix.Rdata'), collapse = '_'))
  } else{
    if (inherits(genome, "FaFile")) {
      return(paste0(c('FaFile',
                      strsplit(basename(genome$path), '\\.')[[1]][1],
                      gsub(' ', '_', tags$tax_group),
                      tags$collection,
                      'cutoff_motifs_matrix.Rdata'), collapse = '_'))
      } else {
        stop("The used genome should be FaFile or BSgenome.")
      }
    }
}

path.connect <- function(path, file){
  if(substr(path, nchar(path), nchar(path)) == .Platform$file.sep){
    path <- substr(path, 1, nchar(path) - 1)
  }
  file.path(path, file)
}

get_bg <- function(bg_method, subject, genome){
  if (bg_method == "subject"){
    bg <- getNucFreqs(subject)
  } else if (bg_method == "genome"){
    if (is.null(genome))
      stop("If bg is genome, then a genome argument must ",
           "be provided!")
    genome <- validate_genome_input(genome)
    bg <- getNucFreqs(genome)
  } else{
    bg <- rep(0.25,4)
  }
  return(bg)
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

# find the distance of start location from the nearest TSS
distance2nearestTss <- function(start.loc, genes, distance.cutoff=10000){
  if (!is.wholenumber(start.loc)) {
    stop('start location should be integer')
  }
  if (class(genes) != 'GRanges') {
    stop('genes should be GRanges')
  }
  genes.tss <- promoters(genes, upstream = 0, downstream = 1)
  dis_tss <- start.loc - genes.tss@ranges@start
  if (any(abs(dis_tss) < distance.cutoff)) {
    dis_neartest <- dis_tss[which(abs(dis_tss) == min(abs(dis_tss)))]
    if (genes.tss[which(dis_tss == dis_neartest)]@strand@values == '-') {
      dis_neartest <- -dis_neartest
    }
  }else{
    dis_neartest <- NULL
  }
  return(dis_neartest)
}

# generate random start from gene TSS location
randomStartFromGene <- function(genes, distance2Tss, random.seed){
  set.seed(random.seed)
  if (class(genes) != 'GRanges') {
    stop('genes should be GRanges')
  }
  if (!is.wholenumber(distance2Tss)) {
    stop('distance2Tss should be integer')
  }
  random_gene <- sample(genes ,1)
  random_gene_tss <- promoters(random_gene, upstream = 0, downstream = 1)
  if (random_gene@strand@values == '-') {
    start <- random_gene_tss@ranges@start - distance2Tss
  } else{
    start <- random_gene_tss@ranges@start + distance2Tss
  }
  return(start)
}


#' pwmType
#'
#' Determines type of PWM
#' @param pwm PWMatrix object
#' @return 'log','log2', or 'frequency' depending on type of pwm
#' @export
#' @keywords internal
#'
pwmType <- function(pwm) {
  # Determine whether un-logged, natural log, or log2
  if (isTRUE(all.equal(colSums(as.matrix(pwm)),
                       rep(1, length(pwm)),
                       check.attributes = FALSE))) {
    return("frequency")
  } else if (isTRUE(all.equal(colSums(2^(as.matrix(pwm)) *
                                      matrix(bg(pwm),
                                             byrow = FALSE,
                                             ncol = length(pwm),
                                             nrow = 4)),
                              rep(1, length(pwm)), tolerance = 10^-5,
                              check.attributes = FALSE))) {
    return("log2")
  } else if (isTRUE(all.equal(colSums(exp(as.matrix(pwm)) *
                                      matrix(bg(pwm),
                                             byrow = FALSE,
                                             ncol = length(pwm),
                                             nrow = 4)),
                              rep(1, length(pwm)), tolerance = 10^-5,
                              check.attributes = FALSE))) {
    return("log")
  } else {
    stop("Can't determine format of PWM -- should be numeric ",
         "frequency summing to 1 or log or log2 odds ratio")
  }
}


