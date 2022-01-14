#' getJasparMotifs
#'
#' @description function to get motifs from JASPAR database
#'
#' @param species Which species?  use either jaspar code or latin name.
#' default is 'Homo sapiens'
#' @param collection Which collection to use?  default is 'CORE'
#' @param ... additional arguments to opts for
#' \code{\link[TFBSTools]{getMatrixSet}}
#'
#' @details Simply a wrapper function for \code{\link[TFBSTools]{getMatrixSet}}
#'  that calls JASPAR2020 database using \code{\link[JASPAR2020]{JASPAR2020}}
#'
#' @return \code{\link[TFBSTools]{PFMatrixList}}
#'
#' @export
#'
#' @examples
#'
#' motifs <- getJasparMotifs()
#'
#'
getJasparMotifs <- function(species = "Homo sapiens",
                            collection = "CORE", ...) {
  opts <- list()
  opts["species"] <- species
  opts["collection"] <- collection
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out))))
    names(out) <- paste(names(out), TFBSTools::name(out), sep = "_")
  return(out)
}




PWMCutoff <- setClass("PWMCutoff",
                      slots = c(pseudocounts = "numeric",
                                ID = "character",
                                name = "character",
                                matrixClass = "character",
                                strand = "character",
                                bg = "numeric",
                                tags = "list",
                                profileMatrix = "matrix",
                                score = "numeric"))
setClass("PWMCutoffList",
         contains="SimpleList",
         representation(
         ),
         prototype(
           elementType="PWMCutoff"
         )
)

setGeneric("PWMCutoff", function(object){
  standardGeneric("PWMCutoff")
})

setMethod("PWMCutoff", signature(object = "PWMatrix"),
          function(object){
            new("PWMCutoff",
                pseudocounts = object@pseudocounts,
                ID = object@ID,
                name = object@name,
                matrixClass = object@matrixClass,
                strand = object@strand,
                bg = object@bg,
                tags = object@tags,
                profileMatrix = object@profileMatrix)
          })

setMethod("PWMCutoff", signature(object = "PWMCutoff"),
          function(object){
            new("PWMCutoff",
                pseudocounts = object@pseudocounts,
                ID = object@ID,
                name = object@name,
                matrixClass = object@matrixClass,
                strand = object@strand,
                bg = object@bg,
                tags = object@tags,
                profileMatrix = object@profileMatrix,
                score = object@score)
          })

setGeneric("PWMCutoffList", function(object){
  standardGeneric("PWMCutoffList")
})

setMethod("PWMCutoffList", signature(object = "PWMatrixList"),
          function(object){
            S4Vectors:::new_SimpleList_from_list("PWMCutoffList", sapply(object, PWMCutoff))
          })

setMethod("PWMCutoffList", signature(object = "list"),
          function(object){
            S4Vectors:::new_SimpleList_from_list("PWMCutoffList", sapply(object, PWMCutoff))
          })

setGeneric("MaxScore", function(object){standardGeneric("MaxScore")})

setGeneric("MinScore", function(object){standardGeneric("MinScore")})

setGeneric("setCutoff", function(object, cutscore){standardGeneric("setCutoff")})

setMethod("ID", signature(x = "PWMCutoff"), function(x) x@ID)

setMethod("ID", signature(x = "PWMCutoffList"), function(x) lapply(x, ID))

setMethod("name", signature(x = "PWMCutoff"), function(x) x@name)

setMethod("name", signature(x = "PWMCutoffList"), function(x) lapply(x, name))

setMethod("tags", signature(x = "PWMCutoff"), function(x) x@tags)

setMethod("tags", signature(x = "PWMCutoffList"), function(x) lapply(x, tags))

setMethod("as.matrix", signature(x = "PWMCutoff"), function(x) x@profileMatrix)

setMethod("as.matrix", signature(x = "RangedSummarizedExperiment"), function(x) {
  out_mat <- x@assays@data[[1]]
  colnames(out_mat) <- x@colData@listData$name
  return(out_mat)})

setMethod("length", signature(x = "PWMCutoff"), function(x) ncol(x@profileMatrix))

setMethod("MaxScore", signature(object = "PWMCutoff"),
          function(object){
            sum(apply(object@profileMatrix, 2, max))
          })

setMethod("MinScore", signature(object = "PWMCutoff"),
          function(object){
            sum(apply(object@profileMatrix, 2, min))
          })

setMethod("setCutoff", signature(object = "PWMCutoffList",
                                 cutscore = "matrix"),
          function(object, cutscore){
            n_bits <- min(c(nchar(as.character(dim(cutscore)[2])), 7))
            pwm_list <- PWMCutoffList(lapply(object@listData, function(pwm){
              pwm_name <- paste0(c(pwm@ID, pwm@name), collapse = '_')
              scores <- cutscore[match(pwm_name, names(object@listData)),]
              cutoffs <- unlist(lapply(2:(n_bits - 1), function(n){quantile(scores, 1 - .1^n, type = 1)}))
              names(cutoffs) <- as.numeric(paste0('1e-', 2:(n_bits - 1)))
              pwm@score <- cutoffs
              pwm
              }))
            return(pwm_list)
          })
