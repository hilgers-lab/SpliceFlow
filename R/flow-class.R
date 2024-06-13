### RegionsFlow Class Definition ###
#' S4 class for intron annotation data for a specific annotation version
#'
#' @slot reference.introns A GRanges object. The intron ranges annotated with the
#'   promoter information.
#' @slot introns A GRanges object. intron coordinates
#' @slot exonFlow A GRanges object. exon coordinates
#' @slot metadata data.frame specifying the links between the two.
#'   gene id and internal promoter state
#' @importFrom GenomicRanges GRanges
#'
#' @name RegionsFlow-class
#' @rdname RegionsFlow-class
#' @exportClass RegionsFlow

setClass(
  "RegionsFlow",
  slots = c(
    reference.introns = "GRanges",
    intronFlow = "GRanges",
    exonFlow = "GRanges",
    metadata = "data.frame"
  ),
  prototype = list(
    reference.introns = GenomicRanges::GRanges(),
    intronFlow = GenomicRanges::GRanges(),
    exonFlow = GenomicRanges::GRanges(),
    metadata= data.frame()
  )
)
###################
### Constructor ###

#' @param reference.introns A GRanges object containing annotated intron ranges
#' @param intronFlow A data.frame containing mapping between transcript,
#'   TSS, promoter and gene ids
#' @param exonFlow A GRanges object containing promoter coordinates
#' @param metadata A GRanges object containing promoter coordinates

#' @name RegionsFlow
#' @rdname RegionsFlow-class
#'
#' @importFrom methods new
#' @importFrom GenomicRanges GRanges
#'
#' @export
#' @return A promoter annotation object with three slots: intronRanges,
#'   promoterIdMapping and promoter Coordinates
RegionsFlow <-
  function(reference.introns = GenomicRanges::GRanges(),
           intronFlow = GenomicRanges::GRanges(),
           exonFlow = GenomicRanges::GRanges(),
           metadata = data.frame()) {
    new(
      "RegionsFlow",
      reference.introns = reference.introns,
      intronFlow = intronFlow,
      exonFlow = exonFlow,
      metadata = metadata)
  }

setValidity("RegionsFlow", function(object) {
  check <- TRUE
  if (is(object@reference.introns, 'GRanges') == FALSE) {
    check <- FALSE
  }
  if (is(object@intronFlow, 'GRanges') == FALSE) {
    check <- FALSE
  }
  if (is(object@exonFlow, 'GRanges') == FALSE) {
    check <- FALSE
  }
  return(check)
})

###############
### Getters ###

#' @param x A RegionsFlow object
#'
#' @describeIn RegionsFlow-class Getter for reference.introns
#' @exportMethod reference.introns

setGeneric("reference.introns", function(x) standardGeneric("reference.introns"))

#' @describeIn intronFlow-class Getter for reference.introns
#' @aliases reference.introns,intronFlow-method

setMethod("reference.introns", "RegionsFlow", function(x) x@reference.introns)

#' @describeIn intronFlow-class Getter for intronFlow
#' @exportMethod intronFlow

setGeneric("intronFlow",
           function(x) standardGeneric("intronFlow"))

#' @describeIn intronFlow-class Getter for intronFlow
#' @aliases intronFlow,intronFlow-method

setMethod("intronFlow", "RegionsFlow",
          function(x) x@intronFlow)

#' @describeIn intronFlow-class Getter for intronFlow
#' @exportMethod intronFlow

setGeneric("exonFlow",
           function(x) standardGeneric("exonFlow"))

#' @describeIn intronFlow-class Getter for exonFlow
#' @aliases exonFlow,intronFlow-method

setMethod("exonFlow", "RegionsFlow",
          function(x) x@exonFlow)


###############
### Setters ###

#' @param value refintronFlow, intron or exon to
#'   be assigned
#'
#' @describeIn refintronFlow-class Setter for intron
#' @exportMethod 'refintronFlow<-'
#' @importFrom methods validObject

setGeneric("reference.introns<-",
           function(x, value) standardGeneric("reference.introns<-"))

#' @describeIn RegionsFlow-class Setter for intronRanges
#' @aliases 'refintronFlow<-',RegionsFlow-method

setMethod("reference.introns<-", "RegionsFlow", function(x, value) {
  x@reference.introns <- value
  validObject(x)
  x
})

#' @describeIn RegionsFlow-class Setter for intronFlow
#' @exportMethod 'intronFlow<-'
#' @importFrom methods validObject

setGeneric("intronFlow<-",
           function(x, value) standardGeneric("intronFlow<-"))

#' @describeIn RegionsFlow-class Setter for intronFlow
#' @aliases 'intronFlow<-',RegionsFlow-method

setMethod("intronFlow<-", "RegionsFlow", function(x, value) {
  x@intronFlow <- value
  validObject(x)
  x
})

#' @describeIn RegionsFlow-class Setter for intronFlow
#' @exportMethod 'intronFlow<-'
#' @importFrom methods validObject

setGeneric("exonFlow<-",
           function(x, value) standardGeneric("exonFlow<-"))

#' @describeIn intronFlow-class Setter for intronFlow
#' @aliases 'intronFlow<-',RegionsFlow-method

setMethod("exonFlow<-", "RegionsFlow", function(x, value) {
  x@exonFlow <- value
  validObject(x)
  x
})

#' @describeIn RegionsFlow-class Setter for intronFlow
#' @exportMethod 'intronFlow<-'
#' @importFrom methods validObject

setGeneric("metadata<-",
           function(x, value) standardGeneric("metadata<-"))

#' @describeIn intronFlow-class Setter for intronFlow
#' @aliases 'intronFlow<-',RegionsFlow-method

setMethod("metadata<-", "RegionsFlow", function(x, value) {
  x@metadata <- value
  validObject(x)
  x
})

#' @describeIn RegionsFlow-class Getter for metadata
#' @exportMethod intronFlow

setGeneric("metadata",
           function(x) standardGeneric("metadata"))

#' @describeIn RegionsFlow-class Getter for metadata
#' @aliases metadata,regionFlow-method

setMethod("metadata", "RegionsFlow",
          function(x) x@metadata)

#' @rdname RegionsFlow-class
#' @aliases show
setMethod("show", "RegionsFlow", function(object) {
  cat("RegionsFlow object\n")
  cat("Introns per isoform:\n")
  print(object@reference.introns)
  cat("\nIntron Flow:\n")
  print(object@intronFlow)
  cat("\nExon Flow:\n")
  print(object@exonFlow)
  cat("\nMetadata:\n")
  cat("Intron-Exon pairs that will be quantified \n")
  print(head(object@metadata))
})


