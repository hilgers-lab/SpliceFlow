#' Retrieve introns from annotation
#'
#' @param refAnnotation.path path to reference annotation in gtf format
#'
#' @return regionFlow-class
#' @export
#' @import dplyr
#' @import GenomicRanges
#' @import IRanges
#' @import BiocGenerics
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom rtracklayer import.gff
#' @importFrom txdbmaker makeTxDbFromGFF
#' @examples
retrieveIntrons <- function(refAnnotation.path, junctions) {
  txdb <- txdbmaker::makeTxDbFromGFF(file = refAnnotation.path, format = "gtf")
  message("Retrieving intron annotation")
  introns <- GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)
  genes <- rtracklayer::import.gff(refAnnotation.path)
  # remove introns overlapping exons exactly
  utrs.annot <- genes[genes$type %in% c("five_prime_utr", "three_prime_utr")]
  exons.annot <- genes[genes$type %in% c("exon")]
  message("Cleaning intron reference")
  # remove exons with full match of UTRs
  exons <- IRanges::subsetByOverlaps(exons.annot, utrs.annot, invert = TRUE)
  # remove introns fully overlapping exons
  intron.annotation <- IRanges::subsetByOverlaps(introns, exons, invert = TRUE)
  # remove introns of less than 10 nt
  intron.annotation <- intron.annotation[GenomicRanges::width(intron.annotation) >= 10]
  intron.annotation <- unlist(intron.annotation)
  intron.annotation$transcript_id <- names(intron.annotation)
  names(intron.annotation) <- NULL
  intron.annotation <- addGeneName(genes, intron.annotation) %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  message("Create intron database")
  # subset intron annotation using junctions from the tissue/sample used.
  intron.annotation <- IRanges::subsetByOverlaps(intron.annotation,
                                        junctions,
                                        maxgap = 1)
  # create intron_id
  introns.db <- sort(intron.annotation) %>%
    as.data.frame(.) %>%
    group_by(gene_id) %>%
    mutate(intron_id = paste0(gene_id, ":", "intron" , row_number())) %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  message("Create exon database")
  # create exon id
  exons.db <- sort(exons) %>%
    as.data.frame(.) %>%
    group_by(gene_id) %>%
    mutate(exon_id = paste0(gene_id, ":", "exon" , row_number())) %>%
    filter(gene_biotype == "protein_coding") %>%
    dplyr::select(c(
      seqnames,
      start,
      end,
      strand,
      transcript_id,
      gene_id,
      gene_name,
      exon_id
    )) %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  message("get Exon region to quantify")
  exon_starts <- getExonStart(exons.db)
  intron_ends <- getIntronFraction(introns.db)
  message("Assign exon to intron")
  exon.annotation <- assignExon(BiocGenerics::unique(exon_starts), BiocGenerics::unique(intron_ends))
  message("Assign intron to exon")
  intron.annotation <- assignIntron(BiocGenerics::unique(exon_starts), BiocGenerics::unique(intron_ends))
  #show(exon.annotation)
  #show(intron.annotation)
  features_combined <- c(exon.annotation, intron.annotation)
  message("Identify redudant splice regions")
  features_combined <- as.data.frame(features_combined)
  features_combined <- features_combined %>%
    dplyr::mutate(exon_intron_id = paste0(intron_id, "|", exon_id))
  splicing_region <- GenomicRanges::makeGRangesFromDataFrame(features_combined, keep.extra.columns= TRUE)
  splicing_region.list <- GenomicRanges::split(splicing_region, f=splicing_region$exon_intron_id)
  splicing_region.list <- GenomicRanges::reduce(splicing_region.list)
  splicing_region.list <- BiocGenerics::unlist(splicing_region.list)
  splicing_region <- unique(splicing_region.list)
  # handle same exon multiple intron assignments
  # take intron with 5 start closest to exon
  message("fix redudant exon assignment ")
  show(splicing_region)
  splicing_region$region_id <- names(splicing_region)
  splicing_region$intron_id <- stringr::str_split_fixed(splicing_region$region_id, pattern = "\\|", n = 2)[,1]
  splicing_region$exon_id <- stringr::str_split_fixed(splicing_region$region_id, pattern = "\\|", n = 2)[,2]
  rownames(splicing_region) <- NULL
  show(splicing_region)
  splicing_region <- splicing_region %>%
    as.data.frame(.) %>%
    group_by(exon_id) %>%
    filter(width == min(width)) %>%
    dplyr::ungroup() %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns= TRUE)
  intron_features <- list()
  intron_features$exon <- exon.annotation
  intron_features$intron <- intron.annotation
  intron_features$all <- splicing_region
  return(intron_features)
}

#' Title
#'
#' @param exons
#' @param introns
#'
#' @return
#' @export
#' @import GenomicRanges dplyr IRanges
#' @examples
assignIntron <- function(exons, introns) {
  # annotate exon to intron
  hits.intron.exon <-
    GenomicRanges::findOverlaps(introns, exons, maxgap = 1)
  trimmed.intron.to_exon <- introns[S4Vectors::queryHits(hits.intron.exon)]
  trimmed.intron.to_exon$exon_start <- GenomicRanges::start(exons[S4Vectors::subjectHits(hits.intron.exon)])
  trimmed.intron.to_exon$exon_id <- exons[S4Vectors::subjectHits(hits.intron.exon)]$exon_id
  trimmed.intron.to_exon$exon_start_neg <- GenomicRanges::end(exons[S4Vectors::subjectHits(hits.intron.exon)])
  trimmed.intron.to_exon <- trimmed.intron.to_exon %>%
    as.data.frame() %>%
    mutate(distance_to_exon =
             ifelse(strand == "+", exon_start - end, exon_start_neg - start)) %>%
  filter(strand == "+" & distance_to_exon %in% c(1) |
           strand == "-" & distance_to_exon %in% c(-1)) %>%
    dplyr::select(c(
      seqnames,
      start,
      end,
      strand,
      gene_id,
      transcript_id,
      intron_id,
      exon_id
    )) %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  return(trimmed.intron.to_exon)
}
#' Title
#'
#' @param exons
#' @param introns
#'
#' @return
#' @export
#' @import GenomicRanges dplyr IRanges
#' @examples
assignExon <- function(exons, introns) {
  # annotate exon to intron
  hits.intron.exon <-
    GenomicRanges::findOverlaps(introns, exons, maxgap = 1)
  trimmed.exon.to_intron <-
    exons[S4Vectors::subjectHits(hits.intron.exon)]
  trimmed.exon.to_intron$intron_end <-
    GenomicRanges::end(introns[S4Vectors::queryHits(hits.intron.exon)])
  trimmed.exon.to_intron$intron_id <-
    introns[S4Vectors::queryHits(hits.intron.exon)]$intron_id
  trimmed.exon.to_intron$intron_end_neg <-
    GenomicRanges::start(introns[S4Vectors::queryHits(hits.intron.exon)])
  exons <- trimmed.exon.to_intron %>%
    as.data.frame() %>%
    mutate(distance_to_exon =
             ifelse(strand == "+", intron_end - start, intron_end_neg - end)) %>%
    filter(distance_to_exon %in% c(-1,1)) %>%
    dplyr::select(c(
      seqnames,
      start,
      end,
      strand,
      gene_id,
      transcript_id,
      intron_id,
      exon_id
    )) %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  return(exons)
}



