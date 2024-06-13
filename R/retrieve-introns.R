#' Retrieve introns from annotation
#'
#' @param refAnnotation.path path to reference annotation in gtf format
#'
#' @return
#' @export
#' @import dplyr
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
  intron.annotation <- subsetByOverlaps(intron.annotation, junctions, maxgap = 1)
  # create intron_id
  introns.db <- sort(intron.annotation) %>%
    as.data.frame(.) %>%
    group_by(gene_id) %>%
    mutate(intron_id = paste0(gene_id, ":", "intron" , row_number())) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
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
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  message("Assign exon to intron")
  exon.annotation <- assignExon(exons.db, introns.db)
  message("Assign intron to exon")
  intron.annotation <- assignIntron(exons.db, introns.db)
  intron_features <- list()
  intron_features$exon <- exons.db
  intron_features$intron <- introns.db
  intron_features$all <- intron.annotation
  return(intron_features)
}
assignIntron <- function(exons, introns) {
  # annotate exon to intron
  hits.intron.exon <-
    findOverlaps(introns, exons, maxgap = 1)
  trimmed.intron.to_exon <- introns[queryHits(hits.intron.exon)]
  trimmed.intron.to_exon$exon_start <-
    start(exons[subjectHits(hits.intron.exon)])
  trimmed.intron.to_exon$exon_id <-
    exons[subjectHits(hits.intron.exon)]$exon_id
  trimmed.intron.to_exon$exon_start_neg <-
    end(exons[subjectHits(hits.intron.exon)])
  trimmed.intron.to_exon <-
    trimmed.intron.to_exon %>% as.data.frame() %>%
    mutate(distance_to_exon =
             ifelse(strand == "+", exon_start - end, exon_start_neg - start)) %>%
    filter(distance_to_exon %in% c(-1, 1)) %>%
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
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  return(trimmed.intron.to_exon)
}
assignExon <- function(exons, introns) {
  # annotate exon to intron
  hits.intron.exon <-
    findOverlaps(introns, exons, maxgap = 1)
  trimmed.exon.to_intron <-
    exons[subjectHits(hits.intron.exon)]
  trimmed.exon.to_intron$intron_end <-
    end(introns[queryHits(hits.intron.exon)])
  trimmed.exon.to_intron$intron_id <-
    introns[queryHits(hits.intron.exon)]$intron_id
  trimmed.exon.to_intron$intron_end_neg <-
    start(introns[queryHits(hits.intron.exon)])
  exons <- trimmed.exon.to_intron %>%
    as.data.frame() %>%
    mutate(distance_to_exon =
             ifelse(strand == "+", intron_end - start, intron_end_neg - end)) %>%
    filter(distance_to_exon %in% c(-1, 1)) %>%
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
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  return(exons)
}
