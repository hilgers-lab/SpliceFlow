#' Add gene_id and gene_name
#'
#' @param annotation ENSEMBL annotation in granges
#'
#' @return trimmed exon
#' @export
#' @import dplyr
#' @examples
addGeneName <- function(annotation, intron_database) {
  gene_to_gene_name_map <- annotation[annotation$type == "transcript"] %>%
    as.data.frame(.) %>%
    dplyr::select(c(gene_id, gene_name, transcript_id))
  names(intron_database) <- NULL
  intron_database <- intron_database %>%
    as.data.frame(.) %>%
    left_join(., gene_to_gene_name_map,
              by = "transcript_id")
  return(intron_database)
}
#' Retrieve intron 3SS last 0.3 nucleotide fraction
#'
#' @param grIntron GenomicRanges intron
#'
#' @return trimmed intron
#' @export
#' @import dplyr
#' @examples
getIntronFraction <- function(grIntron) {
  pos <- grIntron[GenomicRanges::strand(grIntron) == "+"]
  intron_width <- GenomicRanges::width(pos)
  nucleotides_to_take <- intron_width * 0.3
  GenomicRanges::start(pos) <- GenomicRanges::end(pos) - nucleotides_to_take
  neg <- grIntron[!GenomicRanges::strand(grIntron) == "+"]
  intron_width <- GenomicRanges::width(neg)
  nucleotides_to_take <- intron_width * 0.3
  GenomicRanges::end(neg) <- GenomicRanges::start(neg) + nucleotides_to_take
  trimmed_intron <- c(pos, neg)
  return(trimmed_intron)
}
#' Retrieve Exon 5SS first 25 nucleotides
#'
#' @param grExon GenomicRanges exon
#'
#' @return trimmed exon
#' @export
#' @import dplyr
#' @examples
getExonStart <- function(x) {
  pos <- x[GenomicRanges::strand(x) == "+"]
  GenomicRanges::end(pos) <- GenomicRanges::start(pos) + 25
  neg <- x[!GenomicRanges::strand(x) == "+"]
  GenomicRanges::start(neg) <- GenomicRanges::end(neg) - 25
  trimmed_exon <- c(pos, neg)
  return(trimmed_exon)
}

