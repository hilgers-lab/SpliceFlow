#' SpliceFlow calculate splicing efficiencies
#'
#' @param se_counts : counts from exons and introns
#' @param regionsFlow: regionFlow-object containing exons and pairs
#'
#' @return SpliceFlow calculations
#' @export
#' @importFrom DEXSeq DEXSeqDataSetFromSE
#' @importFrom DEXSeq DEXSeq
#' @import dplyr
#' @examples
SpliceFlow <- function(se_counts, regionsFlow){
  dxd = DEXSeqDataSetFromSE(se_counts, design = ~ sample + exon + condition:exon)
  dxr = DEXSeq(dxd)
  # prepare dataset
  # exons
  hits.exons <- GenomicRanges::findOverlaps(dxr, regionsFlow@exonFlow, type = "equal")
  tx_exons <- dxr[S4Vectors::queryHits(hits.exons),]
  tx_exons <- as.data.frame(tx_exons)
  tx_exons$type <- "exon"
  tx_exons$region_id <- regionsFlow@exonFlow[S4Vectors::subjectHits(hits.exons),]$region_id
  tx_exons <- tx_exons %>%
    dplyr::select(c(1,2,8,9,11,12,13,15,
                    starts_with("countData.", ),
                    region_id)) %>%
    dplyr::rename_with(.fn = ~ paste0("exon_", .),
                .cols = -c(groupID,
                           featureID,
                           region_id))
  # reformat columns
  # introns
  hits.introns <- GenomicRanges::findOverlaps(dxr, regionsFlow@intronFlow, type = "equal")
  tx_introns <- dxr[S4Vectors::queryHits(hits.introns),]
  tx_introns <- as.data.frame(tx_introns)
  tx_introns$type <- "introns"
  tx_introns$region_id <- regionsFlow@intronFlow[S4Vectors::subjectHits(hits.introns),]$region_id
  tx_introns <- tx_introns %>%
    dplyr::select(c(8,9,11,12,13,15,
                    starts_with("countData.", ),
                    region_id)) %>%
    dplyr::rename_with(.fn = ~ paste0("intron_", .),
                .cols = -c(region_id))
  splice_region_data <- dplyr::inner_join(tx_exons,tx_introns, by = "region_id")
  # build summarized experiment
  exon_counts <- splice_region_data %>%
    dplyr::select(contains("exon_countData."), region_id)  %>%
    dplyr::rename_with(.fn = ~ gsub("exon_", "", .),
                .cols = -c(region_id))
  rownames(exon_counts) <- exon_counts$region_id
  exon_counts$region_id <- NULL
  intron_counts <- splice_region_data %>%
    dplyr::select(contains("intron_countData."), region_id)  %>%
    dplyr::rename_with(.fn = ~ gsub("intron_", "", .),
                .cols = -c(region_id))
  rownames(intron_counts) <- intron_counts$region_id
  intron_counts$region_id <- NULL
  # splicing effiency calculations from normalized data per condition
  splicing_calculations <- splice_region_data %>%
    dplyr::select(c(,3,14,4,15)) %>%
    dplyr::mutate(control_splicing.efficiency = 1 - .[[2]] / .[[1]]) %>%
    dplyr::mutate(experimental_splicing.efficiency = 1 - .[[4]] / .[[3]])
  result <- SummarizedExperiment::SummarizedExperiment(assays = list(
    exonCounts = exon_counts,
    intronCounts = intron_counts
  ))
  rowData(result) <- splicing_calculations
  rowData(result)$intron <- splice_region_data %>%
    dplyr::select(contains("intron_"), region_id) %>%
    dplyr::rename_with(.fn = ~ gsub("intron_", "", .),
                .cols = -c(region_id)) %>%
    GenomicRanges::makeGRangesFromDataFrame(.)
  rowData(result)$exon <- splice_region_data %>%
    dplyr::select(contains("exon_"), region_id) %>%
    dplyr::rename_with(.fn = ~ gsub("exon_", "", .),
                .cols = -c(region_id)) %>%
    GenomicRanges::makeGRangesFromDataFrame(.)
  return(result)
}
