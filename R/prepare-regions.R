prepareFlowRegions <- function(refAnnotation.path, junctions){
  intron_features<- retrieveIntrons(refAnnotation.path, junctions)
  regionsFlow <- RegionsFlow()
  # Make pair regions and remove redundancy
  # prepare metadata column
  message("preparing regions to quantify quantFlow")
  # prepare quantFlow regions
  # regions to quantify
  exonFlow <- intron_features$exon %>%
    as.data.frame(.) %>%
    dplyr::inner_join(., intron_features$all %>% as.data.frame(.) %>%
                dplyr::select(c(exon_id, intron_id ,region_id))) %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  intronFlow <- intron_features$intron %>%
    as.data.frame(.) %>%
    dplyr::inner_join(., intron_features$all %>% as.data.frame(.) %>%
                dplyr::select(c(exon_id, intron_id ,region_id))) %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  regionsToQuantify <-
    c(exonFlow, intronFlow)
  regionsToQuantify$id <- seq(1,length(regionsToQuantify))
  regionsToQuantify$former_exon_id <- regionsToQuantify$exon_id
  regionsToQuantify$exon_id <- regionsToQuantify$id
  regionsToQuantify$transcript_id <- regionsToQuantify$region_id

  message("output regionsFlow object")
  FlowRegion(regionsFlow) <- intron_features$all
  intronFlow(regionsFlow) <- intronFlow
  exonFlow(regionsFlow) <- exonFlow
  quantFlow(regionsFlow) <- regionsToQuantify
  metadata(regionsFlow) <- intron_features$all %>% as.data.frame(.)
  return(regionsFlow)
}


