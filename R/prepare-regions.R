prepareFlowRegions <- function(refAnnotation.path, junctions){
  intron_features<- retrieveIntrons(refAnnotation.path, junctions)
  regionsFlow <- RegionsFlow()
  # trim exons
  exon_starts <- getExonStart(intron_features$exon)
  # trim introns
  intron_ends <- getIntronFraction(intron_features$intron)
  # Make pair regions and remove redundancy
  regions <- test$all %>% as.data.frame(.) %>%
    mutate(exon_intron_id = paste0(intron_id, "|", exon_id)) %>%
    makeGRangesFromDataFrame(., keep.extra.columns= TRUE) %>%
    split(., f=.$exon_intron_id) %>% range() %>%
    unlist(.) %>% unique(.)
  regions$region_id <- names(regions)
  regions$intron_id <- stringr::str_split_fixed(regions$region_id, pattern = "\\|", n = 2)[,1]
  regions$exon_id <- stringr::str_split_fixed(regions$region_id, pattern = "\\|", n = 2)[,2]
  # prepare metadata column
  metadata <- regions %>%
    as.data.frame(.) %>%
    mutate(exon_intron_region = paste0(seqnames, ":", start, "-", end, "|", strand)) %>%
    dplyr::select(c(region_id, intron_id, exon_id, exon_intron_region)) %>%
    mutate(gene_id = stringr::str_split_fixed(region_id, pattern = "[|:]", n = 4)[,1]) %>%
    mutate(region_id = paste0(stringr::str_split_fixed(region_id, pattern = "[|:]", n = 4)[,1],":",
                              stringr::str_split_fixed(region_id, pattern = "[|:]", n = 4)[,2],
                              "-",
                              stringr::str_split_fixed(region_id, pattern = "[|:]", n = 4)[,4]) )
  rownames(metadata) <- NULL
  exon_starts <- unique(exon_starts[exon_starts$exon_id %in% regions$exon_id ])
  intron_ends <- unique(intron_ends[intron_ends$intron_id %in% regions$intron_id ])
  reference.introns(regionsFlow) <- intron_features$all
  intronFlow(regionsFlow) <- intron_ends
  exonFlow(regionsFlow) <- exon_starts
  metadata(regionsFlow) <- metadata
  return(regionsFlow)
}


