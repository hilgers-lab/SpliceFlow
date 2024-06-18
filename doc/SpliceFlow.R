## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
require(dplyr)

## ----setup--------------------------------------------------------------------
library(SpliceFlow)
library(GenomicAlignments)
library(Rsamtools)
library(SummarizedExperiment)

## -----------------------------------------------------------------------------
ref <- system.file("data", "protein_coding.genes.gtf.gz", package = "SpliceFlow")
data("neuronal.junctions", package = "SpliceFlow")
regionsFlow <- prepareFlowRegions(ref, neuronal_junctions )
show(regionsFlow)

## -----------------------------------------------------------------------------
pathBams <- list.files(system.file("data", package = "SpliceFlow"), 
                       pattern = "\\.bam$", full.names = TRUE)
bamFiles = Rsamtools::BamFileList(pathBams)
names(bamFiles) <- c("nascent_R1", "nascent_R2", "polyA_R1", "polyA_R2")
# regions to quantify
# counts 
se = GenomicAlignments::summarizeOverlaps(regionsFlow@quantFlow,
                       bamFiles,
                       singleEnd = TRUE,
                       ignore.strand = TRUE)
rowData(se)$tx_name <- paste0( rowData(se)$gene_id, ":E", rowData(se)$exon_id) 
rowData(se)$transcript_id <- rowData(se)$exon_intron_region
colData(se)$condition = factor(c("nascent", "nascent", "mrna", "mrna"))
result <- SpliceFlow(se, regionsFlow)
# results stored in rowData
rowData(result)

