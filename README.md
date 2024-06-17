# <img src="https://github.com/hilgers-lab/SpliceFlow/blob/main/data/logo.webp" alt="Logo" width="20%" align="left"> SpliceFlow

## Quantification of splicing efficiencies from total RNA-seq data
<!-- badges: start -->

  ![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/hilgers-lab/LASER)
[![Maintained?](https://img.shields.io/badge/Maintained%3F-Yes-brightgreen)](https://github.com/hilgers-lab/LASER/graphs/contributors)
[![Install](https://img.shields.io/badge/Install-Github-brightgreen)](#installation)
  ![GitHub](https://img.shields.io/github/license/hilgers-lab/LASER)
  [![DOI](https://zenodo.org/badge/580128861.svg)](https://zenodo.org/badge/latestdoi/580128861)

  <!-- badges: end -->
  
  Splice flow determines the efficiency of splicing of invidividual exons. 


  ### Installation

  ```
  install.packages("devtools")
  devtools::install_github("hilgers-lab/SpliceFlow")
  ```
  ### Input files
  * Genome Alignment bam files [STAR]
  * Reference annotation in gtf format. Example file [here](https://github.com/hilgers-lab/LASER/blob/master/inst/exdata/dm6.annot.gtf.gz)
  * Short read sequencing SJ.out files [STAR](https://github.com/alexdobin/STAR). Example file in [here](https://github.com/hilgers-lab/LASER/blob/master/inst/exdata/short_read_junctions.SJ.out.tab). We recommend to pull SJ.out into a single SJ.out from many experiments and filter by min counts.

  ### Usage
  A step by step guide for data processing and identification of exon-couplings can be found
  [here](https://hilgers-lab.github.io/LASER/docs/LASER.html).

  ```
  library(SpliceFlow)
  vignette("SpliceFlow")
  ```


  # Release

  Initial Release 0.1.0

  Release date: Xth X 202X
  This release corresponds to the SpliceFlow version used by [Alfonso-Gonzalez et al. 202X]()

  ## Contact

  Developer Carlos Alfonso-Gonzalez. For questions or feedback you can contact:

    alfonso@ie-freiburg.mpg.de
