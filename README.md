[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
# DE_RNAseq_crawler

A differential expression analysis pipeline for bulk RNAseq data with Gene of Interest (GOI) data fetching into ENSEMBLE, UNIPROT & Biogrid databases for protein data and pathway analysis. 

This is NOT reccomemded for any serious DE analysis as it lacks sufficient robustness and statistical rigor. Please use a more standard FASTQC -> cutadapt -> STAR -> FeatureCounts -> [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)/[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)/[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) pipeline. 
 
## Description

Differential gene expression analysis is a method used to identify genes whose expression levels vary significantly between different biological conditions or sample groups. Here we provide a differential gene expression analysis pipeline for bulk RNAseq data, with an example from _Trypanosoma congolense_ RNAseq in fastq format, where information about the GOI's protein product and possible pathways are fetched from ENSEMBLE, UNIPROT & Biogrid databases.


## Getting Started



### Adding new deconvolution methods


### Dependencies


### Installing




## Authors

[Jay Chow (Chi Lung)](https://github.com/jaychowcl/benchdeconv/)


## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the MIT License


## References
