[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

![Picture3](https://github.com/user-attachments/assets/328450ff-97ba-4669-b497-7e24c4105b03)

A differential expression analysis pipeline for bulk RNAseq data with Gene of Interest (GOI) data fetching into ENSEMBLE, UNIPROT & Biogrid databases for protein data and pathway analysis. 

This is NOT reccomemded for any serious DE analysis as it lacks sufficient robustness and statistical rigor. Please use a more standard FASTQC -> cutadapt -> STAR -> FeatureCounts -> [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)/[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)/[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) pipeline. 
 
## Description

Differential gene expression analysis is a method used to identify genes whose expression levels vary significantly between different biological conditions or sample groups. Here we provide a differential gene expression analysis pipeline for bulk RNAseq data, with an example from _Trypanosoma congolense_ RNAseq in fastq format, where information about the GOI's protein product and possible pathways are fetched from ENSEMBLE, UNIPROT & Biogrid databases.

### General Overview

1. Imports scRNA data
2. Quality checking FastQC reports are generated to inspect sequencing quality
3. Filters bad quality reads
4. Aligns read pairs to reference genome
5. Generate gene counts
6. Generate replicate means and condition fold changes (Will process in multi factor manner)
7. Filter genes by significance threshold
8. Queries UNIPROT ID Mapping API to gather protein accessions
9. Queries UNIPROT, ENSEMBLE, Biogrid databases for relevant information
10. Export data to MySQL database

## Getting Started

1. git clone the repo dawg
2. In ./seqdata place bulk RNAseq data in .fq.gz format. Also provide .fqfiles for metadata (example .fqfiles given for _Trypanosoma congolense_).
3. In ./refseqdata place reference genome in .bed format. Can be downloaded from [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)
4. Set significant DE fold change threshold (Default: >=2)
6. Run foldchange.sh to gather gene expression fold changes. This will also automatically run GOI_crawler.R to gather protein and pathway data from different databases to place into a MySQL database.

(Note: GOI_crawler can also be run independently on a BED file with the foldchange in the last column.)



### Dependencies

[FastqQC](https://github.com/s-andrews/FastQC)

[bedtools](https://github.com/arq5x/bedtools2)

[samtools](https://github.com/samtools/samtools)

[bowtie2](https://github.com/BenLangmead/bowtie2)

library(httr)
library(jsonlite)
library(xml2)
library(RCurl)
library(curl)
library(xml2)
library(rentrez)
library(RMySQL)
library(queryup)
library(argparser)


## Authors

[Jay Chow (Chi Lung)](https://github.com/jaychowcl/benchdeconv/)


## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the MIT License
