

#import libraries
library(httr)
library(jsonlite)
library(xml2)
#install.packages("RCurl")
library(RCurl)
library(curl)
library(xml2)
library(rentrez)
library(RMySQL)
# install.packages("queryup")
library(queryup)

#parse args for argparser
library(argparser)
input_args <- arg_parser(GOI_crawler - Info fetching from databases on GOIs for DE_RNAseq_crawler)
input_args <- add_argument(input_args, "--goi", help="Fold change output file from foldchange.sh",
                           type="character",
                           default = "./RESULTS/FoldChanges/foldchange.txt")
input_args <- add_argument(input_args, "--thresh", help="Fold change output file from foldchange.sh",
                           type="character",
                           default = "./RESULTS/FoldChanges/foldchange.txt")

argv <- parse_args(input_args)


print("Crawling for GOIs with threshold ")
#Set Sys.sleep() time
sleeptime <- 0.1

################# FUNCTIONS #################
#capitalize first letter (for species renaming for biogrid)
capitalize <- function(x){
  paste(toupper(substring(x, 1, 1)),
        tolower(substring(x, 2, nchar(x))),
        sep = "")
}


################# Import Query Genes  #################
#original genes
# querygenes <- c("ENSMUSG00000015002","ENSMUSG00000017548","ENSMUSG00000032333","ENSMUSG00000036202","ENSMUSG00000041272","ENSMUSG00000050953","ENSMUSG00000058589","ENSMUSG00000064722","ENSMUSG00000074830","ENSMUSG00000078592")

#read query from csv file and convert to list (UNCOMMENT TO USE THIS FUNCTION)
querygenes <-read.csv("./GeneIdentifiers.csv", header=FALSE)
querygenes <-as.character(querygenes)


################# Creating Dataframe Columns  #################
#creating intial vectors to create dataframe later
##ensembl lookup
commonname <- c()
genedescription <- c()
species <- c()
strand <- c()
biotype <- c()
objecttype <- c()
chromosome <- c()
startpos <- c()
endpos <- c()
##ensembl seq
seq <- c()
##uniprot id mapping
uniprotid <- rep(NA ,length(querygenes))
##uniprot proteins
proteinname <- rep(NA, length(querygenes))
proteinseq <- rep(NA, length(querygenes))
keywords <- rep(NA, length(querygenes))
##biogrid tax id
taxonid <- rep(NA, length(querygenes))
##biogrid interactions
interactswith <- rep(NA, length(querygenes))
##uniprot data
proteinlength<- rep(NA ,length(querygenes))
proteinmass<- rep(NA ,length(querygenes))
proteinfunction <- rep(NA ,length(querygenes))
proteinfeatures <- rep(NA ,length(querygenes))
proteinUniParc <- rep(NA ,length(querygenes))
proteinsubunit <- rep(NA ,length(querygenes))
proteindevstage <- rep(NA ,length(querygenes))
proteintissuespec <- rep(NA ,length(querygenes))
proteingo <- rep(NA ,length(querygenes))
proteingoid <- rep(NA ,length(querygenes))
proteinsubcellloci <- rep(NA ,length(querygenes))
proteinmodres <- rep(NA ,length(querygenes))
proteinptm <- rep(NA ,length(querygenes))
proteinfam <- rep(NA ,length(querygenes))
proteinrefseq <- rep(NA ,length(querygenes))
proteinpdb <- rep(NA ,length(querygenes))
proteinbiogrid <- rep(NA ,length(querygenes))
proteinstring <- rep(NA ,length(querygenes))


################# Querying ENSEMBL API  #################
#connect to ensembl api for lookup for each gene in query
ensemblapibase <- "https://rest.ensembl.org/lookup/id/"
ensemblapicontent <- "?content-type=application/json"
ensemblapiseqbase <- "https://rest.ensembl.org/sequence/id/"
ensemblapiseqcontent <- "?content-type=text/plain"

#for loop for querying every gene in query to ensembl api and appending desired results to initial list
for (query in querygenes){
  #query ensembl-lookup 
  cat("Gathering Ensembl lookup data for: ", query, "\n")
  ensemblquery <- paste(ensemblapibase, query, ensemblapicontent, sep="")
  ensembllookup <- GET(ensemblquery)
  ensembljson <- fromJSON(toJSON(content(ensembllookup)))
  ensembljson <- as.data.frame(t(ensembljson))
  
  #sleep for 0.1s to avoid being rate limited
  Sys.sleep(sleeptime)
  
  #extracting desired data
  commonname <- c(commonname, as.character(ensembljson$display_name))
  genedescription <- c(genedescription, as.character(ensembljson$description))
  species <- c(species, as.character(ensembljson$species))
  strand <- c(strand, as.character(ensembljson$strand))
  biotype <- c(biotype, as.character(ensembljson$biotype))
  objecttype <- c(objecttype, as.character(ensembljson$object_type))
  chromosome <- c(chromosome, as.character(ensembljson$seq_region_name))
  startpos <- c(startpos, as.character(ensembljson$start))
  endpos <- c(endpos, as.character(ensembljson$end))
  
  #query ensembl-sequence
  cat("Gathering Ensembl sequence data for: ", query, "\n")
  ensemblquery1 <- paste(ensemblapiseqbase, query, ensemblapiseqcontent, sep="")
  ensembllookup1 <- GET(ensemblquery1)
  ensemblseq <- as.character(ensembllookup1)
  seq <- c(seq, ensemblseq)
  
  #sleep for 0.1s to avoid being rate limited
  Sys.sleep(sleeptime)
}


################# Querying UNIPROT ID Mapping API  #################
#get uniprot ID for query genes via id mapping api
count <- 1
for (query in querygenes){
  #query uniprot ID mapping to get jobID
  cat("Submitting ID mapping job to Uniprot for: ", query, "\n")
  idmappost <- postForm(
    "https://rest.uniprot.org/idmapping/run",
    from="Ensembl",
    to="UniProtKB",
    ids=query
  )
  
  #sleep to wait for job completion
  Sys.sleep(1)
  
  #take id mapping jobID to get results
  cat("Getting Uniprot ID from jobs for : ", query, "\n")
  idmapjobid <- as.character(fromJSON(idmappost))
  idmapapi <- paste("https://rest.uniprot.org/idmapping/results/", idmapjobid, sep="")
  h <- basicTextGatherer()
  curlPerform(url=idmapapi, writefunction=h$update)
  idmapped <- h$value()
  idmappedjson <- fromJSON(idmapped)
  
  if(is.null(idmappedjson$failedIds) == TRUE){
    idmappedjson <- as.data.frame(idmappedjson)
    uniprotid[count] <- idmappedjson[1,2]
  } else {
    print("No protein mapped to gene")
  }
  
  #sleep to avoid rate limit
  Sys.sleep(sleeptime)
  count <- count + 1
}


################# Querying UNIPROT PROTEINS API #################
#query uniprot proteins API using uniprot ids and gather info
uniprotapibase <- "www.ebi.ac.uk/proteins/api/proteins/"

count=1
for (query in querygenes){
  #check if protein coding
  if(is.na(uniprotid[count]) == FALSE){
    #call api and import data
    cat("Getting Uniprot Proteins API data from", query, ":\t", uniprotid[count], "\n")
    uniprotapicall <- paste(uniprotapibase, uniprotid[count], sep="")
    uniprotget <- GET(uniprotapicall)
    uniprotget1 <- fromJSON(toJSON(content(uniprotget)))
    
    if (is.null(uniprotget1$protein$submittedName$fullName$value) == FALSE){
      proteinname[count] <- as.character(uniprotget1$protein$submittedName$fullName$value)
    } else {
      proteinname[count] <- as.character(uniprotget1$protein$recommendedName$fullName$value)
    }
    
    #protein sequence
    proteinseq[count] <- uniprotget1$sequence$sequence
    
    #keywords
    keywordconcat <- c()
    for (keywordin in uniprotget1$keywords$value) {
      keywordconcat <- paste(keywordconcat, keywordin, sep=",")
    }
    
    keywords[count] <- keywordconcat
    
  } else {
    cat("No protein data for ", query, ": \t", uniprotid[count], "\n")
  }

  #sleep to avoid rate limit
  Sys.sleep(sleeptime)
  count = count + 1
}


################# Querying BIOGRID TaxonID API  #################
#biogrid api to get taxID
print("Gathering TaxonID from BiogridAPI")
biogridtaxidbase <- "https://webservice.thebiogrid.org/organisms/?accesskey="
biogridtaxidtoken <- "34e9134bdaaefe8eeb2012f1099ed599"
biogridtaxidparams <- "" #&format=json
biogridtaxidapi <- paste(biogridtaxidbase, biogridtaxidtoken, biogridtaxidparams, sep ="")
biogridtaxidget <- GET(biogridtaxidapi)
biogridtaxidgetcontent <- content(biogridtaxidget)
biogridtaxidjson <- read.table(text=biogridtaxidgetcontent, sep="\t")

#biogrid format speciesname of species list
count = 1
for(entry in species){
  species[count] <- gsub("_", " ", capitalize(species[count]))
  count = count + 1
}

#get taxonid for each
count=1
for (queryspecies in species){
  taxcount = 1
  cat("Gathering TaxonID for: ", queryspecies, "\n")
  
  for (taxid in biogridtaxidjson[,2]){
    if(taxid == queryspecies){
      print(taxid)
      taxonid[count] <- biogridtaxidjson[taxcount, 1]
    }
    
    taxcount = taxcount + 1
  }
  
  count = count + 1
}


################# Querying BIOGRID API  #################
#biogrid access token and base api 
biogridtoken <- "&accesskey=34e9134bdaaefe8eeb2012f1099ed599"
biogridbaseapi <- "https://webservice.thebiogrid.org/interactions/"
biogridapiparams <- "?additionalIdentifierTypes=OFFICIAL_SYMBOL&includeInteractors=true&format=jsonExtended&geneList="
biogridapiparamstaxid <- "&taxId="

#query biogrid api to extract information
count=1
for (query in commonname){
  
  biogridapiquery <- paste(biogridbaseapi, biogridapiparams, query, biogridapiparamstaxid, taxonid[count], biogridtoken, sep = "")
  biogridget <- GET(biogridapiquery)
  biogridjson <- fromJSON(toJSON(content(biogridget)))
  cat("Biogrid query for: ", query, "\n")
  
  biogridinteracts <- c()
  for(i in 1:length(biogridjson)){
    
    if (length(biogridjson) == 0){
      biogridinteracts <- NA
    }else {
      biogridinteracts <- c(biogridinteracts, as.character(biogridjson[[i]][8])) 
    }
    
  }
  
  interactswith[count] <- paste(unlist(unique(biogridinteracts)), collapse = ",")
  count = count + 1
  Sys.sleep(sleeptime)
}


################# Querying UNIPROT API  #################
#query Uniprot api
count=1
for (query in commonname){
  
  if(is.na(uniprotid[count]) == FALSE ){
  uniprotquery1 <- query_uniprot(
    query=list("gene_exact"=c(query), "organism_id"=c(taxonid[count])),
    base_url="https://rest.uniprot.org/uniprotkb/",
    columns=c("accession", "gene_primary", "protein_name", 
              "length", "mass", "cc_activity_regulation", "ft_binding", "cc_catalytic_activity", 
              "ft_dna_bind", "ec", "cc_function", "feature_count", "keywordid", "protein_existence",
              "uniparc_id", "cc_subunit", "cc_developmental_stage", "cc_tissue_specificity",
              "go_p", "go_c", "go", "go_f", "go_id", "cc_subcellular_location", "ft_chain",
              "ft_mod_res", "cc_ptm", "protein_families", "xref_refseq", "xref_pdb", "xref_biogrid", "xref_string"),
    max_keys=200,
    updateProgress = NULL,
    show_progress=TRUE
  )
  
  cat("Uniprot query for: ", query, "\n")
  proteinlength[count] <- uniprotquery1$Length[1]
  proteinmass[count] <- uniprotquery1$Mass[1]
  proteinfeatures[count] <- uniprotquery1$Features[1]
  proteinUniParc[count] <- uniprotquery1$UniParc[1]
  proteinsubunit[count] <- uniprotquery1$`Subunit structure`[1]
  proteindevstage[count] <- uniprotquery1$`Developmental stage`[1]
  proteintissuespec[count] <- uniprotquery1$`Tissue specificity`[1]
  proteingo[count] <- uniprotquery1$`Gene Ontology (GO)`[1]
  proteingoid[count] <- uniprotquery1$`Gene Ontology IDs`[1]
  proteinsubcellloci[count] <- uniprotquery1$`Subcellular location [CC]`[1]
  proteinmodres[count] <- uniprotquery1$`Modified residue`[1]
  proteinptm[count] <- uniprotquery1$`Post-translational modification`[1]
  proteinfam[count] <- uniprotquery1$`Protein families`[1]
  proteinrefseq[count] <- uniprotquery1$RefSeq[1]
  proteinpdb[count] <- uniprotquery1$PDB[1]
  proteinbiogrid[count] <- uniprotquery1$BioGRID[1]
  proteinstring[count] <- uniprotquery1$STRING[1]
  proteinfunction[count] <- uniprotquery1$`Function [CC]`[1]
  } else {
    cat("No uniprot entry for ", query, "\n")
  }
  
  count = count + 1
  Sys.sleep(sleeptime)
}

################# Final Dataframe Generation  #################
#Renaming cols
EnsemblID <- querygenes
UniprotID <- uniprotid
GeneName <- commonname
ProteinName <- proteinname
UniparcID <- proteinUniParc
RefseqID <- proteinrefseq
pdbID <- proteinpdb
BiogridID <- proteinbiogrid
StringID <- proteinstring
goID <- proteingoid
goDescription <- proteingo
EnsemblObjectType <- objecttype
GeneDescription <- genedescription
Species <- species
TaxonID <- taxonid
GeneProduct <- biotype
DNASeq <- seq
ChromosomeLoci <- chromosome
Strand <- strand
StartPos <- startpos
EndPos <- endpos
Keywords <- keywords
ProteinSeq <- proteinseq
ProteinFamily <- proteinfam
ProteinLength <- proteinlength
ProteinMass <- proteinmass
ProteinFunction <- proteinfunction
ProteinFeatures <- proteinfeatures
ProteinSubunits <- proteinsubunit
ProteinDevStage <- proteindevstage
ProteinTissueLoci <- proteintissuespec
ProteinSubCellLoci <- proteinsubcellloci
Interactions <- interactswith
ProteinResidueMod <- proteinmodres
ProteinPTM <- proteinptm

#final df
finaldf <- data.frame(EnsemblID,
  UniprotID,
  GeneName,
  GeneProduct,
  ProteinName,
  UniparcID,
  RefseqID,
  pdbID,
  BiogridID,
  StringID,
  goID,
  goDescription,
  EnsemblObjectType,
  GeneDescription,
  Species,
  TaxonID,
  DNASeq,
  ChromosomeLoci,
  Strand,
  StartPos,
  EndPos,
  Keywords,
  ProteinFamily,
  ProteinSeq,
  ProteinLength,
  ProteinMass,
  ProteinFunction,
  ProteinFeatures,
  ProteinSubunits,
  ProteinDevStage,
  ProteinTissueLoci,
  ProteinSubCellLoci,
  ProteinResidueMod,
  ProteinPTM,
  Interactions
)


################# Exporting to MySql  #################

#exporting finaldf to sql

db <-dbConnect(MySQL(),
               user='s2600569',
               password='KlslXlQX',
               dbname='s2600569')

dbWriteTable(db,name='Summary',value=finaldf, overwrite = TRUE)
dbClearResult(dbListResults(db)[[1]])

ac <- dbListConnections(MySQL())
for(con in ac){
  dbDisconnect(con)
}
dbListConnections(MySQL())

