### Information

# This script is used in R to taxonomically assign MOTUs in the Gillet primer dataset for bats and shrews
# This particular write up assesses the MOTUs clustered at 98%

### 1) Set up ###

library(dplyr)
library(phyloseq)
library(knitr)
# the seqRFLP package is required for turning a table into a fasta file
library(seqRFLP)
library(taxize)

### 2) Load data and remove singletons. These are MOTUs with a single read in the entire dataset and are therefore unreliable ###
# Data loaded is output from Step_1__Raw_Sequences_to_MOTUs_table.txt
# Phyloseq is a great package for filtering data, so we enter .csv output from the 'Raw_Sequences_to_MOTU_table.txt' and a sample sheet containing sample information to create a phyloseq object
# For more information on creating a phyloseq object, see this useful tutorial: https://vaulot.github.io/tutorials/Phyloseq_tutorial.html

# Load .csv and create an OTU table

otu_mat <- read.csv("./gillet.sumaclust98.csv") # Read in OTU abundance table
otu_mat <- otu_mat %>% select (-c("total_reads", "sequence")) # remove columns 'total_reads' and 'sequence'
row.names(otu_mat) <- otu_mat$id        # Name rows by OTU ID
otu_mat <- otu_mat %>% select (-id)     # Since we have rownames, We no longer need the 'id' column
otu_mat <- as.matrix(otu_mat)  # Change from dataframe format to Matrix (phyloseq needs it in matrix format)


# Load .csv and create a taxa table

tax_mat <- read.csv("./gillet.sumaclust98.csv") # Read in Taxa table
tax_mat <- tax_mat %>% select (c("id", "sequence")) # remove columns 'id' and 'sequence'
row.names(tax_mat) <- tax_mat$id        # Name rows by OTU ID
tax_mat <- tax_mat %>% select (-id)     # Since we have rownames, We no longer need the 'id' column
tax_mat <- as.matrix(tax_mat) # Change from dataframe format to Matrix (phyloseq needs it in matrix format)

# Load and prepare sample sheet

samples_df <- read.csv("./bat_shrew_trial_sample_sheet.csv")
row.names(samples_df) <- samples_df$Sample   # Name rows by sample ID
samples_df <- samples_df %>% select (-Sample)  # Since we have rownames, We no longer need the 'sample' column

# Create Phyloseq object

OTU = otu_table(otu_mat, taxa_are_rows = TRUE) # Formats otu_table for phyloseq
TAX = tax_table(tax_mat)                       # Formats taxa_table for phyloseq
samples = sample_data(samples_df)              # # Formats sample_sheet for phyloseq

phylo <- phyloseq(OTU, TAX, samples) # create phyloseq object

phylo # just entering the phyloseq object will give its details

# Save phyloseq object

save(phylo, file="./gillet_clust98_unfiltered_notaxa.RData")

# Clean up

rm(list =ls())


### 3) Taxonomic Assignment ###
### Components of this script are taken and altered from:
# https://github.com/tobiasgf/lulu
# Frøslev, T. G., Kjøller, R., Bruun, H. H., Ejrnæs, R., Brunbjerg, A. K., Pietroni, C., & Hansen, A. J. (2017). Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates. Nature Communications, 8(1), 1188.

# Reload phyloseq object for taxonomic assignment

load("./gillet_clust98_unfiltered_notaxa.RData")

# remove sequences with less than 2 reads (i.e. singletons)
ps <- prune_taxa(taxa_sums(phylo) > 1, phylo)
ps

# isolate the taxa table from the phyloseq object

taxadf <- as.data.frame(as.matrix(tax_table(ps)[,"sequence"]))

# Keep only the MOTU identifier column and the sequence column - for blasting
motudf <- as.data.frame(cbind(MOTU = row.names(taxadf), sequence = as.character(taxadf[,"sequence"])))

head(motudf)


# convert table into a fasta file
seqRFLP::dataframe2fas(motudf, file="./MOTU_list_to_blast_98.fasta")

### 3.2) Blast ###
## Next section is done outside of R, but activated from within R
## we use blastn

query_file <- "./MOTU_list_to_blast_98.fasta"
blastn <- "blastn"

system2(blastn,
        c("-db /extrastorage/home/refdbs/blastdb/nt/nt", # path to blastn database
        "-num_threads 23", # how many cores to use
        "-max_target_seqs 25", # maximum sequence targets found to be reported
        "-outfmt '6 std qlen ssciname staxid'", # what information to extract
        "-out clust98.motus.blasthits.txt", # name of output file
        "-qcov_hsp_perc 90", # minimum coverage %
        "-perc_identity 80", # minimum identity %
        "-query", query_file)) # the input/query file


# Use blast results and create tables with desired information for each MOTU

path <- "./"
IDtable_name <- file.path(path,"clust98.motus.blasthits.txt")

IDtable=read.csv(IDtable_name,
                 sep='\t',
                 header=F,
                 as.is=TRUE)

names(IDtable) <- c("qseqid","sseqid","pident",
                    "length","mismatch","gapopen",
                    "qstart","qend","sstart",
                    "send","evalue","bitscore",
                    "qlen","ssciname","staxid")

margin <- 0.51
new_IDtable <- IDtable[0,] # prepare filtered matchlist
ids <- names(table(IDtable$qseqid))
i=1
o=length(ids)
for (name in ids){
  print(paste0("progress: ", round(((i/o) * 100),0) ,"%")) # make a progressline
  test <- IDtable[which(IDtable$qseqid == name),] # select all lines for a query
  max <- max(test$pident)
  test <- test[which(test$pident > (max-margin)),] # select all lines for a query
  #These lines can be included if analysing a taxonomic group with a lot of
  #"unassigned" sequences in GenBank, to exclude those from further evaluation.
  #test2 <- test[!grepl("uncultured eukaryote",
  #          test$truncated_ssciname,ignore.case = TRUE),]
  #if (nrow(test2) > 1) {test <- test2}
  #test <- test[!grepl("Environmental",
  #          test$truncated_ssciname,ignore.case = TRUE),]
  if (nrow(test) > 0 ) { test$string <- toString(names(table(test$ssciname))) }
  new_IDtable = rbind(new_IDtable,test) # add this row to the filtered IDtable
  i=i+1
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

new_IDtable$majority_taxid <-  with(new_IDtable, ave(staxid, qseqid , FUN=Mode))
IDtable2 = new_IDtable[!duplicated(new_IDtable[c(1,17)]),]

### 3.3) Begin using taxize package to use taxid's from blast to retrieve all taxonomic levels ###

all_staxids <- names(table(IDtable2$staxid)) # get all taxids for table
all_classifications <- list() # prepare list for taxize output
o=length(all_staxids) # number of taxids

Start_from <- 1 # change if loop needs to be restarted due to time-out

for (cl in Start_from:o){ # the taxize command "classification" can be run on
  #the all_staxids vector in one line, but often there is
  #a timeout command, therefor this loop workaround.

  #make a progressline (indicating the index the loops needs to be
  #restarted from if it quits)
  print(paste0("processing: ", cl , " of ", o , " taxids"))
  all_classifications[cl] <- classification(all_staxids[cl], db = "ncbi")
}

### 4) Tidy up results ###

output <- data.frame(staxid=character(),taxpath=character(),
                     stringsAsFactors=FALSE)
totalnames <- length(all_staxids)
for (curpart in seq(1:totalnames)){
  print(paste0("progress: ", round(((curpart/totalnames)
                                    * 100),0) ,"%")) # make a progressline
  currenttaxon <- all_classifications[curpart][[1]]
  if ( !is.na(currenttaxon)) {
    spec <- all_staxids[curpart]
    gen <- currenttaxon[which(currenttaxon$rank == "genus"),"name"]
    fam <- currenttaxon[which(currenttaxon$rank == "family"),"name"]
    ord <- currenttaxon[which(currenttaxon$rank == "order"),"name"]
    cla <- currenttaxon[which(currenttaxon$rank == "class"),"name"]
    phy <- currenttaxon[which(currenttaxon$rank == "phylum"),"name"]
    kin <- currenttaxon[which(currenttaxon$rank == "kingdom"),"name"]
    spe <- currenttaxon[which(currenttaxon$rank == "species"),"name"]
    currentpath <- gsub(" ", "_",
                        #paste0("k__",kin,";p__",phy,";c__",cla,";o__",ord,";f__",fam,";g__",gen,";s__",spe))
                        paste0(kin,";",phy,";",cla,";",ord,";",fam,";",gen,";",spe))
    output[curpart,"staxid"] <-  spec # add row to the filtered IDtable
    output[curpart,"taxpath"] <-  currentpath # add row to the filtered IDtable
  }
}


taxonomic_info <- merge(IDtable2,output,by = "staxid", all=TRUE)
tbname <- file.path(path,"Table_otu_clust98_taxonomy1.txt")
{write.table(taxonomic_info, tbname, sep="\t",quote=FALSE, col.names = NA)}


#Split taxonomic string into levels for the OTU data
tab_name <- file.path(path,"Table_otu_clust98_taxonomy1.txt")
otutaxonomy <- read.table(tab_name, sep="\t", header=TRUE, as.is=TRUE)
library(stringr)
otulevels <- str_split_fixed(otutaxonomy$taxpath, ";", 7)
otulevels <- gsub(".__","",otulevels)
otulevels <- as.data.frame(otulevels)
names(otulevels) <- c("kingdom","phylum","class","order","family","genus",
                      "species")
otutaxlevels <- cbind(otutaxonomy,otulevels)
tab_name <- file.path(path,"Table_otu_clust98_taxonomy_levels1.txt")

# Save table with taxonomic assignment
{write.table(otutaxlevels, tab_name, sep="\t",quote=FALSE, col.names = NA)}

### Clean up taxonomic table and remove restrict assignment of MOTUs with low percentage of identity
# Everything has been taxonomically assigned
# All information is in the file "Table_otu_clust*_taxonomy_levels1.txt"

tax.ass.tab <- read.csv("./Table_otu_clust98_taxonomy_levels1.txt", sep="\t", header = TRUE, as.is=TRUE)

taxonomy_table <- tax.ass.tab[,c("qseqid",  # the query sequence ID - i.e. the MOTU value/ID
                                    "kingdom",
                                    "phylum",
                                    "class",
                                    "order",
                                    "family",
                                    "genus",
                                    "species",
                                    "pident",
                                    "evalue")]



newtaxatable <- as.data.frame(tax_table(ps))
newtaxatable$qseqid <- rownames(newtaxatable)

newtaxatable <- merge(taxonomy_table, newtaxatable, by= "qseqid", all=TRUE)
rownames(newtaxatable) <- newtaxatable$qseqid
newtaxatable <- newtaxatable %>% select (-qseqid)

library(stringr)

tax.clean <- as.data.frame(newtaxatable)
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- "none" # if the values in this dataframe is NA, replace with blank i.e. ""

for (i in 1:nrow(tax.clean)){

  #Fill in missing taxonomy

  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

# save(tax.clean, file="./tax.clean.clust95.preidenfilt.csv")

### The extra filter to restrict taxa assignment of MOTUs with <98% identity

tax.clean$pident <- as.numeric(tax.clean$pident)

tax.clean[is.na(tax.clean)] <- 0

for (i in 1:nrow(tax.clean)){

  # identity percentage filter
   if (tax.clean[i,"pident"] == 0){
    tax.clean[i, 1:7] <- "none"
   } else if (tax.clean[i,"pident"] < 90){
    class <- paste("Class_", tax.clean[i,"class"], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,"pident"] < 93){
    order <- paste("Order_", tax.clean[i,"order"], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,"pident"] < 95){
    family <- paste("Family_", tax.clean[i,"family"], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,"pident"] < 98){
    genus <- paste("Genus_", tax.clean[i,"genus"], sep = "")
    tax.clean[i, "species"] <- genus
  }
}


### 5) Merge new, cleaned taxonomic table to phyloseq object created at the start of this script ###

new.tax.mat <- as.matrix(tax.clean)
gillet.phylo <- phyloseq(tax_table(new.tax.mat),
               sample_data(sample_data(ps)),
               otu_table(otu_table(ps),
                         taxa_are_rows = TRUE))
# Check new phyloseq object
gillet.phylo

# Save Phyloseq Object
save(gillet.phylo, file="./phyloseq_object_clust_iden98_taxa_assignednone_no_singletons.RData")


### 6) Empirically comparing taxonomic assignment from different clustering thresholds (95% - 98% in this case) is recommended
# See Step_3__clustering_threshold_comparison.html for script on comparing these results
